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

const static char cvsid[] = "$Id: g2gPoint.cc,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

#include "g2gPoint.h"
#include <iostream>

using namespace pfft;

/**********************************************************************
 * operator < --
 **********************************************************************/
bool 
pfft::operator < (
		  const G2GPoint& p1, 
		  const G2GPoint& p2) 
{
  return p1.pos() < p2.pos(); 
}


/**********************************************************************
 * operator != --
 **********************************************************************/
bool 
pfft::operator != (
		   const G2GPoint& p1, 
		   const G2GPoint& p2) 
{
  return p1.pos() != p2.pos(); 
}

/**********************************************************************
 * operator == --
 **********************************************************************/
bool 
pfft::operator == (
		   const G2GPoint& p1, 
		   const G2GPoint& p2) 
{
  return p1.pos() == p2.pos(); 
}


/**********************************************************************
 * output --
 **********************************************************************/
std::ostream& 
operator << (
	     std::ostream& os, 
	     const G2GPoint& p) 
{
  os <<  "{" << "(" << p.pos().x() << "," << p.pos().y() << "," << p.pos().z() << ")"
     << ", " << p.directStencilPointIndex()
     << ", " << p.interpStencilPointIndex() 
     << "}";
  return os;
}
