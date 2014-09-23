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

  const static char cvsid[] = "$Id: meshOtherShape.h,v 1.2 2002/07/18 15:05:45 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _MESH_OTHER_SHAPE_H_
#define _MESH_OTHER_SHAPE_H_

#include "meshCond.h"

namespace mesh {

  class MeshOtherShape : public MeshCond {

  public:

    MeshOtherShape(void) { }
    MeshOtherShape(const int numNode, const int numPanel): 
      MeshCond(numNode, numPanel) {}

  private:

  };

} //namespace mesh

#endif
