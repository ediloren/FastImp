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

  const static char cvsid[] = "$Id: meshWire.h,v 1.4 2002/12/01 15:50:06 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _MESH_WIRE_H_
#define _MESH_WIRE_H_

#include "meshCond.h"
#include "wire.h"

namespace mesh {

  class MeshWire : public MeshCond {

  public:

    MeshWire(void) {}
    MeshWire(const surf::Wire* WirePtr, const bool useNonUniformMesh,
	     const bool useQuadPanel, const bool aspectRatioWarning);

  private:
    void mapToCube (const surf::Wire* WirePtr);

  };

} // namespace mesh

#endif
