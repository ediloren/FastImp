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

  const static char cvsid[] = "$Id: meshRing.h,v 1.4 2002/12/01 15:50:06 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _MESH_RING_H_
#define _MESH_RING_H_

#include "meshCond.h"
#include "ring.h"

namespace mesh {

  class MeshRing : public MeshCond {

  public:

    MeshRing(void) {}
    MeshRing(const surf::Ring* RingPtr, const bool useNonUniformMesh,
	     const bool useQuadPanel, const bool aspectRatioWarning);

  private:
    void mapToCube (const surf::Ring* RingPtr);
    void mapPolarToCard (void);

  };

} // namespace mesh

#endif
