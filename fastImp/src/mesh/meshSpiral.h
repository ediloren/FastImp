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

  const static char cvsid[] = "$Id: meshSpiral.h,v 1.4 2002/12/01 15:50:06 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _MESH_SPIRAL_H_
#define _MESH_SPIRAL_H_

#include "meshCond.h"
#include "spiral.h"

namespace mesh {

  class MeshSpiral : public MeshCond {

  public:

    MeshSpiral(void) {}
    MeshSpiral(const surf::Spiral* SpiralPtr, const bool useNonUniformMesh,
	       const bool useQuadPanel, const bool aspectRatioWarning);

  private:
    void mapToCube (const surf::Spiral* SpiralPtr);
    void mapPolarToCard (const surf::Spiral* SpiralPtr);

  };

} // namespace mesh

#endif
