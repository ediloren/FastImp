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

  const static char cvsid[] = "$Id: projectMat.h,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _PROJECT_MAT_H_
#define _PROJECT_MAT_H_

#include "kernelIntegration.h"
#include "stencil.h"
#include "gridElement.h"
#include "grid.h"
#include "spColMat.h"
#include <vector>

namespace pfft {

  class ProjectMat {

  public:

    ProjectMat(void) {};
    ProjectMat(
	       Stencil* const projectStencil,
	       GridElement* const srcGridElement,
	       KernelIntegration* const integration,
	       Grid* const grid)
      : projectStencilPtr(projectStencil), srcGridElementPtr(srcGridElement), 
      integrationPtr(integration), gridPtr(grid)
    {
      matList.resize(integrationPtr->numInnerOperator());     
    }

    size_t memoryEstimation(void) const;
    void setup(void);
    const size_t numMat(void) const { return matList.size(); }
    const SpColMat<double>& operator[] (size_t matIndex) const { 
      return matList[matIndex]; }  
    SpColMat<double>& operator[] (size_t matIndex) { return matList[matIndex]; }  

  private:

    Stencil* projectStencilPtr;
    GridElement* srcGridElementPtr;
    KernelIntegration* integrationPtr;
    Grid* gridPtr;
    std::vector<SpColMat<double> > matList;

    void setupEachMatrix (const DifferentialOperator, const size_t);
  };

} //namespace pfft

#endif
