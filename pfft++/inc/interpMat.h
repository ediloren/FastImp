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

  const static char cvsid[] = "$Id: interpMat.h,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _INTERP_MAT_H_
#define _INTERP_MAT_H_

#include "stencil.h"
#include "gridElement.h"
#include "kernelIntegration.h"
#include "grid.h"
#include "spRowMat.h"
#include <vector>

namespace pfft {

  class InterpMat {

  public:

    InterpMat(void) {}
    InterpMat(
	      Stencil* const interpStencilIn,
	      GridElement* const evalGridElementIn,
	      KernelIntegration* const integrationIn,
	      Grid* const gridIn)
      : interpStencilPtr(interpStencilIn), evalGridElementPtr(evalGridElementIn), 
      integrationPtr(integrationIn), gridPtr(gridIn)
    {
      matList.resize(integrationPtr->numOuterOperator());
    }

    size_t memoryEstimation(void) const;
    void setup(void);
    const size_t numMat(void) const { return matList.size(); }
    const SpRowMat<double>& operator[] (size_t matIndex) const { 
      return matList[matIndex]; }  
    SpRowMat<double>& operator[] (size_t matIndex) { return matList[matIndex]; }  

  private:
    Stencil* interpStencilPtr;
    GridElement* evalGridElementPtr;
    KernelIntegration* integrationPtr;
    Grid* gridPtr;
    std::vector<SpRowMat<double> > matList;

    void setupEachMatrix (const DifferentialOperator, const size_t);

  };

} //namespace pfft 

#endif
