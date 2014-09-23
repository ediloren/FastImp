/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Zhenhai Zhu, Ben Song
  
  Description:

  Resources:

  See also:

  const static char cvsid[] = "$Id: gridData.h,v 1.4 2002/10/01 21:17:37 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _GRID_DATA_H_
#define _GRID_DATA_H_

#include "TArray.h"
#include "x3dconv_rr.h"
#include "x3dconv_rc.h"
#include "x3dconv_cc.h"
#include "spRowMat.h"
#include "spColMat.h"

namespace pfft {

  template <class KernelType, class DataType, class GreenFunc>
  class GridData : public Fast3DConv<KernelType, DataType, GreenFunc>  {
    
  public:
    GridData(void) {}
    GridData(const GreenFunc& greenFuncIn) : greenFunc_(greenFuncIn) {}

    size_t memoryUsage (size_t, size_t, size_t) const;
    void fillKernel(void) { setupKernel(gridKernel_); }
    void fftKernel(void) { fftOnKernel(gridKernel_); }
    void conv3D(void) { operator() (gridKernel_, gridData_, gridData_); }
    void allocate(const size_t nx, const size_t ny, const size_t nz, 
		  const double gridStep)
    {
      gridKernel_.resize(nx, ny, nz);
      gridData_.resize(nx, ny, nz);
      initialize(nx, ny, nz, gridStep, gridStep, gridStep, greenFunc_);
    }

    void outputGridData(std::ostream&) const;

    template <class ProjectMat, class SrcVec>
    void projectOnGrid(const ProjectMat& projectMat, const SrcVec& src) {
      matMultVec(gridData_, projectMat, src);
    }

    template <class InterpMat, class DesVec>
    void interpFromGrid(const InterpMat& interpMat, DesVec& des) const {
      matMultVec(des, interpMat, gridData_);
    }

    const typename TypeSelect<DataType>::TYPE& operator[] (size_t i) const {
      return gridData_[i];
    }

  private:
    GreenFunc greenFunc_;
    TArray<typename TypeSelect<KernelType>::TYPE, KernelTag> gridKernel_;
    TArray<typename TypeSelect<DataType>::TYPE, GridDataTag> gridData_; 

  };
    
  /**********************************************************************
   * memoryEstimate --
   **********************************************************************/
  template <class KernelType, class DataType, class GreenFunc>
  size_t
  GridData<KernelType, DataType, GreenFunc>::
  memoryUsage (
	       size_t numPointX, 
	       size_t numPointY, 
	       size_t numPointZ) const
  {
    size_t nx = 2 * numPointX;
    size_t ny = 2 * numPointY;
    size_t nz = 2 * numPointZ;
    size_t gridDataUsage = sizeof(KernelType) * (nx * ny * numPointZ);
    size_t kernelUsage = sizeof(DataType) * (nx * ny * nz) * 1;
      
    return gridDataUsage + kernelUsage;
  }
  
  /**********************************************************************
   * outputGridData --
   **********************************************************************/
  template <class KernelType, class DataType, class GreenFunc>
  void
  GridData<KernelType, DataType, GreenFunc>::
  outputGridData (std::ostream& out) const
  {
    out<<" \n Grid Data \n"; 
    for (size_t ii=0; ii < gridData_.size(); ++ii) {
      out << gridData_[ii] << std::endl;
    }
    out << std::endl;
  }

} //namespace pfft

#endif








