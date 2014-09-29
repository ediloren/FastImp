/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Ben Song
  
  Description:

  Resources:

  See also:

  const static char cvsid[] = "$Id: x3dconv_rc.h,v 1.4 2003/02/11 03:09:24 zhzhu Exp $";

  ==========================================================================
*/
#ifndef _X3DCONV_RC_H_
#define _X3DCONV_RC_H_

#include "x3dconv.h"
#include "TArray.h"
#include "rfft3d.h"

//Enrico
extern double mydrand(void);

namespace pfft{
  
  template<class T3>
  class Fast3DConv<Real, Cplx, T3>{
    
   public:
    inline Fast3DConv(): isKernelSetup_(false), isKernelTransformed_(false),
			 isGridDataSetup_(false), isCalcKernelSetup_(false),
			 calcKernelPtr_(static_cast<T3*>(0)) {}
    
    inline void initialize(const size_t, const size_t, const size_t,
			   const double, const double, const double,
			   const T3&);

     // test module:
     template<class Vector>
     void setupTestGridData(Vector&, bool = true);
  
     template<class VectorA, class VectorB>
     void getDirectPotential(const VectorA, VectorB&);

     template<class VectorA, class VectorB>
     bool compare(const VectorA&, const VectorB&);
    
    // pfft module:
     template<class KernelContainer>
     void setupKernel(KernelContainer&);

     template<class KernelContainer>
     void fftOnKernel(KernelContainer&);

     template<class Vector, class GridDataContainer>
     void setupGridData(const Vector&, GridDataContainer&);
     
     template<class GridDataContainer, class Vector>
     void extractFastPotential(const GridDataContainer&, Vector&);

     template<class KernelContainer, class GridDataContainer>
     void operator() (const KernelContainer&,
		      GridDataContainer&,
		      GridDataContainer&);
		       
   protected:
     // state variables
     size_t nx_, ny_, nz_;
     Real dx_, dy_, dz_;

     bool isKernelSetup_, isKernelTransformed_;
     bool isGridDataSetup_, isCalcKernelSetup_;

     const T3 *calcKernelPtr_;
     
     // working variables
     size_t nx2_, ny2_, nz2_; // 2 means Real 
     size_t dim_[3], dim2_[3];
     Real pace_[3];
     size_t length1D_, length1D8_, length1D4_;
   };

   /**********************************************************************
    * Fast3DConv<Real, Cplx, T3>::initialize --
    **********************************************************************/
   template<class T3>
   void Fast3DConv<Real, Cplx, T3>
   ::initialize(const size_t n1, const size_t n2, const size_t n3,
		const double d1, const double d2, const double d3,
		const T3& calcKernel)
   {
     nx_ = n1; ny_ = n2; nz_ = n3;
     dx_ = d1; dy_ = d2; dz_ = d3;

     calcKernelPtr_ = &calcKernel;
     isCalcKernelSetup_ = true;
     
     pace_[0] = d1; pace_[1] = d2; pace_[2] = d3;

     nx2_ = nx_ * 2; ny2_ = ny_ * 2; nz2_ = nz_ * 2; 
     dim_[0] = nx_; dim_[1] = ny_; dim_[2] = nz_;
     dim2_[0] = nx2_; dim2_[1] = ny2_; dim2_[2] = nz2_;
     length1D_ = nx_ * ny_ * nz_;
     length1D8_ = nx2_ * ny2_ * nz2_;
     length1D4_ = length1D8_ / 2;
  }
  
   /**********************************************************************
    * Fast3DConv<Real, Cplx, T3>::SetupKernel --
    **********************************************************************/
   template<class T3>
   template<class KernelContainer>
   void Fast3DConv<Real, Cplx, T3>
   ::setupKernel (KernelContainer& x)
   {
     assert(x.size() == length1D8_);
     assert(isCalcKernelSetup_);
     
     pfft::setupKernelTopelitz(dim2_, pace_, x, *calcKernelPtr_);
     isKernelSetup_ = true;
   }
   
  /**********************************************************************
   * Fast3DConv<Real, Cplx, T3>::Fast3DConv --
   **********************************************************************/
   template<class T3>
   template<class KernelContainer>
   void Fast3DConv<Real, Cplx, T3>
   ::fftOnKernel(KernelContainer& x)
   {
     assert(isKernelSetup_ == true);
     
     rfft3DGeneral(dim2_, x, FFTW_BACKWARD, NO_SCALE);
     
     isKernelTransformed_ = true;
   };

  /**********************************************************************
   * Fast3DConv<Real, Cplx, T3>::Fast3DConv --
   **********************************************************************/
   template<class T3>
   template<class VectorA, class VectorB>
   void Fast3DConv<Real, Cplx, T3>
   ::setupGridData(const VectorA& gridData,
		   VectorB& gridDataPadded)
   {
     assert(gridData.size() == length1D_);
     assert(gridDataPadded.size() == length1D4_);
     
     for(size_t ii=0; ii<nx_; ++ii){
       for(size_t jj=0; jj<ny2_; ++jj){
	 for(size_t kk=0; kk<nz2_; ++kk){
	   int paddedIndex = kk + nz2_ * (jj + ny2_ * ii);
	   if(jj<ny_ && kk<nz_){
	     int index = kk + nz_ * (jj + ny_ * ii);
	     gridDataPadded[paddedIndex].re = gridData[index].real();
	     gridDataPadded[paddedIndex].im = gridData[index].imag();
	   } else {
	     gridDataPadded[paddedIndex].re = 0.;
	    gridDataPadded[paddedIndex].im = 0.;
	   }
	 }
       }
     } // the grid data is zero padded, all these zeroes should
     // be ensured at its definition !!
     
     isGridDataSetup_ = true;
   }
   
   /**********************************************************************
    * Fast3DConv<Real, Cplx, T3>::Fast3DConv --
    **********************************************************************/
   template<class T3>
   template<class GridDataContainer, class Vector>
   void Fast3DConv<Real, Cplx, T3>
   ::extractFastPotential(const GridDataContainer& gridDataPadded,
			  Vector& gridData)
   {
     assert(gridData.size() == length1D_); 
     pfft::extractFastPotential(dim2_, gridDataPadded, gridData);
   }
   
   /**********************************************************************
   * Fast3DConv<Real, Cplx, T3>::Fast3DConv --
   **********************************************************************/
   template<class T3>
   template<class KernelContainer, class GridDataContainer>
   void Fast3DConv<Real, Cplx, T3>
   ::operator()(const KernelContainer& kernel,
		GridDataContainer& grid,
		GridDataContainer& result)
   {
     assert(isKernelTransformed_ == true);
     assert(result.size() == length1D4_);
     
     int stride, dist, times;
     
     if(&grid != &result){
       copy(grid.begin(), grid.end(), result.begin());
     }
     
     fftw_plan plan;
     
     // transform along z-direction !
     times = ny_;
     stride = 1; dist = nz2_;
     for(size_t ii=0; ii<nx_; ++ii){
       plan = ::fftw_create_plan_specific(nz2_, FFTW_BACKWARD,
					  FFTW_ESTIMATE|FFTW_IN_PLACE,
					  result.p_ + ii * (ny2_ * nz2_),
					  stride,
					  static_cast<fftw_complex*>(0),
					  stride);

       ::fftw(plan, times, result.p_ + ii*(ny2_ * nz2_),
	      stride, dist, static_cast<fftw_complex*>(0), stride, dist);

       fftw_destroy_plan(plan);
     }
    
     times = nz2_;
     stride = nz2_; dist = 1;
     for(size_t ii=0; ii<nx_; ++ii){
       plan = ::fftw_create_plan_specific(ny2_, FFTW_BACKWARD,
					  FFTW_ESTIMATE|FFTW_IN_PLACE,
					  result.p_ + ii * (ny2_ * nz2_),
					  stride,
					  static_cast<fftw_complex*>(0),
					  stride);
       
       ::fftw(plan, times, result.p_ + ii*(ny2_ * nz2_),
	      stride, dist, static_cast<fftw_complex*>(0), stride, dist);
       
       fftw_destroy_plan(plan);
     }
     
     // 3-in-1
     // combined the last dimension transform(1d), multiplication, inverse transform,
     // then just put the half of the result back.
     fftw_complex *work = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
     assert(work != 0);
     
     for(size_t ii=0; ii<nx2_; ++ii){
       work[ii].re = 0.;
       work[ii].im = 0.;
     }
     
     fftw_plan planPlus = fftw_create_plan(nx2_, FFTW_BACKWARD,
					   FFTW_ESTIMATE|FFTW_IN_PLACE);
     fftw_plan planMinus = fftw_create_plan(nx2_, FFTW_FORWARD,
					    FFTW_ESTIMATE|FFTW_IN_PLACE);
     
     stride = ny2_ * nz2_;
     size_t index2D, idx0, idxr, idxi;
     for(size_t jj=0; jj<ny2_; ++jj){
       for(size_t kk=0; kk<nz2_; ++kk){
	 index2D = kk + nz2_ * jj;
	 // step 1: copy into work vector
	 for(size_t ii=0; ii<nx_; ++ii){
	   work[ii] = result.p_[index2D + ii*stride];
	 }
	 for(size_t ii=0; ii<nx_; ++ii){
	   work[ii + nx_].re = 0.;
	   work[ii + nx_].im = 0.;
	 }
	 // step 2: transform it
	 fftw_one(planPlus, work, 0);
	 // step 3: multiplication with kernel spectrum
	 for(size_t ii=0; ii<nx2_; ++ii){
	   if(kk < nz_){
	     idxr = (kk<<1) + nz2_*(jj + ny2_*ii);
	     idxi = ((kk<<1)+1) + nz2_*(jj + ny2_*ii);
	     work[ii] *= Cplx(kernel.p_[idxr], kernel.p_[idxi]);
	   } else if(kk == nz_) {
	     idxr = (jj<<1) + (2*ny2_)*ii;
	     idxi = ((jj<<1)+1) + (2*ny2_)*ii;
	     work[ii] *= Cplx(kernel.speq_[idxr], kernel.speq_[idxi]);
	   } else { // kk > nz_
	     idxr = ((nz2_ - kk)<<1) + nz2_ * ((jj == 0 ? 0:ny2_ - jj)
					       + ny2_ * (ii == 0 ? 0:nx2_ - ii));
	    
	     idxi = (((nz2_ - kk)<<1)+1) + nz2_ * ((jj == 0 ? 0:ny2_ - jj)
						   + ny2_*(ii == 0 ? 0:nx2_ - ii));
	  
	     work[ii] *= Cplx(kernel.p_[idxr], kernel.p_[idxi]);
	   }
	 }
	 // step 4: inverse transform
	 fftw_one(planMinus, work, 0);
	 // step 5: copy half back
	 for(size_t ii=0; ii<nx_; ++ii){
	  result.p_[index2D + ii*stride] = work[ii];
	 } 
       }
     }
     fftw_destroy_plan(planPlus);
     fftw_destroy_plan(planMinus);
     
     times = nz2_;
     stride = nz2_; dist = 1;
     for(size_t ii=0; ii<nx_; ++ii){
       plan = ::fftw_create_plan_specific(ny2_, FFTW_FORWARD,
					  FFTW_ESTIMATE|FFTW_IN_PLACE,
					  result.p_ + ii * (ny2_ * nz2_),
					  stride, 
					  static_cast<fftw_complex*>(0),
					  stride);
      
       ::fftw(plan, times, result.p_ + ii*(ny2_ * nz2_),
	      stride, dist, static_cast<fftw_complex*>(0), stride, dist);
       
       fftw_destroy_plan(plan);
     }

     times = ny_;
     stride = 1; dist = nz2_;
     for(size_t ii=0; ii<nx_; ++ii){
       plan = ::fftw_create_plan_specific(nz2_, FFTW_FORWARD,
					  FFTW_ESTIMATE|FFTW_IN_PLACE,
					  result.p_ + ii * (ny2_ * nz2_),
					  stride,
					  static_cast<fftw_complex*>(0),
					  stride);

       ::fftw(plan, times, result.p_ + ii*(ny2_ * nz2_),
	      stride, dist, static_cast<fftw_complex*>(0), stride, dist);
       
       fftw_destroy_plan(plan);
     }
    
     Real scale_factor = 1./length1D8_;
     for(size_t ii=0; ii<length1D4_; ++ii){
       result.p_[ii] *= scale_factor;
     }
     
     ::operator delete(work);
   }   
   
  /**********************************************************************
   * Fast3DConv<Real, Cplx, T3>
   * ::SetupTestGridData<Vector> --
   **********************************************************************/
   template<class T3>
   template<class Vector>
   void Fast3DConv<Real, Cplx, T3>
   ::setupTestGridData (Vector& gridDataVec, bool randomOrNot)
   {
     assert(gridDataVec.size() == length1D_);
     for(size_t ii=0; ii<length1D_; ++ii){
       if(randomOrNot){
	 gridDataVec[ii] = Cplx(mydrand(), mydrand());
       } else {
	 gridDataVec[ii] = Cplx(1., 0.);
       }
     } 
   }

  /**********************************************************************
   * Fast3DConv<Real, Cplx, T3>::GetDirectPotential --
   **********************************************************************/
   template<class T3>
   template<class VectorA, class VectorB>
   void Fast3DConv<Real, Cplx, T3> 
   ::getDirectPotential (const VectorA gridData,
			 VectorB& directPotential)
   {
     assert(directPotential.size() == length1D_);
     assert(isCalcKernelSetup_);
     
     int tIndex, sIndex;
     for(int tx=0; tx<nx_; ++tx){
       for(int ty=0; ty<ny_; ++ty){
	 for(int tz=0; tz<nz_; ++tz){
	   int tIndex = tz + nz_ * (ty + ny_ * tx);
	   directPotential[tIndex] = 0.;
	   for(int sx=0; sx<nx_; ++sx){
	     for(int sy=0; sy<ny_; ++sy){
	       for(int sz=0; sz<nz_; ++sz){
		 int sIndex = sz + nz_ * (sy + ny_ * sx);
		 directPotential[tIndex] 
		   += (*calcKernelPtr_)((tx-sx)*pace_[0], (ty-sy)*pace_[1], (tz-sz)*pace_[0]) * gridData[sIndex];
	       }
	     }
	   }
	 }
       }
     }
   }
   
  
   /**********************************************************************
    * compare --
    **********************************************************************/
   template<class T3>
   template<class VectorA, class VectorB>
   bool Fast3DConv<Real, Cplx, T3>
   ::compare (
	      const VectorA& directPotential,
	      const VectorB& fastPotential)
	     
   {
     assert(directPotential.size() == fastPotential.size());
     
     Real TOLERANCE = 1.e-6;
     Real reErr;
     
     typename VectorA::const_iterator iterA = directPotential.begin();
     typename VectorA::const_iterator iterB = fastPotential.begin();
     while(iterA != directPotential.end()){
       reErr = abs(*iterA - *iterB) / abs(*iterA);
       if(reErr > TOLERANCE){
	 return false;
       }
       ++iterA; ++iterB;
     }
     
     return true;
   }

 } // namespace pfft

#endif




