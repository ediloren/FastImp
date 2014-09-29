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

  When both kernel and data are real, we could save memory by puut them into
  complex arrays.  This involves messy scramble of indexes.

  Resources:

  See also:

  const static char cvsid[] = "$Id: x3dconv_rr.h,v 1.4 2003/02/11 03:09:24 zhzhu Exp $";

  ==========================================================================
*/
#ifndef _X3DCONV_RR_H_
#define _X3DCONV_RR_H_

#include "x3dconv.h"
#include "rfft3d.h"
#include "TArray.h"
#include <vector>
#include <complex>

//Enrico
extern double mydrand(void);

namespace pfft{
  
  template<class T3>
  class Fast3DConv<Real, Real, T3>{

  public:
    inline Fast3DConv(): isKernelSetup_(false), isKernelTransformed_(false),
			 isGridDataSetup_(false), isCalcKernelSetup_(false),
			 calcKernelPtr_(static_cast<T3*>(0)){}

    inline void initialize(const size_t, const size_t, const size_t,
			   const double, const double, const double,
			   const T3&);
    // test module:
    template<class Vector>
    void setupTestGridData(Vector&, bool randomOrNot = true);
  
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
    void setupGridData(const Vector&,GridDataContainer&);
    
    template<class GridDataContainer, class Vector>
    void extractFastPotential(const GridDataContainer&, Vector&);

    template<class KernelContainer, class GridDataContainer>
    void operator() (const KernelContainer&, GridDataContainer&, GridDataContainer&);
		       
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
   * Fast3DConv<Real, Real>::initialize --
   **********************************************************************/
  template<class T3>
  void Fast3DConv<Real, Real, T3>
  ::initialize(const size_t n1, const size_t n2, const size_t n3,
	       const double d1, const double d2, const double d3,
	       const T3 &calcKernel)
  {
    calcKernelPtr_ = &calcKernel;
    isCalcKernelSetup_ = true;

    nx_ = n1; ny_ = n2; nz_ = n3;
    dx_ = d1; dy_ = d2; dz_ = d3;
    pace_[0] = d1; pace_[1] = d2; pace_[2] = d3;

    nx2_ = nx_ * 2; ny2_ = ny_ * 2; nz2_ = nz_ * 2; 
    dim_[0] = nx_; dim_[1] = ny_; dim_[2] = nz_;
    dim2_[0] = nx2_; dim2_[1] = ny2_; dim2_[2] = nz2_;
    length1D_ = nx_ * ny_ * nz_;
    length1D8_ = nx2_ * ny2_ * nz2_;
    length1D4_ = length1D8_ / 2;
  }
  
  /**********************************************************************
   * Fast3DConv<Real, Real>::SetupKernel --
   **********************************************************************/
  template<class T3>
  template<class KernelContainer>
  void Fast3DConv<Real, Real, T3>
  ::setupKernel (KernelContainer& x)
  {
    assert(x.size() == length1D8_);
    assert(isCalcKernelSetup_);

    pfft::setupKernelTopelitz(dim2_, pace_, x, *calcKernelPtr_);
    isKernelSetup_ = true;
  }
		 
  /**********************************************************************
   * Fast3DConv<Real, Real>::Fast3DConv --
   **********************************************************************/
  template<class T3>
  template<class KernelContainer>
  void Fast3DConv<Real, Real, T3>
  ::fftOnKernel(KernelContainer& x)
  {
    assert(isKernelSetup_ == true);
    
    rfft3DGeneral(dim2_, x, FFTW_BACKWARD, NO_SCALE);

    isKernelTransformed_ = true;
  };

  /**********************************************************************
   * Fast3DConv<Real, Real>::Fast3DConv --
   **********************************************************************/
  template<class T3>
  template<class Vector, class GridDataContainer>
  void Fast3DConv<Real, Real, T3>
  ::setupGridData(const Vector& gridData,
		  GridDataContainer& gridDataPadded)
  {
    assert(gridData.size() == length1D_);
    assert(gridDataPadded.size() == length1D4_);

    for(size_t ii=0; ii<nx_; ++ii){
      for(size_t jj=0; jj<ny2_; ++jj){
	for(size_t kk=0; kk<nz2_; ++kk){
	  int paddedIndex = kk + nz2_ * (jj + ny2_ * ii);
	  if(jj<ny_ && kk<nz_){
	    int index = kk + nz_ * (jj + ny_ * ii);
	    gridDataPadded[paddedIndex] = gridData[index];
	  } else {
	    gridDataPadded[paddedIndex] = 0.;
	  }
	}
      }
    } // the grid data is zero padded, all these zeroes should
    // be ensured at its definition !!

    isGridDataSetup_ = true;
  }
 
  /**********************************************************************
   * Fast3DConv<Real, Real>::Fast3DConv --
   **********************************************************************/
  template<class T3>
  template<class GridDataContainer, class Vector>
  void Fast3DConv<Real, Real, T3>
  ::extractFastPotential(const GridDataContainer& gridDataPadded,
			 Vector& gridData)
  {
    assert(gridData.size() == length1D_);
    
    pfft::extractFastPotential(dim2_, gridDataPadded, gridData);
  }  
  
  /**********************************************************************
   * Fast3DConv<Real, Real>::Fast3DConv --
   **********************************************************************/
  template<class T3>
  template<class KernelContainer, class GridDataContainer>
  void Fast3DConv<Real, Real, T3>
  ::operator()(const KernelContainer& kernel,
	       GridDataContainer& grid,
	       GridDataContainer& result)
  {
    assert(isKernelTransformed_ == true);
    // assert(isGridDataSetup_ == true);
    assert(result.size() == length1D4_);

    int stride, dist, times;

    if(&result != &grid){
      copy(grid.begin(), grid.end(), result.begin());
    }

    fftw_plan plan;
 
    fftw_complex* out
      = reinterpret_cast<fftw_complex*>(const_cast<Real*>(result.p_));
    fftw_complex* core
      = reinterpret_cast<fftw_complex*>(const_cast<Real*>(kernel.p_));
    fftw_complex* kspeq
      = reinterpret_cast<fftw_complex*>(const_cast<Real*>(kernel.speq_));

    // transform along x-direction !
    size_t step = ny2_ * nz_;
    times = ny_;
    stride = 1; dist = nz_;

    for(size_t ii=0; ii<nx_; ++ii){
      plan = ::fftw_create_plan_specific(nz_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 out + ii * step,
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);
      ::fftw(plan, times, out + ii*step,
	   stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plan);
    }

    // transform along y-direction !
    times = nz_;
    stride = nz_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
     
      plan = ::fftw_create_plan_specific(ny2_, FFTW_BACKWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 out + ii * step, 
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);
            
      ::fftw(plan, times, out + ii*step,
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plan);
    }
    

    // transform along z-direction !
    size_t j1, j2, j3;

    const Real c1 = 0.5;
    const Real c1_ = 0.5;

    const Real c2 = -0.5;
    const Real c2_ = +0.5;

    const Real theta = 2.0 * PAI / nz2_;
    const Real theta_ = -theta;
    
    Real wtmp = sin(0.5 * theta);
    Real wtmp_ = sin(0.5 * theta_);

    const Real wpr = -2.0 * wtmp * wtmp;
    const Real wpr_ = -2.0 * wtmp_ * wtmp_;

    const Real wpi = sin(theta);
    const Real wpi_ = sin(theta_);

    Real h1r, h1i, h2i, h2r;

    Real *wr = static_cast<Real*> (::operator new( ((nz_>>1)+1) * sizeof(Real), nothrow) );
    assert(wr != 0);

    Real *wi = static_cast<Real*> (::operator new( ((nz_>>1)+1) * sizeof(Real), nothrow) );
    assert(wi != 0);

    Real *wr_ = static_cast<Real*> (::operator new( ((nz_>>1)+1) * sizeof(Real), nothrow) );
    assert(wr_ != 0);

    Real *wi_ = static_cast<Real*> (::operator new( ((nz_>>1)+1) * sizeof(Real), nothrow) );
    assert(wi_ != 0);


    wr[0] = 1.0; wi[0] = 0.0;
    for(size_t kk=1; kk<=(nz_>>1); ++kk){
      wr[kk] = wr[kk-1] * wpr - wi[kk-1]*wpi + wr[kk-1];
      wi[kk] = wi[kk-1] * wpr + wr[kk-1]*wpi + wi[kk-1];
    }
   
    wr_[0] = 1.0; wi_[0] = 0.0;
    for(size_t kk=1; kk<=(nz_>>1); ++kk){
      wr_[kk] = wr_[kk-1] * wpr_ - wi_[kk-1]*wpi_ + wr_[kk-1];
      wi_[kk] = wi_[kk-1] * wpr_ + wr_[kk-1]*wpi_ + wi_[kk-1];
    }

    fftw_complex *x = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(x != 0);

    fftw_complex *y = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    assert(y != 0);


    fftw_plan planPlus = fftw_create_plan(nx2_, FFTW_BACKWARD,
					  FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_plan planMinus = fftw_create_plan(nx2_, FFTW_FORWARD,
					   FFTW_ESTIMATE|FFTW_IN_PLACE);

    fftw_complex *dataLower = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    fftw_complex *dataHigher = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    fftw_complex *speqLower = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));
    fftw_complex *speqHigher = static_cast<fftw_complex*>(::operator new(nx2_ * sizeof(fftw_complex), nothrow));

    stride = nz_ * ny2_;
    size_t xidx2d, yidx2d;
    Real tmpr1, tmpr2, tmpr3, tmpr4;
    Real tmpi1, tmpi2, tmpi3, tmpi4;

    //for(size_t ii3=0, i3=0; i3<(nz2_>>2)+1; i3++, ii3+=2){
    for(size_t i3=0; i3<=(nz_>>1); ++i3){
      j3 = (i3 != 0 ? (nz_ - i3) : 0); 
      for(size_t i2=0; i2<ny2_; ++i2){
	j2 = (i2 != 0 ? (ny2_-i2) : 0);
	
	xidx2d = i3 + nz_*i2; // data 
	yidx2d = j3 + nz_*j2; // data
	
	// step 1: copy to working vector
	for(size_t i1 = 0; i1<nx_; ++i1){ //tmp nx2_ !
	  x[i1] = out[xidx2d + stride * i1];
	  y[i1] = out[yidx2d + stride * i1];
	} // never forget to padding it with zero!!
	for(size_t i1 = 0; i1<nx_; ++i1){
	  x[i1 + nx_].re = 0.0; x[i1 + nx_].im = 0.0;
	  y[i1 + nx_].re = 0.0; y[i1 + nx_].im = 0.0;
	}
 
	// step 2: do the fft.
	// pay attention to special case of i3
	if(i3 != 0 && i3 != (nz_>>1)){
	  fftw_one(planPlus, x, 0);
	  fftw_one(planPlus, y, 0);
	} else {
	  if(i2 <= j2){
	    fftw_one(planPlus, x, 0);
	    fftw_one(planPlus, y, 0);
	  }
	}
	
	// step 3: unscramble the data using the symmetry 
	for(size_t i1=0; i1<nx2_; ++i1){
	  j1 = ((i1 != 0) ? (nx2_ - i1) : 0);
	    
	  if(i3 != 0){
	    h1r = c1 * (x[i1].re + y[j1].re);
	    h1i = c1 * (x[i1].im - y[j1].im);	    
	    h2i = c2 * (x[i1].re - y[j1].re);	    
	    h2r =-c2 * (x[i1].im + y[j1].im);
	  
	    x[i1].re = h1r + wr[i3]*h2r - wi[i3]*h2i; 
	    x[i1].im = h1i + wr[i3]*h2i + wi[i3]*h2r; 
	    y[j1].re = h1r - wr[i3]*h2r + wi[i3]*h2i; 
	    y[j1].im =-h1i + wr[i3]*h2i + wi[i3]*h2r; 
	  
	  } else { // i3 == 0 !
	    if(i2 <= j2){ 
	      h1r = c1 * (x[i1].re + y[j1].re );
	      h1i = c1 * (x[i1].im - y[j1].im);
	      h2i = c2 * (x[i1].re - y[j1].re );
	      h2r = -c2 * (x[i1].im + y[j1].im);

	      dataLower[i1].re = h1r + h2r;
	      dataLower[i1].im = h1i + h2i;
	      
	      speqHigher[j1].re = h1r - h2r;
	      speqHigher[j1].im = h2i - h1i;
	      
	      if(i2 != j2 && i2 != 0){
		h1r = c1 * (y[i1].re + x[j1].re );  
		h1i = c1 * (y[i1].im - x[j1].im);	    
		h2i = c2 * (y[i1].re - x[j1].re );	    
		h2r = -c2 * (y[i1].im + x[j1].im);
	      
		dataHigher[i1].re = h1r + h2r;
		dataHigher[i1].im = h1i + h2i;
	      
		speqLower[j1].re = h1r - h2r;
		speqLower[j1].im = h2i - h1i;
	      }
	    }
	  }
	}
	
	// step 4: multipication in spectral domain
	if(i3 == 0){
	  if(i2 <= j2){
	    for(size_t i1=0; i1<nx2_; ++i1){
	      dataLower[i1] *= core[xidx2d + stride*i1];
	      speqHigher[i1] *= kspeq[j2 + ny2_ * i1];
	      if(i2 != j2 && i2 != 0){
		dataHigher[i1] *= core[yidx2d + stride*i1];
		speqLower[i1] *= kspeq[i2 + ny2_ * i1];
	      }
	    }
	  } 
	}else { // i3 != 0
	  if(i3 != (nz_>>1)){
	    for(size_t i1=0; i1<nx2_; ++i1){
	      x[i1] *= core[xidx2d + stride * i1];
	      y[i1] *= core[yidx2d + stride * i1];
	    }
	  } else {
	    if(i2 <= j2){
	      for(size_t i1=0; i1<nx2_; ++i1){
		x[i1] *= core[xidx2d + stride * i1];
		y[i1] *= core[yidx2d + stride * i1];
	      } 
	    }
	  }
	}
	
	// step 5: unscramble the data again !
	if(i3 == 0){
	  if(i2 <= j2) {
	    for(size_t i1=0; i1<nx2_; ++i1){
	      j1 = ((i1 != 0) ? nx2_ - i1 : 0);
 
	      h1r = c1_ * (dataLower[i1].re + speqHigher[j1].re );
	      h1i = c1_* (dataLower[i1].im - speqHigher[j1].im);
	      h2i = c2_* (dataLower[i1].re - speqHigher[j1].re );
	      h2r = -c2_* (dataLower[i1].im + speqHigher[j1].im);

	      dataLower[i1].re = h1r + h2r;
	      dataLower[i1].im = h1i + h2i;
	      
	      speqHigher[j1].re = h1r - h2r;
	      speqHigher[j1].im = h2i - h1i;
	      
	      if(i2 != j2 && i2 != 0){
		h1r = c1_* (dataHigher[i1].re + speqLower[j1].re );  
		h1i = c1_* (dataHigher[i1].im - speqLower[j1].im);	    
		h2i = c2_* (dataHigher[i1].re - speqLower[j1].re );	    
		h2r = -c2_* (dataHigher[i1].im + speqLower[j1].im);
	      
		dataHigher[i1].re = h1r + h2r;
		dataHigher[i1].im = h1i + h2i;
	      
		speqLower[j1].re = h1r - h2r;
		speqLower[j1].im = h2i - h1i;
	      }
	    }
	  }
	} else { // i3 != 0
	  for(size_t i1=0; i1<nx2_; ++i1){
	    j1 = ((i1 != 0) ? (nx2_ - i1) : 0);
	    
	    h1r = c1_ * (x[i1].re + y[j1].re);
	    h1i = c1_ * (x[i1].im - y[j1].im);	    
	    h2i = c2_ * (x[i1].re - y[j1].re);	    
	    h2r =-c2_ * (x[i1].im + y[j1].im);
	  
	    x[i1].re = h1r + wr_[i3]*h2r - wi_[i3]*h2i; 
	    x[i1].im = h1i + wr_[i3]*h2i + wi_[i3]*h2r; 
	    y[j1].re = h1r - wr_[i3]*h2r + wi_[i3]*h2i; 
	    y[j1].im =-h1i + wr_[i3]*h2i + wi_[i3]*h2r; 
	  }
	}

	// step 6: inverse fft! 
	if(i3 == 0){
	  if(i2 <= j2){
	    fftw_one(planMinus, dataLower, 0);
	    if(i2 != j2 && i2 != 0){
	      fftw_one(planMinus, dataHigher, 0);
	    }
	  }
	} else { // != 0
	  if(i3 != (nz_>>1)){
	    fftw_one(planMinus, x, 0);
	    fftw_one(planMinus, y, 0);
	  } else {
	    if(i2 <= j2){
	      fftw_one(planMinus, x, 0);
	      fftw_one(planMinus, y, 0);
	    }
	  }
	}

	// step 7: copy back
	if(i3 == 0){
	  if(i2 <= j2){
	    for(size_t i1 = 0; i1<nx_; ++i1){ //tmp nx2_ !
	      out[xidx2d + stride * i1] = dataLower[i1];
	      if(i2 != j2 && i2 != 0){
		out[yidx2d + stride * i1] = dataHigher[i1];
	      }
	    }
	  }
	} else {
	  for(size_t i1 = 0; i1<nx_; ++i1){ //tmp nx2_ !
	    out[xidx2d + stride * i1] = x[i1];
	    out[yidx2d + stride * i1] = y[i1];
	  }
	}
      }
    }

    ::operator delete(wr);
    ::operator delete(wi);
    ::operator delete(wr_);
    ::operator delete(wi_);
    ::operator delete(x);
    ::operator delete(y);
    ::operator delete(dataLower);
    ::operator delete(dataHigher);
    ::operator delete(speqLower);
    ::operator delete(speqHigher);

    times = nz_;
    stride = nz_; dist = 1;
    for(size_t ii=0; ii<nx_; ++ii){
      plan = ::fftw_create_plan_specific(ny2_, FFTW_FORWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 out + ii * (ny2_ * nz_),
					 stride,
					 static_cast<fftw_complex*>(0),
					 stride);

      ::fftw(plan, times, out + ii*(ny2_ * nz_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plan);
    }

    times = ny_;
    stride = 1; dist = nz_;
    for(size_t ii=0; ii<nx_; ++ii){
      plan = ::fftw_create_plan_specific(nz_, FFTW_FORWARD,
					 FFTW_ESTIMATE|FFTW_IN_PLACE,
					 out + ii * (ny2_ * nz_),
					 stride,
					 static_cast<fftw_complex*>(0), stride);

      ::fftw(plan, times, out + ii*(ny2_ * nz_),
	     stride, dist, static_cast<fftw_complex*>(0), stride, dist);

      fftw_destroy_plan(plan);
    }

    Real scale_factor = 2./length1D8_;
    size_t length = length1D8_ >> 2;
    for(size_t ii=0; ii<length; ++ii){
      out[ii] *= scale_factor;
    }
  }

  /**********************************************************************
   * Fast3DConv<Real, Real>
   * ::SetupTestGridData<Vector> --
   **********************************************************************/
  template<class T3>
  template<class Vector>
  void Fast3DConv<Real, Real, T3>
  ::setupTestGridData (Vector& gridDataVec, 
		       bool randomOrNot)
  {
    assert(gridDataVec.size() == length1D_);
    for(size_t ii=0; ii<length1D_; ++ii){
      if(randomOrNot){
	gridDataVec[ii] = mydrand();
      } else {
	gridDataVec[ii] = 1.0;
      }
    } 
  }

  /**********************************************************************
   * Fast3DConv<Real, Real>::GetDirectPotential --
   **********************************************************************/
  template<class T3>
  template<class VectorA, class VectorB>
  void Fast3DConv<Real, Real, T3> 
  ::getDirectPotential (const VectorA gridData,
			VectorB& directPotential)
  {
    assert(directPotential.size() == length1D4_);
    assert(isCalcKernelSetup_);

    int tIndex, sIndex;
    for(int tx=0; tx<nx_; ++tx){
      for(int ty=0; ty<ny2_; ++ty){
	for(int tz=0; tz<nz2_; ++tz){
	  tIndex = tz + nz2_ * (ty + ny2_ * tx);
	  directPotential[tIndex] = 0.;
	  if (ty < ny_ && tz < nz_) {
	    for(int sx=0; sx<nx_; ++sx){
	      for(int sy=0; sy<ny2_; ++sy){
		for(int sz=0; sz<nz2_; ++sz){
		  if (sy < ny_ && sz < nz_) {
		    sIndex = sz + nz2_ * (sy + ny2_ * sx);
		    directPotential[tIndex] 
		      += (*calcKernelPtr_)(fabs(tx-sx)*pace_[0], fabs(ty-sy)*pace_[1], fabs(tz-sz)*pace_[2]) * gridData[sIndex];
		  }
		}
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
  bool Fast3DConv<Real, Real, T3>
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









