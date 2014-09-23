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

  These are shared function by x3dconv_rr, x3dconv_rc and x3dconv_cc

  Resources:

  See also:

  const static char cvsid[] = "$Id: x3dconv.h,v 1.4 2003/02/11 03:09:24 zhzhu Exp $";

  ==========================================================================
*/
#ifndef _X3DCONV_H_
#define _X3DCONV_H_

#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <complex>
#include <string>
#include <cmath>
#include <functional>
#include "fftwInterface.h"
#include "utils.h"

namespace pfft{

  using namespace std;

  typedef complex<double> Cplx; 
  typedef double Real;

  template<class T>
  struct TypeSelect{
    typedef T TYPE;
  };

  template<>
  struct TypeSelect<Cplx>{
    typedef ::fftw_complex TYPE;
  };
 
  //  const double PAI = 3.1415926535897932385;
  const bool DO_SCALE = true;
  const bool NO_SCALE = false;

  struct KernelTag{};
  struct GridDataTag{};

  // T1: Kernel Type
  // T2: GridData Type
  // T3: CalcKernel Type
  template<class T1, class T2, class T3>
  class Fast3DConv;

  template<class T, class Tag>  class TArray;

  /**********************************************************************
   * circulantIndex --
   **********************************************************************/
  inline int circulantIndex(int index, int maxSize)
  {
    if(index == maxSize/2 ){
      return 0;
    } else {
      return ( (index < maxSize/2) ? index : (maxSize-index) );
    }
  }

  /**********************************************************************
   * operator*= --
   **********************************************************************/
  struct fftw_complex2Cplx : public unary_function<fftw_complex, Cplx>{
    Cplx operator()(fftw_complex x){return Cplx(x.re, x.im);}
  };

  
  /**********************************************************************
   * operator*= --
   **********************************************************************/
  struct Cplx2fftw_complex : public unary_function<Cplx, fftw_complex>{
    fftw_complex operator()(Cplx x){
	fftw_complex a;
	a.re = x.real();
	a.im = x.imag();
	return a;
    }
  };

  /**********************************************************************
   * conj --
   **********************************************************************/
  inline fftw_complex conj(fftw_complex& a)
  {
    fftw_complex b;
    b.re = a.re; b.im = -a.im;
    return b;
  }

  inline double conj(const double& a)
  {
    return a;
  }

  /**********************************************************************
   * operator * --
   **********************************************************************/
  inline fftw_complex operator* (const fftw_complex& x, const Cplx& y)
  {
    fftw_complex a;
    a.re = x.re * y.real() - x.im * y.imag();
    a.im = x.re * y.imag() + x.im * y.real();
    return a;
  }

/**********************************************************************
   * operator * --
   **********************************************************************/
  inline fftw_complex operator* (const double& op1, const fftw_complex& op2)
  {
    fftw_complex a;
    a.re = op1 * op2.re;
    a.im = op1 * op2.im;
    return a;
  }

  /**********************************************************************
   * operator += --
   **********************************************************************/
  inline fftw_complex& operator+= (fftw_complex& lhs, const fftw_complex& rhs)
  {
    lhs.re += rhs.re;
    lhs.im += rhs.im;
    return lhs;
  }

  /**********************************************************************
   * operator += --
   **********************************************************************/
  inline fftw_complex& operator+= (fftw_complex& lhs, const std::complex<double>& rhs)
  {
    lhs.re += rhs.real();
    lhs.im += rhs.imag();
    return lhs;
  }

  /**********************************************************************
   * operator += --
   **********************************************************************/
  inline Cplx& operator+= (Cplx& lhs, const fftw_complex& rhs)
  {
    lhs += Cplx(rhs.re, rhs.im);
    return lhs;
  }

  /**********************************************************************
   * operator*= --
   **********************************************************************/
  inline fftw_complex& operator*= (fftw_complex& lhs, const fftw_complex& rhs)
  {
    double tmpr = lhs.re; 
    double tmpi = lhs.im;
    lhs.re = tmpr*rhs.re - tmpi*rhs.im;
    lhs.im = tmpr*rhs.im + tmpi*rhs.re;
    return lhs;
  }

  inline fftw_complex& operator*= (fftw_complex& lhs, const complex<double>& rhs)
  {
    double tmpr = lhs.re; 
    double tmpi = lhs.im;
    lhs.re = tmpr*rhs.real() - tmpi*rhs.imag();
    lhs.im = tmpr*rhs.imag() + tmpi*rhs.real();
    return lhs;
  }

  /**********************************************************************
   * operator*= --
   **********************************************************************/
  inline fftw_complex& operator*= (fftw_complex& lhs, const double& rhs)
  {
    lhs.re *= rhs;
    lhs.im *= rhs;
    return lhs;
  }

  /**********************************************************************
   * operator* --
   **********************************************************************/
  inline void multiply(const fftw_complex& a, const std::complex<double>& b,
		       fftw_complex& result)
  {
    result.re = a.re * b.real() - a.im * b.imag();
    result.im = a.re * b.imag() + a.im * b.real();
  }

  /**********************************************************************
   * operator<< --
   **********************************************************************/
  inline std::ostream& operator<< (std::ostream& out, const fftw_complex& x)
  {
    out<<"("<<x.re<<","<<x.im<<")";
    return out;
  }

  
  inline double TConvert (const double& rhs) { return rhs; }
  inline fftw_complex TConvert (const std::complex<double>& rhs) {
    fftw_complex a;
    a.re = rhs.real(); a.im = rhs.imag();
    return a;
  }

  inline void TSetZero (double& x) { x = 0.; }
  inline void TSetZero (std::complex<double>& x) { x = 0.; }
  inline void TSetZero (fftw_complex& x) { x.re = 0.; x.im = 0.; }

  template<class T1, class T2>
  inline void TCopy(T1& lhs, const T2& rhs) { lhs = rhs; }
    
  template<>
  inline void TCopy(Cplx& lhs, const fftw_complex& rhs)
  {
    lhs = Cplx(rhs.re, rhs.im);
  }

  template<>
  inline void TCopy(fftw_complex& lhs, const Cplx& rhs)
  {
    lhs.re = rhs.real();
    lhs.im = rhs.imag();
  }

  /**********************************************************************
   * template<class VectorA, class VectorB> ExtractFastPotential --
   **********************************************************************/
  template<class VectorA, class VectorB>
  void extractFastPotential(const size_t *dim2,
			    // 2 means double in all directions
			    VectorA& src,
			    VectorB& obj)
  {
    
    const size_t nx2 = dim2[0];
    const size_t ny2 = dim2[1];
    const size_t nz2 = dim2[2];

    const size_t nx = nx2 / 2;
    const size_t ny = ny2 / 2;
    const size_t nz = nz2 / 2;

    for(size_t ii=0; ii<nx; ++ii){
      for(size_t jj=0; jj<ny2; ++jj){
	for(size_t kk=0; kk<nz2; ++kk){
	  size_t index2 = kk + nz2*(jj + ny2*ii);
	  if(jj < ny && kk < nz){
	    size_t index = kk + nz*(jj + ny*ii);
	    TCopy(obj[index], src[index2]);
	  }
	}
      }
    }
  }

  /**********************************************************************
   * template<class Vector, class CalcKernel>
   * SetupKernelTopelitz --
   **********************************************************************/
  template<class Vector, class CalcKernel>
  void setupKernelTopelitz (const size_t *dim2,
			    const double *pace,
			    Vector& kernel,
			    const CalcKernel& f)
  {
    const size_t nx2 = dim2[0];
    const size_t ny2 = dim2[1];
    const size_t nz2 = dim2[2];

    const double dx = pace[0];
    const double dy = pace[1];
    const double dz = pace[2];
    
    //    const int length1D2 = nx2 * ny2 * nz2;
    
    for(size_t ii=0; ii<nx2; ++ii){
      for(size_t jj=0; jj<ny2; ++jj){
	for(size_t kk=0; kk<nz2; ++kk){
	  int index2 = kk + nz2*(jj + ny2*ii);
	  TCopy( kernel[index2], 
		 f(circulantIndex(ii, nx2) * dx,
		   circulantIndex(jj, ny2) * dy,
		   circulantIndex(kk, nz2) * dz) );
	}
      }
    }
  }



} // namespace x3dconv

#endif

















