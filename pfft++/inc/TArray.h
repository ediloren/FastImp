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

  const static char cvsid[] = "$Id: TArray.h,v 1.2 2002/08/01 14:02:38 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _TARRAY_H_
#define _TARRAY_H_

#include <iostream>
#include "myArray.h"
#include "x3dconv.h"

namespace pfft{
    
  template<class T>
  class TArray<T, GridDataTag> : public MyArray<T>{
  public: 
  
    TArray(): MyArray<T>() {}

    TArray(size_t n1, size_t n2, size_t n3)
      : MyArray<T>(n1*n2*n3*4), nx_(n1), ny_(n2), nz_(n3) {}
    
    TArray(size_t n1, size_t n2, size_t n3, const T& t)
      : MyArray<T>(n1*n2*n3*4, t), nx_(n1), ny_(n2), nz_(n3) {}

    void resize(size_t n1, size_t n2, size_t n3) {
      this->MyArray<T>::resize(n1 * n2 * n3 * 4);
      nx_ = n1; ny_ = n2; nz_ = n3;
    }
      
    template<class T1, class T2, class T3>
    friend class Fast3DConv;
  
  private:
    size_t nx_, ny_, nz_;
  };

  
  template<>
  class TArray<Real, KernelTag> : public MyArray<Real> {
  public:
    TArray() : MyArray<Real>(), speq_(static_cast<Real*>(0)) {}

    TArray(size_t n1, size_t n2, size_t n3)
      : MyArray<Real>(n1*n2*n3*8), nx_(n1), ny_(n2), nz_(n3) {
      speq_ = static_cast<Real*>(::operator new( 2*(n1*2)*(n2*2) * sizeof(Real), nothrow));
      assert(speq_ != 0);
    }

    TArray(size_t n1, size_t n2, size_t n3, const Real& t)
      : MyArray<Real>(n1*n2*n3*8, t), nx_(n1), ny_(n2), nz_(n3) {
      speq_ = static_cast<Real*>(::operator new( 2*(n1*2)*(n2*2) * sizeof(Real), nothrow));
      assert(speq_ != 0);
      for(size_t ii=0; ii<2*(n1*2)*(n2*2); ++ii){
	speq_[ii] = t;
      }
    }

    void resize(size_t n1, size_t n2, size_t n3) {
      this->MyArray<Real>::resize(n1 * n2 * n3 * 8);
      speq_ = static_cast<Real*>(::operator new( 2*(n1*2)*(n2*2) * sizeof(Real), nothrow));
      assert(speq_ != 0);
      nx_ = n1, ny_ = n2, nz_ = n3;
    }

    ~TArray(){::operator delete(speq_);}
 
    inline size_t SizeOfSpeq() const { return 2 * (nx_*2) * (ny_*2); }
    
    template<class T1, class T2, class T3>
    friend class Fast3DConv;
    
    template<class KernelContainer>
    friend void rfft3DGeneral(const size_t*, KernelContainer&, fftw_direction, bool);

  private:
    size_t nx_, ny_, nz_;
    Real *speq_;

  };
    
  template<>
  class TArray<fftw_complex, KernelTag> : public MyArray<fftw_complex> {
  public:
    TArray() : MyArray<fftw_complex>() {}
    TArray(size_t n1, size_t n2, size_t n3)
      : MyArray<fftw_complex>(n1*n2*n3*8), nx_(n1), ny_(n2), nz_(n3) {}

    TArray(size_t n1, size_t n2, size_t n3, const fftw_complex& z)
      : MyArray<fftw_complex>(n1*n2*n3*8, z), nx_(n1), ny_(n2), nz_(n3) {}

    void resize(size_t n1, size_t n2, size_t n3) {
      this->MyArray<fftw_complex>::resize(n1 * n2 * n3 * 8);
      nx_ = n1, ny_ = n2, nz_ = n3;
    }

    template<class T1, class T2, class T3>
    friend class Fast3DConv;

  private:
    size_t nx_, ny_, nz_;
  };

}//namespace pfft

#endif









