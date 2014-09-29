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

  numerical vector (0-based [i] AND 1-based (i) indexing )

  Resources:

  See also:

  const static char cvsid[] = "$Id: vec.h,v 1.8 2003/03/03 19:51:44 zhzhu Exp $";

  ==========================================================================
*/



#ifndef __VEC_H_
#define __VEC_H_

#include "subscript.h"
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <complex>
#include <cfloat> // for DBL_MAX
#include <vector> // for function unique()

namespace TNT {

  template <class T>
  class Vector {

  public:
    typedef Subscript   size_type;
    typedef         T   value_type;
    typedef         T   element_type;
    typedef         T*  pointer;
    typedef         T*  iterator;
    typedef         T&  reference;
    typedef const   T*  const_iterator;
    typedef const   T&  const_reference;

    Subscript lbound() const { return 1;}

    Vector() : v_(0), vm1_(0), n_(0)  {};
    Vector(const Vector<T> &A) : v_(0), vm1_(0), n_(0) {
      initialize(A.n_);
      copy(A.v_);
    }
    template <class STL_Container>
    Vector(const STL_Container &A) : v_(0), vm1_(0), n_(0) {
      initialize(A.size());
      for (Subscript i=0; i< n_; i++) v_[i] = A[i];
    }
    /* This interfere with the above constructor. Compilor could not 
       differentiate Subscript and STL_Container.
       So I use two specilizations, shown below.
    Vector(Subscript N) :  v_(0), vm1_(0), n_(0) {
      initialize(N);
    }
    */
    Vector(int N) :  v_(0), vm1_(0), n_(0) {
      initialize(N);
    }
    Vector(size_t N) :  v_(0), vm1_(0), n_(0) {
      initialize(N);
    }
    Vector(Subscript N, const T& value) :  v_(0), vm1_(0), n_(0) {
      initialize(N);
      set(value);
    }
    Vector(Subscript N, const T* v) :  v_(0), vm1_(0), n_(0) {
      initialize(N);
      copy(v);
    }
    //    Vector(Subscript N, char *s) :  v_(0), vm1_(0), n_(0){
    //      initialize(N);
    //      std::istrstream ins(s);
    //      for (Subscript i=0; i<N; i++) ins >> v_[i];
    //    }
    ~Vector() { destroy(); }

    iterator begin() { return v_;}
    iterator end()   { return v_ + n_; }
    const iterator begin() const { return v_;}
    const iterator end() const  { return v_ + n_; }

    Vector<T>& newsize(Subscript N) {
      if (n_ == N) return *this;
      destroy();
      initialize(N);
      return *this;
    }

    Vector<T>& operator=(const Vector<T> &A) {
      if (v_ == A.v_) return *this;
      if (n_ == A.n_) {
	copy(A.v_);
      } else {
	destroy();
	initialize(A.n_);
	copy(A.v_);
      }
      return *this;
    }
        
    Vector<T>& operator=(const T& scalar) { 
      set(scalar);  
      return *this;
    }

    Subscript dim() const { return n_; }
    Subscript size() const { return  n_; }
    reference operator[](Subscript i) { return v_[i]; }
    const_reference operator[](Subscript i) const { return v_[i]; }
    reference first(void) { return v_[0]; }
    reference last(void) { return v_[n_-1]; }

    Vector<T>& operator += (const Vector<T>& A ) {
      for (Subscript i=0; i < n_; i++) v_[i] += A.v_[i];
      return *this;
    }

    Vector<T>& operator -= (const Vector<T>& A ) {
      for (Subscript i=0; i < n_; i++) v_[i] -= A.v_[i];
      return *this;
    }

    template <class T1>
    Vector<T>& operator *= (const T1& scalar ) {
      for (Subscript i=0; i < n_; i++) v_[i] *= scalar;
      return *this;
    }

    template <class T1>
    Vector<T>& operator /= (const T1& scalar ) {
      for (Subscript i=0; i < n_; i++) v_[i] /= scalar;
      return *this;
    }
 
    const T min(void) const {
      T minimum = DBL_MAX;
      for (Subscript i=0; i < n_; i++) minimum = minimum < v_[i] ? minimum : v_[i];
      return minimum;
    }

    const T max(void) const {
      T maximum = DBL_MIN;
      for (Subscript i=0; i < n_; i++) maximum = maximum > v_[i] ? maximum : v_[i];
      return maximum;
    }

    void unique (void) {
      if (n_ == 0) return;
      std::vector<T> tmp;
      tmp.push_back(v_[0]);
      for (Subscript i=1; i < n_; i++) {
	if (v_[i] != v_[i-1]) tmp.push_back(v_[i]);
      }
      destroy();
      initialize(tmp.size());
      for (Subscript i=0; i < n_; i++) v_[i] = tmp[i];
    }

    void outputToFile(const char* fileName) const {
      std::ofstream fout(fileName);
      for (Subscript i=0; i < n_; i++) {
	fout << v_[i] << std::endl;
      }
      fout.close();
    }

  protected:
    T* v_;                  
    T* vm1_;        // pointer adjustment for optimzied 1-offset indexing
    Subscript n_;

    void initialize (Subscript N) {
      // adjust pointers so that they are 1-offset:
      // v_[] is the internal contiguous array, it is still 0-offset
      assert(v_ == NULL);
      v_ = new T[N];
      assert(v_  != NULL);
      vm1_ = v_-1;
      n_ = N;
    }
    void copy(const T* v) { for (Subscript i=0; i< n_; i++) v_[i] = v[i]; }
    void set(const T& val) { for (Subscript i=0; i< n_; i++) v_[i] = val; }
    void destroy(void) {     
      if (v_ == NULL) return ;
      delete [] (v_);     
      v_ = NULL;
      vm1_ = NULL;
    }

  };

  /**********************************************************************
   * operator<< --
   **********************************************************************/
  template <class T>
  std::ostream& 
  operator<< (
	      std::ostream &s, 
	      const Vector<T> &A)
  {
    Subscript N=A.dim();
  
    //    s <<  "length = " << N << endl;
    for (Subscript i=0; i<N; i++)
      // Enrico, modified to prevent gcc error
      s   << A[i] << std::endl;
      //s   << A[i] << endl;
  
    return s;
  }

  /**********************************************************************
   * operator>> --
   **********************************************************************/
  template <class T>
  std::istream & 
  operator >> (
	       std::istream &s, 
	       Vector<T> &A)
  {
    Subscript N;
  
    s >> N;
    if ( !(N == A.size() )) {
      A.newsize(N);
    }

    for (Subscript i=0; i<N; i++)
      s >>  A[i];

    return s;
  }

  /**********************************************************************
   * operator + --
   **********************************************************************/
  template <class T1, class T2>
  Vector<T1> 
  operator + (
	      const Vector<T1> &A, 
	      const Vector<T2> &B)
  {
    assert(A.dim()==B.dim());
    Vector<T1> tmp(A);
    tmp += B;
    return tmp;
  }

  /**********************************************************************
   * operator - --
   **********************************************************************/
  template <class T1, class T2>
  Vector<T1> 
  operator - (
	      const Vector<T1> &A, 
	      const Vector<T2> &B)
  {
    assert(A.dim()==B.dim());
    Vector<T1> tmp(A);
    tmp -= B;
    return tmp;
  }

  /**********************************************************************
   * operator * --
   **********************************************************************/
  template <class T1, class T2>
  Vector<T1> 
  operator * (
	      const Vector<T1> &A, 
	      const Vector<T2> &B)
  {
    Subscript N = A.dim();

    assert(N==B.dim());

    Vector<T1> tmp(N);
    for (Subscript i=0; i<N; i++)
      tmp[i] = A[i] * B[i];

    return tmp;
  }

  template <class T1, class T2>
  Vector<T1> 
  operator * (
	      const Vector<T1> &A, 
	      const T2& scalar )
  {
    Vector<T1> tmp(A);
    tmp *= scalar;
    return tmp;
  }

  // The following function is ambiguously overloaded as operator 
  // in spColMat. Hence not used
  /*
  template <class T1,  class T2>
  Vector<T1> 
  operator * (
	      const T2& scalar,
	      const Vector<T1> &A)
  {
    Vector<T1> tmp(A);
    tmp *= scalar;
    return tmp;
  }
  */

  /**********************************************************************
   * operator / --
   **********************************************************************/
  template <class T1, class T2>
  Vector<T1> 
  operator / (
	      const Vector<T1> &A, 
	      const T2& scalar )
  {
    Vector<T1> tmp(A);
    tmp /= scalar;
    return tmp;
  }

  /**********************************************************************
   * inner_prod --
   **********************************************************************/
  template <class T>
  T 
  inner_prod ( 
	      const Vector<T> &A, 
	      const Vector<T> &B)
  {
    Subscript N = A.dim();
    assert(N == B.dim());

    T sum = 0;
    for (Subscript i=0; i<N; i++)
      sum += A[i] * conj(B[i]);

    return sum;
  }
    
  template <class T1, class T2, class T3>
  void
  inner_prod ( 
	      T3& res,
	      const Vector<T1> &A, 
	      const Vector<T2> &B)
  {
    Subscript N = A.dim();
    assert(N == B.dim());

    res = 0.;
    for (Subscript i=0; i<N; i++)
      res += A[i] * B[i];

    return res;
  }

  /*
  std::complex<double> 
  inner_prod ( 
	      const Vector<std::complex<double> > &A, 
	      const Vector<std::complex<double> > &B)
  {
    Subscript N = A.dim();
    assert(N == B.dim());

    std::complex<double> sum = 0;
    for (Subscript i=0; i<N; i++)
      sum += A[i] * conj(B[i]);

    return sum;
  }
  */
 
  /**********************************************************************
   * two_norm --
   **********************************************************************/
  template <class T>
  double
  two_norm ( 
	    const Vector<T>& A)
  {
    double norm = 0.;
    for (Subscript i=0; i < A.size(); i++) {
      norm += pow(abs(A[i]), 2);
    }
    return sqrt(norm);
  }

  /**********************************************************************
   * infinite_norm --
   **********************************************************************/
  template <class T>
  double
  infinite_norm ( 
		 const Vector<T>& A)
  {
    double norm = 0.;
    for (Subscript i=0; i < A.size(); i++) {
      norm = max(norm, abs(A[i]));
    }
    return norm;
  }

  /**********************************************************************
   * separate_two_norm --
   **********************************************************************/
  template <class T>
  void
  separate_two_norm ( 
		     const Vector<std::complex<T> >& A,
		     double& rnorm,
		     double& inorm)
  {
    rnorm = 0.;
    inorm = 0.;
    for (Subscript i=0; i < A.size(); i++) {
      rnorm += pow(real(A[i]), 2);
      inorm += pow(imag(A[i]), 2);
    }
    rnorm = sqrt(rnorm);
    inorm = sqrt(inorm);
  }

  /**********************************************************************
   * separate_infinite_norm --
   **********************************************************************/
  template <class T>
  void
  separate_infinite_norm ( 
			  const Vector<std::complex<T> >& A,
			  double& rnorm,
			  double& inorm)
  {
    rnorm = 0.;
    inorm = 0.;
    for (Subscript i=0; i < A.size(); i++) {
      rnorm = max(rnorm, abs(real(A[i])));
      inorm = max(inorm, abs(imag(A[i])));
    }
  }

  /**********************************************************************
   * sum --
   **********************************************************************/
  template <class T>
  T
  sum ( 
       const Vector<T>& A)
  {
    T sum = static_cast<T>(0);
    for (Subscript i=0; i < A.size(); i++) {
      sum += A[i];
    }
    return sum;
  }

}   /* namespace TNT */

#endif
// VEC_H
