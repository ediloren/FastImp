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

  const static char cvsid[] = "$Id: myArray.h,v 1.2 2002/10/01 21:17:37 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _MY_ARRAY_H_
#define _MY_ARRAY_H_

#include <new>
#include <cassert>
#include <iostream>
#include <iterator>

namespace pfft{
 
// Enrico, prototype declaration
template<class T>
  class MyArray;

// Enrico, prototype declaration
template<class T>
  std::ostream& operator<< (std::ostream& out, const MyArray<T>& x);

  template<class T>
  class MyArray{
  public:
    typedef T* iterator;
    typedef const T* const_iterator;

    MyArray(): N_(0), p_(0) {}
    MyArray(size_t n); // non-initialization
    MyArray(size_t n, const T&); // initialization 
    ~MyArray();

    void resize(size_t n);
    size_t size(void) const { return N_; }

    iterator begin() { return p_; }
    const_iterator begin() const { return p_; }

    iterator end() { return p_ + N_; }
    const_iterator end() const { return p_ + N_; }

    inline T& operator[] (size_t n);
    inline const T& operator[] (size_t n) const;

    // friend classes & functions:
    // add pfft

    // I/O
    friend std::ostream& operator<< <> (std::ostream&, const MyArray<T>&);

  protected:
    size_t N_;
    T *p_;
  };

  template<class T>
  MyArray<T>::MyArray(size_t n) : N_(n)
  {
    if(n == 0){
      p_ = 0;
    } else {  

      // Enrico, corrected to avoid gcc error
      p_ = static_cast<T*>(::operator new(N_ * sizeof(T), std::nothrow));
      //p_ = static_cast<T*>(::operator new(N_ * sizeof(T), nothrow));
      assert(p_ != 0);

    }
  }

  template<class T>
  MyArray<T>::MyArray(size_t n, const T& t): N_(n)
  {
    if(n == 0){
      p_ = 0;
    } else {
      // Enrico, corrected to avoid gcc error
      p_ = static_cast<T*>(::operator new(N_ * sizeof(T), std::nothrow));
      //p_ = static_cast<T*>(::operator new(N_ * sizeof(T), nothrow));
      assert(p_ != 0);
      for(size_t ii=0; ii<N_; ++ii){
	p_[ii] = t;
      }
    }
  }
  
  template<class T>
  void MyArray<T>::resize(size_t n)
  {
    assert(n > 0);
    if(n > N_){
      ::operator delete(p_);
     // Enrico, corrected to avoid gcc error
      p_ = static_cast<T*>(::operator new(n * sizeof(T), std::nothrow));
      //p_ = static_cast<T*>(::operator new(n * sizeof(T), nothrow));
      assert(p_ != 0);
    }
    N_ = n;
      
  }

  /*
  template<class T>
  template<class Vin<Elt>>
  MyArray<T>::MyArray(const Vin& v)
  {
    assert(!v.empty());

    if(this != &v){
      if(v.size() > N_){
	::operator delete(p_);
	N_ = v.size();
        p_ = static_cast<T*>(::operator new(N_ * sizeof(T))

  */

  template<class T>
  T& MyArray<T>::operator[] (size_t n)
  {
    assert(n <= N_);
    return p_[n];
  }

  template<class T>
  const T& MyArray<T>::operator[] (size_t n) const
  {
    assert(n <= N_);
    return p_[n];
  }

  /**********************************************************************
   * ~MyArray --
   **********************************************************************/
  template<class T>
  MyArray<T>::~MyArray()
  {
    if(p_ != 0){
      ::operator delete(p_);
    }
  }

  
  /**********************************************************************
   * operator<< --
   **********************************************************************/
  template<class T>
  std::ostream& operator<< (std::ostream& out, const MyArray<T>& x)
  {
    assert(x.p_ != 0);
    copy(x.p_, x.p_ + x.N_, std::ostream_iterator<T>(out, "\t"));
    return out;
  }

} // namespace pfft

#endif 

  









