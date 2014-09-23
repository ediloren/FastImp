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

  const static char cvsid[] = "$Id: spVecElement.h,v 1.3 2002/07/18 14:32:19 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _SP_VEC_ELEMENT_H_
#define _SP_VEC_ELEMENT_H_

namespace pfft {

  template <class T> 
  class SpVecElement {
  public:
    SpVecElement(void) : index_(0), value_(0) {};
    SpVecElement(size_t index, T value) : index_(index), value_(value) {};

    friend bool operator< <> (const SpVecElement<T>&, const SpVecElement<T>&);
    SpVecElement<T>& operator+= (const T s) { value_ += s; return *this; }
    SpVecElement<T>& operator-= (const T s) { value_ -= s; return *this; }
    SpVecElement<T>& operator*= (const T s) { value_ *= s; return *this; }
    SpVecElement<T>& operator/= (const T s) { value_ /= s; return *this; }

    const T value(void) const { return value_; }
    T& value(void) { return value_; }
    size_t index(void) const { return index_; }

  private:
    size_t index_;
    T value_;
    
  };

  /**********************************************************************
   * output --
   **********************************************************************/
  template <class T> 
  std::ostream& 
  operator << (
	       std::ostream& os, 
	       const SpVecElement<T>& ele) 
  {
    os << "(" << ele.index() << "," << ele.value() << ")";
    return os;
  }

  /**********************************************************************
   * operator< --
   **********************************************************************/
  template<class T>
  bool operator< (
		  const SpVecElement<T>& ele1, 
		  const SpVecElement<T>& ele2) 
  {
    return (ele1.index_ < ele2.index_);
  }

} // namespace pfft


#endif

