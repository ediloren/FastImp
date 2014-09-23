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

  const static char cvsid[] = "$Id: spVec.h,v 1.4 2003/02/11 02:56:04 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _SP_VEC_H_
#define _SP_VEC_H_

#include <vector>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <float.h>
#include <cfloat> // for DBL_MIN
#include "spVecElement.h"
#include "x3dconv.h"

namespace pfft {

  template <class T> 

  class SpVec {

  public:
    typedef typename std::vector<SpVecElement<T> >::const_iterator SpVecConstIterator;
    typedef typename std::vector<SpVecElement<T> >::iterator SpVecIterator;

    friend std::ostream& operator << <> (std::ostream& os, const SpVec<T>& sv);     
	        
    SpVec(void) {};
    SpVec(size_t numEle) { v.resize(numEle); }
    SpVec<T>& operator*= (const T& s);
    void push_back(const SpVecElement<T> & ele) { v.push_back(ele); }
    void reserve(size_t num) { v.reserve(num); }
    const size_t capacity(void) const { return v.capacity(); }
    void resize(size_t num) { v.resize(num); }
    const size_t size(void) const { return v.size(); }
    void clear(void) { return v.clear(); }
    const SpVecConstIterator begin(void) const { return v.begin(); }
    const SpVecConstIterator end(void) const { return v.end(); }
    const SpVecElement<T> & operator[] (size_t i) const { return v[i]; }
    SpVecElement<T> & operator[] (size_t i) { return v[i]; }
    const size_t largestNonZeroIndex (void) const;
    const size_t index(size_t i) const { return v[i].index(); }
    const T value(size_t i) const { return v[i].value(); }

  private:
    std::vector<SpVecElement<T> > v;
  };

  /**********************************************************************
   * largestNonZeroIndex --
   **********************************************************************/
  template <class T> 
  const size_t
  SpVec<T>::largestNonZeroIndex (
				 void) const
  {
    size_t largest = 0;
    for (SpVecConstIterator ele = v.begin(); ele != v.end(); ++ele) {
      largest = max(largest, ele->index());
    }
    return largest;
  }

  /**********************************************************************
   * operator*= --
   **********************************************************************/
  template <class T> 
  SpVec<T>&
  SpVec<T>::operator*= (
			const T& s)
  {
    for (SpVecIterator ele = v.begin(); ele != v.end(); ++ele) {
      ele->value() *= s;
    }
    return *this;
  }

  /**********************************************************************
   * operator* --
   **********************************************************************/
  template <class T> 
  SpVec<T>
  operator* (
	     SpVec<T>& sv,
	     const T s)
  {
    SpVec<T> ans = sv;
    ans *= s;
    return ans;
  }

  /**********************************************************************
   * dotProd --
   * SpVec * SpVec --
   **********************************************************************/
  template <class T1, class T2, class T3> 
  void
  dotProd (
	   T1& ans,
	   const SpVec<T2>& sv1, 
	   const SpVec<T3>& sv2)
  {
    std::vector<T2> fv(sv1.largestNonZeroIndex());
    scatter(sv1, fv);
    ans = 0;
    for (size_t i = 0; i < sv2.size(); i++) { 
      ans += sv2.value(i) * fv[sv2.index(i)];
    }
  }

  /**********************************************************************
   * dotProd --
   * SpVec * FullVector --
   **********************************************************************/
  template <class T1, class T2, class FullVector> 
  void
  dotProd (
	   T1& ans,
	   const SpVec<T2>& sv, 
	   const FullVector& fv)
  {
    ans = 0;
    for (size_t i = 0; i < sv.size(); i++) {
      ans += sv.value(i) * fv[ sv.index(i) ];
    }
  }

  /**********************************************************************
   * dotProd --
   * FullVector * SpVec --
   **********************************************************************/
  template <class T1, class T2, class FullVector> 
  void
  dotProd (
	   T1& ans,
	   const FullVector& fv,
	   const SpVec<T2>& sv)
  {
    dotProd(ans, sv, fv);
  }

  /**********************************************************************
   * compactDotProd --
   * SpVec * FullVector --
   **********************************************************************/
  template <class T1, class T2, class FullVector> 
  void
  compactDotProd (
		  T1& ans,
		  const SpVec<T2>& sv, 
		  const FullVector& fv)
  {
    ans = 0;
    for (size_t i = 0; i < sv.size(); i++) {
      ans += sv.value(i) * fv[i];
    }
  }

  /**********************************************************************
   * compactDotProd --
   * FullVector * SpVec --
   **********************************************************************/
  template <class T1, class T2, class FullVector> 
  void
  compactDotProd (
		  T1& ans,
		  const FullVector& fv,
		  const SpVec<T2>& sv)
  {
    compactDotProd(ans, sv, fv);
  }

  /**********************************************************************
   * selectDotProd --
   * SpVec * FullVector --
   **********************************************************************/
  template <class T1, class T2, class FullVector> 
  void
  selectDotProd (
		 T1& ans,
		 const SpVec<T2>& sv, 
		 const FullVector& fv, 
		 const std::vector<int>& indexSelect)
  {
    ans = 0;
    for (size_t i = 0; i < sv.size(); i++) {
      ans += sv.value(i) * fv[indexSelect[i]];
    }
  }

  /**********************************************************************
   * matMultVec --
   * Assumption:
   * Memory of the full vector has been allocated before calling this 
   * function
   **********************************************************************/
  template <class T, class FullVector, class DenseMat>
  void
  matMultVec (
	      FullVector& fv,
	      const DenseMat& mat,
	      const SpVec<T>& sv)
  {
    for (int row = 0; row < fv.size(); row++) {
      fv[row] = 0.;
      //      const T* rowi = mat[row];
      for (size_t col = 0; col < sv.size(); col++) {
	fv[row] += mat[row][col] * sv.value(col);
	//	fv[row] += rowi[col] * sv.value(col);
      }
    }
  }

  /**********************************************************************
   * scatter --
   **********************************************************************/
  template <class T, class FullVector> 
  void 
  scatter (
	   const SpVec<T> sv,
	   FullVector& fv)
  {
    for (typename SpVec<T>::SpVecConstIterator ele = sv.begin(); 
	 ele != sv.end(); ++ele) {
      fv[ele->index()] = ele->value();
    }    
  }

  /**********************************************************************
   * output --
   **********************************************************************/
  template <class T> 
  std::ostream& 
  operator << (
	       std::ostream& os, 
	       const SpVec<T>& sv) 
  {
    for (typename SpVec<T>::SpVecConstIterator ele = sv.v.begin();
	 ele != sv.v.end(); ++ele) {
      os << *ele;
    }
    os << std::endl;
    return os;
  }

  /**********************************************************************
   * copyFullVectorToSpVec --
   **********************************************************************/
  template <class T, class FullVectorIndexIterator, class FullVectorValueIterator> 
  void 
  copyFullVectorToSpVec (
			 const FullVectorIndexIterator indexBegin,
			 const FullVectorIndexIterator indexEnd, 
			 const FullVectorValueIterator valueBegin, 
			 SpVec<T>& spVec)
  {
    FullVectorValueIterator value = valueBegin;
    for (FullVectorIndexIterator index = indexBegin; index != indexEnd; 
	 ++index, ++value) {
      if (abs(*value) <= DBL_MIN) continue; // skip zeros
      spVec.push_back(SpVecElement<T>(*index, TConvert(*value)));
    }
  }

  /**********************************************************************
   * copySpVecToFullVector --
   **********************************************************************/
  template <class SpVecIterator, class FullVectorIterator> 
  void 
  copySpVecToFullVector (
			 const SpVecIterator svBegin,
			 const SpVecIterator svEnd,
			 FullVectorIterator fvBegin)
  {
    FullVectorIterator fv = fvBegin; 
    for (SpVecIterator sv = svBegin; sv != svEnd; ++sv, ++fv) {
      *fv = sv->value();
    }    
  }

} // namespace pfft

#endif
