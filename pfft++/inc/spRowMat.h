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

  Compressed row scheme is used in this class to store the sparse matrix entries
  So this class could be used to store rectangular matrix as well.

  Resources:

  See also:

  const static char cvsid[] = "$Id: spRowMat.h,v 1.4 2003/07/17 02:46:16 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _SP_ROW_MAT_H_
#define _SP_ROW_MAT_H_

#include <vector>
#include <float.h>
#include <cfloat> // for DBL_EPSILON
#include "spVec.h"

namespace pfft {

  template <class T>
  class SpRowMat {

  public:
    typedef SpVec<T> SpRow;

    friend std::ostream& operator << <> (std::ostream& os,
					 const SpRowMat<T>& sv); 

    SpRowMat() : numCol_(0), numNonZero_(0) {};
    SpRowMat(size_t numRow, size_t numCol=0);
    const SpRow& operator[] (size_t row) const { return mat[row]; }  
    SpRow& operator[] (size_t row) { return mat[row]; }  
    const size_t numRow(void) const { return mat.size(); }
    const size_t numCol(void) const { return numCol_; }
    size_t memoryEstimate (size_t, size_t) const;
    size_t numNonZero(void) { if (!numNonZero_) countNumNonZero(); return numNonZero_; }
    void insertElement(const size_t row, const size_t col, const T& value) {
      mat[row].push_back(SpVecElement<T>(col, value)); }

    typename std::vector<SpRow>::const_iterator begin(void) const { return mat.begin(); }  
    typename std::vector<SpRow>::const_iterator end(void) const { return mat.end(); }  

  private:
    std::vector<SpRow> mat;
    size_t numCol_;  //This means the numCol if the matrix is full
    size_t numNonZero_;

    void countNumNonZero(void);
  };

  /**********************************************************************
   * SpRowMat --
   **********************************************************************/
  template <class T>
  SpRowMat<T>::SpRowMat (
			 size_t numRow, 
			 size_t numCol)
    : numCol_(numCol), numNonZero_(0) 
  {
    mat.resize(numRow);
    if (! numCol_)
      numCol_ = numRow;
  }

  /**********************************************************************
   * fullVec1 = SpColMat * fullVec2 --
   **********************************************************************/
  template <class T, class FullVector>
  FullVector
  operator* (
	     const SpRowMat<T>& mat, 
	     const FullVector& fv)
  {
    FullVector ans(mat.numRow());
    matMultVec(ans, mat, fv);
    return ans;
  }

  template <class T, class FullVector1, class FullVector2>
  void
  matMultVec (
	      FullVector1& ans,
	      const SpRowMat<T>& mat, 
	      const FullVector2& fv)
  {
    for (size_t row = 0; row < mat.numRow(); row++) {
      ans[row] = 0.;
      for (size_t col = 0; col < mat[row].size(); col++) {
	ans[row] += mat[row].value(col) * fv[ mat[row].index(col) ];
      }
    }    
  }

  // the result is put into part of a large vector. So the input
  // ans is just a pointer to the starting position.
  template <class T, class FullVector1Ptr, class FullVector2>
  void
  matMultVecPartial (
		     FullVector1Ptr ansPtr,
		     const SpRowMat<T>& mat, 
		     const FullVector2& fv)
  {
    for (size_t row = 0; row < mat.numRow(); row++) {
      ansPtr[row] = 0.;
      for (size_t col = 0; col < mat[row].size(); col++) {
	ansPtr[row] += mat[row].value(col) * fv[ mat[row].index(col) ];
      }
    }    
  }

  /**********************************************************************
   * fullVec1 += SpRowMat * fullVec2 --
   **********************************************************************/
  template <class T, class FullVector1, class FullVector2>
  void
  matMultVecAdd (
		 FullVector1& ans,
		 const SpRowMat<T>& mat, 
		 const FullVector2& fv)
  {
    for (size_t row = 0; row < mat.numRow(); row++) {
      for (size_t col = 0; col < mat[row].size(); col++) {
	ans[row] += mat[row].value(col) * fv[ mat[row].index(col) ];
      }
    }    
  }

  /**********************************************************************
   * fullVec1 = fullVec2 * SpRowMat --
   **********************************************************************/
  template <class T, class FullVector>
  FullVector
  operator* (
	     const FullVector& fv,
	     const SpRowMat<T>& mat)
  {
    FullVector ans(mat.numCol());
    vecMultMat(ans, fv, mat);
    return ans;
  }

  template <class T, class FullVector1, class FullVector2>
  void
  vecMultMat (
	      FullVector1& ans,
	      const FullVector2& fv,
	      const SpRowMat<T>& mat)
  {
    for (size_t i = 0; i < mat.numCol(); i++) {
      ans[i] = 0.;
    }

    for (size_t row = 0; row < mat.numRow(); row++) {
      for (size_t col = 0; col < mat[row].size(); col++) {
	ans[mat[row].index(col)] += mat[row].value(col) * vec[row];
      }
    }    
  }

  /**********************************************************************
   * fullVec1 += fullVec2 * SpRowMat --
   **********************************************************************/
  template <class T, class FullVector1, class FullVector2>
  void
  vecMultMatAdd (
		 FullVector1& ans,
		 const FullVector2& fv,
		 const SpRowMat<T>& mat)
  {
    for (size_t row = 0; row < mat.numRow(); row++) {
      for (size_t col = 0; col < mat[row].size(); col++) {
	ans[mat[row].index(col)] += mat[row].value(col) * fv[row];
      }
    }    
  }

  /**********************************************************************
   * SpRowMat * SpVec --
   **********************************************************************/
  template <class T, class FullVector>
  void
  matMultVec (
	      FullVector& ans,
	      const SpRowMat<T>& mat, 
	      const SpVec<T>& sv)
  {
    //    ans = FullVector(mat.numRow());
    for (size_t row = 0; row < mat.numRow(); row++) {
      dotProd(ans[row], mat[row], sv);
    }
  }

  /**********************************************************************
   * SpVec * SpRowMat --
   **********************************************************************/
  template <class T, class FullVector>
  void
  vecMultMat (
	      FullVector& ans,
	      const SpVec<T>& sv,
	      const SpRowMat<T>& mat)
  {
    for (size_t i = 0; i < mat.numCol(); i++) {
      ans[i] = 0.;
    }
    for (size_t i = 0; i < sv.size(); i++) {
      size_t row = sv.index(i);
      for (size_t col = 0; col < mat[row].size(); col++) {
	ans[mat[row].index(col)] += mat[row].value(col) * sv.value(i);
      }
    }    
  }

  /**********************************************************************
   * countNumNonZero --
   **********************************************************************/
  template <class T>
  void
  SpRowMat<T>::countNumNonZero(
			       void)
  {
    numNonZero_ = 0;
    for (size_t row = 0; row < numRow(); row++) {
      numNonZero_ += mat[row].size();
    } 
  }

  /**********************************************************************
   * memoryEstimate --
   **********************************************************************/
  template <class T>
  size_t
  SpRowMat<T>::memoryEstimate (
			       size_t numEntry,
			       size_t numRow) const
  {
    size_t entry = (sizeof(size_t) + sizeof(T)) * numEntry;
    size_t row = numRow * sizeof(SpRow);
    return entry + row;
  }

  /**********************************************************************
   * operator << --
   **********************************************************************/
  template<class T>
  std::ostream& operator<< (std::ostream& out,
			    const SpRowMat<T>& a)
  {
    size_t ii = 0;
    for(; ii < a.numRow(); ++ii) {
      std::vector<SpVecElement<T> > tmp;
      for(size_t jj=0; jj<a.mat[ii].size(); ++jj) {
	// only print out values bigger than machine precision
	if (abs(a.mat[ii][jj].value()) > DBL_EPSILON) {
	  tmp.push_back(SpVecElement<T>(a.mat[ii][jj].index(), a.mat[ii][jj].value()));
	}
      }
      sort(tmp.begin(), tmp.end());
      out<<" row # "<<ii<<" size: "<<tmp.size()<<std::endl;
      for(size_t jj=0; jj<tmp.size(); ++jj) {
	out<<"  "<<tmp[jj].index()<<", "<<tmp[jj].value()<<"  ";
      } out <<std::endl;
    }
    out<<std::endl;
    return out;
  }

} //namespace pfft 

#endif




