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

  Compressed col scheme is used in this class to store the sparse matrix entries
  So this class could be used to store rectangular matrix as well.

  Resources:

  See also:

  const static char cvsid[] = "$Id: spColMat.h,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _SP_COL_MAT_H_
#define _SP_COL_MAT_H_

#include <vector>
#include <stdexcept>
#include <float.h>
#include "x3dconv.h"
#include "spVec.h"

namespace pfft {

  // Enrico
  template <class T>
  class SpColMat;

  // Enrico
  template<class T>
  std::ostream& operator<< (std::ostream& out,
                          const SpColMat<T>& a);

  template <class T>
  class SpColMat {

  public:
    typedef SpVec<T> SpCol;

    friend std::ostream& operator << <> (std::ostream& os,
					 const SpColMat<T>& sv);     

    SpColMat() : numRow_(0), numNonZero_(0) {};
    SpColMat(size_t numCol, size_t numRow=0);
    const SpCol& operator[] (size_t col) const { return mat[col]; }  
    SpCol& operator[] (size_t col) { return mat[col]; }  
    const size_t numCol(void) const { return mat.size(); }
    const size_t numRow(void) const { return numRow_; }
    size_t memoryEstimate (size_t, size_t) const;
    size_t numNonZero(void) { if (!numNonZero_) countNumNonZero(); return numNonZero_; }
    void insertElement(const size_t row, const size_t col, const T& value) {
      mat[col].push_back(SpVecElement<T>(row, value)); }

  private:
    std::vector<SpCol> mat;
    size_t numRow_;  //This means the numRow if the matrix is full
    size_t numNonZero_;

    void countNumNonZero(void);
  };

  
  /**********************************************************************
   * SpColMat --
   **********************************************************************/
  template <class T>
  SpColMat<T>::SpColMat (
			 size_t numCol, 
			 size_t numRow)
    : numRow_(numRow), numNonZero_(0) 
  {
    mat.resize(numCol);
    if (! numRow_)
      numRow_ = numCol;    
  }

  /**********************************************************************
   * operator <<  --
   **********************************************************************/
  template<class T>
  std::ostream& operator<< (std::ostream& out,
			    const SpColMat<T>& a)
  {
    
    size_t ii = 0;
    for(; ii < a.numCol(); ++ii) {
      out<<" col # "<<ii<<std::endl;
      out<<a.mat[ii]<<std::endl;
    }
    out<<std::endl;
    return out;
  }

  /**********************************************************************
   * fullVec1 = SpColMat * fullVec2 --
   **********************************************************************/
  template <class T, class FullVector>
  FullVector
  operator* (
	     const SpColMat<T>& mat, 
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
	      const SpColMat<T>& mat, 
	      const FullVector2& fv)
  {
    for (size_t row = 0; row < mat.numRow(); row++) {
      TSetZero(ans[row]);
    }
    for (size_t col = 0; col < mat.numCol(); col++) {
      for (size_t row = 0; row < mat[col].size(); row++) {
	ans[mat[col].index(row)] += mat[col].value(row) * fv[col];
      }
    }    
  }

  /**********************************************************************
   * fullVec1 += SpColMat * fullVec2 --
   **********************************************************************/
  template <class T1, class T2, class FullVector1, class FullVector2>
  void
  matMultVecAdd (
		 FullVector1& ans,
		 const SpColMat<T1>& mat, 
		 const FullVector2& fv,
		 const T2& fvSample)
  {
    for (size_t col = 0; col < mat.numCol(); col++) {
      const SpVec<T1>& matCol = mat[col];
      T2 fve = fv[col];
      for (size_t row = 0; row < mat[col].size(); row++) {
	ans[matCol[row].index()] += matCol[row].value() * fve;
      }
    }    
  }

  template <class T, class FullVector1, class FullVector2>
  void
  matMultVecAdd (
		 FullVector1& ans,
		 const SpColMat<T>& mat, 
		 const FullVector2& fv)
  {
    for (size_t col = 0; col < mat.numCol(); col++) {
      const SpVec<T>& matCol = mat[col];
      for (size_t row = 0; row < matCol.size(); row++) {
	ans[matCol[row].index()] += matCol[row].value() * fv[col];
      }
    }    
  }

  /**********************************************************************
   * fullVec1 = fullVec2 * SpColMat --
   **********************************************************************/
  template <class T, class FullVector>
  FullVector
  operator* (
	     const FullVector& fv,
	     const SpColMat<T>& mat)
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
	      const SpColMat<T>& mat)
  {
    for (size_t col = 0; col < mat.numCol(); col++) {
      ans[col] = 0.;
      for (size_t row = 0; row < mat[col].size(); row++) {
	ans[col] += mat[col].value(row) * fv[ mat[col].index(row) ];
      }
    }    
  }

  /**********************************************************************
   * fullVec1 += fullVec2 * SpColMat --
   **********************************************************************/
  template <class T, class FullVector1, class FullVector2>
  void
  vecMultMatAdd (
		 FullVector1& ans,
		 const FullVector2& fv,
		 const SpColMat<T>& mat)
  {
    for (size_t col = 0; col < mat.numCol(); col++) {
      for (size_t row = 0; row < mat[col].size(); row++) {
	ans[col] += mat[col].value(row) * fv[ mat[col].index(row) ];
      }
    }    
  }

  /**********************************************************************
   * SpColMat * SpVec --
   **********************************************************************/
  template <class T, class FullVector>
  void
  matMultVec (
	      FullVector& ans,
	      const SpColMat<T>& mat, 
	      const SpVec<T>& sv)
  {
    for (size_t row = 0; row < mat.numRow(); row++) {
      ans[row] = 0.;
    }

    for (size_t i = 0; i < sv.size(); i++) {
      size_t col = sv.index(i);
      for (size_t row = 0; row < mat[col].size(); row++) {
	ans[mat[col].index(row)] += mat[col].value(row) * sv.value(i);
      }
    }    
  }

  /**********************************************************************
   * SpVec * SpColMat --
   **********************************************************************/
  template <class T, class FullVector>
  void
  vecMultMat (
	      FullVector& ans,
	      const SpVec<T>& sv,
	      const SpColMat<T>& mat)
  {
    //    ans = FullVector(mat.numCol());
    for (size_t col = 0; col < mat.numCol(); col++) {
      dotProd(ans[col], mat[col], sv);
    }
  }

  /**********************************************************************
   * countNumNonZero --
   **********************************************************************/
  template <class T>
  void
  SpColMat<T>::countNumNonZero(
			       void)
  {
    numNonZero_ = 0;
    for (size_t col = 0; col < numCol(); col++) {
      numNonZero_ += mat[col].size();
    } 
  }

  /**********************************************************************
   * memoryEstimate --
   **********************************************************************/
  template <class T>
  size_t
  SpColMat<T>::memoryEstimate (
			       size_t numEntry,
			       size_t numCol) const
  {
    size_t entry = (sizeof(size_t) + sizeof(T)) * numEntry;
    size_t col = numCol * sizeof(SpCol);
    return entry + col;
  }


} //namespace pfft 

#endif




