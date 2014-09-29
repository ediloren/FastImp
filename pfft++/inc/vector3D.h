/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Zhenhai Zhu, Jing Wang
  
  Description:

  Resources:

  See also:

  const static char cvsid[] = "$Id: vector3D.h,v 1.6 2003/02/11 02:56:42 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _VECTOR_3D_H_
#define _VECTOR_3D_H_

#include <complex>
#include <iostream>

namespace pfft {

// Enrico, adding the template for the friend function, otherwise gcc will throw an error
// 'declaration of 'operator>>' as non-function.
// Removing the <> in the program produces some useful information - and then adding
// function prototypes makes the error go away
// (of course, since the friend function does use the class in which it is declared friend,
// you must add the template for the class as well, in this case 'vector3D'
// Explanation in the gcc ticket:
// Wolfgang Bangerth 2004-08-23 14:07:29 UTC
// I think the error message is actually quite clear: you try to make
// something a friend which isn't known: for a friend declaration, if
// name lookup doesn't find an existing name, and if it is a declaration
// of either a function or a function template, then this injects the
// name into the surrounding namespace. However, this doesn't apply for
// a function template specialization, as in your case, if the actual
// template declaration isn't available. Thus the compiler complains
// that you are trying to declare something you can't declare (a non-function,
// non-template). It should be obvious that the compiler can't do anything
// about it unless you provide a declaration of the general template.
template <class T>
class vector3D;

template <class T>
std::istream&
operator >> (
    std::istream& is,
    vector3D<T>& vec);

  template <class T>

  class vector3D {

    friend std::istream& operator>><> (std::istream&, vector3D<T>&);

  public:
    vector3D(void) {}
    vector3D(const T x, const T y, const T z) : x_(x), y_(y), z_(z) {}
    
    vector3D& operator= (const vector3D& rhs) {
      if(this != &rhs) {
	x_ = rhs.x_; y_ = rhs.y_; z_ = rhs.z_;
      } 
      return *this;
    }
	
    const T x(void) const { return x_; }
    const T y(void) const { return y_; }
    const T z(void) const { return z_; }
    
    T& x(void) { return x_; }
    T& y(void) { return y_; }
    T& z(void) { return z_; }

    vector3D<T> & operator += (const vector3D<T>& vec);
    vector3D<T> & operator -= (const vector3D<T>& vec);
    vector3D<T> & operator *= (const T& s);
    vector3D<T> & operator /= (const T& s);
    void normalize (void) ;
    void transferLocalToGlobalCoord (
				     const vector3D<T>& origin,
				     const vector3D<T>& X,
				     const vector3D<T>& Y,
				     const vector3D<T>& Z);
    void transferGlobalToLocalCoord (
				     const vector3D<T>& origin,
				     const vector3D<T>& X,
				     const vector3D<T>& Y,
				     const vector3D<T>& Z);
    template <class T2> void rotateCoord (
					  const vector3D<T2>& X, 
					  const vector3D<T2>& Y, 
					  const vector3D<T2>& Z);
    void shiftOrigin (const vector3D<T>& newOrigin);

  private:
    T x_, y_, z_;

  };

  typedef vector3D<double> point3D;

  /**********************************************************************
   * += --
   **********************************************************************/
  template <class T>
  vector3D<T> & 
  vector3D<T>::operator += (
			    const vector3D<T>& vec) 
  { 
    x_ += vec.x(); 
    y_ += vec.y(); 
    z_ += vec.z(); 
    return *this;
  }  

  /**********************************************************************
   * -= --
   **********************************************************************/
  template <class T>
  vector3D<T> & 
  vector3D<T>::operator -= (
			    const vector3D<T>& vec) 
  { 
    x_ -= vec.x(); 
    y_ -= vec.y(); 
    z_ -= vec.z(); 
    return *this;
  }  

  /**********************************************************************
   * *= --
   **********************************************************************/
  template <class T>
  vector3D<T> & 
  vector3D<T>::operator *= (
			    const T& s) 
  { 
    x_ *= s; 
    y_ *= s; 
    z_ *= s; 
    return *this;
  }  

  /**********************************************************************
   * /= --
   **********************************************************************/
  template <class T>
  vector3D<T> & 
  vector3D<T>::operator /= (
			    const T& s) 
  { 
    x_ /= s; 
    y_ /= s; 
    z_ /= s; 
    return *this;
  }  

  /**********************************************************************
   * normalize --
   **********************************************************************/
  template <class T>
  void 
  vector3D<T>::normalize (
			  void) 
  { 
    double len = length(*this); 
    if (len != 0.) {
      x_ /= static_cast<T>(len); 
      y_ /= static_cast<T>(len); 
      z_ /= static_cast<T>(len); 
    }
  }

  /**********************************************************************
   * normalize --
   **********************************************************************/
  template <class T>
  vector3D<T>
  normalize (
	     const vector3D<T>& v) 
  { 
    double len = length(v); 
    if (len != 0.) {
      T x = v.x() / static_cast<T>(len); 
      T y = v.y() / static_cast<T>(len); 
      T z = v.z() / static_cast<T>(len); 
      return vector3D<T>(x, y, z);
    }

    return v;
  }

  /**********************************************************************
   * operator* --
   **********************************************************************/
  template<class T1, class T2>
  vector3D<T1> 
  operator * (
	      const vector3D<T1>& vec, 
	      const T2& s)
  {
    vector3D<T1> vt = vec; 
    vt *= s; 
    return vt;
  }

  template<class T1, class T2>
  vector3D<T1> 
  operator * (
	      const T2& s, 
	      const vector3D<T1>& vec) 
  {
    vector3D<T1> vt = vec; 
    vt *= s; 
    return vt;
  }

  /**********************************************************************
   * operator/ --
   **********************************************************************/
  template <class T1, class T2>
  vector3D<T1>
  operator / (
	      const vector3D<T1>& vec, 
	      const T2& s) 
  {
    vector3D<T1> vt = vec; 
    vt /= s; 
    return vt;
  }

  /**********************************************************************
   * operator+ --
   **********************************************************************/
  template <class T>
  vector3D<T>
  operator + (
	      const vector3D<T>& v1, 
	      const vector3D<T>& v2)
  {
    vector3D<T> vt = v1; 
    vt += v2; 
    return vt;
  }

  /**********************************************************************
   * operator- --
   **********************************************************************/
  template <class T>
  vector3D<T>
  operator - (
	      const vector3D<T>& v1, 
	      const vector3D<T>& v2)
  {
    vector3D<T> vt = v1; 
    vt -= v2; 
    return vt;
  }

  /**********************************************************************
   * dot prod --
   **********************************************************************/
  template <class T>
  T 
  operator * (
	      const vector3D<T>& v1, 
	      const vector3D<T>& v2) 
  {
    return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
  }

  template <class T1, class T2>
  T1 
  operator * (
	      const vector3D<T1>& v1, 
	      const vector3D<T2>& v2) 
  {
    return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
  }

  template <class T1, class T2, class T3>
  void
  dotProd (
	   T3& ans,
	   const vector3D<T1>& v1, 
	   const vector3D<T2>& v2) 
  {
    ans = v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
  }

  /**********************************************************************
   * cross product --
   **********************************************************************/
  template <class T>
  vector3D<T> 
  crossProd (
	     const vector3D<T>& v1, 
	     const vector3D<T>& v2)
  {
    T x = v1.y()*v2.z() - v1.z()*v2.y();
    T y = v1.z()*v2.x() - v1.x()*v2.z();
    T z = v1.x()*v2.y() - v1.y()*v2.x();
    return vector3D<T>(x, y, z);
  }

  template <class T1, class T2, class T3>
  void
  crossProd (
	     vector3D<T3>& ans,
	     const vector3D<T1>& v1, 
	     const vector3D<T2>& v2)
  {
    T3 x = v1.y()*v2.z() - v1.z()*v2.y();
    T3 y = v1.z()*v2.x() - v1.x()*v2.z();
    T3 z = v1.x()*v2.y() - v1.y()*v2.x();
    ans = vector3D<T3>(x, y, z);
  }

  /**********************************************************************
   * comparison --
   **********************************************************************/
  template <class T>
  bool
  operator == (
	       const vector3D<T>& v1, 
	       const vector3D<T>& v2)
  {
    return ( (v1.x()==v2.x()) && (v1.y()==v2.y()) && (v1.z()==v2.z()) );
  }

  template <class T>
  bool
  operator != (
	       const vector3D<T>& v1, 
	       const vector3D<T>& v2)
  {
    return ( (v1.x()!=v2.x()) || (v1.y()!=v2.y()) || (v1.z()!=v2.z()) );
  }

  template <class T>
  bool
  operator < (
	       const vector3D<T>& v1, 
	       const vector3D<T>& v2)
  {
    if (v1.x() < v2.x()) {
      return true;
    } else if (v1.x() > v2.x()) {
      return false;
    } else {
      if (v1.y() < v2.y()) {
	return true;
      } else if (v1.y() > v2.y()) {
	return false;
      } else {
	if (v1.z() < v2.z()) {
	  return true;
	} else {
	  return false;
	}
      }
    }
  }

  template <class T>
  bool
  operator > (
	       const vector3D<T>& v1, 
	       const vector3D<T>& v2)
  {
    if (v1.x() > v2.x()) {
      return true;
    } else if (v1.x() < v2.x()) {
      return false;
    } else {
      if (v1.y() > v2.y()) {
	return true;
      } else if (v1.y() < v2.y()) {
	return false;
      } else {
	if (v1.z() > v2.z()) {
	  return true;
	} else {
	  return false;
	}
      }
    }
  }

  /**********************************************************************
   * input --
   **********************************************************************/
  template <class T>
  std::istream& 
  operator >> (
	       std::istream& is, 
	       vector3D<T>& vec)
  {
    T x, y, z;
    is >> x >> y >> z;
    vec = vector3D<T>(x, y, z);
    return is;
  }

  /**********************************************************************
   * output --
   **********************************************************************/
  template <class T>
  std::ostream& 
  operator << (
	       std::ostream& os, 
	       const vector3D<T>& vec) 
  {
    os <<  "(" << vec.x() 
       << ", " << vec.y()
       << ", " << vec.z() 
       << ")";
      // << ", length = "
      //<< length(vec)
    //       << std::endl;
    return os;
  }

  /**********************************************************************
   * length --
   **********************************************************************/
  // for non-complex 
  template <class T>
  T
  length(
	 const vector3D<T>& vec)
  {
    return sqrt(vec * vec);
  }

  // for complex vector only
  template <class T>
  T
  length(
	 const vector3D<std::complex<T> >& vec)
  {
    T len = abs(vec.x()) * abs(vec.x()) +
      abs(vec.y()) * abs(vec.y()) +
      abs(vec.z()) * abs(vec.z());
    return sqrt(len);
  }

  /**********************************************************************
   * conj --
   **********************************************************************/
  template <class T>
  vector3D<std::complex<T> >
  conj (
	const vector3D<std::complex<T> >& vec)
  { 
    return vector3D<std::complex<T> >(conj(vec.x()), 
				      conj(vec.y()), conj(vec.z())); 
  }

  /**********************************************************************
   * transferLocalToGlobalCoord --
   * Vectors origin, X, Y and Z are in global coordinate system.
   * Vector vec is in the local coordinate system consists of origin
   * X, Y and Z. 
   * This function find the global coordinates of vec
   **********************************************************************/
  template <class T1>
  vector3D<T1>
  transferLocalToGlobalCoord (
			      const vector3D<T1>& vec,
			      const vector3D<T1>& origin,
			      const vector3D<T1>& X,
			      const vector3D<T1>& Y,
			      const vector3D<T1>& Z)
  {
    vector3D<T1> ans(vec);
    ans.transferLocalToGlobalCoord(origin, X, Y, Z);
    return ans;
  }

  /**********************************************************************
   * transferLocalToGlobalCoord --
   * Vectors origin, X, Y and Z are in global coordinate system.
   * Vector *this is in the local coordinate system consists of origin
   * X, Y and Z. 
   * This function find the global coordinates of *this
   **********************************************************************/
  template <class T1>
  void
  vector3D<T1>::transferLocalToGlobalCoord (
					    const vector3D<T1>& origin,
					    const vector3D<T1>& X,
					    const vector3D<T1>& Y,
					    const vector3D<T1>& Z)
  {
    vector3D<T1> rotation(X.x(), Y.x(), Z.x());
    T1 x = (*this) * rotation + origin.x();

    rotation = vector3D<T1>(X.y(), Y.y(), Z.y());
    T1 y = (*this) * rotation + origin.y();

    rotation = vector3D<T1>(X.z(), Y.z(), Z.z());
    T1 z = (*this) * rotation + origin.z();

    x_ = x;
    y_ = y;
    z_ = z;
  }

  /**********************************************************************
   * transferGlobalToLocalCoord --
   * Vectors vec, origin, X, Y and Z are in global coordinate system.
   * This function find the local coordinates of vector vec.
   * This local coordinate system consists of origin, X, Y and Z. 
   **********************************************************************/
  template <class T1>
  vector3D<T1>
  transferGlobalToLocalCoord (
			      const vector3D<T1>& vec,
			      const vector3D<T1>& origin,
			      const vector3D<T1>& X,
			      const vector3D<T1>& Y,
			      const vector3D<T1>& Z)
  {
    vector3D<T1> ans(vec);
    ans.transferGlobalToLocalCoord(origin, X, Y, Z);
    return ans;
  }

  /**********************************************************************
   * transferGlobalToLocalCoord --
   * Vectors *this, origin, X, Y and Z are in global coordinate system.
   * This function find the local coordinates of vector *this.
   * This local coordinate system consists of origin, X, Y and Z. 
   **********************************************************************/
  template <class T1>
  void
  vector3D<T1>::transferGlobalToLocalCoord (
					    const vector3D<T1>& origin,
					    const vector3D<T1>& X,
					    const vector3D<T1>& Y,
					    const vector3D<T1>& Z)
  {
    T1 x = (*this - origin) * X;
    T1 y = (*this - origin) * Y;
    T1 z = (*this - origin) * Z;

    x_ = x;
    y_ = y;
    z_ = z;
  }

  /**********************************************************************
   * rotateCoord --
   * Vector vec is in the local coordinate system which is consists of 
   * X, Y and Z.
   * X, Y and Z are vectors described in the global coordinate system. 
   * This function transfer vector vec to the global coordinate system.
   * Since two systems share the same origin, only rotation is performed. 
   **********************************************************************/
  template <class T1, class T2>
  vector3D<T1>
  rotateCoord (
	       const vector3D<T1>& vec,
	       const vector3D<T2>& X,
	       const vector3D<T2>& Y,
	       const vector3D<T2>& Z)
  {
    vector3D<T1> ans(vec);
    ans.rotateCoord(X, Y, Z);
    return ans;
  }

  /**********************************************************************
   * rotateCoord --
   * Vector *this is in the local coordinate system which is consists of 
   * X, Y and Z.
   * X, Y and Z are vectors described in the global coordinate system. 
   * This function transfer vector *this' to the global coordinate system.
   * Since two systems share the same origin, only rotation is performed. 
   * 
   * Note:
   * The difference between this function and transferLocalToGlobalcoord() 
   * is not just the origin is 0. More important thing is the value type of
   * x_, y_ and z_ could be different from that of X, Y, Z. 
   **********************************************************************/
  template <class T1>
  template <class T2>
  void
  vector3D<T1>::rotateCoord (
			     const vector3D<T2>& X,
			     const vector3D<T2>& Y,
			     const vector3D<T2>& Z)
  {
    vector3D<T2> rotation(X.x(), Y.x(), Z.x());
    T1 x = (*this) * rotation;

    rotation = vector3D<T2>(X.y(), Y.y(), Z.y());
    T1 y = (*this) * rotation;

    rotation = vector3D<T2>(X.z(), Y.z(), Z.z());
    T1 z = (*this) * rotation;

    x_ = x;
    y_ = y;
    z_ = z;
  }

  /**********************************************************************
   * shiftOrigin --
   **********************************************************************/
  template <class T>
  void
  vector3D<T>::shiftOrigin (
			    const vector3D<T>& newOrigin)
  {
    x_ -= newOrigin.x();
    y_ -= newOrigin.y();
    z_ -= newOrigin.z();
  }

} // namespace pfft 

#endif

