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

  ==========================================================================
*/

#include "vector3D.h"
#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>

using namespace pfft;

static void
testDouble1 (
	     void);

static void
testDouble (
	    void);

static void
testComplex (
	    void);

static void
testIO (
	void);

static void
testComparision (
		 void);

template <class T>
static void
printVec (
	  const string& name,
	  vector3D<T>& vec);

/**********************************************************************
 * main --
 **********************************************************************/
int
main (
      void)
{
  //  testDouble();
  //  testDouble1();
  testComplex();
  //  testIO ();
  //  testComparision ();

}

/**********************************************************************
 * testComparision --
 **********************************************************************/
static void
testComparision (
		 void)
{
  vector3D<double> v1(2.1, 2.1, 0.1);
  vector3D<double> v2(0., 0., 1.);
  if (v1==v2)
    cout << "true" << endl;
  else if (v1!=v2)
    cout << "false" << endl;

}

/**********************************************************************
 * testIO --
 **********************************************************************/
static void
testIO (
	void)
{
  ifstream dataSource("data2.dat");
  vector< vector3D<double> > vec;
  vector3D<double> tmp;

  if (!dataSource) {
    
  } else {
    while (dataSource >> tmp) {
      vec.push_back(tmp);
    }
  }

  for (vector<vector3D<double> >::const_iterator it=vec.begin(); it!=vec.end(); it++) {
    cout << *it;
  }
}


/**********************************************************************
 * testDouble1 --
 **********************************************************************/
static void
testDouble1 (
	     void)
{
  vector3D<double> v1(2.1, 2.1, 0.1);
  vector3D<double> v2(0., 0., 1.);
  vector3D<double> v3(1., 1., 2.);
  vector3D<double> v4;
  v4 = (v1+v2+v3)/3.;
  printVec<double>("v4 = ", v4);

  v4 = v1*2;
  printVec<double>("v4 = ", v4);

  v4 = 2*v1;
  printVec<double>("v4 = ", v4);

  v4 = v1;
  v4 *= 2;
  printVec<double>("v4 = ", v4);
}


/**********************************************************************
 * testDouble --
 **********************************************************************/
static void
testDouble (
	    void)
{
  vector3D<double> v1(1., 2., 3.);
  vector3D<double> v2(v1);
  vector3D<double> v3;
  double x;

  printVec<double>("v1 = ", v1);

  x = v1*v2;
  cout << x << endl;

  v2 *= 2.;
  printVec<double>("v2 = ", v2);

  v3 = v1 * 2.;
  printVec<double>("v3 = ", v3);

  v3 = 2. * v1;
  printVec<double>("v3 = ", v3);

  v3 += v1;
  printVec<double>("v3 = ", v3);

  v3 = v3 - v1;
  printVec<double>("v3 = ", v3);

  v3 = v3 + v1;
  printVec<double>("v3 = ", v3);

  v2.normalize();
  printVec<double>("v2 = ", v2);

  vector3D<double> v4(0., 0., 0.);
  v4.normalize();
  printVec<double>("v4 = ", v4);

  vector3D<double> v5(1., 0., 0.);
  v3 = crossProd(v1, v5);
  printVec<double>("v3 = ", v3);

  v5.x() = 2.;
  v3 = crossProd(v1, v5);
  printVec<double>("v3 = ", v3);
}

/**********************************************************************
 * testComplex --
 **********************************************************************/
static void
testComplex (
	    void)
{
  vector3D<complex<double> > v1(1., 2., 3.);
  vector3D<complex<double> > v2(v1);
  vector3D<complex<double> > v3;
  complex<double>  x;
  complex<double>  con(2.,0.);
  
  printVec<complex<double> >("v1 = ", v1);

  x = v1*v2;
  cout << x << endl;

  v2 *= 2.;
  printVec<complex<double> >("v2 = ", v2);

  v3 = v1 * con;
  printVec<complex<double> >("v3 = ", v3);

  v3 = con * v1;
  printVec<complex<double> >("v3 = ", v3);

  v3 += v1;
  printVec<complex<double> >("v3 = ", v3);

  v3 = v3 - v1;
  printVec<complex<double> >("v3 = ", v3);

  v3 = v3 + v1;
  printVec<complex<double> >("v3 = ", v3);

  v2.normalize();
  printVec<complex<double> >("v2 = ", v2);

  vector3D<complex<double> > v4(0., 0., 0.);
  v4.normalize();
  printVec<complex<double> >("v4 = ", v4);

  vector3D<complex<double> > v5(1., 0., 0.);
  v3 = crossProd(v1, v5);
  printVec<complex<double> >("v3 = ", v3);

  v5 = vector3D<complex<double> >(2., v5.y(), v5.z());
  v3 = crossProd(v1, v5);
  printVec<complex<double> >("v3 = ", v3);

  complex<double> c1(1.,0.), c2(2.,0.), c3(0.,3.);
  vector3D<complex<double> > v10(c1, c2, c3);
  printVec<complex<double> >("v10 = ", v10);
  
  vector3D<complex<double> > v11(conj(v10));
  printVec<complex<double> >("v11 = ", v11);

  vector3D<double> v12(1., 2., 3.);
  printVec<double>("v12 = ", v12);

  complex<float> c4(1.,0.), c5(2.,0.), c6(0.,3.);
  vector3D<complex<float> > v13(c4, c5, c6);
  printVec<complex<float> >("v13 = ", v13);
  cout << v13 << endl;

}

/**********************************************************************
 * printVec --
 **********************************************************************/
template <class T>
static void
printVec (
	  const string& name,
	  vector3D<T>& vec)
{
  cout << name << vec << " " << length(vec) << endl;
}

