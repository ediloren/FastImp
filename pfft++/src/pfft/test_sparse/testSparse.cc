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

  const static char cvsid[] = "$Id: testSparse.cc,v 1.1 2002/03/20 03:48:48 bsong Exp $";

  ==========================================================================
*/

#include "spVec.h"
#include "spVecElement.h"
#include "spRowMat.h"
#include "spColMat.h"
#include <iostream>
#include <complex>

using namespace std;
using namespace pfft;

void testSpColMat(void);
void testSpColMat(void);

/**********************************************************************
 * testSpRowMat --
 **********************************************************************/
void testSpRowMat (
		   void)
{
  const size_t numRow = 4;
  SpRowMat<complex<double> > A(numRow);
  SpVec<complex<double> > spRow;

  // fill the sparse matrix: method I:
  // form a compressed row then put it in.
  spRow.push_back(SpVecElement<complex<double> >(0, 1.));
  spRow.push_back(SpVecElement<complex<double> >(1, 2.));
  spRow.push_back(SpVecElement<complex<double> >(3, 3.));
  A[0] = spRow;

  spRow.clear();
  spRow.clear();
  spRow.push_back(SpVecElement<complex<double> >(0, 10.));
  spRow.push_back(SpVecElement<complex<double> >(2, 14.));
  A[1] = spRow;

  spRow.clear();
  spRow.push_back(SpVecElement<complex<double> >(2, .2));
  spRow.push_back(SpVecElement<complex<double> >(3, .5));
  A[2] = spRow;

  spRow.clear();
  spRow.push_back(SpVecElement<complex<double> >(0, 100));
  spRow.push_back(SpVecElement<complex<double> >(3, 200));
  A[3] = spRow;

  // form the matrix: method II
  // randomly insert element.
  A.insertElement(0, 5, complex<double>(1, -1));
  A.insertElement(3, 9, complex<double>(1, -1));

  // output the sparse matrix
  cout << "<--------- A: a sparse matrix based on compressed row ----------->" <<endl<<endl;
  cout << A;
  cout << endl << endl;

  // matrix By a full vector
  vector<complex<double> > ans (A.numRow());
  vector<double> v (A.numCol());
  cout << "<---------------  v: a full vector  -----------------> "<<endl<<endl;
  for (size_t ii = 0; ii<A.numCol(); ii++) {
    v[ii] = double(ii) / double(A.numCol());
    cout <<"v["<<ii<<"] = "<<v[ii] <<endl;
  }

  matMultVec(ans, A, v);
  cout << "<---------------  ans = A x v = ? --------------------> " <<endl<<endl; 
  for (size_t ii=0; ii<ans.size(); ii++) {
    cout << " ans["<<ii<<"]  = "<< ans[ii] << endl;
  }
  cout << endl;
  
  // v2 += ans * mat
  vector<complex<double> > vv (A.numCol());
  for (size_t ii=0; ii<vv.size(); ii++) {
    vv[ii] = v[ii];
  }

  vecMultMatAdd(vv, ans, A);
  cout << "<-----------------  vv(=v) += ans * mat ?  ---------------->" <<endl<<endl; 
  for (size_t ii=0; ii<ans.size(); ii++) {
    cout << " ans["<<ii<<"]  = "<< ans[ii] << endl;
  }
  cout << endl;
  
}

/**********************************************************************
 * testSpColMat --
 **********************************************************************/
void testSpColMat (
		   void)
{
  const size_t numCol = 4;
  SpColMat<complex<double> > A(numCol);
  SpVec<complex<double> > spCol;

  // fill the sparse matrix: method I:
  // form a compressed row then put it in.
  spCol.push_back(SpVecElement<complex<double> >(0, 1.));
  spCol.push_back(SpVecElement<complex<double> >(1, 2.));
  spCol.push_back(SpVecElement<complex<double> >(3, 3.));
  A[0] = spCol;

  spCol.clear();
  spCol.clear();
  spCol.push_back(SpVecElement<complex<double> >(0, 10.));
  spCol.push_back(SpVecElement<complex<double> >(2, 14.));
  A[1] = spCol;

  spCol.clear();
  spCol.push_back(SpVecElement<complex<double> >(2, .2));
  spCol.push_back(SpVecElement<complex<double> >(3, .5));
  A[2] = spCol;

  spCol.clear();
  spCol.push_back(SpVecElement<complex<double> >(0, 100));
  spCol.push_back(SpVecElement<complex<double> >(3, 200));
  A[3] = spCol;

  // form the matrix: method II
  // randomly insert element.
  A.insertElement(0, 5, complex<double>(1., -1.));
  A.insertElement(3, 9, complex<double>(1., -1.));

  // output the sparse matrix
  cout << "<--------- A: a sparse matrix based on compressed col ----------->" <<endl<<endl;
  cout << A;
  cout << endl << endl;

  // matrix By a full vector
  vector<complex<double> > ans (A.numRow());
  vector<double> v (A.numCol());
  cout << "<---------------  v: a full vector  -----------------> "<<endl<<endl;
  for (size_t ii = 0; ii<A.numCol(); ii++) {
    v[ii] = double(ii) / double(A.numCol());
    cout <<"v["<<ii<<"] = "<<v[ii] <<endl;
  }

  matMultVecAdd(ans, A, v);
  cout << "<---------------  ans += A x v = ? --------------------> " <<endl<<endl; 
  for (size_t ii=0; ii<ans.size(); ii++) {
    cout << " ans["<<ii<<"]  = "<< ans[ii] << endl;
  }
  cout << endl;
  
  // vv = ans * mat
  vector<complex<double> > vv (A.numCol());
  for (size_t ii=0; ii<vv.size(); ii++) {
    vv[ii] = v[ii];
  }

  vecMultMatAdd(vv, ans, A);
  cout << "<-----------------  vv(=v) = ans * mat ?  ---------------->" <<endl<<endl; 
  for (size_t ii=0; ii<ans.size(); ii++) {
    cout << " ans["<<ii<<"]  = "<< ans[ii] << endl;
  }
  cout << endl;
  
}

/**********************************************************************
 * main --
 **********************************************************************/
int main (
	  void)
{
  testSpRowMat();
  testSpColMat();
  return 0;
}
 


