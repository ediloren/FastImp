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

  const static char cvsid[] = "$Id: gmres.h,v 1.10 2003/02/05 21:14:34 zhzhu Exp $";

  ==========================================================================
*/


#ifndef _GMRES_H_
#define _GMRES_H_

#include <complex>
#include "cmat.h"
#include "vec.h"

/**********************************************************************
 * outputResid --
 **********************************************************************/
static void
outputResid (
	     const int outIt,
	     const int inIt,
	     const double resid)
{
#ifdef PRINT_GMRES_RESID
  if (outIt == 0) {
    std::cout << std::endl;
  }
  std::cout << inIt << "\t" 
	    << outIt << "th iter, residual = " << resid << std::endl;

#else
  const int printPoint = 10; 
  if (outIt == 0) {
    std::cout << std::endl;
  } else if ((outIt != 0) && (inIt == 0)) {    
    std::cout << "\touter resid = " << resid << std::endl;
  } else {    
    //    if (inIt == 1)  std::cout << std::endl;
    std::cout << outIt << " ";
    if (outIt % printPoint == 0) {
      std::cout << " inner resid = " << resid << std::endl;
    }
    flush(std::cout);
  }
#endif
}

/**********************************************************************
 * givens --
 **********************************************************************/
template <class T>
static void
givens (
	const T a, 
	const T b, 
	T& c, 
	T& s)
{
  if (b == 0.) {
    c = 1.;
    s = 0.;
  } else {
    double r1 = abs(a);
    double r2 = abs(b);
    double length = sqrt(r1*r1 + r2*r2);
    c = a / length;
    s = -b / length; 
  }
}

/**********************************************************************
 * givensRotate --
 **********************************************************************/
template <class T>
static void
givensRotate (
	      const T c, 
	      const T s, 
	      T& a, 
	      T& b)
{
  T a1 = a;
  T b1 = b;
  a = a1 * conj(c) - b1 * conj(s);
  b = c * b1 + s * a1;
}

/**********************************************************************
 * givensRotate --
 **********************************************************************/
static void
givensRotate (
	      const double c, 
	      const double s, 
	      double& a, 
	      double& b)
{
  double a1 = a;
  double b1 = b;
  a = a1 * c - b1 * s;
  b = c * b1 + s * a1;
}

/**********************************************************************
 * gmres --
 * 
 * On return flag: 
 * true -- converged, false -- diverged
 **********************************************************************/
template <class Matrix, class Preconditioner, class T>
bool
gmres (
       Matrix &A,
       TNT::Vector<T> &x0, // initial guess and solution.
       const TNT::Vector<T> &b,
       Preconditioner &M,
       int maxiters = 300,
       int restart = 20,
       double tol = 1e-6)
{
  TNT::Vector<T> c(restart + 1);
  TNT::Vector<T> s(restart + 1);
  TNT::Vector<T> g(restart + 1);
  TNT::Vector<T> y(restart + 1);
  TNT::Matrix<T> H(restart+1, restart+1);

  int vectorLength = b.size();

  TNT::Vector<T> AP(vectorLength);
  TNT::Vector<T> P(vectorLength);
  std::vector<TNT::Vector<T> > bv(restart+1, TNT::Vector<T>(vectorLength));
  TNT::Vector<T> xLast(x0);
  
  if (M.isLeftPreconditioner()) {
    // M^(-1) * ( A * x0 = b )
    A.matMultVec(P, x0); // P = A * x0;
    P = b - P;
    M.solve(P, AP); // AP = precondition(P);
  } else {
    // x0 = M^(-1) * y, where y is the working unknown vector in the 
    // following iteration. Hence r = b - A * M^(-1) * y = b - A*x0
    A.matMultVec(AP, x0); // AP = A * x0;
    AP = b - AP;
  } 

  double outer_rnorm = two_norm(AP); 
  double bnorm = outer_rnorm; // assuming the initial guess is zero

  int outIt = 0;
  outputResid(outIt, 0, outer_rnorm/bnorm);

  while (outIt < maxiters && outer_rnorm/bnorm > tol) {

    g[0] = outer_rnorm;
    AP /= outer_rnorm;
    bv[0] = AP;

    int inIt = 0;
    double inner_rnorm;
    do {
      if (M.isLeftPreconditioner()) {
	// M^(-1) * A * bv[i]
	A.matMultVec(P, bv[inIt]);
	M.solve(P, AP);
      } else if (M.isRightPreconditioner()) {
	// A * M^(-1) * bv[i]
	M.solve(bv[inIt], P);
	A.matMultVec(AP, P);
      } else {
	A.matMultVec(AP, bv[inIt]);
      }

      for (int j = 0; j <= inIt; j++) {
	H[inIt][j] = inner_prod(AP, bv[j]);
	AP -= bv[j] * H[inIt][j];
      }
      H[inIt][inIt+1] = two_norm(AP);   
      bv[inIt+1] = AP / H[inIt][inIt+1];

      for (int i = 0; i < inIt; i++) {
	givensRotate(c[i], s[i], H[inIt][i], H[inIt][i+1]);
      }
      givens(H[inIt][inIt], H[inIt][inIt+1], c[inIt], s[inIt]);
      givensRotate(c[inIt], s[inIt], H[inIt][inIt], H[inIt][inIt+1]);
      g[inIt+1] = 0.; // this should be zero
      givensRotate(c[inIt], s[inIt], g[inIt], g[inIt+1]);
      inner_rnorm = abs(g[inIt+1]);

      inIt++;
      outIt++;

      outputResid(outIt, inIt, inner_rnorm/bnorm);

    } while (inIt < restart && inner_rnorm/bnorm > tol);

    // it is incremented by one before exit. Adjustment is needed here
    inIt --;
    
    // Compute solution, note, H is H[col][row]. 
    // backsolve for y: H*y = g
    for (int k = 0; k <= inIt; ++k) {
      y[k] = g[k];
    }
    for(int i = inIt; i >= 0; i--) {
      y[i] /= H[i][i];
      for(int j = i-1; j >= 0; j--) {
	y[j] -= H[i][j] * y[i]; 
      }
    }
  
    // x = x0 + bv * y
    for(int j=0; j <= inIt; j++) {
      x0 +=  bv[j] * y[j];
    }

    if (M.isLeftPreconditioner()) {
      // M^(-1) * ( A * x0 = b )
      A.matMultVec(P, x0); // P = A * x0;
      P = b - P;
      M.solve(P, AP); // AP = precondition(P);
    } else if (M.isRightPreconditioner()) {
      // A * M^(-1) * x0 = b
      M.solve(x0, P); // P = precondition(x0);
      A.matMultVec(AP, P); // AP = A * P;
      AP = b - AP;
    } else {
      A.matMultVec(AP, x0); // AP = A * x0;
      AP = b - AP;
    }

    outer_rnorm = two_norm(AP); 
    outputResid(outIt, 0, outer_rnorm/bnorm);
  }

  std::cout << "\n\t Number of Gmres iterations := " << outIt << std::endl;
  std::cout << "\t Final residual := " << outer_rnorm/bnorm << std::endl;

  if (M.isRightPreconditioner()) {
    // M^(-1) * x0 = real solution
    M.solve(x0, x0);
  }

  if (outer_rnorm/bnorm < tol) {
    return true;
  } else {
    return false;
  }

}


#endif

    
