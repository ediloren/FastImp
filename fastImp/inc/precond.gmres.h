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

  const static char cvsid[] = "$Id: precond.gmres.h,v 1.1.1.1 2002/04/15 21:41:01 bsong Exp $";

  ==========================================================================
*/


#ifndef _GMRES_H_
#define _GMRES_H_

#include <complex>
#include "cmat.h"
#include "vec.h"

typedef TNT::Vector<std::complex<double> > ComplexVector;
typedef TNT::Matrix<std::complex<double> >::size_type Index;

/**********************************************************************
 * outputResid --
 **********************************************************************/
void
outputResid (
	     const Index outIt,
	     const Index inIt,
	     const double resid)
{
#ifdef PRINT_GMRES_RESID
  if (outIt == 0) {
    std::cout << std::endl;
  }
  std::cout << inIt << "\t" 
	    << outIt << "th iter, residual = " << resid << std::endl;

#else
  if (outIt == 0) {
    std::cout << std::endl;
  } else {    
    if (inIt == 1)  std::cout << std::endl;
    std::cout << outIt << " ";
  }
#endif
}

/**********************************************************************
 * givens --
 **********************************************************************/
void
givens (
	const std::complex<double> a, 
	const std::complex<double> b, 
	std::complex<double>& c, 
	std::complex<double>& s)
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
void
givensRotate (
	      const std::complex<double> c, 
	      const std::complex<double> s, 
	      std::complex<double>& a, 
	      std::complex<double>& b)
{
  std::complex<double> a1 = a;
  std::complex<double> b1 = b;
  a = a1 * conj(c) - b1 * conj(s);
  b = c * b1 + s * a1;
}

/**********************************************************************
 * printVec --
 **********************************************************************/
template <class Vector>
void
printVector (
	     const char* fileName,
	     const Vector& v)
{
  std::ofstream fout(fileName);
  for (size_t i=0; i < v.size(); i++) {
    fout << v[i] << endl;
  }
  fout.close();
}

void
printVector (
	     const char* fileName,
	     const std::complex<double>* v,
	     const size_t length)
{
  std::ofstream fout(fileName);
  for (size_t i=0; i < length; i++) {
    fout << v[i] << endl;
  }
  fout.close();
}

/**********************************************************************
 * gmres --
 * 
 * On return flag: 
 * true -- converged, false -- diverged
 **********************************************************************/
template <class Matrix, class Vector, class VectorRHS, class Preconditioner>
bool
gmres (
       Matrix &A,
       Vector &x0, // initial guess and solution.
       const VectorRHS &b,
       Preconditioner &M,
       int maxiters = 300,
       int restart = 20,
       double tol = 1e-6)
{
  ComplexVector c(restart + 1);
  ComplexVector s(restart + 1);
  ComplexVector g(restart + 1);
  ComplexVector y(restart + 1);
  TNT::Matrix<std::complex<double> > H(restart+1, restart+1);

  Index vectorLength = b.size();

  ComplexVector AP(vectorLength);
  ComplexVector P(vectorLength);
  std::vector<ComplexVector> bv(restart+1, ComplexVector(vectorLength));
  
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

  double rnorm = two_norm(AP); 
  double bnorm = rnorm; // assuming the initial guess is zero

  Index outIt = 0;
  outputResid(outIt, 0, rnorm/bnorm);

  while (outIt < maxiters && rnorm/bnorm > tol) {

    g[0] = rnorm;
    AP /= rnorm;
    bv[0] = AP;

#ifdef DEBUG_GMRES
    printVector("q0.tmp", bv[0]);
#endif

    Index inIt = 0;
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

      for (Index j = 0; j <= inIt; j++) {
	H[inIt][j] = inner_prod(AP, bv[j]);
	AP -= bv[j] * H[inIt][j];
      }
      H[inIt][inIt+1] = two_norm(AP);   
      bv[inIt+1] = AP / H[inIt][inIt+1];

#ifdef DEBUG_GMRES
      if (inIt == 0) {
	printVector("q1.tmp", bv[inIt+1]);
      } else if (inIt == 1) {
	printVector("q2.tmp", bv[inIt+1]);
      } else if (inIt == 2) {
	printVector("q3.tmp", bv[inIt+1]);
      } else if (inIt == 3) {
	printVector("q4.tmp", bv[inIt+1]);
      }
#endif

      for (Index i = 0; i < inIt; i++) {
	givensRotate(c[i], s[i], H[inIt][i], H[inIt][i+1]);
      }
      givens(H[inIt][inIt], H[inIt][inIt+1], c[inIt], s[inIt]);
      givensRotate(c[inIt], s[inIt], H[inIt][inIt], H[inIt][inIt+1]);
      g[inIt+1] = 0.; // this should be zero
      givensRotate(c[inIt], s[inIt], g[inIt], g[inIt+1]);
      rnorm = abs(g[inIt+1]);

      inIt++;
      outIt++;

      outputResid(outIt, inIt, rnorm/bnorm);

    } while (inIt < restart && rnorm/bnorm > tol);

    // it is incremented by one before exit. Adjustment is needed here
    inIt --;
    
    // Compute solution, note, H is H[col][row]. 
    // backsolve for y: H*y = g
    for (Index k = 0; k <= inIt; ++k) {
      y[k] = g[k];
    }
    for(Index i = inIt; i >= 0; i--) {
      y[i] /= H[i][i];
      for(Index j = i-1; j >= 0; j--) {
	y[j] -= H[i][j] * y[i]; 
      }
    }
  
    // x = x0 + bv * y
    for(Index j=0; j <= inIt; j++) {
      x0 += y[j] * bv[j];
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

    rnorm = two_norm(AP); 
  }

  std::cout << "\n\t Number of Gmres iterations := " << outIt << std::endl;
  std::cout << "\t Final residual := " << rnorm/bnorm << std::endl;

  if (M.isRightPreconditioner()) {
    // M^(-1) * x0 = real solution
    M.solve(x0, x0);
  }

#ifdef DEBUG_GMRES
  printVector("x0.tmp", x0);
  printVector("y.tmp", y);
#endif

  if (rnorm/bnorm < tol) {
    return true;
  } else {
    return false;
  }

}

#endif

    
