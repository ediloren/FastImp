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

  This module extract impedance from the solution of the linear equations

  Resources:

  See also:

  const static char cvsid[] = "$Id: extractZ.cc,v 1.21 2003/07/16 02:07:57 zhzhu Exp $";

  ==========================================================================
*/

#include <fstream> // for outputZ
#include <complex>
#include <vector>
#include "extractZ.h"
#include "gmres.h"
#include "dense.h" // for pseudo_inverse()

using namespace std;
using namespace surf;
using namespace mesh;
using namespace pfft;

/**********************************************************************
 * extractZ --
 **********************************************************************/
void
surf::extractZ (
		const CurrentCompMode currentCompMode,
		const PreconditionerType preCondType,
		Formulation& formulation,
		Preconditioner& preconditioner,
		const vector<CondInfo*>& condInfoPtrList,
		const Mesh& condMesh,
		const SurfConst& surfConst,
		const SolverType solverType,
		const SimulationType simuType)
{
  if (simuType == MQS) {
    extract_MQS_Z(currentCompMode, preCondType, formulation, preconditioner,
		  condInfoPtrList, condMesh, surfConst, solverType);
  } else {
    extract_EMQS_fullwave_Z(currentCompMode, preCondType, formulation, preconditioner,
			    condInfoPtrList, condMesh, surfConst, solverType);
  }
}

/**********************************************************************
 * extract_MQS_Z --
 * Without displacement current, the conduction current at one end of a 
 * conductor is equal to the conduction current at the other end. Hence
 * this two ends can be regarded as two terminals of a single port.
 * Only one excitation per conductor is needed.
 **********************************************************************/
void
surf::extract_MQS_Z (
		     const CurrentCompMode currentCompMode,
		     const PreconditionerType preCondType,
		     Formulation& formulation,
		     Preconditioner& preconditioner,
		     const vector<CondInfo*>& condInfoPtrList,
		     const Mesh& condMesh,
		     const SurfConst& surfConst,
		     const SolverType solverType)
{
  ComplexVector RHS(formulation.totalNumUnknown());
  ComplexVector solution(formulation.totalNumUnknown());
  size_t numCond = condMesh.numCond();
  TNT::Matrix<std::complex<double> > Ym(numCond, numCond, 0.);

  double leftContactVolt = 1;
  double rightContactVolt = 0;  
  double voltDiff = leftContactVolt - rightContactVolt;
  PortCurrent portCurrent(currentCompMode, condMesh, condInfoPtrList, 
			  formulation, formulation.frequency(), voltDiff,
			  surfConst); 
  TNT::stopwatch time;
  size_t colIndex = 0;
  for (size_t sCondIndex = 0; sCondIndex < numCond; sCondIndex++) {
    if ( condInfoPtrList[sCondIndex]->condType() != CondInfo::GROUND) {
      formulation.setupRHS(sCondIndex, leftContactVolt, rightContactVolt, RHS);

      time.reset();  time.start();
      if (solverType == DIRECT) {
	formulation.directSolve(RHS, solution);
	formulation.checkCondNum();
	//	formulation.checkSolution(RHS, solution);
      } else {
	iterativeSolver(preCondType, formulation, preconditioner, RHS, solution);    
      }
      timeReport("Time for solving the system := ", time.read());

      size_t rowIndex = 0;
      for (size_t tCondIndex = 0; tCondIndex < numCond; tCondIndex++) {
	if ( condInfoPtrList[tCondIndex]->condType() != CondInfo::GROUND) {
	  complex<double> leftCurrent, rightCurrent;
	  portCurrent(solution, sCondIndex, tCondIndex, leftCurrent, rightCurrent);
	  // take the average to eliminate the round off error. left and right 
	  // current should be the same for MQS.
	  complex<double> current = (leftCurrent - rightCurrent) / 2.;
	  Ym[rowIndex][colIndex] = current / voltDiff;
	  rowIndex ++;
	}
      }
      colIndex ++;      
    }
  }

  if (colIndex != numCond) {
    size_t matSize = colIndex;
    Ym.resize(matSize, matSize);
  }
  /* now Ym has been overwritten by Zm */
  invert(Ym);

  outputZm(Ym, formulation.frequency());
}

/**********************************************************************
 * extract_EMQS_fullwave_Z --
 * With displacement current, the conduction current at one end of a 
 * conductor is NOT equal to the conduction current at the other end. 
 * Hence this two ends can NOT be regarded as two terminals of a single port.
 * I have to treat this conductor as a two port circuit at first. 
 * This means the excitation has to be posed individually on each contact.
 * The input impedance looking into the final port, which consists of the
 * two ends of this conductor could only be obtained through circuit
 * theory. 
 * Two excitation per conductor are needed, if the conductor only has two contacts. 
 * So it is more time consuming than MQS analysis.
 **********************************************************************/
void
surf::extract_EMQS_fullwave_Z (
			       const CurrentCompMode currentCompMode,
			       const PreconditionerType preCondType,
			       Formulation& formulation,
			       Preconditioner& preconditioner,
			       const vector<CondInfo*>& condInfoPtrList,
			       const Mesh& condMesh,
			       const SurfConst& surfConst,
			       const SolverType solverType)
{
  ComplexVector RHS(formulation.totalNumUnknown());
  ComplexVector solution(formulation.totalNumUnknown());
  size_t numCond = condMesh.numCond();
  TNT::Matrix<std::complex<double> > Ym(2*numCond, 2*numCond, 0.);

  vector<double> leftContactVolt(2);
  vector<double> rightContactVolt(2);  
  leftContactVolt[0] = 1;
  rightContactVolt[0] = 0;
  leftContactVolt[1] = 0;
  rightContactVolt[1] = 1;

  TNT::stopwatch time;
  size_t colIndex = 0;
  for (size_t sCondIndex = 0; sCondIndex < numCond; sCondIndex++) {
    if ( condInfoPtrList[sCondIndex]->condType() != CondInfo::GROUND) {

      for (size_t ci = 0; ci < 2; ci++) {
	// excite each contact separately
	double voltDiff = leftContactVolt[ci] - rightContactVolt[ci];
	PortCurrent portCurrent(currentCompMode, condMesh, condInfoPtrList, 
				formulation, formulation.frequency(), voltDiff,
				surfConst); 
	formulation.setupRHS(sCondIndex, leftContactVolt[ci], rightContactVolt[ci], RHS);
	
	time.reset();  time.start();
	if (solverType == DIRECT) {
	  formulation.directSolve(RHS, solution);
	  formulation.checkCondNum();
	  //	formulation.checkSolution(RHS, solution);
	} else {
	  iterativeSolver(preCondType, formulation, preconditioner, RHS, solution);    
	}
	timeReport("Time for solving the system := ", time.read());

	size_t rowIndex = 0;
	for (size_t tCondIndex = 0; tCondIndex < numCond; tCondIndex++) {
	  if ( condInfoPtrList[tCondIndex]->condType() != CondInfo::GROUND) {
	    complex<double> leftCurrent, rightCurrent;
	    portCurrent(solution, sCondIndex, tCondIndex, leftCurrent, rightCurrent);
	    Ym[rowIndex][colIndex] = leftCurrent / abs(voltDiff);
	    rowIndex ++;
	    Ym[rowIndex][colIndex] = rightCurrent / abs(voltDiff);
	    rowIndex ++;
	  }
	}
	colIndex ++;      
      }
    }
  }

  size_t matSize = colIndex;
  if (colIndex != 2*numCond) {
    Ym.resize(matSize, matSize);
  }

  /*
  complex<double> **tmp;
  tmp = new complex<double>*[matSize];
  for (size_t i=0; i<matSize; i++) {
    tmp[i] = new complex<double>[matSize];
    for (size_t j=0; j<matSize; j++) {
      tmp[i][j] = Ym[i][j];
    }    
  }
  pseudo_inverse(tmp, matSize, matSize);
  for (size_t i=0; i<matSize; i++) {
    for (size_t j=0; j<matSize; j++) {
      Ym[i][j] = tmp[i][j];
    }    
  }
  for (size_t i=0; i<matSize; i++) {
    delete tmp[i];
  }
  delete tmp;
  */

  // At low frequency, when leftCurrent ~= rightCurrent, the matrix Ym is
  // nearly singular. The right way is to do pseudo inverse. 
  // For the moment, user has to know which mode to pick. If he uses
  // EMQS to calculate low-f Z, the results won't be good.
  /* now Ym has been overwritten by Zm */
  invert(Ym);

  TNT::Matrix<std::complex<double> > T(colIndex/2, colIndex, 0.);
  for (size_t i=0; i < colIndex/2; i++) {
    T[i][2*i] = 1.;
    T[i][2*i+1] = -1.;
  }
  Ym = T * Ym * transpose(T);

  outputZm(Ym, formulation.frequency());
}

/**********************************************************************
 * outputZm --
 **********************************************************************/
void
surf::outputZm (
		const TNT::Matrix<std::complex<double> >& Zm,
		const double freq)
{ 
  ofstream RFile("R.dat", std::ios::app);
  ofstream LFile("L.dat", std::ios::app);
  ofstream ZmFile("Zm.dat", std::ios::app);
  for (size_t row = 0; row < Zm.num_rows(); row++) {
    for (size_t col = 0; col < Zm.num_cols(); col++) {
      double R = real(Zm[row][col]);
      double L = imag(Zm[row][col]);
      if (freq != 0) {
	L = L/(2*PI*freq);
      } else {
	L = 0;
      }
      cout << endl 
	   << "\t Zm[" << row << ", " << col << "] = " << Zm[row][col] 
	   << endl;
      cout << "\t f = " << freq << "\tR = " << R << "\tL = " << L << endl;
      ZmFile << freq 
	     << "\t Zm[" << row << ", " << col << "] = " << Zm[row][col] 
	     << endl;
      RFile << freq << "\t" << R << endl;
      LFile << freq << "\t" << L << endl;
    }
  }
}

/**********************************************************************
 * extractTransLineZ --
 **********************************************************************/
void
surf::extractTransLineZ (
			 const CurrentCompMode currentCompMode,
			 const PreconditionerType preCondType,
			 Formulation& formulation,
			 Preconditioner& preconditioner,
			 const vector<CondInfo*>& condInfoPtrList,
			 const Mesh& condMesh,
			 const SurfConst& surfConst,
			 const SolverType solverType,
			 const SimulationType simuType)
{
  if (simuType == MQS) {
    extractTransLineZ_MQS(currentCompMode, preCondType, formulation, preconditioner,
			  condInfoPtrList, condMesh, surfConst, solverType);
  } else {
    extractTransLineZ_EMQS_fullwave(currentCompMode, preCondType, formulation, preconditioner,
				    condInfoPtrList, condMesh, surfConst, solverType);
  }
}

/**********************************************************************
 * extractTransLineZ_MQS --
 **********************************************************************/
void
surf::extractTransLineZ_MQS (
			     const CurrentCompMode currentCompMode,
			     const PreconditionerType preCondType,
			     Formulation& formulation,
			     Preconditioner& preconditioner,
			     const vector<CondInfo*>& condInfoPtrList,
			     const Mesh& condMesh,
			     const SurfConst& surfConst,
			     const SolverType solverType)
{
  ComplexVector RHS(formulation.totalNumUnknown());
  ComplexVector solution(formulation.totalNumUnknown());
  size_t numCond = condMesh.numCond();

  vector<double> leftContactVolt(numCond, 0);
  vector<double> rightContactVolt(numCond, 0);
  int sCondIndex = -1;
  int tCondIndex = -1;
  for (size_t condIndex=0; condIndex < numCond; condIndex++) {
    if (condInfoPtrList[condIndex]->condType() == CondInfo::GROUND)  continue;
    // only excite one conductor and extract one column of [Z].
    // The other column is symmetric
    if (sCondIndex == -1) {
      leftContactVolt[condIndex] = 1.;
      rightContactVolt[condIndex] = -1;
      sCondIndex = condIndex;
    } else {
      tCondIndex = condIndex;
      break;
    }
  }
  PortCurrent portCurrent(currentCompMode, condMesh, condInfoPtrList, 
			  formulation, formulation.frequency(), 
			  1., surfConst); 

  TNT::stopwatch time;   time.start();
  formulation.setupRHS(leftContactVolt, rightContactVolt, RHS);
  if (solverType == DIRECT) {
    formulation.directSolve(RHS, solution);
    formulation.checkCondNum();
  } else {
    iterativeSolver(preCondType, formulation, preconditioner, RHS, solution);    
  }
  timeReport("Time for solving the system := ", time.read());

  complex<double> leftCurrent, rightCurrent;
  portCurrent(solution, sCondIndex, sCondIndex, leftCurrent, rightCurrent);
  complex<double> I0 = (leftCurrent - rightCurrent) / 2.;
  portCurrent(solution, sCondIndex, tCondIndex, leftCurrent, rightCurrent);
  complex<double> I1 = (leftCurrent - rightCurrent) / 2.;
#ifdef DEBUG_TRANSLINE
  cout << "I0 := " << I0 << "  " << "I1 := " << I1 << endl;
#endif

  TNT::Matrix<std::complex<double> > Y(numCond, numCond, 0.);
  double voltDiff = leftContactVolt[sCondIndex] - rightContactVolt[sCondIndex];
  Y[0][0] = I0 / abs(voltDiff);
  Y[0][1] = I1 / abs(voltDiff);
  Y[1][0] = Y[0][1];
  Y[1][1] = Y[0][0];
  invert(Y);

#ifdef DEBUG_TRANSLINE
  cout << "z00 := " << Y[0][0] << "  " << "z01 := " << Y[0][1] << endl;
#endif
  complex<double> z = 2.*(Y[0][0]-Y[0][1]);

  double freq = formulation.frequency();
  cout << endl << "f = " << freq
       << "\tZ = " << z << "\t|Z| = " << abs(z) << endl;

  ofstream ZFile("z.dat", std::ios::app);
  ZFile << freq << "\t" << z << "\t" << abs(z) << endl;  

  ofstream ZPlotFile("zplot.dat", std::ios::app);
  ZPlotFile << freq << "\t" << abs(z) << endl;  
}

/**********************************************************************
 * extractTransLineZ_EMQS_fullwave --
 **********************************************************************/
void
surf::extractTransLineZ_EMQS_fullwave (
				       const CurrentCompMode currentCompMode,
				       const PreconditionerType preCondType,
				       Formulation& formulation,
				       Preconditioner& preconditioner,
				       const vector<CondInfo*>& condInfoPtrList,
				       const Mesh& condMesh,
				       const SurfConst& surfConst,
				       const SolverType solverType)
{
  ComplexVector RHS(formulation.totalNumUnknown());
  ComplexVector solution(formulation.totalNumUnknown());
  size_t numCond = condMesh.numCond();

  int sCondIndex = -1;
  int tCondIndex = -1;
  for (size_t condIndex=0; condIndex < numCond; condIndex++) {
    if (condInfoPtrList[condIndex]->condType() == CondInfo::GROUND)  continue;
    if (sCondIndex == -1) {
      sCondIndex = condIndex;
    } else {
      tCondIndex = condIndex;
      break;
    }
  }

  double leftContactVolt = 1;
  double rightContactVolt = 0;
  // excite one contact only due to symmetric structure
  double voltDiff = leftContactVolt - rightContactVolt;
  PortCurrent portCurrent(currentCompMode, condMesh, condInfoPtrList, 
			  formulation, formulation.frequency(), voltDiff,
			  surfConst); 

  TNT::stopwatch time;   time.start();
  formulation.setupRHS(sCondIndex, leftContactVolt, rightContactVolt, RHS);
  if (solverType == DIRECT) {
    formulation.directSolve(RHS, solution);
    formulation.checkCondNum();
  } else {
    iterativeSolver(preCondType, formulation, preconditioner, RHS, solution);    
  }
  timeReport("Time for solving the system := ", time.read());


   // so we got the impedance matrix
   // Here we assume that port indices are ordered strictly
   // as port0 of 1st cond, port 1 of 1st cond, port 0 of 2nd
   // cond, port1 of 2nd cond, ...
  TNT::Matrix<std::complex<double> > Ym(4, 4, 0.);
  complex<double> leftCurrent, rightCurrent;

  portCurrent(solution, sCondIndex, sCondIndex, leftCurrent, rightCurrent);
  Ym[0][0] = leftCurrent / abs(voltDiff);
  Ym[1][0] = rightCurrent / abs(voltDiff);

  portCurrent(solution, sCondIndex, tCondIndex, leftCurrent, rightCurrent);
  Ym[2][0] = leftCurrent / abs(voltDiff);
  Ym[3][0] = rightCurrent / abs(voltDiff);

  Ym[0][1] = Ym[1][0];
  Ym[1][1] = Ym[0][0];
  Ym[2][1] = Ym[3][0];
  Ym[3][1] = Ym[2][0];
  
  Ym[0][2] = Ym[2][0];
  Ym[1][2] = Ym[3][0];
  Ym[2][2] = Ym[0][0];
  Ym[3][2] = Ym[1][0];

  Ym[0][3] = Ym[3][0];
  Ym[1][3] = Ym[2][0];
  Ym[2][3] = Ym[1][0];
  Ym[3][3] = Ym[0][0];

  // Treat left contacts of two conductors as one port, and right contacts as
  // the other port.
  // [1 0 -1 0] * IYm [1  0] [I1] = [V1]
  // [0 1 0 -1]       [0  1] [I2]   [V2]
  //                  [-1 0]
  //                  [0 -1]
  TNT::Matrix<std::complex<double> > T(2, 4, 0.);
  T[0][0] = 1.;
  T[0][2] = -1.;
  T[1][1] = 1.;
  T[1][3] = -1.;

  TNT::Matrix<std::complex<double> > Z(2, 2);
  invert(Ym);
  // This is the impedance matrix for the two-port T-line
  Z = T * Ym * transpose(T);
  // Now I have the admittance matrix for the two-port T-line
  invert(Z);
  // The admittance for the shorted T-line is just Y = Z[0][0] and its impedance
  // is obviously 1/Y;
  complex<double> z = 1./Z[0][0];

  double freq = formulation.frequency();
  cout << endl << "f = " << freq
       << "\tZ = " << z << "\t|Z| = " << abs(z) << endl;

  ofstream ZFile("z.dat", std::ios::app);
  ZFile << freq << "\t" << z << "\t" << abs(z) << endl;  

  ofstream ZPlotFile("zplot.dat", std::ios::app);
  ZPlotFile << freq << "\t" << abs(z) << endl;  
}

/**********************************************************************
 * extractOneColY --
 **********************************************************************/
void
surf::extractOneColY (
		      const CurrentCompMode currentCompMode,
		      const PreconditionerType preCondType,
		      Formulation& formulation,
		      Preconditioner& preconditioner,
		      const vector<CondInfo*>& condInfoPtrList,
		      const Mesh& condMesh,
		      const int oneColIndex,
		      const SurfConst& surfConst,
		      const SolverType solverType)
{
  ComplexVector RHS(formulation.totalNumUnknown());
  ComplexVector solution(formulation.totalNumUnknown());
  size_t numCond = condMesh.numCond();
  TNT::Matrix<std::complex<double> > Ym(numCond, numCond, 0.);

  double leftContactVolt = 1;
  double rightContactVolt = -1;  
  double voltDiff = leftContactVolt - rightContactVolt;
  PortCurrent portCurrent(currentCompMode, condMesh, condInfoPtrList, 
			  formulation, formulation.frequency(), voltDiff, 
			  surfConst); 

  TNT::stopwatch time;   time.start();
  size_t sCondIndex = oneColIndex;
  formulation.setupRHS(sCondIndex, leftContactVolt, rightContactVolt, RHS);
  if (solverType == DIRECT) {
    formulation.directSolve(RHS, solution);
  } else {
    iterativeSolver(preCondType, formulation, preconditioner, RHS, solution);    
  }
  timeReport("Time for solving the system := ", time.read());

  size_t rowIndex = 0;
  for (size_t tCondIndex = 0; tCondIndex < numCond; tCondIndex++) {
    if ( condInfoPtrList[tCondIndex]->condType() != CondInfo::GROUND) {
      complex<double> leftCurrent, rightCurrent;
      portCurrent(solution, sCondIndex, tCondIndex, leftCurrent, rightCurrent);
      complex<double> current = (leftCurrent - rightCurrent) / 2.;;
      Ym[rowIndex][0] = current / abs(voltDiff);
      rowIndex ++;
    }
  }

  if (rowIndex != numCond) {
    Ym.resize(rowIndex, rowIndex);
  }

  outputZm(Ym, formulation.frequency());
}

/**********************************************************************
 * iterativeSolver --
 **********************************************************************/
void
surf::iterativeSolver (
		       const PreconditionerType preCondType,
		       Formulation& formulation,
		       Preconditioner& preconditioner,
		       const ComplexVector& RHS,
		       ComplexVector& solution)
{
  for (size_t i = 0; i < formulation.totalNumUnknown(); i++) {
    solution[i] = 0.;
  }

#ifdef HIGH_ACCURACY
  double tol = 1e-6;
#else
  double tol = 1e-3;
#endif

  if (preconditioner.isLeftPreconditioner()) {
    // left pre-conditioner amplifies the residual
    tol *= 1e-3;
  }
#ifdef DEBUG_ITERATIVE
  size_t restart = 100;
#else
  size_t restart = formulation.getRestart();
#endif
  size_t max_iter = 10 * restart;
  bool converged;
  try {
    converged = gmres(formulation, solution, RHS, 
		      preconditioner, max_iter, restart, tol);
  }
  catch (std::domain_error e) {
    cout << e.what() << endl;
    errorMessage("extracZ.cc::iterativeSolver()",  "Bug in gmres");
  }

  if (!converged) {
    warningMessage("extracZ.cc::iterativeSolver()", "gmres fails to converge! ");
  }
}

