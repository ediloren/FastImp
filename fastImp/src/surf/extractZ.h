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

  const static char cvsid[] = "$Id: extractZ.h,v 1.7 2003/05/01 18:04:00 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _EXTRACT_Z_H_
#define _EXTRACT_Z_H_

#include <vector>
#include <complex>
#include "formulation.h"
#include "preconditioner.h"
#include "portCurrent.h"
#include "condInfo.h"
#include "mesh.h"
#include "surfConst.h" // for ComplexVector
#include "cmat.h"

namespace surf {

  void
  extractZ (
	    const CurrentCompMode currentCompMode,
	    const PreconditionerType preCondType,
	    Formulation& formulation,
	    Preconditioner& preconditioner,
	    const std::vector<CondInfo*>& condInfoPtrList,
	    const mesh::Mesh& condMesh,
	    const SurfConst& surfConst,
	    const SolverType solverType,
	    const SimulationType simuType);

  void
  extract_MQS_Z (
		 const CurrentCompMode currentCompMode,
		 const PreconditionerType preCondType,
		 Formulation& formulation,
		 Preconditioner& preconditioner,
		 const std::vector<CondInfo*>& condInfoPtrList,
		 const mesh::Mesh& condMesh,
		 const SurfConst& surfConst,
		 const SolverType solverType);
  void
  extract_EMQS_fullwave_Z (
			   const CurrentCompMode currentCompMode,
			   const PreconditionerType preCondType,
			   Formulation& formulation,
			   Preconditioner& preconditioner,
			   const std::vector<CondInfo*>& condInfoPtrList,
			   const mesh::Mesh& condMesh,
			   const SurfConst& surfConst,
			   const SolverType solverType);
  void 
  extractOneColY (
		  const CurrentCompMode currentCompMode,
		  const PreconditionerType preCondType,
		  Formulation& formulation,
		  Preconditioner& preconditioner,
		  const std::vector<CondInfo*>& condInfoPtrList,
		  const mesh::Mesh& condMesh,
		  const int oneColIndex,
		  const SurfConst& surfConst,
		  const SolverType solverType);
  
  void
  extractTransLineZ (
		     const CurrentCompMode currentCompMode,
		     const PreconditionerType preCondType,
		     Formulation& formulation,
		     Preconditioner& preconditioner,
		     const std::vector<CondInfo*>& condInfoPtrList,
		     const mesh::Mesh& condMesh,
		     const SurfConst& surfConst,
		     const SolverType solverType,
		     const SimulationType simuType);

  void
  extractTransLineZ_MQS (
			 const CurrentCompMode currentCompMode,
			 const PreconditionerType preCondType,
			 Formulation& formulation,
			 Preconditioner& preconditioner,
			 const std::vector<CondInfo*>& condInfoPtrList,
			 const mesh::Mesh& condMesh,
			 const SurfConst& surfConst,
			 const SolverType solverType);

  void
  extractTransLineZ_EMQS_fullwave (
				   const CurrentCompMode currentCompMode,
				   const PreconditionerType preCondType,
				   Formulation& formulation,
				   Preconditioner& preconditioner,
				   const std::vector<CondInfo*>& condInfoPtrList,
				   const mesh::Mesh& condMesh,
				   const SurfConst& surfConst,
				   const SolverType solverType);
  
  void
  iterativeSolver (
		   const PreconditionerType preCondType,
		   Formulation& formulation,
		   Preconditioner& preconditioner,
		   const ComplexVector& RHS,
		   ComplexVector& solution);
  void
  outputZm (
	    const TNT::Matrix<std::complex<double> >& Zm,
	    const double freq);

  void
  outputTransLineZ (
		    const std::vector<std::vector<std::complex<double> > >& Zm,
		    const double freq);

} //namespace surf

#endif
