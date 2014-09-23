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

#ifndef __SURF_CONST_H_
#define __SURF_CONST_H_

#include <string>
#include <vector>
#include <cfloat> // for DBL_EPSILON
#include "vec.h" // for TNT::Vector
#include "service.h"

namespace surf {

  typedef TNT::Vector<std::complex<double> > ComplexVector;

  const double PI = 3.14159265358979;

  enum SimulationType {MQS, EMQS, FULL_WAVE, NUM_SIMULATION_TYPE };
  enum FormulationType { OLD, NEW, NUM_FORMULATION_TYPE };
  enum PreconditionerType { NO_PRE_CONDITIONER, LEFT_PRE_CONDITIONER,
			    RIGHT_PRE_CONDITIONER, NUM_PRE_CONDITIONER_TYPE };
  enum CurrentCompMode {AUTO_MODE, LOW_FREQUENCY_MODE, HIGH_FREQUENCY_MODE,
			TWO_DIM_MODE, ENERGY_MODE, NUM_CURRENT_COMP_MODE };
  enum SolverType {DIRECT, ITERATIVE, ITERATIVE_WITH_PFFT, NUM_SOLVER_TYPE };

  typedef struct SurfConst {
    double EPSILON0;
    double MU0;
    double normFactor;
  } SurfConst;

} // namespace surf

#endif


