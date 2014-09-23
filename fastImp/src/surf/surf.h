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

#ifndef __SURF_H_
#define __SURF_H_

#include <complex>
#include <string> // for exit()
#include "vec.h" // for TNT::Vector
#include "service.h"
#include "surfConst.h"
#include "frequencySample.h"
#include "condInfo.h"
#include "mesh.h"
#include "element.h"
#include "formulation.h"
#include "preconditioner.h"
#include "extractZ.h"

namespace surf {

  class Surf {

  public:
    Surf(
	 char* inputStructFileIn,
	 char* outputMeshFileIn,
	 char* inputMeshFileIn,
	 const SimulationType simuTypeIn,
	 const PreconditionerType preCondTypeIn,
	 const bool useGlobalCoord,
	 const FormulationType formulationTypeIn,
	 const CurrentCompMode currentCompModeIn,
	 const SolverType solverType,
	 const bool useNonUniformMeshIn,
	 const bool aspectRatioWarningIn,
	 const bool isTransmissionLineIn,
	 const bool useQuadPanelIn,
	 const FrequencySample& freqSamIn,
	 const int oneColIndexIn,
	 const std::vector<surf::CondInfo*>& condInfoPtrList);
    void extract(void);

  private:
    std::vector<surf::CondInfo*> condInfoPtrList;
    mesh::Mesh condMesh;

    char* inputStructFile;
    char* outputMeshFile;
    char* inputMeshFile;
    SimulationType simuType;
    PreconditionerType preCondType;
    bool useGlobalCoord;
    FormulationType formulationType;
    CurrentCompMode currentCompMode;
    SolverType solverType;
    bool useNonUniformMesh;
    bool checkPanelSize_;
    bool isTransmissionLine;
    bool useQuadPanel;
    FrequencySample freqSam;
    int oneColIndex;
    bool isMeshGenerated_;

    SampleMethod samMethod;
    SurfConst surfConst_;
   
    void statusReport (const double freq) const;
    void checkPanelSize(const double f, const Formulation& formulation) const;
    void checkPanelSize(const double f) const;
  };

} // namespace surf

#endif


