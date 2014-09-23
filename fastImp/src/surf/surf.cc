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

  const static char cvsid[] = "$Id: surf.cc,v 1.23 2003/07/16 15:32:33 zhzhu Exp $";

  ==========================================================================
*/

#include <iostream>
#include <stdio.h>
#include "surf.h"

using namespace surf;
using namespace std;
using namespace mesh;
using namespace pfft;

const char* simuTypeName[3] = {"MQS", "EMQS", "FULL_WAVE"};
const char* formulationTypeName[2] = {"OLD", "NEW"};
const char* freqSamMethodName[2] = {"LOGRITHMIC", "LINEAR"};
const char* preCondTypeName[3] = {"NONE", "LEFT_PRE_CONDITIONER", 
				  "RIGHT_PRE_CONDITIONER"}; 
const char* accelarartorName[3] = {"NONE", "pFFT"}; 
const char* currentCompModeName[5] = {"AUTO", "LOW_FREQUENCY_MODE", 
				      "HIGH_FREQUENCY_MODE", "2D_MODE",
				      "ENERGY_MODE" };

/**********************************************************************
 * Surf --
 **********************************************************************/
Surf::Surf (
	    char* inputStructFileIn,
	    char* outputMeshFileIn,
	    char* inputMeshFileIn,
	    const SimulationType simuTypeIn,
	    const PreconditionerType preCondTypeIn,
	    const bool useGlobalCoordIn,
	    const FormulationType formulationTypeIn,
	    const CurrentCompMode currentCompModeIn,
	    const SolverType solverTypeIn,
	    const bool useNonUniformMeshIn,
	    const bool checkPanelSizeIn,
	    const bool isTransmissionLineIn,
	    const bool useQuadPanelIn,
	    const FrequencySample& freqSamIn,
	    const int oneColIndexIn,
	    const vector<surf::CondInfo*>& condInfoPtrListIn)
  :inputStructFile(inputStructFileIn), outputMeshFile(outputMeshFileIn), 
  inputMeshFile(inputMeshFileIn), simuType(simuTypeIn), 
  preCondType(preCondTypeIn), useGlobalCoord(useGlobalCoordIn),
  formulationType(formulationTypeIn),
  currentCompMode(currentCompModeIn), solverType(solverTypeIn), 
  useNonUniformMesh(useNonUniformMeshIn), checkPanelSize_(checkPanelSizeIn),
  isTransmissionLine(isTransmissionLineIn), useQuadPanel(useQuadPanelIn), 
  freqSam(freqSamIn), oneColIndex(oneColIndexIn), 
  condInfoPtrList(condInfoPtrListIn)
{
  if (inputMeshFile) {
    condMesh.readMeshFile(condInfoPtrList, inputMeshFile);
    isMeshGenerated_ = false;
  } else {
    if (outputMeshFile) checkPanelSize_ = false;
    condMesh.generateMesh(condInfoPtrList, useNonUniformMesh, 
			  useQuadPanel, checkPanelSize_);
    isMeshGenerated_ = true;
  }

  if (outputMeshFile) {
    if (inputMeshFile) { // we already have patran format, so output fastcap
      condMesh.outputMesh(outputMeshFile, FASTCAP);
    } else {
      condMesh.outputMesh(outputMeshFile, PATRAN);
    }
    exit(0);
  }

  // normalize the mesh size to meter range to avoid bad scale in system matrix
#ifdef DISABLE_SCALING
  surfConst_.normFactor = 1.;  // without scaling
#else
  surfConst_.normFactor = 1. / condMesh.unit();   // with scaling
#endif

  double er = 1.;
  surfConst_.EPSILON0 = er * 1e-9 / (36*PI) / surfConst_.normFactor;
  surfConst_.MU0 = 4e-7 * PI / surfConst_.normFactor;

  if (isMeshGenerated_) {
    for (int i=0; i < freqSam.numSam(); i++) {
      if (checkPanelSize_) checkPanelSize(freqSam.point(i));
    }
  }
}

/**********************************************************************
 * extract --
 **********************************************************************/
void
Surf::extract (
	       void)
{
  Formulation formulation(condMesh, condInfoPtrList, simuType, formulationType,
			  solverType, surfConst_, useGlobalCoord);
  Preconditioner preconditioner(formulation, preCondType, 
				surfConst_, useGlobalCoord);

  TNT::stopwatch totalTime;
  totalTime.start();
  TNT::stopwatch OneFreqTime;
  for (int i=0; i < freqSam.numSam(); i++) {
    statusReport(freqSam.point(i));
    if (checkPanelSize_) checkPanelSize(freqSam.point(i), formulation);

    OneFreqTime.reset();  OneFreqTime.start();
    formulation.setupSystem(freqSam.point(i));
    if (solverType != DIRECT) {
      preconditioner.setup(freqSam.point(i));
    }
    if (oneColIndex >= 0) {
      extractOneColY(currentCompMode, preCondType, formulation, preconditioner,
      		     condInfoPtrList, condMesh, oneColIndex, surfConst_, solverType);
    } else if (isTransmissionLine) {
      extractTransLineZ(currentCompMode, preCondType, formulation, preconditioner,
      			condInfoPtrList, condMesh, surfConst_, solverType, simuType);
    } else {
      extractZ(currentCompMode, preCondType, formulation, preconditioner,
	       condInfoPtrList, condMesh, surfConst_, solverType, simuType);
    }
    surf::timeReport("Total time for extraction at one frequency point := ", 
		     OneFreqTime.read());
  }
  surf::timeReport("Total run time of fastImp := ", totalTime.read());
  std::cout << "\n\n\tAll impedance matrices dumped to file Zm.dat\n" 
	    << "\tReal and imaginary parts dumped to file R.dat and L.dat\n";
}

/**********************************************************************
 * statusReport --
 **********************************************************************/
void
Surf::statusReport (
		    const double freq) const
{
  cout << endl
       << "\t Sampling frequency := " << freq << endl
       << "\t Simulation Type is " << simuTypeName[simuType] << endl
       << "\t Formulation Type used is " << formulationTypeName[formulationType] << endl
       << "\t Current Computation mode is " << currentCompModeName[currentCompMode] << endl;

  if (solverType == DIRECT) {
    cout << "\t A direct solver is used" << endl;
  } else {
    if (solverType == ITERATIVE) {
      cout << "\t An iterative solver is used" << endl;
    } else {
      cout << "\t An iterative solver with pFFT is used" << endl;
    }

    cout << "\t Pre-condtioner is " << preCondTypeName[preCondType] << endl;
    if (useGlobalCoord) {
      cout << "\t Pre-conditioner is generated in global coord." << endl;
    } else {
      cout << "\t Pre-conditioner is generated in local coord." << endl;
    }
  }

  if (isTransmissionLine) {
    cout << "\t The two conductors are shorted at one end" << endl;
  }

  if (oneColIndex >= 0) {
    cout << "\t Only " << oneColIndex 
	 << "th column of admmitance matrix will be computed" << endl;
  } 
}

/**********************************************************************
 * checkPanelSize --
 * This has to be called after the formulation object has been created.
 **********************************************************************/
void
Surf::checkPanelSize(
		     const double f,
		     const Formulation& formulation) const
{
  const double threshold_1 = 8.;
  const double threshold_2 = 0.5;

  if (simuType != MQS) {

    double waveLen = 3e8 / f;
    double maxPanelSize = formulation.maxElementSize() / surfConst_.normFactor;
    if (waveLen / maxPanelSize < threshold_1) {
      cout << endl 
	   << "\t wave length := " << waveLen 
	   << ", max panel size := " << maxPanelSize;
      surf::warningMessage("surf.cc : checkPanelSize",
			   "Less than 8 panels per wave length!! Please modefy your input file. ");
    }

  } else {

    // If the mesh is generated by my own mesher, it wass already checked
    if (! isMeshGenerated_) {
      for (size_t i = 0; i < condMesh.numCond(); i++) {
	double sigma = condInfoPtrList[i]->conductivity();
	double w = 2*PI*f;
	double skinDepth = sqrt(2. / (w * surfConst_.MU0 * sigma)) / surfConst_.normFactor;
	double minPanelSize;
	minPanelSize = formulation.minElementSize(i) / surfConst_.normFactor;
	if (skinDepth / minPanelSize < threshold_2) {
	  cout << endl 
	       << "\t skinDepth := " << skinDepth
	       << ", min panel size := " << minPanelSize;
	  surf::warningMessage("surf.cc : checkPanelSize",
			       "min panel size is bigger than skinDepth!! Please check the input file. ");
	}
      }
    }
  }
}

/**********************************************************************
 * checkPanelSize --
 * This can be called before the formulation object has been created
 * It is intended to be used for the mesh generated by my own mesher.
 **********************************************************************/
void
Surf::checkPanelSize(
		     const double f) const
{
  const double threshold_2 = 0.5;
  for (size_t i = 0; i < condMesh.numCond(); i++) {
    double sigma = condInfoPtrList[i]->conductivity();
    double w = 2*PI*f;
    double skinDepth = sqrt(2. / (w * surfConst_.MU0 * sigma)) / surfConst_.normFactor;
    double minPanelSize = condMesh.minStep(i) / surfConst_.normFactor;
    if (skinDepth / minPanelSize < threshold_2) {
      cout << endl 
	   << "\t skinDepth := " << skinDepth << " at " << f << " Hz"
	   << ", min panel size := " << minPanelSize;
      surf::errorMessage("surf.cc : checkPanelSize",
			 "min panel size is bigger than skinDepth!! Please check the input file. ");
    }
  }
}

  
  
