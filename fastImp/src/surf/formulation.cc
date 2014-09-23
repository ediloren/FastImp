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

  const static char cvsid[] = "$Id: formulation.cc,v 1.28 2003/07/10 19:59:16 zhzhu Exp $";

  ==========================================================================
*/

#include <stdexcept>
#include "formulation.h"

using namespace surf;
using namespace std;
using namespace mesh;
using namespace pfft;

/**********************************************************************
 * Formulation --
 **********************************************************************/
Formulation::Formulation (
			  const Mesh& condMeshIn, 
			  const vector<CondInfo*>& condInfoPtrListIn,
			  const SimulationType simuType,
			  const FormulationType formulationType,
			  const SolverType solverTypeIn,
			  const SurfConst& surfConst,
			  const bool useGlobalCoord)
  :condMesh(condMeshIn),condInfoPtrList(condInfoPtrListIn), simuType_(simuType),
   formulationType_(formulationType), solverType_(solverTypeIn),
   useGlobalCoord_(useGlobalCoord),
   numCond_(condMesh.numCond()), totalNumE_(condMeshIn.totalNumPanel()), 
   totalNumPhi_(condMeshIn.totalNumNode()), 
   totalNumRho_(condMeshIn.totalNumPanel()),
   printMatMultVecTime_(true), surfConst_(surfConst), 
   LUFactorDone(false), formOneSystemMatDone(false)
{
  setupElementList();
  setupDualPanelList();

  setupGlobalStartRowIndex();
  setupGlobalStartColIndex();
  setupNonIntegralEqua();

  setupIntegralOperatorList();
  if (solverType_ == ITERATIVE_WITH_PFFT) {
    setupIntegralEqua_frequencyIndependentPart();
  }
  allocateWorkUnit();

  if (! useGlobalCoord_) {
    setupLocalStartRowIndex();
    setupLocalStartColIndex();
    setupLocalToGlobalMat();
    setupGlobalToLocalMat();

#ifdef DEBUG_PRE_COND
    std::ofstream fout("g2l.tmp");
    fout << globalToLocalMat;
    std::ofstream fout1("l2g.tmp");
    fout1 << localToGlobalMat;
#endif
  }

  checkMemoryUsage();
}

/**********************************************************************
 * setupSystem --
 **********************************************************************/
void 
Formulation::setupSystem (
			  const double freq)
{
  TNT::stopwatch time;  time.start();
  frequency_ = freq;
  compWaveNumber(freq, K0, KC);
  cout << "\n\t The whole system has " << totalNumRow_ << " unknowns" << endl;
  if (solverType_ == ITERATIVE_WITH_PFFT) {
    setupIntegralEqua_frequencyDependentPart();    
    timeReport("Total time for setting up frequency-dependent part of pfft := ", time.read());
  } else {
    fillIntegralEqua();
    timeReport("Total time for forming system matrix := ", time.read());
    LUFactorDone = false; // need to be reset for each frequ point
    formOneSystemMatDone = false; // need to be reset for each frequ point
  }
}

/**********************************************************************
 * setupIntegralOperatorList --
 **********************************************************************/
void
Formulation::setupIntegralOperatorList (
					void)
{
  //
  // Equa_1 
  //
  if (formulationType_ == OLD) {
    // single layer integral operator
    outerOperatorList1.push_back(NONE);
    innerOperatorList1.push_back(NONE);
    // double layer integral operator
    outerOperatorList1.push_back(NONE);
    innerOperatorList1.push_back(D_DN);
  } else if (formulationType_ == NEW) {
    // d_dn of single layer integral operator
    outerOperatorList1.push_back(D_DN);
    innerOperatorList1.push_back(NONE);
    // d_dn of double layer integral operator
    outerOperatorList1.push_back(D_DN);
    innerOperatorList1.push_back(D_DN);
  }
  calcp1List.resize(numCond_);
  if (solverType_ == ITERATIVE_WITH_PFFT) {
    kernel1List.resize(numCond_);
    pfft1List.resize(numCond_);
  }

  //
  // Equa_2a
  // 
  // single layer operator
  outerOperatorList2.push_back(NONE);
  innerOperatorList2.push_back(NONE);
  // double layer operator
  outerOperatorList2.push_back(NONE);
  innerOperatorList2.push_back(D_DN);

  //
  // Equa_3
  // 
  if (simuType_ != MQS) {
    outerOperatorList3.push_back(NONE);
    innerOperatorList3.push_back(NONE);
  }
}

/**********************************************************************
 * setupIntegralEqua_frequencyIndependentPart --
 **********************************************************************/
void
Formulation::setupIntegralEqua_frequencyIndependentPart (
							 void)
{
  TNT::stopwatch time;  time.start();

  setupEqua1_frequencyIndependentPart();
  setupEqua2a_frequencyIndependentPart();
  if (simuType_ != MQS)
    setupEqua3_frequencyIndependentPart();

  timeReport("Total time for setting up frequency-independent part of pfft := ", time.read());

}

/**********************************************************************
 * setupIntegralEqua_frequencyDependentPart --
 **********************************************************************/
void
Formulation::setupIntegralEqua_frequencyDependentPart (
						       void)
{
  setupEqua1_frequencyDependentPart();
  setupEqua2a_frequencyDependentPart();
  if (simuType_ != MQS)
    setupEqua3_frequencyDependentPart();
}

/**********************************************************************
 * setupEqua1_frequencyIndependentPart --
   Old:    Int (G1 dE/dn') - Int (E dG1/dn') = E 
   New:    d/dn { Int (G1 dE/dn') - Int (E dG1/dn') = E }
       Enforced on each conductor seperately
 **********************************************************************/
void
Formulation::setupEqua1_frequencyIndependentPart (
						  void)
{
  progressReport("setup kernel Independent part of pfft for each conductor surface");

  size_t directStencilSize, projectStencilSize, interpStencilSize;
#ifdef HIGH_ACCURACY
  directStencilSize = 6;
  projectStencilSize = 2;
  interpStencilSize = 2;
#else
  directStencilSize = 3;
  projectStencilSize = 1;
  interpStencilSize = 1;
#endif

  for (size_t i = 0; i < numCond_; i++) {
    // this is a dummy calcp with right operator lists
    calcp1List[i] = Calcp1(outerOperatorList1, innerOperatorList1, 0.);
    pfft1List[i] = Pfft1(condElementList[i], condElementList[i], calcp1List[i], 
			 directStencilSize, projectStencilSize, interpStencilSize);
  }
}

/**********************************************************************
 * setupEqua1_frequencyDependentPart --
   Old:    Int (G1 dE/dn') - Int (E dG1/dn') = E 
   New:    d/dn { Int (G1 dE/dn') - Int (E dG1/dn') = E }
       Enforced on each conductor seperately
 **********************************************************************/
void
Formulation::setupEqua1_frequencyDependentPart (
						void)
{
  progressReport("setup kernel dependent part of pfft for each conductor surface");

  TNT::stopwatch time; 
  for (size_t i = 0; i < numCond_; i++) {
    time.reset();   time.start();
    kernel1List[i] = Kernel1(KC[i]);
#ifdef DEBUG_MEMORY
    int numQudraturePoint = 0;
#else
    int numQudraturePoint = 24;
#endif
    calcp1List[i] = Calcp1(outerOperatorList1, innerOperatorList1, KC[i], numQudraturePoint);
    pfft1List[i].constructKernelDependent(kernel1List[i], calcp1List[i]);    
    timeReport("time for the pfft setup for one conductor:= ", time.read());
  }
}

/**********************************************************************
 * setupEqua2a_frequencyIndependentPart --
     t dot { Int (G0 dE/dn) - Int(E dG0/dn) } + t dot grad(phi) = 0
     Enforced on union of non-contact panels
 **********************************************************************/
void
Formulation::setupEqua2a_frequencyIndependentPart (
						   void)
{
  progressReport("setup kernel independent part of pfft for the union of all non-contact conductor surface");

  size_t directStencilSize, projectStencilSize, interpStencilSize;
#ifdef HIGH_ACCURACY
  directStencilSize = 6;
  projectStencilSize = 2;
  interpStencilSize = 2;
#else
  directStencilSize = 3;
  projectStencilSize = 1;
  interpStencilSize = 1;
#endif

  // this is a dummy calcp with right operator lists
  if (simuType_ == FULL_WAVE) {
    calcp2_FW = Calcp2_FW(outerOperatorList2, innerOperatorList2, 0.);
    pfft2_FW = Pfft2_FW(allElementList, allNonContactElementList, calcp2_FW, 
			directStencilSize, projectStencilSize, interpStencilSize);
  } else {
    calcp2_EMQS = Calcp2_EMQS(outerOperatorList2, innerOperatorList2);
    pfft2_EMQS = Pfft2_EMQS(allElementList, allNonContactElementList, calcp2_EMQS, 
			    directStencilSize, projectStencilSize, interpStencilSize);
  }
}

/**********************************************************************
 * setupEqua2a_frequencyDependentPart --
     t dot { Int (G0 dE/dn) - Int(E dG0/dn) } + t dot grad(phi) = 0
     Enforced on union of non-contact panels
 **********************************************************************/
void
Formulation::setupEqua2a_frequencyDependentPart (
						 void)
{
  progressReport("setup kernel dependent part of pfft for the union of all non-contact conductor surface");

  TNT::stopwatch time;  time.start();

#ifdef DEBUG_MEMORY
  int numQudraturePoint = 0;
  bool skipPanelIntegration = true;
#else
  int numQudraturePoint = 24;
  bool skipPanelIntegration = false;
#endif
  if (simuType_ == FULL_WAVE) {
    kernel2_FW = Kernel2_FW(K0);
    calcp2_FW = Calcp2_FW(outerOperatorList2, innerOperatorList2, K0, numQudraturePoint);
    pfft2_FW.constructKernelDependent(kernel2_FW, calcp2_FW);
  } else {
    calcp2_EMQS = Calcp2_EMQS(outerOperatorList2, innerOperatorList2, skipPanelIntegration);
    pfft2_EMQS.constructKernelDependent(kernel2_EMQS, calcp2_EMQS);
  }

  timeReport("time for the pfft setup := ", time.read());
}

/**********************************************************************
 * setupEqua3_frequencyIndependentPart --
 Int (G0 rho / epsilon0) = 4*PI*phi
 This equation is not used for MQS analysis
 **********************************************************************/
void
Formulation::setupEqua3_frequencyIndependentPart (
						  void)
{
  progressReport("setup kernel independent part of pfft for the union of all conductor surface");

  size_t directStencilSize, projectStencilSize, interpStencilSize;
#ifdef HIGH_ACCURACY
  directStencilSize = 6;
  projectStencilSize = 2;
  interpStencilSize = 2;
#else
  directStencilSize = 3;
  projectStencilSize = 1;
  interpStencilSize = 1;
#endif

  if (simuType_ == FULL_WAVE) {
    calcp3_FW = Calcp3_FW(outerOperatorList3, innerOperatorList3, 0);
#ifdef DEBUG_RHO
    pfft3_FW = Pfft3_FW(allElementList, allNonContactElementList, calcp3_FW,
			directStencilSize, projectStencilSize, interpStencilSize);
#else
    pfft3_FW = Pfft3_FW(allElementList, allElementList, calcp3_FW,
			directStencilSize, projectStencilSize, interpStencilSize);
#endif
  } else {
    calcp3_EMQS = Calcp3_EMQS(outerOperatorList3, innerOperatorList3);
#ifdef DEBUG_RHO
    pfft3_EMQS = Pfft3_EMQS(allElementList, allNonContactElementList, calcp3_EMQS,
			    directStencilSize, projectStencilSize, interpStencilSize);
#else
    pfft3_EMQS = Pfft3_EMQS(allElementList, allElementList, calcp3_EMQS,
			    directStencilSize, projectStencilSize, interpStencilSize);
#endif
  }
}

/**********************************************************************
 * setupEqua3_frequencyDependentPart --
 Int (G0 rho / epsilon0) = 4*PI*phi
 This equation is not used for MQS analysis
 **********************************************************************/
void
Formulation::setupEqua3_frequencyDependentPart (
						void)
{
  progressReport("setup kernel dependent part of pfft for the union of all conductor surface");

  TNT::stopwatch time;  time.start();
  
#ifdef DEBUG_MEMORY
  int numQudraturePoint = 0;
  bool skipPanelIntegration = true;
#else
  int numQudraturePoint = 24;
  bool skipPanelIntegration = false;
#endif
  if (simuType_ == FULL_WAVE) {
    kernel3_FW = Kernel3_FW(K0);
    calcp3_FW = Calcp3_FW(outerOperatorList3, innerOperatorList3, K0, numQudraturePoint);
    pfft3_FW.constructKernelDependent(kernel3_FW, calcp3_FW);
  } else {
    calcp3_EMQS = Calcp3_EMQS(outerOperatorList3, innerOperatorList3, skipPanelIntegration);
    pfft3_EMQS.constructKernelDependent(kernel3_EMQS, calcp3_EMQS);
  }
  timeReport("time for the pfft setup := ", time.read());
}

/**********************************************************************
 * fillIntegralEqua --
 **********************************************************************/
void
Formulation::fillIntegralEqua (
			       void)
{
  for (size_t i = 0; i < numCond_; i++) {
    calcp1List[i] = Calcp1(outerOperatorList1, innerOperatorList1, KC[i]);
  }

#ifdef DEBUG_MEMORY
  int numQudraturePoint = 0;
  bool skipPanelIntegration = true;
#else
  int numQudraturePoint = 24;
  bool skipPanelIntegration = false;
#endif

  if (simuType_ == FULL_WAVE) {
    calcp2_FW = Calcp2_FW(outerOperatorList2, innerOperatorList2, K0, numQudraturePoint);
    calcp3_FW = Calcp3_FW(outerOperatorList3, innerOperatorList3, K0, numQudraturePoint);
  } else if (simuType_ == EMQS) {
    calcp2_EMQS = Calcp2_EMQS(outerOperatorList2, innerOperatorList2, skipPanelIntegration);
    calcp3_EMQS = Calcp3_EMQS(outerOperatorList3, innerOperatorList3, skipPanelIntegration);
  } else if (simuType_ == MQS) {
    calcp2_EMQS = Calcp2_EMQS(outerOperatorList2, innerOperatorList2, skipPanelIntegration);
  }

  equa1Mat = SpRowMat<complex<double> >(3*totalNumE_, 6*totalNumE_);
  equa2aMat1 = SpRowMat<complex<double> >(condMesh.totalNumNonContactPanel(), 
					  6*totalNumE_);
  equa2aMat2 = SpRowMat<complex<double> >(condMesh.totalNumNonContactPanel(), 
					  6*totalNumE_);
  if (simuType_ != MQS) {
#ifdef DEBUG_RHO
    equa3Mat = SpRowMat<complex<double> >(condMesh.totalNumNonContactPanel(), totalNumUnknown_);
#else
    equa3Mat = SpRowMat<complex<double> >(totalNumE_, totalNumUnknown_);
#endif
  }

  size_t rowIndex_equa2a = 0;
  for (size_t tCondIndex=0; tCondIndex < numCond_; tCondIndex++) {
    for (size_t tPanelIndex=0; tPanelIndex < condMesh.numPanel(tCondIndex); 
	 tPanelIndex++) {
      for (size_t sCondIndex=0; sCondIndex < numCond_; sCondIndex++) {
	for (size_t sPanelIndex=0; sPanelIndex < condMesh.numPanel(sCondIndex); 
	     sPanelIndex++) {

	  if (sCondIndex == tCondIndex) {
	    if (formulationType_ == OLD) {
	      fillEqua_1_old(sCondIndex, sPanelIndex, tPanelIndex);
	    } else if (formulationType_ == NEW) {
	      fillEqua_1_new(sCondIndex, sPanelIndex, tPanelIndex);      
	    }
	  }
	  
	  complex<double> singleLayer;
	  if (condMesh.isNonContact(tCondIndex, tPanelIndex)) {
	    singleLayer = fillEqua_2a(sCondIndex, sPanelIndex, tCondIndex, 
				      tPanelIndex, rowIndex_equa2a);
	  }
	  
	  if (simuType_ != MQS) {
	    // Single layer potential is same as equation 2. So it is reused here.
#ifdef DEBUG_RHO
	    if (condMesh.isNonContact(tCondIndex, tPanelIndex)) {
	      fillEqua_3(sCondIndex, sPanelIndex, tCondIndex, tPanelIndex, singleLayer);
	    }
#else
	    fillEqua_3(sCondIndex, sPanelIndex, tCondIndex, tPanelIndex, singleLayer);
#endif
	  }

	}
      }

      // The index of equa2a here is not necessary the same as the panel 
      // index of which it is enforced. I have to compress it into 
      // two blocks of sparse matrices. Hence a separate index counter
      // is used here. The results of matrix vector product, however,
      // will be put into proper places in formulation.h:equa2a_directPart
      if (condMesh.isNonContact(tCondIndex, tPanelIndex)) {
	rowIndex_equa2a++;
      }

    }
  }

#ifdef DEBUG_ITERATIVE
  std::ofstream fout1("equa1.tmp");
  fout1 << equa1Mat;
  std::ofstream fout2("equa21.tmp");
  fout2 << equa2aMat1;
  std::ofstream fout3("equa22.tmp");
  fout3 << equa2aMat2;
  if (simuType_ != MQS) {
    std::ofstream fout4("equa3.tmp");
    fout4 << equa3Mat;
  }
#endif

}

/**********************************************************************
 * checkMemoryUsage --
 **********************************************************************/
void 
Formulation::checkMemoryUsage (
			       void) 
{
  const size_t MB = 1024*1024;
  size_t memoryUsedBySystem;
  if (solverType_ == ITERATIVE_WITH_PFFT) {
    memoryUsedBySystem = estimateMemoryUsage_pfft()/MB;
    memoryReport("pfft", memoryUsedBySystem);
  } else {
    memoryUsedBySystem = estimateMemoryUsage_systemMat()/MB;
    memoryReport("the System Matrix", memoryUsedBySystem);
  }

  size_t totalMemoryUsage;
  if (solverType_ != DIRECT) {
    const int numNonZeroEachRow = 15;
    size_t memoryUsedByPreCond = 
      equa1Mat.memoryEstimate(numNonZeroEachRow * totalNumUnknown_, 
			      totalNumUnknown_) / MB;
    memoryReport("the pre-conditioner Matrix", memoryUsedByPreCond);
    messageWithOneNumber("Assumed number of non-zero per row after LU factorization:= ", numNonZeroEachRow);
    
    int gmresRestart = getRestart(); // assume restart is used
    size_t numWorkVector = 4; // RHS, x, P and AP
    size_t memoryForGmres = (gmresRestart + numWorkVector) * 
      totalNumUnknown_ * sizeof(complex<double>) / MB;
    memoryReport("the Hessenberg Matrix in Gmres", memoryForGmres);
    messageWithOneNumber("Assumed restart in Gmres := ", gmresRestart);
  
    totalMemoryUsage = memoryUsedByPreCond + memoryUsedBySystem + memoryForGmres;
  } else {
    const size_t expansionFactor = 10;
    totalMemoryUsage = expansionFactor * memoryUsedBySystem;
    memoryReport("the LU factored system Matrix", totalMemoryUsage);
    messageWithOneNumber("Assume #nonzeros after LU factorization increases by a factor of ", 
			 expansionFactor);
  }
  memoryReport("the whole program", totalMemoryUsage);

  if (totalMemoryUsage > 2.*1024) {
    warningMessage("formulation.cc: checkMemoryUsage()",
		   "Estimated total memory usage has exceeded 2Gb!");
  }
}

/**********************************************************************
 * estimateMemoryUsage_pfft --
 * generate dummy pfft objects and use them to call memoryEstimate().
 * These objects will be override later by real data.
 **********************************************************************/
size_t
Formulation::estimateMemoryUsage_pfft (
				       void) 
{
  size_t memoryEstimate = 0;
  for (size_t i = 0; i < numCond_; i++) {
    memoryEstimate += pfft1List[i].memoryEstimate();
  }
  
  if (simuType_ == FULL_WAVE) {
    memoryEstimate += pfft2_FW.memoryEstimate();
  } else {
    memoryEstimate += pfft2_EMQS.memoryEstimate();
  }

  if (simuType_ == FULL_WAVE) {
    memoryEstimate += pfft3_FW.memoryEstimate();
  } else if (simuType_ == EMQS) {
    memoryEstimate += pfft3_EMQS.memoryEstimate();
  }

  return memoryEstimate;
}

/**********************************************************************
 * estimateMemoryUsage_systemMat --
 **********************************************************************/
size_t
Formulation::estimateMemoryUsage_systemMat (
					    void) const
{
  size_t memoryUsage = equa1Mat.memoryEstimate(2*totalNumE_*3*totalNumE_, 
					       3*totalNumE_);
  memoryUsage += 
    equa2aMat1.memoryEstimate(6*totalNumE_*condMesh.totalNumNonContactPanel(),
			      condMesh.totalNumNonContactPanel());
  memoryUsage += 
    equa2aMat2.memoryEstimate(6*totalNumE_*condMesh.totalNumNonContactPanel(),
			      condMesh.totalNumNonContactPanel());
  // equa_2a_sparse
  memoryUsage +=
    2*nonIntegralEquaMat.memoryEstimate(4*condMesh.totalNumNonContactPanel(),
					condMesh.totalNumNonContactPanel());
  // equa_2b
  memoryUsage +=
    2*nonIntegralEquaMat.memoryEstimate(3*condMesh.totalNumContactPanel(),
					condMesh.totalNumContactPanel());
  // equa_3
  if (simuType_ != MQS) {
    memoryUsage += equa3Mat.memoryEstimate(totalNumRho_*totalNumE_, totalNumE_);
  }

  // equa4
  memoryUsage += 
    nonIntegralEquaMat.memoryEstimate(3*condMesh.totalNumNonContactPanel(),
				      condMesh.totalNumNonContactPanel());
  // equa5
  memoryUsage += 
    nonIntegralEquaMat.memoryEstimate(3*condMesh.totalNumContactPanel(),
				      condMesh.totalNumContactPanel());
  // equa6
  memoryUsage += 
    nonIntegralEquaMat.memoryEstimate(24*condMesh.totalNumNonContactPanel(),
				      condMesh.totalNumNonContactPanel());
  // equa7
  memoryUsage += 
    nonIntegralEquaMat.memoryEstimate(1*condMesh.totalNumContactPanel(),
				      condMesh.totalNumContactPanel());

  return memoryUsage;
}

/**********************************************************************
 * fillEqua_1_old --
 *   Old:    Int (G1 dE/dn') - Int (E dG1/dn') = 4*PI*E 
 *   Enforced on each conductor seperately
 **********************************************************************/
void
Formulation::fillEqua_1_old (
			     const size_t condIndex,
			     const size_t sPanelIndex,
			     const size_t tPanelIndex)
{
  size_t sGlobalPanelIndex = condMesh.globalPanelIndex(condIndex, sPanelIndex);
  size_t tGlobalPanelIndex = condMesh.globalPanelIndex(condIndex, tPanelIndex);

  try {
    calcp1List[condIndex](allElementList[sGlobalPanelIndex], 
			  allElementList[tGlobalPanelIndex]);
  }
  catch (std::domain_error e) {
    cout << e.what() << endl;
    errorMessage("formulation.cc::fillEqua_1_old()",  "Bug in calcp");
  }

  complex<double> singleLayer = calcp1List[condIndex].result(0);
  complex<double> doubleLayer = calcp1List[condIndex].result(1);
  if (tPanelIndex == sPanelIndex) {
    // Int dG1/dn + 4*PI = 2 * PAI
    doubleLayer = 2 * PI;
  }
    
  for (size_t m = 0; m < 3; m++) {
    // m = 0:  x equation
    // m = 1:  y equation
    // m = 2:  z equation
    size_t rowIndex = tGlobalPanelIndex + m*totalNumE_;

    // Int G1
    size_t colIndex = sGlobalPanelIndex + globalStartColIndex.dEdn[m];
    equa1Mat.insertElement(rowIndex, colIndex, singleLayer);
    
    // -Int dG1/dn
    colIndex = sGlobalPanelIndex + globalStartColIndex.E[m];
    equa1Mat.insertElement(rowIndex, colIndex, -doubleLayer);
  }
}

/**********************************************************************
 * fillEqua_1_new --
 *   New:    d/dn { Int (G1 dE/dn') - Int (E dG1/dn') = 4*PI*E }
 *   Enforced on each conductor seperately
 **********************************************************************/
void
Formulation::fillEqua_1_new (
			     const size_t condIndex,
			     const size_t sPanelIndex,
			     const size_t tPanelIndex)
{
  size_t sGlobalPanelIndex = condMesh.globalPanelIndex(condIndex, sPanelIndex);
  size_t tGlobalPanelIndex = condMesh.globalPanelIndex(condIndex, tPanelIndex);

  calcp1List[condIndex](allElementList[sGlobalPanelIndex], 
			allElementList[tGlobalPanelIndex]);
  complex<double> dn_singleLayer = calcp1List[condIndex].result(0);
  complex<double> dn_doubleLayer = calcp1List[condIndex].result(1);

  if (tPanelIndex == sPanelIndex) {
    /* When target and source panel are same, n = n', d/dn = -d/dn'
       d/dn Int(G) = - Int (dG/dn') = - Intp (dG/dn') + 2*PAI, 
       where Intp means the singular point has been excluded from 
       the integration area. 
       Calcp shoule give me this result. But I want to make sure.
    */
    if (abs(real(dn_singleLayer)/2./PI - 1.) > 1e-3) {
      errorMessage("formulation.c:fillEqua_1_new",
		   "bug in calcp for hyper-singular case");
    }
    // d/dn Int G1 - 4*PI = -2*PI
    dn_singleLayer = -2*PI;
  }

  for (size_t m = 0; m < 3; m++) {
    // m = 0:  x equation
    // m = 1:  y equation
    // m = 2:  z equation
    size_t rowIndex = tGlobalPanelIndex + m*totalNumE_;

    // d/dn Int G1
    size_t colIndex = sGlobalPanelIndex + globalStartColIndex.dEdn[m];
    equa1Mat.insertElement(rowIndex, colIndex, dn_singleLayer);
    
    // -d/dn Int dG1/dn
    colIndex = sGlobalPanelIndex + globalStartColIndex.E[m];
    equa1Mat.insertElement(rowIndex, colIndex, -dn_doubleLayer);
  }
}

/**********************************************************************
 * fillEqua_2a --
 * t dot { Int (G0 dE/dn') - Int(E dG0/dn') } + t dot grad(phi) * 4 * PAI = 0
 * the equation is enforced only along two tangent directions
 * on the non-contact target panel 
 *
 * 4*PAI before phi is because calcp does not include 4*PAI in the denominator 
 * of the kernel exp(ikr)/r
 *
 * Note:
 * 1) The single layer potential is returned for reuse in equation 
 * Int G rho = phi
 **********************************************************************/
std::complex<double>
Formulation::fillEqua_2a (
			  const size_t sCondIndex,
			  const size_t sPanelIndex,
			  const size_t tCondIndex,
			  const size_t tPanelIndex,
			  const size_t rowIndex)
{
  size_t sGlobalPanelIndex = condMesh.globalPanelIndex(sCondIndex, 
						       sPanelIndex);
  size_t tGlobalPanelIndex = condMesh.globalPanelIndex(tCondIndex, 
						       tPanelIndex);

  try {
    if (simuType_ == FULL_WAVE) {
      calcp2_FW(allElementList[sGlobalPanelIndex], allElementList[tGlobalPanelIndex]);
    } else {
      calcp2_EMQS(allElementList[sGlobalPanelIndex], allElementList[tGlobalPanelIndex]);
    }
  }
  catch (std::domain_error e) {
    cout << e.what() << endl;
  }

  // the best way
  // TemplateCalcp2::ValueType singleLayer = calcp2.result(0);
  complex<double> singleLayer;
  complex<double> doubleLayer;
  
  if (simuType_ == FULL_WAVE) {
    singleLayer = calcp2_FW.result(0);
    doubleLayer = calcp2_FW.result(1);
  } else {
    singleLayer = calcp2_EMQS.result(0);
    doubleLayer = calcp2_EMQS.result(1);
  }

  for (size_t m = 0; m < 3; m++) {
    size_t colIndex1 = sGlobalPanelIndex + globalStartColIndex.dEdn[m];
    size_t colIndex2 = sGlobalPanelIndex + globalStartColIndex.E[m];

    // t1 Int G0 
    equa2aMat1.insertElement(rowIndex, colIndex1, 
			     singleLayer * allElementList[tGlobalPanelIndex].tangent1(m));
    // -t1 Int dG0/dn
    equa2aMat1.insertElement(rowIndex, colIndex2, 
			     -doubleLayer * allElementList[tGlobalPanelIndex].tangent1(m));

    // t2 Int G0 
    equa2aMat2.insertElement(rowIndex, colIndex1, 
			     singleLayer * allElementList[tGlobalPanelIndex].tangent2(m));    
    // -t2 Int dG0/dn
    equa2aMat2.insertElement(rowIndex, colIndex2, 
			     -doubleLayer * allElementList[tGlobalPanelIndex].tangent2(m));
  }

  return singleLayer;
}

/**********************************************************************
 * fillEqua_3 --
 **********************************************************************/
void
Formulation::fillEqua_3 (
			 const size_t sCondIndex,
			 const size_t sPanelIndex,
			 const size_t tCondIndex,
			 const size_t tPanelIndex,
			 std::complex<double>& singleLayer)
{
  size_t sGlobalPanelIndex = condMesh.globalPanelIndex(sCondIndex, 
						       sPanelIndex);
  size_t tGlobalPanelIndex = condMesh.globalPanelIndex(tCondIndex, 
						       tPanelIndex);

  if (condMesh.isContact(tCondIndex, tPanelIndex)) {
    // then we could not reuse the single-layer computed in fillEqua2a()
    if (simuType_ == FULL_WAVE) {
      calcp3_FW(allElementList[sGlobalPanelIndex], allElementList[tGlobalPanelIndex]);
      singleLayer = calcp3_FW.result(0);
    } else {
      calcp3_EMQS(allElementList[sGlobalPanelIndex], allElementList[tGlobalPanelIndex]);
      singleLayer = calcp3_EMQS.result(0);
    }
  }
    
  size_t rowIndex = tGlobalPanelIndex;
  size_t colIndex = globalStartColIndex.rho + sGlobalPanelIndex;
  equa3Mat.insertElement(rowIndex, colIndex, singleLayer);	
}

/**********************************************************************
 * setupNonIntegralEqua --
 * equation 2a:    4*PI* t dot grad Phi term   on the non-contact panels
 * equation 2b:    t dot E = 0               on the contact panels
 * equation 3:     -4*PI*average * phi term  on all panels
 * equation 4                                on non-contact panels
 *  n dot E = j*2*PI*freq*epsilon0/sigma * rho    EMQS or full-wave analysis
 *  n dot E = 0                                   MQS analysis
 *  note:  rho = rho / epsilon0, to be consistant with equa3
 * equation 5     dEn/dn = 0               on contact panels.	 
 * equation 6     div E = 0     on the vertices of non-contact panels
 * equation 7     phi = constant  on vertices of contact panels.
 * 
 * Without or With pFFT accelaration,
 * these equations are treated as a single sparse matrix here. 
 * This matrix is frequency independent. So it should 
 * be setup only once. The frequency term in equa4 is handled in
 * function nonIntegralEqua() in formulation.h.
 **********************************************************************/
void
Formulation::setupNonIntegralEqua (
				   void)
{
  nonIntegralEquaMat = SpRowMat<double>(numRowOfNonIntegralMat_, 
					totalNumUnknown_);

  // pFFT takes care of all integral operators. The remaining equations 
  // are represented by a single sparse matrix. So the start row index
  // for each equation should be adjusted accordingly.
  size_t rowIndex_equa2 = 0;
  size_t rowIndex_equa3 = globalStartRowIndex.equa_3 - globalStartRowIndex.equa_2;
  size_t rowIndex_equa45 = globalStartRowIndex.equa_4_5 - globalStartRowIndex.equa_2;
  size_t rowIndex_equa67 = globalStartRowIndex.equa_6_7 - globalStartRowIndex.equa_2;

  for (size_t condIndex = 0; condIndex < condMesh.numCond(); condIndex++) {
    for (size_t panelIndex = 0; panelIndex < condMesh.numPanel(condIndex); 
	 panelIndex++) {
      size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);

      if (condMesh.isNonContact(condIndex, panelIndex)) {
	fillEqua_2a_sparsePart(globalPanelIndex, condIndex, panelIndex, 
			       rowIndex_equa2);	
      } else {
	fillEqua_2b(globalPanelIndex, rowIndex_equa2);
      }
      rowIndex_equa2 ++;

      if (simuType_ != MQS) {
#ifdef DEBUG_RHO
	if (condMesh.isNonContact(condIndex, panelIndex)) {
	  fillEqua_3_sparsePart(condIndex, panelIndex, rowIndex_equa3);
	} else {
	  fillEqua_3_contact(globalPanelIndex, rowIndex_equa3);
	}
#else
	fillEqua_3_sparsePart(condIndex, panelIndex, rowIndex_equa3);
#endif
	rowIndex_equa3 ++;
      }

      // The row indices of equa 4 and 5 are mixed for my convenience.
      if (condMesh.isNonContact(condIndex, panelIndex)) {
	fillEqua_4(globalPanelIndex, rowIndex_equa45);
      } else {
	fillEqua_5(globalPanelIndex, rowIndex_equa45);
      }
      rowIndex_equa45 ++;
    }
  }

  // The row indices of equa 6 and 7 are mixed for my convenience.
  for (size_t globalNodeIndex = 0; globalNodeIndex < condMesh.totalNumNode(); 
       globalNodeIndex++) {
    if ( (condMesh.vertexType(globalNodeIndex) == NON_CONTACT) ||
	 (condMesh.vertexType(globalNodeIndex) == BUFFER) ) {
      fillEqua_6(globalNodeIndex, rowIndex_equa67);      
    } else {
      fillEqua_7(globalNodeIndex, rowIndex_equa67);
    }
    rowIndex_equa67++;
  }

}

/**********************************************************************
 * fillEqua_2a_sparsePart --
 * The sparse part of the equa_2a, the 4*PI*grad (phi) term. 
 **********************************************************************/
void 
Formulation::fillEqua_2a_sparsePart (
				     const size_t globalPanelIndex, 
				     const size_t condIndex, 
				     const size_t panelIndex, 
				     const size_t rowIndex)
{
  if (allElementList[globalPanelIndex].shape() == 4) {    
    compGradPhiForQuad(globalPanelIndex, condIndex, panelIndex, rowIndex);
  } else if (allElementList[globalPanelIndex].shape() == 3) {
    compGradPhiForTriangle(globalPanelIndex, condIndex, panelIndex, rowIndex);
  } else {
    errorMessage("formulation.cc::fillEqua_2a_sparsePart", 
		 "Illegal panel shape!");
  }
}

/**********************************************************************
 * compGradPhiForTriangle --
 * phi(t1, t2) = phi1*L1(t1, t2) + phi2*L2(t1, t2) + phi3*L3(t1, t2)
 * where L1, L2 and L3 are the linear shape functions defined as 
 * L_j(t1, t2) = a_j * t1 + b_j * t2 + c
 * t1 dot grad(phi) = d/dt1 phi = sum_j {a_j phi_j}, j = 1...3
 * t2 dot grad(phi) = d/dt2 phi = sum_j {b_j phi_j}, j = 1...3
 **********************************************************************/
void 
Formulation::compGradPhiForTriangle (
				     const size_t globalPanelIndex, 
				     const size_t condIndex, 
				     const size_t panelIndex, 
				     const size_t rowIndex)
{
  for (size_t vertexIndex = 0; vertexIndex < 3; vertexIndex++) {
    size_t globalNodeIndex = condMesh.globalNodeIndex(condIndex, panelIndex, vertexIndex);
    size_t col = globalStartColIndex.phi + globalNodeIndex;

    // t1 dot grad(phi) = d/dt1 phi = sum_j {a_j phi_j}, j = 1...3
    double a = allElementList[globalPanelIndex].linearShapeFuncCoe(vertexIndex,  0);
    nonIntegralEquaMat.insertElement(rowIndex, col, 4*PI*a);

    // t2 dot grad(phi) = d/dt2 phi = sum_j {b_j phi_j}, j = 1...3
    double b = allElementList[globalPanelIndex].linearShapeFuncCoe(vertexIndex,  1);
    nonIntegralEquaMat.insertElement(rowIndex + totalNumE_, col, 4*PI*b);
  }
}

/**********************************************************************
 * compGradPhiForQuad --
 * Finite Difference scheme is used here 
 * Assumption: 
 * 1)  t1 and t2 are aligned with dist1 and dist2
 * This is generally true if the panels are rectangualr or trapezoidal.
 * It is WRONG if the panels are irregular. In this case, the quadratic 
 * interpolation should be used, just as in FEM.
 * 2)  The shape of panel is quad.
 **********************************************************************/
void 
Formulation::compGradPhiForQuad (
				 const size_t globalPanelIndex, 
				 const size_t condIndex, 
				 const size_t panelIndex, 
				 const size_t rowIndex)
{
  // NOTE:
  //  The sign here should be consistent with the sign defined in 
  //  pfft::element.cc: compNormalAndTangentVector() 
  //
  // t1 dot 4*PI*grad(phi) = 4*PI*(phi_2 + phi_3 - phi_1 - phi_0) / (2 * dist1)
  double coe1 = 4*PI/(2. * allElementList[globalPanelIndex].edgeMidPointDist1());
  vector<double> FD1(4);
  FD1[0] = -coe1;
  FD1[1] = -coe1;
  FD1[2] = coe1;
  FD1[3] = coe1;

  // t2 dot 4*PI*grad(phi) = 4*PI*(-phi_0 - phi_3 + phi_1 + phi_2) / (2 * dist2) 
  double coe2 = 4*PI/(2. * allElementList[globalPanelIndex].edgeMidPointDist2());
  vector<double> FD2(4);
  FD2[0] = -coe2;
  FD2[1] = +coe2;
  FD2[2] = +coe2;
  FD2[3] = -coe2;
  
  for (size_t i = 0; i < condMesh.numVertex(condIndex, panelIndex); i++) {
    size_t globalNodeIndex = condMesh.globalNodeIndex(condIndex, panelIndex, i);
    size_t col = globalStartColIndex.phi + globalNodeIndex;
    nonIntegralEquaMat.insertElement(rowIndex, col, FD1[i]);
    nonIntegralEquaMat.insertElement(rowIndex+totalNumE_, col, FD2[i]);
  }
}

/**********************************************************************
 * fillEqua_2b --
 * equation 2b t dot E = 0               on the contact panels
 **********************************************************************/
void
Formulation::fillEqua_2b (
			  const size_t globalPanelIndex, 
			  const size_t rowIndex)
{
  for (size_t m=0; m<3; m++) {
    /* m = 0: Ex column 
       m = 1: Ey column
       m = 2: Ez column
    */
    size_t colIndex = globalPanelIndex + globalStartColIndex.E[m];

    // equation for t1
    nonIntegralEquaMat.insertElement(rowIndex, colIndex, 
				     allElementList[globalPanelIndex].tangent1(m));

    // equation for t2
    nonIntegralEquaMat.insertElement(rowIndex+totalNumE_, colIndex,
				     allElementList[globalPanelIndex].tangent2(m));
  }
}

/**********************************************************************
 * fillEqua_3_sparsePart --
 * The sparse part of equa_3, the -4*PI * average(phi) term
 **********************************************************************/
void 
Formulation::fillEqua_3_sparsePart (
				    const size_t condIndex, 
				    const size_t panelIndex, 
				    const size_t rowIndex)
{
  double value = -4.*PI / condMesh.numVertex(condIndex, panelIndex);
  for (size_t li = 0; li < condMesh.numVertex(condIndex, panelIndex); li++) {
    size_t colIndex = globalStartColIndex.phi + 
      condMesh.globalNodeIndex(condIndex, panelIndex, li);
    nonIntegralEquaMat.insertElement(rowIndex, colIndex, value);
  }
}

#ifdef DEBUG_RHO
/**********************************************************************
 * fillEqua_3_contact --
 * rho = 0               on the contact panels
 **********************************************************************/
void
Formulation::fillEqua_3_contact (
				 const size_t globalPanelIndex, 
				 const size_t rowIndex)
{
  size_t colIndex = globalPanelIndex + globalStartColIndex.rho;
  nonIntegralEquaMat.insertElement(rowIndex, colIndex, 1.);
}
#endif

/**********************************************************************
 * fillEqua_4 --
 * equation 4                                on non-contact panels
 *  n dot E = j*2*PI*freq*epsilon0/sigma * rho    EMQS or full-wave analysis
 *  n dot E = 0                                   MQS analysis
 * Note: 
 * The rho is not taken care of here. So, the equation is simply:
 *         n dot E = 0                                   
 **********************************************************************/
void
Formulation::fillEqua_4 (
			 const size_t globalPanelIndex, 
			 const size_t rowIndex)
{
  for (size_t m=0; m<3; m++) {
    /* m = 0: Ex column 
       m = 1: Ey column
       m = 2: Ez column
    */  
    size_t col = globalStartColIndex.E[m] + globalPanelIndex;
    nonIntegralEquaMat.insertElement(rowIndex, col, 
				     allElementList[globalPanelIndex].normal(m));
  }
}

/**********************************************************************
 * fillEqua_5 --
 * n dot dE/dn = 0
 **********************************************************************/
void
Formulation::fillEqua_5 (
			 const size_t globalPanelIndex, 
			 const size_t rowIndex)
{
  for (size_t m=0; m<3; m++) {
    /* m = 0: dEx/dn column 
       m = 1: dEy/dn column
       m = 2: dEz/dn column
    */  
    size_t col = globalStartColIndex.dEdn[m] + globalPanelIndex;
    nonIntegralEquaMat.insertElement(rowIndex, col, 
				     allElementList[globalPanelIndex].normal(m));
  }
}

/**********************************************************************
 * fillEqua_6 --
 * div E = 0     on the vertices of non-contact panels
 * Int Et dot (n cross l) dl - Int dE/dn ds = 0
 * The equation is enforced on the dual panel around each vertex. 
 * So the number of equations in this block is equal to the number 
 * of unknown scalar potential phi 
 *
 * Note:
 * To balance magnitude of the coefficients and that of single-layer potential
 * such that the pivoting in LU factorization of pre-conditioner will not
 * cause row interchanges and hence will not introduce larger number of 
 * fill-in's, I multiply this equation with the normalization factor.
 * And since it is a homogeneous equation, this manipulation won't change
 * the result.
 **********************************************************************/
void
Formulation::fillEqua_6 (
			 const size_t dualPanelIndex, 
			 const size_t rowIndex)
{
  size_t numSubPanel = dualPanelList[dualPanelIndex].numSubPanel();
  for (size_t subPanelIndex=0; subPanelIndex<numSubPanel; subPanelIndex++) {
    size_t parentPanelGlobalIndex = 
      dualPanelList[dualPanelIndex].parentPanelGlobalIndex(subPanelIndex);

    vector3D<double> normal = allElementList[parentPanelGlobalIndex].normal();

    point3D node0 = dualPanelList[dualPanelIndex].subPanelVertex(subPanelIndex, 0);
    point3D node1 = dualPanelList[dualPanelIndex].subPanelVertex(subPanelIndex, 1);
    point3D node2 = dualPanelList[dualPanelIndex].subPanelVertex(subPanelIndex, 2);
    vector3D<double> e1 = node0 - node1;
    vector3D<double> e2 = node1 - node2;
    vector3D<double> en1 = crossProd(normal, e1); // without normalization, the desired
    vector3D<double> en2 = crossProd(normal, e2); // edge length is included here.
    double cx = (en1.x() + en2.x()) / surfConst_.normFactor;
    double cy = (en1.y() + en2.y()) / surfConst_.normFactor;
    double cz = (en1.z() + en2.z()) / surfConst_.normFactor;

    /* Ex */
    size_t col = globalStartColIndex.E[0] + parentPanelGlobalIndex;
    nonIntegralEquaMat.insertElement(rowIndex, col, cx);

    /* Ey */
    col = globalStartColIndex.E[1] + parentPanelGlobalIndex;
    nonIntegralEquaMat.insertElement(rowIndex, col, cy);

    /* Ez */
    col = globalStartColIndex.E[2] + parentPanelGlobalIndex;
    nonIntegralEquaMat.insertElement(rowIndex, col, cz);

    double subPanelArea = 
      dualPanelList[dualPanelIndex].subPanelArea(subPanelIndex) / surfConst_.normFactor;
    double dcx = subPanelArea * normal.x();
    double dcy = subPanelArea * normal.y();
    double dcz = subPanelArea * normal.z();

    /* dEx/dn */
    col = globalStartColIndex.dEdn[0] + parentPanelGlobalIndex;
    nonIntegralEquaMat.insertElement(rowIndex, col, -dcx);

    /* dEy/dn */
    col = globalStartColIndex.dEdn[1] + parentPanelGlobalIndex;
    nonIntegralEquaMat.insertElement(rowIndex, col, -dcy);

    /* dEz/dn */
    col = globalStartColIndex.dEdn[2] + parentPanelGlobalIndex;
    nonIntegralEquaMat.insertElement(rowIndex, col, -dcz);
  }

}

/**********************************************************************
 * fillEqua_7 --
 * phi = constant  
 * This equation is only enforced on vertices of contact panels.
 **********************************************************************/
void
Formulation::fillEqua_7 (
			 const size_t globalNodeIndex, 
			 const size_t rowIndex)
{
  size_t col = globalStartColIndex.phi + globalNodeIndex;
  nonIntegralEquaMat.insertElement(rowIndex, col, 1.);
}

/**********************************************************************
 * setupElementList --
 * One element list for each conductor, and then one long element list
 * for the union of all conductors
 * Note: for one conductor structure, the longer list is redundant.
 **********************************************************************/
void
Formulation::setupElementList(
			      void)
{
  condElementList.resize(numCond_);

  for (size_t condIndex = 0; condIndex < numCond_; condIndex++) {
    for (size_t eleIndex = 0; eleIndex < condMesh.numPanel(condIndex); 
	 eleIndex++) {      
      condElementList[condIndex].push_back(element(condMesh.vertex(condIndex, 
								   eleIndex), 
						   eleIndex, condIndex));

    }
  }

  minElementSizeList_.resize(numCond_);
  for (size_t condIndex = 0; condIndex < numCond_; condIndex++) {
    minElementSizeList_[condIndex] = DBL_MAX;
    for (size_t eleIndex = 0; eleIndex < condMesh.numPanel(condIndex); 
	 eleIndex++) {      
      allElementList.push_back(condElementList[condIndex][eleIndex]);
      minElementSizeList_[condIndex] = 
	min(minElementSizeList_[condIndex], 
	    2.*condElementList[condIndex][eleIndex].boundingSphereRadius());
    }
  }

  minElementSize_ = DBL_MAX;
  for (size_t condIndex = 0; condIndex < numCond_; condIndex++) {
    minElementSize_ = min(minElementSize_, minElementSizeList_[condIndex]);
  }

  maxElementSize_ = 0.;
  for (size_t i = 0; i < allElementList.size(); i++) {
    maxElementSize_ = max(maxElementSize_, 
			  2.*allElementList[i].boundingSphereRadius());
  }

  // equa_2 is enforced only on non-contact panels, so I have to construct 
  // a separate list of panels here.
  for (size_t condIndex = 0; condIndex < numCond_; condIndex++) {
    for (size_t eleIndex = 0; eleIndex < condMesh.numPanel(condIndex); 
	 eleIndex++) {      
      if (condMesh.isNonContact(condIndex, eleIndex)) {
	allNonContactElementList.push_back(condElementList[condIndex][eleIndex]);
      }
    }
  }

  //  cout << endl << "\t The element list has been generated for " 
  //       << numCond_ << " conductor(s)" << endl;

}

/**********************************************************************
 * compWaveNumber --
 **********************************************************************/
void
Formulation::compWaveNumber (
			     const double f,
			     double& K0,
			     std::vector<std::complex<double> >& KC) const
{
  if (simuType_ == FULL_WAVE) {
    /*   k^2 = w^2*epsilon0*mu  */
    K0 = 2. * PI * f * sqrt(surfConst_.EPSILON0 * surfConst_.MU0);
  } else {
    K0 = 0.;
  }

  KC.clear();
  for (size_t i=0; i < numCond_; i++) {
    double sigma = condInfoPtrList[i]->conductivity();
    // KC^2 = w^2*epsilon0*mu - j*w*mu*sigma = K^2 - j*w*mu*sigma
    complex<double> kc(K0*K0, - 2. * PI * f * surfConst_.MU0 * sigma);
    kc = sqrt(kc);
    if (real(kc) > 0) {
      kc *= -1;
    }
    KC.push_back(kc);
  }
}

/**********************************************************************
 * allocateWorkUnit --
 **********************************************************************/
void
Formulation::allocateWorkUnit (
			       void)
{
  w1.resize(condMesh.totalNumPanel());
  w2.resize(condMesh.totalNumPanel());
  xw.resize(condMesh.totalNumPanel());
}

/**********************************************************************
 * setupGlobalStartColIndex --
 **********************************************************************/
void
Formulation::setupGlobalStartColIndex (
				       void)
{
  if (formulationType_ == OLD) {
    // make sure the double-layer is on the diagonal because it is always 2*PI
    // This option is better because it is more 
    // diagonally dominant, hence less row interchanges if partial 
    // pivoting is used for LU factoring the pre-conditioner.
    globalStartColIndex.E[0] = 0 * totalNumE_;
    globalStartColIndex.E[1] = 1 * totalNumE_;
    globalStartColIndex.E[2] = 2 * totalNumE_;
    globalStartColIndex.dEdn[0] = 3 * totalNumE_;
    globalStartColIndex.dEdn[1] = 4 * totalNumE_;
    globalStartColIndex.dEdn[2] = 5 * totalNumE_;
  } else if (formulationType_ == NEW) {
    // make sure the d/dn single-layer is on the diagonal because it is always 2*PI
    globalStartColIndex.dEdn[0] = 0 * totalNumE_;
    globalStartColIndex.dEdn[1] = 1 * totalNumE_;
    globalStartColIndex.dEdn[2] = 2 * totalNumE_;
    globalStartColIndex.E[0] = 3 * totalNumE_;
    globalStartColIndex.E[1] = 4 * totalNumE_;
    globalStartColIndex.E[2] = 5 * totalNumE_;
  }
  globalStartColIndex.phi = 6 * totalNumE_;
  globalStartColIndex.rho = 6 * totalNumE_ + totalNumPhi_;

  totalNumUnknown_ = 6 * totalNumE_ + totalNumPhi_;
  if (simuType_ != MQS) 
    totalNumUnknown_ += totalNumRho_;
}

/**********************************************************************
 * setupGlobalStartRowIndex --
 **********************************************************************/
void
Formulation::setupGlobalStartRowIndex (
				       void)
{
  int totalNumPanel = condMesh.totalNumPanel();
  int totalNumNonContactPanel = condMesh.totalNumNonContactPanel();
  int totalNumContactPanel = totalNumPanel - totalNumNonContactPanel;
  int totalNumNode = condMesh.totalNumNode();

  totalNumRow_ = 0;

  globalStartRowIndex.equa_1 = totalNumRow_;
  int numRowInThisBlock = 3*totalNumPanel;
  totalNumRow_ += numRowInThisBlock;

  globalStartRowIndex.equa_2 = totalNumRow_;
  numRowInThisBlock = 2*totalNumPanel;
  totalNumRow_ += numRowInThisBlock;

  if (simuType_ != MQS) {
    globalStartRowIndex.equa_3 = totalNumRow_;
    numRowInThisBlock = totalNumPanel;
    totalNumRow_ += numRowInThisBlock;
  }

  globalStartRowIndex.equa_4_5 = totalNumRow_;
  numRowInThisBlock = totalNumPanel;
  totalNumRow_ += numRowInThisBlock;

  globalStartRowIndex.equa_6_7 = totalNumRow_;
  numRowInThisBlock = totalNumNode; 
  totalNumRow_ += numRowInThisBlock;

  numRowOfNonIntegralMat_ = totalNumRow_ - globalStartRowIndex.equa_2;
  globalStartRowIndex.nonIntegralEqua = globalStartRowIndex.equa_2;
}

/**********************************************************************
 * setupDualPanelList --
 **********************************************************************/
void
Formulation::setupDualPanelList (
				 void)
{
  dualPanelList.resize(condMesh.totalNumNode());
  vector<point3D> node(3);

  for (size_t condIndex = 0; condIndex < condMesh.numCond(); condIndex++) {
    for (size_t localNodeIndex = 0; localNodeIndex < condMesh.numNode(condIndex); 
	 localNodeIndex++) {
      size_t centerNodeIndex = condMesh.globalNodeIndex(condIndex, localNodeIndex);
      point3D center = condMesh.vertexPos(centerNodeIndex);
      size_t numSubPanel = condMesh.numSharedPanel(centerNodeIndex);
      
      for (size_t spi = 0; spi < numSubPanel; spi++) {
	size_t parentPanelLocalIndex = 
	  condMesh.sharedPanelIndex(centerNodeIndex, spi);
	size_t parentPanelGlobalIndex = 
	  condMesh.globalPanelIndex(condIndex, parentPanelLocalIndex);
	node[1] = allElementList[parentPanelGlobalIndex].centroid();
	int numVertex = allElementList[parentPanelGlobalIndex].numVertex();
	for (size_t i = 0; i < numVertex; i++) {
	  if (condMesh.globalNodeIndex(condIndex, parentPanelLocalIndex, i) == 
	      centerNodeIndex) {
	    int j = (i - 1 + numVertex) % numVertex;
	    node[2] = 0.5 * (center + 
			     allElementList[parentPanelGlobalIndex].vertex(j));
	  
	    j = (i + 1) % numVertex;
	    node[0] = 0.5 * (center + 
			     allElementList[parentPanelGlobalIndex].vertex(j));
	  }
	}
	double area = compTriangleArea(center, node[0], node[1]) + 
	  compTriangleArea(center, node[1], node[2]);

	dualPanelList[centerNodeIndex].addSubPanel(SubPanel(area, 
							    parentPanelGlobalIndex, node));
      }
    }
  }
}

/**********************************************************************
 * compTriangleArea --
 **********************************************************************/
double
Formulation::compTriangleArea (
			       const point3D& node1,
			       const point3D& node2,
			       const point3D& node3) const
{
  vector3D<double> edge1 = node2 - node1;
  vector3D<double> edge2 = node2 - node3;
  vector3D<double> vec = crossProd(edge1, edge2);
  return 0.5 * length(vec);
}

/**********************************************************************
 * memoryReport --
 **********************************************************************/
void
Formulation::memoryReport (
			   const std::string object,
			   const size_t memoryUsage) const
{
#ifdef PRINT_MEMORY
  std::streamsize prec = std::cout.precision();
  std::cout << setprecision(3) << endl 
	    << "\tEstimated memory usage by " << object 
	    <<  " := "<< memoryUsage << " (MB)"
	    << setprecision(prec) << std::endl;
#endif
}

/**********************************************************************
 * setupLocalStartColIndex --
 * the sequence of unknowns in local coordinate system
 **********************************************************************/
void
Formulation::setupLocalStartColIndex (
				      void)
{
  if (formulationType_ == OLD) {
    // make sure the double-layer is on the diagonal because it is always 2*PI
    localStartColIndex.E[0] = 0; // Et1
    localStartColIndex.E[1] = totalNumE_; // Et2
    localStartColIndex.E[2] = 2 * totalNumE_; // En
    localStartColIndex.dEdn[0] = 3 * totalNumE_; // dEt1
    localStartColIndex.dEdn[1] = 4 * totalNumE_; // dEt2
    localStartColIndex.dEdn[2] = 5 * totalNumE_; // dEn
  } else if (formulationType_ == NEW) {
    // make sure the d/dn single-layer is on the diagonal because it is always 2*PI
    localStartColIndex.dEdn[0] = 0 * totalNumE_; // dEt1
    localStartColIndex.dEdn[1] = 1 * totalNumE_; // dEt2
    localStartColIndex.dEdn[2] = 2 * totalNumE_; // dEn
    localStartColIndex.E[0] = 3 * totalNumE_; // Et1
    localStartColIndex.E[1] = 4 * totalNumE_; // Et2
    localStartColIndex.E[2] = 5 * totalNumE_; // En
  }
  localStartColIndex.phi = 6 * totalNumE_;
  localStartColIndex.rho = 6 * totalNumE_ + totalNumPhi_;
}

/**********************************************************************
 * setupLocalStartRowIndex --
 * the order of equations enforced in the local coordinate system
 **********************************************************************/
void
Formulation::setupLocalStartRowIndex (
				      void)
{
  int totalNumPanel = condMesh.totalNumPanel();
  int totalNumNode = condMesh.totalNumNode();
  int totalNumRow = 0;

  localStartRowIndex.equa_1 = totalNumRow;
  int numRowInThisBlock = 3*totalNumPanel;
  totalNumRow += numRowInThisBlock;

  localStartRowIndex.equa_2 = totalNumRow;
  numRowInThisBlock = 2*totalNumPanel;
  totalNumRow += numRowInThisBlock;

  // putting equa3 before 45 and 67 turns out to be better
  // because it generates fewer fill-in's
  if (simuType_ != MQS) {
    localStartRowIndex.equa_3 = totalNumRow;
    numRowInThisBlock = totalNumPanel;
    totalNumRow += numRowInThisBlock;
  }

  localStartRowIndex.equa_4_5 = totalNumRow;
  numRowInThisBlock = totalNumPanel;
  totalNumRow += numRowInThisBlock;

  localStartRowIndex.equa_6_7 = totalNumRow;
  numRowInThisBlock = totalNumNode; 
  totalNumRow += numRowInThisBlock;

  localStartRowIndex.nonIntegralEqua = localStartRowIndex.equa_2;
}

/**********************************************************************
 * setupLocalToGlobalMat --
 **********************************************************************/
void
Formulation::setupLocalToGlobalMat (
				    void)
{
  localToGlobalMat = SpRowMat<double>(totalNumUnknown_, totalNumUnknown_);

  size_t numElement = allElementList.size();
  for (size_t elementIndex = 0; elementIndex < numElement; elementIndex++) {
    for (size_t m = 0; m < 3; m++) {
      // Et1, Et2, En to Ex, Ey, Ez
      size_t rowIndex = elementIndex + globalStartColIndex.E[m];
      size_t colIndex = elementIndex + localStartColIndex.E[0];
      localToGlobalMat.insertElement(rowIndex, colIndex, 
				     allElementList[elementIndex].tangent1(m));
      colIndex = elementIndex + localStartColIndex.E[1];
      localToGlobalMat.insertElement(rowIndex, colIndex, 
				     allElementList[elementIndex].tangent2(m));
      colIndex = elementIndex + localStartColIndex.E[2];
      localToGlobalMat.insertElement(rowIndex, colIndex, 
				     allElementList[elementIndex].normal(m));

      // dEt1, dEt2, dEn to dEx, dEy, dEz
      rowIndex = elementIndex + globalStartColIndex.dEdn[m];
      colIndex = elementIndex + localStartColIndex.dEdn[0];
      localToGlobalMat.insertElement(rowIndex, colIndex, 
				     allElementList[elementIndex].tangent1(m));
      colIndex = elementIndex + localStartColIndex.dEdn[1];
      localToGlobalMat.insertElement(rowIndex, colIndex, 
				     allElementList[elementIndex].tangent2(m));
      colIndex = elementIndex + localStartColIndex.dEdn[2];
      localToGlobalMat.insertElement(rowIndex, colIndex, 
				     allElementList[elementIndex].normal(m));
    }
    
    // rho to rho
    if (simuType_ != MQS) {
      size_t rowIndex = elementIndex + globalStartColIndex.rho;
      size_t colIndex = elementIndex + localStartColIndex.rho;
      localToGlobalMat.insertElement(rowIndex, colIndex, 1.);
    }
  }

  // phi to phi
  size_t numNode = condMesh.totalNumNode();
  for (size_t nodeIndex = 0; nodeIndex < numNode; nodeIndex++) {
    size_t rowIndex = nodeIndex + globalStartColIndex.phi;
    size_t colIndex = nodeIndex + localStartColIndex.phi;
    localToGlobalMat.insertElement(rowIndex, colIndex, 1.);
  }
}

/**********************************************************************
 * setupGlobalToLocalMat --
 **********************************************************************/
void
Formulation::setupGlobalToLocalMat (
				    void)
{
  globalToLocalMat = SpRowMat<double>(totalNumUnknown_, totalNumUnknown_);

  size_t numElement = allElementList.size();
  for (size_t elementIndex = 0; elementIndex < numElement; elementIndex++) {
    for (size_t k = 0; k < 3; k++) {
      // k = 0: t1
      // k = 1: t2
      // k = 3: n
      size_t rowIndex = elementIndex + k*numElement +localStartRowIndex.equa_1;
      for (size_t m = 0; m < 3; m++) {
	size_t colIndex = elementIndex + m*numElement + globalStartRowIndex.equa_1;
	double projection;
	if (k == 0) {
	  // equa1 x, y, z to t1
	  projection = allElementList[elementIndex].tangent1(m);
	} else if (k == 1) {
	  // equa1 x, y, z to t2
	  projection = allElementList[elementIndex].tangent2(m);
	} else {
	  // equa1 x, y, z to n
	  projection = allElementList[elementIndex].normal(m);
	}
	//globalToLocalMat.insertElement(rowIndex, colIndex, abs(projection));
	globalToLocalMat.insertElement(rowIndex, colIndex, projection);
      }
    }

    // equa2 in both coordinate system is the same
    for (size_t m = 0; m < 2; m++) {
      // m=0: t1
      // m=1: t2
      size_t rowIndex = elementIndex + m*numElement +localStartRowIndex.equa_2;
      size_t colIndex = elementIndex + m*numElement + globalStartRowIndex.equa_2;
      globalToLocalMat.insertElement(rowIndex, colIndex, 1.);
    }

    // equa3
    if (simuType_ != MQS) {
      size_t rowIndex = elementIndex + localStartRowIndex.equa_3;
      size_t colIndex = elementIndex + globalStartRowIndex.equa_3;
      globalToLocalMat.insertElement(rowIndex, colIndex, 1.);
    }

    // equa4_5
    size_t rowIndex = elementIndex + localStartRowIndex.equa_4_5;
    size_t colIndex = elementIndex + globalStartRowIndex.equa_4_5;
    globalToLocalMat.insertElement(rowIndex, colIndex, 1.);
  }

  // equa6_7
  size_t numNode = condMesh.totalNumNode();
  for (size_t nodeIndex = 0; nodeIndex < numNode; nodeIndex++) {
    size_t rowIndex = nodeIndex + localStartRowIndex.equa_6_7;
    size_t colIndex = nodeIndex + globalStartRowIndex.equa_6_7;
    globalToLocalMat.insertElement(rowIndex, colIndex, 1.);
  }
}

/**********************************************************************
 * getRestart --
 * Assume the total memory is 1Gb, use at most 500Mb for Hessenberg matrix
 * in Gmres. From this I could figure out the optimal restart.
 * The key idea is:
 * For small problems, use a fixed restart
 * for large problems, use all the available memory to maximize restart
 **********************************************************************/
size_t
Formulation::getRestart (
			 void)
{
  //  const double memoryQuota = 500e6; // 500Mb
  //  size_t oneIterMemoryUsage = totalNumUnknown_ * sizeof(complex<double>);
  //  double upperBound = floor(memoryQuota / oneIterMemoryUsage);

  double restartTableLookUp;
  if (totalNumUnknown_ <= 100) {
    restartTableLookUp = 10;
  } else if ((totalNumUnknown_ > 100) && (totalNumUnknown_ <= 1e3) ) {
    restartTableLookUp = 20;
  } else if ((totalNumUnknown_ > 1e3) && (totalNumUnknown_ <= 1e4) ) {
    restartTableLookUp = 30;
  } else if ((totalNumUnknown_ > 1e4) && (totalNumUnknown_ <= 1e5) ) {
    restartTableLookUp = 40;
  } else if ((totalNumUnknown_ > 1e5) && (totalNumUnknown_ <= 1e6) ) {
    restartTableLookUp = 50;
  } else {
    restartTableLookUp = 5e-5 * totalNumUnknown_;
  }

  //  double restart = min(upperBound, restartTableLookUp);
  //  return static_cast<size_t>(restart);
  return static_cast<size_t>(restartTableLookUp);
}

/**********************************************************************
 * LUFactor --
 **********************************************************************/
void
Formulation::LUFactor (
		       void)
{
  if (!formOneSystemMatDone) {
    formOneSystemMat();
  }
  superLU.setup(systemMat);
  TNT::stopwatch time;
  time.start();
  superLU.factor();
  timeReport("Time used for LU factorizing the system matrix := ", time.read());
  LUFactorDone = true;
}

/**********************************************************************
 * checkCondNum --
 **********************************************************************/
void
Formulation::checkCondNum (
			   void)
{
  if (!LUFactorDone) {
    LUFactor();
  }
  messageWithOneNumber("Condition number -- one norm := ", 
		       superLU.estimateCondNum_OneNorm());
  messageWithOneNumber("Condition number -- Infinite norm := ", 
		       superLU.estimateCondNum_InfNorm());
}

/**********************************************************************
 * formOneSystemMat --
 * assemble equa1, 2 and 3 and the nonIntegralEqua into one matrix.
 * Also transfer format from compressed row to compressed col.
 **********************************************************************/
void
Formulation::formOneSystemMat (
			       void)
{
  systemMat = SpColMat<complex<double> >(totalNumRow_, totalNumUnknown_);

  // This is to filter some of the zeros in the original system so that superLU 
  // only need to handle a sparser but equivalent matrix.
  const double threshold = 1e-16; 

  StartRowIndex startRowIndex;
  if (! useGlobalCoord_) {
    startRowIndex = localStartRowIndex;
  } else {
    startRowIndex = globalStartRowIndex;
  }

  copyToOneSystemMat(startRowIndex.equa_1, equa1Mat, threshold);
  copyToOneSystemMat(startRowIndex.nonIntegralEqua, nonIntegralEquaMat, threshold);
  if (simuType_ != MQS) {
    copyToOneSystemMat(startRowIndex.equa_3, equa3Mat, threshold);

    for (size_t condIndex = 0; condIndex < condMesh.numCond(); condIndex++) {
      complex<double> coe(0., -1.);
      coe *= 2 * PI * frequency_ * surfConst_.EPSILON0 / 
	condInfoPtrList[condIndex]->conductivity();
      for (size_t panelIndex = 0; panelIndex < condMesh.numPanel(condIndex); panelIndex++) {
	if (condMesh.isNonContact(condIndex, panelIndex)) {
	  size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
	  size_t index = startRowIndex.equa_4_5 + globalPanelIndex;
	  systemMat.insertElement(index, index, coe);
	}
      }
    }
  }

  for (size_t i = 0; i < allNonContactElementList.size(); i++) {
    // the row indices of 2a and 2b are mixed so I have to 
    // copy the rows accordingly
    size_t condIndex = allNonContactElementList[i].boundaryIndex();
    size_t panelIndex = allNonContactElementList[i].index();
    size_t globalElementIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
    // t1 dot (Int (G0 dE/dn) - Int (E dG0/dn))
    size_t rowIndex1 = startRowIndex.equa_2 + globalElementIndex;
    for (size_t colIndex = 0; colIndex < equa2aMat1[i].size(); colIndex++) {
      complex<double> value = equa2aMat1[i].value(colIndex);
      if (abs(value) > threshold) {
	systemMat.insertElement(rowIndex1, equa2aMat1[i].index(colIndex), value);
      }
    }

    // t2 dot (Int (G0 dE/dn) - Int (E dG0/dn))
    size_t rowIndex2 = rowIndex1 + totalNumE_;
    for (size_t colIndex = 0; colIndex < equa2aMat2[i].size(); colIndex++) {
      complex<double> value = equa2aMat2[i].value(colIndex);
      if (abs(value) > threshold) {
	systemMat.insertElement(rowIndex2, equa2aMat2[i].index(colIndex), value);
      }
    }
  }

  formOneSystemMatDone = true;
}

/**********************************************************************
 * copyToOneSystemMat --
 **********************************************************************/
//see formulation.h
