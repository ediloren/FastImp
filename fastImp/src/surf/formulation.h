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

  const static char cvsid[] = "$Id: formulation.h,v 1.11 2003/07/11 02:48:43 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _FORMULATION_H_
#define _FORMULATION_H_

#include <complex>
#include "surfConst.h"
#include "condInfo.h"
#include "mesh.h"
#include "element.h"
#include "calcpForEikrOverR.h"
#include "oneOverR.h"
#include "eikrOverR.h"
#include "staticCollocation.h"
#include "fullwaveCollocation.h"
#include "pfft.h"
#include "spRowMat.h"
#include "dualPanel.h"
#include "superLU.h"

namespace surf {

  class Formulation {

    friend class Preconditioner;

  public:
    //    Formulation(void) {}
    Formulation(const mesh::Mesh&,  const std::vector<CondInfo*>&, 
		const SimulationType, const FormulationType, 
		const SolverType solverType, const SurfConst&, const bool useGlobalCoord);

    void setupSystem(const double freq);

    template <class VecRHS> void setupRHS (const std::vector<double>& leftContactVolt, 
					   const std::vector<double>& rightContactVolt, 
					   VecRHS& RHS);

    template <class VecRHS> void setupRHS (const size_t condIndex, 
					   const double leftContactVolt, 
					   const double rightContactVolt, 
					   VecRHS& RHS);

    template <class VecA, class VecB> 
    void matMultVec (VecA& ans, VecB& x);
    void LUFactor(void);
    void checkCondNum(void);

    template <class VecRHS, class VecX> 
    void directSolve (const VecRHS& RHS, VecX& x) {  
      if (!LUFactorDone) {
	LUFactor();
      }
      superLU.solve(RHS, x);
    }

    template <class VecRHS, class VecX> 
    void checkSolution (const VecRHS& RHS, VecX& x) {  
      // Enrico, added std:: to cout and endl
      std::cout << "\tThe 2-norm of x := " << two_norm(x)
           << "\tThe 2-norm of RHS := " << two_norm(RHS) << std::endl;
      std::cout << "\tThe 2-norm of b - A * x:= "
           << two_norm(RHS - systemMat * x) << std::endl;
      std::cout << "\tThe infinite-norm of b - A * x:= "
           << infinite_norm(RHS - systemMat * x) << std::endl;
    }

    template <class VecA, class VecB> 
    void operator () (VecA& ans, const VecB& x) { matMultVec(ans, x); }
    
    const size_t numCond(void) const { return numCond_;}
    const size_t totalNumE(void) const { return totalNumE_; }
    const size_t totalNumPhi(void) const { return totalNumPhi_; }
    const size_t totalNumRho(void) const { return totalNumRho_; }
    const size_t totalNumUnknown(void) const { return totalNumUnknown_; }
    const size_t totalNumRow(void) const { return totalNumRow_; }
    const double panelArea(const size_t condIndex, 
			   const size_t panelIndex) const {
      size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, 
							  panelIndex);
      return allElementList[globalPanelIndex].area();
    }
    const double maxElementSize(void) const { return maxElementSize_; }
    const double minElementSize(void) const { return minElementSize_; }
    const double minElementSize(const size_t i) const { 
      return minElementSizeList_[i]; }
    const pfft::point3D panelCentroid(const size_t condIndex, 
				      const size_t panelIndex) const {
      size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
      return allElementList[globalPanelIndex].centroid();
    }
    const SimulationType simulationType(void) const { return simuType_; }
    const double frequency(void) const { return frequency_; }

    template<class Vec> std::complex<double> 
    getPhi(const size_t, const size_t, const Vec&) const;

    template<class Vec> std::complex<double> 
    getRho(const size_t, const size_t, const Vec&) const;

    template<class Vec> std::complex<double> 
    getEn(const size_t, const size_t, const Vec&) const;
    template<class Vec> std::complex<double> 
    getEnInLocalCoord(const size_t, const size_t, const Vec&) const;
    template<class Vec> std::complex<double> 
    getEnInGlobalCoord(const size_t, const size_t, const Vec&) const;

    template<class Vec> std::complex<double> 
    getEt1(const size_t, const size_t, const Vec&) const;
    template<class Vec> std::complex<double> 
    getEt1InLocalCoord(const size_t, const size_t, const Vec&) const;
    template<class Vec> std::complex<double> 
    getEt1InGlobalCoord(const size_t, const size_t, const Vec&) const;

    template<class Vec> std::complex<double> 
    getEt2(const size_t, const size_t, const Vec&) const;
    template<class Vec> std::complex<double> 
    getEt2InLocalCoord(const size_t, const size_t, const Vec&) const;
    template<class Vec> std::complex<double> 
    getEt2InGlobalCoord(const size_t, const size_t, const Vec&) const;

    template<class Vec> pfft::vector3D<std::complex<double> > 
    getE(const size_t, const size_t, const Vec&) const;
    template<class Vec> pfft::vector3D<std::complex<double> > 
    getEInLocalCoord(const size_t, const size_t, const Vec&) const;
    template<class Vec> pfft::vector3D<std::complex<double> > 
    getEInGlobalCoord(const size_t, const size_t, const Vec&) const;

    template<class Vec> pfft::vector3D<std::complex<double> > 
    getdEdn(const size_t, const size_t, const Vec&) const;
    template<class Vec> pfft::vector3D<std::complex<double> > 
    getdEdnInLocalCoord(const size_t, const size_t, const Vec&) const;
    template<class Vec> pfft::vector3D<std::complex<double> > 
    getdEdnInGlobalCoord(const size_t, const size_t, const Vec&) const;

    template<class Vec> std::complex<double> 
    getPhi_deNormalized(const size_t, const size_t, const Vec&) const;

    template<class Vec> std::complex<double> 
    getRho_deNormalized(const size_t, const size_t, const Vec&) const;

    template<class Vec> std::complex<double> 
    getEn_deNormalized(const size_t, const size_t, const Vec&) const;

    template<class Vec> pfft::vector3D<std::complex<double> > 
    getE_deNormalized(const size_t, const size_t, const Vec&) const;

    template<class Vec> pfft::vector3D<std::complex<double> > 
    getdEdn_deNormalized(const size_t, const size_t, const Vec&) const;

    void compWaveNumber(const double f, 
			double& K0,
			std::vector<std::complex<double> >& KC) const;

    size_t getRestart(void);

  private:
    // reference must be used here otherwise meshCond will be copied and
    // later deallocated once in formulation and once in meshCond.
    // Since pointers are involved here, this will cause segmentation fault.
    // The side effect is that formulation can not have default constructor.
    // This is fine since formulation is not used as a private member
    // data in surfMain.cc
    const mesh::Mesh& condMesh; 
    const std::vector<CondInfo*> condInfoPtrList;
    SimulationType simuType_;
    FormulationType formulationType_;
    SolverType solverType_;
    bool useGlobalCoord_;
    std::vector<std::vector<pfft::element> > condElementList;
    std::vector<pfft::element> allElementList;
    std::vector<pfft::element> allNonContactElementList;
    std::vector<DualPanel> dualPanelList;
    std::vector<std::complex<double> > w1;
    std::vector<std::complex<double> > w2;
    std::vector<std::complex<double> > xw;
    size_t numCond_;
    size_t totalNumE_;
    size_t totalNumPhi_;
    size_t totalNumRho_;
    size_t totalNumUnknown_;
    size_t totalNumRow_;
    size_t numRowOfNonIntegralMat_;
    double frequency_;
    double maxElementSize_;
    std::vector<double> minElementSizeList_;
    double minElementSize_;
    SurfConst surfConst_;

    // for equa_1
    typedef pfft::EikrOverR<std::complex<double> > Kernel1;
    typedef pfft::FullwaveCollocation<std::complex<double>, pfft::element> Calcp1;
    typedef pfft::Pfft<std::complex<double>, std::complex<double>, 
      Kernel1, Calcp1> Pfft1;
    std::vector<std::complex<double> > KC;
    std::vector<pfft::DifferentialOperator> innerOperatorList1;
    std::vector<pfft::DifferentialOperator> outerOperatorList1;
    std::vector<Kernel1> kernel1List;
    std::vector<Calcp1> calcp1List;
    std::vector<Pfft1> pfft1List;

    // for equa_2a
    std::vector<pfft::DifferentialOperator> innerOperatorList2;
    std::vector<pfft::DifferentialOperator> outerOperatorList2;

    typedef pfft::OneOverR Kernel2_EMQS;
    typedef pfft::StaticCollocation<pfft::element> Calcp2_EMQS;
    typedef pfft::Pfft<double, std::complex<double>, 
      Kernel2_EMQS, Calcp2_EMQS> Pfft2_EMQS;
    Kernel2_EMQS kernel2_EMQS;
    Calcp2_EMQS calcp2_EMQS;
    Pfft2_EMQS pfft2_EMQS;

    typedef pfft::EikrOverR<double> Kernel2_FW;
    typedef pfft::FullwaveCollocation<double, pfft::element> Calcp2_FW;
    typedef pfft::Pfft<std::complex<double>, std::complex<double>, 
      Kernel2_FW, Calcp2_FW> Pfft2_FW;
    double K0;
    Kernel2_FW kernel2_FW;
    Calcp2_FW calcp2_FW;
    Pfft2_FW pfft2_FW;

    // for equa_3 
    std::vector<pfft::DifferentialOperator> innerOperatorList3;
    std::vector<pfft::DifferentialOperator> outerOperatorList3;

    typedef Kernel2_FW Kernel3_FW;
    typedef Calcp2_FW Calcp3_FW;
    typedef Pfft2_FW Pfft3_FW;
    Kernel3_FW kernel3_FW;
    Calcp3_FW calcp3_FW;
    Pfft3_FW pfft3_FW;

    typedef Kernel2_EMQS Kernel3_EMQS;
    typedef Calcp2_EMQS Calcp3_EMQS;
    typedef Pfft2_EMQS Pfft3_EMQS;
    Kernel3_EMQS kernel3_EMQS;
    Calcp3_EMQS calcp3_EMQS;
    Pfft3_EMQS pfft3_EMQS;

    // for non-integral part of the system matrix
    pfft::SpRowMat<double> nonIntegralEquaMat;

    // conversion between global and local coordinate system
    pfft::SpRowMat<double> localToGlobalMat;
    pfft::SpRowMat<double> globalToLocalMat;

    // when pFFT is not used, the part of system matrix corresponding to
    // integral equations are treated as sparse matrices. 
    pfft::SpRowMat<std::complex<double> > equa1Mat;
    pfft::SpRowMat<std::complex<double> > equa2aMat1;
    pfft::SpRowMat<std::complex<double> > equa2aMat2;
    pfft::SpRowMat<std::complex<double> > equa3Mat;

    // for estimate condition number and even directly solve the system
    pfft::SpColMat<std::complex<double> > systemMat;
    SuperLU superLU;
    bool LUFactorDone;
    bool formOneSystemMatDone;

    /* this struct contains the starting column index for each variable 
       index 0 means x component
       index 1 means y component
       index 2 means z component
       The number of unknowns for E, dEdn and rho is the number of panels
       The number of unknowns for phi is number of vertices
    */
    typedef struct StartColIndex {
      size_t E[3];
      size_t dEdn[3];
      size_t phi;
      size_t rho;
    } StartColIndex;
    StartColIndex globalStartColIndex; // Ex, Ey, Ez, dEx, dEy, dEz
    StartColIndex localStartColIndex; // Et1, Et2, En, dEt1, dEt2, dEn

    typedef struct StartRowIndex {
      size_t equa_1;
      size_t equa_2;
      size_t equa_3;
      size_t equa_4_5;
      size_t equa_6_7;
      size_t nonIntegralEqua;
    } StartRowIndex;
    StartRowIndex globalStartRowIndex; // equa_1 in global, equa_2 in local coord.
    StartRowIndex localStartRowIndex; // all in local coord.

    bool printMatMultVecTime_;
    TNT::stopwatch timer;

    void setupElementList(void);
    void setupDualPanelList(void);
    double compTriangleArea(const pfft::point3D&, const pfft::point3D&, 
			    const pfft::point3D&) const;
    void setupGlobalStartRowIndex(void);
    void setupGlobalStartColIndex(void);
    void allocateWorkUnit(void);

    template <class VecA, class VecB> void equa1 (VecA& ans, const VecB& x);
    template <class VecA, class VecB> void equa2a (VecA& ans, const VecB& x);
    template <class VecA, class VecB> void equa2a_direct (VecA& ans, const VecB& x);
    template <class VecA, class VecB> void equa3 (VecA& ans, const VecB& x);
    template <class VecA, class VecB> void nonIntegralEqua (VecA& ans, const VecB& x);

    void setupIntegralOperatorList (void);
    void setupIntegralEqua_frequencyDependentPart(void);
    void setupIntegralEqua_frequencyIndependentPart(void);
    void setupEqua1_frequencyDependentPart(void);
    void setupEqua1_frequencyIndependentPart(void);
    void setupEqua2a_frequencyDependentPart(void);
    void setupEqua2a_frequencyIndependentPart(void);
    void setupEqua3_frequencyDependentPart(void);
    void setupEqua3_frequencyIndependentPart(void);
    void fillIntegralEqua (void);
    void fillEqua_1_old (const size_t, const size_t, const size_t);
    void fillEqua_1_new (const size_t, const size_t, const size_t);
    std::complex<double> fillEqua_2a (const size_t, const size_t, const size_t, const size_t, const size_t);
    void fillEqua_3 (const size_t, const size_t, const size_t, const size_t, 
		     std::complex<double>&);
#ifdef DEBUG_RHO
    void fillEqua_3_contact (const size_t globalPanelIndex, const size_t rowIndex);
#endif

    void setupNonIntegralEqua(void);
    void fillEqua_1_old(void);
    void fillEqua_1_new(void);
    void fillEqua_2a_sparsePart(const size_t, const size_t, 
				const size_t, const size_t);
    void compGradPhiForTriangle(const size_t, const size_t, 
				const size_t, const size_t);
    void compGradPhiForQuad(const size_t, const size_t, 
			    const size_t, const size_t);
    void fillEqua_2b(const size_t, const size_t);
    void fillEqua_3_sparsePart(const size_t, const size_t, const size_t);
    void fillEqua_4(const size_t, const size_t);
    void fillEqua_5(const size_t, const size_t);
    void fillEqua_6(const size_t, const size_t);
    void fillEqua_7(const size_t, const size_t);

    void setupLocalStartRowIndex(void);
    void setupLocalStartColIndex(void);
    void setupLocalToGlobalMat(void);
    void setupGlobalToLocalMat(void);

    void formOneSystemMat(void);
    template<class T> void
    copyToOneSystemMat(const size_t startRow, const pfft::SpRowMat<T>& mat, 
		       const double threshold);

    void checkMemoryUsage (void);
    size_t estimateMemoryUsage_pfft (void);
    size_t estimateMemoryUsage_systemMat (void) const;
    void memoryReport (const std::string object, 
		       const size_t memoryUsage) const;
  };

  /**********************************************************************
   * setupRHS --
   **********************************************************************/
  template <class VecRHS>
  void 
  Formulation::setupRHS (
			 const std::vector<double>& leftContactVolt, 
			 const std::vector<double>& rightContactVolt, 
			 VecRHS& RHS)
  {
    for (size_t i = 0; i < totalNumUnknown_; i++) {
      RHS[i] = 0.;
    }
    
    for (size_t condIndex = 0; condIndex < condMesh.numCond(); condIndex++) {
      for (size_t i = 0; i < condMesh.numNode(condIndex); i++) {
	size_t gi = condMesh.globalNodeIndex(condIndex, i);
	size_t row;
	if (! useGlobalCoord_) {
	  row = localStartRowIndex.equa_6_7 + gi;
	} else {
	  row = globalStartRowIndex.equa_6_7 + gi;
	}
        // Enrico
        //if (condMesh.vertexType(gi) == LEFT_CONTACT) {
        if (condMesh.vertexType(gi) == mesh::LEFT_CONTACT) {
	  RHS[row] = leftContactVolt[condIndex];
        // Enrico
        //} else if (condMesh.vertexType(gi) == RIGHT_CONTACT) {
        } else if (condMesh.vertexType(gi) == mesh::RIGHT_CONTACT) {
	  RHS[row] = rightContactVolt[condIndex];
	}
      }
    }
  }

  /**********************************************************************
   * setupRHS --
   **********************************************************************/
  template <class VecRHS>
  void 
  Formulation::setupRHS (
			 const size_t condIndex, 
			 const double leftContactVolt, 
			 const double rightContactVolt, 
			 VecRHS& RHS)
  {
    for (size_t i = 0; i < totalNumUnknown_; i++) {
      RHS[i] = 0.;
    }
    
    for (size_t i = 0; i < condMesh.numNode(condIndex); i++) {
      size_t gi = condMesh.globalNodeIndex(condIndex, i);
      size_t row;
      if (! useGlobalCoord_) {
	row = localStartRowIndex.equa_6_7 + gi;
      } else {
	row = globalStartRowIndex.equa_6_7 + gi;
      }
      // Enrico
      //if (condMesh.vertexType(gi) == LEFT_CONTACT) {
      if (condMesh.vertexType(gi) == mesh::LEFT_CONTACT) {
	RHS[row] = leftContactVolt;
      // Enrico
      //} else if (condMesh.vertexType(gi) == RIGHT_CONTACT) {
      } else if (condMesh.vertexType(gi) == mesh::RIGHT_CONTACT) {
	RHS[row] = rightContactVolt;
      }
    }
  }
  
  /**********************************************************************
   * matMultVec --
   * The index of unknown vector x is consistant with panel index.
   * the index of vector ans, or the index of equations, is not.
   **********************************************************************/
  template <class VecA, class VecB> 
  void
  Formulation::matMultVec (
			   VecA& ans, 
			   VecB& x)
  {
    if (printMatMultVecTime_)  timer.start();

    if (! useGlobalCoord_) {
      x = localToGlobalMat * x;
    }

    if (solverType_ == ITERATIVE_WITH_PFFT) {
      nonIntegralEqua(ans, x);
      equa1(ans, x);
      equa2a(ans, x);
      if (simuType_ != MQS) {
	equa3(ans, x);
      }
    } else {
      pfft::matMultVec(ans, equa1Mat, x);
      nonIntegralEqua(ans, x);
      equa2a_direct(ans, x);
      if (simuType_ != MQS) {
	equa3(ans, x);
      }
    }

    if (! useGlobalCoord_) {
      ans = globalToLocalMat * ans;
    }

    if (printMatMultVecTime_)  {
      timeReport("Time used for one matrix vector product := ", timer.read());
      printMatMultVecTime_ = false;
    }
  }

  /**********************************************************************
   * equa1 --
   * Old:   Int (G1 dE/dn') - Int (E dG1/dn') = 4 * PI * E 
   * New:   d/dn { Int (G1 dE/dn') - Int (E dG1/dn') = 4 * PI * E }
   * Enforced on each conductor seperately
   **********************************************************************/
  template <class VecA, class VecB> 
  void
  Formulation::equa1 (
		      VecA& ans, 
		      const VecB& x)
  {
    for (size_t m = 0; m < 3; m++) {
      /* m = 0:  x equation
	 m = 1:  y equation
	 m = 2:  z equation
      */

      size_t rowIndexStart = globalStartRowIndex.equa_1 + m*totalNumE_;
      size_t colIndexStart_dEdn = globalStartColIndex.dEdn[m];
      size_t colIndexStart_E = globalStartColIndex.E[m];
      for (size_t condIndex = 0; condIndex < numCond_; condIndex++) {
	size_t numRow = condElementList[condIndex].size();
	size_t numCol = numRow;
      
	// w1 = Int (G1 dE/dn')
	pfft1List[condIndex].IEoperator(w1, outerOperatorList1[0], 
					innerOperatorList1[0], 
					&x[colIndexStart_dEdn]);

	// w2 = Int (E dG1/dn')
	pfft1List[condIndex].IEoperator(w2, outerOperatorList1[1], 
					innerOperatorList1[1], 
					&x[colIndexStart_E]);

	if (formulationType_ == OLD) {
	  for (size_t i = 0; i < numRow; i++) {
	    ans[rowIndexStart + i] = w1[i] - w2[i] - 4*PI*x[colIndexStart_E + i];
	  }
	} else if (formulationType_ == NEW) {
	  for (size_t i = 0; i < numRow; i++) {
	    ans[rowIndexStart + i] = w1[i] - w2[i] - 4*PI*x[colIndexStart_dEdn + i];
	  }
	}

	colIndexStart_E += numCol;
	colIndexStart_dEdn += numCol;
	rowIndexStart += numRow;
      }
    }
  }

  /**********************************************************************
   * equa2a --
   * t dot { Int (G0 dE/dn) - Int(E dG0/dn) } + 4 * PI * t dot grad(phi) = 0
   * Enforced only along two tangent directions on the union of non-contact 
   * panels 
   **********************************************************************/
  template <class VecA, class VecB> 
  void
  Formulation::equa2a (
		       VecA& ans, 
		       const VecB& x)
  {
    for (size_t m = 0; m < 3; m++) {
      /* m = 0: Ex and dEx/dn column 
	 m = 1: Ey and dEy/dn column
	 m = 2: Ez and dEz/dn column
      */

      // Int (G0 dE/dn)
      if (simuType_ == FULL_WAVE) {
	pfft2_FW.IEoperator(
			    w1,
			    outerOperatorList2[0],
			    innerOperatorList2[0], 
			    &x[globalStartColIndex.dEdn[m]]);
      // Int (E dG0/dn)
	pfft2_FW.IEoperator(
			    w2,
			    outerOperatorList2[1],
			    innerOperatorList2[1], 
			    &x[globalStartColIndex.E[m]]);
      } else {
	pfft2_EMQS.IEoperator(
			      w1,
			      outerOperatorList2[0],
			      innerOperatorList2[0],
			      &x[globalStartColIndex.dEdn[m]]);

	pfft2_EMQS.IEoperator(
			      w2,
			      outerOperatorList2[1],
			      innerOperatorList2[1], 
			      &x[globalStartColIndex.E[m]]);
      }

      for (size_t i = 0; i < allNonContactElementList.size(); i++) {
	// the row indices of 2a and 2b are mixed so I have to 
	// fill in the result vector accordingly
	size_t condIndex = allNonContactElementList[i].boundaryIndex();
	size_t panelIndex = allNonContactElementList[i].index();
	size_t globalElementIndex = 
	  condMesh.globalPanelIndex(condIndex, panelIndex);
	// t1 dot (Int (G0 dE/dn) - Int (E dG0/dn))
	ans[globalStartRowIndex.equa_2 + globalElementIndex] += 
	  (w1[i] - w2[i]) * allNonContactElementList[i].tangent1(m);
	// t2 dot (Int (G0 dE/dn) - Int (E dG0/dn))
	ans[globalStartRowIndex.equa_2 + totalNumE_ + globalElementIndex] += 
	  (w1[i] - w2[i]) * allNonContactElementList[i].tangent2(m);
      }
    }
  }

  /**********************************************************************
   * equa2a_direct --
   * t dot { Int (G0 dE/dn) - Int(E dG0/dn) } + 4 * PI * t dot grad(phi) = 0
   * Enforced only along two tangent directions on the union of non-contact 
   * panels 
   **********************************************************************/
  template <class VecA, class VecB> 
  void
  Formulation::equa2a_direct (
			      VecA& ans, 
			      const VecB& x)
  {
    // t1 * { Int (G0 dE/dn) - Int (E dG0/dn) }
    pfft::matMultVec(w1, equa2aMat1, x);
    // t2 * { Int (G0 dE/dn) - Int (E dG0/dn) }
    pfft::matMultVec(w2, equa2aMat2, x);

    for (size_t i = 0; i < allNonContactElementList.size(); i++) {
      // the row indices of 2a and 2b are mixed so I have to 
      // fill in the result vector accordingly
      size_t condIndex = allNonContactElementList[i].boundaryIndex();
      size_t panelIndex = allNonContactElementList[i].index();
      size_t globalElementIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
      // t1 dot (Int (G0 dE/dn) - Int (E dG0/dn))
      ans[globalStartRowIndex.equa_2 + globalElementIndex] += w1[i];
      // t2 dot (Int (G0 dE/dn) - Int (E dG0/dn))
      ans[globalStartRowIndex.equa_2 + totalNumE_ + globalElementIndex] += w2[i];
    }
  }

  /**********************************************************************
   * equa3 --
   *       Int (G0 rho / epsilon0) = 4 * PAI * phi
   * The equation is enforced on the union of all panels
   * This equation is not used for MQS analysis
   *
   * Note: 
   * 1) To avoid bad scaling, let rho = rho / epsilon0, the quation now is
   * Int (G0 rho) = 4 * PAI * phi
   * 2) phi is defined on vertex. So its value at collocation point, 
   * the centroid, is the average of values on vertices of a panel
   **********************************************************************/
  template <class VecA, class VecB> 
  void
  Formulation::equa3 (
		      VecA& ans, 
		      const VecB& x)
  {
    if (solverType_ == ITERATIVE_WITH_PFFT) {
      // Int G0 rho
      if (simuType_ == FULL_WAVE) {
	pfft3_FW.IEoperator(
			    w1,
			    outerOperatorList3[0],
			    innerOperatorList3[0], 
			    &x[globalStartColIndex.rho]);
      } else {
	pfft3_EMQS.IEoperator(
			      w1,
			      outerOperatorList3[0],
			      innerOperatorList3[0], 
			      &x[globalStartColIndex.rho]);
      }
    } else {
      // Int G0 rho
      pfft::matMultVec(w1, equa3Mat, x);
    }

#ifdef DEBUG_RHO
    for (size_t i = 0; i < allNonContactElementList.size(); i++) {
      // the row indices of 3a and 3b are mixed so I have to 
      // fill in the result vector accordingly
      size_t condIndex = allNonContactElementList[i].boundaryIndex();
      size_t panelIndex = allNonContactElementList[i].index();
      size_t globalElementIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
      ans[globalStartRowIndex.equa_3 + globalElementIndex] += w1[i];
    }
#else
    for (size_t i = 0; i < allElementList.size(); i++) {
      ans[globalStartRowIndex.equa_3 + i] += w1[i];
    }
#endif
  }

  /**********************************************************************
   * nonIntegralEqua --
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
   * Note:
   * 1) The only frequency-dependent term in these equations is the one
   * on the right-hand side of equa4. This is specially handled here
   * because the sparse matrix that represents the rest of other equations 
   * are frequency-independent and is setup only once.
   * 2) To make sure this special handled block is a diag matrix, the rowIndex
   * of equa4 and equa5 has to be mixed. This should be consistant
   * with the function setupNonIntegralEqua()
   **********************************************************************/
  template <class VecA, class VecB> 
  void
  Formulation::nonIntegralEqua (
				VecA& ans, 
				const VecB& x)
  {
    pfft::matMultVecPartial(&ans[globalStartRowIndex.nonIntegralEqua], 
			    nonIntegralEquaMat, x);

#ifdef DEBUG_EMQS
#else
    if (simuType_ != MQS) {
      for (size_t condIndex = 0; condIndex < condMesh.numCond(); condIndex++) {
	std::complex<double> coe(0., -1.);
	coe *= 2. * PI * frequency_ * surfConst_.EPSILON0 / 
	  condInfoPtrList[condIndex]->conductivity();
	for (size_t panelIndex = 0; panelIndex < condMesh.numPanel(condIndex); 
	     panelIndex++) {
	  size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex,
							      panelIndex);
	  if (condMesh.isNonContact(condIndex, panelIndex)) {
	    ans[globalStartRowIndex.equa_4_5 + globalPanelIndex] += 
	      coe * x[globalStartColIndex.rho + globalPanelIndex];
	  }
	}
      }
    }
#endif
  }

  /**********************************************************************
   * getE --
   **********************************************************************/
  template<class Vec>
  pfft::vector3D<std::complex<double> >
  Formulation::getE (
		     const size_t condIndex,
		     const size_t panelIndex,
		     const Vec& solution) const 
  {
    if (! useGlobalCoord_) 
      return getEInLocalCoord(condIndex, panelIndex, solution);
    else
      return getEInGlobalCoord(condIndex, panelIndex, solution);
  }

  /**********************************************************************
   * getEInGlobalCoord --
   **********************************************************************/
  template<class Vec>
  pfft::vector3D<std::complex<double> >
  Formulation::getEInGlobalCoord (
				  const size_t condIndex,
				  const size_t panelIndex,
				  const Vec& solution) const 
  {
    size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
    size_t j = globalPanelIndex + globalStartColIndex.E[0];
    std::complex<double> Ex = solution[j];
    j = globalPanelIndex + globalStartColIndex.E[1];
    std::complex<double> Ey = solution[j];
    j = globalPanelIndex + globalStartColIndex.E[2];
    std::complex<double> Ez = solution[j];

    return pfft::vector3D<std::complex<double> >(Ex, Ey, Ez);
  }

  /**********************************************************************
   * getEInLocalCoord --
   * E are all defined in local coordinate system (t1, t2, n). 
   * This function transfers E to global coordiante system.
   **********************************************************************/
  template<class Vec>
  pfft::vector3D<std::complex<double> >
  Formulation::getEInLocalCoord (
				 const size_t condIndex,
				 const size_t panelIndex,
				 const Vec& solution) const 
  {
    size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
    size_t j = globalPanelIndex + localStartColIndex.E[0];
    std::complex<double> Et1 = solution[j];
    j = globalPanelIndex + localStartColIndex.E[1];
    std::complex<double> Et2 = solution[j];
    j = globalPanelIndex + localStartColIndex.E[2];
    std::complex<double> En = solution[j];
    pfft::vector3D<std::complex<double> > E(Et1, Et2, En);

    E.rotateCoord(allElementList[globalPanelIndex].tangent1(), 
		  allElementList[globalPanelIndex].tangent2(),
		  allElementList[globalPanelIndex].normal());
    return E;
  }

  /**********************************************************************
   * getEt1 --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getEt1 (
		       const size_t condIndex,
		       const size_t panelIndex,
		       const Vec& solution) const 
  {
    if (! useGlobalCoord_) 
      return getEt1InLocalCoord(condIndex, panelIndex, solution);
    else
      return getEt1InGlobalCoord(condIndex, panelIndex, solution);
  }

  /**********************************************************************
   * getEt1InLocalCoord --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getEt1InLocalCoord (
				   const size_t condIndex,
				   const size_t panelIndex,
				   const Vec& solution) const 
  {
    size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
    size_t j = globalPanelIndex + localStartColIndex.E[0];
    return solution[j];
  }

  /**********************************************************************
   * getEt1InGlobalCoord --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getEt1InGlobalCoord (
				    const size_t condIndex,
				    const size_t panelIndex,
				    const Vec& solution) const 
  {
    size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
    pfft::vector3D<double> t1 = allElementList[globalPanelIndex].tangent1();
    pfft::vector3D<std::complex<double> > E = getE(condIndex, panelIndex, solution);
    std::complex<double> Et1;
    dotProd(Et1, E, t1);
    return Et1;
  }

  /**********************************************************************
   * getEt2 --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getEt2 (
		       const size_t condIndex,
		       const size_t panelIndex,
		       const Vec& solution) const 
  {
    if (! useGlobalCoord_) 
      return getEt2InLocalCoord(condIndex, panelIndex, solution);
    else
      return getEt2InGlobalCoord(condIndex, panelIndex, solution);
  }

  /**********************************************************************
   * getEt2InLocalCoord --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getEt2InLocalCoord (
				   const size_t condIndex,
				   const size_t panelIndex,
				   const Vec& solution) const 
  {
    size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
    size_t j = globalPanelIndex + localStartColIndex.E[1];
    return solution[j];
  }

  /**********************************************************************
   * getEt2InGlobalCoord --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getEt2InGlobalCoord (
				    const size_t condIndex,
				    const size_t panelIndex,
				    const Vec& solution) const 
  {
    size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
    pfft::vector3D<double> t2 = allElementList[globalPanelIndex].tangent2();
    pfft::vector3D<std::complex<double> > E = getE(condIndex, panelIndex, solution);
    std::complex<double> Et2;
    dotProd(Et2, E, t2);
    return Et2;
  }

  /**********************************************************************
   * getEn --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getEn (
		      const size_t condIndex,
		      const size_t panelIndex,
		      const Vec& solution) const 
  {
    if (! useGlobalCoord_) 
      return getEnInLocalCoord(condIndex, panelIndex, solution);
    else
      return getEnInGlobalCoord(condIndex, panelIndex, solution);
  }

  /**********************************************************************
   * getEnInLocalCoord --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getEnInLocalCoord (
				  const size_t condIndex,
				  const size_t panelIndex,
				  const Vec& solution) const 
  {
    size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
    size_t j = globalPanelIndex + localStartColIndex.E[2];
    return solution[j];
  }

  /**********************************************************************
   * getEnInGlobalCoord --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getEnInGlobalCoord (
				   const size_t condIndex,
				   const size_t panelIndex,
				   const Vec& solution) const 
  {
    size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
    pfft::vector3D<double> panelNormal = allElementList[globalPanelIndex].normal();
    pfft::vector3D<std::complex<double> > E = getE(condIndex, panelIndex, solution);
    std::complex<double> En;
    dotProd(En, E, panelNormal);
    return En;
  }

  /**********************************************************************
   * getdEdn --
   **********************************************************************/
  template<class Vec>
  pfft::vector3D<std::complex<double> >
  Formulation::getdEdn (
			const size_t condIndex,
			const size_t panelIndex,
			const Vec& solution) const
  {
    if (! useGlobalCoord_) 
      return getdEdnInLocalCoord(condIndex, panelIndex, solution);
    else
      return getdEdnInGlobalCoord(condIndex, panelIndex, solution);
  }

  /**********************************************************************
   * getdEdnInLocalCoord --
   * dEdn are all defined in local coordinate system (t1, t2, n). 
   * This function transfers dEdn to global coordiante system.
   **********************************************************************/
  template<class Vec>
  pfft::vector3D<std::complex<double> >
  Formulation::getdEdnInLocalCoord (
				    const size_t condIndex,
				    const size_t panelIndex,
				    const Vec& solution) const
  {
    size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
    size_t j = globalPanelIndex + localStartColIndex.dEdn[0];
    std::complex<double> dEt1dn = solution[j];
    j = globalPanelIndex + localStartColIndex.dEdn[1];
    std::complex<double> dEt2dn = solution[j];
    j = globalPanelIndex + localStartColIndex.dEdn[2];
    std::complex<double> dEndn = solution[j];
    pfft::vector3D<std::complex<double> > dEdn(dEt1dn, dEt2dn, dEndn);

    dEdn.rotateCoord(allElementList[globalPanelIndex].tangent1(), 
		     allElementList[globalPanelIndex].tangent2(),
		     allElementList[globalPanelIndex].normal());
    return dEdn;
  }

  /**********************************************************************
   * getdEdnInGlobalCoord --
   **********************************************************************/
  template<class Vec>
  pfft::vector3D<std::complex<double> >
  Formulation::getdEdnInGlobalCoord (
				     const size_t condIndex,
				     const size_t panelIndex,
				     const Vec& solution) const
  {
    size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
    size_t j = globalPanelIndex + globalStartColIndex.dEdn[0];
    std::complex<double> dExdn = solution[j];
    j = globalPanelIndex + globalStartColIndex.dEdn[1];
    std::complex<double> dEydn = solution[j];
    j = globalPanelIndex + globalStartColIndex.dEdn[2];
    std::complex<double> dEzdn = solution[j];

    return pfft::vector3D<std::complex<double> >(dExdn, dEydn, dEzdn);
  }

  /**********************************************************************
   * getPhi --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getPhi (
		       const size_t condIndex,
		       const size_t nodeIndexOnOneCond,
		       const Vec& solution) const
  {
    size_t globalNodeIndex = condMesh.globalNodeIndex(condIndex, nodeIndexOnOneCond);
    if (! useGlobalCoord_) 
      return solution[globalNodeIndex + localStartColIndex.phi];
    else
      return solution[globalNodeIndex + globalStartColIndex.phi];
  }

  /**********************************************************************
   * getRho --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getRho (
		       const size_t condIndex,
		       const size_t panelIndex,
		       const Vec& solution) const
  {
    size_t globalPanelIndex = condMesh.globalPanelIndex(condIndex, panelIndex);
    if (! useGlobalCoord_) 
      return solution[globalPanelIndex + localStartColIndex.rho];
    else
      return solution[globalPanelIndex + globalStartColIndex.rho];
  }

  /**********************************************************************
   * getE_deNormalized --
   **********************************************************************/
  template<class Vec>
  pfft::vector3D<std::complex<double> >
  Formulation::getE_deNormalized (
				  const size_t condIndex,
				  const size_t panelIndex,
				  const Vec& solution) const 
  {
    return getE(condIndex, panelIndex, solution);
  }

  /**********************************************************************
   * getEn_deNormalized --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getEn_deNormalized (
				  const size_t condIndex,
				  const size_t panelIndex,
				  const Vec& solution) const 
  {
    return getEn(condIndex, panelIndex, solution);
  }

  /**********************************************************************
   * getdEdn_deNormalized --
   **********************************************************************/
  template<class Vec>
  pfft::vector3D<std::complex<double> >
  Formulation::getdEdn_deNormalized (
				     const size_t condIndex,
				     const size_t panelIndex,
				     const Vec& solution) const
  {
    return getdEdn(condIndex, panelIndex, solution) * surfConst_.normFactor;
  }

  /**********************************************************************
   * getPhi_deNormalized --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getPhi_deNormalized (
				    const size_t condIndex,
				    const size_t nodeIndexOnOneCond,
				    const Vec& solution) const
  {
    return getPhi(condIndex, nodeIndexOnOneCond, solution) / surfConst_.normFactor;
  }

  /**********************************************************************
   * getRho_deNormalized --
   **********************************************************************/
  template<class Vec>
  std::complex<double>
  Formulation::getRho_deNormalized (
				    const size_t condIndex,
				    const size_t panelIndexw,
				    const Vec& solution) const
  {
    return getRho(condIndex, condIndex, solution) * surfConst_.EPSILON0;
  }

  /**********************************************************************
   * copyToOneSystemMat --
   **********************************************************************/
  template<class T>
  void
  Formulation::copyToOneSystemMat (
				   const size_t startRow,
				   const pfft::SpRowMat<T>& mat,
				   const double threshold)
  {
    for (size_t rowIndex = 0; rowIndex < mat.numRow(); rowIndex++) {
      for (size_t colIndex = 0; colIndex < mat[rowIndex].size(); colIndex++) {
	T value = mat[rowIndex].value(colIndex);
    // Enrico, calling double std::abs(double) instead of casting to int abs(int)
//	if (abs(value) > threshold) {
	if (std::abs(value) > threshold) {
	  systemMat.insertElement(rowIndex + startRow, 
				  mat[rowIndex].index(colIndex), value);
	}
      }
    }
  }

} // namespace surf

#endif
