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

  const static char cvsid[] = "$Id: preconditioner.cc,v 1.14 2002/12/31 21:30:39 zhzhu Exp $";

  ==========================================================================
*/

#include "preconditioner.h"

using namespace surf;
using namespace mesh;
using namespace std;
using namespace pfft;

/**********************************************************************
 * Preconditioner --
 **********************************************************************/
Preconditioner::Preconditioner (
				const Formulation& formulationIn, 
				const PreconditionerType precondTypeIn,
				const SurfConst& surfConst,
				const bool useGlobalCoord)
  : fp(formulationIn), precondType(precondTypeIn), totalNumE_(fp.totalNumE_), 
  totalNumRow_(fp.totalNumRow_), printPrecondSolveTime_(true), 
  surfConst_(surfConst), useGlobalCoord_(useGlobalCoord)
{
  if (precondType != NO_PRE_CONDITIONER) {
    innerOperatorList1.resize(1);
    outerOperatorList1.resize(1);
    innerOperatorList2.resize(1);
    outerOperatorList2.resize(1);
    calcp1List.resize(fp.numCond_);
  }
}

/**********************************************************************
 * setup --
 **********************************************************************/
void
Preconditioner::setup (
		       const double freq)
{
  if (precondType == NO_PRE_CONDITIONER) return;

  // Instaed of using K0 and KC in formulation, I compute them here so that 
  // I could have the flexibility of picking different frequency points for
  // preconditioner. The ideal case would be that one preconditioner at a 
  // certain frequency point is good for formulation at all frequency points. 
  // This will certainly save LU factorization time for constructing the 
  // pre-conditioner.
  fp.compWaveNumber(freq, K0, KC);
  
  // Only self-interaction is considered in preconditioner. 
  if (fp.formulationType_ == OLD) {  
    // double-layer is always -2*PI, no need to compute it.
    // Only single-layer is to be computed.
    outerOperatorList1[0] = NONE;
    innerOperatorList1[0] = NONE;
  } else if (fp.formulationType_ == NEW) {  
    // d/dn(single-layer) is always 2*PI, no need to compute it.
    // Only d/dn(double-layer) is to be computed.
    outerOperatorList1[0] = D_DN;
    innerOperatorList1[0] = D_DN;
  }
  for (size_t i = 0; i < fp.numCond_; i++) {
    calcp1List[i] = Formulation::Calcp1(outerOperatorList1, 
					innerOperatorList1, KC[i]);
  }

  outerOperatorList2[0] = NONE;
  innerOperatorList2[0] = NONE;

  if (fp.simuType_ == FULL_WAVE) {
    calcp2_FW = Formulation::Calcp2_FW(outerOperatorList2, innerOperatorList2, K0);
  } else {
    calcp2_EMQS = Formulation::Calcp2_EMQS(outerOperatorList2, innerOperatorList2);
  }

  if (fp.simuType_ == FULL_WAVE) {
    calcp3_FW = Formulation::Calcp3_FW(outerOperatorList2, innerOperatorList2, K0);
  } else if (fp.simuType_ == EMQS){
    calcp3_EMQS = Formulation::Calcp3_EMQS(outerOperatorList2, innerOperatorList2);
  }

#ifdef DEBUG_PRE_COND
  rowBookKeeping.resize(totalNumRow_);
  colBookKeeping.resize(totalNumRow_);
  rowSum.resize(totalNumRow_);
  colSum.resize(totalNumRow_);
#endif

  preCondMat = SpColMat<complex<double> >(fp.totalNumUnknown_, fp.totalNumUnknown_);

  TNT::stopwatch time;  time.start();
  fillSelfTermOfIntegralEqua();
  if (! useGlobalCoord_) {
    fillNonIntegralEquaMatInLocalCoord();
  } else {
    copyNonIntegralEquaMat();
  }

#ifdef DEBUG_EMQS
#else
  if (fp.simuType_ != MQS) {
    fillEqua4_freq_depend(freq);
  }
#endif

  superLU.setup(preCondMat);
  timeReport("Time used for forming the pre-conditioner matrix := ", 
	     time.read());
  
#ifdef DEBUG_PRE_COND
  checkRowBookKeeping();
  checkColBookKeeping();
  checkRowSum();
  checkColSum();
  std::ofstream fout("precond.tmp");
  fout << preCondMat;
#endif

  time.reset();   time.start();
  superLU.factor();
  timeReport("Time used for LU factorizing the pre-conditioner matrix := ", 
	     time.read());

}

/**********************************************************************
 * fillSelfTermOfIntegralEqua --
 **********************************************************************/
void
Preconditioner::fillSelfTermOfIntegralEqua (
					    void)
{
  if (fp.formulationType_ == OLD) {  
    fillEqua1_old();
  } else {
    fillEqua1_new();
  }
  fillEqua2a();
  if (fp.simuType_ != MQS) 
    fillEqua3();
}

/**********************************************************************
 * fillEqua1_old --
 **********************************************************************/
void
Preconditioner::fillEqua1_old (
			       void)
{
  if (! useGlobalCoord_) {
    fillEqua1_old_inLocalCoord();
  } else {
    fillEqua1_old_inGlobalCoord();
  }
}

/**********************************************************************
 * fillEqua1_old_inLocalCoord --
 **********************************************************************/
void
Preconditioner::fillEqua1_old_inLocalCoord (
					    void)
{
  for (size_t condIndex = 0; condIndex < fp.numCond_; condIndex++) {
    for (size_t panelIndex = 0; panelIndex < fp.condMesh.numPanel(condIndex); 
	 panelIndex++) {
      size_t globalPanelIndex = fp.condMesh.globalPanelIndex(condIndex, 
							      panelIndex);
      calcp1List[condIndex](fp.allElementList[globalPanelIndex],
			    fp.allElementList[globalPanelIndex]);
      complex<double> singleLayer = calcp1List[condIndex].result(0);
      double doubleLayer = -2*PI;
      for (size_t m = 0; m < 3; m++) {
	size_t rowIndex = 
	  fp.localStartRowIndex.equa_1 + globalPanelIndex + m*totalNumE_;

	// -Int dG1/dn - 4*PI
	size_t colIndex1 = globalPanelIndex + fp.localStartColIndex.E[m];
	insertSpEle(rowIndex, colIndex1, -doubleLayer - 4*PI);

	// Int G1
	size_t colIndex2 = globalPanelIndex + fp.localStartColIndex.dEdn[m];
	insertSpEle(rowIndex, colIndex2, singleLayer);
      }
    }
  }
}

/**********************************************************************
 * fillEqua1_old_inGlobalCoord --
 **********************************************************************/
void
Preconditioner::fillEqua1_old_inGlobalCoord (
					     void)
{
  for (size_t condIndex = 0; condIndex < fp.numCond_; condIndex++) {
    for (size_t panelIndex = 0; panelIndex < fp.condMesh.numPanel(condIndex); 
	 panelIndex++) {
      size_t globalPanelIndex = fp.condMesh.globalPanelIndex(condIndex, 
							      panelIndex);
      calcp1List[condIndex](fp.allElementList[globalPanelIndex],
			    fp.allElementList[globalPanelIndex]);
      complex<double> singleLayer = calcp1List[condIndex].result(0);
      double doubleLayer = -2*PI;
      for (size_t m = 0; m < 3; m++) {
	size_t rowIndex = fp.globalStartRowIndex.equa_1 + globalPanelIndex + m*totalNumE_;

	// Int G1
	size_t colIndex = globalPanelIndex + fp.globalStartColIndex.dEdn[m];
	insertSpEle(rowIndex, colIndex, singleLayer);

	// -Int dG1/dn - 4*PI
	colIndex = globalPanelIndex + fp.globalStartColIndex.E[m];
	insertSpEle(rowIndex, colIndex, -doubleLayer - 4*PI);
      }
    }
  }
}

/**********************************************************************
 * fillEqua1_new --
 **********************************************************************/
void
Preconditioner::fillEqua1_new (
			       void)
{
  if (! useGlobalCoord_) {
    fillEqua1_new_inLocalCoord();
  } else {
    fillEqua1_new_inGlobalCoord();
  }
}

/**********************************************************************
 * fillEqua1_new_inLocalCoord --
 **********************************************************************/
void
Preconditioner::fillEqua1_new_inLocalCoord (
					    void)
{
  for (size_t condIndex = 0; condIndex < fp.numCond_; condIndex++) {
    for (size_t panelIndex = 0; panelIndex < fp.condMesh.numPanel(condIndex); 
	 panelIndex++) {
      size_t globalPanelIndex = fp.condMesh.globalPanelIndex(condIndex, 
							      panelIndex);
      calcp1List[condIndex](fp.allElementList[globalPanelIndex],
			    fp.allElementList[globalPanelIndex]);
      complex<double> dn_doubleLayer = calcp1List[condIndex].result(0);
      complex<double> dn_singleLayer = 2*PI;

      for (size_t m = 0; m < 3; m++) {
	size_t rowIndex = fp.localStartRowIndex.equa_1 + globalPanelIndex + m*totalNumE_;

	// d/dn Int G1 - 4*PI
	size_t colIndex = globalPanelIndex + fp.localStartColIndex.dEdn[m];
	insertSpEle(rowIndex, colIndex, dn_singleLayer - 4*PI);

	// -d/dn Int dG1/dn
	colIndex = globalPanelIndex + fp.localStartColIndex.E[m];
	insertSpEle(rowIndex, colIndex, -dn_doubleLayer);
      }
    }
  }
}

/**********************************************************************
 * fillEqua1_new_inGlobalCoord --
 **********************************************************************/
void
Preconditioner::fillEqua1_new_inGlobalCoord (
					     void)
{
  for (size_t condIndex = 0; condIndex < fp.numCond_; condIndex++) {
    for (size_t panelIndex = 0; panelIndex < fp.condMesh.numPanel(condIndex); 
	 panelIndex++) {
      size_t globalPanelIndex = fp.condMesh.globalPanelIndex(condIndex, 
							      panelIndex);
      calcp1List[condIndex](fp.allElementList[globalPanelIndex],
			    fp.allElementList[globalPanelIndex]);
      complex<double> dn_doubleLayer = calcp1List[condIndex].result(0);
      complex<double> dn_singleLayer = 2*PI;

      for (size_t m = 0; m < 3; m++) {
	size_t rowIndex = fp.globalStartRowIndex.equa_1 + globalPanelIndex + m*totalNumE_;

	// d/dn Int G1 - 4*PI
	size_t colIndex = globalPanelIndex + fp.globalStartColIndex.dEdn[m];
	insertSpEle(rowIndex, colIndex, dn_singleLayer - 4*PI);

	// -d/dn Int dG1/dn
	colIndex = globalPanelIndex + fp.globalStartColIndex.E[m];
	insertSpEle(rowIndex, colIndex, -dn_doubleLayer);
      }
    }
  }
}

/**********************************************************************
 * fillEqua2a --
 **********************************************************************/
void
Preconditioner::fillEqua2a (
			    void)
{
  if (! useGlobalCoord_) {
    fillEqua2a_inLocalCoord();
  } else { 
    fillEqua2a_inGlobalCoord();
  }
}

/**********************************************************************
 * fillEqua2a_inLocalCoord --
 **********************************************************************/
void
Preconditioner::fillEqua2a_inLocalCoord (
					 void)
{
  size_t numElement = fp.allElementList.size();
  for (size_t elementIndex = 0; elementIndex < numElement; elementIndex++) {
    size_t condIndex = fp.allElementList[elementIndex].boundaryIndex();
    size_t localPanelIndex = fp.allElementList[elementIndex].index();
    if (fp.condMesh.isContact(condIndex, localPanelIndex)) continue;

    if (fp.simuType_ == FULL_WAVE) {
      calcp2_FW(fp.allElementList[elementIndex], fp.allElementList[elementIndex]);
    } else {
      calcp2_EMQS(fp.allElementList[elementIndex], fp.allElementList[elementIndex]);
    }

    complex<double> singleLayer;
    if (fp.simuType_ == FULL_WAVE) {
      singleLayer = calcp2_FW.result(0);
    } else {
      singleLayer = calcp2_EMQS.result(0);
    }
    double doubleLayer = -2*PI;

    size_t rowIndex = elementIndex + fp.localStartRowIndex.equa_2;
    // t1*t1 Int G0 
    size_t colIndex = elementIndex + fp.localStartColIndex.dEdn[0];
    insertSpEle(rowIndex, colIndex, singleLayer);
    // t2*t2 Int G0 
    colIndex = elementIndex + fp.localStartColIndex.dEdn[1];
    insertSpEle(rowIndex + numElement, colIndex, singleLayer);
      
#ifndef DIAG_LOCAL_PRECOND
    // -t1*t1 Int dG0/dn
    colIndex = elementIndex + fp.localStartColIndex.E[0];
    insertSpEle(rowIndex, colIndex, -doubleLayer);
    // -t2*t2 Int dG0/dn
    colIndex = elementIndex + fp.localStartColIndex.E[1];
    insertSpEle(rowIndex + numElement, colIndex, -doubleLayer);
#endif
  }
}

/**********************************************************************
 * fillEqua2a_inGlobalCoord --
 **********************************************************************/
void
Preconditioner::fillEqua2a_inGlobalCoord (
					 void)
{
  size_t numElement = fp.allElementList.size();
  for (size_t elementIndex = 0; elementIndex < numElement; elementIndex++) {
    size_t condIndex = fp.allElementList[elementIndex].boundaryIndex();
    size_t localPanelIndex = fp.allElementList[elementIndex].index();
    if (fp.condMesh.isContact(condIndex, localPanelIndex)) continue;

    if (fp.simuType_ == FULL_WAVE) {
      calcp2_FW(fp.allElementList[elementIndex], fp.allElementList[elementIndex]);
    } else {
      calcp2_EMQS(fp.allElementList[elementIndex], fp.allElementList[elementIndex]);
    }

    complex<double> singleLayer;
    if (fp.simuType_ == FULL_WAVE) {
      singleLayer = calcp2_FW.result(0);
    } else {
      singleLayer = calcp2_EMQS.result(0);
    }

    double doubleLayer = -2*PI;
    
    size_t rowIndex = elementIndex + fp.globalStartRowIndex.equa_2;
    for (size_t m = 0; m < 3; m++) {
      size_t colIndex = elementIndex + fp.globalStartColIndex.dEdn[m];
      // t1 Int G0 
      insertSpEle(rowIndex, colIndex, 
		  singleLayer * fp.allElementList[elementIndex].tangent1(m));
      // t2 Int G0 
      insertSpEle(rowIndex + numElement, colIndex, 
		  singleLayer * fp.allElementList[elementIndex].tangent2(m));
      
      // -t1 Int dG0/dn
      colIndex = elementIndex + fp.globalStartColIndex.E[m];
      insertSpEle(rowIndex, colIndex, 
		  -doubleLayer * fp.allElementList[elementIndex].tangent1(m));
      // -t2 Int dG0/dn
      insertSpEle(rowIndex + numElement, colIndex, 
		  -doubleLayer * fp.allElementList[elementIndex].tangent2(m));
    }
  }
}

/**********************************************************************
 * fillEqua3 --
 **********************************************************************/
void
Preconditioner::fillEqua3 (
			   void)
{
  // self term of Int G0
  for (size_t elementIndex = 0; elementIndex < fp.allElementList.size(); 
       elementIndex++) {

#ifdef DEBUG_RHO
    size_t condIndex = fp.allElementList[elementIndex].boundaryIndex();
    size_t localPanelIndex = fp.allElementList[elementIndex].index();
    if (fp.condMesh.isContact(condIndex, localPanelIndex)) continue;
#endif

    if (fp.simuType_ == FULL_WAVE) {
      calcp2_FW(fp.allElementList[elementIndex], fp.allElementList[elementIndex]);
    } else {
      calcp2_EMQS(fp.allElementList[elementIndex], fp.allElementList[elementIndex]);
    }
    
    complex<double> singleLayer;
    if (fp.simuType_ == FULL_WAVE) {
      singleLayer = calcp2_FW.result(0);
    } else {
      singleLayer = calcp2_EMQS.result(0);
    }

    size_t rowIndex, colIndex;
    if (! useGlobalCoord_) {
      rowIndex = fp.localStartRowIndex.equa_3 + elementIndex;
      colIndex = fp.localStartColIndex.rho + elementIndex;
    } else {
      rowIndex = fp.globalStartRowIndex.equa_3 + elementIndex;
      colIndex = fp.globalStartColIndex.rho + elementIndex;
    }

    insertSpEle(rowIndex, colIndex, singleLayer);
  }     
}

/**********************************************************************
 * copyNonIntegralMat --
 **********************************************************************/
void
Preconditioner::copyNonIntegralEquaMat (
					void)
{
  for (size_t srcRow = 0; srcRow < fp.nonIntegralEquaMat.numRow(); srcRow++) {
    size_t rowIndex = srcRow + fp.globalStartRowIndex.nonIntegralEqua;
    for (size_t srcCol = 0; srcCol < fp.nonIntegralEquaMat[srcRow].size(); 
	 srcCol++) {
      insertSpEle(rowIndex, fp.nonIntegralEquaMat[srcRow].index(srcCol), 
		  fp.nonIntegralEquaMat[srcRow].value(srcCol));
    }
  }
}

/**********************************************************************
 * fillNonIntegralMatInLocalCoord --
 * equation 2a:    4*PI* t dot grad E term   on the non-contact panels
 * equation 2b:    t dot E = 0               on the contact panels
 * equation 6     div E = 0     on the vertices of non-contact panels
 * equation 7     phi = constant  on vertices of contact panels.
 * equation 4                                on non-contact panels
 *  n dot E = j*2*PI*freq*epsilon0/sigma * rho    EMQS or full-wave analysis
 *  n dot E = 0                                   MQS analysis
 *  note:  rho = rho / epsilon0, to be consistant with equa3
 * equation 5     dEn/dn = 0               on contact panels.	 
 * equation 3:     -4*PI*average * phi term  on all panels
 * 
 * equa_2b, equa4, equa5 and equa6 are different in pre-cond and the system
 * matrix, which is define is global coordinate system. Other equations
 * are the same. I could still copy these equations from 
 * formulation.nonIntegralEquaMat.
 **********************************************************************/
void
Preconditioner::fillNonIntegralEquaMatInLocalCoord (
						    void)
{
  for (size_t condIndex = 0; condIndex < fp.condMesh.numCond(); condIndex++) {
    for (size_t panelIndex = 0; panelIndex < fp.condMesh.numPanel(condIndex); 
	 panelIndex++) {
      size_t globalPanelIndex = fp.condMesh.globalPanelIndex(condIndex, panelIndex);

      for (size_t m=0; m<2; m++) {
	// m = 0 for t1
	// m = 1 for t2
	size_t rowIndex = globalPanelIndex + m*fp.condMesh.totalNumPanel() +
	  fp.localStartRowIndex.equa_2;
	if (fp.condMesh.isNonContact(condIndex, panelIndex)) {
	  // equa2a sparse part
	  size_t srcRow = globalPanelIndex + m * fp.condMesh.totalNumPanel();
	  for (size_t srcCol = 0; srcCol < fp.nonIntegralEquaMat[srcRow].size(); 
	       srcCol++) {
	    insertSpEle(rowIndex, fp.nonIntegralEquaMat[srcRow].index(srcCol), 
			fp.nonIntegralEquaMat[srcRow].value(srcCol));
	  }
	} else {
	  // equa2b
	  size_t colIndex = globalPanelIndex + m * fp.condMesh.totalNumPanel();
	  insertSpEle(rowIndex, colIndex, 1.);
	}
      }

      size_t rowIndex = globalPanelIndex + fp.localStartRowIndex.equa_4_5;
      if (fp.condMesh.isNonContact(condIndex, panelIndex)) {
	// euqa4
	size_t colIndex = globalPanelIndex + fp.localStartColIndex.E[2];
	insertSpEle(rowIndex, colIndex, 1.);
      } else {
	// euqa5
	size_t colIndex = globalPanelIndex + fp.localStartColIndex.dEdn[2];
	insertSpEle(rowIndex, colIndex, 1.);
      }

      
      if (fp.simuType_ != MQS) {
	// equa3
	rowIndex = globalPanelIndex + fp.localStartRowIndex.equa_3;
	size_t srcRow = globalPanelIndex + 
	  (fp.globalStartRowIndex.equa_3 - fp.globalStartRowIndex.nonIntegralEqua);
	for (size_t srcCol = 0; srcCol < fp.nonIntegralEquaMat[srcRow].size(); 
	     srcCol++) {
	  insertSpEle(rowIndex, fp.nonIntegralEquaMat[srcRow].index(srcCol), 
		      fp.nonIntegralEquaMat[srcRow].value(srcCol));
	}
      }
    }
  }

  for (size_t globalNodeIndex = 0; globalNodeIndex < fp.condMesh.totalNumNode(); 
       globalNodeIndex++) {
    size_t rowIndex = globalNodeIndex + fp.localStartRowIndex.equa_6_7;
    if ( (fp.condMesh.vertexType(globalNodeIndex) == NON_CONTACT) ||
	 (fp.condMesh.vertexType(globalNodeIndex) == BUFFER) ) {
      // euqa6
      fillEqua6_inLocalCoord(globalNodeIndex, rowIndex);      
    } else {
      // euqa7
      size_t colIndex = globalNodeIndex + fp.localStartColIndex.phi;
      insertSpEle(rowIndex, colIndex, 1.);
    }
  }
}

/**********************************************************************
 * fillEqua6_inLocalCoord --
 * div E = 0     on the vertices of non-contact panels
 * Int Et dot (n cross l) dl - Int dE/dn ds = 0
 * The equation is enforced on the dual panel around each vertex. 
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
Preconditioner::fillEqua6_inLocalCoord (
					const size_t dualPanelIndex, 
					const size_t rowIndex)
{
  size_t numSubPanel = fp.dualPanelList[dualPanelIndex].numSubPanel();
  for (size_t subPanelIndex=0; subPanelIndex<numSubPanel; subPanelIndex++) {
    size_t parentPanelGlobalIndex = 
      fp.dualPanelList[dualPanelIndex].parentPanelGlobalIndex(subPanelIndex);

    vector3D<double> normal = fp.allElementList[parentPanelGlobalIndex].normal();
    vector3D<double> t1 = fp.allElementList[parentPanelGlobalIndex].tangent1();
    vector3D<double> t2 = fp.allElementList[parentPanelGlobalIndex].tangent2();

    point3D node0 = fp.dualPanelList[dualPanelIndex].subPanelVertex(subPanelIndex, 0);
    point3D node1 = fp.dualPanelList[dualPanelIndex].subPanelVertex(subPanelIndex, 1);
    point3D node2 = fp.dualPanelList[dualPanelIndex].subPanelVertex(subPanelIndex, 2);
    vector3D<double> e1 = node0 - node1;
    vector3D<double> e2 = node1 - node2;
    vector3D<double> en1 = crossProd(normal, e1); // without normalization, the desired
    vector3D<double> en2 = crossProd(normal, e2); // edge length is included here.
    double ct1 = (en1 + en2) * t1 / surfConst_.normFactor;
    double ct2 = (en1 + en2) * t2 / surfConst_.normFactor;

    /* Et1 */
    size_t col = fp.localStartColIndex.E[0] + parentPanelGlobalIndex;
    insertSpEle(rowIndex, col, ct1);

    /* Et2 */
    col = fp.localStartColIndex.E[1] + parentPanelGlobalIndex;
    insertSpEle(rowIndex, col, ct2);

    /* dEn/dn */
    double subPanelArea = 
      fp.dualPanelList[dualPanelIndex].subPanelArea(subPanelIndex) / surfConst_.normFactor;
    col = fp.localStartColIndex.dEdn[2] + parentPanelGlobalIndex;
    insertSpEle(rowIndex, col, -subPanelArea);
  }
}

/**********************************************************************
 * fillEqua4_freq_depend --
 **********************************************************************/
void
Preconditioner::fillEqua4_freq_depend (
				       double freq)
{
  for (size_t condIndex = 0; condIndex < fp.numCond_; condIndex++) {
    std::complex<double> coe(0., -1.) ;
    coe *= 2 * PI * freq * surfConst_.EPSILON0 / 
      fp.condInfoPtrList[condIndex]->conductivity();
    for (size_t panelIndex = 0; panelIndex < fp.condMesh.numPanel(condIndex); 
	 panelIndex++) {
      size_t globalPanelIndex = fp.condMesh.globalPanelIndex(condIndex,
							     panelIndex);
      if (fp.condMesh.isNonContact(condIndex, panelIndex)) {

	size_t rowIndex, colIndex;
	if (! useGlobalCoord_) {
	  rowIndex = fp.localStartRowIndex.equa_4_5 + globalPanelIndex;
	  colIndex = fp.localStartColIndex.rho + globalPanelIndex;
	} else {
	  rowIndex = fp.globalStartRowIndex.equa_4_5 + globalPanelIndex;
	  colIndex = fp.globalStartColIndex.rho + globalPanelIndex;
	}

	insertSpEle(rowIndex, colIndex, coe);
      }
    }
  }
}

/**********************************************************************
 * insertSpEle --
 **********************************************************************/
void
Preconditioner::insertSpEle (
			     const size_t rowIndex,
			     const size_t colIndex,
			     const double value)
{
#ifdef DEBUG_PRE_COND
  checkIndex(rowIndex);
  checkIndex(colIndex);
  updateRowBookKeeping(rowIndex);
  updateColBookKeeping(colIndex);
  updateRowSum(rowIndex, value);
  updateColSum(colIndex, value);
#endif

  preCondMat.insertElement(rowIndex, colIndex, value);
}

/**********************************************************************
 * insertSpEle --
 **********************************************************************/
void
Preconditioner::insertSpEle (
			     const size_t rowIndex,
			     const size_t colIndex,
			     const complex<double>& value)
{
#ifdef DEBUG_PRE_COND
  checkIndex(rowIndex);
  checkIndex(colIndex);
  updateRowBookKeeping(rowIndex);
  updateColBookKeeping(colIndex);
  updateRowSum(rowIndex, value);
  updateColSum(colIndex, value);
#endif

  preCondMat.insertElement(rowIndex, colIndex, value);
}

/**********************************************************************
 * updateRowBookKeeping --
 **********************************************************************/
void
Preconditioner::updateRowBookKeeping (
				      size_t row)
{
  rowBookKeeping[row] += 1;
}

/**********************************************************************
 * updateColBookKeeping --
 **********************************************************************/
void
Preconditioner::updateColBookKeeping (
				      size_t col)
{
  colBookKeeping[col] += 1;
}

/**********************************************************************
 * updateRowSum --
 **********************************************************************/
void
Preconditioner::updateRowSum (
			      size_t row,
			      const complex<double>& value)
{
  if (row%2 ==0) {
    rowSum[row] += value;
  } else {
    rowSum[row] -= value;
  }
}

void
Preconditioner::updateRowSum (
			      size_t row,
			      const double value)
{
  if (row%2 ==0) {
    rowSum[row] += value;
  } else {
    rowSum[row] -= value;
  }
}

/**********************************************************************
 * updateColSum --
 **********************************************************************/
void
Preconditioner::updateColSum (
			      size_t col,
			      const complex<double>& value)
{
  if (col%2 ==0) {
    colSum[col] += value;
  } else {
    colSum[col] -= value;
  }
}

void
Preconditioner::updateColSum (
			      size_t col,
			      const double value)
{
  if (col%2 ==0) {
    colSum[col] += value;
  } else {
    colSum[col] -= value;
  }
}

/**********************************************************************
 * checkIndex --
 **********************************************************************/
void
Preconditioner::checkIndex (
			    size_t index)
{
  if (index >= totalNumRow_) {
    errorMessage("preconditioner.cc : checkIndex",
		 "Illegal matrix index, there must be a bug");
  }
}

/**********************************************************************
 * checkRowBookKeeping --
 **********************************************************************/
void
Preconditioner::checkRowBookKeeping (
				     void)
{
  std::ofstream fout("row.tmp");

  for (size_t i=0; i<totalNumRow_; i++) {
    if (rowBookKeeping[i] == 0) {
      cout << "\n\t the problematic row is the " << i << "th row." 
	   << "\n\t total number of rows is" << totalNumRow_ 
	   << "\n\t globalStartRowIndex.equa_1 = " << fp.globalStartRowIndex.equa_1
	   << "\n\t globalStartRowIndex.equa_2 = " << fp.globalStartRowIndex.equa_2
	   << "\n\t globalStartRowIndex.equa_3 = " << fp.globalStartRowIndex.equa_3
	   << "\n\t globalStartRowIndex.equa_4_5 = " << fp.globalStartRowIndex.equa_4_5
	   << "\n\t globalStartRowIndex.equa_6_7 = " << fp.globalStartRowIndex.equa_6_7
	   << endl;
      errorMessage("precondition.cc : checkRowBookKeeping",
		   "Empty row, there must be a bug");
    }
    fout << rowBookKeeping[i] << std::endl;    
  }
}

/**********************************************************************
 * checkColBookKeeping --
 **********************************************************************/
void
Preconditioner::checkColBookKeeping (
				     void)
{
  std::ofstream fout("col.tmp");

  for (size_t i=0; i<totalNumRow_; i++) {
    if (colBookKeeping[i] == 0) {
      cout << "\n\t the problematic col is the " << i << "th col." 
	   << "\n\t total number of cols is" << totalNumRow_ << endl;
      errorMessage("precondition.cc : checkColBookKeeping",
		   "Empty row, there must be a bug");
    }
    fout << colBookKeeping[i] << std::endl;    
  }
}

/**********************************************************************
 * checkRowSum --
 **********************************************************************/
void
Preconditioner::checkRowSum (
			     void)
{
  std::ofstream fout("rowSum.tmp");
  for (size_t i=0; i<totalNumRow_; i++) {
    fout << rowSum[i] << std::endl;    
  }
}

/**********************************************************************
 * checkColSum --
 **********************************************************************/
void
Preconditioner::checkColSum (
			     void)
{
  std::ofstream fout("colSum.tmp");
  for (size_t i=0; i<totalNumRow_; i++) {
    fout << colSum[i] << std::endl;    
  }
}
