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

  This module computes the contact current. The assumed current direction is
  into the contact for all contacts. This way it is straight
  forward to define the [Z].

  This definition is crucial in deciding the impedance matrix. All subsequent
  calculation should be based on this assumption. 

  Resources:

  See also:

  const static char cvsid[] = "$Id: portCurrent.cc,v 1.23 2003/07/17 16:12:34 zhzhu Exp $";

  ==========================================================================
*/

#include <fstream> 
#include "portCurrent.h"

using namespace std;
using namespace surf;
using namespace mesh;
using namespace pfft;

/**********************************************************************
 * operator () --
 **********************************************************************/
void
PortCurrent::operator () (
			  const ComplexVector& sol,
			  const size_t sCondIndex,
			  const size_t tCondIndex,
			  std::complex<double>& leftCurrent,
			  std::complex<double>& rightCurrent)
{
  switch (currentCompMode) {
  case AUTO_MODE:
    switch (checkFrequencyRange(tCondIndex)) {
    case HIGH_F:
      highFrequencyModeCurrent(sol, sCondIndex, tCondIndex, leftCurrent, rightCurrent);
      break;
    case LOW_F:
      lowFrequencyModeCurrent(sol, sCondIndex, tCondIndex, leftCurrent, rightCurrent);
      break;
    case GREY_WINDOW:
#ifdef NO_GRAY_WINDOW
      lowFrequencyModeCurrent(sol, sCondIndex, tCondIndex, leftCurrent, rightCurrent);
#else
      {
	std::complex<double> lc0, lc1, rc0, rc1;
	lowFrequencyModeCurrent(sol, sCondIndex, tCondIndex, lc0, rc0);
	highFrequencyModeCurrent(sol, sCondIndex, tCondIndex, lc1, rc1);
	leftCurrent = 0.5 * (lc0 + lc1);
	rightCurrent = 0.5 * (rc0 + rc1);
      }
#endif
      break;
    default: 
      errorMessage("portCurrent.cc", "Unknow frequency range !");
    }
    break;
  case LOW_FREQUENCY_MODE:
    lowFrequencyModeCurrent(sol, sCondIndex, tCondIndex, leftCurrent, rightCurrent);
    break;
  case HIGH_FREQUENCY_MODE:
    highFrequencyModeCurrent(sol, sCondIndex, tCondIndex, leftCurrent, rightCurrent);
    break;
  case ENERGY_MODE:
    leftCurrent = energyModeCurrent(sol, sCondIndex, tCondIndex);
    rightCurrent = leftCurrent;
    break;
  case TWO_DIM_MODE:
    //    current = surfCompCurrent2DMode(sol, surfPara, sCondIndex, tCondIndex, voltDiff);
    //    progressReport("Terminal current computation uses 2D mode");
    errorMessage("portCurrent.cc", "2D mode has not been implemented yet !");
    break;
  default: 
    errorMessage("portCurrent.cc", "Illegal currentCompMode !");
  }

#ifdef DEBUG_CONTACT
  outputContactE (tCondIndex, sol);
#endif

#ifdef DEBUG_BUFFER
  outputEAlongLength (tCondIndex, sol);
#endif

}

/**********************************************************************
 * checkFrequencyRange --
 * The true criteria to determine if a frequency is high enough to justify
 * the non-obvious way of computing the current, as shown in function
 * compCurrent_highFrequency, is that the skin depth is on the same order
 * of magnitude as the size of the cross-section of the conductor,
 * such as the width or thickness if the cross-section is rectangular.
 **********************************************************************/
PortCurrent::FrequencyRegime
PortCurrent::checkFrequencyRange (
				  const int condIndex) const
{
  vector<pfft::point3D> leftContactCorner, rightContactCorner;
  findContactCorner(condIndex, leftContactCorner, rightContactCorner);
  double leftMinSize = getMinContactSize(leftContactCorner);
  double rightMinSize = getMinContactSize(rightContactCorner);
  double minSize = min(leftMinSize, rightMinSize);
  double skinDepth = compSkinDepth(condIndex);
  double ratio = minSize / skinDepth;

  // Because each contact has at least four sides. SkinDepth < minSize
  // does not mean the current distribuation is already non-uniform.
  // For a rectangular contact, 0.5<ratio<1.5 might still get a uniform 
  // distribuation. Hence this is a grey period.
  if ( ratio > 1.5) {
    return HIGH_F;
  } if ( ratio < 0.5) {
    return LOW_F;
  } else {
    return GREY_WINDOW;
  }
}

/**********************************************************************
 * skinDepth --
 **********************************************************************/
double
PortCurrent::compSkinDepth (
			    const int condIndex) const 
{
  return sqrt(1./ (surf::PI * condInfoPtrList[condIndex]->conductivity()
		   * frequency * surfConst_.MU0) );
}


/**********************************************************************
 * findContactCorner --
 **********************************************************************/
void
PortCurrent::findContactCorner(
			       const int condIndex,
			       vector<pfft::point3D>& leftContactCorner,
			       vector<pfft::point3D>& rightContactCorner) const
{
  for (int panelIndex=0; panelIndex < condMesh.numPanel(condIndex); panelIndex++) {
    if (IsCornerPanel(condIndex, panelIndex)) {
      if ( condMesh.panelType(condIndex, panelIndex) == LEFT_CONTACT) {
	getCornerVertex(condIndex, panelIndex, leftContactCorner);
      } else if (condMesh.panelType(condIndex, panelIndex) == RIGHT_CONTACT) {
	getCornerVertex(condIndex, panelIndex, rightContactCorner);
      }
    }
  }
}

/**********************************************************************
 * IsCornerPanel --
 **********************************************************************/
bool
PortCurrent::IsCornerPanel (
			    int condIndex, 
			    int panelIndex) const
{
  if (condMesh.isContact(condIndex, panelIndex)) {
    for (int localVertexIndex = 0; 
	 localVertexIndex < condMesh.numVertex(condIndex, panelIndex); 
	 localVertexIndex ++) {
      if ( IsCornerVertex(condIndex, panelIndex, localVertexIndex) ) {
	return true;
      }
    }
    return false;
  } else {
    return false;
  }
}

/**********************************************************************
 * IsCornerVertex --
 **********************************************************************/
bool
PortCurrent::IsCornerVertex (
			     const int condIndex, 
			     const int panelIndex, 
			     const int localVertexIndex) const 
{
  int globalVertexIndex = condMesh.globalNodeIndex(condIndex, panelIndex, 
						   localVertexIndex);
  if ((condMesh.numSharedPanel(globalVertexIndex) == 3) && 
      ( (condMesh.vertexType(globalVertexIndex) == LEFT_CONTACT) ||
        (condMesh.vertexType(globalVertexIndex) == RIGHT_CONTACT) ) ) {
      return true;
  }
  return false;
}

/**********************************************************************
 * getCornerVertex --
 **********************************************************************/
void
PortCurrent::getCornerVertex (
			      const int condIndex, 
			      const int panelIndex,
			      std::vector<pfft::point3D>& corner) const
{
  for (int localVertexIndex = 0; 
       localVertexIndex < condMesh.numVertex(condIndex, panelIndex); 
       localVertexIndex ++) {
    if ( IsCornerVertex(condIndex, panelIndex, localVertexIndex) ) {
      corner.push_back(condMesh.vertex(condIndex, panelIndex, localVertexIndex));
    }
  }  

  if (corner.empty()) {
    errorMessage("surfCompCurrent.c : getCornerVertex",
		 "No corner vertex has been found. There must be a bug here!");
  }
}

/**********************************************************************
 * getMinContactSize --
 **********************************************************************/
double
PortCurrent::getMinContactSize(
			       const std::vector<pfft::point3D>& corner) const
{
  double dist1 = length(corner[0] - corner[1]);
  double dist2 = length(corner[0] - corner[2]);
  double dist3 = length(corner[0] - corner[3]);

  return min(min(dist1, dist2), dist3);
}

/**********************************************************************
 * lowFrequencyModeCurrent --
 * I = Sum_j { sigma * (n dot E) * panelArea_j }
 **********************************************************************/
void
PortCurrent::lowFrequencyModeCurrent (
				      const ComplexVector& sol,
				      const size_t sCondIndex,
				      const size_t tCondIndex,
				      std::complex<double>& leftCurrent,
				      std::complex<double>& rightCurrent) const
{
  leftCurrent = complex<double>(0., 0.);
  rightCurrent = complex<double>(0., 0.);
  size_t numPanel = static_cast<size_t>(condMesh.numPanel(tCondIndex));
  for (size_t panelIndex = 0; panelIndex < numPanel; panelIndex++){
    if (condMesh.panelType(tCondIndex, panelIndex) == LEFT_CONTACT) {
      leftCurrent += compContactPanelCurrent(tCondIndex, panelIndex, sol);
    }
    if (condMesh.panelType(tCondIndex, panelIndex) == RIGHT_CONTACT) {
      rightCurrent += compContactPanelCurrent(tCondIndex, panelIndex, sol);
    }
  }
  adjustImagPartOfcurrent(leftCurrent);
  adjustImagPartOfcurrent(rightCurrent);

  // In calculating low-freq current, the assumed direction for current is
  // out the contact on both left and right contacts. 
  // The assumed universal current direction for each contact is into the contact.
  // So I have to reverse the direction here
  leftCurrent *= -1.;
  rightCurrent *= -1.;
  /*
  if (voltDiff > 0.) {
    // The direction for left current should be reversed if left voltage is higher
    leftCurrent *= -1.;
  } else {
    // The direction for right current should be reversed if right voltage is higher
    rightCurrent *= -1.;
  }
  */
#ifdef DEBUG_CURRENT
  message("Terminal current computation uses low-frequency mode (left and right) :=", 
	  leftCurrent, rightCurrent);
#else
  progressReport("Terminal current computation uses low-frequency mode");
#endif
}

/**********************************************************************
 * adjustImagPartOfcurrent --
 * For certain small structures, at low frequency, we will have R >> wL.
 * So the imag part of current basically becomes the computational noise.
 * There is no point to use this noise to do somethimg useful, for example, 
 * extracting L from an impedance. So I just enforce it here to be zero. 
 * To calculate inductance accurately, the working frequency should be set
 * at middle range such that R and wL are comparable.
 **********************************************************************/
void
PortCurrent::adjustImagPartOfcurrent (
				      std::complex<double>& current) const
{
  const double threshold = 1e-8;
  if ( abs(imag(current) / real(current)) <= threshold ) {
    current = complex<double>(real(current), 0.);
  }
}

/**********************************************************************
 * compContactPanelCurrent --
 **********************************************************************/
std::complex<double>
PortCurrent::compContactPanelCurrent (
				      const size_t condIndex,
				      const size_t panelIndex,
				      const ComplexVector& sol) const
{
  /* I_j = sigma * (n dot E) * panelArea_j */
  complex<double> En = formulation.getEn(condIndex, panelIndex, sol);
  double sigma = condInfoPtrList[condIndex]->conductivity();
  double area = formulation.panelArea(condIndex, panelIndex);
  return En * sigma * area;
}

/**********************************************************************
 * highFrequencyModeCurrent --
 **********************************************************************/
void
PortCurrent::highFrequencyModeCurrent (
				       const ComplexVector& sol,
				       const size_t sCondIndex,
				       const size_t tCondIndex,
				       complex<double>& leftCurrent,
				       complex<double>& rightCurrent)
{
  setupTerminalPanel(tCondIndex);

  leftCurrent = 
    compTerminalCurrent(tCondIndex, leftTerminalPanelIndexList, 
			nextToLeftTerminalPanelIndexList, 
			leftSharedEdgeLengthList, sol);
  rightCurrent = 
    compTerminalCurrent(tCondIndex, rightTerminalPanelIndexList, 
			nextToRightTerminalPanelIndexList, 
			rightSharedEdgeLengthList, sol);

  // In calculating high-freq current, the assumed direction for current is
  // into the contact on both left and right contacts. This is the universal
  // current direction. So no adjustment is needed.
  /*
  if (voltDiff > 0.) {
    // The direction for right current should be reversed if left voltage is higher
    rightCurrent *= -1.;
  } else {
    // The direction for left current should be reversed if right voltage is higher
    leftCurrent *= -1.;
  }
  */

#ifdef DEBUG_CURRENT
  message("Terminal current computation uses high-frequency mode (left and right):=", 
	  leftCurrent, rightCurrent);
#else
  progressReport("Terminal current computation uses high-frequency mode");
#endif
}

/**********************************************************************
 * compTerminalCurrent --
 **********************************************************************/
std::complex<double> 
PortCurrent::compTerminalCurrent (
				  const size_t condIndex, 
				  const vector<int>& terminalPanelIndex, 
				  const vector<int>& nextToTerminalPanelIndex, 
				  const vector<double>& sharedEdgeLength,
				  const ComplexVector& sol) const
{
  double omega = 2*PI*formulation.frequency();
  complex<double> coe(0., 1./(omega * surfConst_.MU0));
  if (formulation.simulationType() == FULL_WAVE) {
    double sigma = condInfoPtrList[condIndex]->conductivity();
    coe /= complex<double>(1, omega * surfConst_.EPSILON0 / sigma);
  } 

  complex<double> current(0., 0.);
  for (size_t pairIndex = 0; pairIndex < terminalPanelIndex.size(); 
       pairIndex++) {
    int panelIndex1 = terminalPanelIndex[pairIndex];
    int panelIndex2 = nextToTerminalPanelIndex[pairIndex];
    pfft::point3D centroid1 = formulation.panelCentroid(condIndex, panelIndex1);
    pfft::point3D centroid2 = formulation.panelCentroid(condIndex, panelIndex2);
    pfft::vector3D<double> tangentForward = centroid2 -centroid1;
    double df = length(tangentForward);
    tangentForward.normalize();

    complex<double> dEdn_tf;
    //    dotProd(dEdn_tf, formulation.getdEdn(condIndex, panelIndex1, sol),
    //    	    tangentForward);
    dotProd(dEdn_tf, formulation.getdEdn(condIndex, panelIndex2, sol),
    	    tangentForward);

    complex<double> dEn_dtf(0., 0.);
    if (formulation.simulationType() != MQS) {
      // For MQS, En = charge = 0.
      complex<double> En1 = formulation.getEn(condIndex, panelIndex1, sol);
      complex<double> En2 = formulation.getEn(condIndex, panelIndex2, sol);
      dEn_dtf = (En2 - En1) / df;
    }
    
    current += (dEn_dtf - dEdn_tf) * sharedEdgeLength[pairIndex] * coe;

#ifdef DEBUG_CURRENT
    cout << dEn_dtf << " " << dEdn_tf << " " << coe 
	 << " " << sharedEdgeLength[pairIndex] << endl;
    cout << current << endl;
#endif
  }

  return current;
}

/**********************************************************************
 * setupTerminalPanel --
 **********************************************************************/
void
PortCurrent::setupTerminalPanel (
				 const size_t condIndex)
{
  int numPanel = condMesh.numPanel(condIndex);
  std::vector<TerminalPanelType> terminalPanelTypeList(numPanel, NORMAL_PANEL);
  setupTerminalPanelTypeList(condIndex, terminalPanelTypeList);

  leftTerminalPanelIndexList.clear();
  rightTerminalPanelIndexList.clear();
  nextToLeftTerminalPanelIndexList.clear();
  nextToRightTerminalPanelIndexList.clear();

  for (int panelIndex = 0; panelIndex < numPanel; panelIndex++) {
    switch (terminalPanelTypeList[panelIndex]) {
    case LEFT_TERMINAL:
      checkPanelShape(condIndex, panelIndex);
      leftTerminalPanelIndexList.push_back(panelIndex);
      break;
    case NEXT_TO_LEFT_TERMINAL:
      checkPanelShape(condIndex, panelIndex);
      nextToLeftTerminalPanelIndexList.push_back(panelIndex);
      break;
    case RIGHT_TERMINAL:
      checkPanelShape(condIndex, panelIndex);
      rightTerminalPanelIndexList.push_back(panelIndex);
      break;
    case NEXT_TO_RIGHT_TERMINAL:
      checkPanelShape(condIndex, panelIndex);
      nextToRightTerminalPanelIndexList.push_back(panelIndex);
      break;
    }
  }  

  if (nextToLeftTerminalPanelIndexList.empty()) {
    // this means the number of panels along the current direction is only 2.
    // So Next-To-Left-Panel is same as right-panel. 
    nextToLeftTerminalPanelIndexList = rightTerminalPanelIndexList;
  }
  if (nextToRightTerminalPanelIndexList.empty()) {
    // this means the number of panels along the current direction is only 2.
    // So Next-To-Right-Panel is same as left-panel
    nextToRightTerminalPanelIndexList = leftTerminalPanelIndexList;
  }

  checkNumPanel();

  leftSharedEdgeLengthList.clear();
  for (size_t i = 0; i < leftTerminalPanelIndexList.size(); i++) {
    double sharedEdgeLength = 
      condMesh.sharedPanelEdgeLength(condIndex,
				     leftTerminalPanelIndexList[i],
				     nextToLeftTerminalPanelIndexList[i]);
    leftSharedEdgeLengthList.push_back(sharedEdgeLength);
  }
  
  rightSharedEdgeLengthList.clear();
  for (size_t i = 0; i < rightTerminalPanelIndexList.size(); i++) {
    double sharedEdgeLength = 
      condMesh.sharedPanelEdgeLength(condIndex,
				     rightTerminalPanelIndexList[i],
				     nextToRightTerminalPanelIndexList[i]);
    rightSharedEdgeLengthList.push_back(sharedEdgeLength);
  }
}

/**********************************************************************
 * checkPanelShape --
 **********************************************************************/
void
PortCurrent::checkPanelShape (
			      const size_t condIndex,
			      const size_t panelIndex)
{
  if (condMesh.panelShape(condIndex, panelIndex) != QUAD) {
    errorMessage("portCurrent.cc",
		 "Sorry, can not handle non-quadrilateral panel!");
  }
}

/**********************************************************************
 * setupTerminalPanelTypeList --
 **********************************************************************/
void
PortCurrent::setupTerminalPanelTypeList(
					const size_t condIndex,
					std::vector<TerminalPanelType>& typeList) const
{
  /* find the panels that touch the contact */
  for (int panelIndex = 0; panelIndex < condMesh.numPanel(condIndex); panelIndex++) {
    if (condMesh.isBuffer(condIndex, panelIndex)) {
      for (int localVertexIndex = 0; 
	   localVertexIndex < condMesh.numVertex(condIndex, panelIndex); 
	   localVertexIndex++) {
	int globalNodeIndex = 
	  condMesh.globalNodeIndex(condIndex, panelIndex, localVertexIndex);
	if (condMesh.vertexType(globalNodeIndex) == LEFT_CONTACT) {
	  typeList[panelIndex] = LEFT_TERMINAL;
	} else if (condMesh.vertexType(globalNodeIndex) == RIGHT_CONTACT) {
	  typeList[panelIndex] = RIGHT_TERMINAL;
	}
      }
    }
  }  

  /* find the panels that touch the terminalPanels */
  for (int panelIndex = 0; panelIndex < condMesh.numPanel(condIndex); panelIndex++) {
    if ( (condMesh.isBuffer(condIndex, panelIndex)) && 
	 (typeList[panelIndex] == NORMAL_PANEL) ) {

      for (int localVertexIndex = 0; 
	   localVertexIndex < condMesh.numVertex(condIndex, panelIndex); 
	   localVertexIndex++) {
	int globalNodeIndex = 
	  condMesh.globalNodeIndex(condIndex, panelIndex, localVertexIndex);
	
	int numSharePanel = condMesh.numSharedPanel(globalNodeIndex);
	for (int i = 0; i < numSharePanel; i++) {
	  int panelIndex2 = condMesh.sharedPanelIndex(globalNodeIndex, i);
	  if (typeList[panelIndex2] == LEFT_TERMINAL) {
	    typeList[panelIndex] = NEXT_TO_LEFT_TERMINAL;
	  } else if (typeList[panelIndex2] == RIGHT_TERMINAL) {
	    typeList[panelIndex] = NEXT_TO_RIGHT_TERMINAL;
	  }	    	    
	}
      }
    }
  }
}

/**********************************************************************
 * checkNumPanel --
 **********************************************************************/
void
PortCurrent::checkNumPanel (
			    void) const
{
  const size_t reasonableNumber = 4;
 
  if ( (leftTerminalPanelIndexList.size() < reasonableNumber) || 
       (rightTerminalPanelIndexList.size() < reasonableNumber) ) {
    errorMessage("portCurrent.cc : checkNumPanel",
		 "Wrong number of terminal panels. There could be a bug here!");
  }

  if ( (nextToLeftTerminalPanelIndexList.size() < reasonableNumber) || 
       (nextToRightTerminalPanelIndexList.size() < reasonableNumber) ) {
    errorMessage("portCurrent.cc : checkNumPanel",
		 "Wrong number of next-to-terminal panels. There could be a bug here!");
  }

  if (leftTerminalPanelIndexList.size() != nextToLeftTerminalPanelIndexList.size()) {
    errorMessage("portCurrent.cc : checkNumPanel",
		 "Bug: Number of left-terminal and next-to-left-terminal panels should be same. This could be caused when the number of panels along current direction is 3");
  }
  
  if (rightTerminalPanelIndexList.size() != nextToRightTerminalPanelIndexList.size()) {
    errorMessage("portCurrent.cc : checkNumPanel",
		 "Bug: Number of right-terminal and next-to-right-terminal panels should be same. This could be caused when the number of panels along current direction is 3");
  }

}

/**********************************************************************
 * energyModeCurrent --
 * Only comupte the equivalent dissipation current. 
 **********************************************************************/
complex<double>
PortCurrent::energyModeCurrent (
				const ComplexVector& sol,
				const size_t sCondIndex,
				const size_t tCondIndex)
{
  if (checkFrequencyRange(tCondIndex) == LOW_F) {
    errorMessage("PortCurrent::energyModeCurrent",
		 "Energy Mode for current computation is invalid because frequency is not high enough.");      
  }

  double powerLoss = 0.;
  for (size_t panelIndex = 0; panelIndex < condMesh.numPanel(tCondIndex); panelIndex++) {
    complex<double> Et1 = formulation.getEt1(tCondIndex, panelIndex, sol);
    complex<double> Et2 = formulation.getEt2(tCondIndex, panelIndex, sol);
    double area = formulation.panelArea(tCondIndex, panelIndex);
    powerLoss += (abs(Et1)*abs(Et1) + abs(Et2)*abs(Et2)) * area;

    complex<double> En = formulation.getEn(tCondIndex, panelIndex, sol);
    vector3D<complex<double> > dEdn = formulation.getdEdn(tCondIndex, panelIndex, sol);
#ifdef DEBUG_CURRENT
    cout << Et1 << Et2 << En << dEdn << endl;
#endif
  }    

  double skinDepth = compSkinDepth(tCondIndex);
  double sigma = condInfoPtrList[tCondIndex]->conductivity();
  double current = powerLoss * skinDepth * sigma / 2. / voltDiff;
  current *= 1e15 / surfConst_.normFactor / surfConst_.normFactor;

  progressReport("Terminal current computation uses energy mode");
  return complex<double>(current);
}

/**********************************************************************
 * outputContactE --
 * Ef_buf : electric field component along the longitudinal direction
 * En_buf : electric field component along the normal direction
 * dEf/dn_buf : 
 * En : E filed normal to the left contact panels
 **********************************************************************/
void
PortCurrent::outputContactE (
			     const int condIndex,
			     const ComplexVector& sol)
{
  // setting up panel pair at terminal again
  setupTerminalPanel(condIndex);

  //  ofstream Ef_File("E.tmp", std::ios::app);
  //  ofstream dEfdn_File("dEdn.tmp", std::ios::app);
  //  ofstream En_File("En.tmp", std::ios::app);
  ofstream Ef_File("Ef_buf.tmp");
  ofstream dEfdn_File("dEdn_buf.tmp");
  ofstream En_File("En_buf.tmp");
  for (size_t pairIndex = 0; pairIndex < leftTerminalPanelIndexList.size(); 
       pairIndex++) {
    int panelIndex1 = leftTerminalPanelIndexList[pairIndex];
    int panelIndex2 = nextToLeftTerminalPanelIndexList[pairIndex];
    pfft::point3D centroid1 = formulation.panelCentroid(condIndex, panelIndex1);
    pfft::point3D centroid2 = formulation.panelCentroid(condIndex, panelIndex2);
    pfft::vector3D<double> tangentForward = centroid2 - centroid1;
    double df = length(tangentForward);
    tangentForward.normalize();

    complex<double> dEdn_tf;
    dotProd(dEdn_tf, formulation.getdEdn_deNormalized(condIndex, panelIndex1, sol),
	    tangentForward);
    dEfdn_File << dEdn_tf << "  " << abs( dEdn_tf) << endl;

    complex<double> E_tf;
    dotProd(E_tf, formulation.getE_deNormalized(condIndex, panelIndex1, sol), 
	    tangentForward);
    Ef_File << E_tf << "  " << abs( E_tf) << endl;

    complex<double> En = formulation.getEn_deNormalized(condIndex, panelIndex1, sol);
    En_File << En << "  " << abs( En) << endl;
  }

  ofstream Jn_File("En_contact.tmp");
  size_t numPanel = static_cast<size_t>(condMesh.numPanel(condIndex));
  for (size_t panelIndex = 0; panelIndex < numPanel; panelIndex++){
    if (condMesh.panelType(condIndex, panelIndex) == LEFT_CONTACT) {
      complex<double> En = formulation.getEn_deNormalized(condIndex, panelIndex, sol);
      double area = formulation.panelArea(condIndex, panelIndex);
      Jn_File << En << "  " << abs( En) << "  " << area << endl;
    }
  }
}

/**********************************************************************
 * outputEAlongLength --
 * En : electric field component along the normal direction
 * dEf/dn : component of dE/dn along the longitudinal direction
 **********************************************************************/
void
PortCurrent::outputEAlongLength (
				 const int condIndex,
				 const ComplexVector& sol)
{
  // setting up panel pair at terminal again
  setupTerminalPanel(condIndex);

  ofstream dEfdn_File("dEdn_len.tmp");
  ofstream En_File("En_len.tmp");
  ofstream phi_File("phi_len.tmp");
  ofstream rho_File("rho_len.tmp");
  size_t numBufferLayer = 5;
  vector<int> panelIndex(numBufferLayer);
  for (size_t pairIndex = 0; pairIndex < leftTerminalPanelIndexList.size(); 
       pairIndex++) {
    panelIndex[0] = leftTerminalPanelIndexList[pairIndex];
    panelIndex[1] = nextToLeftTerminalPanelIndexList[pairIndex];
    // panelIndex[i] is along the length direction
    for (size_t i = 2; i < numBufferLayer; i++) {
      panelIndex[i] = panelIndex[i-1] + panelIndex[i-1] - panelIndex[i-2];
    }

    for (size_t i=0; i < numBufferLayer; i++) {
      size_t j1, j2;
      if ((i==0) || (i==1) ) {
	j1 = 0;
	j2 = 1;
      } else {
	j1 = i - 1;
	j2 = i;
      }
      point3D centroid1 = formulation.panelCentroid(condIndex, panelIndex[j1]);
      point3D centroid2 = formulation.panelCentroid(condIndex, panelIndex[j2]);
      vector3D<double> tangentForward = centroid2 - centroid1;
      tangentForward.normalize();

      complex<double> dEdn_f;
      dotProd(dEdn_f, formulation.getdEdn_deNormalized(condIndex, 
						       panelIndex[i], 
						       sol),
	      tangentForward);
      dEfdn_File << dEdn_f << "  " << abs(dEdn_f) << endl;

      complex<double> En = formulation.getEn_deNormalized(condIndex, 
							  panelIndex[i], sol);
      En_File << En << "  " << abs(En) << endl;

      complex<double> phi = formulation.getPhi(condIndex, panelIndex[i], sol);
      phi_File << phi << "  " << abs(phi) << endl;

      if (formulation.simulationType() != MQS) {
	complex<double> rho = formulation.getRho(condIndex, panelIndex[i], sol);
	rho_File << rho << "  " << abs(rho) << endl;
      }
    }
    dEfdn_File << endl;
    En_File << endl;
    phi_File << endl;
    rho_File << endl;
  }
}
