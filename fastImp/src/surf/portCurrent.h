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

  const static char cvsid[] = "$Id: portCurrent.h,v 1.9 2003/05/01 18:04:00 zhzhu Exp $";

  ==========================================================================
*/


#ifndef _PORT_CURRENT_H_
#define _PORT_CURRENT_H_

#include <complex>
//#include "mtl/mtl.h"
#include "surfConst.h"
#include "condInfo.h"
#include "mesh.h"
#include "vector3D.h"
#include "formulation.h"

namespace surf {

  class PortCurrent {
    
  public:

    // PortCurrent(void) {;}
    PortCurrent(
		const CurrentCompMode currentCompModeIn,
		const mesh::Mesh& condMeshIn,
		const std::vector<CondInfo*>& condInfoPtrListIn,
		const Formulation& formulationIn,
		const double frequencyIn,
		const double voltDiffIn,
		const SurfConst& surfConst) :  
      currentCompMode(currentCompModeIn), condMesh(condMeshIn), 
      condInfoPtrList(condInfoPtrListIn), formulation(formulationIn),
      frequency(frequencyIn), voltDiff(voltDiffIn), surfConst_(surfConst)
    {};
    
    void operator()(
		    const ComplexVector& sol,
		    const size_t sCondIndex,
		    const size_t tCondIndex,
		    std::complex<double>& leftCurrent,
		    std::complex<double>& rightCurrent);

  private:
    CurrentCompMode currentCompMode;
    // I have to use reference here. Otherwise, copy will
    // be performed in the constructor. And destructors of these
    // classes will be called when the PortCurrent is out of the scope.
    // Since pointers are involved here, this will cause trouble
    // becasue the memory they are pointing to have been freed already.
    const mesh::Mesh& condMesh;
    const std::vector<CondInfo*>& condInfoPtrList;
    const Formulation& formulation;
    double frequency;
    double voltDiff;
    SurfConst surfConst_;

    enum TerminalPanelType { NORMAL_PANEL, 
			    LEFT_TERMINAL, NEXT_TO_LEFT_TERMINAL, 
			    RIGHT_TERMINAL, NEXT_TO_RIGHT_TERMINAL };
    std::vector<int> leftTerminalPanelIndexList;
    std::vector<int> nextToLeftTerminalPanelIndexList;
    std::vector<double> leftSharedEdgeLengthList;
    std::vector<int> rightTerminalPanelIndexList;
    std::vector<int> nextToRightTerminalPanelIndexList;
    std::vector<double> rightSharedEdgeLengthList;
    enum FrequencyRegime { HIGH_F, LOW_F, GREY_WINDOW};
    FrequencyRegime checkFrequencyRange(const int condIndex) const;
    double compSkinDepth (const int condIndex) const;
    void findContactCorner(
			   const int condIndex,
			   std::vector<pfft::point3D>& leftCorner,
			   std::vector<pfft::point3D>& rightCorner) const;
    bool IsCornerPanel (
			int condIndex, 
			int panelIndex) const;
    bool IsCornerVertex (
			 const int condIndex, 
			 const int panelIndex, 
			 const int localVertexIndex) const;
    void getCornerVertex (
			  const int condIndex, 
			  const int panelIndex,
			  std::vector<pfft::point3D>& corner) const;
    double getMinContactSize(const std::vector<pfft::point3D>& corner) const;
    void lowFrequencyModeCurrent (
				  const ComplexVector& sol,
				  const size_t sCondIndex,
				  const size_t tCondIndex,
				  std::complex<double>& leftCurrent,
				  std::complex<double>& rightCurrent) const;
    std::complex<double> 
    compContactPanelCurrent (
			     const size_t condIndex,
			     const size_t panelIndex,
			     const ComplexVector& sol) const;
    void highFrequencyModeCurrent (
				   const ComplexVector& sol,
				   const size_t sCondIndex,
				   const size_t tCondIndex,
				   std::complex<double>& leftCurrent,
				   std::complex<double>& rightCurrent);
    void setupTerminalPanel (
			     const size_t condIndex);
    void setupTerminalPanelTypeList(
				    const size_t condIndex,
				    std::vector<TerminalPanelType>& typeList) const;
    void checkNumPanel (void) const;
    std::complex<double>  
    compTerminalCurrent (
			 const size_t condIndex, 
			 const std::vector<int>& terminalPanelIndex, 
			 const std::vector<int>& nextToTerminalPanelIndex, 
			 const std::vector<double>& sharedEdgeLength,
			 const ComplexVector& sol) const;

    void adjustImagPartOfcurrent (
				  std::complex<double>& current) const;

    std::complex<double>
    PortCurrent::energyModeCurrent(
				   const ComplexVector& sol,
				   const size_t sCondIndex,
				   const size_t tCondIndex);

    void outputContactE (const int condIndex, const ComplexVector& sol);
    void outputEAlongLength (const int condIndex, const ComplexVector& sol);
    void checkPanelShape(const size_t condIndex, const size_t panelIndex);

  };

} // namespace surf

#endif
