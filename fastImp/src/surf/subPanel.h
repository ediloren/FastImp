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

  const static char cvsid[] = "$Id: subPanel.h,v 1.3 2003/05/27 19:13:09 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _SUB_PANEL_H_
#define _SUB_PANEL_H_

#include "vector3D.h"
#include <vector>

namespace surf {

  class SubPanel {

  public:
    SubPanel(void) {}
    SubPanel(const double area, const size_t parentPanelGlobalIndex, 
	     const std::vector<pfft::point3D>& nodeList) 
      : area_(area), parentPanelGlobalIndex_(parentPanelGlobalIndex), 
	       nodeList_(nodeList) {}

    double area(void) const { return area_; }
    size_t parentPanelGlobalIndex(void) const { return parentPanelGlobalIndex_; }
    pfft::point3D vertex(size_t ni) const { return nodeList_[ni]; }

  private:
      double area_; 

      /* Each sub-panel is part of a parent panel, this is the
	 global index of this parent panel 
      */
      size_t parentPanelGlobalIndex_;  

      /* Subpanel is always quadrilateral. One vertex is the center of the dual panel.
	 One vertex is the center of one of the panels sharing the center of the dual
	 panel. The other two vertices are the middle points of the edges of the panel 
	 whose center I just mentioned. 
	 To save memory, only three vertices of the sub-panel are stored. The center 
	 of the dual panel is not stored. These three vertices are in the clock-wise 
	 sequence
      */
      std::vector<pfft::point3D> nodeList_;  
  };

} // namespace surf

#endif
