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

  const static char cvsid[] = "$Id: dualPanel.h,v 1.3 2003/05/27 19:12:59 zhzhu Exp $";

  ==========================================================================
*/


#ifndef _DUAL_PANEL_H_
#define _DUAL_PANEL_H_

#include "vector3D.h"
#include "subPanel.h"

namespace surf {

  class DualPanel {

  public:
    DualPanel(void) {}
    DualPanel(const std::vector<SubPanel> subPanelListIn) 
      : subPanelList(subPanelListIn) {}

    size_t numSubPanel(void) const { return subPanelList.size(); }
    double subPanelArea(size_t i) const { return subPanelList[i].area(); }
    size_t parentPanelGlobalIndex(size_t i) const { 
      return subPanelList[i].parentPanelGlobalIndex(); }
    pfft::point3D subPanelVertex(size_t pi, size_t vi) const { 
      return subPanelList[pi].vertex(vi); }
    void addSubPanel(const SubPanel& sp) { subPanelList.push_back(sp); }

  private:
    /* Around each panel vertex, there is a dual panel. The vertices of this dual 
       panel are the centroids of the panels sharing this specific panel vertex.

       each dual panel is devided by the edges of origianl panels into a few 
       sub-panels. For quadrilateral panels, it usually is 4. At corner it is 3. For 
       trangle panels it could be more than 4. This number is equal to the number of 
       panels sharing this vertex
    */
    std::vector<SubPanel> subPanelList;

  };

} // namespace surf

#endif
