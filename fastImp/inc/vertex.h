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

  const static char cvsid[] = "$Id: vertex.h,v 1.1.1.1 2002/04/15 21:41:01 bsong Exp $";

  ==========================================================================
*/

#ifndef _VERTEX_H_
#define _VERTEX_H_

#include "vector3D.h"
#include "meshConst.h"

namespace mesh {

  class DualPanel;

  class Vertex {

  public:
    Vertex(void) {}
    Vertex(const pfft::point3D& coord, 
	   const std::vector<size_t>& sharePanelIndex,
	   const mesh::PanelType type) 
      : coordinate_(coord), sharePanelIndex_(sharePanelIndex),
	type_(type) {}

    pfft::point3D coordinate(void) const { return coordinate_; }
    size_t numSharePanel(void) const { return sharePanelIndex_.size(); }
    size_t sharePanelIndex(size_t i) const { return sharePanelIndex_[i]; }
    mesh::PanelType type(void) const { return type_; }

  private:
    pfft::point3D coordinate_;

    // indices of the panels sharing this node, these indices are local to each cond
    std::vector<size_t> sharePanelIndex_; 

    /* this flag shows the type of panel this vertex belongs to.
       A vertex is usually shared by several panels. If one of 
       these panels is not NON_CONTACT type, this flag is same
       as that panel type. So the CANTACT type always supercedes
       the NON_CONTACT type. 
    */
    mesh::PanelType type_;  

  };

} // namespace mesh

#endif
