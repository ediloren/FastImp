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

  const static char cvsid[] = "$Id: panel.h,v 1.3 2002/07/22 01:46:35 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _PANEL_H_
#define _PANEL_H_

#include "vector3D.h"
#include "meshConst.h"
#include <vector>

namespace mesh {

  class Panel {

  public:

    Panel(void) {}
    Panel(const std::vector<int>& nodeIndex, 
	  const PanelType type = NON_CONTACT, const PanelShape shape = QUAD ) 
      : shape_(shape), type_(type), globalNodeIndex_(nodeIndex) {}
    Panel(const int n1, const int n2, const int n3, 
	  const PanelType type = NON_CONTACT) 
      : shape_(TRIANGLE), type_(type) {
      globalNodeIndex_.resize(3);
      globalNodeIndex_[0] = n1; 
      globalNodeIndex_[1] = n2; 
      globalNodeIndex_[2] = n3;
    }

    PanelShape shape(void) const { return shape_; }
    PanelType type(void) const { return type_; }
    size_t numVertex(void) const { return globalNodeIndex_.size(); }
    size_t globalNodeIndex(size_t localIndex) const { return globalNodeIndex_[localIndex]; }

    bool isNonContact(void) const {
      return ((type_ == NON_CONTACT) || (type_ == BUFFER) );
    }

    bool isContact(void) const {
      return ((type_ == LEFT_CONTACT) || (type_ == RIGHT_CONTACT) || (type_ == CONTACT));
    }

    bool isNonBuffer(void) const {
      return ( (type_ != BUFFER) && (type_ != LEFT_CONTACT) &&
	       (type_ != RIGHT_CONTACT) && (type_ != CONTACT) );
    }

    bool isBuffer(void) const { return (type_ == BUFFER); }

  private:
    std::vector<int> globalNodeIndex_; 
    PanelType type_;
    PanelShape shape_;

  };

} //namespace mesh

#endif
