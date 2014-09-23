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

  const static char cvsid[] = "$Id";

  ==========================================================================
*/

#ifndef _MESH_CONST_H_
#define _MESH_CONST_H_

namespace mesh {

  enum PanelShape { TRIANGLE, QUAD };
  enum PanelType { NON_CONTACT, CONTACT, LEFT_CONTACT, RIGHT_CONTACT, BUFFER } ;
  enum MeshFormat { PATRAN, FASTCAP };

} //namespace mesh

#endif
