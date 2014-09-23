/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Ben Song
  
  Description:

  Resources:

  See also:

  const static char cvsid[] = "$Id: rectSpiral.cc,v 1.1 2002/08/26 01:55:22 bsong Exp $";

  ==========================================================================
*/

#include "rectSpiral.h"
#include "fileIo.h"

using namespace std;
using namespace surf;

/**********************************************************************
 * RectSpiral::readInputFile --
 **********************************************************************/
void 
RectSpiral::readInputFile (void)
{
  int iFullTurns;
  int iHalfTurns;
  int iQuadTurns;

  Rd_Int(&iFullTurns); fullTurns_ = iFullTurns;
  Rd_Int(&iHalfTurns); halfTurns_ = iHalfTurns;
  Rd_Int(&iQuadTurns); quadTurns_ = iQuadTurns;
  
  width_ = readSize(unit_);
  spacing_ = readSize(unit_);
  thickness_ = readSize(unit_);
  startLineLength_ = readSize(unit_);

  int nw, nt, nl;
  Rd_Int(&nw); numPanelAlongWidth_ = nw;
  Rd_Int(&nt); numPanelAlongThickness_ = nt;
  Rd_Int(&nl); numPanelAlongShortestLength_ = nl;

  Rd_Double(&conductivity_);

  head_ = readPoint(unit_);
  startLineUnitY_ = readPoint(unit_);
  planeNormal_ = readPoint(unit_);

  startLineUnitY_.normalize();
  planeNormal_.normalize();
  
}

