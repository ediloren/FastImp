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

  const static char cvsid[] = "$Id: wirePath.cc,v 1.1 2002/08/26 01:55:22 bsong Exp $";

  ==========================================================================
*/
#include <fstream>
#include "wirePath.h"


using namespace pfft;
using namespace std;
using namespace mesh;

/**********************************************************************
 * WirePath::prnInfo --
 **********************************************************************/
void WirePath::prnInfo (void)
{
  std::ofstream fout("wirepath.dat");
  fout << " <-------   wire path info  --------> ";
  fout << std::endl << std::endl;
  fout << " width of wire path --> " 
	    << this->width_ 
	    << std::endl << std::endl;
  fout << " thickness of wire path --> "
	    << this->thickness_
	    << std::endl << std::endl;
  fout << " # of vtx along wire path --> "
	    << this->vtxArray_.size()
	    << std::endl << std::endl;
  for(size_t i = 0; i< this->vtxArray_.size(); i ++) {
    fout << vtxArray_[i] << std::endl;
  } 
  fout << " # of segments --> "
	    << this->lengthArray_.size()
	    << std::endl << std::endl;
  for (size_t i = 0; i< this->vtxArray_.size() - 1; i ++) {
    fout << " (" << i <<")th segment" << std::endl;
    
    fout << " length --> " 
	      << this->lengthArray_[i]
	      << std::endl;
    fout << " unix Y --> " 
	      << this->localUnitYArray_[i]
	      << std::endl;
    fout << " unix X --> " 
	      << this->localUnitXArray_[i]
	      << std::endl;
  }
  fout << std::endl << std::endl;
  fout.close();
}


