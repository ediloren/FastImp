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

  const static char cvsid[] = "$Id: condInfo.cc,v 1.8 2003/07/16 15:32:33 zhzhu Exp $";

  ==========================================================================
*/

#include "condInfo.h"
#include "ring.h"
#include "spiral.h"
#include "wire.h"
#include "ground.h"
#include "rectSpiral.h"
#include "otherShape.h"
#include "fileIo.h"
#include "service.h" // for function errorMessage

using namespace std;
using namespace surf;

/**********************************************************************
 * readCondInfo --
 **********************************************************************/
void
surf::readCondInfo (
		    char* inputFileName,
		    vector<CondInfo*>& condInfoPtrList)
{
  if (! Rd_OpenFile(inputFileName)) {
    errorMessage ("condInfo.cc:",  "Fail to open input file!");
  }

  double unit;
  if (Rd_Double(&unit) == FILE_IO_FALSE) {
    errorMessage ("condInfo.cc:", 
		  "Fail to read unit defined in the input file!");
  }

  int numCond;
  if (Rd_Int(&numCond) == FILE_IO_FALSE) {
    errorMessage ("condInfo.cc:", 
		  "Fail to read numCond defined in the input file!");
  }

  condInfoPtrList.resize(numCond);
  for (int condIndex=0; condIndex < numCond; condIndex++) {

    Rd_BlockStart(""); {

      int condTypeInt;
      if (Rd_Int(&condTypeInt) == FILE_IO_FALSE) {
	errorMessage ("condInfo.cc:", 
		      "Fail to read condType defined in the input file!");
      }

      CondInfo::CondType condType = static_cast<CondInfo::CondType>(condTypeInt);
      switch (condType) {
      case CondInfo::RING: {
	condInfoPtrList[condIndex] = new Ring(condType, unit);
	break;
      }
      case CondInfo::WIRE: {
	condInfoPtrList[condIndex] = new Wire(condType, unit);
	break;
      }
      case CondInfo::SPIRAL: {
	condInfoPtrList[condIndex] = new Spiral(condType, unit);
	break;
      }
      case CondInfo::GROUND: {
	condInfoPtrList[condIndex] = new Ground(condType, unit);
	break;
      }
      case CondInfo::RECTSPIRAL: {
	condInfoPtrList[condIndex] = new RectSpiral(condType, unit);
	break;
      }
      case CondInfo::OTHER_SHAPE: {
	condInfoPtrList[condIndex] = new OtherShape(condType, unit);
	break;
      }
      case CondInfo::OTHER_SHAPE_GROUND: {
	condType = CondInfo::GROUND;
	condInfoPtrList[condIndex] = new OtherShape(condType, unit);
	break;
      }
      default:
	errorMessage ("condInfo.cc:", 
		      "Illegal structure type in input file");
	break;
      }

      condInfoPtrList[condIndex]->readInputFile();

    } Rd_BlockEnd("");
  }

  Rd_CloseFile();
}

/**********************************************************************
 * readSize --
 **********************************************************************/
double
surf::readSize (
		const double unit) 
{
  double size;
  Rd_Double(&size);
#ifdef DISABLE_SCALING
  return (size * unit); // without scaling
#else
  return size; // with scaling
#endif

}

/**********************************************************************
 * readPoint --
 **********************************************************************/
pfft::point3D
surf::readPoint (
		 const double unit) 
{
  double x, y, z;

  Rd_BlockStart(""); {
    x = readSize(unit);
    y = readSize(unit);
    z = readSize(unit);
  } Rd_BlockEnd("");

  return pfft::point3D(x, y, z);
}

