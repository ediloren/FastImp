/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: zhenhai zhu
  
  Description:

  Resources:

  See also:

  const static char cvsid[] = "$Id: driver.h,v 1.4 2003/02/11 03:14:26 zhzhu Exp $";

  ==========================================================================
*/

#ifndef __DRIVER_H_
#define __DRIVER_H_

#include <iostream>
#include <stdio.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <string>
#include <complex>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <float.h>
#include <vector>

#include "pfft.h"
#include "element.h"
#include "oneOverR.h"
#include "eikrOverR.h"
#include "calcpForOneOverR.h"
#include "calcpForEikrOverR.h"
#include "staticCollocation.h"
#include "fullwaveCollocation.h"

typedef pfft::OneOverR StaticGreenFunc;
typedef pfft::StaticCollocation<pfft::element> CalcpOne;
typedef pfft::Pfft<double, double, StaticGreenFunc, CalcpOne> PfftOne;

typedef pfft::EikrOverR<double> DynamicGreenFunc;
typedef pfft::FullwaveCollocation<double, pfft::element> CalcpTwo;
typedef pfft::Pfft<std::complex<double>, std::complex<double>, DynamicGreenFunc, CalcpTwo> PfftTwo;

typedef enum {
  INPUT_FILE_NAME         = 'i',
  KERNEL_TYPE             = 'k',
  HELP                    = 'h'
} cmdLineOptions;

typedef enum {SINGLE_LAYER=1, DOUBLE_LAYER, NUMBER_OF_KERNEL} KernelType;

void cmdLineParsor (int argc, char *argv[]);
void printUsage (void);

/**********************************************************************
 * fillDenseMatrix
 **********************************************************************/
template <class Calcp, class T>
void fillDenseMatrix (
		      Calcp& calcp,
		      const std::vector<pfft::element>& srcElementList,
		      const std::vector<pfft::element>& evalElementList,
		      TNT::Matrix<T>& mat)
{
  const size_t nRow = evalElementList.size();
  const size_t nCol = srcElementList.size();
  for (size_t rowIdx = 0; rowIdx < nRow; rowIdx ++) {
    for (size_t colIdx = 0; colIdx < nCol; colIdx ++) {
      calcp(srcElementList[colIdx], evalElementList[rowIdx]);
      mat[rowIdx][colIdx] = calcp.result(0);
    }
  }
}

/**********************************************************************
 * randomSourceGenerator --
 **********************************************************************/
template <class T>
void 
randomSourceGenerator (
		       TNT::Vector<T>& v)
{
  for (size_t ii = 0; ii < v.size(); ii ++) {
    v[ii] = static_cast<T>(rand() / double(RAND_MAX));
  }
}

#endif
