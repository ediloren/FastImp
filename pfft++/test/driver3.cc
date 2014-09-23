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

  A simple driver to pfft++. The kernel is just e^{ikr}/r or d/dn(e^{ikr}/r).
  This drive is to demonstrate how to do frequency sweep efficiently.
  The key here is that you should only construct the frequency-independent
  parts only once.

  When using this driver, please make sure the panel size is at most one tenth
  of a wave length. Input files hfSphere48.qui and hfSphere588.qui 
  satisfy this criteria up to 10GHz.

  This drive only does matrix vector product. It is not a fast solver.

  Resources:

  See also:

  const static char cvsid[] = "$Id: driver3.cc,v 1.4 2003/02/11 03:14:26 zhzhu Exp $";

  ==========================================================================
*/

#include "driver.h"

using namespace std;
using namespace pfft;
using namespace TNT;

char * inputFilename;
KernelType kernelType = SINGLE_LAYER; //default type

/**********************************************************************
 * main --
 **********************************************************************/
int main(
	 int argc,
	 char *argv[]) 
{
  cmdLineParsor (argc, argv);

  vector<element> srcElementList;
  readElementFile(inputFilename, srcElementList);    
  vector<element>& evalElementList = srcElementList;

  vector<pfft::DifferentialOperator> innerOperatorList;
  vector<pfft::DifferentialOperator> outerOperatorList;
  if (kernelType == SINGLE_LAYER) {
    // e^{ikr}/r
    outerOperatorList.push_back(pfft::NONE);
    innerOperatorList.push_back(pfft::NONE);
  } else if (kernelType == DOUBLE_LAYER) {
    // d/dn(e^{ikr}/r)
    outerOperatorList.push_back(pfft::NONE);
    innerOperatorList.push_back(pfft::D_DN);
  }

  TNT::Vector<complex<double> > randomSource(srcElementList.size());
  randomSourceGenerator(randomSource);
  TNT::Matrix<complex<double> > mat(evalElementList.size(), srcElementList.size());
  TNT::Vector<complex<double> > y1(evalElementList.size());
  TNT::Vector<complex<double> > y2(evalElementList.size());

  // this is a dummy calcp with right operator lists
  CalcpTwo calcp(outerOperatorList, innerOperatorList, 0.);
  // This constructs the frequency independent parts.
  PfftTwo pfft(srcElementList, evalElementList, calcp, 4, 2,2);

  const double pi = 3.1415926535897932;
  const double lightspeed = 3.0e8;
  for (double freq = 1e1; freq <= 1e10; freq *= 1e9) {
    // k = w/c; 
    double K0 = (2.0*pi*freq)/ lightspeed;
    DynamicGreenFunc kernel(K0);
    calcp = CalcpTwo(outerOperatorList, innerOperatorList, K0);

    // direct matrix vector product
    fillDenseMatrix(calcp, srcElementList, evalElementList, mat);
    y1 = mat * randomSource;

    // using pfft
    pfft.constructKernelDependent(kernel, calcp);
    pfft.IEoperator(y2, outerOperatorList[0], innerOperatorList[0], randomSource);
    cout << endl 
	 << "frequency =" << freq 
	 << "\t relative error = " << two_norm(y1-y2) / two_norm(y1) 
	 << endl;
  }
}

/**********************************************************************
 * cmdLineParsor --
 **********************************************************************/
void
cmdLineParsor (
	       int argc,
	       char *argv[])
{
  if (argc == 1) {
    printUsage();	
    exit(0);
  }	

  opterr = 0;  /* Disable stderr Error Message by getopt */
  int option;
  while ((option = getopt(argc, argv, "i:k:h")) != EOF) {
    switch (option) {
    case INPUT_FILE_NAME:
      inputFilename = optarg;
      break;
    case KERNEL_TYPE:
      kernelType = KernelType(atoi(optarg));
      break;
    case HELP:
      printUsage();
      exit(0);
      break;
    case '?': /* unknown option */
      printf("\n\n\t Warning in driver3:"
	     "\n\t   Wrong command line options"
	     "\n\t   please check.\n");
      printUsage();
      exit(1);
      break;
    default: /* should never happen */
      printUsage();
      exit(2);
    }
  }
}

/**********************************************************************
 * printUsage --
 **********************************************************************/
void
printUsage (
	    void)
{
  printf("\n\tUsage:  driver3 <options> \n"
	 "\tWhere options are:\n"
	 "\t -i <mesh file name> \n"
	 "\t -k <kernel type> 1: single_layer,   2: double_layer\n"
	 "\t -h Print this message and Exit\n"
	 );
}


