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
  When using this driver, please make sure the panel size is at most one tenth
  of a wave length. Input files hfSphere48.qui and hfSphere588.qui 
  satisfy this criteria up to 10GHz.

  This drive only does matrix vector product. It is not a fast solver.

  Resources:

  See also:

  const static char cvsid[] = "$Id: driver2.cc,v 1.6 2003/02/11 03:14:26 zhzhu Exp $";

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

  // k = w/c; 
  double freq = 1e10;
  const double pi = 3.1415926535897932;
  const double lightspeed = 3.0e8;
  const double K0 = (2.0*pi*freq)/ lightspeed;
  DynamicGreenFunc kernel(K0);
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
  CalcpTwo calcp(outerOperatorList, innerOperatorList, K0);

  // direct matrix vector product
  TNT::Vector<complex<double> > randomSource(srcElementList.size());
  randomSourceGenerator(randomSource);
  TNT::Matrix<complex<double> > mat(evalElementList.size(), srcElementList.size());
  cout << "\t Fill the dense system matrix. It might take a while" << endl; 
  fillDenseMatrix(calcp, srcElementList, evalElementList, mat);
  TNT::Vector<complex<double> > y1(evalElementList.size());
  cout << "\t Calculate matrix vector directly" << endl;
  y1 = mat * randomSource;

  // using pfft
  TNT::Vector<complex<double> > y2(evalElementList.size());
  ofstream fout("error.dat");
  size_t min_size = 1;
  size_t max_size = 2;
  for (size_t size = min_size; size <= max_size; size++) {
    size_t interpStencilSize = size;
    size_t projectStencilSize = size;
    size_t directStencilSize = 2*size;
    cout << endl << "\t Calculate matrix vector product indirectly using pfft" << endl;
    PfftTwo pfft(srcElementList, evalElementList, kernel, calcp,
		 directStencilSize, projectStencilSize, interpStencilSize);
    pfft.IEoperator(y2, outerOperatorList[0], innerOperatorList[0], randomSource);
    cout << endl 
	 << "stencil size = " << directStencilSize
	 << ", " << interpStencilSize
	 << ", " << projectStencilSize << endl
	 << "\t relative error: two_norm = " << two_norm(y1-y2) / two_norm(y1) << endl
	 << endl ;
    fout << endl 
	 << "stencil size = " << directStencilSize
	 << ", " << interpStencilSize
	 << ", " << projectStencilSize << endl
	 << "\t relative error: two_norm = " << two_norm(y1-y2) / two_norm(y1) << endl
	 << endl ;
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


