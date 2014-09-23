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

  A comprehensive driver used for generating those plots on the paper 
  we have published. To get an idea of how to use the pfft++ quickly, please
  refer to driver1.cc and driver2.cc.

  Resources:

  See also:

  const static char cvsid[] = "$Id: testAccuracy.cc,v 1.7 2002/12/26 19:19:20 bsong Exp $";

  ==========================================================================
*/

#include <iostream>
#include <stdio.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <string>
#include <complex>
#include <iomanip>
#include <cstdlib>
//#include "sphereGen.h"
#include "testAccuracy.h"
#include "pfft.h"
#include "element.h"
#include "oneOverR.h"
#include "eikrOverR.h"
#include "staticCollocation.h"
#include "fullwaveCollocation.h"

using namespace std;
using namespace pfft;

static double freqBegin; // -b
static double freqEnd; // -e
static size_t numFreq; // -n
static double maxDirectMem; // -m
static StencilTestType stencilTestType; // -s
static TestType testType; // -t 
static PotentialType testPotentialType; // -p
static char * elementFilename; // -i

static vector<double> freqPnt;
static vector<PotentialType> potentialType; 

const size_t BYTES_PER_MB = 1024 * 1024;

int main(
	 int argc,
	 char *argv[]) 
{
  cmdLineParsor (argc, argv);

  setupPotentialType();
 
  if (testType == ONCEDLP) {
    dynamicTestOnDlpOnly(elementFilename);
    exit(1);
  } else if (testType == MIXED) {
    mixedTest(elementFilename);
  } else {
    if(testType == STATIC) {
      staticTest(elementFilename);
    } else { // testType == DYNAMIC
      setupFreqPnt();
      for (size_t ii=0; ii<numFreq; ii++) {
	dynamicTest (freqPnt[ii], elementFilename);
      }
    }
  }
}

/**********************************************************************
 * staticTest --
 **********************************************************************/
void staticTest (
		 const char* elementFilename)
{
  size_t nPanel;

  // from inner to outer
  // the 1st fold: # of tests with different nearby rule,
  // the 2nd fold: # of tests on slp or dlp,
  vector<vector<double> > reError(potentialType.size());
  vector<vector<double> > mempFFT(potentialType.size());
  vector<vector<double> > timepFFTSetup(potentialType.size());
  vector<vector<double> > timepFFTOneIter(potentialType.size());

  vector<double> memDirect(potentialType.size());
  vector<double> timeDirectSetup(potentialType.size());
  vector<double> timeDirectOneIter(potentialType.size());

  cout << "!!! this is test for 1/r !!!" << endl << endl;

  cout << " Input file = ";
  cout << elementFilename << endl << endl;

  vector<element> srcElementList;
  try {
    string title;
    size_t nCond;      
    readElementFile(elementFilename,
		    srcElementList, title, nCond);    
  } catch (const std::domain_error error) {
    cerr << error.what() << endl;
    cout << "wrong command line option! illegal input file name" << endl;
    exit(1);
  }
  
  nPanel = srcElementList.size();
    
  vector<element>& evalElementList = srcElementList;
    
  accuracyTestForOneOverR (
			   potentialType,
			   srcElementList,
			   evalElementList,
			   memDirect,
			   timeDirectSetup,
			   timeDirectOneIter,
			   mempFFT, 
			   timepFFTSetup,
			   timepFFTOneIter,
			   reError);
    
  
  dumpResultForMatlab (
		       potentialType,
		       "static",
		       nPanel,
		       memDirect,
		       timeDirectSetup,
		       timeDirectOneIter,
		       mempFFT,
		       timepFFTSetup,
		       timepFFTOneIter,
		       reError);
  
}

/**********************************************************************
 * mixedTest --
 **********************************************************************/
void mixedTest (
		const char* elementFilename)
{
  const size_t directStencilSize = 3;
  const size_t interpStencilSize = 1;
  const size_t projStencilSize = 1;

  vector<element> srcElementList;
  try {
    string title;
    size_t nCond;
    readElementFile (
		     elementFilename,
		     srcElementList,
		     title,
		     nCond);
  } catch (
	   const std::domain_error error) {
    cerr << error.what() << endl;
    cout << " wrong command line option! illegal input file name" << endl;
    exit(1);
  }
  
  vector<element>& evalElementList = srcElementList;

  vector<double> realsource;

  const size_t probSize = 2 * srcElementList.size();
  cerr << " # of panels = " << srcElementList.size() << endl;

  char s;
  cerr << " pls key in a choice for test type: " << endl;
  cerr << " u: unit vector " << endl;
  cerr << " r: random vecotr " << endl;
  cerr << " o: all ones vector " << endl;
  
  cin >> s;
  switch(s) {
  case 'u':
    size_t col;
    cerr << " please input the col # " << endl;
    cin >> col;
    realsource.resize(probSize);
    realsource[col] = 1.0;
    break;
  case 'r':
    randomSourceGenerator(probSize, realsource);
    break;
  case 'o':
    realsource.resize(probSize, 1.);
    break;
  default:
    cerr << " unknown type " << endl;
    exit(1);
    break;
  }
  
  vector<DifferentialOperator> outerOperatorList;
  vector<DifferentialOperator> innerOperatorList;
  
  outerOperatorList.push_back(pfft::NONE);
  innerOperatorList.push_back(pfft::NONE);
  
  outerOperatorList.push_back(pfft::NONE);
  innerOperatorList.push_back(pfft::D_DN);
  
  typedef EikrOverR<complex<double> > EikrOverR;
  typedef FullwaveCollocation<complex<double>, element> CalcpEikrOverR;
  typedef Pfft<complex<double>, complex<double>, EikrOverR, CalcpEikrOverR> DynamicPFFT;
  
  vector<complex<double> > cplxsource(probSize);
  for (size_t ii = 0; ii < realsource.size(); ii ++) {
    cplxsource[ii] = realsource[ii];
  }

  const double Pi = 3.1415926535897932;

  const double freq = 1.0e8;
  const double omega = 2 * Pi * freq;

  cerr << " freq = " << freq << endl;
  cerr << " omega = " << omega << endl;

  const double lightspeed= 3.0e8;
  const double K0 = omega / lightspeed;
  const double mu0 = 4.0e-7 * Pi;
  const double sigma = 5.8e7;

  complex<double> KC(K0*K0, -omega * mu0 * sigma);
  KC = sqrt(KC);
  if (real(KC) > 0) {
    KC *= -1;
  } 
  cerr << " K0 = " << K0 << endl;
  cerr << " KC = " << KC << endl;

  EikrOverR eikrOverR(KC);
  CalcpEikrOverR calcpEikrOverR(outerOperatorList, innerOperatorList, KC);
    
  // exp(ikr)/r: direct mv product
  vector<vector<complex<double> > > cplxCoefMatrix;
  vector<complex<double> > dynamicResultDirect (evalElementList.size());
  fillMixedCoefMatrix(calcpEikrOverR, srcElementList, evalElementList, cplxCoefMatrix);
  directProd(cplxsource, cplxCoefMatrix, dynamicResultDirect);
    
  // exp(ikr)/r: pFFT mv product
  DynamicPFFT pfftEikrOverR (
			     srcElementList,
			     evalElementList,
			     eikrOverR,
			     calcpEikrOverR,
			     directStencilSize,
			     projStencilSize,
			     interpStencilSize);
    
  vector<complex<double> > tmp1 (dynamicResultDirect.size());
  vector<complex<double> > tmp2 (dynamicResultDirect.size());
  vector<complex<double> > dynamicResultpFFT(dynamicResultDirect.size());

  pfftEikrOverR.IEoperator(
			   tmp1,
			   outerOperatorList[0],
			   innerOperatorList[0],
			   &cplxsource[0]);

  pfftEikrOverR.IEoperator(
			   tmp2,
			   outerOperatorList[1],
			   innerOperatorList[1],
			   &cplxsource[srcElementList.size()]);

  for (size_t i = 0; i < evalElementList.size(); i++) {
    dynamicResultpFFT[i] = tmp1[i] + tmp2[i];
  }

  double err_mixed
    = relativenorm(dynamicResultDirect, dynamicResultpFFT);
  
  cout << " exp(ikr)/r direct   -   exp(ikr)/r pFFT ";
  cout << err_mixed << endl << endl;

  ofstream fout("err_mixed.dat");
  for (size_t rowIdx = 0; rowIdx < evalElementList.size(); rowIdx ++) {
    double err_mixed
      = abs(dynamicResultDirect[rowIdx] - dynamicResultpFFT[rowIdx]);
    if (abs(dynamicResultDirect[rowIdx]) != 0) {
      err_mixed /= abs(dynamicResultDirect[rowIdx]);
    } else {
      err_mixed = abs(dynamicResultDirect[rowIdx]);
    }

    const char* dlm = "  ";
    
    fout << endl;
    fout << rowIdx << dlm;
    fout << dynamicResultDirect[rowIdx].real() << dlm;
    fout << dynamicResultDirect[rowIdx].imag() << dlm;
    fout << dynamicResultpFFT[rowIdx].real() << dlm;
    fout << dynamicResultpFFT[rowIdx].imag() << dlm;
    fout << err_mixed << dlm;
    fout << endl;
      
  }
  fout.close();
}
  



/**********************************************************************
 * dynamicTestOnDlpOnly --
 **********************************************************************/
void dynamicTestOnDlpOnly (
			   const char* elementFilename)
{
  cout << " this is test for dyanmic kernel only on dlp \n";
  const size_t directStencilSize = 3;
  const size_t interpStencilSize = 1;
  const size_t projStencilSize = 1;

  vector<element> srcElementList;
  try {
    string title;
    size_t nCond;
    readElementFile (
		     elementFilename,
		     srcElementList,
		     title,
		     nCond);
  } catch (
	   const std::domain_error error) {
    cerr << error.what() << endl;
    cout << " wrong command line option! illegal input file name" << endl;
    exit(1);
  }
  
  vector<element>& evalElementList = srcElementList;

  vector<double> realsource;
    
  char s;
  cerr << " pls key in a choice for test type: " << endl;
  cerr << " u: unit vector " << endl;
  cerr << " r: random vecotr " << endl;
  cerr << " o: all ones vector " << endl;
  
  cin >> s;
  switch(s) {
  case 'u':
    size_t col;
    cerr << " please input the col # " << endl;
    cin >> col;
    realsource.resize(srcElementList.size());
    realsource[col] = 1.0;
    break;
  case 'r':
    randomSourceGenerator(srcElementList.size(), realsource);
    break;
  case 'o':
    realsource.resize(srcElementList.size(), 1.);
    break;
  default:
    cerr << " unknown type " << endl;
    exit(1);
    break;
  }
    
  for (size_t ii = 0; ii < potentialType.size(); ii ++) {  
  
    vector<DifferentialOperator> outerOperatorList;
    vector<DifferentialOperator> innerOperatorList;
    
    if (potentialType[ii] == SLP) {
      outerOperatorList.push_back(pfft::NONE);
      innerOperatorList.push_back(pfft::NONE);
    } else {
      outerOperatorList.push_back(pfft::NONE);
      innerOperatorList.push_back(pfft::D_DN);
    }
    
    typedef StaticCollocation<pfft::element> CalcpOneOverR;
    typedef Pfft<double, double, OneOverR, CalcpOneOverR> StaticPfft;
    
    
    OneOverR oneOverR;
    CalcpOneOverR calcpOneOverR(outerOperatorList, innerOperatorList);
    
    // 1/r: direct mv product
    vector<vector<double> > realCoefMatrix;
    vector<double> staticResultDirect (evalElementList.size());
    fillSingleCoefMatrix(calcpOneOverR, srcElementList, evalElementList, 0, realCoefMatrix);
    dumpCoefMatrix("realmat.dat",realCoefMatrix);
    directProd(realsource, realCoefMatrix, staticResultDirect);
    
    // 1/r: pFFT mv product
    StaticPfft pfftOneOverR (
			     srcElementList,
			     evalElementList,
			     oneOverR,
			     calcpOneOverR,
			     directStencilSize,
			     projStencilSize,
			     interpStencilSize);
    
    vector<double> staticResultpFFT (staticResultDirect.size());
    pfftOneOverR.IEoperator (
			     staticResultpFFT,
			     outerOperatorList[0],
			     innerOperatorList[0],
			     realsource);
    
    typedef EikrOverR<double> EikrOverR;
    typedef FullwaveCollocation<complex<double>, element> CalcpEikrOverR;
    typedef Pfft<complex<double>, complex<double>, EikrOverR, CalcpEikrOverR> DynamicPFFT;
    
    vector<complex<double> > cplxsource(realsource.size());
    for (size_t ii = 0; ii < realsource.size(); ii ++) {
      cout << staticResultDirect[ii] << "   " << staticResultpFFT[ii] << endl;
      cplxsource[ii] = realsource[ii];
    }
    
    const double K0 = 0.;
    EikrOverR eikrOverR(K0);
    CalcpEikrOverR calcpEikrOverR(outerOperatorList, innerOperatorList, K0);
    
    // exp(ikr)/r: direct mv product
    vector<vector<complex<double> > > cplxCoefMatrix;
    vector<complex<double> > dynamicResultDirect (evalElementList.size());
    fillSingleCoefMatrix(calcpEikrOverR, srcElementList, evalElementList, 0, cplxCoefMatrix);
    directProd(cplxsource, cplxCoefMatrix, dynamicResultDirect);
    
    // exp(ikr)/r: pFFT mv product
    DynamicPFFT pfftEikrOverR (
			       srcElementList,
			       evalElementList,
			       eikrOverR,
			       calcpEikrOverR,
			       directStencilSize,
			       projStencilSize,
			       interpStencilSize);
    
    vector<complex<double> > dynamicResultpFFT (dynamicResultDirect.size());
    pfftEikrOverR.IEoperator(
			     dynamicResultpFFT,
			     outerOperatorList[0],
			     innerOperatorList[0],
			     cplxsource);
    
    double err_real_direct_pFFT
      = relativenorm(staticResultDirect, staticResultpFFT);
    
    double err_real_direct_cplx_direct
      = relativenorm(staticResultDirect, dynamicResultDirect);
    
    double err_cplx_direct_pFFT 
      = relativenorm(dynamicResultDirect, dynamicResultpFFT);
    
    double err_real_pFFT_cplx_pFFT
      = relativenorm(staticResultpFFT, dynamicResultpFFT);
    
    cout << " 1/r direct   -   1/r pFFT ";
    cout << err_real_direct_pFFT << endl << endl;
    
    cout << " 1/r direct   -   exp(ikr)/r direct ";
    cout << err_real_direct_cplx_direct << endl << endl;
    
    cout << " exp(ikr)/r direct   -   exp(ikr)/r pFFT ";
    cout << err_cplx_direct_pFFT << endl << endl;

    cout << " 1/r pFFT   - exp(ikr)/r pFFT ";
    cout << err_real_pFFT_cplx_pFFT << endl << endl;
    
    ofstream fout("err_real_direct_pFFT.dat");
    for (size_t rowIdx = 0; rowIdx < evalElementList.size(); rowIdx ++) {
      double err_real_direct_pFFT
	= abs(staticResultDirect[rowIdx] - staticResultpFFT[rowIdx]);
      err_real_direct_pFFT /= abs(staticResultDirect[rowIdx]);
      
      double err_cplx_direct_pFFT
	= abs(dynamicResultDirect[rowIdx] - dynamicResultpFFT[rowIdx]);
      err_cplx_direct_pFFT /= abs(dynamicResultDirect[rowIdx]);
      
      fout << endl;
      fout << rowIdx << " ";
      fout << staticResultDirect[rowIdx] << "   ";
      fout << staticResultpFFT[rowIdx] << "   ";
      fout << err_real_direct_pFFT << "   ";
      fout << dynamicResultDirect[rowIdx].real() << "   ";
      fout << dynamicResultDirect[rowIdx].imag() << "   ";
      fout << dynamicResultpFFT[rowIdx].real() << "   ";
      fout << dynamicResultpFFT[rowIdx].imag() << "   ";
      fout << err_cplx_direct_pFFT << "   ";
      fout << endl;
      
      double err1 = abs(realsource[rowIdx] - cplxsource[rowIdx]);
      if (err1 > 1.e-8) {
	cout << " realsource != cplxsource @ ii = " << rowIdx << endl;
	cout << " err1 = " << err1 << endl; 
      } 

      for (size_t colIdx = 0; colIdx < cplxsource.size(); colIdx++) {
	double err2
	  = abs(realCoefMatrix[rowIdx][colIdx] - cplxCoefMatrix[rowIdx][colIdx]) 
	  / abs(cplxCoefMatrix[rowIdx][colIdx]);
	
	if(err2 > 1.e-8 &&
	   abs(cplxCoefMatrix[rowIdx][colIdx]) > 1.e-14) {
	  cout << " realCoefMatrix != cplxCoefMatrix @ (";
	  cout << rowIdx << "," << colIdx <<")";
	  cout << endl << "realCoefMatrix = " << realCoefMatrix[rowIdx][colIdx];
	  cout << "  cplxCoefMatrix = " << cplxCoefMatrix[rowIdx][colIdx] << endl;
	}
      }
    }
    fout.close();
  }
}

/**********************************************************************
 * dynamicTest --
 **********************************************************************/
void dynamicTest (
		  const double freq,
		  const char* elementFilename)
{
  size_t nPanel;

  // from inner to outer
  // the 1st fold: # of tests with different nearby rule,
  // the 2nd fold: # of tests on slp or dlp,
  vector<vector<double> > errorInReal(potentialType.size());
  vector<vector<double> > errorInImag(potentialType.size());
  vector<vector<double> > errorInAmp(potentialType.size());
  vector<vector<double> > errorInPhase(potentialType.size());

  vector<vector<double> > mempFFT(potentialType.size());
  vector<vector<double> > timepFFTSetup(potentialType.size());
  vector<vector<double> > timepFFTOneIter(potentialType.size());
  
  // the 1st fold: # of freq points;
  // the 2nd fold: # of test on different input datum.
  // the 3rd fold: # of test type; at most 2
  vector<double> memDirect(potentialType.size());
  vector<double> timeDirectSetup(potentialType.size());
  vector<double> timeDirectOneIter(potentialType.size());

  cout << "!!! this is test for exp(ikr)/r !!!" << endl << endl;

  cout << " Input filename = ";
  cout << elementFilename << endl << endl;

  string title;
  size_t nCond;
  vector<element> srcElementList;
  
  cout << " now testing with input file = ";
  cout << elementFilename << endl << endl;

  try {
    readElementFile(elementFilename,
		    srcElementList, title, nCond);    
  } catch (const std::domain_error error) {
    cerr << error.what() << endl;
    cout << "wrong command line option! illegal input file name" << endl;
    exit(1);
  }
  
  nPanel = srcElementList.size();
      
  vector<element>& evalElementList = srcElementList;
      
  accuracyTestForEikrOverR (
			    freq,
			    potentialType,
			    srcElementList,
			    evalElementList,
			    memDirect,
			    timeDirectSetup,
			    timeDirectOneIter,
			    mempFFT,
			    timepFFTSetup,
			    timepFFTOneIter,
			    errorInReal,
			    errorInImag,
			    errorInAmp,
			    errorInPhase);
      
  
  dumpResultForMatlab (
		       freq,
		       potentialType,
		       "dynamic_",
		       nPanel,
		       memDirect,
		       timeDirectSetup,
		       timeDirectOneIter,
		       mempFFT,
		       timepFFTSetup,
		       timepFFTOneIter,
		       errorInReal,
		       errorInImag,
		       errorInAmp,
		       errorInPhase);
}

/**********************************************************************
 * dumpResultForMatlab --
 **********************************************************************/
void dumpResultForMatlab (
			  const vector<PotentialType>& potentialType,
			  const string& filenamePrefix,
			  const size_t nPanel,
			  const vector<double>& memDirect,
			  const vector<double>& timeDirectSetup,
			  const vector<double>& timeDirectOneIter,
			  const vector<vector<double> >& mempFFT,
			  const vector<vector<double> >& timepFFTSetup,
			  const vector<vector<double> >& timepFFTOneIter,
			  const vector<vector<double> >& reError)
{
  const string dlm = " ";

  size_t n = 0;
  for (size_t ii = 0; ii < potentialType.size() ; ii ++) {
    string filename(filenamePrefix);
    if (potentialType[ii] == SLP) {
	filename += "_slp.dat";
    } else {
      filename += "_dlp.dat";
    }
    
    ofstream fout(filename.c_str(), std::ios::app);
    if (! fout) {
      // if NOT, abort
      cerr << " fail to open file \"" << filename << "\"" << endl;
      exit (EXIT_FAILURE);
    }

    fout << nPanel << dlm;
    streamsize original_precision = fout.precision(10);
    fout << memDirect[ii] << dlm;
    fout << fabs(timeDirectSetup[ii]) << dlm;
    fout << fabs(timeDirectOneIter[ii]) << dlm;
    for (size_t jj = 0; jj < mempFFT[ii].size(); jj ++) {
      fout << mempFFT[ii][jj] << dlm;
      fout << fabs(timepFFTSetup[ii][jj]) << dlm;
      fout << fabs(timepFFTOneIter[ii][jj]) << dlm;
      fout << reError[ii][jj] << dlm;
    }
    fout << endl;
    fout.precision(original_precision);
    fout.close();
  }
}


/**********************************************************************
 * dumpResultForMatlab --
 **********************************************************************/
void dumpResultForMatlab (
			  const double freq,
			  const vector<PotentialType>& potentialType,
			  const string& filenamePrefix,
			  const size_t nPanel,
			  const vector<double>& memDirect,
			  const vector<double>& timeDirectSetup,
			  const vector<double>& timeDirectOneIter,
			  const vector<vector<double> >& mempFFT,
			  const vector<vector<double> >& timepFFTSetup,
			  const vector<vector<double> >& timepFFTOneIter,
			  const vector<vector<double> >& errorInReal,
			  const vector<vector<double> >& errorInImag,
			  const vector<vector<double> >& errorInAmp,
			  const vector<vector<double> >& errorInPhase)
{
  const string dlm = " ";

  size_t n = 0;
  for (size_t ii = 0; ii < potentialType.size() ; ii ++) {
    string filename(filenamePrefix);
    if (potentialType[ii] == SLP) {
      filename += "_slp.dat";
    } else {
	filename += "_dlp.dat";
    }
  
    ofstream fout(filename.c_str(), std::ios::app);
    if (! fout) {
      // if NOT, abort
      cerr << " fail to open file \"" << filename << "\"" << endl;
      exit (EXIT_FAILURE);
    }

    fout << freq << dlm;
    fout << nPanel << dlm;
    streamsize original_precision = fout.precision(10);
    fout << memDirect[ii] << dlm;
    fout << fabs(timeDirectSetup[ii]) << dlm;
    fout << fabs(timeDirectOneIter[ii]) << dlm;
    for (size_t jj = 0; jj < mempFFT[ii].size(); jj ++) {
      fout << mempFFT[ii][jj] << dlm;
      fout << fabs(timepFFTSetup[ii][jj]) << dlm;
      fout << fabs(timepFFTOneIter[ii][jj]) << dlm;
      fout << errorInReal[ii][jj] << dlm;
      fout << errorInImag[ii][jj] << dlm;
      fout << errorInAmp[ii][jj] << dlm;
      fout << errorInPhase[ii][jj] << dlm;
    }
    fout << endl;
    fout.precision(original_precision);
    fout.close();
  }
}

/**********************************************************************
 * accuracyTestForOneOverR --
 **********************************************************************/
void accuracyTestForOneOverR (
			      const vector<PotentialType>& potentialType,
			      vector<pfft::element>& srcElementList,
			      vector<pfft::element>& evalElementList,
			      vector<double>& memDirect,
			      vector<double>& timeDirectSetup,
			      vector<double>& timeDirectOneIter,
			      vector<vector<double> >& mempFFT,
			      vector<vector<double> >& timepFFTSetup,
			      vector<vector<double> >& timepFFTOneIter,
			      vector<vector<double> >& reError)
{
  // typedef list
  typedef OneOverR GreenFunc;
  typedef StaticCollocation<pfft::element> Calcp;
  typedef Pfft<double, double, GreenFunc, Calcp> PfftType;

  // generate the random source !  
  const size_t nPanels = srcElementList.size();
  vector<double> randomSource(nPanels);
  randomSourceGenerator(nPanels, randomSource);

  bool isTooLargeForDirect;
  for (size_t ii = 0; ii < potentialType.size(); ii++ ) {
    isTooLargeForDirect = memDirectEstimate (nPanels, sizeof(double), memDirect[ii]);
  }

  // generate test stencil size settings.
  vector<size_t> interpStencilSizeVector;
  vector<size_t> projStencilSizeVector;
  vector<size_t> directStencilSizeVector;

  size_t n2;
  // test different directStecilSize
  if (stencilTestType == DIRECT) {
    size_t minStencilSize = 3;
    size_t maxStencilSize = 6;
    n2 = genStencilSize_2 (
			   minStencilSize,
			   maxStencilSize,
			   interpStencilSizeVector,
			   projStencilSizeVector,
			   directStencilSizeVector);
  } else if (stencilTestType == PROJECT) {
      size_t minStencilSize = 1;
      size_t maxStencilSize = 3;
      
      n2 = genStencilSize_3 (
			     minStencilSize,
			     maxStencilSize,
			     interpStencilSizeVector,
			     projStencilSizeVector,
			     directStencilSizeVector);
  } else if (stencilTestType == REGULAR) {
    n2 = 1;
    interpStencilSizeVector.push_back(1);
    projStencilSizeVector.push_back(1);
    directStencilSizeVector.push_back(3);
  } else {
    cerr << " unknown stencil type ! exit ... " << endl;
    exit(1);
  }


  vector<pfft::DifferentialOperator> innerOperatorList;
  vector<pfft::DifferentialOperator> outerOperatorList;
  
  for (size_t ii = 0; ii < potentialType.size(); ii ++) { 
    // operator list 
    

    if (potentialType[ii] == SLP) { // slp
      outerOperatorList.push_back(pfft::NONE);
      innerOperatorList.push_back(pfft::NONE);
    } else { // dlp
      outerOperatorList.push_back(pfft::NONE);
      innerOperatorList.push_back(pfft::D_DN);
    }

    // define kerenl & integration
    GreenFunc kernel;
    Calcp calcp(outerOperatorList, innerOperatorList);

    // direct matrix-vector product
    vector<vector<double> > coefmat;
    vector<double> resultDirect(evalElementList.size());

    if (isTooLargeForDirect) {
      TimeDirectEstimate<double> timeDirectEstimate;
      timeDirectEstimate(nPanels, timeDirectSetup[ii], timeDirectOneIter[ii]);
    } else {
      TNT::stopwatch timer;
      timer.start();

      // 
      fillSingleCoefMatrix(calcp, srcElementList, evalElementList, ii, coefmat);
      timer.stop();
      timeDirectSetup[ii] = timer.read();

      timer.reset();
      timer.start();
      directProd(randomSource, coefmat, resultDirect);
      timer.stop();
      timeDirectOneIter[ii] = timer.read();
    }
    //  mv product accelerated by pFFT
    mempFFT[ii].resize(n2);
    timepFFTOneIter[ii].resize(n2);
    timepFFTSetup[ii].resize(n2);
    reError[ii].resize(n2);


    // loop on different stencil size settings.
    for (size_t jj = 0; jj < n2; jj ++) {
      
      cout << endl 
	   << "directStencil = " << directStencilSizeVector[jj] << endl
	   << "projectStencil = " << projStencilSizeVector[jj] << endl
	   << "interpStencil = " << interpStencilSizeVector[jj] << endl;

      TNT::stopwatch timepFFT;
      timepFFT.reset();
      timepFFT.start();
      PfftType pfft(
		    srcElementList,
		    evalElementList,
		    kernel,
		    calcp,
		    directStencilSizeVector[jj],
		    projStencilSizeVector[jj],
		    interpStencilSizeVector[jj]);

      timepFFT.stop();
      timepFFTSetup[ii][jj] = timepFFT.read();

      mempFFT[ii][jj] = pfft.memoryEstimate() / double(BYTES_PER_MB);
      testEngine (
		  isTooLargeForDirect,
		  pfft,
		  outerOperatorList[ii],
		  innerOperatorList[ii],
		  randomSource,
		  resultDirect,
		  timepFFTOneIter[ii][jj],
		  reError[ii][jj]);

    }
  }

  showResult (
	      potentialType,
	      memDirect,
	      mempFFT,
	      timeDirectSetup,
	      timeDirectOneIter,
	      timepFFTOneIter,
	      timepFFTSetup,
	      interpStencilSizeVector,
	      projStencilSizeVector,
	      directStencilSizeVector,
	      reError,
	      "static test: fixed proj/interp stencil size");
  
  // because zhzhu hardcode the basic type and number 
  // (refer to findBasisFuncTypeAndNumber in Stencil.cc)
  // we cannot do the test 1.
  /* for test 1,2 seperatedly
  vector<vector<double> > reError_1(2); 
  vector<vector<double> > mempFFT_1(2);
  vector<vector<double> > timepFFTOneIter_1(2);
  vector<vector<double> > timepFFTSetup_1(2);

  vector<vector<double> > reError_2(2);
  vector<vector<double> > mempFFT_2(2);
  vector<vector<double> > timepFFTOneIter_2(2);
  vector<vector<double> > timepFFTSetup_2(2);


  size_t minStencilSize;
  size_t maxStencilSize;
  
  vector<size_t> interpStencilSizeVector_1;
  vector<size_t> projStencilSizeVector_1;
  vector<size_t> directStencilSizeVector_1;

  minStencilSize = 1;
  maxStencilSize = 4;
  const size_t n1 = genStencilSize_1 (
				      minStencilSize,
				      maxStencilSize,
				      interpStencilSizeVector_1,
				      projStencilSizeVector_1,
				      directStencilSizeVector_1);
  vector<size_t> interpStencilSizeVector_2;
  vector<size_t> projStencilSizeVector_2;
  vector<size_t> directStencilSizeVector_2;
  
  minStencilSize = 3;
  maxStencilSize = 6;
  const size_t n2 = genStencilSize_2 (minStencilSize,
				      maxStencilSize,
				      interpStencilSizeVector_2,
				      projStencilSizeVector_2,
				      directStencilSizeVector_2);


  // ii = 0: slp; ii = 1: dlp
  for (size_t ii = 0; ii < 2; ii ++) { // 0:slp 1:dlp

    // test 1: fixed the nearby rule
    // ...
    mempFFT_1[ii].resize(n1);
    timepFFTOneIter_1[ii].resize(n1);
    timepFFTSetup_1[ii].resize(n1);
    reError_1[ii].resize(n1);

    for (size_t jj = 0; jj < n1; jj ++) {
      TNT::stopwatch timepFFT;
      cout << " direct stencil size " << directStencilSizeVector_1[jj];
      cout << " interp stencil size " << interpStencilSizeVector_1[jj];
      cout << " proj stencil size " << projStencilSizeVector_1[jj] << endl;

      PfftType pfft(
		    srcElementList,
		    evalElementList,
		    kernel,
		    calcp,
		    directStencilSizeVector_1[jj],
		    projStencilSizeVector_1[jj],
		    interpStencilSizeVector_1[jj]);
      timepFFT.stop();
      timepFFTSetup_1[ii][jj] = timepFFT.read();
      
      mempFFT_1[ii][jj] = pfft.memoryEstimate();

      testEngine (
		  pfft,
		  outerOperatorList[ii],
		  innerOperatorList[ii],
		  randomSource,
		  resultDirect[ii],
		  timepFFTOneIter_1[ii][jj],
		  timepFFTSetup_1[ii][jj],
		  reError_1[ii][jj]);

    }
  
    // test 2: fixed the nearby rule
    // ...
    mempFFT_2[ii].resize(n2);
    timepFFTOneIter_2[ii].resize(n2);
    timepFFTSetup_2[ii].resize(n2);
    reError_2[ii].resize(n1);

    for (size_t jj = 0; jj < n2; jj ++) {
      TNT::stopwatch timepFFT;
      PfftType pfft(
		    srcElementList,
		    evalElementList,
		    kernel,
		    calcp,
		    directStencilSizeVector_2[jj],
		    projStencilSizeVector_2[jj],
		    interpStencilSizeVector_2[jj]);

      timepFFT.stop();
      timepFFTSetup_2[ii][jj] = timepFFT.read();

      mempFFT_2[ii][jj] = pfft.memoryEstimate();
      testEngine (
		  pfft,
		  outerOperatorList[ii],
		  innerOperatorList[ii],
		  randomSource,
		  resultDirect[ii],
		  timepFFTOneIter_2[ii][jj],
		  timepFFTSetup_2[ii][jj],
		  reError_2[ii][jj]);

    }
  }
  */
}

/**********************************************************************
 * accuracyTestForEikrOverR --
 **********************************************************************/
void accuracyTestForEikrOverR (
			       const double freq,
			       const vector<PotentialType>& potentialType,
			       vector<pfft::element>& srcElementList,
			       vector<pfft::element>& evalElementList,
			       vector<double>& memDirect,
			       vector<double>& timeDirectSetup,
			       vector<double>& timeDirectOneIter,
			       vector<vector<double> >& mempFFT,
			       vector<vector<double> >& timepFFTSetup,
			       vector<vector<double> >& timepFFTOneIter,
			       vector<vector<double> >& errorInReal,
			       vector<vector<double> >& errorInImag,
			       vector<vector<double> >& errorInAmp,
			       vector<vector<double> >& errorInPhase)
{
  // comp wavenumber
  // k = w/c; 
  const double pi = 3.1415926535897932;
  const double lightspeed = 3.0e8;
  const double K0 = (2.0*pi*freq)/ lightspeed;

  // typedef list
  typedef EikrOverR<double> GreenFunc;
  typedef FullwaveCollocation<double, element> Calcp;
  typedef Pfft<complex<double>, complex<double>, GreenFunc, Calcp> PfftType;

  // generate the random source !  
  const size_t nPanels = srcElementList.size();
  vector<complex<double> > randomSource(nPanels);
  randomSourceGenerator(nPanels, randomSource);

  bool isTooLargeForDirect;
  for (size_t ii = 0; ii < potentialType.size(); ii++ ) {
    isTooLargeForDirect
      = memDirectEstimate (nPanels, sizeof(complex<double>), memDirect[ii]);
  }

  // generate test stencil size settings.
  vector<size_t> interpStencilSizeVector;
  vector<size_t> projStencilSizeVector;
  vector<size_t> directStencilSizeVector;

  size_t n2;

  // test different directStecilSize
  if (stencilTestType == DIRECT) {
    size_t minStencilSize = 3;
    size_t maxStencilSize = 6;
    n2 = genStencilSize_2 (
			   minStencilSize,
			   maxStencilSize,
			   interpStencilSizeVector,
			   projStencilSizeVector,
			   directStencilSizeVector);
  } else if (stencilTestType == PROJECT) {
    size_t minStencilSize = 1;
    size_t maxStencilSize = 3;
    n2 = genStencilSize_3 (
			   minStencilSize,
			   maxStencilSize,
			   interpStencilSizeVector,
			   projStencilSizeVector,
			   directStencilSizeVector);
  } else if (stencilTestType == REGULAR) {
    n2 = 1;
    interpStencilSizeVector.push_back(1);
    projStencilSizeVector.push_back(1);
    directStencilSizeVector.push_back(3);
  } else {
    cerr << " unknown stencil type ! exit ... " << endl;
    exit(1);
  }

  vector<pfft::DifferentialOperator> innerOperatorList;
  vector<pfft::DifferentialOperator> outerOperatorList;
  for (size_t ii = 0; ii <potentialType.size(); ii ++) { 
    // operator list
    
    if (potentialType[ii] == SLP) { // slp
      outerOperatorList.push_back(pfft::NONE);
      innerOperatorList.push_back(pfft::NONE);
    } else if (potentialType[ii] == DLP) { // dlp
      outerOperatorList.push_back(pfft::NONE);
      innerOperatorList.push_back(pfft::D_DN);
    } else {
      cerr << " unknown potential type ! exit ..." << endl;
      exit(1);
    }

    // define kerenl & integration
    GreenFunc kernel(K0);
    Calcp calcp(outerOperatorList, innerOperatorList, K0);

    // direct m-v product !
    vector<vector<complex<double> > > coefmat;
    vector<complex<double> > resultDirect(evalElementList.size());
    
    if (isTooLargeForDirect) {
      TimeDirectEstimate<complex<double> > timeDirectEstimate;
      timeDirectEstimate(nPanels, timeDirectSetup[ii], timeDirectOneIter[ii]);
    } else {
      TNT::stopwatch timer;
      timer.start();
      fillSingleCoefMatrix(calcp, srcElementList, evalElementList, ii, coefmat);
      timer.stop();
      timeDirectSetup[ii] = timer.read();

      timer.reset();
      timer.start();
      directProd(randomSource, coefmat, resultDirect);
      timer.stop();
      timeDirectOneIter[ii] = timer.read();
    }

    //  mv product accelerated by pFFT
    mempFFT[ii].resize(n2);
    timepFFTSetup[ii].resize(n2);
    timepFFTOneIter[ii].resize(n2);

    errorInReal[ii].resize(n2);
    errorInImag[ii].resize(n2);
    errorInAmp[ii].resize(n2);
    errorInPhase[ii].resize(n2);

    // loop on different stencil size settings.
    for (size_t jj = 0; jj < n2; jj ++) {
      
      TNT::stopwatch timepFFT;
      timepFFT.reset();
      timepFFT.start();
      PfftType pfft(
		    srcElementList,
		    evalElementList,
		    kernel,
		    calcp,
		    directStencilSizeVector[jj],
		    projStencilSizeVector[jj],
		    interpStencilSizeVector[jj]);

      timepFFT.stop();
      timepFFTSetup[ii][jj] = timepFFT.read();

      mempFFT[ii][jj] = pfft.memoryEstimate() / double(BYTES_PER_MB);
      testEngine (
		  isTooLargeForDirect,
		  pfft,
		  outerOperatorList[ii],
		  innerOperatorList[ii],
		  randomSource,
		  resultDirect,
		  timepFFTOneIter[ii][jj],
		  errorInReal[ii][jj],
		  errorInImag[ii][jj],
		  errorInAmp[ii][jj],
		  errorInPhase[ii][jj]);
    }
  }
  
  showResult (
	      freq,
	      potentialType,
	      memDirect,
	      mempFFT,
	      timeDirectSetup,
	      timeDirectOneIter,
	      timepFFTOneIter,
	      timepFFTSetup,
	      interpStencilSizeVector,
	      projStencilSizeVector,
	      directStencilSizeVector,
	      errorInReal,
	      errorInImag,
	      errorInAmp,
	      errorInPhase,
	      " dynamic test: fixed proj/interp stencil size");
}

/**********************************************************************
 * showResult --
 **********************************************************************/
void showResult (
		 const vector<PotentialType>& potentialType,
		 const vector<double>& memDirect,
		 const vector<vector<double> >& mempFFT,
		 const vector<double>& timeDirectSetup,
		 const vector<double>& timeDirectOneIter,
		 const vector<vector<double> >& timepFFTOneIter,
		 const vector<vector<double> >& timepFFTSetup,
		 const vector<size_t>& interpStencilSizeVector,
		 const vector<size_t>& projStencilSizeVector,
		 const vector<size_t>& directStencilSizeVector,
		 const vector<vector<double> >& reError,
		 const string& title)
{
  cout << endl;
  cout << " <---------------------- ";
  cout << title;
  cout << " ----------------------> " << endl;

  const size_t n = mempFFT[0].size();
  
  for (size_t ii = 0; ii < potentialType.size(); ii ++) {
    if (potentialType[ii] == SLP) {
      cout << endl;
      cout << " <---------------------- ";
      cout << " single layer potential ";
      cout << " -----------------------> " << endl << endl;
    } else {
      cout << " <---------------------- ";
      cout << " double layer potential ";
      cout << " -----------------------> " << endl << endl;
    }

    cout << endl;
    cout << " <<<<<<<<<<<<<<<<<<<   ";
    cout << " direct m-v product ";
    cout << "    >>>>>>>>>>>>>>>>>> " << endl << endl;
    
    streamsize original_precision = cout.precision(10);
    cout << " memory requirement: " << memDirect[ii] << " (MB)" << endl;
    cout << " cpu time requirement for setup coef matrix: ";
    cout << timeDirectSetup[ii] << " (s)" << endl;
    cout << " cpu time requirement for one iteration: ";
    cout << timeDirectOneIter[ii] << " (s)" << endl;
    
    cout << endl;
    cout << " <<<<<<<<<<<<<<<<<<<   ";
    cout << " pFFT accelerated m-v product ";
    cout << "    >>>>>>>>>>>>>>>>>> " << endl << endl;
    
    for (size_t jj = 0; jj < n; jj++) {
      cout << endl;
      cout << " stencil size [iterp, proj, direct]: ";
      cout << interpStencilSizeVector[jj] << " , ";
      cout << projStencilSizeVector[jj] << " , ";
      cout << directStencilSizeVector[jj] << endl;
      cout << " relative error: " << reError[ii][jj] << endl;   
      cout << " memory requirement: "<< mempFFT[ii][jj] << " (MB)" << endl;
      cout << " cpu time requirement for setup pFFT: ";
      cout << timepFFTSetup[ii][jj] << endl;
      cout << " cpu time requirement for one iteration: ";
      cout << timepFFTOneIter[ii][jj] << endl;
      cout << endl;
    }

    cout.precision(original_precision);
    cout << endl << endl;
  }
}

/**********************************************************************
 * showResult --
 **********************************************************************/
void showResult (
		 const double freq,
		 const vector<PotentialType>& potentialType,
		 const vector<double>& memDirect,
		 const vector<vector<double> >& mempFFT,
		 const vector<double>& timeDirectSetup,
		 const vector<double>& timeDirectOneIter,
		 const vector<vector<double> >& timepFFTOneIter,
		 const vector<vector<double> >& timepFFTSetup,
		 const vector<size_t>& interpStencilSizeVector,
		 const vector<size_t>& projStencilSizeVector,
		 const vector<size_t>& directStencilSizeVector,
		 const vector<vector<double> >& errorInReal,
		 const vector<vector<double> >& errorInImag,
		 const vector<vector<double> >& errorInAmp,
		 const vector<vector<double> >& errorInPhase,
		 const string& title)
{
  cout << endl;
  cout << " <---------------------- ";
  cout << title;
  cout << " ----------------------> " << endl;

  cout << " freq = " << freq << endl << endl;
  const size_t n = mempFFT[0].size();
  
  for (size_t ii = 0; ii < potentialType.size(); ii ++) {
    if (potentialType[ii] == SLP) {
      cout << endl;
      cout << " <---------------------- ";
      cout << " single layer potential ";
      cout << " -----------------------> " << endl << endl;
    } else {
      cout << " <---------------------- ";
      cout << " double layer potential ";
      cout << " -----------------------> " << endl << endl;
    }

    cout << endl;
    cout << " <<<<<<<<<<<<<<<<<<<   ";
    cout << " direct m-v product ";
    cout << "    >>>>>>>>>>>>>>>>>> " << endl << endl;

    streamsize original_precision = cout.precision(10);

    cout << " memory requirement: " << memDirect[ii] << " (MB)" << endl;
    cout << " cpu time requirement for setup coef matrix: ";
    cout << timeDirectSetup[ii] << " (s)" << endl;
    cout << " cpu time requirement for one iteration: ";
    cout << timeDirectOneIter[ii] << " (s)" << endl;
    
    cout << endl;
    cout << " <<<<<<<<<<<<<<<<<<<   ";
    cout << " pFFT accelerated m-v product ";
    cout << "    >>>>>>>>>>>>>>>>>> " << endl << endl;
    
    for (size_t jj = 0; jj < n; jj++) {
      cout << endl;
      cout << " stencil size [iterp, proj, direct]: ";
      cout << interpStencilSizeVector[jj] << " , ";
      cout << projStencilSizeVector[jj] << " , ";
      cout << directStencilSizeVector[jj] << endl;
      cout << " relative error in real part: " << errorInReal[ii][jj] << endl;
      cout << " relative error in imag part: " << errorInImag[ii][jj] << endl;
      cout << " relative error in amplitude: " << errorInAmp[ii][jj] << endl;
      cout << " relative error in angle: " << errorInPhase[ii][jj] << endl;
      cout << " memory requirement: "<< mempFFT[ii][jj] << " (MB)" << endl;
      cout << " cpu time requirement for setup pFFT: ";
      cout << timepFFTSetup[ii][jj] << endl;
      cout << " cpu time requirement for one iteration: ";
      cout << timepFFTOneIter[ii][jj] << endl;
      cout << endl;
    }

    cout.precision(original_precision);
    cout << endl << endl;
  }
}

/**********************************************************************
 * mydrand --
 **********************************************************************/
inline double mydrand (void) 
{
  //int seed = int(TNT::seconds());
  //void srand(seed);
  return rand()/double(RAND_MAX);
}

/**********************************************************************
 * memDirectEstimate --
 **********************************************************************/
bool memDirectEstimate (
			const size_t n,
			const size_t unitSize,
			double& memInMB)
{
  const size_t BytesPerMB = 1024 * 1024;
  memInMB = n*n * unitSize / double(BytesPerMB);
  bool isTooLargeForDirect = false;
  if (memInMB > maxDirectMem) {
    cout << " warning: expected mem: " 
  	 << memInMB 
  	 << " large than "
	 << maxDirectMem
 	 << endl;
    cout << " skip direct calculation ! " << endl;
    isTooLargeForDirect = true;
  }
      
  // if (memInMB > 30720) {
  //  cout << " warning: expected mem: " <<memInMB<< " large than 30G " << endl;
  //  isTooLargeForDirect = true;
  // } 
  return isTooLargeForDirect;
}

/**********************************************************************
 * genStencilSize_1 --
 **********************************************************************/
size_t genStencilSize_1 (
			 const size_t minStencilSize,
			 const size_t maxStencilSize,
			 std::vector<size_t>& interpStencilSizeVector,
			 std::vector<size_t>& projStencilSizeVector,
			 std::vector<size_t>& directStencilSizeVector)
{
  const size_t n1 = maxStencilSize - minStencilSize + 1;
  
  interpStencilSizeVector.resize(n1);
  projStencilSizeVector.resize(n1);
  directStencilSizeVector.resize(n1);

  for (size_t ii=0; ii<n1; ii++) {
    interpStencilSizeVector[ii] = minStencilSize + ii;
    projStencilSizeVector[ii] = minStencilSize + ii;
    directStencilSizeVector[ii] = interpStencilSizeVector[ii] * 3;
  }
  return n1;
}

/**********************************************************************
 * genStencilSize_2 --
 **********************************************************************/
size_t genStencilSize_2 (
			 const size_t minStencilSize,
			 const size_t maxStencilSize,
			 std::vector<size_t>& interpStencilSizeVector,
			 std::vector<size_t>& projStencilSizeVector,
			 std::vector<size_t>& directStencilSizeVector)
{
  const size_t n2 = maxStencilSize - minStencilSize + 1;
  
  interpStencilSizeVector.resize(n2);
  projStencilSizeVector.resize(n2);
  directStencilSizeVector.resize(n2);

  for (size_t ii=0; ii<n2; ii++) {
    interpStencilSizeVector[ii] = 1;
    projStencilSizeVector[ii] = 1;
    directStencilSizeVector[ii] = minStencilSize + ii;
  }
  
  return n2;
}

/**********************************************************************
 * genStencilSize_3 --
 **********************************************************************/
size_t genStencilSize_3 (
			 const size_t minStencilSize,
			 const size_t maxStencilSize,
			 std::vector<size_t>& interpStencilSizeVector,
			 std::vector<size_t>& projStencilSizeVector,
			 std::vector<size_t>& directStencilSizeVector)
{
  const size_t n2 = maxStencilSize - minStencilSize + 1;
  
  interpStencilSizeVector.resize(n2);
  projStencilSizeVector.resize(n2);
  directStencilSizeVector.resize(n2);

  for (size_t ii=0; ii<n2; ii++) {
    interpStencilSizeVector[ii] = minStencilSize + ii;
    projStencilSizeVector[ii] = minStencilSize + ii;
    directStencilSizeVector[ii] = 3 * projStencilSizeVector[ii];
  }
  
  return n2;
}

/**********************************************************************
 * compareCplxResult -- norm(v1-v2)/norm(v1)
 **********************************************************************/
void compareCplxResult (
			double& errorInReal,
			double& errorInImag,
			double& errorInAmp,
			double& errorInPhase,
			const std::vector<std::complex<double> >& v1,
			const std::vector<std::complex<double> >& v2)
{
  std::vector<double> v1real(v1.size());
  std::vector<double> v1imag(v1.size());
  std::vector<double> v1amp(v1.size());
  std::vector<double> v1phase(v1.size());

  std::vector<double> v2real(v2.size());
  std::vector<double> v2imag(v2.size());
  std::vector<double> v2amp(v2.size());
  std::vector<double> v2phase(v2.size());

  for (size_t i = 0; i < v1.size(); i++) {
    v1real[i] = v1[i].real();
    v2real[i] = v2[i].real();

    v1imag[i] = v1[i].imag();
    v2imag[i] = v2[i].imag();

    v1amp[i] = abs(v1[i]);
    v2amp[i] = abs(v2[i]);
    
    v1phase[i] = arg(v1[i]);
    v2phase[i] = arg(v2[i]);
  }
  
  errorInReal = relativenorm(v1real, v2real);
  errorInImag = relativenorm(v1imag, v2imag);
  errorInAmp = relativenorm(v1amp, v2amp);
  errorInPhase = relativenorm(v1phase, v2phase);
}

/**********************************************************************
 * cmdLineParsor --
 **********************************************************************/
static void
cmdLineParsor (
	       int argc,
	       char *argv[])
{
  int option;

  if (argc == 1) {
    printUsage();	
    exit(0);
  }	

  opterr = 0;  /* Disable stderr Error Message by getopt */
  while ((option = getopt(argc, argv, "i:s:t:p:b:e:n:m:h")) != EOF) {
    switch (option) {
    case INPUT_FILE_NAME:
      elementFilename = optarg;
      break;
    case STENCIL_TYPE:
      stencilTestType = StencilTestType(atoi(optarg));
      break;
    case TEST_TYPE:
      testType = TestType(atoi(optarg));
      break;
    case TEST_POTENTIAL_TYPE:
      testPotentialType = PotentialType(atoi(optarg));
      break;
    case FREQ_BEGIN:
      sscanf(optarg,"%lf", &freqBegin);
      break;
    case FREQ_END:
      sscanf(optarg,"%lf", &(freqEnd));
      break;
    case FREQ_PNT_NUM:
      numFreq = atoi(optarg);
      break;
    case MAX_DIRECT_MEM:
      sscanf(optarg,"%lf", &maxDirectMem);
      break;
    case HELP:
      printUsage();
      exit(0);
      break;
    case '?': /* unknown option */
      printf("\n\n\t Warning in pfft:"
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
static void
printUsage (
	    void)
{
  printf("\n\tUsage:  pfft <options> \n"
	 "\tWhere options are:\n"
	 "\t -i <filename> \n"
	 "\t -t <testtype> 0: static, 1, dynamic, 2 once \n"
	 "\t -s <stenciltype> 0: direct, 1: project, 2: regular \n"
	 "\t -p <potential_type> 0: slp, 1: dlp, 2: both \n"
	 "\t -b <freqBegin> \n"
	 "\t -e <freqEnd> \n"
	 "\t -n <freqNum> \n"
	 "\t -h Print this message and Exit\n"
	 );
}

/**********************************************************************
 * setupPotentialType --
 **********************************************************************/
static inline void setupPotentialType (void)
{
  potentialType.clear();
  if (testPotentialType == BOTH) {
    potentialType.push_back(SLP);
    potentialType.push_back(DLP);
  } else {
    potentialType.push_back(testPotentialType);
  }
}


/**********************************************************************
 * setupFreqPnt --
 **********************************************************************/
static void setupFreqPnt (void)
{
  double fbl = freqBegin;
  double fel = freqEnd;
  double step;
  if (numFreq >= 2) {
    fbl = log10(freqBegin);
    fel = log10(freqEnd);
    step = (fel - fbl)/(numFreq-1);
    step = pow(10, step);
  }

  freqPnt.clear();
  freqPnt.reserve(numFreq);

  size_t n = 0;
  double f;
  while (n < numFreq) {
    if (n == 0) {
      f = pow(10, fbl);
    } else {
      f *= step;
    }
    freqPnt.push_back(f);
    n ++;
  }
}
