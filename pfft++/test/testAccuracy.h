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

  const static char cvsid[] = "$Id: testAccuracy.h,v 1.6 2002/12/26 19:19:20 bsong Exp $";

  ==========================================================================
*/

#ifndef __TEST_ACCURACY_H_
#define __TEST_ACCURACY_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include "pfft.h"
#include "element.h"
#include "oneOverR.h"
#include "eikrOverR.h"
#include "calcpForOneOverR.h"
#include "calcpForEikrOverR.h"
#include "staticCollocation.h"
#include "fullwaveCollocation.h"
#include "stpwatch.h"

typedef enum {
  INPUT_FILE_NAME         = 'i',
  STENCIL_TYPE            = 's',
  TEST_TYPE               = 't',
  TEST_POTENTIAL_TYPE     = 'p',
  FREQ_BEGIN              = 'b',
  FREQ_END                = 'e',
  FREQ_PNT_NUM            = 'n',
  MAX_DIRECT_MEM          = 'm',
  HELP                    = 'h'
} cmdLineOptions;


typedef enum {
  STATIC = 0,
  DYNAMIC = 1,
  ONCEDLP = 2,
  MIXED = 3
} TestType;

typedef enum {
  DIRECT = 0,
  PROJECT = 1,
  REGULAR = 2 // (interp/proj/direct) = (1/1/3) !!
} StencilTestType;

typedef enum {
  SLP = 0,
  DLP = 1,
  BOTH = 2,
} PotentialType;

template<class T>
class TimeDirectEstimate;

template<>
class TimeDirectEstimate<double> {
public:
  TimeDirectEstimate() {};
  inline double operator() (const size_t n,
			    double& timeSetup,
			    double& timeOneIter) {
    double ratio = (double(n)/double(base_length));
    timeSetup = pow(ratio,2) * base_time_setup;
    timeOneIter = pow(ratio,2) * base_time_oneIter;
  }
private:
  static const size_t base_length = 10800;
  static const double base_time_setup = 637.94; // for dahlquist
  static const double base_time_oneIter = 1.85; // for dahlquist
};

template<>
const size_t TimeDirectEstimate<double>::base_length;

template<>
const double TimeDirectEstimate<double>::base_time_setup;

template<>
const double TimeDirectEstimate<double>::base_time_oneIter;

template<>
class TimeDirectEstimate<std::complex<double> > {
public:
  TimeDirectEstimate() {};
  inline double operator() (const size_t n,
			    double& timeSetup,
			    double& timeOneIter) {
    double ratio = (double(n)/double(base_length));
    timeSetup = pow(ratio,2) * base_time_setup;
    timeOneIter = pow(ratio,2) * base_time_oneIter;
  }
private:
  static const size_t base_length = 10800;
  static const double base_time_setup = 21225.8; // for dahlquist
  static const double base_time_oneIter = 20.902; // for dahlquist
};

template<>
const size_t TimeDirectEstimate<std::complex<double> >::base_length;

template<>
const double TimeDirectEstimate<std::complex<double> >::base_time_setup;

template<>
const double TimeDirectEstimate<std::complex<double> >::base_time_oneIter;

inline double mydrand(void);

static void cmdLineParsor(
			  int argc, 
			  char *argv[]);
static void
printUsage (
	    void);

static void setupFreqPnt (void);
static inline void setupPotentialType(void);

bool memDirectEstimate (
			const size_t n,
			const size_t unitSize,
			double& memDirect);

void compareCplxResult (
			double& errorInReal,
			double& errorInImag,
			double& errorInAmp,
			double& errorInPhase,
			const std::vector<std::complex<double> >& v1,
			const std::vector<std::complex<double> >& v2);

void dynamicTestOnDlpOnly (
			   const char* elementFileNameList);

size_t genStencilSize_3 (
			 const size_t minStencilSize,
			 const size_t maxStencilSize,
			 std::vector<size_t>& interpStencilSizeVector,
			 std::vector<size_t>& projStencilSizeVector,
			 std::vector<size_t>& directStencilSizeVector);

size_t genStencilSize_2 (
			 const size_t minStencilSize,
			 const size_t maxStencilSize,
			 std::vector<size_t>& interpStencilSizeVector,
			 std::vector<size_t>& projStencilSizeVector,
			 std::vector<size_t>& directStencilSizeVector);

size_t genStencilSize_1 (
			 const size_t minStencilSize,
			 const size_t maxStencilSize,
			 std::vector<size_t>& interpStencilSizeVector,
			 std::vector<size_t>& projStencilSizeVector,
			 std::vector<size_t>& directStencilSizeVector);


void staticTest (
		 const char* elementFileNameList);

void mixedTest (
		 const char* elementFileNameList);

void dynamicTest (
		  const double freq,
		  const char* elementFileNameList);

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
			  const vector<vector<double> >& reError);

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
			  const vector<vector<double> >& errorInPhase);

void accuracyTestForOneOverR (
			      const vector<PotentialType>& potentialType,
			      std::vector<pfft::element>& srcElementList,
			      std::vector<pfft::element>& evalElementList,
			      vector<double>& memDirect,
			      vector<double>& timeDirectSetup,
			      vector<double>& timeDirectOneIter,
			      vector<vector<double> >& mempFFT,
			      vector<vector<double> >& timepFFTSetup,
			      vector<vector<double> >& timepFFTOneIter,
			      vector<vector<double> >& reError);

void accuracyTestForEikrOverR (
			       const double freq,
			       const vector<PotentialType>& potentialType,
			       std::vector<pfft::element>& srcElementList,
			       std::vector<pfft::element>& evalElementList,
			       vector<double>& memDirect,
			       vector<double>& timeDirectSetup,
			       vector<double>& timeDirectOneIter,
			       vector<vector<double> >& mempFFT,
			       vector<vector<double> >& timepFFTSetup,
			       vector<vector<double> >& timepFFTOneIter,
			       vector<vector<double> >& errorInReal,
			       vector<vector<double> >& errorInImag,
			       vector<vector<double> >& errorInAmp,
			       vector<vector<double> >& errorInPhase);

void showResult (
		 const vector<PotentialType>& potentialType,
		 const vector<double>& memDirect,
		 const vector<vector<double> >& memPFFT,
		 const vector<double>& timeDirectSetup,
		 const vector<double>& timeDirectOneIter,
		 const vector<vector<double> >& timePFFTMin,
		 const vector<vector<double> >& timePFFTMax,
		 const vector<size_t>& interpStencilSizeVector,
		 const vector<size_t>& projStencilSizeVector,
		 const vector<size_t>& directStencilSizeVector,
		 const vector<vector<double> >& reError,
		 const string& title);

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
		 const string& title);

/**********************************************************************
 * randomSourceGenerator --
 **********************************************************************/
template <class T>
void randomSourceGenerator (
			    const size_t n,
			    std::vector<T>& v);

template <>
void randomSourceGenerator (
			    const size_t n,
			    std::vector<double>& v)
{
  v.resize(n);
  for (size_t ii = 0; ii < n; ii ++) {
    v[ii] = mydrand();
  }
}

template <>
void randomSourceGenerator (
			    const size_t n,
			    std::vector<std::complex<double> >& v)
{
  v.resize(n);
  for (size_t ii=0; ii < n; ii++) {
    v[ii] = std::complex<double>(mydrand(), mydrand());
  }
}

/**********************************************************************
 * fillMixedCoefMatrix -- fill coef matrix for direct m-v product
 **********************************************************************/
template <
class T,
class CalcpType
>
void fillMixedCoefMatrix (
			  CalcpType& calcp,
			  const std::vector<pfft::element>& srcElementList,
			  const std::vector<pfft::element>& evalElementList,
			  std::vector<std::vector<T> >& coefmat)
{
  const size_t nRow = evalElementList.size();
  const size_t nCol = srcElementList.size();
  coefmat.resize(nRow);
  for (size_t rowIdx = 0; rowIdx < nRow; rowIdx ++) {
    coefmat[rowIdx].resize(2 * nCol);
    for (size_t colIdx = 0; colIdx < nCol; colIdx ++) {
      calcp(srcElementList[colIdx], evalElementList[rowIdx]);
      coefmat[rowIdx][colIdx] = calcp.result(0);
      coefmat[rowIdx][colIdx + nCol] = calcp.result(1);
    }
  }
}

/**********************************************************************
 * fillCoefMatrix -- fill coef matrix for direct m-v product
 **********************************************************************/
template <
class T,
class CalcpType
>
void fillSingleCoefMatrix (
			   CalcpType& calcp,
			   const std::vector<pfft::element>& srcElementList,
			   const std::vector<pfft::element>& evalElementList,
			   const size_t integralIndex,
			   std::vector<std::vector<T> >& coefmat)
{
  const size_t nRow = evalElementList.size();
  const size_t nCol = srcElementList.size();
  coefmat.resize(nRow);
  for (size_t rowIdx = 0; rowIdx < nRow; rowIdx ++) {
    coefmat[rowIdx].resize(nCol);
    for (size_t colIdx = 0; colIdx < nCol; colIdx ++) {
      calcp(srcElementList[colIdx], evalElementList[rowIdx]);
      coefmat[rowIdx][colIdx] = calcp.result(integralIndex);
    }
  }
}

/**********************************************************************
 * directProd -- direct matrix-vector productor
 **********************************************************************/
template <class T>
void directProd (
		 const std::vector<T>& source,
		 const std::vector<vector<T> >& coefmat,
		 std::vector<T>& result)
{
  const size_t nRow = coefmat.size();
  const size_t nCol = coefmat[0].size();
  for (size_t row = 0; row < nRow; row ++) {
    result[row] = T(0.);
    for (size_t col = 0; col < nCol; col ++) {
      result[row] += source[col] * coefmat[row][col];
    }
  }
}

/**********************************************************************
 * relativenorm -- norm(v1-v2)/norm(v1)
 **********************************************************************/
template <class T1, class T2>
inline double relativenorm (
			    const T1& x,
			    const T2& y)
{
  return abs(x - y)/abs(x);
}

/**********************************************************************
 * relativenorm -- norm(v1-v2)/norm(v1)
 **********************************************************************/
template <class T1, class T2>
double relativenorm (
		     const std::vector<T1>& v1,
		     const std::vector<T2>& v2)
{
  if (v1.size() != v2.size()) {
    cerr << " v1 , v2 have different size " << endl;
    exit(1);
  }
  return twonorm(v1,v2) / twonorm(v1);
}

/**********************************************************************
 * twonorm(v1,v2) --
 **********************************************************************/
template <class T1, class T2>
double twonorm (
		const std::vector<T1>& v1,
		const std::vector<T2>& v2)
{
  if (v1.size() != v2.size()) {
    cerr << " v1,v2 have different size " << std::endl;
    exit(1);
  }
  double s = 0;
  for (size_t i = 0; i<v1.size(); i++) {
    s += pow(abs(v1[i] - v2[i]),2);
  }
  return sqrt(s);
}

/**********************************************************************
 * twonorm --
 **********************************************************************/
template <class T>
double twonorm (
		const std::vector<T>& v)
{
  double s = 0;
  for (size_t i = 0; i<v.size(); i++) {
    s += pow(abs(v[i]),2);
  }
  return sqrt(s);
}

/**********************************************************************
 * compareResult -- norm(v1-v2)/norm(v1)
 **********************************************************************/
template <class T1, class T2>
double compareResult (
		      const std::vector<T1>& v1,
		      const std::vector<T2>& v2)
{
  
  double normDiff = 0;
  double normExact = 0;
  
  for (size_t ii=0; ii < v1.size(); ii++) {
    normDiff += pow(abs(v1[ii]-v2[ii]),2);
    normExact += pow(abs(v1[ii]),2);
  }

  return sqrt(normDiff/normExact);
  
}


/**********************************************************************
 * testEngine --
 **********************************************************************/
template <class KernelType, class PfftType>
void testEngine (
		 bool isTooLargeForDirect,
		 PfftType& pfft0,
		 const pfft::DifferentialOperator& outerOperator,
		 const pfft::DifferentialOperator& innerOperator,
		 const std::vector<KernelType>& source,
		 const std::vector<KernelType>& resultDirect,
		 double& timepFFTOneIter,
		 double& reError)
{
  std::vector<typename PfftType::DataType> resultpFFT(resultDirect.size());
  
  // maximum order of proj/interp stencil size is 4^3 = 64
  TNT::stopwatch timer;
  timer.reset();
  timer.start();
  pfft0.IEoperator (
		    resultpFFT,
		    outerOperator, 
		    innerOperator,
		    source);
  timer.stop();

  timepFFTOneIter = timer.read();

  if (isTooLargeForDirect){
    reError = 1.0e-5;
  } else {
    reError = relativenorm(resultDirect, resultpFFT);
  }
}

/**********************************************************************
 * testEngine --
 **********************************************************************/
template <class KernelType, class PfftType>
void testEngine (
		 bool isTooLargeForDirect,
		 PfftType& pfft0,
		 const pfft::DifferentialOperator& outerOperator,
		 const pfft::DifferentialOperator& innerOperator,
		 const std::vector<KernelType>& source,
		 const std::vector<KernelType>& resultDirect,
		 double& timepFFTOneIter,
		 double& errorInReal,
		 double& errorInImag,
		 double& errorInAmp,
		 double& errorInPhase)
{
  std::vector<typename PfftType::DataType> resultpFFT(resultDirect.size());
  
  // maximum order of proj/interp stencil size is 4^3 = 64
  TNT::stopwatch timer;
  timer.reset();
  timer.start();
  pfft0.IEoperator (
		    resultpFFT,
		    outerOperator, 
		    innerOperator,
		    source);
  timer.stop();

  timepFFTOneIter = timer.read();
  
  if (isTooLargeForDirect) {
    errorInReal = errorInImag = errorInPhase = errorInAmp = 1e-6;
  } else {
    compareCplxResult(
		      errorInReal,
		      errorInImag,
		      errorInAmp,
		      errorInPhase,
		      resultDirect,
		      resultpFFT);
  }
}

/**********************************************************************
 * dumpCoefMatrix --
 **********************************************************************/
template<class T>
void dumpCoefMatrix (
		     const char* filename,
		     const std::vector<std::vector<T> >& a)
{
  std::ofstream fout(filename);
  for (size_t ii = 0; ii < a.size(); ii ++) {
    for (size_t jj = 0; jj < a[ii].size(); jj ++) {
      if (jj != a[ii].size() - 1) {
	fout << a[ii][jj] << " ";
      } else {
	fout << a[ii][jj] << std::endl;
      }
    }
  }
  fout.close(); 
}


#endif




