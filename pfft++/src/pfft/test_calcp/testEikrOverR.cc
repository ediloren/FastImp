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

  const static char cvsid[] = "$Id: testEikrOverR.cc,v 1.5 2003/07/16 18:36:41 zhzhu Exp $";

  ==========================================================================
*/
#include <vector>
#include <complex>
#include <string>
#include "calcpForEikrOverR.h"
#include "calcpForOneOverR.h"
#include "element.h"
#include "vector3D.h"
#include "utils.h" // for PAI

using namespace std;
using namespace pfft;

void 
setSrcPanelEvalPnt (
		    element& srcPanel,
		    vector<point3D>& evalPntVector,
		    vector<string>& evalPntDesList,
		    const double skinDepth);
std::complex<double>
compWaveNumber1 (
		 const double f,
		 double& delta);

double
compWaveNumber2 (
		 const double f);

/**********************************************************************
 * main --
 **********************************************************************/
int
main (
      void)
{
  element srcPanel;
  vector<point3D> evalPntList;
  vector<string> evalPntDesList;

  double freq = 1e10;
  double delta;
#ifdef TEST_FREESPACE
  double k = compWaveNumber2(freq);
  cout << "k := " << k << endl;
#else
  complex<double> k = compWaveNumber1(freq, delta);
  cout << "k := " << k << endl;
#endif
  setSrcPanelEvalPnt(srcPanel, evalPntList, evalPntDesList, delta);

  bool needSlp = true;
  bool needDlp = true;
  bool needGradSlp = true;
  int numQuadPoint = 24;
#ifdef TEST_FREESPACE
  CalcpForEikrOverR<double, element> 
#else
  CalcpForEikrOverR<complex<double>, element> 
#endif
    calcp(k, FOLLOW_NORMAL_DIRECTION, 
	  needSlp, needDlp, needGradSlp, numQuadPoint);

  complex<double> slp = 0;
  complex<double> dlp = 0;
  vector3D<complex<double> > gradSlp(0., 0., 0.);
  for (size_t i=0; i<evalPntList.size(); i++) {
    calcp(srcPanel, evalPntList[i], slp, dlp, gradSlp);
    cout << endl
	 << evalPntDesList[i] << endl
	 << "\tslp = " << slp << endl
	 << "\tdlp = " << dlp << endl
	 << "\tgradSlp = " << gradSlp
	 << endl;
  }
  return 1;
}

/**********************************************************************
 * setSrcPanelEvalPnt --
 **********************************************************************/
void 
setSrcPanelEvalPnt (
		    element& srcPanel,
		    vector<point3D>& evalPntVector,
		    vector<string>& evalPntDesList,
		    const double skinDepth)
{
#ifdef TEST_FREESPACE
  double a = 0.1e-3;
  double b = a*0.1;
  double h = 0.;
#else
  double a = 0.1e-3;
  //  double a = 14*skinDepth;
  double b = 1.*skinDepth;
  double h = 0.;
#endif

  point3D v1(0., 0., 0.);
  point3D v2(a, 0., 0.);
  point3D v3(a, b, 0.);
  point3D v4(0., b ,0.);  
  srcPanel = element(v1, v2, v3, v4);
  cout << "panel centroid := " << srcPanel.centroid() << endl;

  evalPntVector.push_back(point3D(srcPanel.centroid())); 
  evalPntDesList.push_back(" center ");

  evalPntVector.push_back(point3D(1.5*a, b/2, 0.));
  evalPntDesList.push_back(" 1.5 side length away ");

  evalPntVector.push_back(point3D(a, b/2, b/2));
  evalPntDesList.push_back(" above edge ");
}

/**********************************************************************
 * compWaveNumber1 --
 **********************************************************************/
std::complex<double>
compWaveNumber1 (
		 const double f,
		 double& skinDepth)
{
  double sigma = 5.8e7; // that of copper
  double mu = 4e-7 * PAI;
  double w = 2. * PAI * f;
  skinDepth = sqrt(2./(w*mu*sigma));
  complex<double> kc = complex<double>(-1., 1.) / skinDepth;
  return kc;
}

/**********************************************************************
 * compWaveNumber2 --
 **********************************************************************/
double
compWaveNumber2 (
		 const double f)
{
  double w = 2. * PAI * f;
  double k = w / 3e8;
  return k;
}
