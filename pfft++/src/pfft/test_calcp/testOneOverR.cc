/*-*-c++-*-*/
#include <vector>
#include <complex>
#include <string>
#include "calcpForEikrOverR.h"
#include "calcpForOneOverR.h"
#include "element.h"
#include "vector3D.h"

using namespace std;
using namespace pfft;

void 
setSrcPanelEvalPnt (
		    element& srcPanel,
		    vector<point3D>& evalPntVector,
		    vector<string>& evalPntDesVector);

/**********************************************************************
 * main --
 **********************************************************************/
int
main (
      void)
{
  element srcPanel;
  vector<point3D> evalPntList;
  vector<string> evalPntDesVector;
  setSrcPanelEvalPnt(srcPanel, evalPntList, evalPntDesVector);

  bool needSlp = true;
  bool needDlp = true;
  bool needGradSlp = true;
  bool needGradDlp = true;
  CalcpForOneOverR<element> calcp(needSlp, needDlp, needGradSlp, needGradDlp);

  double slp = 0;
  double dlp = 0;
  vector3D<double> gradSlp(0., 0., 0.);
  vector3D<double> gradDlp(0., 0., 0.);
  for (size_t i=0; i<evalPntList.size(); i++) {
    calcp(srcPanel, evalPntList[i], slp, dlp, gradSlp, gradDlp);
    cout << endl
	 << evalPntDesVector[i] << endl
	 << "\tslp = " << slp << endl
	 << "\tdlp = " << dlp << endl
	 << "\tgradSlp = " << gradSlp
	 << "\tgradDlp = " << gradDlp;
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
		    vector<string>& evalPntDesVector)
{
  double a = 0.1e-3;
  double b = a*0.1;
  double h = 0.;

  point3D v1(0., 0., 0.);
  point3D v2(a, 0., 0.);
  point3D v3(a, b, 0.);
  point3D v4(0., b ,0.);  
  srcPanel = element(v1, v2, v3, v4);

  evalPntVector.push_back(point3D(a/2, b/2, 0.)); 
  evalPntDesVector.push_back(" center ");

  evalPntVector.push_back(point3D(1.5*a, b/2, 0.));
  evalPntDesVector.push_back(" 1.5 side length away ");

  evalPntVector.push_back(point3D(a, b/2, b/2));
  evalPntDesVector.push_back(" above edge ");
}





