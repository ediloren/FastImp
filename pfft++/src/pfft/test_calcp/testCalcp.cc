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

void compWaveNumber (
		     const double f,
		     vector<double>& K0,
		     vector<complex<double> >& KC,
		     double& wavelen);

void compareTwoCalcp (
		      double wavelen);

void testCalcpForEikrOverR (
			    const vector<double>& K0,
			    const vector<complex<double> >& KC,
			    const element& srcPanel,
			    const vector<point3D>& evalPntVector,
			    const vector<string>& evalPntDesVector,
			    const double wavelen,
			    const double freq);

void showSrcPanelEvalPnt (
			  const element& srcPanel,
			  const vector<point3D>& evalPntVector,
			  const vector<string>& evalPntDesVector);

template <class T>
void showResult (
		 const vector<T>& slp,
		 const vector<T>& dlp,
		 const vector<vector3D<T> >& grad);
/*
void setSrcPanelEvalPnt (
			 element& srcPanel,
			 vector<point3D>& evalPntVector,
			 vector<string>& evalPntDesVector);
*/

void setSrcPanelEvalPntRectangle (
				  element& srcPanel,
				  vector<point3D>& evalPntVector,
				  vector<string>& evalPntDesVector,
				  double wavelen);			  

void setSrcPanelEvalPntTriangle (
				 element& srcPanel,
				 vector<point3D>& evalPntVector,
				 vector<string>& evalPntDesVector,
				 double wavelen,
				 double alpha);

void compute (
	      CalcpForOneOverR<element>& calcp,
	      const element* srcPanelPtr,
	      const vector<point3D>& evalPntVector,
	      vector<double>& slp,
	      vector<double>& dlp,
	      vector<vector3D<double> >& gradSlp);

template <class KT>
void compute (
	      CalcpForEikrOverR<KT, element>& calcp,
	      element srcPanel,
	      vector<point3D> evalPntVector,
	      vector<complex<double> >& slp,
	      vector<complex<double> >& dlp,
	      vector<vector3D<complex<double> > >& grad);

void compareResult (
		    const vector<double>& slp1,
		    const vector<double>& dlp1,
		    const vector<vector3D<double> >& grad1,
		    const vector<complex<double> >& slp2,
		    const vector<complex<double> >& dlp2,
		    const vector<vector3D<complex<double> > >& grad2,
		    const element& srcPanel,
		    const vector<point3D>& evalPntVector,
		    const vector<string>& evalPntDesVector);
/*


/**********************************************************************
 * main --
 **********************************************************************/
int
main (
      void)
{
  double freq;
  double alpha;
  double wavelen;

  vector<double> K0(2);
  vector<complex<double> > KC(2);
  
  element srcPanel;
  vector<point3D> evalPntVector;
  vector<string> evalPntDesVector;
  
  freq = 1.0e10;
  compWaveNumber (freq, K0, KC, wavelen);

  cout <<" #######################  Triangular Panel #################### "<<endl<<endl;

  // triangular panel: coarse mesh 
  alpha = 1.;
  cout << " <------------------ coarse mesh -------------------> ";
  cout << " all edges are around " << alpha <<" x wavelen " << endl << endl;
  setSrcPanelEvalPntTriangle (srcPanel,
			      evalPntVector, evalPntDesVector,
			      wavelen, alpha);

  testCalcpForEikrOverR(K0, KC, srcPanel,
			evalPntVector, evalPntDesVector,
			wavelen, freq);
  
  // triangular panel: finer mesh
  evalPntVector.clear();
  evalPntDesVector.clear();

  alpha = 0.5;
  cout << " <------------------ fine mesh -------------------> ";
  cout << " all edges are around " << alpha <<" x wavelen " << endl << endl;
  setSrcPanelEvalPntTriangle(srcPanel,
			     evalPntVector, evalPntDesVector,
			     wavelen, alpha);

  testCalcpForEikrOverR(K0, KC, srcPanel,
			evalPntVector, evalPntDesVector,
			wavelen, freq);

  // rectangular panel
  cout << "<-------------------- Rectangular Panel -------------------->" << endl<<endl;
  evalPntVector.clear();
  evalPntDesVector.clear();

  setSrcPanelEvalPntRectangle(srcPanel, evalPntVector, evalPntDesVector, wavelen);
  testCalcpForEikrOverR(K0, KC, srcPanel,
			evalPntVector, evalPntDesVector,
			wavelen, freq);

  compareTwoCalcp(wavelen);
  
  return 0;
}

/**********************************************************************
 * compareTwoCalcp --
 * -------------------------------------------------------------------
 * compare CalcpForOneOverR with CalcpForEikrOver @f=0 
 * those two calcp should have the same results.
 **********************************************************************/
void compareTwoCalcp (double wavelen)
{
  CalcpEvalPntDirection evalPointDirection = FOLLOW_NORMAL_DIRECTION;

  bool slpWanted = true;
  bool dlpWanted = true;
  bool gradOfSlpWanted = true;
  bool gradOfDlpWanted = false;

  typedef CalcpForOneOverR<element> CalcpForOneOverR;
  typedef CalcpForEikrOverR<double, element> CalcpForEikrOverR;

  CalcpForOneOverR calcpForOneOverR;
  CalcpForEikrOverR calcpForEikrOverR;

  double K0 = 0.;

  calcpForOneOverR.setupBasicParameters(slpWanted, dlpWanted,
					gradOfSlpWanted, gradOfDlpWanted);

  calcpForEikrOverR = CalcpForEikrOverR(K0, evalPointDirection,
					slpWanted, dlpWanted, gradOfSlpWanted);

  element srcPanel;
  vector<point3D> evalPntVector;
  vector<string> evalPntDesVector;

  setSrcPanelEvalPntRectangle (srcPanel, evalPntVector, evalPntDesVector, wavelen);

  calcpForOneOverR.findPanelVertexOrder(srcPanel);

  vector<double> slpOneOverR;
  vector<double> dlpOneOverR;
  vector<vector3D<double> > gradOneOverR;
  compute(calcpForOneOverR, &srcPanel, evalPntVector,
	  slpOneOverR, dlpOneOverR, gradOneOverR);

  vector<complex<double> > slpEikrOverR;
  vector<complex<double> > dlpEikrOverR;
  vector<vector3D<complex<double> > > gradEikrOverR;
  compute(calcpForEikrOverR, srcPanel, evalPntVector,
	  slpEikrOverR, dlpEikrOverR, gradEikrOverR);
  
  compareResult(slpOneOverR, dlpOneOverR, gradOneOverR,
		slpEikrOverR, dlpEikrOverR, gradEikrOverR,
		srcPanel, evalPntVector, evalPntDesVector);

}

/**********************************************************************
 * testCalcpForEikrOverR --
 **********************************************************************/
void testCalcpForEikrOverR (
			    const vector<double>& K0,
			    const vector<complex<double> >& KC,
			    const element& srcPanel,
			    const vector<point3D>& evalPntVector,
			    const vector<string>& evalPntDesVector,
			    const double wavelen,
			    const double freq)
{
  CalcpEvalPntDirection evalPointDirection = FOLLOW_NORMAL_DIRECTION;
  
  bool slpWanted = true;
  bool dlpWanted = true;
  bool gradWanted = true;
  
  CalcpForEikrOverR<complex<double>, element> calcp1;
  CalcpForEikrOverR<double, element> calcp2;

  vector<complex<double> > slp1, slp2;
  vector<complex<double> > dlp1, dlp2;
  vector<vector3D<complex<double> > > grad1, grad2;

  showSrcPanelEvalPnt (srcPanel, evalPntVector, evalPntDesVector);

  cout << endl <<" <---------------- freq = " << freq;
  cout << " ------------------> " << endl << endl;
  cout << endl << " <---------------- EMQS ----------------> " << endl;
  cout << " K0 := " << K0[0] << endl << endl;
  calcp2 = 
    CalcpForEikrOverR<double, element>(K0[0], evalPointDirection, 
				       slpWanted, dlpWanted, 
				       gradWanted);

  compute (calcp2, srcPanel, evalPntVector, slp2, dlp2, grad2);
  showResult(slp2, dlp2, grad2);

  cout << " KC := " << KC[0] << endl << endl;
  calcp1 = 
    CalcpForEikrOverR<complex<double>, element>(KC[0], evalPointDirection, 
						slpWanted, dlpWanted, 
						gradWanted);

  compute (calcp1, srcPanel, evalPntVector, slp1, dlp1, grad1);
  showResult(slp1, dlp1, grad1);

  cout << endl << " <----------------- Full Wave --------------> " << endl;
  cout << " K0 := " << K0[1] << endl << endl;
  calcp2 = 
    CalcpForEikrOverR<double, element>(K0[1], evalPointDirection, 
				       slpWanted, dlpWanted, 
				       gradWanted);

  compute (calcp2, srcPanel, evalPntVector, slp2, dlp2, grad2);
  showResult(slp2, dlp2, grad2);

  cout << " KC := " << KC[1] << endl <<endl;
  calcp1 = 
    CalcpForEikrOverR<complex<double>, element>(KC[1], evalPointDirection, 
						slpWanted, dlpWanted, 
						gradWanted);
  compute (calcp1, srcPanel, evalPntVector, slp1, dlp1, grad1);
  showResult(slp1, dlp1, grad1);

}

/**********************************************************************
 * showSrcPanelEvalPnt --
 **********************************************************************/
void showSrcPanelEvalPnt (
			  const element& srcPanel,
			  const vector<point3D>& evalPntVector,
			  const vector<string>& evalPntDesVector)
{
  cout << "           source panel is :            " <<endl;
  cout << " -------------------------------------- " <<endl<<endl;
  cout << srcPanel <<endl;

  for (size_t ii=0; ii<evalPntVector.size(); ii++) {
    cout << " eval pnt " << ii << ": " << evalPntDesVector[ii] << endl;
    cout << " ( " << evalPntVector[ii].x() 
	 << " , " << evalPntVector[ii].y()
	 << " , " << evalPntVector[ii].z() << " )" << endl << endl;
  } 
}

/**********************************************************************
 * showResult --
 **********************************************************************/
template <class T>
void showResult (
		 const vector<T>& slp,
		 const vector<T>& dlp,
		 const vector<vector3D<T> >& grad)
{
  for (size_t ii=0; ii<slp.size(); ii++) {
    cout << " pnt "<<ii<<"-> slp: "<<slp[ii]<<" dlp: "<<dlp[ii]<<endl;
    cout << "         grad ("<<grad[ii].x()<<" , "<<grad[ii].y()<<" , "<<grad[ii].z()<<endl<<endl;
  }				
}

/**********************************************************************
 * setSrcPanelEvalPnt --
 * Des - Description
 **********************************************************************/
/*
void setSrcPanelEvalPnt (
			 element& srcPanel,
			 vector<point3D>& evalPntVector,
			 vector<string>& evalPntDesVector)
{
  double a = 0.5e-3;
  double b = 2*a;
  double h = a/2.;
  double threshold = 1.e-5;

  point3D v1(0., 0., 0.);
  point3D v2(a, 0., 0.);
  point3D v3(a, b, 0.);
  point3D v4(0., b ,0.);  
  srcPanel = element(v1, v2, v3, v4);

  evalPntVector.reserve(3);
  evalPntDesVector.reserve(3);

  evalPntVector.push_back(point3D(a/2, b/2, 0.)); 
  evalPntDesVector.push_back(" center ");

  evalPntVector.push_back(point3D(1.5*a, b/2, 0.));
  evalPntDesVector.push_back(" 1.5 side length away ");

  evalPntVector.push_back(point3D(a, b/2, 1.5*a));
  evalPntDesVector.push_back(" above edge ");
}
*/

/**********************************************************************
 * setSrcPanelEvalPntTriangle --
 **********************************************************************/
void setSrcPanelEvalPntTriangle (
				 element& srcPanel,
				 vector<point3D>& evalPntVector,
				 vector<string>& evalPntDesVector,
				 const double wavelen,
				 const double alpha)
{
  const double FAR_FIELD_FACTOR = 200.0;
  double a = alpha * wavelen / 2.0;
  //double h = .1 * wavelen;
  double h = 2.0 * wavelen;

  point3D v1(0., 0., 0.);
  point3D v2(2*a, 0., 0.);
  point3D v3(a, 2*a, 0.);

  point3D centroid = (v1 + v2 + v3) / 3.0;
  point3D centerOfEdge01 = (v1 + v2) / 2.0;

  point3D unitZ(0., 0., 1.);
  point3D direction = centroid - v1;
  direction.normalize();

  srcPanel = element(v1, v2, v3);
  
  evalPntVector.reserve(9);
  evalPntDesVector.reserve(9);
  
  evalPntVector.push_back(centroid);
  evalPntDesVector.push_back(" right on the center ");

  evalPntVector.push_back(centroid + h * unitZ);
  evalPntDesVector.push_back(" right above center ");

  evalPntVector.push_back(centerOfEdge01 + h * unitZ);
  evalPntDesVector.push_back(" right above the edge ");

  evalPntVector.push_back(v3 + h * unitZ);
  evalPntDesVector.push_back(" right above the vertex ");

  evalPntVector.push_back(v1 + 0.05 * wavelen * direction);
  evalPntDesVector.push_back(" 1/20 wavelen away from centroid along 1c ");

  evalPntVector.push_back(v1 + 0.1 * wavelen * direction); 
  evalPntDesVector.push_back(" 1/10 wavelen away from centroid along 1c ");

  evalPntVector.push_back(v1 + 0.25 * wavelen * direction);
  evalPntDesVector.push_back(" 1/4 wavelen away from centroid along 1c ");

  evalPntVector.push_back(v1 + .5 * wavelen * direction);
  evalPntDesVector.push_back(" 1/2 wavelen away from centroid along 1c ");
  
  evalPntVector.push_back(centroid +
			  FAR_FIELD_FACTOR * vector3D<double>(1.0, 0.5, 0.7));
  evalPntDesVector.push_back(" far field approximation should be trigged ");

}


/**********************************************************************
 * setSrcPanelEvalPntRectangle --
 * ...................................................................
 * setup src panel and eval pnts for the test of calcpForEikrOverR. 
 * the correct results are:
 * for f = 1.0e10, 
 * EMQS case: K0 = 0. KC = 
 * FULLWAVE case: K0 = 
 *   
 * 
 **********************************************************************/
void setSrcPanelEvalPntRectangle (
				  element& srcPanel,
				  vector<point3D>& evalPntVector,
				  vector<string>& evalPntDesVector,
				  const double wavelen)
{
  const double FAR_FIELD_FACTOR = 200.0;

  double a = 0.5e-3;
  double b = 2*a;
  double h = a/2.0;
  double threshold = 1.e-5;

  point3D v1(0., 0., 0.);
  point3D v2(0., b, 0.);
  point3D v3(a, b, 0.);
  point3D v4(a, 0.,0.);  
  
  srcPanel = element(v1, v2, v3, v4);
  
  evalPntVector.reserve(9);
  evalPntDesVector.reserve(9);
  
  evalPntVector.push_back(point3D(a/2, b/2, 0.));
  evalPntDesVector.push_back(" right on the center ");

  evalPntVector.push_back(point3D(a/2, b/2, h));
  evalPntDesVector.push_back(" right above center ");

  evalPntVector.push_back(point3D(a, b/2, h));
  evalPntDesVector.push_back(" above the edge ");

  evalPntVector.push_back(point3D(a, b, h));
  evalPntDesVector.push_back(" right above one vertex ");

  evalPntVector.push_back(point3D(.05*wavelen, b/2, 0));
  evalPntDesVector.push_back(" 1/20 wave length away ");

  evalPntVector.push_back(point3D(.1*wavelen, b/2, 0));
  evalPntDesVector.push_back(" 1/10 wave length away");

  evalPntVector.push_back(point3D(.25*wavelen, b/2, 0.));
  evalPntDesVector.push_back(" 1/4 wave length away ");

  evalPntVector.push_back(point3D(.5*wavelen, b/2, 0.));
  evalPntDesVector.push_back(" 1/2 wave length away ");

  evalPntVector.push_back(point3D(a/2, b/2, 0.)
			  + FAR_FIELD_FACTOR * point3D(1.0, 0.5, 0.7));
  evalPntDesVector.push_back(" far field ");

}

/**********************************************************************
 * compute --
 **********************************************************************/
void compute (
	      CalcpForOneOverR<element>& calcp,
	      const element* srcPanelPtr,
	      const vector<point3D>& evalPntVector,
	      vector<double>& slp,
	      vector<double>& dlp,
	      vector<vector3D<double> >& gradSlp)
{
  size_t n = evalPntVector.size();
  
  vector<vector3D<double> > gradDlp;

  slp.resize(n);
  dlp.resize(n);
  gradSlp.resize(n);
  gradDlp.resize(n);

  for (size_t ii=0; ii<n; ii++) {
    calcp(srcPanelPtr, evalPntVector[ii],
	  slp[ii], dlp[ii], gradSlp[ii], gradDlp[ii]);
  }
}

/**********************************************************************
 * compute --
 **********************************************************************/
template <class KT>
void compute (
	      CalcpForEikrOverR<KT, element>& calcp,
	      element srcPanel,
	      vector<point3D> evalPntVector,
	      vector<complex<double> >& slp,
	      vector<complex<double> >& dlp,
	      vector<vector3D<complex<double> > >& gradSlp)
{
  size_t n = evalPntVector.size();
  
  slp.resize(n);
  dlp.resize(n);
  gradSlp.resize(n);
  
  for (size_t ii=0; ii<n; ii++) {
    calcp(srcPanel, evalPntVector[ii],
	  slp[ii], dlp[ii], gradSlp[ii]);
  }
}

/**********************************************************************
 * compareResult --
 **********************************************************************/
//template <class T1, class T2>
void compareResult (
		    const vector<double>& slp1,
		    const vector<double>& dlp1,
		    const vector<vector3D<double> >& grad1,
		    const vector<complex<double> >& slp2,
		    const vector<complex<double> >& dlp2,
		    const vector<vector3D<complex<double> > >& grad2,
		    const element& srcPanel,
		    const vector<point3D>& evalPntVector,
		    const vector<string>& evalPntDesVector)
{
  cout << "                    source panel "<<endl;
  cout << " --------------------------------------------------- "<<endl;
  cout << "  " <<srcPanel<< "  " <<endl<<endl;
  
  size_t n = evalPntVector.size();
  for (size_t i=0; i<n; i++) {
    cout << " evalPnt 1: " << evalPntDesVector[i] << endl;
    cout << " " << evalPntVector[i] << endl;
    cout << " Slp 1: "<<slp1[i]<<"\tSlp 2: "<<slp2[i]<<endl;
    cout << " Dlp 1: "<<dlp1[i]<<"\tDlp 2: "<<dlp2[i]<<endl;
    cout << " Grad 1.x: "<<grad1[i].x()<<"\t Grad 2.x "<<grad2[i].x()<<endl;
    cout << " Grad 1.y: "<<grad1[i].y()<<"\t Grad 2.y "<<grad2[i].y()<<endl;
    cout << " Grad 1.z: "<<grad1[i].z()<<"\t Grad 2.z "<<grad2[i].z()<<endl;
  }
}
  
/**********************************************************************
 * compWaveNumber --
 **********************************************************************/
void
compWaveNumber (
		const double f,
		vector<double>& K0,
		vector<complex<double> >& KC,
		double& wavelen)
{
  // 0: EMQS; 1: FULL WAVE
  K0[0] = 0.;
  K0[1] = 2. * PAI * f / 3.e8;

  wavelen = 3.e8 / f;

  const double sigma = 5.8e7;
  const double mu0 = 4e-7 * PAI;
  
  for (size_t ii=0; ii<2; ii++) {
    KC[ii] = complex<double>(K0[ii]*K0[ii], - 2. * PAI * f * mu0 * sigma);
    KC[ii] = sqrt(KC[ii]);
    if (real(KC[ii]) > 0) {
      KC[ii] *= -1;
    }
  }
}


/**********************************************************************
 * computeMore --
 **********************************************************************/
/*
template <class Calcp>
void
computeMore (
	     Calcp calcp)
{
  complex<double> slp;
  complex<double> dlp;
  vector3D<complex<double> > grad;

  double a = 0.5e-3;
  double b = 2*a;
  double h = a/2.;
  point3D v1(0., 0., 0.);
  point3D v2(0., b, 0.);
  point3D v3(a, b, 0.);
  point3D v4(a, 0.,0.);  
  element panel(v1, v2, v3, v4);
  double threshold = 1e-5;

  point3D evalPnt3(a/2, b/2, h);
  calcp(panel, evalPnt3, slp, dlp, grad);
  cout << "above center==" << slp << "\t" << dlp << endl;

  point3D evalPnt4(a, b/2, 0.);
  calcp(panel, evalPnt4, slp, dlp, grad);
  cout << "on the edge==" << slp << "\t" << dlp << endl;
}
*/







