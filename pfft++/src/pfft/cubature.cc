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

  ==========================================================================
*/

const static char cvsid[] = "$Id: cubature.cc,v 1.5 2002/10/19 01:03:51 bsong Exp $";

#include "cubature.h"
#include "utils.h"
#include <stdexcept>
#include <assert.h>

using namespace pfft;
using namespace std;

/**********************************************************************
 * cubature --
 **********************************************************************/
Cubature::Cubature (
		    const int degree, 
		    const int rule, 
		    const Shape shape)
  : degree_(degree), rule_(rule), shape_(shape), currentPointIndex_(0)
{
  if (shape_ == TRIANGLE) {
    setupTriangleCubtature();
  } else if (shape_ == SQUARE) {
    setupSquareCubtature();
  } else {
    throw domain_error("unknown cubature shape");
  }
}

/**********************************************************************
 * point --
* Barycentric coordinates (x,y,z) are used to represent points in the 
* plane of the triangle. If the corners of the triangle are denoted 
* by (V1,V2,V3), then the Cartesian coordinates of a point are 
* x*V1 + y*V2 + z*V3.
 **********************************************************************/
point3D 
Cubature::point (
		 const vector<point3D>& v, 
		 const size_t pointIndex) const 
{
  if (shape_ == TRIANGLE) {
    return point_[pointIndex].x() * v[0] + 
           point_[pointIndex].y() * v[1] +
           point_[pointIndex].z() * v[2];
  } else if (shape_ == SQUARE) {
    return point3D(0., 0., 0.);
  } else {
    throw domain_error("unknown cubature shape");
  }
}

/**********************************************************************
 * weight --
 * The triangle cubature rules are set up for a triangle with area 1/2. 
 * Thus, the Jacobian of the transfomation from this "unit" 
 * triangle to the present element is equal to TWICE the 
 * area of the element.                                              
 **********************************************************************/
double 
Cubature::weight (
		  const double area, 
		  const size_t pointIndex) const 
{
  if (shape_ == SQUARE)
    return area*2. * weight_[pointIndex];
  else 
    return area * weight_[pointIndex];
}

/**********************************************************************
 * setupTriangleCubtature --
  Cubature formulae for integration over a triangle with corners
  (0,0) (0,1) and (0,1). This is a two-dimensional simplex.
  Returns points (x,y) and weights for the chosen scheme.
  The schemes have been found from the on-line 
  "Encyclopaedia of Cubature Formulas" at
  URL: http://www.cs.kuleuven.ac.be/~nines/research/ecf/ecf.html
  Rules have been copied using "cut and paste" to reduce the risk of 
  typos. Any rule of order N integrates exactly all polynomials of 
  degree less than or equal to N.
  
  degree: Degree of the cubature rule to use.
  rule:  Number of rule of this degree to use (normally use irule=1).
  
  (a,b,c) Barycentric coordinates of cubture points (see below).
  w: Weights of the cubature rule. w[i] is the weight to use 
  at point (x[i],y[i]), i.e. the cubature rule is:

  INT \approx SUM_{i=0}^{order-1} w[i]*f(x[i],y[i]) 

  Since the area of the domain of integration is one half, the weights 
  of each scheme should add up to one half.

  Note that a+b+c=1 (the Barycentric coordinates contain redundant 
  information). This is exploited in the present routine, so that 
  only two of the three coordinates are given at any time. The last 
  coordinate is found by subtracting from unity the values of the other 
  two coordinates. 
 **********************************************************************/
void
Cubature::setupTriangleCubtature (
				void)
{
  findNumPointInTriangleCubature();
  point_.resize(numPoint_);
  weight_.resize(numPoint_);
  findPointAndWeightInTriangle();
}

/**********************************************************************
 * findNumPointInTriangleCubature --
 **********************************************************************/
void
Cubature::findNumPointInTriangleCubature (
					  void)
{
  switch (100*degree_ + rule_) {
  case  101: numPoint_= 1; break;
  case  201: numPoint_= 3; break;
  case  301: numPoint_= 4; break; /* Contains negative weights */
  case  302: numPoint_= 6; break;
  case  401: numPoint_= 6; break;
  case  501: numPoint_= 7; break;
  case  601: numPoint_=12; break;
  case  602: numPoint_=12; break;
  case  701: numPoint_=12; break;
  case  801: numPoint_=16; break;
  case  901: numPoint_=19; break;
  case 1001: numPoint_=25; break;
  case 1002: numPoint_=25; break;
  case 1101: numPoint_=28; break;
  default:
    warningMessage("cubature.cc : findNumPointInTriangleCubature",
		   "Unknown rule or degree, use the default highest degree cubature instead");
    numPoint_=28; break;
  }
}

/**********************************************************************
 * findPointAndWeightInTriangle --
 **********************************************************************/
void
Cubature::findPointAndWeightInTriangle (
					void)
{
  switch (100*degree_ + rule_) {
  case 101:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 1 
     * Points: 1 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 0 0 
     * REF: [Str71] */
    /*npoints=1;*/
    /* Generator: [ Fully symmetric ] (= single node in this case) */
    SIMPLEX2_FULLY1(0.33333333333333333,0.33333333333333333, 0.50000000000000000);

    currentPointIndex_ = 0;
    break;

  case 201:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 2 
     * Points: 3 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 0 1 0 
     * REF: [Str71] */
    /*npoints=3;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.16666666666666666,0.16666666666666666, 0.16666666666666666);
    currentPointIndex_ = 0;
    break;

  case 301:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 3 
     * Points: 4 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 1 0 
     * REF: [Str71] 
     *      NOTE: This rule contains negative weights! */
    /*npoints=4;*/
    /* Generator: [ Fully symmetric ] (= single node in this case) */
    SIMPLEX2_FULLY1(0.33333333333333333,0.33333333333333333, -0.28125000000000000);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.20000000000000000,0.20000000000000000,  0.26041666666666667);
    currentPointIndex_ = 0;
    break;

  case 302:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 3 
     * Points: 6 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 0 0 1 
     * REF: [Str71] */
    /*npoints=6;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.10903900907287721,0.23193336855303057, 0.083333333333333333);
    break;

  case 401:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 4 
     * Points: 6 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 0 2 0 
     * REF: [cow73][dun85b][lg78][lj75][moa74][sf73] */
    /*npoints=6;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.091576213509770743,0.091576213509770743, 0.054975871827660933);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.44594849091596488, 0.44594849091596488,  0.11169079483900573);

    currentPointIndex_ = 0;
    break;

  case 501:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 5 
     * Points: 7 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 2 0 
     * REF: [str71] */
    /*npoints=7;*/
    /* Generator: [ Fully symmetric ]  (= single node in this case) */
    SIMPLEX2_FULLY1(0.33333333333333333,0.33333333333333333, 0.11250000000000000);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.10128650732345633,0.10128650732345633, 0.062969590272413576);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.47014206410511508,0.47014206410511508, 0.066197076394253090);

    currentPointIndex_ = 0;
    break;

  case 601:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 6 
     * Points: 12 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 0 2 1 
     * REF: [cow73][dun85b][lg78][lj75][moa74][sf73] */
    /*npoints=12;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.063089014491502228,0.063089014491502228, 0.025422453185103408);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.24928674517091042, 0.24928674517091042,  0.058393137863189683);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.053145049844816947,0.31035245103378440,  0.041425537809186787);

    currentPointIndex_ = 0;
    break;

  case 602:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 6 
     * Points: 12 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 0 2 1 
     * REF: [bec87] */
    /*npoints=12;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.21942998254978296,0.21942998254978296, 0.085666562076490515);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.48013796411221504,0.48013796411221504, 0.040365544796515489);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.83900925971479105,0.14161901592396815, 0.020317279896830331);

    currentPointIndex_ = 0;
    break;

  case 701:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 7 
     * Points: 12 
     * Structure: Ro3 invariant 
     * Rule struct: 0 0 0 0 0 4 
     * REF: [gat88] */
    /*npoints=12;*/
    /* Generator: [ Ro3 invariant ] (Cyclic rotations) */
    SIMPLEX2_ROTATIONAL3(0.062382265094402118,0.067517867073916085, 0.026517028157436251);
    /* Generator: [ Ro3 invariant ] */
    SIMPLEX2_ROTATIONAL3(0.055225456656926611,0.32150249385198182,  0.043881408714446055);
    /* Generator: [ Ro3 invariant ] */
    SIMPLEX2_ROTATIONAL3(0.034324302945097146,0.66094919618673565,  0.028775042784981585);
    /* Generator: [ Ro3 invariant ] */
    SIMPLEX2_ROTATIONAL3(0.51584233435359177,0.27771616697639178, 0.067493187009802774);

    currentPointIndex_ = 0;
    break;

  case 801:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 8 
     * Points: 16 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 3 1 
     * REF: [lj75][dun85b][lg78] */
    /*npoints=16;*/
    /* Generator: [ Fully symmetric ] (= single node in this case) */
    SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.072157803838893584);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.17056930775176020, 0.17056930775176020,  0.051608685267359125);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.050547228317030975,0.050547228317030975, 0.016229248811599040);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.45929258829272315, 0.45929258829272315,  0.047545817133642312);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.72849239295540428, 0.26311282963463811,  0.013615157087217497);

    currentPointIndex_ = 0;
    break;

  case 901:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 9 
     * Points: 19 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 4 1 
     * REF: [lj75][dun85b][lg78] */
    /*npoints=19;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.048567898141399416);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.48968251919873762, 0.48968251919873762,  0.015667350113569535);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.43708959149293663, 0.43708959149293663,  0.038913770502387139);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.18820353561903273, 0.18820353561903273,  0.039823869463605126);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.044729513394452709,0.044729513394452709, 0.012788837829349015);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.74119859878449802, 0.036838412054736283, 0.021641769688644688);

    currentPointIndex_ = 0;
    break;

  case 1001:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 10 
     * Points: 25 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 2 3 
     * REF [lg78] */
    /*npoints=25;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.039947252370619853);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.42508621060209057, 0.42508621060209057,  0.035561901116188667);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.023308867510000190,0.023308867510000190, 4.1119093452320977e-3);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.62830740021349255, 0.22376697357697300,  0.022715296148085009);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.61131382618139764, 0.35874014186443146,  0.018679928117152638);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.82107206998562937, 0.14329537042686714,  0.015443328442281994);

    currentPointIndex_ = 0;
    break;

  case 1002:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 10 
     * Points: 25 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 0 1 2 3 
     * REF: [lg78] */
    /*npoints=25;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.040871664573142983);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.14216110105656438, 0.14216110105656438,  0.022978981802372364);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.032055373216943512,0.032055373216943512, 6.6764844065747831e-3);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.53005411892734402, 0.32181299528883542,  0.031952453198212022);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.60123332868345924, 0.36914678182781098,  0.017092324081479714);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.80793060092287906, 0.16370173373718249,  0.012648878853644192);

    currentPointIndex_ = 0;
    break;

  case 1101:
    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 11 
     * Points: 28 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 1 1 5 1 
     * REF: [lj75] */
    /*npoints=28;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.85887028128263670, 0.14112971871736329,  3.6811918916502771e-3);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.043988650581116119);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.025989140928287395,0.025989140928287395, 4.3721557768680115e-3);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.094287502647922495,0.094287502647922495, 0.019040785996967468);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.49463677501721381, 0.49463677501721381,  9.4277240280656460e-3); 
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.20734338261451133, 0.20734338261451133,  0.036079848772369763);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.43890780570049209, 0.43890780570049209,  0.034664569352767949);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.67793765488259040, 0.044841677589130443, 0.020528157714644283);

    currentPointIndex_ = 0;
    break;

  default:
    warningMessage("cubature.cc : findPointAndWeightInTriangle",
		   "Unknown rule or degree, use the default highest degree cubature instead");

    /* Region: Simplex 
     * Dimension: 2 
     * Degree: 11 
     * Points: 28 
     * Structure: Fully symmetric 
     * Rule struct: 0 0 1 1 5 1 
     * REF: [lj75] */
    /*npoints=28;*/
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.85887028128263670, 0.14112971871736329,  3.6811918916502771e-3);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY1(0.33333333333333333, 0.33333333333333333,  0.043988650581116119);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.025989140928287395,0.025989140928287395, 4.3721557768680115e-3);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.094287502647922495,0.094287502647922495, 0.019040785996967468);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.49463677501721381, 0.49463677501721381,  9.4277240280656460e-3); 
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.20734338261451133, 0.20734338261451133,  0.036079848772369763);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY3(0.43890780570049209, 0.43890780570049209,  0.034664569352767949);
    /* Generator: [ Fully symmetric ] */
    SIMPLEX2_FULLY6(0.67793765488259040, 0.044841677589130443, 0.020528157714644283);

    currentPointIndex_ = 0;
    break;
    //    throw domain_error("Unknown rule or degree");
  } 
}

/**********************************************************************
 * SIMPLEX2_FULLY1 --
 * Fully symmetric, 1-point resultant: must be (1/3,1/3,1/3) 
 **********************************************************************/
void
Cubature::SIMPLEX2_FULLY1 (
			   double coord1,
			   double coord2,
			   double weight)
{
  assert(coord1==0.33333333333333333); 
  assert(coord1==coord2); 
  double coord3 = 1.0 - coord1 - coord2;
  point_[currentPointIndex_] = vector3D<double>(coord1, coord2, coord3);
  weight_[currentPointIndex_++] = weight;
}


/**********************************************************************
 * SIMPLEX2_FULLY3 --
 * Fully symmetric, 3-point resultant: can be written so that coord1=coord2. 
 **********************************************************************/
void
Cubature::SIMPLEX2_FULLY3 (
			   double coord1, 
			   double coord2,
			   double weight)
{
  assert(coord1==coord2); 
  double coord3 = 1.0 - coord1 - coord2;
  point_[currentPointIndex_] = vector3D<double>(coord1, coord2, coord3); 
  weight_[currentPointIndex_ ++] = weight;

  point_[currentPointIndex_] = vector3D<double>(coord3, coord1, coord2);
  weight_[currentPointIndex_ ++] = weight;

  point_[currentPointIndex_] = vector3D<double>(coord2, coord3, coord1);
  weight_[currentPointIndex_ ++] = weight;

}

/**********************************************************************
 * SIMPLEX2_FULLY3 --
 * Fully symmetric, 6-point resultant 
 **********************************************************************/
void
Cubature::SIMPLEX2_FULLY6 (
			   double coord1,
			   double coord2,
			   double weight)
{ 
  assert(coord1!=coord2); 
  double coord3 = 1.0 - coord1 - coord2;
  point_[currentPointIndex_] = vector3D<double>(coord1, coord2, coord3);
  weight_[currentPointIndex_ ++] = weight;

  point_[currentPointIndex_] = vector3D<double>(coord1, coord3, coord2); 
  weight_[currentPointIndex_ ++] = weight;

  point_[currentPointIndex_] = vector3D<double>(coord2, coord1, coord3); 
  weight_[currentPointIndex_ ++] = weight;

  point_[currentPointIndex_] = vector3D<double>(coord3, coord1, coord2); 
  weight_[currentPointIndex_ ++] = weight;

  point_[currentPointIndex_] = vector3D<double>(coord2, coord3, coord1); 
  weight_[currentPointIndex_ ++] = weight;

  point_[currentPointIndex_] = vector3D<double>(coord3, coord2, coord1); 
  weight_[currentPointIndex_ ++] = weight;

}

/**********************************************************************
 * SIMPLEX2_FULLY3 --
 * Rotational symmetric (Ro3 invariant) 3-point resultant
 **********************************************************************/
void
Cubature::SIMPLEX2_ROTATIONAL3 (
				double coord1,
				double coord2,
				double weight)
{ 
  double coord3 = 1.0 - coord1 - coord2;
  point_[currentPointIndex_] = vector3D<double>(coord1, coord2, coord3); 
  weight_[currentPointIndex_ ++] = weight;

  point_[currentPointIndex_] = vector3D<double>(coord3, coord1, coord2); 
  weight_[currentPointIndex_ ++] = weight;

  point_[currentPointIndex_] = vector3D<double>(coord2, coord3, coord1); 
  weight_[currentPointIndex_ ++] = weight;
}

/**********************************************************************
 * setupSquareCubtature --
 **********************************************************************/
void
Cubature::setupSquareCubtature (
				void)
{

}











