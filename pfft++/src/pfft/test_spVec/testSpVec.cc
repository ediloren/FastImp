#include "spVec.h"
#include <iostream>
#include <complex>
#include <valarray>

using namespace pfft;

/**********************************************************************
 * main --
 **********************************************************************/
int
main (
      void)
{
  SpVec<double> sv;
  sv.reserve(10);
  sv.push_back(SpVecElement<double>(1, 2.));
  sv.push_back(SpVecElement<double>(2, 3.));
  sv.push_back(SpVecElement<double>(3, 4.));
  sv.push_back(SpVecElement<double>(4, 5.));
  cout << sv;

  cout << sv.index(0) << sv.value(0) << endl;
  cout << sv[0].index() << sv[0].value() << endl;

  SpVec<double> sv1 = sv;
  cout << sv1;

  sv1 *= 2;
  cout << sv1;

  SpVec<double> sv2(sv);
  cout << sv2;

  double ans;
  dotProd(ans, sv2, sv1);
  cout << ans << endl;

  sv2[0] += 10.;
  cout << sv2[0] << std::endl;

  sv2.clear();
  cout << sv2;

}

