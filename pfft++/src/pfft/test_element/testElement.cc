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

#include "element.h"
#include "vector3D.h"
#include <fstream>

using namespace pfft;
void
testIO_1 (
	  void);

/**********************************************************************
 * main --
 **********************************************************************/
int
main (
      void)
{
  //  char *fileName = "woven01n01.qui";
  char *fileName = "sphere108.qui";
  vector<element> elementList;
  std::string title;
  size_t numOfConductor;

  //    testIO_1();
  readElementFile(fileName, elementList, title, numOfConductor);
  cout << title << endl << endl;
  cout << " numOfConductor "<< numOfConductor << endl;
  for (size_t ii=0; ii<numOfConductor; ii++) {
    cout << " <---------------  conductor " << ii << "  ----------------->" << endl;
    cout << elementList[ii] << endl;
  }

  return 1;
}

/**********************************************************************
 * testIO_1 --
 **********************************************************************/
void
testIO_1 (
	  void)
{
  vector3D<double> v1(0., 0., 0.);
  vector3D<double> v2(1., 0., 0.);
  vector3D<double> v3(0., 1., 0.);

  element panel1(v1, v2, v3);
  cout << panel1;
  
  element panel2(panel1);
  cout << panel2;
}


