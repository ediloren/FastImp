#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "sphereGen.h"

int main(void) 
{
   const size_t n = 2;
   vector<const char*> filenameList;
   vector<double> radiusList;
   vector<size_t> ndenseList;

/* for the generation of static testing cases. 
   the size of the sphere is fixed at 1. 
*/
   ndenseList.push_back(2);
   ndenseList.push_back(3);
   ndenseList.push_back(7);
   ndenseList.push_back(9);
   ndenseList.push_back(20);
   ndenseList.push_back(30);
   ndenseList.push_back(65);

   radiusList.push_back(1.0);
   radiusList.push_back(1.0);
   radiusList.push_back(1.0);
   radiusList.push_back(1.0);
   radiusList.push_back(1.0);
   radiusList.push_back(1.0);
   radiusList.push_back(1.0);
   
   filenameList.push_back("sphere48.qui");
   filenameList.push_back("sphere108.qui");
   filenameList.push_back("sphere588.qui");
   filenameList.push_back("sphere972.qui");
   filenameList.push_back("sphere4800.qui");
   filenameList.push_back("sphere10800.qui");
   filenameList.push_back("sphere50700.qui");

/*   for the generation of high frequency testing elements.
     the maximum size of the elements is only a fraction 
     (1/8 or 1/10) of the wavelength. 
   ndenseList.push_back(2);
   ndenseList.push_back(3);
   ndenseList.push_back(7);
   ndenseList.push_back(9);
   ndenseList.push_back(20);
   ndenseList.push_back(30);

   radiusList.push_back(0.0084);
   radiusList.push_back(0.0093);
   radiusList.push_back(0.0194);
   radiusList.push_back(0.0248);
   radiusList.push_back(0.0530);
   radiusList.push_back(0.0795);
   
   filenameList.push_back("hfSphere48.qui");
   filenameList.push_back("hfSphere108.qui");
   filenameList.push_back("hfSphere588.qui");
   filenameList.push_back("hfSphere972.qui");
   filenameList.push_back("hfSphere4800.qui");
   filenameList.push_back("hfSphere10800.qui");
 */

   for (size_t ii = 0; ii < radiusList.size(); ii ++) {
	   sphereGen(filenameList[ii], radiusList[ii], ndenseList[ii]);
   }
	 
   return 0;
}

