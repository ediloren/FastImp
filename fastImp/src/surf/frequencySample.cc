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

  const static char cvsid[] = "$Id: frequencySample.cc,v 1.3 2002/07/18 15:05:45 zhzhu Exp $";

  ==========================================================================
*/

#include "frequencySample.h"
#include "service.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace surf;
using namespace std;

/**********************************************************************
 * FrequencySample --
 **********************************************************************/
FrequencySample::FrequencySample (
				  const double firstPoint, 
				  const double lastPoint, 
				  const int numSam,
				  const SampleMethod samMethod)
  : firstPoint_(firstPoint), lastPoint_(lastPoint), 
  numSam_(numSam), samMethod_(samMethod)
{
  if (samMethod == LOG) {
    setupLogSamPoint();
  } else if (samMethod == LINEAR) {
    setupLinearSamPoint();
  } else {
    errorMessage("frequencySample.cc : frequencySample",
		 "Illegal frequency sampling method");
  }

  try {
    checkRange();
  } catch (const domain_error error) {
    cerr << error.what() << endl;
    errorMessage("frequencySample.cc : frequencySample",
		 "Illegal frequency sampling data");
  }

  if (numSam_ >= 3) {
    size_t midPointIndex = numSam_/2;
    middlePoint_ = samPointList_[midPointIndex];
  } else {
    middlePoint_ = firstPoint_;
  }
}

/**********************************************************************
 * setupLogSamPoint --
 **********************************************************************/
void
FrequencySample::setupLogSamPoint (
				   void)
{
  double step;
  if (numSam_ >= 2) {
    step = (log10(lastPoint_) - log10(firstPoint_))/(numSam_-1);
    step = pow(10, step);
  }

  double f = firstPoint_;
  samPointList_.reserve(numSam_);
  samPointList_.push_back(f);
  for (int i = 1; i < numSam_; i++) {
    f *= step;
    samPointList_.push_back(f);
  }
}

/**********************************************************************
 * setupLinearSamPoint --
 **********************************************************************/
void
FrequencySample::setupLinearSamPoint (
				      void)
{
  samPointList_.resize(numSam_);
  if (numSam_ >= 2) {
    double step = (lastPoint_ - firstPoint_) / (numSam_-1);
    for (size_t i = 0; i < numSam_; i++) {
      samPointList_[i] = firstPoint_ + i * step;
    }
  } else {
    samPointList_[0] = firstPoint_;
  }
}

/**********************************************************************
 * checkRange --
 **********************************************************************/
void
FrequencySample::checkRange (
			     void)
{
  if (firstPoint_ < 0) {
    throw domain_error("Begin of the frequency range should be positive!");
  }

  if (lastPoint_ < 0) {
    throw domain_error("End of the frequency range should be positive!");
  }

  if (lastPoint_ < firstPoint_) {
    throw domain_error("Begin of the frequency range should be less than end of the frequency range!");
  }

  if ( ((firstPoint_ != lastPoint_) && (numSam_ == 1) )   ||
       ((firstPoint_ == lastPoint_) && (numSam_ != 1) ) ) {
    throw domain_error("Inconsistant frequency range and sampling number!");
  } 

  if (numSam_ <= 0) {
    throw domain_error("Number of frequency sampling points should be positive!");
  }

  if ( (numSam_ > 1) && (lastPoint_ < firstPoint_) ) {
    throw domain_error("End of the frequency range should be larger than the begin the range!");
  }
}

