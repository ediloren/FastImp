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

  const static char cvsid[] = "$Id: frequencySample.h,v 1.1.1.1 2002/04/15 21:41:01 bsong Exp $";

  ==========================================================================
*/

#ifndef _FREQUENCY_SAMPLE_H_
#define _FREQUENCY_SAMPLE_H_

#include <vector>

namespace surf {

  enum SampleMethod { LOG, LINEAR, NUM_SAMPLE_METHOD };

  class FrequencySample {

  public:
    FrequencySample(const double firstPoint = 1e6, 
		    const double lastPoint = 1e6, 
		    const int numSam = 1,
		    const SampleMethod samMethod = LOG);
    double point(int i) const { return samPointList_[i]; }
    int numSam(void) const { return numSam_; }
    double firstPoint(void) const { return firstPoint_; }
    double lastPoint(void) const { return lastPoint_; }
    SampleMethod samMethod(void) const { return samMethod_; }
    double middlePoint(void) const { return middlePoint_; }

  private:
    double firstPoint_;
    double lastPoint_;
    double middlePoint_;
    int numSam_;
    SampleMethod samMethod_;
    std::vector<double> samPointList_;

    void setupLogSamPoint (void);
    void setupLinearSamPoint (void);
    void checkRange(void);

  };

} // namespace surf

#endif
