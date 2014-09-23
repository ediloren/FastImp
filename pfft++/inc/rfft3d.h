/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Zhenhai Zhu, Ben Song
  
  Description:

  Resources:

  See also:

  const static char cvsid[] = "$Id: rfft3d.h,v 1.3 2002/08/01 14:02:38 zhzhu Exp $";

  ==========================================================================
*/
#ifndef _RFFT3D_H_
#define _RFFT3D_H_

#include <iostream>
#include <vector>
#include "fftwInterface.h"
#include "x3dconv.h"

namespace pfft{

  /**********************************************************************
   * rfft3DGeneral --
   **********************************************************************/
  template<class KernelContainer>
  void rfft3DGeneral(const size_t *dim,
		     KernelContainer&  a,
		     fftw_direction dir,
		     bool SCALE_FLAG)
  {
    const size_t n1 = dim[0];
    const size_t n2 = dim[1];
    const size_t n3 = dim[2];
    const size_t dimHalf[] = {n1, n2, n3>>1};
    
    int isign;
    isign = (dir == FFTW_BACKWARD) ? +1 : -1;
    
    const double c1 = 0.5;
    const double c2 = -0.5 * isign;
    const double theta = isign * 2.0 * PAI / static_cast<double>(n3);
    
    double wr, wi;
    double wtmp = sin(0.5 * theta);

    const double wpr = -2.0 * wtmp * wtmp;
    const double wpi = sin(theta);
    
    size_t i1, i2, i3, j1, j2, j3, ii3;
    double h1r, h1i, h2r, h2i;

    if(isign == 1){
      // fft3d view real a as packed complex array!
      fftw_complex *in = reinterpret_cast<fftw_complex*>(a.p_);
     
      ::fftwnd_plan plan 
	  = ::fftw3d_create_plan(dimHalf[0], dimHalf[1], dimHalf[2],
				 FFTW_BACKWARD,
				 FFTW_ESTIMATE|FFTW_IN_PLACE);

      // this function always causes segmentation fault on ia64-linux
      // might be a bug for fftw.
      /*
      fftwnd_plan plan
	= fftwnd_create_plan(3, 
			     reinterpret_cast<int*>(const_cast<size_t*>(dimHalf)),
			     dir,
			     FFTW_ESTIMATE|FFTW_IN_PLACE);
      */

      ::fftwnd_one(plan, in, static_cast<fftw_complex*>(0));
      ::fftwnd_destroy_plan(plan);
      for(i1=0; i1<n1; ++i1){
	for(i2=0, j2=0; i2<n2; ++i2){
	  a.speq_[(j2++) + (2*n2)*i1] = a.p_[0 + n3*(i2+n2*i1)]; 
	  a.speq_[(j2++) + (2*n2)*i1] = a.p_[1 + n3*(i2+n2*i1)];
	}
      }
    }    

    for(i1=0; i1<n1; ++i1){
      j1 = (i1 != 0? n1 - i1 : 0);
      //Zero frequency is its own reflection, otherwise locate corresponding
      //negative frequency in wrap-around order.
      wr = 1.0;
      wi = 0.0;
      for(ii3=0, i3=0; i3<(n3>>2)+1; i3++, ii3+=2){
	for(i2=0; i2<n2; ++i2){
	  if(i3 == 0){
	    j2 = (i2 != 0 ? ((n2-i2)<<1) : 0);
	    
	    h1r = c1 * (a.p_[0 + n3*(i2 + n2*i1)] +
			a.speq_[j2 + (2*n2)*j1] );
	  
	    h1i = c1 * (a.p_[1 + n3*(i2 + n2*i1)] -
			a.speq_[j2+1 + (2*n2)*j1]);
	    
	    h2i = c2 * (a.p_[0 + n3*(i2 + n2*i1)] -
			a.speq_[j2 + (2*n2)*j1] );
	    
	    h2r = -c2 * (a.p_[1 + n3*(i2 + n2*i1)] +
			 a.speq_[j2+1 + (2*n2)*j1]);

	    a.p_[0 + n3*(i2 + n2*i1)] = h1r + h2r; 
	    a.p_[1 + n3*(i2 + n2*i1)] = h1i + h2i; 
	    
	    a.speq_[j2 + (2*n2)*j1] = h1r - h2r; 
	    a.speq_[j2+1 + (2*n2)*j1] = h2i - h1i;
	    
	  } else {
	    j2 = (i2 != 0 ? n2-i2: 0);
	    j3 = n3 - (i3<<1);
	    
	    h1r = c1 * (a.p_[ii3 + n3*(i2 + n2*i1)] +
			a.p_[j3 + n3*(j2 + n2*j1)]);

	    h1i = c1 * (a.p_[ii3+1 + n3*(i2 + n2*i1)] -
			a.p_[ j3+1 + n3*(j2 + n2*j1)]);
	    
	    h2i = c2 * (a.p_[ii3 + n3*(i2 + n2*i1)] -
			a.p_[ j3 + n3*(j2 + n2*j1)]);

	    h2r =-c2 * (a.p_[ii3+1 + n3*(i2 + n2*i1)] +
			a.p_[ j3+1 + n3*(j2 + n2*j1)]);

	    a.p_[ii3   + n3*(i2 + n2*i1)] = h1r + wr*h2r - wi*h2i; 
	    a.p_[ii3+1 + n3*(i2 + n2*i1)] = h1i + wr*h2i + wi*h2r; 
	    a.p_[j3    + n3*(j2 + n2*j1)] = h1r - wr*h2r + wi*h2i; 
	    a.p_[j3+1  + n3*(j2 + n2*j1)] =-h1i + wr*h2i + wi*h2r; 

	  }
	}
	
	wr = (wtmp = wr) * wpr - wi*wpi + wr;
	wi = wi*wpr + wtmp*wpi + wi;
      }
    }

    if(isign == -1){
	// fft3d view real a as packed complex array!
      fftw_complex *in = reinterpret_cast<fftw_complex*>(a.p_);

      ::fftwnd_plan plan 
	  = ::fftw3d_create_plan(dimHalf[0], dimHalf[1], dimHalf[2],
				 FFTW_BACKWARD,
				 FFTW_ESTIMATE|FFTW_IN_PLACE);

      // this function always causes segmentation fault on ia64-linux
      // might be a bug for fftw.
      /*
      fftwnd_plan plan = fftwnd_create_plan(3, reinterpret_cast<int*>(const_cast<size_t*>(dimHalf)), dir,
					    FFTW_ESTIMATE|FFTW_IN_PLACE);
      */
      ::fftwnd_one(plan, in, static_cast<fftw_complex*>(0));
      fftwnd_destroy_plan(plan);
    }

    if(SCALE_FLAG){
      double scale_factor = 2./(n1*n2*n3);
      for(double *ptr = a.p_; ptr != a.p_ + a.size(); ++ptr){
	(*ptr) *= scale_factor;
      }
      for(double *ptr = a.speq_; ptr != a.speq_ + a.SizeOfSpeq(); ++ptr){
	(*ptr) *= scale_factor;
      }
    }
  }


  /**********************************************************************
   * checkSymmetry --
   **********************************************************************/
  template<class VectorA>
  int checkSymmetry(const size_t& n1,
		    const size_t& n2,
		    const size_t& n3,
		    const VectorA& a)
  {
    const size_t length1D = n1*n2*n3;
    size_t differentElementCnts = 0;
    size_t conjugateElementCnts = 0;
    size_t symmetryElementCnts = 0;
    size_t negativeSymmetryCnts = 0;
    double EPSILON = 1e-6;

    //assert(a.size() == length1D);
 
    for(size_t i3=0; i3<n3; ++i3){
      size_t j3 = (i3 == 0) ? 0:n3-i3;
      for(size_t i2=0; i2<n2; ++i2){
	size_t j2 = (i2 == 0) ? 0:n2-i2;
	size_t i1 = 0;
	while(i1 < n1){
	  size_t j1 = (i1 == 0) ? 0:n1-i1;
	  if(i3 != j3 || i2 != j2 || i1 != j1){
	    ++differentElementCnts;
	    size_t idx0 = i3 + n3*(i2 + n2*i1);
	    size_t idx1 = j3 + n3*(j2 + n2*j1);
	    if(fabs(a[idx0].re - a[idx1].re) < EPSILON &&
	       fabs(a[idx0].im - a[idx1].im) < EPSILON ){
	      ++symmetryElementCnts;
	    } else if(fabs(a[idx0].re - a[idx1].re) < EPSILON &&
		      fabs(a[idx0].im + a[idx1].im) < EPSILON ){
	      ++conjugateElementCnts;
	    } else if(fabs(a[idx0].re + a[idx1].re) < EPSILON &&
		      fabs(a[idx0].im + a[idx1].im) < EPSILON ){
	      ++negativeSymmetryCnts;
	    } else {
	      cout<<"a["<<i1<<"]["<<i2<<"]["<<i3<<"] = ("<<a[idx0].re;
	      cout<<","<<a[idx0].im<<")"<<endl;
	      cout<<"a["<<j1<<"]["<<j2<<"]["<<j3<<"]=  ("<<a[idx1].re;
	      cout<<","<<a[idx1].im<<")"<<endl;
	      return 0;
	    }
	  }
	  ++ i1;
	}
      }
    }

    if(symmetryElementCnts == differentElementCnts){
      return +1;
    } else if(conjugateElementCnts == differentElementCnts){
      return -1;
    } else if(negativeSymmetryCnts == differentElementCnts){
      return +2;
    } else {
      return 0;
    }
  }

  template<class VectorA>
  void print3DMatrix(const size_t& n1, 
		     const size_t& n2,
		     const size_t& n3,
		     const VectorA& a)
  {
    for(size_t i1=0; i1<n1; ++i1){
      for(size_t i2=0; i2<n2; ++i2){
	for(size_t i3=0; i3<n3; ++i3){
	  size_t idx = i3 + n3*(i2 + n2*i1);
	  cout<<"a["<<i1<<"]["<<i2<<"]["<<i3<<"] = ("<<a[idx].re;
	  cout<<","<<a[idx].im<<") ";
	} cout<<endl;
      } cout<<endl;
    }
    
  }
} // namespace

#endif










