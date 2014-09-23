/*-*-c-*-*/
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

  const static char cvsid[] = "$Id: dense.h,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

  ==========================================================================
*/

#ifndef _DENSE_H_
#define _DENSE_H_

#ifdef __cplusplus
extern "C" {
#endif

int 
pseudo_inverse(
	       double **A, 
	       int rows, 
	       int cols);


#ifdef __cplusplus
} /* End of extern "C" { */
#endif

#endif
