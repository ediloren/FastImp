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

  CLAPACK routines

  Resources:

  See also:

  const static char cvsid[] = "$Id: dgesvd.h,v 1.1.1.1 2002/03/12 20:08:56 zhzhu Exp $";

  ==========================================================================
*/

int dgesvd_(char *jobu, char *jobvt, long int  *m, long int  *n, 
	    double *a, long int  *lda, double *s, 
	    double *u, long int  *ldu, double *vt, long int  *ldvt, 
	    double *work, long int *lwork, long int  *info);
