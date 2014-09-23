#include "f2c.h"

/* Subroutine */ int dlasr_(char *side, char *pivot, char *direct, integer *m,
	 integer *n, doublereal *c, doublereal *s, doublereal *a, integer *
	lda)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASR   performs the transformation   

       A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )   

       A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )   

    where A is an m by n real matrix and P is an orthogonal matrix,   
    consisting of a sequence of plane rotations determined by the   
    parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l' 
  
    and z = n when SIDE = 'R' or 'r' ):   

    When  DIRECT = 'F' or 'f'  ( Forward sequence ) then   

       P = P( z - 1 )*...*P( 2 )*P( 1 ),   

    and when DIRECT = 'B' or 'b'  ( Backward sequence ) then   

       P = P( 1 )*P( 2 )*...*P( z - 1 ),   

    where  P( k ) is a plane rotation matrix for the following planes:   

       when  PIVOT = 'V' or 'v'  ( Variable pivot ),   
          the plane ( k, k + 1 )   

       when  PIVOT = 'T' or 't'  ( Top pivot ),   
          the plane ( 1, k + 1 )   

       when  PIVOT = 'B' or 'b'  ( Bottom pivot ),   
          the plane ( k, z )   

    c( k ) and s( k )  must contain the  cosine and sine that define the 
  
    matrix  P( k ).  The two by two plane rotation part of the matrix   
    P( k ), R( k ), is assumed to be of the form   

       R( k ) = (  c( k )  s( k ) ).   
                ( -s( k )  c( k ) )   

    This version vectorises across rows of the array A when SIDE = 'L'.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            Specifies whether the plane rotation matrix P is applied to   
            A on the left or the right.   
            = 'L':  Left, compute A := P*A   
            = 'R':  Right, compute A:= A*P'   

    DIRECT  (input) CHARACTER*1   
            Specifies whether P is a forward or backward sequence of   
            plane rotations.   
            = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )   
            = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )   

    PIVOT   (input) CHARACTER*1   
            Specifies the plane for which P(k) is a plane rotation   
            matrix.   
            = 'V':  Variable pivot, the plane (k,k+1)   
            = 'T':  Top pivot, the plane (1,k+1)   
            = 'B':  Bottom pivot, the plane (k,z)   

    M       (input) INTEGER   
            The number of rows of the matrix A.  If m <= 1, an immediate 
  
            return is effected.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  If n <= 1, an   
            immediate return is effected.   

    C, S    (input) DOUBLE PRECISION arrays, dimension   
                    (M-1) if SIDE = 'L'   
                    (N-1) if SIDE = 'R'   
            c(k) and s(k) contain the cosine and sine that define the   
            matrix P(k).  The two by two plane rotation part of the   
            matrix P(k), R(k), is assumed to be of the form   
            R( k ) = (  c( k )  s( k ) ).   
                     ( -s( k )  c( k ) )   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            The m by n matrix A.  On exit, A is overwritten by P*A if   
            SIDE = 'R' or by A*P' if SIDE = 'L'.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static doublereal ctemp, stemp;
    extern /* Subroutine */ int xerbla_(char *, integer *);


#define C(I) c[(I)-1]
#define S(I) s[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! (lsame_(side, "L") || lsame_(side, "R"))) {
	info = 1;
    } else if (! (lsame_(pivot, "V") || lsame_(pivot, "T") || 
	    lsame_(pivot, "B"))) {
	info = 2;
    } else if (! (lsame_(direct, "F") || lsame_(direct, "B")))
	     {
	info = 3;
    } else if (*m < 0) {
	info = 4;
    } else if (*n < 0) {
	info = 5;
    } else if (*lda < max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("DLASR ", &info);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }
    if (lsame_(side, "L")) {

/*        Form  P * A */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= *m-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i = 1; i <= *n; ++i) {
			    temp = A(j+1,i);
			    A(j+1,i) = ctemp * temp - stemp * A(j,i);
			    A(j,i) = stemp * temp + ctemp * A(j,i);
/* L10: */
			}
		    }
/* L20: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i = 1; i <= *n; ++i) {
			    temp = A(j+1,i);
			    A(j+1,i) = ctemp * temp - stemp * A(j,i);
			    A(j,i) = stemp * temp + ctemp * A(j,i);
/* L30: */
			}
		    }
/* L40: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m;
		for (j = 2; j <= *m; ++j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i = 1; i <= *n; ++i) {
			    temp = A(j,i);
			    A(j,i) = ctemp * temp - stemp * A(1,i);
			    A(1,i) = stemp * temp + ctemp * A(1,i);
/* L50: */
			}
		    }
/* L60: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m; j >= 2; --j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i = 1; i <= *n; ++i) {
			    temp = A(j,i);
			    A(j,i) = ctemp * temp - stemp * A(1,i);
			    A(1,i) = stemp * temp + ctemp * A(1,i);
/* L70: */
			}
		    }
/* L80: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= *m-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i = 1; i <= *n; ++i) {
			    temp = A(j,i);
			    A(j,i) = stemp * A(*m,i) + 
				    ctemp * temp;
			    A(*m,i) = ctemp * A(*m,i) - 
				    stemp * temp;
/* L90: */
			}
		    }
/* L100: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i = 1; i <= *n; ++i) {
			    temp = A(j,i);
			    A(j,i) = stemp * A(*m,i) + 
				    ctemp * temp;
			    A(*m,i) = ctemp * A(*m,i) - 
				    stemp * temp;
/* L110: */
			}
		    }
/* L120: */
		}
	    }
	}
    } else if (lsame_(side, "R")) {

/*        Form A * P' */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= *n-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i = 1; i <= *m; ++i) {
			    temp = A(i,j+1);
			    A(i,j+1) = ctemp * temp - stemp * 
				    A(i,j);
			    A(i,j) = stemp * temp + ctemp * A(i,j);
/* L130: */
			}
		    }
/* L140: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i = 1; i <= *m; ++i) {
			    temp = A(i,j+1);
			    A(i,j+1) = ctemp * temp - stemp * 
				    A(i,j);
			    A(i,j) = stemp * temp + ctemp * A(i,j);
/* L150: */
			}
		    }
/* L160: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n;
		for (j = 2; j <= *n; ++j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i = 1; i <= *m; ++i) {
			    temp = A(i,j);
			    A(i,j) = ctemp * temp - stemp * A(i,1);
			    A(i,1) = stemp * temp + ctemp * A(i,1);
/* L170: */
			}
		    }
/* L180: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n; j >= 2; --j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i = 1; i <= *m; ++i) {
			    temp = A(i,j);
			    A(i,j) = ctemp * temp - stemp * A(i,1);
			    A(i,1) = stemp * temp + ctemp * A(i,1);
/* L190: */
			}
		    }
/* L200: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= *n-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i = 1; i <= *m; ++i) {
			    temp = A(i,j);
			    A(i,j) = stemp * A(i,*n) + 
				    ctemp * temp;
			    A(i,*n) = ctemp * A(i,*n) - 
				    stemp * temp;
/* L210: */
			}
		    }
/* L220: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i = 1; i <= *m; ++i) {
			    temp = A(i,j);
			    A(i,j) = stemp * A(i,*n) + 
				    ctemp * temp;
			    A(i,*n) = ctemp * A(i,*n) - 
				    stemp * temp;
/* L230: */
			}
		    }
/* L240: */
		}
	    }
	}
    }

    return 0;

/*     End of DLASR */

} /* dlasr_ */

