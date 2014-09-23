#include "f2c.h"

/* Subroutine */ int dlasq4_(integer *n, doublereal *q, doublereal *e, 
	doublereal *tau, doublereal *sup)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


       Purpose   
       =======   

       DLASQ4 estimates TAU, the smallest eigenvalue of a matrix. This   
       routine improves the input value of SUP which is an upper bound   
       for the smallest eigenvalue for this matrix .   

       Arguments   
       =========   

    N       (input) INTEGER   
            On entry, N specifies the number of rows and columns   
            in the matrix. N must be at least 0.   

    Q       (input) DOUBLE PRECISION array, dimension (N)   
            Q array   

    E       (input) DOUBLE PRECISION array, dimension (N)   
            E array   

    TAU     (output) DOUBLE PRECISION   
            Estimate of the shift   

    SUP     (input/output) DOUBLE PRECISION   
            Upper bound for the smallest singular value   

    ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublereal c_b4 = .7;
    
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    /* Builtin functions */
    double pow_di(doublereal *, integer *);
    /* Local variables */
    static doublereal xinf, d;
    static integer i;
    static doublereal dm;
    static integer ifl;



#define E(I) e[(I)-1]
#define Q(I) q[(I)-1]


    ifl = 1;
/* Computing MIN */
    d__1 = min(*sup,Q(1)), d__1 = min(d__1,Q(2)), d__1 = min(d__1,Q(3)), d__2 
	    = Q(*n), d__1 = min(d__1,d__2), d__2 = Q(*n - 1), d__1 = min(d__1,
	    d__2), d__2 = Q(*n - 2);
    *sup = min(d__1,d__2);
    *tau = *sup * .9999;
    xinf = 0.;
L10:
    if (ifl == 5) {
	*tau = xinf;
	return 0;
    }
    d = Q(1) - *tau;
    dm = d;
    i__1 = *n - 2;
    for (i = 1; i <= *n-2; ++i) {
	d = d / (d + E(i)) * Q(i + 1) - *tau;
	if (dm > d) {
	    dm = d;
	}
	if (d < 0.) {
	    *sup = *tau;
/* Computing MAX */
	    d__1 = *sup * pow_di(&c_b4, &ifl), d__2 = d + *tau;
	    *tau = max(d__1,d__2);
	    ++ifl;
	    goto L10;
	}
/* L20: */
    }
    d = d / (d + E(*n - 1)) * Q(*n) - *tau;
    if (dm > d) {
	dm = d;
    }
    if (d < 0.) {
	*sup = *tau;
/* Computing MAX */
	d__1 = xinf, d__2 = d + *tau;
	xinf = max(d__1,d__2);
	if (*sup * pow_di(&c_b4, &ifl) <= xinf) {
	    *tau = xinf;
	} else {
	    *tau = *sup * pow_di(&c_b4, &ifl);
	    ++ifl;
	    goto L10;
	}
    } else {
/* Computing MIN */
	d__1 = *sup, d__2 = dm + *tau;
	*sup = min(d__1,d__2);
    }
    return 0;

/*     End of DLASQ4 */

} /* dlasq4_ */

