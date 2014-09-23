#include "f2c.h"

/* Subroutine */ int dlasq1_(integer *n, doublereal *d, doublereal *e, 
	doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


       Purpose   
       =======   

       DLASQ1 computes the singular values of a real N-by-N bidiagonal   
       matrix with diagonal D and off-diagonal E. The singular values are 
  
       computed to high relative accuracy, barring over/underflow or   
       denormalization. The algorithm is described in   

       "Accurate singular values and differential qd algorithms," by   
       K. V. Fernando and B. N. Parlett,   
       Numer. Math., Vol-67, No. 2, pp. 191-230,1994.   

       See also   
       "Implementation of differential qd algorithms," by   
       K. V. Fernando and B. N. Parlett, Technical Report,   
       Department of Mathematics, University of California at Berkeley,   
       1994 (Under preparation).   

       Arguments   
       =========   

    N       (input) INTEGER   
            The number of rows and columns in the matrix. N >= 0.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, D contains the diagonal elements of the   
            bidiagonal matrix whose SVD is desired. On normal exit,   
            D contains the singular values in decreasing order.   

    E       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, elements E(1:N-1) contain the off-diagonal elements 
  
            of the bidiagonal matrix whose SVD is desired.   
            On exit, E is overwritten.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the algorithm did not converge;  i   
                  specifies how many superdiagonals did not converge.   

    ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublereal c_b8 = .125;
    static integer c__1 = 1;
    static integer c__0 = 0;
    
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;
    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);
    /* Local variables */
    static integer kend, ierr;
    extern /* Subroutine */ int dlas2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    static integer i, j, m;
    static doublereal sfmin, sigmn;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal sigmx;
    extern /* Subroutine */ int dlasq2_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *);
    static doublereal small2;
    static integer ke;
    static doublereal dm;
    extern doublereal dlamch_(char *);
    static doublereal dx;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *);
    static integer ny;
    extern /* Subroutine */ int xerbla_(char *, integer *), dlasrt_(
	    char *, integer *, doublereal *, integer *);
    static doublereal thresh, tolmul;
    static logical restrt;
    static doublereal scl, eps, tol, sig1, sig2, tol2;



#define WORK(I) work[(I)-1]
#define E(I) e[(I)-1]
#define D(I) d[(I)-1]


    *info = 0;
    if (*n < 0) {
	*info = -2;
	i__1 = -(*info);
	xerbla_("DLASQ1", &i__1);
	return 0;
    } else if (*n == 0) {
	return 0;
    } else if (*n == 1) {
	D(1) = abs(D(1));
	return 0;
    } else if (*n == 2) {
	dlas2_(&D(1), &E(1), &D(2), &sigmn, &sigmx);
	D(1) = sigmx;
	D(2) = sigmn;
	return 0;
    }

/*     Estimate the largest singular value */

    sigmx = 0.;
    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
/* Computing MAX */
	d__2 = sigmx, d__3 = (d__1 = E(i), abs(d__1));
	sigmx = max(d__2,d__3);
/* L10: */
    }

/*     Early return if sigmx is zero (matrix is already diagonal) */

    if (sigmx == 0.) {
	goto L70;
    }

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	D(i) = (d__1 = D(i), abs(d__1));
/* Computing MAX */
	d__1 = sigmx, d__2 = D(i);
	sigmx = max(d__1,d__2);
/* L20: */
    }

/*     Get machine parameters */

    eps = dlamch_("EPSILON");
    sfmin = dlamch_("SAFE MINIMUM");

/*     Compute singular values to relative accuracy TOL   
       It is assumed that tol**2 does not underflow.   

   Computing MAX   
   Computing MIN */
    d__3 = 100., d__4 = pow_dd(&eps, &c_b8);
    d__1 = 10., d__2 = min(d__3,d__4);
    tolmul = max(d__1,d__2);
    tol = tolmul * eps;
/* Computing 2nd power */
    d__1 = tol;
    tol2 = d__1 * d__1;

    thresh = sigmx * sqrt(sfmin) * tol;

/*     Scale matrix so the square of the largest element is   
       1 / ( 256 * SFMIN ) */

    scl = sqrt(1. / (sfmin * 256.));
/* Computing 2nd power */
    d__1 = tolmul;
    small2 = 1. / (d__1 * d__1 * 256.);
    dcopy_(n, &D(1), &c__1, &WORK(1), &c__1);
    i__1 = *n - 1;
    dcopy_(&i__1, &E(1), &c__1, &WORK(*n + 1), &c__1);
    dlascl_("G", &c__0, &c__0, &sigmx, &scl, n, &c__1, &WORK(1), n, &ierr)
	    ;
    i__1 = *n - 1;
    i__2 = *n - 1;
    dlascl_("G", &c__0, &c__0, &sigmx, &scl, &i__1, &c__1, &WORK(*n + 1), &
	    i__2, &ierr);

/*     Square D and E (the input for the qd algorithm) */

    i__1 = (*n << 1) - 1;
    for (j = 1; j <= (*n<<1)-1; ++j) {
/* Computing 2nd power */
	d__1 = WORK(j);
	WORK(j) = d__1 * d__1;
/* L30: */
    }

/*     Apply qd algorithm */

    m = 0;
    E(*n) = 0.;
    dx = WORK(1);
    dm = dx;
    ke = 0;
    restrt = FALSE_;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if ((d__1 = E(i), abs(d__1)) <= thresh || WORK(*n + i) <= tol2 * (dm /
		 (doublereal) (i - m))) {
	    ny = i - m;
	    if (ny == 1) {
		goto L50;
	    } else if (ny == 2) {
		dlas2_(&D(m + 1), &E(m + 1), &D(m + 2), &sig1, &sig2);
		D(m + 1) = sig1;
		D(m + 2) = sig2;
	    } else {
		kend = ke + 1 - m;
		dlasq2_(&ny, &D(m + 1), &E(m + 1), &WORK(m + 1), &WORK(m + *n 
			+ 1), &eps, &tol2, &small2, &dm, &kend, info);

/*                 Return, INFO = number of unconverged superd
iagonals */

		if (*info != 0) {
		    *info += i;
		    return 0;
		}

/*                 Undo scaling */

		i__2 = m + ny;
		for (j = m + 1; j <= m+ny; ++j) {
		    D(j) = sqrt(D(j));
/* L40: */
		}
		dlascl_("G", &c__0, &c__0, &scl, &sigmx, &ny, &c__1, &D(m + 1)
			, &ny, &ierr);
	    }
L50:
	    m = i;
	    if (i != *n) {
		dx = WORK(i + 1);
		dm = dx;
		ke = i;
		restrt = TRUE_;
	    }
	}
	if (i != *n && ! restrt) {
	    dx = WORK(i + 1) * (dx / (dx + WORK(*n + i)));
	    if (dm > dx) {
		dm = dx;
		ke = i;
	    }
	}
	restrt = FALSE_;
/* L60: */
    }
    kend = ke + 1;

/*     Sort the singular values into decreasing order */

L70:
    dlasrt_("D", n, &D(1), info);
    return 0;

/*     End of DLASQ1 */

} /* dlasq1_ */

