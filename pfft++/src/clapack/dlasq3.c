#include "f2c.h"

/* Subroutine */ int dlasq3_(integer *n, doublereal *q, doublereal *e, 
	doublereal *qq, doublereal *ee, doublereal *sup, doublereal *sigma, 
	integer *kend, integer *off, integer *iphase, integer *iconv, 
	doublereal *eps, doublereal *tol2, doublereal *small2)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


       Purpose   
       =======   

       DLASQ3 is the workhorse of the whole bidiagonal SVD algorithm.   
       This can be described as the differential qd with shifts.   

       Arguments   
       =========   

    N       (input/output) INTEGER   
            On entry, N specifies the number of rows and columns   
            in the matrix. N must be at least 3.   
            On exit N is non-negative and less than the input value.   

    Q       (input/output) DOUBLE PRECISION array, dimension (N)   
            Q array in ping (see IPHASE below)   

    E       (input/output) DOUBLE PRECISION array, dimension (N)   
            E array in ping (see IPHASE below)   

    QQ      (input/output) DOUBLE PRECISION array, dimension (N)   
            Q array in pong (see IPHASE below)   

    EE      (input/output) DOUBLE PRECISION array, dimension (N)   
            E array in pong (see IPHASE below)   

    SUP     (input/output) DOUBLE PRECISION   
            Upper bound for the smallest eigenvalue   

    SIGMA   (input/output) DOUBLE PRECISION   
            Accumulated shift for the present submatrix   

    KEND    (input/output) INTEGER   
            Index where minimum D(i) occurs in recurrence for   
            splitting criterion   

    OFF     (input/output) INTEGER   
            Offset for arrays   

    IPHASE  (input/output) INTEGER   
            If IPHASE = 1 (ping) then data is in Q and E arrays   
            If IPHASE = 2 (pong) then data is in QQ and EE arrays   

    ICONV   (input) INTEGER   
            If ICONV = 0 a bottom part of a matrix (with a split)   
            If ICONV =-3 a top part of a matrix (with a split)   

    EPS     (input) DOUBLE PRECISION   
            Machine epsilon   

    TOL2    (input) DOUBLE PRECISION   
            Square of the relative tolerance TOL as defined in DLASQ1   

    SMALL2  (input) DOUBLE PRECISION   
            A threshold value as defined in DLASQ1   

    ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static logical ldef;
    static integer icnt;
    static doublereal tolx, toly, tolz;
    static integer k1end, k2end;
    static doublereal d;
    static integer i;
    static doublereal qemax;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer maxit, n1, n2;
    static doublereal t1;
    extern /* Subroutine */ int dlasq4_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static integer ic, ke;
    static doublereal dm;
    static integer ip, ks;
    static doublereal xx, yy;
    static logical lsplit;
    static integer ifl;
    static doublereal tau;
    static integer isp;



#define EE(I) ee[(I)-1]
#define QQ(I) qq[(I)-1]
#define E(I) e[(I)-1]
#define Q(I) q[(I)-1]


    icnt = 0;
    tau = 0.;
    dm = *sup;
    tolx = *sigma * *tol2;
    tolz = max(*small2,*sigma) * *tol2;

/*     Set maximum number of iterations */

    maxit = *n * 100;

/*     Flipping */

    ic = 2;
    if (*n > 3) {
	if (*iphase == 1) {
	    i__1 = *n - 2;
	    for (i = 1; i <= *n-2; ++i) {
		if (Q(i) > Q(i + 1)) {
		    ++ic;
		}
		if (E(i) > E(i + 1)) {
		    ++ic;
		}
/* L10: */
	    }
	    if (Q(*n - 1) > Q(*n)) {
		++ic;
	    }
	    if (ic < *n) {
		dcopy_(n, &Q(1), &c__1, &QQ(1), &c_n1);
		i__1 = *n - 1;
		dcopy_(&i__1, &E(1), &c__1, &EE(1), &c_n1);
		if (*kend != 0) {
		    *kend = *n - *kend + 1;
		}
		*iphase = 2;
	    }
	} else {
	    i__1 = *n - 2;
	    for (i = 1; i <= *n-2; ++i) {
		if (QQ(i) > QQ(i + 1)) {
		    ++ic;
		}
		if (EE(i) > EE(i + 1)) {
		    ++ic;
		}
/* L20: */
	    }
	    if (QQ(*n - 1) > QQ(*n)) {
		++ic;
	    }
	    if (ic < *n) {
		dcopy_(n, &QQ(1), &c__1, &Q(1), &c_n1);
		i__1 = *n - 1;
		dcopy_(&i__1, &EE(1), &c__1, &E(1), &c_n1);
		if (*kend != 0) {
		    *kend = *n - *kend + 1;
		}
		*iphase = 1;
	    }
	}
    }
    if (*iconv == -3) {
	if (*iphase == 1) {
	    goto L180;
	} else {
	    goto L80;
	}
    }
    if (*iphase == 2) {
	goto L130;
    }

/*     The ping section of the code */

L30:
    ifl = 0;

/*     Compute the shift */

    if (*kend == 0 || *sup == 0.) {
	tau = 0.;
    } else if (icnt > 0 && dm <= tolz) {
	tau = 0.;
    } else {
/* Computing MAX */
	i__1 = 5, i__2 = *n / 32;
	ip = max(i__1,i__2);
	n2 = (ip << 1) + 1;
	if (n2 >= *n) {
	    n1 = 1;
	    n2 = *n;
	} else if (*kend + ip > *n) {
	    n1 = *n - (ip << 1);
	} else if (*kend - ip < 1) {
	    n1 = 1;
	} else {
	    n1 = *kend - ip;
	}
	dlasq4_(&n2, &Q(n1), &E(n1), &tau, sup);
    }
L40:
    ++icnt;
    if (icnt > maxit) {
	*sup = -1.;
	return 0;
    }
    if (tau == 0.) {

/*     dqd algorithm */

	d = Q(1);
	dm = d;
	ke = 0;
	i__1 = *n - 3;
	for (i = 1; i <= *n-3; ++i) {
	    QQ(i) = d + E(i);
	    d = d / QQ(i) * Q(i + 1);
	    if (dm > d) {
		dm = d;
		ke = i;
	    }
/* L50: */
	}
	++ke;

/*     Penultimate dqd step (in ping) */

	k2end = ke;
	QQ(*n - 2) = d + E(*n - 2);
	d = d / QQ(*n - 2) * Q(*n - 1);
	if (dm > d) {
	    dm = d;
	    ke = *n - 1;
	}

/*     Final dqd step (in ping) */

	k1end = ke;
	QQ(*n - 1) = d + E(*n - 1);
	d = d / QQ(*n - 1) * Q(*n);
	if (dm > d) {
	    dm = d;
	    ke = *n;
	}
	QQ(*n) = d;
    } else {

/*     The dqds algorithm (in ping) */

	d = Q(1) - tau;
	dm = d;
	ke = 0;
	if (d < 0.) {
	    goto L120;
	}
	i__1 = *n - 3;
	for (i = 1; i <= *n-3; ++i) {
	    QQ(i) = d + E(i);
	    d = d / QQ(i) * Q(i + 1) - tau;
	    if (dm > d) {
		dm = d;
		ke = i;
		if (d < 0.) {
		    goto L120;
		}
	    }
/* L60: */
	}
	++ke;

/*     Penultimate dqds step (in ping) */

	k2end = ke;
	QQ(*n - 2) = d + E(*n - 2);
	d = d / QQ(*n - 2) * Q(*n - 1) - tau;
	if (dm > d) {
	    dm = d;
	    ke = *n - 1;
	    if (d < 0.) {
		goto L120;
	    }
	}

/*     Final dqds step (in ping) */

	k1end = ke;
	QQ(*n - 1) = d + E(*n - 1);
	d = d / QQ(*n - 1) * Q(*n) - tau;
	if (dm > d) {
	    dm = d;
	    ke = *n;
	}
	QQ(*n) = d;
    }

/*        Convergence when QQ(N) is small (in ping) */

    if ((d__1 = QQ(*n), abs(d__1)) <= *sigma * *tol2) {
	QQ(*n) = 0.;
	dm = 0.;
	ke = *n;
    }
    if (QQ(*n) < 0.) {
	goto L120;
    }

/*     Non-negative qd array: Update the e's */

    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
	EE(i) = E(i) / QQ(i) * Q(i + 1);
/* L70: */
    }

/*     Updating sigma and iphase in ping */

    *sigma += tau;
    *iphase = 2;
L80:
    tolx = *sigma * *tol2;
    toly = *sigma * *eps;
    tolz = max(*sigma,*small2) * *tol2;

/*     Checking for deflation and convergence (in ping) */

L90:
    if (*n <= 2) {
	return 0;
    }

/*        Deflation: bottom 1x1 (in ping) */

    ldef = FALSE_;
    if (EE(*n - 1) <= tolz) {
	ldef = TRUE_;
    } else if (*sigma > 0.) {
	if (EE(*n - 1) <= *eps * (*sigma + QQ(*n))) {
	    if (EE(*n - 1) * (QQ(*n) / (QQ(*n) + *sigma)) <= *tol2 * (QQ(*n) 
		    + *sigma)) {
		ldef = TRUE_;
	    }
	}
    } else {
	if (EE(*n - 1) <= QQ(*n) * *tol2) {
	    ldef = TRUE_;
	}
    }
    if (ldef) {
	Q(*n) = QQ(*n) + *sigma;
	--(*n);
	++(*iconv);
	goto L90;
    }

/*        Deflation: bottom 2x2 (in ping) */

    ldef = FALSE_;
    if (EE(*n - 2) <= tolz) {
	ldef = TRUE_;
    } else if (*sigma > 0.) {
	t1 = *sigma + EE(*n - 1) * (*sigma / (*sigma + QQ(*n)));
	if (EE(*n - 2) * (t1 / (QQ(*n - 1) + t1)) <= toly) {
	    if (EE(*n - 2) * (QQ(*n - 1) / (QQ(*n - 1) + t1)) <= tolx) {
		ldef = TRUE_;
	    }
	}
    } else {
	if (EE(*n - 2) <= QQ(*n) / (QQ(*n) + EE(*n - 1) + QQ(*n - 1)) * QQ(*n 
		- 1) * *tol2) {
	    ldef = TRUE_;
	}
    }
    if (ldef) {
/* Computing MAX */
	d__1 = QQ(*n), d__2 = QQ(*n - 1), d__1 = max(d__1,d__2), d__2 = EE(*n 
		- 1);
	qemax = max(d__1,d__2);
	if (qemax != 0.) {
	    if (qemax == QQ(*n - 1)) {
/* Computing 2nd power */
		d__1 = (QQ(*n) - QQ(*n - 1) + EE(*n - 1)) / qemax;
		xx = (QQ(*n) + QQ(*n - 1) + EE(*n - 1) + qemax * sqrt(d__1 * 
			d__1 + EE(*n - 1) * 4. / qemax)) * .5;
	    } else if (qemax == QQ(*n)) {
/* Computing 2nd power */
		d__1 = (QQ(*n - 1) - QQ(*n) + EE(*n - 1)) / qemax;
		xx = (QQ(*n) + QQ(*n - 1) + EE(*n - 1) + qemax * sqrt(d__1 * 
			d__1 + EE(*n - 1) * 4. / qemax)) * .5;
	    } else {
/* Computing 2nd power */
		d__1 = (QQ(*n) - QQ(*n - 1) + EE(*n - 1)) / qemax;
		xx = (QQ(*n) + QQ(*n - 1) + EE(*n - 1) + qemax * sqrt(d__1 * 
			d__1 + QQ(*n - 1) * 4. / qemax)) * .5;
	    }
/* Computing MAX */
	    d__1 = QQ(*n), d__2 = QQ(*n - 1);
/* Computing MIN */
	    d__3 = QQ(*n), d__4 = QQ(*n - 1);
	    yy = max(d__1,d__2) / xx * min(d__3,d__4);
	} else {
	    xx = 0.;
	    yy = 0.;
	}
	Q(*n - 1) = *sigma + xx;
	Q(*n) = yy + *sigma;
	*n += -2;
	*iconv += 2;
	goto L90;
    }

/*     Updating bounds before going to pong */

    if (*iconv == 0) {
	*kend = ke;
/* Computing MIN */
	d__1 = dm, d__2 = *sup - tau;
	*sup = min(d__1,d__2);
    } else if (*iconv > 0) {
/* Computing MIN */
	d__1 = QQ(*n), d__2 = QQ(*n - 1), d__1 = min(d__1,d__2), d__2 = QQ(*n 
		- 2), d__1 = min(d__1,d__2), d__1 = min(d__1,QQ(1)), d__1 = 
		min(d__1,QQ(2));
	*sup = min(d__1,QQ(3));
	if (*iconv == 1) {
	    *kend = k1end;
	} else if (*iconv == 2) {
	    *kend = k2end;
	} else {
	    *kend = *n;
	}
	icnt = 0;
	maxit = *n * 100;
    }

/*     Checking for splitting in ping */

    lsplit = FALSE_;
    for (ks = *n - 3; ks >= 3; --ks) {
	if (EE(ks) <= toly) {
/* Computing MIN */
	    d__1 = QQ(ks + 1), d__2 = QQ(ks);
/* Computing MIN */
	    d__3 = QQ(ks + 1), d__4 = QQ(ks);
	    if (EE(ks) * (min(d__1,d__2) / (min(d__3,d__4) + *sigma)) <= tolx)
		     {
		lsplit = TRUE_;
		goto L110;
	    }
	}
/* L100: */
    }

    ks = 2;
    if (EE(2) <= tolz) {
	lsplit = TRUE_;
    } else if (*sigma > 0.) {
	t1 = *sigma + EE(1) * (*sigma / (*sigma + QQ(1)));
	if (EE(2) * (t1 / (QQ(1) + t1)) <= toly) {
	    if (EE(2) * (QQ(1) / (QQ(1) + t1)) <= tolx) {
		lsplit = TRUE_;
	    }
	}
    } else {
	if (EE(2) <= QQ(1) / (QQ(1) + EE(1) + QQ(2)) * QQ(2) * *tol2) {
	    lsplit = TRUE_;
	}
    }
    if (lsplit) {
	goto L110;
    }

    ks = 1;
    if (EE(1) <= tolz) {
	lsplit = TRUE_;
    } else if (*sigma > 0.) {
	if (EE(1) <= *eps * (*sigma + QQ(1))) {
	    if (EE(1) * (QQ(1) / (QQ(1) + *sigma)) <= *tol2 * (QQ(1) + *sigma)
		    ) {
		lsplit = TRUE_;
	    }
	}
    } else {
	if (EE(1) <= QQ(1) * *tol2) {
	    lsplit = TRUE_;
	}
    }

L110:
    if (lsplit) {
/* Computing MIN */
	d__1 = QQ(*n), d__2 = QQ(*n - 1), d__1 = min(d__1,d__2), d__2 = QQ(*n 
		- 2);
	*sup = min(d__1,d__2);
	isp = -(*off + 1);
	*off += ks;
	*n -= ks;
/* Computing MAX */
	i__1 = 1, i__2 = *kend - ks;
	*kend = max(i__1,i__2);
	E(ks) = *sigma;
	EE(ks) = (doublereal) isp;
	*iconv = 0;
	return 0;
    }

/*     Coincidence */

    if (tau == 0. && dm <= tolz && *kend != *n && *iconv == 0 && icnt > 0) {
	i__1 = *n - ke;
	dcopy_(&i__1, &E(ke), &c__1, &QQ(ke), &c__1);
	QQ(*n) = 0.;
	i__1 = *n - ke;
	dcopy_(&i__1, &Q(ke + 1), &c__1, &EE(ke), &c__1);
	*sup = 0.;
    }
    *iconv = 0;
    goto L130;

/*     A new shift when the previous failed (in ping) */

L120:
    ++ifl;
    *sup = tau;

/*     SUP is small or   
       Too many bad shifts (ping) */

    if (*sup <= tolz || ifl >= 2) {
	tau = 0.;
	goto L40;

/*     The asymptotic shift (in ping) */

    } else {
/* Computing MAX */
	d__1 = tau + d;
	tau = max(d__1,0.);
	if (tau <= tolz) {
	    tau = 0.;
	}
	goto L40;
    }

/*     the pong section of the code */

L130:
    ifl = 0;

/*     Compute the shift (in pong) */

    if (*kend == 0 && *sup == 0.) {
	tau = 0.;
    } else if (icnt > 0 && dm <= tolz) {
	tau = 0.;
    } else {
/* Computing MAX */
	i__1 = 5, i__2 = *n / 32;
	ip = max(i__1,i__2);
	n2 = (ip << 1) + 1;
	if (n2 >= *n) {
	    n1 = 1;
	    n2 = *n;
	} else if (*kend + ip > *n) {
	    n1 = *n - (ip << 1);
	} else if (*kend - ip < 1) {
	    n1 = 1;
	} else {
	    n1 = *kend - ip;
	}
	dlasq4_(&n2, &QQ(n1), &EE(n1), &tau, sup);
    }
L140:
    ++icnt;
    if (icnt > maxit) {
	*sup = -(*sup);
	return 0;
    }
    if (tau == 0.) {

/*     The dqd algorithm (in pong) */

	d = QQ(1);
	dm = d;
	ke = 0;
	i__1 = *n - 3;
	for (i = 1; i <= *n-3; ++i) {
	    Q(i) = d + EE(i);
	    d = d / Q(i) * QQ(i + 1);
	    if (dm > d) {
		dm = d;
		ke = i;
	    }
/* L150: */
	}
	++ke;

/*     Penultimate dqd step (in pong) */

	k2end = ke;
	Q(*n - 2) = d + EE(*n - 2);
	d = d / Q(*n - 2) * QQ(*n - 1);
	if (dm > d) {
	    dm = d;
	    ke = *n - 1;
	}

/*     Final dqd step (in pong) */

	k1end = ke;
	Q(*n - 1) = d + EE(*n - 1);
	d = d / Q(*n - 1) * QQ(*n);
	if (dm > d) {
	    dm = d;
	    ke = *n;
	}
	Q(*n) = d;
    } else {

/*     The dqds algorithm (in pong) */

	d = QQ(1) - tau;
	dm = d;
	ke = 0;
	if (d < 0.) {
	    goto L220;
	}
	i__1 = *n - 3;
	for (i = 1; i <= *n-3; ++i) {
	    Q(i) = d + EE(i);
	    d = d / Q(i) * QQ(i + 1) - tau;
	    if (dm > d) {
		dm = d;
		ke = i;
		if (d < 0.) {
		    goto L220;
		}
	    }
/* L160: */
	}
	++ke;

/*     Penultimate dqds step (in pong) */

	k2end = ke;
	Q(*n - 2) = d + EE(*n - 2);
	d = d / Q(*n - 2) * QQ(*n - 1) - tau;
	if (dm > d) {
	    dm = d;
	    ke = *n - 1;
	    if (d < 0.) {
		goto L220;
	    }
	}

/*     Final dqds step (in pong) */

	k1end = ke;
	Q(*n - 1) = d + EE(*n - 1);
	d = d / Q(*n - 1) * QQ(*n) - tau;
	if (dm > d) {
	    dm = d;
	    ke = *n;
	}
	Q(*n) = d;
    }

/*        Convergence when is small (in pong) */

    if ((d__1 = Q(*n), abs(d__1)) <= *sigma * *tol2) {
	Q(*n) = 0.;
	dm = 0.;
	ke = *n;
    }
    if (Q(*n) < 0.) {
	goto L220;
    }

/*     Non-negative qd array: Update the e's */

    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
	E(i) = EE(i) / Q(i) * QQ(i + 1);
/* L170: */
    }

/*     Updating sigma and iphase in pong */

    *sigma += tau;
L180:
    *iphase = 1;
    tolx = *sigma * *tol2;
    toly = *sigma * *eps;

/*     Checking for deflation and convergence (in pong) */

L190:
    if (*n <= 2) {
	return 0;
    }

/*        Deflation: bottom 1x1 (in pong) */

    ldef = FALSE_;
    if (E(*n - 1) <= tolz) {
	ldef = TRUE_;
    } else if (*sigma > 0.) {
	if (E(*n - 1) <= *eps * (*sigma + Q(*n))) {
	    if (E(*n - 1) * (Q(*n) / (Q(*n) + *sigma)) <= *tol2 * (Q(*n) + *
		    sigma)) {
		ldef = TRUE_;
	    }
	}
    } else {
	if (E(*n - 1) <= Q(*n) * *tol2) {
	    ldef = TRUE_;
	}
    }
    if (ldef) {
	Q(*n) += *sigma;
	--(*n);
	++(*iconv);
	goto L190;
    }

/*        Deflation: bottom 2x2 (in pong) */

    ldef = FALSE_;
    if (E(*n - 2) <= tolz) {
	ldef = TRUE_;
    } else if (*sigma > 0.) {
	t1 = *sigma + E(*n - 1) * (*sigma / (*sigma + Q(*n)));
	if (E(*n - 2) * (t1 / (Q(*n - 1) + t1)) <= toly) {
	    if (E(*n - 2) * (Q(*n - 1) / (Q(*n - 1) + t1)) <= tolx) {
		ldef = TRUE_;
	    }
	}
    } else {
	if (E(*n - 2) <= Q(*n) / (Q(*n) + EE(*n - 1) + Q(*n - 1)) * Q(*n - 1) 
		* *tol2) {
	    ldef = TRUE_;
	}
    }
    if (ldef) {
/* Computing MAX */
	d__1 = Q(*n), d__2 = Q(*n - 1), d__1 = max(d__1,d__2), d__2 = E(*n - 
		1);
	qemax = max(d__1,d__2);
	if (qemax != 0.) {
	    if (qemax == Q(*n - 1)) {
/* Computing 2nd power */
		d__1 = (Q(*n) - Q(*n - 1) + E(*n - 1)) / qemax;
		xx = (Q(*n) + Q(*n - 1) + E(*n - 1) + qemax * sqrt(d__1 * 
			d__1 + E(*n - 1) * 4. / qemax)) * .5;
	    } else if (qemax == Q(*n)) {
/* Computing 2nd power */
		d__1 = (Q(*n - 1) - Q(*n) + E(*n - 1)) / qemax;
		xx = (Q(*n) + Q(*n - 1) + E(*n - 1) + qemax * sqrt(d__1 * 
			d__1 + E(*n - 1) * 4. / qemax)) * .5;
	    } else {
/* Computing 2nd power */
		d__1 = (Q(*n) - Q(*n - 1) + E(*n - 1)) / qemax;
		xx = (Q(*n) + Q(*n - 1) + E(*n - 1) + qemax * sqrt(d__1 * 
			d__1 + Q(*n - 1) * 4. / qemax)) * .5;
	    }
/* Computing MAX */
	    d__1 = Q(*n), d__2 = Q(*n - 1);
/* Computing MIN */
	    d__3 = Q(*n), d__4 = Q(*n - 1);
	    yy = max(d__1,d__2) / xx * min(d__3,d__4);
	} else {
	    xx = 0.;
	    yy = 0.;
	}
	Q(*n - 1) = *sigma + xx;
	Q(*n) = yy + *sigma;
	*n += -2;
	*iconv += 2;
	goto L190;
    }

/*     Updating bounds before going to pong */

    if (*iconv == 0) {
	*kend = ke;
/* Computing MIN */
	d__1 = dm, d__2 = *sup - tau;
	*sup = min(d__1,d__2);
    } else if (*iconv > 0) {
/* Computing MIN */
	d__1 = Q(*n), d__2 = Q(*n - 1), d__1 = min(d__1,d__2), d__2 = Q(*n - 
		2), d__1 = min(d__1,d__2), d__1 = min(d__1,Q(1)), d__1 = min(
		d__1,Q(2));
	*sup = min(d__1,Q(3));
	if (*iconv == 1) {
	    *kend = k1end;
	} else if (*iconv == 2) {
	    *kend = k2end;
	} else {
	    *kend = *n;
	}
	icnt = 0;
	maxit = *n * 100;
    }

/*     Checking for splitting in pong */

    lsplit = FALSE_;
    for (ks = *n - 3; ks >= 3; --ks) {
	if (E(ks) <= toly) {
/* Computing MIN */
	    d__1 = Q(ks + 1), d__2 = Q(ks);
/* Computing MIN */
	    d__3 = Q(ks + 1), d__4 = Q(ks);
	    if (E(ks) * (min(d__1,d__2) / (min(d__3,d__4) + *sigma)) <= tolx) 
		    {
		lsplit = TRUE_;
		goto L210;
	    }
	}
/* L200: */
    }

    ks = 2;
    if (E(2) <= tolz) {
	lsplit = TRUE_;
    } else if (*sigma > 0.) {
	t1 = *sigma + E(1) * (*sigma / (*sigma + Q(1)));
	if (E(2) * (t1 / (Q(1) + t1)) <= toly) {
	    if (E(2) * (Q(1) / (Q(1) + t1)) <= tolx) {
		lsplit = TRUE_;
	    }
	}
    } else {
	if (E(2) <= Q(1) / (Q(1) + E(1) + Q(2)) * Q(2) * *tol2) {
	    lsplit = TRUE_;
	}
    }
    if (lsplit) {
	goto L210;
    }

    ks = 1;
    if (E(1) <= tolz) {
	lsplit = TRUE_;
    } else if (*sigma > 0.) {
	if (E(1) <= *eps * (*sigma + Q(1))) {
	    if (E(1) * (Q(1) / (Q(1) + *sigma)) <= *tol2 * (Q(1) + *sigma)) {
		lsplit = TRUE_;
	    }
	}
    } else {
	if (E(1) <= Q(1) * *tol2) {
	    lsplit = TRUE_;
	}
    }

L210:
    if (lsplit) {
/* Computing MIN */
	d__1 = Q(*n), d__2 = Q(*n - 1), d__1 = min(d__1,d__2), d__2 = Q(*n - 
		2);
	*sup = min(d__1,d__2);
	isp = *off + 1;
	*off += ks;
/* Computing MAX */
	i__1 = 1, i__2 = *kend - ks;
	*kend = max(i__1,i__2);
	*n -= ks;
	E(ks) = *sigma;
	EE(ks) = (doublereal) isp;
	*iconv = 0;
	return 0;
    }

/*     Coincidence */

    if (tau == 0. && dm <= tolz && *kend != *n && *iconv == 0 && icnt > 0) {
	i__1 = *n - ke;
	dcopy_(&i__1, &EE(ke), &c__1, &Q(ke), &c__1);
	Q(*n) = 0.;
	i__1 = *n - ke;
	dcopy_(&i__1, &QQ(ke + 1), &c__1, &E(ke), &c__1);
	*sup = 0.;
    }
    *iconv = 0;
    goto L30;

/*     Computation of a new shift when the previous failed (in pong) */

L220:
    ++ifl;
    *sup = tau;

/*     SUP is small or   
       Too many bad shifts (in pong) */

    if (*sup <= tolz || ifl >= 2) {
	tau = 0.;
	goto L140;

/*     The asymptotic shift (in pong) */

    } else {
/* Computing MAX */
	d__1 = tau + d;
	tau = max(d__1,0.);
	if (tau <= tolz) {
	    tau = 0.;
	}
	goto L140;
    }

/*     End of DLASQ3 */

    return 0;
} /* dlasq3_ */

