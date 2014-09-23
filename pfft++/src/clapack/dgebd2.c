#include "f2c.h"

/* Subroutine */ int dgebd2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *d, doublereal *e, doublereal *tauq, doublereal *taup,
	 doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DGEBD2 reduces a real general m by n matrix A to upper or lower   
    bidiagonal form B by an orthogonal transformation: Q' * A * P = B.   

    If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows in the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns in the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the m by n general matrix to be reduced.   
            On exit,   
            if m >= n, the diagonal and the first superdiagonal are   
              overwritten with the upper bidiagonal matrix B; the   
              elements below the diagonal, with the array TAUQ, represent 
  
              the orthogonal matrix Q as a product of elementary   
              reflectors, and the elements above the first superdiagonal, 
  
              with the array TAUP, represent the orthogonal matrix P as   
              a product of elementary reflectors;   
            if m < n, the diagonal and the first subdiagonal are   
              overwritten with the lower bidiagonal matrix B; the   
              elements below the first subdiagonal, with the array TAUQ, 
  
              represent the orthogonal matrix Q as a product of   
              elementary reflectors, and the elements above the diagonal, 
  
              with the array TAUP, represent the orthogonal matrix P as   
              a product of elementary reflectors.   
            See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    D       (output) DOUBLE PRECISION array, dimension (min(M,N))   
            The diagonal elements of the bidiagonal matrix B:   
            D(i) = A(i,i).   

    E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)   
            The off-diagonal elements of the bidiagonal matrix B:   
            if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;   
            if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.   

    TAUQ    (output) DOUBLE PRECISION array dimension (min(M,N))   
            The scalar factors of the elementary reflectors which   
            represent the orthogonal matrix Q. See Further Details.   

    TAUP    (output) DOUBLE PRECISION array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors which   
            represent the orthogonal matrix P. See Further Details.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (max(M,N))   

    INFO    (output) INTEGER   
            = 0: successful exit.   
            < 0: if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    The matrices Q and P are represented as products of elementary   
    reflectors:   

    If m >= n,   

       Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)   

    Each H(i) and G(i) has the form:   

       H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'   

    where tauq and taup are real scalars, and v and u are real vectors;   
    v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i); 
  
    u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n); 
  
    tauq is stored in TAUQ(i) and taup in TAUP(i).   

    If m < n,   

       Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)   

    Each H(i) and G(i) has the form:   

       H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'   

    where tauq and taup are real scalars, and v and u are real vectors;   
    v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i); 
  
    u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n); 
  
    tauq is stored in TAUQ(i) and taup in TAUP(i).   

    The contents of A on exit are illustrated by the following examples: 
  

    m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):   

      (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )   
      (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )   
      (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )   
      (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )   
      (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )   
      (  v1  v2  v3  v4  v5 )   

    where d and e denote diagonal and off-diagonal elements of B, vi   
    denotes an element of the vector defining H(i), and ui an element of 
  
    the vector defining G(i).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    static integer i;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *), dlarfg_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *), xerbla_(char *, integer *);



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define TAUQ(I) tauq[(I)-1]
#define TAUP(I) taup[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    }
    if (*info < 0) {
	i__1 = -(*info);
	xerbla_("DGEBD2", &i__1);
	return 0;
    }

    if (*m >= *n) {

/*        Reduce to upper bidiagonal form */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {

/*           Generate elementary reflector H(i) to annihilate A(i+
1:m,i) */

	    i__2 = *m - i + 1;
/* Computing MIN */
	    i__3 = i + 1;
	    dlarfg_(&i__2, &A(i,i), &A(min(i+1,*m),i), 
		    &c__1, &TAUQ(i));
	    D(i) = A(i,i);
	    A(i,i) = 1.;

/*           Apply H(i) to A(i:m,i+1:n) from the left */

	    i__2 = *m - i + 1;
	    i__3 = *n - i;
	    dlarf_("Left", &i__2, &i__3, &A(i,i), &c__1, &TAUQ(i), 
		    &A(i,i+1), lda, &WORK(1));
	    A(i,i) = D(i);

	    if (i < *n) {

/*              Generate elementary reflector G(i) to annihila
te   
                A(i,i+2:n) */

		i__2 = *n - i;
/* Computing MIN */
		i__3 = i + 2;
		dlarfg_(&i__2, &A(i,i+1), &A(i,min(i+2,*n)), lda, &TAUP(i));
		E(i) = A(i,i+1);
		A(i,i+1) = 1.;

/*              Apply G(i) to A(i+1:m,i+1:n) from the right */

		i__2 = *m - i;
		i__3 = *n - i;
		dlarf_("Right", &i__2, &i__3, &A(i,i+1), lda, &
			TAUP(i), &A(i+1,i+1), lda, &WORK(1));
		A(i,i+1) = E(i);
	    } else {
		TAUP(i) = 0.;
	    }
/* L10: */
	}
    } else {

/*        Reduce to lower bidiagonal form */

	i__1 = *m;
	for (i = 1; i <= *m; ++i) {

/*           Generate elementary reflector G(i) to annihilate A(i,
i+1:n) */

	    i__2 = *n - i + 1;
/* Computing MIN */
	    i__3 = i + 1;
	    dlarfg_(&i__2, &A(i,i), &A(i,min(i+1,*n)), 
		    lda, &TAUP(i));
	    D(i) = A(i,i);
	    A(i,i) = 1.;

/*           Apply G(i) to A(i+1:m,i:n) from the right */

	    i__2 = *m - i;
	    i__3 = *n - i + 1;
/* Computing MIN */
	    i__4 = i + 1;
	    dlarf_("Right", &i__2, &i__3, &A(i,i), lda, &TAUP(i), &
		    A(min(i+1,*m),i), lda, &WORK(1));
	    A(i,i) = D(i);

	    if (i < *m) {

/*              Generate elementary reflector H(i) to annihila
te   
                A(i+2:m,i) */

		i__2 = *m - i;
/* Computing MIN */
		i__3 = i + 2;
		dlarfg_(&i__2, &A(i+1,i), &A(min(i+2,*m),i), &c__1, &TAUQ(i));
		E(i) = A(i+1,i);
		A(i+1,i) = 1.;

/*              Apply H(i) to A(i+1:m,i+1:n) from the left */

		i__2 = *m - i;
		i__3 = *n - i;
		dlarf_("Left", &i__2, &i__3, &A(i+1,i), &c__1, &
			TAUQ(i), &A(i+1,i+1), lda, &WORK(1));
		A(i+1,i) = E(i);
	    } else {
		TAUQ(i) = 0.;
	    }
/* L20: */
	}
    }
    return 0;

/*     End of DGEBD2 */

} /* dgebd2_ */

