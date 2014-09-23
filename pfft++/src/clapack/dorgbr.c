#include "f2c.h"

/* Subroutine */ int dorgbr_(char *vect, integer *m, integer *n, integer *k, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORGBR generates one of the real orthogonal matrices Q or P**T   
    determined by DGEBRD when reducing a real matrix A to bidiagonal   
    form: A = Q * B * P**T.  Q and P**T are defined as products of   
    elementary reflectors H(i) or G(i) respectively.   

    If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q   
    is of order M:   
    if m >= k, Q = H(1) H(2) . . . H(k) and DORGBR returns the first n   
    columns of Q, where m >= n >= k;   
    if m < k, Q = H(1) H(2) . . . H(m-1) and DORGBR returns Q as an   
    M-by-M matrix.   

    If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T   
    is of order N:   
    if k < n, P**T = G(k) . . . G(2) G(1) and DORGBR returns the first m 
  
    rows of P**T, where n >= m >= k;   
    if k >= n, P**T = G(n-1) . . . G(2) G(1) and DORGBR returns P**T as   
    an N-by-N matrix.   

    Arguments   
    =========   

    VECT    (input) CHARACTER*1   
            Specifies whether the matrix Q or the matrix P**T is   
            required, as defined in the transformation applied by DGEBRD: 
  
            = 'Q':  generate Q;   
            = 'P':  generate P**T.   

    M       (input) INTEGER   
            The number of rows of the matrix Q or P**T to be returned.   
            M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q or P**T to be returned. 
  
            N >= 0.   
            If VECT = 'Q', M >= N >= min(M,K);   
            if VECT = 'P', N >= M >= min(N,K).   

    K       (input) INTEGER   
            If VECT = 'Q', the number of columns in the original M-by-K   
            matrix reduced by DGEBRD.   
            If VECT = 'P', the number of rows in the original K-by-N   
            matrix reduced by DGEBRD.   
            K >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the vectors which define the elementary reflectors, 
  
            as returned by DGEBRD.   
            On exit, the M-by-N matrix Q or P**T.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,M).   

    TAU     (input) DOUBLE PRECISION array, dimension   
                                  (min(M,K)) if VECT = 'Q'   
                                  (min(N,K)) if VECT = 'P'   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i) or G(i), which determines Q or P**T, as   
            returned by DGEBRD in its array argument TAUQ or TAUP.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= max(1,min(M,N)).   
            For optimum performance LWORK >= min(M,N)*NB, where NB   
            is the optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer iinfo;
    static logical wantq;
    extern /* Subroutine */ int xerbla_(char *, integer *), dorglq_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), dorgqr_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);


#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    wantq = lsame_(vect, "Q");
    if (! wantq && ! lsame_(vect, "P")) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*n < 0 || wantq && (*n > *m || *n < min(*m,*k)) || ! wantq && (
	    *m > *n || *m < min(*n,*k))) {
	*info = -3;
    } else if (*k < 0) {
	*info = -4;
    } else if (*lda < max(1,*m)) {
	*info = -6;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = min(*m,*n);
	if (*lwork < max(i__1,i__2)) {
	    *info = -9;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORGBR", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	WORK(1) = 1.;
	return 0;
    }

    if (wantq) {

/*        Form Q, determined by a call to DGEBRD to reduce an m-by-k 
  
          matrix */

	if (*m >= *k) {

/*           If m >= k, assume m >= n >= k */

	    dorgqr_(m, n, k, &A(1,1), lda, &TAU(1), &WORK(1), lwork, &
		    iinfo);

	} else {

/*           If m < k, assume m = n   

             Shift the vectors which define the elementary reflect
ors one   
             column to the right, and set the first row and column
 of Q   
             to those of the unit matrix */

	    for (j = *m; j >= 2; --j) {
		A(1,j) = 0.;
		i__1 = *m;
		for (i = j + 1; i <= *m; ++i) {
		    A(i,j) = A(i,j-1);
/* L10: */
		}
/* L20: */
	    }
	    A(1,1) = 1.;
	    i__1 = *m;
	    for (i = 2; i <= *m; ++i) {
		A(i,1) = 0.;
/* L30: */
	    }
	    if (*m > 1) {

/*              Form Q(2:m,2:m) */

		i__1 = *m - 1;
		i__2 = *m - 1;
		i__3 = *m - 1;
		dorgqr_(&i__1, &i__2, &i__3, &A(2,2), lda, &TAU(
			1), &WORK(1), lwork, &iinfo);
	    }
	}
    } else {

/*        Form P', determined by a call to DGEBRD to reduce a k-by-n 
  
          matrix */

	if (*k < *n) {

/*           If k < n, assume k <= m <= n */

	    dorglq_(m, n, k, &A(1,1), lda, &TAU(1), &WORK(1), lwork, &
		    iinfo);

	} else {

/*           If k >= n, assume m = n   

             Shift the vectors which define the elementary reflect
ors one   
             row downward, and set the first row and column of P' 
to   
             those of the unit matrix */

	    A(1,1) = 1.;
	    i__1 = *n;
	    for (i = 2; i <= *n; ++i) {
		A(i,1) = 0.;
/* L40: */
	    }
	    i__1 = *n;
	    for (j = 2; j <= *n; ++j) {
		for (i = j - 1; i >= 2; --i) {
		    A(i,j) = A(i-1,j);
/* L50: */
		}
		A(1,j) = 0.;
/* L60: */
	    }
	    if (*n > 1) {

/*              Form P'(2:n,2:n) */

		i__1 = *n - 1;
		i__2 = *n - 1;
		i__3 = *n - 1;
		dorglq_(&i__1, &i__2, &i__3, &A(2,2), lda, &TAU(
			1), &WORK(1), lwork, &iinfo);
	    }
	}
    }
    return 0;

/*     End of DORGBR */

} /* dorgbr_ */

