/* cholesky.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine cholesky  --  modified Cholesky linear solver  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "cholesky" uses a modified Cholesky method to solve the linear */
/*     system Ax = b, returning "x" in "b"; "A" is assumed to be a */
/*     real symmetric positive definite matrix with its diagonal and */
/*     upper triangle stored by rows */

/*     literature reference: */

/*     R. S. Martin, G. Peters and J. H. Wilkinson, Numerische */
/*     Mathematik, 7, 362-383 (1965) */


/* Subroutine */ int cholesky_(integer *nvar, doublereal *a, doublereal *b)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal r__, s, t;
    static integer ii, ij, ik, ki, kk, im, jk, jm;



/*     Cholesky factorization to reduce "A" to (L)(D)(L transpose) */
/*     "L" has a unit diagonal; store 1.0/D on the diagonal of "A" */

    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    ii = 1;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	im = i__ - 1;
	if (i__ != 1) {
	    ij = i__;
	    i__2 = im;
	    for (j = 1; j <= i__2; ++j) {
		r__ = a[ij];
		if (j != 1) {
		    ik = i__;
		    jk = j;
		    jm = j - 1;
		    i__3 = jm;
		    for (k = 1; k <= i__3; ++k) {
			r__ -= a[ik] * a[jk];
			ik = *nvar - k + ik;
			jk = *nvar - k + jk;
		    }
		}
		a[ij] = r__;
		ij = *nvar - j + ij;
	    }
	}
	r__ = a[ii];
	if (i__ != 1) {
	    kk = 1;
	    ik = i__;
	    i__2 = im;
	    for (k = 1; k <= i__2; ++k) {
		s = a[ik];
		t = s * a[kk];
		a[ik] = t;
		r__ -= s * t;
		ik = *nvar - k + ik;
		kk = *nvar - k + 1 + kk;
	    }
	}
	a[ii] = 1. / r__;
	ii = *nvar - i__ + 1 + ii;
    }

/*     solve linear equations; first solve Ly = b for y */

    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ != 1) {
	    ik = i__;
	    im = i__ - 1;
	    r__ = b[i__];
	    i__2 = im;
	    for (k = 1; k <= i__2; ++k) {
		r__ -= b[k] * a[ik];
		ik = *nvar - k + ik;
	    }
	    b[i__] = r__;
	}
    }

/*     finally, solve (D)(L transpose)(x) = y for x */

    ii = *nvar * (*nvar + 1) / 2;
    i__1 = *nvar;
    for (j = 1; j <= i__1; ++j) {
	i__ = *nvar + 1 - j;
	r__ = b[i__] * a[ii];
	if (j != 1) {
	    im = i__ + 1;
	    ki = ii + 1;
	    i__2 = *nvar;
	    for (k = im; k <= i__2; ++k) {
		r__ -= a[ki] * b[k];
		++ki;
	    }
	}
	b[i__] = r__;
	ii = ii - j - 1;
    }
    return 0;
} /* cholesky_ */

