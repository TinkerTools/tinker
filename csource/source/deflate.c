/* deflate.f -- translated by f2c (version 20050501).
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

/* Common Block Declarations */

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

/* Table of constant values */

static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine deflate  --  eigenvalues by method of deflation  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "deflate" uses the power method with deflation to compute the */
/*     few largest eigenvalues and eigenvectors of a symmetric matrix */

/*     n     logical dimension of the matrix to be diagonalized */
/*     np    physical dimension of the matrix storage area */
/*     nv    number of largest eigenvalues to be extracted */
/*     a     input with the matrix to be diagonalized; only */
/*              the lower triangle and diagonal are required */
/*     ev    returned with the eigenvalues in descending order */
/*     vec   returned with the eigenvectors of the matrix */
/*     work  temporary work vector */


/* Subroutine */ int deflate_(integer *n, integer *np, integer *nv, 
	doublereal *a, doublereal *ev, doublereal *vec, doublereal *work)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 DEFLATE  --  Eigenvalue\002,i3,\002 not "
	    "Fully Converged\002)";

    /* System generated locals */
    integer vec_dim1, vec_offset, a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k;
    static doublereal eps, dot1, dot2;
    static integer iter;
    static doublereal ratio;
    extern doublereal random_(void);
    static integer maxiter;

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 0, 0, fmt_10, 0 };



#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define vec_ref(a_1,a_2) vec[(a_2)*vec_dim1 + a_1]



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     initialize number of iterations and convergence criteria */

    /* Parameter adjustments */
    --work;
    a_dim1 = *np;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    vec_dim1 = *np;
    vec_offset = 1 + vec_dim1;
    vec -= vec_offset;
    --ev;

    /* Function Body */
    maxiter = 500;
    eps = 1e-6;

/*     use identity vector as initial guess for eigenvectors */

    i__1 = *nv;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    vec_ref(i__, j) = 1.;
	}
    }

/*     find the few largest eigenvalues and eigenvectors */

    i__1 = *nv;
    for (k = 1; k <= i__1; ++k) {
	ev[k] = 0.;
	dot1 = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[i__] = 0.;
	    i__3 = i__ - 1;
	    for (j = 1; j <= i__3; ++j) {
		work[i__] += a_ref(i__, j) * vec_ref(j, k);
	    }
	    i__3 = *n;
	    for (j = i__; j <= i__3; ++j) {
		work[i__] += a_ref(j, i__) * vec_ref(j, k);
	    }
/* Computing 2nd power */
	    d__1 = work[i__];
	    dot1 += d__1 * d__1;
	}

/*     if in or near null space, use random guess as eigenvector */

	if (dot1 <= eps * 100. * (doublereal) (*n)) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		work[i__] = random_();
	    }
	}

/*     find the current eigenvalue by iterating to convergence; */
/*     first multiply vector by matrix and compute dot products */

	i__2 = maxiter;
	for (iter = 1; iter <= i__2; ++iter) {
	    dot1 = 0.;
	    dot2 = 0.;
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		vec_ref(i__, k) = 0.;
		i__4 = i__ - 1;
		for (j = 1; j <= i__4; ++j) {
		    vec_ref(i__, k) = vec_ref(i__, k) + a_ref(i__, j) * work[
			    j];
		}
		i__4 = *n;
		for (j = i__; j <= i__4; ++j) {
		    vec_ref(i__, k) = vec_ref(i__, k) + a_ref(j, i__) * work[
			    j];
		}
/* Computing 2nd power */
		d__1 = vec_ref(i__, k);
		dot1 += d__1 * d__1;
		dot2 += vec_ref(i__, k) * work[i__];
	    }

/*     normalize new eigenvector and substitute for old one */

	    ratio = (d__1 = (ev[k] - dot2) / dot2, abs(d__1));
	    ev[k] = dot2;
	    dot1 = sqrt(dot1);
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		vec_ref(i__, k) = vec_ref(i__, k) / dot1;
		work[i__] = vec_ref(i__, k);
	    }
	    if (ratio < eps) {
		goto L20;
	    }
	}
	io___10.ciunit = iounit_1.iout;
	s_wsfe(&io___10);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	e_wsfe();

/*     eliminate the current eigenvalue from the matrix */

L20:
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *n;
	    for (j = i__; j <= i__3; ++j) {
		a_ref(j, i__) = a_ref(j, i__) - ev[k] * vec_ref(i__, k) * 
			vec_ref(j, k);
	    }
	}
    }
    return 0;
} /* deflate_ */

#undef vec_ref
#undef a_ref


