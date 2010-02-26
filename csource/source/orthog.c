/* orthog.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine orthog  --  Gram-Schmidt orthogonalization  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "orthog" performs an orthogonalization of an input matrix */
/*     via the modified Gram-Schmidt algorithm */

/*     variables and parameters: */

/*     m     first logical dimension of matrix to orthogonalize */
/*     mp    first physical dimension of matrix storage area */
/*     n     second logical dimension of matrix to orthogonalize */
/*     a     matrix to orthogonalize; contains result on exit */


/* Subroutine */ int orthog_(integer *m, integer *mp, integer *n, doublereal *
	a)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal rkj, rkk;


#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]



/*     compute the modified Gram-Schmidt orthogonalization */

    /* Parameter adjustments */
    a_dim1 = *mp;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	rkk = 0.;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    d__1 = a_ref(i__, k);
	    rkk += d__1 * d__1;
	}
	rkk = sqrt(rkk);
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    a_ref(i__, k) = a_ref(i__, k) / rkk;
	}
	i__2 = *n;
	for (j = k + 1; j <= i__2; ++j) {
	    rkj = 0.;
	    i__3 = *m;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		rkj += a_ref(i__, k) * a_ref(i__, j);
	    }
	    i__3 = *m;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) - a_ref(i__, k) * rkj;
	    }
	}
    }
    return 0;
} /* orthog_ */

#undef a_ref


