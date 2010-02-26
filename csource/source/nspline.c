/* nspline.f -- translated by f2c (version 20050501).
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



/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Chuanjie Wu and Jay William Ponder  ## */
/*     ##                    All Rights Reserved                     ## */
/*     ################################################################ */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##   subroutine nspline  --  nonperiodic natural cubic spline   ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "nspline" computes coefficients for an nonperiodic cubic spline */
/*     with natural boundary conditions where the first and last second */
/*     derivatives are already known */


/* Subroutine */ int nspline_(integer *n, integer *np, doublereal *x0, 
	doublereal *y0, doublereal *s1, doublereal *s2, doublereal *h__, 
	doublereal *g, doublereal *dy, doublereal *dla, doublereal *dmu)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal t, y21, y2n;



/*     set first and last second deriviatives to zero */

    y21 = 0.;
    y2n = 0.;

/*     find the intervals to be used */

    i__1 = *n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	h__[i__] = x0[i__ + 1] - x0[i__];
	dy[i__] = (y0[i__ + 1] - y0[i__]) / h__[i__];
    }

/*     calculate the spline coeffcients */

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dla[i__] = h__[i__] / (h__[i__] + h__[i__ - 1]);
	dmu[i__] = 1. - dla[i__];
	g[i__] = (dla[i__] * dy[i__ - 1] + dmu[i__] * dy[i__]) * 3.;
    }

/*     set the initial value via natural boundary condition */

    dla[*n] = 1.;
    dla[0] = 0.;
    dmu[*n] = 0.;
    dmu[0] = 1.;
    g[0] = dy[0] * 3. - h__[0] * .5 * y21;
    g[*n] = dy[*n - 1] * 3. + h__[*n - 1] * .5 * y2n;

/*     solve the triagonal system of linear equations */

    dmu[0] *= .5;
    g[0] *= .5;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = 2. - dmu[i__ - 1] * dla[i__];
	dmu[i__] /= t;
	g[i__] = (g[i__] - g[i__ - 1] * dla[i__]) / t;
    }
    for (i__ = *n - 1; i__ >= 0; --i__) {
	g[i__] -= dmu[i__] * g[i__ + 1];
    }

/*     get the first derivative at each grid point */

    i__1 = *n;
    for (i__ = 0; i__ <= i__1; ++i__) {
	s1[i__] = g[i__];
    }

/*     get the second derivative at each grid point */

    s2[0] = y21;
    s2[*n] = y2n;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2[i__] = (y0[i__ + 1] - y0[i__]) * 6. / (h__[i__] * h__[i__]) - s1[
		i__] * 4. / h__[i__] - s1[i__ + 1] * 2. / h__[i__];
    }
    return 0;
} /* nspline_ */

