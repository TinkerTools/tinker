/* bspline.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ######################################################### */
/*     ##                                                     ## */
/*     ##  subroutine bspline  --  get B-spline coefficients  ## */
/*     ##                                                     ## */
/*     ######################################################### */


/*     "bspline" calculates the coefficients for an n-th order */
/*     B-spline approximation */


/* Subroutine */ int bspline_(doublereal *x, integer *n, doublereal *c__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k;
    static doublereal denom;



/*     initialize the B-spline as the linear case */

    /* Parameter adjustments */
    --c__;

    /* Function Body */
    c__[1] = 1. - *x;
    c__[2] = *x;

/*     compute standard B-spline recursion to n-th order */

    i__1 = *n;
    for (k = 3; k <= i__1; ++k) {
	denom = 1. / (doublereal) (k - 1);
	c__[k] = *x * c__[k - 1] * denom;
	i__2 = k - 2;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    c__[k - i__] = ((*x + (doublereal) i__) * c__[k - i__ - 1] + ((
		    doublereal) (k - i__) - *x) * c__[k - i__]) * denom;
	}
	c__[1] = (1. - *x) * c__[1] * denom;
    }
    return 0;
} /* bspline_ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine bspline1  --  get B-spline coeffs & derivs  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "bspline1" calculates the coefficients and derivative */
/*     coefficients for an n-th order B-spline approximation */


/* Subroutine */ int bspline1_(doublereal *x, integer *n, doublereal *c__, 
	doublereal *d__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k;
    static doublereal denom;



/*     initialize the B-spline as the linear case */

    /* Parameter adjustments */
    --d__;
    --c__;

    /* Function Body */
    c__[1] = 1. - *x;
    c__[2] = *x;

/*     compute standard B-spline recursion to n-1-th order */

    i__1 = *n - 1;
    for (k = 3; k <= i__1; ++k) {
	denom = 1. / (doublereal) (k - 1);
	c__[k] = *x * c__[k - 1] * denom;
	i__2 = k - 2;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    c__[k - i__] = ((*x + (doublereal) i__) * c__[k - i__ - 1] + ((
		    doublereal) (k - i__) - *x) * c__[k - i__]) * denom;
	}
	c__[1] = (1. - *x) * c__[1] * denom;
    }

/*     get the derivative from n-1-th order coefficients */

    d__[1] = -c__[1];
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	d__[i__] = c__[i__ - 1] - c__[i__];
    }
    d__[*n] = c__[*n - 1];

/*     use one final recursion to get n-th order coefficients */

    denom = 1. / (doublereal) (*n - 1);
    c__[*n] = *x * c__[*n - 1] * denom;
    i__1 = *n - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[*n - i__] = ((*x + (doublereal) i__) * c__[*n - i__ - 1] + ((
		doublereal) (*n - i__) - *x) * c__[*n - i__]) * denom;
    }
    c__[1] = (1. - *x) * c__[1] * denom;
    return 0;
} /* bspline1_ */

