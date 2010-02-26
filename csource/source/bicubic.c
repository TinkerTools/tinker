/* bicubic.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine bcuint  --  bicubic interpolation of function  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "bcuint" performs a bicubic interpolation of the function */
/*     value on a 2D spline grid */


/* Subroutine */ int bcuint_(doublereal *y, doublereal *y1, doublereal *y2, 
	doublereal *y12, doublereal *x1l, doublereal *x1u, doublereal *x2l, 
	doublereal *x2u, doublereal *x1, doublereal *x2, doublereal *ansy)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal c__[16]	/* was [4][4] */;
    static integer i__;
    static doublereal t, u;
    extern /* Subroutine */ int bcucof_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);


#define c___ref(a_1,a_2) c__[(a_2)*4 + a_1 - 5]



/*     get coefficients, then perform bicubic interpolation */

    /* Parameter adjustments */
    --y12;
    --y2;
    --y1;
    --y;

    /* Function Body */
    d__1 = *x1u - *x1l;
    d__2 = *x2u - *x2l;
    bcucof_(&y[1], &y1[1], &y2[1], &y12[1], &d__1, &d__2, c__);
    t = (*x1 - *x1l) / (*x1u - *x1l);
    u = (*x2 - *x2l) / (*x2u - *x2l);
    *ansy = 0.;
    for (i__ = 4; i__ >= 1; --i__) {
	*ansy = t * *ansy + ((c___ref(i__, 4) * u + c___ref(i__, 3)) * u + 
		c___ref(i__, 2)) * u + c___ref(i__, 1);
    }
    return 0;
} /* bcuint_ */

#undef c___ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine bcuint1  --  bicubic interpolation of gradient  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "bcuint1" performs a bicubic interpolation of the function */
/*     value and gradient along the directions of a 2D spline grid */

/*     literature reference: */

/*     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. */
/*     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge */
/*     University Press, 1992, Section 3.6 */


/* Subroutine */ int bcuint1_(doublereal *y, doublereal *y1, doublereal *y2, 
	doublereal *y12, doublereal *x1l, doublereal *x1u, doublereal *x2l, 
	doublereal *x2u, doublereal *x1, doublereal *x2, doublereal *ansy, 
	doublereal *ansy1, doublereal *ansy2)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal c__[16]	/* was [4][4] */;
    static integer i__;
    static doublereal t, u;
    extern /* Subroutine */ int bcucof_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);


#define c___ref(a_1,a_2) c__[(a_2)*4 + a_1 - 5]



/*     get coefficients, then perform bicubic interpolation */

    /* Parameter adjustments */
    --y12;
    --y2;
    --y1;
    --y;

    /* Function Body */
    d__1 = *x1u - *x1l;
    d__2 = *x2u - *x2l;
    bcucof_(&y[1], &y1[1], &y2[1], &y12[1], &d__1, &d__2, c__);
    t = (*x1 - *x1l) / (*x1u - *x1l);
    u = (*x2 - *x2l) / (*x2u - *x2l);
    *ansy = 0.;
    *ansy1 = 0.;
    *ansy2 = 0.;
    for (i__ = 4; i__ >= 1; --i__) {
	*ansy = t * *ansy + ((c___ref(i__, 4) * u + c___ref(i__, 3)) * u + 
		c___ref(i__, 2)) * u + c___ref(i__, 1);
	*ansy1 = u * *ansy1 + (c___ref(4, i__) * 3. * t + c___ref(3, i__) * 
		2.) * t + c___ref(2, i__);
	*ansy2 = t * *ansy2 + (c___ref(i__, 4) * 3. * u + c___ref(i__, 3) * 
		2.) * u + c___ref(i__, 2);
    }
    *ansy1 /= *x1u - *x1l;
    *ansy2 /= *x2u - *x2l;
    return 0;
} /* bcuint1_ */

#undef c___ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine bcuint2  --  bicubic interpolation of Hessian  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "bcuint2" performs a bicubic interpolation of the function value, */
/*     gradient and Hessian along the directions of a 2D spline grid */


/* Subroutine */ int bcuint2_(doublereal *y, doublereal *y1, doublereal *y2, 
	doublereal *y12, doublereal *x1l, doublereal *x1u, doublereal *x2l, 
	doublereal *x2u, doublereal *x1, doublereal *x2, doublereal *ansy, 
	doublereal *ansy1, doublereal *ansy2, doublereal *ansy12, doublereal *
	ansy11, doublereal *ansy22)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal c__[16]	/* was [4][4] */;
    static integer i__;
    static doublereal t, u;
    extern /* Subroutine */ int bcucof_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);


#define c___ref(a_1,a_2) c__[(a_2)*4 + a_1 - 5]



/*     get coefficients, then perform bicubic interpolation */

    /* Parameter adjustments */
    --y12;
    --y2;
    --y1;
    --y;

    /* Function Body */
    d__1 = *x1u - *x1l;
    d__2 = *x2u - *x2l;
    bcucof_(&y[1], &y1[1], &y2[1], &y12[1], &d__1, &d__2, c__);
    t = (*x1 - *x1l) / (*x1u - *x1l);
    u = (*x2 - *x2l) / (*x2u - *x2l);
    *ansy = 0.;
    *ansy1 = 0.;
    *ansy2 = 0.;
    *ansy11 = 0.;
    *ansy22 = 0.;
    for (i__ = 4; i__ >= 1; --i__) {
	*ansy = t * *ansy + ((c___ref(i__, 4) * u + c___ref(i__, 3)) * u + 
		c___ref(i__, 2)) * u + c___ref(i__, 1);
	*ansy1 = u * *ansy1 + (c___ref(4, i__) * 3. * t + c___ref(3, i__) * 
		2.) * t + c___ref(2, i__);
	*ansy2 = t * *ansy2 + (c___ref(i__, 4) * 3. * u + c___ref(i__, 3) * 
		2.) * u + c___ref(i__, 2);
	*ansy11 = u * *ansy11 + c___ref(4, i__) * 6. * t + c___ref(3, i__) * 
		2.;
	*ansy22 = t * *ansy22 + c___ref(i__, 4) * 6. * u + c___ref(i__, 3) * 
		2.;
    }
    *ansy12 = t * 3. * t * ((c___ref(4, 4) * 3. * u + c___ref(4, 3) * 2.) * u 
	    + c___ref(4, 2)) + t * 2. * ((c___ref(3, 4) * 3. * u + c___ref(3, 
	    3) * 2.) * u + c___ref(3, 2)) + (c___ref(2, 4) * 3. * u + c___ref(
	    2, 3) * 2.) * u + c___ref(2, 2);
    *ansy1 /= *x1u - *x1l;
    *ansy2 /= *x2u - *x2l;
    *ansy11 /= (*x1u - *x1l) * (*x1u - *x1l);
    *ansy22 /= (*x2u - *x2l) * (*x2u - *x2l);
    *ansy12 /= (*x1u - *x1l) * (*x2u - *x2l);
    return 0;
} /* bcuint2_ */

#undef c___ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine bcucof  --  bicubic interpolation coefficients  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "bcucof" determines the coefficient matrix needed for bicubic */
/*     interpolation of a function, gradients and cross derivatives */

/*     literature reference: */

/*     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. */
/*     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge */
/*     University Press, 1992, Section 3.6 */


/* Subroutine */ int bcucof_(doublereal *y, doublereal *y1, doublereal *y2, 
	doublereal *y12, doublereal *d1, doublereal *d2, doublereal *c__)
{
    /* Initialized data */

    static doublereal wt[256]	/* was [16][16] */ = { 1.,0.,-3.,2.,0.,0.,0.,
	    0.,-3.,0.,9.,-6.,2.,0.,-6.,4.,0.,0.,0.,0.,0.,0.,0.,0.,3.,0.,-9.,
	    6.,-2.,0.,6.,-4.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,9.,-6.,0.,0.,-6.,
	    4.,0.,0.,3.,-2.,0.,0.,0.,0.,0.,0.,-9.,6.,0.,0.,6.,-4.,0.,0.,0.,0.,
	    1.,0.,-3.,2.,-2.,0.,6.,-4.,1.,0.,-3.,2.,0.,0.,0.,0.,0.,0.,0.,0.,
	    -1.,0.,3.,-2.,1.,0.,-3.,2.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-3.,2.,
	    0.,0.,3.,-2.,0.,0.,0.,0.,0.,0.,3.,-2.,0.,0.,-6.,4.,0.,0.,3.,-2.,
	    0.,1.,-2.,1.,0.,0.,0.,0.,0.,-3.,6.,-3.,0.,2.,-4.,2.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,3.,-6.,3.,0.,-2.,4.,-2.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	    0.,-3.,3.,0.,0.,2.,-2.,0.,0.,-1.,1.,0.,0.,0.,0.,0.,0.,3.,-3.,0.,
	    0.,-2.,2.,0.,0.,0.,0.,0.,1.,-2.,1.,0.,-2.,4.,-2.,0.,1.,-2.,1.,0.,
	    0.,0.,0.,0.,0.,0.,0.,0.,-1.,2.,-1.,0.,1.,-2.,1.,0.,0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,1.,-1.,0.,0.,-1.,1.,0.,0.,0.,0.,0.,0.,-1.,1.,0.,0.,2.,
	    -2.,0.,0.,-1.,1. };

    static integer i__, j, k;
    static doublereal x[16], cl[16], xx, d1d2;


#define c___ref(a_1,a_2) c__[(a_2)*4 + a_1]
#define wt_ref(a_1,a_2) wt[(a_2)*16 + a_1 - 17]

    /* Parameter adjustments */
    c__ -= 5;
    --y12;
    --y2;
    --y1;
    --y;

    /* Function Body */


/*     pack a temporary vector of corner values */

    d1d2 = *d1 * *d2;
    for (i__ = 1; i__ <= 4; ++i__) {
	x[i__ - 1] = y[i__];
	x[i__ + 3] = y1[i__] * *d1;
	x[i__ + 7] = y2[i__] * *d2;
	x[i__ + 11] = y12[i__] * d1d2;
    }

/*     matrix multiply by the stored weight table */

    for (i__ = 1; i__ <= 16; ++i__) {
	xx = 0.;
	for (k = 1; k <= 16; ++k) {
	    xx += wt_ref(i__, k) * x[k - 1];
	}
	cl[i__ - 1] = xx;
    }

/*     unpack the result into the coefficient table */

    j = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	for (k = 1; k <= 4; ++k) {
	    ++j;
	    c___ref(i__, k) = cl[j - 1];
	}
    }
    return 0;
} /* bcucof_ */

#undef wt_ref
#undef c___ref


