/* erf.f -- translated by f2c (version 20050501).
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



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function erf  --  evaluate the standard error function  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "erf" computes a numerical approximation to the value of */
/*     the error function via a Chebyshev approximation */


doublereal erf_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer mode;
    static doublereal result;
    extern /* Subroutine */ int erfcore_(doublereal *, doublereal *, integer *
	    );



/*     compute the error function via Chebyshev fitting */

    mode = 0;
    erfcore_(x, &result, &mode);
    ret_val = result;
    return ret_val;
} /* erf_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function erfc  --  evaluate complementary error function  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "erfc" computes a numerical approximation to the value of the */
/*     complementary error function via a Chebyshev approximation */


doublereal erfc_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer mode;
    static doublereal result;
    extern /* Subroutine */ int erfcore_(doublereal *, doublereal *, integer *
	    );



/*     get the complementary error function via Chebyshev fitting */

    mode = 1;
    erfcore_(x, &result, &mode);
    ret_val = result;
    return ret_val;
} /* erfc_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine erfcore  --  erf and erfc via Chebyshev approx  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "erfcore" evaluates erf(x) or erfc(x) for a real argument x; */
/*     when called with mode set to 0 it returns erf, a mode of 1 */
/*     returns erfc; uses rational functions that approximate erf(x) */
/*     and erfc(x) to at least 18 significant decimal digits */

/*     literature reference: */

/*     W. J. Cody, "Rational Chebyshev Approximations for the Error */
/*     Function", Mathematics of Computation, 631-638, 1969 */

/*     adapted from an original program written by W. J. Cody, */
/*     Mathematics and Computer Science Division, Argonne National */
/*     Laboratory, Argonne, IL 60439 */

/*     machine-dependent constants: */

/*     xtiny   argument below which erf(x) may be represented by */
/*             2*x/sqrt(pi) and above which x*x won't underflow; */
/*             a conservative value is the largest machine number */
/*             X such that 1.0 + X = 1.0 to machine precision */

/*     xbig    largest argument acceptable for erfc; solution to */
/*             the equation:  W(x) * (1-0.5/x**2) = XMIN, where */
/*             W(x) = exp(-x*x)/[x*sqrt(pi)] */


/* Subroutine */ int erfcore_(doublereal *arg, doublereal *result, integer *
	mode)
{
    /* Initialized data */

    static doublereal c__[9] = { .564188496988670089,8.88314979438837594,
	    66.1191906371416295,298.635138197400131,881.95222124176909,
	    1712.04761263407058,2051.07837782607147,1230.33935479799725,
	    2.15311535474403846e-8 };
    static doublereal d__[8] = { 15.7449261107098347,117.693950891312499,
	    537.181101862009858,1621.38957456669019,3290.79923573345963,
	    4362.61909014324716,3439.36767414372164,1230.33935480374942 };
    static doublereal p[6] = { .305326634961232344,.360344899949804439,
	    .125781726111229246,.0160837851487422766,6.58749161529837803e-4,
	    .0163153871373020978 };
    static doublereal q[5] = { 2.56852019228982242,1.87295284992346047,
	    .527905102951428412,.0605183413124413191,.00233520497626869185 };
    static doublereal sqrpi = .56418958354775628695;
    static doublereal xtiny = 1.11e-16;
    static doublereal xbig = 26.543;
    static doublereal a[5] = { 3.1611237438705656,113.864154151050156,
	    377.485237685302021,3209.37758913846947,.185777706184603153 };
    static doublereal b[4] = { 23.6012909523441209,244.024637934444173,
	    1282.61652607737228,2844.23683343917062 };

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double d_int(doublereal *), exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal x, y, del, ysq, xden, xnum;


/*     mathematical and machine-dependent constants */


/*     coefficients for approximation to erf in first interval */


/*     coefficients for approximation to erfc in second interval */


/*     coefficients for approximation to erfc in third interval */



/*     store the argument and its absolute value */

    x = *arg;
    y = abs(x);

/*     evaluate error function for |x| less than 0.46875 */

    if (y <= .46875) {
	ysq = 0.;
	if (y > xtiny) {
	    ysq = y * y;
	}
	xnum = a[4] * ysq;
	xden = ysq;
	for (i__ = 1; i__ <= 3; ++i__) {
	    xnum = (xnum + a[i__ - 1]) * ysq;
	    xden = (xden + b[i__ - 1]) * ysq;
	}
	*result = x * (xnum + a[3]) / (xden + b[3]);
	if (*mode != 0) {
	    *result = 1. - *result;
	}

/*     get complementary error function for 0.46875 <= |x| <= 4.0 */

    } else if (y <= 4.) {
	xnum = c__[8] * y;
	xden = y;
	for (i__ = 1; i__ <= 7; ++i__) {
	    xnum = (xnum + c__[i__ - 1]) * y;
	    xden = (xden + d__[i__ - 1]) * y;
	}
	*result = (xnum + c__[7]) / (xden + d__[7]);
	d__1 = y * 16.;
	ysq = d_int(&d__1) / 16.;
	del = (y - ysq) * (y + ysq);
/*        result = exp(-ysq*ysq) * exp(-del) * result */
	*result = exp(-ysq * ysq - del) * *result;
	if (*mode == 0) {
	    *result = 1. - *result;
	    if (x < 0.) {
		*result = -(*result);
	    }
	} else {
	    if (x < 0.) {
		*result = 2. - *result;
	    }
	}

/*     get complementary error function for |x| greater than 4.0 */

    } else {
	*result = 0.;
	if (y < xbig) {
	    ysq = 1. / (y * y);
	    xnum = p[5] * ysq;
	    xden = ysq;
	    for (i__ = 1; i__ <= 4; ++i__) {
		xnum = (xnum + p[i__ - 1]) * ysq;
		xden = (xden + q[i__ - 1]) * ysq;
	    }
	    *result = ysq * (xnum + p[4]) / (xden + q[4]);
	    *result = (sqrpi - *result) / y;
	    d__1 = y * 16.;
	    ysq = d_int(&d__1) / 16.;
	    del = (y - ysq) * (y + ysq);
/*           result = exp(-ysq*ysq) * exp(-del) * result */
	    *result = exp(-ysq * ysq - del) * *result;
	}
	if (*mode == 0) {
	    *result = 1. - *result;
	    if (x < 0.) {
		*result = -(*result);
	    }
	} else {
	    if (x < 0.) {
		*result = 2. - *result;
	    }
	}
    }
    return 0;
} /* erfcore_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function erfinv  --  evaluate the error function inverse  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "erfinv" evaluates the inverse of the error function for */
/*     an argument in the range (-1,1) using a rational function */
/*     approximation followed by cycles of Newton-Raphson correction */

/*     adapted from the pseudocode for the Matlab function of the */
/*     same name; Matlab, version 4.2c, March 1995 */


doublereal erfinv_(doublereal *x)
{
    /* Initialized data */

    static doublereal a[4] = { .886226899,-1.645349621,.914624893,-.140543331 
	    };
    static doublereal b[4] = { -2.118377725,1.442710462,-.329097515,
	    .012229801 };
    static doublereal c__[4] = { -1.970840454,-1.624906493,3.429567803,
	    1.641345311 };
    static doublereal d__[2] = { 3.5438892,1.6370678 };

    /* Format strings */
    static char fmt_10[] = "(/,\002 ERFINV  --  Illegal Argument to Invers"
	    "e\002,\002 Error Function\002)";

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);
    double exp(doublereal);

    /* Local variables */
    static doublereal y, z__;
    extern doublereal erf_(doublereal *);
    extern /* Subroutine */ int fatal_(void);

    /* Fortran I/O blocks */
    static cilist io___27 = { 0, 0, 0, fmt_10, 0 };




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  math.i  --  mathematical and geometrical constants  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     radian   conversion factor from radians to degrees */
/*     pi       numerical value of the geometric constant */
/*     sqrtpi   numerical value of the square root of Pi */
/*     logten   numerical value of the natural log of ten */
/*     sqrttwo  numerical value of the square root of two */
/*     twosix   numerical value of the sixth root of two */



/*     coefficients for approximation to erfinv in central range */


/*     coefficients for approximation to erfinv near endpoints */



/*     get an initial estimate for the inverse error function */

    if (abs(*x) <= .7) {
	y = *x * *x;
	z__ = *x * (((a[3] * y + a[2]) * y + a[1]) * y + a[0]) / ((((b[3] * y 
		+ b[2]) * y + b[1]) * y + b[0]) * y + 1.);
    } else if (*x > .7 && *x < 1.) {
	y = sqrt(-log((1. - *x) / 2.));
	z__ = (((c__[3] * y + c__[2]) * y + c__[1]) * y + c__[0]) / ((d__[1] *
		 y + d__[0]) * y + 1.);
    } else if (*x < -.7 && *x > -1.) {
	y = sqrt(-log((*x + 1.) / 2.));
	z__ = -(((c__[3] * y + c__[2]) * y + c__[1]) * y + c__[0]) / ((d__[1] 
		* y + d__[0]) * y + 1.);
    } else {
	io___27.ciunit = iounit_1.iout;
	s_wsfe(&io___27);
	e_wsfe();
	fatal_();
    }

/*     use two steps of Newton-Raphson correction to increase accuracy */

    z__ -= (erf_(&z__) - *x) / (exp(-z__ * z__) * 1.1283791670955126);
    z__ -= (erf_(&z__) - *x) / (exp(-z__ * z__) * 1.1283791670955126);
    ret_val = z__;
    return ret_val;
} /* erfinv_ */

