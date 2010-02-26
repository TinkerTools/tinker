/* invbeta.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function invbeta  --  inverse Beta distribution function  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "invbeta" computes the inverse Beta distribution function */
/*     via a combination of Newton iteration and bisection search */

/*     literature reference: */

/*     K. L. Majumder and G. P. Bhattacharjee, "Inverse of the */
/*     Incomplete Beta Function Ratio", Applied Statistics, 22, */
/*     411-414 (1973) */


doublereal invbeta_(doublereal *a, doublereal *b, doublereal *y)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal), pow_dd(doublereal *, doublereal 
	    *);

    /* Local variables */
    static doublereal x, x0, x1, beta, mean;
    static logical done;
    static doublereal aexp, bexp;
    extern doublereal betai_(doublereal *, doublereal *, doublereal *);
    static doublereal slope, stdev, error;
    extern doublereal gammln_(doublereal *);



/*     use limiting values when input argument is out of range */

    done = FALSE_;
    if (*y <= 0.) {
	x = 0.;
	done = TRUE_;
    } else if (*y >= 1.) {
	x = 1.;
	done = TRUE_;
    }

/*     initial guess from mean and variance of probability function */

    if (! done) {
	aexp = *a - 1.;
	bexp = *b - 1.;
	d__1 = *a + *b;
	beta = exp(gammln_(a) + gammln_(b) - gammln_(&d__1));
	mean = *a / (*a + *b);
/* Computing 2nd power */
	d__1 = *a + *b;
	stdev = sqrt(*a * *b / ((*a + *b + 1.) * (d__1 * d__1)));
	if (*y > 0. && *y <= .167) {
	    x = mean + (*y / .167 - 2.) * stdev;
	} else if (*y > .167 && *y < .833) {
	    x = mean + (*y / .333 - 1.5) * stdev;
	} else if (*y >= .833 && *y < 1.) {
	    x = mean + (*y / .167 - 4.) * stdev;
	}
/* Computing MAX */
	d__1 = 1e-5, d__2 = min(.99999000000000005,x);
	x = max(d__1,d__2);
    }

/*     refine inverse distribution value via Newton iteration */

    while(! done) {
	d__1 = 1. - x;
	slope = pow_dd(&x, &aexp) * pow_dd(&d__1, &bexp) / beta;
	error = betai_(a, b, &x) - *y;
	x -= error / slope;
	if (abs(error) < 1e-5) {
	    done = TRUE_;
	}
	if (x < 0. || x > 1.) {
	    done = TRUE_;
	}
    }

/*     try bisection search if Newton iteration moved out of range */

    if (x < 0. || x > 1.) {
	x0 = 0.;
	x1 = 1.;
	done = FALSE_;
    }

/*     refine inverse distribution value via bisection search */

    while(! done) {
	x = (x0 + x1) * .5;
	error = betai_(a, b, &x) - *y;
	if (error > 0.) {
	    x1 = x;
	}
	if (error < 0.) {
	    x0 = x;
	}
	if (abs(error) < 1e-5) {
	    done = TRUE_;
	}
    }

/*     return best estimate of the inverse beta distribution value */

    ret_val = x;
    return ret_val;
} /* invbeta_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  function betai  --  cumulative Beta distribution function  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "betai" evaluates the cumulative Beta distribution function */
/*     as the probability that a random variable from a distribution */
/*     with Beta parameters "a" and "b" will be less than "x" */


doublereal betai_(doublereal *a, doublereal *b, doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal bt;
    extern doublereal betacf_(doublereal *, doublereal *, doublereal *), 
	    gammln_(doublereal *);



/*     get cumulative distribution directly or via reflection */

    if (*x <= 0.) {
	ret_val = 0.;
    } else if (*x >= 1.) {
	ret_val = 1.;
    } else {
	d__1 = *a + *b;
	bt = exp(gammln_(&d__1) - gammln_(a) - gammln_(b) + *a * log(*x) + *b 
		* log(1. - *x));
	if (*x < (*a + 1.) / (*a + *b + 2.)) {
	    ret_val = bt / *a * betacf_(a, b, x);
	} else {
	    d__1 = 1. - *x;
	    ret_val = 1. - bt / *b * betacf_(b, a, &d__1);
	}
    }
    return ret_val;
} /* betai_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  function betacf  --  continued fraction routine for betai  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "betacf" computes a rapidly convergent continued fraction needed */
/*     by routine "betai" to evaluate the cumulative Beta distribution */


doublereal betacf_(doublereal *a, doublereal *b, doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal c__, d__, h__;
    static integer i__;
    static doublereal m, m2, aa, qab, del, qam, qap;



/*     establish an initial guess for the Beta continued fraction */

    qab = *a + *b;
    qap = *a + 1.;
    qam = *a - 1.;
    c__ = 1.;
    d__ = 1. - qab * *x / qap;
    if (abs(d__) < 1e-30) {
	d__ = 1e-30;
    }
    d__ = 1. / d__;
    h__ = d__;

/*     iteratively improve the continued fraction to convergence */

    for (i__ = 1; i__ <= 100; ++i__) {
	m = (doublereal) i__;
	m2 = m * 2.;
	aa = m * (*b - m) * *x / ((qam + m2) * (*a + m2));
	d__ = aa * d__ + 1.;
	if (abs(d__) < 1e-30) {
	    d__ = 1e-30;
	}
	c__ = aa / c__ + 1.;
	if (abs(c__) < 1e-30) {
	    c__ = 1e-30;
	}
	d__ = 1. / d__;
	h__ = h__ * d__ * c__;
	aa = -(*a + m) * (qab + m) * *x / ((*a + m2) * (qap + m2));
	d__ = aa * d__ + 1.;
	if (abs(d__) < 1e-30) {
	    d__ = 1e-30;
	}
	c__ = aa / c__ + 1.;
	if (abs(c__) < 1e-30) {
	    c__ = 1e-30;
	}
	d__ = 1. / d__;
	del = d__ * c__;
	h__ *= del;
	if ((d__1 = del - 1., abs(d__1)) < 1e-10) {
	    goto L10;
	}
    }
L10:
    ret_val = h__;
    return ret_val;
} /* betacf_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function gammln  --  natural log of the Gamma function  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "gammln" uses a series expansion due to Lanczos to compute */
/*     the natural logarithm of the Gamma function at "x" in [0,1] */


doublereal gammln_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal temp, series;



/*     get the natural log of Gamma via a series expansion */

    temp = *x + 5.5;
    temp = (*x + .5) * log(temp) - temp;
    series = 76.18009172947146 / (*x + 1.) + 1.000000000190015 + 
	    -86.50532032941677 / (*x + 2.) + 24.01409824083091 / (*x + 3.) + 
	    -1.231739572450155 / (*x + 4.) + .001208650973866179 / (*x + 5.) 
	    + -5.395239384953e-6 / (*x + 6.);
    ret_val = temp + log(series * 2.5066282746310005 / *x);
    return ret_val;
} /* gammln_ */

