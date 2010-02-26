/* cspline.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine cspline  --  periodic interpolating cube spline  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "cspline" computes the coefficients for a periodic interpolating */
/*     cubic spline */

/*     literature reference: */

/*     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran, */
/*     Springer Verlag, 1996, Section 10.1.2  [see routine "isplpe"] */


/* Subroutine */ int cspline_(integer *n, integer *np, doublereal *xn, 
	doublereal *fn, doublereal *b, doublereal *c__, doublereal *d__, 
	doublereal *h__, doublereal *du, doublereal *dm, doublereal *rc, 
	doublereal *rs)
{
    /* Format strings */
    static char fmt_10[] = "(\002 CSPLINE  --  Warning, Input Values are No"
	    "t\002,\002 Periodic\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__;
    static doublereal eps, temp1, temp2;
    static integer iflag;
    extern /* Subroutine */ int cytsy_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;
    static doublereal average;

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, fmt_10, 0 };




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




/*     check the periodicity of fn, and for subsequent call */

    eps = 1e-6;
    if ((d__1 = fn[*n] - fn[0], abs(d__1)) > eps) {
	io___2.ciunit = iounit_1.iout;
	s_wsfe(&io___2);
	e_wsfe();
    }
    average = (fn[0] + fn[*n]) * .5;
    fn[0] = average;
    fn[*n] = average;

/*     get auxiliary variables and matrix elements on first call */

    i__1 = *n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	h__[i__] = xn[i__ + 1] - xn[i__];
    }
    h__[*n] = h__[0];
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	du[i__] = h__[i__];
    }
    du[*n] = h__[0];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dm[i__] = (h__[i__ - 1] + h__[i__]) * 2.;
    }

/*     compute the right hand side */

    temp1 = (fn[1] - fn[0]) / h__[0];
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp2 = (fn[i__ + 1] - fn[i__]) / h__[i__];
	rs[i__] = (temp2 - temp1) * 3.;
	temp1 = temp2;
    }
    rs[*n] = ((fn[1] - fn[0]) / h__[0] - temp1) * 3.;

/*     solve the linear system with factorization */

    cytsy_(n, np, dm, du, rc, rs, c__, &iflag);
    if (iflag != 1) {
	return 0;
    }

/*     compute remaining spline coefficients */

    c__[0] = c__[*n];
    i__1 = *n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	b[i__] = (fn[i__ + 1] - fn[i__]) / h__[i__] - h__[i__] / 3. * (c__[
		i__ + 1] + c__[i__] * 2.);
	d__[i__] = (c__[i__ + 1] - c__[i__]) / (h__[i__] * 3.);
    }
    b[*n] = (fn[1] - fn[*n]) / h__[*n] - h__[*n] / 3. * (c__[1] + c__[*n] * 
	    2.);
    return 0;
} /* cspline_ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine cytsy  --  solve cyclic tridiagonal system  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "cytsy" solves a system of linear equations for a cyclically */
/*     tridiagonal, symmetric, positive definite matrix */

/*     literature reference: */

/*     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran, */
/*     Springer Verlag, 1996, Section 4.11.2 */


/* Subroutine */ int cytsy_(integer *n, integer *np, doublereal *dm, 
	doublereal *du, doublereal *cr, doublereal *rs, doublereal *x, 
	integer *iflag)
{
    extern /* Subroutine */ int cytsyp_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *), cytsys_(integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *);



/*     factorization of the input matrix */

    *iflag = -2;
    if (*n < 3) {
	return 0;
    }
    cytsyp_(n, np, dm, du, cr, iflag);

/*     update and back substitute as necessary */

    if (*iflag == 1) {
	cytsys_(n, np, dm, du, cr, rs, x);
    }
    return 0;
} /* cytsy_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine cytsyp  --  tridiagonal Cholesky factorization  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "cytsyp" finds the Cholesky factors of a cyclically tridiagonal */
/*     symmetric, positive definite matrix given by two vectors */

/*     literature reference: */

/*     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran, */
/*     Springer Verlag, 1996, Section 4.11.2 */


/* Subroutine */ int cytsyp_(integer *n, integer *np, doublereal *dm, 
	doublereal *du, doublereal *cr, integer *iflag)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal d__;
    static integer i__;
    static doublereal eps, row, temp1, temp2;



/*     set error bound and test for condition n greater than 2 */

    eps = 1e-8;
    *iflag = -2;
    if (*n < 3) {
	return 0;
    }

/*     checking to see if matrix is positive definite */

    row = abs(dm[1]) + abs(du[1]) + (d__1 = du[*n], abs(d__1));
    if (row == 0.) {
	*iflag = 0;
	return 0;
    }
    d__ = 1. / row;
    if (dm[1] < 0.) {
	*iflag = -1;
	return 0;
    } else if (abs(dm[1]) * d__ <= eps) {
	*iflag = 0;
	return 0;
    }

/*     factoring a while checking for a positive definite and strong */
/*     nonsingular matrix a */

    temp1 = du[1];
    du[1] /= dm[1];
    cr[1] = du[*n] / dm[1];
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	row = (d__1 = dm[i__], abs(d__1)) + (d__2 = du[i__], abs(d__2)) + abs(
		temp1);
	if (row == 0.) {
	    *iflag = 0;
	    return 0;
	}
	d__ = 1. / row;
	dm[i__] -= temp1 * du[i__ - 1];
	if (dm[i__] < 0.) {
	    *iflag = -1;
	    return 0;
	} else if ((d__1 = dm[i__], abs(d__1)) * d__ <= eps) {
	    *iflag = 0;
	    return 0;
	}
	if (i__ < *n - 1) {
	    cr[i__] = -temp1 * cr[i__ - 1] / dm[i__];
	    temp1 = du[i__];
	    du[i__] /= dm[i__];
	} else {
	    temp2 = du[i__];
	    du[i__] = (du[i__] - temp1 * cr[i__ - 1]) / dm[i__];
	}
    }
    row = (d__1 = du[*n], abs(d__1)) + (d__2 = dm[*n], abs(d__2)) + abs(temp2)
	    ;
    if (row == 0.) {
	*iflag = 0;
	return 0;
    }
    d__ = 1. / row;
    dm[*n] -= dm[*n - 1] * du[*n - 1] * du[*n - 1];
    temp1 = 0.;
    i__1 = *n - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp1 += dm[i__] * cr[i__] * cr[i__];
    }
    dm[*n] -= temp1;
    if (dm[*n] < 0.) {
	*iflag = -1;
	return 0;
    } else if ((d__1 = dm[*n], abs(d__1)) * d__ <= eps) {
	*iflag = 0;
	return 0;
    }
    *iflag = 1;
    return 0;
} /* cytsyp_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine cytsys  --  tridiagonal solution from factors  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "cytsys" solves a cyclically tridiagonal linear system */
/*     given the Cholesky factors */

/*     literature reference: */

/*     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran, */
/*     Springer Verlag, 1996, Section 4.11.2 */


/* Subroutine */ int cytsys_(integer *n, integer *np, doublereal *dm, 
	doublereal *du, doublereal *cr, doublereal *rs, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal sum, temp;



/*     updating phase */

    temp = rs[1];
    rs[1] = temp / dm[1];
    sum = cr[1] * temp;
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	temp = rs[i__] - du[i__ - 1] * temp;
	rs[i__] = temp / dm[i__];
	if (i__ != *n - 1) {
	    sum += cr[i__] * temp;
	}
    }
    temp = rs[*n] - du[*n - 1] * temp;
    temp -= sum;
    rs[*n] = temp / dm[*n];

/*     backsubstitution phase */

    x[*n] = rs[*n];
    x[*n - 1] = rs[*n - 1] - du[*n - 1] * x[*n];
    for (i__ = *n - 2; i__ >= 1; --i__) {
	x[i__] = rs[i__] - du[i__] * x[i__ + 1] - cr[i__] * x[*n];
    }
    return 0;
} /* cytsys_ */

