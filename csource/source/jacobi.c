/* jacobi.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine jacobi  --  jacobi matrix diagonalization  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "jacobi" performs a matrix diagonalization of a real */
/*     symmetric matrix by the method of Jacobi rotations */

/*     variables and parameters: */

/*     n     logical dimension of the matrix to be diagonalized */
/*     np    physical dimension of the matrix storage area */
/*     a     input with the matrix to be diagonalized; only */
/*              the upper triangle and diagonal are required */
/*     d     returned with the eigenvalues in ascending order */
/*     v     returned with the eigenvectors of the matrix */
/*     b     temporary work vector */
/*     z     temporary work vector */


/* Subroutine */ int jacobi_(integer *n, integer *np, doublereal *a, 
	doublereal *d__, doublereal *v, doublereal *b, doublereal *z__)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 JACOBI  --  Matrix Diagonalization not C"
	    "onverged\002)";

    /* System generated locals */
    integer a_dim1, a_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static doublereal c__, g, h__;
    static integer i__, j, k;
    static doublereal p, s, t;
    static integer ip, iq;
    static doublereal sm, tau;
    static integer nrot;
    static doublereal theta, tresh;
    static integer maxrot;

    /* Fortran I/O blocks */
    static cilist io___16 = { 0, 0, 0, fmt_20, 0 };



#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define v_ref(a_1,a_2) v[(a_2)*v_dim1 + a_1]



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




/*     setup and initialization */

    /* Parameter adjustments */
    --z__;
    --b;
    v_dim1 = *np;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --d__;
    a_dim1 = *np;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    maxrot = 100;
    nrot = 0;
    i__1 = *n;
    for (ip = 1; ip <= i__1; ++ip) {
	i__2 = *n;
	for (iq = 1; iq <= i__2; ++iq) {
	    v_ref(ip, iq) = 0.;
	}
	v_ref(ip, ip) = 1.;
    }
    i__1 = *n;
    for (ip = 1; ip <= i__1; ++ip) {
	b[ip] = a_ref(ip, ip);
	d__[ip] = b[ip];
	z__[ip] = 0.;
    }

/*     perform the jacobi rotations */

    i__1 = maxrot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sm = 0.;
	i__2 = *n - 1;
	for (ip = 1; ip <= i__2; ++ip) {
	    i__3 = *n;
	    for (iq = ip + 1; iq <= i__3; ++iq) {
		sm += (d__1 = a_ref(ip, iq), abs(d__1));
	    }
	}
	if (sm == 0.) {
	    goto L10;
	}
	if (i__ < 4) {
/* Computing 2nd power */
	    i__2 = *n;
	    tresh = sm * .2 / (i__2 * i__2);
	} else {
	    tresh = 0.;
	}
	i__2 = *n - 1;
	for (ip = 1; ip <= i__2; ++ip) {
	    i__3 = *n;
	    for (iq = ip + 1; iq <= i__3; ++iq) {
		g = (d__1 = a_ref(ip, iq), abs(d__1)) * 100.;
		if (i__ > 4 && (d__1 = d__[ip], abs(d__1)) + g == (d__2 = d__[
			ip], abs(d__2)) && (d__3 = d__[iq], abs(d__3)) + g == 
			(d__4 = d__[iq], abs(d__4))) {
		    a_ref(ip, iq) = 0.;
		} else if ((d__1 = a_ref(ip, iq), abs(d__1)) > tresh) {
		    h__ = d__[iq] - d__[ip];
		    if (abs(h__) + g == abs(h__)) {
			t = a_ref(ip, iq) / h__;
		    } else {
			theta = h__ * .5 / a_ref(ip, iq);
/* Computing 2nd power */
			d__1 = theta;
			t = 1. / (abs(theta) + sqrt(d__1 * d__1 + 1.));
			if (theta < 0.) {
			    t = -t;
			}
		    }
/* Computing 2nd power */
		    d__1 = t;
		    c__ = 1. / sqrt(d__1 * d__1 + 1.);
		    s = t * c__;
		    tau = s / (c__ + 1.);
		    h__ = t * a_ref(ip, iq);
		    z__[ip] -= h__;
		    z__[iq] += h__;
		    d__[ip] -= h__;
		    d__[iq] += h__;
		    a_ref(ip, iq) = 0.;
		    i__4 = ip - 1;
		    for (j = 1; j <= i__4; ++j) {
			g = a_ref(j, ip);
			h__ = a_ref(j, iq);
			a_ref(j, ip) = g - s * (h__ + g * tau);
			a_ref(j, iq) = h__ + s * (g - h__ * tau);
		    }
		    i__4 = iq - 1;
		    for (j = ip + 1; j <= i__4; ++j) {
			g = a_ref(ip, j);
			h__ = a_ref(j, iq);
			a_ref(ip, j) = g - s * (h__ + g * tau);
			a_ref(j, iq) = h__ + s * (g - h__ * tau);
		    }
		    i__4 = *n;
		    for (j = iq + 1; j <= i__4; ++j) {
			g = a_ref(ip, j);
			h__ = a_ref(iq, j);
			a_ref(ip, j) = g - s * (h__ + g * tau);
			a_ref(iq, j) = h__ + s * (g - h__ * tau);
		    }
		    i__4 = *n;
		    for (j = 1; j <= i__4; ++j) {
			g = v_ref(j, ip);
			h__ = v_ref(j, iq);
			v_ref(j, ip) = g - s * (h__ + g * tau);
			v_ref(j, iq) = h__ + s * (g - h__ * tau);
		    }
		    ++nrot;
		}
	    }
	}
	i__2 = *n;
	for (ip = 1; ip <= i__2; ++ip) {
	    b[ip] += z__[ip];
	    d__[ip] = b[ip];
	    z__[ip] = 0.;
	}
    }

/*     print warning if not converged */

L10:
    if (nrot == maxrot) {
	io___16.ciunit = iounit_1.iout;
	s_wsfe(&io___16);
	e_wsfe();
    }

/*     sort the eigenvalues and vectors */

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = i__;
	p = d__[i__];
	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    if (d__[j] < p) {
		k = j;
		p = d__[j];
	    }
	}
	if (k != i__) {
	    d__[k] = d__[i__];
	    d__[i__] = p;
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		p = v_ref(j, i__);
		v_ref(j, i__) = v_ref(j, k);
		v_ref(j, k) = p;
	    }
	}
    }
    return 0;
} /* jacobi_ */

#undef v_ref
#undef a_ref


