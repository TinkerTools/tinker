/* invert.f -- translated by f2c (version 20050501).
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
/*     ##  subroutine invert  --  gauss-jordan matrix inversion  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "invert" inverts a matrix using the Gauss-Jordan method */

/*     variables and parameters: */

/*     n     logical dimension of the matrix to be inverted */
/*     np    physical dimension of the matrix storage area */
/*     a     matrix to invert; contains inverse on exit */


/* Subroutine */ int invert_(integer *n, integer *np, doublereal *a)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 INVERT  --  Matrix Too Large; Increase M"
	    "AXINV\002)";
    static char fmt_20[] = "(/,\002 INVERT  --  Cannot Invert\002,\002 a Sin"
	    "gular Matrix\002)";
    static char fmt_30[] = "(/,\002 INVERT  --  Cannot Invert a Singular Mat"
	    "rix\002)";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k;
    static doublereal big;
    static integer icol;
    static doublereal temp;
    static integer irow;
    extern /* Subroutine */ int fatal_(void);
    static integer indxc[100], indxr[100];
    static doublereal pivot;
    static integer ipivot[100];

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_30, 0 };



#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]



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




/*     check to see if the matrix is too large to handle */

    /* Parameter adjustments */
    a_dim1 = *np;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    if (*n > 100) {
	io___1.ciunit = iounit_1.iout;
	s_wsfe(&io___1);
	e_wsfe();
	fatal_();
    }

/*     perform matrix inversion via the Gauss-Jordan algorithm */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ipivot[i__ - 1] = 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	big = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    if (ipivot[j - 1] != 1) {
		i__3 = *n;
		for (k = 1; k <= i__3; ++k) {
		    if (ipivot[k - 1] == 0) {
			if ((d__1 = a_ref(j, k), abs(d__1)) >= big) {
			    big = (d__1 = a_ref(j, k), abs(d__1));
			    irow = j;
			    icol = k;
			}
		    } else if (ipivot[k - 1] > 1) {
			io___9.ciunit = iounit_1.iout;
			s_wsfe(&io___9);
			e_wsfe();
			fatal_();
		    }
		}
	    }
	}
	++ipivot[icol - 1];
	if (irow != icol) {
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		temp = a_ref(irow, j);
		a_ref(irow, j) = a_ref(icol, j);
		a_ref(icol, j) = temp;
	    }
	}
	indxr[i__ - 1] = irow;
	indxc[i__ - 1] = icol;
	if (a_ref(icol, icol) == 0.) {
	    io___13.ciunit = iounit_1.iout;
	    s_wsfe(&io___13);
	    e_wsfe();
	    fatal_();
	}
	pivot = a_ref(icol, icol);
	a_ref(icol, icol) = 1.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    a_ref(icol, j) = a_ref(icol, j) / pivot;
	}
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    if (j != icol) {
		temp = a_ref(j, icol);
		a_ref(j, icol) = 0.;
		i__3 = *n;
		for (k = 1; k <= i__3; ++k) {
		    a_ref(j, k) = a_ref(j, k) - a_ref(icol, k) * temp;
		}
	    }
	}
    }
    for (i__ = *n; i__ >= 1; --i__) {
	if (indxr[i__ - 1] != indxc[i__ - 1]) {
	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {
		temp = a_ref(k, indxr[i__ - 1]);
		a_ref(k, indxr[i__ - 1]) = a_ref(k, indxc[i__ - 1]);
		a_ref(k, indxc[i__ - 1]) = temp;
	    }
	}
    }
    return 0;
} /* invert_ */

#undef a_ref


