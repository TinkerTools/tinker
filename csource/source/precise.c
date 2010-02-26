/* precise.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ######################################################### */
/*     ##                                                     ## */
/*     ##  function precise  --  determine machine precision  ## */
/*     ##                                                     ## */
/*     ######################################################### */


/*     "precise" finds a machine precision value as selected by */
/*     the input argument: (1) the smallest positive floating */
/*     point value, (2) the smallest relative floating point */
/*     spacing, (3) the largest relative floating point spacing */


doublereal precise_(integer *i__)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal one, zero, delta, value;



/*     set values for zero, one and multiplicative factor */

    zero = 0.;
    one = 1.;
    delta = 1.1;
    ret_val = one;

/*     find the smallest positive floating point value; */
/*     minimum of 0.24x10-307 is a patch for some SGI's, */
/*     for Sparc cpu's under Linux, etc. */

    if (*i__ == 1) {
/*        dowhile (precise .ne. zero) */
	while(ret_val >= 2.4e-308) {
	    value = ret_val;
	    ret_val /= delta;
	}
	ret_val = value;

/*     find the smallest relative floating point spacing */

    } else if (*i__ == 2) {
	while(one + ret_val != one) {
	    value = ret_val;
	    ret_val /= delta;
	}
	ret_val = value;

/*     find the largest relative floating point spacing */

    } else if (*i__ == 3) {
	while(one + ret_val != ret_val) {
	    value = ret_val;
	    ret_val *= delta;
	}
	ret_val = value;
    }
    return ret_val;
} /* precise_ */

