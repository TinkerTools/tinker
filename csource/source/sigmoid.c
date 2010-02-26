/* sigmoid.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  function sigmoid  --  general sigmoidal functional form  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "sigmoid" implements a normalized sigmoidal function on the */
/*     interval [0,1]; the curves connect (0,0) to (1,1) and have */
/*     a cooperativity controlled by beta, they approach a straight */
/*     line as beta -> 0 and get more nonlinear as beta increases */


doublereal sigmoid_(doublereal *beta, doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal expmin, expmax, expterm;



/*     compute the value of the normalized sigmoidal function */

    if (*beta == 0.) {
	ret_val = *x;
    } else {
	expmax = 1. / (exp(-(*beta)) + 1.);
	expmin = 1. / (exp(*beta) + 1.);
	expterm = 1. / (exp(*beta * (*x * 2. - 1.)) + 1.);
	ret_val = (expmax - expterm) / (expmax - expmin);
    }
    return ret_val;
} /* sigmoid_ */

