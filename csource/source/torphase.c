/* torphase.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine torphase  --  torsional amplitude and phase  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "torphase" sets the n-fold amplitude and phase values */
/*     for each torsion via sorting of the input parameters */


/* Subroutine */ int torphase_(integer *ft, doublereal *vt, doublereal *st)
{
    static integer i__, k;
    static doublereal phase[6], ampli[6];



/*     copy the input fold, amplitude and phase angles */

    /* Parameter adjustments */
    --st;
    --vt;
    --ft;

    /* Function Body */
    for (i__ = 1; i__ <= 6; ++i__) {
	ampli[i__ - 1] = vt[i__];
	phase[i__ - 1] = st[i__];
	vt[i__] = 0.;
	st[i__] = 0.;
    }

/*     shift the phase angles into the standard range */

    for (i__ = 1; i__ <= 6; ++i__) {
	while(phase[i__ - 1] < -180.) {
	    phase[i__ - 1] += 360.;
	}
	while(phase[i__ - 1] > 180.) {
	    phase[i__ - 1] += -360.;
	}
    }

/*     convert input torsional parameters to storage format */

    for (i__ = 1; i__ <= 6; ++i__) {
	k = ft[i__];
	if (k == 0) {
	    goto L10;
	} else if (k <= 6) {
	    vt[k] = ampli[i__ - 1];
	    st[k] = phase[i__ - 1];
	}
    }
L10:
    return 0;
} /* torphase_ */

