/* nexttext.f -- translated by f2c (version 20050501).
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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  function nexttext  --  find next non-blank character  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "nexttext" finds and returns the location of the first */
/*     non-blank character within an input text string; zero */
/*     is returned if no such character is found */


integer nexttext_(char *string, ftnlen string_len)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen);

    /* Local variables */
    static integer i__, size;



/*     move forward through the string, one character */
/*     at a time, looking for first non-blank character */

    ret_val = 0;
    size = i_len(string, string_len);
    i__1 = size;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*(unsigned char *)&string[i__ - 1] > ' ') {
	    ret_val = i__;
	    goto L10;
	}
    }
L10:
    return ret_val;
} /* nexttext_ */

