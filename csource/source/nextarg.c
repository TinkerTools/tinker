/* nextarg.f -- translated by f2c (version 20050501).
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
    integer narg;
    logical listarg[21];
    char arg[2520];
} argue_;

#define argue_1 argue_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine nextarg  --  find next command line argument  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "nextarg" finds the next unused command line argument */
/*     and returns it in the input character string */


/* Subroutine */ int nextarg_(char *string, logical *exist, ftnlen string_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_len(char *, ftnlen);

    /* Local variables */
    static integer i__, length;


#define arg_ref(a_0,a_1) &argue_1.arg[(a_1)*120 + a_0 - 0]



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  argue.i  --  command line arguments at program startup  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     maxarg    maximum number of command line arguments */

/*     narg      number of command line arguments to the program */
/*     listarg   flag to mark available command line arguments */
/*     arg       strings containing the command line arguments */




/*     initialize the command argument as a blank string */

    s_copy(string, "          ", string_len, (ftnlen)10);
    *exist = FALSE_;

/*     get the next command line argument and mark it as used */

    if (argue_1.narg != 0) {
/* Computing MIN */
	i__1 = i_len(string, string_len), i__2 = i_len(arg_ref(0, 20), (
		ftnlen)120);
	length = min(i__1,i__2);
	i__1 = argue_1.narg;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (argue_1.listarg[i__]) {
		argue_1.listarg[i__] = FALSE_;
		s_copy(string, arg_ref(0, i__), string_len, length);
		*exist = TRUE_;
		goto L10;
	    }
	}
L10:
	;
    }
    return 0;
} /* nextarg_ */

#undef arg_ref


