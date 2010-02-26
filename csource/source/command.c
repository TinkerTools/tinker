/* command.f -- translated by f2c (version 20050501).
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

/* Table of constant values */

static integer c__3 = 3;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine command  --  get any command line arguments  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "command" uses the standard Unix-like iargc/getarg routines */
/*     to get the number and values of arguments specified on the */
/*     command line at program runtime */


/* Subroutine */ int command_(void)
{
    /* System generated locals */
    address a__1[3];
    integer i__1[3], i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i__;
    extern integer iargc_(void);
    static char blank[20];
    extern /* Subroutine */ int getarg_(integer *, char *, ftnlen), upcase_(
	    char *, ftnlen);
    static char letter[1];


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




/*     initialize command line arguments as blank strings */

    s_copy(blank, "                    ", (ftnlen)20, (ftnlen)20);
    for (i__ = 0; i__ <= 20; ++i__) {
/* Writing concatenation */
	i__1[0] = 20, a__1[0] = blank;
	i__1[1] = 20, a__1[1] = blank;
	i__1[2] = 20, a__1[2] = blank;
	s_cat(arg_ref(0, i__), a__1, i__1, &c__3, (ftnlen)120);
    }

/*     get the number of arguments and store each in a string */

    argue_1.narg = iargc_();
    if (argue_1.narg > 20) {
	argue_1.narg = 20;
    }
    i__2 = argue_1.narg;
    for (i__ = 0; i__ <= i__2; ++i__) {
	getarg_(&i__, arg_ref(0, i__), (ftnlen)120);
    }

/*     mark the command line options as unuseable for input */

    argue_1.listarg[0] = FALSE_;
    i__2 = argue_1.narg;
    for (i__ = 1; i__ <= i__2; ++i__) {
	argue_1.listarg[i__] = TRUE_;
    }
    i__2 = argue_1.narg;
    for (i__ = 1; i__ <= i__2; ++i__) {
	*(unsigned char *)letter = *(unsigned char *)arg_ref(0, i__);
	if (*(unsigned char *)letter == '-') {
	    *(unsigned char *)letter = *(unsigned char *)arg_ref(1, i__);
	    upcase_(letter, (ftnlen)1);
	    if (*(unsigned char *)letter >= 'A' && *(unsigned char *)letter <=
		     'Z') {
		argue_1.listarg[i__] = FALSE_;
		argue_1.listarg[i__ + 1] = FALSE_;
	    }
	}
    }
    return 0;
} /* command_ */

#undef arg_ref


