/* promo.f -- translated by f2c (version 20050501).
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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine promo  --  copywrite notice and version info  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "promo" writes a short message containing information */
/*     about the TINKER version number and the copyright notice */


/* Subroutine */ int promo_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 \002,78(\002#\002),/,\002 \002,78(\002"
	    "#\002),/,\002 ##\002,74x,\002##\002,/,\002 ##\002,13x,\002TINKER"
	    "  ---  Software Tools for\002,\002 Molecular Design\002,13x,\002"
	    "##\002,/,\002 ##\002,74x,\002##\002,/,\002 ##\002,24x,\002Versio"
	    "n 5.1  February 2010\002,24x,\002##\002,/,\002 ##\002,74x,\002#"
	    "#\002,/,\002 ##\002,15x,\002Copyright (c)  Jay William Ponder"
	    "\002,\002  1990-2010\002,15x,\002##\002,/,\002 ##\002,28x,\002Al"
	    "l Rights Reserved\002,27x,\002##\002,/,\002 ##\002,74x,\002##"
	    "\002,/,\002 \002,78(\002#\002),/,\002 \002,78(\002#\002),/)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };




/*     print out the informational header message */



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


    io___1.ciunit = iounit_1.iout;
    s_wsfe(&io___1);
    e_wsfe();
    return 0;
} /* promo_ */

