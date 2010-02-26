/* freeunit.f -- translated by f2c (version 20050501).
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
/*     ##  function freeunit  --  gets an unopened logical unit  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "freeunit" finds an unopened Fortran I/O unit and returns */
/*     its numerical value from 1 to 99; the units already assigned */
/*     to "input" and "iout" (usually 5 and 6) are skipped since */
/*     they have special meaning as the default I/O units */


integer freeunit_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 FREEUNIT  --  No Available Fortran\002"
	    ",\002 I/O Units\002)";

    /* System generated locals */
    integer ret_val;
    inlist ioin__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), f_inqu(inlist *);

    /* Local variables */
    static logical used;
    extern /* Subroutine */ int fatal_(void);

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




/*     try each logical unit until an unopened one is found */

    ret_val = 0;
    used = TRUE_;
    while(used) {
	++ret_val;
	if (ret_val != iounit_1.input && ret_val != iounit_1.iout) {
	    if (ret_val > 99) {
		io___2.ciunit = iounit_1.iout;
		s_wsfe(&io___2);
		e_wsfe();
		fatal_();
	    }
	    ioin__1.inerr = 0;
	    ioin__1.inunit = ret_val;
	    ioin__1.infile = 0;
	    ioin__1.inex = 0;
	    ioin__1.inopen = &used;
	    ioin__1.innum = 0;
	    ioin__1.innamed = 0;
	    ioin__1.inname = 0;
	    ioin__1.inacc = 0;
	    ioin__1.inseq = 0;
	    ioin__1.indir = 0;
	    ioin__1.infmt = 0;
	    ioin__1.inform = 0;
	    ioin__1.inunf = 0;
	    ioin__1.inrecl = 0;
	    ioin__1.innrec = 0;
	    ioin__1.inblank = 0;
	    f_inqu(&ioin__1);
	}
    }
    return ret_val;
} /* freeunit_ */

