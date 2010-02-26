/* number.f -- translated by f2c (version 20050501).
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
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  function number  --  convert text string to number  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     "number" converts a text numeral into an integer value; */
/*     the input string must contain only numeric characters */


integer number_(char *string, ftnlen string_len)
{
    /* Initialized data */

    static integer place[10] = { 1,10,100,1000,10000,100000,1000000,10000000,
	    100000000,1000000000 };

    /* Format strings */
    static char fmt_10[] = "(\002 NUMBER  --  Input Text String is Too Lon"
	    "g\002)";
    static char fmt_30[] = "(/,\002 NUMBER  --  Non-Numeric Characters Foun"
	    "d\002,\002 in Numeral String\002)";

    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j, last, digit, first;
    static char letter[1];

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_30, 0 };




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  inform.i  --  control values for I/O and program flow  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     digits    decimal places output for energy and coordinates */
/*     iprint    steps between status printing (0=no printing) */
/*     iwrite    steps between coordinate dumps (0=no dumps) */
/*     isend     steps between socket communication (0=no sockets) */
/*     verbose   logical flag to turn on extra information */
/*     debug     logical flag to turn on full debug printing */
/*     holdup    logical flag to wait for carriage return on exit */
/*     abort     logical flag to stop execution at next chance */




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




/*     initialize the integer value of number to zero */

    ret_val = 0;

/*     get the first and last nonblank characters */

    last = trimtext_(string, string_len);
    if (last > 10) {
	io___3.ciunit = iounit_1.iout;
	s_wsfe(&io___3);
	e_wsfe();
	return ret_val;
    }
    first = 1;
    i__1 = last;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)letter = *(unsigned char *)&string[i__ - 1];
	if (*(unsigned char *)letter != ' ') {
	    first = i__;
	    goto L20;
	}
    }
L20:

/*     convert the text numeral into an integer number */

    j = 0;
    i__1 = first;
    for (i__ = last; i__ >= i__1; --i__) {
	++j;
	*(unsigned char *)letter = *(unsigned char *)&string[i__ - 1];
	if (*(unsigned char *)letter == '0') {
	    digit = 0;
	} else if (*(unsigned char *)letter == '1') {
	    digit = 1;
	} else if (*(unsigned char *)letter == '2') {
	    digit = 2;
	} else if (*(unsigned char *)letter == '3') {
	    digit = 3;
	} else if (*(unsigned char *)letter == '4') {
	    digit = 4;
	} else if (*(unsigned char *)letter == '5') {
	    digit = 5;
	} else if (*(unsigned char *)letter == '6') {
	    digit = 6;
	} else if (*(unsigned char *)letter == '7') {
	    digit = 7;
	} else if (*(unsigned char *)letter == '8') {
	    digit = 8;
	} else if (*(unsigned char *)letter == '9') {
	    digit = 9;
	} else {
	    if (inform_1.debug) {
		io___9.ciunit = iounit_1.iout;
		s_wsfe(&io___9);
		e_wsfe();
	    }
	    ret_val = 0;
	    goto L40;
	}
	ret_val += digit * place[j - 1];
    }
L40:
    return ret_val;
} /* number_ */

