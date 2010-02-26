/* trimtext.f -- translated by f2c (version 20050501).
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
/*     ##  function trimtext  --  find last non-blank character  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "trimtext" finds and returns the location of the last */
/*     non-blank character before the first null character in */
/*     an input text string; the function returns zero if no */
/*     such character is found */


integer trimtext_(char *string, ftnlen string_len)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen);

    /* Local variables */
    static integer i__, last, size;
    static char null[1];



/*     move forward through the string, one character */
/*     at a time, looking for first null character */

    ret_val = 0;
    size = i_len(string, string_len);
    *(unsigned char *)null = '\0';
    last = size;
    i__1 = size;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*(unsigned char *)&string[i__ - 1] == *(unsigned char *)null) {
	    last = i__ - 1;
	    goto L10;
	}
    }
L10:

/*     move backward through the string, one character */
/*     at a time, looking for first non-blank character */

    for (i__ = last; i__ >= 1; --i__) {
	if (*(unsigned char *)&string[i__ - 1] > ' ') {
	    ret_val = i__;
	    goto L20;
	}
    }
L20:
    return ret_val;
} /* trimtext_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine justify  --  convert string to right justified  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "justify" converts a text string to right justified format */
/*     with leading blank spaces */


/* Subroutine */ int justify_(char *string, ftnlen string_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen);

    /* Local variables */
    static integer i__, k, last, size;
    static char null[1], letter[1];



/*     move backward through the string, one character */
/*     at a time, looking for first non-blank character */

    size = i_len(string, string_len);
    *(unsigned char *)null = '\0';
    last = 0;
    for (i__ = size; i__ >= 1; --i__) {
	*(unsigned char *)letter = *(unsigned char *)&string[i__ - 1];
	if (*(unsigned char *)letter != ' ' && *(unsigned char *)letter != *(
		unsigned char *)null) {
	    last = i__;
	    goto L10;
	}
    }
L10:

/*     move string to the right and pad with leading blanks */

    for (i__ = last; i__ >= 1; --i__) {
	k = i__ + size - last;
	*(unsigned char *)&string[k - 1] = *(unsigned char *)&string[i__ - 1];
    }
    i__1 = size - last;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)&string[i__ - 1] = ' ';
    }
    return 0;
} /* justify_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine upcase  --  convert string to all upper case  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "upcase" converts a text string to all upper case letters */


/* Subroutine */ int upcase_(char *string, ftnlen string_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen);

    /* Local variables */
    static integer i__, code, size;
    static char letter[1];



/*     convert lower case to upper case one letter at a time */

    size = i_len(string, string_len);
    i__1 = size;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)letter = *(unsigned char *)&string[i__ - 1];
	code = *(unsigned char *)letter;
	if (*(unsigned char *)letter >= 'a' && *(unsigned char *)letter <= 
		'z') {
	    *(unsigned char *)&string[i__ - 1] = (char) (code - 32);
	}
    }
    return 0;
} /* upcase_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine lowcase  --  convert string to all lower case  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "lowcase" converts a text string to all lower case letters */


/* Subroutine */ int lowcase_(char *string, ftnlen string_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen);

    /* Local variables */
    static integer i__, code, size;
    static char letter[1];



/*     convert upper case to lower case one letter at a time */

    size = i_len(string, string_len);
    i__1 = size;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)letter = *(unsigned char *)&string[i__ - 1];
	code = *(unsigned char *)letter;
	if (*(unsigned char *)letter >= 'A' && *(unsigned char *)letter <= 
		'Z') {
	    *(unsigned char *)&string[i__ - 1] = (char) (code + 32);
	}
    }
    return 0;
} /* lowcase_ */

