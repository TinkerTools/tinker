/* getnumb.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine getnumb  --  extract an integer from a string  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "getnumb" searchs an input string from left to right for an */
/*     integer and puts the numeric value in "number"; returns zero */
/*     with "next" unchanged if no integer value is found */

/*     variables and parameters: */

/*     string    input character string to be searched */
/*     number    output with the first integer in the string */
/*     next      input with first position of search string; */
/*                 output with the position following the number */


/* Subroutine */ int getnumb_(char *string, integer *number, integer *next, 
	ftnlen string_len)
{
    /* Initialized data */

    static integer place[10] = { 1,10,100,1000,10000,100000,1000000,10000000,
	    100000000,1000000000 };

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j, code, last, final, digit, first;
    static logical negate;
    static integer length;
    static char letter[1];
    static integer initial;
    static logical numeral;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  2004  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  ascii.i  --  selected values of ASCII character codes  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     null         decimal value of ASCII code for null (0) */
/*     tab          decimal value of ASCII code for tab (9) */
/*     linefeed     decimal value of ASCII code for linefeed (10) */
/*     formfeed     decimal value of ASCII code for formfeed (12) */
/*     carriage     decimal value of ASCII code for carriage return (13) */
/*     escape       decimal value of ASCII code for escape (27) */
/*     space        decimal value of ASCII code for blank space (32) */
/*     exclamation  decimal value of ASCII code for exclamation (33) */
/*     quote        decimal value of ASCII code for double quote (34) */
/*     pound        decimal value of ASCII code for pound sign (35) */
/*     dollar       decimal value of ASCII code for dollar sign (36) */
/*     percent      decimal value of ASCII code for percent sign (37) */
/*     ampersand    decimal value of ASCII code for ampersand (38) */
/*     apostrophe   decimal value of ASCII code for single quote (39) */
/*     asterisk     decimal value of ASCII code for asterisk (42) */
/*     plus         decimal value of ASCII code for plus sign (43) */
/*     comma        decimal value of ASCII code for comma (44) */
/*     minus        decimal value of ASCII code for minus sign (45) */
/*     period       decimal value of ASCII code for period (46) */
/*     frontslash   decimal value of ASCII codd for frontslash (47) */
/*     colon        decimal value of ASCII code for colon (58) */
/*     semicolon    decimal value of ASCII code for semicolon (59) */
/*     equal        decimal value of ASCII code for equal sign (61) */
/*     question     decimal value of ASCII code for question mark (63) */
/*     atsign       decimal value of ASCII code for at sign (64) */
/*     backslash    decimal value of ASCII code for backslash (92) */
/*     caret        decimal value of ASCII code for caret (94) */
/*     underbar     decimal value of ASCII code for underbar (95) */
/*     vertical     decimal value of ASCII code for vertical bar (124) */
/*     tilde        decimal value of ASCII code for tilde (126) */




/*     initialize number and get the input text string length */

    *number = 0;
    negate = FALSE_;
    numeral = FALSE_;
    length = trimtext_(string + (*next - 1), string_len - (*next - 1));

/*     move through the string one character at a time, */
/*     searching for the first run of numeric characters */

    first = *next;
    last = 0;
    initial = *next;
    final = *next + length - 1;
    i__1 = final;
    for (i__ = initial; i__ <= i__1; ++i__) {
	*(unsigned char *)letter = *(unsigned char *)&string[i__ - 1];
	code = *(unsigned char *)letter;
	if (*(unsigned char *)letter >= '0' && *(unsigned char *)letter <= 
		'9') {
	    if (! numeral) {
		numeral = TRUE_;
		first = i__;
	    }
	    if (i__ == final) {
		last = final;
		*next = i__ + 1;
	    }
	} else if (code == 45 && ! negate) {
	    negate = TRUE_;
	} else if (numeral) {
	    if (code == 32 || code == 9 || code == 44 || code == 59 || code ==
		     58 || code == 95) {
		last = i__ - 1;
		*next = i__;
	    } else {
		numeral = FALSE_;
	    }
	    goto L10;
	} else if (negate) {
	    numeral = FALSE_;
	    goto L10;
	} else if (code != 32 && code != 9) {
	    numeral = FALSE_;
	    goto L10;
	}
    }
L10:

/*     trim the actual number if it is too big to return */

    if (! numeral) {
	*next = initial;
    }
/* Computing MIN */
    i__1 = last, i__2 = first + 9;
    last = min(i__1,i__2);

/*     convert the text numeral into an integer number */

    j = 0;
    i__1 = first;
    for (i__ = last; i__ >= i__1; --i__) {
	++j;
	if (*(unsigned char *)&string[i__ - 1] == '0') {
	    digit = 0;
	} else if (*(unsigned char *)&string[i__ - 1] == '1') {
	    digit = 1;
	} else if (*(unsigned char *)&string[i__ - 1] == '2') {
	    digit = 2;
	} else if (*(unsigned char *)&string[i__ - 1] == '3') {
	    digit = 3;
	} else if (*(unsigned char *)&string[i__ - 1] == '4') {
	    digit = 4;
	} else if (*(unsigned char *)&string[i__ - 1] == '5') {
	    digit = 5;
	} else if (*(unsigned char *)&string[i__ - 1] == '6') {
	    digit = 6;
	} else if (*(unsigned char *)&string[i__ - 1] == '7') {
	    digit = 7;
	} else if (*(unsigned char *)&string[i__ - 1] == '8') {
	    digit = 8;
	} else if (*(unsigned char *)&string[i__ - 1] == '9') {
	    digit = 9;
	}
	*number += digit * place[j - 1];
    }
    if (negate) {
	*number = -(*number);
    }
    return 0;
} /* getnumb_ */

