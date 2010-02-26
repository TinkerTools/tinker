/* getstring.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine getstring  --  extract single quoted string  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "getstring" searchs for a quoted text string within an input */
/*     character string; the region between the first and second */
/*     quotes is returned as the "text"; if the actual text is too */
/*     long, only the first part is returned */

/*     variables and parameters: */

/*     string    input character string to be searched */
/*     text      the quoted text found in the input string */
/*     next      input with first position of search string; */
/*                 output with the position following text */


/* Subroutine */ int getstring_(char *string, char *text, integer *next, 
	ftnlen string_len, ftnlen text_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer i_len(char *, ftnlen);

    /* Local variables */
    static integer i__, j, code, last, size, final, first, length, extent, 
	    initial;



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




/*     get the length of input string and output text */

    length = i_len(string + (*next - 1), string_len - (*next - 1));
    size = i_len(text, text_len);

/*     move through the string one character at a time, */
/*     searching for the quoted text string characters */

    first = *next;
    last = 0;
    initial = *next;
    final = *next + length - 1;
    i__1 = final;
    for (i__ = initial; i__ <= i__1; ++i__) {
	code = *(unsigned char *)&string[i__ - 1];
	if (code == 34) {
	    first = i__ + 1;
	    i__2 = final;
	    for (j = first; j <= i__2; ++j) {
		code = *(unsigned char *)&string[j - 1];
		if (code == 34) {
		    last = j - 1;
		    *next = j + 1;
		    goto L10;
		}
	    }
	}
    }
L10:

/*     trim the actual word if it is too long to return */

    extent = last - first + 1;
    final = first + size - 1;
    if (extent > size) {
	last = final;
    }

/*     transfer the text into the return string */

    j = 0;
    i__1 = last;
    for (i__ = first; i__ <= i__1; ++i__) {
	++j;
	*(unsigned char *)&text[j - 1] = *(unsigned char *)&string[i__ - 1];
    }
    i__1 = final;
    for (i__ = last + 1; i__ <= i__1; ++i__) {
	++j;
	*(unsigned char *)&text[j - 1] = ' ';
    }
    return 0;
} /* getstring_ */

