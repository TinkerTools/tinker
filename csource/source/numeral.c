/* numeral.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine numeral  --  convert number to text string  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "numeral" converts an input integer number into the */
/*     corresponding right- or left-justified text numeral */

/*     number  integer value of the number to be transformed */
/*     string  text string to be filled with corresponding numeral */
/*     size    on input, the minimal acceptable numeral length, if */
/*               zero then output will be right justified, if */
/*               nonzero then numeral is left-justified and padded */
/*               with leading zeros as necessary; upon output, the */
/*               number of non-blank characters in the numeral */


/* Subroutine */ int numeral_(integer *number, char *string, integer *size, 
	ftnlen string_len)
{
    /* Initialized data */

    static char digit[1*10] = "0" "1" "2" "3" "4" "5" "6" "7" "8" "9";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len(char *, ftnlen);

    /* Local variables */
    static logical negative;
    static integer thousand, i__, pos, ones, tens;
    static logical right;
    static integer multi, length, hundred, million, minsize, tenthou, hunthou;



/*     set justification and size bounds for numeral string */

    if (*size == 0) {
	right = TRUE_;
	*size = 1;
    } else {
	right = FALSE_;
    }
    minsize = *size;
    length = i_len(string, string_len);

/*     test the sign of the original number */

    if (*number >= 0) {
	negative = FALSE_;
    } else {
	negative = TRUE_;
	*number = -(*number);
    }

/*     use modulo arithmetic to find place-holding digits */

    million = *number / 1000000;
    multi = million * 1000000;
    hunthou = (*number - multi) / 100000;
    multi += hunthou * 100000;
    tenthou = (*number - multi) / 10000;
    multi += tenthou * 10000;
    thousand = (*number - multi) / 1000;
    multi += thousand * 1000;
    hundred = (*number - multi) / 100;
    multi += hundred * 100;
    tens = (*number - multi) / 10;
    multi += tens * 10;
    ones = *number - multi;

/*     find the correct length to be used for the numeral */

    if (million != 0) {
	*size = 7;
    } else if (hunthou != 0) {
	*size = 6;
    } else if (tenthou != 0) {
	*size = 5;
    } else if (thousand != 0) {
	*size = 4;
    } else if (hundred != 0) {
	*size = 3;
    } else if (tens != 0) {
	*size = 2;
    } else {
	*size = 1;
    }
    *size = min(*size,length);
    *size = max(*size,minsize);

/*     convert individual digits to a string of numerals */

    if (*size == 7) {
	*(unsigned char *)string = *(unsigned char *)&digit[million];
	*(unsigned char *)&string[1] = *(unsigned char *)&digit[hunthou];
	*(unsigned char *)&string[2] = *(unsigned char *)&digit[tenthou];
	*(unsigned char *)&string[3] = *(unsigned char *)&digit[thousand];
	*(unsigned char *)&string[4] = *(unsigned char *)&digit[hundred];
	*(unsigned char *)&string[5] = *(unsigned char *)&digit[tens];
	*(unsigned char *)&string[6] = *(unsigned char *)&digit[ones];
    } else if (*size == 6) {
	*(unsigned char *)string = *(unsigned char *)&digit[hunthou];
	*(unsigned char *)&string[1] = *(unsigned char *)&digit[tenthou];
	*(unsigned char *)&string[2] = *(unsigned char *)&digit[thousand];
	*(unsigned char *)&string[3] = *(unsigned char *)&digit[hundred];
	*(unsigned char *)&string[4] = *(unsigned char *)&digit[tens];
	*(unsigned char *)&string[5] = *(unsigned char *)&digit[ones];
    } else if (*size == 5) {
	*(unsigned char *)string = *(unsigned char *)&digit[tenthou];
	*(unsigned char *)&string[1] = *(unsigned char *)&digit[thousand];
	*(unsigned char *)&string[2] = *(unsigned char *)&digit[hundred];
	*(unsigned char *)&string[3] = *(unsigned char *)&digit[tens];
	*(unsigned char *)&string[4] = *(unsigned char *)&digit[ones];
    } else if (*size == 4) {
	*(unsigned char *)string = *(unsigned char *)&digit[thousand];
	*(unsigned char *)&string[1] = *(unsigned char *)&digit[hundred];
	*(unsigned char *)&string[2] = *(unsigned char *)&digit[tens];
	*(unsigned char *)&string[3] = *(unsigned char *)&digit[ones];
    } else if (*size == 3) {
	*(unsigned char *)string = *(unsigned char *)&digit[hundred];
	*(unsigned char *)&string[1] = *(unsigned char *)&digit[tens];
	*(unsigned char *)&string[2] = *(unsigned char *)&digit[ones];
    } else if (*size == 2) {
	*(unsigned char *)string = *(unsigned char *)&digit[tens];
	*(unsigned char *)&string[1] = *(unsigned char *)&digit[ones];
    } else {
	*(unsigned char *)string = *(unsigned char *)&digit[ones];
    }

/*     right-justify if desired, and pad with blanks */

    if (right) {
	for (i__ = *size; i__ >= 1; --i__) {
	    pos = length - *size + i__;
	    *(unsigned char *)&string[pos - 1] = *(unsigned char *)&string[
		    i__ - 1];
	}
	i__1 = length - *size;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    *(unsigned char *)&string[i__ - 1] = ' ';
	}
    } else {
	i__1 = length;
	for (i__ = *size + 1; i__ <= i__1; ++i__) {
	    *(unsigned char *)&string[i__ - 1] = ' ';
	}
    }
    return 0;
} /* numeral_ */

