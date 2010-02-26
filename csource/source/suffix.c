/* suffix.f -- translated by f2c (version 20050501).
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

/* Table of constant values */

static integer c__3 = 3;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine suffix  --  test for default file extension  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "suffix" checks a filename for the presence of an */
/*     extension, and appends an extension if none is found */


/* Subroutine */ int suffix_(char *filename, char *extension, ftnlen 
	filename_len, ftnlen extension_len)
{
    /* System generated locals */
    address a__1[3];
    integer i__1, i__2[3];
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer i__, leng, last, lext;
    static logical exist;
    static char letter[1];



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




/*     get the length of the current filename */

    leng = trimtext_(filename, filename_len);
    lext = trimtext_(extension, extension_len);

/*     check for an extension on the current filename */

    last = leng;
    i__1 = leng;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)letter = *(unsigned char *)&filename[i__ - 1];
	if (*(unsigned char *)letter == '/') {
	    last = leng;
	}
/*        if (letter .eq. '\')  last = leng */
	if (*(unsigned char *)letter == 92) {
	    last = leng;
	}
	if (*(unsigned char *)letter == ']') {
	    last = leng;
	}
	if (*(unsigned char *)letter == ':') {
	    last = leng;
	}
	if (*(unsigned char *)letter == '~') {
	    last = leng;
	}
	if (*(unsigned char *)letter == '.') {
	    last = i__ - 1;
	}
    }
    if (last != leng) {
	return 0;
    }

/*     append extension if current name does not exist */

    exist = FALSE_;
    if (leng != 0) {
	ioin__1.inerr = 0;
	ioin__1.infilen = leng;
	ioin__1.infile = filename;
	ioin__1.inex = &exist;
	ioin__1.inopen = 0;
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
    if (! exist) {
/* Writing concatenation */
	i__2[0] = leng, a__1[0] = filename;
	i__2[1] = 1, a__1[1] = ".";
	i__2[2] = lext, a__1[2] = extension;
	s_cat(filename, a__1, i__2, &c__3, filename_len);
    }
    return 0;
} /* suffix_ */

