/* basefile.f -- translated by f2c (version 20050501).
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
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine basefile  --  get base prefix from a filename  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "basefile" extracts from an input filename the portion */
/*     consisting of any directory name and the base filename */


/* Subroutine */ int basefile_(char *string, ftnlen string_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer i__, k;
    extern /* Subroutine */ int getkey_(void);
    static char letter[1];
    extern /* Subroutine */ int control_(void);



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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  files.i  --  name and number of current structure files  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     nprior     number of previously existing cycle files */
/*     ldir       length in characters of the directory name */
/*     leng       length in characters of the base filename */
/*     filename   base filename used by default for all files */
/*     outfile    output filename used for intermediate results */




/*     store the input filename and find its full length */

    s_copy(files_1.filename, string, (ftnlen)120, (ftnlen)120);
    files_1.leng = trimtext_(string, (ftnlen)120);

/*     count the number of characters prior to any extension */

    k = files_1.leng;
    i__1 = files_1.leng;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)letter = *(unsigned char *)&string[i__ - 1];
	if (*(unsigned char *)letter == '/') {
	    k = files_1.leng;
	}
/*        if (letter .eq. '\')  k = leng */
	if (*(unsigned char *)letter == 92) {
	    k = files_1.leng;
	}
	if (*(unsigned char *)letter == ']') {
	    k = files_1.leng;
	}
	if (*(unsigned char *)letter == ':') {
	    k = files_1.leng;
	}
	if (*(unsigned char *)letter == '~') {
	    k = files_1.leng;
	}
	if (*(unsigned char *)letter == '.') {
	    k = i__ - 1;
	}
    }
    files_1.leng = min(files_1.leng,k);

/*     find the length of any directory name prefix */

    k = 0;
    for (i__ = files_1.leng; i__ >= 1; --i__) {
	*(unsigned char *)letter = *(unsigned char *)&string[i__ - 1];
	if (*(unsigned char *)letter == '/') {
	    k = i__;
	}
/*        if (letter .eq. '\')  k = i */
	if (*(unsigned char *)letter == 92) {
	    k = i__;
	}
	if (*(unsigned char *)letter == ']') {
	    k = i__;
	}
	if (*(unsigned char *)letter == ':') {
	    k = i__;
	}
	if (*(unsigned char *)letter == '~') {
	    k = i__;
	}
/*        if (letter .eq. '.')  k = i */
	if (k != 0) {
	    goto L10;
	}
    }
L10:
    files_1.ldir = k;

/*     read and store the keywords from the keyfile */

    getkey_();

/*     get the information level and output style */

    control_();
    return 0;
} /* basefile_ */

