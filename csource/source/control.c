/* control.f -- translated by f2c (version 20050501).
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

struct {
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

struct {
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine control  --  set information and output types  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "control" gets initial values for parameters that determine */
/*     the output style and information level provided by TINKER */


/* Subroutine */ int control_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);

    /* Local variables */
    static integer i__, next;
    static logical exist;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120], keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___7 = { 1, string, 0, 0, 120, 1 };



#define arg_ref(a_0,a_1) &argue_1.arg[(a_1)*120 + a_0 - 0]
#define keyline_ref(a_0,a_1) &keys_1.keyline[(a_1)*120 + a_0 - 120]



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  sizes.i  --  parameter values to set array dimensions  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "sizes.i" sets values for critical array dimensions used */
/*     throughout the software; these parameters will fix the size */
/*     of the largest systems that can be handled; values too large */
/*     for the computer's memory and/or swap space to accomodate */
/*     will result in poor performance or outright failure */

/*     parameter:      maximum allowed number of: */

/*     maxatm          atoms in the molecular system */
/*     maxval          atoms directly bonded to an atom */
/*     maxgrp          user-defined groups of atoms */
/*     maxref          stored reference molecular systems */
/*     maxtyp          force field atom type definitions */
/*     maxclass        force field atom class definitions */
/*     maxprm          lines in the parameter file */
/*     maxkey          lines in the keyword file */
/*     maxrot          bonds for torsional rotation */
/*     maxvar          optimization variables (vector storage) */
/*     maxopt          optimization variables (matrix storage) */
/*     maxhess         off-diagonal Hessian elements */
/*     maxlight        sites for method of lights neighbors */
/*     maxvlst         atom neighbors in van der Waals pair list */
/*     maxelst         atom neighbors in electrostatics pair list */
/*     maxfft          grid points in each FFT dimension */
/*     maxfix          geometric constraints and restraints */
/*     maxvib          vibrational frequencies */
/*     maxgeo          distance geometry points */
/*     maxcell         unit cells in replicated crystal */
/*     maxring         3-, 4-, or 5-membered rings */
/*     maxbio          biopolymer atom definitions */
/*     maxres          residues in the macromolecule */
/*     maxamino        amino acid residue types */
/*     maxnuc          nucleic acid residue types */
/*     maxbnd          covalent bonds in molecular system */
/*     maxang          bond angles in molecular system */
/*     maxtors         torsional angles in molecular system */
/*     maxbitor        bitorsions in molecular system */
/*     maxpi           atoms in conjugated pisystem */
/*     maxpib          covalent bonds involving pisystem */
/*     maxpit          torsional angles involving pisystem */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  keys.i  --  contents of current keyword parameter file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nkey      number of nonblank lines in the keyword file */
/*     keyline   contents of each individual keyword file line */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  output.i  --  control of coordinate output file format  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     archive    logical flag to save structures in an archive */
/*     noversion  logical flag governing use of filename versions */
/*     overwrite  logical flag to overwrite intermediate files inplace */
/*     cyclesave  logical flag to mark use of numbered cycle files */
/*     coordtype  selects Cartesian, internal, rigid body or none */




/*     set default values for information and output variables */

    inform_1.digits = 4;
    inform_1.abort = FALSE_;
    inform_1.verbose = FALSE_;
    inform_1.debug = FALSE_;
    inform_1.holdup = FALSE_;
    output_1.archive = FALSE_;
    output_1.noversion = FALSE_;
    output_1.overwrite = FALSE_;
    output_1.cyclesave = FALSE_;

/*     check for control parameters on the command line */

    exist = FALSE_;
    i__1 = argue_1.narg - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(string, arg_ref(0, i__), (ftnlen)120, (ftnlen)120);
	upcase_(string, (ftnlen)120);
	if (s_cmp(string, "-V", (ftnlen)2, (ftnlen)2) == 0) {
	    inform_1.verbose = TRUE_;
	} else if (s_cmp(string, "-D", (ftnlen)2, (ftnlen)2) == 0) {
	    inform_1.debug = TRUE_;
	    inform_1.verbose = TRUE_;
	}
    }

/*     search keywords for various control parameters */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "DIGITS ", (ftnlen)7, (ftnlen)7) == 0) {
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___7);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&inform_1.digits, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "VERBOSE ", (ftnlen)8, (ftnlen)8) == 0) {
	    inform_1.verbose = TRUE_;
	} else if (s_cmp(keyword, "DEBUG ", (ftnlen)6, (ftnlen)6) == 0) {
	    inform_1.debug = TRUE_;
	    inform_1.verbose = TRUE_;
	} else if (s_cmp(keyword, "EXIT-PAUSE ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    inform_1.holdup = TRUE_;
	} else if (s_cmp(keyword, "ARCHIVE ", (ftnlen)8, (ftnlen)8) == 0) {
	    output_1.archive = TRUE_;
	} else if (s_cmp(keyword, "NOVERSION ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    output_1.noversion = TRUE_;
	} else if (s_cmp(keyword, "OVERWRITE ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    output_1.overwrite = TRUE_;
	} else if (s_cmp(keyword, "SAVE-CYCLE ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    output_1.cyclesave = TRUE_;
	}
L10:
	;
    }
    return 0;
} /* control_ */

#undef keyline_ref
#undef arg_ref


