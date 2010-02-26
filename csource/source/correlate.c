/* correlate.f -- translated by f2c (version 20050501).
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
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

struct {
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

struct {
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  program correlate  --  time correlation of a property  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "correlate" computes the time correlation function of some */
/*     user-supplied property from individual snapshot frames taken */
/*     from a molecular dynamics or other trajectory */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter Base Name of Coordinate Cycle File"
	    "s :  \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_40[] = "(/,\002 Numbers of First & Last File and Step"
	    "\002,\002 Increment :  \002,$)";
    static char fmt_50[] = "(a120)";
    static char fmt_80[] = "(/,\002 Maximum Frame Separation to be Used i"
	    "n\002,\002 Correlation [ALL] :  \002,$)";
    static char fmt_90[] = "(a120)";
    static char fmt_110[] = "(/,\002 Correlation Function Computed using\002"
	    ",i5,\002 Blocks of\002,i6,\002 Frames\002)";
    static char fmt_120[] = "(/,3x,\002Correlation within Frame Block :   "
	    " \002,i8)";
    static char fmt_130[] = "(3x,\002Correlation between Frame Blocks :  "
	    "\002,2i8)";
    static char fmt_140[] = "(/,3x,\002Separation\002,7x,\002Samples\002,8x"
	    ",\002Average Value\002,7x,\002Normalized\002,/)";
    static char fmt_150[] = "(i9,6x,i10,6x,2f17.6)";
    static char fmt_160[] = "(/,3x,\002Separation\002,7x,\002Samples\002,8x"
	    ",\002Average Value\002,/)";
    static char fmt_170[] = "(i9,6x,i10,6x,f17.6)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    extern doublereal property_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *)
	    ;
    static integer i__, j, k, m, n1, n2, t1[1000], t2[1000];
    static doublereal x1[1000000]	/* was [1000][1000] */, y1[1000000]	
	    /* was [1000][1000] */, z1[1000000]	/* was [1000][1000] */, x2[
	    1000000]	/* was [1000][1000] */, y2[1000000]	/* was [1000][
	    1000] */, z2[1000000]	/* was [1000][1000] */;
    static integer dt, last, step, stop;
    extern /* Subroutine */ int final_(void);
    static integer nfile;
    static doublereal value;
    static integer icorr[100001], first;
    static doublereal ucorr[100001], vcorr[100001];
    static integer start;
    static logical exist, query;
    static integer blkgap, nblock, maxgap;
    static logical normal;
    static char letter[1], string[120];
    static integer blkdiff;
    extern /* Subroutine */ int readblk_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *), 
	    initial_(void);
    static integer blksize;
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___3 = { 0, 0, 0, fmt_20, 0 };
    static icilist io___11 = { 1, string, 1, 0, 120, 1 };
    static icilist io___12 = { 1, string, 1, 0, 120, 1 };
    static icilist io___13 = { 1, string, 1, 0, 120, 1 };
    static cilist io___14 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___16 = { 1, string, 1, 0, 120, 1 };
    static icilist io___18 = { 1, string, 1, 0, 120, 1 };
    static cilist io___19 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_90, 0 };
    static icilist io___21 = { 1, string, 1, 0, 120, 1 };
    static cilist io___26 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___48 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___52 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___54 = { 0, 0, 0, fmt_170, 0 };




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
/*     ##  atoms.i  --  number, position and type of current atoms  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     x       current x-coordinate for each atom in the system */
/*     y       current y-coordinate for each atom in the system */
/*     z       current z-coordinate for each atom in the system */
/*     n       total number of atoms in the current system */
/*     type    atom type number for each atom in the system */




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




/*     get the base name of user specified input structures */

    initial_();
    nextarg_(files_1.filename, &exist, (ftnlen)120);
    if (! exist) {
	io___2.ciunit = iounit_1.iout;
	s_wsfe(&io___2);
	e_wsfe();
	io___3.ciunit = iounit_1.input;
	s_rsfe(&io___3);
	do_fio(&c__1, files_1.filename, (ftnlen)120);
	e_rsfe();
    }

/*     remove any extension from the filename */

    files_1.leng = trimtext_(files_1.filename, (ftnlen)120);
    last = files_1.leng;
    i__1 = files_1.leng;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)letter = *(unsigned char *)&files_1.filename[i__ - 
		1];
	if (*(unsigned char *)letter == '/') {
	    last = files_1.leng;
	}
/*        if (letter .eq. '\')  last = leng */
	if (*(unsigned char *)letter == 92) {
	    last = files_1.leng;
	}
	if (*(unsigned char *)letter == ']') {
	    last = files_1.leng;
	}
	if (*(unsigned char *)letter == ':') {
	    last = files_1.leng;
	}
	if (*(unsigned char *)letter == '~') {
	    last = files_1.leng;
	}
	if (*(unsigned char *)letter == '.') {
	    last = i__ - 1;
	}
    }
    files_1.leng = min(files_1.leng,last);

/*     set first and last snapshot frames and step increment */

    first = 0;
    last = 0;
    step = 0;
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___11);
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&first, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L30;
	}
	query = FALSE_;
    }
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___12);
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&last, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L30;
	}
    }
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___13);
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&step, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L30;
	}
    }
L30:
    if (query) {
	io___14.ciunit = iounit_1.iout;
	s_wsfe(&io___14);
	e_wsfe();
	io___15.ciunit = iounit_1.input;
	s_rsfe(&io___15);
	do_fio(&c__1, string, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___16);
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&first, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&last, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&step, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L60;
	}
L60:
	;
    }
    if (last == 0) {
	last = first;
    }
    if (step == 0) {
	step = 1;
    }

/*     set the maximum frame separation to be used for correlation */

    maxgap = last - first;
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___18);
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&maxgap, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L70;
	}
	query = FALSE_;
    }
L70:
    if (query) {
	io___19.ciunit = iounit_1.iout;
	s_wsfe(&io___19);
	e_wsfe();
	io___20.ciunit = iounit_1.input;
	s_rsfe(&io___20);
	do_fio(&c__1, string, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___21);
	if (i__1 != 0) {
	    goto L100;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&maxgap, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L100;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L100;
	}
L100:
	;
    }
    if (maxgap == 0) {
	maxgap = last - first;
    }

/*     get the number of file blocks from the total files */

    nfile = (last - first) / step + 1;
    nblock = (nfile - 1) / 1000 + 1;
    blksize = step * 1000;
    blkgap = (maxgap - 1) / blksize + 1;
    io___26.ciunit = iounit_1.iout;
    s_wsfe(&io___26);
    do_fio(&c__1, (char *)&nblock, (ftnlen)sizeof(integer));
    i__1 = min(nfile,1000);
    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
    e_wsfe();

/*     zero out the time correlation function cumulative values */

    i__1 = maxgap;
    for (i__ = 0; i__ <= i__1; ++i__) {
	icorr[i__] = 0;
	vcorr[i__] = 0.;
    }

/*     cycle over all pairs of snapshot frame blocks */

    i__1 = nblock;
    for (i__ = 1; i__ <= i__1; ++i__) {
	start = first + (i__ - 1) * blksize;
	stop = start + blksize - step;
	stop = min(last,stop);
	readblk_(&start, &stop, &step, &n1, t1, x1, y1, z1);
	io___36.ciunit = iounit_1.iout;
	s_wsfe(&io___36);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();

/*     compute time correlation for frames within single block */

	i__2 = n1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = n1;
	    for (m = k; m <= i__3; ++m) {
		dt = t1[m - 1] - t1[k - 1];
		if (dt <= maxgap) {
		    value = property_(&k, x1, y1, z1, &m, x1, y1, z1);
		    ++icorr[dt];
		    vcorr[dt] += value;
		}
	    }
	}

/*     compute time correlation for frames between two blocks */

/* Computing MIN */
	i__3 = i__ + blkgap;
	i__2 = min(i__3,nblock);
	for (j = i__ + 1; j <= i__2; ++j) {
	    start = first + (j - 1) * blksize;
	    stop = start + blksize - step;
	    stop = min(last,stop);
	    blkdiff = (j - i__) * 1000;
	    readblk_(&start, &stop, &step, &n2, t2, x2, y2, z2);
	    io___48.ciunit = iounit_1.iout;
	    s_wsfe(&io___48);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	    e_wsfe();
	    i__3 = n1;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = n2;
		for (m = 1; m <= i__4; ++m) {
		    dt = t2[m - 1] - t1[k - 1] + blkdiff;
		    if (dt <= maxgap) {
			value = property_(&k, x1, y1, z1, &m, x2, y2, z2);
			++icorr[dt];
			vcorr[dt] += value;
		    }
		}
	    }
	}
    }

/*     compute the average correlation function values */

    i__1 = maxgap;
    for (i__ = 0; i__ <= i__1; ++i__) {
	if (icorr[i__] != 0) {
	    vcorr[i__] /= (doublereal) icorr[i__];
	}
    }

/*     get the normalized correlation function if applicable */

    normal = FALSE_;
    if (vcorr[0] != 0.) {
	normal = TRUE_;
    }
    if (normal) {
	i__1 = maxgap;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    ucorr[i__] = vcorr[i__] / vcorr[0];
	}
    }

/*     print the final values of the correlation function */

    if (normal) {
	io___51.ciunit = iounit_1.iout;
	s_wsfe(&io___51);
	e_wsfe();
	i__1 = maxgap;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    if (icorr[i__] != 0) {
		io___52.ciunit = iounit_1.iout;
		s_wsfe(&io___52);
		i__2 = i__ * step;
		do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&icorr[i__], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&vcorr[i__], (ftnlen)sizeof(doublereal))
			;
		do_fio(&c__1, (char *)&ucorr[i__], (ftnlen)sizeof(doublereal))
			;
		e_wsfe();
	    }
	}
    } else {
	io___53.ciunit = iounit_1.iout;
	s_wsfe(&io___53);
	e_wsfe();
	i__1 = maxgap;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    if (icorr[i__] != 0) {
		io___54.ciunit = iounit_1.iout;
		s_wsfe(&io___54);
		i__2 = i__ * step;
		do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&icorr[i__], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&vcorr[i__], (ftnlen)sizeof(doublereal))
			;
		e_wsfe();
	    }
	}
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine readblk  --  read a block of snapshot frames  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "readblk" reads in a set of snapshot frames and transfers */
/*     the values to internal arrays for use in the computation */
/*     of time correlation functions */


/* Subroutine */ int readblk_(integer *start, integer *stop, integer *step, 
	integer *nb, integer *tb, doublereal *xb, doublereal *yb, doublereal *
	zb)
{
    /* Format strings */
    static char fmt_10[] = "(a120)";
    static char fmt_20[] = "(/,\002 READBLK  --  Too many Correlation Sites"
	    ";\002,\002 Increase MAXSITE\002)";
    static char fmt_40[] = "(a120)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2, i__3[3], i__4;
    olist o__1;
    cllist cl__1;
    inlist ioin__1;

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_inqu(inlist *), f_open(olist *), s_rsfe(cilist *), do_fio(
	    integer *, char *, ftnlen), e_rsfe(void), s_rsli(icilist *), 
	    do_lio(integer *, integer *, char *, ftnlen), e_rsli(void), 
	    s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_clos(cllist *);

    /* Local variables */
    extern integer freeunit_(void);
    static integer i__, k, nt;
    static char ext[7];
    static integer lext, next, ixyz, label;
    extern /* Subroutine */ int fatal_(void);
    static logical exist;
    static char record[120], string[120];
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    , getword_(char *, char *, integer *, ftnlen, ftnlen);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static cilist io___62 = { 0, 0, 0, fmt_10, 0 };
    static icilist io___64 = { 0, record, 0, 0, 120, 1 };
    static cilist io___65 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___68 = { 0, 0, 0, fmt_40, 0 };
    static icilist io___69 = { 0, record, 0, 0, 120, 1 };
    static icilist io___72 = { 0, string, 0, 0, 120, 1 };



#define xb_ref(a_1,a_2) xb[(a_2)*1000 + a_1]
#define yb_ref(a_1,a_2) yb[(a_2)*1000 + a_1]
#define zb_ref(a_1,a_2) zb[(a_2)*1000 + a_1]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]



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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  atmtyp.i  --  atomic properties for each current atom  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     mass      atomic weight for each atom in the system */
/*     tag       integer atom labels from input coordinates file */
/*     class     atom class number for each atom in the system */
/*     atomic    atomic number for each atom in the system */
/*     valence   valence number for each atom in the system */
/*     name      atom name for each atom in the system */
/*     story     descriptive type for each atom in system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  atoms.i  --  number, position and type of current atoms  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     x       current x-coordinate for each atom in the system */
/*     y       current y-coordinate for each atom in the system */
/*     z       current z-coordinate for each atom in the system */
/*     n       total number of atoms in the current system */
/*     type    atom type number for each atom in the system */




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




/*     initialize the number of files and the numeral size */

    /* Parameter adjustments */
    zb -= 1001;
    yb -= 1001;
    xb -= 1001;
    --tb;

    /* Function Body */
    nt = 0;
    *nb = 0;
    lext = 3;

/*     cycle over all snapshot frames in the block of files */

    i__1 = *stop;
    i__2 = *step;
    for (i__ = *start; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	++nt;
	numeral_(&i__, ext, &lext, (ftnlen)7);
/* Writing concatenation */
	i__3[0] = files_1.leng, a__1[0] = files_1.filename;
	i__3[1] = 1, a__1[1] = ".";
	i__3[2] = lext, a__1[2] = ext;
	s_cat(xyzfile, a__1, i__3, &c__3, (ftnlen)120);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = xyzfile;
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

/*     add file to the current block and get number of atoms */

	if (exist) {
	    ++(*nb);
	    tb[*nb] = nt;
	    ixyz = freeunit_();
	    o__1.oerr = 0;
	    o__1.ounit = ixyz;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = xyzfile;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	    io___62.ciunit = ixyz;
	    s_rsfe(&io___62);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    s_rsli(&io___64);
	    do_lio(&c__3, &c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
	    e_rsli();

/*     check for too many correlation sites in the frame */

	    if (atoms_1.n > 1000) {
		io___65.ciunit = iounit_1.iout;
		s_wsfe(&io___65);
		e_wsfe();
		fatal_();
	    }

/*     read the frame in the TINKER-generated coordinate format; */
/*     this is fast, but assumes the fixed format shown below */

/*           do k = 1, n */
/*              read (ixyz,30)  name(k),xb(k,nb),yb(k,nb),zb(k,nb) */
/*  30          format (8x,a3,3f12.6) */
/*           end do */

/*     alternatively, get each frame from a free formated file; */
/*     this is slow, but correctly handles any valid TINKER file */

	    i__4 = atoms_1.n;
	    for (k = 1; k <= i__4; ++k) {
		next = 1;
		io___68.ciunit = ixyz;
		s_rsfe(&io___68);
		do_fio(&c__1, record, (ftnlen)120);
		e_rsfe();
		s_rsli(&io___69);
		do_lio(&c__3, &c__1, (char *)&label, (ftnlen)sizeof(integer));
		e_rsli();
		getword_(record, name___ref(0, k), &next, (ftnlen)120, (
			ftnlen)3);
		s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next 
			- 1));
		s_rsli(&io___72);
		do_lio(&c__5, &c__1, (char *)&xb_ref(k, *nb), (ftnlen)sizeof(
			doublereal));
		do_lio(&c__5, &c__1, (char *)&yb_ref(k, *nb), (ftnlen)sizeof(
			doublereal));
		do_lio(&c__5, &c__1, (char *)&zb_ref(k, *nb), (ftnlen)sizeof(
			doublereal));
		e_rsli();
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = ixyz;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	}
    }
    return 0;
} /* readblk_ */

#undef name___ref
#undef zb_ref
#undef yb_ref
#undef xb_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  function property  --  compute correlation property value  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "property" takes two input snapshot frames and computes the */
/*     value of the property for which the correlation function is */
/*     being accumulated */

/*     this version of "property" finds the velocity autocorrelation */
/*     or the rms fit as a function of time, and is merely provided */
/*     as an example; the user will need to write a similar custom */
/*     function to compute other properties to be correlated */


doublereal property_(integer *i__, doublereal *xi, doublereal *yi, doublereal 
	*zi, integer *k, doublereal *xk, doublereal *yk, doublereal *zk)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer j;
    static doublereal x1[1000], y1[1000], z1[1000], x2[1000], y2[1000], z2[
	    1000], value;
    extern /* Subroutine */ int impose_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *);


#define xi_ref(a_1,a_2) xi[(a_2)*1000 + a_1]
#define yi_ref(a_1,a_2) yi[(a_2)*1000 + a_1]
#define zi_ref(a_1,a_2) zi[(a_2)*1000 + a_1]
#define xk_ref(a_1,a_2) xk[(a_2)*1000 + a_1]
#define yk_ref(a_1,a_2) yk[(a_2)*1000 + a_1]
#define zk_ref(a_1,a_2) zk[(a_2)*1000 + a_1]



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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  atoms.i  --  number, position and type of current atoms  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     x       current x-coordinate for each atom in the system */
/*     y       current y-coordinate for each atom in the system */
/*     z       current z-coordinate for each atom in the system */
/*     n       total number of atoms in the current system */
/*     type    atom type number for each atom in the system */




/*     transfer the input trajectory frames to local vectors */

    /* Parameter adjustments */
    zk -= 1001;
    yk -= 1001;
    xk -= 1001;
    zi -= 1001;
    yi -= 1001;
    xi -= 1001;

    /* Function Body */
    value = 0.;
    i__1 = atoms_1.n;
    for (j = 1; j <= i__1; ++j) {
	x1[j - 1] = xi_ref(j, *i__);
	y1[j - 1] = yi_ref(j, *i__);
	z1[j - 1] = zi_ref(j, *i__);
	x2[j - 1] = xk_ref(j, *k);
	y2[j - 1] = yk_ref(j, *k);
	z2[j - 1] = zk_ref(j, *k);
    }

/*     sample code to find the velocity autocorrelation function */

/*     do j = 1, n */
/*        value = value + x1(j)*x2(j) + y1(j)*y2(j) + z1(j)*z2(j) */
/*     end do */

/*     sample code to find the rms deviation upon superposition */

    impose_(&atoms_1.n, x1, y1, z1, &atoms_1.n, x2, y2, z2, &value);

/*     set property value to be returned for this frame pair */

    ret_val = value;
    return ret_val;
} /* property_ */

#undef zk_ref
#undef yk_ref
#undef xk_ref
#undef zi_ref
#undef yi_ref
#undef xi_ref


/* Main program alias */ int correlate_ () { MAIN__ (); return 0; }
