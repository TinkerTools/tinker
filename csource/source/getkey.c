/* getkey.f -- translated by f2c (version 20050501).
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
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine getkey  --  find and store contents of keyfile  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "getkey" finds a valid keyfile and stores its contents as */
/*     line images for subsequent keyword parameter searching */


/* Subroutine */ int getkey_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 GETKEY  --  Keyfile Specified\002,\002 o"
	    "n Command Line was not Found\002)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(/,\002 GETKEY  --  Keyfile Too Large;\002,\002 "
	    "Increase MAXKEY\002)";
    static char fmt_50[] = "()";
    static char fmt_60[] = "()";
    static char fmt_70[] = "(a)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2];
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), f_inqu(inlist *), s_wsfe(
	    cilist *), e_wsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_rew(alist *), s_rsfe(cilist *), do_fio(integer 
	    *, char *, ftnlen), e_rsfe(void), f_clos(cllist *);

    /* Local variables */
    extern integer freeunit_(void), trimtext_(char *, ftnlen);
    static integer i__, ikey, next;
    extern /* Subroutine */ int fatal_(void);
    static logical exist, header;
    static char record[120];
    static integer length;
    extern /* Subroutine */ int upcase_(char *, ftnlen), suffix_(char *, char 
	    *, ftnlen, ftnlen);
    static char string[120], keyfile[120], comment[120], keyword[20];
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen), 
	    gettext_(char *, char *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___7 = { 1, 0, 1, fmt_20, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_70, 0 };



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




/*     check for a keyfile specified on the command line */

    exist = FALSE_;
    i__1 = argue_1.narg - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(string, arg_ref(0, i__), (ftnlen)120, (ftnlen)120);
	upcase_(string, (ftnlen)120);
	if (s_cmp(string, "-K", (ftnlen)2, (ftnlen)2) == 0) {
	    s_copy(keyfile, arg_ref(0, i__ + 1), (ftnlen)120, (ftnlen)120);
	    suffix_(keyfile, "key", (ftnlen)120, (ftnlen)3);
	    version_(keyfile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = keyfile;
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
	    if (! exist) {
		io___5.ciunit = iounit_1.iout;
		s_wsfe(&io___5);
		e_wsfe();
		fatal_();
	    }
	}
    }

/*     try to get keyfile from base name of current system */

    if (! exist) {
/* Writing concatenation */
	i__2[0] = files_1.leng, a__1[0] = files_1.filename;
	i__2[1] = 4, a__1[1] = ".key";
	s_cat(keyfile, a__1, i__2, &c__2, (ftnlen)120);
	version_(keyfile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = keyfile;
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

/*     check for the existence of a generic keyfile */

    if (! exist) {
	if (files_1.ldir == 0) {
	    s_copy(keyfile, "tinker.key", (ftnlen)120, (ftnlen)10);
	} else {
/* Writing concatenation */
	    i__2[0] = files_1.ldir, a__1[0] = files_1.filename;
	    i__2[1] = 10, a__1[1] = "tinker.key";
	    s_cat(keyfile, a__1, i__2, &c__2, (ftnlen)120);
	}
	version_(keyfile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = keyfile;
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

/*     read the keyfile and store it for latter use */

    keys_1.nkey = 0;
    if (exist) {
	ikey = freeunit_();
	o__1.oerr = 0;
	o__1.ounit = ikey;
	o__1.ofnmlen = 120;
	o__1.ofnm = keyfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = ikey;
	f_rew(&al__1);
	while(TRUE_) {
	    io___7.ciunit = ikey;
	    i__1 = s_rsfe(&io___7);
	    if (i__1 != 0) {
		goto L40;
	    }
	    i__1 = do_fio(&c__1, record, (ftnlen)120);
	    if (i__1 != 0) {
		goto L40;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L40;
	    }
	    ++keys_1.nkey;
	    s_copy(keyline_ref(0, keys_1.nkey), record, (ftnlen)120, (ftnlen)
		    120);
	    if (keys_1.nkey >= 10000) {
		io___9.ciunit = iounit_1.iout;
		s_wsfe(&io___9);
		e_wsfe();
		fatal_();
	    }
	}
L40:
	cl__1.cerr = 0;
	cl__1.cunit = ikey;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     check for comment lines to be echoed to the output */

    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "ECHO ", (ftnlen)5, (ftnlen)5) == 0) {
	    s_copy(comment, record + (next - 1), (ftnlen)120, 120 - (next - 1)
		    );
	    length = trimtext_(comment, (ftnlen)120);
	    if (header) {
		header = FALSE_;
		io___15.ciunit = iounit_1.iout;
		s_wsfe(&io___15);
		e_wsfe();
	    }
	    if (length == 0) {
		io___16.ciunit = iounit_1.iout;
		s_wsfe(&io___16);
		e_wsfe();
	    } else {
		io___17.ciunit = iounit_1.iout;
		s_wsfe(&io___17);
		do_fio(&c__1, comment, length);
		e_wsfe();
	    }
	}
    }
    return 0;
} /* getkey_ */

#undef keyline_ref
#undef arg_ref


