/* getprm.f -- translated by f2c (version 20050501).
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

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

struct {
    integer nprm;
    char prmline[3000000];
} params_;

#define params_1 params_

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine getprm  --  get force field parameter file  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "getprm" finds the potential energy parameter file */
/*     and then opens and reads the parameters */


/* Subroutine */ int getprm_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter Potential Parameter File Name : "
	    " \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(a120)";
    static char fmt_40[] = "(/,\002 GETPRM  --  Parameter File Too Large;"
	    "\002,\002 Increase MAXPRM\002)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen),
	     s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), f_inqu(inlist *), s_wsfe(
	    cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *, char 
	    *, ftnlen), e_rsfe(void), f_open(olist *), f_rew(alist *), f_clos(
	    cllist *);

    /* Local variables */
    extern integer freeunit_(void);
    static integer i__;
    extern /* Subroutine */ int getstring_(char *, char *, integer *, ftnlen, 
	    ftnlen);
    static char none[4];
    static integer iprm, next;
    extern /* Subroutine */ int fatal_(void);
    static logical exist;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), suffix_(char *, char 
	    *, ftnlen, ftnlen);
    static char string[120];
    static logical useprm;
    extern /* Subroutine */ int readprm_(void);
    static char prmfile[120];
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), initprm_(
	    void);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___13 = { 1, 0, 1, fmt_30, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_40, 0 };



#define keyline_ref(a_0,a_1) &keys_1.keyline[(a_1)*120 + a_0 - 120]
#define prmline_ref(a_0,a_1) &params_1.prmline[(a_1)*120 + a_0 - 120]



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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  params.i  --  contents of force field parameter file  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     nprm      number of nonblank lines in the parameter file */
/*     prmline   contents of each individual parameter file line */




/*     set the default name for the parameter file */

    useprm = TRUE_;
/* Writing concatenation */
    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
    i__1[1] = 4, a__1[1] = ".prm";
    s_cat(prmfile, a__1, i__1, &c__2, (ftnlen)120);

/*     search the keyword list for the parameter filename */

    i__2 = keys_1.nkey;
    for (i__ = 1; i__ <= i__2; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "PARAMETERS ", (ftnlen)11, (ftnlen)11) == 0) {
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    next = 1;
	    getstring_(string, prmfile, &next, (ftnlen)120, (ftnlen)120);
	    if (next == 1) {
		gettext_(string, prmfile, &next, (ftnlen)120, (ftnlen)120);
	    }
	}
    }

/*     check existence of default or specified parameter file */

    suffix_(prmfile, "prm", (ftnlen)120, (ftnlen)3);
    version_(prmfile, "old", (ftnlen)120, (ftnlen)3);
    ioin__1.inerr = 0;
    ioin__1.infilen = 120;
    ioin__1.infile = prmfile;
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

/*     test for user specified absence of a parameter file */

    if (! exist) {
	s_copy(none, prmfile, (ftnlen)4, (ftnlen)4);
	upcase_(none, (ftnlen)4);
	if (s_cmp(none, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    exist = TRUE_;
	    useprm = FALSE_;
	}
    }

/*     try to get a parameter filename from the command line */

    if (! exist) {
	nextarg_(prmfile, &exist, (ftnlen)120);
	if (exist) {
	    suffix_(prmfile, "prm", (ftnlen)120, (ftnlen)3);
	    version_(prmfile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = prmfile;
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
    }

/*     if necessary, ask for the parameter filename */

    while(! exist) {
	io___10.ciunit = iounit_1.iout;
	s_wsfe(&io___10);
	e_wsfe();
	io___11.ciunit = iounit_1.input;
	s_rsfe(&io___11);
	do_fio(&c__1, prmfile, (ftnlen)120);
	e_rsfe();
	suffix_(prmfile, "prm", (ftnlen)120, (ftnlen)3);
	version_(prmfile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = prmfile;
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

/*     initialize force field control and parameter values */

    initprm_();

/*     read the parameter file and store it for latter use */

    params_1.nprm = 0;
    if (useprm) {
	iprm = freeunit_();
	o__1.oerr = 0;
	o__1.ounit = iprm;
	o__1.ofnmlen = 120;
	o__1.ofnm = prmfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = iprm;
	f_rew(&al__1);
	while(TRUE_) {
	    io___13.ciunit = iprm;
	    i__2 = s_rsfe(&io___13);
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_fio(&c__1, record, (ftnlen)120);
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = e_rsfe();
	    if (i__2 != 0) {
		goto L50;
	    }
	    ++params_1.nprm;
	    s_copy(prmline_ref(0, params_1.nprm), record, (ftnlen)120, (
		    ftnlen)120);
	    if (params_1.nprm >= 25000) {
		io___14.ciunit = iounit_1.iout;
		s_wsfe(&io___14);
		e_wsfe();
		fatal_();
	    }
	}
L50:
	cl__1.cerr = 0;
	cl__1.cunit = iprm;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     get control and parameter values from the parameter file */

    if (useprm) {
	readprm_();
    }
    return 0;
} /* getprm_ */

#undef prmline_ref
#undef keyline_ref


