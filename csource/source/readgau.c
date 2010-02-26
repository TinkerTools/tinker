/* readgau.f -- translated by f2c (version 20050501).
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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal gx[25000], gy[25000], gz[25000], gforce[75000]	/* was [3][
	    25000] */, gh[1000000], gfreq[1000];
    integer ngatom;
} qmstuf_;

#define qmstuf_1 qmstuf_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;



/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine readgau  --  read data from G03 output file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "readgau" reads an ab initio optimized structure, forces, */
/*     Hessian and frequencies from a Gaussian 03 output file */


/* Subroutine */ int readgau_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter the Name of the Gaussian Output Fi"
	    "le :  \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(a120)";
    static char fmt_40[] = "(a120)";
    static char fmt_60[] = "(a120)";
    static char fmt_70[] = "(a120)";
    static char fmt_80[] = "(a120)";

    /* System generated locals */
    integer i__1, i__2;
    icilist ici__1;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_inqu(inlist *), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *)
	    , do_fio(integer *, char *, ftnlen), e_rsfe(void), f_open(olist *)
	    , f_rew(alist *), s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(
	    icilist *), do_lio(integer *, integer *, char *, ftnlen), e_rsli(
	    void), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int basefile_(char *, ftnlen);
    extern integer freeunit_(void);
    static char arcstart[4];
    static doublereal hessunit;
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j;
    extern /* Subroutine */ int readarcword_(integer *, char *, char *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer code, igau, itmp, jtmp, ktmp;
    static char word[120];
    static integer next, nfreq;
    static logical exist;
    static char record[120];
    static integer length, nghess;
    extern /* Subroutine */ int upcase_(char *, ftnlen), suffix_(char *, char 
	    *, ftnlen, ftnlen);
    static char string[120], gaufile[120];
    static doublereal frcunit;
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen);
    static char keyword[120];
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___8 = { 1, 0, 1, fmt_30, 0 };
    static cilist io___13 = { 1, 0, 1, fmt_30, 0 };
    static cilist io___15 = { 1, 0, 1, fmt_40, 0 };
    static icilist io___16 = { 1, record, 1, 0, 120, 1 };
    static cilist io___21 = { 1, 0, 1, fmt_60, 0 };
    static cilist io___22 = { 1, 0, 1, fmt_70, 0 };
    static cilist io___23 = { 1, 0, 1, fmt_80, 0 };
    static icilist io___24 = { 1, record, 1, 0, 120, 1 };
    static icilist io___26 = { 1, string, 1, 0, 120, 1 };



#define gforce_ref(a_1,a_2) qmstuf_1.gforce[(a_2)*3 + a_1 - 4]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2009 by Chuanjie Wu and Jay William Ponder  ## */
/*     ##                    All Rights Reserved                     ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  qmstuf.i  --  quantum data from Gaussian 03 calculation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     gx       x-coordinate of each atom in the QM data file */
/*     gy       y-coordinate of each atom in the QM data file */
/*     gz       z-coordinate of each atom in the QM data file */
/*     gforce   force components on each atom from QM data */
/*     gh       Hessian maxtrix elements from QM data */
/*     gfreq    calculated vibrational frequencies from QM data */
/*     ngatom   number of atoms in the QM data file */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  units.i  --  physical constants and unit conversions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     literature reference: */

/*     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended */
/*     Values of the Fundamental Physical Constants: 2006", Reviews of */
/*     Modern Physics, 80, 633-730 (2008) */

/*     The "2006 CODATA Recommended Values" are also available from */
/*     the NIST Reference on Constants, Units, and Uncertainty site */
/*     at http://physics.nist.gov/cuu/index.html */

/*     Most values below are derived from the 2006 CODATA reference */

/*     The conversion from calorie to Joule is the definition of the */
/*     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992) */

/*     The "coulomb" energy conversion factor is found by dimensional */
/*     analysis of Coulomb's Law, ie, by dividing the square of the */
/*     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is */
/*     the permittivity of vacuum (the "electric constant"); note that */
/*     eps0 is typically given in F/m, equivalent to C**2/(J-m) */

/*     The approximate value used for the Debye, 3.33564 x 10-30 C-m, */
/*     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997) */

/*     The value of "prescon" is based on definition of 1 atmosphere */
/*     as 101325 Pa set by the 10th Conference Generale des Poids et */
/*     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3 */

/*     avogadro    Avogadro's number (N) in particles/mole */
/*     lightspd    speed of light in vacuum (c) in cm/ps */
/*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K */
/*     gasconst    ideal gas constant (R) in kcal/mole/K */
/*     emass       mass of an electron in atomic mass units */
/*     joule       conversion from calories to joules */
/*     convert     conversion from kcal to g*Ang**2/ps**2 */
/*     bohr        conversion from Bohrs to Angstroms */
/*     hartree     conversion from Hartree to kcal/mole */
/*     evolt       conversion from Hartree to electron-volts */
/*     efreq       conversion from Hartree to cm-1 */
/*     coulomb     conversion from electron**2/Ang to kcal/mole */
/*     debye       conversion from electron-Ang to Debyes */
/*     prescon     conversion from kcal/mole/Ang**3 to Atm */




/*     initialize some values prior to opening the log file */

    exist = FALSE_;
    qmstuf_1.ngatom = 0;
    nfreq = 0;
    s_copy(arcstart, "1\\1\\", (ftnlen)4, (ftnlen)4);

/*     specify and open the Gaussian 03 output log file */

    nextarg_(gaufile, &exist, (ftnlen)120);
    if (exist) {
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = gaufile;
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
	igau = freeunit_();
	basefile_(gaufile, (ftnlen)120);
	suffix_(gaufile, "log", (ftnlen)120, (ftnlen)3);
	version_(gaufile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = gaufile;
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
	    basefile_(gaufile, (ftnlen)120);
	    suffix_(gaufile, "out", (ftnlen)120, (ftnlen)3);
	    version_(gaufile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = gaufile;
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
    while(! exist) {
	io___6.ciunit = iounit_1.iout;
	s_wsfe(&io___6);
	e_wsfe();
	io___7.ciunit = iounit_1.input;
	s_rsfe(&io___7);
	do_fio(&c__1, gaufile, (ftnlen)120);
	e_rsfe();
	igau = freeunit_();
	basefile_(gaufile, (ftnlen)120);
	suffix_(gaufile, "log", (ftnlen)120, (ftnlen)3);
	version_(gaufile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = gaufile;
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
	    basefile_(gaufile, (ftnlen)120);
	    suffix_(gaufile, "out", (ftnlen)120, (ftnlen)3);
	    version_(gaufile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = gaufile;
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
    o__1.oerr = 0;
    o__1.ounit = igau;
    o__1.ofnmlen = 120;
    o__1.ofnm = gaufile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    al__1.aerr = 0;
    al__1.aunit = igau;
    f_rew(&al__1);

/*     read each line of Gaussian output and find nonblank string */

    while(TRUE_) {
	io___8.ciunit = igau;
	i__1 = s_rsfe(&io___8);
	if (i__1 != 0) {
	    goto L120;
	}
	i__1 = do_fio(&c__1, keyword, (ftnlen)120);
	if (i__1 != 0) {
	    goto L120;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L120;
	}
	j = 0;
	for (i__ = 120; i__ >= 1; --i__) {
	    if (*(unsigned char *)&keyword[i__ - 1] != ' ') {
		j = i__;
	    }
	}
	--j;
	i__1 = 120 - j;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = i__ + j - 1;
	    s_copy(keyword + (i__ - 1), keyword + i__2, (ftnlen)1, i__ + j - 
		    i__2);
	}
	for (i__ = 121 - j; i__ <= 120; ++i__) {
	    *(unsigned char *)&keyword[i__ - 1] = ' ';
	}
	length = trimtext_(keyword, (ftnlen)120);
	upcase_(keyword, (ftnlen)120);

/*     get structure, forces and frequencies from Gaussian output */

	if (s_cmp(keyword, "STANDARD ORIENTATION", (ftnlen)20, (ftnlen)20) == 
		0) {
	    for (i__ = 1; i__ <= 4; ++i__) {
		io___13.ciunit = igau;
		i__1 = s_rsfe(&io___13);
		if (i__1 != 0) {
		    goto L120;
		}
		i__1 = do_fio(&c__1, record, (ftnlen)120);
		if (i__1 != 0) {
		    goto L120;
		}
		i__1 = e_rsfe();
		if (i__1 != 0) {
		    goto L120;
		}
	    }
	    i__ = 1;
	    while(TRUE_) {
		io___15.ciunit = igau;
		i__1 = s_rsfe(&io___15);
		if (i__1 != 0) {
		    goto L120;
		}
		i__1 = do_fio(&c__1, record, (ftnlen)120);
		if (i__1 != 0) {
		    goto L120;
		}
		i__1 = e_rsfe();
		if (i__1 != 0) {
		    goto L120;
		}
		i__1 = s_rsli(&io___16);
		if (i__1 != 0) {
		    goto L50;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&itmp, (ftnlen)sizeof(
			integer));
		if (i__1 != 0) {
		    goto L50;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&jtmp, (ftnlen)sizeof(
			integer));
		if (i__1 != 0) {
		    goto L50;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&ktmp, (ftnlen)sizeof(
			integer));
		if (i__1 != 0) {
		    goto L50;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&qmstuf_1.gx[i__ - 1], (
			ftnlen)sizeof(doublereal));
		if (i__1 != 0) {
		    goto L50;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&qmstuf_1.gy[i__ - 1], (
			ftnlen)sizeof(doublereal));
		if (i__1 != 0) {
		    goto L50;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&qmstuf_1.gz[i__ - 1], (
			ftnlen)sizeof(doublereal));
		if (i__1 != 0) {
		    goto L50;
		}
		i__1 = e_rsli();
		if (i__1 != 0) {
		    goto L50;
		}
		if (jtmp <= 0) {
		    goto L50;
		}
		++i__;
	    }
L50:
	    qmstuf_1.ngatom = i__ - 1;
	} else if (s_cmp(keyword + 36, "FORCES (HARTREES/BOHR)", (ftnlen)22, (
		ftnlen)22) == 0) {
	    frcunit = 1185.8210418245483;
	    io___21.ciunit = igau;
	    i__1 = s_rsfe(&io___21);
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = do_fio(&c__1, record, (ftnlen)120);
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L120;
	    }
	    io___22.ciunit = igau;
	    i__1 = s_rsfe(&io___22);
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = do_fio(&c__1, record, (ftnlen)120);
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = qmstuf_1.ngatom;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___23.ciunit = igau;
		i__2 = s_rsfe(&io___23);
		if (i__2 != 0) {
		    goto L120;
		}
		i__2 = do_fio(&c__1, record, (ftnlen)120);
		if (i__2 != 0) {
		    goto L120;
		}
		i__2 = e_rsfe();
		if (i__2 != 0) {
		    goto L120;
		}
		i__2 = s_rsli(&io___24);
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__3, &c__1, (char *)&itmp, (ftnlen)sizeof(
			integer));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__3, &c__1, (char *)&jtmp, (ftnlen)sizeof(
			integer));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&gforce_ref(1, i__), (
			ftnlen)sizeof(doublereal));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&gforce_ref(2, i__), (
			ftnlen)sizeof(doublereal));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&gforce_ref(3, i__), (
			ftnlen)sizeof(doublereal));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L90;
		}
		for (j = 1; j <= 3; ++j) {
		    gforce_ref(j, i__) = gforce_ref(j, i__) * frcunit;
		}
L90:
		;
	    }
	} else if (s_cmp(keyword, "FREQUENCIES --", (ftnlen)14, (ftnlen)14) ==
		 0) {
	    s_copy(string, keyword + 14, (ftnlen)120, (ftnlen)106);
	    i__1 = s_rsli(&io___26);
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&qmstuf_1.gfreq[nfreq], (
		    ftnlen)sizeof(doublereal));
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&qmstuf_1.gfreq[nfreq + 1], (
		    ftnlen)sizeof(doublereal));
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&qmstuf_1.gfreq[nfreq + 2], (
		    ftnlen)sizeof(doublereal));
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L100;
	    }
L100:
	    nfreq += 3;

/*     get the Hessian from archive section at bottom of output */

	} else if (s_cmp(keyword, arcstart, (ftnlen)4, (ftnlen)4) == 0) {
	    hessunit = 2240.8770116615283;
	    next = 1;
	    while(TRUE_) {
		readarcword_(&igau, record, word, &length, &next, (ftnlen)120,
			 (ftnlen)120);
		if (s_cmp(word, "NImag", (ftnlen)5, (ftnlen)5) == 0) {
		    for (i__ = 1; i__ <= 4; ++i__) {
			readarcword_(&igau, record, word, &length, &next, (
				ftnlen)120, (ftnlen)120);
		    }
		    nghess = qmstuf_1.ngatom * 3 * (qmstuf_1.ngatom * 3 + 1) /
			     2;
		    i__1 = nghess;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			readarcword_(&igau, record, word, &length, &next, (
				ftnlen)120, (ftnlen)120);
			ici__1.icierr = 0;
			ici__1.iciend = 0;
			ici__1.icirnum = 1;
			ici__1.icirlen = length;
			ici__1.iciunit = word;
			ici__1.icifmt = 0;
			s_rsli(&ici__1);
			do_lio(&c__5, &c__1, (char *)&qmstuf_1.gh[i__ - 1], (
				ftnlen)sizeof(doublereal));
			e_rsli();
			qmstuf_1.gh[i__ - 1] *= hessunit;
		    }
		    goto L110;
		}
		code = *(unsigned char *)word;
		if (code == 64) {
		    goto L110;
		}
	    }
L110:
	    ;
	}
    }
L120:
    cl__1.cerr = 0;
    cl__1.cunit = igau;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     zero out the frequencies if none were in Gaussian output */

    if (nfreq == 0) {
	i__1 = qmstuf_1.ngatom * 3;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    qmstuf_1.gfreq[i__ - 1] = 0.;
	}
    }
    return 0;
} /* readgau_ */

#undef gforce_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine getarcword  --  read Gaussian archive section  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "getarcword" reads data from Gaussian archive section; each */
/*     entry is terminated with a backslash symbol */

/*     igau     file unit of the Gaussian output file */
/*     word     information to be read */
/*     length   length of the word */


/* Subroutine */ int readarcword_(integer *igau, char *string, char *word, 
	integer *length, integer *next, ftnlen string_len, ftnlen word_len)
{
    /* Format strings */
    static char fmt_20[] = "(a120)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void);

    /* Local variables */
    static integer i__, code;
    static char letter[1];

    /* Fortran I/O blocks */
    static cilist io___35 = { 1, 0, 1, fmt_20, 0 };




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




/*     initialize some values prior to parsing the test string */

    *length = 1;
    *(unsigned char *)letter = ' ';
    for (i__ = 1; i__ <= 120; ++i__) {
	*(unsigned char *)&word[i__ - 1] = ' ';
    }

/*     attempt to read a text word entry from the input string */

    *(unsigned char *)letter = *(unsigned char *)&string[*next - 1];
    code = *(unsigned char *)letter;
    if (code == 92 || code == 61) {
	*(unsigned char *)word = *(unsigned char *)letter;
	++(*next);
	*length = 1;
	return 0;
    }
L10:
    for (i__ = *next; i__ <= 75; ++i__) {
	if (code == 92 || code == 61) {
	    return 0;
	}
	if (*next > 70) {
	    io___35.ciunit = *igau;
	    i__1 = s_rsfe(&io___35);
	    if (i__1 != 0) {
		goto L30;
	    }
	    i__1 = do_fio(&c__1, string, (ftnlen)120);
	    if (i__1 != 0) {
		goto L30;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L30;
	    }
	    *next = 1;
	    goto L10;
	}
	if (code == 44) {
	    ++(*next);
	    return 0;
	}
	if (code == 92 || code == 61) {
	    return 0;
	}
	*(unsigned char *)&word[*length - 1] = *(unsigned char *)letter;
	++(*next);
	*(unsigned char *)letter = *(unsigned char *)&string[*next - 1];
	code = *(unsigned char *)letter;
	++(*length);
    }
    if (code == 64) {
	*(unsigned char *)word = *(unsigned char *)letter;
	*length = 1;
    }
L30:
    return 0;
} /* readarcword_ */

