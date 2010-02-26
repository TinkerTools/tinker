/* pss.f -- translated by f2c (version 20050501).
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
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal dihed[1000];
    integer nomega, iomega[2000]	/* was [2][1000] */, zline[1000];
} omega_;

#define omega_1 omega_

struct {
    doublereal xref[250000]	/* was [25000][10] */, yref[250000]	/* 
	    was [25000][10] */, zref[250000]	/* was [25000][10] */;
    integer nref[10], reftyp[250000]	/* was [25000][10] */, n12ref[250000]	
	    /* was [25000][10] */, i12ref[2000000]	/* was [8][25000][10] 
	    */, refleng[10], refltitle[10];
    char refnam[750000]	/* was [25000][10] */, reffile[1200], reftitle[1200];
} refer_;

#define refer_1 refer_

struct {
    doublereal etree, ilevel[501];
    integer nlevel;
} tree_;

#define tree_1 tree_

struct {
    doublereal m2[25000], deform, difft, diffv, diffc;
    logical use_smooth__, use_dem__, use_gda__, use_tophat__, use_stophat__;
} warp_;

#define warp_1 warp_

struct {
    doublereal hesscut;
} hescut_;

#define hescut_1 hescut_

struct {
    doublereal zbond[25000], zang[25000], ztors[25000];
    integer iz[100000]	/* was [4][25000] */;
} zcoord_;

#define zcoord_1 zcoord_

struct {
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__5 = 5;
static doublereal c_b70 = 12.;
static integer c__1000 = 1000;



/*     ############################################################## */
/*     ##  COPYRIGHT (C) 1997 by Rohit Pappu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  program pss  --  Cartesian potential smoothing & search  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "pss" implements the potential smoothing plus search method */
/*     for global optimization in Cartesian coordinate space with */
/*     local searches performed in Cartesian or torsional space */

/*     literature reference: */

/*     J. Kostrowicki and H. A. Scheraga, "Application of the Diffusion */
/*     Equation Method for Global Optimization to Oligopeptides", Journal */
/*     of Physical Chemistry, 96, 7442-7449 (1992) */

/*     S. Nakamura, H. Hirose, M. Ikeguchi and J. Doi, "Conformational */
/*     Energy Minimization Using a Two-Stage Method", Journal of Physical */
/*     Chemistry, 99, 8374-8378 (1995) */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Enter the Number of Steps for Smoothing "
	    "Schedule\002,\002 [100] :  \002,$)";
    static char fmt_30[] = "(i10)";
    static char fmt_40[] = "(/,\002 Perform Forward Smoothing from Input Str"
	    "ucture\002,\002 [Y] :  \002,$)";
    static char fmt_50[] = "(a120)";
    static char fmt_60[] = "(/,\002 Use Quadratic, Cubic or Sigmoidal Schedu"
	    "le\002,\002 (Q [C] or S) :  \002,$)";
    static char fmt_70[] = "(a120)";
    static char fmt_80[] = "(/,\002 Local Search Type - Cartesian, Torsional"
	    " or None\002,\002 (C T or [N]) :  \002,$)";
    static char fmt_90[] = "(a120)";
    static char fmt_110[] = "(/,\002 Enter the Range of Local Search Directi"
	    "ons\002,\002 (1=Highest Freq) :  \002,$)";
    static char fmt_120[] = "(a120)";
    static char fmt_140[] = "(/,\002 Enter the Largest Smoothing Level fo"
	    "r\002,\002 Local Search [5.0] :  \002,$)";
    static char fmt_150[] = "(f20.0)";
    static char fmt_160[] = "(/,\002 Restrict Local Search to Children of In"
	    "put\002,\002 Structure [N] :  \002,$)";
    static char fmt_170[] = "(a120)";
    static char fmt_190[] = "(/,\002 Enter RMS Gradient per Atom Criterio"
	    "n\002,\002 [0.0001] :  \002,$)";
    static char fmt_200[] = "(f20.0)";
    static char fmt_210[] = "(/,\002 Final Function Value and Deformation "
	    ":\002,2f15.4)";
    static char fmt_220[] = "(/,\002 Final Function Value and Deformation "
	    ":\002,2f15.4)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void), modecart_(integer *, integer 
	    *, doublereal *, doublereal *, logical *);
    static logical use_cart__;
    extern /* Subroutine */ int modetors_(integer *, integer *, doublereal *, 
	    doublereal *, logical *);
    static logical use_tors__;
    extern /* Subroutine */ int localxyz_(doublereal *, doublereal *), 
	    psswrite_(integer *);
    static integer i__;
    static doublereal rms;
    static logical use_forward__;
    static integer next, stop;
    static logical check;
    extern /* Subroutine */ int final_(void);
    static integer range;
    static doublereal ratio;
    static logical exist;
    static integer start;
    extern /* Subroutine */ int active_(void);
    static char record[120];
    static doublereal grdmin;
    static char answer[1], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), impose_(integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *), getxyz_(void), 
	    makeref_(integer *), makeint_(integer *), initial_(void);
    extern doublereal sigmoid_(doublereal *, doublereal *);
    static doublereal srchmax;
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen);
    static doublereal minimum;
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), initrot_(void);
    static char formtyp[1];

    /* Fortran I/O blocks */
    static icilist io___3 = { 1, string, 1, 0, 120, 1 };
    static cilist io___4 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_90, 0 };
    static icilist io___21 = { 1, string, 1, 0, 120, 1 };
    static icilist io___22 = { 1, string, 1, 0, 120, 1 };
    static cilist io___23 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_120, 0 };
    static icilist io___25 = { 0, record, 0, 0, 120, 1 };
    static icilist io___28 = { 1, string, 1, 0, 120, 1 };
    static cilist io___29 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_170, 0 };
    static icilist io___35 = { 1, string, 1, 0, 120, 1 };
    static cilist io___36 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_220, 0 };




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  refer.i  --  storage of reference atomic coordinate set  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xref        reference x-coordinates for atoms in each system */
/*     yref        reference y-coordinates for atoms in each system */
/*     zref        reference z-coordinates for atoms in each system */
/*     nref        total number of atoms in each reference system */
/*     reftyp      atom types of the atoms in each reference system */
/*     n12ref      number of atoms bonded to each reference atom */
/*     i12ref      atom numbers of atoms 1-2 connected to each atom */
/*     refleng     length in characters of each reference filename */
/*     refltitle   length in characters of each reference title line */
/*     refnam      atom names of the atoms in each reference system */
/*     reffile     base filename for each reference system */
/*     reftitle    title used to describe each reference system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  tree.i  --  potential smoothing and search tree levels  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     maxpss   maximum number of potential smoothing levels */

/*     etree    energy reference value at the top of the tree */
/*     ilevel   smoothing deformation value at each tree level */
/*     nlevel   number of levels of potential smoothing used */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  warp.i  --  parameters for potential surface smoothing  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     m2           second moment of the GDA gaussian for each atom */
/*     deform       value of the smoothing deformation parameter */
/*     difft        diffusion coefficient for torsional potential */
/*     diffv        diffusion coefficient for van der Waals potential */
/*     diffc        diffusion coefficient for charge-charge potential */
/*     use_smooth   flag to use a potential energy smoothing method */
/*     use_dem      flag to use diffusion equation method potential */
/*     use_gda      flag to use gaussian density annealing potential */
/*     use_tophat   flag to use analytical tophat smoothed potential */
/*     use_stophat  flag to use shifted tophat smoothed potential */




/*     set up the structure, mechanics calculation and smoothing */

    initial_();
    getxyz_();
    warp_1.use_smooth__ = TRUE_;
    warp_1.use_dem__ = TRUE_;
    mechanic_();
    inform_1.iwrite = 0;

/*     get the number of points along the deformation schedule */

    tree_1.nlevel = -1;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___3);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&tree_1.nlevel, (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    }
L10:
    if (tree_1.nlevel < 0) {
	io___4.ciunit = iounit_1.iout;
	s_wsfe(&io___4);
	e_wsfe();
	io___5.ciunit = iounit_1.input;
	s_rsfe(&io___5);
	do_fio(&c__1, (char *)&tree_1.nlevel, (ftnlen)sizeof(integer));
	e_rsfe();
	if (tree_1.nlevel <= 0) {
	    tree_1.nlevel = 100;
	}
    }

/*     decide whether to use forward smoothing of initial structure */

    use_forward__ = TRUE_;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___8.ciunit = iounit_1.iout;
	s_wsfe(&io___8);
	e_wsfe();
	io___9.ciunit = iounit_1.input;
	s_rsfe(&io___9);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'N') {
	use_forward__ = FALSE_;
    }

/*     get the functional form for the deformation schedule */

    *(unsigned char *)formtyp = 'C';
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___13.ciunit = iounit_1.iout;
	s_wsfe(&io___13);
	e_wsfe();
	io___14.ciunit = iounit_1.input;
	s_rsfe(&io___14);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'Q') {
	*(unsigned char *)formtyp = *(unsigned char *)answer;
    }
    if (*(unsigned char *)answer == 'S') {
	*(unsigned char *)formtyp = *(unsigned char *)answer;
    }

/*     decide which type of local search procedure to use */

    use_cart__ = FALSE_;
    use_tors__ = FALSE_;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___17.ciunit = iounit_1.iout;
	s_wsfe(&io___17);
	e_wsfe();
	io___18.ciunit = iounit_1.input;
	s_rsfe(&io___18);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'C') {
	use_cart__ = TRUE_;
    }
    if (*(unsigned char *)answer == 'T') {
	use_tors__ = TRUE_;
    }

/*     get the rotatable bonds for torsional local search */

    if (use_tors__) {
	makeint_(&c__0);
	initrot_();
	active_();
    }

/*     get the number of eigenvectors to use for local search */

    if (use_cart__ || use_tors__) {
	start = -1;
	stop = -1;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___21);
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&start, (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L100;
	    }
	}
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___22);
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&stop, (ftnlen)sizeof(integer)
		    );
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L100;
	    }
	}
L100:
	if (stop <= 0) {
	    io___23.ciunit = iounit_1.iout;
	    s_wsfe(&io___23);
	    e_wsfe();
	    io___24.ciunit = iounit_1.input;
	    s_rsfe(&io___24);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    s_rsli(&io___25);
	    do_lio(&c__3, &c__1, (char *)&start, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&stop, (ftnlen)sizeof(integer));
	    e_rsli();
	    range = (i__1 = stop - start, abs(i__1));
	    start = min(start,stop);
	    stop = start + range;
	}
	if (use_cart__) {
/* Computing MIN */
	    i__1 = stop, i__2 = atoms_1.n * 3 - 6;
	    stop = min(i__1,i__2);
	}
	if (use_tors__) {
	    stop = min(stop,omega_1.nomega);
	}
    }

/*     get the maximal smoothing level for use of local search */

    if (use_cart__ || use_tors__) {
	srchmax = -1.;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___28);
	    if (i__1 != 0) {
		goto L130;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&srchmax, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L130;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L130;
	    }
	}
L130:
	if (srchmax < 0.) {
	    io___29.ciunit = iounit_1.iout;
	    s_wsfe(&io___29);
	    e_wsfe();
	    io___30.ciunit = iounit_1.input;
	    s_rsfe(&io___30);
	    do_fio(&c__1, (char *)&srchmax, (ftnlen)sizeof(doublereal));
	    e_rsfe();
	    if (srchmax < 0.) {
		srchmax = 5.;
	    }
	}
    }

/*     decide whether to use forward smoothing of initial structure */

    check = FALSE_;
    if ((use_cart__ || use_tors__) && ! use_forward__) {
	nextarg_(answer, &exist, (ftnlen)1);
	if (! exist) {
	    io___32.ciunit = iounit_1.iout;
	    s_wsfe(&io___32);
	    e_wsfe();
	    io___33.ciunit = iounit_1.input;
	    s_rsfe(&io___33);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    next = 1;
	    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
	}
	upcase_(answer, (ftnlen)1);
	if (*(unsigned char *)answer == 'Y') {
	    check = TRUE_;
	}
    }

/*     get the termination criterion as RMS gradient per atom */

    grdmin = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___35);
	if (i__1 != 0) {
	    goto L180;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&grdmin, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L180;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L180;
	}
    }
L180:
    if (grdmin <= 0.) {
	io___36.ciunit = iounit_1.iout;
	s_wsfe(&io___36);
	e_wsfe();
	io___37.ciunit = iounit_1.input;
	s_rsfe(&io___37);
	do_fio(&c__1, (char *)&grdmin, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (grdmin <= 0.) {
	grdmin = 1e-4;
    }

/*     compute the smoothing levels for the desired protocol */

    i__1 = tree_1.nlevel;
    for (i__ = 0; i__ <= i__1; ++i__) {
	ratio = 1. - (doublereal) (tree_1.nlevel - i__) / (doublereal) 
		tree_1.nlevel;
	if (*(unsigned char *)formtyp == 'Q') {
/* Computing 2nd power */
	    d__1 = ratio;
	    tree_1.ilevel[i__] = warp_1.deform * (d__1 * d__1);
	} else if (*(unsigned char *)formtyp == 'C') {
/* Computing 3rd power */
	    d__1 = ratio;
	    tree_1.ilevel[i__] = warp_1.deform * (d__1 * (d__1 * d__1));
	} else if (*(unsigned char *)formtyp == 'S') {
	    tree_1.ilevel[i__] = warp_1.deform * sigmoid_(&c_b70, &ratio);
	}
    }

/*     perform forward PSS by looping over smoothed surfaces */

    if (use_forward__) {
	i__1 = tree_1.nlevel - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    warp_1.deform = tree_1.ilevel[i__];
	    makeref_(&c__1);
	    inform_1.iprint = 1;
	    localxyz_(&minimum, &grdmin);
	    impose_(&atoms_1.n, refer_1.xref, refer_1.yref, refer_1.zref, &
		    atoms_1.n, atoms_1.x, atoms_1.y, atoms_1.z__, &rms);
	    psswrite_(&i__);
	    io___42.ciunit = iounit_1.iout;
	    s_wsfe(&io___42);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&warp_1.deform, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     perform PSS reversal by looping over smoothed surfaces */

    for (i__ = tree_1.nlevel; i__ >= 0; --i__) {
	warp_1.deform = tree_1.ilevel[i__];
	makeref_(&c__1);
	inform_1.iprint = 1;
	localxyz_(&minimum, &grdmin);
	impose_(&atoms_1.n, refer_1.xref, refer_1.yref, refer_1.zref, &
		atoms_1.n, atoms_1.x, atoms_1.y, atoms_1.z__, &rms);
	if (i__ == tree_1.nlevel) {
	    tree_1.etree = minimum;
	}
	if (warp_1.deform <= srchmax) {
	    if (use_cart__) {
		modecart_(&start, &stop, &minimum, &grdmin, &check);
	    } else if (use_tors__) {
		modetors_(&start, &stop, &minimum, &grdmin, &check);
	    }
	}
	if (use_forward__) {
	    i__1 = (tree_1.nlevel << 1) - i__;
	    psswrite_(&i__1);
	} else {
	    i__1 = tree_1.nlevel - i__;
	    psswrite_(&i__1);
	}
	io___43.ciunit = iounit_1.iout;
	s_wsfe(&io___43);
	do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&warp_1.deform, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  function pss1  --  energy and gradient values for PSS  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "pss1" is a service routine that computes the energy */
/*     and gradient during PSS global optimization in Cartesian */
/*     coordinate space */


doublereal pss1_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static doublereal e;
    static integer i__, nvar;
    static doublereal derivs[75000]	/* was [3][25000] */;


#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1 - 4]



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




/*     translate optimization parameters to atomic coordinates */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	atoms_1.x[i__ - 1] = xx[nvar];
	++nvar;
	atoms_1.y[i__ - 1] = xx[nvar];
	++nvar;
	atoms_1.z__[i__ - 1] = xx[nvar];
    }

/*     compute and store the energy and gradient */

    gradient_(&e, derivs);
    ret_val = e;

/*     store Cartesian gradient as optimization gradient */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	g[nvar] = derivs_ref(1, i__);
	++nvar;
	g[nvar] = derivs_ref(2, i__);
	++nvar;
	g[nvar] = derivs_ref(3, i__);
    }
    return ret_val;
} /* pss1_ */

#undef derivs_ref




/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  subroutine pss2  --  Hessian matrix values for PSS  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     "pss2" is a service routine that computes the sparse */
/*     matrix Hessian elements during PSS global optimization */
/*     in Cartesian coordinate space */


/* Subroutine */ int pss2_(char *mode, doublereal *xx, doublereal *h__, 
	integer *hinit, integer *hstop, integer *hindex, doublereal *hdiag, 
	ftnlen mode_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, nvar;
    extern /* Subroutine */ int hessian_(doublereal *, integer *, integer *, 
	    integer *, doublereal *);



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




/*     translate optimization parameters to atomic coordinates */

    /* Parameter adjustments */
    --hdiag;
    --hindex;
    --hstop;
    --hinit;
    --h__;
    --xx;

    /* Function Body */
    if (s_cmp(mode, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	return 0;
    }
    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	atoms_1.x[i__ - 1] = xx[nvar];
	++nvar;
	atoms_1.y[i__ - 1] = xx[nvar];
	++nvar;
	atoms_1.z__[i__ - 1] = xx[nvar];
    }

/*     compute and store the Hessian elements */

    hessian_(&h__[1], &hinit[1], &hstop[1], &hindex[1], &hdiag[1]);
    return 0;
} /* pss2_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine modecart  --  Cartesian local search for PSS  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/* Subroutine */ int modecart_(integer *start, integer *stop, doublereal *
	minimum, doublereal *grdmin, logical *check)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Cartesian Mode Search :\002,5x,\002Itera"
	    "tion\002,i4,6x,\002Energy\002,f12.4,/)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int climbxyz_(integer *, doublereal *, doublereal 
	    *, doublereal *, logical *), eigenxyz_(doublereal *, doublereal *)
	    ;
    static integer i__, j, k;
    static doublereal eps, rms;
    static logical done;
    static doublereal size, step[3000]	/* was [3][1000] */, eigen[1000];
    static integer niter;
    static doublereal vects[1000000]	/* was [1000][1000] */, xbest[25000], 
	    ybest[25000], zbest[25000];
    extern /* Subroutine */ int getref_(integer *);
    static doublereal minref;
    extern /* Subroutine */ int impose_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *), makeref_(integer *);
    static integer nsearch;
    static doublereal minbest;

    /* Fortran I/O blocks */
    static cilist io___55 = { 0, 0, 0, fmt_10, 0 };



#define step_ref(a_1,a_2) step[(a_2)*3 + a_1 - 4]
#define vects_ref(a_1,a_2) vects[(a_2)*1000 + a_1 - 1001]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  refer.i  --  storage of reference atomic coordinate set  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xref        reference x-coordinates for atoms in each system */
/*     yref        reference y-coordinates for atoms in each system */
/*     zref        reference z-coordinates for atoms in each system */
/*     nref        total number of atoms in each reference system */
/*     reftyp      atom types of the atoms in each reference system */
/*     n12ref      number of atoms bonded to each reference atom */
/*     i12ref      atom numbers of atoms 1-2 connected to each atom */
/*     refleng     length in characters of each reference filename */
/*     refltitle   length in characters of each reference title line */
/*     refnam      atom names of the atoms in each reference system */
/*     reffile     base filename for each reference system */
/*     reftitle    title used to describe each reference system */




/*     store the current coordinates as the reference set */

    makeref_(&c__1);

/*     set parameters related to the local search procedure */

    done = FALSE_;
    eps = 1e-4;
    minref = *minimum;
    minbest = *minimum;
    niter = 0;

/*     find local minimum along each of the steepest directions */

    while(! done) {
	++niter;
	io___55.ciunit = iounit_1.iout;
	s_wsfe(&io___55);
	do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&minref, (ftnlen)sizeof(doublereal));
	e_wsfe();
	eigenxyz_(eigen, vects);

/*     search both directions along each eigenvector in turn */

	nsearch = 0;
	i__1 = *stop;
	for (i__ = *start; i__ <= i__1; ++i__) {
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		j = (k - 1) * 3;
		size = 1. / sqrt((d__1 = eigen[atoms_1.n * 3 - i__], abs(d__1)
			));
		step_ref(1, k) = size * vects_ref(j + 1, atoms_1.n * 3 - i__ 
			+ 1);
		step_ref(2, k) = size * vects_ref(j + 2, atoms_1.n * 3 - i__ 
			+ 1);
		step_ref(3, k) = size * vects_ref(j + 3, atoms_1.n * 3 - i__ 
			+ 1);
	    }
	    ++nsearch;
	    getref_(&c__1);
	    climbxyz_(&nsearch, minimum, step, grdmin, check);
	    if (*minimum < minbest) {
		minbest = *minimum;
		i__2 = atoms_1.n;
		for (k = 1; k <= i__2; ++k) {
		    xbest[k - 1] = atoms_1.x[k - 1];
		    ybest[k - 1] = atoms_1.y[k - 1];
		    zbest[k - 1] = atoms_1.z__[k - 1];
		}
	    }
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		step_ref(1, k) = -step_ref(1, k);
		step_ref(2, k) = -step_ref(2, k);
		step_ref(3, k) = -step_ref(3, k);
	    }
	    ++nsearch;
	    getref_(&c__1);
	    climbxyz_(&nsearch, minimum, step, grdmin, check);
	    if (*minimum < minbest) {
		minbest = *minimum;
		i__2 = atoms_1.n;
		for (k = 1; k <= i__2; ++k) {
		    xbest[k - 1] = atoms_1.x[k - 1];
		    ybest[k - 1] = atoms_1.y[k - 1];
		    zbest[k - 1] = atoms_1.z__[k - 1];
		}
	    }
	}

/*     check for convergence of the local search procedure */

	if (minbest < minref - eps) {
	    done = FALSE_;
	    minref = minbest;
	    impose_(&atoms_1.n, refer_1.xref, refer_1.yref, refer_1.zref, &
		    atoms_1.n, xbest, ybest, zbest, &rms);
	    i__1 = atoms_1.n;
	    for (k = 1; k <= i__1; ++k) {
		atoms_1.x[k - 1] = xbest[k - 1];
		atoms_1.y[k - 1] = ybest[k - 1];
		atoms_1.z__[k - 1] = zbest[k - 1];
	    }
	    makeref_(&c__1);
	} else {
	    done = TRUE_;
	    *minimum = minref;
	    getref_(&c__1);
	}
    }
    return 0;
} /* modecart_ */

#undef vects_ref
#undef step_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine modetors  --  torsional local search for PSS  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/* Subroutine */ int modetors_(integer *start, integer *stop, doublereal *
	minimum, doublereal *grdmin, logical *check)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Torsional Mode Search :\002,5x,\002Itera"
	    "tion\002,i4,6x,\002Energy\002,f12.4,/)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int eigentor_(doublereal *, doublereal *), 
	    climbtor_(integer *, doublereal *, doublereal *, doublereal *, 
	    logical *);
    static integer i__, k;
    static doublereal eps, rms;
    static logical done;
    static doublereal step[1000], eigen[1000];
    static integer niter;
    static doublereal vects[1000000]	/* was [1000][1000] */, xbest[25000], 
	    ybest[25000], zbest[25000];
    extern /* Subroutine */ int getref_(integer *);
    static doublereal minref;
    extern /* Subroutine */ int impose_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *), makeref_(integer *);
    static integer nsearch;
    extern /* Subroutine */ int makeint_(integer *);
    static doublereal minbest;

    /* Fortran I/O blocks */
    static cilist io___73 = { 0, 0, 0, fmt_10, 0 };



#define vects_ref(a_1,a_2) vects[(a_2)*1000 + a_1 - 1001]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  refer.i  --  storage of reference atomic coordinate set  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xref        reference x-coordinates for atoms in each system */
/*     yref        reference y-coordinates for atoms in each system */
/*     zref        reference z-coordinates for atoms in each system */
/*     nref        total number of atoms in each reference system */
/*     reftyp      atom types of the atoms in each reference system */
/*     n12ref      number of atoms bonded to each reference atom */
/*     i12ref      atom numbers of atoms 1-2 connected to each atom */
/*     refleng     length in characters of each reference filename */
/*     refltitle   length in characters of each reference title line */
/*     refnam      atom names of the atoms in each reference system */
/*     reffile     base filename for each reference system */
/*     reftitle    title used to describe each reference system */




/*     store the current coordinates as the reference set */

    makeref_(&c__1);

/*     set parameters related to the local search procedure */

    done = FALSE_;
    eps = 1e-4;
    minref = *minimum;
    minbest = *minimum;
    niter = 0;

/*     find local minimum along each of the steepest directions */

    while(! done) {
	++niter;
	io___73.ciunit = iounit_1.iout;
	s_wsfe(&io___73);
	do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&minref, (ftnlen)sizeof(doublereal));
	e_wsfe();
	makeint_(&c__0);
	eigentor_(eigen, vects);

/*     search both directions along each eigenvector in turn */

	nsearch = 0;
	i__1 = *stop;
	for (i__ = *start; i__ <= i__1; ++i__) {
	    i__2 = omega_1.nomega;
	    for (k = 1; k <= i__2; ++k) {
		step[k - 1] = vects_ref(k, omega_1.nomega - i__ + 1);
	    }
	    ++nsearch;
	    climbtor_(&nsearch, minimum, step, grdmin, check);
	    if (*minimum < minbest) {
		minbest = *minimum;
		i__2 = atoms_1.n;
		for (k = 1; k <= i__2; ++k) {
		    xbest[k - 1] = atoms_1.x[k - 1];
		    ybest[k - 1] = atoms_1.y[k - 1];
		    zbest[k - 1] = atoms_1.z__[k - 1];
		}
	    }
	    i__2 = omega_1.nomega;
	    for (k = 1; k <= i__2; ++k) {
		step[k - 1] = -step[k - 1];
	    }
	    ++nsearch;
	    climbtor_(&nsearch, minimum, step, grdmin, check);
	    if (*minimum < minbest) {
		minbest = *minimum;
		i__2 = atoms_1.n;
		for (k = 1; k <= i__2; ++k) {
		    xbest[k - 1] = atoms_1.x[k - 1];
		    ybest[k - 1] = atoms_1.y[k - 1];
		    zbest[k - 1] = atoms_1.z__[k - 1];
		}
	    }
	}

/*     check for convergence of the local search procedure */

	if (minbest < minref - eps) {
	    done = FALSE_;
	    minref = minbest;
	    impose_(&atoms_1.n, refer_1.xref, refer_1.yref, refer_1.zref, &
		    atoms_1.n, xbest, ybest, zbest, &rms);
	    i__1 = atoms_1.n;
	    for (k = 1; k <= i__1; ++k) {
		atoms_1.x[k - 1] = xbest[k - 1];
		atoms_1.y[k - 1] = ybest[k - 1];
		atoms_1.z__[k - 1] = zbest[k - 1];
	    }
	    makeref_(&c__1);
	} else {
	    done = TRUE_;
	    *minimum = minref;
	    getref_(&c__1);
	}
    }
    return 0;
} /* modetors_ */

#undef vects_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine eigenxyz  --  Cartesian Hessian eigenvectors  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/* Subroutine */ int eigenxyz_(doublereal *eigen, doublereal *vects)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal a[1000], b[1000], h__[1000000];
    static integer i__, j, k;
    static doublereal p[1000], w[1000], ta[1000], tb[1000], ty[1000], hdiag[
	    75000]	/* was [3][25000] */;
    extern /* Subroutine */ int diagq_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer ihess, nfreq, hinit[75000]	/* was [3][25000] */, hstop[
	    75000]	/* was [3][25000] */, hindex[1000000];
    static doublereal matrix[500500];
    extern /* Subroutine */ int hessian_(doublereal *, integer *, integer *, 
	    integer *, doublereal *);


#define hdiag_ref(a_1,a_2) hdiag[(a_2)*3 + a_1 - 4]
#define hinit_ref(a_1,a_2) hinit[(a_2)*3 + a_1 - 4]
#define hstop_ref(a_1,a_2) hstop[(a_2)*3 + a_1 - 4]



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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  hescut.i  --  cutoff value for Hessian matrix elements  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     hesscut   magnitude of smallest allowed Hessian element */




/*     compute the Hessian matrix in Cartesian space */

    /* Parameter adjustments */
    vects -= 1001;
    --eigen;

    /* Function Body */
    hescut_1.hesscut = 0.;
    hessian_(h__, hinit, hstop, hindex, hdiag);

/*     place Hessian elements into triangular form */

    ihess = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    ++ihess;
	    matrix[ihess - 1] = hdiag_ref(j, i__);
	    i__2 = hstop_ref(j, i__);
	    for (k = hinit_ref(j, i__); k <= i__2; ++k) {
		++ihess;
		matrix[ihess - 1] = h__[k - 1];
	    }
	}
    }

/*     diagonalize the Hessian to obtain eigenvalues */

    nfreq = atoms_1.n * 3;
    diagq_(&nfreq, &c__1000, &nfreq, matrix, &eigen[1], &vects[1001], a, b, p,
	     w, ta, tb, ty);
    return 0;
} /* eigenxyz_ */

#undef hstop_ref
#undef hinit_ref
#undef hdiag_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine eigentor  --  torsional Hessian eigenvectors  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/* Subroutine */ int eigentor_(doublereal *eigen, doublereal *vects)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal a[1000], b[1000];
    static integer i__, j;
    static doublereal p[1000], w[1000], ta[1000], tb[1000], ty[1000], hrot[
	    1000000]	/* was [1000][1000] */;
    extern /* Subroutine */ int diagq_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer ihess;
    static doublereal matrix[500500];
    extern /* Subroutine */ int hessrot_(char *, doublereal *, ftnlen);


#define hrot_ref(a_1,a_2) hrot[(a_2)*1000 + a_1 - 1001]



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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




/*     compute the Hessian in torsional space */

    /* Parameter adjustments */
    vects -= 1001;
    --eigen;

    /* Function Body */
    hessrot_("FULL", hrot, (ftnlen)4);

/*     place Hessian elements into triangular form */

    ihess = 0;
    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = omega_1.nomega;
	for (j = i__; j <= i__2; ++j) {
	    ++ihess;
	    matrix[ihess - 1] = hrot_ref(i__, j);
	}
    }

/*     diagonalize the Hessian to obtain eigenvalues */

    diagq_(&omega_1.nomega, &c__1000, &omega_1.nomega, matrix, &eigen[1], &
	    vects[1001], a, b, p, w, ta, tb, ty);
    return 0;
} /* eigentor_ */

#undef hrot_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine climbxyz  --  Cartesian local search direction  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/* Subroutine */ int climbxyz_(integer *nsearch, doublereal *minimum, 
	doublereal *step, doublereal *grdmin, logical *check)
{
    /* Format strings */
    static char fmt_10[] = "(4x,\002Search Direction\002,i4,10x,\002Step\002"
	    ",i6,10x,2f12.4)";
    static char fmt_20[] = "(4x,\002Search Direction\002,i4,10x,\002Step\002"
	    ",i6,10x,f12.4)";
    static char fmt_30[] = "(4x,\002Search Direction\002,i4,36x,\002-----"
	    "-\002)";
    static char fmt_40[] = "(4x,\002Search Direction\002,i4,36x,\002-----"
	    "-\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int localxyz_(doublereal *, doublereal *);
    static integer i__;
    static doublereal big;
    static logical keep, done;
    static doublereal estep[501];
    static integer kstep, nstep;
    extern /* Subroutine */ int getref_(integer *);
    static doublereal parent;
    extern doublereal energy_(void);
    extern /* Subroutine */ int chktree_(doublereal *, doublereal *, logical *
	    );

    /* Fortran I/O blocks */
    static cilist io___122 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___123 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___124 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___125 = { 0, 0, 0, fmt_40, 0 };



#define step_ref(a_1,a_2) step[(a_2)*3 + a_1]



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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  refer.i  --  storage of reference atomic coordinate set  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xref        reference x-coordinates for atoms in each system */
/*     yref        reference y-coordinates for atoms in each system */
/*     zref        reference z-coordinates for atoms in each system */
/*     nref        total number of atoms in each reference system */
/*     reftyp      atom types of the atoms in each reference system */
/*     n12ref      number of atoms bonded to each reference atom */
/*     i12ref      atom numbers of atoms 1-2 connected to each atom */
/*     refleng     length in characters of each reference filename */
/*     refltitle   length in characters of each reference title line */
/*     refnam      atom names of the atoms in each reference system */
/*     reffile     base filename for each reference system */
/*     reftitle    title used to describe each reference system */




/*     convert current reference coordinates to a Z-matrix */

    /* Parameter adjustments */
    step -= 4;

    /* Function Body */
    getref_(&c__1);

/*     set the maximum number of steps and the step size */

    done = FALSE_;
    keep = TRUE_;
    inform_1.iprint = 0;
    big = 1e6;
    *minimum = big;
    kstep = 0;
    nstep = 65;

/*     scan the search direction for a minimization candidate */

    while(! done) {
	if (kstep != 0) {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		atoms_1.x[i__ - 1] += step_ref(1, i__);
		atoms_1.y[i__ - 1] += step_ref(2, i__);
		atoms_1.z__[i__ - 1] += step_ref(3, i__);
	    }
	}
	estep[kstep] = energy_();
	if (kstep >= 2) {
	    if (estep[kstep] < estep[kstep - 2] && estep[kstep - 1] < estep[
		    kstep - 2]) {
		done = TRUE_;
		i__1 = atoms_1.n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    atoms_1.x[i__ - 1] -= step_ref(1, i__);
		    atoms_1.y[i__ - 1] -= step_ref(2, i__);
		    atoms_1.z__[i__ - 1] -= step_ref(3, i__);
		}
		localxyz_(minimum, grdmin);
		parent = *minimum;
		if (*check) {
		    chktree_(&parent, grdmin, &keep);
		}
		if (*minimum >= -big) {
		    if (*check) {
			io___122.ciunit = iounit_1.iout;
			s_wsfe(&io___122);
			do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(
				integer));
			i__1 = kstep - 1;
			do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&(*minimum), (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&parent, (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    } else {
			io___123.ciunit = iounit_1.iout;
			s_wsfe(&io___123);
			do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(
				integer));
			i__1 = kstep - 1;
			do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&(*minimum), (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    }
		} else {
		    *minimum = big;
		    io___124.ciunit = iounit_1.iout;
		    s_wsfe(&io___124);
		    do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(integer)
			    );
		    e_wsfe();
		}
		if (! keep) {
		    *minimum = big;
		}
	    }
	}
	if (kstep >= nstep && ! done) {
	    done = TRUE_;
	    io___125.ciunit = iounit_1.iout;
	    s_wsfe(&io___125);
	    do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	++kstep;
    }
    return 0;
} /* climbxyz_ */

#undef step_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine climbtor  --  torsional local search direction  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/* Subroutine */ int climbtor_(integer *nsearch, doublereal *minimum, 
	doublereal *step, doublereal *grdmin, logical *check)
{
    /* Format strings */
    static char fmt_10[] = "(4x,\002Search Direction\002,i4,10x,\002Step\002"
	    ",i6,10x,2f12.4)";
    static char fmt_20[] = "(4x,\002Search Direction\002,i4,10x,\002Step\002"
	    ",i6,10x,f12.4)";
    static char fmt_30[] = "(4x,\002Search Direction\002,i4,36x,\002-----"
	    "-\002)";
    static char fmt_40[] = "(4x,\002Search Direction\002,i4,36x,\002-----"
	    "-\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int localxyz_(doublereal *, doublereal *);
    static integer i__;
    static doublereal big;
    static logical keep, done;
    static doublereal size, estep[501];
    static integer kstep, nstep;
    extern /* Subroutine */ int getref_(integer *);
    static doublereal parent;
    extern doublereal energy_(void);
    extern /* Subroutine */ int chktree_(doublereal *, doublereal *, logical *
	    ), makeint_(integer *), makexyz_(void);

    /* Fortran I/O blocks */
    static cilist io___135 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___136 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___137 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___138 = { 0, 0, 0, fmt_40, 0 };




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  math.i  --  mathematical and geometrical constants  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     radian   conversion factor from radians to degrees */
/*     pi       numerical value of the geometric constant */
/*     sqrtpi   numerical value of the square root of Pi */
/*     logten   numerical value of the natural log of ten */
/*     sqrttwo  numerical value of the square root of two */
/*     twosix   numerical value of the sixth root of two */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     zbond   bond length used to define each Z-matrix atom */
/*     zang    bond angle used to define each Z-matrix atom */
/*     ztors   angle or torsion used to define Z-matrix atom */
/*     iz      defining atom numbers for each Z-matrix atom */




/*     convert current reference coordinates to a Z-matrix */

    /* Parameter adjustments */
    --step;

    /* Function Body */
    getref_(&c__1);
    makeint_(&c__0);

/*     set the maximum number of steps and the step size */

    done = FALSE_;
    keep = TRUE_;
    inform_1.iprint = 0;
    big = 1e6;
    *minimum = big;
    kstep = 0;
    nstep = 65;
    size = 5.729577951308233;
    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	step[i__] = size * step[i__];
    }

/*     scan the search direction for a minimization candidate */

    while(! done) {
	if (kstep != 0) {
	    i__1 = omega_1.nomega;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] += step[i__];
	    }
	}
	makexyz_();
	estep[kstep] = energy_();
	if (kstep >= 2) {
	    if (estep[kstep] < estep[kstep - 2] && estep[kstep - 1] < estep[
		    kstep - 2]) {
		done = TRUE_;
		i__1 = omega_1.nomega;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] -= step[i__];
		}
		makexyz_();
		localxyz_(minimum, grdmin);
		parent = *minimum;
		if (*check) {
		    chktree_(&parent, grdmin, &keep);
		}
		if (*minimum >= -big) {
		    if (*check) {
			io___135.ciunit = iounit_1.iout;
			s_wsfe(&io___135);
			do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(
				integer));
			i__1 = kstep - 1;
			do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&(*minimum), (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&parent, (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    } else {
			io___136.ciunit = iounit_1.iout;
			s_wsfe(&io___136);
			do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(
				integer));
			i__1 = kstep - 1;
			do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&(*minimum), (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    }
		} else {
		    *minimum = big;
		    io___137.ciunit = iounit_1.iout;
		    s_wsfe(&io___137);
		    do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(integer)
			    );
		    e_wsfe();
		}
		if (! keep) {
		    *minimum = big;
		}
	    }
	}
	if (kstep >= nstep && ! done) {
	    done = TRUE_;
	    io___138.ciunit = iounit_1.iout;
	    s_wsfe(&io___138);
	    do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	++kstep;
    }
    return 0;
} /* climbtor_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine localxyz  --  PSS local search optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "localxyz" is used during the potential smoothing and search */
/*     procedure to perform a local optimization at the current */
/*     smoothing level */


/* Subroutine */ int localxyz_(doublereal *minimum, doublereal *grdmin)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal xx[75000];
    extern doublereal pss1_(doublereal *, doublereal *);
    extern /* Subroutine */ int pss2_(char *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, ftnlen);
    static char mode[6];
    extern /* Subroutine */ int tncg_(char *, char *, integer *, doublereal *,
	     doublereal *, doublereal *, D_fp, S_fp, U_fp, ftnlen, ftnlen);
    static integer nvar;
    static char method[6];
    static logical oldverb;
    extern /* Subroutine */ int optsave_();



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




/*     translate the coordinates of each atom */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	xx[nvar - 1] = atoms_1.x[i__ - 1];
	++nvar;
	xx[nvar - 1] = atoms_1.y[i__ - 1];
	++nvar;
	xx[nvar - 1] = atoms_1.z__[i__ - 1];
    }

/*     make the call to the optimization routine */

    oldverb = inform_1.verbose;
    inform_1.verbose = FALSE_;
    s_copy(mode, "AUTO", (ftnlen)6, (ftnlen)4);
    s_copy(method, "AUTO", (ftnlen)6, (ftnlen)4);
    tncg_(mode, method, &nvar, xx, minimum, grdmin, (D_fp)pss1_, (S_fp)pss2_, 
	    (U_fp)optsave_, (ftnlen)6, (ftnlen)6);
    inform_1.verbose = oldverb;

/*     untranslate the final coordinates for each atom */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	atoms_1.x[i__ - 1] = xx[nvar - 1];
	++nvar;
	atoms_1.y[i__ - 1] = xx[nvar - 1];
	++nvar;
	atoms_1.z__[i__ - 1] = xx[nvar - 1];
    }
    return 0;
} /* localxyz_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine chktree  --  check for legitimacy of branch  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "chktree" tests a minimum energy structure to see if it */
/*     belongs to the correct progenitor in the existing map */


/* Subroutine */ int chktree_(doublereal *parent, doublereal *grdmin, logical 
	*keep)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ int localxyz_(doublereal *, doublereal *);
    static integer i__;
    static doublereal x0[25000], y0[25000], z0[25000], eps, deform0;



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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  tree.i  --  potential smoothing and search tree levels  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     maxpss   maximum number of potential smoothing levels */

/*     etree    energy reference value at the top of the tree */
/*     ilevel   smoothing deformation value at each tree level */
/*     nlevel   number of levels of potential smoothing used */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  warp.i  --  parameters for potential surface smoothing  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     m2           second moment of the GDA gaussian for each atom */
/*     deform       value of the smoothing deformation parameter */
/*     difft        diffusion coefficient for torsional potential */
/*     diffv        diffusion coefficient for van der Waals potential */
/*     diffc        diffusion coefficient for charge-charge potential */
/*     use_smooth   flag to use a potential energy smoothing method */
/*     use_dem      flag to use diffusion equation method potential */
/*     use_gda      flag to use gaussian density annealing potential */
/*     use_tophat   flag to use analytical tophat smoothed potential */
/*     use_stophat  flag to use shifted tophat smoothed potential */



/*     store the current smoothing level and coordinates */

    deform0 = warp_1.deform;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x0[i__ - 1] = atoms_1.x[i__ - 1];
	y0[i__ - 1] = atoms_1.y[i__ - 1];
	z0[i__ - 1] = atoms_1.z__[i__ - 1];
    }

/*     forward smoothing optimizations back to highest level */

    i__1 = tree_1.nlevel;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (warp_1.deform < tree_1.ilevel[i__]) {
	    warp_1.deform = tree_1.ilevel[i__];
	    localxyz_(parent, grdmin);
	}
    }

/*     compare energy to reference value for this tree branch */

    eps = 1e-4;
    *keep = FALSE_;
    if ((d__1 = *parent - tree_1.etree, abs(d__1)) < eps) {
	*keep = TRUE_;
    }

/*     restore the original smoothing level and coordinates */

    warp_1.deform = deform0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = x0[i__ - 1];
	atoms_1.y[i__ - 1] = y0[i__ - 1];
	atoms_1.z__[i__ - 1] = z0[i__ - 1];
    }
    return 0;
} /* chktree_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine psswrite  --  output structures on PSS path  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/* Subroutine */ int psswrite_(integer *i__)
{
    /* System generated locals */
    address a__1[3];
    integer i__1[3];
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern integer freeunit_(void);
    static char ext[7];
    static integer lext, ixyz;
    extern /* Subroutine */ int prtxyz_(integer *), numeral_(integer *, char *
	    , integer *, ftnlen), version_(char *, char *, ftnlen, ftnlen);
    static char xyzfile[120];



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




/*     write the coordinates of the current minimum to a file */

    lext = 3;
    numeral_(i__, ext, &lext, (ftnlen)7);
    ixyz = freeunit_();
/* Writing concatenation */
    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
    i__1[1] = 1, a__1[1] = ".";
    i__1[2] = lext, a__1[2] = ext;
    s_cat(xyzfile, a__1, i__1, &c__3, (ftnlen)120);
    version_(xyzfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = ixyz;
    o__1.ofnmlen = 120;
    o__1.ofnm = xyzfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    prtxyz_(&ixyz);
    cl__1.cerr = 0;
    cl__1.cunit = ixyz;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* psswrite_ */

/* Main program alias */ int pss_ () { MAIN__ (); return 0; }
