/* testhess.f -- translated by f2c (version 20050501).
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
    doublereal hesscut;
} hescut_;

#define hescut_1 hescut_

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
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__1 = 1;
static integer c__5 = 5;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  program testhess  --  Hessian matrix test; cart. version  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "testhess" computes and compares the analytical and numerical */
/*     Hessian matrices of the potential energy function with respect */
/*     to Cartesian coordinates */


/* Main program */ int MAIN__(void)
{
    /* Initialized data */

    static char axis[1*3] = "X" "Y" "Z";

    /* Format strings */
    static char fmt_10[] = "(/,\002 Compute Analytical Hessian Matrix [Y] "
	    ":  \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(/,\002 Compute Numerical Hessian Matrix [Y] :"
	    "   \002,$)";
    static char fmt_40[] = "(a120)";
    static char fmt_50[] = "(/,\002 Numerical Hessian from Gradient\002,\002"
	    " or Function [G] :  \002,$)";
    static char fmt_60[] = "(a120)";
    static char fmt_80[] = "(/,\002 Enter a Numerical Stepsize [\002,d7.1"
	    ",\002 Ang] :  \002,$)";
    static char fmt_90[] = "(f20.0)";
    static char fmt_100[] = "(/,\002 List Individual Hessian Components [N] "
	    ":   \002,$)";
    static char fmt_110[] = "(a120)";
    static char fmt_120[] = "(/,\002 Comparison of Analytical and\002,\002 N"
	    "umerical Hessian Elements :\002,//,3x,\0021st Atom\002,4x,\0022n"
	    "d Atom\002,10x,\002Analytical\002,8x,\002Numerical\002,6x,\002Di"
	    "fference\002,/)";
    static char fmt_130[] = "(1x,i6,\002 (\002,a1,\002) \002,1x,i6,\002 ("
	    "\002,a1,\002) \002,2x,2f17.8,f16.8)";
    static char fmt_140[] = "(1x,i6,\002 (\002,a1,\002) \002,1x,i6,\002 ("
	    "\002,a1,\002) \002,2x,2f17.6,f16.6)";
    static char fmt_150[] = "(1x,i6,\002 (\002,a1,\002) \002,1x,i6,\002 ("
	    "\002,a1,\002) \002,2x,2f17.4,f16.4)";
    static char fmt_160[] = "(/,\002 Comparison of Analytical and\002,\002 N"
	    "umerical Hessian Elements :\002,//,3x,\0021st Atom\002,4x,\0022n"
	    "d Atom\002,10x,\002Analytical\002,8x,\002Numerical\002,6x,\002Di"
	    "fference\002,/)";
    static char fmt_170[] = "(1x,i6,\002 (\002,a1,\002) \002,1x,i6,\002 ("
	    "\002,a1,\002) \002,2x,2f17.8,f16.8)";
    static char fmt_180[] = "(1x,i6,\002 (\002,a1,\002) \002,1x,i6,\002 ("
	    "\002,a1,\002) \002,2x,2f17.6,f16.6)";
    static char fmt_190[] = "(1x,i6,\002 (\002,a1,\002) \002,1x,i6,\002 ("
	    "\002,a1,\002) \002,2x,2f17.4,f16.4)";
    static char fmt_200[] = "(/,\002 Analytical and Numerical Hessian Elemen"
	    "ts\002,\002 are Identical\002)";
    static char fmt_210[] = "(/,\002 Diagonal Hessian Elements for Each Atom"
	    " :\002,//,6x,\002Atom\002,21x,\002X\002,19x,\002Y\002,19x,\002"
	    "Z\002,/)";
    static char fmt_220[] = "(/,\002 Diagonal Hessian Elements for Each Atom"
	    " :\002,//,6x,\002Atom\002,19x,\002X\002,17x,\002Y\002,17x,\002"
	    "Z\002,/)";
    static char fmt_230[] = "(/,\002 Diagonal Hessian Elements for Each Atom"
	    " :\002,//,6x,\002Atom\002,17x,\002X\002,15x,\002Y\002,15x,\002"
	    "Z\002,/)";
    static char fmt_240[] = "(i10,5x,3f20.8)";
    static char fmt_250[] = "(i10,5x,3f18.6)";
    static char fmt_260[] = "(i10,5x,3f16.4)";
    static char fmt_270[] = "(/,\002 Sum of Diagonal Hessian Elements :\002,"
	    "6x,f20.8)";
    static char fmt_280[] = "(/,\002 Sum of Diagonal Hessian Elements :\002,"
	    "6x,f18.6)";
    static char fmt_290[] = "(/,\002 Sum of Diagonal Hessian Elements :\002,"
	    "6x,f16.4)";
    static char fmt_300[] = "(/,\002 3x3 Hessian Block for Atoms :\002,3x,2i"
	    "8,/)";
    static char fmt_310[] = "(\002 Numer\002,5x,3f20.8)";
    static char fmt_320[] = "(\002 Numer\002,5x,3f18.6)";
    static char fmt_330[] = "(\002 Numer\002,5x,3f16.4)";
    static char fmt_340[] = "(/,\002 Hessian Matrix written to File :  \002,"
	    "a40)";
    static char fmt_350[] = "(/,\002 Diagonal Hessian Elements  (3 per Atom"
	    ")\002,/)";
    static char fmt_360[] = "(4f16.8)";
    static char fmt_370[] = "(5f14.6)";
    static char fmt_380[] = "(6f12.4)";
    static char fmt_390[] = "(/,\002 Off-diagonal Hessian Elements for Ato"
	    "m\002,i6,1x,a1,/)";
    static char fmt_400[] = "(4f16.8)";
    static char fmt_410[] = "(5f14.6)";
    static char fmt_420[] = "(6f12.4)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2];
    doublereal d__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void), gradient_(doublereal *, 
	    doublereal *);
    static char hessfile[120];
    static logical doanalyt;
    extern integer freeunit_(void);
    static doublereal e, g[75000]	/* was [3][25000] */, h__[1000000];
    static integer i__, j, k, m;
    static logical identical;
    static doublereal g0[75000]	/* was [3][25000] */;
    static integer ii, jj;
    static doublereal old, eps, sum, eps0, diff;
    static integer ihes, next;
    static doublereal hdiag[75000]	/* was [3][25000] */, delta;
    extern /* Subroutine */ int final_(void);
    static integer index, hinit[75000]	/* was [3][25000] */;
    static doublereal nhess[810000]	/* was [3][300][3][300] */;
    static logical exist;
    static integer hstop[75000]	/* was [3][25000] */;
    static logical query, dograd;
    static char record[120];
    static integer hindex[1000000];
    extern doublereal energy_();
    static logical dofull;
    static char answer[1], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), getxyz_(void), 
	    initial_(void), hessian_(doublereal *, integer *, integer *, 
	    integer *, doublereal *), numgrad_(D_fp, doublereal *, doublereal 
	    *);
    static logical donumer;
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), gettext_(
	    char *, char *, integer *, ftnlen, ftnlen), version_(char *, char 
	    *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_60, 0 };
    static icilist io___20 = { 1, string, 1, 0, 120, 1 };
    static cilist io___21 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___22 = { 1, 0, 0, fmt_90, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___48 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___49 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___52 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___54 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___55 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___56 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___57 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___58 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___60 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___64 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___66 = { 0, 0, 0, fmt_320, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_330, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_340, 0 };
    static cilist io___71 = { 0, 0, 0, fmt_350, 0 };
    static cilist io___72 = { 0, 0, 0, fmt_360, 0 };
    static cilist io___73 = { 0, 0, 0, fmt_370, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_380, 0 };
    static cilist io___75 = { 0, 0, 0, fmt_390, 0 };
    static cilist io___76 = { 0, 0, 0, fmt_400, 0 };
    static cilist io___77 = { 0, 0, 0, fmt_410, 0 };
    static cilist io___78 = { 0, 0, 0, fmt_420, 0 };



#define g_ref(a_1,a_2) g[(a_2)*3 + a_1 - 4]
#define g0_ref(a_1,a_2) g0[(a_2)*3 + a_1 - 4]
#define hdiag_ref(a_1,a_2) hdiag[(a_2)*3 + a_1 - 4]
#define hinit_ref(a_1,a_2) hinit[(a_2)*3 + a_1 - 4]
#define nhess_ref(a_1,a_2,a_3,a_4) nhess[(((a_4)*3 + (a_3))*300 + (a_2))*\
3 + a_1 - 3604]
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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  hescut.i  --  cutoff value for Hessian matrix elements  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     hesscut   magnitude of smallest allowed Hessian element */




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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     set up the structure and mechanics calculation */

    initial_();
    getxyz_();
    mechanic_();

/*     set difference threshhold via the energy precision */

    delta = 1e-4;
    if (inform_1.digits >= 6) {
	delta = 1e-6;
    }
    if (inform_1.digits >= 8) {
	delta = 1e-8;
    }

/*     decide whether to do an analytical Hessian calculation */

    doanalyt = TRUE_;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___6.ciunit = iounit_1.iout;
	s_wsfe(&io___6);
	e_wsfe();
	io___7.ciunit = iounit_1.input;
	s_rsfe(&io___7);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'N') {
	doanalyt = FALSE_;
    }

/*     decide whether to do a numerical Hessian calculation */

    donumer = FALSE_;
    if (atoms_1.n <= 300) {
	donumer = TRUE_;
	nextarg_(answer, &exist, (ftnlen)1);
	if (! exist) {
	    io___11.ciunit = iounit_1.iout;
	    s_wsfe(&io___11);
	    e_wsfe();
	    io___12.ciunit = iounit_1.input;
	    s_rsfe(&io___12);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    next = 1;
	    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
	}
	upcase_(answer, (ftnlen)1);
	if (*(unsigned char *)answer == 'N') {
	    donumer = FALSE_;
	}
    }

/*     get numerical Hessian from either gradient or energy */

    if (donumer) {
	dograd = TRUE_;
	nextarg_(answer, &exist, (ftnlen)1);
	if (! exist) {
	    io___14.ciunit = iounit_1.iout;
	    s_wsfe(&io___14);
	    e_wsfe();
	    io___15.ciunit = iounit_1.input;
	    s_rsfe(&io___15);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    next = 1;
	    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
	}
	upcase_(answer, (ftnlen)1);
	if (*(unsigned char *)answer == 'F') {
	    dograd = FALSE_;
	}

/*     get the stepsize for numerical Hessian calculation */

	eps = -1.;
	eps0 = .001;
	if (dograd) {
	    eps0 = 1e-5;
	}
	query = TRUE_;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___20);
	    if (i__1 != 0) {
		goto L70;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&eps, (ftnlen)sizeof(
		    doublereal));
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
	    io___21.ciunit = iounit_1.iout;
	    s_wsfe(&io___21);
	    do_fio(&c__1, (char *)&eps0, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    io___22.ciunit = iounit_1.input;
	    i__1 = s_rsfe(&io___22);
	    if (i__1 != 0) {
		goto L70;
	    }
	    i__1 = do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
	    if (i__1 != 0) {
		goto L70;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L70;
	    }
	}
	if (eps <= 0.) {
	    eps = eps0;
	}
    }

/*     decide whether to output results by Hessian component */

    dofull = FALSE_;
    if (atoms_1.n <= 20 && donumer) {
	nextarg_(answer, &exist, (ftnlen)1);
	if (! exist) {
	    io___24.ciunit = iounit_1.iout;
	    s_wsfe(&io___24);
	    e_wsfe();
	    io___25.ciunit = iounit_1.input;
	    s_rsfe(&io___25);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    next = 1;
	    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
	}
	upcase_(answer, (ftnlen)1);
	if (*(unsigned char *)answer == 'Y') {
	    dofull = TRUE_;
	}
    }

/*     get the analytical Hessian matrix elements */

    identical = TRUE_;
    if (doanalyt) {
	hescut_1.hesscut = 0.;
	hessian_(h__, hinit, hstop, hindex, hdiag);
	hessian_(h__, hinit, hstop, hindex, hdiag);
	hessian_(h__, hinit, hstop, hindex, hdiag);
    }

/*     get the two-sided numerical Hessian matrix elements */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (donumer && usage_1.use[i__ - 1]) {
	    old = atoms_1.x[i__ - 1];
	    atoms_1.x[i__ - 1] -= eps * .5;
	    if (dograd) {
		gradient_(&e, g);
	    } else {
		numgrad_((D_fp)energy_, g, &eps);
	    }
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		for (j = 1; j <= 3; ++j) {
		    g0_ref(j, k) = g_ref(j, k);
		}
	    }
	    atoms_1.x[i__ - 1] += eps;
	    if (dograd) {
		gradient_(&e, g);
	    } else {
		numgrad_((D_fp)energy_, g, &eps);
	    }
	    atoms_1.x[i__ - 1] = old;
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		for (j = 1; j <= 3; ++j) {
		    nhess_ref(j, k, 1, i__) = (g_ref(j, k) - g0_ref(j, k)) / 
			    eps;
		}
	    }
	    old = atoms_1.y[i__ - 1];
	    atoms_1.y[i__ - 1] -= eps * .5;
	    if (dograd) {
		gradient_(&e, g);
	    } else {
		numgrad_((D_fp)energy_, g, &eps);
	    }
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		for (j = 1; j <= 3; ++j) {
		    g0_ref(j, k) = g_ref(j, k);
		}
	    }
	    atoms_1.y[i__ - 1] += eps;
	    if (dograd) {
		gradient_(&e, g);
	    } else {
		numgrad_((D_fp)energy_, g, &eps);
	    }
	    atoms_1.y[i__ - 1] = old;
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		for (j = 1; j <= 3; ++j) {
		    nhess_ref(j, k, 2, i__) = (g_ref(j, k) - g0_ref(j, k)) / 
			    eps;
		}
	    }
	    old = atoms_1.z__[i__ - 1];
	    atoms_1.z__[i__ - 1] -= eps * .5;
	    if (dograd) {
		gradient_(&e, g);
	    } else {
		numgrad_((D_fp)energy_, g, &eps);
	    }
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		for (j = 1; j <= 3; ++j) {
		    g0_ref(j, k) = g_ref(j, k);
		}
	    }
	    atoms_1.z__[i__ - 1] += eps;
	    if (dograd) {
		gradient_(&e, g);
	    } else {
		numgrad_((D_fp)energy_, g, &eps);
	    }
	    atoms_1.z__[i__ - 1] = old;
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		for (j = 1; j <= 3; ++j) {
		    nhess_ref(j, k, 3, i__) = (g_ref(j, k) - g0_ref(j, k)) / 
			    eps;
		}
	    }
	}

/*     compare the analytical and numerical diagonal elements */

	if (doanalyt && donumer) {
	    for (j = 1; j <= 3; ++j) {
		diff = (d__1 = hdiag_ref(j, i__) - nhess_ref(j, i__, j, i__), 
			abs(d__1));
		if (diff > delta) {
		    if (identical) {
			identical = FALSE_;
			io___41.ciunit = iounit_1.iout;
			s_wsfe(&io___41);
			e_wsfe();
		    }
		    if (inform_1.digits >= 8) {
			io___42.ciunit = iounit_1.iout;
			s_wsfe(&io___42);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			do_fio(&c__1, (char *)&hdiag_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&nhess_ref(j, i__, j, i__), (
				ftnlen)sizeof(doublereal));
			d__1 = hdiag_ref(j, i__) - nhess_ref(j, i__, j, i__);
			do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    } else if (inform_1.digits >= 6) {
			io___43.ciunit = iounit_1.iout;
			s_wsfe(&io___43);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			do_fio(&c__1, (char *)&hdiag_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&nhess_ref(j, i__, j, i__), (
				ftnlen)sizeof(doublereal));
			d__1 = hdiag_ref(j, i__) - nhess_ref(j, i__, j, i__);
			do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    } else {
			io___44.ciunit = iounit_1.iout;
			s_wsfe(&io___44);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			do_fio(&c__1, (char *)&hdiag_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&nhess_ref(j, i__, j, i__), (
				ftnlen)sizeof(doublereal));
			d__1 = hdiag_ref(j, i__) - nhess_ref(j, i__, j, i__);
			do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    }
		}

/*     compare the analytical and numerical off-diagonal elements */

		i__2 = hstop_ref(j, i__);
		for (k = hinit_ref(j, i__); k <= i__2; ++k) {
		    index = hindex[k - 1];
		    jj = index % 3;
		    if (jj == 0) {
			jj = 3;
		    }
		    ii = (index + 2) / 3;
		    diff = (d__1 = h__[k - 1] - nhess_ref(jj, ii, j, i__), 
			    abs(d__1));
		    if (diff > delta) {
			if (identical) {
			    identical = FALSE_;
			    io___48.ciunit = iounit_1.iout;
			    s_wsfe(&io___48);
			    e_wsfe();
			}
			if (inform_1.digits >= 8) {
			    io___49.ciunit = iounit_1.iout;
			    s_wsfe(&io___49);
			    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer)
				    );
			    do_fio(&c__1, axis + (jj - 1), (ftnlen)1);
			    do_fio(&c__1, (char *)&h__[k - 1], (ftnlen)sizeof(
				    doublereal));
			    do_fio(&c__1, (char *)&nhess_ref(jj, ii, j, i__), 
				    (ftnlen)sizeof(doublereal));
			    d__1 = h__[k - 1] - nhess_ref(jj, ii, j, i__);
			    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				    doublereal));
			    e_wsfe();
			} else if (inform_1.digits >= 6) {
			    io___50.ciunit = iounit_1.iout;
			    s_wsfe(&io___50);
			    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer)
				    );
			    do_fio(&c__1, axis + (jj - 1), (ftnlen)1);
			    do_fio(&c__1, (char *)&h__[k - 1], (ftnlen)sizeof(
				    doublereal));
			    do_fio(&c__1, (char *)&nhess_ref(jj, ii, j, i__), 
				    (ftnlen)sizeof(doublereal));
			    d__1 = h__[k - 1] - nhess_ref(jj, ii, j, i__);
			    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				    doublereal));
			    e_wsfe();
			} else {
			    io___51.ciunit = iounit_1.iout;
			    s_wsfe(&io___51);
			    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer)
				    );
			    do_fio(&c__1, axis + (jj - 1), (ftnlen)1);
			    do_fio(&c__1, (char *)&h__[k - 1], (ftnlen)sizeof(
				    doublereal));
			    do_fio(&c__1, (char *)&nhess_ref(jj, ii, j, i__), 
				    (ftnlen)sizeof(doublereal));
			    d__1 = h__[k - 1] - nhess_ref(jj, ii, j, i__);
			    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				    doublereal));
			    e_wsfe();
			}
		    }
		}
	    }
	}
    }

/*     success if the analytical and numerical elements are the same */

    if (doanalyt && donumer) {
	if (identical) {
	    io___52.ciunit = iounit_1.iout;
	    s_wsfe(&io___52);
	    e_wsfe();
	}
    }

/*     write out the diagonal Hessian elements for each atom */

    if (doanalyt) {
	if (inform_1.digits >= 8) {
	    io___53.ciunit = iounit_1.iout;
	    s_wsfe(&io___53);
	    e_wsfe();
	} else if (inform_1.digits >= 6) {
	    io___54.ciunit = iounit_1.iout;
	    s_wsfe(&io___54);
	    e_wsfe();
	} else {
	    io___55.ciunit = iounit_1.iout;
	    s_wsfe(&io___55);
	    e_wsfe();
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (inform_1.digits >= 8) {
		io___56.ciunit = iounit_1.iout;
		s_wsfe(&io___56);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&hdiag_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    } else if (inform_1.digits >= 6) {
		io___57.ciunit = iounit_1.iout;
		s_wsfe(&io___57);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&hdiag_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    } else {
		io___58.ciunit = iounit_1.iout;
		s_wsfe(&io___58);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&hdiag_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	}
    }

/*     write out the Hessian trace as sum of diagonal elements */

    if (doanalyt) {
	sum = 0.;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		sum += hdiag_ref(j, i__);
	    }
	}
	if (inform_1.digits >= 8) {
	    io___60.ciunit = iounit_1.iout;
	    s_wsfe(&io___60);
	    do_fio(&c__1, (char *)&sum, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (inform_1.digits >= 6) {
	    io___61.ciunit = iounit_1.iout;
	    s_wsfe(&io___61);
	    do_fio(&c__1, (char *)&sum, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___62.ciunit = iounit_1.iout;
	    s_wsfe(&io___62);
	    do_fio(&c__1, (char *)&sum, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     write out the full matrix of numerical Hessian elements */

    if (dofull && donumer) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		io___63.ciunit = iounit_1.iout;
		s_wsfe(&io___63);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		e_wsfe();
		for (j = 1; j <= 3; ++j) {
		    if (inform_1.digits >= 8) {
			io___64.ciunit = iounit_1.iout;
			s_wsfe(&io___64);
			for (m = 1; m <= 3; ++m) {
			    do_fio(&c__1, (char *)&nhess_ref(m, i__, j, k), (
				    ftnlen)sizeof(doublereal));
			}
			e_wsfe();
		    } else if (inform_1.digits >= 6) {
			io___66.ciunit = iounit_1.iout;
			s_wsfe(&io___66);
			for (m = 1; m <= 3; ++m) {
			    do_fio(&c__1, (char *)&nhess_ref(m, i__, j, k), (
				    ftnlen)sizeof(doublereal));
			}
			e_wsfe();
		    } else {
			io___67.ciunit = iounit_1.iout;
			s_wsfe(&io___67);
			for (m = 1; m <= 3; ++m) {
			    do_fio(&c__1, (char *)&nhess_ref(m, i__, j, k), (
				    ftnlen)sizeof(doublereal));
			}
			e_wsfe();
		    }
		}
	    }
	}
    }

/*     write out the full matrix of analytical Hessian elements */

    if (doanalyt && ! donumer) {
	ihes = freeunit_();
/* Writing concatenation */
	i__3[0] = files_1.leng, a__1[0] = files_1.filename;
	i__3[1] = 4, a__1[1] = ".hes";
	s_cat(hessfile, a__1, i__3, &c__2, (ftnlen)120);
	version_(hessfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = ihes;
	o__1.ofnmlen = 120;
	o__1.ofnm = hessfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	io___70.ciunit = iounit_1.iout;
	s_wsfe(&io___70);
	do_fio(&c__1, hessfile, (ftnlen)120);
	e_wsfe();
	io___71.ciunit = ihes;
	s_wsfe(&io___71);
	e_wsfe();
	if (inform_1.digits >= 8) {
	    io___72.ciunit = ihes;
	    s_wsfe(&io___72);
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&hdiag_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
	    }
	    e_wsfe();
	} else if (inform_1.digits >= 6) {
	    io___73.ciunit = ihes;
	    s_wsfe(&io___73);
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&hdiag_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
	    }
	    e_wsfe();
	} else {
	    io___74.ciunit = ihes;
	    s_wsfe(&io___74);
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&hdiag_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
	    }
	    e_wsfe();
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		if (hinit_ref(j, i__) <= hstop_ref(j, i__)) {
		    io___75.ciunit = ihes;
		    s_wsfe(&io___75);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, axis + (j - 1), (ftnlen)1);
		    e_wsfe();
		    if (inform_1.digits >= 8) {
			io___76.ciunit = ihes;
			s_wsfe(&io___76);
			i__2 = hstop_ref(j, i__);
			for (k = hinit_ref(j, i__); k <= i__2; ++k) {
			    do_fio(&c__1, (char *)&h__[k - 1], (ftnlen)sizeof(
				    doublereal));
			}
			e_wsfe();
		    } else if (inform_1.digits >= 6) {
			io___77.ciunit = ihes;
			s_wsfe(&io___77);
			i__2 = hstop_ref(j, i__);
			for (k = hinit_ref(j, i__); k <= i__2; ++k) {
			    do_fio(&c__1, (char *)&h__[k - 1], (ftnlen)sizeof(
				    doublereal));
			}
			e_wsfe();
		    } else {
			io___78.ciunit = ihes;
			s_wsfe(&io___78);
			i__2 = hstop_ref(j, i__);
			for (k = hinit_ref(j, i__); k <= i__2; ++k) {
			    do_fio(&c__1, (char *)&h__[k - 1], (ftnlen)sizeof(
				    doublereal));
			}
			e_wsfe();
		    }
		}
	    }
	}
	cl__1.cerr = 0;
	cl__1.cunit = ihes;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef hstop_ref
#undef nhess_ref
#undef hinit_ref
#undef hdiag_ref
#undef g0_ref
#undef g_ref


/* Main program alias */ int testhess_ () { MAIN__ (); return 0; }
