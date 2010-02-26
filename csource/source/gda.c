/* gda.f -- translated by f2c (version 20050501).
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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal fctmin, hguess;
    integer maxiter, nextiter;
} minima_;

#define minima_1 minima_

struct {
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

struct {
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_

struct {
    doublereal m2[25000], deform, difft, diffv, diffc;
    logical use_smooth__, use_dem__, use_gda__, use_tophat__, use_stophat__;
} warp_;

#define warp_1 warp_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static doublereal c_b15 = .33333333333333331;
static integer c__5 = 5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  program gda  --  simulated annealing on gaussian density  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "gda" implements Gaussian Density Annealing (GDA) algorithm */
/*     for global optimization via simulated annealing */

/*     literature reference: */

/*     J. Ma and J. E. Straub, "Simulated Annealing using the */
/*     Classical Density Distribution", Journal of Chemical Physics, */
/*     101, 533-541 (1994) */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Enter Number of Annealing Trials [1] : "
	    " \002,$)";
    static char fmt_30[] = "(i10)";
    static char fmt_40[] = "(/,\002 Use Randomized Initial Coordinates [N] :"
	    "  \002,$)";
    static char fmt_50[] = "(a120)";
    static char fmt_70[] = "(/,\002 Enter Initial and Final Beta [0.01, 10**"
	    "10] :  \002,$)";
    static char fmt_80[] = "(a120)";
    static char fmt_100[] = "(//,\002 Gaussian Density Annealing Global Opti"
	    "mization :\002,//,\002 BS Step\002,5x,\002Log(Beta)\002,6x,\002E"
	    "nergy\002,9x,\002Rg\002,8x,\002Log(M2)\002,7x,\002Status\002,/)";
    static char fmt_110[] = "(/,\002 Global Energy Minimum for Trial\002,i4"
	    ",\002 :\002,f15.4)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2[3], i__3;
    doublereal d__1;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);
    double pow_dd(doublereal *, doublereal *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_rew(alist *);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    extern integer freeunit_(void);
    static integer i__;
    static doublereal h1;
    static logical randomize;
    static doublereal xx[100000];
    static integer nok;
    static doublereal eps, xcm, ycm, zcm;
    static char ext[7];
    extern /* Subroutine */ int gda1_();
    extern doublereal gda2_();
    extern /* Subroutine */ int gda3_();
    static integer igda, nbad;
    static char mode[6];
    static doublereal hmin;
    extern /* Subroutine */ int tncg_(char *, char *, integer *, doublereal *,
	     doublereal *, doublereal *, D_fp, U_fp, U_fp, ftnlen, ftnlen);
    static integer nvar, lext, next;
    extern /* Subroutine */ int final_(void);
    static doublereal bstop;
    static integer nstep;
    static logical exist;
    static doublereal m2init[25000];
    extern /* Subroutine */ int diffeq_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *, U_fp);
    static char record[120];
    static integer itrial, ntrial;
    extern doublereal random_(void);
    static doublereal bstart, grdmin;
    static char answer[1], method[6], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char status[7];
    extern /* Subroutine */ int getxyz_(void), prtxyz_(integer *);
    static char gdafile[120];
    extern /* Subroutine */ int gdastat_(integer *, doublereal *, doublereal *
	    , char *, ftnlen), initial_(void), numeral_(integer *, char *, 
	    integer *, ftnlen), nextarg_(char *, logical *, ftnlen);
    static doublereal minimum;
    extern /* Subroutine */ int optsave_();
    static doublereal boxsize;
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___6 = { 1, string, 1, 0, 120, 1 };
    static cilist io___7 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___18 = { 1, string, 1, 0, 120, 1 };
    static icilist io___19 = { 1, string, 1, 0, 120, 1 };
    static cilist io___20 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_80, 0 };
    static icilist io___22 = { 1, record, 1, 0, 120, 1 };
    static cilist io___35 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_110, 0 };




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
/*     ##  minima.i  --  general parameters for minimizations  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     fctmin    value below which function is deemed optimized */
/*     hguess    initial value for the H-matrix diagonal elements */
/*     maxiter   maximum number of iterations during optimization */
/*     nextiter  iteration number to use for the first iteration */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  potent.i  --  usage of each potential energy component  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     use_bond    logical flag governing use of bond stretch potential */
/*     use_angle   logical flag governing use of angle bend potential */
/*     use_strbnd  logical flag governing use of stretch-bend potential */
/*     use_urey    logical flag governing use of Urey-Bradley potential */
/*     use_angang  logical flag governing use of angle-angle cross term */
/*     use_opbend  logical flag governing use of out-of-plane bend term */
/*     use_opdist  logical flag governing use of out-of-plane distance */
/*     use_improp  logical flag governing use of improper dihedral term */
/*     use_imptor  logical flag governing use of improper torsion term */
/*     use_tors    logical flag governing use of torsional potential */
/*     use_pitors  logical flag governing use of pi-orbital torsion term */
/*     use_strtor  logical flag governing use of stretch-torsion term */
/*     use_tortor  logical flag governing use of torsion-torsion term */
/*     use_vdw     logical flag governing use of vdw der Waals potential */
/*     use_charge  logical flag governing use of charge-charge potential */
/*     use_chgdpl  logical flag governing use of charge-dipole potential */
/*     use_dipole  logical flag governing use of dipole-dipole potential */
/*     use_mpole   logical flag governing use of multipole potential */
/*     use_polar   logical flag governing use of polarization term */
/*     use_rxnfld  logical flag governing use of reaction field term */
/*     use_solv    logical flag governing use of continuum solvation */
/*     use_metal   logical flag governing use of ligand field term */
/*     use_geom    logical flag governing use of geometric restraints */
/*     use_extra   logical flag governing use of extra potential term */
/*     use_born    logical flag governing use of Born radii values */
/*     use_orbit   logical flag governing use of pisystem computation */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  vdwpot.i  --  specifics of van der Waals functional form  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     abuck      value of "A" constant in Buckingham vdw potential */
/*     bbuck      value of "B" constant in Buckingham vdw potential */
/*     cbuck      value of "C" constant in Buckingham vdw potential */
/*     ghal       value of "gamma" in buffered 14-7 vdw potential */
/*     dhal       value of "delta" in buffered 14-7 vdw potential */
/*     v2scale    factor by which 1-2 vdw interactions are scaled */
/*     v3scale    factor by which 1-3 vdw interactions are scaled */
/*     v4scale    factor by which 1-4 vdw interactions are scaled */
/*     v5scale    factor by which 1-5 vdw interactions are scaled */
/*     igauss     coefficients of Gaussian fit to vdw potential */
/*     ngauss     number of Gaussians used in fit to vdw potential */
/*     vdwindex   indexing mode (atom type or class) for vdw parameters */
/*     vdwtyp     type of van der Waals potential energy function */
/*     radtyp     type of parameter (sigma or R-min) for atomic size */
/*     radsiz     atomic size provided as radius or diameter */
/*     radrule    combining rule for atomic size parameters */
/*     epsrule    combining rule for vdw well depth parameters */
/*     gausstyp   type of Gaussian fit to van der Waals potential */




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
    warp_1.use_gda__ = TRUE_;
    mechanic_();

/*     store the initial values of the squared mean Gaussian width */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m2init[i__ - 1] = warp_1.m2[0];
    }

/*     get the number of optimized structures to be constructed */

    ntrial = 0;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___6);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ntrial, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    }
L10:
    if (ntrial <= 0) {
	io___7.ciunit = iounit_1.iout;
	s_wsfe(&io___7);
	e_wsfe();
	io___8.ciunit = iounit_1.input;
	s_rsfe(&io___8);
	do_fio(&c__1, (char *)&ntrial, (ftnlen)sizeof(integer));
	e_rsfe();
    }
    if (ntrial <= 0) {
	ntrial = 1;
    }

/*     see if random coordinates are desired as starting structures */

    randomize = TRUE_;
    if (ntrial == 1) {
	randomize = FALSE_;
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
	if (*(unsigned char *)answer == 'Y') {
	    randomize = TRUE_;
	}
    }
    if (randomize) {
	d__1 = (doublereal) atoms_1.n;
	boxsize = pow_dd(&d__1, &c_b15) * 10.;
    }

/*     get initial and final values of inverse temperature */

    bstart = -1.;
    bstop = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___18);
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&bstart, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L60;
	}
    }
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___19);
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&bstop, (ftnlen)sizeof(doublereal)
		);
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L60;
	}
    }
L60:
    if (bstart <= 0. || bstop <= 0.) {
	io___20.ciunit = iounit_1.iout;
	s_wsfe(&io___20);
	e_wsfe();
	io___21.ciunit = iounit_1.input;
	s_rsfe(&io___21);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___22);
	if (i__1 != 0) {
	    goto L90;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&bstart, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L90;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&bstop, (ftnlen)sizeof(doublereal)
		);
	if (i__1 != 0) {
	    goto L90;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L90;
	}
L90:
	;
    }
    if (bstart <= 0.) {
	bstart = .01;
    }
    if (bstop <= 0.) {
	bstop = 1e10;
    }

/*     write out a copy of coordinates for later update */

    i__1 = ntrial;
    for (itrial = 1; itrial <= i__1; ++itrial) {
	lext = 3;
	numeral_(&itrial, ext, &lext, (ftnlen)7);
/* Writing concatenation */
	i__2[0] = files_1.leng, a__1[0] = files_1.filename;
	i__2[1] = 1, a__1[1] = ".";
	i__2[2] = lext, a__1[2] = ext;
	s_cat(gdafile, a__1, i__2, &c__3, (ftnlen)120);
	version_(gdafile, "new", (ftnlen)120, (ftnlen)3);
	igda = freeunit_();
	o__1.oerr = 0;
	o__1.ounit = igda;
	o__1.ofnmlen = 120;
	o__1.ofnm = gdafile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	prtxyz_(&igda);
	cl__1.cerr = 0;
	cl__1.cunit = igda;
	cl__1.csta = 0;
	f_clos(&cl__1);
	s_copy(files_1.outfile, gdafile, (ftnlen)120, (ftnlen)120);

/*     set an initial box size and generate random coordinates */

	if (randomize) {
	    i__3 = atoms_1.n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		atoms_1.x[i__ - 1] = boxsize * random_();
		atoms_1.y[i__ - 1] = boxsize * random_();
		atoms_1.z__[i__ - 1] = boxsize * random_();
	    }
	}

/*     translate coordinates and M2's to optimization parameters */

	nvar = 0;
	i__3 = atoms_1.n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    ++nvar;
	    xx[nvar - 1] = atoms_1.x[i__ - 1];
	    ++nvar;
	    xx[nvar - 1] = atoms_1.y[i__ - 1];
	    ++nvar;
	    xx[nvar - 1] = atoms_1.z__[i__ - 1];
	}
	i__3 = atoms_1.n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    ++nvar;
	    xx[nvar - 1] = m2init[i__ - 1];
	}

/*     make changes to the potential to use potential smoothing */

	warp_1.use_smooth__ = TRUE_;
	potent_1.use_geom__ = TRUE_;
	s_copy(vdwpot_1.vdwtyp, "GAUSSIAN", (ftnlen)13, (ftnlen)8);

/*     make the call to the Bulirsch-Stoer integration routine */

	nstep = 0;
	s_copy(status, "       ", (ftnlen)7, (ftnlen)7);
	eps = 1e-8;
	h1 = .01;
	hmin = 0.;
	io___35.ciunit = iounit_1.iout;
	s_wsfe(&io___35);
	e_wsfe();
	gdastat_(&nstep, &bstart, xx, status, (ftnlen)7);
	diffeq_(&nvar, xx, &bstart, &bstop, &eps, &h1, &hmin, &nok, &nbad, (
		U_fp)gda1_);
	nstep = nok + nbad;

/*     make changes to the potential for standard optimization */

	warp_1.use_smooth__ = FALSE_;
	potent_1.use_geom__ = FALSE_;
	s_copy(vdwpot_1.vdwtyp, "LENNARD-JONES", (ftnlen)13, (ftnlen)13);

/*     make the call to the energy minimization routine */

	s_copy(mode, "DTNCG", (ftnlen)6, (ftnlen)5);
	s_copy(method, "AUTO", (ftnlen)6, (ftnlen)4);
	nvar = atoms_1.n * 3;
	grdmin = 1e-4;
	minima_1.nextiter = nstep + 1;
	tncg_(mode, method, &nvar, xx, &minimum, &grdmin, (D_fp)gda2_, (U_fp)
		gda3_, (U_fp)optsave_, (ftnlen)6, (ftnlen)6);
/*        call lbfgs (nvar,xx,minimum,grdmin,gda2,optsave) */
	io___42.ciunit = iounit_1.iout;
	s_wsfe(&io___42);
	do_fio(&c__1, (char *)&itrial, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	e_wsfe();

/*     translate optimization parameters into atomic coordinates */

	nvar = 0;
	i__3 = atoms_1.n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    ++nvar;
	    atoms_1.x[i__ - 1] = xx[nvar - 1];
	    ++nvar;
	    atoms_1.y[i__ - 1] = xx[nvar - 1];
	    ++nvar;
	    atoms_1.z__[i__ - 1] = xx[nvar - 1];
	}

/*     move the center of mass to the origin */

	xcm = 0.;
	ycm = 0.;
	zcm = 0.;
	i__3 = atoms_1.n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    xcm += atoms_1.x[i__ - 1];
	    ycm += atoms_1.y[i__ - 1];
	    zcm += atoms_1.z__[i__ - 1];
	}
	xcm /= (doublereal) atoms_1.n;
	ycm /= (doublereal) atoms_1.n;
	zcm /= (doublereal) atoms_1.n;
	i__3 = atoms_1.n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    atoms_1.x[i__ - 1] -= xcm;
	    atoms_1.y[i__ - 1] -= xcm;
	    atoms_1.z__[i__ - 1] -= xcm;
	}

/*     write the final coordinates into a file */

	igda = freeunit_();
	o__1.oerr = 0;
	o__1.ounit = igda;
	o__1.ofnmlen = 120;
	o__1.ofnm = gdafile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = igda;
	f_rew(&al__1);
	prtxyz_(&igda);
	cl__1.cerr = 0;
	cl__1.cunit = igda;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine gda1  --  gaussian density annealing gradient  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/* Subroutine */ int gda1_(doublereal *beta, doublereal *xx, doublereal *g)
{
    /* Format strings */
    static char fmt_10[] = "(\002 GDA1  --  Warning, Negative M2 at Atom\002"
	    ",i6,\002 with Value\002,d12.4)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static doublereal e, h__[1000000];
    static integer i__;
    static doublereal sum;
    static integer nvar;
    static doublereal hdiag[75000];
    static integer hinit[75000], hstop[75000], hindex[1000000];
    static doublereal derivs[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int hessian_(doublereal *, integer *, integer *, 
	    integer *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___48 = { 0, 0, 0, fmt_10, 0 };



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




/*     translate optimization parameters to coordinates and M2's */

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
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	warp_1.m2[i__ - 1] = xx[nvar];
	if (warp_1.m2[i__ - 1] < 0.) {
	    io___48.ciunit = iounit_1.iout;
	    s_wsfe(&io___48);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&warp_1.m2[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    warp_1.m2[i__ - 1] = -warp_1.m2[i__ - 1];
	}
    }

/*     compute and store the Cartesian energy gradient vector */

    gradient_(&e, derivs);

/*     translate the energy gradient into a dr/dbeta vector */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	g[nvar] = -(warp_1.m2[i__ - 1] / 3.) * derivs_ref(1, i__);
	++nvar;
	g[nvar] = -(warp_1.m2[i__ - 1] / 3.) * derivs_ref(2, i__);
	++nvar;
	g[nvar] = -(warp_1.m2[i__ - 1] / 3.) * derivs_ref(3, i__);
    }

/*     compute and store the Hessian elements */

    hessian_(h__, hinit, hstop, hindex, hdiag);

/*     translate the Hessian diagonal into a dM2/dbeta vector */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	sum = hdiag[i__ * 3 - 3] + hdiag[i__ * 3 - 2] + hdiag[i__ * 3 - 1];
/* Computing 2nd power */
	d__1 = warp_1.m2[i__ - 1] / 3.;
	g[nvar] = -(d__1 * d__1) * sum;
    }
    return 0;
} /* gda1_ */

#undef derivs_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function gda2  --  energy/gradient for TNCG optimization  ## */
/*     ##                                                            ## */
/*     ################################################################ */


doublereal gda2_(doublereal *xx, doublereal *g)
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
} /* gda2_ */

#undef derivs_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine gda3  --  Hessian values for TNCG optimization  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/* Subroutine */ int gda3_(char *mode, doublereal *xx, doublereal *h__, 
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
} /* gda3_ */

/* Main program alias */ int gda_ () { MAIN__ (); return 0; }
