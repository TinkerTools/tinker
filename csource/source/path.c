/* path.f -- translated by f2c (version 20050501).
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
    doublereal wfit[25000];
    integer nfit, ifit[50000]	/* was [2][25000] */;
} align_;

#define align_1 align_

struct {
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

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
    doublereal stpmin, stpmax, cappa, slpmax, angmax;
    integer intmax;
} linmin_;

#define linmin_1 linmin_

struct {
    doublereal fctmin, hguess;
    integer maxiter, nextiter;
} minima_;

#define minima_1 minima_

struct {
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

struct {
    doublereal p0[75000], p1[75000], pmid[75000], pvect[75000], pstep[75000], 
	    pzet[75000], pnorm, acoeff[49]	/* was [7][7] */, gc[525000]	
	    /* was [75000][7] */;
} paths_;

#define paths_1 paths_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__75000 = 75000;
static integer c__7 = 7;



/*     ############################################################### */
/*     ##  COPYRIGHT (C) 1991 by Shawn Huston & Jay William Ponder  ## */
/*     ##                    All Rights Reserved                    ## */
/*     ############################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  program path  --  conformational interconversion pathway  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "path" locates a series of structures equally spaced along */
/*     a conformational pathway connecting the input reactant and */
/*     product structures; a series of constrained optimizations */
/*     orthogonal to the path is done via Lagrangian multipliers */

/*     literature reference: */

/*     R. Czerminski and R. Elber, "Reaction Path Study of */
/*     Conformational Transitions in Flexible Systems: Applications */
/*     to Peptides", Journal of Chemical Physics, 92, 5580-5601 (1990) */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Enter Number of Path Points to Generate "
	    "[9] :  \002,$)";
    static char fmt_30[] = "(i10)";
    static char fmt_50[] = "(/,\002 Enter RMS Gradient per Atom Criterion"
	    "\002,\002 [0.1] :  \002,$)";
    static char fmt_60[] = "(f20.0)";
    static char fmt_70[] = "(/,\002 RMS Fit for Reactant and Product :\002,f"
	    "12.4)";
    static char fmt_80[] = "(/,\002 Reactant Potential Energy :\002,f12.4,/"
	    ",\002 Product Potential Energy : \002,f12.4)";
    static char fmt_90[] = "(/,\002 Path Point :\002,i12)";
    static char fmt_100[] = "(\002 Initial Point :\002,12x,f12.4)";
    static char fmt_110[] = "(\002 Optimized Point :\002,10x,f12.4)";
    static char fmt_120[] = "(\002 Target-Energy Difference :\002,d13.3)";
    static char fmt_130[] = "(\002 Gradient along Path :\002,6x,f12.4)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    static doublereal rmsvalue, g[75000];
    static integer i__, j, k;
    static doublereal p[75000], ge[75000];
    static integer ix, iy, iz;
    static doublereal sum, temp[525000]	/* was [75000][7] */;
    static integer nvar;
    static doublereal epot, etot, xtmp[25000], ytmp[25000], ztmp[25000];
    extern doublereal path1_();
    static doublereal epot0, epot1;
    extern /* Subroutine */ int final_(void), lbfgs_(integer *, doublereal *, 
	    doublereal *, doublereal *, D_fp, S_fp);
    static integer ipath, npath;
    static doublereal rplen;
    static logical exist;
    static doublereal grdmin;
    extern doublereal potnrg_(doublereal *, doublereal *);
    static char string[120];
    extern /* Subroutine */ int impose_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *), orthog_(integer *, integer *, integer *, 
	    doublereal *), invert_(integer *, integer *, doublereal *), 
	    getxyz_(void), initial_(void);
    static doublereal project;
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), optsave_(
	    integer *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static icilist io___9 = { 1, string, 1, 0, 120, 1 };
    static cilist io___10 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_30, 0 };
    static icilist io___13 = { 1, string, 1, 0, 120, 1 };
    static cilist io___14 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_130, 0 };



#define gc_ref(a_1,a_2) paths_1.gc[(a_2)*75000 + a_1 - 75001]
#define ifit_ref(a_1,a_2) align_1.ifit[(a_2)*2 + a_1 - 3]
#define temp_ref(a_1,a_2) temp[(a_2)*75000 + a_1 - 75001]
#define acoeff_ref(a_1,a_2) paths_1.acoeff[(a_2)*7 + a_1 - 8]



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  align.i  --  information for superposition of structures  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     wfit    weights assigned to atom pairs during superposition */
/*     nfit    number of atoms to use in superimposing two structures */
/*     ifit    atom numbers of pairs of atoms to be superimposed */




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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  linmin.i  --  parameters for line search minimization  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     stpmin   minimum step length in current line search direction */
/*     stpmax   maximum step length in current line search direction */
/*     cappa    stringency of line search (0=tight < cappa < 1=loose) */
/*     slpmax   projected gradient above which stepsize is reduced */
/*     angmax   maximum angle between search direction and -gradient */
/*     intmax   maximum number of interpolations during line search */




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
/*     ##  output.i  --  control of coordinate output file format  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     archive    logical flag to save structures in an archive */
/*     noversion  logical flag governing use of filename versions */
/*     overwrite  logical flag to overwrite intermediate files inplace */
/*     cyclesave  logical flag to mark use of numbered cycle files */
/*     coordtype  selects Cartesian, internal, rigid body or none */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  paths.i  --  parameters for Elber reaction path method  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     p0       reactant Cartesian coordinates as variables */
/*     p1       product Cartesian coordinates as variables */
/*     pmid     midpoint between the reactant and product */
/*     pvect    vector connecting the reactant and product */
/*     pstep    step per cycle along reactant-product vector */
/*     pzet     current projection on reactant-product vector */
/*     pnorm    length of the reactant-product vector */
/*     acoeff   transformation matrix 'A' from Elber paper */
/*     gc       gradients of the path constraints */




/*     initialize some constants and variables */

    initial_();

/*     get and store the initial structure coordinates */

    getxyz_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	paths_1.p0[i__ * 3 - 3] = atoms_1.x[i__ - 1];
	paths_1.p0[i__ * 3 - 2] = atoms_1.y[i__ - 1];
	paths_1.p0[i__ * 3 - 1] = atoms_1.z__[i__ - 1];
	xtmp[i__ - 1] = atoms_1.x[i__ - 1];
	ytmp[i__ - 1] = atoms_1.y[i__ - 1];
	ztmp[i__ - 1] = atoms_1.z__[i__ - 1];
    }

/*     get the coordinates for the final structure */

    getxyz_();
    mechanic_();

/*     set default values for some control variables */

    nvar = atoms_1.n * 3;
    output_1.cyclesave = TRUE_;
    linmin_1.stpmax = 1.;
    inform_1.iwrite = 0;
    if (inform_1.verbose) {
	inform_1.iprint = 1;
    } else {
	inform_1.iprint = 0;
    }

/*     get the number of path points to be generated */

    npath = -1;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___9);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&npath, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    }
L10:
    if (npath <= 0) {
	io___10.ciunit = iounit_1.iout;
	s_wsfe(&io___10);
	e_wsfe();
	io___11.ciunit = iounit_1.input;
	s_rsfe(&io___11);
	do_fio(&c__1, (char *)&npath, (ftnlen)sizeof(integer));
	e_rsfe();
    }
    if (npath <= 0) {
	npath = 9;
    }

/*     get the termination criterion as RMS gradient along path */

    grdmin = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___13);
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&grdmin, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L40;
	}
    }
L40:
    if (grdmin <= 0.) {
	io___14.ciunit = iounit_1.iout;
	s_wsfe(&io___14);
	e_wsfe();
	io___15.ciunit = iounit_1.input;
	s_rsfe(&io___15);
	do_fio(&c__1, (char *)&grdmin, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (grdmin <= 0.) {
	grdmin = .1;
    }

/*     superimpose the reactant and product structures */

    align_1.nfit = atoms_1.n;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ifit_ref(1, i__) = i__;
	ifit_ref(2, i__) = i__;
	align_1.wfit[i__ - 1] = atmtyp_1.mass[i__ - 1];
    }
    impose_(&atoms_1.n, xtmp, ytmp, ztmp, &atoms_1.n, atoms_1.x, atoms_1.y, 
	    atoms_1.z__, &rmsvalue);
    io___17.ciunit = iounit_1.iout;
    s_wsfe(&io___17);
    do_fio(&c__1, (char *)&rmsvalue, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     store the coordinates for the superimposed product */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	paths_1.p1[i__ * 3 - 3] = atoms_1.x[i__ - 1];
	paths_1.p1[i__ * 3 - 2] = atoms_1.y[i__ - 1];
	paths_1.p1[i__ * 3 - 1] = atoms_1.z__[i__ - 1];
    }

/*     write out the starting potential energy values */

    epot0 = potnrg_(paths_1.p0, g);
    epot1 = potnrg_(paths_1.p1, g);
    io___21.ciunit = iounit_1.iout;
    s_wsfe(&io___21);
    do_fio(&c__1, (char *)&epot0, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&epot1, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     construct step vector for getting */
/*     optimization-initial coordinates */

    rplen = (doublereal) (npath + 1);
    paths_1.pnorm = 0.;
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	paths_1.pvect[i__ - 1] = paths_1.p1[i__ - 1] - paths_1.p0[i__ - 1];
	paths_1.pstep[i__ - 1] = paths_1.pvect[i__ - 1] / rplen;
/* Computing 2nd power */
	d__1 = paths_1.pvect[i__ - 1];
	paths_1.pnorm += d__1 * d__1;
    }
    paths_1.pnorm = sqrt(paths_1.pnorm);

/*     set the gradient of constraints array */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ix = (i__ - 1) * 3 + 1;
	iy = ix + 1;
	iz = iy + 1;
	gc_ref(ix, 1) = paths_1.pvect[ix - 1];
	gc_ref(iy, 1) = paths_1.pvect[iy - 1];
	gc_ref(iz, 1) = paths_1.pvect[iz - 1];
	gc_ref(ix, 2) = atmtyp_1.mass[i__ - 1];
	gc_ref(iy, 2) = 0.;
	gc_ref(iz, 2) = 0.;
	gc_ref(ix, 3) = 0.;
	gc_ref(iy, 3) = atmtyp_1.mass[i__ - 1];
	gc_ref(iz, 3) = 0.;
	gc_ref(ix, 4) = 0.;
	gc_ref(iy, 4) = 0.;
	gc_ref(iz, 4) = atmtyp_1.mass[i__ - 1];
	gc_ref(ix, 5) = 0.;
	gc_ref(iy, 5) = atmtyp_1.mass[i__ - 1] * paths_1.p0[iz - 1];
	gc_ref(iz, 5) = -atmtyp_1.mass[i__ - 1] * paths_1.p0[iy - 1];
	gc_ref(ix, 6) = -atmtyp_1.mass[i__ - 1] * paths_1.p0[iz - 1];
	gc_ref(iy, 6) = 0.;
	gc_ref(iz, 6) = atmtyp_1.mass[i__ - 1] * paths_1.p0[ix - 1];
	gc_ref(ix, 7) = atmtyp_1.mass[i__ - 1] * paths_1.p0[iy - 1];
	gc_ref(iy, 7) = -atmtyp_1.mass[i__ - 1] * paths_1.p0[ix - 1];
	gc_ref(iz, 7) = 0.;
    }

/*     copy to temporary storage and orthogonalize */

    for (i__ = 1; i__ <= 7; ++i__) {
	i__1 = nvar;
	for (j = 1; j <= i__1; ++j) {
	    temp_ref(j, i__) = gc_ref(j, i__);
	}
    }
    orthog_(&nvar, &c__75000, &c__7, paths_1.gc);

/*     set the A matrix to transform sigma into C space */

    for (i__ = 1; i__ <= 7; ++i__) {
	for (k = 1; k <= 7; ++k) {
	    sum = 0.;
	    i__1 = nvar;
	    for (j = 1; j <= i__1; ++j) {
		sum += temp_ref(j, i__) * gc_ref(j, k);
	    }
	    acoeff_ref(i__, k) = sum;
	}
    }

/*     perform the matrix inversion to get A matrix */
/*     which transforms C into sigma space */

    invert_(&c__7, &c__7, paths_1.acoeff);

/*     set the current path point to be the reactant */

    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__ - 1] = paths_1.p0[i__ - 1];
    }

/*     loop over structures along path to be optimized */

    i__1 = npath;
    for (ipath = 1; ipath <= i__1; ++ipath) {
	io___32.ciunit = iounit_1.iout;
	s_wsfe(&io___32);
	do_fio(&c__1, (char *)&ipath, (ftnlen)sizeof(integer));
	e_wsfe();

/*     get r(zeta), set initial path point and energy */

	i__2 = nvar;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    paths_1.pzet[i__ - 1] = paths_1.p0[i__ - 1] + ipath * 
		    paths_1.pstep[i__ - 1];
	    p[i__ - 1] += paths_1.pstep[i__ - 1];
	}
	epot = potnrg_(p, ge);
	io___35.ciunit = iounit_1.iout;
	s_wsfe(&io___35);
	do_fio(&c__1, (char *)&epot, (ftnlen)sizeof(doublereal));
	e_wsfe();

/*     call optimizer to get constrained minimum */

	lbfgs_(&nvar, p, &etot, &grdmin, (D_fp)path1_, (S_fp)optsave_);
/*        call ocvm (nvar,p,etot,grdmin,path1,optsave) */

/*     print energy and constraint value at the minimum */

	epot = potnrg_(p, ge);
	io___37.ciunit = iounit_1.iout;
	s_wsfe(&io___37);
	do_fio(&c__1, (char *)&epot, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___38.ciunit = iounit_1.iout;
	s_wsfe(&io___38);
	d__1 = etot - epot;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfe();

/*     write coordinates of the current path point */

	optsave_(&ipath, &epot, p);

/*     find projection of the gradient along path direction */

	project = 0.;
	i__2 = nvar;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    project += ge[i__ - 1] * paths_1.pvect[i__ - 1] / paths_1.pnorm;
	}
	io___40.ciunit = iounit_1.iout;
	s_wsfe(&io___40);
	do_fio(&c__1, (char *)&project, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef acoeff_ref
#undef temp_ref
#undef ifit_ref
#undef gc_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  function path1  --  value and gradient of target function  ## */
/*     ##                                                             ## */
/*     ################################################################# */


doublereal path1_(doublereal *p, doublereal *gt)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, j;
    static doublereal ge[75000];
    static integer ix, iy, iz;
    static doublereal xx, yy, zz;
    static integer nvar;
    static doublereal cnst[7], gamma[7], sigma[7], cterm;
    extern doublereal potnrg_(doublereal *, doublereal *);


#define gc_ref(a_1,a_2) paths_1.gc[(a_2)*75000 + a_1 - 75001]
#define acoeff_ref(a_1,a_2) paths_1.acoeff[(a_2)*7 + a_1 - 8]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  paths.i  --  parameters for Elber reaction path method  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     p0       reactant Cartesian coordinates as variables */
/*     p1       product Cartesian coordinates as variables */
/*     pmid     midpoint between the reactant and product */
/*     pvect    vector connecting the reactant and product */
/*     pstep    step per cycle along reactant-product vector */
/*     pzet     current projection on reactant-product vector */
/*     pnorm    length of the reactant-product vector */
/*     acoeff   transformation matrix 'A' from Elber paper */
/*     gc       gradients of the path constraints */




/*     get the value of the potential energy */

    /* Parameter adjustments */
    --gt;
    --p;

    /* Function Body */
    nvar = atoms_1.n * 3;
    ret_val = potnrg_(&p[1], ge);

/*     construct the Lagrangian multipliers */

    for (i__ = 1; i__ <= 7; ++i__) {
	gamma[i__ - 1] = 0.;
	i__1 = nvar;
	for (j = 1; j <= i__1; ++j) {
	    gamma[i__ - 1] -= ge[j - 1] * gc_ref(j, i__);
	}
    }

/*     set the path value, translation and rotation constraints */

    for (i__ = 1; i__ <= 7; ++i__) {
	cnst[i__ - 1] = 0.;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ix = (i__ - 1) * 3 + 1;
	iy = ix + 1;
	iz = iy + 1;
	xx = p[ix] - paths_1.pzet[ix - 1];
	yy = p[iy] - paths_1.pzet[iy - 1];
	zz = p[iz] - paths_1.pzet[iz - 1];
	cnst[0] = cnst[0] + xx * paths_1.pvect[ix - 1] + yy * paths_1.pvect[
		iy - 1] + zz * paths_1.pvect[iz - 1];
	cnst[1] += atmtyp_1.mass[i__ - 1] * (p[ix] - paths_1.p0[ix - 1]);
	cnst[2] += atmtyp_1.mass[i__ - 1] * (p[iy] - paths_1.p0[iy - 1]);
	cnst[3] += atmtyp_1.mass[i__ - 1] * (p[iz] - paths_1.p0[iz - 1]);
	cnst[4] += atmtyp_1.mass[i__ - 1] * (p[iy] * paths_1.p0[iz - 1] - p[
		iz] * paths_1.p0[iy - 1]);
	cnst[5] += atmtyp_1.mass[i__ - 1] * (p[iz] * paths_1.p0[ix - 1] - p[
		ix] * paths_1.p0[iz - 1]);
	cnst[6] += atmtyp_1.mass[i__ - 1] * (p[ix] * paths_1.p0[iy - 1] - p[
		iy] * paths_1.p0[ix - 1]);
    }

/*     construct the orthonormal "sigma" constraints */

    for (i__ = 1; i__ <= 7; ++i__) {
	sigma[i__ - 1] = 0.;
	for (j = 1; j <= 7; ++j) {
	    sigma[i__ - 1] += acoeff_ref(i__, j) * cnst[j - 1];
	}
    }

/*     find the target function value */

    cterm = 0.;
    for (i__ = 1; i__ <= 7; ++i__) {
	cterm += gamma[i__ - 1] * sigma[i__ - 1];
    }
    ret_val += cterm;

/*     construct the gradient of the target function */

    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gt[i__] = ge[i__ - 1];
	for (j = 1; j <= 7; ++j) {
	    gt[i__] += gamma[j - 1] * gc_ref(i__, j);
	}
    }
    return ret_val;
} /* path1_ */

#undef acoeff_ref
#undef gc_ref




/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  function potnrg  --  potential energy and gradient  ## */
/*     ##                                                      ## */
/*     ########################################################## */


doublereal potnrg_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static integer i__;
    static doublereal deriv[75000]	/* was [3][25000] */, energy;


#define deriv_ref(a_1,a_2) deriv[(a_2)*3 + a_1 - 4]



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
/*     ##  paths.i  --  parameters for Elber reaction path method  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     p0       reactant Cartesian coordinates as variables */
/*     p1       product Cartesian coordinates as variables */
/*     pmid     midpoint between the reactant and product */
/*     pvect    vector connecting the reactant and product */
/*     pstep    step per cycle along reactant-product vector */
/*     pzet     current projection on reactant-product vector */
/*     pnorm    length of the reactant-product vector */
/*     acoeff   transformation matrix 'A' from Elber paper */
/*     gc       gradients of the path constraints */




/*     copy position vector into atomic coordinates */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xx[i__ * 3 - 2];
	atoms_1.y[i__ - 1] = xx[i__ * 3 - 1];
	atoms_1.z__[i__ - 1] = xx[i__ * 3];
    }

/*     compute potential energy and Cartesian derivatives */

    gradient_(&energy, deriv);

/*     set the energy value and gradient vector */

    ret_val = energy;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__ * 3 - 2] = deriv_ref(1, i__);
	g[i__ * 3 - 1] = deriv_ref(2, i__);
	g[i__ * 3] = deriv_ref(3, i__);
    }
    return ret_val;
} /* potnrg_ */

#undef deriv_ref


/* Main program alias */ int path_ () { MAIN__ (); return 0; }
