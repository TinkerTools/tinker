/* sniffer.f -- translated by f2c (version 20050501).
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
    doublereal scale[75000];
    logical set_scale__;
} scales_;

#define scales_1 scales_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  program sniffer  --  discrete generalized descent search  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "sniffer" performs a global energy minimization using a */
/*     discrete version of Griewank's global search trajectory */

/*     literature references: */

/*     A. O. Griewank, "Generalized Descent for Global Optimization", */
/*     Journal of Optimization Theory and Applications, 34, 11-39 (1981) */

/*     R. A. R. Butler and E. E. Slaminka, "An Evaluation of the Sniffer */
/*     Global Optimization Algorithm Using Standard Test Functions", */
/*     Journal of Computational Physics, 99, 28-32 (1992) */

/*     J. W. Rogers, Jr. and R. A. Donnelly, "Potential Transformation */
/*     Methods for Large-Scale Global Optimization", SIAM Journal of */
/*     Optimization, 5, 871-891 (1995) */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Enter Number of Steps in the Initial Se"
	    "t\002,\002 [100] :  \002,$)";
    static char fmt_30[] = "(f20.0)";
    static char fmt_50[] = "(/,\002 Enter Target Energy for the Global Minim"
	    "um\002,\002 [0.0] :  \002,$)";
    static char fmt_60[] = "(f20.0)";
    static char fmt_80[] = "(/,\002 Enter RMS Gradient per Atom Criterion [1"
	    ".0] :  \002,$)";
    static char fmt_90[] = "(f20.0)";
    static char fmt_100[] = "(/,\002 Discrete Generalized Descent Globa"
	    "l\002,\002 Optimization :\002)";
    static char fmt_110[] = "(/,4x,\002Iter\002,11x,\002F Value\002,13x,\002"
	    "G RMS\002,8x,\002X Move\002,9x,\002Angle\002,/)";
    static char fmt_120[] = "(i8,2f18.4,2f14.4)";
    static char fmt_130[] = "(i8,2d18.4,2f14.4)";
    static char fmt_140[] = "(/,\002 SNIFFER  --  Incomplete Convergence due"
	    " to \002,a9)";
    static char fmt_150[] = "(/,\002 SNIFFER  --  Normal Termination due to"
	    " \002,a9)";
    static char fmt_160[] = "(/,\002 Final Function Value :\002,f15.4,/,\002"
	    " Final RMS Gradient :  \002,f15.4,/,\002 Final Gradient Norm :"
	    " \002,f15.4)";
    static char fmt_170[] = "(/,\002 Final Function Value :\002,f15.4,/,\002"
	    " Final RMS Gradient :  \002,d15.4,/,\002 Final Gradient Norm :"
	    " \002,d15.4)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2], i__3;
    doublereal d__1, d__2;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);
    integer f_rew(alist *);
    double acos(doublereal);
    integer i_dnnt(doublereal *), s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void), gradient_(doublereal *, 
	    doublereal *);
    extern integer freeunit_(void);
    static doublereal d__[75000], f, g[75000];
    static integer i__, j;
    static doublereal mu, xx[75000], dot, eps, rms;
    static logical done;
    static doublereal fmin, gmin;
    static integer imin, nvar;
    static doublereal grms, size;
    static integer stop;
    static doublereal alpha, angle;
    extern /* Subroutine */ int final_(void);
    static doublereal mufac, dnorm;
    static integer niter;
    static doublereal gnorm;
    static integer istep;
    static logical exist;
    static integer start;
    static doublereal epsfac, scaler;
    extern /* Subroutine */ int getref_(integer *);
    static doublereal grdmin, cosine, derivs[75000]	/* was [3][25000] */;
    static char string[120], status[9];
    extern /* Subroutine */ int getxyz_(void), prtxyz_(integer *), makeref_(
	    integer *);
    static char minfile[120];
    static doublereal stepfac;
    extern /* Subroutine */ int initial_(void), nextarg_(char *, logical *, 
	    ftnlen);
    static doublereal minimum;
    static integer maxstep;
    extern /* Subroutine */ int optsave_(integer *, doublereal *, doublereal *
	    ), version_(char *, char *, ftnlen, ftnlen);
    extern doublereal sniffer1_(doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static icilist io___4 = { 1, string, 1, 0, 120, 1 };
    static cilist io___5 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_30, 0 };
    static icilist io___7 = { 1, string, 1, 0, 120, 1 };
    static cilist io___8 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_60, 0 };
    static icilist io___11 = { 1, string, 1, 0, 120, 1 };
    static cilist io___12 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___48 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___49 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___54 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___55 = { 0, 0, 0, fmt_170, 0 };



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
/*     ##  scales.i  --  parameter scale factors for optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     scale      multiplicative factor for each optimization parameter */
/*     set_scale  logical flag to show if scale factors have been set */




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
    makeref_(&c__1);

/*     get the number of steps in the initial block */

    maxstep = -1;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___4);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&maxstep, (ftnlen)sizeof(integer))
		;
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    }
L10:
    if (maxstep <= 0) {
	io___5.ciunit = iounit_1.iout;
	s_wsfe(&io___5);
	e_wsfe();
	io___6.ciunit = iounit_1.input;
	s_rsfe(&io___6);
	do_fio(&c__1, (char *)&maxstep, (ftnlen)sizeof(integer));
	e_rsfe();
    }
    if (maxstep <= 0) {
	maxstep = 100;
    }

/*     get the target value for the global energy minimum */

    minima_1.fctmin = 1e6;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___7);
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&minima_1.fctmin, (ftnlen)sizeof(
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
    if (minima_1.fctmin >= 1e6) {
	io___8.ciunit = iounit_1.iout;
	s_wsfe(&io___8);
	e_wsfe();
	io___9.ciunit = iounit_1.input;
	s_rsfe(&io___9);
	do_fio(&c__1, (char *)&minima_1.fctmin, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (minima_1.fctmin >= 1e6) {
	minima_1.fctmin = 0.;
    }

/*     get termination criterion as RMS gradient per atom */

    grdmin = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___11);
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&grdmin, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L70;
	}
    }
L70:
    if (grdmin <= 0.) {
	io___12.ciunit = iounit_1.iout;
	s_wsfe(&io___12);
	e_wsfe();
	io___13.ciunit = iounit_1.input;
	s_rsfe(&io___13);
	do_fio(&c__1, (char *)&grdmin, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (grdmin <= 0.) {
	grdmin = 1.;
    }

/*     write out a copy of coordinates for later update */

    imin = freeunit_();
/* Writing concatenation */
    i__2[0] = files_1.leng, a__1[0] = files_1.filename;
    i__2[1] = 4, a__1[1] = ".xyz";
    s_cat(minfile, a__1, i__2, &c__2, (ftnlen)120);
    version_(minfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = imin;
    o__1.ofnmlen = 120;
    o__1.ofnm = minfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    prtxyz_(&imin);
    cl__1.cerr = 0;
    cl__1.cunit = imin;
    cl__1.csta = 0;
    f_clos(&cl__1);
    s_copy(files_1.outfile, minfile, (ftnlen)120, (ftnlen)120);
    s_copy(output_1.coordtype, "CARTESIAN", (ftnlen)9, (ftnlen)9);

/*     set scaling parameter for function and derivative values; */
/*     use square root of median eigenvalue of a typical Hessian */

    scales_1.set_scale__ = TRUE_;
    scaler = 1.;
    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    for (j = 1; j <= 3; ++j) {
		++nvar;
		scales_1.scale[nvar - 1] = scaler;
	    }
	}
    }

/*     scale the coordinates of each active atom */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    xx[nvar - 1] = atoms_1.x[i__ - 1] * scales_1.scale[nvar - 1];
	    ++nvar;
	    xx[nvar - 1] = atoms_1.y[i__ - 1] * scales_1.scale[nvar - 1];
	    ++nvar;
	    xx[nvar - 1] = atoms_1.z__[i__ - 1] * scales_1.scale[nvar - 1];
	}
    }

/*     set initial values for the control parameters */

    epsfac = 1.1;
    mufac = 1.7;
    stepfac = 1.1;

/*     set initial values for optimization parameters */

    inform_1.iprint = 1;
    inform_1.iwrite = 100;
    rms = sqrt((doublereal) atoms_1.n);
    start = 0;
    stop = start + maxstep;
    eps = 1.;
    mu = 1.;
    linmin_1.stpmax = rms * .1;
    linmin_1.stpmin = .001;

/*     initialize unit direction vector along negative gradient */

    f = sniffer1_(xx, g);
    gnorm = 0.;
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = g[i__ - 1];
	gnorm += d__1 * d__1;
    }
    gnorm = sqrt(gnorm);
    grms = gnorm / rms;
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__ - 1] = -g[i__ - 1] / gnorm;
    }
    fmin = f;
    gmin = grms;

/*     tests of the successful termination criteria */

    if (fmin <= minima_1.fctmin) {
	s_copy(status, "TargetVal", (ftnlen)9, (ftnlen)9);
	done = TRUE_;
    } else if (gmin <= grdmin) {
	s_copy(status, "SmallGrad", (ftnlen)9, (ftnlen)9);
	done = TRUE_;
    } else {
	done = FALSE_;
    }

/*     print header information prior to iterations */

    if (inform_1.iprint > 0) {
	io___38.ciunit = iounit_1.iout;
	s_wsfe(&io___38);
	e_wsfe();
    }

/*     perform a set of basic sniffer search steps */

    niter = 0;
    while(! done) {
	io___40.ciunit = iounit_1.iout;
	s_wsfe(&io___40);
	e_wsfe();
	i__1 = stop;
	for (istep = start; istep <= i__1; ++istep) {

/*     get the current energy and gradient values */

	    f = sniffer1_(xx, g);

/*     if current energy is lowest yet, save the coordinates */

	    if (f < fmin) {
		nvar = 0;
		i__3 = atoms_1.n;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    if (usage_1.use[i__ - 1]) {
			++nvar;
			atoms_1.x[i__ - 1] = xx[nvar - 1] / scales_1.scale[
				nvar - 1];
			++nvar;
			atoms_1.y[i__ - 1] = xx[nvar - 1] / scales_1.scale[
				nvar - 1];
			++nvar;
			atoms_1.z__[i__ - 1] = xx[nvar - 1] / scales_1.scale[
				nvar - 1];
		    }
		}
		makeref_(&c__1);
		imin = freeunit_();
		o__1.oerr = 0;
		o__1.ounit = imin;
		o__1.ofnmlen = 120;
		o__1.ofnm = minfile;
		o__1.orl = 0;
		o__1.osta = "old";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
		al__1.aerr = 0;
		al__1.aunit = imin;
		f_rew(&al__1);
		prtxyz_(&imin);
		cl__1.cerr = 0;
		cl__1.cunit = imin;
		cl__1.csta = 0;
		f_clos(&cl__1);
	    }

/*     get rms gradient and dot product with search direction */

	    gnorm = 0.;
	    dot = 0.;
	    i__3 = nvar;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		gnorm += g[i__ - 1] * g[i__ - 1];
		dot += d__[i__ - 1] * g[i__ - 1];
	    }
	    gnorm = sqrt(gnorm);
	    grms = gnorm / (scaler * rms);

/*     compute the next direction vector and its length */

/* Computing MAX */
	    d__1 = 0., d__2 = (eps + 1.) * dot + 1.;
	    alpha = max(d__1,d__2);
	    dnorm = 0.;
	    i__3 = nvar;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		d__[i__ - 1] = -eps * g[i__ - 1] + alpha * d__[i__ - 1];
		dnorm += d__[i__ - 1] * d__[i__ - 1];
	    }
	    dnorm = sqrt(dnorm);

/*     normalize direction and get angle with negative gradient */

	    dot = 0.;
	    i__3 = nvar;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		d__[i__ - 1] /= dnorm;
		dot += d__[i__ - 1] * g[i__ - 1];
	    }
	    cosine = -dot / gnorm;
/* Computing MIN */
	    d__1 = 1., d__2 = max(-1.,cosine);
	    cosine = min(d__1,d__2);
	    angle = acos(cosine) * 57.29577951308232088;

/*     move atomic positions along the direction vector */

/* Computing MIN */
	    d__1 = linmin_1.stpmax, d__2 = mu * (f - minima_1.fctmin);
	    size = min(d__1,d__2);
	    i__3 = nvar;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		xx[i__ - 1] += size * d__[i__ - 1];
	    }

/*     compute the size of the step taken */

/* Computing MIN */
	    d__1 = linmin_1.stpmax, d__2 = mu * (f - minima_1.fctmin);
	    size = min(d__1,d__2);
	    size /= rms;

/*     update the best value and gradient found so far */

	    fmin = min(fmin,f);
	    gmin = min(gmin,grms);

/*     print intermediate results every few iterations */

	    if (inform_1.iprint > 0) {
		if (done || niter % inform_1.iprint == 0) {
		    if (f < 1e12 && f > -1e11 && grms < 1e12) {
			io___48.ciunit = iounit_1.iout;
			s_wsfe(&io___48);
			do_fio(&c__1, (char *)&istep, (ftnlen)sizeof(integer))
				;
			do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
			do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&size, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&angle, (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    } else {
			io___49.ciunit = iounit_1.iout;
			s_wsfe(&io___49);
			do_fio(&c__1, (char *)&istep, (ftnlen)sizeof(integer))
				;
			do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
			do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&size, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&angle, (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    }
		}
	    }
	}

/*     tests of the various termination and error criteria */

	if (fmin <= minima_1.fctmin) {
	    s_copy(status, "TargetVal", (ftnlen)9, (ftnlen)9);
	    done = TRUE_;
	} else if (gmin <= grdmin) {
	    s_copy(status, "SmallGrad", (ftnlen)9, (ftnlen)9);
	    done = TRUE_;
	} else if (size <= linmin_1.stpmin) {
	    s_copy(status, "SmallMove", (ftnlen)9, (ftnlen)9);
	    done = TRUE_;
	}

/*     write the final coordinates for this set of steps */

	++niter;
	if (output_1.cyclesave) {
	    optsave_(&niter, &fmin, xx);
	}

/*     update the optimization parameters for the next set */

	eps *= epsfac;
	mu /= mufac;
	d__1 = (doublereal) maxstep * stepfac;
	maxstep = i_dnnt(&d__1);
	start = stop + 1;
	stop = start + maxstep;
    }

/*     write message about satisfaction of termination criteria */

    if (s_cmp(status, "SmallMove", (ftnlen)9, (ftnlen)9) == 0) {
	io___50.ciunit = iounit_1.iout;
	s_wsfe(&io___50);
	do_fio(&c__1, status, (ftnlen)9);
	e_wsfe();
    } else {
	io___51.ciunit = iounit_1.iout;
	s_wsfe(&io___51);
	do_fio(&c__1, status, (ftnlen)9);
	e_wsfe();
    }

/*     use lowest energy structure as global minimum estimate */

    getref_(&c__1);

/*     write out final function value and gradient */

    gradient_(&minimum, derivs);
    gnorm = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    for (j = 1; j <= 3; ++j) {
/* Computing 2nd power */
		d__1 = derivs_ref(j, i__);
		gnorm += d__1 * d__1;
	    }
	}
    }
    gnorm = sqrt(gnorm);
    grms = gnorm / rms;
    if (grms > 1e-4) {
	io___54.ciunit = iounit_1.iout;
	s_wsfe(&io___54);
	do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else {
	io___55.ciunit = iounit_1.iout;
	s_wsfe(&io___55);
	do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     write the final coordinates into a file */

    imin = freeunit_();
    o__1.oerr = 0;
    o__1.ounit = imin;
    o__1.ofnmlen = 120;
    o__1.ofnm = minfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    al__1.aerr = 0;
    al__1.aunit = imin;
    f_rew(&al__1);
    prtxyz_(&imin);
    cl__1.cerr = 0;
    cl__1.cunit = imin;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef derivs_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function sniffer1  --  energy and gradient for sniffer  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "sniffer1" is a service routine that computes the energy */
/*     and gradient for the Sniffer global optimization method */


doublereal sniffer1_(doublereal *xx, doublereal *g)
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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  scales.i  --  parameter scale factors for optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     scale      multiplicative factor for each optimization parameter */
/*     set_scale  logical flag to show if scale factors have been set */




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




/*     translate optimization parameters to atomic coordinates */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    atoms_1.x[i__ - 1] = xx[nvar] / scales_1.scale[nvar - 1];
	    ++nvar;
	    atoms_1.y[i__ - 1] = xx[nvar] / scales_1.scale[nvar - 1];
	    ++nvar;
	    atoms_1.z__[i__ - 1] = xx[nvar] / scales_1.scale[nvar - 1];
	}
    }

/*     compute and store the energy and gradient */

    gradient_(&e, derivs);
    ret_val = e;

/*     store Cartesian gradient as optimization gradient */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    g[nvar] = derivs_ref(1, i__) / scales_1.scale[nvar - 1];
	    ++nvar;
	    g[nvar] = derivs_ref(2, i__) / scales_1.scale[nvar - 1];
	    ++nvar;
	    g[nvar] = derivs_ref(3, i__) / scales_1.scale[nvar - 1];
	}
    }
    return ret_val;
} /* sniffer1_ */

#undef derivs_ref


/* Main program alias */ int sniffer_ () { MAIN__ (); return 0; }
