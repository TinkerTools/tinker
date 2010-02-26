/* tncg.f -- translated by f2c (version 20050501).
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
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

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
    integer norbit, iorbit[100], reorbit, piperp[300]	/* was [3][100] */, 
	    nbpi, ibpi[600]	/* was [3][200] */, ntpi, itpi[800]	/* 
	    was [2][400] */;
    logical listpi[25000];
} piorbs_;

#define piorbs_1 piorbs_

struct {
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__3 = 3;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine tncg  --  truncated newton optimization method  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "tncg" implements a truncated Newton optimization algorithm */
/*     in which a preconditioned linear conjugate gradient method is */
/*     used to approximately solve Newton's equations; special features */
/*     include use of an explicit sparse Hessian or finite-difference */
/*     gradient-Hessian products within the PCG iteration; the exact */
/*     Newton search directions can be used optionally; by default the */
/*     algorithm checks for negative curvature to prevent convergence */
/*     to a stationary point having negative eigenvalues; if a saddle */
/*     point is desired this test can be removed by disabling "negtest" */

/*     literature references: */

/*     J. W. Ponder and F. M Richards, "An Efficient Newton-like */
/*     Method for Molecular Mechanics Energy Minimization of */
/*     Large Molecules", Journal of Computational Chemistry, */
/*     8, 1016-1024 (1987) */

/*     R. S. Dembo and T. Steihaug, "Truncated-Newton Algorithms */
/*     for Large-Scale Unconstrained Optimization", Mathematical */
/*     Programming, 26, 190-212 (1983) */

/*     variables and parameters: */

/*     mode       determines optimization method; choice of */
/*                  Newton's method, truncated Newton, or */
/*                  truncated Newton with finite differencing */
/*     method     determines which type of preconditioning will */
/*                  be used on the Newton equations; choice */
/*                  of none, diagonal, 3x3 block diagonal, */
/*                  SSOR or incomplete Cholesky preconditioning */
/*     nvar       number of parameters in the objective function */
/*     minimum    upon return contains the best value of the */
/*                  function found during the optimization */
/*     f          contains current best value of function */
/*     x          contains starting point upon input, upon */
/*                  return contains the best point found */
/*     g          contains gradient of current best point */
/*     h          contains the Hessian matrix values in an */
/*                  indexed linear array */
/*     h_mode     controls amount of Hessian matrix computed; */
/*                  either the full matrix, diagonal or none */
/*     h_init     points to the first Hessian matrix element */
/*                  associated with each parameter */
/*     h_stop     points to the last Hessian matrix element */
/*                  associated with each parameter */
/*     h_index    contains second parameter involved in each */
/*                  element of the Hessian array */
/*     h_diag     contains diagonal of the Hessian matrix */
/*     p          search direction resulting from pcg iteration */
/*     f_move     function decrease over last tncg iteration */
/*     f_old      function value at end of last iteration */
/*     x_move     rms movement per atom over last tn iteration */
/*     x_old      parameters value at end of last tn iteration */
/*     g_norm     Euclidian norm of the gradient vector */
/*     g_rms      root mean square gradient value */
/*     fg_call    cumulative number of function/gradient calls */
/*     grdmin     termination criterion based on RMS gradient */
/*     iprint     print iteration results every iprint iterations */
/*     iwrite     call user-supplied output every iwrite iterations */
/*     newhess    number of iterations between the computation */
/*                  of new Hessian matrix values */
/*     negtest    determines whether test for negative curvature */
/*                  is performed during the PCG iterations */
/*     maxiter    maximum number of tncg iterations to attempt */

/*     parameters used in the line search: */

/*     cappa      accuarcy of line search control  (0 < cappa < 1) */
/*     stpmin     minimum allowed line search step size */
/*     stpmax     maximum allowed line search step size */
/*     angmax     maximum angle between search and -grad directions */
/*     intmax     maximum number of interpolations in line search */

/*     required external routines: */

/*     fgvalue    function to evaluate function and gradient values */
/*     hmatrix    subroutine which evaluates Hessian diagonal */
/*                  and large off-diagonal matrix elements */
/*     optsave    subroutine to write out info about current status */


/* Subroutine */ int tncg_(char *mode, char *method, integer *nvar, 
	doublereal *x, doublereal *minimum, doublereal *grdmin, D_fp fgvalue, 
	S_fp hmatrix, S_fp optsave, ftnlen mode_len, ftnlen method_len)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 TNCG  --  Too many Parameters,\002,\002 "
	    "Increase the Value of MAXVAR\002)";
    static char fmt_30[] = "(/,\002 Full-Newton Conjugate-Gradient\002,\002 "
	    "Optimization :\002)";
    static char fmt_40[] = "(/,\002 Truncated-Newton Conjugate-Gradient\002"
	    ",\002 Optimization :\002)";
    static char fmt_50[] = "(/,\002 Finite-Difference Truncated-Newton\002"
	    ",\002 Conjugate-Gradient Optimization :\002)";
    static char fmt_60[] = "(/,\002 Variable-Mode Truncated-Newton\002,\002 "
	    "Conjugate-Gradient Optimization :\002)";
    static char fmt_70[] = "(/,\002 Algorithm : \002,a6,5x,\002Preconditioni"
	    "ng : \002,a6,5x,\002 RMS Grad : \002,d8.2)";
    static char fmt_80[] = "(/,\002 TN Iter   F Value       G RMS     F Move"
	    "  \002,\002  X Move   CG Iter   Solve   FG Call\002,/)";
    static char fmt_90[] = "(i6,f12.4,f12.4,41x,i7)";
    static char fmt_100[] = "(i6,d12.4,d12.4,41x,i7)";
    static char fmt_110[] = "(/,\002 TNCG  --  Normal Termination due to Sma"
	    "llGrad\002)";
    static char fmt_120[] = "(/,\002 TNCG  --  Normal Termination due to Sma"
	    "llFct\002)";
    static char fmt_130[] = "(/,\002 TNCG  --  Incomplete Convergence\002"
	    ",\002 due to IterLimit\002)";
    static char fmt_140[] = "(i6,f12.4,f12.4,f11.4,f10.4,i8,3x,a9,i7)";
    static char fmt_150[] = "(i6,d12.4,d12.4,d11.4,f10.4,i8,3x,a9,i7)";
    static char fmt_160[] = "(/,\002 TNCG  --  Normal Termination due to "
	    "\002,a9)";
    static char fmt_170[] = "(/,\002 TNCG  --  Incomplete Convergence\002"
	    ",\002 due to \002,a9)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    double sqrt(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static logical automode;
    static doublereal f, g[75000], h__[1000000];
    static integer i__;
    static doublereal p[75000];
    static logical automatic;
    static char info_solve__[9];
    static doublereal rms;
    static char info_search__[9];
    static logical done;
    static integer nerr, next;
    static doublereal f_old__, angle;
    extern /* Subroutine */ int piscf_(void);
    static doublereal x_old__[75000], g_rms__, h_diag__[75000];
    static char h_mode__[4];
    extern /* Subroutine */ int search_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , D_fp, char *, ftnlen);
    static integer h_init__[75000];
    static doublereal f_move__;
    static char record[120];
    static integer maxerr, h_stop__[75000];
    static doublereal x_move__, g_norm__;
    static char status[9], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static integer fg_call__, iter_cg__, h_index__[1000000], iter_tn__;
    static logical negtest;
    static integer newhess;
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), tnsolve_(char *, char *, logical *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, doublereal *, integer *, integer *, 
	    integer *, D_fp, char *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };
    static icilist io___15 = { 1, string, 1, 0, 120, 1 };
    static icilist io___16 = { 1, string, 1, 0, 120, 1 };
    static icilist io___17 = { 1, string, 1, 0, 120, 1 };
    static icilist io___18 = { 1, string, 1, 0, 120, 1 };
    static icilist io___19 = { 1, string, 1, 0, 120, 1 };
    static icilist io___20 = { 1, string, 1, 0, 120, 1 };
    static icilist io___21 = { 1, string, 1, 0, 120, 1 };
    static icilist io___22 = { 1, string, 1, 0, 120, 1 };
    static icilist io___23 = { 1, string, 1, 0, 120, 1 };
    static icilist io___24 = { 1, string, 1, 0, 120, 1 };
    static cilist io___28 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___46 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___58 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___59 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___60 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_170, 0 };



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  piorbs.i  --  conjugated system in the current structure  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     norbit    total number of pisystem orbitals in the system */
/*     iorbit    numbers of the atoms containing pisystem orbitals */
/*     reorbit   number of evaluations between orbital updates */
/*     piperp    atoms defining a normal plane to each orbital */
/*     nbpi      total number of bonds affected by the pisystem */
/*     ibpi      bond and piatom numbers for each pisystem bond */
/*     ntpi      total number of torsions affected by the pisystem */
/*     itpi      torsion and pibond numbers for each pisystem torsion */
/*     listpi    atom list indicating whether each atom has an orbital */




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




/*     check number of variables and get type of optimization */

    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*nvar > 75000) {
	io___1.ciunit = iounit_1.iout;
	s_wsfe(&io___1);
	e_wsfe();
	return 0;
    }
    rms = sqrt((doublereal) (*nvar));
    if (s_cmp(output_1.coordtype, "CARTESIAN", (ftnlen)9, (ftnlen)9) == 0) {
	rms /= sqrt(3.);
    } else if (s_cmp(output_1.coordtype, "RIGIDBODY", (ftnlen)9, (ftnlen)9) ==
	     0) {
	rms /= sqrt(6.);
    }

/*     set default parameters for the optimization */

    if (minima_1.fctmin == 0.) {
	minima_1.fctmin = -1e6;
    }
    if (inform_1.iwrite < 0) {
	inform_1.iwrite = 1;
    }
    if (inform_1.iprint < 0) {
	inform_1.iprint = 1;
    }
    if (minima_1.maxiter == 0) {
	minima_1.maxiter = 1000;
    }
    if (minima_1.nextiter == 0) {
	minima_1.nextiter = 1;
    }
    newhess = 1;
    maxerr = 3;
    done = FALSE_;
    s_copy(status, "          ", (ftnlen)9, (ftnlen)10);
    negtest = TRUE_;
    automode = FALSE_;
    automatic = FALSE_;
    if (s_cmp(mode, "AUTO", (ftnlen)6, (ftnlen)4) == 0) {
	automode = TRUE_;
    }
    if (s_cmp(method, "AUTO", (ftnlen)6, (ftnlen)4) == 0) {
	automatic = TRUE_;
    }

/*     set default parameters for the line search */

    if (linmin_1.cappa == 0.) {
	linmin_1.cappa = .1;
    }
    if (linmin_1.stpmin == 0.) {
	linmin_1.stpmin = 1e-20;
    }
    if (linmin_1.stpmax == 0.) {
	linmin_1.stpmax = 5.;
    }
    linmin_1.slpmax = 1e4;
    linmin_1.angmax = 180.;
    linmin_1.intmax = 8;

/*     search each line of the keyword file for options */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "FCTMIN ", (ftnlen)7, (ftnlen)7) == 0) {
	    i__2 = s_rsli(&io___15);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&minima_1.fctmin, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "MAXITER ", (ftnlen)8, (ftnlen)8) == 0) {
	    i__2 = s_rsli(&io___16);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&minima_1.maxiter, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "NEXTITER ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___17);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&minima_1.nextiter, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "NEWHESS ", (ftnlen)8, (ftnlen)8) == 0) {
	    i__2 = s_rsli(&io___18);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&newhess, (ftnlen)sizeof(
		    integer));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "SADDLEPOINT ", (ftnlen)12, (ftnlen)12) == 
		0) {
	    negtest = FALSE_;
	} else if (s_cmp(keyword, "STEPMIN ", (ftnlen)8, (ftnlen)8) == 0) {
	    i__2 = s_rsli(&io___19);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&linmin_1.stpmin, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "STEPMAX ", (ftnlen)8, (ftnlen)8) == 0) {
	    i__2 = s_rsli(&io___20);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&linmin_1.stpmax, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "CAPPA ", (ftnlen)6, (ftnlen)6) == 0) {
	    i__2 = s_rsli(&io___21);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&linmin_1.cappa, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "SLOPEMAX ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___22);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&linmin_1.slpmax, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "ANGMAX ", (ftnlen)7, (ftnlen)7) == 0) {
	    i__2 = s_rsli(&io___23);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&linmin_1.angmax, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "INTMAX ", (ftnlen)7, (ftnlen)7) == 0) {
	    i__2 = s_rsli(&io___24);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&linmin_1.intmax, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	}
L20:
	;
    }

/*     initialize the function call and iteration counters */

    fg_call__ = 0;
    nerr = 0;
    iter_tn__ = minima_1.nextiter - 1;
    minima_1.maxiter = iter_tn__ + minima_1.maxiter;

/*     print header information about the method used */

    if (inform_1.iprint > 0) {
	if (s_cmp(mode, "NEWTON", (ftnlen)6, (ftnlen)6) == 0) {
	    io___28.ciunit = iounit_1.iout;
	    s_wsfe(&io___28);
	    e_wsfe();
	} else if (s_cmp(mode, "TNCG", (ftnlen)6, (ftnlen)4) == 0) {
	    io___29.ciunit = iounit_1.iout;
	    s_wsfe(&io___29);
	    e_wsfe();
	} else if (s_cmp(mode, "DTNCG", (ftnlen)6, (ftnlen)5) == 0) {
	    io___30.ciunit = iounit_1.iout;
	    s_wsfe(&io___30);
	    e_wsfe();
	} else if (s_cmp(mode, "AUTO", (ftnlen)6, (ftnlen)4) == 0) {
	    io___31.ciunit = iounit_1.iout;
	    s_wsfe(&io___31);
	    e_wsfe();
	}
	io___32.ciunit = iounit_1.iout;
	s_wsfe(&io___32);
	do_fio(&c__1, mode, (ftnlen)6);
	do_fio(&c__1, method, (ftnlen)6);
	do_fio(&c__1, (char *)&(*grdmin), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___33.ciunit = iounit_1.iout;
	s_wsfe(&io___33);
	e_wsfe();
    }

/*     evaluate the function and get the initial gradient */

    iter_cg__ = 0;
    ++fg_call__;
    f = (*fgvalue)(&x[1], g);
    f_old__ = f;
    g_norm__ = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x_old__[i__ - 1] = x[i__];
/* Computing 2nd power */
	d__1 = g[i__ - 1];
	g_norm__ += d__1 * d__1;
    }
    g_norm__ = sqrt(g_norm__);
    f_move__ = linmin_1.stpmax * .5 * g_norm__;
    g_rms__ = g_norm__ / rms;

/*     print initial information prior to first iteration */

    if (inform_1.iprint > 0) {
	if (f < 1e6 && f > -1e5 && g_rms__ < 1e6) {
	    io___42.ciunit = iounit_1.iout;
	    s_wsfe(&io___42);
	    do_fio(&c__1, (char *)&iter_tn__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g_rms__, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&fg_call__, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___43.ciunit = iounit_1.iout;
	    s_wsfe(&io___43);
	    do_fio(&c__1, (char *)&iter_tn__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g_rms__, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&fg_call__, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }

/*     write initial intermediate prior to first iteration */

    if (inform_1.iwrite > 0) {
	(*optsave)(&iter_tn__, &f, &x[1]);
    }

/*     check for termination criteria met by initial point */

    if (g_rms__ <= *grdmin) {
	done = TRUE_;
	*minimum = f;
	if (inform_1.iprint > 0) {
	    io___44.ciunit = iounit_1.iout;
	    s_wsfe(&io___44);
	    e_wsfe();
	}
    } else if (f <= minima_1.fctmin) {
	done = TRUE_;
	*minimum = f;
	if (inform_1.iprint > 0) {
	    io___45.ciunit = iounit_1.iout;
	    s_wsfe(&io___45);
	    e_wsfe();
	}
    } else if (iter_tn__ >= minima_1.maxiter) {
	done = TRUE_;
	*minimum = f;
	if (inform_1.iprint > 0) {
	    io___46.ciunit = iounit_1.iout;
	    s_wsfe(&io___46);
	    e_wsfe();
	}
    }

/*     beginning of the outer truncated Newton iteration */

    while(! done) {
	++iter_tn__;

/*     if pisystem is present, update the molecular orbitals */

	if (potent_1.use_orbit__) {
	    piorbs_1.reorbit = 1;
	    piscf_();
	    ++fg_call__;
	    f = (*fgvalue)(&x[1], g);
	    piorbs_1.reorbit = 0;
	}

/*     choose the optimization mode based on the gradient value */

	if (automode) {
	    if (g_rms__ >= 3.) {
		s_copy(mode, "TNCG", (ftnlen)6, (ftnlen)4);
	    } else {
		s_copy(mode, "DTNCG", (ftnlen)6, (ftnlen)5);
	    }
	}

/*     decide on an optimal preconditioning based on the gradient */

	if (automatic) {
	    if (*nvar < 10) {
		s_copy(method, "DIAG", (ftnlen)6, (ftnlen)4);
		hescut_1.hesscut = 0.;
	    } else if (g_rms__ >= 10.) {
		s_copy(method, "DIAG", (ftnlen)6, (ftnlen)4);
		hescut_1.hesscut = 1.;
	    } else if (g_rms__ >= 1.) {
		s_copy(method, "ICCG", (ftnlen)6, (ftnlen)4);
		hescut_1.hesscut = *nvar * .001;
		if (hescut_1.hesscut > 1.) {
		    hescut_1.hesscut = 1.;
		}
	    } else {
		s_copy(method, "ICCG", (ftnlen)6, (ftnlen)4);
		hescut_1.hesscut = *nvar * .001;
		if (hescut_1.hesscut > .1) {
		    hescut_1.hesscut = .1;
		}
	    }
	}

/*     compute needed portions of the Hessian matrix */

	s_copy(h_mode__, "FULL", (ftnlen)4, (ftnlen)4);
	if ((iter_tn__ - 1) % newhess != 0) {
	    s_copy(h_mode__, "NONE", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(mode, "DTNCG", (ftnlen)6, (ftnlen)5) == 0 && s_cmp(method, 
		"NONE", (ftnlen)6, (ftnlen)4) == 0) {
	    s_copy(h_mode__, "NONE", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(mode, "DTNCG", (ftnlen)6, (ftnlen)5) == 0 && s_cmp(method, 
		"DIAG", (ftnlen)6, (ftnlen)4) == 0) {
	    s_copy(h_mode__, "DIAG", (ftnlen)4, (ftnlen)4);
	}
	(*hmatrix)(h_mode__, &x[1], h__, h_init__, h_stop__, h_index__, 
		h_diag__, (ftnlen)4);

/*     find the next approximate Newton search direction */

	tnsolve_(mode, method, &negtest, nvar, p, &x[1], g, h__, h_init__, 
		h_stop__, h_index__, h_diag__, &iter_tn__, &iter_cg__, &
		fg_call__, (D_fp)fgvalue, info_solve__, (ftnlen)6, (ftnlen)6, 
		(ftnlen)9);

/*     perform a line search in the chosen direction */

	s_copy(info_search__, "         ", (ftnlen)9, (ftnlen)9);
	search_(nvar, &f, g, &x[1], p, &f_move__, &angle, &fg_call__, (D_fp)
		fgvalue, info_search__, (ftnlen)9);
	if (s_cmp(info_search__, " Success ", (ftnlen)9, (ftnlen)9) != 0) {
	    s_copy(info_solve__, info_search__, (ftnlen)9, (ftnlen)9);
	}

/*     update variables to reflect this iteration */

	f_move__ = f_old__ - f;
	f_old__ = f;
	x_move__ = 0.;
	g_norm__ = 0.;
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = x[i__] - x_old__[i__ - 1];
	    x_move__ += d__1 * d__1;
	    x_old__[i__ - 1] = x[i__];
/* Computing 2nd power */
	    d__1 = g[i__ - 1];
	    g_norm__ += d__1 * d__1;
	}
	x_move__ = sqrt(x_move__);
	x_move__ /= rms;
	if (s_cmp(output_1.coordtype, "INTERNAL", (ftnlen)9, (ftnlen)8) == 0) 
		{
	    x_move__ *= 57.29577951308232088;
	}
	g_norm__ = sqrt(g_norm__);
	g_rms__ = g_norm__ / rms;

/*     quit if the maximum number of iterations is exceeded */

	if (iter_tn__ >= minima_1.maxiter) {
	    done = TRUE_;
	    s_copy(status, "IterLimit", (ftnlen)9, (ftnlen)9);
	}

/*     quit if the function value did not change */

	if (f_move__ == 0.) {
	    done = TRUE_;
	    s_copy(status, "NoMotion ", (ftnlen)9, (ftnlen)9);
	}

/*     quit if either of the normal termination tests are met */

	if (g_rms__ <= *grdmin) {
	    done = TRUE_;
	    s_copy(status, "SmallGrad", (ftnlen)9, (ftnlen)9);
	} else if (f <= minima_1.fctmin) {
	    done = TRUE_;
	    s_copy(status, "SmallFct ", (ftnlen)9, (ftnlen)9);
	}

/*     quit if the line search encounters successive problems */

	if (s_cmp(info_search__, "BadIntpln", (ftnlen)9, (ftnlen)9) == 0 || 
		s_cmp(info_search__, "IntplnErr", (ftnlen)9, (ftnlen)9) == 0) 
		{
	    ++nerr;
	    if (nerr >= maxerr) {
		done = TRUE_;
		s_copy(status, info_search__, (ftnlen)9, (ftnlen)9);
	    }
	} else {
	    nerr = 0;
	}

/*     print intermediate results for the current iteration */

	if (inform_1.iprint > 0) {
	    if (done || iter_tn__ % inform_1.iprint == 0) {
		if (f < 1e6 && f > -1e5 && g_rms__ < 1e6 && f_move__ < 1e5) {
		    io___58.ciunit = iounit_1.iout;
		    s_wsfe(&io___58);
		    do_fio(&c__1, (char *)&iter_tn__, (ftnlen)sizeof(integer))
			    ;
		    do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&g_rms__, (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&f_move__, (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&x_move__, (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&iter_cg__, (ftnlen)sizeof(integer))
			    ;
		    do_fio(&c__1, info_solve__, (ftnlen)9);
		    do_fio(&c__1, (char *)&fg_call__, (ftnlen)sizeof(integer))
			    ;
		    e_wsfe();
		} else {
		    io___59.ciunit = iounit_1.iout;
		    s_wsfe(&io___59);
		    do_fio(&c__1, (char *)&iter_tn__, (ftnlen)sizeof(integer))
			    ;
		    do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&g_rms__, (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&f_move__, (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&x_move__, (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&iter_cg__, (ftnlen)sizeof(integer))
			    ;
		    do_fio(&c__1, info_solve__, (ftnlen)9);
		    do_fio(&c__1, (char *)&fg_call__, (ftnlen)sizeof(integer))
			    ;
		    e_wsfe();
		}
	    }
	}

/*     write intermediate results for the current iteration */

	if (inform_1.iwrite > 0) {
	    if (done || iter_tn__ % inform_1.iwrite == 0) {
		(*optsave)(&iter_tn__, &f, &x[1]);
	    }
	}

/*     print the reason for terminating the optimization */

	if (done) {
	    *minimum = f;
	    if (inform_1.iprint > 0) {
		if (g_rms__ <= *grdmin || f <= minima_1.fctmin) {
		    io___60.ciunit = iounit_1.iout;
		    s_wsfe(&io___60);
		    do_fio(&c__1, status, (ftnlen)9);
		    e_wsfe();
		} else {
		    io___61.ciunit = iounit_1.iout;
		    s_wsfe(&io___61);
		    do_fio(&c__1, status, (ftnlen)9);
		    e_wsfe();
		}
	    }
	}
    }
    return 0;
} /* tncg_ */

#undef keyline_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine tnsolve  --  approx linear equation solution   ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "tnsolve" uses a linear conjugate gradient method to find */
/*     an approximate solution to the set of linear equations */
/*     represented in matrix form by Hp = -g (Newton's equations) */

/*     status codes upon return: */

/*     TruncNewt    convergence to (truncated) Newton criterion */
/*     NegCurve     termination upon detecting negative curvature */
/*     OverLimit    maximum number of CG iterations exceeded */


/* Subroutine */ int tnsolve_(char *mode, char *method, logical *negtest, 
	integer *nvar, doublereal *p, doublereal *x, doublereal *g, 
	doublereal *h__, integer *h_init__, integer *h_stop__, integer *
	h_index__, doublereal *h_diag__, integer *cycle, integer *iter_cg__, 
	integer *fg_call__, D_fp fgvalue, char *status, ftnlen mode_len, 
	ftnlen method_len, ftnlen status_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);
    integer i_dnnt(doublereal *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal converge, d__[75000];
    static integer i__, j, k;
    static doublereal m[75000], q[75000], r__[75000], s[75000], dd, gg, hj, 
	    dq, rr, rs, eps, beta;
    static integer iter;
    static doublereal alpha, delta, sigma, g_rms__, g_norm__, rs_new__, 
	    r_norm__, f_sigma__, g_sigma__[75000], x_sigma__[75000];
    extern /* Subroutine */ int precond_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, ftnlen);
    static integer maxiter;



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




/*     transformation using exact Hessian diagonal */

    /* Parameter adjustments */
    --h_diag__;
    --h_index__;
    --h_stop__;
    --h_init__;
    --h__;
    --g;
    --x;
    --p;

    /* Function Body */
    if (s_cmp(mode, "DTNCG", (ftnlen)6, (ftnlen)5) != 0 && s_cmp(method, 
	    "NONE", (ftnlen)6, (ftnlen)4) != 0) {
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    m[i__ - 1] = 1. / sqrt((d__1 = h_diag__[i__], abs(d__1)));
	}
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    g[i__] *= m[i__ - 1];
	    h_diag__[i__] = h_diag__[i__] * m[i__ - 1] * m[i__ - 1];
	    i__2 = h_stop__[i__];
	    for (j = h_init__[i__]; j <= i__2; ++j) {
		k = h_index__[j];
		h__[j] = h__[j] * m[i__ - 1] * m[k - 1];
	    }
	}
    }

/*     setup prior to linear conjugate gradient iterations */

    iter = 0;
    gg = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = 0.;
	r__[i__ - 1] = -g[i__];
	gg += g[i__] * g[i__];
    }
    g_norm__ = sqrt(gg);
    precond_(method, &iter, nvar, s, r__, &h__[1], &h_init__[1], &h_stop__[1],
	     &h_index__[1], &h_diag__[1], (ftnlen)6);
    rs = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__ - 1] = s[i__ - 1];
	rs += r__[i__ - 1] * s[i__ - 1];
    }
    if (s_cmp(mode, "NEWTON", (ftnlen)6, (ftnlen)6) == 0) {
	eps = 1e-10;
	maxiter = *nvar;
    } else if (s_cmp(mode, "TNCG", (ftnlen)6, (ftnlen)4) == 0 || s_cmp(mode, 
	    "DTNCG", (ftnlen)6, (ftnlen)5) == 0) {
	delta = 1.;
	eps = delta / (doublereal) (*cycle);
	g_rms__ = g_norm__ / sqrt((doublereal) (*nvar));
	eps = min(eps,g_rms__);
	converge = 1.;
	eps = pow_dd(&eps, &converge);
	d__1 = sqrt((doublereal) (*nvar)) * 10.;
	maxiter = i_dnnt(&d__1);
    }
    iter = 1;

/*     evaluate or estimate the matrix-vector product */

    while(TRUE_) {
	if (s_cmp(mode, "TNCG", (ftnlen)6, (ftnlen)4) == 0 || s_cmp(mode, 
		"NEWTON", (ftnlen)6, (ftnlen)6) == 0) {
	    i__1 = *nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		q[i__ - 1] = 0.;
	    }
	    i__1 = *nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		q[i__ - 1] += h_diag__[i__] * d__[i__ - 1];
		i__2 = h_stop__[i__];
		for (j = h_init__[i__]; j <= i__2; ++j) {
		    k = h_index__[j];
		    hj = h__[j];
		    q[i__ - 1] += hj * d__[k - 1];
		    q[k - 1] += hj * d__[i__ - 1];
		}
	    }
	} else if (s_cmp(mode, "DTNCG", (ftnlen)6, (ftnlen)5) == 0) {
	    dd = 0.;
	    i__1 = *nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dd += d__[i__ - 1] * d__[i__ - 1];
	    }
	    sigma = 1e-7 / sqrt(dd);
	    if (s_cmp(output_1.coordtype, "INTERNAL", (ftnlen)9, (ftnlen)8) ==
		     0) {
		sigma = 1e-4 / sqrt(dd);
	    }
	    i__1 = *nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		x_sigma__[i__ - 1] = x[i__] + sigma * d__[i__ - 1];
	    }
	    ++(*fg_call__);
	    f_sigma__ = (*fgvalue)(x_sigma__, g_sigma__);
	    i__1 = *nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		q[i__ - 1] = (g_sigma__[i__ - 1] - g[i__]) / sigma;
	    }
	}

/*     check for a direction of negative curvature */

	dq = 0.;
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dq += d__[i__ - 1] * q[i__ - 1];
	}
	if (*negtest) {
	    if (dq <= 0.) {
		if (iter == 1) {
		    i__1 = *nvar;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			p[i__] = d__[i__ - 1];
		    }
		}
		s_copy(status, " NegCurve", (ftnlen)9, (ftnlen)9);
		goto L10;
	    }
	}

/*     test the truncated Newton termination criterion */

	alpha = rs / dq;
	rr = 0.;
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    p[i__] += alpha * d__[i__ - 1];
	    r__[i__ - 1] -= alpha * q[i__ - 1];
	    rr += r__[i__ - 1] * r__[i__ - 1];
	}
	r_norm__ = sqrt(rr);
	if (r_norm__ / g_norm__ <= eps) {
	    s_copy(status, "TruncNewt", (ftnlen)9, (ftnlen)9);
	    goto L10;
	}

/*     solve the preconditioning equations */

	precond_(method, &iter, nvar, s, r__, &h__[1], &h_init__[1], &
		h_stop__[1], &h_index__[1], &h_diag__[1], (ftnlen)6);

/*     update the truncated Newton direction */

	rs_new__ = 0.;
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    rs_new__ += r__[i__ - 1] * s[i__ - 1];
	}
	beta = rs_new__ / rs;
	rs = rs_new__;
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d__[i__ - 1] = s[i__ - 1] + beta * d__[i__ - 1];
	}

/*     check for overlimit, then begin next iteration */

	if (iter >= maxiter) {
	    s_copy(status, "OverLimit", (ftnlen)9, (ftnlen)9);
	    goto L10;
	}
	++iter;
    }

/*     retransform and increment total iterations, then terminate */

L10:
    if (s_cmp(mode, "DTNCG", (ftnlen)6, (ftnlen)5) != 0 && s_cmp(method, 
	    "NONE", (ftnlen)6, (ftnlen)4) != 0) {
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    p[i__] *= m[i__ - 1];
	    g[i__] /= m[i__ - 1];
	}
    }
    *iter_cg__ += iter;
    return 0;
} /* tnsolve_ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine precond  --  precondition linear CG method  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "precond" solves a simplified version of the Newton equations */
/*     Ms = r, and uses the result to precondition linear conjugate */
/*     gradient iterations on the full Newton equations in "tnsolve" */

/*     reference for incomplete Cholesky factorization : */

/*     T. A. Manteuffel, "An Incomplete Factorization Technique */
/*     for Positive Definite Linear Systems", Mathematics of */
/*     Computation, 34, 473-497 (1980); the present method is */
/*     based upon the SICCG(0) method described in this paper */

/*     types of preconditioning methods : */

/*     none     use no preconditioning at all */
/*     diag     exact Hessian diagonal preconditioning */
/*     block    3x3 block diagonal preconditioning */
/*     ssor     symmetric successive over-relaxation */
/*     iccg     shifted incomplete Cholesky factorization */


/* Subroutine */ int precond_(char *method, integer *iter, integer *nvar, 
	doublereal *s, doublereal *r__, doublereal *h__, integer *h_init__, 
	integer *h_stop__, integer *h_index__, doublereal *h_diag__, ftnlen 
	method_len)
{
    /* Format strings */
    static char fmt_20[] = "(\002 PRECOND  --  Incomplete Cholesky is\002"
	    ",\002 Unstable, using Diagonal Method\002)";
    static char fmt_30[] = "(\002 PRECOND  --  Incomplete Cholesky\002,i12"
	    ",\002 Operations\002,f8.3,\002 Alpha Value\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal maxalpha;
    extern /* Subroutine */ int cholesky_(integer *, doublereal *, doublereal 
	    *);
    static doublereal a[6], b[3], f[1000000];
    static integer i__, j, k, ii, kk, ix, iy, iz;
    static doublereal f_i__, f_k__;
    static integer iii, kkk;
    static doublereal diag[75000], alpha, omega, f_diag__[75000];
    static integer c_init__[75000], nblock;
    static logical stable;
    static doublereal factor;
    static integer c_stop__[75000], icount;
    extern /* Subroutine */ int column_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static integer c_index__[1000000], c_value__[1000000];

    /* Fortran I/O blocks */
    static cilist io___111 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___120 = { 0, 0, 0, fmt_30, 0 };




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




/*     use no preconditioning, using M = identity matrix */

    /* Parameter adjustments */
    --h_diag__;
    --h_index__;
    --h_stop__;
    --h_init__;
    --h__;
    --r__;
    --s;

    /* Function Body */
    if (s_cmp(method, "NONE", (ftnlen)6, (ftnlen)4) == 0) {
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s[i__] = r__[i__];
	}
    }

/*     diagonal preconditioning, using M = abs(Hessian diagonal) */

    if (s_cmp(method, "DIAG", (ftnlen)6, (ftnlen)4) == 0) {
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s[i__] = r__[i__] / (d__1 = h_diag__[i__], abs(d__1));
	}
    }

/*     block diagonal preconditioning with exact atom blocks */
/*     (using M = 3x3 blocks from diagonal of full Hessian) */

    if (s_cmp(method, "BLOCK", (ftnlen)6, (ftnlen)5) == 0) {
	nblock = 3;
	i__1 = *nvar / 3;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iz = i__ * 3;
	    iy = iz - 1;
	    ix = iz - 2;
	    a[0] = h_diag__[ix];
	    if (h_index__[h_init__[ix]] == iy) {
		a[1] = h__[h_init__[ix]];
	    } else {
		a[1] = 0.;
	    }
	    if (h_index__[h_init__[ix] + 1] == iz) {
		a[2] = h__[h_init__[ix] + 1];
	    } else {
		a[2] = 0.;
	    }
	    a[3] = h_diag__[iy];
	    if (h_index__[h_init__[iy]] == iz) {
		a[4] = h__[h_init__[iy]];
	    } else {
		a[4] = 0.;
	    }
	    a[5] = h_diag__[iz];
	    b[0] = r__[ix];
	    b[1] = r__[iy];
	    b[2] = r__[iz];
	    cholesky_(&nblock, a, b);
	    s[ix] = b[0];
	    s[iy] = b[1];
	    s[iz] = b[2];
	}
    }

/*     symmetric successive over-relaxation (SSOR) preconditioning */
/*     (using M = (D/w+U)T * (D/w)-1 * (D/w+U) with 0 < w < 2) */

    if (s_cmp(method, "SSOR", (ftnlen)6, (ftnlen)4) == 0) {
	omega = 1.;
	factor = 2. - omega;
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s[i__] = r__[i__] * factor;
	    diag[i__ - 1] = h_diag__[i__] / omega;
	}
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s[i__] /= diag[i__ - 1];
	    i__2 = h_stop__[i__];
	    for (j = h_init__[i__]; j <= i__2; ++j) {
		k = h_index__[j];
		s[k] -= h__[j] * s[i__];
	    }
	}
	for (i__ = *nvar; i__ >= 1; --i__) {
	    s[i__] *= diag[i__ - 1];
	    i__1 = h_stop__[i__];
	    for (j = h_init__[i__]; j <= i__1; ++j) {
		k = h_index__[j];
		s[i__] -= h__[j] * s[k];
	    }
	    s[i__] /= diag[i__ - 1];
	}
    }

/*     factorization phase of incomplete cholesky preconditioning */

    if (s_cmp(method, "ICCG", (ftnlen)6, (ftnlen)4) == 0 && *iter == 0) {
	column_(nvar, &h_init__[1], &h_stop__[1], &h_index__[1], c_init__, 
		c_stop__, c_index__, c_value__);
	stable = TRUE_;
	icount = 0;
	maxalpha = 2.1;
	alpha = -.001;
L10:
	if (alpha <= 0.) {
	    alpha += .001;
	} else {
	    alpha *= 2.;
	}
	if (alpha > maxalpha) {
	    stable = FALSE_;
	    if (inform_1.verbose) {
		io___111.ciunit = iounit_1.iout;
		s_wsfe(&io___111);
		e_wsfe();
	    }
	} else {
	    factor = alpha + 1.;
	    i__1 = *nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		f_diag__[i__ - 1] = factor * h_diag__[i__];
		i__2 = c_stop__[i__ - 1];
		for (j = c_init__[i__ - 1]; j <= i__2; ++j) {
		    k = c_index__[j - 1];
		    f_i__ = f[c_value__[j - 1] - 1];
		    f_diag__[i__ - 1] -= f_i__ * f_i__ * f_diag__[k - 1];
		    ++icount;
		}
		if (f_diag__[i__ - 1] <= 0.) {
		    goto L10;
		}
		if (f_diag__[i__ - 1] < 1e-7) {
		    f_diag__[i__ - 1] = 1e-7;
		}
		f_diag__[i__ - 1] = 1. / f_diag__[i__ - 1];
		i__2 = h_stop__[i__];
		for (j = h_init__[i__]; j <= i__2; ++j) {
		    k = h_index__[j];
		    f[j - 1] = h__[j];
		    ii = c_init__[i__ - 1];
		    kk = c_init__[k - 1];
		    while(ii <= c_stop__[i__ - 1] && kk <= c_stop__[k - 1]) {
			iii = c_index__[ii - 1];
			kkk = c_index__[kk - 1];
			if (iii < kkk) {
			    ++ii;
			} else if (kkk < iii) {
			    ++kk;
			} else {
			    f_i__ = f[c_value__[ii - 1] - 1];
			    f_k__ = f[c_value__[kk - 1] - 1];
			    f[j - 1] -= f_i__ * f_k__ * f_diag__[iii - 1];
			    ++ii;
			    ++kk;
			    ++icount;
			}
		    }
		}
	    }
	    if (inform_1.verbose) {
		io___120.ciunit = iounit_1.iout;
		s_wsfe(&io___120);
		do_fio(&c__1, (char *)&icount, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&alpha, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
    }

/*     solution phase of incomplete cholesky preconditioning */

    if (s_cmp(method, "ICCG", (ftnlen)6, (ftnlen)4) == 0) {
	if (stable) {
	    i__1 = *nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		s[i__] = r__[i__];
	    }
	    i__1 = *nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		s[i__] *= f_diag__[i__ - 1];
		i__2 = h_stop__[i__];
		for (j = h_init__[i__]; j <= i__2; ++j) {
		    k = h_index__[j];
		    s[k] -= f[j - 1] * s[i__];
		}
	    }
	    for (i__ = *nvar; i__ >= 1; --i__) {
		s[i__] /= f_diag__[i__ - 1];
		i__1 = h_stop__[i__];
		for (j = h_init__[i__]; j <= i__1; ++j) {
		    k = h_index__[j];
		    s[i__] -= f[j - 1] * s[k];
		}
		s[i__] *= f_diag__[i__ - 1];
	    }
	} else {
	    i__1 = *nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		s[i__] = r__[i__] / (d__1 = h_diag__[i__], abs(d__1));
	    }
	}
    }
    return 0;
} /* precond_ */

