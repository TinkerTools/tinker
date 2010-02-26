/* ocvm.f -- translated by f2c (version 20050501).
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
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

struct {
    doublereal scale[75000];
    logical set_scale__;
} scales_;

#define scales_1 scales_

/* Table of constant values */

static integer c__2 = 2;
static integer c__5 = 5;
static integer c__1 = 1;
static integer c__3 = 3;
static real c_b96 = 0.f;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine ocvm  --  variable metric optimization method  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "ocvm" is an optimally conditioned variable metric nonlinear */
/*     optimization routine without line searches */

/*     literature references: */

/*     W. C. Davidon, "Optimally Conditioned Optimization Algorithms */
/*     Without Line Searches", Mathematical Programming, 9, 1-30 (1975) */

/*     D. F. Shanno and K-H. Phua, "Matrix Conditioning and Nonlinear */
/*     Optimization", Mathematical Programming, 14, 149-16 (1977) */

/*     D. F. Shanno and K-H. Phua, "Numerical Comparison of Several */
/*     Variable-Metric Algorithms", Journal of Optimization Theory */
/*     and Applications, 25, 507-518 (1978) */

/*     variables and parameters: */

/*     nvar      number of parameters in the objective function */
/*     x0        contains starting point upon input, upon return */
/*                 contains the best point found */
/*     f0        during optimization contains best current function */
/*                 value; returns final best function value */
/*     grdmin    normal exit if rms gradient gets below this value */
/*     ncalls    total number of function/gradient evaluations */

/*     required external routines: */

/*     fgvalue    function to evaluate function and gradient values */
/*     optsave    subroutine to write out info about current status */


/* Subroutine */ int ocvm_(integer *nvar, doublereal *x0, doublereal *f0, 
	doublereal *grdmin, D_fp fgvalue, S_fp optsave)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 OCVM  --  Too many Parameters,\002,\002 "
	    "Increase the Value of MAXOPT\002)";
    static char fmt_30[] = "(/,\002 Optimally Conditioned Variable Metric"
	    "\002,\002 Optimization :\002)";
    static char fmt_40[] = "(/,\002 VM Iter    F Value       G RMS     F M"
	    "ove\002,\002    X Move      Angle   FG Call\002,/)";
    static char fmt_50[] = "(i6,f13.4,f12.4,32x,i9)";
    static char fmt_60[] = "(i6,d13.4,d12.4,32x,i9)";
    static char fmt_70[] = "(i6,f13.4,f12.4,f11.4,f10.4,f11.4,i9)";
    static char fmt_80[] = "(i6,d13.4,d12.4,d11.4,f10.4,f11.4,i9)";
    static char fmt_90[] = "(i6,f13.4,f12.4,f11.4,f10.4,f11.4,i9)";
    static char fmt_100[] = "(i6,d13.4,d12.4,d11.4,f10.4,f11.4,i9)";
    static char fmt_110[] = "(/,\002 OCVM  --  Incomplete Convergence\002"
	    ",\002 due to \002,a9)";
    static char fmt_120[] = "(/,\002 OCVM  --  Normal Termination\002,\002 d"
	    "ue to \002,a9)";
    static char fmt_140[] = "(i6,f13.4,f12.4,f11.4,f10.4,f11.4,i9)";
    static char fmt_150[] = "(i6,d13.4,d12.4,f11.4,f10.4,f11.4,i9)";
    static char fmt_160[] = "(/,\002 OCVM  --  Incomplete Convergence\002"
	    ",\002 due to \002,a9)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    double sqrt(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), do_fio(integer *, char *, ftnlen);
    double acos(doublereal);

    /* Local variables */
    static doublereal srchnorm, a, b, c__, f, g[1000], h__[1000000]	/* 
	    was [1000][1000] */;
    static integer i__, j;
    static doublereal k[1000], m[1000], n[1000], p[1000], q[1000], s[1000], u[
	    1000], v, w[1000], x[1000], b0, k0[1000], m2, n2, u2, sg, hq[1000]
	    , mw, us, qk0, eps, rms;
    static integer nbig;
    static logical done;
    static doublereal zeta;
    static integer mvar;
    static doublereal grms;
    static integer next;
    static doublereal f0old, x0old[1000], gamma, alpha, delta, fmove;
    static integer niter;
    static doublereal gnorm;
    static integer nstep;
    static doublereal snorm, xmove, search[1000];
    static integer maxbig;
    static doublereal cosang;
    static integer ncalls;
    static char record[120];
    static doublereal fprime, micron;
    static char string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char status[9];
    static doublereal f0prime, sgangle;
    extern doublereal precise_(integer *);
    static integer maxstep;
    static logical restart;
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };
    static icilist io___14 = { 1, string, 1, 0, 120, 1 };
    static icilist io___15 = { 1, string, 1, 0, 120, 1 };
    static icilist io___16 = { 1, string, 1, 0, 120, 1 };
    static icilist io___17 = { 1, string, 1, 0, 120, 1 };
    static icilist io___18 = { 1, string, 1, 0, 120, 1 };
    static icilist io___19 = { 1, string, 1, 0, 120, 1 };
    static cilist io___25 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___56 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___57 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___58 = { 0, 0, 0, fmt_160, 0 };



#define h___ref(a_1,a_2) h__[(a_2)*1000 + a_1 - 1001]
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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  scales.i  --  parameter scale factors for optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     scale      multiplicative factor for each optimization parameter */
/*     set_scale  logical flag to show if scale factors have been set */




/*     initialization and set-up for the optimization */

    /* Parameter adjustments */
    --x0;

    /* Function Body */
    if (*nvar > 1000) {
	io___1.ciunit = iounit_1.iout;
	s_wsfe(&io___1);
	e_wsfe();
	return 0;
    }
    mvar = *nvar;
    rms = sqrt((doublereal) (*nvar));
    if (s_cmp(output_1.coordtype, "CARTESIAN", (ftnlen)9, (ftnlen)9) == 0) {
	rms /= sqrt(3.);
    } else if (s_cmp(output_1.coordtype, "RIGIDBODY", (ftnlen)9, (ftnlen)9) ==
	     0) {
	rms /= sqrt(6.);
    }
    maxbig = 2;
    maxstep = 10;
    eps = precise_(&c__2);
    restart = TRUE_;
    done = FALSE_;

/*     set default values for variable scale factors */

    if (! scales_1.set_scale__) {
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (scales_1.scale[i__ - 1] == 0.) {
		scales_1.scale[i__ - 1] = 1.;
	    }
	}
    }

/*     set default parameters for the optimization */

    if (minima_1.fctmin == 0.) {
	minima_1.fctmin = -1e6;
    }
    if (minima_1.maxiter == 0) {
	minima_1.maxiter = 1000000;
    }
    if (minima_1.nextiter == 0) {
	minima_1.nextiter = 1;
    }
    if (inform_1.iprint < 0) {
	inform_1.iprint = 1;
    }
    if (inform_1.iwrite < 0) {
	inform_1.iwrite = 1;
    }
    if (linmin_1.stpmax == 0.) {
	linmin_1.stpmax = 5.;
    }
    if (minima_1.hguess == 0.) {
	minima_1.hguess = .4;
    }
    linmin_1.angmax = 180.;

/*     search the keywords for optimization parameters */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "FCTMIN ", (ftnlen)7, (ftnlen)7) == 0) {
	    i__2 = s_rsli(&io___14);
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
	    i__2 = s_rsli(&io___15);
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
	    i__2 = s_rsli(&io___16);
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
	} else if (s_cmp(keyword, "HGUESS ", (ftnlen)7, (ftnlen)7) == 0) {
	    i__2 = s_rsli(&io___17);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&minima_1.hguess, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "STEPMAX ", (ftnlen)8, (ftnlen)8) == 0) {
	    i__2 = s_rsli(&io___18);
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
	} else if (s_cmp(keyword, "ANGMAX ", (ftnlen)7, (ftnlen)7) == 0) {
	    i__2 = s_rsli(&io___19);
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
	}
L20:
	;
    }

/*     evaluate the function and get the initial gradient */

    niter = minima_1.nextiter - 1;
    minima_1.maxiter = niter + minima_1.maxiter;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x0old[i__ - 1] = x0[i__];
    }
    ncalls = 1;
    *f0 = (*fgvalue)(&x0[1], g);
    f0old = *f0;

/*     print initial information prior to first iteration */

    if (inform_1.iprint > 0) {
	io___25.ciunit = iounit_1.iout;
	s_wsfe(&io___25);
	e_wsfe();
	io___26.ciunit = iounit_1.iout;
	s_wsfe(&io___26);
	e_wsfe();
    }

/*     set the "h" matrix to a diagonal upon restarting */

    while(! done) {
	if (restart) {
	    i__1 = mvar;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *nvar;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    h___ref(i__, j) = 0.;
		}
	    }
	    i__1 = mvar;
	    for (j = 1; j <= i__1; ++j) {
		h___ref(j, j) = minima_1.hguess;
	    }
	    i__1 = mvar;
	    for (j = 1; j <= i__1; ++j) {
		k0[j - 1] = 0.;
		i__2 = *nvar;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    k0[j - 1] += h___ref(i__, j) * g[i__ - 1];
		}
		w[j - 1] = k0[j - 1];
	    }
	    restart = FALSE_;
	}

/*     start the next iteration using either an updated "h" */
/*     matrix or the "h" matrix from the previous iteration */

	gnorm = 0.;
	grms = 0.;
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = g[i__ - 1];
	    gnorm += d__1 * d__1;
/* Computing 2nd power */
	    d__1 = g[i__ - 1] * scales_1.scale[i__ - 1];
	    grms += d__1 * d__1;
	}
	gnorm = sqrt(gnorm);
	grms = sqrt(grms) / rms;
	xmove = 0.;
	if (niter != 0) {
	    i__1 = *nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		d__1 = (x0[i__] - x0old[i__ - 1]) / scales_1.scale[i__ - 1];
		xmove += d__1 * d__1;
		x0old[i__ - 1] = x0[i__];
	    }
	    xmove = sqrt(xmove) / rms;
	    if (s_cmp(output_1.coordtype, "INTERNAL", (ftnlen)9, (ftnlen)8) ==
		     0) {
		xmove *= 57.29577951308232088;
	    }
	    fmove = f0old - *f0;
	    f0old = *f0;
	}

/*     print intermediate results for the current iteration */

	if (inform_1.iprint > 0) {
	    if (niter == 0) {
		if (*f0 < 1e7 && *f0 > -1e6 && grms < 1e6) {
		    io___35.ciunit = iounit_1.iout;
		    s_wsfe(&io___35);
		    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*f0), (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
		    e_wsfe();
		} else {
		    io___36.ciunit = iounit_1.iout;
		    s_wsfe(&io___36);
		    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*f0), (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
		    e_wsfe();
		}
	    } else if (niter % inform_1.iprint == 0) {
		if (*f0 < 1e7 && *f0 > -1e6 && grms < 1e6 && fmove < 1e5) {
		    io___37.ciunit = iounit_1.iout;
		    s_wsfe(&io___37);
		    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*f0), (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&fmove, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&xmove, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&sgangle, (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
		    e_wsfe();
		} else {
		    io___39.ciunit = iounit_1.iout;
		    s_wsfe(&io___39);
		    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*f0), (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&fmove, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&xmove, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&sgangle, (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
		    e_wsfe();
		}
	    }
	}

/*     write intermediate results for the current iteration */

	if (inform_1.iwrite > 0) {
	    if (niter % inform_1.iwrite == 0) {
		(*optsave)(&niter, f0, x);
	    }
	}

/*     before starting the next iteration, check to see whether */
/*     the gradient norm, function decrease or iteration limit */
/*     termination criteria have been satisfied */

	if (grms < *grdmin || *f0 < minima_1.fctmin || niter >= 
		minima_1.maxiter) {
	    if (inform_1.iprint > 0) {
		if (niter != 0 && niter % inform_1.iprint != 0) {
		    if (*f0 < 1e7 && *f0 > -1e6 && grms < 1e6 && fmove < 1e5) 
			    {
			io___41.ciunit = iounit_1.iout;
			s_wsfe(&io___41);
			do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer))
				;
			do_fio(&c__1, (char *)&(*f0), (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&fmove, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&xmove, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&sgangle, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer)
				);
			e_wsfe();
		    } else {
			io___42.ciunit = iounit_1.iout;
			s_wsfe(&io___42);
			do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer))
				;
			do_fio(&c__1, (char *)&(*f0), (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&fmove, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&xmove, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&sgangle, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer)
				);
			e_wsfe();
		    }
		}
		if (niter >= minima_1.maxiter) {
		    s_copy(status, "IterLimit", (ftnlen)9, (ftnlen)9);
		}
		if (*f0 < minima_1.fctmin) {
		    s_copy(status, "SmallFct ", (ftnlen)9, (ftnlen)9);
		}
		if (grms < *grdmin) {
		    s_copy(status, "SmallGrad", (ftnlen)9, (ftnlen)9);
		}
		if (s_cmp(status, "IterLimit", (ftnlen)9, (ftnlen)9) == 0) {
		    io___44.ciunit = iounit_1.iout;
		    s_wsfe(&io___44);
		    do_fio(&c__1, status, (ftnlen)9);
		    e_wsfe();
		} else {
		    io___45.ciunit = iounit_1.iout;
		    s_wsfe(&io___45);
		    do_fio(&c__1, status, (ftnlen)9);
		    e_wsfe();
		}
	    }
	    if (inform_1.iwrite > 0) {
		if (niter % inform_1.iwrite != 0) {
		    (*optsave)(&niter, f0, x);
		}
	    }
	    done = TRUE_;
	    goto L170;
	}

/*     start of the next iteration */

	++niter;
	sg = 0.;
	snorm = 0.;
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
	    s[j - 1] = -k0[j - 1];
/* Computing 2nd power */
	    d__1 = s[j - 1];
	    snorm += d__1 * d__1;
	    sg -= s[j - 1] * g[j - 1];
	}
	f0prime = -snorm;
	snorm = sqrt(snorm);
	cosang = sg / (snorm * gnorm);
/* Computing MIN */
	d__1 = 1., d__2 = max(-1.,cosang);
	cosang = min(d__1,d__2);
	sgangle = acos(cosang) * 57.29577951308232088;
	if (sgangle > linmin_1.angmax) {
	    ++nbig;
	} else {
	    nbig = 0;
	}
	zeta = 2.;
	if ((*f0 - minima_1.fctmin) * 4. < -f0prime) {
	    i__1 = mvar;
	    for (j = 1; j <= i__1; ++j) {
		s[j - 1] = -s[j - 1] * ((*f0 - minima_1.fctmin) * 4. / 
			f0prime);
	    }
	    f0prime = (*f0 - minima_1.fctmin) * -4.;
	}

/*     location of the next starting point */

	nstep = 0;
L130:
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    search[i__ - 1] = 0.;
	}
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *nvar;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		search[i__ - 1] += h___ref(i__, j) * s[j - 1];
	    }
	}
	srchnorm = 0.;
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = search[i__ - 1];
	    srchnorm += d__1 * d__1;
	}
	srchnorm = sqrt(srchnorm);
	if (srchnorm > linmin_1.stpmax) {
	    i__1 = mvar;
	    for (j = 1; j <= i__1; ++j) {
		s[j - 1] = linmin_1.stpmax / srchnorm * s[j - 1];
	    }
	    i__1 = *nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		search[i__ - 1] = linmin_1.stpmax / srchnorm * search[i__ - 1]
			;
	    }
	    f0prime = linmin_1.stpmax / srchnorm * f0prime;
	    zeta = .5;
	}

/*     invoke abnormal termination if -f0prime is too small */

	if (-f0prime < eps) {
	    if (inform_1.iprint > 0) {
		if (niter != 0 && niter % inform_1.iprint != 0) {
		    if (*f0 < 1e7 && *f0 > -1e6 && grms < 1e6) {
			io___56.ciunit = iounit_1.iout;
			s_wsfe(&io___56);
			do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer))
				;
			do_fio(&c__1, (char *)&(*f0), (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&c_b96, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c_b96, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&sgangle, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer)
				);
			e_wsfe();
		    } else {
			io___57.ciunit = iounit_1.iout;
			s_wsfe(&io___57);
			do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer))
				;
			do_fio(&c__1, (char *)&(*f0), (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&c_b96, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c_b96, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&sgangle, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer)
				);
			e_wsfe();
		    }
		}
		s_copy(status, "SmallMove", (ftnlen)9, (ftnlen)9);
		io___58.ciunit = iounit_1.iout;
		s_wsfe(&io___58);
		do_fio(&c__1, status, (ftnlen)9);
		e_wsfe();
	    }
	    if (inform_1.iwrite > 0) {
		if (niter % inform_1.iwrite != 0) {
		    (*optsave)(&niter, f0, x);
		}
	    }
	    done = TRUE_;
	    goto L170;
	}
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__ - 1] = x0[i__] + search[i__ - 1];
	}
	++ncalls;
	f = (*fgvalue)(x, g);
	if (f >= *f0) {
	    i__1 = mvar;
	    for (j = 1; j <= i__1; ++j) {
		s[j - 1] *= .5;
	    }
	    f0prime *= .5;
	    zeta = .5;
	    goto L130;
	}

/*     decide whether to update or take another step */

	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
	    k[j - 1] = 0.;
	    i__2 = *nvar;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		k[j - 1] += h___ref(i__, j) * g[i__ - 1];
	    }
	}
	fprime = 0.;
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
	    fprime += k[j - 1] * s[j - 1];
	}
	b0 = fprime - f0prime;
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
	    m[j - 1] = s[j - 1] + k0[j - 1] - k[j - 1];
	    k0[j - 1] = k[j - 1];
	}
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x0[i__] = x[i__ - 1];
	}
	*f0 = f;
	f0prime = fprime;
	if (b0 < eps) {
	    ++nstep;
	    if (nstep >= maxstep) {
		restart = TRUE_;
		goto L170;
	    }
	    i__1 = mvar;
	    for (j = 1; j <= i__1; ++j) {
		s[j - 1] *= zeta;
	    }
	    f0prime *= zeta;
	    goto L130;
	}

/*     check to see if we need to update */

	if (nbig >= maxbig) {
	    restart = TRUE_;
	    goto L170;
	}
	m2 = 0.;
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	    d__1 = m[j - 1];
	    m2 += d__1 * d__1;
	}
	if (m2 < eps) {
	    goto L170;
	}
	v = 0.;
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
	    v += m[j - 1] * s[j - 1];
	}
	micron = v - m2;
	mw = 0.;
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
	    mw += m[j - 1] * w[j - 1];
	}
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
	    u[j - 1] = w[j - 1] - m[j - 1] * (mw / m2);
	}
	u2 = 0.;
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	    d__1 = u[j - 1];
	    u2 += d__1 * d__1;
	}
	if (m2 * u2 >= eps) {
	    us = 0.;
	    i__1 = mvar;
	    for (j = 1; j <= i__1; ++j) {
		us += u[j - 1] * s[j - 1];
	    }
	    i__1 = mvar;
	    for (j = 1; j <= i__1; ++j) {
		n[j - 1] = u[j - 1] * (us / u2);
	    }
	    n2 = us * us / u2;
	} else {
	    i__1 = mvar;
	    for (j = 1; j <= i__1; ++j) {
		n[j - 1] = 0.;
	    }
	    n2 = 0.;
	}

/*     test inner product of projected s and del-g */

	b = n2 + micron * v / m2;
	if (b < eps) {
	    i__1 = mvar;
	    for (j = 1; j <= i__1; ++j) {
		n[j - 1] = s[j - 1] - m[j - 1] * (v / m2);
	    }
	    n2 = b0 - micron * v / m2;
	    b = b0;
	}

/*     set "gamma" and "delta" for the update */

	if (micron * v >= m2 * n2) {
	    gamma = 0.;
	    delta = sqrt(v / micron);
	} else {
	    a = b - micron;
	    c__ = b + v;
	    gamma = sqrt((1. - micron * v / (m2 * n2)) / (a * b));
	    delta = sqrt(c__ / a);
	    if (c__ < a) {
		gamma = -gamma;
	    }
	}

/*     perform the update of the "h" matrix */

	alpha = v + micron * delta + m2 * n2 * gamma;
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
	    p[j - 1] = m[j - 1] * (delta - n2 * gamma) + n[j - 1] * (gamma * 
		    v);
	    q[j - 1] = m[j - 1] * ((n2 * gamma + 1.) / alpha) - n[j - 1] * (
		    gamma * micron / alpha);
	    w[j - 1] = m[j - 1] * (n2 * (gamma * micron * v + 1.) / alpha) - 
		    n[j - 1] * ((delta + 1.) * micron * v / alpha);
	}
	qk0 = 0.;
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
	    qk0 += q[j - 1] * k0[j - 1];
	}
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
	    k0[j - 1] += p[j - 1] * qk0;
	}
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    hq[i__ - 1] = 0.;
	}
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *nvar;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		hq[i__ - 1] += h___ref(i__, j) * q[j - 1];
	    }
	}
	i__1 = mvar;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *nvar;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		h___ref(i__, j) = h___ref(i__, j) + hq[i__ - 1] * p[j - 1];
	    }
	}
	if (n2 <= 0.) {
	    i__1 = mvar;
	    for (j = 1; j <= i__1; ++j) {
		w[j - 1] = k0[j - 1];
	    }
	}
L170:
	;
    }
    return 0;
} /* ocvm_ */

#undef keyline_ref
#undef h___ref


