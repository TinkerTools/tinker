/* lbfgs.f -- translated by f2c (version 20050501).
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
    doublereal scale[75000];
    logical set_scale__;
} scales_;

#define scales_1 scales_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine lbfgs  --  limited memory BFGS optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "lbfgs" is a limited memory BFGS quasi-newton nonlinear */
/*     optimization routine */

/*     literature references: */

/*     J. Nocedal, "Updating Quasi-Newton Matrices with Limited */
/*     Storage", Mathematics of Computation, 35, 773-782 (1980) */

/*     D. C. Lui and J. Nocedal, "On the Limited Memory BFGS Method */
/*     for Large Scale Optimization", Mathematical Programming, */
/*     45, 503-528 (1989) */

/*     J. Nocedal and S. J. Wright, "Numerical Optimization", */
/*     Springer-Verlag New York, 1999, Section 9.1 */

/*     variables and parameters: */

/*     nvar      number of parameters in the objective function */
/*     x         contains starting point upon input, upon return */
/*                 contains the best point found */
/*     minimum   during optimization contains best current function */
/*                 value; returns final best function value */
/*     grdmin    normal exit if rms gradient gets below this value */
/*     ncalls    total number of function/gradient evaluations */

/*     required external routines: */

/*     fgvalue    function to evaluate function and gradient values */
/*     optsave    subroutine to write out info about current status */


/* Subroutine */ int lbfgs_(integer *nvar, doublereal *x, doublereal *minimum,
	 doublereal *grdmin, D_fp fgvalue, S_fp optsave)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 LBFGS  --  Too many Parameters,\002,\002"
	    " Increase the Value of MAXVAR\002)";
    static char fmt_30[] = "(/,\002 LBFGS  --  Number of Saved Vectors Set"
	    " to\002,\002 the Maximum of\002,i5)";
    static char fmt_40[] = "(/,\002 Steepest Descent Gradient Optimization "
	    ":\002)";
    static char fmt_50[] = "(/,\002 SD Iter    F Value      G RMS     F Mov"
	    "e\002,\002    X Move    Angle  FG Call  Comment\002,/)";
    static char fmt_60[] = "(/,\002 Limited Memory BFGS Quasi-Newton\002,"
	    "\002 Optimization :\002)";
    static char fmt_70[] = "(/,\002 QN Iter    F Value      G RMS     F Mov"
	    "e\002,\002    X Move    Angle  FG Call  Comment\002,/)";
    static char fmt_80[] = "(i6,f13.4,f11.4,30x,i7)";
    static char fmt_90[] = "(i6,d13.4,d11.4,30x,i7)";
    static char fmt_100[] = "(i6,f13.4,f11.4,f11.4,f10.4,f9.2,i7,3x,a9)";
    static char fmt_110[] = "(i6,d13.4,d11.4,d11.4,f10.4,f9.2,i7,3x,a9)";
    static char fmt_120[] = "(/,\002 LBFGS  --  Normal Termination due to"
	    " \002,a9)";
    static char fmt_130[] = "(/,\002 LBFGS  --  Incomplete Convergence due t"
	    "o \002,a9)";

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
    static doublereal f, g[75000];
    static integer i__, j, k, m;
    static doublereal p[75000], q[75000], r__[75000], s[1500000]	/* 
	    was [75000][20] */, y[1500000]	/* was [75000][20] */, h0[
	    75000], ys, yy, rho[20], rms, beta;
    static logical done;
    static integer msav, muse, nerr, next;
    static doublereal gamma, f_old__, g_old__[75000], angle, alpha[20];
    static char blank[9];
    static doublereal x_old__[75000], g_rms__;
    static integer niter;
    extern /* Subroutine */ int search_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , D_fp, char *, ftnlen);
    static integer ncalls;
    static doublereal f_move__;
    static integer maxerr;
    static doublereal x_move__, g_norm__;
    static char record[120], status[9], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
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
    static icilist io___20 = { 1, string, 1, 0, 120, 1 };
    static icilist io___21 = { 1, string, 1, 0, 120, 1 };
    static icilist io___22 = { 1, string, 1, 0, 120, 1 };
    static cilist io___23 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___58 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___59 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___60 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_130, 0 };



#define s_ref(a_1,a_2) s[(a_2)*75000 + a_1 - 75001]
#define y_ref(a_1,a_2) y[(a_2)*75000 + a_1 - 75001]
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
/*     ##  scales.i  --  parameter scale factors for optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     scale      multiplicative factor for each optimization parameter */
/*     set_scale  logical flag to show if scale factors have been set */




/*     initialize some values to be used below */

    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*nvar > 75000) {
	io___1.ciunit = iounit_1.iout;
	s_wsfe(&io___1);
	e_wsfe();
	return 0;
    }
    ncalls = 0;
    rms = sqrt((doublereal) (*nvar));
    if (s_cmp(output_1.coordtype, "CARTESIAN", (ftnlen)9, (ftnlen)9) == 0) {
	rms /= sqrt(3.);
    } else if (s_cmp(output_1.coordtype, "RIGIDBODY", (ftnlen)9, (ftnlen)9) ==
	     0) {
	rms /= sqrt(6.);
    }
    s_copy(blank, "         ", (ftnlen)9, (ftnlen)9);
    done = FALSE_;
    nerr = 0;
    maxerr = 2;

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

    msav = min(*nvar,20);
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

/*     set default parameters for the line search */

    if (linmin_1.stpmin == 0.) {
	linmin_1.stpmin = 1e-16;
    }
    if (linmin_1.stpmax == 0.) {
	linmin_1.stpmax = 5.;
    }
    if (linmin_1.cappa == 0.) {
	linmin_1.cappa = .9;
    }
    linmin_1.slpmax = 1e4;
    linmin_1.angmax = 180.;
    linmin_1.intmax = 5;

/*     search the keywords for optimization parameters */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "LBFGS-VECTORS ", (ftnlen)14, (ftnlen)14) == 0) {
	    i__2 = s_rsli(&io___14);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&msav, (ftnlen)sizeof(integer)
		    );
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "STEEPEST-DESCENT ", (ftnlen)17, (ftnlen)17)
		 == 0) {
	    msav = 0;
	} else if (s_cmp(keyword, "FCTMIN ", (ftnlen)7, (ftnlen)7) == 0) {
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
	} else if (s_cmp(keyword, "STEPMAX ", (ftnlen)8, (ftnlen)8) == 0) {
	    i__2 = s_rsli(&io___17);
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
	} else if (s_cmp(keyword, "STEPMIN ", (ftnlen)8, (ftnlen)8) == 0) {
	    i__2 = s_rsli(&io___18);
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
	} else if (s_cmp(keyword, "CAPPA ", (ftnlen)6, (ftnlen)6) == 0) {
	    i__2 = s_rsli(&io___19);
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
	    i__2 = s_rsli(&io___20);
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
	    i__2 = s_rsli(&io___21);
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
	    i__2 = s_rsli(&io___22);
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

/*     check the number of saved correction vectors */

    if (msav < 0 || msav > min(*nvar,20)) {
	msav = min(*nvar,20);
	io___23.ciunit = iounit_1.iout;
	s_wsfe(&io___23);
	do_fio(&c__1, (char *)&msav, (ftnlen)sizeof(integer));
	e_wsfe();
    }

/*     evaluate the function and get the initial gradient */

    niter = minima_1.nextiter - 1;
    minima_1.maxiter = niter + minima_1.maxiter;
    ++ncalls;
    f = (*fgvalue)(&x[1], g);
    f_old__ = f;
    m = 0;
    gamma = 1.;
    g_norm__ = 0.;
    g_rms__ = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g_norm__ += g[i__ - 1] * g[i__ - 1];
/* Computing 2nd power */
	d__1 = g[i__ - 1] * scales_1.scale[i__ - 1];
	g_rms__ += d__1 * d__1;
    }
    g_norm__ = sqrt(g_norm__);
    g_rms__ = sqrt(g_rms__) / rms;
    f_move__ = linmin_1.stpmax * .5 * g_norm__;

/*     print initial information prior to first iteration */

    if (inform_1.iprint > 0) {
	if (msav == 0) {
	    io___33.ciunit = iounit_1.iout;
	    s_wsfe(&io___33);
	    e_wsfe();
	    io___34.ciunit = iounit_1.iout;
	    s_wsfe(&io___34);
	    e_wsfe();
	} else {
	    io___35.ciunit = iounit_1.iout;
	    s_wsfe(&io___35);
	    e_wsfe();
	    io___36.ciunit = iounit_1.iout;
	    s_wsfe(&io___36);
	    e_wsfe();
	}
	if (f < 1e7 && f > -1e6 && g_rms__ < 1e5) {
	    io___37.ciunit = iounit_1.iout;
	    s_wsfe(&io___37);
	    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g_rms__, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___38.ciunit = iounit_1.iout;
	    s_wsfe(&io___38);
	    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g_rms__, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }

/*     write initial intermediate prior to first iteration */

    if (inform_1.iwrite > 0) {
	(*optsave)(&niter, &f, &x[1]);
    }

/*     tests of the various termination criteria */

    if (niter >= minima_1.maxiter) {
	s_copy(status, "IterLimit", (ftnlen)9, (ftnlen)9);
	done = TRUE_;
    }
    if (f <= minima_1.fctmin) {
	s_copy(status, "SmallFct ", (ftnlen)9, (ftnlen)9);
	done = TRUE_;
    }
    if (g_rms__ <= *grdmin) {
	s_copy(status, "SmallGrad", (ftnlen)9, (ftnlen)9);
	done = TRUE_;
    }

/*     start of a new limited memory BFGS iteration */

    while(! done) {
	++niter;
/* Computing MIN */
	i__1 = niter - 1;
	muse = min(i__1,msav);
	++m;
	if (m > msav) {
	    m = 1;
	}

/*     estimate Hessian diagonal and compute the Hg product */

	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    h0[i__ - 1] = gamma;
	    q[i__ - 1] = g[i__ - 1];
	}
	k = m;
	i__1 = muse;
	for (j = 1; j <= i__1; ++j) {
	    --k;
	    if (k == 0) {
		k = msav;
	    }
	    alpha[k - 1] = 0.;
	    i__2 = *nvar;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		alpha[k - 1] += s_ref(i__, k) * q[i__ - 1];
	    }
	    alpha[k - 1] *= rho[k - 1];
	    i__2 = *nvar;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		q[i__ - 1] -= alpha[k - 1] * y_ref(i__, k);
	    }
	}
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    r__[i__ - 1] = h0[i__ - 1] * q[i__ - 1];
	}
	i__1 = muse;
	for (j = 1; j <= i__1; ++j) {
	    beta = 0.;
	    i__2 = *nvar;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		beta += y_ref(i__, k) * r__[i__ - 1];
	    }
	    beta *= rho[k - 1];
	    i__2 = *nvar;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		r__[i__ - 1] += s_ref(i__, k) * (alpha[k - 1] - beta);
	    }
	    ++k;
	    if (k > msav) {
		k = 1;
	    }
	}

/*     set search direction and store current point and gradient */

	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    p[i__ - 1] = -r__[i__ - 1];
	    x_old__[i__ - 1] = x[i__];
	    g_old__[i__ - 1] = g[i__ - 1];
	}

/*     perform line search along the new conjugate direction */

	s_copy(status, blank, (ftnlen)9, (ftnlen)9);
	search_(nvar, &f, g, &x[1], p, &f_move__, &angle, &ncalls, (D_fp)
		fgvalue, status, (ftnlen)9);

/*     update variables based on results of this iteration */

	ys = 0.;
	yy = 0.;
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s_ref(i__, m) = x[i__] - x_old__[i__ - 1];
	    y_ref(i__, m) = g[i__ - 1] - g_old__[i__ - 1];
	    ys += y_ref(i__, m) * s_ref(i__, m);
	    yy += y_ref(i__, m) * y_ref(i__, m);
	}
	gamma = (d__1 = ys / yy, abs(d__1));
	rho[m - 1] = 1. / ys;

/*     get the sizes of the moves made during this iteration */

	f_move__ = f_old__ - f;
	f_old__ = f;
	x_move__ = 0.;
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = (x[i__] - x_old__[i__ - 1]) / scales_1.scale[i__ - 1];
	    x_move__ += d__1 * d__1;
	}
	x_move__ = sqrt(x_move__) / rms;
	if (s_cmp(output_1.coordtype, "INTERNAL", (ftnlen)9, (ftnlen)8) == 0) 
		{
	    x_move__ *= 57.29577951308232088;
	}

/*     compute the rms gradient per optimization parameter */

	g_rms__ = 0.;
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = g[i__ - 1] * scales_1.scale[i__ - 1];
	    g_rms__ += d__1 * d__1;
	}
	g_rms__ = sqrt(g_rms__) / rms;

/*     test for error due to line search problems */

	if (s_cmp(status, "BadIntpln", (ftnlen)9, (ftnlen)9) == 0 || s_cmp(
		status, "IntplnErr", (ftnlen)9, (ftnlen)9) == 0) {
	    ++nerr;
	    if (nerr >= maxerr) {
		done = TRUE_;
	    }
	} else {
	    nerr = 0;
	}

/*     test for too many total iterations */

	if (niter >= minima_1.maxiter) {
	    s_copy(status, "IterLimit", (ftnlen)9, (ftnlen)9);
	    done = TRUE_;
	}

/*     test the normal termination criteria */

	if (f <= minima_1.fctmin) {
	    s_copy(status, "SmallFct ", (ftnlen)9, (ftnlen)9);
	    done = TRUE_;
	}
	if (g_rms__ <= *grdmin) {
	    s_copy(status, "SmallGrad", (ftnlen)9, (ftnlen)9);
	    done = TRUE_;
	}

/*     print intermediate results for the current iteration */

	if (inform_1.iprint > 0) {
	    if (done || niter % inform_1.iprint == 0) {
		if (f < 1e7 && f > -1e6 && g_rms__ < 1e5 && f_move__ < 1e5) {
		    io___58.ciunit = iounit_1.iout;
		    s_wsfe(&io___58);
		    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&g_rms__, (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&f_move__, (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&x_move__, (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&angle, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
		    do_fio(&c__1, status, (ftnlen)9);
		    e_wsfe();
		} else {
		    io___59.ciunit = iounit_1.iout;
		    s_wsfe(&io___59);
		    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&g_rms__, (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&f_move__, (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&x_move__, (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&angle, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
		    do_fio(&c__1, status, (ftnlen)9);
		    e_wsfe();
		}
	    }
	}

/*     write intermediate results for the current iteration */

	if (inform_1.iwrite > 0) {
	    if (done || niter % inform_1.iwrite == 0) {
		(*optsave)(&niter, &f, &x[1]);
	    }
	}
    }

/*     set final value of the objective function */

    *minimum = f;
    if (inform_1.iprint > 0) {
	if (s_cmp(status, "SmallGrad", (ftnlen)9, (ftnlen)9) == 0 || s_cmp(
		status, "SmallFct ", (ftnlen)9, (ftnlen)9) == 0) {
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
    return 0;
} /* lbfgs_ */

#undef keyline_ref
#undef y_ref
#undef s_ref


