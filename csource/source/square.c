/* square.f -- translated by f2c (version 20050501).
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
    doublereal fctmin, hguess;
    integer maxiter, nextiter;
} minima_;

#define minima_1 minima_

/* Table of constant values */

static integer c__2 = 2;
static doublereal c_b7 = 10.;
static doublereal c_b9 = .33333333333333331;
static doublereal c_b10 = .66666666666666663;
static integer c__5 = 5;
static integer c__1 = 1;
static integer c__3 = 3;
static logical c_true = TRUE_;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine square  --  nonlinear least squares with bounds  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "square" is a nonlinear least squares routine derived from */
/*     the IMSL routine BCLSF and More's Minpack routine LMDER; the */
/*     Jacobian is estimated by finite differences and bounds can */
/*     be specified for the variables to be refined */

/*     arguments and variables : */

/*     n        number of variables to optimize */
/*     m        number of residual functions */
/*     xlo      vector of length n containing the lower bounds for */
/*                the variables */
/*     xhi      vector of length n containing the upper bounds for */
/*                the variables */
/*     xscale   vector of length n containing the diagonal scaling */
/*                matrix for the variables */
/*     xc       vector of length n containing the variable values */
/*                at the approximate solution */
/*     fc       vector of length m containing the residuals at the */
/*                approximate solution */
/*     fp       vector of length m containing the updated residual */
/*     xp       vector of length n containing the updated point */
/*     sc       vector of length n containing the last step taken */
/*     gc       vector of length n containing an estimate of */
/*                the gradient at the approximate solution */
/*     fjac     real m by n matrix containing an estimate of */
/*                the Jacobian at the approximate solution */
/*     mdim     leading dimension of fjac exactly as specified in */
/*                the dimension statement of the calling program */
/*     iactive  integer vector of length n indicating if xc(i) had */
/*                to be moved to an upper or lower bound */
/*     ipvt     vector of length n containing the permutation matrix */
/*                used in the QR factorization of the Jacobian at */
/*                the approximate solution */
/*     stpmax   real scalar containing the maximum allowed step size */
/*     delta    real scalar containing the trust region radius */

/*     required external routines : */

/*     rsdvalue   subroutine to evaluate residual function values */
/*     lsqwrite   subroutine to write out info about current status */


/* Subroutine */ int square_(integer *m, integer *n, doublereal *xlo, 
	doublereal *xhi, doublereal *xc, doublereal *fc, doublereal *gc, 
	doublereal *fjac, integer *mdim, doublereal *grdmin, S_fp rsdvalue, 
	S_fp lsqwrite)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 SQUARE  --  Too many Parameters,\002,"
	    "\002 Increase the Value of MAXLSQ\002)";
    static char fmt_20[] = "(/,\002 SQUARE  --  Too many Residuals,\002,\002"
	    " Increase the Value of MAXRSD\002)";
    static char fmt_40[] = "(/,\002 Levenberg-Marquardt Nonlinear Least Squa"
	    "res :\002)";
    static char fmt_50[] = "(/,\002 LS Iter    F Value      Total G     Acti"
	    "ve G\002,\002    N Active   F Calls\002,/)";
    static char fmt_60[] = "(i6,3f13.4,2i10)";
    static char fmt_70[] = "(i6,3d13.4,2i10)";
    static char fmt_100[] = "(i6,3f13.4,2i10)";
    static char fmt_110[] = "(i6,3d13.4,2i10)";
    static char fmt_130[] = "(/,\002 SQUARE  --  Normal Termination of Least"
	    " Squares\002)";
    static char fmt_140[] = "(/,\002 SQUARE  --  Maximum Number of Allowed I"
	    "terations\002)";
    static char fmt_150[] = "(/,\002 SQUARE  --  Relative Function Convergen"
	    "ce\002,//,\002 Both the scaled actual and predicted\002,\002 red"
	    "uctions in the function\002,/,\002 are less than or equal to the"
	    " relative\002,\002 convergence tolerance\002)";
    static char fmt_160[] = "(/,\002 SQUARE  --  Possible False Convergenc"
	    "e\002,//,\002 The iterates appear to be converging to\002,\002 a"
	    " noncritical point due\002,/,\002 to bad gradient information, d"
	    "iscontinuous\002,\002 function, or stopping\002,/,\002 tolerance"
	    "s being too tight\002)";
    static char fmt_170[] = "(/,\002 SQUARE  --  Five Consecutive Maximum"
	    "\002,\002 Length Steps\002,//,\002 Either the function is unboun"
	    "ded below,\002,\002 or has a finite\002,/,\002 asymptote in some"
	    " direction, or STEPMAX\002,\002 is too small\002)";

    /* System generated locals */
    integer fjac_dim1, fjac_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    double pow_di(doublereal *, integer *), pow_dd(doublereal *, doublereal *)
	    , sqrt(doublereal);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), do_fio(
	    integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k;
    static doublereal ga[50], sa[50], fp[100], sc[50], gs[50], xp[50], amu, 
	    eps, qtf[50], xsa[50], sum;
    static logical done;
    static doublereal temp;
    static integer next, ipvt[50];
    static doublereal work[50];
    static integer icode;
    static doublereal rdiag[50], delta, ftemp[100];
    static integer niter;
    static logical gauss;
    static doublereal rftol;
    static logical first;
    static doublereal xtemp;
    extern /* Subroutine */ int trust_(S_fp, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, logical *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, logical *, integer *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *);
    static integer ncalls, ndigit;
    static doublereal epsfcn, fcnorm, fpnorm, gcnorm, ganorm, stpmax, stpmin, 
	    xscale[50];
    static logical bigstp;
    static char record[120], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static doublereal stepsz;
    extern /* Subroutine */ int qrfact_(integer *, integer *, doublereal *, 
	    integer *, logical *, integer *, doublereal *, doublereal *), 
	    lmstep_(integer *, doublereal *, doublereal *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *, doublereal *, logical *);
    static integer iactive[50];
    static doublereal faketol;
    static integer nactive;
    extern doublereal precise_(integer *);
    static integer nbigstp;
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);
    static doublereal stpnorm;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___2 = { 0, 0, 0, fmt_20, 0 };
    static icilist io___20 = { 1, string, 1, 0, 120, 1 };
    static icilist io___21 = { 1, string, 1, 0, 120, 1 };
    static icilist io___22 = { 1, string, 1, 0, 120, 1 };
    static icilist io___23 = { 1, string, 1, 0, 120, 1 };
    static cilist io___36 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___58 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___59 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___64 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___65 = { 0, 0, 0, fmt_170, 0 };



#define fjac_ref(a_1,a_2) fjac[(a_2)*fjac_dim1 + a_1]
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

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  minima.i  --  general parameters for minimizations  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     fctmin    value below which function is deemed optimized */
/*     hguess    initial value for the H-matrix diagonal elements */
/*     maxiter   maximum number of iterations during optimization */
/*     nextiter  iteration number to use for the first iteration */




/*     check for too many variables or residuals */

    /* Parameter adjustments */
    --xlo;
    --xhi;
    --xc;
    --fc;
    --gc;
    fjac_dim1 = *mdim;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;

    /* Function Body */
    if (*n > 50) {
	io___1.ciunit = iounit_1.iout;
	s_wsfe(&io___1);
	e_wsfe();
	return 0;
    } else if (*m > 100) {
	io___2.ciunit = iounit_1.iout;
	s_wsfe(&io___2);
	e_wsfe();
	return 0;
    }

/*     initialize various counters and status code */

    niter = 0;
    ncalls = 0;
    nbigstp = 0;
    done = FALSE_;

/*     setup the default tolerances and parameter values */

    eps = precise_(&c__2);
    ndigit = 10;
    i__1 = -ndigit;
    if (eps < pow_di(&c_b7, &i__1)) {
	i__1 = -ndigit;
	eps = pow_di(&c_b7, &i__1);
    }
    if (minima_1.maxiter == 0) {
	minima_1.maxiter = 100;
    }
    if (inform_1.iprint < 0) {
	inform_1.iprint = 1;
    }
    if (inform_1.iwrite < 0) {
	inform_1.iwrite = 1;
    }
    if (minima_1.fctmin == 0.) {
	minima_1.fctmin = eps;
    }
    if (*grdmin == 0.) {
	*grdmin = pow_dd(&eps, &c_b9);
    }
    epsfcn = sqrt(eps);
    delta = 0.;
    stpmax = sqrt((doublereal) (*n)) * 1e3;
    stpmin = pow_dd(&eps, &c_b10);
    rftol = pow_dd(&eps, &c_b10);
    faketol = eps * 100.;

/*     search each line of the keyword file for options */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "FCTMIN ", (ftnlen)7, (ftnlen)7) == 0) {
	    i__2 = s_rsli(&io___20);
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&minima_1.fctmin, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L30;
	    }
	} else if (s_cmp(keyword, "MAXITER ", (ftnlen)8, (ftnlen)8) == 0) {
	    i__2 = s_rsli(&io___21);
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&minima_1.maxiter, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L30;
	    }
	} else if (s_cmp(keyword, "PRINTOUT ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___22);
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&inform_1.iprint, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L30;
	    }
	} else if (s_cmp(keyword, "WRITEOUT ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___23);
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&inform_1.iwrite, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L30;
	    }
	}
L30:
	;
    }

/*     check feasibility of variables and use bounds if needed */

    nactive = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (xc[j] < xlo[j]) {
	    xc[j] = xlo[j];
	    iactive[j - 1] = -1;
	} else if (xc[j] > xhi[j]) {
	    xc[j] = xhi[j];
	    iactive[j - 1] = 1;
	} else {
	    ++nactive;
	    iactive[j - 1] = 0;
	}
    }

/*     evaluate the function at the initial point */

    ++ncalls;
    (*rsdvalue)(m, n, &xc[1], &fc[1]);
    fcnorm = 0.;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = fc[i__];
	fcnorm += d__1 * d__1;
    }
    fcnorm *= .5;

/*     evaluate the Jacobian at the initial point by finite */
/*     differences; replace loop with user routine if desired */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	stepsz = epsfcn * (d__1 = xc[j], abs(d__1));
	if (stepsz < epsfcn) {
	    stepsz = epsfcn;
	}
	if (xc[j] < 0.) {
	    stepsz = -stepsz;
	}
	xtemp = xc[j];
	xc[j] = xtemp + stepsz;
	++ncalls;
	(*rsdvalue)(m, n, &xc[1], ftemp);
	xc[j] = xtemp;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    fjac_ref(i__, j) = (ftemp[i__ - 1] - fc[i__]) / stepsz;
	}
    }

/*     compute More's adaptive variable scale factors */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = 0.;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    d__1 = fjac_ref(i__, j);
	    temp += d__1 * d__1;
	}
	xscale[j - 1] = sqrt(temp);
	if (xscale[j - 1] == 0.) {
	    xscale[j - 1] = 1.;
	}
    }

/*     compute the total gradient vector for all variables */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	gc[j] = 0.;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    gc[j] += fjac_ref(i__, j) * fc[i__];
	}
    }

/*     compute the norm of the scaled total gradient */
/*     and the scaled gradient for active variables */

    gcnorm = 0.;
    ganorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__2 = (d__1 = xc[j], abs(d__1)), d__3 = 1. / xscale[j - 1];
	gs[j - 1] = gc[j] * max(d__2,d__3);
/* Computing 2nd power */
	d__1 = gs[j - 1];
	gcnorm += d__1 * d__1;
	if (iactive[j - 1] == 0) {
/* Computing 2nd power */
	    d__1 = gs[j - 1];
	    ganorm += d__1 * d__1;
	}
    }
    gcnorm = sqrt(gcnorm / *n);
    if (nactive != 0) {
	ganorm = sqrt(ganorm / nactive);
    }

/*     print out information about initial conditions */

    if (inform_1.iprint > 0) {
	io___36.ciunit = iounit_1.iout;
	s_wsfe(&io___36);
	e_wsfe();
	io___37.ciunit = iounit_1.iout;
	s_wsfe(&io___37);
	e_wsfe();
	if (max(fcnorm,gcnorm) < 1e7) {
	    io___38.ciunit = iounit_1.iout;
	    s_wsfe(&io___38);
	    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&fcnorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gcnorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ganorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nactive, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___39.ciunit = iounit_1.iout;
	    s_wsfe(&io___39);
	    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&fcnorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gcnorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ganorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nactive, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }

/*     write out the parameters, derivatives and residuals */

    if (inform_1.iwrite != 0) {
	(*lsqwrite)(&niter, &xc[1], gs, m, &fc[1]);
    }

/*     check stopping criteria at the initial point; test the */
/*     absolute function value and gradient norm for termination */

    if (fcnorm <= minima_1.fctmin) {
	return 0;
    }
    if (ganorm <= *grdmin) {
	return 0;
    }

/*     start of the main body of least squares iteration */

L80:
    ++niter;

/*     repack the Jacobian to include only active variables */

    if (nactive != *n) {
	k = 0;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (iactive[j - 1] != 0) {
		if (k == 0) {
		    k = j;
		}
	    } else {
		if (k != 0) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			fjac_ref(i__, k) = fjac_ref(i__, j);
		    }
		    ++k;
		}
	    }
	}
    }

/*     repack scale factors and gradient for active variables */

    k = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (iactive[j - 1] == 0) {
	    ++k;
	    xsa[k - 1] = xscale[j - 1];
	    ga[k - 1] = gc[j];
	}
    }

/*     compute the QR factorization of the Jacobian */

    qrfact_(m, &nactive, &fjac[fjac_offset], mdim, &c_true, ipvt, rdiag, work)
	    ;

/*     compute the vector Q(transpose) * residuals */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	qtf[i__ - 1] = fc[i__];
    }
    i__1 = nactive;
    for (j = 1; j <= i__1; ++j) {
	if (fjac_ref(j, j) != 0.) {
	    sum = 0.;
	    i__2 = *m;
	    for (i__ = j; i__ <= i__2; ++i__) {
		sum += fjac_ref(i__, j) * qtf[i__ - 1];
	    }
	    temp = -sum / fjac_ref(j, j);
	    i__2 = *m;
	    for (i__ = j; i__ <= i__2; ++i__) {
		qtf[i__ - 1] += fjac_ref(i__, j) * temp;
	    }
	}
	fjac_ref(j, j) = rdiag[j - 1];
    }

/*     compute the Levenberg-Marquardt step */

    icode = 6;
    first = TRUE_;
    while(icode >= 4) {
	lmstep_(&nactive, ga, &fjac[fjac_offset], mdim, ipvt, xsa, qtf, &
		stpmax, &delta, &amu, &first, sa, &gauss);

/*     unpack the step vector to include all variables */

	k = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (iactive[i__ - 1] != 0) {
		sc[i__ - 1] = 0.;
	    } else {
		++k;
		sc[i__ - 1] = sa[k - 1];
	    }
	}

/*     check new point and update the trust region */

	trust_((S_fp)rsdvalue, m, n, &xc[1], &fcnorm, &gc[1], &fjac[
		fjac_offset], mdim, ipvt, sc, sa, xscale, &gauss, &stpmax, &
		delta, &icode, xp, &fc[1], fp, &fpnorm, &bigstp, &ncalls, &
		xlo[1], &xhi[1], &nactive, &stpmin, &rftol, &faketol);
    }
    if (icode == 1) {
	done = TRUE_;
    }

/*     update to the new variables and residuals */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	xc[j] = xp[j - 1];
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	fc[i__] = fp[i__ - 1];
    }
    fcnorm = fpnorm;

/*     update the active vs inactive status of the variables; */
/*     in a true active set strategy, at most one constraint */
/*     is added to the active set per iteration */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (iactive[j - 1] == 0) {
	    if ((d__1 = xc[j] - xlo[j], abs(d__1)) <= eps) {
		--nactive;
		iactive[j - 1] = -1;
/*              goto 110 */
	    } else if ((d__1 = xc[j] - xhi[j], abs(d__1)) <= eps) {
		--nactive;
		iactive[j - 1] = 1;
/*              goto 110 */
	    }
	}
    }
/* L90: */

/*     evaluate the Jacobian at the new point using finite */
/*     differences; replace loop with user routine if desired */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__2 = (d__1 = xc[j], abs(d__1)), d__3 = 1. / xscale[j - 1];
	stepsz = epsfcn * max(d__2,d__3);
	if (xc[j] < 0.) {
	    stepsz = -stepsz;
	}
	xtemp = xc[j];
	xc[j] = xtemp + stepsz;
	++ncalls;
	(*rsdvalue)(m, n, &xc[1], ftemp);
	xc[j] = xtemp;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    fjac_ref(i__, j) = (ftemp[i__ - 1] - fc[i__]) / stepsz;
	}
    }

/*     compute More's adaptive variable scale factors */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = 0.;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    d__1 = fjac_ref(i__, j);
	    temp += d__1 * d__1;
	}
/* Computing MAX */
	d__1 = xscale[j - 1], d__2 = sqrt(temp);
	xscale[j - 1] = max(d__1,d__2);
    }

/*     compute the total gradient vector for all variables */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	gc[j] = 0.;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    gc[j] += fjac_ref(i__, j) * fc[i__];
	}
    }

/*     compute the norm of the scaled total gradient */
/*     and the scaled gradient for active variables */

    gcnorm = 0.;
    ganorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__2 = (d__1 = xc[j], abs(d__1)), d__3 = 1. / xscale[j - 1];
	gs[j - 1] = gc[j] * max(d__2,d__3);
/* Computing 2nd power */
	d__1 = gs[j - 1];
	gcnorm += d__1 * d__1;
	if (iactive[j - 1] == 0) {
/* Computing 2nd power */
	    d__1 = gs[j - 1];
	    ganorm += d__1 * d__1;
	}
    }
    gcnorm = sqrt(gcnorm / *n);
    if (nactive != 0) {
	ganorm = sqrt(ganorm / nactive);
    }

/*     print out information about current iteration */

    if (inform_1.iprint != 0 && niter % inform_1.iprint == 0) {
	if (max(fcnorm,gcnorm) < 1e7) {
	    io___58.ciunit = iounit_1.iout;
	    s_wsfe(&io___58);
	    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&fcnorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gcnorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ganorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nactive, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___59.ciunit = iounit_1.iout;
	    s_wsfe(&io___59);
	    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&fcnorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gcnorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ganorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nactive, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }

/*     check stopping criteria at the new point; test the absolute */
/*     function value, the gradient norm and step for termination */

    if (fcnorm <= minima_1.fctmin) {
	done = TRUE_;
    }
    if (ganorm <= *grdmin) {
	done = TRUE_;
    }
    stpnorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__2 = (d__1 = xc[j], abs(d__1)), d__3 = 1. / xscale[j - 1];
	temp = max(d__2,d__3);
/* Computing 2nd power */
	d__1 = sc[j - 1] / temp;
	stpnorm += d__1 * d__1;
    }
    stpnorm = sqrt(stpnorm / *n);
    if (stpnorm <= stpmin) {
	done = TRUE_;
    }

/*     check for inactive variables that can be made active; */
/*     in a true active set strategy, variables are released */
/*     one at a time at a minimum of the current active set */

/*     if (done) then */
    if (nactive != *n) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (iactive[j - 1] == -1 && gc[j] < 0.) {
		++nactive;
		iactive[j - 1] = 0;
		done = FALSE_;
/*              goto 140 */
	    } else if (iactive[j - 1] == 1 && gc[j] > 0.) {
		++nactive;
		iactive[j - 1] = 0;
		done = FALSE_;
/*              goto 140 */
	    }
	}
/* L120: */
    }
/*     end if */

/*     if still done, then normal termination has been achieved */

    if (done) {
	io___61.ciunit = iounit_1.iout;
	s_wsfe(&io___61);
	e_wsfe();

/*     check the limit on the number of iterations */

    } else if (niter >= minima_1.maxiter) {
	done = TRUE_;
	io___62.ciunit = iounit_1.iout;
	s_wsfe(&io___62);
	e_wsfe();

/*     check for termination due to relative function convergence */

    } else if (icode == 2) {
	done = TRUE_;
	io___63.ciunit = iounit_1.iout;
	s_wsfe(&io___63);
	e_wsfe();

/*     check for termination due to false convergence */

    } else if (icode == 3) {
	done = TRUE_;
	io___64.ciunit = iounit_1.iout;
	s_wsfe(&io___64);
	e_wsfe();

/*     check for several consecutive maximum steps taken */

    } else if (bigstp) {
	++nbigstp;
	if (nbigstp == 5) {
	    done = TRUE_;
	    io___65.ciunit = iounit_1.iout;
	    s_wsfe(&io___65);
	    e_wsfe();
	}

/*     no reason to quit, so prepare to take another step */

    } else {
	nbigstp = 0;
    }

/*     write out the parameters, derivatives and residuals */

    if (inform_1.iwrite != 0 && niter % inform_1.iwrite == 0) {
	if (! done) {
	    (*lsqwrite)(&niter, &xc[1], gs, m, &fc[1]);
	}
    }

/*     continue with the next iteration if not finished */

    if (! done) {
	goto L80;
    }
    return 0;
} /* square_ */

#undef keyline_ref
#undef fjac_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine lmstep  --  finds Levenberg-Marquardt step  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "lmstep" computes the Levenberg-Marquardt step during a */
/*     nonlinear least squares calculation; this version is based */
/*     upon ideas from the Minpack routine LMPAR together with */
/*     with the internal doubling strategy of Dennis and Schnabel */

/*     arguments and variables : */

/*     n        dimension of the problem */
/*     ga       real vector of length n containing the gradient */
/*                of the residual vector */
/*     a        real n by n array which on input contains in the full */
/*                upper triangle the upper triangle of the matrix r */
/*                resulting from the QR factorization of the Jacobian; */
/*                on output the full upper triangle is unaltered, and */
/*                the strict lower triangle contains the strict lower */
/*                triangle of the lower triangular matrix l which is */
/*                equal to the Cholesky factor of (j**t)*j + amu*xscale */
/*     lda      leading dimension of a exactly as specified in the */
/*                dimension statement of the calling program */
/*     ipvt     integer array of length n containing the pivoting */
/*                infomation from the QR factorization routine */
/*     xscale   real vector of length n containing the diagonal */
/*                scaling matrix for the variables */
/*     qtf      real vector containing the first n elements of */
/*                Q(transpose) * (the scaled residual vector) */
/*     amu      scalar containing an initial estimate of the */
/*                Levenberg-Marquardt parameter on input; on */
/*                output, amu contains the final estimate of */
/*                the Levenberg-Marquardt parameter */
/*     first    logical variable set true only if this is the first */
/*                call to this routine in this iteration */
/*     sa       real vector of length n containing the */
/*                Levenberg-Marquardt step */
/*     gnstep   real vector of length n containing the */
/*                Gauss-Newton step */
/*     gauss    logical variable which is true if the Gauss-Newton */
/*                step is acceptable, and is false otherwise */
/*     diag     vector of length n containing the diagonal elements */
/*                of the Cholesky factor of (j**t)*j + amu*xscale */


/* Subroutine */ int lmstep_(integer *n, doublereal *ga, doublereal *a, 
	integer *lda, integer *ipvt, doublereal *xscale, doublereal *qtf, 
	doublereal *stpmax, doublereal *delta, doublereal *amu, logical *
	first, doublereal *sa, logical *gauss)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal phi, sum, diag[50], beta, high;
    static logical done;
    static doublereal phip, alow, temp, tiny, work1[50], work2[50], alpha, 
	    amuhi, small, phipi;
    static integer nsing;
    static doublereal deltap, gnleng, gnstep[50], amulow, sgnorm, stplen;
    extern doublereal precise_(integer *);
    extern /* Subroutine */ int qrsolve_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *);


#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]



/*     set smallest floating point magnitude and spacing */

    /* Parameter adjustments */
    --ga;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --xscale;
    --qtf;
    --sa;

    /* Function Body */
    tiny = precise_(&c__1);
    small = precise_(&c__2);

/*     if initial trust region is not provided by the user, */
/*     compute and use the length of the Cauchy step given */
/*     by beta = norm2(r*trans(p)*d**(-2)*g)**2 */

    if (*delta == 0.) {
	*amu = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work1[i__ - 1] = ga[i__] / xscale[i__];
	}
	alpha = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = work1[i__ - 1];
	    alpha += d__1 * d__1;
	}
	beta = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    temp = 0.;
	    i__2 = *n;
	    for (j = i__; j <= i__2; ++j) {
		k = ipvt[j];
/* Computing 2nd power */
		d__1 = xscale[k];
		temp += a_ref(i__, j) * ga[k] / (d__1 * d__1);
	    }
/* Computing 2nd power */
	    d__1 = temp;
	    beta += d__1 * d__1;
	}
	if (beta <= tiny) {
	    *delta = alpha * sqrt(alpha);
	} else {
	    *delta = alpha * sqrt(alpha) / beta;
	}
	*delta = min(*delta,*stpmax);
    }

/*     the following is done only on the first time through */
/*     this iteration: (1) compute the Gauss-Newton step; */
/*     if the Jacobian is rank-deficient, obtain a least */
/*     squares solution, (2) compute the length of the scaled */
/*     Gauss-Newton step, (3) compute the norm of the scaled */
/*     gradient used in computing an upper bound for "amu" */

    if (*first) {
	nsing = *n;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (a_ref(j, j) == 0. && nsing == *n) {
		nsing = j - 1;
	    }
	    if (nsing < *n) {
		work1[j - 1] = 0.;
	    }
	}
	work1[nsing - 1] = qtf[nsing] / a_ref(nsing, nsing);
	for (j = nsing - 1; j >= 1; --j) {
	    sum = 0.;
	    i__1 = nsing;
	    for (i__ = j + 1; i__ <= i__1; ++i__) {
		sum += a_ref(j, i__) * work1[i__ - 1];
	    }
	    work1[j - 1] = (qtf[j] - sum) / a_ref(j, j);
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    gnstep[ipvt[j] - 1] = -work1[j - 1];
	}

/*     find the length of scaled Gauss-Newton step */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    work1[j - 1] = xscale[j] * gnstep[j - 1];
	}
	gnleng = 0.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	    d__1 = work1[j - 1];
	    gnleng += d__1 * d__1;
	}
	gnleng = sqrt(gnleng);

/*     find the length of the scaled gradient */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    work1[j - 1] = ga[j] / xscale[j];
	}
	sgnorm = 0.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	    d__1 = work1[j - 1];
	    sgnorm += d__1 * d__1;
	}
	sgnorm = sqrt(sgnorm);
    }

/*     set the bounds on the computed step */

    high = 1.5;
    alow = .75;

/*     check to see if the Gauss-Newton step is acceptable */

    if (gnleng <= high * *delta) {
	*gauss = TRUE_;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sa[j] = gnstep[j - 1];
	}
	*amu = 0.;
	*delta = min(*delta,gnleng);

/*     the Gauss-Newton step is rejected, find a nontrivial step; */
/*     first compute a starting value of "amu" if previous step */
/*     was not a Gauss-Newton step */

    } else {
	*gauss = FALSE_;
	if (*amu > 0.) {
	    *amu -= (phi + deltap) / *delta * ((deltap - *delta + phi) / phip)
		    ;
	}
	phi = gnleng - *delta;

/*     if the Jacobian is not rank deficient, the Newton step */
/*     provides a lower bound for "amu"; else set bound to zero */

	if (nsing == *n) {
	    if (*first) {
		*first = FALSE_;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    k = ipvt[j];
/* Computing 2nd power */
		    d__1 = xscale[k];
		    work1[j - 1] = gnstep[k - 1] * (d__1 * d__1);
		}

/*     obtain trans(r**-1)*(trans(p)*s) by solving the */
/*     system of equations trans(r)*work1 = work1 */

		work1[*n - 1] /= a_ref(*n, *n);
		for (j = *n - 1; j >= 1; --j) {
		    sum = 0.;
		    i__1 = *n;
		    for (i__ = j + 1; i__ <= i__1; ++i__) {
			sum += a_ref(j, i__) * work1[i__ - 1];
		    }
		    work1[j - 1] = (work1[j - 1] - sum) / a_ref(j, j);
		}
		phipi = 0.;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
		    d__1 = work1[j - 1];
		    phipi -= d__1 * d__1;
		}
		phipi /= gnleng;
	    }
	    amulow = -phi / phipi;
	} else {
	    *first = FALSE_;
	    amulow = 0.;
	}
	amuhi = sgnorm / *delta;

/*     iterate until a satisfactory "amu" is generated */

	done = FALSE_;
	while(! done) {
	    if (*amu < amulow || *amu > amuhi) {
/* Computing MAX */
		d__1 = sqrt(amulow * amuhi), d__2 = amuhi * .001;
		*amu = max(d__1,d__2);
	    }
	    temp = sqrt(*amu);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		work1[j - 1] = temp * xscale[j];
	    }

/*     solve the damped least squares system for the value of the */
/*     Levenberg-Marquardt step using More's Minpack technique */

	    qrsolve_(n, &a[a_offset], lda, &ipvt[1], work1, &qtf[1], &sa[1], 
		    diag, work2);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		sa[j] = -sa[j];
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		work2[j - 1] = xscale[j] * sa[j];
	    }
	    stplen = 0.;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
		d__1 = work2[j - 1];
		stplen += d__1 * d__1;
	    }
	    stplen = sqrt(stplen);
	    phi = stplen - *delta;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		k = ipvt[j];
		work1[j - 1] = xscale[k] * work2[k - 1];
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if ((d__1 = diag[j - 1], abs(d__1)) >= tiny) {
		    work1[j - 1] /= diag[j - 1];
		}
		if (j < *n) {
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			work1[i__ - 1] -= work1[j - 1] * a_ref(i__, j);
		    }
		}
	    }
	    phip = 0.;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
		d__1 = work1[j - 1];
		phip -= d__1 * d__1;
	    }
	    phip /= stplen;

/*     check to see if the step is acceptable; if not, */
/*     update amulow, amuhi and amu for next iteration */

	    if (stplen >= alow * *delta && stplen <= high * *delta || amuhi - 
		    amulow <= small) {
		done = TRUE_;
	    } else {
/* Computing MAX */
		d__1 = amulow, d__2 = *amu - phi / phip;
		amulow = max(d__1,d__2);
		if (phi < 0.) {
		    amuhi = *amu;
		}
		*amu -= stplen / *delta * (phi / phip);
	    }
	}
    }
    deltap = *delta;
    return 0;
} /* lmstep_ */

#undef a_ref




/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine qrfact  --  QR factorization of a matrix  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "qrfact" performs Householder transformations with column */
/*     pivoting (optional) to compute a QR factorization of the */
/*     m by n matrix a; the routine determines an orthogonal */
/*     matrix q, a permutation matrix p, and an upper trapezoidal */
/*     matrix r with diagonal elements of nonincreasing magnitude, */
/*     such that a*p = q*r; the Householder transformation for */
/*     column k, k = 1,2,...,min(m,n), is of the form */

/*               i - (1/u(k))*u*u(transpose) */

/*     where u has zeros in the first k-1 positions */

/*     arguments and variables : */

/*     m        positive integer input variable set to */
/*                the number of rows of a */
/*     n        positive integer input variable set to */
/*                the number of columns of a */
/*     a        m by n array; on input a contains the matrix for */
/*                which the QR factorization is to be computed; on */
/*                output the strict upper trapezoidal part contains */
/*                the strict upper trapezoidal part of r, the lower */
/*                trapezoidal part contains a factored form of q */
/*                (non-trivial elements of u vectors described above) */
/*     lda      positive integer input variable not less than m */
/*                which specifies the leading dimension of array a */
/*     pivot    logical input variable; if pivot is set true, */
/*                then column pivoting is enforced; if pivot is */
/*                set false, then no column pivoting is done */
/*     ipvt     integer output array which defines the permutation */
/*                matrix p such that a*p = q*r; column j of p is */
/*                column ipvt(j) of the identity matrix */
/*     rdiag    an output array of length n which contains */
/*                the diagonal elements of r */


/* Subroutine */ int qrfact_(integer *m, integer *n, doublereal *a, integer *
	lda, logical *pivot, integer *ipvt, doublereal *rdiag, doublereal *
	work)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, jmax;
    static doublereal temp;
    static integer minmn, itemp;
    static doublereal aknorm;


#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]



/*     find the initial column norms and initialize some arrays */

    /* Parameter adjustments */
    --work;
    --rdiag;
    --ipvt;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = 0.;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    d__1 = a_ref(i__, j);
	    temp += d__1 * d__1;
	}
	rdiag[j] = sqrt(temp);
	work[j] = rdiag[j];
	if (*pivot) {
	    ipvt[j] = j;
	}
    }

/*     reduce the matrix with Householder transformations */

    minmn = min(*m,*n);
    i__1 = minmn;
    for (k = 1; k <= i__1; ++k) {

/*     bring the column of largest norm into the pivot position */

	if (*pivot) {
	    jmax = k;
	    i__2 = *n;
	    for (j = k; j <= i__2; ++j) {
		if (rdiag[j] > rdiag[jmax]) {
		    jmax = j;
		}
	    }
	    if (jmax != k) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = a_ref(i__, k);
		    a_ref(i__, k) = a_ref(i__, jmax);
		    a_ref(i__, jmax) = temp;
		}
		rdiag[jmax] = rdiag[k];
		work[jmax] = work[k];
		itemp = ipvt[k];
		ipvt[k] = ipvt[jmax];
		ipvt[jmax] = itemp;
	    }
	}

/*     compute the Householder transformation to reduce the */
/*     k-th column of a to a multiple of the k-th unit vector */

	aknorm = 0.;
	i__2 = *m;
	for (i__ = k; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    d__1 = a_ref(i__, k);
	    aknorm += d__1 * d__1;
	}
	aknorm = sqrt(aknorm);
	if (aknorm != 0.) {
	    if (a_ref(k, k) < 0.) {
		aknorm = -aknorm;
	    }
	    i__2 = *m;
	    for (i__ = k; i__ <= i__2; ++i__) {
		a_ref(i__, k) = a_ref(i__, k) / aknorm;
	    }
	    a_ref(k, k) = a_ref(k, k) + 1.;

/*     apply the transformation to the remaining columns */
/*     and update the column norms */

	    if (*n >= k + 1) {
		i__2 = *n;
		for (j = k + 1; j <= i__2; ++j) {
		    temp = 0.;
		    i__3 = *m;
		    for (i__ = k; i__ <= i__3; ++i__) {
			temp += a_ref(i__, k) * a_ref(i__, j);
		    }
		    temp /= a_ref(k, k);
		    i__3 = *m;
		    for (i__ = k; i__ <= i__3; ++i__) {
			a_ref(i__, j) = a_ref(i__, j) - temp * a_ref(i__, k);
		    }
		    if (*pivot && rdiag[j] != 0.) {
			temp = a_ref(k, j) / rdiag[j];
			if (abs(temp) < 1.) {
/* Computing 2nd power */
			    d__1 = temp;
			    rdiag[j] *= sqrt(1. - d__1 * d__1);
			} else {
			    temp = 0.;
			    i__3 = *m;
			    for (i__ = k + 1; i__ <= i__3; ++i__) {
/* Computing 2nd power */
				d__1 = a_ref(i__, j);
				temp += d__1 * d__1;
			    }
			    rdiag[j] = sqrt(temp);
			    work[j] = rdiag[j];
			}
		    }
		}
	    }
	}
	rdiag[k] = -aknorm;
    }
    return 0;
} /* qrfact_ */

#undef a_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine qrsolve  --  get least squares via QR factors  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "qrsolve" solves a*x=b and d*x=0 in the least squares sense; */
/*     normally used in combination with routine "qrfact" to solve */
/*     least squares problems */

/*     arguments and variables : */

/*     n        number of rows and columns in the matrix r */
/*     r        an n by n array containing the upper triangular */
/*                matrix r; on output the full triangle is unaltered, */
/*                and the strict lower triangle contains the transpose */
/*                of the strict upper triangular matrix s */
/*     ldr      leading dimension of r exactly as specified in */
/*                the dimension statement of the calling program */
/*     ipvt     vector of length n which defines the permutation */
/*                matrix p such that a*p = q*r; column j of p is */
/*                column ipvt(j) of the identity matrix */
/*     diag     vector of length n containing the diagonal elements */
/*                of the matrix d */
/*     qtb      vector of length n containing the first n elements */
/*                of the vector q(transpose)*b */
/*     x        vector of length n containing the least squares */
/*                solution of the systems a*x = b, d*x = 0 */
/*     sdiag    vector of length n containing the diagonal elements */
/*                of the upper triangular matrix s */


/* Subroutine */ int qrsolve_(integer *n, doublereal *r__, integer *ldr, 
	integer *ipvt, doublereal *diag, doublereal *qtb, doublereal *x, 
	doublereal *sdiag, doublereal *work)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal cotangent;
    static integer jj;
    static doublereal sine, temp;
    static integer nsing;
    static doublereal qtbpj, cosine, tangent;


#define r___ref(a_1,a_2) r__[(a_2)*r_dim1 + a_1]



/*     copy r and (q transpose)*b to preserve input and */
/*     initialize s; in particular, save the diagonal */
/*     elements of r in x */

    /* Parameter adjustments */
    --work;
    --sdiag;
    --x;
    --qtb;
    --diag;
    --ipvt;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;

    /* Function Body */
    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (k = j + 1; k <= i__2; ++k) {
	    r___ref(k, j) = r___ref(j, k);
	}
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] = r___ref(j, j);
	work[j] = qtb[j];
    }

/*     eliminate the diagonal matrix d using a Givens rotation */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {

/*     prepare the row of d to be eliminated, locating */
/*     the diagonal element using p from the QR factorization */

	jj = ipvt[j];
	if (diag[jj] != 0.) {
	    i__2 = *n;
	    for (k = j; k <= i__2; ++k) {
		sdiag[k] = 0.;
	    }
	    sdiag[j] = diag[jj];

/*     the transformations to eliminate the row of d modify */
/*     only a single element of (q transpose)*b beyond the */
/*     first n, which is initially zero */

	    qtbpj = 0.;
	    i__2 = *n;
	    for (k = j; k <= i__2; ++k) {

/*     determine a Givens rotation which eliminates the */
/*     appropriate element in the current row of d */

		if (sdiag[k] != 0.) {
		    if ((d__1 = r___ref(k, k), abs(d__1)) < (d__2 = sdiag[k], 
			    abs(d__2))) {
			cotangent = r___ref(k, k) / sdiag[k];
/* Computing 2nd power */
			d__1 = cotangent;
			sine = .5 / sqrt(d__1 * d__1 * .25 + .25);
			cosine = sine * cotangent;
		    } else {
			tangent = sdiag[k] / r___ref(k, k);
/* Computing 2nd power */
			d__1 = tangent;
			cosine = .5 / sqrt(d__1 * d__1 * .25 + .25);
			sine = cosine * tangent;
		    }

/*     compute the modified diagonal element of r */
/*     and the modified element of ((q transpose)*b,0) */

		    r___ref(k, k) = cosine * r___ref(k, k) + sine * sdiag[k];
		    temp = cosine * work[k] + sine * qtbpj;
		    qtbpj = -sine * work[k] + cosine * qtbpj;
		    work[k] = temp;

/*     accumulate the tranformation in the row of s */

		    if (*n >= k + 1) {
			i__3 = *n;
			for (i__ = k + 1; i__ <= i__3; ++i__) {
			    temp = cosine * r___ref(i__, k) + sine * sdiag[
				    i__];
			    sdiag[i__] = -sine * r___ref(i__, k) + cosine * 
				    sdiag[i__];
			    r___ref(i__, k) = temp;
			}
		    }
		}
	    }
	}

/*     store the diagonal element of s and restore */
/*     the corresponding diagonal element of r */

	sdiag[j] = r___ref(j, j);
	r___ref(j, j) = x[j];
    }

/*     solve the triangular system for z; if the system */
/*     is singular, then obtain a least squares solution */

    nsing = *n;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (sdiag[j] == 0. && nsing == *n) {
	    nsing = j - 1;
	}
	if (nsing < *n) {
	    work[j] = 0.;
	}
    }
    if (nsing >= 1) {
	i__1 = nsing;
	for (k = 1; k <= i__1; ++k) {
	    j = nsing - k + 1;
	    temp = 0.;
	    if (nsing >= j + 1) {
		i__2 = nsing;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    temp += r___ref(i__, j) * work[i__];
		}
	    }
	    work[j] = (work[j] - temp) / sdiag[j];
	}
    }

/*     permute the components of z back to components of x */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	k = ipvt[j];
	x[k] = work[j];
    }
    return 0;
} /* qrsolve_ */

#undef r___ref




/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine trust  --  update the model trust region  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "trust" updates the model trust region for a nonlinear */
/*     least squares calculation; this version is based on the */
/*     ideas found in NL2SOL and in Dennis and Schnabel's book */

/*     arguments and variables : */

/*     m        number of functions */
/*     n        number of variables */
/*     xc       vector of length n containing the current iterate */
/*     fcnorm   real scalar containing the norm of f(xc) */
/*     gc       vector of length n containing the gradient at xc */
/*     a        real m by n matrix containing the upper triangular */
/*                matrix r from the QR factorization of the current */
/*                Jacobian in the upper triangle */
/*     lda      leading dimension of A exactly as specified in */
/*                the dimension statement of the calling program */
/*     ipvt     integer vector of length n containing the permutation */
/*                matrix from QR factorization of the Jacobian */
/*     sc       vector of length n containing the Newton step */
/*     sa       vector of length n containing current step */
/*     xscale   vector of length n containing the diagonal */
/*                scaling matrix for x */
/*     gauss    logical variable equal is true when the Gauss-Newton */
/*                step is taken */
/*     stpmax   maximum allowable step size */
/*     delta    trust region radius with value retained between calls */
/*     icode    return code set upon exit */
/*                0  means xp accepted as next iterate, delta */
/*                     is trust region for next iteration */
/*                1  means the algorithm was unable to find a */
/*                     satisfactory xp sufficiently distinct from xc */
/*                2  means both the scaled actual and predicted */
/*                     function reductions are smaller than rftol */
/*                3  means that false convergence is detected */
/*                4  means fpnorm is too large, current iteration is */
/*                     continued with a new, reduced trust region */
/*                5  means fpnorm is sufficiently small, but the */
/*                     chance of taking a longer successful step */
/*                     seems good that the current iteration is to */
/*                     be continued with a new, doubled trust region */
/*     xpprev   vector of length n containing the value of xp */
/*                at the previous call within this iteration */
/*     fpprev   vector of length m containing f(xpprev) */
/*     xp       vector of length n containing the new iterate */
/*     fp       vector of length m containing the functions at xp */
/*     fpnorm   scalar containing the norm of f(xp) */
/*     bigstp   logical variable; true, if maximum step length was taken, */
/*                false  otherwise */
/*     ncalls   number of function evaluations used */
/*     xlo      vector of length n containing the lower bounds */
/*     xhi      vector of length n containing the upper bounds */
/*     nactive  number of columns in the active Jacobian */

/*     required external routines : */

/*     rsdvalue   subroutine to evaluate residual function values */


/* Subroutine */ int trust_(S_fp rsdvalue, integer *m, integer *n, doublereal 
	*xc, doublereal *fcnorm, doublereal *gc, doublereal *a, integer *lda, 
	integer *ipvt, doublereal *sc, doublereal *sa, doublereal *xscale, 
	logical *gauss, doublereal *stpmax, doublereal *delta, integer *icode,
	 doublereal *xp, doublereal *fc, doublereal *fp, doublereal *fpnorm, 
	logical *bigstp, integer *ncalls, doublereal *xlo, doublereal *xhi, 
	integer *nactive, doublereal *stpmin, doublereal *rftol, doublereal *
	faketol)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static logical feas;
    static doublereal temp, tiny, alpha;
    static logical ltemp;
    static doublereal slope, reduce, rellen, fpnrmp, fpprev[100], stplen, 
	    xpprev[50], predict;
    extern doublereal precise_(integer *);


#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]



/*     set value of alpha, logicals and step length */

    /* Parameter adjustments */
    --xc;
    --gc;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --sc;
    --sa;
    --xscale;
    --xp;
    --fc;
    --fp;
    --xlo;
    --xhi;

    /* Function Body */
    alpha = 1e-4;
    *bigstp = FALSE_;
    feas = TRUE_;
    stplen = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = xscale[i__] * sc[i__];
	stplen += d__1 * d__1;
    }
    stplen = sqrt(stplen);

/*     compute new trial point and new function values */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xp[i__] = xc[i__] + sc[i__];
	if (xp[i__] > xhi[i__]) {
	    sc[i__] = xhi[i__] - xc[i__];
	    xp[i__] = xhi[i__];
	    feas = FALSE_;
	} else if (xp[i__] < xlo[i__]) {
	    sc[i__] = xlo[i__] - xc[i__];
	    xp[i__] = xlo[i__];
	    feas = FALSE_;
	}
    }
    ++(*ncalls);
    (*rsdvalue)(m, n, &xp[1], &fp[1]);
    *fpnorm = 0.;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = fp[i__];
	*fpnorm += d__1 * d__1;
    }
    *fpnorm *= .5;
    reduce = *fpnorm - *fcnorm;
    slope = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	slope += gc[i__] * sc[i__];
    }
    if (*icode != 5) {
	fpnrmp = 0.;
    }

/*     internal doubling no good; reset to previous and quit */

    if (*icode == 5 && (*fpnorm >= fpnrmp || reduce > alpha * slope)) {
	*icode = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xp[i__] = xpprev[i__ - 1];
	}
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    fp[i__] = fpprev[i__ - 1];
	}
	*fpnorm = fpnrmp;
	*delta *= .5;

/*     fpnorm is too large; the step is unacceptable */

    } else if (reduce >= alpha * slope) {
	rellen = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__3 = (d__2 = xp[i__], abs(d__2)), d__4 = 1. / xscale[i__];
	    temp = (d__1 = sc[i__], abs(d__1)) / max(d__3,d__4);
	    rellen = max(rellen,temp);
	}

/*     magnitude of (xp-xc) is too small, end the global step */

	if (rellen < *stpmin) {
	    *icode = 1;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		xp[i__] = xc[i__];
	    }
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		fp[i__] = fc[i__];
	    }

/*     quadratic interpolation step; reduce delta and continue */

	} else {
	    *icode = 4;
	    tiny = precise_(&c__1);
	    if ((d__1 = reduce - slope, abs(d__1)) > tiny) {
		temp = -slope * stplen / ((reduce - slope) * 2.);
	    } else {
		temp = -slope * stplen / 2.;
	    }
	    if (temp < *delta * .1) {
		*delta *= .1;
	    } else if (temp > *delta * .5) {
		*delta *= .5;
	    } else {
		*delta = temp;
	    }
	}

/*     fpnorm is sufficiently small; the step is acceptable compute */
/*     the predicted reduction as predict = g(T)*s + (1/2)*s(T)*h*s */
/*     with h = p * r**t * r * p**t */

    } else {
	predict = slope;
	i__1 = *nactive;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = ipvt[i__];
	    temp = 0.;
	    i__2 = *nactive;
	    for (j = i__; j <= i__2; ++j) {
		temp += sa[k] * a_ref(i__, j);
	    }
/* Computing 2nd power */
	    d__1 = temp;
	    predict += d__1 * d__1 * .5;
	}
	ltemp = (d__1 = predict - reduce, abs(d__1)) <= abs(reduce) * .1;

/*     if reduce and predict agree to within relative error of 0.1 */
/*     or if negative curvature is indicated, and a longer step is */
/*     possible and delta has not been decreased this iteration, */
/*     then double trust region and continue global step */

	if (*icode != 4 && (ltemp || reduce <= slope) && feas && ! (*gauss) &&
		 *delta <= *stpmax * .99) {
	    *icode = 5;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		xpprev[i__ - 1] = xp[i__];
	    }
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		fpprev[i__ - 1] = fp[i__];
	    }
	    fpnrmp = *fpnorm;
/* Computing MIN */
	    d__1 = *delta * 2.;
	    *delta = min(d__1,*stpmax);

/*     accept the point; choose new trust region for next iteration */

	} else {
	    *icode = 0;
	    if (stplen > *stpmax * .99) {
		*bigstp = TRUE_;
	    }
	    if (reduce >= predict * .1) {
		*delta *= .5;
	    } else if (reduce <= predict * .75) {
/* Computing MIN */
		d__1 = *delta * 2.;
		*delta = min(d__1,*stpmax);
	    }
	}

/*     check relative function convergence and false convergence */

	if (reduce <= predict * 2.) {
	    if (abs(reduce) <= *rftol * abs(*fcnorm) && abs(predict) <= *
		    rftol * abs(*fcnorm)) {
		*icode = 2;
	    }
	} else {
	    rellen = 0.;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
		d__3 = (d__2 = xp[i__], abs(d__2)), d__4 = 1. / xscale[i__];
		temp = (d__1 = sc[i__], abs(d__1)) / max(d__3,d__4);
		rellen = max(rellen,temp);
	    }
	    if (rellen < *faketol) {
		*icode = 3;
	    }
	}
    }
    return 0;
} /* trust_ */

#undef a_ref


