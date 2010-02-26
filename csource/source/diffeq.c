/* diffeq.f -- translated by f2c (version 20050501).
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
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

struct {
    doublereal m2[25000], deform, difft, diffv, diffc;
    logical use_smooth__, use_dem__, use_gda__, use_tophat__, use_stophat__;
} warp_;

#define warp_1 warp_

/* Table of constant values */

static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine diffeq  --  differential equation integration  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "diffeq" performs the numerical integration of an ordinary */
/*     differential equation using an adaptive stepsize method to */
/*     solve the corresponding coupled first-order equations of the */
/*     general form dyi/dx = f(x,y1,...,yn) for yi = y1,...,yn */

/*     variables and parameters : */

/*     nvar      number of coupled first-order differential equations */
/*     y         contains the values of the dependent variables */
/*     x1        value of the beginning integration limit */
/*     x2        value of the ending integration limit */
/*     eps       relative accuracy required of the integration steps */
/*     h1        initial guess for the first integration stepsize */
/*     hmin      minimum allowed integration stepsize */
/*     nok       number of initially successful integration steps */
/*     nbad      number of integration steps that required retry */

/*     required external routines : */

/*     gvalue    subroutine to find the right-hand side of the */
/*                  first-order differential equations */


/* Subroutine */ int diffeq_(integer *nvar, doublereal *y, doublereal *x1, 
	doublereal *x2, doublereal *eps, doublereal *h1, doublereal *hmin, 
	integer *nok, integer *nbad, S_fp gvalue)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 DIFFEQ  --  Normal Termination\002,\002 "
	    "at Integration Limit\002)";
    static char fmt_20[] = "(/,\002 DIFFEQ  --  Incomplete Integration\002"
	    ",\002 due to SmallStep\002)";
    static char fmt_30[] = "(/,\002 DIFFEQ  --  Incomplete Integration\002"
	    ",\002 due to IterLimit\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static doublereal h__;
    static integer i__;
    static doublereal x;
    static logical terminate;
    static doublereal hdid, dydx[100000], yscal[100000], hnext;
    static integer nstep;
    extern /* Subroutine */ int bsstep_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, S_fp);
    static char status[7];
    extern /* Subroutine */ int gdastat_(integer *, doublereal *, doublereal *
	    , char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_30, 0 };




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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     initialize starting limit, step size and status counters */

    /* Parameter adjustments */
    --y;

    /* Function Body */
    terminate = FALSE_;
    x = *x1;
    d__1 = *x2 - *x1;
    h__ = d_sign(h1, &d__1);
    nstep = 0;
    *nok = 0;
    *nbad = 0;

/*     perform a series of individual integration steps */

    while(! terminate) {
	(*gvalue)(&x, &y[1], dydx);
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    yscal[i__ - 1] = (d__1 = y[i__], abs(d__1)) + (d__2 = h__ * dydx[
		    i__ - 1], abs(d__2)) + 1e-30;
	}

/*     set the final step to stop at the integration limit */

	if ((x + h__ - *x2) * (x + h__ - *x1) > 0.) {
	    h__ = *x2 - x;
	}

/*     take a Bulirsch-Stoer integration step */

	bsstep_(nvar, &x, dydx, &y[1], &h__, eps, yscal, &hdid, &hnext, (S_fp)
		gvalue);

/*     mark the current step as either good or bad */

	if (hdid == h__) {
	    ++(*nok);
	    s_copy(status, "Success", (ftnlen)7, (ftnlen)7);
	} else {
	    ++(*nbad);
	    s_copy(status, " Retry ", (ftnlen)7, (ftnlen)7);
	}

/*     update stepsize and get information about the current step */

	h__ = hnext;
	++nstep;
	gdastat_(&nstep, &x, &y[1], status, (ftnlen)7);

/*     test for convergence to the final integration limit */

	if ((x - *x2) * (*x2 - *x1) >= 0.) {
	    io___11.ciunit = iounit_1.iout;
	    s_wsfe(&io___11);
	    e_wsfe();
	    terminate = TRUE_;
	}

/*     test for a trial stepsize that is too small */

	if (abs(hnext) < *hmin) {
	    io___12.ciunit = iounit_1.iout;
	    s_wsfe(&io___12);
	    e_wsfe();
	    terminate = TRUE_;
	}

/*     test for too many total integration steps */

	if (nstep >= 1000) {
	    io___13.ciunit = iounit_1.iout;
	    s_wsfe(&io___13);
	    e_wsfe();
	    terminate = TRUE_;
	}
    }
    return 0;
} /* diffeq_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine bsstep  --  Bulirsch-Stoer integration step  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "bsstep" takes a single Bulirsch-Stoer step with monitoring */
/*     of local truncation error to ensure accuracy */

/*     literature reference: */

/*     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. */
/*     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge */
/*     University Press, 1992, Section 16.4 */


/* Subroutine */ int bsstep_(integer *nv, doublereal *x, doublereal *dydx, 
	doublereal *y, doublereal *htry, doublereal *eps, doublereal *yscal, 
	doublereal *hdid, doublereal *hnext, S_fp gvalue)
{
    /* Initialized data */

    static logical first = TRUE_;
    static doublereal epsold = -1.;
    static integer nseq[9] = { 2,4,6,8,10,12,14,16,18 };

    /* Format strings */
    static char fmt_30[] = "(\002 BSSTEP  --  Underflow of Step Size\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static doublereal a[9], h__;
    static integer i__, k, kk, km, iq;
    static doublereal alf[64]	/* was [8][8] */, red, err[8], eps1, fact;
    extern /* Subroutine */ int mmid_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, S_fp);
    static integer kmax, kopt;
    static doublereal xnew, work, xest, yerr[100000], ysav[100000], yseq[
	    100000], scale;
    extern /* Subroutine */ int fatal_(void);
    static logical reduct;
    static doublereal errmax, wrkmin;
    extern /* Subroutine */ int pzextr_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___29 = { 0, 0, 0, fmt_30, 0 };



#define alf_ref(a_1,a_2) alf[(a_2)*8 + a_1 - 9]



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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */


    /* Parameter adjustments */
    --yscal;
    --y;
    --dydx;

    /* Function Body */


    if (*eps != epsold) {
	*hnext = -1e29;
	xnew = -1e29;
	eps1 = *eps * .25;
	a[0] = (doublereal) nseq[0] + 1.;
	for (k = 1; k <= 8; ++k) {
	    a[k] = a[k - 1] + (doublereal) nseq[k];
	}
	for (iq = 2; iq <= 8; ++iq) {
	    i__1 = iq - 1;
	    for (k = 1; k <= i__1; ++k) {
		d__1 = (a[k] - a[iq]) / ((a[iq] - a[0] + 1.) * ((k << 1) + 1))
			;
		alf_ref(k, iq) = pow_dd(&eps1, &d__1);
	    }
	}
	epsold = *eps;
	for (kopt = 2; kopt <= 7; ++kopt) {
	    if (a[kopt] > a[kopt - 1] * alf_ref(kopt - 1, kopt)) {
		goto L10;
	    }
	}
L10:
	kmax = kopt;
    }
    h__ = *htry;
    i__1 = *nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ysav[i__ - 1] = y[i__];
    }
    if (h__ != *hnext || *x != xnew) {
	first = TRUE_;
	kopt = kmax;
    }
    reduct = FALSE_;
L20:
    i__1 = kmax;
    for (k = 1; k <= i__1; ++k) {
	xnew = *x + h__;
	if (xnew == *x) {
	    io___29.ciunit = iounit_1.iout;
	    s_wsfe(&io___29);
	    e_wsfe();
	    fatal_();
	}
	mmid_(&nseq[k - 1], &h__, nv, x, &dydx[1], ysav, yseq, (S_fp)gvalue);
/* Computing 2nd power */
	d__1 = h__ / (doublereal) nseq[k - 1];
	xest = d__1 * d__1;
	pzextr_(&k, nv, &xest, yseq, &y[1], yerr);
	if (k != 1) {
	    errmax = 1e-30;
	    i__2 = *nv;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
		d__2 = errmax, d__3 = (d__1 = yerr[i__ - 1] / yscal[i__], abs(
			d__1));
		errmax = max(d__2,d__3);
	    }
	    errmax /= *eps;
	    km = k - 1;
	    d__1 = errmax / .25;
	    d__2 = 1. / ((km << 1) + 1);
	    err[km - 1] = pow_dd(&d__1, &d__2);
	}
	if (k != 1 && (k >= kopt - 1 || first)) {
	    if (errmax < 1.) {
		goto L50;
	    }
	    if (k == kmax || k == kopt + 1) {
		red = .7 / err[km - 1];
		goto L40;
	    } else if (k == kopt) {
		if (alf_ref(kopt - 1, kopt) < err[km - 1]) {
		    red = 1. / err[km - 1];
		    goto L40;
		}
	    } else if (kopt == kmax) {
		if (alf_ref(km, kmax - 1) < err[km - 1]) {
		    red = alf_ref(km, kmax - 1) * .7 / err[km - 1];
		    goto L40;
		}
	    } else if (alf_ref(km, kopt) < err[km - 1]) {
		red = alf_ref(km, kopt - 1) / err[km - 1];
		goto L40;
	    }
	}
    }
L40:
    red = min(red,.7);
    red = max(red,1e-5);
    h__ *= red;
    reduct = TRUE_;
    goto L20;
L50:
    *x = xnew;
    *hdid = h__;
    first = FALSE_;
    wrkmin = 1e35;
    i__1 = km;
    for (kk = 1; kk <= i__1; ++kk) {
/* Computing MAX */
	d__1 = err[kk - 1];
	fact = max(d__1,.1);
	work = fact * a[kk];
	if (work < wrkmin) {
	    scale = fact;
	    wrkmin = work;
	    kopt = kk + 1;
	}
    }
    *hnext = h__ / scale;
    if (kopt >= k && kopt != kmax && ! reduct) {
/* Computing MAX */
	d__1 = scale / alf_ref(kopt - 1, kopt);
	fact = max(d__1,.1);
	if (a[kopt] * fact <= wrkmin) {
	    *hnext = h__ / fact;
	    ++kopt;
	}
    }
    return 0;
} /* bsstep_ */

#undef alf_ref




/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine mmid  --  takes a modified midpoint step  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "mmid" implements a modified midpoint method to advance the */
/*     integration of a set of first order differential equations */


/* Subroutine */ int mmid_(integer *nstep, doublereal *htot, integer *nvar, 
	doublereal *xs, doublereal *dydx, doublereal *y, doublereal *yout, 
	S_fp gvalue)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal h__;
    static integer i__, k;
    static doublereal x, h2, ym[100000], yn[100000], temp;



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




/*     set substep size based on number of steps to be taken */

    /* Parameter adjustments */
    --yout;
    --y;
    --dydx;

    /* Function Body */
    h__ = *htot / (doublereal) (*nstep);
    h2 = h__ * 2.;

/*     take the first substep and get values at ends of step */

    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ym[i__ - 1] = y[i__];
	yn[i__ - 1] = y[i__] + h__ * dydx[i__];
    }
    x = *xs + h__;
    (*gvalue)(&x, yn, &yout[1]);

/*     take the second and subsequent substeps */

    i__1 = *nstep;
    for (k = 2; k <= i__1; ++k) {
	i__2 = *nvar;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    temp = ym[i__ - 1] + h2 * yout[i__];
	    ym[i__ - 1] = yn[i__ - 1];
	    yn[i__ - 1] = temp;
	}
	x += h__;
	(*gvalue)(&x, yn, &yout[1]);
    }

/*     complete the update of values for the last substep */

    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yout[i__] = (ym[i__ - 1] + yn[i__ - 1] + h__ * yout[i__]) * .5;
    }
    return 0;
} /* mmid_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine pzextr  --  polynomial extrapolation method  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "pzextr" is a polynomial extrapolation routine used during */
/*     Bulirsch-Stoer integration of ordinary differential equations */


/* Subroutine */ int pzextr_(integer *iest, integer *nvar, doublereal *xest, 
	doublereal *yest, doublereal *yz, doublereal *dy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal d__[100000];
    static integer i__, j;
    static doublereal q, x[13], f1, f2, qcol[1300000]	/* was [100000][13] */
	    , delta;


#define qcol_ref(a_1,a_2) qcol[(a_2)*100000 + a_1 - 100001]



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




    /* Parameter adjustments */
    --dy;
    --yz;
    --yest;

    /* Function Body */
    x[*iest - 1] = *xest;
    i__1 = *nvar;
    for (j = 1; j <= i__1; ++j) {
	dy[j] = yest[j];
	yz[j] = yest[j];
    }
    if (*iest == 1) {
	i__1 = *nvar;
	for (j = 1; j <= i__1; ++j) {
	    qcol_ref(j, 1) = yest[j];
	}
    } else {
	i__1 = *nvar;
	for (j = 1; j <= i__1; ++j) {
	    d__[j - 1] = yest[j];
	}
	i__1 = *iest - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    delta = 1. / (x[*iest - i__ - 1] - *xest);
	    f1 = *xest * delta;
	    f2 = x[*iest - i__ - 1] * delta;
	    i__2 = *nvar;
	    for (j = 1; j <= i__2; ++j) {
		q = qcol_ref(j, i__);
		qcol_ref(j, i__) = dy[j];
		delta = d__[j - 1] - q;
		dy[j] = f1 * delta;
		d__[j - 1] = f2 * delta;
		yz[j] += dy[j];
	    }
	}
	i__1 = *nvar;
	for (j = 1; j <= i__1; ++j) {
	    qcol_ref(j, *iest) = dy[j];
	}
    }
    return 0;
} /* pzextr_ */

#undef qcol_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine gdastat  --  results for GDA integration step  ## */
/*     ##                                                            ## */
/*     ################################################################ */

/*     "gdastat" finds the energy, radius of gyration, and average M2 */
/*     for a GDA integration step; also saves the coordinates */


/* Subroutine */ int gdastat_(integer *nstep, doublereal *beta, doublereal *
	xx, char *status, ftnlen status_len)
{
    /* Format strings */
    static char fmt_10[] = "(i6,2x,4f13.4,6x,a7)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *);
    double log(doublereal);
    integer do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal e;
    static integer i__;
    static doublereal rg;
    static integer nvar;
    static doublereal m2ave;
    extern doublereal energy_(void);
    extern /* Subroutine */ int gyrate_(doublereal *), optsave_(integer *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___64 = { 0, 0, 0, fmt_10, 0 };




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
	warp_1.m2[i__ - 1] = (d__1 = xx[nvar], abs(d__1));
    }

/*     get some info about the current integration step */

    e = energy_();
    gyrate_(&rg);
    m2ave = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m2ave += warp_1.m2[i__ - 1];
    }
    m2ave /= (doublereal) atoms_1.n;
    io___64.ciunit = iounit_1.iout;
    s_wsfe(&io___64);
    do_fio(&c__1, (char *)&(*nstep), (ftnlen)sizeof(integer));
    d__1 = log(*beta) / 2.302585092994045684;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&rg, (ftnlen)sizeof(doublereal));
    d__2 = log(m2ave) / 2.302585092994045684;
    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, status, (ftnlen)7);
    e_wsfe();

/*     save the current coordinates to a disk file */

    optsave_(nstep, &e, &xx[1]);
    return 0;
} /* gdastat_ */

