/* search.f -- translated by f2c (version 20050501).
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
    doublereal stpmin, stpmax, cappa, slpmax, angmax;
    integer intmax;
} linmin_;

#define linmin_1 linmin_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine search  --  perform unidimensional line search  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "search" is a unidimensional line search based upon parabolic */
/*     extrapolation and cubic interpolation using both function and */
/*     gradient values */

/*     variables used by the routine : */

/*     f       function value at the best line search point */
/*     x       current values of variables during line search */
/*     g       gradient at the current point during line search */
/*     p       initial search vector, unchanged by this routine */
/*     s       scaled search vector at current line search point */
/*     angle   angle between search and negative gradient vector */

/*     parameters used by the routine : */

/*     stpmin   minimum step length in current line search direction */
/*     stpmax   maximum step length in current line search direction */
/*     cappa    stringency of line search (0=tight < cappa < 1=loose) */
/*     slpmax   projected gradient above which stepsize is reduced */
/*     angmax   maximum angle between search direction and -gradient */
/*     intmax   maximum number of interpolations during line search */

/*     status codes upon return : */

/*     Success     normal termination after satisfying "cappa" test */
/*     ScaleStep   normal termination after a step size rescaling */
/*     ReSearch    normal termination after a reinterpolation */
/*     WideAngle   large angle between search direction and -gradient */
/*     BadIntpln   unsatisfied "cappa" test after two searches */
/*     IntplnErr   function value increase or serious gradient error */


/* Subroutine */ int search_(integer *nvar, doublereal *f, doublereal *g, 
	doublereal *x, doublereal *p, doublereal *f_move__, doublereal *angle,
	 integer *ncalls, D_fp fgvalue, char *status, ftnlen status_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal), acos(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal s[75000], f_0__, f_1__, x_0__[75000], f_a__, f_b__, 
	    f_c__, sss, ttt, sg_0__, sg_1__, sg_a__, sg_b__, sg_c__, cube, 
	    step, parab;
    static char blank[9];
    static doublereal cosang, g_norm__, s_norm__;
    static integer intpln;
    static doublereal cubstp;
    static logical restart;



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




/*     use default parameters for the line search if needed */

    /* Parameter adjustments */
    --p;
    --x;
    --g;

    /* Function Body */
    s_copy(blank, "         ", (ftnlen)9, (ftnlen)9);
    if (linmin_1.stpmin == 0.) {
	linmin_1.stpmin = 1e-20;
    }
    if (linmin_1.stpmax == 0.) {
	linmin_1.stpmax = 2.;
    }
    if (linmin_1.cappa == 0.) {
	linmin_1.cappa = .1;
    }
    if (linmin_1.slpmax == 0.) {
	linmin_1.slpmax = 1e4;
    }
    if (linmin_1.angmax == 0.) {
	linmin_1.angmax = 180.;
    }
    if (linmin_1.intmax == 0) {
	linmin_1.intmax = 5;
    }

/*     copy the search direction into a new vector */

    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s[i__ - 1] = p[i__];
    }

/*     compute the length of gradient and search direction */

    g_norm__ = 0.;
    s_norm__ = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g_norm__ += g[i__] * g[i__];
	s_norm__ += s[i__ - 1] * s[i__ - 1];
    }
    g_norm__ = sqrt(g_norm__);
    s_norm__ = sqrt(s_norm__);

/*     store initial function, then normalize the */
/*     search vector and find projected gradient */

    f_0__ = *f;
    sg_0__ = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x_0__[i__ - 1] = x[i__];
	s[i__ - 1] /= s_norm__;
	sg_0__ += s[i__ - 1] * g[i__];
    }

/*     check the angle between the search direction */
/*     and the negative gradient vector */

    cosang = -sg_0__ / g_norm__;
/* Computing MIN */
    d__1 = 1., d__2 = max(-1.,cosang);
    cosang = min(d__1,d__2);
    *angle = acos(cosang) * 57.29577951308232088;
    if (*angle > linmin_1.angmax) {
	s_copy(status, "WideAngle", (ftnlen)9, (ftnlen)9);
	return 0;
    }

/*     set the initial stepsize to the length of the passed */
/*     search vector, or based on previous function decrease */

    step = (d__1 = *f_move__ / sg_0__, abs(d__1)) * 2.;
    step = min(step,s_norm__);
    if (step > linmin_1.stpmax) {
	step = linmin_1.stpmax;
    }
    if (step < linmin_1.stpmin) {
	step = linmin_1.stpmin;
    }

/*     beginning of the parabolic extrapolation procedure */

L10:
    restart = TRUE_;
    intpln = 0;
    f_b__ = f_0__;
    sg_b__ = sg_0__;

/*     replace last point by latest and take another step */

L20:
    f_a__ = f_b__;
    sg_a__ = sg_b__;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] += step * s[i__ - 1];
    }

/*     get new function and projected gradient following a step */

    ++(*ncalls);
    f_b__ = (*fgvalue)(&x[1], &g[1]);
    sg_b__ = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sg_b__ += s[i__ - 1] * g[i__];
    }

/*     scale stepsize if initial gradient change is too large */

    if ((d__1 = sg_b__ / sg_a__, abs(d__1)) >= linmin_1.slpmax && restart) {
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] = x_0__[i__ - 1];
	}
	step /= 10.;
	s_copy(status, "ScaleStep", (ftnlen)9, (ftnlen)9);
	goto L10;
    }
    restart = FALSE_;

/*     return if the gradient is small and function decreases */

    if ((d__1 = sg_b__ / sg_0__, abs(d__1)) <= linmin_1.cappa && f_b__ < 
	    f_a__) {
	*f = f_b__;
	if (s_cmp(status, blank, (ftnlen)9, (ftnlen)9) == 0) {
	    s_copy(status, " Success ", (ftnlen)9, (ftnlen)9);
	}
	return 0;
    }

/*     interpolate if gradient changes sign or function increases */

    if (sg_b__ * sg_a__ < 0. || f_b__ > f_a__) {
	goto L30;
    }

/*     if the finite difference curvature is negative double the step; */
/*     or if  step < parabolic estimate < 4*step  use this estimate, */
/*     otherwise truncate to step or 4*step, respectively */

    step *= 2.;
    if (sg_b__ > sg_a__) {
	parab = (f_a__ - f_b__) / (sg_b__ - sg_a__);
	if (parab > step * 2.) {
	    parab = step * 2.;
	}
	if (parab < step * .5) {
	    parab = step * .5;
	}
	step = parab;
    }
    if (step > linmin_1.stpmax) {
	step = linmin_1.stpmax;
    }
    goto L20;

/*     beginning of the cubic interpolation procedure */

L30:
    ++intpln;
    sss = (f_b__ - f_a__) * 3. / step - sg_a__ - sg_b__;
    ttt = sss * sss - sg_a__ * sg_b__;
    if (ttt < 0.) {
	*f = f_b__;
	s_copy(status, "IntplnErr", (ftnlen)9, (ftnlen)9);
	return 0;
    }
    ttt = sqrt(ttt);
    cube = step * (sg_b__ + ttt + sss) / (sg_b__ - sg_a__ + ttt * 2.);
    if (cube < 0. || cube > step) {
	*f = f_b__;
	s_copy(status, "IntplnErr", (ftnlen)9, (ftnlen)9);
	return 0;
    }
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] -= cube * s[i__ - 1];
    }

/*     get new function and gradient, then test for termination */

    ++(*ncalls);
    f_c__ = (*fgvalue)(&x[1], &g[1]);
    sg_c__ = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sg_c__ += s[i__ - 1] * g[i__];
    }
    if ((d__1 = sg_c__ / sg_0__, abs(d__1)) <= linmin_1.cappa) {
	*f = f_c__;
	if (s_cmp(status, blank, (ftnlen)9, (ftnlen)9) == 0) {
	    s_copy(status, " Success ", (ftnlen)9, (ftnlen)9);
	}
	return 0;
    }

/*     get the next pair of bracketing points by replacing one */
/*     of the current brackets with the interpolated point */

    if (f_c__ <= f_a__ || f_c__ <= f_b__) {
/* Computing MIN */
	d__2 = abs(cube), d__3 = (d__1 = step - cube, abs(d__1));
	cubstp = min(d__2,d__3);
	if (cubstp >= linmin_1.stpmin && intpln < linmin_1.intmax) {

/*     if the current brackets have slopes of opposite sign, */
/*     then substitute the interpolated point for the bracket */
/*     point with slope of same sign as the interpolated point */

	    if (sg_a__ * sg_b__ < 0.) {
		if (sg_a__ * sg_c__ < 0.) {
		    f_b__ = f_c__;
		    sg_b__ = sg_c__;
		    step -= cube;
		} else {
		    f_a__ = f_c__;
		    sg_a__ = sg_c__;
		    step = cube;
		    i__1 = *nvar;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			x[i__] += cube * s[i__ - 1];
		    }
		}

/*     if current brackets have slope of same sign, then replace */
/*     the far bracket if the interpolated point has a slope of */
/*     the opposite sign or a lower function value than the near */
/*     bracket, otherwise replace the far bracket point */

	    } else {
		if (sg_a__ * sg_c__ < 0. || f_a__ <= f_c__) {
		    f_b__ = f_c__;
		    sg_b__ = sg_c__;
		    step -= cube;
		} else {
		    f_a__ = f_c__;
		    sg_a__ = sg_c__;
		    step = cube;
		    i__1 = *nvar;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			x[i__] += cube * s[i__ - 1];
		    }
		}
	    }
	    goto L30;
	}
    }

/*     interpolation has failed, reset to best current point */

/* Computing MIN */
    d__1 = min(f_a__,f_b__);
    f_1__ = min(d__1,f_c__);
    if (f_1__ == f_a__) {
	sg_1__ = sg_a__;
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] += (cube - step) * s[i__ - 1];
	}
    } else if (f_1__ == f_b__) {
	sg_1__ = sg_b__;
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] += cube * s[i__ - 1];
	}
    } else if (f_1__ == f_c__) {
	sg_1__ = sg_c__;
    }

/*     try to restart from best point with smaller stepsize */

    if (f_1__ > f_0__) {
	++(*ncalls);
	*f = (*fgvalue)(&x[1], &g[1]);
	s_copy(status, "IntplnErr", (ftnlen)9, (ftnlen)9);
	return 0;
    }
    f_0__ = f_1__;
    sg_0__ = sg_1__;
    if (sg_1__ > 0.) {
	i__1 = *nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s[i__ - 1] = -s[i__ - 1];
	}
	sg_0__ = -sg_1__;
    }
/* Computing MAX */
    d__1 = cube, d__2 = step - cube;
    step = max(d__1,d__2) / 10.;
    if (step < linmin_1.stpmin) {
	step = linmin_1.stpmin;
    }

/*     if already restarted once, then return with best point */

    if (s_cmp(status, " ReSearch", (ftnlen)9, (ftnlen)9) == 0) {
	++(*ncalls);
	*f = (*fgvalue)(&x[1], &g[1]);
	s_copy(status, "BadIntpln", (ftnlen)9, (ftnlen)9);
	return 0;
    } else {
	s_copy(status, " ReSearch", (ftnlen)9, (ftnlen)9);
	goto L10;
    }
    return 0;
} /* search_ */

