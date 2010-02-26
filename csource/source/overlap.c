/* overlap.f -- translated by f2c (version 20050501).
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

/* Table of constant values */

static integer c__2 = 2;
static integer c__5 = 5;
static integer c__4 = 4;
static doublereal c_b14 = -1.;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine overlap  --  p-orbital overlap for pisystem  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "overlap" computes the overlap for two parallel p-orbitals */
/*     given the atomic numbers and distance of separation */


/* Subroutine */ int overlap_(integer *atmnum1, integer *atmnum2, doublereal *
	rang, doublereal *ovlap)
{
    /* Initialized data */

    static doublereal zeta[18] = { 1.,1.7,.65,.975,1.3,1.625,1.95,2.275,2.6,
	    2.925,.733,.95,1.167,1.383,1.6,1.817,2.033,2.25 };

    static doublereal s[3];
    static integer la, lb, na, nb;
    static doublereal za, zb, rbohr;
    extern /* Subroutine */ int slater_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *);



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  units.i  --  physical constants and unit conversions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     literature reference: */

/*     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended */
/*     Values of the Fundamental Physical Constants: 2006", Reviews of */
/*     Modern Physics, 80, 633-730 (2008) */

/*     The "2006 CODATA Recommended Values" are also available from */
/*     the NIST Reference on Constants, Units, and Uncertainty site */
/*     at http://physics.nist.gov/cuu/index.html */

/*     Most values below are derived from the 2006 CODATA reference */

/*     The conversion from calorie to Joule is the definition of the */
/*     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992) */

/*     The "coulomb" energy conversion factor is found by dimensional */
/*     analysis of Coulomb's Law, ie, by dividing the square of the */
/*     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is */
/*     the permittivity of vacuum (the "electric constant"); note that */
/*     eps0 is typically given in F/m, equivalent to C**2/(J-m) */

/*     The approximate value used for the Debye, 3.33564 x 10-30 C-m, */
/*     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997) */

/*     The value of "prescon" is based on definition of 1 atmosphere */
/*     as 101325 Pa set by the 10th Conference Generale des Poids et */
/*     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3 */

/*     avogadro    Avogadro's number (N) in particles/mole */
/*     lightspd    speed of light in vacuum (c) in cm/ps */
/*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K */
/*     gasconst    ideal gas constant (R) in kcal/mole/K */
/*     emass       mass of an electron in atomic mass units */
/*     joule       conversion from calories to joules */
/*     convert     conversion from kcal to g*Ang**2/ps**2 */
/*     bohr        conversion from Bohrs to Angstroms */
/*     hartree     conversion from Hartree to kcal/mole */
/*     evolt       conversion from Hartree to electron-volts */
/*     efreq       conversion from Hartree to cm-1 */
/*     coulomb     conversion from electron**2/Ang to kcal/mole */
/*     debye       conversion from electron-Ang to Debyes */
/*     prescon     conversion from kcal/mole/Ang**3 to Atm */



/*     Slater orbital exponents for hydrogen through argon */



/*     principal quantum number from atomic number */

    na = 2;
    nb = 2;
    if (*atmnum1 > 10) {
	na = 3;
    }
    if (*atmnum2 > 10) {
	nb = 3;
    }

/*     azimuthal quantum number for p-orbitals */

    la = 1;
    lb = 1;

/*     orbital exponent from stored ideal values */

    za = zeta[*atmnum1 - 1];
    zb = zeta[*atmnum2 - 1];

/*     convert interatomic distance to bohrs */

    rbohr = *rang / .52917720859;

/*     get pi-overlap via generic overlap integral routine */

    slater_(&na, &la, &za, &nb, &lb, &zb, &rbohr, s);
    *ovlap = s[1];
    return 0;
} /* overlap_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine slater  --  find overlap integrals for STO's  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "slater" is a general routine for computing the overlap */
/*     integrals between two Slater-type orbitals */

/*     literature reference: */

/*     D. B. Cook, "Structures and Approximations for Electrons in */
/*     Molecules", Ellis Horwood Limited, Sussex, England, 1978, */
/*     adapted from the code in Chapter 7 */

/*     variables and parameters: */

/*     na   principle quantum number for first orbital */
/*     la   azimuthal quantum number for first orbital */
/*     za   orbital exponent for the first orbital */
/*     nb   principle quantum number for second orbital */
/*     lb   azimuthal quantum number for second orbital */
/*     zb   orbital exponent for the second orbital */
/*     r    interatomic distance in atomic units */
/*     s    vector containing the sigma-sigma, pi-pi */
/*            and delta-delta overlaps upon output */


/* Subroutine */ int slater_(integer *na, integer *la, doublereal *za, 
	integer *nb, integer *lb, doublereal *zb, doublereal *r__, doublereal 
	*s)
{
    /* Initialized data */

    static integer icosa[2] = { 0,1 };
    static integer icosb[2] = { 0,1 };
    static doublereal cosa[2] = { 1.,1. };
    static doublereal cosb[2] = { -1.,1. };
    static integer idsga[5] = { 0,1,2,2,0 };
    static integer idsgb[5] = { 0,1,2,0,2 };
    static doublereal dsiga[5] = { 3.,4.,3.,-1.,-1. };
    static doublereal dsigb[5] = { 3.,-4.,3.,-1.,-1. };
    static integer isina[4] = { 0,2,0,2 };
    static integer isinb[4] = { 0,0,2,2 };
    static doublereal sinab[4] = { -1.,1.,1.,-1. };
    static doublereal theta[6] = { .7071068,1.224745,.8660254,.7905694,
	    1.9364916,.9682458 };
    static doublereal fact[15] = { 1.,1.,2.,6.,24.,120.,720.,5040.,40320.,
	    362880.,3628800.,39916800.,479001600.,6227020800.,87178291200. };

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_di(doublereal *, integer *), sqrt(doublereal);

    /* Local variables */
    static doublereal a[20], b[20], c__[200];
    static integer j, k, m;
    static doublereal p;
    static integer ia[200], ja, jb, ib[200];
    static doublereal an;
    static integer nn;
    static doublereal pt, ana, anb, anr;
    static integer max__;
    static doublereal coef;
    extern doublereal cjkm_(integer *, integer *, integer *);
    static logical done;
    extern /* Subroutine */ int aset_(doublereal *, integer *, doublereal *), 
	    bset_(doublereal *, integer *, doublereal *);
    static integer novi, maxx;
    static doublereal cbase[20], rhalf;
    extern /* Subroutine */ int polyp_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *);

    /* Parameter adjustments */
    --s;

    /* Function Body */


/*     zero out the overlap integrals */

    done = FALSE_;
    s[1] = 0.;
    s[2] = 0.;
    s[3] = 0.;
    d__1 = *za * 2.;
    i__1 = (*na << 1) + 1;
    ana = pow_di(&d__1, &i__1) / fact[*na * 2];
    d__1 = *zb * 2.;
    i__1 = (*nb << 1) + 1;
    anb = pow_di(&d__1, &i__1) / fact[*nb * 2];

/*     orbitals are on the same atomic center */

    if (*r__ < 1e-6) {
	anr = 1.;
	j = *na + *nb + 1;
	d__1 = *za + *zb;
	s[1] = fact[j - 1] / pow_di(&d__1, &j);
	an = sqrt(ana * anb);
	for (novi = 1; novi <= 3; ++novi) {
	    s[novi] = s[novi] * an * anr;
	}
	return 0;
    }

/*     compute overlap integrals for general case */

    rhalf = *r__ * .5;
    p = rhalf * (*za + *zb);
    pt = rhalf * (*za - *zb);
    nn = *na + *nb;
    aset_(&p, &nn, a);
    bset_(&pt, &nn, b);
    k = *na - *la;
    m = *nb - *lb;
    max__ = k + m + 1;
    i__1 = max__;
    for (j = 1; j <= i__1; ++j) {
	ia[j - 1] = j - 1;
	ib[j - 1] = max__ - j;
	i__2 = j - 1;
	cbase[j - 1] = cjkm_(&i__2, &k, &m);
	c__[j - 1] = cbase[j - 1];
    }
    maxx = max__;
    if (*la == 1) {
	polyp_(c__, ia, ib, &maxx, cosa, icosa, icosb, &c__2);
    } else if (*la == 2) {
	polyp_(c__, ia, ib, &maxx, dsiga, idsga, idsgb, &c__5);
    }
    if (*lb == 1) {
	polyp_(c__, ia, ib, &maxx, cosb, icosa, icosb, &c__2);
    } else if (*lb == 2) {
	polyp_(c__, ia, ib, &maxx, dsigb, idsga, idsgb, &c__5);
    }
    novi = 1;
    while(! done) {
	i__1 = maxx;
	for (j = 1; j <= i__1; ++j) {
	    ja = ia[j - 1] + 1;
	    jb = ib[j - 1] + 1;
	    coef = c__[j - 1];
	    if (abs(coef) >= 1e-8) {
		s[novi] += coef * a[ja - 1] * b[jb - 1];
	    }
	}
	ja = *la * (*la + 1) / 2 + novi;
	jb = *lb * (*lb + 1) / 2 + novi;
	s[novi] = s[novi] * theta[ja - 1] * theta[jb - 1];
	if (novi == 1 && *la != 0 && *lb != 0) {
	    maxx = max__;
	    i__1 = maxx;
	    for (j = 1; j <= i__1; ++j) {
		c__[j - 1] = cbase[j - 1];
	    }
	    polyp_(c__, ia, ib, &maxx, sinab, isina, isinb, &c__4);
	    if (*la == 2) {
		polyp_(c__, ia, ib, &maxx, cosa, icosa, icosb, &c__2);
	    }
	    if (*lb == 2) {
		polyp_(c__, ia, ib, &maxx, cosb, icosa, icosb, &c__2);
	    }
	    novi = 2;
	} else if (novi == 2 && *la == 2 && *lb == 2) {
	    maxx = max__;
	    i__1 = maxx;
	    for (j = 1; j <= i__1; ++j) {
		c__[j - 1] = cbase[j - 1];
	    }
	    polyp_(c__, ia, ib, &maxx, sinab, isina, isinb, &c__4);
	    polyp_(c__, ia, ib, &maxx, sinab, isina, isinb, &c__4);
	    novi = 3;
	} else {
	    i__1 = *na + *nb + 1;
	    anr = pow_di(&rhalf, &i__1);
	    an = sqrt(ana * anb);
	    for (novi = 1; novi <= 3; ++novi) {
		s[novi] = s[novi] * an * anr;
	    }
	    done = TRUE_;
	}
    }
    return 0;
} /* slater_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine polyp  --  polynomial product for STO overlap  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "polyp" is a polynomial product routine that multiplies two */
/*     algebraic forms */


/* Subroutine */ int polyp_(doublereal *c__, integer *ia, integer *ib, 
	integer *max__, doublereal *d__, integer *iaa, integer *ibb, integer *
	n)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, m;



    /* Parameter adjustments */
    --c__;
    --ia;
    --ib;
    --ibb;
    --iaa;
    --d__;

    /* Function Body */
    i__1 = *max__;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {
	    i__ = *n - k + 1;
	    m = (i__ - 1) * *max__ + j;
	    c__[m] = c__[j] * d__[i__];
	    ia[m] = ia[j] + iaa[i__];
	    ib[m] = ib[j] + ibb[i__];
	}
    }
    *max__ = *n * *max__;
    return 0;
} /* polyp_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function cjkm  --  coefficients of spherical harmonics  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "cjkm" computes the coefficients of spherical harmonics */
/*     expressed in prolate spheroidal coordinates */


doublereal cjkm_(integer *j, integer *k, integer *m)
{
    /* Initialized data */

    static doublereal fact[15] = { 1.,1.,2.,6.,24.,120.,720.,5040.,40320.,
	    362880.,3628800.,39916800.,479001600.,6227020800.,87178291200. };

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static integer i__;
    static doublereal b1, b2;
    static integer id, ip1, idd, min__, max__;
    static doublereal sum;



    min__ = 1;
    if (*j > *m) {
	min__ = *j - *m + 1;
    }
    max__ = *j + 1;
    if (*k < *j) {
	max__ = *k + 1;
    }
    sum = 0.;
    i__1 = max__;
    for (ip1 = min__; ip1 <= i__1; ++ip1) {
	i__ = ip1 - 1;
	id = *k - i__ + 1;
	b1 = fact[*k] / (fact[i__] * fact[id - 1]);
	if (*j < i__) {
	    b2 = 1.;
	} else {
	    id = *m - (*j - i__) + 1;
	    idd = *j - i__ + 1;
	    b2 = fact[*m] / (fact[idd - 1] * fact[id - 1]);
	}
	sum += b1 * b2 * pow_di(&c_b14, &i__);
    }
    i__1 = *m - *j;
    ret_val = sum * pow_di(&c_b14, &i__1);
    return ret_val;
} /* cjkm_ */



/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine aset  --  get "A" functions by recursion  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "aset" computes by recursion the A functions used in the */
/*     evaluation of Slater-type (STO) overlap integrals */


/* Subroutine */ int aset_(doublereal *alpha, integer *n, doublereal *a)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal alp;



    /* Parameter adjustments */
    --a;

    /* Function Body */
    alp = 1. / *alpha;
    a[1] = exp(-(*alpha)) * alp;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__ + 1] = a[1] + (doublereal) i__ * a[i__] * alp;
    }
    return 0;
} /* aset_ */



/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine bset  --  get "B" functions by recursion  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "bset" computes by downward recursion the B functions used */
/*     in the evaluation of Slater-type (STO) overlap integrals */


/* Subroutine */ int bset_(doublereal *beta, integer *n, doublereal *b)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal d1, d2;
    extern doublereal bmax_(doublereal *, integer *);
    static doublereal betam;



    /* Parameter adjustments */
    --b;

    /* Function Body */
    if (abs(*beta) < 1e-6) {
	i__1 = *n + 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    b[i__] = 2. / (doublereal) i__;
	    if (i__ / 2 << 1 == i__) {
		b[i__] = 0.;
	    }
	}
    } else if (abs(*beta) > (doublereal) (*n) / 2.3) {
	d1 = exp(*beta);
	d2 = 1. / d1;
	betam = 1. / *beta;
	b[1] = (d1 - d2) * betam;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d1 = -d1;
	    b[i__ + 1] = (d1 - d2 + (doublereal) i__ * b[i__]) * betam;
	}
    } else {
	b[*n + 1] = bmax_(beta, n);
	d1 = exp(*beta);
	d2 = 1. / d1;
	if (*n / 2 << 1 != *n) {
	    d1 = -d1;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j = *n - i__ + 1;
	    d1 = -d1;
	    b[j] = (d1 + d2 + *beta * b[j + 1]) / (doublereal) j;
	}
    }
    return 0;
} /* bset_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function bmax  --  find maximum order of "B" functions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "bmax" computes the maximum order of the B functions needed */
/*     for evaluation of Slater-type (STO) overlap integrals */


doublereal bmax_(doublereal *beta, integer *n)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal b, fi, bot, top, sum;
    static logical done;
    static doublereal sign, term;



    done = FALSE_;
/* Computing 2nd power */
    d__1 = *beta;
    b = d__1 * d__1;
    top = (doublereal) (*n) + 1.;
    sum = 1. / top;
    fi = 2.;
    sign = 2.;
    if (*n / 2 << 1 != *n) {
	top += 1.;
	sum = *beta / top;
	fi += 1.;
	sign = -2.;
    }
    term = sum;
    while(! done) {
	bot = top + 2.;
	term = term * b * top / (fi * (fi - 1.) * bot);
	sum += term;
	if (abs(term) <= 1e-7) {
	    done = TRUE_;
	} else {
	    fi += 2.;
	    top = bot;
	}
    }
    ret_val = sign * sum;
    return ret_val;
} /* bmax_ */

