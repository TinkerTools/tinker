/* random.f -- translated by f2c (version 20050501).
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

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static doublereal c_b13 = 6.2831853071795862;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  function random  --  portable random number generator  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "random" generates a random number on [0,1] via a long */
/*     period generator due to L'Ecuyer with Bays-Durham shuffle */

/*     literature references: */

/*     P. L'Ecuyer, Communications of the ACM, 31, 742-774 (1988) */

/*     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. */
/*     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge */
/*     University Press, 1992, Section 7.1 */


doublereal random_(void)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* Format strings */
    static char fmt_20[] = "(/,\002 RANDOM  --  Initialized with SEED of\002"
	    ",i12)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int calendar_(integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    static integer ishuffle[32], i__, k, iy, day, seed, year, hour, next, 
	    seed2, month, second;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static integer minute;
    static char string[120], keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___14 = { 1, string, 0, 0, 120, 1 };
    static cilist io___15 = { 0, 0, 0, fmt_20, 0 };



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




/*     random number seed is first set to a big number, */
/*     then incremented by the seconds elapsed this decade */

    if (first) {
	first = FALSE_;
	seed = 141803398;
	calendar_(&year, &month, &day, &hour, &minute, &second);
	year %= 10;
	seed = seed + year * 32140800 + (month - 1) * 2678400;
	seed = seed + (day - 1) * 86400 + hour * 3600;
	seed = seed + minute * 60 + second;

/*     search the keywords for a random number seed */

	i__1 = keys_1.nkey;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    next = 1;
	    s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	    gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	    upcase_(keyword, (ftnlen)20);
	    if (s_cmp(keyword, "RANDOMSEED ", (ftnlen)11, (ftnlen)11) == 0) {
		s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next 
			- 1));
		i__2 = s_rsli(&io___14);
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = do_lio(&c__3, &c__1, (char *)&seed, (ftnlen)sizeof(
			integer));
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L10;
		}
		seed = max(1,seed);
	    }
L10:
	    ;
	}

/*     print the value used for the random number seed */

	if (inform_1.verbose) {
	    io___15.ciunit = iounit_1.iout;
	    s_wsfe(&io___15);
	    do_fio(&c__1, (char *)&seed, (ftnlen)sizeof(integer));
	    e_wsfe();
	}

/*     warm up and then load the shuffling table */

	seed2 = seed;
	for (i__ = 40; i__ >= 1; --i__) {
	    k = seed / 53668;
	    seed = (seed - k * 53668) * 40014 - k * 12211;
	    if (seed < 0) {
		seed += 2147483563;
	    }
	    if (i__ <= 32) {
		ishuffle[i__ - 1] = seed;
	    }
	}
	iy = ishuffle[0];
    }

/*     get a new random number value each call */

    k = seed / 53668;
    seed = (seed - k * 53668) * 40014 - k * 12211;
    if (seed < 0) {
	seed += 2147483563;
    }
    k = seed2 / 52774;
    seed2 = (seed2 - k * 52774) * 40692 - k * 3791;
    if (seed2 < 0) {
	seed2 += 2147483399;
    }
    i__ = iy / 67108862 + 1;
    iy = ishuffle[i__ - 1] - seed2;
    ishuffle[i__ - 1] = seed;
    if (iy < 1) {
	iy += 2147483562;
    }
    ret_val = iy * 4.6566130573917691e-10;

/*     print the value of the current random number */

/*     if (debug) then */
/*        write (iout,30)  random */
/*  30    format (' RANDOM  --  The Random Number Value is',f12.8) */
/*     end if */
    return ret_val;
} /* random_ */

#undef keyline_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  function normal  --  random number from normal curve  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "normal" generates a random number from a normal Gaussian */
/*     distribution with a mean of zero and a variance of one */


doublereal normal_(void)
{
    /* Initialized data */

    static logical compute = TRUE_;

    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal v1, v2, rsq, store, factor;
    extern doublereal random_(void);



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




/*     get a pair of random values from the distribution */

    if (compute) {
L10:
	v1 = random_() * 2. - 1.;
	v2 = random_() * 2. - 1.;
/* Computing 2nd power */
	d__1 = v1;
/* Computing 2nd power */
	d__2 = v2;
	rsq = d__1 * d__1 + d__2 * d__2;
	if (rsq >= 1.) {
	    goto L10;
	}
	factor = sqrt(log(rsq) * -2. / rsq);
	store = v1 * factor;
	ret_val = v2 * factor;
	compute = FALSE_;

/*     use the second random value computed at the last call */

    } else {
	ret_val = store;
	compute = TRUE_;
    }

/*     print the value of the current random number */

/*     if (debug) then */
/*        write (iout,20)  normal */
/*  20    format (' NORMAL  --  The Random Number Value is',f12.8) */
/*     end if */
    return ret_val;
} /* normal_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine ranvec  --  unit vector in random direction  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "ranvec" generates a unit vector in 3-dimensional */
/*     space with uniformly distributed random orientation */

/*     literature references: */

/*     G. Marsaglia, Ann. Math. Stat., 43, 645 (1972) */

/*     R. C. Rapaport, The Art of Molecular Dynamics Simulation, */
/*     2nd Edition, Cambridge University Press, 2004, Section 18.4 */


/* Subroutine */ int ranvec_(doublereal *vector)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal s, x, y;
    extern doublereal random_(void);



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




/*     get a pair of appropriate components in the plane */

    /* Parameter adjustments */
    --vector;

    /* Function Body */
    s = 2.;
    while(s >= 1.) {
	x = random_() * 2. - 1.;
	y = random_() * 2. - 1.;
/* Computing 2nd power */
	d__1 = x;
/* Computing 2nd power */
	d__2 = y;
	s = d__1 * d__1 + d__2 * d__2;
    }

/*     construct the 3-dimensional random unit vector */

    vector[3] = 1. - s * 2.;
    s = sqrt(1. - s) * 2.;
    vector[2] = s * y;
    vector[1] = s * x;

/*     print the components of the random unit vector */

/*     if (debug) then */
/*        write (iout,10)  vector(1),vector(2),vector(3) */
/*  10    format (' RANVEC  --  The Random Vector is',3f10.4) */
/*     end if */
    return 0;
} /* ranvec_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine sphere  --  uniform set of points on sphere  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "sphere" finds a specified number of uniformly distributed */
/*     points on a sphere of unit radius centered at the origin */

/*     literature reference: */

/*     E. B. Saff and A. B. J. Kuijlaars, "Distributing Many */
/*     Points on a Sphere", The Mathematical Intelligencer, */
/*     19, 5-11 (1997) */


/* Subroutine */ int sphere_(integer *ndot, doublereal *dot)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double acos(doublereal), sqrt(doublereal), d_mod(doublereal *, doublereal 
	    *), sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal h__;
    static integer i__;
    static doublereal phi, tot, tot1, theta, phiold;


#define dot_ref(a_1,a_2) dot[(a_2)*3 + a_1]



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




/*     find spherical coordinates then convert to Cartesian */

    /* Parameter adjustments */
    dot -= 4;

    /* Function Body */
    tot = (doublereal) (*ndot);
    tot1 = (doublereal) (*ndot - 1);
    i__1 = *ndot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__ = (doublereal) (i__ - 1) * 2. / tot1 - 1.;
	theta = acos(h__);
	if (i__ == 1 || i__ == *ndot) {
	    phi = 0.;
	} else {
	    d__1 = phiold + 3.6 / sqrt(tot * (1. - h__ * h__));
	    phi = d_mod(&d__1, &c_b13);
	}
	dot_ref(1, i__) = sin(theta) * cos(phi);
	dot_ref(2, i__) = sin(theta) * sin(phi);
	dot_ref(3, i__) = cos(theta);
	phiold = phi;
    }
    return 0;
} /* sphere_ */

#undef dot_ref


