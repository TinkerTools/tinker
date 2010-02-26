/* surface.f -- translated by f2c (version 20050501).
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
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine surface  --  find the accessible surface area  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "surface" performs an analytical computation of the weighted */
/*     solvent accessible surface area of each atom and the first */
/*     derivatives of the area with respect to Cartesian coordinates */

/*     literature references: */

/*     T. J. Richmond, "Solvent Accessible Surface Area and */
/*     Excluded Volume in Proteins", Journal of Molecular Biology, */
/*     178, 63-89 (1984) */

/*     L. Wesson and D. Eisenberg, "Atomic Solvation Parameters */
/*     Applied to Molecular Dynamics of Proteins in Solution", */
/*     Protein Science, 1, 227-235 (1992) */

/*     variables and parameters: */

/*     total    total surface area of the whole structure */
/*     area     accessible surface area of each atom */
/*     radius   radii of the individual atoms */
/*     weight   weight assigned to each atom's area; if set to */
/*                1.0, return is actual area in square Angstroms */
/*     probe    radius of the probe sphere */
/*     delta    tolerance used in the tests for sphere overlaps */
/*                and for colinearity */
/*     rmove    connectivity errors can usually be avoided if the */
/*                offending atom is shifted by this small amount */


/* Subroutine */ int surface_(doublereal *total, doublereal *area, doublereal 
	*radius, doublereal *weight, doublereal *probe)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 SURFACE  --  Increase the Value of MAX"
	    "ARC\002)";
    static char fmt_70[] = "(/,\002 SURFACE  --  Increase the Value\002,\002"
	    " of MAXARC\002)";
    static char fmt_80[] = "(/,\002 SURFACE  --  Increase the Value\002,\002"
	    " of MAXARC\002)";
    static char fmt_140[] = "(/,\002 SURFACE  --  Connectivity Error at Ato"
	    "m\002,i6)";
    static char fmt_170[] = "(/,\002 SURFACE  --  Negative Area at Atom\002,"
	    "i6)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);
    double asin(doublereal), acos(doublereal);
    integer do_fio(integer *, char *, ftnlen);
    double d_mod(doublereal *, doublereal *);

    /* Local variables */
    static doublereal b[1000];
    static integer i__, j, k, l, m;
    static doublereal r__[25000], t, b1[1000], t1, cc, bg[1000];
    static integer ib, jb;
    static doublereal bk, dk, gi;
    static integer ii;
    static doublereal gk;
    static integer mi, ni, io;
    static doublereal td;
    static integer ir;
    static doublereal tb, ti, tf, ri[1000];
    static integer lt[1000];
    static doublereal ex[1000], gr[1000], xc[1000], rr, yc[1000], tr, zc[1000]
	    , tt, xr, yr, tx, ty, tz, zr, ux[1000], uy[1000], uz[1000], xc1[
	    1000], yc1[1000], zc1[1000], tk1, tk2, tr2, the, rik, bsq[1000], 
	    eps;
    static integer key[1000];
    static doublereal dsq[1000], txb, tyb, axx, axy, axz, ayx, ayy, azx, azy, 
	    azz, uxl, uyl, uzl, txk, tyk, txr, tyr, tzk, txl, tyl, tzl;
    static logical top;
    static doublereal pid2, bsq1[1000], dsq1[1000], pix2, pix4, rrx2, arcf[
	    1000], arci[1000];
    static integer narc;
    static doublereal thec, ccsq, bsqk;
    static integer kent[1000];
    static doublereal ther[1000], dsql;
    static logical skip[25000], omit[1000];
    static doublereal wght, risq[1000];
    static integer kout[1000];
    static doublereal rrsq, xysq;
    extern /* Subroutine */ int sort2_(integer *, doublereal *, integer *), 
	    fatal_(void);
    static doublereal delta;
    static integer intag[1000];
    static doublereal exang;
    static logical moved;
    static doublereal therk;
    static logical komit;
    static doublereal rmove, risqk, rplus, delta2;
    static integer intag1[1000];
    static doublereal arclen, cosine, arcsum, rminus;

    /* Fortran I/O blocks */
    static cilist io___42 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___104 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___116 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___121 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___122 = { 0, 0, 0, fmt_170, 0 };




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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     zero the area and derivatives, and set the sphere radii */

    /* Parameter adjustments */
    --weight;
    --radius;
    --area;

    /* Function Body */
    *total = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	area[i__] = 0.;
	r__[i__ - 1] = radius[i__];
	if (r__[i__ - 1] != 0.) {
	    r__[i__ - 1] += *probe;
	}
    }

/*     set pi multiples, overlap criterion and tolerances */

    pix2 = 6.2831853071795862;
    pix4 = 12.566370614359172;
    pid2 = 1.5707963267948966;
    delta = 1e-8;
/* Computing 2nd power */
    d__1 = delta;
    delta2 = d__1 * d__1;
    eps = 1e-8;
    rmove = 1e-8;

/*     set the "skip" array to exclude all inactive atoms */
/*     that do not overlap any of the current active atoms */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	skip[i__ - 1] = TRUE_;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    xr = atoms_1.x[i__ - 1];
	    yr = atoms_1.y[i__ - 1];
	    zr = atoms_1.z__[i__ - 1];
	    rr = r__[i__ - 1];
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
		d__1 = rr + r__[k - 1];
		rplus = d__1 * d__1;
/* Computing 2nd power */
		d__1 = atoms_1.x[k - 1] - xr;
/* Computing 2nd power */
		d__2 = atoms_1.y[k - 1] - yr;
/* Computing 2nd power */
		d__3 = atoms_1.z__[k - 1] - zr;
		ccsq = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
		if (ccsq <= rplus) {
		    skip[k - 1] = FALSE_;
		}
	    }
	}
    }

/*     compute the area and derivatives of current "ir" sphere */

    i__1 = atoms_1.n;
    for (ir = 1; ir <= i__1; ++ir) {
	if (skip[ir - 1]) {
	    goto L180;
	}
	xr = atoms_1.x[ir - 1];
	yr = atoms_1.y[ir - 1];
	zr = atoms_1.z__[ir - 1];
	rr = r__[ir - 1];
	rrx2 = rr * 2.;
/* Computing 2nd power */
	d__1 = rr;
	rrsq = d__1 * d__1;
	wght = weight[ir];
	moved = FALSE_;

/*     initialize some counters and sums for the "ir" sphere */

L10:
	io = 0;
	jb = 0;
	ib = 0;
	arclen = 0.;
	exang = 0.;

/*     test each sphere to see if it overlaps the "ir" sphere */

	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (i__ == ir) {
		goto L30;
	    }
	    rplus = rr + r__[i__ - 1];
	    tx = atoms_1.x[i__ - 1] - xr;
	    if (abs(tx) >= rplus) {
		goto L30;
	    }
	    ty = atoms_1.y[i__ - 1] - yr;
	    if (abs(ty) >= rplus) {
		goto L30;
	    }
	    tz = atoms_1.z__[i__ - 1] - zr;
	    if (abs(tz) >= rplus) {
		goto L30;
	    }

/*     check for overlap of spheres by testing center to */
/*     center distance against sum and difference of radii */

/* Computing 2nd power */
	    d__1 = tx;
/* Computing 2nd power */
	    d__2 = ty;
	    xysq = d__1 * d__1 + d__2 * d__2;
	    if (xysq < delta2) {
		tx = delta;
		ty = 0.;
		xysq = delta2;
	    }
/* Computing 2nd power */
	    d__1 = tz;
	    ccsq = xysq + d__1 * d__1;
	    cc = sqrt(ccsq);
	    if (rplus - cc <= delta) {
		goto L30;
	    }
	    rminus = rr - r__[i__ - 1];

/*     check for a completely buried "ir" sphere */

	    if (cc - abs(rminus) <= delta) {
		if (rminus <= 0.) {
		    goto L180;
		}
		goto L30;
	    }

/*     calculate overlap parameters between "i" and "ir" sphere */

	    ++io;
	    xc1[io - 1] = tx;
	    yc1[io - 1] = ty;
	    zc1[io - 1] = tz;
	    dsq1[io - 1] = xysq;
	    bsq1[io - 1] = ccsq;
	    b1[io - 1] = cc;
	    gr[io - 1] = (ccsq + rplus * rminus) / (rrx2 * b1[io - 1]);
	    intag1[io - 1] = i__;
	    if (io > 1000) {
		io___42.ciunit = iounit_1.iout;
		s_wsfe(&io___42);
		e_wsfe();
		fatal_();
	    }
L30:
	    ;
	}

/*     case where no other spheres overlap the current sphere */

	if (io == 0) {
	    area[ir] = pix4;
	    goto L160;
	}

/*     case where only one sphere overlaps the current sphere */

	if (io == 1) {
	    k = 1;
	    txk = xc1[0];
	    tyk = yc1[0];
	    tzk = zc1[0];
	    bsqk = bsq1[0];
	    bk = b1[0];
	    intag[0] = intag1[0];
	    arcsum = pix2;
	    ++ib;
	    arclen += gr[k - 1] * arcsum;
	    goto L150;
	}

/*     general case where more than one sphere intersects the */
/*     current sphere; sort intersecting spheres by their degree */
/*     of overlap with the current main sphere */

	sort2_(&io, gr, key);
	i__2 = io;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    k = key[i__ - 1];
	    intag[i__ - 1] = intag1[k - 1];
	    xc[i__ - 1] = xc1[k - 1];
	    yc[i__ - 1] = yc1[k - 1];
	    zc[i__ - 1] = zc1[k - 1];
	    dsq[i__ - 1] = dsq1[k - 1];
	    b[i__ - 1] = b1[k - 1];
	    bsq[i__ - 1] = bsq1[k - 1];
	    omit[i__ - 1] = FALSE_;
	}

/*     radius of the each circle on the surface of the "ir" sphere */

	i__2 = io;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    gi = gr[i__ - 1] * rr;
	    bg[i__ - 1] = b[i__ - 1] * gi;
/* Computing 2nd power */
	    d__1 = gi;
	    risq[i__ - 1] = rrsq - d__1 * d__1;
	    ri[i__ - 1] = sqrt(risq[i__ - 1]);
/* Computing MIN */
/* Computing MAX */
	    d__3 = -1., d__4 = gr[i__ - 1];
	    d__1 = 1., d__2 = max(d__3,d__4);
	    ther[i__ - 1] = pid2 - asin((min(d__1,d__2)));
	}

/*     find boundary of inaccessible area on "ir" sphere */

	i__2 = io - 1;
	for (k = 1; k <= i__2; ++k) {
	    if (! omit[k - 1]) {
		txk = xc[k - 1];
		tyk = yc[k - 1];
		tzk = zc[k - 1];
		bk = b[k - 1];
		therk = ther[k - 1];
		i__3 = io;
		for (j = k + 1; j <= i__3; ++j) {
		    if (omit[j - 1]) {
			goto L60;
		    }

/*     check to see if J circle is intersecting K circle; */
/*     get distance between circle centers and sum of radii */

		    cc = (txk * xc[j - 1] + tyk * yc[j - 1] + tzk * zc[j - 1])
			     / (bk * b[j - 1]);
/* Computing MIN */
		    d__1 = 1., d__2 = max(-1.,cc);
		    cc = acos((min(d__1,d__2)));
		    td = therk + ther[j - 1];

/*     check to see if circles enclose separate regions */

		    if (cc >= td) {
			goto L60;
		    }

/*     check for circle J completely inside circle K */

		    if (cc + ther[j - 1] < therk) {
			goto L40;
		    }

/*     check for circles essentially parallel */

		    if (cc > delta) {
			goto L50;
		    }
L40:
		    omit[j - 1] = TRUE_;
		    goto L60;

/*     check for "ir" sphere completely buried */

L50:
		    if (pix2 - cc <= td) {
			goto L180;
		    }
L60:
		    ;
		}
	    }
	}

/*     find T value of circle intersections */

	i__2 = io;
	for (k = 1; k <= i__2; ++k) {
	    if (omit[k - 1]) {
		goto L110;
	    }
	    komit = omit[k - 1];
	    omit[k - 1] = TRUE_;
	    narc = 0;
	    top = FALSE_;
	    txk = xc[k - 1];
	    tyk = yc[k - 1];
	    tzk = zc[k - 1];
	    dk = sqrt(dsq[k - 1]);
	    bsqk = bsq[k - 1];
	    bk = b[k - 1];
	    gk = gr[k - 1] * rr;
	    risqk = risq[k - 1];
	    rik = ri[k - 1];
	    therk = ther[k - 1];

/*     rotation matrix elements */

	    t1 = tzk / (bk * dk);
	    axx = txk * t1;
	    axy = tyk * t1;
	    axz = dk / bk;
	    ayx = tyk / dk;
	    ayy = txk / dk;
	    azx = txk / bk;
	    azy = tyk / bk;
	    azz = tzk / bk;
	    i__3 = io;
	    for (l = 1; l <= i__3; ++l) {
		if (! omit[l - 1]) {
		    txl = xc[l - 1];
		    tyl = yc[l - 1];
		    tzl = zc[l - 1];

/*     rotate spheres so K vector colinear with z-axis */

		    uxl = txl * axx + tyl * axy - tzl * axz;
		    uyl = tyl * ayy - txl * ayx;
		    uzl = txl * azx + tyl * azy + tzl * azz;
/* Computing MIN */
/* Computing MAX */
		    d__3 = -1., d__4 = uzl / b[l - 1];
		    d__1 = 1., d__2 = max(d__3,d__4);
		    cosine = min(d__1,d__2);
		    if (acos(cosine) < therk + ther[l - 1]) {
/* Computing 2nd power */
			d__1 = uxl;
/* Computing 2nd power */
			d__2 = uyl;
			dsql = d__1 * d__1 + d__2 * d__2;
			tb = uzl * gk - bg[l - 1];
			txb = uxl * tb;
			tyb = uyl * tb;
			td = rik * dsql;
/* Computing 2nd power */
			d__1 = tb;
			tr2 = risqk * dsql - d__1 * d__1;
			tr2 = max(eps,tr2);
			tr = sqrt(tr2);
			txr = uxl * tr;
			tyr = uyl * tr;

/*     get T values of intersection for K circle */

			tb = (txb + tyr) / td;
/* Computing MIN */
			d__1 = 1., d__2 = max(-1.,tb);
			tb = min(d__1,d__2);
			tk1 = acos(tb);
			if (tyb - txr < 0.) {
			    tk1 = pix2 - tk1;
			}
			tb = (txb - tyr) / td;
/* Computing MIN */
			d__1 = 1., d__2 = max(-1.,tb);
			tb = min(d__1,d__2);
			tk2 = acos(tb);
			if (tyb + txr < 0.) {
			    tk2 = pix2 - tk2;
			}
			thec = (rrsq * uzl - gk * bg[l - 1]) / (rik * ri[l - 
				1] * b[l - 1]);
			if (abs(thec) < 1.) {
			    the = -acos(thec);
			} else if (thec >= 1.) {
			    the = 0.;
			} else if (thec <= -1.) {
			    the = -3.141592653589793238;
			}

/*     see if "tk1" is entry or exit point; check t=0 point; */
/*     "ti" is exit point, "tf" is entry point */

/* Computing MIN */
/* Computing MAX */
			d__3 = -1., d__4 = (uzl * gk - uxl * rik) / (b[l - 1] 
				* rr);
			d__1 = 1., d__2 = max(d__3,d__4);
			cosine = min(d__1,d__2);
			if ((acos(cosine) - ther[l - 1]) * (tk2 - tk1) <= 0.) 
				{
			    ti = tk2;
			    tf = tk1;
			} else {
			    ti = tk1;
			    tf = tk2;
			}
			++narc;
			if (narc >= 1000) {
			    io___104.ciunit = iounit_1.iout;
			    s_wsfe(&io___104);
			    e_wsfe();
			    fatal_();
			}
			if (tf <= ti) {
			    arcf[narc - 1] = tf;
			    arci[narc - 1] = 0.;
			    tf = pix2;
			    lt[narc - 1] = l;
			    ex[narc - 1] = the;
			    top = TRUE_;
			    ++narc;
			}
			arcf[narc - 1] = tf;
			arci[narc - 1] = ti;
			lt[narc - 1] = l;
			ex[narc - 1] = the;
			ux[l - 1] = uxl;
			uy[l - 1] = uyl;
			uz[l - 1] = uzl;
		    }
		}
	    }
	    omit[k - 1] = komit;

/*     special case; K circle without intersections */

	    if (narc <= 0) {
		goto L90;
	    }

/*     general case; sum up arclength and set connectivity code */

	    sort2_(&narc, arci, key);
	    arcsum = arci[0];
	    mi = key[0];
	    t = arcf[mi - 1];
	    ni = mi;
	    if (narc > 1) {
		i__3 = narc;
		for (j = 2; j <= i__3; ++j) {
		    m = key[j - 1];
		    if (t < arci[j - 1]) {
			arcsum = arcsum + arci[j - 1] - t;
			exang += ex[ni - 1];
			++jb;
			if (jb >= 1000) {
			    io___116.ciunit = iounit_1.iout;
			    s_wsfe(&io___116);
			    e_wsfe();
			    fatal_();
			}
			l = lt[ni - 1];
			kent[jb - 1] = l * 1000 + k;
			l = lt[m - 1];
			kout[jb - 1] = k * 1000 + l;
		    }
		    tt = arcf[m - 1];
		    if (tt >= t) {
			t = tt;
			ni = m;
		    }
		}
	    }
	    arcsum = arcsum + pix2 - t;
	    if (! top) {
		exang += ex[ni - 1];
		++jb;
		l = lt[ni - 1];
		kent[jb - 1] = l * 1000 + k;
		l = lt[mi - 1];
		kout[jb - 1] = k * 1000 + l;
	    }
	    goto L100;
L90:
	    arcsum = pix2;
	    ++ib;
L100:
	    arclen += gr[k - 1] * arcsum;
L110:
	    ;
	}
	if (arclen == 0.) {
	    goto L180;
	}
	if (jb == 0) {
	    goto L150;
	}

/*     find number of independent boundaries and check connectivity */

	j = 0;
	i__2 = jb;
	for (k = 1; k <= i__2; ++k) {
	    if (kout[k - 1] != 0) {
		i__ = k;
L120:
		m = kout[i__ - 1];
		kout[i__ - 1] = 0;
		++j;
		i__3 = jb;
		for (ii = 1; ii <= i__3; ++ii) {
		    if (m == kent[ii - 1]) {
			if (ii == k) {
			    ++ib;
			    if (j == jb) {
				goto L150;
			    }
			    goto L130;
			}
			i__ = ii;
			goto L120;
		    }
		}
L130:
		;
	    }
	}
	++ib;

/*     attempt to fix connectivity error by moving atom slightly */

	if (moved) {
	    io___121.ciunit = iounit_1.iout;
	    s_wsfe(&io___121);
	    do_fio(&c__1, (char *)&ir, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    moved = TRUE_;
	    xr += rmove;
	    yr += rmove;
	    zr += rmove;
	    goto L10;
	}

/*     form the accessible area for the current atom */

L150:
	area[ir] = ib * pix2 + exang + arclen;
	area[ir] = d_mod(&area[ir], &pix4);
L160:
	area[ir] *= rrsq;

/*     attempt to fix negative area by moving atom slightly */

	if (area[ir] < 0.) {
	    if (moved) {
		io___122.ciunit = iounit_1.iout;
		s_wsfe(&io___122);
		do_fio(&c__1, (char *)&ir, (ftnlen)sizeof(integer));
		e_wsfe();
	    } else {
		moved = TRUE_;
		xr += rmove;
		yr += rmove;
		zr += rmove;
		goto L10;
	    }
	}

/*     weight the accessible area by the scale factor */

	area[ir] *= wght;
	*total += area[ir];
L180:
	;
    }

/*     print out the surface area values for each atom */

/*     if (debug) then */
/*        write (iout,190) */
/* 190    format (/,' Weighted Atomic Surface Areas Values :', */
/*    &           //,4x,'Atom',7x,'Area Term',6x,'Weight',/) */
/*        do i = 1, n */
/*           if (.not. skip(i)) then */
/*              write (iout,200)  i,area(i),weight(i) */
/* 200          format (i8,4x,2f12.4) */
/*           end if */
/*        end do */
/*        write (iout,210)  total */
/* 210    format (/,' Total Weighted Surface Area :',5x,f16.4) */
/*     end if */
    return 0;
} /* surface_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine surface1  --  accessible surface area & derivs  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "surface1" performs an analytical computation of the weighted */
/*     solvent accessible surface area of each atom and the first */
/*     derivatives of the area with respect to Cartesian coordinates */

/*     literature references: */

/*     T. J. Richmond, "Solvent Accessible Surface Area and */
/*     Excluded Volume in Proteins", Journal of Molecular Biology, */
/*     178, 63-89 (1984) */

/*     L. Wesson and D. Eisenberg, "Atomic Solvation Parameters */
/*     Applied to Molecular Dynamics of Proteins in Solution", */
/*     Protein Science, 1, 227-235 (1992) */

/*     variables and parameters: */

/*     total    total surface area of the whole structure */
/*     area     accessible surface area of each atom */
/*     darea    x,y,z components of the gradient of the area of */
/*                the molecule with respect to atomic coordinates */
/*     radius   radii of the individual atoms */
/*     weight   weight assigned to each atom's area; if set to */
/*                1.0, return is actual area in square Angstroms */
/*     probe    radius of the probe sphere */
/*     delta    tolerance used in the tests for sphere overlaps */
/*                and for colinearity */
/*     rmove    connectivity errors can usually be avoided if the */
/*                offending atom is shifted by this small amount */


/* Subroutine */ int surface1_(doublereal *total, doublereal *area, 
	doublereal *darea, doublereal *radius, doublereal *weight, doublereal 
	*probe)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 SURFACE1  --  Increase the Value of MAXA"
	    "RC\002)";
    static char fmt_70[] = "(/,\002 SURFACE1  --  Increase the Value\002,"
	    "\002 of MAXARC\002)";
    static char fmt_80[] = "(/,\002 SURFACE1  --  Increase the Value\002,"
	    "\002 of MAXARC\002)";
    static char fmt_140[] = "(/,\002 SURFACE1  --  Connectivity Error at A"
	    "tom\002,i6)";
    static char fmt_170[] = "(/,\002 SURFACE1  --  Negative Area at Atom\002"
	    ",i6)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);
    double asin(doublereal), acos(doublereal);
    integer do_fio(integer *, char *, ftnlen);
    double d_mod(doublereal *, doublereal *);

    /* Local variables */
    static doublereal b[1000];
    static integer i__, j, k, l, m;
    static doublereal p, r__[25000], s, t, v, b1[1000];
    static integer sign_yder__[1000];
    static doublereal t1, t2, cc, bg[1000];
    static integer ib, jb;
    static doublereal bk, dk, gi;
    static integer ii;
    static doublereal gl, gk;
    static integer mi, in, io, ni, ir;
    static doublereal td, tb, ti, tf;
    static integer lt[1000];
    static doublereal ri[1000], ex[1000], gr[1000], rr, xc[1000], tr, yc[1000]
	    , tt, zc[1000], xr, yr, tx, ty, tz, zr, ux[1000], uy[1000], uz[
	    1000], xc1[1000], yc1[1000], zc1[1000], tk1, tk2, tr2, bgl, dax, 
	    day, daz, the, rcn, rik, bsq[1000], eps;
    static integer key[1000];
    static doublereal dsq[1000], txb, tyb, axx, axy, axz, ayx, ayy, azx, azy, 
	    azz, uxl, uyl, wxl, uzl, txk, txr, tyr, tyk, tzk, txl, tyl, tzl;
    static logical top;
    static doublereal pid2, bsq1[1000], dsq1[1000], pix2, pix4, faca, facb, 
	    facc, rrx2, gaca, gacb, deal, decl, arcf[1000], arci[1000];
    static integer narc, ider[1000];
    static doublereal thec, ccsq, bsqk;
    static integer kent[1000];
    static doublereal bsql, dsql, ther[1000];
    static logical skip[25000], omit[1000];
    static doublereal wght, risq[1000];
    static integer kout[1000];
    static doublereal rrsq, xysq;
    extern /* Subroutine */ int sort2_(integer *, doublereal *, integer *), 
	    fatal_(void);
    static doublereal delta, dtkal, dtlal, dtkcl;
    static integer intag[1000];
    static doublereal exang, dtlcl;
    static logical moved;
    static doublereal therk;
    static logical komit;
    static doublereal rmove, risqk, risql, rplus, delta2, wxlsq;
    static integer intag1[1000];
    static doublereal arclen, cosine, arcsum, rminus;

    /* Fortran I/O blocks */
    static cilist io___166 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___229 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___241 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___271 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___272 = { 0, 0, 0, fmt_170, 0 };



#define darea_ref(a_1,a_2) darea[(a_2)*3 + a_1]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     zero the area and derivatives, and set the sphere radii */

    /* Parameter adjustments */
    --weight;
    --radius;
    darea -= 4;
    --area;

    /* Function Body */
    *total = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	area[i__] = 0.;
	darea_ref(1, i__) = 0.;
	darea_ref(2, i__) = 0.;
	darea_ref(3, i__) = 0.;
	r__[i__ - 1] = radius[i__];
	if (r__[i__ - 1] != 0.) {
	    r__[i__ - 1] += *probe;
	}
    }

/*     set pi multiples, overlap criterion and tolerances */

    pix2 = 6.2831853071795862;
    pix4 = 12.566370614359172;
    pid2 = 1.5707963267948966;
    delta = 1e-8;
/* Computing 2nd power */
    d__1 = delta;
    delta2 = d__1 * d__1;
    eps = 1e-8;
    rmove = 1e-8;
    for (i__ = 1; i__ <= 1000; ++i__) {
	ider[i__ - 1] = 0;
	sign_yder__[i__ - 1] = 0;
    }

/*     set the "skip" array to exclude all inactive atoms */
/*     that do not overlap any of the current active atoms */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	skip[i__ - 1] = TRUE_;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    xr = atoms_1.x[i__ - 1];
	    yr = atoms_1.y[i__ - 1];
	    zr = atoms_1.z__[i__ - 1];
	    rr = r__[i__ - 1];
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
		d__1 = rr + r__[k - 1];
		rplus = d__1 * d__1;
/* Computing 2nd power */
		d__1 = atoms_1.x[k - 1] - xr;
/* Computing 2nd power */
		d__2 = atoms_1.y[k - 1] - yr;
/* Computing 2nd power */
		d__3 = atoms_1.z__[k - 1] - zr;
		ccsq = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
		if (ccsq <= rplus) {
		    skip[k - 1] = FALSE_;
		}
	    }
	}
    }

/*     compute the area and derivatives of current "ir" sphere */

    i__1 = atoms_1.n;
    for (ir = 1; ir <= i__1; ++ir) {
	if (skip[ir - 1]) {
	    goto L180;
	}
	xr = atoms_1.x[ir - 1];
	yr = atoms_1.y[ir - 1];
	zr = atoms_1.z__[ir - 1];
	rr = r__[ir - 1];
	rrx2 = rr * 2.;
/* Computing 2nd power */
	d__1 = rr;
	rrsq = d__1 * d__1;
	wght = weight[ir];
	moved = FALSE_;

/*     initialize some counters and sums for the "ir" sphere */

L10:
	io = 0;
	jb = 0;
	ib = 0;
	arclen = 0.;
	exang = 0.;

/*     test each sphere to see if it overlaps the "ir" sphere */

	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (i__ == ir) {
		goto L30;
	    }
	    rplus = rr + r__[i__ - 1];
	    tx = atoms_1.x[i__ - 1] - xr;
	    if (abs(tx) >= rplus) {
		goto L30;
	    }
	    ty = atoms_1.y[i__ - 1] - yr;
	    if (abs(ty) >= rplus) {
		goto L30;
	    }
	    tz = atoms_1.z__[i__ - 1] - zr;
	    if (abs(tz) >= rplus) {
		goto L30;
	    }

/*     check for overlap of spheres by testing center to */
/*     center distance against sum and difference of radii */

/* Computing 2nd power */
	    d__1 = tx;
/* Computing 2nd power */
	    d__2 = ty;
	    xysq = d__1 * d__1 + d__2 * d__2;
	    if (xysq < delta2) {
		tx = delta;
		ty = 0.;
		xysq = delta2;
	    }
/* Computing 2nd power */
	    d__1 = tz;
	    ccsq = xysq + d__1 * d__1;
	    cc = sqrt(ccsq);
	    if (rplus - cc <= delta) {
		goto L30;
	    }
	    rminus = rr - r__[i__ - 1];

/*     check for a completely buried "ir" sphere */

	    if (cc - abs(rminus) <= delta) {
		if (rminus <= 0.) {
		    goto L180;
		}
		goto L30;
	    }

/*     calculate overlap parameters between "i" and "ir" sphere */

	    ++io;
	    xc1[io - 1] = tx;
	    yc1[io - 1] = ty;
	    zc1[io - 1] = tz;
	    dsq1[io - 1] = xysq;
	    bsq1[io - 1] = ccsq;
	    b1[io - 1] = cc;
	    gr[io - 1] = (ccsq + rplus * rminus) / (rrx2 * b1[io - 1]);
	    intag1[io - 1] = i__;
	    if (io > 1000) {
		io___166.ciunit = iounit_1.iout;
		s_wsfe(&io___166);
		e_wsfe();
		fatal_();
	    }
L30:
	    ;
	}

/*     case where no other spheres overlap the current sphere */

	if (io == 0) {
	    area[ir] = pix4;
	    goto L160;
	}

/*     case where only one sphere overlaps the current sphere */

	if (io == 1) {
	    k = 1;
	    txk = xc1[0];
	    tyk = yc1[0];
	    tzk = zc1[0];
	    bsqk = bsq1[0];
	    bk = b1[0];
	    intag[0] = intag1[0];
	    arcsum = pix2;
	    ++ib;
	    arclen += gr[k - 1] * arcsum;
	    if (! moved) {
		in = intag[k - 1];
/* Computing 2nd power */
		d__1 = r__[in - 1];
		t1 = arcsum * rrsq * (bsqk - rrsq + d__1 * d__1) / (rrx2 * 
			bsqk * bk);
		darea_ref(1, ir) = darea_ref(1, ir) - txk * t1 * wght;
		darea_ref(2, ir) = darea_ref(2, ir) - tyk * t1 * wght;
		darea_ref(3, ir) = darea_ref(3, ir) - tzk * t1 * wght;
		darea_ref(1, in) = darea_ref(1, in) + txk * t1 * wght;
		darea_ref(2, in) = darea_ref(2, in) + tyk * t1 * wght;
		darea_ref(3, in) = darea_ref(3, in) + tzk * t1 * wght;
	    }
	    goto L150;
	}

/*     general case where more than one sphere intersects the */
/*     current sphere; sort intersecting spheres by their degree */
/*     of overlap with the current main sphere */

	sort2_(&io, gr, key);
	i__2 = io;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    k = key[i__ - 1];
	    intag[i__ - 1] = intag1[k - 1];
	    xc[i__ - 1] = xc1[k - 1];
	    yc[i__ - 1] = yc1[k - 1];
	    zc[i__ - 1] = zc1[k - 1];
	    dsq[i__ - 1] = dsq1[k - 1];
	    b[i__ - 1] = b1[k - 1];
	    bsq[i__ - 1] = bsq1[k - 1];
	    omit[i__ - 1] = FALSE_;
	}

/*     radius of the each circle on the surface of the "ir" sphere */

	i__2 = io;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    gi = gr[i__ - 1] * rr;
	    bg[i__ - 1] = b[i__ - 1] * gi;
/* Computing 2nd power */
	    d__1 = gi;
	    risq[i__ - 1] = rrsq - d__1 * d__1;
	    ri[i__ - 1] = sqrt(risq[i__ - 1]);
/* Computing MIN */
/* Computing MAX */
	    d__3 = -1., d__4 = gr[i__ - 1];
	    d__1 = 1., d__2 = max(d__3,d__4);
	    ther[i__ - 1] = pid2 - asin((min(d__1,d__2)));
	}

/*     find boundary of inaccessible area on "ir" sphere */

	i__2 = io - 1;
	for (k = 1; k <= i__2; ++k) {
	    if (! omit[k - 1]) {
		txk = xc[k - 1];
		tyk = yc[k - 1];
		tzk = zc[k - 1];
		bk = b[k - 1];
		therk = ther[k - 1];
		i__3 = io;
		for (j = k + 1; j <= i__3; ++j) {
		    if (omit[j - 1]) {
			goto L60;
		    }

/*     check to see if J circle is intersecting K circle; */
/*     get distance between circle centers and sum of radii */

		    cc = (txk * xc[j - 1] + tyk * yc[j - 1] + tzk * zc[j - 1])
			     / (bk * b[j - 1]);
/* Computing MIN */
		    d__1 = 1., d__2 = max(-1.,cc);
		    cc = acos((min(d__1,d__2)));
		    td = therk + ther[j - 1];

/*     check to see if circles enclose separate regions */

		    if (cc >= td) {
			goto L60;
		    }

/*     check for circle J completely inside circle K */

		    if (cc + ther[j - 1] < therk) {
			goto L40;
		    }

/*     check for circles essentially parallel */

		    if (cc > delta) {
			goto L50;
		    }
L40:
		    omit[j - 1] = TRUE_;
		    goto L60;

/*     check for "ir" sphere completely buried */

L50:
		    if (pix2 - cc <= td) {
			goto L180;
		    }
L60:
		    ;
		}
	    }
	}

/*     find T value of circle intersections */

	i__2 = io;
	for (k = 1; k <= i__2; ++k) {
	    if (omit[k - 1]) {
		goto L110;
	    }
	    komit = omit[k - 1];
	    omit[k - 1] = TRUE_;
	    narc = 0;
	    top = FALSE_;
	    txk = xc[k - 1];
	    tyk = yc[k - 1];
	    tzk = zc[k - 1];
	    dk = sqrt(dsq[k - 1]);
	    bsqk = bsq[k - 1];
	    bk = b[k - 1];
	    gk = gr[k - 1] * rr;
	    risqk = risq[k - 1];
	    rik = ri[k - 1];
	    therk = ther[k - 1];

/*     rotation matrix elements */

	    t1 = tzk / (bk * dk);
	    axx = txk * t1;
	    axy = tyk * t1;
	    axz = dk / bk;
	    ayx = tyk / dk;
	    ayy = txk / dk;
	    azx = txk / bk;
	    azy = tyk / bk;
	    azz = tzk / bk;
	    i__3 = io;
	    for (l = 1; l <= i__3; ++l) {
		if (! omit[l - 1]) {
		    txl = xc[l - 1];
		    tyl = yc[l - 1];
		    tzl = zc[l - 1];

/*     rotate spheres so K vector colinear with z-axis */

		    uxl = txl * axx + tyl * axy - tzl * axz;
		    uyl = tyl * ayy - txl * ayx;
		    uzl = txl * azx + tyl * azy + tzl * azz;
/* Computing MIN */
/* Computing MAX */
		    d__3 = -1., d__4 = uzl / b[l - 1];
		    d__1 = 1., d__2 = max(d__3,d__4);
		    cosine = min(d__1,d__2);
		    if (acos(cosine) < therk + ther[l - 1]) {
/* Computing 2nd power */
			d__1 = uxl;
/* Computing 2nd power */
			d__2 = uyl;
			dsql = d__1 * d__1 + d__2 * d__2;
			tb = uzl * gk - bg[l - 1];
			txb = uxl * tb;
			tyb = uyl * tb;
			td = rik * dsql;
/* Computing 2nd power */
			d__1 = tb;
			tr2 = risqk * dsql - d__1 * d__1;
			tr2 = max(eps,tr2);
			tr = sqrt(tr2);
			txr = uxl * tr;
			tyr = uyl * tr;

/*     get T values of intersection for K circle */

			tb = (txb + tyr) / td;
/* Computing MIN */
			d__1 = 1., d__2 = max(-1.,tb);
			tb = min(d__1,d__2);
			tk1 = acos(tb);
			if (tyb - txr < 0.) {
			    tk1 = pix2 - tk1;
			}
			tb = (txb - tyr) / td;
/* Computing MIN */
			d__1 = 1., d__2 = max(-1.,tb);
			tb = min(d__1,d__2);
			tk2 = acos(tb);
			if (tyb + txr < 0.) {
			    tk2 = pix2 - tk2;
			}
			thec = (rrsq * uzl - gk * bg[l - 1]) / (rik * ri[l - 
				1] * b[l - 1]);
			if (abs(thec) < 1.) {
			    the = -acos(thec);
			} else if (thec >= 1.) {
			    the = 0.;
			} else if (thec <= -1.) {
			    the = -3.141592653589793238;
			}

/*     see if "tk1" is entry or exit point; check t=0 point; */
/*     "ti" is exit point, "tf" is entry point */

/* Computing MIN */
/* Computing MAX */
			d__3 = -1., d__4 = (uzl * gk - uxl * rik) / (b[l - 1] 
				* rr);
			d__1 = 1., d__2 = max(d__3,d__4);
			cosine = min(d__1,d__2);
			if ((acos(cosine) - ther[l - 1]) * (tk2 - tk1) <= 0.) 
				{
			    ti = tk2;
			    tf = tk1;
			} else {
			    ti = tk1;
			    tf = tk2;
			}
			++narc;
			if (narc >= 1000) {
			    io___229.ciunit = iounit_1.iout;
			    s_wsfe(&io___229);
			    e_wsfe();
			    fatal_();
			}
			if (tf <= ti) {
			    arcf[narc - 1] = tf;
			    arci[narc - 1] = 0.;
			    tf = pix2;
			    lt[narc - 1] = l;
			    ex[narc - 1] = the;
			    top = TRUE_;
			    ++narc;
			}
			arcf[narc - 1] = tf;
			arci[narc - 1] = ti;
			lt[narc - 1] = l;
			ex[narc - 1] = the;
			ux[l - 1] = uxl;
			uy[l - 1] = uyl;
			uz[l - 1] = uzl;
		    }
		}
	    }
	    omit[k - 1] = komit;

/*     special case; K circle without intersections */

	    if (narc <= 0) {
		goto L90;
	    }

/*     general case; sum up arclength and set connectivity code */

	    sort2_(&narc, arci, key);
	    arcsum = arci[0];
	    mi = key[0];
	    t = arcf[mi - 1];
	    ni = mi;
	    if (narc > 1) {
		i__3 = narc;
		for (j = 2; j <= i__3; ++j) {
		    m = key[j - 1];
		    if (t < arci[j - 1]) {
			arcsum = arcsum + arci[j - 1] - t;
			exang += ex[ni - 1];
			++jb;
			if (jb >= 1000) {
			    io___241.ciunit = iounit_1.iout;
			    s_wsfe(&io___241);
			    e_wsfe();
			    fatal_();
			}
			l = lt[ni - 1];
			++ider[l - 1];
			++sign_yder__[l - 1];
			kent[jb - 1] = l * 1000 + k;
			l = lt[m - 1];
			++ider[l - 1];
			--sign_yder__[l - 1];
			kout[jb - 1] = k * 1000 + l;
		    }
		    tt = arcf[m - 1];
		    if (tt >= t) {
			t = tt;
			ni = m;
		    }
		}
	    }
	    arcsum = arcsum + pix2 - t;
	    if (! top) {
		exang += ex[ni - 1];
		++jb;
		l = lt[ni - 1];
		++ider[l - 1];
		++sign_yder__[l - 1];
		kent[jb - 1] = l * 1000 + k;
		l = lt[mi - 1];
		++ider[l - 1];
		--sign_yder__[l - 1];
		kout[jb - 1] = k * 1000 + l;
	    }

/*     calculate the surface area derivatives */

	    i__3 = io;
	    for (l = 1; l <= i__3; ++l) {
		if (ider[l - 1] != 0) {
		    rcn = ider[l - 1] * rrsq;
		    ider[l - 1] = 0;
		    uzl = uz[l - 1];
		    gl = gr[l - 1] * rr;
		    bgl = bg[l - 1];
		    bsql = bsq[l - 1];
		    risql = risq[l - 1];
/* Computing 2nd power */
		    d__1 = uzl;
		    wxlsq = bsql - d__1 * d__1;
		    wxl = sqrt(wxlsq);
		    p = bgl - gk * uzl;
/* Computing 2nd power */
		    d__1 = p;
		    v = risqk * wxlsq - d__1 * d__1;
		    v = max(eps,v);
		    v = sqrt(v);
		    t1 = rr * (gk * (bgl - bsql) + uzl * (bgl - rrsq)) / (v * 
			    risql * bsql);
		    deal = -wxl * t1;
		    decl = -uzl * t1 - rr / v;
		    dtkal = (wxlsq - p) / (wxl * v);
		    dtkcl = (uzl - gk) / v;
		    s = gk * b[l - 1] - gl * uzl;
		    t1 = gk * 2. - uzl;
		    t2 = rrsq - bgl;
		    dtlal = -(risql * wxlsq * b[l - 1] * t1 - s * (wxlsq * t2 
			    + risql * bsql)) / (risql * wxl * bsql * v);
		    dtlcl = -(risql * b[l - 1] * (uzl * t1 - bgl) - uzl * t2 *
			     s) / (risql * bsql * v);
		    gaca = rcn * (deal - (gk * dtkal - gl * dtlal) / rr) / 
			    wxl;
		    gacb = (gk - uzl * gl / b[l - 1]) * sign_yder__[l - 1] * 
			    rr / wxlsq;
		    sign_yder__[l - 1] = 0;
		    if (! moved) {
			faca = ux[l - 1] * gaca - uy[l - 1] * gacb;
			facb = uy[l - 1] * gaca + ux[l - 1] * gacb;
			facc = rcn * (decl - (gk * dtkcl - gl * dtlcl) / rr);
			dax = axx * faca - ayx * facb + azx * facc;
			day = axy * faca + ayy * facb + azy * facc;
			daz = azz * facc - axz * faca;
			in = intag[l - 1];
			darea_ref(1, ir) = darea_ref(1, ir) + dax * wght;
			darea_ref(2, ir) = darea_ref(2, ir) + day * wght;
			darea_ref(3, ir) = darea_ref(3, ir) + daz * wght;
			darea_ref(1, in) = darea_ref(1, in) - dax * wght;
			darea_ref(2, in) = darea_ref(2, in) - day * wght;
			darea_ref(3, in) = darea_ref(3, in) - daz * wght;
		    }
		}
	    }
	    goto L100;
L90:
	    arcsum = pix2;
	    ++ib;
L100:
	    arclen += gr[k - 1] * arcsum;
	    if (! moved) {
		in = intag[k - 1];
/* Computing 2nd power */
		d__1 = r__[in - 1];
		t1 = arcsum * rrsq * (bsqk - rrsq + d__1 * d__1) / (rrx2 * 
			bsqk * bk);
		darea_ref(1, ir) = darea_ref(1, ir) - txk * t1 * wght;
		darea_ref(2, ir) = darea_ref(2, ir) - tyk * t1 * wght;
		darea_ref(3, ir) = darea_ref(3, ir) - tzk * t1 * wght;
		darea_ref(1, in) = darea_ref(1, in) + txk * t1 * wght;
		darea_ref(2, in) = darea_ref(2, in) + tyk * t1 * wght;
		darea_ref(3, in) = darea_ref(3, in) + tzk * t1 * wght;
	    }
L110:
	    ;
	}
	if (arclen == 0.) {
	    goto L180;
	}
	if (jb == 0) {
	    goto L150;
	}

/*     find number of independent boundaries and check connectivity */

	j = 0;
	i__2 = jb;
	for (k = 1; k <= i__2; ++k) {
	    if (kout[k - 1] != 0) {
		i__ = k;
L120:
		m = kout[i__ - 1];
		kout[i__ - 1] = 0;
		++j;
		i__3 = jb;
		for (ii = 1; ii <= i__3; ++ii) {
		    if (m == kent[ii - 1]) {
			if (ii == k) {
			    ++ib;
			    if (j == jb) {
				goto L150;
			    }
			    goto L130;
			}
			i__ = ii;
			goto L120;
		    }
		}
L130:
		;
	    }
	}
	++ib;

/*     attempt to fix connectivity error by moving atom slightly */

	if (moved) {
	    io___271.ciunit = iounit_1.iout;
	    s_wsfe(&io___271);
	    do_fio(&c__1, (char *)&ir, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    moved = TRUE_;
	    xr += rmove;
	    yr += rmove;
	    zr += rmove;
	    goto L10;
	}

/*     form the accessible area for the current atom */

L150:
	area[ir] = ib * pix2 + exang + arclen;
	area[ir] = d_mod(&area[ir], &pix4);
L160:
	area[ir] *= rrsq;

/*     attempt to fix negative area by moving atom slightly */

	if (area[ir] < 0.) {
	    if (moved) {
		io___272.ciunit = iounit_1.iout;
		s_wsfe(&io___272);
		do_fio(&c__1, (char *)&ir, (ftnlen)sizeof(integer));
		e_wsfe();
	    } else {
		moved = TRUE_;
		xr += rmove;
		yr += rmove;
		zr += rmove;
		goto L10;
	    }
	}

/*     weight the accessible area by the scale factor */

	area[ir] *= wght;
	*total += area[ir];
L180:
	;
    }

/*     zero out the area derivatives for the inactive atoms */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (! usage_1.use[i__ - 1]) {
	    darea_ref(1, i__) = 0.;
	    darea_ref(2, i__) = 0.;
	    darea_ref(3, i__) = 0.;
	}
    }

/*     print out the surface area and derivatives for each atom */

/*     if (debug) then */
/*        write (iout,190) */
/* 190    format (/,' Weighted Atomic Surface Areas and Derivatives :', */
/*    &           //,4x,'Atom',7x,'Area Term',10x,'dA/dx', */
/*    &              7x,'dA/dy',7x,'dA/dz',6x,'Weight',/) */
/*        do i = 1, n */
/*           if (.not. skip(i)) then */
/*              write (iout,200)  i,area(i),(darea(j,i),j=1,3),weight(i) */
/* 200          format (i8,4x,f12.4,3x,3f12.4,f12.4) */
/*           end if */
/*        end do */
/*        write (iout,210)  total */
/* 210    format (/,' Total Weighted Surface Area :',5x,f16.4) */
/*     end if */
    return 0;
} /* surface1_ */

#undef darea_ref


