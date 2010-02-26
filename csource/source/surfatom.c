/* surfatom.f -- translated by f2c (version 20050501).
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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

/* Table of constant values */

static doublereal c_b7 = 12.566370614359172;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine surfatom  --  exposed surface area of an atom  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "surfatom" performs an analytical computation of the surface */
/*     area of a specified atom; a simplified version of "surface" */

/*     literature references: */

/*     T. J. Richmond, "Solvent Accessible Surface Area and */
/*     Excluded Volume in Proteins", Journal of Molecular Biology, */
/*     178, 63-89 (1984) */

/*     L. Wesson and D. Eisenberg, "Atomic Solvation Parameters */
/*     Applied to Molecular Dynamics of Proteins in Solution", */
/*     Protein Science, 1, 227-235 (1992) */

/*     variables and parameters: */

/*     ir       number of atom for which area is desired */
/*     area     accessible surface area of the atom */
/*     radius   radii of each of the individual atoms */


/* Subroutine */ int surfatom_(integer *ir, doublereal *area, doublereal *
	radius)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 SURFATOM  --  Increase the Value of MAXA"
	    "RC\002)";
    static char fmt_70[] = "(/,\002 SURFATOM  --  Increase the Value\002,"
	    "\002 of MAXARC\002)";
    static char fmt_80[] = "(/,\002 SURFATOM  --  Increase the Value\002,"
	    "\002 of MAXARC\002)";
    static char fmt_140[] = "(/,\002 SURFATOM  --  Connectivity Error at A"
	    "tom\002,i6)";
    static char fmt_160[] = "(/,\002 SURFATOM  --  Negative Area at Atom\002"
	    ",i6)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);
    double d_mod(doublereal *, doublereal *), asin(doublereal), acos(
	    doublereal);
    integer do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal b[1000];
    static integer i__, j, k, m;
    static doublereal t, b1[1000], t1, cc, bg[1000];
    static integer ib, jb;
    static doublereal bk, dk, gi;
    static integer ii;
    static doublereal gk;
    static integer mi, ni, io;
    static doublereal td, tb, ti, tf, ri[1000];
    static integer lt[1000];
    static doublereal ex[1000], gr[1000], xc[1000], rr, yc[1000], tr, zc[1000]
	    , tt, xr, yr, zr, tx, ty, tz, ux[1000], uy[1000], uz[1000], xc1[
	    1000], yc1[1000], zc1[1000], tk1, tk2, tr2, the, rik, bsq[1000], 
	    eps;
    static integer key[1000];
    static doublereal dsq[1000], txb, tyb, axx, axy, axz, ayx, ayy, azx, azy, 
	    azz, uxj, uyj, uzj, txk, tyk, txr, tyr, tzk, txj, tyj, tzj;
    static logical top;
    static doublereal bsq1[1000], dsq1[1000], pix2, arcf[1000], arci[1000];
    static integer narc;
    static doublereal thec, ccsq, bsqk;
    static integer kent[1000];
    static doublereal dsqj, ther[1000];
    static logical omit[1000];
    static doublereal risq[1000];
    static integer kout[1000];
    static doublereal rrsq, xysq;
    extern /* Subroutine */ int sort2_(integer *, doublereal *, integer *), 
	    fatal_(void);
    static doublereal delta;
    static integer intag[1000];
    static doublereal exang;
    static logical moved;
    static doublereal therk, rmove, risqk, rplus, delta2;
    static integer intag1[1000];
    static doublereal arclen, cosine, arcsum, rminus;

    /* Fortran I/O blocks */
    static cilist io___26 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___94 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___107 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___112 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___113 = { 0, 0, 0, fmt_160, 0 };




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




/*     zero out the surface area for the sphere of interest */

    /* Parameter adjustments */
    --radius;

    /* Function Body */
    *area = 0.;
    if (radius[*ir] == 0.) {
	return 0;
    }

/*     set the overlap significance and connectivity shift */

    pix2 = 6.2831853071795862;
    delta = 1e-8;
    delta2 = delta * delta;
    eps = 1e-8;
    moved = FALSE_;
    rmove = 1e-8;

/*     store coordinates and radius of the sphere of interest */

    xr = atoms_1.x[*ir - 1];
    yr = atoms_1.y[*ir - 1];
    zr = atoms_1.z__[*ir - 1];
    rr = radius[*ir];
    rrsq = rr * rr;

/*     initialize values of some counters and summations */

L10:
    io = 0;
    jb = 0;
    ib = 0;
    arclen = 0.;
    exang = 0.;

/*     test each sphere to see if it overlaps the sphere of interest */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == *ir || radius[i__] == 0.) {
	    goto L30;
	}
	rplus = rr + radius[i__];
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

/*     check for sphere overlap by testing distance against radii */

	xysq = tx * tx + ty * ty;
	if (xysq < delta2) {
	    tx = delta;
	    ty = 0.;
	    xysq = delta2;
	}
	ccsq = xysq + tz * tz;
	cc = sqrt(ccsq);
	if (rplus - cc <= delta) {
	    goto L30;
	}
	rminus = rr - radius[i__];

/*     check to see if sphere of interest is completely buried */

	if (cc - abs(rminus) <= delta) {
	    if (rminus <= 0.) {
		goto L170;
	    }
	    goto L30;
	}

/*     check for too many overlaps with sphere of interest */

	if (io >= 1000) {
	    io___26.ciunit = iounit_1.iout;
	    s_wsfe(&io___26);
	    e_wsfe();
	    fatal_();
	}

/*     get overlap between current sphere and sphere of interest */

	++io;
	xc1[io - 1] = tx;
	yc1[io - 1] = ty;
	zc1[io - 1] = tz;
	dsq1[io - 1] = xysq;
	bsq1[io - 1] = ccsq;
	b1[io - 1] = cc;
	gr[io - 1] = (ccsq + rplus * rminus) / (rr * 2. * b1[io - 1]);
	intag1[io - 1] = i__;
	omit[io - 1] = FALSE_;
L30:
	;
    }

/*     case where no other spheres overlap the sphere of interest */

    if (io == 0) {
	*area = rrsq * 12.566370614359172;
	return 0;
    }

/*     case where only one sphere overlaps the sphere of interest */

    if (io == 1) {
	*area = pix2 * (gr[0] + 1.);
	*area = d_mod(area, &c_b7) * rrsq;
	return 0;
    }

/*     case where many spheres intersect the sphere of interest; */
/*     sort the intersecting spheres by their degree of overlap */

    sort2_(&io, gr, key);
    i__1 = io;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = key[i__ - 1];
	intag[i__ - 1] = intag1[k - 1];
	xc[i__ - 1] = xc1[k - 1];
	yc[i__ - 1] = yc1[k - 1];
	zc[i__ - 1] = zc1[k - 1];
	dsq[i__ - 1] = dsq1[k - 1];
	b[i__ - 1] = b1[k - 1];
	bsq[i__ - 1] = bsq1[k - 1];
    }

/*     get radius of each overlap circle on surface of the sphere */

    i__1 = io;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gi = gr[i__ - 1] * rr;
	bg[i__ - 1] = b[i__ - 1] * gi;
	risq[i__ - 1] = rrsq - gi * gi;
	ri[i__ - 1] = sqrt(risq[i__ - 1]);
/* Computing MIN */
/* Computing MAX */
	d__3 = -1., d__4 = gr[i__ - 1];
	d__1 = 1., d__2 = max(d__3,d__4);
	ther[i__ - 1] = 1.5707963267948966 - asin((min(d__1,d__2)));
    }

/*     find boundary of inaccessible area on sphere of interest */

    i__1 = io - 1;
    for (k = 1; k <= i__1; ++k) {
	if (! omit[k - 1]) {
	    txk = xc[k - 1];
	    tyk = yc[k - 1];
	    tzk = zc[k - 1];
	    bk = b[k - 1];
	    therk = ther[k - 1];

/*     check to see if J circle is intersecting K circle; */
/*     get distance between circle centers and sum of radii */

	    i__2 = io;
	    for (j = k + 1; j <= i__2; ++j) {
		if (omit[j - 1]) {
		    goto L60;
		}
		cc = (txk * xc[j - 1] + tyk * yc[j - 1] + tzk * zc[j - 1]) / (
			bk * b[j - 1]);
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

/*     check for circles that are essentially parallel */

		if (cc > delta) {
		    goto L50;
		}
L40:
		omit[j - 1] = TRUE_;
		goto L60;

/*     check to see if sphere of interest is completely buried */

L50:
		if (pix2 - cc <= td) {
		    goto L170;
		}
L60:
		;
	    }
	}
    }

/*     find T value of circle intersections */

    i__1 = io;
    for (k = 1; k <= i__1; ++k) {
	if (omit[k - 1]) {
	    goto L110;
	}
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
	i__2 = io;
	for (j = 1; j <= i__2; ++j) {
	    if (! omit[j - 1]) {
		txj = xc[j - 1];
		tyj = yc[j - 1];
		tzj = zc[j - 1];

/*     rotate spheres so K vector colinear with z-axis */

		uxj = txj * axx + tyj * axy - tzj * axz;
		uyj = tyj * ayy - txj * ayx;
		uzj = txj * azx + tyj * azy + tzj * azz;
/* Computing MIN */
/* Computing MAX */
		d__3 = -1., d__4 = uzj / b[j - 1];
		d__1 = 1., d__2 = max(d__3,d__4);
		cosine = min(d__1,d__2);
		if (acos(cosine) < therk + ther[j - 1]) {
		    dsqj = uxj * uxj + uyj * uyj;
		    tb = uzj * gk - bg[j - 1];
		    txb = uxj * tb;
		    tyb = uyj * tb;
		    td = rik * dsqj;
		    tr2 = risqk * dsqj - tb * tb;
		    tr2 = max(eps,tr2);
		    tr = sqrt(tr2);
		    txr = uxj * tr;
		    tyr = uyj * tr;

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
		    thec = (rrsq * uzj - gk * bg[j - 1]) / (rik * ri[j - 1] * 
			    b[j - 1]);
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
		    d__3 = -1., d__4 = (uzj * gk - uxj * rik) / (b[j - 1] * 
			    rr);
		    d__1 = 1., d__2 = max(d__3,d__4);
		    cosine = min(d__1,d__2);
		    if ((acos(cosine) - ther[j - 1]) * (tk2 - tk1) <= 0.) {
			ti = tk2;
			tf = tk1;
		    } else {
			ti = tk2;
			tf = tk1;
		    }
		    ++narc;
		    if (narc >= 1000) {
			io___94.ciunit = iounit_1.iout;
			s_wsfe(&io___94);
			e_wsfe();
			fatal_();
		    }
		    if (tf <= ti) {
			arcf[narc - 1] = tf;
			arci[narc - 1] = 0.;
			tf = pix2;
			lt[narc - 1] = j;
			ex[narc - 1] = the;
			top = TRUE_;
			++narc;
		    }
		    arcf[narc - 1] = tf;
		    arci[narc - 1] = ti;
		    lt[narc - 1] = j;
		    ex[narc - 1] = the;
		    ux[j - 1] = uxj;
		    uy[j - 1] = uyj;
		    uz[j - 1] = uzj;
		}
	    }
	}
	omit[k - 1] = FALSE_;

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
	    i__2 = narc;
	    for (j = 2; j <= i__2; ++j) {
		m = key[j - 1];
		if (t < arci[j - 1]) {
		    arcsum = arcsum + arci[j - 1] - t;
		    exang += ex[ni - 1];
		    ++jb;
		    if (jb >= 1000) {
			io___107.ciunit = iounit_1.iout;
			s_wsfe(&io___107);
			e_wsfe();
			fatal_();
		    }
		    i__ = lt[ni - 1];
		    kent[jb - 1] = i__ * 1000 + k;
		    i__ = lt[m - 1];
		    kout[jb - 1] = k * 1000 + i__;
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
	    i__ = lt[ni - 1];
	    kent[jb - 1] = i__ * 1000 + k;
	    i__ = lt[mi - 1];
	    kout[jb - 1] = k * 1000 + i__;
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
	goto L170;
    }
    if (jb == 0) {
	goto L150;
    }

/*     find number of independent boundaries and check connectivity */

    j = 0;
    i__1 = jb;
    for (k = 1; k <= i__1; ++k) {
	if (kout[k - 1] != 0) {
	    i__ = k;
L120:
	    m = kout[i__ - 1];
	    kout[i__ - 1] = 0;
	    ++j;
	    i__2 = jb;
	    for (ii = 1; ii <= i__2; ++ii) {
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
	io___112.ciunit = iounit_1.iout;
	s_wsfe(&io___112);
	do_fio(&c__1, (char *)&(*ir), (ftnlen)sizeof(integer));
	e_wsfe();
    } else {
	moved = TRUE_;
	xr += rmove;
	yr += rmove;
	zr += rmove;
	goto L10;
    }

/*     compute the exposed surface area for the sphere of interest */

L150:
    *area = ib * pix2 + exang + arclen;
    *area = d_mod(area, &c_b7) * rrsq;

/*     attempt to fix negative area by moving atom slightly */

    if (*area < 0.) {
	if (moved) {
	    io___113.ciunit = iounit_1.iout;
	    s_wsfe(&io___113);
	    do_fio(&c__1, (char *)&(*ir), (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    moved = TRUE_;
	    xr += rmove;
	    yr += rmove;
	    zr += rmove;
	    goto L10;
	}
    }
L170:
    return 0;
} /* surfatom_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine surfatom1  --  surface area and derivs of atom  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "surfatom1" performs an analytical computation of the surface */
/*     area and first derivatives with respect to Cartesian coordinates */
/*     of a specified atom */


/* Subroutine */ int surfatom1_(integer *ir, doublereal *area, doublereal *
	darea, doublereal *radius)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 SURFATOM  --  Increase the Value of MAXA"
	    "RC\002)";
    static char fmt_70[] = "(/,\002 SURFATOM  --  Increase the Value\002,"
	    "\002 of MAXARC\002)";
    static char fmt_80[] = "(/,\002 SURFATOM  --  Increase the Value\002,"
	    "\002 of MAXARC\002)";
    static char fmt_140[] = "(/,\002 SURFATOM  --  Connectivity Error at A"
	    "tom\002,i6)";
    static char fmt_160[] = "(/,\002 SURFATOM  --  Negative Area at Atom\002"
	    ",i6)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);
    double d_mod(doublereal *, doublereal *), asin(doublereal), acos(
	    doublereal);
    integer do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal b[1000];
    static integer i__, j, k, m;
    static doublereal p, s, t, v, b1[1000];
    static integer sign_yder__[1000];
    static doublereal t1, t2, cc, bg[1000];
    static integer ib, jb;
    static doublereal bk, dk, gi;
    static integer ii;
    static doublereal gl, gk;
    static integer mi, in, io, ni;
    static doublereal td, tb, ti, tf;
    static integer lt[1000];
    static doublereal ri[1000], ex[1000], gr[1000], rr, xc[1000], tr, yc[1000]
	    , tt, zc[1000], xr, yr, zr, tx, ty, tz, ux[1000], uy[1000], uz[
	    1000], xc1[1000], yc1[1000], zc1[1000], tk1, tk2, tr2, bgl, dax, 
	    day, daz, the, rcn, rik, bsq[1000], eps;
    static integer key[1000];
    static doublereal rin, dsq[1000], txb, tyb, axx, axy, axz, ayx, ayy, azx, 
	    azy, azz, uxj, uyj, wxl, uzj, txk, txr, tyr, tyk, tzk, txj, tyj, 
	    tzj, uzl;
    static logical top;
    static doublereal bsq1[1000], dsq1[1000], pix2, faca, facb, facc, rrx2, 
	    gaca, gacb, deal, decl, arcf[1000], arci[1000];
    static integer narc, ider[1000];
    static doublereal thec, ccsq, bsqk;
    static integer kent[1000];
    static doublereal bsql, dsqj, ther[1000];
    static logical omit[1000];
    static doublereal risq[1000];
    static integer kout[1000];
    static doublereal rrsq, xysq;
    extern /* Subroutine */ int sort2_(integer *, doublereal *, integer *), 
	    fatal_(void);
    static doublereal delta, dtkal, dtlal, dtkcl;
    static integer intag[1000];
    static doublereal exang, dtlcl;
    static logical moved;
    static doublereal therk, rmove, risqk, risql, rplus, delta2, wxlsq;
    static integer intag1[1000];
    static doublereal arclen, cosine, arcsum, rminus;

    /* Fortran I/O blocks */
    static cilist io___142 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___213 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___225 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___256 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___257 = { 0, 0, 0, fmt_160, 0 };



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




/*     zero out the area and derivatives for sphere of interest */

    /* Parameter adjustments */
    --radius;
    darea -= 4;

    /* Function Body */
    *area = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	darea_ref(1, i__) = 0.;
	darea_ref(2, i__) = 0.;
	darea_ref(3, i__) = 0.;
    }
    if (radius[*ir] == 0.) {
	return 0;
    }

/*     set the overlap significance and connectivity shift */

    pix2 = 6.2831853071795862;
    delta = 1e-8;
    delta2 = delta * delta;
    eps = 1e-8;
    moved = FALSE_;
    rmove = 1e-8;
    for (i__ = 1; i__ <= 1000; ++i__) {
	ider[i__ - 1] = 0;
	sign_yder__[i__ - 1] = 0;
    }

/*     store coordinates and radius of the sphere of interest */

    xr = atoms_1.x[*ir - 1];
    yr = atoms_1.y[*ir - 1];
    zr = atoms_1.z__[*ir - 1];
    rr = radius[*ir];
    rrx2 = rr * 2.;
    rrsq = rr * rr;

/*     initialize values of some counters and summations */

L10:
    io = 0;
    jb = 0;
    ib = 0;
    arclen = 0.;
    exang = 0.;

/*     test each sphere to see if it overlaps the sphere of interest */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == *ir || radius[i__] == 0.) {
	    goto L30;
	}
	rplus = rr + radius[i__];
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

/*     check for sphere overlap by testing distance against radii */

	xysq = tx * tx + ty * ty;
	if (xysq < delta2) {
	    tx = delta;
	    ty = 0.;
	    xysq = delta2;
	}
	ccsq = xysq + tz * tz;
	cc = sqrt(ccsq);
	if (rplus - cc <= delta) {
	    goto L30;
	}
	rminus = rr - radius[i__];

/*     check to see if sphere of interest is completely buried */

	if (cc - abs(rminus) <= delta) {
	    if (rminus <= 0.) {
		goto L170;
	    }
	    goto L30;
	}

/*     check for too many overlaps with sphere of interest */

	if (io >= 1000) {
	    io___142.ciunit = iounit_1.iout;
	    s_wsfe(&io___142);
	    e_wsfe();
	    fatal_();
	}

/*     get overlap between current sphere and sphere of interest */

	++io;
	xc1[io - 1] = tx;
	yc1[io - 1] = ty;
	zc1[io - 1] = tz;
	dsq1[io - 1] = xysq;
	bsq1[io - 1] = ccsq;
	b1[io - 1] = cc;
	gr[io - 1] = (ccsq + rplus * rminus) / (rr * 2. * b1[io - 1]);
	intag1[io - 1] = i__;
	omit[io - 1] = FALSE_;
L30:
	;
    }

/*     case where no other spheres overlap the sphere of interest */

    if (io == 0) {
	*area = rrsq * 12.566370614359172;
	return 0;
    }

/*     case where only one sphere overlaps the sphere of interest */

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
	    rin = radius[in];
	    t1 = arcsum * rrsq * (bsqk - rrsq + rin * rin) / (rrx2 * bsqk * 
		    bk);
	    darea_ref(1, *ir) = darea_ref(1, *ir) - txk * t1;
	    darea_ref(2, *ir) = darea_ref(2, *ir) - tyk * t1;
	    darea_ref(3, *ir) = darea_ref(3, *ir) - tzk * t1;
	    darea_ref(1, in) = darea_ref(1, in) + txk * t1;
	    darea_ref(2, in) = darea_ref(2, in) + tyk * t1;
	    darea_ref(3, in) = darea_ref(3, in) + tzk * t1;
	}
	*area = pix2 * (gr[0] + 1.);
	*area = d_mod(area, &c_b7) * rrsq;
	return 0;
    }

/*     case where many spheres intersect the sphere of interest; */
/*     sort the intersecting spheres by their degree of overlap */

    sort2_(&io, gr, key);
    i__1 = io;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = key[i__ - 1];
	intag[i__ - 1] = intag1[k - 1];
	xc[i__ - 1] = xc1[k - 1];
	yc[i__ - 1] = yc1[k - 1];
	zc[i__ - 1] = zc1[k - 1];
	dsq[i__ - 1] = dsq1[k - 1];
	b[i__ - 1] = b1[k - 1];
	bsq[i__ - 1] = bsq1[k - 1];
    }

/*     get radius of each overlap circle on surface of the sphere */

    i__1 = io;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gi = gr[i__ - 1] * rr;
	bg[i__ - 1] = b[i__ - 1] * gi;
	risq[i__ - 1] = rrsq - gi * gi;
	ri[i__ - 1] = sqrt(risq[i__ - 1]);
/* Computing MIN */
/* Computing MAX */
	d__3 = -1., d__4 = gr[i__ - 1];
	d__1 = 1., d__2 = max(d__3,d__4);
	ther[i__ - 1] = 1.5707963267948966 - asin((min(d__1,d__2)));
    }

/*     find boundary of inaccessible area on sphere of interest */

    i__1 = io - 1;
    for (k = 1; k <= i__1; ++k) {
	if (! omit[k - 1]) {
	    txk = xc[k - 1];
	    tyk = yc[k - 1];
	    tzk = zc[k - 1];
	    bk = b[k - 1];
	    therk = ther[k - 1];

/*     check to see if J circle is intersecting K circle; */
/*     get distance between circle centers and sum of radii */

	    i__2 = io;
	    for (j = k + 1; j <= i__2; ++j) {
		if (omit[j - 1]) {
		    goto L60;
		}
		cc = (txk * xc[j - 1] + tyk * yc[j - 1] + tzk * zc[j - 1]) / (
			bk * b[j - 1]);
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

/*     check for circles that are essentially parallel */

		if (cc > delta) {
		    goto L50;
		}
L40:
		omit[j - 1] = TRUE_;
		goto L60;

/*     check to see if sphere of interest is completely buried */

L50:
		if (pix2 - cc <= td) {
		    goto L170;
		}
L60:
		;
	    }
	}
    }

/*     find T value of circle intersections */

    i__1 = io;
    for (k = 1; k <= i__1; ++k) {
	if (omit[k - 1]) {
	    goto L110;
	}
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
	i__2 = io;
	for (j = 1; j <= i__2; ++j) {
	    if (! omit[j - 1]) {
		txj = xc[j - 1];
		tyj = yc[j - 1];
		tzj = zc[j - 1];

/*     rotate spheres so K vector colinear with z-axis */

		uxj = txj * axx + tyj * axy - tzj * axz;
		uyj = tyj * ayy - txj * ayx;
		uzj = txj * azx + tyj * azy + tzj * azz;
/* Computing MIN */
/* Computing MAX */
		d__3 = -1., d__4 = uzj / b[j - 1];
		d__1 = 1., d__2 = max(d__3,d__4);
		cosine = min(d__1,d__2);
		if (acos(cosine) < therk + ther[j - 1]) {
		    dsqj = uxj * uxj + uyj * uyj;
		    tb = uzj * gk - bg[j - 1];
		    txb = uxj * tb;
		    tyb = uyj * tb;
		    td = rik * dsqj;
		    tr2 = risqk * dsqj - tb * tb;
		    tr2 = max(eps,tr2);
		    tr = sqrt(tr2);
		    txr = uxj * tr;
		    tyr = uyj * tr;

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
		    thec = (rrsq * uzj - gk * bg[j - 1]) / (rik * ri[j - 1] * 
			    b[j - 1]);
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
		    d__3 = -1., d__4 = (uzj * gk - uxj * rik) / (b[j - 1] * 
			    rr);
		    d__1 = 1., d__2 = max(d__3,d__4);
		    cosine = min(d__1,d__2);
		    if ((acos(cosine) - ther[j - 1]) * (tk2 - tk1) <= 0.) {
			ti = tk2;
			tf = tk1;
		    } else {
			ti = tk2;
			tf = tk1;
		    }
		    ++narc;
		    if (narc >= 1000) {
			io___213.ciunit = iounit_1.iout;
			s_wsfe(&io___213);
			e_wsfe();
			fatal_();
		    }
		    if (tf <= ti) {
			arcf[narc - 1] = tf;
			arci[narc - 1] = 0.;
			tf = pix2;
			lt[narc - 1] = j;
			ex[narc - 1] = the;
			top = TRUE_;
			++narc;
		    }
		    arcf[narc - 1] = tf;
		    arci[narc - 1] = ti;
		    lt[narc - 1] = j;
		    ex[narc - 1] = the;
		    ux[j - 1] = uxj;
		    uy[j - 1] = uyj;
		    uz[j - 1] = uzj;
		}
	    }
	}
	omit[k - 1] = FALSE_;

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
	    i__2 = narc;
	    for (j = 2; j <= i__2; ++j) {
		m = key[j - 1];
		if (t < arci[j - 1]) {
		    arcsum = arcsum + arci[j - 1] - t;
		    exang += ex[ni - 1];
		    ++jb;
		    if (jb >= 1000) {
			io___225.ciunit = iounit_1.iout;
			s_wsfe(&io___225);
			e_wsfe();
			fatal_();
		    }
		    i__ = lt[ni - 1];
		    ++ider[i__ - 1];
		    ++sign_yder__[i__ - 1];
		    kent[jb - 1] = i__ * 1000 + k;
		    i__ = lt[m - 1];
		    ++ider[i__ - 1];
		    --sign_yder__[i__ - 1];
		    kout[jb - 1] = k * 1000 + i__;
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
	    i__ = lt[ni - 1];
	    ++ider[i__ - 1];
	    ++sign_yder__[i__ - 1];
	    kent[jb - 1] = i__ * 1000 + k;
	    i__ = lt[mi - 1];
	    ++ider[i__ - 1];
	    --sign_yder__[i__ - 1];
	    kout[jb - 1] = k * 1000 + i__;
	}

/*     calculate the surface area derivatives */

	i__2 = io;
	for (j = 1; j <= i__2; ++j) {
	    if (ider[j - 1] != 0) {
		rcn = ider[j - 1] * rrsq;
		ider[j - 1] = 0;
		uzl = uz[j - 1];
		gl = gr[j - 1] * rr;
		bgl = bg[j - 1];
		bsql = bsq[j - 1];
		risql = risq[j - 1];
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
		s = gk * b[j - 1] - gl * uzl;
		t1 = gk * 2. - uzl;
		t2 = rrsq - bgl;
		dtlal = -(risql * wxlsq * b[j - 1] * t1 - s * (wxlsq * t2 + 
			risql * bsql)) / (risql * wxl * bsql * v);
		dtlcl = -(risql * b[j - 1] * (uzl * t1 - bgl) - uzl * t2 * s) 
			/ (risql * bsql * v);
		gaca = rcn * (deal - (gk * dtkal - gl * dtlal) / rr) / wxl;
		gacb = (gk - uzl * gl / b[j - 1]) * sign_yder__[j - 1] * rr / 
			wxlsq;
		sign_yder__[j - 1] = 0;
		if (! moved) {
		    faca = ux[j - 1] * gaca - uy[j - 1] * gacb;
		    facb = uy[j - 1] * gaca + ux[j - 1] * gacb;
		    facc = rcn * (decl - (gk * dtkcl - gl * dtlcl) / rr);
		    dax = axx * faca - ayx * facb + azx * facc;
		    day = axy * faca + ayy * facb + azy * facc;
		    daz = azz * facc - axz * faca;
		    in = intag[j - 1];
		    darea_ref(1, *ir) = darea_ref(1, *ir) + dax;
		    darea_ref(2, *ir) = darea_ref(2, *ir) + day;
		    darea_ref(3, *ir) = darea_ref(3, *ir) + daz;
		    darea_ref(1, in) = darea_ref(1, in) - dax;
		    darea_ref(2, in) = darea_ref(2, in) - day;
		    darea_ref(3, in) = darea_ref(3, in) - daz;
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
	    rin = radius[in];
	    t1 = arcsum * rrsq * (bsqk - rrsq + rin * rin) / (rrx2 * bsqk * 
		    bk);
	    darea_ref(1, *ir) = darea_ref(1, *ir) - txk * t1;
	    darea_ref(2, *ir) = darea_ref(2, *ir) - tyk * t1;
	    darea_ref(3, *ir) = darea_ref(3, *ir) - tzk * t1;
	    darea_ref(1, in) = darea_ref(1, in) + txk * t1;
	    darea_ref(2, in) = darea_ref(2, in) + tyk * t1;
	    darea_ref(3, in) = darea_ref(3, in) + tzk * t1;
	}
L110:
	;
    }
    if (arclen == 0.) {
	goto L170;
    }
    if (jb == 0) {
	goto L150;
    }

/*     find number of independent boundaries and check connectivity */

    j = 0;
    i__1 = jb;
    for (k = 1; k <= i__1; ++k) {
	if (kout[k - 1] != 0) {
	    i__ = k;
L120:
	    m = kout[i__ - 1];
	    kout[i__ - 1] = 0;
	    ++j;
	    i__2 = jb;
	    for (ii = 1; ii <= i__2; ++ii) {
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
	io___256.ciunit = iounit_1.iout;
	s_wsfe(&io___256);
	do_fio(&c__1, (char *)&(*ir), (ftnlen)sizeof(integer));
	e_wsfe();
    } else {
	moved = TRUE_;
	xr += rmove;
	yr += rmove;
	zr += rmove;
	goto L10;
    }

/*     compute the exposed surface area for the sphere of interest */

L150:
    *area = ib * pix2 + exang + arclen;
    *area = d_mod(area, &c_b7) * rrsq;

/*     attempt to fix negative area by moving atom slightly */

    if (*area < 0.) {
	if (moved) {
	    io___257.ciunit = iounit_1.iout;
	    s_wsfe(&io___257);
	    do_fio(&c__1, (char *)&(*ir), (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    moved = TRUE_;
	    xr += rmove;
	    yr += rmove;
	    zr += rmove;
	    goto L10;
	}
    }
L170:
    return 0;
} /* surfatom1_ */

#undef darea_ref


