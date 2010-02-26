/* xyzatm.f -- translated by f2c (version 20050501).
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

/* Table of constant values */

static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine xyzatm  --  single atom internal to Cartesian  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "xyzatm" computes the Cartesian coordinates of a single */
/*     atom from its defining internal coordinate values */


/* Subroutine */ int xyzatm_(integer *i__, integer *ia, doublereal *bond, 
	integer *ib, doublereal *angle1, integer *ic, doublereal *angle2, 
	integer *chiral)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 XYZATM  --  Undefined Dihedral\002,\002 "
	    "Angle at Atom\002,i6)";
    static char fmt_20[] = "(/,\002 XYZATM  --  Defining Atoms Colinear\002"
	    ",\002 at Atom\002,i6)";
    static char fmt_30[] = "(/,\002 XYZATM  --  Sum of Bond Angles\002,\002 "
	    "Too Large at Atom\002,i6)";

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal), sqrt(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal a, b, c__, xt, yt, zt, xu, yu, zu, rab, rba, rbc, rac, 
	    xab, yab, zab, xba, yba, zba, xbc, ybc, zbc, xac, yac, zac, eps, 
	    rad1, rad2, cos1, cos2, sin1, sin2, cosb, sinb, cosg, sine, sing, 
	    xtmp, ztmp, sine2, cosine;

    /* Fortran I/O blocks */
    static cilist io___27 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_30, 0 };




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




/*     convert angles to radians, and get their sines and cosines */

    eps = 1e-8;
    rad1 = *angle1 / 57.29577951308232088;
    rad2 = *angle2 / 57.29577951308232088;
    sin1 = sin(rad1);
    cos1 = cos(rad1);
    sin2 = sin(rad2);
    cos2 = cos(rad2);

/*     if no second site given, place the atom at the origin */

    if (*ia == 0) {
	atoms_1.x[*i__ - 1] = 0.;
	atoms_1.y[*i__ - 1] = 0.;
	atoms_1.z__[*i__ - 1] = 0.;

/*     if no third site given, place the atom along the z-axis */

    } else if (*ib == 0) {
	atoms_1.x[*i__ - 1] = atoms_1.x[*ia - 1];
	atoms_1.y[*i__ - 1] = atoms_1.y[*ia - 1];
	atoms_1.z__[*i__ - 1] = atoms_1.z__[*ia - 1] + *bond;

/*     if no fourth site given, place the atom in the x,z-plane */

    } else if (*ic == 0) {
	xab = atoms_1.x[*ia - 1] - atoms_1.x[*ib - 1];
	yab = atoms_1.y[*ia - 1] - atoms_1.y[*ib - 1];
	zab = atoms_1.z__[*ia - 1] - atoms_1.z__[*ib - 1];
/* Computing 2nd power */
	d__1 = xab;
/* Computing 2nd power */
	d__2 = yab;
/* Computing 2nd power */
	d__3 = zab;
	rab = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	xab /= rab;
	yab /= rab;
	zab /= rab;
	cosb = zab;
/* Computing 2nd power */
	d__1 = xab;
/* Computing 2nd power */
	d__2 = yab;
	sinb = sqrt(d__1 * d__1 + d__2 * d__2);
	if (sinb == 0.) {
	    cosg = 1.;
	    sing = 0.;
	} else {
	    cosg = yab / sinb;
	    sing = xab / sinb;
	}
	xtmp = *bond * sin1;
	ztmp = rab - *bond * cos1;
	atoms_1.x[*i__ - 1] = atoms_1.x[*ib - 1] + xtmp * cosg + ztmp * sing *
		 sinb;
	atoms_1.y[*i__ - 1] = atoms_1.y[*ib - 1] - xtmp * sing + ztmp * cosg *
		 sinb;
	atoms_1.z__[*i__ - 1] = atoms_1.z__[*ib - 1] + ztmp * cosb;

/*     general case where the second angle is a dihedral angle */

    } else if (*chiral == 0) {
	xab = atoms_1.x[*ia - 1] - atoms_1.x[*ib - 1];
	yab = atoms_1.y[*ia - 1] - atoms_1.y[*ib - 1];
	zab = atoms_1.z__[*ia - 1] - atoms_1.z__[*ib - 1];
/* Computing 2nd power */
	d__1 = xab;
/* Computing 2nd power */
	d__2 = yab;
/* Computing 2nd power */
	d__3 = zab;
	rab = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	xab /= rab;
	yab /= rab;
	zab /= rab;
	xbc = atoms_1.x[*ib - 1] - atoms_1.x[*ic - 1];
	ybc = atoms_1.y[*ib - 1] - atoms_1.y[*ic - 1];
	zbc = atoms_1.z__[*ib - 1] - atoms_1.z__[*ic - 1];
/* Computing 2nd power */
	d__1 = xbc;
/* Computing 2nd power */
	d__2 = ybc;
/* Computing 2nd power */
	d__3 = zbc;
	rbc = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	xbc /= rbc;
	ybc /= rbc;
	zbc /= rbc;
	xt = zab * ybc - yab * zbc;
	yt = xab * zbc - zab * xbc;
	zt = yab * xbc - xab * ybc;
	cosine = xab * xbc + yab * ybc + zab * zbc;
/* Computing MAX */
/* Computing 2nd power */
	d__2 = cosine;
	d__1 = 1. - d__2 * d__2;
	sine = sqrt((max(d__1,eps)));
	if (abs(cosine) >= 1.) {
	    io___27.ciunit = iounit_1.iout;
	    s_wsfe(&io___27);
	    do_fio(&c__1, (char *)&(*i__), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	xt /= sine;
	yt /= sine;
	zt /= sine;
	xu = yt * zab - zt * yab;
	yu = zt * xab - xt * zab;
	zu = xt * yab - yt * xab;
	atoms_1.x[*i__ - 1] = atoms_1.x[*ia - 1] + *bond * (xu * sin1 * cos2 
		+ xt * sin1 * sin2 - xab * cos1);
	atoms_1.y[*i__ - 1] = atoms_1.y[*ia - 1] + *bond * (yu * sin1 * cos2 
		+ yt * sin1 * sin2 - yab * cos1);
	atoms_1.z__[*i__ - 1] = atoms_1.z__[*ia - 1] + *bond * (zu * sin1 * 
		cos2 + zt * sin1 * sin2 - zab * cos1);

/*     general case where the second angle is a bond angle */

    } else if (abs(*chiral) == 1) {
	xba = atoms_1.x[*ib - 1] - atoms_1.x[*ia - 1];
	yba = atoms_1.y[*ib - 1] - atoms_1.y[*ia - 1];
	zba = atoms_1.z__[*ib - 1] - atoms_1.z__[*ia - 1];
/* Computing 2nd power */
	d__1 = xba;
/* Computing 2nd power */
	d__2 = yba;
/* Computing 2nd power */
	d__3 = zba;
	rba = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	xba /= rba;
	yba /= rba;
	zba /= rba;
	xac = atoms_1.x[*ia - 1] - atoms_1.x[*ic - 1];
	yac = atoms_1.y[*ia - 1] - atoms_1.y[*ic - 1];
	zac = atoms_1.z__[*ia - 1] - atoms_1.z__[*ic - 1];
/* Computing 2nd power */
	d__1 = xac;
/* Computing 2nd power */
	d__2 = yac;
/* Computing 2nd power */
	d__3 = zac;
	rac = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	xac /= rac;
	yac /= rac;
	zac /= rac;
	xt = zba * yac - yba * zac;
	yt = xba * zac - zba * xac;
	zt = yba * xac - xba * yac;
	cosine = xba * xac + yba * yac + zba * zac;
/* Computing MAX */
/* Computing 2nd power */
	d__2 = cosine;
	d__1 = 1. - d__2 * d__2;
	sine2 = max(d__1,eps);
	if (abs(cosine) >= 1.) {
	    io___40.ciunit = iounit_1.iout;
	    s_wsfe(&io___40);
	    do_fio(&c__1, (char *)&(*i__), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	a = (-cos2 - cosine * cos1) / sine2;
	b = (cos1 + cosine * cos2) / sine2;
	c__ = (a * cos2 + 1. - b * cos1) / sine2;
	if (c__ > eps) {
	    c__ = *chiral * sqrt(c__);
	} else if (c__ < -eps) {
/* Computing 2nd power */
	    d__1 = a * xac + b * xba;
/* Computing 2nd power */
	    d__2 = a * yac + b * yba;
/* Computing 2nd power */
	    d__3 = a * zac + b * zba;
	    c__ = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	    a /= c__;
	    b /= c__;
	    c__ = 0.;
	    if (inform_1.debug) {
		io___44.ciunit = iounit_1.iout;
		s_wsfe(&io___44);
		do_fio(&c__1, (char *)&(*ia), (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	} else {
	    c__ = 0.;
	}
	atoms_1.x[*i__ - 1] = atoms_1.x[*ia - 1] + *bond * (a * xac + b * xba 
		+ c__ * xt);
	atoms_1.y[*i__ - 1] = atoms_1.y[*ia - 1] + *bond * (a * yac + b * yba 
		+ c__ * yt);
	atoms_1.z__[*i__ - 1] = atoms_1.z__[*ia - 1] + *bond * (a * zac + b * 
		zba + c__ * zt);
    }
    return 0;
} /* xyzatm_ */

