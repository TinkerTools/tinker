/* geometry.f -- translated by f2c (version 20050501).
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



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function geometry  --  evaluate distance, angle, torsion  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "geometry" finds the value of the interatomic distance, angle */
/*     or dihedral angle defined by two to four input atoms */


doublereal geometry_(integer *ia, integer *ib, integer *ic, integer *id)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal);

    /* Local variables */
    static doublereal xt, yt, zt, xu, yu, zu, rt2, ru2, xab, yab, zab, xba, 
	    yba, zba, xcb, ycb, zcb, xdc, ydc, zdc, rab2, rcb2, rabc, sign, 
	    rtru, cosine;



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




/*     set default in case atoms are coincident or colinear */

    ret_val = 0.;

/*     compute the value of the distance in angstroms */

    if (*ic == 0) {
	xab = atoms_1.x[*ia - 1] - atoms_1.x[*ib - 1];
	yab = atoms_1.y[*ia - 1] - atoms_1.y[*ib - 1];
	zab = atoms_1.z__[*ia - 1] - atoms_1.z__[*ib - 1];
	ret_val = sqrt(xab * xab + yab * yab + zab * zab);

/*     compute the value of the angle in degrees */

    } else if (*id == 0) {
	xab = atoms_1.x[*ia - 1] - atoms_1.x[*ib - 1];
	yab = atoms_1.y[*ia - 1] - atoms_1.y[*ib - 1];
	zab = atoms_1.z__[*ia - 1] - atoms_1.z__[*ib - 1];
	xcb = atoms_1.x[*ic - 1] - atoms_1.x[*ib - 1];
	ycb = atoms_1.y[*ic - 1] - atoms_1.y[*ib - 1];
	zcb = atoms_1.z__[*ic - 1] - atoms_1.z__[*ib - 1];
	rab2 = xab * xab + yab * yab + zab * zab;
	rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
	rabc = sqrt(rab2 * rcb2);
	if (rabc != 0.) {
	    cosine = (xab * xcb + yab * ycb + zab * zcb) / rabc;
/* Computing MIN */
	    d__1 = 1., d__2 = max(-1.,cosine);
	    cosine = min(d__1,d__2);
	    ret_val = acos(cosine) * 57.29577951308232088;
	}

/*     compute the value of the dihedral angle in degrees */

    } else {
	xba = atoms_1.x[*ib - 1] - atoms_1.x[*ia - 1];
	yba = atoms_1.y[*ib - 1] - atoms_1.y[*ia - 1];
	zba = atoms_1.z__[*ib - 1] - atoms_1.z__[*ia - 1];
	xcb = atoms_1.x[*ic - 1] - atoms_1.x[*ib - 1];
	ycb = atoms_1.y[*ic - 1] - atoms_1.y[*ib - 1];
	zcb = atoms_1.z__[*ic - 1] - atoms_1.z__[*ib - 1];
	xdc = atoms_1.x[*id - 1] - atoms_1.x[*ic - 1];
	ydc = atoms_1.y[*id - 1] - atoms_1.y[*ic - 1];
	zdc = atoms_1.z__[*id - 1] - atoms_1.z__[*ic - 1];
	xt = yba * zcb - ycb * zba;
	yt = xcb * zba - xba * zcb;
	zt = xba * ycb - xcb * yba;
	xu = ycb * zdc - ydc * zcb;
	yu = xdc * zcb - xcb * zdc;
	zu = xcb * ydc - xdc * ycb;
	rt2 = xt * xt + yt * yt + zt * zt;
	ru2 = xu * xu + yu * yu + zu * zu;
	rtru = sqrt(rt2 * ru2);
	if (rtru != 0.) {
	    cosine = (xt * xu + yt * yu + zt * zu) / rtru;
/* Computing MIN */
	    d__1 = 1., d__2 = max(-1.,cosine);
	    cosine = min(d__1,d__2);
	    ret_val = acos(cosine) * 57.29577951308232088;
	    sign = xba * xu + yba * yu + zba * zu;
	    if (sign < 0.) {
		ret_val = -ret_val;
	    }
	}
    }
    return ret_val;
} /* geometry_ */

