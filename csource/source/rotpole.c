/* rotpole.f -- translated by f2c (version 20050501).
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
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

struct {
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_



/*     ############################################################ */
/*     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ## */
/*     ##                  All Rights Reserved                   ## */
/*     ############################################################ */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine rotpole  --  rotate multipoles to global frame  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "rotpole" constructs the set of atomic multipoles in the global */
/*     frame by applying the correct rotation matrix for each site */


/* Subroutine */ int rotpole_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal a[9]	/* was [3][3] */;
    static integer i__;
    extern /* Subroutine */ int rotmat_(integer *, doublereal *), rotsite_(
	    integer *, doublereal *);



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
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




/*     rotate the atomic multipoles at each site in turn */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rotmat_(&i__, a);
	rotsite_(&i__, a);
    }
    return 0;
} /* rotpole_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine rotmat  --  find global frame rotation matrix  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "rotmat" finds the rotation matrix that converts from the local */
/*     coordinate system to the global frame at a multipole site */


/* Subroutine */ int rotmat_(integer *i__, doublereal *a)
{
    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static doublereal r__;
    static integer ii;
    static doublereal dx, dy, dz;
    static integer ix, iy, iz;
    static doublereal xi, yi, zi, dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3,
	     dot;


#define a_ref(a_1,a_2) a[(a_2)*3 + a_1]
#define polaxe_ref(a_0,a_1) &mpole_1.polaxe[(a_1)*8 + a_0 - 8]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




/*     get coordinates and frame definition for the multipole site */

    /* Parameter adjustments */
    a -= 4;

    /* Function Body */
    ii = mpole_1.ipole[*i__ - 1];
    xi = atoms_1.x[ii - 1];
    yi = atoms_1.y[ii - 1];
    zi = atoms_1.z__[ii - 1];
    ix = mpole_1.xaxis[*i__ - 1];
    iy = mpole_1.yaxis[*i__ - 1];
    iz = mpole_1.zaxis[*i__ - 1];

/*     use the identity matrix as the default rotation matrix */

    a_ref(1, 1) = 1.;
    a_ref(2, 1) = 0.;
    a_ref(3, 1) = 0.;
    a_ref(1, 3) = 0.;
    a_ref(2, 3) = 0.;
    a_ref(3, 3) = 1.;

/*     z-then-x method rotation matrix elements for z- and x-axes */

    if (s_cmp(polaxe_ref(0, *i__), "Z-then-X", (ftnlen)8, (ftnlen)8) == 0) {
	dx = atoms_1.x[iz - 1] - xi;
	dy = atoms_1.y[iz - 1] - yi;
	dz = atoms_1.z__[iz - 1] - zi;
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	a_ref(1, 3) = dx / r__;
	a_ref(2, 3) = dy / r__;
	a_ref(3, 3) = dz / r__;
	dx = atoms_1.x[ix - 1] - xi;
	dy = atoms_1.y[ix - 1] - yi;
	dz = atoms_1.z__[ix - 1] - zi;
	dot = dx * a_ref(1, 3) + dy * a_ref(2, 3) + dz * a_ref(3, 3);
	dx -= dot * a_ref(1, 3);
	dy -= dot * a_ref(2, 3);
	dz -= dot * a_ref(3, 3);
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	a_ref(1, 1) = dx / r__;
	a_ref(2, 1) = dy / r__;
	a_ref(3, 1) = dz / r__;

/*     bisector method rotation matrix elements for z- and x-axes */

    } else if (s_cmp(polaxe_ref(0, *i__), "Bisector", (ftnlen)8, (ftnlen)8) ==
	     0) {
	dx = atoms_1.x[iz - 1] - xi;
	dy = atoms_1.y[iz - 1] - yi;
	dz = atoms_1.z__[iz - 1] - zi;
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	dx1 = dx / r__;
	dy1 = dy / r__;
	dz1 = dz / r__;
	dx = atoms_1.x[ix - 1] - xi;
	dy = atoms_1.y[ix - 1] - yi;
	dz = atoms_1.z__[ix - 1] - zi;
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	dx2 = dx / r__;
	dy2 = dy / r__;
	dz2 = dz / r__;
	dx = dx1 + dx2;
	dy = dy1 + dy2;
	dz = dz1 + dz2;
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	a_ref(1, 3) = dx / r__;
	a_ref(2, 3) = dy / r__;
	a_ref(3, 3) = dz / r__;
	dot = dx2 * a_ref(1, 3) + dy2 * a_ref(2, 3) + dz2 * a_ref(3, 3);
	dx = dx2 - dot * a_ref(1, 3);
	dy = dy2 - dot * a_ref(2, 3);
	dz = dz2 - dot * a_ref(3, 3);
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	a_ref(1, 1) = dx / r__;
	a_ref(2, 1) = dy / r__;
	a_ref(3, 1) = dz / r__;

/*     z-bisect method rotation matrix elements for z- and x-axes */

    } else if (s_cmp(polaxe_ref(0, *i__), "Z-Bisect", (ftnlen)8, (ftnlen)8) ==
	     0) {
	dx = atoms_1.x[iz - 1] - xi;
	dy = atoms_1.y[iz - 1] - yi;
	dz = atoms_1.z__[iz - 1] - zi;
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	a_ref(1, 3) = dx / r__;
	a_ref(2, 3) = dy / r__;
	a_ref(3, 3) = dz / r__;
	dx = atoms_1.x[ix - 1] - xi;
	dy = atoms_1.y[ix - 1] - yi;
	dz = atoms_1.z__[ix - 1] - zi;
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	dx1 = dx / r__;
	dy1 = dy / r__;
	dz1 = dz / r__;
	dx = atoms_1.x[iy - 1] - xi;
	dy = atoms_1.y[iy - 1] - yi;
	dz = atoms_1.z__[iy - 1] - zi;
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	dx2 = dx / r__;
	dy2 = dy / r__;
	dz2 = dz / r__;
	dx = dx1 + dx2;
	dy = dy1 + dy2;
	dz = dz1 + dz2;
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	dx /= r__;
	dy /= r__;
	dz /= r__;
	dot = dx * a_ref(1, 3) + dy * a_ref(2, 3) + dz * a_ref(3, 3);
	dx -= dot * a_ref(1, 3);
	dy -= dot * a_ref(2, 3);
	dz -= dot * a_ref(3, 3);
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	a_ref(1, 1) = dx / r__;
	a_ref(2, 1) = dy / r__;
	a_ref(3, 1) = dz / r__;

/*     3-fold method rotation matrix elements for z- and x-axes */

    } else if (s_cmp(polaxe_ref(0, *i__), "3-Fold", (ftnlen)8, (ftnlen)6) == 
	    0) {
	dx = atoms_1.x[iz - 1] - xi;
	dy = atoms_1.y[iz - 1] - yi;
	dz = atoms_1.z__[iz - 1] - zi;
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	dx1 = dx / r__;
	dy1 = dy / r__;
	dz1 = dz / r__;
	dx = atoms_1.x[ix - 1] - xi;
	dy = atoms_1.y[ix - 1] - yi;
	dz = atoms_1.z__[ix - 1] - zi;
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	dx2 = dx / r__;
	dy2 = dy / r__;
	dz2 = dz / r__;
	dx = atoms_1.x[iy - 1] - xi;
	dy = atoms_1.y[iy - 1] - yi;
	dz = atoms_1.z__[iy - 1] - zi;
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	dx3 = dx / r__;
	dy3 = dy / r__;
	dz3 = dz / r__;
	dx = dx1 + dx2 + dx3;
	dy = dy1 + dy2 + dy3;
	dz = dz1 + dz2 + dz3;
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	a_ref(1, 3) = dx / r__;
	a_ref(2, 3) = dy / r__;
	a_ref(3, 3) = dz / r__;
	dot = dx2 * a_ref(1, 3) + dy2 * a_ref(2, 3) + dz2 * a_ref(3, 3);
	dx = dx2 - dot * a_ref(1, 3);
	dy = dy2 - dot * a_ref(2, 3);
	dz = dz2 - dot * a_ref(3, 3);
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	a_ref(1, 1) = dx / r__;
	a_ref(2, 1) = dy / r__;
	a_ref(3, 1) = dz / r__;
    }

/*     finally, find rotation matrix elements for the y-axis */

    a_ref(1, 2) = a_ref(3, 1) * a_ref(2, 3) - a_ref(2, 1) * a_ref(3, 3);
    a_ref(2, 2) = a_ref(1, 1) * a_ref(3, 3) - a_ref(3, 1) * a_ref(1, 3);
    a_ref(3, 2) = a_ref(2, 1) * a_ref(1, 3) - a_ref(1, 1) * a_ref(2, 3);
    return 0;
} /* rotmat_ */

#undef polaxe_ref
#undef a_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine rotsite  --  rotate multipoles at single site  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "rotsite" computes the atomic multipoles at a specified site */
/*     in the global coordinate frame by applying a rotation matrix */


/* Subroutine */ int rotsite_(integer *isite, doublereal *a)
{
    static integer i__, j, k, m;
    static doublereal m2[9]	/* was [3][3] */, r2[9]	/* was [3][3] */;


#define a_ref(a_1,a_2) a[(a_2)*3 + a_1]
#define m2_ref(a_1,a_2) m2[(a_2)*3 + a_1 - 4]
#define r2_ref(a_1,a_2) r2[(a_2)*3 + a_1 - 4]
#define pole_ref(a_1,a_2) mpole_1.pole[(a_2)*13 + a_1 - 14]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




/*     monopoles have the same value in any coordinate frame */

    /* Parameter adjustments */
    a -= 4;

    /* Function Body */
    rpole_ref(1, *isite) = pole_ref(1, *isite);

/*     rotate the dipoles to the global coordinate frame */

    for (i__ = 2; i__ <= 4; ++i__) {
	rpole_ref(i__, *isite) = 0.;
	for (j = 2; j <= 4; ++j) {
	    rpole_ref(i__, *isite) = rpole_ref(i__, *isite) + pole_ref(j, *
		    isite) * a_ref(i__ - 1, j - 1);
	}
    }

/*     rotate the quadrupoles to the global coordinate frame */

    k = 5;
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    m2_ref(i__, j) = pole_ref(k, *isite);
	    r2_ref(i__, j) = 0.;
	    ++k;
	}
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    if (j < i__) {
		r2_ref(i__, j) = r2_ref(j, i__);
	    } else {
		for (k = 1; k <= 3; ++k) {
		    for (m = 1; m <= 3; ++m) {
			r2_ref(i__, j) = r2_ref(i__, j) + a_ref(i__, k) * 
				a_ref(j, m) * m2_ref(k, m);
		    }
		}
	    }
	}
    }
    k = 5;
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    rpole_ref(k, *isite) = r2_ref(i__, j);
	    ++k;
	}
    }
    return 0;
} /* rotsite_ */

#undef rpole_ref
#undef pole_ref
#undef r2_ref
#undef m2_ref
#undef a_ref


