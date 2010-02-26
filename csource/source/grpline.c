/* grpline.f -- translated by f2c (version 20050501).
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
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

struct {
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

struct {
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    doublereal vcm[3000]	/* was [3][1000] */, wcm[3000]	/* was [3][
	    1000] */, lm[3000]	/* was [3][1000] */, vc[3000]	/* was [3][
	    1000] */, wc[3000]	/* was [3][1000] */;
    logical linear[1000];
} rgddyn_;

#define rgddyn_1 rgddyn_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine grpline  --  test atom groups for linearity  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "grpline" tests each atom group for linearity of the sites */
/*     contained in the group */


/* Subroutine */ int grpline_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal x2, y2, z2, xx, yy, zz, det, rcm[3], eps, xcm[25000], 
	    ycm[25000], zcm[25000];
    static integer size, stop;
    static doublereal weigh, inert[6];
    static integer start;


#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]



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
/*     ##  atmtyp.i  --  atomic properties for each current atom  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     mass      atomic weight for each atom in the system */
/*     tag       integer atom labels from input coordinates file */
/*     class     atom class number for each atom in the system */
/*     atomic    atomic number for each atom in the system */
/*     valence   valence number for each atom in the system */
/*     name      atom name for each atom in the system */
/*     story     descriptive type for each atom in system */




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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  group.i  --  partitioning of system into atom groups  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     grpmass     total mass of all the atoms in each group */
/*     wgrp        weight for each set of group-group interactions */
/*     ngrp        total number of atom groups in the system */
/*     kgrp        contiguous list of the atoms in each group */
/*     igrp        first and last atom of each group in the list */
/*     grplist     number of the group to which each atom belongs */
/*     use_group   flag to use partitioning of system into groups */
/*     use_intra   flag to include only intragroup interactions */
/*     use_inter   flag to include only intergroup interactions */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  rgddyn.i  --  velocities and momenta for rigid body MD  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     vcm     current translational velocity of each rigid body */
/*     wcm     current angular velocity of each rigid body */
/*     lm      current angular momentum of each rigid body */
/*     vc      half-step translational velocity for kinetic energy */
/*     wc      half-step angular velocity for kinetic energy */
/*     linear  logical flag to mark group as linear or nonlinear */




/*     compute the center of mass coordinates of each group */

    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	start = igrp_ref(1, i__);
	stop = igrp_ref(2, i__);
	for (j = 1; j <= 3; ++j) {
	    rcm[j - 1] = 0.;
	}
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    weigh = atmtyp_1.mass[k - 1];
	    rcm[0] += atoms_1.x[k - 1] * weigh;
	    rcm[1] += atoms_1.y[k - 1] * weigh;
	    rcm[2] += atoms_1.z__[k - 1] * weigh;
	}
/* Computing MAX */
	d__1 = 1., d__2 = group_1.grpmass[i__ - 1];
	weigh = max(d__1,d__2);
	for (j = 1; j <= 3; ++j) {
	    rcm[j - 1] /= weigh;
	}
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    xcm[k - 1] = atoms_1.x[k - 1] - rcm[0];
	    ycm[k - 1] = atoms_1.y[k - 1] - rcm[1];
	    zcm[k - 1] = atoms_1.z__[k - 1] - rcm[2];
	}
    }

/*     use the moments of inertia to check for linearity */

    eps = 1e-8;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	size = igrp_ref(2, i__) - igrp_ref(1, i__) + 1;
	rgddyn_1.linear[i__ - 1] = FALSE_;
	if (size == 2) {
	    rgddyn_1.linear[i__ - 1] = TRUE_;
	} else if (size > 2) {
	    for (j = 1; j <= 6; ++j) {
		inert[j - 1] = 0.;
	    }
	    i__2 = igrp_ref(2, i__);
	    for (j = igrp_ref(1, i__); j <= i__2; ++j) {
		k = group_1.kgrp[j - 1];
		xx = xcm[k - 1];
		yy = ycm[k - 1];
		zz = zcm[k - 1];
		x2 = xx * xx;
		y2 = yy * yy;
		z2 = zz * zz;
		weigh = atmtyp_1.mass[k - 1];
		inert[0] += weigh * (y2 + z2);
		inert[1] -= weigh * xx * yy;
		inert[2] += weigh * (x2 + z2);
		inert[3] -= weigh * xx * zz;
		inert[4] -= weigh * yy * zz;
		inert[5] += weigh * (x2 + y2);
	    }
	    det = inert[0] * inert[2] * inert[5] + inert[1] * 2. * inert[4] * 
		    inert[3] - inert[2] * inert[3] * inert[3] - inert[0] * 
		    inert[4] * inert[4] - inert[1] * inert[1] * inert[5];
	    if (abs(det) < eps) {
		rgddyn_1.linear[i__ - 1] = TRUE_;
	    }
	}
    }
    return 0;
} /* grpline_ */

#undef igrp_ref


