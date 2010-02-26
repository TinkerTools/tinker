/* gradrgd.f -- translated by f2c (version 20050501).
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
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    doublereal xrb[25000], yrb[25000], zrb[25000], rbc[6000]	/* was [6][
	    1000] */;
    logical use_rigid__;
} rigid_;

#define rigid_1 rigid_



/*     ############################################################## */
/*     ##  COPYRIGHT (C) 1997 by Rohit Pappu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine gradrgd  --  energy & gradient of rigid body  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "gradrgd" calls subroutines to calculate the potential energy */
/*     and first derivatives with respect to rigid body coordinates */


/* Subroutine */ int gradrgd_(doublereal *energy, doublereal *derivs)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static doublereal g[75000]	/* was [3][25000] */;
    static integer i__, j, k;
    static doublereal phi, xcm, ycm, zcm, tau[3], cphi, ephi[3], epsi[3];
    static integer init;
    static doublereal sphi;
    static integer stop;
    static doublereal theta, xterm, yterm, zterm, ctheta, etheta[3], stheta;


#define g_ref(a_1,a_2) g[(a_2)*3 + a_1 - 4]
#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]
#define derivs_ref(a_1,a_2) derivs[(a_2)*6 + a_1]



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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  rigid.i  --  rigid body coordinates for atom groups  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     xrb         rigid body reference x-coordinate for each atom */
/*     yrb         rigid body reference y-coordinate for each atom */
/*     zrb         rigid body reference z-coordinate for each atom */
/*     rbc         current rigid body coordinates for each group */
/*     use_rigid   flag to mark use of rigid body coordinate system */




/*     zero out the total of rigid body derivative components */

    /* Parameter adjustments */
    derivs -= 7;

    /* Function Body */
    for (i__ = 1; i__ <= 1000; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    derivs_ref(j, i__) = 0.;
	}
    }

/*     calculate the energy and Cartesian first derivatives */

    gradient_(energy, g);

/*     compute the rigid body gradient components for each group */

    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	init = igrp_ref(1, i__);
	stop = igrp_ref(2, i__);
	xcm = rbc_ref(1, i__);
	ycm = rbc_ref(2, i__);
	zcm = rbc_ref(3, i__);
	phi = rbc_ref(4, i__);
	theta = rbc_ref(5, i__);
	cphi = cos(phi);
	sphi = sin(phi);
	ctheta = cos(theta);
	stheta = sin(theta);

/*     get unit vectors along the phi, theta and psi rotation axes */

	ephi[0] = 0.;
	ephi[1] = 0.;
	ephi[2] = 1.;
	etheta[0] = -sphi;
	etheta[1] = cphi;
	etheta[2] = 0.;
	epsi[0] = ctheta * cphi;
	epsi[1] = ctheta * sphi;
	epsi[2] = -stheta;

/*     first, get the rigid body gradients for translations */

	i__2 = stop;
	for (j = init; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    derivs_ref(1, i__) = derivs_ref(1, i__) + g_ref(1, k);
	    derivs_ref(2, i__) = derivs_ref(2, i__) + g_ref(2, k);
	    derivs_ref(3, i__) = derivs_ref(3, i__) + g_ref(3, k);
	}

/*     accumulate the moment arm along each axis of rotation */

	for (j = 1; j <= 3; ++j) {
	    tau[j - 1] = 0.;
	}
	i__2 = stop;
	for (j = init; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    xterm = atoms_1.x[k - 1] - xcm;
	    yterm = atoms_1.y[k - 1] - ycm;
	    zterm = atoms_1.z__[k - 1] - zcm;
	    tau[0] = tau[0] + yterm * g_ref(3, k) - zterm * g_ref(2, k);
	    tau[1] = tau[1] + zterm * g_ref(1, k) - xterm * g_ref(3, k);
	    tau[2] = tau[2] + xterm * g_ref(2, k) - yterm * g_ref(1, k);
	}

/*     now, set the rigid body gradients for rotations */

	for (j = 1; j <= 3; ++j) {
	    derivs_ref(4, i__) = derivs_ref(4, i__) + tau[j - 1] * ephi[j - 1]
		    ;
	    derivs_ref(5, i__) = derivs_ref(5, i__) + tau[j - 1] * etheta[j - 
		    1];
	    derivs_ref(6, i__) = derivs_ref(6, i__) + tau[j - 1] * epsi[j - 1]
		    ;
	}
    }
    return 0;
} /* gradrgd_ */

#undef derivs_ref
#undef igrp_ref
#undef rbc_ref
#undef g_ref


