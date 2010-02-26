/* hessrgd.f -- translated by f2c (version 20050501).
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



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine hessrgd  --  rigid body Hessian elements  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "hessrgd" computes the numerical Hessian elements with */
/*     respect to rigid body coordinates via 6*ngroup+1 gradient */
/*     evaluations */


/* Subroutine */ int hessrgd_(doublereal *hrigid)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int rigidxyz_(void);
    static doublereal e, g[6000]	/* was [6][1000] */;
    static integer i__, j, k, m;
    static doublereal g0[6000]	/* was [6][1000] */, old, eps;
    static integer nvar;
    extern /* Subroutine */ int gradrgd_(doublereal *, doublereal *);


#define g_ref(a_1,a_2) g[(a_2)*6 + a_1 - 7]
#define g0_ref(a_1,a_2) g0[(a_2)*6 + a_1 - 7]
#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]
#define hrigid_ref(a_1,a_2) hrigid[(a_2)*6000 + a_1]



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




/*     calculate base values for the rigid body gradient */

    /* Parameter adjustments */
    hrigid -= 6001;

    /* Function Body */
    eps = 1e-5;
    gradrgd_(&e, g0);

/*     compute one-sided numerical Hessian from gradient values; */
/*     set off-diagonal elements to the average symmetric value */

    nvar = group_1.ngrp * 6;
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = (i__ - 1) / 6 + 1;
	k = (i__ - 1) % 6 + 1;
	old = rbc_ref(k, j);
	rbc_ref(k, j) = rbc_ref(k, j) + eps;
	rigidxyz_();
	gradrgd_(&e, g);
	rbc_ref(k, j) = old;
	i__2 = nvar;
	for (m = 1; m <= i__2; ++m) {
	    j = (m - 1) / 6 + 1;
	    k = (m - 1) % 6 + 1;
	    hrigid_ref(m, i__) = (g_ref(k, j) - g0_ref(k, j)) / eps;
	}
	i__2 = i__ - 1;
	for (m = 1; m <= i__2; ++m) {
	    hrigid_ref(m, i__) = (hrigid_ref(m, i__) + hrigid_ref(i__, m)) * 
		    .5;
	    hrigid_ref(i__, m) = hrigid_ref(m, i__);
	}
    }

/*     restore the Cartesian coordinates to original values */

    rigidxyz_();
    return 0;
} /* hessrgd_ */

#undef hrigid_ref
#undef rbc_ref
#undef g0_ref
#undef g_ref


