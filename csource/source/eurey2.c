/* eurey2.f -- translated by f2c (version 20050501).
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
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    doublereal hessx[75000]	/* was [3][25000] */, hessy[75000]	/* 
	    was [3][25000] */, hessz[75000]	/* was [3][25000] */;
} hessn_;

#define hessn_1 hessn_

struct {
    doublereal uk[75000], ul[75000];
    integer nurey, iury[225000]	/* was [3][75000] */;
} urey_;

#define urey_1 urey_

struct {
    doublereal cury, qury, ureyunit;
} urypot_;

#define urypot_1 urypot_

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine eurey2  --  atom-by-atom Urey-Bradley Hessian  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "eurey2" calculates second derivatives of the Urey-Bradley */
/*     interaction energy for a single atom at a time */


/* Subroutine */ int eurey2_(integer *i__)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j;
    static doublereal de;
    static integer ia, ic;
    static doublereal dt, d2e[9]	/* was [3][3] */, dt2, rac, xac, yac, 
	    zac, rac2, fgrp, term, ideal;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal deddt, force;
    static integer iurey;
    static doublereal termx, termy, termz, d2eddt2;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;


#define d2e_ref(a_1,a_2) d2e[(a_2)*3 + a_1 - 4]
#define iury_ref(a_1,a_2) urey_1.iury[(a_2)*3 + a_1 - 4]
#define hessx_ref(a_1,a_2) hessn_1.hessx[(a_2)*3 + a_1 - 4]
#define hessy_ref(a_1,a_2) hessn_1.hessy[(a_2)*3 + a_1 - 4]
#define hessz_ref(a_1,a_2) hessn_1.hessz[(a_2)*3 + a_1 - 4]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  bound.i  --  control of periodic boundary conditions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     polycut       cutoff distance for infinite polymer nonbonds */
/*     polycut2      square of infinite polymer nonbond cutoff */
/*     use_bounds    flag to use periodic boundary conditions */
/*     use_replica   flag to use replicates for periodic system */
/*     use_polymer   flag to mark presence of infinite polymer */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  couple.i  --  near-neighbor atom connectivity lists   ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxn13   maximum number of atoms 1-3 connected to an atom */
/*     maxn14   maximum number of atoms 1-4 connected to an atom */
/*     maxn15   maximum number of atoms 1-5 connected to an atom */

/*     n12      number of atoms directly bonded to each atom */
/*     i12      atom numbers of atoms 1-2 connected to each atom */
/*     n13      number of atoms in a 1-3 relation to each atom */
/*     i13      atom numbers of atoms 1-3 connected to each atom */
/*     n14      number of atoms in a 1-4 relation to each atom */
/*     i14      atom numbers of atoms 1-4 connected to each atom */
/*     n15      number of atoms in a 1-5 relation to each atom */
/*     i15      atom numbers of atoms 1-5 connected to each atom */




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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  hessn.i  --  Cartesian Hessian elements for a single atom  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     hessx   Hessian elements for x-component of current atom */
/*     hessy   Hessian elements for y-component of current atom */
/*     hessz   Hessian elements for z-component of current atom */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  urey.i  --  Urey-Bradley interactions in the structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     uk      Urey-Bradley force constants (kcal/mole/Ang**2) */
/*     ul      ideal 1-3 distance values in Angstroms */
/*     nurey   total number of Urey-Bradley terms in the system */
/*     iury    numbers of the atoms in each Urey-Bradley interaction */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  urypot.i  --  specifics of Urey-Bradley functional form  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     cury       cubic coefficient in Urey-Bradley potential */
/*     qury       quartic coefficient in Urey-Bradley potential */
/*     ureyunit   convert Urey-Bradley energy to kcal/mole */




/*     compute the Hessian elements of the Urey-Bradley energy */

    i__1 = urey_1.nurey;
    for (iurey = 1; iurey <= i__1; ++iurey) {
	ia = iury_ref(1, iurey);
	ic = iury_ref(2, iurey);
	ideal = urey_1.ul[iurey - 1];
	force = urey_1.uk[iurey - 1];

/*     decide whether to compute the current interaction */

	proceed = *i__ == ia || *i__ == ic;
	if (proceed && group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ic, &c__0, &c__0, &c__0, &c__0);
	}

/*     compute the value of the 1-3 distance deviation */

	if (proceed) {
	    if (*i__ == ic) {
		ic = ia;
		ia = *i__;
	    }
	    xac = atoms_1.x[ia - 1] - atoms_1.x[ic - 1];
	    yac = atoms_1.y[ia - 1] - atoms_1.y[ic - 1];
	    zac = atoms_1.z__[ia - 1] - atoms_1.z__[ic - 1];
	    if (bound_1.use_polymer__) {
		image_(&xac, &yac, &zac);
	    }
	    rac2 = xac * xac + yac * yac + zac * zac;
	    rac = sqrt(rac2);
	    dt = rac - ideal;
	    dt2 = dt * dt;
	    deddt = urypot_1.ureyunit * 2. * force * dt * (urypot_1.cury * 
		    1.5 * dt + 1. + urypot_1.qury * 2. * dt2);
	    d2eddt2 = urypot_1.ureyunit * 2. * force * (urypot_1.cury * 3. * 
		    dt + 1. + urypot_1.qury * 6. * dt2);

/*     scale the interaction based on its group membership */

	    if (group_1.use_group__) {
		deddt *= fgrp;
		d2eddt2 *= fgrp;
	    }

/*     set the chain rule terms for the Hessian elements */

	    de = deddt / rac;
	    term = (d2eddt2 - de) / rac2;
	    termx = term * xac;
	    termy = term * yac;
	    termz = term * zac;
	    d2e_ref(1, 1) = termx * xac + de;
	    d2e_ref(1, 2) = termx * yac;
	    d2e_ref(1, 3) = termx * zac;
	    d2e_ref(2, 1) = d2e_ref(1, 2);
	    d2e_ref(2, 2) = termy * yac + de;
	    d2e_ref(2, 3) = termy * zac;
	    d2e_ref(3, 1) = d2e_ref(1, 3);
	    d2e_ref(3, 2) = d2e_ref(2, 3);
	    d2e_ref(3, 3) = termz * zac + de;

/*     increment diagonal and non-diagonal Hessian elements */

	    for (j = 1; j <= 3; ++j) {
		hessx_ref(j, ia) = hessx_ref(j, ia) + d2e_ref(1, j);
		hessy_ref(j, ia) = hessy_ref(j, ia) + d2e_ref(2, j);
		hessz_ref(j, ia) = hessz_ref(j, ia) + d2e_ref(3, j);
		hessx_ref(j, ic) = hessx_ref(j, ic) - d2e_ref(1, j);
		hessy_ref(j, ic) = hessy_ref(j, ic) - d2e_ref(2, j);
		hessz_ref(j, ic) = hessz_ref(j, ic) - d2e_ref(3, j);
	    }
	}
    }
    return 0;
} /* eurey2_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef iury_ref
#undef d2e_ref


