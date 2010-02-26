/* ebond2.f -- translated by f2c (version 20050501).
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
    integer bndlist[200000]	/* was [8][25000] */, anglist[700000]	/* 
	    was [28][25000] */;
} atmlst_;

#define atmlst_1 atmlst_

struct {
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

struct {
    doublereal cbnd, qbnd, bndunit;
    char bndtyp[8];
} bndpot_;

#define bndpot_1 bndpot_

struct {
    doublereal bk[50000], bl[50000];
    integer nbond, ibnd[100000]	/* was [2][50000] */;
} bond_;

#define bond_1 bond_

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

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine ebond2  --  atom-by-atom bond stretch Hessian  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "ebond2" calculates second derivatives of the bond */
/*     stretching energy for a single atom at a time */


/* Subroutine */ int ebond2_(integer *i__)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double exp(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal de;
    static integer ia, ib;
    static doublereal dt, d2e[9]	/* was [3][3] */, dt2, bde, rab, xab, 
	    yab, zab, rab2, fgrp, term, ideal;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal deddt, force, termx, termy, termz, d2eddt2;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;
    static doublereal expterm;


#define d2e_ref(a_1,a_2) d2e[(a_2)*3 + a_1 - 4]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define hessx_ref(a_1,a_2) hessn_1.hessx[(a_2)*3 + a_1 - 4]
#define hessy_ref(a_1,a_2) hessn_1.hessy[(a_2)*3 + a_1 - 4]
#define hessz_ref(a_1,a_2) hessn_1.hessz[(a_2)*3 + a_1 - 4]
#define bndlist_ref(a_1,a_2) atmlst_1.bndlist[(a_2)*8 + a_1 - 9]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  atmlst.i  --  local geometry terms involving each atom  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     bndlist   list of the bond numbers involving each atom */
/*     anglist   list of the angle numbers centered on each atom */




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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  bndpot.i  --  specifics of bond stretch functional forms  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     cbnd      cubic coefficient in bond stretch potential */
/*     qbnd      quartic coefficient in bond stretch potential */
/*     bndunit   convert bond stretch energy to kcal/mole */
/*     bndtyp    type of bond stretch potential energy function */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  bond.i  --  covalent bonds in the current structure  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     bk      bond stretch force constants (kcal/mole/Ang**2) */
/*     bl      ideal bond length values in Angstroms */
/*     nbond   total number of bond stretches in the system */
/*     ibnd    numbers of the atoms in each bond stretch */




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




/*     compute the Hessian elements of the bond stretch energy */

    ia = *i__;
    i__1 = couple_1.n12[ia - 1];
    for (k = 1; k <= i__1; ++k) {
	j = bndlist_ref(k, ia);
	if (ibnd_ref(1, j) == ia) {
	    ib = ibnd_ref(2, j);
	} else {
	    ib = ibnd_ref(1, j);
	}
	ideal = bond_1.bl[j - 1];
	force = bond_1.bk[j - 1];

/*     decide whether to compute the current interaction */

	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &c__0, &c__0, &c__0, &c__0);
	}

/*     compute the value of the bond length deviation */

	if (proceed) {
	    xab = atoms_1.x[ia - 1] - atoms_1.x[ib - 1];
	    yab = atoms_1.y[ia - 1] - atoms_1.y[ib - 1];
	    zab = atoms_1.z__[ia - 1] - atoms_1.z__[ib - 1];
	    if (bound_1.use_polymer__) {
		image_(&xab, &yab, &zab);
	    }
	    rab2 = xab * xab + yab * yab + zab * zab;
	    rab = sqrt(rab2);
	    dt = rab - ideal;

/*     harmonic potential uses Taylor expansion of Morse potential */
/*     through the fourth power of the bond length deviation */

	    if (s_cmp(bndpot_1.bndtyp, "HARMONIC", (ftnlen)8, (ftnlen)8) == 0)
		     {
		dt2 = dt * dt;
		deddt = bndpot_1.bndunit * 2. * force * dt * (bndpot_1.cbnd * 
			1.5 * dt + 1. + bndpot_1.qbnd * 2. * dt2);
		d2eddt2 = bndpot_1.bndunit * 2. * force * (bndpot_1.cbnd * 3. 
			* dt + 1. + bndpot_1.qbnd * 6. * dt2);

/*     Morse potential uses energy = BDE * (1 - e**(-alpha*dt))**2) */
/*     with the approximations alpha = sqrt(ForceConst/BDE) = -2 */
/*     and BDE = Bond Dissociation Energy = ForceConst/alpha**2 */

	    } else if (s_cmp(bndpot_1.bndtyp, "MORSE", (ftnlen)8, (ftnlen)5) 
		    == 0) {
		expterm = exp(dt * -2.);
		bde = bndpot_1.bndunit * .25 * force;
		deddt = bde * 4. * (1. - expterm) * expterm;
		d2eddt2 = bde * -8. * (1. - expterm * 2.) * expterm;
	    }

/*     scale the interaction based on its group membership */

	    if (group_1.use_group__) {
		deddt *= fgrp;
		d2eddt2 *= fgrp;
	    }

/*     set the chain rule terms for the Hessian elements */

	    if (rab2 == 0.) {
		de = 0.;
		term = 0.;
	    } else {
		de = deddt / rab;
		term = (d2eddt2 - de) / rab2;
	    }
	    termx = term * xab;
	    termy = term * yab;
	    termz = term * zab;
	    d2e_ref(1, 1) = termx * xab + de;
	    d2e_ref(1, 2) = termx * yab;
	    d2e_ref(1, 3) = termx * zab;
	    d2e_ref(2, 1) = d2e_ref(1, 2);
	    d2e_ref(2, 2) = termy * yab + de;
	    d2e_ref(2, 3) = termy * zab;
	    d2e_ref(3, 1) = d2e_ref(1, 3);
	    d2e_ref(3, 2) = d2e_ref(2, 3);
	    d2e_ref(3, 3) = termz * zab + de;

/*     increment diagonal and non-diagonal Hessian elements */

	    for (j = 1; j <= 3; ++j) {
		hessx_ref(j, ia) = hessx_ref(j, ia) + d2e_ref(1, j);
		hessy_ref(j, ia) = hessy_ref(j, ia) + d2e_ref(2, j);
		hessz_ref(j, ia) = hessz_ref(j, ia) + d2e_ref(3, j);
		hessx_ref(j, ib) = hessx_ref(j, ib) - d2e_ref(1, j);
		hessy_ref(j, ib) = hessy_ref(j, ib) - d2e_ref(2, j);
		hessz_ref(j, ib) = hessz_ref(j, ib) - d2e_ref(3, j);
	    }
	}
    }
    return 0;
} /* ebond2_ */

#undef bndlist_ref
#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef ibnd_ref
#undef d2e_ref


