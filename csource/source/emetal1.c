/* emetal1.f -- translated by f2c (version 20050501).
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
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal desum[75000]	/* was [3][25000] */, deb[75000]	/* 
	    was [3][25000] */, dea[75000]	/* was [3][25000] */, deba[
	    75000]	/* was [3][25000] */, deub[75000]	/* was [3][
	    25000] */, deaa[75000]	/* was [3][25000] */, deopb[75000]	
	    /* was [3][25000] */, deopd[75000]	/* was [3][25000] */, deid[
	    75000]	/* was [3][25000] */, deit[75000]	/* was [3][
	    25000] */, det[75000]	/* was [3][25000] */, dept[75000]	
	    /* was [3][25000] */, debt[75000]	/* was [3][25000] */, dett[
	    75000]	/* was [3][25000] */, dev[75000]	/* was [3][
	    25000] */, dec[75000]	/* was [3][25000] */, decd[75000]	
	    /* was [3][25000] */, ded[75000]	/* was [3][25000] */, dem[
	    75000]	/* was [3][25000] */, dep[75000]	/* was [3][
	    25000] */, der[75000]	/* was [3][25000] */, des[75000]	
	    /* was [3][25000] */, delf[75000]	/* was [3][25000] */, deg[
	    75000]	/* was [3][25000] */, dex[75000]	/* was [3][
	    25000] */;
} deriv_;

#define deriv_1 deriv_

struct {
    doublereal esum, eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, 
	    ebt, ett, ev, ec, ecd, ed, em, ep, er, es, elf, eg, ex;
} energi_;

#define energi_1 energi_

struct {
    doublereal chg[5000];
} kchrge_;

#define kchrge_1 kchrge_



/*     ################################################################## */
/*     ##  COPYRIGHT (C) 2001 by Anders Carlsson & Jay William Ponder  ## */
/*     ##                     All Rights Reserved                      ## */
/*     ################################################################## */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine emetal1  --  ligand field energy & derivatives  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "emetal1" calculates the transition metal ligand field energy */
/*     and its first derivatives with respect to Cartesian coordinates */


/* Subroutine */ int emetal1_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal dedrback[30]	/* was [3][10] */;
    static integer neighnum[10];
    static doublereal argument, e;
    static integer i__, j, k;
    static doublereal e0;
    static integer k0;
    static doublereal r2, de[30]	/* was [3][10] */;
    static integer jj;
    static doublereal cij, dot, elf0, r0cu, beta, rcut, xmet, ymet, zmet, 
	    rback[30]	/* was [3][10] */, alpha;
    extern /* Subroutine */ int fatal_(void);
    static doublereal kappa, demet[3], rback2, angfac[100]	/* was [10][
	    10] */, delfdh[10];
    static integer ineigh, jneigh;
    static doublereal expfac[10];
    static integer nneigh;
    static doublereal rneigh[10], xneigh[10], yneigh[10], zneigh[10], cosmat[
	    100]	/* was [10][10] */, esqtet, facback[10], dangfac[100]	
	    /* was [10][10] */, detback[30]	/* was [3][10] */, dftback[30]
	    	/* was [3][10] */, delfrad[10], dfacback[30]	/* was [3][10]
	     */;


#define dedrback_ref(a_1,a_2) dedrback[(a_2)*3 + a_1 - 4]
#define de_ref(a_1,a_2) de[(a_2)*3 + a_1 - 4]
#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define delf_ref(a_1,a_2) deriv_1.delf[(a_2)*3 + a_1 - 4]
#define rback_ref(a_1,a_2) rback[(a_2)*3 + a_1 - 4]
#define angfac_ref(a_1,a_2) angfac[(a_2)*10 + a_1 - 11]
#define cosmat_ref(a_1,a_2) cosmat[(a_2)*10 + a_1 - 11]
#define dangfac_ref(a_1,a_2) dangfac[(a_2)*10 + a_1 - 11]
#define detback_ref(a_1,a_2) detback[(a_2)*3 + a_1 - 4]
#define dftback_ref(a_1,a_2) dftback[(a_2)*3 + a_1 - 4]
#define dfacback_ref(a_1,a_2) dfacback[(a_2)*3 + a_1 - 4]



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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  deriv.i  --  Cartesian coordinate derivative components  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     desum   total energy Cartesian coordinate derivatives */
/*     deb     bond stretch Cartesian coordinate derivatives */
/*     dea     angle bend Cartesian coordinate derivatives */
/*     deba    stretch-bend Cartesian coordinate derivatives */
/*     deub    Urey-Bradley Cartesian coordinate derivatives */
/*     deaa    angle-angle Cartesian coordinate derivatives */
/*     deopb   out-of-plane bend Cartesian coordinate derivatives */
/*     deopd   out-of-plane distance Cartesian coordinate derivatives */
/*     deid    improper dihedral Cartesian coordinate derivatives */
/*     deit    improper torsion Cartesian coordinate derivatives */
/*     det     torsional Cartesian coordinate derivatives */
/*     dept    pi-orbital torsion Cartesian coordinate derivatives */
/*     debt    stretch-torsion Cartesian coordinate derivatives */
/*     dett    torsion-torsion Cartesian coordinate derivatives */
/*     dev     van der Waals Cartesian coordinate derivatives */
/*     dec     charge-charge Cartesian coordinate derivatives */
/*     decd    charge-dipole Cartesian coordinate derivatives */
/*     ded     dipole-dipole Cartesian coordinate derivatives */
/*     dem     multipole Cartesian coordinate derivatives */
/*     dep     polarization Cartesian coordinate derivatives */
/*     der     reaction field Cartesian coordinate derivatives */
/*     des     solvation Cartesian coordinate derivatives */
/*     delf    metal ligand field Cartesian coordinate derivatives */
/*     deg     geometric restraint Cartesian coordinate derivatives */
/*     dex     extra energy term Cartesian coordinate derivatives */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  kchrge.i  --  forcefield parameters for partial charges  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     chg   partial charge parameters for each atom type */




/*     zero out metal ligand field energy and first derivatives */

    energi_1.elf = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	delf_ref(1, i__) = 0.;
	delf_ref(2, i__) = 0.;
	delf_ref(3, i__) = 0.;
    }

/*     begin ligand field calculation; for now, only for Cu+2 */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (atmtyp_1.atomic[i__ - 1] != 29 || kchrge_1.chg[atoms_1.type__[i__ 
		- 1] - 1] != 2.) {
	    goto L30;
	}
	nneigh = 0;
	xmet = atoms_1.x[i__ - 1];
	ymet = atoms_1.y[i__ - 1];
	zmet = atoms_1.z__[i__ - 1];
	i__2 = atoms_1.n;
	for (j = 1; j <= i__2; ++j) {
	    if (j == i__) {
		goto L10;
	    }
	    if (atmtyp_1.atomic[j - 1] != 7 && atmtyp_1.atomic[j - 1] != 8) {
		goto L10;
	    }

/*     next are standard distance, decay factor and splitting energy */

	    r0cu = 2.06;
	    kappa = 1.;
	    esqtet = 1.64;

/*     semiclassical method obtains only 65% of sq-tet difference */

	    elf0 = esqtet * 1.78 / .78;
	    elf0 *= .65;
	    e0 = elf0 * 23.05 / 2.608;
	    rcut = 2.5;
	    alpha = 0.;
	    beta = -1.;
/* Computing 2nd power */
	    d__1 = atoms_1.x[i__ - 1] - atoms_1.x[j - 1];
/* Computing 2nd power */
	    d__2 = atoms_1.y[i__ - 1] - atoms_1.y[j - 1];
/* Computing 2nd power */
	    d__3 = atoms_1.z__[i__ - 1] - atoms_1.z__[j - 1];
	    r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
/* Computing 2nd power */
	    d__1 = rcut;
	    if (r2 > d__1 * d__1) {
		goto L10;
	    }
	    ++nneigh;
	    xneigh[nneigh - 1] = atoms_1.x[j - 1];
	    yneigh[nneigh - 1] = atoms_1.y[j - 1];
	    zneigh[nneigh - 1] = atoms_1.z__[j - 1];
	    neighnum[nneigh - 1] = j;
	    if (couple_1.n12[j - 1] <= 0) {
		fatal_();
	    }

/*     set average of back neighbors relative to j, */
/*     notice that it is defined as a difference */

	    rback_ref(1, nneigh) = 0.;
	    rback_ref(2, nneigh) = 0.;
	    rback_ref(3, nneigh) = 0.;
	    i__3 = couple_1.n12[j - 1];
	    for (k0 = 1; k0 <= i__3; ++k0) {
		k = i12_ref(k0, j);
		rback_ref(1, nneigh) = rback_ref(1, nneigh) + (atoms_1.x[k - 
			1] - atoms_1.x[j - 1]) / (doublereal) couple_1.n12[j 
			- 1];
		rback_ref(2, nneigh) = rback_ref(2, nneigh) + (atoms_1.y[k - 
			1] - atoms_1.y[j - 1]) / (doublereal) couple_1.n12[j 
			- 1];
		rback_ref(3, nneigh) = rback_ref(3, nneigh) + (atoms_1.z__[k 
			- 1] - atoms_1.z__[j - 1]) / (doublereal) 
			couple_1.n12[j - 1];
	    }
	    facback[nneigh - 1] = 0.;
	    dot = rback_ref(1, nneigh) * (atoms_1.x[i__ - 1] - atoms_1.x[j - 
		    1]) + rback_ref(2, nneigh) * (atoms_1.y[i__ - 1] - 
		    atoms_1.y[j - 1]) + rback_ref(3, nneigh) * (atoms_1.z__[
		    i__ - 1] - atoms_1.z__[j - 1]);
/* Computing 2nd power */
	    d__1 = rback_ref(1, nneigh);
/* Computing 2nd power */
	    d__2 = rback_ref(2, nneigh);
/* Computing 2nd power */
	    d__3 = rback_ref(3, nneigh);
	    rback2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
/* Computing 2nd power */
	    d__1 = alpha * sqrt(rback2) * sqrt(r2) + beta * dot;
	    facback[nneigh - 1] = d__1 * d__1 / (rback2 * r2);

/*     dfacback is derivative of back factor with respect to center */
/*     of gravity of back atoms; detback is corresponding energy */

	    dfacback_ref(1, nneigh) = (alpha * sqrt(r2) * rback_ref(1, nneigh)
		     / sqrt(rback2) + beta * (atoms_1.x[i__ - 1] - atoms_1.x[
		    j - 1])) * 2. * (alpha * sqrt(rback2 * r2) + beta * dot) /
		     (rback2 * r2) - facback[nneigh - 1] * 2. * rback_ref(1, 
		    nneigh) / rback2;
	    dfacback_ref(2, nneigh) = (alpha * sqrt(r2) * rback_ref(2, nneigh)
		     / sqrt(rback2) + beta * (atoms_1.y[i__ - 1] - atoms_1.y[
		    j - 1])) * 2. * (alpha * sqrt(rback2 * r2) + beta * dot) /
		     (rback2 * r2) - facback[nneigh - 1] * 2. * rback_ref(2, 
		    nneigh) / rback2;
	    dfacback_ref(3, nneigh) = (alpha * sqrt(r2) * rback_ref(3, nneigh)
		     / sqrt(rback2) + beta * (atoms_1.z__[i__ - 1] - 
		    atoms_1.z__[j - 1])) * 2. * (alpha * sqrt(rback2 * r2) + 
		    beta * dot) / (rback2 * r2) - facback[nneigh - 1] * 2. * 
		    rback_ref(3, nneigh) / rback2;
	    dftback_ref(1, nneigh) = (alpha * sqrt(rback2) * (atoms_1.x[i__ - 
		    1] - atoms_1.x[j - 1]) / sqrt(r2) + beta * rback_ref(1, 
		    nneigh)) * 2. * (alpha * sqrt(rback2 * r2) + beta * dot) /
		     (rback2 * r2) - facback[nneigh - 1] * 2. * (atoms_1.x[
		    i__ - 1] - atoms_1.x[j - 1]) / r2;
	    dftback_ref(2, nneigh) = (alpha * sqrt(rback2) * (atoms_1.y[i__ - 
		    1] - atoms_1.y[j - 1]) / sqrt(r2) + beta * rback_ref(2, 
		    nneigh)) * 2. * (alpha * sqrt(rback2 * r2) + beta * dot) /
		     (rback2 * r2) - facback[nneigh - 1] * 2. * (atoms_1.y[
		    i__ - 1] - atoms_1.y[j - 1]) / r2;
	    dftback_ref(3, nneigh) = (alpha * sqrt(rback2) * (atoms_1.z__[i__ 
		    - 1] - atoms_1.z__[j - 1]) / sqrt(r2) + beta * rback_ref(
		    3, nneigh)) * 2. * (alpha * sqrt(rback2 * r2) + beta * 
		    dot) / (rback2 * r2) - facback[nneigh - 1] * 2. * (
		    atoms_1.z__[i__ - 1] - atoms_1.z__[j - 1]) / r2;
L10:
	    ;
	}

/*     calculate the energy and derivatives for current interaction */

	i__2 = nneigh;
	for (ineigh = 1; ineigh <= i__2; ++ineigh) {
/* Computing 2nd power */
	    d__1 = xneigh[ineigh - 1] - xmet;
/* Computing 2nd power */
	    d__2 = yneigh[ineigh - 1] - ymet;
/* Computing 2nd power */
	    d__3 = zneigh[ineigh - 1] - zmet;
	    rneigh[ineigh - 1] = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	    rneigh[ineigh - 1] = sqrt(rneigh[ineigh - 1]);
	    expfac[ineigh - 1] = exp(-kappa * (rneigh[ineigh - 1] - r0cu));
	    jj = neighnum[ineigh - 1];
	    if (atmtyp_1.atomic[jj - 1] == 8) {
		expfac[ineigh - 1] *= .4;
	    }
	}
	i__2 = nneigh;
	for (ineigh = 1; ineigh <= i__2; ++ineigh) {
	    jj = neighnum[ineigh - 1];
	}
	i__2 = nneigh;
	for (ineigh = 1; ineigh <= i__2; ++ineigh) {
	    i__3 = nneigh;
	    for (jneigh = 1; jneigh <= i__3; ++jneigh) {
		dot = (xneigh[ineigh - 1] - xmet) * (xneigh[jneigh - 1] - 
			xmet) + (yneigh[ineigh - 1] - ymet) * (yneigh[jneigh 
			- 1] - ymet) + (zneigh[ineigh - 1] - zmet) * (zneigh[
			jneigh - 1] - zmet);
		cij = dot / (rneigh[ineigh - 1] * rneigh[jneigh - 1]);
		cosmat_ref(ineigh, jneigh) = cij;
/* Computing 4th power */
		d__1 = cij, d__1 *= d__1;
/* Computing 2nd power */
		d__2 = cij;
		angfac_ref(ineigh, jneigh) = d__1 * d__1 * 2.25 - d__2 * d__2 
			* 1.5 + .005;
/* Computing 3rd power */
		d__1 = cij;
		dangfac_ref(ineigh, jneigh) = d__1 * (d__1 * d__1) * 9. - cij 
			* 3.;
	    }
	}
	argument = 0.;
	i__2 = nneigh;
	for (ineigh = 1; ineigh <= i__2; ++ineigh) {
	    i__3 = nneigh;
	    for (jneigh = 1; jneigh <= i__3; ++jneigh) {
		argument += expfac[ineigh - 1] * expfac[jneigh - 1] * 
			angfac_ref(ineigh, jneigh) * facback[ineigh - 1] * 
			facback[jneigh - 1];
	    }
	}
	e = 0.;
	i__2 = nneigh;
	for (ineigh = 1; ineigh <= i__2; ++ineigh) {
	    de_ref(1, ineigh) = 0.;
	    de_ref(2, ineigh) = 0.;
	    de_ref(3, ineigh) = 0.;
	}
	if (argument <= 0.) {
	    goto L20;
	}
	if (argument > 0.) {
	    e = -e0 * sqrt(argument);
	}

/*     set up radial derivatives of energy */

	i__2 = nneigh;
	for (ineigh = 1; ineigh <= i__2; ++ineigh) {
	    delfrad[ineigh - 1] = 0.;
	    i__3 = nneigh;
	    for (jneigh = 1; jneigh <= i__3; ++jneigh) {
		delfrad[ineigh - 1] += expfac[jneigh - 1] * angfac_ref(ineigh,
			 jneigh) * facback[jneigh - 1];
	    }
	    delfdh[ineigh - 1] = delfrad[ineigh - 1] * (e0 / e);

/*     note two minus signs above cancel */

/* Computing 2nd power */
	    d__1 = e0;
	    delfrad[ineigh - 1] = -delfrad[ineigh - 1] * (d__1 * d__1 / e) * 
		    kappa * expfac[ineigh - 1] * facback[ineigh - 1];

/*     note cancelling factors of two from square root and product */

	}

/*     below does angular derivatives */

	i__2 = nneigh;
	for (ineigh = 1; ineigh <= i__2; ++ineigh) {
	    de_ref(1, ineigh) = 0.;
	    de_ref(2, ineigh) = 0.;
	    de_ref(3, ineigh) = 0.;
	    i__3 = nneigh;
	    for (jneigh = 1; jneigh <= i__3; ++jneigh) {
/* Computing 2nd power */
		d__1 = rneigh[ineigh - 1];
		de_ref(1, ineigh) = de_ref(1, ineigh) + expfac[jneigh - 1] * 
			facback[jneigh - 1] * dangfac_ref(ineigh, jneigh) * ((
			xneigh[jneigh - 1] - xmet) / (rneigh[ineigh - 1] * 
			rneigh[jneigh - 1]) - (xneigh[ineigh - 1] - xmet) * 
			cosmat_ref(ineigh, jneigh) / (d__1 * d__1));
/* Computing 2nd power */
		d__1 = rneigh[ineigh - 1];
		de_ref(2, ineigh) = de_ref(2, ineigh) + expfac[jneigh - 1] * 
			facback[jneigh - 1] * dangfac_ref(ineigh, jneigh) * ((
			yneigh[jneigh - 1] - ymet) / (rneigh[ineigh - 1] * 
			rneigh[jneigh - 1]) - (yneigh[ineigh - 1] - ymet) * 
			cosmat_ref(ineigh, jneigh) / (d__1 * d__1));
/* Computing 2nd power */
		d__1 = rneigh[ineigh - 1];
		de_ref(3, ineigh) = de_ref(3, ineigh) + expfac[jneigh - 1] * 
			facback[jneigh - 1] * dangfac_ref(ineigh, jneigh) * ((
			zneigh[jneigh - 1] - zmet) / (rneigh[ineigh - 1] * 
			rneigh[jneigh - 1]) - (zneigh[ineigh - 1] - zmet) * 
			cosmat_ref(ineigh, jneigh) / (d__1 * d__1));
	    }
	    de_ref(1, ineigh) = de_ref(1, ineigh) * e0 * e0 * expfac[ineigh - 
		    1] * facback[ineigh - 1] / e;
	    de_ref(2, ineigh) = de_ref(2, ineigh) * e0 * e0 * expfac[ineigh - 
		    1] * facback[ineigh - 1] / e;
	    de_ref(3, ineigh) = de_ref(3, ineigh) * e0 * e0 * expfac[ineigh - 
		    1] * facback[ineigh - 1] / e;
	}
	i__2 = nneigh;
	for (ineigh = 1; ineigh <= i__2; ++ineigh) {
	    de_ref(1, ineigh) = de_ref(1, ineigh) + delfrad[ineigh - 1] * (
		    xneigh[ineigh - 1] - xmet) / rneigh[ineigh - 1];
	    de_ref(2, ineigh) = de_ref(2, ineigh) + delfrad[ineigh - 1] * (
		    yneigh[ineigh - 1] - ymet) / rneigh[ineigh - 1];
	    de_ref(3, ineigh) = de_ref(3, ineigh) + delfrad[ineigh - 1] * (
		    zneigh[ineigh - 1] - zmet) / rneigh[ineigh - 1];
	}
	i__2 = nneigh;
	for (ineigh = 1; ineigh <= i__2; ++ineigh) {
	    dedrback_ref(1, ineigh) = dfacback_ref(1, ineigh) * e0 * expfac[
		    ineigh - 1] * delfdh[ineigh - 1];
	    dedrback_ref(2, ineigh) = dfacback_ref(2, ineigh) * e0 * expfac[
		    ineigh - 1] * delfdh[ineigh - 1];
	    dedrback_ref(3, ineigh) = dfacback_ref(3, ineigh) * e0 * expfac[
		    ineigh - 1] * delfdh[ineigh - 1];
	    detback_ref(1, ineigh) = dftback_ref(1, ineigh) * e0 * expfac[
		    ineigh - 1] * delfdh[ineigh - 1];
	    detback_ref(2, ineigh) = dftback_ref(2, ineigh) * e0 * expfac[
		    ineigh - 1] * delfdh[ineigh - 1];
	    detback_ref(3, ineigh) = dftback_ref(3, ineigh) * e0 * expfac[
		    ineigh - 1] * delfdh[ineigh - 1];
	}
L20:
	demet[0] = 0.;
	demet[1] = 0.;
	demet[2] = 0.;
	i__2 = nneigh;
	for (ineigh = 1; ineigh <= i__2; ++ineigh) {
	    demet[0] = demet[0] - de_ref(1, ineigh) + detback_ref(1, ineigh);
	    demet[1] = demet[1] - de_ref(2, ineigh) + detback_ref(2, ineigh);
	    demet[2] = demet[2] - de_ref(3, ineigh) + detback_ref(3, ineigh);
	}
	energi_1.elf += e;
	delf_ref(1, i__) = delf_ref(1, i__) + demet[0];
	delf_ref(2, i__) = delf_ref(2, i__) + demet[1];
	delf_ref(3, i__) = delf_ref(3, i__) + demet[2];
	i__2 = nneigh;
	for (ineigh = 1; ineigh <= i__2; ++ineigh) {
	    j = neighnum[ineigh - 1];
	    delf_ref(1, j) = delf_ref(1, j) + de_ref(1, ineigh) - 
		    dedrback_ref(1, ineigh) - detback_ref(1, ineigh);
	    delf_ref(2, j) = delf_ref(2, j) + de_ref(2, ineigh) - 
		    dedrback_ref(2, ineigh) - detback_ref(2, ineigh);
	    delf_ref(3, j) = delf_ref(3, j) + de_ref(3, ineigh) - 
		    dedrback_ref(3, ineigh) - detback_ref(3, ineigh);
	}
	i__2 = nneigh;
	for (ineigh = 1; ineigh <= i__2; ++ineigh) {
	    j = neighnum[ineigh - 1];
	    if (couple_1.n12[j - 1] <= 0) {
		fatal_();
	    }
	    i__3 = couple_1.n12[j - 1];
	    for (k0 = 1; k0 <= i__3; ++k0) {
		k = i12_ref(k0, j);
		delf_ref(1, k) = delf_ref(1, k) + dedrback_ref(1, ineigh) / (
			doublereal) couple_1.n12[j - 1];
		delf_ref(2, k) = delf_ref(2, k) + dedrback_ref(2, ineigh) / (
			doublereal) couple_1.n12[j - 1];
		delf_ref(3, k) = delf_ref(3, k) + dedrback_ref(3, ineigh) / (
			doublereal) couple_1.n12[j - 1];
	    }
	}
L30:
	;
    }
    return 0;
} /* emetal1_ */

#undef dfacback_ref
#undef dftback_ref
#undef detback_ref
#undef dangfac_ref
#undef cosmat_ref
#undef angfac_ref
#undef rback_ref
#undef delf_ref
#undef i12_ref
#undef de_ref
#undef dedrback_ref


