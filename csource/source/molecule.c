/* molecule.f -- translated by f2c (version 20050501).
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
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  subroutine molecule  --  assign atoms to molecules  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     "molecule" counts the molecules, assigns each atom to */
/*     its molecule and computes the mass of each molecule */


/* Subroutine */ int molecule_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, mi, mj, mk, list[25000];
    extern /* Subroutine */ int sort_(integer *, integer *), sort3_(integer *,
	     integer *, integer *);
    static integer iattach;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define imol_ref(a_1,a_2) molcul_1.imol[(a_2)*2 + a_1 - 3]



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  molcul.i  --  individual molecules within current system  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     molmass   molecular weight for each molecule in the system */
/*     totmass   total weight of all the molecules in the system */
/*     nmol      total number of separate molecules in the system */
/*     kmol      contiguous list of the atoms in each molecule */
/*     imol      first and last atom of each molecule in the list */
/*     molcule   number of the molecule to which each atom belongs */




/*     zero number of molecules and molecule membership list */

    molcul_1.nmol = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	molcul_1.molcule[i__ - 1] = 0;
    }

/*     assign each atom to its respective molecule */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (molcul_1.molcule[i__ - 1] == 0) {
	    ++molcul_1.nmol;
	    molcul_1.molcule[i__ - 1] = molcul_1.nmol;
	}
	mi = molcul_1.molcule[i__ - 1];
	i__2 = couple_1.n12[i__ - 1];
	for (iattach = 1; iattach <= i__2; ++iattach) {
	    j = i12_ref(iattach, i__);
	    mj = molcul_1.molcule[j - 1];
	    if (mj == 0) {
		molcul_1.molcule[j - 1] = mi;
	    } else if (mi < mj) {
		--molcul_1.nmol;
		i__3 = atoms_1.n;
		for (k = 1; k <= i__3; ++k) {
		    mk = molcul_1.molcule[k - 1];
		    if (mk == mj) {
			molcul_1.molcule[k - 1] = mi;
		    } else if (mk > mj) {
			molcul_1.molcule[k - 1] = mk - 1;
		    }
		}
	    } else if (mi > mj) {
		--molcul_1.nmol;
		i__3 = atoms_1.n;
		for (k = 1; k <= i__3; ++k) {
		    mk = molcul_1.molcule[k - 1];
		    if (mk == mi) {
			molcul_1.molcule[k - 1] = mj;
		    } else if (mk > mi) {
			molcul_1.molcule[k - 1] = mk - 1;
		    }
		}
		mi = mj;
	    }
	}
    }

/*     pack atoms of each molecule into a contiguous indexed list */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	list[i__ - 1] = molcul_1.molcule[i__ - 1];
    }
    sort3_(&atoms_1.n, list, molcul_1.kmol);

/*     find the first and last atom in each molecule */

    k = 1;
    imol_ref(1, 1) = 1;
    i__1 = atoms_1.n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	j = list[i__ - 1];
	if (j != k) {
	    imol_ref(2, k) = i__ - 1;
	    k = j;
	    imol_ref(1, k) = i__;
	}
    }
    imol_ref(2, molcul_1.nmol) = atoms_1.n;

/*     sort the list of atoms in each molecule by atom number */

    i__1 = molcul_1.nmol;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = imol_ref(2, i__) - imol_ref(1, i__) + 1;
	sort_(&k, &molcul_1.kmol[imol_ref(1, i__) - 1]);
    }

/*     if all atomic masses are zero, set them all to unity */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (atmtyp_1.mass[i__ - 1] != 0.) {
	    goto L10;
	}
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atmtyp_1.mass[i__ - 1] = 1.;
    }
L10:

/*     compute the mass of each molecule and the total mass */

    molcul_1.totmass = 0.;
    i__1 = molcul_1.nmol;
    for (i__ = 1; i__ <= i__1; ++i__) {
	molcul_1.molmass[i__ - 1] = 0.;
	i__2 = imol_ref(2, i__);
	for (k = imol_ref(1, i__); k <= i__2; ++k) {
	    molcul_1.molmass[i__ - 1] += atmtyp_1.mass[molcul_1.kmol[k - 1] - 
		    1];
	}
	molcul_1.totmass += molcul_1.molmass[i__ - 1];
    }
    return 0;
} /* molecule_ */

#undef imol_ref
#undef i12_ref


