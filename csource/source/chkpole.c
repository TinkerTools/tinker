/* chkpole.f -- translated by f2c (version 20050501).
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
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine chkpole  --  check multipoles at chiral sites  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "chkpole" inverts atomic multipole moments as necessary */
/*     at sites with chiral local reference frame definitions */


/* Subroutine */ int chkpole_(void)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, k;
    static doublereal c1, c2, c3;
    static integer ia, ib, ic, id;
    static doublereal xad, yad, zad, xbd, ybd, zbd, xcd, ycd, zcd, vol;
    static logical check;


#define pole_ref(a_1,a_2) mpole_1.pole[(a_2)*13 + a_1 - 14]
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




/*     loop over multipole sites testing for chirality inversion */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	check = TRUE_;
	if (s_cmp(polaxe_ref(0, i__), "Z-then-X", (ftnlen)8, (ftnlen)8) != 0) 
		{
	    check = FALSE_;
	}
	if (mpole_1.yaxis[i__ - 1] == 0) {
	    check = FALSE_;
	}
	if (check) {
	    k = mpole_1.yaxis[i__ - 1];
	    ia = mpole_1.ipole[i__ - 1];
	    ib = mpole_1.zaxis[i__ - 1];
	    ic = mpole_1.xaxis[i__ - 1];
	    id = abs(k);

/*     compute the signed parallelpiped volume at chiral site */

	    xad = atoms_1.x[ia - 1] - atoms_1.x[id - 1];
	    yad = atoms_1.y[ia - 1] - atoms_1.y[id - 1];
	    zad = atoms_1.z__[ia - 1] - atoms_1.z__[id - 1];
	    xbd = atoms_1.x[ib - 1] - atoms_1.x[id - 1];
	    ybd = atoms_1.y[ib - 1] - atoms_1.y[id - 1];
	    zbd = atoms_1.z__[ib - 1] - atoms_1.z__[id - 1];
	    xcd = atoms_1.x[ic - 1] - atoms_1.x[id - 1];
	    ycd = atoms_1.y[ic - 1] - atoms_1.y[id - 1];
	    zcd = atoms_1.z__[ic - 1] - atoms_1.z__[id - 1];
	    c1 = ybd * zcd - zbd * ycd;
	    c2 = ycd * zad - zcd * yad;
	    c3 = yad * zbd - zad * ybd;
	    vol = xad * c1 + xbd * c2 + xcd * c3;

/*     invert atomic multipole components involving the y-axis */

	    if (k < 0 && vol > 0. || k > 0 && vol < 0.) {
		mpole_1.yaxis[i__ - 1] = -k;
		pole_ref(3, i__) = -pole_ref(3, i__);
		pole_ref(6, i__) = -pole_ref(6, i__);
		pole_ref(8, i__) = -pole_ref(8, i__);
		pole_ref(10, i__) = -pole_ref(10, i__);
		pole_ref(12, i__) = -pole_ref(12, i__);
	    }
	}
    }
    return 0;
} /* chkpole_ */

#undef polaxe_ref
#undef pole_ref


