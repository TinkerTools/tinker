/* numgrad.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine numgrad  --  numerical gradient of a function  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "numgrad" computes the gradient of the objective function */
/*     "fvalue" with respect to Cartesian coordinates of the atoms */
/*     via a one-sided or two-sided numerical differentiation */


/* Subroutine */ int numgrad_(D_fp fvalue, doublereal *g, doublereal *eps)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static logical twosided;
    static doublereal f;
    static integer i__;
    static doublereal f0, old;


#define g_ref(a_1,a_2) g[(a_2)*3 + a_1]



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




/*     chose between use of one-sided or two-sided gradient */

    /* Parameter adjustments */
    g -= 4;

    /* Function Body */
    twosided = TRUE_;
    if (! twosided) {
	f0 = (*fvalue)();
    }

/*     compute the numerical gradient from function values */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	old = atoms_1.x[i__ - 1];
	if (twosided) {
	    atoms_1.x[i__ - 1] -= *eps * .5;
	    f0 = (*fvalue)();
	}
	atoms_1.x[i__ - 1] += *eps;
	f = (*fvalue)();
	atoms_1.x[i__ - 1] = old;
	g_ref(1, i__) = (f - f0) / *eps;
	old = atoms_1.y[i__ - 1];
	if (twosided) {
	    atoms_1.y[i__ - 1] -= *eps * .5;
	    f0 = (*fvalue)();
	}
	atoms_1.y[i__ - 1] += *eps;
	f = (*fvalue)();
	atoms_1.y[i__ - 1] = old;
	g_ref(2, i__) = (f - f0) / *eps;
	old = atoms_1.z__[i__ - 1];
	if (twosided) {
	    atoms_1.z__[i__ - 1] -= *eps * .5;
	    f0 = (*fvalue)();
	}
	atoms_1.z__[i__ - 1] += *eps;
	f = (*fvalue)();
	atoms_1.z__[i__ - 1] = old;
	g_ref(3, i__) = (f - f0) / *eps;
    }
    return 0;
} /* numgrad_ */

#undef g_ref


