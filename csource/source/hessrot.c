/* hessrot.f -- translated by f2c (version 20050501).
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
    doublereal dihed[1000];
    integer nomega, iomega[2000]	/* was [2][1000] */, zline[1000];
} omega_;

#define omega_1 omega_

struct {
    doublereal zbond[25000], zang[25000], ztors[25000];
    integer iz[100000]	/* was [4][25000] */;
} zcoord_;

#define zcoord_1 zcoord_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  subroutine hessrot  --  torsional Hessian elements  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     "hessrot" computes the numerical Hessian elements with */
/*     respect to torsional angles; either the full matrix or */
/*     just the diagonal can be calculated; the full matrix */
/*     needs nomega+1 gradient evaluations while the diagonal */
/*     requires just two gradient calls */


/* Subroutine */ int hessrot_(char *mode, doublereal *hrot, ftnlen mode_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal e, g[1000];
    static integer i__, j;
    static doublereal g0[1000], old[1000], eps;
    static integer line;
    extern /* Subroutine */ int gradrot_(doublereal *, doublereal *), 
	    makexyz_(void);


#define hrot_ref(a_1,a_2) hrot[(a_2)*1000 + a_1]



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
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  math.i  --  mathematical and geometrical constants  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     radian   conversion factor from radians to degrees */
/*     pi       numerical value of the geometric constant */
/*     sqrtpi   numerical value of the square root of Pi */
/*     logten   numerical value of the natural log of ten */
/*     sqrttwo  numerical value of the square root of two */
/*     twosix   numerical value of the sixth root of two */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     zbond   bond length used to define each Z-matrix atom */
/*     zang    bond angle used to define each Z-matrix atom */
/*     ztors   angle or torsion used to define Z-matrix atom */
/*     iz      defining atom numbers for each Z-matrix atom */




/*     calculate base values for the torsional gradient */

    /* Parameter adjustments */
    hrot -= 1001;

    /* Function Body */
    eps = 1e-4;
    gradrot_(&e, g0);

/*     compute one-sided numerical Hessian from gradient values; */
/*     set off-diagonal elements to the average symmetric value */

    if (s_cmp(mode, "FULL", (ftnlen)4, (ftnlen)4) == 0) {
	i__1 = omega_1.nomega;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    line = omega_1.zline[i__ - 1];
	    old[i__ - 1] = zcoord_1.ztors[line - 1];
	    zcoord_1.ztors[line - 1] += eps * 57.29577951308232088;
	    makexyz_();
	    gradrot_(&e, g);
	    zcoord_1.ztors[line - 1] = old[i__ - 1];
	    i__2 = omega_1.nomega;
	    for (j = 1; j <= i__2; ++j) {
		hrot_ref(j, i__) = (g[j - 1] - g0[j - 1]) / eps;
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		hrot_ref(j, i__) = (hrot_ref(j, i__) + hrot_ref(i__, j)) * .5;
		hrot_ref(i__, j) = hrot_ref(j, i__);
	    }
	}

/*     compute numerical Hessian diagonal from gradient values */

    } else if (s_cmp(mode, "DIAG", (ftnlen)4, (ftnlen)4) == 0) {
	i__1 = omega_1.nomega;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    line = omega_1.zline[i__ - 1];
	    old[i__ - 1] = zcoord_1.ztors[line - 1];
	    zcoord_1.ztors[line - 1] += eps * 57.29577951308232088;
	}
	makexyz_();
	gradrot_(&e, g);
	i__1 = omega_1.nomega;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    hrot_ref(i__, i__) = (g[i__ - 1] - g0[i__ - 1]) / eps;
	    line = omega_1.zline[i__ - 1];
	    zcoord_1.ztors[line - 1] = old[i__ - 1];
	}
    }

/*     restore the Cartesian coordinates to original values */

    makexyz_();
    return 0;
} /* hessrot_ */

#undef hrot_ref


