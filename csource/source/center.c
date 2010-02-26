/* center.f -- translated by f2c (version 20050501).
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
    doublereal wfit[25000];
    integer nfit, ifit[50000]	/* was [2][25000] */;
} align_;

#define align_1 align_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine center  --  superimpose structure centroids  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "center" moves the weighted centroid of each coordinate */
/*     set to the origin during least squares superposition */


/* Subroutine */ int center_(integer *n1, doublereal *x1, doublereal *y1, 
	doublereal *z1, integer *n2, doublereal *x2, doublereal *y2, 
	doublereal *z2, doublereal *xmid, doublereal *ymid, doublereal *zmid)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k;
    static doublereal norm, weigh;


#define ifit_ref(a_1,a_2) align_1.ifit[(a_2)*2 + a_1 - 3]



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  align.i  --  information for superposition of structures  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     wfit    weights assigned to atom pairs during superposition */
/*     nfit    number of atoms to use in superimposing two structures */
/*     ifit    atom numbers of pairs of atoms to be superimposed */




/*     find the weighted centroid of the second */
/*     structure and translate it to the origin */

    /* Parameter adjustments */
    --z2;
    --y2;
    --x2;
    --z1;
    --y1;
    --x1;

    /* Function Body */
    *xmid = 0.;
    *ymid = 0.;
    *zmid = 0.;
    norm = 0.;
    i__1 = align_1.nfit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = ifit_ref(2, i__);
	weigh = align_1.wfit[i__ - 1];
	*xmid += x2[k] * weigh;
	*ymid += y2[k] * weigh;
	*zmid += z2[k] * weigh;
	norm += weigh;
    }
    *xmid /= norm;
    *ymid /= norm;
    *zmid /= norm;
    i__1 = *n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x2[i__] -= *xmid;
	y2[i__] -= *ymid;
	z2[i__] -= *zmid;
    }

/*     now repeat for the first structure, note */
/*     that this centroid position gets returned */

    *xmid = 0.;
    *ymid = 0.;
    *zmid = 0.;
    norm = 0.;
    i__1 = align_1.nfit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = ifit_ref(1, i__);
	weigh = align_1.wfit[i__ - 1];
	*xmid += x1[k] * weigh;
	*ymid += y1[k] * weigh;
	*zmid += z1[k] * weigh;
	norm += weigh;
    }
    *xmid /= norm;
    *ymid /= norm;
    *zmid /= norm;
    i__1 = *n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x1[i__] -= *xmid;
	y1[i__] -= *ymid;
	z1[i__] -= *zmid;
    }
    return 0;
} /* center_ */

#undef ifit_ref


