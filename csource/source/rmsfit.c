/* rmsfit.f -- translated by f2c (version 20050501).
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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  function rmsfit  --  rms deviation for paired atoms  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "rmsfit" computes the rms fit of two coordinate sets */


doublereal rmsfit_(doublereal *x1, doublereal *y1, doublereal *z1, doublereal 
	*x2, doublereal *y2, doublereal *z2)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, i1, i2;
    static doublereal xr, yr, zr, norm, dist2, weigh, rmsterm;


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




/*     compute the rms fit over superimposed atom pairs */

    /* Parameter adjustments */
    --z2;
    --y2;
    --x2;
    --z1;
    --y1;
    --x1;

    /* Function Body */
    ret_val = 0.;
    norm = 0.;
    i__1 = align_1.nfit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = ifit_ref(1, i__);
	i2 = ifit_ref(2, i__);
	weigh = align_1.wfit[i__ - 1];
	xr = x1[i1] - x2[i2];
	yr = y1[i1] - y2[i2];
	zr = z1[i1] - z2[i2];
/* Computing 2nd power */
	d__1 = xr;
/* Computing 2nd power */
	d__2 = yr;
/* Computing 2nd power */
	d__3 = zr;
	dist2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	norm += weigh;
	rmsterm = dist2 * weigh;
	ret_val += rmsterm;
    }
    ret_val = sqrt(ret_val / norm);
    return ret_val;
} /* rmsfit_ */

#undef ifit_ref


