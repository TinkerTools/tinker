/* quatfit.f -- translated by f2c (version 20050501).
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

/* Table of constant values */

static integer c__4 = 4;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine quatfit  --  quaternion superposition of coords  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "quatfit" uses a quaternion-based method to achieve the best */
/*     fit superposition of two sets of coordinates */

/*     literature reference: */

/*     S. K. Kearsley, "On the Orthogonal Transformation Used for */
/*     Structural Comparisons", Acta Crystallographica Section A, */
/*     45, 208-210 (1989) */

/*     adapted from an original program written by D. J. Heisterberg, */
/*     Ohio Supercomputer Center, Columbus, OH */


/* Subroutine */ int quatfit_(integer *n1, doublereal *x1, doublereal *y1, 
	doublereal *z1, integer *n2, doublereal *x2, doublereal *y2, 
	doublereal *z2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static doublereal c__[16]	/* was [4][4] */, d__[4];
    static integer i__;
    static doublereal q[4], v[16]	/* was [4][4] */;
    static integer i1, i2;
    static doublereal rot[9]	/* was [3][3] */, xrot, yrot, zrot, xxyx, 
	    xxyy, xxyz, xyyx, xyyy, xyyz, xzyx, xzyy, xzyz, work1[4], work2[4]
	    , weigh;
    extern /* Subroutine */ int jacobi_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


#define c___ref(a_1,a_2) c__[(a_2)*4 + a_1 - 5]
#define v_ref(a_1,a_2) v[(a_2)*4 + a_1 - 5]
#define rot_ref(a_1,a_2) rot[(a_2)*3 + a_1 - 4]
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




/*     build the upper triangle of the quadratic form matrix */

    /* Parameter adjustments */
    --z2;
    --y2;
    --x2;
    --z1;
    --y1;
    --x1;

    /* Function Body */
    xxyx = 0.;
    xxyy = 0.;
    xxyz = 0.;
    xyyx = 0.;
    xyyy = 0.;
    xyyz = 0.;
    xzyx = 0.;
    xzyy = 0.;
    xzyz = 0.;
    i__1 = align_1.nfit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = ifit_ref(1, i__);
	i2 = ifit_ref(2, i__);
	weigh = align_1.wfit[i__ - 1];
	xxyx += weigh * x1[i1] * x2[i2];
	xxyy += weigh * y1[i1] * x2[i2];
	xxyz += weigh * z1[i1] * x2[i2];
	xyyx += weigh * x1[i1] * y2[i2];
	xyyy += weigh * y1[i1] * y2[i2];
	xyyz += weigh * z1[i1] * y2[i2];
	xzyx += weigh * x1[i1] * z2[i2];
	xzyy += weigh * y1[i1] * z2[i2];
	xzyz += weigh * z1[i1] * z2[i2];
    }
    c___ref(1, 1) = xxyx + xyyy + xzyz;
    c___ref(1, 2) = xzyy - xyyz;
    c___ref(2, 2) = xxyx - xyyy - xzyz;
    c___ref(1, 3) = xxyz - xzyx;
    c___ref(2, 3) = xxyy + xyyx;
    c___ref(3, 3) = xyyy - xzyz - xxyx;
    c___ref(1, 4) = xyyx - xxyy;
    c___ref(2, 4) = xzyx + xxyz;
    c___ref(3, 4) = xyyz + xzyy;
    c___ref(4, 4) = xzyz - xxyx - xyyy;

/*     diagonalize the quadratic form matrix */

    jacobi_(&c__4, &c__4, c__, d__, v, work1, work2);

/*     extract the desired quaternion */

    q[0] = v_ref(1, 4);
    q[1] = v_ref(2, 4);
    q[2] = v_ref(3, 4);
    q[3] = v_ref(4, 4);

/*     assemble rotation matrix that superimposes the molecules */

/* Computing 2nd power */
    d__1 = q[0];
/* Computing 2nd power */
    d__2 = q[1];
/* Computing 2nd power */
    d__3 = q[2];
/* Computing 2nd power */
    d__4 = q[3];
    rot_ref(1, 1) = d__1 * d__1 + d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
    rot_ref(2, 1) = (q[1] * q[2] - q[0] * q[3]) * 2.;
    rot_ref(3, 1) = (q[1] * q[3] + q[0] * q[2]) * 2.;
    rot_ref(1, 2) = (q[2] * q[1] + q[0] * q[3]) * 2.;
/* Computing 2nd power */
    d__1 = q[0];
/* Computing 2nd power */
    d__2 = q[1];
/* Computing 2nd power */
    d__3 = q[2];
/* Computing 2nd power */
    d__4 = q[3];
    rot_ref(2, 2) = d__1 * d__1 - d__2 * d__2 + d__3 * d__3 - d__4 * d__4;
    rot_ref(3, 2) = (q[2] * q[3] - q[0] * q[1]) * 2.;
    rot_ref(1, 3) = (q[3] * q[1] - q[0] * q[2]) * 2.;
    rot_ref(2, 3) = (q[3] * q[2] + q[0] * q[1]) * 2.;
/* Computing 2nd power */
    d__1 = q[0];
/* Computing 2nd power */
    d__2 = q[1];
/* Computing 2nd power */
    d__3 = q[2];
/* Computing 2nd power */
    d__4 = q[3];
    rot_ref(3, 3) = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 + d__4 * d__4;

/*     rotate second molecule to best fit with first molecule */

    i__1 = *n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xrot = x2[i__] * rot_ref(1, 1) + y2[i__] * rot_ref(1, 2) + z2[i__] * 
		rot_ref(1, 3);
	yrot = x2[i__] * rot_ref(2, 1) + y2[i__] * rot_ref(2, 2) + z2[i__] * 
		rot_ref(2, 3);
	zrot = x2[i__] * rot_ref(3, 1) + y2[i__] * rot_ref(3, 2) + z2[i__] * 
		rot_ref(3, 3);
	x2[i__] = xrot;
	y2[i__] = yrot;
	z2[i__] = zrot;
    }
    return 0;
} /* quatfit_ */

#undef ifit_ref
#undef rot_ref
#undef v_ref
#undef c___ref


