/* bounds.f -- translated by f2c (version 20050501).
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
    doublereal xbox, ybox, zbox, alpha, beta, gamma, xbox2, ybox2, zbox2, 
	    box34, lvec[9]	/* was [3][3] */, recip[9]	/* was [3][3] 
	    */, volbox, beta_sin__, beta_cos__, gamma_sin__, gamma_cos__, 
	    beta_term__, gamma_term__;
    logical orthogonal, monoclinic, triclinic, octahedron;
    char spacegrp[10];
} boxes_;

#define boxes_1 boxes_

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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine bounds  --  check periodic boundary conditions  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "bounds" finds the center of mass of each molecule and */
/*     translates any stray molecules back into the periodic box */


/* Subroutine */ int bounds_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal xmid;
    static integer init;
    static doublereal ymid, zmid, xoff[25000], yoff[25000], zoff[25000];
    static integer stop;
    static doublereal weigh, xfrac, yfrac, zfrac;


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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  boxes.i  --  parameters for periodic boundary conditions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     xbox        length of a-axis of periodic box in Angstroms */
/*     ybox        length of b-axis of periodic box in Angstroms */
/*     zbox        length of c-axis of periodic box in Angstroms */
/*     alpha       angle between b- and c-axes of box in degrees */
/*     beta        angle between a- and c-axes of box in degrees */
/*     gamma       angle between a- and b-axes of box in degrees */
/*     xbox2       half of the a-axis length of periodic box */
/*     ybox2       half of the b-axis length of periodic box */
/*     zbox2       half of the c-axis length of periodic box */
/*     box34       three-fourths axis length of truncated octahedron */
/*     lvec        real space lattice vectors as matrix rows */
/*     recip       reciprocal lattice vectors as matrix columns */
/*     volbox      volume in Ang**3 of the periodic box */
/*     beta_sin    sine of the beta periodic box angle */
/*     beta_cos    cosine of the beta periodic box angle */
/*     gamma_sin   sine of the gamma periodic box angle */
/*     gamma_cos   cosine of the gamma periodic box angle */
/*     beta_term   term used in generating triclinic box */
/*     gamma_term  term used in generating triclinic box */
/*     orthogonal  flag to mark periodic box as orthogonal */
/*     monoclinic  flag to mark periodic box as monoclinic */
/*     triclinic   flag to mark periodic box as triclinic */
/*     octahedron  flag to mark box as truncated octahedron */
/*     spacegrp    space group symbol for the unitcell type */




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




/*     locate the center of mass of each molecule */

    i__1 = molcul_1.nmol;
    for (i__ = 1; i__ <= i__1; ++i__) {
	init = imol_ref(1, i__);
	stop = imol_ref(2, i__);
	xmid = 0.;
	ymid = 0.;
	zmid = 0.;
	i__2 = stop;
	for (j = init; j <= i__2; ++j) {
	    k = molcul_1.kmol[j - 1];
	    weigh = atmtyp_1.mass[k - 1];
	    xmid += atoms_1.x[k - 1] * weigh;
	    ymid += atoms_1.y[k - 1] * weigh;
	    zmid += atoms_1.z__[k - 1] * weigh;
	}
	weigh = molcul_1.molmass[i__ - 1];
	xmid /= weigh;
	ymid /= weigh;
	zmid /= weigh;

/*     save atomic coordinates relative to center of mass */

	i__2 = stop;
	for (j = init; j <= i__2; ++j) {
	    k = molcul_1.kmol[j - 1];
	    xoff[k - 1] = atoms_1.x[k - 1] - xmid;
	    yoff[k - 1] = atoms_1.y[k - 1] - ymid;
	    zoff[k - 1] = atoms_1.z__[k - 1] - zmid;
	}

/*     get fractional coordinates of center of mass */

	if (boxes_1.orthogonal || boxes_1.octahedron) {
	    zfrac = zmid;
	    yfrac = ymid;
	    xfrac = xmid;
	} else if (boxes_1.monoclinic) {
	    zfrac = zmid / boxes_1.beta_sin__;
	    yfrac = ymid;
	    xfrac = xmid - zfrac * boxes_1.beta_cos__;
	} else if (boxes_1.triclinic) {
	    zfrac = zmid / boxes_1.gamma_term__;
	    yfrac = (ymid - zfrac * boxes_1.beta_term__) / 
		    boxes_1.gamma_sin__;
	    xfrac = xmid - yfrac * boxes_1.gamma_cos__ - zfrac * 
		    boxes_1.beta_cos__;
	}

/*     translate center of mass into the periodic box */

	while(xfrac > boxes_1.xbox2) {
	    xfrac -= boxes_1.xbox;
	}
	while(xfrac < -boxes_1.xbox2) {
	    xfrac += boxes_1.xbox;
	}
	while(yfrac > boxes_1.ybox2) {
	    yfrac -= boxes_1.ybox;
	}
	while(yfrac < -boxes_1.ybox2) {
	    yfrac += boxes_1.ybox;
	}
	while(zfrac > boxes_1.zbox2) {
	    zfrac -= boxes_1.zbox;
	}
	while(zfrac < -boxes_1.zbox2) {
	    zfrac += boxes_1.zbox;
	}

/*     truncated octahedron needs to have corners removed */

	if (boxes_1.octahedron) {
	    if (abs(xfrac) + abs(yfrac) + abs(zfrac) > boxes_1.box34) {
		xfrac -= d_sign(&boxes_1.xbox2, &xfrac);
		yfrac -= d_sign(&boxes_1.ybox2, &yfrac);
		zfrac -= d_sign(&boxes_1.zbox2, &zfrac);
	    }
	}

/*     convert fractional center of mass back to Cartesian */

	if (boxes_1.orthogonal || boxes_1.octahedron) {
	    xmid = xfrac;
	    ymid = yfrac;
	    zmid = zfrac;
	} else if (boxes_1.monoclinic) {
	    xmid = xfrac + zfrac * boxes_1.beta_cos__;
	    ymid = yfrac;
	    zmid = zfrac * boxes_1.beta_sin__;
	} else if (boxes_1.triclinic) {
	    xmid = xfrac + yfrac * boxes_1.gamma_cos__ + zfrac * 
		    boxes_1.beta_cos__;
	    ymid = yfrac * boxes_1.gamma_sin__ + zfrac * boxes_1.beta_term__;
	    zmid = zfrac * boxes_1.gamma_term__;
	}

/*     translate coordinates via offset from center of mass */

	i__2 = stop;
	for (j = init; j <= i__2; ++j) {
	    k = molcul_1.kmol[j - 1];
	    atoms_1.x[k - 1] = xoff[k - 1] + xmid;
	    atoms_1.y[k - 1] = yoff[k - 1] + ymid;
	    atoms_1.z__[k - 1] = zoff[k - 1] + zmid;
	}
    }
    return 0;
} /* bounds_ */

#undef imol_ref


