/* image.f -- translated by f2c (version 20050501).
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
    doublereal xbox, ybox, zbox, alpha, beta, gamma, xbox2, ybox2, zbox2, 
	    box34, lvec[9]	/* was [3][3] */, recip[9]	/* was [3][3] 
	    */, volbox, beta_sin__, beta_cos__, gamma_sin__, gamma_cos__, 
	    beta_term__, gamma_term__;
    logical orthogonal, monoclinic, triclinic, octahedron;
    char spacegrp[10];
} boxes_;

#define boxes_1 boxes_

struct {
    doublereal xcell, ycell, zcell, xcell2, ycell2, zcell2;
    integer ncell, icell[30000]	/* was [3][10000] */;
} cell_;

#define cell_1 cell_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine image  --  pairwise minimum image distance  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "image" takes the components of pairwise distance between */
/*     two points in a periodic box and converts to the components */
/*     of the minimum image distance */


/* Subroutine */ int image_(doublereal *xr, doublereal *yr, doublereal *zr)
{
    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal xf, yf, zf;



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cell.i  --  periodic boundaries using replicated cells  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xcell    length of the a-axis of the complete replicated cell */
/*     ycell    length of the b-axis of the complete replicated cell */
/*     zcell    length of the c-axis of the complete replicated cell */
/*     xcell2   half the length of the a-axis of the replicated cell */
/*     ycell2   half the length of the b-axis of the replicated cell */
/*     zcell2   half the length of the c-axis of the replicated cell */
/*     ncell    total number of cell replicates for periodic boundaries */
/*     icell    offset along axes for each replicate periodic cell */




/*     for orthogonal lattice, find the desired image directly */

    if (boxes_1.orthogonal) {
	while(abs(*xr) > cell_1.xcell2) {
	    *xr -= d_sign(&cell_1.xcell, xr);
	}
	while(abs(*yr) > cell_1.ycell2) {
	    *yr -= d_sign(&cell_1.ycell, yr);
	}
	while(abs(*zr) > cell_1.zcell2) {
	    *zr -= d_sign(&cell_1.zcell, zr);
	}

/*     for monoclinic lattice, convert "xr" and "zr" to */
/*     fractional coordinates, find desired image and then */
/*     translate fractional coordinates back to Cartesian */

    } else if (boxes_1.monoclinic) {
	zf = *zr / boxes_1.beta_sin__;
	xf = *xr - zf * boxes_1.beta_cos__;
	while(abs(xf) > cell_1.xcell2) {
	    xf -= d_sign(&cell_1.xcell, &xf);
	}
	while(abs(*yr) > cell_1.ycell2) {
	    *yr -= d_sign(&cell_1.ycell, yr);
	}
	while(abs(zf) > cell_1.zcell2) {
	    zf -= d_sign(&cell_1.zcell, &zf);
	}
	*xr = xf + zf * boxes_1.beta_cos__;
	*zr = zf * boxes_1.beta_sin__;

/*     for triclinic lattice, convert pairwise components to */
/*     fractional coordinates, find desired image and then */
/*     translate fractional coordinates back to Cartesian */

    } else if (boxes_1.triclinic) {
	zf = *zr / boxes_1.gamma_term__;
	yf = (*yr - zf * boxes_1.beta_term__) / boxes_1.gamma_sin__;
	xf = *xr - yf * boxes_1.gamma_cos__ - zf * boxes_1.beta_cos__;
	while(abs(xf) > cell_1.xcell2) {
	    xf -= d_sign(&cell_1.xcell, &xf);
	}
	while(abs(yf) > cell_1.ycell2) {
	    yf -= d_sign(&cell_1.ycell, &yf);
	}
	while(abs(zf) > cell_1.zcell2) {
	    zf -= d_sign(&cell_1.zcell, &zf);
	}
	*xr = xf + yf * boxes_1.gamma_cos__ + zf * boxes_1.beta_cos__;
	*yr = yf * boxes_1.gamma_sin__ + zf * boxes_1.beta_term__;
	*zr = zf * boxes_1.gamma_term__;

/*     for truncated octahedron, use orthogonal box equations, */
/*     then perform extra tests to remove corner pieces */

    } else if (boxes_1.octahedron) {
	while(abs(*xr) > boxes_1.xbox2) {
	    *xr -= d_sign(&boxes_1.xbox, xr);
	}
	while(abs(*yr) > boxes_1.ybox2) {
	    *yr -= d_sign(&boxes_1.ybox, yr);
	}
	while(abs(*zr) > boxes_1.zbox2) {
	    *zr -= d_sign(&boxes_1.zbox, zr);
	}
	if (abs(*xr) + abs(*yr) + abs(*zr) > boxes_1.box34) {
	    *xr -= d_sign(&boxes_1.xbox2, xr);
	    *yr -= d_sign(&boxes_1.ybox2, yr);
	    *zr -= d_sign(&boxes_1.zbox2, zr);
	}
    }
    return 0;
} /* image_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine imager  --  replicate minimum image distance  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "imager" takes the components of pairwise distance between */
/*     two points in the same or neighboring periodic boxes and */
/*     converts to the components of the minimum image distance */


/* Subroutine */ int imager_(doublereal *xr, doublereal *yr, doublereal *zr, 
	integer *i__)
{
    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal xf, yf, zf, xmove, ymove, zmove, xsize, ysize, zsize, 
	    xsize2, ysize2, zsize2;


#define icell_ref(a_1,a_2) cell_1.icell[(a_2)*3 + a_1 - 4]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cell.i  --  periodic boundaries using replicated cells  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xcell    length of the a-axis of the complete replicated cell */
/*     ycell    length of the b-axis of the complete replicated cell */
/*     zcell    length of the c-axis of the complete replicated cell */
/*     xcell2   half the length of the a-axis of the replicated cell */
/*     ycell2   half the length of the b-axis of the replicated cell */
/*     zcell2   half the length of the c-axis of the replicated cell */
/*     ncell    total number of cell replicates for periodic boundaries */
/*     icell    offset along axes for each replicate periodic cell */




/*     set dimensions for either single box or replicated cell */

    if (*i__ >= 0) {
	xsize = cell_1.xcell;
	ysize = cell_1.ycell;
	zsize = cell_1.zcell;
	xsize2 = cell_1.xcell2;
	ysize2 = cell_1.ycell2;
	zsize2 = cell_1.zcell2;
    } else {
	xsize = boxes_1.xbox;
	ysize = boxes_1.ybox;
	zsize = boxes_1.zbox;
	xsize2 = boxes_1.xbox2;
	ysize2 = boxes_1.ybox2;
	zsize2 = boxes_1.zbox2;
    }

/*     compute the distance to translate along each cell axis */

    if (*i__ <= 0) {
	xmove = 0.;
	ymove = 0.;
	zmove = 0.;
    } else {
	xmove = icell_ref(1, *i__) * boxes_1.xbox;
	ymove = icell_ref(2, *i__) * boxes_1.ybox;
	zmove = icell_ref(3, *i__) * boxes_1.zbox;
    }

/*     for orthogonal lattice, find the desired image directly */

    if (boxes_1.orthogonal) {
	*xr += xmove;
	while(abs(*xr) > xsize2) {
	    *xr -= d_sign(&xsize, xr);
	}
	*yr += ymove;
	while(abs(*yr) > ysize2) {
	    *yr -= d_sign(&ysize, yr);
	}
	*zr += zmove;
	while(abs(*zr) > zsize2) {
	    *zr -= d_sign(&zsize, zr);
	}

/*     for monoclinic lattice, convert "xr" and "zr" to */
/*     fractional coordinates, find desired image and then */
/*     translate fractional coordinates back to Cartesian */

    } else if (boxes_1.monoclinic) {
	zf = *zr / boxes_1.beta_sin__;
	xf = *xr - zf * boxes_1.beta_cos__;
	xf += xmove;
	while(abs(xf) > xsize2) {
	    xf -= d_sign(&xsize, &xf);
	}
	*yr += ymove;
	while(abs(*yr) > ysize2) {
	    *yr -= d_sign(&ysize, yr);
	}
	zf += zmove;
	while(abs(zf) > zsize2) {
	    zf -= d_sign(&zsize, &zf);
	}
	*xr = xf + zf * boxes_1.beta_cos__;
	*zr = zf * boxes_1.beta_sin__;

/*     for triclinic lattice, convert pairwise components to */
/*     fractional coordinates, find desired image and then */
/*     translate fractional coordinates back to Cartesian */

    } else if (boxes_1.triclinic) {
	zf = *zr / boxes_1.gamma_term__;
	yf = (*yr - zf * boxes_1.beta_term__) / boxes_1.gamma_sin__;
	xf = *xr - yf * boxes_1.gamma_cos__ - zf * boxes_1.beta_cos__;
	xf += xmove;
	while(abs(xf) > xsize2) {
	    xf -= d_sign(&xsize, &xf);
	}
	yf += ymove;
	while(abs(yf) > ysize2) {
	    yf -= d_sign(&ysize, &yf);
	}
	zf += zmove;
	while(abs(zf) > zsize2) {
	    zf -= d_sign(&zsize, &zf);
	}
	*xr = xf + yf * boxes_1.gamma_cos__ + zf * boxes_1.beta_cos__;
	*yr = yf * boxes_1.gamma_sin__ + zf * boxes_1.beta_term__;
	*zr = zf * boxes_1.gamma_term__;

/*     for truncated octahedron, use orthogonal box equations, */
/*     then perform extra tests to remove corner pieces */

    } else if (boxes_1.octahedron) {
	while(abs(*xr) > boxes_1.xbox2) {
	    *xr -= d_sign(&boxes_1.xbox, xr);
	}
	while(abs(*yr) > boxes_1.ybox2) {
	    *yr -= d_sign(&boxes_1.ybox, yr);
	}
	while(abs(*zr) > boxes_1.zbox2) {
	    *zr -= d_sign(&boxes_1.zbox, zr);
	}
	if (abs(*xr) + abs(*yr) + abs(*zr) > boxes_1.box34) {
	    *xr -= d_sign(&boxes_1.xbox2, xr);
	    *yr -= d_sign(&boxes_1.ybox2, yr);
	    *zr -= d_sign(&boxes_1.zbox2, zr);
	}
    }
    return 0;
} /* imager_ */

#undef icell_ref


