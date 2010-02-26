/* lattice.f -- translated by f2c (version 20050501).
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

struct {
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine lattice  --  setup periodic boundary conditions  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "lattice" stores the periodic box dimensions and sets angle */
/*     values to be used in computing fractional coordinates */


/* Subroutine */ int lattice_(void)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal alpha_cos__, ar1, ar2, ar3, br1, br2, br3, cr1, cr2, 
	    cr3;


#define lvec_ref(a_1,a_2) boxes_1.lvec[(a_2)*3 + a_1 - 4]
#define recip_ref(a_1,a_2) boxes_1.recip[(a_2)*3 + a_1 - 4]



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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  inform.i  --  control values for I/O and program flow  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     digits    decimal places output for energy and coordinates */
/*     iprint    steps between status printing (0=no printing) */
/*     iwrite    steps between coordinate dumps (0=no dumps) */
/*     isend     steps between socket communication (0=no sockets) */
/*     verbose   logical flag to turn on extra information */
/*     debug     logical flag to turn on full debug printing */
/*     holdup    logical flag to wait for carriage return on exit */
/*     abort     logical flag to stop execution at next chance */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




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




/*     compute and store the half box length values */

    boxes_1.xbox2 = boxes_1.xbox * .5;
    boxes_1.ybox2 = boxes_1.ybox * .5;
    boxes_1.zbox2 = boxes_1.zbox * .5;
    if (boxes_1.octahedron) {
	boxes_1.box34 = boxes_1.xbox * .75;
    }

/*     set replicated cell dimensions equal to the unitcell */

    cell_1.xcell = boxes_1.xbox;
    cell_1.ycell = boxes_1.ybox;
    cell_1.zcell = boxes_1.zbox;
    cell_1.xcell2 = boxes_1.xbox2;
    cell_1.ycell2 = boxes_1.ybox2;
    cell_1.zcell2 = boxes_1.zbox2;

/*     get values needed for fractional coordinate computations */

    if (boxes_1.orthogonal || boxes_1.octahedron) {
	alpha_cos__ = 0.;
	boxes_1.beta_sin__ = 1.;
	boxes_1.beta_cos__ = 0.;
	boxes_1.gamma_sin__ = 1.;
	boxes_1.gamma_cos__ = 0.;
	boxes_1.beta_term__ = 0.;
	boxes_1.gamma_term__ = 1.;
    } else if (boxes_1.monoclinic) {
	alpha_cos__ = 0.;
	boxes_1.beta_sin__ = sin(boxes_1.beta / 57.29577951308232088);
	boxes_1.beta_cos__ = cos(boxes_1.beta / 57.29577951308232088);
	boxes_1.gamma_sin__ = 1.;
	boxes_1.gamma_cos__ = 0.;
	boxes_1.beta_term__ = 0.;
	boxes_1.gamma_term__ = boxes_1.beta_sin__;
    } else if (boxes_1.triclinic) {
	alpha_cos__ = cos(boxes_1.alpha / 57.29577951308232088);
	boxes_1.beta_sin__ = sin(boxes_1.beta / 57.29577951308232088);
	boxes_1.beta_cos__ = cos(boxes_1.beta / 57.29577951308232088);
	boxes_1.gamma_sin__ = sin(boxes_1.gamma / 57.29577951308232088);
	boxes_1.gamma_cos__ = cos(boxes_1.gamma / 57.29577951308232088);
	boxes_1.beta_term__ = (alpha_cos__ - boxes_1.beta_cos__ * 
		boxes_1.gamma_cos__) / boxes_1.gamma_sin__;
/* Computing 2nd power */
	d__1 = boxes_1.beta_sin__;
/* Computing 2nd power */
	d__2 = boxes_1.beta_term__;
	boxes_1.gamma_term__ = sqrt(d__1 * d__1 - d__2 * d__2);
    }

/*     determine the volume of the parent periodic box */

    boxes_1.volbox = 0.;
    if (boxes_1.orthogonal || boxes_1.octahedron) {
	boxes_1.volbox = boxes_1.xbox * boxes_1.ybox * boxes_1.zbox;
    } else if (boxes_1.monoclinic) {
	boxes_1.volbox = boxes_1.beta_sin__ * boxes_1.xbox * boxes_1.ybox * 
		boxes_1.zbox;
    } else if (boxes_1.triclinic) {
	boxes_1.volbox = boxes_1.gamma_sin__ * boxes_1.gamma_term__ * 
		boxes_1.xbox * boxes_1.ybox * boxes_1.zbox;
    }

/*     compute and store real space lattice vectors as rows */

    ar1 = boxes_1.xbox;
    ar2 = 0.;
    ar3 = 0.;
    br1 = boxes_1.ybox * boxes_1.gamma_cos__;
    br2 = boxes_1.ybox * boxes_1.gamma_sin__;
    br3 = 0.;
    cr1 = boxes_1.zbox * boxes_1.beta_cos__;
    cr2 = boxes_1.zbox * boxes_1.beta_term__;
    cr3 = boxes_1.zbox * boxes_1.gamma_term__;
    lvec_ref(1, 1) = ar1;
    lvec_ref(1, 2) = ar2;
    lvec_ref(1, 3) = ar3;
    lvec_ref(2, 1) = br1;
    lvec_ref(2, 2) = br2;
    lvec_ref(2, 3) = br3;
    lvec_ref(3, 1) = cr1;
    lvec_ref(3, 2) = cr2;
    lvec_ref(3, 3) = cr3;

/*     compute and store reciprocal lattice vectors as columns */

    if (boxes_1.volbox != 0.) {
	recip_ref(1, 1) = (br2 * cr3 - cr2 * br3) / boxes_1.volbox;
	recip_ref(2, 1) = (br3 * cr1 - cr3 * br1) / boxes_1.volbox;
	recip_ref(3, 1) = (br1 * cr2 - cr1 * br2) / boxes_1.volbox;
	recip_ref(1, 2) = (cr2 * ar3 - ar2 * cr3) / boxes_1.volbox;
	recip_ref(2, 2) = (cr3 * ar1 - ar3 * cr1) / boxes_1.volbox;
	recip_ref(3, 2) = (cr1 * ar2 - ar1 * cr2) / boxes_1.volbox;
	recip_ref(1, 3) = (ar2 * br3 - br2 * ar3) / boxes_1.volbox;
	recip_ref(2, 3) = (ar3 * br1 - br3 * ar1) / boxes_1.volbox;
	recip_ref(3, 3) = (ar1 * br2 - br1 * ar2) / boxes_1.volbox;
    }

/*     volume of truncated octahedron is half of cubic parent */

    if (boxes_1.octahedron) {
	boxes_1.volbox *= .5;
    }
    return 0;
} /* lattice_ */

#undef recip_ref
#undef lvec_ref


