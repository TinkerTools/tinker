/* replica.f -- translated by f2c (version 20050501).
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
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

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

/* Table of constant values */

static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine replica  --  periodicity via cell replication  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "replica" decides between images and replicates for generation */
/*     of periodic boundary conditions, and sets the cell replicate */
/*     list if the replicates method is to be used */


/* Subroutine */ int replica_(doublereal *cutoff)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 REPLICA  --  Truncated Octahedron\002"
	    ",\002 cannot be Replicated\002)";
    static char fmt_20[] = "(/,\002 REPLICA  --  Increase MAXCELL or Decre"
	    "ase\002,\002 the Interaction Cutoffs\002)";
    static char fmt_30[] = "(/,\002 REPLICA  --  Period Boundary via\002,i3"
	    ",\002 x\002,i3,\002 x\002,i3,\002 Set of Cell Replicates\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal maximage;
    static integer i__, j, k, nx, ny, nz;
    extern /* Subroutine */ int fatal_(void);
    static doublereal xlimit, ylimit, zlimit;

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_30, 0 };



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  bound.i  --  control of periodic boundary conditions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     polycut       cutoff distance for infinite polymer nonbonds */
/*     polycut2      square of infinite polymer nonbond cutoff */
/*     use_bounds    flag to use periodic boundary conditions */
/*     use_replica   flag to use replicates for periodic system */
/*     use_polymer   flag to mark presence of infinite polymer */




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




/*     only necessary if periodic boundaries are in use */

    if (! bound_1.use_bounds__) {
	return 0;
    }

/*     find the maximum sphere radius inscribed in periodic box */

    if (boxes_1.orthogonal) {
	xlimit = boxes_1.xbox2;
	ylimit = boxes_1.ybox2;
	zlimit = boxes_1.zbox2;
    } else if (boxes_1.monoclinic) {
	xlimit = boxes_1.xbox2 * boxes_1.beta_sin__;
	ylimit = boxes_1.ybox2;
	zlimit = boxes_1.zbox2 * boxes_1.beta_sin__;
    } else if (boxes_1.triclinic) {
	xlimit = boxes_1.xbox2 * boxes_1.beta_sin__ * boxes_1.gamma_sin__;
	ylimit = boxes_1.ybox2 * boxes_1.gamma_sin__;
	zlimit = boxes_1.zbox2 * boxes_1.beta_sin__;
    } else if (boxes_1.octahedron) {
	xlimit = sqrt(3.) / 4. * boxes_1.xbox;
	ylimit = xlimit;
	zlimit = xlimit;
    }
/* Computing MIN */
    d__1 = min(xlimit,ylimit);
    maximage = min(d__1,zlimit);

/*     use replicate method to handle cutoffs too large for images */

    if (*cutoff <= maximage) {
	bound_1.use_replica__ = FALSE_;
    } else {
	bound_1.use_replica__ = TRUE_;
    }

/*     truncated octahedron cannot use the replicates method */

    if (boxes_1.octahedron && bound_1.use_replica__) {
	io___5.ciunit = iounit_1.iout;
	s_wsfe(&io___5);
	e_wsfe();
	fatal_();
    }

/*     find the number of replicates needed based on cutoff */

    nx = (integer) (*cutoff / xlimit);
    ny = (integer) (*cutoff / ylimit);
    nz = (integer) (*cutoff / zlimit);
    if (*cutoff > (doublereal) nx * xlimit) {
	++nx;
    }
    if (*cutoff > (doublereal) ny * ylimit) {
	++ny;
    }
    if (*cutoff > (doublereal) nz * zlimit) {
	++nz;
    }
    if (nx < 1) {
	nx = 1;
    }
    if (ny < 1) {
	ny = 1;
    }
    if (nz < 1) {
	nz = 1;
    }

/*     set the replicated cell length and the half width */

    cell_1.xcell = (doublereal) nx * boxes_1.xbox;
    cell_1.ycell = (doublereal) ny * boxes_1.ybox;
    cell_1.zcell = (doublereal) nz * boxes_1.zbox;
    cell_1.xcell2 = cell_1.xcell * .5;
    cell_1.ycell2 = cell_1.ycell * .5;
    cell_1.zcell2 = cell_1.zcell * .5;

/*     check the total number of replicated unit cells */

    cell_1.ncell = nx * ny * nz - 1;
    if (cell_1.ncell > 10000) {
	io___9.ciunit = iounit_1.iout;
	s_wsfe(&io___9);
	e_wsfe();
	fatal_();
    }

/*     assign indices to the required cell replicates */

    cell_1.ncell = 0;
    i__1 = nz - 1;
    for (k = 0; k <= i__1; ++k) {
	i__2 = ny - 1;
	for (j = 0; j <= i__2; ++j) {
	    i__3 = nx - 1;
	    for (i__ = 0; i__ <= i__3; ++i__) {
		if (k != 0 || j != 0 || i__ != 0) {
		    ++cell_1.ncell;
		    icell_ref(1, cell_1.ncell) = i__;
		    icell_ref(2, cell_1.ncell) = j;
		    icell_ref(3, cell_1.ncell) = k;
		}
	    }
	}
    }

/*     print a message indicating the number of replicates used */

    if (inform_1.debug && cell_1.ncell != 0) {
	io___13.ciunit = iounit_1.iout;
	s_wsfe(&io___13);
	do_fio(&c__1, (char *)&nx, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ny, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nz, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    return 0;
} /* replica_ */

#undef icell_ref


