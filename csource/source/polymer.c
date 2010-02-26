/* polymer.f -- translated by f2c (version 20050501).
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
    doublereal bk[50000], bl[50000];
    integer nbond, ibnd[100000]	/* was [2][50000] */;
} bond_;

#define bond_1 bond_

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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine polymer  --  check for an infinite polymer  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "polymer" tests for the presence of an infinite polymer */
/*     extending across periodic boundaries */


/* Subroutine */ int polymer_(void)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 POLYMER  --  Image Conflicts for Infin"
	    "ite\002,\002 Polymer in Small Cell\002)";
    static char fmt_40[] = "(/,\002 POLYMER  --  Warning, Infinite Polyme"
	    "r\002,\002 Cutoff may be Too Small\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static doublereal maximage;
    static integer i__, j, ia, ib;
    static doublereal xr, yr, zr, xab, yab, zab, eps;
    static integer next;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *), fatal_(void);
    static doublereal delta;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static doublereal xlimit, ylimit, zlimit;
    static char string[120], keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___6 = { 1, string, 1, 0, 120, 1 };
    static cilist io___22 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_40, 0 };



#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define keyline_ref(a_0,a_1) &keys_1.keyline[(a_1)*120 + a_0 - 120]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  bond.i  --  covalent bonds in the current structure  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     bk      bond stretch force constants (kcal/mole/Ang**2) */
/*     bl      ideal bond length values in Angstroms */
/*     nbond   total number of bond stretches in the system */
/*     ibnd    numbers of the atoms in each bond stretch */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  keys.i  --  contents of current keyword parameter file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nkey      number of nonblank lines in the keyword file */
/*     keyline   contents of each individual keyword file line */




/*     set defaults of infinite polymer usage and cutoff distance */

    bound_1.use_polymer__ = FALSE_;
    bound_1.polycut = 5.5;

/*     get any keywords containing infinite polymer cutoff parameters */

    i__1 = keys_1.nkey;
    for (j = 1; j <= i__1; ++j) {
	next = 1;
	s_copy(record, keyline_ref(0, j), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "POLYMER-CUTOFF ", (ftnlen)15, (ftnlen)15) == 0) {
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___6);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&bound_1.polycut, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
L10:
	    ;
	}
    }

/*     see if any bond connections require a minimum image */

    if (bound_1.use_bounds__) {
	eps = 1e-4;
	i__1 = bond_1.nbond;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = ibnd_ref(1, i__);
	    ib = ibnd_ref(2, i__);
	    xab = atoms_1.x[ia - 1] - atoms_1.x[ib - 1];
	    yab = atoms_1.y[ia - 1] - atoms_1.y[ib - 1];
	    zab = atoms_1.z__[ia - 1] - atoms_1.z__[ib - 1];
	    xr = xab;
	    yr = yab;
	    zr = zab;
	    image_(&xr, &yr, &zr);
	    delta = (d__1 = xr - xab, abs(d__1)) + (d__2 = yr - yab, abs(d__2)
		    ) + (d__3 = zr - zab, abs(d__3));
	    if (delta > eps) {
		bound_1.use_polymer__ = TRUE_;
		goto L20;
	    }
	}
L20:
	;
    }

/*     find the maximum sphere radius inscribed in periodic box */

    if (bound_1.use_polymer__) {
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

/*     check for too large or small an infinite polymer cutoff */

	if (bound_1.polycut > maximage) {
	    io___22.ciunit = iounit_1.iout;
	    s_wsfe(&io___22);
	    e_wsfe();
	    fatal_();
	} else if (bound_1.polycut < 5.5) {
	    io___23.ciunit = iounit_1.iout;
	    s_wsfe(&io___23);
	    e_wsfe();
	}
    }

/*     set square of cutoff distance for use with nonbonded terms */

    bound_1.polycut2 = bound_1.polycut * bound_1.polycut;
    return 0;
} /* polymer_ */

#undef keyline_ref
#undef ibnd_ref


