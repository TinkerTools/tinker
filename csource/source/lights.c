/* lights.f -- translated by f2c (version 20050501).
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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nlight, kbx[25000], kby[25000], kbz[25000], kex[25000], key[25000]
	    , kez[25000], locx[200000], locy[200000], locz[200000], rgx[
	    200000], rgy[200000], rgz[200000];
} light_;

#define light_1 light_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine lights  --  get neighbors via method of lights  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "lights" computes the set of nearest neighbor interactions */
/*     using the method of lights algorithm */

/*     literature reference: */

/*     F. Sullivan, R. D. Mountain and J. O'Connell, "Molecular */
/*     Dynamics on Vector Computers", Journal of Computational */
/*     Physics, 61, 138-153 (1985) */


/* Subroutine */ int lights_(doublereal *cutoff, integer *nsite, doublereal *
	xsort, doublereal *ysort, doublereal *zsort)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 LIGHTS  --  Number of Replicas is Too"
	    "\002,\002 Large for Method of Lights\002)";
    static char fmt_20[] = "(/,\002 LIGHTS  --  Truncated Octahedron no"
	    "t\002,\002 Supported by Method of Lights\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal box, xcut, ycut, zcut;
    extern /* Subroutine */ int sort2_(integer *, doublereal *, integer *), 
	    sort5_(integer *, integer *, integer *), fatal_(void);
    static doublereal xfrac[25000], yfrac[25000], zfrac[25000], xmove, ymove, 
	    zmove;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___2 = { 0, 0, 0, fmt_20, 0 };



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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  light.i  --  indices for method of lights pair neighbors  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     nlight  total number of sites for method of lights calculation */
/*     kbx     low index of neighbors of each site in the x-sorted list */
/*     kby     low index of neighbors of each site in the y-sorted list */
/*     kbz     low index of neighbors of each site in the z-sorted list */
/*     kex     high index of neighbors of each site in the x-sorted list */
/*     key     high index of neighbors of each site in the y-sorted list */
/*     kez     high index of neighbors of each site in the z-sorted list */
/*     locx    pointer from x-sorted list into original interaction list */
/*     locy    pointer from y-sorted list into original interaction list */
/*     locz    pointer from z-sorted list into original interaction list */
/*     rgx     pointer from original interaction list into x-sorted list */
/*     rgy     pointer from original interaction list into y-sorted list */
/*     rgz     pointer from original interaction list into z-sorted list */




/*     check that maximum number of replicates is not exceeded */

    /* Parameter adjustments */
    --zsort;
    --ysort;
    --xsort;

    /* Function Body */
    if (bound_1.use_replica__) {
	if (cell_1.xcell2 > boxes_1.xbox || cell_1.ycell2 > boxes_1.ybox || 
		cell_1.zcell2 > boxes_1.zbox) {
	    io___1.ciunit = iounit_1.iout;
	    s_wsfe(&io___1);
	    e_wsfe();
	    fatal_();
	}
    }

/*     truncated octahedron periodicity is not handled at present */

    if (bound_1.use_bounds__) {
	if (boxes_1.octahedron) {
	    io___2.ciunit = iounit_1.iout;
	    s_wsfe(&io___2);
	    e_wsfe();
	    fatal_();
	}
    }

/*     set the light width based on input distance cutoff */

    xcut = *cutoff;
    ycut = *cutoff;
    zcut = *cutoff;
    if (bound_1.use_bounds__) {
	if (boxes_1.monoclinic) {
	    zcut /= boxes_1.beta_sin__;
	    xcut += zcut * abs(boxes_1.beta_cos__);
	} else if (boxes_1.triclinic) {
	    zcut /= boxes_1.gamma_term__;
	    ycut = (ycut + zcut * abs(boxes_1.beta_term__)) / 
		    boxes_1.gamma_sin__;
	    xcut = xcut + ycut * abs(boxes_1.gamma_cos__) + zcut * abs(
		    boxes_1.beta_cos__);
	}
	xcut = min(xcut,cell_1.xcell2);
	ycut = min(ycut,cell_1.ycell2);
	zcut = min(zcut,cell_1.zcell2);
    }

/*     find fractional coordinates for the unitcell atoms */

    if (bound_1.use_bounds__) {
	if (boxes_1.orthogonal) {
	    i__1 = *nsite;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		zfrac[i__ - 1] = zsort[i__];
		yfrac[i__ - 1] = ysort[i__];
		xfrac[i__ - 1] = xsort[i__];
	    }
	} else if (boxes_1.monoclinic) {
	    i__1 = *nsite;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		zfrac[i__ - 1] = zsort[i__] / boxes_1.beta_sin__;
		yfrac[i__ - 1] = ysort[i__];
		xfrac[i__ - 1] = xsort[i__] - zfrac[i__ - 1] * 
			boxes_1.beta_cos__;
	    }
	} else if (boxes_1.triclinic) {
	    i__1 = *nsite;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		zfrac[i__ - 1] = zsort[i__] / boxes_1.gamma_term__;
		yfrac[i__ - 1] = (ysort[i__] - zfrac[i__ - 1] * 
			boxes_1.beta_term__) / boxes_1.gamma_sin__;
		xfrac[i__ - 1] = xsort[i__] - yfrac[i__ - 1] * 
			boxes_1.gamma_cos__ - zfrac[i__ - 1] * 
			boxes_1.beta_cos__;
	    }
	}
    }

/*     use images to move coordinates into periodic cell */

    if (bound_1.use_bounds__) {
	i__1 = *nsite;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xsort[i__] = xfrac[i__ - 1];
	    ysort[i__] = yfrac[i__ - 1];
	    zsort[i__] = zfrac[i__ - 1];
	    while((d__1 = xsort[i__], abs(d__1)) > cell_1.xcell2) {
		xsort[i__] -= d_sign(&cell_1.xcell, &xsort[i__]);
	    }
	    while((d__1 = ysort[i__], abs(d__1)) > cell_1.ycell2) {
		ysort[i__] -= d_sign(&cell_1.ycell, &ysort[i__]);
	    }
	    while((d__1 = zsort[i__], abs(d__1)) > cell_1.zcell2) {
		zsort[i__] -= d_sign(&cell_1.zcell, &zsort[i__]);
	    }
	}
    }

/*     generate the replica coordinates for the sort arrays */

    if (bound_1.use_replica__) {
	k = *nsite;
	i__1 = cell_1.ncell;
	for (j = 1; j <= i__1; ++j) {
	    xmove = icell_ref(1, j) * boxes_1.xbox;
	    ymove = icell_ref(2, j) * boxes_1.ybox;
	    zmove = icell_ref(3, j) * boxes_1.zbox;
	    i__2 = *nsite;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		++k;
		xsort[k] = xfrac[i__ - 1] + xmove;
		ysort[k] = yfrac[i__ - 1] + ymove;
		zsort[k] = zfrac[i__ - 1] + zmove;
		while((d__1 = xsort[k], abs(d__1)) > cell_1.xcell2) {
		    xsort[k] -= d_sign(&cell_1.xcell, &xsort[k]);
		}
		while((d__1 = ysort[k], abs(d__1)) > cell_1.ycell2) {
		    ysort[k] -= d_sign(&cell_1.ycell, &ysort[k]);
		}
		while((d__1 = zsort[k], abs(d__1)) > cell_1.zcell2) {
		    zsort[k] -= d_sign(&cell_1.zcell, &zsort[k]);
		}
	    }
	}
    }

/*     sort the coordinate components into ascending order */

    light_1.nlight = (cell_1.ncell + 1) * *nsite;
    sort2_(&light_1.nlight, &xsort[1], light_1.locx);
    sort2_(&light_1.nlight, &ysort[1], light_1.locy);
    sort2_(&light_1.nlight, &zsort[1], light_1.locz);

/*     use of replicates requires secondary sorting along x-axis */

    if (bound_1.use_replica__) {
	j = 1;
	i__1 = light_1.nlight - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (xsort[i__ + 1] != xsort[i__]) {
		i__2 = i__ - j + 1;
		sort5_(&i__2, &light_1.locx[j - 1], nsite);
		j = i__ + 1;
	    }
	}
	i__1 = light_1.nlight - j + 1;
	sort5_(&i__1, &light_1.locx[j - 1], nsite);
    }

/*     index the position of each atom in the sorted coordinates */

    i__1 = light_1.nlight;
    for (i__ = 1; i__ <= i__1; ++i__) {
	light_1.rgx[light_1.locx[i__ - 1] - 1] = i__;
	light_1.rgy[light_1.locy[i__ - 1] - 1] = i__;
	light_1.rgz[light_1.locz[i__ - 1] - 1] = i__;
    }

/*     find the negative x-coordinate boundary for each atom */

    for (i__ = light_1.nlight; i__ >= 1; --i__) {
	k = light_1.locx[i__ - 1];
	if (k <= *nsite) {
	    light_1.kbx[k - 1] = i__;
	}
    }

/*     find the positive x-coordinate boundary for each atom */

    j = 1;
    box = 0.;
    i__1 = light_1.nlight;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = light_1.locx[i__ - 1];
	if (k <= *nsite) {
	    while(xsort[j] - xsort[i__] + box < xcut) {
		if (j == light_1.nlight) {
		    if (bound_1.use_bounds__) {
			j = 0;
			box = cell_1.xcell;
		    }
		}
		++j;
		if (j > light_1.nlight) {
		    goto L30;
		}
	    }
L30:
	    --j;
	    if (j < 1) {
		j = light_1.nlight;
		box = 0.;
	    }
	    light_1.kex[k - 1] = j;
	}
    }

/*     find the negative y-coordinate boundary for each atom */

    j = light_1.nlight;
    box = 0.;
    for (i__ = light_1.nlight; i__ >= 1; --i__) {
	k = light_1.locy[i__ - 1];
	if (k <= *nsite) {
	    while(ysort[i__] - ysort[j] + box <= ycut) {
		if (j == 1) {
		    if (bound_1.use_bounds__) {
			j = light_1.nlight + 1;
			box = cell_1.ycell;
		    }
		}
		--j;
		if (j < 1) {
		    goto L40;
		}
	    }
L40:
	    ++j;
	    if (j > light_1.nlight) {
		j = 1;
		box = 0.;
	    }
	    light_1.kby[k - 1] = j;
	}
    }

/*     find the positive y-coordinate boundary for each atom */

    j = 1;
    box = 0.;
    i__1 = light_1.nlight;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = light_1.locy[i__ - 1];
	if (k <= *nsite) {
	    while(ysort[j] - ysort[i__] + box < ycut) {
		if (j == light_1.nlight) {
		    if (bound_1.use_bounds__) {
			j = 0;
			box = cell_1.ycell;
		    }
		}
		++j;
		if (j > light_1.nlight) {
		    goto L50;
		}
	    }
L50:
	    --j;
	    if (j < 1) {
		j = light_1.nlight;
		box = 0.;
	    }
	    light_1.key[k - 1] = j;
	}
    }

/*     find the negative z-coordinate boundary for each atom */

    j = light_1.nlight;
    box = 0.;
    for (i__ = light_1.nlight; i__ >= 1; --i__) {
	k = light_1.locz[i__ - 1];
	if (k <= *nsite) {
	    while(zsort[i__] - zsort[j] + box <= zcut) {
		if (j == 1) {
		    if (bound_1.use_bounds__) {
			j = light_1.nlight + 1;
			box = cell_1.zcell;
		    }
		}
		--j;
		if (j < 1) {
		    goto L60;
		}
	    }
L60:
	    ++j;
	    if (j > light_1.nlight) {
		j = 1;
		box = 0.;
	    }
	    light_1.kbz[k - 1] = j;
	}
    }

/*     find the positive z-coordinate boundary for each atom */

    j = 1;
    box = 0.;
    i__1 = light_1.nlight;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = light_1.locz[i__ - 1];
	if (k <= *nsite) {
	    while(zsort[j] - zsort[i__] + box < zcut) {
		if (j == light_1.nlight) {
		    if (bound_1.use_bounds__) {
			j = 0;
			box = cell_1.zcell;
		    }
		}
		++j;
		if (j > light_1.nlight) {
		    goto L70;
		}
	    }
L70:
	    --j;
	    if (j < 1) {
		j = light_1.nlight;
		box = 0.;
	    }
	    light_1.kez[k - 1] = j;
	}
    }
    return 0;
} /* lights_ */

#undef icell_ref


