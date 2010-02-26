/* chkring.f -- translated by f2c (version 20050501).
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
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  2004  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine chkring  --  check atom set for small rings  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "chkring" tests an atom or a set of connected atoms for */
/*     their presence within a single 3- to 6-membered ring */


/* Subroutine */ int chkring_(integer *iring, integer *ia, integer *ib, 
	integer *ic, integer *id)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, m, p, q, r__, nset;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]



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
/*     ##  couple.i  --  near-neighbor atom connectivity lists   ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxn13   maximum number of atoms 1-3 connected to an atom */
/*     maxn14   maximum number of atoms 1-4 connected to an atom */
/*     maxn15   maximum number of atoms 1-5 connected to an atom */

/*     n12      number of atoms directly bonded to each atom */
/*     i12      atom numbers of atoms 1-2 connected to each atom */
/*     n13      number of atoms in a 1-3 relation to each atom */
/*     i13      atom numbers of atoms 1-3 connected to each atom */
/*     n14      number of atoms in a 1-4 relation to each atom */
/*     i14      atom numbers of atoms 1-4 connected to each atom */
/*     n15      number of atoms in a 1-5 relation to each atom */
/*     i15      atom numbers of atoms 1-5 connected to each atom */




/*     initialize the ring size and number of atoms to test */

    *iring = 0;
    nset = 0;
    if (*ia > 0) {
	nset = 1;
    }
    if (*ib > 0) {
	nset = 2;
    }
    if (*ic > 0) {
	nset = 3;
    }
    if (*id > 0) {
	nset = 4;
    }

/*     cannot be in a ring if the terminal atoms are univalent */

    if (nset == 1) {
	if (couple_1.n12[*ia - 1] <= 1) {
	    nset = 0;
	}
    } else if (nset == 2) {
/* Computing MIN */
	i__1 = couple_1.n12[*ia - 1], i__2 = couple_1.n12[*ib - 1];
	if (min(i__1,i__2) <= 1) {
	    nset = 0;
	}
    } else if (nset == 3) {
/* Computing MIN */
	i__1 = couple_1.n12[*ia - 1], i__2 = couple_1.n12[*ic - 1];
	if (min(i__1,i__2) <= 1) {
	    nset = 0;
	}
    } else if (nset == 4) {
/* Computing MIN */
	i__1 = couple_1.n12[*ia - 1], i__2 = couple_1.n12[*id - 1];
	if (min(i__1,i__2) <= 1) {
	    nset = 0;
	}
    }

/*     check the input atoms for sequential connectivity */

    if (nset > 1) {
	i__1 = couple_1.n12[*ia - 1];
	for (j = 1; j <= i__1; ++j) {
	    i__ = i12_ref(j, *ia);
	    if (*ib == i__) {
		if (nset == 2) {
		    goto L10;
		}
		i__2 = couple_1.n12[*ib - 1];
		for (k = 1; k <= i__2; ++k) {
		    m = i12_ref(k, *ib);
		    if (*ic == m) {
			if (nset == 3) {
			    goto L10;
			}
			i__3 = couple_1.n12[*ic - 1];
			for (p = 1; p <= i__3; ++p) {
			    q = i12_ref(p, *ic);
			    if (*id == q) {
				goto L10;
			    }
			}
		    }
		}
	    }
	}
	nset = 0;
L10:
	;
    }

/*     check for an atom contained inside a small ring */

    if (nset == 1) {
	i__1 = couple_1.n12[*ia - 1] - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__ = i12_ref(j, *ia);
	    i__2 = couple_1.n12[*ia - 1];
	    for (k = j + 1; k <= i__2; ++k) {
		m = i12_ref(k, *ia);
		i__3 = couple_1.n12[i__ - 1];
		for (p = 1; p <= i__3; ++p) {
		    if (m == i12_ref(p, i__)) {
			*iring = 3;
			goto L20;
		    }
		}
	    }
	}
	i__1 = couple_1.n12[*ia - 1] - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__ = i12_ref(j, *ia);
	    i__2 = couple_1.n12[*ia - 1];
	    for (k = j + 1; k <= i__2; ++k) {
		m = i12_ref(k, *ia);
		i__3 = couple_1.n12[i__ - 1];
		for (p = 1; p <= i__3; ++p) {
		    r__ = i12_ref(p, i__);
		    if (r__ != *ia) {
			i__4 = couple_1.n12[m - 1];
			for (q = 1; q <= i__4; ++q) {
			    if (r__ == i12_ref(q, m)) {
				*iring = 4;
				goto L20;
			    }
			}
		    }
		}
	    }
	}
	i__1 = couple_1.n13[*ia - 1] - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__ = i13_ref(j, *ia);
	    i__2 = couple_1.n13[*ia - 1];
	    for (k = j + 1; k <= i__2; ++k) {
		m = i13_ref(k, *ia);
		i__3 = couple_1.n12[i__ - 1];
		for (p = 1; p <= i__3; ++p) {
		    if (m == i12_ref(p, i__)) {
			*iring = 5;
			goto L20;
		    }
		}
		i__3 = couple_1.n13[i__ - 1];
		for (p = 1; p <= i__3; ++p) {
		    if (m == i13_ref(p, i__)) {
			*iring = 6;
			goto L20;
		    }
		}
	    }
	}
L20:

/*     check for a bond contained inside a small ring */

	;
    } else if (nset == 2) {
	i__1 = couple_1.n12[*ia - 1];
	for (j = 1; j <= i__1; ++j) {
	    i__ = i12_ref(j, *ia);
	    i__2 = couple_1.n12[*ib - 1];
	    for (k = 1; k <= i__2; ++k) {
		if (i__ == i12_ref(k, *ib)) {
		    *iring = 3;
		    goto L30;
		}
	    }
	}
	i__1 = couple_1.n12[*ia - 1];
	for (j = 1; j <= i__1; ++j) {
	    i__ = i12_ref(j, *ia);
	    if (*ib != i__) {
		i__2 = couple_1.n12[*ib - 1];
		for (k = 1; k <= i__2; ++k) {
		    m = i12_ref(k, *ib);
		    if (*ia != m) {
			i__3 = couple_1.n12[i__ - 1];
			for (p = 1; p <= i__3; ++p) {
			    if (m == i12_ref(p, i__)) {
				*iring = 4;
				goto L30;
			    }
			}
		    }
		}
	    }
	}
	i__1 = couple_1.n13[*ia - 1];
	for (j = 1; j <= i__1; ++j) {
	    i__ = i13_ref(j, *ia);
	    i__2 = couple_1.n13[*ib - 1];
	    for (k = 1; k <= i__2; ++k) {
		if (i__ == i13_ref(k, *ib)) {
		    *iring = 5;
		    goto L30;
		}
	    }
	}
	i__1 = couple_1.n12[*ia - 1];
	for (j = 1; j <= i__1; ++j) {
	    i__ = i12_ref(j, *ia);
	    if (*ib != i__) {
		i__2 = couple_1.n13[*ib - 1];
		for (k = 1; k <= i__2; ++k) {
		    m = i13_ref(k, *ib);
		    i__3 = couple_1.n13[i__ - 1];
		    for (p = 1; p <= i__3; ++p) {
			if (m == i13_ref(p, i__)) {
			    *iring = 6;
			    i__4 = couple_1.n12[*ia - 1];
			    for (q = 1; q <= i__4; ++q) {
				if (m == i12_ref(q, *ia)) {
				    *iring = 0;
				}
			    }
			    if (*iring == 6) {
				goto L30;
			    }
			}
		    }
		}
	    }
	}
L30:

/*     check for an angle contained inside a small ring */

	;
    } else if (nset == 3) {
	i__1 = couple_1.n12[*ia - 1];
	for (j = 1; j <= i__1; ++j) {
	    if (*ic == i12_ref(j, *ia)) {
		*iring = 3;
		goto L40;
	    }
	}
	i__1 = couple_1.n12[*ia - 1];
	for (j = 1; j <= i__1; ++j) {
	    i__ = i12_ref(j, *ia);
	    if (*ib != i__) {
		i__2 = couple_1.n12[*ic - 1];
		for (k = 1; k <= i__2; ++k) {
		    if (i__ == i12_ref(k, *ic)) {
			*iring = 4;
			goto L40;
		    }
		}
	    }
	}
	i__1 = couple_1.n12[*ia - 1];
	for (j = 1; j <= i__1; ++j) {
	    i__ = i12_ref(j, *ia);
	    if (*ib != i__) {
		i__2 = couple_1.n13[*ic - 1];
		for (k = 1; k <= i__2; ++k) {
		    if (i__ == i13_ref(k, *ic)) {
			*iring = 5;
			goto L40;
		    }
		}
	    }
	}
	i__1 = couple_1.n13[*ia - 1];
	for (j = 1; j <= i__1; ++j) {
	    i__ = i13_ref(j, *ia);
	    if (*ic != i__) {
		i__2 = couple_1.n13[*ic - 1];
		for (k = 1; k <= i__2; ++k) {
		    if (i__ == i13_ref(k, *ic)) {
			*iring = 6;
			goto L40;
		    }
		}
	    }
	}
L40:

/*     check for a torsion contained inside a small ring */

	;
    } else if (nset == 4) {
	i__1 = couple_1.n12[*ia - 1];
	for (j = 1; j <= i__1; ++j) {
	    if (*id == i12_ref(j, *ia)) {
		*iring = 4;
		goto L50;
	    }
	}
	i__1 = couple_1.n12[*ia - 1];
	for (j = 1; j <= i__1; ++j) {
	    i__ = i12_ref(j, *ia);
	    if (*ib != i__) {
		i__2 = couple_1.n12[*id - 1];
		for (k = 1; k <= i__2; ++k) {
		    if (i__ == i12_ref(k, *id)) {
			*iring = 5;
			goto L50;
		    }
		}
	    }
	}
	i__1 = couple_1.n12[*ia - 1];
	for (j = 1; j <= i__1; ++j) {
	    i__ = i12_ref(j, *ia);
	    if (*ib != i__) {
		i__2 = couple_1.n13[*id - 1];
		for (k = 1; k <= i__2; ++k) {
		    if (i__ == i13_ref(k, *id)) {
			*iring = 6;
			goto L50;
		    }
		}
	    }
	}
L50:
	;
    }
    return 0;
} /* chkring_ */

#undef i13_ref
#undef i12_ref


