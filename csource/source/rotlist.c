/* rotlist.f -- translated by f2c (version 20050501).
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
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

struct {
    integer nrot, rot[25000];
    logical use_short__;
} rotate_;

#define rotate_1 rotate_

struct {
    integer nadd, iadd[50000]	/* was [2][25000] */, ndel, idel[50000]	/* 
	    was [2][25000] */;
} zclose_;

#define zclose_1 zclose_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine rotlist  --  find atoms on one side of a bond  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "rotlist" generates the minimum list of all the atoms lying */
/*     to one side of a pair of directly bonded atoms; optionally */
/*     finds the minimal list by choosing the side with fewer atoms */


/* Subroutine */ int rotlist_(integer *base, integer *partner)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 ROTLIST  --  Maximum Valence Exceeded"
	    ";\002,\002 Increase MAXVAL\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, k, ia, ib, mark, swap, list[25001], test;
    extern /* Subroutine */ int fatal_(void);
    static logical bonded;
    static integer nattach;

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 0, 0, fmt_10, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define iadd_ref(a_1,a_2) zclose_1.iadd[(a_2)*2 + a_1 - 3]
#define idel_ref(a_1,a_2) zclose_1.idel[(a_2)*2 + a_1 - 3]



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
/*     ##  molcul.i  --  individual molecules within current system  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     molmass   molecular weight for each molecule in the system */
/*     totmass   total weight of all the molecules in the system */
/*     nmol      total number of separate molecules in the system */
/*     kmol      contiguous list of the atoms in each molecule */
/*     imol      first and last atom of each molecule in the list */
/*     molcule   number of the molecule to which each atom belongs */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  rotate.i  --  molecule partitions for rotation of a bond  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     nrot        total number of atoms moving when bond rotates */
/*     rot         atom numbers of atoms moving when bond rotates */
/*     use_short   logical flag governing use of shortest atom list */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  zclose.i  --  ring openings and closures for Z-matrix  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     nadd   number of added bonds between Z-matrix atoms */
/*     iadd   numbers of the atom pairs defining added bonds */
/*     ndel   number of bonds between Z-matrix bonds to delete */
/*     idel   numbers of the atom pairs defining deleted bonds */




/*     initialize the number of atoms to one side of the bond */

    rotate_1.nrot = 0;

/*     remove any bonds needed for intramolecular ring closures */

    i__1 = zclose_1.nadd;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iadd_ref(1, i__);
	ib = iadd_ref(2, i__);
	if (molcul_1.molcule[ia - 1] == molcul_1.molcule[ib - 1]) {
	    i__2 = couple_1.n12[ia - 1];
	    for (k = 1; k <= i__2; ++k) {
		if (i12_ref(k, ia) == ib) {
		    i12_ref(k, ia) = 0;
		}
	    }
	    i__2 = couple_1.n12[ib - 1];
	    for (k = 1; k <= i__2; ++k) {
		if (i12_ref(k, ib) == ia) {
		    i12_ref(k, ib) = 0;
		}
	    }
	}
    }

/*     add any links needed to make intermolecular connections */

    i__1 = zclose_1.ndel;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = idel_ref(1, i__);
	ib = idel_ref(2, i__);
	if (molcul_1.molcule[ia - 1] != molcul_1.molcule[ib - 1]) {
	    if (couple_1.n12[ia - 1] == 8 || couple_1.n12[ib - 1] == 8) {
		io___5.ciunit = iounit_1.iout;
		s_wsfe(&io___5);
		e_wsfe();
		fatal_();
	    }
	    ++couple_1.n12[ia - 1];
	    i12_ref(couple_1.n12[ia - 1], ia) = ib;
	    ++couple_1.n12[ib - 1];
	    i12_ref(couple_1.n12[ib - 1], ib) = ia;
	}
    }

/*     check to see if the two atoms are still directly bonded */

    bonded = FALSE_;
    i__1 = couple_1.n12[*base - 1];
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i12_ref(i__, *base) == *partner) {
	    bonded = TRUE_;
	}
    }

/*     make a list of atoms to one side of this pair of atoms, */
/*     taking note of any rings in which the atom pair resides */

    if (bonded) {
	list[0] = 1;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    rotate_1.rot[i__ - 1] = 0;
	}
L20:
	rotate_1.nrot = 0;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    list[i__] = 0;
	}
	list[*base] = 1;
	list[*partner] = 1;
	nattach = couple_1.n12[*base - 1];
	i__1 = nattach;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    test = i12_ref(i__, *base);
	    if (list[test] == 0) {
		++rotate_1.nrot;
		if (rotate_1.use_short__ && rotate_1.nrot >= atoms_1.n / 2) {
		    goto L30;
		}
		rotate_1.rot[rotate_1.nrot - 1] = test;
		list[test] = 1;
	    }
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    mark = rotate_1.rot[i__ - 1];
	    if (mark == 0) {
		goto L40;
	    }
	    nattach = couple_1.n12[mark - 1];
	    if (nattach > 1) {
		i__2 = nattach;
		for (k = 1; k <= i__2; ++k) {
		    test = i12_ref(k, mark);
		    if (list[test] == 0) {
			++rotate_1.nrot;
			if (rotate_1.use_short__ && rotate_1.nrot >= 
				atoms_1.n / 2) {
			    goto L30;
			}
			rotate_1.rot[rotate_1.nrot - 1] = test;
			list[test] = 1;
		    }
		}
	    }
	}

/*     the list contains over half the total number of atoms, */
/*     so reverse the base and partner atoms, then start over */

L30:
	swap = *base;
	*base = *partner;
	*partner = swap;
	i__1 = rotate_1.nrot;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    rotate_1.rot[i__ - 1] = 0;
	}
	goto L20;
    }

/*     remove links added to make intermolecular connections */

L40:
    i__1 = zclose_1.ndel;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = idel_ref(1, i__);
	ib = idel_ref(2, i__);
	if (molcul_1.molcule[ia - 1] != molcul_1.molcule[ib - 1]) {
	    --couple_1.n12[ia - 1];
	    --couple_1.n12[ib - 1];
	}
    }

/*     add any bonds required for intramolecular ring closures */

    i__1 = zclose_1.nadd;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iadd_ref(1, i__);
	ib = iadd_ref(2, i__);
	if (molcul_1.molcule[ia - 1] == molcul_1.molcule[ib - 1]) {
	    i__2 = couple_1.n12[ia - 1];
	    for (k = 1; k <= i__2; ++k) {
		if (i12_ref(k, ia) == 0) {
		    i12_ref(k, ia) = ib;
		    goto L50;
		}
	    }
L50:
	    i__2 = couple_1.n12[ib - 1];
	    for (k = 1; k <= i__2; ++k) {
		if (i12_ref(k, ib) == 0) {
		    i12_ref(k, ib) = ia;
		    goto L60;
		}
	    }
L60:
	    ;
	}
    }
    return 0;
} /* rotlist_ */

#undef idel_ref
#undef iadd_ref
#undef i12_ref


