/* attach.f -- translated by f2c (version 20050501).
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

/* Table of constant values */

static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine attach  --  setup of connectivity arrays  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "attach" generates lists of 1-3, 1-4 and 1-5 connectivities */
/*     starting from the previously determined list of attached */
/*     atoms (ie, 1-2 connectivity) */


/* Subroutine */ int attach_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 ATTACH  --  Too many 1-3 Connected Atom"
	    "s\002,\002 Attached to Atom\002,i6)";
    static char fmt_40[] = "(/,\002 ATTACH  --  Too many 1-4 Connected Atom"
	    "s\002,\002 Attached to Atom\002,i6)";
    static char fmt_60[] = "(/,\002 ATTACH  --  Too many 1-5 Connected Atom"
	    "s\002,\002 Attached to Atom\002,i6)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k, m, jj, kk;
    extern /* Subroutine */ int sort_(integer *, integer *), fatal_(void);

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_60, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]



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




/*     loop over all atoms finding all the 1-3 relationships; */
/*     note "n12" and "i12" have already been setup elsewhere */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	couple_1.n13[i__ - 1] = 0;
	i__2 = couple_1.n12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    jj = i12_ref(j, i__);
	    i__3 = couple_1.n12[jj - 1];
	    for (k = 1; k <= i__3; ++k) {
		kk = i12_ref(k, jj);
		if (kk == i__) {
		    goto L10;
		}
		i__4 = couple_1.n12[i__ - 1];
		for (m = 1; m <= i__4; ++m) {
		    if (kk == i12_ref(m, i__)) {
			goto L10;
		    }
		}
		++couple_1.n13[i__ - 1];
		i13_ref(couple_1.n13[i__ - 1], i__) = kk;
L10:
		;
	    }
	}
	if (couple_1.n13[i__ - 1] > 24) {
	    io___7.ciunit = iounit_1.iout;
	    s_wsfe(&io___7);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    e_wsfe();
	    fatal_();
	}
	sort_(&couple_1.n13[i__ - 1], &i13_ref(1, i__));
    }

/*     loop over all atoms finding all the 1-4 relationships */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	couple_1.n14[i__ - 1] = 0;
	i__2 = couple_1.n13[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    jj = i13_ref(j, i__);
	    i__3 = couple_1.n12[jj - 1];
	    for (k = 1; k <= i__3; ++k) {
		kk = i12_ref(k, jj);
		if (kk == i__) {
		    goto L30;
		}
		i__4 = couple_1.n12[i__ - 1];
		for (m = 1; m <= i__4; ++m) {
		    if (kk == i12_ref(m, i__)) {
			goto L30;
		    }
		}
		i__4 = couple_1.n13[i__ - 1];
		for (m = 1; m <= i__4; ++m) {
		    if (kk == i13_ref(m, i__)) {
			goto L30;
		    }
		}
		++couple_1.n14[i__ - 1];
		i14_ref(couple_1.n14[i__ - 1], i__) = kk;
L30:
		;
	    }
	}
	if (couple_1.n14[i__ - 1] > 72) {
	    io___8.ciunit = iounit_1.iout;
	    s_wsfe(&io___8);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    e_wsfe();
	    fatal_();
	}
	sort_(&couple_1.n14[i__ - 1], &i14_ref(1, i__));
    }

/*     loop over all atoms finding all the 1-5 relationships */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	couple_1.n15[i__ - 1] = 0;
	i__2 = couple_1.n14[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    jj = i14_ref(j, i__);
	    i__3 = couple_1.n12[jj - 1];
	    for (k = 1; k <= i__3; ++k) {
		kk = i12_ref(k, jj);
		if (kk == i__) {
		    goto L50;
		}
		i__4 = couple_1.n12[i__ - 1];
		for (m = 1; m <= i__4; ++m) {
		    if (kk == i12_ref(m, i__)) {
			goto L50;
		    }
		}
		i__4 = couple_1.n13[i__ - 1];
		for (m = 1; m <= i__4; ++m) {
		    if (kk == i13_ref(m, i__)) {
			goto L50;
		    }
		}
		i__4 = couple_1.n14[i__ - 1];
		for (m = 1; m <= i__4; ++m) {
		    if (kk == i14_ref(m, i__)) {
			goto L50;
		    }
		}
		++couple_1.n15[i__ - 1];
		i15_ref(couple_1.n15[i__ - 1], i__) = kk;
L50:
		;
	    }
	}
	if (couple_1.n15[i__ - 1] > 216) {
	    io___9.ciunit = iounit_1.iout;
	    s_wsfe(&io___9);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    e_wsfe();
	    fatal_();
	}
	sort_(&couple_1.n15[i__ - 1], &i15_ref(1, i__));
    }
    return 0;
} /* attach_ */

#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref


