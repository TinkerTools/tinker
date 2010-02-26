/* mutate.f -- translated by f2c (version 20050501).
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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal weight[5000];
    integer atmcls[5000], atmnum[5000], ligand[5000];
    char symbol[15000], describe[120000];
} katoms_;

#define katoms_1 katoms_

struct {
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

struct {
    doublereal lambda, vlambda, clambda, dlambda, mlambda, plambda;
    integer nmut, imut[25000], type0[25000], type1[25000], class0[25000], 
	    class1[25000];
    logical mut[25000];
} mutant_;

#define mutant_1 mutant_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__3 = 3;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine mutate  --  set the chimeric atoms and values  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "mutate" sets the list of chimeric atoms and mutational */
/*     parameters that are used during a free energy calculation */


/* Subroutine */ int mutate_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);

    /* Local variables */
    static integer i__, j, k, ia, it0, it1, list[50], next;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120], keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___7 = { 1, string, 0, 0, 120, 1 };
    static icilist io___8 = { 1, string, 0, 0, 120, 1 };
    static icilist io___9 = { 1, string, 0, 0, 120, 1 };
    static icilist io___10 = { 1, string, 0, 0, 120, 1 };
    static icilist io___11 = { 1, string, 0, 0, 120, 1 };
    static icilist io___12 = { 1, string, 0, 0, 120, 1 };
    static icilist io___13 = { 1, string, 0, 0, 120, 1 };
    static icilist io___18 = { 1, string, 1, 0, 120, 1 };



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
/*     ##  katoms.i  --  forcefield parameters for the atom types  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     weight     average atomic mass of each atom type */
/*     atmcls     atom class number for each of the atom types */
/*     atmnum     atomic number for each of the atom types */
/*     ligand     number of atoms to be attached to each atom type */
/*     symbol     modified atomic symbol for each atom type */
/*     describe   string identifying each of the atom types */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




/*     zero out the number of mutated atoms and atom list */

    mutant_1.nmut = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mutant_1.mut[i__ - 1] = FALSE_;
    }

/*     search keywords for free energy perturbation options */

    i__1 = keys_1.nkey;
    for (j = 1; j <= i__1; ++j) {
	next = 1;
	s_copy(record, keyline_ref(0, j), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "LAMBDA ", (ftnlen)7, (ftnlen)7) == 0) {
	    mutant_1.lambda = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___7);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&mutant_1.lambda, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	    mutant_1.vlambda = mutant_1.lambda;
	    mutant_1.clambda = mutant_1.lambda;
	    mutant_1.dlambda = mutant_1.lambda;
	    mutant_1.mlambda = mutant_1.lambda;
	    mutant_1.plambda = mutant_1.lambda;
	} else if (s_cmp(keyword, "VDW-LAMBDA ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___8);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&mutant_1.vlambda, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "CHG-LAMBDA ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___9);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&mutant_1.clambda, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "DPL-LAMBDA ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___10);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&mutant_1.dlambda, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "MPOLE-LAMBDA ", (ftnlen)13, (ftnlen)13) == 
		0) {
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___11);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&mutant_1.mlambda, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "POLAR-LAMBDA ", (ftnlen)13, (ftnlen)13) == 
		0) {
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___12);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&mutant_1.plambda, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "MUTATE ", (ftnlen)7, (ftnlen)7) == 0) {
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___13);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&it0, (ftnlen)sizeof(integer))
		    ;
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&it1, (ftnlen)sizeof(integer))
		    ;
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	    ++mutant_1.nmut;
	    mutant_1.imut[mutant_1.nmut - 1] = ia;
	    mutant_1.type0[mutant_1.nmut - 1] = it0;
	    mutant_1.type1[mutant_1.nmut - 1] = it1;
	    mutant_1.class0[mutant_1.nmut - 1] = katoms_1.atmcls[it0 - 1];
	    mutant_1.class1[mutant_1.nmut - 1] = katoms_1.atmcls[it1 - 1];
	    mutant_1.mut[ia - 1] = TRUE_;
	} else if (s_cmp(keyword, "LIGAND ", (ftnlen)7, (ftnlen)7) == 0) {
	    for (i__ = 1; i__ <= 50; ++i__) {
		list[i__ - 1] = 0;
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___18);
	    if (i__2 != 0) {
		goto L10;
	    }
	    for (i__ = 1; i__ <= 50; ++i__) {
		i__2 = do_lio(&c__3, &c__1, (char *)&list[i__ - 1], (ftnlen)
			sizeof(integer));
		if (i__2 != 0) {
		    goto L10;
		}
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
L10:
	    k = 1;
	    while(list[k - 1] != 0) {
		ia = list[k - 1];
		if (list[k - 1] > 0) {
		    if (! mutant_1.mut[ia - 1]) {
			++mutant_1.nmut;
			mutant_1.imut[mutant_1.nmut - 1] = ia;
			mutant_1.type0[mutant_1.nmut - 1] = 0;
			mutant_1.type1[mutant_1.nmut - 1] = atoms_1.type__[ia 
				- 1];
			mutant_1.class0[mutant_1.nmut - 1] = 0;
			mutant_1.class1[mutant_1.nmut - 1] = atmtyp_1.class__[
				ia - 1];
			mutant_1.mut[ia - 1] = TRUE_;
		    }
		    ++k;
		} else {
		    i__4 = (i__3 = list[k], abs(i__3));
		    for (i__ = (i__2 = list[k - 1], abs(i__2)); i__ <= i__4; 
			    ++i__) {
			if (! mutant_1.mut[i__ - 1]) {
			    ++mutant_1.nmut;
			    mutant_1.imut[mutant_1.nmut - 1] = i__;
			    mutant_1.type0[mutant_1.nmut - 1] = 0;
			    mutant_1.type1[mutant_1.nmut - 1] = 
				    atoms_1.type__[i__ - 1];
			    mutant_1.class0[mutant_1.nmut - 1] = 0;
			    mutant_1.class1[mutant_1.nmut - 1] = 
				    atmtyp_1.class__[i__ - 1];
			    mutant_1.mut[i__ - 1] = TRUE_;
			}
		    }
		    k += 2;
		}
	    }
	}
L20:
	;
    }
    return 0;
} /* mutate_ */

#undef keyline_ref


