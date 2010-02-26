/* vibrot.f -- translated by f2c (version 20050501).
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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal dihed[1000];
    integer nomega, iomega[2000]	/* was [2][1000] */, zline[1000];
} omega_;

#define omega_1 omega_

/* Table of constant values */

static integer c__1 = 1;
static integer c__1000 = 1000;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  program vibrot  --  vibrational analysis over torsions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "vibrot" computes the eigenvalues and eigenvectors of the */
/*     torsional Hessian matrix */

/*     literature reference: */

/*     M. Levitt, C. Sander and P. S. Stern, "Protein Normal-mode */
/*     Dynamics: Trypsin Inhibitor, Crambin, Ribonuclease and Lysozyme", */
/*     Journal of Molecular Biology, 181, 423-447 (1985) */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Diagonal of the Torsional Hessian :\002,"
	    "/)";
    static char fmt_20[] = "(4(i8,f11.3))";
    static char fmt_30[] = "(/,\002 Torsional Hessian Matrix Elements :\002)";
    static char fmt_40[] = "()";
    static char fmt_50[] = "(6f13.4)";
    static char fmt_60[] = "(/,\002 Eigenvalues of the Hessian Matrix :\002,"
	    "/)";
    static char fmt_70[] = "(4(i8,f11.3))";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    static doublereal a[1000], b[1000];
    static integer i__, j;
    static doublereal p[1000], w[1000], ta[1000], tb[1000], ty[1000], hrot[
	    1000000]	/* was [1000][1000] */;
    extern /* Subroutine */ int diagq_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal eigen[1000];
    extern /* Subroutine */ int final_(void);
    static integer ihess;
    static doublereal vects[1000000]	/* was [1000][1000] */;
    extern /* Subroutine */ int getint_(void);
    static doublereal matrix[500500];
    extern /* Subroutine */ int initial_(void), hessrot_(char *, doublereal *,
	     ftnlen), initrot_(void);

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___3 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_70, 0 };



#define hrot_ref(a_1,a_2) hrot[(a_2)*1000 + a_1 - 1001]



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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




/*     compute torsional Hessian matrix elements */

    initial_();
    getint_();
    mechanic_();
    initrot_();
    hessrot_("FULL", hrot, (ftnlen)4);

/*     write out the torsional Hessian diagonal */

    io___2.ciunit = iounit_1.iout;
    s_wsfe(&io___2);
    e_wsfe();
    io___3.ciunit = iounit_1.iout;
    s_wsfe(&io___3);
    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&hrot_ref(i__, i__), (ftnlen)sizeof(doublereal))
		;
    }
    e_wsfe();

/*     write out the torsional Hessian elements */

    if (omega_1.nomega <= 30) {
	io___5.ciunit = iounit_1.iout;
	s_wsfe(&io___5);
	e_wsfe();
	i__1 = omega_1.nomega;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___6.ciunit = iounit_1.iout;
	    s_wsfe(&io___6);
	    e_wsfe();
	    io___7.ciunit = iounit_1.iout;
	    s_wsfe(&io___7);
	    i__2 = omega_1.nomega;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&hrot_ref(j, i__), (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	}
    }

/*     place Hessian elements into triangular form */

    ihess = 0;
    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = omega_1.nomega;
	for (j = i__; j <= i__2; ++j) {
	    ++ihess;
	    matrix[ihess - 1] = hrot_ref(i__, j);
	}
    }

/*     perform diagonalization to get Hessian eigenvalues */

    diagq_(&omega_1.nomega, &c__1000, &omega_1.nomega, matrix, eigen, vects, 
	    a, b, p, w, ta, tb, ty);
    io___20.ciunit = iounit_1.iout;
    s_wsfe(&io___20);
    e_wsfe();
    io___21.ciunit = iounit_1.iout;
    s_wsfe(&io___21);
    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&eigen[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef hrot_ref


/* Main program alias */ int vibrot_ () { MAIN__ (); return 0; }
