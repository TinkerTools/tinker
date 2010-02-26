/* fft3d.f -- translated by f2c (version 20050501).
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
    doublereal bsmod1[100], bsmod2[100], bsmod3[100], table[1200]	/* 
	    was [400][3] */, qgrid[2000000]	/* was [2][100][100][100] */, 
	    qfac[1000000]	/* was [100][100][100] */, thetai1[1000000]	
	    /* was [4][10][25000] */, thetai2[1000000]	/* was [4][10][25000] 
	    */, thetai3[1000000]	/* was [4][10][25000] */;
    integer nfft1, nfft2, nfft3, bsorder, iprime[45]	/* was [15][3] */, 
	    igrid[75000]	/* was [3][25000] */;
} pme_;

#define pme_1 pme_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ####################################################### */
/*     ##                                                   ## */
/*     ##  subroutine fftsetup  --  3-D FFT initialization  ## */
/*     ##                                                   ## */
/*     ####################################################### */


/*     "fftsetup" does the initialization for a 3-D FFT via */
/*     three separate 1-D initializations */


/* Subroutine */ int fftsetup_(void)
{
    extern /* Subroutine */ int cffti_(integer *, doublereal *, integer *);


#define table_ref(a_1,a_2) pme_1.table[(a_2)*400 + a_1 - 401]
#define iprime_ref(a_1,a_2) pme_1.iprime[(a_2)*15 + a_1 - 16]



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




/*     perform initialization along X, Y and Z directions */



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  pme.i  --  values for particle mesh Ewald summation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxorder   maximum order of the B-spline approximation */
/*     maxprime   maximum number of prime factors of FFT dimension */
/*     maxtable   maximum size of the FFT table array */

/*     bsmod1     B-spline moduli along the a-axis direction */
/*     bsmod2     B-spline moduli along the b-axis direction */
/*     bsmod3     B-spline moduli along the c-axis direction */
/*     table      intermediate array used by the FFT calculation */
/*     qgrid      values on the particle mesh Ewald charge grid */
/*     qfac       prefactors for particle mesh Ewald charge grid */
/*     thetai1    B-spline coefficients along the a-axis */
/*     thetai2    B-spline coefficients along the b-axis */
/*     thetai3    B-spline coefficients along the c-axis */
/*     nfft1      number of grid points along the a-axis direction */
/*     nfft2      number of grid points along the b-axis direction */
/*     nfft3      number of grid points along the c-axis direction */
/*     bsorder    order of the PME B-spline approximation */
/*     iprime     prime factorization of each FFT dimension */
/*     igrid      initial Ewald charge grid values for B-spline */


    cffti_(&pme_1.nfft1, &table_ref(1, 1), &iprime_ref(1, 1));
    cffti_(&pme_1.nfft2, &table_ref(1, 2), &iprime_ref(1, 2));
    cffti_(&pme_1.nfft3, &table_ref(1, 3), &iprime_ref(1, 3));
    return 0;
} /* fftsetup_ */

#undef iprime_ref
#undef table_ref




/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  subroutine fftfront  --  3-D FFT forward transform  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     "fftfront" does a 3-D FFT forward transform via three */
/*     separate 1-D transformations */


/* Subroutine */ int fftfront_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal work[200]	/* was [2][100] */;
    extern /* Subroutine */ int cfftf_(integer *, doublereal *, doublereal *, 
	    integer *);


#define work_ref(a_1,a_2) work[(a_2)*2 + a_1 - 3]
#define table_ref(a_1,a_2) pme_1.table[(a_2)*400 + a_1 - 401]
#define qgrid_ref(a_1,a_2,a_3,a_4) pme_1.qgrid[(((a_4)*100 + (a_3))*100 + \
(a_2))*2 + a_1 - 20203]
#define iprime_ref(a_1,a_2) pme_1.iprime[(a_2)*15 + a_1 - 16]



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
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  pme.i  --  values for particle mesh Ewald summation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxorder   maximum order of the B-spline approximation */
/*     maxprime   maximum number of prime factors of FFT dimension */
/*     maxtable   maximum size of the FFT table array */

/*     bsmod1     B-spline moduli along the a-axis direction */
/*     bsmod2     B-spline moduli along the b-axis direction */
/*     bsmod3     B-spline moduli along the c-axis direction */
/*     table      intermediate array used by the FFT calculation */
/*     qgrid      values on the particle mesh Ewald charge grid */
/*     qfac       prefactors for particle mesh Ewald charge grid */
/*     thetai1    B-spline coefficients along the a-axis */
/*     thetai2    B-spline coefficients along the b-axis */
/*     thetai3    B-spline coefficients along the c-axis */
/*     nfft1      number of grid points along the a-axis direction */
/*     nfft2      number of grid points along the b-axis direction */
/*     nfft3      number of grid points along the c-axis direction */
/*     bsorder    order of the PME B-spline approximation */
/*     iprime     prime factorization of each FFT dimension */
/*     igrid      initial Ewald charge grid values for B-spline */




/*     perform forward transform along X, Y and Z directions */

    i__1 = pme_1.nfft3;
    for (k = 1; k <= i__1; ++k) {
	i__2 = pme_1.nfft2;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = pme_1.nfft1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		work_ref(1, i__) = qgrid_ref(1, i__, j, k);
		work_ref(2, i__) = qgrid_ref(2, i__, j, k);
	    }
	    cfftf_(&pme_1.nfft1, work, &table_ref(1, 1), &iprime_ref(1, 1));
	    i__3 = pme_1.nfft1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		qgrid_ref(1, i__, j, k) = work_ref(1, i__);
		qgrid_ref(2, i__, j, k) = work_ref(2, i__);
	    }
	}
    }
    i__1 = pme_1.nfft3;
    for (k = 1; k <= i__1; ++k) {
	i__2 = pme_1.nfft1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = pme_1.nfft2;
	    for (j = 1; j <= i__3; ++j) {
		work_ref(1, j) = qgrid_ref(1, i__, j, k);
		work_ref(2, j) = qgrid_ref(2, i__, j, k);
	    }
	    cfftf_(&pme_1.nfft2, work, &table_ref(1, 2), &iprime_ref(1, 2));
	    i__3 = pme_1.nfft2;
	    for (j = 1; j <= i__3; ++j) {
		qgrid_ref(1, i__, j, k) = work_ref(1, j);
		qgrid_ref(2, i__, j, k) = work_ref(2, j);
	    }
	}
    }
    i__1 = pme_1.nfft1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = pme_1.nfft2;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = pme_1.nfft3;
	    for (k = 1; k <= i__3; ++k) {
		work_ref(1, k) = qgrid_ref(1, i__, j, k);
		work_ref(2, k) = qgrid_ref(2, i__, j, k);
	    }
	    cfftf_(&pme_1.nfft3, work, &table_ref(1, 3), &iprime_ref(1, 3));
	    i__3 = pme_1.nfft3;
	    for (k = 1; k <= i__3; ++k) {
		qgrid_ref(1, i__, j, k) = work_ref(1, k);
		qgrid_ref(2, i__, j, k) = work_ref(2, k);
	    }
	}
    }
    return 0;
} /* fftfront_ */

#undef iprime_ref
#undef qgrid_ref
#undef table_ref
#undef work_ref




/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  subroutine fftback  --  3-D FFT backward transform  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     "fftback" does a 3-D FFT backward transform via three */
/*     separate 1-D transformations */


/* Subroutine */ int fftback_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal work[200]	/* was [2][100] */;
    extern /* Subroutine */ int cfftb_(integer *, doublereal *, doublereal *, 
	    integer *);


#define work_ref(a_1,a_2) work[(a_2)*2 + a_1 - 3]
#define table_ref(a_1,a_2) pme_1.table[(a_2)*400 + a_1 - 401]
#define qgrid_ref(a_1,a_2,a_3,a_4) pme_1.qgrid[(((a_4)*100 + (a_3))*100 + \
(a_2))*2 + a_1 - 20203]
#define iprime_ref(a_1,a_2) pme_1.iprime[(a_2)*15 + a_1 - 16]



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
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  pme.i  --  values for particle mesh Ewald summation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxorder   maximum order of the B-spline approximation */
/*     maxprime   maximum number of prime factors of FFT dimension */
/*     maxtable   maximum size of the FFT table array */

/*     bsmod1     B-spline moduli along the a-axis direction */
/*     bsmod2     B-spline moduli along the b-axis direction */
/*     bsmod3     B-spline moduli along the c-axis direction */
/*     table      intermediate array used by the FFT calculation */
/*     qgrid      values on the particle mesh Ewald charge grid */
/*     qfac       prefactors for particle mesh Ewald charge grid */
/*     thetai1    B-spline coefficients along the a-axis */
/*     thetai2    B-spline coefficients along the b-axis */
/*     thetai3    B-spline coefficients along the c-axis */
/*     nfft1      number of grid points along the a-axis direction */
/*     nfft2      number of grid points along the b-axis direction */
/*     nfft3      number of grid points along the c-axis direction */
/*     bsorder    order of the PME B-spline approximation */
/*     iprime     prime factorization of each FFT dimension */
/*     igrid      initial Ewald charge grid values for B-spline */




/*     perform backward transform along X, Y and Z directions */

    i__1 = pme_1.nfft3;
    for (k = 1; k <= i__1; ++k) {
	i__2 = pme_1.nfft2;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = pme_1.nfft1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		work_ref(1, i__) = qgrid_ref(1, i__, j, k);
		work_ref(2, i__) = qgrid_ref(2, i__, j, k);
	    }
	    cfftb_(&pme_1.nfft1, work, &table_ref(1, 1), &iprime_ref(1, 1));
	    i__3 = pme_1.nfft1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		qgrid_ref(1, i__, j, k) = work_ref(1, i__);
		qgrid_ref(2, i__, j, k) = work_ref(2, i__);
	    }
	}
    }
    i__1 = pme_1.nfft3;
    for (k = 1; k <= i__1; ++k) {
	i__2 = pme_1.nfft1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = pme_1.nfft2;
	    for (j = 1; j <= i__3; ++j) {
		work_ref(1, j) = qgrid_ref(1, i__, j, k);
		work_ref(2, j) = qgrid_ref(2, i__, j, k);
	    }
	    cfftb_(&pme_1.nfft2, work, &table_ref(1, 2), &iprime_ref(1, 2));
	    i__3 = pme_1.nfft2;
	    for (j = 1; j <= i__3; ++j) {
		qgrid_ref(1, i__, j, k) = work_ref(1, j);
		qgrid_ref(2, i__, j, k) = work_ref(2, j);
	    }
	}
    }
    i__1 = pme_1.nfft1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = pme_1.nfft2;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = pme_1.nfft3;
	    for (k = 1; k <= i__3; ++k) {
		work_ref(1, k) = qgrid_ref(1, i__, j, k);
		work_ref(2, k) = qgrid_ref(2, i__, j, k);
	    }
	    cfftb_(&pme_1.nfft3, work, &table_ref(1, 3), &iprime_ref(1, 3));
	    i__3 = pme_1.nfft3;
	    for (k = 1; k <= i__3; ++k) {
		qgrid_ref(1, i__, j, k) = work_ref(1, k);
		qgrid_ref(2, i__, j, k) = work_ref(2, k);
	    }
	}
    }
    return 0;
} /* fftback_ */

#undef iprime_ref
#undef qgrid_ref
#undef table_ref
#undef work_ref


