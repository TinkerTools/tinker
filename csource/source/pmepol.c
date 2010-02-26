/* pmepol.f -- translated by f2c (version 20050501).
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
    doublereal xbox, ybox, zbox, alpha, beta, gamma, xbox2, ybox2, zbox2, 
	    box34, lvec[9]	/* was [3][3] */, recip[9]	/* was [3][3] 
	    */, volbox, beta_sin__, beta_cos__, gamma_sin__, gamma_cos__, 
	    beta_term__, gamma_term__;
    logical orthogonal, monoclinic, triclinic, octahedron;
    char spacegrp[10];
} boxes_;

#define boxes_1 boxes_

struct {
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

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



/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Thomas Darden & Jay William Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  routines below implement various coordinate and B-spline  ## */
/*     ##  manipulations for particle mesh Ewald summation applied   ## */
/*     ##  to polarizable atomic multipoles; modified from original  ## */
/*     ##  code due to Thomas Darden, NIEHS, Research Triangle, NC   ## */
/*     ##                                                            ## */
/*     ################################################################ */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine bspline_fill  --  set multipole B-spline coeffs  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "bspline_fill" finds B-spline coefficients and derivatives */
/*     for multipole sites along the fractional coordinate axes */


/* Subroutine */ int bspline_fill__(void)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double d_nint(doublereal *);

    /* Local variables */
    static integer i__;
    static doublereal w;
    static integer ii;
    static doublereal fr;
    static integer ifr;
    extern /* Subroutine */ int bspline_gen__(doublereal *, doublereal *);


#define igrid_ref(a_1,a_2) pme_1.igrid[(a_2)*3 + a_1 - 4]
#define recip_ref(a_1,a_2) boxes_1.recip[(a_2)*3 + a_1 - 4]
#define thetai1_ref(a_1,a_2,a_3) pme_1.thetai1[((a_3)*10 + (a_2))*4 + a_1 \
- 45]
#define thetai2_ref(a_1,a_2,a_3) pme_1.thetai2[((a_3)*10 + (a_2))*4 + a_1 \
- 45]
#define thetai3_ref(a_1,a_2,a_3) pme_1.thetai3[((a_3)*10 + (a_2))*4 + a_1 \
- 45]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




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




/*     get the B-spline coefficients for each multipole site */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	w = atoms_1.x[ii - 1] * recip_ref(1, 1) + atoms_1.y[ii - 1] * 
		recip_ref(2, 1) + atoms_1.z__[ii - 1] * recip_ref(3, 1);
	fr = (doublereal) pme_1.nfft1 * (w - d_nint(&w) + .5);
	ifr = (integer) fr;
	w = fr - (doublereal) ifr;
	igrid_ref(1, i__) = ifr - pme_1.bsorder;
	bspline_gen__(&w, &thetai1_ref(1, 1, i__));
	w = atoms_1.x[ii - 1] * recip_ref(1, 2) + atoms_1.y[ii - 1] * 
		recip_ref(2, 2) + atoms_1.z__[ii - 1] * recip_ref(3, 2);
	fr = (doublereal) pme_1.nfft2 * (w - d_nint(&w) + .5);
	ifr = (integer) fr;
	w = fr - (doublereal) ifr;
	igrid_ref(2, i__) = ifr - pme_1.bsorder;
	bspline_gen__(&w, &thetai2_ref(1, 1, i__));
	w = atoms_1.x[ii - 1] * recip_ref(1, 3) + atoms_1.y[ii - 1] * 
		recip_ref(2, 3) + atoms_1.z__[ii - 1] * recip_ref(3, 3);
	fr = (doublereal) pme_1.nfft3 * (w - d_nint(&w) + .5);
	ifr = (integer) fr;
	w = fr - (doublereal) ifr;
	igrid_ref(3, i__) = ifr - pme_1.bsorder;
	bspline_gen__(&w, &thetai3_ref(1, 1, i__));
    }
    return 0;
} /* bspline_fill__ */

#undef thetai3_ref
#undef thetai2_ref
#undef thetai1_ref
#undef recip_ref
#undef igrid_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine bspline_gen  --  single-site B-spline coeffs  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "bspline_gen" gets B-spline coefficients and derivatives for */
/*     a single multipole site along a particular direction */


/* Subroutine */ int bspline_gen__(doublereal *w, doublereal *thetai)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal denom, array[100]	/* was [10][10] */;


#define array_ref(a_1,a_2) array[(a_2)*10 + a_1 - 11]
#define thetai_ref(a_1,a_2) thetai[(a_2)*4 + a_1]



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




/*     initialization to get to 2nd order recursion */

    /* Parameter adjustments */
    thetai -= 5;

    /* Function Body */
    array_ref(2, 2) = *w;
    array_ref(2, 1) = 1. - *w;

/*     perform one pass to get to 3rd order recursion */

    array_ref(3, 3) = *w * .5 * array_ref(2, 2);
    array_ref(3, 2) = ((*w + 1.) * array_ref(2, 1) + (2. - *w) * array_ref(2, 
	    2)) * .5;
    array_ref(3, 1) = (1. - *w) * .5 * array_ref(2, 1);

/*     compute standard B-spline recursion to desired order */

    i__1 = pme_1.bsorder;
    for (i__ = 4; i__ <= i__1; ++i__) {
	k = i__ - 1;
	denom = 1. / (doublereal) k;
	array_ref(i__, i__) = denom * *w * array_ref(k, k);
	i__2 = i__ - 2;
	for (j = 1; j <= i__2; ++j) {
	    array_ref(i__, i__ - j) = denom * ((*w + (doublereal) j) * 
		    array_ref(k, i__ - j - 1) + ((doublereal) (i__ - j) - *w) 
		    * array_ref(k, i__ - j));
	}
	array_ref(i__, 1) = denom * (1. - *w) * array_ref(k, 1);
    }

/*     get coefficients for the B-spline first derivative */

    k = pme_1.bsorder - 1;
    array_ref(k, pme_1.bsorder) = array_ref(k, pme_1.bsorder - 1);
    for (i__ = pme_1.bsorder - 1; i__ >= 2; --i__) {
	array_ref(k, i__) = array_ref(k, i__ - 1) - array_ref(k, i__);
    }
    array_ref(k, 1) = -array_ref(k, 1);

/*     get coefficients for the B-spline second derivative */

    k = pme_1.bsorder - 2;
    array_ref(k, pme_1.bsorder - 1) = array_ref(k, pme_1.bsorder - 2);
    for (i__ = pme_1.bsorder - 2; i__ >= 2; --i__) {
	array_ref(k, i__) = array_ref(k, i__ - 1) - array_ref(k, i__);
    }
    array_ref(k, 1) = -array_ref(k, 1);
    array_ref(k, pme_1.bsorder) = array_ref(k, pme_1.bsorder - 1);
    for (i__ = pme_1.bsorder - 1; i__ >= 2; --i__) {
	array_ref(k, i__) = array_ref(k, i__ - 1) - array_ref(k, i__);
    }
    array_ref(k, 1) = -array_ref(k, 1);

/*     get coefficients for the B-spline third derivative */

    k = pme_1.bsorder - 3;
    array_ref(k, pme_1.bsorder - 2) = array_ref(k, pme_1.bsorder - 3);
    for (i__ = pme_1.bsorder - 3; i__ >= 2; --i__) {
	array_ref(k, i__) = array_ref(k, i__ - 1) - array_ref(k, i__);
    }
    array_ref(k, 1) = -array_ref(k, 1);
    array_ref(k, pme_1.bsorder - 1) = array_ref(k, pme_1.bsorder - 2);
    for (i__ = pme_1.bsorder - 2; i__ >= 2; --i__) {
	array_ref(k, i__) = array_ref(k, i__ - 1) - array_ref(k, i__);
    }
    array_ref(k, 1) = -array_ref(k, 1);
    array_ref(k, pme_1.bsorder) = array_ref(k, pme_1.bsorder - 1);
    for (i__ = pme_1.bsorder - 1; i__ >= 2; --i__) {
	array_ref(k, i__) = array_ref(k, i__ - 1) - array_ref(k, i__);
    }
    array_ref(k, 1) = -array_ref(k, 1);

/*     copy coefficients from temporary to permanent storage */

    i__1 = pme_1.bsorder;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    thetai_ref(j, i__) = array_ref(pme_1.bsorder - j + 1, i__);
	}
    }
    return 0;
} /* bspline_gen__ */

#undef thetai_ref
#undef array_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine grid_mpole  --  put multipoles on PME grid  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "grid_mpole" places the fractional atomic multipoles onto */
/*     the particle mesh Ewald grid */


/* Subroutine */ int grid_mpole__(doublereal *fmp)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer i_sign(integer *, integer *);

    /* Local variables */
    static integer i__, j, k, m, i0, j0, k0;
    static doublereal t0, u0, v0, v1, u1, t1, v2, u2, t2;
    static integer it1, it2, it3, igrd0, jgrd0, kgrd0;
    static doublereal term0, term1, term2;


#define fmp_ref(a_1,a_2) fmp[(a_2)*10 + a_1]
#define igrid_ref(a_1,a_2) pme_1.igrid[(a_2)*3 + a_1 - 4]
#define qgrid_ref(a_1,a_2,a_3,a_4) pme_1.qgrid[(((a_4)*100 + (a_3))*100 + \
(a_2))*2 + a_1 - 20203]
#define thetai1_ref(a_1,a_2,a_3) pme_1.thetai1[((a_3)*10 + (a_2))*4 + a_1 \
- 45]
#define thetai2_ref(a_1,a_2,a_3) pme_1.thetai2[((a_3)*10 + (a_2))*4 + a_1 \
- 45]
#define thetai3_ref(a_1,a_2,a_3) pme_1.thetai3[((a_3)*10 + (a_2))*4 + a_1 \
- 45]



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
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




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




/*     zero out the particle mesh Ewald charge grid */

    /* Parameter adjustments */
    fmp -= 11;

    /* Function Body */
    i__1 = pme_1.nfft3;
    for (k = 1; k <= i__1; ++k) {
	i__2 = pme_1.nfft2;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = pme_1.nfft1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		qgrid_ref(1, i__, j, k) = 0.;
		qgrid_ref(2, i__, j, k) = 0.;
	    }
	}
    }

/*     put the permanent multipole moments onto the grid */

    i__1 = mpole_1.npole;
    for (m = 1; m <= i__1; ++m) {
	igrd0 = igrid_ref(1, m);
	jgrd0 = igrid_ref(2, m);
	kgrd0 = igrid_ref(3, m);
	k0 = kgrd0;
	i__2 = pme_1.bsorder;
	for (it3 = 1; it3 <= i__2; ++it3) {
	    ++k0;
	    k = k0 + 1 + (pme_1.nfft3 - i_sign(&pme_1.nfft3, &k0)) / 2;
	    v0 = thetai3_ref(1, it3, m);
	    v1 = thetai3_ref(2, it3, m);
	    v2 = thetai3_ref(3, it3, m);
	    j0 = jgrd0;
	    i__3 = pme_1.bsorder;
	    for (it2 = 1; it2 <= i__3; ++it2) {
		++j0;
		j = j0 + 1 + (pme_1.nfft2 - i_sign(&pme_1.nfft2, &j0)) / 2;
		u0 = thetai2_ref(1, it2, m);
		u1 = thetai2_ref(2, it2, m);
		u2 = thetai2_ref(3, it2, m);
		term0 = fmp_ref(1, m) * u0 * v0 + fmp_ref(3, m) * u1 * v0 + 
			fmp_ref(4, m) * u0 * v1 + fmp_ref(6, m) * u2 * v0 + 
			fmp_ref(7, m) * u0 * v2 + fmp_ref(10, m) * u1 * v1;
		term1 = fmp_ref(2, m) * u0 * v0 + fmp_ref(8, m) * u1 * v0 + 
			fmp_ref(9, m) * u0 * v1;
		term2 = fmp_ref(5, m) * u0 * v0;
		i0 = igrd0;
		i__4 = pme_1.bsorder;
		for (it1 = 1; it1 <= i__4; ++it1) {
		    ++i0;
		    i__ = i0 + 1 + (pme_1.nfft1 - i_sign(&pme_1.nfft1, &i0)) /
			     2;
		    t0 = thetai1_ref(1, it1, m);
		    t1 = thetai1_ref(2, it1, m);
		    t2 = thetai1_ref(3, it1, m);
		    qgrid_ref(1, i__, j, k) = qgrid_ref(1, i__, j, k) + term0 
			    * t0 + term1 * t1 + term2 * t2;
		}
	    }
	}
    }
    return 0;
} /* grid_mpole__ */

#undef thetai3_ref
#undef thetai2_ref
#undef thetai1_ref
#undef qgrid_ref
#undef igrid_ref
#undef fmp_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine grid_uind  --  put induced dipoles on PME grid  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "grid_uind" places the fractional induced dipoles onto the */
/*     particle mesh Ewald grid */


/* Subroutine */ int grid_uind__(doublereal *fuind, doublereal *fuinp)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer i_sign(integer *, integer *);

    /* Local variables */
    static integer i__, j, k, m, i0, j0, k0;
    static doublereal t0, u0, v0, v1, u1, t1;
    static integer it1, it2, it3, igrd0, jgrd0, kgrd0;
    static doublereal term01, term11, term02, term12;


#define igrid_ref(a_1,a_2) pme_1.igrid[(a_2)*3 + a_1 - 4]
#define fuind_ref(a_1,a_2) fuind[(a_2)*3 + a_1]
#define qgrid_ref(a_1,a_2,a_3,a_4) pme_1.qgrid[(((a_4)*100 + (a_3))*100 + \
(a_2))*2 + a_1 - 20203]
#define fuinp_ref(a_1,a_2) fuinp[(a_2)*3 + a_1]
#define thetai1_ref(a_1,a_2,a_3) pme_1.thetai1[((a_3)*10 + (a_2))*4 + a_1 \
- 45]
#define thetai2_ref(a_1,a_2,a_3) pme_1.thetai2[((a_3)*10 + (a_2))*4 + a_1 \
- 45]
#define thetai3_ref(a_1,a_2,a_3) pme_1.thetai3[((a_3)*10 + (a_2))*4 + a_1 \
- 45]



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
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




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




/*     zero out the particle mesh Ewald charge grid */

    /* Parameter adjustments */
    fuinp -= 4;
    fuind -= 4;

    /* Function Body */
    i__1 = pme_1.nfft3;
    for (k = 1; k <= i__1; ++k) {
	i__2 = pme_1.nfft2;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = pme_1.nfft1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		qgrid_ref(1, i__, j, k) = 0.;
		qgrid_ref(2, i__, j, k) = 0.;
	    }
	}
    }

/*     put the induced dipole moments onto the grid */

    i__1 = mpole_1.npole;
    for (m = 1; m <= i__1; ++m) {
	igrd0 = igrid_ref(1, m);
	jgrd0 = igrid_ref(2, m);
	kgrd0 = igrid_ref(3, m);
	k0 = kgrd0;
	i__2 = pme_1.bsorder;
	for (it3 = 1; it3 <= i__2; ++it3) {
	    ++k0;
	    k = k0 + 1 + (pme_1.nfft3 - i_sign(&pme_1.nfft3, &k0)) / 2;
	    v0 = thetai3_ref(1, it3, m);
	    v1 = thetai3_ref(2, it3, m);
	    j0 = jgrd0;
	    i__3 = pme_1.bsorder;
	    for (it2 = 1; it2 <= i__3; ++it2) {
		++j0;
		j = j0 + 1 + (pme_1.nfft2 - i_sign(&pme_1.nfft2, &j0)) / 2;
		u0 = thetai2_ref(1, it2, m);
		u1 = thetai2_ref(2, it2, m);
		term01 = fuind_ref(2, m) * u1 * v0 + fuind_ref(3, m) * u0 * 
			v1;
		term11 = fuind_ref(1, m) * u0 * v0;
		term02 = fuinp_ref(2, m) * u1 * v0 + fuinp_ref(3, m) * u0 * 
			v1;
		term12 = fuinp_ref(1, m) * u0 * v0;
		i0 = igrd0;
		i__4 = pme_1.bsorder;
		for (it1 = 1; it1 <= i__4; ++it1) {
		    ++i0;
		    i__ = i0 + 1 + (pme_1.nfft1 - i_sign(&pme_1.nfft1, &i0)) /
			     2;
		    t0 = thetai1_ref(1, it1, m);
		    t1 = thetai1_ref(2, it1, m);
		    qgrid_ref(1, i__, j, k) = qgrid_ref(1, i__, j, k) + 
			    term01 * t0 + term11 * t1;
		    qgrid_ref(2, i__, j, k) = qgrid_ref(2, i__, j, k) + 
			    term02 * t0 + term12 * t1;
		}
	    }
	}
    }
    return 0;
} /* grid_uind__ */

#undef thetai3_ref
#undef thetai2_ref
#undef thetai1_ref
#undef fuinp_ref
#undef qgrid_ref
#undef fuind_ref
#undef igrid_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine fphi_mpole  --  multipole potential from grid  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "fphi_mpole" extracts the permanent multipole potential from */
/*     the particle mesh Ewald grid */


/* Subroutine */ int fphi_mpole__(doublereal *fphi)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer i_sign(integer *, integer *);

    /* Local variables */
    static integer i__, j, k, m, i0, j0, k0;
    static doublereal t0, u0, v0, v1, v2, v3, u1, u2, u3, t1, t2, t3, tq;
    static integer it1, it2, it3;
    static doublereal tu00, tu10, tu01, tu20, tu11, tu02, tu21, tu12, tu30, 
	    tu03;
    static integer igrd0, jgrd0, kgrd0;
    static doublereal tuv000, tuv100, tuv010, tuv001, tuv200, tuv020, tuv002, 
	    tuv110, tuv101, tuv011, tuv300, tuv030, tuv003, tuv210, tuv201, 
	    tuv120, tuv021, tuv102, tuv012, tuv111;


#define fphi_ref(a_1,a_2) fphi[(a_2)*20 + a_1]
#define igrid_ref(a_1,a_2) pme_1.igrid[(a_2)*3 + a_1 - 4]
#define qgrid_ref(a_1,a_2,a_3,a_4) pme_1.qgrid[(((a_4)*100 + (a_3))*100 + \
(a_2))*2 + a_1 - 20203]
#define thetai1_ref(a_1,a_2,a_3) pme_1.thetai1[((a_3)*10 + (a_2))*4 + a_1 \
- 45]
#define thetai2_ref(a_1,a_2,a_3) pme_1.thetai2[((a_3)*10 + (a_2))*4 + a_1 \
- 45]
#define thetai3_ref(a_1,a_2,a_3) pme_1.thetai3[((a_3)*10 + (a_2))*4 + a_1 \
- 45]



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
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




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




/*     extract the permanent multipole field at each site */

    /* Parameter adjustments */
    fphi -= 21;

    /* Function Body */
    i__1 = mpole_1.npole;
    for (m = 1; m <= i__1; ++m) {
	igrd0 = igrid_ref(1, m);
	jgrd0 = igrid_ref(2, m);
	kgrd0 = igrid_ref(3, m);
	tuv000 = 0.;
	tuv001 = 0.;
	tuv010 = 0.;
	tuv100 = 0.;
	tuv200 = 0.;
	tuv020 = 0.;
	tuv002 = 0.;
	tuv110 = 0.;
	tuv101 = 0.;
	tuv011 = 0.;
	tuv300 = 0.;
	tuv030 = 0.;
	tuv003 = 0.;
	tuv210 = 0.;
	tuv201 = 0.;
	tuv120 = 0.;
	tuv021 = 0.;
	tuv102 = 0.;
	tuv012 = 0.;
	tuv111 = 0.;
	k0 = kgrd0;
	i__2 = pme_1.bsorder;
	for (it3 = 1; it3 <= i__2; ++it3) {
	    ++k0;
	    k = k0 + 1 + (pme_1.nfft3 - i_sign(&pme_1.nfft3, &k0)) / 2;
	    v0 = thetai3_ref(1, it3, m);
	    v1 = thetai3_ref(2, it3, m);
	    v2 = thetai3_ref(3, it3, m);
	    v3 = thetai3_ref(4, it3, m);
	    tu00 = 0.;
	    tu10 = 0.;
	    tu01 = 0.;
	    tu20 = 0.;
	    tu11 = 0.;
	    tu02 = 0.;
	    tu30 = 0.;
	    tu21 = 0.;
	    tu12 = 0.;
	    tu03 = 0.;
	    j0 = jgrd0;
	    i__3 = pme_1.bsorder;
	    for (it2 = 1; it2 <= i__3; ++it2) {
		++j0;
		j = j0 + 1 + (pme_1.nfft2 - i_sign(&pme_1.nfft2, &j0)) / 2;
		u0 = thetai2_ref(1, it2, m);
		u1 = thetai2_ref(2, it2, m);
		u2 = thetai2_ref(3, it2, m);
		u3 = thetai2_ref(4, it2, m);
		t0 = 0.;
		t1 = 0.;
		t2 = 0.;
		t3 = 0.;
		i0 = igrd0;
		i__4 = pme_1.bsorder;
		for (it1 = 1; it1 <= i__4; ++it1) {
		    ++i0;
		    i__ = i0 + 1 + (pme_1.nfft1 - i_sign(&pme_1.nfft1, &i0)) /
			     2;
		    tq = qgrid_ref(1, i__, j, k);
		    t0 += tq * thetai1_ref(1, it1, m);
		    t1 += tq * thetai1_ref(2, it1, m);
		    t2 += tq * thetai1_ref(3, it1, m);
		    t3 += tq * thetai1_ref(4, it1, m);
		}
		tu00 += t0 * u0;
		tu10 += t1 * u0;
		tu01 += t0 * u1;
		tu20 += t2 * u0;
		tu11 += t1 * u1;
		tu02 += t0 * u2;
		tu30 += t3 * u0;
		tu21 += t2 * u1;
		tu12 += t1 * u2;
		tu03 += t0 * u3;
	    }
	    tuv000 += tu00 * v0;
	    tuv100 += tu10 * v0;
	    tuv010 += tu01 * v0;
	    tuv001 += tu00 * v1;
	    tuv200 += tu20 * v0;
	    tuv020 += tu02 * v0;
	    tuv002 += tu00 * v2;
	    tuv110 += tu11 * v0;
	    tuv101 += tu10 * v1;
	    tuv011 += tu01 * v1;
	    tuv300 += tu30 * v0;
	    tuv030 += tu03 * v0;
	    tuv003 += tu00 * v3;
	    tuv210 += tu21 * v0;
	    tuv201 += tu20 * v1;
	    tuv120 += tu12 * v0;
	    tuv021 += tu02 * v1;
	    tuv102 += tu10 * v2;
	    tuv012 += tu01 * v2;
	    tuv111 += tu11 * v1;
	}
	fphi_ref(1, m) = tuv000;
	fphi_ref(2, m) = tuv100;
	fphi_ref(3, m) = tuv010;
	fphi_ref(4, m) = tuv001;
	fphi_ref(5, m) = tuv200;
	fphi_ref(6, m) = tuv020;
	fphi_ref(7, m) = tuv002;
	fphi_ref(8, m) = tuv110;
	fphi_ref(9, m) = tuv101;
	fphi_ref(10, m) = tuv011;
	fphi_ref(11, m) = tuv300;
	fphi_ref(12, m) = tuv030;
	fphi_ref(13, m) = tuv003;
	fphi_ref(14, m) = tuv210;
	fphi_ref(15, m) = tuv201;
	fphi_ref(16, m) = tuv120;
	fphi_ref(17, m) = tuv021;
	fphi_ref(18, m) = tuv102;
	fphi_ref(19, m) = tuv012;
	fphi_ref(20, m) = tuv111;
    }
    return 0;
} /* fphi_mpole__ */

#undef thetai3_ref
#undef thetai2_ref
#undef thetai1_ref
#undef qgrid_ref
#undef igrid_ref
#undef fphi_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine fphi_uind  --  induced potential from grid  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "fphi_uind" extracts the induced dipole potential from */
/*     the particle mesh Ewald grid */


/* Subroutine */ int fphi_uind__(doublereal *fdip_phi1__, doublereal *
	fdip_phi2__, doublereal *fdip_sum_phi__)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer i_sign(integer *, integer *);

    /* Local variables */
    static integer i__, j, k, m, i0, j0, k0;
    static doublereal t0, u0, v0, v1, v2, v3, u1, u2, u3, t1, t2, t3;
    static integer it1, it2, it3;
    static doublereal t0_1__, t0_2__, t1_1__, t1_2__, t2_1__, t2_2__, tu00, 
	    tu10, tu01, tu20, tu11, tu02, tu30, tu21, tu12, tu03, tq_1__, 
	    tq_2__;
    static integer igrd0, jgrd0, kgrd0;
    static doublereal tu00_1__, tu01_1__, tu10_1__, tu00_2__, tu01_2__, 
	    tu10_2__, tu20_1__, tu11_1__, tu02_1__, tu20_2__, tu11_2__, 
	    tu02_2__, tuv000, tuv100, tuv010, tuv001, tuv200, tuv020, tuv002, 
	    tuv110, tuv101, tuv011, tuv300, tuv030, tuv003, tuv210, tuv201, 
	    tuv120, tuv021, tuv102, tuv012, tuv111, tuv100_1__, tuv010_1__, 
	    tuv001_1__, tuv100_2__, tuv010_2__, tuv001_2__, tuv200_1__, 
	    tuv020_1__, tuv002_1__, tuv110_1__, tuv101_1__, tuv011_1__, 
	    tuv200_2__, tuv020_2__, tuv002_2__, tuv110_2__, tuv101_2__, 
	    tuv011_2__;


#define fdip_phi1___ref(a_1,a_2) fdip_phi1__[(a_2)*10 + a_1]
#define fdip_phi2___ref(a_1,a_2) fdip_phi2__[(a_2)*10 + a_1]
#define fdip_sum_phi___ref(a_1,a_2) fdip_sum_phi__[(a_2)*20 + a_1]
#define igrid_ref(a_1,a_2) pme_1.igrid[(a_2)*3 + a_1 - 4]
#define qgrid_ref(a_1,a_2,a_3,a_4) pme_1.qgrid[(((a_4)*100 + (a_3))*100 + \
(a_2))*2 + a_1 - 20203]
#define thetai1_ref(a_1,a_2,a_3) pme_1.thetai1[((a_3)*10 + (a_2))*4 + a_1 \
- 45]
#define thetai2_ref(a_1,a_2,a_3) pme_1.thetai2[((a_3)*10 + (a_2))*4 + a_1 \
- 45]
#define thetai3_ref(a_1,a_2,a_3) pme_1.thetai3[((a_3)*10 + (a_2))*4 + a_1 \
- 45]



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
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




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




/*     extract the induced dipole field at each site */

    /* Parameter adjustments */
    fdip_sum_phi__ -= 21;
    fdip_phi2__ -= 11;
    fdip_phi1__ -= 11;

    /* Function Body */
    i__1 = mpole_1.npole;
    for (m = 1; m <= i__1; ++m) {
	igrd0 = igrid_ref(1, m);
	jgrd0 = igrid_ref(2, m);
	kgrd0 = igrid_ref(3, m);
	tuv100_1__ = 0.;
	tuv010_1__ = 0.;
	tuv001_1__ = 0.;
	tuv200_1__ = 0.;
	tuv020_1__ = 0.;
	tuv002_1__ = 0.;
	tuv110_1__ = 0.;
	tuv101_1__ = 0.;
	tuv011_1__ = 0.;
	tuv100_2__ = 0.;
	tuv010_2__ = 0.;
	tuv001_2__ = 0.;
	tuv200_2__ = 0.;
	tuv020_2__ = 0.;
	tuv002_2__ = 0.;
	tuv110_2__ = 0.;
	tuv101_2__ = 0.;
	tuv011_2__ = 0.;
	tuv000 = 0.;
	tuv001 = 0.;
	tuv010 = 0.;
	tuv100 = 0.;
	tuv200 = 0.;
	tuv020 = 0.;
	tuv002 = 0.;
	tuv110 = 0.;
	tuv101 = 0.;
	tuv011 = 0.;
	tuv300 = 0.;
	tuv030 = 0.;
	tuv003 = 0.;
	tuv210 = 0.;
	tuv201 = 0.;
	tuv120 = 0.;
	tuv021 = 0.;
	tuv102 = 0.;
	tuv012 = 0.;
	tuv111 = 0.;
	k0 = kgrd0;
	i__2 = pme_1.bsorder;
	for (it3 = 1; it3 <= i__2; ++it3) {
	    ++k0;
	    k = k0 + 1 + (pme_1.nfft3 - i_sign(&pme_1.nfft3, &k0)) / 2;
	    v0 = thetai3_ref(1, it3, m);
	    v1 = thetai3_ref(2, it3, m);
	    v2 = thetai3_ref(3, it3, m);
	    v3 = thetai3_ref(4, it3, m);
	    tu00_1__ = 0.;
	    tu01_1__ = 0.;
	    tu10_1__ = 0.;
	    tu20_1__ = 0.;
	    tu11_1__ = 0.;
	    tu02_1__ = 0.;
	    tu00_2__ = 0.;
	    tu01_2__ = 0.;
	    tu10_2__ = 0.;
	    tu20_2__ = 0.;
	    tu11_2__ = 0.;
	    tu02_2__ = 0.;
	    tu00 = 0.;
	    tu10 = 0.;
	    tu01 = 0.;
	    tu20 = 0.;
	    tu11 = 0.;
	    tu02 = 0.;
	    tu30 = 0.;
	    tu21 = 0.;
	    tu12 = 0.;
	    tu03 = 0.;
	    j0 = jgrd0;
	    i__3 = pme_1.bsorder;
	    for (it2 = 1; it2 <= i__3; ++it2) {
		++j0;
		j = j0 + 1 + (pme_1.nfft2 - i_sign(&pme_1.nfft2, &j0)) / 2;
		u0 = thetai2_ref(1, it2, m);
		u1 = thetai2_ref(2, it2, m);
		u2 = thetai2_ref(3, it2, m);
		u3 = thetai2_ref(4, it2, m);
		t0_1__ = 0.;
		t1_1__ = 0.;
		t2_1__ = 0.;
		t0_2__ = 0.;
		t1_2__ = 0.;
		t2_2__ = 0.;
		t3 = 0.;
		i0 = igrd0;
		i__4 = pme_1.bsorder;
		for (it1 = 1; it1 <= i__4; ++it1) {
		    ++i0;
		    i__ = i0 + 1 + (pme_1.nfft1 - i_sign(&pme_1.nfft1, &i0)) /
			     2;
		    tq_1__ = qgrid_ref(1, i__, j, k);
		    tq_2__ = qgrid_ref(2, i__, j, k);
		    t0_1__ += tq_1__ * thetai1_ref(1, it1, m);
		    t1_1__ += tq_1__ * thetai1_ref(2, it1, m);
		    t2_1__ += tq_1__ * thetai1_ref(3, it1, m);
		    t0_2__ += tq_2__ * thetai1_ref(1, it1, m);
		    t1_2__ += tq_2__ * thetai1_ref(2, it1, m);
		    t2_2__ += tq_2__ * thetai1_ref(3, it1, m);
		    t3 += (tq_1__ + tq_2__) * thetai1_ref(4, it1, m);
		}
		tu00_1__ += t0_1__ * u0;
		tu10_1__ += t1_1__ * u0;
		tu01_1__ += t0_1__ * u1;
		tu20_1__ += t2_1__ * u0;
		tu11_1__ += t1_1__ * u1;
		tu02_1__ += t0_1__ * u2;
		tu00_2__ += t0_2__ * u0;
		tu10_2__ += t1_2__ * u0;
		tu01_2__ += t0_2__ * u1;
		tu20_2__ += t2_2__ * u0;
		tu11_2__ += t1_2__ * u1;
		tu02_2__ += t0_2__ * u2;
		t0 = t0_1__ + t0_2__;
		t1 = t1_1__ + t1_2__;
		t2 = t2_1__ + t2_2__;
		tu00 += t0 * u0;
		tu10 += t1 * u0;
		tu01 += t0 * u1;
		tu20 += t2 * u0;
		tu11 += t1 * u1;
		tu02 += t0 * u2;
		tu30 += t3 * u0;
		tu21 += t2 * u1;
		tu12 += t1 * u2;
		tu03 += t0 * u3;
	    }
	    tuv100_1__ += tu10_1__ * v0;
	    tuv010_1__ += tu01_1__ * v0;
	    tuv001_1__ += tu00_1__ * v1;
	    tuv200_1__ += tu20_1__ * v0;
	    tuv020_1__ += tu02_1__ * v0;
	    tuv002_1__ += tu00_1__ * v2;
	    tuv110_1__ += tu11_1__ * v0;
	    tuv101_1__ += tu10_1__ * v1;
	    tuv011_1__ += tu01_1__ * v1;
	    tuv100_2__ += tu10_2__ * v0;
	    tuv010_2__ += tu01_2__ * v0;
	    tuv001_2__ += tu00_2__ * v1;
	    tuv200_2__ += tu20_2__ * v0;
	    tuv020_2__ += tu02_2__ * v0;
	    tuv002_2__ += tu00_2__ * v2;
	    tuv110_2__ += tu11_2__ * v0;
	    tuv101_2__ += tu10_2__ * v1;
	    tuv011_2__ += tu01_2__ * v1;
	    tuv000 += tu00 * v0;
	    tuv100 += tu10 * v0;
	    tuv010 += tu01 * v0;
	    tuv001 += tu00 * v1;
	    tuv200 += tu20 * v0;
	    tuv020 += tu02 * v0;
	    tuv002 += tu00 * v2;
	    tuv110 += tu11 * v0;
	    tuv101 += tu10 * v1;
	    tuv011 += tu01 * v1;
	    tuv300 += tu30 * v0;
	    tuv030 += tu03 * v0;
	    tuv003 += tu00 * v3;
	    tuv210 += tu21 * v0;
	    tuv201 += tu20 * v1;
	    tuv120 += tu12 * v0;
	    tuv021 += tu02 * v1;
	    tuv102 += tu10 * v2;
	    tuv012 += tu01 * v2;
	    tuv111 += tu11 * v1;
	}
	fdip_phi1___ref(2, m) = tuv100_1__;
	fdip_phi1___ref(3, m) = tuv010_1__;
	fdip_phi1___ref(4, m) = tuv001_1__;
	fdip_phi1___ref(5, m) = tuv200_1__;
	fdip_phi1___ref(6, m) = tuv020_1__;
	fdip_phi1___ref(7, m) = tuv002_1__;
	fdip_phi1___ref(8, m) = tuv110_1__;
	fdip_phi1___ref(9, m) = tuv101_1__;
	fdip_phi1___ref(10, m) = tuv011_1__;
	fdip_phi2___ref(2, m) = tuv100_2__;
	fdip_phi2___ref(3, m) = tuv010_2__;
	fdip_phi2___ref(4, m) = tuv001_2__;
	fdip_phi2___ref(5, m) = tuv200_2__;
	fdip_phi2___ref(6, m) = tuv020_2__;
	fdip_phi2___ref(7, m) = tuv002_2__;
	fdip_phi2___ref(8, m) = tuv110_2__;
	fdip_phi2___ref(9, m) = tuv101_2__;
	fdip_phi2___ref(10, m) = tuv011_2__;
	fdip_sum_phi___ref(1, m) = tuv000;
	fdip_sum_phi___ref(2, m) = tuv100;
	fdip_sum_phi___ref(3, m) = tuv010;
	fdip_sum_phi___ref(4, m) = tuv001;
	fdip_sum_phi___ref(5, m) = tuv200;
	fdip_sum_phi___ref(6, m) = tuv020;
	fdip_sum_phi___ref(7, m) = tuv002;
	fdip_sum_phi___ref(8, m) = tuv110;
	fdip_sum_phi___ref(9, m) = tuv101;
	fdip_sum_phi___ref(10, m) = tuv011;
	fdip_sum_phi___ref(11, m) = tuv300;
	fdip_sum_phi___ref(12, m) = tuv030;
	fdip_sum_phi___ref(13, m) = tuv003;
	fdip_sum_phi___ref(14, m) = tuv210;
	fdip_sum_phi___ref(15, m) = tuv201;
	fdip_sum_phi___ref(16, m) = tuv120;
	fdip_sum_phi___ref(17, m) = tuv021;
	fdip_sum_phi___ref(18, m) = tuv102;
	fdip_sum_phi___ref(19, m) = tuv012;
	fdip_sum_phi___ref(20, m) = tuv111;
    }
    return 0;
} /* fphi_uind__ */

#undef thetai3_ref
#undef thetai2_ref
#undef thetai1_ref
#undef qgrid_ref
#undef igrid_ref
#undef fdip_sum_phi___ref
#undef fdip_phi2___ref
#undef fdip_phi1___ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine cmp_to_fmp  --  transformation of multipoles  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "cmp_to_fmp" transforms the atomic multipoles from Cartesian */
/*     to fractional coordinates */


/* Subroutine */ int cmp_to_fmp__(doublereal *cmp, doublereal *fmp)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal ctf[100]	/* was [10][10] */;
    extern /* Subroutine */ int cart_to_frac__(doublereal *);


#define ctf_ref(a_1,a_2) ctf[(a_2)*10 + a_1 - 11]
#define cmp_ref(a_1,a_2) cmp[(a_2)*10 + a_1]
#define fmp_ref(a_1,a_2) fmp[(a_2)*10 + a_1]



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
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




/*     find the matrix to convert Cartesian to fractional */

    /* Parameter adjustments */
    fmp -= 11;
    cmp -= 11;

    /* Function Body */
    cart_to_frac__(ctf);

/*     apply the transformation to get the fractional multipoles */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	fmp_ref(1, i__) = ctf_ref(1, 1) * cmp_ref(1, i__);
	for (j = 2; j <= 4; ++j) {
	    fmp_ref(j, i__) = 0.;
	    for (k = 2; k <= 4; ++k) {
		fmp_ref(j, i__) = fmp_ref(j, i__) + ctf_ref(j, k) * cmp_ref(k,
			 i__);
	    }
	}
	for (j = 5; j <= 10; ++j) {
	    fmp_ref(j, i__) = 0.;
	    for (k = 5; k <= 10; ++k) {
		fmp_ref(j, i__) = fmp_ref(j, i__) + ctf_ref(j, k) * cmp_ref(k,
			 i__);
	    }
	}
    }
    return 0;
} /* cmp_to_fmp__ */

#undef fmp_ref
#undef cmp_ref
#undef ctf_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine cart_to_frac  --  Cartesian to fractional  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "cart_to_frac" computes a transformation matrix to convert */
/*     a multipole object in Cartesian coordinates to fractional */

/*     note the multipole components are stored in the condensed */
/*     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz) */


/* Subroutine */ int cart_to_frac__(doublereal *ctf)
{
    /* Initialized data */

    static integer qi1[6] = { 1,2,3,1,1,2 };
    static integer qi2[6] = { 1,2,3,2,3,3 };

    static doublereal a[9]	/* was [3][3] */;
    static integer i__, j, k, m, i1, i2;


#define a_ref(a_1,a_2) a[(a_2)*3 + a_1 - 4]
#define ctf_ref(a_1,a_2) ctf[(a_2)*10 + a_1]
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


    /* Parameter adjustments */
    ctf -= 11;

    /* Function Body */


/*     set the reciprocal vector transformation matrix */

    for (i__ = 1; i__ <= 3; ++i__) {
	a_ref(1, i__) = (doublereal) pme_1.nfft1 * recip_ref(i__, 1);
	a_ref(2, i__) = (doublereal) pme_1.nfft2 * recip_ref(i__, 2);
	a_ref(3, i__) = (doublereal) pme_1.nfft3 * recip_ref(i__, 3);
    }

/*     get the Cartesian to fractional conversion matrix */

    for (i__ = 1; i__ <= 10; ++i__) {
	for (j = 1; j <= 10; ++j) {
	    ctf_ref(j, i__) = 0.;
	}
    }
    ctf_ref(1, 1) = 1.;
    for (i__ = 2; i__ <= 4; ++i__) {
	for (j = 2; j <= 4; ++j) {
	    ctf_ref(i__, j) = a_ref(i__ - 1, j - 1);
	}
    }
    for (i1 = 1; i1 <= 3; ++i1) {
	k = qi1[i1 - 1];
	for (i2 = 1; i2 <= 6; ++i2) {
	    i__ = qi1[i2 - 1];
	    j = qi2[i2 - 1];
	    ctf_ref(i1 + 4, i2 + 4) = a_ref(k, i__) * a_ref(k, j);
	}
    }
    for (i1 = 4; i1 <= 6; ++i1) {
	k = qi1[i1 - 1];
	m = qi2[i1 - 1];
	for (i2 = 1; i2 <= 6; ++i2) {
	    i__ = qi1[i2 - 1];
	    j = qi2[i2 - 1];
	    ctf_ref(i1 + 4, i2 + 4) = a_ref(k, i__) * a_ref(m, j) + a_ref(k, 
		    j) * a_ref(m, i__);
	}
    }
    return 0;
} /* cart_to_frac__ */

#undef recip_ref
#undef ctf_ref
#undef a_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine fphi_to_cphi  --  transformation of potential  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "fphi_to_cphi" transforms the reciprocal space potential from */
/*     fractional to Cartesian coordinates */


/* Subroutine */ int fphi_to_cphi__(doublereal *fphi, doublereal *cphi)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal ftc[100]	/* was [10][10] */;
    extern /* Subroutine */ int frac_to_cart__(doublereal *);


#define ftc_ref(a_1,a_2) ftc[(a_2)*10 + a_1 - 11]
#define cphi_ref(a_1,a_2) cphi[(a_2)*10 + a_1]
#define fphi_ref(a_1,a_2) fphi[(a_2)*20 + a_1]



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
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




/*     find the matrix to convert fractional to Cartesian */

    /* Parameter adjustments */
    cphi -= 11;
    fphi -= 21;

    /* Function Body */
    frac_to_cart__(ftc);

/*     apply the transformation to get the Cartesian potential */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cphi_ref(1, i__) = ftc_ref(1, 1) * fphi_ref(1, i__);
	for (j = 2; j <= 4; ++j) {
	    cphi_ref(j, i__) = 0.;
	    for (k = 2; k <= 4; ++k) {
		cphi_ref(j, i__) = cphi_ref(j, i__) + ftc_ref(j, k) * 
			fphi_ref(k, i__);
	    }
	}
	for (j = 5; j <= 10; ++j) {
	    cphi_ref(j, i__) = 0.;
	    for (k = 5; k <= 10; ++k) {
		cphi_ref(j, i__) = cphi_ref(j, i__) + ftc_ref(j, k) * 
			fphi_ref(k, i__);
	    }
	}
    }
    return 0;
} /* fphi_to_cphi__ */

#undef fphi_ref
#undef cphi_ref
#undef ftc_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine frac_to_cart  --  fractional to Cartesian  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "frac_to_cart" computes a transformation matrix to convert */
/*     a multipole object in fraction coordinates to Cartesian */

/*     note the multipole components are stored in the condensed */
/*     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz) */


/* Subroutine */ int frac_to_cart__(doublereal *ftc)
{
    /* Initialized data */

    static integer qi1[6] = { 1,2,3,1,1,2 };
    static integer qi2[6] = { 1,2,3,2,3,3 };

    static doublereal a[9]	/* was [3][3] */;
    static integer i__, j, k, m, i1, i2;


#define a_ref(a_1,a_2) a[(a_2)*3 + a_1 - 4]
#define ftc_ref(a_1,a_2) ftc[(a_2)*10 + a_1]
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


    /* Parameter adjustments */
    ftc -= 11;

    /* Function Body */


/*     set the reciprocal vector transformation matrix */

    for (i__ = 1; i__ <= 3; ++i__) {
	a_ref(i__, 1) = (doublereal) pme_1.nfft1 * recip_ref(i__, 1);
	a_ref(i__, 2) = (doublereal) pme_1.nfft2 * recip_ref(i__, 2);
	a_ref(i__, 3) = (doublereal) pme_1.nfft3 * recip_ref(i__, 3);
    }

/*     get the fractional to Cartesian conversion matrix */

    for (i__ = 1; i__ <= 10; ++i__) {
	for (j = 1; j <= 10; ++j) {
	    ftc_ref(j, i__) = 0.;
	}
    }
    ftc_ref(1, 1) = 1.;
    for (i__ = 2; i__ <= 4; ++i__) {
	for (j = 2; j <= 4; ++j) {
	    ftc_ref(i__, j) = a_ref(i__ - 1, j - 1);
	}
    }
    for (i1 = 1; i1 <= 3; ++i1) {
	k = qi1[i1 - 1];
	for (i2 = 1; i2 <= 3; ++i2) {
	    i__ = qi1[i2 - 1];
	    ftc_ref(i1 + 4, i2 + 4) = a_ref(k, i__) * a_ref(k, i__);
	}
	for (i2 = 4; i2 <= 6; ++i2) {
	    i__ = qi1[i2 - 1];
	    j = qi2[i2 - 1];
	    ftc_ref(i1 + 4, i2 + 4) = a_ref(k, i__) * 2. * a_ref(k, j);
	}
    }
    for (i1 = 4; i1 <= 6; ++i1) {
	k = qi1[i1 - 1];
	m = qi2[i1 - 1];
	for (i2 = 1; i2 <= 3; ++i2) {
	    i__ = qi1[i2 - 1];
	    ftc_ref(i1 + 4, i2 + 4) = a_ref(k, i__) * a_ref(m, i__);
	}
	for (i2 = 4; i2 <= 6; ++i2) {
	    i__ = qi1[i2 - 1];
	    j = qi2[i2 - 1];
	    ftc_ref(i1 + 4, i2 + 4) = a_ref(k, i__) * a_ref(m, j) + a_ref(m, 
		    i__) * a_ref(k, j);
	}
    }
    return 0;
} /* frac_to_cart__ */

#undef recip_ref
#undef ftc_ref
#undef a_ref


