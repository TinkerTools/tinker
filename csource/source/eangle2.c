/* eangle2.f -- translated by f2c (version 20050501).
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
    doublereal ak[75000], anat[75000], afld[75000];
    integer nangle, iang[300000]	/* was [4][75000] */;
} angle_;

#define angle_1 angle_

struct {
    doublereal angunit, stbnunit, aaunit, opbunit, opdunit, cang, qang, pang, 
	    sang, copb, qopb, popb, sopb, copd, qopd, popd, sopd;
    char angtyp[600000], opbtyp[8];
} angpot_;

#define angpot_1 angpot_

struct {
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

struct {
    doublereal desum[75000]	/* was [3][25000] */, deb[75000]	/* 
	    was [3][25000] */, dea[75000]	/* was [3][25000] */, deba[
	    75000]	/* was [3][25000] */, deub[75000]	/* was [3][
	    25000] */, deaa[75000]	/* was [3][25000] */, deopb[75000]	
	    /* was [3][25000] */, deopd[75000]	/* was [3][25000] */, deid[
	    75000]	/* was [3][25000] */, deit[75000]	/* was [3][
	    25000] */, det[75000]	/* was [3][25000] */, dept[75000]	
	    /* was [3][25000] */, debt[75000]	/* was [3][25000] */, dett[
	    75000]	/* was [3][25000] */, dev[75000]	/* was [3][
	    25000] */, dec[75000]	/* was [3][25000] */, decd[75000]	
	    /* was [3][25000] */, ded[75000]	/* was [3][25000] */, dem[
	    75000]	/* was [3][25000] */, dep[75000]	/* was [3][
	    25000] */, der[75000]	/* was [3][25000] */, des[75000]	
	    /* was [3][25000] */, delf[75000]	/* was [3][25000] */, deg[
	    75000]	/* was [3][25000] */, dex[75000]	/* was [3][
	    25000] */;
} deriv_;

#define deriv_1 deriv_

struct {
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    doublereal hessx[75000]	/* was [3][25000] */, hessy[75000]	/* 
	    was [3][25000] */, hessz[75000]	/* was [3][25000] */;
} hessn_;

#define hessn_1 hessn_

struct {
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine eangle2  --  atom-by-atom angle bend Hessian  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "eangle2" calculates second derivatives of the angle bending */
/*     energy for a single atom using a mixture of analytical and */
/*     finite difference methods; projected in-plane angles at trigonal */
/*     centers, special linear or Fourier angle bending terms are */
/*     optionally used */


/* Subroutine */ int eangle2_(integer *i__)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static logical twosided;
    static integer j, k;
    static doublereal d0[75000]	/* was [3][25000] */;
    static integer ia, ib, ic, id;
    static doublereal old, eps, fgrp, term;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;
    extern /* Subroutine */ int eangle2a_(integer *), eangle2b_(integer *);


#define d0_ref(a_1,a_2) d0[(a_2)*3 + a_1 - 4]
#define dea_ref(a_1,a_2) deriv_1.dea[(a_2)*3 + a_1 - 4]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define hessx_ref(a_1,a_2) hessn_1.hessx[(a_2)*3 + a_1 - 4]
#define hessy_ref(a_1,a_2) hessn_1.hessy[(a_2)*3 + a_1 - 4]
#define hessz_ref(a_1,a_2) hessn_1.hessz[(a_2)*3 + a_1 - 4]
#define angtyp_ref(a_0,a_1) &angpot_1.angtyp[(a_1)*8 + a_0 - 8]



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
/*     ##  angle.i  --  bond angles within the current structure  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     ak       harmonic angle force constant (kcal/mole/rad**2) */
/*     anat     ideal bond angle or phase shift angle (degrees) */
/*     afld     periodicity for Fourier bond angle term */
/*     nangle   total number of bond angles in the system */
/*     iang     numbers of the atoms in each bond angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  angpot.i  --  specifics of angle bend functional forms  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     angunit    convert angle bending energy to kcal/mole */
/*     stbnunit   convert stretch-bend energy to kcal/mole */
/*     aaunit     convert angle-angle energy to kcal/mole */
/*     opbunit    convert out-of-plane bend energy to kcal/mole */
/*     opdunit    convert out-of-plane distance energy to kcal/mole */
/*     cang       cubic coefficient in angle bending potential */
/*     qang       quartic coefficient in angle bending potential */
/*     pang       quintic coefficient in angle bending potential */
/*     sang       sextic coefficient in angle bending potential */
/*     copb       cubic coefficient in out-of-plane bend potential */
/*     qopb       quartic coefficient in out-of-plane bend potential */
/*     popb       quintic coefficient in out-of-plane bend potential */
/*     sopb       sextic coefficient in out-of-plane bend potential */
/*     copd       cubic coefficient in out-of-plane distance potential */
/*     qopd       quartic coefficient in out-of-plane distance potential */
/*     popd       quintic coefficient in out-of-plane distance potential */
/*     sopd       sextic coefficient in out-of-plane distance potential */
/*     angtyp     type of angle bending function for each bond angle */
/*     opbtyp     type of out-of-plane bend potential energy function */




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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  deriv.i  --  Cartesian coordinate derivative components  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     desum   total energy Cartesian coordinate derivatives */
/*     deb     bond stretch Cartesian coordinate derivatives */
/*     dea     angle bend Cartesian coordinate derivatives */
/*     deba    stretch-bend Cartesian coordinate derivatives */
/*     deub    Urey-Bradley Cartesian coordinate derivatives */
/*     deaa    angle-angle Cartesian coordinate derivatives */
/*     deopb   out-of-plane bend Cartesian coordinate derivatives */
/*     deopd   out-of-plane distance Cartesian coordinate derivatives */
/*     deid    improper dihedral Cartesian coordinate derivatives */
/*     deit    improper torsion Cartesian coordinate derivatives */
/*     det     torsional Cartesian coordinate derivatives */
/*     dept    pi-orbital torsion Cartesian coordinate derivatives */
/*     debt    stretch-torsion Cartesian coordinate derivatives */
/*     dett    torsion-torsion Cartesian coordinate derivatives */
/*     dev     van der Waals Cartesian coordinate derivatives */
/*     dec     charge-charge Cartesian coordinate derivatives */
/*     decd    charge-dipole Cartesian coordinate derivatives */
/*     ded     dipole-dipole Cartesian coordinate derivatives */
/*     dem     multipole Cartesian coordinate derivatives */
/*     dep     polarization Cartesian coordinate derivatives */
/*     der     reaction field Cartesian coordinate derivatives */
/*     des     solvation Cartesian coordinate derivatives */
/*     delf    metal ligand field Cartesian coordinate derivatives */
/*     deg     geometric restraint Cartesian coordinate derivatives */
/*     dex     extra energy term Cartesian coordinate derivatives */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  group.i  --  partitioning of system into atom groups  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     grpmass     total mass of all the atoms in each group */
/*     wgrp        weight for each set of group-group interactions */
/*     ngrp        total number of atom groups in the system */
/*     kgrp        contiguous list of the atoms in each group */
/*     igrp        first and last atom of each group in the list */
/*     grplist     number of the group to which each atom belongs */
/*     use_group   flag to use partitioning of system into groups */
/*     use_intra   flag to include only intragroup interactions */
/*     use_inter   flag to include only intergroup interactions */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  hessn.i  --  Cartesian Hessian elements for a single atom  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     hessx   Hessian elements for x-component of current atom */
/*     hessy   Hessian elements for y-component of current atom */
/*     hessz   Hessian elements for z-component of current atom */




/*     compute analytical angle bending Hessian elements */

    eangle2a_(i__);

/*     set stepsize for derivatives and default group weight */

    eps = 1e-5;
    fgrp = 1.;
    twosided = FALSE_;
    if (atoms_1.n <= 50) {
	twosided = TRUE_;
    }

/*     compute numerical in-plane bend Hessian for current atom */

    i__1 = angle_1.nangle;
    for (k = 1; k <= i__1; ++k) {
	proceed = FALSE_;
	if (s_cmp(angtyp_ref(0, k), "IN-PLANE", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = iang_ref(1, k);
	    ib = iang_ref(2, k);
	    ic = iang_ref(3, k);
	    id = iang_ref(4, k);
	    proceed = *i__ == ia || *i__ == ib || *i__ == ic || *i__ == id;
	    if (proceed && group_1.use_group__) {
		groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &c__0, &c__0);
	    }
	}
	if (proceed) {
	    term = fgrp / eps;

/*     find first derivatives for the base structure */

	    if (! twosided) {
		eangle2b_(&k);
		for (j = 1; j <= 3; ++j) {
		    d0_ref(j, ia) = dea_ref(j, ia);
		    d0_ref(j, ib) = dea_ref(j, ib);
		    d0_ref(j, ic) = dea_ref(j, ic);
		    d0_ref(j, id) = dea_ref(j, id);
		}
	    }

/*     find numerical x-components via perturbed structures */

	    old = atoms_1.x[*i__ - 1];
	    if (twosided) {
		atoms_1.x[*i__ - 1] -= eps * .5;
		eangle2b_(&k);
		for (j = 1; j <= 3; ++j) {
		    d0_ref(j, ia) = dea_ref(j, ia);
		    d0_ref(j, ib) = dea_ref(j, ib);
		    d0_ref(j, ic) = dea_ref(j, ic);
		    d0_ref(j, id) = dea_ref(j, id);
		}
	    }
	    atoms_1.x[*i__ - 1] += eps;
	    eangle2b_(&k);
	    atoms_1.x[*i__ - 1] = old;
	    for (j = 1; j <= 3; ++j) {
		hessx_ref(j, ia) = hessx_ref(j, ia) + term * (dea_ref(j, ia) 
			- d0_ref(j, ia));
		hessx_ref(j, ib) = hessx_ref(j, ib) + term * (dea_ref(j, ib) 
			- d0_ref(j, ib));
		hessx_ref(j, ic) = hessx_ref(j, ic) + term * (dea_ref(j, ic) 
			- d0_ref(j, ic));
		hessx_ref(j, id) = hessx_ref(j, id) + term * (dea_ref(j, id) 
			- d0_ref(j, id));
	    }

/*     find numerical y-components via perturbed structures */

	    old = atoms_1.y[*i__ - 1];
	    if (twosided) {
		atoms_1.y[*i__ - 1] -= eps * .5;
		eangle2b_(&k);
		for (j = 1; j <= 3; ++j) {
		    d0_ref(j, ia) = dea_ref(j, ia);
		    d0_ref(j, ib) = dea_ref(j, ib);
		    d0_ref(j, ic) = dea_ref(j, ic);
		    d0_ref(j, id) = dea_ref(j, id);
		}
	    }
	    atoms_1.y[*i__ - 1] += eps;
	    eangle2b_(&k);
	    atoms_1.y[*i__ - 1] = old;
	    for (j = 1; j <= 3; ++j) {
		hessy_ref(j, ia) = hessy_ref(j, ia) + term * (dea_ref(j, ia) 
			- d0_ref(j, ia));
		hessy_ref(j, ib) = hessy_ref(j, ib) + term * (dea_ref(j, ib) 
			- d0_ref(j, ib));
		hessy_ref(j, ic) = hessy_ref(j, ic) + term * (dea_ref(j, ic) 
			- d0_ref(j, ic));
		hessy_ref(j, id) = hessy_ref(j, id) + term * (dea_ref(j, id) 
			- d0_ref(j, id));
	    }

/*     find numerical z-components via perturbed structures */

	    old = atoms_1.z__[*i__ - 1];
	    if (twosided) {
		atoms_1.z__[*i__ - 1] -= eps * .5;
		eangle2b_(&k);
		for (j = 1; j <= 3; ++j) {
		    d0_ref(j, ia) = dea_ref(j, ia);
		    d0_ref(j, ib) = dea_ref(j, ib);
		    d0_ref(j, ic) = dea_ref(j, ic);
		    d0_ref(j, id) = dea_ref(j, id);
		}
	    }
	    atoms_1.z__[*i__ - 1] += eps;
	    eangle2b_(&k);
	    atoms_1.z__[*i__ - 1] = old;
	    for (j = 1; j <= 3; ++j) {
		hessz_ref(j, ia) = hessz_ref(j, ia) + term * (dea_ref(j, ia) 
			- d0_ref(j, ia));
		hessz_ref(j, ib) = hessz_ref(j, ib) + term * (dea_ref(j, ib) 
			- d0_ref(j, ib));
		hessz_ref(j, ic) = hessz_ref(j, ic) + term * (dea_ref(j, ic) 
			- d0_ref(j, ic));
		hessz_ref(j, id) = hessz_ref(j, id) + term * (dea_ref(j, id) 
			- d0_ref(j, id));
	    }
	}
    }
    return 0;
} /* eangle2_ */

#undef angtyp_ref
#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef iang_ref
#undef dea_ref
#undef d0_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine eangle2a  --  angle bending Hessian; analytical  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "eangle2a" calculates bond angle bending potential energy */
/*     second derivatives with respect to Cartesian coordinates */


/* Subroutine */ int eangle2a_(integer *iatom)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal), acos(doublereal), cos(doublereal), sin(
	    doublereal);

    /* Local variables */
    static doublereal dziaxic, dziayic, dziazic;
    static logical proceed;
    static integer i__, ia, ib, ic;
    static doublereal dt, rp, xp, yp, zp, dt2, dt3, dt4, rp2, xab, yab, zab, 
	    xcb, ycb, zcb, xia, yia, zia, xib, yib, dot, zib, xic, yic, zic, 
	    xpo, ypo, zpo, rab2, rcb2, fold, xabp, yabp, xrab, yrab, sine, 
	    fgrp, zrab, xrcb, yrcb, zrcb, zabp, xcbp, ycbp, zcbp, ideal;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal deddt, angle, force, terma, termc, d2eddt2;
    static logical linear;
    static doublereal factor, cosine;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal ddtdxia, ddtdyia, ddtdzia, ddtdxib, ddtdyib, ddtdzib, 
	    ddtdxic, ddtdyic, ddtdzic, dxiaxia, dxiayia, dxiazia, dxibxib, 
	    dxibyib, dxibzib, dxicxic, dxicyic, dxiczic, dyiayia, dyiazia, 
	    dziazia, dyibyib, dyibzib, dzibzib, dyicyic, dyiczic, dziczic, 
	    dxibxia, dxibyia, dxibzia, dyibxia, dyibyia, dyibzia, dzibxia, 
	    dzibyia, dzibzia, dxibxic, dxibyic, dxibzic, dyibxic, dyibyic, 
	    dyibzic, dzibxic, dzibyic, dzibzic, dxiaxic, dxiayic, dxiazic, 
	    dyiaxic, dyiayic, dyiazic;


#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define hessx_ref(a_1,a_2) hessn_1.hessx[(a_2)*3 + a_1 - 4]
#define hessy_ref(a_1,a_2) hessn_1.hessy[(a_2)*3 + a_1 - 4]
#define hessz_ref(a_1,a_2) hessn_1.hessz[(a_2)*3 + a_1 - 4]
#define angtyp_ref(a_0,a_1) &angpot_1.angtyp[(a_1)*8 + a_0 - 8]



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
/*     ##  angle.i  --  bond angles within the current structure  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     ak       harmonic angle force constant (kcal/mole/rad**2) */
/*     anat     ideal bond angle or phase shift angle (degrees) */
/*     afld     periodicity for Fourier bond angle term */
/*     nangle   total number of bond angles in the system */
/*     iang     numbers of the atoms in each bond angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  angpot.i  --  specifics of angle bend functional forms  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     angunit    convert angle bending energy to kcal/mole */
/*     stbnunit   convert stretch-bend energy to kcal/mole */
/*     aaunit     convert angle-angle energy to kcal/mole */
/*     opbunit    convert out-of-plane bend energy to kcal/mole */
/*     opdunit    convert out-of-plane distance energy to kcal/mole */
/*     cang       cubic coefficient in angle bending potential */
/*     qang       quartic coefficient in angle bending potential */
/*     pang       quintic coefficient in angle bending potential */
/*     sang       sextic coefficient in angle bending potential */
/*     copb       cubic coefficient in out-of-plane bend potential */
/*     qopb       quartic coefficient in out-of-plane bend potential */
/*     popb       quintic coefficient in out-of-plane bend potential */
/*     sopb       sextic coefficient in out-of-plane bend potential */
/*     copd       cubic coefficient in out-of-plane distance potential */
/*     qopd       quartic coefficient in out-of-plane distance potential */
/*     popd       quintic coefficient in out-of-plane distance potential */
/*     sopd       sextic coefficient in out-of-plane distance potential */
/*     angtyp     type of angle bending function for each bond angle */
/*     opbtyp     type of out-of-plane bend potential energy function */




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
/*     ##  bound.i  --  control of periodic boundary conditions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     polycut       cutoff distance for infinite polymer nonbonds */
/*     polycut2      square of infinite polymer nonbond cutoff */
/*     use_bounds    flag to use periodic boundary conditions */
/*     use_replica   flag to use replicates for periodic system */
/*     use_polymer   flag to mark presence of infinite polymer */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  group.i  --  partitioning of system into atom groups  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     grpmass     total mass of all the atoms in each group */
/*     wgrp        weight for each set of group-group interactions */
/*     ngrp        total number of atom groups in the system */
/*     kgrp        contiguous list of the atoms in each group */
/*     igrp        first and last atom of each group in the list */
/*     grplist     number of the group to which each atom belongs */
/*     use_group   flag to use partitioning of system into groups */
/*     use_intra   flag to include only intragroup interactions */
/*     use_inter   flag to include only intergroup interactions */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  hessn.i  --  Cartesian Hessian elements for a single atom  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     hessx   Hessian elements for x-component of current atom */
/*     hessy   Hessian elements for y-component of current atom */
/*     hessz   Hessian elements for z-component of current atom */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  math.i  --  mathematical and geometrical constants  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     radian   conversion factor from radians to degrees */
/*     pi       numerical value of the geometric constant */
/*     sqrtpi   numerical value of the square root of Pi */
/*     logten   numerical value of the natural log of ten */
/*     sqrttwo  numerical value of the square root of two */
/*     twosix   numerical value of the sixth root of two */




/*     calculate the bond angle bending energy term */

    i__1 = angle_1.nangle;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	ideal = angle_1.anat[i__ - 1];
	force = angle_1.ak[i__ - 1];

/*     decide whether to compute the current interaction */

	proceed = *iatom == ia || *iatom == ib || *iatom == ic;
	if (proceed && group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &c__0, &c__0, &c__0);
	}

/*     get the coordinates of the atoms in the angle */

	if (proceed) {
	    xia = atoms_1.x[ia - 1];
	    yia = atoms_1.y[ia - 1];
	    zia = atoms_1.z__[ia - 1];
	    xib = atoms_1.x[ib - 1];
	    yib = atoms_1.y[ib - 1];
	    zib = atoms_1.z__[ib - 1];
	    xic = atoms_1.x[ic - 1];
	    yic = atoms_1.y[ic - 1];
	    zic = atoms_1.z__[ic - 1];

/*     compute the bond angle bending Hessian elements */

	    if (s_cmp(angtyp_ref(0, i__), "IN-PLANE", (ftnlen)8, (ftnlen)8) !=
		     0) {
		xab = xia - xib;
		yab = yia - yib;
		zab = zia - zib;
		xcb = xic - xib;
		ycb = yic - yib;
		zcb = zic - zib;
		if (bound_1.use_polymer__) {
		    image_(&xab, &yab, &zab);
		    image_(&xcb, &ycb, &zcb);
		}
		rab2 = xab * xab + yab * yab + zab * zab;
		rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
		if (rab2 != 0. && rcb2 != 0.) {
		    xp = ycb * zab - zcb * yab;
		    yp = zcb * xab - xcb * zab;
		    zp = xcb * yab - ycb * xab;
		    rp = sqrt(xp * xp + yp * yp + zp * zp);
		    dot = xab * xcb + yab * ycb + zab * zcb;
		    cosine = dot / sqrt(rab2 * rcb2);
/* Computing MIN */
		    d__1 = 1., d__2 = max(-1.,cosine);
		    cosine = min(d__1,d__2);
		    angle = acos(cosine) * 57.29577951308232088;

/*     get the master chain rule terms for derivatives */

		    if (s_cmp(angtyp_ref(0, i__), "HARMONIC", (ftnlen)8, (
			    ftnlen)8) == 0) {
			dt = angle - ideal;
			dt2 = dt * dt;
			dt3 = dt2 * dt;
			dt4 = dt3 * dt;
			deddt = angpot_1.angunit * force * dt * 
				57.29577951308232088 * (angpot_1.cang * 3. * 
				dt + 2. + angpot_1.qang * 4. * dt2 + 
				angpot_1.pang * 5. * dt3 + angpot_1.sang * 6. 
				* dt4);
			d2eddt2 = angpot_1.angunit * force * 
				3282.8063500117441 * (angpot_1.cang * 6. * dt 
				+ 2. + angpot_1.qang * 12. * dt2 + 
				angpot_1.pang * 20. * dt3 + angpot_1.sang * 
				30. * dt4);
		    } else if (s_cmp(angtyp_ref(0, i__), "LINEAR", (ftnlen)8, 
			    (ftnlen)6) == 0) {
			factor = angpot_1.angunit * 2. * 3282.8063500117441;
			sine = sqrt(1. - cosine * cosine);
			deddt = -factor * force * sine;
			d2eddt2 = -factor * force * cosine;
		    } else if (s_cmp(angtyp_ref(0, i__), "FOURIER", (ftnlen)8,
			     (ftnlen)7) == 0) {
			fold = angle_1.afld[i__ - 1];
			factor = angpot_1.angunit * 2. * (3282.8063500117441 /
				 fold);
			cosine = cos((fold * angle - ideal) / 
				57.29577951308232088);
			sine = sin((fold * angle - ideal) / 
				57.29577951308232088);
			deddt = -factor * force * sine;
			d2eddt2 = -factor * force * fold * cosine;
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			deddt *= fgrp;
			d2eddt2 *= fgrp;
		    }

/*     construct an orthogonal direction for linear angles */

		    linear = FALSE_;
		    if (rp < 1e-6) {
			linear = TRUE_;
			if (xab != 0. && yab != 0.) {
			    xp = -yab;
			    yp = xab;
			    zp = 0.;
			} else if (xab == 0. && yab == 0.) {
			    xp = 1.;
			    yp = 0.;
			    zp = 0.;
			} else if (xab != 0. && yab == 0.) {
			    xp = 0.;
			    yp = 1.;
			    zp = 0.;
			} else if (xab == 0. && yab != 0.) {
			    xp = 1.;
			    yp = 0.;
			    zp = 0.;
			}
			rp = sqrt(xp * xp + yp * yp + zp * zp);
		    }

/*     first derivatives of bond angle with respect to coordinates */

L10:
		    terma = -1. / (rab2 * rp);
		    termc = 1. / (rcb2 * rp);
		    ddtdxia = terma * (yab * zp - zab * yp);
		    ddtdyia = terma * (zab * xp - xab * zp);
		    ddtdzia = terma * (xab * yp - yab * xp);
		    ddtdxic = termc * (ycb * zp - zcb * yp);
		    ddtdyic = termc * (zcb * xp - xcb * zp);
		    ddtdzic = termc * (xcb * yp - ycb * xp);
		    ddtdxib = -ddtdxia - ddtdxic;
		    ddtdyib = -ddtdyia - ddtdyic;
		    ddtdzib = -ddtdzia - ddtdzic;

/*     abbreviations used in defining chain rule terms */

		    xrab = xab * 2. / rab2;
		    yrab = yab * 2. / rab2;
		    zrab = zab * 2. / rab2;
		    xrcb = xcb * 2. / rcb2;
		    yrcb = ycb * 2. / rcb2;
		    zrcb = zcb * 2. / rcb2;
		    rp2 = 1. / (rp * rp);
		    xabp = (yab * zp - zab * yp) * rp2;
		    yabp = (zab * xp - xab * zp) * rp2;
		    zabp = (xab * yp - yab * xp) * rp2;
		    xcbp = (ycb * zp - zcb * yp) * rp2;
		    ycbp = (zcb * xp - xcb * zp) * rp2;
		    zcbp = (xcb * yp - ycb * xp) * rp2;

/*     chain rule terms for second derivative components */

		    dxiaxia = terma * (xab * xcb - dot) + ddtdxia * (xcbp - 
			    xrab);
		    dxiayia = terma * (zp + yab * xcb) + ddtdxia * (ycbp - 
			    yrab);
		    dxiazia = terma * (zab * xcb - yp) + ddtdxia * (zcbp - 
			    zrab);
		    dyiayia = terma * (yab * ycb - dot) + ddtdyia * (ycbp - 
			    yrab);
		    dyiazia = terma * (xp + zab * ycb) + ddtdyia * (zcbp - 
			    zrab);
		    dziazia = terma * (zab * zcb - dot) + ddtdzia * (zcbp - 
			    zrab);
		    dxicxic = termc * (dot - xab * xcb) - ddtdxic * (xabp + 
			    xrcb);
		    dxicyic = termc * (zp - ycb * xab) - ddtdxic * (yabp + 
			    yrcb);
		    dxiczic = -termc * (yp + zcb * xab) - ddtdxic * (zabp + 
			    zrcb);
		    dyicyic = termc * (dot - yab * ycb) - ddtdyic * (yabp + 
			    yrcb);
		    dyiczic = termc * (xp - zcb * yab) - ddtdyic * (zabp + 
			    zrcb);
		    dziczic = termc * (dot - zab * zcb) - ddtdzic * (zabp + 
			    zrcb);
		    dxiaxic = terma * (yab * yab + zab * zab) - ddtdxia * 
			    xabp;
		    dxiayic = -terma * xab * yab - ddtdxia * yabp;
		    dxiazic = -terma * xab * zab - ddtdxia * zabp;
		    dyiaxic = -terma * xab * yab - ddtdyia * xabp;
		    dyiayic = terma * (xab * xab + zab * zab) - ddtdyia * 
			    yabp;
		    dyiazic = -terma * yab * zab - ddtdyia * zabp;
		    dziaxic = -terma * xab * zab - ddtdzia * xabp;
		    dziayic = -terma * yab * zab - ddtdzia * yabp;
		    dziazic = terma * (xab * xab + yab * yab) - ddtdzia * 
			    zabp;

/*     get some second derivative chain rule terms by difference */

		    dxibxia = -dxiaxia - dxiaxic;
		    dxibyia = -dxiayia - dyiaxic;
		    dxibzia = -dxiazia - dziaxic;
		    dyibxia = -dxiayia - dxiayic;
		    dyibyia = -dyiayia - dyiayic;
		    dyibzia = -dyiazia - dziayic;
		    dzibxia = -dxiazia - dxiazic;
		    dzibyia = -dyiazia - dyiazic;
		    dzibzia = -dziazia - dziazic;
		    dxibxic = -dxicxic - dxiaxic;
		    dxibyic = -dxicyic - dxiayic;
		    dxibzic = -dxiczic - dxiazic;
		    dyibxic = -dxicyic - dyiaxic;
		    dyibyic = -dyicyic - dyiayic;
		    dyibzic = -dyiczic - dyiazic;
		    dzibxic = -dxiczic - dziaxic;
		    dzibyic = -dyiczic - dziayic;
		    dzibzic = -dziczic - dziazic;
		    dxibxib = -dxibxia - dxibxic;
		    dxibyib = -dxibyia - dxibyic;
		    dxibzib = -dxibzia - dxibzic;
		    dyibyib = -dyibyia - dyibyic;
		    dyibzib = -dyibzia - dyibzic;
		    dzibzib = -dzibzia - dzibzic;

/*     increment diagonal and off-diagonal Hessian elements */

		    if (ia == *iatom) {
			hessx_ref(1, ia) = hessx_ref(1, ia) + deddt * dxiaxia 
				+ d2eddt2 * ddtdxia * ddtdxia;
			hessx_ref(2, ia) = hessx_ref(2, ia) + deddt * dxiayia 
				+ d2eddt2 * ddtdxia * ddtdyia;
			hessx_ref(3, ia) = hessx_ref(3, ia) + deddt * dxiazia 
				+ d2eddt2 * ddtdxia * ddtdzia;
			hessy_ref(1, ia) = hessy_ref(1, ia) + deddt * dxiayia 
				+ d2eddt2 * ddtdyia * ddtdxia;
			hessy_ref(2, ia) = hessy_ref(2, ia) + deddt * dyiayia 
				+ d2eddt2 * ddtdyia * ddtdyia;
			hessy_ref(3, ia) = hessy_ref(3, ia) + deddt * dyiazia 
				+ d2eddt2 * ddtdyia * ddtdzia;
			hessz_ref(1, ia) = hessz_ref(1, ia) + deddt * dxiazia 
				+ d2eddt2 * ddtdzia * ddtdxia;
			hessz_ref(2, ia) = hessz_ref(2, ia) + deddt * dyiazia 
				+ d2eddt2 * ddtdzia * ddtdyia;
			hessz_ref(3, ia) = hessz_ref(3, ia) + deddt * dziazia 
				+ d2eddt2 * ddtdzia * ddtdzia;
			hessx_ref(1, ib) = hessx_ref(1, ib) + deddt * dxibxia 
				+ d2eddt2 * ddtdxia * ddtdxib;
			hessx_ref(2, ib) = hessx_ref(2, ib) + deddt * dyibxia 
				+ d2eddt2 * ddtdxia * ddtdyib;
			hessx_ref(3, ib) = hessx_ref(3, ib) + deddt * dzibxia 
				+ d2eddt2 * ddtdxia * ddtdzib;
			hessy_ref(1, ib) = hessy_ref(1, ib) + deddt * dxibyia 
				+ d2eddt2 * ddtdyia * ddtdxib;
			hessy_ref(2, ib) = hessy_ref(2, ib) + deddt * dyibyia 
				+ d2eddt2 * ddtdyia * ddtdyib;
			hessy_ref(3, ib) = hessy_ref(3, ib) + deddt * dzibyia 
				+ d2eddt2 * ddtdyia * ddtdzib;
			hessz_ref(1, ib) = hessz_ref(1, ib) + deddt * dxibzia 
				+ d2eddt2 * ddtdzia * ddtdxib;
			hessz_ref(2, ib) = hessz_ref(2, ib) + deddt * dyibzia 
				+ d2eddt2 * ddtdzia * ddtdyib;
			hessz_ref(3, ib) = hessz_ref(3, ib) + deddt * dzibzia 
				+ d2eddt2 * ddtdzia * ddtdzib;
			hessx_ref(1, ic) = hessx_ref(1, ic) + deddt * dxiaxic 
				+ d2eddt2 * ddtdxia * ddtdxic;
			hessx_ref(2, ic) = hessx_ref(2, ic) + deddt * dxiayic 
				+ d2eddt2 * ddtdxia * ddtdyic;
			hessx_ref(3, ic) = hessx_ref(3, ic) + deddt * dxiazic 
				+ d2eddt2 * ddtdxia * ddtdzic;
			hessy_ref(1, ic) = hessy_ref(1, ic) + deddt * dyiaxic 
				+ d2eddt2 * ddtdyia * ddtdxic;
			hessy_ref(2, ic) = hessy_ref(2, ic) + deddt * dyiayic 
				+ d2eddt2 * ddtdyia * ddtdyic;
			hessy_ref(3, ic) = hessy_ref(3, ic) + deddt * dyiazic 
				+ d2eddt2 * ddtdyia * ddtdzic;
			hessz_ref(1, ic) = hessz_ref(1, ic) + deddt * dziaxic 
				+ d2eddt2 * ddtdzia * ddtdxic;
			hessz_ref(2, ic) = hessz_ref(2, ic) + deddt * dziayic 
				+ d2eddt2 * ddtdzia * ddtdyic;
			hessz_ref(3, ic) = hessz_ref(3, ic) + deddt * dziazic 
				+ d2eddt2 * ddtdzia * ddtdzic;
		    } else if (ib == *iatom) {
			hessx_ref(1, ib) = hessx_ref(1, ib) + deddt * dxibxib 
				+ d2eddt2 * ddtdxib * ddtdxib;
			hessx_ref(2, ib) = hessx_ref(2, ib) + deddt * dxibyib 
				+ d2eddt2 * ddtdxib * ddtdyib;
			hessx_ref(3, ib) = hessx_ref(3, ib) + deddt * dxibzib 
				+ d2eddt2 * ddtdxib * ddtdzib;
			hessy_ref(1, ib) = hessy_ref(1, ib) + deddt * dxibyib 
				+ d2eddt2 * ddtdyib * ddtdxib;
			hessy_ref(2, ib) = hessy_ref(2, ib) + deddt * dyibyib 
				+ d2eddt2 * ddtdyib * ddtdyib;
			hessy_ref(3, ib) = hessy_ref(3, ib) + deddt * dyibzib 
				+ d2eddt2 * ddtdyib * ddtdzib;
			hessz_ref(1, ib) = hessz_ref(1, ib) + deddt * dxibzib 
				+ d2eddt2 * ddtdzib * ddtdxib;
			hessz_ref(2, ib) = hessz_ref(2, ib) + deddt * dyibzib 
				+ d2eddt2 * ddtdzib * ddtdyib;
			hessz_ref(3, ib) = hessz_ref(3, ib) + deddt * dzibzib 
				+ d2eddt2 * ddtdzib * ddtdzib;
			hessx_ref(1, ia) = hessx_ref(1, ia) + deddt * dxibxia 
				+ d2eddt2 * ddtdxib * ddtdxia;
			hessx_ref(2, ia) = hessx_ref(2, ia) + deddt * dxibyia 
				+ d2eddt2 * ddtdxib * ddtdyia;
			hessx_ref(3, ia) = hessx_ref(3, ia) + deddt * dxibzia 
				+ d2eddt2 * ddtdxib * ddtdzia;
			hessy_ref(1, ia) = hessy_ref(1, ia) + deddt * dyibxia 
				+ d2eddt2 * ddtdyib * ddtdxia;
			hessy_ref(2, ia) = hessy_ref(2, ia) + deddt * dyibyia 
				+ d2eddt2 * ddtdyib * ddtdyia;
			hessy_ref(3, ia) = hessy_ref(3, ia) + deddt * dyibzia 
				+ d2eddt2 * ddtdyib * ddtdzia;
			hessz_ref(1, ia) = hessz_ref(1, ia) + deddt * dzibxia 
				+ d2eddt2 * ddtdzib * ddtdxia;
			hessz_ref(2, ia) = hessz_ref(2, ia) + deddt * dzibyia 
				+ d2eddt2 * ddtdzib * ddtdyia;
			hessz_ref(3, ia) = hessz_ref(3, ia) + deddt * dzibzia 
				+ d2eddt2 * ddtdzib * ddtdzia;
			hessx_ref(1, ic) = hessx_ref(1, ic) + deddt * dxibxic 
				+ d2eddt2 * ddtdxib * ddtdxic;
			hessx_ref(2, ic) = hessx_ref(2, ic) + deddt * dxibyic 
				+ d2eddt2 * ddtdxib * ddtdyic;
			hessx_ref(3, ic) = hessx_ref(3, ic) + deddt * dxibzic 
				+ d2eddt2 * ddtdxib * ddtdzic;
			hessy_ref(1, ic) = hessy_ref(1, ic) + deddt * dyibxic 
				+ d2eddt2 * ddtdyib * ddtdxic;
			hessy_ref(2, ic) = hessy_ref(2, ic) + deddt * dyibyic 
				+ d2eddt2 * ddtdyib * ddtdyic;
			hessy_ref(3, ic) = hessy_ref(3, ic) + deddt * dyibzic 
				+ d2eddt2 * ddtdyib * ddtdzic;
			hessz_ref(1, ic) = hessz_ref(1, ic) + deddt * dzibxic 
				+ d2eddt2 * ddtdzib * ddtdxic;
			hessz_ref(2, ic) = hessz_ref(2, ic) + deddt * dzibyic 
				+ d2eddt2 * ddtdzib * ddtdyic;
			hessz_ref(3, ic) = hessz_ref(3, ic) + deddt * dzibzic 
				+ d2eddt2 * ddtdzib * ddtdzic;
		    } else if (ic == *iatom) {
			hessx_ref(1, ic) = hessx_ref(1, ic) + deddt * dxicxic 
				+ d2eddt2 * ddtdxic * ddtdxic;
			hessx_ref(2, ic) = hessx_ref(2, ic) + deddt * dxicyic 
				+ d2eddt2 * ddtdxic * ddtdyic;
			hessx_ref(3, ic) = hessx_ref(3, ic) + deddt * dxiczic 
				+ d2eddt2 * ddtdxic * ddtdzic;
			hessy_ref(1, ic) = hessy_ref(1, ic) + deddt * dxicyic 
				+ d2eddt2 * ddtdyic * ddtdxic;
			hessy_ref(2, ic) = hessy_ref(2, ic) + deddt * dyicyic 
				+ d2eddt2 * ddtdyic * ddtdyic;
			hessy_ref(3, ic) = hessy_ref(3, ic) + deddt * dyiczic 
				+ d2eddt2 * ddtdyic * ddtdzic;
			hessz_ref(1, ic) = hessz_ref(1, ic) + deddt * dxiczic 
				+ d2eddt2 * ddtdzic * ddtdxic;
			hessz_ref(2, ic) = hessz_ref(2, ic) + deddt * dyiczic 
				+ d2eddt2 * ddtdzic * ddtdyic;
			hessz_ref(3, ic) = hessz_ref(3, ic) + deddt * dziczic 
				+ d2eddt2 * ddtdzic * ddtdzic;
			hessx_ref(1, ib) = hessx_ref(1, ib) + deddt * dxibxic 
				+ d2eddt2 * ddtdxic * ddtdxib;
			hessx_ref(2, ib) = hessx_ref(2, ib) + deddt * dyibxic 
				+ d2eddt2 * ddtdxic * ddtdyib;
			hessx_ref(3, ib) = hessx_ref(3, ib) + deddt * dzibxic 
				+ d2eddt2 * ddtdxic * ddtdzib;
			hessy_ref(1, ib) = hessy_ref(1, ib) + deddt * dxibyic 
				+ d2eddt2 * ddtdyic * ddtdxib;
			hessy_ref(2, ib) = hessy_ref(2, ib) + deddt * dyibyic 
				+ d2eddt2 * ddtdyic * ddtdyib;
			hessy_ref(3, ib) = hessy_ref(3, ib) + deddt * dzibyic 
				+ d2eddt2 * ddtdyic * ddtdzib;
			hessz_ref(1, ib) = hessz_ref(1, ib) + deddt * dxibzic 
				+ d2eddt2 * ddtdzic * ddtdxib;
			hessz_ref(2, ib) = hessz_ref(2, ib) + deddt * dyibzic 
				+ d2eddt2 * ddtdzic * ddtdyib;
			hessz_ref(3, ib) = hessz_ref(3, ib) + deddt * dzibzic 
				+ d2eddt2 * ddtdzic * ddtdzib;
			hessx_ref(1, ia) = hessx_ref(1, ia) + deddt * dxiaxic 
				+ d2eddt2 * ddtdxic * ddtdxia;
			hessx_ref(2, ia) = hessx_ref(2, ia) + deddt * dyiaxic 
				+ d2eddt2 * ddtdxic * ddtdyia;
			hessx_ref(3, ia) = hessx_ref(3, ia) + deddt * dziaxic 
				+ d2eddt2 * ddtdxic * ddtdzia;
			hessy_ref(1, ia) = hessy_ref(1, ia) + deddt * dxiayic 
				+ d2eddt2 * ddtdyic * ddtdxia;
			hessy_ref(2, ia) = hessy_ref(2, ia) + deddt * dyiayic 
				+ d2eddt2 * ddtdyic * ddtdyia;
			hessy_ref(3, ia) = hessy_ref(3, ia) + deddt * dziayic 
				+ d2eddt2 * ddtdyic * ddtdzia;
			hessz_ref(1, ia) = hessz_ref(1, ia) + deddt * dxiazic 
				+ d2eddt2 * ddtdzic * ddtdxia;
			hessz_ref(2, ia) = hessz_ref(2, ia) + deddt * dyiazic 
				+ d2eddt2 * ddtdzic * ddtdyia;
			hessz_ref(3, ia) = hessz_ref(3, ia) + deddt * dziazic 
				+ d2eddt2 * ddtdzic * ddtdzia;
		    }

/*     construct a second orthogonal direction for linear angles */

		    if (linear) {
			linear = FALSE_;
			xpo = xp;
			ypo = yp;
			zpo = zp;
			xp = ypo * zab - zpo * yab;
			yp = zpo * xab - xpo * zab;
			zp = xpo * yab - ypo * xab;
			rp = sqrt(xp * xp + yp * yp + zp * zp);
			goto L10;
		    }
		}
	    }
	}
    }
    return 0;
} /* eangle2a_ */

#undef angtyp_ref
#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef iang_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine eangle2b  --  in-plane bend Hessian; numerical  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "eangle2b" computes projected in-plane bending first derivatives */
/*     for a single angle with respect to Cartesian coordinates; */
/*     used in computation of finite difference second derivatives */


/* Subroutine */ int eangle2b_(integer *i__)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal);

    /* Local variables */
    static integer ia, ib, ic, id;
    static doublereal dt, rm, xm, ym, zm, xt, yt, zt, dt2, dt3, dt4, rt2, xad,
	     yad, xia, yia, zia, xib, yib, dot, zib, xic, yic, zic, xid, yid, 
	    zid, zad, xbd, ybd, zbd, xcd, ycd, zcd, xip, yip, zip, xap, yap, 
	    zap, xcp, ycp, zcp, rap2, rcp2, term, ptrt2, ideal;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal deddt, angle, delta, force, terma, termc, delta2, 
	    dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, dedxic, dedyic, 
	    dedzic, dedxid, dedyid, dedzid, dedxip, cosine, dedyip, dedzip, 
	    dpdxia, dpdyia, dpdzia, dpdxic, dpdyic, dpdzic;


#define dea_ref(a_1,a_2) deriv_1.dea[(a_2)*3 + a_1 - 4]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]



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
/*     ##  angle.i  --  bond angles within the current structure  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     ak       harmonic angle force constant (kcal/mole/rad**2) */
/*     anat     ideal bond angle or phase shift angle (degrees) */
/*     afld     periodicity for Fourier bond angle term */
/*     nangle   total number of bond angles in the system */
/*     iang     numbers of the atoms in each bond angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  angpot.i  --  specifics of angle bend functional forms  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     angunit    convert angle bending energy to kcal/mole */
/*     stbnunit   convert stretch-bend energy to kcal/mole */
/*     aaunit     convert angle-angle energy to kcal/mole */
/*     opbunit    convert out-of-plane bend energy to kcal/mole */
/*     opdunit    convert out-of-plane distance energy to kcal/mole */
/*     cang       cubic coefficient in angle bending potential */
/*     qang       quartic coefficient in angle bending potential */
/*     pang       quintic coefficient in angle bending potential */
/*     sang       sextic coefficient in angle bending potential */
/*     copb       cubic coefficient in out-of-plane bend potential */
/*     qopb       quartic coefficient in out-of-plane bend potential */
/*     popb       quintic coefficient in out-of-plane bend potential */
/*     sopb       sextic coefficient in out-of-plane bend potential */
/*     copd       cubic coefficient in out-of-plane distance potential */
/*     qopd       quartic coefficient in out-of-plane distance potential */
/*     popd       quintic coefficient in out-of-plane distance potential */
/*     sopd       sextic coefficient in out-of-plane distance potential */
/*     angtyp     type of angle bending function for each bond angle */
/*     opbtyp     type of out-of-plane bend potential energy function */




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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  deriv.i  --  Cartesian coordinate derivative components  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     desum   total energy Cartesian coordinate derivatives */
/*     deb     bond stretch Cartesian coordinate derivatives */
/*     dea     angle bend Cartesian coordinate derivatives */
/*     deba    stretch-bend Cartesian coordinate derivatives */
/*     deub    Urey-Bradley Cartesian coordinate derivatives */
/*     deaa    angle-angle Cartesian coordinate derivatives */
/*     deopb   out-of-plane bend Cartesian coordinate derivatives */
/*     deopd   out-of-plane distance Cartesian coordinate derivatives */
/*     deid    improper dihedral Cartesian coordinate derivatives */
/*     deit    improper torsion Cartesian coordinate derivatives */
/*     det     torsional Cartesian coordinate derivatives */
/*     dept    pi-orbital torsion Cartesian coordinate derivatives */
/*     debt    stretch-torsion Cartesian coordinate derivatives */
/*     dett    torsion-torsion Cartesian coordinate derivatives */
/*     dev     van der Waals Cartesian coordinate derivatives */
/*     dec     charge-charge Cartesian coordinate derivatives */
/*     decd    charge-dipole Cartesian coordinate derivatives */
/*     ded     dipole-dipole Cartesian coordinate derivatives */
/*     dem     multipole Cartesian coordinate derivatives */
/*     dep     polarization Cartesian coordinate derivatives */
/*     der     reaction field Cartesian coordinate derivatives */
/*     des     solvation Cartesian coordinate derivatives */
/*     delf    metal ligand field Cartesian coordinate derivatives */
/*     deg     geometric restraint Cartesian coordinate derivatives */
/*     dex     extra energy term Cartesian coordinate derivatives */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  math.i  --  mathematical and geometrical constants  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     radian   conversion factor from radians to degrees */
/*     pi       numerical value of the geometric constant */
/*     sqrtpi   numerical value of the square root of Pi */
/*     logten   numerical value of the natural log of ten */
/*     sqrttwo  numerical value of the square root of two */
/*     twosix   numerical value of the sixth root of two */




/*     set the atom numbers and parameters for this angle */

    ia = iang_ref(1, *i__);
    ib = iang_ref(2, *i__);
    ic = iang_ref(3, *i__);
    id = iang_ref(4, *i__);
    ideal = angle_1.anat[*i__ - 1];
    force = angle_1.ak[*i__ - 1];

/*     get the coordinates of the atoms in the angle */

    xia = atoms_1.x[ia - 1];
    yia = atoms_1.y[ia - 1];
    zia = atoms_1.z__[ia - 1];
    xib = atoms_1.x[ib - 1];
    yib = atoms_1.y[ib - 1];
    zib = atoms_1.z__[ib - 1];
    xic = atoms_1.x[ic - 1];
    yic = atoms_1.y[ic - 1];
    zic = atoms_1.z__[ic - 1];
    xid = atoms_1.x[id - 1];
    yid = atoms_1.y[id - 1];
    zid = atoms_1.z__[id - 1];

/*     zero out the first derivative components */

    dea_ref(1, ia) = 0.;
    dea_ref(2, ia) = 0.;
    dea_ref(3, ia) = 0.;
    dea_ref(1, ib) = 0.;
    dea_ref(2, ib) = 0.;
    dea_ref(3, ib) = 0.;
    dea_ref(1, ic) = 0.;
    dea_ref(2, ic) = 0.;
    dea_ref(3, ic) = 0.;
    dea_ref(1, id) = 0.;
    dea_ref(2, id) = 0.;
    dea_ref(3, id) = 0.;

/*     compute the projected in-plane angle gradient */

    xad = xia - xid;
    yad = yia - yid;
    zad = zia - zid;
    xbd = xib - xid;
    ybd = yib - yid;
    zbd = zib - zid;
    xcd = xic - xid;
    ycd = yic - yid;
    zcd = zic - zid;
    if (bound_1.use_polymer__) {
	image_(&xad, &yad, &zad);
	image_(&xbd, &ybd, &zbd);
	image_(&xcd, &ycd, &zcd);
    }
    xt = yad * zcd - zad * ycd;
    yt = zad * xcd - xad * zcd;
    zt = xad * ycd - yad * xcd;
    rt2 = xt * xt + yt * yt + zt * zt;
    delta = -(xt * xbd + yt * ybd + zt * zbd) / rt2;
    xip = xib + xt * delta;
    yip = yib + yt * delta;
    zip = zib + zt * delta;
    xap = xia - xip;
    yap = yia - yip;
    zap = zia - zip;
    xcp = xic - xip;
    ycp = yic - yip;
    zcp = zic - zip;
    if (bound_1.use_polymer__) {
	image_(&xap, &yap, &zap);
	image_(&xcp, &ycp, &zcp);
    }
    rap2 = xap * xap + yap * yap + zap * zap;
    rcp2 = xcp * xcp + ycp * ycp + zcp * zcp;
    if (rap2 != 0. && rcp2 != 0.) {
	xm = ycp * zap - zcp * yap;
	ym = zcp * xap - xcp * zap;
	zm = xcp * yap - ycp * xap;
	rm = sqrt(xm * xm + ym * ym + zm * zm);
	rm = max(rm,1e-6);
	dot = xap * xcp + yap * ycp + zap * zcp;
	cosine = dot / sqrt(rap2 * rcp2);
/* Computing MIN */
	d__1 = 1., d__2 = max(-1.,cosine);
	cosine = min(d__1,d__2);
	angle = acos(cosine) * 57.29577951308232088;

/*     get the master chain rule term for derivatives */

	dt = angle - ideal;
	dt2 = dt * dt;
	dt3 = dt2 * dt;
	dt4 = dt2 * dt2;
	deddt = angpot_1.angunit * force * dt * 57.29577951308232088 * (
		angpot_1.cang * 3. * dt + 2. + angpot_1.qang * 4. * dt2 + 
		angpot_1.pang * 5. * dt3 + angpot_1.sang * 6. * dt4);

/*     chain rule terms for first derivative components */

	terma = -deddt / (rap2 * rm);
	termc = deddt / (rcp2 * rm);
	dedxia = terma * (yap * zm - zap * ym);
	dedyia = terma * (zap * xm - xap * zm);
	dedzia = terma * (xap * ym - yap * xm);
	dedxic = termc * (ycp * zm - zcp * ym);
	dedyic = termc * (zcp * xm - xcp * zm);
	dedzic = termc * (xcp * ym - ycp * xm);
	dedxip = -dedxia - dedxic;
	dedyip = -dedyia - dedyic;
	dedzip = -dedzia - dedzic;

/*     chain rule components for the projection of the central atom */

	delta2 = delta * 2.;
	ptrt2 = (dedxip * xt + dedyip * yt + dedzip * zt) / rt2;
	term = zcd * ybd - ycd * zbd + delta2 * (yt * zcd - zt * ycd);
	dpdxia = delta * (ycd * dedzip - zcd * dedyip) + term * ptrt2;
	term = xcd * zbd - zcd * xbd + delta2 * (zt * xcd - xt * zcd);
	dpdyia = delta * (zcd * dedxip - xcd * dedzip) + term * ptrt2;
	term = ycd * xbd - xcd * ybd + delta2 * (xt * ycd - yt * xcd);
	dpdzia = delta * (xcd * dedyip - ycd * dedxip) + term * ptrt2;
	term = yad * zbd - zad * ybd + delta2 * (zt * yad - yt * zad);
	dpdxic = delta * (zad * dedyip - yad * dedzip) + term * ptrt2;
	term = zad * xbd - xad * zbd + delta2 * (xt * zad - zt * xad);
	dpdyic = delta * (xad * dedzip - zad * dedxip) + term * ptrt2;
	term = xad * ybd - yad * xbd + delta2 * (yt * xad - xt * yad);
	dpdzic = delta * (yad * dedxip - xad * dedyip) + term * ptrt2;

/*     compute derivative components for this interaction */

	dedxia += dpdxia;
	dedyia += dpdyia;
	dedzia += dpdzia;
	dedxib = dedxip;
	dedyib = dedyip;
	dedzib = dedzip;
	dedxic += dpdxic;
	dedyic += dpdyic;
	dedzic += dpdzic;
	dedxid = -dedxia - dedxib - dedxic;
	dedyid = -dedyia - dedyib - dedyic;
	dedzid = -dedzia - dedzib - dedzic;

/*     set the in-plane angle bending first derivatives */

	dea_ref(1, ia) = dedxia;
	dea_ref(2, ia) = dedyia;
	dea_ref(3, ia) = dedzia;
	dea_ref(1, ib) = dedxib;
	dea_ref(2, ib) = dedyib;
	dea_ref(3, ib) = dedzib;
	dea_ref(1, ic) = dedxic;
	dea_ref(2, ic) = dedyic;
	dea_ref(3, ic) = dedzic;
	dea_ref(1, id) = dedxid;
	dea_ref(2, id) = dedyid;
	dea_ref(3, id) = dedzid;
    }
    return 0;
} /* eangle2b_ */

#undef iang_ref
#undef dea_ref


