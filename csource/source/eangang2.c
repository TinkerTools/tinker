/* eangang2.f -- translated by f2c (version 20050501).
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
    doublereal kaa[100000];
    integer nangang, iaa[200000]	/* was [2][100000] */;
} angang_;

#define angang_1 angang_

struct {
    doublereal ak[75000], anat[75000], afld[75000];
    integer nangle, iang[300000]	/* was [4][75000] */;
} angle_;

#define angle_1 angle_

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
    doublereal angunit, stbnunit, aaunit, opbunit, opdunit, cang, qang, pang, 
	    sang, copb, qopb, popb, sopb, copd, qopd, popd, sopd;
    char angtyp[600000], opbtyp[8];
} angpot_;

#define angpot_1 angpot_

struct {
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine eangang2  --  angle-angle Hessian; numerical  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "eangang2" calculates the angle-angle potential energy */
/*     second derivatives with respect to Cartesian coordinates */
/*     using finite difference methods */


/* Subroutine */ int eangang2_(integer *i__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static logical twosided;
    extern /* Subroutine */ int eangang2a_(integer *);
    static integer j, k;
    static doublereal d0[75000]	/* was [3][25000] */;
    static integer ia, ib, ic, id, ie;
    static doublereal old, eps, fgrp, term;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static integer iangang;
    static logical proceed;


#define d0_ref(a_1,a_2) d0[(a_2)*3 + a_1 - 4]
#define iaa_ref(a_1,a_2) angang_1.iaa[(a_2)*2 + a_1 - 3]
#define deaa_ref(a_1,a_2) deriv_1.deaa[(a_2)*3 + a_1 - 4]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define hessx_ref(a_1,a_2) hessn_1.hessx[(a_2)*3 + a_1 - 4]
#define hessy_ref(a_1,a_2) hessn_1.hessy[(a_2)*3 + a_1 - 4]
#define hessz_ref(a_1,a_2) hessn_1.hessz[(a_2)*3 + a_1 - 4]



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
/*     ##  angang.i  --  angle-angle terms in current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     kaa       force constant for angle-angle cross terms */
/*     nangang   total number of angle-angle interactions */
/*     iaa       angle numbers used in each angle-angle term */




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




/*     set stepsize for derivatives and default group weight */

    eps = 1e-5;
    fgrp = 1.;
    twosided = FALSE_;
    if (atoms_1.n <= 50) {
	twosided = TRUE_;
    }

/*     compute numerical angle-angle Hessian for current atom */

    i__1 = angang_1.nangang;
    for (iangang = 1; iangang <= i__1; ++iangang) {
	j = iaa_ref(1, iangang);
	k = iaa_ref(2, iangang);
	ia = iang_ref(1, j);
	ib = iang_ref(2, j);
	ic = iang_ref(3, j);
	id = iang_ref(1, k);
	ie = iang_ref(3, k);

/*     decide whether to compute the current interaction */

	proceed = *i__ == ia || *i__ == ib || *i__ == ic || *i__ == id || *
		i__ == ie;
	if (proceed && group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &ie, &c__0);
	}

/*     eliminate any duplicate atoms in the pair of angles */

	if (proceed) {
	    if (id == ia || id == ic) {
		id = ie;
		ie = 0;
	    } else if (ie == ia || ie == ic) {
		ie = 0;
	    }
	    term = fgrp / eps;

/*     find first derivatives for the base structure */

	    if (! twosided) {
		eangang2a_(&iangang);
		for (j = 1; j <= 3; ++j) {
		    d0_ref(j, ia) = deaa_ref(j, ia);
		    d0_ref(j, ib) = deaa_ref(j, ib);
		    d0_ref(j, ic) = deaa_ref(j, ic);
		    d0_ref(j, id) = deaa_ref(j, id);
		    if (ie != 0) {
			d0_ref(j, ie) = deaa_ref(j, ie);
		    }
		}
	    }

/*     find numerical x-components via perturbed structures */

	    old = atoms_1.x[*i__ - 1];
	    if (twosided) {
		atoms_1.x[*i__ - 1] -= eps * .5;
		eangang2a_(&iangang);
		for (j = 1; j <= 3; ++j) {
		    d0_ref(j, ia) = deaa_ref(j, ia);
		    d0_ref(j, ib) = deaa_ref(j, ib);
		    d0_ref(j, ic) = deaa_ref(j, ic);
		    d0_ref(j, id) = deaa_ref(j, id);
		    if (ie != 0) {
			d0_ref(j, ie) = deaa_ref(j, ie);
		    }
		}
	    }
	    atoms_1.x[*i__ - 1] += eps;
	    eangang2a_(&iangang);
	    atoms_1.x[*i__ - 1] = old;
	    for (j = 1; j <= 3; ++j) {
		hessx_ref(j, ia) = hessx_ref(j, ia) + term * (deaa_ref(j, ia) 
			- d0_ref(j, ia));
		hessx_ref(j, ib) = hessx_ref(j, ib) + term * (deaa_ref(j, ib) 
			- d0_ref(j, ib));
		hessx_ref(j, ic) = hessx_ref(j, ic) + term * (deaa_ref(j, ic) 
			- d0_ref(j, ic));
		hessx_ref(j, id) = hessx_ref(j, id) + term * (deaa_ref(j, id) 
			- d0_ref(j, id));
		if (ie != 0) {
		    hessx_ref(j, ie) = hessx_ref(j, ie) + term * (deaa_ref(j, 
			    ie) - d0_ref(j, ie));
		}
	    }

/*     find numerical y-components via perturbed structures */

	    old = atoms_1.y[*i__ - 1];
	    if (twosided) {
		atoms_1.y[*i__ - 1] -= eps * .5;
		eangang2a_(&iangang);
		for (j = 1; j <= 3; ++j) {
		    d0_ref(j, ia) = deaa_ref(j, ia);
		    d0_ref(j, ib) = deaa_ref(j, ib);
		    d0_ref(j, ic) = deaa_ref(j, ic);
		    d0_ref(j, id) = deaa_ref(j, id);
		    if (ie != 0) {
			d0_ref(j, ie) = deaa_ref(j, ie);
		    }
		}
	    }
	    atoms_1.y[*i__ - 1] += eps;
	    eangang2a_(&iangang);
	    atoms_1.y[*i__ - 1] = old;
	    for (j = 1; j <= 3; ++j) {
		hessy_ref(j, ia) = hessy_ref(j, ia) + term * (deaa_ref(j, ia) 
			- d0_ref(j, ia));
		hessy_ref(j, ib) = hessy_ref(j, ib) + term * (deaa_ref(j, ib) 
			- d0_ref(j, ib));
		hessy_ref(j, ic) = hessy_ref(j, ic) + term * (deaa_ref(j, ic) 
			- d0_ref(j, ic));
		hessy_ref(j, id) = hessy_ref(j, id) + term * (deaa_ref(j, id) 
			- d0_ref(j, id));
		if (ie != 0) {
		    hessy_ref(j, ie) = hessy_ref(j, ie) + term * (deaa_ref(j, 
			    ie) - d0_ref(j, ie));
		}
	    }

/*     find numerical z-components via perturbed structures */

	    old = atoms_1.z__[*i__ - 1];
	    if (twosided) {
		atoms_1.z__[*i__ - 1] -= eps * .5;
		eangang2a_(&iangang);
		for (j = 1; j <= 3; ++j) {
		    d0_ref(j, ia) = deaa_ref(j, ia);
		    d0_ref(j, ib) = deaa_ref(j, ib);
		    d0_ref(j, ic) = deaa_ref(j, ic);
		    d0_ref(j, id) = deaa_ref(j, id);
		    if (ie != 0) {
			d0_ref(j, ie) = deaa_ref(j, ie);
		    }
		}
	    }
	    atoms_1.z__[*i__ - 1] += eps;
	    eangang2a_(&iangang);
	    atoms_1.z__[*i__ - 1] = old;
	    for (j = 1; j <= 3; ++j) {
		hessz_ref(j, ia) = hessz_ref(j, ia) + term * (deaa_ref(j, ia) 
			- d0_ref(j, ia));
		hessz_ref(j, ib) = hessz_ref(j, ib) + term * (deaa_ref(j, ib) 
			- d0_ref(j, ib));
		hessz_ref(j, ic) = hessz_ref(j, ic) + term * (deaa_ref(j, ic) 
			- d0_ref(j, ic));
		hessz_ref(j, id) = hessz_ref(j, id) + term * (deaa_ref(j, id) 
			- d0_ref(j, id));
		if (ie != 0) {
		    hessz_ref(j, ie) = hessz_ref(j, ie) + term * (deaa_ref(j, 
			    ie) - d0_ref(j, ie));
		}
	    }
	}
    }
    return 0;
} /* eangang2_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef iang_ref
#undef deaa_ref
#undef iaa_ref
#undef d0_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine eangang2a  --  angle-angle interaction derivs  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "eangang2a" calculates the angle-angle first derivatives for */
/*     a single interaction with respect to Cartesian coordinates; */
/*     used in computation of finite difference second derivatives */


/* Subroutine */ int eangang2a_(integer *i__)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal);

    /* Local variables */
    static integer j, k, ia, ib, ic, id, ie;
    static doublereal rp, rq, xp, yp, zp, xq, yq, zq, dt1, dt2, xab, yab, xia,
	     yia, zia, xib, yib, dot, zib, xic, yic, zic, xid, yid, zid, xie, 
	    yie, zie, zab, xcb, ycb, zcb, xdb, ydb, zdb, xeb, yeb, zeb, rab2, 
	    rcb2, rdb2, reb2, term;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal angle, terma, termc, termd, terme, deddt1, deddt2, 
	    dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, dedxic, dedyic, 
	    dedzic, dedxid, dedyid, dedzid, dedxie, cosine, dedyie, dedzie;


#define iaa_ref(a_1,a_2) angang_1.iaa[(a_2)*2 + a_1 - 3]
#define deaa_ref(a_1,a_2) deriv_1.deaa[(a_2)*3 + a_1 - 4]
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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  angang.i  --  angle-angle terms in current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     kaa       force constant for angle-angle cross terms */
/*     nangang   total number of angle-angle interactions */
/*     iaa       angle numbers used in each angle-angle term */




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




/*     set the coordinates of the involved atoms */

    j = iaa_ref(1, *i__);
    k = iaa_ref(2, *i__);
    ia = iang_ref(1, j);
    ib = iang_ref(2, j);
    ic = iang_ref(3, j);
    id = iang_ref(1, k);
    ie = iang_ref(3, k);
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
    xie = atoms_1.x[ie - 1];
    yie = atoms_1.y[ie - 1];
    zie = atoms_1.z__[ie - 1];

/*     zero out the first derivative components */

    deaa_ref(1, ia) = 0.;
    deaa_ref(2, ia) = 0.;
    deaa_ref(3, ia) = 0.;
    deaa_ref(1, ib) = 0.;
    deaa_ref(2, ib) = 0.;
    deaa_ref(3, ib) = 0.;
    deaa_ref(1, ic) = 0.;
    deaa_ref(2, ic) = 0.;
    deaa_ref(3, ic) = 0.;
    deaa_ref(1, id) = 0.;
    deaa_ref(2, id) = 0.;
    deaa_ref(3, id) = 0.;
    deaa_ref(1, ie) = 0.;
    deaa_ref(2, ie) = 0.;
    deaa_ref(3, ie) = 0.;

/*     compute the values of the two bond angles */

    xab = xia - xib;
    yab = yia - yib;
    zab = zia - zib;
    xcb = xic - xib;
    ycb = yic - yib;
    zcb = zic - zib;
    xdb = xid - xib;
    ydb = yid - yib;
    zdb = zid - zib;
    xeb = xie - xib;
    yeb = yie - yib;
    zeb = zie - zib;
    if (bound_1.use_polymer__) {
	image_(&xab, &yab, &zab);
	image_(&xcb, &ycb, &zcb);
	image_(&xdb, &ydb, &zdb);
	image_(&xeb, &yeb, &zeb);
    }
    rab2 = xab * xab + yab * yab + zab * zab;
    rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
    rdb2 = xdb * xdb + ydb * ydb + zdb * zdb;
    reb2 = xeb * xeb + yeb * yeb + zeb * zeb;
    xp = ycb * zab - zcb * yab;
    yp = zcb * xab - xcb * zab;
    zp = xcb * yab - ycb * xab;
    xq = yeb * zdb - zeb * ydb;
    yq = zeb * xdb - xeb * zdb;
    zq = xeb * ydb - yeb * xdb;
    rp = sqrt(xp * xp + yp * yp + zp * zp);
    rq = sqrt(xq * xq + yq * yq + zq * zq);
    if (rp * rq != 0.) {
	dot = xab * xcb + yab * ycb + zab * zcb;
	cosine = dot / sqrt(rab2 * rcb2);
/* Computing MIN */
	d__1 = 1., d__2 = max(-1.,cosine);
	cosine = min(d__1,d__2);
	angle = acos(cosine) * 57.29577951308232088;
	dt1 = angle - angle_1.anat[j - 1];
	dot = xdb * xeb + ydb * yeb + zdb * zeb;
	cosine = dot / sqrt(rdb2 * reb2);
/* Computing MIN */
	d__1 = 1., d__2 = max(-1.,cosine);
	cosine = min(d__1,d__2);
	angle = acos(cosine) * 57.29577951308232088;
	dt2 = angle - angle_1.anat[k - 1];

/*     get the energy and master chain rule terms for derivatives */

	term = angpot_1.aaunit * 57.29577951308232088 * angang_1.kaa[*i__ - 1]
		;
	deddt1 = term * dt2;
	deddt2 = term * dt1;

/*     find chain rule terms for the first bond angle deviation */

	terma = -deddt1 / (rab2 * rp);
	termc = deddt1 / (rcb2 * rp);
	dedxia = terma * (yab * zp - zab * yp);
	dedyia = terma * (zab * xp - xab * zp);
	dedzia = terma * (xab * yp - yab * xp);
	dedxic = termc * (ycb * zp - zcb * yp);
	dedyic = termc * (zcb * xp - xcb * zp);
	dedzic = termc * (xcb * yp - ycb * xp);

/*     find chain rule terms for the second bond angle deviation */

	termd = -deddt2 / (rdb2 * rq);
	terme = deddt2 / (reb2 * rq);
	dedxid = termd * (ydb * zq - zdb * yq);
	dedyid = termd * (zdb * xq - xdb * zq);
	dedzid = termd * (xdb * yq - ydb * xq);
	dedxie = terme * (yeb * zq - zeb * yq);
	dedyie = terme * (zeb * xq - xeb * zq);
	dedzie = terme * (xeb * yq - yeb * xq);

/*     get the central atom derivative terms by difference */

	dedxib = -dedxia - dedxic - dedxid - dedxie;
	dedyib = -dedyia - dedyic - dedyid - dedyie;
	dedzib = -dedzia - dedzic - dedzid - dedzie;

/*     set the angle-angle interaction first derivatives */

	deaa_ref(1, ia) = deaa_ref(1, ia) + dedxia;
	deaa_ref(2, ia) = deaa_ref(2, ia) + dedyia;
	deaa_ref(3, ia) = deaa_ref(3, ia) + dedzia;
	deaa_ref(1, ib) = deaa_ref(1, ib) + dedxib;
	deaa_ref(2, ib) = deaa_ref(2, ib) + dedyib;
	deaa_ref(3, ib) = deaa_ref(3, ib) + dedzib;
	deaa_ref(1, ic) = deaa_ref(1, ic) + dedxic;
	deaa_ref(2, ic) = deaa_ref(2, ic) + dedyic;
	deaa_ref(3, ic) = deaa_ref(3, ic) + dedzic;
	deaa_ref(1, id) = deaa_ref(1, id) + dedxid;
	deaa_ref(2, id) = deaa_ref(2, id) + dedyid;
	deaa_ref(3, id) = deaa_ref(3, id) + dedzid;
	deaa_ref(1, ie) = deaa_ref(1, ie) + dedxie;
	deaa_ref(2, ie) = deaa_ref(2, ie) + dedyie;
	deaa_ref(3, ie) = deaa_ref(3, ie) + dedzie;
    }
    return 0;
} /* eangang2a_ */

#undef iang_ref
#undef deaa_ref
#undef iaa_ref


