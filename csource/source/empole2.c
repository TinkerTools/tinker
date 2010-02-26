/* empole2.f -- translated by f2c (version 20050501).
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
    doublereal hessx[75000]	/* was [3][25000] */, hessy[75000]	/* 
	    was [3][25000] */, hessz[75000]	/* was [3][25000] */;
} hessn_;

#define hessn_1 hessn_

struct {
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

struct {
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

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
    doublereal electric, dielec, ebuffer, c2scale, c3scale, c4scale, c5scale;
    logical neutnbr, neutcut;
} chgpot_;

#define chgpot_1 chgpot_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal vdwcut, chgcut, dplcut, mpolecut, vdwtaper, chgtaper, dpltaper,
	     mpoletaper, ewaldcut;
    logical use_ewald__, use_lights__, use_list__, use_vlist__, use_clist__, 
	    use_mlist__;
} cutoff_;

#define cutoff_1 cutoff_

struct {
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

struct {
    doublereal m2scale, m3scale, m4scale, m5scale;
} mplpot_;

#define mplpot_1 mplpot_

struct {
    doublereal polarity[25000], thole[25000], pdamp[25000], uind[75000]	/* 
	    was [3][25000] */, uinp[75000]	/* was [3][25000] */, uinds[
	    75000]	/* was [3][25000] */, uinps[75000]	/* was [3][
	    25000] */;
    integer npolar;
} polar_;

#define polar_1 polar_

struct {
    integer np11[25000], ip11[2500000]	/* was [100][25000] */, np12[25000], 
	    ip12[1250000]	/* was [50][25000] */, np13[25000], ip13[
	    1250000]	/* was [50][25000] */, np14[25000], ip14[1250000]	
	    /* was [50][25000] */;
} polgrp_;

#define polgrp_1 polgrp_

struct {
    doublereal poleps, polsor, p2scale, p3scale, p4scale, p5scale, d1scale, 
	    d2scale, d3scale, d4scale, u1scale, u2scale, u3scale, u4scale;
    char poltyp[6];
} polpot_;

#define polpot_1 polpot_

struct {
    doublereal off, off2, cut, cut2, c0, c1, c2, c3, c4, c5, f0, f1, f2, f3, 
	    f4, f5, f6, f7;
} shunt_;

#define shunt_1 shunt_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__0 = 0;



/*     ############################################################# */
/*     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine empole2  --  multipole & polarization Hessian  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "empole2" calculates second derivatives of the multipole and */
/*     dipole polarization energy for a single atom at a time */

/*     since polarization effects are many-body, it is incorrect to */
/*     neglect interactions with atoms not directly involved as the */
/*     multipole site or its local axis definition; to get better */
/*     accuracy, "list" should include all atoms, an option available */
/*     by setting "biglist" to "true" */

/*     also, the "reinduce" flag controls whether the induced dipoles */
/*     are recomputed every time an atom is moved during computation */
/*     of the numerical Hessian; setting the flag to "true" produces */
/*     a much slower calculation, but can greatly aid convergence of */
/*     minimizations, accuracy of vibrational frequencies, etc. */

/*     the "twosided" flag controls use of one-sided vs. two-sided */
/*     numerical derivatives; setting the flag to "true" gives a more */
/*     accurate Hessian at the expense of 50% additional computation */

/*     in the current version, all of the above accuracy improvements */
/*     are turned on for systems containing 50 atoms or fewer */


/* Subroutine */ int empole2_(integer *i__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static logical reinduce, twosided;
    static integer j, k;
    static doublereal d0[75000]	/* was [3][25000] */, old, eps;
    static integer list[25000], nlist;
    static logical biglist;
    extern /* Subroutine */ int empole2a_(integer *, integer *, logical *);


#define d0_ref(a_1,a_2) d0[(a_2)*3 + a_1 - 4]
#define dem_ref(a_1,a_2) deriv_1.dem[(a_2)*3 + a_1 - 4]
#define dep_ref(a_1,a_2) deriv_1.dep[(a_2)*3 + a_1 - 4]
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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  potent.i  --  usage of each potential energy component  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     use_bond    logical flag governing use of bond stretch potential */
/*     use_angle   logical flag governing use of angle bend potential */
/*     use_strbnd  logical flag governing use of stretch-bend potential */
/*     use_urey    logical flag governing use of Urey-Bradley potential */
/*     use_angang  logical flag governing use of angle-angle cross term */
/*     use_opbend  logical flag governing use of out-of-plane bend term */
/*     use_opdist  logical flag governing use of out-of-plane distance */
/*     use_improp  logical flag governing use of improper dihedral term */
/*     use_imptor  logical flag governing use of improper torsion term */
/*     use_tors    logical flag governing use of torsional potential */
/*     use_pitors  logical flag governing use of pi-orbital torsion term */
/*     use_strtor  logical flag governing use of stretch-torsion term */
/*     use_tortor  logical flag governing use of torsion-torsion term */
/*     use_vdw     logical flag governing use of vdw der Waals potential */
/*     use_charge  logical flag governing use of charge-charge potential */
/*     use_chgdpl  logical flag governing use of charge-dipole potential */
/*     use_dipole  logical flag governing use of dipole-dipole potential */
/*     use_mpole   logical flag governing use of multipole potential */
/*     use_polar   logical flag governing use of polarization term */
/*     use_rxnfld  logical flag governing use of reaction field term */
/*     use_solv    logical flag governing use of continuum solvation */
/*     use_metal   logical flag governing use of ligand field term */
/*     use_geom    logical flag governing use of geometric restraints */
/*     use_extra   logical flag governing use of extra potential term */
/*     use_born    logical flag governing use of Born radii values */
/*     use_orbit   logical flag governing use of pisystem computation */




/*     set the default stepsize and flag for induced dipoles */

    eps = 1e-5;
    biglist = FALSE_;
    reinduce = FALSE_;
    twosided = FALSE_;
    if (atoms_1.n <= 50) {
	biglist = TRUE_;
	reinduce = TRUE_;
	twosided = TRUE_;
    }

/*     find the multipole definitions involving the current atom; */
/*     results in a faster but approximate Hessian calculation */

    nlist = 0;
    i__1 = mpole_1.npole;
    for (k = 1; k <= i__1; ++k) {
	if (biglist || mpole_1.ipole[k - 1] == *i__ || mpole_1.zaxis[k - 1] ==
		 *i__ || mpole_1.xaxis[k - 1] == *i__) {
	    ++nlist;
	    list[nlist - 1] = k;
	}
    }

/*     get multipole first derivatives for the base structure */

    if (! twosided) {
	empole2a_(&nlist, list, &reinduce);
	i__1 = atoms_1.n;
	for (k = 1; k <= i__1; ++k) {
	    for (j = 1; j <= 3; ++j) {
		d0_ref(j, k) = dem_ref(j, k) + dep_ref(j, k);
	    }
	}
    }

/*     find numerical x-components via perturbed structures */

    old = atoms_1.x[*i__ - 1];
    if (twosided) {
	atoms_1.x[*i__ - 1] -= eps * .5;
	empole2a_(&nlist, list, &reinduce);
	i__1 = atoms_1.n;
	for (k = 1; k <= i__1; ++k) {
	    for (j = 1; j <= 3; ++j) {
		d0_ref(j, k) = dem_ref(j, k) + dep_ref(j, k);
	    }
	}
    }
    atoms_1.x[*i__ - 1] += eps;
    empole2a_(&nlist, list, &reinduce);
    atoms_1.x[*i__ - 1] = old;
    i__1 = atoms_1.n;
    for (k = 1; k <= i__1; ++k) {
	for (j = 1; j <= 3; ++j) {
	    hessx_ref(j, k) = hessx_ref(j, k) + (dem_ref(j, k) + dep_ref(j, k)
		     - d0_ref(j, k)) / eps;
	}
    }

/*     find numerical y-components via perturbed structures */

    old = atoms_1.y[*i__ - 1];
    if (twosided) {
	atoms_1.y[*i__ - 1] -= eps * .5;
	empole2a_(&nlist, list, &reinduce);
	i__1 = atoms_1.n;
	for (k = 1; k <= i__1; ++k) {
	    for (j = 1; j <= 3; ++j) {
		d0_ref(j, k) = dem_ref(j, k) + dep_ref(j, k);
	    }
	}
    }
    atoms_1.y[*i__ - 1] += eps;
    empole2a_(&nlist, list, &reinduce);
    atoms_1.y[*i__ - 1] = old;
    i__1 = atoms_1.n;
    for (k = 1; k <= i__1; ++k) {
	for (j = 1; j <= 3; ++j) {
	    hessy_ref(j, k) = hessy_ref(j, k) + (dem_ref(j, k) + dep_ref(j, k)
		     - d0_ref(j, k)) / eps;
	}
    }

/*     find numerical z-components via perturbed structures */

    old = atoms_1.z__[*i__ - 1];
    if (twosided) {
	atoms_1.z__[*i__ - 1] -= eps * .5;
	empole2a_(&nlist, list, &reinduce);
	i__1 = atoms_1.n;
	for (k = 1; k <= i__1; ++k) {
	    for (j = 1; j <= 3; ++j) {
		d0_ref(j, k) = dem_ref(j, k) + dep_ref(j, k);
	    }
	}
    }
    atoms_1.z__[*i__ - 1] += eps;
    empole2a_(&nlist, list, &reinduce);
    atoms_1.z__[*i__ - 1] = old;
    i__1 = atoms_1.n;
    for (k = 1; k <= i__1; ++k) {
	for (j = 1; j <= 3; ++j) {
	    hessz_ref(j, k) = hessz_ref(j, k) + (dem_ref(j, k) + dep_ref(j, k)
		     - d0_ref(j, k)) / eps;
	}
    }
    return 0;
} /* empole2_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef dep_ref
#undef dem_ref
#undef d0_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine empole2a  --  mpole & polar Hessian; numerical  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "empole2a" computes multipole and dipole polarization first */
/*     derivatives for a single atom with respect to Cartesian */
/*     coordinates; used to get finite difference second derivatives */


/* Subroutine */ int empole2a_(integer *nlist, integer *list, logical *
	reinduce)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal f;
    static integer i__, j, k;
    static doublereal r__, r2, ci;
    static integer ii, kk;
    static doublereal di[3], ck, dk[3], qi[9], gl[9], qk[9], sc[10], gf[7];
    static integer ix, iz, kx, kz;
    static doublereal xr, yr, zr, rr1, rr3, rr5, rr7, rr9, gfd, gfi[6];
    static integer iii;
    static doublereal pdi, rr11, pti, qir[3], qkr[3], gli[7], sci[8], gti[6], 
	    dsc3, dsc5, dsc7, psc3, ftm2[3], psc5, psc7, ttm2[3], ttm3[3], 
	    damp, fdir[3], qidk[3], qkdi[3], glip[7], fgrp, qiuk[3], qkui[3], 
	    dixr[3], dkxr[3], scip[8];
    static logical usei, usek;
    static doublereal ddsc3[3], ddsc5[3], ddsc7[3], ftm2i[3], temp3, temp5, 
	    temp7, ttm2i[3], ttm3i[3];
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal dixdk[3], frcxi[3], frcxk[3], frcyi[3], frcyk[3], frczi[
	    3], frczk[3], dixuk[3], dkxui[3], qiukp[3], qkuip[3], qiqkr[3], 
	    qkqir[3], qixqk[3], rxqir[3], rxqkr[3], scale3, scale5, scale7, 
	    dscale[25000], pgamma, mscale[25000], pscale[25000];
    extern /* Subroutine */ int induce_(void);
    static doublereal uscale[25000], findmp[3], fridmp[3], dixukp[3], dkxuip[
	    3], dixqkr[3], scale3i, scale5i, scale7i, uixqkr[3], ukxqir[3], 
	    rxqiuk[3], rxqkui[3], dkxqir[3], rxqikr[3], rxqkir[3], rxqidk[3], 
	    rxqkdi[3];
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *), torque_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static logical proceed;
    extern /* Subroutine */ int chkpole_(void);
    static doublereal expdamp;
    extern /* Subroutine */ int rotpole_(void);
    static doublereal qkrxqir[3], uixqkrp[3], ukxqirp[3], rxqiukp[3], rxqkuip[
	    3];


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define dem_ref(a_1,a_2) deriv_1.dem[(a_2)*3 + a_1 - 4]
#define dep_ref(a_1,a_2) deriv_1.dep[(a_2)*3 + a_1 - 4]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define ip12_ref(a_1,a_2) polgrp_1.ip12[(a_2)*50 + a_1 - 51]
#define ip13_ref(a_1,a_2) polgrp_1.ip13[(a_2)*50 + a_1 - 51]
#define ip14_ref(a_1,a_2) polgrp_1.ip14[(a_2)*50 + a_1 - 51]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define uinp_ref(a_1,a_2) polar_1.uinp[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  chgpot.i  --  specifics of charge-charge functional form  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     electric   energy factor in kcal/mole for current force field */
/*     dielec     dielectric constant for electrostatic interactions */
/*     ebuffer    electrostatic buffering constant added to distance */
/*     c2scale    factor by which 1-2 charge interactions are scaled */
/*     c3scale    factor by which 1-3 charge interactions are scaled */
/*     c4scale    factor by which 1-4 charge interactions are scaled */
/*     c5scale    factor by which 1-5 charge interactions are scaled */
/*     neutnbr    logical flag governing use of neutral group neighbors */
/*     neutcut    logical flag governing use of neutral group cutoffs */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cutoff.i  --  cutoff distances for energy interactions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     vdwcut      cutoff distance for van der Waals interactions */
/*     chgcut      cutoff distance for charge-charge interactions */
/*     dplcut      cutoff distance for dipole-dipole interactions */
/*     mpolecut    cutoff distance for atomic multipole interactions */
/*     vdwtaper    distance at which van der Waals switching begins */
/*     chgtaper    distance at which charge-charge switching begins */
/*     dpltaper    distance at which dipole-dipole switching begins */
/*     mpoletaper  distance at which atomic multipole switching begins */
/*     ewaldcut    cutoff distance for direct space Ewald summation */
/*     use_ewald   logical flag governing use of Ewald summation */
/*     use_lights  logical flag governing use of method of lights */
/*     use_list    logical flag governing use of any neighbor lists */
/*     use_vlist   logical flag governing use of vdw neighbor lists */
/*     use_clist   logical flag governing use of charge neighbor lists */
/*     use_mlist   logical flag governing use of multipole neighbor lists */




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
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mplpot.i  --  specifics of atomic multipole functions  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     m2scale   factor by which 1-2 multipole interactions are scaled */
/*     m3scale   factor by which 1-3 multipole interactions are scaled */
/*     m4scale   factor by which 1-4 multipole interactions are scaled */
/*     m5scale   factor by which 1-5 multipole interactions are scaled */




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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  polar.i  --  polarizabilities and induced dipole moments  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     polarity  dipole polarizability for each multipole site (Ang**3) */
/*     thole     Thole polarizability damping value for each site */
/*     pdamp     value of polarizability scale factor for each site */
/*     uind      induced dipole components at each multipole site */
/*     uinp      induced dipoles in field used for energy interactions */
/*     uinds     GK or PB induced dipoles at each multipole site */
/*     uinps     induced dipoles in field used for GK or PB energy */
/*     npolar    total number of polarizable sites in the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  polgrp.i  --  polarizable site group connectivity lists   ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxp11   maximum number of atoms in a polarization group */
/*     maxp12   maximum number of atoms in groups 1-2 to an atom */
/*     maxp13   maximum number of atoms in groups 1-3 to an atom */
/*     maxp14   maximum number of atoms in groups 1-4 to an atom */

/*     np11     number of atoms in polarization group of each atom */
/*     ip11     atom numbers of atoms in same group as each atom */
/*     np12     number of atoms in groups 1-2 to each atom */
/*     ip12     atom numbers of atoms in groups 1-2 to each atom */
/*     np13     number of atoms in groups 1-3 to each atom */
/*     ip13     atom numbers of atoms in groups 1-3 to each atom */
/*     np14     number of atoms in groups 1-4 to each atom */
/*     ip14     atom numbers of atoms in groups 1-4 to each atom */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  polpot.i  --  specifics of polarization functional form  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     poleps    induced dipole convergence criterion (rms Debyes/atom) */
/*     polsor    induced dipole SOR convergence acceleration factor */
/*     p2scale   field 1-2 scale factor for energy evaluations */
/*     p3scale   field 1-3 scale factor for energy evaluations */
/*     p4scale   field 1-4 scale factor for energy evaluations */
/*     p5scale   field 1-5 scale factor for energy evaluations */
/*     d1scale   field intra-group scale factor for direct induced */
/*     d2scale   field 1-2 group scale factor for direct induced */
/*     d3scale   field 1-3 group scale factor for direct induced */
/*     d4scale   field 1-4 group scale factor for direct induced */
/*     u1scale   field intra-group scale factor for mutual induced */
/*     u2scale   field 1-2 group scale factor for mutual induced */
/*     u3scale   field 1-3 group scale factor for mutual induced */
/*     u4scale   field 1-4 group scale factor for mutual induced */
/*     poltyp    type of polarization potential (direct or mutual) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  potent.i  --  usage of each potential energy component  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     use_bond    logical flag governing use of bond stretch potential */
/*     use_angle   logical flag governing use of angle bend potential */
/*     use_strbnd  logical flag governing use of stretch-bend potential */
/*     use_urey    logical flag governing use of Urey-Bradley potential */
/*     use_angang  logical flag governing use of angle-angle cross term */
/*     use_opbend  logical flag governing use of out-of-plane bend term */
/*     use_opdist  logical flag governing use of out-of-plane distance */
/*     use_improp  logical flag governing use of improper dihedral term */
/*     use_imptor  logical flag governing use of improper torsion term */
/*     use_tors    logical flag governing use of torsional potential */
/*     use_pitors  logical flag governing use of pi-orbital torsion term */
/*     use_strtor  logical flag governing use of stretch-torsion term */
/*     use_tortor  logical flag governing use of torsion-torsion term */
/*     use_vdw     logical flag governing use of vdw der Waals potential */
/*     use_charge  logical flag governing use of charge-charge potential */
/*     use_chgdpl  logical flag governing use of charge-dipole potential */
/*     use_dipole  logical flag governing use of dipole-dipole potential */
/*     use_mpole   logical flag governing use of multipole potential */
/*     use_polar   logical flag governing use of polarization term */
/*     use_rxnfld  logical flag governing use of reaction field term */
/*     use_solv    logical flag governing use of continuum solvation */
/*     use_metal   logical flag governing use of ligand field term */
/*     use_geom    logical flag governing use of geometric restraints */
/*     use_extra   logical flag governing use of extra potential term */
/*     use_born    logical flag governing use of Born radii values */
/*     use_orbit   logical flag governing use of pisystem computation */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  shunt.i  --  polynomial switching function coefficients  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     off    distance at which the potential energy goes to zero */
/*     off2   square of distance at which the potential goes to zero */
/*     cut    distance at which switching of the potential begins */
/*     cut2   square of distance at which the switching begins */
/*     c0     zeroth order coefficient of multiplicative switch */
/*     c1     first order coefficient of multiplicative switch */
/*     c2     second order coefficient of multiplicative switch */
/*     c3     third order coefficient of multiplicative switch */
/*     c4     fourth order coefficient of multiplicative switch */
/*     c5     fifth order coefficient of multiplicative switch */
/*     f0     zeroth order coefficient of additive switch function */
/*     f1     first order coefficient of additive switch function */
/*     f2     second order coefficient of additive switch function */
/*     f3     third order coefficient of additive switch function */
/*     f4     fourth order coefficient of additive switch function */
/*     f5     fifth order coefficient of additive switch function */
/*     f6     sixth order coefficient of additive switch function */
/*     f7     seventh order coefficient of additive switch function */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     zero out the multipole and polarization first derivatives */

    /* Parameter adjustments */
    --list;

    /* Function Body */
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    dem_ref(j, i__) = 0.;
	    dep_ref(j, i__) = 0.;
	}
    }
    if (*nlist == 0) {
	return 0;
    }

/*     set arrays needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mscale[i__ - 1] = 1.;
	pscale[i__ - 1] = 1.;
	dscale[i__ - 1] = 1.;
	uscale[i__ - 1] = 1.;
    }

/*     set conversion factor, cutoff and switching coefficients */

    f = chgpot_1.electric / chgpot_1.dielec;
    switch_("MPOLE", (ftnlen)5);

/*     check the sign of multipole components at chiral sites */

    if (*reinduce) {
	chkpole_();

/*     rotate the multipole components into the global frame */

	rotpole_();

/*     compute the induced dipoles at each polarizable atom */

	induce_();
    }

/*     set scale factors for permanent multipole and induced terms */

    i__1 = *nlist;
    for (iii = 1; iii <= i__1; ++iii) {
	i__ = list[iii];
	ii = mpole_1.ipole[i__ - 1];
	iz = mpole_1.zaxis[i__ - 1];
	ix = mpole_1.xaxis[i__ - 1];
	pdi = polar_1.pdamp[i__ - 1];
	pti = polar_1.thole[i__ - 1];
	ci = rpole_ref(1, i__);
	di[0] = rpole_ref(2, i__);
	di[1] = rpole_ref(3, i__);
	di[2] = rpole_ref(4, i__);
	qi[0] = rpole_ref(5, i__);
	qi[1] = rpole_ref(6, i__);
	qi[2] = rpole_ref(7, i__);
	qi[3] = rpole_ref(8, i__);
	qi[4] = rpole_ref(9, i__);
	qi[5] = rpole_ref(10, i__);
	qi[6] = rpole_ref(11, i__);
	qi[7] = rpole_ref(12, i__);
	qi[8] = rpole_ref(13, i__);
	usei = usage_1.use[ii - 1] || usage_1.use[iz - 1] || usage_1.use[ix - 
		1];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = mplpot_1.m2scale;
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = mplpot_1.m3scale;
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = mplpot_1.m4scale;
	    pscale[i14_ref(j, ii) - 1] = polpot_1.p4scale;
	    i__3 = polgrp_1.np11[ii - 1];
	    for (k = 1; k <= i__3; ++k) {
		if (i14_ref(j, ii) == ip11_ref(k, ii)) {
		    pscale[i14_ref(j, ii) - 1] = pscale[i14_ref(j, ii) - 1] * 
			    .5;
		}
	    }
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i15_ref(j, ii) - 1] = mplpot_1.m5scale;
	    pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = polpot_1.d1scale;
	    uscale[ip11_ref(j, ii) - 1] = polpot_1.u1scale;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = polpot_1.d2scale;
	    uscale[ip12_ref(j, ii) - 1] = polpot_1.u2scale;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = polpot_1.d3scale;
	    uscale[ip13_ref(j, ii) - 1] = polpot_1.u3scale;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = polpot_1.d4scale;
	    uscale[ip14_ref(j, ii) - 1] = polpot_1.u4scale;
	}
	i__2 = mpole_1.npole;
	for (k = i__ + 1; k <= i__2; ++k) {
	    kk = mpole_1.ipole[k - 1];
	    kz = mpole_1.zaxis[k - 1];
	    kx = mpole_1.xaxis[k - 1];
	    usek = usage_1.use[kk - 1] || usage_1.use[kz - 1] || usage_1.use[
		    kx - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &ii, &kk, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (! group_1.use_intra__) {
		proceed = TRUE_;
	    }
	    if (proceed) {
		proceed = usei || usek;
	    }
	    if (! proceed) {
		goto L10;
	    }
	    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
	    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
	    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
	    image_(&xr, &yr, &zr);
	    r2 = xr * xr + yr * yr + zr * zr;
	    if (r2 <= shunt_1.off2) {
		r__ = sqrt(r2);
		ck = rpole_ref(1, k);
		dk[0] = rpole_ref(2, k);
		dk[1] = rpole_ref(3, k);
		dk[2] = rpole_ref(4, k);
		qk[0] = rpole_ref(5, k);
		qk[1] = rpole_ref(6, k);
		qk[2] = rpole_ref(7, k);
		qk[3] = rpole_ref(8, k);
		qk[4] = rpole_ref(9, k);
		qk[5] = rpole_ref(10, k);
		qk[6] = rpole_ref(11, k);
		qk[7] = rpole_ref(12, k);
		qk[8] = rpole_ref(13, k);

/*     apply Thole polarization damping to scale factors */

		rr1 = 1. / r__;
		rr3 = rr1 / r2;
		rr5 = rr3 * 3. / r2;
		rr7 = rr5 * 5. / r2;
		rr9 = rr7 * 7. / r2;
		rr11 = rr9 * 9. / r2;
		scale3 = 1.;
		scale5 = 1.;
		scale7 = 1.;
		for (j = 1; j <= 3; ++j) {
		    ddsc3[j - 1] = 0.;
		    ddsc5[j - 1] = 0.;
		    ddsc7[j - 1] = 0.;
		}
		damp = pdi * polar_1.pdamp[k - 1];
		if (damp != 0.) {
/* Computing MIN */
		    d__1 = pti, d__2 = polar_1.thole[k - 1];
		    pgamma = min(d__1,d__2);
/* Computing 3rd power */
		    d__1 = r__ / damp;
		    damp = -pgamma * (d__1 * (d__1 * d__1));
		    if (damp > -50.) {
			expdamp = exp(damp);
			scale3 = 1. - expdamp;
			scale5 = 1. - (1. - damp) * expdamp;
/* Computing 2nd power */
			d__1 = damp;
			scale7 = 1. - (1. - damp + d__1 * d__1 * .6) * 
				expdamp;
			temp3 = damp * -3. * expdamp / r2;
			temp5 = -damp;
			temp7 = -.2 - damp * .6;
			ddsc3[0] = temp3 * xr;
			ddsc3[1] = temp3 * yr;
			ddsc3[2] = temp3 * zr;
			ddsc5[0] = temp5 * ddsc3[0];
			ddsc5[1] = temp5 * ddsc3[1];
			ddsc5[2] = temp5 * ddsc3[2];
			ddsc7[0] = temp7 * ddsc5[0];
			ddsc7[1] = temp7 * ddsc5[1];
			ddsc7[2] = temp7 * ddsc5[2];
		    }
		}
		scale3i = scale3 * uscale[kk - 1];
		scale5i = scale5 * uscale[kk - 1];
		scale7i = scale7 * uscale[kk - 1];
		dsc3 = scale3 * dscale[kk - 1];
		dsc5 = scale5 * dscale[kk - 1];
		dsc7 = scale7 * dscale[kk - 1];
		psc3 = scale3 * pscale[kk - 1];
		psc5 = scale5 * pscale[kk - 1];
		psc7 = scale7 * pscale[kk - 1];

/*     construct necessary auxiliary vectors */

		dixdk[0] = di[1] * dk[2] - di[2] * dk[1];
		dixdk[1] = di[2] * dk[0] - di[0] * dk[2];
		dixdk[2] = di[0] * dk[1] - di[1] * dk[0];
		dixuk[0] = di[1] * uind_ref(3, k) - di[2] * uind_ref(2, k);
		dixuk[1] = di[2] * uind_ref(1, k) - di[0] * uind_ref(3, k);
		dixuk[2] = di[0] * uind_ref(2, k) - di[1] * uind_ref(1, k);
		dkxui[0] = dk[1] * uind_ref(3, i__) - dk[2] * uind_ref(2, i__)
			;
		dkxui[1] = dk[2] * uind_ref(1, i__) - dk[0] * uind_ref(3, i__)
			;
		dkxui[2] = dk[0] * uind_ref(2, i__) - dk[1] * uind_ref(1, i__)
			;
		dixukp[0] = di[1] * uinp_ref(3, k) - di[2] * uinp_ref(2, k);
		dixukp[1] = di[2] * uinp_ref(1, k) - di[0] * uinp_ref(3, k);
		dixukp[2] = di[0] * uinp_ref(2, k) - di[1] * uinp_ref(1, k);
		dkxuip[0] = dk[1] * uinp_ref(3, i__) - dk[2] * uinp_ref(2, 
			i__);
		dkxuip[1] = dk[2] * uinp_ref(1, i__) - dk[0] * uinp_ref(3, 
			i__);
		dkxuip[2] = dk[0] * uinp_ref(2, i__) - dk[1] * uinp_ref(1, 
			i__);
		dixr[0] = di[1] * zr - di[2] * yr;
		dixr[1] = di[2] * xr - di[0] * zr;
		dixr[2] = di[0] * yr - di[1] * xr;
		dkxr[0] = dk[1] * zr - dk[2] * yr;
		dkxr[1] = dk[2] * xr - dk[0] * zr;
		dkxr[2] = dk[0] * yr - dk[1] * xr;
		qir[0] = qi[0] * xr + qi[3] * yr + qi[6] * zr;
		qir[1] = qi[1] * xr + qi[4] * yr + qi[7] * zr;
		qir[2] = qi[2] * xr + qi[5] * yr + qi[8] * zr;
		qkr[0] = qk[0] * xr + qk[3] * yr + qk[6] * zr;
		qkr[1] = qk[1] * xr + qk[4] * yr + qk[7] * zr;
		qkr[2] = qk[2] * xr + qk[5] * yr + qk[8] * zr;
		qiqkr[0] = qi[0] * qkr[0] + qi[3] * qkr[1] + qi[6] * qkr[2];
		qiqkr[1] = qi[1] * qkr[0] + qi[4] * qkr[1] + qi[7] * qkr[2];
		qiqkr[2] = qi[2] * qkr[0] + qi[5] * qkr[1] + qi[8] * qkr[2];
		qkqir[0] = qk[0] * qir[0] + qk[3] * qir[1] + qk[6] * qir[2];
		qkqir[1] = qk[1] * qir[0] + qk[4] * qir[1] + qk[7] * qir[2];
		qkqir[2] = qk[2] * qir[0] + qk[5] * qir[1] + qk[8] * qir[2];
		qixqk[0] = qi[1] * qk[2] + qi[4] * qk[5] + qi[7] * qk[8] - qi[
			2] * qk[1] - qi[5] * qk[4] - qi[8] * qk[7];
		qixqk[1] = qi[2] * qk[0] + qi[5] * qk[3] + qi[8] * qk[6] - qi[
			0] * qk[2] - qi[3] * qk[5] - qi[6] * qk[8];
		qixqk[2] = qi[0] * qk[1] + qi[3] * qk[4] + qi[6] * qk[7] - qi[
			1] * qk[0] - qi[4] * qk[3] - qi[7] * qk[6];
		rxqir[0] = yr * qir[2] - zr * qir[1];
		rxqir[1] = zr * qir[0] - xr * qir[2];
		rxqir[2] = xr * qir[1] - yr * qir[0];
		rxqkr[0] = yr * qkr[2] - zr * qkr[1];
		rxqkr[1] = zr * qkr[0] - xr * qkr[2];
		rxqkr[2] = xr * qkr[1] - yr * qkr[0];
		rxqikr[0] = yr * qiqkr[2] - zr * qiqkr[1];
		rxqikr[1] = zr * qiqkr[0] - xr * qiqkr[2];
		rxqikr[2] = xr * qiqkr[1] - yr * qiqkr[0];
		rxqkir[0] = yr * qkqir[2] - zr * qkqir[1];
		rxqkir[1] = zr * qkqir[0] - xr * qkqir[2];
		rxqkir[2] = xr * qkqir[1] - yr * qkqir[0];
		qkrxqir[0] = qkr[1] * qir[2] - qkr[2] * qir[1];
		qkrxqir[1] = qkr[2] * qir[0] - qkr[0] * qir[2];
		qkrxqir[2] = qkr[0] * qir[1] - qkr[1] * qir[0];
		qidk[0] = qi[0] * dk[0] + qi[3] * dk[1] + qi[6] * dk[2];
		qidk[1] = qi[1] * dk[0] + qi[4] * dk[1] + qi[7] * dk[2];
		qidk[2] = qi[2] * dk[0] + qi[5] * dk[1] + qi[8] * dk[2];
		qkdi[0] = qk[0] * di[0] + qk[3] * di[1] + qk[6] * di[2];
		qkdi[1] = qk[1] * di[0] + qk[4] * di[1] + qk[7] * di[2];
		qkdi[2] = qk[2] * di[0] + qk[5] * di[1] + qk[8] * di[2];
		qiuk[0] = qi[0] * uind_ref(1, k) + qi[3] * uind_ref(2, k) + 
			qi[6] * uind_ref(3, k);
		qiuk[1] = qi[1] * uind_ref(1, k) + qi[4] * uind_ref(2, k) + 
			qi[7] * uind_ref(3, k);
		qiuk[2] = qi[2] * uind_ref(1, k) + qi[5] * uind_ref(2, k) + 
			qi[8] * uind_ref(3, k);
		qkui[0] = qk[0] * uind_ref(1, i__) + qk[3] * uind_ref(2, i__) 
			+ qk[6] * uind_ref(3, i__);
		qkui[1] = qk[1] * uind_ref(1, i__) + qk[4] * uind_ref(2, i__) 
			+ qk[7] * uind_ref(3, i__);
		qkui[2] = qk[2] * uind_ref(1, i__) + qk[5] * uind_ref(2, i__) 
			+ qk[8] * uind_ref(3, i__);
		qiukp[0] = qi[0] * uinp_ref(1, k) + qi[3] * uinp_ref(2, k) + 
			qi[6] * uinp_ref(3, k);
		qiukp[1] = qi[1] * uinp_ref(1, k) + qi[4] * uinp_ref(2, k) + 
			qi[7] * uinp_ref(3, k);
		qiukp[2] = qi[2] * uinp_ref(1, k) + qi[5] * uinp_ref(2, k) + 
			qi[8] * uinp_ref(3, k);
		qkuip[0] = qk[0] * uinp_ref(1, i__) + qk[3] * uinp_ref(2, i__)
			 + qk[6] * uinp_ref(3, i__);
		qkuip[1] = qk[1] * uinp_ref(1, i__) + qk[4] * uinp_ref(2, i__)
			 + qk[7] * uinp_ref(3, i__);
		qkuip[2] = qk[2] * uinp_ref(1, i__) + qk[5] * uinp_ref(2, i__)
			 + qk[8] * uinp_ref(3, i__);
		dixqkr[0] = di[1] * qkr[2] - di[2] * qkr[1];
		dixqkr[1] = di[2] * qkr[0] - di[0] * qkr[2];
		dixqkr[2] = di[0] * qkr[1] - di[1] * qkr[0];
		dkxqir[0] = dk[1] * qir[2] - dk[2] * qir[1];
		dkxqir[1] = dk[2] * qir[0] - dk[0] * qir[2];
		dkxqir[2] = dk[0] * qir[1] - dk[1] * qir[0];
		uixqkr[0] = uind_ref(2, i__) * qkr[2] - uind_ref(3, i__) * 
			qkr[1];
		uixqkr[1] = uind_ref(3, i__) * qkr[0] - uind_ref(1, i__) * 
			qkr[2];
		uixqkr[2] = uind_ref(1, i__) * qkr[1] - uind_ref(2, i__) * 
			qkr[0];
		ukxqir[0] = uind_ref(2, k) * qir[2] - uind_ref(3, k) * qir[1];
		ukxqir[1] = uind_ref(3, k) * qir[0] - uind_ref(1, k) * qir[2];
		ukxqir[2] = uind_ref(1, k) * qir[1] - uind_ref(2, k) * qir[0];
		uixqkrp[0] = uinp_ref(2, i__) * qkr[2] - uinp_ref(3, i__) * 
			qkr[1];
		uixqkrp[1] = uinp_ref(3, i__) * qkr[0] - uinp_ref(1, i__) * 
			qkr[2];
		uixqkrp[2] = uinp_ref(1, i__) * qkr[1] - uinp_ref(2, i__) * 
			qkr[0];
		ukxqirp[0] = uinp_ref(2, k) * qir[2] - uinp_ref(3, k) * qir[1]
			;
		ukxqirp[1] = uinp_ref(3, k) * qir[0] - uinp_ref(1, k) * qir[2]
			;
		ukxqirp[2] = uinp_ref(1, k) * qir[1] - uinp_ref(2, k) * qir[0]
			;
		rxqidk[0] = yr * qidk[2] - zr * qidk[1];
		rxqidk[1] = zr * qidk[0] - xr * qidk[2];
		rxqidk[2] = xr * qidk[1] - yr * qidk[0];
		rxqkdi[0] = yr * qkdi[2] - zr * qkdi[1];
		rxqkdi[1] = zr * qkdi[0] - xr * qkdi[2];
		rxqkdi[2] = xr * qkdi[1] - yr * qkdi[0];
		rxqiuk[0] = yr * qiuk[2] - zr * qiuk[1];
		rxqiuk[1] = zr * qiuk[0] - xr * qiuk[2];
		rxqiuk[2] = xr * qiuk[1] - yr * qiuk[0];
		rxqkui[0] = yr * qkui[2] - zr * qkui[1];
		rxqkui[1] = zr * qkui[0] - xr * qkui[2];
		rxqkui[2] = xr * qkui[1] - yr * qkui[0];
		rxqiukp[0] = yr * qiukp[2] - zr * qiukp[1];
		rxqiukp[1] = zr * qiukp[0] - xr * qiukp[2];
		rxqiukp[2] = xr * qiukp[1] - yr * qiukp[0];
		rxqkuip[0] = yr * qkuip[2] - zr * qkuip[1];
		rxqkuip[1] = zr * qkuip[0] - xr * qkuip[2];
		rxqkuip[2] = xr * qkuip[1] - yr * qkuip[0];

/*     calculate scalar products for permanent components */

		sc[1] = di[0] * dk[0] + di[1] * dk[1] + di[2] * dk[2];
		sc[2] = di[0] * xr + di[1] * yr + di[2] * zr;
		sc[3] = dk[0] * xr + dk[1] * yr + dk[2] * zr;
		sc[4] = qir[0] * xr + qir[1] * yr + qir[2] * zr;
		sc[5] = qkr[0] * xr + qkr[1] * yr + qkr[2] * zr;
		sc[6] = qir[0] * dk[0] + qir[1] * dk[1] + qir[2] * dk[2];
		sc[7] = qkr[0] * di[0] + qkr[1] * di[1] + qkr[2] * di[2];
		sc[8] = qir[0] * qkr[0] + qir[1] * qkr[1] + qir[2] * qkr[2];
		sc[9] = qi[0] * qk[0] + qi[1] * qk[1] + qi[2] * qk[2] + qi[3] 
			* qk[3] + qi[4] * qk[4] + qi[5] * qk[5] + qi[6] * qk[
			6] + qi[7] * qk[7] + qi[8] * qk[8];

/*     calculate scalar products for induced components */

		sci[0] = uind_ref(1, i__) * dk[0] + uind_ref(2, i__) * dk[1] 
			+ uind_ref(3, i__) * dk[2] + di[0] * uind_ref(1, k) + 
			di[1] * uind_ref(2, k) + di[2] * uind_ref(3, k);
		sci[1] = uind_ref(1, i__) * uind_ref(1, k) + uind_ref(2, i__) 
			* uind_ref(2, k) + uind_ref(3, i__) * uind_ref(3, k);
		sci[2] = uind_ref(1, i__) * xr + uind_ref(2, i__) * yr + 
			uind_ref(3, i__) * zr;
		sci[3] = uind_ref(1, k) * xr + uind_ref(2, k) * yr + uind_ref(
			3, k) * zr;
		sci[6] = qir[0] * uind_ref(1, k) + qir[1] * uind_ref(2, k) + 
			qir[2] * uind_ref(3, k);
		sci[7] = qkr[0] * uind_ref(1, i__) + qkr[1] * uind_ref(2, i__)
			 + qkr[2] * uind_ref(3, i__);
		scip[0] = uinp_ref(1, i__) * dk[0] + uinp_ref(2, i__) * dk[1] 
			+ uinp_ref(3, i__) * dk[2] + di[0] * uinp_ref(1, k) + 
			di[1] * uinp_ref(2, k) + di[2] * uinp_ref(3, k);
		scip[1] = uind_ref(1, i__) * uinp_ref(1, k) + uind_ref(2, i__)
			 * uinp_ref(2, k) + uind_ref(3, i__) * uinp_ref(3, k) 
			+ uinp_ref(1, i__) * uind_ref(1, k) + uinp_ref(2, i__)
			 * uind_ref(2, k) + uinp_ref(3, i__) * uind_ref(3, k);
		scip[2] = uinp_ref(1, i__) * xr + uinp_ref(2, i__) * yr + 
			uinp_ref(3, i__) * zr;
		scip[3] = uinp_ref(1, k) * xr + uinp_ref(2, k) * yr + 
			uinp_ref(3, k) * zr;
		scip[6] = qir[0] * uinp_ref(1, k) + qir[1] * uinp_ref(2, k) + 
			qir[2] * uinp_ref(3, k);
		scip[7] = qkr[0] * uinp_ref(1, i__) + qkr[1] * uinp_ref(2, 
			i__) + qkr[2] * uinp_ref(3, i__);

/*     calculate the gl functions for permanent components */

		gl[0] = ci * ck;
		gl[1] = ck * sc[2] - ci * sc[3];
		gl[2] = ci * sc[5] + ck * sc[4] - sc[2] * sc[3];
		gl[3] = sc[2] * sc[5] - sc[3] * sc[4];
		gl[4] = sc[4] * sc[5];
		gl[5] = sc[8] * -4.;
		gl[6] = sc[1];
		gl[7] = (sc[6] - sc[7]) * 2.;
		gl[8] = sc[9] * 2.;

/*     calculate the gl functions for induced components */

		gli[0] = ck * sci[2] - ci * sci[3];
		gli[1] = -sc[2] * sci[3] - sci[2] * sc[3];
		gli[2] = sci[2] * sc[5] - sci[3] * sc[4];
		gli[5] = sci[0];
		gli[6] = (sci[6] - sci[7]) * 2.;
		glip[0] = ck * scip[2] - ci * scip[3];
		glip[1] = -sc[2] * scip[3] - scip[2] * sc[3];
		glip[2] = scip[2] * sc[5] - scip[3] * sc[4];
		glip[5] = scip[0];
		glip[6] = (scip[6] - scip[7]) * 2.;

/*     intermediate variables for the permanent components */

		gf[0] = rr3 * gl[0] + rr5 * (gl[1] + gl[6]) + rr7 * (gl[2] + 
			gl[7] + gl[8]) + rr9 * (gl[3] + gl[5]) + rr11 * gl[4];
		gf[1] = -ck * rr3 + sc[3] * rr5 - sc[5] * rr7;
		gf[2] = ci * rr3 + sc[2] * rr5 + sc[4] * rr7;
		gf[3] = rr5 * 2.;
		gf[4] = (-ck * rr5 + sc[3] * rr7 - sc[5] * rr9) * 2.;
		gf[5] = (-ci * rr5 - sc[2] * rr7 - sc[4] * rr9) * 2.;
		gf[6] = rr7 * 4.;

/*     intermediate variables for the induced components */

		gfi[0] = rr5 * .5 * ((gli[0] + gli[5]) * psc3 + (glip[0] + 
			glip[5]) * dsc3 + scip[1] * scale3i) + rr7 * .5 * ((
			gli[6] + gli[1]) * psc5 + (glip[6] + glip[1]) * dsc5 
			- (sci[2] * scip[3] + scip[2] * sci[3]) * scale5i) + 
			rr9 * .5 * (gli[2] * psc7 + glip[2] * dsc7);
		gfi[1] = -rr3 * ck + rr5 * sc[3] - rr7 * sc[5];
		gfi[2] = rr3 * ci + rr5 * sc[2] + rr7 * sc[4];
		gfi[3] = rr5 * 2.;
		gfi[4] = rr7 * (sci[3] * psc7 + scip[3] * dsc7);
		gfi[5] = -rr7 * (sci[2] * psc7 + scip[2] * dsc7);

/*     get the permanent force components */

		ftm2[0] = gf[0] * xr + gf[1] * di[0] + gf[2] * dk[0] + gf[3] *
			 (qkdi[0] - qidk[0]) + gf[4] * qir[0] + gf[5] * qkr[0]
			 + gf[6] * (qiqkr[0] + qkqir[0]);
		ftm2[1] = gf[0] * yr + gf[1] * di[1] + gf[2] * dk[1] + gf[3] *
			 (qkdi[1] - qidk[1]) + gf[4] * qir[1] + gf[5] * qkr[1]
			 + gf[6] * (qiqkr[1] + qkqir[1]);
		ftm2[2] = gf[0] * zr + gf[1] * di[2] + gf[2] * dk[2] + gf[3] *
			 (qkdi[2] - qidk[2]) + gf[4] * qir[2] + gf[5] * qkr[2]
			 + gf[6] * (qiqkr[2] + qkqir[2]);

/*     get the induced force components */

		ftm2i[0] = gfi[0] * xr + (-rr3 * ck * (uind_ref(1, i__) * 
			psc3 + uinp_ref(1, i__) * dsc3) + rr5 * sc[3] * (
			uind_ref(1, i__) * psc5 + uinp_ref(1, i__) * dsc5) - 
			rr7 * sc[5] * (uind_ref(1, i__) * psc7 + uinp_ref(1, 
			i__) * dsc7)) * .5 + (rr3 * ci * (uind_ref(1, k) * 
			psc3 + uinp_ref(1, k) * dsc3) + rr5 * sc[2] * (
			uind_ref(1, k) * psc5 + uinp_ref(1, k) * dsc5) + rr7 *
			 sc[4] * (uind_ref(1, k) * psc7 + uinp_ref(1, k) * 
			dsc7)) * .5 + rr5 * scale5i * (sci[3] * uinp_ref(1, 
			i__) + scip[3] * uind_ref(1, i__) + sci[2] * uinp_ref(
			1, k) + scip[2] * uind_ref(1, k)) * .5 + (sci[3] * 
			psc5 + scip[3] * dsc5) * .5 * rr5 * di[0] + (sci[2] * 
			psc5 + scip[2] * dsc5) * .5 * rr5 * dk[0] + gfi[3] * 
			.5 * ((qkui[0] - qiuk[0]) * psc5 + (qkuip[0] - qiukp[
			0]) * dsc5) + gfi[4] * qir[0] + gfi[5] * qkr[0];
		ftm2i[1] = gfi[0] * yr + (-rr3 * ck * (uind_ref(2, i__) * 
			psc3 + uinp_ref(2, i__) * dsc3) + rr5 * sc[3] * (
			uind_ref(2, i__) * psc5 + uinp_ref(2, i__) * dsc5) - 
			rr7 * sc[5] * (uind_ref(2, i__) * psc7 + uinp_ref(2, 
			i__) * dsc7)) * .5 + (rr3 * ci * (uind_ref(2, k) * 
			psc3 + uinp_ref(2, k) * dsc3) + rr5 * sc[2] * (
			uind_ref(2, k) * psc5 + uinp_ref(2, k) * dsc5) + rr7 *
			 sc[4] * (uind_ref(2, k) * psc7 + uinp_ref(2, k) * 
			dsc7)) * .5 + rr5 * scale5i * (sci[3] * uinp_ref(2, 
			i__) + scip[3] * uind_ref(2, i__) + sci[2] * uinp_ref(
			2, k) + scip[2] * uind_ref(2, k)) * .5 + (sci[3] * 
			psc5 + scip[3] * dsc5) * .5 * rr5 * di[1] + (sci[2] * 
			psc5 + scip[2] * dsc5) * .5 * rr5 * dk[1] + gfi[3] * 
			.5 * ((qkui[1] - qiuk[1]) * psc5 + (qkuip[1] - qiukp[
			1]) * dsc5) + gfi[4] * qir[1] + gfi[5] * qkr[1];
		ftm2i[2] = gfi[0] * zr + (-rr3 * ck * (uind_ref(3, i__) * 
			psc3 + uinp_ref(3, i__) * dsc3) + rr5 * sc[3] * (
			uind_ref(3, i__) * psc5 + uinp_ref(3, i__) * dsc5) - 
			rr7 * sc[5] * (uind_ref(3, i__) * psc7 + uinp_ref(3, 
			i__) * dsc7)) * .5 + (rr3 * ci * (uind_ref(3, k) * 
			psc3 + uinp_ref(3, k) * dsc3) + rr5 * sc[2] * (
			uind_ref(3, k) * psc5 + uinp_ref(3, k) * dsc5) + rr7 *
			 sc[4] * (uind_ref(3, k) * psc7 + uinp_ref(3, k) * 
			dsc7)) * .5 + rr5 * scale5i * (sci[3] * uinp_ref(3, 
			i__) + scip[3] * uind_ref(3, i__) + sci[2] * uinp_ref(
			3, k) + scip[2] * uind_ref(3, k)) * .5 + (sci[3] * 
			psc5 + scip[3] * dsc5) * .5 * rr5 * di[2] + (sci[2] * 
			psc5 + scip[2] * dsc5) * .5 * rr5 * dk[2] + gfi[3] * 
			.5 * ((qkui[2] - qiuk[2]) * psc5 + (qkuip[2] - qiukp[
			2]) * dsc5) + gfi[4] * qir[2] + gfi[5] * qkr[2];

/*     account for partially excluded induced interactions */

		temp3 = rr3 * .5 * ((gli[0] + gli[5]) * pscale[kk - 1] + (
			glip[0] + glip[5]) * dscale[kk - 1]);
		temp5 = rr5 * .5 * ((gli[1] + gli[6]) * pscale[kk - 1] + (
			glip[1] + glip[6]) * dscale[kk - 1]);
		temp7 = rr7 * .5 * (gli[2] * pscale[kk - 1] + glip[2] * 
			dscale[kk - 1]);
		fridmp[0] = temp3 * ddsc3[0] + temp5 * ddsc5[0] + temp7 * 
			ddsc7[0];
		fridmp[1] = temp3 * ddsc3[1] + temp5 * ddsc5[1] + temp7 * 
			ddsc7[1];
		fridmp[2] = temp3 * ddsc3[2] + temp5 * ddsc5[2] + temp7 * 
			ddsc7[2];

/*     find some scaling terms for induced-induced force */

		temp3 = rr3 * .5 * uscale[kk - 1] * scip[1];
		temp5 = rr5 * -.5 * uscale[kk - 1] * (sci[2] * scip[3] + scip[
			2] * sci[3]);
		findmp[0] = temp3 * ddsc3[0] + temp5 * ddsc5[0];
		findmp[1] = temp3 * ddsc3[1] + temp5 * ddsc5[1];
		findmp[2] = temp3 * ddsc3[2] + temp5 * ddsc5[2];

/*     modify induced force for partially excluded interactions */

		ftm2i[0] = ftm2i[0] - fridmp[0] - findmp[0];
		ftm2i[1] = ftm2i[1] - fridmp[1] - findmp[1];
		ftm2i[2] = ftm2i[2] - fridmp[2] - findmp[2];

/*     correction to convert mutual to direct polarization force */

		if (s_cmp(polpot_1.poltyp, "DIRECT", (ftnlen)6, (ftnlen)6) == 
			0) {
		    gfd = (rr5 * scip[1] * scale3i - rr7 * (scip[2] * sci[3] 
			    + sci[2] * scip[3]) * scale5i) * .5;
		    temp5 = rr5 * .5 * scale5i;
		    fdir[0] = gfd * xr + temp5 * (sci[3] * uinp_ref(1, i__) + 
			    scip[3] * uind_ref(1, i__) + sci[2] * uinp_ref(1, 
			    k) + scip[2] * uind_ref(1, k));
		    fdir[1] = gfd * yr + temp5 * (sci[3] * uinp_ref(2, i__) + 
			    scip[3] * uind_ref(2, i__) + sci[2] * uinp_ref(2, 
			    k) + scip[2] * uind_ref(2, k));
		    fdir[2] = gfd * zr + temp5 * (sci[3] * uinp_ref(3, i__) + 
			    scip[3] * uind_ref(3, i__) + sci[2] * uinp_ref(3, 
			    k) + scip[2] * uind_ref(3, k));
		    ftm2i[0] = ftm2i[0] - fdir[0] + findmp[0];
		    ftm2i[1] = ftm2i[1] - fdir[1] + findmp[1];
		    ftm2i[2] = ftm2i[2] - fdir[2] + findmp[2];
		}

/*     intermediate variables for induced torque components */

		gti[1] = rr5 * .5 * (sci[3] * psc5 + scip[3] * dsc5);
		gti[2] = rr5 * .5 * (sci[2] * psc5 + scip[2] * dsc5);
		gti[3] = gfi[3];
		gti[4] = gfi[4];
		gti[5] = gfi[5];

/*     get the permanent torque components */

		ttm2[0] = -rr3 * dixdk[0] + gf[1] * dixr[0] - gf[4] * rxqir[0]
			 + gf[3] * (dixqkr[0] + dkxqir[0] + rxqidk[0] - qixqk[
			0] * 2.) - gf[6] * (rxqikr[0] + qkrxqir[0]);
		ttm2[1] = -rr3 * dixdk[1] + gf[1] * dixr[1] - gf[4] * rxqir[1]
			 + gf[3] * (dixqkr[1] + dkxqir[1] + rxqidk[1] - qixqk[
			1] * 2.) - gf[6] * (rxqikr[1] + qkrxqir[1]);
		ttm2[2] = -rr3 * dixdk[2] + gf[1] * dixr[2] - gf[4] * rxqir[2]
			 + gf[3] * (dixqkr[2] + dkxqir[2] + rxqidk[2] - qixqk[
			2] * 2.) - gf[6] * (rxqikr[2] + qkrxqir[2]);
		ttm3[0] = rr3 * dixdk[0] + gf[2] * dkxr[0] - gf[5] * rxqkr[0] 
			- gf[3] * (dixqkr[0] + dkxqir[0] + rxqkdi[0] - qixqk[
			0] * 2.) - gf[6] * (rxqkir[0] - qkrxqir[0]);
		ttm3[1] = rr3 * dixdk[1] + gf[2] * dkxr[1] - gf[5] * rxqkr[1] 
			- gf[3] * (dixqkr[1] + dkxqir[1] + rxqkdi[1] - qixqk[
			1] * 2.) - gf[6] * (rxqkir[1] - qkrxqir[1]);
		ttm3[2] = rr3 * dixdk[2] + gf[2] * dkxr[2] - gf[5] * rxqkr[2] 
			- gf[3] * (dixqkr[2] + dkxqir[2] + rxqkdi[2] - qixqk[
			2] * 2.) - gf[6] * (rxqkir[2] - qkrxqir[2]);

/*     get the induced torque components */

		ttm2i[0] = -rr3 * (dixuk[0] * psc3 + dixukp[0] * dsc3) * .5 + 
			gti[1] * dixr[0] + gti[3] * ((ukxqir[0] + rxqiuk[0]) *
			 psc5 + (ukxqirp[0] + rxqiukp[0]) * dsc5) * .5 - gti[
			4] * rxqir[0];
		ttm2i[1] = -rr3 * (dixuk[1] * psc3 + dixukp[1] * dsc3) * .5 + 
			gti[1] * dixr[1] + gti[3] * ((ukxqir[1] + rxqiuk[1]) *
			 psc5 + (ukxqirp[1] + rxqiukp[1]) * dsc5) * .5 - gti[
			4] * rxqir[1];
		ttm2i[2] = -rr3 * (dixuk[2] * psc3 + dixukp[2] * dsc3) * .5 + 
			gti[1] * dixr[2] + gti[3] * ((ukxqir[2] + rxqiuk[2]) *
			 psc5 + (ukxqirp[2] + rxqiukp[2]) * dsc5) * .5 - gti[
			4] * rxqir[2];
		ttm3i[0] = -rr3 * (dkxui[0] * psc3 + dkxuip[0] * dsc3) * .5 + 
			gti[2] * dkxr[0] - gti[3] * ((uixqkr[0] + rxqkui[0]) *
			 psc5 + (uixqkrp[0] + rxqkuip[0]) * dsc5) * .5 - gti[
			5] * rxqkr[0];
		ttm3i[1] = -rr3 * (dkxui[1] * psc3 + dkxuip[1] * dsc3) * .5 + 
			gti[2] * dkxr[1] - gti[3] * ((uixqkr[1] + rxqkui[1]) *
			 psc5 + (uixqkrp[1] + rxqkuip[1]) * dsc5) * .5 - gti[
			5] * rxqkr[1];
		ttm3i[2] = -rr3 * (dkxui[2] * psc3 + dkxuip[2] * dsc3) * .5 + 
			gti[2] * dkxr[2] - gti[3] * ((uixqkr[2] + rxqkui[2]) *
			 psc5 + (uixqkrp[2] + rxqkuip[2]) * dsc5) * .5 - gti[
			5] * rxqkr[2];

/*     handle the case where scaling is used */

		for (j = 1; j <= 3; ++j) {
		    ftm2[j - 1] = f * ftm2[j - 1] * mscale[kk - 1];
		    ftm2i[j - 1] = f * ftm2i[j - 1];
		    ttm2[j - 1] = f * ttm2[j - 1] * mscale[kk - 1];
		    ttm2i[j - 1] = f * ttm2i[j - 1];
		    ttm3[j - 1] = f * ttm3[j - 1] * mscale[kk - 1];
		    ttm3i[j - 1] = f * ttm3i[j - 1];
		}

/*     increment gradient due to force and torque on first site */

		dem_ref(1, ii) = dem_ref(1, ii) + ftm2[0];
		dem_ref(2, ii) = dem_ref(2, ii) + ftm2[1];
		dem_ref(3, ii) = dem_ref(3, ii) + ftm2[2];
		dep_ref(1, ii) = dep_ref(1, ii) + ftm2i[0];
		dep_ref(2, ii) = dep_ref(2, ii) + ftm2i[1];
		dep_ref(3, ii) = dep_ref(3, ii) + ftm2i[2];
		torque_(&i__, ttm2, ttm2i, frcxi, frcyi, frczi);

/*     increment gradient due to force and torque on second site */

		dem_ref(1, kk) = dem_ref(1, kk) - ftm2[0];
		dem_ref(2, kk) = dem_ref(2, kk) - ftm2[1];
		dem_ref(3, kk) = dem_ref(3, kk) - ftm2[2];
		dep_ref(1, kk) = dep_ref(1, kk) - ftm2i[0];
		dep_ref(2, kk) = dep_ref(2, kk) - ftm2i[1];
		dep_ref(3, kk) = dep_ref(3, kk) - ftm2i[2];
		torque_(&k, ttm3, ttm3i, frcxk, frcyk, frczk);
	    }
L10:
	    ;
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = 1.;
	    pscale[i12_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = 1.;
	    pscale[i13_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = 1.;
	    pscale[i14_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i15_ref(j, ii) - 1] = 1.;
	    pscale[i15_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = 1.;
	    uscale[ip11_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = 1.;
	    uscale[ip12_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = 1.;
	    uscale[ip13_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = 1.;
	    uscale[ip14_ref(j, ii) - 1] = 1.;
	}
    }

/*     rezero out the derivatives for terms that are not used */

    if (! potent_1.use_mpole__) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		dem_ref(j, i__) = 0.;
	    }
	}
    }
    if (! potent_1.use_polar__) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		dep_ref(j, i__) = 0.;
	    }
	}
    }
    return 0;
} /* empole2a_ */

#undef rpole_ref
#undef uinp_ref
#undef uind_ref
#undef ip14_ref
#undef ip13_ref
#undef ip12_ref
#undef ip11_ref
#undef dep_ref
#undef dem_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref


