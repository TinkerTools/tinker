/* empole1.f -- translated by f2c (version 20050501).
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
    doublereal vdwcut, chgcut, dplcut, mpolecut, vdwtaper, chgtaper, dpltaper,
	     mpoletaper, ewaldcut;
    logical use_ewald__, use_lights__, use_list__, use_vlist__, use_clist__, 
	    use_mlist__;
} cutoff_;

#define cutoff_1 cutoff_

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
    doublereal esum, eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, 
	    ebt, ett, ev, ec, ecd, ed, em, ep, er, es, elf, eg, ex;
} energi_;

#define energi_1 energi_

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
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

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
    doublereal xcell, ycell, zcell, xcell2, ycell2, zcell2;
    integer ncell, icell[30000]	/* was [3][10000] */;
} cell_;

#define cell_1 cell_

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
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    doublereal einter;
} inter_;

#define inter_1 inter_

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

struct {
    doublereal vir[9]	/* was [3][3] */;
} virial_;

#define virial_1 virial_

struct {
    doublereal lbuffer, lbuf2, vbuf2, cbuf2, mbuf2;
    integer nvlst[25000], vlst[45000000]	/* was [1800][25000] */, 
	    nelst[25000], elst[30000000]	/* was [1200][25000] */;
    logical dovlst, doclst, domlst;
} neigh_;

#define neigh_1 neigh_

struct {
    doublereal aewald;
    char boundary[7];
} ewald_;

#define ewald_1 ewald_

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

/* Table of constant values */

static integer c__0 = 0;



/*     ############################################################# */
/*     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine empole1  --  mpole/polar energy & derivatives  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "empole1" calculates the multipole and dipole polarization */
/*     energy and derivatives with respect to Cartesian coordinates */


/* Subroutine */ int empole1_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, ii;
    extern /* Subroutine */ int empole1a_(void), empole1b_(void), empole1c_(
	    void), empole1d_(void);


#define dem_ref(a_1,a_2) deriv_1.dem[(a_2)*3 + a_1 - 4]
#define dep_ref(a_1,a_2) deriv_1.dep[(a_2)*3 + a_1 - 4]



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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




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




/*     choose the method for summing over multipole interactions */

    if (cutoff_1.use_ewald__) {
	if (cutoff_1.use_mlist__) {
	    empole1d_();
	} else {
	    empole1c_();
	}
    } else {
	if (cutoff_1.use_mlist__) {
	    empole1b_();
	} else {
	    empole1a_();
	}
    }

/*     zero out energy and derivative terms which are not in use */

    if (! potent_1.use_mpole__) {
	energi_1.em = 0.;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    for (j = 1; j <= 3; ++j) {
		dem_ref(j, ii) = 0.;
	    }
	}
    }
    if (! potent_1.use_polar__) {
	energi_1.ep = 0.;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    for (j = 1; j <= 3; ++j) {
		dep_ref(j, ii) = 0.;
	    }
	}
    }
    return 0;
} /* empole1_ */

#undef dep_ref
#undef dem_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine empole1a  --  double loop multipole derivatives  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "empole1a" calculates the multipole and dipole polarization */
/*     energy and derivatives with respect to Cartesian coordinates */
/*     using a pairwise double loop */


/* Subroutine */ int empole1a_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__, r2;
    static integer ii, kk;
    static doublereal ei, ci;
    static integer ix;
    static doublereal di[3];
    static integer iz, kx, kz;
    static doublereal qi[9], ck, dk[3], qk[9], xr, yr, zr, gl[9], sc[10], gf[
	    7], rr1, rr3, rr5, rr7, rr9, gfd, gfi[6];
    static integer iax, iay, iaz, kax, kay, kaz;
    static doublereal pdi, pti, xix, yix, zix, xiy, yiy, ziy, xiz, yiz, ziz, 
	    xkx, ykx, zkx, xky, yky, zky, xkz, ykz, zkz, rr11, vxx, dsc3, 
	    dsc5, vyy, dsc7, vzz, vyx, vzx, vzy, qir[3], qkr[3], gli[7], psc3,
	     ftm2[3], psc5, sci[8], psc7, gti[6], ttm2[3], ttm3[3], damp, 
	    fdir[3], qidk[3], qkdi[3], glip[7], fgrp, qiuk[3], qkui[3], dixr[
	    3], dkxr[3], scip[8];
    static logical usei, usek;
    static doublereal ddsc3[3], ddsc5[3], ddsc7[3], ftm2i[3], temp3, temp5, 
	    temp7, ttm2i[3], ttm3i[3];
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static integer jcell;
    static doublereal dixdk[3], frcxi[3], frcxk[3], frcyi[3], frcyk[3], frczi[
	    3], frczk[3], dixuk[3], dkxui[3], qiukp[3], qkuip[3], qiqkr[3], 
	    qkqir[3], qixqk[3], rxqir[3], scale3, rxqkr[3], scale5, scale7, 
	    dscale[25000], pgamma, mscale[25000];
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal pscale[25000];
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
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cell.i  --  periodic boundaries using replicated cells  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xcell    length of the a-axis of the complete replicated cell */
/*     ycell    length of the b-axis of the complete replicated cell */
/*     zcell    length of the c-axis of the complete replicated cell */
/*     xcell2   half the length of the a-axis of the replicated cell */
/*     ycell2   half the length of the b-axis of the replicated cell */
/*     zcell2   half the length of the c-axis of the replicated cell */
/*     ncell    total number of cell replicates for periodic boundaries */
/*     icell    offset along axes for each replicate periodic cell */




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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  inter.i  --  sum of intermolecular energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     einter   total intermolecular potential energy */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  virial.i  --  components of internal virial tensor  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     vir    total internal virial Cartesian tensor components */




/*     zero out multipole and polarization energy and derivatives */

    energi_1.em = 0.;
    energi_1.ep = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    dem_ref(j, i__) = 0.;
	    dep_ref(j, i__) = 0.;
	}
    }

/*     check the sign of multipole components at chiral sites */

    chkpole_();

/*     rotate the multipole components into the global frame */

    rotpole_();

/*     compute the induced dipoles at each polarizable atom */

    induce_();

/*     set arrays needed to scale connected atom interactions */

    if (mpole_1.npole == 0) {
	return 0;
    }
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

/*     set scale factors for permanent multipole and induced terms */

    i__1 = mpole_1.npole - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
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

/*     compute the energy contributions for this interaction */

		e = rr1 * gl[0] + rr3 * (gl[1] + gl[6]) + rr5 * (gl[2] + gl[7]
			 + gl[8]) + rr7 * (gl[3] + gl[5]) + rr9 * gl[4];
		ei = (rr3 * (gli[0] + gli[5]) * psc3 + rr5 * (gli[1] + gli[6])
			 * psc5 + rr7 * gli[2] * psc7) * .5;
		e = f * mscale[kk - 1] * e;
		ei = f * ei;
		energi_1.em += e;
		energi_1.ep += ei;

/*     increment the total intermolecular energy */

		if (molcul_1.molcule[ii - 1] != molcul_1.molcule[kk - 1]) {
		    inter_1.einter = inter_1.einter + e + ei;
		}

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

/*     intermediate terms for induced torque on multipoles */

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

/*     increment the internal virial tensor components */

		iaz = mpole_1.zaxis[i__ - 1];
		iax = mpole_1.xaxis[i__ - 1];
		iay = mpole_1.yaxis[i__ - 1];
		kaz = mpole_1.zaxis[k - 1];
		kax = mpole_1.xaxis[k - 1];
		kay = mpole_1.yaxis[k - 1];
		if (iaz == 0) {
		    iaz = ii;
		}
		if (iax == 0) {
		    iax = ii;
		}
		if (iay == 0) {
		    iay = ii;
		}
		if (kaz == 0) {
		    kaz = kk;
		}
		if (kax == 0) {
		    kax = kk;
		}
		if (kay == 0) {
		    kay = kk;
		}
		xiz = atoms_1.x[iaz - 1] - atoms_1.x[ii - 1];
		yiz = atoms_1.y[iaz - 1] - atoms_1.y[ii - 1];
		ziz = atoms_1.z__[iaz - 1] - atoms_1.z__[ii - 1];
		xix = atoms_1.x[iax - 1] - atoms_1.x[ii - 1];
		yix = atoms_1.y[iax - 1] - atoms_1.y[ii - 1];
		zix = atoms_1.z__[iax - 1] - atoms_1.z__[ii - 1];
		xiy = atoms_1.x[iay - 1] - atoms_1.x[ii - 1];
		yiy = atoms_1.y[iay - 1] - atoms_1.y[ii - 1];
		ziy = atoms_1.z__[iay - 1] - atoms_1.z__[ii - 1];
		xkz = atoms_1.x[kaz - 1] - atoms_1.x[kk - 1];
		ykz = atoms_1.y[kaz - 1] - atoms_1.y[kk - 1];
		zkz = atoms_1.z__[kaz - 1] - atoms_1.z__[kk - 1];
		xkx = atoms_1.x[kax - 1] - atoms_1.x[kk - 1];
		ykx = atoms_1.y[kax - 1] - atoms_1.y[kk - 1];
		zkx = atoms_1.z__[kax - 1] - atoms_1.z__[kk - 1];
		xky = atoms_1.x[kay - 1] - atoms_1.x[kk - 1];
		yky = atoms_1.y[kay - 1] - atoms_1.y[kk - 1];
		zky = atoms_1.z__[kay - 1] - atoms_1.z__[kk - 1];
		vxx = -xr * (ftm2[0] + ftm2i[0]) + xix * frcxi[0] + xiy * 
			frcyi[0] + xiz * frczi[0] + xkx * frcxk[0] + xky * 
			frcyk[0] + xkz * frczk[0];
		vyx = -yr * (ftm2[0] + ftm2i[0]) + yix * frcxi[0] + yiy * 
			frcyi[0] + yiz * frczi[0] + ykx * frcxk[0] + yky * 
			frcyk[0] + ykz * frczk[0];
		vzx = -zr * (ftm2[0] + ftm2i[0]) + zix * frcxi[0] + ziy * 
			frcyi[0] + ziz * frczi[0] + zkx * frcxk[0] + zky * 
			frcyk[0] + zkz * frczk[0];
		vyy = -yr * (ftm2[1] + ftm2i[1]) + yix * frcxi[1] + yiy * 
			frcyi[1] + yiz * frczi[1] + ykx * frcxk[1] + yky * 
			frcyk[1] + ykz * frczk[1];
		vzy = -zr * (ftm2[1] + ftm2i[1]) + zix * frcxi[1] + ziy * 
			frcyi[1] + ziz * frczi[1] + zkx * frcxk[1] + zky * 
			frcyk[1] + zkz * frczk[1];
		vzz = -zr * (ftm2[2] + ftm2i[2]) + zix * frcxi[2] + ziy * 
			frcyi[2] + ziz * frczi[2] + zkx * frcxk[2] + zky * 
			frcyk[2] + zkz * frczk[2];
		vir_ref(1, 1) = vir_ref(1, 1) + vxx;
		vir_ref(2, 1) = vir_ref(2, 1) + vyx;
		vir_ref(3, 1) = vir_ref(3, 1) + vzx;
		vir_ref(1, 2) = vir_ref(1, 2) + vyx;
		vir_ref(2, 2) = vir_ref(2, 2) + vyy;
		vir_ref(3, 2) = vir_ref(3, 2) + vzy;
		vir_ref(1, 3) = vir_ref(1, 3) + vzx;
		vir_ref(2, 3) = vir_ref(2, 3) + vzy;
		vir_ref(3, 3) = vir_ref(3, 3) + vzz;
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

/*     for periodic boundary conditions with large cutoffs */
/*     neighbors must be found by the replicates method */

    if (bound_1.use_replica__) {

/*     calculate interaction with other unit cells */

	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
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
	    usei = usage_1.use[ii - 1] || usage_1.use[iz - 1] || usage_1.use[
		    ix - 1];

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
	    for (k = i__; k <= i__2; ++k) {
		kk = mpole_1.ipole[k - 1];
		kz = mpole_1.zaxis[k - 1];
		kx = mpole_1.xaxis[k - 1];
		usek = usage_1.use[kk - 1] || usage_1.use[kz - 1] || 
			usage_1.use[kx - 1];
		if (group_1.use_group__) {
		    groups_(&proceed, &fgrp, &ii, &kk, &c__0, &c__0, &c__0, &
			    c__0);
		}
		proceed = TRUE_;
		if (proceed) {
		    proceed = usei || usek;
		}
		if (! proceed) {
		    goto L20;
		}
		i__3 = cell_1.ncell;
		for (jcell = 1; jcell <= i__3; ++jcell) {
		    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
		    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
		    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
		    imager_(&xr, &yr, &zr, &jcell);
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
			scale3i = scale3;
			scale5i = scale5;
			scale7i = scale7;
			dsc3 = scale3;
			dsc5 = scale5;
			dsc7 = scale7;
			psc3 = scale3;
			psc5 = scale5;
			psc7 = scale7;
			if (bound_1.use_polymer__) {
			    if (r2 <= bound_1.polycut2) {
				scale3i *= uscale[kk - 1];
				scale5i *= uscale[kk - 1];
				scale7i *= uscale[kk - 1];
				dsc3 *= dscale[kk - 1];
				dsc5 *= dscale[kk - 1];
				dsc7 *= dscale[kk - 1];
				psc3 *= pscale[kk - 1];
				psc5 *= pscale[kk - 1];
				psc7 *= pscale[kk - 1];
			    }
			}

/*     construct necessary auxiliary vectors */

			dixdk[0] = di[1] * dk[2] - di[2] * dk[1];
			dixdk[1] = di[2] * dk[0] - di[0] * dk[2];
			dixdk[2] = di[0] * dk[1] - di[1] * dk[0];
			dixuk[0] = di[1] * uind_ref(3, k) - di[2] * uind_ref(
				2, k);
			dixuk[1] = di[2] * uind_ref(1, k) - di[0] * uind_ref(
				3, k);
			dixuk[2] = di[0] * uind_ref(2, k) - di[1] * uind_ref(
				1, k);
			dkxui[0] = dk[1] * uind_ref(3, i__) - dk[2] * 
				uind_ref(2, i__);
			dkxui[1] = dk[2] * uind_ref(1, i__) - dk[0] * 
				uind_ref(3, i__);
			dkxui[2] = dk[0] * uind_ref(2, i__) - dk[1] * 
				uind_ref(1, i__);
			dixukp[0] = di[1] * uinp_ref(3, k) - di[2] * uinp_ref(
				2, k);
			dixukp[1] = di[2] * uinp_ref(1, k) - di[0] * uinp_ref(
				3, k);
			dixukp[2] = di[0] * uinp_ref(2, k) - di[1] * uinp_ref(
				1, k);
			dkxuip[0] = dk[1] * uinp_ref(3, i__) - dk[2] * 
				uinp_ref(2, i__);
			dkxuip[1] = dk[2] * uinp_ref(1, i__) - dk[0] * 
				uinp_ref(3, i__);
			dkxuip[2] = dk[0] * uinp_ref(2, i__) - dk[1] * 
				uinp_ref(1, i__);
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
			qiqkr[0] = qi[0] * qkr[0] + qi[3] * qkr[1] + qi[6] * 
				qkr[2];
			qiqkr[1] = qi[1] * qkr[0] + qi[4] * qkr[1] + qi[7] * 
				qkr[2];
			qiqkr[2] = qi[2] * qkr[0] + qi[5] * qkr[1] + qi[8] * 
				qkr[2];
			qkqir[0] = qk[0] * qir[0] + qk[3] * qir[1] + qk[6] * 
				qir[2];
			qkqir[1] = qk[1] * qir[0] + qk[4] * qir[1] + qk[7] * 
				qir[2];
			qkqir[2] = qk[2] * qir[0] + qk[5] * qir[1] + qk[8] * 
				qir[2];
			qixqk[0] = qi[1] * qk[2] + qi[4] * qk[5] + qi[7] * qk[
				8] - qi[2] * qk[1] - qi[5] * qk[4] - qi[8] * 
				qk[7];
			qixqk[1] = qi[2] * qk[0] + qi[5] * qk[3] + qi[8] * qk[
				6] - qi[0] * qk[2] - qi[3] * qk[5] - qi[6] * 
				qk[8];
			qixqk[2] = qi[0] * qk[1] + qi[3] * qk[4] + qi[6] * qk[
				7] - qi[1] * qk[0] - qi[4] * qk[3] - qi[7] * 
				qk[6];
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
			qidk[0] = qi[0] * dk[0] + qi[3] * dk[1] + qi[6] * dk[
				2];
			qidk[1] = qi[1] * dk[0] + qi[4] * dk[1] + qi[7] * dk[
				2];
			qidk[2] = qi[2] * dk[0] + qi[5] * dk[1] + qi[8] * dk[
				2];
			qkdi[0] = qk[0] * di[0] + qk[3] * di[1] + qk[6] * di[
				2];
			qkdi[1] = qk[1] * di[0] + qk[4] * di[1] + qk[7] * di[
				2];
			qkdi[2] = qk[2] * di[0] + qk[5] * di[1] + qk[8] * di[
				2];
			qiuk[0] = qi[0] * uind_ref(1, k) + qi[3] * uind_ref(2,
				 k) + qi[6] * uind_ref(3, k);
			qiuk[1] = qi[1] * uind_ref(1, k) + qi[4] * uind_ref(2,
				 k) + qi[7] * uind_ref(3, k);
			qiuk[2] = qi[2] * uind_ref(1, k) + qi[5] * uind_ref(2,
				 k) + qi[8] * uind_ref(3, k);
			qkui[0] = qk[0] * uind_ref(1, i__) + qk[3] * uind_ref(
				2, i__) + qk[6] * uind_ref(3, i__);
			qkui[1] = qk[1] * uind_ref(1, i__) + qk[4] * uind_ref(
				2, i__) + qk[7] * uind_ref(3, i__);
			qkui[2] = qk[2] * uind_ref(1, i__) + qk[5] * uind_ref(
				2, i__) + qk[8] * uind_ref(3, i__);
			qiukp[0] = qi[0] * uinp_ref(1, k) + qi[3] * uinp_ref(
				2, k) + qi[6] * uinp_ref(3, k);
			qiukp[1] = qi[1] * uinp_ref(1, k) + qi[4] * uinp_ref(
				2, k) + qi[7] * uinp_ref(3, k);
			qiukp[2] = qi[2] * uinp_ref(1, k) + qi[5] * uinp_ref(
				2, k) + qi[8] * uinp_ref(3, k);
			qkuip[0] = qk[0] * uinp_ref(1, i__) + qk[3] * 
				uinp_ref(2, i__) + qk[6] * uinp_ref(3, i__);
			qkuip[1] = qk[1] * uinp_ref(1, i__) + qk[4] * 
				uinp_ref(2, i__) + qk[7] * uinp_ref(3, i__);
			qkuip[2] = qk[2] * uinp_ref(1, i__) + qk[5] * 
				uinp_ref(2, i__) + qk[8] * uinp_ref(3, i__);
			dixqkr[0] = di[1] * qkr[2] - di[2] * qkr[1];
			dixqkr[1] = di[2] * qkr[0] - di[0] * qkr[2];
			dixqkr[2] = di[0] * qkr[1] - di[1] * qkr[0];
			dkxqir[0] = dk[1] * qir[2] - dk[2] * qir[1];
			dkxqir[1] = dk[2] * qir[0] - dk[0] * qir[2];
			dkxqir[2] = dk[0] * qir[1] - dk[1] * qir[0];
			uixqkr[0] = uind_ref(2, i__) * qkr[2] - uind_ref(3, 
				i__) * qkr[1];
			uixqkr[1] = uind_ref(3, i__) * qkr[0] - uind_ref(1, 
				i__) * qkr[2];
			uixqkr[2] = uind_ref(1, i__) * qkr[1] - uind_ref(2, 
				i__) * qkr[0];
			ukxqir[0] = uind_ref(2, k) * qir[2] - uind_ref(3, k) *
				 qir[1];
			ukxqir[1] = uind_ref(3, k) * qir[0] - uind_ref(1, k) *
				 qir[2];
			ukxqir[2] = uind_ref(1, k) * qir[1] - uind_ref(2, k) *
				 qir[0];
			uixqkrp[0] = uinp_ref(2, i__) * qkr[2] - uinp_ref(3, 
				i__) * qkr[1];
			uixqkrp[1] = uinp_ref(3, i__) * qkr[0] - uinp_ref(1, 
				i__) * qkr[2];
			uixqkrp[2] = uinp_ref(1, i__) * qkr[1] - uinp_ref(2, 
				i__) * qkr[0];
			ukxqirp[0] = uinp_ref(2, k) * qir[2] - uinp_ref(3, k) 
				* qir[1];
			ukxqirp[1] = uinp_ref(3, k) * qir[0] - uinp_ref(1, k) 
				* qir[2];
			ukxqirp[2] = uinp_ref(1, k) * qir[1] - uinp_ref(2, k) 
				* qir[0];
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
			sc[6] = qir[0] * dk[0] + qir[1] * dk[1] + qir[2] * dk[
				2];
			sc[7] = qkr[0] * di[0] + qkr[1] * di[1] + qkr[2] * di[
				2];
			sc[8] = qir[0] * qkr[0] + qir[1] * qkr[1] + qir[2] * 
				qkr[2];
			sc[9] = qi[0] * qk[0] + qi[1] * qk[1] + qi[2] * qk[2] 
				+ qi[3] * qk[3] + qi[4] * qk[4] + qi[5] * qk[
				5] + qi[6] * qk[6] + qi[7] * qk[7] + qi[8] * 
				qk[8];

/*     calculate scalar products for induced components */

			sci[0] = uind_ref(1, i__) * dk[0] + uind_ref(2, i__) *
				 dk[1] + uind_ref(3, i__) * dk[2] + di[0] * 
				uind_ref(1, k) + di[1] * uind_ref(2, k) + di[
				2] * uind_ref(3, k);
			sci[1] = uind_ref(1, i__) * uind_ref(1, k) + uind_ref(
				2, i__) * uind_ref(2, k) + uind_ref(3, i__) * 
				uind_ref(3, k);
			sci[2] = uind_ref(1, i__) * xr + uind_ref(2, i__) * 
				yr + uind_ref(3, i__) * zr;
			sci[3] = uind_ref(1, k) * xr + uind_ref(2, k) * yr + 
				uind_ref(3, k) * zr;
			sci[6] = qir[0] * uind_ref(1, k) + qir[1] * uind_ref(
				2, k) + qir[2] * uind_ref(3, k);
			sci[7] = qkr[0] * uind_ref(1, i__) + qkr[1] * 
				uind_ref(2, i__) + qkr[2] * uind_ref(3, i__);
			scip[0] = uinp_ref(1, i__) * dk[0] + uinp_ref(2, i__) 
				* dk[1] + uinp_ref(3, i__) * dk[2] + di[0] * 
				uinp_ref(1, k) + di[1] * uinp_ref(2, k) + di[
				2] * uinp_ref(3, k);
			scip[1] = uind_ref(1, i__) * uinp_ref(1, k) + 
				uind_ref(2, i__) * uinp_ref(2, k) + uind_ref(
				3, i__) * uinp_ref(3, k) + uinp_ref(1, i__) * 
				uind_ref(1, k) + uinp_ref(2, i__) * uind_ref(
				2, k) + uinp_ref(3, i__) * uind_ref(3, k);
			scip[2] = uinp_ref(1, i__) * xr + uinp_ref(2, i__) * 
				yr + uinp_ref(3, i__) * zr;
			scip[3] = uinp_ref(1, k) * xr + uinp_ref(2, k) * yr + 
				uinp_ref(3, k) * zr;
			scip[6] = qir[0] * uinp_ref(1, k) + qir[1] * uinp_ref(
				2, k) + qir[2] * uinp_ref(3, k);
			scip[7] = qkr[0] * uinp_ref(1, i__) + qkr[1] * 
				uinp_ref(2, i__) + qkr[2] * uinp_ref(3, i__);

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

/*     compute the energy contributions for this interaction */

			e = rr1 * gl[0] + rr3 * (gl[1] + gl[6]) + rr5 * (gl[2]
				 + gl[7] + gl[8]) + rr7 * (gl[3] + gl[5]) + 
				rr9 * gl[4];
			ei = (rr3 * (gli[0] + gli[5]) * psc3 + rr5 * (gli[1] 
				+ gli[6]) * psc5 + rr7 * gli[2] * psc7) * .5;
			e = f * e;
			ei = f * ei;
			if (bound_1.use_polymer__) {
			    if (r2 <= bound_1.polycut2) {
				e *= mscale[kk - 1];
			    }
			}
			if (group_1.use_group__) {
			    e *= fgrp;
/*                 ei = ei * fgrp */
			}
			if (ii == kk) {
			    e *= .5;
			    ei *= .5;
			}
			energi_1.em += e;
			energi_1.ep += ei;

/*     increment the total intermolecular energy */

			if (molcul_1.molcule[ii - 1] != molcul_1.molcule[kk - 
				1]) {
			    inter_1.einter = inter_1.einter + e + ei;
			}

/*     intermediate variables for the permanent components */

			gf[0] = rr3 * gl[0] + rr5 * (gl[1] + gl[6]) + rr7 * (
				gl[2] + gl[7] + gl[8]) + rr9 * (gl[3] + gl[5])
				 + rr11 * gl[4];
			gf[1] = -ck * rr3 + sc[3] * rr5 - sc[5] * rr7;
			gf[2] = ci * rr3 + sc[2] * rr5 + sc[4] * rr7;
			gf[3] = rr5 * 2.;
			gf[4] = (-ck * rr5 + sc[3] * rr7 - sc[5] * rr9) * 2.;
			gf[5] = (-ci * rr5 - sc[2] * rr7 - sc[4] * rr9) * 2.;
			gf[6] = rr7 * 4.;

/*     intermediate variables for the induced components */

			gfi[0] = rr5 * .5 * ((gli[0] + gli[5]) * psc3 + (glip[
				0] + glip[5]) * dsc3 + scip[1] * scale3i) + 
				rr7 * .5 * ((gli[6] + gli[1]) * psc5 + (glip[
				6] + glip[1]) * dsc5 - (sci[2] * scip[3] + 
				scip[2] * sci[3]) * scale5i) + rr9 * .5 * (
				gli[2] * psc7 + glip[2] * dsc7);
			gfi[1] = -rr3 * ck + rr5 * sc[3] - rr7 * sc[5];
			gfi[2] = rr3 * ci + rr5 * sc[2] + rr7 * sc[4];
			gfi[3] = rr5 * 2.;
			gfi[4] = rr7 * (sci[3] * psc7 + scip[3] * dsc7);
			gfi[5] = -rr7 * (sci[2] * psc7 + scip[2] * dsc7);

/*     get the permanent force components */

			ftm2[0] = gf[0] * xr + gf[1] * di[0] + gf[2] * dk[0] 
				+ gf[3] * (qkdi[0] - qidk[0]) + gf[4] * qir[0]
				 + gf[5] * qkr[0] + gf[6] * (qiqkr[0] + qkqir[
				0]);
			ftm2[1] = gf[0] * yr + gf[1] * di[1] + gf[2] * dk[1] 
				+ gf[3] * (qkdi[1] - qidk[1]) + gf[4] * qir[1]
				 + gf[5] * qkr[1] + gf[6] * (qiqkr[1] + qkqir[
				1]);
			ftm2[2] = gf[0] * zr + gf[1] * di[2] + gf[2] * dk[2] 
				+ gf[3] * (qkdi[2] - qidk[2]) + gf[4] * qir[2]
				 + gf[5] * qkr[2] + gf[6] * (qiqkr[2] + qkqir[
				2]);

/*     get the induced force components */

			ftm2i[0] = gfi[0] * xr + (-rr3 * ck * (uind_ref(1, 
				i__) * psc3 + uinp_ref(1, i__) * dsc3) + rr5 *
				 sc[3] * (uind_ref(1, i__) * psc5 + uinp_ref(
				1, i__) * dsc5) - rr7 * sc[5] * (uind_ref(1, 
				i__) * psc7 + uinp_ref(1, i__) * dsc7)) * .5 
				+ (rr3 * ci * (uind_ref(1, k) * psc3 + 
				uinp_ref(1, k) * dsc3) + rr5 * sc[2] * (
				uind_ref(1, k) * psc5 + uinp_ref(1, k) * dsc5)
				 + rr7 * sc[4] * (uind_ref(1, k) * psc7 + 
				uinp_ref(1, k) * dsc7)) * .5 + rr5 * scale5i *
				 (sci[3] * uinp_ref(1, i__) + scip[3] * 
				uind_ref(1, i__) + sci[2] * uinp_ref(1, k) + 
				scip[2] * uind_ref(1, k)) * .5 + (sci[3] * 
				psc5 + scip[3] * dsc5) * .5 * rr5 * di[0] + (
				sci[2] * psc5 + scip[2] * dsc5) * .5 * rr5 * 
				dk[0] + gfi[3] * .5 * ((qkui[0] - qiuk[0]) * 
				psc5 + (qkuip[0] - qiukp[0]) * dsc5) + gfi[4] 
				* qir[0] + gfi[5] * qkr[0];
			ftm2i[1] = gfi[0] * yr + (-rr3 * ck * (uind_ref(2, 
				i__) * psc3 + uinp_ref(2, i__) * dsc3) + rr5 *
				 sc[3] * (uind_ref(2, i__) * psc5 + uinp_ref(
				2, i__) * dsc5) - rr7 * sc[5] * (uind_ref(2, 
				i__) * psc7 + uinp_ref(2, i__) * dsc7)) * .5 
				+ (rr3 * ci * (uind_ref(2, k) * psc3 + 
				uinp_ref(2, k) * dsc3) + rr5 * sc[2] * (
				uind_ref(2, k) * psc5 + uinp_ref(2, k) * dsc5)
				 + rr7 * sc[4] * (uind_ref(2, k) * psc7 + 
				uinp_ref(2, k) * dsc7)) * .5 + rr5 * scale5i *
				 (sci[3] * uinp_ref(2, i__) + scip[3] * 
				uind_ref(2, i__) + sci[2] * uinp_ref(2, k) + 
				scip[2] * uind_ref(2, k)) * .5 + (sci[3] * 
				psc5 + scip[3] * dsc5) * .5 * rr5 * di[1] + (
				sci[2] * psc5 + scip[2] * dsc5) * .5 * rr5 * 
				dk[1] + gfi[3] * .5 * ((qkui[1] - qiuk[1]) * 
				psc5 + (qkuip[1] - qiukp[1]) * dsc5) + gfi[4] 
				* qir[1] + gfi[5] * qkr[1];
			ftm2i[2] = gfi[0] * zr + (-rr3 * ck * (uind_ref(3, 
				i__) * psc3 + uinp_ref(3, i__) * dsc3) + rr5 *
				 sc[3] * (uind_ref(3, i__) * psc5 + uinp_ref(
				3, i__) * dsc5) - rr7 * sc[5] * (uind_ref(3, 
				i__) * psc7 + uinp_ref(3, i__) * dsc7)) * .5 
				+ (rr3 * ci * (uind_ref(3, k) * psc3 + 
				uinp_ref(3, k) * dsc3) + rr5 * sc[2] * (
				uind_ref(3, k) * psc5 + uinp_ref(3, k) * dsc5)
				 + rr7 * sc[4] * (uind_ref(3, k) * psc7 + 
				uinp_ref(3, k) * dsc7)) * .5 + rr5 * scale5i *
				 (sci[3] * uinp_ref(3, i__) + scip[3] * 
				uind_ref(3, i__) + sci[2] * uinp_ref(3, k) + 
				scip[2] * uind_ref(3, k)) * .5 + (sci[3] * 
				psc5 + scip[3] * dsc5) * .5 * rr5 * di[2] + (
				sci[2] * psc5 + scip[2] * dsc5) * .5 * rr5 * 
				dk[2] + gfi[3] * .5 * ((qkui[2] - qiuk[2]) * 
				psc5 + (qkuip[2] - qiukp[2]) * dsc5) + gfi[4] 
				* qir[2] + gfi[5] * qkr[2];

/*     account for partially excluded induced interactions */

			temp3 = rr3 * .5 * ((gli[0] + gli[5]) * pscale[kk - 1]
				 + (glip[0] + glip[5]) * dscale[kk - 1]);
			temp5 = rr5 * .5 * ((gli[1] + gli[6]) * pscale[kk - 1]
				 + (glip[1] + glip[6]) * dscale[kk - 1]);
			temp7 = rr7 * .5 * (gli[2] * pscale[kk - 1] + glip[2] 
				* dscale[kk - 1]);
			fridmp[0] = temp3 * ddsc3[0] + temp5 * ddsc5[0] + 
				temp7 * ddsc7[0];
			fridmp[1] = temp3 * ddsc3[1] + temp5 * ddsc5[1] + 
				temp7 * ddsc7[1];
			fridmp[2] = temp3 * ddsc3[2] + temp5 * ddsc5[2] + 
				temp7 * ddsc7[2];

/*     find some scaling terms for induced-induced force */

			temp3 = rr3 * .5 * uscale[kk - 1] * scip[1];
			temp5 = rr5 * -.5 * uscale[kk - 1] * (sci[2] * scip[3]
				 + scip[2] * sci[3]);
			findmp[0] = temp3 * ddsc3[0] + temp5 * ddsc5[0];
			findmp[1] = temp3 * ddsc3[1] + temp5 * ddsc5[1];
			findmp[2] = temp3 * ddsc3[2] + temp5 * ddsc5[2];

/*     modify induced force for partially excluded interactions */

			ftm2i[0] = ftm2i[0] - fridmp[0] - findmp[0];
			ftm2i[1] = ftm2i[1] - fridmp[1] - findmp[1];
			ftm2i[2] = ftm2i[2] - fridmp[2] - findmp[2];

/*     correction to convert mutual to direct polarization force */

			if (s_cmp(polpot_1.poltyp, "DIRECT", (ftnlen)6, (
				ftnlen)6) == 0) {
			    gfd = (rr5 * scip[1] * scale3i - rr7 * (scip[2] * 
				    sci[3] + sci[2] * scip[3]) * scale5i) * 
				    .5;
			    temp5 = rr5 * .5 * scale5i;
			    fdir[0] = gfd * xr + temp5 * (sci[3] * uinp_ref(1,
				     i__) + scip[3] * uind_ref(1, i__) + sci[
				    2] * uinp_ref(1, k) + scip[2] * uind_ref(
				    1, k));
			    fdir[1] = gfd * yr + temp5 * (sci[3] * uinp_ref(2,
				     i__) + scip[3] * uind_ref(2, i__) + sci[
				    2] * uinp_ref(2, k) + scip[2] * uind_ref(
				    2, k));
			    fdir[2] = gfd * zr + temp5 * (sci[3] * uinp_ref(3,
				     i__) + scip[3] * uind_ref(3, i__) + sci[
				    2] * uinp_ref(3, k) + scip[2] * uind_ref(
				    3, k));
			    ftm2i[0] = ftm2i[0] - fdir[0] + findmp[0];
			    ftm2i[1] = ftm2i[1] - fdir[1] + findmp[1];
			    ftm2i[2] = ftm2i[2] - fdir[2] + findmp[2];
			}

/*     intermediate terms for induced torque on multipoles */

			gti[1] = rr5 * .5 * (sci[3] * psc5 + scip[3] * dsc5);
			gti[2] = rr5 * .5 * (sci[2] * psc5 + scip[2] * dsc5);
			gti[3] = gfi[3];
			gti[4] = gfi[4];
			gti[5] = gfi[5];

/*     get the permanent torque components */

			ttm2[0] = -rr3 * dixdk[0] + gf[1] * dixr[0] - gf[4] * 
				rxqir[0] + gf[3] * (dixqkr[0] + dkxqir[0] + 
				rxqidk[0] - qixqk[0] * 2.) - gf[6] * (rxqikr[
				0] + qkrxqir[0]);
			ttm2[1] = -rr3 * dixdk[1] + gf[1] * dixr[1] - gf[4] * 
				rxqir[1] + gf[3] * (dixqkr[1] + dkxqir[1] + 
				rxqidk[1] - qixqk[1] * 2.) - gf[6] * (rxqikr[
				1] + qkrxqir[1]);
			ttm2[2] = -rr3 * dixdk[2] + gf[1] * dixr[2] - gf[4] * 
				rxqir[2] + gf[3] * (dixqkr[2] + dkxqir[2] + 
				rxqidk[2] - qixqk[2] * 2.) - gf[6] * (rxqikr[
				2] + qkrxqir[2]);
			ttm3[0] = rr3 * dixdk[0] + gf[2] * dkxr[0] - gf[5] * 
				rxqkr[0] - gf[3] * (dixqkr[0] + dkxqir[0] + 
				rxqkdi[0] - qixqk[0] * 2.) - gf[6] * (rxqkir[
				0] - qkrxqir[0]);
			ttm3[1] = rr3 * dixdk[1] + gf[2] * dkxr[1] - gf[5] * 
				rxqkr[1] - gf[3] * (dixqkr[1] + dkxqir[1] + 
				rxqkdi[1] - qixqk[1] * 2.) - gf[6] * (rxqkir[
				1] - qkrxqir[1]);
			ttm3[2] = rr3 * dixdk[2] + gf[2] * dkxr[2] - gf[5] * 
				rxqkr[2] - gf[3] * (dixqkr[2] + dkxqir[2] + 
				rxqkdi[2] - qixqk[2] * 2.) - gf[6] * (rxqkir[
				2] - qkrxqir[2]);

/*     get the induced torque components */

			ttm2i[0] = -rr3 * (dixuk[0] * psc3 + dixukp[0] * dsc3)
				 * .5 + gti[1] * dixr[0] + gti[3] * ((ukxqir[
				0] + rxqiuk[0]) * psc5 + (ukxqirp[0] + 
				rxqiukp[0]) * dsc5) * .5 - gti[4] * rxqir[0];
			ttm2i[1] = -rr3 * (dixuk[1] * psc3 + dixukp[1] * dsc3)
				 * .5 + gti[1] * dixr[1] + gti[3] * ((ukxqir[
				1] + rxqiuk[1]) * psc5 + (ukxqirp[1] + 
				rxqiukp[1]) * dsc5) * .5 - gti[4] * rxqir[1];
			ttm2i[2] = -rr3 * (dixuk[2] * psc3 + dixukp[2] * dsc3)
				 * .5 + gti[1] * dixr[2] + gti[3] * ((ukxqir[
				2] + rxqiuk[2]) * psc5 + (ukxqirp[2] + 
				rxqiukp[2]) * dsc5) * .5 - gti[4] * rxqir[2];
			ttm3i[0] = -rr3 * (dkxui[0] * psc3 + dkxuip[0] * dsc3)
				 * .5 + gti[2] * dkxr[0] - gti[3] * ((uixqkr[
				0] + rxqkui[0]) * psc5 + (uixqkrp[0] + 
				rxqkuip[0]) * dsc5) * .5 - gti[5] * rxqkr[0];
			ttm3i[1] = -rr3 * (dkxui[1] * psc3 + dkxuip[1] * dsc3)
				 * .5 + gti[2] * dkxr[1] - gti[3] * ((uixqkr[
				1] + rxqkui[1]) * psc5 + (uixqkrp[1] + 
				rxqkuip[1]) * dsc5) * .5 - gti[5] * rxqkr[1];
			ttm3i[2] = -rr3 * (dkxui[2] * psc3 + dkxuip[2] * dsc3)
				 * .5 + gti[2] * dkxr[2] - gti[3] * ((uixqkr[
				2] + rxqkui[2]) * psc5 + (uixqkrp[2] + 
				rxqkuip[2]) * dsc5) * .5 - gti[5] * rxqkr[2];

/*     handle the case where scaling is used */

			for (j = 1; j <= 3; ++j) {
			    ftm2[j - 1] = f * ftm2[j - 1];
			    ftm2i[j - 1] = f * ftm2i[j - 1];
			    ttm2[j - 1] = f * ttm2[j - 1];
			    ttm2i[j - 1] = f * ttm2i[j - 1];
			    ttm3[j - 1] = f * ttm3[j - 1];
			    ttm3i[j - 1] = f * ttm3i[j - 1];
			}
			if (bound_1.use_polymer__) {
			    if (r2 <= bound_1.polycut2) {
				for (j = 1; j <= 3; ++j) {
				    ftm2[j - 1] *= mscale[kk - 1];
				    ttm2[j - 1] *= mscale[kk - 1];
				    ttm3[j - 1] *= mscale[kk - 1];
				}
			    }
			}
			if (group_1.use_group__) {
			    for (j = 1; j <= 3; ++j) {
				ftm2[j - 1] *= fgrp;
				ttm2[j - 1] *= fgrp;
				ttm3[j - 1] *= fgrp;
/*                    ftm2i(j) = ftm2i(j) * fgrp */
/*                    ttm2i(j) = ttm2i(j) * fgrp */
/*                    ttm3i(j) = ttm3i(j) * fgrp */
			    }
			}
			if (ii == kk) {
			    for (j = 1; j <= 3; ++j) {
				ftm2[j - 1] *= .5;
				ftm2i[j - 1] *= .5;
				ttm2[j - 1] *= .5;
				ttm2i[j - 1] *= .5;
				ttm3[j - 1] *= .5;
				ttm3i[j - 1] *= .5;
			    }
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

/*     increment the internal virial tensor components */

			iaz = mpole_1.zaxis[i__ - 1];
			iax = mpole_1.xaxis[i__ - 1];
			iay = mpole_1.yaxis[i__ - 1];
			kaz = mpole_1.zaxis[k - 1];
			kax = mpole_1.xaxis[k - 1];
			kay = mpole_1.yaxis[k - 1];
			if (iaz == 0) {
			    iaz = ii;
			}
			if (iax == 0) {
			    iax = ii;
			}
			if (iay == 0) {
			    iay = ii;
			}
			if (kaz == 0) {
			    kaz = kk;
			}
			if (kax == 0) {
			    kax = kk;
			}
			if (kay == 0) {
			    kay = kk;
			}
			xiz = atoms_1.x[iaz - 1] - atoms_1.x[ii - 1];
			yiz = atoms_1.y[iaz - 1] - atoms_1.y[ii - 1];
			ziz = atoms_1.z__[iaz - 1] - atoms_1.z__[ii - 1];
			xix = atoms_1.x[iax - 1] - atoms_1.x[ii - 1];
			yix = atoms_1.y[iax - 1] - atoms_1.y[ii - 1];
			zix = atoms_1.z__[iax - 1] - atoms_1.z__[ii - 1];
			xiy = atoms_1.x[iay - 1] - atoms_1.x[ii - 1];
			yiy = atoms_1.y[iay - 1] - atoms_1.y[ii - 1];
			ziy = atoms_1.z__[iay - 1] - atoms_1.z__[ii - 1];
			xkz = atoms_1.x[kaz - 1] - atoms_1.x[kk - 1];
			ykz = atoms_1.y[kaz - 1] - atoms_1.y[kk - 1];
			zkz = atoms_1.z__[kaz - 1] - atoms_1.z__[kk - 1];
			xkx = atoms_1.x[kax - 1] - atoms_1.x[kk - 1];
			ykx = atoms_1.y[kax - 1] - atoms_1.y[kk - 1];
			zkx = atoms_1.z__[kax - 1] - atoms_1.z__[kk - 1];
			xky = atoms_1.x[kay - 1] - atoms_1.x[kk - 1];
			yky = atoms_1.y[kay - 1] - atoms_1.y[kk - 1];
			zky = atoms_1.z__[kay - 1] - atoms_1.z__[kk - 1];
			vxx = -xr * (ftm2[0] + ftm2i[0]) + xix * frcxi[0] + 
				xiy * frcyi[0] + xiz * frczi[0] + xkx * frcxk[
				0] + xky * frcyk[0] + xkz * frczk[0];
			vyx = -yr * (ftm2[0] + ftm2i[0]) + yix * frcxi[0] + 
				yiy * frcyi[0] + yiz * frczi[0] + ykx * frcxk[
				0] + yky * frcyk[0] + ykz * frczk[0];
			vzx = -zr * (ftm2[0] + ftm2i[0]) + zix * frcxi[0] + 
				ziy * frcyi[0] + ziz * frczi[0] + zkx * frcxk[
				0] + zky * frcyk[0] + zkz * frczk[0];
			vyy = -yr * (ftm2[1] + ftm2i[1]) + yix * frcxi[1] + 
				yiy * frcyi[1] + yiz * frczi[1] + ykx * frcxk[
				1] + yky * frcyk[1] + ykz * frczk[1];
			vzy = -zr * (ftm2[1] + ftm2i[1]) + zix * frcxi[1] + 
				ziy * frcyi[1] + ziz * frczi[1] + zkx * frcxk[
				1] + zky * frcyk[1] + zkz * frczk[1];
			vzz = -zr * (ftm2[2] + ftm2i[2]) + zix * frcxi[2] + 
				ziy * frcyi[2] + ziz * frczi[2] + zkx * frcxk[
				2] + zky * frcyk[2] + zkz * frczk[2];
			vir_ref(1, 1) = vir_ref(1, 1) + vxx;
			vir_ref(2, 1) = vir_ref(2, 1) + vyx;
			vir_ref(3, 1) = vir_ref(3, 1) + vzx;
			vir_ref(1, 2) = vir_ref(1, 2) + vyx;
			vir_ref(2, 2) = vir_ref(2, 2) + vyy;
			vir_ref(3, 2) = vir_ref(3, 2) + vzy;
			vir_ref(1, 3) = vir_ref(1, 3) + vzx;
			vir_ref(2, 3) = vir_ref(2, 3) + vzy;
			vir_ref(3, 3) = vir_ref(3, 3) + vzz;
		    }
		}
L20:
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
    }
    return 0;
} /* empole1a_ */

#undef rpole_ref
#undef uinp_ref
#undef uind_ref
#undef vir_ref
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




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine empole1b  --  neighbor list multipole derivs  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "empole1b" calculates the multipole and dipole polarization */
/*     energy and derivatives with respect to Cartesian coordinates */
/*     using a neighbor list */


/* Subroutine */ int empole1b_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__, r2;
    static integer ii, kk;
    static doublereal ei, ci;
    static integer ix;
    static doublereal di[3];
    static integer iz, kx, kz;
    static doublereal qi[9], ck, dk[3], qk[9], xr, yr, zr, gl[9], sc[10], gf[
	    7], rr1, rr3, rr5, rr7, rr9, gfd, gfi[6];
    static integer kkk, iax, iay, iaz, kax, kay, kaz;
    static doublereal pdi, pti, xix, yix, zix, xiy, yiy, ziy, xiz, yiz, ziz, 
	    xkx, ykx, zkx, xky, yky, zky, xkz, ykz, zkz, rr11, dsc3, dsc5, 
	    vxx, dsc7, vyy, vzz, vyx, vzx, vzy, qir[3], qkr[3], psc3, ftm2[3],
	     psc5, gli[7], psc7, sci[8], gti[6], ttm2[3], ttm3[3], damp, fdir[
	    3], qidk[3], qkdi[3], glip[7], fgrp, qiuk[3], qkui[3], dixr[3], 
	    dkxr[3], scip[8];
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
    static doublereal uscale[25000], findmp[3], fridmp[3], scale3i, scale5i, 
	    scale7i, dixukp[3], dkxuip[3], uixqkr[3], ukxqir[3], rxqiuk[3], 
	    rxqkui[3], dixqkr[3], dkxqir[3], rxqikr[3], rxqkir[3], rxqidk[3], 
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
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define elst_ref(a_1,a_2) neigh_1.elst[(a_2)*1200 + a_1 - 1201]
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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  inter.i  --  sum of intermolecular energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     einter   total intermolecular potential energy */




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




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  neigh.i  --  pairwise neighbor list indices and storage  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     lbuffer     width of the neighbor list buffer region */
/*     lbuf2       motion squared needed to trigger list rebuild */
/*     vbuf2       square of vdw cutoff plus neighbor list buffer */
/*     cbuf2       square of charge cutoff plus neighbor list buffer */
/*     mbuf2       square of multipole cutoff plus neighbor list buffer */
/*     nvlst       number of sites in list for each vdw site */
/*     vlst        site numbers in neighbor list of each vdw site */
/*     nelst       number of sites in list for each electrostatic site */
/*     elst        site numbers in list of each electrostatic site */
/*     dovlst      logical flag to rebuild vdw neighbor list */
/*     doclst      logical flag to rebuild charge neighbor list */
/*     domlst      logical flag to rebuild multipole neighbor list */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  virial.i  --  components of internal virial tensor  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     vir    total internal virial Cartesian tensor components */




/*     zero out multipole and polarization energy and derivatives */

    energi_1.em = 0.;
    energi_1.ep = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    dem_ref(j, i__) = 0.;
	    dep_ref(j, i__) = 0.;
	}
    }

/*     check the sign of multipole components at chiral sites */

    chkpole_();

/*     rotate the multipole components into the global frame */

    rotpole_();

/*     compute the induced dipoles at each polarizable atom */

    induce_();

/*     set arrays needed to scale connected atom interactions */

    if (mpole_1.npole == 0) {
	return 0;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mscale[i__ - 1] = 1.;
	pscale[i__ - 1] = 1.;
	dscale[i__ - 1] = 1.;
	uscale[i__ - 1] = 1.;
    }

/*     set conversion factor, cutoff and scaling coefficients */

    f = chgpot_1.electric / chgpot_1.dielec;
    switch_("MPOLE", (ftnlen)5);

/*     set scale factors for permanent multipole and induced terms */

    i__1 = mpole_1.npole - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
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
	i__2 = neigh_1.nelst[i__ - 1];
	for (kkk = 1; kkk <= i__2; ++kkk) {
	    k = elst_ref(kkk, i__);
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

/*     compute the energy contributions for this interaction */

		e = rr1 * gl[0] + rr3 * (gl[1] + gl[6]) + rr5 * (gl[2] + gl[7]
			 + gl[8]) + rr7 * (gl[3] + gl[5]) + rr9 * gl[4];
		ei = (rr3 * (gli[0] + gli[5]) * psc3 + rr5 * (gli[1] + gli[6])
			 * psc5 + rr7 * gli[2] * psc7) * .5;
		e = f * mscale[kk - 1] * e;
		ei = f * ei;
		energi_1.em += e;
		energi_1.ep += ei;

/*     increment the total intermolecular energy */

		if (molcul_1.molcule[ii - 1] != molcul_1.molcule[kk - 1]) {
		    inter_1.einter = inter_1.einter + e + ei;
		}

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

/*     increment the internal virial tensor components */

		iaz = mpole_1.zaxis[i__ - 1];
		iax = mpole_1.xaxis[i__ - 1];
		iay = mpole_1.yaxis[i__ - 1];
		kaz = mpole_1.zaxis[k - 1];
		kax = mpole_1.xaxis[k - 1];
		kay = mpole_1.yaxis[k - 1];
		if (iaz == 0) {
		    iaz = ii;
		}
		if (iax == 0) {
		    iax = ii;
		}
		if (iay == 0) {
		    iay = ii;
		}
		if (kaz == 0) {
		    kaz = kk;
		}
		if (kax == 0) {
		    kax = kk;
		}
		if (kay == 0) {
		    kay = kk;
		}
		xiz = atoms_1.x[iaz - 1] - atoms_1.x[ii - 1];
		yiz = atoms_1.y[iaz - 1] - atoms_1.y[ii - 1];
		ziz = atoms_1.z__[iaz - 1] - atoms_1.z__[ii - 1];
		xix = atoms_1.x[iax - 1] - atoms_1.x[ii - 1];
		yix = atoms_1.y[iax - 1] - atoms_1.y[ii - 1];
		zix = atoms_1.z__[iax - 1] - atoms_1.z__[ii - 1];
		xiy = atoms_1.x[iay - 1] - atoms_1.x[ii - 1];
		yiy = atoms_1.y[iay - 1] - atoms_1.y[ii - 1];
		ziy = atoms_1.z__[iay - 1] - atoms_1.z__[ii - 1];
		xkz = atoms_1.x[kaz - 1] - atoms_1.x[kk - 1];
		ykz = atoms_1.y[kaz - 1] - atoms_1.y[kk - 1];
		zkz = atoms_1.z__[kaz - 1] - atoms_1.z__[kk - 1];
		xkx = atoms_1.x[kax - 1] - atoms_1.x[kk - 1];
		ykx = atoms_1.y[kax - 1] - atoms_1.y[kk - 1];
		zkx = atoms_1.z__[kax - 1] - atoms_1.z__[kk - 1];
		xky = atoms_1.x[kay - 1] - atoms_1.x[kk - 1];
		yky = atoms_1.y[kay - 1] - atoms_1.y[kk - 1];
		zky = atoms_1.z__[kay - 1] - atoms_1.z__[kk - 1];
		vxx = -xr * (ftm2[0] + ftm2i[0]) + xix * frcxi[0] + xiy * 
			frcyi[0] + xiz * frczi[0] + xkx * frcxk[0] + xky * 
			frcyk[0] + xkz * frczk[0];
		vyx = -yr * (ftm2[0] + ftm2i[0]) + yix * frcxi[0] + yiy * 
			frcyi[0] + yiz * frczi[0] + ykx * frcxk[0] + yky * 
			frcyk[0] + ykz * frczk[0];
		vzx = -zr * (ftm2[0] + ftm2i[0]) + zix * frcxi[0] + ziy * 
			frcyi[0] + ziz * frczi[0] + zkx * frcxk[0] + zky * 
			frcyk[0] + zkz * frczk[0];
		vyy = -yr * (ftm2[1] + ftm2i[1]) + yix * frcxi[1] + yiy * 
			frcyi[1] + yiz * frczi[1] + ykx * frcxk[1] + yky * 
			frcyk[1] + ykz * frczk[1];
		vzy = -zr * (ftm2[1] + ftm2i[1]) + zix * frcxi[1] + ziy * 
			frcyi[1] + ziz * frczi[1] + zkx * frcxk[1] + zky * 
			frcyk[1] + zkz * frczk[1];
		vzz = -zr * (ftm2[2] + ftm2i[2]) + zix * frcxi[2] + ziy * 
			frcyi[2] + ziz * frczi[2] + zkx * frcxk[2] + zky * 
			frcyk[2] + zkz * frczk[2];
		vir_ref(1, 1) = vir_ref(1, 1) + vxx;
		vir_ref(2, 1) = vir_ref(2, 1) + vyx;
		vir_ref(3, 1) = vir_ref(3, 1) + vzx;
		vir_ref(1, 2) = vir_ref(1, 2) + vyx;
		vir_ref(2, 2) = vir_ref(2, 2) + vyy;
		vir_ref(3, 2) = vir_ref(3, 2) + vzy;
		vir_ref(1, 3) = vir_ref(1, 3) + vzx;
		vir_ref(2, 3) = vir_ref(2, 3) + vzy;
		vir_ref(3, 3) = vir_ref(3, 3) + vzz;
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
    return 0;
} /* empole1b_ */

#undef rpole_ref
#undef uinp_ref
#undef uind_ref
#undef elst_ref
#undef vir_ref
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




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine empole1c  --  Ewald multipole derivs via loop  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "empole1c" calculates the multipole and dipole polarization */
/*     energy and derivatives with respect to Cartesian coordinates */
/*     using particle mesh Ewald summation and a double loop */


/* Subroutine */ int empole1c_(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j;
    static doublereal ci, ei;
    static integer ii;
    static doublereal xd, yd, zd, xq, yq, zq, xu, yu, zu, xv, yv, zv, cii, 
	    dii, qii, dix, diy, uii, diz, uix, uiy, uiz, trq[3], xup, yup, 
	    zup, frcx[3], frcy[3], frcz[3], term, trqi[3], qixx, qixy, qixz, 
	    qiyy, qiyz, qizz, fterm, vterm;
    extern /* Subroutine */ int induce_(void);
    static doublereal eintra;
    extern /* Subroutine */ int ereal1c_(doublereal *), torque_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal xdfield, ydfield, zdfield;
    extern /* Subroutine */ int chkpole_(void);
    static doublereal xufield, yufield, zufield;
    extern /* Subroutine */ int rotpole_(void), emrecip1_(void);


#define dem_ref(a_1,a_2) deriv_1.dem[(a_2)*3 + a_1 - 4]
#define dep_ref(a_1,a_2) deriv_1.dep[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  ewald.i  --  parameters and options for Ewald summation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     aewald     Ewald convergence coefficient value (Ang-1) */
/*     boundary   Ewald boundary condition; none, tinfoil or vacuum */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  inter.i  --  sum of intermolecular energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     einter   total intermolecular potential energy */




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

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  virial.i  --  components of internal virial tensor  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     vir    total internal virial Cartesian tensor components */




/*     zero out multipole and polarization energy and derivatives */

    energi_1.em = 0.;
    energi_1.ep = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    dem_ref(j, i__) = 0.;
	    dep_ref(j, i__) = 0.;
	}
    }

/*     set the energy unit conversion factor */

    f = chgpot_1.electric / chgpot_1.dielec;

/*     check the sign of multipole components at chiral sites */

    chkpole_();

/*     rotate the multipole components into the global frame */

    rotpole_();

/*     compute the induced dipole moment at each atom */

    induce_();

/*     compute the reciprocal space part of the Ewald summation */

    emrecip1_();

/*     compute the real space part of the Ewald summation */

    ereal1c_(&eintra);

/*     compute the Ewald self-energy term over all the atoms */

    term = ewald_1.aewald * 2. * ewald_1.aewald;
    fterm = -f * ewald_1.aewald / 1.772453850905516027;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ci = rpole_ref(1, i__);
	dix = rpole_ref(2, i__);
	diy = rpole_ref(3, i__);
	diz = rpole_ref(4, i__);
	qixx = rpole_ref(5, i__);
	qixy = rpole_ref(6, i__);
	qixz = rpole_ref(7, i__);
	qiyy = rpole_ref(9, i__);
	qiyz = rpole_ref(10, i__);
	qizz = rpole_ref(13, i__);
	uix = uind_ref(1, i__);
	uiy = uind_ref(2, i__);
	uiz = uind_ref(3, i__);
	cii = ci * ci;
	dii = dix * dix + diy * diy + diz * diz;
	qii = qixx * qixx + qiyy * qiyy + qizz * qizz + (qixy * qixy + qixz * 
		qixz + qiyz * qiyz) * 2.;
	uii = dix * uix + diy * uiy + diz * uiz;
	e = fterm * (cii + term * (dii / 3. + term * 2. * qii / 5.));
	ei = fterm * term * uii / 3.;
	energi_1.em += e;
	energi_1.ep += ei;
    }

/*     compute the self-energy torque term due to induced dipole */

    trq[0] = 0.;
    trq[1] = 0.;
    trq[2] = 0.;
/* Computing 3rd power */
    d__1 = ewald_1.aewald;
    term = f * 1.3333333333333333 * (d__1 * (d__1 * d__1)) / 
	    1.772453850905516027;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dix = rpole_ref(2, i__);
	diy = rpole_ref(3, i__);
	diz = rpole_ref(4, i__);
	uix = (uind_ref(1, i__) + uinp_ref(1, i__)) * .5;
	uiy = (uind_ref(2, i__) + uinp_ref(2, i__)) * .5;
	uiz = (uind_ref(3, i__) + uinp_ref(3, i__)) * .5;
	trqi[0] = term * (diy * uiz - diz * uiy);
	trqi[1] = term * (diz * uix - dix * uiz);
	trqi[2] = term * (dix * uiy - diy * uix);
	torque_(&i__, trq, trqi, frcx, frcy, frcz);
    }

/*     compute the cell dipole boundary correction term */

    if (s_cmp(ewald_1.boundary, "VACUUM", (ftnlen)7, (ftnlen)6) == 0) {
	xd = 0.;
	yd = 0.;
	zd = 0.;
	xu = 0.;
	yu = 0.;
	zu = 0.;
	xup = 0.;
	yup = 0.;
	zup = 0.;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    xd = xd + rpole_ref(2, i__) + rpole_ref(1, i__) * atoms_1.x[ii - 
		    1];
	    yd = yd + rpole_ref(3, i__) + rpole_ref(1, i__) * atoms_1.y[ii - 
		    1];
	    zd = zd + rpole_ref(4, i__) + rpole_ref(1, i__) * atoms_1.z__[ii 
		    - 1];
	    xu += uind_ref(1, i__);
	    yu += uind_ref(2, i__);
	    zu += uind_ref(3, i__);
	    xup += uinp_ref(1, i__);
	    yup += uinp_ref(2, i__);
	    zup += uinp_ref(3, i__);
	}
	term = f * .66666666666666663 * (3.141592653589793238 / 
		boxes_1.volbox);
	energi_1.em += term * (xd * xd + yd * yd + zd * zd);
	energi_1.ep += term * (xd * xu + yd * yu + zd * zu);
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    dem_ref(1, ii) = dem_ref(1, ii) + term * 2. * rpole_ref(1, i__) * 
		    xd;
	    dem_ref(2, ii) = dem_ref(2, ii) + term * 2. * rpole_ref(1, i__) * 
		    yd;
	    dem_ref(3, ii) = dem_ref(3, ii) + term * 2. * rpole_ref(1, i__) * 
		    zd;
	    dep_ref(1, ii) = dep_ref(1, ii) + term * rpole_ref(1, i__) * (xu 
		    + xup);
	    dep_ref(2, ii) = dep_ref(2, ii) + term * rpole_ref(1, i__) * (yu 
		    + yup);
	    dep_ref(3, ii) = dep_ref(3, ii) + term * rpole_ref(1, i__) * (zu 
		    + zup);
	}
	xdfield = term * -2. * xd;
	ydfield = term * -2. * yd;
	zdfield = term * -2. * zd;
	xufield = -term * (xu + xup);
	yufield = -term * (yu + yup);
	zufield = -term * (zu + zup);
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    trq[0] = rpole_ref(3, i__) * zdfield - rpole_ref(4, i__) * 
		    ydfield;
	    trq[1] = rpole_ref(4, i__) * xdfield - rpole_ref(2, i__) * 
		    zdfield;
	    trq[2] = rpole_ref(2, i__) * ydfield - rpole_ref(3, i__) * 
		    xdfield;
	    trqi[0] = rpole_ref(3, i__) * zufield - rpole_ref(4, i__) * 
		    yufield;
	    trqi[1] = rpole_ref(4, i__) * xufield - rpole_ref(2, i__) * 
		    zufield;
	    trqi[2] = rpole_ref(2, i__) * yufield - rpole_ref(3, i__) * 
		    xufield;
	    torque_(&i__, trq, trqi, frcx, frcy, frcz);
	}

/*     boundary correction to virial due to overall cell dipole */

	xd = 0.;
	yd = 0.;
	zd = 0.;
	xq = 0.;
	yq = 0.;
	zq = 0.;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    xd += rpole_ref(2, i__);
	    yd += rpole_ref(3, i__);
	    zd += rpole_ref(4, i__);
	    xq += rpole_ref(1, i__) * atoms_1.x[ii - 1];
	    yq += rpole_ref(1, i__) * atoms_1.y[ii - 1];
	    zq += rpole_ref(1, i__) * atoms_1.z__[ii - 1];
	}
	xv = xq * (xd + (xu + xup) * .5);
	yv = yq * (yd + (yu + yup) * .5);
	zv = zq * (zd + (zu + zup) * .5);
	vterm = term * (xq * xq + yq * yq + zq * zq + (xv + yv + zv) * 2. + 
		xu * xup + yu * yup + zu * zup + xd * (xd + xu + xup) + yd * (
		yd + yu + yup) + zd * (zd + zu + zup));
	vir_ref(1, 1) = vir_ref(1, 1) + term * 2. * (xq * xq + xv) + vterm;
	vir_ref(2, 1) = vir_ref(2, 1) + term * 2. * (xq * yq + xv);
	vir_ref(3, 1) = vir_ref(3, 1) + term * 2. * (xq * zq + xv);
	vir_ref(1, 2) = vir_ref(1, 2) + term * 2. * (yq * xq + yv);
	vir_ref(2, 2) = vir_ref(2, 2) + term * 2. * (yq * yq + yv) + vterm;
	vir_ref(3, 2) = vir_ref(3, 2) + term * 2. * (yq * zq + yv);
	vir_ref(1, 3) = vir_ref(1, 3) + term * 2. * (zq * xq + zv);
	vir_ref(2, 3) = vir_ref(2, 3) + term * 2. * (zq * yq + zv);
	vir_ref(3, 3) = vir_ref(3, 3) + term * 2. * (zq * zq + zv) + vterm;
	if (s_cmp(polpot_1.poltyp, "DIRECT", (ftnlen)6, (ftnlen)6) == 0) {
	    vterm = term * (xu * xup + yu * yup + zu * zup);
	    vir_ref(1, 1) = vir_ref(1, 1) + vterm;
	    vir_ref(2, 2) = vir_ref(2, 2) + vterm;
	    vir_ref(3, 3) = vir_ref(3, 3) + vterm;
	}
    }

/*     intermolecular energy is total minus intramolecular part */

    inter_1.einter = inter_1.einter + energi_1.em + energi_1.ep - eintra;
    return 0;
} /* empole1c_ */

#undef rpole_ref
#undef uinp_ref
#undef uind_ref
#undef vir_ref
#undef dep_ref
#undef dem_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine ereal1c  --  ewald real space derivs via loop  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "ereal1c" evaluates the real space portion of the regular Ewald */
/*     summation energy and gradient due to atomic multipole interactions */
/*     and dipole polarizability */


/* Subroutine */ int ereal1c_(doublereal *eintra)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__, r2;
    static integer ii, kk;
    static doublereal ei, ci, di[3], qi[9], ck, dk[3], qk[9], bn[6], sc[10], 
	    gl[9], xr, yr, zr, gf[7], rr1, rr3, rr5, rr7, rr9, gfd, gfi[6];
    static integer iax, iay, iaz, kax, kay, kaz;
    static doublereal pdi, pti, rr11, xix, yix, zix, xiy, yiy, ziy, xiz, yiz, 
	    ziz, xkx, ykx, zkx, xky, yky, zky, xkz, ykz, zkz, erl, vxx, dsc3, 
	    dsc5, vyy, dsc7, vzz, vyx, vzx, vzy, qir[3], qkr[3], sci[8], psc3,
	     ftm2[3], psc5, gli[7], psc7, usc3, usc5, gfr[7], gti[6], ttm2[3],
	     ttm3[3], bfac;
    extern doublereal erfc_(doublereal *);
    static doublereal damp, gfdr, fdir[3], qidk[3], qkdi[3], erli, glip[7], 
	    scip[8], dixr[3], gfri[6], dkxr[3], qiuk[3], qkui[3], gtri[6], 
	    ddsc3[3], ddsc5[3], ddsc7[3], exp2a, ftm2i[3], alsq2, temp3, 
	    temp5, ftm2r[3], temp7, ttm2i[3], ttm3i[3], ttm2r[3], ttm3r[3];
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static integer jcell;
    static doublereal dixdk[3], frcxi[3], frcxk[3], frcyi[3], frcyk[3], frczi[
	    3], frczk[3], dkxui[3], dixuk[3], qiukp[3], qkuip[3], qiqkr[3], 
	    qkqir[3], qixqk[3], rxqir[3], scale3, rxqkr[3], scale5, scale7, 
	    alsq2n, ftm2ri[3], ttm2ri[3], ttm3ri[3], dscale[25000], pgamma, 
	    mscale[25000];
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal ralpha, pscale[25000], uscale[25000], findmp[3], fridmp[
	    3], dixukp[3], dkxuip[3], dixqkr[3], uixqkr[3], ukxqir[3], rxqiuk[
	    3], rxqkui[3], dkxqir[3], rxqikr[3], rxqkir[3], rxqidk[3], rxqkdi[
	    3];
    extern /* Subroutine */ int switch_(char *, ftnlen), torque_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal expdamp, qkrxqir[3], uixqkrp[3], ukxqirp[3], rxqiukp[3],
	     rxqkuip[3];


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
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cell.i  --  periodic boundaries using replicated cells  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xcell    length of the a-axis of the complete replicated cell */
/*     ycell    length of the b-axis of the complete replicated cell */
/*     zcell    length of the c-axis of the complete replicated cell */
/*     xcell2   half the length of the a-axis of the replicated cell */
/*     ycell2   half the length of the b-axis of the replicated cell */
/*     zcell2   half the length of the c-axis of the replicated cell */
/*     ncell    total number of cell replicates for periodic boundaries */
/*     icell    offset along axes for each replicate periodic cell */




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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  ewald.i  --  parameters and options for Ewald summation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     aewald     Ewald convergence coefficient value (Ang-1) */
/*     boundary   Ewald boundary condition; none, tinfoil or vacuum */




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

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  virial.i  --  components of internal virial tensor  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     vir    total internal virial Cartesian tensor components */




/*     zero out the intramolecular portion of the Ewald energy */

    *eintra = 0.;
    if (mpole_1.npole == 0) {
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
    switch_("EWALD", (ftnlen)5);

/*     set the permanent multipole and induced dipole values */

    i__1 = mpole_1.npole - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
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

/*     calculate the real space error function terms */

		ralpha = ewald_1.aewald * r__;
		bn[0] = erfc_(&ralpha) / r__;
/* Computing 2nd power */
		d__1 = ewald_1.aewald;
		alsq2 = d__1 * d__1 * 2.;
		alsq2n = 0.;
		if (ewald_1.aewald > 0.) {
		    alsq2n = 1. / (ewald_1.aewald * 1.772453850905516027);
		}
/* Computing 2nd power */
		d__1 = ralpha;
		exp2a = exp(-(d__1 * d__1));
		for (j = 1; j <= 5; ++j) {
		    bfac = (doublereal) ((j << 1) - 1);
		    alsq2n = alsq2 * alsq2n;
		    bn[j] = (bfac * bn[j - 1] + alsq2n * exp2a) / r2;
		}

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
		dsc3 = 1. - scale3 * dscale[kk - 1];
		dsc5 = 1. - scale5 * dscale[kk - 1];
		dsc7 = 1. - scale7 * dscale[kk - 1];
		psc3 = 1. - scale3 * pscale[kk - 1];
		psc5 = 1. - scale5 * pscale[kk - 1];
		psc7 = 1. - scale7 * pscale[kk - 1];
		usc3 = 1. - scale3 * uscale[kk - 1];
		usc5 = 1. - scale5 * uscale[kk - 1];

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

/*     calculate the scalar products for permanent components */

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

/*     calculate the scalar products for induced components */

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

/*     compute the energy contributions for this interaction */

		e = bn[0] * gl[0] + bn[1] * (gl[1] + gl[6]) + bn[2] * (gl[2] 
			+ gl[7] + gl[8]) + bn[3] * (gl[3] + gl[5]) + bn[4] * 
			gl[4];
		ei = (bn[1] * (gli[0] + gli[5]) + bn[2] * (gli[1] + gli[6]) + 
			bn[3] * gli[2]) * .5;

/*     get the real energy without any screening function */

		erl = rr1 * gl[0] + rr3 * (gl[1] + gl[6]) + rr5 * (gl[2] + gl[
			7] + gl[8]) + rr7 * (gl[3] + gl[5]) + rr9 * gl[4];
		erli = (rr3 * (gli[0] + gli[5]) * psc3 + rr5 * (gli[1] + gli[
			6]) * psc5 + rr7 * gli[2] * psc7) * .5;
		e -= (1. - mscale[kk - 1]) * erl;
		ei -= erli;
		e = f * e;
		ei = f * ei;
		energi_1.em += e;
		energi_1.ep += ei;

/*     increment the total intramolecular energy; assumes */
/*     intramolecular distances are less than half of cell */
/*     length and less than the ewald cutoff */

		if (molcul_1.molcule[ii - 1] == molcul_1.molcule[kk - 1]) {
		    *eintra += mscale[kk - 1] * erl * f;
		    *eintra += pscale[kk - 1] * .5 * (rr3 * (gli[0] + gli[5]) 
			    * scale3 + rr5 * (gli[1] + gli[6]) * scale5 + rr7 
			    * gli[2] * scale7);
		}

/*     intermediate variables for permanent force terms */

		gf[0] = bn[1] * gl[0] + bn[2] * (gl[1] + gl[6]) + bn[3] * (gl[
			2] + gl[7] + gl[8]) + bn[4] * (gl[3] + gl[5]) + bn[5] 
			* gl[4];
		gf[1] = -ck * bn[1] + sc[3] * bn[2] - sc[5] * bn[3];
		gf[2] = ci * bn[1] + sc[2] * bn[2] + sc[4] * bn[3];
		gf[3] = bn[2] * 2.;
		gf[4] = (-ck * bn[2] + sc[3] * bn[3] - sc[5] * bn[4]) * 2.;
		gf[5] = (-ci * bn[2] - sc[2] * bn[3] - sc[4] * bn[4]) * 2.;
		gf[6] = bn[3] * 4.;
		gfr[0] = rr3 * gl[0] + rr5 * (gl[1] + gl[6]) + rr7 * (gl[2] + 
			gl[7] + gl[8]) + rr9 * (gl[3] + gl[5]) + rr11 * gl[4];
		gfr[1] = -ck * rr3 + sc[3] * rr5 - sc[5] * rr7;
		gfr[2] = ci * rr3 + sc[2] * rr5 + sc[4] * rr7;
		gfr[3] = rr5 * 2.;
		gfr[4] = (-ck * rr5 + sc[3] * rr7 - sc[5] * rr9) * 2.;
		gfr[5] = (-ci * rr5 - sc[2] * rr7 - sc[4] * rr9) * 2.;
		gfr[6] = rr7 * 4.;

/*     intermediate variables for induced force terms */

		gfi[0] = bn[2] * .5 * (gli[0] + glip[0] + gli[5] + glip[5]) + 
			bn[2] * .5 * scip[1] + bn[3] * .5 * (gli[1] + glip[1] 
			+ gli[6] + glip[6]) - bn[3] * .5 * (sci[2] * scip[3] 
			+ scip[2] * sci[3]) + bn[4] * .5 * (gli[2] + glip[2]);
		gfi[1] = -ck * bn[1] + sc[3] * bn[2] - sc[5] * bn[3];
		gfi[2] = ci * bn[1] + sc[2] * bn[2] + sc[4] * bn[3];
		gfi[3] = bn[2] * 2.;
		gfi[4] = bn[3] * (sci[3] + scip[3]);
		gfi[5] = -bn[3] * (sci[2] + scip[2]);
		gfri[0] = rr5 * .5 * ((gli[0] + gli[5]) * psc3 + (glip[0] + 
			glip[5]) * dsc3 + scip[1] * usc3) + rr7 * .5 * ((gli[
			6] + gli[1]) * psc5 + (glip[6] + glip[1]) * dsc5 - (
			sci[2] * scip[3] + scip[2] * sci[3]) * usc5) + rr9 * 
			.5 * (gli[2] * psc7 + glip[2] * dsc7);
		gfri[1] = -rr3 * ck + rr5 * sc[3] - rr7 * sc[5];
		gfri[2] = rr3 * ci + rr5 * sc[2] + rr7 * sc[4];
		gfri[3] = rr5 * 2.;
		gfri[4] = rr7 * (sci[3] * psc7 + scip[3] * dsc7);
		gfri[5] = -rr7 * (sci[2] * psc7 + scip[2] * dsc7);

/*     get the permanent force with screening */

		ftm2[0] = gf[0] * xr + gf[1] * di[0] + gf[2] * dk[0] + gf[3] *
			 (qkdi[0] - qidk[0]) + gf[4] * qir[0] + gf[5] * qkr[0]
			 + gf[6] * (qiqkr[0] + qkqir[0]);
		ftm2[1] = gf[0] * yr + gf[1] * di[1] + gf[2] * dk[1] + gf[3] *
			 (qkdi[1] - qidk[1]) + gf[4] * qir[1] + gf[5] * qkr[1]
			 + gf[6] * (qiqkr[1] + qkqir[1]);
		ftm2[2] = gf[0] * zr + gf[1] * di[2] + gf[2] * dk[2] + gf[3] *
			 (qkdi[2] - qidk[2]) + gf[4] * qir[2] + gf[5] * qkr[2]
			 + gf[6] * (qiqkr[2] + qkqir[2]);

/*     get the permanent force without screening */

		ftm2r[0] = gfr[0] * xr + gfr[1] * di[0] + gfr[2] * dk[0] + 
			gfr[3] * (qkdi[0] - qidk[0]) + gfr[4] * qir[0] + gfr[
			5] * qkr[0] + gfr[6] * (qiqkr[0] + qkqir[0]);
		ftm2r[1] = gfr[0] * yr + gfr[1] * di[1] + gfr[2] * dk[1] + 
			gfr[3] * (qkdi[1] - qidk[1]) + gfr[4] * qir[1] + gfr[
			5] * qkr[1] + gfr[6] * (qiqkr[1] + qkqir[1]);
		ftm2r[2] = gfr[0] * zr + gfr[1] * di[2] + gfr[2] * dk[2] + 
			gfr[3] * (qkdi[2] - qidk[2]) + gfr[4] * qir[2] + gfr[
			5] * qkr[2] + gfr[6] * (qiqkr[2] + qkqir[2]);

/*     get the induced force with screening */

		ftm2i[0] = gfi[0] * xr + (gfi[1] * (uind_ref(1, i__) + 
			uinp_ref(1, i__)) + bn[2] * (sci[3] * uinp_ref(1, i__)
			 + scip[3] * uind_ref(1, i__)) + gfi[2] * (uind_ref(1,
			 k) + uinp_ref(1, k)) + bn[2] * (sci[2] * uinp_ref(1, 
			k) + scip[2] * uind_ref(1, k)) + (sci[3] + scip[3]) * 
			bn[2] * di[0] + (sci[2] + scip[2]) * bn[2] * dk[0] + 
			gfi[3] * (qkui[0] + qkuip[0] - qiuk[0] - qiukp[0])) * 
			.5 + gfi[4] * qir[0] + gfi[5] * qkr[0];
		ftm2i[1] = gfi[0] * yr + (gfi[1] * (uind_ref(2, i__) + 
			uinp_ref(2, i__)) + bn[2] * (sci[3] * uinp_ref(2, i__)
			 + scip[3] * uind_ref(2, i__)) + gfi[2] * (uind_ref(2,
			 k) + uinp_ref(2, k)) + bn[2] * (sci[2] * uinp_ref(2, 
			k) + scip[2] * uind_ref(2, k)) + (sci[3] + scip[3]) * 
			bn[2] * di[1] + (sci[2] + scip[2]) * bn[2] * dk[1] + 
			gfi[3] * (qkui[1] + qkuip[1] - qiuk[1] - qiukp[1])) * 
			.5 + gfi[4] * qir[1] + gfi[5] * qkr[1];
		ftm2i[2] = gfi[0] * zr + (gfi[1] * (uind_ref(3, i__) + 
			uinp_ref(3, i__)) + bn[2] * (sci[3] * uinp_ref(3, i__)
			 + scip[3] * uind_ref(3, i__)) + gfi[2] * (uind_ref(3,
			 k) + uinp_ref(3, k)) + bn[2] * (sci[2] * uinp_ref(3, 
			k) + scip[2] * uind_ref(3, k)) + (sci[3] + scip[3]) * 
			bn[2] * di[2] + (sci[2] + scip[2]) * bn[2] * dk[2] + 
			gfi[3] * (qkui[2] + qkuip[2] - qiuk[2] - qiukp[2])) * 
			.5 + gfi[4] * qir[2] + gfi[5] * qkr[2];

/*     get the induced force without screening */

		ftm2ri[0] = gfri[0] * xr + (-rr3 * ck * (uind_ref(1, i__) * 
			psc3 + uinp_ref(1, i__) * dsc3) + rr5 * sc[3] * (
			uind_ref(1, i__) * psc5 + uinp_ref(1, i__) * dsc5) - 
			rr7 * sc[5] * (uind_ref(1, i__) * psc7 + uinp_ref(1, 
			i__) * dsc7)) * .5 + (rr3 * ci * (uind_ref(1, k) * 
			psc3 + uinp_ref(1, k) * dsc3) + rr5 * sc[2] * (
			uind_ref(1, k) * psc5 + uinp_ref(1, k) * dsc5) + rr7 *
			 sc[4] * (uind_ref(1, k) * psc7 + uinp_ref(1, k) * 
			dsc7)) * .5 + rr5 * usc5 * (sci[3] * uinp_ref(1, i__) 
			+ scip[3] * uind_ref(1, i__) + sci[2] * uinp_ref(1, k)
			 + scip[2] * uind_ref(1, k)) * .5 + (sci[3] * psc5 + 
			scip[3] * dsc5) * .5 * rr5 * di[0] + (sci[2] * psc5 + 
			scip[2] * dsc5) * .5 * rr5 * dk[0] + gfri[3] * .5 * ((
			qkui[0] - qiuk[0]) * psc5 + (qkuip[0] - qiukp[0]) * 
			dsc5) + gfri[4] * qir[0] + gfri[5] * qkr[0];
		ftm2ri[1] = gfri[0] * yr + (-rr3 * ck * (uind_ref(2, i__) * 
			psc3 + uinp_ref(2, i__) * dsc3) + rr5 * sc[3] * (
			uind_ref(2, i__) * psc5 + uinp_ref(2, i__) * dsc5) - 
			rr7 * sc[5] * (uind_ref(2, i__) * psc7 + uinp_ref(2, 
			i__) * dsc7)) * .5 + (rr3 * ci * (uind_ref(2, k) * 
			psc3 + uinp_ref(2, k) * dsc3) + rr5 * sc[2] * (
			uind_ref(2, k) * psc5 + uinp_ref(2, k) * dsc5) + rr7 *
			 sc[4] * (uind_ref(2, k) * psc7 + uinp_ref(2, k) * 
			dsc7)) * .5 + rr5 * usc5 * (sci[3] * uinp_ref(2, i__) 
			+ scip[3] * uind_ref(2, i__) + sci[2] * uinp_ref(2, k)
			 + scip[2] * uind_ref(2, k)) * .5 + (sci[3] * psc5 + 
			scip[3] * dsc5) * .5 * rr5 * di[1] + (sci[2] * psc5 + 
			scip[2] * dsc5) * .5 * rr5 * dk[1] + gfri[3] * .5 * ((
			qkui[1] - qiuk[1]) * psc5 + (qkuip[1] - qiukp[1]) * 
			dsc5) + gfri[4] * qir[1] + gfri[5] * qkr[1];
		ftm2ri[2] = gfri[0] * zr + (-rr3 * ck * (uind_ref(3, i__) * 
			psc3 + uinp_ref(3, i__) * dsc3) + rr5 * sc[3] * (
			uind_ref(3, i__) * psc5 + uinp_ref(3, i__) * dsc5) - 
			rr7 * sc[5] * (uind_ref(3, i__) * psc7 + uinp_ref(3, 
			i__) * dsc7)) * .5 + (rr3 * ci * (uind_ref(3, k) * 
			psc3 + uinp_ref(3, k) * dsc3) + rr5 * sc[2] * (
			uind_ref(3, k) * psc5 + uinp_ref(3, k) * dsc5) + rr7 *
			 sc[4] * (uind_ref(3, k) * psc7 + uinp_ref(3, k) * 
			dsc7)) * .5 + rr5 * usc5 * (sci[3] * uinp_ref(3, i__) 
			+ scip[3] * uind_ref(3, i__) + sci[2] * uinp_ref(3, k)
			 + scip[2] * uind_ref(3, k)) * .5 + (sci[3] * psc5 + 
			scip[3] * dsc5) * .5 * rr5 * di[2] + (sci[2] * psc5 + 
			scip[2] * dsc5) * .5 * rr5 * dk[2] + gfri[3] * .5 * ((
			qkui[2] - qiuk[2]) * psc5 + (qkuip[2] - qiukp[2]) * 
			dsc5) + gfri[4] * qir[2] + gfri[5] * qkr[2];

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

/*     modify the forces for partially excluded interactions */

		ftm2i[0] = ftm2i[0] - fridmp[0] - findmp[0];
		ftm2i[1] = ftm2i[1] - fridmp[1] - findmp[1];
		ftm2i[2] = ftm2i[2] - fridmp[2] - findmp[2];

/*     correction to convert mutual to direct polarization force */

		if (s_cmp(polpot_1.poltyp, "DIRECT", (ftnlen)6, (ftnlen)6) == 
			0) {
		    gfd = (bn[2] * scip[1] - bn[3] * (scip[2] * sci[3] + sci[
			    2] * scip[3])) * .5;
		    gfdr = (rr5 * scip[1] * usc3 - rr7 * (scip[2] * sci[3] + 
			    sci[2] * scip[3]) * usc5) * .5;
		    ftm2i[0] = ftm2i[0] - gfd * xr - bn[2] * .5 * (sci[3] * 
			    uinp_ref(1, i__) + scip[3] * uind_ref(1, i__) + 
			    sci[2] * uinp_ref(1, k) + scip[2] * uind_ref(1, k)
			    );
		    ftm2i[1] = ftm2i[1] - gfd * yr - bn[2] * .5 * (sci[3] * 
			    uinp_ref(2, i__) + scip[3] * uind_ref(2, i__) + 
			    sci[2] * uinp_ref(2, k) + scip[2] * uind_ref(2, k)
			    );
		    ftm2i[2] = ftm2i[2] - gfd * zr - bn[2] * .5 * (sci[3] * 
			    uinp_ref(3, i__) + scip[3] * uind_ref(3, i__) + 
			    sci[2] * uinp_ref(3, k) + scip[2] * uind_ref(3, k)
			    );
		    fdir[0] = gfdr * xr + usc5 * .5 * rr5 * (sci[3] * 
			    uinp_ref(1, i__) + scip[3] * uind_ref(1, i__) + 
			    sci[2] * uinp_ref(1, k) + scip[2] * uind_ref(1, k)
			    );
		    fdir[1] = gfdr * yr + usc5 * .5 * rr5 * (sci[3] * 
			    uinp_ref(2, i__) + scip[3] * uind_ref(2, i__) + 
			    sci[2] * uinp_ref(2, k) + scip[2] * uind_ref(2, k)
			    );
		    fdir[2] = gfdr * zr + usc5 * .5 * rr5 * (sci[3] * 
			    uinp_ref(3, i__) + scip[3] * uind_ref(3, i__) + 
			    sci[2] * uinp_ref(3, k) + scip[2] * uind_ref(3, k)
			    );
		    ftm2i[0] = ftm2i[0] + fdir[0] + findmp[0];
		    ftm2i[1] = ftm2i[1] + fdir[1] + findmp[1];
		    ftm2i[2] = ftm2i[2] + fdir[2] + findmp[2];
		}

/*     intermediate variables for induced torque terms */

		gti[1] = bn[2] * .5 * (sci[3] + scip[3]);
		gti[2] = bn[2] * .5 * (sci[2] + scip[2]);
		gti[3] = gfi[3];
		gti[4] = gfi[4];
		gti[5] = gfi[5];
		gtri[1] = rr5 * .5 * (sci[3] * psc5 + scip[3] * dsc5);
		gtri[2] = rr5 * .5 * (sci[2] * psc5 + scip[2] * dsc5);
		gtri[3] = gfri[3];
		gtri[4] = gfri[4];
		gtri[5] = gfri[5];

/*     get the permanent torque with screening */

		ttm2[0] = -bn[1] * dixdk[0] + gf[1] * dixr[0] + gf[3] * (
			dixqkr[0] + dkxqir[0] + rxqidk[0] - qixqk[0] * 2.) - 
			gf[4] * rxqir[0] - gf[6] * (rxqikr[0] + qkrxqir[0]);
		ttm2[1] = -bn[1] * dixdk[1] + gf[1] * dixr[1] + gf[3] * (
			dixqkr[1] + dkxqir[1] + rxqidk[1] - qixqk[1] * 2.) - 
			gf[4] * rxqir[1] - gf[6] * (rxqikr[1] + qkrxqir[1]);
		ttm2[2] = -bn[1] * dixdk[2] + gf[1] * dixr[2] + gf[3] * (
			dixqkr[2] + dkxqir[2] + rxqidk[2] - qixqk[2] * 2.) - 
			gf[4] * rxqir[2] - gf[6] * (rxqikr[2] + qkrxqir[2]);
		ttm3[0] = bn[1] * dixdk[0] + gf[2] * dkxr[0] - gf[3] * (
			dixqkr[0] + dkxqir[0] + rxqkdi[0] - qixqk[0] * 2.) - 
			gf[5] * rxqkr[0] - gf[6] * (rxqkir[0] - qkrxqir[0]);
		ttm3[1] = bn[1] * dixdk[1] + gf[2] * dkxr[1] - gf[3] * (
			dixqkr[1] + dkxqir[1] + rxqkdi[1] - qixqk[1] * 2.) - 
			gf[5] * rxqkr[1] - gf[6] * (rxqkir[1] - qkrxqir[1]);
		ttm3[2] = bn[1] * dixdk[2] + gf[2] * dkxr[2] - gf[3] * (
			dixqkr[2] + dkxqir[2] + rxqkdi[2] - qixqk[2] * 2.) - 
			gf[5] * rxqkr[2] - gf[6] * (rxqkir[2] - qkrxqir[2]);

/*     get the permanent torque without screening */

		ttm2r[0] = -rr3 * dixdk[0] + gfr[1] * dixr[0] - gfr[4] * 
			rxqir[0] + gfr[3] * (dixqkr[0] + dkxqir[0] + rxqidk[0]
			 - qixqk[0] * 2.) - gfr[6] * (rxqikr[0] + qkrxqir[0]);
		ttm2r[1] = -rr3 * dixdk[1] + gfr[1] * dixr[1] - gfr[4] * 
			rxqir[1] + gfr[3] * (dixqkr[1] + dkxqir[1] + rxqidk[1]
			 - qixqk[1] * 2.) - gfr[6] * (rxqikr[1] + qkrxqir[1]);
		ttm2r[2] = -rr3 * dixdk[2] + gfr[1] * dixr[2] - gfr[4] * 
			rxqir[2] + gfr[3] * (dixqkr[2] + dkxqir[2] + rxqidk[2]
			 - qixqk[2] * 2.) - gfr[6] * (rxqikr[2] + qkrxqir[2]);
		ttm3r[0] = rr3 * dixdk[0] + gfr[2] * dkxr[0] - gfr[5] * rxqkr[
			0] - gfr[3] * (dixqkr[0] + dkxqir[0] + rxqkdi[0] - 
			qixqk[0] * 2.) - gfr[6] * (rxqkir[0] - qkrxqir[0]);
		ttm3r[1] = rr3 * dixdk[1] + gfr[2] * dkxr[1] - gfr[5] * rxqkr[
			1] - gfr[3] * (dixqkr[1] + dkxqir[1] + rxqkdi[1] - 
			qixqk[1] * 2.) - gfr[6] * (rxqkir[1] - qkrxqir[1]);
		ttm3r[2] = rr3 * dixdk[2] + gfr[2] * dkxr[2] - gfr[5] * rxqkr[
			2] - gfr[3] * (dixqkr[2] + dkxqir[2] + rxqkdi[2] - 
			qixqk[2] * 2.) - gfr[6] * (rxqkir[2] - qkrxqir[2]);

/*     get the induced torque with screening */

		ttm2i[0] = -bn[1] * (dixuk[0] + dixukp[0]) * .5 + gti[1] * 
			dixr[0] + gti[3] * (ukxqir[0] + rxqiuk[0] + ukxqirp[0]
			 + rxqiukp[0]) * .5 - gti[4] * rxqir[0];
		ttm2i[1] = -bn[1] * (dixuk[1] + dixukp[1]) * .5 + gti[1] * 
			dixr[1] + gti[3] * (ukxqir[1] + rxqiuk[1] + ukxqirp[1]
			 + rxqiukp[1]) * .5 - gti[4] * rxqir[1];
		ttm2i[2] = -bn[1] * (dixuk[2] + dixukp[2]) * .5 + gti[1] * 
			dixr[2] + gti[3] * (ukxqir[2] + rxqiuk[2] + ukxqirp[2]
			 + rxqiukp[2]) * .5 - gti[4] * rxqir[2];
		ttm3i[0] = -bn[1] * (dkxui[0] + dkxuip[0]) * .5 + gti[2] * 
			dkxr[0] - gti[3] * (uixqkr[0] + rxqkui[0] + uixqkrp[0]
			 + rxqkuip[0]) * .5 - gti[5] * rxqkr[0];
		ttm3i[1] = -bn[1] * (dkxui[1] + dkxuip[1]) * .5 + gti[2] * 
			dkxr[1] - gti[3] * (uixqkr[1] + rxqkui[1] + uixqkrp[1]
			 + rxqkuip[1]) * .5 - gti[5] * rxqkr[1];
		ttm3i[2] = -bn[1] * (dkxui[2] + dkxuip[2]) * .5 + gti[2] * 
			dkxr[2] - gti[3] * (uixqkr[2] + rxqkui[2] + uixqkrp[2]
			 + rxqkuip[2]) * .5 - gti[5] * rxqkr[2];

/*     get the induced torque without screening */

		ttm2ri[0] = -rr3 * (dixuk[0] * psc3 + dixukp[0] * dsc3) * .5 
			+ gtri[1] * dixr[0] + gtri[3] * ((ukxqir[0] + rxqiuk[
			0]) * psc5 + (ukxqirp[0] + rxqiukp[0]) * dsc5) * .5 - 
			gtri[4] * rxqir[0];
		ttm2ri[1] = -rr3 * (dixuk[1] * psc3 + dixukp[1] * dsc3) * .5 
			+ gtri[1] * dixr[1] + gtri[3] * ((ukxqir[1] + rxqiuk[
			1]) * psc5 + (ukxqirp[1] + rxqiukp[1]) * dsc5) * .5 - 
			gtri[4] * rxqir[1];
		ttm2ri[2] = -rr3 * (dixuk[2] * psc3 + dixukp[2] * dsc3) * .5 
			+ gtri[1] * dixr[2] + gtri[3] * ((ukxqir[2] + rxqiuk[
			2]) * psc5 + (ukxqirp[2] + rxqiukp[2]) * dsc5) * .5 - 
			gtri[4] * rxqir[2];
		ttm3ri[0] = -rr3 * (dkxui[0] * psc3 + dkxuip[0] * dsc3) * .5 
			+ gtri[2] * dkxr[0] - gtri[3] * ((uixqkr[0] + rxqkui[
			0]) * psc5 + (uixqkrp[0] + rxqkuip[0]) * dsc5) * .5 - 
			gtri[5] * rxqkr[0];
		ttm3ri[1] = -rr3 * (dkxui[1] * psc3 + dkxuip[1] * dsc3) * .5 
			+ gtri[2] * dkxr[1] - gtri[3] * ((uixqkr[1] + rxqkui[
			1]) * psc5 + (uixqkrp[1] + rxqkuip[1]) * dsc5) * .5 - 
			gtri[5] * rxqkr[1];
		ttm3ri[2] = -rr3 * (dkxui[2] * psc3 + dkxuip[2] * dsc3) * .5 
			+ gtri[2] * dkxr[2] - gtri[3] * ((uixqkr[2] + rxqkui[
			2]) * psc5 + (uixqkrp[2] + rxqkuip[2]) * dsc5) * .5 - 
			gtri[5] * rxqkr[2];

/*     handle the case where scaling is used */

		for (j = 1; j <= 3; ++j) {
		    ftm2[j - 1] = f * (ftm2[j - 1] - (1. - mscale[kk - 1]) * 
			    ftm2r[j - 1]);
		    ftm2i[j - 1] = f * (ftm2i[j - 1] - ftm2ri[j - 1]);
		    ttm2[j - 1] = f * (ttm2[j - 1] - (1. - mscale[kk - 1]) * 
			    ttm2r[j - 1]);
		    ttm2i[j - 1] = f * (ttm2i[j - 1] - ttm2ri[j - 1]);
		    ttm3[j - 1] = f * (ttm3[j - 1] - (1. - mscale[kk - 1]) * 
			    ttm3r[j - 1]);
		    ttm3i[j - 1] = f * (ttm3i[j - 1] - ttm3ri[j - 1]);
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

/*     increment the internal virial tensor components */

		iaz = mpole_1.zaxis[i__ - 1];
		iax = mpole_1.xaxis[i__ - 1];
		iay = mpole_1.yaxis[i__ - 1];
		kaz = mpole_1.zaxis[k - 1];
		kax = mpole_1.xaxis[k - 1];
		kay = mpole_1.yaxis[k - 1];
		if (iaz == 0) {
		    iaz = ii;
		}
		if (iax == 0) {
		    iax = ii;
		}
		if (iay == 0) {
		    iay = ii;
		}
		if (kaz == 0) {
		    kaz = kk;
		}
		if (kax == 0) {
		    kax = kk;
		}
		if (kay == 0) {
		    kay = kk;
		}
		xiz = atoms_1.x[iaz - 1] - atoms_1.x[ii - 1];
		yiz = atoms_1.y[iaz - 1] - atoms_1.y[ii - 1];
		ziz = atoms_1.z__[iaz - 1] - atoms_1.z__[ii - 1];
		xix = atoms_1.x[iax - 1] - atoms_1.x[ii - 1];
		yix = atoms_1.y[iax - 1] - atoms_1.y[ii - 1];
		zix = atoms_1.z__[iax - 1] - atoms_1.z__[ii - 1];
		xiy = atoms_1.x[iay - 1] - atoms_1.x[ii - 1];
		yiy = atoms_1.y[iay - 1] - atoms_1.y[ii - 1];
		ziy = atoms_1.z__[iay - 1] - atoms_1.z__[ii - 1];
		xkz = atoms_1.x[kaz - 1] - atoms_1.x[kk - 1];
		ykz = atoms_1.y[kaz - 1] - atoms_1.y[kk - 1];
		zkz = atoms_1.z__[kaz - 1] - atoms_1.z__[kk - 1];
		xkx = atoms_1.x[kax - 1] - atoms_1.x[kk - 1];
		ykx = atoms_1.y[kax - 1] - atoms_1.y[kk - 1];
		zkx = atoms_1.z__[kax - 1] - atoms_1.z__[kk - 1];
		xky = atoms_1.x[kay - 1] - atoms_1.x[kk - 1];
		yky = atoms_1.y[kay - 1] - atoms_1.y[kk - 1];
		zky = atoms_1.z__[kay - 1] - atoms_1.z__[kk - 1];
		vxx = -xr * (ftm2[0] + ftm2i[0]) + xix * frcxi[0] + xiy * 
			frcyi[0] + xiz * frczi[0] + xkx * frcxk[0] + xky * 
			frcyk[0] + xkz * frczk[0];
		vyx = -yr * (ftm2[0] + ftm2i[0]) + yix * frcxi[0] + yiy * 
			frcyi[0] + yiz * frczi[0] + ykx * frcxk[0] + yky * 
			frcyk[0] + ykz * frczk[0];
		vzx = -zr * (ftm2[0] + ftm2i[0]) + zix * frcxi[0] + ziy * 
			frcyi[0] + ziz * frczi[0] + zkx * frcxk[0] + zky * 
			frcyk[0] + zkz * frczk[0];
		vyy = -yr * (ftm2[1] + ftm2i[1]) + yix * frcxi[1] + yiy * 
			frcyi[1] + yiz * frczi[1] + ykx * frcxk[1] + yky * 
			frcyk[1] + ykz * frczk[1];
		vzy = -zr * (ftm2[1] + ftm2i[1]) + zix * frcxi[1] + ziy * 
			frcyi[1] + ziz * frczi[1] + zkx * frcxk[1] + zky * 
			frcyk[1] + zkz * frczk[1];
		vzz = -zr * (ftm2[2] + ftm2i[2]) + zix * frcxi[2] + ziy * 
			frcyi[2] + ziz * frczi[2] + zkx * frcxk[2] + zky * 
			frcyk[2] + zkz * frczk[2];
		vir_ref(1, 1) = vir_ref(1, 1) + vxx;
		vir_ref(2, 1) = vir_ref(2, 1) + vyx;
		vir_ref(3, 1) = vir_ref(3, 1) + vzx;
		vir_ref(1, 2) = vir_ref(1, 2) + vyx;
		vir_ref(2, 2) = vir_ref(2, 2) + vyy;
		vir_ref(3, 2) = vir_ref(3, 2) + vzy;
		vir_ref(1, 3) = vir_ref(1, 3) + vzx;
		vir_ref(2, 3) = vir_ref(2, 3) + vzy;
		vir_ref(3, 3) = vir_ref(3, 3) + vzz;
	    }
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

/*     for periodic boundary conditions with large cutoffs */
/*     neighbors must be found by the replicates method */

    if (bound_1.use_replica__) {

/*     calculate interaction with other unit cells */

	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
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
	    for (k = i__; k <= i__2; ++k) {
		kk = mpole_1.ipole[k - 1];
		i__3 = cell_1.ncell;
		for (jcell = 1; jcell <= i__3; ++jcell) {
		    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
		    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
		    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
		    imager_(&xr, &yr, &zr, &jcell);
		    r2 = xr * xr + yr * yr + zr * zr;
		    if (! (bound_1.use_polymer__ && r2 <= bound_1.polycut2)) {
			mscale[kk - 1] = 1.;
			pscale[kk - 1] = 1.;
			dscale[kk - 1] = 1.;
			uscale[kk - 1] = 1.;
		    }
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

/*     calculate the real space error function terms */

			ralpha = ewald_1.aewald * r__;
			bn[0] = erfc_(&ralpha) / r__;
/* Computing 2nd power */
			d__1 = ewald_1.aewald;
			alsq2 = d__1 * d__1 * 2.;
			alsq2n = 0.;
			if (ewald_1.aewald > 0.) {
			    alsq2n = 1. / (ewald_1.aewald * 
				    1.772453850905516027);
			}
/* Computing 2nd power */
			d__1 = ralpha;
			exp2a = exp(-(d__1 * d__1));
			for (j = 1; j <= 5; ++j) {
			    bfac = (doublereal) ((j << 1) - 1);
			    alsq2n = alsq2 * alsq2n;
			    bn[j] = (bfac * bn[j - 1] + alsq2n * exp2a) / r2;
			}

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
			dsc3 = 1. - scale3 * dscale[kk - 1];
			dsc5 = 1. - scale5 * dscale[kk - 1];
			dsc7 = 1. - scale7 * dscale[kk - 1];
			psc3 = 1. - scale3 * pscale[kk - 1];
			psc5 = 1. - scale5 * pscale[kk - 1];
			psc7 = 1. - scale7 * pscale[kk - 1];
			usc3 = 1. - scale3 * uscale[kk - 1];
			usc5 = 1. - scale5 * uscale[kk - 1];

/*     construct necessary auxiliary vectors */

			dixdk[0] = di[1] * dk[2] - di[2] * dk[1];
			dixdk[1] = di[2] * dk[0] - di[0] * dk[2];
			dixdk[2] = di[0] * dk[1] - di[1] * dk[0];
			dixuk[0] = di[1] * uind_ref(3, k) - di[2] * uind_ref(
				2, k);
			dixuk[1] = di[2] * uind_ref(1, k) - di[0] * uind_ref(
				3, k);
			dixuk[2] = di[0] * uind_ref(2, k) - di[1] * uind_ref(
				1, k);
			dkxui[0] = dk[1] * uind_ref(3, i__) - dk[2] * 
				uind_ref(2, i__);
			dkxui[1] = dk[2] * uind_ref(1, i__) - dk[0] * 
				uind_ref(3, i__);
			dkxui[2] = dk[0] * uind_ref(2, i__) - dk[1] * 
				uind_ref(1, i__);
			dixukp[0] = di[1] * uinp_ref(3, k) - di[2] * uinp_ref(
				2, k);
			dixukp[1] = di[2] * uinp_ref(1, k) - di[0] * uinp_ref(
				3, k);
			dixukp[2] = di[0] * uinp_ref(2, k) - di[1] * uinp_ref(
				1, k);
			dkxuip[0] = dk[1] * uinp_ref(3, i__) - dk[2] * 
				uinp_ref(2, i__);
			dkxuip[1] = dk[2] * uinp_ref(1, i__) - dk[0] * 
				uinp_ref(3, i__);
			dkxuip[2] = dk[0] * uinp_ref(2, i__) - dk[1] * 
				uinp_ref(1, i__);
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
			qiqkr[0] = qi[0] * qkr[0] + qi[3] * qkr[1] + qi[6] * 
				qkr[2];
			qiqkr[1] = qi[1] * qkr[0] + qi[4] * qkr[1] + qi[7] * 
				qkr[2];
			qiqkr[2] = qi[2] * qkr[0] + qi[5] * qkr[1] + qi[8] * 
				qkr[2];
			qkqir[0] = qk[0] * qir[0] + qk[3] * qir[1] + qk[6] * 
				qir[2];
			qkqir[1] = qk[1] * qir[0] + qk[4] * qir[1] + qk[7] * 
				qir[2];
			qkqir[2] = qk[2] * qir[0] + qk[5] * qir[1] + qk[8] * 
				qir[2];
			qixqk[0] = qi[1] * qk[2] + qi[4] * qk[5] + qi[7] * qk[
				8] - qi[2] * qk[1] - qi[5] * qk[4] - qi[8] * 
				qk[7];
			qixqk[1] = qi[2] * qk[0] + qi[5] * qk[3] + qi[8] * qk[
				6] - qi[0] * qk[2] - qi[3] * qk[5] - qi[6] * 
				qk[8];
			qixqk[2] = qi[0] * qk[1] + qi[3] * qk[4] + qi[6] * qk[
				7] - qi[1] * qk[0] - qi[4] * qk[3] - qi[7] * 
				qk[6];
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
			qidk[0] = qi[0] * dk[0] + qi[3] * dk[1] + qi[6] * dk[
				2];
			qidk[1] = qi[1] * dk[0] + qi[4] * dk[1] + qi[7] * dk[
				2];
			qidk[2] = qi[2] * dk[0] + qi[5] * dk[1] + qi[8] * dk[
				2];
			qkdi[0] = qk[0] * di[0] + qk[3] * di[1] + qk[6] * di[
				2];
			qkdi[1] = qk[1] * di[0] + qk[4] * di[1] + qk[7] * di[
				2];
			qkdi[2] = qk[2] * di[0] + qk[5] * di[1] + qk[8] * di[
				2];
			qiuk[0] = qi[0] * uind_ref(1, k) + qi[3] * uind_ref(2,
				 k) + qi[6] * uind_ref(3, k);
			qiuk[1] = qi[1] * uind_ref(1, k) + qi[4] * uind_ref(2,
				 k) + qi[7] * uind_ref(3, k);
			qiuk[2] = qi[2] * uind_ref(1, k) + qi[5] * uind_ref(2,
				 k) + qi[8] * uind_ref(3, k);
			qkui[0] = qk[0] * uind_ref(1, i__) + qk[3] * uind_ref(
				2, i__) + qk[6] * uind_ref(3, i__);
			qkui[1] = qk[1] * uind_ref(1, i__) + qk[4] * uind_ref(
				2, i__) + qk[7] * uind_ref(3, i__);
			qkui[2] = qk[2] * uind_ref(1, i__) + qk[5] * uind_ref(
				2, i__) + qk[8] * uind_ref(3, i__);
			qiukp[0] = qi[0] * uinp_ref(1, k) + qi[3] * uinp_ref(
				2, k) + qi[6] * uinp_ref(3, k);
			qiukp[1] = qi[1] * uinp_ref(1, k) + qi[4] * uinp_ref(
				2, k) + qi[7] * uinp_ref(3, k);
			qiukp[2] = qi[2] * uinp_ref(1, k) + qi[5] * uinp_ref(
				2, k) + qi[8] * uinp_ref(3, k);
			qkuip[0] = qk[0] * uinp_ref(1, i__) + qk[3] * 
				uinp_ref(2, i__) + qk[6] * uinp_ref(3, i__);
			qkuip[1] = qk[1] * uinp_ref(1, i__) + qk[4] * 
				uinp_ref(2, i__) + qk[7] * uinp_ref(3, i__);
			qkuip[2] = qk[2] * uinp_ref(1, i__) + qk[5] * 
				uinp_ref(2, i__) + qk[8] * uinp_ref(3, i__);
			dixqkr[0] = di[1] * qkr[2] - di[2] * qkr[1];
			dixqkr[1] = di[2] * qkr[0] - di[0] * qkr[2];
			dixqkr[2] = di[0] * qkr[1] - di[1] * qkr[0];
			dkxqir[0] = dk[1] * qir[2] - dk[2] * qir[1];
			dkxqir[1] = dk[2] * qir[0] - dk[0] * qir[2];
			dkxqir[2] = dk[0] * qir[1] - dk[1] * qir[0];
			uixqkr[0] = uind_ref(2, i__) * qkr[2] - uind_ref(3, 
				i__) * qkr[1];
			uixqkr[1] = uind_ref(3, i__) * qkr[0] - uind_ref(1, 
				i__) * qkr[2];
			uixqkr[2] = uind_ref(1, i__) * qkr[1] - uind_ref(2, 
				i__) * qkr[0];
			ukxqir[0] = uind_ref(2, k) * qir[2] - uind_ref(3, k) *
				 qir[1];
			ukxqir[1] = uind_ref(3, k) * qir[0] - uind_ref(1, k) *
				 qir[2];
			ukxqir[2] = uind_ref(1, k) * qir[1] - uind_ref(2, k) *
				 qir[0];
			uixqkrp[0] = uinp_ref(2, i__) * qkr[2] - uinp_ref(3, 
				i__) * qkr[1];
			uixqkrp[1] = uinp_ref(3, i__) * qkr[0] - uinp_ref(1, 
				i__) * qkr[2];
			uixqkrp[2] = uinp_ref(1, i__) * qkr[1] - uinp_ref(2, 
				i__) * qkr[0];
			ukxqirp[0] = uinp_ref(2, k) * qir[2] - uinp_ref(3, k) 
				* qir[1];
			ukxqirp[1] = uinp_ref(3, k) * qir[0] - uinp_ref(1, k) 
				* qir[2];
			ukxqirp[2] = uinp_ref(1, k) * qir[1] - uinp_ref(2, k) 
				* qir[0];
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

/*     calculate the scalar products for permanent components */

			sc[1] = di[0] * dk[0] + di[1] * dk[1] + di[2] * dk[2];
			sc[2] = di[0] * xr + di[1] * yr + di[2] * zr;
			sc[3] = dk[0] * xr + dk[1] * yr + dk[2] * zr;
			sc[4] = qir[0] * xr + qir[1] * yr + qir[2] * zr;
			sc[5] = qkr[0] * xr + qkr[1] * yr + qkr[2] * zr;
			sc[6] = qir[0] * dk[0] + qir[1] * dk[1] + qir[2] * dk[
				2];
			sc[7] = qkr[0] * di[0] + qkr[1] * di[1] + qkr[2] * di[
				2];
			sc[8] = qir[0] * qkr[0] + qir[1] * qkr[1] + qir[2] * 
				qkr[2];
			sc[9] = qi[0] * qk[0] + qi[1] * qk[1] + qi[2] * qk[2] 
				+ qi[3] * qk[3] + qi[4] * qk[4] + qi[5] * qk[
				5] + qi[6] * qk[6] + qi[7] * qk[7] + qi[8] * 
				qk[8];

/*     calculate the scalar products for induced components */

			sci[0] = uind_ref(1, i__) * dk[0] + uind_ref(2, i__) *
				 dk[1] + uind_ref(3, i__) * dk[2] + di[0] * 
				uind_ref(1, k) + di[1] * uind_ref(2, k) + di[
				2] * uind_ref(3, k);
			sci[1] = uind_ref(1, i__) * uind_ref(1, k) + uind_ref(
				2, i__) * uind_ref(2, k) + uind_ref(3, i__) * 
				uind_ref(3, k);
			sci[2] = uind_ref(1, i__) * xr + uind_ref(2, i__) * 
				yr + uind_ref(3, i__) * zr;
			sci[3] = uind_ref(1, k) * xr + uind_ref(2, k) * yr + 
				uind_ref(3, k) * zr;
			sci[6] = qir[0] * uind_ref(1, k) + qir[1] * uind_ref(
				2, k) + qir[2] * uind_ref(3, k);
			sci[7] = qkr[0] * uind_ref(1, i__) + qkr[1] * 
				uind_ref(2, i__) + qkr[2] * uind_ref(3, i__);
			scip[0] = uinp_ref(1, i__) * dk[0] + uinp_ref(2, i__) 
				* dk[1] + uinp_ref(3, i__) * dk[2] + di[0] * 
				uinp_ref(1, k) + di[1] * uinp_ref(2, k) + di[
				2] * uinp_ref(3, k);
			scip[1] = uind_ref(1, i__) * uinp_ref(1, k) + 
				uind_ref(2, i__) * uinp_ref(2, k) + uind_ref(
				3, i__) * uinp_ref(3, k) + uinp_ref(1, i__) * 
				uind_ref(1, k) + uinp_ref(2, i__) * uind_ref(
				2, k) + uinp_ref(3, i__) * uind_ref(3, k);
			scip[2] = uinp_ref(1, i__) * xr + uinp_ref(2, i__) * 
				yr + uinp_ref(3, i__) * zr;
			scip[3] = uinp_ref(1, k) * xr + uinp_ref(2, k) * yr + 
				uinp_ref(3, k) * zr;
			scip[6] = qir[0] * uinp_ref(1, k) + qir[1] * uinp_ref(
				2, k) + qir[2] * uinp_ref(3, k);
			scip[7] = qkr[0] * uinp_ref(1, i__) + qkr[1] * 
				uinp_ref(2, i__) + qkr[2] * uinp_ref(3, i__);

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

/*     compute the energy contributions for this interaction */

			e = bn[0] * gl[0] + bn[1] * (gl[1] + gl[6]) + bn[2] * 
				(gl[2] + gl[7] + gl[8]) + bn[3] * (gl[3] + gl[
				5]) + bn[4] * gl[4];
			ei = (bn[1] * (gli[0] + gli[5]) + bn[2] * (gli[1] + 
				gli[6]) + bn[3] * gli[2]) * .5;

/*     get the real energy without any screening function */

			erl = rr1 * gl[0] + rr3 * (gl[1] + gl[6]) + rr5 * (gl[
				2] + gl[7] + gl[8]) + rr7 * (gl[3] + gl[5]) + 
				rr9 * gl[4];
			erli = (rr3 * (gli[0] + gli[5]) * psc3 + rr5 * (gli[1]
				 + gli[6]) * psc5 + rr7 * gli[2] * psc7) * .5;
			if (bound_1.use_polymer__ && r2 <= bound_1.polycut2) {
			    e -= (1. - mscale[kk - 1]) * erl;
			}
			ei -= erli;
			e = f * e;
			ei = f * ei;
			if (ii == kk) {
			    e *= .5;
			    ei *= .5;
			}
			energi_1.em += e;
			energi_1.ep += ei;

/*     increment the total intramolecular energy; assumes */
/*     intramolecular distances are less than half of cell */
/*     length and less than the ewald cutoff */

			if (molcul_1.molcule[ii - 1] == molcul_1.molcule[kk - 
				1]) {
			    *eintra += mscale[kk - 1] * erl * f;
			    *eintra += pscale[kk - 1] * .5 * (rr3 * (gli[0] + 
				    gli[5]) * scale3 + rr5 * (gli[1] + gli[6])
				     * scale5 + rr7 * gli[2] * scale7);
			}

/*     intermediate variables for permanent force terms */

			gf[0] = bn[1] * gl[0] + bn[2] * (gl[1] + gl[6]) + bn[
				3] * (gl[2] + gl[7] + gl[8]) + bn[4] * (gl[3] 
				+ gl[5]) + bn[5] * gl[4];
			gf[1] = -ck * bn[1] + sc[3] * bn[2] - sc[5] * bn[3];
			gf[2] = ci * bn[1] + sc[2] * bn[2] + sc[4] * bn[3];
			gf[3] = bn[2] * 2.;
			gf[4] = (-ck * bn[2] + sc[3] * bn[3] - sc[5] * bn[4]) 
				* 2.;
			gf[5] = (-ci * bn[2] - sc[2] * bn[3] - sc[4] * bn[4]) 
				* 2.;
			gf[6] = bn[3] * 4.;
			gfr[0] = rr3 * gl[0] + rr5 * (gl[1] + gl[6]) + rr7 * (
				gl[2] + gl[7] + gl[8]) + rr9 * (gl[3] + gl[5])
				 + rr11 * gl[4];
			gfr[1] = -ck * rr3 + sc[3] * rr5 - sc[5] * rr7;
			gfr[2] = ci * rr3 + sc[2] * rr5 + sc[4] * rr7;
			gfr[3] = rr5 * 2.;
			gfr[4] = (-ck * rr5 + sc[3] * rr7 - sc[5] * rr9) * 2.;
			gfr[5] = (-ci * rr5 - sc[2] * rr7 - sc[4] * rr9) * 2.;
			gfr[6] = rr7 * 4.;

/*     intermediate variables for induced force terms */

			gfi[0] = bn[2] * .5 * (gli[0] + glip[0] + gli[5] + 
				glip[5]) + bn[2] * .5 * scip[1] + bn[3] * .5 *
				 (gli[1] + glip[1] + gli[6] + glip[6]) - bn[3]
				 * .5 * (sci[2] * scip[3] + scip[2] * sci[3]) 
				+ bn[4] * .5 * (gli[2] + glip[2]);
			gfi[1] = -ck * bn[1] + sc[3] * bn[2] - sc[5] * bn[3];
			gfi[2] = ci * bn[1] + sc[2] * bn[2] + sc[4] * bn[3];
			gfi[3] = bn[2] * 2.;
			gfi[4] = bn[3] * (sci[3] + scip[3]);
			gfi[5] = -bn[3] * (sci[2] + scip[2]);
			gfri[0] = rr5 * .5 * ((gli[0] + gli[5]) * psc3 + (
				glip[0] + glip[5]) * dsc3 + scip[1] * usc3) + 
				rr7 * .5 * ((gli[6] + gli[1]) * psc5 + (glip[
				6] + glip[1]) * dsc5 - (sci[2] * scip[3] + 
				scip[2] * sci[3]) * usc5) + rr9 * .5 * (gli[2]
				 * psc7 + glip[2] * dsc7);
			gfri[1] = -rr3 * ck + rr5 * sc[3] - rr7 * sc[5];
			gfri[2] = rr3 * ci + rr5 * sc[2] + rr7 * sc[4];
			gfri[3] = rr5 * 2.;
			gfri[4] = rr7 * (sci[3] * psc7 + scip[3] * dsc7);
			gfri[5] = -rr7 * (sci[2] * psc7 + scip[2] * dsc7);

/*     get the permanent force with screening */

			ftm2[0] = gf[0] * xr + gf[1] * di[0] + gf[2] * dk[0] 
				+ gf[3] * (qkdi[0] - qidk[0]) + gf[4] * qir[0]
				 + gf[5] * qkr[0] + gf[6] * (qiqkr[0] + qkqir[
				0]);
			ftm2[1] = gf[0] * yr + gf[1] * di[1] + gf[2] * dk[1] 
				+ gf[3] * (qkdi[1] - qidk[1]) + gf[4] * qir[1]
				 + gf[5] * qkr[1] + gf[6] * (qiqkr[1] + qkqir[
				1]);
			ftm2[2] = gf[0] * zr + gf[1] * di[2] + gf[2] * dk[2] 
				+ gf[3] * (qkdi[2] - qidk[2]) + gf[4] * qir[2]
				 + gf[5] * qkr[2] + gf[6] * (qiqkr[2] + qkqir[
				2]);

/*     get the permanent force without screening */

			ftm2r[0] = gfr[0] * xr + gfr[1] * di[0] + gfr[2] * dk[
				0] + gfr[3] * (qkdi[0] - qidk[0]) + gfr[4] * 
				qir[0] + gfr[5] * qkr[0] + gfr[6] * (qiqkr[0] 
				+ qkqir[0]);
			ftm2r[1] = gfr[0] * yr + gfr[1] * di[1] + gfr[2] * dk[
				1] + gfr[3] * (qkdi[1] - qidk[1]) + gfr[4] * 
				qir[1] + gfr[5] * qkr[1] + gfr[6] * (qiqkr[1] 
				+ qkqir[1]);
			ftm2r[2] = gfr[0] * zr + gfr[1] * di[2] + gfr[2] * dk[
				2] + gfr[3] * (qkdi[2] - qidk[2]) + gfr[4] * 
				qir[2] + gfr[5] * qkr[2] + gfr[6] * (qiqkr[2] 
				+ qkqir[2]);

/*     get the induced force with screening */

			ftm2i[0] = gfi[0] * xr + (gfi[1] * (uind_ref(1, i__) 
				+ uinp_ref(1, i__)) + bn[2] * (sci[3] * 
				uinp_ref(1, i__) + scip[3] * uind_ref(1, i__))
				 + gfi[2] * (uind_ref(1, k) + uinp_ref(1, k)) 
				+ bn[2] * (sci[2] * uinp_ref(1, k) + scip[2] *
				 uind_ref(1, k)) + (sci[3] + scip[3]) * bn[2] 
				* di[0] + (sci[2] + scip[2]) * bn[2] * dk[0] 
				+ gfi[3] * (qkui[0] + qkuip[0] - qiuk[0] - 
				qiukp[0])) * .5 + gfi[4] * qir[0] + gfi[5] * 
				qkr[0];
			ftm2i[1] = gfi[0] * yr + (gfi[1] * (uind_ref(2, i__) 
				+ uinp_ref(2, i__)) + bn[2] * (sci[3] * 
				uinp_ref(2, i__) + scip[3] * uind_ref(2, i__))
				 + gfi[2] * (uind_ref(2, k) + uinp_ref(2, k)) 
				+ bn[2] * (sci[2] * uinp_ref(2, k) + scip[2] *
				 uind_ref(2, k)) + (sci[3] + scip[3]) * bn[2] 
				* di[1] + (sci[2] + scip[2]) * bn[2] * dk[1] 
				+ gfi[3] * (qkui[1] + qkuip[1] - qiuk[1] - 
				qiukp[1])) * .5 + gfi[4] * qir[1] + gfi[5] * 
				qkr[1];
			ftm2i[2] = gfi[0] * zr + (gfi[1] * (uind_ref(3, i__) 
				+ uinp_ref(3, i__)) + bn[2] * (sci[3] * 
				uinp_ref(3, i__) + scip[3] * uind_ref(3, i__))
				 + gfi[2] * (uind_ref(3, k) + uinp_ref(3, k)) 
				+ bn[2] * (sci[2] * uinp_ref(3, k) + scip[2] *
				 uind_ref(3, k)) + (sci[3] + scip[3]) * bn[2] 
				* di[2] + (sci[2] + scip[2]) * bn[2] * dk[2] 
				+ gfi[3] * (qkui[2] + qkuip[2] - qiuk[2] - 
				qiukp[2])) * .5 + gfi[4] * qir[2] + gfi[5] * 
				qkr[2];

/*     get the induced force without screening */

			ftm2ri[0] = gfri[0] * xr + (-rr3 * ck * (uind_ref(1, 
				i__) * psc3 + uinp_ref(1, i__) * dsc3) + rr5 *
				 sc[3] * (uind_ref(1, i__) * psc5 + uinp_ref(
				1, i__) * dsc5) - rr7 * sc[5] * (uind_ref(1, 
				i__) * psc7 + uinp_ref(1, i__) * dsc7)) * .5 
				+ (rr3 * ci * (uind_ref(1, k) * psc3 + 
				uinp_ref(1, k) * dsc3) + rr5 * sc[2] * (
				uind_ref(1, k) * psc5 + uinp_ref(1, k) * dsc5)
				 + rr7 * sc[4] * (uind_ref(1, k) * psc7 + 
				uinp_ref(1, k) * dsc7)) * .5 + rr5 * usc5 * (
				sci[3] * uinp_ref(1, i__) + scip[3] * 
				uind_ref(1, i__) + sci[2] * uinp_ref(1, k) + 
				scip[2] * uind_ref(1, k)) * .5 + (sci[3] * 
				psc5 + scip[3] * dsc5) * .5 * rr5 * di[0] + (
				sci[2] * psc5 + scip[2] * dsc5) * .5 * rr5 * 
				dk[0] + gfri[3] * .5 * ((qkui[0] - qiuk[0]) * 
				psc5 + (qkuip[0] - qiukp[0]) * dsc5) + gfri[4]
				 * qir[0] + gfri[5] * qkr[0];
			ftm2ri[1] = gfri[0] * yr + (-rr3 * ck * (uind_ref(2, 
				i__) * psc3 + uinp_ref(2, i__) * dsc3) + rr5 *
				 sc[3] * (uind_ref(2, i__) * psc5 + uinp_ref(
				2, i__) * dsc5) - rr7 * sc[5] * (uind_ref(2, 
				i__) * psc7 + uinp_ref(2, i__) * dsc7)) * .5 
				+ (rr3 * ci * (uind_ref(2, k) * psc3 + 
				uinp_ref(2, k) * dsc3) + rr5 * sc[2] * (
				uind_ref(2, k) * psc5 + uinp_ref(2, k) * dsc5)
				 + rr7 * sc[4] * (uind_ref(2, k) * psc7 + 
				uinp_ref(2, k) * dsc7)) * .5 + rr5 * usc5 * (
				sci[3] * uinp_ref(2, i__) + scip[3] * 
				uind_ref(2, i__) + sci[2] * uinp_ref(2, k) + 
				scip[2] * uind_ref(2, k)) * .5 + (sci[3] * 
				psc5 + scip[3] * dsc5) * .5 * rr5 * di[1] + (
				sci[2] * psc5 + scip[2] * dsc5) * .5 * rr5 * 
				dk[1] + gfri[3] * .5 * ((qkui[1] - qiuk[1]) * 
				psc5 + (qkuip[1] - qiukp[1]) * dsc5) + gfri[4]
				 * qir[1] + gfri[5] * qkr[1];
			ftm2ri[2] = gfri[0] * zr + (-rr3 * ck * (uind_ref(3, 
				i__) * psc3 + uinp_ref(3, i__) * dsc3) + rr5 *
				 sc[3] * (uind_ref(3, i__) * psc5 + uinp_ref(
				3, i__) * dsc5) - rr7 * sc[5] * (uind_ref(3, 
				i__) * psc7 + uinp_ref(3, i__) * dsc7)) * .5 
				+ (rr3 * ci * (uind_ref(3, k) * psc3 + 
				uinp_ref(3, k) * dsc3) + rr5 * sc[2] * (
				uind_ref(3, k) * psc5 + uinp_ref(3, k) * dsc5)
				 + rr7 * sc[4] * (uind_ref(3, k) * psc7 + 
				uinp_ref(3, k) * dsc7)) * .5 + rr5 * usc5 * (
				sci[3] * uinp_ref(3, i__) + scip[3] * 
				uind_ref(3, i__) + sci[2] * uinp_ref(3, k) + 
				scip[2] * uind_ref(3, k)) * .5 + (sci[3] * 
				psc5 + scip[3] * dsc5) * .5 * rr5 * di[2] + (
				sci[2] * psc5 + scip[2] * dsc5) * .5 * rr5 * 
				dk[2] + gfri[3] * .5 * ((qkui[2] - qiuk[2]) * 
				psc5 + (qkuip[2] - qiukp[2]) * dsc5) + gfri[4]
				 * qir[2] + gfri[5] * qkr[2];

/*     account for partially excluded induced interactions */

			temp3 = rr3 * .5 * ((gli[0] + gli[5]) * pscale[kk - 1]
				 + (glip[0] + glip[5]) * dscale[kk - 1]);
			temp5 = rr5 * .5 * ((gli[1] + gli[6]) * pscale[kk - 1]
				 + (glip[1] + glip[6]) * dscale[kk - 1]);
			temp7 = rr7 * .5 * (gli[2] * pscale[kk - 1] + glip[2] 
				* dscale[kk - 1]);
			fridmp[0] = temp3 * ddsc3[0] + temp5 * ddsc5[0] + 
				temp7 * ddsc7[0];
			fridmp[1] = temp3 * ddsc3[1] + temp5 * ddsc5[1] + 
				temp7 * ddsc7[1];
			fridmp[2] = temp3 * ddsc3[2] + temp5 * ddsc5[2] + 
				temp7 * ddsc7[2];

/*     find some scaling terms for induced-induced force */

			temp3 = rr3 * .5 * uscale[kk - 1] * scip[1];
			temp5 = rr5 * -.5 * uscale[kk - 1] * (sci[2] * scip[3]
				 + scip[2] * sci[3]);
			findmp[0] = temp3 * ddsc3[0] + temp5 * ddsc5[0];
			findmp[1] = temp3 * ddsc3[1] + temp5 * ddsc5[1];
			findmp[2] = temp3 * ddsc3[2] + temp5 * ddsc5[2];

/*     modify the forces for partially excluded interactions */

			ftm2i[0] = ftm2i[0] - fridmp[0] - findmp[0];
			ftm2i[1] = ftm2i[1] - fridmp[1] - findmp[1];
			ftm2i[2] = ftm2i[2] - fridmp[2] - findmp[2];

/*     correction to convert mutual to direct polarization force */

			if (s_cmp(polpot_1.poltyp, "DIRECT", (ftnlen)6, (
				ftnlen)6) == 0) {
			    gfd = (bn[2] * scip[1] - bn[3] * (scip[2] * sci[3]
				     + sci[2] * scip[3])) * .5;
			    gfdr = (rr5 * scip[1] * usc3 - rr7 * (scip[2] * 
				    sci[3] + sci[2] * scip[3]) * usc5) * .5;
			    ftm2i[0] = ftm2i[0] - gfd * xr - bn[2] * .5 * (
				    sci[3] * uinp_ref(1, i__) + scip[3] * 
				    uind_ref(1, i__) + sci[2] * uinp_ref(1, k)
				     + scip[2] * uind_ref(1, k));
			    ftm2i[1] = ftm2i[1] - gfd * yr - bn[2] * .5 * (
				    sci[3] * uinp_ref(2, i__) + scip[3] * 
				    uind_ref(2, i__) + sci[2] * uinp_ref(2, k)
				     + scip[2] * uind_ref(2, k));
			    ftm2i[2] = ftm2i[2] - gfd * zr - bn[2] * .5 * (
				    sci[3] * uinp_ref(3, i__) + scip[3] * 
				    uind_ref(3, i__) + sci[2] * uinp_ref(3, k)
				     + scip[2] * uind_ref(3, k));
			    fdir[0] = gfdr * xr + usc5 * .5 * rr5 * (sci[3] * 
				    uinp_ref(1, i__) + scip[3] * uind_ref(1, 
				    i__) + sci[2] * uinp_ref(1, k) + scip[2] *
				     uind_ref(1, k));
			    fdir[1] = gfdr * yr + usc5 * .5 * rr5 * (sci[3] * 
				    uinp_ref(2, i__) + scip[3] * uind_ref(2, 
				    i__) + sci[2] * uinp_ref(2, k) + scip[2] *
				     uind_ref(2, k));
			    fdir[2] = gfdr * zr + usc5 * .5 * rr5 * (sci[3] * 
				    uinp_ref(3, i__) + scip[3] * uind_ref(3, 
				    i__) + sci[2] * uinp_ref(3, k) + scip[2] *
				     uind_ref(3, k));
			    ftm2i[0] = ftm2i[0] + fdir[0] + findmp[0];
			    ftm2i[1] = ftm2i[1] + fdir[1] + findmp[1];
			    ftm2i[2] = ftm2i[2] + fdir[2] + findmp[2];
			}

/*     intermediate variables for induced torque terms */

			gti[1] = bn[2] * .5 * (sci[3] + scip[3]);
			gti[2] = bn[2] * .5 * (sci[2] + scip[2]);
			gti[3] = gfi[3];
			gti[4] = gfi[4];
			gti[5] = gfi[5];
			gtri[1] = rr5 * .5 * (sci[3] * psc5 + scip[3] * dsc5);
			gtri[2] = rr5 * .5 * (sci[2] * psc5 + scip[2] * dsc5);
			gtri[3] = gfri[3];
			gtri[4] = gfri[4];
			gtri[5] = gfri[5];

/*     get the permanent torque with screening */

			ttm2[0] = -bn[1] * dixdk[0] + gf[1] * dixr[0] + gf[3] 
				* (dixqkr[0] + dkxqir[0] + rxqidk[0] - qixqk[
				0] * 2.) - gf[4] * rxqir[0] - gf[6] * (rxqikr[
				0] + qkrxqir[0]);
			ttm2[1] = -bn[1] * dixdk[1] + gf[1] * dixr[1] + gf[3] 
				* (dixqkr[1] + dkxqir[1] + rxqidk[1] - qixqk[
				1] * 2.) - gf[4] * rxqir[1] - gf[6] * (rxqikr[
				1] + qkrxqir[1]);
			ttm2[2] = -bn[1] * dixdk[2] + gf[1] * dixr[2] + gf[3] 
				* (dixqkr[2] + dkxqir[2] + rxqidk[2] - qixqk[
				2] * 2.) - gf[4] * rxqir[2] - gf[6] * (rxqikr[
				2] + qkrxqir[2]);
			ttm3[0] = bn[1] * dixdk[0] + gf[2] * dkxr[0] - gf[3] *
				 (dixqkr[0] + dkxqir[0] + rxqkdi[0] - qixqk[0]
				 * 2.) - gf[5] * rxqkr[0] - gf[6] * (rxqkir[0]
				 - qkrxqir[0]);
			ttm3[1] = bn[1] * dixdk[1] + gf[2] * dkxr[1] - gf[3] *
				 (dixqkr[1] + dkxqir[1] + rxqkdi[1] - qixqk[1]
				 * 2.) - gf[5] * rxqkr[1] - gf[6] * (rxqkir[1]
				 - qkrxqir[1]);
			ttm3[2] = bn[1] * dixdk[2] + gf[2] * dkxr[2] - gf[3] *
				 (dixqkr[2] + dkxqir[2] + rxqkdi[2] - qixqk[2]
				 * 2.) - gf[5] * rxqkr[2] - gf[6] * (rxqkir[2]
				 - qkrxqir[2]);

/*     get the permanent torque without screening */

			ttm2r[0] = -rr3 * dixdk[0] + gfr[1] * dixr[0] - gfr[4]
				 * rxqir[0] + gfr[3] * (dixqkr[0] + dkxqir[0] 
				+ rxqidk[0] - qixqk[0] * 2.) - gfr[6] * (
				rxqikr[0] + qkrxqir[0]);
			ttm2r[1] = -rr3 * dixdk[1] + gfr[1] * dixr[1] - gfr[4]
				 * rxqir[1] + gfr[3] * (dixqkr[1] + dkxqir[1] 
				+ rxqidk[1] - qixqk[1] * 2.) - gfr[6] * (
				rxqikr[1] + qkrxqir[1]);
			ttm2r[2] = -rr3 * dixdk[2] + gfr[1] * dixr[2] - gfr[4]
				 * rxqir[2] + gfr[3] * (dixqkr[2] + dkxqir[2] 
				+ rxqidk[2] - qixqk[2] * 2.) - gfr[6] * (
				rxqikr[2] + qkrxqir[2]);
			ttm3r[0] = rr3 * dixdk[0] + gfr[2] * dkxr[0] - gfr[5] 
				* rxqkr[0] - gfr[3] * (dixqkr[0] + dkxqir[0] 
				+ rxqkdi[0] - qixqk[0] * 2.) - gfr[6] * (
				rxqkir[0] - qkrxqir[0]);
			ttm3r[1] = rr3 * dixdk[1] + gfr[2] * dkxr[1] - gfr[5] 
				* rxqkr[1] - gfr[3] * (dixqkr[1] + dkxqir[1] 
				+ rxqkdi[1] - qixqk[1] * 2.) - gfr[6] * (
				rxqkir[1] - qkrxqir[1]);
			ttm3r[2] = rr3 * dixdk[2] + gfr[2] * dkxr[2] - gfr[5] 
				* rxqkr[2] - gfr[3] * (dixqkr[2] + dkxqir[2] 
				+ rxqkdi[2] - qixqk[2] * 2.) - gfr[6] * (
				rxqkir[2] - qkrxqir[2]);

/*     get the induced torque with screening */

			ttm2i[0] = -bn[1] * (dixuk[0] + dixukp[0]) * .5 + gti[
				1] * dixr[0] + gti[3] * (ukxqir[0] + rxqiuk[0]
				 + ukxqirp[0] + rxqiukp[0]) * .5 - gti[4] * 
				rxqir[0];
			ttm2i[1] = -bn[1] * (dixuk[1] + dixukp[1]) * .5 + gti[
				1] * dixr[1] + gti[3] * (ukxqir[1] + rxqiuk[1]
				 + ukxqirp[1] + rxqiukp[1]) * .5 - gti[4] * 
				rxqir[1];
			ttm2i[2] = -bn[1] * (dixuk[2] + dixukp[2]) * .5 + gti[
				1] * dixr[2] + gti[3] * (ukxqir[2] + rxqiuk[2]
				 + ukxqirp[2] + rxqiukp[2]) * .5 - gti[4] * 
				rxqir[2];
			ttm3i[0] = -bn[1] * (dkxui[0] + dkxuip[0]) * .5 + gti[
				2] * dkxr[0] - gti[3] * (uixqkr[0] + rxqkui[0]
				 + uixqkrp[0] + rxqkuip[0]) * .5 - gti[5] * 
				rxqkr[0];
			ttm3i[1] = -bn[1] * (dkxui[1] + dkxuip[1]) * .5 + gti[
				2] * dkxr[1] - gti[3] * (uixqkr[1] + rxqkui[1]
				 + uixqkrp[1] + rxqkuip[1]) * .5 - gti[5] * 
				rxqkr[1];
			ttm3i[2] = -bn[1] * (dkxui[2] + dkxuip[2]) * .5 + gti[
				2] * dkxr[2] - gti[3] * (uixqkr[2] + rxqkui[2]
				 + uixqkrp[2] + rxqkuip[2]) * .5 - gti[5] * 
				rxqkr[2];

/*     get the induced torque without screening */

			ttm2ri[0] = -rr3 * (dixuk[0] * psc3 + dixukp[0] * 
				dsc3) * .5 + gtri[1] * dixr[0] + gtri[3] * ((
				ukxqir[0] + rxqiuk[0]) * psc5 + (ukxqirp[0] + 
				rxqiukp[0]) * dsc5) * .5 - gtri[4] * rxqir[0];
			ttm2ri[1] = -rr3 * (dixuk[1] * psc3 + dixukp[1] * 
				dsc3) * .5 + gtri[1] * dixr[1] + gtri[3] * ((
				ukxqir[1] + rxqiuk[1]) * psc5 + (ukxqirp[1] + 
				rxqiukp[1]) * dsc5) * .5 - gtri[4] * rxqir[1];
			ttm2ri[2] = -rr3 * (dixuk[2] * psc3 + dixukp[2] * 
				dsc3) * .5 + gtri[1] * dixr[2] + gtri[3] * ((
				ukxqir[2] + rxqiuk[2]) * psc5 + (ukxqirp[2] + 
				rxqiukp[2]) * dsc5) * .5 - gtri[4] * rxqir[2];
			ttm3ri[0] = -rr3 * (dkxui[0] * psc3 + dkxuip[0] * 
				dsc3) * .5 + gtri[2] * dkxr[0] - gtri[3] * ((
				uixqkr[0] + rxqkui[0]) * psc5 + (uixqkrp[0] + 
				rxqkuip[0]) * dsc5) * .5 - gtri[5] * rxqkr[0];
			ttm3ri[1] = -rr3 * (dkxui[1] * psc3 + dkxuip[1] * 
				dsc3) * .5 + gtri[2] * dkxr[1] - gtri[3] * ((
				uixqkr[1] + rxqkui[1]) * psc5 + (uixqkrp[1] + 
				rxqkuip[1]) * dsc5) * .5 - gtri[5] * rxqkr[1];
			ttm3ri[2] = -rr3 * (dkxui[2] * psc3 + dkxuip[2] * 
				dsc3) * .5 + gtri[2] * dkxr[2] - gtri[3] * ((
				uixqkr[2] + rxqkui[2]) * psc5 + (uixqkrp[2] + 
				rxqkuip[2]) * dsc5) * .5 - gtri[5] * rxqkr[2];

/*     handle the case where scaling is used */

			if (bound_1.use_polymer__ && r2 <= bound_1.polycut2) {
			    for (j = 1; j <= 3; ++j) {
				ftm2[j - 1] = f * (ftm2[j - 1] - (1. - mscale[
					kk - 1]) * ftm2r[j - 1]);
				ftm2i[j - 1] = f * (ftm2i[j - 1] - ftm2ri[j - 
					1]);
				ttm2[j - 1] = f * (ttm2[j - 1] - (1. - mscale[
					kk - 1]) * ttm2r[j - 1]);
				ttm2i[j - 1] = f * (ttm2i[j - 1] - ttm2ri[j - 
					1]);
				ttm3[j - 1] = f * (ttm3[j - 1] - (1. - mscale[
					kk - 1]) * ttm3r[j - 1]);
				ttm3i[j - 1] = f * (ttm3i[j - 1] - ttm3ri[j - 
					1]);
			    }
			} else {
			    for (j = 1; j <= 3; ++j) {
				ftm2[j - 1] = f * ftm2[j - 1];
				ftm2i[j - 1] = f * (ftm2i[j - 1] - ftm2ri[j - 
					1]);
				ttm2[j - 1] = f * ttm2[j - 1];
				ttm2i[j - 1] = f * (ttm2i[j - 1] - ttm2ri[j - 
					1]);
				ttm3[j - 1] = f * ttm3[j - 1];
				ttm3i[j - 1] = f * (ttm3i[j - 1] - ttm3ri[j - 
					1]);
			    }
			}
			if (ii == kk) {
			    for (j = 1; j <= 3; ++j) {
				ftm2[j - 1] *= .5;
				ftm2i[j - 1] *= .5;
				ttm2[j - 1] *= .5;
				ttm2i[j - 1] *= .5;
				ttm3[j - 1] *= .5;
				ttm3i[j - 1] *= .5;
			    }
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

/*     increment the internal virial tensor components */

			iaz = mpole_1.zaxis[i__ - 1];
			iax = mpole_1.xaxis[i__ - 1];
			iay = mpole_1.yaxis[i__ - 1];
			kaz = mpole_1.zaxis[k - 1];
			kax = mpole_1.xaxis[k - 1];
			kay = mpole_1.yaxis[k - 1];
			if (iaz == 0) {
			    iaz = ii;
			}
			if (iax == 0) {
			    iax = ii;
			}
			if (iay == 0) {
			    iay = ii;
			}
			if (kaz == 0) {
			    kaz = kk;
			}
			if (kax == 0) {
			    kax = kk;
			}
			if (kay == 0) {
			    kay = kk;
			}
			xiz = atoms_1.x[iaz - 1] - atoms_1.x[ii - 1];
			yiz = atoms_1.y[iaz - 1] - atoms_1.y[ii - 1];
			ziz = atoms_1.z__[iaz - 1] - atoms_1.z__[ii - 1];
			xix = atoms_1.x[iax - 1] - atoms_1.x[ii - 1];
			yix = atoms_1.y[iax - 1] - atoms_1.y[ii - 1];
			zix = atoms_1.z__[iax - 1] - atoms_1.z__[ii - 1];
			xiy = atoms_1.x[iay - 1] - atoms_1.x[ii - 1];
			yiy = atoms_1.y[iay - 1] - atoms_1.y[ii - 1];
			ziy = atoms_1.z__[iay - 1] - atoms_1.z__[ii - 1];
			xkz = atoms_1.x[kaz - 1] - atoms_1.x[kk - 1];
			ykz = atoms_1.y[kaz - 1] - atoms_1.y[kk - 1];
			zkz = atoms_1.z__[kaz - 1] - atoms_1.z__[kk - 1];
			xkx = atoms_1.x[kax - 1] - atoms_1.x[kk - 1];
			ykx = atoms_1.y[kax - 1] - atoms_1.y[kk - 1];
			zkx = atoms_1.z__[kax - 1] - atoms_1.z__[kk - 1];
			xky = atoms_1.x[kay - 1] - atoms_1.x[kk - 1];
			yky = atoms_1.y[kay - 1] - atoms_1.y[kk - 1];
			zky = atoms_1.z__[kay - 1] - atoms_1.z__[kk - 1];
			vxx = -xr * (ftm2[0] + ftm2i[0]) + xix * frcxi[0] + 
				xiy * frcyi[0] + xiz * frczi[0] + xkx * frcxk[
				0] + xky * frcyk[0] + xkz * frczk[0];
			vyx = -yr * (ftm2[0] + ftm2i[0]) + yix * frcxi[0] + 
				yiy * frcyi[0] + yiz * frczi[0] + ykx * frcxk[
				0] + yky * frcyk[0] + ykz * frczk[0];
			vzx = -zr * (ftm2[0] + ftm2i[0]) + zix * frcxi[0] + 
				ziy * frcyi[0] + ziz * frczi[0] + zkx * frcxk[
				0] + zky * frcyk[0] + zkz * frczk[0];
			vyy = -yr * (ftm2[1] + ftm2i[1]) + yix * frcxi[1] + 
				yiy * frcyi[1] + yiz * frczi[1] + ykx * frcxk[
				1] + yky * frcyk[1] + ykz * frczk[1];
			vzy = -zr * (ftm2[1] + ftm2i[1]) + zix * frcxi[1] + 
				ziy * frcyi[1] + ziz * frczi[1] + zkx * frcxk[
				1] + zky * frcyk[1] + zkz * frczk[1];
			vzz = -zr * (ftm2[2] + ftm2i[2]) + zix * frcxi[2] + 
				ziy * frcyi[2] + ziz * frczi[2] + zkx * frcxk[
				2] + zky * frcyk[2] + zkz * frczk[2];
			vir_ref(1, 1) = vir_ref(1, 1) + vxx;
			vir_ref(2, 1) = vir_ref(2, 1) + vyx;
			vir_ref(3, 1) = vir_ref(3, 1) + vzx;
			vir_ref(1, 2) = vir_ref(1, 2) + vyx;
			vir_ref(2, 2) = vir_ref(2, 2) + vyy;
			vir_ref(3, 2) = vir_ref(3, 2) + vzy;
			vir_ref(1, 3) = vir_ref(1, 3) + vzx;
			vir_ref(2, 3) = vir_ref(2, 3) + vzy;
			vir_ref(3, 3) = vir_ref(3, 3) + vzz;
		    }
		}
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
    }
    return 0;
} /* ereal1c_ */

#undef rpole_ref
#undef uinp_ref
#undef uind_ref
#undef vir_ref
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




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine empole1d  --  Ewald multipole derivs via list  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "empole1d" calculates the multipole and dipole polarization */
/*     energy and derivatives with respect to Cartesian coordinates */
/*     using particle mesh Ewald summation and a neighbor list */


/* Subroutine */ int empole1d_(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j;
    static doublereal ci, ei;
    static integer ii;
    static doublereal xd, yd, zd, xq, yq, zq, xu, yu, zu, xv, yv, zv, cii, 
	    dii, qii, dix, diy, uii, diz, uix, uiy, uiz, trq[3], xup, yup, 
	    zup, frcx[3], frcy[3], frcz[3], term, trqi[3], qixx, qixy, qixz, 
	    qiyy, qiyz, qizz, fterm, vterm;
    extern /* Subroutine */ int induce_(void);
    static doublereal eintra;
    extern /* Subroutine */ int ereal1d_(doublereal *), torque_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal xdfield, ydfield, zdfield;
    extern /* Subroutine */ int chkpole_(void);
    static doublereal xufield, yufield, zufield;
    extern /* Subroutine */ int rotpole_(void), emrecip1_(void);


#define dem_ref(a_1,a_2) deriv_1.dem[(a_2)*3 + a_1 - 4]
#define dep_ref(a_1,a_2) deriv_1.dep[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  ewald.i  --  parameters and options for Ewald summation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     aewald     Ewald convergence coefficient value (Ang-1) */
/*     boundary   Ewald boundary condition; none, tinfoil or vacuum */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  inter.i  --  sum of intermolecular energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     einter   total intermolecular potential energy */




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

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  virial.i  --  components of internal virial tensor  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     vir    total internal virial Cartesian tensor components */




/*     zero out multipole and polarization energy and derivatives */

    energi_1.em = 0.;
    energi_1.ep = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    dem_ref(j, i__) = 0.;
	    dep_ref(j, i__) = 0.;
	}
    }

/*     set the energy unit conversion factor */

    f = chgpot_1.electric / chgpot_1.dielec;

/*     check the sign of multipole components at chiral sites */

    chkpole_();

/*     rotate the multipole components into the global frame */

    rotpole_();

/*     compute the induced dipole moment at each atom */

    induce_();

/*     compute the reciprocal space part of the Ewald summation */

    emrecip1_();

/*     compute the real space part of the Ewald summation */

    ereal1d_(&eintra);

/*     compute the Ewald self-energy term over all the atoms */

    term = ewald_1.aewald * 2. * ewald_1.aewald;
    fterm = -f * ewald_1.aewald / 1.772453850905516027;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ci = rpole_ref(1, i__);
	dix = rpole_ref(2, i__);
	diy = rpole_ref(3, i__);
	diz = rpole_ref(4, i__);
	qixx = rpole_ref(5, i__);
	qixy = rpole_ref(6, i__);
	qixz = rpole_ref(7, i__);
	qiyy = rpole_ref(9, i__);
	qiyz = rpole_ref(10, i__);
	qizz = rpole_ref(13, i__);
	uix = uind_ref(1, i__);
	uiy = uind_ref(2, i__);
	uiz = uind_ref(3, i__);
	cii = ci * ci;
	dii = dix * dix + diy * diy + diz * diz;
	qii = qixx * qixx + qiyy * qiyy + qizz * qizz + (qixy * qixy + qixz * 
		qixz + qiyz * qiyz) * 2.;
	uii = dix * uix + diy * uiy + diz * uiz;
	e = fterm * (cii + term * (dii / 3. + term * 2. * qii / 5.));
	ei = fterm * term * uii / 3.;
	energi_1.em += e;
	energi_1.ep += ei;
    }

/*     compute the self-energy torque term due to induced dipole */

    trq[0] = 0.;
    trq[1] = 0.;
    trq[2] = 0.;
/* Computing 3rd power */
    d__1 = ewald_1.aewald;
    term = f * 1.3333333333333333 * (d__1 * (d__1 * d__1)) / 
	    1.772453850905516027;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dix = rpole_ref(2, i__);
	diy = rpole_ref(3, i__);
	diz = rpole_ref(4, i__);
	uix = (uind_ref(1, i__) + uinp_ref(1, i__)) * .5;
	uiy = (uind_ref(2, i__) + uinp_ref(2, i__)) * .5;
	uiz = (uind_ref(3, i__) + uinp_ref(3, i__)) * .5;
	trqi[0] = term * (diy * uiz - diz * uiy);
	trqi[1] = term * (diz * uix - dix * uiz);
	trqi[2] = term * (dix * uiy - diy * uix);
	torque_(&i__, trq, trqi, frcx, frcy, frcz);
    }

/*     compute the cell dipole boundary correction term */

    if (s_cmp(ewald_1.boundary, "VACUUM", (ftnlen)7, (ftnlen)6) == 0) {
	xd = 0.;
	yd = 0.;
	zd = 0.;
	xu = 0.;
	yu = 0.;
	zu = 0.;
	xup = 0.;
	yup = 0.;
	zup = 0.;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    xd = xd + rpole_ref(2, i__) + rpole_ref(1, i__) * atoms_1.x[ii - 
		    1];
	    yd = yd + rpole_ref(3, i__) + rpole_ref(1, i__) * atoms_1.y[ii - 
		    1];
	    zd = zd + rpole_ref(4, i__) + rpole_ref(1, i__) * atoms_1.z__[ii 
		    - 1];
	    xu += uind_ref(1, i__);
	    yu += uind_ref(2, i__);
	    zu += uind_ref(3, i__);
	    xup += uinp_ref(1, i__);
	    yup += uinp_ref(2, i__);
	    zup += uinp_ref(3, i__);
	}
	term = f * .66666666666666663 * (3.141592653589793238 / 
		boxes_1.volbox);
	energi_1.em += term * (xd * xd + yd * yd + zd * zd);
	energi_1.ep += term * (xd * xu + yd * yu + zd * zu);
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    dem_ref(1, ii) = dem_ref(1, ii) + term * 2. * rpole_ref(1, i__) * 
		    xd;
	    dem_ref(2, ii) = dem_ref(2, ii) + term * 2. * rpole_ref(1, i__) * 
		    yd;
	    dem_ref(3, ii) = dem_ref(3, ii) + term * 2. * rpole_ref(1, i__) * 
		    zd;
	    dep_ref(1, ii) = dep_ref(1, ii) + term * rpole_ref(1, i__) * (xu 
		    + xup);
	    dep_ref(2, ii) = dep_ref(2, ii) + term * rpole_ref(1, i__) * (yu 
		    + yup);
	    dep_ref(3, ii) = dep_ref(3, ii) + term * rpole_ref(1, i__) * (zu 
		    + zup);
	}
	xdfield = term * -2. * xd;
	ydfield = term * -2. * yd;
	zdfield = term * -2. * zd;
	xufield = -term * (xu + xup);
	yufield = -term * (yu + yup);
	zufield = -term * (zu + zup);
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    trq[0] = rpole_ref(3, i__) * zdfield - rpole_ref(4, i__) * 
		    ydfield;
	    trq[1] = rpole_ref(4, i__) * xdfield - rpole_ref(2, i__) * 
		    zdfield;
	    trq[2] = rpole_ref(2, i__) * ydfield - rpole_ref(3, i__) * 
		    xdfield;
	    trqi[0] = rpole_ref(3, i__) * zufield - rpole_ref(4, i__) * 
		    yufield;
	    trqi[1] = rpole_ref(4, i__) * xufield - rpole_ref(2, i__) * 
		    zufield;
	    trqi[2] = rpole_ref(2, i__) * yufield - rpole_ref(3, i__) * 
		    xufield;
	    torque_(&i__, trq, trqi, frcx, frcy, frcz);
	}

/*     boundary correction to virial due to overall cell dipole */

	xd = 0.;
	yd = 0.;
	zd = 0.;
	xq = 0.;
	yq = 0.;
	zq = 0.;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    xd += rpole_ref(2, i__);
	    yd += rpole_ref(3, i__);
	    zd += rpole_ref(4, i__);
	    xq += rpole_ref(1, i__) * atoms_1.x[ii - 1];
	    yq += rpole_ref(1, i__) * atoms_1.y[ii - 1];
	    zq += rpole_ref(1, i__) * atoms_1.z__[ii - 1];
	}
	xv = xq * (xd + (xu + xup) * .5);
	yv = yq * (yd + (yu + yup) * .5);
	zv = zq * (zd + (zu + zup) * .5);
	vterm = term * (xq * xq + yq * yq + zq * zq + (xv + yv + zv) * 2. + 
		xu * xup + yu * yup + zu * zup + xd * (xd + xu + xup) + yd * (
		yd + yu + yup) + zd * (zd + zu + zup));
	vir_ref(1, 1) = vir_ref(1, 1) + term * 2. * (xq * xq + xv) + vterm;
	vir_ref(2, 1) = vir_ref(2, 1) + term * 2. * (xq * yq + xv);
	vir_ref(3, 1) = vir_ref(3, 1) + term * 2. * (xq * zq + xv);
	vir_ref(1, 2) = vir_ref(1, 2) + term * 2. * (yq * xq + yv);
	vir_ref(2, 2) = vir_ref(2, 2) + term * 2. * (yq * yq + yv) + vterm;
	vir_ref(3, 2) = vir_ref(3, 2) + term * 2. * (yq * zq + yv);
	vir_ref(1, 3) = vir_ref(1, 3) + term * 2. * (zq * xq + zv);
	vir_ref(2, 3) = vir_ref(2, 3) + term * 2. * (zq * yq + zv);
	vir_ref(3, 3) = vir_ref(3, 3) + term * 2. * (zq * zq + zv) + vterm;
	if (s_cmp(polpot_1.poltyp, "DIRECT", (ftnlen)6, (ftnlen)6) == 0) {
	    vterm = term * (xu * xup + yu * yup + zu * zup);
	    vir_ref(1, 1) = vir_ref(1, 1) + vterm;
	    vir_ref(2, 2) = vir_ref(2, 2) + vterm;
	    vir_ref(3, 3) = vir_ref(3, 3) + vterm;
	}
    }

/*     intermolecular energy is total minus intramolecular part */

    inter_1.einter = inter_1.einter + energi_1.em + energi_1.ep - eintra;
    return 0;
} /* empole1d_ */

#undef rpole_ref
#undef uinp_ref
#undef uind_ref
#undef vir_ref
#undef dep_ref
#undef dem_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine ereal1d  --  ewald real space derivs via list  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "ereal1d" evaluates the real space portion of the regular Ewald */
/*     summation energy and gradient due to atomic multipole interactions */
/*     and dipole polarizability */


/* Subroutine */ int ereal1d_(doublereal *eintra)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__, r2;
    static integer ii, kk;
    static doublereal ei, ci, di[3], qi[9], ck, dk[3], qk[9], bn[6], sc[10], 
	    gl[9], xr, yr, zr, gf[7], rr1, rr3, rr5, rr7, rr9, gfd, gfi[6];
    static integer kkk, iax, iay, iaz, kax, kay, kaz;
    static doublereal pdi, pti, xix, yix, zix, xiy, yiy, ziy, xiz, yiz, ziz, 
	    xkx, ykx, zkx, xky, yky, zky, xkz, ykz, zkz, rr11, erl, dsc3, 
	    dsc5, vxx, dsc7, vyy, vzz, vyx, vzx, vzy, qir[3], qkr[3], psc3, 
	    ftm2[3], psc5, sci[8], psc7, usc3, usc5, gli[7], gfr[7], gti[6], 
	    ttm2[3], ttm3[3], bfac;
    extern doublereal erfc_(doublereal *);
    static doublereal damp, gfdr, fdir[3], qidk[3], qkdi[3], erli, glip[7], 
	    scip[8], dixr[3], dkxr[3], qiuk[3], qkui[3], gfri[6], gtri[6];
    static logical dorl;
    static doublereal ddsc3[3], ddsc5[3], ddsc7[3], exp2a, ftm2i[3], alsq2, 
	    temp3, temp5, ftm2r[3], temp7, ttm2i[3], ttm3i[3], ttm2r[3], 
	    ttm3r[3];
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal dixdk[3], frcxi[3], frcxk[3], frcyi[3], frcyk[3], frczi[
	    3], frczk[3], dkxui[3], dixuk[3], qiukp[3], qkuip[3], qiqkr[3], 
	    qkqir[3], scale3, qixqk[3], scale5, rxqir[3], scale7, rxqkr[3];
    static logical dorli;
    static doublereal alsq2n, ftm2ri[3], ttm2ri[3], ttm3ri[3], dscale[25000], 
	    pgamma, mscale[25000], ralpha, pscale[25000], uscale[25000], 
	    findmp[3], fridmp[3], dixukp[3], dkxuip[3], uixqkr[3], ukxqir[3], 
	    rxqiuk[3], rxqkui[3], dixqkr[3], dkxqir[3], rxqikr[3], rxqkir[3], 
	    rxqidk[3], rxqkdi[3];
    extern /* Subroutine */ int switch_(char *, ftnlen), torque_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal expdamp, qkrxqir[3], uixqkrp[3], ukxqirp[3], rxqiukp[3],
	     rxqkuip[3];


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
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define elst_ref(a_1,a_2) neigh_1.elst[(a_2)*1200 + a_1 - 1201]
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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cell.i  --  periodic boundaries using replicated cells  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xcell    length of the a-axis of the complete replicated cell */
/*     ycell    length of the b-axis of the complete replicated cell */
/*     zcell    length of the c-axis of the complete replicated cell */
/*     xcell2   half the length of the a-axis of the replicated cell */
/*     ycell2   half the length of the b-axis of the replicated cell */
/*     zcell2   half the length of the c-axis of the replicated cell */
/*     ncell    total number of cell replicates for periodic boundaries */
/*     icell    offset along axes for each replicate periodic cell */




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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  ewald.i  --  parameters and options for Ewald summation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     aewald     Ewald convergence coefficient value (Ang-1) */
/*     boundary   Ewald boundary condition; none, tinfoil or vacuum */




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




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  neigh.i  --  pairwise neighbor list indices and storage  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     lbuffer     width of the neighbor list buffer region */
/*     lbuf2       motion squared needed to trigger list rebuild */
/*     vbuf2       square of vdw cutoff plus neighbor list buffer */
/*     cbuf2       square of charge cutoff plus neighbor list buffer */
/*     mbuf2       square of multipole cutoff plus neighbor list buffer */
/*     nvlst       number of sites in list for each vdw site */
/*     vlst        site numbers in neighbor list of each vdw site */
/*     nelst       number of sites in list for each electrostatic site */
/*     elst        site numbers in list of each electrostatic site */
/*     dovlst      logical flag to rebuild vdw neighbor list */
/*     doclst      logical flag to rebuild charge neighbor list */
/*     domlst      logical flag to rebuild multipole neighbor list */




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

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  virial.i  --  components of internal virial tensor  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     vir    total internal virial Cartesian tensor components */




/*     zero out the intramolecular portion of the Ewald energy */

    *eintra = 0.;
    if (mpole_1.npole == 0) {
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
    switch_("EWALD", (ftnlen)5);

/*     set the permanent multipole and induced dipole values */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
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
	i__2 = neigh_1.nelst[i__ - 1];
	for (kkk = 1; kkk <= i__2; ++kkk) {
	    k = elst_ref(kkk, i__);
	    kk = mpole_1.ipole[k - 1];
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

/*     calculate the real space error function terms */

		ralpha = ewald_1.aewald * r__;
		bn[0] = erfc_(&ralpha) / r__;
/* Computing 2nd power */
		d__1 = ewald_1.aewald;
		alsq2 = d__1 * d__1 * 2.;
		alsq2n = 0.;
		if (ewald_1.aewald > 0.) {
		    alsq2n = 1. / (ewald_1.aewald * 1.772453850905516027);
		}
/* Computing 2nd power */
		d__1 = ralpha;
		exp2a = exp(-(d__1 * d__1));
		for (j = 1; j <= 5; ++j) {
		    bfac = (doublereal) ((j << 1) - 1);
		    alsq2n = alsq2 * alsq2n;
		    bn[j] = (bfac * bn[j - 1] + alsq2n * exp2a) / r2;
		}

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
		dsc3 = 1. - scale3 * dscale[kk - 1];
		dsc5 = 1. - scale5 * dscale[kk - 1];
		dsc7 = 1. - scale7 * dscale[kk - 1];
		psc3 = 1. - scale3 * pscale[kk - 1];
		psc5 = 1. - scale5 * pscale[kk - 1];
		psc7 = 1. - scale7 * pscale[kk - 1];
		usc3 = 1. - scale3 * uscale[kk - 1];
		usc5 = 1. - scale5 * uscale[kk - 1];

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

/*     calculate the scalar products for permanent components */

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

/*     calculate the scalar products for induced components */

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

/*     compute the energy contributions for this interaction */

		e = bn[0] * gl[0] + bn[1] * (gl[1] + gl[6]) + bn[2] * (gl[2] 
			+ gl[7] + gl[8]) + bn[3] * (gl[3] + gl[5]) + bn[4] * 
			gl[4];
		ei = (bn[1] * (gli[0] + gli[5]) + bn[2] * (gli[1] + gli[6]) + 
			bn[3] * gli[2]) * .5;

/*     get the real energy without any screening function */

		erl = rr1 * gl[0] + rr3 * (gl[1] + gl[6]) + rr5 * (gl[2] + gl[
			7] + gl[8]) + rr7 * (gl[3] + gl[5]) + rr9 * gl[4];
		erl *= 1. - mscale[kk - 1];
		erli = (rr3 * (gli[0] + gli[5]) * psc3 + rr5 * (gli[1] + gli[
			6]) * psc5 + rr7 * gli[2] * psc7) * .5;
		e -= erl;
		ei -= erli;
		e = f * e;
		ei = f * ei;
		energi_1.em += e;
		energi_1.ep += ei;

/*     increment the total intramolecular energy; assumes */
/*     intramolecular distances are less than half of cell */
/*     length and less than the ewald cutoff */

		if (molcul_1.molcule[ii - 1] == molcul_1.molcule[kk - 1]) {
		    *eintra += mscale[kk - 1] * erl * f;
		    *eintra += pscale[kk - 1] * .5 * (rr3 * (gli[0] + gli[5]) 
			    * scale3 + rr5 * (gli[1] + gli[6]) * scale5 + rr7 
			    * gli[2] * scale7);
		}

/*     set flags to compute components without screening */

		dorl = FALSE_;
		dorli = FALSE_;
		if (mscale[kk - 1] != 1.) {
		    dorl = TRUE_;
		}
		if (psc3 != 0.) {
		    dorli = TRUE_;
		}
		if (dsc3 != 0.) {
		    dorli = TRUE_;
		}
		if (usc3 != 0.) {
		    dorli = TRUE_;
		}

/*     zero out force and torque components without screening */

		for (j = 1; j <= 3; ++j) {
		    ftm2r[j - 1] = 0.;
		    ftm2ri[j - 1] = 0.;
		    ttm2r[j - 1] = 0.;
		    ttm2ri[j - 1] = 0.;
		    ttm3r[j - 1] = 0.;
		    ttm3ri[j - 1] = 0.;
		}

/*     get the permanent force with screening */

		gf[0] = bn[1] * gl[0] + bn[2] * (gl[1] + gl[6]) + bn[3] * (gl[
			2] + gl[7] + gl[8]) + bn[4] * (gl[3] + gl[5]) + bn[5] 
			* gl[4];
		gf[1] = -ck * bn[1] + sc[3] * bn[2] - sc[5] * bn[3];
		gf[2] = ci * bn[1] + sc[2] * bn[2] + sc[4] * bn[3];
		gf[3] = bn[2] * 2.;
		gf[4] = (-ck * bn[2] + sc[3] * bn[3] - sc[5] * bn[4]) * 2.;
		gf[5] = (-ci * bn[2] - sc[2] * bn[3] - sc[4] * bn[4]) * 2.;
		gf[6] = bn[3] * 4.;
		ftm2[0] = gf[0] * xr + gf[1] * di[0] + gf[2] * dk[0] + gf[3] *
			 (qkdi[0] - qidk[0]) + gf[4] * qir[0] + gf[5] * qkr[0]
			 + gf[6] * (qiqkr[0] + qkqir[0]);
		ftm2[1] = gf[0] * yr + gf[1] * di[1] + gf[2] * dk[1] + gf[3] *
			 (qkdi[1] - qidk[1]) + gf[4] * qir[1] + gf[5] * qkr[1]
			 + gf[6] * (qiqkr[1] + qkqir[1]);
		ftm2[2] = gf[0] * zr + gf[1] * di[2] + gf[2] * dk[2] + gf[3] *
			 (qkdi[2] - qidk[2]) + gf[4] * qir[2] + gf[5] * qkr[2]
			 + gf[6] * (qiqkr[2] + qkqir[2]);

/*     get the permanent force without screening */

		if (dorl) {
		    gfr[0] = rr3 * gl[0] + rr5 * (gl[1] + gl[6]) + rr7 * (gl[
			    2] + gl[7] + gl[8]) + rr9 * (gl[3] + gl[5]) + 
			    rr11 * gl[4];
		    gfr[1] = -ck * rr3 + sc[3] * rr5 - sc[5] * rr7;
		    gfr[2] = ci * rr3 + sc[2] * rr5 + sc[4] * rr7;
		    gfr[3] = rr5 * 2.;
		    gfr[4] = (-ck * rr5 + sc[3] * rr7 - sc[5] * rr9) * 2.;
		    gfr[5] = (-ci * rr5 - sc[2] * rr7 - sc[4] * rr9) * 2.;
		    gfr[6] = rr7 * 4.;
		    ftm2r[0] = gfr[0] * xr + gfr[1] * di[0] + gfr[2] * dk[0] 
			    + gfr[3] * (qkdi[0] - qidk[0]) + gfr[4] * qir[0] 
			    + gfr[5] * qkr[0] + gfr[6] * (qiqkr[0] + qkqir[0])
			    ;
		    ftm2r[1] = gfr[0] * yr + gfr[1] * di[1] + gfr[2] * dk[1] 
			    + gfr[3] * (qkdi[1] - qidk[1]) + gfr[4] * qir[1] 
			    + gfr[5] * qkr[1] + gfr[6] * (qiqkr[1] + qkqir[1])
			    ;
		    ftm2r[2] = gfr[0] * zr + gfr[1] * di[2] + gfr[2] * dk[2] 
			    + gfr[3] * (qkdi[2] - qidk[2]) + gfr[4] * qir[2] 
			    + gfr[5] * qkr[2] + gfr[6] * (qiqkr[2] + qkqir[2])
			    ;
		}

/*     get the induced force with screening */

		gfi[0] = bn[2] * .5 * (gli[0] + glip[0] + gli[5] + glip[5]) + 
			bn[2] * .5 * scip[1] + bn[3] * .5 * (gli[1] + glip[1] 
			+ gli[6] + glip[6]) - bn[3] * .5 * (sci[2] * scip[3] 
			+ scip[2] * sci[3]) + bn[4] * .5 * (gli[2] + glip[2]);
		gfi[1] = -ck * bn[1] + sc[3] * bn[2] - sc[5] * bn[3];
		gfi[2] = ci * bn[1] + sc[2] * bn[2] + sc[4] * bn[3];
		gfi[3] = bn[2] * 2.;
		gfi[4] = bn[3] * (sci[3] + scip[3]);
		gfi[5] = -bn[3] * (sci[2] + scip[2]);
		ftm2i[0] = gfi[0] * xr + (gfi[1] * (uind_ref(1, i__) + 
			uinp_ref(1, i__)) + bn[2] * (sci[3] * uinp_ref(1, i__)
			 + scip[3] * uind_ref(1, i__)) + gfi[2] * (uind_ref(1,
			 k) + uinp_ref(1, k)) + bn[2] * (sci[2] * uinp_ref(1, 
			k) + scip[2] * uind_ref(1, k)) + (sci[3] + scip[3]) * 
			bn[2] * di[0] + (sci[2] + scip[2]) * bn[2] * dk[0] + 
			gfi[3] * (qkui[0] + qkuip[0] - qiuk[0] - qiukp[0])) * 
			.5 + gfi[4] * qir[0] + gfi[5] * qkr[0];
		ftm2i[1] = gfi[0] * yr + (gfi[1] * (uind_ref(2, i__) + 
			uinp_ref(2, i__)) + bn[2] * (sci[3] * uinp_ref(2, i__)
			 + scip[3] * uind_ref(2, i__)) + gfi[2] * (uind_ref(2,
			 k) + uinp_ref(2, k)) + bn[2] * (sci[2] * uinp_ref(2, 
			k) + scip[2] * uind_ref(2, k)) + (sci[3] + scip[3]) * 
			bn[2] * di[1] + (sci[2] + scip[2]) * bn[2] * dk[1] + 
			gfi[3] * (qkui[1] + qkuip[1] - qiuk[1] - qiukp[1])) * 
			.5 + gfi[4] * qir[1] + gfi[5] * qkr[1];
		ftm2i[2] = gfi[0] * zr + (gfi[1] * (uind_ref(3, i__) + 
			uinp_ref(3, i__)) + bn[2] * (sci[3] * uinp_ref(3, i__)
			 + scip[3] * uind_ref(3, i__)) + gfi[2] * (uind_ref(3,
			 k) + uinp_ref(3, k)) + bn[2] * (sci[2] * uinp_ref(3, 
			k) + scip[2] * uind_ref(3, k)) + (sci[3] + scip[3]) * 
			bn[2] * di[2] + (sci[2] + scip[2]) * bn[2] * dk[2] + 
			gfi[3] * (qkui[2] + qkuip[2] - qiuk[2] - qiukp[2])) * 
			.5 + gfi[4] * qir[2] + gfi[5] * qkr[2];

/*     get the induced force without screening */

		if (dorli) {
		    gfri[0] = rr5 * .5 * ((gli[0] + gli[5]) * psc3 + (glip[0] 
			    + glip[5]) * dsc3 + scip[1] * usc3) + rr7 * .5 * (
			    (gli[6] + gli[1]) * psc5 + (glip[6] + glip[1]) * 
			    dsc5 - (sci[2] * scip[3] + scip[2] * sci[3]) * 
			    usc5) + rr9 * .5 * (gli[2] * psc7 + glip[2] * 
			    dsc7);
		    gfri[1] = -rr3 * ck + rr5 * sc[3] - rr7 * sc[5];
		    gfri[2] = rr3 * ci + rr5 * sc[2] + rr7 * sc[4];
		    gfri[3] = rr5 * 2.;
		    gfri[4] = rr7 * (sci[3] * psc7 + scip[3] * dsc7);
		    gfri[5] = -rr7 * (sci[2] * psc7 + scip[2] * dsc7);
		    ftm2ri[0] = gfri[0] * xr + (-rr3 * ck * (uind_ref(1, i__) 
			    * psc3 + uinp_ref(1, i__) * dsc3) + rr5 * sc[3] * 
			    (uind_ref(1, i__) * psc5 + uinp_ref(1, i__) * 
			    dsc5) - rr7 * sc[5] * (uind_ref(1, i__) * psc7 + 
			    uinp_ref(1, i__) * dsc7)) * .5 + (rr3 * ci * (
			    uind_ref(1, k) * psc3 + uinp_ref(1, k) * dsc3) + 
			    rr5 * sc[2] * (uind_ref(1, k) * psc5 + uinp_ref(1,
			     k) * dsc5) + rr7 * sc[4] * (uind_ref(1, k) * 
			    psc7 + uinp_ref(1, k) * dsc7)) * .5 + rr5 * usc5 *
			     (sci[3] * uinp_ref(1, i__) + scip[3] * uind_ref(
			    1, i__) + sci[2] * uinp_ref(1, k) + scip[2] * 
			    uind_ref(1, k)) * .5 + (sci[3] * psc5 + scip[3] * 
			    dsc5) * .5 * rr5 * di[0] + (sci[2] * psc5 + scip[
			    2] * dsc5) * .5 * rr5 * dk[0] + gfri[3] * .5 * ((
			    qkui[0] - qiuk[0]) * psc5 + (qkuip[0] - qiukp[0]) 
			    * dsc5) + gfri[4] * qir[0] + gfri[5] * qkr[0];
		    ftm2ri[1] = gfri[0] * yr + (-rr3 * ck * (uind_ref(2, i__) 
			    * psc3 + uinp_ref(2, i__) * dsc3) + rr5 * sc[3] * 
			    (uind_ref(2, i__) * psc5 + uinp_ref(2, i__) * 
			    dsc5) - rr7 * sc[5] * (uind_ref(2, i__) * psc7 + 
			    uinp_ref(2, i__) * dsc7)) * .5 + (rr3 * ci * (
			    uind_ref(2, k) * psc3 + uinp_ref(2, k) * dsc3) + 
			    rr5 * sc[2] * (uind_ref(2, k) * psc5 + uinp_ref(2,
			     k) * dsc5) + rr7 * sc[4] * (uind_ref(2, k) * 
			    psc7 + uinp_ref(2, k) * dsc7)) * .5 + rr5 * usc5 *
			     (sci[3] * uinp_ref(2, i__) + scip[3] * uind_ref(
			    2, i__) + sci[2] * uinp_ref(2, k) + scip[2] * 
			    uind_ref(2, k)) * .5 + (sci[3] * psc5 + scip[3] * 
			    dsc5) * .5 * rr5 * di[1] + (sci[2] * psc5 + scip[
			    2] * dsc5) * .5 * rr5 * dk[1] + gfri[3] * .5 * ((
			    qkui[1] - qiuk[1]) * psc5 + (qkuip[1] - qiukp[1]) 
			    * dsc5) + gfri[4] * qir[1] + gfri[5] * qkr[1];
		    ftm2ri[2] = gfri[0] * zr + (-rr3 * ck * (uind_ref(3, i__) 
			    * psc3 + uinp_ref(3, i__) * dsc3) + rr5 * sc[3] * 
			    (uind_ref(3, i__) * psc5 + uinp_ref(3, i__) * 
			    dsc5) - rr7 * sc[5] * (uind_ref(3, i__) * psc7 + 
			    uinp_ref(3, i__) * dsc7)) * .5 + (rr3 * ci * (
			    uind_ref(3, k) * psc3 + uinp_ref(3, k) * dsc3) + 
			    rr5 * sc[2] * (uind_ref(3, k) * psc5 + uinp_ref(3,
			     k) * dsc5) + rr7 * sc[4] * (uind_ref(3, k) * 
			    psc7 + uinp_ref(3, k) * dsc7)) * .5 + rr5 * usc5 *
			     (sci[3] * uinp_ref(3, i__) + scip[3] * uind_ref(
			    3, i__) + sci[2] * uinp_ref(3, k) + scip[2] * 
			    uind_ref(3, k)) * .5 + (sci[3] * psc5 + scip[3] * 
			    dsc5) * .5 * rr5 * di[2] + (sci[2] * psc5 + scip[
			    2] * dsc5) * .5 * rr5 * dk[2] + gfri[3] * .5 * ((
			    qkui[2] - qiuk[2]) * psc5 + (qkuip[2] - qiukp[2]) 
			    * dsc5) + gfri[4] * qir[2] + gfri[5] * qkr[2];
		}

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

/*     modify the forces for partially excluded interactions */

		ftm2i[0] = ftm2i[0] - fridmp[0] - findmp[0];
		ftm2i[1] = ftm2i[1] - fridmp[1] - findmp[1];
		ftm2i[2] = ftm2i[2] - fridmp[2] - findmp[2];

/*     correction to convert mutual to direct polarization force */

		if (s_cmp(polpot_1.poltyp, "DIRECT", (ftnlen)6, (ftnlen)6) == 
			0) {
		    gfd = (bn[2] * scip[1] - bn[3] * (scip[2] * sci[3] + sci[
			    2] * scip[3])) * .5;
		    gfdr = (rr5 * scip[1] * usc3 - rr7 * (scip[2] * sci[3] + 
			    sci[2] * scip[3]) * usc5) * .5;
		    ftm2i[0] = ftm2i[0] - gfd * xr - bn[2] * .5 * (sci[3] * 
			    uinp_ref(1, i__) + scip[3] * uind_ref(1, i__) + 
			    sci[2] * uinp_ref(1, k) + scip[2] * uind_ref(1, k)
			    );
		    ftm2i[1] = ftm2i[1] - gfd * yr - bn[2] * .5 * (sci[3] * 
			    uinp_ref(2, i__) + scip[3] * uind_ref(2, i__) + 
			    sci[2] * uinp_ref(2, k) + scip[2] * uind_ref(2, k)
			    );
		    ftm2i[2] = ftm2i[2] - gfd * zr - bn[2] * .5 * (sci[3] * 
			    uinp_ref(3, i__) + scip[3] * uind_ref(3, i__) + 
			    sci[2] * uinp_ref(3, k) + scip[2] * uind_ref(3, k)
			    );
		    fdir[0] = gfdr * xr + usc5 * .5 * rr5 * (sci[3] * 
			    uinp_ref(1, i__) + scip[3] * uind_ref(1, i__) + 
			    sci[2] * uinp_ref(1, k) + scip[2] * uind_ref(1, k)
			    );
		    fdir[1] = gfdr * yr + usc5 * .5 * rr5 * (sci[3] * 
			    uinp_ref(2, i__) + scip[3] * uind_ref(2, i__) + 
			    sci[2] * uinp_ref(2, k) + scip[2] * uind_ref(2, k)
			    );
		    fdir[2] = gfdr * zr + usc5 * .5 * rr5 * (sci[3] * 
			    uinp_ref(3, i__) + scip[3] * uind_ref(3, i__) + 
			    sci[2] * uinp_ref(3, k) + scip[2] * uind_ref(3, k)
			    );
		    ftm2i[0] = ftm2i[0] + fdir[0] + findmp[0];
		    ftm2i[1] = ftm2i[1] + fdir[1] + findmp[1];
		    ftm2i[2] = ftm2i[2] + fdir[2] + findmp[2];
		}

/*     get the permanent torque with screening */

		ttm2[0] = -bn[1] * dixdk[0] + gf[1] * dixr[0] + gf[3] * (
			dixqkr[0] + dkxqir[0] + rxqidk[0] - qixqk[0] * 2.) - 
			gf[4] * rxqir[0] - gf[6] * (rxqikr[0] + qkrxqir[0]);
		ttm2[1] = -bn[1] * dixdk[1] + gf[1] * dixr[1] + gf[3] * (
			dixqkr[1] + dkxqir[1] + rxqidk[1] - qixqk[1] * 2.) - 
			gf[4] * rxqir[1] - gf[6] * (rxqikr[1] + qkrxqir[1]);
		ttm2[2] = -bn[1] * dixdk[2] + gf[1] * dixr[2] + gf[3] * (
			dixqkr[2] + dkxqir[2] + rxqidk[2] - qixqk[2] * 2.) - 
			gf[4] * rxqir[2] - gf[6] * (rxqikr[2] + qkrxqir[2]);
		ttm3[0] = bn[1] * dixdk[0] + gf[2] * dkxr[0] - gf[3] * (
			dixqkr[0] + dkxqir[0] + rxqkdi[0] - qixqk[0] * 2.) - 
			gf[5] * rxqkr[0] - gf[6] * (rxqkir[0] - qkrxqir[0]);
		ttm3[1] = bn[1] * dixdk[1] + gf[2] * dkxr[1] - gf[3] * (
			dixqkr[1] + dkxqir[1] + rxqkdi[1] - qixqk[1] * 2.) - 
			gf[5] * rxqkr[1] - gf[6] * (rxqkir[1] - qkrxqir[1]);
		ttm3[2] = bn[1] * dixdk[2] + gf[2] * dkxr[2] - gf[3] * (
			dixqkr[2] + dkxqir[2] + rxqkdi[2] - qixqk[2] * 2.) - 
			gf[5] * rxqkr[2] - gf[6] * (rxqkir[2] - qkrxqir[2]);

/*     get the permanent torque without screening */

		if (dorl) {
		    ttm2r[0] = -rr3 * dixdk[0] + gfr[1] * dixr[0] + gfr[3] * (
			    dixqkr[0] + dkxqir[0] + rxqidk[0] - qixqk[0] * 2.)
			     - gfr[4] * rxqir[0] - gfr[6] * (rxqikr[0] + 
			    qkrxqir[0]);
		    ttm2r[1] = -rr3 * dixdk[1] + gfr[1] * dixr[1] + gfr[3] * (
			    dixqkr[1] + dkxqir[1] + rxqidk[1] - qixqk[1] * 2.)
			     - gfr[4] * rxqir[1] - gfr[6] * (rxqikr[1] + 
			    qkrxqir[1]);
		    ttm2r[2] = -rr3 * dixdk[2] + gfr[1] * dixr[2] + gfr[3] * (
			    dixqkr[2] + dkxqir[2] + rxqidk[2] - qixqk[2] * 2.)
			     - gfr[4] * rxqir[2] - gfr[6] * (rxqikr[2] + 
			    qkrxqir[2]);
		    ttm3r[0] = rr3 * dixdk[0] + gfr[2] * dkxr[0] - gfr[3] * (
			    dixqkr[0] + dkxqir[0] + rxqkdi[0] - qixqk[0] * 2.)
			     - gfr[5] * rxqkr[0] - gfr[6] * (rxqkir[0] - 
			    qkrxqir[0]);
		    ttm3r[1] = rr3 * dixdk[1] + gfr[2] * dkxr[1] - gfr[3] * (
			    dixqkr[1] + dkxqir[1] + rxqkdi[1] - qixqk[1] * 2.)
			     - gfr[5] * rxqkr[1] - gfr[6] * (rxqkir[1] - 
			    qkrxqir[1]);
		    ttm3r[2] = rr3 * dixdk[2] + gfr[2] * dkxr[2] - gfr[3] * (
			    dixqkr[2] + dkxqir[2] + rxqkdi[2] - qixqk[2] * 2.)
			     - gfr[5] * rxqkr[2] - gfr[6] * (rxqkir[2] - 
			    qkrxqir[2]);
		}

/*     get the induced torque with screening */

		gti[1] = bn[2] * .5 * (sci[3] + scip[3]);
		gti[2] = bn[2] * .5 * (sci[2] + scip[2]);
		gti[3] = gfi[3];
		gti[4] = gfi[4];
		gti[5] = gfi[5];
		ttm2i[0] = bn[1] * -.5 * (dixuk[0] + dixukp[0]) + gti[1] * 
			dixr[0] - gti[4] * rxqir[0] + gti[3] * .5 * (ukxqir[0]
			 + rxqiuk[0] + ukxqirp[0] + rxqiukp[0]);
		ttm2i[1] = bn[1] * -.5 * (dixuk[1] + dixukp[1]) + gti[1] * 
			dixr[1] - gti[4] * rxqir[1] + gti[3] * .5 * (ukxqir[1]
			 + rxqiuk[1] + ukxqirp[1] + rxqiukp[1]);
		ttm2i[2] = bn[1] * -.5 * (dixuk[2] + dixukp[2]) + gti[1] * 
			dixr[2] - gti[4] * rxqir[2] + gti[3] * .5 * (ukxqir[2]
			 + rxqiuk[2] + ukxqirp[2] + rxqiukp[2]);
		ttm3i[0] = bn[1] * -.5 * (dkxui[0] + dkxuip[0]) + gti[2] * 
			dkxr[0] - gti[5] * rxqkr[0] - gti[3] * .5 * (uixqkr[0]
			 + rxqkui[0] + uixqkrp[0] + rxqkuip[0]);
		ttm3i[1] = bn[1] * -.5 * (dkxui[1] + dkxuip[1]) + gti[2] * 
			dkxr[1] - gti[5] * rxqkr[1] - gti[3] * .5 * (uixqkr[1]
			 + rxqkui[1] + uixqkrp[1] + rxqkuip[1]);
		ttm3i[2] = bn[1] * -.5 * (dkxui[2] + dkxuip[2]) + gti[2] * 
			dkxr[2] - gti[5] * rxqkr[2] - gti[3] * .5 * (uixqkr[2]
			 + rxqkui[2] + uixqkrp[2] + rxqkuip[2]);

/*     get the induced torque without screening */

		if (dorli) {
		    gtri[1] = rr5 * .5 * (sci[3] * psc5 + scip[3] * dsc5);
		    gtri[2] = rr5 * .5 * (sci[2] * psc5 + scip[2] * dsc5);
		    gtri[3] = gfri[3];
		    gtri[4] = gfri[4];
		    gtri[5] = gfri[5];
		    ttm2ri[0] = -rr3 * (dixuk[0] * psc3 + dixukp[0] * dsc3) * 
			    .5 + gtri[1] * dixr[0] - gtri[4] * rxqir[0] + 
			    gtri[3] * ((ukxqir[0] + rxqiuk[0]) * psc5 + (
			    ukxqirp[0] + rxqiukp[0]) * dsc5) * .5;
		    ttm2ri[1] = -rr3 * (dixuk[1] * psc3 + dixukp[1] * dsc3) * 
			    .5 + gtri[1] * dixr[1] - gtri[4] * rxqir[1] + 
			    gtri[3] * ((ukxqir[1] + rxqiuk[1]) * psc5 + (
			    ukxqirp[1] + rxqiukp[1]) * dsc5) * .5;
		    ttm2ri[2] = -rr3 * (dixuk[2] * psc3 + dixukp[2] * dsc3) * 
			    .5 + gtri[1] * dixr[2] - gtri[4] * rxqir[2] + 
			    gtri[3] * ((ukxqir[2] + rxqiuk[2]) * psc5 + (
			    ukxqirp[2] + rxqiukp[2]) * dsc5) * .5;
		    ttm3ri[0] = -rr3 * (dkxui[0] * psc3 + dkxuip[0] * dsc3) * 
			    .5 + gtri[2] * dkxr[0] - gtri[5] * rxqkr[0] - 
			    gtri[3] * ((uixqkr[0] + rxqkui[0]) * psc5 + (
			    uixqkrp[0] + rxqkuip[0]) * dsc5) * .5;
		    ttm3ri[1] = -rr3 * (dkxui[1] * psc3 + dkxuip[1] * dsc3) * 
			    .5 + gtri[2] * dkxr[1] - gtri[5] * rxqkr[1] - 
			    gtri[3] * ((uixqkr[1] + rxqkui[1]) * psc5 + (
			    uixqkrp[1] + rxqkuip[1]) * dsc5) * .5;
		    ttm3ri[2] = -rr3 * (dkxui[2] * psc3 + dkxuip[2] * dsc3) * 
			    .5 + gtri[2] * dkxr[2] - gtri[5] * rxqkr[2] - 
			    gtri[3] * ((uixqkr[2] + rxqkui[2]) * psc5 + (
			    uixqkrp[2] + rxqkuip[2]) * dsc5) * .5;
		}

/*     handle the case where scaling is used */

		for (j = 1; j <= 3; ++j) {
		    ftm2[j - 1] = f * (ftm2[j - 1] - (1. - mscale[kk - 1]) * 
			    ftm2r[j - 1]);
		    ftm2i[j - 1] = f * (ftm2i[j - 1] - ftm2ri[j - 1]);
		    ttm2[j - 1] = f * (ttm2[j - 1] - (1. - mscale[kk - 1]) * 
			    ttm2r[j - 1]);
		    ttm2i[j - 1] = f * (ttm2i[j - 1] - ttm2ri[j - 1]);
		    ttm3[j - 1] = f * (ttm3[j - 1] - (1. - mscale[kk - 1]) * 
			    ttm3r[j - 1]);
		    ttm3i[j - 1] = f * (ttm3i[j - 1] - ttm3ri[j - 1]);
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

/*     increment the internal virial tensor components */

		iaz = mpole_1.zaxis[i__ - 1];
		iax = mpole_1.xaxis[i__ - 1];
		iay = mpole_1.yaxis[i__ - 1];
		kaz = mpole_1.zaxis[k - 1];
		kax = mpole_1.xaxis[k - 1];
		kay = mpole_1.yaxis[k - 1];
		if (iaz == 0) {
		    iaz = ii;
		}
		if (iax == 0) {
		    iax = ii;
		}
		if (iay == 0) {
		    iay = ii;
		}
		if (kaz == 0) {
		    kaz = kk;
		}
		if (kax == 0) {
		    kax = kk;
		}
		if (kay == 0) {
		    kay = kk;
		}
		xiz = atoms_1.x[iaz - 1] - atoms_1.x[ii - 1];
		yiz = atoms_1.y[iaz - 1] - atoms_1.y[ii - 1];
		ziz = atoms_1.z__[iaz - 1] - atoms_1.z__[ii - 1];
		xix = atoms_1.x[iax - 1] - atoms_1.x[ii - 1];
		yix = atoms_1.y[iax - 1] - atoms_1.y[ii - 1];
		zix = atoms_1.z__[iax - 1] - atoms_1.z__[ii - 1];
		xiy = atoms_1.x[iay - 1] - atoms_1.x[ii - 1];
		yiy = atoms_1.y[iay - 1] - atoms_1.y[ii - 1];
		ziy = atoms_1.z__[iay - 1] - atoms_1.z__[ii - 1];
		xkz = atoms_1.x[kaz - 1] - atoms_1.x[kk - 1];
		ykz = atoms_1.y[kaz - 1] - atoms_1.y[kk - 1];
		zkz = atoms_1.z__[kaz - 1] - atoms_1.z__[kk - 1];
		xkx = atoms_1.x[kax - 1] - atoms_1.x[kk - 1];
		ykx = atoms_1.y[kax - 1] - atoms_1.y[kk - 1];
		zkx = atoms_1.z__[kax - 1] - atoms_1.z__[kk - 1];
		xky = atoms_1.x[kay - 1] - atoms_1.x[kk - 1];
		yky = atoms_1.y[kay - 1] - atoms_1.y[kk - 1];
		zky = atoms_1.z__[kay - 1] - atoms_1.z__[kk - 1];
		vxx = -xr * (ftm2[0] + ftm2i[0]) + xix * frcxi[0] + xiy * 
			frcyi[0] + xiz * frczi[0] + xkx * frcxk[0] + xky * 
			frcyk[0] + xkz * frczk[0];
		vyx = -yr * (ftm2[0] + ftm2i[0]) + yix * frcxi[0] + yiy * 
			frcyi[0] + yiz * frczi[0] + ykx * frcxk[0] + yky * 
			frcyk[0] + ykz * frczk[0];
		vzx = -zr * (ftm2[0] + ftm2i[0]) + zix * frcxi[0] + ziy * 
			frcyi[0] + ziz * frczi[0] + zkx * frcxk[0] + zky * 
			frcyk[0] + zkz * frczk[0];
		vyy = -yr * (ftm2[1] + ftm2i[1]) + yix * frcxi[1] + yiy * 
			frcyi[1] + yiz * frczi[1] + ykx * frcxk[1] + yky * 
			frcyk[1] + ykz * frczk[1];
		vzy = -zr * (ftm2[1] + ftm2i[1]) + zix * frcxi[1] + ziy * 
			frcyi[1] + ziz * frczi[1] + zkx * frcxk[1] + zky * 
			frcyk[1] + zkz * frczk[1];
		vzz = -zr * (ftm2[2] + ftm2i[2]) + zix * frcxi[2] + ziy * 
			frcyi[2] + ziz * frczi[2] + zkx * frcxk[2] + zky * 
			frcyk[2] + zkz * frczk[2];
		vir_ref(1, 1) = vir_ref(1, 1) + vxx;
		vir_ref(2, 1) = vir_ref(2, 1) + vyx;
		vir_ref(3, 1) = vir_ref(3, 1) + vzx;
		vir_ref(1, 2) = vir_ref(1, 2) + vyx;
		vir_ref(2, 2) = vir_ref(2, 2) + vyy;
		vir_ref(3, 2) = vir_ref(3, 2) + vzy;
		vir_ref(1, 3) = vir_ref(1, 3) + vzx;
		vir_ref(2, 3) = vir_ref(2, 3) + vzy;
		vir_ref(3, 3) = vir_ref(3, 3) + vzz;
	    }
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
    return 0;
} /* ereal1d_ */

#undef rpole_ref
#undef uinp_ref
#undef elst_ref
#undef uind_ref
#undef vir_ref
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




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine emrecip1  --  mpole Ewald recip energy & derivs  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "emrecip1" evaluates the reciprocal space portion of the particle */
/*     mesh Ewald summation energy and gradient due to atomic multipole */
/*     interactions and dipole polarizability */

/*     literature reference: */

/*     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate */
/*     Representation of Electrostatics in Classical Force Fields: */
/*     Efficient Implementation of Multipolar Interactions in */
/*     Biomolecular Simulations", Journal of Chemical Physics, 120, */
/*     73-87 (2004) */

/*     modifications for nonperiodic systems suggested by Tom Darden */
/*     during May 2007 */


/* Subroutine */ int emrecip1_(void)
{
    /* Initialized data */

    static integer deriv1[10] = { 2,5,8,9,11,16,18,14,15,20 };
    static integer deriv2[10] = { 3,8,6,10,14,12,19,16,20,17 };
    static integer deriv3[10] = { 4,9,10,7,15,17,13,20,18,19 };

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal), cos(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int grid_mpole__(doublereal *), fphi_mpole__(
	    doublereal *), cmp_to_fmp__(doublereal *, doublereal *);
    static doublereal a[9]	/* was [3][3] */, e;
    static integer i__, j, k;
    static doublereal f1, f2, h1, h2;
    static integer j1, j2, j3, k1, k2, k3, m1, m2, m3;
    static doublereal r1, r2, r3, h3, f3;
    static integer ii;
    extern /* Subroutine */ int frac_to_cart__(doublereal *), fphi_to_cphi__(
	    doublereal *, doublereal *);
    static integer nf1, nf2, nf3;
    extern /* Subroutine */ int bspline_fill__(void);
    static integer nff;
    static doublereal ftc[100]	/* was [10][10] */, frc[75000]	/* was [3][
	    25000] */, cmp[250000]	/* was [10][25000] */, fmp[250000]	
	    /* was [10][25000] */, hsq, trq[75000]	/* was [3][25000] */, 
	    vxx, vyx, vzx, vyy, vzy, vzz, cphi[250000]	/* was [10][25000] */,
	     fphi[500000]	/* was [20][25000] */, term;
    static integer ntot;
    static doublereal cphid[4], fphid[250000]	/* was [10][25000] */, cphim[
	    4], denom, cphip[4], fuind[75000]	/* was [3][25000] */, fphip[
	    250000]	/* was [10][25000] */, eterm, fuinp[75000]	/* 
	    was [3][25000] */, qgrip[2000000]	/* was [2][100][100][100] */, 
	    pterm, vterm, struc2, fphidp[500000]	/* was [20][25000] */;
    extern /* Subroutine */ int fftback_(void), torque2_(doublereal *, 
	    doublereal *);
    static doublereal expterm, volterm;
    extern /* Subroutine */ int fftfront_(void), grid_uind__(doublereal *, 
	    doublereal *), fphi_uind__(doublereal *, doublereal *, doublereal 
	    *);


#define a_ref(a_1,a_2) a[(a_2)*3 + a_1 - 4]
#define dem_ref(a_1,a_2) deriv_1.dem[(a_2)*3 + a_1 - 4]
#define dep_ref(a_1,a_2) deriv_1.dep[(a_2)*3 + a_1 - 4]
#define ftc_ref(a_1,a_2) ftc[(a_2)*10 + a_1 - 11]
#define frc_ref(a_1,a_2) frc[(a_2)*3 + a_1 - 4]
#define cmp_ref(a_1,a_2) cmp[(a_2)*10 + a_1 - 11]
#define fmp_ref(a_1,a_2) fmp[(a_2)*10 + a_1 - 11]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define trq_ref(a_1,a_2) trq[(a_2)*3 + a_1 - 4]
#define qfac_ref(a_1,a_2,a_3) pme_1.qfac[((a_3)*100 + (a_2))*100 + a_1 - \
10101]
#define cphi_ref(a_1,a_2) cphi[(a_2)*10 + a_1 - 11]
#define fphi_ref(a_1,a_2) fphi[(a_2)*20 + a_1 - 21]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define uinp_ref(a_1,a_2) polar_1.uinp[(a_2)*3 + a_1 - 4]
#define fphid_ref(a_1,a_2) fphid[(a_2)*10 + a_1 - 11]
#define recip_ref(a_1,a_2) boxes_1.recip[(a_2)*3 + a_1 - 4]
#define fuind_ref(a_1,a_2) fuind[(a_2)*3 + a_1 - 4]
#define qgrid_ref(a_1,a_2,a_3,a_4) pme_1.qgrid[(((a_4)*100 + (a_3))*100 + \
(a_2))*2 + a_1 - 20203]
#define fphip_ref(a_1,a_2) fphip[(a_2)*10 + a_1 - 11]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define fuinp_ref(a_1,a_2) fuinp[(a_2)*3 + a_1 - 4]
#define qgrip_ref(a_1,a_2,a_3,a_4) qgrip[(((a_4)*100 + (a_3))*100 + (a_2))\
*2 + a_1 - 20203]
#define fphidp_ref(a_1,a_2) fphidp[(a_2)*20 + a_1 - 21]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  ewald.i  --  parameters and options for Ewald summation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     aewald     Ewald convergence coefficient value (Ang-1) */
/*     boundary   Ewald boundary condition; none, tinfoil or vacuum */




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

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  virial.i  --  components of internal virial tensor  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     vir    total internal virial Cartesian tensor components */



/*     derivative indices into the fphi and fdip_phi arrays */



/*     return if the Ewald coefficient is zero */

    if (ewald_1.aewald < 1e-6) {
	return 0;
    }

/*     zero out the temporary virial accumulation variables */

    vxx = 0.;
    vyx = 0.;
    vzx = 0.;
    vyy = 0.;
    vzy = 0.;
    vzz = 0.;

/*     copy multipole moments and coordinates to local storage */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cmp_ref(1, i__) = rpole_ref(1, i__);
	cmp_ref(2, i__) = rpole_ref(2, i__);
	cmp_ref(3, i__) = rpole_ref(3, i__);
	cmp_ref(4, i__) = rpole_ref(4, i__);
	cmp_ref(5, i__) = rpole_ref(5, i__);
	cmp_ref(6, i__) = rpole_ref(9, i__);
	cmp_ref(7, i__) = rpole_ref(13, i__);
	cmp_ref(8, i__) = rpole_ref(6, i__) * 2.;
	cmp_ref(9, i__) = rpole_ref(7, i__) * 2.;
	cmp_ref(10, i__) = rpole_ref(10, i__) * 2.;
    }

/*     get the fractional to Cartesian transformation matrix */

    frac_to_cart__(ftc);

/*     compute the arrays of B-spline coefficients */

    if (! potent_1.use_polar__) {
	bspline_fill__();
    }

/*     assign permanent and induced multipoles to PME grid */
/*     and perform the 3-D FFT forward transformation */

    if (potent_1.use_polar__) {
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 2; j <= 4; ++j) {
		cmp_ref(j, i__) = cmp_ref(j, i__) + uinp_ref(j - 1, i__);
	    }
	}
	cmp_to_fmp__(cmp, fmp);
	grid_mpole__(fmp);
	fftfront_();
	i__1 = pme_1.nfft3;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = pme_1.nfft2;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = pme_1.nfft1;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    qgrip_ref(1, i__, j, k) = qgrid_ref(1, i__, j, k);
		    qgrip_ref(2, i__, j, k) = qgrid_ref(2, i__, j, k);
		}
	    }
	}
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 2; j <= 4; ++j) {
		cmp_ref(j, i__) = cmp_ref(j, i__) + uind_ref(j - 1, i__) - 
			uinp_ref(j - 1, i__);
	    }
	}
	cmp_to_fmp__(cmp, fmp);
	grid_mpole__(fmp);
	fftfront_();
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 2; j <= 4; ++j) {
		cmp_ref(j, i__) = cmp_ref(j, i__) - uind_ref(j - 1, i__);
	    }
	}
    } else {
	cmp_to_fmp__(cmp, fmp);
	grid_mpole__(fmp);
	fftfront_();
	i__1 = pme_1.nfft3;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = pme_1.nfft2;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = pme_1.nfft1;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    qgrip_ref(1, i__, j, k) = qgrid_ref(1, i__, j, k);
		    qgrip_ref(2, i__, j, k) = qgrid_ref(2, i__, j, k);
		}
	    }
	}
    }

/*     make the scalar summation over reciprocal lattice */

    ntot = pme_1.nfft1 * pme_1.nfft2 * pme_1.nfft3;
/* Computing 2nd power */
    d__1 = 3.141592653589793238 / ewald_1.aewald;
    pterm = d__1 * d__1;
    volterm = boxes_1.volbox * 3.141592653589793238;
    nff = pme_1.nfft1 * pme_1.nfft2;
    nf1 = (pme_1.nfft1 + 1) / 2;
    nf2 = (pme_1.nfft2 + 1) / 2;
    nf3 = (pme_1.nfft3 + 1) / 2;
    i__1 = ntot - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k3 = i__ / nff + 1;
	j = i__ - (k3 - 1) * nff;
	k2 = j / pme_1.nfft1 + 1;
	k1 = j - (k2 - 1) * pme_1.nfft1 + 1;
	m1 = k1 - 1;
	m2 = k2 - 1;
	m3 = k3 - 1;
	if (k1 > nf1) {
	    m1 -= pme_1.nfft1;
	}
	if (k2 > nf2) {
	    m2 -= pme_1.nfft2;
	}
	if (k3 > nf3) {
	    m3 -= pme_1.nfft3;
	}
	r1 = (doublereal) m1;
	r2 = (doublereal) m2;
	r3 = (doublereal) m3;
	h1 = recip_ref(1, 1) * r1 + recip_ref(1, 2) * r2 + recip_ref(1, 3) * 
		r3;
	h2 = recip_ref(2, 1) * r1 + recip_ref(2, 2) * r2 + recip_ref(2, 3) * 
		r3;
	h3 = recip_ref(3, 1) * r1 + recip_ref(3, 2) * r2 + recip_ref(3, 3) * 
		r3;
	hsq = h1 * h1 + h2 * h2 + h3 * h3;
	term = -pterm * hsq;
	expterm = 0.;
	if (term > -50.) {
	    denom = volterm * hsq * pme_1.bsmod1[k1 - 1] * pme_1.bsmod2[k2 - 
		    1] * pme_1.bsmod3[k3 - 1];
	    expterm = exp(term) / denom;
	    if (! bound_1.use_bounds__) {
		expterm *= 1. - cos(boxes_1.xbox * 3.141592653589793238 * 
			sqrt(hsq));
	    } else if (boxes_1.octahedron) {
		if ((m1 + m2 + m3) % 2 != 0) {
		    expterm = 0.;
		}
	    }
	    struc2 = qgrid_ref(1, k1, k2, k3) * qgrip_ref(1, k1, k2, k3) + 
		    qgrid_ref(2, k1, k2, k3) * qgrip_ref(2, k1, k2, k3);
	    eterm = chgpot_1.electric * .5 * expterm * struc2;
	    vterm = 2. / hsq * (1. - term) * eterm;
	    vxx = vxx + h1 * h1 * vterm - eterm;
	    vyx += h2 * h1 * vterm;
	    vzx += h3 * h1 * vterm;
	    vyy = vyy + h2 * h2 * vterm - eterm;
	    vzy += h3 * h2 * vterm;
	    vzz = vzz + h3 * h3 * vterm - eterm;
	}
	qfac_ref(k1, k2, k3) = expterm;
    }

/*     transform permanent multipoles without induced dipoles */

    if (potent_1.use_polar__) {
	cmp_to_fmp__(cmp, fmp);
	grid_mpole__(fmp);
	fftfront_();
    }

/*     account for the zeroth grid point for a finite system */

    qfac_ref(1, 1, 1) = 0.;
    if (! bound_1.use_bounds__) {
	expterm = 1.5707963267948966 / boxes_1.xbox;
/* Computing 2nd power */
	d__1 = qgrid_ref(1, 1, 1, 1);
/* Computing 2nd power */
	d__2 = qgrid_ref(2, 1, 1, 1);
	struc2 = d__1 * d__1 + d__2 * d__2;
	e = expterm * .5 * struc2;
	qfac_ref(1, 1, 1) = expterm;
    }

/*     complete the transformation of the PME grid */

    i__1 = pme_1.nfft3;
    for (k = 1; k <= i__1; ++k) {
	i__2 = pme_1.nfft2;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = pme_1.nfft1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		term = qfac_ref(i__, j, k);
		qgrid_ref(1, i__, j, k) = term * qgrid_ref(1, i__, j, k);
		qgrid_ref(2, i__, j, k) = term * qgrid_ref(2, i__, j, k);
	    }
	}
    }

/*     perform 3-D FFT backward transform and get potential */

    fftback_();
    fphi_mpole__(fphi);
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 20; ++j) {
	    fphi_ref(j, i__) = chgpot_1.electric * fphi_ref(j, i__);
	}
    }
    fphi_to_cphi__(fphi, cphi);

/*     increment the permanent multipole energy and gradient */

    e = 0.;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f1 = 0.;
	f2 = 0.;
	f3 = 0.;
	for (k = 1; k <= 10; ++k) {
	    e += fmp_ref(k, i__) * fphi_ref(k, i__);
	    f1 += fmp_ref(k, i__) * fphi_ref(deriv1[k - 1], i__);
	    f2 += fmp_ref(k, i__) * fphi_ref(deriv2[k - 1], i__);
	    f3 += fmp_ref(k, i__) * fphi_ref(deriv3[k - 1], i__);
	}
	f1 = (doublereal) pme_1.nfft1 * f1;
	f2 = (doublereal) pme_1.nfft2 * f2;
	f3 = (doublereal) pme_1.nfft3 * f3;
	frc_ref(1, i__) = recip_ref(1, 1) * f1 + recip_ref(1, 2) * f2 + 
		recip_ref(1, 3) * f3;
	frc_ref(2, i__) = recip_ref(2, 1) * f1 + recip_ref(2, 2) * f2 + 
		recip_ref(2, 3) * f3;
	frc_ref(3, i__) = recip_ref(3, 1) * f1 + recip_ref(3, 2) * f2 + 
		recip_ref(3, 3) * f3;
    }
    e *= .5;
    energi_1.em += e;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	dem_ref(1, ii) = dem_ref(1, ii) + frc_ref(1, i__);
	dem_ref(2, ii) = dem_ref(2, ii) + frc_ref(2, i__);
	dem_ref(3, ii) = dem_ref(3, ii) + frc_ref(3, i__);
    }

/*     distribute torques into the permanent multipole gradient */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	trq_ref(1, i__) = cmp_ref(4, i__) * cphi_ref(3, i__) - cmp_ref(3, i__)
		 * cphi_ref(4, i__) + (cmp_ref(7, i__) - cmp_ref(6, i__)) * 
		2. * cphi_ref(10, i__) + cmp_ref(9, i__) * cphi_ref(8, i__) + 
		cmp_ref(10, i__) * cphi_ref(6, i__) - cmp_ref(8, i__) * 
		cphi_ref(9, i__) - cmp_ref(10, i__) * cphi_ref(7, i__);
	trq_ref(2, i__) = cmp_ref(2, i__) * cphi_ref(4, i__) - cmp_ref(4, i__)
		 * cphi_ref(2, i__) + (cmp_ref(5, i__) - cmp_ref(7, i__)) * 
		2. * cphi_ref(9, i__) + cmp_ref(8, i__) * cphi_ref(10, i__) + 
		cmp_ref(9, i__) * cphi_ref(7, i__) - cmp_ref(9, i__) * 
		cphi_ref(5, i__) - cmp_ref(10, i__) * cphi_ref(8, i__);
	trq_ref(3, i__) = cmp_ref(3, i__) * cphi_ref(2, i__) - cmp_ref(2, i__)
		 * cphi_ref(3, i__) + (cmp_ref(6, i__) - cmp_ref(5, i__)) * 
		2. * cphi_ref(8, i__) + cmp_ref(8, i__) * cphi_ref(5, i__) + 
		cmp_ref(10, i__) * cphi_ref(9, i__) - cmp_ref(8, i__) * 
		cphi_ref(6, i__) - cmp_ref(9, i__) * cphi_ref(10, i__);
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	frc_ref(1, i__) = 0.;
	frc_ref(2, i__) = 0.;
	frc_ref(3, i__) = 0.;
    }
    torque2_(trq, frc);
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dem_ref(1, i__) = dem_ref(1, i__) + frc_ref(1, i__);
	dem_ref(2, i__) = dem_ref(2, i__) + frc_ref(2, i__);
	dem_ref(3, i__) = dem_ref(3, i__) + frc_ref(3, i__);
    }

/*     permanent torque contribution to the internal virial */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vxx = vxx - cmp_ref(2, i__) * cphi_ref(2, i__) - cmp_ref(5, i__) * 2. 
		* cphi_ref(5, i__) - cmp_ref(8, i__) * cphi_ref(8, i__) - 
		cmp_ref(9, i__) * cphi_ref(9, i__);
	vyx = vyx - (cmp_ref(3, i__) * cphi_ref(2, i__) + cmp_ref(2, i__) * 
		cphi_ref(3, i__)) * .5 - (cmp_ref(5, i__) + cmp_ref(6, i__)) *
		 cphi_ref(8, i__) - cmp_ref(8, i__) * .5 * (cphi_ref(5, i__) 
		+ cphi_ref(6, i__)) - (cmp_ref(9, i__) * cphi_ref(10, i__) + 
		cmp_ref(10, i__) * cphi_ref(9, i__)) * .5;
	vzx = vzx - (cmp_ref(4, i__) * cphi_ref(2, i__) + cmp_ref(2, i__) * 
		cphi_ref(4, i__)) * .5 - (cmp_ref(5, i__) + cmp_ref(7, i__)) *
		 cphi_ref(9, i__) - cmp_ref(9, i__) * .5 * (cphi_ref(5, i__) 
		+ cphi_ref(7, i__)) - (cmp_ref(8, i__) * cphi_ref(10, i__) + 
		cmp_ref(10, i__) * cphi_ref(8, i__)) * .5;
	vyy = vyy - cmp_ref(3, i__) * cphi_ref(3, i__) - cmp_ref(6, i__) * 2. 
		* cphi_ref(6, i__) - cmp_ref(8, i__) * cphi_ref(8, i__) - 
		cmp_ref(10, i__) * cphi_ref(10, i__);
	vzy = vzy - (cmp_ref(4, i__) * cphi_ref(3, i__) + cmp_ref(3, i__) * 
		cphi_ref(4, i__)) * .5 - (cmp_ref(6, i__) + cmp_ref(7, i__)) *
		 cphi_ref(10, i__) - cmp_ref(10, i__) * .5 * (cphi_ref(6, i__)
		 + cphi_ref(7, i__)) - (cmp_ref(8, i__) * cphi_ref(9, i__) + 
		cmp_ref(9, i__) * cphi_ref(8, i__)) * .5;
	vzz = vzz - cmp_ref(4, i__) * cphi_ref(4, i__) - cmp_ref(7, i__) * 2. 
		* cphi_ref(7, i__) - cmp_ref(9, i__) * cphi_ref(9, i__) - 
		cmp_ref(10, i__) * cphi_ref(10, i__);
    }

/*     convert Cartesian induced dipoles to fractional coordinates */

    if (potent_1.use_polar__) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    a_ref(1, i__) = (doublereal) pme_1.nfft1 * recip_ref(i__, 1);
	    a_ref(2, i__) = (doublereal) pme_1.nfft2 * recip_ref(i__, 2);
	    a_ref(3, i__) = (doublereal) pme_1.nfft3 * recip_ref(i__, 3);
	}
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		fuind_ref(j, i__) = a_ref(j, 1) * uind_ref(1, i__) + a_ref(j, 
			2) * uind_ref(2, i__) + a_ref(j, 3) * uind_ref(3, i__)
			;
		fuinp_ref(j, i__) = a_ref(j, 1) * uinp_ref(1, i__) + a_ref(j, 
			2) * uinp_ref(2, i__) + a_ref(j, 3) * uinp_ref(3, i__)
			;
	    }
	}

/*     assign PME grid and perform 3-D FFT forward transform */

	grid_uind__(fuind, fuinp);
	fftfront_();

/*     account for the zeroth grid point for a finite system */

	if (! bound_1.use_bounds__) {
	    expterm = 1.5707963267948966 / boxes_1.xbox;
/* Computing 2nd power */
	    d__1 = qgrid_ref(1, 1, 1, 1);
/* Computing 2nd power */
	    d__2 = qgrid_ref(2, 1, 1, 1);
	    struc2 = d__1 * d__1 + d__2 * d__2;
	    e = expterm * .5 * struc2;
	}

/*     complete the transformation of the PME grid */

	i__1 = pme_1.nfft3;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = pme_1.nfft2;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = pme_1.nfft1;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    term = qfac_ref(i__, j, k);
		    qgrid_ref(1, i__, j, k) = term * qgrid_ref(1, i__, j, k);
		    qgrid_ref(2, i__, j, k) = term * qgrid_ref(2, i__, j, k);
		}
	    }
	}

/*     perform 3-D FFT backward transform and get potential */

	fftback_();
	fphi_uind__(fphid, fphip, fphidp);
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 10; ++j) {
		fphid_ref(j, i__) = chgpot_1.electric * fphid_ref(j, i__);
		fphip_ref(j, i__) = chgpot_1.electric * fphip_ref(j, i__);
	    }
	    for (j = 1; j <= 20; ++j) {
		fphidp_ref(j, i__) = chgpot_1.electric * fphidp_ref(j, i__);
	    }
	}

/*     increment the induced dipole energy and gradient */

	e = 0.;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    f1 = 0.;
	    f2 = 0.;
	    f3 = 0.;
	    for (k = 1; k <= 3; ++k) {
		j1 = deriv1[k];
		j2 = deriv2[k];
		j3 = deriv3[k];
		e += fuind_ref(k, i__) * fphi_ref(k + 1, i__);
		f1 = f1 + (fuind_ref(k, i__) + fuinp_ref(k, i__)) * fphi_ref(
			j1, i__) + fuind_ref(k, i__) * fphip_ref(j1, i__) + 
			fuinp_ref(k, i__) * fphid_ref(j1, i__);
		f2 = f2 + (fuind_ref(k, i__) + fuinp_ref(k, i__)) * fphi_ref(
			j2, i__) + fuind_ref(k, i__) * fphip_ref(j2, i__) + 
			fuinp_ref(k, i__) * fphid_ref(j2, i__);
		f3 = f3 + (fuind_ref(k, i__) + fuinp_ref(k, i__)) * fphi_ref(
			j3, i__) + fuind_ref(k, i__) * fphip_ref(j3, i__) + 
			fuinp_ref(k, i__) * fphid_ref(j3, i__);
		if (s_cmp(polpot_1.poltyp, "DIRECT", (ftnlen)6, (ftnlen)6) == 
			0) {
		    f1 = f1 - fuind_ref(k, i__) * fphip_ref(j1, i__) - 
			    fuinp_ref(k, i__) * fphid_ref(j1, i__);
		    f2 = f2 - fuind_ref(k, i__) * fphip_ref(j2, i__) - 
			    fuinp_ref(k, i__) * fphid_ref(j2, i__);
		    f3 = f3 - fuind_ref(k, i__) * fphip_ref(j3, i__) - 
			    fuinp_ref(k, i__) * fphid_ref(j3, i__);
		}
	    }
	    for (k = 1; k <= 10; ++k) {
		f1 += fmp_ref(k, i__) * fphidp_ref(deriv1[k - 1], i__);
		f2 += fmp_ref(k, i__) * fphidp_ref(deriv2[k - 1], i__);
		f3 += fmp_ref(k, i__) * fphidp_ref(deriv3[k - 1], i__);
	    }
	    f1 = (doublereal) pme_1.nfft1 * .5 * f1;
	    f2 = (doublereal) pme_1.nfft2 * .5 * f2;
	    f3 = (doublereal) pme_1.nfft3 * .5 * f3;
	    frc_ref(1, i__) = recip_ref(1, 1) * f1 + recip_ref(1, 2) * f2 + 
		    recip_ref(1, 3) * f3;
	    frc_ref(2, i__) = recip_ref(2, 1) * f1 + recip_ref(2, 2) * f2 + 
		    recip_ref(2, 3) * f3;
	    frc_ref(3, i__) = recip_ref(3, 1) * f1 + recip_ref(3, 2) * f2 + 
		    recip_ref(3, 3) * f3;
	}
	e *= .5;
	energi_1.ep += e;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    dep_ref(1, ii) = dep_ref(1, ii) + frc_ref(1, i__);
	    dep_ref(2, ii) = dep_ref(2, ii) + frc_ref(2, i__);
	    dep_ref(3, ii) = dep_ref(3, ii) + frc_ref(3, i__);
	}

/*     set the potential to be the induced dipole average */

	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (k = 1; k <= 10; ++k) {
		fphidp_ref(k, i__) = fphidp_ref(k, i__) * .5;
	    }
	}
	fphi_to_cphi__(fphidp, cphi);

/*     distribute torques into the induced dipole gradient */

	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    trq_ref(1, i__) = cmp_ref(4, i__) * cphi_ref(3, i__) - cmp_ref(3, 
		    i__) * cphi_ref(4, i__) + (cmp_ref(7, i__) - cmp_ref(6, 
		    i__)) * 2. * cphi_ref(10, i__) + cmp_ref(9, i__) * 
		    cphi_ref(8, i__) + cmp_ref(10, i__) * cphi_ref(6, i__) - 
		    cmp_ref(8, i__) * cphi_ref(9, i__) - cmp_ref(10, i__) * 
		    cphi_ref(7, i__);
	    trq_ref(2, i__) = cmp_ref(2, i__) * cphi_ref(4, i__) - cmp_ref(4, 
		    i__) * cphi_ref(2, i__) + (cmp_ref(5, i__) - cmp_ref(7, 
		    i__)) * 2. * cphi_ref(9, i__) + cmp_ref(8, i__) * 
		    cphi_ref(10, i__) + cmp_ref(9, i__) * cphi_ref(7, i__) - 
		    cmp_ref(9, i__) * cphi_ref(5, i__) - cmp_ref(10, i__) * 
		    cphi_ref(8, i__);
	    trq_ref(3, i__) = cmp_ref(3, i__) * cphi_ref(2, i__) - cmp_ref(2, 
		    i__) * cphi_ref(3, i__) + (cmp_ref(6, i__) - cmp_ref(5, 
		    i__)) * 2. * cphi_ref(8, i__) + cmp_ref(8, i__) * 
		    cphi_ref(5, i__) + cmp_ref(10, i__) * cphi_ref(9, i__) - 
		    cmp_ref(8, i__) * cphi_ref(6, i__) - cmp_ref(9, i__) * 
		    cphi_ref(10, i__);
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    frc_ref(1, i__) = 0.;
	    frc_ref(2, i__) = 0.;
	    frc_ref(3, i__) = 0.;
	}
	torque2_(trq, frc);
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dep_ref(1, i__) = dep_ref(1, i__) + frc_ref(1, i__);
	    dep_ref(2, i__) = dep_ref(2, i__) + frc_ref(2, i__);
	    dep_ref(3, i__) = dep_ref(3, i__) + frc_ref(3, i__);
	}

/*     induced torque contribution to the internal virial */

	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 2; j <= 4; ++j) {
		cphim[j - 1] = 0.;
		cphid[j - 1] = 0.;
		cphip[j - 1] = 0.;
		for (k = 2; k <= 4; ++k) {
		    cphim[j - 1] += ftc_ref(j, k) * fphi_ref(k, i__);
		    cphid[j - 1] += ftc_ref(j, k) * fphid_ref(k, i__);
		    cphip[j - 1] += ftc_ref(j, k) * fphip_ref(k, i__);
		}
	    }
	    vxx = vxx - cphi_ref(2, i__) * cmp_ref(2, i__) - (cphim[1] * (
		    uind_ref(1, i__) + uinp_ref(1, i__)) + cphid[1] * 
		    uinp_ref(1, i__) + cphip[1] * uind_ref(1, i__)) * .5;
	    vyx = vyx - (cphi_ref(2, i__) * cmp_ref(3, i__) + cphi_ref(3, i__)
		     * cmp_ref(2, i__)) * .5 - (cphim[1] * (uind_ref(2, i__) 
		    + uinp_ref(2, i__)) + cphim[2] * (uind_ref(1, i__) + 
		    uinp_ref(1, i__)) + cphid[1] * uinp_ref(2, i__) + cphip[1]
		     * uind_ref(2, i__) + cphid[2] * uinp_ref(1, i__) + cphip[
		    2] * uind_ref(1, i__)) * .25;
	    vzx = vzx - (cphi_ref(2, i__) * cmp_ref(4, i__) + cphi_ref(4, i__)
		     * cmp_ref(2, i__)) * .5 - (cphim[1] * (uind_ref(3, i__) 
		    + uinp_ref(3, i__)) + cphim[3] * (uind_ref(1, i__) + 
		    uinp_ref(1, i__)) + cphid[1] * uinp_ref(3, i__) + cphip[1]
		     * uind_ref(3, i__) + cphid[3] * uinp_ref(1, i__) + cphip[
		    3] * uind_ref(1, i__)) * .25;
	    vyy = vyy - cphi_ref(3, i__) * cmp_ref(3, i__) - (cphim[2] * (
		    uind_ref(2, i__) + uinp_ref(2, i__)) + cphid[2] * 
		    uinp_ref(2, i__) + cphip[2] * uind_ref(2, i__)) * .5;
	    vzy = vzy - (cphi_ref(3, i__) * cmp_ref(4, i__) + cphi_ref(4, i__)
		     * cmp_ref(3, i__)) * .5 - (cphim[2] * (uind_ref(3, i__) 
		    + uinp_ref(3, i__)) + cphim[3] * (uind_ref(2, i__) + 
		    uinp_ref(2, i__)) + cphid[2] * uinp_ref(3, i__) + cphip[2]
		     * uind_ref(3, i__) + cphid[3] * uinp_ref(2, i__) + cphip[
		    3] * uind_ref(2, i__)) * .25;
	    vzz = vzz - cphi_ref(4, i__) * cmp_ref(4, i__) - (cphim[3] * (
		    uind_ref(3, i__) + uinp_ref(3, i__)) + cphid[3] * 
		    uinp_ref(3, i__) + cphip[3] * uind_ref(3, i__)) * .5;
	    vxx = vxx - cmp_ref(5, i__) * 2. * cphi_ref(5, i__) - cmp_ref(8, 
		    i__) * cphi_ref(8, i__) - cmp_ref(9, i__) * cphi_ref(9, 
		    i__);
	    vyx = vyx - (cmp_ref(5, i__) + cmp_ref(6, i__)) * cphi_ref(8, i__)
		     - (cmp_ref(8, i__) * (cphi_ref(6, i__) + cphi_ref(5, i__)
		    ) + cmp_ref(9, i__) * cphi_ref(10, i__) + cmp_ref(10, i__)
		     * cphi_ref(9, i__)) * .5;
	    vzx = vzx - (cmp_ref(5, i__) + cmp_ref(7, i__)) * cphi_ref(9, i__)
		     - (cmp_ref(9, i__) * (cphi_ref(5, i__) + cphi_ref(7, i__)
		    ) + cmp_ref(8, i__) * cphi_ref(10, i__) + cmp_ref(10, i__)
		     * cphi_ref(8, i__)) * .5;
	    vyy = vyy - cmp_ref(6, i__) * 2. * cphi_ref(6, i__) - cmp_ref(8, 
		    i__) * cphi_ref(8, i__) - cmp_ref(10, i__) * cphi_ref(10, 
		    i__);
	    vzy = vzy - (cmp_ref(6, i__) + cmp_ref(7, i__)) * cphi_ref(10, 
		    i__) - (cmp_ref(10, i__) * (cphi_ref(6, i__) + cphi_ref(7,
		     i__)) + cmp_ref(8, i__) * cphi_ref(9, i__) + cmp_ref(9, 
		    i__) * cphi_ref(8, i__)) * .5;
	    vzz = vzz - cmp_ref(7, i__) * 2. * cphi_ref(7, i__) - cmp_ref(9, 
		    i__) * cphi_ref(9, i__) - cmp_ref(10, i__) * cphi_ref(10, 
		    i__);
	}
    }

/*     increment the internal virial tensor components */

    vir_ref(1, 1) = vir_ref(1, 1) + vxx;
    vir_ref(2, 1) = vir_ref(2, 1) + vyx;
    vir_ref(3, 1) = vir_ref(3, 1) + vzx;
    vir_ref(1, 2) = vir_ref(1, 2) + vyx;
    vir_ref(2, 2) = vir_ref(2, 2) + vyy;
    vir_ref(3, 2) = vir_ref(3, 2) + vzy;
    vir_ref(1, 3) = vir_ref(1, 3) + vzx;
    vir_ref(2, 3) = vir_ref(2, 3) + vzy;
    vir_ref(3, 3) = vir_ref(3, 3) + vzz;
    return 0;
} /* emrecip1_ */

#undef fphidp_ref
#undef qgrip_ref
#undef fuinp_ref
#undef rpole_ref
#undef fphip_ref
#undef qgrid_ref
#undef fuind_ref
#undef recip_ref
#undef fphid_ref
#undef uinp_ref
#undef uind_ref
#undef fphi_ref
#undef cphi_ref
#undef qfac_ref
#undef trq_ref
#undef vir_ref
#undef fmp_ref
#undef cmp_ref
#undef frc_ref
#undef ftc_ref
#undef dep_ref
#undef dem_ref
#undef a_ref


