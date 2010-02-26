/* echarge.f -- translated by f2c (version 20050501).
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
    doublereal m2[25000], deform, difft, diffv, diffc;
    logical use_smooth__, use_dem__, use_gda__, use_tophat__, use_stophat__;
} warp_;

#define warp_1 warp_

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
    doublereal xcell, ycell, zcell, xcell2, ycell2, zcell2;
    integer ncell, icell[30000]	/* was [3][10000] */;
} cell_;

#define cell_1 cell_

struct {
    doublereal pchg[25000];
    integer nion, iion[25000], jion[25000], kion[25000], chglist[25000];
} charge_;

#define charge_1 charge_

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
    doublereal esum, eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, 
	    ebt, ett, ev, ec, ecd, ed, em, ep, er, es, elf, eg, ex;
} energi_;

#define energi_1 energi_

struct {
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

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
    doublereal xbox, ybox, zbox, alpha, beta, gamma, xbox2, ybox2, zbox2, 
	    box34, lvec[9]	/* was [3][3] */, recip[9]	/* was [3][3] 
	    */, volbox, beta_sin__, beta_cos__, gamma_sin__, gamma_cos__, 
	    beta_term__, gamma_term__;
    logical orthogonal, monoclinic, triclinic, octahedron;
    char spacegrp[10];
} boxes_;

#define boxes_1 boxes_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nlight, kbx[25000], kby[25000], kbz[25000], kex[25000], key[25000]
	    , kez[25000], locx[200000], locy[200000], locz[200000], rgx[
	    200000], rgy[200000], rgz[200000];
} light_;

#define light_1 light_

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



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine echarge  --  charge-charge potential energy  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "echarge" calculates the charge-charge interaction energy */


/* Subroutine */ int echarge_(void)
{
    extern /* Subroutine */ int echarge0a_(void), echarge0b_(void), 
	    echarge0c_(void), echarge0d_(void), echarge0e_(void), echarge0f_(
	    void), echarge0g_(void);



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




/*     choose the method for summing over pairwise interactions */



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  warp.i  --  parameters for potential surface smoothing  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     m2           second moment of the GDA gaussian for each atom */
/*     deform       value of the smoothing deformation parameter */
/*     difft        diffusion coefficient for torsional potential */
/*     diffv        diffusion coefficient for van der Waals potential */
/*     diffc        diffusion coefficient for charge-charge potential */
/*     use_smooth   flag to use a potential energy smoothing method */
/*     use_dem      flag to use diffusion equation method potential */
/*     use_gda      flag to use gaussian density annealing potential */
/*     use_tophat   flag to use analytical tophat smoothed potential */
/*     use_stophat  flag to use shifted tophat smoothed potential */


    if (warp_1.use_smooth__) {
	echarge0g_();
    } else if (cutoff_1.use_ewald__) {
	if (cutoff_1.use_clist__) {
	    echarge0f_();
	} else if (cutoff_1.use_lights__) {
	    echarge0e_();
	} else {
	    echarge0d_();
	}
    } else if (cutoff_1.use_clist__) {
	echarge0c_();
    } else if (cutoff_1.use_lights__) {
	echarge0b_();
    } else {
	echarge0a_();
    }
    return 0;
} /* echarge_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine echarge0a  --  charge energy via double loop  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "echarge0a" calculates the charge-charge interaction energy */
/*     using a pairwise double loop */


/* Subroutine */ int echarge0a_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__, r2;
    static integer ii, kk, in, kn, ic, kc;
    static doublereal rb, fi, xi, yi, zi, xc, yc, zc, xr, yr, zr, rc, rc2, 
	    rc3, rc4, rc5, rc6, rc7, fik, xic, yic, zic, fgrp;
    static logical usei;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal shift, taper, trans, cscale[25000];
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *), switch_(char *, ftnlen), groups_(
	    logical *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *);
    static logical proceed;


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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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




/*     zero out the charge interaction energy */

    energi_1.ec = 0.;
    if (charge_1.nion == 0) {
	return 0;
    }

/*     set array needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cscale[i__ - 1] = 1.;
    }

/*     set conversion factor, cutoff and switching coefficients */

    f = chgpot_1.electric / chgpot_1.dielec;
    switch_("CHARGE", (ftnlen)6);

/*     calculate the charge interaction energy term */

    i__1 = charge_1.nion - 1;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = charge_1.iion[ii - 1];
	in = charge_1.jion[ii - 1];
	ic = charge_1.kion[ii - 1];
	xic = atoms_1.x[ic - 1];
	yic = atoms_1.y[ic - 1];
	zic = atoms_1.z__[ic - 1];
	xi = atoms_1.x[i__ - 1] - xic;
	yi = atoms_1.y[i__ - 1] - yic;
	zi = atoms_1.z__[i__ - 1] - zic;
	fi = f * charge_1.pchg[ii - 1];
	usei = usage_1.use[i__ - 1] || usage_1.use[ic - 1];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = chgpot_1.c2scale;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = chgpot_1.c3scale;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = chgpot_1.c4scale;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = chgpot_1.c5scale;
	}

/*     decide whether to compute the current interaction */

	i__2 = charge_1.nion;
	for (kk = ii + 1; kk <= i__2; ++kk) {
	    k = charge_1.iion[kk - 1];
	    kn = charge_1.jion[kk - 1];
	    kc = charge_1.kion[kk - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1] || usage_1.use[kc - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xc = xic - atoms_1.x[kc - 1];
		yc = yic - atoms_1.y[kc - 1];
		zc = zic - atoms_1.z__[kc - 1];
		image_(&xc, &yc, &zc);
		rc2 = xc * xc + yc * yc + zc * zc;
		if (rc2 <= shunt_1.off2) {
		    xr = xc + xi - atoms_1.x[k - 1] + atoms_1.x[kc - 1];
		    yr = yc + yi - atoms_1.y[k - 1] + atoms_1.y[kc - 1];
		    zr = zc + zi - atoms_1.z__[k - 1] + atoms_1.z__[kc - 1];
		    r2 = xr * xr + yr * yr + zr * zr;
		    r__ = sqrt(r2);
		    rb = r__ + chgpot_1.ebuffer;
		    fik = fi * charge_1.pchg[kk - 1] * cscale[kn - 1];
		    e = fik / rb;

/*     use shifted energy switching if near the cutoff distance */

		    shift = fik / ((shunt_1.off + shunt_1.cut) * .5);
		    e -= shift;
		    if (rc2 > shunt_1.cut2) {
			rc = sqrt(rc2);
			rc3 = rc2 * rc;
			rc4 = rc2 * rc2;
			rc5 = rc2 * rc3;
			rc6 = rc3 * rc3;
			rc7 = rc3 * rc4;
			taper = shunt_1.c5 * rc5 + shunt_1.c4 * rc4 + 
				shunt_1.c3 * rc3 + shunt_1.c2 * rc2 + 
				shunt_1.c1 * rc + shunt_1.c0;
			trans = fik * (shunt_1.f7 * rc7 + shunt_1.f6 * rc6 + 
				shunt_1.f5 * rc5 + shunt_1.f4 * rc4 + 
				shunt_1.f3 * rc3 + shunt_1.f2 * rc2 + 
				shunt_1.f1 * rc + shunt_1.f0);
			e = e * taper + trans;
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
		    }

/*     increment the overall charge-charge energy component */

		    energi_1.ec += e;
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = 1.;
	}
    }

/*     for periodic boundary conditions with large cutoffs */
/*     neighbors must be found by the replicates method */

    if (! bound_1.use_replica__) {
	return 0;
    }

/*     calculate interaction energy with other unit cells */

    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = charge_1.iion[ii - 1];
	in = charge_1.jion[ii - 1];
	ic = charge_1.kion[ii - 1];
	usei = usage_1.use[i__ - 1] || usage_1.use[ic - 1];
	xic = atoms_1.x[ic - 1];
	yic = atoms_1.y[ic - 1];
	zic = atoms_1.z__[ic - 1];
	xi = atoms_1.x[i__ - 1] - xic;
	yi = atoms_1.y[i__ - 1] - yic;
	zi = atoms_1.z__[i__ - 1] - zic;
	fi = f * charge_1.pchg[ii - 1];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = chgpot_1.c2scale;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = chgpot_1.c3scale;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = chgpot_1.c4scale;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = chgpot_1.c5scale;
	}

/*     decide whether to compute the current interaction */

	i__2 = charge_1.nion;
	for (kk = ii; kk <= i__2; ++kk) {
	    k = charge_1.iion[kk - 1];
	    kn = charge_1.jion[kk - 1];
	    kc = charge_1.kion[kk - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1] || usage_1.use[kc - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		i__3 = cell_1.ncell;
		for (j = 1; j <= i__3; ++j) {
		    xc = xic - atoms_1.x[kc - 1];
		    yc = yic - atoms_1.y[kc - 1];
		    zc = zic - atoms_1.z__[kc - 1];
		    imager_(&xc, &yc, &zc, &j);
		    rc2 = xc * xc + yc * yc + zc * zc;
		    if (rc2 <= shunt_1.off2) {
			xr = xc + xi - atoms_1.x[k - 1] + atoms_1.x[kc - 1];
			yr = yc + yi - atoms_1.y[k - 1] + atoms_1.y[kc - 1];
			zr = zc + zi - atoms_1.z__[k - 1] + atoms_1.z__[kc - 
				1];
			r2 = xr * xr + yr * yr + zr * zr;
			r__ = sqrt(r2);
			rb = r__ + chgpot_1.ebuffer;
			fik = fi * charge_1.pchg[kk - 1];
			if (bound_1.use_polymer__) {
			    if (r2 <= bound_1.polycut2) {
				fik *= cscale[kn - 1];
			    }
			}
			e = fik / rb;

/*     use shifted energy switching if near the cutoff distance */

			shift = fik / ((shunt_1.off + shunt_1.cut) * .5);
			e -= shift;
			if (rc2 > shunt_1.cut2) {
			    rc = sqrt(rc2);
			    rc3 = rc2 * rc;
			    rc4 = rc2 * rc2;
			    rc5 = rc2 * rc3;
			    rc6 = rc3 * rc3;
			    rc7 = rc3 * rc4;
			    taper = shunt_1.c5 * rc5 + shunt_1.c4 * rc4 + 
				    shunt_1.c3 * rc3 + shunt_1.c2 * rc2 + 
				    shunt_1.c1 * rc + shunt_1.c0;
			    trans = fik * (shunt_1.f7 * rc7 + shunt_1.f6 * 
				    rc6 + shunt_1.f5 * rc5 + shunt_1.f4 * rc4 
				    + shunt_1.f3 * rc3 + shunt_1.f2 * rc2 + 
				    shunt_1.f1 * rc + shunt_1.f0);
			    e = e * taper + trans;
			}

/*     scale the interaction based on its group membership */

			if (group_1.use_group__) {
			    e *= fgrp;
			}

/*     increment the overall charge-charge energy component */

			if (i__ == k) {
			    e *= .5;
			}
			energi_1.ec += e;
		    }
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = 1.;
	}
    }
    return 0;
} /* echarge0a_ */

#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine echarge0b  --  method of lights charge energy  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "echarge0b" calculates the charge-charge interaction energy */
/*     using the method of lights */


/* Subroutine */ int echarge0b_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__, r2;
    static integer ii, kk, in, ic, kn, kc;
    static doublereal rb, fi, xi, yi, zi, xc, yc, zc, xr, yr, zr, rc, rc2, 
	    rc3, rc4, rc5, rc6, rc7, fik, xic;
    static integer kgy, kgz;
    static doublereal yic, zic;
    static integer kmap;
    static doublereal fgrp;
    static logical usei;
    static integer stop;
    static doublereal shift, taper;
    static logical prime;
    static doublereal trans;
    static integer start;
    static doublereal xsort[200000], ysort[200000], zsort[200000], cscale[
	    25000];
    static logical repeat;
    extern /* Subroutine */ int lights_(doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *), switch_(char *, ftnlen), groups_(
	    logical *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *);
    static logical proceed;


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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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
/*     ##  light.i  --  indices for method of lights pair neighbors  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     nlight  total number of sites for method of lights calculation */
/*     kbx     low index of neighbors of each site in the x-sorted list */
/*     kby     low index of neighbors of each site in the y-sorted list */
/*     kbz     low index of neighbors of each site in the z-sorted list */
/*     kex     high index of neighbors of each site in the x-sorted list */
/*     key     high index of neighbors of each site in the y-sorted list */
/*     kez     high index of neighbors of each site in the z-sorted list */
/*     locx    pointer from x-sorted list into original interaction list */
/*     locy    pointer from y-sorted list into original interaction list */
/*     locz    pointer from z-sorted list into original interaction list */
/*     rgx     pointer from original interaction list into x-sorted list */
/*     rgy     pointer from original interaction list into y-sorted list */
/*     rgz     pointer from original interaction list into z-sorted list */




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




/*     zero out the charge interaction energy */

    energi_1.ec = 0.;
    if (charge_1.nion == 0) {
	return 0;
    }

/*     set array needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cscale[i__ - 1] = 1.;
    }

/*     set conversion factor, cutoff and switching coefficients */

    f = chgpot_1.electric / chgpot_1.dielec;
    switch_("CHARGE", (ftnlen)6);

/*     transfer the interaction site coordinates to sorting arrays */

    i__1 = charge_1.nion;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = charge_1.kion[i__ - 1];
	xsort[i__ - 1] = atoms_1.x[k - 1];
	ysort[i__ - 1] = atoms_1.y[k - 1];
	zsort[i__ - 1] = atoms_1.z__[k - 1];
    }

/*     use the method of lights to generate neighbors */

    lights_(&shunt_1.off, &charge_1.nion, xsort, ysort, zsort);

/*     loop over all atoms computing the interactions */

    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = charge_1.iion[ii - 1];
	in = charge_1.jion[ii - 1];
	ic = charge_1.kion[ii - 1];
	xic = xsort[light_1.rgx[ii - 1] - 1];
	yic = ysort[light_1.rgy[ii - 1] - 1];
	zic = zsort[light_1.rgz[ii - 1] - 1];
	xi = atoms_1.x[i__ - 1] - atoms_1.x[ic - 1];
	yi = atoms_1.y[i__ - 1] - atoms_1.y[ic - 1];
	zi = atoms_1.z__[i__ - 1] - atoms_1.z__[ic - 1];
	fi = f * charge_1.pchg[ii - 1];
	usei = usage_1.use[i__ - 1] || usage_1.use[ic - 1];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = chgpot_1.c2scale;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = chgpot_1.c3scale;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = chgpot_1.c4scale;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = chgpot_1.c5scale;
	}

/*     loop over method of lights neighbors of current atom */

	if (light_1.kbx[ii - 1] <= light_1.kex[ii - 1]) {
	    repeat = FALSE_;
	    start = light_1.kbx[ii - 1] + 1;
	    stop = light_1.kex[ii - 1];
	} else {
	    repeat = TRUE_;
	    start = 1;
	    stop = light_1.kex[ii - 1];
	}
L10:
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    kk = light_1.locx[j - 1];
	    kgy = light_1.rgy[kk - 1];
	    if (light_1.kby[ii - 1] <= light_1.key[ii - 1]) {
		if (kgy < light_1.kby[ii - 1] || kgy > light_1.key[ii - 1]) {
		    goto L20;
		}
	    } else {
		if (kgy < light_1.kby[ii - 1] && kgy > light_1.key[ii - 1]) {
		    goto L20;
		}
	    }
	    kgz = light_1.rgz[kk - 1];
	    if (light_1.kbz[ii - 1] <= light_1.kez[ii - 1]) {
		if (kgz < light_1.kbz[ii - 1] || kgz > light_1.kez[ii - 1]) {
		    goto L20;
		}
	    } else {
		if (kgz < light_1.kbz[ii - 1] && kgz > light_1.kez[ii - 1]) {
		    goto L20;
		}
	    }
	    kmap = kk - (kk - 1) / charge_1.nion * charge_1.nion;
	    k = charge_1.iion[kmap - 1];
	    kn = charge_1.jion[kmap - 1];
	    kc = charge_1.kion[kmap - 1];
	    prime = kk <= charge_1.nion;

/*     decide whether to compute the current interaction */

	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1] || usage_1.use[kc - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xc = xic - xsort[j - 1];
		yc = yic - ysort[kgy - 1];
		zc = zic - zsort[kgz - 1];
		if (bound_1.use_bounds__) {
		    if (abs(xc) > cell_1.xcell2) {
			xc -= d_sign(&cell_1.xcell, &xc);
		    }
		    if (abs(yc) > cell_1.ycell2) {
			yc -= d_sign(&cell_1.ycell, &yc);
		    }
		    if (abs(zc) > cell_1.zcell2) {
			zc -= d_sign(&cell_1.zcell, &zc);
		    }
		    if (boxes_1.monoclinic) {
			xc += zc * boxes_1.beta_cos__;
			zc *= boxes_1.beta_sin__;
		    } else if (boxes_1.triclinic) {
			xc = xc + yc * boxes_1.gamma_cos__ + zc * 
				boxes_1.beta_cos__;
			yc = yc * boxes_1.gamma_sin__ + zc * 
				boxes_1.beta_term__;
			zc *= boxes_1.gamma_term__;
		    }
		}
		rc2 = xc * xc + yc * yc + zc * zc;
		if (rc2 <= shunt_1.off2) {
		    xr = xc + xi - atoms_1.x[k - 1] + atoms_1.x[kc - 1];
		    yr = yc + yi - atoms_1.y[k - 1] + atoms_1.y[kc - 1];
		    zr = zc + zi - atoms_1.z__[k - 1] + atoms_1.z__[kc - 1];
		    r2 = xr * xr + yr * yr + zr * zr;
		    r__ = sqrt(r2);
		    rb = r__ + chgpot_1.ebuffer;
		    fik = fi * charge_1.pchg[kmap - 1];
		    if (prime) {
			fik *= cscale[kn - 1];
		    }
		    if (bound_1.use_polymer__) {
			if (r2 > bound_1.polycut2) {
			    fik = fi * charge_1.pchg[kmap - 1];
			}
		    }
		    e = fik / rb;

/*     use shifted energy switching if near the cutoff distance */

		    shift = fik / ((shunt_1.off + shunt_1.cut) * .5);
		    e -= shift;
		    if (rc2 > shunt_1.cut2) {
			rc = sqrt(rc2);
			rc3 = rc2 * rc;
			rc4 = rc2 * rc2;
			rc5 = rc2 * rc3;
			rc6 = rc3 * rc3;
			rc7 = rc3 * rc4;
			taper = shunt_1.c5 * rc5 + shunt_1.c4 * rc4 + 
				shunt_1.c3 * rc3 + shunt_1.c2 * rc2 + 
				shunt_1.c1 * rc + shunt_1.c0;
			trans = fik * (shunt_1.f7 * rc7 + shunt_1.f6 * rc6 + 
				shunt_1.f5 * rc5 + shunt_1.f4 * rc4 + 
				shunt_1.f3 * rc3 + shunt_1.f2 * rc2 + 
				shunt_1.f1 * rc + shunt_1.f0);
			e = e * taper + trans;
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
		    }

/*     increment the overall charge-charge energy component */

		    energi_1.ec += e;
		}
	    }
L20:
	    ;
	}
	if (repeat) {
	    repeat = FALSE_;
	    start = light_1.kbx[ii - 1] + 1;
	    stop = light_1.nlight;
	    goto L10;
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = 1.;
	}
    }
    return 0;
} /* echarge0b_ */

#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine echarge0c  --  charge energy via neighbor list  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "echarge0c" calculates the charge-charge interaction energy */
/*     using a pairwise neighbor list */


/* Subroutine */ int echarge0c_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__, r2;
    static integer ii, kk, in, kn, ic, kc;
    static doublereal rb, fi, xi, yi, zi, xc, yc, zc, xr, yr, zr, rc, rc2, 
	    rc3, rc4, rc5, rc6, rc7, fik;
    static integer kkk;
    static doublereal xic, yic, zic, fgrp;
    static logical usei;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal taper, shift, trans, cscale[25000];
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *);
    static logical proceed;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define elst_ref(a_1,a_2) neigh_1.elst[(a_2)*1200 + a_1 - 1201]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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




/*     zero out the charge interaction energy */

    energi_1.ec = 0.;
    if (charge_1.nion == 0) {
	return 0;
    }

/*     set array needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cscale[i__ - 1] = 1.;
    }

/*     set conversion factor, cutoff and switching coefficients */

    f = chgpot_1.electric / chgpot_1.dielec;
    switch_("CHARGE", (ftnlen)6);

/*     calculate the charge interaction energy term */

    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = charge_1.iion[ii - 1];
	in = charge_1.jion[ii - 1];
	ic = charge_1.kion[ii - 1];
	xic = atoms_1.x[ic - 1];
	yic = atoms_1.y[ic - 1];
	zic = atoms_1.z__[ic - 1];
	xi = atoms_1.x[i__ - 1] - xic;
	yi = atoms_1.y[i__ - 1] - yic;
	zi = atoms_1.z__[i__ - 1] - zic;
	fi = f * charge_1.pchg[ii - 1];
	usei = usage_1.use[i__ - 1] || usage_1.use[ic - 1];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = chgpot_1.c2scale;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = chgpot_1.c3scale;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = chgpot_1.c4scale;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = chgpot_1.c5scale;
	}

/*     decide whether to compute the current interaction */

	i__2 = neigh_1.nelst[ii - 1];
	for (kkk = 1; kkk <= i__2; ++kkk) {
	    kk = elst_ref(kkk, ii);
	    k = charge_1.iion[kk - 1];
	    kn = charge_1.jion[kk - 1];
	    kc = charge_1.kion[kk - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1] || usage_1.use[kc - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xc = xic - atoms_1.x[kc - 1];
		yc = yic - atoms_1.y[kc - 1];
		zc = zic - atoms_1.z__[kc - 1];
		image_(&xc, &yc, &zc);
		rc2 = xc * xc + yc * yc + zc * zc;
		if (rc2 <= shunt_1.off2) {
		    xr = xc + xi - atoms_1.x[k - 1] + atoms_1.x[kc - 1];
		    yr = yc + yi - atoms_1.y[k - 1] + atoms_1.y[kc - 1];
		    zr = zc + zi - atoms_1.z__[k - 1] + atoms_1.z__[kc - 1];
		    r2 = xr * xr + yr * yr + zr * zr;
		    r__ = sqrt(r2);
		    rb = r__ + chgpot_1.ebuffer;
		    fik = fi * charge_1.pchg[kk - 1] * cscale[kn - 1];
		    e = fik / rb;

/*     use shifted energy switching if near the cutoff distance */

		    shift = fik / ((shunt_1.off + shunt_1.cut) * .5);
		    e -= shift;
		    if (rc2 > shunt_1.cut2) {
			rc = sqrt(rc2);
			rc3 = rc2 * rc;
			rc4 = rc2 * rc2;
			rc5 = rc2 * rc3;
			rc6 = rc3 * rc3;
			rc7 = rc3 * rc4;
			taper = shunt_1.c5 * rc5 + shunt_1.c4 * rc4 + 
				shunt_1.c3 * rc3 + shunt_1.c2 * rc2 + 
				shunt_1.c1 * rc + shunt_1.c0;
			trans = fik * (shunt_1.f7 * rc7 + shunt_1.f6 * rc6 + 
				shunt_1.f5 * rc5 + shunt_1.f4 * rc4 + 
				shunt_1.f3 * rc3 + shunt_1.f2 * rc2 + 
				shunt_1.f1 * rc + shunt_1.f0);
			e = e * taper + trans;
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
		    }

/*     increment the overall charge-charge energy component */

		    energi_1.ec += e;
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = 1.;
	}
    }
    return 0;
} /* echarge0c_ */

#undef elst_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine echarge0d  --  double loop Ewald charge energy  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "echarge0d" calculates the charge-charge interaction energy */
/*     using a particle mesh Ewald summation */


/* Subroutine */ int echarge0d_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__, scaleterm, r2;
    static integer ii, kk, in, kn;
    static doublereal fs, fi, rb, xi, yi, zi, xd, yd, zd, xr, yr, zr, fik, 
	    rew;
    extern doublereal erfc_(doublereal *);
    static doublereal fgrp;
    static logical usei;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal scale, cscale[25000];
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *), switch_(char *, ftnlen), groups_(
	    logical *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *), ecrecip_(void);
    static logical proceed;
    static doublereal erfterm;


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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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




/*     zero out the Ewald charge interaction energy */

    energi_1.ec = 0.;
    if (charge_1.nion == 0) {
	return 0;
    }

/*     set array needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cscale[i__ - 1] = 1.;
    }

/*     set conversion factor, cutoff and switching coefficients */

    f = chgpot_1.electric / chgpot_1.dielec;
    switch_("EWALD", (ftnlen)5);

/*     compute the Ewald self-energy term over all the atoms */

    fs = -f * ewald_1.aewald / 1.772453850905516027;
    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
/* Computing 2nd power */
	d__1 = charge_1.pchg[ii - 1];
	e = fs * (d__1 * d__1);
	energi_1.ec += e;
    }

/*     compute the cell dipole boundary correction term */

    if (s_cmp(ewald_1.boundary, "VACUUM", (ftnlen)7, (ftnlen)6) == 0) {
	xd = 0.;
	yd = 0.;
	zd = 0.;
	i__1 = charge_1.nion;
	for (ii = 1; ii <= i__1; ++ii) {
	    i__ = charge_1.iion[ii - 1];
	    xd += charge_1.pchg[ii - 1] * atoms_1.x[i__ - 1];
	    yd += charge_1.pchg[ii - 1] * atoms_1.y[i__ - 1];
	    zd += charge_1.pchg[ii - 1] * atoms_1.z__[i__ - 1];
	}
	e = f * .66666666666666663 * (3.141592653589793238 / boxes_1.volbox) *
		 (xd * xd + yd * yd + zd * zd);
	energi_1.ec += e;
    }

/*     compute the reciprocal space part of the Ewald summation */

    ecrecip_();

/*     compute the real space portion of the Ewald summation */

    i__1 = charge_1.nion - 1;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = charge_1.iion[ii - 1];
	in = charge_1.jion[ii - 1];
	usei = usage_1.use[i__ - 1];
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	fi = f * charge_1.pchg[ii - 1];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = chgpot_1.c2scale;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = chgpot_1.c3scale;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = chgpot_1.c4scale;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = chgpot_1.c5scale;
	}

/*     decide whether to compute the current interaction */

	i__2 = charge_1.nion;
	for (kk = ii + 1; kk <= i__2; ++kk) {
	    k = charge_1.iion[kk - 1];
	    kn = charge_1.jion[kk - 1];
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    proceed = TRUE_;
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xr = xi - atoms_1.x[k - 1];
		yr = yi - atoms_1.y[k - 1];
		zr = zi - atoms_1.z__[k - 1];
		image_(&xr, &yr, &zr);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    r__ = sqrt(r2);
		    rb = r__ + chgpot_1.ebuffer;
		    fik = fi * charge_1.pchg[kk - 1];
		    rew = ewald_1.aewald * r__;
		    erfterm = erfc_(&rew);
		    scale = cscale[kn - 1];
		    if (group_1.use_group__) {
			scale *= fgrp;
		    }
		    scaleterm = scale - 1.;
		    e = fik / rb * (erfterm + scaleterm);

/*     increment the overall charge-charge energy component */

		    energi_1.ec += e;
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = 1.;
	}
    }

/*     for periodic boundary conditions with large cutoffs */
/*     neighbors must be found by the replicates method */

    if (! bound_1.use_replica__) {
	return 0;
    }

/*     calculate real space portion involving other unit cells */

    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = charge_1.iion[ii - 1];
	in = charge_1.jion[ii - 1];
	usei = usage_1.use[i__ - 1];
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	fi = f * charge_1.pchg[ii - 1];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = chgpot_1.c2scale;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = chgpot_1.c3scale;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = chgpot_1.c4scale;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = chgpot_1.c5scale;
	}

/*     decide whether to compute the current interaction */

	i__2 = charge_1.nion;
	for (kk = ii; kk <= i__2; ++kk) {
	    k = charge_1.iion[kk - 1];
	    kn = charge_1.jion[kk - 1];
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    proceed = TRUE_;
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		i__3 = cell_1.ncell;
		for (j = 1; j <= i__3; ++j) {
		    xr = xi - atoms_1.x[k - 1];
		    yr = yi - atoms_1.y[k - 1];
		    zr = zi - atoms_1.z__[k - 1];
		    imager_(&xr, &yr, &zr, &j);
		    r2 = xr * xr + yr * yr + zr * zr;
		    if (r2 <= shunt_1.off2) {
			r__ = sqrt(r2);
			rb = r__ + chgpot_1.ebuffer;
			fik = fi * charge_1.pchg[kk - 1];
			rew = ewald_1.aewald * r__;
			erfterm = erfc_(&rew);
			scale = 1.;
			if (group_1.use_group__) {
			    scale *= fgrp;
			}
			if (bound_1.use_polymer__) {
			    if (r2 <= bound_1.polycut2) {
				scale *= cscale[kn - 1];
			    }
			}
			scaleterm = scale - 1.;
			e = fik / rb * (erfterm + scaleterm);

/*     increment the overall charge-charge energy component; */
/*     interaction of an atom with its own image counts half */

			if (i__ == k) {
			    e *= .5;
			}
			energi_1.ec += e;
		    }
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = 1.;
	}
    }
    return 0;
} /* echarge0d_ */

#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine echarge0e  --  Ewald charge energy via lights  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "echarge0e" calculates the charge-charge interaction energy */
/*     using a particle mesh Ewald summation and the method of lights */


/* Subroutine */ int echarge0e_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__, scaleterm, r2;
    static integer ii, kk, in, kn;
    static doublereal fs, fi, rb, xi, yi, zi, xd, yd, zd, xr, yr, zr, fik;
    static integer kgy, kgz;
    static doublereal rew;
    extern doublereal erfc_(doublereal *);
    static integer kmap;
    static doublereal fgrp;
    static logical usei;
    static integer stop;
    static doublereal scale;
    static logical prime;
    static integer start;
    static doublereal xsort[200000], ysort[200000], zsort[200000], cscale[
	    25000];
    static logical repeat;
    extern /* Subroutine */ int lights_(doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *), switch_(char *, ftnlen), groups_(
	    logical *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *), ecrecip_(void);
    static logical proceed;
    static doublereal erfterm;


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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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
/*     ##  light.i  --  indices for method of lights pair neighbors  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     nlight  total number of sites for method of lights calculation */
/*     kbx     low index of neighbors of each site in the x-sorted list */
/*     kby     low index of neighbors of each site in the y-sorted list */
/*     kbz     low index of neighbors of each site in the z-sorted list */
/*     kex     high index of neighbors of each site in the x-sorted list */
/*     key     high index of neighbors of each site in the y-sorted list */
/*     kez     high index of neighbors of each site in the z-sorted list */
/*     locx    pointer from x-sorted list into original interaction list */
/*     locy    pointer from y-sorted list into original interaction list */
/*     locz    pointer from z-sorted list into original interaction list */
/*     rgx     pointer from original interaction list into x-sorted list */
/*     rgy     pointer from original interaction list into y-sorted list */
/*     rgz     pointer from original interaction list into z-sorted list */




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




/*     zero out the Ewald charge interaction energy */

    energi_1.ec = 0.;
    if (charge_1.nion == 0) {
	return 0;
    }

/*     set array needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cscale[i__ - 1] = 1.;
    }

/*     set conversion factor, cutoff and switching coefficients */

    f = chgpot_1.electric / chgpot_1.dielec;
    switch_("EWALD", (ftnlen)5);

/*     compute the Ewald self-energy term over all the atoms */

    fs = -f * ewald_1.aewald / 1.772453850905516027;
    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
/* Computing 2nd power */
	d__1 = charge_1.pchg[ii - 1];
	e = fs * (d__1 * d__1);
	energi_1.ec += e;
    }

/*     compute the cell dipole boundary correction term */

    if (s_cmp(ewald_1.boundary, "VACUUM", (ftnlen)7, (ftnlen)6) == 0) {
	xd = 0.;
	yd = 0.;
	zd = 0.;
	i__1 = charge_1.nion;
	for (ii = 1; ii <= i__1; ++ii) {
	    i__ = charge_1.iion[ii - 1];
	    xd += charge_1.pchg[ii - 1] * atoms_1.x[i__ - 1];
	    yd += charge_1.pchg[ii - 1] * atoms_1.y[i__ - 1];
	    zd += charge_1.pchg[ii - 1] * atoms_1.z__[i__ - 1];
	}
	e = f * .66666666666666663 * (3.141592653589793238 / boxes_1.volbox) *
		 (xd * xd + yd * yd + zd * zd);
	energi_1.ec += e;
    }

/*     compute the reciprocal space part of the Ewald summation */

    ecrecip_();

/*     compute the real space portion of the Ewald summation; */
/*     transfer the interaction site coordinates to sorting arrays */

    i__1 = charge_1.nion;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = charge_1.iion[i__ - 1];
	xsort[i__ - 1] = atoms_1.x[k - 1];
	ysort[i__ - 1] = atoms_1.y[k - 1];
	zsort[i__ - 1] = atoms_1.z__[k - 1];
    }

/*     use the method of lights to generate neighbors */

    lights_(&shunt_1.off, &charge_1.nion, xsort, ysort, zsort);

/*     loop over all atoms computing the interactions */

    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = charge_1.iion[ii - 1];
	in = charge_1.jion[ii - 1];
	xi = xsort[light_1.rgx[ii - 1] - 1];
	yi = ysort[light_1.rgy[ii - 1] - 1];
	zi = zsort[light_1.rgz[ii - 1] - 1];
	fi = f * charge_1.pchg[ii - 1];
	usei = usage_1.use[i__ - 1];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = chgpot_1.c2scale;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = chgpot_1.c3scale;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = chgpot_1.c4scale;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = chgpot_1.c5scale;
	}

/*     loop over method of lights neighbors of current atom */

	if (light_1.kbx[ii - 1] <= light_1.kex[ii - 1]) {
	    repeat = FALSE_;
	    start = light_1.kbx[ii - 1] + 1;
	    stop = light_1.kex[ii - 1];
	} else {
	    repeat = TRUE_;
	    start = 1;
	    stop = light_1.kex[ii - 1];
	}
L10:
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    kk = light_1.locx[j - 1];
	    kgy = light_1.rgy[kk - 1];
	    if (light_1.kby[ii - 1] <= light_1.key[ii - 1]) {
		if (kgy < light_1.kby[ii - 1] || kgy > light_1.key[ii - 1]) {
		    goto L20;
		}
	    } else {
		if (kgy < light_1.kby[ii - 1] && kgy > light_1.key[ii - 1]) {
		    goto L20;
		}
	    }
	    kgz = light_1.rgz[kk - 1];
	    if (light_1.kbz[ii - 1] <= light_1.kez[ii - 1]) {
		if (kgz < light_1.kbz[ii - 1] || kgz > light_1.kez[ii - 1]) {
		    goto L20;
		}
	    } else {
		if (kgz < light_1.kbz[ii - 1] && kgz > light_1.kez[ii - 1]) {
		    goto L20;
		}
	    }
	    kmap = kk - (kk - 1) / charge_1.nion * charge_1.nion;
	    k = charge_1.iion[kmap - 1];
	    kn = charge_1.jion[kmap - 1];
	    prime = kk <= charge_1.nion;

/*     decide whether to compute the current interaction */

	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    proceed = TRUE_;
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xr = xi - xsort[j - 1];
		yr = yi - ysort[kgy - 1];
		zr = zi - zsort[kgz - 1];
		if (bound_1.use_bounds__) {
		    if (abs(xr) > cell_1.xcell2) {
			xr -= d_sign(&cell_1.xcell, &xr);
		    }
		    if (abs(yr) > cell_1.ycell2) {
			yr -= d_sign(&cell_1.ycell, &yr);
		    }
		    if (abs(zr) > cell_1.zcell2) {
			zr -= d_sign(&cell_1.zcell, &zr);
		    }
		    if (boxes_1.monoclinic) {
			xr += zr * boxes_1.beta_cos__;
			zr *= boxes_1.beta_sin__;
		    } else if (boxes_1.triclinic) {
			xr = xr + yr * boxes_1.gamma_cos__ + zr * 
				boxes_1.beta_cos__;
			yr = yr * boxes_1.gamma_sin__ + zr * 
				boxes_1.beta_term__;
			zr *= boxes_1.gamma_term__;
		    }
		}
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    r__ = sqrt(r2);
		    rb = r__ + chgpot_1.ebuffer;
		    rew = ewald_1.aewald * r__;
		    erfterm = erfc_(&rew);
		    scale = 1.;
		    if (prime) {
			scale = cscale[kn - 1];
		    }
		    if (group_1.use_group__) {
			scale *= fgrp;
		    }
		    fik = fi * charge_1.pchg[kmap - 1];
		    if (bound_1.use_polymer__) {
			if (r2 > bound_1.polycut2) {
			    fik = fi * charge_1.pchg[kmap - 1];
			}
		    }
		    scaleterm = scale - 1.;
		    e = fik / rb * (erfterm + scaleterm);

/*     increment the overall charge-charge energy component */

		    energi_1.ec += e;
		}
	    }
L20:
	    ;
	}
	if (repeat) {
	    repeat = FALSE_;
	    start = light_1.kbx[ii - 1] + 1;
	    stop = light_1.nlight;
	    goto L10;
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = 1.;
	}
    }
    return 0;
} /* echarge0e_ */

#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine echarge0f  --  Ewald charge energy via list  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "echarge0f" calculates the charge-charge interaction energy */
/*     using a particle mesh Ewald summation and a neighbor list */


/* Subroutine */ int echarge0f_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__, scaleterm, r2;
    static integer ii, kk, in, kn;
    static doublereal fs, fi, rb, xi, yi, zi, xd, yd, zd, xr, yr, zr, fik;
    static integer kkk;
    static doublereal rew;
    extern doublereal erfc_(doublereal *);
    static doublereal fgrp;
    static logical usei;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal scale, cscale[25000];
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *), ecrecip_(void);
    static logical proceed;
    static doublereal erfterm;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define elst_ref(a_1,a_2) neigh_1.elst[(a_2)*1200 + a_1 - 1201]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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




/*     zero out the Ewald charge interaction energy */

    energi_1.ec = 0.;
    if (charge_1.nion == 0) {
	return 0;
    }

/*     set array needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cscale[i__ - 1] = 1.;
    }

/*     set conversion factor, cutoff and switching coefficients */

    f = chgpot_1.electric / chgpot_1.dielec;
    switch_("EWALD", (ftnlen)5);

/*     compute the Ewald self-energy term over all the atoms */

    fs = -f * ewald_1.aewald / 1.772453850905516027;
    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
/* Computing 2nd power */
	d__1 = charge_1.pchg[ii - 1];
	e = fs * (d__1 * d__1);
	energi_1.ec += e;
    }

/*     compute the cell dipole boundary correction term */

    if (s_cmp(ewald_1.boundary, "VACUUM", (ftnlen)7, (ftnlen)6) == 0) {
	xd = 0.;
	yd = 0.;
	zd = 0.;
	i__1 = charge_1.nion;
	for (ii = 1; ii <= i__1; ++ii) {
	    i__ = charge_1.iion[ii - 1];
	    xd += charge_1.pchg[ii - 1] * atoms_1.x[i__ - 1];
	    yd += charge_1.pchg[ii - 1] * atoms_1.y[i__ - 1];
	    zd += charge_1.pchg[ii - 1] * atoms_1.z__[i__ - 1];
	}
	e = f * .66666666666666663 * (3.141592653589793238 / boxes_1.volbox) *
		 (xd * xd + yd * yd + zd * zd);
	energi_1.ec += e;
    }

/*     compute the reciprocal space part of the Ewald summation */

    ecrecip_();

/*     compute the real space portion of the Ewald summation */

    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = charge_1.iion[ii - 1];
	in = charge_1.jion[ii - 1];
	usei = usage_1.use[i__ - 1];
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	fi = f * charge_1.pchg[ii - 1];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = chgpot_1.c2scale;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = chgpot_1.c3scale;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = chgpot_1.c4scale;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = chgpot_1.c5scale;
	}

/*     decide whether to compute the current interaction */

	i__2 = neigh_1.nelst[ii - 1];
	for (kkk = 1; kkk <= i__2; ++kkk) {
	    kk = elst_ref(kkk, ii);
	    k = charge_1.iion[kk - 1];
	    kn = charge_1.jion[kk - 1];
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    proceed = TRUE_;
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xr = xi - atoms_1.x[k - 1];
		yr = yi - atoms_1.y[k - 1];
		zr = zi - atoms_1.z__[k - 1];
		image_(&xr, &yr, &zr);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    r__ = sqrt(r2);
		    rb = r__ + chgpot_1.ebuffer;
		    fik = fi * charge_1.pchg[kk - 1];
		    rew = ewald_1.aewald * r__;
		    erfterm = erfc_(&rew);
		    scale = cscale[kn - 1];
		    if (group_1.use_group__) {
			scale *= fgrp;
		    }
		    scaleterm = scale - 1.;
		    e = fik / rb * (erfterm + scaleterm);

/*     increment the overall charge-charge energy component */

		    energi_1.ec += e;
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = 1.;
	}
    }
    return 0;
} /* echarge0f_ */

#undef elst_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine echarge0g  --  charge energy for smoothing  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "echarge0g" calculates the charge-charge interaction energy */
/*     for use with potential smoothing methods */


/* Subroutine */ int echarge0g_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__, r2;
    static integer ii, kk, in, kn;
    static doublereal rb, fi, xi, yi, zi, xr, yr, zr, rb2, fik;
    extern doublereal erf_(doublereal *);
    static doublereal fgrp;
    static logical usei;
    static doublereal width, wterm, width2, width3, cscale[25000];
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;


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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  warp.i  --  parameters for potential surface smoothing  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     m2           second moment of the GDA gaussian for each atom */
/*     deform       value of the smoothing deformation parameter */
/*     difft        diffusion coefficient for torsional potential */
/*     diffv        diffusion coefficient for van der Waals potential */
/*     diffc        diffusion coefficient for charge-charge potential */
/*     use_smooth   flag to use a potential energy smoothing method */
/*     use_dem      flag to use diffusion equation method potential */
/*     use_gda      flag to use gaussian density annealing potential */
/*     use_tophat   flag to use analytical tophat smoothed potential */
/*     use_stophat  flag to use shifted tophat smoothed potential */




/*     zero out the charge interaction energy */

    energi_1.ec = 0.;
    if (charge_1.nion == 0) {
	return 0;
    }

/*     set array needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cscale[i__ - 1] = 1.;
    }

/*     set the energy units conversion factor */

    f = chgpot_1.electric / chgpot_1.dielec;

/*     set the extent of smoothing to be performed */

    width = warp_1.deform * warp_1.diffc;
    if (warp_1.use_dem__) {
	if (width > 0.) {
	    width = .5 / sqrt(width);
	}
    } else if (warp_1.use_gda__) {
	wterm = sqrt(3. / (warp_1.diffc * 2.));
    }
    width2 = width * width;
    width3 = width * width2;

/*     calculate the charge interaction energy term */

    i__1 = charge_1.nion - 1;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = charge_1.iion[ii - 1];
	in = charge_1.jion[ii - 1];
	usei = usage_1.use[i__ - 1];
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	fi = f * charge_1.pchg[ii - 1];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = chgpot_1.c2scale;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = chgpot_1.c3scale;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = chgpot_1.c4scale;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = chgpot_1.c5scale;
	}

/*     decide whether to compute the current interaction */

	i__2 = charge_1.nion;
	for (kk = ii + 1; kk <= i__2; ++kk) {
	    k = charge_1.iion[kk - 1];
	    kn = charge_1.jion[kk - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xr = xi - atoms_1.x[k - 1];
		yr = yi - atoms_1.y[k - 1];
		zr = zi - atoms_1.z__[k - 1];
		r2 = xr * xr + yr * yr + zr * zr;
		r__ = sqrt(r2);
		rb = r__ + chgpot_1.ebuffer;
		fik = fi * charge_1.pchg[kk - 1] * cscale[kn - 1];
		e = fik / rb;

/*     transform the potential function via smoothing */

		if (warp_1.use_dem__) {
		    if (width > 0.) {
			d__1 = width * rb;
			e *= erf_(&d__1);
		    }
		} else if (warp_1.use_gda__) {
		    width = warp_1.m2[i__ - 1] + warp_1.m2[k - 1];
		    if (width > 0.) {
			width = wterm / sqrt(width);
			d__1 = width * rb;
			e *= erf_(&d__1);
		    }
		} else if (warp_1.use_tophat__) {
		    if (width > rb) {
			rb2 = rb * rb;
			e = fik * (width2 * 3. - rb2) / (width3 * 2.);
		    }
		} else if (warp_1.use_stophat__) {
		    e = fik / (rb + width);
		}

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    e *= fgrp;
		}

/*     increment the overall charge-charge energy component */

		energi_1.ec += e;
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i12_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n13[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i13_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n14[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i14_ref(j, in) - 1] = 1.;
	}
	i__2 = couple_1.n15[in - 1];
	for (j = 1; j <= i__2; ++j) {
	    cscale[i15_ref(j, in) - 1] = 1.;
	}
    }
    return 0;
} /* echarge0g_ */

#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine ecrecip  --  PME reciprocal space charge energy  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "ecrecip" evaluates the reciprocal space portion of the particle */
/*     mesh Ewald energy due to partial charges */

/*     literature reference: */

/*     U. Essmann, L. Perera, M. L Berkowitz, T. Darden, H. Lee and */
/*     L. G. Pedersen, "A Smooth Particle Mesh Ewald Method", Journal */
/*     of Chemical Physics, 103, 8577-8593 (1995) */

/*     modifications for nonperiodic systems suggested by Tom Darden */
/*     during May 2007 */


/* Subroutine */ int ecrecip_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_nint(doublereal *);
    integer i_sign(integer *, integer *);
    double exp(doublereal), sqrt(doublereal), cos(doublereal);

    /* Local variables */
    extern /* Subroutine */ int fftfront_(void);
    static doublereal e, f;
    static integer i__, j, k, m;
    static doublereal w;
    static integer i0, j0, k0, k1, k2, k3, m1, m2, m3;
    static doublereal h1, h2, t1, t2, t3, h3, r1, r2, r3;
    static integer ii;
    static doublereal fr, xi, yi, zi;
    static integer nf1, nf2, nf3, it1, it2, it3, nff, ifr;
    static doublereal hsq, term, denom, pterm, theta1[250000]	/* was [10][
	    25000] */, theta2[250000]	/* was [10][25000] */, theta3[250000]	
	    /* was [10][25000] */, struc2;
    static integer npoint;
    extern /* Subroutine */ int bspline_(doublereal *, integer *, doublereal *
	    );
    static doublereal expterm, volterm;


#define igrid_ref(a_1,a_2) pme_1.igrid[(a_2)*3 + a_1 - 4]
#define recip_ref(a_1,a_2) boxes_1.recip[(a_2)*3 + a_1 - 4]
#define qgrid_ref(a_1,a_2,a_3,a_4) pme_1.qgrid[(((a_4)*100 + (a_3))*100 + \
(a_2))*2 + a_1 - 20203]
#define theta1_ref(a_1,a_2) theta1[(a_2)*10 + a_1 - 11]
#define theta2_ref(a_1,a_2) theta2[(a_2)*10 + a_1 - 11]
#define theta3_ref(a_1,a_2) theta3[(a_2)*10 + a_1 - 11]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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

/*     get the B-spline coefficients for each charge site */

    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = charge_1.iion[ii - 1];
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	w = xi * recip_ref(1, 1) + yi * recip_ref(2, 1) + zi * recip_ref(3, 1)
		;
	fr = (doublereal) pme_1.nfft1 * (w - d_nint(&w) + .5);
	ifr = (integer) fr;
	w = fr - (doublereal) ifr;
	igrid_ref(1, i__) = ifr - pme_1.bsorder;
	bspline_(&w, &pme_1.bsorder, &theta1_ref(1, i__));
	w = xi * recip_ref(1, 2) + yi * recip_ref(2, 2) + zi * recip_ref(3, 2)
		;
	fr = (doublereal) pme_1.nfft2 * (w - d_nint(&w) + .5);
	ifr = (integer) fr;
	w = fr - (doublereal) ifr;
	igrid_ref(2, i__) = ifr - pme_1.bsorder;
	bspline_(&w, &pme_1.bsorder, &theta2_ref(1, i__));
	w = xi * recip_ref(1, 3) + yi * recip_ref(2, 3) + zi * recip_ref(3, 3)
		;
	fr = (doublereal) pme_1.nfft3 * (w - d_nint(&w) + .5);
	ifr = (integer) fr;
	w = fr - (doublereal) ifr;
	igrid_ref(3, i__) = ifr - pme_1.bsorder;
	bspline_(&w, &pme_1.bsorder, &theta3_ref(1, i__));
    }

/*     put the atomic partial charges onto the grid */

    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
	m = charge_1.iion[ii - 1];
	k0 = igrid_ref(3, m);
	i__2 = pme_1.bsorder;
	for (it3 = 1; it3 <= i__2; ++it3) {
	    ++k0;
	    k = k0 + 1 + (pme_1.nfft3 - i_sign(&pme_1.nfft3, &k0)) / 2;
	    t3 = theta3_ref(it3, m) * charge_1.pchg[ii - 1];
	    j0 = igrid_ref(2, m);
	    i__3 = pme_1.bsorder;
	    for (it2 = 1; it2 <= i__3; ++it2) {
		++j0;
		j = j0 + 1 + (pme_1.nfft2 - i_sign(&pme_1.nfft2, &j0)) / 2;
		t2 = theta2_ref(it2, m);
		term = t2 * t3;
		i0 = igrid_ref(1, m);
		i__4 = pme_1.bsorder;
		for (it1 = 1; it1 <= i__4; ++it1) {
		    ++i0;
		    i__ = i0 + 1 + (pme_1.nfft1 - i_sign(&pme_1.nfft1, &i0)) /
			     2;
		    t1 = theta1_ref(it1, m);
		    qgrid_ref(1, i__, j, k) = qgrid_ref(1, i__, j, k) + term *
			     t1;
		}
	    }
	}
    }

/*     perform the 3-D FFT forward transformation */

    fftfront_();

/*     use scalar sum to get the reciprocal space energy */

    f = chgpot_1.electric * .5 / chgpot_1.dielec;
    npoint = pme_1.nfft1 * pme_1.nfft2 * pme_1.nfft3;
/* Computing 2nd power */
    d__1 = 3.141592653589793238 / ewald_1.aewald;
    pterm = d__1 * d__1;
    volterm = boxes_1.volbox * 3.141592653589793238;
    nff = pme_1.nfft1 * pme_1.nfft2;
    nf1 = (pme_1.nfft1 + 1) / 2;
    nf2 = (pme_1.nfft2 + 1) / 2;
    nf3 = (pme_1.nfft3 + 1) / 2;
    i__1 = npoint - 1;
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
/* Computing 2nd power */
	    d__1 = qgrid_ref(1, k1, k2, k3);
/* Computing 2nd power */
	    d__2 = qgrid_ref(2, k1, k2, k3);
	    struc2 = d__1 * d__1 + d__2 * d__2;
	    e = f * expterm * struc2;
	    energi_1.ec += e;
	}
    }

/*     account for the zeroth grid point for a finite system */

    if (! bound_1.use_bounds__) {
	expterm = 1.5707963267948966 / boxes_1.xbox;
/* Computing 2nd power */
	d__1 = qgrid_ref(1, 1, 1, 1);
/* Computing 2nd power */
	d__2 = qgrid_ref(2, 1, 1, 1);
	struc2 = d__1 * d__1 + d__2 * d__2;
	e = f * expterm * struc2;
	energi_1.ec += e;
    }
    return 0;
} /* ecrecip_ */

#undef theta3_ref
#undef theta2_ref
#undef theta1_ref
#undef qgrid_ref
#undef recip_ref
#undef igrid_ref


