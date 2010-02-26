/* echarge2.f -- translated by f2c (version 20050501).
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
    doublereal off, off2, cut, cut2, c0, c1, c2, c3, c4, c5, f0, f1, f2, f3, 
	    f4, f5, f6, f7;
} shunt_;

#define shunt_1 shunt_

struct {
    doublereal aewald;
    char boundary[7];
} ewald_;

#define ewald_1 ewald_

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine echarge2  --  atomwise charge-charge Hessian  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "echarge2" calculates second derivatives of the */
/*     charge-charge interaction energy for a single atom */


/* Subroutine */ int echarge2_(integer *i__)
{
    extern /* Subroutine */ int echarge2a_(integer *), echarge2b_(integer *), 
	    echarge2c_(integer *);



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




/*     choose the method for summing over pairwise interactions */

    if (warp_1.use_smooth__) {
	echarge2c_(i__);
    } else if (cutoff_1.use_ewald__) {
	echarge2b_(i__);
    } else {
	echarge2a_(i__);
    }
    return 0;
} /* echarge2_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine echarge2a  --  charge Hessian via double loop  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "echarge2a" calculates second derivatives of the charge-charge */
/*     interaction energy for a single atom using a pairwise double loop */


/* Subroutine */ int echarge2a_(integer *i__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e;
    static integer j, k;
    static doublereal r__, r2, r3, r4, r5, r6, r7, de, fi, rb;
    static integer kk, in, kn;
    static doublereal xi, yi, zi, xr, yr, zr, d2e, rb2, fik, fgrp, term[9]	
	    /* was [3][3] */, d2edx, d2edy, d2edz;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static integer jcell;
    static doublereal shift, taper, trans, cscale[25000];
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal dtaper, dtrans;
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *);
    static doublereal d2taper, d2trans;
    static logical proceed;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define term_ref(a_1,a_2) term[(a_2)*3 + a_1 - 4]
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




/*     first see if the atom of interest carries a charge */

    i__1 = charge_1.nion;
    for (k = 1; k <= i__1; ++k) {
	if (charge_1.iion[k - 1] == *i__) {
	    fi = chgpot_1.electric * charge_1.pchg[k - 1] / chgpot_1.dielec;
	    in = charge_1.jion[k - 1];
	    goto L10;
	}
    }
    return 0;
L10:

/*     store the coordinates of the atom of interest */

    xi = atoms_1.x[*i__ - 1];
    yi = atoms_1.y[*i__ - 1];
    zi = atoms_1.z__[*i__ - 1];

/*     set array needed to scale connected atom interactions */

    i__1 = charge_1.nion;
    for (j = 1; j <= i__1; ++j) {
	cscale[charge_1.iion[j - 1] - 1] = 1.;
    }
    i__1 = couple_1.n12[in - 1];
    for (j = 1; j <= i__1; ++j) {
	cscale[i12_ref(j, in) - 1] = chgpot_1.c2scale;
    }
    i__1 = couple_1.n13[in - 1];
    for (j = 1; j <= i__1; ++j) {
	cscale[i13_ref(j, in) - 1] = chgpot_1.c3scale;
    }
    i__1 = couple_1.n14[in - 1];
    for (j = 1; j <= i__1; ++j) {
	cscale[i14_ref(j, in) - 1] = chgpot_1.c4scale;
    }
    i__1 = couple_1.n15[in - 1];
    for (j = 1; j <= i__1; ++j) {
	cscale[i15_ref(j, in) - 1] = chgpot_1.c5scale;
    }

/*     set cutoff distances and switching function coefficients */

    switch_("CHARGE", (ftnlen)6);

/*     calculate the charge interaction energy Hessian elements */

    i__1 = charge_1.nion;
    for (kk = 1; kk <= i__1; ++kk) {
	k = charge_1.iion[kk - 1];
	kn = charge_1.jion[kk - 1];
	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, i__, &k, &c__0, &c__0, &c__0, &c__0);
	}
	if (proceed) {
	    proceed = kn != *i__;
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
		rb2 = rb * rb;
		fik = fi * charge_1.pchg[kk - 1] * cscale[kn - 1];

/*     compute chain rule terms for Hessian matrix elements */

		de = -fik / rb2;
		d2e = de * -2. / rb;

/*     use shifted energy switching if near the cutoff distance */

		if (r2 > shunt_1.cut2) {
		    e = fik / r__;
		    shift = fik / ((shunt_1.off + shunt_1.cut) * .5);
		    e -= shift;
		    r3 = r2 * r__;
		    r4 = r2 * r2;
		    r5 = r2 * r3;
		    r6 = r3 * r3;
		    r7 = r3 * r4;
		    taper = shunt_1.c5 * r5 + shunt_1.c4 * r4 + shunt_1.c3 * 
			    r3 + shunt_1.c2 * r2 + shunt_1.c1 * r__ + 
			    shunt_1.c0;
		    dtaper = shunt_1.c5 * 5. * r4 + shunt_1.c4 * 4. * r3 + 
			    shunt_1.c3 * 3. * r2 + shunt_1.c2 * 2. * r__ + 
			    shunt_1.c1;
		    d2taper = shunt_1.c5 * 20. * r3 + shunt_1.c4 * 12. * r2 + 
			    shunt_1.c3 * 6. * r__ + shunt_1.c2 * 2.;
		    trans = fik * (shunt_1.f7 * r7 + shunt_1.f6 * r6 + 
			    shunt_1.f5 * r5 + shunt_1.f4 * r4 + shunt_1.f3 * 
			    r3 + shunt_1.f2 * r2 + shunt_1.f1 * r__ + 
			    shunt_1.f0);
		    dtrans = fik * (shunt_1.f7 * 7. * r6 + shunt_1.f6 * 6. * 
			    r5 + shunt_1.f5 * 5. * r4 + shunt_1.f4 * 4. * r3 
			    + shunt_1.f3 * 3. * r2 + shunt_1.f2 * 2. * r__ + 
			    shunt_1.f1);
		    d2trans = fik * (shunt_1.f7 * 42. * r5 + shunt_1.f6 * 30. 
			    * r4 + shunt_1.f5 * 20. * r3 + shunt_1.f4 * 12. * 
			    r2 + shunt_1.f3 * 6. * r__ + shunt_1.f2 * 2.);
		    d2e = e * d2taper + de * 2. * dtaper + d2e * taper + 
			    d2trans;
		    de = e * dtaper + de * taper + dtrans;
		}

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    de *= fgrp;
		    d2e *= fgrp;
		}

/*     form the individual Hessian element components */

		de /= r__;
		d2e = (d2e - de) / r2;
		d2edx = d2e * xr;
		d2edy = d2e * yr;
		d2edz = d2e * zr;
		term_ref(1, 1) = d2edx * xr + de;
		term_ref(1, 2) = d2edx * yr;
		term_ref(1, 3) = d2edx * zr;
		term_ref(2, 1) = term_ref(1, 2);
		term_ref(2, 2) = d2edy * yr + de;
		term_ref(2, 3) = d2edy * zr;
		term_ref(3, 1) = term_ref(1, 3);
		term_ref(3, 2) = term_ref(2, 3);
		term_ref(3, 3) = d2edz * zr + de;

/*     increment diagonal and non-diagonal Hessian elements */

		for (j = 1; j <= 3; ++j) {
		    hessx_ref(j, *i__) = hessx_ref(j, *i__) + term_ref(1, j);
		    hessy_ref(j, *i__) = hessy_ref(j, *i__) + term_ref(2, j);
		    hessz_ref(j, *i__) = hessz_ref(j, *i__) + term_ref(3, j);
		    hessx_ref(j, k) = hessx_ref(j, k) - term_ref(1, j);
		    hessy_ref(j, k) = hessy_ref(j, k) - term_ref(2, j);
		    hessz_ref(j, k) = hessz_ref(j, k) - term_ref(3, j);
		}
	    }
	}
    }

/*     for periodic boundary conditions with large cutoffs */
/*     neighbors must be found by the replicates method */

    if (! bound_1.use_replica__) {
	return 0;
    }

/*     calculate interaction energy with other unit cells */

    i__1 = charge_1.nion;
    for (kk = 1; kk <= i__1; ++kk) {
	k = charge_1.iion[kk - 1];
	kn = charge_1.jion[kk - 1];
	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, i__, &k, &c__0, &c__0, &c__0, &c__0);
	}

/*     compute the energy contribution for this interaction */

	if (proceed) {
	    i__2 = cell_1.ncell;
	    for (jcell = 1; jcell <= i__2; ++jcell) {
		xr = xi - atoms_1.x[k - 1];
		yr = yi - atoms_1.y[k - 1];
		zr = zi - atoms_1.z__[k - 1];
		imager_(&xr, &yr, &zr, &jcell);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    r__ = sqrt(r2);
		    rb = r__ + chgpot_1.ebuffer;
		    rb2 = rb * rb;
		    fik = fi * charge_1.pchg[kk - 1];
		    if (bound_1.use_polymer__) {
			if (r2 <= bound_1.polycut2) {
			    fik *= cscale[kn - 1];
			}
		    }

/*     compute chain rule terms for Hessian matrix elements */

		    de = -fik / rb2;
		    d2e = de * -2. / rb;

/*     use shifted energy switching if near the cutoff distance */

		    if (r2 > shunt_1.cut2) {
			e = fik / r__;
			shift = fik / ((shunt_1.off + shunt_1.cut) * .5);
			e -= shift;
			r3 = r2 * r__;
			r4 = r2 * r2;
			r5 = r2 * r3;
			r6 = r3 * r3;
			r7 = r3 * r4;
			taper = shunt_1.c5 * r5 + shunt_1.c4 * r4 + 
				shunt_1.c3 * r3 + shunt_1.c2 * r2 + 
				shunt_1.c1 * r__ + shunt_1.c0;
			dtaper = shunt_1.c5 * 5. * r4 + shunt_1.c4 * 4. * r3 
				+ shunt_1.c3 * 3. * r2 + shunt_1.c2 * 2. * 
				r__ + shunt_1.c1;
			d2taper = shunt_1.c5 * 20. * r3 + shunt_1.c4 * 12. * 
				r2 + shunt_1.c3 * 6. * r__ + shunt_1.c2 * 2.;
			trans = fik * (shunt_1.f7 * r7 + shunt_1.f6 * r6 + 
				shunt_1.f5 * r5 + shunt_1.f4 * r4 + 
				shunt_1.f3 * r3 + shunt_1.f2 * r2 + 
				shunt_1.f1 * r__ + shunt_1.f0);
			dtrans = fik * (shunt_1.f7 * 7. * r6 + shunt_1.f6 * 
				6. * r5 + shunt_1.f5 * 5. * r4 + shunt_1.f4 * 
				4. * r3 + shunt_1.f3 * 3. * r2 + shunt_1.f2 * 
				2. * r__ + shunt_1.f1);
			d2trans = fik * (shunt_1.f7 * 42. * r5 + shunt_1.f6 * 
				30. * r4 + shunt_1.f5 * 20. * r3 + shunt_1.f4 
				* 12. * r2 + shunt_1.f3 * 6. * r__ + 
				shunt_1.f2 * 2.);
			d2e = e * d2taper + de * 2. * dtaper + d2e * taper + 
				d2trans;
			de = e * dtaper + de * taper + dtrans;
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			de *= fgrp;
			d2e *= fgrp;
		    }

/*     form the individual Hessian element components */

		    de /= r__;
		    d2e = (d2e - de) / r2;
		    d2edx = d2e * xr;
		    d2edy = d2e * yr;
		    d2edz = d2e * zr;
		    term_ref(1, 1) = d2edx * xr + de;
		    term_ref(1, 2) = d2edx * yr;
		    term_ref(1, 3) = d2edx * zr;
		    term_ref(2, 1) = term_ref(1, 2);
		    term_ref(2, 2) = d2edy * yr + de;
		    term_ref(2, 3) = d2edy * zr;
		    term_ref(3, 1) = term_ref(1, 3);
		    term_ref(3, 2) = term_ref(2, 3);
		    term_ref(3, 3) = d2edz * zr + de;

/*     increment diagonal and non-diagonal Hessian elements */

		    for (j = 1; j <= 3; ++j) {
			hessx_ref(j, *i__) = hessx_ref(j, *i__) + term_ref(1, 
				j);
			hessy_ref(j, *i__) = hessy_ref(j, *i__) + term_ref(2, 
				j);
			hessz_ref(j, *i__) = hessz_ref(j, *i__) + term_ref(3, 
				j);
			hessx_ref(j, k) = hessx_ref(j, k) - term_ref(1, j);
			hessy_ref(j, k) = hessy_ref(j, k) - term_ref(2, j);
			hessz_ref(j, k) = hessz_ref(j, k) - term_ref(3, j);
		    }
		}
	    }
	}
    }
    return 0;
} /* echarge2a_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef term_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine echarge2b  --  Ewald summation charge Hessian  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "echarge2b" calculates second derivatives of the charge-charge */
/*     interaction energy for a single atom using a particle mesh */
/*     Ewald summation */

/*     note this version does not incorporate the reciprocal space */
/*     contribution to the Hessian calculation */


/* Subroutine */ int echarge2b_(integer *i__)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal r__, scaleterm, r2, de, fi, rb;
    static integer kk, in, kn;
    static doublereal xi, yi, zi, xr, yr, zr, d2e, rb2, fik, rew, cut2;
    extern doublereal erfc_(doublereal *);
    static doublereal fgrp, term[9]	/* was [3][3] */, d2edx, d2edy, d2edz;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal scale;
    static integer jcell;
    static doublereal cscale[25000];
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *), groups_(logical *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *);
    static logical proceed;
    static doublereal erfterm, expterm;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define term_ref(a_1,a_2) term[(a_2)*3 + a_1 - 4]
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




/*     first see if the atom of interest carries a charge */

    i__1 = charge_1.nion;
    for (k = 1; k <= i__1; ++k) {
	if (charge_1.iion[k - 1] == *i__) {
	    fi = chgpot_1.electric * charge_1.pchg[k - 1] / chgpot_1.dielec;
	    in = charge_1.jion[k - 1];
	    cut2 = cutoff_1.ewaldcut * cutoff_1.ewaldcut;
	    goto L10;
	}
    }
    return 0;
L10:

/*     store the coordinates of the atom of interest */

    xi = atoms_1.x[*i__ - 1];
    yi = atoms_1.y[*i__ - 1];
    zi = atoms_1.z__[*i__ - 1];

/*     set array needed to scale connected atom interactions */

    i__1 = charge_1.nion;
    for (j = 1; j <= i__1; ++j) {
	cscale[charge_1.iion[j - 1] - 1] = 1.;
    }
    i__1 = couple_1.n12[in - 1];
    for (j = 1; j <= i__1; ++j) {
	cscale[i12_ref(j, in) - 1] = chgpot_1.c2scale;
    }
    i__1 = couple_1.n13[in - 1];
    for (j = 1; j <= i__1; ++j) {
	cscale[i13_ref(j, in) - 1] = chgpot_1.c3scale;
    }
    i__1 = couple_1.n14[in - 1];
    for (j = 1; j <= i__1; ++j) {
	cscale[i14_ref(j, in) - 1] = chgpot_1.c4scale;
    }
    i__1 = couple_1.n15[in - 1];
    for (j = 1; j <= i__1; ++j) {
	cscale[i15_ref(j, in) - 1] = chgpot_1.c5scale;
    }

/*     calculate the real space Ewald interaction Hessian elements */

    i__1 = charge_1.nion;
    for (kk = 1; kk <= i__1; ++kk) {
	k = charge_1.iion[kk - 1];
	kn = charge_1.jion[kk - 1];
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, i__, &k, &c__0, &c__0, &c__0, &c__0);
	}
	proceed = TRUE_;
	if (proceed) {
	    proceed = kn != *i__;
	}

/*     compute the energy contribution for this interaction */

	if (proceed) {
	    xr = xi - atoms_1.x[k - 1];
	    yr = yi - atoms_1.y[k - 1];
	    zr = zi - atoms_1.z__[k - 1];
	    image_(&xr, &yr, &zr);
	    r2 = xr * xr + yr * yr + zr * zr;
	    if (r2 <= cut2) {
		r__ = sqrt(r2);
		rb = r__ + chgpot_1.ebuffer;
		rb2 = rb * rb;
		fik = fi * charge_1.pchg[kk - 1] * cscale[kn - 1];
		rew = ewald_1.aewald * r__;
		erfterm = erfc_(&rew);
/* Computing 2nd power */
		d__1 = rew;
		expterm = exp(-(d__1 * d__1));
		scale = cscale[kn - 1];
		if (group_1.use_group__) {
		    scale *= fgrp;
		}
		scaleterm = scale - 1.;

/*     compute chain rule terms for Hessian matrix elements */

		de = -fik * ((erfterm + scaleterm) / rb2 + ewald_1.aewald * 
			2. / 1.772453850905516027 * expterm / r__);
/* Computing 3rd power */
		d__1 = ewald_1.aewald;
		d2e = de * -2. / rb + fik / (rb * rb2) * 2. * scaleterm + fik 
			* 4. * (d__1 * (d__1 * d__1)) / 1.772453850905516027 *
			 expterm + fik / (rb * rb2) * 2. * scaleterm;

/*     form the individual Hessian element components */

		de /= r__;
		d2e = (d2e - de) / r2;
		d2edx = d2e * xr;
		d2edy = d2e * yr;
		d2edz = d2e * zr;
		term_ref(1, 1) = d2edx * xr + de;
		term_ref(1, 2) = d2edx * yr;
		term_ref(1, 3) = d2edx * zr;
		term_ref(2, 1) = term_ref(1, 2);
		term_ref(2, 2) = d2edy * yr + de;
		term_ref(2, 3) = d2edy * zr;
		term_ref(3, 1) = term_ref(1, 3);
		term_ref(3, 2) = term_ref(2, 3);
		term_ref(3, 3) = d2edz * zr + de;

/*     increment diagonal and non-diagonal Hessian elements */

		for (j = 1; j <= 3; ++j) {
		    hessx_ref(j, *i__) = hessx_ref(j, *i__) + term_ref(1, j);
		    hessy_ref(j, *i__) = hessy_ref(j, *i__) + term_ref(2, j);
		    hessz_ref(j, *i__) = hessz_ref(j, *i__) + term_ref(3, j);
		    hessx_ref(j, k) = hessx_ref(j, k) - term_ref(1, j);
		    hessy_ref(j, k) = hessy_ref(j, k) - term_ref(2, j);
		    hessz_ref(j, k) = hessz_ref(j, k) - term_ref(3, j);
		}
	    }
	}
    }

/*     for periodic boundary conditions with large cutoffs */
/*     neighbors must be found by the replicates method */

    if (! bound_1.use_replica__) {
	return 0;
    }

/*     calculate interaction energy with other unit cells */

    i__1 = charge_1.nion;
    for (kk = 1; kk <= i__1; ++kk) {
	k = charge_1.iion[kk - 1];
	kn = charge_1.jion[kk - 1];
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, i__, &k, &c__0, &c__0, &c__0, &c__0);
	}
	proceed = TRUE_;

/*     compute the energy contribution for this interaction */

	if (proceed) {
	    i__2 = cell_1.ncell;
	    for (jcell = 1; jcell <= i__2; ++jcell) {
		xr = xi - atoms_1.x[k - 1];
		yr = yi - atoms_1.y[k - 1];
		zr = zi - atoms_1.z__[k - 1];
		imager_(&xr, &yr, &zr, &jcell);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= cut2) {
		    r__ = sqrt(r2);
		    rb = r__ + chgpot_1.ebuffer;
		    rb2 = rb * rb;
		    fik = fi * charge_1.pchg[kk - 1];
		    rew = ewald_1.aewald * r__;
		    erfterm = erfc_(&rew);
/* Computing 2nd power */
		    d__1 = rew;
		    expterm = exp(-(d__1 * d__1));
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

/*     compute chain rule terms for Hessian matrix elements */

/* Computing 2nd power */
		    d__1 = rew;
		    de = -fik * ((erfterm + scaleterm) / rb2 + ewald_1.aewald 
			    * 2. / 1.772453850905516027 * exp(-(d__1 * d__1)) 
			    / r__);
/* Computing 3rd power */
		    d__1 = ewald_1.aewald;
		    d2e = de * -2. / rb + fik / (rb * rb2) * 2. * scaleterm + 
			    fik * 4. * (d__1 * (d__1 * d__1)) / 
			    1.772453850905516027 * expterm + fik / (rb * rb2) 
			    * 2. * scaleterm;

/*     form the individual Hessian element components */

		    de /= r__;
		    d2e = (d2e - de) / r2;
		    d2edx = d2e * xr;
		    d2edy = d2e * yr;
		    d2edz = d2e * zr;
		    term_ref(1, 1) = d2edx * xr + de;
		    term_ref(1, 2) = d2edx * yr;
		    term_ref(1, 3) = d2edx * zr;
		    term_ref(2, 1) = term_ref(1, 2);
		    term_ref(2, 2) = d2edy * yr + de;
		    term_ref(2, 3) = d2edy * zr;
		    term_ref(3, 1) = term_ref(1, 3);
		    term_ref(3, 2) = term_ref(2, 3);
		    term_ref(3, 3) = d2edz * zr + de;

/*     increment diagonal and non-diagonal Hessian elements */

		    for (j = 1; j <= 3; ++j) {
			hessx_ref(j, *i__) = hessx_ref(j, *i__) + term_ref(1, 
				j);
			hessy_ref(j, *i__) = hessy_ref(j, *i__) + term_ref(2, 
				j);
			hessz_ref(j, *i__) = hessz_ref(j, *i__) + term_ref(3, 
				j);
			hessx_ref(j, k) = hessx_ref(j, k) - term_ref(1, j);
			hessy_ref(j, k) = hessy_ref(j, k) - term_ref(2, j);
			hessz_ref(j, k) = hessz_ref(j, k) - term_ref(3, j);
		    }
		}
	    }
	}
    }
    return 0;
} /* echarge2b_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef term_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine echarge2c  --  charge Hessian for smoothing  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "echarge2c" calculates second derivatives of the charge-charge */
/*     interaction energy for a single atom for use with potential */
/*     smoothing methods */


/* Subroutine */ int echarge2c_(integer *i__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal r__, r2, de, fi, rb;
    static integer kk, in, kn;
    static doublereal xi, yi, zi, xr, yr, zr, d2e, rb2, fik;
    extern doublereal erf_(doublereal *);
    static doublereal fgrp, term[9]	/* was [3][3] */, d2edx, d2edy, d2edz,
	     width, wterm, width2, width3, cscale[25000], expcut;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;
    static doublereal erfterm, expterm;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define term_ref(a_1,a_2) term[(a_2)*3 + a_1 - 4]
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




/*     first see if the atom of interest carries a charge */

    i__1 = charge_1.nion;
    for (k = 1; k <= i__1; ++k) {
	if (charge_1.iion[k - 1] == *i__) {
	    fi = chgpot_1.electric * charge_1.pchg[k - 1] / chgpot_1.dielec;
	    in = charge_1.jion[k - 1];
	    goto L10;
	}
    }
    return 0;
L10:

/*     store the coordinates of the atom of interest */

    xi = atoms_1.x[*i__ - 1];
    yi = atoms_1.y[*i__ - 1];
    zi = atoms_1.z__[*i__ - 1];

/*     set array needed to scale connected atom interactions */

    i__1 = charge_1.nion;
    for (j = 1; j <= i__1; ++j) {
	cscale[charge_1.iion[j - 1] - 1] = 1.;
    }
    i__1 = couple_1.n12[in - 1];
    for (j = 1; j <= i__1; ++j) {
	cscale[i12_ref(j, in) - 1] = chgpot_1.c2scale;
    }
    i__1 = couple_1.n13[in - 1];
    for (j = 1; j <= i__1; ++j) {
	cscale[i13_ref(j, in) - 1] = chgpot_1.c3scale;
    }
    i__1 = couple_1.n14[in - 1];
    for (j = 1; j <= i__1; ++j) {
	cscale[i14_ref(j, in) - 1] = chgpot_1.c4scale;
    }
    i__1 = couple_1.n15[in - 1];
    for (j = 1; j <= i__1; ++j) {
	cscale[i15_ref(j, in) - 1] = chgpot_1.c5scale;
    }

/*     set the smallest exponential terms to be calculated */

    expcut = -50.;

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

/*     calculate the charge interaction energy Hessian elements */

    i__1 = charge_1.nion;
    for (kk = 1; kk <= i__1; ++kk) {
	k = charge_1.iion[kk - 1];
	kn = charge_1.jion[kk - 1];
	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, i__, &k, &c__0, &c__0, &c__0, &c__0);
	}
	if (proceed) {
	    proceed = kn != *i__;
	}

/*     compute the energy contribution for this interaction */

	if (proceed) {
	    xr = xi - atoms_1.x[k - 1];
	    yr = yi - atoms_1.y[k - 1];
	    zr = zi - atoms_1.z__[k - 1];
	    r2 = xr * xr + yr * yr + zr * zr;
	    r__ = sqrt(r2);
	    rb = r__ + chgpot_1.ebuffer;
	    rb2 = rb * rb;
	    fik = fi * charge_1.pchg[kk - 1] * cscale[kn - 1];

/*     compute chain rule terms for Hessian matrix elements */

	    de = -fik / rb2;
	    d2e = de * -2. / rb;

/*     transform the potential function via smoothing */

	    if (warp_1.use_dem__) {
		if (width > 0.) {
		    d__1 = width * rb;
		    erfterm = erf_(&d__1);
		    expterm = -rb2 * width2;
		    if (expterm > expcut) {
			expterm = fik * 2. * width * exp(expterm) / (rb * 
				1.772453850905516027);
		    } else {
			expterm = 0.;
		    }
		    de = de * erfterm + expterm;
		    d2e = (de / rb + expterm * rb * width2) * -2.;
		}
	    } else if (warp_1.use_gda__) {
		width = warp_1.m2[*i__ - 1] + warp_1.m2[k - 1];
		if (width > 0.) {
		    width = wterm / sqrt(width);
		    width2 = width * width;
		    d__1 = width * rb;
		    erfterm = erf_(&d__1);
		    expterm = -rb2 * width2;
		    if (expterm > expcut) {
			expterm = fik * 2. * width * exp(expterm) / (rb * 
				1.772453850905516027);
		    } else {
			expterm = 0.;
		    }
		    de = de * erfterm + expterm;
		    d2e = (de / rb + expterm * r__ * width2) * -2.;
		}
	    } else if (warp_1.use_tophat__) {
		if (width > rb) {
		    d2e = -fik / width3;
		    de = d2e * rb;
		}
	    } else if (warp_1.use_stophat__) {
		wterm = rb + width;
		de = -fik / (wterm * wterm);
		d2e = de * -2. / wterm;
	    }

/*     scale the interaction based on its group membership */

	    if (group_1.use_group__) {
		de *= fgrp;
		d2e *= fgrp;
	    }

/*     form the individual Hessian element components */

	    de /= r__;
	    d2e = (d2e - de) / r2;
	    d2edx = d2e * xr;
	    d2edy = d2e * yr;
	    d2edz = d2e * zr;
	    term_ref(1, 1) = d2edx * xr + de;
	    term_ref(1, 2) = d2edx * yr;
	    term_ref(1, 3) = d2edx * zr;
	    term_ref(2, 1) = term_ref(1, 2);
	    term_ref(2, 2) = d2edy * yr + de;
	    term_ref(2, 3) = d2edy * zr;
	    term_ref(3, 1) = term_ref(1, 3);
	    term_ref(3, 2) = term_ref(2, 3);
	    term_ref(3, 3) = d2edz * zr + de;

/*     increment diagonal and non-diagonal Hessian elements */

	    for (j = 1; j <= 3; ++j) {
		hessx_ref(j, *i__) = hessx_ref(j, *i__) + term_ref(1, j);
		hessy_ref(j, *i__) = hessy_ref(j, *i__) + term_ref(2, j);
		hessz_ref(j, *i__) = hessz_ref(j, *i__) + term_ref(3, j);
		hessx_ref(j, k) = hessx_ref(j, k) - term_ref(1, j);
		hessy_ref(j, k) = hessy_ref(j, k) - term_ref(2, j);
		hessz_ref(j, k) = hessz_ref(j, k) - term_ref(3, j);
	    }
	}
    }
    return 0;
} /* echarge2c_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef term_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref


