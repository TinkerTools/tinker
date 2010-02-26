/* esolv2.f -- translated by f2c (version 20050501).
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
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

struct {
    doublereal rsolv[25000], asolv[25000], rborn[25000], drb[25000], drbp[
	    25000], drobc[25000], doffset, p1, p2, p3, p4, p5, gpol[25000], 
	    shct[25000], aobc[25000], bobc[25000], gobc[25000], vsolv[25000], 
	    wace[1000000]	/* was [1000][1000] */, s2ace[1000000]	/* 
	    was [1000][1000] */, uace[1000000]	/* was [1000][1000] */;
    char solvtyp[8], borntyp[8];
} solute_;

#define solute_1 solute_

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
    doublereal hessx[75000]	/* was [3][25000] */, hessy[75000]	/* 
	    was [3][25000] */, hessz[75000]	/* was [3][25000] */;
} hessn_;

#define hessn_1 hessn_

struct {
    doublereal off, off2, cut, cut2, c0, c1, c2, c3, c4, c5, f0, f1, f2, f3, 
	    f4, f5, f6, f7;
} shunt_;

#define shunt_1 shunt_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine esolv2  --  atom-by-atom solvation Hessian  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "esolv2" calculates second derivatives of the continuum */
/*     solvation energy for surface area, generalized Born, */
/*     generalized Kirkwood and Poisson-Boltzmann solvation models */

/*     note this version does not contain the chain rule terms */
/*     for derivatives of Born radii with respect to coordinates */


/* Subroutine */ int esolv2_(integer *i__)
{
    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int egb2a_(integer *), egb2b_(integer *);
    static doublereal probe;



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
/*     ##  solute.i  --  parameters for continuum solvation models  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     rsolv     atomic radius of each atom for continuum solvation */
/*     asolv     atomic surface area solvation parameters */
/*     rborn     Born radius of each atom for GB/SA solvation */
/*     drb       solvation derivatives with respect to Born radii */
/*     drbp      GK polarization derivatives with respect to Born radii */
/*     drobc     chain rule term for Onufriev-Bashford-Case radii */
/*     doffset   dielectric offset to continuum solvation atomic radii */
/*     p1        single-atom scale factor for analytical Still radii */
/*     p2        1-2 interaction scale factor for analytical Still radii */
/*     p3        1-3 interaction scale factor for analytical Still radii */
/*     p4        nonbonded scale factor for analytical Still radii */
/*     p5        soft cutoff parameter for analytical Still radii */
/*     gpol      polarization self-energy values for each atom */
/*     shct      overlap scale factors for Hawkins-Cramer-Truhlar radii */
/*     aobc      alpha values for Onufriev-Bashford-Case radii */
/*     bobc      beta values for Onufriev-Bashford-Case radii */
/*     gobc      gamma values for Onufriev-Bashford-Case radii */
/*     vsolv     atomic volume of each atom for use with ACE */
/*     wace      "omega" values for atom class pairs for use with ACE */
/*     s2ace     "sigma^2" values for atom class pairs for use with ACE */
/*     uace      "mu" values for atom class pairs for use with ACE */
/*     solvtyp   type of continuum solvation energy model in use */
/*     borntyp   method to be used for the Born radius computation */




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


/*     real*8 aes(maxatm) */
/*     real*8 des(3,maxatm) */


/*     set a value for the solvent molecule probe radius */

    probe = 1.4;

/*     compute the surface area-based solvation energy term */

/*     call surface1 (es,aes,des,rsolv,asolv,probe) */

/*     get the generalized Born term for GB/SA solvation */

    if (potent_1.use_born__ && s_cmp(solute_1.solvtyp, "GK", (ftnlen)8, (
	    ftnlen)2) != 0 && s_cmp(solute_1.solvtyp, "PB", (ftnlen)8, (
	    ftnlen)2) != 0) {
	if (warp_1.use_smooth__) {
	    egb2b_(i__);
	} else {
	    egb2a_(i__);
	}
    }
    return 0;
} /* esolv2_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine egb2a  --  atom-by-atom GB solvation Hessian  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "egb2a" calculates second derivatives of the generalized */
/*     Born energy term for the GB/SA solvation models */


/* Subroutine */ int egb2a_(integer *i__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal e;
    static integer j, k;
    static doublereal r__, r2, r3, r4, r5, r6, r7, de, fi;
    static integer kk;
    static doublereal xi, yi, zi, xr, yr, zr, d2e, rb2, rm2, fgb, fik, fgb2, 
	    dfgb, term[9]	/* was [3][3] */, dfgb2, d2fgb, d2edx, d2edy, 
	    d2edz, taper, shift, trans, dtaper, dwater, dtrans;
    extern /* Subroutine */ int switch_(char *, ftnlen);
    static doublereal d2taper, d2trans, expterm;


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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  solute.i  --  parameters for continuum solvation models  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     rsolv     atomic radius of each atom for continuum solvation */
/*     asolv     atomic surface area solvation parameters */
/*     rborn     Born radius of each atom for GB/SA solvation */
/*     drb       solvation derivatives with respect to Born radii */
/*     drbp      GK polarization derivatives with respect to Born radii */
/*     drobc     chain rule term for Onufriev-Bashford-Case radii */
/*     doffset   dielectric offset to continuum solvation atomic radii */
/*     p1        single-atom scale factor for analytical Still radii */
/*     p2        1-2 interaction scale factor for analytical Still radii */
/*     p3        1-3 interaction scale factor for analytical Still radii */
/*     p4        nonbonded scale factor for analytical Still radii */
/*     p5        soft cutoff parameter for analytical Still radii */
/*     gpol      polarization self-energy values for each atom */
/*     shct      overlap scale factors for Hawkins-Cramer-Truhlar radii */
/*     aobc      alpha values for Onufriev-Bashford-Case radii */
/*     bobc      beta values for Onufriev-Bashford-Case radii */
/*     gobc      gamma values for Onufriev-Bashford-Case radii */
/*     vsolv     atomic volume of each atom for use with ACE */
/*     wace      "omega" values for atom class pairs for use with ACE */
/*     s2ace     "sigma^2" values for atom class pairs for use with ACE */
/*     uace      "mu" values for atom class pairs for use with ACE */
/*     solvtyp   type of continuum solvation energy model in use */
/*     borntyp   method to be used for the Born radius computation */




/*     first see if the atom of interest carries a charge */

    i__1 = charge_1.nion;
    for (k = 1; k <= i__1; ++k) {
	if (charge_1.iion[k - 1] == *i__) {
	    fi = charge_1.pchg[k - 1];
	    goto L10;
	}
    }
    return 0;
L10:

/*     store the coordinates of the atom of interest */

    xi = atoms_1.x[*i__ - 1];
    yi = atoms_1.y[*i__ - 1];
    zi = atoms_1.z__[*i__ - 1];

/*     set the solvent dielectric and energy conversion factor */

    dwater = 78.3;
    fi = -chgpot_1.electric * (1. - 1. / dwater) * fi;

/*     set cutoff distances and switching function coefficients */

    switch_("CHARGE", (ftnlen)6);

/*     calculate GB polarization energy Hessian elements */

    i__1 = charge_1.nion;
    for (kk = 1; kk <= i__1; ++kk) {
	k = charge_1.iion[kk - 1];
	if (*i__ != k) {
	    xr = xi - atoms_1.x[k - 1];
	    yr = yi - atoms_1.y[k - 1];
	    zr = zi - atoms_1.z__[k - 1];
	    r2 = xr * xr + yr * yr + zr * zr;
	    if (r2 <= shunt_1.off2) {
		r__ = sqrt(r2);
		fik = fi * charge_1.pchg[kk - 1];

/*     compute chain rule terms for Hessian matrix elements */

		rb2 = solute_1.rborn[*i__ - 1] * solute_1.rborn[k - 1];
		expterm = exp(r2 * -.25 / rb2);
		fgb2 = r2 + rb2 * expterm;
		fgb = sqrt(fgb2);
		dfgb = (1. - expterm * .25) * r__ / fgb;
		dfgb2 = dfgb * dfgb;
		d2fgb = -dfgb2 / fgb + dfgb / r__ + r2 / rb2 * .125 * expterm 
			/ fgb;
		de = -fik * dfgb / fgb2;
		d2e = -fik * (d2fgb - dfgb2 * 2. / fgb) / fgb2;

/*     use energy switching if near the cutoff distance */

		if (r2 > shunt_1.cut2) {
		    e = fik / fgb;
/* Computing 2nd power */
		    d__1 = (shunt_1.off + shunt_1.cut) * .5;
		    rm2 = d__1 * d__1;
		    shift = fik / sqrt(rm2 + rb2 * exp(rm2 * -.25 / rb2));
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
    return 0;
} /* egb2a_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef term_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine egb2b  --  GB solvation Hessian for smoothing  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "egb2b" calculates second derivatives of the generalized */
/*     Born energy term for the GB/SA solvation models for use with */
/*     potential smoothing methods */


/* Subroutine */ int egb2b_(integer *i__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal r__, r2, de, fi;
    static integer kk;
    static doublereal xi, yi, zi, xr, yr, zr, d2e, rb2, fgb, fik;
    extern doublereal erf_(doublereal *);
    static doublereal fgb2, dfgb, term[9]	/* was [3][3] */, dfgb2, 
	    d2fgb, d2edx, d2edy, d2edz, width, sterm, dwater, erfterm, 
	    expterm;


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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  solute.i  --  parameters for continuum solvation models  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     rsolv     atomic radius of each atom for continuum solvation */
/*     asolv     atomic surface area solvation parameters */
/*     rborn     Born radius of each atom for GB/SA solvation */
/*     drb       solvation derivatives with respect to Born radii */
/*     drbp      GK polarization derivatives with respect to Born radii */
/*     drobc     chain rule term for Onufriev-Bashford-Case radii */
/*     doffset   dielectric offset to continuum solvation atomic radii */
/*     p1        single-atom scale factor for analytical Still radii */
/*     p2        1-2 interaction scale factor for analytical Still radii */
/*     p3        1-3 interaction scale factor for analytical Still radii */
/*     p4        nonbonded scale factor for analytical Still radii */
/*     p5        soft cutoff parameter for analytical Still radii */
/*     gpol      polarization self-energy values for each atom */
/*     shct      overlap scale factors for Hawkins-Cramer-Truhlar radii */
/*     aobc      alpha values for Onufriev-Bashford-Case radii */
/*     bobc      beta values for Onufriev-Bashford-Case radii */
/*     gobc      gamma values for Onufriev-Bashford-Case radii */
/*     vsolv     atomic volume of each atom for use with ACE */
/*     wace      "omega" values for atom class pairs for use with ACE */
/*     s2ace     "sigma^2" values for atom class pairs for use with ACE */
/*     uace      "mu" values for atom class pairs for use with ACE */
/*     solvtyp   type of continuum solvation energy model in use */
/*     borntyp   method to be used for the Born radius computation */




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
	    fi = charge_1.pchg[k - 1];
	    goto L10;
	}
    }
    return 0;
L10:

/*     store the coordinates of the atom of interest */

    xi = atoms_1.x[*i__ - 1];
    yi = atoms_1.y[*i__ - 1];
    zi = atoms_1.z__[*i__ - 1];

/*     set the solvent dielectric and energy conversion factor */

    dwater = 78.3;
    fi = -chgpot_1.electric * (1. - 1. / dwater) * fi;

/*     set the extent of smoothing to be performed */

    sterm = .5 / sqrt(warp_1.diffc);

/*     calculate GB polarization energy Hessian elements */

    i__1 = charge_1.nion;
    for (kk = 1; kk <= i__1; ++kk) {
	k = charge_1.iion[kk - 1];
	if (*i__ != k) {
	    xr = xi - atoms_1.x[k - 1];
	    yr = yi - atoms_1.y[k - 1];
	    zr = zi - atoms_1.z__[k - 1];
	    r2 = xr * xr + yr * yr + zr * zr;
	    r__ = sqrt(r2);
	    fik = fi * charge_1.pchg[kk - 1];

/*     compute chain rule terms for Hessian matrix elements */

	    rb2 = solute_1.rborn[*i__ - 1] * solute_1.rborn[k - 1];
	    expterm = exp(r2 * -.25 / rb2);
	    fgb2 = r2 + rb2 * expterm;
	    fgb = sqrt(fgb2);
	    dfgb = (1. - expterm * .25) * r__ / fgb;
	    dfgb2 = dfgb * dfgb;
	    d2fgb = -dfgb2 / fgb + dfgb / r__ + r2 / rb2 * .125 * expterm / 
		    fgb;
	    de = -fik * dfgb / fgb2;
	    d2e = -fik * (d2fgb - dfgb2 * 2. / fgb) / fgb2;

/*     use a smoothable GB analogous to the Coulomb solution */

	    if (warp_1.deform > 0.) {
		width = warp_1.deform + rb2 * .15 * exp(rb2 * -.006 / 
			warp_1.deform);
		width = sterm / sqrt(width);
		d__1 = width * fgb;
		erfterm = erf_(&d__1);
/* Computing 2nd power */
		d__1 = width * fgb;
		expterm = width * exp(-(d__1 * d__1)) / 1.772453850905516027;
		de *= erfterm - expterm * 2. * fgb;
/* Computing 2nd power */
		d__1 = width;
		d2e = d2e * erfterm + fik * 2. * expterm * (d2fgb / fgb - 
			dfgb2 * 2. * (1. / fgb2 + d__1 * d__1));
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
} /* egb2b_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef term_ref


