/* echgdpl2.f -- translated by f2c (version 20050501).
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
    doublereal bdpl[50000], sdpl[50000];
    integer ndipole, idpl[100000]	/* was [2][50000] */;
} dipole_;

#define dipole_1 dipole_

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

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine echgdpl2  --  atomwise charge-dipole Hessian  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "echgdpl2" calculates second derivatives of the */
/*     charge-dipole interaction energy for a single atom */


/* Subroutine */ int echgdpl2_(integer *i__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal dtxkdxi1, dtykdxi1, dtxkdxk1, dtxkdxk2, dtykdxk1, 
	    dtykdxk2, dtzkdxi1, dtzkdxk1, dtzkdxk2, dtxkdyi1, dtxkdyk1, 
	    dtxkdyk2, dtykdyi1, dtykdyk1, dtykdyk2, dtzkdyi1, dtzkdyk1, 
	    dtzkdyk2, dtxkdzi1, dtxkdzk1, dtxkdzk2, dtykdzi1, dtykdzk1, 
	    dtykdzk2, dtzkdzi1, dtzkdzk1, dtzkdzk2, e, f;
    static integer k;
    static doublereal r__, d2taperxx, d2taperxy, d2taperyy, d2taperxz, 
	    d2taperzz, d2taperyz;
    static integer i1, k1, k2;
    static doublereal r2, r3, r4, r5, fi;
    static integer ii;
    static doublereal fk, xi, yi, zi, xk, yk, zk, xq, xr, yr, zr, yq, zq, sk1,
	     sk2, rk2, fik, rkr3, xrr2, yrr2, zrr2, fgrp;
    static integer skip[25000], omit[25000];
    static doublereal dotk, term, part, dotk2, part2, term2, xkrk2, ykrk2, 
	    zkrk2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static integer jcell;
    static doublereal taper, termx, termy, termz, dedxi1, dedyi1, dedzi1, 
	    dedxk1, dedyk1, dedzk1, dedxk2, dedyk2, dedzk2, dtdxi1, dtdyi1, 
	    dtdzi1, dtdxk1, dtdyk1, dtdzk1, dtdxk2, dtdyk2, dotkr2, dtdzk2;
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal factor, dtaper;
    extern /* Subroutine */ int switch_(char *, ftnlen);
    static doublereal termxk, termyk, termzk;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal d2taper, dotkrk2, dtxdxi1, dtydxi1, dtxdxk1, dtxdxk2, 
	    dtydxk1, dtydxk2, dtzdxi1, dtzdxk1, dtzdxk2, dtxdyi1, dtxdyk1, 
	    dtxdyk2, dtydyi1, dtydyk1, dtydyk2, dtzdyi1, dtzdyk1, dtzdyk2, 
	    dtxdzi1, dtxdzk1, dtxdzk2, dtydzi1, dtydzk1, dtydzk2, dtzdzi1, 
	    dtzdzk1, dtzdzk2;
    static logical proceed;
    static doublereal factork, dtaperx, dtapery, dtaperz;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define idpl_ref(a_1,a_2) dipole_1.idpl[(a_2)*2 + a_1 - 3]
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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  dipole.i  --  atom & bond dipoles for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     bdpl      magnitude of each of the dipoles (Debyes) */
/*     sdpl      position of each dipole between defining atoms */
/*     ndipole   total number of dipoles in the system */
/*     idpl      numbers of atoms that define each dipole */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  units.i  --  physical constants and unit conversions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     literature reference: */

/*     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended */
/*     Values of the Fundamental Physical Constants: 2006", Reviews of */
/*     Modern Physics, 80, 633-730 (2008) */

/*     The "2006 CODATA Recommended Values" are also available from */
/*     the NIST Reference on Constants, Units, and Uncertainty site */
/*     at http://physics.nist.gov/cuu/index.html */

/*     Most values below are derived from the 2006 CODATA reference */

/*     The conversion from calorie to Joule is the definition of the */
/*     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992) */

/*     The "coulomb" energy conversion factor is found by dimensional */
/*     analysis of Coulomb's Law, ie, by dividing the square of the */
/*     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is */
/*     the permittivity of vacuum (the "electric constant"); note that */
/*     eps0 is typically given in F/m, equivalent to C**2/(J-m) */

/*     The approximate value used for the Debye, 3.33564 x 10-30 C-m, */
/*     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997) */

/*     The value of "prescon" is based on definition of 1 atmosphere */
/*     as 101325 Pa set by the 10th Conference Generale des Poids et */
/*     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3 */

/*     avogadro    Avogadro's number (N) in particles/mole */
/*     lightspd    speed of light in vacuum (c) in cm/ps */
/*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K */
/*     gasconst    ideal gas constant (R) in kcal/mole/K */
/*     emass       mass of an electron in atomic mass units */
/*     joule       conversion from calories to joules */
/*     convert     conversion from kcal to g*Ang**2/ps**2 */
/*     bohr        conversion from Bohrs to Angstroms */
/*     hartree     conversion from Hartree to kcal/mole */
/*     evolt       conversion from Hartree to electron-volts */
/*     efreq       conversion from Hartree to cm-1 */
/*     coulomb     conversion from electron**2/Ang to kcal/mole */
/*     debye       conversion from electron-Ang to Debyes */
/*     prescon     conversion from kcal/mole/Ang**3 to Atm */




/*     zero out the lists of atoms to be skipped */

    if (dipole_1.ndipole == 0 || charge_1.nion == 0) {
	return 0;
    }
    i__1 = atoms_1.n;
    for (k = 1; k <= i__1; ++k) {
	skip[k - 1] = 0;
	omit[k - 1] = 0;
    }

/*     set conversion factor and switching function coefficients */

    f = -chgpot_1.electric / (chgpot_1.dielec * 4.80321);
    switch_("CHGDPL", (ftnlen)6);

/*     first see if the atom of interest carries a charge */

    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
	i1 = charge_1.iion[ii - 1];
	if (i1 != *i__) {
	    goto L10;
	}
	skip[i1 - 1] = i1;
	i__2 = couple_1.n12[i1 - 1];
	for (k = 1; k <= i__2; ++k) {
	    skip[i12_ref(k, i1) - 1] = i1;
	}
	xi = atoms_1.x[i1 - 1];
	yi = atoms_1.y[i1 - 1];
	zi = atoms_1.z__[i1 - 1];
	fi = f * charge_1.pchg[ii - 1];

/*     decide whether to compute the current interaction */

	i__2 = dipole_1.ndipole;
	for (k = 1; k <= i__2; ++k) {
	    k1 = idpl_ref(1, k);
	    k2 = idpl_ref(2, k);
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i1, &k1, &k2, &c__0, &c__0, &c__0);
	    }
	    if (proceed) {
		proceed = skip[k1 - 1] != i1 && skip[k2 - 1] != i1;
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		sk1 = 1. - dipole_1.sdpl[k - 1];
		sk2 = dipole_1.sdpl[k - 1];
		xk = atoms_1.x[k1 - 1] - atoms_1.x[k2 - 1];
		yk = atoms_1.y[k1 - 1] - atoms_1.y[k2 - 1];
		zk = atoms_1.z__[k1 - 1] - atoms_1.z__[k2 - 1];
		xr = xi - atoms_1.x[k1 - 1] + xk * sk2;
		yr = yi - atoms_1.y[k1 - 1] + yk * sk2;
		zr = zi - atoms_1.z__[k1 - 1] + zk * sk2;
		image_(&xr, &yr, &zr);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    rk2 = xk * xk + yk * yk + zk * zk;
		    rkr3 = sqrt(rk2 * r2) * r2;
		    dotk = xk * xr + yk * yr + zk * zr;
		    fik = -fi * dipole_1.bdpl[k - 1];

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			fik *= fgrp;
		    }

/*     some abbreviations used in various chain rule terms */

		    xrr2 = xr / r2;
		    yrr2 = yr / r2;
		    zrr2 = zr / r2;
		    xkrk2 = xk / rk2;
		    ykrk2 = yk / rk2;
		    zkrk2 = zk / rk2;
		    dotk2 = dotk * 2.;
		    dotkr2 = dotk / r2;

/*     form the chain rule terms for first derivatives */

		    term = fik / rkr3;
		    term2 = dotk * -3.;
		    termx = term * (xk + xrr2 * term2);
		    termy = term * (yk + yrr2 * term2);
		    termz = term * (zk + zrr2 * term2);
		    termxk = term * (xr - dotk * xkrk2);
		    termyk = term * (yr - dotk * ykrk2);
		    termzk = term * (zr - dotk * zkrk2);

/*     use energy switching if near the cutoff distance */

		    if (r2 > shunt_1.cut2) {
			e = fik * dotk / rkr3;
			dedxi1 = termx;
			dedyi1 = termy;
			dedzi1 = termz;
			dedxk1 = -sk1 * termx + termxk;
			dedyk1 = -sk1 * termy + termyk;
			dedzk1 = -sk1 * termz + termzk;
			dedxk2 = -sk2 * termx - termxk;
			dedyk2 = -sk2 * termy - termyk;
			dedzk2 = -sk2 * termz - termzk;
			r__ = sqrt(r2);
			r3 = r2 * r__;
			r4 = r2 * r2;
			r5 = r2 * r3;
			taper = shunt_1.c5 * r5 + shunt_1.c4 * r4 + 
				shunt_1.c3 * r3 + shunt_1.c2 * r2 + 
				shunt_1.c1 * r__ + shunt_1.c0;
			dtaper = shunt_1.c5 * 5. * r4 + shunt_1.c4 * 4. * r3 
				+ shunt_1.c3 * 3. * r2 + shunt_1.c2 * 2. * 
				r__ + shunt_1.c1;
			d2taper = shunt_1.c5 * 20. * r3 + shunt_1.c4 * 12. * 
				r2 + shunt_1.c3 * 6. * r__ + shunt_1.c2 * 2.;
			dtaper /= r__;
			dtaperx = xr * dtaper;
			dtapery = yr * dtaper;
			dtaperz = zr * dtaper;
			d2taper = e * (d2taper - dtaper);
			dtaper = e * dtaper;
			d2taperxx = xr * xrr2 * d2taper + dtaper;
			d2taperxy = xr * yrr2 * d2taper;
			d2taperxz = xr * zrr2 * d2taper;
			d2taperyy = yr * yrr2 * d2taper + dtaper;
			d2taperyz = yr * zrr2 * d2taper;
			d2taperzz = zr * zrr2 * d2taper + dtaper;
			term *= taper;
			termx *= taper;
			termy *= taper;
			termz *= taper;
			termxk *= taper;
			termyk *= taper;
			termzk *= taper;
		    }

/*     chain rule terms for second derivative components */

		    dtdxi1 = xrr2 * -3.;
		    part = xk - dotk2 * xrr2;
		    factor = (dotkr2 + xrr2 * part) * -3.;
		    factork = 1. - xk * xkrk2;
		    dtxdxi1 = dtdxi1 * termx + term * factor;
		    dtxkdxi1 = dtdxi1 * termxk + term * factork;
		    factor = yrr2 * -3. * part;
		    factork = -yk * xkrk2;
		    dtydxi1 = dtdxi1 * termy + term * factor;
		    dtykdxi1 = dtdxi1 * termyk + term * factork;
		    factor = zrr2 * -3. * part;
		    factork = -zk * xkrk2;
		    dtzdxi1 = dtdxi1 * termz + term * factor;
		    dtzkdxi1 = dtdxi1 * termzk + term * factork;
		    dtdyi1 = yrr2 * -3.;
		    part = yk - dotk2 * yrr2;
		    factor = xrr2 * -3. * part;
		    factork = -xk * ykrk2;
		    dtxdyi1 = dtdyi1 * termx + term * factor;
		    dtxkdyi1 = dtdyi1 * termxk + term * factork;
		    factor = (dotkr2 + yrr2 * part) * -3.;
		    factork = 1. - yk * ykrk2;
		    dtydyi1 = dtdyi1 * termy + term * factor;
		    dtykdyi1 = dtdyi1 * termyk + term * factork;
		    factor = zrr2 * -3. * part;
		    factork = -zk * ykrk2;
		    dtzdyi1 = dtdyi1 * termz + term * factor;
		    dtzkdyi1 = dtdyi1 * termzk + term * factork;
		    dtdzi1 = zrr2 * -3.;
		    part = zk - dotk2 * zrr2;
		    factor = xrr2 * -3. * part;
		    factork = -xk * zkrk2;
		    dtxdzi1 = dtdzi1 * termx + term * factor;
		    dtxkdzi1 = dtdzi1 * termxk + term * factork;
		    factor = yrr2 * -3. * part;
		    factork = -yk * zkrk2;
		    dtydzi1 = dtdzi1 * termy + term * factor;
		    dtykdzi1 = dtdzi1 * termyk + term * factork;
		    factor = (dotkr2 + zrr2 * part) * -3.;
		    factork = 1. - zk * zkrk2;
		    dtzdzi1 = dtdzi1 * termz + term * factor;
		    dtzkdzi1 = dtdzi1 * termzk + term * factork;

/*     increment diagonal and off-diagonal Hessian elements */

		    hessx_ref(1, i1) = hessx_ref(1, i1) + dtxdxi1;
		    hessx_ref(2, i1) = hessx_ref(2, i1) + dtydxi1;
		    hessx_ref(3, i1) = hessx_ref(3, i1) + dtzdxi1;
		    hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * dtxdxi1 + 
			    dtxkdxi1;
		    hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * dtydxi1 + 
			    dtykdxi1;
		    hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * dtzdxi1 + 
			    dtzkdxi1;
		    hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * dtxdxi1 - 
			    dtxkdxi1;
		    hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * dtydxi1 - 
			    dtykdxi1;
		    hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * dtzdxi1 - 
			    dtzkdxi1;
		    hessy_ref(1, i1) = hessy_ref(1, i1) + dtxdyi1;
		    hessy_ref(2, i1) = hessy_ref(2, i1) + dtydyi1;
		    hessy_ref(3, i1) = hessy_ref(3, i1) + dtzdyi1;
		    hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * dtxdyi1 + 
			    dtxkdyi1;
		    hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * dtydyi1 + 
			    dtykdyi1;
		    hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * dtzdyi1 + 
			    dtzkdyi1;
		    hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * dtxdyi1 - 
			    dtxkdyi1;
		    hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * dtydyi1 - 
			    dtykdyi1;
		    hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * dtzdyi1 - 
			    dtzkdyi1;
		    hessz_ref(1, i1) = hessz_ref(1, i1) + dtxdzi1;
		    hessz_ref(2, i1) = hessz_ref(2, i1) + dtydzi1;
		    hessz_ref(3, i1) = hessz_ref(3, i1) + dtzdzi1;
		    hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * dtxdzi1 + 
			    dtxkdzi1;
		    hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * dtydzi1 + 
			    dtykdzi1;
		    hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * dtzdzi1 + 
			    dtzkdzi1;
		    hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * dtxdzi1 - 
			    dtxkdzi1;
		    hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * dtydzi1 - 
			    dtykdzi1;
		    hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * dtzdzi1 - 
			    dtzkdzi1;

/*     more energy switching if near the cutoff distance */

		    if (r2 > shunt_1.cut2) {
			hessx_ref(1, i1) = hessx_ref(1, i1) + dtaperx * 
				dedxi1 + dtaperx * dedxi1 + d2taperxx;
			hessx_ref(2, i1) = hessx_ref(2, i1) + dtaperx * 
				dedyi1 + dtapery * dedxi1 + d2taperxy;
			hessx_ref(3, i1) = hessx_ref(3, i1) + dtaperx * 
				dedzi1 + dtaperz * dedxi1 + d2taperxz;
			hessx_ref(1, k1) = hessx_ref(1, k1) + dtaperx * 
				dedxk1 - sk1 * (dtaperx * dedxi1 + d2taperxx);
			hessx_ref(2, k1) = hessx_ref(2, k1) + dtaperx * 
				dedyk1 - sk1 * (dtapery * dedxi1 + d2taperxy);
			hessx_ref(3, k1) = hessx_ref(3, k1) + dtaperx * 
				dedzk1 - sk1 * (dtaperz * dedxi1 + d2taperxz);
			hessx_ref(1, k2) = hessx_ref(1, k2) + dtaperx * 
				dedxk2 - sk2 * (dtaperx * dedxi1 + d2taperxx);
			hessx_ref(2, k2) = hessx_ref(2, k2) + dtaperx * 
				dedyk2 - sk2 * (dtapery * dedxi1 + d2taperxy);
			hessx_ref(3, k2) = hessx_ref(3, k2) + dtaperx * 
				dedzk2 - sk2 * (dtaperz * dedxi1 + d2taperxz);
			hessy_ref(1, i1) = hessy_ref(1, i1) + dtapery * 
				dedxi1 + dtaperx * dedyi1 + d2taperxy;
			hessy_ref(2, i1) = hessy_ref(2, i1) + dtapery * 
				dedyi1 + dtapery * dedyi1 + d2taperyy;
			hessy_ref(3, i1) = hessy_ref(3, i1) + dtapery * 
				dedzi1 + dtaperz * dedyi1 + d2taperyz;
			hessy_ref(1, k1) = hessy_ref(1, k1) + dtapery * 
				dedxk1 - sk1 * (dtaperx * dedyi1 + d2taperxy);
			hessy_ref(2, k1) = hessy_ref(2, k1) + dtapery * 
				dedyk1 - sk1 * (dtapery * dedyi1 + d2taperyy);
			hessy_ref(3, k1) = hessy_ref(3, k1) + dtapery * 
				dedzk1 - sk1 * (dtaperz * dedyi1 + d2taperyz);
			hessy_ref(1, k2) = hessy_ref(1, k2) + dtapery * 
				dedxk2 - sk2 * (dtaperx * dedyi1 + d2taperxy);
			hessy_ref(2, k2) = hessy_ref(2, k2) + dtapery * 
				dedyk2 - sk2 * (dtapery * dedyi1 + d2taperyy);
			hessy_ref(3, k2) = hessy_ref(3, k2) + dtapery * 
				dedzk2 - sk2 * (dtaperz * dedyi1 + d2taperyz);
			hessz_ref(1, i1) = hessz_ref(1, i1) + dtaperz * 
				dedxi1 + dtaperx * dedzi1 + d2taperxz;
			hessz_ref(2, i1) = hessz_ref(2, i1) + dtaperz * 
				dedyi1 + dtapery * dedzi1 + d2taperyz;
			hessz_ref(3, i1) = hessz_ref(3, i1) + dtaperz * 
				dedzi1 + dtaperz * dedzi1 + d2taperzz;
			hessz_ref(1, k1) = hessz_ref(1, k1) + dtaperz * 
				dedxk1 - sk1 * (dtaperx * dedzi1 + d2taperxz);
			hessz_ref(2, k1) = hessz_ref(2, k1) + dtaperz * 
				dedyk1 - sk1 * (dtapery * dedzi1 + d2taperyz);
			hessz_ref(3, k1) = hessz_ref(3, k1) + dtaperz * 
				dedzk1 - sk1 * (dtaperz * dedzi1 + d2taperzz);
			hessz_ref(1, k2) = hessz_ref(1, k2) + dtaperz * 
				dedxk2 - sk2 * (dtaperx * dedzi1 + d2taperxz);
			hessz_ref(2, k2) = hessz_ref(2, k2) + dtaperz * 
				dedyk2 - sk2 * (dtapery * dedzi1 + d2taperyz);
			hessz_ref(3, k2) = hessz_ref(3, k2) + dtaperz * 
				dedzk2 - sk2 * (dtaperz * dedzi1 + d2taperzz);
		    }
		}
	    }
	}
L10:
	;
    }

/*     see if the atom of interest is part of a dipole */

    i__1 = dipole_1.ndipole;
    for (k = 1; k <= i__1; ++k) {
	k1 = idpl_ref(1, k);
	k2 = idpl_ref(2, k);
	if (k1 != *i__ && k2 != *i__) {
	    goto L20;
	}
	i__2 = couple_1.n12[k1 - 1];
	for (ii = 1; ii <= i__2; ++ii) {
	    omit[i12_ref(ii, k1) - 1] = k;
	}
	i__2 = couple_1.n12[k2 - 1];
	for (ii = 1; ii <= i__2; ++ii) {
	    omit[i12_ref(ii, k2) - 1] = k;
	}
	sk1 = 1. - dipole_1.sdpl[k - 1];
	sk2 = dipole_1.sdpl[k - 1];
	xk = atoms_1.x[k1 - 1] - atoms_1.x[k2 - 1];
	yk = atoms_1.y[k1 - 1] - atoms_1.y[k2 - 1];
	zk = atoms_1.z__[k1 - 1] - atoms_1.z__[k2 - 1];
	rk2 = xk * xk + yk * yk + zk * zk;
	xq = atoms_1.x[k1 - 1] - xk * sk2;
	yq = atoms_1.y[k1 - 1] - yk * sk2;
	zq = atoms_1.z__[k1 - 1] - zk * sk2;
	fk = -f * dipole_1.bdpl[k - 1];

/*     decide whether to compute the current interaction */

	i__2 = charge_1.nion;
	for (ii = 1; ii <= i__2; ++ii) {
	    i1 = charge_1.iion[ii - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i1, &k1, &k2, &c__0, &c__0, &c__0);
	    }
	    if (proceed) {
		proceed = omit[i1 - 1] != k;
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xr = atoms_1.x[i1 - 1] - xq;
		yr = atoms_1.y[i1 - 1] - yq;
		zr = atoms_1.z__[i1 - 1] - zq;
		image_(&xr, &yr, &zr);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    rkr3 = sqrt(rk2 * r2) * r2;
		    dotk = xk * xr + yk * yr + zk * zr;
		    fik = fk * charge_1.pchg[ii - 1];

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			fik *= fgrp;
		    }

/*     some abbreviations used in various chain rule terms */

		    xrr2 = xr / r2;
		    yrr2 = yr / r2;
		    zrr2 = zr / r2;
		    xkrk2 = xk / rk2;
		    ykrk2 = yk / rk2;
		    zkrk2 = zk / rk2;
		    dotk2 = dotk * 2.;
		    dotkr2 = dotk / r2;
		    dotkrk2 = dotk / rk2;

/*     form the chain rule terms for first derivatives */

		    term = fik / rkr3;
		    term2 = dotk * -3.;
		    termx = term * (xk + xrr2 * term2);
		    termy = term * (yk + yrr2 * term2);
		    termz = term * (zk + zrr2 * term2);
		    termxk = term * (xr - dotk * xkrk2);
		    termyk = term * (yr - dotk * ykrk2);
		    termzk = term * (zr - dotk * zkrk2);

/*     use energy switching if near the cutoff distance */

		    if (r2 > shunt_1.cut2) {
			e = fik * dotk / rkr3;
			dedxi1 = termx;
			dedyi1 = termy;
			dedzi1 = termz;
			dedxk1 = -sk1 * termx + termxk;
			dedyk1 = -sk1 * termy + termyk;
			dedzk1 = -sk1 * termz + termzk;
			dedxk2 = -sk2 * termx - termxk;
			dedyk2 = -sk2 * termy - termyk;
			dedzk2 = -sk2 * termz - termzk;
			r__ = sqrt(r2);
			r3 = r2 * r__;
			r4 = r2 * r2;
			r5 = r2 * r3;
			taper = shunt_1.c5 * r5 + shunt_1.c4 * r4 + 
				shunt_1.c3 * r3 + shunt_1.c2 * r2 + 
				shunt_1.c1 * r__ + shunt_1.c0;
			dtaper = shunt_1.c5 * 5. * r4 + shunt_1.c4 * 4. * r3 
				+ shunt_1.c3 * 3. * r2 + shunt_1.c2 * 2. * 
				r__ + shunt_1.c1;
			d2taper = shunt_1.c5 * 20. * r3 + shunt_1.c4 * 12. * 
				r2 + shunt_1.c3 * 6. * r__ + shunt_1.c2 * 2.;
			dtaper /= r__;
			dtaperx = xr * dtaper;
			dtapery = yr * dtaper;
			dtaperz = zr * dtaper;
			d2taper = e * (d2taper - dtaper);
			dtaper = e * dtaper;
			d2taperxx = xr * xrr2 * d2taper + dtaper;
			d2taperxy = xr * yrr2 * d2taper;
			d2taperxz = xr * zrr2 * d2taper;
			d2taperyy = yr * yrr2 * d2taper + dtaper;
			d2taperyz = yr * zrr2 * d2taper;
			d2taperzz = zr * zrr2 * d2taper + dtaper;
			term *= taper;
			termx *= taper;
			termy *= taper;
			termz *= taper;
			termxk *= taper;
			termyk *= taper;
			termzk *= taper;
		    }

/*     chain rule terms for second derivative components */

		    if (k1 == *i__) {
			dtdxk1 = sk1 * 3. * xrr2 - xkrk2;
			part = sk1 * xk - xr;
			part2 = sk1 * dotk2 * xrr2 - part;
			factor = 1. - xrr2 * 3. * part2 + sk1 * 3. * dotkr2;
			factork = -sk1 + dotk2 * xkrk2 * xkrk2 + xkrk2 * part 
				- dotkrk2;
			dtxdxk1 = dtdxk1 * termx + term * factor;
			dtxkdxk1 = dtdxk1 * termxk + term * factork;
			factor = yrr2 * -3. * part2;
			factork = dotk2 * ykrk2 * xkrk2 + ykrk2 * part;
			dtydxk1 = dtdxk1 * termy + term * factor;
			dtykdxk1 = dtdxk1 * termyk + term * factork;
			factor = zrr2 * -3. * part2;
			factork = dotk2 * zkrk2 * xkrk2 + zkrk2 * part;
			dtzdxk1 = dtdxk1 * termz + term * factor;
			dtzkdxk1 = dtdxk1 * termzk + term * factork;
			dtdyk1 = sk1 * 3. * yrr2 - ykrk2;
			part = sk1 * yk - yr;
			part2 = sk1 * dotk2 * yrr2 - part;
			factor = xrr2 * -3. * part2;
			factork = dotk2 * xkrk2 * ykrk2 + xkrk2 * part;
			dtxdyk1 = dtdyk1 * termx + term * factor;
			dtxkdyk1 = dtdyk1 * termxk + term * factork;
			factor = 1. - yrr2 * 3. * part2 + sk1 * 3. * dotkr2;
			factork = -sk1 + dotk2 * ykrk2 * ykrk2 + ykrk2 * part 
				- dotkrk2;
			dtydyk1 = dtdyk1 * termy + term * factor;
			dtykdyk1 = dtdyk1 * termyk + term * factork;
			factor = zrr2 * -3. * part2;
			factork = dotk2 * zkrk2 * ykrk2 + zkrk2 * part;
			dtzdyk1 = dtdyk1 * termz + term * factor;
			dtzkdyk1 = dtdyk1 * termzk + term * factork;
			dtdzk1 = sk1 * 3. * zrr2 - zkrk2;
			part = sk1 * zk - zr;
			part2 = sk1 * dotk2 * zrr2 - part;
			factor = xrr2 * -3. * part2;
			factork = dotk2 * xkrk2 * zkrk2 + xkrk2 * part;
			dtxdzk1 = dtdzk1 * termx + term * factor;
			dtxkdzk1 = dtdzk1 * termxk + term * factork;
			factor = yrr2 * -3. * part2;
			factork = dotk2 * ykrk2 * zkrk2 + ykrk2 * part;
			dtydzk1 = dtdzk1 * termy + term * factor;
			dtykdzk1 = dtdzk1 * termyk + term * factork;
			factor = 1. - zrr2 * 3. * part2 + sk1 * 3. * dotkr2;
			factork = -sk1 + dotk2 * zkrk2 * zkrk2 + zkrk2 * part 
				- dotkrk2;
			dtzdzk1 = dtdzk1 * termz + term * factor;
			dtzkdzk1 = dtdzk1 * termzk + term * factork;
		    } else if (k2 == *i__) {
			dtdxk2 = sk2 * 3. * xrr2 + xkrk2;
			part = sk2 * xk + xr;
			part2 = sk2 * dotk2 * xrr2 - part;
			factor = -1. - xrr2 * 3. * part2 + sk2 * 3. * dotkr2;
			factork = -sk2 - dotk2 * xkrk2 * xkrk2 + xkrk2 * part 
				+ dotkrk2;
			dtxdxk2 = dtdxk2 * termx + term * factor;
			dtxkdxk2 = dtdxk2 * termxk + term * factork;
			factor = yrr2 * -3. * part2;
			factork = -dotk2 * ykrk2 * xkrk2 + ykrk2 * part;
			dtydxk2 = dtdxk2 * termy + term * factor;
			dtykdxk2 = dtdxk2 * termyk + term * factork;
			factor = zrr2 * -3. * part2;
			factork = -dotk2 * zkrk2 * xkrk2 + zkrk2 * part;
			dtzdxk2 = dtdxk2 * termz + term * factor;
			dtzkdxk2 = dtdxk2 * termzk + term * factork;
			dtdyk2 = sk2 * 3. * yrr2 + ykrk2;
			part = sk2 * yk + yr;
			part2 = sk2 * dotk2 * yrr2 - part;
			factor = xrr2 * -3. * part2;
			factork = -dotk2 * xkrk2 * ykrk2 + xkrk2 * part;
			dtxdyk2 = dtdyk2 * termx + term * factor;
			dtxkdyk2 = dtdyk2 * termxk + term * factork;
			factor = -1. - yrr2 * 3. * part2 + sk2 * 3. * dotkr2;
			factork = -sk2 - dotk2 * ykrk2 * ykrk2 + ykrk2 * part 
				+ dotkrk2;
			dtydyk2 = dtdyk2 * termy + term * factor;
			dtykdyk2 = dtdyk2 * termyk + term * factork;
			factor = zrr2 * -3. * part2;
			factork = -dotk2 * zkrk2 * ykrk2 + zkrk2 * part;
			dtzdyk2 = dtdyk2 * termz + term * factor;
			dtzkdyk2 = dtdyk2 * termzk + term * factork;
			dtdzk2 = sk2 * 3. * zrr2 + zkrk2;
			part = sk2 * zk + zr;
			part2 = sk2 * dotk2 * zrr2 - part;
			factor = xrr2 * -3. * part2;
			factork = -dotk2 * xkrk2 * zkrk2 + xkrk2 * part;
			dtxdzk2 = dtdzk2 * termx + term * factor;
			dtxkdzk2 = dtdzk2 * termxk + term * factork;
			factor = yrr2 * -3. * part2;
			factork = -dotk2 * ykrk2 * zkrk2 + ykrk2 * part;
			dtydzk2 = dtdzk2 * termy + term * factor;
			dtykdzk2 = dtdzk2 * termyk + term * factork;
			factor = -1. - zrr2 * 3. * part2 + sk2 * 3. * dotkr2;
			factork = -sk2 - dotk2 * zkrk2 * zkrk2 + zkrk2 * part 
				+ dotkrk2;
			dtzdzk2 = dtdzk2 * termz + term * factor;
			dtzkdzk2 = dtdzk2 * termzk + term * factork;
		    }

/*     increment diagonal and off-diagonal Hessian elements */

		    if (*i__ == k1) {
			hessx_ref(1, i1) = hessx_ref(1, i1) + dtxdxk1;
			hessx_ref(2, i1) = hessx_ref(2, i1) + dtydxk1;
			hessx_ref(3, i1) = hessx_ref(3, i1) + dtzdxk1;
			hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * dtxdxk1 + 
				dtxkdxk1;
			hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * dtydxk1 + 
				dtykdxk1;
			hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * dtzdxk1 + 
				dtzkdxk1;
			hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * dtxdxk1 - 
				dtxkdxk1;
			hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * dtydxk1 - 
				dtykdxk1;
			hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * dtzdxk1 - 
				dtzkdxk1;
			hessy_ref(1, i1) = hessy_ref(1, i1) + dtxdyk1;
			hessy_ref(2, i1) = hessy_ref(2, i1) + dtydyk1;
			hessy_ref(3, i1) = hessy_ref(3, i1) + dtzdyk1;
			hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * dtxdyk1 + 
				dtxkdyk1;
			hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * dtydyk1 + 
				dtykdyk1;
			hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * dtzdyk1 + 
				dtzkdyk1;
			hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * dtxdyk1 - 
				dtxkdyk1;
			hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * dtydyk1 - 
				dtykdyk1;
			hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * dtzdyk1 - 
				dtzkdyk1;
			hessz_ref(1, i1) = hessz_ref(1, i1) + dtxdzk1;
			hessz_ref(2, i1) = hessz_ref(2, i1) + dtydzk1;
			hessz_ref(3, i1) = hessz_ref(3, i1) + dtzdzk1;
			hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * dtxdzk1 + 
				dtxkdzk1;
			hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * dtydzk1 + 
				dtykdzk1;
			hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * dtzdzk1 + 
				dtzkdzk1;
			hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * dtxdzk1 - 
				dtxkdzk1;
			hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * dtydzk1 - 
				dtykdzk1;
			hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * dtzdzk1 - 
				dtzkdzk1;
		    } else if (*i__ == k2) {
			hessx_ref(1, i1) = hessx_ref(1, i1) + dtxdxk2;
			hessx_ref(2, i1) = hessx_ref(2, i1) + dtydxk2;
			hessx_ref(3, i1) = hessx_ref(3, i1) + dtzdxk2;
			hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * dtxdxk2 + 
				dtxkdxk2;
			hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * dtydxk2 + 
				dtykdxk2;
			hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * dtzdxk2 + 
				dtzkdxk2;
			hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * dtxdxk2 - 
				dtxkdxk2;
			hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * dtydxk2 - 
				dtykdxk2;
			hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * dtzdxk2 - 
				dtzkdxk2;
			hessy_ref(1, i1) = hessy_ref(1, i1) + dtxdyk2;
			hessy_ref(2, i1) = hessy_ref(2, i1) + dtydyk2;
			hessy_ref(3, i1) = hessy_ref(3, i1) + dtzdyk2;
			hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * dtxdyk2 + 
				dtxkdyk2;
			hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * dtydyk2 + 
				dtykdyk2;
			hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * dtzdyk2 + 
				dtzkdyk2;
			hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * dtxdyk2 - 
				dtxkdyk2;
			hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * dtydyk2 - 
				dtykdyk2;
			hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * dtzdyk2 - 
				dtzkdyk2;
			hessz_ref(1, i1) = hessz_ref(1, i1) + dtxdzk2;
			hessz_ref(2, i1) = hessz_ref(2, i1) + dtydzk2;
			hessz_ref(3, i1) = hessz_ref(3, i1) + dtzdzk2;
			hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * dtxdzk2 + 
				dtxkdzk2;
			hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * dtydzk2 + 
				dtykdzk2;
			hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * dtzdzk2 + 
				dtzkdzk2;
			hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * dtxdzk2 - 
				dtxkdzk2;
			hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * dtydzk2 - 
				dtykdzk2;
			hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * dtzdzk2 - 
				dtzkdzk2;
		    }

/*     more energy switching if near the cutoff distance */

		    if (r2 > shunt_1.cut2 && *i__ == k1) {
			hessx_ref(1, i1) = hessx_ref(1, i1) - sk1 * dtaperx * 
				dedxi1 + dtaperx * dedxk1 - sk1 * d2taperxx;
			hessx_ref(2, i1) = hessx_ref(2, i1) - sk1 * dtaperx * 
				dedyi1 + dtapery * dedxk1 - sk1 * d2taperxy;
			hessx_ref(3, i1) = hessx_ref(3, i1) - sk1 * dtaperx * 
				dedzi1 + dtaperz * dedxk1 - sk1 * d2taperxz;
			hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * dtaperx * 
				dedxk1 - sk1 * dtaperx * dedxk1 + sk1 * sk1 * 
				d2taperxx;
			hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * dtaperx * 
				dedyk1 - sk1 * dtapery * dedxk1 + sk1 * sk1 * 
				d2taperxy;
			hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * dtaperx * 
				dedzk1 - sk1 * dtaperz * dedxk1 + sk1 * sk1 * 
				d2taperxz;
			hessx_ref(1, k2) = hessx_ref(1, k2) - sk1 * dtaperx * 
				dedxk2 - sk2 * dtaperx * dedxk1 + sk1 * sk2 * 
				d2taperxx;
			hessx_ref(2, k2) = hessx_ref(2, k2) - sk1 * dtaperx * 
				dedyk2 - sk2 * dtapery * dedxk1 + sk1 * sk2 * 
				d2taperxy;
			hessx_ref(3, k2) = hessx_ref(3, k2) - sk1 * dtaperx * 
				dedzk2 - sk2 * dtaperz * dedxk1 + sk1 * sk2 * 
				d2taperxz;
			hessy_ref(1, i1) = hessy_ref(1, i1) - sk1 * dtapery * 
				dedxi1 + dtaperx * dedyk1 - sk1 * d2taperxy;
			hessy_ref(2, i1) = hessy_ref(2, i1) - sk1 * dtapery * 
				dedyi1 + dtapery * dedyk1 - sk1 * d2taperyy;
			hessy_ref(3, i1) = hessy_ref(3, i1) - sk1 * dtapery * 
				dedzi1 + dtaperz * dedyk1 - sk1 * d2taperyz;
			hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * dtapery * 
				dedxk1 - sk1 * dtaperx * dedyk1 + sk1 * sk1 * 
				d2taperxy;
			hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * dtapery * 
				dedyk1 - sk1 * dtapery * dedyk1 + sk1 * sk1 * 
				d2taperyy;
			hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * dtapery * 
				dedzk1 - sk1 * dtaperz * dedyk1 + sk1 * sk1 * 
				d2taperyz;
			hessy_ref(1, k2) = hessy_ref(1, k2) - sk1 * dtapery * 
				dedxk2 - sk2 * dtaperx * dedyk1 + sk1 * sk2 * 
				d2taperxy;
			hessy_ref(2, k2) = hessy_ref(2, k2) - sk1 * dtapery * 
				dedyk2 - sk2 * dtapery * dedyk1 + sk1 * sk2 * 
				d2taperyy;
			hessy_ref(3, k2) = hessy_ref(3, k2) - sk1 * dtapery * 
				dedzk2 - sk2 * dtaperz * dedyk1 + sk1 * sk2 * 
				d2taperyz;
			hessz_ref(1, i1) = hessz_ref(1, i1) - sk1 * dtaperz * 
				dedxi1 + dtaperx * dedzk1 - sk1 * d2taperxz;
			hessz_ref(2, i1) = hessz_ref(2, i1) - sk1 * dtaperz * 
				dedyi1 + dtapery * dedzk1 - sk1 * d2taperyz;
			hessz_ref(3, i1) = hessz_ref(3, i1) - sk1 * dtaperz * 
				dedzi1 + dtaperz * dedzk1 - sk1 * d2taperzz;
			hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * dtaperz * 
				dedxk1 - sk1 * dtaperx * dedzk1 + sk1 * sk1 * 
				d2taperxz;
			hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * dtaperz * 
				dedyk1 - sk1 * dtapery * dedzk1 + sk1 * sk1 * 
				d2taperyz;
			hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * dtaperz * 
				dedzk1 - sk1 * dtaperz * dedzk1 + sk1 * sk1 * 
				d2taperzz;
			hessz_ref(1, k2) = hessz_ref(1, k2) - sk1 * dtaperz * 
				dedxk2 - sk2 * dtaperx * dedzk1 + sk1 * sk2 * 
				d2taperxz;
			hessz_ref(2, k2) = hessz_ref(2, k2) - sk1 * dtaperz * 
				dedyk2 - sk2 * dtapery * dedzk1 + sk1 * sk2 * 
				d2taperyz;
			hessz_ref(3, k2) = hessz_ref(3, k2) - sk1 * dtaperz * 
				dedzk2 - sk2 * dtaperz * dedzk1 + sk1 * sk2 * 
				d2taperzz;
		    } else if (r2 > shunt_1.cut2 && *i__ == k2) {
			hessx_ref(1, i1) = hessx_ref(1, i1) - sk2 * dtaperx * 
				dedxi1 + dtaperx * dedxk2 - sk2 * d2taperxx;
			hessx_ref(2, i1) = hessx_ref(2, i1) - sk2 * dtaperx * 
				dedyi1 + dtapery * dedxk2 - sk2 * d2taperxy;
			hessx_ref(3, i1) = hessx_ref(3, i1) - sk2 * dtaperx * 
				dedzi1 + dtaperz * dedxk2 - sk2 * d2taperxz;
			hessx_ref(1, k1) = hessx_ref(1, k1) - sk2 * dtaperx * 
				dedxk1 - sk1 * dtaperx * dedxk2 + sk1 * sk2 * 
				d2taperxx;
			hessx_ref(2, k1) = hessx_ref(2, k1) - sk2 * dtaperx * 
				dedyk1 - sk1 * dtapery * dedxk2 + sk1 * sk2 * 
				d2taperxy;
			hessx_ref(3, k1) = hessx_ref(3, k1) - sk2 * dtaperx * 
				dedzk1 - sk1 * dtaperz * dedxk2 + sk1 * sk2 * 
				d2taperxz;
			hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * dtaperx * 
				dedxk2 - sk2 * dtaperx * dedxk2 + sk2 * sk2 * 
				d2taperxx;
			hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * dtaperx * 
				dedyk2 - sk2 * dtapery * dedxk2 + sk2 * sk2 * 
				d2taperxy;
			hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * dtaperx * 
				dedzk2 - sk2 * dtaperz * dedxk2 + sk2 * sk2 * 
				d2taperxz;
			hessy_ref(1, i1) = hessy_ref(1, i1) - sk2 * dtapery * 
				dedxi1 + dtaperx * dedyk2 - sk2 * d2taperxy;
			hessy_ref(2, i1) = hessy_ref(2, i1) - sk2 * dtapery * 
				dedyi1 + dtapery * dedyk2 - sk2 * d2taperyy;
			hessy_ref(3, i1) = hessy_ref(3, i1) - sk2 * dtapery * 
				dedzi1 + dtaperz * dedyk2 - sk2 * d2taperyz;
			hessy_ref(1, k1) = hessy_ref(1, k1) - sk2 * dtapery * 
				dedxk1 - sk1 * dtaperx * dedyk2 + sk1 * sk2 * 
				d2taperxy;
			hessy_ref(2, k1) = hessy_ref(2, k1) - sk2 * dtapery * 
				dedyk1 - sk1 * dtapery * dedyk2 + sk1 * sk2 * 
				d2taperyy;
			hessy_ref(3, k1) = hessy_ref(3, k1) - sk2 * dtapery * 
				dedzk1 - sk1 * dtaperz * dedyk2 + sk1 * sk2 * 
				d2taperyz;
			hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * dtapery * 
				dedxk2 - sk2 * dtaperx * dedyk2 + sk2 * sk2 * 
				d2taperxy;
			hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * dtapery * 
				dedyk2 - sk2 * dtapery * dedyk2 + sk2 * sk2 * 
				d2taperyy;
			hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * dtapery * 
				dedzk2 - sk2 * dtaperz * dedyk2 + sk2 * sk2 * 
				d2taperyz;
			hessz_ref(1, i1) = hessz_ref(1, i1) - sk2 * dtaperz * 
				dedxi1 + dtaperx * dedzk2 - sk2 * d2taperxz;
			hessz_ref(2, i1) = hessz_ref(2, i1) - sk2 * dtaperz * 
				dedyi1 + dtapery * dedzk2 - sk2 * d2taperyz;
			hessz_ref(3, i1) = hessz_ref(3, i1) - sk2 * dtaperz * 
				dedzi1 + dtaperz * dedzk2 - sk2 * d2taperzz;
			hessz_ref(1, k1) = hessz_ref(1, k1) - sk2 * dtaperz * 
				dedxk1 - sk1 * dtaperx * dedzk2 + sk1 * sk2 * 
				d2taperxz;
			hessz_ref(2, k1) = hessz_ref(2, k1) - sk2 * dtaperz * 
				dedyk1 - sk1 * dtapery * dedzk2 + sk1 * sk2 * 
				d2taperyz;
			hessz_ref(3, k1) = hessz_ref(3, k1) - sk2 * dtaperz * 
				dedzk1 - sk1 * dtaperz * dedzk2 + sk1 * sk2 * 
				d2taperzz;
			hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * dtaperz * 
				dedxk2 - sk2 * dtaperx * dedzk2 + sk2 * sk2 * 
				d2taperxz;
			hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * dtaperz * 
				dedyk2 - sk2 * dtapery * dedzk2 + sk2 * sk2 * 
				d2taperyz;
			hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * dtaperz * 
				dedzk2 - sk2 * dtaperz * dedzk2 + sk2 * sk2 * 
				d2taperzz;
		    }
		}
	    }
	}
L20:
	;
    }

/*     for periodic boundary conditions with large cutoffs */
/*     neighbors must be found by the replicates method */

    if (! bound_1.use_replica__) {
	return 0;
    }

/*     calculate interaction energy with other unit cells */

    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
	i1 = charge_1.iion[ii - 1];
	if (i1 != *i__) {
	    goto L30;
	}
	skip[i1 - 1] = i1;
	i__2 = couple_1.n12[i1 - 1];
	for (k = 1; k <= i__2; ++k) {
	    skip[i12_ref(k, i1) - 1] = i1;
	}
	xi = atoms_1.x[i1 - 1];
	yi = atoms_1.y[i1 - 1];
	zi = atoms_1.z__[i1 - 1];
	fi = f * charge_1.pchg[ii - 1];

/*     decide whether to compute the current interaction */

	i__2 = dipole_1.ndipole;
	for (k = 1; k <= i__2; ++k) {
	    k1 = idpl_ref(1, k);
	    k2 = idpl_ref(2, k);
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i1, &k1, &k2, &c__0, &c__0, &c__0);
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		sk1 = 1. - dipole_1.sdpl[k - 1];
		sk2 = dipole_1.sdpl[k - 1];
		i__3 = cell_1.ncell;
		for (jcell = 1; jcell <= i__3; ++jcell) {
		    xk = atoms_1.x[k1 - 1] - atoms_1.x[k2 - 1];
		    yk = atoms_1.y[k1 - 1] - atoms_1.y[k2 - 1];
		    zk = atoms_1.z__[k1 - 1] - atoms_1.z__[k2 - 1];
		    xr = xi - atoms_1.x[k1 - 1] + xk * sk2;
		    yr = yi - atoms_1.y[k1 - 1] + yk * sk2;
		    zr = zi - atoms_1.z__[k1 - 1] + zk * sk2;
		    imager_(&xr, &yr, &zr, &jcell);
		    r2 = xr * xr + yr * yr + zr * zr;
		    if (r2 <= shunt_1.off2) {
			rk2 = xk * xk + yk * yk + zk * zk;
			rkr3 = sqrt(rk2 * r2) * r2;
			dotk = xk * xr + yk * yr + zk * zr;
			fik = -fi * dipole_1.bdpl[k - 1];
			if (bound_1.use_polymer__) {
			    if (r2 < bound_1.polycut2) {
				if (skip[k1 - 1] == i1 || skip[k2 - 1] != i1) 
					{
				    fik = 0.;
				}
			    }
			}

/*     scale the interaction based on its group membership */

			if (group_1.use_group__) {
			    fik *= fgrp;
			}

/*     some abbreviations used in various chain rule terms */

			xrr2 = xr / r2;
			yrr2 = yr / r2;
			zrr2 = zr / r2;
			xkrk2 = xk / rk2;
			ykrk2 = yk / rk2;
			zkrk2 = zk / rk2;
			dotk2 = dotk * 2.;
			dotkr2 = dotk / r2;

/*     form the chain rule terms for first derivatives */

			term = fik / rkr3;
			term2 = dotk * -3.;
			termx = term * (xk + xrr2 * term2);
			termy = term * (yk + yrr2 * term2);
			termz = term * (zk + zrr2 * term2);
			termxk = term * (xr - dotk * xkrk2);
			termyk = term * (yr - dotk * ykrk2);
			termzk = term * (zr - dotk * zkrk2);

/*     use energy switching if near the cutoff distance */

			if (r2 > shunt_1.cut2) {
			    e = fik * dotk / rkr3;
			    dedxi1 = termx;
			    dedyi1 = termy;
			    dedzi1 = termz;
			    dedxk1 = -sk1 * termx + termxk;
			    dedyk1 = -sk1 * termy + termyk;
			    dedzk1 = -sk1 * termz + termzk;
			    dedxk2 = -sk2 * termx - termxk;
			    dedyk2 = -sk2 * termy - termyk;
			    dedzk2 = -sk2 * termz - termzk;
			    r__ = sqrt(r2);
			    r3 = r2 * r__;
			    r4 = r2 * r2;
			    r5 = r2 * r3;
			    taper = shunt_1.c5 * r5 + shunt_1.c4 * r4 + 
				    shunt_1.c3 * r3 + shunt_1.c2 * r2 + 
				    shunt_1.c1 * r__ + shunt_1.c0;
			    dtaper = shunt_1.c5 * 5. * r4 + shunt_1.c4 * 4. * 
				    r3 + shunt_1.c3 * 3. * r2 + shunt_1.c2 * 
				    2. * r__ + shunt_1.c1;
			    d2taper = shunt_1.c5 * 20. * r3 + shunt_1.c4 * 
				    12. * r2 + shunt_1.c3 * 6. * r__ + 
				    shunt_1.c2 * 2.;
			    dtaper /= r__;
			    dtaperx = xr * dtaper;
			    dtapery = yr * dtaper;
			    dtaperz = zr * dtaper;
			    d2taper = e * (d2taper - dtaper);
			    dtaper = e * dtaper;
			    d2taperxx = xr * xrr2 * d2taper + dtaper;
			    d2taperxy = xr * yrr2 * d2taper;
			    d2taperxz = xr * zrr2 * d2taper;
			    d2taperyy = yr * yrr2 * d2taper + dtaper;
			    d2taperyz = yr * zrr2 * d2taper;
			    d2taperzz = zr * zrr2 * d2taper + dtaper;
			    term *= taper;
			    termx *= taper;
			    termy *= taper;
			    termz *= taper;
			    termxk *= taper;
			    termyk *= taper;
			    termzk *= taper;
			}

/*     chain rule terms for second derivative components */

			dtdxi1 = xrr2 * -3.;
			part = xk - dotk2 * xrr2;
			factor = (dotkr2 + xrr2 * part) * -3.;
			factork = 1. - xk * xkrk2;
			dtxdxi1 = dtdxi1 * termx + term * factor;
			dtxkdxi1 = dtdxi1 * termxk + term * factork;
			factor = yrr2 * -3. * part;
			factork = -yk * xkrk2;
			dtydxi1 = dtdxi1 * termy + term * factor;
			dtykdxi1 = dtdxi1 * termyk + term * factork;
			factor = zrr2 * -3. * part;
			factork = -zk * xkrk2;
			dtzdxi1 = dtdxi1 * termz + term * factor;
			dtzkdxi1 = dtdxi1 * termzk + term * factork;
			dtdyi1 = yrr2 * -3.;
			part = yk - dotk2 * yrr2;
			factor = xrr2 * -3. * part;
			factork = -xk * ykrk2;
			dtxdyi1 = dtdyi1 * termx + term * factor;
			dtxkdyi1 = dtdyi1 * termxk + term * factork;
			factor = (dotkr2 + yrr2 * part) * -3.;
			factork = 1. - yk * ykrk2;
			dtydyi1 = dtdyi1 * termy + term * factor;
			dtykdyi1 = dtdyi1 * termyk + term * factork;
			factor = zrr2 * -3. * part;
			factork = -zk * ykrk2;
			dtzdyi1 = dtdyi1 * termz + term * factor;
			dtzkdyi1 = dtdyi1 * termzk + term * factork;
			dtdzi1 = zrr2 * -3.;
			part = zk - dotk2 * zrr2;
			factor = xrr2 * -3. * part;
			factork = -xk * zkrk2;
			dtxdzi1 = dtdzi1 * termx + term * factor;
			dtxkdzi1 = dtdzi1 * termxk + term * factork;
			factor = yrr2 * -3. * part;
			factork = -yk * zkrk2;
			dtydzi1 = dtdzi1 * termy + term * factor;
			dtykdzi1 = dtdzi1 * termyk + term * factork;
			factor = (dotkr2 + zrr2 * part) * -3.;
			factork = 1. - zk * zkrk2;
			dtzdzi1 = dtdzi1 * termz + term * factor;
			dtzkdzi1 = dtdzi1 * termzk + term * factork;

/*     increment diagonal and off-diagonal Hessian elements */

			hessx_ref(1, i1) = hessx_ref(1, i1) + dtxdxi1;
			hessx_ref(2, i1) = hessx_ref(2, i1) + dtydxi1;
			hessx_ref(3, i1) = hessx_ref(3, i1) + dtzdxi1;
			hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * dtxdxi1 + 
				dtxkdxi1;
			hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * dtydxi1 + 
				dtykdxi1;
			hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * dtzdxi1 + 
				dtzkdxi1;
			hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * dtxdxi1 - 
				dtxkdxi1;
			hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * dtydxi1 - 
				dtykdxi1;
			hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * dtzdxi1 - 
				dtzkdxi1;
			hessy_ref(1, i1) = hessy_ref(1, i1) + dtxdyi1;
			hessy_ref(2, i1) = hessy_ref(2, i1) + dtydyi1;
			hessy_ref(3, i1) = hessy_ref(3, i1) + dtzdyi1;
			hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * dtxdyi1 + 
				dtxkdyi1;
			hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * dtydyi1 + 
				dtykdyi1;
			hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * dtzdyi1 + 
				dtzkdyi1;
			hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * dtxdyi1 - 
				dtxkdyi1;
			hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * dtydyi1 - 
				dtykdyi1;
			hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * dtzdyi1 - 
				dtzkdyi1;
			hessz_ref(1, i1) = hessz_ref(1, i1) + dtxdzi1;
			hessz_ref(2, i1) = hessz_ref(2, i1) + dtydzi1;
			hessz_ref(3, i1) = hessz_ref(3, i1) + dtzdzi1;
			hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * dtxdzi1 + 
				dtxkdzi1;
			hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * dtydzi1 + 
				dtykdzi1;
			hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * dtzdzi1 + 
				dtzkdzi1;
			hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * dtxdzi1 - 
				dtxkdzi1;
			hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * dtydzi1 - 
				dtykdzi1;
			hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * dtzdzi1 - 
				dtzkdzi1;

/*     more energy switching if near the cutoff distance */

			if (r2 > shunt_1.cut2) {
			    hessx_ref(1, i1) = hessx_ref(1, i1) + dtaperx * 
				    dedxi1 + dtaperx * dedxi1 + d2taperxx;
			    hessx_ref(2, i1) = hessx_ref(2, i1) + dtaperx * 
				    dedyi1 + dtapery * dedxi1 + d2taperxy;
			    hessx_ref(3, i1) = hessx_ref(3, i1) + dtaperx * 
				    dedzi1 + dtaperz * dedxi1 + d2taperxz;
			    hessx_ref(1, k1) = hessx_ref(1, k1) + dtaperx * 
				    dedxk1 - sk1 * (dtaperx * dedxi1 + 
				    d2taperxx);
			    hessx_ref(2, k1) = hessx_ref(2, k1) + dtaperx * 
				    dedyk1 - sk1 * (dtapery * dedxi1 + 
				    d2taperxy);
			    hessx_ref(3, k1) = hessx_ref(3, k1) + dtaperx * 
				    dedzk1 - sk1 * (dtaperz * dedxi1 + 
				    d2taperxz);
			    hessx_ref(1, k2) = hessx_ref(1, k2) + dtaperx * 
				    dedxk2 - sk2 * (dtaperx * dedxi1 + 
				    d2taperxx);
			    hessx_ref(2, k2) = hessx_ref(2, k2) + dtaperx * 
				    dedyk2 - sk2 * (dtapery * dedxi1 + 
				    d2taperxy);
			    hessx_ref(3, k2) = hessx_ref(3, k2) + dtaperx * 
				    dedzk2 - sk2 * (dtaperz * dedxi1 + 
				    d2taperxz);
			    hessy_ref(1, i1) = hessy_ref(1, i1) + dtapery * 
				    dedxi1 + dtaperx * dedyi1 + d2taperxy;
			    hessy_ref(2, i1) = hessy_ref(2, i1) + dtapery * 
				    dedyi1 + dtapery * dedyi1 + d2taperyy;
			    hessy_ref(3, i1) = hessy_ref(3, i1) + dtapery * 
				    dedzi1 + dtaperz * dedyi1 + d2taperyz;
			    hessy_ref(1, k1) = hessy_ref(1, k1) + dtapery * 
				    dedxk1 - sk1 * (dtaperx * dedyi1 + 
				    d2taperxy);
			    hessy_ref(2, k1) = hessy_ref(2, k1) + dtapery * 
				    dedyk1 - sk1 * (dtapery * dedyi1 + 
				    d2taperyy);
			    hessy_ref(3, k1) = hessy_ref(3, k1) + dtapery * 
				    dedzk1 - sk1 * (dtaperz * dedyi1 + 
				    d2taperyz);
			    hessy_ref(1, k2) = hessy_ref(1, k2) + dtapery * 
				    dedxk2 - sk2 * (dtaperx * dedyi1 + 
				    d2taperxy);
			    hessy_ref(2, k2) = hessy_ref(2, k2) + dtapery * 
				    dedyk2 - sk2 * (dtapery * dedyi1 + 
				    d2taperyy);
			    hessy_ref(3, k2) = hessy_ref(3, k2) + dtapery * 
				    dedzk2 - sk2 * (dtaperz * dedyi1 + 
				    d2taperyz);
			    hessz_ref(1, i1) = hessz_ref(1, i1) + dtaperz * 
				    dedxi1 + dtaperx * dedzi1 + d2taperxz;
			    hessz_ref(2, i1) = hessz_ref(2, i1) + dtaperz * 
				    dedyi1 + dtapery * dedzi1 + d2taperyz;
			    hessz_ref(3, i1) = hessz_ref(3, i1) + dtaperz * 
				    dedzi1 + dtaperz * dedzi1 + d2taperzz;
			    hessz_ref(1, k1) = hessz_ref(1, k1) + dtaperz * 
				    dedxk1 - sk1 * (dtaperx * dedzi1 + 
				    d2taperxz);
			    hessz_ref(2, k1) = hessz_ref(2, k1) + dtaperz * 
				    dedyk1 - sk1 * (dtapery * dedzi1 + 
				    d2taperyz);
			    hessz_ref(3, k1) = hessz_ref(3, k1) + dtaperz * 
				    dedzk1 - sk1 * (dtaperz * dedzi1 + 
				    d2taperzz);
			    hessz_ref(1, k2) = hessz_ref(1, k2) + dtaperz * 
				    dedxk2 - sk2 * (dtaperx * dedzi1 + 
				    d2taperxz);
			    hessz_ref(2, k2) = hessz_ref(2, k2) + dtaperz * 
				    dedyk2 - sk2 * (dtapery * dedzi1 + 
				    d2taperyz);
			    hessz_ref(3, k2) = hessz_ref(3, k2) + dtaperz * 
				    dedzk2 - sk2 * (dtaperz * dedzi1 + 
				    d2taperzz);
			}
		    }
		}
	    }
	}
L30:
	;
    }

/*     see if the atom of interest is part of a dipole */

    i__1 = dipole_1.ndipole;
    for (k = 1; k <= i__1; ++k) {
	k1 = idpl_ref(1, k);
	k2 = idpl_ref(2, k);
	if (k1 != *i__ && k2 != *i__) {
	    goto L40;
	}
	i__2 = couple_1.n12[k1 - 1];
	for (ii = 1; ii <= i__2; ++ii) {
	    omit[i12_ref(ii, k1) - 1] = k;
	}
	i__2 = couple_1.n12[k2 - 1];
	for (ii = 1; ii <= i__2; ++ii) {
	    omit[i12_ref(ii, k2) - 1] = k;
	}
	sk1 = 1. - dipole_1.sdpl[k - 1];
	sk2 = dipole_1.sdpl[k - 1];
	xk = atoms_1.x[k1 - 1] - atoms_1.x[k2 - 1];
	yk = atoms_1.y[k1 - 1] - atoms_1.y[k2 - 1];
	zk = atoms_1.z__[k1 - 1] - atoms_1.z__[k2 - 1];
	rk2 = xk * xk + yk * yk + zk * zk;
	xq = atoms_1.x[k1 - 1] - xk * sk2;
	yq = atoms_1.y[k1 - 1] - yk * sk2;
	zq = atoms_1.z__[k1 - 1] - zk * sk2;
	fk = -f * dipole_1.bdpl[k - 1];

/*     decide whether to compute the current interaction */

	i__2 = charge_1.nion;
	for (ii = 1; ii <= i__2; ++ii) {
	    i1 = charge_1.iion[ii - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i1, &k1, &k2, &c__0, &c__0, &c__0);
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		i__3 = cell_1.ncell;
		for (jcell = 1; jcell <= i__3; ++jcell) {
		    xr = atoms_1.x[i1 - 1] - xq;
		    yr = atoms_1.y[i1 - 1] - yq;
		    zr = atoms_1.z__[i1 - 1] - zq;
		    imager_(&xr, &yr, &zr, &jcell);
		    r2 = xr * xr + yr * yr + zr * zr;
		    if (r2 <= shunt_1.off2) {
			rkr3 = sqrt(rk2 * r2) * r2;
			dotk = xk * xr + yk * yr + zk * zr;
			fik = fk * charge_1.pchg[ii - 1];
			if (bound_1.use_polymer__) {
			    if (r2 < bound_1.polycut2) {
				if (omit[i1 - 1] != k) {
				    fik = 0.;
				}
			    }
			}

/*     scale the interaction based on its group membership */

			if (group_1.use_group__) {
			    fik *= fgrp;
			}

/*     some abbreviations used in various chain rule terms */

			xrr2 = xr / r2;
			yrr2 = yr / r2;
			zrr2 = zr / r2;
			xkrk2 = xk / rk2;
			ykrk2 = yk / rk2;
			zkrk2 = zk / rk2;
			dotk2 = dotk * 2.;
			dotkr2 = dotk / r2;
			dotkrk2 = dotk / rk2;

/*     form the chain rule terms for first derivatives */

			term = fik / rkr3;
			term2 = dotk * -3.;
			termx = term * (xk + xrr2 * term2);
			termy = term * (yk + yrr2 * term2);
			termz = term * (zk + zrr2 * term2);
			termxk = term * (xr - dotk * xkrk2);
			termyk = term * (yr - dotk * ykrk2);
			termzk = term * (zr - dotk * zkrk2);

/*     use energy switching if near the cutoff distance */

			if (r2 > shunt_1.cut2) {
			    e = fik * dotk / rkr3;
			    dedxi1 = termx;
			    dedyi1 = termy;
			    dedzi1 = termz;
			    dedxk1 = -sk1 * termx + termxk;
			    dedyk1 = -sk1 * termy + termyk;
			    dedzk1 = -sk1 * termz + termzk;
			    dedxk2 = -sk2 * termx - termxk;
			    dedyk2 = -sk2 * termy - termyk;
			    dedzk2 = -sk2 * termz - termzk;
			    r__ = sqrt(r2);
			    r3 = r2 * r__;
			    r4 = r2 * r2;
			    r5 = r2 * r3;
			    taper = shunt_1.c5 * r5 + shunt_1.c4 * r4 + 
				    shunt_1.c3 * r3 + shunt_1.c2 * r2 + 
				    shunt_1.c1 * r__ + shunt_1.c0;
			    dtaper = shunt_1.c5 * 5. * r4 + shunt_1.c4 * 4. * 
				    r3 + shunt_1.c3 * 3. * r2 + shunt_1.c2 * 
				    2. * r__ + shunt_1.c1;
			    d2taper = shunt_1.c5 * 20. * r3 + shunt_1.c4 * 
				    12. * r2 + shunt_1.c3 * 6. * r__ + 
				    shunt_1.c2 * 2.;
			    dtaper /= r__;
			    dtaperx = xr * dtaper;
			    dtapery = yr * dtaper;
			    dtaperz = zr * dtaper;
			    d2taper = e * (d2taper - dtaper);
			    dtaper = e * dtaper;
			    d2taperxx = xr * xrr2 * d2taper + dtaper;
			    d2taperxy = xr * yrr2 * d2taper;
			    d2taperxz = xr * zrr2 * d2taper;
			    d2taperyy = yr * yrr2 * d2taper + dtaper;
			    d2taperyz = yr * zrr2 * d2taper;
			    d2taperzz = zr * zrr2 * d2taper + dtaper;
			    term *= taper;
			    termx *= taper;
			    termy *= taper;
			    termz *= taper;
			    termxk *= taper;
			    termyk *= taper;
			    termzk *= taper;
			}

/*     chain rule terms for second derivative components */

			if (k1 == *i__) {
			    dtdxk1 = sk1 * 3. * xrr2 - xkrk2;
			    part = sk1 * xk - xr;
			    part2 = sk1 * dotk2 * xrr2 - part;
			    factor = 1. - xrr2 * 3. * part2 + sk1 * 3. * 
				    dotkr2;
			    factork = -sk1 + dotk2 * xkrk2 * xkrk2 + xkrk2 * 
				    part - dotkrk2;
			    dtxdxk1 = dtdxk1 * termx + term * factor;
			    dtxkdxk1 = dtdxk1 * termxk + term * factork;
			    factor = yrr2 * -3. * part2;
			    factork = dotk2 * ykrk2 * xkrk2 + ykrk2 * part;
			    dtydxk1 = dtdxk1 * termy + term * factor;
			    dtykdxk1 = dtdxk1 * termyk + term * factork;
			    factor = zrr2 * -3. * part2;
			    factork = dotk2 * zkrk2 * xkrk2 + zkrk2 * part;
			    dtzdxk1 = dtdxk1 * termz + term * factor;
			    dtzkdxk1 = dtdxk1 * termzk + term * factork;
			    dtdyk1 = sk1 * 3. * yrr2 - ykrk2;
			    part = sk1 * yk - yr;
			    part2 = sk1 * dotk2 * yrr2 - part;
			    factor = xrr2 * -3. * part2;
			    factork = dotk2 * xkrk2 * ykrk2 + xkrk2 * part;
			    dtxdyk1 = dtdyk1 * termx + term * factor;
			    dtxkdyk1 = dtdyk1 * termxk + term * factork;
			    factor = 1. - yrr2 * 3. * part2 + sk1 * 3. * 
				    dotkr2;
			    factork = -sk1 + dotk2 * ykrk2 * ykrk2 + ykrk2 * 
				    part - dotkrk2;
			    dtydyk1 = dtdyk1 * termy + term * factor;
			    dtykdyk1 = dtdyk1 * termyk + term * factork;
			    factor = zrr2 * -3. * part2;
			    factork = dotk2 * zkrk2 * ykrk2 + zkrk2 * part;
			    dtzdyk1 = dtdyk1 * termz + term * factor;
			    dtzkdyk1 = dtdyk1 * termzk + term * factork;
			    dtdzk1 = sk1 * 3. * zrr2 - zkrk2;
			    part = sk1 * zk - zr;
			    part2 = sk1 * dotk2 * zrr2 - part;
			    factor = xrr2 * -3. * part2;
			    factork = dotk2 * xkrk2 * zkrk2 + xkrk2 * part;
			    dtxdzk1 = dtdzk1 * termx + term * factor;
			    dtxkdzk1 = dtdzk1 * termxk + term * factork;
			    factor = yrr2 * -3. * part2;
			    factork = dotk2 * ykrk2 * zkrk2 + ykrk2 * part;
			    dtydzk1 = dtdzk1 * termy + term * factor;
			    dtykdzk1 = dtdzk1 * termyk + term * factork;
			    factor = 1. - zrr2 * 3. * part2 + sk1 * 3. * 
				    dotkr2;
			    factork = -sk1 + dotk2 * zkrk2 * zkrk2 + zkrk2 * 
				    part - dotkrk2;
			    dtzdzk1 = dtdzk1 * termz + term * factor;
			    dtzkdzk1 = dtdzk1 * termzk + term * factork;
			} else if (k2 == *i__) {
			    dtdxk2 = sk2 * 3. * xrr2 + xkrk2;
			    part = sk2 * xk + xr;
			    part2 = sk2 * dotk2 * xrr2 - part;
			    factor = -1. - xrr2 * 3. * part2 + sk2 * 3. * 
				    dotkr2;
			    factork = -sk2 - dotk2 * xkrk2 * xkrk2 + xkrk2 * 
				    part + dotkrk2;
			    dtxdxk2 = dtdxk2 * termx + term * factor;
			    dtxkdxk2 = dtdxk2 * termxk + term * factork;
			    factor = yrr2 * -3. * part2;
			    factork = -dotk2 * ykrk2 * xkrk2 + ykrk2 * part;
			    dtydxk2 = dtdxk2 * termy + term * factor;
			    dtykdxk2 = dtdxk2 * termyk + term * factork;
			    factor = zrr2 * -3. * part2;
			    factork = -dotk2 * zkrk2 * xkrk2 + zkrk2 * part;
			    dtzdxk2 = dtdxk2 * termz + term * factor;
			    dtzkdxk2 = dtdxk2 * termzk + term * factork;
			    dtdyk2 = sk2 * 3. * yrr2 + ykrk2;
			    part = sk2 * yk + yr;
			    part2 = sk2 * dotk2 * yrr2 - part;
			    factor = xrr2 * -3. * part2;
			    factork = -dotk2 * xkrk2 * ykrk2 + xkrk2 * part;
			    dtxdyk2 = dtdyk2 * termx + term * factor;
			    dtxkdyk2 = dtdyk2 * termxk + term * factork;
			    factor = -1. - yrr2 * 3. * part2 + sk2 * 3. * 
				    dotkr2;
			    factork = -sk2 - dotk2 * ykrk2 * ykrk2 + ykrk2 * 
				    part + dotkrk2;
			    dtydyk2 = dtdyk2 * termy + term * factor;
			    dtykdyk2 = dtdyk2 * termyk + term * factork;
			    factor = zrr2 * -3. * part2;
			    factork = -dotk2 * zkrk2 * ykrk2 + zkrk2 * part;
			    dtzdyk2 = dtdyk2 * termz + term * factor;
			    dtzkdyk2 = dtdyk2 * termzk + term * factork;
			    dtdzk2 = sk2 * 3. * zrr2 + zkrk2;
			    part = sk2 * zk + zr;
			    part2 = sk2 * dotk2 * zrr2 - part;
			    factor = xrr2 * -3. * part2;
			    factork = -dotk2 * xkrk2 * zkrk2 + xkrk2 * part;
			    dtxdzk2 = dtdzk2 * termx + term * factor;
			    dtxkdzk2 = dtdzk2 * termxk + term * factork;
			    factor = yrr2 * -3. * part2;
			    factork = -dotk2 * ykrk2 * zkrk2 + ykrk2 * part;
			    dtydzk2 = dtdzk2 * termy + term * factor;
			    dtykdzk2 = dtdzk2 * termyk + term * factork;
			    factor = -1. - zrr2 * 3. * part2 + sk2 * 3. * 
				    dotkr2;
			    factork = -sk2 - dotk2 * zkrk2 * zkrk2 + zkrk2 * 
				    part + dotkrk2;
			    dtzdzk2 = dtdzk2 * termz + term * factor;
			    dtzkdzk2 = dtdzk2 * termzk + term * factork;
			}

/*     increment diagonal and off-diagonal Hessian elements */

			if (*i__ == k1) {
			    hessx_ref(1, i1) = hessx_ref(1, i1) + dtxdxk1;
			    hessx_ref(2, i1) = hessx_ref(2, i1) + dtydxk1;
			    hessx_ref(3, i1) = hessx_ref(3, i1) + dtzdxk1;
			    hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * 
				    dtxdxk1 + dtxkdxk1;
			    hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * 
				    dtydxk1 + dtykdxk1;
			    hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * 
				    dtzdxk1 + dtzkdxk1;
			    hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * 
				    dtxdxk1 - dtxkdxk1;
			    hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * 
				    dtydxk1 - dtykdxk1;
			    hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * 
				    dtzdxk1 - dtzkdxk1;
			    hessy_ref(1, i1) = hessy_ref(1, i1) + dtxdyk1;
			    hessy_ref(2, i1) = hessy_ref(2, i1) + dtydyk1;
			    hessy_ref(3, i1) = hessy_ref(3, i1) + dtzdyk1;
			    hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * 
				    dtxdyk1 + dtxkdyk1;
			    hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * 
				    dtydyk1 + dtykdyk1;
			    hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * 
				    dtzdyk1 + dtzkdyk1;
			    hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * 
				    dtxdyk1 - dtxkdyk1;
			    hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * 
				    dtydyk1 - dtykdyk1;
			    hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * 
				    dtzdyk1 - dtzkdyk1;
			    hessz_ref(1, i1) = hessz_ref(1, i1) + dtxdzk1;
			    hessz_ref(2, i1) = hessz_ref(2, i1) + dtydzk1;
			    hessz_ref(3, i1) = hessz_ref(3, i1) + dtzdzk1;
			    hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * 
				    dtxdzk1 + dtxkdzk1;
			    hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * 
				    dtydzk1 + dtykdzk1;
			    hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * 
				    dtzdzk1 + dtzkdzk1;
			    hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * 
				    dtxdzk1 - dtxkdzk1;
			    hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * 
				    dtydzk1 - dtykdzk1;
			    hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * 
				    dtzdzk1 - dtzkdzk1;
			} else if (*i__ == k2) {
			    hessx_ref(1, i1) = hessx_ref(1, i1) + dtxdxk2;
			    hessx_ref(2, i1) = hessx_ref(2, i1) + dtydxk2;
			    hessx_ref(3, i1) = hessx_ref(3, i1) + dtzdxk2;
			    hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * 
				    dtxdxk2 + dtxkdxk2;
			    hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * 
				    dtydxk2 + dtykdxk2;
			    hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * 
				    dtzdxk2 + dtzkdxk2;
			    hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * 
				    dtxdxk2 - dtxkdxk2;
			    hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * 
				    dtydxk2 - dtykdxk2;
			    hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * 
				    dtzdxk2 - dtzkdxk2;
			    hessy_ref(1, i1) = hessy_ref(1, i1) + dtxdyk2;
			    hessy_ref(2, i1) = hessy_ref(2, i1) + dtydyk2;
			    hessy_ref(3, i1) = hessy_ref(3, i1) + dtzdyk2;
			    hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * 
				    dtxdyk2 + dtxkdyk2;
			    hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * 
				    dtydyk2 + dtykdyk2;
			    hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * 
				    dtzdyk2 + dtzkdyk2;
			    hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * 
				    dtxdyk2 - dtxkdyk2;
			    hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * 
				    dtydyk2 - dtykdyk2;
			    hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * 
				    dtzdyk2 - dtzkdyk2;
			    hessz_ref(1, i1) = hessz_ref(1, i1) + dtxdzk2;
			    hessz_ref(2, i1) = hessz_ref(2, i1) + dtydzk2;
			    hessz_ref(3, i1) = hessz_ref(3, i1) + dtzdzk2;
			    hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * 
				    dtxdzk2 + dtxkdzk2;
			    hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * 
				    dtydzk2 + dtykdzk2;
			    hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * 
				    dtzdzk2 + dtzkdzk2;
			    hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * 
				    dtxdzk2 - dtxkdzk2;
			    hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * 
				    dtydzk2 - dtykdzk2;
			    hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * 
				    dtzdzk2 - dtzkdzk2;
			}

/*     more energy switching if near the cutoff distance */

			if (r2 > shunt_1.cut2 && *i__ == k1) {
			    hessx_ref(1, i1) = hessx_ref(1, i1) - sk1 * 
				    dtaperx * dedxi1 + dtaperx * dedxk1 - sk1 
				    * d2taperxx;
			    hessx_ref(2, i1) = hessx_ref(2, i1) - sk1 * 
				    dtaperx * dedyi1 + dtapery * dedxk1 - sk1 
				    * d2taperxy;
			    hessx_ref(3, i1) = hessx_ref(3, i1) - sk1 * 
				    dtaperx * dedzi1 + dtaperz * dedxk1 - sk1 
				    * d2taperxz;
			    hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * 
				    dtaperx * dedxk1 - sk1 * dtaperx * dedxk1 
				    + sk1 * sk1 * d2taperxx;
			    hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * 
				    dtaperx * dedyk1 - sk1 * dtapery * dedxk1 
				    + sk1 * sk1 * d2taperxy;
			    hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * 
				    dtaperx * dedzk1 - sk1 * dtaperz * dedxk1 
				    + sk1 * sk1 * d2taperxz;
			    hessx_ref(1, k2) = hessx_ref(1, k2) - sk1 * 
				    dtaperx * dedxk2 - sk2 * dtaperx * dedxk1 
				    + sk1 * sk2 * d2taperxx;
			    hessx_ref(2, k2) = hessx_ref(2, k2) - sk1 * 
				    dtaperx * dedyk2 - sk2 * dtapery * dedxk1 
				    + sk1 * sk2 * d2taperxy;
			    hessx_ref(3, k2) = hessx_ref(3, k2) - sk1 * 
				    dtaperx * dedzk2 - sk2 * dtaperz * dedxk1 
				    + sk1 * sk2 * d2taperxz;
			    hessy_ref(1, i1) = hessy_ref(1, i1) - sk1 * 
				    dtapery * dedxi1 + dtaperx * dedyk1 - sk1 
				    * d2taperxy;
			    hessy_ref(2, i1) = hessy_ref(2, i1) - sk1 * 
				    dtapery * dedyi1 + dtapery * dedyk1 - sk1 
				    * d2taperyy;
			    hessy_ref(3, i1) = hessy_ref(3, i1) - sk1 * 
				    dtapery * dedzi1 + dtaperz * dedyk1 - sk1 
				    * d2taperyz;
			    hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * 
				    dtapery * dedxk1 - sk1 * dtaperx * dedyk1 
				    + sk1 * sk1 * d2taperxy;
			    hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * 
				    dtapery * dedyk1 - sk1 * dtapery * dedyk1 
				    + sk1 * sk1 * d2taperyy;
			    hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * 
				    dtapery * dedzk1 - sk1 * dtaperz * dedyk1 
				    + sk1 * sk1 * d2taperyz;
			    hessy_ref(1, k2) = hessy_ref(1, k2) - sk1 * 
				    dtapery * dedxk2 - sk2 * dtaperx * dedyk1 
				    + sk1 * sk2 * d2taperxy;
			    hessy_ref(2, k2) = hessy_ref(2, k2) - sk1 * 
				    dtapery * dedyk2 - sk2 * dtapery * dedyk1 
				    + sk1 * sk2 * d2taperyy;
			    hessy_ref(3, k2) = hessy_ref(3, k2) - sk1 * 
				    dtapery * dedzk2 - sk2 * dtaperz * dedyk1 
				    + sk1 * sk2 * d2taperyz;
			    hessz_ref(1, i1) = hessz_ref(1, i1) - sk1 * 
				    dtaperz * dedxi1 + dtaperx * dedzk1 - sk1 
				    * d2taperxz;
			    hessz_ref(2, i1) = hessz_ref(2, i1) - sk1 * 
				    dtaperz * dedyi1 + dtapery * dedzk1 - sk1 
				    * d2taperyz;
			    hessz_ref(3, i1) = hessz_ref(3, i1) - sk1 * 
				    dtaperz * dedzi1 + dtaperz * dedzk1 - sk1 
				    * d2taperzz;
			    hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * 
				    dtaperz * dedxk1 - sk1 * dtaperx * dedzk1 
				    + sk1 * sk1 * d2taperxz;
			    hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * 
				    dtaperz * dedyk1 - sk1 * dtapery * dedzk1 
				    + sk1 * sk1 * d2taperyz;
			    hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * 
				    dtaperz * dedzk1 - sk1 * dtaperz * dedzk1 
				    + sk1 * sk1 * d2taperzz;
			    hessz_ref(1, k2) = hessz_ref(1, k2) - sk1 * 
				    dtaperz * dedxk2 - sk2 * dtaperx * dedzk1 
				    + sk1 * sk2 * d2taperxz;
			    hessz_ref(2, k2) = hessz_ref(2, k2) - sk1 * 
				    dtaperz * dedyk2 - sk2 * dtapery * dedzk1 
				    + sk1 * sk2 * d2taperyz;
			    hessz_ref(3, k2) = hessz_ref(3, k2) - sk1 * 
				    dtaperz * dedzk2 - sk2 * dtaperz * dedzk1 
				    + sk1 * sk2 * d2taperzz;
			} else if (r2 > shunt_1.cut2 && *i__ == k2) {
			    hessx_ref(1, i1) = hessx_ref(1, i1) - sk2 * 
				    dtaperx * dedxi1 + dtaperx * dedxk2 - sk2 
				    * d2taperxx;
			    hessx_ref(2, i1) = hessx_ref(2, i1) - sk2 * 
				    dtaperx * dedyi1 + dtapery * dedxk2 - sk2 
				    * d2taperxy;
			    hessx_ref(3, i1) = hessx_ref(3, i1) - sk2 * 
				    dtaperx * dedzi1 + dtaperz * dedxk2 - sk2 
				    * d2taperxz;
			    hessx_ref(1, k1) = hessx_ref(1, k1) - sk2 * 
				    dtaperx * dedxk1 - sk1 * dtaperx * dedxk2 
				    + sk1 * sk2 * d2taperxx;
			    hessx_ref(2, k1) = hessx_ref(2, k1) - sk2 * 
				    dtaperx * dedyk1 - sk1 * dtapery * dedxk2 
				    + sk1 * sk2 * d2taperxy;
			    hessx_ref(3, k1) = hessx_ref(3, k1) - sk2 * 
				    dtaperx * dedzk1 - sk1 * dtaperz * dedxk2 
				    + sk1 * sk2 * d2taperxz;
			    hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * 
				    dtaperx * dedxk2 - sk2 * dtaperx * dedxk2 
				    + sk2 * sk2 * d2taperxx;
			    hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * 
				    dtaperx * dedyk2 - sk2 * dtapery * dedxk2 
				    + sk2 * sk2 * d2taperxy;
			    hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * 
				    dtaperx * dedzk2 - sk2 * dtaperz * dedxk2 
				    + sk2 * sk2 * d2taperxz;
			    hessy_ref(1, i1) = hessy_ref(1, i1) - sk2 * 
				    dtapery * dedxi1 + dtaperx * dedyk2 - sk2 
				    * d2taperxy;
			    hessy_ref(2, i1) = hessy_ref(2, i1) - sk2 * 
				    dtapery * dedyi1 + dtapery * dedyk2 - sk2 
				    * d2taperyy;
			    hessy_ref(3, i1) = hessy_ref(3, i1) - sk2 * 
				    dtapery * dedzi1 + dtaperz * dedyk2 - sk2 
				    * d2taperyz;
			    hessy_ref(1, k1) = hessy_ref(1, k1) - sk2 * 
				    dtapery * dedxk1 - sk1 * dtaperx * dedyk2 
				    + sk1 * sk2 * d2taperxy;
			    hessy_ref(2, k1) = hessy_ref(2, k1) - sk2 * 
				    dtapery * dedyk1 - sk1 * dtapery * dedyk2 
				    + sk1 * sk2 * d2taperyy;
			    hessy_ref(3, k1) = hessy_ref(3, k1) - sk2 * 
				    dtapery * dedzk1 - sk1 * dtaperz * dedyk2 
				    + sk1 * sk2 * d2taperyz;
			    hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * 
				    dtapery * dedxk2 - sk2 * dtaperx * dedyk2 
				    + sk2 * sk2 * d2taperxy;
			    hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * 
				    dtapery * dedyk2 - sk2 * dtapery * dedyk2 
				    + sk2 * sk2 * d2taperyy;
			    hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * 
				    dtapery * dedzk2 - sk2 * dtaperz * dedyk2 
				    + sk2 * sk2 * d2taperyz;
			    hessz_ref(1, i1) = hessz_ref(1, i1) - sk2 * 
				    dtaperz * dedxi1 + dtaperx * dedzk2 - sk2 
				    * d2taperxz;
			    hessz_ref(2, i1) = hessz_ref(2, i1) - sk2 * 
				    dtaperz * dedyi1 + dtapery * dedzk2 - sk2 
				    * d2taperyz;
			    hessz_ref(3, i1) = hessz_ref(3, i1) - sk2 * 
				    dtaperz * dedzi1 + dtaperz * dedzk2 - sk2 
				    * d2taperzz;
			    hessz_ref(1, k1) = hessz_ref(1, k1) - sk2 * 
				    dtaperz * dedxk1 - sk1 * dtaperx * dedzk2 
				    + sk1 * sk2 * d2taperxz;
			    hessz_ref(2, k1) = hessz_ref(2, k1) - sk2 * 
				    dtaperz * dedyk1 - sk1 * dtapery * dedzk2 
				    + sk1 * sk2 * d2taperyz;
			    hessz_ref(3, k1) = hessz_ref(3, k1) - sk2 * 
				    dtaperz * dedzk1 - sk1 * dtaperz * dedzk2 
				    + sk1 * sk2 * d2taperzz;
			    hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * 
				    dtaperz * dedxk2 - sk2 * dtaperx * dedzk2 
				    + sk2 * sk2 * d2taperxz;
			    hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * 
				    dtaperz * dedyk2 - sk2 * dtapery * dedzk2 
				    + sk2 * sk2 * d2taperyz;
			    hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * 
				    dtaperz * dedzk2 - sk2 * dtaperz * dedzk2 
				    + sk2 * sk2 * d2taperzz;
			}
		    }
		}
	    }
	}
L40:
	;
    }
    return 0;
} /* echgdpl2_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef idpl_ref
#undef i12_ref


