/* edipole.f -- translated by f2c (version 20050501).
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
    doublereal electric, dielec, ebuffer, c2scale, c3scale, c4scale, c5scale;
    logical neutnbr, neutcut;
} chgpot_;

#define chgpot_1 chgpot_

struct {
    doublereal bdpl[50000], sdpl[50000];
    integer ndipole, idpl[100000]	/* was [2][50000] */;
} dipole_;

#define dipole_1 dipole_

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

/* Table of constant values */

static integer c_n1 = -1;
static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine edipole  --  dipole-dipole potential energy  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "edipole" calculates the dipole-dipole interaction energy */


/* Subroutine */ int edipole_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__;
    static integer i1, i2, k1, k2;
    static doublereal r2, r3, r4, r5, fi, xi, yi, zi, xk, yk, zk, xq, yq, zq, 
	    xr, yr, zr, ri2, rk2, fik, fgrp, doti, dotk, dotp;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal taper, rirkr3;
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *), switch_(char *, ftnlen), groups_(
	    logical *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *);
    static logical proceed;


#define idpl_ref(a_1,a_2) dipole_1.idpl[(a_2)*2 + a_1 - 3]



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
/*     ##  dipole.i  --  atom & bond dipoles for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     bdpl      magnitude of each of the dipoles (Debyes) */
/*     sdpl      position of each dipole between defining atoms */
/*     ndipole   total number of dipoles in the system */
/*     idpl      numbers of atoms that define each dipole */




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




/*     zero out the overall dipole interaction energy */
/*     and set up the constants for the calculation */

    energi_1.ed = 0.;
    if (dipole_1.ndipole == 0) {
	return 0;
    }

/*     set conversion factor and switching function coefficients */

    f = chgpot_1.electric / (chgpot_1.dielec * 23.070826304099999);
    switch_("DIPOLE", (ftnlen)6);

/*     calculate the pairwise dipole interaction energy term */

    i__1 = dipole_1.ndipole - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = idpl_ref(1, i__);
	i2 = idpl_ref(2, i__);
	xi = atoms_1.x[i2 - 1] - atoms_1.x[i1 - 1];
	yi = atoms_1.y[i2 - 1] - atoms_1.y[i1 - 1];
	zi = atoms_1.z__[i2 - 1] - atoms_1.z__[i1 - 1];
	if (bound_1.use_polymer__) {
	    imager_(&xi, &yi, &zi, &c_n1);
	}
	ri2 = xi * xi + yi * yi + zi * zi;
	xq = atoms_1.x[i1 - 1] + xi * dipole_1.sdpl[i__ - 1];
	yq = atoms_1.y[i1 - 1] + yi * dipole_1.sdpl[i__ - 1];
	zq = atoms_1.z__[i1 - 1] + zi * dipole_1.sdpl[i__ - 1];
	fi = f * dipole_1.bdpl[i__ - 1];

/*     decide whether to compute the current interaction */

	i__2 = dipole_1.ndipole;
	for (k = i__ + 1; k <= i__2; ++k) {
	    k1 = idpl_ref(1, k);
	    k2 = idpl_ref(2, k);
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i1, &i2, &k1, &k2, &c__0, &c__0);
	    }
	    if (proceed) {
		proceed = usage_1.use[i1 - 1] || usage_1.use[i2 - 1] || 
			usage_1.use[k1 - 1] || usage_1.use[k2 - 1];
	    }
	    if (proceed) {
		proceed = k1 != i1 && k1 != i2 && k2 != i1 && k2 != i2;
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xk = atoms_1.x[k2 - 1] - atoms_1.x[k1 - 1];
		yk = atoms_1.y[k2 - 1] - atoms_1.y[k1 - 1];
		zk = atoms_1.z__[k2 - 1] - atoms_1.z__[k1 - 1];
		if (bound_1.use_polymer__) {
		    imager_(&xk, &yk, &zk, &c_n1);
		}
		xr = xq - atoms_1.x[k1 - 1] - xk * dipole_1.sdpl[k - 1];
		yr = yq - atoms_1.y[k1 - 1] - yk * dipole_1.sdpl[k - 1];
		zr = zq - atoms_1.z__[k1 - 1] - zk * dipole_1.sdpl[k - 1];
		image_(&xr, &yr, &zr);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    rk2 = xk * xk + yk * yk + zk * zk;
		    rirkr3 = sqrt(ri2 * rk2 * r2) * r2;
		    dotp = xi * xk + yi * yk + zi * zk;
		    doti = xi * xr + yi * yr + zi * zr;
		    dotk = xk * xr + yk * yr + zk * zr;
		    fik = fi * dipole_1.bdpl[k - 1];
		    e = fik * (dotp - doti * 3. * dotk / r2) / rirkr3;

/*     use energy switching if near the cutoff distance */

		    if (r2 > shunt_1.cut2) {
			r__ = sqrt(r2);
			r3 = r2 * r__;
			r4 = r2 * r2;
			r5 = r2 * r3;
			taper = shunt_1.c5 * r5 + shunt_1.c4 * r4 + 
				shunt_1.c3 * r3 + shunt_1.c2 * r2 + 
				shunt_1.c1 * r__ + shunt_1.c0;
			e *= taper;
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
		    }

/*     increment the overall dipole-dipole energy component */

		    energi_1.ed += e;
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

    i__1 = dipole_1.ndipole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = idpl_ref(1, i__);
	i2 = idpl_ref(2, i__);
	xi = atoms_1.x[i2 - 1] - atoms_1.x[i1 - 1];
	yi = atoms_1.y[i2 - 1] - atoms_1.y[i1 - 1];
	zi = atoms_1.z__[i2 - 1] - atoms_1.z__[i1 - 1];
	if (bound_1.use_polymer__) {
	    imager_(&xi, &yi, &zi, &c_n1);
	}
	ri2 = xi * xi + yi * yi + zi * zi;
	xq = atoms_1.x[i1 - 1] + xi * dipole_1.sdpl[i__ - 1];
	yq = atoms_1.y[i1 - 1] + yi * dipole_1.sdpl[i__ - 1];
	zq = atoms_1.z__[i1 - 1] + zi * dipole_1.sdpl[i__ - 1];
	fi = f * dipole_1.bdpl[i__ - 1];

/*     decide whether to compute the current interaction */

	i__2 = dipole_1.ndipole;
	for (k = i__; k <= i__2; ++k) {
	    k1 = idpl_ref(1, k);
	    k2 = idpl_ref(2, k);
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i1, &i2, &k1, &k2, &c__0, &c__0);
	    }
	    if (proceed) {
		proceed = usage_1.use[i1 - 1] || usage_1.use[i2 - 1] || 
			usage_1.use[k1 - 1] || usage_1.use[k2 - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		i__3 = cell_1.ncell;
		for (j = 1; j <= i__3; ++j) {
		    xk = atoms_1.x[k2 - 1] - atoms_1.x[k1 - 1];
		    yk = atoms_1.y[k2 - 1] - atoms_1.y[k1 - 1];
		    zk = atoms_1.z__[k2 - 1] - atoms_1.z__[k1 - 1];
		    if (bound_1.use_polymer__) {
			imager_(&xk, &yk, &zk, &c_n1);
		    }
		    xr = xq - atoms_1.x[k1 - 1] - xk * dipole_1.sdpl[k - 1];
		    yr = yq - atoms_1.y[k1 - 1] - yk * dipole_1.sdpl[k - 1];
		    zr = zq - atoms_1.z__[k1 - 1] - zk * dipole_1.sdpl[k - 1];
		    imager_(&xr, &yr, &zr, &j);
		    r2 = xr * xr + yr * yr + zr * zr;
		    if (r2 <= shunt_1.off2) {
			rk2 = xk * xk + yk * yk + zk * zk;
			rirkr3 = sqrt(ri2 * rk2 * r2) * r2;
			dotp = xi * xk + yi * yk + zi * zk;
			doti = xi * xr + yi * yr + zi * zr;
			dotk = xk * xr + yk * yr + zk * zr;
			fik = fi * dipole_1.bdpl[k - 1];
			if (bound_1.use_polymer__) {
			    if (r2 < bound_1.polycut2) {
				if (k1 == i1 || k1 == i2 || k2 == i1 || k2 == 
					i2) {
				    fik = 0.;
				}
			    }
			}
			e = fik * (dotp - doti * 3. * dotk / r2) / rirkr3;

/*     use energy switching if near the cutoff distance */

			if (r2 > shunt_1.cut2) {
			    r__ = sqrt(r2);
			    r3 = r2 * r__;
			    r4 = r2 * r2;
			    r5 = r2 * r3;
			    taper = shunt_1.c5 * r5 + shunt_1.c4 * r4 + 
				    shunt_1.c3 * r3 + shunt_1.c2 * r2 + 
				    shunt_1.c1 * r__ + shunt_1.c0;
			    e *= taper;
			}

/*     scale the interaction based on its group membership */

			if (group_1.use_group__) {
			    e *= fgrp;
			}

/*     increment the overall dipole-dipole energy component */

			if (i__ == k) {
			    e *= .5;
			}
			energi_1.ed += e;
		    }
		}
	    }
	}
    }
    return 0;
} /* edipole_ */

#undef idpl_ref


