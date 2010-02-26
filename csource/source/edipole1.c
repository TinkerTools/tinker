/* edipole1.f -- translated by f2c (version 20050501).
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

/* Table of constant values */

static integer c_n1 = -1;
static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine edipole1  --  dipole-dipole energy & derivs  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "edipole1" calculates the dipole-dipole interaction energy */
/*     and first derivatives with respect to Cartesian coordinates */


/* Subroutine */ int edipole1_(void)
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
    static doublereal r2, r3, r4, r5, de, fi, dedrirkri2, dedrirkrk2, xi, yi, 
	    zi, xk, yk, zk, xq, yq, zq, xr, yr, zr, ri2, si1, rk2, si2, sk1, 
	    sk2, xq1, yq1, zq1, xq2, yq2, zq2, fik, vxx, vyx, vyy, vzx, vzz, 
	    vzy, dedr, fgrp, doti, dotk, dotp;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal taper, termx, termy, termz, dedxi1, dedyi1, dedzi1, 
	    dedxi2, dedyi2, dedzi2, dedxk1, dedyk1, dedzk1, dedxk2, dedyk2, 
	    dedzk2, rirkr3;
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal dtaper;
    extern /* Subroutine */ int switch_(char *, ftnlen);
    static doublereal termxi, termyi, termzi, termxk, termyk, termzk;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal deddoti, deddotk;
    static logical proceed;
    static doublereal deddotp, dedrirk, dtaperx, dtapery, dtaperz;


#define ded_ref(a_1,a_2) deriv_1.ded[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
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




/*     zero out the overall dipole interaction energy and derivs, */
/*     then set up the constants for the calculation */

    energi_1.ed = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ded_ref(1, i__) = 0.;
	ded_ref(2, i__) = 0.;
	ded_ref(3, i__) = 0.;
    }
    if (dipole_1.ndipole == 0) {
	return 0;
    }

/*     set conversion factor and switching function coefficients */

    f = chgpot_1.electric / (chgpot_1.dielec * 23.070826304099999);
    switch_("DIPOLE", (ftnlen)6);

/*     compute the dipole interaction energy and first derivatives */

    i__1 = dipole_1.ndipole - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = idpl_ref(1, i__);
	i2 = idpl_ref(2, i__);
	si1 = 1. - dipole_1.sdpl[i__ - 1];
	si2 = dipole_1.sdpl[i__ - 1];
	xi = atoms_1.x[i2 - 1] - atoms_1.x[i1 - 1];
	yi = atoms_1.y[i2 - 1] - atoms_1.y[i1 - 1];
	zi = atoms_1.z__[i2 - 1] - atoms_1.z__[i1 - 1];
	if (bound_1.use_polymer__) {
	    imager_(&xi, &yi, &zi, &c_n1);
	}
	ri2 = xi * xi + yi * yi + zi * zi;
	xq = atoms_1.x[i1 - 1] + xi * si2;
	yq = atoms_1.y[i1 - 1] + yi * si2;
	zq = atoms_1.z__[i1 - 1] + zi * si2;
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
		sk1 = 1. - dipole_1.sdpl[k - 1];
		sk2 = dipole_1.sdpl[k - 1];
		xk = atoms_1.x[k2 - 1] - atoms_1.x[k1 - 1];
		yk = atoms_1.y[k2 - 1] - atoms_1.y[k1 - 1];
		zk = atoms_1.z__[k2 - 1] - atoms_1.z__[k1 - 1];
		if (bound_1.use_polymer__) {
		    imager_(&xk, &yk, &zk, &c_n1);
		}
		xr = xq - atoms_1.x[k1 - 1] - xk * sk2;
		yr = yq - atoms_1.y[k1 - 1] - yk * sk2;
		zr = zq - atoms_1.z__[k1 - 1] - zk * sk2;
		image_(&xr, &yr, &zr);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    rk2 = xk * xk + yk * yk + zk * zk;
		    rirkr3 = sqrt(ri2 * rk2 * r2) * r2;
		    dotp = xi * xk + yi * yk + zi * zk;
		    doti = xi * xr + yi * yr + zi * zr;
		    dotk = xk * xr + yk * yr + zk * zr;
		    fik = fi * dipole_1.bdpl[k - 1];

/*     form the energy and master chain rule term for derivatives */

		    e = fik * (dotp - doti * 3. * dotk / r2) / rirkr3;
		    de = -fik / (rirkr3 * r2);

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
			de *= fgrp;
		    }

/*     secondary chain rule terms for derivative expressions */

		    deddotp = -de * r2;
		    deddoti = de * 3. * dotk;
		    deddotk = de * 3. * doti;
		    dedr = de * (dotp * 3. - doti * 15. * dotk / r2);
		    dedrirk = -e;
		    dedrirkri2 = dedrirk / ri2;
		    dedrirkrk2 = dedrirk / rk2;

/*     more chain rule terms for derivative expressions */

		    termx = dedr * xr + deddoti * xi + deddotk * xk;
		    termy = dedr * yr + deddoti * yi + deddotk * yk;
		    termz = dedr * zr + deddoti * zi + deddotk * zk;
		    termxi = dedrirkri2 * xi + deddotp * xk + deddoti * xr;
		    termyi = dedrirkri2 * yi + deddotp * yk + deddoti * yr;
		    termzi = dedrirkri2 * zi + deddotp * zk + deddoti * zr;
		    termxk = dedrirkrk2 * xk + deddotp * xi + deddotk * xr;
		    termyk = dedrirkrk2 * yk + deddotp * yi + deddotk * yr;
		    termzk = dedrirkrk2 * zk + deddotp * zi + deddotk * zr;

/*     finally, the individual first derivative components */

		    dedxi1 = si1 * termx - termxi;
		    dedyi1 = si1 * termy - termyi;
		    dedzi1 = si1 * termz - termzi;
		    dedxi2 = si2 * termx + termxi;
		    dedyi2 = si2 * termy + termyi;
		    dedzi2 = si2 * termz + termzi;
		    dedxk1 = -sk1 * termx - termxk;
		    dedyk1 = -sk1 * termy - termyk;
		    dedzk1 = -sk1 * termz - termzk;
		    dedxk2 = -sk2 * termx + termxk;
		    dedyk2 = -sk2 * termy + termyk;
		    dedzk2 = -sk2 * termz + termzk;

/*     use energy switching if near the cutoff distance */

		    if (r2 > shunt_1.cut2) {
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
			dtaper = dtaper * e / r__;
			dtaperx = xr * dtaper;
			dtapery = yr * dtaper;
			dtaperz = zr * dtaper;
			e *= taper;
			dedxi1 = dedxi1 * taper + si1 * dtaperx;
			dedyi1 = dedyi1 * taper + si1 * dtapery;
			dedzi1 = dedzi1 * taper + si1 * dtaperz;
			dedxi2 = dedxi2 * taper + si2 * dtaperx;
			dedyi2 = dedyi2 * taper + si2 * dtapery;
			dedzi2 = dedzi2 * taper + si2 * dtaperz;
			dedxk1 = dedxk1 * taper - sk1 * dtaperx;
			dedyk1 = dedyk1 * taper - sk1 * dtapery;
			dedzk1 = dedzk1 * taper - sk1 * dtaperz;
			dedxk2 = dedxk2 * taper - sk2 * dtaperx;
			dedyk2 = dedyk2 * taper - sk2 * dtapery;
			dedzk2 = dedzk2 * taper - sk2 * dtaperz;
		    }

/*     increment the overall energy and derivative expressions */

		    energi_1.ed += e;
		    ded_ref(1, i1) = ded_ref(1, i1) + dedxi1;
		    ded_ref(2, i1) = ded_ref(2, i1) + dedyi1;
		    ded_ref(3, i1) = ded_ref(3, i1) + dedzi1;
		    ded_ref(1, i2) = ded_ref(1, i2) + dedxi2;
		    ded_ref(2, i2) = ded_ref(2, i2) + dedyi2;
		    ded_ref(3, i2) = ded_ref(3, i2) + dedzi2;
		    ded_ref(1, k1) = ded_ref(1, k1) + dedxk1;
		    ded_ref(2, k1) = ded_ref(2, k1) + dedyk1;
		    ded_ref(3, k1) = ded_ref(3, k1) + dedzk1;
		    ded_ref(1, k2) = ded_ref(1, k2) + dedxk2;
		    ded_ref(2, k2) = ded_ref(2, k2) + dedyk2;
		    ded_ref(3, k2) = ded_ref(3, k2) + dedzk2;

/*     increment the internal virial tensor components */

		    xq1 = atoms_1.x[k1 - 1] - xq;
		    yq1 = atoms_1.y[k1 - 1] - yq;
		    zq1 = atoms_1.z__[k1 - 1] - zq;
		    xq2 = atoms_1.x[k2 - 1] - xq;
		    yq2 = atoms_1.y[k2 - 1] - yq;
		    zq2 = atoms_1.z__[k2 - 1] - zq;
		    vxx = xq1 * dedxk1 + xq2 * dedxk2;
		    vyx = (yq1 * dedxk1 + yq2 * dedxk2 + xq1 * dedyk1 + xq2 * 
			    dedyk2) * .5;
		    vzx = (zq1 * dedxk1 + zq2 * dedxk2 + xq1 * dedzk1 + xq2 * 
			    dedzk2) * .5;
		    vyy = yq1 * dedyk1 + yq2 * dedyk2;
		    vzy = (zq1 * dedyk1 + zq2 * dedyk2 + yq1 * dedzk1 + yq2 * 
			    dedzk2) * .5;
		    vzz = zq1 * dedzk1 + zq2 * dedzk2;
		    vir_ref(1, 1) = vir_ref(1, 1) + vxx;
		    vir_ref(2, 1) = vir_ref(2, 1) + vyx;
		    vir_ref(3, 1) = vir_ref(3, 1) + vzx;
		    vir_ref(1, 2) = vir_ref(1, 2) + vyx;
		    vir_ref(2, 2) = vir_ref(2, 2) + vyy;
		    vir_ref(3, 2) = vir_ref(3, 2) + vzy;
		    vir_ref(1, 3) = vir_ref(1, 3) + vzx;
		    vir_ref(2, 3) = vir_ref(2, 3) + vzy;
		    vir_ref(3, 3) = vir_ref(3, 3) + vzz;

/*     increment the total intermolecular energy */

		    if (molcul_1.molcule[i1 - 1] != molcul_1.molcule[k1 - 1]) 
			    {
			inter_1.einter += e;
		    }
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
	si1 = 1. - dipole_1.sdpl[i__ - 1];
	si2 = dipole_1.sdpl[i__ - 1];
	xi = atoms_1.x[i2 - 1] - atoms_1.x[i1 - 1];
	yi = atoms_1.y[i2 - 1] - atoms_1.y[i1 - 1];
	zi = atoms_1.z__[i2 - 1] - atoms_1.z__[i1 - 1];
	if (bound_1.use_polymer__) {
	    imager_(&xi, &yi, &zi, &c_n1);
	}
	ri2 = xi * xi + yi * yi + zi * zi;
	xq = atoms_1.x[i1 - 1] + xi * si2;
	yq = atoms_1.y[i1 - 1] + yi * si2;
	zq = atoms_1.z__[i1 - 1] + zi * si2;
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
		sk1 = 1. - dipole_1.sdpl[k - 1];
		sk2 = dipole_1.sdpl[k - 1];
		i__3 = cell_1.ncell;
		for (j = 1; j <= i__3; ++j) {
		    xk = atoms_1.x[k2 - 1] - atoms_1.x[k1 - 1];
		    yk = atoms_1.y[k2 - 1] - atoms_1.y[k1 - 1];
		    zk = atoms_1.z__[k2 - 1] - atoms_1.z__[k1 - 1];
		    if (bound_1.use_polymer__) {
			imager_(&xk, &yk, &zk, &c_n1);
		    }
		    xr = xq - atoms_1.x[k1 - 1] - xk * sk2;
		    yr = yq - atoms_1.y[k1 - 1] - yk * sk2;
		    zr = zq - atoms_1.z__[k1 - 1] - zk * sk2;
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

/*     form the energy and master chain rule term for derivatives */

			e = fik * (dotp - doti * 3. * dotk / r2) / rirkr3;
			de = -fik / (rirkr3 * r2);

/*     scale the interaction based on its group membership */

			if (group_1.use_group__) {
			    e *= fgrp;
			    de *= fgrp;
			}

/*     secondary chain rule terms for derivative expressions */

			deddotp = -de * r2;
			deddoti = de * 3. * dotk;
			deddotk = de * 3. * doti;
			dedr = de * (dotp * 3. - doti * 15. * dotk / r2);
			dedrirk = -e;
			dedrirkri2 = dedrirk / ri2;
			dedrirkrk2 = dedrirk / rk2;

/*     more chain rule terms for derivative expressions */

			termx = dedr * xr + deddoti * xi + deddotk * xk;
			termy = dedr * yr + deddoti * yi + deddotk * yk;
			termz = dedr * zr + deddoti * zi + deddotk * zk;
			termxi = dedrirkri2 * xi + deddotp * xk + deddoti * 
				xr;
			termyi = dedrirkri2 * yi + deddotp * yk + deddoti * 
				yr;
			termzi = dedrirkri2 * zi + deddotp * zk + deddoti * 
				zr;
			termxk = dedrirkrk2 * xk + deddotp * xi + deddotk * 
				xr;
			termyk = dedrirkrk2 * yk + deddotp * yi + deddotk * 
				yr;
			termzk = dedrirkrk2 * zk + deddotp * zi + deddotk * 
				zr;

/*     finally, the individual first derivative components */

			dedxi1 = si1 * termx - termxi;
			dedyi1 = si1 * termy - termyi;
			dedzi1 = si1 * termz - termzi;
			dedxi2 = si2 * termx + termxi;
			dedyi2 = si2 * termy + termyi;
			dedzi2 = si2 * termz + termzi;
			dedxk1 = -sk1 * termx - termxk;
			dedyk1 = -sk1 * termy - termyk;
			dedzk1 = -sk1 * termz - termzk;
			dedxk2 = -sk2 * termx + termxk;
			dedyk2 = -sk2 * termy + termyk;
			dedzk2 = -sk2 * termz + termzk;

/*     use energy switching if near the cutoff distance */

			if (r2 > shunt_1.cut2) {
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
			    dtaper = dtaper * e / r__;
			    dtaperx = xr * dtaper;
			    dtapery = yr * dtaper;
			    dtaperz = zr * dtaper;
			    e *= taper;
			    dedxi1 = dedxi1 * taper + si1 * dtaperx;
			    dedyi1 = dedyi1 * taper + si1 * dtapery;
			    dedzi1 = dedzi1 * taper + si1 * dtaperz;
			    dedxi2 = dedxi2 * taper + si2 * dtaperx;
			    dedyi2 = dedyi2 * taper + si2 * dtapery;
			    dedzi2 = dedzi2 * taper + si2 * dtaperz;
			    dedxk1 = dedxk1 * taper - sk1 * dtaperx;
			    dedyk1 = dedyk1 * taper - sk1 * dtapery;
			    dedzk1 = dedzk1 * taper - sk1 * dtaperz;
			    dedxk2 = dedxk2 * taper - sk2 * dtaperx;
			    dedyk2 = dedyk2 * taper - sk2 * dtapery;
			    dedzk2 = dedzk2 * taper - sk2 * dtaperz;
			}

/*     increment the overall energy and derivative expressions */

			if (i__ == k) {
			    e *= .5;
			}
			energi_1.ed += e;
			ded_ref(1, i1) = ded_ref(1, i1) + dedxi1;
			ded_ref(2, i1) = ded_ref(2, i1) + dedyi1;
			ded_ref(3, i1) = ded_ref(3, i1) + dedzi1;
			ded_ref(1, i2) = ded_ref(1, i2) + dedxi2;
			ded_ref(2, i2) = ded_ref(2, i2) + dedyi2;
			ded_ref(3, i2) = ded_ref(3, i2) + dedzi2;
			if (i__ != k) {
			    ded_ref(1, k1) = ded_ref(1, k1) + dedxk1;
			    ded_ref(2, k1) = ded_ref(2, k1) + dedyk1;
			    ded_ref(3, k1) = ded_ref(3, k1) + dedzk1;
			    ded_ref(1, k2) = ded_ref(1, k2) + dedxk2;
			    ded_ref(2, k2) = ded_ref(2, k2) + dedyk2;
			    ded_ref(3, k2) = ded_ref(3, k2) + dedzk2;
			}

/*     increment the internal virial tensor components */

			xq1 = atoms_1.x[k1 - 1] - xq;
			yq1 = atoms_1.y[k1 - 1] - yq;
			zq1 = atoms_1.z__[k1 - 1] - zq;
			xq2 = atoms_1.x[k2 - 1] - xq;
			yq2 = atoms_1.y[k2 - 1] - yq;
			zq2 = atoms_1.z__[k2 - 1] - zq;
			vxx = xq1 * dedxk1 + xq2 * dedxk2;
			vyx = (yq1 * dedxk1 + yq2 * dedxk2 + xq1 * dedyk1 + 
				xq2 * dedyk2) * .5;
			vzx = (zq1 * dedxk1 + zq2 * dedxk2 + xq1 * dedzk1 + 
				xq2 * dedzk2) * .5;
			vyy = yq1 * dedyk1 + yq2 * dedyk2;
			vzy = (zq1 * dedyk1 + zq2 * dedyk2 + yq1 * dedzk1 + 
				yq2 * dedzk2) * .5;
			vzz = zq1 * dedzk1 + zq2 * dedzk2;
			vir_ref(1, 1) = vir_ref(1, 1) + vxx;
			vir_ref(2, 1) = vir_ref(2, 1) + vyx;
			vir_ref(3, 1) = vir_ref(3, 1) + vzx;
			vir_ref(1, 2) = vir_ref(1, 2) + vyx;
			vir_ref(2, 2) = vir_ref(2, 2) + vyy;
			vir_ref(3, 2) = vir_ref(3, 2) + vzy;
			vir_ref(1, 3) = vir_ref(1, 3) + vzx;
			vir_ref(2, 3) = vir_ref(2, 3) + vzy;
			vir_ref(3, 3) = vir_ref(3, 3) + vzz;

/*     increment the total intermolecular energy */

			inter_1.einter += e;
		    }
		}
	    }
	}
    }
    return 0;
} /* edipole1_ */

#undef idpl_ref
#undef vir_ref
#undef ded_ref


