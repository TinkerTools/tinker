/* edipole2.f -- translated by f2c (version 20050501).
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

static integer c_n1 = -1;
static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine edipole2  --  atomwise dipole-dipole Hessian  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "edipole2" calculates second derivatives of the */
/*     dipole-dipole interaction energy for a single atom */


/* Subroutine */ int edipole2_(integer *i__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal dtxidxi1, dtxidxi2, dtxkdxi1, dtxkdxi2, dtyidxi1, 
	    dtykdxi1, dtyidxi2, dtykdxi2, dtzidxi1, dtzkdxi1, dtzidxi2, 
	    dtzkdxi2, dtxidyi1, dtxkdyi1, dtxidyi2, dtxkdyi2, dtyidyi1, 
	    dtykdyi1, dtyidyi2, dtykdyi2, dtzidyi1, dtzkdyi1, dtzidyi2, 
	    dtzkdyi2, dtxidzi1, dtxkdzi1, dtxidzi2, dtxkdzi2, dtyidzi1, 
	    dtykdzi1, dtyidzi2, dtykdzi2, dtzidzi1, dtzkdzi1, dtzidzi2, 
	    dtzkdzi2, e, f, r__, d2taperxx, d2taperxy, d2taperyy, d2taperxz, 
	    d2taperzz, d2taperyz;
    static integer i1, i2, k1, k2;
    static doublereal r2, r3, r4, r5, de, fi, xi, yi, zi, xk, yk, zk, xq, yq, 
	    zq, xr, yr, zr, ri2, si1, rk2, si2, sk1, sk2, fik, xrr2, yrr2, 
	    zrr2, dedr, fgrp, doti, dotk, enum__, dotp, part, xixk, xiyk, 
	    xizk, yixk, yiyk, yizk, zixk, xixr, xiyr, xizr, yixr, yiyr, yizr, 
	    zixr, ziyr, zizr, xkxr, xkyr, xkzr, ykxr, ykyr, ykzr, zkxr, zkyr, 
	    zkzr, ziyk, zizk, xrxr, xryr, xrzr, yryr, yrzr, zrzr, xiri2, 
	    yiri2, ziri2, r2inv, xkrk2, ykrk2, zkrk2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static integer jcell;
    static doublereal dotik, taper, termx, termy, termz, dedxi1, dedyi1, 
	    dedzi1, dedxi2, dedyi2, dedzi2, dedxk1, dedyk1, dedzk1, dedxk2, 
	    dedyk2, dedzk2, dtdxi1, dtdyi1, dtdzi1, dtdxi2, dtdyi2, dtdzi2, 
	    ri2inv, rirkr3;
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal factor, dtaper, partik, xidotk, yidotk, zidotk, xkdoti, 
	    termxi, termyi, termzi, termxk, termyk, termzk, ykdoti, zkdoti;
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *);
    static doublereal d2taper, dtxdxi1, dtxdxi2, dtydxi1, dtydxi2, dtzdxi1, 
	    dtzdxi2, dtxdyi1, dtxdyi2, dtydyi1, dtydyi2, dtzdyi1, dtzdyi2, 
	    dtxdzi1, dtxdzi2, dtydzi1, dtydzi2, dtzdzi1, dtzdzi2, deddoti, 
	    deddotk;
    static logical proceed;
    static doublereal deddotp, dedrirk;
    static integer idipole, kdipole;
    static doublereal factori, factork, dtaperx, dtapery, dtaperz;


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




/*     set conversion factor and switching function coefficients */

    if (dipole_1.ndipole == 0) {
	return 0;
    }
    f = chgpot_1.electric / (chgpot_1.dielec * 23.070826304099999);
    switch_("DIPOLE", (ftnlen)6);

/*     calculate the dipole interaction energy Hessian elements */

    i__1 = dipole_1.ndipole;
    for (idipole = 1; idipole <= i__1; ++idipole) {
	i1 = idpl_ref(1, idipole);
	i2 = idpl_ref(2, idipole);
	si1 = 1. - dipole_1.sdpl[idipole - 1];
	si2 = dipole_1.sdpl[idipole - 1];
	if (i1 != *i__ && i2 != *i__) {
	    goto L10;
	}
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
	fi = f * dipole_1.bdpl[idipole - 1];

/*     decide whether to compute the current interaction */

	i__2 = dipole_1.ndipole;
	for (kdipole = 1; kdipole <= i__2; ++kdipole) {
	    k1 = idpl_ref(1, kdipole);
	    k2 = idpl_ref(2, kdipole);
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i1, &i2, &k1, &k2, &c__0, &c__0);
	    }
	    if (proceed) {
		proceed = k1 != i1 && k1 != i2 && k2 != i1 && k2 != i2;
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		sk1 = 1. - dipole_1.sdpl[kdipole - 1];
		sk2 = dipole_1.sdpl[kdipole - 1];
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
		    fik = fi * dipole_1.bdpl[kdipole - 1];

/*     some abbreviations used in various chain rule terms */

		    dotik = doti * dotk;
		    enum__ = dotp * r2 - dotik * 3.;
		    r2inv = 15. / r2;
		    ri2inv = 1. / ri2;
		    xrr2 = xr / r2;
		    yrr2 = yr / r2;
		    zrr2 = zr / r2;
		    xiri2 = xi / ri2;
		    yiri2 = yi / ri2;
		    ziri2 = zi / ri2;
		    xkrk2 = xk / rk2;
		    ykrk2 = yk / rk2;
		    zkrk2 = zk / rk2;
		    xixr = xi * xr;
		    xiyr = xi * yr;
		    xizr = xi * zr;
		    yixr = yi * xr;
		    yiyr = yi * yr;
		    yizr = yi * zr;
		    zixr = zi * xr;
		    ziyr = zi * yr;
		    zizr = zi * zr;
		    xkxr = xk * xr;
		    xkyr = xk * yr;
		    xkzr = xk * zr;
		    ykxr = yk * xr;
		    ykyr = yk * yr;
		    ykzr = yk * zr;
		    zkxr = zk * xr;
		    zkyr = zk * yr;
		    zkzr = zk * zr;
		    xixk = xi * xk;
		    xiyk = xi * yk;
		    xizk = xi * zk;
		    yixk = yi * xk;
		    yiyk = yi * yk;
		    yizk = yi * zk;
		    zixk = zi * xk;
		    ziyk = zi * yk;
		    zizk = zi * zk;
		    xrxr = xr * 3. * xr;
		    xryr = xr * 3. * yr;
		    xrzr = xr * 3. * zr;
		    yryr = yr * 3. * yr;
		    yrzr = yr * 3. * zr;
		    zrzr = zr * 3. * zr;
		    xidotk = xi * dotk;
		    yidotk = yi * dotk;
		    zidotk = zi * dotk;
		    xkdoti = xk * doti;
		    ykdoti = yk * doti;
		    zkdoti = zk * doti;

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			fik *= fgrp;
		    }

/*     form the master chain rule term for derivatives */

		    de = -fik / (rirkr3 * r2);

/*     form the chain rule terms for first derivatives */

		    deddotp = -de * r2;
		    deddoti = de * 3. * dotk;
		    deddotk = de * 3. * doti;
		    dedr = de * (dotp * 3. - dotik * 15. / r2);
		    dedrirk = de * enum__;

/*     more first derivative chain rule expressions */

		    termx = dedr * xr + deddoti * xi + deddotk * xk;
		    termy = dedr * yr + deddoti * yi + deddotk * yk;
		    termz = dedr * zr + deddoti * zi + deddotk * zk;
		    termxi = dedrirk * xiri2 + deddotp * xk + deddoti * xr;
		    termyi = dedrirk * yiri2 + deddotp * yk + deddoti * yr;
		    termzi = dedrirk * ziri2 + deddotp * zk + deddoti * zr;
		    termxk = dedrirk * xkrk2 + deddotp * xi + deddotk * xr;
		    termyk = dedrirk * ykrk2 + deddotp * yi + deddotk * yr;
		    termzk = dedrirk * zkrk2 + deddotp * zi + deddotk * zr;

/*     use energy switching if near the cutoff distance */

		    if (r2 > shunt_1.cut2) {
			e = fik * (dotp - doti * 3. * dotk / r2) / rirkr3;
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
			de *= taper;
			termx *= taper;
			termy *= taper;
			termz *= taper;
			termxi *= taper;
			termyi *= taper;
			termzi *= taper;
			termxk *= taper;
			termyk *= taper;
			termzk *= taper;
		    }

/*     chain rule terms for second derivative components */

		    if (*i__ == i1) {
			dtdxi1 = si1 * -5. * xrr2 + xiri2;
			part = si1 * xkdoti - dotk * xr + si1 * xidotk - si1 *
				 2. * dotik * xrr2;
			partik = -xk * r2 + si1 * 2. * dotp * xr - si1 * 3. * 
				xkdoti + xr * 3. * dotk - si1 * 3. * xidotk;
			factor = si1 * 3. * dotp - xkxr * 6. + si1 * 6. * 
				xixk - dotk * 3. - r2inv * (xr * part + si1 * 
				dotik);
			factori = si1 * 3. * dotk + si1 * xkxr + xiri2 * 
				partik - enum__ * (ri2inv - xiri2 * 2. * 
				xiri2);
			factork = r2 + si1 * 3. * doti + si1 * xixr - xrxr + 
				xkrk2 * partik;
			dtxdxi1 = dtdxi1 * termx + de * factor;
			dtxidxi1 = dtdxi1 * termxi + de * factori;
			dtxkdxi1 = dtdxi1 * termxk + de * factork;
			factor = xkyr * -3. - ykxr * 3. + si1 * 3. * xiyk + 
				si1 * 3. * yixk - r2inv * yr * part;
			factori = si1 * -2. * ykxr + si1 * 3. * xkyr + yiri2 *
				 partik + enum__ * 2. * yiri2 * xiri2;
			factork = si1 * -2. * yixr - xryr + si1 * 3. * xiyr + 
				ykrk2 * partik;
			dtydxi1 = dtdxi1 * termy + de * factor;
			dtyidxi1 = dtdxi1 * termyi + de * factori;
			dtykdxi1 = dtdxi1 * termyk + de * factork;
			factor = xkzr * -3. - zkxr * 3. + si1 * 3. * xizk + 
				si1 * 3. * zixk - r2inv * zr * part;
			factori = si1 * -2. * zkxr + si1 * 3. * xkzr + ziri2 *
				 partik + enum__ * 2. * ziri2 * xiri2;
			factork = si1 * -2. * zixr - xrzr + si1 * 3. * xizr + 
				zkrk2 * partik;
			dtzdxi1 = dtdxi1 * termz + de * factor;
			dtzidxi1 = dtdxi1 * termzi + de * factori;
			dtzkdxi1 = dtdxi1 * termzk + de * factork;
			dtdyi1 = si1 * -5. * yrr2 + yiri2;
			part = si1 * ykdoti - dotk * yr + si1 * yidotk - si1 *
				 2. * dotik * yrr2;
			partik = -yk * r2 + si1 * 2. * dotp * yr - si1 * 3. * 
				ykdoti + yr * 3. * dotk - si1 * 3. * yidotk;
			factor = ykxr * -3. - xkyr * 3. + si1 * 3. * yixk + 
				si1 * 3. * xiyk - r2inv * xr * part;
			factori = si1 * -2. * xkyr + si1 * 3. * ykxr + xiri2 *
				 partik + enum__ * 2. * xiri2 * yiri2;
			factork = si1 * -2. * xiyr - xryr + si1 * 3. * yixr + 
				xkrk2 * partik;
			dtxdyi1 = dtdyi1 * termx + de * factor;
			dtxidyi1 = dtdyi1 * termxi + de * factori;
			dtxkdyi1 = dtdyi1 * termxk + de * factork;
			factor = si1 * 3. * dotp - ykyr * 6. + si1 * 6. * 
				yiyk - dotk * 3. - r2inv * (yr * part + si1 * 
				dotik);
			factori = si1 * 3. * dotk + si1 * ykyr + yiri2 * 
				partik - enum__ * (ri2inv - yiri2 * 2. * 
				yiri2);
			factork = r2 + si1 * 3. * doti + si1 * yiyr - yryr + 
				ykrk2 * partik;
			dtydyi1 = dtdyi1 * termy + de * factor;
			dtyidyi1 = dtdyi1 * termyi + de * factori;
			dtykdyi1 = dtdyi1 * termyk + de * factork;
			factor = ykzr * -3. - zkyr * 3. + si1 * 3. * yizk + 
				si1 * 3. * ziyk - r2inv * zr * part;
			factori = si1 * -2. * zkyr + si1 * 3. * ykzr + ziri2 *
				 partik + enum__ * 2. * ziri2 * yiri2;
			factork = si1 * -2. * ziyr - yrzr + si1 * 3. * yizr + 
				zkrk2 * partik;
			dtzdyi1 = dtdyi1 * termz + de * factor;
			dtzidyi1 = dtdyi1 * termzi + de * factori;
			dtzkdyi1 = dtdyi1 * termzk + de * factork;
			dtdzi1 = si1 * -5. * zrr2 + ziri2;
			part = si1 * zkdoti - dotk * zr + si1 * zidotk - si1 *
				 2. * dotik * zrr2;
			partik = -zk * r2 + si1 * 2. * dotp * zr - si1 * 3. * 
				zkdoti + zr * 3. * dotk - si1 * 3. * zidotk;
			factor = zkxr * -3. - xkzr * 3. + si1 * 3. * zixk + 
				si1 * 3. * xizk - r2inv * xr * part;
			factori = si1 * -2. * xkzr + si1 * 3. * zkxr + xiri2 *
				 partik + enum__ * 2. * xiri2 * ziri2;
			factork = si1 * -2. * xizr - xrzr + si1 * 3. * zixr + 
				xkrk2 * partik;
			dtxdzi1 = dtdzi1 * termx + de * factor;
			dtxidzi1 = dtdzi1 * termxi + de * factori;
			dtxkdzi1 = dtdzi1 * termxk + de * factork;
			factor = zkyr * -3. - ykzr * 3. + si1 * 3. * ziyk + 
				si1 * 3. * yizk - r2inv * yr * part;
			factori = si1 * -2. * ykzr + si1 * 3. * zkyr + yiri2 *
				 partik + enum__ * 2. * yiri2 * ziri2;
			factork = si1 * -2. * yizr - yrzr + si1 * 3. * ziyr + 
				ykrk2 * partik;
			dtydzi1 = dtdzi1 * termy + de * factor;
			dtyidzi1 = dtdzi1 * termyi + de * factori;
			dtykdzi1 = dtdzi1 * termyk + de * factork;
			factor = si1 * 3. * dotp - zkzr * 6. + si1 * 6. * 
				zizk - dotk * 3. - r2inv * (zr * part + si1 * 
				dotik);
			factori = si1 * 3. * dotk + si1 * zkzr + ziri2 * 
				partik - enum__ * (ri2inv - ziri2 * 2. * 
				ziri2);
			factork = r2 + si1 * 3. * doti + si1 * zizr - zrzr + 
				zkrk2 * partik;
			dtzdzi1 = dtdzi1 * termz + de * factor;
			dtzidzi1 = dtdzi1 * termzi + de * factori;
			dtzkdzi1 = dtdzi1 * termzk + de * factork;
		    } else if (*i__ == i2) {
			dtdxi2 = si2 * -5. * xrr2 - xiri2;
			part = si2 * xkdoti + dotk * xr + si2 * xidotk - si2 *
				 2. * dotik * xrr2;
			partik = xk * r2 + si2 * 2. * dotp * xr - si2 * 3. * 
				xkdoti - xr * 3. * dotk - si2 * 3. * xidotk;
			factor = si2 * 3. * dotp + xkxr * 6. + si2 * 6. * 
				xixk + dotk * 3. - r2inv * (xr * part + si2 * 
				dotik);
			factori = si2 * 3. * dotk + si2 * xkxr + xiri2 * 
				partik + enum__ * (ri2inv - xiri2 * 2. * 
				xiri2);
			factork = -r2 + si2 * 3. * doti + si2 * xixr + xrxr + 
				xkrk2 * partik;
			dtxdxi2 = dtdxi2 * termx + de * factor;
			dtxidxi2 = dtdxi2 * termxi + de * factori;
			dtxkdxi2 = dtdxi2 * termxk + de * factork;
			factor = xkyr * 3. + ykxr * 3. + si2 * 3. * xiyk + 
				si2 * 3. * yixk - r2inv * yr * part;
			factori = si2 * -2. * ykxr + si2 * 3. * xkyr + yiri2 *
				 partik - enum__ * 2. * yiri2 * xiri2;
			factork = si2 * -2. * yixr + xryr + si2 * 3. * xiyr + 
				ykrk2 * partik;
			dtydxi2 = dtdxi2 * termy + de * factor;
			dtyidxi2 = dtdxi2 * termyi + de * factori;
			dtykdxi2 = dtdxi2 * termyk + de * factork;
			factor = xkzr * 3. + zkxr * 3. + si2 * 3. * xizk + 
				si2 * 3. * zixk - r2inv * zr * part;
			factori = si2 * -2. * zkxr + si2 * 3. * xkzr + ziri2 *
				 partik - enum__ * 2. * ziri2 * xiri2;
			factork = si2 * -2. * zixr + xrzr + si2 * 3. * xizr + 
				zkrk2 * partik;
			dtzdxi2 = dtdxi2 * termz + de * factor;
			dtzidxi2 = dtdxi2 * termzi + de * factori;
			dtzkdxi2 = dtdxi2 * termzk + de * factork;
			dtdyi2 = si2 * -5. * yrr2 - yiri2;
			part = si2 * ykdoti + dotk * yr + si2 * yidotk - si2 *
				 2. * dotik * yrr2;
			partik = yk * r2 + si2 * 2. * dotp * yr - si2 * 3. * 
				ykdoti - yr * 3. * dotk - si2 * 3. * yidotk;
			factor = ykxr * 3. + xkyr * 3. + si2 * 3. * yixk + 
				si2 * 3. * xiyk - r2inv * xr * part;
			factori = si2 * -2. * xkyr + si2 * 3. * ykxr + xiri2 *
				 partik - enum__ * 2. * xiri2 * yiri2;
			factork = si2 * -2. * xiyr + xryr + si2 * 3. * yixr + 
				xkrk2 * partik;
			dtxdyi2 = dtdyi2 * termx + de * factor;
			dtxidyi2 = dtdyi2 * termxi + de * factori;
			dtxkdyi2 = dtdyi2 * termxk + de * factork;
			factor = si2 * 3. * dotp + ykyr * 6. + si2 * 6. * 
				yiyk + dotk * 3. - r2inv * (yr * part + si2 * 
				dotik);
			factori = si2 * 3. * dotk + si2 * ykyr + yiri2 * 
				partik + enum__ * (ri2inv - yiri2 * 2. * 
				yiri2);
			factork = -r2 + si2 * 3. * doti + si2 * yiyr + yryr + 
				ykrk2 * partik;
			dtydyi2 = dtdyi2 * termy + de * factor;
			dtyidyi2 = dtdyi2 * termyi + de * factori;
			dtykdyi2 = dtdyi2 * termyk + de * factork;
			factor = ykzr * 3. + zkyr * 3. + si2 * 3. * yizk + 
				si2 * 3. * ziyk - r2inv * zr * part;
			factori = si2 * -2. * zkyr + si2 * 3. * ykzr + ziri2 *
				 partik - enum__ * 2. * ziri2 * yiri2;
			factork = si2 * -2. * ziyr + yrzr + si2 * 3. * yizr + 
				zkrk2 * partik;
			dtzdyi2 = dtdyi2 * termz + de * factor;
			dtzidyi2 = dtdyi2 * termzi + de * factori;
			dtzkdyi2 = dtdyi2 * termzk + de * factork;
			dtdzi2 = si2 * -5. * zrr2 - ziri2;
			part = si2 * zkdoti + dotk * zr + si2 * zidotk - si2 *
				 2. * dotik * zrr2;
			partik = zk * r2 + si2 * 2. * dotp * zr - si2 * 3. * 
				zkdoti - zr * 3. * dotk - si2 * 3. * zidotk;
			factor = zkxr * 3. + xkzr * 3. + si2 * 3. * zixk + 
				si2 * 3. * xizk - r2inv * xr * part;
			factori = si2 * -2. * xkzr + si2 * 3. * zkxr + xiri2 *
				 partik - enum__ * 2. * xiri2 * ziri2;
			factork = si2 * -2. * xizr + xrzr + si2 * 3. * zixr + 
				xkrk2 * partik;
			dtxdzi2 = dtdzi2 * termx + de * factor;
			dtxidzi2 = dtdzi2 * termxi + de * factori;
			dtxkdzi2 = dtdzi2 * termxk + de * factork;
			factor = zkyr * 3. + ykzr * 3. + si2 * 3. * ziyk + 
				si2 * 3. * yizk - r2inv * yr * part;
			factori = si2 * -2. * ykzr + si2 * 3. * zkyr + yiri2 *
				 partik - enum__ * 2. * yiri2 * ziri2;
			factork = si2 * -2. * yizr + yrzr + si2 * 3. * ziyr + 
				ykrk2 * partik;
			dtydzi2 = dtdzi2 * termy + de * factor;
			dtyidzi2 = dtdzi2 * termyi + de * factori;
			dtykdzi2 = dtdzi2 * termyk + de * factork;
			factor = si2 * 3. * dotp + zkzr * 6. + si2 * 6. * 
				zizk + dotk * 3. - r2inv * (zr * part + si2 * 
				dotik);
			factori = si2 * 3. * dotk + si2 * zkzr + ziri2 * 
				partik + enum__ * (ri2inv - ziri2 * 2. * 
				ziri2);
			factork = -r2 + si2 * 3. * doti + si2 * zizr + zrzr + 
				zkrk2 * partik;
			dtzdzi2 = dtdzi2 * termz + de * factor;
			dtzidzi2 = dtdzi2 * termzi + de * factori;
			dtzkdzi2 = dtdzi2 * termzk + de * factork;
		    }

/*     increment diagonal and off-diagonal Hessian elements */

		    if (*i__ == i1) {
			hessx_ref(1, i1) = hessx_ref(1, i1) + si1 * dtxdxi1 - 
				dtxidxi1;
			hessx_ref(2, i1) = hessx_ref(2, i1) + si1 * dtydxi1 - 
				dtyidxi1;
			hessx_ref(3, i1) = hessx_ref(3, i1) + si1 * dtzdxi1 - 
				dtzidxi1;
			hessx_ref(1, i2) = hessx_ref(1, i2) + si2 * dtxdxi1 + 
				dtxidxi1;
			hessx_ref(2, i2) = hessx_ref(2, i2) + si2 * dtydxi1 + 
				dtyidxi1;
			hessx_ref(3, i2) = hessx_ref(3, i2) + si2 * dtzdxi1 + 
				dtzidxi1;
			hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * dtxdxi1 - 
				dtxkdxi1;
			hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * dtydxi1 - 
				dtykdxi1;
			hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * dtzdxi1 - 
				dtzkdxi1;
			hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * dtxdxi1 + 
				dtxkdxi1;
			hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * dtydxi1 + 
				dtykdxi1;
			hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * dtzdxi1 + 
				dtzkdxi1;
			hessy_ref(1, i1) = hessy_ref(1, i1) + si1 * dtxdyi1 - 
				dtxidyi1;
			hessy_ref(2, i1) = hessy_ref(2, i1) + si1 * dtydyi1 - 
				dtyidyi1;
			hessy_ref(3, i1) = hessy_ref(3, i1) + si1 * dtzdyi1 - 
				dtzidyi1;
			hessy_ref(1, i2) = hessy_ref(1, i2) + si2 * dtxdyi1 + 
				dtxidyi1;
			hessy_ref(3, i2) = hessy_ref(3, i2) + si2 * dtzdyi1 + 
				dtzidyi1;
			hessy_ref(2, i2) = hessy_ref(2, i2) + si2 * dtydyi1 + 
				dtyidyi1;
			hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * dtxdyi1 - 
				dtxkdyi1;
			hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * dtydyi1 - 
				dtykdyi1;
			hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * dtzdyi1 - 
				dtzkdyi1;
			hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * dtxdyi1 + 
				dtxkdyi1;
			hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * dtydyi1 + 
				dtykdyi1;
			hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * dtzdyi1 + 
				dtzkdyi1;
			hessz_ref(1, i1) = hessz_ref(1, i1) + si1 * dtxdzi1 - 
				dtxidzi1;
			hessz_ref(2, i1) = hessz_ref(2, i1) + si1 * dtydzi1 - 
				dtyidzi1;
			hessz_ref(3, i1) = hessz_ref(3, i1) + si1 * dtzdzi1 - 
				dtzidzi1;
			hessz_ref(1, i2) = hessz_ref(1, i2) + si2 * dtxdzi1 + 
				dtxidzi1;
			hessz_ref(2, i2) = hessz_ref(2, i2) + si2 * dtydzi1 + 
				dtyidzi1;
			hessz_ref(3, i2) = hessz_ref(3, i2) + si2 * dtzdzi1 + 
				dtzidzi1;
			hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * dtxdzi1 - 
				dtxkdzi1;
			hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * dtydzi1 - 
				dtykdzi1;
			hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * dtzdzi1 - 
				dtzkdzi1;
			hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * dtxdzi1 + 
				dtxkdzi1;
			hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * dtydzi1 + 
				dtykdzi1;
			hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * dtzdzi1 + 
				dtzkdzi1;
		    } else if (*i__ == i2) {
			hessx_ref(1, i1) = hessx_ref(1, i1) + si1 * dtxdxi2 - 
				dtxidxi2;
			hessx_ref(2, i1) = hessx_ref(2, i1) + si1 * dtydxi2 - 
				dtyidxi2;
			hessx_ref(3, i1) = hessx_ref(3, i1) + si1 * dtzdxi2 - 
				dtzidxi2;
			hessx_ref(1, i2) = hessx_ref(1, i2) + si2 * dtxdxi2 + 
				dtxidxi2;
			hessx_ref(2, i2) = hessx_ref(2, i2) + si2 * dtydxi2 + 
				dtyidxi2;
			hessx_ref(3, i2) = hessx_ref(3, i2) + si2 * dtzdxi2 + 
				dtzidxi2;
			hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * dtxdxi2 - 
				dtxkdxi2;
			hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * dtydxi2 - 
				dtykdxi2;
			hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * dtzdxi2 - 
				dtzkdxi2;
			hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * dtxdxi2 + 
				dtxkdxi2;
			hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * dtydxi2 + 
				dtykdxi2;
			hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * dtzdxi2 + 
				dtzkdxi2;
			hessy_ref(1, i1) = hessy_ref(1, i1) + si1 * dtxdyi2 - 
				dtxidyi2;
			hessy_ref(2, i1) = hessy_ref(2, i1) + si1 * dtydyi2 - 
				dtyidyi2;
			hessy_ref(3, i1) = hessy_ref(3, i1) + si1 * dtzdyi2 - 
				dtzidyi2;
			hessy_ref(1, i2) = hessy_ref(1, i2) + si2 * dtxdyi2 + 
				dtxidyi2;
			hessy_ref(2, i2) = hessy_ref(2, i2) + si2 * dtydyi2 + 
				dtyidyi2;
			hessy_ref(3, i2) = hessy_ref(3, i2) + si2 * dtzdyi2 + 
				dtzidyi2;
			hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * dtxdyi2 - 
				dtxkdyi2;
			hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * dtydyi2 - 
				dtykdyi2;
			hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * dtzdyi2 - 
				dtzkdyi2;
			hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * dtxdyi2 + 
				dtxkdyi2;
			hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * dtydyi2 + 
				dtykdyi2;
			hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * dtzdyi2 + 
				dtzkdyi2;
			hessz_ref(1, i1) = hessz_ref(1, i1) + si1 * dtxdzi2 - 
				dtxidzi2;
			hessz_ref(2, i1) = hessz_ref(2, i1) + si1 * dtydzi2 - 
				dtyidzi2;
			hessz_ref(3, i1) = hessz_ref(3, i1) + si1 * dtzdzi2 - 
				dtzidzi2;
			hessz_ref(1, i2) = hessz_ref(1, i2) + si2 * dtxdzi2 + 
				dtxidzi2;
			hessz_ref(2, i2) = hessz_ref(2, i2) + si2 * dtydzi2 + 
				dtyidzi2;
			hessz_ref(3, i2) = hessz_ref(3, i2) + si2 * dtzdzi2 + 
				dtzidzi2;
			hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * dtxdzi2 - 
				dtxkdzi2;
			hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * dtydzi2 - 
				dtykdzi2;
			hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * dtzdzi2 - 
				dtzkdzi2;
			hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * dtxdzi2 + 
				dtxkdzi2;
			hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * dtydzi2 + 
				dtykdzi2;
			hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * dtzdzi2 + 
				dtzkdzi2;
		    }

/*     more energy switching if near the cutoff distance */

		    if (r2 > shunt_1.cut2) {
			if (*i__ == i1) {
			    hessx_ref(1, i1) = hessx_ref(1, i1) + si1 * 
				    dtaperx * dedxi1 + si1 * dtaperx * dedxi1 
				    + si1 * si1 * d2taperxx;
			    hessx_ref(2, i1) = hessx_ref(2, i1) + si1 * 
				    dtapery * dedxi1 + si1 * dtaperx * dedyi1 
				    + si1 * si1 * d2taperxy;
			    hessx_ref(3, i1) = hessx_ref(3, i1) + si1 * 
				    dtaperz * dedxi1 + si1 * dtaperx * dedzi1 
				    + si1 * si1 * d2taperxz;
			    hessx_ref(1, i2) = hessx_ref(1, i2) + si2 * 
				    dtaperx * dedxi1 + si1 * dtaperx * dedxi2 
				    + si2 * si1 * d2taperxx;
			    hessx_ref(2, i2) = hessx_ref(2, i2) + si2 * 
				    dtapery * dedxi1 + si1 * dtaperx * dedyi2 
				    + si2 * si1 * d2taperxy;
			    hessx_ref(3, i2) = hessx_ref(3, i2) + si2 * 
				    dtaperz * dedxi1 + si1 * dtaperx * dedzi2 
				    + si2 * si1 * d2taperxz;
			    hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * 
				    dtaperx * dedxi1 + si1 * dtaperx * dedxk1 
				    - sk1 * si1 * d2taperxx;
			    hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * 
				    dtapery * dedxi1 + si1 * dtaperx * dedyk1 
				    - sk1 * si1 * d2taperxy;
			    hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * 
				    dtaperz * dedxi1 + si1 * dtaperx * dedzk1 
				    - sk1 * si1 * d2taperxz;
			    hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * 
				    dtaperx * dedxi1 + si1 * dtaperx * dedxk2 
				    - sk2 * si1 * d2taperxx;
			    hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * 
				    dtapery * dedxi1 + si1 * dtaperx * dedyk2 
				    - sk2 * si1 * d2taperxy;
			    hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * 
				    dtaperz * dedxi1 + si1 * dtaperx * dedzk2 
				    - sk2 * si1 * d2taperxz;
			    hessy_ref(1, i1) = hessy_ref(1, i1) + si1 * 
				    dtaperx * dedyi1 + si1 * dtapery * dedxi1 
				    + si1 * si1 * d2taperxy;
			    hessy_ref(2, i1) = hessy_ref(2, i1) + si1 * 
				    dtapery * dedyi1 + si1 * dtapery * dedyi1 
				    + si1 * si1 * d2taperyy;
			    hessy_ref(3, i1) = hessy_ref(3, i1) + si1 * 
				    dtaperz * dedyi1 + si1 * dtapery * dedzi1 
				    + si1 * si1 * d2taperyz;
			    hessy_ref(1, i2) = hessy_ref(1, i2) + si2 * 
				    dtaperx * dedyi1 + si1 * dtapery * dedxi2 
				    + si2 * si1 * d2taperxy;
			    hessy_ref(2, i2) = hessy_ref(2, i2) + si2 * 
				    dtapery * dedyi1 + si1 * dtapery * dedyi2 
				    + si2 * si1 * d2taperyy;
			    hessy_ref(3, i2) = hessy_ref(3, i2) + si2 * 
				    dtaperz * dedyi1 + si1 * dtapery * dedzi2 
				    + si2 * si1 * d2taperyz;
			    hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * 
				    dtaperx * dedyi1 + si1 * dtapery * dedxk1 
				    - sk1 * si1 * d2taperxy;
			    hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * 
				    dtapery * dedyi1 + si1 * dtapery * dedyk1 
				    - sk1 * si1 * d2taperyy;
			    hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * 
				    dtaperz * dedyi1 + si1 * dtapery * dedzk1 
				    - sk1 * si1 * d2taperyz;
			    hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * 
				    dtaperx * dedyi1 + si1 * dtapery * dedxk2 
				    - sk2 * si1 * d2taperxy;
			    hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * 
				    dtapery * dedyi1 + si1 * dtapery * dedyk2 
				    - sk2 * si1 * d2taperyy;
			    hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * 
				    dtaperz * dedyi1 + si1 * dtapery * dedzk2 
				    - sk2 * si1 * d2taperyz;
			    hessz_ref(1, i1) = hessz_ref(1, i1) + si1 * 
				    dtaperx * dedzi1 + si1 * dtaperz * dedxi1 
				    + si1 * si1 * d2taperxz;
			    hessz_ref(2, i1) = hessz_ref(2, i1) + si1 * 
				    dtapery * dedzi1 + si1 * dtaperz * dedyi1 
				    + si1 * si1 * d2taperyz;
			    hessz_ref(3, i1) = hessz_ref(3, i1) + si1 * 
				    dtaperz * dedzi1 + si1 * dtaperz * dedzi1 
				    + si1 * si1 * d2taperzz;
			    hessz_ref(1, i2) = hessz_ref(1, i2) + si2 * 
				    dtaperx * dedzi1 + si1 * dtaperz * dedxi2 
				    + si2 * si1 * d2taperxz;
			    hessz_ref(2, i2) = hessz_ref(2, i2) + si2 * 
				    dtapery * dedzi1 + si1 * dtaperz * dedyi2 
				    + si2 * si1 * d2taperyz;
			    hessz_ref(3, i2) = hessz_ref(3, i2) + si2 * 
				    dtaperz * dedzi1 + si1 * dtaperz * dedzi2 
				    + si2 * si1 * d2taperzz;
			    hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * 
				    dtaperx * dedzi1 + si1 * dtaperz * dedxk1 
				    - sk1 * si1 * d2taperxz;
			    hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * 
				    dtapery * dedzi1 + si1 * dtaperz * dedyk1 
				    - sk1 * si1 * d2taperyz;
			    hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * 
				    dtaperz * dedzi1 + si1 * dtaperz * dedzk1 
				    - sk1 * si1 * d2taperzz;
			    hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * 
				    dtaperx * dedzi1 + si1 * dtaperz * dedxk2 
				    - sk2 * si1 * d2taperxz;
			    hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * 
				    dtapery * dedzi1 + si1 * dtaperz * dedyk2 
				    - sk2 * si1 * d2taperyz;
			    hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * 
				    dtaperz * dedzi1 + si1 * dtaperz * dedzk2 
				    - sk2 * si1 * d2taperzz;
			} else if (*i__ == i2) {
			    hessx_ref(1, i1) = hessx_ref(1, i1) + si1 * 
				    dtaperx * dedxi2 + si2 * dtaperx * dedxi1 
				    + si1 * si2 * d2taperxx;
			    hessx_ref(2, i1) = hessx_ref(2, i1) + si1 * 
				    dtapery * dedxi2 + si2 * dtaperx * dedyi1 
				    + si1 * si2 * d2taperxy;
			    hessx_ref(3, i1) = hessx_ref(3, i1) + si1 * 
				    dtaperz * dedxi2 + si2 * dtaperx * dedzi1 
				    + si1 * si2 * d2taperxz;
			    hessx_ref(1, i2) = hessx_ref(1, i2) + si2 * 
				    dtaperx * dedxi2 + si2 * dtaperx * dedxi2 
				    + si2 * si2 * d2taperxx;
			    hessx_ref(2, i2) = hessx_ref(2, i2) + si2 * 
				    dtapery * dedxi2 + si2 * dtaperx * dedyi2 
				    + si2 * si2 * d2taperxy;
			    hessx_ref(3, i2) = hessx_ref(3, i2) + si2 * 
				    dtaperz * dedxi2 + si2 * dtaperx * dedzi2 
				    + si2 * si2 * d2taperxz;
			    hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * 
				    dtaperx * dedxi2 + si2 * dtaperx * dedxk1 
				    - sk1 * si2 * d2taperxx;
			    hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * 
				    dtapery * dedxi2 + si2 * dtaperx * dedyk1 
				    - sk1 * si2 * d2taperxy;
			    hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * 
				    dtaperz * dedxi2 + si2 * dtaperx * dedzk1 
				    - sk1 * si2 * d2taperxz;
			    hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * 
				    dtaperx * dedxi2 + si2 * dtaperx * dedxk2 
				    - sk2 * si2 * d2taperxx;
			    hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * 
				    dtapery * dedxi2 + si2 * dtaperx * dedyk2 
				    - sk2 * si2 * d2taperxy;
			    hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * 
				    dtaperz * dedxi2 + si2 * dtaperx * dedzk2 
				    - sk2 * si2 * d2taperxz;
			    hessy_ref(1, i1) = hessy_ref(1, i1) + si1 * 
				    dtaperx * dedyi2 + si2 * dtapery * dedxi1 
				    + si1 * si2 * d2taperxy;
			    hessy_ref(2, i1) = hessy_ref(2, i1) + si1 * 
				    dtapery * dedyi2 + si2 * dtapery * dedyi1 
				    + si1 * si2 * d2taperyy;
			    hessy_ref(3, i1) = hessy_ref(3, i1) + si1 * 
				    dtaperz * dedyi2 + si2 * dtapery * dedzi1 
				    + si1 * si2 * d2taperyz;
			    hessy_ref(1, i2) = hessy_ref(1, i2) + si2 * 
				    dtaperx * dedyi2 + si2 * dtapery * dedxi2 
				    + si2 * si2 * d2taperxy;
			    hessy_ref(2, i2) = hessy_ref(2, i2) + si2 * 
				    dtapery * dedyi2 + si2 * dtapery * dedyi2 
				    + si2 * si2 * d2taperyy;
			    hessy_ref(3, i2) = hessy_ref(3, i2) + si2 * 
				    dtaperz * dedyi2 + si2 * dtapery * dedzi2 
				    + si2 * si2 * d2taperyz;
			    hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * 
				    dtaperx * dedyi2 + si2 * dtapery * dedxk1 
				    - sk1 * si2 * d2taperxy;
			    hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * 
				    dtapery * dedyi2 + si2 * dtapery * dedyk1 
				    - sk1 * si2 * d2taperyy;
			    hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * 
				    dtaperz * dedyi2 + si2 * dtapery * dedzk1 
				    - sk1 * si2 * d2taperyz;
			    hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * 
				    dtaperx * dedyi2 + si2 * dtapery * dedxk2 
				    - sk2 * si2 * d2taperxy;
			    hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * 
				    dtapery * dedyi2 + si2 * dtapery * dedyk2 
				    - sk2 * si2 * d2taperyy;
			    hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * 
				    dtaperz * dedyi2 + si2 * dtapery * dedzk2 
				    - sk2 * si2 * d2taperyz;
			    hessz_ref(1, i1) = hessz_ref(1, i1) + si1 * 
				    dtaperx * dedzi2 + si2 * dtaperz * dedxi1 
				    + si1 * si2 * d2taperxz;
			    hessz_ref(2, i1) = hessz_ref(2, i1) + si1 * 
				    dtapery * dedzi2 + si2 * dtaperz * dedyi1 
				    + si1 * si2 * d2taperyz;
			    hessz_ref(3, i1) = hessz_ref(3, i1) + si1 * 
				    dtaperz * dedzi2 + si2 * dtaperz * dedzi1 
				    + si1 * si2 * d2taperzz;
			    hessz_ref(1, i2) = hessz_ref(1, i2) + si2 * 
				    dtaperx * dedzi2 + si2 * dtaperz * dedxi2 
				    + si2 * si2 * d2taperxz;
			    hessz_ref(2, i2) = hessz_ref(2, i2) + si2 * 
				    dtapery * dedzi2 + si2 * dtaperz * dedyi2 
				    + si2 * si2 * d2taperyz;
			    hessz_ref(3, i2) = hessz_ref(3, i2) + si2 * 
				    dtaperz * dedzi2 + si2 * dtaperz * dedzi2 
				    + si2 * si2 * d2taperzz;
			    hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * 
				    dtaperx * dedzi2 + si2 * dtaperz * dedxk1 
				    - sk1 * si2 * d2taperxz;
			    hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * 
				    dtapery * dedzi2 + si2 * dtaperz * dedyk1 
				    - sk1 * si2 * d2taperyz;
			    hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * 
				    dtaperz * dedzi2 + si2 * dtaperz * dedzk1 
				    - sk1 * si2 * d2taperzz;
			    hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * 
				    dtaperx * dedzi2 + si2 * dtaperz * dedxk2 
				    - sk2 * si2 * d2taperxz;
			    hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * 
				    dtapery * dedzi2 + si2 * dtaperz * dedyk2 
				    - sk2 * si2 * d2taperyz;
			    hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * 
				    dtaperz * dedzi2 + si2 * dtaperz * dedzk2 
				    - sk2 * si2 * d2taperzz;
			}
		    }
		}
	    }
	}
L10:
	;
    }

/*     for periodic boundary conditions with large cutoffs */
/*     neighbors must be found by the replicates method */

    if (! bound_1.use_replica__) {
	return 0;
    }

/*     calculate interaction energy with other unit cells */

    i__1 = dipole_1.ndipole;
    for (idipole = 1; idipole <= i__1; ++idipole) {
	i1 = idpl_ref(1, idipole);
	i2 = idpl_ref(2, idipole);
	si1 = 1. - dipole_1.sdpl[idipole - 1];
	si2 = dipole_1.sdpl[idipole - 1];
	if (i1 != *i__ && i2 != *i__) {
	    goto L30;
	}
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
	fi = f * dipole_1.bdpl[idipole - 1];

/*     decide whether to compute the current interaction */

	i__2 = dipole_1.ndipole;
	for (kdipole = 1; kdipole <= i__2; ++kdipole) {
	    k1 = idpl_ref(1, kdipole);
	    k2 = idpl_ref(2, kdipole);
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i1, &i2, &k1, &k2, &c__0, &c__0);
	    }
	    if (! proceed) {
		goto L20;
	    }

/*     compute the energy contribution for this interaction */

	    sk1 = 1. - dipole_1.sdpl[kdipole - 1];
	    sk2 = dipole_1.sdpl[kdipole - 1];
	    i__3 = cell_1.ncell;
	    for (jcell = 1; jcell <= i__3; ++jcell) {
		xk = atoms_1.x[k2 - 1] - atoms_1.x[k1 - 1];
		yk = atoms_1.y[k2 - 1] - atoms_1.y[k1 - 1];
		zk = atoms_1.z__[k2 - 1] - atoms_1.z__[k1 - 1];
		if (bound_1.use_polymer__) {
		    imager_(&xk, &yk, &zk, &c_n1);
		}
		xr = xq - atoms_1.x[k1 - 1] - xk * sk2;
		yr = yq - atoms_1.y[k1 - 1] - yk * sk2;
		zr = zq - atoms_1.z__[k1 - 1] - zk * sk2;
		imager_(&xr, &yr, &zr, &jcell);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    rk2 = xk * xk + yk * yk + zk * zk;
		    rirkr3 = sqrt(ri2 * rk2 * r2) * r2;
		    dotp = xi * xk + yi * yk + zi * zk;
		    doti = xi * xr + yi * yr + zi * zr;
		    dotk = xk * xr + yk * yr + zk * zr;
		    fik = fi * dipole_1.bdpl[kdipole - 1];
		    if (bound_1.use_polymer__) {
			if (r2 < bound_1.polycut2) {
			    if (k1 == i1 || k1 == i2 || k2 == i1 || k2 == i2) 
				    {
				fik = 0.;
			    }
			}
		    }

/*     some abbreviations used in various chain rule terms */

		    dotik = doti * dotk;
		    enum__ = dotp * r2 - dotik * 3.;
		    r2inv = 15. / r2;
		    ri2inv = 1. / ri2;
		    xrr2 = xr / r2;
		    yrr2 = yr / r2;
		    zrr2 = zr / r2;
		    xiri2 = xi / ri2;
		    yiri2 = yi / ri2;
		    ziri2 = zi / ri2;
		    xkrk2 = xk / rk2;
		    ykrk2 = yk / rk2;
		    zkrk2 = zk / rk2;
		    xixr = xi * xr;
		    xiyr = xi * yr;
		    xizr = xi * zr;
		    yixr = yi * xr;
		    yiyr = yi * yr;
		    yizr = yi * zr;
		    zixr = zi * xr;
		    ziyr = zi * yr;
		    zizr = zi * zr;
		    xkxr = xk * xr;
		    xkyr = xk * yr;
		    xkzr = xk * zr;
		    ykxr = yk * xr;
		    ykyr = yk * yr;
		    ykzr = yk * zr;
		    zkxr = zk * xr;
		    zkyr = zk * yr;
		    zkzr = zk * zr;
		    xixk = xi * xk;
		    xiyk = xi * yk;
		    xizk = xi * zk;
		    yixk = yi * xk;
		    yiyk = yi * yk;
		    yizk = yi * zk;
		    zixk = zi * xk;
		    ziyk = zi * yk;
		    zizk = zi * zk;
		    xrxr = xr * 3. * xr;
		    xryr = xr * 3. * yr;
		    xrzr = xr * 3. * zr;
		    yryr = yr * 3. * yr;
		    yrzr = yr * 3. * zr;
		    zrzr = zr * 3. * zr;
		    xidotk = xi * dotk;
		    yidotk = yi * dotk;
		    zidotk = zi * dotk;
		    xkdoti = xk * doti;
		    ykdoti = yk * doti;
		    zkdoti = zk * doti;

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			fik *= fgrp;
		    }

/*     form the master chain rule term for derivatives */

		    de = -fik / (rirkr3 * r2);

/*     form the chain rule terms for first derivatives */

		    deddotp = -de * r2;
		    deddoti = de * 3. * dotk;
		    deddotk = de * 3. * doti;
		    dedr = de * (dotp * 3. - dotik * 15. / r2);
		    dedrirk = de * enum__;

/*     more first derivative chain rule expressions */

		    termx = dedr * xr + deddoti * xi + deddotk * xk;
		    termy = dedr * yr + deddoti * yi + deddotk * yk;
		    termz = dedr * zr + deddoti * zi + deddotk * zk;
		    termxi = dedrirk * xiri2 + deddotp * xk + deddoti * xr;
		    termyi = dedrirk * yiri2 + deddotp * yk + deddoti * yr;
		    termzi = dedrirk * ziri2 + deddotp * zk + deddoti * zr;
		    termxk = dedrirk * xkrk2 + deddotp * xi + deddotk * xr;
		    termyk = dedrirk * ykrk2 + deddotp * yi + deddotk * yr;
		    termzk = dedrirk * zkrk2 + deddotp * zi + deddotk * zr;

/*     use energy switching if near the cutoff distance */

		    if (r2 > shunt_1.cut2) {
			e = fik * (dotp - doti * 3. * dotk / r2) / rirkr3;
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
			de *= taper;
			termx *= taper;
			termy *= taper;
			termz *= taper;
			termxi *= taper;
			termyi *= taper;
			termzi *= taper;
			termxk *= taper;
			termyk *= taper;
			termzk *= taper;
		    }

/*     chain rule terms for second derivative components */

		    if (*i__ == i1) {
			dtdxi1 = si1 * -5. * xrr2 + xiri2;
			part = si1 * xkdoti - dotk * xr + si1 * xidotk - si1 *
				 2. * dotik * xrr2;
			partik = -xk * r2 + si1 * 2. * dotp * xr - si1 * 3. * 
				xkdoti + xr * 3. * dotk - si1 * 3. * xidotk;
			factor = si1 * 3. * dotp - xkxr * 6. + si1 * 6. * 
				xixk - dotk * 3. - r2inv * (xr * part + si1 * 
				dotik);
			factori = si1 * 3. * dotk + si1 * xkxr + xiri2 * 
				partik - enum__ * (ri2inv - xiri2 * 2. * 
				xiri2);
			factork = r2 + si1 * 3. * doti + si1 * xixr - xrxr + 
				xkrk2 * partik;
			dtxdxi1 = dtdxi1 * termx + de * factor;
			dtxidxi1 = dtdxi1 * termxi + de * factori;
			dtxkdxi1 = dtdxi1 * termxk + de * factork;
			factor = xkyr * -3. - ykxr * 3. + si1 * 3. * xiyk + 
				si1 * 3. * yixk - r2inv * yr * part;
			factori = si1 * -2. * ykxr + si1 * 3. * xkyr + yiri2 *
				 partik + enum__ * 2. * yiri2 * xiri2;
			factork = si1 * -2. * yixr - xryr + si1 * 3. * xiyr + 
				ykrk2 * partik;
			dtydxi1 = dtdxi1 * termy + de * factor;
			dtyidxi1 = dtdxi1 * termyi + de * factori;
			dtykdxi1 = dtdxi1 * termyk + de * factork;
			factor = xkzr * -3. - zkxr * 3. + si1 * 3. * xizk + 
				si1 * 3. * zixk - r2inv * zr * part;
			factori = si1 * -2. * zkxr + si1 * 3. * xkzr + ziri2 *
				 partik + enum__ * 2. * ziri2 * xiri2;
			factork = si1 * -2. * zixr - xrzr + si1 * 3. * xizr + 
				zkrk2 * partik;
			dtzdxi1 = dtdxi1 * termz + de * factor;
			dtzidxi1 = dtdxi1 * termzi + de * factori;
			dtzkdxi1 = dtdxi1 * termzk + de * factork;
			dtdyi1 = si1 * -5. * yrr2 + yiri2;
			part = si1 * ykdoti - dotk * yr + si1 * yidotk - si1 *
				 2. * dotik * yrr2;
			partik = -yk * r2 + si1 * 2. * dotp * yr - si1 * 3. * 
				ykdoti + yr * 3. * dotk - si1 * 3. * yidotk;
			factor = ykxr * -3. - xkyr * 3. + si1 * 3. * yixk + 
				si1 * 3. * xiyk - r2inv * xr * part;
			factori = si1 * -2. * xkyr + si1 * 3. * ykxr + xiri2 *
				 partik + enum__ * 2. * xiri2 * yiri2;
			factork = si1 * -2. * xiyr - xryr + si1 * 3. * yixr + 
				xkrk2 * partik;
			dtxdyi1 = dtdyi1 * termx + de * factor;
			dtxidyi1 = dtdyi1 * termxi + de * factori;
			dtxkdyi1 = dtdyi1 * termxk + de * factork;
			factor = si1 * 3. * dotp - ykyr * 6. + si1 * 6. * 
				yiyk - dotk * 3. - r2inv * (yr * part + si1 * 
				dotik);
			factori = si1 * 3. * dotk + si1 * ykyr + yiri2 * 
				partik - enum__ * (ri2inv - yiri2 * 2. * 
				yiri2);
			factork = r2 + si1 * 3. * doti + si1 * yiyr - yryr + 
				ykrk2 * partik;
			dtydyi1 = dtdyi1 * termy + de * factor;
			dtyidyi1 = dtdyi1 * termyi + de * factori;
			dtykdyi1 = dtdyi1 * termyk + de * factork;
			factor = ykzr * -3. - zkyr * 3. + si1 * 3. * yizk + 
				si1 * 3. * ziyk - r2inv * zr * part;
			factori = si1 * -2. * zkyr + si1 * 3. * ykzr + ziri2 *
				 partik + enum__ * 2. * ziri2 * yiri2;
			factork = si1 * -2. * ziyr - yrzr + si1 * 3. * yizr + 
				zkrk2 * partik;
			dtzdyi1 = dtdyi1 * termz + de * factor;
			dtzidyi1 = dtdyi1 * termzi + de * factori;
			dtzkdyi1 = dtdyi1 * termzk + de * factork;
			dtdzi1 = si1 * -5. * zrr2 + ziri2;
			part = si1 * zkdoti - dotk * zr + si1 * zidotk - si1 *
				 2. * dotik * zrr2;
			partik = -zk * r2 + si1 * 2. * dotp * zr - si1 * 3. * 
				zkdoti + zr * 3. * dotk - si1 * 3. * zidotk;
			factor = zkxr * -3. - xkzr * 3. + si1 * 3. * zixk + 
				si1 * 3. * xizk - r2inv * xr * part;
			factori = si1 * -2. * xkzr + si1 * 3. * zkxr + xiri2 *
				 partik + enum__ * 2. * xiri2 * ziri2;
			factork = si1 * -2. * xizr - xrzr + si1 * 3. * zixr + 
				xkrk2 * partik;
			dtxdzi1 = dtdzi1 * termx + de * factor;
			dtxidzi1 = dtdzi1 * termxi + de * factori;
			dtxkdzi1 = dtdzi1 * termxk + de * factork;
			factor = zkyr * -3. - ykzr * 3. + si1 * 3. * ziyk + 
				si1 * 3. * yizk - r2inv * yr * part;
			factori = si1 * -2. * ykzr + si1 * 3. * zkyr + yiri2 *
				 partik + enum__ * 2. * yiri2 * ziri2;
			factork = si1 * -2. * yizr - yrzr + si1 * 3. * ziyr + 
				ykrk2 * partik;
			dtydzi1 = dtdzi1 * termy + de * factor;
			dtyidzi1 = dtdzi1 * termyi + de * factori;
			dtykdzi1 = dtdzi1 * termyk + de * factork;
			factor = si1 * 3. * dotp - zkzr * 6. + si1 * 6. * 
				zizk - dotk * 3. - r2inv * (zr * part + si1 * 
				dotik);
			factori = si1 * 3. * dotk + si1 * zkzr + ziri2 * 
				partik - enum__ * (ri2inv - ziri2 * 2. * 
				ziri2);
			factork = r2 + si1 * 3. * doti + si1 * zizr - zrzr + 
				zkrk2 * partik;
			dtzdzi1 = dtdzi1 * termz + de * factor;
			dtzidzi1 = dtdzi1 * termzi + de * factori;
			dtzkdzi1 = dtdzi1 * termzk + de * factork;
		    } else if (*i__ == i2) {
			dtdxi2 = si2 * -5. * xrr2 - xiri2;
			part = si2 * xkdoti + dotk * xr + si2 * xidotk - si2 *
				 2. * dotik * xrr2;
			partik = xk * r2 + si2 * 2. * dotp * xr - si2 * 3. * 
				xkdoti - xr * 3. * dotk - si2 * 3. * xidotk;
			factor = si2 * 3. * dotp + xkxr * 6. + si2 * 6. * 
				xixk + dotk * 3. - r2inv * (xr * part + si2 * 
				dotik);
			factori = si2 * 3. * dotk + si2 * xkxr + xiri2 * 
				partik + enum__ * (ri2inv - xiri2 * 2. * 
				xiri2);
			factork = -r2 + si2 * 3. * doti + si2 * xixr + xrxr + 
				xkrk2 * partik;
			dtxdxi2 = dtdxi2 * termx + de * factor;
			dtxidxi2 = dtdxi2 * termxi + de * factori;
			dtxkdxi2 = dtdxi2 * termxk + de * factork;
			factor = xkyr * 3. + ykxr * 3. + si2 * 3. * xiyk + 
				si2 * 3. * yixk - r2inv * yr * part;
			factori = si2 * -2. * ykxr + si2 * 3. * xkyr + yiri2 *
				 partik - enum__ * 2. * yiri2 * xiri2;
			factork = si2 * -2. * yixr + xryr + si2 * 3. * xiyr + 
				ykrk2 * partik;
			dtydxi2 = dtdxi2 * termy + de * factor;
			dtyidxi2 = dtdxi2 * termyi + de * factori;
			dtykdxi2 = dtdxi2 * termyk + de * factork;
			factor = xkzr * 3. + zkxr * 3. + si2 * 3. * xizk + 
				si2 * 3. * zixk - r2inv * zr * part;
			factori = si2 * -2. * zkxr + si2 * 3. * xkzr + ziri2 *
				 partik - enum__ * 2. * ziri2 * xiri2;
			factork = si2 * -2. * zixr + xrzr + si2 * 3. * xizr + 
				zkrk2 * partik;
			dtzdxi2 = dtdxi2 * termz + de * factor;
			dtzidxi2 = dtdxi2 * termzi + de * factori;
			dtzkdxi2 = dtdxi2 * termzk + de * factork;
			dtdyi2 = si2 * -5. * yrr2 - yiri2;
			part = si2 * ykdoti + dotk * yr + si2 * yidotk - si2 *
				 2. * dotik * yrr2;
			partik = yk * r2 + si2 * 2. * dotp * yr - si2 * 3. * 
				ykdoti - yr * 3. * dotk - si2 * 3. * yidotk;
			factor = ykxr * 3. + xkyr * 3. + si2 * 3. * yixk + 
				si2 * 3. * xiyk - r2inv * xr * part;
			factori = si2 * -2. * xkyr + si2 * 3. * ykxr + xiri2 *
				 partik - enum__ * 2. * xiri2 * yiri2;
			factork = si2 * -2. * xiyr + xryr + si2 * 3. * yixr + 
				xkrk2 * partik;
			dtxdyi2 = dtdyi2 * termx + de * factor;
			dtxidyi2 = dtdyi2 * termxi + de * factori;
			dtxkdyi2 = dtdyi2 * termxk + de * factork;
			factor = si2 * 3. * dotp + ykyr * 6. + si2 * 6. * 
				yiyk + dotk * 3. - r2inv * (yr * part + si2 * 
				dotik);
			factori = si2 * 3. * dotk + si2 * ykyr + yiri2 * 
				partik + enum__ * (ri2inv - yiri2 * 2. * 
				yiri2);
			factork = -r2 + si2 * 3. * doti + si2 * yiyr + yryr + 
				ykrk2 * partik;
			dtydyi2 = dtdyi2 * termy + de * factor;
			dtyidyi2 = dtdyi2 * termyi + de * factori;
			dtykdyi2 = dtdyi2 * termyk + de * factork;
			factor = ykzr * 3. + zkyr * 3. + si2 * 3. * yizk + 
				si2 * 3. * ziyk - r2inv * zr * part;
			factori = si2 * -2. * zkyr + si2 * 3. * ykzr + ziri2 *
				 partik - enum__ * 2. * ziri2 * yiri2;
			factork = si2 * -2. * ziyr + yrzr + si2 * 3. * yizr + 
				zkrk2 * partik;
			dtzdyi2 = dtdyi2 * termz + de * factor;
			dtzidyi2 = dtdyi2 * termzi + de * factori;
			dtzkdyi2 = dtdyi2 * termzk + de * factork;
			dtdzi2 = si2 * -5. * zrr2 - ziri2;
			part = si2 * zkdoti + dotk * zr + si2 * zidotk - si2 *
				 2. * dotik * zrr2;
			partik = zk * r2 + si2 * 2. * dotp * zr - si2 * 3. * 
				zkdoti - zr * 3. * dotk - si2 * 3. * zidotk;
			factor = zkxr * 3. + xkzr * 3. + si2 * 3. * zixk + 
				si2 * 3. * xizk - r2inv * xr * part;
			factori = si2 * -2. * xkzr + si2 * 3. * zkxr + xiri2 *
				 partik - enum__ * 2. * xiri2 * ziri2;
			factork = si2 * -2. * xizr + xrzr + si2 * 3. * zixr + 
				xkrk2 * partik;
			dtxdzi2 = dtdzi2 * termx + de * factor;
			dtxidzi2 = dtdzi2 * termxi + de * factori;
			dtxkdzi2 = dtdzi2 * termxk + de * factork;
			factor = zkyr * 3. + ykzr * 3. + si2 * 3. * ziyk + 
				si2 * 3. * yizk - r2inv * yr * part;
			factori = si2 * -2. * ykzr + si2 * 3. * zkyr + yiri2 *
				 partik - enum__ * 2. * yiri2 * ziri2;
			factork = si2 * -2. * yizr + yrzr + si2 * 3. * ziyr + 
				ykrk2 * partik;
			dtydzi2 = dtdzi2 * termy + de * factor;
			dtyidzi2 = dtdzi2 * termyi + de * factori;
			dtykdzi2 = dtdzi2 * termyk + de * factork;
			factor = si2 * 3. * dotp + zkzr * 6. + si2 * 6. * 
				zizk + dotk * 3. - r2inv * (zr * part + si2 * 
				dotik);
			factori = si2 * 3. * dotk + si2 * zkzr + ziri2 * 
				partik + enum__ * (ri2inv - ziri2 * 2. * 
				ziri2);
			factork = -r2 + si2 * 3. * doti + si2 * zizr + zrzr + 
				zkrk2 * partik;
			dtzdzi2 = dtdzi2 * termz + de * factor;
			dtzidzi2 = dtdzi2 * termzi + de * factori;
			dtzkdzi2 = dtdzi2 * termzk + de * factork;
		    }

/*     increment diagonal and off-diagonal Hessian elements */

		    if (*i__ == i1) {
			hessx_ref(1, i1) = hessx_ref(1, i1) + si1 * dtxdxi1 - 
				dtxidxi1;
			hessx_ref(2, i1) = hessx_ref(2, i1) + si1 * dtydxi1 - 
				dtyidxi1;
			hessx_ref(3, i1) = hessx_ref(3, i1) + si1 * dtzdxi1 - 
				dtzidxi1;
			hessx_ref(1, i2) = hessx_ref(1, i2) + si2 * dtxdxi1 + 
				dtxidxi1;
			hessx_ref(2, i2) = hessx_ref(2, i2) + si2 * dtydxi1 + 
				dtyidxi1;
			hessx_ref(3, i2) = hessx_ref(3, i2) + si2 * dtzdxi1 + 
				dtzidxi1;
			hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * dtxdxi1 - 
				dtxkdxi1;
			hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * dtydxi1 - 
				dtykdxi1;
			hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * dtzdxi1 - 
				dtzkdxi1;
			hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * dtxdxi1 + 
				dtxkdxi1;
			hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * dtydxi1 + 
				dtykdxi1;
			hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * dtzdxi1 + 
				dtzkdxi1;
			hessy_ref(1, i1) = hessy_ref(1, i1) + si1 * dtxdyi1 - 
				dtxidyi1;
			hessy_ref(2, i1) = hessy_ref(2, i1) + si1 * dtydyi1 - 
				dtyidyi1;
			hessy_ref(3, i1) = hessy_ref(3, i1) + si1 * dtzdyi1 - 
				dtzidyi1;
			hessy_ref(1, i2) = hessy_ref(1, i2) + si2 * dtxdyi1 + 
				dtxidyi1;
			hessy_ref(3, i2) = hessy_ref(3, i2) + si2 * dtzdyi1 + 
				dtzidyi1;
			hessy_ref(2, i2) = hessy_ref(2, i2) + si2 * dtydyi1 + 
				dtyidyi1;
			hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * dtxdyi1 - 
				dtxkdyi1;
			hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * dtydyi1 - 
				dtykdyi1;
			hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * dtzdyi1 - 
				dtzkdyi1;
			hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * dtxdyi1 + 
				dtxkdyi1;
			hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * dtydyi1 + 
				dtykdyi1;
			hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * dtzdyi1 + 
				dtzkdyi1;
			hessz_ref(1, i1) = hessz_ref(1, i1) + si1 * dtxdzi1 - 
				dtxidzi1;
			hessz_ref(2, i1) = hessz_ref(2, i1) + si1 * dtydzi1 - 
				dtyidzi1;
			hessz_ref(3, i1) = hessz_ref(3, i1) + si1 * dtzdzi1 - 
				dtzidzi1;
			hessz_ref(1, i2) = hessz_ref(1, i2) + si2 * dtxdzi1 + 
				dtxidzi1;
			hessz_ref(2, i2) = hessz_ref(2, i2) + si2 * dtydzi1 + 
				dtyidzi1;
			hessz_ref(3, i2) = hessz_ref(3, i2) + si2 * dtzdzi1 + 
				dtzidzi1;
			hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * dtxdzi1 - 
				dtxkdzi1;
			hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * dtydzi1 - 
				dtykdzi1;
			hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * dtzdzi1 - 
				dtzkdzi1;
			hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * dtxdzi1 + 
				dtxkdzi1;
			hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * dtydzi1 + 
				dtykdzi1;
			hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * dtzdzi1 + 
				dtzkdzi1;
		    } else if (*i__ == i2) {
			hessx_ref(1, i1) = hessx_ref(1, i1) + si1 * dtxdxi2 - 
				dtxidxi2;
			hessx_ref(2, i1) = hessx_ref(2, i1) + si1 * dtydxi2 - 
				dtyidxi2;
			hessx_ref(3, i1) = hessx_ref(3, i1) + si1 * dtzdxi2 - 
				dtzidxi2;
			hessx_ref(1, i2) = hessx_ref(1, i2) + si2 * dtxdxi2 + 
				dtxidxi2;
			hessx_ref(2, i2) = hessx_ref(2, i2) + si2 * dtydxi2 + 
				dtyidxi2;
			hessx_ref(3, i2) = hessx_ref(3, i2) + si2 * dtzdxi2 + 
				dtzidxi2;
			hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * dtxdxi2 - 
				dtxkdxi2;
			hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * dtydxi2 - 
				dtykdxi2;
			hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * dtzdxi2 - 
				dtzkdxi2;
			hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * dtxdxi2 + 
				dtxkdxi2;
			hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * dtydxi2 + 
				dtykdxi2;
			hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * dtzdxi2 + 
				dtzkdxi2;
			hessy_ref(1, i1) = hessy_ref(1, i1) + si1 * dtxdyi2 - 
				dtxidyi2;
			hessy_ref(2, i1) = hessy_ref(2, i1) + si1 * dtydyi2 - 
				dtyidyi2;
			hessy_ref(3, i1) = hessy_ref(3, i1) + si1 * dtzdyi2 - 
				dtzidyi2;
			hessy_ref(1, i2) = hessy_ref(1, i2) + si2 * dtxdyi2 + 
				dtxidyi2;
			hessy_ref(2, i2) = hessy_ref(2, i2) + si2 * dtydyi2 + 
				dtyidyi2;
			hessy_ref(3, i2) = hessy_ref(3, i2) + si2 * dtzdyi2 + 
				dtzidyi2;
			hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * dtxdyi2 - 
				dtxkdyi2;
			hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * dtydyi2 - 
				dtykdyi2;
			hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * dtzdyi2 - 
				dtzkdyi2;
			hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * dtxdyi2 + 
				dtxkdyi2;
			hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * dtydyi2 + 
				dtykdyi2;
			hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * dtzdyi2 + 
				dtzkdyi2;
			hessz_ref(1, i1) = hessz_ref(1, i1) + si1 * dtxdzi2 - 
				dtxidzi2;
			hessz_ref(2, i1) = hessz_ref(2, i1) + si1 * dtydzi2 - 
				dtyidzi2;
			hessz_ref(3, i1) = hessz_ref(3, i1) + si1 * dtzdzi2 - 
				dtzidzi2;
			hessz_ref(1, i2) = hessz_ref(1, i2) + si2 * dtxdzi2 + 
				dtxidzi2;
			hessz_ref(2, i2) = hessz_ref(2, i2) + si2 * dtydzi2 + 
				dtyidzi2;
			hessz_ref(3, i2) = hessz_ref(3, i2) + si2 * dtzdzi2 + 
				dtzidzi2;
			hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * dtxdzi2 - 
				dtxkdzi2;
			hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * dtydzi2 - 
				dtykdzi2;
			hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * dtzdzi2 - 
				dtzkdzi2;
			hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * dtxdzi2 + 
				dtxkdzi2;
			hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * dtydzi2 + 
				dtykdzi2;
			hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * dtzdzi2 + 
				dtzkdzi2;
		    }

/*     more energy switching if near the cutoff distance */

		    if (r2 > shunt_1.cut2) {
			if (*i__ == i1) {
			    hessx_ref(1, i1) = hessx_ref(1, i1) + si1 * 
				    dtaperx * dedxi1 + si1 * dtaperx * dedxi1 
				    + si1 * si1 * d2taperxx;
			    hessx_ref(2, i1) = hessx_ref(2, i1) + si1 * 
				    dtapery * dedxi1 + si1 * dtaperx * dedyi1 
				    + si1 * si1 * d2taperxy;
			    hessx_ref(3, i1) = hessx_ref(3, i1) + si1 * 
				    dtaperz * dedxi1 + si1 * dtaperx * dedzi1 
				    + si1 * si1 * d2taperxz;
			    hessx_ref(1, i2) = hessx_ref(1, i2) + si2 * 
				    dtaperx * dedxi1 + si1 * dtaperx * dedxi2 
				    + si2 * si1 * d2taperxx;
			    hessx_ref(2, i2) = hessx_ref(2, i2) + si2 * 
				    dtapery * dedxi1 + si1 * dtaperx * dedyi2 
				    + si2 * si1 * d2taperxy;
			    hessx_ref(3, i2) = hessx_ref(3, i2) + si2 * 
				    dtaperz * dedxi1 + si1 * dtaperx * dedzi2 
				    + si2 * si1 * d2taperxz;
			    hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * 
				    dtaperx * dedxi1 + si1 * dtaperx * dedxk1 
				    - sk1 * si1 * d2taperxx;
			    hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * 
				    dtapery * dedxi1 + si1 * dtaperx * dedyk1 
				    - sk1 * si1 * d2taperxy;
			    hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * 
				    dtaperz * dedxi1 + si1 * dtaperx * dedzk1 
				    - sk1 * si1 * d2taperxz;
			    hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * 
				    dtaperx * dedxi1 + si1 * dtaperx * dedxk2 
				    - sk2 * si1 * d2taperxx;
			    hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * 
				    dtapery * dedxi1 + si1 * dtaperx * dedyk2 
				    - sk2 * si1 * d2taperxy;
			    hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * 
				    dtaperz * dedxi1 + si1 * dtaperx * dedzk2 
				    - sk2 * si1 * d2taperxz;
			    hessy_ref(1, i1) = hessy_ref(1, i1) + si1 * 
				    dtaperx * dedyi1 + si1 * dtapery * dedxi1 
				    + si1 * si1 * d2taperxy;
			    hessy_ref(2, i1) = hessy_ref(2, i1) + si1 * 
				    dtapery * dedyi1 + si1 * dtapery * dedyi1 
				    + si1 * si1 * d2taperyy;
			    hessy_ref(3, i1) = hessy_ref(3, i1) + si1 * 
				    dtaperz * dedyi1 + si1 * dtapery * dedzi1 
				    + si1 * si1 * d2taperyz;
			    hessy_ref(1, i2) = hessy_ref(1, i2) + si2 * 
				    dtaperx * dedyi1 + si1 * dtapery * dedxi2 
				    + si2 * si1 * d2taperxy;
			    hessy_ref(2, i2) = hessy_ref(2, i2) + si2 * 
				    dtapery * dedyi1 + si1 * dtapery * dedyi2 
				    + si2 * si1 * d2taperyy;
			    hessy_ref(3, i2) = hessy_ref(3, i2) + si2 * 
				    dtaperz * dedyi1 + si1 * dtapery * dedzi2 
				    + si2 * si1 * d2taperyz;
			    hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * 
				    dtaperx * dedyi1 + si1 * dtapery * dedxk1 
				    - sk1 * si1 * d2taperxy;
			    hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * 
				    dtapery * dedyi1 + si1 * dtapery * dedyk1 
				    - sk1 * si1 * d2taperyy;
			    hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * 
				    dtaperz * dedyi1 + si1 * dtapery * dedzk1 
				    - sk1 * si1 * d2taperyz;
			    hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * 
				    dtaperx * dedyi1 + si1 * dtapery * dedxk2 
				    - sk2 * si1 * d2taperxy;
			    hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * 
				    dtapery * dedyi1 + si1 * dtapery * dedyk2 
				    - sk2 * si1 * d2taperyy;
			    hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * 
				    dtaperz * dedyi1 + si1 * dtapery * dedzk2 
				    - sk2 * si1 * d2taperyz;
			    hessz_ref(1, i1) = hessz_ref(1, i1) + si1 * 
				    dtaperx * dedzi1 + si1 * dtaperz * dedxi1 
				    + si1 * si1 * d2taperxz;
			    hessz_ref(2, i1) = hessz_ref(2, i1) + si1 * 
				    dtapery * dedzi1 + si1 * dtaperz * dedyi1 
				    + si1 * si1 * d2taperyz;
			    hessz_ref(3, i1) = hessz_ref(3, i1) + si1 * 
				    dtaperz * dedzi1 + si1 * dtaperz * dedzi1 
				    + si1 * si1 * d2taperzz;
			    hessz_ref(1, i2) = hessz_ref(1, i2) + si2 * 
				    dtaperx * dedzi1 + si1 * dtaperz * dedxi2 
				    + si2 * si1 * d2taperxz;
			    hessz_ref(2, i2) = hessz_ref(2, i2) + si2 * 
				    dtapery * dedzi1 + si1 * dtaperz * dedyi2 
				    + si2 * si1 * d2taperyz;
			    hessz_ref(3, i2) = hessz_ref(3, i2) + si2 * 
				    dtaperz * dedzi1 + si1 * dtaperz * dedzi2 
				    + si2 * si1 * d2taperzz;
			    hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * 
				    dtaperx * dedzi1 + si1 * dtaperz * dedxk1 
				    - sk1 * si1 * d2taperxz;
			    hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * 
				    dtapery * dedzi1 + si1 * dtaperz * dedyk1 
				    - sk1 * si1 * d2taperyz;
			    hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * 
				    dtaperz * dedzi1 + si1 * dtaperz * dedzk1 
				    - sk1 * si1 * d2taperzz;
			    hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * 
				    dtaperx * dedzi1 + si1 * dtaperz * dedxk2 
				    - sk2 * si1 * d2taperxz;
			    hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * 
				    dtapery * dedzi1 + si1 * dtaperz * dedyk2 
				    - sk2 * si1 * d2taperyz;
			    hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * 
				    dtaperz * dedzi1 + si1 * dtaperz * dedzk2 
				    - sk2 * si1 * d2taperzz;
			} else if (*i__ == i2) {
			    hessx_ref(1, i1) = hessx_ref(1, i1) + si1 * 
				    dtaperx * dedxi2 + si2 * dtaperx * dedxi1 
				    + si1 * si2 * d2taperxx;
			    hessx_ref(2, i1) = hessx_ref(2, i1) + si1 * 
				    dtapery * dedxi2 + si2 * dtaperx * dedyi1 
				    + si1 * si2 * d2taperxy;
			    hessx_ref(3, i1) = hessx_ref(3, i1) + si1 * 
				    dtaperz * dedxi2 + si2 * dtaperx * dedzi1 
				    + si1 * si2 * d2taperxz;
			    hessx_ref(1, i2) = hessx_ref(1, i2) + si2 * 
				    dtaperx * dedxi2 + si2 * dtaperx * dedxi2 
				    + si2 * si2 * d2taperxx;
			    hessx_ref(2, i2) = hessx_ref(2, i2) + si2 * 
				    dtapery * dedxi2 + si2 * dtaperx * dedyi2 
				    + si2 * si2 * d2taperxy;
			    hessx_ref(3, i2) = hessx_ref(3, i2) + si2 * 
				    dtaperz * dedxi2 + si2 * dtaperx * dedzi2 
				    + si2 * si2 * d2taperxz;
			    hessx_ref(1, k1) = hessx_ref(1, k1) - sk1 * 
				    dtaperx * dedxi2 + si2 * dtaperx * dedxk1 
				    - sk1 * si2 * d2taperxx;
			    hessx_ref(2, k1) = hessx_ref(2, k1) - sk1 * 
				    dtapery * dedxi2 + si2 * dtaperx * dedyk1 
				    - sk1 * si2 * d2taperxy;
			    hessx_ref(3, k1) = hessx_ref(3, k1) - sk1 * 
				    dtaperz * dedxi2 + si2 * dtaperx * dedzk1 
				    - sk1 * si2 * d2taperxz;
			    hessx_ref(1, k2) = hessx_ref(1, k2) - sk2 * 
				    dtaperx * dedxi2 + si2 * dtaperx * dedxk2 
				    - sk2 * si2 * d2taperxx;
			    hessx_ref(2, k2) = hessx_ref(2, k2) - sk2 * 
				    dtapery * dedxi2 + si2 * dtaperx * dedyk2 
				    - sk2 * si2 * d2taperxy;
			    hessx_ref(3, k2) = hessx_ref(3, k2) - sk2 * 
				    dtaperz * dedxi2 + si2 * dtaperx * dedzk2 
				    - sk2 * si2 * d2taperxz;
			    hessy_ref(1, i1) = hessy_ref(1, i1) + si1 * 
				    dtaperx * dedyi2 + si2 * dtapery * dedxi1 
				    + si1 * si2 * d2taperxy;
			    hessy_ref(2, i1) = hessy_ref(2, i1) + si1 * 
				    dtapery * dedyi2 + si2 * dtapery * dedyi1 
				    + si1 * si2 * d2taperyy;
			    hessy_ref(3, i1) = hessy_ref(3, i1) + si1 * 
				    dtaperz * dedyi2 + si2 * dtapery * dedzi1 
				    + si1 * si2 * d2taperyz;
			    hessy_ref(1, i2) = hessy_ref(1, i2) + si2 * 
				    dtaperx * dedyi2 + si2 * dtapery * dedxi2 
				    + si2 * si2 * d2taperxy;
			    hessy_ref(2, i2) = hessy_ref(2, i2) + si2 * 
				    dtapery * dedyi2 + si2 * dtapery * dedyi2 
				    + si2 * si2 * d2taperyy;
			    hessy_ref(3, i2) = hessy_ref(3, i2) + si2 * 
				    dtaperz * dedyi2 + si2 * dtapery * dedzi2 
				    + si2 * si2 * d2taperyz;
			    hessy_ref(1, k1) = hessy_ref(1, k1) - sk1 * 
				    dtaperx * dedyi2 + si2 * dtapery * dedxk1 
				    - sk1 * si2 * d2taperxy;
			    hessy_ref(2, k1) = hessy_ref(2, k1) - sk1 * 
				    dtapery * dedyi2 + si2 * dtapery * dedyk1 
				    - sk1 * si2 * d2taperyy;
			    hessy_ref(3, k1) = hessy_ref(3, k1) - sk1 * 
				    dtaperz * dedyi2 + si2 * dtapery * dedzk1 
				    - sk1 * si2 * d2taperyz;
			    hessy_ref(1, k2) = hessy_ref(1, k2) - sk2 * 
				    dtaperx * dedyi2 + si2 * dtapery * dedxk2 
				    - sk2 * si2 * d2taperxy;
			    hessy_ref(2, k2) = hessy_ref(2, k2) - sk2 * 
				    dtapery * dedyi2 + si2 * dtapery * dedyk2 
				    - sk2 * si2 * d2taperyy;
			    hessy_ref(3, k2) = hessy_ref(3, k2) - sk2 * 
				    dtaperz * dedyi2 + si2 * dtapery * dedzk2 
				    - sk2 * si2 * d2taperyz;
			    hessz_ref(1, i1) = hessz_ref(1, i1) + si1 * 
				    dtaperx * dedzi2 + si2 * dtaperz * dedxi1 
				    + si1 * si2 * d2taperxz;
			    hessz_ref(2, i1) = hessz_ref(2, i1) + si1 * 
				    dtapery * dedzi2 + si2 * dtaperz * dedyi1 
				    + si1 * si2 * d2taperyz;
			    hessz_ref(3, i1) = hessz_ref(3, i1) + si1 * 
				    dtaperz * dedzi2 + si2 * dtaperz * dedzi1 
				    + si1 * si2 * d2taperzz;
			    hessz_ref(1, i2) = hessz_ref(1, i2) + si2 * 
				    dtaperx * dedzi2 + si2 * dtaperz * dedxi2 
				    + si2 * si2 * d2taperxz;
			    hessz_ref(2, i2) = hessz_ref(2, i2) + si2 * 
				    dtapery * dedzi2 + si2 * dtaperz * dedyi2 
				    + si2 * si2 * d2taperyz;
			    hessz_ref(3, i2) = hessz_ref(3, i2) + si2 * 
				    dtaperz * dedzi2 + si2 * dtaperz * dedzi2 
				    + si2 * si2 * d2taperzz;
			    hessz_ref(1, k1) = hessz_ref(1, k1) - sk1 * 
				    dtaperx * dedzi2 + si2 * dtaperz * dedxk1 
				    - sk1 * si2 * d2taperxz;
			    hessz_ref(2, k1) = hessz_ref(2, k1) - sk1 * 
				    dtapery * dedzi2 + si2 * dtaperz * dedyk1 
				    - sk1 * si2 * d2taperyz;
			    hessz_ref(3, k1) = hessz_ref(3, k1) - sk1 * 
				    dtaperz * dedzi2 + si2 * dtaperz * dedzk1 
				    - sk1 * si2 * d2taperzz;
			    hessz_ref(1, k2) = hessz_ref(1, k2) - sk2 * 
				    dtaperx * dedzi2 + si2 * dtaperz * dedxk2 
				    - sk2 * si2 * d2taperxz;
			    hessz_ref(2, k2) = hessz_ref(2, k2) - sk2 * 
				    dtapery * dedzi2 + si2 * dtaperz * dedyk2 
				    - sk2 * si2 * d2taperyz;
			    hessz_ref(3, k2) = hessz_ref(3, k2) - sk2 * 
				    dtaperz * dedzi2 + si2 * dtaperz * dedzk2 
				    - sk2 * si2 * d2taperzz;
			}
		    }
		}
	    }
L20:
	    ;
	}
L30:
	;
    }
    return 0;
} /* edipole2_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef idpl_ref


