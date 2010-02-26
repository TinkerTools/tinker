/* moments.f -- translated by f2c (version 20050501).
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
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

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
    doublereal bdpl[50000], sdpl[50000];
    integer ndipole, idpl[100000]	/* was [2][50000] */;
} dipole_;

#define dipole_1 dipole_

struct {
    doublereal netchg, netdpl, netqdp[3], xdpl, ydpl, zdpl, xxqdp, xyqdp, 
	    xzqdp, yxqdp, yyqdp, yzqdp, zxqdp, zyqdp, zzqdp;
} moment_;

#define moment_1 moment_

struct {
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

struct {
    doublereal polarity[25000], thole[25000], pdamp[25000], uind[75000]	/* 
	    was [3][25000] */, uinp[75000]	/* was [3][25000] */, uinds[
	    75000]	/* was [3][25000] */, uinps[75000]	/* was [3][
	    25000] */;
    integer npolar;
} polar_;

#define polar_1 polar_

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
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__3 = 3;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine moments  --  total electric multipole moments  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "moments" computes the total electric charge, dipole and */
/*     quadrupole moments for the active atoms as a sum over the */
/*     partial charges, bond dipoles and atomic multipole moments */


/* Subroutine */ int moments_(void)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal a[9]	/* was [3][3] */, b[9]	/* was [3][3] */;
    static integer i__, j, k;
    static doublereal xc, yc, zc, ri, xi, yi, zi, xcm[25000], ycm[25000], zcm[
	    25000], xbnd, qave, ybnd, xmid, ymid, zmid, zbnd;
    extern /* Subroutine */ int born_(void);
    static doublereal work1[3], work2[3], weigh;
    extern /* Subroutine */ int jacobi_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), induce_(
	    void), chkpole_(void), rotpole_(void);


#define a_ref(a_1,a_2) a[(a_2)*3 + a_1 - 4]
#define idpl_ref(a_1,a_2) dipole_1.idpl[(a_2)*2 + a_1 - 3]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define uinds_ref(a_1,a_2) polar_1.uinds[(a_2)*3 + a_1 - 4]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  atmtyp.i  --  atomic properties for each current atom  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     mass      atomic weight for each atom in the system */
/*     tag       integer atom labels from input coordinates file */
/*     class     atom class number for each atom in the system */
/*     atomic    atomic number for each atom in the system */
/*     valence   valence number for each atom in the system */
/*     name      atom name for each atom in the system */
/*     story     descriptive type for each atom in system */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  moment.i  --  components of electric multipole moments  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     netchg   net electric charge for the total system */
/*     netdpl   dipole moment magnitude for the total system */
/*     netqdp   diagonal quadrupole (Qxx, Qyy, Qzz) for system */
/*     xdpl     dipole vector x-component in the global frame */
/*     ydpl     dipole vector y-component in the global frame */
/*     zdpl     dipole vector z-component in the global frame */
/*     xxqdp    quadrupole tensor xx-component in global frame */
/*     xyqdp    quadrupole tensor xy-component in global frame */
/*     xzqdp    quadrupole tensor xz-component in global frame */
/*     yxqdp    quadrupole tensor yx-component in global frame */
/*     yyqdp    quadrupole tensor yy-component in global frame */
/*     yzqdp    quadrupole tensor yz-component in global frame */
/*     zxqdp    quadrupole tensor zx-component in global frame */
/*     zyqdp    quadrupole tensor zy-component in global frame */
/*     zzqdp    quadrupole tensor zz-component in global frame */




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




/*     zero out total charge, dipole and quadrupole components */

    moment_1.netchg = 0.;
    moment_1.netdpl = 0.;
    moment_1.netqdp[0] = 0.;
    moment_1.netqdp[1] = 0.;
    moment_1.netqdp[2] = 0.;
    moment_1.xdpl = 0.;
    moment_1.ydpl = 0.;
    moment_1.zdpl = 0.;
    moment_1.xxqdp = 0.;
    moment_1.xyqdp = 0.;
    moment_1.xzqdp = 0.;
    moment_1.yxqdp = 0.;
    moment_1.yyqdp = 0.;
    moment_1.yzqdp = 0.;
    moment_1.zxqdp = 0.;
    moment_1.zyqdp = 0.;
    moment_1.zzqdp = 0.;

/*     find the center of mass of the set of active atoms */

    weigh = 0.;
    xmid = 0.;
    ymid = 0.;
    zmid = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    weigh += atmtyp_1.mass[i__ - 1];
	    xmid += atoms_1.x[i__ - 1] * atmtyp_1.mass[i__ - 1];
	    ymid += atoms_1.y[i__ - 1] * atmtyp_1.mass[i__ - 1];
	    zmid += atoms_1.z__[i__ - 1] * atmtyp_1.mass[i__ - 1];
	}
    }
    if (weigh != 0.) {
	xmid /= weigh;
	ymid /= weigh;
	zmid /= weigh;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xcm[i__ - 1] = atoms_1.x[i__ - 1] - xmid;
	ycm[i__ - 1] = atoms_1.y[i__ - 1] - ymid;
	zcm[i__ - 1] = atoms_1.z__[i__ - 1] - zmid;
    }

/*     set the multipole moment components due to partial charges */

    i__1 = charge_1.nion;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = charge_1.iion[i__ - 1];
	if (usage_1.use[k - 1]) {
	    moment_1.netchg += charge_1.pchg[i__ - 1];
	    moment_1.xdpl += xcm[k - 1] * charge_1.pchg[i__ - 1];
	    moment_1.ydpl += ycm[k - 1] * charge_1.pchg[i__ - 1];
	    moment_1.zdpl += zcm[k - 1] * charge_1.pchg[i__ - 1];
	    moment_1.xxqdp += xcm[k - 1] * xcm[k - 1] * charge_1.pchg[i__ - 1]
		    ;
	    moment_1.xyqdp += xcm[k - 1] * ycm[k - 1] * charge_1.pchg[i__ - 1]
		    ;
	    moment_1.xzqdp += xcm[k - 1] * zcm[k - 1] * charge_1.pchg[i__ - 1]
		    ;
	    moment_1.yxqdp += ycm[k - 1] * xcm[k - 1] * charge_1.pchg[i__ - 1]
		    ;
	    moment_1.yyqdp += ycm[k - 1] * ycm[k - 1] * charge_1.pchg[i__ - 1]
		    ;
	    moment_1.yzqdp += ycm[k - 1] * zcm[k - 1] * charge_1.pchg[i__ - 1]
		    ;
	    moment_1.zxqdp += zcm[k - 1] * xcm[k - 1] * charge_1.pchg[i__ - 1]
		    ;
	    moment_1.zyqdp += zcm[k - 1] * ycm[k - 1] * charge_1.pchg[i__ - 1]
		    ;
	    moment_1.zzqdp += zcm[k - 1] * zcm[k - 1] * charge_1.pchg[i__ - 1]
		    ;
	}
    }

/*     set the multipole moment components due to bond dipoles */

    i__1 = dipole_1.ndipole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = idpl_ref(1, i__);
	k = idpl_ref(2, i__);
	if (usage_1.use[j - 1] || usage_1.use[k - 1]) {
	    xi = atoms_1.x[j - 1] - atoms_1.x[k - 1];
	    yi = atoms_1.y[j - 1] - atoms_1.y[k - 1];
	    zi = atoms_1.z__[j - 1] - atoms_1.z__[k - 1];
	    ri = sqrt(xi * xi + yi * yi + zi * zi);
	    xbnd = dipole_1.bdpl[i__ - 1] * (xi / ri) / 4.80321;
	    ybnd = dipole_1.bdpl[i__ - 1] * (yi / ri) / 4.80321;
	    zbnd = dipole_1.bdpl[i__ - 1] * (zi / ri) / 4.80321;
	    xc = atoms_1.x[j - 1] - xi * dipole_1.sdpl[i__ - 1];
	    yc = atoms_1.y[j - 1] - yi * dipole_1.sdpl[i__ - 1];
	    zc = atoms_1.z__[j - 1] - zi * dipole_1.sdpl[i__ - 1];
	    moment_1.xdpl += xbnd;
	    moment_1.ydpl += ybnd;
	    moment_1.zdpl += zbnd;
	    moment_1.xxqdp += xc * 2. * xbnd;
	    moment_1.xyqdp = moment_1.xyqdp + xc * ybnd + yc * xbnd;
	    moment_1.xzqdp = moment_1.xzqdp + xc * zbnd + zc * xbnd;
	    moment_1.yxqdp = moment_1.yxqdp + yc * xbnd + xc * ybnd;
	    moment_1.yyqdp += yc * 2. * ybnd;
	    moment_1.yzqdp = moment_1.yzqdp + yc * zbnd + zc * ybnd;
	    moment_1.zxqdp = moment_1.zxqdp + zc * xbnd + xc * zbnd;
	    moment_1.zyqdp = moment_1.zyqdp + zc * ybnd + yc * zbnd;
	    moment_1.zzqdp += zc * 2. * zbnd;
	}
    }

/*     find atomic multipoles and induced dipoles in global frame */

    chkpole_();
    rotpole_();
    if (potent_1.use_born__) {
	born_();
    }
    induce_();
    if (s_cmp(solute_1.solvtyp, "GK", (ftnlen)8, (ftnlen)2) == 0 || s_cmp(
	    solute_1.solvtyp, "PB", (ftnlen)8, (ftnlen)2) == 0) {
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    rpole_ref(2, i__) = rpole_ref(2, i__) + uinds_ref(1, i__);
	    rpole_ref(3, i__) = rpole_ref(3, i__) + uinds_ref(2, i__);
	    rpole_ref(4, i__) = rpole_ref(4, i__) + uinds_ref(3, i__);
	}
    } else {
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    rpole_ref(2, i__) = rpole_ref(2, i__) + uind_ref(1, i__);
	    rpole_ref(3, i__) = rpole_ref(3, i__) + uind_ref(2, i__);
	    rpole_ref(4, i__) = rpole_ref(4, i__) + uind_ref(3, i__);
	}
    }

/*     set the multipole moment components due to atomic multipoles */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = mpole_1.ipole[i__ - 1];
	if (usage_1.use[k - 1]) {
	    moment_1.netchg += rpole_ref(1, i__);
	    moment_1.xdpl = moment_1.xdpl + xcm[k - 1] * rpole_ref(1, i__) + 
		    rpole_ref(2, i__);
	    moment_1.ydpl = moment_1.ydpl + ycm[k - 1] * rpole_ref(1, i__) + 
		    rpole_ref(3, i__);
	    moment_1.zdpl = moment_1.zdpl + zcm[k - 1] * rpole_ref(1, i__) + 
		    rpole_ref(4, i__);
	    moment_1.xxqdp = moment_1.xxqdp + xcm[k - 1] * xcm[k - 1] * 
		    rpole_ref(1, i__) + xcm[k - 1] * 2. * rpole_ref(2, i__);
	    moment_1.xyqdp = moment_1.xyqdp + xcm[k - 1] * ycm[k - 1] * 
		    rpole_ref(1, i__) + xcm[k - 1] * rpole_ref(3, i__) + ycm[
		    k - 1] * rpole_ref(2, i__);
	    moment_1.xzqdp = moment_1.xzqdp + xcm[k - 1] * zcm[k - 1] * 
		    rpole_ref(1, i__) + xcm[k - 1] * rpole_ref(4, i__) + zcm[
		    k - 1] * rpole_ref(2, i__);
	    moment_1.yxqdp = moment_1.yxqdp + ycm[k - 1] * xcm[k - 1] * 
		    rpole_ref(1, i__) + ycm[k - 1] * rpole_ref(2, i__) + xcm[
		    k - 1] * rpole_ref(3, i__);
	    moment_1.yyqdp = moment_1.yyqdp + ycm[k - 1] * ycm[k - 1] * 
		    rpole_ref(1, i__) + ycm[k - 1] * 2. * rpole_ref(3, i__);
	    moment_1.yzqdp = moment_1.yzqdp + ycm[k - 1] * zcm[k - 1] * 
		    rpole_ref(1, i__) + ycm[k - 1] * rpole_ref(4, i__) + zcm[
		    k - 1] * rpole_ref(3, i__);
	    moment_1.zxqdp = moment_1.zxqdp + zcm[k - 1] * xcm[k - 1] * 
		    rpole_ref(1, i__) + zcm[k - 1] * rpole_ref(2, i__) + xcm[
		    k - 1] * rpole_ref(4, i__);
	    moment_1.zyqdp = moment_1.zyqdp + zcm[k - 1] * ycm[k - 1] * 
		    rpole_ref(1, i__) + zcm[k - 1] * rpole_ref(3, i__) + ycm[
		    k - 1] * rpole_ref(4, i__);
	    moment_1.zzqdp = moment_1.zzqdp + zcm[k - 1] * zcm[k - 1] * 
		    rpole_ref(1, i__) + zcm[k - 1] * 2. * rpole_ref(4, i__);
	}
    }

/*     convert the quadrupole from traced to traceless form */

    qave = (moment_1.xxqdp + moment_1.yyqdp + moment_1.zzqdp) / 3.;
    moment_1.xxqdp = (moment_1.xxqdp - qave) * 1.5;
    moment_1.xyqdp *= 1.5;
    moment_1.xzqdp *= 1.5;
    moment_1.yxqdp *= 1.5;
    moment_1.yyqdp = (moment_1.yyqdp - qave) * 1.5;
    moment_1.yzqdp *= 1.5;
    moment_1.zxqdp *= 1.5;
    moment_1.zyqdp *= 1.5;
    moment_1.zzqdp = (moment_1.zzqdp - qave) * 1.5;

/*     add the traceless atomic quadrupoles to total quadrupole */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = mpole_1.ipole[i__ - 1];
	if (usage_1.use[k - 1]) {
	    moment_1.xxqdp += rpole_ref(5, i__) * 3.;
	    moment_1.xyqdp += rpole_ref(6, i__) * 3.;
	    moment_1.xzqdp += rpole_ref(7, i__) * 3.;
	    moment_1.yxqdp += rpole_ref(8, i__) * 3.;
	    moment_1.yyqdp += rpole_ref(9, i__) * 3.;
	    moment_1.yzqdp += rpole_ref(10, i__) * 3.;
	    moment_1.zxqdp += rpole_ref(11, i__) * 3.;
	    moment_1.zyqdp += rpole_ref(12, i__) * 3.;
	    moment_1.zzqdp += rpole_ref(13, i__) * 3.;
	}
    }

/*     convert dipole to Debyes and quadrupole to Buckinghams */

    moment_1.xdpl *= 4.80321;
    moment_1.ydpl *= 4.80321;
    moment_1.zdpl *= 4.80321;
    moment_1.xxqdp *= 4.80321;
    moment_1.xyqdp *= 4.80321;
    moment_1.xzqdp *= 4.80321;
    moment_1.yxqdp *= 4.80321;
    moment_1.yyqdp *= 4.80321;
    moment_1.yzqdp *= 4.80321;
    moment_1.zxqdp *= 4.80321;
    moment_1.zyqdp *= 4.80321;
    moment_1.zzqdp *= 4.80321;

/*     get dipole magnitude and diagonalize quadrupole tensor */

    moment_1.netdpl = sqrt(moment_1.xdpl * moment_1.xdpl + moment_1.ydpl * 
	    moment_1.ydpl + moment_1.zdpl * moment_1.zdpl);
    a_ref(1, 1) = moment_1.xxqdp;
    a_ref(1, 2) = moment_1.xyqdp;
    a_ref(1, 3) = moment_1.xzqdp;
    a_ref(2, 1) = moment_1.yxqdp;
    a_ref(2, 2) = moment_1.yyqdp;
    a_ref(2, 3) = moment_1.yzqdp;
    a_ref(3, 1) = moment_1.zxqdp;
    a_ref(3, 2) = moment_1.zyqdp;
    a_ref(3, 3) = moment_1.zzqdp;
    jacobi_(&c__3, &c__3, a, moment_1.netqdp, b, work1, work2);
    return 0;
} /* moments_ */

#undef uinds_ref
#undef rpole_ref
#undef uind_ref
#undef idpl_ref
#undef a_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine momfull  --  multipole moments for full system  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "momfull"  computes the electric moments for the full system */
/*     as a sum over the partial charges, bond dipoles and atomic */
/*     multipole moments */


/* Subroutine */ int momfull_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static logical temp[25000];
    extern /* Subroutine */ int moments_(void);



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     store active atom list, and make all atoms active */

    if (usage_1.nuse != atoms_1.n) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    temp[i__ - 1] = usage_1.use[i__ - 1];
	    usage_1.use[i__ - 1] = TRUE_;
	}
    }

/*     compute the electric multipole moments for the system */

    moments_();

/*     revert to the original set of active atoms */

    if (usage_1.nuse != atoms_1.n) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    usage_1.use[i__ - 1] = temp[i__ - 1];
	}
    }
    return 0;
} /* momfull_ */

