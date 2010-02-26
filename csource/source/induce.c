/* induce.f -- translated by f2c (version 20050501).
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
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

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
    doublereal rsolv[25000], asolv[25000], rborn[25000], drb[25000], drbp[
	    25000], drobc[25000], doffset, p1, p2, p3, p4, p5, gpol[25000], 
	    shct[25000], aobc[25000], bobc[25000], gobc[25000], vsolv[25000], 
	    wace[1000000]	/* was [1000][1000] */, s2ace[1000000]	/* 
	    was [1000][1000] */, uace[1000000]	/* was [1000][1000] */;
    char solvtyp[8], borntyp[8];
} solute_;

#define solute_1 solute_

struct {
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

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
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

struct {
    doublereal off, off2, cut, cut2, c0, c1, c2, c3, c4, c5, f0, f1, f2, f3, 
	    f4, f5, f6, f7;
} shunt_;

#define shunt_1 shunt_

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

struct {
    doublereal gkr[5000], gkc;
} gk_;

#define gk_1 gk_

struct {
    doublereal pbe, apbe[25000], pbr[25000], pbep[75000]	/* was [3][
	    25000] */, pbfp[75000]	/* was [3][25000] */, pbtp[75000]	
	    /* was [3][25000] */, pbeuind[75000]	/* was [3][25000] */, 
	    pbeuinp[75000]	/* was [3][25000] */, grid[3], gcent[3], 
	    cgrid[3], cgcent[3], fgrid[3], fgcent[3], ionr[10], ionc[10], 
	    pdie, sdie, srad, swin, sdens, smin;
    integer ionn, dime[3], ionq[10];
    char pbtyp[20], pbsoln[20], bcfl[20], srfm[20], chgm[20];
} pb_;

#define pb_1 pb_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;



/*     ############################################################# */
/*     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine induce  --  evaluate induced dipole moments  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "induce" computes the induced dipole moments at polarizable */
/*     sites due to direct or mutual polarization */

/*     assumes multipole components have already been rotated into */
/*     the global coordinate frame; computes induced dipoles based */
/*     on full system, use of active or inactive atoms is ignored */


/* Subroutine */ int induce_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Vacuum Induced Dipole Moments\002,\002 ("
	    "Debyes) :\002)";
    static char fmt_20[] = "(/,\002 Induced Dipole Moments (Debyes) :\002)";
    static char fmt_30[] = "(/,4x,\002Atom\002,14x,\002X\002,15x,\002Y\002,1"
	    "5x,\002Z\002,15x,\002Total\002,/)";
    static char fmt_40[] = "(/,4x,\002Atom\002,14x,\002X\002,13x,\002Y\002,1"
	    "3x,\002Z\002,12x,\002Total\002,/)";
    static char fmt_50[] = "(/,4x,\002Atom\002,14x,\002X\002,11x,\002Y\002,1"
	    "1x,\002Z\002,9x,\002Total\002,/)";
    static char fmt_60[] = "(i8,3x,4f16.8)";
    static char fmt_70[] = "(i8,4x,4f14.6)";
    static char fmt_80[] = "(i8,5x,4f12.4)";
    static char fmt_90[] = "(/,\002 SCRF Induced Dipole Moments\002,\002 (De"
	    "byes) :\002)";
    static char fmt_100[] = "(/,4x,\002Atom\002,14x,\002X\002,15x,\002Y\002,"
	    "15x,\002Z\002,15x,\002Total\002,/)";
    static char fmt_110[] = "(/,4x,\002Atom\002,14x,\002X\002,13x,\002Y\002,"
	    "13x,\002Z\002,12x,\002Total\002,/)";
    static char fmt_120[] = "(/,4x,\002Atom\002,14x,\002X\002,11x,\002Y\002,"
	    "11x,\002Z\002,9x,\002Total\002,/)";
    static char fmt_130[] = "(i8,3x,4f16.8)";
    static char fmt_140[] = "(i8,4x,4f14.6)";
    static char fmt_150[] = "(i8,5x,4f12.4)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void);
    double sqrt(doublereal);
    integer do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k;
    static doublereal norm;
    static logical header;
    extern /* Subroutine */ int induce0a_(void), induce0b_(void), induce0c_(
	    void), induce0d_(void), induce0e_(void), induce0f_(void);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_150, 0 };



#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  inform.i  --  control values for I/O and program flow  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     digits    decimal places output for energy and coordinates */
/*     iprint    steps between status printing (0=no printing) */
/*     iwrite    steps between coordinate dumps (0=no dumps) */
/*     isend     steps between socket communication (0=no sockets) */
/*     verbose   logical flag to turn on extra information */
/*     debug     logical flag to turn on full debug printing */
/*     holdup    logical flag to wait for carriage return on exit */
/*     abort     logical flag to stop execution at next chance */




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




/*     choose the method for summing over multipole field */

    if (s_cmp(solute_1.solvtyp, "PB", (ftnlen)8, (ftnlen)2) == 0) {
	induce0f_();
    } else if (s_cmp(solute_1.solvtyp, "GK", (ftnlen)8, (ftnlen)2) == 0) {
	induce0e_();
    } else if (cutoff_1.use_ewald__) {
	if (cutoff_1.use_mlist__) {
	    induce0d_();
	} else {
	    induce0c_();
	}
    } else {
	if (cutoff_1.use_mlist__) {
	    induce0b_();
	} else {
	    induce0a_();
	}
    }

/*     print out a list of the final induced dipole moments */

    if (inform_1.debug) {
	header = TRUE_;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (polar_1.polarity[i__ - 1] != 0.) {
		if (header) {
		    header = FALSE_;
		    if (s_cmp(solute_1.solvtyp, "GK", (ftnlen)8, (ftnlen)2) ==
			     0 || s_cmp(solute_1.solvtyp, "PB", (ftnlen)8, (
			    ftnlen)2) == 0) {
			io___3.ciunit = iounit_1.iout;
			s_wsfe(&io___3);
			e_wsfe();
		    } else {
			io___4.ciunit = iounit_1.iout;
			s_wsfe(&io___4);
			e_wsfe();
		    }
		    if (inform_1.digits >= 8) {
			io___5.ciunit = iounit_1.iout;
			s_wsfe(&io___5);
			e_wsfe();
		    } else if (inform_1.digits >= 6) {
			io___6.ciunit = iounit_1.iout;
			s_wsfe(&io___6);
			e_wsfe();
		    } else {
			io___7.ciunit = iounit_1.iout;
			s_wsfe(&io___7);
			e_wsfe();
		    }
		}
		k = mpole_1.ipole[i__ - 1];
/* Computing 2nd power */
		d__1 = uind_ref(1, i__);
/* Computing 2nd power */
		d__2 = uind_ref(2, i__);
/* Computing 2nd power */
		d__3 = uind_ref(3, i__);
		norm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
		if (inform_1.digits >= 8) {
		    io___10.ciunit = iounit_1.iout;
		    s_wsfe(&io___10);
		    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		    for (j = 1; j <= 3; ++j) {
			d__1 = uind_ref(j, i__) * 4.80321;
			do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				doublereal));
		    }
		    d__2 = norm * 4.80321;
		    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else if (inform_1.digits >= 6) {
		    io___12.ciunit = iounit_1.iout;
		    s_wsfe(&io___12);
		    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		    for (j = 1; j <= 3; ++j) {
			d__1 = uind_ref(j, i__) * 4.80321;
			do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				doublereal));
		    }
		    d__2 = norm * 4.80321;
		    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else {
		    io___13.ciunit = iounit_1.iout;
		    s_wsfe(&io___13);
		    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		    for (j = 1; j <= 3; ++j) {
			d__1 = uind_ref(j, i__) * 4.80321;
			do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				doublereal));
		    }
		    d__2 = norm * 4.80321;
		    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	}
	header = TRUE_;
	if (s_cmp(solute_1.solvtyp, "GK", (ftnlen)8, (ftnlen)2) == 0 || s_cmp(
		solute_1.solvtyp, "PB", (ftnlen)8, (ftnlen)2) == 0) {
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (polar_1.polarity[i__ - 1] != 0.) {
		    if (header) {
			header = FALSE_;
			io___14.ciunit = iounit_1.iout;
			s_wsfe(&io___14);
			e_wsfe();
			if (inform_1.digits >= 8) {
			    io___15.ciunit = iounit_1.iout;
			    s_wsfe(&io___15);
			    e_wsfe();
			} else if (inform_1.digits >= 6) {
			    io___16.ciunit = iounit_1.iout;
			    s_wsfe(&io___16);
			    e_wsfe();
			} else {
			    io___17.ciunit = iounit_1.iout;
			    s_wsfe(&io___17);
			    e_wsfe();
			}
		    }
		    k = mpole_1.ipole[i__ - 1];
/* Computing 2nd power */
		    d__1 = uinds_ref(1, i__);
/* Computing 2nd power */
		    d__2 = uinds_ref(2, i__);
/* Computing 2nd power */
		    d__3 = uinds_ref(3, i__);
		    norm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
		    if (inform_1.digits >= 8) {
			io___18.ciunit = iounit_1.iout;
			s_wsfe(&io___18);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			for (j = 1; j <= 3; ++j) {
			    d__1 = uinds_ref(j, i__) * 4.80321;
			    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				    doublereal));
			}
			d__2 = norm * 4.80321;
			do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    } else if (inform_1.digits >= 6) {
			io___19.ciunit = iounit_1.iout;
			s_wsfe(&io___19);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			for (j = 1; j <= 3; ++j) {
			    d__1 = uinds_ref(j, i__) * 4.80321;
			    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				    doublereal));
			}
			d__2 = norm * 4.80321;
			do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    } else {
			io___20.ciunit = iounit_1.iout;
			s_wsfe(&io___20);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			for (j = 1; j <= 3; ++j) {
			    d__1 = uinds_ref(j, i__) * 4.80321;
			    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				    doublereal));
			}
			d__2 = norm * 4.80321;
			do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    }
		}
	    }
	}
    }
    return 0;
} /* induce_ */

#undef uinds_ref
#undef uind_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine induce0a  --  induced dipoles via double loop  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "induce0a" computes the induced dipole moments at polarizable */
/*     sites using a pairwise double loop */


/* Subroutine */ int induce0a_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Determination of Induced Dipole\002,\002"
	    " Moments :\002,//,4x,\002Iter\002,8x,\002RMS Change (Debyes)\002"
	    ",/)";
    static char fmt_20[] = "(i8,7x,f16.10)";
    static char fmt_30[] = "(/,\002 INDUCE  --  Warning, Induced Dipole"
	    "s\002,\002 are not Converged\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal r__, r2, ci;
    static integer ii, kk;
    static doublereal ck, xr, yr, zr, rr3, rr5, rr7, fid[3], fkd[3], dir, dix,
	     diy, diz, dkx, dky, dkz, dkr, qir, qkr, eps, pdi, pti, qix, qiy, 
	    qiz, qkx, qky, qkz, fip[3], fkp[3], damp;
    static logical done;
    static doublereal epsd, fgrp;
    static integer iter;
    static doublereal duir, dukr, epsp, duix, duiy, duiz, dukx, duky, puix, 
	    puiy, puiz, dukz, qixx, qixy, qixz, qiyy, qiyz, qizz, pukx, puky, 
	    pukz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, puir, pukr, udir[75000]	
	    /* was [3][25000] */, uold[75000]	/* was [3][25000] */, field[
	    75000]	/* was [3][25000] */;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *), fatal_(void);
    static doublereal udirp[75000]	/* was [3][25000] */, uoldp[75000]	
	    /* was [3][25000] */, scale3, scale5, scale7, dscale[25000], 
	    pgamma, fieldp[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal pscale[25000], epsold;
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *), prterr_(void);
    static logical proceed;
    static doublereal expdamp;
    static integer maxiter;

    /* Fortran I/O blocks */
    static cilist io___110 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___111 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___112 = { 0, 0, 0, fmt_30, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define ip12_ref(a_1,a_2) polgrp_1.ip12[(a_2)*50 + a_1 - 51]
#define ip13_ref(a_1,a_2) polgrp_1.ip13[(a_2)*50 + a_1 - 51]
#define ip14_ref(a_1,a_2) polgrp_1.ip14[(a_2)*50 + a_1 - 51]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define uinp_ref(a_1,a_2) polar_1.uinp[(a_2)*3 + a_1 - 4]
#define udir_ref(a_1,a_2) udir[(a_2)*3 + a_1 - 4]
#define uold_ref(a_1,a_2) uold[(a_2)*3 + a_1 - 4]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define udirp_ref(a_1,a_2) udirp[(a_2)*3 + a_1 - 4]
#define uoldp_ref(a_1,a_2) uoldp[(a_2)*3 + a_1 - 4]
#define fieldp_ref(a_1,a_2) fieldp[(a_2)*3 + a_1 - 4]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  inform.i  --  control values for I/O and program flow  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     digits    decimal places output for energy and coordinates */
/*     iprint    steps between status printing (0=no printing) */
/*     iwrite    steps between coordinate dumps (0=no dumps) */
/*     isend     steps between socket communication (0=no sockets) */
/*     verbose   logical flag to turn on extra information */
/*     debug     logical flag to turn on full debug printing */
/*     holdup    logical flag to wait for carriage return on exit */
/*     abort     logical flag to stop execution at next chance */




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




/*     zero out the induced dipole and the field at each site */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    uind_ref(j, i__) = 0.;
	    uinp_ref(j, i__) = 0.;
	    field_ref(j, i__) = 0.;
	    fieldp_ref(j, i__) = 0.;
	}
    }
    if (! potent_1.use_polar__) {
	return 0;
    }

/*     set the switching function coefficients */

    switch_("MPOLE", (ftnlen)5);

/*     compute the direct induced dipole moment at each atom */

    i__1 = mpole_1.npole - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	pdi = polar_1.pdamp[i__ - 1];
	pti = polar_1.thole[i__ - 1];
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
	i__2 = mpole_1.npole;
	for (j = i__ + 1; j <= i__2; ++j) {
	    dscale[mpole_1.ipole[j - 1] - 1] = 1.;
	    pscale[mpole_1.ipole[j - 1] - 1] = 1.;
	}
	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
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
	    pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = polpot_1.d1scale;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = polpot_1.d2scale;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = polpot_1.d3scale;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = polpot_1.d4scale;
	}
	i__2 = mpole_1.npole;
	for (k = i__ + 1; k <= i__2; ++k) {
	    kk = mpole_1.ipole[k - 1];
	    proceed = TRUE_;
	    if (group_1.use_intra__) {
		groups_(&proceed, &fgrp, &ii, &kk, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
		yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
		zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
		image_(&xr, &yr, &zr);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    r__ = sqrt(r2);
		    ck = rpole_ref(1, k);
		    dkx = rpole_ref(2, k);
		    dky = rpole_ref(3, k);
		    dkz = rpole_ref(4, k);
		    qkxx = rpole_ref(5, k);
		    qkxy = rpole_ref(6, k);
		    qkxz = rpole_ref(7, k);
		    qkyy = rpole_ref(9, k);
		    qkyz = rpole_ref(10, k);
		    qkzz = rpole_ref(13, k);
		    scale3 = 1.;
		    scale5 = 1.;
		    scale7 = 1.;
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
			    scale5 = 1. - expdamp * (1. - damp);
/* Computing 2nd power */
			    d__1 = damp;
			    scale7 = 1. - expdamp * (1. - damp + d__1 * d__1 *
				     .6);
			}
		    }
		    rr3 = scale3 / (r__ * r2);
		    rr5 = scale5 * 3. / (r__ * r2 * r2);
		    rr7 = scale7 * 15. / (r__ * r2 * r2 * r2);
		    dir = dix * xr + diy * yr + diz * zr;
		    qix = qixx * xr + qixy * yr + qixz * zr;
		    qiy = qixy * xr + qiyy * yr + qiyz * zr;
		    qiz = qixz * xr + qiyz * yr + qizz * zr;
		    qir = qix * xr + qiy * yr + qiz * zr;
		    dkr = dkx * xr + dky * yr + dkz * zr;
		    qkx = qkxx * xr + qkxy * yr + qkxz * zr;
		    qky = qkxy * xr + qkyy * yr + qkyz * zr;
		    qkz = qkxz * xr + qkyz * yr + qkzz * zr;
		    qkr = qkx * xr + qky * yr + qkz * zr;
		    fid[0] = -xr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * 
			    dkx + rr5 * 2. * qkx;
		    fid[1] = -yr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * 
			    dky + rr5 * 2. * qky;
		    fid[2] = -zr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * 
			    dkz + rr5 * 2. * qkz;
		    fkd[0] = xr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * 
			    dix - rr5 * 2. * qix;
		    fkd[1] = yr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * 
			    diy - rr5 * 2. * qiy;
		    fkd[2] = zr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * 
			    diz - rr5 * 2. * qiz;
		    for (j = 1; j <= 3; ++j) {
			field_ref(j, i__) = field_ref(j, i__) + fid[j - 1] * 
				dscale[kk - 1];
			field_ref(j, k) = field_ref(j, k) + fkd[j - 1] * 
				dscale[kk - 1];
			fieldp_ref(j, i__) = fieldp_ref(j, i__) + fid[j - 1] *
				 pscale[kk - 1];
			fieldp_ref(j, k) = fieldp_ref(j, k) + fkd[j - 1] * 
				pscale[kk - 1];
		    }
		}
	    }
	}
    }

/*     periodic boundary for large cutoffs via replicates method */

    if (bound_1.use_replica__) {
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    pdi = polar_1.pdamp[i__ - 1];
	    pti = polar_1.thole[i__ - 1];
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
	    i__2 = mpole_1.npole;
	    for (j = i__; j <= i__2; ++j) {
		dscale[mpole_1.ipole[j - 1] - 1] = 1.;
		pscale[mpole_1.ipole[j - 1] - 1] = 1.;
	    }
	    i__2 = couple_1.n12[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	    }
	    i__2 = couple_1.n13[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	    }
	    i__2 = couple_1.n14[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[i14_ref(j, ii) - 1] = polpot_1.p4scale;
	    }
	    i__2 = couple_1.n15[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	    }
	    i__2 = polgrp_1.np11[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip11_ref(j, ii) - 1] = polpot_1.d1scale;
	    }
	    i__2 = polgrp_1.np12[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip12_ref(j, ii) - 1] = polpot_1.d2scale;
	    }
	    i__2 = polgrp_1.np13[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip13_ref(j, ii) - 1] = polpot_1.d3scale;
	    }
	    i__2 = polgrp_1.np14[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip14_ref(j, ii) - 1] = polpot_1.d4scale;
	    }
	    i__2 = mpole_1.npole;
	    for (k = i__; k <= i__2; ++k) {
		kk = mpole_1.ipole[k - 1];
		ck = rpole_ref(1, k);
		dkx = rpole_ref(2, k);
		dky = rpole_ref(3, k);
		dkz = rpole_ref(4, k);
		qkxx = rpole_ref(5, k);
		qkxy = rpole_ref(6, k);
		qkxz = rpole_ref(7, k);
		qkyy = rpole_ref(9, k);
		qkyz = rpole_ref(10, k);
		qkzz = rpole_ref(13, k);
		i__3 = cell_1.ncell;
		for (m = 1; m <= i__3; ++m) {
		    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
		    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
		    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
		    imager_(&xr, &yr, &zr, &m);
		    r2 = xr * xr + yr * yr + zr * zr;
		    if (r2 <= shunt_1.off2) {
			r__ = sqrt(r2);
			scale3 = 1.;
			scale5 = 1.;
			scale7 = 1.;
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
				scale5 = 1. - expdamp * (1. - damp);
/* Computing 2nd power */
				d__1 = damp;
				scale7 = 1. - expdamp * (1. - damp + d__1 * 
					d__1 * .6);
			    }
			}
			rr3 = scale3 / (r__ * r2);
			rr5 = scale5 * 3. / (r__ * r2 * r2);
			rr7 = scale7 * 15. / (r__ * r2 * r2 * r2);
			dir = dix * xr + diy * yr + diz * zr;
			qix = qixx * xr + qixy * yr + qixz * zr;
			qiy = qixy * xr + qiyy * yr + qiyz * zr;
			qiz = qixz * xr + qiyz * yr + qizz * zr;
			qir = qix * xr + qiy * yr + qiz * zr;
			dkr = dkx * xr + dky * yr + dkz * zr;
			qkx = qkxx * xr + qkxy * yr + qkxz * zr;
			qky = qkxy * xr + qkyy * yr + qkyz * zr;
			qkz = qkxz * xr + qkyz * yr + qkzz * zr;
			qkr = qkx * xr + qky * yr + qkz * zr;
			fid[0] = -xr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - 
				rr3 * dkx + rr5 * 2. * qkx;
			fid[1] = -yr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - 
				rr3 * dky + rr5 * 2. * qky;
			fid[2] = -zr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - 
				rr3 * dkz + rr5 * 2. * qkz;
			fkd[0] = xr * (rr3 * ci + rr5 * dir + rr7 * qir) - 
				rr3 * dix - rr5 * 2. * qix;
			fkd[1] = yr * (rr3 * ci + rr5 * dir + rr7 * qir) - 
				rr3 * diy - rr5 * 2. * qiy;
			fkd[2] = zr * (rr3 * ci + rr5 * dir + rr7 * qir) - 
				rr3 * diz - rr5 * 2. * qiz;
			for (j = 1; j <= 3; ++j) {
			    fip[j - 1] = fid[j - 1];
			    fkp[j - 1] = fkd[j - 1];
			}
			if (bound_1.use_polymer__ && r2 <= bound_1.polycut2) {
			    for (j = 1; j <= 3; ++j) {
				fid[j - 1] *= dscale[kk - 1];
				fip[j - 1] *= pscale[kk - 1];
				fkd[j - 1] *= dscale[kk - 1];
				fkp[j - 1] *= pscale[kk - 1];
			    }
			}
			for (j = 1; j <= 3; ++j) {
			    field_ref(j, i__) = field_ref(j, i__) + fid[j - 1]
				    ;
			    fieldp_ref(j, i__) = fieldp_ref(j, i__) + fip[j - 
				    1];
			    if (ii != kk) {
				field_ref(j, k) = field_ref(j, k) + fkd[j - 1]
					;
				fieldp_ref(j, k) = fieldp_ref(j, k) + fkp[j - 
					1];
			    }
			}
		    }
		}
	    }
	}
    }

/*     set induced dipoles to polarizability times direct field */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    udir_ref(j, i__) = polar_1.polarity[i__ - 1] * field_ref(j, i__);
	    udirp_ref(j, i__) = polar_1.polarity[i__ - 1] * fieldp_ref(j, i__)
		    ;
	    uind_ref(j, i__) = udir_ref(j, i__);
	    uinp_ref(j, i__) = udirp_ref(j, i__);
	}
    }

/*     set tolerances for computation of mutual induced dipoles */

    if (s_cmp(polpot_1.poltyp, "MUTUAL", (ftnlen)6, (ftnlen)6) == 0) {
	done = FALSE_;
	maxiter = 500;
	iter = 0;
	eps = 100.;

/*     compute mutual induced dipole moments by an iterative method */

	while(! done) {
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    field_ref(j, i__) = 0.;
		    fieldp_ref(j, i__) = 0.;
		}
	    }
	    i__1 = mpole_1.npole - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ii = mpole_1.ipole[i__ - 1];
		pdi = polar_1.pdamp[i__ - 1];
		pti = polar_1.thole[i__ - 1];
		duix = uind_ref(1, i__);
		duiy = uind_ref(2, i__);
		duiz = uind_ref(3, i__);
		puix = uinp_ref(1, i__);
		puiy = uinp_ref(2, i__);
		puiz = uinp_ref(3, i__);
		i__2 = mpole_1.npole;
		for (j = i__ + 1; j <= i__2; ++j) {
		    dscale[mpole_1.ipole[j - 1] - 1] = 1.;
		}
		i__2 = polgrp_1.np11[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip11_ref(j, ii) - 1] = polpot_1.u1scale;
		}
		i__2 = polgrp_1.np12[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip12_ref(j, ii) - 1] = polpot_1.u2scale;
		}
		i__2 = polgrp_1.np13[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip13_ref(j, ii) - 1] = polpot_1.u3scale;
		}
		i__2 = polgrp_1.np14[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip14_ref(j, ii) - 1] = polpot_1.u4scale;
		}
		i__2 = mpole_1.npole;
		for (k = i__ + 1; k <= i__2; ++k) {
		    kk = mpole_1.ipole[k - 1];
		    proceed = TRUE_;
		    if (group_1.use_intra__) {
			groups_(&proceed, &fgrp, &ii, &kk, &c__0, &c__0, &
				c__0, &c__0);
		    }
		    if (proceed) {
			xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
			yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
			zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
			image_(&xr, &yr, &zr);
			r2 = xr * xr + yr * yr + zr * zr;
			if (r2 <= shunt_1.off2) {
			    r__ = sqrt(r2);
			    dukx = uind_ref(1, k);
			    duky = uind_ref(2, k);
			    dukz = uind_ref(3, k);
			    pukx = uinp_ref(1, k);
			    puky = uinp_ref(2, k);
			    pukz = uinp_ref(3, k);
			    scale3 = dscale[kk - 1];
			    scale5 = dscale[kk - 1];
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
				    scale3 *= 1. - expdamp;
				    scale5 *= 1. - expdamp * (1. - damp);
				}
			    }
			    rr3 = scale3 / (r__ * r2);
			    rr5 = scale5 * 3. / (r__ * r2 * r2);
			    duir = xr * duix + yr * duiy + zr * duiz;
			    dukr = xr * dukx + yr * duky + zr * dukz;
			    puir = xr * puix + yr * puiy + zr * puiz;
			    pukr = xr * pukx + yr * puky + zr * pukz;
			    fid[0] = -rr3 * dukx + rr5 * dukr * xr;
			    fid[1] = -rr3 * duky + rr5 * dukr * yr;
			    fid[2] = -rr3 * dukz + rr5 * dukr * zr;
			    fkd[0] = -rr3 * duix + rr5 * duir * xr;
			    fkd[1] = -rr3 * duiy + rr5 * duir * yr;
			    fkd[2] = -rr3 * duiz + rr5 * duir * zr;
			    fip[0] = -rr3 * pukx + rr5 * pukr * xr;
			    fip[1] = -rr3 * puky + rr5 * pukr * yr;
			    fip[2] = -rr3 * pukz + rr5 * pukr * zr;
			    fkp[0] = -rr3 * puix + rr5 * puir * xr;
			    fkp[1] = -rr3 * puiy + rr5 * puir * yr;
			    fkp[2] = -rr3 * puiz + rr5 * puir * zr;
			    for (j = 1; j <= 3; ++j) {
				field_ref(j, i__) = field_ref(j, i__) + fid[j 
					- 1];
				field_ref(j, k) = field_ref(j, k) + fkd[j - 1]
					;
				fieldp_ref(j, i__) = fieldp_ref(j, i__) + fip[
					j - 1];
				fieldp_ref(j, k) = fieldp_ref(j, k) + fkp[j - 
					1];
			    }
			}
		    }
		}
	    }

/*     periodic boundary for large cutoffs via replicates method */

	    if (bound_1.use_replica__) {
		i__1 = mpole_1.npole;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    ii = mpole_1.ipole[i__ - 1];
		    pdi = polar_1.pdamp[i__ - 1];
		    pti = polar_1.thole[i__ - 1];
		    duix = uind_ref(1, i__);
		    duiy = uind_ref(2, i__);
		    duiz = uind_ref(3, i__);
		    puix = uinp_ref(1, i__);
		    puiy = uinp_ref(2, i__);
		    puiz = uinp_ref(3, i__);
		    i__2 = mpole_1.npole;
		    for (j = i__; j <= i__2; ++j) {
			dscale[mpole_1.ipole[j - 1] - 1] = 1.;
		    }
		    i__2 = polgrp_1.np11[ii - 1];
		    for (j = 1; j <= i__2; ++j) {
			dscale[ip11_ref(j, ii) - 1] = polpot_1.u1scale;
		    }
		    i__2 = polgrp_1.np12[ii - 1];
		    for (j = 1; j <= i__2; ++j) {
			dscale[ip12_ref(j, ii) - 1] = polpot_1.u2scale;
		    }
		    i__2 = polgrp_1.np13[ii - 1];
		    for (j = 1; j <= i__2; ++j) {
			dscale[ip13_ref(j, ii) - 1] = polpot_1.u3scale;
		    }
		    i__2 = polgrp_1.np14[ii - 1];
		    for (j = 1; j <= i__2; ++j) {
			dscale[ip14_ref(j, ii) - 1] = polpot_1.u4scale;
		    }
		    i__2 = mpole_1.npole;
		    for (k = i__; k <= i__2; ++k) {
			kk = mpole_1.ipole[k - 1];
			dukx = uind_ref(1, k);
			duky = uind_ref(2, k);
			dukz = uind_ref(3, k);
			pukx = uinp_ref(1, k);
			puky = uinp_ref(2, k);
			pukz = uinp_ref(3, k);
			proceed = TRUE_;
			i__3 = cell_1.ncell;
			for (m = 1; m <= i__3; ++m) {
			    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
			    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
			    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
			    imager_(&xr, &yr, &zr, &m);
			    r2 = xr * xr + yr * yr + zr * zr;
			    if (r2 <= shunt_1.off2) {
				r__ = sqrt(r2);
				scale3 = 1.;
				scale5 = 1.;
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
					scale5 = 1. - expdamp * (1. - damp);
				    }
				}
				rr3 = scale3 / (r__ * r2);
				rr5 = scale5 * 3. / (r__ * r2 * r2);
				duir = xr * duix + yr * duiy + zr * duiz;
				dukr = xr * dukx + yr * duky + zr * dukz;
				puir = xr * puix + yr * puiy + zr * puiz;
				pukr = xr * pukx + yr * puky + zr * pukz;
				fid[0] = -rr3 * dukx + rr5 * dukr * xr;
				fid[1] = -rr3 * duky + rr5 * dukr * yr;
				fid[2] = -rr3 * dukz + rr5 * dukr * zr;
				fkd[0] = -rr3 * duix + rr5 * duir * xr;
				fkd[1] = -rr3 * duiy + rr5 * duir * yr;
				fkd[2] = -rr3 * duiz + rr5 * duir * zr;
				fip[0] = -rr3 * pukx + rr5 * pukr * xr;
				fip[1] = -rr3 * puky + rr5 * pukr * yr;
				fip[2] = -rr3 * pukz + rr5 * pukr * zr;
				fkp[0] = -rr3 * puix + rr5 * puir * xr;
				fkp[1] = -rr3 * puiy + rr5 * puir * yr;
				fkp[2] = -rr3 * puiz + rr5 * puir * zr;
				if (bound_1.use_polymer__) {
				    if (r2 <= bound_1.polycut2) {
					for (j = 1; j <= 3; ++j) {
					    fid[j - 1] *= dscale[kk - 1];
					    fkd[j - 1] *= dscale[kk - 1];
					    fip[j - 1] *= dscale[kk - 1];
					    fkp[j - 1] *= dscale[kk - 1];
					}
				    }
				}
				for (j = 1; j <= 3; ++j) {
				    field_ref(j, i__) = field_ref(j, i__) + 
					    fid[j - 1];
				    fieldp_ref(j, i__) = fieldp_ref(j, i__) + 
					    fip[j - 1];
				    if (ii != kk) {
					field_ref(j, k) = field_ref(j, k) + 
						fkd[j - 1];
					fieldp_ref(j, k) = fieldp_ref(j, k) + 
						fkp[j - 1];
				    }
				}
			    }
			}
		    }
		}
	    }

/*     check to see if the mutual induced dipoles have converged */

	    ++iter;
	    epsold = eps;
	    epsd = 0.;
	    epsp = 0.;
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    uold_ref(j, i__) = uind_ref(j, i__);
		    uoldp_ref(j, i__) = uinp_ref(j, i__);
		    uind_ref(j, i__) = udir_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * field_ref(j, i__);
		    uinp_ref(j, i__) = udirp_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * fieldp_ref(j, i__);
		    uind_ref(j, i__) = uold_ref(j, i__) + polpot_1.polsor * (
			    uind_ref(j, i__) - uold_ref(j, i__));
		    uinp_ref(j, i__) = uoldp_ref(j, i__) + polpot_1.polsor * (
			    uinp_ref(j, i__) - uoldp_ref(j, i__));
/* Computing 2nd power */
		    d__1 = uind_ref(j, i__) - uold_ref(j, i__);
		    epsd += d__1 * d__1;
/* Computing 2nd power */
		    d__1 = uinp_ref(j, i__) - uoldp_ref(j, i__);
		    epsp += d__1 * d__1;
		}
	    }
	    eps = max(epsd,epsp);
	    eps = sqrt(eps / (doublereal) polar_1.npolar) * 4.80321;
	    if (inform_1.debug) {
		if (iter == 1) {
		    io___110.ciunit = iounit_1.iout;
		    s_wsfe(&io___110);
		    e_wsfe();
		}
		io___111.ciunit = iounit_1.iout;
		s_wsfe(&io___111);
		do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (eps < polpot_1.poleps) {
		done = TRUE_;
	    }
	    if (eps > epsold) {
		done = TRUE_;
	    }
	    if (iter >= maxiter) {
		done = TRUE_;
	    }
	}

/*     terminate the calculation if dipoles failed to converge */

	if (eps > polpot_1.poleps) {
	    io___112.ciunit = iounit_1.iout;
	    s_wsfe(&io___112);
	    e_wsfe();
	    prterr_();
	    fatal_();
	}
    }
    return 0;
} /* induce0a_ */

#undef fieldp_ref
#undef uoldp_ref
#undef udirp_ref
#undef rpole_ref
#undef field_ref
#undef uold_ref
#undef udir_ref
#undef uinp_ref
#undef uind_ref
#undef ip14_ref
#undef ip13_ref
#undef ip12_ref
#undef ip11_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine induce0b  --  induced dipoles via neighbor list  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "induce0b" computes the induced dipole moments at polarizable */
/*     sites using a neighbor list */


/* Subroutine */ int induce0b_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Determination of Induced Dipole\002,\002"
	    " Moments :\002,//,4x,\002Iter\002,8x,\002RMS Change (Debyes)\002"
	    ",/)";
    static char fmt_20[] = "(i8,7x,f16.10)";
    static char fmt_30[] = "(/,\002 INDUCE  --  Warning, Induced Dipole"
	    "s\002,\002 are not Converged\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k;
    static doublereal r__, r2, ci;
    static integer ii, kk;
    static doublereal ck, xr, yr, zr, rr3, rr5, rr7, fid[3], fkd[3];
    static integer kkk;
    static doublereal dix, diy, diz, dkx, dky, dkz, dir, dkr, qir, qkr, eps, 
	    pdi, qix, qiy, qiz, qkx, qky, qkz, pti, fip[3], fkp[3], damp;
    static logical done;
    static doublereal epsd, fgrp;
    static integer iter;
    static doublereal duir, dukr, duix, duiy, duiz, dukx, duky, puix, puiy, 
	    puiz, qixx, qixy, qixz, qiyy, qiyz, qizz, dukz, pukx, puky, pukz, 
	    qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, puir, pukr, epsp, udir[75000]	
	    /* was [3][25000] */, uold[75000]	/* was [3][25000] */, field[
	    75000]	/* was [3][25000] */;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *), fatal_(void);
    static doublereal udirp[75000]	/* was [3][25000] */, uoldp[75000]	
	    /* was [3][25000] */, scale3, scale5, scale7, dscale[25000], 
	    pgamma, fieldp[75000]	/* was [3][25000] */, pscale[25000], 
	    epsold;
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *), prterr_(void);
    static logical proceed;
    static doublereal expdamp;
    static integer maxiter;

    /* Fortran I/O blocks */
    static cilist io___202 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___203 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___204 = { 0, 0, 0, fmt_30, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define ip12_ref(a_1,a_2) polgrp_1.ip12[(a_2)*50 + a_1 - 51]
#define ip13_ref(a_1,a_2) polgrp_1.ip13[(a_2)*50 + a_1 - 51]
#define ip14_ref(a_1,a_2) polgrp_1.ip14[(a_2)*50 + a_1 - 51]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define elst_ref(a_1,a_2) neigh_1.elst[(a_2)*1200 + a_1 - 1201]
#define uinp_ref(a_1,a_2) polar_1.uinp[(a_2)*3 + a_1 - 4]
#define udir_ref(a_1,a_2) udir[(a_2)*3 + a_1 - 4]
#define uold_ref(a_1,a_2) uold[(a_2)*3 + a_1 - 4]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define udirp_ref(a_1,a_2) udirp[(a_2)*3 + a_1 - 4]
#define uoldp_ref(a_1,a_2) uoldp[(a_2)*3 + a_1 - 4]
#define fieldp_ref(a_1,a_2) fieldp[(a_2)*3 + a_1 - 4]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  inform.i  --  control values for I/O and program flow  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     digits    decimal places output for energy and coordinates */
/*     iprint    steps between status printing (0=no printing) */
/*     iwrite    steps between coordinate dumps (0=no dumps) */
/*     isend     steps between socket communication (0=no sockets) */
/*     verbose   logical flag to turn on extra information */
/*     debug     logical flag to turn on full debug printing */
/*     holdup    logical flag to wait for carriage return on exit */
/*     abort     logical flag to stop execution at next chance */




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




/*     zero out the induced dipole and the field at each site */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    uind_ref(j, i__) = 0.;
	    uinp_ref(j, i__) = 0.;
	    field_ref(j, i__) = 0.;
	    fieldp_ref(j, i__) = 0.;
	}
    }
    if (! potent_1.use_polar__) {
	return 0;
    }

/*     set the switching function coefficients */

    switch_("MPOLE", (ftnlen)5);

/*     compute the direct induced dipole moment at each atom */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	pdi = polar_1.pdamp[i__ - 1];
	pti = polar_1.thole[i__ - 1];
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
	i__2 = neigh_1.nelst[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[mpole_1.ipole[elst_ref(j, i__) - 1] - 1] = 1.;
	    pscale[mpole_1.ipole[elst_ref(j, i__) - 1] - 1] = 1.;
	}
	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
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
	    pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = polpot_1.d1scale;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = polpot_1.d2scale;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = polpot_1.d3scale;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = polpot_1.d4scale;
	}
	i__2 = neigh_1.nelst[i__ - 1];
	for (kkk = 1; kkk <= i__2; ++kkk) {
	    k = elst_ref(kkk, i__);
	    kk = mpole_1.ipole[k - 1];
	    proceed = TRUE_;
	    if (group_1.use_intra__) {
		groups_(&proceed, &fgrp, &ii, &kk, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
		yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
		zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
		image_(&xr, &yr, &zr);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    r__ = sqrt(r2);
		    ck = rpole_ref(1, k);
		    dkx = rpole_ref(2, k);
		    dky = rpole_ref(3, k);
		    dkz = rpole_ref(4, k);
		    qkxx = rpole_ref(5, k);
		    qkxy = rpole_ref(6, k);
		    qkxz = rpole_ref(7, k);
		    qkyy = rpole_ref(9, k);
		    qkyz = rpole_ref(10, k);
		    qkzz = rpole_ref(13, k);
		    scale3 = 1.;
		    scale5 = 1.;
		    scale7 = 1.;
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
			    scale5 = 1. - expdamp * (1. - damp);
/* Computing 2nd power */
			    d__1 = damp;
			    scale7 = 1. - expdamp * (1. - damp + d__1 * d__1 *
				     .6);
			}
		    }
		    rr3 = scale3 / (r__ * r2);
		    rr5 = scale5 * 3. / (r__ * r2 * r2);
		    rr7 = scale7 * 15. / (r__ * r2 * r2 * r2);
		    dir = dix * xr + diy * yr + diz * zr;
		    qix = qixx * xr + qixy * yr + qixz * zr;
		    qiy = qixy * xr + qiyy * yr + qiyz * zr;
		    qiz = qixz * xr + qiyz * yr + qizz * zr;
		    qir = qix * xr + qiy * yr + qiz * zr;
		    dkr = dkx * xr + dky * yr + dkz * zr;
		    qkx = qkxx * xr + qkxy * yr + qkxz * zr;
		    qky = qkxy * xr + qkyy * yr + qkyz * zr;
		    qkz = qkxz * xr + qkyz * yr + qkzz * zr;
		    qkr = qkx * xr + qky * yr + qkz * zr;
		    fid[0] = -xr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * 
			    dkx + rr5 * 2. * qkx;
		    fid[1] = -yr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * 
			    dky + rr5 * 2. * qky;
		    fid[2] = -zr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * 
			    dkz + rr5 * 2. * qkz;
		    fkd[0] = xr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * 
			    dix - rr5 * 2. * qix;
		    fkd[1] = yr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * 
			    diy - rr5 * 2. * qiy;
		    fkd[2] = zr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * 
			    diz - rr5 * 2. * qiz;
		    for (j = 1; j <= 3; ++j) {
			field_ref(j, i__) = field_ref(j, i__) + fid[j - 1] * 
				dscale[kk - 1];
			field_ref(j, k) = field_ref(j, k) + fkd[j - 1] * 
				dscale[kk - 1];
			fieldp_ref(j, i__) = fieldp_ref(j, i__) + fid[j - 1] *
				 pscale[kk - 1];
			fieldp_ref(j, k) = fieldp_ref(j, k) + fkd[j - 1] * 
				pscale[kk - 1];
		    }
		}
	    }
	}
    }

/*     set induced dipoles to polarizability times direct field */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    udir_ref(j, i__) = polar_1.polarity[i__ - 1] * field_ref(j, i__);
	    udirp_ref(j, i__) = polar_1.polarity[i__ - 1] * fieldp_ref(j, i__)
		    ;
	    uind_ref(j, i__) = udir_ref(j, i__);
	    uinp_ref(j, i__) = udirp_ref(j, i__);
	}
    }

/*     set tolerances for computation of mutual induced dipoles */

    if (s_cmp(polpot_1.poltyp, "MUTUAL", (ftnlen)6, (ftnlen)6) == 0) {
	done = FALSE_;
	maxiter = 500;
	iter = 0;
	eps = 100.;

/*     compute mutual induced dipole moments by an iterative method */

	while(! done) {
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    field_ref(j, i__) = 0.;
		    fieldp_ref(j, i__) = 0.;
		}
	    }
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ii = mpole_1.ipole[i__ - 1];
		pdi = polar_1.pdamp[i__ - 1];
		pti = polar_1.thole[i__ - 1];
		duix = uind_ref(1, i__);
		duiy = uind_ref(2, i__);
		duiz = uind_ref(3, i__);
		puix = uinp_ref(1, i__);
		puiy = uinp_ref(2, i__);
		puiz = uinp_ref(3, i__);
		i__2 = neigh_1.nelst[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[mpole_1.ipole[elst_ref(j, i__) - 1] - 1] = 1.;
		}
		i__2 = polgrp_1.np11[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip11_ref(j, ii) - 1] = polpot_1.u1scale;
		}
		i__2 = polgrp_1.np12[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip12_ref(j, ii) - 1] = polpot_1.u2scale;
		}
		i__2 = polgrp_1.np13[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip13_ref(j, ii) - 1] = polpot_1.u3scale;
		}
		i__2 = polgrp_1.np14[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip14_ref(j, ii) - 1] = polpot_1.u4scale;
		}
		i__2 = neigh_1.nelst[i__ - 1];
		for (kkk = 1; kkk <= i__2; ++kkk) {
		    k = elst_ref(kkk, i__);
		    kk = mpole_1.ipole[k - 1];
		    proceed = TRUE_;
		    if (group_1.use_intra__) {
			groups_(&proceed, &fgrp, &ii, &kk, &c__0, &c__0, &
				c__0, &c__0);
		    }
		    if (proceed) {
			xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
			yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
			zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
			image_(&xr, &yr, &zr);
			r2 = xr * xr + yr * yr + zr * zr;
			if (r2 <= shunt_1.off2) {
			    r__ = sqrt(r2);
			    dukx = uind_ref(1, k);
			    duky = uind_ref(2, k);
			    dukz = uind_ref(3, k);
			    pukx = uinp_ref(1, k);
			    puky = uinp_ref(2, k);
			    pukz = uinp_ref(3, k);
			    scale3 = dscale[kk - 1];
			    scale5 = dscale[kk - 1];
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
				    scale3 *= 1. - expdamp;
				    scale5 *= 1. - (1. - damp) * expdamp;
				}
			    }
			    rr3 = scale3 / (r__ * r2);
			    rr5 = scale5 * 3. / (r__ * r2 * r2);
			    duir = xr * duix + yr * duiy + zr * duiz;
			    dukr = xr * dukx + yr * duky + zr * dukz;
			    puir = xr * puix + yr * puiy + zr * puiz;
			    pukr = xr * pukx + yr * puky + zr * pukz;
			    fid[0] = -rr3 * dukx + rr5 * dukr * xr;
			    fid[1] = -rr3 * duky + rr5 * dukr * yr;
			    fid[2] = -rr3 * dukz + rr5 * dukr * zr;
			    fkd[0] = -rr3 * duix + rr5 * duir * xr;
			    fkd[1] = -rr3 * duiy + rr5 * duir * yr;
			    fkd[2] = -rr3 * duiz + rr5 * duir * zr;
			    fip[0] = -rr3 * pukx + rr5 * pukr * xr;
			    fip[1] = -rr3 * puky + rr5 * pukr * yr;
			    fip[2] = -rr3 * pukz + rr5 * pukr * zr;
			    fkp[0] = -rr3 * puix + rr5 * puir * xr;
			    fkp[1] = -rr3 * puiy + rr5 * puir * yr;
			    fkp[2] = -rr3 * puiz + rr5 * puir * zr;
			    for (j = 1; j <= 3; ++j) {
				field_ref(j, i__) = field_ref(j, i__) + fid[j 
					- 1];
				field_ref(j, k) = field_ref(j, k) + fkd[j - 1]
					;
				fieldp_ref(j, i__) = fieldp_ref(j, i__) + fip[
					j - 1];
				fieldp_ref(j, k) = fieldp_ref(j, k) + fkp[j - 
					1];
			    }
			}
		    }
		}
	    }

/*     check to see if the mutual induced dipoles have converged */

	    ++iter;
	    epsold = eps;
	    epsd = 0.;
	    epsp = 0.;
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    uold_ref(j, i__) = uind_ref(j, i__);
		    uoldp_ref(j, i__) = uinp_ref(j, i__);
		    uind_ref(j, i__) = udir_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * field_ref(j, i__);
		    uinp_ref(j, i__) = udirp_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * fieldp_ref(j, i__);
		    uind_ref(j, i__) = uold_ref(j, i__) + polpot_1.polsor * (
			    uind_ref(j, i__) - uold_ref(j, i__));
		    uinp_ref(j, i__) = uoldp_ref(j, i__) + polpot_1.polsor * (
			    uinp_ref(j, i__) - uoldp_ref(j, i__));
/* Computing 2nd power */
		    d__1 = uind_ref(j, i__) - uold_ref(j, i__);
		    epsd += d__1 * d__1;
/* Computing 2nd power */
		    d__1 = uinp_ref(j, i__) - uoldp_ref(j, i__);
		    epsp += d__1 * d__1;
		}
	    }
	    eps = max(epsd,epsp);
	    eps = sqrt(eps / (doublereal) polar_1.npolar) * 4.80321;
	    if (inform_1.debug) {
		if (iter == 1) {
		    io___202.ciunit = iounit_1.iout;
		    s_wsfe(&io___202);
		    e_wsfe();
		}
		io___203.ciunit = iounit_1.iout;
		s_wsfe(&io___203);
		do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (eps < polpot_1.poleps) {
		done = TRUE_;
	    }
	    if (eps > epsold) {
		done = TRUE_;
	    }
	    if (iter >= maxiter) {
		done = TRUE_;
	    }
	}

/*     terminate the calculation if dipoles failed to converge */

	if (eps > polpot_1.poleps) {
	    io___204.ciunit = iounit_1.iout;
	    s_wsfe(&io___204);
	    e_wsfe();
	    prterr_();
	    fatal_();
	}
    }
    return 0;
} /* induce0b_ */

#undef fieldp_ref
#undef uoldp_ref
#undef udirp_ref
#undef rpole_ref
#undef field_ref
#undef uold_ref
#undef udir_ref
#undef uinp_ref
#undef elst_ref
#undef uind_ref
#undef ip14_ref
#undef ip13_ref
#undef ip12_ref
#undef ip11_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine induce0c  --  Ewald induced dipoles via loop  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "induce0c" computes the induced dipole moments at polarizable */
/*     sites using particle mesh Ewald summation a double loop */


/* Subroutine */ int induce0c_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Determination of Induced Dipole\002,\002"
	    " Moments :\002,//,4x,\002Iter\002,8x,\002RMS Change (Debyes)\002"
	    ",/)";
    static char fmt_20[] = "(i8,7x,f16.10)";
    static char fmt_30[] = "(/,\002 INDUCE  --  Warning, Induced Dipole"
	    "s\002,\002 are not Converged\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int umutual1_(doublereal *, doublereal *), 
	    udirect2a_(doublereal *, doublereal *);
    static integer i__, j;
    extern /* Subroutine */ int umutual2a_(doublereal *, doublereal *);
    static integer ii;
    static doublereal eps;
    static logical done;
    static doublereal epsd;
    static integer iter;
    static doublereal udir[75000]	/* was [3][25000] */, uold[75000]	
	    /* was [3][25000] */, term, epsp, field[75000]	/* was [3][
	    25000] */;
    extern /* Subroutine */ int fatal_(void);
    static doublereal ucell[3], udirp[75000]	/* was [3][25000] */, uoldp[
	    75000]	/* was [3][25000] */, fieldp[75000]	/* was [3][
	    25000] */, ucellp[3], epsold;
    extern /* Subroutine */ int prterr_(void);
    static integer maxiter;
    extern /* Subroutine */ int udirect1_(doublereal *);

    /* Fortran I/O blocks */
    static cilist io___224 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___225 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___226 = { 0, 0, 0, fmt_30, 0 };



#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define udir_ref(a_1,a_2) udir[(a_2)*3 + a_1 - 4]
#define uold_ref(a_1,a_2) uold[(a_2)*3 + a_1 - 4]
#define uinp_ref(a_1,a_2) polar_1.uinp[(a_2)*3 + a_1 - 4]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define udirp_ref(a_1,a_2) udirp[(a_2)*3 + a_1 - 4]
#define uoldp_ref(a_1,a_2) uoldp[(a_2)*3 + a_1 - 4]
#define fieldp_ref(a_1,a_2) fieldp[(a_2)*3 + a_1 - 4]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  inform.i  --  control values for I/O and program flow  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     digits    decimal places output for energy and coordinates */
/*     iprint    steps between status printing (0=no printing) */
/*     iwrite    steps between coordinate dumps (0=no dumps) */
/*     isend     steps between socket communication (0=no sockets) */
/*     verbose   logical flag to turn on extra information */
/*     debug     logical flag to turn on full debug printing */
/*     holdup    logical flag to wait for carriage return on exit */
/*     abort     logical flag to stop execution at next chance */




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




/*     zero out the induced dipole and the field at each site */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    uind_ref(j, i__) = 0.;
	    uinp_ref(j, i__) = 0.;
	    field_ref(j, i__) = 0.;
	    fieldp_ref(j, i__) = 0.;
	}
    }
    if (! potent_1.use_polar__) {
	return 0;
    }

/*     get the reciprical space part of the electrostatic field */

    udirect1_(field);

/*     get the real space portion of the electrostatic field */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    fieldp_ref(j, i__) = field_ref(j, i__);
	}
    }
    udirect2a_(field, fieldp);

/*     get the self-energy portion of the electrostatic field */

/* Computing 3rd power */
    d__1 = ewald_1.aewald;
    term = d__1 * (d__1 * d__1) * 1.3333333333333333 / 1.772453850905516027;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    field_ref(j, i__) = field_ref(j, i__) + term * rpole_ref(j + 1, 
		    i__);
	    fieldp_ref(j, i__) = fieldp_ref(j, i__) + term * rpole_ref(j + 1, 
		    i__);
	}
    }

/*     compute the cell dipole boundary correction to field */

    if (s_cmp(ewald_1.boundary, "VACUUM", (ftnlen)7, (ftnlen)6) == 0) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    ucell[i__ - 1] = 0.;
	}
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    ucell[0] = ucell[0] + rpole_ref(2, i__) + rpole_ref(1, i__) * 
		    atoms_1.x[ii - 1];
	    ucell[1] = ucell[1] + rpole_ref(3, i__) + rpole_ref(1, i__) * 
		    atoms_1.y[ii - 1];
	    ucell[2] = ucell[2] + rpole_ref(4, i__) + rpole_ref(1, i__) * 
		    atoms_1.z__[ii - 1];
	}
	term = 4.1887902047863905 / boxes_1.volbox;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		field_ref(j, i__) = field_ref(j, i__) - term * ucell[j - 1];
		fieldp_ref(j, i__) = fieldp_ref(j, i__) - term * ucell[j - 1];
	    }
	}
    }

/*     set induced dipoles to polarizability times direct field */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    udir_ref(j, i__) = polar_1.polarity[i__ - 1] * field_ref(j, i__);
	    udirp_ref(j, i__) = polar_1.polarity[i__ - 1] * fieldp_ref(j, i__)
		    ;
	    uind_ref(j, i__) = udir_ref(j, i__);
	    uinp_ref(j, i__) = udirp_ref(j, i__);
	}
    }

/*     set tolerances for computation of mutual induced dipoles */

    if (s_cmp(polpot_1.poltyp, "MUTUAL", (ftnlen)6, (ftnlen)6) == 0) {
	done = FALSE_;
	maxiter = 500;
	iter = 0;
	eps = 100.;

/*     compute mutual induced dipole moments by an iterative method */

	while(! done) {
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    field_ref(j, i__) = 0.;
		    fieldp_ref(j, i__) = 0.;
		}
	    }
	    umutual1_(field, fieldp);
	    umutual2a_(field, fieldp);
/* Computing 3rd power */
	    d__1 = ewald_1.aewald;
	    term = d__1 * (d__1 * d__1) * 1.3333333333333333 / 
		    1.772453850905516027;
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    field_ref(j, i__) = field_ref(j, i__) + term * uind_ref(j,
			     i__);
		    fieldp_ref(j, i__) = fieldp_ref(j, i__) + term * uinp_ref(
			    j, i__);
		}
	    }
	    if (s_cmp(ewald_1.boundary, "VACUUM", (ftnlen)7, (ftnlen)6) == 0) 
		    {
		for (i__ = 1; i__ <= 3; ++i__) {
		    ucell[i__ - 1] = 0.;
		    ucellp[i__ - 1] = 0.;
		}
		i__1 = mpole_1.npole;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    for (j = 1; j <= 3; ++j) {
			ucell[j - 1] += uind_ref(j, i__);
			ucellp[j - 1] += uinp_ref(j, i__);
		    }
		}
		term = 4.1887902047863905 / boxes_1.volbox;
		i__1 = mpole_1.npole;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    for (j = 1; j <= 3; ++j) {
			field_ref(j, i__) = field_ref(j, i__) - term * ucell[
				j - 1];
			fieldp_ref(j, i__) = fieldp_ref(j, i__) - term * 
				ucellp[j - 1];
		    }
		}
	    }

/*     check to see if the mutual induced dipoles have converged */

	    ++iter;
	    epsold = eps;
	    epsd = 0.;
	    epsp = 0.;
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    uold_ref(j, i__) = uind_ref(j, i__);
		    uoldp_ref(j, i__) = uinp_ref(j, i__);
		    uind_ref(j, i__) = udir_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * field_ref(j, i__);
		    uinp_ref(j, i__) = udirp_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * fieldp_ref(j, i__);
		    uind_ref(j, i__) = uold_ref(j, i__) + polpot_1.polsor * (
			    uind_ref(j, i__) - uold_ref(j, i__));
		    uinp_ref(j, i__) = uoldp_ref(j, i__) + polpot_1.polsor * (
			    uinp_ref(j, i__) - uoldp_ref(j, i__));
/* Computing 2nd power */
		    d__1 = uind_ref(j, i__) - uold_ref(j, i__);
		    epsd += d__1 * d__1;
/* Computing 2nd power */
		    d__1 = uinp_ref(j, i__) - uoldp_ref(j, i__);
		    epsp += d__1 * d__1;
		}
	    }
	    eps = max(epsd,epsp);
	    eps = sqrt(eps / (doublereal) polar_1.npolar) * 4.80321;
	    if (inform_1.debug) {
		if (iter == 1) {
		    io___224.ciunit = iounit_1.iout;
		    s_wsfe(&io___224);
		    e_wsfe();
		}
		io___225.ciunit = iounit_1.iout;
		s_wsfe(&io___225);
		do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (eps < polpot_1.poleps) {
		done = TRUE_;
	    }
	    if (eps > epsold) {
		done = TRUE_;
	    }
	    if (iter >= maxiter) {
		done = TRUE_;
	    }
	}

/*     terminate the calculation if dipoles failed to converge */

	if (eps > polpot_1.poleps) {
	    io___226.ciunit = iounit_1.iout;
	    s_wsfe(&io___226);
	    e_wsfe();
	    prterr_();
	    fatal_();
	}
    }
    return 0;
} /* induce0c_ */

#undef fieldp_ref
#undef uoldp_ref
#undef udirp_ref
#undef rpole_ref
#undef field_ref
#undef uinp_ref
#undef uold_ref
#undef udir_ref
#undef uind_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine induce0d  --  Ewald induced dipoles via list  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "induce0d" computes the induced dipole moments at polarizable */
/*     sites using particle mesh Ewald summation and a neighbor list */


/* Subroutine */ int induce0d_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Determination of Induced Dipole\002,\002"
	    " Moments :\002,//,4x,\002Iter\002,8x,\002RMS Change (Debyes)\002"
	    ",/)";
    static char fmt_20[] = "(i8,7x,f16.10)";
    static char fmt_30[] = "(/,\002 INDUCE  --  Warning, Induced Dipole"
	    "s\002,\002 are not Converged\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int umutual1_(doublereal *, doublereal *), 
	    udirect2b_(doublereal *, doublereal *);
    static integer i__, j;
    extern /* Subroutine */ int umutual2b_(doublereal *, doublereal *);
    static integer ii;
    static doublereal eps;
    static logical done;
    static doublereal epsd;
    static integer iter;
    static doublereal udir[75000]	/* was [3][25000] */, uold[75000]	
	    /* was [3][25000] */, term, epsp, field[75000]	/* was [3][
	    25000] */;
    extern /* Subroutine */ int fatal_(void);
    static doublereal ucell[3], udirp[75000]	/* was [3][25000] */, uoldp[
	    75000]	/* was [3][25000] */, fieldp[75000]	/* was [3][
	    25000] */, ucellp[3], epsold;
    extern /* Subroutine */ int prterr_(void);
    static integer maxiter;
    extern /* Subroutine */ int udirect1_(doublereal *);

    /* Fortran I/O blocks */
    static cilist io___246 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___247 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___248 = { 0, 0, 0, fmt_30, 0 };



#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define udir_ref(a_1,a_2) udir[(a_2)*3 + a_1 - 4]
#define uold_ref(a_1,a_2) uold[(a_2)*3 + a_1 - 4]
#define uinp_ref(a_1,a_2) polar_1.uinp[(a_2)*3 + a_1 - 4]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define udirp_ref(a_1,a_2) udirp[(a_2)*3 + a_1 - 4]
#define uoldp_ref(a_1,a_2) uoldp[(a_2)*3 + a_1 - 4]
#define fieldp_ref(a_1,a_2) fieldp[(a_2)*3 + a_1 - 4]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  inform.i  --  control values for I/O and program flow  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     digits    decimal places output for energy and coordinates */
/*     iprint    steps between status printing (0=no printing) */
/*     iwrite    steps between coordinate dumps (0=no dumps) */
/*     isend     steps between socket communication (0=no sockets) */
/*     verbose   logical flag to turn on extra information */
/*     debug     logical flag to turn on full debug printing */
/*     holdup    logical flag to wait for carriage return on exit */
/*     abort     logical flag to stop execution at next chance */




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




/*     zero out the induced dipole and the field at each site */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    uind_ref(j, i__) = 0.;
	    uinp_ref(j, i__) = 0.;
	    field_ref(j, i__) = 0.;
	    fieldp_ref(j, i__) = 0.;
	}
    }
    if (! potent_1.use_polar__) {
	return 0;
    }

/*     get the reciprical space part of the electrostatic field */

    udirect1_(field);

/*     get the real space portion of the electrostatic field */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    fieldp_ref(j, i__) = field_ref(j, i__);
	}
    }
    udirect2b_(field, fieldp);

/*     get the self-energy portion of the electrostatic field */

/* Computing 3rd power */
    d__1 = ewald_1.aewald;
    term = d__1 * (d__1 * d__1) * 1.3333333333333333 / 1.772453850905516027;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    field_ref(j, i__) = field_ref(j, i__) + term * rpole_ref(j + 1, 
		    i__);
	    fieldp_ref(j, i__) = fieldp_ref(j, i__) + term * rpole_ref(j + 1, 
		    i__);
	}
    }

/*     compute the cell dipole boundary correction to field */

    if (s_cmp(ewald_1.boundary, "VACUUM", (ftnlen)7, (ftnlen)6) == 0) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    ucell[i__ - 1] = 0.;
	}
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    ucell[0] = ucell[0] + rpole_ref(2, i__) + rpole_ref(1, i__) * 
		    atoms_1.x[ii - 1];
	    ucell[1] = ucell[1] + rpole_ref(3, i__) + rpole_ref(1, i__) * 
		    atoms_1.y[ii - 1];
	    ucell[2] = ucell[2] + rpole_ref(4, i__) + rpole_ref(1, i__) * 
		    atoms_1.z__[ii - 1];
	}
	term = 4.1887902047863905 / boxes_1.volbox;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		field_ref(j, i__) = field_ref(j, i__) - term * ucell[j - 1];
		fieldp_ref(j, i__) = fieldp_ref(j, i__) - term * ucell[j - 1];
	    }
	}
    }

/*     set induced dipoles to polarizability times direct field */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    udir_ref(j, i__) = polar_1.polarity[i__ - 1] * field_ref(j, i__);
	    udirp_ref(j, i__) = polar_1.polarity[i__ - 1] * fieldp_ref(j, i__)
		    ;
	    uind_ref(j, i__) = udir_ref(j, i__);
	    uinp_ref(j, i__) = udirp_ref(j, i__);
	}
    }

/*     set tolerances for computation of mutual induced dipoles */

    if (s_cmp(polpot_1.poltyp, "MUTUAL", (ftnlen)6, (ftnlen)6) == 0) {
	done = FALSE_;
	maxiter = 500;
	iter = 0;
	eps = 100.;

/*     compute mutual induced dipole moments by an iterative method */

	while(! done) {
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    field_ref(j, i__) = 0.;
		    fieldp_ref(j, i__) = 0.;
		}
	    }
	    umutual1_(field, fieldp);
	    umutual2b_(field, fieldp);
/* Computing 3rd power */
	    d__1 = ewald_1.aewald;
	    term = d__1 * (d__1 * d__1) * 1.3333333333333333 / 
		    1.772453850905516027;
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    field_ref(j, i__) = field_ref(j, i__) + term * uind_ref(j,
			     i__);
		    fieldp_ref(j, i__) = fieldp_ref(j, i__) + term * uinp_ref(
			    j, i__);
		}
	    }
	    if (s_cmp(ewald_1.boundary, "VACUUM", (ftnlen)7, (ftnlen)6) == 0) 
		    {
		for (i__ = 1; i__ <= 3; ++i__) {
		    ucell[i__ - 1] = 0.;
		    ucellp[i__ - 1] = 0.;
		}
		i__1 = mpole_1.npole;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    for (j = 1; j <= 3; ++j) {
			ucell[j - 1] += uind_ref(j, i__);
			ucellp[j - 1] += uinp_ref(j, i__);
		    }
		}
		term = 4.1887902047863905 / boxes_1.volbox;
		i__1 = mpole_1.npole;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    for (j = 1; j <= 3; ++j) {
			field_ref(j, i__) = field_ref(j, i__) - term * ucell[
				j - 1];
			fieldp_ref(j, i__) = fieldp_ref(j, i__) - term * 
				ucellp[j - 1];
		    }
		}
	    }

/*     check to see if the mutual induced dipoles have converged */

	    ++iter;
	    epsold = eps;
	    epsd = 0.;
	    epsp = 0.;
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    uold_ref(j, i__) = uind_ref(j, i__);
		    uoldp_ref(j, i__) = uinp_ref(j, i__);
		    uind_ref(j, i__) = udir_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * field_ref(j, i__);
		    uinp_ref(j, i__) = udirp_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * fieldp_ref(j, i__);
		    uind_ref(j, i__) = uold_ref(j, i__) + polpot_1.polsor * (
			    uind_ref(j, i__) - uold_ref(j, i__));
		    uinp_ref(j, i__) = uoldp_ref(j, i__) + polpot_1.polsor * (
			    uinp_ref(j, i__) - uoldp_ref(j, i__));
/* Computing 2nd power */
		    d__1 = uind_ref(j, i__) - uold_ref(j, i__);
		    epsd += d__1 * d__1;
/* Computing 2nd power */
		    d__1 = uinp_ref(j, i__) - uoldp_ref(j, i__);
		    epsp += d__1 * d__1;
		}
	    }
	    eps = max(epsd,epsp);
	    eps = sqrt(eps / (doublereal) polar_1.npolar) * 4.80321;
	    if (inform_1.debug) {
		if (iter == 1) {
		    io___246.ciunit = iounit_1.iout;
		    s_wsfe(&io___246);
		    e_wsfe();
		}
		io___247.ciunit = iounit_1.iout;
		s_wsfe(&io___247);
		do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (eps < polpot_1.poleps) {
		done = TRUE_;
	    }
	    if (eps > epsold) {
		done = TRUE_;
	    }
	    if (iter >= maxiter) {
		done = TRUE_;
	    }
	}

/*     terminate the calculation if dipoles failed to converge */

	if (eps > polpot_1.poleps) {
	    io___248.ciunit = iounit_1.iout;
	    s_wsfe(&io___248);
	    e_wsfe();
	    prterr_();
	    fatal_();
	}
    }
    return 0;
} /* induce0d_ */

#undef fieldp_ref
#undef uoldp_ref
#undef udirp_ref
#undef rpole_ref
#undef field_ref
#undef uinp_ref
#undef uold_ref
#undef udir_ref
#undef uind_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine udirect1  --  Ewald recip direct induced field  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "udirect1" computes the reciprocal space contribution of the */
/*     permanent atomic multipole moments to the field */


/* Subroutine */ int udirect1_(doublereal *field)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal), cos(doublereal);

    /* Local variables */
    extern /* Subroutine */ int fftfront_(void);
    static integer i__, j, k;
    static doublereal h1, h2, h3;
    static integer k1, k2, k3, m1, m2, m3;
    static doublereal r1, r2, r3;
    extern /* Subroutine */ int grid_mpole__(doublereal *), fphi_mpole__(
	    doublereal *), cmp_to_fmp__(doublereal *, doublereal *);
    static integer nf1, nf2, nf3, nff;
    static doublereal cmp[250000]	/* was [10][25000] */, fmp[250000]	
	    /* was [10][25000] */, hsq, cphi[250000]	/* was [10][25000] */,
	     fphi[500000]	/* was [20][25000] */, term;
    static integer ntot;
    extern /* Subroutine */ int fphi_to_cphi__(doublereal *, doublereal *), 
	    bspline_fill__(void);
    static doublereal denom, pterm;
    extern /* Subroutine */ int fftback_(void);
    static doublereal expterm, volterm;


#define cmp_ref(a_1,a_2) cmp[(a_2)*10 + a_1 - 11]
#define qfac_ref(a_1,a_2,a_3) pme_1.qfac[((a_3)*100 + (a_2))*100 + a_1 - \
10101]
#define cphi_ref(a_1,a_2) cphi[(a_2)*10 + a_1 - 11]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1]
#define recip_ref(a_1,a_2) boxes_1.recip[(a_2)*3 + a_1 - 4]
#define qgrid_ref(a_1,a_2,a_3,a_4) pme_1.qgrid[(((a_4)*100 + (a_3))*100 + \
(a_2))*2 + a_1 - 20203]
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




/*     return if the Ewald coefficient is zero */

    /* Parameter adjustments */
    field -= 4;

    /* Function Body */
    if (ewald_1.aewald < 1e-6) {
	return 0;
    }

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

/*     compute the arrays of B-spline coefficients */

    bspline_fill__();

/*     convert Cartesian multipoles to fractional coordinates */

    cmp_to_fmp__(cmp, fmp);

/*     assign PME grid and perform 3-D FFT forward transform */

    grid_mpole__(fmp);
    fftfront_();

/*     make the scalar summation over reciprocal lattice */

    qfac_ref(1, 1, 1) = 0.;
/* Computing 2nd power */
    d__1 = 3.141592653589793238 / ewald_1.aewald;
    pterm = d__1 * d__1;
    volterm = boxes_1.volbox * 3.141592653589793238;
    nff = pme_1.nfft1 * pme_1.nfft2;
    nf1 = (pme_1.nfft1 + 1) / 2;
    nf2 = (pme_1.nfft2 + 1) / 2;
    nf3 = (pme_1.nfft3 + 1) / 2;
    ntot = pme_1.nfft1 * pme_1.nfft2 * pme_1.nfft3;
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
	}
	qfac_ref(k1, k2, k3) = expterm;
    }

/*     account for the zeroth grid point for a finite system */

    qfac_ref(1, 1, 1) = 0.;
    if (! bound_1.use_bounds__) {
	expterm = 1.5707963267948966 / boxes_1.xbox;
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

/*     perform 3-D FFT backward transform and get field */

    fftback_();
    fphi_mpole__(fphi);

/*     convert the field from fractional to Cartesian */

    fphi_to_cphi__(fphi, cphi);

/*     increment the field at each multipole site */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	field_ref(1, i__) = field_ref(1, i__) - cphi_ref(2, i__);
	field_ref(2, i__) = field_ref(2, i__) - cphi_ref(3, i__);
	field_ref(3, i__) = field_ref(3, i__) - cphi_ref(4, i__);
    }
    return 0;
} /* udirect1_ */

#undef rpole_ref
#undef qgrid_ref
#undef recip_ref
#undef field_ref
#undef cphi_ref
#undef qfac_ref
#undef cmp_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine udirect2a  --  Ewald real direct field via loop  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "udirect2a" computes the real space contribution of the permanent */
/*     atomic multipole moments to the field via a double loop */


/* Subroutine */ int udirect2a_(doublereal *field, doublereal *fieldp)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal r__, r2, ci;
    static integer ii, kk;
    static doublereal ck, bn[4], xr, yr, zr, fid[3], fkd[3], dir, dix, diy, 
	    diz, dkx, dky, dkz, dkr, qir, qkr, pdi, pti, fim[3], qix, qiy, 
	    qiz, qkx, qky, qkz, fkm[3], fip[3], fkp[3], dsc3, dsc5, dsc7, 
	    psc3, drr3, psc5, drr5, drr7, psc7, prr3, prr5, prr7, bfac;
    extern doublereal erfc_(doublereal *);
    static doublereal damp, qixx, qixy, qiyy, qixz, qizz, qiyz, qkxx, qkyy, 
	    qkzz, qkxy, qkxz, qkyz, exp2a, alsq2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal scale3, scale5, scale7, alsq2n, dscale[25000], pgamma;
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal ralpha, pscale[25000];
    extern /* Subroutine */ int switch_(char *, ftnlen);
    static doublereal expdamp;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define ip12_ref(a_1,a_2) polgrp_1.ip12[(a_2)*50 + a_1 - 51]
#define ip13_ref(a_1,a_2) polgrp_1.ip13[(a_2)*50 + a_1 - 51]
#define ip14_ref(a_1,a_2) polgrp_1.ip14[(a_2)*50 + a_1 - 51]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define fieldp_ref(a_1,a_2) fieldp[(a_2)*3 + a_1]



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




/*     check for multipoles and set cutoff coefficients */

    /* Parameter adjustments */
    fieldp -= 4;
    field -= 4;

    /* Function Body */
    if (mpole_1.npole == 0) {
	return 0;
    }
    switch_("EWALD", (ftnlen)5);

/*     set arrays needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pscale[i__ - 1] = 1.;
	dscale[i__ - 1] = 1.;
    }

/*     compute the real space portion of the Ewald summation */

    i__1 = mpole_1.npole - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	pdi = polar_1.pdamp[i__ - 1];
	pti = polar_1.thole[i__ - 1];
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
	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
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
	    pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = polpot_1.d1scale;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = polpot_1.d2scale;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = polpot_1.d3scale;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = polpot_1.d4scale;
	}
	i__2 = mpole_1.npole;
	for (k = i__ + 1; k <= i__2; ++k) {
	    kk = mpole_1.ipole[k - 1];
	    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
	    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
	    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
	    image_(&xr, &yr, &zr);
	    r2 = xr * xr + yr * yr + zr * zr;
	    if (r2 <= shunt_1.cut2) {
		r__ = sqrt(r2);
		ck = rpole_ref(1, k);
		dkx = rpole_ref(2, k);
		dky = rpole_ref(3, k);
		dkz = rpole_ref(4, k);
		qkxx = rpole_ref(5, k);
		qkxy = rpole_ref(6, k);
		qkxz = rpole_ref(7, k);
		qkyy = rpole_ref(9, k);
		qkyz = rpole_ref(10, k);
		qkzz = rpole_ref(13, k);

/*     calculate the error function damping terms */

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
		for (j = 1; j <= 3; ++j) {
		    bfac = (doublereal) (j + j - 1);
		    alsq2n = alsq2 * alsq2n;
		    bn[j] = (bfac * bn[j - 1] + alsq2n * exp2a) / r2;
		}

/*     compute the error function scaled and unscaled terms */

		scale3 = 1.;
		scale5 = 1.;
		scale7 = 1.;
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
			scale5 = 1. - expdamp * (1. - damp);
/* Computing 2nd power */
			d__1 = damp;
			scale7 = 1. - expdamp * (1. - damp + d__1 * d__1 * .6)
				;
		    }
		}
		dsc3 = scale3 * dscale[kk - 1];
		dsc5 = scale5 * dscale[kk - 1];
		dsc7 = scale7 * dscale[kk - 1];
		psc3 = scale3 * pscale[kk - 1];
		psc5 = scale5 * pscale[kk - 1];
		psc7 = scale7 * pscale[kk - 1];
		drr3 = (1. - dsc3) / (r__ * r2);
		drr5 = (1. - dsc5) * 3. / (r__ * r2 * r2);
		drr7 = (1. - dsc7) * 15. / (r__ * r2 * r2 * r2);
		prr3 = (1. - psc3) / (r__ * r2);
		prr5 = (1. - psc5) * 3. / (r__ * r2 * r2);
		prr7 = (1. - psc7) * 15. / (r__ * r2 * r2 * r2);
		dir = dix * xr + diy * yr + diz * zr;
		qix = qixx * xr + qixy * yr + qixz * zr;
		qiy = qixy * xr + qiyy * yr + qiyz * zr;
		qiz = qixz * xr + qiyz * yr + qizz * zr;
		qir = qix * xr + qiy * yr + qiz * zr;
		dkr = dkx * xr + dky * yr + dkz * zr;
		qkx = qkxx * xr + qkxy * yr + qkxz * zr;
		qky = qkxy * xr + qkyy * yr + qkyz * zr;
		qkz = qkxz * xr + qkyz * yr + qkzz * zr;
		qkr = qkx * xr + qky * yr + qkz * zr;
		fim[0] = -xr * (bn[1] * ck - bn[2] * dkr + bn[3] * qkr) - bn[
			1] * dkx + bn[2] * 2. * qkx;
		fim[1] = -yr * (bn[1] * ck - bn[2] * dkr + bn[3] * qkr) - bn[
			1] * dky + bn[2] * 2. * qky;
		fim[2] = -zr * (bn[1] * ck - bn[2] * dkr + bn[3] * qkr) - bn[
			1] * dkz + bn[2] * 2. * qkz;
		fkm[0] = xr * (bn[1] * ci + bn[2] * dir + bn[3] * qir) - bn[1]
			 * dix - bn[2] * 2. * qix;
		fkm[1] = yr * (bn[1] * ci + bn[2] * dir + bn[3] * qir) - bn[1]
			 * diy - bn[2] * 2. * qiy;
		fkm[2] = zr * (bn[1] * ci + bn[2] * dir + bn[3] * qir) - bn[1]
			 * diz - bn[2] * 2. * qiz;
		fid[0] = -xr * (drr3 * ck - drr5 * dkr + drr7 * qkr) - drr3 * 
			dkx + drr5 * 2. * qkx;
		fid[1] = -yr * (drr3 * ck - drr5 * dkr + drr7 * qkr) - drr3 * 
			dky + drr5 * 2. * qky;
		fid[2] = -zr * (drr3 * ck - drr5 * dkr + drr7 * qkr) - drr3 * 
			dkz + drr5 * 2. * qkz;
		fkd[0] = xr * (drr3 * ci + drr5 * dir + drr7 * qir) - drr3 * 
			dix - drr5 * 2. * qix;
		fkd[1] = yr * (drr3 * ci + drr5 * dir + drr7 * qir) - drr3 * 
			diy - drr5 * 2. * qiy;
		fkd[2] = zr * (drr3 * ci + drr5 * dir + drr7 * qir) - drr3 * 
			diz - drr5 * 2. * qiz;
		fip[0] = -xr * (prr3 * ck - prr5 * dkr + prr7 * qkr) - prr3 * 
			dkx + prr5 * 2. * qkx;
		fip[1] = -yr * (prr3 * ck - prr5 * dkr + prr7 * qkr) - prr3 * 
			dky + prr5 * 2. * qky;
		fip[2] = -zr * (prr3 * ck - prr5 * dkr + prr7 * qkr) - prr3 * 
			dkz + prr5 * 2. * qkz;
		fkp[0] = xr * (prr3 * ci + prr5 * dir + prr7 * qir) - prr3 * 
			dix - prr5 * 2. * qix;
		fkp[1] = yr * (prr3 * ci + prr5 * dir + prr7 * qir) - prr3 * 
			diy - prr5 * 2. * qiy;
		fkp[2] = zr * (prr3 * ci + prr5 * dir + prr7 * qir) - prr3 * 
			diz - prr5 * 2. * qiz;

/*     increment the field at each site due to this interaction */

		for (j = 1; j <= 3; ++j) {
		    field_ref(j, i__) = field_ref(j, i__) + fim[j - 1] - fid[
			    j - 1];
		    field_ref(j, k) = field_ref(j, k) + fkm[j - 1] - fkd[j - 
			    1];
		    fieldp_ref(j, i__) = fieldp_ref(j, i__) + fim[j - 1] - 
			    fip[j - 1];
		    fieldp_ref(j, k) = fieldp_ref(j, k) + fkm[j - 1] - fkp[j 
			    - 1];
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i12_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i13_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i14_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i15_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = 1.;
	}
    }

/*     periodic boundary for large cutoffs via replicates method */

    if (bound_1.use_replica__) {
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    pdi = polar_1.pdamp[i__ - 1];
	    pti = polar_1.thole[i__ - 1];
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
	    i__2 = couple_1.n12[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	    }
	    i__2 = couple_1.n13[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	    }
	    i__2 = couple_1.n14[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[i14_ref(j, ii) - 1] = polpot_1.p4scale;
	    }
	    i__2 = couple_1.n15[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	    }
	    i__2 = polgrp_1.np11[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip11_ref(j, ii) - 1] = polpot_1.d1scale;
	    }
	    i__2 = polgrp_1.np12[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip12_ref(j, ii) - 1] = polpot_1.d2scale;
	    }
	    i__2 = polgrp_1.np13[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip13_ref(j, ii) - 1] = polpot_1.d3scale;
	    }
	    i__2 = polgrp_1.np14[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip14_ref(j, ii) - 1] = polpot_1.d4scale;
	    }
	    i__2 = mpole_1.npole;
	    for (k = i__; k <= i__2; ++k) {
		kk = mpole_1.ipole[k - 1];
		ck = rpole_ref(1, k);
		dkx = rpole_ref(2, k);
		dky = rpole_ref(3, k);
		dkz = rpole_ref(4, k);
		qkxx = rpole_ref(5, k);
		qkxy = rpole_ref(6, k);
		qkxz = rpole_ref(7, k);
		qkyy = rpole_ref(9, k);
		qkyz = rpole_ref(10, k);
		qkzz = rpole_ref(13, k);
		i__3 = cell_1.ncell;
		for (m = 1; m <= i__3; ++m) {
		    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
		    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
		    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
		    imager_(&xr, &yr, &zr, &m);
		    r2 = xr * xr + yr * yr + zr * zr;

/*     calculate the error function damping terms */

		    if (r2 <= shunt_1.cut2) {
			r__ = sqrt(r2);
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
			for (j = 1; j <= 3; ++j) {
			    bfac = (doublereal) (j + j - 1);
			    alsq2n = alsq2 * alsq2n;
			    bn[j] = (bfac * bn[j - 1] + alsq2n * exp2a) / r2;
			}

/*     compute the error function scaled and unscaled terms */

			scale3 = 1.;
			scale5 = 1.;
			scale7 = 1.;
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
				scale5 = 1. - expdamp * (1. - damp);
/* Computing 2nd power */
				d__1 = damp;
				scale7 = 1. - expdamp * (1. - damp + d__1 * 
					d__1 * .6);
			    }
			}
			dsc3 = scale3;
			dsc5 = scale5;
			dsc7 = scale7;
			psc3 = scale3;
			psc5 = scale5;
			psc7 = scale7;
			if (bound_1.use_polymer__) {
			    if (r2 <= bound_1.polycut2) {
				dsc3 = scale3 * dscale[kk - 1];
				dsc5 = scale5 * dscale[kk - 1];
				dsc7 = scale7 * dscale[kk - 1];
				psc3 = scale3 * pscale[kk - 1];
				psc5 = scale5 * pscale[kk - 1];
				psc7 = scale7 * pscale[kk - 1];
			    }
			}
			drr3 = (1. - dsc3) / (r__ * r2);
			drr5 = (1. - dsc5) * 3. / (r__ * r2 * r2);
			drr7 = (1. - dsc7) * 15. / (r__ * r2 * r2 * r2);
			prr3 = (1. - psc3) / (r__ * r2);
			prr5 = (1. - psc5) * 3. / (r__ * r2 * r2);
			prr7 = (1. - psc7) * 15. / (r__ * r2 * r2 * r2);
			dir = dix * xr + diy * yr + diz * zr;
			qix = qixx * xr + qixy * yr + qixz * zr;
			qiy = qixy * xr + qiyy * yr + qiyz * zr;
			qiz = qixz * xr + qiyz * yr + qizz * zr;
			qir = qix * xr + qiy * yr + qiz * zr;
			dkr = dkx * xr + dky * yr + dkz * zr;
			qkx = qkxx * xr + qkxy * yr + qkxz * zr;
			qky = qkxy * xr + qkyy * yr + qkyz * zr;
			qkz = qkxz * xr + qkyz * yr + qkzz * zr;
			qkr = qkx * xr + qky * yr + qkz * zr;
			fim[0] = -xr * (bn[1] * ck - bn[2] * dkr + bn[3] * 
				qkr) - bn[1] * dkx + bn[2] * 2. * qkx;
			fim[1] = -yr * (bn[1] * ck - bn[2] * dkr + bn[3] * 
				qkr) - bn[1] * dky + bn[2] * 2. * qky;
			fim[2] = -zr * (bn[1] * ck - bn[2] * dkr + bn[3] * 
				qkr) - bn[1] * dkz + bn[2] * 2. * qkz;
			fkm[0] = xr * (bn[1] * ci + bn[2] * dir + bn[3] * qir)
				 - bn[1] * dix - bn[2] * 2. * qix;
			fkm[1] = yr * (bn[1] * ci + bn[2] * dir + bn[3] * qir)
				 - bn[1] * diy - bn[2] * 2. * qiy;
			fkm[2] = zr * (bn[1] * ci + bn[2] * dir + bn[3] * qir)
				 - bn[1] * diz - bn[2] * 2. * qiz;
			fid[0] = -xr * (drr3 * ck - drr5 * dkr + drr7 * qkr) 
				- drr3 * dkx + drr5 * 2. * qkx;
			fid[1] = -yr * (drr3 * ck - drr5 * dkr + drr7 * qkr) 
				- drr3 * dky + drr5 * 2. * qky;
			fid[2] = -zr * (drr3 * ck - drr5 * dkr + drr7 * qkr) 
				- drr3 * dkz + drr5 * 2. * qkz;
			fkd[0] = xr * (drr3 * ci + drr5 * dir + drr7 * qir) - 
				drr3 * dix - drr5 * 2. * qix;
			fkd[1] = yr * (drr3 * ci + drr5 * dir + drr7 * qir) - 
				drr3 * diy - drr5 * 2. * qiy;
			fkd[2] = zr * (drr3 * ci + drr5 * dir + drr7 * qir) - 
				drr3 * diz - drr5 * 2. * qiz;
			fip[0] = -xr * (prr3 * ck - prr5 * dkr + prr7 * qkr) 
				- prr3 * dkx + prr5 * 2. * qkx;
			fip[1] = -yr * (prr3 * ck - prr5 * dkr + prr7 * qkr) 
				- prr3 * dky + prr5 * 2. * qky;
			fip[2] = -zr * (prr3 * ck - prr5 * dkr + prr7 * qkr) 
				- prr3 * dkz + prr5 * 2. * qkz;
			fkp[0] = xr * (prr3 * ci + prr5 * dir + prr7 * qir) - 
				prr3 * dix - prr5 * 2. * qix;
			fkp[1] = yr * (prr3 * ci + prr5 * dir + prr7 * qir) - 
				prr3 * diy - prr5 * 2. * qiy;
			fkp[2] = zr * (prr3 * ci + prr5 * dir + prr7 * qir) - 
				prr3 * diz - prr5 * 2. * qiz;

/*     increment the field at each site due to this interaction */

			for (j = 1; j <= 3; ++j) {
			    field_ref(j, i__) = field_ref(j, i__) + fim[j - 1]
				     - fid[j - 1];
			    fieldp_ref(j, i__) = fieldp_ref(j, i__) + fim[j - 
				    1] - fip[j - 1];
			    if (ii != kk) {
				field_ref(j, k) = field_ref(j, k) + fkm[j - 1]
					 - fkd[j - 1];
				fieldp_ref(j, k) = fieldp_ref(j, k) + fkm[j - 
					1] - fkp[j - 1];
			    }
			}
		    }
		}
	    }

/*     reset interaction scaling coefficients for connected atoms */

	    i__2 = couple_1.n12[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[i12_ref(j, ii) - 1] = 1.;
	    }
	    i__2 = couple_1.n13[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[i13_ref(j, ii) - 1] = 1.;
	    }
	    i__2 = couple_1.n14[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[i14_ref(j, ii) - 1] = 1.;
	    }
	    i__2 = couple_1.n15[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[i15_ref(j, ii) - 1] = 1.;
	    }
	    i__2 = polgrp_1.np11[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip11_ref(j, ii) - 1] = 1.;
	    }
	    i__2 = polgrp_1.np12[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip12_ref(j, ii) - 1] = 1.;
	    }
	    i__2 = polgrp_1.np13[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip13_ref(j, ii) - 1] = 1.;
	    }
	    i__2 = polgrp_1.np14[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip14_ref(j, ii) - 1] = 1.;
	    }
	}
    }
    return 0;
} /* udirect2a_ */

#undef fieldp_ref
#undef rpole_ref
#undef field_ref
#undef ip14_ref
#undef ip13_ref
#undef ip12_ref
#undef ip11_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine udirect2b  --  Ewald real direct field via list  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "udirect2b" computes the real space contribution of the permanent */
/*     atomic multipole moments to the field via a neighbor list */


/* Subroutine */ int udirect2b_(doublereal *field, doublereal *fieldp)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal r__, r2, ci;
    static integer ii, kk;
    static doublereal ck, bn[4], xr, yr, zr, fid[3], fkd[3];
    static integer kkk;
    static doublereal dix, diy, diz, dkx, dky, dkz, dir, dkr, qir, qkr, pdi, 
	    pti, qix, qiy, qiz, qkx, qky, qkz, fim[3], fkm[3], fip[3], fkp[3],
	     dsc3, dsc5, dsc7, psc3, drr3, psc5, drr5, drr7, psc7, prr3, prr5,
	     prr7, bfac;
    extern doublereal erfc_(doublereal *);
    static doublereal damp, qixx, qixy, qiyy, qixz, qizz, qiyz, qkxx, qkyy, 
	    qkzz, qkxy, qkxz, qkyz, exp2a, alsq2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal scale3, scale5, scale7, alsq2n, dscale[25000], pgamma, 
	    ralpha, pscale[25000];
    extern /* Subroutine */ int switch_(char *, ftnlen);
    static doublereal expdamp;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define ip12_ref(a_1,a_2) polgrp_1.ip12[(a_2)*50 + a_1 - 51]
#define ip13_ref(a_1,a_2) polgrp_1.ip13[(a_2)*50 + a_1 - 51]
#define ip14_ref(a_1,a_2) polgrp_1.ip14[(a_2)*50 + a_1 - 51]
#define elst_ref(a_1,a_2) neigh_1.elst[(a_2)*1200 + a_1 - 1201]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define fieldp_ref(a_1,a_2) fieldp[(a_2)*3 + a_1]



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




/*     check for multipoles and set cutoff coefficients */

    /* Parameter adjustments */
    fieldp -= 4;
    field -= 4;

    /* Function Body */
    if (mpole_1.npole == 0) {
	return 0;
    }
    switch_("EWALD", (ftnlen)5);

/*     set arrays needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pscale[i__ - 1] = 1.;
	dscale[i__ - 1] = 1.;
    }

/*     compute the real space portion of the Ewald summation */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	pdi = polar_1.pdamp[i__ - 1];
	pti = polar_1.thole[i__ - 1];
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
	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
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
	    pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = polpot_1.d1scale;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = polpot_1.d2scale;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = polpot_1.d3scale;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = polpot_1.d4scale;
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
	    if (r2 <= shunt_1.cut2) {
		r__ = sqrt(r2);
		ck = rpole_ref(1, k);
		dkx = rpole_ref(2, k);
		dky = rpole_ref(3, k);
		dkz = rpole_ref(4, k);
		qkxx = rpole_ref(5, k);
		qkxy = rpole_ref(6, k);
		qkxz = rpole_ref(7, k);
		qkyy = rpole_ref(9, k);
		qkyz = rpole_ref(10, k);
		qkzz = rpole_ref(13, k);

/*     calculate the error function damping terms */

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
		for (j = 1; j <= 3; ++j) {
		    bfac = (doublereal) (j + j - 1);
		    alsq2n = alsq2 * alsq2n;
		    bn[j] = (bfac * bn[j - 1] + alsq2n * exp2a) / r2;
		}

/*     compute the error function scaled and unscaled terms */

		scale3 = 1.;
		scale5 = 1.;
		scale7 = 1.;
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
			scale5 = 1. - expdamp * (1. - damp);
/* Computing 2nd power */
			d__1 = damp;
			scale7 = 1. - expdamp * (1. - damp + d__1 * d__1 * .6)
				;
		    }
		}
		dsc3 = scale3 * dscale[kk - 1];
		dsc5 = scale5 * dscale[kk - 1];
		dsc7 = scale7 * dscale[kk - 1];
		psc3 = scale3 * pscale[kk - 1];
		psc5 = scale5 * pscale[kk - 1];
		psc7 = scale7 * pscale[kk - 1];
		drr3 = (1. - dsc3) / (r__ * r2);
		drr5 = (1. - dsc5) * 3. / (r__ * r2 * r2);
		drr7 = (1. - dsc7) * 15. / (r__ * r2 * r2 * r2);
		prr3 = (1. - psc3) / (r__ * r2);
		prr5 = (1. - psc5) * 3. / (r__ * r2 * r2);
		prr7 = (1. - psc7) * 15. / (r__ * r2 * r2 * r2);
		dir = dix * xr + diy * yr + diz * zr;
		qix = qixx * xr + qixy * yr + qixz * zr;
		qiy = qixy * xr + qiyy * yr + qiyz * zr;
		qiz = qixz * xr + qiyz * yr + qizz * zr;
		qir = qix * xr + qiy * yr + qiz * zr;
		dkr = dkx * xr + dky * yr + dkz * zr;
		qkx = qkxx * xr + qkxy * yr + qkxz * zr;
		qky = qkxy * xr + qkyy * yr + qkyz * zr;
		qkz = qkxz * xr + qkyz * yr + qkzz * zr;
		qkr = qkx * xr + qky * yr + qkz * zr;
		fim[0] = -xr * (bn[1] * ck - bn[2] * dkr + bn[3] * qkr) - bn[
			1] * dkx + bn[2] * 2. * qkx;
		fim[1] = -yr * (bn[1] * ck - bn[2] * dkr + bn[3] * qkr) - bn[
			1] * dky + bn[2] * 2. * qky;
		fim[2] = -zr * (bn[1] * ck - bn[2] * dkr + bn[3] * qkr) - bn[
			1] * dkz + bn[2] * 2. * qkz;
		fkm[0] = xr * (bn[1] * ci + bn[2] * dir + bn[3] * qir) - bn[1]
			 * dix - bn[2] * 2. * qix;
		fkm[1] = yr * (bn[1] * ci + bn[2] * dir + bn[3] * qir) - bn[1]
			 * diy - bn[2] * 2. * qiy;
		fkm[2] = zr * (bn[1] * ci + bn[2] * dir + bn[3] * qir) - bn[1]
			 * diz - bn[2] * 2. * qiz;
		fid[0] = -xr * (drr3 * ck - drr5 * dkr + drr7 * qkr) - drr3 * 
			dkx + drr5 * 2. * qkx;
		fid[1] = -yr * (drr3 * ck - drr5 * dkr + drr7 * qkr) - drr3 * 
			dky + drr5 * 2. * qky;
		fid[2] = -zr * (drr3 * ck - drr5 * dkr + drr7 * qkr) - drr3 * 
			dkz + drr5 * 2. * qkz;
		fkd[0] = xr * (drr3 * ci + drr5 * dir + drr7 * qir) - drr3 * 
			dix - drr5 * 2. * qix;
		fkd[1] = yr * (drr3 * ci + drr5 * dir + drr7 * qir) - drr3 * 
			diy - drr5 * 2. * qiy;
		fkd[2] = zr * (drr3 * ci + drr5 * dir + drr7 * qir) - drr3 * 
			diz - drr5 * 2. * qiz;
		fip[0] = -xr * (prr3 * ck - prr5 * dkr + prr7 * qkr) - prr3 * 
			dkx + prr5 * 2. * qkx;
		fip[1] = -yr * (prr3 * ck - prr5 * dkr + prr7 * qkr) - prr3 * 
			dky + prr5 * 2. * qky;
		fip[2] = -zr * (prr3 * ck - prr5 * dkr + prr7 * qkr) - prr3 * 
			dkz + prr5 * 2. * qkz;
		fkp[0] = xr * (prr3 * ci + prr5 * dir + prr7 * qir) - prr3 * 
			dix - prr5 * 2. * qix;
		fkp[1] = yr * (prr3 * ci + prr5 * dir + prr7 * qir) - prr3 * 
			diy - prr5 * 2. * qiy;
		fkp[2] = zr * (prr3 * ci + prr5 * dir + prr7 * qir) - prr3 * 
			diz - prr5 * 2. * qiz;

/*     increment the field at each site due to this interaction */

		for (j = 1; j <= 3; ++j) {
		    field_ref(j, i__) = field_ref(j, i__) + fim[j - 1] - fid[
			    j - 1];
		    field_ref(j, k) = field_ref(j, k) + fkm[j - 1] - fkd[j - 
			    1];
		    fieldp_ref(j, i__) = fieldp_ref(j, i__) + fim[j - 1] - 
			    fip[j - 1];
		    fieldp_ref(j, k) = fieldp_ref(j, k) + fkm[j - 1] - fkp[j 
			    - 1];
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i12_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i13_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i14_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i15_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = 1.;
	}
    }
    return 0;
} /* udirect2b_ */

#undef fieldp_ref
#undef rpole_ref
#undef field_ref
#undef elst_ref
#undef ip14_ref
#undef ip13_ref
#undef ip12_ref
#undef ip11_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine umutual1  --  Ewald recip mutual induced field  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "umutual1" computes the reciprocal space contribution of the */
/*     induced atomic dipole moments to the field */


/* Subroutine */ int umutual1_(doublereal *field, doublereal *fieldp)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int fftfront_(void);
    static doublereal dipfield1[75000]	/* was [3][25000] */, dipfield2[75000]
	    	/* was [3][25000] */, fdip_phi1__[250000]	/* was [10][
	    25000] */, fdip_phi2__[250000]	/* was [10][25000] */, a[9]	
	    /* was [3][3] */;
    static integer i__, j, k;
    extern /* Subroutine */ int grid_uind__(doublereal *, doublereal *), 
	    fphi_uind__(doublereal *, doublereal *, doublereal *);
    static doublereal term, fdip_sum_phi__[500000]	/* was [20][25000] */,
	     fuind[75000]	/* was [3][25000] */, fuinp[75000]	/* 
	    was [3][25000] */;
    extern /* Subroutine */ int fftback_(void);


#define dipfield1_ref(a_1,a_2) dipfield1[(a_2)*3 + a_1 - 4]
#define dipfield2_ref(a_1,a_2) dipfield2[(a_2)*3 + a_1 - 4]
#define fdip_phi1___ref(a_1,a_2) fdip_phi1__[(a_2)*10 + a_1 - 11]
#define fdip_phi2___ref(a_1,a_2) fdip_phi2__[(a_2)*10 + a_1 - 11]
#define a_ref(a_1,a_2) a[(a_2)*3 + a_1 - 4]
#define qfac_ref(a_1,a_2,a_3) pme_1.qfac[((a_3)*100 + (a_2))*100 + a_1 - \
10101]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define uinp_ref(a_1,a_2) polar_1.uinp[(a_2)*3 + a_1 - 4]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1]
#define recip_ref(a_1,a_2) boxes_1.recip[(a_2)*3 + a_1 - 4]
#define fuind_ref(a_1,a_2) fuind[(a_2)*3 + a_1 - 4]
#define qgrid_ref(a_1,a_2,a_3,a_4) pme_1.qgrid[(((a_4)*100 + (a_3))*100 + \
(a_2))*2 + a_1 - 20203]
#define fuinp_ref(a_1,a_2) fuinp[(a_2)*3 + a_1 - 4]
#define fieldp_ref(a_1,a_2) fieldp[(a_2)*3 + a_1]



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




/*     return if the Ewald coefficient is zero */

    /* Parameter adjustments */
    fieldp -= 4;
    field -= 4;

    /* Function Body */
    if (ewald_1.aewald < 1e-6) {
	return 0;
    }

/*     convert Cartesian dipoles to fractional coordinates */

    for (i__ = 1; i__ <= 3; ++i__) {
	a_ref(1, i__) = (doublereal) pme_1.nfft1 * recip_ref(i__, 1);
	a_ref(2, i__) = (doublereal) pme_1.nfft2 * recip_ref(i__, 2);
	a_ref(3, i__) = (doublereal) pme_1.nfft3 * recip_ref(i__, 3);
    }
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (k = 1; k <= 3; ++k) {
	    fuind_ref(k, i__) = a_ref(k, 1) * uind_ref(1, i__) + a_ref(k, 2) *
		     uind_ref(2, i__) + a_ref(k, 3) * uind_ref(3, i__);
	    fuinp_ref(k, i__) = a_ref(k, 1) * uinp_ref(1, i__) + a_ref(k, 2) *
		     uinp_ref(2, i__) + a_ref(k, 3) * uinp_ref(3, i__);
	}
    }

/*     assign PME grid and perform 3-D FFT forward transform */

    grid_uind__(fuind, fuinp);
    fftfront_();

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

/*     perform 3-D FFT backward transform and get field */

    fftback_();
    fphi_uind__(fdip_phi1__, fdip_phi2__, fdip_sum_phi__);

/*     convert the dipole fields from fractional to Cartesian */

    for (i__ = 1; i__ <= 3; ++i__) {
	a_ref(i__, 1) = (doublereal) pme_1.nfft1 * recip_ref(i__, 1);
	a_ref(i__, 2) = (doublereal) pme_1.nfft2 * recip_ref(i__, 2);
	a_ref(i__, 3) = (doublereal) pme_1.nfft3 * recip_ref(i__, 3);
    }
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (k = 1; k <= 3; ++k) {
	    dipfield1_ref(k, i__) = a_ref(k, 1) * fdip_phi1___ref(2, i__) + 
		    a_ref(k, 2) * fdip_phi1___ref(3, i__) + a_ref(k, 3) * 
		    fdip_phi1___ref(4, i__);
	    dipfield2_ref(k, i__) = a_ref(k, 1) * fdip_phi2___ref(2, i__) + 
		    a_ref(k, 2) * fdip_phi2___ref(3, i__) + a_ref(k, 3) * 
		    fdip_phi2___ref(4, i__);
	}
    }

/*     increment the field at each multipole site */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (k = 1; k <= 3; ++k) {
	    field_ref(k, i__) = field_ref(k, i__) - dipfield1_ref(k, i__);
	    fieldp_ref(k, i__) = fieldp_ref(k, i__) - dipfield2_ref(k, i__);
	}
    }
    return 0;
} /* umutual1_ */

#undef fieldp_ref
#undef fuinp_ref
#undef qgrid_ref
#undef fuind_ref
#undef recip_ref
#undef field_ref
#undef uinp_ref
#undef uind_ref
#undef qfac_ref
#undef a_ref
#undef fdip_phi2___ref
#undef fdip_phi1___ref
#undef dipfield2_ref
#undef dipfield1_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine umutual2a  --  Ewald real mutual field via loop  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "umutual2a" computes the real space contribution of the induced */
/*     atomic dipole moments to the field via a double loop */


/* Subroutine */ int umutual2a_(doublereal *field, doublereal *fieldp)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal r__, r2, bn[3];
    static integer ii, kk;
    static doublereal xr, yr, zr, rr3, rr5, fid[3], fkd[3], pdi, fip[3], fkp[
	    3], pti, bfac;
    extern doublereal erfc_(doublereal *);
    static doublereal fimd[3], damp, fkmd[3], fimp[3], fkmp[3], duir, dukr, 
	    duix, duiy, duiz, dukx, duky, puir, pukr, puix, puiy, puiz, dukz, 
	    pukx, puky, pukz, exp2a, alsq2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal scale3, scale5, alsq2n, dscale[25000], pgamma;
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal ralpha;
    extern /* Subroutine */ int switch_(char *, ftnlen);
    static doublereal expdamp;


#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define ip12_ref(a_1,a_2) polgrp_1.ip12[(a_2)*50 + a_1 - 51]
#define ip13_ref(a_1,a_2) polgrp_1.ip13[(a_2)*50 + a_1 - 51]
#define ip14_ref(a_1,a_2) polgrp_1.ip14[(a_2)*50 + a_1 - 51]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define uinp_ref(a_1,a_2) polar_1.uinp[(a_2)*3 + a_1 - 4]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1]
#define fieldp_ref(a_1,a_2) fieldp[(a_2)*3 + a_1]



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




/*     check for multipoles and set cutoff coefficients */

    /* Parameter adjustments */
    fieldp -= 4;
    field -= 4;

    /* Function Body */
    if (mpole_1.npole == 0) {
	return 0;
    }
    switch_("EWALD", (ftnlen)5);

/*     set array needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dscale[i__ - 1] = 1.;
    }

/*     compute the real space portion of the Ewald summation */

    i__1 = mpole_1.npole - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	pdi = polar_1.pdamp[i__ - 1];
	pti = polar_1.thole[i__ - 1];
	duix = uind_ref(1, i__);
	duiy = uind_ref(2, i__);
	duiz = uind_ref(3, i__);
	puix = uinp_ref(1, i__);
	puiy = uinp_ref(2, i__);
	puiz = uinp_ref(3, i__);
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = polpot_1.u1scale;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = polpot_1.u2scale;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = polpot_1.u3scale;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = polpot_1.u4scale;
	}
	i__2 = mpole_1.npole;
	for (k = i__ + 1; k <= i__2; ++k) {
	    kk = mpole_1.ipole[k - 1];
	    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
	    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
	    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
	    image_(&xr, &yr, &zr);
	    r2 = xr * xr + yr * yr + zr * zr;
	    if (r2 <= shunt_1.cut2) {
		r__ = sqrt(r2);
		dukx = uind_ref(1, k);
		duky = uind_ref(2, k);
		dukz = uind_ref(3, k);
		pukx = uinp_ref(1, k);
		puky = uinp_ref(2, k);
		pukz = uinp_ref(3, k);

/*     calculate the error function damping terms */

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
		for (j = 1; j <= 2; ++j) {
		    bfac = (doublereal) (j + j - 1);
		    alsq2n = alsq2 * alsq2n;
		    bn[j] = (bfac * bn[j - 1] + alsq2n * exp2a) / r2;
		}

/*     compute the error function scaled and unscaled terms */

		scale3 = dscale[kk - 1];
		scale5 = dscale[kk - 1];
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
			scale3 *= 1. - expdamp;
			scale5 *= 1. - (1. - damp) * expdamp;
		    }
		}
		rr3 = (1. - scale3) / (r__ * r2);
		rr5 = (1. - scale5) * 3. / (r__ * r2 * r2);
		duir = xr * duix + yr * duiy + zr * duiz;
		dukr = xr * dukx + yr * duky + zr * dukz;
		puir = xr * puix + yr * puiy + zr * puiz;
		pukr = xr * pukx + yr * puky + zr * pukz;
		fimd[0] = -bn[1] * dukx + bn[2] * dukr * xr;
		fimd[1] = -bn[1] * duky + bn[2] * dukr * yr;
		fimd[2] = -bn[1] * dukz + bn[2] * dukr * zr;
		fkmd[0] = -bn[1] * duix + bn[2] * duir * xr;
		fkmd[1] = -bn[1] * duiy + bn[2] * duir * yr;
		fkmd[2] = -bn[1] * duiz + bn[2] * duir * zr;
		fimp[0] = -bn[1] * pukx + bn[2] * pukr * xr;
		fimp[1] = -bn[1] * puky + bn[2] * pukr * yr;
		fimp[2] = -bn[1] * pukz + bn[2] * pukr * zr;
		fkmp[0] = -bn[1] * puix + bn[2] * puir * xr;
		fkmp[1] = -bn[1] * puiy + bn[2] * puir * yr;
		fkmp[2] = -bn[1] * puiz + bn[2] * puir * zr;
		fid[0] = -rr3 * dukx + rr5 * dukr * xr;
		fid[1] = -rr3 * duky + rr5 * dukr * yr;
		fid[2] = -rr3 * dukz + rr5 * dukr * zr;
		fkd[0] = -rr3 * duix + rr5 * duir * xr;
		fkd[1] = -rr3 * duiy + rr5 * duir * yr;
		fkd[2] = -rr3 * duiz + rr5 * duir * zr;
		fip[0] = -rr3 * pukx + rr5 * pukr * xr;
		fip[1] = -rr3 * puky + rr5 * pukr * yr;
		fip[2] = -rr3 * pukz + rr5 * pukr * zr;
		fkp[0] = -rr3 * puix + rr5 * puir * xr;
		fkp[1] = -rr3 * puiy + rr5 * puir * yr;
		fkp[2] = -rr3 * puiz + rr5 * puir * zr;

/*     increment the field at each site due to this interaction */

		for (j = 1; j <= 3; ++j) {
		    field_ref(j, i__) = field_ref(j, i__) + fimd[j - 1] - fid[
			    j - 1];
		    field_ref(j, k) = field_ref(j, k) + fkmd[j - 1] - fkd[j - 
			    1];
		    fieldp_ref(j, i__) = fieldp_ref(j, i__) + fimp[j - 1] - 
			    fip[j - 1];
		    fieldp_ref(j, k) = fieldp_ref(j, k) + fkmp[j - 1] - fkp[j 
			    - 1];
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = 1.;
	}
    }

/*     periodic boundary for large cutoffs via replicates method */

    if (bound_1.use_replica__) {
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    pdi = polar_1.pdamp[i__ - 1];
	    pti = polar_1.thole[i__ - 1];
	    duix = uind_ref(1, i__);
	    duiy = uind_ref(2, i__);
	    duiz = uind_ref(3, i__);
	    puix = uinp_ref(1, i__);
	    puiy = uinp_ref(2, i__);
	    puiz = uinp_ref(3, i__);
	    i__2 = polgrp_1.np11[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip11_ref(j, ii) - 1] = polpot_1.u1scale;
	    }
	    i__2 = polgrp_1.np12[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip12_ref(j, ii) - 1] = polpot_1.u2scale;
	    }
	    i__2 = polgrp_1.np13[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip13_ref(j, ii) - 1] = polpot_1.u3scale;
	    }
	    i__2 = polgrp_1.np14[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip14_ref(j, ii) - 1] = polpot_1.u4scale;
	    }
	    i__2 = mpole_1.npole;
	    for (k = i__; k <= i__2; ++k) {
		kk = mpole_1.ipole[k - 1];
		dukx = uind_ref(1, k);
		duky = uind_ref(2, k);
		dukz = uind_ref(3, k);
		pukx = uinp_ref(1, k);
		puky = uinp_ref(2, k);
		pukz = uinp_ref(3, k);
		i__3 = cell_1.ncell;
		for (m = 1; m <= i__3; ++m) {
		    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
		    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
		    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
		    imager_(&xr, &yr, &zr, &m);
		    r2 = xr * xr + yr * yr + zr * zr;

/*     calculate the error function damping terms */

		    if (r2 <= shunt_1.cut2) {
			r__ = sqrt(r2);
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
			for (j = 1; j <= 2; ++j) {
			    bfac = (doublereal) (j + j - 1);
			    alsq2n = alsq2 * alsq2n;
			    bn[j] = (bfac * bn[j - 1] + alsq2n * exp2a) / r2;
			}

/*     compute the error function scaled and unscaled terms */

			scale3 = 1.;
			scale5 = 1.;
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
			    }
			}
			if (bound_1.use_polymer__) {
			    if (r2 <= bound_1.polycut2) {
				scale3 *= dscale[kk - 1];
				scale5 *= dscale[kk - 1];
			    }
			}
			rr3 = (1. - scale3) / (r__ * r2);
			rr5 = (1. - scale5) * 3. / (r__ * r2 * r2);
			duir = xr * duix + yr * duiy + zr * duiz;
			dukr = xr * dukx + yr * duky + zr * dukz;
			puir = xr * puix + yr * puiy + zr * puiz;
			pukr = xr * pukx + yr * puky + zr * pukz;
			fimd[0] = -bn[1] * dukx + bn[2] * dukr * xr;
			fimd[1] = -bn[1] * duky + bn[2] * dukr * yr;
			fimd[2] = -bn[1] * dukz + bn[2] * dukr * zr;
			fkmd[0] = -bn[1] * duix + bn[2] * duir * xr;
			fkmd[1] = -bn[1] * duiy + bn[2] * duir * yr;
			fkmd[2] = -bn[1] * duiz + bn[2] * duir * zr;
			fimp[0] = -bn[1] * pukx + bn[2] * pukr * xr;
			fimp[1] = -bn[1] * puky + bn[2] * pukr * yr;
			fimp[2] = -bn[1] * pukz + bn[2] * pukr * zr;
			fkmp[0] = -bn[1] * puix + bn[2] * puir * xr;
			fkmp[1] = -bn[1] * puiy + bn[2] * puir * yr;
			fkmp[2] = -bn[1] * puiz + bn[2] * puir * zr;
			fid[0] = -rr3 * dukx + rr5 * dukr * xr;
			fid[1] = -rr3 * duky + rr5 * dukr * yr;
			fid[2] = -rr3 * dukz + rr5 * dukr * zr;
			fkd[0] = -rr3 * duix + rr5 * duir * xr;
			fkd[1] = -rr3 * duiy + rr5 * duir * yr;
			fkd[2] = -rr3 * duiz + rr5 * duir * zr;
			fip[0] = -rr3 * pukx + rr5 * pukr * xr;
			fip[1] = -rr3 * puky + rr5 * pukr * yr;
			fip[2] = -rr3 * pukz + rr5 * pukr * zr;
			fkp[0] = -rr3 * puix + rr5 * puir * xr;
			fkp[1] = -rr3 * puiy + rr5 * puir * yr;
			fkp[2] = -rr3 * puiz + rr5 * puir * zr;

/*     increment the field at each site due to this interaction */

			for (j = 1; j <= 3; ++j) {
			    field_ref(j, i__) = field_ref(j, i__) + fimd[j - 
				    1] - fid[j - 1];
			    fieldp_ref(j, i__) = fieldp_ref(j, i__) + fimp[j 
				    - 1] - fip[j - 1];
			    if (ii != kk) {
				field_ref(j, k) = field_ref(j, k) + fkmd[j - 
					1] - fkd[j - 1];
				fieldp_ref(j, k) = fieldp_ref(j, k) + fkmp[j 
					- 1] - fkp[j - 1];
			    }
			}
		    }
		}
	    }

/*     reset interaction scaling coefficients for connected atoms */

	    i__2 = polgrp_1.np11[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip11_ref(j, ii) - 1] = 1.;
	    }
	    i__2 = polgrp_1.np12[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip12_ref(j, ii) - 1] = 1.;
	    }
	    i__2 = polgrp_1.np13[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip13_ref(j, ii) - 1] = 1.;
	    }
	    i__2 = polgrp_1.np14[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		dscale[ip14_ref(j, ii) - 1] = 1.;
	    }
	}
    }
    return 0;
} /* umutual2a_ */

#undef fieldp_ref
#undef field_ref
#undef uinp_ref
#undef uind_ref
#undef ip14_ref
#undef ip13_ref
#undef ip12_ref
#undef ip11_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine umutual2b  --  Ewald real mutual field via list  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "umutual2b" computes the real space contribution of the induced */
/*     atomic dipole moments to the field via a neighbor list */


/* Subroutine */ int umutual2b_(doublereal *field, doublereal *fieldp)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal r__, r2, bn[3];
    static integer ii, kk;
    static doublereal xr, yr, zr, rr3, rr5, fid[3], fkd[3];
    static integer kkk;
    static doublereal pdi, fip[3], fkp[3], pti, bfac;
    extern doublereal erfc_(doublereal *);
    static doublereal fimd[3], damp, fkmd[3], fimp[3], fkmp[3], duir, dukr, 
	    duix, duiy, duiz, dukx, duky, puir, pukr, puix, puiy, puiz, dukz, 
	    pukx, puky, pukz, exp2a, alsq2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal scale3, scale5, alsq2n, dscale[25000], pgamma, ralpha;
    extern /* Subroutine */ int switch_(char *, ftnlen);
    static doublereal expdamp;


#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define ip12_ref(a_1,a_2) polgrp_1.ip12[(a_2)*50 + a_1 - 51]
#define ip13_ref(a_1,a_2) polgrp_1.ip13[(a_2)*50 + a_1 - 51]
#define ip14_ref(a_1,a_2) polgrp_1.ip14[(a_2)*50 + a_1 - 51]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define elst_ref(a_1,a_2) neigh_1.elst[(a_2)*1200 + a_1 - 1201]
#define uinp_ref(a_1,a_2) polar_1.uinp[(a_2)*3 + a_1 - 4]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1]
#define fieldp_ref(a_1,a_2) fieldp[(a_2)*3 + a_1]



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




/*     check for multipoles and set cutoff coefficients */

    /* Parameter adjustments */
    fieldp -= 4;
    field -= 4;

    /* Function Body */
    if (mpole_1.npole == 0) {
	return 0;
    }
    switch_("EWALD", (ftnlen)5);

/*     set array needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dscale[i__ - 1] = 1.;
    }

/*     compute the real space portion of the Ewald summation */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	pdi = polar_1.pdamp[i__ - 1];
	pti = polar_1.thole[i__ - 1];
	duix = uind_ref(1, i__);
	duiy = uind_ref(2, i__);
	duiz = uind_ref(3, i__);
	puix = uinp_ref(1, i__);
	puiy = uinp_ref(2, i__);
	puiz = uinp_ref(3, i__);
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = polpot_1.u1scale;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = polpot_1.u2scale;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = polpot_1.u3scale;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = polpot_1.u4scale;
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
	    if (r2 <= shunt_1.cut2) {
		r__ = sqrt(r2);
		dukx = uind_ref(1, k);
		duky = uind_ref(2, k);
		dukz = uind_ref(3, k);
		pukx = uinp_ref(1, k);
		puky = uinp_ref(2, k);
		pukz = uinp_ref(3, k);

/*     calculate the error function damping terms */

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
		for (j = 1; j <= 2; ++j) {
		    bfac = (doublereal) (j + j - 1);
		    alsq2n = alsq2 * alsq2n;
		    bn[j] = (bfac * bn[j - 1] + alsq2n * exp2a) / r2;
		}

/*     compute the error function scaled and unscaled terms */

		scale3 = dscale[kk - 1];
		scale5 = dscale[kk - 1];
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
			scale3 *= 1. - expdamp;
			scale5 *= 1. - (1. - damp) * expdamp;
		    }
		}
		rr3 = (1. - scale3) / (r__ * r2);
		rr5 = (1. - scale5) * 3. / (r__ * r2 * r2);
		duir = xr * duix + yr * duiy + zr * duiz;
		dukr = xr * dukx + yr * duky + zr * dukz;
		puir = xr * puix + yr * puiy + zr * puiz;
		pukr = xr * pukx + yr * puky + zr * pukz;
		fimd[0] = -bn[1] * dukx + bn[2] * dukr * xr;
		fimd[1] = -bn[1] * duky + bn[2] * dukr * yr;
		fimd[2] = -bn[1] * dukz + bn[2] * dukr * zr;
		fkmd[0] = -bn[1] * duix + bn[2] * duir * xr;
		fkmd[1] = -bn[1] * duiy + bn[2] * duir * yr;
		fkmd[2] = -bn[1] * duiz + bn[2] * duir * zr;
		fimp[0] = -bn[1] * pukx + bn[2] * pukr * xr;
		fimp[1] = -bn[1] * puky + bn[2] * pukr * yr;
		fimp[2] = -bn[1] * pukz + bn[2] * pukr * zr;
		fkmp[0] = -bn[1] * puix + bn[2] * puir * xr;
		fkmp[1] = -bn[1] * puiy + bn[2] * puir * yr;
		fkmp[2] = -bn[1] * puiz + bn[2] * puir * zr;
		fid[0] = -rr3 * dukx + rr5 * dukr * xr;
		fid[1] = -rr3 * duky + rr5 * dukr * yr;
		fid[2] = -rr3 * dukz + rr5 * dukr * zr;
		fkd[0] = -rr3 * duix + rr5 * duir * xr;
		fkd[1] = -rr3 * duiy + rr5 * duir * yr;
		fkd[2] = -rr3 * duiz + rr5 * duir * zr;
		fip[0] = -rr3 * pukx + rr5 * pukr * xr;
		fip[1] = -rr3 * puky + rr5 * pukr * yr;
		fip[2] = -rr3 * pukz + rr5 * pukr * zr;
		fkp[0] = -rr3 * puix + rr5 * puir * xr;
		fkp[1] = -rr3 * puiy + rr5 * puir * yr;
		fkp[2] = -rr3 * puiz + rr5 * puir * zr;

/*     increment the field at each site due to this interaction */

		for (j = 1; j <= 3; ++j) {
		    field_ref(j, i__) = field_ref(j, i__) + fimd[j - 1] - fid[
			    j - 1];
		    field_ref(j, k) = field_ref(j, k) + fkmd[j - 1] - fkd[j - 
			    1];
		    fieldp_ref(j, i__) = fieldp_ref(j, i__) + fimp[j - 1] - 
			    fip[j - 1];
		    fieldp_ref(j, k) = fieldp_ref(j, k) + fkmp[j - 1] - fkp[j 
			    - 1];
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = 1.;
	}
    }
    return 0;
} /* umutual2b_ */

#undef fieldp_ref
#undef field_ref
#undef uinp_ref
#undef elst_ref
#undef uind_ref
#undef ip14_ref
#undef ip13_ref
#undef ip12_ref
#undef ip11_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine induce0e  --  Kirkwood SCRF induced dipoles  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "induce0e" computes the induced dipole moments at polarizable */
/*     sites for generalized Kirkwood SCRF and vacuum environments */


/* Subroutine */ int induce0e_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Determination of Induced Dipole\002,\002"
	    " Moments :\002,//,4x,\002Iter\002,8x,\002RMS Change (Debyes)\002"
	    ",/)";
    static char fmt_20[] = "(i8,7x,f16.10)";
    static char fmt_30[] = "(/,\002 INDUCE  --  Warning, Induced Dipole"
	    "s\002,\002 are not Converged\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal gkfieldp[75000]	/* was [3][25000] */, a[12]	/* 
	    was [4][3] */;
    static integer i__, j, k;
    static doublereal r__, r2, expcdexpc, fc, fd, gc[4], ci;
    static integer ii, kk;
    static doublereal ck, fq, gf, xr, yr, zr, gf2, gf3, gf5, gf7, rb2, rr3, 
	    rr5, rr7, xr2, yr2, zr2, fid[3], fkd[3], dir, dkr, pdi, eps, rbi, 
	    rbk, fip[3], qir, qkr, pti, fkp[3], qix, qiy, qiz, qkx, uxi, uyi, 
	    uzi, uxk, uyk, uzk, qky, qkz, gux[10], guy[10], guz[10], damp, 
	    fids[3];
    static logical done;
    static doublereal fkds[3], epsd, fgrp;
    static integer iter;
    static doublereal duir, dukr, duix, duiy, duiz, dukx, duky, dukz, puir, 
	    pukr, puix, puiy, puiz, pukx, qxxi, qxyi, qxzi, qyyi, qyzi, qzzi, 
	    qxxk, qxyk, qxzk, qyyk, qyzk, qzzk, puky, pukz, epsp, expc, gqxx[
	    4], gqxy[4], gqxz[4], gqyy[4], gqyz[4], gqzz[4], expc1, fips[3], 
	    fkps[3], udir[75000]	/* was [3][25000] */, uold[75000]	
	    /* was [3][25000] */, field[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int fatal_(void);
    static doublereal dexpc, duirs, dukrs, duixs, duiys, duizs, dukxs, dukys, 
	    dukzs, puixs, puiys, puizs, pukxs, pukys, pukzs, puirs, pukrs, 
	    epsds, epsps, scale3, scale5, scale7, udirp[75000]	/* was [3][
	    25000] */, udirs[75000]	/* was [3][25000] */, uoldp[75000]	
	    /* was [3][25000] */, uolds[75000]	/* was [3][25000] */, dscale[
	    25000], pgamma, fieldp[75000]	/* was [3][25000] */, fields[
	    75000]	/* was [3][25000] */, pscale[25000], epsold, dwater, 
	    udirps[75000]	/* was [3][25000] */, uoldps[75000]	/* 
	    was [3][25000] */;
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *), prterr_(void);
    static doublereal gkfield[75000]	/* was [3][25000] */;
    static logical proceed;
    static doublereal fieldps[75000]	/* was [3][25000] */, expdamp;
    static integer maxiter;
    static doublereal expterm;

    /* Fortran I/O blocks */
    static cilist io___692 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___693 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___694 = { 0, 0, 0, fmt_30, 0 };



#define gkfieldp_ref(a_1,a_2) gkfieldp[(a_2)*3 + a_1 - 4]
#define a_ref(a_1,a_2) a[(a_2)*4 + a_1 - 0]
#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define ip12_ref(a_1,a_2) polgrp_1.ip12[(a_2)*50 + a_1 - 51]
#define ip13_ref(a_1,a_2) polgrp_1.ip13[(a_2)*50 + a_1 - 51]
#define ip14_ref(a_1,a_2) polgrp_1.ip14[(a_2)*50 + a_1 - 51]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define uinp_ref(a_1,a_2) polar_1.uinp[(a_2)*3 + a_1 - 4]
#define udir_ref(a_1,a_2) udir[(a_2)*3 + a_1 - 4]
#define uold_ref(a_1,a_2) uold[(a_2)*3 + a_1 - 4]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define uinds_ref(a_1,a_2) polar_1.uinds[(a_2)*3 + a_1 - 4]
#define uinps_ref(a_1,a_2) polar_1.uinps[(a_2)*3 + a_1 - 4]
#define udirp_ref(a_1,a_2) udirp[(a_2)*3 + a_1 - 4]
#define udirs_ref(a_1,a_2) udirs[(a_2)*3 + a_1 - 4]
#define uoldp_ref(a_1,a_2) uoldp[(a_2)*3 + a_1 - 4]
#define uolds_ref(a_1,a_2) uolds[(a_2)*3 + a_1 - 4]
#define fieldp_ref(a_1,a_2) fieldp[(a_2)*3 + a_1 - 4]
#define fields_ref(a_1,a_2) fields[(a_2)*3 + a_1 - 4]
#define udirps_ref(a_1,a_2) udirps[(a_2)*3 + a_1 - 4]
#define uoldps_ref(a_1,a_2) uoldps[(a_2)*3 + a_1 - 4]
#define gkfield_ref(a_1,a_2) gkfield[(a_2)*3 + a_1 - 4]
#define fieldps_ref(a_1,a_2) fieldps[(a_2)*3 + a_1 - 4]



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




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  gk.i  --  parameters for generalized Kirkwood solvation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     gkr      generalized Kirkwood cavity radii for atom types */
/*     gkc      tuning parameter exponent in the f(GB) function */




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
/*     ##  inform.i  --  control values for I/O and program flow  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     digits    decimal places output for energy and coordinates */
/*     iprint    steps between status printing (0=no printing) */
/*     iwrite    steps between coordinate dumps (0=no dumps) */
/*     isend     steps between socket communication (0=no sockets) */
/*     verbose   logical flag to turn on extra information */
/*     debug     logical flag to turn on full debug printing */
/*     holdup    logical flag to wait for carriage return on exit */
/*     abort     logical flag to stop execution at next chance */




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




/*     zero out the induced dipoles and the fields at each site; */
/*     uind & uinp are vacuum dipoles, uinds & uinps are SCRF dipoles, */
/*     and field & fieldp are solute fields, gkfield & gkfieldp */
/*     are reaction fields */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    uind_ref(j, i__) = 0.;
	    uinp_ref(j, i__) = 0.;
	    uinds_ref(j, i__) = 0.;
	    uinps_ref(j, i__) = 0.;
	    field_ref(j, i__) = 0.;
	    fieldp_ref(j, i__) = 0.;
	    gkfield_ref(j, i__) = 0.;
	    gkfieldp_ref(j, i__) = 0.;
	}
    }
    if (! potent_1.use_polar__ && ! potent_1.use_solv__) {
	return 0;
    }
    dwater = 78.3;
    fc = (1. - dwater) * 1. / (dwater * 1. + 0.);
    fd = (1. - dwater) * 2. / (dwater * 2. + 1.);
    fq = (1. - dwater) * 3. / (dwater * 3. + 2.);

/*     set the switching function coefficients */

    switch_("MPOLE", (ftnlen)5);

/*     set arrays needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pscale[i__ - 1] = 1.;
	dscale[i__ - 1] = 1.;
    }

/*     compute the direct induced dipole moment at each atom, and */
/*     another set that also includes RF due to permanent multipoles */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	pdi = polar_1.pdamp[i__ - 1];
	pti = polar_1.thole[i__ - 1];
	rbi = solute_1.rborn[ii - 1];
	ci = rpole_ref(1, i__);
	uxi = rpole_ref(2, i__);
	uyi = rpole_ref(3, i__);
	uzi = rpole_ref(4, i__);
	qxxi = rpole_ref(5, i__);
	qxyi = rpole_ref(6, i__);
	qxzi = rpole_ref(7, i__);
	qyyi = rpole_ref(9, i__);
	qyzi = rpole_ref(10, i__);
	qzzi = rpole_ref(13, i__);
	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
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
	    pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = polpot_1.d1scale;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = polpot_1.d2scale;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = polpot_1.d3scale;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = polpot_1.d4scale;
	}
	i__2 = mpole_1.npole;
	for (k = i__; k <= i__2; ++k) {
	    kk = mpole_1.ipole[k - 1];
	    rbk = solute_1.rborn[kk - 1];
	    proceed = TRUE_;
	    if (group_1.use_intra__) {
		groups_(&proceed, &fgrp, &ii, &kk, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
		yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
		zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
		xr2 = xr * xr;
		yr2 = yr * yr;
		zr2 = zr * zr;
		r2 = xr2 + yr2 + zr2;
		if (r2 <= shunt_1.off2) {
		    r__ = sqrt(r2);
		    ck = rpole_ref(1, k);
		    uxk = rpole_ref(2, k);
		    uyk = rpole_ref(3, k);
		    uzk = rpole_ref(4, k);
		    qxxk = rpole_ref(5, k);
		    qxyk = rpole_ref(6, k);
		    qxzk = rpole_ref(7, k);
		    qyyk = rpole_ref(9, k);
		    qyzk = rpole_ref(10, k);
		    qzzk = rpole_ref(13, k);

/*     self-interactions for the solute field are skipped */

		    if (i__ != k) {
			scale3 = 1.;
			scale5 = 1.;
			scale7 = 1.;
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
				scale5 = 1. - expdamp * (1. - damp);
/* Computing 2nd power */
				d__1 = damp;
				scale7 = 1. - expdamp * (1. - damp + d__1 * 
					d__1 * .6);
			    }
			}
			rr3 = scale3 / (r__ * r2);
			rr5 = scale5 * 3. / (r__ * r2 * r2);
			rr7 = scale7 * 15. / (r__ * r2 * r2 * r2);
			dir = uxi * xr + uyi * yr + uzi * zr;
			qix = qxxi * xr + qxyi * yr + qxzi * zr;
			qiy = qxyi * xr + qyyi * yr + qyzi * zr;
			qiz = qxzi * xr + qyzi * yr + qzzi * zr;
			qir = qix * xr + qiy * yr + qiz * zr;
			dkr = uxk * xr + uyk * yr + uzk * zr;
			qkx = qxxk * xr + qxyk * yr + qxzk * zr;
			qky = qxyk * xr + qyyk * yr + qyzk * zr;
			qkz = qxzk * xr + qyzk * yr + qzzk * zr;
			qkr = qkx * xr + qky * yr + qkz * zr;
			fid[0] = -xr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - 
				rr3 * uxk + rr5 * 2. * qkx;
			fid[1] = -yr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - 
				rr3 * uyk + rr5 * 2. * qky;
			fid[2] = -zr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - 
				rr3 * uzk + rr5 * 2. * qkz;
			fkd[0] = xr * (rr3 * ci + rr5 * dir + rr7 * qir) - 
				rr3 * uxi - rr5 * 2. * qix;
			fkd[1] = yr * (rr3 * ci + rr5 * dir + rr7 * qir) - 
				rr3 * uyi - rr5 * 2. * qiy;
			fkd[2] = zr * (rr3 * ci + rr5 * dir + rr7 * qir) - 
				rr3 * uzi - rr5 * 2. * qiz;
			for (j = 1; j <= 3; ++j) {
			    field_ref(j, i__) = field_ref(j, i__) + fid[j - 1]
				     * dscale[kk - 1];
			    field_ref(j, k) = field_ref(j, k) + fkd[j - 1] * 
				    dscale[kk - 1];
			    fieldp_ref(j, i__) = fieldp_ref(j, i__) + fid[j - 
				    1] * pscale[kk - 1];
			    fieldp_ref(j, k) = fieldp_ref(j, k) + fkd[j - 1] *
				     pscale[kk - 1];
			}
		    }
		    rb2 = rbi * rbk;
		    expterm = exp(-r2 / (gk_1.gkc * rb2));
		    expc = expterm / gk_1.gkc;
		    dexpc = -2. / (gk_1.gkc * rb2);
		    gf2 = 1. / (r2 + rb2 * expterm);
		    gf = sqrt(gf2);
		    gf3 = gf2 * gf;
		    gf5 = gf3 * gf2;
		    gf7 = gf5 * gf2;

/*     reaction potential auxiliary terms */

		    a_ref(0, 0) = gf;
		    a_ref(1, 0) = -gf3;
		    a_ref(2, 0) = gf5 * 3.;
		    a_ref(3, 0) = gf7 * -15.;

/*     reaction potential gradient auxiliary terms */

		    expc1 = 1. - expc;
		    a_ref(0, 1) = expc1 * a_ref(1, 0);
		    a_ref(1, 1) = expc1 * a_ref(2, 0);
		    a_ref(2, 1) = expc1 * a_ref(3, 0);

/*     dipole second reaction potential gradient auxiliary term */

		    expcdexpc = -expc * dexpc;
		    a_ref(1, 2) = expc1 * a_ref(2, 1) + expcdexpc * a_ref(2, 
			    0);

/*     multiply the auxillary terms by dielectric functions */

		    a_ref(0, 1) = fc * a_ref(0, 1);
		    a_ref(1, 0) = fd * a_ref(1, 0);
		    a_ref(1, 1) = fd * a_ref(1, 1);
		    a_ref(1, 2) = fd * a_ref(1, 2);
		    a_ref(2, 0) = fq * a_ref(2, 0);
		    a_ref(2, 1) = fq * a_ref(2, 1);

/*     unweighted dipole reaction potential tensor */

		    gux[0] = xr * a_ref(1, 0);
		    guy[0] = yr * a_ref(1, 0);
		    guz[0] = zr * a_ref(1, 0);

/*     unweighted reaction potential gradient tensor */

		    gc[1] = xr * a_ref(0, 1);
		    gc[2] = yr * a_ref(0, 1);
		    gc[3] = zr * a_ref(0, 1);
		    gux[1] = a_ref(1, 0) + xr2 * a_ref(1, 1);
		    gux[2] = xr * yr * a_ref(1, 1);
		    gux[3] = xr * zr * a_ref(1, 1);
		    guy[1] = gux[2];
		    guy[2] = a_ref(1, 0) + yr2 * a_ref(1, 1);
		    guy[3] = yr * zr * a_ref(1, 1);
		    guz[1] = gux[3];
		    guz[2] = guy[3];
		    guz[3] = a_ref(1, 0) + zr2 * a_ref(1, 1);
		    gqxx[1] = xr * (a_ref(2, 0) * 2. + xr2 * a_ref(2, 1));
		    gqxx[2] = yr * xr2 * a_ref(2, 1);
		    gqxx[3] = zr * xr2 * a_ref(2, 1);
		    gqyy[1] = xr * yr2 * a_ref(2, 1);
		    gqyy[2] = yr * (a_ref(2, 0) * 2. + yr2 * a_ref(2, 1));
		    gqyy[3] = zr * yr2 * a_ref(2, 1);
		    gqzz[1] = xr * zr2 * a_ref(2, 1);
		    gqzz[2] = yr * zr2 * a_ref(2, 1);
		    gqzz[3] = zr * (a_ref(2, 0) * 2. + zr2 * a_ref(2, 1));
		    gqxy[1] = yr * (a_ref(2, 0) + xr2 * a_ref(2, 1));
		    gqxy[2] = xr * (a_ref(2, 0) + yr2 * a_ref(2, 1));
		    gqxy[3] = zr * xr * yr * a_ref(2, 1);
		    gqxz[1] = zr * (a_ref(2, 0) + xr2 * a_ref(2, 1));
		    gqxz[2] = gqxy[3];
		    gqxz[3] = xr * (a_ref(2, 0) + zr2 * a_ref(2, 1));
		    gqyz[1] = gqxy[3];
		    gqyz[2] = zr * (a_ref(2, 0) + yr2 * a_ref(2, 1));
		    gqyz[3] = yr * (a_ref(2, 0) + zr2 * a_ref(2, 1));

/*     unweighted dipole second reaction potential gradient tensor */

		    gux[4] = xr * (a_ref(1, 1) * 3. + xr2 * a_ref(1, 2));
		    gux[5] = yr * (a_ref(1, 1) + xr2 * a_ref(1, 2));
		    gux[6] = zr * (a_ref(1, 1) + xr2 * a_ref(1, 2));
		    gux[7] = xr * (a_ref(1, 1) + yr2 * a_ref(1, 2));
		    gux[8] = zr * xr * yr * a_ref(1, 2);
		    gux[9] = xr * (a_ref(1, 1) + zr2 * a_ref(1, 2));
		    guy[4] = yr * (a_ref(1, 1) + xr2 * a_ref(1, 2));
		    guy[5] = xr * (a_ref(1, 1) + yr2 * a_ref(1, 2));
		    guy[6] = gux[8];
		    guy[7] = yr * (a_ref(1, 1) * 3. + yr2 * a_ref(1, 2));
		    guy[8] = zr * (a_ref(1, 1) + yr2 * a_ref(1, 2));
		    guy[9] = yr * (a_ref(1, 1) + zr2 * a_ref(1, 2));
		    guz[4] = zr * (a_ref(1, 1) + xr2 * a_ref(1, 2));
		    guz[5] = gux[8];
		    guz[6] = xr * (a_ref(1, 1) + zr2 * a_ref(1, 2));
		    guz[7] = zr * (a_ref(1, 1) + yr2 * a_ref(1, 2));
		    guz[8] = yr * (a_ref(1, 1) + zr2 * a_ref(1, 2));
		    guz[9] = zr * (a_ref(1, 1) * 3. + zr2 * a_ref(1, 2));

/*     generalized Kirkwood permanent reaction field */

		    fid[0] = uxk * gux[1] + uyk * gux[2] + uzk * gux[3] + (ck 
			    * gux[0] + qxxk * gux[4] + qyyk * gux[7] + qzzk * 
			    gux[9] + (qxyk * gux[5] + qxzk * gux[6] + qyzk * 
			    gux[8]) * 2.) * .5 + (ck * gc[1] + qxxk * gqxx[1] 
			    + qyyk * gqyy[1] + qzzk * gqzz[1] + (qxyk * gqxy[
			    1] + qxzk * gqxz[1] + qyzk * gqyz[1]) * 2.) * .5;
		    fid[1] = uxk * guy[1] + uyk * guy[2] + uzk * guy[3] + (ck 
			    * guy[0] + qxxk * guy[4] + qyyk * guy[7] + qzzk * 
			    guy[9] + (qxyk * guy[5] + qxzk * guy[6] + qyzk * 
			    guy[8]) * 2.) * .5 + (ck * gc[2] + qxxk * gqxx[2] 
			    + qyyk * gqyy[2] + qzzk * gqzz[2] + (qxyk * gqxy[
			    2] + qxzk * gqxz[2] + qyzk * gqyz[2]) * 2.) * .5;
		    fid[2] = uxk * guz[1] + uyk * guz[2] + uzk * guz[3] + (ck 
			    * guz[0] + qxxk * guz[4] + qyyk * guz[7] + qzzk * 
			    guz[9] + (qxyk * guz[5] + qxzk * guz[6] + qyzk * 
			    guz[8]) * 2.) * .5 + (ck * gc[3] + qxxk * gqxx[3] 
			    + qyyk * gqyy[3] + qzzk * gqzz[3] + (qxyk * gqxy[
			    3] + qxzk * gqxz[3] + qyzk * gqyz[3]) * 2.) * .5;
		    fkd[0] = uxi * gux[1] + uyi * gux[2] + uzi * gux[3] - (ci 
			    * gux[0] + qxxi * gux[4] + qyyi * gux[7] + qzzi * 
			    gux[9] + (qxyi * gux[5] + qxzi * gux[6] + qyzi * 
			    gux[8]) * 2.) * .5 - (ci * gc[1] + qxxi * gqxx[1] 
			    + qyyi * gqyy[1] + qzzi * gqzz[1] + (qxyi * gqxy[
			    1] + qxzi * gqxz[1] + qyzi * gqyz[1]) * 2.) * .5;
		    fkd[1] = uxi * guy[1] + uyi * guy[2] + uzi * guy[3] - (ci 
			    * guy[0] + qxxi * guy[4] + qyyi * guy[7] + qzzi * 
			    guy[9] + (qxyi * guy[5] + qxzi * guy[6] + qyzi * 
			    guy[8]) * 2.) * .5 - (ci * gc[2] + qxxi * gqxx[2] 
			    + qyyi * gqyy[2] + qzzi * gqzz[2] + (qxyi * gqxy[
			    2] + qxzi * gqxz[2] + qyzi * gqyz[2]) * 2.) * .5;
		    fkd[2] = uxi * guz[1] + uyi * guz[2] + uzi * guz[3] - (ci 
			    * guz[0] + qxxi * guz[4] + qyyi * guz[7] + qzzi * 
			    guz[9] + (qxyi * guz[5] + qxzi * guz[6] + qyzi * 
			    guz[8]) * 2.) * .5 - (ci * gc[3] + qxxi * gqxx[3] 
			    + qyyi * gqyy[3] + qzzi * gqzz[3] + (qxyi * gqxy[
			    3] + qxzi * gqxz[3] + qyzi * gqyz[3]) * 2.) * .5;

/*     scale the self-field by half, such that it sums to one below */

		    if (i__ == k) {
			for (j = 1; j <= 3; ++j) {
			    fid[j - 1] *= .5;
			    fkd[j - 1] *= .5;
			}
		    }
		    for (j = 1; j <= 3; ++j) {
			gkfield_ref(j, i__) = gkfield_ref(j, i__) + fid[j - 1]
				;
			gkfield_ref(j, k) = gkfield_ref(j, k) + fkd[j - 1];
			gkfieldp_ref(j, i__) = gkfieldp_ref(j, i__) + fid[j - 
				1];
			gkfieldp_ref(j, k) = gkfieldp_ref(j, k) + fkd[j - 1];
		    }
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i12_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i13_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i14_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i15_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = 1.;
	}
    }

/*     set vacuum induced dipoles to polarizability times direct field; */
/*     set SCRF induced dipoles to polarizability times direct field */
/*     plus the GK reaction field due to permanent multipoles */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    udir_ref(j, i__) = polar_1.polarity[i__ - 1] * field_ref(j, i__);
	    udirp_ref(j, i__) = polar_1.polarity[i__ - 1] * fieldp_ref(j, i__)
		    ;
	    udirs_ref(j, i__) = polar_1.polarity[i__ - 1] * (field_ref(j, i__)
		     + gkfield_ref(j, i__));
	    udirps_ref(j, i__) = polar_1.polarity[i__ - 1] * (fieldp_ref(j, 
		    i__) + gkfieldp_ref(j, i__));
	    uind_ref(j, i__) = udir_ref(j, i__);
	    uinp_ref(j, i__) = udirp_ref(j, i__);
	    uinds_ref(j, i__) = udirs_ref(j, i__);
	    uinps_ref(j, i__) = udirps_ref(j, i__);
	}
    }

/*     set tolerances for computation of mutual induced dipoles */

    if (s_cmp(polpot_1.poltyp, "MUTUAL", (ftnlen)6, (ftnlen)6) == 0) {
	done = FALSE_;
	maxiter = 500;
	iter = 0;
	eps = 100.;

/*     compute mutual induced dipole moments by an iterative method */

	while(! done) {
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    field_ref(j, i__) = 0.;
		    fieldp_ref(j, i__) = 0.;
		    fields_ref(j, i__) = 0.;
		    fieldps_ref(j, i__) = 0.;
		    gkfield_ref(j, i__) = 0.;
		    gkfieldp_ref(j, i__) = 0.;
		}
	    }
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ii = mpole_1.ipole[i__ - 1];
		pdi = polar_1.pdamp[i__ - 1];
		pti = polar_1.thole[i__ - 1];
		rbi = solute_1.rborn[ii - 1];
		duix = uind_ref(1, i__);
		duiy = uind_ref(2, i__);
		duiz = uind_ref(3, i__);
		puix = uinp_ref(1, i__);
		puiy = uinp_ref(2, i__);
		puiz = uinp_ref(3, i__);
		duixs = uinds_ref(1, i__);
		duiys = uinds_ref(2, i__);
		duizs = uinds_ref(3, i__);
		puixs = uinps_ref(1, i__);
		puiys = uinps_ref(2, i__);
		puizs = uinps_ref(3, i__);
		i__2 = polgrp_1.np11[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip11_ref(j, ii) - 1] = polpot_1.u1scale;
		}
		i__2 = polgrp_1.np12[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip12_ref(j, ii) - 1] = polpot_1.u2scale;
		}
		i__2 = polgrp_1.np13[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip13_ref(j, ii) - 1] = polpot_1.u3scale;
		}
		i__2 = polgrp_1.np14[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip14_ref(j, ii) - 1] = polpot_1.u4scale;
		}
		i__2 = mpole_1.npole;
		for (k = i__; k <= i__2; ++k) {
		    kk = mpole_1.ipole[k - 1];
		    rbk = solute_1.rborn[kk - 1];
		    proceed = TRUE_;
		    if (group_1.use_intra__) {
			groups_(&proceed, &fgrp, &ii, &kk, &c__0, &c__0, &
				c__0, &c__0);
		    }
		    if (proceed) {
			xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
			yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
			zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
			xr2 = xr * xr;
			yr2 = yr * yr;
			zr2 = zr * zr;
			r2 = xr2 + yr2 + zr2;
			if (r2 <= shunt_1.off2) {
			    r__ = sqrt(r2);
			    dukx = uind_ref(1, k);
			    duky = uind_ref(2, k);
			    dukz = uind_ref(3, k);
			    pukx = uinp_ref(1, k);
			    puky = uinp_ref(2, k);
			    pukz = uinp_ref(3, k);
			    dukxs = uinds_ref(1, k);
			    dukys = uinds_ref(2, k);
			    dukzs = uinds_ref(3, k);
			    pukxs = uinps_ref(1, k);
			    pukys = uinps_ref(2, k);
			    pukzs = uinps_ref(3, k);
			    if (i__ != k) {
				scale3 = dscale[kk - 1];
				scale5 = dscale[kk - 1];
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
					scale3 *= 1. - expdamp;
					scale5 *= 1. - (1. - damp) * expdamp;
				    }
				}
				rr3 = scale3 / (r__ * r2);
				rr5 = scale5 * 3. / (r__ * r2 * r2);
				duir = xr * duix + yr * duiy + zr * duiz;
				dukr = xr * dukx + yr * duky + zr * dukz;
				puir = xr * puix + yr * puiy + zr * puiz;
				pukr = xr * pukx + yr * puky + zr * pukz;
				duirs = xr * duixs + yr * duiys + zr * duizs;
				dukrs = xr * dukxs + yr * dukys + zr * dukzs;
				puirs = xr * puixs + yr * puiys + zr * puizs;
				pukrs = xr * pukxs + yr * pukys + zr * pukzs;
				fid[0] = -rr3 * dukx + rr5 * dukr * xr;
				fid[1] = -rr3 * duky + rr5 * dukr * yr;
				fid[2] = -rr3 * dukz + rr5 * dukr * zr;
				fkd[0] = -rr3 * duix + rr5 * duir * xr;
				fkd[1] = -rr3 * duiy + rr5 * duir * yr;
				fkd[2] = -rr3 * duiz + rr5 * duir * zr;
				fip[0] = -rr3 * pukx + rr5 * pukr * xr;
				fip[1] = -rr3 * puky + rr5 * pukr * yr;
				fip[2] = -rr3 * pukz + rr5 * pukr * zr;
				fkp[0] = -rr3 * puix + rr5 * puir * xr;
				fkp[1] = -rr3 * puiy + rr5 * puir * yr;
				fkp[2] = -rr3 * puiz + rr5 * puir * zr;
				fids[0] = -rr3 * dukxs + rr5 * dukrs * xr;
				fids[1] = -rr3 * dukys + rr5 * dukrs * yr;
				fids[2] = -rr3 * dukzs + rr5 * dukrs * zr;
				fkds[0] = -rr3 * duixs + rr5 * duirs * xr;
				fkds[1] = -rr3 * duiys + rr5 * duirs * yr;
				fkds[2] = -rr3 * duizs + rr5 * duirs * zr;
				fips[0] = -rr3 * pukxs + rr5 * pukrs * xr;
				fips[1] = -rr3 * pukys + rr5 * pukrs * yr;
				fips[2] = -rr3 * pukzs + rr5 * pukrs * zr;
				fkps[0] = -rr3 * puixs + rr5 * puirs * xr;
				fkps[1] = -rr3 * puiys + rr5 * puirs * yr;
				fkps[2] = -rr3 * puizs + rr5 * puirs * zr;
				for (j = 1; j <= 3; ++j) {
				    field_ref(j, i__) = field_ref(j, i__) + 
					    fid[j - 1];
				    field_ref(j, k) = field_ref(j, k) + fkd[j 
					    - 1];
				    fieldp_ref(j, i__) = fieldp_ref(j, i__) + 
					    fip[j - 1];
				    fieldp_ref(j, k) = fieldp_ref(j, k) + fkp[
					    j - 1];
				    fields_ref(j, i__) = fields_ref(j, i__) + 
					    fids[j - 1];
				    fields_ref(j, k) = fields_ref(j, k) + 
					    fkds[j - 1];
				    fieldps_ref(j, i__) = fieldps_ref(j, i__) 
					    + fips[j - 1];
				    fieldps_ref(j, k) = fieldps_ref(j, k) + 
					    fkps[j - 1];
				}
			    }
			    rb2 = rbi * rbk;
			    expterm = exp(-r2 / (gk_1.gkc * rb2));
			    expc = expterm / gk_1.gkc;
			    dexpc = -2. / (gk_1.gkc * rbi * rbk);
			    gf2 = 1. / (r2 + rb2 * expterm);
			    gf = sqrt(gf2);
			    gf3 = gf2 * gf;
			    gf5 = gf3 * gf2;

/*     reaction potential auxiliary terms */

			    a_ref(1, 0) = -gf3;
			    a_ref(2, 0) = gf5 * 3.;

/*     reaction potential gradient auxiliary terms */

			    expc1 = 1. - expc;
			    a_ref(1, 1) = expc1 * a_ref(2, 0);

/*     unweighted dipole reaction potential gradient tensor */

			    gux[1] = fd * (a_ref(1, 0) + xr2 * a_ref(1, 1));
			    gux[2] = fd * xr * yr * a_ref(1, 1);
			    gux[3] = fd * xr * zr * a_ref(1, 1);
			    guy[1] = gux[2];
			    guy[2] = fd * (a_ref(1, 0) + yr2 * a_ref(1, 1));
			    guy[3] = fd * yr * zr * a_ref(1, 1);
			    guz[1] = gux[3];
			    guz[2] = guy[3];
			    guz[3] = fd * (a_ref(1, 0) + zr2 * a_ref(1, 1));
			    fids[0] = dukxs * gux[1] + dukys * guy[1] + dukzs 
				    * guz[1];
			    fids[1] = dukxs * gux[2] + dukys * guy[2] + dukzs 
				    * guz[2];
			    fids[2] = dukxs * gux[3] + dukys * guy[3] + dukzs 
				    * guz[3];
			    fkds[0] = duixs * gux[1] + duiys * guy[1] + duizs 
				    * guz[1];
			    fkds[1] = duixs * gux[2] + duiys * guy[2] + duizs 
				    * guz[2];
			    fkds[2] = duixs * gux[3] + duiys * guy[3] + duizs 
				    * guz[3];
			    fips[0] = pukxs * gux[1] + pukys * guy[1] + pukzs 
				    * guz[1];
			    fips[1] = pukxs * gux[2] + pukys * guy[2] + pukzs 
				    * guz[2];
			    fips[2] = pukxs * gux[3] + pukys * guy[3] + pukzs 
				    * guz[3];
			    fkps[0] = puixs * gux[1] + puiys * guy[1] + puizs 
				    * guz[1];
			    fkps[1] = puixs * gux[2] + puiys * guy[2] + puizs 
				    * guz[2];
			    fkps[2] = puixs * gux[3] + puiys * guy[3] + puizs 
				    * guz[3];
			    if (i__ == k) {
				for (j = 1; j <= 3; ++j) {
				    fids[j - 1] *= .5;
				    fkds[j - 1] *= .5;
				    fips[j - 1] *= .5;
				    fkps[j - 1] *= .5;
				}
			    }
			    for (j = 1; j <= 3; ++j) {
				gkfield_ref(j, i__) = gkfield_ref(j, i__) + 
					fids[j - 1];
				gkfield_ref(j, k) = gkfield_ref(j, k) + fkds[
					j - 1];
				gkfieldp_ref(j, i__) = gkfieldp_ref(j, i__) + 
					fips[j - 1];
				gkfieldp_ref(j, k) = gkfieldp_ref(j, k) + 
					fkps[j - 1];
			    }
			}
		    }
		}

/*     reset interaction scaling coefficients for connected atoms */

		i__2 = polgrp_1.np11[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip11_ref(j, ii) - 1] = 1.;
		}
		i__2 = polgrp_1.np12[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip12_ref(j, ii) - 1] = 1.;
		}
		i__2 = polgrp_1.np13[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip13_ref(j, ii) - 1] = 1.;
		}
		i__2 = polgrp_1.np14[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip14_ref(j, ii) - 1] = 1.;
		}
	    }

/*     check to see if the mutual induced dipoles have converged */

	    ++iter;
	    epsold = eps;
	    epsd = 0.;
	    epsp = 0.;
	    epsds = 0.;
	    epsps = 0.;
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    uold_ref(j, i__) = uind_ref(j, i__);
		    uoldp_ref(j, i__) = uinp_ref(j, i__);
		    uolds_ref(j, i__) = uinds_ref(j, i__);
		    uoldps_ref(j, i__) = uinps_ref(j, i__);
		    uind_ref(j, i__) = udir_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * field_ref(j, i__);
		    uinp_ref(j, i__) = udirp_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * fieldp_ref(j, i__);
		    uinds_ref(j, i__) = udirs_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * (fields_ref(j, i__) + gkfield_ref(j, 
			    i__));
		    uinps_ref(j, i__) = udirps_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * (fieldps_ref(j, i__) + gkfieldp_ref(j, 
			    i__));
		    uind_ref(j, i__) = uold_ref(j, i__) + polpot_1.polsor * (
			    uind_ref(j, i__) - uold_ref(j, i__));
		    uinp_ref(j, i__) = uoldp_ref(j, i__) + polpot_1.polsor * (
			    uinp_ref(j, i__) - uoldp_ref(j, i__));
		    uinds_ref(j, i__) = uolds_ref(j, i__) + polpot_1.polsor * 
			    (uinds_ref(j, i__) - uolds_ref(j, i__));
		    uinps_ref(j, i__) = uoldps_ref(j, i__) + polpot_1.polsor *
			     (uinps_ref(j, i__) - uoldps_ref(j, i__));
/* Computing 2nd power */
		    d__1 = uind_ref(j, i__) - uold_ref(j, i__);
		    epsd += d__1 * d__1;
/* Computing 2nd power */
		    d__1 = uinp_ref(j, i__) - uoldp_ref(j, i__);
		    epsp += d__1 * d__1;
/* Computing 2nd power */
		    d__1 = uinds_ref(j, i__) - uolds_ref(j, i__);
		    epsds += d__1 * d__1;
/* Computing 2nd power */
		    d__1 = uinps_ref(j, i__) - uoldps_ref(j, i__);
		    epsps += d__1 * d__1;
		}
	    }
/* Computing MAX */
	    d__1 = max(epsd,epsp), d__1 = max(d__1,epsds);
	    eps = max(d__1,epsps);
	    eps = sqrt(eps / (doublereal) polar_1.npolar) * 4.80321;
	    if (inform_1.debug) {
		if (iter == 1) {
		    io___692.ciunit = iounit_1.iout;
		    s_wsfe(&io___692);
		    e_wsfe();
		}
		io___693.ciunit = iounit_1.iout;
		s_wsfe(&io___693);
		do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (eps < polpot_1.poleps) {
		done = TRUE_;
	    }
	    if (eps > epsold) {
		done = TRUE_;
	    }
	    if (iter >= maxiter) {
		done = TRUE_;
	    }
	}

/*     terminate the calculation if dipoles failed to converge */

	if (eps > polpot_1.poleps) {
	    io___694.ciunit = iounit_1.iout;
	    s_wsfe(&io___694);
	    e_wsfe();
	    prterr_();
	    fatal_();
	}
    }
    return 0;
} /* induce0e_ */

#undef fieldps_ref
#undef gkfield_ref
#undef uoldps_ref
#undef udirps_ref
#undef fields_ref
#undef fieldp_ref
#undef uolds_ref
#undef uoldp_ref
#undef udirs_ref
#undef udirp_ref
#undef uinps_ref
#undef uinds_ref
#undef rpole_ref
#undef field_ref
#undef uold_ref
#undef udir_ref
#undef uinp_ref
#undef uind_ref
#undef ip14_ref
#undef ip13_ref
#undef ip12_ref
#undef ip11_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref
#undef a_ref
#undef gkfieldp_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine induce0f  --  Poisson-Boltzmann induced dipoles  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "induce0f" computes the induced dipole moments at polarizable */
/*     sites for Poisson-Boltzmann SCRF and vacuum environments */


/* Subroutine */ int induce0f_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Determination of Induced Dipole\002,\002"
	    " Moments :\002,//,4x,\002Iter\002,8x,\002RMS Change (Debyes)\002"
	    ",/)";
    static char fmt_20[] = "(i8,7x,f16.10)";
    static char fmt_30[] = "(/,\002 INDUCE  --  Warning, Induced Dipole"
	    "s\002,\002 are not Converged\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int pbempole_(void);
    static integer i__, j, k;
    static doublereal r__, r2, ci;
    static integer ii, kk;
    static doublereal ck, xr, yr, zr;
    extern /* Subroutine */ int apbsinduce_(doublereal *, doublereal *);
    static doublereal rr3, rr5, rr7, xr2, yr2, zr2, fid[3], fkd[3], dir, dkr, 
	    eps, pdi, fip[3], fkp[3], qir, qkr, pti, qix, qiy, qiz, qkx, uxi, 
	    uyi, uzi, uxk, uyk, uzk, qky, qkz, damp, epsd, fids[3], fgrp;
    static integer iter;
    static doublereal duix, duiy, duiz, puix, puiy, puiz, dukx, qxxi, qxyi, 
	    qxzi, qyyi, qyzi, qzzi, qxxk, qxyk, qxzk, qyyk, qyzk, qzzk, duky, 
	    dukz, pukx, puky, pukz, duir, puir, dukr, pukr, epsp, fkds[3], 
	    fips[3], fkps[3], udir[75000]	/* was [3][25000] */, uold[
	    75000]	/* was [3][25000] */;
    static logical done;
    extern /* Subroutine */ int apbsnlinduce_(doublereal *, doublereal *);
    static doublereal field[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int fatal_(void);
    static doublereal duirs, duixs, duiys, duizs, dukxs, dukys, puixs, puiys, 
	    puizs, dukzs, pukxs, pukys, pukzs, puirs, dukrs, pukrs, epsds, 
	    epsps, scale3, scale5, scale7, udirp[75000]	/* was [3][25000] */, 
	    udirs[75000]	/* was [3][25000] */, uoldp[75000]	/* 
	    was [3][25000] */, uolds[75000]	/* was [3][25000] */, dscale[
	    25000], pgamma, fieldp[75000]	/* was [3][25000] */, fields[
	    75000]	/* was [3][25000] */, pscale[25000], epsold, udirps[
	    75000]	/* was [3][25000] */, uoldps[75000]	/* was [3][
	    25000] */;
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *), prterr_(void);
    static logical proceed;
    static doublereal fieldps[75000]	/* was [3][25000] */, indpole[75000]	
	    /* was [3][25000] */, expdamp;
    static integer maxiter;
    static doublereal inppole[75000]	/* was [3][25000] */;

    /* Fortran I/O blocks */
    static cilist io___816 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___817 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___818 = { 0, 0, 0, fmt_30, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define ip12_ref(a_1,a_2) polgrp_1.ip12[(a_2)*50 + a_1 - 51]
#define ip13_ref(a_1,a_2) polgrp_1.ip13[(a_2)*50 + a_1 - 51]
#define ip14_ref(a_1,a_2) polgrp_1.ip14[(a_2)*50 + a_1 - 51]
#define pbep_ref(a_1,a_2) pb_1.pbep[(a_2)*3 + a_1 - 4]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define uinp_ref(a_1,a_2) polar_1.uinp[(a_2)*3 + a_1 - 4]
#define udir_ref(a_1,a_2) udir[(a_2)*3 + a_1 - 4]
#define uold_ref(a_1,a_2) uold[(a_2)*3 + a_1 - 4]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define uinds_ref(a_1,a_2) polar_1.uinds[(a_2)*3 + a_1 - 4]
#define uinps_ref(a_1,a_2) polar_1.uinps[(a_2)*3 + a_1 - 4]
#define udirp_ref(a_1,a_2) udirp[(a_2)*3 + a_1 - 4]
#define udirs_ref(a_1,a_2) udirs[(a_2)*3 + a_1 - 4]
#define uoldp_ref(a_1,a_2) uoldp[(a_2)*3 + a_1 - 4]
#define uolds_ref(a_1,a_2) uolds[(a_2)*3 + a_1 - 4]
#define fieldp_ref(a_1,a_2) fieldp[(a_2)*3 + a_1 - 4]
#define fields_ref(a_1,a_2) fields[(a_2)*3 + a_1 - 4]
#define udirps_ref(a_1,a_2) udirps[(a_2)*3 + a_1 - 4]
#define uoldps_ref(a_1,a_2) uoldps[(a_2)*3 + a_1 - 4]
#define pbeuind_ref(a_1,a_2) pb_1.pbeuind[(a_2)*3 + a_1 - 4]
#define fieldps_ref(a_1,a_2) fieldps[(a_2)*3 + a_1 - 4]
#define indpole_ref(a_1,a_2) indpole[(a_2)*3 + a_1 - 4]
#define pbeuinp_ref(a_1,a_2) pb_1.pbeuinp[(a_2)*3 + a_1 - 4]
#define inppole_ref(a_1,a_2) inppole[(a_2)*3 + a_1 - 4]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  inform.i  --  control values for I/O and program flow  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     digits    decimal places output for energy and coordinates */
/*     iprint    steps between status printing (0=no printing) */
/*     iwrite    steps between coordinate dumps (0=no dumps) */
/*     isend     steps between socket communication (0=no sockets) */
/*     verbose   logical flag to turn on extra information */
/*     debug     logical flag to turn on full debug printing */
/*     holdup    logical flag to wait for carriage return on exit */
/*     abort     logical flag to stop execution at next chance */




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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  pb.i  --  parameters for Poisson-Boltzmann solvation  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     pbe      Poisson-Boltzman permanent multipole solvation energy */
/*     apbe     Poisson-Boltzman permanent multipole energy over atoms */
/*     pbr      Poisson-Boltzman cavity radii for atom types */
/*     pbep     Poisson-Boltzman energies on permanent multipoles */
/*     pbfp     Poisson-Boltzman forces on permanent multipoles */
/*     pbtp     Poisson-Boltzman torques on permanent multipoles */
/*     pbeuind  Poisson-Boltzman field due to induced dipoles */
/*     pbeuinp  Poisson-Boltzman field due to non-local induced dipoles */

/*     APBS configuration parameters (see APBS documentation for more details) */
/*     In the column on the right are possible values for each variable, with */
/*     the default values given in brackets. Note that only a subset of APBS */
/*     options are supported and/or are appropriate for use with AMOEBA. */

/*     pbtyp                                   lpbe */

/*     At some point AMOEBA with the non-linear PBE could be supported, but */
/*     there is only have theory for energies (no gradients). */

/*     pbsoln                                  mg-auto, [mg-manual] */

/*     Currently there is only limited support for focusing calculations, */
/*     which is a powerful feature of APBS. The current requirement is */
/*     that energies and forces must all be calculated using the finest */
/*     solution. */

/*     bcfl     boundary conditions            zero, sdh, [mdh] */
/*     chgm     multipole discretization       spl4 */

/*     other charge discretization methods are not appropriate for AMOEBA */

/*     srfm     surface method                 mol, smol, [spl4] */

/*     spl4 is required for forces calculations, although mol is useful for */
/*     comparison with generalized Kirkwood */

/*     dime     number of grid points          [65, 65, 65] */
/*     grid     grid spacing (mg-manual)       fxn of "dime" */
/*     cgrid    coarse grid spacing            fxn of "dime" */
/*     fgrid    fine grid spacing              cgrid / 2 */

/*     stable results require grid spacing to be fine enough to keep */
/*     multipoles inside the dielectric boundary (2.5 * grid < PBR) */

/*     gcent    grid center  (mg-manual)       center of mass */
/*     cgcent   coarse grid center             center of mass */
/*     fgcent   fine grid center               center of mass */
/*     pdie     solute/homogeneous dieletric   [1.0] */
/*     sdie     solvent dieletric              [78.3] */
/*     ionn     number of ion species          [0] */
/*     ionc     ion concentration (M)          [0.0] */
/*     ionq     ion charge (electrons)         [1.0] */
/*     ionr     ion radius (A)                 [2.0] */
/*     srad     solvent probe radius (A)       [1.4] */
/*     swin     surface spline window width    [0.3] */
/*     sdens    density of surface points      [10.0] */

/*     additional parameter to facilitate default grid setup */

/*     smin     minimum distance between an    [10.0] */
/*              atomic center and the grid */
/*              boundary (A) */




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  gk.i  --  parameters for generalized Kirkwood solvation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     gkr      generalized Kirkwood cavity radii for atom types */
/*     gkc      tuning parameter exponent in the f(GB) function */




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




/*     zero out the induced dipoles and the fields at each site; */
/*     uind & uinp are vacuum dipoles, uinds & uinps are SCRF dipoles, */
/*     and field & fieldp are solute fields */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    uind_ref(j, i__) = 0.;
	    uinp_ref(j, i__) = 0.;
	    uinds_ref(j, i__) = 0.;
	    uinps_ref(j, i__) = 0.;
	    field_ref(j, i__) = 0.;
	    fieldp_ref(j, i__) = 0.;
	}
    }
    if (! potent_1.use_polar__ || ! potent_1.use_solv__) {
	return 0;
    }

/*     set the switching function coefficients */

    switch_("MPOLE", (ftnlen)5);

/*     set arrays needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pscale[i__ - 1] = 1.;
	dscale[i__ - 1] = 1.;
    }

/*     compute the direct induced dipole moment at each atom, and */
/*     another set that also includes RF due to permanent multipoles */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	pdi = polar_1.pdamp[i__ - 1];
	pti = polar_1.thole[i__ - 1];
	ci = rpole_ref(1, i__);
	uxi = rpole_ref(2, i__);
	uyi = rpole_ref(3, i__);
	uzi = rpole_ref(4, i__);
	qxxi = rpole_ref(5, i__);
	qxyi = rpole_ref(6, i__);
	qxzi = rpole_ref(7, i__);
	qyyi = rpole_ref(9, i__);
	qyzi = rpole_ref(10, i__);
	qzzi = rpole_ref(13, i__);
	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
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
	    pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = polpot_1.d1scale;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = polpot_1.d2scale;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = polpot_1.d3scale;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = polpot_1.d4scale;
	}
	i__2 = mpole_1.npole;
	for (k = i__ + 1; k <= i__2; ++k) {
	    kk = mpole_1.ipole[k - 1];
	    proceed = TRUE_;
	    if (group_1.use_intra__) {
		groups_(&proceed, &fgrp, &ii, &kk, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
		yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
		zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
		xr2 = xr * xr;
		yr2 = yr * yr;
		zr2 = zr * zr;
		r2 = xr2 + yr2 + zr2;
		if (r2 <= shunt_1.off2) {
		    r__ = sqrt(r2);
		    ck = rpole_ref(1, k);
		    uxk = rpole_ref(2, k);
		    uyk = rpole_ref(3, k);
		    uzk = rpole_ref(4, k);
		    qxxk = rpole_ref(5, k);
		    qxyk = rpole_ref(6, k);
		    qxzk = rpole_ref(7, k);
		    qyyk = rpole_ref(9, k);
		    qyzk = rpole_ref(10, k);
		    qzzk = rpole_ref(13, k);

/*     self-interactions for the solute field are skipped */

		    scale3 = 1.;
		    scale5 = 1.;
		    scale7 = 1.;
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
			    scale5 = 1. - expdamp * (1. - damp);
/* Computing 2nd power */
			    d__1 = damp;
			    scale7 = 1. - expdamp * (1. - damp + d__1 * d__1 *
				     .6);
			}
		    }
		    rr3 = scale3 / (r__ * r2);
		    rr5 = scale5 * 3. / (r__ * r2 * r2);
		    rr7 = scale7 * 15. / (r__ * r2 * r2 * r2);
		    dir = uxi * xr + uyi * yr + uzi * zr;
		    qix = qxxi * xr + qxyi * yr + qxzi * zr;
		    qiy = qxyi * xr + qyyi * yr + qyzi * zr;
		    qiz = qxzi * xr + qyzi * yr + qzzi * zr;
		    qir = qix * xr + qiy * yr + qiz * zr;
		    dkr = uxk * xr + uyk * yr + uzk * zr;
		    qkx = qxxk * xr + qxyk * yr + qxzk * zr;
		    qky = qxyk * xr + qyyk * yr + qyzk * zr;
		    qkz = qxzk * xr + qyzk * yr + qzzk * zr;
		    qkr = qkx * xr + qky * yr + qkz * zr;
		    fid[0] = -xr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * 
			    uxk + rr5 * 2. * qkx;
		    fid[1] = -yr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * 
			    uyk + rr5 * 2. * qky;
		    fid[2] = -zr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * 
			    uzk + rr5 * 2. * qkz;
		    fkd[0] = xr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * 
			    uxi - rr5 * 2. * qix;
		    fkd[1] = yr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * 
			    uyi - rr5 * 2. * qiy;
		    fkd[2] = zr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * 
			    uzi - rr5 * 2. * qiz;
		    for (j = 1; j <= 3; ++j) {
			field_ref(j, i__) = field_ref(j, i__) + fid[j - 1] * 
				dscale[kk - 1];
			field_ref(j, k) = field_ref(j, k) + fkd[j - 1] * 
				dscale[kk - 1];
			fieldp_ref(j, i__) = fieldp_ref(j, i__) + fid[j - 1] *
				 pscale[kk - 1];
			fieldp_ref(j, k) = fieldp_ref(j, k) + fkd[j - 1] * 
				pscale[kk - 1];
		    }
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i12_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i13_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i14_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i15_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = 1.;
	}
    }

/*     find Poisson-Boltzmann reaction field due to permanent multipoles */

    pbempole_();

/*     set vacuum induced dipoles to polarizability times direct field; */
/*     SCRF induced dipoles are polarizability times direct field */
/*     plus the reaction field due to permanent multipoles */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	for (j = 1; j <= 3; ++j) {
	    udir_ref(j, i__) = polar_1.polarity[i__ - 1] * field_ref(j, i__);
	    udirp_ref(j, i__) = polar_1.polarity[i__ - 1] * fieldp_ref(j, i__)
		    ;
	    udirs_ref(j, i__) = polar_1.polarity[i__ - 1] * (field_ref(j, i__)
		     + pbep_ref(j, ii));
	    udirps_ref(j, i__) = polar_1.polarity[i__ - 1] * (fieldp_ref(j, 
		    i__) + pbep_ref(j, ii));
	    uind_ref(j, i__) = udir_ref(j, i__);
	    uinp_ref(j, i__) = udirp_ref(j, i__);
	    uinds_ref(j, i__) = udirs_ref(j, i__);
	    uinps_ref(j, i__) = udirps_ref(j, i__);
	}
    }

/*     set tolerances for computation of mutual induced dipoles */

    if (s_cmp(polpot_1.poltyp, "MUTUAL", (ftnlen)6, (ftnlen)6) == 0) {
	done = FALSE_;
	maxiter = 500;
	iter = 0;
	eps = 100.;

/*     get mutual induced dipole moments by an iterative method */

	while(! done) {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    indpole_ref(j, i__) = 0.;
		    inppole_ref(j, i__) = 0.;
		    pbeuind_ref(j, i__) = 0.;
		    pbeuinp_ref(j, i__) = 0.;
		}
	    }
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ii = mpole_1.ipole[i__ - 1];
		for (j = 1; j <= 3; ++j) {
		    indpole_ref(j, ii) = uinds_ref(j, i__);
		    inppole_ref(j, ii) = uinps_ref(j, i__);
		}
	    }
	    apbsinduce_(indpole, pb_1.pbeuind);
	    apbsnlinduce_(inppole, pb_1.pbeuinp);

/*     compute intra-solute mutual polarization */

	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    field_ref(j, i__) = 0.;
		    fieldp_ref(j, i__) = 0.;
		    fields_ref(j, i__) = 0.;
		    fieldps_ref(j, i__) = 0.;
		}
	    }
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ii = mpole_1.ipole[i__ - 1];
		pdi = polar_1.pdamp[i__ - 1];
		pti = polar_1.thole[i__ - 1];
		duix = uind_ref(1, i__);
		duiy = uind_ref(2, i__);
		duiz = uind_ref(3, i__);
		puix = uinp_ref(1, i__);
		puiy = uinp_ref(2, i__);
		puiz = uinp_ref(3, i__);
		duixs = uinds_ref(1, i__);
		duiys = uinds_ref(2, i__);
		duizs = uinds_ref(3, i__);
		puixs = uinps_ref(1, i__);
		puiys = uinps_ref(2, i__);
		puizs = uinps_ref(3, i__);
		i__2 = polgrp_1.np11[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip11_ref(j, ii) - 1] = polpot_1.u1scale;
		}
		i__2 = polgrp_1.np12[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip12_ref(j, ii) - 1] = polpot_1.u2scale;
		}
		i__2 = polgrp_1.np13[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip13_ref(j, ii) - 1] = polpot_1.u3scale;
		}
		i__2 = polgrp_1.np14[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip14_ref(j, ii) - 1] = polpot_1.u4scale;
		}
		i__2 = mpole_1.npole;
		for (k = i__ + 1; k <= i__2; ++k) {
		    kk = mpole_1.ipole[k - 1];
		    proceed = TRUE_;
		    if (group_1.use_intra__) {
			groups_(&proceed, &fgrp, &ii, &kk, &c__0, &c__0, &
				c__0, &c__0);
		    }
		    if (proceed) {
			xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
			yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
			zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
			xr2 = xr * xr;
			yr2 = yr * yr;
			zr2 = zr * zr;
			r2 = xr2 + yr2 + zr2;
			if (r2 <= shunt_1.off2) {
			    r__ = sqrt(r2);
			    dukx = uind_ref(1, k);
			    duky = uind_ref(2, k);
			    dukz = uind_ref(3, k);
			    pukx = uinp_ref(1, k);
			    puky = uinp_ref(2, k);
			    pukz = uinp_ref(3, k);
			    dukxs = uinds_ref(1, k);
			    dukys = uinds_ref(2, k);
			    dukzs = uinds_ref(3, k);
			    pukxs = uinps_ref(1, k);
			    pukys = uinps_ref(2, k);
			    pukzs = uinps_ref(3, k);
			    scale3 = dscale[kk - 1];
			    scale5 = dscale[kk - 1];
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
				    scale3 *= 1. - expdamp;
				    scale5 *= 1. - (1. - damp) * expdamp;
				}
			    }
			    rr3 = scale3 / (r__ * r2);
			    rr5 = scale5 * 3. / (r__ * r2 * r2);
			    duir = xr * duix + yr * duiy + zr * duiz;
			    dukr = xr * dukx + yr * duky + zr * dukz;
			    puir = xr * puix + yr * puiy + zr * puiz;
			    pukr = xr * pukx + yr * puky + zr * pukz;
			    duirs = xr * duixs + yr * duiys + zr * duizs;
			    dukrs = xr * dukxs + yr * dukys + zr * dukzs;
			    puirs = xr * puixs + yr * puiys + zr * puizs;
			    pukrs = xr * pukxs + yr * pukys + zr * pukzs;
			    fid[0] = -rr3 * dukx + rr5 * dukr * xr;
			    fid[1] = -rr3 * duky + rr5 * dukr * yr;
			    fid[2] = -rr3 * dukz + rr5 * dukr * zr;
			    fkd[0] = -rr3 * duix + rr5 * duir * xr;
			    fkd[1] = -rr3 * duiy + rr5 * duir * yr;
			    fkd[2] = -rr3 * duiz + rr5 * duir * zr;
			    fip[0] = -rr3 * pukx + rr5 * pukr * xr;
			    fip[1] = -rr3 * puky + rr5 * pukr * yr;
			    fip[2] = -rr3 * pukz + rr5 * pukr * zr;
			    fkp[0] = -rr3 * puix + rr5 * puir * xr;
			    fkp[1] = -rr3 * puiy + rr5 * puir * yr;
			    fkp[2] = -rr3 * puiz + rr5 * puir * zr;
			    fids[0] = -rr3 * dukxs + rr5 * dukrs * xr;
			    fids[1] = -rr3 * dukys + rr5 * dukrs * yr;
			    fids[2] = -rr3 * dukzs + rr5 * dukrs * zr;
			    fkds[0] = -rr3 * duixs + rr5 * duirs * xr;
			    fkds[1] = -rr3 * duiys + rr5 * duirs * yr;
			    fkds[2] = -rr3 * duizs + rr5 * duirs * zr;
			    fips[0] = -rr3 * pukxs + rr5 * pukrs * xr;
			    fips[1] = -rr3 * pukys + rr5 * pukrs * yr;
			    fips[2] = -rr3 * pukzs + rr5 * pukrs * zr;
			    fkps[0] = -rr3 * puixs + rr5 * puirs * xr;
			    fkps[1] = -rr3 * puiys + rr5 * puirs * yr;
			    fkps[2] = -rr3 * puizs + rr5 * puirs * zr;
			    for (j = 1; j <= 3; ++j) {
				field_ref(j, i__) = field_ref(j, i__) + fid[j 
					- 1];
				field_ref(j, k) = field_ref(j, k) + fkd[j - 1]
					;
				fieldp_ref(j, i__) = fieldp_ref(j, i__) + fip[
					j - 1];
				fieldp_ref(j, k) = fieldp_ref(j, k) + fkp[j - 
					1];
				fields_ref(j, i__) = fields_ref(j, i__) + 
					fids[j - 1];
				fields_ref(j, k) = fields_ref(j, k) + fkds[j 
					- 1];
				fieldps_ref(j, i__) = fieldps_ref(j, i__) + 
					fips[j - 1];
				fieldps_ref(j, k) = fieldps_ref(j, k) + fkps[
					j - 1];
			    }
			}
		    }
		}

/*     reset interaction scaling coefficients for connected atoms */

		i__2 = polgrp_1.np11[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip11_ref(j, ii) - 1] = 1.;
		}
		i__2 = polgrp_1.np12[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip12_ref(j, ii) - 1] = 1.;
		}
		i__2 = polgrp_1.np13[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip13_ref(j, ii) - 1] = 1.;
		}
		i__2 = polgrp_1.np14[ii - 1];
		for (j = 1; j <= i__2; ++j) {
		    dscale[ip14_ref(j, ii) - 1] = 1.;
		}
	    }

/*     check to see if the mutual induced dipoles have converged */

	    ++iter;
	    epsold = eps;
	    epsd = 0.;
	    epsp = 0.;
	    epsds = 0.;
	    epsps = 0.;
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ii = mpole_1.ipole[i__ - 1];
		for (j = 1; j <= 3; ++j) {
		    uold_ref(j, i__) = uind_ref(j, i__);
		    uoldp_ref(j, i__) = uinp_ref(j, i__);
		    uolds_ref(j, i__) = uinds_ref(j, i__);
		    uoldps_ref(j, i__) = uinps_ref(j, i__);
		    uind_ref(j, i__) = udir_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * field_ref(j, i__);
		    uinp_ref(j, i__) = udirp_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * fieldp_ref(j, i__);
		    uinds_ref(j, i__) = udirs_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * (fields_ref(j, i__) + pbeuind_ref(j, 
			    ii));
		    uinps_ref(j, i__) = udirps_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * (fieldps_ref(j, i__) + pbeuinp_ref(j, 
			    ii));
		    uind_ref(j, i__) = uold_ref(j, i__) + polpot_1.polsor * (
			    uind_ref(j, i__) - uold_ref(j, i__));
		    uinp_ref(j, i__) = uoldp_ref(j, i__) + polpot_1.polsor * (
			    uinp_ref(j, i__) - uoldp_ref(j, i__));
		    uinds_ref(j, i__) = uolds_ref(j, i__) + polpot_1.polsor * 
			    (uinds_ref(j, i__) - uolds_ref(j, i__));
		    uinps_ref(j, i__) = uoldps_ref(j, i__) + polpot_1.polsor *
			     (uinps_ref(j, i__) - uoldps_ref(j, i__));
/* Computing 2nd power */
		    d__1 = uind_ref(j, i__) - uold_ref(j, i__);
		    epsd += d__1 * d__1;
/* Computing 2nd power */
		    d__1 = uinp_ref(j, i__) - uoldp_ref(j, i__);
		    epsp += d__1 * d__1;
/* Computing 2nd power */
		    d__1 = uinds_ref(j, i__) - uolds_ref(j, i__);
		    epsds += d__1 * d__1;
/* Computing 2nd power */
		    d__1 = uinps_ref(j, i__) - uoldps_ref(j, i__);
		    epsps += d__1 * d__1;
		}
	    }
/* Computing MAX */
	    d__1 = max(epsd,epsp), d__1 = max(d__1,epsds);
	    eps = max(d__1,epsps);
	    eps = sqrt(eps / (doublereal) polar_1.npolar) * 4.80321;
	    if (inform_1.debug) {
		if (iter == 1) {
		    io___816.ciunit = iounit_1.iout;
		    s_wsfe(&io___816);
		    e_wsfe();
		}
		io___817.ciunit = iounit_1.iout;
		s_wsfe(&io___817);
		do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (eps < polpot_1.poleps) {
		done = TRUE_;
	    }
	    if (eps > epsold) {
		done = TRUE_;
	    }
	    if (iter >= maxiter) {
		done = TRUE_;
	    }
	}

/*     terminate the calculation if dipoles failed to converge */

	if (eps > polpot_1.poleps) {
	    io___818.ciunit = iounit_1.iout;
	    s_wsfe(&io___818);
	    e_wsfe();
	    prterr_();
	    fatal_();
	}
    }
    return 0;
} /* induce0f_ */

#undef inppole_ref
#undef pbeuinp_ref
#undef indpole_ref
#undef fieldps_ref
#undef pbeuind_ref
#undef uoldps_ref
#undef udirps_ref
#undef fields_ref
#undef fieldp_ref
#undef uolds_ref
#undef uoldp_ref
#undef udirs_ref
#undef udirp_ref
#undef uinps_ref
#undef uinds_ref
#undef rpole_ref
#undef field_ref
#undef uold_ref
#undef udir_ref
#undef uinp_ref
#undef uind_ref
#undef pbep_ref
#undef ip14_ref
#undef ip13_ref
#undef ip12_ref
#undef ip11_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref


