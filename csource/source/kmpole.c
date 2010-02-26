/* kmpole.f -- translated by f2c (version 20050501).
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
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

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
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

struct {
    doublereal multip[26000]	/* was [13][2000] */;
    char mpaxis[16000], kmp[32000];
} kmulti_;

#define kmulti_1 kmulti_

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
    integer np11[25000], ip11[2500000]	/* was [100][25000] */, np12[25000], 
	    ip12[1250000]	/* was [50][25000] */, np13[25000], ip13[
	    1250000]	/* was [50][25000] */, np14[25000], ip14[1250000]	
	    /* was [50][25000] */;
} polgrp_;

#define polgrp_1 polgrp_

struct {
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__4 = 4;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine kmpole  --  multipole parameter assignment  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "kmpole" assigns atomic multipole moments to the atoms of */
/*     the structure and processes any new or changed values */


/* Subroutine */ int kmpole_(void)
{
    /* Format strings */
    static char fmt_40[] = "(/,\002 KMPOLE  --  Too many Atomic Multipole"
	    "\002,\002 Parameters\002)";
    static char fmt_70[] = "(/,\002 Additional Atomic Multipole Parameters "
	    ":\002,//,5x,\002Atom Type\002,5x,\002Coordinate Frame\002,\002 D"
	    "efinition\002,8x,\002Multipole Moments\002)";
    static char fmt_80[] = "(/,4x,i6,5x,i6,1x,i6,1x,i6,3x,a8,2x,f9.5,/,48x,3"
	    "f9.5,/,48x,f9.5,/,48x,2f9.5,/,48x,3f9.5)";
    static char fmt_130[] = "(/,\002 Additional Atomic Multipoles\002,\002 f"
	    "or Specific Atoms :\002,//,6x,\002Atom\002,9x,\002Coordinate Fra"
	    "me\002,\002 Definition\002,8x,\002Multipole Moments\002)";
    static char fmt_140[] = "(/,4x,i6,5x,i6,1x,i6,1x,i6,3x,a8,2x,f9.5,/,48x,"
	    "3f9.5,/,48x,f9.5,/,48x,2f9.5,/,48x,3f9.5)";

    /* System generated locals */
    address a__1[4];
    integer i__1, i__2, i__3[4], i__4, i__5, i__6, i__7;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i__, j, k, l, m, ji, ki, li;
    static char pa[4], pb[4], pc[4], pd[4];
    static integer it, jt, kt, lt, kx, ky, kz;
    static char pt[16];
    static integer big, imp;
    static doublereal mpl[13];
    static integer nmp;
    static char axt[8];
    static integer mpt[2000], mpx[2000], mpy[2000], mpz[2000];
    static logical path;
    static integer size, next, xtyp, ytyp, ztyp;
    static char blank[16];
    static logical header;
    static char record[120];
    extern doublereal random_(void);
    extern integer number_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), chkpole_(void), 
	    numeral_(integer *, char *, integer *, ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___10 = { 1, string, 1, 0, 120, 1 };
    static icilist io___15 = { 1, string, 1, 0, 120, 1 };
    static icilist io___16 = { 1, record, 1, 0, 120, 1 };
    static icilist io___17 = { 1, record, 1, 0, 120, 1 };
    static icilist io___18 = { 1, record, 1, 0, 120, 1 };
    static icilist io___19 = { 1, record, 1, 0, 120, 1 };
    static cilist io___20 = { 0, 0, 0, fmt_40, 0 };
    static icilist io___25 = { 1, string, 1, 0, 120, 1 };
    static icilist io___26 = { 1, string, 1, 0, 120, 1 };
    static icilist io___27 = { 1, record, 1, 0, 120, 1 };
    static icilist io___28 = { 1, record, 1, 0, 120, 1 };
    static icilist io___29 = { 1, record, 1, 0, 120, 1 };
    static icilist io___30 = { 1, record, 1, 0, 120, 1 };
    static cilist io___31 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_80, 0 };
    static icilist io___55 = { 1, string, 1, 0, 120, 1 };
    static icilist io___56 = { 1, string, 1, 0, 120, 1 };
    static icilist io___57 = { 1, record, 1, 0, 120, 1 };
    static icilist io___58 = { 1, record, 1, 0, 120, 1 };
    static icilist io___59 = { 1, record, 1, 0, 120, 1 };
    static icilist io___60 = { 1, record, 1, 0, 120, 1 };
    static cilist io___61 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_140, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define kmp_ref(a_0,a_1) &kmulti_1.kmp[(a_1)*16 + a_0 - 16]
#define pole_ref(a_1,a_2) mpole_1.pole[(a_2)*13 + a_1 - 14]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define uinp_ref(a_1,a_2) polar_1.uinp[(a_2)*3 + a_1 - 4]
#define uinds_ref(a_1,a_2) polar_1.uinds[(a_2)*3 + a_1 - 4]
#define uinps_ref(a_1,a_2) polar_1.uinps[(a_2)*3 + a_1 - 4]
#define polaxe_ref(a_0,a_1) &mpole_1.polaxe[(a_1)*8 + a_0 - 8]
#define mpaxis_ref(a_0,a_1) &kmulti_1.mpaxis[(a_1)*8 + a_0 - 8]
#define multip_ref(a_1,a_2) kmulti_1.multip[(a_2)*13 + a_1 - 14]
#define keyline_ref(a_0,a_1) &keys_1.keyline[(a_1)*120 + a_0 - 120]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  keys.i  --  contents of current keyword parameter file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nkey      number of nonblank lines in the keyword file */
/*     keyline   contents of each individual keyword file line */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kmulti.i  --  forcefield parameters for atomic multipoles  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxnmp   maximum number of atomic multipole parameter entries */

/*     multip   atomic monopole, dipole and quadrupole values */
/*     mpaxis   type of local axis definition for atomic multipoles */
/*     kmp      string of atom types for atomic multipoles */




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




/*     count the number of existing multipole parameters */

    s_copy(blank, "                ", (ftnlen)16, (ftnlen)16);
    nmp = 2000;
    for (i__ = 2000; i__ >= 1; --i__) {
	if (s_cmp(kmp_ref(0, i__), blank, (ftnlen)16, (ftnlen)16) == 0) {
	    nmp = i__ - 1;
	}
    }

/*     find and count new multipole parameters in the keyfile */

    imp = 0;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "MULTIPOLE ", (ftnlen)10, (ftnlen)10) == 0) {
	    k = 0;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___10);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&kz, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&kx, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ky, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&mpl[0], (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	    goto L20;
L10:
	    i__2 = s_rsli(&io___15);
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&kz, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&kx, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&mpl[0], (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L30;
	    }
L20:
	    if (k > 0) {
		s_copy(record, keyline_ref(0, i__ + 1), (ftnlen)120, (ftnlen)
			120);
		i__2 = s_rsli(&io___16);
		if (i__2 != 0) {
		    goto L30;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[1], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L30;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[2], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L30;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[3], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L30;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L30;
		}
		s_copy(record, keyline_ref(0, i__ + 2), (ftnlen)120, (ftnlen)
			120);
		i__2 = s_rsli(&io___17);
		if (i__2 != 0) {
		    goto L30;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[4], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L30;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L30;
		}
		s_copy(record, keyline_ref(0, i__ + 3), (ftnlen)120, (ftnlen)
			120);
		i__2 = s_rsli(&io___18);
		if (i__2 != 0) {
		    goto L30;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[7], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L30;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[8], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L30;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L30;
		}
		s_copy(record, keyline_ref(0, i__ + 4), (ftnlen)120, (ftnlen)
			120);
		i__2 = s_rsli(&io___19);
		if (i__2 != 0) {
		    goto L30;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[10], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L30;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[11], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L30;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[12], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L30;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L30;
		}
		++imp;
	    }
L30:
	    ;
	}
    }

/*     check for too many combined parameter values */

    nmp += imp;
    if (nmp > 2000) {
	io___20.ciunit = iounit_1.iout;
	s_wsfe(&io___20);
	e_wsfe();
	inform_1.abort = TRUE_;
    }

/*     move existing parameters to make room for new values */

    if (imp != 0) {
	i__1 = imp + 1;
	for (j = nmp; j >= i__1; --j) {
	    k = j - imp;
	    s_copy(kmp_ref(0, j), kmp_ref(0, k), (ftnlen)16, (ftnlen)16);
	    s_copy(mpaxis_ref(0, j), mpaxis_ref(0, k), (ftnlen)8, (ftnlen)8);
	    for (m = 1; m <= 13; ++m) {
		multip_ref(m, j) = multip_ref(m, k);
	    }
	}
    }

/*     process keywords containing atomic multipole parameters */

    imp = 0;
    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "MULTIPOLE ", (ftnlen)10, (ftnlen)10) == 0) {
	    k = 0;
	    kz = 0;
	    kx = 0;
	    ky = 0;
	    s_copy(axt, "Z-then-X", (ftnlen)8, (ftnlen)8);
	    for (j = 1; j <= 13; ++j) {
		mpl[j - 1] = 0.;
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___25);
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&kz, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&kx, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ky, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&mpl[0], (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L50;
	    }
	    goto L60;
L50:
	    ky = 0;
	    i__2 = s_rsli(&io___26);
	    if (i__2 != 0) {
		goto L90;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L90;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&kz, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L90;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&kx, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L90;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&mpl[0], (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L90;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L90;
	    }
L60:
	    if (k > 0) {
		if (kz < 0 || kx < 0) {
		    s_copy(axt, "Bisector", (ftnlen)8, (ftnlen)8);
		}
		if (kx < 0 && ky < 0) {
		    s_copy(axt, "Z-Bisect", (ftnlen)8, (ftnlen)8);
		}
/* Computing MAX */
		i__2 = max(kz,kx);
		if (max(i__2,ky) < 0) {
		    s_copy(axt, "3-Fold", (ftnlen)8, (ftnlen)6);
		}
		kz = abs(kz);
		kx = abs(kx);
		ky = abs(ky);
		s_copy(record, keyline_ref(0, i__ + 1), (ftnlen)120, (ftnlen)
			120);
		i__2 = s_rsli(&io___27);
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[1], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[2], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[3], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L90;
		}
		s_copy(record, keyline_ref(0, i__ + 2), (ftnlen)120, (ftnlen)
			120);
		i__2 = s_rsli(&io___28);
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[4], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L90;
		}
		s_copy(record, keyline_ref(0, i__ + 3), (ftnlen)120, (ftnlen)
			120);
		i__2 = s_rsli(&io___29);
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[7], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[8], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L90;
		}
		s_copy(record, keyline_ref(0, i__ + 4), (ftnlen)120, (ftnlen)
			120);
		i__2 = s_rsli(&io___30);
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[10], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[11], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[12], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L90;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L90;
		}
		mpl[5] = mpl[7];
		mpl[6] = mpl[10];
		mpl[9] = mpl[11];
		if (header) {
		    header = FALSE_;
		    io___31.ciunit = iounit_1.iout;
		    s_wsfe(&io___31);
		    e_wsfe();
		}
		io___32.ciunit = iounit_1.iout;
		s_wsfe(&io___32);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kz, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kx, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ky, (ftnlen)sizeof(integer));
		do_fio(&c__1, axt, (ftnlen)8);
		for (j = 1; j <= 5; ++j) {
		    do_fio(&c__1, (char *)&mpl[j - 1], (ftnlen)sizeof(
			    doublereal));
		}
		do_fio(&c__1, (char *)&mpl[7], (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&mpl[8], (ftnlen)sizeof(doublereal));
		for (j = 11; j <= 13; ++j) {
		    do_fio(&c__1, (char *)&mpl[j - 1], (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
		size = 4;
		numeral_(&k, pa, &size, (ftnlen)4);
		numeral_(&kz, pb, &size, (ftnlen)4);
		numeral_(&kx, pc, &size, (ftnlen)4);
		numeral_(&ky, pd, &size, (ftnlen)4);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
		++imp;
		s_copy(kmp_ref(0, imp), pt, (ftnlen)16, (ftnlen)16);
		s_copy(mpaxis_ref(0, imp), axt, (ftnlen)8, (ftnlen)8);
		for (j = 1; j <= 13; ++j) {
		    multip_ref(j, imp) = mpl[j - 1];
		}
	    }
L90:
	    ;
	}
    }

/*     zero out local axes, multipoles and polarization attachments */

    mpole_1.npole = atoms_1.n;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mpole_1.pollist[i__ - 1] = 0;
	mpole_1.zaxis[i__ - 1] = 0;
	mpole_1.xaxis[i__ - 1] = 0;
	mpole_1.yaxis[i__ - 1] = 0;
	s_copy(polaxe_ref(0, i__), "        ", (ftnlen)8, (ftnlen)8);
	for (j = 1; j <= 13; ++j) {
	    pole_ref(j, i__) = 0.;
	}
	polgrp_1.np11[i__ - 1] = 0;
	polgrp_1.np12[i__ - 1] = 0;
	polgrp_1.np13[i__ - 1] = 0;
	polgrp_1.np14[i__ - 1] = 0;
    }

/*     store the atom types associated with each parameter */

    i__1 = nmp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mpt[i__ - 1] = number_(kmp_ref(0, i__), (ftnlen)4);
	mpz[i__ - 1] = number_(kmp_ref(4, i__), (ftnlen)4);
	mpx[i__ - 1] = number_(kmp_ref(8, i__), (ftnlen)4);
	mpy[i__ - 1] = number_(kmp_ref(12, i__), (ftnlen)4);
    }

/*     assign multipole parameters via only 1-2 connected atoms */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	it = atoms_1.type__[i__ - 1];
	i__2 = nmp;
	for (imp = 1; imp <= i__2; ++imp) {
	    if (it == mpt[imp - 1]) {
		ztyp = mpz[imp - 1];
		xtyp = mpx[imp - 1];
		ytyp = mpy[imp - 1];
		i__4 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__4; ++j) {
		    ji = i12_ref(j, i__);
		    jt = atoms_1.type__[ji - 1];
		    if (jt == ztyp) {
			i__5 = couple_1.n12[i__ - 1];
			for (k = 1; k <= i__5; ++k) {
			    ki = i12_ref(k, i__);
			    kt = atoms_1.type__[ki - 1];
			    if (kt == xtyp && ki != ji) {
				if (ytyp == 0) {
				    mpole_1.zaxis[i__ - 1] = ji;
				    mpole_1.xaxis[i__ - 1] = ki;
				    s_copy(polaxe_ref(0, i__), mpaxis_ref(0, 
					    imp), (ftnlen)8, (ftnlen)8);
				    for (m = 1; m <= 13; ++m) {
					pole_ref(m, i__) = multip_ref(m, imp);
				    }
				    goto L100;
				}
				i__6 = couple_1.n12[i__ - 1];
				for (l = 1; l <= i__6; ++l) {
				    li = i12_ref(l, i__);
				    lt = atoms_1.type__[li - 1];
				    if (lt == ytyp && li != ji && li != ki) {
					mpole_1.zaxis[i__ - 1] = ji;
					mpole_1.xaxis[i__ - 1] = ki;
					mpole_1.yaxis[i__ - 1] = li;
					s_copy(polaxe_ref(0, i__), mpaxis_ref(
						0, imp), (ftnlen)8, (ftnlen)8)
						;
					for (m = 1; m <= 13; ++m) {
					    pole_ref(m, i__) = multip_ref(m, 
						    imp);
					}
					goto L100;
				    }
				}
			    }
			}
		    }
		}
	    }
	}

/*     assign multipole parameters via 1-2 and 1-3 connected atoms */

	i__2 = nmp;
	for (imp = 1; imp <= i__2; ++imp) {
	    if (it == mpt[imp - 1]) {
		ztyp = mpz[imp - 1];
		xtyp = mpx[imp - 1];
		ytyp = mpy[imp - 1];
		i__4 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__4; ++j) {
		    ji = i12_ref(j, i__);
		    jt = atoms_1.type__[ji - 1];
		    if (jt == ztyp) {
			i__5 = couple_1.n13[i__ - 1];
			for (k = 1; k <= i__5; ++k) {
			    ki = i13_ref(k, i__);
			    kt = atoms_1.type__[ki - 1];
			    path = FALSE_;
			    i__6 = couple_1.n12[ki - 1];
			    for (m = 1; m <= i__6; ++m) {
				if (i12_ref(m, ki) == ji) {
				    path = TRUE_;
				}
			    }
			    if (kt == xtyp && path) {
				if (ytyp == 0) {
				    mpole_1.zaxis[i__ - 1] = ji;
				    mpole_1.xaxis[i__ - 1] = ki;
				    s_copy(polaxe_ref(0, i__), mpaxis_ref(0, 
					    imp), (ftnlen)8, (ftnlen)8);
				    for (m = 1; m <= 13; ++m) {
					pole_ref(m, i__) = multip_ref(m, imp);
				    }
				    goto L100;
				}
				i__6 = couple_1.n13[i__ - 1];
				for (l = 1; l <= i__6; ++l) {
				    li = i13_ref(l, i__);
				    lt = atoms_1.type__[li - 1];
				    path = FALSE_;
				    i__7 = couple_1.n12[li - 1];
				    for (m = 1; m <= i__7; ++m) {
					if (i12_ref(m, li) == ji) {
					    path = TRUE_;
					}
				    }
				    if (lt == ytyp && li != ki && path) {
					mpole_1.zaxis[i__ - 1] = ji;
					mpole_1.xaxis[i__ - 1] = ki;
					mpole_1.yaxis[i__ - 1] = li;
					s_copy(polaxe_ref(0, i__), mpaxis_ref(
						0, imp), (ftnlen)8, (ftnlen)8)
						;
					for (m = 1; m <= 13; ++m) {
					    pole_ref(m, i__) = multip_ref(m, 
						    imp);
					}
					goto L100;
				    }
				}
			    }
			}
		    }
		}
	    }
	}

/*     assign multipole parameters via only a z-defining atom */

	i__2 = nmp;
	for (imp = 1; imp <= i__2; ++imp) {
	    if (it == mpt[imp - 1]) {
		ztyp = mpz[imp - 1];
		xtyp = mpx[imp - 1];
		ytyp = mpy[imp - 1];
		i__4 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__4; ++j) {
		    ji = i12_ref(j, i__);
		    jt = atoms_1.type__[ji - 1];
		    if (jt == ztyp) {
			if (xtyp == 0) {
			    mpole_1.zaxis[i__ - 1] = ji;
			    mpole_1.xaxis[i__ - 1] = atoms_1.n + 1;
			    s_copy(polaxe_ref(0, i__), mpaxis_ref(0, imp), (
				    ftnlen)8, (ftnlen)8);
			    for (m = 1; m <= 13; ++m) {
				pole_ref(m, i__) = multip_ref(m, imp);
			    }
			    goto L100;
			}
		    }
		}
	    }
	}

/*     assign multipole parameters via no connected atoms */

	i__2 = nmp;
	for (imp = 1; imp <= i__2; ++imp) {
	    if (it == mpt[imp - 1]) {
		ztyp = mpz[imp - 1];
		xtyp = mpx[imp - 1];
		ytyp = mpy[imp - 1];
		if (ztyp == 0) {
		    mpole_1.zaxis[i__ - 1] = atoms_1.n + 1;
		    mpole_1.xaxis[i__ - 1] = atoms_1.n + 2;
		    s_copy(polaxe_ref(0, i__), mpaxis_ref(0, imp), (ftnlen)8, 
			    (ftnlen)8);
		    for (m = 1; m <= 13; ++m) {
			pole_ref(m, i__) = multip_ref(m, imp);
		    }
		    goto L100;
		}
	    }
	}
L100:
	;
    }

/*     process keywords with multipole parameters for specific atoms */

    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "MULTIPOLE ", (ftnlen)10, (ftnlen)10) == 0) {
	    k = 0;
	    kz = 0;
	    kx = 0;
	    ky = 0;
	    s_copy(axt, "Z-then-X", (ftnlen)8, (ftnlen)8);
	    for (j = 1; j <= 13; ++j) {
		mpl[j - 1] = 0.;
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___55);
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&kz, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&kx, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ky, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&mpl[0], (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L110;
	    }
	    goto L120;
L110:
	    ky = 0;
	    i__2 = s_rsli(&io___56);
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&kz, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&kx, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&mpl[0], (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L150;
	    }
L120:
	    if (k < 0 && k >= -atoms_1.n) {
		k = -k;
		if (kz < 0 || kx < 0) {
		    s_copy(axt, "Bisector", (ftnlen)8, (ftnlen)8);
		}
		if (kx < 0 && ky < 0) {
		    s_copy(axt, "Z-Bisect", (ftnlen)8, (ftnlen)8);
		}
/* Computing MAX */
		i__2 = max(kz,kx);
		if (max(i__2,ky) < 0) {
		    s_copy(axt, "3-Fold", (ftnlen)8, (ftnlen)6);
		}
		kz = abs(kz);
		kx = abs(kx);
		ky = abs(ky);
		s_copy(record, keyline_ref(0, i__ + 1), (ftnlen)120, (ftnlen)
			120);
		i__2 = s_rsli(&io___57);
		if (i__2 != 0) {
		    goto L150;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[1], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L150;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[2], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L150;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[3], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L150;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L150;
		}
		s_copy(record, keyline_ref(0, i__ + 2), (ftnlen)120, (ftnlen)
			120);
		i__2 = s_rsli(&io___58);
		if (i__2 != 0) {
		    goto L150;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[4], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L150;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L150;
		}
		s_copy(record, keyline_ref(0, i__ + 3), (ftnlen)120, (ftnlen)
			120);
		i__2 = s_rsli(&io___59);
		if (i__2 != 0) {
		    goto L150;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[7], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L150;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[8], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L150;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L150;
		}
		s_copy(record, keyline_ref(0, i__ + 4), (ftnlen)120, (ftnlen)
			120);
		i__2 = s_rsli(&io___60);
		if (i__2 != 0) {
		    goto L150;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[10], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L150;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[11], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L150;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mpl[12], (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L150;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L150;
		}
		mpl[5] = mpl[7];
		mpl[6] = mpl[10];
		mpl[9] = mpl[11];
		if (header) {
		    header = FALSE_;
		    io___61.ciunit = iounit_1.iout;
		    s_wsfe(&io___61);
		    e_wsfe();
		}
		io___62.ciunit = iounit_1.iout;
		s_wsfe(&io___62);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kz, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kx, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ky, (ftnlen)sizeof(integer));
		do_fio(&c__1, axt, (ftnlen)8);
		for (j = 1; j <= 5; ++j) {
		    do_fio(&c__1, (char *)&mpl[j - 1], (ftnlen)sizeof(
			    doublereal));
		}
		do_fio(&c__1, (char *)&mpl[7], (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&mpl[8], (ftnlen)sizeof(doublereal));
		for (j = 11; j <= 13; ++j) {
		    do_fio(&c__1, (char *)&mpl[j - 1], (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
		if (kz == 0) {
		    kz = atoms_1.n + 1;
		}
		if (kx == 0) {
		    kx = atoms_1.n + 2;
		}
		mpole_1.zaxis[k - 1] = kz;
		mpole_1.xaxis[k - 1] = kx;
		mpole_1.yaxis[k - 1] = ky;
		s_copy(polaxe_ref(0, k), axt, (ftnlen)8, (ftnlen)8);
		for (j = 1; j <= 13; ++j) {
		    pole_ref(j, k) = mpl[j - 1];
		}
	    }
L150:
	    ;
	}
    }

/*     convert the dipole and quadrupole moments to Angstroms, */
/*     quadrupole divided by 3 for use as traceless values */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (k = 2; k <= 4; ++k) {
	    pole_ref(k, i__) = pole_ref(k, i__) * .52917720859;
	}
	for (k = 5; k <= 13; ++k) {
	    pole_ref(k, i__) = pole_ref(k, i__) * .28002851809110429 / 3.;
	}
    }

/*     get the order of the multipole expansion at each site */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	size = 0;
	for (k = 1; k <= 13; ++k) {
	    if (pole_ref(k, i__) != 0.) {
		size = max(k,size);
	    }
	}
	if (size > 4) {
	    size = 13;
	} else if (size > 1) {
	    size = 4;
	}
	mpole_1.polsiz[i__ - 1] = size;
    }

/*     if needed, get random coordinates for dummy axis defining atoms */

    big = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	i__2 = big, i__4 = mpole_1.zaxis[i__ - 1], i__2 = max(i__2,i__4), 
		i__4 = mpole_1.xaxis[i__ - 1], i__2 = max(i__2,i__4), i__4 = 
		mpole_1.yaxis[i__ - 1];
	big = max(i__2,i__4);
    }
    if (big > atoms_1.n) {
	i__1 = big;
	for (i__ = atoms_1.n + 1; i__ <= i__1; ++i__) {
	    atoms_1.x[i__ - 1] = random_();
	    atoms_1.y[i__ - 1] = random_();
	    atoms_1.z__[i__ - 1] = random_();
	}
    }

/*     if polarization not used, zero out induced dipoles */

    if (! potent_1.use_polar__) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		uind_ref(j, i__) = 0.;
		uinp_ref(j, i__) = 0.;
		uinds_ref(j, i__) = 0.;
		uinps_ref(j, i__) = 0.;
	    }
	}
    }

/*     remove any zero or undefined atomic multipoles */

    if (! potent_1.use_polar__ && ! potent_1.use_solv__) {
	mpole_1.npole = 0;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (mpole_1.polsiz[i__ - 1] != 0) {
		++mpole_1.npole;
		mpole_1.ipole[mpole_1.npole - 1] = i__;
		mpole_1.pollist[i__ - 1] = mpole_1.npole;
		mpole_1.zaxis[mpole_1.npole - 1] = mpole_1.zaxis[i__ - 1];
		mpole_1.xaxis[mpole_1.npole - 1] = mpole_1.xaxis[i__ - 1];
		mpole_1.yaxis[mpole_1.npole - 1] = mpole_1.yaxis[i__ - 1];
		s_copy(polaxe_ref(0, mpole_1.npole), polaxe_ref(0, i__), (
			ftnlen)8, (ftnlen)8);
		for (j = 1; j <= 13; ++j) {
		    pole_ref(j, mpole_1.npole) = pole_ref(j, i__);
		}
	    }
	}

/*     test multipoles at chiral sites and invert if necessary */

	chkpole_();

/*     turn off the atomic multipole potential if it is not used */

	if (mpole_1.npole == 0) {
	    potent_1.use_mpole__ = FALSE_;
	}
    }
    return 0;
} /* kmpole_ */

#undef keyline_ref
#undef multip_ref
#undef mpaxis_ref
#undef polaxe_ref
#undef uinps_ref
#undef uinds_ref
#undef uinp_ref
#undef uind_ref
#undef pole_ref
#undef kmp_ref
#undef i13_ref
#undef i12_ref


