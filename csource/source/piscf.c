/* piscf.f -- translated by f2c (version 20050501).
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
    doublereal bk[50000], bl[50000];
    integer nbond, ibnd[100000]	/* was [2][50000] */;
} bond_;

#define bond_1 bond_

struct {
    doublereal pbpl[200], pnpl[200];
} border_;

#define border_1 border_

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
    doublereal q[100], w[100], em[100];
    integer nfill;
} orbits_;

#define orbits_1 orbits_

struct {
    integer norbit, iorbit[100], reorbit, piperp[300]	/* was [3][100] */, 
	    nbpi, ibpi[600]	/* was [3][200] */, ntpi, itpi[800]	/* 
	    was [2][400] */;
    logical listpi[25000];
} piorbs_;

#define piorbs_1 piorbs_

struct {
    doublereal bkpi[200], blpi[200], kslope[200], lslope[200], torsp2[400];
} pistuf_;

#define pistuf_1 pistuf_

struct {
    doublereal tors1[400000]	/* was [4][100000] */, tors2[400000]	/* 
	    was [4][100000] */, tors3[400000]	/* was [4][100000] */, tors4[
	    400000]	/* was [4][100000] */, tors5[400000]	/* was [4][
	    100000] */, tors6[400000]	/* was [4][100000] */;
    integer ntors, itors[400000]	/* was [4][100000] */;
} tors_;

#define tors_1 tors_

/* Table of constant values */

static integer c__6 = 6;
static integer c__100 = 100;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine piscf  --  scf molecular orbital calculation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "piscf" performs an scf molecular orbital calculation for */
/*     the pisystem using a modified Pariser-Parr-Pople method */


/* Subroutine */ int piscf_(void)
{
    /* Initialized data */

    static logical first = TRUE_;
    static integer ncalls = 0;

    /* Format strings */
    static char fmt_10[] = "(\002 PISCF  --  The SCF-MO Iteration has\002"
	    ",\002 not reached Self-Consistency\002)";
    static char fmt_20[] = "(/,\002 Pi-SCF-MO Calculation for Planar Syste"
	    "m :\002)";
    static char fmt_30[] = "(/,\002 Pi-SCF-MO Calculation for Non-Plana"
	    "r\002,\002 System :\002)";
    static char fmt_40[] = "(/,\002 Total Energy\002,11x,f12.4,/,\002 Conver"
	    "gence\002,12x,d12.4,/,\002 Iterations\002,13x,i12)";
    static char fmt_50[] = "(/,\002 Core Integrals\002,9x,f12.4,/,\002 Coulo"
	    "mb Repulsion\002,6x,f12.4,/,\002 Exchange Repulsion\002,5x,f12.4"
	    ",/,\002 Nuclear Repulsion\002,6x,f12.4)";
    static char fmt_60[] = "(/,\002 Orbital Energies\002)";
    static char fmt_70[] = "(8f9.4)";
    static char fmt_80[] = "(/,\002 Molecular Orbitals\002)";
    static char fmt_90[] = "(8f9.4)";
    static char fmt_100[] = "(/,\002 Fock Matrix\002)";
    static char fmt_110[] = "(8f9.4)";
    static char fmt_120[] = "(/,\002 Electron Densities\002)";
    static char fmt_130[] = "(8f9.4)";
    static char fmt_140[] = "(/,\002 Density Matrix\002)";
    static char fmt_150[] = "(8f9.4)";
    static char fmt_160[] = "(/,\002 H-Core Matrix\002)";
    static char fmt_170[] = "(8f9.4)";
    static char fmt_180[] = "(/,\002 Gamma Matrix\002)";
    static char fmt_190[] = "(8f9.4)";
    static char fmt_200[] = "(/,\002 Pi Bond Orders\002)";
    static char fmt_210[] = "(4x,2i5,3x,f10.4)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal), exp(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal converge;
    static integer i__, j, k, m;
    static doublereal p, v[10000]	/* was [100][100] */, s1, s2, g11, 
	    g12, hc[10000]	/* was [100][100] */, g14, en[100], ip[100], 
	    qi, ed[10000]	/* was [100][100] */, xg, xi, xj, xk, ebb, 
	    ebe, blb, ble, gii, gij, gjk, rij, vij, vik, xij, yij, zij, vmj, 
	    vmk, hcii, hcij, aeth, fock[10000]	/* was [100][100] */;
    static char mode[6];
    static doublereal brij, erij;
    static integer iatn, jatn, iorb, jorb;
    static doublereal g11sq, abnz;
    static integer iter;
    static doublereal work1[100], work2[100], ebeta, gamma[10000]	/* 
	    was [100][100] */, delta, ovlap, total, rijsq;
    extern /* Subroutine */ int jacobi_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal bebond, eebond, covlap, totold, povlap[200];
    extern /* Subroutine */ int pitilt_(doublereal *);
    static doublereal cionize;
    extern /* Subroutine */ int pialter_(void);
    static integer maxiter;
    static doublereal iionize, jionize;
    extern /* Subroutine */ int overlap_(integer *, integer *, doublereal *, 
	    doublereal *);

    /* Fortran I/O blocks */
    static cilist io___68 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___71 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___72 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___73 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___75 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___76 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___77 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___78 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___79 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___80 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___81 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___82 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___83 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___85 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___86 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___87 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___88 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___90 = { 0, 0, 0, fmt_210, 0 };



#define v_ref(a_1,a_2) v[(a_2)*100 + a_1 - 101]
#define hc_ref(a_1,a_2) hc[(a_2)*100 + a_1 - 101]
#define ed_ref(a_1,a_2) ed[(a_2)*100 + a_1 - 101]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define fock_ref(a_1,a_2) fock[(a_2)*100 + a_1 - 101]
#define ibpi_ref(a_1,a_2) piorbs_1.ibpi[(a_2)*3 + a_1 - 4]
#define gamma_ref(a_1,a_2) gamma[(a_2)*100 + a_1 - 101]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  bond.i  --  covalent bonds in the current structure  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     bk      bond stretch force constants (kcal/mole/Ang**2) */
/*     bl      ideal bond length values in Angstroms */
/*     nbond   total number of bond stretches in the system */
/*     ibnd    numbers of the atoms in each bond stretch */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  border.i  --  bond orders for a conjugated pisystem  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     pbpl   pi-bond orders for bonds in "planar" pisystem */
/*     pnpl   pi-bond orders for bonds in "nonplanar" pisystem */




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
/*     ##  orbits.i  --  orbital energies for conjugated pisystem  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     q       number of pi-electrons contributed by each atom */
/*     w       ionization potential of each pisystem atom */
/*     em      repulsion integral for each pisystem atom */
/*     nfill   number of filled pisystem molecular orbitals */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  piorbs.i  --  conjugated system in the current structure  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     norbit    total number of pisystem orbitals in the system */
/*     iorbit    numbers of the atoms containing pisystem orbitals */
/*     reorbit   number of evaluations between orbital updates */
/*     piperp    atoms defining a normal plane to each orbital */
/*     nbpi      total number of bonds affected by the pisystem */
/*     ibpi      bond and piatom numbers for each pisystem bond */
/*     ntpi      total number of torsions affected by the pisystem */
/*     itpi      torsion and pibond numbers for each pisystem torsion */
/*     listpi    atom list indicating whether each atom has an orbital */




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




/*     only needs to be done if pisystem is present */

    if (piorbs_1.norbit == 0) {
	return 0;
    }

/*     increment the number of calls to this routine */

    ++ncalls;
    if (piorbs_1.reorbit == 0 || ncalls < piorbs_1.reorbit) {
	return 0;
    }
    ncalls = 0;

/*     initialize some constants and parameters: */

/*     mode      planar or nonplanar pi-calculation */
/*     maxiter   maximum number of scf iterations */
/*     converge  criterion for scf convergence */
/*     ebeta     value of resonance integral for ethylene */
/*     cionize   ionizaton potential of carbon (hartree) */

    s_copy(mode, "PLANAR", (ftnlen)6, (ftnlen)6);
    maxiter = 50;
    converge = 1e-8;
    ebeta = -.0757;
    cionize = -.41012247144125952;

/*     set the bond energies, alpha values and ideal bond length */
/*     parameter for carbon-carbon pibond type parameters: */

/*     ebe = equilibrium bond energy in ethylene */
/*     ebb = equilibrium bond energy in benzene */
/*     aeth = the P-P-P constant "a" in ethylene */
/*     abnz = the P-P-P constant "a" in benzene */
/*     ble = equilibrium bond length in ethylene */
/*     blb = equilibrium bond length in benzene */

    ebe = 129.37;
    ebb = 117.58;
    aeth = 2.309;
    abnz = 2.142;
    ble = 1.338;
    blb = 1.397;

/*     assign empirical one-center Coulomb integrals, and */
/*     first or second ionization potential depending on */
/*     whether the orbital contributes one or two electrons */

    i__1 = piorbs_1.norbit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gamma_ref(i__, i__) = orbits_1.em[i__ - 1];
	ip[i__ - 1] = orbits_1.w[i__ - 1] + (1. - orbits_1.q[i__ - 1]) * 
		orbits_1.em[i__ - 1];
    }

/*     calculate two-center repulsion integrals */
/*     according to Ohno's semi-empirical formula */

    i__1 = piorbs_1.norbit - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iorb = piorbs_1.iorbit[i__ - 1];
	gii = gamma_ref(i__, i__);
	i__2 = piorbs_1.norbit;
	for (j = i__ + 1; j <= i__2; ++j) {
	    jorb = piorbs_1.iorbit[j - 1];
	    g11 = (gii + gamma_ref(j, j)) * .5;
/* Computing 2nd power */
	    d__1 = g11;
	    g11sq = 1. / (d__1 * d__1);
	    xij = atoms_1.x[iorb - 1] - atoms_1.x[jorb - 1];
	    yij = atoms_1.y[iorb - 1] - atoms_1.y[jorb - 1];
	    zij = atoms_1.z__[iorb - 1] - atoms_1.z__[jorb - 1];
/* Computing 2nd power */
	    d__1 = xij;
/* Computing 2nd power */
	    d__2 = yij;
/* Computing 2nd power */
	    d__3 = zij;
	    rijsq = (d__1 * d__1 + d__2 * d__2 + d__3 * d__3) / 
		    .28002851809110429;
	    g12 = 1. / sqrt(rijsq + g11sq);
	    gamma_ref(i__, j) = g12;
	    gamma_ref(j, i__) = g12;
	}
    }

/*     zero out the resonance integral values */

    i__1 = piorbs_1.norbit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = piorbs_1.norbit;
	for (j = 1; j <= i__2; ++j) {
	    hc_ref(j, i__) = 0.;
	}
    }

/*     the first term in the sum to find alpha is the first */
/*     or second ionization potential, then the two-center */
/*     repulsion integrals are added */

    i__1 = piorbs_1.norbit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	hcii = ip[i__ - 1];
	i__2 = piorbs_1.norbit;
	for (j = 1; j <= i__2; ++j) {
	    if (i__ != j) {
		hcii -= orbits_1.q[j - 1] * gamma_ref(i__, j);
	    }
	}
	hc_ref(i__, i__) = hcii;
    }

/*     get two-center repulsion integrals via Ohno's formula */

    i__1 = piorbs_1.nbpi;
    for (k = 1; k <= i__1; ++k) {
	i__ = ibpi_ref(2, k);
	j = ibpi_ref(3, k);
	iorb = piorbs_1.iorbit[i__ - 1];
	jorb = piorbs_1.iorbit[j - 1];
	iatn = atmtyp_1.atomic[iorb - 1];
	jatn = atmtyp_1.atomic[jorb - 1];
	xij = atoms_1.x[iorb - 1] - atoms_1.x[jorb - 1];
	yij = atoms_1.y[iorb - 1] - atoms_1.y[jorb - 1];
	zij = atoms_1.z__[iorb - 1] - atoms_1.z__[jorb - 1];
/* Computing 2nd power */
	d__1 = xij;
/* Computing 2nd power */
	d__2 = yij;
/* Computing 2nd power */
	d__3 = zij;
	rij = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
	d__1 = rij;
	rijsq = d__1 * d__1 / .28002851809110429;
	g11 = (gamma_ref(i__, i__) + gamma_ref(j, j)) * .5;
/* Computing 2nd power */
	d__1 = g11;
	g11sq = 1. / (d__1 * d__1);
	g12 = gamma_ref(i__, j);

/*     compute the bond energy using a Morse potential */

	erij = aeth * (ble - rij);
	brij = abnz * (blb - rij);
	eebond = (exp(erij) * 2. - exp(erij * 2.)) * ebe / 627.5094688;
	bebond = (exp(brij) * 2. - exp(brij * 2.)) * ebb / 627.5094688;

/*     compute the carbon-carbon resonance integral using */
/*     the Whitehead and Lo formula */

	g14 = 1. / sqrt(rijsq * 4. + g11sq);
	hcij = (bebond - eebond) * 1.5 - g11 * .375 + g12 * 
		.41666666666666669 - g14 / 24.;

/*     if either atom is non-carbon, then factor the resonance */
/*     integral by overlap ratio and ionization potential ratio */

	if (iatn != 6 || jatn != 6) {
	    overlap_(&iatn, &jatn, &rij, &ovlap);
	    overlap_(&c__6, &c__6, &rij, &covlap);
	    hcij *= ovlap / covlap;
	    iionize = ip[i__ - 1];
	    if (orbits_1.q[i__ - 1] != 1.) {
		if (iatn == 7) {
		    iionize *= .595;
		}
		if (iatn == 8) {
		    iionize *= .525;
		}
		if (iatn == 16) {
		    iionize *= .89;
		}
	    }
	    jionize = ip[j - 1];
	    if (orbits_1.q[j - 1] != 1.) {
		if (jatn == 7) {
		    jionize *= .595;
		}
		if (jatn == 8) {
		    jionize *= .525;
		}
		if (jatn == 16) {
		    jionize *= .89;
		}
	    }
	    hcij = hcij * (iionize + jionize) / (cionize * 2.);
	}

/*     set symmetric elements to the same value */

	hc_ref(i__, j) = hcij;
	hc_ref(j, i__) = hcij;
    }

/*     make an initial guess at the Fock matrix if needed */

    if (first) {
	first = FALSE_;
	i__1 = piorbs_1.norbit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = piorbs_1.norbit;
	    for (j = 1; j <= i__2; ++j) {
		fock_ref(j, i__) = hc_ref(j, i__);
	    }
	}
	i__1 = piorbs_1.norbit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    fock_ref(i__, i__) = ip[i__ - 1] * .5;
	}
    }

/*     now, do the scf-mo computation; note that it needs to */
/*     be done twice, initially for the planar analog of the */
/*     actual system; then for the nonplanar (actual) system */

    while(s_cmp(mode, "PLANAR", (ftnlen)6, (ftnlen)6) == 0 || s_cmp(mode, 
	    "NONPLN", (ftnlen)6, (ftnlen)6) == 0) {
	if (s_cmp(mode, "NONPLN", (ftnlen)6, (ftnlen)6) == 0) {
	    pitilt_(povlap);
	    i__1 = piorbs_1.nbpi;
	    for (k = 1; k <= i__1; ++k) {
		i__ = ibpi_ref(2, k);
		j = ibpi_ref(3, k);
		hc_ref(i__, j) = hc_ref(i__, j) * povlap[k - 1];
		hc_ref(j, i__) = hc_ref(i__, j);
	    }
	}

/*     perform scf iterations until convergence is reached; */
/*     diagonalize the Fock matrix "f" to get the mo's, */
/*     then use mo's to form the next "f" matrix assuming */
/*     zero differential overlap except for one-center */
/*     exchange repulsions */

	iter = 0;
	delta = converge * 2.;
	while(delta > converge && iter < maxiter) {
	    ++iter;
	    jacobi_(&piorbs_1.norbit, &c__100, fock, en, v, work1, work2);
	    i__1 = piorbs_1.norbit;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = piorbs_1.norbit;
		for (j = i__; j <= i__2; ++j) {
		    s1 = 0.;
		    s2 = 0.;
		    gij = gamma_ref(i__, j);
		    i__3 = orbits_1.nfill;
		    for (k = 1; k <= i__3; ++k) {
			s2 -= v_ref(i__, k) * v_ref(j, k) * gij;
			if (i__ == j) {
			    i__4 = piorbs_1.norbit;
			    for (m = 1; m <= i__4; ++m) {
/* Computing 2nd power */
				d__1 = v_ref(m, k);
				s1 += gamma_ref(i__, m) * 2. * (d__1 * d__1);
			    }
			}
		    }
		    fock_ref(i__, j) = s1 + s2 + hc_ref(i__, j);
		    fock_ref(j, i__) = fock_ref(i__, j);
		}
	    }

/*     calculate the ground state energy, where "xi" sums the */
/*     molecular core integrals, "xj" sums the molecular coulomb */
/*     repulsion integrals, "xk" sums the molecular exchange */
/*     repulsion integrals, and "xg" is sums the nuclear repulsion */

	    xi = 0.;
	    xj = 0.;
	    xk = 0.;
	    xg = 0.;
	    i__1 = orbits_1.nfill;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = piorbs_1.norbit;
		for (j = 1; j <= i__2; ++j) {
		    vij = v_ref(j, i__);
		    i__3 = piorbs_1.norbit;
		    for (k = 1; k <= i__3; ++k) {
			vik = v_ref(k, i__);
			gjk = gamma_ref(j, k);
			xi += vij * 2. * vik * hc_ref(j, k);
			i__4 = orbits_1.nfill;
			for (m = 1; m <= i__4; ++m) {
			    vmj = v_ref(j, m);
			    vmk = v_ref(k, m);
			    xj += vij * 2. * vij * vmk * vmk * gjk;
			    xk -= vij * vmj * vik * vmk * gjk;
			}
		    }
		}
	    }
	    i__1 = piorbs_1.norbit - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		qi = orbits_1.q[i__ - 1];
		i__2 = piorbs_1.norbit;
		for (j = i__ + 1; j <= i__2; ++j) {
		    xg += qi * orbits_1.q[j - 1] * gamma_ref(i__, j);
		}
	    }
	    total = xi + xj + xk + xg;
	    if (iter != 1) {
		delta = (d__1 = total - totold, abs(d__1));
	    }
	    totold = total;
	}

/*     print warning if scf-mo iteration did not converge */

	if (delta > converge) {
	    io___68.ciunit = iounit_1.iout;
	    s_wsfe(&io___68);
	    e_wsfe();
	}

/*     calculate electron densities from filled mo's */

	i__1 = piorbs_1.norbit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = piorbs_1.norbit;
	    for (j = 1; j <= i__2; ++j) {
		ed_ref(i__, j) = 0.;
		i__3 = orbits_1.nfill;
		for (k = 1; k <= i__3; ++k) {
		    ed_ref(i__, j) = ed_ref(i__, j) + v_ref(i__, k) * 2. * 
			    v_ref(j, k);
		}
	    }
	}

/*     print out results for the scf computation */

	if (inform_1.debug) {
	    if (s_cmp(mode, "PLANAR", (ftnlen)6, (ftnlen)6) == 0) {
		io___70.ciunit = iounit_1.iout;
		s_wsfe(&io___70);
		e_wsfe();
	    } else {
		io___71.ciunit = iounit_1.iout;
		s_wsfe(&io___71);
		e_wsfe();
	    }
	    io___72.ciunit = iounit_1.iout;
	    s_wsfe(&io___72);
	    do_fio(&c__1, (char *)&total, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&delta, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___73.ciunit = iounit_1.iout;
	    s_wsfe(&io___73);
	    do_fio(&c__1, (char *)&xi, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&xj, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&xk, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&xg, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    io___74.ciunit = iounit_1.iout;
	    s_wsfe(&io___74);
	    e_wsfe();
	    io___75.ciunit = iounit_1.iout;
	    s_wsfe(&io___75);
	    i__1 = piorbs_1.norbit;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&en[i__ - 1], (ftnlen)sizeof(doublereal)
			);
	    }
	    e_wsfe();
	    io___76.ciunit = iounit_1.iout;
	    s_wsfe(&io___76);
	    e_wsfe();
	    i__1 = piorbs_1.norbit;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___77.ciunit = iounit_1.iout;
		s_wsfe(&io___77);
		i__2 = piorbs_1.norbit;
		for (j = 1; j <= i__2; ++j) {
		    do_fio(&c__1, (char *)&v_ref(i__, j), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	    io___78.ciunit = iounit_1.iout;
	    s_wsfe(&io___78);
	    e_wsfe();
	    i__1 = piorbs_1.norbit;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___79.ciunit = iounit_1.iout;
		s_wsfe(&io___79);
		i__2 = piorbs_1.norbit;
		for (j = 1; j <= i__2; ++j) {
		    do_fio(&c__1, (char *)&fock_ref(i__, j), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	    io___80.ciunit = iounit_1.iout;
	    s_wsfe(&io___80);
	    e_wsfe();
	    io___81.ciunit = iounit_1.iout;
	    s_wsfe(&io___81);
	    i__1 = piorbs_1.norbit;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&ed_ref(i__, i__), (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	    io___82.ciunit = iounit_1.iout;
	    s_wsfe(&io___82);
	    e_wsfe();
	    i__1 = piorbs_1.norbit;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___83.ciunit = iounit_1.iout;
		s_wsfe(&io___83);
		i__2 = piorbs_1.norbit;
		for (j = 1; j <= i__2; ++j) {
		    do_fio(&c__1, (char *)&ed_ref(i__, j), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	    io___84.ciunit = iounit_1.iout;
	    s_wsfe(&io___84);
	    e_wsfe();
	    i__1 = piorbs_1.norbit;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___85.ciunit = iounit_1.iout;
		s_wsfe(&io___85);
		i__2 = piorbs_1.norbit;
		for (j = 1; j <= i__2; ++j) {
		    do_fio(&c__1, (char *)&hc_ref(i__, j), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	    io___86.ciunit = iounit_1.iout;
	    s_wsfe(&io___86);
	    e_wsfe();
	    i__1 = piorbs_1.norbit;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___87.ciunit = iounit_1.iout;
		s_wsfe(&io___87);
		i__2 = piorbs_1.norbit;
		for (j = 1; j <= i__2; ++j) {
		    do_fio(&c__1, (char *)&gamma_ref(i__, j), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	}

/*     now, get the bond orders (compute p and p*b) */

	if (inform_1.debug) {
	    io___88.ciunit = iounit_1.iout;
	    s_wsfe(&io___88);
	    e_wsfe();
	}
	i__1 = piorbs_1.nbpi;
	for (k = 1; k <= i__1; ++k) {
	    i__ = ibpi_ref(2, k);
	    j = ibpi_ref(3, k);
	    p = 0.;
	    i__2 = orbits_1.nfill;
	    for (m = 1; m <= i__2; ++m) {
		p += v_ref(i__, m) * 2. * v_ref(j, m);
	    }
	    if (s_cmp(mode, "PLANAR", (ftnlen)6, (ftnlen)6) == 0) {
		border_1.pbpl[k - 1] = p * hc_ref(i__, j) / ebeta;
	    } else if (s_cmp(mode, "NONPLN", (ftnlen)6, (ftnlen)6) == 0) {
		border_1.pnpl[k - 1] = p;
	    }
	    if (inform_1.debug) {
		i__ = ibnd_ref(1, ibpi_ref(1, k));
		j = ibnd_ref(2, ibpi_ref(1, k));
		io___90.ciunit = iounit_1.iout;
		s_wsfe(&io___90);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&p, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}

/*     if we have done planar calculation, do the nonplanar; */
/*     when both are complete, alter the pisystem constants */

	if (s_cmp(mode, "PLANAR", (ftnlen)6, (ftnlen)6) == 0) {
	    s_copy(mode, "NONPLN", (ftnlen)6, (ftnlen)6);
	} else if (s_cmp(mode, "NONPLN", (ftnlen)6, (ftnlen)6) == 0) {
	    s_copy(mode, "      ", (ftnlen)6, (ftnlen)6);
	}
    }

/*     alter torsional and bond constants for pisystem */

    pialter_();
    return 0;
} /* piscf_ */

#undef gamma_ref
#undef ibpi_ref
#undef fock_ref
#undef ibnd_ref
#undef ed_ref
#undef hc_ref
#undef v_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine pitilt  --  direction cosines for pisystem  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "pitilt" calculates for each pibond the ratio of the */
/*     actual p-orbital overlap integral to the ideal overlap */
/*     if the same orbitals were perfectly parallel */


/* Subroutine */ int pitilt_(doublereal *povlap)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal a1, b1, c1, a2, b2, c2, x2, y2, z2, x3, y3, z3, xr[8], 
	    yr[8], zr[8], rij, xij, yij, zij;
    static integer iorb, jorb, list[8];
    static doublereal ideal, rnorm, cosine;
    extern /* Subroutine */ int pimove_(integer *, doublereal *, doublereal *,
	     doublereal *), overlap_(integer *, integer *, doublereal *, 
	    doublereal *);


#define ibpi_ref(a_1,a_2) piorbs_1.ibpi[(a_2)*3 + a_1 - 4]
#define piperp_ref(a_1,a_2) piorbs_1.piperp[(a_2)*3 + a_1 - 4]



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  piorbs.i  --  conjugated system in the current structure  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     norbit    total number of pisystem orbitals in the system */
/*     iorbit    numbers of the atoms containing pisystem orbitals */
/*     reorbit   number of evaluations between orbital updates */
/*     piperp    atoms defining a normal plane to each orbital */
/*     nbpi      total number of bonds affected by the pisystem */
/*     ibpi      bond and piatom numbers for each pisystem bond */
/*     ntpi      total number of torsions affected by the pisystem */
/*     itpi      torsion and pibond numbers for each pisystem torsion */
/*     listpi    atom list indicating whether each atom has an orbital */




/*     planes defining each p-orbital are in "piperp"; transform */
/*     coordinates of "iorb", "jorb" and their associated planes */
/*     to put "iorb" at origin and "jorb" along the x-axis */

    /* Parameter adjustments */
    --povlap;

    /* Function Body */
    i__1 = piorbs_1.nbpi;
    for (k = 1; k <= i__1; ++k) {
	i__ = ibpi_ref(2, k);
	j = ibpi_ref(3, k);
	iorb = piorbs_1.iorbit[i__ - 1];
	jorb = piorbs_1.iorbit[j - 1];
	list[0] = iorb;
	list[1] = jorb;
	for (m = 1; m <= 3; ++m) {
	    list[m + 1] = piperp_ref(m, i__);
	    list[m + 4] = piperp_ref(m, j);
	}
	pimove_(list, xr, yr, zr);

/*     check for sp-hybridized carbon in current bond; */
/*     assume perfect overlap for any such pibond */

	if (atmtyp_1.atomic[iorb - 1] == 6 && couple_1.n12[iorb - 1] == 2 || 
		atmtyp_1.atomic[jorb - 1] == 6 && couple_1.n12[jorb - 1] == 2)
		 {
	    povlap[k] = 1.;

/*     find and normalize a vector parallel to first p-orbital */

	} else {
	    x2 = xr[3] - xr[2];
	    y2 = yr[3] - yr[2];
	    z2 = zr[3] - zr[2];
	    x3 = xr[4] - xr[2];
	    y3 = yr[4] - yr[2];
	    z3 = zr[4] - zr[2];
	    a1 = y2 * z3 - y3 * z2;
	    b1 = x3 * z2 - x2 * z3;
	    c1 = x2 * y3 - x3 * y2;
	    rnorm = sqrt(a1 * a1 + b1 * b1 + c1 * c1);
	    a1 /= rnorm;
	    b1 /= rnorm;
	    c1 /= rnorm;

/*     now find vector parallel to the second p-orbital, */
/*     "a2" changes sign to correspond to internuclear axis */

	    x2 = xr[6] - xr[5];
	    y2 = yr[6] - yr[5];
	    z2 = zr[6] - zr[5];
	    x3 = xr[7] - xr[5];
	    y3 = yr[7] - yr[5];
	    z3 = zr[7] - zr[5];
	    a2 = y2 * z3 - y3 * z2;
	    b2 = x3 * z2 - x2 * z3;
	    c2 = x2 * y3 - x3 * y2;
	    rnorm = sqrt(a2 * a2 + b2 * b2 + c2 * c2);
	    a2 = -a2 / rnorm;
	    b2 /= rnorm;
	    c2 /= rnorm;

/*     compute the cosine of the angle between p-orbitals; */
/*     if more than 90 degrees, reverse one of the vectors */

	    cosine = a1 * a2 + b1 * b2 + c1 * c2;
	    if (cosine < 0.) {
		a2 = -a2;
		b2 = -b2;
		c2 = -c2;
	    }

/*     find overlap if the orbitals were perfectly parallel */

	    xij = atoms_1.x[iorb - 1] - atoms_1.x[jorb - 1];
	    yij = atoms_1.y[iorb - 1] - atoms_1.y[jorb - 1];
	    zij = atoms_1.z__[iorb - 1] - atoms_1.z__[jorb - 1];
/* Computing 2nd power */
	    d__1 = xij;
/* Computing 2nd power */
	    d__2 = yij;
/* Computing 2nd power */
	    d__3 = zij;
	    rij = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	    overlap_(&atmtyp_1.atomic[iorb - 1], &atmtyp_1.atomic[jorb - 1], &
		    rij, &ideal);

/*     set ratio of actual to ideal overlap for current pibond */

	    povlap[k] = ideal * a1 * a2 + b1 * b2 + c1 * c2;
	}
    }
    return 0;
} /* pitilt_ */

#undef piperp_ref
#undef ibpi_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine pimove  --  translate & rotate bond vector  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "pimove" rotates the vector between atoms "list(1)" and */
/*     "list(2)" so that atom 1 is at the origin and atom 2 along */
/*     the x-axis; the atoms defining the respective planes are */
/*     also moved and their bond lengths normalized */


/* Subroutine */ int pimove_(integer *list, doublereal *xr, doublereal *yr, 
	doublereal *zr)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal xt, yt, zt, sine, xold, denom, cosine;



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




/*     translate "list" atoms to place atom 1 at origin */

    /* Parameter adjustments */
    --zr;
    --yr;
    --xr;
    --list;

    /* Function Body */
    j = list[1];
    xt = atoms_1.x[j - 1];
    yt = atoms_1.y[j - 1];
    zt = atoms_1.z__[j - 1];
    for (i__ = 1; i__ <= 8; ++i__) {
	j = list[i__];
	xr[i__] = atoms_1.x[j - 1] - xt;
	yr[i__] = atoms_1.y[j - 1] - yt;
	zr[i__] = atoms_1.z__[j - 1] - zt;
    }

/*     rotate "list" atoms to place atom 2 on the x-axis */

/* Computing 2nd power */
    d__1 = xr[2];
/* Computing 2nd power */
    d__2 = yr[2];
    denom = sqrt(d__1 * d__1 + d__2 * d__2);
    if (denom != 0.) {
	sine = yr[2] / denom;
	cosine = xr[2] / denom;
	for (i__ = 1; i__ <= 8; ++i__) {
	    xold = xr[i__];
	    xr[i__] = xr[i__] * cosine + yr[i__] * sine;
	    yr[i__] = yr[i__] * cosine - xold * sine;
	}
    }
/* Computing 2nd power */
    d__1 = xr[2];
/* Computing 2nd power */
    d__2 = zr[2];
    denom = sqrt(d__1 * d__1 + d__2 * d__2);
    if (denom != 0.) {
	sine = zr[2] / denom;
	cosine = xr[2] / denom;
	for (i__ = 1; i__ <= 8; ++i__) {
	    xold = xr[i__];
	    xr[i__] = xr[i__] * cosine + zr[i__] * sine;
	    zr[i__] = zr[i__] * cosine - xold * sine;
	}
    }

/*     normalize the coordinates of atoms defining the plane */
/*     for atom 1 (ie, make all these atoms have unit length to */
/*     atom 1) so that the orbital makes equal angles with the */
/*     atoms rather than simply being perpendicular to the common */
/*     plane of the atoms */

    for (i__ = 3; i__ <= 5; ++i__) {
	if (list[i__] != list[1]) {
/* Computing 2nd power */
	    d__1 = xr[i__];
/* Computing 2nd power */
	    d__2 = yr[i__];
/* Computing 2nd power */
	    d__3 = zr[i__];
	    denom = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	    xr[i__] /= denom;
	    yr[i__] /= denom;
	    zr[i__] /= denom;
	}
    }

/*     normalization of plane defining atoms for atom 2; for the */
/*     x-coordinate we translate back to the origin, normalize */
/*     and then retranslate back along the x-axis */

    for (i__ = 6; i__ <= 8; ++i__) {
	if (list[i__] != list[2]) {
/* Computing 2nd power */
	    d__1 = xr[i__] - xr[2];
/* Computing 2nd power */
	    d__2 = yr[i__];
/* Computing 2nd power */
	    d__3 = zr[i__];
	    denom = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	    xr[i__] = (xr[i__] - xr[2]) / denom + xr[2];
	    yr[i__] /= denom;
	    zr[i__] /= denom;
	}
    }
    return 0;
} /* pimove_ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine pialter  --  modify constants for pisystem  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "pialter" first modifies bond lengths and force constants */
/*     according to the standard bond slope parameters and the */
/*     bond order values stored in "pnpl"; also alters some 2-fold */
/*     torsional parameters based on the bond-order * beta matrix */


/* Subroutine */ int pialter_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Altered Bond Stretching Parameters\002"
	    ",\002 for Pi-System :\002,//,\002 Type\002,13x,\002Atom Names"
	    "\002,15x,\002Initial\002,16x,\002Final\002,/)";
    static char fmt_20[] = "(\002 Bond\002,7x,i5,\002-\002,a3,1x,i5,\002-"
	    "\002,a3,5x,f9.3,f8.4,2x,\002-->\002,f9.3,f8.4)";
    static char fmt_30[] = "(/,\002 Altered 2-Fold Torsional Parameters\002"
	    ",\002 for Pi-System :\002,//,\002 Type\002,23x,\002Atom Names"
	    "\002,17x,\002Initial\002,8x,\002Final\002,/)";
    static char fmt_40[] = "(\002 Torsion\002,4x,i5,\002-\002,a3,3(1x,i5,"
	    "\002-\002,a3),3x,f8.3,2x,\002-->\002,f8.3)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, ia, ib, ic, id;

    /* Fortran I/O blocks */
    static cilist io___129 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___134 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___135 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___139 = { 0, 0, 0, fmt_40, 0 };



#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define ibpi_ref(a_1,a_2) piorbs_1.ibpi[(a_2)*3 + a_1 - 4]
#define itpi_ref(a_1,a_2) piorbs_1.itpi[(a_2)*2 + a_1 - 3]
#define tors2_ref(a_1,a_2) tors_1.tors2[(a_2)*4 + a_1 - 5]
#define itors_ref(a_1,a_2) tors_1.itors[(a_2)*4 + a_1 - 5]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  bond.i  --  covalent bonds in the current structure  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     bk      bond stretch force constants (kcal/mole/Ang**2) */
/*     bl      ideal bond length values in Angstroms */
/*     nbond   total number of bond stretches in the system */
/*     ibnd    numbers of the atoms in each bond stretch */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  border.i  --  bond orders for a conjugated pisystem  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     pbpl   pi-bond orders for bonds in "planar" pisystem */
/*     pnpl   pi-bond orders for bonds in "nonplanar" pisystem */




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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  piorbs.i  --  conjugated system in the current structure  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     norbit    total number of pisystem orbitals in the system */
/*     iorbit    numbers of the atoms containing pisystem orbitals */
/*     reorbit   number of evaluations between orbital updates */
/*     piperp    atoms defining a normal plane to each orbital */
/*     nbpi      total number of bonds affected by the pisystem */
/*     ibpi      bond and piatom numbers for each pisystem bond */
/*     ntpi      total number of torsions affected by the pisystem */
/*     itpi      torsion and pibond numbers for each pisystem torsion */
/*     listpi    atom list indicating whether each atom has an orbital */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  pistuf.i  --  bonds and torsions in the current pisystem  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     bkpi     bond stretch force constants for pi-bond order of 1.0 */
/*     blpi     ideal bond length values for a pi-bond order of 1.0 */
/*     kslope   rate of force constant decrease with bond order decrease */
/*     lslope   rate of bond length increase with a bond order decrease */
/*     torsp2   2-fold torsional energy barrier for pi-bond order of 1.0 */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  tors.i  --  torsional angles within the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     tors1   1-fold amplitude and phase for each torsional angle */
/*     tors2   2-fold amplitude and phase for each torsional angle */
/*     tors3   3-fold amplitude and phase for each torsional angle */
/*     tors4   4-fold amplitude and phase for each torsional angle */
/*     tors5   5-fold amplitude and phase for each torsional angle */
/*     tors6   6-fold amplitude and phase for each torsional angle */
/*     ntors   total number of torsional angles in the system */
/*     itors   numbers of the atoms in each torsional angle */




/*     modify the stretching constants and natural bond lengths */

    if (inform_1.debug && piorbs_1.nbpi != 0) {
	io___129.ciunit = iounit_1.iout;
	s_wsfe(&io___129);
	e_wsfe();
    }
    i__1 = piorbs_1.nbpi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = ibpi_ref(1, i__);
	ia = ibnd_ref(1, k);
	ib = ibnd_ref(2, k);
	bond_1.bk[k - 1] = pistuf_1.bkpi[i__ - 1] - pistuf_1.kslope[i__ - 1] *
		 (1. - border_1.pnpl[i__ - 1]);
	bond_1.bl[k - 1] = pistuf_1.blpi[i__ - 1] + pistuf_1.lslope[i__ - 1] *
		 (1. - border_1.pnpl[i__ - 1]);
	if (inform_1.debug) {
	    io___134.ciunit = iounit_1.iout;
	    s_wsfe(&io___134);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ia), (ftnlen)3);
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ib), (ftnlen)3);
	    do_fio(&c__1, (char *)&pistuf_1.bkpi[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pistuf_1.blpi[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&bond_1.bk[k - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&bond_1.bl[k - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }

/*     modify the 2-fold torsional constants across pibonds */

    if (inform_1.debug && piorbs_1.ntpi != 0) {
	io___135.ciunit = iounit_1.iout;
	s_wsfe(&io___135);
	e_wsfe();
    }
    i__1 = piorbs_1.ntpi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = itpi_ref(1, i__);
	k = itpi_ref(2, i__);
	ia = itors_ref(1, j);
	ib = itors_ref(2, j);
	ic = itors_ref(3, j);
	id = itors_ref(4, j);
	tors2_ref(1, j) = border_1.pbpl[k - 1] * pistuf_1.torsp2[i__ - 1];
	if (inform_1.debug) {
	    io___139.ciunit = iounit_1.iout;
	    s_wsfe(&io___139);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ia), (ftnlen)3);
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ib), (ftnlen)3);
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ic), (ftnlen)3);
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, id), (ftnlen)3);
	    do_fio(&c__1, (char *)&pistuf_1.torsp2[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&tors2_ref(1, j), (ftnlen)sizeof(doublereal)
		    );
	    e_wsfe();
	}
    }
    return 0;
} /* pialter_ */

#undef itors_ref
#undef tors2_ref
#undef itpi_ref
#undef ibpi_ref
#undef name___ref
#undef ibnd_ref


