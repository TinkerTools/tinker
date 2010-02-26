/* korbit.f -- translated by f2c (version 20050501).
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
    doublereal electron[1000], ionize[1000], repulse[1000], sslope[500], 
	    tslope[500], sslope5[200], tslope5[200], sslope4[200], tslope4[
	    200];
    char kpi[4000], kpi5[1600], kpi4[1600];
} korbs_;

#define korbs_1 korbs_

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

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__2 = 2;
static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine korbit  --  pisystem parameter assignment  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "korbit" assigns pi-orbital parameters to conjugated */
/*     systems and processes any new or changed parameters */


/* Subroutine */ int korbit_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Additional Pisystem Atom Parameters :"
	    "\002,//,6x,\002Atom Type\002,3x,\002Electron\002,3x,\002Ionizati"
	    "on\002,3x,\002Repulsion\002,/)";
    static char fmt_30[] = "(8x,i4,3x,f10.3,3x,f10.3,2x,f10.3)";
    static char fmt_60[] = "(/,\002 Additional Pisystem Bond Parameters :"
	    "\002,//,6x,\002Atom Types\002,7x,\002d Force\002,4x,\002d Lengt"
	    "h\002,/)";
    static char fmt_70[] = "(6x,2i4,5x,2f11.3)";
    static char fmt_80[] = "(6x,2i4,5x,2f11.3,3x,a6)";
    static char fmt_90[] = "(/,\002 KORBIT  --  Too many Pisystem Bond\002"
	    ",\002 Type Parameters\002)";
    static char fmt_110[] = "(/,\002 KORBIT  --  Too many 5-Ring Pisystem Bo"
	    "nd\002,\002 Type Parameters\002)";
    static char fmt_130[] = "(/,\002 KORBIT  --  Too many 4-Ring Pisystem Bo"
	    "nd\002,\002 Type Parameters\002)";
    static char fmt_150[] = "(/,\002 Undefined Conjugated Pibond Parameter"
	    "s :\002,//,\002 Type\002,13x,\002Atom Names\002,11x,\002Atom Cla"
	    "sses\002,/)";
    static char fmt_160[] = "(1x,a6,5x,i6,\002-\002,a3,i6,\002-\002,a3,7x,2i"
	    "5)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer i_dnnt(doublereal *);

    /* Local variables */
    static logical use_ring__;
    static integer i__, j, k, ia, ib;
    static char pa[4], pb[4];
    static integer it;
    static char pt[8];
    static integer ita, itb, npi, npi4, npi5, size, next;
    static char label[6], blank[8];
    static doublereal elect;
    static integer iring;
    static doublereal ioniz, sslop, tslop;
    static logical header;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static doublereal repuls;
    static char string[120];
    extern /* Subroutine */ int chkring_(integer *, integer *, integer *, 
	    integer *, integer *), numeral_(integer *, char *, integer *, 
	    ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___12 = { 1, string, 0, 0, 120, 1 };
    static cilist io___13 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_30, 0 };
    static icilist io___19 = { 1, string, 0, 0, 120, 1 };
    static cilist io___20 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_160, 0 };



#define kpi_ref(a_0,a_1) &korbs_1.kpi[(a_1)*8 + a_0 - 8]
#define kpi4_ref(a_0,a_1) &korbs_1.kpi4[(a_1)*8 + a_0 - 8]
#define kpi5_ref(a_0,a_1) &korbs_1.kpi5[(a_1)*8 + a_0 - 8]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define ibpi_ref(a_1,a_2) piorbs_1.ibpi[(a_2)*3 + a_1 - 4]
#define itpi_ref(a_1,a_2) piorbs_1.itpi[(a_2)*2 + a_1 - 3]
#define tors2_ref(a_1,a_2) tors_1.tors2[(a_2)*4 + a_1 - 5]
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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  korbs.i  --  forcefield parameters for pisystem orbitals  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxnpi     maximum number of pisystem bond parameter entries */
/*     maxnpi5    maximum number of 5-membered ring pibond entries */
/*     maxnpi4    maximum number of 4-membered ring pibond entries */

/*     electron   number of pi-electrons for each atom class */
/*     ionize     ionization potential for each atom class */
/*     repulse    repulsion integral value for each atom class */
/*     sslope     slope for bond stretch vs. pi-bond order */
/*     tslope     slope for 2-fold torsion vs. pi-bond order */
/*     sslope5    slope for 5-ring bond stretch vs. pi-bond order */
/*     tslope5    slope for 5-ring 2-fold torsion vs. pi-bond order */
/*     sslope4    slope for 4-ring bond stretch vs. pi-bond order */
/*     tslope4    slope for 4-ring 2-fold torsion vs. pi-bond order */
/*     kpi        string of atom classes for pisystem bonds */
/*     kpi5       string of atom classes for 5-ring pisystem bonds */
/*     kpi4       string of atom classes for 4-ring pisystem bonds */




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




/*     process keywords containing pisystem atom parameters */

    s_copy(blank, "        ", (ftnlen)8, (ftnlen)8);
    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "PIATOM ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    elect = 0.;
	    ioniz = 0.;
	    repuls = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___12);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&elect, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&ioniz, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&repuls, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
L10:
	    if (header) {
		header = FALSE_;
		io___13.ciunit = iounit_1.iout;
		s_wsfe(&io___13);
		e_wsfe();
	    }
	    io___14.ciunit = iounit_1.iout;
	    s_wsfe(&io___14);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&elect, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ioniz, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&repuls, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    if (ia > 0 && ia <= 1000) {
		korbs_1.electron[ia - 1] = elect;
		korbs_1.ionize[ia - 1] = ioniz;
		korbs_1.repulse[ia - 1] = repuls;
	    } else {
/* L40: */
		inform_1.abort = TRUE_;
	    }
	}
    }

/*     process keywords containing pisystem bond parameters */

    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	iring = -1;
	if (s_cmp(keyword, "PIBOND ", (ftnlen)7, (ftnlen)7) == 0) {
	    iring = 0;
	}
	if (s_cmp(keyword, "PIBOND5 ", (ftnlen)8, (ftnlen)8) == 0) {
	    iring = 5;
	}
	if (s_cmp(keyword, "PIBOND4 ", (ftnlen)8, (ftnlen)8) == 0) {
	    iring = 4;
	}
	if (iring >= 0) {
	    ia = 0;
	    ib = 0;
	    sslop = 0.;
	    tslop = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___19);
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&sslop, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&tslop, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L50;
	    }
L50:
	    if (header) {
		header = FALSE_;
		io___20.ciunit = iounit_1.iout;
		s_wsfe(&io___20);
		e_wsfe();
	    }
	    if (iring == 0) {
		io___21.ciunit = iounit_1.iout;
		s_wsfe(&io___21);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&sslop, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&tslop, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		if (iring == 5) {
		    s_copy(label, "5-Ring", (ftnlen)6, (ftnlen)6);
		}
		if (iring == 4) {
		    s_copy(label, "4-Ring", (ftnlen)6, (ftnlen)6);
		}
		io___23.ciunit = iounit_1.iout;
		s_wsfe(&io___23);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&sslop, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&tslop, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, label, (ftnlen)6);
		e_wsfe();
	    }
	    size = 4;
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    if (ia <= ib) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pb;
		i__3[1] = 4, a__1[1] = pa;
		s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
	    }
	    if (iring == 0) {
		for (j = 1; j <= 500; ++j) {
		    if (s_cmp(kpi_ref(0, j), blank, (ftnlen)8, (ftnlen)8) == 
			    0 || s_cmp(kpi_ref(0, j), pt, (ftnlen)8, (ftnlen)
			    8) == 0) {
			s_copy(kpi_ref(0, j), pt, (ftnlen)8, (ftnlen)8);
			korbs_1.sslope[j - 1] = sslop;
			korbs_1.tslope[j - 1] = tslop;
			goto L100;
		    }
		}
		io___29.ciunit = iounit_1.iout;
		s_wsfe(&io___29);
		e_wsfe();
		inform_1.abort = TRUE_;
L100:
		;
	    } else if (iring == 5) {
		for (j = 1; j <= 200; ++j) {
		    if (s_cmp(kpi5_ref(0, j), blank, (ftnlen)8, (ftnlen)8) == 
			    0 || s_cmp(kpi5_ref(0, j), pt, (ftnlen)8, (ftnlen)
			    8) == 0) {
			s_copy(kpi5_ref(0, j), pt, (ftnlen)8, (ftnlen)8);
			korbs_1.sslope5[j - 1] = sslop;
			korbs_1.tslope5[j - 1] = tslop;
			goto L120;
		    }
		}
		io___30.ciunit = iounit_1.iout;
		s_wsfe(&io___30);
		e_wsfe();
		inform_1.abort = TRUE_;
L120:
		;
	    } else if (iring == 4) {
		for (j = 1; j <= 200; ++j) {
		    if (s_cmp(kpi4_ref(0, j), blank, (ftnlen)8, (ftnlen)8) == 
			    0 || s_cmp(kpi4_ref(0, j), pt, (ftnlen)8, (ftnlen)
			    8) == 0) {
			s_copy(kpi4_ref(0, j), pt, (ftnlen)8, (ftnlen)8);
			korbs_1.sslope4[j - 1] = sslop;
			korbs_1.tslope4[j - 1] = tslop;
			goto L140;
		    }
		}
		io___31.ciunit = iounit_1.iout;
		s_wsfe(&io___31);
		e_wsfe();
		inform_1.abort = TRUE_;
L140:
		;
	    }
	}
    }

/*     determine the total number of forcefield parameters */

    npi = 500;
    npi5 = 200;
    npi4 = 200;
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(kpi_ref(0, i__), blank, (ftnlen)8, (ftnlen)8) == 0) {
	    npi = i__ - 1;
	}
    }
    for (i__ = 200; i__ >= 1; --i__) {
	if (s_cmp(kpi5_ref(0, i__), blank, (ftnlen)8, (ftnlen)8) == 0) {
	    npi5 = i__ - 1;
	}
    }
    for (i__ = 200; i__ >= 1; --i__) {
	if (s_cmp(kpi4_ref(0, i__), blank, (ftnlen)8, (ftnlen)8) == 0) {
	    npi4 = i__ - 1;
	}
    }
    use_ring__ = FALSE_;
    if (min(npi5,npi4) != 0) {
	use_ring__ = TRUE_;
    }

/*     assign the values characteristic of the piatom types; */
/*     count the number of filled pi molecular orbitals */

    orbits_1.nfill = 0;
    i__1 = piorbs_1.norbit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	it = atoms_1.type__[piorbs_1.iorbit[i__ - 1] - 1];
	orbits_1.q[i__ - 1] = korbs_1.electron[it - 1];
	orbits_1.w[i__ - 1] = korbs_1.ionize[it - 1] / 27.21138386;
	orbits_1.em[i__ - 1] = korbs_1.repulse[it - 1] / 27.21138386;
	orbits_1.nfill += i_dnnt(&orbits_1.q[i__ - 1]);
    }
    orbits_1.nfill /= 2;

/*     assign parameters for all bonds between piatoms; */
/*     store the original bond lengths and force constants */

    i__1 = piorbs_1.nbpi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = ibpi_ref(1, i__);
	ia = ibpi_ref(2, i__);
	ib = ibpi_ref(3, i__);
	ita = atmtyp_1.class__[piorbs_1.iorbit[ia - 1] - 1];
	itb = atmtyp_1.class__[piorbs_1.iorbit[ib - 1] - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	if (ita <= itb) {
/* Writing concatenation */
	    i__3[0] = 4, a__1[0] = pa;
	    i__3[1] = 4, a__1[1] = pb;
	    s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
	} else {
/* Writing concatenation */
	    i__3[0] = 4, a__1[0] = pb;
	    i__3[1] = 4, a__1[1] = pa;
	    s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
	}

/*     make a check for bonds contained inside small rings */

	iring = 0;
	if (use_ring__) {
	    chkring_(&iring, &ia, &ib, &c__0, &c__0);
	    if (iring == 6) {
		iring = 0;
	    }
	    if (iring == 5 && npi5 == 0) {
		iring = 0;
	    }
	    if (iring == 4 && npi4 == 0) {
		iring = 0;
	    }
	    if (iring == 3) {
		iring = 0;
	    }
	}

/*     assign conjugated bond parameters for each pibond */

	if (iring == 0) {
	    i__2 = npi;
	    for (k = 1; k <= i__2; ++k) {
		if (s_cmp(kpi_ref(0, k), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    pistuf_1.bkpi[i__ - 1] = bond_1.bk[j - 1];
		    pistuf_1.blpi[i__ - 1] = bond_1.bl[j - 1];
		    pistuf_1.kslope[i__ - 1] = korbs_1.sslope[k - 1];
		    pistuf_1.lslope[i__ - 1] = korbs_1.tslope[k - 1];
		    goto L170;
		}
	    }

/*     assign bond parameters for 5-membered ring pibonds */

	} else if (iring == 5) {
	    i__2 = npi5;
	    for (k = 1; k <= i__2; ++k) {
		if (s_cmp(kpi5_ref(0, k), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    pistuf_1.bkpi[i__ - 1] = bond_1.bk[j - 1];
		    pistuf_1.blpi[i__ - 1] = bond_1.bl[j - 1];
		    pistuf_1.kslope[i__ - 1] = korbs_1.sslope5[k - 1];
		    pistuf_1.lslope[i__ - 1] = korbs_1.tslope5[k - 1];
		    goto L170;
		}
	    }

/*     assign bond parameters for 4-membered ring pibonds */

	} else if (iring == 4) {
	    i__2 = npi4;
	    for (k = 1; k <= i__2; ++k) {
		if (s_cmp(kpi4_ref(0, k), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    pistuf_1.bkpi[i__ - 1] = bond_1.bk[j - 1];
		    pistuf_1.blpi[i__ - 1] = bond_1.bl[j - 1];
		    pistuf_1.kslope[i__ - 1] = korbs_1.sslope4[k - 1];
		    pistuf_1.lslope[i__ - 1] = korbs_1.tslope4[k - 1];
		    goto L170;
		}
	    }
	}

/*     warning if suitable conjugated pibond parameters not found */

	inform_1.abort = TRUE_;
	if (header) {
	    header = FALSE_;
	    io___40.ciunit = iounit_1.iout;
	    s_wsfe(&io___40);
	    e_wsfe();
	}
	s_copy(label, "Pibond", (ftnlen)6, (ftnlen)6);
	if (iring == 5) {
	    s_copy(label, "5-Ring", (ftnlen)6, (ftnlen)6);
	}
	if (iring == 4) {
	    s_copy(label, "4-Ring", (ftnlen)6, (ftnlen)6);
	}
	io___41.ciunit = iounit_1.iout;
	s_wsfe(&io___41);
	do_fio(&c__1, label, (ftnlen)6);
	do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	do_fio(&c__1, name___ref(0, ia), (ftnlen)3);
	do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	do_fio(&c__1, name___ref(0, ib), (ftnlen)3);
	do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
	e_wsfe();
L170:
	;
    }

/*     store the original torsional constants across pibonds */

    i__1 = piorbs_1.ntpi;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = itpi_ref(1, i__);
	pistuf_1.torsp2[i__ - 1] = tors2_ref(1, j);
    }
    return 0;
} /* korbit_ */

#undef keyline_ref
#undef tors2_ref
#undef itpi_ref
#undef ibpi_ref
#undef name___ref
#undef kpi5_ref
#undef kpi4_ref
#undef kpi_ref


