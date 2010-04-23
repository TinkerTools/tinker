/* vibbig.f -- translated by f2c (version 20050501).
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
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

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
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal vdwcut, chgcut, dplcut, mpolecut, vdwtaper, chgtaper, dpltaper,
	     mpoletaper, ewaldcut;
    logical use_ewald__, use_lights__, use_list__, use_vlist__, use_clist__, 
	    use_mlist__;
} cutoff_;

#define cutoff_1 cutoff_

struct {
    doublereal hesscut;
} hescut_;

#define hescut_1 hescut_

struct {
    doublereal hessx[75000]	/* was [3][25000] */, hessy[75000]	/* 
	    was [3][25000] */, hessz[75000]	/* was [3][25000] */;
} hessn_;

#define hessn_1 hessn_

struct {
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

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
    doublereal xrb[25000], yrb[25000], zrb[25000], rbc[6000]	/* was [6][
	    1000] */;
    logical use_rigid__;
} rigid_;

#define rigid_1 rigid_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal radmin[1000000]	/* was [1000][1000] */, epsilon[1000000]	
	    /* was [1000][1000] */, radmin4[1000000]	/* was [1000][1000] */
	    , epsilon4[1000000]	/* was [1000][1000] */, radhbnd[1000000]	
	    /* was [1000][1000] */, epshbnd[1000000]	/* was [1000][1000] */
	    , kred[25000];
    integer ired[25000], nvdw, ivdw[25000], jvdw[25000], nvt, ivt[25000], jvt[
	    25000];
} vdw_;

#define vdw_1 vdw_

struct {
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
    logical use_vcorr__;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__2 = 2;
static doublereal c_b75 = 1.;
static integer c__0 = 0;
static integer c__150 = 150;
static integer c__1000 = 1000;



/*     ################################################################# */
/*     ##  COPYRIGHT (C) 2007 by Alexey Kaledin & Jay William Ponder  ## */
/*     ##                     All Rights Reserved                     ## */
/*     ################################################################# */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  program vibbig  --  block iterative vibrational analysis  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "vibbig" performs large-scale vibrational mode analysis using */
/*     only vector storage and gradient evaluations; preconditioning */
/*     is via an approximate inverse from a block diagonal Hessian, */
/*     and a sliding block method is used to converge any number of */
/*     eigenvectors starting from either lowest or highest frequency */

/*     literature references: */

/*     C. Murray, S. C. Racine and E. R. Davidson, "Improved Algorithms */
/*     for the Lowest Few Eigenvalues and Associated Eigenvectors of */
/*     Large Matrices", Journal of Computational Physics, 103, 382-389 */
/*     (1992) */

/*     A. L. Kaledin, "Gradient-Based Direct Normal-Mode Analysis", */
/*     Journal of Chemical Physics, 122, 184106 (2005) */

/*     A. L. Kaledin, M. Kaledin and J. M. Bowman, "All-Atom Calculation */
/*     of the Normal Modes of Bacteriorhodopsin Using a Sliding Block */
/*     Iterative Diagonalization Method", Journal of Chemical Theory */
/*     and Computation, 2, 166-174 (2006) */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Start at Lowest or Highest Frequenc"
	    "y\002,\002 Normal Mode [\002,a1,\002] :  \002,$)";
    static char fmt_30[] = "(a120)";
    static char fmt_50[] = "(/,\002 Enter Desired Frequency Cutoff in cm-"
	    "1\002,\002 [0.0] :  \002,$)";
    static char fmt_60[] = "(f20.0)";
    static char fmt_90[] = "(/,\002 Atom Blocks Used to Subdivide the System"
	    " :\002,/)";
    static char fmt_100[] = "(\002 Block :\002,i7,9x,\002Size :\002,i7,9x"
	    ",\002Atoms :\002,i7,\002  to\002,i7)";
    static char fmt_110[] = "(/,\002 Storage for Preconditioning Array :\002"
	    ",5x,i12)";
    static char fmt_120[] = "(/,\002 VIBBIG  --  Preconditioning Too Large"
	    ";\002,\002 Increase MAXHESS\002)";
    static char fmt_140[] = "(/,\002 Prior Normal Modes Available at Restart"
	    " :\002,i11)";
    static char fmt_170[] = "(/,\002 Number of Converged Normal Modes :\002,"
	    "6x,i12)";
    static char fmt_180[] = "(/,\002 VIBBIG  --  Loss of Root Identity; Plea"
	    "se\002,\002 Try to Restart\002)";
    static char fmt_200[] = "(/,\002 Converged Normal Modes from Iterativ"
	    "e\002,\002 Vibrational Analysis :\002)";
    static char fmt_210[] = "(/,4x,\002Mode\002,7x,\002Frequency\002,8x,\002"
	    "Delta\002,10x,\002R Norm\002,10x,\002Orthog\002)";
    static char fmt_220[] = "()";
    static char fmt_230[] = "(i8,f15.3,3d16.4)";
    static char fmt_240[] = "(/,\002 Iteration\002,i7,11x,\002New Modes\002,"
	    "i6,10x,\002 Total Modes\002,i6,/)";
    static char fmt_250[] = "(4x,\002Mode\002,7x,\002Frequency\002,8x,\002De"
	    "lta\002,10x,\002R Norm\002,10x,\002Orthog\002)";
    static char fmt_260[] = "(i8,f15.3,3d16.4)";
    static char fmt_270[] = "(/,\002 Number of Converged Normal Modes :\002,"
	    "6x,i12)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2], i__4;
    doublereal d__1;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void), 
	    s_rsfe(cilist *), e_rsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_inqu(inlist *), f_open(olist *), i_dnnt(doublereal *), s_rsle(
	    cilist *), e_rsle(void), f_clos(cllist *), s_rsue(cilist *), 
	    do_uio(integer *, char *, ftnlen), e_rsue(void), f_rew(alist *);
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);
    integer s_wsue(cilist *), e_wsue(void);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    static char datafile[120];
    extern integer freeunit_(void);
    static doublereal c__[22500]	/* was [150][150] */, h__[22500]	
	    /* was [150][150] */;
    static integer i__, j, k;
    static doublereal p[75000], u[450000]	/* was [75000][6] */;
    static char blockfile[120];
    static integer i1, i2, k0, k1, k2;
    extern /* Subroutine */ int preconblk_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), transform_(integer *, doublereal *, doublereal *);
    static integer ii;
    static doublereal pk[75000], xe[75000];
    static integer np;
    static doublereal xm[75000], ur[225000]	/* was [75000][3] */, uu[
	    1000000], phi[11250000]	/* was [75000][150] */, sum, uku[
	    75000];
    static integer ivb1, ivb2;
    static doublereal tmp1[150], tmp2[150], uku0[75000];
    static integer iblk[25000], nblk;
    static logical done;
    static integer ivib;
    static doublereal fmax, hmin[75000], freq[150], phik[11250000]	/* 
	    was [75000][150] */;
    static integer iter, nvar;
    static doublereal size;
    static integer next;
    static doublereal wtol;
    extern /* Subroutine */ int fatal_(void);
    static doublereal space, dfreq;
    static integer nlock, npair;
    static doublereal shift;
    static integer iconv, idump;
    static doublereal rcomp;
    static integer nconv;
    static doublereal ratio, funit;
    static logical exist;
    static doublereal rnorm;
    extern /* Subroutine */ int gsort_(integer *, integer *, doublereal *, 
	    doublereal *);
    static integer nroot;
    static logical header;
    static integer iblock, irange;
    static doublereal factor;
    static integer nbasis;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char answer[1], string[120];
    extern /* Subroutine */ int konvec_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *), qonvec_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *), prtvib_(integer *, 
	    doublereal *), getxyz_(void), diagblk_(integer *, integer *, 
	    integer *, doublereal *, doublereal *);
    static integer ifactor;
    extern /* Subroutine */ int initial_(void), hessblk_(doublereal *, 
	    integer *, integer *, integer *, doublereal *);
    static doublereal freqold[150];
    extern /* Subroutine */ int trigger_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *);
    static integer maxiter;
    static doublereal uku_min__, uku_max__;
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), trbasis_(
	    integer *, doublereal *, doublereal *, doublereal *), project_(
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *);
    static char keyword[20];
    static logical restart;
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___15 = { 1, string, 1, 0, 120, 1 };
    static icilist io___16 = { 1, string, 1, 0, 120, 1 };
    static icilist io___17 = { 1, string, 1, 0, 120, 1 };
    static icilist io___18 = { 1, string, 1, 0, 120, 1 };
    static cilist io___22 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_30, 0 };
    static icilist io___25 = { 1, string, 1, 0, 120, 1 };
    static cilist io___26 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___50 = { 1, 0, 1, 0, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___52 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___54 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___55 = { 1, 0, 1, 0, 0 };
    static cilist io___57 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___68 = { 0, 0, 0, 0, 0 };
    static cilist io___70 = { 0, 0, 0, 0, 0 };
    static cilist io___83 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___88 = { 0, 0, 0, 0, 0 };
    static cilist io___92 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___93 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___94 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___95 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___96 = { 0, 0, 0, 0, 0 };
    static cilist io___97 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___98 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___99 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___100 = { 0, 0, 0, 0, 0 };
    static cilist io___101 = { 0, 0, 0, 0, 0 };
    static cilist io___102 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___103 = { 0, 0, 0, 0, 0 };
    static cilist io___104 = { 0, 0, 0, 0, 0 };



#define c___ref(a_1,a_2) c__[(a_2)*150 + a_1 - 151]
#define h___ref(a_1,a_2) h__[(a_2)*150 + a_1 - 151]
#define ur_ref(a_1,a_2) ur[(a_2)*75000 + a_1 - 75001]
#define phi_ref(a_1,a_2) phi[(a_2)*75000 + a_1 - 75001]
#define phik_ref(a_1,a_2) phik[(a_2)*75000 + a_1 - 75001]
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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  files.i  --  name and number of current structure files  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     nprior     number of previously existing cycle files */
/*     ldir       length in characters of the directory name */
/*     leng       length in characters of the base filename */
/*     filename   base filename used by default for all files */
/*     outfile    output filename used for intermediate results */




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




/*     set up the structure and mechanics calculation */

    initial_();
    getxyz_();
    mechanic_();

/*     set default parameters for the normal mode computation */

    nvar = atoms_1.n * 3;
    nroot = 50;
    np = 6;
    iter = 0;
    idump = 10;
    maxiter = 100000;
    wtol = 1e-5;
    header = TRUE_;
    restart = TRUE_;

/*     search the keywords for normal mode parameters */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "MAXITER ", (ftnlen)8, (ftnlen)8) == 0) {
	    i__2 = s_rsli(&io___15);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&maxiter, (ftnlen)sizeof(
		    integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "IDUMP ", (ftnlen)6, (ftnlen)6) == 0) {
	    i__2 = s_rsli(&io___16);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&idump, (ftnlen)sizeof(
		    integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "VIB-ROOTS ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    i__2 = s_rsli(&io___17);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&nroot, (ftnlen)sizeof(
		    integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	    nroot = min(nroot,50);
	} else if (s_cmp(keyword, "VIB-TOLERANCE ", (ftnlen)14, (ftnlen)14) ==
		 0) {
	    i__2 = s_rsli(&io___18);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&wtol, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	}
L10:
	;
    }

/*     find either the lowest or highest normal modes */

    factor = 1.;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	*(unsigned char *)answer = 'L';
	io___22.ciunit = iounit_1.iout;
	s_wsfe(&io___22);
	do_fio(&c__1, answer, (ftnlen)1);
	e_wsfe();
	io___23.ciunit = iounit_1.input;
	s_rsfe(&io___23);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'H') {
	factor = -1.;
    }

/*     find cutoff value for desired extreme frequency */

    fmax = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___25);
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&fmax, (ftnlen)sizeof(doublereal))
		;
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L40;
	}
    }
L40:
    if (fmax <= 0.) {
	io___26.ciunit = iounit_1.iout;
	s_wsfe(&io___26);
	e_wsfe();
	io___27.ciunit = iounit_1.input;
	s_rsfe(&io___27);
	do_fio(&c__1, (char *)&fmax, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (fmax <= 0.) {
	fmax = 0.;
    }

/*     set default values for some additional parameters */

    funit = factor * 219474.6313705 * 5.4857990943e-4;
    ifactor = (integer) factor;
/* Computing MAX */
    i__1 = (1 - ifactor) / 2;
    irange = (nvar - np + 1) * max(i__1,0);
    nconv = nlock;
    npair = nroot << 1;
    nbasis = nroot * 3;

/*     open or create eigenvector file for use during restarts */

    ivb1 = freeunit_();
/* Writing concatenation */
    i__3[0] = files_1.leng, a__1[0] = files_1.filename;
    i__3[1] = 4, a__1[1] = ".vb1";
    s_cat(datafile, a__1, i__3, &c__2, (ftnlen)120);
    version_(datafile, "old", (ftnlen)120, (ftnlen)3);
    ioin__1.inerr = 0;
    ioin__1.infilen = 120;
    ioin__1.infile = datafile;
    ioin__1.inex = &exist;
    ioin__1.inopen = 0;
    ioin__1.innum = 0;
    ioin__1.innamed = 0;
    ioin__1.inname = 0;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    if (exist) {
	o__1.oerr = 0;
	o__1.ounit = ivb1;
	o__1.ofnmlen = 120;
	o__1.ofnm = datafile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = "unformatted";
	o__1.oblnk = 0;
	f_open(&o__1);
    } else {
	o__1.oerr = 0;
	o__1.ounit = ivb1;
	o__1.ofnmlen = 120;
	o__1.ofnm = datafile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = "unformatted";
	o__1.oblnk = 0;
	f_open(&o__1);
    }

/*     open or create basis vector file for use during restarts */

    ivb2 = freeunit_();
/* Writing concatenation */
    i__3[0] = files_1.leng, a__1[0] = files_1.filename;
    i__3[1] = 4, a__1[1] = ".vb2";
    s_cat(datafile, a__1, i__3, &c__2, (ftnlen)120);
    version_(datafile, "old", (ftnlen)120, (ftnlen)3);
    ioin__1.inerr = 0;
    ioin__1.infilen = 120;
    ioin__1.infile = datafile;
    ioin__1.inex = &exist;
    ioin__1.inopen = 0;
    ioin__1.innum = 0;
    ioin__1.innamed = 0;
    ioin__1.inname = 0;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    if (exist) {
	o__1.oerr = 0;
	o__1.ounit = ivb2;
	o__1.ofnmlen = 120;
	o__1.ofnm = datafile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = "unformatted";
	o__1.oblnk = 0;
	f_open(&o__1);
    } else {
	restart = FALSE_;
	o__1.oerr = 0;
	o__1.ounit = ivb2;
	o__1.ofnmlen = 120;
	o__1.ofnm = datafile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = "unformatted";
	o__1.oblnk = 0;
	f_open(&o__1);
    }

/*     store a coordinate vector for each atom */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xe[i__ * 3 - 3] = atoms_1.x[i__ - 1] / .52917720859;
	xe[i__ * 3 - 2] = atoms_1.y[i__ - 1] / .52917720859;
	xe[i__ * 3 - 1] = atoms_1.z__[i__ - 1] / .52917720859;
    }

/*     store atomic mass for each coordinate component */

    k = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atmtyp_1.mass[i__ - 1] *= 5.4857990943e-4;
	for (j = 1; j <= 3; ++j) {
	    ++k;
	    xm[k - 1] = atmtyp_1.mass[i__ - 1];
	}
    }

/*     remove pure translational and rotational modes */

    trbasis_(&np, xe, u, ur);

/*     set number and size of blocks based on storage space */

    space = 9e5;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = (doublereal) atoms_1.n;
	size = d__1 * d__1 * 9. / (doublereal) i__;
	if (size < space) {
	    nblk = i__;
	    goto L70;
	}
    }
L70:
    nblk = max(3,nblk);
    size = (doublereal) atoms_1.n / (doublereal) nblk;
    size = min(size,333.33333333333331);
    i__1 = nblk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = (doublereal) i__ * size;
	iblk[i__ - 1] = i_dnnt(&d__1);
    }
    for (i__ = nblk; i__ >= 2; --i__) {
	iblk[i__ - 1] -= iblk[i__ - 2];
    }

/*     get number and size of blocks from an external file */

    iblock = freeunit_();
/* Writing concatenation */
    i__3[0] = files_1.leng, a__1[0] = files_1.filename;
    i__3[1] = 4, a__1[1] = ".blk";
    s_cat(blockfile, a__1, i__3, &c__2, (ftnlen)120);
    version_(blockfile, "old", (ftnlen)120, (ftnlen)3);
    ioin__1.inerr = 0;
    ioin__1.infilen = 120;
    ioin__1.infile = blockfile;
    ioin__1.inex = &exist;
    ioin__1.inopen = 0;
    ioin__1.innum = 0;
    ioin__1.innamed = 0;
    ioin__1.inname = 0;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    if (exist) {
	o__1.oerr = 0;
	o__1.ounit = iblock;
	o__1.ofnmlen = 120;
	o__1.ofnm = blockfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	i__ = 0;
	while(TRUE_) {
	    ++i__;
	    io___50.ciunit = iblock;
	    i__1 = s_rsle(&io___50);
	    if (i__1 != 0) {
		goto L80;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&iblk[i__ - 1], (ftnlen)
		    sizeof(integer));
	    if (i__1 != 0) {
		goto L80;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L80;
	    }
	}
L80:
	nblk = i__ - 1;
	cl__1.cerr = 0;
	cl__1.cunit = iblock;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     print info about the atom blocks and preconditioning */

    io___51.ciunit = iounit_1.iout;
    s_wsfe(&io___51);
    e_wsfe();
    k = 0;
    i__1 = nblk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___52.ciunit = iounit_1.iout;
	s_wsfe(&io___52);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iblk[i__ - 1], (ftnlen)sizeof(integer));
	i__2 = k + 1;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	i__4 = k + iblk[i__ - 1];
	do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
	e_wsfe();
	k += iblk[i__ - 1];
    }
    k = 0;
    i__1 = nblk;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	i__2 = iblk[i__ - 1];
	k += i__2 * i__2 * 9;
    }
    if (k < 1000000) {
	io___53.ciunit = iounit_1.iout;
	s_wsfe(&io___53);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	e_wsfe();
    } else {
	io___54.ciunit = iounit_1.iout;
	s_wsfe(&io___54);
	e_wsfe();
	fatal_();
    }

/*     determine number of prior modes available at restart */

    nlock = 0;
    while(TRUE_) {
	io___55.ciunit = ivb1;
	i__1 = s_rsue(&io___55);
	if (i__1 != 0) {
	    goto L130;
	}
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = do_uio(&c__1, (char *)&p[k - 1], (ftnlen)sizeof(doublereal)
		    );
	    if (i__1 != 0) {
		goto L130;
	    }
	}
	i__1 = e_rsue();
	if (i__1 != 0) {
	    goto L130;
	}
	++nlock;
    }
L130:
    al__1.aerr = 0;
    al__1.aunit = ivb1;
    f_rew(&al__1);
    if (nlock != 0) {
	io___57.ciunit = iounit_1.iout;
	s_wsfe(&io___57);
	do_fio(&c__1, (char *)&nlock, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    nconv = nlock;

/*     compute and diagonalize the Hessian for each block */

    k0 = 0;
    i1 = 1;
    i__1 = nblk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ > 1) {
/* Computing 2nd power */
	    i__2 = iblk[i__ - 2];
	    k0 += i__2 * i__2 * 9;
	    i1 += iblk[i__ - 2];
	}
	i2 = i1 + iblk[i__ - 1] - 1;
	k1 = i1 * 3 - 2;
	k2 = i2 * 3;
	hessblk_(atmtyp_1.mass, &k0, &i1, &i2, uu);
	i__2 = iblk[i__ - 1] * 3;
	diagblk_(&k0, &k1, &i__2, uu, uku);
    }

/*     use negative of eigenvalues if doing high frequencies */

    i__1 = nvar;
    for (k = 1; k <= i__1; ++k) {
	uku[k - 1] = factor * uku[k - 1];
	uku0[k - 1] = uku[k - 1];
    }
    uku_max__ = uku[0];
    uku_min__ = uku[0];
    i__1 = nvar;
    for (k = 2; k <= i__1; ++k) {
	if (uku[k - 1] > uku_max__) {
	    uku_max__ = uku[k - 1];
	}
	if (uku[k - 1] < uku_min__) {
	    uku_min__ = uku[k - 1];
	}
    }

/*     if restarting, read trial vectors and estimate eigenvalues */

    if (restart) {
	i__1 = npair;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___68.ciunit = ivb2;
	    s_rsue(&io___68);
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		do_uio(&c__1, (char *)&phi_ref(k, i__), (ftnlen)sizeof(
			doublereal));
	    }
	    e_rsue();
	    io___70.ciunit = ivb2;
	    s_rsue(&io___70);
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		do_uio(&c__1, (char *)&phik_ref(k, i__), (ftnlen)sizeof(
			doublereal));
	    }
	    e_rsue();
	}
	i__1 = nroot;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    h___ref(i__, i__) = 0.;
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		h___ref(i__, i__) = h___ref(i__, i__) + phik_ref(k, i__) * 
			phi_ref(k, i__);
	    }
	    freqold[i__ - 1] = d_sign(&c_b75, &h___ref(i__, i__)) * sqrt((
		    d__1 = h___ref(i__, i__), abs(d__1)));
	}
	goto L150;
    }

/*     if not restarting, generate initial guess eigenvectors */

    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	trigger_(&nvar, &nbasis, &np, &ifactor, &nblk, iblk, u, uu, p);
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phi_ref(k, i__) = p[k - 1];
	}
    }

/*     project out locked roots from components of phi */

    project_(&nvar, &nconv, &ivb1, &nroot, &c__0, phi);
    project_(&nvar, &nconv, &ivb1, &nroot, &c__0, phik);

/*     reload and make vector orthonormal to existing basis */

    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    p[k - 1] = phi_ref(k, i__);
	}
	if (i__ == 1) {
	    sum = 0.;
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		sum += p[k - 1] * p[k - 1];
	    }
	    sum = sqrt(sum);
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		p[k - 1] /= sum;
	    }
	} else {
	    i__2 = i__ - 1;
	    gsort_(&nvar, &i__2, phi, p);
	}
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phi_ref(k, i__) = p[k - 1];
	}
    }

/*     store K on phi */

    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    p[k - 1] = phi_ref(k, i__);
	}
	konvec_(&nvar, xm, xe, p, pk);
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phik_ref(k, i__) = factor * pk[k - 1];
	}
    }

/*     make nroot-by-nroot CI matrix */

    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nroot;
	for (j = i__; j <= i__2; ++j) {
	    h___ref(i__, j) = 0.;
	    i__4 = nvar;
	    for (k = 1; k <= i__4; ++k) {
		h___ref(i__, j) = h___ref(i__, j) + phik_ref(k, i__) * 
			phi_ref(k, j);
	    }
	    h___ref(j, i__) = h___ref(i__, j);
	}
    }

/*     diagonalize and use first nroot solutions as starting basis */

    transform_(&nroot, h__, c__);

/*     fill up arrays */

    i__1 = nvar;
    for (k = 1; k <= i__1; ++k) {
	i__2 = nroot;
	for (j = 1; j <= i__2; ++j) {
	    tmp1[j - 1] = 0.;
	    tmp2[j - 1] = 0.;
	    i__4 = nroot;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		tmp1[j - 1] += c___ref(i__, j) * phi_ref(k, i__);
		tmp2[j - 1] += c___ref(i__, j) * phik_ref(k, i__);
	    }
	}
	i__2 = nroot;
	for (j = 1; j <= i__2; ++j) {
	    phi_ref(k, j) = tmp1[j - 1];
	    phik_ref(k, j) = tmp2[j - 1];
	}
    }

/*     residues of guesses */

    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	freq[i__ - 1] = funit * d_sign(&c_b75, &h___ref(i__, i__)) * sqrt((
		d__1 = h___ref(i__, i__), abs(d__1)));
	freqold[i__ - 1] = freq[i__ - 1];
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    pk[k - 1] = phik_ref(k, i__) - h___ref(i__, i__) * phi_ref(k, i__)
		    ;
	}

/*     use Davidson preconditioner if finding low frequencies */

	if (factor > 0.) {
	    preconblk_(&nvar, &nblk, iblk, uku, uu, &h___ref(i__, i__), &hmin[
		    i__ - 1], pk);
	}

/*     project residual onto P-space */

	qonvec_(&nvar, &np, u, pk, p);
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phi_ref(k, i__ + nroot) = p[k - 1];
	}
    }

/*     project out locked roots from components of phi */

    project_(&nvar, &nconv, &ivb1, &nroot, &nroot, phi);
    project_(&nvar, &nconv, &ivb1, &nroot, &nroot, phik);

/*     reload and make vector orthonormal to existing basis */

    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    p[k - 1] = phi_ref(k, i__ + nroot);
	}
	i__2 = nroot + i__ - 1;
	gsort_(&nvar, &i__2, phi, p);
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phi_ref(k, i__ + nroot) = p[k - 1];
	}
    }

/*     store K on phi */

    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    p[k - 1] = phi_ref(k, i__ + nroot);
	}
	konvec_(&nvar, xm, xe, p, pk);
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phik_ref(k, i__ + nroot) = factor * pk[k - 1];
	}
    }

/*     make npair-by-npair CI matrix */

    i__1 = npair;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = npair;
	for (j = i__; j <= i__2; ++j) {
	    h___ref(i__, j) = 0.;
	    i__4 = nvar;
	    for (k = 1; k <= i__4; ++k) {
		h___ref(i__, j) = h___ref(i__, j) + phik_ref(k, i__) * 
			phi_ref(k, j);
	    }
	    h___ref(j, i__) = h___ref(i__, j);
	}
    }

/*     diagonalize and use first nroot solutions as new guess */

    transform_(&npair, h__, c__);
    i__1 = nvar;
    for (k = 1; k <= i__1; ++k) {
	i__2 = nroot;
	for (j = 1; j <= i__2; ++j) {
	    tmp1[j - 1] = 0.;
	    tmp2[j - 1] = 0.;
	    i__4 = npair;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		tmp1[j - 1] += c___ref(i__, j) * phi_ref(k, i__);
		tmp2[j - 1] += c___ref(i__, j) * phik_ref(k, i__);
	    }
	}

/*     old solution fills up 2nd block */

	i__2 = nroot;
	for (j = 1; j <= i__2; ++j) {
	    phi_ref(k, j + nroot) = phi_ref(k, j);
	    phik_ref(k, j + nroot) = phik_ref(k, j);
	}

/*     new solution fills up 1st block */

	i__2 = nroot;
	for (j = 1; j <= i__2; ++j) {
	    phi_ref(k, j) = tmp1[j - 1];
	    phik_ref(k, j) = tmp2[j - 1];
	}

/*     orthogonalize 2nd block to 1st */

	i__2 = nroot;
	for (j = 1; j <= i__2; ++j) {
	    i__4 = nroot;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		phi_ref(k, j + nroot) = phi_ref(k, j + nroot) - c___ref(j, 
			i__) * phi_ref(k, i__);
		phik_ref(k, j + nroot) = phik_ref(k, j + nroot) - c___ref(j, 
			i__) * phik_ref(k, i__);
	    }
	}
    }

/*     orthogonalize 2nd block on itself */

    sum = 0.;
    i__1 = nvar;
    for (k = 1; k <= i__1; ++k) {
	sum += phi_ref(k, nroot + 1) * phi_ref(k, nroot + 1);
    }
    sum = sqrt(sum);

/*     normalize leading vector */

    i__1 = nvar;
    for (k = 1; k <= i__1; ++k) {
	phi_ref(k, nroot + 1) = phi_ref(k, nroot + 1) / sum;
	phik_ref(k, nroot + 1) = phik_ref(k, nroot + 1) / sum;
    }

/*     orthogonalize the rest one-by-one */

    if (nroot > 1) {
	i__1 = max(2,nroot);
	for (i__ = 2; i__ <= i__1; ++i__) {
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		sum = 0.;
		i__4 = nvar;
		for (k = 1; k <= i__4; ++k) {
		    sum += phi_ref(k, i__ + nroot) * phi_ref(k, j + nroot);
		}
		i__4 = nvar;
		for (k = 1; k <= i__4; ++k) {
		    phi_ref(k, i__ + nroot) = phi_ref(k, i__ + nroot) - sum * 
			    phi_ref(k, j + nroot);
		    phik_ref(k, i__ + nroot) = phik_ref(k, i__ + nroot) - sum 
			    * phik_ref(k, j + nroot);
		}
	    }
	    sum = 0.;
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		sum += phi_ref(k, i__ + nroot) * phi_ref(k, i__ + nroot);
	    }
	    sum = sqrt(sum);
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		phi_ref(k, i__ + nroot) = phi_ref(k, i__ + nroot) / sum;
		phik_ref(k, i__ + nroot) = phik_ref(k, i__ + nroot) / sum;
	    }
	}
    }

/*     residue of new solution (if restarting, begin here) */

L150:
    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	freq[i__ - 1] = funit * d_sign(&c_b75, &h___ref(i__, i__)) * sqrt((
		d__1 = h___ref(i__, i__), abs(d__1)));
	freq[i__ + nroot - 1] = funit * d_sign(&c_b75, &h___ref(i__ + nroot, 
		i__ + nroot)) * sqrt((d__1 = h___ref(i__ + nroot, i__ + nroot)
		, abs(d__1)));
	freq[i__ + npair - 1] = funit * d_sign(&c_b75, &h___ref(i__ + npair, 
		i__ + npair)) * sqrt((d__1 = h___ref(i__ + npair, i__ + npair)
		, abs(d__1)));
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    pk[k - 1] = phik_ref(k, i__) - h___ref(i__, i__) * phi_ref(k, i__)
		    ;
	}

/*     use Davidson preconditioner if finding low frequencies */

	if (factor > 0.) {
	    preconblk_(&nvar, &nblk, iblk, uku, uu, &h___ref(i__, i__), &hmin[
		    i__ - 1], pk);
	}

/*     project onto P-space */

	qonvec_(&nvar, &np, u, pk, p);
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phi_ref(k, i__ + npair) = p[k - 1];
	}
    }

/*     project out locked roots from components of phi */

    project_(&nvar, &nconv, &ivb1, &nroot, &npair, phi);
    project_(&nvar, &nconv, &ivb1, &nroot, &npair, phik);

/*     reload and orthogonalize to 1st + 2nd */

    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    p[k - 1] = phi_ref(k, i__ + npair);
	}
	i__2 = npair + i__ - 1;
	gsort_(&nvar, &i__2, phi, p);
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phi_ref(k, i__ + npair) = p[k - 1];
	}
    }

/*     store K on phi */

    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    p[k - 1] = phi_ref(k, i__ + npair);
	}
	konvec_(&nvar, xm, xe, p, pk);
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phik_ref(k, i__ + npair) = factor * pk[k - 1];
	}
    }

/*     beginning of iterations */

    iconv = 0;
L160:
    done = FALSE_;
    ++iter;

/*     make nbasis-by-nbasis CI matrix */

    i__1 = nbasis;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nbasis;
	for (j = i__; j <= i__2; ++j) {
	    h___ref(i__, j) = 0.;
	    i__4 = nvar;
	    for (k = 1; k <= i__4; ++k) {
		h___ref(i__, j) = h___ref(i__, j) + phik_ref(k, i__) * 
			phi_ref(k, j);
	    }
	    h___ref(j, i__) = h___ref(i__, j);
	}
    }

/*     list of previous frequencies */

    i__1 = npair;
    for (i__ = 1; i__ <= i__1; ++i__) {
	freqold[i__ - 1] = freq[i__ - 1];
    }

/*     diagonalize and use first nroot solutions as new guess */

    transform_(&nbasis, h__, c__);

/*     check for collapse based on leading component of ground state */

    if (iconv == 0 && nconv > 0) {
/* Computing 2nd power */
	d__1 = c___ref(1, 1);
	sum = sqrt(1. - d__1 * d__1);
	if (sum > .9) {
	    io___83.ciunit = iounit_1.iout;
	    s_wsfe(&io___83);
	    i__1 = nconv - nlock;
	    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___84.ciunit = iounit_1.iout;
	    s_wsfe(&io___84);
	    e_wsfe();
	    cl__1.cerr = 0;
	    cl__1.cunit = ivb2;
	    cl__1.csta = "delete";
	    f_clos(&cl__1);
	    goto L280;
	}
    }

/*     list of new frequencies */

    i__1 = npair;
    for (i__ = 1; i__ <= i__1; ++i__) {
	freq[i__ - 1] = funit * d_sign(&c_b75, &h___ref(i__, i__)) * sqrt((
		d__1 = h___ref(i__, i__), abs(d__1)));
    }

/*     check if first few have converged */

    iconv = 0;
L190:
    dfreq = freqold[iconv] - freq[iconv];
    if (dfreq * factor > 0. && dfreq * factor < wtol) {
	++iconv;
	goto L190;
    }
    if (iconv > 0) {

/*     shift levels of preconditioner matrix; since the Hessian */
/*     is gradually deflated, reduce effect of the preconditioner */
/*     based on a simple 1/x curve, the uku levels are squeezed */
/*     upwards to eventually lead to a unit operator */

	ratio = (doublereal) (nconv + iconv) / (doublereal) nvar;
	shift = uku_min__ / (1. - ratio);
	shift += h___ref(iconv + nroot, iconv + nroot);

/*     do a regular shift, which also seems to work */

	i__1 = nvar;
	for (k = 1; k <= i__1; ++k) {
	    uku[k - 1] = uku_max__ + (uku0[k - 1] - uku_max__) * (uku_max__ - 
		    shift) / (uku_max__ - uku_min__);
	}

/*     move cursor to end of storage file */

	i__1 = nconv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___88.ciunit = ivb1;
	    s_rsue(&io___88);
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		do_uio(&c__1, (char *)&pk[k - 1], (ftnlen)sizeof(doublereal));
	    }
	    e_rsue();
	}

/*     norm of residual */

	i__1 = iconv;
	for (j = 1; j <= i__1; ++j) {
	    rnorm = 0.;
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		p[k - 1] = 0.;
		pk[k - 1] = 0.;
		i__4 = nbasis;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    p[k - 1] += c___ref(i__, j) * phi_ref(k, i__);
		    pk[k - 1] += c___ref(i__, j) * phik_ref(k, i__);
		}
/* Computing 2nd power */
		d__1 = pk[k - 1] - h___ref(j, j) * p[k - 1];
		rnorm += d__1 * d__1;
	    }
	    rnorm = sqrt(rnorm);

/*     component of root in R-space */

	    for (i__ = 1; i__ <= 3; ++i__) {
		tmp1[i__ - 1] = 0.;
		i__2 = nvar;
		for (k = 1; k <= i__2; ++k) {
		    tmp1[i__ - 1] += ur_ref(k, i__) * p[k - 1];
		}
	    }
	    rcomp = 0.;
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		sum = 0.;
		for (i__ = 1; i__ <= 3; ++i__) {
		    sum += ur_ref(k, i__) * tmp1[i__ - 1];
		}
		rcomp += sum * sum;
	    }
	    rcomp = sqrt(rcomp);

/*     write the converged mode to formatted and binary files */

	    ivib = irange + ifactor * (nconv + j);
	    if ((header || inform_1.verbose) && j == 1) {
		header = FALSE_;
		io___92.ciunit = iounit_1.iout;
		s_wsfe(&io___92);
		e_wsfe();
		io___93.ciunit = iounit_1.iout;
		s_wsfe(&io___93);
		e_wsfe();
		if (! inform_1.verbose) {
		    io___94.ciunit = iounit_1.iout;
		    s_wsfe(&io___94);
		    e_wsfe();
		}
	    }
	    dfreq = freqold[j - 1] - freq[j - 1];
	    io___95.ciunit = iounit_1.iout;
	    s_wsfe(&io___95);
	    do_fio(&c__1, (char *)&ivib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&freq[j - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&dfreq, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&rnorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&rcomp, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    prtvib_(&ivib, p);
	    io___96.ciunit = ivb1;
	    s_wsue(&io___96);
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		do_uio(&c__1, (char *)&p[k - 1], (ftnlen)sizeof(doublereal));
	    }
	    e_wsue();
	}
	al__1.aerr = 0;
	al__1.aunit = ivb1;
	f_rew(&al__1);

/*     update total number of vectors locked on disk */

	nconv += iconv;
	if (freq[iconv - 1] * factor >= fmax * factor) {
	    done = TRUE_;
	    cl__1.cerr = 0;
	    cl__1.cunit = ivb1;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	}
    }

/*     shift frequency arrays by iconv */

    i__1 = npair;
    for (i__ = 1; i__ <= i__1; ++i__) {
	freq[i__ - 1] = freq[i__ + iconv - 1];
	freqold[i__ - 1] = freqold[i__ + iconv - 1];
    }
    i__1 = nvar;
    for (k = 1; k <= i__1; ++k) {
	i__2 = nroot + iconv;
	for (j = 1; j <= i__2; ++j) {
	    tmp1[j - 1] = 0.;
	    tmp2[j - 1] = 0.;
	    i__4 = nbasis;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		tmp1[j - 1] += c___ref(i__, j) * phi_ref(k, i__);
		tmp2[j - 1] += c___ref(i__, j) * phik_ref(k, i__);
	    }
	}

/*     old solution fills up 2nd block */

	i__2 = nroot;
	for (j = 1; j <= i__2; ++j) {
	    phi_ref(k, j + nroot + iconv) = phi_ref(k, j + iconv);
	    phik_ref(k, j + nroot + iconv) = phik_ref(k, j + iconv);
	}

/*     new solution fills up 1st block */

	i__2 = nroot;
	for (j = 1; j <= i__2; ++j) {
	    phi_ref(k, j + iconv) = tmp1[j + iconv - 1];
	    phik_ref(k, j + iconv) = tmp2[j + iconv - 1];
	}

/*     shift index down by iconv */

	i__2 = npair;
	for (j = 1; j <= i__2; ++j) {
	    phi_ref(k, j) = phi_ref(k, j + iconv);
	    phik_ref(k, j) = phik_ref(k, j + iconv);
	}

/*     orthogonalize 2nd block to 1st + iconv roots */

	i__2 = nroot;
	for (j = 1; j <= i__2; ++j) {
	    i__4 = nroot;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		phi_ref(k, j + nroot) = phi_ref(k, j + nroot) - c___ref(j + 
			iconv, i__ + iconv) * phi_ref(k, i__);
		phik_ref(k, j + nroot) = phik_ref(k, j + nroot) - c___ref(j + 
			iconv, i__ + iconv) * phik_ref(k, i__);
	    }
	    i__4 = iconv;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		phi_ref(k, j + nroot) = phi_ref(k, j + nroot) - c___ref(j + 
			iconv, i__) * tmp1[i__ - 1];
		phik_ref(k, j + nroot) = phik_ref(k, j + nroot) - c___ref(j + 
			iconv, i__) * tmp2[i__ - 1];
	    }
	}
    }

/*     orthogonalize 2nd block on itself */

    sum = 0.;
    i__1 = nvar;
    for (k = 1; k <= i__1; ++k) {
	sum += phi_ref(k, nroot + 1) * phi_ref(k, nroot + 1);
    }
    sum = sqrt(sum);

/*     normalize leading vector */

    i__1 = nvar;
    for (k = 1; k <= i__1; ++k) {
	phi_ref(k, nroot + 1) = phi_ref(k, nroot + 1) / sum;
	phik_ref(k, nroot + 1) = phik_ref(k, nroot + 1) / sum;
    }

/*     orthogonalize the rest one-by-one */

    if (nroot > 1) {
	i__1 = max(2,nroot);
	for (i__ = 2; i__ <= i__1; ++i__) {
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		sum = 0.;
		i__4 = nvar;
		for (k = 1; k <= i__4; ++k) {
		    sum += phi_ref(k, i__ + nroot) * phi_ref(k, j + nroot);
		}
		i__4 = nvar;
		for (k = 1; k <= i__4; ++k) {
		    phi_ref(k, i__ + nroot) = phi_ref(k, i__ + nroot) - sum * 
			    phi_ref(k, j + nroot);
		    phik_ref(k, i__ + nroot) = phik_ref(k, i__ + nroot) - sum 
			    * phik_ref(k, j + nroot);
		}
	    }
	    sum = 0.;
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		sum += phi_ref(k, i__ + nroot) * phi_ref(k, i__ + nroot);
	    }
	    sum = sqrt(sum);
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		phi_ref(k, i__ + nroot) = phi_ref(k, i__ + nroot) / sum;
		phik_ref(k, i__ + nroot) = phik_ref(k, i__ + nroot) / sum;
	    }
	}
    }

/*     print a header for the current iteration */

    if (inform_1.verbose) {
	io___97.ciunit = iounit_1.iout;
	s_wsfe(&io___97);
	do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iconv, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nconv, (ftnlen)sizeof(integer));
	e_wsfe();
	io___98.ciunit = iounit_1.iout;
	s_wsfe(&io___98);
	e_wsfe();
    }

/*     norm of residual */

    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rnorm = 0.;
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    d__1 = phik_ref(k, i__) - h___ref(i__ + iconv, i__ + iconv) * 
		    phi_ref(k, i__);
	    rnorm += d__1 * d__1;
	}
	rnorm = sqrt(rnorm);

/*     calculate root's component in R-space */

	for (j = 1; j <= 3; ++j) {
	    tmp1[j - 1] = 0.;
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		tmp1[j - 1] += ur_ref(k, j) * phi_ref(k, i__);
	    }
	}
	rcomp = 0.;
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    sum = 0.;
	    for (j = 1; j <= 3; ++j) {
		sum += ur_ref(k, j) * tmp1[j - 1];
	    }
	    rcomp += sum * sum;
	}
	rcomp = sqrt(rcomp);
	dfreq = freqold[i__ - 1] - freq[i__ - 1];
	if (inform_1.verbose) {
	    io___99.ciunit = iounit_1.iout;
	    s_wsfe(&io___99);
	    i__2 = irange + ifactor * (i__ + nconv);
	    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&freq[i__ - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&dfreq, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&rnorm, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&rcomp, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     save vectors for restart */

    if (iter % idump == 0) {
	al__1.aerr = 0;
	al__1.aunit = ivb2;
	f_rew(&al__1);
	i__1 = npair;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___100.ciunit = ivb2;
	    s_wsue(&io___100);
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		do_uio(&c__1, (char *)&phi_ref(k, i__), (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsue();
	    io___101.ciunit = ivb2;
	    s_wsue(&io___101);
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		do_uio(&c__1, (char *)&phik_ref(k, i__), (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsue();
	}
    }

/*     prepare restart if finished or iterations exhausted */

    if (done || iter == maxiter) {
	io___102.ciunit = iounit_1.iout;
	s_wsfe(&io___102);
	i__1 = nconv - nlock;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	e_wsfe();
	al__1.aerr = 0;
	al__1.aunit = ivb2;
	f_rew(&al__1);
	i__1 = npair;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___103.ciunit = ivb2;
	    s_wsue(&io___103);
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		do_uio(&c__1, (char *)&phi_ref(k, i__), (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsue();
	    io___104.ciunit = ivb2;
	    s_wsue(&io___104);
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		do_uio(&c__1, (char *)&phik_ref(k, i__), (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsue();
	}
	cl__1.cerr = 0;
	cl__1.cunit = ivb2;
	cl__1.csta = 0;
	f_clos(&cl__1);
	goto L280;
    }

/*     as above, make sure no prior roots are mixed in the basis */

    i__1 = npair;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    p[k - 1] = phi_ref(k, i__);
	}
	qonvec_(&nvar, &np, u, p, pk);
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phi_ref(k, i__) = pk[k - 1];
	}
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    p[k - 1] = phik_ref(k, i__);
	}
	qonvec_(&nvar, &np, u, p, pk);
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phik_ref(k, i__) = pk[k - 1];
	}
    }

/*     project out locked roots from components of phi */

    project_(&nvar, &nconv, &ivb1, &npair, &c__0, phi);
    project_(&nvar, &nconv, &ivb1, &npair, &c__0, phik);

/*     setup next iteration; solution residue, Davidson weight */

    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    pk[k - 1] = phik_ref(k, i__) - h___ref(i__ + iconv, i__ + iconv) *
		     phi_ref(k, i__);
	}

/*     use Davidson preconditioner if finding low frequencies */

	ii = i__ + iconv;
	if (factor > 0.) {
	    preconblk_(&nvar, &nblk, iblk, uku, uu, &h___ref(ii, ii), &hmin[
		    i__ - 1], pk);
	}

/*     project residual onto P-space */

	qonvec_(&nvar, &np, u, pk, p);
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phi_ref(k, i__ + npair) = p[k - 1];
	}
    }

/*     project out locked roots from components of phi */

    project_(&nvar, &nconv, &ivb1, &nroot, &npair, phi);

/*     reload and orthogonalize to 1st + 2nd */

    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    p[k - 1] = phi_ref(k, i__ + npair);
	}
	i__2 = npair + i__ - 1;
	gsort_(&nvar, &i__2, phi, p);
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phi_ref(k, i__ + npair) = p[k - 1];
	}
    }

/*     store K on phi */

    i__1 = nroot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    p[k - 1] = phi_ref(k, i__ + npair);
	}
	konvec_(&nvar, xm, xe, p, pk);
	qonvec_(&nvar, &np, u, pk, p);
	i__2 = nvar;
	for (k = 1; k <= i__2; ++k) {
	    phik_ref(k, i__ + npair) = factor * p[k - 1];
	}
    }

/*     project out locked roots from components of phik */

    project_(&nvar, &nconv, &ivb1, &nroot, &npair, phik);
    goto L160;
L280:
    return 0;
} /* MAIN__ */

#undef keyline_ref
#undef phik_ref
#undef phi_ref
#undef ur_ref
#undef h___ref
#undef c___ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine trigger  --  get initial trial eigenvectors  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "trigger" constructs a set of initial trial vectors for */
/*     use during sliding block iterative matrix diagonalization */


/* Subroutine */ int trigger_(integer *nvar, integer *nbasis, integer *np, 
	integer *ifactor, integer *nblk, integer *iblk, doublereal *u, 
	doublereal *uu, doublereal *p)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal w;
    static integer k0, k1, k2;
    static doublereal tmp[75000], sum;
    extern doublereal random_(void);
    extern /* Subroutine */ int qonvec_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static integer nguess;



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




/*     set the number of random guesses */

    /* Parameter adjustments */
    --p;
    --uu;
    u -= 75001;
    --iblk;

    /* Function Body */
    nguess = (integer) ((doublereal) (*nbasis) / (doublereal) (*nblk)) + 1;

/*     zero out the trial vector */

    i__1 = *nvar;
    for (k = 1; k <= i__1; ++k) {
	p[k] = 0.;
    }

/*     create overlap with the entire P-space */

    k0 = 0;
    k1 = 1;
    i__1 = *nblk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ > 1) {
/* Computing 2nd power */
	    i__2 = iblk[i__ - 1];
	    k0 += i__2 * i__2 * 9;
	    k1 += iblk[i__ - 1] * 3;
	}
	k2 = k1 + iblk[i__] * 3 - 1;

/*     scan over rows of the Hessian */

	m = 0;
	i__2 = iblk[i__] * 3;
	for (j = 1; j <= i__2; ++j) {
	    if (*ifactor == 1) {
/* Computing MIN */
		i__3 = nguess, i__4 = iblk[i__] * 3;
		if (j > min(i__3,i__4)) {
		    w = 0.;
		} else {
		    w = random_() - .5;
		}
	    } else {
/* Computing MIN */
		i__3 = nguess, i__4 = iblk[i__] * 3;
		if (j < iblk[i__] * 3 - min(i__3,i__4) + 1) {
		    w = 0.;
		} else {
		    w = random_() - .5;
		}
	    }
	    i__3 = k2;
	    for (k = k1; k <= i__3; ++k) {
		++m;
		p[k] += w * uu[k0 + m];
	    }
	}
    }

/*     project the vector onto P-space */

    qonvec_(nvar, np, &u[75001], &p[1], tmp);

/*     perform a normalization */

    sum = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = tmp[i__ - 1];
	sum += d__1 * d__1;
    }
    sum = sqrt(sum);
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = tmp[i__ - 1] / sum;
    }
    return 0;
} /* trigger_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine trbasis  --  set translation/rotation vectors  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "trbasis" forms translation and rotation basis vectors used */
/*     during vibrational analysis via block iterative diagonalization */


/* Subroutine */ int trbasis_(integer *np, doublereal *xe, doublereal *u, 
	doublereal *ur)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__[9]	/* was [3][3] */, e[9]	/* was [3][3] */;
    static integer i__, j, k;
    static doublereal p[3], t1[3], t2[3], cm[3], ra, pr, rha, sum;
    static integer nvar;
    static doublereal tmass;
    extern /* Subroutine */ int jacobi_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


#define c___ref(a_1,a_2) c__[(a_2)*3 + a_1 - 4]
#define e_ref(a_1,a_2) e[(a_2)*3 + a_1 - 4]
#define u_ref(a_1,a_2) u[(a_2)*75000 + a_1]
#define ur_ref(a_1,a_2) ur[(a_2)*75000 + a_1]



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




/*     zero out the translation and rotation vectors */

    /* Parameter adjustments */
    ur -= 75001;
    u -= 75001;
    --xe;

    /* Function Body */
    nvar = atoms_1.n * 3;
    for (i__ = 1; i__ <= 6; ++i__) {
	i__1 = nvar;
	for (j = 1; j <= i__1; ++j) {
	    u_ref(j, i__) = 0.;
	}
    }

/*     get the total mass of the system */

    tmass = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tmass += atmtyp_1.mass[i__ - 1];
    }

/*     set basis vectors for translations */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	u_ref(i__ * 3 - 2, 1) = sqrt(atmtyp_1.mass[i__ - 1] / tmass);
	u_ref(i__ * 3 - 1, 1) = 0.;
	u_ref(i__ * 3, 1) = 0.;
	u_ref(i__ * 3 - 2, 2) = 0.;
	u_ref(i__ * 3 - 1, 2) = sqrt(atmtyp_1.mass[i__ - 1] / tmass);
	u_ref(i__ * 3, 2) = 0.;
	u_ref(i__ * 3 - 2, 3) = 0.;
	u_ref(i__ * 3 - 1, 3) = 0.;
	u_ref(i__ * 3, 3) = sqrt(atmtyp_1.mass[i__ - 1] / tmass);
    }

/*     move center of mass to origin */

    for (i__ = 1; i__ <= 3; ++i__) {
	cm[i__ - 1] = 0.;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    cm[j - 1] += xe[(i__ - 1) * 3 + j] * atmtyp_1.mass[i__ - 1];
	}
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    xe[(i__ - 1) * 3 + j] -= cm[j - 1] / tmass;
	}
    }

/*     get the moments of inertia */

    for (i__ = 1; i__ <= 3; ++i__) {
	e_ref(i__, i__) = 0.;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = xe[i__ * 3 - 1];
/* Computing 2nd power */
	d__2 = xe[i__ * 3];
	e_ref(1, 1) = e_ref(1, 1) + (d__1 * d__1 + d__2 * d__2) * 
		atmtyp_1.mass[i__ - 1];
/* Computing 2nd power */
	d__1 = xe[i__ * 3 - 2];
/* Computing 2nd power */
	d__2 = xe[i__ * 3];
	e_ref(2, 2) = e_ref(2, 2) + (d__1 * d__1 + d__2 * d__2) * 
		atmtyp_1.mass[i__ - 1];
/* Computing 2nd power */
	d__1 = xe[i__ * 3 - 2];
/* Computing 2nd power */
	d__2 = xe[i__ * 3 - 1];
	e_ref(3, 3) = e_ref(3, 3) + (d__1 * d__1 + d__2 * d__2) * 
		atmtyp_1.mass[i__ - 1];
    }
    for (i__ = 1; i__ <= 2; ++i__) {
	for (j = i__ + 1; j <= 3; ++j) {
	    e_ref(i__, j) = 0.;
	    i__1 = atoms_1.n;
	    for (k = 1; k <= i__1; ++k) {
		e_ref(i__, j) = e_ref(i__, j) - xe[(k - 1) * 3 + i__] * xe[(k 
			- 1) * 3 + j] * atmtyp_1.mass[k - 1];
	    }
	    e_ref(j, i__) = e_ref(i__, j);
	}
    }

/*     diagonalize to get principal axes */

    jacobi_(&c__3, &c__3, e, cm, c__, t1, t2);

/*     construction of principle rotations */

    for (i__ = 1; i__ <= 3; ++i__) {
	i__1 = atoms_1.n;
	for (j = 1; j <= i__1; ++j) {
	    ra = 0.;
	    pr = 0.;
	    for (k = 1; k <= 3; ++k) {
		cm[k - 1] = xe[(j - 1) * 3 + k];
/* Computing 2nd power */
		d__1 = cm[k - 1];
		ra += d__1 * d__1;
		pr += cm[k - 1] * c___ref(k, i__);
	    }
/* Computing 2nd power */
	    d__1 = pr;
	    rha = sqrt(ra - d__1 * d__1);
	    p[0] = c___ref(2, i__) * cm[2] - c___ref(3, i__) * cm[1];
	    p[1] = c___ref(3, i__) * cm[0] - c___ref(1, i__) * cm[2];
	    p[2] = c___ref(1, i__) * cm[1] - c___ref(2, i__) * cm[0];
	    sum = 0.;
	    for (k = 1; k <= 3; ++k) {
/* Computing 2nd power */
		d__1 = p[k - 1];
		sum += d__1 * d__1;
	    }
	    sum = sqrt(sum);
	    for (k = 1; k <= 3; ++k) {
		ur_ref((j - 1) * 3 + k, i__) = sqrt(atmtyp_1.mass[j - 1]) * 
			rha * p[k - 1] / sum;
	    }
	}
	sum = 0.;
	i__1 = nvar;
	for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	    d__1 = ur_ref(j, i__);
	    sum += d__1 * d__1;
	}
	sum = sqrt(sum);
	i__1 = nvar;
	for (j = 1; j <= i__1; ++j) {
	    ur_ref(j, i__) = ur_ref(j, i__) / sum;
	}
    }

/*     set basis vectors for rotation */

    if (*np == 6) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    i__1 = nvar;
	    for (j = 1; j <= i__1; ++j) {
		u_ref(j, i__ + 3) = ur_ref(j, i__);
	    }
	}
    }
    return 0;
} /* trbasis_ */

#undef ur_ref
#undef u_ref
#undef e_ref
#undef c___ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine preconblk  --  precondition atom block Hessian  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "preconblk" applies a preconditioner to an atom block section */
/*     of the Hessian matrix */


/* Subroutine */ int preconblk_(integer *nvar, integer *nblk, integer *iblk, 
	doublereal *uku, doublereal *uu, doublereal *h__, doublereal *hmin, 
	doublereal *pk)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal d__[75000];
    static integer i__, j, k, l, k0, k1, k2, l2;
    static doublereal work[75000];



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




/*     find smallest element of |h-uku| */

    /* Parameter adjustments */
    --pk;
    --uu;
    --uku;
    --iblk;

    /* Function Body */
    *hmin = (d__1 = *h__ - uku[1], abs(d__1));
    i__1 = *nvar;
    for (k = 2; k <= i__1; ++k) {
	if ((d__1 = *h__ - uku[k], abs(d__1)) < *hmin) {
	    *hmin = (d__1 = *h__ - uku[k], abs(d__1));
	}
    }

/*     assign values to temporary array */

    i__1 = *nvar;
    for (k = 1; k <= i__1; ++k) {
	d__[k - 1] = *h__ - uku[k];
    }

/*     invert array via d=hmin/d, where hmin=min{|d(k)|} */

    i__1 = *nvar;
    for (k = 1; k <= i__1; ++k) {
	d__[k - 1] = *hmin / d__[k - 1];
    }

/*     create overlap with the entire pk array */

    k0 = 0;
    k1 = 1;
    i__1 = *nblk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ > 1) {
/* Computing 2nd power */
	    i__2 = iblk[i__ - 1];
	    k0 += i__2 * i__2 * 9;
	    k1 += iblk[i__ - 1] * 3;
	}
	k2 = k1 + iblk[i__] * 3 - 1;

/*    scan over rows of the Hessian, first part */

	l = 0;
	i__2 = iblk[i__] * 3;
	for (j = 1; j <= i__2; ++j) {
	    l2 = k1 + j - 1;
	    work[l2 - 1] = 0.;
	    i__3 = k2;
	    for (k = k1; k <= i__3; ++k) {
		++l;
		work[l2 - 1] += uu[k0 + l] * pk[k];
	    }
	}

/*    zero out the segment */

	i__2 = k2;
	for (k = k1; k <= i__2; ++k) {
	    pk[k] = 0.;
	}

/*    scan over rows of the Hessian, second part */

	l = 0;
	i__2 = iblk[i__] * 3;
	for (j = 1; j <= i__2; ++j) {
	    l2 = k1 + j - 1;
	    i__3 = k2;
	    for (k = k1; k <= i__3; ++k) {
		++l;
		pk[k] += uu[k0 + l] * d__[l2 - 1] * work[l2 - 1];
	    }
	}
    }
    return 0;
} /* preconblk_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine gsort  --  orthogonal vector via Gram-Schmidt  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "gsort" uses the Gram-Schmidt algorithm to build orthogonal */
/*     vectors for sliding block interative matrix diagonalization */


/* Subroutine */ int gsort_(integer *nvar, integer *nb, doublereal *phi, 
	doublereal *p0)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal s[150], sum, proj[75000];


#define phi_ref(a_1,a_2) phi[(a_2)*75000 + a_1]



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




/*     make overlap between two basis sets */

    /* Parameter adjustments */
    --p0;
    phi -= 75001;

    /* Function Body */
    i__1 = *nb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s[i__ - 1] = 0.;
	i__2 = *nvar;
	for (j = 1; j <= i__2; ++j) {
	    s[i__ - 1] += p0[j] * phi_ref(j, i__);
	}
    }

/*     start the Gram-Schmidt procedure */

    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	proj[i__ - 1] = 0.;
    }

/*     construct projector */

    i__1 = *nb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *nvar;
	for (j = 1; j <= i__2; ++j) {
	    proj[j - 1] += s[i__ - 1] * phi_ref(j, i__);
	}
    }

/*     apply projector and normalize new vector */

    sum = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	proj[i__ - 1] = p0[i__] - proj[i__ - 1];
	sum += proj[i__ - 1] * proj[i__ - 1];
    }
    sum = sqrt(sum);
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	proj[i__ - 1] /= sum;
    }

/*     return original array updated */

    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p0[i__] = proj[i__ - 1];
    }
    return 0;
} /* gsort_ */

#undef phi_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine qonvec  --  block iterative vibration utility  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "qonvec" is a vector utility routine used during sliding */
/*     block iterative matrix diagonalization */


/* Subroutine */ int qonvec_(integer *nvar, integer *np, doublereal *u, 
	doublereal *pk, doublereal *p)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal pku[6];


#define u_ref(a_1,a_2) u[(a_2)*75000 + a_1]



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




/*     operate on vector pk with u-transpose */

    /* Parameter adjustments */
    --p;
    --pk;
    u -= 75001;

    /* Function Body */
    i__1 = *np;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pku[i__ - 1] = 0.;
	i__2 = *nvar;
	for (j = 1; j <= i__2; ++j) {
	    pku[i__ - 1] += u_ref(j, i__) * pk[j];
	}
    }

/*     operate with u on the resultant */

    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = 0.;
	i__2 = *np;
	for (j = 1; j <= i__2; ++j) {
	    p[i__] += u_ref(i__, j) * pku[j - 1];
	}
    }

/*     subtract new product from p */

    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = pk[i__] - p[i__];
    }
    return 0;
} /* qonvec_ */

#undef u_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine project  --  remove known vectors from current  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "project" reads locked vectors from a binary file and projects */
/*     them out of the components of the set of trial eigenvectors */
/*     using the relation Y = X - U * U^T * X */


/* Subroutine */ int project_(integer *nvar, integer *nconv, integer *ivb1, 
	integer *ns, integer *m, doublereal *phi)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    alist al__1;

    /* Builtin functions */
    integer s_rsue(cilist *), do_uio(integer *, char *, ftnlen), e_rsue(void),
	     f_rew(alist *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal u[75000], tmp[150], phik[11250000]	/* was [75000]
	    [150] */;

    /* Fortran I/O blocks */
    static cilist io___153 = { 0, 0, 0, 0, 0 };



#define phi_ref(a_1,a_2) phi[(a_2)*75000 + a_1]
#define phik_ref(a_1,a_2) phik[(a_2)*75000 + a_1 - 75001]



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




/*     zero the temporary storage array */

    /* Parameter adjustments */
    phi -= 75001;

    /* Function Body */
    i__1 = *nvar;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    phik_ref(k, i__ + *m) = 0.;
	}
    }

/*     read and scan over the locked eigenvectors */

    i__1 = *nconv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___153.ciunit = *ivb1;
	s_rsue(&io___153);
	i__2 = *nvar;
	for (k = 1; k <= i__2; ++k) {
	    do_uio(&c__1, (char *)&u[k - 1], (ftnlen)sizeof(doublereal));
	}
	e_rsue();
	i__2 = *ns;
	for (j = 1; j <= i__2; ++j) {
	    tmp[j - 1] = 0.;
	    i__3 = *nvar;
	    for (k = 1; k <= i__3; ++k) {
		tmp[j - 1] += u[k - 1] * phi_ref(k, j + *m);
	    }
	}
	i__2 = *ns;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *nvar;
	    for (k = 1; k <= i__3; ++k) {
		phik_ref(k, j + *m) = phik_ref(k, j + *m) + u[k - 1] * tmp[j 
			- 1];
	    }
	}
    }

/*     project locked vectors out of the current set */

    i__1 = *nvar;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    phi_ref(k, i__ + *m) = phi_ref(k, i__ + *m) - phik_ref(k, i__ + *
		    m);
	}
    }
    if (*nconv > 0) {
	al__1.aerr = 0;
	al__1.aunit = *ivb1;
	f_rew(&al__1);
    }
    return 0;
} /* project_ */

#undef phik_ref
#undef phi_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine konvec  --  evaluate Hessian-vector product  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "konvec" finds a Hessian-vector product via finite-difference */
/*     evaluation of the gradient based on atomic displacements */


/* Subroutine */ int konvec_(integer *nvar, doublereal *xm, doublereal *qe, 
	doublereal *uvec, doublereal *kuvec)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static doublereal e;
    static integer i__, j, k;
    static doublereal eps, sum, grd1[75000]	/* was [3][25000] */, grd2[
	    75000]	/* was [3][25000] */, term, delta[75000];


#define grd1_ref(a_1,a_2) grd1[(a_2)*3 + a_1 - 4]
#define grd2_ref(a_1,a_2) grd2[(a_2)*3 + a_1 - 4]



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




/*     estimate displacement based on total average */

    /* Parameter adjustments */
    --kuvec;
    --uvec;
    --qe;
    --xm;

    /* Function Body */
    sum = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum += uvec[i__] * uvec[i__] / xm[i__];
    }

/*     store the coordinate displacements */

    eps = .001 / sqrt(sum);
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	delta[i__ - 1] = eps * uvec[i__] / sqrt(xm[i__]);
    }

/*     compute the forward displacement */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = (i__ - 1) * 3;
	atoms_1.x[i__ - 1] = (qe[k + 1] + delta[k]) * .52917720859;
	atoms_1.y[i__ - 1] = (qe[k + 2] + delta[k + 1]) * .52917720859;
	atoms_1.z__[i__ - 1] = (qe[k + 3] + delta[k + 2]) * .52917720859;
    }
    gradient_(&e, grd1);

/*     compute the backward displacement */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = (i__ - 1) * 3;
	atoms_1.x[i__ - 1] = (qe[k + 1] - delta[k]) * .52917720859;
	atoms_1.y[i__ - 1] = (qe[k + 2] - delta[k + 1]) * .52917720859;
	atoms_1.z__[i__ - 1] = (qe[k + 3] - delta[k + 2]) * .52917720859;
    }
    gradient_(&e, grd2);

/*     update via finite differences */

    term = .26458860429499997 / (eps * 627.5094688);
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = (i__ - 1) * 3;
	for (j = 1; j <= 3; ++j) {
	    kuvec[k + j] = term * (grd1_ref(j, i__) - grd2_ref(j, i__)) / 
		    sqrt(xm[k + j]);
	}
    }
    return 0;
} /* konvec_ */

#undef grd2_ref
#undef grd1_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine transform  --  diagonalize trial basis vectors  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "transform" diagonalizes the current basis vectors to produce */
/*     trial roots for sliding block iterative matrix diagonalization */


/* Subroutine */ int transform_(integer *nb, doublereal *h__, doublereal *c__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal a[150], b[150];
    static integer i__, j, k;
    static doublereal p[150], w[150], y[150], c1[22500]	/* was [150][150] */, 
	    e1[150], h1[11325], ta[150], tb[150];
    extern /* Subroutine */ int diagq_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);


#define c___ref(a_1,a_2) c__[(a_2)*150 + a_1]
#define h___ref(a_1,a_2) h__[(a_2)*150 + a_1]
#define c1_ref(a_1,a_2) c1[(a_2)*150 + a_1 - 151]



/*     pack the upper triangle of matrix */

    /* Parameter adjustments */
    c__ -= 151;
    h__ -= 151;

    /* Function Body */
    k = 0;
    i__1 = *nb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *nb;
	for (j = i__; j <= i__2; ++j) {
	    ++k;
	    h1[k - 1] = h___ref(i__, j);
	}
    }

/*     perform the matrix diagonalization */

    diagq_(nb, &c__150, nb, h1, e1, c1, a, b, p, w, ta, tb, y);

/*     copy values into the return arrays */

    i__1 = *nb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *nb;
	for (j = 1; j <= i__2; ++j) {
	    h___ref(i__, j) = 0.;
	    c___ref(i__, j) = c1_ref(i__, j);
	}
	h___ref(i__, i__) = e1[i__ - 1];
    }
    return 0;
} /* transform_ */

#undef c1_ref
#undef h___ref
#undef c___ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine diagblk  -- diagonalization for atom block  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "diagblk" performs diagonalization of the Hessian for a */
/*     block of atoms within a larger system */


/* Subroutine */ int diagblk_(integer *k0, integer *k1, integer *n, 
	doublereal *vector, doublereal *wres)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal a[1000], b[1000];
    static integer i__, j, k, m;
    static doublereal p[1000], w[1000], y[1000], ta[1000], tb[1000], hvec[
	    1000000]	/* was [1000][1000] */, hval[1000], hres[500500];
    extern /* Subroutine */ int diagq_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);


#define hvec_ref(a_1,a_2) hvec[(a_2)*1000 + a_1 - 1001]



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




/*     pack the upper triangle of matrix */

    /* Parameter adjustments */
    --wres;
    --vector;

    /* Function Body */
    k = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = *k0 + (i__ - 1) * *n;
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    ++k;
	    hres[k - 1] = vector[m + j];
	}
    }

/*     perform the matrix diagonalization */

    diagq_(n, &c__1000, n, hres, hval, hvec, a, b, p, w, ta, tb, y);

/*     copy values into return arrays */

    k = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    ++k;
	    vector[*k0 + k] = hvec_ref(j, i__);
	}
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wres[*k1 + i__ - 1] = hval[i__ - 1];
    }
    return 0;
} /* diagblk_ */

#undef hvec_ref




/*     ######################################################### */
/*     ##                                                     ## */
/*     ##  subroutine prtvib  --  output of vibrational mode  ## */
/*     ##                                                     ## */
/*     ######################################################### */


/*     "prtvib" writes to an external disk file a series of */
/*     coordinate sets representing motion along a vibrational */
/*     normal mode */


/* Subroutine */ int prtvib_(integer *ivib, doublereal *p)
{
    /* System generated locals */
    address a__1[3];
    integer i__1[3], i__2, i__3;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern integer freeunit_(void);
    static integer i__, j, k;
    static char ext[7];
    static doublereal xref[25000], yref[25000], zref[25000];
    static integer lext, ixyz;
    static doublereal ratio;
    static integer nview;
    extern /* Subroutine */ int prtxyz_(integer *), numeral_(integer *, char *
	    , integer *, ftnlen), version_(char *, char *, ftnlen, ftnlen);
    static char xyzfile[120];



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
/*     ##  files.i  --  name and number of current structure files  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     nprior     number of previously existing cycle files */
/*     ldir       length in characters of the directory name */
/*     leng       length in characters of the base filename */
/*     filename   base filename used by default for all files */
/*     outfile    output filename used for intermediate results */




/*     create a name for the vibrational displacement file */

    /* Parameter adjustments */
    --p;

    /* Function Body */
    lext = 3;
    numeral_(ivib, ext, &lext, (ftnlen)7);
/* Writing concatenation */
    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
    i__1[1] = 1, a__1[1] = ".";
    i__1[2] = lext, a__1[2] = ext;
    s_cat(xyzfile, a__1, i__1, &c__3, (ftnlen)120);
    ixyz = freeunit_();
    version_(xyzfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = ixyz;
    o__1.ofnmlen = 120;
    o__1.ofnm = xyzfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);

/*     store the original atomic coordinates */

    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	xref[i__ - 1] = atoms_1.x[i__ - 1];
	yref[i__ - 1] = atoms_1.y[i__ - 1];
	zref[i__ - 1] = atoms_1.z__[i__ - 1];
    }

/*     make file with plus and minus the current vibration */

    nview = 3;
    i__2 = nview;
    for (i__ = -nview; i__ <= i__2; ++i__) {
	ratio = (doublereal) i__ / (doublereal) nview;
	i__3 = atoms_1.n;
	for (k = 1; k <= i__3; ++k) {
	    j = (k - 1) * 3;
	    atoms_1.x[k - 1] = xref[k - 1] + ratio * p[j + 1];
	    atoms_1.y[k - 1] = yref[k - 1] + ratio * p[j + 2];
	    atoms_1.z__[k - 1] = zref[k - 1] + ratio * p[j + 3];
	}
	prtxyz_(&ixyz);
    }
    cl__1.cerr = 0;
    cl__1.cunit = ixyz;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     restore the original atomic coordinates */

    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	atoms_1.x[i__ - 1] = xref[i__ - 1];
	atoms_1.y[i__ - 1] = yref[i__ - 1];
	atoms_1.z__[i__ - 1] = zref[i__ - 1];
    }
    return 0;
} /* prtvib_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine hessblk  --  Hessian elements for atom block  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "hessblk" calls subroutines to calculate the Hessian elements */
/*     for each atom in turn with respect to Cartesian coordinates */


/* Subroutine */ int hessblk_(doublereal *amass, integer *k0, integer *i1, 
	integer *i2, doublereal *vector)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int estrbnd2_(integer *), erxnfld2_(integer *), 
	    eopdist2_(integer *), eimprop2_(integer *), eimptor2_(integer *), 
	    epitors2_(integer *), etortor2_(integer *), estrtor2_(integer *);
    static integer i__, j, k, ii;
    static doublereal ami, rdn;
    extern /* Subroutine */ int elj2_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal amik;
    extern /* Subroutine */ int born_(void);
    static doublereal xred[25000], yred[25000], zred[25000];
    extern /* Subroutine */ int ehal2_(integer *, doublereal *, doublereal *, 
	    doublereal *), piscf_(void), ebond2_(integer *), ebuck2_(integer *
	    , doublereal *, doublereal *, doublereal *), egeom2_(integer *), 
	    extra2_(integer *), esolv2_(integer *), eurey2_(integer *), 
	    etors2_(integer *), emm3hb2_(integer *, doublereal *, doublereal *
	    , doublereal *), induce_(void);
    static doublereal cutoff;
    extern /* Subroutine */ int bounds_(void), nblist_(void), eangle2_(
	    integer *), emetal2_(integer *), empole2_(integer *), egauss2_(
	    integer *, doublereal *, doublereal *, doublereal *), replica_(
	    doublereal *), chkpole_(void), echarge2_(integer *), eangang2_(
	    integer *), rotpole_(void), echgdpl2_(integer *), eopbend2_(
	    integer *), edipole2_(integer *);


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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  hescut.i  --  cutoff value for Hessian matrix elements  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     hesscut   magnitude of smallest allowed Hessian element */




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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  rigid.i  --  rigid body coordinates for atom groups  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     xrb         rigid body reference x-coordinate for each atom */
/*     yrb         rigid body reference y-coordinate for each atom */
/*     zrb         rigid body reference z-coordinate for each atom */
/*     rbc         current rigid body coordinates for each group */
/*     use_rigid   flag to mark use of rigid body coordinate system */




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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  vdw.i  --  van der Waals parameters for current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     radmin     minimum energy distance for each atom class pair */
/*     epsilon    well depth parameter for each atom class pair */
/*     radmin4    minimum energy distance for 1-4 interaction pairs */
/*     epsilon4   well depth parameter for 1-4 interaction pairs */
/*     radhbnd    minimum energy distance for hydrogen bonding pairs */
/*     epshbnd    well depth parameter for hydrogen bonding pairs */
/*     kred       value of reduction factor parameter for each atom */
/*     ired       attached atom from which reduction factor is applied */
/*     nvdw       total number van der Waals active sites in the system */
/*     ivdw       number of the atom for each van der Waals active site */
/*     jvdw       type or class index into vdw parameters for each atom */
/*     nvt        number of distinct vdw types/classes in the system */
/*     ivt        type/class index for each distinct vdw type or class */
/*     jvt        frequency of each vdw type or class in the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  vdwpot.i  --  specifics of van der Waals functional form  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     abuck       value of "A" constant in Buckingham vdw potential */
/*     bbuck       value of "B" constant in Buckingham vdw potential */
/*     cbuck       value of "C" constant in Buckingham vdw potential */
/*     ghal        value of "gamma" in buffered 14-7 vdw potential */
/*     dhal        value of "delta" in buffered 14-7 vdw potential */
/*     v2scale     factor by which 1-2 vdw interactions are scaled */
/*     v3scale     factor by which 1-3 vdw interactions are scaled */
/*     v4scale     factor by which 1-4 vdw interactions are scaled */
/*     v5scale     factor by which 1-5 vdw interactions are scaled */
/*     igauss      coefficients of Gaussian fit to vdw potential */
/*     ngauss      number of Gaussians used in fit to vdw potential */
/*     use_vcorr   flag to use long range vdw der Waals correction */
/*     vdwindex    indexing mode (atom type or class) for vdw parameters */
/*     vdwtyp      type of van der Waals potential energy function */
/*     radtyp      type of parameter (sigma or R-min) for atomic size */
/*     radsiz      atomic size provided as radius or diameter */
/*     radrule     combining rule for atomic size parameters */
/*     epsrule     combining rule for vdw well depth parameters */
/*     gausstyp    type of Gaussian fit to van der Waals potential */




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




/*     maintain any periodic boundary conditions */

    /* Parameter adjustments */
    --vector;
    --amass;

    /* Function Body */
    if (bound_1.use_bounds__ && ! rigid_1.use_rigid__) {
	bounds_();
    }

/*     update the pairwise interaction neighbor lists */

    if (cutoff_1.use_list__) {
	nblist_();
    }

/*     many implicit solvation models require Born radii */

    if (potent_1.use_born__) {
	born_();
    }

/*     alter bond and torsion constants for pisystem */

    if (potent_1.use_orbit__) {
	piscf_();
    }

/*     compute the induced dipoles at polarizable atoms */

    if (potent_1.use_polar__) {
	chkpole_();
	rotpole_();
	induce_();
    }

/*     calculate the "reduced" atomic coordinates */

    if (potent_1.use_vdw__) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = vdw_1.ired[i__ - 1];
	    rdn = vdw_1.kred[i__ - 1];
	    xred[i__ - 1] = rdn * (atoms_1.x[i__ - 1] - atoms_1.x[ii - 1]) + 
		    atoms_1.x[ii - 1];
	    yred[i__ - 1] = rdn * (atoms_1.y[i__ - 1] - atoms_1.y[ii - 1]) + 
		    atoms_1.y[ii - 1];
	    zred[i__ - 1] = rdn * (atoms_1.z__[i__ - 1] - atoms_1.z__[ii - 1])
		     + atoms_1.z__[ii - 1];
	}
    }

/*     zero out the Hessian elements for the current atom */

    ii = 0;
    i__1 = *i2;
    for (i__ = *i1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    i__2 = *i2;
	    for (k = *i1; k <= i__2; ++k) {
		for (j = 1; j <= 3; ++j) {
		    hessx_ref(j, k) = 0.;
		    hessy_ref(j, k) = 0.;
		    hessz_ref(j, k) = 0.;
		}
	    }

/*     remove any previous use of the replicates method */

	    cutoff = 0.;
	    replica_(&cutoff);

/*     call the local geometry Hessian component routines */

	    if (potent_1.use_bond__) {
		ebond2_(&i__);
	    }
	    if (potent_1.use_angle__) {
		eangle2_(&i__);
	    }
	    if (potent_1.use_strbnd__) {
		estrbnd2_(&i__);
	    }
	    if (potent_1.use_urey__) {
		eurey2_(&i__);
	    }
	    if (potent_1.use_angang__) {
		eangang2_(&i__);
	    }
	    if (potent_1.use_opbend__) {
		eopbend2_(&i__);
	    }
	    if (potent_1.use_opdist__) {
		eopdist2_(&i__);
	    }
	    if (potent_1.use_improp__) {
		eimprop2_(&i__);
	    }
	    if (potent_1.use_imptor__) {
		eimptor2_(&i__);
	    }
	    if (potent_1.use_tors__) {
		etors2_(&i__);
	    }
	    if (potent_1.use_pitors__) {
		epitors2_(&i__);
	    }
	    if (potent_1.use_strtor__) {
		estrtor2_(&i__);
	    }
	    if (potent_1.use_tortor__) {
		etortor2_(&i__);
	    }

/*     call the van der Waals Hessian component routines */

	    if (potent_1.use_vdw__) {
		if (s_cmp(vdwpot_1.vdwtyp, "LENNARD-JONES", (ftnlen)13, (
			ftnlen)13) == 0) {
		    elj2_(&i__, xred, yred, zred);
		} else if (s_cmp(vdwpot_1.vdwtyp, "BUCKINGHAM", (ftnlen)13, (
			ftnlen)10) == 0) {
		    ebuck2_(&i__, xred, yred, zred);
		} else if (s_cmp(vdwpot_1.vdwtyp, "MM3-HBOND", (ftnlen)13, (
			ftnlen)9) == 0) {
		    emm3hb2_(&i__, xred, yred, zred);
		} else if (s_cmp(vdwpot_1.vdwtyp, "BUFFERED-14-7", (ftnlen)13,
			 (ftnlen)13) == 0) {
		    ehal2_(&i__, xred, yred, zred);
		} else if (s_cmp(vdwpot_1.vdwtyp, "GAUSSIAN", (ftnlen)13, (
			ftnlen)8) == 0) {
		    egauss2_(&i__, xred, yred, zred);
		}
	    }

/*     call the electrostatic Hessian component routines */

	    if (potent_1.use_charge__) {
		echarge2_(&i__);
	    }
	    if (potent_1.use_chgdpl__) {
		echgdpl2_(&i__);
	    }
	    if (potent_1.use_dipole__) {
		edipole2_(&i__);
	    }
	    if (potent_1.use_mpole__ || potent_1.use_polar__) {
		empole2_(&i__);
	    }
	    if (potent_1.use_rxnfld__) {
		erxnfld2_(&i__);
	    }

/*     call any miscellaneous Hessian component routines */

	    if (potent_1.use_solv__) {
		esolv2_(&i__);
	    }
	    if (potent_1.use_metal__) {
		emetal2_(&i__);
	    }
	    if (potent_1.use_geom__) {
		egeom2_(&i__);
	    }
	    if (potent_1.use_extra__) {
		extra2_(&i__);
	    }

/*     store Hessian for the current atom block as a vector */

	    ami = .28002851809110429 / (sqrt(amass[i__]) * 627.5094688);
	    i__2 = *i2;
	    for (k = *i1; k <= i__2; ++k) {
		amik = ami / sqrt(amass[k]);
		for (j = 1; j <= 3; ++j) {
		    ++ii;
		    vector[*k0 + ii] = hessx_ref(j, k) * amik;
		}
	    }
	    i__2 = *i2;
	    for (k = *i1; k <= i__2; ++k) {
		amik = ami / sqrt(amass[k]);
		for (j = 1; j <= 3; ++j) {
		    ++ii;
		    vector[*k0 + ii] = hessy_ref(j, k) * amik;
		}
	    }
	    i__2 = *i2;
	    for (k = *i1; k <= i__2; ++k) {
		amik = ami / sqrt(amass[k]);
		for (j = 1; j <= 3; ++j) {
		    ++ii;
		    vector[*k0 + ii] = hessz_ref(j, k) * amik;
		}
	    }
	}
    }
    return 0;
} /* hessblk_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref


/* Main program alias */ int vibbig_ () { MAIN__ (); return 0; }
