/* monte.f -- translated by f2c (version 20050501).
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
    doublereal dihed[1000];
    integer nomega, iomega[2000]	/* was [2][1000] */, zline[1000];
} omega_;

#define omega_1 omega_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal zbond[25000], zang[25000], ztors[25000];
    integer iz[100000]	/* was [4][25000] */;
} zcoord_;

#define zcoord_1 zcoord_

struct {
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__5 = 5;
static integer c__2 = 2;



/*     ################################################################ */
/*     ##                   COPYRIGHT (C) 2001 by                    ## */
/*     ##  Michael Schnieders, Alan Grossfield & Jay William Ponder  ## */
/*     ##                    All Rights Reserved                     ## */
/*     ################################################################ */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  program monte  --  Monte Carlo-Minimization search method  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "monte" performs a Monte Carlo-Minimization conformational */
/*     search using Cartesian single atom or torsional move sets */

/*     literature references: */

/*     Z. Li and H. A. Scheraga, "Monte Carlo-Minimization Approach */
/*     to the Multiple-Minima Problem in Protein Folding", Proc. Natl. */
/*     Acad. Sci. USA, 84, 6611-6615 (1987) */

/*     D. J. Wales, "Energy Landscapes with Applications to Clusters, */
/*     Biomolecules and Glasses", Cambridge University Press, 2003, */
/*     Section 6.7.4 */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Number of Monte Carlo Steps [1000] : "
	    "\002,$)";
    static char fmt_30[] = "(i15)";
    static char fmt_40[] = "(/,\002 Use [C]artesian or [T]orsional Moves [C]"
	    " : \002,$)";
    static char fmt_50[] = "(a120)";
    static char fmt_70[] = "(/,\002 Enter Maximum Step in Degrees [180] :"
	    " \002,$)";
    static char fmt_80[] = "(/,\002 Enter Maximum Step in Angstroms [3.0] :"
	    " \002,$)";
    static char fmt_90[] = "(a120)";
    static char fmt_120[] = "(/,\002 Enter the Desired Temperature in Degr"
	    "ees\002,\002 K [500] : \002,$)";
    static char fmt_130[] = "(a120)";
    static char fmt_160[] = "(/,\002 Enter RMS Gradient Criterion [0.01] :"
	    " \002,$)";
    static char fmt_170[] = "(a120)";
    static char fmt_190[] = "(/,\002 Monte Carlo Minimization Global Searc"
	    "h :\002)";
    static char fmt_200[] = "(/,\002 MCM Iter       Current         Global  "
	    "     Temper\002,\002      Ratio      Status\002,/)";
    static char fmt_210[] = "(i8,3x,f12.4)";
    static char fmt_220[] = "(i8,3x,f12.4,3x,f12.4,3x,f10.2,3x,f8.3,6x,a6)";
    static char fmt_230[] = "(i8,9x,\002------\002,3x,f12.4,3x,f10.2,3x,f8.3"
	    ",6x,a6)";
    static char fmt_240[] = "(/,\002 Global Minimum Energy Value :\002,1x,f2"
	    "0.8)";
    static char fmt_250[] = "(/,\002 Global Minimum Energy Value :\002,3x,f1"
	    "8.6)";
    static char fmt_260[] = "(/,\002 Global Minimum Energy Value :\002,5x,f1"
	    "6.4)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2], i__3;
    doublereal d__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);
    double log(doublereal);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double exp(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    extern integer freeunit_(void);
    static doublereal pminimum;
    static logical torsmove;
    static integer i__, k, m;
    static doublereal xg[25000], yg[25000], zg[25000], xi[25000], yi[25000], 
	    zi[25000], xp[25000], yp[25000], zp[25000], big, eps, beta;
    static integer nbig, keep, imin;
    static doublereal size;
    static integer next;
    extern /* Subroutine */ int final_(void);
    static doublereal trial, ratio;
    static logical reset;
    static integer istep, nstep;
    static doublereal boltz;
    static logical exist;
    static doublereal tsize, global, factor, grdmin;
    extern doublereal random_(void);
    static char record[120];
    static doublereal temper, vector[3];
    static char answer[1], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), ranvec_(doublereal *)
	    ;
    static char status[6];
    extern /* Subroutine */ int getxyz_(void), prtxyz_(integer *);
    static char minfile[120];
    extern /* Subroutine */ int makeint_(integer *), initial_(void), nextarg_(
	    char *, logical *, ftnlen);
    static doublereal minimum;
    extern /* Subroutine */ int mcmstep_(doublereal *, doublereal *), 
	    gettext_(char *, char *, integer *, ftnlen, ftnlen), version_(
	    char *, char *, ftnlen, ftnlen), initrot_(void), makexyz_(void);

    /* Fortran I/O blocks */
    static icilist io___9 = { 1, string, 1, 0, 120, 1 };
    static cilist io___10 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___20 = { 1, string, 1, 0, 120, 1 };
    static cilist io___21 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_90, 0 };
    static icilist io___24 = { 1, string, 1, 0, 120, 1 };
    static icilist io___26 = { 1, string, 1, 0, 120, 1 };
    static cilist io___27 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_130, 0 };
    static icilist io___29 = { 1, string, 1, 0, 120, 1 };
    static icilist io___32 = { 1, string, 1, 0, 120, 1 };
    static cilist io___33 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_170, 0 };
    static icilist io___35 = { 1, string, 1, 0, 120, 1 };
    static cilist io___36 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___64 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___65 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___66 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_260, 0 };




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
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     zbond   bond length used to define each Z-matrix atom */
/*     zang    bond angle used to define each Z-matrix atom */
/*     ztors   angle or torsion used to define Z-matrix atom */
/*     iz      defining atom numbers for each Z-matrix atom */




/*     set up the structure and mechanics calculation */

    initial_();
    getxyz_();
    mechanic_();

/*     initialize values of some counters and parameters */

    keep = 0;
    nbig = 0;
    eps = 1e-4;
    big = 1e5;
    reset = FALSE_;

/*     get the desired number of Monte Carlo steps */

    nstep = -1;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___9);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&nstep, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    }
L10:
    if (nstep <= 0) {
	io___10.ciunit = iounit_1.iout;
	s_wsfe(&io___10);
	e_wsfe();
	io___11.ciunit = iounit_1.input;
	s_rsfe(&io___11);
	do_fio(&c__1, (char *)&nstep, (ftnlen)sizeof(integer));
	e_rsfe();
	if (nstep <= 0) {
	    nstep = 1000;
	}
    }

/*     choose either the torsional or single atom move set */

    torsmove = FALSE_;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___14.ciunit = iounit_1.iout;
	s_wsfe(&io___14);
	e_wsfe();
	io___15.ciunit = iounit_1.input;
	s_rsfe(&io___15);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'T') {
	torsmove = TRUE_;
    }

/*     generate the internal coordinates, keep all atoms active */

    if (torsmove) {
	makeint_(&c__0);
	initrot_();
	usage_1.nuse = atoms_1.n;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    usage_1.use[i__ - 1] = TRUE_;
	}
    }

/*     get the desired Cartesian or torsional step size */

    size = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___20);
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&size, (ftnlen)sizeof(doublereal))
		;
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L60;
	}
    }
L60:
    if (size < 0.) {
	if (torsmove) {
	    io___21.ciunit = iounit_1.iout;
	    s_wsfe(&io___21);
	    e_wsfe();
	} else {
	    io___22.ciunit = iounit_1.iout;
	    s_wsfe(&io___22);
	    e_wsfe();
	}
	io___23.ciunit = iounit_1.input;
	s_rsfe(&io___23);
	do_fio(&c__1, string, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___24);
	if (i__1 != 0) {
	    goto L100;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&size, (ftnlen)sizeof(doublereal))
		;
	if (i__1 != 0) {
	    goto L100;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L100;
	}
L100:
	if (size < 0.) {
	    if (torsmove) {
		size = 180.;
	    } else {
		size = 3.;
	    }
	}
	if (torsmove) {
	    size = min(size,180.);
	}
    }

/*     get the desired temperature for Metropolis criterion */

    temper = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___26);
	if (i__1 != 0) {
	    goto L110;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&temper, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L110;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L110;
	}
    }
L110:
    if (temper < 0.) {
	io___27.ciunit = iounit_1.iout;
	s_wsfe(&io___27);
	e_wsfe();
	io___28.ciunit = iounit_1.input;
	s_rsfe(&io___28);
	do_fio(&c__1, string, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___29);
	if (i__1 != 0) {
	    goto L140;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&temper, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L140;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L140;
	}
L140:
	if (temper < 0.) {
	    temper = 500.;
	}
    }
    beta = 1. / (temper * .0019872066);

/*     get the gradient convergence for local minimizations */

    grdmin = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___32);
	if (i__1 != 0) {
	    goto L110;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&grdmin, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L110;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L110;
	}
    }
/* L150: */
    if (grdmin < 0.) {
	io___33.ciunit = iounit_1.iout;
	s_wsfe(&io___33);
	e_wsfe();
	io___34.ciunit = iounit_1.input;
	s_rsfe(&io___34);
	do_fio(&c__1, string, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___35);
	if (i__1 != 0) {
	    goto L180;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&grdmin, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L180;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L180;
	}
L180:
	if (grdmin < 0.) {
	    grdmin = .01f;
	}
    }

/*     print some information prior to initial iteration */

    io___36.ciunit = iounit_1.iout;
    s_wsfe(&io___36);
    e_wsfe();
    io___37.ciunit = iounit_1.iout;
    s_wsfe(&io___37);
    e_wsfe();

/*     store the coordinates, then perform a minimization */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi[i__ - 1] = atoms_1.x[i__ - 1];
	yi[i__ - 1] = atoms_1.y[i__ - 1];
	zi[i__ - 1] = atoms_1.z__[i__ - 1];
    }
    mcmstep_(&minimum, &grdmin);
    pminimum = minimum;
    io___43.ciunit = iounit_1.iout;
    s_wsfe(&io___43);
    do_fio(&c__1, (char *)&c__0, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     save coordinates as the initial global minimum */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xg[i__ - 1] = atoms_1.x[i__ - 1];
	yg[i__ - 1] = atoms_1.y[i__ - 1];
	zg[i__ - 1] = atoms_1.z__[i__ - 1];
    }
    global = minimum;
    imin = freeunit_();
/* Writing concatenation */
    i__2[0] = files_1.leng, a__1[0] = files_1.filename;
    i__2[1] = 4, a__1[1] = ".xyz";
    s_cat(minfile, a__1, i__2, &c__2, (ftnlen)120);
    version_(minfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = imin;
    o__1.ofnmlen = 120;
    o__1.ofnm = minfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    prtxyz_(&imin);
    cl__1.cerr = 0;
    cl__1.cunit = imin;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     optionally reset coordinates to before the minimization */

    if (reset) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    atoms_1.x[i__ - 1] = xi[i__ - 1];
	    atoms_1.y[i__ - 1] = yi[i__ - 1];
	    atoms_1.z__[i__ - 1] = zi[i__ - 1];
	}
    }
    if (torsmove) {
	makeint_(&c__2);
    }

/*     store the prior coordinates to start each MCM iteration */

    i__1 = nstep;
    for (istep = 1; istep <= i__1; ++istep) {
	i__3 = atoms_1.n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    xp[i__ - 1] = atoms_1.x[i__ - 1];
	    yp[i__ - 1] = atoms_1.y[i__ - 1];
	    zp[i__ - 1] = atoms_1.z__[i__ - 1];
	}

/*     generate random angle moves for a few torsions */

	if (torsmove) {
/* Computing MAX */
	    d__1 = random_();
	    m = (integer) (-log((max(d__1,1e-4)))) + 1;
	    i__3 = m;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		k = (integer) (omega_1.nomega * random_()) + 1;
		k = omega_1.zline[k - 1];
		tsize = size * 2. * (random_() - .5);
		zcoord_1.ztors[k - 1] += tsize;
		if (zcoord_1.ztors[k - 1] > 180.) {
		    zcoord_1.ztors[k - 1] += -360.;
		} else if (zcoord_1.ztors[k - 1] < -180.) {
		    zcoord_1.ztors[k - 1] += 360.;
		}
	    }
	    makexyz_();

/*     generate a random Cartesian move for each atom */

	} else {
	    i__3 = atoms_1.n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		if (usage_1.use[i__ - 1]) {
		    ranvec_(vector);
		    factor = size * random_();
		    atoms_1.x[i__ - 1] += factor * vector[0];
		    atoms_1.y[i__ - 1] += factor * vector[1];
		    atoms_1.z__[i__ - 1] += factor * vector[2];
		}
	    }
	}

/*     store the coordinates, then perform a minimization */

	i__3 = atoms_1.n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    xi[i__ - 1] = atoms_1.x[i__ - 1];
	    yi[i__ - 1] = atoms_1.y[i__ - 1];
	    zi[i__ - 1] = atoms_1.z__[i__ - 1];
	}
	mcmstep_(&minimum, &grdmin);

/*     test for an unreasonably low energy at the minimum */

	if (minimum < -big) {
	    minimum = big;
	}

/*     step is probably degenerate if energy is identical */

	if ((d__1 = minimum - pminimum, abs(d__1)) <= eps) {
	    s_copy(status, "Same", (ftnlen)6, (ftnlen)4);
	    pminimum = minimum;

/*     accept the step if the new minimum has lower energy */

	} else if (minimum <= pminimum) {
	    s_copy(status, "Accept", (ftnlen)6, (ftnlen)6);
	    pminimum = minimum;

/*     if the energy increased, apply the Metropolis criterion */

	} else {
	    boltz = exp(-beta * (minimum - pminimum));
	    trial = random_();

/*     reject the step if the energy increase is too large */

	    if (boltz < trial) {
		s_copy(status, "Reject", (ftnlen)6, (ftnlen)6);

/*     accept the step if the energy increase is small enough */

	    } else {
		s_copy(status, "Accept", (ftnlen)6, (ftnlen)6);
		pminimum = minimum;
	    }
	}

/*     save coordinates with the best energy as global minimum */

	if (minimum < global) {
	    i__3 = atoms_1.n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		xg[i__ - 1] = atoms_1.x[i__ - 1];
		yg[i__ - 1] = atoms_1.y[i__ - 1];
		zg[i__ - 1] = atoms_1.z__[i__ - 1];
	    }
	    global = minimum;
	    imin = freeunit_();
/* Writing concatenation */
	    i__2[0] = files_1.leng, a__1[0] = files_1.filename;
	    i__2[1] = 4, a__1[1] = ".xyz";
	    s_cat(minfile, a__1, i__2, &c__2, (ftnlen)120);
	    version_(minfile, "old", (ftnlen)120, (ftnlen)3);
	    o__1.oerr = 0;
	    o__1.ounit = imin;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = minfile;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	    prtxyz_(&imin);
	    cl__1.cerr = 0;
	    cl__1.cunit = imin;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	}

/*     update the overall Monte Carlo acceptance ratio */

	if (s_cmp(status, "Accept", (ftnlen)6, (ftnlen)6) == 0) {
	    ++keep;
	}
	ratio = (doublereal) keep / (doublereal) istep;

/*     print intermediate results for the current iteration */

	if (minimum < big) {
	    nbig = 0;
	    io___63.ciunit = iounit_1.iout;
	    s_wsfe(&io___63);
	    do_fio(&c__1, (char *)&istep, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&global, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&temper, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ratio, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, status, (ftnlen)6);
	    e_wsfe();
	} else {
	    ++nbig;
	    io___64.ciunit = iounit_1.iout;
	    s_wsfe(&io___64);
	    do_fio(&c__1, (char *)&istep, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&global, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&temper, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ratio, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, status, (ftnlen)6);
	    e_wsfe();
	}

/*     restore global minimum after repeated bad iterations */

	if (nbig >= 3) {
	    nbig = 0;
	    i__3 = atoms_1.n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		atoms_1.x[i__ - 1] = xg[i__ - 1];
		atoms_1.y[i__ - 1] = yg[i__ - 1];
		atoms_1.z__[i__ - 1] = zg[i__ - 1];
	    }

/*     optionally reset coordinates to before the minimization */

	} else if (s_cmp(status, "Same", (ftnlen)6, (ftnlen)4) == 0 || s_cmp(
		status, "Accept", (ftnlen)6, (ftnlen)6) == 0) {
	    if (reset) {
		i__3 = atoms_1.n;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    atoms_1.x[i__ - 1] = xi[i__ - 1];
		    atoms_1.y[i__ - 1] = yi[i__ - 1];
		    atoms_1.z__[i__ - 1] = zi[i__ - 1];
		}
	    }

/*     restore coordinates to those from the previous iteration */

	} else if (s_cmp(status, "Reject", (ftnlen)6, (ftnlen)6) == 0) {
	    i__3 = atoms_1.n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		atoms_1.x[i__ - 1] = xp[i__ - 1];
		atoms_1.y[i__ - 1] = yp[i__ - 1];
		atoms_1.z__[i__ - 1] = zp[i__ - 1];
	    }
	}

/*     update internal coordinates if using torsional moves */

	if (torsmove) {
	    makeint_(&c__2);
	}
    }

/*     write out the final global minimum energy value */

    if (inform_1.digits >= 8) {
	io___65.ciunit = iounit_1.iout;
	s_wsfe(&io___65);
	do_fio(&c__1, (char *)&global, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else if (inform_1.digits >= 6) {
	io___66.ciunit = iounit_1.iout;
	s_wsfe(&io___66);
	do_fio(&c__1, (char *)&global, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else {
	io___67.ciunit = iounit_1.iout;
	s_wsfe(&io___67);
	do_fio(&c__1, (char *)&global, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  function mcmstep  --  minimization phase of an MCM step  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "mcmstep" implements the minimization phase of an MCM step */
/*     via Cartesian minimization following a Monte Carlo step */


/* Subroutine */ int mcmstep_(doublereal *minimum, doublereal *grdmin)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal xx[75000];
    extern doublereal mcm1_();
    extern /* Subroutine */ int mcm2_();
    static char mode[6];
    extern /* Subroutine */ int tncg_(char *, char *, integer *, doublereal *,
	     doublereal *, doublereal *, D_fp, U_fp, U_fp, ftnlen, ftnlen);
    static integer nvar;
    static char method[6];
    extern /* Subroutine */ int optsave_();



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  output.i  --  control of coordinate output file format  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     archive    logical flag to save structures in an archive */
/*     noversion  logical flag governing use of filename versions */
/*     overwrite  logical flag to overwrite intermediate files inplace */
/*     cyclesave  logical flag to mark use of numbered cycle files */
/*     coordtype  selects Cartesian, internal, rigid body or none */




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




/*     prepare for the truncated Newton minimization */

    s_copy(mode, "AUTO", (ftnlen)6, (ftnlen)4);
    s_copy(method, "AUTO", (ftnlen)6, (ftnlen)4);
    inform_1.verbose = FALSE_;
    inform_1.iprint = 0;
    inform_1.iwrite = 0;
    s_copy(output_1.coordtype, "CARTESIAN", (ftnlen)9, (ftnlen)9);

/*     translate the coordinates of each active atom */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    xx[nvar - 1] = atoms_1.x[i__ - 1];
	    ++nvar;
	    xx[nvar - 1] = atoms_1.y[i__ - 1];
	    ++nvar;
	    xx[nvar - 1] = atoms_1.z__[i__ - 1];
	}
    }

/*     make the call to the optimization routine */

    tncg_(mode, method, &nvar, xx, minimum, grdmin, (D_fp)mcm1_, (U_fp)mcm2_, 
	    (U_fp)optsave_, (ftnlen)6, (ftnlen)6);

/*     untranslate the final coordinates for active atoms */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    atoms_1.x[i__ - 1] = xx[nvar - 1];
	    ++nvar;
	    atoms_1.y[i__ - 1] = xx[nvar - 1];
	    ++nvar;
	    atoms_1.z__[i__ - 1] = xx[nvar - 1];
	}
    }
    return 0;
} /* mcmstep_ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  function mcm1  --  energy and gradient for MCM search  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "mcm1" is a service routine that computes the energy */
/*     and gradient for truncated Newton optimization in Cartesian */
/*     coordinate space */


doublereal mcm1_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static doublereal e;
    static integer i__, nvar;
    static doublereal derivs[75000]	/* was [3][25000] */;


#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1 - 4]



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




/*     translate optimization parameters to atomic coordinates */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    atoms_1.x[i__ - 1] = xx[nvar];
	    ++nvar;
	    atoms_1.y[i__ - 1] = xx[nvar];
	    ++nvar;
	    atoms_1.z__[i__ - 1] = xx[nvar];
	}
    }

/*     compute and store the energy and gradient */

    gradient_(&e, derivs);
    ret_val = e;

/*     store Cartesian gradient as optimization gradient */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    g[nvar] = derivs_ref(1, i__);
	    ++nvar;
	    g[nvar] = derivs_ref(2, i__);
	    ++nvar;
	    g[nvar] = derivs_ref(3, i__);
	}
    }
    return ret_val;
} /* mcm1_ */

#undef derivs_ref




/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  subroutine mcm2  --  Hessian values for MCM search  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     "mcm2" is a service routine that computes the sparse */
/*     matrix Hessian elements for truncated Newton optimization */
/*     in Cartesian coordinate space */


/* Subroutine */ int mcm2_(char *mode, doublereal *xx, doublereal *h__, 
	integer *hinit, integer *hstop, integer *hindex, doublereal *hdiag, 
	ftnlen mode_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k, hvar[75000], huse[75000], nvar;
    extern /* Subroutine */ int hessian_(doublereal *, integer *, integer *, 
	    integer *, doublereal *);



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




/*     translate optimization parameters to atomic coordinates */

    /* Parameter adjustments */
    --hdiag;
    --hindex;
    --hstop;
    --hinit;
    --h__;
    --xx;

    /* Function Body */
    if (s_cmp(mode, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	return 0;
    }
    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    atoms_1.x[i__ - 1] = xx[nvar];
	    ++nvar;
	    atoms_1.y[i__ - 1] = xx[nvar];
	    ++nvar;
	    atoms_1.z__[i__ - 1] = xx[nvar];
	}
    }

/*     compute and store the Hessian elements */

    hessian_(&h__[1], &hinit[1], &hstop[1], &hindex[1], &hdiag[1]);

/*     transform the sparse Hessian to use only active atoms */

    nvar = 0;
    if (usage_1.nuse != atoms_1.n) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = (i__ - 1) * 3;
	    if (usage_1.use[i__ - 1]) {
		for (j = 1; j <= 3; ++j) {
		    ++nvar;
		    hvar[nvar - 1] = j + k;
		    huse[j + k - 1] = nvar;
		}
	    } else {
		for (j = 1; j <= 3; ++j) {
		    huse[j + k - 1] = 0;
		}
	    }
	}
	i__1 = nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = hvar[i__ - 1];
	    hinit[i__] = hinit[k];
	    hstop[i__] = hstop[k];
	    hdiag[i__] = hdiag[k];
	    i__2 = hstop[i__];
	    for (j = hinit[i__]; j <= i__2; ++j) {
		hindex[j] = huse[hindex[j] - 1];
	    }
	}
    }
    return 0;
} /* mcm2_ */

/* Main program alias */ int monte_ () { MAIN__ (); return 0; }
