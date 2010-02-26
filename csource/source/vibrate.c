/* vibrate.f -- translated by f2c (version 20050501).
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
    doublereal hesscut;
} hescut_;

#define hescut_1 hescut_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__1000 = 1000;
static integer c__1 = 1;
static doublereal c_b12 = 1.;
static integer c__3 = 3;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  program vibrate  --  vibrational analysis and normal modes  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "vibrate" performs a vibrational normal mode analysis; the */
/*     Hessian matrix of second derivatives is determined and then */
/*     diagonalized both directly and after mass weighting; output */
/*     consists of the eigenvalues of the force constant matrix as */
/*     well as the vibrational frequencies and displacements */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 VIBRATE  --  Too many Atoms in the Molec"
	    "ule\002)";
    static char fmt_20[] = "(/,\002 Eigenvalues of the Hessian Matrix :\002,"
	    "/)";
    static char fmt_30[] = "(5(i5,f10.3))";
    static char fmt_40[] = "(/,\002 Vibrational Frequencies (cm-1) :\002,/)";
    static char fmt_50[] = "(5(i5,f10.3))";
    static char fmt_70[] = "(/,\002 Enter Vibrations to Output [List, A=Al"
	    "l\002,\002 or <CR>=Exit] :  \002,$)";
    static char fmt_80[] = "(a120)";
    static char fmt_100[] = "(/,\002 Vibrational Normal Mode\002,i6,\002 wit"
	    "h Frequency\002,f11.3,\002 cm-1\002,//,5x,\002Atom\002,5x,\002De"
	    "lta X\002,5x,\002Delta Y\002,5x,\002Delta Z\002,/)";
    static char fmt_110[] = "(4x,i5,3f12.6)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2, i__3[3], i__4;
    doublereal d__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    double sqrt(doublereal);
    integer do_fio(integer *, char *, ftnlen);
    double d_sign(doublereal *, doublereal *);
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_rsfe(cilist *), e_rsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    extern integer freeunit_(void);
    static doublereal a[1000], b[1000], h__[1000000];
    static integer i__, j, k, m;
    static doublereal p[1000], w[1000], ta[1000], tb[1000];
    static integer iv[1000];
    static doublereal ty[1000];
    static char ext[7];
    static integer ivib, nvib;
    static doublereal xref[25000], yref[25000], zref[25000];
    static integer list[1000], lext, next, ixyz;
    static doublereal mass2[25000], hdiag[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int diagq_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal eigen[1000];
    extern /* Subroutine */ int fatal_(void), final_(void);
    static integer ihess, nfreq, hinit[75000]	/* was [3][25000] */;
    static doublereal ratio;
    static integer ilist;
    static doublereal vects[1000000]	/* was [1000][1000] */;
    static integer nview, nlist;
    static logical exist;
    static integer hstop[75000]	/* was [3][25000] */;
    static doublereal vnorm;
    static logical query;
    static doublereal factor;
    static integer hindex[1000000];
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char letter[1];
    static doublereal matrix[500500];
    static char string[120];
    static integer ndummy;
    extern /* Subroutine */ int getxyz_(void), prtxyz_(integer *), initial_(
	    void), hessian_(doublereal *, integer *, integer *, integer *, 
	    doublereal *), numeral_(integer *, char *, integer *, ftnlen), 
	    nextarg_(char *, logical *, ftnlen), gettext_(char *, char *, 
	    integer *, ftnlen, ftnlen), version_(char *, char *, ftnlen, 
	    ftnlen);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___38 = { 1, string, 1, 0, 120, 1 };
    static cilist io___39 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_80, 0 };
    static icilist io___44 = { 1, record, 1, 0, 120, 1 };
    static cilist io___47 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___48 = { 0, 0, 0, fmt_110, 0 };



#define hdiag_ref(a_1,a_2) hdiag[(a_2)*3 + a_1 - 4]
#define hinit_ref(a_1,a_2) hinit[(a_2)*3 + a_1 - 4]
#define vects_ref(a_1,a_2) vects[(a_2)*1000 + a_1 - 1001]
#define hstop_ref(a_1,a_2) hstop[(a_2)*3 + a_1 - 4]



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




/*     set up the structure and mechanics calculation */

    initial_();
    getxyz_();
    mechanic_();

/*     initialize various things needed for vibrations */

    nfreq = usage_1.nuse * 3;
    ndummy = 0;
    if (nfreq > 1000 || nfreq * nfreq > 1000000) {
	io___3.ciunit = iounit_1.iout;
	s_wsfe(&io___3);
	e_wsfe();
	fatal_();
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1] && atmtyp_1.atomic[i__ - 1] == 0) {
	    ++ndummy;
	    atmtyp_1.mass[i__ - 1] = .001;
	}
	mass2[i__ - 1] = sqrt(atmtyp_1.mass[i__ - 1]);
    }
    nvib = nfreq - ndummy * 3;

/*     calculate the Hessian matrix of second derivatives */

    hescut_1.hesscut = 0.;
    hessian_(h__, hinit, hstop, hindex, hdiag);

/*     store upper triangle of the Hessian in "matrix" */

    ihess = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    for (j = 1; j <= 3; ++j) {
		++ihess;
		matrix[ihess - 1] = hdiag_ref(j, i__);
		i__2 = hstop_ref(j, i__);
		for (k = hinit_ref(j, i__); k <= i__2; ++k) {
		    m = (hindex[k - 1] + 2) / 3;
		    if (usage_1.use[m - 1]) {
			++ihess;
			matrix[ihess - 1] = h__[k - 1];
		    }
		}
	    }
	}
    }

/*     perform diagonalization to get Hessian eigenvalues */

    diagq_(&nfreq, &c__1000, &nfreq, matrix, eigen, vects, a, b, p, w, ta, tb,
	     ty);
    io___26.ciunit = iounit_1.iout;
    s_wsfe(&io___26);
    e_wsfe();
    io___27.ciunit = iounit_1.iout;
    s_wsfe(&io___27);
    i__1 = nvib;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&eigen[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();

/*     store upper triangle of the mass-weighted Hessian matrix */

    ihess = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    for (j = 1; j <= 3; ++j) {
		++ihess;
		matrix[ihess - 1] = hdiag_ref(j, i__) / atmtyp_1.mass[i__ - 1]
			;
		i__2 = hstop_ref(j, i__);
		for (k = hinit_ref(j, i__); k <= i__2; ++k) {
		    m = (hindex[k - 1] + 2) / 3;
		    if (usage_1.use[m - 1]) {
			++ihess;
			matrix[ihess - 1] = h__[k - 1] / (mass2[i__ - 1] * 
				mass2[m - 1]);
		    }
		}
	    }
	}
    }

/*     diagonalize to get vibrational frequencies and normal modes */

    diagq_(&nfreq, &c__1000, &nfreq, matrix, eigen, vects, a, b, p, w, ta, tb,
	     ty);
    factor = sqrt(418.4) / .1883651567308853;
    i__1 = nvib;
    for (i__ = 1; i__ <= i__1; ++i__) {
	eigen[i__ - 1] = factor * d_sign(&c_b12, &eigen[i__ - 1]) * sqrt((
		d__1 = eigen[i__ - 1], abs(d__1)));
    }
    io___29.ciunit = iounit_1.iout;
    s_wsfe(&io___29);
    e_wsfe();
    io___30.ciunit = iounit_1.iout;
    s_wsfe(&io___30);
    i__1 = nvib;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&eigen[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();

/*     form Cartesian coordinate displacements from normal modes */

    i__1 = nvib;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vnorm = 0.;
	i__2 = nfreq;
	for (j = 1; j <= i__2; ++j) {
	    k = usage_1.iuse[(j + 2) / 3 - 1];
	    vects_ref(j, i__) = vects_ref(j, i__) / mass2[k - 1];
/* Computing 2nd power */
	    d__1 = vects_ref(j, i__);
	    vnorm += d__1 * d__1;
	}
	vnorm = sqrt(vnorm);
	i__2 = nfreq;
	for (j = 1; j <= i__2; ++j) {
	    vects_ref(j, i__) = vects_ref(j, i__) / vnorm;
	}
    }

/*     try to get output vibrational modes from command line */

    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	query = FALSE_;
	*(unsigned char *)letter = *(unsigned char *)string;
	upcase_(letter, (ftnlen)1);
	if (*(unsigned char *)letter == 'A') {
	    nlist = nvib;
	    i__1 = nlist;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		list[i__ - 1] = i__;
	    }
	} else {
	    nlist = 0;
	    i__1 = nvib;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = s_rsli(&io___38);
		if (i__2 != 0) {
		    goto L60;
		}
		i__2 = do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(
			integer));
		if (i__2 != 0) {
		    goto L60;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L60;
		}
		if (k >= 1 && k <= nvib) {
		    ++nlist;
		    list[nlist - 1] = k;
		}
		nextarg_(string, &exist, (ftnlen)120);
	    }
L60:
	    ;
	}
    }

/*     ask the user for the vibrational modes to be output */

    if (query) {
	io___39.ciunit = iounit_1.iout;
	s_wsfe(&io___39);
	e_wsfe();
	io___40.ciunit = iounit_1.input;
	s_rsfe(&io___40);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	*(unsigned char *)letter = ' ';
	next = 1;
	gettext_(record, letter, &next, (ftnlen)120, (ftnlen)1);
	upcase_(letter, (ftnlen)1);
	if (*(unsigned char *)letter == ' ') {
	    nlist = 0;
	} else if (*(unsigned char *)letter == 'A') {
	    nlist = nvib;
	    i__1 = nlist;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		list[i__ - 1] = i__;
	    }
	} else {
	    i__1 = nvib;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		iv[i__ - 1] = 0;
	    }
	    i__1 = s_rsli(&io___44);
	    if (i__1 != 0) {
		goto L90;
	    }
	    i__2 = nvib;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__1 = do_lio(&c__3, &c__1, (char *)&iv[i__ - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L90;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L90;
	    }
L90:
	    nlist = 0;
	    i__1 = nvib;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		k = iv[i__ - 1];
		if (k >= 1 && k <= nvib) {
		    ++nlist;
		    list[nlist - 1] = k;
		}
	    }
	}
    }

/*     print the vibrational frequencies and normal modes */

    i__1 = nlist;
    for (ilist = 1; ilist <= i__1; ++ilist) {
	ivib = list[ilist - 1];
	io___47.ciunit = iounit_1.iout;
	s_wsfe(&io___47);
	do_fio(&c__1, (char *)&ivib, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&eigen[ivib - 1], (ftnlen)sizeof(doublereal));
	e_wsfe();
	i__2 = usage_1.nuse;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    j = (i__ - 1) * 3;
	    io___48.ciunit = iounit_1.iout;
	    s_wsfe(&io___48);
	    do_fio(&c__1, (char *)&usage_1.iuse[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&vects_ref(j + 1, ivib), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&vects_ref(j + 2, ivib), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&vects_ref(j + 3, ivib), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}

/*     create a name for the vibrational displacement file */

	lext = 3;
	numeral_(&ivib, ext, &lext, (ftnlen)7);
/* Writing concatenation */
	i__3[0] = files_1.leng, a__1[0] = files_1.filename;
	i__3[1] = 1, a__1[1] = ".";
	i__3[2] = lext, a__1[2] = ext;
	s_cat(xyzfile, a__1, i__3, &c__3, (ftnlen)120);
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
	    i__4 = usage_1.nuse;
	    for (k = 1; k <= i__4; ++k) {
		j = (k - 1) * 3;
		m = usage_1.iuse[k - 1];
		atoms_1.x[m - 1] = xref[m - 1] + ratio * vects_ref(j + 1, 
			ivib);
		atoms_1.y[m - 1] = yref[m - 1] + ratio * vects_ref(j + 2, 
			ivib);
		atoms_1.z__[m - 1] = zref[m - 1] + ratio * vects_ref(j + 3, 
			ivib);
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
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef hstop_ref
#undef vects_ref
#undef hinit_ref
#undef hdiag_ref


/* Main program alias */ int vibrate_ () { MAIN__ (); return 0; }
