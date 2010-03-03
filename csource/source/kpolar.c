/* kpolar.f -- translated by f2c (version 20050501).
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
    doublereal polr[5000], athl[5000];
    integer pgrp[40000]	/* was [8][5000] */;
} kpolr_;

#define kpolr_1 kpolr_

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
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    integer np11[25000], ip11[2500000]	/* was [100][25000] */, np12[25000], 
	    ip12[1250000]	/* was [50][25000] */, np13[25000], ip13[
	    1250000]	/* was [50][25000] */, np14[25000], ip14[1250000]	
	    /* was [50][25000] */;
} polgrp_;

#define polgrp_1 polgrp_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__3 = 3;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine kpolar  --  assign polarizability parameters  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "kpolar" assigns atomic dipole polarizabilities to the atoms */
/*     within the structure and processes any new or changed values */


/* Subroutine */ int kpolar_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Additional Atomic Dipole\002,\002 Polari"
	    "zability Parameters :\002,//,5x,\002Atom Type\002,11x,\002Alph"
	    "a\002,8x,\002Damp\002,5x,\002Group Atom Types\002/)";
    static char fmt_40[] = "(4x,i6,10x,f10.3,2x,f10.3,7x,20i5)";
    static char fmt_50[] = "(/,\002 KPOLAR  --  Too many Dipole\002,\002 Pol"
	    "arizability Parameters\002)";
    static char fmt_70[] = "(/,\002 Additional Dipole Polarizabilities\002"
	    ",\002 for Specific Atoms :\002,//,6x,\002Atom\002,15x,\002Alph"
	    "a\002,8x,\002Damp\002,/)";
    static char fmt_80[] = "(4x,i6,10x,f10.3,2x,f10.3)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern /* Subroutine */ int polargrp_(void);
    static integer i__, j, k, pg[8], npg;
    static doublereal thl, pol;
    static integer next;
    static doublereal sixth;
    static logical header;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int chkpole_(void), getnumb_(char *, integer *, 
	    integer *, ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___12 = { 1, string, 1, 0, 120, 1 };
    static cilist io___13 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___17 = { 1, string, 1, 0, 120, 1 };
    static cilist io___18 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_80, 0 };



#define pole_ref(a_1,a_2) mpole_1.pole[(a_2)*13 + a_1 - 14]
#define pgrp_ref(a_1,a_2) kpolr_1.pgrp[(a_2)*8 + a_1 - 9]
#define polaxe_ref(a_0,a_1) &mpole_1.polaxe[(a_1)*8 + a_0 - 8]
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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  kpolr.i  --  forcefield parameters for polarizability  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     polr   dipole polarizability parameters for each atom type */
/*     athl   Thole polarizability damping value for each atom type */
/*     pgrp   connected types in polarization group of each atom type */




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




/*     process keywords containing polarizability parameters */

    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "POLARIZE ", (ftnlen)9, (ftnlen)9) == 0) {
	    k = 0;
	    pol = 0.;
	    thl = -1.;
	    for (j = 1; j <= 8; ++j) {
		pg[j - 1] = 0;
	    }
	    getnumb_(record, &k, &next, (ftnlen)120);
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___12);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&pol, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&thl, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    for (j = 1; j <= 8; ++j) {
		i__2 = do_lio(&c__3, &c__1, (char *)&pg[j - 1], (ftnlen)
			sizeof(integer));
		if (i__2 != 0) {
		    goto L10;
		}
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
L10:
	    if (k > 0) {
		if (header) {
		    header = FALSE_;
		    io___13.ciunit = iounit_1.iout;
		    s_wsfe(&io___13);
		    e_wsfe();
		}
		if (k <= 5000) {
		    kpolr_1.polr[k - 1] = pol;
		    kpolr_1.athl[k - 1] = thl;
		    for (j = 1; j <= 8; ++j) {
			pgrp_ref(j, k) = pg[j - 1];
			if (pg[j - 1] == 0) {
			    npg = j - 1;
			    goto L30;
			}
		    }
L30:
		    io___15.ciunit = iounit_1.iout;
		    s_wsfe(&io___15);
		    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&pol, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&thl, (ftnlen)sizeof(doublereal));
		    i__2 = npg;
		    for (j = 1; j <= i__2; ++j) {
			do_fio(&c__1, (char *)&pg[j - 1], (ftnlen)sizeof(
				integer));
		    }
		    e_wsfe();
		} else {
		    io___16.ciunit = iounit_1.iout;
		    s_wsfe(&io___16);
		    e_wsfe();
		    inform_1.abort = TRUE_;
		}
	    }
	}
    }

/*     find and store the atomic dipole polarizability parameters */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	polar_1.polarity[i__ - 1] = kpolr_1.polr[atoms_1.type__[i__ - 1] - 1];
	polar_1.thole[i__ - 1] = kpolr_1.athl[atoms_1.type__[i__ - 1] - 1];
    }

/*     process keywords containing atom specific polarizabilities */

    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "POLARIZE ", (ftnlen)9, (ftnlen)9) == 0) {
	    k = 0;
	    pol = 0.;
	    thl = 0.;
	    getnumb_(record, &k, &next, (ftnlen)120);
	    if (k < 0 && k >= -atoms_1.n) {
		k = -k;
		s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next 
			- 1));
		i__2 = s_rsli(&io___17);
		if (i__2 != 0) {
		    goto L60;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&pol, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L60;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&thl, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L60;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L60;
		}
L60:
		if (header) {
		    header = FALSE_;
		    io___18.ciunit = iounit_1.iout;
		    s_wsfe(&io___18);
		    e_wsfe();
		}
		io___19.ciunit = iounit_1.iout;
		s_wsfe(&io___19);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&pol, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&thl, (ftnlen)sizeof(doublereal));
		e_wsfe();
		polar_1.polarity[k - 1] = pol;
		polar_1.thole[k - 1] = thl;
	    }
	}
    }

/*     remove zero and undefined polarizable sites from the list */

    mpole_1.npole = 0;
    polar_1.npolar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (mpole_1.polsiz[i__ - 1] != 0 || polar_1.polarity[i__ - 1] != 0.) {
	    ++mpole_1.npole;
	    mpole_1.ipole[mpole_1.npole - 1] = i__;
	    mpole_1.pollist[i__ - 1] = mpole_1.npole;
	    mpole_1.zaxis[mpole_1.npole - 1] = mpole_1.zaxis[i__ - 1];
	    mpole_1.xaxis[mpole_1.npole - 1] = mpole_1.xaxis[i__ - 1];
	    mpole_1.yaxis[mpole_1.npole - 1] = mpole_1.yaxis[i__ - 1];
	    s_copy(polaxe_ref(0, mpole_1.npole), polaxe_ref(0, i__), (ftnlen)
		    8, (ftnlen)8);
	    for (k = 1; k <= 13; ++k) {
		pole_ref(k, mpole_1.npole) = pole_ref(k, i__);
	    }
	    if (polar_1.polarity[i__ - 1] != 0.) {
		++polar_1.npolar;
	    }
	    polar_1.polarity[mpole_1.npole - 1] = polar_1.polarity[i__ - 1];
	    polar_1.thole[mpole_1.npole - 1] = polar_1.thole[i__ - 1];
	}
    }

/*     set the values used in the scaling of the polarizability */

    sixth = .16666666666666666;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (polar_1.thole[i__ - 1] == 0.) {
	    polar_1.pdamp[i__ - 1] = 0.;
	} else {
	    polar_1.pdamp[i__ - 1] = pow_dd(&polar_1.polarity[i__ - 1], &
		    sixth);
	}
    }

/*     assign polarization group connectivity of each atom */

    polargrp_();

/*     test multipoles at chiral sites and invert if necessary */

    chkpole_();

/*     turn off polarizable multipole potential if it is not used */

    if (mpole_1.npole == 0) {
	potent_1.use_mpole__ = FALSE_;
    }
    if (polar_1.npolar == 0) {
	potent_1.use_polar__ = FALSE_;
    }
    return 0;
} /* kpolar_ */

#undef keyline_ref
#undef polaxe_ref
#undef pgrp_ref
#undef pole_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine polargrp  --  polarization group connectivity  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "polargrp" generates members of the polarization group of */
/*     each atom and separate lists of the 1-2, 1-3 and 1-4 group */
/*     connectivities */


/* Subroutine */ int polargrp_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 POLARGRP  --  Too many Atoms\002,\002 in"
	    " Polarization Group\002)";
    static char fmt_30[] = "(/,\002 POLARGRP  --  Too many Atoms\002,\002 in"
	    " Polarization Group\002)";
    static char fmt_40[] = "(/,\002 POLARGRP  --  Too many Atoms\002,\002 in"
	    " 1-2 Polarization Groups\002)";
    static char fmt_50[] = "(/,\002 POLARGRP  --  Too many Atoms\002,\002 in"
	    " 1-3 Polarization Groups\002)";
    static char fmt_60[] = "(/,\002 POLARGRP  --  Too many Atoms\002,\002 in"
	    " 1-4 Polarization Groups\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k, jj, kk, it, jt, keep[25000];
    static logical done;
    static integer mask[25000], list[25000], stop;
    extern /* Subroutine */ int sort_(integer *, integer *), sort8_(integer *,
	     integer *);
    static integer nkeep, nlist, start;

    /* Fortran I/O blocks */
    static cilist io___28 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_60, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define ip12_ref(a_1,a_2) polgrp_1.ip12[(a_2)*50 + a_1 - 51]
#define ip13_ref(a_1,a_2) polgrp_1.ip13[(a_2)*50 + a_1 - 51]
#define ip14_ref(a_1,a_2) polgrp_1.ip14[(a_2)*50 + a_1 - 51]
#define pgrp_ref(a_1,a_2) kpolr_1.pgrp[(a_2)*8 + a_1 - 9]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  kpolr.i  --  forcefield parameters for polarizability  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     polr   dipole polarizability parameters for each atom type */
/*     athl   Thole polarizability damping value for each atom type */
/*     pgrp   connected types in polarization group of each atom type */




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




/*     find the directly connected group members for each atom */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	polgrp_1.np11[i__ - 1] = 1;
	ip11_ref(1, i__) = i__;
	it = atoms_1.type__[i__ - 1];
	i__2 = couple_1.n12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    jj = i12_ref(j, i__);
	    jt = atoms_1.type__[jj - 1];
	    for (k = 1; k <= 8; ++k) {
		kk = pgrp_ref(k, it);
		if (kk == 0) {
		    goto L20;
		}
		if (pgrp_ref(k, it) == jt) {
		    ++polgrp_1.np11[i__ - 1];
		    if (polgrp_1.np11[i__ - 1] <= 100) {
			ip11_ref(polgrp_1.np11[i__ - 1], i__) = jj;
		    } else {
			io___28.ciunit = iounit_1.iout;
			s_wsfe(&io___28);
			e_wsfe();
			inform_1.abort = TRUE_;
		    }
		}
	    }
L20:
	    ;
	}
    }

/*     find any other group members for each atom in turn */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	list[i__ - 1] = 0;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	start = 1;
	stop = polgrp_1.np11[i__ - 1];
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    jj = ip11_ref(j, i__);
	    if (jj < i__) {
		done = TRUE_;
		polgrp_1.np11[i__ - 1] = polgrp_1.np11[jj - 1];
		i__3 = polgrp_1.np11[i__ - 1];
		for (k = 1; k <= i__3; ++k) {
		    ip11_ref(k, i__) = ip11_ref(k, jj);
		}
	    } else {
		list[jj - 1] = i__;
	    }
	}
	while(! done) {
	    done = TRUE_;
	    i__2 = stop;
	    for (j = start; j <= i__2; ++j) {
		jj = ip11_ref(j, i__);
		i__3 = polgrp_1.np11[jj - 1];
		for (k = 1; k <= i__3; ++k) {
		    kk = ip11_ref(k, jj);
		    if (list[kk - 1] != i__) {
			++polgrp_1.np11[i__ - 1];
			if (polgrp_1.np11[i__ - 1] <= 100) {
			    ip11_ref(polgrp_1.np11[i__ - 1], i__) = kk;
			} else {
			    io___33.ciunit = iounit_1.iout;
			    s_wsfe(&io___33);
			    e_wsfe();
			    inform_1.abort = TRUE_;
			}
			list[kk - 1] = i__;
		    }
		}
	    }
	    if (polgrp_1.np11[i__ - 1] != stop) {
		done = FALSE_;
		start = stop + 1;
		stop = polgrp_1.np11[i__ - 1];
	    }
	}
	sort_(&polgrp_1.np11[i__ - 1], &ip11_ref(1, i__));
    }

/*     loop over atoms finding all the 1-2 group relationships */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mask[i__ - 1] = 0;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nlist = 0;
	i__2 = polgrp_1.np11[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    jj = ip11_ref(j, i__);
	    ++nlist;
	    list[nlist - 1] = jj;
	    mask[jj - 1] = i__;
	}
	nkeep = 0;
	i__2 = nlist;
	for (j = 1; j <= i__2; ++j) {
	    jj = list[j - 1];
	    i__3 = couple_1.n12[jj - 1];
	    for (k = 1; k <= i__3; ++k) {
		kk = i12_ref(k, jj);
		if (mask[kk - 1] != i__) {
		    ++nkeep;
		    keep[nkeep - 1] = kk;
		}
	    }
	}
	nlist = 0;
	i__2 = nkeep;
	for (j = 1; j <= i__2; ++j) {
	    jj = keep[j - 1];
	    i__3 = polgrp_1.np11[jj - 1];
	    for (k = 1; k <= i__3; ++k) {
		kk = ip11_ref(k, jj);
		++nlist;
		list[nlist - 1] = kk;
	    }
	}
	sort8_(&nlist, list);
	if (nlist <= 50) {
	    polgrp_1.np12[i__ - 1] = nlist;
	    i__2 = nlist;
	    for (j = 1; j <= i__2; ++j) {
		ip12_ref(j, i__) = list[j - 1];
	    }
	} else {
	    io___38.ciunit = iounit_1.iout;
	    s_wsfe(&io___38);
	    e_wsfe();
	    inform_1.abort = TRUE_;
	}
    }

/*     loop over atoms finding all the 1-3 group relationships */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mask[i__ - 1] = 0;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = polgrp_1.np11[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    jj = ip11_ref(j, i__);
	    mask[jj - 1] = i__;
	}
	i__2 = polgrp_1.np12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    jj = ip12_ref(j, i__);
	    mask[jj - 1] = i__;
	}
	nlist = 0;
	i__2 = polgrp_1.np12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    jj = ip12_ref(j, i__);
	    i__3 = polgrp_1.np12[jj - 1];
	    for (k = 1; k <= i__3; ++k) {
		kk = ip12_ref(k, jj);
		if (mask[kk - 1] != i__) {
		    ++nlist;
		    list[nlist - 1] = kk;
		}
	    }
	}
	sort8_(&nlist, list);
	if (nlist <= 50) {
	    polgrp_1.np13[i__ - 1] = nlist;
	    i__2 = nlist;
	    for (j = 1; j <= i__2; ++j) {
		ip13_ref(j, i__) = list[j - 1];
	    }
	} else {
	    io___39.ciunit = iounit_1.iout;
	    s_wsfe(&io___39);
	    e_wsfe();
	    inform_1.abort = TRUE_;
	}
    }

/*     loop over atoms finding all the 1-4 group relationships */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mask[i__ - 1] = 0;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = polgrp_1.np11[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    jj = ip11_ref(j, i__);
	    mask[jj - 1] = i__;
	}
	i__2 = polgrp_1.np12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    jj = ip12_ref(j, i__);
	    mask[jj - 1] = i__;
	}
	i__2 = polgrp_1.np13[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    jj = ip13_ref(j, i__);
	    mask[jj - 1] = i__;
	}
	nlist = 0;
	i__2 = polgrp_1.np13[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    jj = ip13_ref(j, i__);
	    i__3 = polgrp_1.np12[jj - 1];
	    for (k = 1; k <= i__3; ++k) {
		kk = ip12_ref(k, jj);
		if (mask[kk - 1] != i__) {
		    ++nlist;
		    list[nlist - 1] = kk;
		}
	    }
	}
	sort8_(&nlist, list);
	if (nlist <= 50) {
	    polgrp_1.np14[i__ - 1] = nlist;
	    i__2 = nlist;
	    for (j = 1; j <= i__2; ++j) {
		ip14_ref(j, i__) = list[j - 1];
	    }
	} else {
	    io___40.ciunit = iounit_1.iout;
	    s_wsfe(&io___40);
	    e_wsfe();
	    inform_1.abort = TRUE_;
	}
    }
    return 0;
} /* polargrp_ */

#undef pgrp_ref
#undef ip14_ref
#undef ip13_ref
#undef ip12_ref
#undef ip11_ref
#undef i12_ref


