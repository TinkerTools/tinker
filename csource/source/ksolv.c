/* ksolv.f -- translated by f2c (version 20050501).
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
    doublereal gkr[5000], gkc;
} gk_;

#define gk_1 gk_

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
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal ak[75000], anat[75000], afld[75000];
    integer nangle, iang[300000]	/* was [4][75000] */;
} angle_;

#define angle_1 angle_

struct {
    integer bndlist[200000]	/* was [8][25000] */, anglist[700000]	/* 
	    was [28][25000] */;
} atmlst_;

#define atmlst_1 atmlst_

struct {
    doublereal bk[50000], bl[50000];
    integer nbond, ibnd[100000]	/* was [2][50000] */;
} bond_;

#define bond_1 bond_

struct {
    doublereal electric, dielec, ebuffer, c2scale, c3scale, c4scale, c5scale;
    logical neutnbr, neutcut;
} chgpot_;

#define chgpot_1 chgpot_

struct {
    doublereal rad[5000], eps[5000], rad4[5000], eps4[5000], reduct[5000];
} kvdws_;

#define kvdws_1 kvdws_

struct {
    doublereal kelvin0, kelvin, atmsph, tautemp, taupres, compress, collide, 
	    xnh[2], vnh[2], qnh[2], gnh[2], volmove;
    integer voltrial;
    logical isothermal, isobaric, anisotrop;
    char thermostat[11], barostat[10], volscale[9];
} bath_;

#define bath_1 bath_

struct {
    doublereal solvprs, surften, spcut, spoff, stcut, stoff, rcav[25000], 
	    rdisp[25000], cdisp[25000];
} npolar_;

#define npolar_1 npolar_

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

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5000 = 5000;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine ksolv  --  solvation parameter assignment  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "ksolv" assigns continuum solvation energy parameters for */
/*     the surface area, generalized Born, generalized Kirkwood */
/*     and Poisson-Boltzmann solvation models */


/* Subroutine */ int ksolv_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Additional Continuum Solvation\002,\002 "
	    "Parameters:\002,//,5x,\002Atom Type\002,10x,\002Radius\002,/)";
    static char fmt_30[] = "(4x,i6,8x,f12.4)";
    static char fmt_40[] = "(/,\002 KSOLV  --  Only Atom Types Through\002,i"
	    "4,\002 are Allowed\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, k;
    static doublereal rd;
    extern /* Subroutine */ int kgb_(void), kpb_(void), kgk_(void), ksa_(void)
	    , knp_(void);
    static integer next;
    static char value[20];
    static logical header;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int getword_(char *, char *, integer *, ftnlen, 
	    ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___8 = { 1, string, 1, 0, 120, 1 };
    static icilist io___11 = { 1, string, 1, 0, 120, 1 };
    static cilist io___12 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_40, 0 };



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




/*     defaults for continuum solvation term and parameters */

    potent_1.use_solv__ = FALSE_;
    potent_1.use_born__ = FALSE_;
    s_copy(solute_1.solvtyp, "        ", (ftnlen)8, (ftnlen)8);
    s_copy(solute_1.borntyp, "        ", (ftnlen)8, (ftnlen)8);
    solute_1.doffset = .09;

/*     search keywords for continuum solvation commands */

    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "SOLVATE ", (ftnlen)8, (ftnlen)8) == 0) {
	    potent_1.use_solv__ = TRUE_;
	    potent_1.use_born__ = FALSE_;
	    s_copy(solute_1.solvtyp, "ASP", (ftnlen)8, (ftnlen)3);
	    getword_(record, value, &next, (ftnlen)120, (ftnlen)20);
	    upcase_(value, (ftnlen)20);
	    if (s_cmp(value, "ASP", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(solute_1.solvtyp, "ASP", (ftnlen)8, (ftnlen)3);
	    } else if (s_cmp(value, "SASA", (ftnlen)4, (ftnlen)4) == 0) {
		s_copy(solute_1.solvtyp, "SASA", (ftnlen)8, (ftnlen)4);
	    } else if (s_cmp(value, "GBSA", (ftnlen)4, (ftnlen)4) == 0) {
		potent_1.use_born__ = TRUE_;
		s_copy(solute_1.solvtyp, "STILL", (ftnlen)8, (ftnlen)5);
	    } else if (s_cmp(value, "ONION", (ftnlen)5, (ftnlen)5) == 0) {
		potent_1.use_born__ = TRUE_;
		s_copy(solute_1.solvtyp, "ONION", (ftnlen)8, (ftnlen)5);
	    } else if (s_cmp(value, "STILL", (ftnlen)5, (ftnlen)5) == 0) {
		potent_1.use_born__ = TRUE_;
		s_copy(solute_1.solvtyp, "STILL", (ftnlen)8, (ftnlen)5);
	    } else if (s_cmp(value, "HCT", (ftnlen)3, (ftnlen)3) == 0) {
		potent_1.use_born__ = TRUE_;
		s_copy(solute_1.solvtyp, "HCT", (ftnlen)8, (ftnlen)3);
	    } else if (s_cmp(value, "OBC", (ftnlen)3, (ftnlen)3) == 0) {
		potent_1.use_born__ = TRUE_;
		s_copy(solute_1.solvtyp, "OBC", (ftnlen)8, (ftnlen)3);
	    } else if (s_cmp(value, "ACE", (ftnlen)3, (ftnlen)3) == 0) {
		potent_1.use_born__ = TRUE_;
		s_copy(solute_1.solvtyp, "ACE", (ftnlen)8, (ftnlen)3);
	    } else if (s_cmp(value, "GK", (ftnlen)2, (ftnlen)2) == 0) {
		potent_1.use_born__ = TRUE_;
		s_copy(solute_1.solvtyp, "GK", (ftnlen)8, (ftnlen)2);
	    } else if (s_cmp(value, "PB", (ftnlen)2, (ftnlen)2) == 0) {
		s_copy(solute_1.solvtyp, "PB", (ftnlen)8, (ftnlen)2);
	    }
	} else if (s_cmp(keyword, "BORN-RADIUS ", (ftnlen)12, (ftnlen)12) == 
		0) {
	    getword_(record, value, &next, (ftnlen)120, (ftnlen)20);
	    upcase_(value, (ftnlen)20);
	    if (s_cmp(value, "ONION", (ftnlen)5, (ftnlen)5) == 0) {
		s_copy(solute_1.borntyp, "ONION", (ftnlen)8, (ftnlen)5);
	    } else if (s_cmp(value, "STILL", (ftnlen)5, (ftnlen)5) == 0) {
		s_copy(solute_1.borntyp, "STILL", (ftnlen)8, (ftnlen)5);
	    } else if (s_cmp(value, "HCT", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(solute_1.borntyp, "HCT", (ftnlen)8, (ftnlen)3);
	    } else if (s_cmp(value, "OBC", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(solute_1.borntyp, "OBC", (ftnlen)8, (ftnlen)3);
	    } else if (s_cmp(value, "ACE", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(solute_1.borntyp, "ACE", (ftnlen)8, (ftnlen)3);
	    } else if (s_cmp(value, "GRYCUK", (ftnlen)4, (ftnlen)6) == 0) {
		s_copy(solute_1.borntyp, "GRYCUK", (ftnlen)8, (ftnlen)6);
	    } else if (s_cmp(value, "PERFECT", (ftnlen)7, (ftnlen)7) == 0) {
		s_copy(solute_1.borntyp, "PERFECT", (ftnlen)8, (ftnlen)7);
	    }
	} else if (s_cmp(keyword, "DIELECTRIC-OFFSET ", (ftnlen)18, (ftnlen)
		18) == 0) {
	    i__2 = s_rsli(&io___8);
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&solute_1.doffset, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L50;
	    }
	    if (solute_1.doffset < 0.) {
		solute_1.doffset = -solute_1.doffset;
	    }
	} else if (s_cmp(keyword, "GKR ", (ftnlen)4, (ftnlen)4) == 0) {
	    k = 0;
	    rd = 0.;
	    i__2 = s_rsli(&io___11);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&rd, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
L10:
	    if (k >= 1 && k <= 1000) {
		if (header) {
		    header = FALSE_;
		    io___12.ciunit = iounit_1.iout;
		    s_wsfe(&io___12);
		    e_wsfe();
		}
		if (rd < 0.) {
		    rd = 0.;
		}
		gk_1.gkr[k - 1] = rd;
		io___13.ciunit = iounit_1.iout;
		s_wsfe(&io___13);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&rd, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else if (k > 5000) {
		io___14.ciunit = iounit_1.iout;
		s_wsfe(&io___14);
		do_fio(&c__1, (char *)&c__5000, (ftnlen)sizeof(integer));
		e_wsfe();
		inform_1.abort = TRUE_;
	    }
	}
L50:
	;
    }

/*     set a default if no Born radius method was assigned */

    if (potent_1.use_born__ && s_cmp(solute_1.borntyp, "       ", (ftnlen)8, (
	    ftnlen)7) == 0) {
	s_copy(solute_1.borntyp, solute_1.solvtyp, (ftnlen)8, (ftnlen)8);
	if (s_cmp(solute_1.solvtyp, "GK", (ftnlen)8, (ftnlen)2) == 0) {
	    s_copy(solute_1.borntyp, "GRYCUK", (ftnlen)8, (ftnlen)6);
	}
    }

/*     invoke the setup needed for specific Born radius models */

    if (s_cmp(solute_1.borntyp, "PERFECT", (ftnlen)8, (ftnlen)7) == 0) {
	kpb_();
    }

/*     invoke the setup needed for specific solvation models */

    if (s_cmp(solute_1.solvtyp, "ASP", (ftnlen)8, (ftnlen)3) == 0 || s_cmp(
	    solute_1.solvtyp, "SASA", (ftnlen)8, (ftnlen)4) == 0) {
	ksa_();
    } else if (s_cmp(solute_1.solvtyp, "GK", (ftnlen)8, (ftnlen)2) == 0) {
	kgk_();
	knp_();
    } else if (s_cmp(solute_1.solvtyp, "PB", (ftnlen)8, (ftnlen)2) == 0) {
	kpb_();
	knp_();
    } else {
	kgb_();
    }
    return 0;
} /* ksolv_ */

#undef keyline_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine ksa  --  set surface area solvation parameters  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "ksa" initializes parameters needed for surface area only */
/*     solvation models including ASP and SASA */

/*     literature references: */

/*     L. Wesson and D. Eisenberg, "Atomic Solvation Parameters */
/*     Applied to Molecular Dynamics of Proteins in Solution", */
/*     Protein Science, 1, 227-235 (1992)  (Eisenberg-McLachlan ASP) */

/*     T. Ooi, M. Oobatake, G. Nemethy and H. A. Scheraga, "Accessible */
/*     Surface Areas as a Measure of the Thermodynamic Parameters of */
/*     Hydration of Peptides", PNAS, 84, 3086-3090 (1987)  (SASA) */


/* Subroutine */ int ksa_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k, atmnum;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]



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




/*     assign the Eisenberg-McLachlan ASP solvation parameters; */
/*     parameters only available for protein-peptide groups */

    if (s_cmp(solute_1.solvtyp, "ASP", (ftnlen)8, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    if (atmnum == 6) {
		solute_1.rsolv[i__ - 1] = 1.9;
		solute_1.asolv[i__ - 1] = .004;
	    } else if (atmnum == 7) {
		solute_1.rsolv[i__ - 1] = 1.7;
		solute_1.asolv[i__ - 1] = -.113;
		if (couple_1.n12[i__ - 1] == 4) {
		    solute_1.asolv[i__ - 1] = -.169;
		}
	    } else if (atmnum == 8) {
		solute_1.rsolv[i__ - 1] = 1.4;
		solute_1.asolv[i__ - 1] = -.113;
		if (couple_1.n12[i__ - 1] == 1 && atmtyp_1.atomic[i12_ref(1, 
			i__) - 1] == 6) {
		    i__2 = couple_1.n13[i__ - 1];
		    for (j = 1; j <= i__2; ++j) {
			k = i13_ref(j, i__);
			if (couple_1.n12[k - 1] == 1 && atmtyp_1.atomic[k - 1]
				 == 8) {
			    solute_1.asolv[i__ - 1] = -.166;
			}
		    }
		}
		i__2 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, i__);
		    if (atmtyp_1.atomic[k - 1] == 15) {
			solute_1.asolv[i__ - 1] = -.14;
		    }
		}
	    } else if (atmnum == 15) {
		solute_1.rsolv[i__ - 1] = 1.9;
		solute_1.asolv[i__ - 1] = -.14;
	    } else if (atmnum == 16) {
		solute_1.rsolv[i__ - 1] = 1.8;
		solute_1.asolv[i__ - 1] = -.017;
	    } else {
		solute_1.rsolv[i__ - 1] = 0.;
		solute_1.asolv[i__ - 1] = 0.;
	    }
	}
    }

/*     assign the Ooi-Scheraga SASA solvation parameters; */
/*     parameters only available for protein-peptide groups */

    if (s_cmp(solute_1.solvtyp, "SASA", (ftnlen)8, (ftnlen)4) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    if (atmnum == 6) {
		solute_1.rsolv[i__ - 1] = 2.;
		solute_1.asolv[i__ - 1] = .008;
		if (couple_1.n12[i__ - 1] == 3) {
		    solute_1.rsolv[i__ - 1] = 1.75;
		    solute_1.asolv[i__ - 1] = -.008;
		    i__2 = couple_1.n12[i__ - 1];
		    for (j = 1; j <= i__2; ++j) {
			k = i12_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 8) {
			    solute_1.rsolv[i__ - 1] = 1.55;
			    solute_1.asolv[i__ - 1] = .427;
			}
		    }
		}
	    } else if (atmnum == 7) {
		solute_1.rsolv[i__ - 1] = 1.55;
		solute_1.asolv[i__ - 1] = -.132;
		if (couple_1.n12[i__ - 1] == 4) {
		    solute_1.asolv[i__ - 1] = -1.212;
		}
	    } else if (atmnum == 8) {
		solute_1.rsolv[i__ - 1] = 1.4;
		if (couple_1.n12[i__ - 1] == 1) {
		    solute_1.asolv[i__ - 1] = -.038;
		    if (atmtyp_1.atomic[i12_ref(1, i__) - 1] == 6) {
			i__2 = couple_1.n13[i__ - 1];
			for (j = 1; j <= i__2; ++j) {
			    k = i13_ref(j, i__);
			    if (couple_1.n12[k - 1] == 1 && atmtyp_1.atomic[k 
				    - 1] == 8) {
				solute_1.asolv[i__ - 1] = -.77;
			    }
			}
		    }
		} else if (couple_1.n12[i__ - 1] == 2) {
		    solute_1.asolv[i__ - 1] = -.172;
		}
		i__2 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, i__);
		    if (atmtyp_1.atomic[k - 1] == 15) {
			solute_1.asolv[i__ - 1] = -.717;
		    }
		}
	    } else if (atmnum == 15) {
		solute_1.rsolv[i__ - 1] = 2.1;
		solute_1.asolv[i__ - 1] = 0.;
	    } else if (atmnum == 16) {
		solute_1.rsolv[i__ - 1] = 2.;
		solute_1.asolv[i__ - 1] = -.021;
	    } else if (atmnum == 17) {
		solute_1.rsolv[i__ - 1] = 2.;
		solute_1.asolv[i__ - 1] = .012;
	    } else {
		solute_1.rsolv[i__ - 1] = 0.;
		solute_1.asolv[i__ - 1] = 0.;
	    }
	}
    }
    return 0;
} /* ksa_ */

#undef i13_ref
#undef i12_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine kgb  --  assign generalized Born parameters  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "kgb" initializes parameters needed for the generalized */
/*     Born solvation models */

/*     literature references: */

/*     M. Schaefer, C. Bartels, F. Leclerc and M. Karplus, "Effective */
/*     Atom Volumes for Implicit Solvent Models: Comparison between */
/*     Voronoi Volumes and Minimum Fluctuations Volumes", Journal of */
/*     Computational Chemistry, 22, 1857-1879 (2001)  (ACE) */


/* Subroutine */ int kgb_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double cos(doublereal);
    integer i_dnnt(doublereal *);
    double sqrt(doublereal), atan(doublereal);

    /* Local variables */
    static doublereal h__;
    static integer i__, j, k, m;
    static doublereal r__, c1, c2, c3, r2, r4;
    static integer ia, ib, nh, kc, mm, ic, id;
    static doublereal ri, rk, vk, pi2, ri2, rk2, rab, rbc, fik, qik, uik, 
	    s2ik, tik2, s3ik, temp, term, prod2, prod4;
    static logical amide;
    static doublereal alpha;
    extern /* Subroutine */ int kbond_(void);
    static doublereal omgik, ratio, width, qterm, alpha2, alpha4;
    extern /* Subroutine */ int kangle_(void);
    static doublereal factor;
    static integer atmmas;
    static doublereal cosine;
    static integer atmnum;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define uace_ref(a_1,a_2) solute_1.uace[(a_2)*1000 + a_1 - 1001]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define wace_ref(a_1,a_2) solute_1.wace[(a_2)*1000 + a_1 - 1001]
#define s2ace_ref(a_1,a_2) solute_1.s2ace[(a_2)*1000 + a_1 - 1001]
#define bndlist_ref(a_1,a_2) atmlst_1.bndlist[(a_2)*8 + a_1 - 9]



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
/*     ##  angle.i  --  bond angles within the current structure  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     ak       harmonic angle force constant (kcal/mole/rad**2) */
/*     anat     ideal bond angle or phase shift angle (degrees) */
/*     afld     periodicity for Fourier bond angle term */
/*     nangle   total number of bond angles in the system */
/*     iang     numbers of the atoms in each bond angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  atmlst.i  --  local geometry terms involving each atom  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     bndlist   list of the bond numbers involving each atom */
/*     anglist   list of the angle numbers centered on each atom */




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




/*     set offset and scaling values for analytical Still method */

    if (s_cmp(solute_1.solvtyp, "STILL", (ftnlen)8, (ftnlen)5) == 0) {
	solute_1.p1 = .073;
	solute_1.p2 = .921;
	solute_1.p3 = 6.211;
	solute_1.p4 = 15.236;
	solute_1.p5 = 1.254;
	if (! potent_1.use_bond__) {
	    kbond_();
	}
	if (! potent_1.use_angle__) {
	    kangle_();
	}
    }

/*     set overlap scale factors for HCT and OBC methods */

    if (s_cmp(solute_1.solvtyp, "HCT", (ftnlen)8, (ftnlen)3) == 0 || s_cmp(
	    solute_1.solvtyp, "OBC", (ftnlen)8, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.shct[i__ - 1] = .8;
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    if (atmnum == 1) {
		solute_1.shct[i__ - 1] = .85;
	    }
	    if (atmnum == 6) {
		solute_1.shct[i__ - 1] = .72;
	    }
	    if (atmnum == 7) {
		solute_1.shct[i__ - 1] = .79;
	    }
	    if (atmnum == 8) {
		solute_1.shct[i__ - 1] = .85;
	    }
	    if (atmnum == 9) {
		solute_1.shct[i__ - 1] = .88;
	    }
	    if (atmnum == 15) {
		solute_1.shct[i__ - 1] = .86;
	    }
	    if (atmnum == 16) {
		solute_1.shct[i__ - 1] = .96;
	    }
	    if (atmnum == 26) {
		solute_1.shct[i__ - 1] = .88;
	    }
	}
    }

/*     set rescaling coefficients for the OBC method */

    if (s_cmp(solute_1.solvtyp, "OBC", (ftnlen)8, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.aobc[i__ - 1] = 1.;
	    solute_1.bobc[i__ - 1] = .8;
	    solute_1.gobc[i__ - 1] = 4.85;
	}
    }

/*     set the Gaussian width factor for the ACE method */

    if (s_cmp(solute_1.solvtyp, "ACE", (ftnlen)8, (ftnlen)3) == 0) {
	width = 1.2;
    }

/*     assign surface area factors for nonpolar solvation */

    if (s_cmp(solute_1.solvtyp, "ONION", (ftnlen)8, (ftnlen)5) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.asolv[i__ - 1] = .0072;
	}
    } else if (s_cmp(solute_1.solvtyp, "STILL", (ftnlen)8, (ftnlen)5) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.asolv[i__ - 1] = .0049;
	}
    } else if (s_cmp(solute_1.solvtyp, "HCT", (ftnlen)8, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.asolv[i__ - 1] = .0054;
	}
    } else if (s_cmp(solute_1.solvtyp, "OBC", (ftnlen)8, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.asolv[i__ - 1] = .0054;
	}
    } else if (s_cmp(solute_1.solvtyp, "ACE", (ftnlen)8, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.asolv[i__ - 1] = .003;
	}
    }

/*     assign standard radii for GB/SA methods other than ACE; */
/*     taken from Macromodel and OPLS-AA, except for hydrogens */

    if (s_cmp(solute_1.solvtyp, "ACE", (ftnlen)8, (ftnlen)3) != 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    if (atmnum == 1) {
		solute_1.rsolv[i__ - 1] = 1.25;
		k = i12_ref(1, i__);
		if (atmtyp_1.atomic[k - 1] == 7) {
		    solute_1.rsolv[i__ - 1] = 1.15;
		}
		if (atmtyp_1.atomic[k - 1] == 8) {
		    solute_1.rsolv[i__ - 1] = 1.05;
		}
	    } else if (atmnum == 3) {
		solute_1.rsolv[i__ - 1] = 1.432;
	    } else if (atmnum == 6) {
		solute_1.rsolv[i__ - 1] = 1.9;
		if (couple_1.n12[i__ - 1] == 3) {
		    solute_1.rsolv[i__ - 1] = 1.875;
		}
		if (couple_1.n12[i__ - 1] == 2) {
		    solute_1.rsolv[i__ - 1] = 1.825;
		}
	    } else if (atmnum == 7) {
		solute_1.rsolv[i__ - 1] = 1.7063;
		if (couple_1.n12[i__ - 1] == 4) {
		    solute_1.rsolv[i__ - 1] = 1.625;
		}
		if (couple_1.n12[i__ - 1] == 1) {
		    solute_1.rsolv[i__ - 1] = 1.6;
		}
	    } else if (atmnum == 8) {
		solute_1.rsolv[i__ - 1] = 1.535;
		if (couple_1.n12[i__ - 1] == 1) {
		    solute_1.rsolv[i__ - 1] = 1.48;
		}
	    } else if (atmnum == 9) {
		solute_1.rsolv[i__ - 1] = 1.47;
	    } else if (atmnum == 10) {
		solute_1.rsolv[i__ - 1] = 1.39;
	    } else if (atmnum == 11) {
		solute_1.rsolv[i__ - 1] = 1.992;
	    } else if (atmnum == 12) {
		solute_1.rsolv[i__ - 1] = 1.7;
	    } else if (atmnum == 14) {
		solute_1.rsolv[i__ - 1] = 1.8;
	    } else if (atmnum == 15) {
		solute_1.rsolv[i__ - 1] = 1.87;
	    } else if (atmnum == 16) {
		solute_1.rsolv[i__ - 1] = 1.775;
	    } else if (atmnum == 17) {
		solute_1.rsolv[i__ - 1] = 1.735;
	    } else if (atmnum == 18) {
		solute_1.rsolv[i__ - 1] = 1.7;
	    } else if (atmnum == 19) {
		solute_1.rsolv[i__ - 1] = 2.123;
	    } else if (atmnum == 20) {
		solute_1.rsolv[i__ - 1] = 1.817;
	    } else if (atmnum == 35) {
		solute_1.rsolv[i__ - 1] = 1.9;
	    } else if (atmnum == 36) {
		solute_1.rsolv[i__ - 1] = 1.812;
	    } else if (atmnum == 37) {
		solute_1.rsolv[i__ - 1] = 2.26;
	    } else if (atmnum == 53) {
		solute_1.rsolv[i__ - 1] = 2.1;
	    } else if (atmnum == 54) {
		solute_1.rsolv[i__ - 1] = 1.967;
	    } else if (atmnum == 55) {
		solute_1.rsolv[i__ - 1] = 2.507;
	    } else if (atmnum == 56) {
		solute_1.rsolv[i__ - 1] = 2.188;
	    } else {
		solute_1.rsolv[i__ - 1] = 2.;
	    }
	}
    }

/*     compute the atomic volumes for the analytical Still method */

    if (s_cmp(solute_1.solvtyp, "STILL", (ftnlen)8, (ftnlen)5) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 3rd power */
	    d__1 = solute_1.rsolv[i__ - 1];
	    solute_1.vsolv[i__ - 1] = d__1 * (d__1 * d__1) * 
		    4.1887902047863905;
	    ri = solute_1.rsolv[i__ - 1];
	    ri2 = ri * ri;
	    i__2 = couple_1.n12[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		k = i12_ref(j, i__);
		rk = solute_1.rsolv[k - 1];
		r__ = bond_1.bl[bndlist_ref(j, i__) - 1] * 1.01;
		ratio = (rk * rk - ri2 - r__ * r__) / (ri * 2. * r__);
		h__ = ri * (ratio + 1.);
		term = h__ * 1.0471975511965976 * h__ * (ri * 3. - h__);
		solute_1.vsolv[i__ - 1] -= term;
	    }
	}

/*     get self-, 1-2 and 1-3 polarization for analytical Still method */

	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.gpol[i__ - 1] = chgpot_1.electric * -.5 / (
		    solute_1.rsolv[i__ - 1] - solute_1.doffset + solute_1.p1);
	}
	i__1 = bond_1.nbond;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = ibnd_ref(1, i__);
	    ib = ibnd_ref(2, i__);
	    r__ = bond_1.bl[i__ - 1];
/* Computing 4th power */
	    d__1 = r__, d__1 *= d__1;
	    r4 = d__1 * d__1;
	    solute_1.gpol[ia - 1] += solute_1.p2 * solute_1.vsolv[ib - 1] / 
		    r4;
	    solute_1.gpol[ib - 1] += solute_1.p2 * solute_1.vsolv[ia - 1] / 
		    r4;
	}
	i__1 = angle_1.nangle;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = iang_ref(1, i__);
	    ib = iang_ref(2, i__);
	    ic = iang_ref(3, i__);
	    factor = 1.;
	    i__2 = couple_1.n12[ia - 1];
	    for (j = 1; j <= i__2; ++j) {
		id = i12_ref(j, ia);
		if (id == ic) {
		    factor = 0.;
		} else if (id != ib) {
		    i__3 = couple_1.n12[ic - 1];
		    for (k = 1; k <= i__3; ++k) {
			if (i12_ref(k, ic) == id) {
			    factor = .5;
			}
		    }
		}
	    }
	    i__2 = couple_1.n12[ib - 1];
	    for (j = 1; j <= i__2; ++j) {
		if (i12_ref(j, ib) == ia) {
		    rab = bond_1.bl[bndlist_ref(j, ib) - 1];
		} else if (i12_ref(j, ib) == ic) {
		    rbc = bond_1.bl[bndlist_ref(j, ib) - 1];
		}
	    }
	    cosine = cos(angle_1.anat[i__ - 1] / 57.29577951308232088);
/* Computing 2nd power */
	    d__1 = rab;
/* Computing 2nd power */
	    d__2 = rbc;
	    r2 = d__1 * d__1 + d__2 * d__2 - rab * 2. * rbc * cosine;
	    r4 = r2 * r2;
	    solute_1.gpol[ia - 1] += factor * solute_1.p3 * solute_1.vsolv[ic 
		    - 1] / r4;
	    solute_1.gpol[ic - 1] += factor * solute_1.p3 * solute_1.vsolv[ia 
		    - 1] / r4;
	}
    }

/*     assign the atomic radii and volumes for the ACE method; */
/*     volumes taken from average Voronoi values with hydrogens */

    if (s_cmp(solute_1.solvtyp, "ACE", (ftnlen)8, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    atmmas = i_dnnt(&atmtyp_1.mass[i__ - 1]);
	    if (atmnum == 1) {
		solute_1.rsolv[i__ - 1] = 1.468;
		solute_1.vsolv[i__ - 1] = 11.;
		k = i12_ref(1, i__);
		if (atmtyp_1.atomic[k - 1] == 6 && couple_1.n12[k - 1] == 4) {
		    solute_1.vsolv[i__ - 1] = 11.895;
		} else if (atmtyp_1.atomic[k - 1] == 6 && couple_1.n12[k - 1] 
			== 3) {
		    solute_1.vsolv[i__ - 1] = 13.242;
		} else if (atmtyp_1.atomic[k - 1] == 7 && couple_1.n12[k - 1] 
			== 4) {
		    solute_1.rsolv[i__ - 1] = .6;
		    solute_1.vsolv[i__ - 1] = 9.138;
		} else if (atmtyp_1.atomic[k - 1] == 7 || atmtyp_1.atomic[k - 
			1] == 8) {
		    solute_1.rsolv[i__ - 1] = .6;
		    solute_1.vsolv[i__ - 1] = 9.901;
		} else if (atmtyp_1.atomic[k - 1] != 16) {
		    solute_1.rsolv[i__ - 1] = 1.468;
		    solute_1.vsolv[i__ - 1] = 13.071;
		}
	    } else if (atmnum == 6) {
		solute_1.rsolv[i__ - 1] = 2.49;
		solute_1.vsolv[i__ - 1] = 7.;
		nh = 0;
		i__2 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, i__);
		    if (atmtyp_1.atomic[k - 1] == 1) {
			++nh;
		    }
		}
		if (couple_1.n12[i__ - 1] == 4) {
		    if (nh == 3) {
			solute_1.vsolv[i__ - 1] = 3.042;
		    } else if (nh == 2) {
			solute_1.vsolv[i__ - 1] = 3.743;
		    } else if (nh == 1) {
			solute_1.vsolv[i__ - 1] = 4.38;
		    }
		} else if (couple_1.n12[i__ - 1] == 3) {
		    if (nh == 1) {
			solute_1.rsolv[i__ - 1] = 2.1;
			solute_1.vsolv[i__ - 1] = 7.482;
		    } else if (nh == 0) {
			solute_1.rsolv[i__ - 1] = 2.1;
			solute_1.vsolv[i__ - 1] = 8.288;
		    }
		    i__2 = couple_1.n12[i__ - 1];
		    for (j = 1; j <= i__2; ++j) {
			k = i12_ref(1, j);
			if (atmtyp_1.atomic[k - 1] == 8 && couple_1.n12[k - 1]
				 == 1) {
			    solute_1.rsolv[i__ - 1] = 2.1;
			    solute_1.vsolv[i__ - 1] = 7.139;
			}
		    }
		}
		if (atmmas == 15) {
		    solute_1.rsolv[i__ - 1] = 2.165;
		    solute_1.vsolv[i__ - 1] = 33.175;
		} else if (atmmas == 14) {
		    solute_1.rsolv[i__ - 1] = 2.235;
		    solute_1.vsolv[i__ - 1] = 20.862;
		} else if (atmmas == 13 && couple_1.n12[i__ - 1] == 2) {
		    solute_1.rsolv[i__ - 1] = 2.1;
		    solute_1.vsolv[i__ - 1] = 20.329;
		} else if (atmmas == 13) {
		    solute_1.rsolv[i__ - 1] = 2.365;
		    solute_1.vsolv[i__ - 1] = 11.784;
		}
	    } else if (atmnum == 7) {
		solute_1.rsolv[i__ - 1] = 1.6;
		solute_1.vsolv[i__ - 1] = 6.;
		nh = 0;
		i__2 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, i__);
		    if (atmtyp_1.atomic[k - 1] == 1) {
			++nh;
		    }
		}
		if (couple_1.n12[i__ - 1] == 4) {
		    if (nh == 3) {
			solute_1.vsolv[i__ - 1] = 2.549;
		    } else if (nh == 2) {
			solute_1.vsolv[i__ - 1] = 3.304;
		    }
		} else if (couple_1.n12[i__ - 1] == 3) {
		    amide = FALSE_;
		    i__2 = couple_1.n12[i__ - 1];
		    for (j = 1; j <= i__2; ++j) {
			m = i12_ref(j, i__);
			if (atmtyp_1.atomic[m - 1] == 6) {
			    i__3 = couple_1.n12[m - 1];
			    for (k = 1; k <= i__3; ++k) {
				mm = i12_ref(k, m);
				if (atmtyp_1.atomic[mm - 1] == 8 && 
					couple_1.n12[mm - 1] == 1) {
				    amide = TRUE_;
				}
			    }
			}
		    }
		    if (amide) {
			if (nh == 0) {
			    solute_1.vsolv[i__ - 1] = 7.189;
			} else if (nh == 1) {
			    solute_1.vsolv[i__ - 1] = 6.03;
			} else if (nh == 2) {
			    solute_1.vsolv[i__ - 1] = 5.693;
			}
		    } else {
			if (nh == 2) {
			    solute_1.vsolv[i__ - 1] = 5.677;
			} else if (nh == 2) {
			    solute_1.vsolv[i__ - 1] = 6.498;
			}
		    }
		}
	    } else if (atmnum == 8) {
		solute_1.rsolv[i__ - 1] = 1.6;
		solute_1.vsolv[i__ - 1] = 12.;
		if (couple_1.n12[i__ - 1] == 1) {
		    solute_1.vsolv[i__ - 1] = 13.532;
		    k = i12_ref(1, i__);
		    if (atmtyp_1.atomic[k - 1] == 15) {
			solute_1.vsolv[i__ - 1] = 17.202;
		    } else {
			i__2 = couple_1.n13[i__ - 1];
			for (j = 1; j <= i__2; ++j) {
			    k = i13_ref(j, i__);
			    if (atmtyp_1.atomic[j - 1] == 8 && couple_1.n12[j 
				    - 1] == 1) {
				solute_1.vsolv[i__ - 1] = 15.4;
			    }
			}
		    }
		} else if (couple_1.n12[i__ - 1] == 2) {
		    solute_1.vsolv[i__ - 1] = 10.642;
		    i__2 = couple_1.n12[i__ - 1];
		    for (j = 1; j <= i__2; ++j) {
			k = i12_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 15) {
			    solute_1.vsolv[i__ - 1] = 11.416;
			}
		    }
		}
	    } else if (atmnum == 12) {
		solute_1.rsolv[i__ - 1] = 1.;
		solute_1.vsolv[i__ - 1] = 15.235;
	    } else if (atmnum == 15) {
		solute_1.rsolv[i__ - 1] = 1.89;
		solute_1.vsolv[i__ - 1] = 6.131;
	    } else if (atmnum == 16) {
		solute_1.rsolv[i__ - 1] = 1.89;
		solute_1.vsolv[i__ - 1] = 17.232;
		i__2 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, i__);
		    if (atmtyp_1.atomic[k - 1] == 16) {
			solute_1.vsolv[i__ - 1] = 18.465;
		    }
		}
	    } else if (atmnum == 26) {
		solute_1.rsolv[i__ - 1] = .65;
		solute_1.vsolv[i__ - 1] = 9.951;
	    } else {
		solute_1.rsolv[i__ - 1] = 0.;
		solute_1.vsolv[i__ - 1] = 0.;
	    }
	}

/*     calculate the pairwise parameters for the ACE method */

	c1 = .42441318157838759;
	c2 = .6681679418714731;
	c3 = 11.136655993663416;
	pi2 = .10132118364233778;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ic = atmtyp_1.class__[i__ - 1];
	    ri = solute_1.rsolv[i__ - 1];
	    ri2 = ri * ri;
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		kc = atmtyp_1.class__[k - 1];
		rk = solute_1.rsolv[kc - 1];
		vk = solute_1.vsolv[kc - 1];
		rk2 = rk * rk;
/* Computing MAX */
		d__1 = width, d__2 = ri / rk;
		alpha = max(d__1,d__2);
		alpha2 = alpha * alpha;
		alpha4 = alpha2 * alpha2;
		prod2 = alpha2 * rk2;
		prod4 = prod2 * prod2;
		ratio = alpha2 * rk2 / ri2;
		tik2 = ratio * 1.5707963267948966;
		temp = 1. / (tik2 * 2. + 1.);
		fik = 2. / (tik2 + 1.) - temp;
		qik = tik2 * sqrt(temp);
		qterm = qik - atan(qik);
		if (k != i__) {
		    omgik = vk * qterm * pi2 / prod4;
		} else {
		    omgik = c1 * qterm / (alpha4 * ri);
		}
		s2ik = qterm * 3. * prod2 / ((fik + 3.) * qik - atan(qik) * 
			4.);
		s3ik = s2ik * sqrt(s2ik);
		uik = c2 * ri / (1. - c3 * s3ik * ri * omgik / vk);
		wace_ref(ic, kc) = omgik;
		s2ace_ref(ic, kc) = s2ik;
		uace_ref(ic, kc) = uik;
	    }
	}
    }
    return 0;
} /* kgb_ */

#undef bndlist_ref
#undef s2ace_ref
#undef wace_ref
#undef iang_ref
#undef uace_ref
#undef ibnd_ref
#undef i13_ref
#undef i12_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine kgk  --  set generalized Kirkwood parameters  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "kgk" initializes parameters needed for the generalized */
/*     Kirkwood solvation model */


/* Subroutine */ int kgk_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);

    /* Local variables */
    static integer i__, j, k, l, m, next;
    static char value[20];
    static doublereal rscale;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static doublereal offset;
    static integer atmnum;
    static char radtyp[10], string[120];
    extern /* Subroutine */ int getword_(char *, char *, integer *, ftnlen, 
	    ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___73 = { 1, string, 1, 0, 120, 1 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
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

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kvdws.i  --  forcefield parameters for van der Waals terms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     rad      van der Waals radius parameter for each atom type */
/*     eps      van der Waals well depth parameter for each atom type */
/*     rad4     van der Waals radius parameter in 1-4 interactions */
/*     eps4     van der Waals well depth parameter in 1-4 interactions */
/*     reduct   van der Waals reduction factor for each atom type */




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




/*     set default value for exponent in the GB/GK function */

    gk_1.gkc = 2.455;
    s_copy(radtyp, "BONDI", (ftnlen)10, (ftnlen)5);

/*     get any altered generalized Kirkwood values from keyfile */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "GKC ", (ftnlen)4, (ftnlen)4) == 0) {
	    i__2 = s_rsli(&io___73);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gk_1.gkc, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "GK-RADII ", (ftnlen)9, (ftnlen)9) == 0) {
	    getword_(record, value, &next, (ftnlen)120, (ftnlen)20);
	    upcase_(value, (ftnlen)20);
	    if (s_cmp(value, "VDW", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(radtyp, "VDW", (ftnlen)10, (ftnlen)3);
	    } else if (s_cmp(value, "MACROMODEL", (ftnlen)10, (ftnlen)10) == 
		    0) {
		s_copy(radtyp, "MACROMODEL", (ftnlen)10, (ftnlen)10);
	    } else if (s_cmp(value, "AMOEBA", (ftnlen)6, (ftnlen)6) == 0) {
		s_copy(radtyp, "AMOEBA", (ftnlen)10, (ftnlen)6);
	    } else if (s_cmp(value, "BONDI", (ftnlen)5, (ftnlen)5) == 0) {
		s_copy(radtyp, "BONDI", (ftnlen)10, (ftnlen)5);
	    } else if (s_cmp(value, "TOMASI", (ftnlen)6, (ftnlen)6) == 0) {
		s_copy(radtyp, "TOMASI", (ftnlen)10, (ftnlen)6);
	    }
	}
L10:
	;
    }

/*     assign base atomic radii from the van der Waals values */

    if (s_cmp(radtyp, "VDW", (ftnlen)10, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.rsolv[i__ - 1] = 2.;
	    if (atmtyp_1.class__[i__ - 1] != 0) {
		solute_1.rsolv[i__ - 1] = kvdws_1.rad[atmtyp_1.class__[i__ - 
			1] - 1];
	    }
	    solute_1.rsolv[i__ - 1] += -.1;
	}

/*     assign standard solvation radii adapted from Macromodel */

    } else if (s_cmp(radtyp, "MACROMODEL", (ftnlen)10, (ftnlen)10) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.rsolv[i__ - 1] = 2.;
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    if (atmnum == 1) {
		solute_1.rsolv[i__ - 1] = 1.25;
		k = i12_ref(1, i__);
		if (atmtyp_1.atomic[k - 1] == 7) {
		    solute_1.rsolv[i__ - 1] = 1.15;
		}
		if (atmtyp_1.atomic[k - 1] == 8) {
		    solute_1.rsolv[i__ - 1] = 1.05;
		}
	    } else if (atmnum == 3) {
		solute_1.rsolv[i__ - 1] = 1.432;
	    } else if (atmnum == 6) {
		solute_1.rsolv[i__ - 1] = 1.9;
		if (couple_1.n12[i__ - 1] == 3) {
		    solute_1.rsolv[i__ - 1] = 1.875;
		}
		if (couple_1.n12[i__ - 1] == 2) {
		    solute_1.rsolv[i__ - 1] = 1.825;
		}
	    } else if (atmnum == 7) {
		solute_1.rsolv[i__ - 1] = 1.7063;
		if (couple_1.n12[i__ - 1] == 4) {
		    solute_1.rsolv[i__ - 1] = 1.625;
		}
		if (couple_1.n12[i__ - 1] == 1) {
		    solute_1.rsolv[i__ - 1] = 1.6;
		}
	    } else if (atmnum == 8) {
		solute_1.rsolv[i__ - 1] = 1.535;
		if (couple_1.n12[i__ - 1] == 1) {
		    solute_1.rsolv[i__ - 1] = 1.48;
		}
	    } else if (atmnum == 9) {
		solute_1.rsolv[i__ - 1] = 1.47;
	    } else if (atmnum == 10) {
		solute_1.rsolv[i__ - 1] = 1.39;
	    } else if (atmnum == 11) {
		solute_1.rsolv[i__ - 1] = 1.992;
	    } else if (atmnum == 12) {
		solute_1.rsolv[i__ - 1] = 1.7;
	    } else if (atmnum == 14) {
		solute_1.rsolv[i__ - 1] = 1.8;
	    } else if (atmnum == 15) {
		solute_1.rsolv[i__ - 1] = 1.87;
	    } else if (atmnum == 16) {
		solute_1.rsolv[i__ - 1] = 1.775;
	    } else if (atmnum == 17) {
		solute_1.rsolv[i__ - 1] = 1.735;
	    } else if (atmnum == 18) {
		solute_1.rsolv[i__ - 1] = 1.7;
	    } else if (atmnum == 19) {
		solute_1.rsolv[i__ - 1] = 2.123;
	    } else if (atmnum == 20) {
		solute_1.rsolv[i__ - 1] = 1.817;
	    } else if (atmnum == 35) {
		solute_1.rsolv[i__ - 1] = 1.9;
	    } else if (atmnum == 36) {
		solute_1.rsolv[i__ - 1] = 1.812;
	    } else if (atmnum == 37) {
		solute_1.rsolv[i__ - 1] = 2.26;
	    } else if (atmnum == 53) {
		solute_1.rsolv[i__ - 1] = 2.1;
	    } else if (atmnum == 54) {
		solute_1.rsolv[i__ - 1] = 1.967;
	    } else if (atmnum == 55) {
		solute_1.rsolv[i__ - 1] = 2.507;
	    } else if (atmnum == 56) {
		solute_1.rsolv[i__ - 1] = 2.188;
	    }
	}

/*     assign base atomic radii from traditional Bondi values */

    } else if (s_cmp(radtyp, "AMOEBA", (ftnlen)10, (ftnlen)6) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.rsolv[i__ - 1] = 2.;
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    if (atmnum == 1) {
		solute_1.rsolv[i__ - 1] = 1.32;
		k = i12_ref(1, i__);
		if (atmtyp_1.atomic[k - 1] == 7) {
		    solute_1.rsolv[i__ - 1] = 1.1;
		}
		if (atmtyp_1.atomic[k - 1] == 8) {
		    solute_1.rsolv[i__ - 1] = 1.05;
		}
	    }
	    if (atmnum == 3) {
		solute_1.rsolv[i__ - 1] = 1.5;
	    }
	    if (atmnum == 6) {
		solute_1.rsolv[i__ - 1] = 2.;
		if (couple_1.n12[i__ - 1] == 3) {
		    solute_1.rsolv[i__ - 1] = 2.05;
		}
		if (couple_1.n12[i__ - 1] == 4) {
		    i__2 = couple_1.n12[i__ - 1];
		    for (j = 1; j <= i__2; ++j) {
			k = i12_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 7) {
			    solute_1.rsolv[i__ - 1] = 1.75;
			}
			if (atmtyp_1.atomic[k - 1] == 8) {
			    solute_1.rsolv[i__ - 1] = 1.75;
			}
		    }
		}
	    }
	    if (atmnum == 7) {
		solute_1.rsolv[i__ - 1] = 1.6;
	    }
	    if (atmnum == 8) {
		solute_1.rsolv[i__ - 1] = 1.55;
		if (couple_1.n12[i__ - 1] == 2) {
		    solute_1.rsolv[i__ - 1] = 1.45;
		}
	    }
	    if (atmnum == 9) {
		solute_1.rsolv[i__ - 1] = 1.54;
	    }
	    if (atmnum == 10) {
		solute_1.rsolv[i__ - 1] = 1.46;
	    }
	    if (atmnum == 11) {
		solute_1.rsolv[i__ - 1] = 2.09;
	    }
	    if (atmnum == 12) {
		solute_1.rsolv[i__ - 1] = 1.79;
	    }
	    if (atmnum == 14) {
		solute_1.rsolv[i__ - 1] = 1.89;
	    }
	    if (atmnum == 15) {
		solute_1.rsolv[i__ - 1] = 1.96;
	    }
	    if (atmnum == 16) {
		solute_1.rsolv[i__ - 1] = 1.86;
	    }
	    if (atmnum == 17) {
		solute_1.rsolv[i__ - 1] = 1.82;
	    }
	    if (atmnum == 18) {
		solute_1.rsolv[i__ - 1] = 1.79;
	    }
	    if (atmnum == 19) {
		solute_1.rsolv[i__ - 1] = 2.23;
	    }
	    if (atmnum == 20) {
		solute_1.rsolv[i__ - 1] = 1.91;
	    }
	    if (atmnum == 35) {
		solute_1.rsolv[i__ - 1] = 2.;
	    }
	    if (atmnum == 36) {
		solute_1.rsolv[i__ - 1] = 1.9;
	    }
	    if (atmnum == 37) {
		solute_1.rsolv[i__ - 1] = 2.26;
	    }
	    if (atmnum == 53) {
		solute_1.rsolv[i__ - 1] = 2.37;
	    }
	    if (atmnum == 54) {
		solute_1.rsolv[i__ - 1] = 2.07;
	    }
	    if (atmnum == 55) {
		solute_1.rsolv[i__ - 1] = 2.63;
	    }
	    if (atmnum == 56) {
		solute_1.rsolv[i__ - 1] = 2.3;
	    }
	}

/*     assign base atomic radii from traditional Bondi values */

    } else {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.rsolv[i__ - 1] = 2.;
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    if (atmnum == 0) {
		solute_1.rsolv[i__ - 1] = 0.;
	    }
	    if (atmnum == 1) {
		solute_1.rsolv[i__ - 1] = 1.2;
	    }
	    if (atmnum == 2) {
		solute_1.rsolv[i__ - 1] = 1.4;
	    }
	    if (atmnum == 5) {
		solute_1.rsolv[i__ - 1] = 1.8;
	    }
	    if (atmnum == 6) {
		solute_1.rsolv[i__ - 1] = 1.7;
	    }
	    if (atmnum == 7) {
		solute_1.rsolv[i__ - 1] = 1.55;
	    }
	    if (atmnum == 8) {
		solute_1.rsolv[i__ - 1] = 1.52;
	    }
	    if (atmnum == 9) {
		solute_1.rsolv[i__ - 1] = 1.47;
	    }
	    if (atmnum == 10) {
		solute_1.rsolv[i__ - 1] = 1.54;
	    }
	    if (atmnum == 14) {
		solute_1.rsolv[i__ - 1] = 2.1;
	    }
	    if (atmnum == 15) {
		solute_1.rsolv[i__ - 1] = 1.8;
	    }
	    if (atmnum == 16) {
		solute_1.rsolv[i__ - 1] = 1.8;
	    }
	    if (atmnum == 17) {
		solute_1.rsolv[i__ - 1] = 1.75;
	    }
	    if (atmnum == 18) {
		solute_1.rsolv[i__ - 1] = 1.88;
	    }
	    if (atmnum == 34) {
		solute_1.rsolv[i__ - 1] = 1.9;
	    }
	    if (atmnum == 35) {
		solute_1.rsolv[i__ - 1] = 1.85;
	    }
	    if (atmnum == 36) {
		solute_1.rsolv[i__ - 1] = 2.02;
	    }
	    if (atmnum == 53) {
		solute_1.rsolv[i__ - 1] = 1.98;
	    }
	    if (atmnum == 54) {
		solute_1.rsolv[i__ - 1] = 2.16;
	    }
	}
    }

/*     make Tomasi-style modifications to the base atomic radii */

    if (s_cmp(radtyp, "TOMASI", (ftnlen)10, (ftnlen)6) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    offset = 0.;
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    if (atmtyp_1.atomic[i__ - 1] == 1) {
		i__2 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, i__);
		    if (atmtyp_1.atomic[k - 1] == 6) {
			i__3 = couple_1.n12[k - 1];
			for (l = 1; l <= i__3; ++l) {
			    m = i12_ref(l, k);
			    if (atmtyp_1.atomic[m - 1] == 7) {
				offset = -.05;
			    }
			    if (atmtyp_1.atomic[m - 1] == 8) {
				offset = -.1;
			    }
			}
		    }
		    if (atmtyp_1.atomic[k - 1] == 7) {
			offset = -.25;
		    }
		    if (atmtyp_1.atomic[k - 1] == 8) {
			offset = -.4;
		    }
		    if (atmtyp_1.atomic[k - 1] == 16) {
			offset = -.1;
		    }
		}
	    } else if (atmtyp_1.atomic[i__ - 1] == 6) {
		if (couple_1.n12[i__ - 1] == 4) {
		    offset = .05;
		}
		if (couple_1.n12[i__ - 1] == 3) {
		    offset = .02;
		}
		if (couple_1.n12[i__ - 1] == 2) {
		    offset = -.03;
		}
		i__2 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, i__);
		    if (atmtyp_1.atomic[k - 1] == 6) {
			offset += -.07;
		    }
		}
		i__2 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, i__);
		    if (atmtyp_1.atomic[k - 1] == 7 && couple_1.n12[k - 1] == 
			    4) {
			offset = -.2;
		    }
		    if (atmtyp_1.atomic[k - 1] == 7 && couple_1.n12[k - 1] == 
			    3) {
			offset = -.25;
		    }
		    if (atmtyp_1.atomic[k - 1] == 8) {
			offset = -.2;
		    }
		}
	    } else if (atmtyp_1.atomic[i__ - 1] == 7) {
		if (couple_1.n12[i__ - 1] == 3) {
		    offset = -.1;
		    i__2 = couple_1.n12[i__ - 1];
		    for (j = 1; j <= i__2; ++j) {
			k = i12_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 6) {
			    offset += -.24;
			}
		    }
		} else {
		    offset = -.2;
		    i__2 = couple_1.n12[i__ - 1];
		    for (j = 1; j <= i__2; ++j) {
			k = i12_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 6) {
			    offset += -.16;
			}
		    }
		}
	    } else if (atmtyp_1.atomic[i__ - 1] == 8) {
		if (couple_1.n12[i__ - 1] == 2) {
		    offset = -.21;
		    i__2 = couple_1.n12[i__ - 1];
		    for (j = 1; j <= i__2; ++j) {
			k = i12_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 6) {
			    offset = -.36;
			}
		    }
		} else {
		    offset = -.25;
		}
	    } else if (atmtyp_1.atomic[i__ - 1] == 16) {
		offset = -.03;
		i__2 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, i__);
		    if (atmtyp_1.atomic[k - 1] == 6) {
			offset += -.1;
		    }
		}
	    }
	    solute_1.rsolv[i__ - 1] += offset;
	}
    }

/*     apply an overall scale factor to the solvation radii */

    rscale = 1.;
    if (s_cmp(radtyp, "VDW", (ftnlen)10, (ftnlen)3) == 0) {
	rscale = 1.;
    }
    if (s_cmp(radtyp, "MACROMODEL", (ftnlen)10, (ftnlen)10) == 0) {
	rscale = 1.;
    }
    if (s_cmp(radtyp, "AMOEBA", (ftnlen)10, (ftnlen)6) == 0) {
	rscale = 1.;
    }
    if (s_cmp(radtyp, "BONDI", (ftnlen)10, (ftnlen)5) == 0) {
	rscale = 1.03;
    }
    if (s_cmp(radtyp, "TOMASI", (ftnlen)10, (ftnlen)6) == 0) {
	rscale = 1.24;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	solute_1.rsolv[i__ - 1] *= rscale;
    }

/*     assign generic value for the HCT overlap scale factor */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	solute_1.shct[i__ - 1] = .69;
    }
    if (s_cmp(radtyp, "MACROMODEL", (ftnlen)10, (ftnlen)10) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.shct[i__ - 1] = .8;
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    if (atmnum == 1) {
		solute_1.shct[i__ - 1] = .85;
	    }
	    if (atmnum == 6) {
		solute_1.shct[i__ - 1] = .72;
	    }
	    if (atmnum == 7) {
		solute_1.shct[i__ - 1] = .79;
	    }
	    if (atmnum == 8) {
		solute_1.shct[i__ - 1] = .85;
	    }
	    if (atmnum == 9) {
		solute_1.shct[i__ - 1] = .88;
	    }
	    if (atmnum == 15) {
		solute_1.shct[i__ - 1] = .86;
	    }
	    if (atmnum == 16) {
		solute_1.shct[i__ - 1] = .96;
	    }
	    if (atmnum == 26) {
		solute_1.shct[i__ - 1] = .88;
	    }
	}
    }
    return 0;
} /* kgk_ */

#undef keyline_ref
#undef i12_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine kpb  --  assign Poisson-Boltzmann parameters  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "kpb" assigns parameters needed for the Poisson-Boltzmann */
/*     solvation model implemented via APBS */


/* Subroutine */ int kpb_(void)
{
    /* Format strings */
    static char fmt_160[] = "(/,\002 APBS Grid Dimensions and Spacing :\002,"
	    "//,10x,3i8,10x,f10.4)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer pbtyplen;
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j, k, l, m, pbsolnlen;
    static doublereal ri, gx, gy, gz;
    static integer nx, ny, nz;
    static doublereal xcm, ycm, zcm;
    extern /* Subroutine */ int apbsinitial_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, char *, integer *, char *, integer *, 
	    char *, integer *, char *, integer *, char *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen, ftnlen);
    static doublereal xmin, xmax;
    static integer next;
    static doublereal ymin, ymax, zmin, zmax, xlen, ylen, zlen, weigh;
    static char value[20];
    static doublereal total, rscale, pbionc;
    static char record[120];
    static integer maxgrd, atmnum, pbionq;
    static doublereal offset, pbionr;
    static char radtyp[10], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static integer bcfllen, chgmlen;
    static doublereal spacing;
    static integer srfmlen;
    extern /* Subroutine */ int getword_(char *, char *, integer *, ftnlen, 
	    ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___108 = { 1, string, 1, 0, 120, 1 };
    static icilist io___110 = { 1, string, 1, 0, 120, 1 };
    static icilist io___111 = { 1, string, 1, 0, 120, 1 };
    static icilist io___112 = { 1, string, 1, 0, 120, 1 };
    static icilist io___113 = { 1, string, 1, 0, 120, 1 };
    static icilist io___114 = { 1, string, 1, 0, 120, 1 };
    static icilist io___115 = { 1, string, 1, 0, 120, 1 };
    static icilist io___119 = { 1, string, 1, 0, 120, 1 };
    static icilist io___120 = { 1, string, 1, 0, 120, 1 };
    static icilist io___122 = { 1, string, 1, 0, 120, 1 };
    static icilist io___126 = { 1, string, 1, 0, 120, 1 };
    static icilist io___127 = { 1, string, 1, 0, 120, 1 };
    static icilist io___128 = { 1, string, 1, 0, 120, 1 };
    static icilist io___129 = { 1, string, 1, 0, 120, 1 };
    static icilist io___130 = { 1, string, 1, 0, 120, 1 };
    static cilist io___142 = { 0, 0, 0, fmt_160, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
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
/*     ##  bath.i  --  temperature and pressure control parameters  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxnose     maximum length of the Nose-Hoover chain */

/*     kelvin0     target value for the system temperature (K) */
/*     kelvin      variable target temperature for thermostat (K) */
/*     atmsph      target value for the system pressure (atm) */
/*     tautemp     time constant for Berendsen thermostat (psec) */
/*     taupres     time constant for Berendsen barostat (psec) */
/*     compress    isothermal compressibility of medium (atm-1) */
/*     collide     collision frequency for Andersen thermostat */
/*     xnh         position of each chained Nose-Hoover thermostat */
/*     vnh         velocity of each chained Nose-Hoover thermostat */
/*     qnh         mass for each chained Nose-Hoover thermostat */
/*     gnh         coupling between chained Nose-Hoover thermostats */
/*     volmove     maximum volume move for Monte Carlo barostat (Ang**3) */
/*     voltrial    mean number of steps between Monte Carlo moves */
/*     isothermal  logical flag governing use of temperature control */
/*     isobaric    logical flag governing use of pressure control */
/*     anisotrop   logical flag governing use of anisotropic pressure */
/*     thermostat  choice of temperature control method to be used */
/*     barostat    choice of pressure control method to be used */
/*     volscale    choice of scaling method for Monte Carlo barostat */




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

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kvdws.i  --  forcefield parameters for van der Waals terms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     rad      van der Waals radius parameter for each atom type */
/*     eps      van der Waals well depth parameter for each atom type */
/*     rad4     van der Waals radius parameter in 1-4 interactions */
/*     eps4     van der Waals well depth parameter in 1-4 interactions */
/*     reduct   van der Waals reduction factor for each atom type */




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




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  npolar.i  --  nonpolar cavity & dispersion parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     epso      water oxygen eps for implicit dispersion term */
/*     epsh      water hydrogen eps for implicit dispersion term */
/*     rmino     water oxygen Rmin for implicit dispersion term */
/*     rminh     water hydrogen Rmin for implicit dispersion term */
/*     awater    water number density at standard temp & pressure */
/*     slevy     enthalpy-to-free energy scale factor for dispersion */

/*     solvprs   limiting microscopic solvent pressure value */
/*     surften   limiting macroscopic surface tension value */
/*     spcut     starting radius for solvent pressure tapering */
/*     spoff     cutoff radius for solvent pressure tapering */
/*     stcut     starting radius for surface tension tapering */
/*     stoff     cutoff radius for surface tension tapering */
/*     rcav      atomic radius of each atom for cavitation energy */
/*     rdisp     atomic radius of each atom for dispersion energy */
/*     cdisp     maximum dispersion energy for each atom */




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




/*     assign some default APBS configuration parameters */

    s_copy(pb_1.pbtyp, "LPBE", (ftnlen)20, (ftnlen)4);
    s_copy(pb_1.pbsoln, "MG-MANUAL", (ftnlen)20, (ftnlen)9);
    s_copy(radtyp, "BONDI", (ftnlen)10, (ftnlen)5);
    s_copy(pb_1.bcfl, "MDH", (ftnlen)20, (ftnlen)3);
    s_copy(pb_1.chgm, "SPL4", (ftnlen)20, (ftnlen)4);
    s_copy(pb_1.srfm, "SPL4", (ftnlen)20, (ftnlen)4);
    bath_1.kelvin = 298.;
    pb_1.pdie = 1.;
    pb_1.sdie = 78.3;
    pb_1.srad = 1.4;
    pb_1.swin = .3;
    pb_1.sdens = 10.;
    pb_1.smin = 3.;
    pb_1.ionn = 0;
    for (i__ = 1; i__ <= 10; ++i__) {
	pb_1.ionc[i__ - 1] = 0.;
	pb_1.ionq[i__ - 1] = 1;
	pb_1.ionr[i__ - 1] = 2.;
    }
    spacing = .5;
    maxgrd = 225;

/*     compute the position of the center of mass */

    total = 0.;
    xcm = 0.;
    ycm = 0.;
    zcm = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	weigh = atmtyp_1.mass[i__ - 1];
	total += weigh;
	xcm += atoms_1.x[i__ - 1] * weigh;
	ycm += atoms_1.y[i__ - 1] * weigh;
	zcm += atoms_1.z__[i__ - 1] * weigh;
    }
    xcm /= total;
    ycm /= total;
    zcm /= total;
    pb_1.gcent[0] = xcm;
    pb_1.gcent[1] = ycm;
    pb_1.gcent[2] = zcm;

/*     set default APBS grid dimension based on system extent */

    xmin = xcm;
    ymin = ycm;
    zmin = zcm;
    xmax = xcm;
    ymax = ycm;
    zmax = zcm;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ri = solute_1.rsolv[i__ - 1];
/* Computing MIN */
	d__1 = xmin, d__2 = atoms_1.x[i__ - 1] - ri;
	xmin = min(d__1,d__2);
/* Computing MIN */
	d__1 = ymin, d__2 = atoms_1.y[i__ - 1] - ri;
	ymin = min(d__1,d__2);
/* Computing MIN */
	d__1 = zmin, d__2 = atoms_1.z__[i__ - 1] - ri;
	zmin = min(d__1,d__2);
/* Computing MAX */
	d__1 = xmax, d__2 = atoms_1.x[i__ - 1] + ri;
	xmax = max(d__1,d__2);
/* Computing MAX */
	d__1 = ymax, d__2 = atoms_1.y[i__ - 1] + ri;
	ymax = max(d__1,d__2);
/* Computing MAX */
	d__1 = zmax, d__2 = atoms_1.z__[i__ - 1] + ri;
	zmax = max(d__1,d__2);
    }
/* Computing MAX */
    d__1 = xcm - xmin, d__2 = xmax - xcm;
    xlen = (max(d__1,d__2) + pb_1.smin) * 2.;
/* Computing MAX */
    d__1 = ycm - ymin, d__2 = ymax - ycm;
    ylen = (max(d__1,d__2) + pb_1.smin) * 2.;
/* Computing MAX */
    d__1 = zcm - zmin, d__2 = zmax - zcm;
    zlen = (max(d__1,d__2) + pb_1.smin) * 2.;
    pb_1.dime[0] = (integer) (xlen / spacing) + 1;
    pb_1.dime[1] = (integer) (ylen / spacing) + 1;
    pb_1.dime[2] = (integer) (zlen / spacing) + 1;

/*     get any altered APBS parameters from the keyfile */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "MG-AUTO ", (ftnlen)8, (ftnlen)8) == 0) {
	    s_copy(pb_1.pbsoln, "MG-AUTO", (ftnlen)20, (ftnlen)7);
	} else if (s_cmp(keyword, "MG-MANUAL ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    s_copy(pb_1.pbsoln, "MG-MANUAL", (ftnlen)20, (ftnlen)9);
	} else if (s_cmp(keyword, "APBS-GRID ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    nx = pb_1.dime[0];
	    ny = pb_1.dime[1];
	    nz = pb_1.dime[2];
	    i__2 = s_rsli(&io___108);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&nx, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ny, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&nz, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
L10:
	    if (nx >= 33) {
		pb_1.dime[0] = nx;
	    }
	    if (ny >= 33) {
		pb_1.dime[1] = ny;
	    }
	    if (nz >= 33) {
		pb_1.dime[2] = nz;
	    }
	} else if (s_cmp(keyword, "PB-RADII ", (ftnlen)9, (ftnlen)9) == 0) {
	    getword_(record, value, &next, (ftnlen)120, (ftnlen)20);
	    upcase_(value, (ftnlen)20);
	    if (s_cmp(value, "VDW", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(radtyp, "VDW", (ftnlen)10, (ftnlen)3);
	    } else if (s_cmp(value, "MACROMODEL", (ftnlen)10, (ftnlen)10) == 
		    0) {
		s_copy(radtyp, "MACROMODEL", (ftnlen)10, (ftnlen)10);
	    } else if (s_cmp(value, "BONDI", (ftnlen)5, (ftnlen)5) == 0) {
		s_copy(radtyp, "BONDI", (ftnlen)10, (ftnlen)5);
	    } else if (s_cmp(value, "TOMASI", (ftnlen)6, (ftnlen)6) == 0) {
		s_copy(radtyp, "TOMASI", (ftnlen)10, (ftnlen)6);
	    }
	} else if (s_cmp(keyword, "SDENS ", (ftnlen)6, (ftnlen)6) == 0) {
	    i__2 = s_rsli(&io___110);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&pb_1.sdens, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
L20:
	    ;
	} else if (s_cmp(keyword, "PDIE ", (ftnlen)5, (ftnlen)5) == 0) {
	    i__2 = s_rsli(&io___111);
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&pb_1.pdie, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L30;
	    }
L30:
	    ;
	} else if (s_cmp(keyword, "SDIE ", (ftnlen)5, (ftnlen)5) == 0) {
	    i__2 = s_rsli(&io___112);
	    if (i__2 != 0) {
		goto L40;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&pb_1.sdie, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L40;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L40;
	    }
L40:
	    ;
	} else if (s_cmp(keyword, "SRAD ", (ftnlen)5, (ftnlen)5) == 0) {
	    i__2 = s_rsli(&io___113);
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&pb_1.srad, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L50;
	    }
L50:
	    ;
	} else if (s_cmp(keyword, "SWIN ", (ftnlen)5, (ftnlen)5) == 0) {
	    i__2 = s_rsli(&io___114);
	    if (i__2 != 0) {
		goto L60;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&pb_1.swin, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L60;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L60;
	    }
L60:
	    ;
	} else if (s_cmp(keyword, "SMIN ", (ftnlen)5, (ftnlen)5) == 0) {
	    i__2 = s_rsli(&io___115);
	    if (i__2 != 0) {
		goto L70;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&pb_1.smin, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L70;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L70;
	    }
L70:
	    ;
	} else if (s_cmp(keyword, "SRFM ", (ftnlen)5, (ftnlen)5) == 0) {
	    getword_(record, value, &next, (ftnlen)120, (ftnlen)20);
	    upcase_(value, (ftnlen)20);
	    if (s_cmp(value, "MOL", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(pb_1.srfm, "MOL", (ftnlen)20, (ftnlen)3);
	    } else if (s_cmp(value, "SMOL", (ftnlen)4, (ftnlen)4) == 0) {
		s_copy(pb_1.srfm, "SMOL", (ftnlen)20, (ftnlen)4);
	    } else if (s_cmp(value, "SPL2", (ftnlen)4, (ftnlen)4) == 0) {
		s_copy(pb_1.srfm, "SPL2", (ftnlen)20, (ftnlen)4);
	    }
	} else if (s_cmp(keyword, "BCFL ", (ftnlen)5, (ftnlen)5) == 0) {
	    getword_(record, value, &next, (ftnlen)120, (ftnlen)20);
	    upcase_(value, (ftnlen)20);
	    if (s_cmp(value, "ZERO", (ftnlen)3, (ftnlen)4) == 0) {
		s_copy(pb_1.bcfl, "ZERO", (ftnlen)20, (ftnlen)4);
	    } else if (s_cmp(value, "MDH", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(pb_1.bcfl, "MDH", (ftnlen)20, (ftnlen)3);
	    } else if (s_cmp(value, "SDH", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(pb_1.bcfl, "SDH", (ftnlen)20, (ftnlen)3);
	    }
	} else if (s_cmp(keyword, "ION ", (ftnlen)4, (ftnlen)4) == 0) {
	    pbionc = 0.;
	    pbionq = 1;
	    pbionr = 2.;
	    i__2 = s_rsli(&io___119);
	    if (i__2 != 0) {
		goto L80;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&pbionq, (ftnlen)sizeof(
		    integer));
	    if (i__2 != 0) {
		goto L80;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&pbionc, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L80;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&pbionr, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L80;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L80;
	    }
L80:
	    if (pbionq != 0 && pbionc >= 0. && pbionr >= 0.) {
		++pb_1.ionn;
		pb_1.ionc[pb_1.ionn - 1] = pbionc;
		pb_1.ionq[pb_1.ionn - 1] = pbionq;
		pb_1.ionr[pb_1.ionn - 1] = pbionr;
	    }
	}
    }

/*     set APBS grid spacing for the chosen grid dimension */

/* Computing MAX */
    d__1 = xcm - xmin, d__2 = xmax - xcm;
    xlen = (max(d__1,d__2) + pb_1.smin) * 2.;
/* Computing MAX */
    d__1 = ycm - ymin, d__2 = ymax - ycm;
    ylen = (max(d__1,d__2) + pb_1.smin) * 2.;
/* Computing MAX */
    d__1 = zcm - zmin, d__2 = zmax - zcm;
    zlen = (max(d__1,d__2) + pb_1.smin) * 2.;
    pb_1.grid[0] = xlen / pb_1.dime[0];
    pb_1.grid[1] = ylen / pb_1.dime[1];
    pb_1.grid[2] = zlen / pb_1.dime[2];

/*     grid spacing must be equal to maintain traceless quadrupoles */

/* Computing MIN */
    d__1 = min(pb_1.grid[0],pb_1.grid[1]);
    pb_1.grid[0] = min(d__1,pb_1.grid[2]);
    pb_1.grid[1] = pb_1.grid[0];
    pb_1.grid[2] = pb_1.grid[0];

/*     set the grid dimensions to the smallest multiples of 32 */

    pb_1.dime[0] = 33;
    pb_1.dime[1] = 33;
    pb_1.dime[2] = 33;
    while(pb_1.grid[0] * pb_1.dime[0] < xlen) {
	pb_1.dime[0] += 32;
    }
    while(pb_1.grid[1] * pb_1.dime[1] < ylen) {
	pb_1.dime[1] += 32;
    }
    while(pb_1.grid[2] * pb_1.dime[2] < zlen) {
	pb_1.dime[2] += 32;
    }

/*     limit the grid dimensions and recompute the grid spacing */

    pb_1.dime[0] = min(pb_1.dime[0],maxgrd);
    pb_1.dime[1] = min(pb_1.dime[1],maxgrd);
    pb_1.dime[2] = min(pb_1.dime[2],maxgrd);
    pb_1.grid[0] = xlen / pb_1.dime[0];
    pb_1.grid[1] = ylen / pb_1.dime[1];
    pb_1.grid[2] = zlen / pb_1.dime[2];

/*     grid spacing must be equal to maintain traceless quadrupoles */

/* Computing MAX */
    d__1 = max(pb_1.grid[0],pb_1.grid[1]);
    pb_1.grid[0] = max(d__1,pb_1.grid[2]);
    pb_1.grid[1] = pb_1.grid[0];
    pb_1.grid[2] = pb_1.grid[0];

/*     if this is an "mg-auto" (focusing) calculation, set the */
/*     fine grid to the default size, and the coarse grid to */
/*     twice its original size. Currently, all energies and */
/*     forces need to be evaluated at the same resolution */

    if (s_cmp(pb_1.pbsoln, "MG-AUTO", (ftnlen)20, (ftnlen)7) == 0) {
	pb_1.fgrid[0] = pb_1.grid[0];
	pb_1.fgrid[1] = pb_1.grid[1];
	pb_1.fgrid[2] = pb_1.grid[2];
	pb_1.fgcent[0] = pb_1.gcent[0];
	pb_1.fgcent[1] = pb_1.gcent[1];
	pb_1.fgcent[2] = pb_1.gcent[2];
	pb_1.cgrid[0] = pb_1.grid[0] * 2.;
	pb_1.cgrid[1] = pb_1.grid[1] * 2.;
	pb_1.cgrid[2] = pb_1.grid[2] * 2.;
    }

/*     get any custom APBS grid parameters from the keyfile */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "DIME ", (ftnlen)5, (ftnlen)5) == 0) {
	    i__2 = s_rsli(&io___120);
	    if (i__2 != 0) {
		goto L90;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&nx, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L90;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ny, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L90;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&nz, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L90;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L90;
	    }
	    pb_1.dime[0] = nx;
	    pb_1.dime[1] = ny;
	    pb_1.dime[2] = nz;
L90:
	    for (j = 1; j <= 3; ++j) {
		if (pb_1.dime[j - 1] % 32 != 1) {
		    pb_1.dime[j - 1] = ((pb_1.dime[j - 1] - 1) / 32 + 1 << 5) 
			    + 1;
		}
	    }
	} else if (s_cmp(keyword, "GRID ", (ftnlen)5, (ftnlen)5) == 0) {
	    i__2 = s_rsli(&io___122);
	    if (i__2 != 0) {
		goto L100;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gx, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L100;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gy, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L100;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gz, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L100;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L100;
	    }
	    pb_1.grid[0] = gx;
	    pb_1.grid[1] = gy;
	    pb_1.grid[2] = gz;
L100:
	    ;
	} else if (s_cmp(keyword, "CGRID ", (ftnlen)6, (ftnlen)6) == 0) {
	    i__2 = s_rsli(&io___126);
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gx, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gy, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gz, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L110;
	    }
	    pb_1.cgrid[0] = gx;
	    pb_1.cgrid[1] = gy;
	    pb_1.cgrid[2] = gz;
L110:
	    ;
	} else if (s_cmp(keyword, "FGRID ", (ftnlen)6, (ftnlen)6) == 0) {
	    i__2 = s_rsli(&io___127);
	    if (i__2 != 0) {
		goto L120;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gx, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L120;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gy, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L120;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gz, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L120;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L120;
	    }
	    pb_1.fgrid[0] = gx;
	    pb_1.fgrid[1] = gy;
	    pb_1.fgrid[2] = gz;
L120:
	    ;
	} else if (s_cmp(keyword, "GCENT ", (ftnlen)6, (ftnlen)6) == 0) {
	    i__2 = s_rsli(&io___128);
	    if (i__2 != 0) {
		goto L130;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gx, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L130;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gy, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L130;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gz, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L130;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L130;
	    }
	    pb_1.gcent[0] = gx;
	    pb_1.gcent[1] = gy;
	    pb_1.gcent[2] = gz;
L130:
	    ;
	} else if (s_cmp(keyword, "CGCENT ", (ftnlen)7, (ftnlen)7) == 0) {
	    i__2 = s_rsli(&io___129);
	    if (i__2 != 0) {
		goto L140;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gx, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L140;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gy, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L140;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gz, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L140;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L140;
	    }
	    pb_1.cgcent[0] = gx;
	    pb_1.cgcent[1] = gy;
	    pb_1.cgcent[2] = gz;
L140:
	    ;
	} else if (s_cmp(keyword, "FGCENT ", (ftnlen)7, (ftnlen)7) == 0) {
	    i__2 = s_rsli(&io___130);
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gx, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gy, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gz, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L150;
	    }
	    pb_1.fgcent[0] = gx;
	    pb_1.fgcent[1] = gy;
	    pb_1.fgcent[2] = gz;
L150:
	    ;
	}
    }

/*     assign base atomic radii from the van der Waals values */

    if (s_cmp(radtyp, "VDW", (ftnlen)10, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.rsolv[i__ - 1] = 2.;
	    if (atmtyp_1.class__[i__ - 1] != 0) {
		solute_1.rsolv[i__ - 1] = kvdws_1.rad[atmtyp_1.class__[i__ - 
			1] - 1];
	    }
	}

/*     assign standard solvation radii adapted from Macromodel */

    } else if (s_cmp(radtyp, "MACROMODEL", (ftnlen)10, (ftnlen)10) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.rsolv[i__ - 1] = 2.;
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    if (atmnum == 1) {
		solute_1.rsolv[i__ - 1] = 1.25;
		k = i12_ref(1, i__);
		if (atmtyp_1.atomic[k - 1] == 7) {
		    solute_1.rsolv[i__ - 1] = 1.15;
		}
		if (atmtyp_1.atomic[k - 1] == 8) {
		    solute_1.rsolv[i__ - 1] = 1.05;
		}
	    } else if (atmnum == 3) {
		solute_1.rsolv[i__ - 1] = 1.432;
	    } else if (atmnum == 6) {
		solute_1.rsolv[i__ - 1] = 1.9;
		if (couple_1.n12[i__ - 1] == 3) {
		    solute_1.rsolv[i__ - 1] = 1.875;
		}
		if (couple_1.n12[i__ - 1] == 2) {
		    solute_1.rsolv[i__ - 1] = 1.825;
		}
	    } else if (atmnum == 7) {
		solute_1.rsolv[i__ - 1] = 1.7063;
		if (couple_1.n12[i__ - 1] == 4) {
		    solute_1.rsolv[i__ - 1] = 1.625;
		}
		if (couple_1.n12[i__ - 1] == 1) {
		    solute_1.rsolv[i__ - 1] = 1.6;
		}
	    } else if (atmnum == 8) {
		solute_1.rsolv[i__ - 1] = 1.535;
		if (couple_1.n12[i__ - 1] == 1) {
		    solute_1.rsolv[i__ - 1] = 1.48;
		}
	    } else if (atmnum == 9) {
		solute_1.rsolv[i__ - 1] = 1.47;
	    } else if (atmnum == 10) {
		solute_1.rsolv[i__ - 1] = 1.39;
	    } else if (atmnum == 11) {
		solute_1.rsolv[i__ - 1] = 1.992;
	    } else if (atmnum == 12) {
		solute_1.rsolv[i__ - 1] = 1.7;
	    } else if (atmnum == 14) {
		solute_1.rsolv[i__ - 1] = 1.8;
	    } else if (atmnum == 15) {
		solute_1.rsolv[i__ - 1] = 1.87;
	    } else if (atmnum == 16) {
		solute_1.rsolv[i__ - 1] = 1.775;
	    } else if (atmnum == 17) {
		solute_1.rsolv[i__ - 1] = 1.735;
	    } else if (atmnum == 18) {
		solute_1.rsolv[i__ - 1] = 1.7;
	    } else if (atmnum == 19) {
		solute_1.rsolv[i__ - 1] = 2.123;
	    } else if (atmnum == 20) {
		solute_1.rsolv[i__ - 1] = 1.817;
	    } else if (atmnum == 35) {
		solute_1.rsolv[i__ - 1] = 1.9;
	    } else if (atmnum == 36) {
		solute_1.rsolv[i__ - 1] = 1.812;
	    } else if (atmnum == 37) {
		solute_1.rsolv[i__ - 1] = 2.26;
	    } else if (atmnum == 53) {
		solute_1.rsolv[i__ - 1] = 2.1;
	    } else if (atmnum == 54) {
		solute_1.rsolv[i__ - 1] = 1.967;
	    } else if (atmnum == 55) {
		solute_1.rsolv[i__ - 1] = 2.507;
	    } else if (atmnum == 56) {
		solute_1.rsolv[i__ - 1] = 2.188;
	    }
	}

/*     assign base atomic radii from traditional Bondi values */

    } else {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.rsolv[i__ - 1] = 2.;
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    if (atmnum == 0) {
		solute_1.rsolv[i__ - 1] = 0.;
	    }
	    if (atmnum == 1) {
		solute_1.rsolv[i__ - 1] = 1.2;
	    }
	    if (atmnum == 2) {
		solute_1.rsolv[i__ - 1] = 1.4;
	    }
	    if (atmnum == 5) {
		solute_1.rsolv[i__ - 1] = 1.8;
	    }
	    if (atmnum == 6) {
		solute_1.rsolv[i__ - 1] = 1.7;
	    }
	    if (atmnum == 7) {
		solute_1.rsolv[i__ - 1] = 1.55;
	    }
	    if (atmnum == 8) {
		solute_1.rsolv[i__ - 1] = 1.52;
	    }
	    if (atmnum == 9) {
		solute_1.rsolv[i__ - 1] = 1.47;
	    }
	    if (atmnum == 10) {
		solute_1.rsolv[i__ - 1] = 1.54;
	    }
	    if (atmnum == 14) {
		solute_1.rsolv[i__ - 1] = 2.1;
	    }
	    if (atmnum == 15) {
		solute_1.rsolv[i__ - 1] = 1.8;
	    }
	    if (atmnum == 16) {
		solute_1.rsolv[i__ - 1] = 1.8;
	    }
	    if (atmnum == 17) {
		solute_1.rsolv[i__ - 1] = 1.75;
	    }
	    if (atmnum == 18) {
		solute_1.rsolv[i__ - 1] = 1.88;
	    }
	    if (atmnum == 34) {
		solute_1.rsolv[i__ - 1] = 1.9;
	    }
	    if (atmnum == 35) {
		solute_1.rsolv[i__ - 1] = 1.85;
	    }
	    if (atmnum == 36) {
		solute_1.rsolv[i__ - 1] = 2.02;
	    }
	    if (atmnum == 53) {
		solute_1.rsolv[i__ - 1] = 1.98;
	    }
	    if (atmnum == 54) {
		solute_1.rsolv[i__ - 1] = 2.16;
	    }
	}
    }

/*     make Tomasi-style modifications to the base atomic radii */

    if (s_cmp(radtyp, "TOMASI", (ftnlen)10, (ftnlen)6) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    offset = 0.;
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    if (atmtyp_1.atomic[i__ - 1] == 1) {
		i__2 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, i__);
		    if (atmtyp_1.atomic[k - 1] == 6) {
			i__3 = couple_1.n12[k - 1];
			for (l = 1; l <= i__3; ++l) {
			    m = i12_ref(l, k);
			    if (atmtyp_1.atomic[m - 1] == 7) {
				offset = -.05;
			    }
			    if (atmtyp_1.atomic[m - 1] == 8) {
				offset = -.1;
			    }
			}
		    }
		    if (atmtyp_1.atomic[k - 1] == 7) {
			offset = -.25;
		    }
		    if (atmtyp_1.atomic[k - 1] == 8) {
			offset = -.4;
		    }
		    if (atmtyp_1.atomic[k - 1] == 16) {
			offset = -.1;
		    }
		}
	    } else if (atmtyp_1.atomic[i__ - 1] == 6) {
		if (couple_1.n12[i__ - 1] == 4) {
		    offset = .05;
		}
		if (couple_1.n12[i__ - 1] == 3) {
		    offset = .02;
		}
		if (couple_1.n12[i__ - 1] == 2) {
		    offset = -.03;
		}
		i__2 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, i__);
		    if (atmtyp_1.atomic[k - 1] == 6) {
			offset += -.07;
		    }
		}
		i__2 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, i__);
		    if (atmtyp_1.atomic[k - 1] == 7 && couple_1.n12[k - 1] == 
			    4) {
			offset = -.2;
		    }
		    if (atmtyp_1.atomic[k - 1] == 7 && couple_1.n12[k - 1] == 
			    3) {
			offset = -.25;
		    }
		    if (atmtyp_1.atomic[k - 1] == 8) {
			offset = -.2;
		    }
		}
	    } else if (atmtyp_1.atomic[i__ - 1] == 7) {
		if (couple_1.n12[i__ - 1] == 3) {
		    offset = -.1;
		    i__2 = couple_1.n12[i__ - 1];
		    for (j = 1; j <= i__2; ++j) {
			k = i12_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 6) {
			    offset += -.24;
			}
		    }
		} else {
		    offset = -.2;
		    i__2 = couple_1.n12[i__ - 1];
		    for (j = 1; j <= i__2; ++j) {
			k = i12_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 6) {
			    offset += -.16;
			}
		    }
		}
	    } else if (atmtyp_1.atomic[i__ - 1] == 8) {
		if (couple_1.n12[i__ - 1] == 2) {
		    offset = -.21;
		    i__2 = couple_1.n12[i__ - 1];
		    for (j = 1; j <= i__2; ++j) {
			k = i12_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 6) {
			    offset = -.36;
			}
		    }
		} else {
		    offset = -.25;
		}
	    } else if (atmtyp_1.atomic[i__ - 1] == 16) {
		offset = -.03;
		i__2 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, i__);
		    if (atmtyp_1.atomic[k - 1] == 6) {
			offset += -.1;
		    }
		}
	    }
	    solute_1.rsolv[i__ - 1] += offset;
	}
    }

/*     apply an overall scale factor to the solvation radii */

    rscale = 1.;
    if (s_cmp(radtyp, "VDW", (ftnlen)10, (ftnlen)3) == 0) {
	rscale = 1.1;
    }
    if (s_cmp(radtyp, "MACROMODEL", (ftnlen)10, (ftnlen)10) == 0) {
	rscale = 1.15;
    }
    if (s_cmp(radtyp, "BONDI", (ftnlen)10, (ftnlen)5) == 0) {
	rscale = 1.21;
    }
    if (s_cmp(radtyp, "TOMASI", (ftnlen)10, (ftnlen)6) == 0) {
	rscale = 1.47;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	solute_1.rsolv[i__ - 1] *= rscale;
    }

/*     assign generic value for the HCT overlap scale factor */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	solute_1.shct[i__ - 1] = .69;
    }

/*     determine the length of the character arguments */

    pbtyplen = trimtext_(pb_1.pbtyp, (ftnlen)20);
    pbsolnlen = trimtext_(pb_1.pbsoln, (ftnlen)20);
    bcfllen = trimtext_(pb_1.bcfl, (ftnlen)20);
    chgmlen = trimtext_(pb_1.chgm, (ftnlen)20);
    srfmlen = trimtext_(pb_1.srfm, (ftnlen)20);

/*     make call needed to initialize the APBS calculation */

    if (inform_1.verbose) {
	io___142.ciunit = iounit_1.iout;
	s_wsfe(&io___142);
	for (i__ = 1; i__ <= 3; ++i__) {
	    do_fio(&c__1, (char *)&pb_1.dime[i__ - 1], (ftnlen)sizeof(integer)
		    );
	}
	do_fio(&c__1, (char *)&pb_1.grid[0], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    apbsinitial_(pb_1.dime, pb_1.grid, pb_1.gcent, pb_1.cgrid, pb_1.cgcent, 
	    pb_1.fgrid, pb_1.fgcent, &pb_1.pdie, &pb_1.sdie, &pb_1.srad, &
	    pb_1.swin, &pb_1.sdens, &bath_1.kelvin, &pb_1.ionn, pb_1.ionc, 
	    pb_1.ionq, pb_1.ionr, pb_1.pbtyp, &pbtyplen, pb_1.pbsoln, &
	    pbsolnlen, pb_1.bcfl, &bcfllen, pb_1.chgm, &chgmlen, pb_1.srfm, &
	    srfmlen, (ftnlen)20, (ftnlen)20, (ftnlen)20, (ftnlen)20, (ftnlen)
	    20);
    return 0;
} /* kpb_ */

#undef keyline_ref
#undef i12_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine knp  --  assign nonpolar solvation parameters  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "knp" initializes parameters needed for the cavitation-plus- */
/*     dispersion nonpolar portion of solvation models */


/* Subroutine */ int knp_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal ah, ao, ri, ri3, ri7, ri11, epsi;
    static integer next;
    static doublereal emixh, rmini, emixo, rmixh, cross, rmixo, rmixh3, 
	    rmixh7, rmixo3, rmixo7, cavoff;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    static doublereal dispoff;
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___150 = { 1, string, 1, 0, 120, 1 };
    static icilist io___151 = { 1, string, 1, 0, 120, 1 };



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
/*     ##  keys.i  --  contents of current keyword parameter file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nkey      number of nonblank lines in the keyword file */
/*     keyline   contents of each individual keyword file line */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kvdws.i  --  forcefield parameters for van der Waals terms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     rad      van der Waals radius parameter for each atom type */
/*     eps      van der Waals well depth parameter for each atom type */
/*     rad4     van der Waals radius parameter in 1-4 interactions */
/*     eps4     van der Waals well depth parameter in 1-4 interactions */
/*     reduct   van der Waals reduction factor for each atom type */




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




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  npolar.i  --  nonpolar cavity & dispersion parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     epso      water oxygen eps for implicit dispersion term */
/*     epsh      water hydrogen eps for implicit dispersion term */
/*     rmino     water oxygen Rmin for implicit dispersion term */
/*     rminh     water hydrogen Rmin for implicit dispersion term */
/*     awater    water number density at standard temp & pressure */
/*     slevy     enthalpy-to-free energy scale factor for dispersion */

/*     solvprs   limiting microscopic solvent pressure value */
/*     surften   limiting macroscopic surface tension value */
/*     spcut     starting radius for solvent pressure tapering */
/*     spoff     cutoff radius for solvent pressure tapering */
/*     stcut     starting radius for surface tension tapering */
/*     stoff     cutoff radius for surface tension tapering */
/*     rcav      atomic radius of each atom for cavitation energy */
/*     rdisp     atomic radius of each atom for dispersion energy */
/*     cdisp     maximum dispersion energy for each atom */




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




/*     set default values for solvent pressure and surface tension */

    npolar_1.solvprs = .0327;
    npolar_1.surften = .08;

/*     set default values for cavity and dispersion radius offsets */

    cavoff = 0.;
    dispoff = .26;

/*     get any altered surface tension value from keyfile */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "SOLVENT-PRESSURE ", (ftnlen)17, (ftnlen)17) == 0) 
		{
	    i__2 = s_rsli(&io___150);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&npolar_1.solvprs, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "SURFACE-TENSION ", (ftnlen)16, (ftnlen)16) 
		== 0) {
	    i__2 = s_rsli(&io___151);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&npolar_1.surften, (ftnlen)
		    sizeof(doublereal));
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

/*     set switching function values for pressure and tension */

    cross = npolar_1.surften * 3. / npolar_1.solvprs;
    npolar_1.spcut = cross - 3.5;
    npolar_1.spoff = cross + 3.5;
    npolar_1.stcut = cross + 3.9;
    npolar_1.stoff = cross - 3.5;

/*     assign surface area factors for nonpolar solvation */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	solute_1.asolv[i__ - 1] = npolar_1.surften;
    }

/*     set cavity and dispersion radii for nonpolar solvation */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	npolar_1.rcav[i__ - 1] = kvdws_1.rad[atmtyp_1.class__[i__ - 1] - 1] + 
		cavoff;
	npolar_1.rdisp[i__ - 1] = kvdws_1.rad[atmtyp_1.class__[i__ - 1] - 1] 
		+ dispoff;
    }

/*     compute maximum dispersion energies for each atom */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	epsi = kvdws_1.eps[atmtyp_1.class__[i__ - 1] - 1];
	if (npolar_1.rdisp[i__ - 1] > 0. && epsi > 0.) {
	    rmini = kvdws_1.rad[atmtyp_1.class__[i__ - 1] - 1];
/* Computing 2nd power */
	    d__1 = sqrt(.11) + sqrt(epsi);
	    emixo = epsi * .44 / (d__1 * d__1);
/* Computing 3rd power */
	    d__1 = rmini;
/* Computing 2nd power */
	    d__2 = rmini;
	    rmixo = (d__1 * (d__1 * d__1) + 4.9347068906249989) * 2. / (d__2 *
		     d__2 + 2.8985062499999996);
/* Computing 3rd power */
	    d__1 = rmixo;
	    rmixo3 = d__1 * (d__1 * d__1);
/* Computing 7th power */
	    d__1 = rmixo, d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
	    rmixo7 = d__2 * (d__1 * d__1);
	    ao = emixo * rmixo7;
/* Computing 2nd power */
	    d__1 = sqrt(.0135) + sqrt(epsi);
	    emixh = epsi * .053999999999999999 / (d__1 * d__1);
/* Computing 3rd power */
	    d__1 = rmini;
/* Computing 2nd power */
	    d__2 = rmini;
	    rmixh = (d__1 * (d__1 * d__1) + 2.3393951718749992) * 2. / (d__2 *
		     d__2 + 1.7622562499999996);
/* Computing 3rd power */
	    d__1 = rmixh;
	    rmixh3 = d__1 * (d__1 * d__1);
/* Computing 7th power */
	    d__1 = rmixh, d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
	    rmixh7 = d__2 * (d__1 * d__1);
	    ah = emixh * rmixh7;
	    ri = npolar_1.rdisp[i__ - 1];
/* Computing 3rd power */
	    d__1 = ri;
	    ri3 = d__1 * (d__1 * d__1);
/* Computing 7th power */
	    d__1 = ri, d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
	    ri7 = d__2 * (d__1 * d__1);
/* Computing 11th power */
	    d__1 = ri, d__2 = d__1, d__1 *= d__1, d__2 *= d__1, d__1 *= d__1;
	    ri11 = d__2 * (d__1 * d__1);
	    if (ri < rmixh) {
		npolar_1.cdisp[i__ - 1] = emixh * -12.566370614359172 * (
			rmixh3 - ri3) / 3.;
		npolar_1.cdisp[i__ - 1] -= emixh * 18. / 11. * rmixh3 * 
			3.141592653589793238;
	    } else {
		npolar_1.cdisp[i__ - 1] = (rmixh7 * 2. - ri7 * 11.) * 
			6.2831853071795862 * ah;
		npolar_1.cdisp[i__ - 1] /= ri11 * 11.;
	    }
	    npolar_1.cdisp[i__ - 1] *= 2.;
	    if (ri < rmixo) {
		npolar_1.cdisp[i__ - 1] -= emixo * 12.566370614359172 * (
			rmixo3 - ri3) / 3.;
		npolar_1.cdisp[i__ - 1] -= emixo * 18. / 11. * rmixo3 * 
			3.141592653589793238;
	    } else {
		npolar_1.cdisp[i__ - 1] += (rmixo7 * 2. - ri7 * 11.) * 
			6.2831853071795862 * ao / (ri11 * 11.);
	    }
	}
	npolar_1.cdisp[i__ - 1] *= .033427999999999999;
    }
    return 0;
} /* knp_ */

#undef keyline_ref


