/* kdipole.f -- translated by f2c (version 20050501).
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
    integer bndlist[200000]	/* was [8][25000] */, anglist[700000]	/* 
	    was [28][25000] */;
} atmlst_;

#define atmlst_1 atmlst_

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
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal bdpl[50000], sdpl[50000];
    integer ndipole, idpl[100000]	/* was [2][50000] */;
} dipole_;

#define dipole_1 dipole_

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
    doublereal dpl[1000], dpl5[500], dpl4[500], dpl3[500], pos[1000], pos5[
	    500], pos4[500], pos3[500];
    char kd[8000], kd5[4000], kd4[4000], kd3[4000];
} kdipol_;

#define kdipol_1 kdipol_

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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine kdipole  --  assign bond dipole parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "kdipole" assigns bond dipoles to the bonds within */
/*     the structure and processes any new or changed values */


/* Subroutine */ int kdipole_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Additional Bond Dipole Moment \002,\002P"
	    "arameters :\002,//,5x,\002Atom Types\002,9x,\002Moment\002,5x"
	    ",\002Position\002,/)";
    static char fmt_30[] = "(6x,2i4,4x,2f12.3)";
    static char fmt_40[] = "(6x,2i4,4x,2f12.3,3x,a6)";
    static char fmt_50[] = "(/,\002 KDIPOLE  --  Too many Bond Dipole\002"
	    ",\002 Moment Parameters\002)";
    static char fmt_60[] = "(/,\002 KDIPOLE  --  Too many 5-Ring Bond\002"
	    ",\002 Dipole Parameters\002)";
    static char fmt_70[] = "(/,\002 KDIPOLE  --  Too many 4-Ring Bond\002"
	    ",\002 Dipole Parameters\002)";
    static char fmt_80[] = "(/,\002 KDIPOLE  --  Too many 3-Ring Bond\002"
	    ",\002 Dipole Parameters\002)";
    static char fmt_120[] = "(/,\002 Additional Bond Dipoles for\002,\002 Sp"
	    "ecific Bonds :\002,//,5x,\002Bonded Atoms\002,7x,\002Moment\002,"
	    "5x,\002Position\002,/)";
    static char fmt_130[] = "(4x,i5,\002 -\002,i5,2x,2f12.3)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static logical use_ring__;
    static integer i__, j, k, ia, ib, nd;
    static doublereal dp;
    static char pa[4], pb[4];
    static doublereal ps;
    static char pt[8];
    static integer nd4, nd5, nd3, ita, itb, size, next;
    static char label[6], blank[8];
    static integer iring;
    static logical header;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int chkring_(integer *, integer *, integer *, 
	    integer *, integer *), numeral_(integer *, char *, integer *, 
	    ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___13 = { 1, string, 1, 0, 120, 1 };
    static cilist io___14 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_80, 0 };
    static icilist io___34 = { 1, string, 1, 0, 120, 1 };
    static cilist io___35 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_130, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define kd_ref(a_0,a_1) &kdipol_1.kd[(a_1)*8 + a_0 - 8]
#define kd3_ref(a_0,a_1) &kdipol_1.kd3[(a_1)*8 + a_0 - 8]
#define kd4_ref(a_0,a_1) &kdipol_1.kd4[(a_1)*8 + a_0 - 8]
#define kd5_ref(a_0,a_1) &kdipol_1.kd5[(a_1)*8 + a_0 - 8]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define idpl_ref(a_1,a_2) dipole_1.idpl[(a_2)*2 + a_1 - 3]
#define bndlist_ref(a_1,a_2) atmlst_1.bndlist[(a_2)*8 + a_1 - 9]
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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  kdipol.i  --  forcefield parameters for bond dipoles  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxnd    maximum number of bond dipole parameter entries */
/*     maxnd5   maximum number of 5-membered ring dipole entries */
/*     maxnd4   maximum number of 4-membered ring dipole entries */
/*     maxnd3   maximum number of 3-membered ring dipole entries */

/*     dpl      dipole moment parameters for bond dipoles */
/*     dpl5     dipole moment parameters for 5-ring dipoles */
/*     dpl4     dipole moment parameters for 4-ring dipoles */
/*     dpl3     dipole moment parameters for 3-ring dipoles */
/*     pos      dipole position parameters for bond dipoles */
/*     pos5     dipole position parameters for 5-ring dipoles */
/*     pos4     dipole position parameters for 4-ring dipoles */
/*     pos3     dipole position parameters for 3-ring dipoles */
/*     kd       string of atom classes for bond dipoles */
/*     kd5      string of atom classes for 5-ring dipoles */
/*     kd4      string of atom classes for 4-ring dipoles */
/*     kd3      string of atom classes for 3-ring dipoles */




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




/*     process keywords containing bond dipole parameters */

    s_copy(blank, "        ", (ftnlen)8, (ftnlen)8);
    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	iring = -1;
	if (s_cmp(keyword, "DIPOLE ", (ftnlen)7, (ftnlen)7) == 0) {
	    iring = 0;
	}
	if (s_cmp(keyword, "DIPOLE5 ", (ftnlen)8, (ftnlen)8) == 0) {
	    iring = 5;
	}
	if (s_cmp(keyword, "DIPOLE4 ", (ftnlen)8, (ftnlen)8) == 0) {
	    iring = 4;
	}
	if (s_cmp(keyword, "DIPOLE3 ", (ftnlen)8, (ftnlen)8) == 0) {
	    iring = 3;
	}
	if (iring >= 0) {
	    ia = 0;
	    ib = 0;
	    dp = 0.;
	    ps = .5;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___13);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&dp, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&ps, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
L10:
	    if (ia > 0 && ib > 0) {
		if (header) {
		    header = FALSE_;
		    io___14.ciunit = iounit_1.iout;
		    s_wsfe(&io___14);
		    e_wsfe();
		}
		if (iring == 0) {
		    io___15.ciunit = iounit_1.iout;
		    s_wsfe(&io___15);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&dp, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&ps, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else {
		    if (iring == 5) {
			s_copy(label, "5-Ring", (ftnlen)6, (ftnlen)6);
		    }
		    if (iring == 4) {
			s_copy(label, "4-Ring", (ftnlen)6, (ftnlen)6);
		    }
		    if (iring == 3) {
			s_copy(label, "3-Ring", (ftnlen)6, (ftnlen)6);
		    }
		    io___17.ciunit = iounit_1.iout;
		    s_wsfe(&io___17);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&dp, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&ps, (ftnlen)sizeof(doublereal));
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
		    for (j = 1; j <= 1000; ++j) {
			if (s_cmp(kd_ref(0, j), blank, (ftnlen)8, (ftnlen)8) 
				== 0 || s_cmp(kd_ref(0, j), pt, (ftnlen)8, (
				ftnlen)8) == 0) {
			    s_copy(kd_ref(0, j), pt, (ftnlen)8, (ftnlen)8);
			    if (ia <= ib) {
				kdipol_1.dpl[j - 1] = dp;
				kdipol_1.pos[j - 1] = ps;
			    } else {
				kdipol_1.dpl[j - 1] = -dp;
				kdipol_1.pos[j - 1] = 1. - ps;
			    }
			    goto L90;
			}
		    }
		    io___23.ciunit = iounit_1.iout;
		    s_wsfe(&io___23);
		    e_wsfe();
		    inform_1.abort = TRUE_;
		} else if (iring == 5) {
		    for (j = 1; j <= 500; ++j) {
			if (s_cmp(kd5_ref(0, j), blank, (ftnlen)8, (ftnlen)8) 
				== 0 || s_cmp(kd5_ref(0, j), pt, (ftnlen)8, (
				ftnlen)8) == 0) {
			    s_copy(kd5_ref(0, j), pt, (ftnlen)8, (ftnlen)8);
			    if (ia <= ib) {
				kdipol_1.dpl5[j - 1] = dp;
				kdipol_1.pos5[j - 1] = ps;
			    } else {
				kdipol_1.dpl5[j - 1] = -dp;
				kdipol_1.pos5[j - 1] = 1. - ps;
			    }
			    goto L90;
			}
		    }
		    io___24.ciunit = iounit_1.iout;
		    s_wsfe(&io___24);
		    e_wsfe();
		    inform_1.abort = TRUE_;
		} else if (iring == 4) {
		    for (j = 1; j <= 500; ++j) {
			if (s_cmp(kd4_ref(0, j), blank, (ftnlen)8, (ftnlen)8) 
				== 0 || s_cmp(kd4_ref(0, j), pt, (ftnlen)8, (
				ftnlen)8) == 0) {
			    s_copy(kd4_ref(0, j), pt, (ftnlen)8, (ftnlen)8);
			    if (ia <= ib) {
				kdipol_1.dpl4[j - 1] = dp;
				kdipol_1.pos4[j - 1] = ps;
			    } else {
				kdipol_1.dpl4[j - 1] = -dp;
				kdipol_1.pos4[j - 1] = 1. - ps;
			    }
			    goto L90;
			}
		    }
		    io___25.ciunit = iounit_1.iout;
		    s_wsfe(&io___25);
		    e_wsfe();
		    inform_1.abort = TRUE_;
		} else if (iring == 3) {
		    for (j = 1; j <= 500; ++j) {
			if (s_cmp(kd3_ref(0, j), blank, (ftnlen)8, (ftnlen)8) 
				== 0 || s_cmp(kd3_ref(0, j), pt, (ftnlen)8, (
				ftnlen)8) == 0) {
			    s_copy(kd3_ref(0, j), pt, (ftnlen)8, (ftnlen)8);
			    if (ia <= ib) {
				kdipol_1.dpl3[j - 1] = dp;
				kdipol_1.pos3[j - 1] = ps;
			    } else {
				kdipol_1.dpl3[j - 1] = -dp;
				kdipol_1.pos3[j - 1] = 1. - ps;
			    }
			    goto L90;
			}
		    }
		    io___26.ciunit = iounit_1.iout;
		    s_wsfe(&io___26);
		    e_wsfe();
		    inform_1.abort = TRUE_;
		}
	    }
L90:
	    ;
	}
    }

/*     determine the total number of forcefield parameters */

    nd = 1000;
    nd5 = 500;
    nd4 = 500;
    nd3 = 500;
    for (i__ = 1000; i__ >= 1; --i__) {
	if (s_cmp(kd_ref(0, i__), blank, (ftnlen)8, (ftnlen)8) == 0) {
	    nd = i__ - 1;
	}
    }
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(kd5_ref(0, i__), blank, (ftnlen)8, (ftnlen)8) == 0) {
	    nd5 = i__ - 1;
	}
    }
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(kd4_ref(0, i__), blank, (ftnlen)8, (ftnlen)8) == 0) {
	    nd4 = i__ - 1;
	}
    }
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(kd3_ref(0, i__), blank, (ftnlen)8, (ftnlen)8) == 0) {
	    nd3 = i__ - 1;
	}
    }
    use_ring__ = FALSE_;
/* Computing MIN */
    i__1 = min(nd5,nd4);
    if (min(i__1,nd3) != 0) {
	use_ring__ = TRUE_;
    }

/*     find and store all the bond dipole moments */

    i__1 = bond_1.nbond;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ibnd_ref(1, i__);
	ib = ibnd_ref(2, i__);
	ita = atoms_1.type__[ia - 1];
	itb = atoms_1.type__[ib - 1];
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
	dipole_1.bdpl[i__ - 1] = 0.;

/*     make a check for bonds contained inside small rings */

	iring = 0;
	if (use_ring__) {
	    chkring_(&iring, &ia, &ib, &c__0, &c__0);
	    if (iring == 6) {
		iring = 0;
	    }
	    if (iring == 5 && nd5 == 0) {
		iring = 0;
	    }
	    if (iring == 4 && nd4 == 0) {
		iring = 0;
	    }
	    if (iring == 3 && nd3 == 0) {
		iring = 0;
	    }
	}

/*     try to assign bond dipole parameters for the bond */

	if (iring == 0) {
	    i__2 = nd;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kd_ref(0, j), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    if (ita <= itb) {
			idpl_ref(1, i__) = ia;
			idpl_ref(2, i__) = ib;
		    } else {
			idpl_ref(1, i__) = ib;
			idpl_ref(2, i__) = ia;
		    }
		    dipole_1.bdpl[i__ - 1] = kdipol_1.dpl[j - 1];
		    dipole_1.sdpl[i__ - 1] = kdipol_1.pos[j - 1];
		    goto L100;
		}
	    }
	} else if (iring == 5) {
	    i__2 = nd5;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kd5_ref(0, j), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    if (ita <= itb) {
			idpl_ref(1, i__) = ia;
			idpl_ref(2, i__) = ib;
		    } else {
			idpl_ref(1, i__) = ib;
			idpl_ref(2, i__) = ia;
		    }
		    dipole_1.bdpl[i__ - 1] = kdipol_1.dpl5[j - 1];
		    dipole_1.sdpl[i__ - 1] = kdipol_1.pos5[j - 1];
		    goto L100;
		}
	    }
	} else if (iring == 4) {
	    i__2 = nd4;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kd4_ref(0, j), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    if (ita <= itb) {
			idpl_ref(1, i__) = ia;
			idpl_ref(2, i__) = ib;
		    } else {
			idpl_ref(1, i__) = ib;
			idpl_ref(2, i__) = ia;
		    }
		    dipole_1.bdpl[i__ - 1] = kdipol_1.dpl4[j - 1];
		    dipole_1.sdpl[i__ - 1] = kdipol_1.pos4[j - 1];
		    goto L100;
		}
	    }
	} else if (iring == 3) {
	    i__2 = nd3;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kd3_ref(0, j), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    if (ita <= itb) {
			idpl_ref(1, i__) = ia;
			idpl_ref(2, i__) = ib;
		    } else {
			idpl_ref(1, i__) = ib;
			idpl_ref(2, i__) = ia;
		    }
		    dipole_1.bdpl[i__ - 1] = kdipol_1.dpl3[j - 1];
		    dipole_1.sdpl[i__ - 1] = kdipol_1.pos3[j - 1];
		    goto L100;
		}
	    }
	}
L100:
	;
    }

/*     process keywords containing bond specific bond dipoles */

    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "DIPOLE ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    dp = 0.;
	    ps = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___34);
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&dp, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&ps, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L110;
	    }
L110:
	    if (ia < 0 && ib < 0) {
		ia = -ia;
		ib = -ib;
		if (header) {
		    header = FALSE_;
		    io___35.ciunit = iounit_1.iout;
		    s_wsfe(&io___35);
		    e_wsfe();
		}
		i__2 = couple_1.n12[ia - 1];
		for (j = 1; j <= i__2; ++j) {
		    if (i12_ref(j, ia) == ib) {
			k = bndlist_ref(j, ia);
			if (ps == 0.) {
			    ps = .5;
			}
			if (idpl_ref(1, k) == ib) {
			    dipole_1.bdpl[k - 1] = dp;
			    dipole_1.sdpl[k - 1] = ps;
			} else {
			    dipole_1.bdpl[k - 1] = -dp;
			    dipole_1.sdpl[k - 1] = 1. - ps;
			}
			io___37.ciunit = iounit_1.iout;
			s_wsfe(&io___37);
			do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&dp, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&ps, (ftnlen)sizeof(doublereal))
				;
			e_wsfe();
			goto L140;
		    }
		}
	    }
L140:
	    ;
	}
    }

/*     remove zero bond dipoles from the list of dipoles */

    dipole_1.ndipole = 0;
    i__1 = bond_1.nbond;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (dipole_1.bdpl[i__ - 1] != 0.) {
	    ++dipole_1.ndipole;
	    idpl_ref(1, dipole_1.ndipole) = idpl_ref(1, i__);
	    idpl_ref(2, dipole_1.ndipole) = idpl_ref(2, i__);
	    dipole_1.bdpl[dipole_1.ndipole - 1] = dipole_1.bdpl[i__ - 1];
	    dipole_1.sdpl[dipole_1.ndipole - 1] = dipole_1.sdpl[i__ - 1];
	}
    }

/*     turn off dipole-dipole and charge-dipole terms if not used */

    if (dipole_1.ndipole == 0) {
	potent_1.use_dipole__ = FALSE_;
	potent_1.use_chgdpl__ = FALSE_;
    }
    return 0;
} /* kdipole_ */

#undef keyline_ref
#undef bndlist_ref
#undef idpl_ref
#undef ibnd_ref
#undef kd5_ref
#undef kd4_ref
#undef kd3_ref
#undef kd_ref
#undef i12_ref


