/* kbond.f -- translated by f2c (version 20050501).
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
    doublereal bcon[2000], blen[2000], bcon5[500], blen5[500], bcon4[500], 
	    blen4[500], bcon3[500], blen3[500], dlen[500];
    char kb[16000], kb5[4000], kb4[4000], kb3[4000], kel[6000];
} kbonds_;

#define kbonds_1 kbonds_

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
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine kbond  --  bond stretch parameter assignment  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "kbond" assigns a force constant and ideal bond length */
/*     to each bond in the structure and processes any new or */
/*     changed parameter values */


/* Subroutine */ int kbond_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Additional Bond Stretching Parameters "
	    ":\002,//,5x,\002Atom Classes\002,9x,\002K(S)\002,6x,\002Lengt"
	    "h\002,/)";
    static char fmt_30[] = "(6x,2i4,4x,f12.3,f12.4)";
    static char fmt_40[] = "(6x,2i4,4x,f12.3,f12.4,3x,a6)";
    static char fmt_50[] = "(/,\002 KBOND  --  Too many Bond Stretching\002"
	    ",\002 Parameters\002)";
    static char fmt_70[] = "(/,\002 KBOND  --  Too many 5-Ring Stretching"
	    "\002,\002 Parameters\002)";
    static char fmt_90[] = "(/,\002 KBOND  --  Too many 4-Ring Stretching"
	    "\002,\002 Parameters\002)";
    static char fmt_110[] = "(/,\002 KBOND  --  Too many 3-Ring Stretchin"
	    "g\002,\002 Parameters\002)";
    static char fmt_140[] = "(/,\002 Undefined Bond Stretching Parameters "
	    ":\002,//,\002 Type\002,13x,\002Atom Names\002,11x,\002Atom Class"
	    "es\002,/)";
    static char fmt_150[] = "(1x,a6,5x,i6,\002-\002,a3,i6,\002-\002,a3,7x,2i"
	    "5)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2], i__4;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static logical use_ring__;
    static integer i__, j;
    static doublereal bd, fc;
    static integer ia, ib, nb;
    static char pa[4], pb[4], pt[8];
    static integer nb4, nb5, nb3, ita, itb;
    static logical done;
    static integer size, next;
    static char label[6], blank[8];
    extern /* Subroutine */ int keneg_(void);
    static integer minat, iring;
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
    static cilist io___24 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_150, 0 };



#define kb_ref(a_0,a_1) &kbonds_1.kb[(a_1)*8 + a_0 - 8]
#define kb3_ref(a_0,a_1) &kbonds_1.kb3[(a_1)*8 + a_0 - 8]
#define kb4_ref(a_0,a_1) &kbonds_1.kb4[(a_1)*8 + a_0 - 8]
#define kb5_ref(a_0,a_1) &kbonds_1.kb5[(a_1)*8 + a_0 - 8]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  kbonds.i  --  forcefield parameters for bond stretching  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxnb   maximum number of bond stretch parameter entries */
/*     maxnb5  maximum number of 5-membered ring bond stretch entries */
/*     maxnb4  maximum number of 4-membered ring bond stretch entries */
/*     maxnb3  maximum number of 3-membered ring bond stretch entries */
/*     maxnel  maximum number of electronegativity bond corrections */

/*     bcon    force constant parameters for harmonic bond stretch */
/*     blen    bond length parameters for harmonic bond stretch */
/*     bcon5   force constant parameters for 5-ring bond stretch */
/*     blen5   bond length parameters for 5-ring bond stretch */
/*     bcon4   force constant parameters for 4-ring bond stretch */
/*     blen4   bond length parameters for 4-ring bond stretch */
/*     bcon3   force constant parameters for 3-ring bond stretch */
/*     blen3   bond length parameters for 3-ring bond stretch */
/*     dlen    electronegativity bond length correction parameters */
/*     kb      string of atom classes for harmonic bond stretch */
/*     kb5     string of atom classes for 5-ring bond stretch */
/*     kb4     string of atom classes for 4-ring bond stretch */
/*     kb3     string of atom classes for 3-ring bond stretch */
/*     kel     string of atom classes for electronegativity corrections */




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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     process keywords containing bond stretch parameters */

    s_copy(blank, "        ", (ftnlen)8, (ftnlen)8);
    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	iring = -1;
	if (s_cmp(keyword, "BOND ", (ftnlen)5, (ftnlen)5) == 0) {
	    iring = 0;
	}
	if (s_cmp(keyword, "BOND5 ", (ftnlen)6, (ftnlen)6) == 0) {
	    iring = 5;
	}
	if (s_cmp(keyword, "BOND4 ", (ftnlen)6, (ftnlen)6) == 0) {
	    iring = 4;
	}
	if (s_cmp(keyword, "BOND3 ", (ftnlen)6, (ftnlen)6) == 0) {
	    iring = 3;
	}
	if (iring >= 0) {
	    ia = 0;
	    ib = 0;
	    fc = 0.;
	    bd = 0.;
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
	    i__2 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&bd, (ftnlen)sizeof(
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
		io___14.ciunit = iounit_1.iout;
		s_wsfe(&io___14);
		e_wsfe();
	    }
	    if (iring == 0) {
		io___15.ciunit = iounit_1.iout;
		s_wsfe(&io___15);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bd, (ftnlen)sizeof(doublereal));
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
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bd, (ftnlen)sizeof(doublereal));
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
		for (j = 1; j <= 2000; ++j) {
		    if (s_cmp(kb_ref(0, j), blank, (ftnlen)8, (ftnlen)8) == 0 
			    || s_cmp(kb_ref(0, j), pt, (ftnlen)8, (ftnlen)8) 
			    == 0) {
			s_copy(kb_ref(0, j), pt, (ftnlen)8, (ftnlen)8);
			kbonds_1.bcon[j - 1] = fc;
			kbonds_1.blen[j - 1] = bd;
			goto L60;
		    }
		}
		io___23.ciunit = iounit_1.iout;
		s_wsfe(&io___23);
		e_wsfe();
		inform_1.abort = TRUE_;
L60:
		;
	    } else if (iring == 5) {
		for (j = 1; j <= 500; ++j) {
		    if (s_cmp(kb5_ref(0, j), blank, (ftnlen)8, (ftnlen)8) == 
			    0 || s_cmp(kb5_ref(0, j), pt, (ftnlen)8, (ftnlen)
			    8) == 0) {
			s_copy(kb5_ref(0, j), pt, (ftnlen)8, (ftnlen)8);
			kbonds_1.bcon5[j - 1] = fc;
			kbonds_1.blen5[j - 1] = bd;
			goto L80;
		    }
		}
		io___24.ciunit = iounit_1.iout;
		s_wsfe(&io___24);
		e_wsfe();
		inform_1.abort = TRUE_;
L80:
		;
	    } else if (iring == 4) {
		for (j = 1; j <= 500; ++j) {
		    if (s_cmp(kb4_ref(0, j), blank, (ftnlen)8, (ftnlen)8) == 
			    0 || s_cmp(kb4_ref(0, j), pt, (ftnlen)8, (ftnlen)
			    8) == 0) {
			s_copy(kb4_ref(0, j), pt, (ftnlen)8, (ftnlen)8);
			kbonds_1.bcon4[j - 1] = fc;
			kbonds_1.blen4[j - 1] = bd;
			goto L100;
		    }
		}
		io___25.ciunit = iounit_1.iout;
		s_wsfe(&io___25);
		e_wsfe();
		inform_1.abort = TRUE_;
L100:
		;
	    } else if (iring == 3) {
		for (j = 1; j <= 500; ++j) {
		    if (s_cmp(kb3_ref(0, j), blank, (ftnlen)8, (ftnlen)8) == 
			    0 || s_cmp(kb3_ref(0, j), pt, (ftnlen)8, (ftnlen)
			    8) == 0) {
			s_copy(kb3_ref(0, j), pt, (ftnlen)8, (ftnlen)8);
			kbonds_1.bcon3[j - 1] = fc;
			kbonds_1.blen3[j - 1] = bd;
			goto L120;
		    }
		}
		io___26.ciunit = iounit_1.iout;
		s_wsfe(&io___26);
		e_wsfe();
		inform_1.abort = TRUE_;
L120:
		;
	    }
	}
    }

/*     determine the total number of forcefield parameters */

    nb = 2000;
    nb5 = 500;
    nb4 = 500;
    nb3 = 500;
    for (i__ = 2000; i__ >= 1; --i__) {
	if (s_cmp(kb_ref(0, i__), blank, (ftnlen)8, (ftnlen)8) == 0) {
	    nb = i__ - 1;
	}
    }
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(kb5_ref(0, i__), blank, (ftnlen)8, (ftnlen)8) == 0) {
	    nb5 = i__ - 1;
	}
    }
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(kb4_ref(0, i__), blank, (ftnlen)8, (ftnlen)8) == 0) {
	    nb4 = i__ - 1;
	}
    }
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(kb3_ref(0, i__), blank, (ftnlen)8, (ftnlen)8) == 0) {
	    nb3 = i__ - 1;
	}
    }
    use_ring__ = FALSE_;
/* Computing MIN */
    i__1 = min(nb5,nb4);
    if (min(i__1,nb3) != 0) {
	use_ring__ = TRUE_;
    }

/*     assign ideal bond length and force constant for each bond */

    header = TRUE_;
    i__1 = bond_1.nbond;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ibnd_ref(1, i__);
	ib = ibnd_ref(2, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
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
	bond_1.bk[i__ - 1] = 0.;
	bond_1.bl[i__ - 1] = 0.;
	done = FALSE_;

/*     make a check for bonds contained inside small rings */

	iring = 0;
	if (use_ring__) {
	    chkring_(&iring, &ia, &ib, &c__0, &c__0);
	    if (iring == 6) {
		iring = 0;
	    }
	    if (iring == 5 && nb5 == 0) {
		iring = 0;
	    }
	    if (iring == 4 && nb4 == 0) {
		iring = 0;
	    }
	    if (iring == 3 && nb3 == 0) {
		iring = 0;
	    }
	}

/*     assign bond stretching parameters for each bond */

	if (iring == 0) {
	    i__2 = nb;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kb_ref(0, j), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    bond_1.bk[i__ - 1] = kbonds_1.bcon[j - 1];
		    bond_1.bl[i__ - 1] = kbonds_1.blen[j - 1];
		    done = TRUE_;
		    goto L130;
		}
	    }

/*     assign stretching parameters for 5-membered ring bonds */

	} else if (iring == 5) {
	    i__2 = nb5;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kb5_ref(0, j), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    bond_1.bk[i__ - 1] = kbonds_1.bcon5[j - 1];
		    bond_1.bl[i__ - 1] = kbonds_1.blen5[j - 1];
		    done = TRUE_;
		    goto L130;
		}
	    }

/*     assign stretching parameters for 4-membered ring bonds */

	} else if (iring == 4) {
	    i__2 = nb4;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kb4_ref(0, j), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    bond_1.bk[i__ - 1] = kbonds_1.bcon4[j - 1];
		    bond_1.bl[i__ - 1] = kbonds_1.blen4[j - 1];
		    done = TRUE_;
		    goto L130;
		}
	    }

/*     assign stretching parameters for 3-membered ring bonds */

	} else if (iring == 3) {
	    i__2 = nb3;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kb3_ref(0, j), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    bond_1.bk[i__ - 1] = kbonds_1.bcon3[j - 1];
		    bond_1.bl[i__ - 1] = kbonds_1.blen3[j - 1];
		    done = TRUE_;
		    goto L130;
		}
	    }
	}

/*     warning if suitable bond stretching parameter not found */

L130:
/* Computing MIN */
	i__2 = atmtyp_1.atomic[ia - 1], i__4 = atmtyp_1.atomic[ib - 1];
	minat = min(i__2,i__4);
	if (minat == 0) {
	    done = TRUE_;
	}
	if (potent_1.use_bond__ && ! done) {
	    if (usage_1.use[ia - 1] || usage_1.use[ib - 1]) {
		inform_1.abort = TRUE_;
	    }
	    if (header) {
		header = FALSE_;
		io___36.ciunit = iounit_1.iout;
		s_wsfe(&io___36);
		e_wsfe();
	    }
	    s_copy(label, "Bond  ", (ftnlen)6, (ftnlen)6);
	    if (iring == 5) {
		s_copy(label, "5-Ring", (ftnlen)6, (ftnlen)6);
	    }
	    if (iring == 4) {
		s_copy(label, "4-Ring", (ftnlen)6, (ftnlen)6);
	    }
	    if (iring == 3) {
		s_copy(label, "3-Ring", (ftnlen)6, (ftnlen)6);
	    }
	    io___37.ciunit = iounit_1.iout;
	    s_wsfe(&io___37);
	    do_fio(&c__1, label, (ftnlen)6);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ia), (ftnlen)3);
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ib), (ftnlen)3);
	    do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }

/*     check for electronegativity bond length corrections */

    keneg_();

/*     turn off the bond stretch potential if it is not used */

    if (bond_1.nbond == 0) {
	potent_1.use_bond__ = FALSE_;
    }
    return 0;
} /* kbond_ */

#undef keyline_ref
#undef name___ref
#undef ibnd_ref
#undef kb5_ref
#undef kb4_ref
#undef kb3_ref
#undef kb_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine keneg  --  assign electronegativity parameters  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "keneg" applies primary and secondary electronegativity bond */
/*     length corrections to applicable bond parameters */

/*     note this version does not scale multiple corrections to the */
/*     same bond by increasing powers of 0.62 as in MM3 */


/* Subroutine */ int keneg_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Additional Electronegativity Parameter"
	    "s :\002,//,5x,\002Atom Classes\002,18x,\002dLength\002,/)";
    static char fmt_30[] = "(4x,3i4,14x,f12.4)";
    static char fmt_40[] = "(/,\002 KENEG  --  Too many Electronegativity"
	    "\002,\002 Parameters\002)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2, i__3[3], i__4;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i__, j, k, m, ia, ib, ic, id;
    static doublereal dl;
    static char pa[4], pb[4], pc[4], pd[4], pt[12], pt1[12], pt2[12];
    static integer ita, nel, itb, itc, itd, size, next;
    static char blank[12];
    static logical header;
    static doublereal factor;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    ;
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___49 = { 1, string, 1, 0, 120, 1 };
    static cilist io___50 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___58 = { 0, 0, 0, fmt_40, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define kel_ref(a_0,a_1) &kbonds_1.kel[(a_1)*12 + a_0 - 12]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define itors_ref(a_1,a_2) tors_1.itors[(a_2)*4 + a_1 - 5]
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
/*     ##  kbonds.i  --  forcefield parameters for bond stretching  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxnb   maximum number of bond stretch parameter entries */
/*     maxnb5  maximum number of 5-membered ring bond stretch entries */
/*     maxnb4  maximum number of 4-membered ring bond stretch entries */
/*     maxnb3  maximum number of 3-membered ring bond stretch entries */
/*     maxnel  maximum number of electronegativity bond corrections */

/*     bcon    force constant parameters for harmonic bond stretch */
/*     blen    bond length parameters for harmonic bond stretch */
/*     bcon5   force constant parameters for 5-ring bond stretch */
/*     blen5   bond length parameters for 5-ring bond stretch */
/*     bcon4   force constant parameters for 4-ring bond stretch */
/*     blen4   bond length parameters for 4-ring bond stretch */
/*     bcon3   force constant parameters for 3-ring bond stretch */
/*     blen3   bond length parameters for 3-ring bond stretch */
/*     dlen    electronegativity bond length correction parameters */
/*     kb      string of atom classes for harmonic bond stretch */
/*     kb5     string of atom classes for 5-ring bond stretch */
/*     kb4     string of atom classes for 4-ring bond stretch */
/*     kb3     string of atom classes for 3-ring bond stretch */
/*     kel     string of atom classes for electronegativity corrections */




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




/*     process keywords containing electronegativity parameters */

    s_copy(blank, "            ", (ftnlen)12, (ftnlen)12);
    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "ELECTNEG ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    dl = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___49);
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
	    i__2 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&dl, (ftnlen)sizeof(
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
		io___50.ciunit = iounit_1.iout;
		s_wsfe(&io___50);
		e_wsfe();
	    }
	    io___51.ciunit = iounit_1.iout;
	    s_wsfe(&io___51);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dl, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    size = 4;
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
/* Writing concatenation */
	    i__3[0] = 4, a__1[0] = pa;
	    i__3[1] = 4, a__1[1] = pb;
	    i__3[2] = 4, a__1[2] = pc;
	    s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	    for (j = 1; j <= 500; ++j) {
		if (s_cmp(kel_ref(0, j), blank, (ftnlen)12, (ftnlen)12) == 0 
			|| s_cmp(kel_ref(0, j), pt, (ftnlen)12, (ftnlen)12) ==
			 0) {
		    s_copy(kel_ref(0, j), pt, (ftnlen)12, (ftnlen)12);
		    kbonds_1.dlen[j - 1] = dl;
		    goto L50;
		}
	    }
	    io___58.ciunit = iounit_1.iout;
	    s_wsfe(&io___58);
	    e_wsfe();
	    inform_1.abort = TRUE_;
L50:
	    ;
	}
    }

/*     determine the total number of forcefield parameters */

    nel = 500;
    for (i__ = 1; i__ <= 500; ++i__) {
	if (s_cmp(kel_ref(0, i__), blank, (ftnlen)12, (ftnlen)12) == 0) {
	    nel = i__ - 1;
	    goto L60;
	}
    }
L60:

/*     check angles for primary electronegativity corrections */

    if (nel != 0) {
	i__1 = angle_1.nangle;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = iang_ref(1, i__);
	    ib = iang_ref(2, i__);
	    ic = iang_ref(3, i__);
	    ita = atmtyp_1.class__[ia - 1];
	    itb = atmtyp_1.class__[ib - 1];
	    itc = atmtyp_1.class__[ic - 1];
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    numeral_(&itc, pc, &size, (ftnlen)4);
/* Writing concatenation */
	    i__3[0] = 4, a__1[0] = pa;
	    i__3[1] = 4, a__1[1] = pb;
	    i__3[2] = 4, a__1[2] = pc;
	    s_cat(pt1, a__1, i__3, &c__3, (ftnlen)12);
/* Writing concatenation */
	    i__3[0] = 4, a__1[0] = pc;
	    i__3[1] = 4, a__1[1] = pb;
	    i__3[2] = 4, a__1[2] = pa;
	    s_cat(pt2, a__1, i__3, &c__3, (ftnlen)12);

/*     search the parameter set for a match to either bond */

	    i__2 = nel;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kel_ref(0, j), pt1, (ftnlen)12, (ftnlen)12) == 0) {
		    i__4 = couple_1.n12[ia - 1];
		    for (k = 1; k <= i__4; ++k) {
			if (i12_ref(k, ia) == ib) {
			    m = bndlist_ref(k, ia);
			    bond_1.bl[m - 1] += kbonds_1.dlen[j - 1];
			}
		    }
		    goto L70;
		} else if (s_cmp(kel_ref(0, j), pt2, (ftnlen)12, (ftnlen)12) 
			== 0) {
		    i__4 = couple_1.n12[ic - 1];
		    for (k = 1; k <= i__4; ++k) {
			if (i12_ref(k, ic) == ib) {
			    m = bndlist_ref(k, ic);
			    bond_1.bl[m - 1] += kbonds_1.dlen[j - 1];
			}
		    }
		    goto L70;
		}
	    }
L70:
	    ;
	}

/*     check torsions for secondary electronegativity corrections */

	factor = .4;
	i__1 = tors_1.ntors;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = itors_ref(1, i__);
	    ib = itors_ref(2, i__);
	    ic = itors_ref(3, i__);
	    id = itors_ref(4, i__);
	    ita = atmtyp_1.class__[ia - 1];
	    itb = atmtyp_1.class__[ib - 1];
	    itc = atmtyp_1.class__[ic - 1];
	    itd = atmtyp_1.class__[id - 1];
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    numeral_(&itc, pc, &size, (ftnlen)4);
	    numeral_(&itd, pd, &size, (ftnlen)4);
/* Writing concatenation */
	    i__3[0] = 4, a__1[0] = pa;
	    i__3[1] = 4, a__1[1] = pb;
	    i__3[2] = 4, a__1[2] = pd;
	    s_cat(pt1, a__1, i__3, &c__3, (ftnlen)12);
/* Writing concatenation */
	    i__3[0] = 4, a__1[0] = pd;
	    i__3[1] = 4, a__1[1] = pc;
	    i__3[2] = 4, a__1[2] = pa;
	    s_cat(pt2, a__1, i__3, &c__3, (ftnlen)12);

/*     turn off electronegativity effect for attached hydrogens */

	    if (atmtyp_1.atomic[id - 1] <= 1) {
		s_copy(pt1, blank, (ftnlen)12, (ftnlen)12);
	    }
	    if (atmtyp_1.atomic[ia - 1] <= 1) {
		s_copy(pt2, blank, (ftnlen)12, (ftnlen)12);
	    }

/*     search the parameter set for a match to either bond */

	    i__2 = nel;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kel_ref(0, j), pt1, (ftnlen)12, (ftnlen)12) == 0) {
		    i__4 = couple_1.n12[ia - 1];
		    for (k = 1; k <= i__4; ++k) {
			if (i12_ref(k, ia) == ib) {
			    m = bndlist_ref(k, ia);
			    bond_1.bl[m - 1] += factor * kbonds_1.dlen[j - 1];
			}
		    }
		    goto L80;
		} else if (s_cmp(kel_ref(0, j), pt2, (ftnlen)12, (ftnlen)12) 
			== 0) {
		    i__4 = couple_1.n12[id - 1];
		    for (k = 1; k <= i__4; ++k) {
			if (i12_ref(k, id) == ic) {
			    m = bndlist_ref(k, id);
			    bond_1.bl[m - 1] += factor * kbonds_1.dlen[j - 1];
			}
		    }
		    goto L80;
		}
	    }
L80:
	    ;
	}
    }
    return 0;
} /* keneg_ */

#undef keyline_ref
#undef bndlist_ref
#undef itors_ref
#undef iang_ref
#undef kel_ref
#undef i12_ref


