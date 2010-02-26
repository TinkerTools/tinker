/* kstrtor.f -- translated by f2c (version 20050501).
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
    doublereal btcon[1500]	/* was [3][500] */;
    char kbt[8000];
} ksttor_;

#define ksttor_1 ksttor_

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
    doublereal kst[300000]	/* was [3][100000] */;
    integer nstrtor, ist[200000]	/* was [2][100000] */;
} strtor_;

#define strtor_1 strtor_

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
static integer c__4 = 4;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine kstrtor  --  find stretch-torsion parameters  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "kstrtor" assigns stretch-torsion parameters to torsions */
/*     needing them, and processes any new or changed values */


/* Subroutine */ int kstrtor_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Additional Stretch-Torsion Parameters "
	    ":\002,//,5x,\002Atom Classes\002,7x,\0021-Fold\002,6x,\0022-Fol"
	    "d\002,6x,\0023-Fold\002,/)";
    static char fmt_30[] = "(1x,4i4,1x,3f12.3)";
    static char fmt_40[] = "(/,\002 KSTRTOR  --  Too many Stretch-Torsion"
	    "\002,\002 Parameters\002)";

    /* System generated locals */
    address a__1[4], a__2[3];
    integer i__1, i__2, i__3[4], i__4[3], i__5;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i__, j, k, ia, ib, ic, id;
    static char pa[4], pb[4], pc[4], pd[4], pt[16];
    static doublereal bt1, bt2, bt3;
    static char pt0[16];
    static integer ita, itb, itc, itd, nbt, size, next;
    static char blank[16], zeros[4];
    static logical header;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    ;
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___16 = { 1, string, 1, 0, 120, 1 };
    static cilist io___17 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_40, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define kbt_ref(a_0,a_1) &ksttor_1.kbt[(a_1)*16 + a_0 - 16]
#define ist_ref(a_1,a_2) strtor_1.ist[(a_2)*2 + a_1 - 3]
#define kst_ref(a_1,a_2) strtor_1.kst[(a_2)*3 + a_1 - 4]
#define btcon_ref(a_1,a_2) ksttor_1.btcon[(a_2)*3 + a_1 - 4]
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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  ksttor.i  --  forcefield parameters for stretch-torsions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxnbt   maximum number of stretch-torsion parameter entries */

/*     btcon    force constant parameters for stretch-torsion */
/*     kbt      string of atom classes for stretch-torsion terms */




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
/*     ##  strtor.i  --  stretch-torsions in the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     kst       1-, 2- and 3-fold stretch-torsion force constants */
/*     nstrtor   total number of stretch-torsion interactions */
/*     ist       torsion and bond numbers used in stretch-torsion */




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




/*     process keywords containing stretch-torsion parameters */

    s_copy(blank, "                ", (ftnlen)16, (ftnlen)16);
    s_copy(zeros, "0000", (ftnlen)4, (ftnlen)4);
    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "STRTORS ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    bt1 = 0.;
	    bt2 = 0.;
	    bt3 = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___16);
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
	    i__2 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&bt1, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&bt2, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&bt3, (ftnlen)sizeof(
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
		io___17.ciunit = iounit_1.iout;
		s_wsfe(&io___17);
		e_wsfe();
	    }
	    io___18.ciunit = iounit_1.iout;
	    s_wsfe(&io___18);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&bt1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&bt2, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&bt3, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    numeral_(&itc, pc, &size, (ftnlen)4);
	    numeral_(&itd, pd, &size, (ftnlen)4);
	    if (itb < itc) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (itc < itb) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (ita <= itd) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (itd < ita) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    }
	    for (j = 1; j <= 500; ++j) {
		if (s_cmp(kbt_ref(0, j), blank, (ftnlen)16, (ftnlen)16) == 0 
			|| s_cmp(kbt_ref(0, j), pt, (ftnlen)16, (ftnlen)16) ==
			 0) {
		    s_copy(kbt_ref(0, j), pt, (ftnlen)16, (ftnlen)16);
		    btcon_ref(1, j) = bt1;
		    btcon_ref(2, j) = bt2;
		    btcon_ref(3, j) = bt3;
		    goto L50;
		}
	    }
	    io___30.ciunit = iounit_1.iout;
	    s_wsfe(&io___30);
	    e_wsfe();
	    inform_1.abort = TRUE_;
L50:
	    ;
	}
    }

/*     determine the total number of forcefield parameters */

    nbt = 500;
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(kbt_ref(0, i__), blank, (ftnlen)16, (ftnlen)16) == 0) {
	    nbt = i__ - 1;
	}
    }

/*     assign the stretch-torsion parameters for each torsion */

    strtor_1.nstrtor = 0;
    if (nbt != 0) {
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
	    if (itb < itc) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (itc < itb) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (ita <= itd) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (itd < ita) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    }
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = zeros;
	    i__4[1] = 8, a__2[1] = pt + 4;
	    i__4[2] = 4, a__2[2] = zeros;
	    s_cat(pt0, a__2, i__4, &c__3, (ftnlen)16);
	    i__2 = nbt;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kbt_ref(0, j), pt, (ftnlen)16, (ftnlen)16) == 0) {
		    ++strtor_1.nstrtor;
		    kst_ref(1, strtor_1.nstrtor) = btcon_ref(1, j);
		    kst_ref(2, strtor_1.nstrtor) = btcon_ref(2, j);
		    kst_ref(3, strtor_1.nstrtor) = btcon_ref(3, j);
		    ist_ref(1, strtor_1.nstrtor) = i__;
		    i__5 = couple_1.n12[ib - 1];
		    for (k = 1; k <= i__5; ++k) {
			if (i12_ref(k, ib) == ic) {
			    ist_ref(2, strtor_1.nstrtor) = bndlist_ref(k, ib);
			    goto L60;
			}
		    }
		}
	    }
	    i__2 = nbt;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kbt_ref(0, j), pt0, (ftnlen)16, (ftnlen)16) == 0) {
		    ++strtor_1.nstrtor;
		    kst_ref(1, strtor_1.nstrtor) = btcon_ref(1, j);
		    kst_ref(2, strtor_1.nstrtor) = btcon_ref(2, j);
		    kst_ref(3, strtor_1.nstrtor) = btcon_ref(3, j);
		    ist_ref(1, strtor_1.nstrtor) = i__;
		    i__5 = couple_1.n12[ib - 1];
		    for (k = 1; k <= i__5; ++k) {
			if (i12_ref(k, ib) == ic) {
			    ist_ref(2, strtor_1.nstrtor) = bndlist_ref(k, ib);
			    goto L60;
			}
		    }
		}
	    }
L60:
	    ;
	}
    }

/*     turn off the stretch-torsion potential if it is not used */

    if (strtor_1.nstrtor == 0) {
	potent_1.use_strtor__ = FALSE_;
    }
    return 0;
} /* kstrtor_ */

#undef keyline_ref
#undef bndlist_ref
#undef itors_ref
#undef btcon_ref
#undef kst_ref
#undef ist_ref
#undef kbt_ref
#undef i12_ref


