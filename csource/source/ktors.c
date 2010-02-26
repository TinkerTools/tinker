/* ktors.f -- translated by f2c (version 20050501).
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
    doublereal t1[4000]	/* was [2][2000] */, t2[4000]	/* was [2][2000] */, 
	    t3[4000]	/* was [2][2000] */, t4[4000]	/* was [2][2000] */, 
	    t5[4000]	/* was [2][2000] */, t6[4000]	/* was [2][2000] */, 
	    t15[1000]	/* was [2][500] */, t25[1000]	/* was [2][500] */, 
	    t35[1000]	/* was [2][500] */, t45[1000]	/* was [2][500] */, 
	    t55[1000]	/* was [2][500] */, t65[1000]	/* was [2][500] */, 
	    t14[1000]	/* was [2][500] */, t24[1000]	/* was [2][500] */, 
	    t34[1000]	/* was [2][500] */, t44[1000]	/* was [2][500] */, 
	    t54[1000]	/* was [2][500] */, t64[1000]	/* was [2][500] */;
    char kt[32000], kt5[8000], kt4[8000];
} ktorsn_;

#define ktorsn_1 ktorsn_

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
    doublereal tors1[400000]	/* was [4][100000] */, tors2[400000]	/* 
	    was [4][100000] */, tors3[400000]	/* was [4][100000] */, tors4[
	    400000]	/* was [4][100000] */, tors5[400000]	/* was [4][
	    100000] */, tors6[400000]	/* was [4][100000] */;
    integer ntors, itors[400000]	/* was [4][100000] */;
} tors_;

#define tors_1 tors_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__4 = 4;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine ktors  --  torsional parameter assignment  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "ktors" assigns torsional parameters to each torsion in */
/*     the structure and processes any new or changed values */


/* Subroutine */ int ktors_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Additional Torsional Parameters :\002,//"
	    ",5x,\002Atom Classes\002,4x,\0021-Fold\002,4x,\0022-Fold\002,4x"
	    ",\0023-Fold\002,4x,\0024-Fold\002,4x,\0025-Fold\002,4x,\0026-Fold"
	    "\002,/)";
    static char fmt_30[] = "(1x,4i4,1x,6(f6.2,i4))";
    static char fmt_40[] = "(1x,4i4,1x,6(f6.2,i4),3x,a6)";
    static char fmt_50[] = "(/,\002 KTORS  --  Too many Torsional Angle\002"
	    ",\002 Parameters\002)";
    static char fmt_70[] = "(/,\002 KTORS  --  Too many 5-Ring Torsiona"
	    "l\002,\002 Parameters\002)";
    static char fmt_90[] = "(/,\002 KTORS  --  Too many 4-Ring Torsiona"
	    "l\002,\002 Parameters\002)";
    static char fmt_120[] = "(/,\002 Undefined Torsional Parameters :\002,"
	    "//,\002 Type\002,24x,\002Atom Names\002,24x,\002Atom Classes\002"
	    ",/)";
    static char fmt_130[] = "(1x,a7,4x,4(i6,\002-\002,a3),5x,4i5)";

    /* System generated locals */
    address a__1[4], a__2[2], a__3[3];
    integer i__1, i__2, i__3[4], i__4[2], i__5[3], i__6;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen),
	     i_dnnt(doublereal *);
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static logical use_ring__;
    extern /* Subroutine */ int torphase_(integer *, doublereal *, doublereal 
	    *);
    static integer i__, j, ia, ib, ic, id;
    static char pa[4], pb[4], pc[4];
    static integer nt, ft[6];
    static doublereal st[6], vt[6];
    static char pd[4], pt[16], pt0[16];
    static integer nt4, nt5;
    static char pt1[16], pt2[16];
    static integer ita, itb, itc, itd;
    static logical done;
    static integer size, next;
    static char label[7];
    static doublereal angle;
    static char blank[16];
    static integer iring, minat, ilist;
    static char klist[16*2000];
    static integer nlist;
    static char zeros[4];
    static logical header;
    static char record[120];
    static integer kindex[2000];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int chkring_(integer *, integer *, integer *, 
	    integer *, integer *), numeral_(integer *, char *, integer *, 
	    ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___18 = { 1, string, 1, 0, 120, 1 };
    static cilist io___25 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___49 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_130, 0 };



#define t1_ref(a_1,a_2) ktorsn_1.t1[(a_2)*2 + a_1 - 3]
#define t2_ref(a_1,a_2) ktorsn_1.t2[(a_2)*2 + a_1 - 3]
#define t3_ref(a_1,a_2) ktorsn_1.t3[(a_2)*2 + a_1 - 3]
#define t4_ref(a_1,a_2) ktorsn_1.t4[(a_2)*2 + a_1 - 3]
#define t5_ref(a_1,a_2) ktorsn_1.t5[(a_2)*2 + a_1 - 3]
#define t6_ref(a_1,a_2) ktorsn_1.t6[(a_2)*2 + a_1 - 3]
#define t14_ref(a_1,a_2) ktorsn_1.t14[(a_2)*2 + a_1 - 3]
#define t15_ref(a_1,a_2) ktorsn_1.t15[(a_2)*2 + a_1 - 3]
#define t25_ref(a_1,a_2) ktorsn_1.t25[(a_2)*2 + a_1 - 3]
#define t35_ref(a_1,a_2) ktorsn_1.t35[(a_2)*2 + a_1 - 3]
#define t45_ref(a_1,a_2) ktorsn_1.t45[(a_2)*2 + a_1 - 3]
#define t55_ref(a_1,a_2) ktorsn_1.t55[(a_2)*2 + a_1 - 3]
#define t65_ref(a_1,a_2) ktorsn_1.t65[(a_2)*2 + a_1 - 3]
#define t24_ref(a_1,a_2) ktorsn_1.t24[(a_2)*2 + a_1 - 3]
#define t34_ref(a_1,a_2) ktorsn_1.t34[(a_2)*2 + a_1 - 3]
#define t44_ref(a_1,a_2) ktorsn_1.t44[(a_2)*2 + a_1 - 3]
#define t54_ref(a_1,a_2) ktorsn_1.t54[(a_2)*2 + a_1 - 3]
#define t64_ref(a_1,a_2) ktorsn_1.t64[(a_2)*2 + a_1 - 3]
#define kt_ref(a_0,a_1) &ktorsn_1.kt[(a_1)*16 + a_0 - 16]
#define kt4_ref(a_0,a_1) &ktorsn_1.kt4[(a_1)*16 + a_0 - 16]
#define kt5_ref(a_0,a_1) &ktorsn_1.kt5[(a_1)*16 + a_0 - 16]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define tors1_ref(a_1,a_2) tors_1.tors1[(a_2)*4 + a_1 - 5]
#define tors2_ref(a_1,a_2) tors_1.tors2[(a_2)*4 + a_1 - 5]
#define tors3_ref(a_1,a_2) tors_1.tors3[(a_2)*4 + a_1 - 5]
#define tors4_ref(a_1,a_2) tors_1.tors4[(a_2)*4 + a_1 - 5]
#define tors5_ref(a_1,a_2) tors_1.tors5[(a_2)*4 + a_1 - 5]
#define tors6_ref(a_1,a_2) tors_1.tors6[(a_2)*4 + a_1 - 5]
#define klist_ref(a_0,a_1) &klist[(a_1)*16 + a_0 - 16]
#define itors_ref(a_1,a_2) tors_1.itors[(a_2)*4 + a_1 - 5]
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
/*     ##  ktorsn.i  --  forcefield parameters for torsional angles  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxnt    maximum number of torsional angle parameter entries */
/*     maxnt5   maximum number of 5-membered ring torsion entries */
/*     maxnt4   maximum number of 4-membered ring torsion entries */

/*     t1       torsional parameters for standard 1-fold rotation */
/*     t2       torsional parameters for standard 2-fold rotation */
/*     t3       torsional parameters for standard 3-fold rotation */
/*     t4       torsional parameters for standard 4-fold rotation */
/*     t5       torsional parameters for standard 5-fold rotation */
/*     t6       torsional parameters for standard 6-fold rotation */
/*     t15      torsional parameters for 1-fold rotation in 5-ring */
/*     t25      torsional parameters for 2-fold rotation in 5-ring */
/*     t35      torsional parameters for 3-fold rotation in 5-ring */
/*     t45      torsional parameters for 4-fold rotation in 5-ring */
/*     t55      torsional parameters for 5-fold rotation in 5-ring */
/*     t65      torsional parameters for 6-fold rotation in 5-ring */
/*     t14      torsional parameters for 1-fold rotation in 4-ring */
/*     t24      torsional parameters for 2-fold rotation in 4-ring */
/*     t34      torsional parameters for 3-fold rotation in 4-ring */
/*     t44      torsional parameters for 4-fold rotation in 4-ring */
/*     t54      torsional parameters for 5-fold rotation in 4-ring */
/*     t64      torsional parameters for 6-fold rotation in 4-ring */
/*     kt       string of atom classes for torsional angles */
/*     kt5      string of atom classes for 5-ring torsions */
/*     kt4      string of atom classes for 4-ring torsions */




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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     process keywords containing torsional angle parameters */

    s_copy(blank, "                ", (ftnlen)16, (ftnlen)16);
    s_copy(zeros, "0000", (ftnlen)4, (ftnlen)4);
    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	iring = -1;
	if (s_cmp(keyword, "TORSION ", (ftnlen)8, (ftnlen)8) == 0) {
	    iring = 0;
	}
	if (s_cmp(keyword, "TORSION5 ", (ftnlen)9, (ftnlen)9) == 0) {
	    iring = 5;
	}
	if (s_cmp(keyword, "TORSION4 ", (ftnlen)9, (ftnlen)9) == 0) {
	    iring = 4;
	}
	if (iring >= 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    for (j = 1; j <= 6; ++j) {
		vt[j - 1] = 0.;
		st[j - 1] = 0.;
		ft[j - 1] = 0;
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___18);
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
	    for (j = 1; j <= 6; ++j) {
		i__2 = do_lio(&c__5, &c__1, (char *)&vt[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&st[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = do_lio(&c__3, &c__1, (char *)&ft[j - 1], (ftnlen)
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
	    size = 4;
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    numeral_(&id, pd, &size, (ftnlen)4);
	    if (ib < ic) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (ic < ib) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (ia <= id) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (id < ia) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    }
	    torphase_(ft, vt, st);
	    if (header) {
		header = FALSE_;
		io___25.ciunit = iounit_1.iout;
		s_wsfe(&io___25);
		e_wsfe();
	    }
	    if (iring == 0) {
		io___26.ciunit = iounit_1.iout;
		s_wsfe(&io___26);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		for (j = 1; j <= 6; ++j) {
		    do_fio(&c__1, (char *)&vt[j - 1], (ftnlen)sizeof(
			    doublereal));
		    i__2 = i_dnnt(&st[j - 1]);
		    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
		}
		e_wsfe();
	    } else {
		if (iring == 5) {
		    s_copy(label, "5-Ring ", (ftnlen)7, (ftnlen)7);
		}
		if (iring == 4) {
		    s_copy(label, "4-Ring ", (ftnlen)7, (ftnlen)7);
		}
		io___28.ciunit = iounit_1.iout;
		s_wsfe(&io___28);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		for (j = 1; j <= 6; ++j) {
		    do_fio(&c__1, (char *)&vt[j - 1], (ftnlen)sizeof(
			    doublereal));
		    i__2 = i_dnnt(&st[j - 1]);
		    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
		}
		do_fio(&c__1, label, (ftnlen)6);
		e_wsfe();
	    }
	    if (iring == 0) {
		for (j = 1; j <= 2000; ++j) {
		    if (s_cmp(kt_ref(0, j), blank, (ftnlen)16, (ftnlen)16) == 
			    0 || s_cmp(kt_ref(0, j), pt, (ftnlen)16, (ftnlen)
			    16) == 0) {
			s_copy(kt_ref(0, j), pt, (ftnlen)16, (ftnlen)16);
			t1_ref(1, j) = vt[0];
			t1_ref(2, j) = st[0];
			t2_ref(1, j) = vt[1];
			t2_ref(2, j) = st[1];
			t3_ref(1, j) = vt[2];
			t3_ref(2, j) = st[2];
			t4_ref(1, j) = vt[3];
			t4_ref(2, j) = st[3];
			t5_ref(1, j) = vt[4];
			t5_ref(2, j) = st[4];
			t6_ref(1, j) = vt[5];
			t6_ref(2, j) = st[5];
			goto L60;
		    }
		}
		io___29.ciunit = iounit_1.iout;
		s_wsfe(&io___29);
		e_wsfe();
		inform_1.abort = TRUE_;
L60:
		;
	    } else if (iring == 5) {
		for (j = 1; j <= 500; ++j) {
		    if (s_cmp(kt5_ref(0, j), blank, (ftnlen)16, (ftnlen)16) ==
			     0 || s_cmp(kt5_ref(0, j), pt, (ftnlen)16, (
			    ftnlen)16) == 0) {
			s_copy(kt5_ref(0, j), pt, (ftnlen)16, (ftnlen)16);
			t15_ref(1, j) = vt[0];
			t15_ref(2, j) = st[0];
			t25_ref(1, j) = vt[1];
			t25_ref(2, j) = st[1];
			t35_ref(1, j) = vt[2];
			t35_ref(2, j) = st[2];
			t45_ref(1, j) = vt[3];
			t45_ref(2, j) = st[3];
			t55_ref(1, j) = vt[4];
			t55_ref(2, j) = st[4];
			t65_ref(1, j) = vt[5];
			t65_ref(2, j) = st[5];
			goto L80;
		    }
		}
		io___30.ciunit = iounit_1.iout;
		s_wsfe(&io___30);
		e_wsfe();
		inform_1.abort = TRUE_;
L80:
		;
	    } else if (iring == 4) {
		for (j = 1; j <= 500; ++j) {
		    if (s_cmp(kt4_ref(0, j), blank, (ftnlen)16, (ftnlen)16) ==
			     0 || s_cmp(kt4_ref(0, j), pt, (ftnlen)16, (
			    ftnlen)16) == 0) {
			s_copy(kt4_ref(0, j), pt, (ftnlen)16, (ftnlen)16);
			t14_ref(1, j) = vt[0];
			t14_ref(2, j) = st[0];
			t24_ref(1, j) = vt[1];
			t24_ref(2, j) = st[1];
			t34_ref(1, j) = vt[2];
			t34_ref(2, j) = st[2];
			t44_ref(1, j) = vt[3];
			t44_ref(2, j) = st[3];
			t54_ref(1, j) = vt[4];
			t54_ref(2, j) = st[4];
			t64_ref(1, j) = vt[5];
			t64_ref(2, j) = st[5];
			goto L100;
		    }
		}
		io___31.ciunit = iounit_1.iout;
		s_wsfe(&io___31);
		e_wsfe();
		inform_1.abort = TRUE_;
L100:
		;
	    }
	}
    }

/*     determine the total number of forcefield parameters */

    nt = 2000;
    nt5 = 500;
    nt4 = 500;
    for (i__ = 2000; i__ >= 1; --i__) {
	if (s_cmp(kt_ref(0, i__), blank, (ftnlen)16, (ftnlen)16) == 0) {
	    nt = i__ - 1;
	}
    }
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(kt5_ref(0, i__), blank, (ftnlen)16, (ftnlen)16) == 0) {
	    nt5 = i__ - 1;
	}
    }
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(kt4_ref(0, i__), blank, (ftnlen)16, (ftnlen)16) == 0) {
	    nt4 = i__ - 1;
	}
    }
    use_ring__ = FALSE_;
    if (min(nt5,nt4) != 0) {
	use_ring__ = TRUE_;
    }

/*     assign torsional parameters for each torsional angle */
/*     by putting the parameter values into the "tors" arrays */

    header = TRUE_;
    nlist = 0;
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
	i__4[1] = 12, a__2[1] = pt + 4;
	s_cat(pt2, a__2, i__4, &c__2, (ftnlen)16);
/* Writing concatenation */
	i__4[0] = 12, a__2[0] = pt;
	i__4[1] = 4, a__2[1] = zeros;
	s_cat(pt1, a__2, i__4, &c__2, (ftnlen)16);
/* Writing concatenation */
	i__5[0] = 4, a__3[0] = zeros;
	i__5[1] = 8, a__3[1] = pt + 4;
	i__5[2] = 4, a__3[2] = zeros;
	s_cat(pt0, a__3, i__5, &c__3, (ftnlen)16);
	tors1_ref(1, i__) = 0.;
	tors1_ref(2, i__) = 0.;
	tors2_ref(1, i__) = 0.;
	tors2_ref(2, i__) = 0.;
	tors3_ref(1, i__) = 0.;
	tors3_ref(2, i__) = 0.;
	tors4_ref(1, i__) = 0.;
	tors4_ref(2, i__) = 0.;
	tors5_ref(1, i__) = 0.;
	tors5_ref(2, i__) = 0.;
	tors6_ref(1, i__) = 0.;
	tors6_ref(2, i__) = 0.;
	done = FALSE_;

/*     make a check for torsions inside small rings */

	iring = 0;
	if (use_ring__) {
	    chkring_(&iring, &ia, &ib, &ic, &id);
	    if (iring == 6) {
		iring = 0;
	    }
	    if (iring == 5 && nt5 == 0) {
		iring = 0;
	    }
	    if (iring == 4 && nt4 == 0) {
		iring = 0;
	    }
	}

/*     find parameters for this torsion; first check "klist" */
/*     to save time for angle types already located */

	if (iring == 0) {
	    i__2 = nlist;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(klist_ref(0, j), pt, (ftnlen)16, (ftnlen)16) == 0) {
		    ilist = kindex[j - 1];
		    tors1_ref(1, i__) = tors1_ref(1, ilist);
		    tors1_ref(2, i__) = tors1_ref(2, ilist);
		    tors2_ref(1, i__) = tors2_ref(1, ilist);
		    tors2_ref(2, i__) = tors2_ref(2, ilist);
		    tors3_ref(1, i__) = tors3_ref(1, ilist);
		    tors3_ref(2, i__) = tors3_ref(2, ilist);
		    tors4_ref(1, i__) = tors4_ref(1, ilist);
		    tors4_ref(2, i__) = tors4_ref(2, ilist);
		    tors5_ref(1, i__) = tors5_ref(1, ilist);
		    tors5_ref(2, i__) = tors5_ref(2, ilist);
		    tors6_ref(1, i__) = tors6_ref(1, ilist);
		    tors6_ref(2, i__) = tors6_ref(2, ilist);
		    done = TRUE_;
		    goto L110;
		}
	    }
	    i__2 = nt;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kt_ref(0, j), pt, (ftnlen)16, (ftnlen)16) == 0) {
		    ++nlist;
		    s_copy(klist_ref(0, nlist), pt, (ftnlen)16, (ftnlen)16);
		    kindex[nlist - 1] = i__;
		    tors1_ref(1, i__) = t1_ref(1, j);
		    tors1_ref(2, i__) = t1_ref(2, j);
		    tors2_ref(1, i__) = t2_ref(1, j);
		    tors2_ref(2, i__) = t2_ref(2, j);
		    tors3_ref(1, i__) = t3_ref(1, j);
		    tors3_ref(2, i__) = t3_ref(2, j);
		    tors4_ref(1, i__) = t4_ref(1, j);
		    tors4_ref(2, i__) = t4_ref(2, j);
		    tors5_ref(1, i__) = t5_ref(1, j);
		    tors5_ref(2, i__) = t5_ref(2, j);
		    tors6_ref(1, i__) = t6_ref(1, j);
		    tors6_ref(2, i__) = t6_ref(2, j);
		    done = TRUE_;
		    goto L110;
		}
	    }
	    i__2 = nt;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kt_ref(0, j), pt1, (ftnlen)16, (ftnlen)16) == 0 || 
			s_cmp(kt_ref(0, j), pt2, (ftnlen)16, (ftnlen)16) == 0)
			 {
		    tors1_ref(1, i__) = t1_ref(1, j);
		    tors1_ref(2, i__) = t1_ref(2, j);
		    tors2_ref(1, i__) = t2_ref(1, j);
		    tors2_ref(2, i__) = t2_ref(2, j);
		    tors3_ref(1, i__) = t3_ref(1, j);
		    tors3_ref(2, i__) = t3_ref(2, j);
		    tors4_ref(1, i__) = t4_ref(1, j);
		    tors4_ref(2, i__) = t4_ref(2, j);
		    tors5_ref(1, i__) = t5_ref(1, j);
		    tors5_ref(2, i__) = t5_ref(2, j);
		    tors6_ref(1, i__) = t6_ref(1, j);
		    tors6_ref(2, i__) = t6_ref(2, j);
		    done = TRUE_;
		    goto L110;
		}
	    }
	    i__2 = nt;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kt_ref(0, j), pt0, (ftnlen)16, (ftnlen)16) == 0) {
		    tors1_ref(1, i__) = t1_ref(1, j);
		    tors1_ref(2, i__) = t1_ref(2, j);
		    tors2_ref(1, i__) = t2_ref(1, j);
		    tors2_ref(2, i__) = t2_ref(2, j);
		    tors3_ref(1, i__) = t3_ref(1, j);
		    tors3_ref(2, i__) = t3_ref(2, j);
		    tors4_ref(1, i__) = t4_ref(1, j);
		    tors4_ref(2, i__) = t4_ref(2, j);
		    tors5_ref(1, i__) = t5_ref(1, j);
		    tors5_ref(2, i__) = t5_ref(2, j);
		    tors6_ref(1, i__) = t6_ref(1, j);
		    tors6_ref(2, i__) = t6_ref(2, j);
		    done = TRUE_;
		    goto L110;
		}
	    }

/*     find the parameters for a 5-ring torsion */

	} else if (iring == 5) {
	    i__2 = nt5;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kt5_ref(0, j), pt, (ftnlen)16, (ftnlen)16) == 0) {
		    tors1_ref(1, i__) = t15_ref(1, j);
		    tors1_ref(2, i__) = t15_ref(2, j);
		    tors2_ref(1, i__) = t25_ref(1, j);
		    tors2_ref(2, i__) = t25_ref(2, j);
		    tors3_ref(1, i__) = t35_ref(1, j);
		    tors3_ref(2, i__) = t35_ref(2, j);
		    tors4_ref(1, i__) = t45_ref(1, j);
		    tors4_ref(2, i__) = t45_ref(2, j);
		    tors5_ref(1, i__) = t55_ref(1, j);
		    tors5_ref(2, i__) = t55_ref(2, j);
		    tors6_ref(1, i__) = t65_ref(1, j);
		    tors6_ref(2, i__) = t65_ref(2, j);
		    done = TRUE_;
		    goto L110;
		}
	    }
	    i__2 = nt5;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kt5_ref(0, j), pt1, (ftnlen)16, (ftnlen)16) == 0 || 
			s_cmp(kt5_ref(0, j), pt2, (ftnlen)16, (ftnlen)16) == 
			0) {
		    tors1_ref(1, i__) = t15_ref(1, j);
		    tors1_ref(2, i__) = t15_ref(2, j);
		    tors2_ref(1, i__) = t25_ref(1, j);
		    tors2_ref(2, i__) = t25_ref(2, j);
		    tors3_ref(1, i__) = t35_ref(1, j);
		    tors3_ref(2, i__) = t35_ref(2, j);
		    tors4_ref(1, i__) = t45_ref(1, j);
		    tors4_ref(2, i__) = t45_ref(2, j);
		    tors5_ref(1, i__) = t55_ref(1, j);
		    tors5_ref(2, i__) = t55_ref(2, j);
		    tors6_ref(1, i__) = t65_ref(1, j);
		    tors6_ref(2, i__) = t65_ref(2, j);
		    done = TRUE_;
		    goto L110;
		}
	    }
	    i__2 = nt5;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kt5_ref(0, j), pt0, (ftnlen)16, (ftnlen)16) == 0) {
		    tors1_ref(1, i__) = t15_ref(1, j);
		    tors1_ref(2, i__) = t15_ref(2, j);
		    tors2_ref(1, i__) = t25_ref(1, j);
		    tors2_ref(2, i__) = t25_ref(2, j);
		    tors3_ref(1, i__) = t35_ref(1, j);
		    tors3_ref(2, i__) = t35_ref(2, j);
		    tors4_ref(1, i__) = t45_ref(1, j);
		    tors4_ref(2, i__) = t45_ref(2, j);
		    tors5_ref(1, i__) = t55_ref(1, j);
		    tors5_ref(2, i__) = t55_ref(2, j);
		    tors6_ref(1, i__) = t65_ref(1, j);
		    tors6_ref(2, i__) = t65_ref(2, j);
		    done = TRUE_;
		    goto L110;
		}
	    }

/*     find the parameters for a 4-ring torsion */

	} else if (iring == 4) {
	    i__2 = nt4;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kt4_ref(0, j), pt, (ftnlen)16, (ftnlen)16) == 0) {
		    tors1_ref(1, i__) = t14_ref(1, j);
		    tors1_ref(2, i__) = t14_ref(2, j);
		    tors2_ref(1, i__) = t24_ref(1, j);
		    tors2_ref(2, i__) = t24_ref(2, j);
		    tors3_ref(1, i__) = t34_ref(1, j);
		    tors3_ref(2, i__) = t34_ref(2, j);
		    tors4_ref(1, i__) = t44_ref(1, j);
		    tors4_ref(2, i__) = t44_ref(2, j);
		    tors5_ref(1, i__) = t54_ref(1, j);
		    tors5_ref(2, i__) = t54_ref(2, j);
		    tors6_ref(1, i__) = t64_ref(1, j);
		    tors6_ref(2, i__) = t64_ref(2, j);
		    done = TRUE_;
		    goto L110;
		}
	    }
	    i__2 = nt4;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kt4_ref(0, j), pt1, (ftnlen)16, (ftnlen)16) == 0 || 
			s_cmp(kt4_ref(0, j), pt2, (ftnlen)16, (ftnlen)16) == 
			0) {
		    tors1_ref(1, i__) = t14_ref(1, j);
		    tors1_ref(2, i__) = t14_ref(2, j);
		    tors2_ref(1, i__) = t24_ref(1, j);
		    tors2_ref(2, i__) = t24_ref(2, j);
		    tors3_ref(1, i__) = t34_ref(1, j);
		    tors3_ref(2, i__) = t34_ref(2, j);
		    tors4_ref(1, i__) = t44_ref(1, j);
		    tors4_ref(2, i__) = t44_ref(2, j);
		    tors5_ref(1, i__) = t54_ref(1, j);
		    tors5_ref(2, i__) = t54_ref(2, j);
		    tors6_ref(1, i__) = t64_ref(1, j);
		    tors6_ref(2, i__) = t64_ref(2, j);
		    done = TRUE_;
		    goto L110;
		}
	    }
	    i__2 = nt4;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kt4_ref(0, j), pt0, (ftnlen)16, (ftnlen)16) == 0) {
		    tors1_ref(1, i__) = t14_ref(1, j);
		    tors1_ref(2, i__) = t14_ref(2, j);
		    tors2_ref(1, i__) = t24_ref(1, j);
		    tors2_ref(2, i__) = t24_ref(2, j);
		    tors3_ref(1, i__) = t34_ref(1, j);
		    tors3_ref(2, i__) = t34_ref(2, j);
		    tors4_ref(1, i__) = t44_ref(1, j);
		    tors4_ref(2, i__) = t44_ref(2, j);
		    tors5_ref(1, i__) = t54_ref(1, j);
		    tors5_ref(2, i__) = t54_ref(2, j);
		    tors6_ref(1, i__) = t64_ref(1, j);
		    tors6_ref(2, i__) = t64_ref(2, j);
		    done = TRUE_;
		    goto L110;
		}
	    }
	}

/*     warning if suitable torsional parameter not found */

L110:
/* Computing MIN */
	i__2 = atmtyp_1.atomic[ia - 1], i__6 = atmtyp_1.atomic[ib - 1], i__2 =
		 min(i__2,i__6), i__6 = atmtyp_1.atomic[ic - 1], i__2 = min(
		i__2,i__6), i__6 = atmtyp_1.atomic[id - 1];
	minat = min(i__2,i__6);
	if (minat == 0) {
	    done = TRUE_;
	}
	if (! done) {
	    if (potent_1.use_tors__) {
		if (usage_1.use[ia - 1] || usage_1.use[ib - 1] || usage_1.use[
			ic - 1] || usage_1.use[id - 1]) {
		    inform_1.abort = TRUE_;
		}
	    }
	    if (header) {
		header = FALSE_;
		io___49.ciunit = iounit_1.iout;
		s_wsfe(&io___49);
		e_wsfe();
	    }
	    s_copy(label, "Torsion", (ftnlen)7, (ftnlen)7);
	    if (iring == 5) {
		s_copy(label, "5-Ring ", (ftnlen)7, (ftnlen)7);
	    }
	    if (iring == 4) {
		s_copy(label, "4-Ring ", (ftnlen)7, (ftnlen)7);
	    }
	    io___50.ciunit = iounit_1.iout;
	    s_wsfe(&io___50);
	    do_fio(&c__1, label, (ftnlen)7);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ia), (ftnlen)3);
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ib), (ftnlen)3);
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ic), (ftnlen)3);
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, id), (ftnlen)3);
	    do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }

/*     find the cosine and sine of phase angle for each torsion */

    i__1 = tors_1.ntors;
    for (i__ = 1; i__ <= i__1; ++i__) {
	angle = tors1_ref(2, i__) / 57.29577951308232088;
	tors1_ref(3, i__) = cos(angle);
	tors1_ref(4, i__) = sin(angle);
	angle = tors2_ref(2, i__) / 57.29577951308232088;
	tors2_ref(3, i__) = cos(angle);
	tors2_ref(4, i__) = sin(angle);
	angle = tors3_ref(2, i__) / 57.29577951308232088;
	tors3_ref(3, i__) = cos(angle);
	tors3_ref(4, i__) = sin(angle);
	angle = tors4_ref(2, i__) / 57.29577951308232088;
	tors4_ref(3, i__) = cos(angle);
	tors4_ref(4, i__) = sin(angle);
	angle = tors5_ref(2, i__) / 57.29577951308232088;
	tors5_ref(3, i__) = cos(angle);
	tors5_ref(4, i__) = sin(angle);
	angle = tors6_ref(2, i__) / 57.29577951308232088;
	tors6_ref(3, i__) = cos(angle);
	tors6_ref(4, i__) = sin(angle);
    }

/*     turn off the torsional potential if it is not used */

    if (tors_1.ntors == 0) {
	potent_1.use_tors__ = FALSE_;
    }
    return 0;
} /* ktors_ */

#undef keyline_ref
#undef itors_ref
#undef klist_ref
#undef tors6_ref
#undef tors5_ref
#undef tors4_ref
#undef tors3_ref
#undef tors2_ref
#undef tors1_ref
#undef name___ref
#undef kt5_ref
#undef kt4_ref
#undef kt_ref
#undef t64_ref
#undef t54_ref
#undef t44_ref
#undef t34_ref
#undef t24_ref
#undef t65_ref
#undef t55_ref
#undef t45_ref
#undef t35_ref
#undef t25_ref
#undef t15_ref
#undef t14_ref
#undef t6_ref
#undef t5_ref
#undef t4_ref
#undef t3_ref
#undef t2_ref
#undef t1_ref


