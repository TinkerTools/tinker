/* kgeom.f -- translated by f2c (version 20050501).
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
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

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
    doublereal xpfix[25000], ypfix[25000], zpfix[25000], pfix[50000]	/* 
	    was [2][25000] */, dfix[75000]	/* was [3][25000] */, afix[
	    75000]	/* was [3][25000] */, tfix[75000]	/* was [3][
	    25000] */, gfix[75000]	/* was [3][25000] */, chir[75000]	
	    /* was [3][25000] */, depth, width, rwall;
    integer npfix, ipfix[25000], kpfix[75000]	/* was [3][25000] */, ndfix, 
	    idfix[50000]	/* was [2][25000] */, nafix, iafix[75000]	
	    /* was [3][25000] */, ntfix, itfix[100000]	/* was [4][25000] */, 
	    ngfix, igfix[50000]	/* was [2][25000] */, nchir, ichir[100000]	
	    /* was [4][25000] */;
    logical use_basin__, use_wall__;
} kgeoms_;

#define kgeoms_1 kgeoms_

struct {
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

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
static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine kgeom  --  restraint term parameter assignment  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "kgeom" asisgns parameters for geometric restraint terms */
/*     to be included in the potential energy calculation */


/* Subroutine */ int kgeom_(void)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 KGEOM  --  Too many Position Restraints"
	    ";\002,\002 Increase MAXFIX\002)";
    static char fmt_60[] = "(/,\002 KGEOM  --  Too many Distance Restraints"
	    ";\002,\002 Increase MAXFIX\002)";
    static char fmt_90[] = "(/,\002 KGEOM  --  Too many Angle Restraints;"
	    "\002,\002 Increase MAXFIX\002)";
    static char fmt_120[] = "(/,\002 KGEOM  --  Too many Torsion Restraints"
	    ";\002,\002 Increase MAXFIX\002)";
    static char fmt_150[] = "(/,\002 KGEOM  --  Too many Group Restraints"
	    ";\002,\002 Increase MAXFIX\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    static logical intermol;
    extern doublereal geometry_(integer *, integer *, integer *, integer *);
    static integer i__, j, k;
    static doublereal a1, a2, a3, d1, d2, d3, g1, g2, g3, c1, c2, c3, p1, p2, 
	    p3, p4, p5, t1, t2, t3;
    static integer ia, ib, ic, id, ip;
    static doublereal xr, yr, zr, xad, yad, zad, xbd, ybd, zbd, xcd, ycd, zcd,
	     xcm, ycm, zcm, vol;
    static logical keep;
    static integer next;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *), fatal_(void);
    static doublereal weigh, ratio;
    static logical exist;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char letter[1], string[120];
    extern /* Subroutine */ int getword_(char *, char *, integer *, ftnlen, 
	    ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___12 = { 1, string, 1, 0, 120, 1 };
    static icilist io___14 = { 1, string, 1, 0, 120, 1 };
    static icilist io___15 = { 1, string, 1, 0, 120, 1 };
    static cilist io___16 = { 0, 0, 0, fmt_30, 0 };
    static icilist io___21 = { 1, string, 1, 0, 120, 1 };
    static icilist io___24 = { 1, string, 1, 0, 120, 1 };
    static cilist io___29 = { 0, 0, 0, fmt_60, 0 };
    static icilist io___33 = { 1, string, 1, 0, 120, 1 };
    static icilist io___35 = { 1, string, 1, 0, 120, 1 };
    static cilist io___36 = { 0, 0, 0, fmt_90, 0 };
    static icilist io___40 = { 1, string, 1, 0, 120, 1 };
    static icilist io___42 = { 1, string, 1, 0, 120, 1 };
    static cilist io___43 = { 0, 0, 0, fmt_120, 0 };
    static icilist io___47 = { 1, string, 1, 0, 120, 1 };
    static icilist io___48 = { 1, string, 1, 0, 120, 1 };
    static cilist io___55 = { 0, 0, 0, fmt_150, 0 };
    static icilist io___71 = { 1, string, 1, 0, 120, 1 };
    static icilist io___72 = { 1, string, 1, 0, 120, 1 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define chir_ref(a_1,a_2) kgeoms_1.chir[(a_2)*3 + a_1 - 4]
#define afix_ref(a_1,a_2) kgeoms_1.afix[(a_2)*3 + a_1 - 4]
#define dfix_ref(a_1,a_2) kgeoms_1.dfix[(a_2)*3 + a_1 - 4]
#define gfix_ref(a_1,a_2) kgeoms_1.gfix[(a_2)*3 + a_1 - 4]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]
#define pfix_ref(a_1,a_2) kgeoms_1.pfix[(a_2)*2 + a_1 - 3]
#define tfix_ref(a_1,a_2) kgeoms_1.tfix[(a_2)*3 + a_1 - 4]
#define ichir_ref(a_1,a_2) kgeoms_1.ichir[(a_2)*4 + a_1 - 5]
#define iafix_ref(a_1,a_2) kgeoms_1.iafix[(a_2)*3 + a_1 - 4]
#define idfix_ref(a_1,a_2) kgeoms_1.idfix[(a_2)*2 + a_1 - 3]
#define igfix_ref(a_1,a_2) kgeoms_1.igfix[(a_2)*2 + a_1 - 3]
#define kpfix_ref(a_1,a_2) kgeoms_1.kpfix[(a_2)*3 + a_1 - 4]
#define itfix_ref(a_1,a_2) kgeoms_1.itfix[(a_2)*4 + a_1 - 5]
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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  group.i  --  partitioning of system into atom groups  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     grpmass     total mass of all the atoms in each group */
/*     wgrp        weight for each set of group-group interactions */
/*     ngrp        total number of atom groups in the system */
/*     kgrp        contiguous list of the atoms in each group */
/*     igrp        first and last atom of each group in the list */
/*     grplist     number of the group to which each atom belongs */
/*     use_group   flag to use partitioning of system into groups */
/*     use_intra   flag to include only intragroup interactions */
/*     use_inter   flag to include only intergroup interactions */




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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  kgeoms.i  --  parameters for the geometrical restraints  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xpfix      x-coordinate target for each restrained position */
/*     ypfix      y-coordinate target for each restrained position */
/*     zpfix      z-coordinate target for each restrained position */
/*     pfix       force constant and flat-well range for each position */
/*     dfix       force constant and target range for each distance */
/*     afix       force constant and target range for each angle */
/*     tfix       force constant and target range for each torsion */
/*     gfix       force constant and target range for each group distance */
/*     chir       force constant and target range for chiral centers */
/*     depth      depth of shallow Gaussian basin restraint */
/*     width      exponential width coefficient of Gaussian basin */
/*     rwall      radius of spherical droplet boundary restraint */
/*     npfix      number of position restraints to be applied */
/*     ipfix      atom number involved in each position restraint */
/*     kpfix      flags to use x-, y-, z-coordinate position restraints */
/*     ndfix      number of distance restraints to be applied */
/*     idfix      atom numbers defining each distance restraint */
/*     nafix      number of angle restraints to be applied */
/*     iafix      atom numbers defining each angle restraint */
/*     ntfix      number of torsional restraints to be applied */
/*     itfix      atom numbers defining each torsional restraint */
/*     ngfix      number of group distance restraints to be applied */
/*     igfix      group numbers defining each group distance restraint */
/*     nchir      number of chirality restraints to be applied */
/*     ichir      atom numbers defining each chirality restraint */
/*     use_basin  logical flag governing use of Gaussian basin */
/*     use_wall   logical flag governing use of droplet boundary */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  molcul.i  --  individual molecules within current system  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     molmass   molecular weight for each molecule in the system */
/*     totmass   total weight of all the molecules in the system */
/*     nmol      total number of separate molecules in the system */
/*     kmol      contiguous list of the atoms in each molecule */
/*     imol      first and last atom of each molecule in the list */
/*     molcule   number of the molecule to which each atom belongs */




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




/*     set the default values for the restraint variables */

    kgeoms_1.npfix = 0;
    kgeoms_1.ndfix = 0;
    kgeoms_1.nafix = 0;
    kgeoms_1.ntfix = 0;
    kgeoms_1.ngfix = 0;
    kgeoms_1.nchir = 0;
    kgeoms_1.depth = 0.;
    kgeoms_1.width = 0.;
    kgeoms_1.use_basin__ = FALSE_;
    kgeoms_1.use_wall__ = FALSE_;

/*     search the keywords for restraint parameters */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));

/*     get atom restrained to a specified position range */

	if (s_cmp(keyword, "RESTRAIN-POSITION ", (ftnlen)18, (ftnlen)18) == 0)
		 {
	    p1 = 0.;
	    p2 = 0.;
	    p3 = 0.;
	    p4 = 0.;
	    p5 = 0.;
	    next = 1;
	    getword_(string, letter, &next, (ftnlen)120, (ftnlen)1);
	    if (*(unsigned char *)letter == ' ') {
		i__2 = s_rsli(&io___12);
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = do_lio(&c__3, &c__1, (char *)&ip, (ftnlen)sizeof(
			integer));
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&p1, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&p2, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&p3, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&p4, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&p5, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L10;
		}
L10:
		if (p4 == 0.) {
		    p4 = 100.;
		}
		++kgeoms_1.npfix;
		kgeoms_1.ipfix[kgeoms_1.npfix - 1] = ip;
		kpfix_ref(1, kgeoms_1.npfix) = 1;
		kpfix_ref(2, kgeoms_1.npfix) = 1;
		kpfix_ref(3, kgeoms_1.npfix) = 1;
		kgeoms_1.xpfix[kgeoms_1.npfix - 1] = p1;
		kgeoms_1.ypfix[kgeoms_1.npfix - 1] = p2;
		kgeoms_1.zpfix[kgeoms_1.npfix - 1] = p3;
		pfix_ref(1, kgeoms_1.npfix) = p4;
		pfix_ref(2, kgeoms_1.npfix) = p5;
	    } else {
		upcase_(letter, (ftnlen)1);
		i__2 = s_rsli(&io___14);
		if (i__2 != 0) {
		    goto L20;
		}
		i__2 = do_lio(&c__3, &c__1, (char *)&ip, (ftnlen)sizeof(
			integer));
		if (i__2 != 0) {
		    goto L20;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L20;
		}
		s_copy(string, string + (next - 1), (ftnlen)120, 120 - (next 
			- 1));
		i__2 = s_rsli(&io___15);
		if (i__2 != 0) {
		    goto L20;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&p1, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L20;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&p2, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L20;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&p3, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L20;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L20;
		}
L20:
		if (p2 == 0.) {
		    p2 = 100.;
		}
		++kgeoms_1.npfix;
		kgeoms_1.ipfix[kgeoms_1.npfix - 1] = ip;
		kpfix_ref(1, kgeoms_1.npfix) = 0;
		kpfix_ref(2, kgeoms_1.npfix) = 0;
		kpfix_ref(3, kgeoms_1.npfix) = 0;
		if (*(unsigned char *)letter == 'X') {
		    kpfix_ref(1, kgeoms_1.npfix) = 1;
		    kgeoms_1.xpfix[kgeoms_1.npfix - 1] = p1;
		} else if (*(unsigned char *)letter == 'Y') {
		    kpfix_ref(2, kgeoms_1.npfix) = 1;
		    kgeoms_1.ypfix[kgeoms_1.npfix - 1] = p1;
		} else if (*(unsigned char *)letter == 'Z') {
		    kpfix_ref(3, kgeoms_1.npfix) = 1;
		    kgeoms_1.zpfix[kgeoms_1.npfix - 1] = p1;
		}
		pfix_ref(1, kgeoms_1.npfix) = p2;
		pfix_ref(2, kgeoms_1.npfix) = p3;
	    }
	    if (kgeoms_1.npfix > 25000) {
		io___16.ciunit = iounit_1.iout;
		s_wsfe(&io___16);
		e_wsfe();
		fatal_();
	    }

/*     get atoms restrained to a specified distance range */

	} else if (s_cmp(keyword, "RESTRAIN-DISTANCE ", (ftnlen)18, (ftnlen)
		18) == 0) {
	    d1 = 0.;
	    d2 = 0.;
	    d3 = 0.;
	    exist = FALSE_;
	    i__2 = s_rsli(&io___21);
	    if (i__2 != 0) {
		goto L40;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L40;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L40;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&d1, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L40;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&d2, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L40;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L40;
	    }
	    exist = TRUE_;
L40:
	    i__2 = s_rsli(&io___24);
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
	    i__2 = do_lio(&c__5, &c__1, (char *)&d1, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&d2, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&d3, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L50;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L50;
	    }
L50:
	    if (d1 == 0.) {
		d1 = 100.;
	    }
	    if (! exist) {
		xr = atoms_1.x[ia - 1] - atoms_1.x[ib - 1];
		yr = atoms_1.y[ia - 1] - atoms_1.y[ib - 1];
		zr = atoms_1.z__[ia - 1] - atoms_1.z__[ib - 1];
		intermol = molcul_1.molcule[ia - 1] != molcul_1.molcule[ib - 
			1];
		if (bound_1.use_bounds__ && intermol) {
		    image_(&xr, &yr, &zr);
		}
		d2 = sqrt(xr * xr + yr * yr + zr * zr);
	    }
	    if (d3 == 0.) {
		d3 = d2;
	    }
	    ++kgeoms_1.ndfix;
	    idfix_ref(1, kgeoms_1.ndfix) = ia;
	    idfix_ref(2, kgeoms_1.ndfix) = ib;
	    dfix_ref(1, kgeoms_1.ndfix) = d1;
	    dfix_ref(2, kgeoms_1.ndfix) = d2;
	    dfix_ref(3, kgeoms_1.ndfix) = d3;
	    if (kgeoms_1.ndfix > 25000) {
		io___29.ciunit = iounit_1.iout;
		s_wsfe(&io___29);
		e_wsfe();
		fatal_();
	    }

/*     get atoms restrained to a specified angle range */

	} else if (s_cmp(keyword, "RESTRAIN-ANGLE ", (ftnlen)15, (ftnlen)15) 
		== 0) {
	    a1 = 0.;
	    a2 = 0.;
	    a3 = 0.;
	    exist = FALSE_;
	    i__2 = s_rsli(&io___33);
	    if (i__2 != 0) {
		goto L70;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L70;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L70;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L70;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&a1, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L70;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&a2, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L70;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L70;
	    }
	    exist = TRUE_;
L70:
	    i__2 = s_rsli(&io___35);
	    if (i__2 != 0) {
		goto L80;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L80;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L80;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L80;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&a1, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L80;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&a2, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L80;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&a3, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L80;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L80;
	    }
L80:
	    if (a1 == 0.) {
		a1 = 10.;
	    }
	    if (! exist) {
		a2 = geometry_(&ia, &ib, &ic, &c__0);
	    }
	    if (a3 == 0.) {
		a3 = a2;
	    }
	    ++kgeoms_1.nafix;
	    iafix_ref(1, kgeoms_1.nafix) = ia;
	    iafix_ref(2, kgeoms_1.nafix) = ib;
	    iafix_ref(3, kgeoms_1.nafix) = ic;
	    afix_ref(1, kgeoms_1.nafix) = a1;
	    afix_ref(2, kgeoms_1.nafix) = a2;
	    afix_ref(3, kgeoms_1.nafix) = a3;
	    if (kgeoms_1.nafix > 25000) {
		io___36.ciunit = iounit_1.iout;
		s_wsfe(&io___36);
		e_wsfe();
		fatal_();
	    }

/*     get atoms restrained to a specified torsion range */

	} else if (s_cmp(keyword, "RESTRAIN-TORSION ", (ftnlen)17, (ftnlen)17)
		 == 0) {
	    t1 = 0.;
	    t2 = 0.;
	    t3 = 0.;
	    exist = FALSE_;
	    i__2 = s_rsli(&io___40);
	    if (i__2 != 0) {
		goto L100;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L100;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L100;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L100;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L100;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&t1, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L100;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&t2, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L100;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L100;
	    }
	    exist = TRUE_;
L100:
	    i__2 = s_rsli(&io___42);
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
	    i__2 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&t1, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&t2, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&t3, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L110;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L110;
	    }
	    exist = TRUE_;
L110:
	    if (t1 == 0.) {
		t1 = 1.;
	    }
	    if (! exist) {
		t2 = geometry_(&ia, &ib, &ic, &id);
	    }
	    if (t3 == 0.) {
		t3 = t2;
	    }
	    while(t2 > 180.) {
		t2 += -360.;
	    }
	    while(t2 < -180.) {
		t2 += 360.;
	    }
	    while(t3 > 180.) {
		t3 += -360.;
	    }
	    while(t3 < -180.) {
		t3 += 360.;
	    }
	    ++kgeoms_1.ntfix;
	    itfix_ref(1, kgeoms_1.ntfix) = ia;
	    itfix_ref(2, kgeoms_1.ntfix) = ib;
	    itfix_ref(3, kgeoms_1.ntfix) = ic;
	    itfix_ref(4, kgeoms_1.ntfix) = id;
	    tfix_ref(1, kgeoms_1.ntfix) = t1;
	    tfix_ref(2, kgeoms_1.ntfix) = t2;
	    tfix_ref(3, kgeoms_1.ntfix) = t3;
	    if (kgeoms_1.ntfix > 25000) {
		io___43.ciunit = iounit_1.iout;
		s_wsfe(&io___43);
		e_wsfe();
		fatal_();
	    }

/*     get groups restrained to a specified distance range */

	} else if (s_cmp(keyword, "RESTRAIN-GROUPS ", (ftnlen)16, (ftnlen)16) 
		== 0) {
	    g1 = 0.;
	    g2 = 0.;
	    g3 = 0.;
	    exist = FALSE_;
	    i__2 = s_rsli(&io___47);
	    if (i__2 != 0) {
		goto L130;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L130;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L130;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&g1, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L130;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&g2, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L130;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L130;
	    }
	    exist = TRUE_;
L130:
	    i__2 = s_rsli(&io___48);
	    if (i__2 != 0) {
		goto L140;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L140;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L140;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&g1, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L140;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&g2, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L140;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&g3, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L140;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L140;
	    }
L140:
	    if (g1 == 0.) {
		g1 = 100.;
	    }
	    if (! exist) {
		xcm = 0.;
		ycm = 0.;
		zcm = 0.;
		i__2 = igrp_ref(2, ia);
		for (j = igrp_ref(1, ia); j <= i__2; ++j) {
		    k = group_1.kgrp[j - 1];
		    weigh = atmtyp_1.mass[k - 1];
		    xcm += atoms_1.x[k - 1] * weigh;
		    ycm += atoms_1.y[k - 1] * weigh;
		    zcm += atoms_1.z__[k - 1] * weigh;
		}
/* Computing MAX */
		d__1 = 1., d__2 = group_1.grpmass[ia - 1];
		weigh = max(d__1,d__2);
		xr = xcm / weigh;
		yr = ycm / weigh;
		zr = zcm / weigh;
		xcm = 0.;
		ycm = 0.;
		zcm = 0.;
		i__2 = igrp_ref(2, ib);
		for (j = igrp_ref(1, ib); j <= i__2; ++j) {
		    k = group_1.kgrp[j - 1];
		    weigh = atmtyp_1.mass[k - 1];
		    xcm += atoms_1.x[k - 1] * weigh;
		    ycm += atoms_1.y[k - 1] * weigh;
		    zcm += atoms_1.z__[k - 1] * weigh;
		}
/* Computing MAX */
		d__1 = 1., d__2 = group_1.grpmass[ib - 1];
		weigh = max(d__1,d__2);
		xr -= xcm / weigh;
		yr -= ycm / weigh;
		zr -= zcm / weigh;
		intermol = molcul_1.molcule[group_1.kgrp[igrp_ref(1, ia) - 1] 
			- 1] != molcul_1.molcule[group_1.kgrp[igrp_ref(1, ib) 
			- 1] - 1];
		if (bound_1.use_bounds__ && intermol) {
		    image_(&xr, &yr, &zr);
		}
		g2 = sqrt(xr * xr + yr * yr + zr * zr);
	    }
	    if (g3 == 0.) {
		g3 = g2;
	    }
	    ++kgeoms_1.ngfix;
	    igfix_ref(1, kgeoms_1.ngfix) = ia;
	    igfix_ref(2, kgeoms_1.ngfix) = ib;
	    gfix_ref(1, kgeoms_1.ngfix) = g1;
	    gfix_ref(2, kgeoms_1.ngfix) = g2;
	    gfix_ref(3, kgeoms_1.ngfix) = g3;
	    if (kgeoms_1.ngfix > 25000) {
		io___55.ciunit = iounit_1.iout;
		s_wsfe(&io___55);
		e_wsfe();
		fatal_();
	    }

/*     maintain chirality as found in the original input structure */

	} else if (s_cmp(keyword, "ENFORCE-CHIRALITY ", (ftnlen)18, (ftnlen)
		18) == 0) {
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		if (couple_1.n12[j - 1] == 4) {
		    ia = i12_ref(1, j);
		    ib = i12_ref(2, j);
		    ic = i12_ref(3, j);
		    id = i12_ref(4, j);
		    keep = TRUE_;
		    if (couple_1.n12[ia - 1] == 1) {
			if (atoms_1.type__[ia - 1] == atoms_1.type__[ib - 1]) 
				{
			    keep = FALSE_;
			}
			if (atoms_1.type__[ia - 1] == atoms_1.type__[ic - 1]) 
				{
			    keep = FALSE_;
			}
			if (atoms_1.type__[ia - 1] == atoms_1.type__[id - 1]) 
				{
			    keep = FALSE_;
			}
		    } else if (couple_1.n12[ib - 1] == 1) {
			if (atoms_1.type__[ib - 1] == atoms_1.type__[ic - 1]) 
				{
			    keep = FALSE_;
			}
			if (atoms_1.type__[ib - 1] == atoms_1.type__[id - 1]) 
				{
			    keep = FALSE_;
			}
		    } else if (couple_1.n12[ic - 1] == 1) {
			if (atoms_1.type__[ic - 1] == atoms_1.type__[id - 1]) 
				{
			    keep = FALSE_;
			}
		    }
		    if (keep) {
			++kgeoms_1.nchir;
			ichir_ref(1, kgeoms_1.nchir) = ia;
			ichir_ref(2, kgeoms_1.nchir) = ib;
			ichir_ref(3, kgeoms_1.nchir) = ic;
			ichir_ref(4, kgeoms_1.nchir) = id;
			xad = atoms_1.x[ia - 1] - atoms_1.x[id - 1];
			yad = atoms_1.y[ia - 1] - atoms_1.y[id - 1];
			zad = atoms_1.z__[ia - 1] - atoms_1.z__[id - 1];
			xbd = atoms_1.x[ib - 1] - atoms_1.x[id - 1];
			ybd = atoms_1.y[ib - 1] - atoms_1.y[id - 1];
			zbd = atoms_1.z__[ib - 1] - atoms_1.z__[id - 1];
			xcd = atoms_1.x[ic - 1] - atoms_1.x[id - 1];
			ycd = atoms_1.y[ic - 1] - atoms_1.y[id - 1];
			zcd = atoms_1.z__[ic - 1] - atoms_1.z__[id - 1];
			c1 = ybd * zcd - zbd * ycd;
			c2 = ycd * zad - zcd * yad;
			c3 = yad * zbd - zad * ybd;
			vol = xad * c1 + xbd * c2 + xcd * c3;
			ratio = (d__1 = vol / (xad * xbd * xcd), abs(d__1));
			chir_ref(1, kgeoms_1.nchir) = 10.;
			if (ratio > .1) {
			    chir_ref(2, kgeoms_1.nchir) = vol * .5;
			    chir_ref(3, kgeoms_1.nchir) = vol * 2.;
			} else {
			    chir_ref(2, kgeoms_1.nchir) = abs(vol) * -2.;
			    chir_ref(3, kgeoms_1.nchir) = abs(vol) * 2.;
			}
		    }
		}
	    }

/*     setup any shallow Gaussian basin restraint between atoms */

	} else if (s_cmp(keyword, "BASIN ", (ftnlen)6, (ftnlen)6) == 0) {
	    i__2 = s_rsli(&io___71);
	    if (i__2 != 0) {
		goto L160;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&kgeoms_1.depth, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L160;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&kgeoms_1.width, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L160;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L160;
	    }
L160:
	    kgeoms_1.use_basin__ = TRUE_;
	    if (kgeoms_1.depth == 0.) {
		kgeoms_1.use_basin__ = FALSE_;
	    }
	    if (kgeoms_1.width == 0.) {
		kgeoms_1.use_basin__ = FALSE_;
	    }
	    if (kgeoms_1.depth > 0.) {
		kgeoms_1.depth = -kgeoms_1.depth;
	    }

/*     setup any spherical droplet restraint between atoms */

	} else if (s_cmp(keyword, "WALL ", (ftnlen)5, (ftnlen)5) == 0) {
	    i__2 = s_rsli(&io___72);
	    if (i__2 != 0) {
		goto L170;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&kgeoms_1.rwall, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L170;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L170;
	    }
L170:
	    if (kgeoms_1.rwall > 0.) {
		kgeoms_1.use_wall__ = TRUE_;
	    }
	}
    }

/*     turn on the geometric restraint potential if it is used */

    potent_1.use_geom__ = FALSE_;
    if (kgeoms_1.npfix != 0) {
	potent_1.use_geom__ = TRUE_;
    }
    if (kgeoms_1.ndfix != 0) {
	potent_1.use_geom__ = TRUE_;
    }
    if (kgeoms_1.nafix != 0) {
	potent_1.use_geom__ = TRUE_;
    }
    if (kgeoms_1.ntfix != 0) {
	potent_1.use_geom__ = TRUE_;
    }
    if (kgeoms_1.ngfix != 0) {
	potent_1.use_geom__ = TRUE_;
    }
    if (kgeoms_1.nchir != 0) {
	potent_1.use_geom__ = TRUE_;
    }
    if (kgeoms_1.use_basin__) {
	potent_1.use_geom__ = TRUE_;
    }
    if (kgeoms_1.use_wall__) {
	potent_1.use_geom__ = TRUE_;
    }
    return 0;
} /* kgeom_ */

#undef keyline_ref
#undef itfix_ref
#undef kpfix_ref
#undef igfix_ref
#undef idfix_ref
#undef iafix_ref
#undef ichir_ref
#undef tfix_ref
#undef pfix_ref
#undef igrp_ref
#undef gfix_ref
#undef dfix_ref
#undef afix_ref
#undef chir_ref
#undef i12_ref


