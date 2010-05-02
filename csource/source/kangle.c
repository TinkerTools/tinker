/* kangle.f -- translated by f2c (version 20050501).
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
    doublereal ak[75000], anat[75000], afld[75000];
    integer nangle, iang[300000]	/* was [4][75000] */;
} angle_;

#define angle_1 angle_

struct {
    doublereal angunit, stbnunit, aaunit, opbunit, opdunit, cang, qang, pang, 
	    sang, copb, qopb, popb, sopb, copd, qopd, popd, sopd;
    char angtyp[600000], opbtyp[8];
} angpot_;

#define angpot_1 angpot_

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
    doublereal acon[2000], acon5[500], acon4[500], acon3[500], aconf[500], 
	    ang[6000]	/* was [3][2000] */, ang5[1500]	/* was [3][500] */, 
	    ang4[1500]	/* was [3][500] */, ang3[1500]	/* was [3][500] */, 
	    angf[1000]	/* was [2][500] */;
    char ka[24000], ka5[6000], ka4[6000], ka3[6000], kaf[6000];
} kangs_;

#define kangs_1 kangs_

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

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine kangle  --  angle bend parameter assignment  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "kangle" assigns the force constants and ideal angles for */
/*     the bond angles; also processes new or changed parameters */


/* Subroutine */ int kangle_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Additional Angle Bending Parameters :"
	    "\002,//,5x,\002Atom Classes\002,9x,\002K(B)\002,7x,\002Angle\002"
	    ",/)";
    static char fmt_30[] = "(4x,3i4,2x,2f12.3)";
    static char fmt_40[] = "(4x,3i4,2x,2f12.3,3x,\0020-H's\002)";
    static char fmt_50[] = "(4x,3i4,2x,2f12.3,3x,\0021-H's\002)";
    static char fmt_60[] = "(4x,3i4,2x,2f12.3,3x,\0022-H's\002)";
    static char fmt_70[] = "(4x,3i4,2x,2f12.3,3x,a6)";
    static char fmt_80[] = "(4x,3i4,2x,2f12.3,3x,a6,3x,\0020-H's\002)";
    static char fmt_90[] = "(4x,3i4,2x,2f12.3,3x,a6,3x,\0021-H's\002)";
    static char fmt_100[] = "(4x,3i4,2x,2f12.3,3x,a6,3x,\0022-H's\002)";
    static char fmt_110[] = "(/,\002 KANGLE  --  Too many Bond Angle\002,"
	    "\002 Bending Parameters\002)";
    static char fmt_120[] = "(/,\002 KANGLE  --  Too many 5-Ring Angle\002"
	    ",\002 Bending Parameters\002)";
    static char fmt_130[] = "(/,\002 KANGLE  --  Too many 4-Ring Angle\002"
	    ",\002 Bending Parameters\002)";
    static char fmt_140[] = "(/,\002 KANGLE  --  Too many 3-Ring Angle\002"
	    ",\002 Bending Parameters\002)";
    static char fmt_170[] = "(/,\002 Additional Fourier Angle Bending\002"
	    ",\002 Parameters :\002,//,5x,\002Atom Classes\002,9x,\002K(B)"
	    "\002,7x,\002Shift\002,6x,\002Period\002,/)";
    static char fmt_180[] = "(4x,3i4,2x,3f12.3)";
    static char fmt_190[] = "(/,\002 KANGLE  --  Too many Fourier Angle\002"
	    ",\002 Bending Parameters\002)";
    static char fmt_220[] = "(/,\002 Undefined Angle Bending Parameters :"
	    "\002,//,\002 Type\002,18x,\002Atom Names\002,19x,\002Atom Classes"
	    "\002,/)";
    static char fmt_230[] = "(1x,a6,5x,3(i6,\002-\002,a3),7x,3i5)";

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
    static logical use_ring__;
    static integer i__, j;
    static doublereal fc;
    static integer ia, ib, ic, na, ih, nh;
    static doublereal an;
    static char pa[4], pb[4], pc[4];
    static doublereal pr;
    static char pt[12];
    static integer na4, na5, na3;
    static doublereal an1, an2, an3;
    static integer naf, jen, ita, itb, itc;
    static logical done;
    static integer size, next;
    static char label[6], blank[12];
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
    static icilist io___17 = { 1, string, 1, 0, 120, 1 };
    static cilist io___18 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_140, 0 };
    static icilist io___40 = { 1, string, 1, 0, 120, 1 };
    static cilist io___41 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___57 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___58 = { 0, 0, 0, fmt_230, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define ka_ref(a_0,a_1) &kangs_1.ka[(a_1)*12 + a_0 - 12]
#define ka3_ref(a_0,a_1) &kangs_1.ka3[(a_1)*12 + a_0 - 12]
#define ka4_ref(a_0,a_1) &kangs_1.ka4[(a_1)*12 + a_0 - 12]
#define ka5_ref(a_0,a_1) &kangs_1.ka5[(a_1)*12 + a_0 - 12]
#define kaf_ref(a_0,a_1) &kangs_1.kaf[(a_1)*12 + a_0 - 12]
#define ang_ref(a_1,a_2) kangs_1.ang[(a_2)*3 + a_1 - 4]
#define ang3_ref(a_1,a_2) kangs_1.ang3[(a_2)*3 + a_1 - 4]
#define ang4_ref(a_1,a_2) kangs_1.ang4[(a_2)*3 + a_1 - 4]
#define ang5_ref(a_1,a_2) kangs_1.ang5[(a_2)*3 + a_1 - 4]
#define angf_ref(a_1,a_2) kangs_1.angf[(a_2)*2 + a_1 - 3]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define angtyp_ref(a_0,a_1) &angpot_1.angtyp[(a_1)*8 + a_0 - 8]
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
/*     ##  angpot.i  --  specifics of angle bend functional forms  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     angunit    convert angle bending energy to kcal/mole */
/*     stbnunit   convert stretch-bend energy to kcal/mole */
/*     aaunit     convert angle-angle energy to kcal/mole */
/*     opbunit    convert out-of-plane bend energy to kcal/mole */
/*     opdunit    convert out-of-plane distance energy to kcal/mole */
/*     cang       cubic coefficient in angle bending potential */
/*     qang       quartic coefficient in angle bending potential */
/*     pang       quintic coefficient in angle bending potential */
/*     sang       sextic coefficient in angle bending potential */
/*     copb       cubic coefficient in out-of-plane bend potential */
/*     qopb       quartic coefficient in out-of-plane bend potential */
/*     popb       quintic coefficient in out-of-plane bend potential */
/*     sopb       sextic coefficient in out-of-plane bend potential */
/*     copd       cubic coefficient in out-of-plane distance potential */
/*     qopd       quartic coefficient in out-of-plane distance potential */
/*     popd       quintic coefficient in out-of-plane distance potential */
/*     sopd       sextic coefficient in out-of-plane distance potential */
/*     angtyp     type of angle bending function for each bond angle */
/*     opbtyp     type of out-of-plane bend potential energy function */




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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kangs.i  --  forcefield parameters for bond angle bending  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxna    maximum number of harmonic angle bend parameter entries */
/*     maxna5   maximum number of 5-membered ring angle bend entries */
/*     maxna4   maximum number of 4-membered ring angle bend entries */
/*     maxna3   maximum number of 3-membered ring angle bend entries */
/*     maxnaf   maximum number of Fourier angle bend parameter entries */

/*     acon     force constant parameters for harmonic angle bends */
/*     acon5    force constant parameters for 5-ring angle bends */
/*     acon4    force constant parameters for 4-ring angle bends */
/*     acon3    force constant parameters for 3-ring angle bends */
/*     aconf    force constant parameters for Fourier angle bends */
/*     ang      bond angle parameters for harmonic angle bends */
/*     ang5     bond angle parameters for 5-ring angle bends */
/*     ang4     bond angle parameters for 4-ring angle bends */
/*     ang3     bond angle parameters for 3-ring angle bends */
/*     angf     phase shift angle and periodicity for Fourier bends */
/*     ka       string of atom classes for harmonic angle bends */
/*     ka5      string of atom classes for 5-ring angle bends */
/*     ka4      string of atom classes for 4-ring angle bends */
/*     ka3      string of atom classes for 3-ring angle bends */
/*     kaf      string of atom classes for Fourier angle bends */




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




/*     process keywords containing angle bending parameters */

    s_copy(blank, "         ", (ftnlen)12, (ftnlen)9);
    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	iring = -1;
	if (s_cmp(keyword, "ANGLE ", (ftnlen)6, (ftnlen)6) == 0) {
	    iring = 0;
	}
	if (s_cmp(keyword, "ANGLE5 ", (ftnlen)7, (ftnlen)7) == 0) {
	    iring = 5;
	}
	if (s_cmp(keyword, "ANGLE4 ", (ftnlen)7, (ftnlen)7) == 0) {
	    iring = 4;
	}
	if (s_cmp(keyword, "ANGLE3 ", (ftnlen)7, (ftnlen)7) == 0) {
	    iring = 3;
	}
	if (iring >= 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    an1 = 0.;
	    an2 = 0.;
	    an3 = 0.;
	    jen = 0;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___17);
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
	    i__2 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&an1, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&an2, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&an3, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
L10:
	    if (an2 != 0. || an3 != 0.) {
		jen = 1;
	    }
	    if (header) {
		header = FALSE_;
		io___18.ciunit = iounit_1.iout;
		s_wsfe(&io___18);
		e_wsfe();
	    }
	    if (iring == 0) {
		if (jen == 0) {
		    io___19.ciunit = iounit_1.iout;
		    s_wsfe(&io___19);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else if (an1 != 0.) {
		    io___20.ciunit = iounit_1.iout;
		    s_wsfe(&io___20);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
		if (an2 != 0.) {
		    io___21.ciunit = iounit_1.iout;
		    s_wsfe(&io___21);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an2, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
		if (an3 != 0.) {
		    io___22.ciunit = iounit_1.iout;
		    s_wsfe(&io___22);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an3, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
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
		if (jen == 0) {
		    io___24.ciunit = iounit_1.iout;
		    s_wsfe(&io___24);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, label, (ftnlen)6);
		    e_wsfe();
		} else if (an1 != 0.) {
		    io___25.ciunit = iounit_1.iout;
		    s_wsfe(&io___25);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, label, (ftnlen)6);
		    e_wsfe();
		}
		if (an2 != 0.) {
		    io___26.ciunit = iounit_1.iout;
		    s_wsfe(&io___26);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, label, (ftnlen)6);
		    e_wsfe();
		}
		if (an3 != 0.) {
		    io___27.ciunit = iounit_1.iout;
		    s_wsfe(&io___27);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an3, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, label, (ftnlen)6);
		    e_wsfe();
		}
	    }
	    size = 4;
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    if (ia <= ic) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pc;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pa;
		s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	    }
	    if (iring == 0) {
		for (j = 1; j <= 2000; ++j) {
		    if (s_cmp(ka_ref(0, j), blank, (ftnlen)12, (ftnlen)12) == 
			    0 || s_cmp(ka_ref(0, j), pt, (ftnlen)12, (ftnlen)
			    12) == 0) {
			s_copy(ka_ref(0, j), pt, (ftnlen)12, (ftnlen)12);
			kangs_1.acon[j - 1] = fc;
			ang_ref(1, j) = an1;
			ang_ref(2, j) = an2;
			ang_ref(3, j) = an3;
			goto L150;
		    }
		}
		io___34.ciunit = iounit_1.iout;
		s_wsfe(&io___34);
		e_wsfe();
		inform_1.abort = TRUE_;
	    } else if (iring == 5) {
		for (j = 1; j <= 500; ++j) {
		    if (s_cmp(ka5_ref(0, j), blank, (ftnlen)12, (ftnlen)12) ==
			     0 || s_cmp(ka5_ref(0, j), pt, (ftnlen)12, (
			    ftnlen)12) == 0) {
			s_copy(ka5_ref(0, j), pt, (ftnlen)12, (ftnlen)12);
			kangs_1.acon5[j - 1] = fc;
			ang5_ref(1, j) = an1;
			ang5_ref(2, j) = an2;
			ang5_ref(3, j) = an3;
			goto L150;
		    }
		}
		io___35.ciunit = iounit_1.iout;
		s_wsfe(&io___35);
		e_wsfe();
		inform_1.abort = TRUE_;
	    } else if (iring == 4) {
		for (j = 1; j <= 500; ++j) {
		    if (s_cmp(ka4_ref(0, j), blank, (ftnlen)12, (ftnlen)12) ==
			     0 || s_cmp(ka4_ref(0, j), pt, (ftnlen)12, (
			    ftnlen)12) == 0) {
			s_copy(ka4_ref(0, j), pt, (ftnlen)12, (ftnlen)12);
			kangs_1.acon4[j - 1] = fc;
			ang4_ref(1, j) = an1;
			ang4_ref(2, j) = an2;
			ang4_ref(3, j) = an3;
			goto L150;
		    }
		}
		io___36.ciunit = iounit_1.iout;
		s_wsfe(&io___36);
		e_wsfe();
		inform_1.abort = TRUE_;
	    } else if (iring == 3) {
		for (j = 1; j <= 500; ++j) {
		    if (s_cmp(ka3_ref(0, j), blank, (ftnlen)12, (ftnlen)12) ==
			     0 || s_cmp(ka3_ref(0, j), pt, (ftnlen)12, (
			    ftnlen)12) == 0) {
			s_copy(ka3_ref(0, j), pt, (ftnlen)12, (ftnlen)12);
			kangs_1.acon3[j - 1] = fc;
			ang3_ref(1, j) = an1;
			ang3_ref(2, j) = an2;
			ang3_ref(3, j) = an3;
			goto L150;
		    }
		}
		io___37.ciunit = iounit_1.iout;
		s_wsfe(&io___37);
		e_wsfe();
		inform_1.abort = TRUE_;
	    }
L150:
	    ;
	}
    }

/*     process keywords containing Fourier angle bending parameters */

    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	iring = -1;
	if (s_cmp(keyword, "ANGLEF ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    an = 0.;
	    pr = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___40);
	    if (i__2 != 0) {
		goto L160;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L160;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L160;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L160;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L160;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&an, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L160;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&pr, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L160;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L160;
	    }
L160:
	    if (header) {
		header = FALSE_;
		io___41.ciunit = iounit_1.iout;
		s_wsfe(&io___41);
		e_wsfe();
	    }
	    io___42.ciunit = iounit_1.iout;
	    s_wsfe(&io___42);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&an, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&pr, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    size = 4;
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    if (ia <= ic) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pc;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pa;
		s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	    }
	    for (j = 1; j <= 500; ++j) {
		if (s_cmp(kaf_ref(0, j), blank, (ftnlen)12, (ftnlen)12) == 0 
			|| s_cmp(kaf_ref(0, j), pt, (ftnlen)12, (ftnlen)12) ==
			 0) {
		    s_copy(kaf_ref(0, j), pt, (ftnlen)12, (ftnlen)12);
		    kangs_1.aconf[j - 1] = fc;
		    angf_ref(1, j) = an;
		    angf_ref(2, j) = pr;
		    goto L200;
		}
	    }
	    io___43.ciunit = iounit_1.iout;
	    s_wsfe(&io___43);
	    e_wsfe();
	    inform_1.abort = TRUE_;
L200:
	    ;
	}
    }

/*     determine the total number of forcefield parameters */

    na = 2000;
    na5 = 500;
    na4 = 500;
    na3 = 500;
    naf = 500;
    for (i__ = 2000; i__ >= 1; --i__) {
	if (s_cmp(ka_ref(0, i__), blank, (ftnlen)12, (ftnlen)12) == 0) {
	    na = i__ - 1;
	}
    }
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(ka5_ref(0, i__), blank, (ftnlen)12, (ftnlen)12) == 0) {
	    na5 = i__ - 1;
	}
    }
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(ka4_ref(0, i__), blank, (ftnlen)12, (ftnlen)12) == 0) {
	    na4 = i__ - 1;
	}
    }
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(ka3_ref(0, i__), blank, (ftnlen)12, (ftnlen)12) == 0) {
	    na3 = i__ - 1;
	}
    }
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(kaf_ref(0, i__), blank, (ftnlen)12, (ftnlen)12) == 0) {
	    naf = i__ - 1;
	}
    }
    use_ring__ = FALSE_;
/* Computing MIN */
    i__1 = min(na5,na4);
    if (min(i__1,na3) != 0) {
	use_ring__ = TRUE_;
    }

/*     set generic parameters for use with any number of hydrogens */

    i__1 = na;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ang_ref(2, i__) == 0. && ang_ref(3, i__) == 0.) {
	    ang_ref(2, i__) = ang_ref(1, i__);
	    ang_ref(3, i__) = ang_ref(1, i__);
	}
    }
    i__1 = na5;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ang5_ref(2, i__) == 0. && ang5_ref(3, i__) == 0.) {
	    ang5_ref(2, i__) = ang5_ref(1, i__);
	    ang5_ref(3, i__) = ang5_ref(1, i__);
	}
    }
    i__1 = na4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ang4_ref(2, i__) == 0. && ang4_ref(3, i__) == 0.) {
	    ang4_ref(2, i__) = ang4_ref(1, i__);
	    ang4_ref(3, i__) = ang4_ref(1, i__);
	}
    }
    i__1 = na3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ang3_ref(2, i__) == 0. && ang3_ref(3, i__) == 0.) {
	    ang3_ref(2, i__) = ang3_ref(1, i__);
	    ang3_ref(3, i__) = ang3_ref(1, i__);
	}
    }

/*     assign ideal bond angle and force constant for each angle */

    header = TRUE_;
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
	if (ita <= itc) {
/* Writing concatenation */
	    i__3[0] = 4, a__1[0] = pa;
	    i__3[1] = 4, a__1[1] = pb;
	    i__3[2] = 4, a__1[2] = pc;
	    s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	} else {
/* Writing concatenation */
	    i__3[0] = 4, a__1[0] = pc;
	    i__3[1] = 4, a__1[1] = pb;
	    i__3[2] = 4, a__1[2] = pa;
	    s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	}
	angle_1.ak[i__ - 1] = 0.;
	angle_1.anat[i__ - 1] = 0.;
	angle_1.afld[i__ - 1] = 0.;
	s_copy(angtyp_ref(0, i__), "HARMONIC", (ftnlen)8, (ftnlen)8);
	done = FALSE_;

/*     count number of non-angle hydrogens on the central atom */

	nh = 1;
	i__2 = couple_1.n12[ib - 1];
	for (j = 1; j <= i__2; ++j) {
	    ih = i12_ref(j, ib);
	    if (ih != ia && ih != ic && atmtyp_1.atomic[ih - 1] == 1) {
		++nh;
	    }
	}

/*     make a check for bond angles contained inside small rings */

	iring = 0;
	if (use_ring__) {
	    chkring_(&iring, &ia, &ib, &ic, &c__0);
	    if (iring == 6) {
		iring = 0;
	    }
	    if (iring == 5 && na5 == 0) {
		iring = 0;
	    }
	    if (iring == 4 && na4 == 0) {
		iring = 0;
	    }
	    if (iring == 3 && na3 == 0) {
		iring = 0;
	    }
	}

/*     assign angle bending parameters for bond angles */

	if (iring == 0) {
	    i__2 = na;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(ka_ref(0, j), pt, (ftnlen)12, (ftnlen)12) == 0 && 
			ang_ref(nh, j) != 0.) {
		    angle_1.ak[i__ - 1] = kangs_1.acon[j - 1];
		    angle_1.anat[i__ - 1] = ang_ref(nh, j);
		    done = TRUE_;
		    goto L210;
		}
	    }

/*     assign bending parameters for 5-membered ring angles */

	} else if (iring == 5) {
	    i__2 = na5;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(ka5_ref(0, j), pt, (ftnlen)12, (ftnlen)12) == 0 && 
			ang5_ref(nh, j) != 0.) {
		    angle_1.ak[i__ - 1] = kangs_1.acon5[j - 1];
		    angle_1.anat[i__ - 1] = ang5_ref(nh, j);
		    done = TRUE_;
		    goto L210;
		}
	    }

/*     assign bending parameters for 4-membered ring angles */

	} else if (iring == 4) {
	    i__2 = na4;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(ka4_ref(0, j), pt, (ftnlen)12, (ftnlen)12) == 0 && 
			ang4_ref(nh, j) != 0.) {
		    angle_1.ak[i__ - 1] = kangs_1.acon4[j - 1];
		    angle_1.anat[i__ - 1] = ang4_ref(nh, j);
		    done = TRUE_;
		    goto L210;
		}
	    }

/*     assign bending parameters for 3-membered ring angles */

	} else if (iring == 3) {
	    i__2 = na3;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(ka3_ref(0, j), pt, (ftnlen)12, (ftnlen)12) == 0 && 
			ang3_ref(nh, j) != 0.) {
		    angle_1.ak[i__ - 1] = kangs_1.acon3[j - 1];
		    angle_1.anat[i__ - 1] = ang3_ref(nh, j);
		    done = TRUE_;
		    goto L210;
		}
	    }
	}

/*     assign Fourier angle bending parameters for bond angles */

	if (! done) {
	    i__2 = naf;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(kaf_ref(0, j), pt, (ftnlen)12, (ftnlen)12) == 0) {
		    angle_1.ak[i__ - 1] = kangs_1.aconf[j - 1];
		    angle_1.anat[i__ - 1] = angf_ref(1, j);
		    angle_1.afld[i__ - 1] = angf_ref(2, j);
		    s_copy(angtyp_ref(0, i__), "FOURIER", (ftnlen)8, (ftnlen)
			    7);
		    done = TRUE_;
		    goto L210;
		}
	    }
	}

/*     warning if suitable angle bending parameter not found */

L210:
/* Computing MIN */
	i__2 = atmtyp_1.atomic[ia - 1], i__4 = atmtyp_1.atomic[ib - 1], i__2 =
		 min(i__2,i__4), i__4 = atmtyp_1.atomic[ic - 1];
	minat = min(i__2,i__4);
	if (minat == 0) {
	    done = TRUE_;
	}
	if (potent_1.use_angle__ && ! done) {
	    if (usage_1.use[ia - 1] || usage_1.use[ib - 1] || usage_1.use[ic 
		    - 1]) {
		inform_1.abort = TRUE_;
	    }
	    if (header) {
		header = FALSE_;
		io___57.ciunit = iounit_1.iout;
		s_wsfe(&io___57);
		e_wsfe();
	    }
	    s_copy(label, "Angle ", (ftnlen)6, (ftnlen)6);
	    if (iring == 5) {
		s_copy(label, "5-Ring", (ftnlen)6, (ftnlen)6);
	    }
	    if (iring == 4) {
		s_copy(label, "4-Ring", (ftnlen)6, (ftnlen)6);
	    }
	    if (iring == 3) {
		s_copy(label, "3-Ring", (ftnlen)6, (ftnlen)6);
	    }
	    io___58.ciunit = iounit_1.iout;
	    s_wsfe(&io___58);
	    do_fio(&c__1, label, (ftnlen)6);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ia), (ftnlen)3);
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ib), (ftnlen)3);
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ic), (ftnlen)3);
	    do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }

/*     turn off the angle bending potential if it is not used */

    if (angle_1.nangle == 0) {
	potent_1.use_angle__ = FALSE_;
    }
    return 0;
} /* kangle_ */

#undef keyline_ref
#undef angtyp_ref
#undef name___ref
#undef iang_ref
#undef angf_ref
#undef ang5_ref
#undef ang4_ref
#undef ang3_ref
#undef ang_ref
#undef kaf_ref
#undef ka5_ref
#undef ka4_ref
#undef ka3_ref
#undef ka_ref
#undef i12_ref


