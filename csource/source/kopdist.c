/* kopdist.f -- translated by f2c (version 20050501).
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
    doublereal opds[500];
    char kopd[8000];
} kopdst_;

#define kopdst_1 kopdst_

struct {
    doublereal opdk[25000];
    integer nopdist, iopd[100000]	/* was [4][25000] */;
} opdist_;

#define opdist_1 opdist_

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
static integer c__4 = 4;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine kopdist  --  out-of-plane distance parameters  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "kopdist" assigns the force constants for out-of-plane */
/*     distance at trigonal centers via the central atom height; */
/*     also processes any new or changed parameter values */


/* Subroutine */ int kopdist_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Additional Out-of-Plane Distance\002,"
	    "\002 Parameters :\002,//,5x,\002Atom Classes\002,19x,\002K(OPD"
	    ")\002,/)";
    static char fmt_30[] = "(4x,4i4,10x,2f12.3)";
    static char fmt_40[] = "(/,\002 KOPDIST  --  Too many Out-of-Plane Dista"
	    "nce\002,\002 Parameters\002)";

    /* System generated locals */
    address a__1[4], a__2[2];
    integer i__1, i__2, i__3[4], i__4[2];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, ia, ib, ic, id;
    static char pa[4], pb[4], pc[4], pd[4], pt[16], pt0[16];
    static integer ita, itb, itc, itd;
    static doublereal fopd;
    static integer nopd, size, next;
    static char blank[16];
    static integer itmin;
    static char zeros[12];
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
    static icilist io___14 = { 1, string, 1, 0, 120, 1 };
    static cilist io___25 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_40, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define iopd_ref(a_1,a_2) opdist_1.iopd[(a_2)*4 + a_1 - 5]
#define kopd_ref(a_0,a_1) &kopdst_1.kopd[(a_1)*16 + a_0 - 16]
#define angtyp_ref(a_0,a_1) &angpot_1.angtyp[(a_1)*8 + a_0 - 8]
#define keyline_ref(a_0,a_1) &keys_1.keyline[(a_1)*120 + a_0 - 120]
#define anglist_ref(a_1,a_2) atmlst_1.anglist[(a_2)*28 + a_1 - 29]



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
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kopdst.i  --  forcefield parameters for out-plane distance  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     maxnopd   maximum number of out-of-plane distance entries */

/*     opds      force constant parameters for out-of-plane distance */
/*     kopd      string of atom classes for out-of-plane distance */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  opdist.i  --  out-of-plane distances in current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     opdk      force constant values for out-of-plane distance */
/*     nopdist   total number of out-of-plane distances in the system */
/*     iopb      numbers of the atoms in each out-of-plane distance */




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




/*     process keywords containing out-of-plane distance parameters */

    s_copy(blank, "                ", (ftnlen)16, (ftnlen)16);
    s_copy(zeros, "000000000000", (ftnlen)12, (ftnlen)12);
    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "OPDIST ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    fopd = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___14);
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
	    i__2 = do_lio(&c__5, &c__1, (char *)&fopd, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
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
/* Computing MIN */
	    i__2 = min(itb,itc);
	    itmin = min(i__2,itd);
	    if (itb == itmin) {
		if (itc <= itd) {
/* Writing concatenation */
		    i__3[0] = 4, a__1[0] = pa;
		    i__3[1] = 4, a__1[1] = pb;
		    i__3[2] = 4, a__1[2] = pc;
		    i__3[3] = 4, a__1[3] = pd;
		    s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
		} else {
/* Writing concatenation */
		    i__3[0] = 4, a__1[0] = pa;
		    i__3[1] = 4, a__1[1] = pb;
		    i__3[2] = 4, a__1[2] = pd;
		    i__3[3] = 4, a__1[3] = pc;
		    s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
		}
	    } else if (itc == itmin) {
		if (itb <= itd) {
/* Writing concatenation */
		    i__3[0] = 4, a__1[0] = pa;
		    i__3[1] = 4, a__1[1] = pc;
		    i__3[2] = 4, a__1[2] = pb;
		    i__3[3] = 4, a__1[3] = pd;
		    s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
		} else {
/* Writing concatenation */
		    i__3[0] = 4, a__1[0] = pa;
		    i__3[1] = 4, a__1[1] = pc;
		    i__3[2] = 4, a__1[2] = pd;
		    i__3[3] = 4, a__1[3] = pb;
		    s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
		}
	    } else if (itd == itmin) {
		if (itb <= itc) {
/* Writing concatenation */
		    i__3[0] = 4, a__1[0] = pa;
		    i__3[1] = 4, a__1[1] = pd;
		    i__3[2] = 4, a__1[2] = pb;
		    i__3[3] = 4, a__1[3] = pc;
		    s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
		} else {
/* Writing concatenation */
		    i__3[0] = 4, a__1[0] = pa;
		    i__3[1] = 4, a__1[1] = pd;
		    i__3[2] = 4, a__1[2] = pc;
		    i__3[3] = 4, a__1[3] = pb;
		    s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
		}
	    }
	    if (header) {
		header = FALSE_;
		io___25.ciunit = iounit_1.iout;
		s_wsfe(&io___25);
		e_wsfe();
	    }
	    io___26.ciunit = iounit_1.iout;
	    s_wsfe(&io___26);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&fopd, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    for (j = 1; j <= 500; ++j) {
		if (s_cmp(kopd_ref(0, j), blank, (ftnlen)16, (ftnlen)16) == 0 
			|| s_cmp(kopd_ref(0, j), pt, (ftnlen)16, (ftnlen)16) 
			== 0) {
		    s_copy(kopd_ref(0, j), pt, (ftnlen)16, (ftnlen)16);
		    kopdst_1.opds[j - 1] = fopd;
		    goto L50;
		}
	    }
	    io___28.ciunit = iounit_1.iout;
	    s_wsfe(&io___28);
	    e_wsfe();
	    inform_1.abort = TRUE_;
L50:
	    ;
	}
    }

/*     determine the total number of forcefield parameters */

    nopd = 500;
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(kopd_ref(0, i__), blank, (ftnlen)16, (ftnlen)16) == 0) {
	    nopd = i__ - 1;
	}
    }

/*     assign out-of-plane distance parameters for trigonal sites */

    opdist_1.nopdist = 0;
    if (nopd != 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (couple_1.n12[i__ - 1] == 3) {
		ia = i__;
		ib = i12_ref(1, i__);
		ic = i12_ref(2, i__);
		id = i12_ref(3, i__);
		ita = atmtyp_1.class__[ia - 1];
		itb = atmtyp_1.class__[ib - 1];
		itc = atmtyp_1.class__[ic - 1];
		itd = atmtyp_1.class__[id - 1];
		size = 4;
		numeral_(&ita, pa, &size, (ftnlen)4);
		numeral_(&itb, pb, &size, (ftnlen)4);
		numeral_(&itc, pc, &size, (ftnlen)4);
		numeral_(&itd, pd, &size, (ftnlen)4);
/* Computing MIN */
		i__2 = min(itb,itc);
		itmin = min(i__2,itd);
		if (itb == itmin) {
		    if (itc <= itd) {
/* Writing concatenation */
			i__3[0] = 4, a__1[0] = pa;
			i__3[1] = 4, a__1[1] = pb;
			i__3[2] = 4, a__1[2] = pc;
			i__3[3] = 4, a__1[3] = pd;
			s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
		    } else {
/* Writing concatenation */
			i__3[0] = 4, a__1[0] = pa;
			i__3[1] = 4, a__1[1] = pb;
			i__3[2] = 4, a__1[2] = pd;
			i__3[3] = 4, a__1[3] = pc;
			s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
		    }
		} else if (itc == itmin) {
		    if (itb <= itd) {
/* Writing concatenation */
			i__3[0] = 4, a__1[0] = pa;
			i__3[1] = 4, a__1[1] = pc;
			i__3[2] = 4, a__1[2] = pb;
			i__3[3] = 4, a__1[3] = pd;
			s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
		    } else {
/* Writing concatenation */
			i__3[0] = 4, a__1[0] = pa;
			i__3[1] = 4, a__1[1] = pc;
			i__3[2] = 4, a__1[2] = pd;
			i__3[3] = 4, a__1[3] = pb;
			s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
		    }
		} else if (itd == itmin) {
		    if (itb <= itc) {
/* Writing concatenation */
			i__3[0] = 4, a__1[0] = pa;
			i__3[1] = 4, a__1[1] = pd;
			i__3[2] = 4, a__1[2] = pb;
			i__3[3] = 4, a__1[3] = pc;
			s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
		    } else {
/* Writing concatenation */
			i__3[0] = 4, a__1[0] = pa;
			i__3[1] = 4, a__1[1] = pd;
			i__3[2] = 4, a__1[2] = pc;
			i__3[3] = 4, a__1[3] = pb;
			s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
		    }
		}
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pa;
		i__4[1] = 12, a__2[1] = zeros;
		s_cat(pt0, a__2, i__4, &c__2, (ftnlen)16);
		i__2 = nopd;
		for (j = 1; j <= i__2; ++j) {
		    if (s_cmp(kopd_ref(0, j), pt, (ftnlen)16, (ftnlen)16) == 
			    0) {
			++opdist_1.nopdist;
			iopd_ref(1, opdist_1.nopdist) = ia;
			iopd_ref(2, opdist_1.nopdist) = ib;
			iopd_ref(3, opdist_1.nopdist) = ic;
			iopd_ref(4, opdist_1.nopdist) = id;
			opdist_1.opdk[opdist_1.nopdist - 1] = kopdst_1.opds[j 
				- 1];
			goto L60;
		    }
		}
		i__2 = nopd;
		for (j = 1; j <= i__2; ++j) {
		    if (s_cmp(kopd_ref(0, j), pt0, (ftnlen)16, (ftnlen)16) == 
			    0) {
			++opdist_1.nopdist;
			iopd_ref(1, opdist_1.nopdist) = ia;
			iopd_ref(2, opdist_1.nopdist) = ib;
			iopd_ref(3, opdist_1.nopdist) = ic;
			iopd_ref(4, opdist_1.nopdist) = id;
			opdist_1.opdk[opdist_1.nopdist - 1] = kopdst_1.opds[j 
				- 1];
			goto L60;
		    }
		}
L60:
		;
	    }
	}
    }

/*     mark angles at trigonal sites to use projected in-plane values */

    i__1 = opdist_1.nopdist;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iopd_ref(1, i__);
	for (j = 1; j <= 3; ++j) {
	    k = anglist_ref(j, ia);
	    if (s_cmp(angtyp_ref(0, k), "HARMONIC", (ftnlen)8, (ftnlen)8) == 
		    0) {
		s_copy(angtyp_ref(0, k), "IN-PLANE", (ftnlen)8, (ftnlen)8);
	    }
	}
    }

/*     turn off out-of-plane distance potential if it is not used */

    if (opdist_1.nopdist == 0) {
	potent_1.use_opdist__ = FALSE_;
    }
    return 0;
} /* kopdist_ */

#undef anglist_ref
#undef keyline_ref
#undef angtyp_ref
#undef kopd_ref
#undef iopd_ref
#undef i12_ref


