/* kurey.f -- translated by f2c (version 20050501).
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
    doublereal ucon[2000], dst13[2000];
    char ku[24000];
} kurybr_;

#define kurybr_1 kurybr_

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
    doublereal uk[75000], ul[75000];
    integer nurey, iury[225000]	/* was [3][75000] */;
} urey_;

#define urey_1 urey_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine kurey  --  Urey-Bradley parameter assignment  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "kurey" assigns the force constants and ideal distances */
/*     for the Urey-Bradley 1-3 interactions; also processes any */
/*     new or changed parameter values */


/* Subroutine */ int kurey_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Additional Urey-Bradley Parameters :\002"
	    ",//,5x,\002Atom Classes\002,8x,\002K(UB)\002,5x,\002Distance\002"
	    ",/)";
    static char fmt_30[] = "(4x,3i4,2x,f12.3,f12.4)";
    static char fmt_40[] = "(/,\002 KUREY  --  Too many Urey-Bradley\002,"
	    "\002 Interaction Parameters\002)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2, i__3[3];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i__, j;
    static doublereal bb;
    static integer ia, ib, ic;
    static char pa[4], pb[4], pc[4];
    static integer nu;
    static char pt[12];
    static doublereal tt;
    static integer ita, itb, itc, size, next;
    static char blank[12];
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
    static icilist io___13 = { 1, string, 1, 0, 120, 1 };
    static cilist io___14 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_40, 0 };



#define ku_ref(a_0,a_1) &kurybr_1.ku[(a_1)*12 + a_0 - 12]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define iury_ref(a_1,a_2) urey_1.iury[(a_2)*3 + a_1 - 4]
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
/*     ##  kurybr.i  --  forcefield parameters for Urey-Bradley terms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     maxnu   maximum number of Urey-Bradley parameter entries */

/*     ucon    force constant parameters for Urey-Bradley terms */
/*     dst13   ideal 1-3 distance parameters for Urey-Bradley terms */
/*     ku      string of atom classes for Urey-Bradley terms */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  urey.i  --  Urey-Bradley interactions in the structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     uk      Urey-Bradley force constants (kcal/mole/Ang**2) */
/*     ul      ideal 1-3 distance values in Angstroms */
/*     nurey   total number of Urey-Bradley terms in the system */
/*     iury    numbers of the atoms in each Urey-Bradley interaction */




/*     process keywords containing Urey-Bradley parameters */

    s_copy(blank, "            ", (ftnlen)12, (ftnlen)12);
    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "UREYBRAD ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    bb = 0.;
	    tt = 0.;
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
	    i__2 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&bb, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&tt, (ftnlen)sizeof(
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
	    io___15.ciunit = iounit_1.iout;
	    s_wsfe(&io___15);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&bb, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&tt, (ftnlen)sizeof(doublereal));
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
	    for (j = 1; j <= 2000; ++j) {
		if (s_cmp(ku_ref(0, j), blank, (ftnlen)12, (ftnlen)12) == 0 ||
			 s_cmp(ku_ref(0, j), pt, (ftnlen)12, (ftnlen)12) == 0)
			 {
		    s_copy(ku_ref(0, j), pt, (ftnlen)12, (ftnlen)12);
		    kurybr_1.ucon[j - 1] = bb;
		    kurybr_1.dst13[j - 1] = tt;
		    goto L50;
		}
	    }
	    io___22.ciunit = iounit_1.iout;
	    s_wsfe(&io___22);
	    e_wsfe();
	    inform_1.abort = TRUE_;
L50:
	    ;
	}
    }

/*     determine the total number of forcefield parameters */

    nu = 2000;
    for (i__ = 2000; i__ >= 1; --i__) {
	if (s_cmp(ku_ref(0, i__), blank, (ftnlen)12, (ftnlen)12) == 0) {
	    nu = i__ - 1;
	}
    }

/*     assign the Urey-Bradley parameters for each angle */

    urey_1.nurey = 0;
    if (nu != 0) {
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
	    i__2 = nu;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(ku_ref(0, j), pt, (ftnlen)12, (ftnlen)12) == 0) {
		    ++urey_1.nurey;
		    iury_ref(1, urey_1.nurey) = ia;
		    iury_ref(2, urey_1.nurey) = ib;
		    iury_ref(3, urey_1.nurey) = ic;
		    urey_1.uk[urey_1.nurey - 1] = kurybr_1.ucon[j - 1];
		    urey_1.ul[urey_1.nurey - 1] = kurybr_1.dst13[j - 1];
		    goto L60;
		}
	    }
L60:
	    ;
	}
    }

/*     turn off the Urey-Bradley potential if it is not used */

    if (urey_1.nurey == 0) {
	potent_1.use_urey__ = FALSE_;
    }
    return 0;
} /* kurey_ */

#undef keyline_ref
#undef iury_ref
#undef iang_ref
#undef ku_ref


