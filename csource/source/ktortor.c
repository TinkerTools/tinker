/* ktortor.f -- translated by f2c (version 20050501).
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
    integer nbitor, ibitor[500000]	/* was [5][100000] */;
} bitor_;

#define bitor_1 bitor_

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
    doublereal ttx[3000]	/* was [30][100] */, tty[3000]	/* was [30][
	    100] */, tbf[90000]	/* was [900][100] */, tbx[90000]	/* 
	    was [900][100] */, tby[90000]	/* was [900][100] */, tbxy[
	    90000]	/* was [900][100] */;
    integer tnx[100], tny[100];
    char ktt[2000];
} ktrtor_;

#define ktrtor_1 ktrtor_

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
    integer ntortor, itt[300000]	/* was [3][100000] */;
} tortor_;

#define tortor_1 tortor_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__30 = 30;



/*     ############################################################# */
/*     ##  COPYRIGHT (C) 2003 by Pengyu Ren & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine ktortor  --  tors-tors parameter assignment  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "ktortor" assigns torsion-torsion parameters to adjacent */
/*     torsion pairs and processes any new or changed values */


/* Subroutine */ int ktortor_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Additional Torsion-Torsion Parameters "
	    ":\002,//,5x,\002Atom Classes\002,12x,\002GridSize1\002,5x,\002Gr"
	    "idSize2\002,/)";
    static char fmt_30[] = "(1x,5i4,6x,i8,6x,i8)";
    static char fmt_40[] = "(/,\002 KTORTOR  --  Too many Torsion-Torsion"
	    "\002,\002 Parameters\002)";

    /* System generated locals */
    address a__1[5];
    integer i__1, i__2, i__3, i__4[5];
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i__, j, k, m, ia, ib, ic, id, ie;
    static char pa[4], pb[4], pc[4], pd[4];
    static doublereal bs[31], cs[31], ds[31];
    static char pe[4];
    static doublereal tf[900];
    static char pt[20];
    static integer nx, ny;
    static doublereal tx[900], ty[900];
    static char pt1[20], pt2[20];
    static integer ita, itb, itc, itd, ite;
    static doublereal eps;
    static integer ntt, nxy;
    static doublereal tmp1[31], tmp2[31], tmp3[31], tmp4[31], tmp5[31], tmp6[
	    31], tmp7[31];
    static integer size, next;
    extern /* Subroutine */ int sort9_(integer *, doublereal *);
    static char blank[20];
    static logical header, cyclic;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int cspline_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), numeral_(integer *, char *, integer *, ftnlen), 
	    nspline_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___20 = { 1, string, 1, 0, 120, 1 };
    static icilist io___21 = { 1, record, 1, 0, 120, 1 };
    static cilist io___22 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_40, 0 };



#define tbf_ref(a_1,a_2) ktrtor_1.tbf[(a_2)*900 + a_1 - 901]
#define tbx_ref(a_1,a_2) ktrtor_1.tbx[(a_2)*900 + a_1 - 901]
#define tby_ref(a_1,a_2) ktrtor_1.tby[(a_2)*900 + a_1 - 901]
#define itt_ref(a_1,a_2) tortor_1.itt[(a_2)*3 + a_1 - 4]
#define ktt_ref(a_0,a_1) &ktrtor_1.ktt[(a_1)*20 + a_0 - 20]
#define ttx_ref(a_1,a_2) ktrtor_1.ttx[(a_2)*30 + a_1 - 31]
#define tty_ref(a_1,a_2) ktrtor_1.tty[(a_2)*30 + a_1 - 31]
#define tbxy_ref(a_1,a_2) ktrtor_1.tbxy[(a_2)*900 + a_1 - 901]
#define ibitor_ref(a_1,a_2) bitor_1.ibitor[(a_2)*5 + a_1 - 6]
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
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  bitor.i  --  bitorsions within the current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     nbitor  total number of bitorsions in the system */
/*     ibitor  numbers of the atoms in each bitorsion */




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
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  ktrtor.i  --  forcefield parameters for torsion-torsions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxntt    maximum number of torsion-torsion parameter entries */
/*     maxtgrd   maximum dimension of torsion-torsion spline grid */
/*     maxtgrd2  maximum number of torsion-torsion spline grid points */

/*     ttx       angle values for first torsion of spline grid */
/*     tty       angle values for second torsion of spline grid */
/*     tbf       function values at points on spline grid */
/*     tbx       gradient over first torsion of spline grid */
/*     tby       gradient over second torsion of spline grid */
/*     tbxy      Hessian cross components over spline grid */
/*     tnx       number of columns in torsion-torsion spline grid */
/*     tny       number of rows in torsion-torsion spline grid */
/*     ktt       string of torsion-torsion atom classes */




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
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  tortor.i  --  torsion-torsions in the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     ntortor   total number of torsion-torsion interactions */
/*     itt       atoms and parameter indices for torsion-torsion */




/*     process keywords containing torsion-torsion parameters */

    s_copy(blank, "                    ", (ftnlen)20, (ftnlen)20);
    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "TORTORS ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    ie = 0;
	    nx = 0;
	    ny = 0;
	    nxy = 0;
	    for (j = 1; j <= 900; ++j) {
		tx[j - 1] = 0.;
		ty[j - 1] = 0.;
		tf[j - 1] = 0.;
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___20);
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
	    i__2 = do_lio(&c__3, &c__1, (char *)&ie, (ftnlen)sizeof(integer));
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
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	    nxy = nx * ny;
	    i__2 = nxy;
	    for (j = 1; j <= i__2; ++j) {
		s_copy(record, keyline_ref(0, i__ + j), (ftnlen)120, (ftnlen)
			120);
		i__3 = s_rsli(&io___21);
		if (i__3 != 0) {
		    goto L10;
		}
		i__3 = do_lio(&c__5, &c__1, (char *)&tx[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__3 != 0) {
		    goto L10;
		}
		i__3 = do_lio(&c__5, &c__1, (char *)&ty[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__3 != 0) {
		    goto L10;
		}
		i__3 = do_lio(&c__5, &c__1, (char *)&tf[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__3 != 0) {
		    goto L10;
		}
		i__3 = e_rsli();
		if (i__3 != 0) {
		    goto L10;
		}
	    }
L10:
	    if (header) {
		header = FALSE_;
		io___22.ciunit = iounit_1.iout;
		s_wsfe(&io___22);
		e_wsfe();
	    }
	    io___23.ciunit = iounit_1.iout;
	    s_wsfe(&io___23);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ie, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nx, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ny, (ftnlen)sizeof(integer));
	    e_wsfe();
	    size = 4;
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    numeral_(&id, pd, &size, (ftnlen)4);
	    numeral_(&ie, pe, &size, (ftnlen)4);
/* Writing concatenation */
	    i__4[0] = 4, a__1[0] = pa;
	    i__4[1] = 4, a__1[1] = pb;
	    i__4[2] = 4, a__1[2] = pc;
	    i__4[3] = 4, a__1[3] = pd;
	    i__4[4] = 4, a__1[4] = pe;
	    s_cat(pt, a__1, i__4, &c__5, (ftnlen)20);
	    for (j = 1; j <= 100; ++j) {
		if (s_cmp(ktt_ref(0, j), blank, (ftnlen)20, (ftnlen)20) == 0 
			|| s_cmp(ktt_ref(0, j), pt, (ftnlen)20, (ftnlen)20) ==
			 0) {
		    s_copy(ktt_ref(0, j), pt, (ftnlen)20, (ftnlen)20);
		    nx = nxy;
		    sort9_(&nx, tx);
		    ny = nxy;
		    sort9_(&ny, ty);
		    ktrtor_1.tnx[j - 1] = nx;
		    ktrtor_1.tny[j - 1] = ny;
		    i__2 = nx;
		    for (k = 1; k <= i__2; ++k) {
			ttx_ref(k, j) = tx[k - 1];
		    }
		    i__2 = ny;
		    for (k = 1; k <= i__2; ++k) {
			tty_ref(k, j) = ty[k - 1];
		    }
		    i__2 = nxy;
		    for (k = 1; k <= i__2; ++k) {
			tbf_ref(k, j) = tf[k - 1];
		    }
		    goto L50;
		}
	    }
	    io___32.ciunit = iounit_1.iout;
	    s_wsfe(&io___32);
	    e_wsfe();
	    inform_1.abort = TRUE_;
L50:
	    ;
	}
    }

/*     determine the total number of forcefield parameters */

    ntt = 100;
    for (i__ = 100; i__ >= 1; --i__) {
	if (s_cmp(ktt_ref(0, i__), blank, (ftnlen)20, (ftnlen)20) == 0) {
	    ntt = i__ - 1;
	}
    }

/*     check whether each torsion-torsion parameter is periodic; */
/*     assumes the "tbf" array is sorted with both indices in */
/*     increasing order and the first index changing most rapidly */

    i__1 = ntt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cyclic = TRUE_;
	eps = 1e-4;
	nx = ktrtor_1.tnx[i__ - 1] - 1;
	ny = ktrtor_1.tny[i__ - 1] - 1;
	if ((d__1 = ttx_ref(1, i__) - ttx_ref(ktrtor_1.tnx[i__ - 1], i__), 
		abs(d__1)) - 360. > eps) {
	    cyclic = FALSE_;
	}
	if ((d__1 = tty_ref(1, i__) - tty_ref(ktrtor_1.tny[i__ - 1], i__), 
		abs(d__1)) - 360. > eps) {
	    cyclic = FALSE_;
	}
	if (cyclic) {
	    i__2 = ktrtor_1.tny[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		k = (j - 1) * ktrtor_1.tnx[i__ - 1] + 1;
		if ((d__1 = tbf_ref(k, i__) - tbf_ref(k + nx, i__), abs(d__1))
			 > eps) {
		    cyclic = FALSE_;
		}
	    }
	    k = ny * ktrtor_1.tnx[i__ - 1];
	    i__2 = ktrtor_1.tnx[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		if ((d__1 = tbf_ref(j, i__) - tbf_ref(j + k, i__), abs(d__1)) 
			> eps) {
		    cyclic = FALSE_;
		}
	    }
	}

/*     spline fit the derivatives about the first torsion */

	i__2 = ktrtor_1.tnx[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    tmp1[j - 1] = ttx_ref(j, i__);
	}
	m = 0;
	i__2 = ktrtor_1.tny[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    i__3 = ktrtor_1.tnx[i__ - 1];
	    for (k = 1; k <= i__3; ++k) {
		tmp2[k - 1] = tbf_ref(m + k, i__);
	    }
	    if (cyclic) {
		cspline_(&nx, &c__30, tmp1, tmp2, bs, cs, ds, tmp3, tmp4, 
			tmp5, tmp6, tmp7);
	    } else {
		nspline_(&nx, &c__30, tmp1, tmp2, bs, cs, tmp3, tmp4, tmp5, 
			tmp6, tmp7);
	    }
	    i__3 = ktrtor_1.tnx[i__ - 1];
	    for (k = 1; k <= i__3; ++k) {
		tbx_ref(m + k, i__) = bs[k - 1];
	    }
	    m += ktrtor_1.tnx[i__ - 1];
	}

/*     spline fit the derivatives about the second torsion */

	i__2 = ktrtor_1.tny[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    tmp1[j - 1] = tty_ref(j, i__);
	}
	m = 1;
	i__2 = ktrtor_1.tnx[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    i__3 = ktrtor_1.tny[i__ - 1];
	    for (k = 1; k <= i__3; ++k) {
		tmp2[k - 1] = tbf_ref(m + (k - 1) * ktrtor_1.tnx[i__ - 1], 
			i__);
	    }
	    if (cyclic) {
		cspline_(&ny, &c__30, tmp1, tmp2, bs, cs, ds, tmp3, tmp4, 
			tmp5, tmp6, tmp7);
	    } else {
		nspline_(&ny, &c__30, tmp1, tmp2, bs, cs, tmp3, tmp4, tmp5, 
			tmp6, tmp7);
	    }
	    i__3 = ktrtor_1.tny[i__ - 1];
	    for (k = 1; k <= i__3; ++k) {
		tby_ref(m + (k - 1) * ktrtor_1.tnx[i__ - 1], i__) = bs[k - 1];
	    }
	    ++m;
	}

/*     spline fit the cross derivatives about both torsions */

	m = 1;
	i__2 = ktrtor_1.tnx[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    i__3 = ktrtor_1.tny[i__ - 1];
	    for (k = 1; k <= i__3; ++k) {
		tmp2[k - 1] = tbx_ref(m + (k - 1) * ktrtor_1.tnx[i__ - 1], 
			i__);
	    }
	    if (cyclic) {
		cspline_(&ny, &c__30, tmp1, tmp2, bs, cs, ds, tmp3, tmp4, 
			tmp5, tmp6, tmp7);
	    } else {
		nspline_(&ny, &c__30, tmp1, tmp2, bs, cs, tmp3, tmp4, tmp5, 
			tmp6, tmp7);
	    }
	    i__3 = ktrtor_1.tny[i__ - 1];
	    for (k = 1; k <= i__3; ++k) {
		tbxy_ref(m + (k - 1) * ktrtor_1.tnx[i__ - 1], i__) = bs[k - 1]
			;
	    }
	    ++m;
	}
    }

/*     assign torsion-torsion parameters for each bitorsion */

    tortor_1.ntortor = 0;
    i__1 = bitor_1.nbitor;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ibitor_ref(1, i__);
	ib = ibitor_ref(2, i__);
	ic = ibitor_ref(3, i__);
	id = ibitor_ref(4, i__);
	ie = ibitor_ref(5, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	itd = atmtyp_1.class__[id - 1];
	ite = atmtyp_1.class__[ie - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	numeral_(&itd, pd, &size, (ftnlen)4);
	numeral_(&ite, pe, &size, (ftnlen)4);
/* Writing concatenation */
	i__4[0] = 4, a__1[0] = pa;
	i__4[1] = 4, a__1[1] = pb;
	i__4[2] = 4, a__1[2] = pc;
	i__4[3] = 4, a__1[3] = pd;
	i__4[4] = 4, a__1[4] = pe;
	s_cat(pt1, a__1, i__4, &c__5, (ftnlen)20);
/* Writing concatenation */
	i__4[0] = 4, a__1[0] = pe;
	i__4[1] = 4, a__1[1] = pd;
	i__4[2] = 4, a__1[2] = pc;
	i__4[3] = 4, a__1[3] = pb;
	i__4[4] = 4, a__1[4] = pa;
	s_cat(pt2, a__1, i__4, &c__5, (ftnlen)20);

/*     find parameters for this torsion-torsion interaction */

	i__2 = ntt;
	for (j = 1; j <= i__2; ++j) {
	    if (s_cmp(ktt_ref(0, j), pt1, (ftnlen)20, (ftnlen)20) == 0) {
		++tortor_1.ntortor;
		itt_ref(1, tortor_1.ntortor) = i__;
		itt_ref(2, tortor_1.ntortor) = j;
		itt_ref(3, tortor_1.ntortor) = 1;
		goto L60;
	    } else if (s_cmp(ktt_ref(0, j), pt2, (ftnlen)20, (ftnlen)20) == 0)
		     {
		++tortor_1.ntortor;
		itt_ref(1, tortor_1.ntortor) = i__;
		itt_ref(2, tortor_1.ntortor) = j;
		itt_ref(3, tortor_1.ntortor) = -1;
		goto L60;
	    }
	}
L60:
	;
    }

/*     turn off the torsion-torsion potential if it is not used */

    if (tortor_1.ntortor == 0) {
	potent_1.use_tortor__ = FALSE_;
    }
    return 0;
} /* ktortor_ */

#undef keyline_ref
#undef ibitor_ref
#undef tbxy_ref
#undef tty_ref
#undef ttx_ref
#undef ktt_ref
#undef itt_ref
#undef tby_ref
#undef tbx_ref
#undef tbf_ref


