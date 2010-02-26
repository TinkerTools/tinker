/* initrot.f -- translated by f2c (version 20050501).
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
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

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
    doublereal dihed[1000];
    integer nomega, iomega[2000]	/* was [2][1000] */, zline[1000];
} omega_;

#define omega_1 omega_

struct {
    integer nrot, rot[25000];
    logical use_short__;
} rotate_;

#define rotate_1 rotate_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal zbond[25000], zang[25000], ztors[25000];
    integer iz[100000]	/* was [4][25000] */;
} zcoord_;

#define zcoord_1 zcoord_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine initrot  --  set bonds for dihedral rotation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "initrot" sets the torsional angles which are to be rotated */
/*     in subsequent computation, by default automatically selects */
/*     all rotatable single bonds; assumes internal coordinates have */
/*     already been setup */


/* Subroutine */ int initrot_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Selection of Torsional Angles for Rotati"
	    "on :\002,//,\002    0  - Automatic Selection of Torsional Angle"
	    "s\002,/,\002    1  - Manual Selection of Angles to Rotate\002,/"
	    ",\002    2  - Manual Selection of Angles to Freeze\002,//,\002 E"
	    "nter the Method of Choice [0] :  \002,$)";
    static char fmt_30[] = "(i10)";
    static char fmt_40[] = "(/,\002 Enter Atoms in Rotatable Bond\002,i5,"
	    "\002 :  \002,$)";
    static char fmt_50[] = "(a120)";
    static char fmt_60[] = "(/,\002 INITROT  --  Bond between Atoms\002,2i6"
	    ",\002 is not Rotatable\002)";
    static char fmt_90[] = "(/,\002 Enter Atoms in Frozen Bond\002,i5,\002 :"
	    "  \002,$)";
    static char fmt_100[] = "(a120)";
    static char fmt_130[] = "(/,\002 INITROT  --  Rotation about\002,2i6,"
	    "\002 occurs more than once in Z-matrix\002)";
    static char fmt_140[] = "(/,\002 List of Active Atoms for Torsional\002"
	    ",\002 Calculations :\002,/)";
    static char fmt_150[] = "(3x,10i7)";
    static char fmt_160[] = "(/,\002 INITROT  --  No Torsions for Subsequen"
	    "t\002,\002 Computation\002)";
    static char fmt_170[] = "(/,\002 INITROT  --  Too many Torsions;\002,"
	    "\002 Increase MAXROT\002)";
    static char fmt_180[] = "(/,\002 Number of Torsions Used in Derivativ"
	    "e\002,\002 Computation :\002,i6)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);

    /* Local variables */
    extern logical rotcheck_(integer *, integer *);
    static integer i__, j, j1, j2, mode, list[25000], bond1, bond2;
    extern /* Subroutine */ int fatal_(void);
    static integer iring, nlist;
    static logical exist, query;
    static integer ifixed[100000]	/* was [2][50000] */, nfixed;
    static char record[120];
    static logical rotate;
    static char string[120];
    static integer attach1, attach2;
    extern /* Subroutine */ int chkring_(integer *, integer *, integer *, 
	    integer *, integer *), nextarg_(char *, logical *, ftnlen), 
	    rotlist_(integer *, integer *);

    /* Fortran I/O blocks */
    static icilist io___5 = { 1, string, 1, 0, 120, 1 };
    static cilist io___6 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___13 = { 1, record, 1, 0, 120, 1 };
    static cilist io___19 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_100, 0 };
    static icilist io___24 = { 1, record, 1, 0, 120, 1 };
    static cilist io___28 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_180, 0 };



#define iz_ref(a_1,a_2) zcoord_1.iz[(a_2)*4 + a_1 - 5]
#define iomega_ref(a_1,a_2) omega_1.iomega[(a_2)*2 + a_1 - 3]
#define ifixed_ref(a_1,a_2) ifixed[(a_2)*2 + a_1 - 3]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  rotate.i  --  molecule partitions for rotation of a bond  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     nrot        total number of atoms moving when bond rotates */
/*     rot         atom numbers of atoms moving when bond rotates */
/*     use_short   logical flag governing use of shortest atom list */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     zbond   bond length used to define each Z-matrix atom */
/*     zang    bond angle used to define each Z-matrix atom */
/*     ztors   angle or torsion used to define Z-matrix atom */
/*     iz      defining atom numbers for each Z-matrix atom */




/*     initialize the number of rotatable torsional angles */

    omega_1.nomega = 0;

/*     use shortest rotlist if there is no absolute coordinate frame */

    rotate_1.use_short__ = TRUE_;
    if (group_1.use_group__) {
	rotate_1.use_short__ = FALSE_;
    }
    if (kgeoms_1.npfix != 0) {
	rotate_1.use_short__ = FALSE_;
    }

/*     choose automatic or manual selection of torsional angles */

    mode = 0;
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___5);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&mode, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	query = FALSE_;
    }
L10:
    if (query) {
	io___6.ciunit = iounit_1.iout;
	s_wsfe(&io___6);
	e_wsfe();
	io___7.ciunit = iounit_1.input;
	s_rsfe(&io___7);
	do_fio(&c__1, (char *)&mode, (ftnlen)sizeof(integer));
	e_rsfe();
    }
    if (mode != 1 && mode != 2) {
	mode = 0;
    }

/*     manual selection of the torsional angles to be rotated */

    if (mode == 1) {
	while(omega_1.nomega < 1000) {
	    ++omega_1.nomega;
	    j1 = 0;
	    j2 = 0;
	    io___10.ciunit = iounit_1.iout;
	    s_wsfe(&io___10);
	    do_fio(&c__1, (char *)&omega_1.nomega, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___11.ciunit = iounit_1.input;
	    s_rsfe(&io___11);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    i__1 = s_rsli(&io___13);
	    if (i__1 != 0) {
		goto L80;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&j1, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L80;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&j2, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L80;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L80;
	    }
	    if (j1 == 0 && j2 == 0) {
		goto L80;
	    }
	    i__1 = atoms_1.n;
	    for (i__ = 4; i__ <= i__1; ++i__) {
		if (iz_ref(4, i__) == 0) {
		    bond1 = iz_ref(1, i__);
		    bond2 = iz_ref(2, i__);
		    attach1 = couple_1.n12[bond1 - 1];
		    attach2 = couple_1.n12[bond2 - 1];
		    if (attach1 > 1 && attach2 > 1) {
			if (bond1 == j1 && bond2 == j2 || bond1 == j2 && 
				bond2 == j1) {
			    if (rotcheck_(&bond1, &bond2)) {
				iomega_ref(1, omega_1.nomega) = bond1;
				iomega_ref(2, omega_1.nomega) = bond2;
				omega_1.dihed[omega_1.nomega - 1] = 
					zcoord_1.ztors[i__ - 1] / 
					57.29577951308232088;
				omega_1.zline[omega_1.nomega - 1] = i__;
				goto L70;
			    }
			}
		    }
		}
	    }
	    --omega_1.nomega;
	    io___19.ciunit = iounit_1.iout;
	    s_wsfe(&io___19);
	    do_fio(&c__1, (char *)&j1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&j2, (ftnlen)sizeof(integer));
	    e_wsfe();
L70:
	    ;
	}
L80:
	--omega_1.nomega;
    }

/*     manual selection of the torsional angles to be frozen */

    nfixed = 0;
    if (mode == 2) {
	for (i__ = 1; i__ <= 1000; ++i__) {
	    ifixed_ref(1, i__) = 0;
	    ifixed_ref(2, i__) = 0;
	    io___22.ciunit = iounit_1.iout;
	    s_wsfe(&io___22);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___23.ciunit = iounit_1.input;
	    s_rsfe(&io___23);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    i__1 = s_rsli(&io___24);
	    if (i__1 != 0) {
		goto L110;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ifixed_ref(1, i__), (ftnlen)
		    sizeof(integer));
	    if (i__1 != 0) {
		goto L110;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ifixed_ref(2, i__), (ftnlen)
		    sizeof(integer));
	    if (i__1 != 0) {
		goto L110;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L110;
	    }
	    if (ifixed_ref(1, i__) == 0 && ifixed_ref(2, i__) == 0) {
		goto L110;
	    }
	    ++nfixed;
	}
L110:
	;
    }

/*     perform the automatic selection of torsional angles to rotate */

    if (mode == 0 || mode == 2) {
	i__1 = atoms_1.n;
	for (i__ = 4; i__ <= i__1; ++i__) {
	    if (iz_ref(4, i__) == 0) {
		rotate = TRUE_;
		bond1 = iz_ref(1, i__);
		bond2 = iz_ref(2, i__);

/*     do not rotate a bond if either bonded atom is univalent */

		attach1 = couple_1.n12[bond1 - 1];
		attach2 = couple_1.n12[bond2 - 1];
		if (attach1 <= 1 || attach2 <= 1) {
		    rotate = FALSE_;
		}

/*     do not rotate a bond contained within a small ring */

		iring = 0;
		chkring_(&iring, &bond1, &bond2, &c__0, &c__0);
		if (iring != 0) {
		    rotate = FALSE_;
		}

/*     do not rotate bonds explicitly frozen by the user */

		if (mode == 2 && rotate) {
		    i__2 = nfixed;
		    for (j = 1; j <= i__2; ++j) {
			j1 = ifixed_ref(1, j);
			j2 = ifixed_ref(2, j);
			if (bond1 == j1 && bond2 == j2 || bond1 == j2 && 
				bond2 == j1) {
			    rotate = FALSE_;
			    goto L120;
			}
		    }
		}
L120:

/*     do not rotate bonds with inactive atoms on both sides */

		if (rotate) {
		    if (! rotcheck_(&bond1, &bond2)) {
			rotate = FALSE_;
		    }
		}

/*     check for possible duplication of rotatable bonds */

		if (rotate) {
		    i__2 = omega_1.nomega;
		    for (j = 1; j <= i__2; ++j) {
			j1 = iomega_ref(1, j);
			j2 = iomega_ref(2, j);
			if (bond1 == j1 && bond2 == j2 || bond1 == j2 && 
				bond2 == j1) {
			    io___28.ciunit = iounit_1.iout;
			    s_wsfe(&io___28);
			    do_fio(&c__1, (char *)&bond1, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&bond2, (ftnlen)sizeof(
				    integer));
			    e_wsfe();
			    fatal_();
			}
		    }
		    ++omega_1.nomega;
		    iomega_ref(1, omega_1.nomega) = bond1;
		    iomega_ref(2, omega_1.nomega) = bond2;
		    omega_1.dihed[omega_1.nomega - 1] = zcoord_1.ztors[i__ - 
			    1] / 57.29577951308232088;
		    omega_1.zline[omega_1.nomega - 1] = i__;
		}
	    }
	}
    }

/*     make inactive the atoms not rotatable via any torsion */

    if (usage_1.nuse == atoms_1.n) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    usage_1.use[i__ - 1] = FALSE_;
	}
	i__1 = omega_1.nomega;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    bond1 = iomega_ref(1, i__);
	    bond2 = iomega_ref(2, i__);
	    rotlist_(&bond1, &bond2);
	    i__2 = rotate_1.nrot;
	    for (j = 1; j <= i__2; ++j) {
		usage_1.use[rotate_1.rot[j - 1] - 1] = TRUE_;
	    }
	}
	usage_1.nuse = 0;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (usage_1.use[i__ - 1]) {
		++usage_1.nuse;
	    }
	}
	if (inform_1.debug && usage_1.nuse > 0 && usage_1.nuse < atoms_1.n) {
	    nlist = 0;
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (usage_1.use[i__ - 1]) {
		    ++nlist;
		    list[nlist - 1] = i__;
		}
	    }
	    io___31.ciunit = iounit_1.iout;
	    s_wsfe(&io___31);
	    e_wsfe();
	    io___32.ciunit = iounit_1.iout;
	    s_wsfe(&io___32);
	    i__1 = nlist;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&list[i__ - 1], (ftnlen)sizeof(integer))
			;
	    }
	    e_wsfe();
	}
    }

/*     write out the number of rotatable torsions to be used */

    if (omega_1.nomega == 0) {
	io___33.ciunit = iounit_1.iout;
	s_wsfe(&io___33);
	e_wsfe();
	fatal_();
    } else if (omega_1.nomega > 1000) {
	io___34.ciunit = iounit_1.iout;
	s_wsfe(&io___34);
	e_wsfe();
	fatal_();
    }
    io___35.ciunit = iounit_1.iout;
    s_wsfe(&io___35);
    do_fio(&c__1, (char *)&omega_1.nomega, (ftnlen)sizeof(integer));
    e_wsfe();
    return 0;
} /* initrot_ */

#undef ifixed_ref
#undef iomega_ref
#undef iz_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function rotcheck  --  check for fixed atoms across bond  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "rotcheck" tests a specified candidate rotatable bond for */
/*     the disallowed case where inactive atoms are found on both */
/*     sides of the candidate bond */


logical rotcheck_(integer *base, integer *partner)
{
    /* System generated locals */
    integer i__1;
    logical ret_val;

    /* Local variables */
    static integer i__;
    static logical list[25000], value;
    extern /* Subroutine */ int rotlist_(integer *, integer *);



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  rotate.i  --  molecule partitions for rotation of a bond  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     nrot        total number of atoms moving when bond rotates */
/*     rot         atom numbers of atoms moving when bond rotates */
/*     use_short   logical flag governing use of shortest atom list */




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




/*     initialize status and find atoms on short side of the bond */

    value = TRUE_;
    rotlist_(base, partner);

/*     rotation is allowed if all atoms on one side are active */

    i__1 = rotate_1.nrot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (! usage_1.use[rotate_1.rot[i__ - 1] - 1]) {
	    value = FALSE_;
	    goto L10;
	}
    }
L10:

/*     if short side had inactive atoms, check the other side */

    if (! value) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    list[i__ - 1] = TRUE_;
	}
	i__1 = rotate_1.nrot;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    list[rotate_1.rot[i__ - 1] - 1] = FALSE_;
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (list[i__ - 1] && ! usage_1.use[i__ - 1]) {
		goto L20;
	    }
	}
	value = TRUE_;
L20:
	;
    }

/*     set the final return value of the function */

    ret_val = value;
    return ret_val;
} /* rotcheck_ */

