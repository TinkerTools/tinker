/* cluster.f -- translated by f2c (version 20050501).
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
    doublereal vdwcut, chgcut, dplcut, mpolecut, vdwtaper, chgtaper, dpltaper,
	     mpoletaper, ewaldcut;
    logical use_ewald__, use_lights__, use_list__, use_vlist__, use_clist__, 
	    use_mlist__;
} cutoff_;

#define cutoff_1 cutoff_

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
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

struct {
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine cluster  --  set user-defined groups of atoms  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "cluster" gets the partitioning of the system into groups */
/*     and stores a list of the group to which each atom belongs */


/* Subroutine */ int cluster_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 CLUSTER  --  Too many Atom Groups;\002"
	    ",\002 Increase MAXGRP\002)";
    static char fmt_30[] = "(/,\002 CLUSTER  --  Too many Atom Groups;\002"
	    ",\002 Increase MAXGRP\002)";
    static char fmt_50[] = "(/,\002 List of Atoms Contained in Group\002,i3"
	    ",\002 :\002,/)";
    static char fmt_60[] = "(3x,10i7)";
    static char fmt_70[] = "(/,\002 Active Sets of Intra- and InterGrou"
	    "p\002,\002 Interactions :\002,//,11x,\002Groups\002,15x,\002Typ"
	    "e\002,14x,\002Weight\002,/)";
    static char fmt_80[] = "(5x,2i6,12x,\002IntraGroup\002,5x,f12.4)";
    static char fmt_90[] = "(5x,2i6,12x,\002InterGroup\002,5x,f12.4)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), s_rsli(icilist *), do_lio(integer *, integer *, char *, 
	    ftnlen), e_rsli(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, ga, gb;
    static doublereal wg;
    static integer gnum, size[1001], list[25000], next;
    extern /* Subroutine */ int sort_(integer *, integer *), sort3_(integer *,
	     integer *, integer *), fatal_(void);
    static logical header;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int getnumb_(char *, integer *, integer *, ftnlen)
	    , cutoffs_(void);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static cilist io___8 = { 0, 0, 0, fmt_10, 0 };
    static icilist io___10 = { 1, string, 1, 0, 120, 1 };
    static cilist io___12 = { 0, 0, 0, fmt_30, 0 };
    static icilist io___16 = { 1, string, 1, 0, 120, 1 };
    static cilist io___18 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_90, 0 };



#define imol_ref(a_1,a_2) molcul_1.imol[(a_2)*2 + a_1 - 3]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]
#define wgrp_ref(a_1,a_2) group_1.wgrp[(a_2)*1001 + a_1 - 0]
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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cutoff.i  --  cutoff distances for energy interactions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     vdwcut      cutoff distance for van der Waals interactions */
/*     chgcut      cutoff distance for charge-charge interactions */
/*     dplcut      cutoff distance for dipole-dipole interactions */
/*     mpolecut    cutoff distance for atomic multipole interactions */
/*     vdwtaper    distance at which van der Waals switching begins */
/*     chgtaper    distance at which charge-charge switching begins */
/*     dpltaper    distance at which dipole-dipole switching begins */
/*     mpoletaper  distance at which atomic multipole switching begins */
/*     ewaldcut    cutoff distance for direct space Ewald summation */
/*     use_ewald   logical flag governing use of Ewald summation */
/*     use_lights  logical flag governing use of method of lights */
/*     use_list    logical flag governing use of any neighbor lists */
/*     use_vlist   logical flag governing use of vdw neighbor lists */
/*     use_clist   logical flag governing use of charge neighbor lists */
/*     use_mlist   logical flag governing use of multipole neighbor lists */




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
/*     ##  molcul.i  --  individual molecules within current system  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     molmass   molecular weight for each molecule in the system */
/*     totmass   total weight of all the molecules in the system */
/*     nmol      total number of separate molecules in the system */
/*     kmol      contiguous list of the atoms in each molecule */
/*     imol      first and last atom of each molecule in the list */
/*     molcule   number of the molecule to which each atom belongs */




/*     set defaults for the group atom list and weight options */

    group_1.use_group__ = FALSE_;
    group_1.use_intra__ = FALSE_;
    group_1.use_inter__ = FALSE_;
    group_1.ngrp = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	group_1.kgrp[i__ - 1] = 0;
	group_1.grplist[i__ - 1] = 0;
    }
    for (i__ = 1; i__ <= 1000; ++i__) {
	igrp_ref(1, i__) = 1;
	igrp_ref(2, i__) = 0;
    }
    for (i__ = 0; i__ <= 1000; ++i__) {
	for (j = 0; j <= 1000; ++j) {
	    wgrp_ref(j, i__) = 1.;
	}
    }

/*     get any keywords containing atom group definitions */

    i__1 = keys_1.nkey;
    for (j = 1; j <= i__1; ++j) {
	next = 1;
	s_copy(record, keyline_ref(0, j), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "GROUP ", (ftnlen)6, (ftnlen)6) == 0) {
	    group_1.use_group__ = TRUE_;
	    gnum = 0;
	    i__2 = atoms_1.n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		list[i__ - 1] = 0;
	    }
	    getnumb_(record, &gnum, &next, (ftnlen)120);
	    if (gnum > 1000) {
		io___8.ciunit = iounit_1.iout;
		s_wsfe(&io___8);
		e_wsfe();
		fatal_();
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___10);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__3 = atoms_1.n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__2 = do_lio(&c__3, &c__1, (char *)&list[i__ - 1], (ftnlen)
			sizeof(integer));
		if (i__2 != 0) {
		    goto L20;
		}
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
L20:
	    i__ = 1;
	    while(list[i__ - 1] != 0) {
		if (list[i__ - 1] > 0) {
		    group_1.grplist[list[i__ - 1] - 1] = gnum;
		    ++i__;
		} else {
		    i__4 = (i__3 = list[i__], abs(i__3));
		    for (k = (i__2 = list[i__ - 1], abs(i__2)); k <= i__4; 
			    ++k) {
			group_1.grplist[k - 1] = gnum;
		    }
		    i__ += 2;
		}
	    }

/*     get any keywords with weights for group interactions */

	} else if (s_cmp(keyword, "GROUP-MOLECULE ", (ftnlen)15, (ftnlen)15) 
		== 0) {
	    group_1.use_group__ = TRUE_;
	    group_1.use_inter__ = TRUE_;
	    group_1.use_intra__ = FALSE_;
	    if (molcul_1.nmol > 1000) {
		io___12.ciunit = iounit_1.iout;
		s_wsfe(&io___12);
		e_wsfe();
		fatal_();
	    }
	    i__4 = molcul_1.nmol;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__2 = imol_ref(2, i__);
		for (k = imol_ref(1, i__); k <= i__2; ++k) {
		    group_1.grplist[molcul_1.kmol[k - 1] - 1] = i__;
		}
	    }

/*     get any keywords with weights for group interactions */

	} else if (s_cmp(keyword, "GROUP-SELECT ", (ftnlen)13, (ftnlen)13) == 
		0) {
	    ga = 0;
	    gb = 0;
	    wg = -1.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__4 = s_rsli(&io___16);
	    if (i__4 != 0) {
		goto L40;
	    }
	    i__4 = do_lio(&c__3, &c__1, (char *)&ga, (ftnlen)sizeof(integer));
	    if (i__4 != 0) {
		goto L40;
	    }
	    i__4 = do_lio(&c__3, &c__1, (char *)&gb, (ftnlen)sizeof(integer));
	    if (i__4 != 0) {
		goto L40;
	    }
	    i__4 = do_lio(&c__5, &c__1, (char *)&wg, (ftnlen)sizeof(
		    doublereal));
	    if (i__4 != 0) {
		goto L40;
	    }
	    i__4 = e_rsli();
	    if (i__4 != 0) {
		goto L40;
	    }
L40:
	    if (wg < 0.) {
		wg = 1.;
	    }
	    wgrp_ref(ga, gb) = wg;
	    wgrp_ref(gb, ga) = wg;
	    group_1.use_inter__ = FALSE_;

/*     get keywords to select common sets of group interactions */

	} else if (s_cmp(keyword, "GROUP-INTRA ", (ftnlen)12, (ftnlen)12) == 
		0) {
	    group_1.use_intra__ = TRUE_;
	    group_1.use_inter__ = FALSE_;
	} else if (s_cmp(keyword, "GROUP-INTER ", (ftnlen)12, (ftnlen)12) == 
		0) {
	    group_1.use_inter__ = TRUE_;
	    group_1.use_intra__ = FALSE_;
	}
    }

/*     pack atoms of each group into a contiguous indexed list */

    if (group_1.use_group__) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    list[i__ - 1] = group_1.grplist[i__ - 1];
	}
	sort3_(&atoms_1.n, list, group_1.kgrp);

/*     find the first and last atom in each of the groups */

	k = list[0];
	igrp_ref(1, k) = 1;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j = list[i__ - 1];
	    if (j != k) {
		igrp_ref(2, k) = i__ - 1;
		igrp_ref(1, j) = i__;
		k = j;
	    }
	    group_1.ngrp = max(j,group_1.ngrp);
	}
	igrp_ref(2, j) = atoms_1.n;

/*     sort the list of atoms in each group by atom number */

	for (i__ = 0; i__ <= 1000; ++i__) {
	    size[i__] = igrp_ref(2, i__) - igrp_ref(1, i__) + 1;
	    sort_(&size[i__], &group_1.kgrp[igrp_ref(1, i__) - 1]);
	}
    }

/*     use only intragroup or intergroup interactions if selected */

    if (group_1.use_intra__) {
	i__1 = group_1.ngrp;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    i__4 = group_1.ngrp;
	    for (j = 0; j <= i__4; ++j) {
		wgrp_ref(j, i__) = 0.;
	    }
	    wgrp_ref(i__, i__) = 1.;
	}
    }
    if (group_1.use_inter__) {
	i__1 = group_1.ngrp;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    i__4 = group_1.ngrp;
	    for (j = 0; j <= i__4; ++j) {
		wgrp_ref(j, i__) = 1.;
	    }
	    wgrp_ref(i__, i__) = 0.;
	}
    }

/*     disable consideration of interactions with any empty groups */

    i__1 = group_1.ngrp;
    for (i__ = 0; i__ <= i__1; ++i__) {
	if (size[i__] == 0) {
	    i__4 = group_1.ngrp;
	    for (j = 0; j <= i__4; ++j) {
		wgrp_ref(j, i__) = 0.;
		wgrp_ref(i__, j) = 0.;
	    }
	}
    }

/*     turn off bounds and replicas for intragroup calculations */

    if (group_1.use_intra__) {
	bound_1.use_bounds__ = FALSE_;
	bound_1.use_replica__ = FALSE_;
	cutoffs_();
    }

/*     compute the total mass of all atoms in each group */

    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	group_1.grpmass[i__ - 1] = 0.;
	i__4 = igrp_ref(2, i__);
	for (j = igrp_ref(1, i__); j <= i__4; ++j) {
	    group_1.grpmass[i__ - 1] += atmtyp_1.mass[group_1.kgrp[j - 1] - 1]
		    ;
	}
    }

/*     output the final list of atoms in each group */

    if (inform_1.debug && group_1.use_group__) {
	i__1 = group_1.ngrp;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (size[i__] != 0) {
		io___18.ciunit = iounit_1.iout;
		s_wsfe(&io___18);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		e_wsfe();
		io___19.ciunit = iounit_1.iout;
		s_wsfe(&io___19);
		i__4 = igrp_ref(2, i__);
		for (j = igrp_ref(1, i__); j <= i__4; ++j) {
		    do_fio(&c__1, (char *)&group_1.kgrp[j - 1], (ftnlen)
			    sizeof(integer));
		}
		e_wsfe();
	    }
	}
    }

/*     output the weights for intragroup and intergroup interactions */

    if (inform_1.debug && group_1.use_group__) {
	header = TRUE_;
	i__1 = group_1.ngrp;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__4 = group_1.ngrp;
	    for (j = i__; j <= i__4; ++j) {
		if (wgrp_ref(j, i__) != 0.) {
		    if (header) {
			header = FALSE_;
			io___21.ciunit = iounit_1.iout;
			s_wsfe(&io___21);
			e_wsfe();
		    }
		    if (i__ == j) {
			io___22.ciunit = iounit_1.iout;
			s_wsfe(&io___22);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&wgrp_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			e_wsfe();
		    } else {
			io___23.ciunit = iounit_1.iout;
			s_wsfe(&io___23);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&wgrp_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			e_wsfe();
		    }
		}
	    }
	}
    }
    return 0;
} /* cluster_ */

#undef keyline_ref
#undef wgrp_ref
#undef igrp_ref
#undef imol_ref


