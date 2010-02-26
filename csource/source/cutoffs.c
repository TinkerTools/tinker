/* cutoffs.f -- translated by f2c (version 20050501).
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
    doublereal hesscut;
} hescut_;

#define hescut_1 hescut_

struct {
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

struct {
    doublereal lbuffer, lbuf2, vbuf2, cbuf2, mbuf2;
    integer nvlst[25000], vlst[45000000]	/* was [1800][25000] */, 
	    nelst[25000], elst[30000000]	/* was [1200][25000] */;
    logical dovlst, doclst, domlst;
} neigh_;

#define neigh_1 neigh_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine cutoffs  --  set distance and Hessian cutoffs  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "cutoffs" initializes and stores spherical energy cutoff */
/*     distance windows, Hessian element and Ewald sum cutoffs, */
/*     and the pairwise neighbor generation method */


/* Subroutine */ int cutoffs_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);

    /* Local variables */
    static logical truncate;
    static integer i__;
    static doublereal big;
    static integer next;
    static doublereal value;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120], keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___8 = { 1, string, 1, 0, 120, 1 };
    static icilist io___9 = { 1, string, 1, 0, 120, 1 };
    static icilist io___10 = { 1, string, 1, 0, 120, 1 };
    static icilist io___12 = { 1, string, 1, 0, 120, 1 };
    static icilist io___13 = { 1, string, 1, 0, 120, 1 };
    static icilist io___14 = { 1, string, 1, 0, 120, 1 };
    static icilist io___15 = { 1, string, 1, 0, 120, 1 };
    static icilist io___16 = { 1, string, 1, 0, 120, 1 };
    static icilist io___17 = { 1, string, 1, 0, 120, 1 };
    static icilist io___18 = { 1, string, 1, 0, 120, 1 };
    static icilist io___19 = { 1, string, 1, 0, 120, 1 };
    static icilist io___20 = { 1, string, 1, 0, 120, 1 };
    static icilist io___21 = { 1, string, 1, 0, 120, 1 };



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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  hescut.i  --  cutoff value for Hessian matrix elements  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     hesscut   magnitude of smallest allowed Hessian element */




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




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  neigh.i  --  pairwise neighbor list indices and storage  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     lbuffer     width of the neighbor list buffer region */
/*     lbuf2       motion squared needed to trigger list rebuild */
/*     vbuf2       square of vdw cutoff plus neighbor list buffer */
/*     cbuf2       square of charge cutoff plus neighbor list buffer */
/*     mbuf2       square of multipole cutoff plus neighbor list buffer */
/*     nvlst       number of sites in list for each vdw site */
/*     vlst        site numbers in neighbor list of each vdw site */
/*     nelst       number of sites in list for each electrostatic site */
/*     elst        site numbers in list of each electrostatic site */
/*     dovlst      logical flag to rebuild vdw neighbor list */
/*     doclst      logical flag to rebuild charge neighbor list */
/*     domlst      logical flag to rebuild multipole neighbor list */




/*     set defaults for spherical energy cutoff distances */

    big = 1e12;
    if (bound_1.use_bounds__) {
	cutoff_1.vdwcut = 9.;
	cutoff_1.chgcut = 9.;
	cutoff_1.dplcut = 9.;
	cutoff_1.mpolecut = 9.;
    } else {
	cutoff_1.vdwcut = big;
	cutoff_1.chgcut = big;
	cutoff_1.dplcut = big;
	cutoff_1.mpolecut = big;
    }
    cutoff_1.ewaldcut = 7.;

/*     set defaults for tapering, Hessian cutoff and neighbor buffer */

    cutoff_1.vdwtaper = .9;
    cutoff_1.chgtaper = .65;
    cutoff_1.dpltaper = .75;
    cutoff_1.mpoletaper = .65;
    hescut_1.hesscut = 0.;
    neigh_1.lbuffer = 2.;

/*     set defaults for Ewald sum, tapering style and neighbor method */

    cutoff_1.use_ewald__ = FALSE_;
    truncate = FALSE_;
    cutoff_1.use_lights__ = FALSE_;
    cutoff_1.use_list__ = FALSE_;
    cutoff_1.use_vlist__ = FALSE_;
    cutoff_1.use_clist__ = FALSE_;
    cutoff_1.use_mlist__ = FALSE_;

/*     search the keywords for various cutoff parameters */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));

/*     get values related to use of Ewald summation */

	if (s_cmp(keyword, "EWALD ", (ftnlen)6, (ftnlen)6) == 0) {
	    cutoff_1.use_ewald__ = TRUE_;
	} else if (s_cmp(keyword, "EWALD-CUTOFF ", (ftnlen)13, (ftnlen)13) == 
		0) {
	    i__2 = s_rsli(&io___8);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&cutoff_1.ewaldcut, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }

/*     get values for the tapering style and neighbor method */

	} else if (s_cmp(keyword, "TRUNCATE ", (ftnlen)9, (ftnlen)9) == 0) {
	    truncate = TRUE_;
	} else if (s_cmp(keyword, "LIGHTS ", (ftnlen)7, (ftnlen)7) == 0) {
	    cutoff_1.use_lights__ = TRUE_;
	} else if (s_cmp(keyword, "NEIGHBOR-LIST ", (ftnlen)14, (ftnlen)14) ==
		 0) {
	    cutoff_1.use_list__ = TRUE_;
	    cutoff_1.use_vlist__ = TRUE_;
	    cutoff_1.use_clist__ = TRUE_;
	    cutoff_1.use_mlist__ = TRUE_;
	} else if (s_cmp(keyword, "VDW-LIST ", (ftnlen)9, (ftnlen)9) == 0) {
	    cutoff_1.use_list__ = TRUE_;
	    cutoff_1.use_vlist__ = TRUE_;
	} else if (s_cmp(keyword, "CHG-LIST ", (ftnlen)9, (ftnlen)9) == 0) {
	    cutoff_1.use_list__ = TRUE_;
	    cutoff_1.use_clist__ = TRUE_;
	} else if (s_cmp(keyword, "MPOLE-LIST ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    cutoff_1.use_list__ = TRUE_;
	    cutoff_1.use_mlist__ = TRUE_;

/*     get cutoff for the magnitude of Hessian elements */

	} else if (s_cmp(keyword, "HESS-CUTOFF ", (ftnlen)12, (ftnlen)12) == 
		0) {
	    i__2 = s_rsli(&io___9);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&hescut_1.hesscut, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }

/*     get the cutoff radii for potential energy functions */

	} else if (s_cmp(keyword, "CUTOFF ", (ftnlen)7, (ftnlen)7) == 0) {
	    i__2 = s_rsli(&io___10);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&value, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	    cutoff_1.vdwcut = value;
	    cutoff_1.chgcut = value;
	    cutoff_1.dplcut = value;
	    cutoff_1.mpolecut = value;
	    cutoff_1.ewaldcut = value;
	} else if (s_cmp(keyword, "VDW-CUTOFF ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    i__2 = s_rsli(&io___12);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&cutoff_1.vdwcut, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "CHG-CUTOFF ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    i__2 = s_rsli(&io___13);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&cutoff_1.chgcut, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "DPL-CUTOFF ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    i__2 = s_rsli(&io___14);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&cutoff_1.dplcut, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "MPOLE-CUTOFF ", (ftnlen)13, (ftnlen)13) == 
		0) {
	    i__2 = s_rsli(&io___15);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&cutoff_1.mpolecut, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }

/*     get distance for initialization of energy switching */

	} else if (s_cmp(keyword, "TAPER ", (ftnlen)6, (ftnlen)6) == 0) {
	    i__2 = s_rsli(&io___16);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&value, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	    cutoff_1.vdwtaper = value;
	    cutoff_1.chgtaper = value;
	    cutoff_1.dpltaper = value;
	    cutoff_1.mpoletaper = value;
	} else if (s_cmp(keyword, "VDW-TAPER ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    i__2 = s_rsli(&io___17);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&cutoff_1.vdwtaper, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "CHG-TAPER ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    i__2 = s_rsli(&io___18);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&cutoff_1.chgtaper, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "DPL-TAPER ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    i__2 = s_rsli(&io___19);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&cutoff_1.dpltaper, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "MPOLE-TAPER ", (ftnlen)12, (ftnlen)12) == 
		0) {
	    i__2 = s_rsli(&io___20);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&cutoff_1.mpoletaper, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }

/*     get buffer width for pairwise nonbonded neighbor lists */

	} else if (s_cmp(keyword, "LIST-BUFFER ", (ftnlen)12, (ftnlen)12) == 
		0) {
	    i__2 = s_rsli(&io___21);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&neigh_1.lbuffer, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	}
L10:
	;
    }

/*     convert any distance percentages to absolute distances */

    if (cutoff_1.vdwtaper < 1.) {
	cutoff_1.vdwtaper *= cutoff_1.vdwcut;
    }
    if (cutoff_1.chgtaper < 1.) {
	cutoff_1.chgtaper *= cutoff_1.chgcut;
    }
    if (cutoff_1.dpltaper < 1.) {
	cutoff_1.dpltaper *= cutoff_1.dplcut;
    }
    if (cutoff_1.mpoletaper < 1.) {
	cutoff_1.mpoletaper *= cutoff_1.mpolecut;
    }

/*     apply truncation cutoffs if they were requested */

    if (truncate) {
	cutoff_1.vdwtaper = big;
	cutoff_1.chgtaper = big;
	cutoff_1.dpltaper = big;
	cutoff_1.mpoletaper = big;
    }

/*     set buffer region limits for pairwise neighbor lists */

/* Computing 2nd power */
    d__1 = neigh_1.lbuffer * .5;
    neigh_1.lbuf2 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = cutoff_1.vdwcut + neigh_1.lbuffer;
    neigh_1.vbuf2 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = cutoff_1.chgcut + neigh_1.lbuffer;
    neigh_1.cbuf2 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = cutoff_1.mpolecut + neigh_1.lbuffer;
    neigh_1.mbuf2 = d__1 * d__1;
    if (cutoff_1.use_ewald__) {
/* Computing 2nd power */
	d__1 = cutoff_1.ewaldcut + neigh_1.lbuffer;
	neigh_1.cbuf2 = d__1 * d__1;
/* Computing 2nd power */
	d__1 = cutoff_1.ewaldcut + neigh_1.lbuffer;
	neigh_1.mbuf2 = d__1 * d__1;
    }
    return 0;
} /* cutoffs_ */

#undef keyline_ref


