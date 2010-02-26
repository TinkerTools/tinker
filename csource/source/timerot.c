/* timerot.f -- translated by f2c (version 20050501).
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
    doublereal vdwcut, chgcut, dplcut, mpolecut, vdwtaper, chgtaper, dpltaper,
	     mpoletaper, ewaldcut;
    logical use_ewald__, use_lights__, use_list__, use_vlist__, use_clist__, 
	    use_mlist__;
} cutoff_;

#define cutoff_1 cutoff_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  program timerot  --  timer for torsional energy terms  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "timerot" measures the CPU time required for file reading */
/*     and parameter assignment, potential energy computation, */
/*     energy and gradient over torsions, and torsional angle */
/*     Hessian matrix evaluation */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Enter Desired Number of Repetitions [1] "
	    ":  \002,$)";
    static char fmt_30[] = "(i10)";
    static char fmt_40[] = "(/,\002 Include Timing for Hessian Evaluations ["
	    "Y] :  \002,$)";
    static char fmt_50[] = "(a120)";
    static char fmt_60[] = "(/,\002 Computation Set-up :\002,f15.3,\002 Seco"
	    "nds\002)";
    static char fmt_70[] = "(/,\002 Potential Energy :  \002,f15.3,\002 Seco"
	    "nds for\002,i6,\002 Evaluations\002)";
    static char fmt_80[] = "(/,\002 Energy & Gradient : \002,f15.3,\002 Seco"
	    "nds for\002,i6,\002 Evaluations\002)";
    static char fmt_90[] = "(/,\002 Hessian Matrix :    \002,f15.3,\002 Seco"
	    "nds for\002,i6,\002 Evaluations\002)";
    static char fmt_100[] = "(/,\002 Potential Energy :  \002,f15.3,\002 Sec"
	    "onds for\002,i6,\002 Evaluations\002)";
    static char fmt_110[] = "(/,\002 Energy & Gradient : \002,f15.3,\002 Sec"
	    "onds for\002,i6,\002 Evaluations\002)";
    static char fmt_120[] = "(/,\002 Hessian Matrix :    \002,f15.3,\002 Sec"
	    "onds for\002,i6,\002 Evaluations\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    static integer i__;
    static logical dohessian;
    static doublereal hrot[1000000]	/* was [1000][1000] */;
    static integer next;
    extern /* Subroutine */ int final_(void);
    static doublereal value;
    static logical exist, query;
    extern /* Subroutine */ int getime_(doublereal *);
    static integer ncalls;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    extern doublereal energy_(void);
    extern /* Subroutine */ int getint_(void);
    static doublereal derivs[1000];
    static char answer[1], string[120];
    extern /* Subroutine */ int setime_(void), nblist_(void);
    static doublereal elapsed;
    extern /* Subroutine */ int initial_(void), gradrot_(doublereal *, 
	    doublereal *), nextarg_(char *, logical *, ftnlen), gettext_(char 
	    *, char *, integer *, ftnlen, ftnlen), hessrot_(char *, 
	    doublereal *, ftnlen), initrot_(void);

    /* Fortran I/O blocks */
    static icilist io___6 = { 1, string, 1, 0, 120, 1 };
    static cilist io___7 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_120, 0 };




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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     read in the molecular system to be timed */

    initial_();
    getint_();

/*     get the timing for setup of the calculation */

    setime_();
    mechanic_();
    initrot_();
    if (cutoff_1.use_list__) {
	nblist_();
    }
    getime_(&elapsed);

/*     get the number of calculation cycles to perform */

    ncalls = 0;
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___6);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
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
	io___7.ciunit = iounit_1.iout;
	s_wsfe(&io___7);
	e_wsfe();
	io___8.ciunit = iounit_1.input;
	s_rsfe(&io___8);
	do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
	e_rsfe();
    }
    if (ncalls == 0) {
	ncalls = 1;
    }

/*     decide whether to include timing of Hessian evaluations */

    dohessian = TRUE_;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___11.ciunit = iounit_1.iout;
	s_wsfe(&io___11);
	e_wsfe();
	io___12.ciunit = iounit_1.input;
	s_rsfe(&io___12);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'N') {
	dohessian = FALSE_;
    }

/*     print the time required for the computation setup */

    io___15.ciunit = iounit_1.iout;
    s_wsfe(&io___15);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     run the potential energy only timing experiment */

    setime_();
    i__1 = ncalls;
    for (i__ = 1; i__ <= i__1; ++i__) {
	value = energy_();
    }
    getime_(&elapsed);
    io___18.ciunit = iounit_1.iout;
    s_wsfe(&io___18);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
    e_wsfe();

/*     run the energy and gradient timing experiment */

    setime_();
    i__1 = ncalls;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gradrot_(&value, derivs);
    }
    getime_(&elapsed);
    io___20.ciunit = iounit_1.iout;
    s_wsfe(&io___20);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
    e_wsfe();

/*     run the Hessian matrix only timing experiment */

    if (dohessian) {
	setime_();
	i__1 = ncalls;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    hessrot_("FULL", hrot, (ftnlen)4);
	}
	getime_(&elapsed);
	io___22.ciunit = iounit_1.iout;
	s_wsfe(&io___22);
	do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
	e_wsfe();
    }

/*     repeat the potential energy only timing experiment */

    setime_();
    i__1 = ncalls;
    for (i__ = 1; i__ <= i__1; ++i__) {
	value = energy_();
    }
    getime_(&elapsed);
    io___23.ciunit = iounit_1.iout;
    s_wsfe(&io___23);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
    e_wsfe();

/*     repeat the energy and gradient timing experiment */

    setime_();
    i__1 = ncalls;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gradrot_(&value, derivs);
    }
    getime_(&elapsed);
    io___24.ciunit = iounit_1.iout;
    s_wsfe(&io___24);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
    e_wsfe();

/*     repeat the Hessian matrix only timing experiment */

    if (dohessian) {
	setime_();
	i__1 = ncalls;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    hessrot_("FULL", hrot, (ftnlen)4);
	}
	getime_(&elapsed);
	io___25.ciunit = iounit_1.iout;
	s_wsfe(&io___25);
	do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
	e_wsfe();
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

/* Main program alias */ int timerot_ () { MAIN__ (); return 0; }
