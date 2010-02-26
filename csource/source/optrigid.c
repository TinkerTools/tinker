/* optrigid.f -- translated by f2c (version 20050501).
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
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

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
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

struct {
    doublereal xrb[25000], yrb[25000], zrb[25000], rbc[6000]	/* was [6][
	    1000] */;
    logical use_rigid__;
} rigid_;

#define rigid_1 rigid_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  program optrigid  --  variable metric rigid body optimizer  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "optrigid" performs an energy minimization of rigid body atom */
/*     groups using an optimally conditioned variable metric method */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 OPTRIGID  --  Too many Parameters,\002"
	    ",\002 Increase the Value of MAXOPT\002)";
    static char fmt_40[] = "(/,\002 Enter RMS Gradient per Rigid Body Criter"
	    "ion\002,\002 [0.01] :  \002,$)";
    static char fmt_50[] = "(f20.0)";
    static char fmt_60[] = "(/,\002 OPTRIGID  --  Too many Parameters,\002"
	    ",\002 Increase the Value of MAXOPT\002)";
    static char fmt_70[] = "(/,\002 Final Function Value :\002,2x,f20.8,/"
	    ",\002 Final RMS Gradient :\002,4x,f20.8,/,\002 Final Gradient No"
	    "rm :\002,3x,f20.8)";
    static char fmt_80[] = "(/,\002 Final Function Value :\002,2x,f20.8,/"
	    ",\002 Final RMS Gradient :\002,4x,d20.8,/,\002 Final Gradient No"
	    "rm :\002,3x,d20.8)";
    static char fmt_90[] = "(/,\002 Final Function Value :\002,2x,f18.6,/"
	    ",\002 Final RMS Gradient :\002,4x,f18.6,/,\002 Final Gradient No"
	    "rm :\002,3x,f18.6)";
    static char fmt_100[] = "(/,\002 Final Function Value :\002,2x,f18.6,/"
	    ",\002 Final RMS Gradient :\002,4x,d18.6,/,\002 Final Gradient No"
	    "rm :\002,3x,d18.6)";
    static char fmt_110[] = "(/,\002 Final Function Value :\002,2x,f16.4,/"
	    ",\002 Final RMS Gradient :\002,4x,f16.4,/,\002 Final Gradient No"
	    "rm :\002,3x,f16.4)";
    static char fmt_120[] = "(/,\002 Final Function Value :\002,2x,f16.4,/"
	    ",\002 Final RMS Gradient :\002,4x,d16.4,/,\002 Final Gradient No"
	    "rm :\002,3x,d16.4)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2];
    doublereal d__1;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_rsfe(
	    cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);
    double sqrt(doublereal);
    integer f_rew(alist *);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    extern integer freeunit_(void);
    static integer i__, j;
    extern doublereal optrigid1_();
    static doublereal xx[1000];
    static integer imin;
    extern /* Subroutine */ int ocvm_(integer *, doublereal *, doublereal *, 
	    doublereal *, D_fp, U_fp);
    static integer nvar;
    static doublereal grms;
    static integer next;
    extern /* Subroutine */ int fatal_(void), final_(void);
    static doublereal gnorm;
    static logical exist;
    static char record[120];
    static doublereal grdmin;
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static doublereal derivs[6000]	/* was [6][1000] */;
    static char string[120];
    extern /* Subroutine */ int orient_(void), getxyz_(void), prtxyz_(integer 
	    *), gradrgd_(doublereal *, doublereal *);
    static char minfile[120];
    extern /* Subroutine */ int initial_(void), nextarg_(char *, logical *, 
	    ftnlen);
    static doublereal minimum;
    extern /* Subroutine */ int optsave_();
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };
    static icilist io___7 = { 1, string, 1, 0, 120, 1 };
    static icilist io___8 = { 1, string, 1, 0, 120, 1 };
    static icilist io___11 = { 1, string, 1, 0, 120, 1 };
    static cilist io___12 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_120, 0 };



#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]
#define derivs_ref(a_1,a_2) derivs[(a_2)*6 + a_1 - 7]
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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  files.i  --  name and number of current structure files  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     nprior     number of previously existing cycle files */
/*     ldir       length in characters of the directory name */
/*     leng       length in characters of the base filename */
/*     filename   base filename used by default for all files */
/*     outfile    output filename used for intermediate results */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  output.i  --  control of coordinate output file format  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     archive    logical flag to save structures in an archive */
/*     noversion  logical flag governing use of filename versions */
/*     overwrite  logical flag to overwrite intermediate files inplace */
/*     cyclesave  logical flag to mark use of numbered cycle files */
/*     coordtype  selects Cartesian, internal, rigid body or none */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  rigid.i  --  rigid body coordinates for atom groups  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     xrb         rigid body reference x-coordinate for each atom */
/*     yrb         rigid body reference y-coordinate for each atom */
/*     zrb         rigid body reference z-coordinate for each atom */
/*     rbc         current rigid body coordinates for each group */
/*     use_rigid   flag to mark use of rigid body coordinate system */




/*     set up the molecular mechanics calculation */

    initial_();
    getxyz_();
    mechanic_();

/*     set up the use of rigid body coordinate system */

    rigid_1.use_rigid__ = TRUE_;
    orient_();

/*     check for too many parameters to be optimized */

    if (group_1.ngrp * 6 > 1000) {
	io___1.ciunit = iounit_1.iout;
	s_wsfe(&io___1);
	e_wsfe();
	fatal_();
    }

/*     search the keywords for output frequency parameters */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "PRINTOUT ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___7);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&inform_1.iprint, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "WRITEOUT ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___8);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&inform_1.iwrite, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	}
L20:
	;
    }

/*     get termination criterion as RMS rigid body gradient */

    grdmin = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___11);
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&grdmin, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L30;
	}
    }
L30:
    if (grdmin <= 0.) {
	io___12.ciunit = iounit_1.iout;
	s_wsfe(&io___12);
	e_wsfe();
	io___13.ciunit = iounit_1.input;
	s_rsfe(&io___13);
	do_fio(&c__1, (char *)&grdmin, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (grdmin == 0.) {
	grdmin = .01;
    }

/*     write out a copy of coordinates for later update */

    imin = freeunit_();
/* Writing concatenation */
    i__3[0] = files_1.leng, a__1[0] = files_1.filename;
    i__3[1] = 4, a__1[1] = ".xyz";
    s_cat(minfile, a__1, i__3, &c__2, (ftnlen)120);
    version_(minfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = imin;
    o__1.ofnmlen = 120;
    o__1.ofnm = minfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    prtxyz_(&imin);
    cl__1.cerr = 0;
    cl__1.cunit = imin;
    cl__1.csta = 0;
    f_clos(&cl__1);
    s_copy(files_1.outfile, minfile, (ftnlen)120, (ftnlen)120);

/*     transfer rigid body coordinates to optimization parameters */

    nvar = 0;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    ++nvar;
	    xx[nvar - 1] = rbc_ref(j, i__);
	}
    }

/*     check for too many parameters to be optimized */

    if (nvar > 1000) {
	io___19.ciunit = iounit_1.iout;
	s_wsfe(&io___19);
	e_wsfe();
	fatal_();
    }

/*     make the call to the optimization routine */

    s_copy(output_1.coordtype, "RIGIDBODY", (ftnlen)9, (ftnlen)9);
    ocvm_(&nvar, xx, &minimum, &grdmin, (D_fp)optrigid1_, (U_fp)optsave_);

/*     transfer optimization parameters to rigid body coordinates */

    nvar = 0;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    ++nvar;
	    rbc_ref(j, i__) = xx[nvar - 1];
	}
    }

/*     compute the final function and RMS gradient values */

    gradrgd_(&minimum, derivs);
    gnorm = 0.;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 6; ++j) {
/* Computing 2nd power */
	    d__1 = derivs_ref(j, i__);
	    gnorm += d__1 * d__1;
	}
    }
    gnorm = sqrt(gnorm);
    grms = gnorm / sqrt((doublereal) group_1.ngrp);

/*     write out the final function and gradient values */

    if (inform_1.digits >= 8) {
	if (grms > 1e-8) {
	    io___24.ciunit = iounit_1.iout;
	    s_wsfe(&io___24);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___25.ciunit = iounit_1.iout;
	    s_wsfe(&io___25);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    } else if (inform_1.digits >= 6) {
	if (grms > 1e-6) {
	    io___26.ciunit = iounit_1.iout;
	    s_wsfe(&io___26);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___27.ciunit = iounit_1.iout;
	    s_wsfe(&io___27);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    } else {
	if (grms > 1e-4) {
	    io___28.ciunit = iounit_1.iout;
	    s_wsfe(&io___28);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___29.ciunit = iounit_1.iout;
	    s_wsfe(&io___29);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     write the final coordinates into a file */

    imin = freeunit_();
    o__1.oerr = 0;
    o__1.ounit = imin;
    o__1.ofnmlen = 120;
    o__1.ofnm = minfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    al__1.aerr = 0;
    al__1.aunit = imin;
    f_rew(&al__1);
    prtxyz_(&imin);
    cl__1.cerr = 0;
    cl__1.cunit = imin;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef keyline_ref
#undef derivs_ref
#undef rbc_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function optrigid1  --  energy and gradient for optrigid  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "optrigid1" is a service routine that computes the energy */
/*     and gradient for optimally conditioned variable metric */
/*     optimization of rigid bodies */


doublereal optrigid1_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int rigidxyz_(void);
    static doublereal e;
    static integer i__, j, nvar;
    static doublereal derivs[6000]	/* was [6][1000] */;
    extern /* Subroutine */ int gradrgd_(doublereal *, doublereal *);


#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]
#define derivs_ref(a_1,a_2) derivs[(a_2)*6 + a_1 - 7]



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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  rigid.i  --  rigid body coordinates for atom groups  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     xrb         rigid body reference x-coordinate for each atom */
/*     yrb         rigid body reference y-coordinate for each atom */
/*     zrb         rigid body reference z-coordinate for each atom */
/*     rbc         current rigid body coordinates for each group */
/*     use_rigid   flag to mark use of rigid body coordinate system */




/*     translate optimization parameters to rigid body coordinates */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    nvar = 0;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    ++nvar;
	    rbc_ref(j, i__) = xx[nvar];
	}
    }

/*     compute and store the energy and gradient */

    rigidxyz_();
    gradrgd_(&e, derivs);
    ret_val = e;

/*     store rigid body gradient as optimization gradient */

    nvar = 0;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    ++nvar;
	    g[nvar] = derivs_ref(j, i__);
	}
    }
    return ret_val;
} /* optrigid1_ */

#undef derivs_ref
#undef rbc_ref


/* Main program alias */ int optrigid_ () { MAIN__ (); return 0; }
