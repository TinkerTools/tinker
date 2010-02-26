/* optimize.f -- translated by f2c (version 20050501).
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
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

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
    doublereal scale[75000];
    logical set_scale__;
} scales_;

#define scales_1 scales_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  program optimize  --  variable metric Cartesian optimizer  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "optimize" performs energy minimization in Cartesian coordinate */
/*     space using an optimally conditioned variable metric method */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 Enter RMS Gradient per Atom Criterion"
	    "\002,\002 [0.01] :  \002,$)";
    static char fmt_40[] = "(f20.0)";
    static char fmt_50[] = "(/,\002 OPTIMIZE  --  Too many Parameters,\002"
	    ",\002 Increase the Value of MAXOPT\002)";
    static char fmt_60[] = "(/,\002 Final Function Value :\002,2x,f20.8,/"
	    ",\002 Final RMS Gradient :\002,4x,f20.8,/,\002 Final Gradient No"
	    "rm :\002,3x,f20.8)";
    static char fmt_70[] = "(/,\002 Final Function Value :\002,2x,f20.8,/"
	    ",\002 Final RMS Gradient :\002,4x,d20.8,/,\002 Final Gradient No"
	    "rm :\002,3x,d20.8)";
    static char fmt_80[] = "(/,\002 Final Function Value :\002,2x,f18.6,/"
	    ",\002 Final RMS Gradient :\002,4x,f18.6,/,\002 Final Gradient No"
	    "rm :\002,3x,f18.6)";
    static char fmt_90[] = "(/,\002 Final Function Value :\002,2x,f18.6,/"
	    ",\002 Final RMS Gradient :\002,4x,d18.6,/,\002 Final Gradient No"
	    "rm :\002,3x,d18.6)";
    static char fmt_100[] = "(/,\002 Final Function Value :\002,2x,f16.4,/"
	    ",\002 Final RMS Gradient :\002,4x,f16.4,/,\002 Final Gradient No"
	    "rm :\002,3x,f16.4)";
    static char fmt_110[] = "(/,\002 Final Function Value :\002,2x,f16.4,/"
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
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *, char 
	    *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);
    double sqrt(doublereal);
    integer f_rew(alist *);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    extern doublereal optimiz1_();
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static logical analytic;
    extern integer freeunit_(void);
    static integer i__, j;
    static doublereal xx[1000], eps;
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
    extern doublereal energy_(void);
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static doublereal derivs[75000]	/* was [3][25000] */;
    static char string[120];
    extern /* Subroutine */ int getxyz_(void), prtxyz_(integer *);
    static char minfile[120];
    extern /* Subroutine */ int initial_(void), numgrad_(D_fp, doublereal *, 
	    doublereal *), nextarg_(char *, logical *, ftnlen);
    static doublereal minimum;
    extern /* Subroutine */ int optsave_();
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___8 = { 1, string, 1, 0, 120, 1 };
    static icilist io___9 = { 1, string, 1, 0, 120, 1 };
    static icilist io___12 = { 1, string, 1, 0, 120, 1 };
    static cilist io___13 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_110, 0 };



#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1 - 4]
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
/*     ##  scales.i  --  parameter scale factors for optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     scale      multiplicative factor for each optimization parameter */
/*     set_scale  logical flag to show if scale factors have been set */




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




/*     set up the structure and mechanics calculation */

    initial_();
    getxyz_();
    mechanic_();

/*     use either analytical or numerical gradients */

    analytic = TRUE_;
    eps = 1e-5;

/*     search the keywords for output frequency parameters */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "PRINTOUT ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___8);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&inform_1.iprint, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "WRITEOUT ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___9);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&inform_1.iwrite, (ftnlen)
		    sizeof(integer));
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

/*     get termination criterion as RMS gradient per atom */

    grdmin = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___12);
	if (i__1 != 0) {
	    goto L20;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&grdmin, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L20;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L20;
	}
    }
L20:
    if (grdmin <= 0.) {
	io___13.ciunit = iounit_1.iout;
	s_wsfe(&io___13);
	e_wsfe();
	io___14.ciunit = iounit_1.input;
	s_rsfe(&io___14);
	do_fio(&c__1, (char *)&grdmin, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (grdmin <= 0.) {
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

/*     set scaling parameter for function and derivative values; */
/*     use square root of median eigenvalue of typical Hessian */

    scales_1.set_scale__ = TRUE_;
    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    for (j = 1; j <= 3; ++j) {
		++nvar;
		scales_1.scale[nvar - 1] = 12.;
	    }
	}
    }

/*     check for too many parameters to be optimized */

    if (nvar > 1000) {
	io___19.ciunit = iounit_1.iout;
	s_wsfe(&io___19);
	e_wsfe();
	fatal_();
    }

/*     scale the coordinates of each active atom */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    xx[nvar - 1] = atoms_1.x[i__ - 1] * scales_1.scale[nvar - 1];
	    ++nvar;
	    xx[nvar - 1] = atoms_1.y[i__ - 1] * scales_1.scale[nvar - 1];
	    ++nvar;
	    xx[nvar - 1] = atoms_1.z__[i__ - 1] * scales_1.scale[nvar - 1];
	}
    }

/*     make the call to the optimization routine */

    ocvm_(&nvar, xx, &minimum, &grdmin, (D_fp)optimiz1_, (U_fp)optsave_);

/*     unscale the final coordinates for active atoms */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    atoms_1.x[i__ - 1] = xx[nvar - 1] / scales_1.scale[nvar - 1];
	    ++nvar;
	    atoms_1.y[i__ - 1] = xx[nvar - 1] / scales_1.scale[nvar - 1];
	    ++nvar;
	    atoms_1.z__[i__ - 1] = xx[nvar - 1] / scales_1.scale[nvar - 1];
	}
    }

/*     compute the final function and RMS gradient values */

    if (analytic) {
	gradient_(&minimum, derivs);
    } else {
	minimum = energy_();
	numgrad_((D_fp)energy_, derivs, &eps);
    }
    gnorm = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    for (j = 1; j <= 3; ++j) {
/* Computing 2nd power */
		d__1 = derivs_ref(j, i__);
		gnorm += d__1 * d__1;
	    }
	}
    }
    gnorm = sqrt(gnorm);
    grms = gnorm / sqrt((doublereal) (nvar / 3));

/*     write out the final function and gradient values */

    if (inform_1.digits >= 8) {
	if (grms > 1e-8) {
	    io___25.ciunit = iounit_1.iout;
	    s_wsfe(&io___25);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___26.ciunit = iounit_1.iout;
	    s_wsfe(&io___26);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    } else if (inform_1.digits >= 6) {
	if (grms > 1e-6) {
	    io___27.ciunit = iounit_1.iout;
	    s_wsfe(&io___27);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___28.ciunit = iounit_1.iout;
	    s_wsfe(&io___28);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    } else {
	if (grms > 1e-4) {
	    io___29.ciunit = iounit_1.iout;
	    s_wsfe(&io___29);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___30.ciunit = iounit_1.iout;
	    s_wsfe(&io___30);
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




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  function optimiz1  --  energy and gradient for optimize  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "optimiz1" is a service routine that computes the energy and */
/*     gradient for optimally conditioned variable metric optimization */
/*     in Cartesian coordinate space */


doublereal optimiz1_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static logical analytic;
    static doublereal e;
    static integer i__;
    static doublereal eps;
    static integer nvar;
    extern doublereal energy_(void);
    static doublereal derivs[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int numgrad_(D_fp, doublereal *, doublereal *);


#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1 - 4]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  scales.i  --  parameter scale factors for optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     scale      multiplicative factor for each optimization parameter */
/*     set_scale  logical flag to show if scale factors have been set */




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




/*     use either analytical or numerical gradients */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    analytic = TRUE_;
    eps = 1e-5;

/*     translate optimization parameters to atomic coordinates */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    atoms_1.x[i__ - 1] = xx[nvar] / scales_1.scale[nvar - 1];
	    ++nvar;
	    atoms_1.y[i__ - 1] = xx[nvar] / scales_1.scale[nvar - 1];
	    ++nvar;
	    atoms_1.z__[i__ - 1] = xx[nvar] / scales_1.scale[nvar - 1];
	}
    }

/*     compute and store the energy and gradient */

    if (analytic) {
	gradient_(&e, derivs);
    } else {
	e = energy_();
	numgrad_((D_fp)energy_, derivs, &eps);
    }
    ret_val = e;

/*     store Cartesian gradient as optimization gradient */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    g[nvar] = derivs_ref(1, i__) / scales_1.scale[nvar - 1];
	    ++nvar;
	    g[nvar] = derivs_ref(2, i__) / scales_1.scale[nvar - 1];
	    ++nvar;
	    g[nvar] = derivs_ref(3, i__) / scales_1.scale[nvar - 1];
	}
    }
    return ret_val;
} /* optimiz1_ */

#undef derivs_ref


/* Main program alias */ int optimize_ () { MAIN__ (); return 0; }
