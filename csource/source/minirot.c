/* minirot.f -- translated by f2c (version 20050501).
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
    doublereal dihed[1000];
    integer nomega, iomega[2000]	/* was [2][1000] */, zline[1000];
} omega_;

#define omega_1 omega_

struct {
    doublereal scale[75000];
    logical set_scale__;
} scales_;

#define scales_1 scales_

struct {
    doublereal zbond[25000], zang[25000], ztors[25000];
    integer iz[100000]	/* was [4][25000] */;
} zcoord_;

#define zcoord_1 zcoord_

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
/*     ##  program minirot  --  low storage BFGS torsional optimizer  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "minirot" performs an energy minimization in torsional */
/*     angle space using a low storage BFGS nonlinear optimization */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 Enter RMS Gradient per Torsion Criterio"
	    "n\002,\002 [0.01] :  \002,$)";
    static char fmt_40[] = "(f20.0)";
    static char fmt_50[] = "(/,\002 Final Function Value :\002,2x,f20.8,/"
	    ",\002 Final RMS Gradient :\002,4x,f20.8,/,\002 Final Gradient No"
	    "rm :\002,3x,f20.8)";
    static char fmt_60[] = "(/,\002 Final Function Value :\002,2x,f20.8,/"
	    ",\002 Final RMS Gradient :\002,4x,d20.8,/,\002 Final Gradient No"
	    "rm :\002,3x,d20.8)";
    static char fmt_70[] = "(/,\002 Final Function Value :\002,2x,f18.6,/"
	    ",\002 Final RMS Gradient :\002,4x,f18.6,/,\002 Final Gradient No"
	    "rm :\002,3x,f18.6)";
    static char fmt_80[] = "(/,\002 Final Function Value :\002,2x,f18.6,/"
	    ",\002 Final RMS Gradient :\002,4x,d18.6,/,\002 Final Gradient No"
	    "rm :\002,3x,d18.6)";
    static char fmt_90[] = "(/,\002 Final Function Value :\002,2x,f16.4,/"
	    ",\002 Final RMS Gradient :\002,4x,f16.4,/,\002 Final Gradient No"
	    "rm :\002,3x,f16.4)";
    static char fmt_100[] = "(/,\002 Final Function Value :\002,2x,f16.4,/"
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
    extern doublereal minirot1_();
    extern /* Subroutine */ int mechanic_(void);
    extern integer freeunit_(void);
    static integer i__;
    static doublereal xx[75000];
    static integer imin;
    static doublereal grms;
    static integer next;
    extern /* Subroutine */ int final_(void), lbfgs_(integer *, doublereal *, 
	    doublereal *, doublereal *, D_fp, U_fp);
    static doublereal gnorm;
    static logical exist;
    static char record[120];
    static doublereal grdmin;
    extern /* Subroutine */ int getint_(void);
    static doublereal derivs[1000];
    static char string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), prtint_(integer *);
    static char minfile[120];
    extern /* Subroutine */ int initial_(void), gradrot_(doublereal *, 
	    doublereal *), nextarg_(char *, logical *, ftnlen);
    static doublereal minimum;
    extern /* Subroutine */ int optsave_();
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen), initrot_(void);

    /* Fortran I/O blocks */
    static icilist io___6 = { 1, string, 1, 0, 120, 1 };
    static icilist io___7 = { 1, string, 1, 0, 120, 1 };
    static icilist io___10 = { 1, string, 1, 0, 120, 1 };
    static cilist io___11 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_100, 0 };



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     zbond   bond length used to define each Z-matrix atom */
/*     zang    bond angle used to define each Z-matrix atom */
/*     ztors   angle or torsion used to define Z-matrix atom */
/*     iz      defining atom numbers for each Z-matrix atom */




/*     set up the molecular mechanics calculation */

    initial_();
    getint_();
    mechanic_();
    initrot_();

/*     search the keywords for output frequency parameters */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "PRINTOUT ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___6);
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
	    i__2 = s_rsli(&io___7);
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

/*     get termination criterion as RMS torsional gradient */

    grdmin = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___10);
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
	io___11.ciunit = iounit_1.iout;
	s_wsfe(&io___11);
	e_wsfe();
	io___12.ciunit = iounit_1.input;
	s_rsfe(&io___12);
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
    i__3[1] = 4, a__1[1] = ".int";
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
    prtint_(&imin);
    cl__1.cerr = 0;
    cl__1.cunit = imin;
    cl__1.csta = 0;
    f_clos(&cl__1);
    s_copy(files_1.outfile, minfile, (ftnlen)120, (ftnlen)120);

/*     set scaling parameter for function and derivative values; */
/*     use square root of median eigenvalue of typical Hessian */

    scales_1.set_scale__ = TRUE_;
    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	scales_1.scale[i__ - 1] = 5.;
    }

/*     get optimization parameters as scaled dihedral angles */

    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xx[i__ - 1] = omega_1.dihed[i__ - 1] * scales_1.scale[i__ - 1];
    }

/*     make the call to the optimization routine */

    lbfgs_(&omega_1.nomega, xx, &minimum, &grdmin, (D_fp)minirot1_, (U_fp)
	    optsave_);

/*     unscale the final dihedral angle values */

    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	omega_1.dihed[i__ - 1] = xx[i__ - 1] / scales_1.scale[i__ - 1];
	zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] = omega_1.dihed[i__ - 1] * 
		57.29577951308232088;
    }

/*     compute the final function and RMS gradient values */

    gradrot_(&minimum, derivs);
    gnorm = 0.;
    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = derivs[i__ - 1];
	gnorm += d__1 * d__1;
    }
    gnorm = sqrt(gnorm);
    grms = gnorm / sqrt((doublereal) omega_1.nomega);

/*     write out the final function and gradient values */

    if (inform_1.digits >= 8) {
	if (grms > 1e-8) {
	    io___20.ciunit = iounit_1.iout;
	    s_wsfe(&io___20);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___21.ciunit = iounit_1.iout;
	    s_wsfe(&io___21);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    } else if (inform_1.digits >= 6) {
	if (grms > 1e-6) {
	    io___22.ciunit = iounit_1.iout;
	    s_wsfe(&io___22);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___23.ciunit = iounit_1.iout;
	    s_wsfe(&io___23);
	    do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    } else {
	if (grms > 1e-4) {
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
    prtint_(&imin);
    cl__1.cerr = 0;
    cl__1.cunit = imin;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef keyline_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function minirot1  --  energy and gradient for minirot  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "minirot1" is a service routine that computes the energy */
/*     and gradient for a low storage BFGS nonlinear optimization */
/*     in torsional angle space */


doublereal minirot1_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static doublereal e;
    static integer i__;
    static doublereal derivs[1000];
    extern /* Subroutine */ int gradrot_(doublereal *, doublereal *), 
	    makexyz_(void);



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     zbond   bond length used to define each Z-matrix atom */
/*     zang    bond angle used to define each Z-matrix atom */
/*     ztors   angle or torsion used to define Z-matrix atom */
/*     iz      defining atom numbers for each Z-matrix atom */




/*     translate optimization variables into dihedrals */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	omega_1.dihed[i__ - 1] = xx[i__] / scales_1.scale[i__ - 1];
	zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] = omega_1.dihed[i__ - 1] * 
		57.29577951308232088;
    }

/*     get coordinates, then compute energy and gradient */

    makexyz_();
    gradrot_(&e, derivs);
    ret_val = e;

/*     store torsional gradient as optimization gradient */

    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = derivs[i__ - 1] / scales_1.scale[i__ - 1];
    }
    return ret_val;
} /* minirot1_ */

/* Main program alias */ int minirot_ () { MAIN__ (); return 0; }
