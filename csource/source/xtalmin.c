/* xtalmin.f -- translated by f2c (version 20050501).
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
    doublereal xbox, ybox, zbox, alpha, beta, gamma, xbox2, ybox2, zbox2, 
	    box34, lvec[9]	/* was [3][3] */, recip[9]	/* was [3][3] 
	    */, volbox, beta_sin__, beta_cos__, gamma_sin__, gamma_cos__, 
	    beta_term__, gamma_term__;
    logical orthogonal, monoclinic, triclinic, octahedron;
    char spacegrp[10];
} boxes_;

#define boxes_1 boxes_

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

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__2 = 2;



/*     ############################################################# */
/*     ##  COPYRIGHT (C) 2004 by Pengyu Ren & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  program xtalmin  --  full lattice crystal minimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "xtalmin" performs a full crystal energy minimization by */
/*     optimizing over fractional atomic coordinates and the six */
/*     lattice lengths and angles */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 Enter RMS Gradient per Atom Criterion"
	    "\002,\002 [0.01] :  \002,$)";
    static char fmt_40[] = "(f20.0)";
    static char fmt_50[] = "(/,\002 Initial Lattice Dimensions :    a   \002"
	    ",f12.4,/,\002                                 b   \002,f12.4,/"
	    ",\002                                 c   \002,f12.4,/,\002     "
	    "                           Alpha\002,f12.4,/,\002               "
	    "                 Beta \002,f12.4,/,\002                         "
	    "       Gamma\002,f12.4)";
    static char fmt_60[] = "(/,\002 Final Potential Function Value :\002,f16"
	    ".4,/,\002 Final RMS Coordinate Gradient : \002,f16.4,/,\002 Fina"
	    "l Coordinate Gradient Norm :\002,f16.4,/,\002 Final RMS Lattice "
	    "Gradient :    \002,f16.4,/,\002 Final Lattice Gradient Norm :  "
	    " \002,f16.4)";
    static char fmt_70[] = "(/,\002 Final Potential Function Value :\002,f16"
	    ".4,/,\002 Final RMS Coordinate Gradient : \002,d16.4,/,\002 Fina"
	    "l Coordinate Gradient Norm :\002,d16.4,/,\002 Final RMS Lattice "
	    "Gradient :    \002,d16.4,/,\002 Final Lattice Gradient Norm :  "
	    " \002,d16.4)";
    static char fmt_80[] = "(/,\002 Final Lattice Dimensions :      a   \002"
	    ",f12.4,/,\002                                 b   \002,f12.4,/"
	    ",\002                                 c   \002,f12.4,/,\002     "
	    "                           Alpha\002,f12.4,/,\002               "
	    "                 Beta \002,f12.4,/,\002                         "
	    "       Gamma\002,f12.4)";

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
    extern doublereal xtalmin1_(doublereal *, doublereal *);
    extern /* Subroutine */ int mechanic_(void), gradient_(doublereal *, 
	    doublereal *);
    extern integer freeunit_(void);
    static doublereal e;
    static integer i__, j;
    static doublereal xf[25000], yf[25000], zf[25000], xx[75000], glat[75000];
    static integer imin;
    extern /* Subroutine */ int ocvm_(integer *, doublereal *, doublereal *, 
	    doublereal *, D_fp, U_fp);
    static integer nvar;
    static doublereal grms;
    static integer next;
    extern /* Subroutine */ int final_(void), lbfgs_(integer *, doublereal *, 
	    doublereal *, doublereal *, D_fp, U_fp);
    static doublereal gnorm, glrms;
    static logical exist;
    static char record[120];
    static doublereal grdmin;
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static doublereal glnorm, derivs[75000]	/* was [3][25000] */;
    static char string[120];
    extern /* Subroutine */ int getxyz_(void), prtxyz_(integer *);
    static char minfile[120];
    extern /* Subroutine */ int lattice_(void), initial_(void), nextarg_(char 
	    *, logical *, ftnlen);
    static doublereal minimum;
    extern /* Subroutine */ int optsave_();
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___6 = { 1, string, 1, 0, 120, 1 };
    static icilist io___7 = { 1, string, 1, 0, 120, 1 };
    static icilist io___10 = { 1, string, 1, 0, 120, 1 };
    static cilist io___11 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_80, 0 };



#define lvec_ref(a_1,a_2) boxes_1.lvec[(a_2)*3 + a_1 - 4]
#define recip_ref(a_1,a_2) boxes_1.recip[(a_2)*3 + a_1 - 4]
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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  boxes.i  --  parameters for periodic boundary conditions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     xbox        length of a-axis of periodic box in Angstroms */
/*     ybox        length of b-axis of periodic box in Angstroms */
/*     zbox        length of c-axis of periodic box in Angstroms */
/*     alpha       angle between b- and c-axes of box in degrees */
/*     beta        angle between a- and c-axes of box in degrees */
/*     gamma       angle between a- and b-axes of box in degrees */
/*     xbox2       half of the a-axis length of periodic box */
/*     ybox2       half of the b-axis length of periodic box */
/*     zbox2       half of the c-axis length of periodic box */
/*     box34       three-fourths axis length of truncated octahedron */
/*     lvec        real space lattice vectors as matrix rows */
/*     recip       reciprocal lattice vectors as matrix columns */
/*     volbox      volume in Ang**3 of the periodic box */
/*     beta_sin    sine of the beta periodic box angle */
/*     beta_cos    cosine of the beta periodic box angle */
/*     gamma_sin   sine of the gamma periodic box angle */
/*     gamma_cos   cosine of the gamma periodic box angle */
/*     beta_term   term used in generating triclinic box */
/*     gamma_term  term used in generating triclinic box */
/*     orthogonal  flag to mark periodic box as orthogonal */
/*     monoclinic  flag to mark periodic box as monoclinic */
/*     triclinic   flag to mark periodic box as triclinic */
/*     octahedron  flag to mark box as truncated octahedron */
/*     spacegrp    space group symbol for the unitcell type */




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




/*     set up the structure and mechanics calculation */

    initial_();
    getxyz_();
    mechanic_();

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

/*     get termination criterion as RMS gradient per atom */

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

/*     write out the initial values of the lattice parameters */

    io___15.ciunit = iounit_1.iout;
    s_wsfe(&io___15);
    do_fio(&c__1, (char *)&boxes_1.xbox, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.ybox, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.zbox, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.alpha, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.beta, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.gamma, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     set scale factors to apply to optimization variables */

    scales_1.set_scale__ = TRUE_;
    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	scales_1.scale[nvar - 1] = boxes_1.xbox * 12.;
	++nvar;
	scales_1.scale[nvar - 1] = boxes_1.ybox * 12.;
	++nvar;
	scales_1.scale[nvar - 1] = boxes_1.zbox * 12.;
    }
    scales_1.scale[nvar] = sqrt(boxes_1.xbox) * 4.;
    scales_1.scale[nvar + 1] = sqrt(boxes_1.ybox) * 4.;
    scales_1.scale[nvar + 2] = sqrt(boxes_1.zbox) * 4.;
    scales_1.scale[nvar + 3] = sqrt(boxes_1.volbox) * .02;
    scales_1.scale[nvar + 4] = sqrt(boxes_1.volbox) * .02;
    scales_1.scale[nvar + 5] = sqrt(boxes_1.volbox) * .02;

/*     compute the fractional coordinates for each atom */

    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = i__ * 3 - 3;
	xx[j] = atoms_1.x[i__ - 1] * recip_ref(1, 1) + atoms_1.y[i__ - 1] * 
		recip_ref(2, 1) + atoms_1.z__[i__ - 1] * recip_ref(3, 1);
	xx[j + 1] = atoms_1.x[i__ - 1] * recip_ref(1, 2) + atoms_1.y[i__ - 1] 
		* recip_ref(2, 2) + atoms_1.z__[i__ - 1] * recip_ref(3, 2);
	xx[j + 2] = atoms_1.x[i__ - 1] * recip_ref(1, 3) + atoms_1.y[i__ - 1] 
		* recip_ref(2, 3) + atoms_1.z__[i__ - 1] * recip_ref(3, 3);
    }

/*     scale the fractional coordinates and lattice parameters */

    nvar = atoms_1.n * 3;
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xx[i__ - 1] *= scales_1.scale[i__ - 1];
    }
    xx[nvar] = boxes_1.xbox * scales_1.scale[nvar];
    xx[nvar + 1] = boxes_1.ybox * scales_1.scale[nvar + 1];
    xx[nvar + 2] = boxes_1.zbox * scales_1.scale[nvar + 2];
    xx[nvar + 3] = boxes_1.alpha * scales_1.scale[nvar + 3];
    xx[nvar + 4] = boxes_1.beta * scales_1.scale[nvar + 4];
    xx[nvar + 5] = boxes_1.gamma * scales_1.scale[nvar + 5];
    nvar += 6;

/*     make the call to the optimization routine */

    if (nvar <= 1000) {
	ocvm_(&nvar, xx, &minimum, &grdmin, (D_fp)xtalmin1_, (U_fp)optsave_);
    } else {
	lbfgs_(&nvar, xx, &minimum, &grdmin, (D_fp)xtalmin1_, (U_fp)optsave_);
    }

/*     unscale fractional coordinates and get atomic coordinates */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = i__ * 3 - 3;
	xf[i__ - 1] = xx[j] / scales_1.scale[j];
	yf[i__ - 1] = xx[j + 1] / scales_1.scale[j + 1];
	zf[i__ - 1] = xx[j + 2] / scales_1.scale[j + 2];
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }

/*     compute final energy value and coordinate RMS gradient */

    gradient_(&e, derivs);
    gnorm = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* Computing 2nd power */
	    d__1 = derivs_ref(j, i__);
	    gnorm += d__1 * d__1;
	}
    }
    gnorm = sqrt(gnorm);
    nvar = atoms_1.n * 3;
    grms = gnorm / sqrt((doublereal) (nvar / 3));

/*     compute the final RMS gradient for lattice parameters */

    minimum = xtalmin1_(xx, glat);
    glnorm = 0.;
    i__1 = nvar + 6;
    for (i__ = nvar + 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = scales_1.scale[i__ - 1] * glat[i__ - 1];
	glnorm += d__1 * d__1;
    }
    glnorm = sqrt(glnorm);
    glrms = glnorm / sqrt(6.);

/*     write out the final energy and coordinate gradients */

    if (grms > 1e-4 && glrms > 1e-4) {
	io___30.ciunit = iounit_1.iout;
	s_wsfe(&io___30);
	do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&glrms, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&glnorm, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else {
	io___31.ciunit = iounit_1.iout;
	s_wsfe(&io___31);
	do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&gnorm, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&glrms, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&glnorm, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     write out the final values of the lattice parameters */

    io___32.ciunit = iounit_1.iout;
    s_wsfe(&io___32);
    do_fio(&c__1, (char *)&boxes_1.xbox, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.ybox, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.zbox, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.alpha, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.beta, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.gamma, (ftnlen)sizeof(doublereal));
    e_wsfe();

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
#undef recip_ref
#undef lvec_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function xtalmin1  --  energy and gradient for lattice  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "xtalmin1" is a service routine that computes the energy and */
/*     gradient with respect to fractional coordinates and lattice */
/*     dimensions for a crystal energy minimization */


doublereal xtalmin1_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static doublereal e;
    static integer i__, j;
    static doublereal e0, xf[25000], yf[25000], zf[25000], old, eps;
    extern doublereal energy_(void);
    static doublereal derivs[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int lattice_(void);


#define lvec_ref(a_1,a_2) boxes_1.lvec[(a_2)*3 + a_1 - 4]
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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  boxes.i  --  parameters for periodic boundary conditions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     xbox        length of a-axis of periodic box in Angstroms */
/*     ybox        length of b-axis of periodic box in Angstroms */
/*     zbox        length of c-axis of periodic box in Angstroms */
/*     alpha       angle between b- and c-axes of box in degrees */
/*     beta        angle between a- and c-axes of box in degrees */
/*     gamma       angle between a- and b-axes of box in degrees */
/*     xbox2       half of the a-axis length of periodic box */
/*     ybox2       half of the b-axis length of periodic box */
/*     zbox2       half of the c-axis length of periodic box */
/*     box34       three-fourths axis length of truncated octahedron */
/*     lvec        real space lattice vectors as matrix rows */
/*     recip       reciprocal lattice vectors as matrix columns */
/*     volbox      volume in Ang**3 of the periodic box */
/*     beta_sin    sine of the beta periodic box angle */
/*     beta_cos    cosine of the beta periodic box angle */
/*     gamma_sin   sine of the gamma periodic box angle */
/*     gamma_cos   cosine of the gamma periodic box angle */
/*     beta_term   term used in generating triclinic box */
/*     gamma_term  term used in generating triclinic box */
/*     orthogonal  flag to mark periodic box as orthogonal */
/*     monoclinic  flag to mark periodic box as monoclinic */
/*     triclinic   flag to mark periodic box as triclinic */
/*     octahedron  flag to mark box as truncated octahedron */
/*     spacegrp    space group symbol for the unitcell type */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  scales.i  --  parameter scale factors for optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     scale      multiplicative factor for each optimization parameter */
/*     set_scale  logical flag to show if scale factors have been set */




/*     translate optimization variables to fractional coordinates */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = i__ * 3 - 3;
	xf[i__ - 1] = xx[j + 1] / scales_1.scale[j];
	yf[i__ - 1] = xx[j + 2] / scales_1.scale[j + 1];
	zf[i__ - 1] = xx[j + 3] / scales_1.scale[j + 2];
    }

/*     translate optimization variables to lattice parameters */

    boxes_1.xbox = xx[atoms_1.n * 3 + 1] / scales_1.scale[atoms_1.n * 3];
    boxes_1.ybox = xx[atoms_1.n * 3 + 2] / scales_1.scale[atoms_1.n * 3 + 1];
    boxes_1.zbox = xx[atoms_1.n * 3 + 3] / scales_1.scale[atoms_1.n * 3 + 2];
    boxes_1.alpha = xx[atoms_1.n * 3 + 4] / scales_1.scale[atoms_1.n * 3 + 3];
    boxes_1.beta = xx[atoms_1.n * 3 + 5] / scales_1.scale[atoms_1.n * 3 + 4];
    boxes_1.gamma = xx[atoms_1.n * 3 + 6] / scales_1.scale[atoms_1.n * 3 + 5];

/*     update current atomic coordinates based on optimization values */

    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }

/*     find energy and fractional coordinates deriviatives */

    gradient_(&e, derivs);
    ret_val = e;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = i__ * 3 - 3;
	g[j + 1] = derivs_ref(1, i__) * lvec_ref(1, 1) + derivs_ref(2, i__) * 
		lvec_ref(1, 2) + derivs_ref(3, i__) * lvec_ref(1, 3);
	g[j + 2] = derivs_ref(1, i__) * lvec_ref(2, 1) + derivs_ref(2, i__) * 
		lvec_ref(2, 2) + derivs_ref(3, i__) * lvec_ref(2, 3);
	g[j + 3] = derivs_ref(1, i__) * lvec_ref(3, 1) + derivs_ref(2, i__) * 
		lvec_ref(3, 2) + derivs_ref(3, i__) * lvec_ref(3, 3);
    }

/*     find derivative with respect to lattice a-axis length */

    eps = 1e-4;
    old = boxes_1.xbox;
    boxes_1.xbox -= eps * .5;
    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }
    e0 = energy_();
    boxes_1.xbox += eps;
    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }
    e = energy_();
    g[atoms_1.n * 3 + 1] = (e - e0) / eps;
    boxes_1.xbox = old;

/*     find derivative with respect to lattice b-axis length */

    old = boxes_1.ybox;
    boxes_1.ybox -= eps * .5;
    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }
    e0 = energy_();
    boxes_1.ybox += eps;
    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }
    e = energy_();
    g[atoms_1.n * 3 + 2] = (e - e0) / eps;
    boxes_1.ybox = old;

/*     find derivative with respect to lattice c-axis length */

    old = boxes_1.zbox;
    boxes_1.zbox -= eps * .5;
    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }
    e0 = energy_();
    boxes_1.zbox += eps;
    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }
    e = energy_();
    g[atoms_1.n * 3 + 3] = (e - e0) / eps;
    boxes_1.zbox = old;

/*     find derivative with respect to lattice alpha angle */

    eps *= 57.29577951308232088;
    old = boxes_1.alpha;
    boxes_1.alpha -= eps * .5;
    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }
    e0 = energy_();
    boxes_1.alpha += eps;
    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }
    e = energy_();
    g[atoms_1.n * 3 + 4] = (e - e0) / eps;
    boxes_1.alpha = old;

/*     find derivative with respect to lattice beta angle */

    old = boxes_1.beta;
    boxes_1.beta -= eps * .5;
    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }
    e0 = energy_();
    boxes_1.beta += eps;
    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }
    e = energy_();
    g[atoms_1.n * 3 + 5] = (e - e0) / eps;
    boxes_1.beta = old;

/*     find derivative with respect to lattice gamma angle */

    old = boxes_1.gamma;
    boxes_1.gamma -= eps * .5;
    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }
    e0 = energy_();
    boxes_1.gamma += eps;
    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }
    e = energy_();
    g[atoms_1.n * 3 + 6] = (e - e0) / eps;
    boxes_1.gamma = old;

/*     revert to the original atomic coordinate values */

    lattice_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 1) + yf[i__ - 1] * 
		lvec_ref(2, 1) + zf[i__ - 1] * lvec_ref(3, 1);
	atoms_1.y[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 2) + yf[i__ - 1] * 
		lvec_ref(2, 2) + zf[i__ - 1] * lvec_ref(3, 2);
	atoms_1.z__[i__ - 1] = xf[i__ - 1] * lvec_ref(1, 3) + yf[i__ - 1] * 
		lvec_ref(2, 3) + zf[i__ - 1] * lvec_ref(3, 3);
    }

/*     apply scale factors to the coordinate and lattice gradient */

    i__1 = atoms_1.n * 3 + 6;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] /= scales_1.scale[i__ - 1];
    }
    return ret_val;
} /* xtalmin1_ */

#undef derivs_ref
#undef lvec_ref


/* Main program alias */ int xtalmin_ () { MAIN__ (); return 0; }
