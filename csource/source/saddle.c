/* saddle.f -- translated by f2c (version 20050501).
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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

struct {
    doublereal stpmin, stpmax, cappa, slpmax, angmax;
    integer intmax;
} linmin_;

#define linmin_1 linmin_

struct {
    doublereal t, pm, xmin1[75000], xmin2[75000], xm[75000];
} syntrn_;

#define syntrn_1 syntrn_

struct {
    integer ltitle;
    char title[120];
} titles_;

#define titles_1 titles_

struct {
    doublereal zbond[25000], zang[25000], ztors[25000];
    integer iz[100000]	/* was [4][25000] */;
} zcoord_;

#define zcoord_1 zcoord_

struct {
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    doublereal fctmin, hguess;
    integer maxiter, nextiter;
} minima_;

#define minima_1 minima_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  program saddle  --  find conformational transition state  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "saddle" finds a transition state between two conformational */
/*     minima using a combination of ideas from the synchronous transit */
/*     (Halgren-Lipscomb) and quadratic path (Bell-Crighton) methods */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 Enter RMS Gradient per Atom Criterion [0"
	    ".1] :  \002,$)";
    static char fmt_40[] = "(f20.0)";
    static char fmt_50[] = "(/,\002 Perform Synchronous Transit Pathway Sc"
	    "ans\002,\002 [N] :  \002,$)";
    static char fmt_60[] = "(a120)";
    static char fmt_70[] = "(/,\002 RMS Fit for All Atoms of both Structures"
	    " :\002,f10.4)";
    static char fmt_80[] = "(/,\002 Energy Value for Endpoint Structure 1 "
	    ":\002,f13.4,/,\002 Energy Value for Endpoint Structure 2 :\002,f"
	    "13.4)";
    static char fmt_90[] = "(/,\002 Using TSTATE.XYZ as the Transition State"
	    " Estimate\002)";
    static char fmt_110[] = "(/,\002 Search for a Maximum along Synchronous "
	    "Transit :\002,/\002 ST Iter    F Value       Path      RMS G\002,"
	    "\002      G Tan      Gamma   FG Call\002,/)";
    static char fmt_120[] = "(i6,f13.4,f11.4,f11.4,f11.4,f11.5,i8)";
    static char fmt_140[] = "(i6,f13.4,f11.4,f11.4,f11.4,f11.5,i8)";
    static char fmt_160[] = "(/,\002 SADDLE  --  Termination due to Los"
	    "s\002,\002 of Negative Curvature\002)";
    static char fmt_170[] = "(/,\002 Search for a Minimum in Conjugate Direc"
	    "tions :\002,/,\002 CG Iter    F Value      RMS G     F Move\002"
	    ",\002    X Move    Angle   FG Call  Comment\002,/)";
    static char fmt_180[] = "(i6,f13.4,f11.4,30x,i7)";
    static char fmt_190[] = "(/,\002 SADDLE  --  Normal Termination at\002"
	    ",\002 Transition State\002)";
    static char fmt_200[] = "(i6,f13.4,f11.4,f11.4,f10.4,f9.2,i7,3x,a9)";
    static char fmt_210[] = "(/,\002 SADDLE  --  Normal Termination at\002"
	    ",\002 Transition State\002)";
    static char fmt_220[] = "(/,\002 SADDLE  --  Termination due to Maximu"
	    "m\002,\002 Iteration Limit\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *, char 
	    *, ftnlen), e_rsfe(void), f_inqu(inlist *), f_open(olist *), 
	    f_rew(alist *), f_clos(cllist *);
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    static doublereal gammamin;
    extern /* Subroutine */ int pathscan_(integer *, doublereal *, doublereal 
	    *, integer *);
    static integer maxcycle;
    static logical newcycle;
    extern integer freeunit_(void);
    static integer maxinner;
    static doublereal rmsvalue;
    static integer maxouter;
    static doublereal f, g[75000];
    static integer i__;
    static doublereal p, s[75000], h0[75000], g2, s0[75000];
    static logical terminate;
    static doublereal x1[25000], x2[25000], y1[25000], y2[25000], z1[25000], 
	    z2[25000], hg, sg, xx[75000], f_0__, f_1__, f_2__, f_3__, sg0, 
	    tan__[75000];
    static integer its;
    static doublereal beta, dgdt[75000];
    static logical scan, done;
    static integer nvar, next;
    static doublereal zang1[25000], zang2[25000], gamma, f_old__, g_old__[
	    75000], angle;
    extern /* Subroutine */ int fatal_(void);
    static doublereal g_tan__, delta;
    extern /* Subroutine */ int final_(void);
    static doublereal x_old__[75000], g_rms__;
    static integer niter;
    static logical exist;
    static doublereal g2_old__, zbond1[25000], zbond2[25000], ztors1[25000], 
	    ztors2[25000];
    extern /* Subroutine */ int search_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , D_fp, char *, ftnlen);
    static doublereal reduce;
    static integer ncalls, ncycle;
    static doublereal grdmin, f_move__;
    static integer ninner;
    static char tsfile[120];
    static doublereal x_move__;
    static char answer[1], record[120];
    static integer nouter;
    extern doublereal saddle1_(doublereal *, doublereal *);
    static char string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), impose_(integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *);
    static char status[9];
    extern /* Subroutine */ int getxyz_(void);
    static doublereal energy1, energy2;
    extern /* Subroutine */ int prtxyz_(integer *);
    static doublereal diverge;
    static logical spanned;
    extern /* Subroutine */ int initial_(void), makeint_(integer *), pathval_(
	    integer *, doublereal *), tangent_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal epsilon;
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), pathpnt_(
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    ;
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen), readxyz_(
	    integer *), makexyz_(void);

    /* Fortran I/O blocks */
    static icilist io___22 = { 1, string, 1, 0, 120, 1 };
    static icilist io___23 = { 1, string, 1, 0, 120, 1 };
    static icilist io___24 = { 1, string, 1, 0, 120, 1 };
    static icilist io___27 = { 1, string, 1, 0, 120, 1 };
    static cilist io___28 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___46 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___69 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___79 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___80 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___81 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___89 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___90 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___91 = { 0, 0, 0, fmt_220, 0 };



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  linmin.i  --  parameters for line search minimization  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     stpmin   minimum step length in current line search direction */
/*     stpmax   maximum step length in current line search direction */
/*     cappa    stringency of line search (0=tight < cappa < 1=loose) */
/*     slpmax   projected gradient above which stepsize is reduced */
/*     angmax   maximum angle between search direction and -gradient */
/*     intmax   maximum number of interpolations during line search */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  syntrn.i  --  definition of synchronous transit path  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     t       value of the path coordinate (0=reactant, 1=product) */
/*     pm      path coordinate for extra point in quadratic transit */
/*     xmin1   reactant coordinates as array of optimization variables */
/*     xmin2   product coordinates as array of optimization variables */
/*     xm      extra coordinate set for quadratic synchronous transit */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  titles.i  --  title for the current molecular system  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     ltitle   length in characters of the nonblank title string */
/*     title    title used to describe the current structure */




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




/*     set default parameters for the saddle point run */

    initial_();
    terminate = FALSE_;
    ncalls = 0;
    nouter = 0;
    maxouter = 100;
    maxinner = 50;
    maxcycle = 4;
    epsilon = .5;
    gammamin = 1e-5;
    diverge = .005;
    reduce = 0.;

/*     get coordinates for the first endpoint structure */

    getxyz_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x1[i__ - 1] = atoms_1.x[i__ - 1];
	y1[i__ - 1] = atoms_1.y[i__ - 1];
	z1[i__ - 1] = atoms_1.z__[i__ - 1];
    }

/*     get coordinates for the second endpoint structure */

    getxyz_();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x2[i__ - 1] = atoms_1.x[i__ - 1];
	y2[i__ - 1] = atoms_1.y[i__ - 1];
	z2[i__ - 1] = atoms_1.z__[i__ - 1];
    }

/*     setup for the subsequent energy computations */

    mechanic_();

/*     get any altered values from the keyword file */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "DIVERGE ", (ftnlen)8, (ftnlen)8) == 0) {
	    i__2 = s_rsli(&io___22);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&diverge, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "REDUCE ", (ftnlen)7, (ftnlen)7) == 0) {
	    i__2 = s_rsli(&io___23);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&reduce, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "GAMMAMIN ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___24);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&gammamin, (ftnlen)sizeof(
		    doublereal));
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
	i__1 = s_rsli(&io___27);
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
	io___28.ciunit = iounit_1.iout;
	s_wsfe(&io___28);
	e_wsfe();
	io___29.ciunit = iounit_1.input;
	s_rsfe(&io___29);
	do_fio(&c__1, (char *)&grdmin, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (grdmin <= 0.) {
	grdmin = .1;
    }

/*     find out whether syncronous transit scans are desired */

    scan = FALSE_;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___32.ciunit = iounit_1.iout;
	s_wsfe(&io___32);
	e_wsfe();
	io___33.ciunit = iounit_1.input;
	s_rsfe(&io___33);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'Y') {
	scan = TRUE_;
    }

/*     superimpose the two conformational endpoints */

    impose_(&atoms_1.n, x1, y1, z1, &atoms_1.n, x2, y2, z2, &rmsvalue);
    io___35.ciunit = iounit_1.iout;
    s_wsfe(&io___35);
    do_fio(&c__1, (char *)&rmsvalue, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     copy the superimposed structures into vectors */

    nvar = atoms_1.n * 3;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	syntrn_1.xmin1[i__ * 3 - 3] = x1[i__ - 1];
	syntrn_1.xmin1[i__ * 3 - 2] = y1[i__ - 1];
	syntrn_1.xmin1[i__ * 3 - 1] = z1[i__ - 1];
	syntrn_1.xmin2[i__ * 3 - 3] = x2[i__ - 1];
	syntrn_1.xmin2[i__ * 3 - 2] = y2[i__ - 1];
	syntrn_1.xmin2[i__ * 3 - 1] = z2[i__ - 1];
    }

/*     get and store internal coordinates for first endpoint */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = x1[i__ - 1];
	atoms_1.y[i__ - 1] = y1[i__ - 1];
	atoms_1.z__[i__ - 1] = z1[i__ - 1];
    }
    makeint_(&c__0);
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zbond1[i__ - 1] = zcoord_1.zbond[i__ - 1];
	zang1[i__ - 1] = zcoord_1.zang[i__ - 1];
	ztors1[i__ - 1] = zcoord_1.ztors[i__ - 1];
    }

/*     get and store internal coordinates for second endpoint */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = x2[i__ - 1];
	atoms_1.y[i__ - 1] = y2[i__ - 1];
	atoms_1.z__[i__ - 1] = z2[i__ - 1];
    }
    makeint_(&c__2);
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zbond2[i__ - 1] = zcoord_1.zbond[i__ - 1];
	zang2[i__ - 1] = zcoord_1.zang[i__ - 1];
	ztors2[i__ - 1] = zcoord_1.ztors[i__ - 1];
	if (ztors1[i__ - 1] - ztors2[i__ - 1] > 180.) {
	    ztors2[i__ - 1] += 360.;
	} else if (ztors1[i__ - 1] - ztors2[i__ - 1] < -180.) {
	    ztors1[i__ - 1] += 360.;
	}
    }

/*     get the energies for the two endpoint structures */

    ncalls += 2;
    energy1 = saddle1_(syntrn_1.xmin1, g);
    energy2 = saddle1_(syntrn_1.xmin2, g);
    io___46.ciunit = iounit_1.iout;
    s_wsfe(&io___46);
    do_fio(&c__1, (char *)&energy1, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&energy2, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     make a guess at the transition state structure; */
/*     or use the current guess if one is around */

    ioin__1.inerr = 0;
    ioin__1.infilen = 10;
    ioin__1.infile = "tstate.xyz";
    ioin__1.inex = &exist;
    ioin__1.inopen = 0;
    ioin__1.innum = 0;
    ioin__1.innamed = 0;
    ioin__1.inname = 0;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    if (exist) {
	io___47.ciunit = iounit_1.iout;
	s_wsfe(&io___47);
	e_wsfe();
	its = freeunit_();
	s_copy(tsfile, "tstate.xyz", (ftnlen)120, (ftnlen)10);
	version_(tsfile, "old", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = its;
	o__1.ofnmlen = 120;
	o__1.ofnm = tsfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = its;
	f_rew(&al__1);
	readxyz_(&its);
	cl__1.cerr = 0;
	cl__1.cunit = its;
	cl__1.csta = 0;
	f_clos(&cl__1);
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xx[i__ * 3 - 3] = atoms_1.x[i__ - 1];
	    xx[i__ * 3 - 2] = atoms_1.y[i__ - 1];
	    xx[i__ * 3 - 1] = atoms_1.z__[i__ - 1];
	}
    } else {
	syntrn_1.t = .5;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    zcoord_1.zbond[i__ - 1] = (1. - syntrn_1.t) * zbond1[i__ - 1] + 
		    syntrn_1.t * zbond2[i__ - 1];
	    zcoord_1.zang[i__ - 1] = (1. - syntrn_1.t) * zang1[i__ - 1] + 
		    syntrn_1.t * zang2[i__ - 1];
	    zcoord_1.ztors[i__ - 1] = (1. - syntrn_1.t) * ztors1[i__ - 1] + 
		    syntrn_1.t * ztors2[i__ - 1];
	}
	makexyz_();
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xx[i__ * 3 - 3] = atoms_1.x[i__ - 1];
	    xx[i__ * 3 - 2] = atoms_1.y[i__ - 1];
	    xx[i__ * 3 - 1] = atoms_1.z__[i__ - 1];
	}
    }

/*     save the initial estimate of the transition state */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xx[i__ * 3 - 3];
	atoms_1.y[i__ - 1] = xx[i__ * 3 - 2];
	atoms_1.z__[i__ - 1] = xx[i__ * 3 - 1];
    }
    if (! exist) {
	s_copy(titles_1.title, "Transition State Structure", (ftnlen)120, (
		ftnlen)26);
	titles_1.ltitle = 26;
    }
    its = freeunit_();
    s_copy(tsfile, "tstate.xyz", (ftnlen)120, (ftnlen)10);
    version_(tsfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = its;
    o__1.ofnmlen = 120;
    o__1.ofnm = tsfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    prtxyz_(&its);
    cl__1.cerr = 0;
    cl__1.cunit = its;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     start of the major loop for transition state location; */
/*     first, find the value of the transit path coordinate */

L100:
    ++nouter;
    pathval_(&nvar, xx);

/*     make a scan along the synchronous transit pathway */

    if (scan) {
	pathscan_(&nvar, syntrn_1.xmin1, syntrn_1.xmin2, &ncalls);
    }

/*     set parameters for use in quadratic line maximization */

    done = FALSE_;
    niter = 1;
    ncycle = 1;
    delta = .01;

/*     compute initial point for quadratic line maximization */

    syntrn_1.t = syntrn_1.pm;
    pathpnt_(&nvar, &syntrn_1.t, xx, syntrn_1.xmin1, syntrn_1.xmin2);
    ncalls += 3;
    f = saddle1_(xx, g);
    tangent_(&nvar, xx, g, &g_rms__, tan__, &g_tan__, &gamma, dgdt);
    io___61.ciunit = iounit_1.iout;
    s_wsfe(&io___61);
    e_wsfe();
    io___62.ciunit = iounit_1.iout;
    s_wsfe(&io___62);
    do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&syntrn_1.t, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&g_rms__, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&g_tan__, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&gamma, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
    e_wsfe();

/*     make an iterative search for quadratic line maximum */

    while(! done) {
	f_0__ = f;
	syntrn_1.t += delta;
	pathpnt_(&nvar, &syntrn_1.t, xx, syntrn_1.xmin1, syntrn_1.xmin2);
	++ncalls;
	f_1__ = saddle1_(xx, g);
	syntrn_1.t -= delta * 2.;
	pathpnt_(&nvar, &syntrn_1.t, xx, syntrn_1.xmin1, syntrn_1.xmin2);
	syntrn_1.t += delta;
	++ncalls;
	f_2__ = saddle1_(xx, g);
	if (f_1__ > f_0__ && f_2__ > f_0__) {
	    goto L150;
	} else if (f_1__ > f_0__) {
	    syntrn_1.t += delta;
	    p = 1.;
	} else if (f_2__ > f_0__) {
	    syntrn_1.t -= delta;
	    p = -1.;
	    f_1__ = f_2__;
	} else {
	    syntrn_1.t += delta * .5 * (f_2__ - f_1__) / (f_1__ - f_0__ * 2. 
		    + f_2__);
	    goto L130;
	}
	spanned = FALSE_;
	while(! spanned) {
	    p *= 2.;
	    syntrn_1.t += p * delta;
	    if (syntrn_1.t <= 0.) {
		syntrn_1.t = 0.;
		f_2__ = energy1;
	    } else if (syntrn_1.t >= 1.) {
		syntrn_1.t = 1.;
		f_2__ = energy2;
	    } else {
		pathpnt_(&nvar, &syntrn_1.t, xx, syntrn_1.xmin1, 
			syntrn_1.xmin2);
		++ncalls;
		f_2__ = saddle1_(xx, g);
	    }
	    if (f_2__ > f_1__) {
		f_0__ = f_1__;
		f_1__ = f_2__;
	    } else {
		spanned = TRUE_;
	    }
	}
	p *= .5;
	syntrn_1.t -= p * delta;
	if (syntrn_1.t <= 0.) {
	    syntrn_1.t = 0.;
	    f_3__ = energy1;
	} else if (syntrn_1.t >= 1.) {
	    syntrn_1.t = 1.;
	    f_3__ = energy2;
	} else {
	    pathpnt_(&nvar, &syntrn_1.t, xx, syntrn_1.xmin1, syntrn_1.xmin2);
	    ++ncalls;
	    f_3__ = saddle1_(xx, g);
	}
	if (f_3__ > f_1__) {
	    syntrn_1.t += abs(p) * .5 * delta * (f_1__ - f_2__) / (f_2__ - 
		    f_3__ * 2. + f_1__);
	} else {
	    syntrn_1.t -= p * delta;
	    syntrn_1.t += abs(p) * .5 * delta * (f_0__ - f_3__) / (f_3__ - 
		    f_1__ * 2. + f_0__);
	}
L130:
	++niter;
	pathpnt_(&nvar, &syntrn_1.t, xx, syntrn_1.xmin1, syntrn_1.xmin2);
	ncalls += 3;
	f = saddle1_(xx, g);
	tangent_(&nvar, xx, g, &g_rms__, tan__, &g_tan__, &gamma, dgdt);
	io___69.ciunit = iounit_1.iout;
	s_wsfe(&io___69);
	do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&syntrn_1.t, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&g_rms__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&g_tan__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&gamma, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
	e_wsfe();
	if (ncycle >= maxcycle || gamma < gammamin) {
	    done = TRUE_;
	}
L150:
	++ncycle;
	delta *= epsilon;
    }

/*     if the path maximum is too near to an endpoint, */
/*     then negative curvature has probably been lost */

    if (syntrn_1.t <= .05 || syntrn_1.t >= .95) {
	if (! scan) {
	    pathscan_(&nvar, syntrn_1.xmin1, syntrn_1.xmin2, &ncalls);
	}
	io___70.ciunit = iounit_1.iout;
	s_wsfe(&io___70);
	e_wsfe();
	fatal_();
    }

/*     save the current maximum as the transition state estimate */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xx[i__ * 3 - 3];
	atoms_1.y[i__ - 1] = xx[i__ * 3 - 2];
	atoms_1.z__[i__ - 1] = xx[i__ * 3 - 1];
    }
    its = freeunit_();
    s_copy(tsfile, "tstate.xyz", (ftnlen)120, (ftnlen)10);
    version_(tsfile, "old", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = its;
    o__1.ofnmlen = 120;
    o__1.ofnm = tsfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    al__1.aerr = 0;
    al__1.aunit = its;
    f_rew(&al__1);
    prtxyz_(&its);
    cl__1.cerr = 0;
    cl__1.cunit = its;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     the maximum is located, get ready for minimization */

    sg = 0.;
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s0[i__ - 1] = tan__[i__ - 1];
	sg += s0[i__ - 1] * dgdt[i__ - 1];
    }
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h0[i__ - 1] = dgdt[i__ - 1] / sg;
    }

/*     set the initial conjugate direction for minimization */

    ninner = 0;
    g2 = 0.;
    hg = 0.;
    f_move__ = 1e6;
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = g[i__ - 1];
	g2 += d__1 * d__1;
	hg += h0[i__ - 1] * g[i__ - 1];
    }
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s[i__ - 1] = -g[i__ - 1] + hg * s0[i__ - 1];
    }
    g_rms__ = sqrt(g2 / (doublereal) atoms_1.n);
    io___79.ciunit = iounit_1.iout;
    s_wsfe(&io___79);
    e_wsfe();
    io___80.ciunit = iounit_1.iout;
    s_wsfe(&io___80);
    do_fio(&c__1, (char *)&ninner, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&g_rms__, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
    e_wsfe();

/*     check the termination criterion */

    if (g_rms__ < grdmin) {
	terminate = TRUE_;
	io___81.ciunit = iounit_1.iout;
	s_wsfe(&io___81);
	e_wsfe();
    }

/*     line search to find minimum in conjugate direction */

    while(! terminate) {
	++ninner;
	f_old__ = f;
	g2_old__ = g2;
	i__1 = nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x_old__[i__ - 1] = xx[i__ - 1];
	    g_old__[i__ - 1] = g[i__ - 1];
	}
	s_copy(status, "         ", (ftnlen)9, (ftnlen)9);
	linmin_1.angmax = 90.;
	search_(&nvar, &f, g, xx, s, &f_move__, &angle, &ncalls, (D_fp)
		saddle1_, status, (ftnlen)9);

/*     if search direction points uphill, use its negative */

	if (s_cmp(status, "WideAngle", (ftnlen)9, (ftnlen)9) == 0) {
	    i__1 = nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		s[i__ - 1] = -s[i__ - 1];
	    }
	    search_(&nvar, &f, g, xx, s, &f_move__, &angle, &ncalls, (D_fp)
		    saddle1_, status, (ftnlen)9);
	}

/*     compute movement and gradient following line search */

	f_move__ = f_old__ - f;
	x_move__ = 0.;
	g2 = 0.;
	i__1 = nvar;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = xx[i__ - 1] - x_old__[i__ - 1];
	    x_move__ += d__1 * d__1;
/* Computing 2nd power */
	    d__1 = g[i__ - 1];
	    g2 += d__1 * d__1;
	}
	x_move__ = sqrt(x_move__ / (doublereal) atoms_1.n);
	g_rms__ = sqrt(g2 / (doublereal) atoms_1.n);
	io___89.ciunit = iounit_1.iout;
	s_wsfe(&io___89);
	do_fio(&c__1, (char *)&ninner, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&f, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&g_rms__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&f_move__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&x_move__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&angle, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
	do_fio(&c__1, status, (ftnlen)9);
	e_wsfe();

/*     check the termination criteria */

	if (g_rms__ < grdmin) {
	    terminate = TRUE_;
	    io___90.ciunit = iounit_1.iout;
	    s_wsfe(&io___90);
	    e_wsfe();
	} else if (nouter >= maxouter) {
	    terminate = TRUE_;
	    io___91.ciunit = iounit_1.iout;
	    s_wsfe(&io___91);
	    e_wsfe();
	}

/*     check to see if another maximization is needed */

	if (! terminate) {
	    sg0 = 0.;
	    i__1 = nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		sg0 += s0[i__ - 1] * g[i__ - 1];
	    }
	    newcycle = FALSE_;
	    if (ninner >= maxinner) {
		newcycle = TRUE_;
	    }
	    if (sg0 * sg0 / g2 > diverge) {
		newcycle = TRUE_;
	    }
	    if (s_cmp(status, " Success ", (ftnlen)9, (ftnlen)9) != 0) {
		newcycle = TRUE_;
	    }

/*     unfortunately, a new maximization is needed; first save */
/*     the current minimum as the transition state estimate */

	    if (newcycle) {
		i__1 = atoms_1.n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    atoms_1.x[i__ - 1] = xx[i__ * 3 - 3];
		    atoms_1.y[i__ - 1] = xx[i__ * 3 - 2];
		    atoms_1.z__[i__ - 1] = xx[i__ * 3 - 1];
		}
		its = freeunit_();
		s_copy(tsfile, "tstate.xyz", (ftnlen)120, (ftnlen)10);
		version_(tsfile, "old", (ftnlen)120, (ftnlen)3);
		o__1.oerr = 0;
		o__1.ounit = its;
		o__1.ofnmlen = 120;
		o__1.ofnm = tsfile;
		o__1.orl = 0;
		o__1.osta = "old";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
		al__1.aerr = 0;
		al__1.aunit = its;
		f_rew(&al__1);
		prtxyz_(&its);
		cl__1.cerr = 0;
		cl__1.cunit = its;
		cl__1.csta = 0;
		f_clos(&cl__1);

/*     move the path endpoints toward current transition state; */
/*     then jump to the start of the next maximization cycle */

		if (reduce != 0.) {
		    pathval_(&nvar, xx);
		    syntrn_1.t = reduce * syntrn_1.pm;
		    pathpnt_(&nvar, &syntrn_1.t, x_old__, syntrn_1.xmin1, 
			    syntrn_1.xmin2);
		    i__1 = nvar;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			syntrn_1.xmin1[i__ - 1] = x_old__[i__ - 1];
		    }
		    ++ncalls;
		    energy1 = saddle1_(syntrn_1.xmin1, g);
		    syntrn_1.t = 1. - reduce * (1. - syntrn_1.pm);
		    pathpnt_(&nvar, &syntrn_1.t, x_old__, syntrn_1.xmin1, 
			    syntrn_1.xmin2);
		    i__1 = nvar;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			syntrn_1.xmin2[i__ - 1] = x_old__[i__ - 1];
		    }
		    ++ncalls;
		    energy2 = saddle1_(syntrn_1.xmin2, g);
		}
		goto L100;
	    }

/*     find the next conjugate search direction to search; */
/*     choice of "beta" is Fletcher-Reeves or Polak-Ribiere */

	    hg = 0.;
	    i__1 = nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		hg += h0[i__ - 1] * g[i__ - 1];
	    }
	    beta = 0.;
	    i__1 = nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*              beta = beta + g(i) * g(i) */
		beta += g[i__ - 1] * (g[i__ - 1] - g_old__[i__ - 1]);
	    }
	    beta /= g2_old__;
	    i__1 = nvar;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		s[i__ - 1] = -g[i__ - 1] + hg * s0[i__ - 1] + beta * s[i__ - 
			1];
	    }
	}
    }

/*     write out the final transition state structure */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xx[i__ * 3 - 3];
	atoms_1.y[i__ - 1] = xx[i__ * 3 - 2];
	atoms_1.z__[i__ - 1] = xx[i__ * 3 - 1];
    }
    its = freeunit_();
    s_copy(tsfile, "tstate.xyz", (ftnlen)120, (ftnlen)10);
    version_(tsfile, "old", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = its;
    o__1.ofnmlen = 120;
    o__1.ofnm = tsfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    al__1.aerr = 0;
    al__1.aunit = its;
    f_rew(&al__1);
    prtxyz_(&its);
    cl__1.cerr = 0;
    cl__1.cunit = its;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef keyline_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine pathval  --  synchronous transit path values  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "pathval" computes the synchronous transit path value for */
/*     the specified structure */


/* Subroutine */ int pathval_(integer *nvar, doublereal *xx)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal rmsvalue;
    static integer i__;
    static doublereal x1[25000], x2[25000], y1[25000], y2[25000], z1[25000], 
	    z2[25000], dp, dr;
    extern /* Subroutine */ int impose_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *);



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
/*     ##  syntrn.i  --  definition of synchronous transit path  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     t       value of the path coordinate (0=reactant, 1=product) */
/*     pm      path coordinate for extra point in quadratic transit */
/*     xmin1   reactant coordinates as array of optimization variables */
/*     xmin2   product coordinates as array of optimization variables */
/*     xm      extra coordinate set for quadratic synchronous transit */




/*     find the value of the transit path coordinate "pm"; */
/*     it is the ratio of the rms fits to the two endpoints */

    /* Parameter adjustments */
    --xx;

    /* Function Body */
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x1[i__ - 1] = syntrn_1.xmin1[i__ * 3 - 3];
	y1[i__ - 1] = syntrn_1.xmin1[i__ * 3 - 2];
	z1[i__ - 1] = syntrn_1.xmin1[i__ * 3 - 1];
	x2[i__ - 1] = xx[i__ * 3 - 2];
	y2[i__ - 1] = xx[i__ * 3 - 1];
	z2[i__ - 1] = xx[i__ * 3];
    }
    impose_(&atoms_1.n, x1, y1, z1, &atoms_1.n, x2, y2, z2, &dr);
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x1[i__ - 1] = syntrn_1.xmin2[i__ * 3 - 3];
	y1[i__ - 1] = syntrn_1.xmin2[i__ * 3 - 2];
	z1[i__ - 1] = syntrn_1.xmin2[i__ * 3 - 1];
	x2[i__ - 1] = xx[i__ * 3 - 2];
	y2[i__ - 1] = xx[i__ * 3 - 1];
	z2[i__ - 1] = xx[i__ * 3];
    }
    impose_(&atoms_1.n, x1, y1, z1, &atoms_1.n, x2, y2, z2, &dp);
    syntrn_1.pm = dr / (dr + dp);

/*     superimpose on linear transit structure of same path value */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x1[i__ - 1] = (1. - syntrn_1.pm) * syntrn_1.xmin1[i__ * 3 - 3] + 
		syntrn_1.pm * syntrn_1.xmin2[i__ * 3 - 3];
	y1[i__ - 1] = (1. - syntrn_1.pm) * syntrn_1.xmin1[i__ * 3 - 2] + 
		syntrn_1.pm * syntrn_1.xmin2[i__ * 3 - 2];
	z1[i__ - 1] = (1. - syntrn_1.pm) * syntrn_1.xmin1[i__ * 3 - 1] + 
		syntrn_1.pm * syntrn_1.xmin2[i__ * 3 - 1];
	x2[i__ - 1] = xx[i__ * 3 - 2];
	y2[i__ - 1] = xx[i__ * 3 - 1];
	z2[i__ - 1] = xx[i__ * 3];
    }
    impose_(&atoms_1.n, x1, y1, z1, &atoms_1.n, x2, y2, z2, &rmsvalue);
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xx[i__ * 3 - 2] = x2[i__ - 1];
	xx[i__ * 3 - 1] = y2[i__ - 1];
	xx[i__ * 3] = z2[i__ - 1];
    }
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	syntrn_1.xm[i__ - 1] = xx[i__];
    }
    return 0;
} /* pathval_ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine pathpnt  --  get coordinates of path point  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "pathpnt" finds a structure on the synchronous transit path */
/*     with the specified path value "t" */


/* Subroutine */ int pathpnt_(integer *nvar, doublereal *t, doublereal *xx, 
	doublereal *x0, doublereal *x1)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int ocvm_(integer *, doublereal *, doublereal *, 
	    doublereal *, D_fp, U_fp);
    static doublereal value, grdmin;
    extern /* Subroutine */ int optsave_();
    extern doublereal transit_();



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

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  minima.i  --  general parameters for minimizations  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     fctmin    value below which function is deemed optimized */
/*     hguess    initial value for the H-matrix diagonal elements */
/*     maxiter   maximum number of iterations during optimization */
/*     nextiter  iteration number to use for the first iteration */




/*     initialize some parameters for the upcoming optimization */

    /* Parameter adjustments */
    --x1;
    --x0;
    --xx;

    /* Function Body */
    if (inform_1.debug) {
	inform_1.iprint = 1;
    } else {
	inform_1.iprint = 0;
    }
    inform_1.iwrite = 0;
    minima_1.maxiter = 1000;
    grdmin = 1e-5;

/*     interpolate coordinates to give initial estimate */

    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xx[i__] = (1. - *t) * x0[i__] + *t * x1[i__];
    }

/*     optimize the synchronous transit function */

/*     call lbfgs (nvar,xx,value,grdmin,transit,optsave) */
    ocvm_(nvar, &xx[1], &value, &grdmin, (D_fp)transit_, (U_fp)optsave_);
    return 0;
} /* pathpnt_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine pathscan  --  scan along the transit pathway  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "pathscan" makes a scan of a synchronous transit pathway by */
/*     computing structures and energies for specific path values */


/* Subroutine */ int pathscan_(integer *nvar, doublereal *x0, doublereal *x1, 
	integer *ncalls)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Scan of the Synchronous Transit Pathwa"
	    "y :\002,/,\002 N Scan     F Value       Path      RMS G\002,\002"
	    "      G Tan      Gamma   FG Call\002,/)";
    static char fmt_20[] = "(i6,f13.4,f11.4,f11.4,f11.4,f11.5,i8)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal g[75000];
    static integer i__;
    static doublereal xx[75000], tan__[75000], dgdt[75000], gamma, g_tan__, 
	    g_rms__, energy;
    extern doublereal saddle1_(doublereal *, doublereal *);
    extern /* Subroutine */ int tangent_(integer *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), pathpnt_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___108 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___118 = { 0, 0, 0, fmt_20, 0 };




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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  syntrn.i  --  definition of synchronous transit path  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     t       value of the path coordinate (0=reactant, 1=product) */
/*     pm      path coordinate for extra point in quadratic transit */
/*     xmin1   reactant coordinates as array of optimization variables */
/*     xmin2   product coordinates as array of optimization variables */
/*     xm      extra coordinate set for quadratic synchronous transit */




/*     make a scan along the synchronous transit pathway */

    /* Parameter adjustments */
    --x1;
    --x0;

    /* Function Body */
    io___108.ciunit = iounit_1.iout;
    s_wsfe(&io___108);
    e_wsfe();
    for (i__ = 0; i__ <= 10; ++i__) {
	syntrn_1.t = (doublereal) i__ * .1;
	pathpnt_(nvar, &syntrn_1.t, xx, &x0[1], &x1[1]);
	*ncalls += 3;
	energy = saddle1_(xx, g);
	tangent_(nvar, xx, g, &g_rms__, tan__, &g_tan__, &gamma, dgdt);
	io___118.ciunit = iounit_1.iout;
	s_wsfe(&io___118);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&energy, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&syntrn_1.t, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&g_rms__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&g_tan__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&gamma, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*ncalls), (ftnlen)sizeof(integer));
	e_wsfe();
    }
    return 0;
} /* pathscan_ */



/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine tangent  --  synchronous transit tangent  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "tangent" finds the projected gradient on the synchronous */
/*     transit path for a point along the transit pathway */


/* Subroutine */ int tangent_(integer *nvar, doublereal *xx, doublereal *g, 
	doublereal *g_rms__, doublereal *tan__, doublereal *g_tan__, 
	doublereal *gamma, doublereal *dgdt)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal tan_norm__;
    static integer i__;
    static doublereal g2, t0, gb[75000], gf[75000], xb[75000], xf[75000], 
	    delta, energy;
    extern doublereal saddle1_(doublereal *, doublereal *);
    extern /* Subroutine */ int pathpnt_(integer *, doublereal *, doublereal *
	    , doublereal *, doublereal *);



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
/*     ##  syntrn.i  --  definition of synchronous transit path  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     t       value of the path coordinate (0=reactant, 1=product) */
/*     pm      path coordinate for extra point in quadratic transit */
/*     xmin1   reactant coordinates as array of optimization variables */
/*     xmin2   product coordinates as array of optimization variables */
/*     xm      extra coordinate set for quadratic synchronous transit */




/*     set the finite difference path increment */

    /* Parameter adjustments */
    --dgdt;
    --tan__;
    --g;
    --xx;

    /* Function Body */
    delta = .01;

/*     store the initial pathpnt and compute gradient norm */

    t0 = syntrn_1.t;
    g2 = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = g[i__];
	g2 += d__1 * d__1;
    }
    *g_rms__ = sqrt(g2 / (doublereal) atoms_1.n);

/*     compute the forward difference */

    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xf[i__ - 1] = xx[i__];
    }
    syntrn_1.t = t0 + delta;
    pathpnt_(nvar, &syntrn_1.t, xf, xf, xf);
    energy = saddle1_(xf, gf);

/*     compute the backward difference */

    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xb[i__ - 1] = xx[i__];
    }
    syntrn_1.t = t0 - delta;
    pathpnt_(nvar, &syntrn_1.t, xb, xb, xb);
    energy = saddle1_(xb, gb);
    syntrn_1.t = t0;

/*     compute tangent to the path, and projected gradient */

    tan_norm__ = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tan__[i__] = xf[i__ - 1] - xb[i__ - 1];
/* Computing 2nd power */
	d__1 = tan__[i__];
	tan_norm__ += d__1 * d__1;
	dgdt[i__] = gf[i__ - 1] - gb[i__ - 1];
    }
    tan_norm__ = sqrt(tan_norm__);
    *g_tan__ = 0.;
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tan__[i__] /= tan_norm__;
	*g_tan__ += g[i__] * tan__[i__];
    }
    *g_tan__ /= sqrt((doublereal) atoms_1.n);
/* Computing 2nd power */
    d__1 = *g_tan__ / *g_rms__;
    *gamma = d__1 * d__1;
    return 0;
} /* tangent_ */



/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  function transit  --  synchronous transit evaluation  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "transit" evaluates the synchronous transit function and */
/*     gradient; linear and quadratic transit paths are available */


doublereal transit_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal r1, r2, rc, rd, wc, ri, wd, rm, wi;
    static integer ix, iy, iz, jx, jy, jz;
    static doublereal tq, pq, x1d, y1d, z1d, x2d, y2d, x1i, y1i, z1i, x2i, 
	    y2i, z2i, z2d, ri4, xcd, ycd, zcd, xci, yci, zci, xmd, ymd, zmd, 
	    xmi, ymi, zmi;
    static char mode[9];
    static integer nvar;
    static doublereal term, gamma, value, termx, termy, termz, cutoff, 
	    cutoff2;



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
/*     ##  syntrn.i  --  definition of synchronous transit path  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     t       value of the path coordinate (0=reactant, 1=product) */
/*     pm      path coordinate for extra point in quadratic transit */
/*     xmin1   reactant coordinates as array of optimization variables */
/*     xmin2   product coordinates as array of optimization variables */
/*     xm      extra coordinate set for quadratic synchronous transit */




/*     zero out the synchronous transit function and gradient */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    value = 0.;
    nvar = atoms_1.n * 3;
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = 0.;
    }
    tq = 1. - syntrn_1.t;

/*     set the cutoff distance for interatomic distances */

    cutoff = 1e3;
/* Computing 2nd power */
    d__1 = cutoff;
    cutoff2 = d__1 * d__1;

/*     set the type of synchronous transit path to be used */

    if (syntrn_1.pm == 0.) {
	s_copy(mode, "LINEAR", (ftnlen)9, (ftnlen)6);
    } else {
	s_copy(mode, "QUADRATIC", (ftnlen)9, (ftnlen)9);
	pq = 1. - syntrn_1.pm;
    }

/*     portion based on interpolated interatomic distances */

    i__1 = atoms_1.n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iz = i__ * 3;
	iy = iz - 1;
	ix = iz - 2;
	xci = xx[ix];
	yci = xx[iy];
	zci = xx[iz];
	x1i = syntrn_1.xmin1[ix - 1];
	y1i = syntrn_1.xmin1[iy - 1];
	z1i = syntrn_1.xmin1[iz - 1];
	x2i = syntrn_1.xmin2[ix - 1];
	y2i = syntrn_1.xmin2[iy - 1];
	z2i = syntrn_1.xmin2[iz - 1];
	if (s_cmp(mode, "QUADRATIC", (ftnlen)9, (ftnlen)9) == 0) {
	    xmi = syntrn_1.xm[ix - 1];
	    ymi = syntrn_1.xm[iy - 1];
	    zmi = syntrn_1.xm[iz - 1];
	}
	i__2 = atoms_1.n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    jz = j * 3;
	    jy = jz - 1;
	    jx = jz - 2;
	    xcd = xci - xx[jx];
	    ycd = yci - xx[jy];
	    zcd = zci - xx[jz];
	    x1d = x1i - syntrn_1.xmin1[jx - 1];
	    y1d = y1i - syntrn_1.xmin1[jy - 1];
	    z1d = z1i - syntrn_1.xmin1[jz - 1];
	    x2d = x2i - syntrn_1.xmin2[jx - 1];
	    y2d = y2i - syntrn_1.xmin2[jy - 1];
	    z2d = z2i - syntrn_1.xmin2[jz - 1];
/* Computing 2nd power */
	    d__1 = xcd;
/* Computing 2nd power */
	    d__2 = ycd;
/* Computing 2nd power */
	    d__3 = zcd;
	    rc = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
/* Computing 2nd power */
	    d__1 = x1d;
/* Computing 2nd power */
	    d__2 = y1d;
/* Computing 2nd power */
	    d__3 = z1d;
	    r1 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
/* Computing 2nd power */
	    d__1 = x2d;
/* Computing 2nd power */
	    d__2 = y2d;
/* Computing 2nd power */
	    d__3 = z2d;
	    r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
/* Computing MIN */
	    d__1 = min(rc,r1);
	    if (min(d__1,r2) < cutoff2) {
		rc = sqrt(rc);
		r1 = sqrt(r1);
		r2 = sqrt(r2);
		ri = tq * r1 + syntrn_1.t * r2;
		if (s_cmp(mode, "QUADRATIC", (ftnlen)9, (ftnlen)9) == 0) {
		    xmd = xmi - syntrn_1.xm[jx - 1];
		    ymd = ymi - syntrn_1.xm[jy - 1];
		    zmd = zmi - syntrn_1.xm[jz - 1];
/* Computing 2nd power */
		    d__1 = xmd;
/* Computing 2nd power */
		    d__2 = ymd;
/* Computing 2nd power */
		    d__3 = zmd;
		    rm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
		    gamma = (rm - pq * r1 - syntrn_1.pm * r2) / (syntrn_1.pm *
			     pq);
		    ri += gamma * syntrn_1.t * tq;
		}
/* Computing 4th power */
		d__1 = ri, d__1 *= d__1;
		ri4 = d__1 * d__1;
		rd = rc - ri;
/* Computing 2nd power */
		d__1 = rd;
		value += d__1 * d__1 / ri4;
		term = rd * 2. / (ri4 * rc);
		termx = term * xcd;
		termy = term * ycd;
		termz = term * zcd;
		g[ix] += termx;
		g[iy] += termy;
		g[iz] += termz;
		g[jx] -= termx;
		g[jy] -= termy;
		g[jz] -= termz;
	    }
	}
    }

/*     portion used to supress rigid rotations and translations */

    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wc = xx[i__];
	wi = tq * syntrn_1.xmin1[i__ - 1] + syntrn_1.t * syntrn_1.xmin2[i__ - 
		1];
	wd = wc - wi;
/* Computing 2nd power */
	d__1 = wd;
	value += d__1 * d__1 * 1e-6;
	g[i__] += wd * 2e-6;
    }
    ret_val = value;
    return ret_val;
} /* transit_ */



/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  function saddle1  --  energy and gradient for saddle  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "saddle1" is a service routine that computes the energy and */
/*     gradient for transition state optimization */


doublereal saddle1_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static doublereal e;
    static integer i__;



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




    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = xx[i__ * 3 - 2];
	atoms_1.y[i__ - 1] = xx[i__ * 3 - 1];
	atoms_1.z__[i__ - 1] = xx[i__ * 3];
    }
    gradient_(&e, &g[1]);
    ret_val = e;
    return ret_val;
} /* saddle1_ */

/* Main program alias */ int saddle_ () { MAIN__ (); return 0; }
