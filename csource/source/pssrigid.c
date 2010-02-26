/* pssrigid.f -- translated by f2c (version 20050501).
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
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

struct {
    doublereal xref[250000]	/* was [25000][10] */, yref[250000]	/* 
	    was [25000][10] */, zref[250000]	/* was [25000][10] */;
    integer nref[10], reftyp[250000]	/* was [25000][10] */, n12ref[250000]	
	    /* was [25000][10] */, i12ref[2000000]	/* was [8][25000][10] 
	    */, refleng[10], refltitle[10];
    char refnam[750000]	/* was [25000][10] */, reffile[1200], reftitle[1200];
} refer_;

#define refer_1 refer_

struct {
    doublereal xrb[25000], yrb[25000], zrb[25000], rbc[6000]	/* was [6][
	    1000] */;
    logical use_rigid__;
} rigid_;

#define rigid_1 rigid_

struct {
    doublereal m2[25000], deform, difft, diffv, diffc;
    logical use_smooth__, use_dem__, use_gda__, use_tophat__, use_stophat__;
} warp_;

#define warp_1 warp_

struct {
    doublereal fctmin, hguess;
    integer maxiter, nextiter;
} minima_;

#define minima_1 minima_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static doublereal c_b41 = 12.;
static integer c__6000 = 6000;



/*     ############################################################## */
/*     ##  COPYRIGHT (C) 1997 by Rohit Pappu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  program pssrigid  --  smoothing & search over rigid bodies  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "pssrigid" implements the potential smoothing plus search method */
/*     for global optimization for a set of rigid bodies */

/*     literature reference: */

/*     J. Kostrowicki and H. A. Scheraga, "Application of the Diffusion */
/*     Equation Method for Global Optimization to Oligopeptides", Journal */
/*     of Physical Chemistry, 96, 7442-7449 (1992) */

/*     S. Nakamura, H. Hirose, M. Ikeguchi and J. Doi, "Conformational */
/*     Energy Minimization Using a Two-Stage Method", Journal of Physical */
/*     Chemistry, 99, 8374-8378 (1995) */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Enter the Number of Steps for the PSS Sc"
	    "hedule\002,\002 [100] :  \002,$)";
    static char fmt_30[] = "(i10)";
    static char fmt_40[] = "(/,\002 Use Local Search to Explore Each Smoothi"
	    "ng Level\002,\002 [N] :  \002,$)";
    static char fmt_50[] = "(a120)";
    static char fmt_70[] = "(/,\002 Enter the Number of Directions for Loca"
	    "l\002,\002 Search [\002,i2,\002] :  \002,$)";
    static char fmt_80[] = "(i10)";
    static char fmt_100[] = "(/,\002 Enter the Largest Smoothing Value for L"
	    "ocal\002,\002 Search [5.0] :  \002,$)";
    static char fmt_110[] = "(f20.0)";
    static char fmt_130[] = "(/,\002 Enter RMS Gradient per Rigid Body Crite"
	    "rion\002,\002 [0.0001] :  \002,$)";
    static char fmt_140[] = "(f20.0)";
    static char fmt_150[] = "(/,\002 Final Function Value and Deformation "
	    ":\002,2f15.4)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2, i__3[3];
    doublereal d__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    extern integer freeunit_(void);
    extern /* Subroutine */ int rigidxyz_(void);
    static integer i__, j, k;
    static logical use_local__;
    static doublereal xx[75000];
    static char ext[7];
    static doublereal rms;
    static integer nvar;
    extern /* Subroutine */ int ocvm_(integer *, doublereal *, doublereal *, 
	    doublereal *, D_fp, U_fp);
    static integer lext, next, ixyz;
    extern /* Subroutine */ int final_(void);
    static doublereal ratio;
    static logical exist;
    static integer neigen;
    static char record[120];
    static doublereal grdmin;
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char answer[1];
    static integer npoint;
    static char string[120];
    extern /* Subroutine */ int orient_(void), impose_(integer *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *), getxyz_(void);
    static doublereal deform0;
    extern /* Subroutine */ int prtxyz_(integer *);
    extern doublereal pssrgd1_();
    extern /* Subroutine */ int makeref_(integer *), initial_(void);
    extern doublereal sigmoid_(doublereal *, doublereal *);
    extern /* Subroutine */ int rgdsrch_(integer *, doublereal *, doublereal *
	    ), numeral_(integer *, char *, integer *, ftnlen);
    static doublereal srchmax;
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen);
    static doublereal minimum;
    extern /* Subroutine */ int optsave_();
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static icilist io___5 = { 1, string, 1, 0, 120, 1 };
    static cilist io___6 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___15 = { 1, string, 1, 0, 120, 1 };
    static cilist io___17 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_80, 0 };
    static icilist io___20 = { 1, string, 1, 0, 120, 1 };
    static cilist io___21 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_110, 0 };
    static icilist io___24 = { 1, string, 1, 0, 120, 1 };
    static cilist io___25 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_150, 0 };



#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]



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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  refer.i  --  storage of reference atomic coordinate set  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xref        reference x-coordinates for atoms in each system */
/*     yref        reference y-coordinates for atoms in each system */
/*     zref        reference z-coordinates for atoms in each system */
/*     nref        total number of atoms in each reference system */
/*     reftyp      atom types of the atoms in each reference system */
/*     n12ref      number of atoms bonded to each reference atom */
/*     i12ref      atom numbers of atoms 1-2 connected to each atom */
/*     refleng     length in characters of each reference filename */
/*     refltitle   length in characters of each reference title line */
/*     refnam      atom names of the atoms in each reference system */
/*     reffile     base filename for each reference system */
/*     reftitle    title used to describe each reference system */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  warp.i  --  parameters for potential surface smoothing  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     m2           second moment of the GDA gaussian for each atom */
/*     deform       value of the smoothing deformation parameter */
/*     difft        diffusion coefficient for torsional potential */
/*     diffv        diffusion coefficient for van der Waals potential */
/*     diffc        diffusion coefficient for charge-charge potential */
/*     use_smooth   flag to use a potential energy smoothing method */
/*     use_dem      flag to use diffusion equation method potential */
/*     use_gda      flag to use gaussian density annealing potential */
/*     use_tophat   flag to use analytical tophat smoothed potential */
/*     use_stophat  flag to use shifted tophat smoothed potential */




/*     set up the structure, mechanics calculation and smoothing */

    initial_();
    getxyz_();
    warp_1.use_smooth__ = TRUE_;
    warp_1.use_dem__ = TRUE_;
    mechanic_();

/*     get rigid body coordinates and save the Cartesian coordinates */

    rigid_1.use_rigid__ = TRUE_;
    orient_();
    makeref_(&c__1);

/*     set maximum deformation value and disable coordinate dumps */

    deform0 = warp_1.deform;
    inform_1.iwrite = 0;

/*     get the number of points along the deformation schedule */

    npoint = -1;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___5);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&npoint, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    }
L10:
    if (npoint < 0) {
	io___6.ciunit = iounit_1.iout;
	s_wsfe(&io___6);
	e_wsfe();
	io___7.ciunit = iounit_1.input;
	s_rsfe(&io___7);
	do_fio(&c__1, (char *)&npoint, (ftnlen)sizeof(integer));
	e_rsfe();
	if (npoint <= 0) {
	    npoint = 100;
	}
    }

/*     decide whether to use the local search procedure */

    use_local__ = FALSE_;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___10.ciunit = iounit_1.iout;
	s_wsfe(&io___10);
	e_wsfe();
	io___11.ciunit = iounit_1.input;
	s_rsfe(&io___11);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'Y') {
	use_local__ = TRUE_;
    }

/*     get the number of eigenvectors to use for the local search */

    if (use_local__) {
	neigen = -1;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___15);
	    if (i__1 != 0) {
		goto L60;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&neigen, (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L60;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L60;
	    }
	}
L60:
	if (neigen <= 0) {
	    nvar = (group_1.ngrp - 1) * 6;
	    io___17.ciunit = iounit_1.iout;
	    s_wsfe(&io___17);
	    do_fio(&c__1, (char *)&nvar, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___18.ciunit = iounit_1.input;
	    s_rsfe(&io___18);
	    do_fio(&c__1, (char *)&neigen, (ftnlen)sizeof(integer));
	    e_rsfe();
	    if (neigen > nvar) {
		neigen = nvar;
	    }
	}
    }

/*     get the maximal smoothing level for use of local search */

    if (use_local__) {
	srchmax = -1.;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___20);
	    if (i__1 != 0) {
		goto L90;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&srchmax, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L90;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L90;
	    }
	}
L90:
	if (srchmax < 0.) {
	    io___21.ciunit = iounit_1.iout;
	    s_wsfe(&io___21);
	    e_wsfe();
	    io___22.ciunit = iounit_1.input;
	    s_rsfe(&io___22);
	    do_fio(&c__1, (char *)&srchmax, (ftnlen)sizeof(doublereal));
	    e_rsfe();
	    if (srchmax < 0.) {
		srchmax = 5.;
	    }
	}
    }

/*     get the termination criterion as RMS gradient per atom */

    grdmin = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___24);
	if (i__1 != 0) {
	    goto L120;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&grdmin, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L120;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L120;
	}
    }
L120:
    if (grdmin <= 0.) {
	io___25.ciunit = iounit_1.iout;
	s_wsfe(&io___25);
	e_wsfe();
	io___26.ciunit = iounit_1.input;
	s_rsfe(&io___26);
	do_fio(&c__1, (char *)&grdmin, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (grdmin <= 0.) {
	grdmin = 1e-4;
    }

/*     perform PSS iteration by looping over smoothed surfaces */

    i__1 = npoint << 1;
    for (k = 0; k <= i__1; ++k) {
	ratio = 1. - (i__2 = npoint - k, (doublereal) abs(i__2)) / (
		doublereal) npoint;
	if (molcul_1.nmol == 1) {
/* Computing 3rd power */
	    d__1 = ratio;
	    warp_1.deform = deform0 * (d__1 * (d__1 * d__1));
	} else {
	    warp_1.deform = deform0 * sigmoid_(&c_b41, &ratio);
	}

/*     transfer rigid body coordinates to optimization parameters */

	nvar = 0;
	i__2 = group_1.ngrp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    for (j = 1; j <= 6; ++j) {
		++nvar;
		xx[nvar - 1] = rbc_ref(j, i__);
	    }
	}

/*     make the call to the variable metric optimization routine */

	inform_1.iprint = 1;
	ocvm_(&nvar, xx, &minimum, &grdmin, (D_fp)pssrgd1_, (U_fp)optsave_);

/*     transfer optimization parameters to rigid body coordinates */

	nvar = 0;
	i__2 = group_1.ngrp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    for (j = 1; j <= 6; ++j) {
		++nvar;
		rbc_ref(j, i__) = xx[nvar - 1];
	    }
	}

/*     use normal mode local search to explore adjacent minima */

	if (use_local__) {
	    if (warp_1.deform <= srchmax && k >= npoint) {
		rgdsrch_(&neigen, &minimum, &grdmin);
	    }
	}

/*     write out final energy function value and smoothing level */

	io___33.ciunit = iounit_1.iout;
	s_wsfe(&io___33);
	do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&warp_1.deform, (ftnlen)sizeof(doublereal));
	e_wsfe();

/*     get Cartesian coordinates and superimpose on reference */

	rigidxyz_();
	if (igrp_ref(1, 1) == 1 && igrp_ref(2, group_1.ngrp) == atoms_1.n) {
	    impose_(&atoms_1.n, refer_1.xref, refer_1.yref, refer_1.zref, &
		    atoms_1.n, atoms_1.x, atoms_1.y, atoms_1.z__, &rms);
	}

/*     write the coordinates of the current minimum to a file */

	lext = 3;
	numeral_(&k, ext, &lext, (ftnlen)7);
	ixyz = freeunit_();
/* Writing concatenation */
	i__3[0] = files_1.leng, a__1[0] = files_1.filename;
	i__3[1] = 1, a__1[1] = ".";
	i__3[2] = lext, a__1[2] = ext;
	s_cat(xyzfile, a__1, i__3, &c__3, (ftnlen)120);
	version_(xyzfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = ixyz;
	o__1.ofnmlen = 120;
	o__1.ofnm = xyzfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	prtxyz_(&ixyz);
	cl__1.cerr = 0;
	cl__1.cunit = ixyz;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef igrp_ref
#undef rbc_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function pssrgd1  --  energy and gradient values for PSS  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "pssrgd1" is a service routine that computes the energy and */
/*     gradient during PSS global optimization over rigid bodies */


doublereal pssrgd1_(doublereal *xx, doublereal *g)
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
} /* pssrgd1_ */

#undef derivs_ref
#undef rbc_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine rgdsrch  --  local search for rigid body PSS  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/* Subroutine */ int rgdsrch_(integer *neigen, doublereal *minimum, 
	doublereal *grdmin)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Normal Mode Search :\002,8x,\002Iterat"
	    "ion\002,i4,6x,\002Energy\002,f12.4,/)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int climbrgd_(integer *, doublereal *, doublereal 
	    *, doublereal *), eigenrgd_(doublereal *, doublereal *), 
	    rigidxyz_(void);
    static integer i__, j, k;
    static doublereal eps;
    static logical done;
    static integer ndoi, nvar;
    static doublereal step[6000], eigen[6000], rbref[6000]	/* was [6][
	    1000] */, vects[36000000]	/* was [6000][6000] */, minref, 
	    rbbest[6000]	/* was [6][1000] */;
    static integer nsearch;
    static doublereal minbest;

    /* Fortran I/O blocks */
    static cilist io___53 = { 0, 0, 0, fmt_10, 0 };



#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]
#define rbref_ref(a_1,a_2) rbref[(a_2)*6 + a_1 - 7]
#define vects_ref(a_1,a_2) vects[(a_2)*6000 + a_1 - 6001]
#define rbbest_ref(a_1,a_2) rbbest[(a_2)*6 + a_1 - 7]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




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




/*     set parameters related to the local search procedure */

    done = FALSE_;
    eps = 1e-4;
    minref = *minimum;
    minbest = *minimum;
    ndoi = 0;
    nvar = group_1.ngrp * 6;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    rbref_ref(j, i__) = rbc_ref(j, i__);
	}
    }

/*     find local minimum along each of the steepest directions */

    while(! done) {
	++ndoi;
	io___53.ciunit = iounit_1.iout;
	s_wsfe(&io___53);
	do_fio(&c__1, (char *)&ndoi, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&minref, (ftnlen)sizeof(doublereal));
	e_wsfe();
	rigidxyz_();
	eigenrgd_(eigen, vects);

/*     search both directions along each eigenvector in turn */

	nsearch = 0;
	i__1 = *neigen;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		step[k - 1] = vects_ref(k, nvar - i__ + 1);
	    }
	    i__2 = group_1.ngrp;
	    for (k = 1; k <= i__2; ++k) {
		for (j = 1; j <= 6; ++j) {
		    rbc_ref(j, k) = rbref_ref(j, k);
		}
	    }
	    ++nsearch;
	    climbrgd_(&nsearch, minimum, step, grdmin);
	    if (*minimum < minbest) {
		minbest = *minimum;
		i__2 = group_1.ngrp;
		for (k = 1; k <= i__2; ++k) {
		    for (j = 1; j <= 6; ++j) {
			rbbest_ref(j, k) = rbc_ref(j, k);
		    }
		}
	    }
	    i__2 = nvar;
	    for (k = 1; k <= i__2; ++k) {
		step[k - 1] = -vects_ref(k, nvar - i__ + 1);
	    }
	    i__2 = group_1.ngrp;
	    for (k = 1; k <= i__2; ++k) {
		for (j = 1; j <= 6; ++j) {
		    rbc_ref(j, k) = rbref_ref(j, k);
		}
	    }
	    ++nsearch;
	    climbrgd_(&nsearch, minimum, step, grdmin);
	    if (*minimum < minbest) {
		minbest = *minimum;
		i__2 = group_1.ngrp;
		for (k = 1; k <= i__2; ++k) {
		    for (j = 1; j <= 6; ++j) {
			rbbest_ref(j, k) = rbc_ref(j, k);
		    }
		}
	    }
	}

/*     check for convergence of the local search procedure */

	if (minbest < minref - eps) {
	    done = FALSE_;
	    minref = minbest;
	    i__1 = group_1.ngrp;
	    for (k = 1; k <= i__1; ++k) {
		for (j = 1; j <= 6; ++j) {
		    rbref_ref(j, k) = rbbest_ref(j, k);
		}
	    }
	} else {
	    done = TRUE_;
	    *minimum = minref;
	    i__1 = group_1.ngrp;
	    for (k = 1; k <= i__1; ++k) {
		for (j = 1; j <= 6; ++j) {
		    rbc_ref(j, k) = rbref_ref(j, k);
		}
	    }
	}
    }
    return 0;
} /* rgdsrch_ */

#undef rbbest_ref
#undef vects_ref
#undef rbref_ref
#undef rbc_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine eigenrgd  --  rigid body Hessian eigenvectors  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/* Subroutine */ int eigenrgd_(doublereal *eigen, doublereal *vects)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal a[6000], b[6000];
    static integer i__, j;
    static doublereal p[6000], w[6000], ta[6000], tb[6000], ty[6000];
    static integer nvar;
    extern /* Subroutine */ int diagq_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer ihess;
    static doublereal vnorm, hrigid[36000000]	/* was [6000][6000] */, 
	    matrix[18003000];
    extern /* Subroutine */ int hessrgd_(doublereal *);


#define vects_ref(a_1,a_2) vects[(a_2)*6000 + a_1]
#define hrigid_ref(a_1,a_2) hrigid[(a_2)*6000 + a_1 - 6001]



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




/*     compute the Hessian for rigid body motion */

    /* Parameter adjustments */
    vects -= 6001;
    --eigen;

    /* Function Body */
    hessrgd_(hrigid);

/*     place Hessian elements into triangular form */

    nvar = group_1.ngrp * 6;
    ihess = 0;
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nvar;
	for (j = i__; j <= i__2; ++j) {
	    ++ihess;
	    matrix[ihess - 1] = hrigid_ref(i__, j);
	}
    }

/*     diagonalize the Hessian to obtain eigenvalues */

    diagq_(&nvar, &c__6000, &nvar, matrix, &eigen[1], &vects[6001], a, b, p, 
	    w, ta, tb, ty);

/*     normalize the rigid body Hessian eigenvectors */

    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vnorm = 0.;
	i__2 = nvar;
	for (j = 1; j <= i__2; ++j) {
/* Computing 2nd power */
	    d__1 = vects_ref(j, i__);
	    vnorm += d__1 * d__1;
	}
	vnorm = sqrt(vnorm);
	i__2 = nvar;
	for (j = 1; j <= i__2; ++j) {
	    vects_ref(j, i__) = vects_ref(j, i__) / vnorm;
	}
    }
    return 0;
} /* eigenrgd_ */

#undef hrigid_ref
#undef vects_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine climbrgd  --  minimum from a PSS local search  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/* Subroutine */ int climbrgd_(integer *nsearch, doublereal *minimum, 
	doublereal *step, doublereal *grdmin)
{
    /* Format strings */
    static char fmt_10[] = "(4x,\002Search Direction\002,i4,10x,\002Step\002"
	    ",i6,10x,f12.4)";
    static char fmt_20[] = "(4x,\002Search Direction\002,i4,36x,\002-----"
	    "-\002)";
    static char fmt_30[] = "(4x,\002Search Direction\002,i4,36x,\002-----"
	    "-\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int localrgd_(doublereal *, doublereal *), 
	    rigidxyz_(void);
    static integer i__, j;
    static doublereal big;
    static logical done;
    static integer nvar;
    static doublereal size, estep[501];
    static integer kstep, nstep;
    extern doublereal energy_(void);

    /* Fortran I/O blocks */
    static cilist io___83 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___85 = { 0, 0, 0, fmt_30, 0 };



#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]



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




/*     set the maximum number of steps and the step size */

    /* Parameter adjustments */
    --step;

    /* Function Body */
    done = FALSE_;
    big = 1e5;
    *minimum = big;
    kstep = 0;
    nstep = 65;
/*     size = 0.1d0 */
    size = 1.;
    nvar = group_1.ngrp * 6;
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	step[i__] = size * step[i__];
    }

/*     scan the search direction for a minimization candidate */

    while(! done) {
	if (kstep != 0) {
	    nvar = 0;
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 6; ++j) {
		    ++nvar;
		    rbc_ref(j, i__) = rbc_ref(j, i__) + step[nvar];
		}
	    }
	}
	rigidxyz_();
	estep[kstep] = energy_();
	if (kstep >= 2 && estep[kstep] <= 1e4) {
	    if (estep[kstep] < estep[kstep - 2] && estep[kstep - 1] < estep[
		    kstep - 2]) {
		done = TRUE_;
		nvar = 0;
		i__1 = group_1.ngrp;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    for (j = 1; j <= 6; ++j) {
			++nvar;
			rbc_ref(j, i__) = rbc_ref(j, i__) - step[nvar];
		    }
		}
		rigidxyz_();
		localrgd_(minimum, grdmin);
		if (*minimum >= -big) {
		    io___83.ciunit = iounit_1.iout;
		    s_wsfe(&io___83);
		    do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(integer)
			    );
		    i__1 = kstep - 1;
		    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*minimum), (ftnlen)sizeof(
			    doublereal));
		    e_wsfe();
		} else {
		    *minimum = big;
		    io___84.ciunit = iounit_1.iout;
		    s_wsfe(&io___84);
		    do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(integer)
			    );
		    e_wsfe();
		}
	    }
	}
	if (kstep >= nstep && ! done) {
	    done = TRUE_;
	    io___85.ciunit = iounit_1.iout;
	    s_wsfe(&io___85);
	    do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	++kstep;
    }
    return 0;
} /* climbrgd_ */

#undef rbc_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine localrgd  --  PSS local search optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "localrgd" is used during the PSS local search procedure */
/*     to perform a rigid body energy minimization */


/* Subroutine */ int localrgd_(doublereal *minimum, doublereal *grdmin)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;
    static doublereal xx[75000];
    extern /* Subroutine */ int ocvm_(integer *, doublereal *, doublereal *, 
	    doublereal *, D_fp, U_fp);
    static integer nvar, oldprt;
    extern doublereal pssrgd1_(doublereal *, doublereal *);
    static logical oldverb;
    extern /* Subroutine */ int optsave_();


#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]



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
/*     ##  minima.i  --  general parameters for minimizations  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     fctmin    value below which function is deemed optimized */
/*     hguess    initial value for the H-matrix diagonal elements */
/*     maxiter   maximum number of iterations during optimization */
/*     nextiter  iteration number to use for the first iteration */




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




/*     transfer rigid body coordinates to optimization parameters */

    nvar = 0;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    ++nvar;
	    xx[nvar - 1] = rbc_ref(j, i__);
	}
    }

/*     make the call to the optimization routine */

    oldverb = inform_1.verbose;
    oldprt = inform_1.iprint;
    inform_1.verbose = FALSE_;
    inform_1.iprint = 0;
    ocvm_(&nvar, xx, minimum, grdmin, (D_fp)pssrgd1_, (U_fp)optsave_);
    inform_1.verbose = oldverb;
    inform_1.iprint = oldprt;

/*     transfer optimization parameters to rigid body coordinates */

    nvar = 0;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    ++nvar;
	    rbc_ref(j, i__) = xx[nvar - 1];
	}
    }
    return 0;
} /* localrgd_ */

#undef rbc_ref


/* Main program alias */ int pssrigid_ () { MAIN__ (); return 0; }
