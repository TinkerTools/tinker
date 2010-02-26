/* pssrot.f -- translated by f2c (version 20050501).
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
    doublereal dihed[1000];
    integer nomega, iomega[2000]	/* was [2][1000] */, zline[1000];
} omega_;

#define omega_1 omega_

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
    doublereal m2[25000], deform, difft, diffv, diffc;
    logical use_smooth__, use_dem__, use_gda__, use_tophat__, use_stophat__;
} warp_;

#define warp_1 warp_

struct {
    doublereal zbond[25000], zang[25000], ztors[25000];
    integer iz[100000]	/* was [4][25000] */;
} zcoord_;

#define zcoord_1 zcoord_

struct {
    doublereal fctmin, hguess;
    integer maxiter, nextiter;
} minima_;

#define minima_1 minima_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static integer c__1000 = 1000;



/*     ############################################################## */
/*     ##  COPYRIGHT (C) 1997 by Rohit Pappu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  program pssrot  --  torsional potential smoothing & search  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "pssrot" implements the potential smoothing plus search method */
/*     for global optimization in torsional space */

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
    static char fmt_40[] = "(/,\002 Use Local Search to Explore the Smoothin"
	    "g Levels\002,\002 [N] :  \002,$)";
    static char fmt_50[] = "(a120)";
    static char fmt_70[] = "(/,\002 Enter the Number of Directions for Loca"
	    "l\002,\002 Search [5] :  \002,$)";
    static char fmt_80[] = "(i10)";
    static char fmt_100[] = "(/,\002 Enter the Largest Smoothing Value for L"
	    "ocal\002,\002 Search [5.0] :  \002,$)";
    static char fmt_110[] = "(f20.0)";
    static char fmt_130[] = "(/,\002 Enter RMS Gradient per Atom Criterio"
	    "n\002,\002 [0.0001] :  \002,$)";
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
    static integer i__, k;
    static logical use_local__;
    static doublereal xx[75000];
    static char ext[7];
    static doublereal rms;
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
    extern /* Subroutine */ int impose_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *), getxyz_(void);
    static doublereal deform0;
    extern /* Subroutine */ int prtxyz_(integer *), makeref_(integer *);
    extern doublereal pssrot1_();
    extern /* Subroutine */ int initial_(void), numeral_(integer *, char *, 
	    integer *, ftnlen);
    static doublereal srchmax;
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen);
    static doublereal minimum;
    extern /* Subroutine */ int moderot_(integer *, doublereal *, doublereal *
	    );
    extern /* Subroutine */ int optsave_();
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen), initrot_(void),
	     makexyz_(void);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static icilist io___5 = { 1, string, 1, 0, 120, 1 };
    static cilist io___6 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___15 = { 1, string, 1, 0, 120, 1 };
    static cilist io___16 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_80, 0 };
    static icilist io___19 = { 1, string, 1, 0, 120, 1 };
    static cilist io___20 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_110, 0 };
    static icilist io___23 = { 1, string, 1, 0, 120, 1 };
    static cilist io___24 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_150, 0 };




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




/*     set up the structure, mechanics calculation and smoothing */

    initial_();
    getxyz_();
    warp_1.use_smooth__ = TRUE_;
    warp_1.use_dem__ = TRUE_;
    mechanic_();
    initrot_();

/*     convert to Cartesian coordinates and save the initial set */

    makexyz_();
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
	    io___16.ciunit = iounit_1.iout;
	    s_wsfe(&io___16);
	    e_wsfe();
	    io___17.ciunit = iounit_1.input;
	    s_rsfe(&io___17);
	    do_fio(&c__1, (char *)&neigen, (ftnlen)sizeof(integer));
	    e_rsfe();
	    if (neigen <= 0) {
		neigen = 5;
	    }
	}
    }

/*     get the maximal smoothing level for use of local search */

    if (use_local__) {
	srchmax = -1.;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___19);
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
	    io___20.ciunit = iounit_1.iout;
	    s_wsfe(&io___20);
	    e_wsfe();
	    io___21.ciunit = iounit_1.input;
	    s_rsfe(&io___21);
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
	i__1 = s_rsli(&io___23);
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
	io___24.ciunit = iounit_1.iout;
	s_wsfe(&io___24);
	e_wsfe();
	io___25.ciunit = iounit_1.input;
	s_rsfe(&io___25);
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
/* Computing 3rd power */
	d__1 = ratio;
	warp_1.deform = deform0 * (d__1 * (d__1 * d__1));

/*     translate the initial coordinates */

	i__2 = omega_1.nomega;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xx[i__ - 1] = omega_1.dihed[i__ - 1];
	}

/*     make the call to the variable metric optimization routine */

	inform_1.iprint = 1;
	ocvm_(&omega_1.nomega, xx, &minimum, &grdmin, (D_fp)pssrot1_, (U_fp)
		optsave_);

/*     untranslate the final coordinates for each atom */

	i__2 = omega_1.nomega;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    omega_1.dihed[i__ - 1] = xx[i__ - 1];
	    zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] = omega_1.dihed[i__ - 
		    1] * 57.29577951308232088;
	}

/*     use normal mode local search to explore adjacent minima */

	if (use_local__) {
	    if (k >= npoint && warp_1.deform <= srchmax) {
		moderot_(&neigen, &minimum, &grdmin);
	    }
	}

/*     write out final energy function value and smoothing level */

	io___31.ciunit = iounit_1.iout;
	s_wsfe(&io___31);
	do_fio(&c__1, (char *)&minimum, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&warp_1.deform, (ftnlen)sizeof(doublereal));
	e_wsfe();

/*     get Cartesian coordinates and superimpose on reference */

	makexyz_();
	impose_(&atoms_1.n, refer_1.xref, refer_1.yref, refer_1.zref, &
		atoms_1.n, atoms_1.x, atoms_1.y, atoms_1.z__, &rms);

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



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function pssrot1  --  energy and gradient values for PSS  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "pssrot1" is a service routine that computes the energy and */
/*     gradient during PSS global optimization in torsional space */


doublereal pssrot1_(doublereal *xx, doublereal *g)
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
	omega_1.dihed[i__ - 1] = xx[i__];
	zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] = omega_1.dihed[i__ - 1] * 
		57.29577951308232088;
    }

/*     compute and store the energy and gradient */

    makexyz_();
    gradrot_(&e, derivs);
    ret_val = e;

/*     store torsional gradient as optimization gradient */

    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = derivs[i__ - 1];
    }
    return ret_val;
} /* pssrot1_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine moderot  --  torsional local search for PSS  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/* Subroutine */ int moderot_(integer *neigen, doublereal *minimum, 
	doublereal *grdmin)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Torsional Mode Search :\002,5x,\002Itera"
	    "tion\002,i4,6x,\002Energy\002,f12.4,/)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int eigenrot_(doublereal *, doublereal *), 
	    climbrot_(integer *, doublereal *, doublereal *, doublereal *);
    static integer i__, k;
    static doublereal eps;
    static logical done;
    static integer ndoi;
    static doublereal step[1000], eigen[1000], vects[1000000]	/* was [1000][
	    1000] */, zbest[1000], zorig[1000], minref;
    static integer nsearch;
    static doublereal minbest;
    extern /* Subroutine */ int makexyz_(void);

    /* Fortran I/O blocks */
    static cilist io___47 = { 0, 0, 0, fmt_10, 0 };



#define vects_ref(a_1,a_2) vects[(a_2)*1000 + a_1 - 1001]



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
/*     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     zbond   bond length used to define each Z-matrix atom */
/*     zang    bond angle used to define each Z-matrix atom */
/*     ztors   angle or torsion used to define Z-matrix atom */
/*     iz      defining atom numbers for each Z-matrix atom */




/*     set parameters related to the local search procedure */

    done = FALSE_;
    eps = 1e-4;
    minref = *minimum;
    minbest = *minimum;
    ndoi = 0;
    i__1 = omega_1.nomega;
    for (k = 1; k <= i__1; ++k) {
	zorig[k - 1] = zcoord_1.ztors[omega_1.zline[k - 1] - 1];
    }

/*     find local minimum along each of the steepest directions */

    while(! done) {
	++ndoi;
	io___47.ciunit = iounit_1.iout;
	s_wsfe(&io___47);
	do_fio(&c__1, (char *)&ndoi, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&minref, (ftnlen)sizeof(doublereal));
	e_wsfe();
	makexyz_();
	eigenrot_(eigen, vects);

/*     search both directions along each eigenvector in turn */

	nsearch = 0;
	i__1 = *neigen;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = omega_1.nomega;
	    for (k = 1; k <= i__2; ++k) {
		step[k - 1] = vects_ref(k, omega_1.nomega - i__ + 1);
		zcoord_1.ztors[omega_1.zline[k - 1] - 1] = zorig[k - 1];
	    }
	    ++nsearch;
	    climbrot_(&nsearch, minimum, step, grdmin);
	    if (*minimum < minbest) {
		minbest = *minimum;
		i__2 = omega_1.nomega;
		for (k = 1; k <= i__2; ++k) {
		    zbest[k - 1] = zcoord_1.ztors[omega_1.zline[k - 1] - 1];
		}
	    }
	    i__2 = omega_1.nomega;
	    for (k = 1; k <= i__2; ++k) {
		step[k - 1] = -vects_ref(k, omega_1.nomega - i__ + 1);
		zcoord_1.ztors[omega_1.zline[k - 1] - 1] = zorig[k - 1];
	    }
	    ++nsearch;
	    climbrot_(&nsearch, minimum, step, grdmin);
	    if (*minimum < minbest) {
		minbest = *minimum;
		i__2 = omega_1.nomega;
		for (k = 1; k <= i__2; ++k) {
		    zbest[k - 1] = zcoord_1.ztors[omega_1.zline[k - 1] - 1];
		}
	    }
	}

/*     check for convergence of the local search procedure */

	if (minbest < minref - eps) {
	    done = FALSE_;
	    minref = minbest;
	    i__1 = omega_1.nomega;
	    for (k = 1; k <= i__1; ++k) {
		zorig[k - 1] = zbest[k - 1];
	    }
	} else {
	    done = TRUE_;
	    *minimum = minref;
	    i__1 = omega_1.nomega;
	    for (k = 1; k <= i__1; ++k) {
		omega_1.dihed[k - 1] = zorig[k - 1] / 57.29577951308232088;
		zcoord_1.ztors[omega_1.zline[k - 1] - 1] = zorig[k - 1];
	    }
	}
    }
    return 0;
} /* moderot_ */

#undef vects_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine eigenrot  --  torsional Hessian eigenvectors  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/* Subroutine */ int eigenrot_(doublereal *eigen, doublereal *vects)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal a[1000], b[1000];
    static integer i__, j;
    static doublereal p[1000], w[1000], ta[1000], tb[1000], ty[1000], hrot[
	    1000000]	/* was [1000][1000] */;
    extern /* Subroutine */ int diagq_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer ihess;
    static doublereal matrix[500500];
    extern /* Subroutine */ int hessrot_(char *, doublereal *, ftnlen);


#define hrot_ref(a_1,a_2) hrot[(a_2)*1000 + a_1 - 1001]



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
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




/*     compute the Hessian in torsional space */

    /* Parameter adjustments */
    vects -= 1001;
    --eigen;

    /* Function Body */
    hessrot_("FULL", hrot, (ftnlen)4);

/*     place Hessian elements into triangular form */

    ihess = 0;
    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = omega_1.nomega;
	for (j = i__; j <= i__2; ++j) {
	    ++ihess;
	    matrix[ihess - 1] = hrot_ref(i__, j);
	}
    }

/*     diagonalize the Hessian to obtain eigenvalues */

    diagq_(&omega_1.nomega, &c__1000, &omega_1.nomega, matrix, &eigen[1], &
	    vects[1001], a, b, p, w, ta, tb, ty);
    return 0;
} /* eigenrot_ */

#undef hrot_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine climbrot  --  minimum from a PSS local search  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/* Subroutine */ int climbrot_(integer *nsearch, doublereal *minimum, 
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
    extern /* Subroutine */ int localrot_(doublereal *, doublereal *);
    static integer i__;
    static doublereal big;
    static logical done;
    static doublereal size, small, estep[501];
    static integer kstep, nstep;
    extern doublereal energy_(void);
    extern /* Subroutine */ int makexyz_(void);

    /* Fortran I/O blocks */
    static cilist io___74 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___75 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___76 = { 0, 0, 0, fmt_30, 0 };




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
/*     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     zbond   bond length used to define each Z-matrix atom */
/*     zang    bond angle used to define each Z-matrix atom */
/*     ztors   angle or torsion used to define Z-matrix atom */
/*     iz      defining atom numbers for each Z-matrix atom */




/*     set the maximum number of steps and the step size */

    /* Parameter adjustments */
    --step;

    /* Function Body */
    done = FALSE_;
    big = 1e10;
    small = -1e5;
    *minimum = big;
    kstep = 0;
    nstep = 65;
    size = 5.729577951308233;
    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	step[i__] = size * step[i__];
    }

/*     scan the search direction for a minimization candidate */

    while(! done) {
	if (kstep != 0) {
	    i__1 = omega_1.nomega;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] += step[i__];
	    }
	}
	makexyz_();
	estep[kstep] = energy_();
	if (kstep >= 2) {
	    if (estep[kstep] < estep[kstep - 2] && estep[kstep - 1] < estep[
		    kstep - 2]) {
		done = TRUE_;
		i__1 = omega_1.nomega;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] -= step[i__];
		}
		makexyz_();
		localrot_(minimum, grdmin);
		if (*minimum >= small) {
		    io___74.ciunit = iounit_1.iout;
		    s_wsfe(&io___74);
		    do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(integer)
			    );
		    i__1 = kstep - 1;
		    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&(*minimum), (ftnlen)sizeof(
			    doublereal));
		    e_wsfe();
		} else {
		    *minimum = big;
		    io___75.ciunit = iounit_1.iout;
		    s_wsfe(&io___75);
		    do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(integer)
			    );
		    e_wsfe();
		}
	    }
	}
	if (kstep >= nstep && ! done) {
	    done = TRUE_;
	    io___76.ciunit = iounit_1.iout;
	    s_wsfe(&io___76);
	    do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	++kstep;
    }
    return 0;
} /* climbrot_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine localrot  --  PSS local search optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "localrot" is used during the PSS local search procedure */
/*     to perform a torsional space energy minimization */


/* Subroutine */ int localrot_(doublereal *minimum, doublereal *grdmin)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal xx[75000];
    extern /* Subroutine */ int ocvm_(integer *, doublereal *, doublereal *, 
	    doublereal *, D_fp, U_fp);
    static integer oldprt;
    extern doublereal pssrot1_(doublereal *, doublereal *);
    static logical oldverb;
    extern /* Subroutine */ int optsave_();



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
/*     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     zbond   bond length used to define each Z-matrix atom */
/*     zang    bond angle used to define each Z-matrix atom */
/*     ztors   angle or torsion used to define Z-matrix atom */
/*     iz      defining atom numbers for each Z-matrix atom */




/*     translate the coordinates of each atom */

    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	omega_1.dihed[i__ - 1] = zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] / 
		57.29577951308232088;
	xx[i__ - 1] = omega_1.dihed[i__ - 1];
    }

/*     make the call to the optimization routine */

    oldverb = inform_1.verbose;
    oldprt = inform_1.iprint;
    inform_1.verbose = FALSE_;
    inform_1.iprint = 0;
    ocvm_(&omega_1.nomega, xx, minimum, grdmin, (D_fp)pssrot1_, (U_fp)
	    optsave_);
    inform_1.verbose = oldverb;
    inform_1.iprint = oldprt;

/*     untranslate the final coordinates for each atom */

    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	omega_1.dihed[i__ - 1] = xx[i__ - 1];
	zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] = omega_1.dihed[i__ - 1] * 
		57.29577951308232088;
    }
    return 0;
} /* localrot_ */

/* Main program alias */ int pssrot_ () { MAIN__ (); return 0; }
