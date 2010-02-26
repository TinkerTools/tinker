/* scan.f -- translated by f2c (version 20050501).
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
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

struct {
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

struct {
    doublereal dihed[1000];
    integer nomega, iomega[2000]	/* was [2][1000] */, zline[1000];
} omega_;

#define omega_1 omega_

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

static integer c__0 = 0;
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__1000 = 1000;



/*     ############################################################## */
/*     ##  COPYRIGHT (C) 1998 by Rohit Pappu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  program scan  --  maps minima on potential energy surface  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "scan" attempts to find all the local minima on a potential */
/*     energy surface via an iterative series of local searches */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Enter the Number Search Directions for L"
	    "ocal\002,\002 Search [5] :  \002,$)";
    static char fmt_30[] = "(i10)";
    static char fmt_50[] = "(/,\002 Enter the Energy Threshold for Local Min"
	    "ima\002,\002 [100.0] :  \002,$)";
    static char fmt_60[] = "(f20.0)";
    static char fmt_80[] = "(/,\002 Enter RMS Gradient per Atom Criterion"
	    "\002,\002 [0.0001] :  \002,$)";
    static char fmt_90[] = "(f20.0)";
    static char fmt_100[] = "(/,\002 Generating Seed Point for Potential Ene"
	    "rgy\002,\002 Surface Scan\002,/)";
    static char fmt_110[] = "(/,\002 Normal Mode Local Search\002,7x,\002Min"
	    "imum\002,i7,/)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2[3];
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void), mapcheck_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), 
	    localmin_(doublereal *, doublereal *), modesrch_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *);
    extern integer freeunit_(void);
    static char ext[7];
    static doublereal emap[100000];
    static integer nmap, lext, ixyz;
    extern /* Subroutine */ int final_(void);
    static doublereal range;
    static integer niter;
    static logical exist;
    static integer neigen;
    extern /* Subroutine */ int active_(void);
    static doublereal grdmin;
    static char string[120];
    extern /* Subroutine */ int getxyz_(void), makeint_(integer *), initial_(
	    void), numeral_(integer *, char *, integer *, ftnlen), nextarg_(
	    char *, logical *, ftnlen);
    static doublereal minimum;
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen), 
	    readxyz_(integer *), initrot_(void);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static icilist io___5 = { 1, string, 1, 0, 120, 1 };
    static cilist io___6 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_30, 0 };
    static icilist io___9 = { 1, string, 1, 0, 120, 1 };
    static cilist io___10 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_60, 0 };
    static icilist io___13 = { 1, string, 1, 0, 120, 1 };
    static cilist io___14 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_110, 0 };




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
/*     ##  output.i  --  control of coordinate output file format  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     archive    logical flag to save structures in an archive */
/*     noversion  logical flag governing use of filename versions */
/*     overwrite  logical flag to overwrite intermediate files inplace */
/*     cyclesave  logical flag to mark use of numbered cycle files */
/*     coordtype  selects Cartesian, internal, rigid body or none */



/*     set up the structure and mechanics calculation */

    initial_();
    getxyz_();
    mechanic_();

/*     initialize the number of minima and coordinate type */

    nmap = 0;
    s_copy(output_1.coordtype, "CARTESIAN", (ftnlen)9, (ftnlen)9);

/*     get the rotatable bonds for torsional local search */

    makeint_(&c__0);
    initrot_();
    active_();

/*     get the number of eigenvectors to use for the local search */

    neigen = -1;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___5);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&neigen, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    }
L10:
    if (neigen <= 0) {
	io___6.ciunit = iounit_1.iout;
	s_wsfe(&io___6);
	e_wsfe();
	io___7.ciunit = iounit_1.input;
	s_rsfe(&io___7);
	do_fio(&c__1, (char *)&neigen, (ftnlen)sizeof(integer));
	e_rsfe();
	if (neigen <= 0) {
	    neigen = 5;
	}
    }

/*     get the energy threshold criterion for map membership */

    range = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___9);
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&range, (ftnlen)sizeof(doublereal)
		);
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L40;
	}
    }
L40:
    if (range <= 0.) {
	io___10.ciunit = iounit_1.iout;
	s_wsfe(&io___10);
	e_wsfe();
	io___11.ciunit = iounit_1.input;
	s_rsfe(&io___11);
	do_fio(&c__1, (char *)&range, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (range <= 0.) {
	range = 100.;
    }

/*     get the termination criterion as RMS gradient per atom */

    grdmin = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___13);
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&grdmin, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L70;
	}
    }
L70:
    if (grdmin <= 0.) {
	io___14.ciunit = iounit_1.iout;
	s_wsfe(&io___14);
	e_wsfe();
	io___15.ciunit = iounit_1.input;
	s_rsfe(&io___15);
	do_fio(&c__1, (char *)&grdmin, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (grdmin <= 0.) {
	grdmin = 1e-4;
    }

/*     set the energy output precision via convergence criterion */

    if (grdmin <= 1e-6) {
	inform_1.digits = 6;
    }
    if (grdmin <= 1e-8) {
	inform_1.digits = 8;
    }

/*     find the first map point from the input structure */

    io___16.ciunit = iounit_1.iout;
    s_wsfe(&io___16);
    e_wsfe();
    localmin_(&minimum, &grdmin);
    mapcheck_(&nmap, emap, &range, &minimum, &grdmin);

/*     use normal mode local search to explore adjacent minima */

    niter = 0;
    while(niter < nmap) {
	++niter;
	io___20.ciunit = iounit_1.iout;
	s_wsfe(&io___20);
	do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	e_wsfe();
	lext = 3;
	numeral_(&niter, ext, &lext, (ftnlen)7);
	ixyz = freeunit_();
/* Writing concatenation */
	i__2[0] = files_1.leng, a__1[0] = files_1.filename;
	i__2[1] = 1, a__1[1] = ".";
	i__2[2] = lext, a__1[2] = ext;
	s_cat(xyzfile, a__1, i__2, &c__3, (ftnlen)120);
	version_(xyzfile, "old", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = ixyz;
	o__1.ofnmlen = 120;
	o__1.ofnm = xyzfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	readxyz_(&ixyz);
	cl__1.cerr = 0;
	cl__1.cunit = ixyz;
	cl__1.csta = 0;
	f_clos(&cl__1);
	modesrch_(&nmap, emap, &range, &neigen, &grdmin);
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine mapcheck  --  addition to local minimum list  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "mapcheck" checks the current minimum energy structure */
/*     for possible addition to the master list of local minima */


/* Subroutine */ int mapcheck_(integer *nmap, doublereal *emap, doublereal *
	range, doublereal *minimum, doublereal *grdmin)
{
    /* Format strings */
    static char fmt_10[] = "(/,4x,\002Potential Surface Map\002,7x,\002Minim"
	    "um\002,i7,6x,f20.8,/)";
    static char fmt_20[] = "(/,4x,\002Potential Surface Map\002,7x,\002Minim"
	    "um\002,i7,6x,f18.6,/)";
    static char fmt_30[] = "(/,4x,\002Potential Surface Map\002,7x,\002Minim"
	    "um\002,i7,6x,f16.4,/)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2[3];
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern integer freeunit_(void);
    static integer i__;
    static doublereal eps;
    static char ext[7];
    static integer lext, ixyz;
    static doublereal delta;
    static logical unique;
    extern /* Subroutine */ int prtxyz_(integer *), numeral_(integer *, char *
	    , integer *, ftnlen), version_(char *, char *, ftnlen, ftnlen);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static cilist io___29 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_30, 0 };




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




/*     check to see if the current minimum was previously found */

    /* Parameter adjustments */
    --emap;

    /* Function Body */
    eps = *grdmin;
    unique = TRUE_;
    i__1 = *nmap;
    for (i__ = 1; i__ <= i__1; ++i__) {
	delta = *minimum - emap[i__];
	if (abs(delta) < eps) {
	    unique = FALSE_;
	}
	if (delta > *range) {
	    unique = FALSE_;
	}
    }

/*     add minimum to master list if it was not previously known */

    if (unique) {
	++(*nmap);
	emap[*nmap] = *minimum;
	if (inform_1.digits >= 8) {
	    io___29.ciunit = iounit_1.iout;
	    s_wsfe(&io___29);
	    do_fio(&c__1, (char *)&(*nmap), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*minimum), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (inform_1.digits >= 6) {
	    io___30.ciunit = iounit_1.iout;
	    s_wsfe(&io___30);
	    do_fio(&c__1, (char *)&(*nmap), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*minimum), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___31.ciunit = iounit_1.iout;
	    s_wsfe(&io___31);
	    do_fio(&c__1, (char *)&(*nmap), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*minimum), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}

/*     write the coordinates of the new minimum to a file */

	lext = 3;
	numeral_(nmap, ext, &lext, (ftnlen)7);
	ixyz = freeunit_();
/* Writing concatenation */
	i__2[0] = files_1.leng, a__1[0] = files_1.filename;
	i__2[1] = 1, a__1[1] = ".";
	i__2[2] = lext, a__1[2] = ext;
	s_cat(xyzfile, a__1, i__2, &c__3, (ftnlen)120);
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
    return 0;
} /* mapcheck_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  function scan1  --  energy and gradient values for scan  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "scan1" is a service routine that computes the energy and */
/*     gradient during exploration of a potential energy surface */
/*     via iterative local search */


doublereal scan1_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static doublereal e;
    static integer i__, nvar;
    static doublereal derivs[75000]	/* was [3][25000] */;


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




/*     translate optimization parameters to atomic coordinates */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	atoms_1.x[i__ - 1] = xx[nvar];
	++nvar;
	atoms_1.y[i__ - 1] = xx[nvar];
	++nvar;
	atoms_1.z__[i__ - 1] = xx[nvar];
    }

/*     compute and store the energy and gradient */

    gradient_(&e, derivs);
    ret_val = e;

/*     store Cartesian gradient as optimization gradient */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	g[nvar] = derivs_ref(1, i__);
	++nvar;
	g[nvar] = derivs_ref(2, i__);
	++nvar;
	g[nvar] = derivs_ref(3, i__);
    }
    return ret_val;
} /* scan1_ */

#undef derivs_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine scan2  --  Hessian matrix values for scan  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "scan2" is a service routine that computes the sparse matrix */
/*     Hessian elements during exploration of a potential energy */
/*     surface via iterative local search */


/* Subroutine */ int scan2_(char *mode, doublereal *xx, doublereal *h__, 
	integer *hinit, integer *hstop, integer *hindex, doublereal *hdiag, 
	ftnlen mode_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, nvar;
    extern /* Subroutine */ int hessian_(doublereal *, integer *, integer *, 
	    integer *, doublereal *);



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




/*     translate optimization parameters to atomic coordinates */

    /* Parameter adjustments */
    --hdiag;
    --hindex;
    --hstop;
    --hinit;
    --h__;
    --xx;

    /* Function Body */
    if (s_cmp(mode, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	return 0;
    }
    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	atoms_1.x[i__ - 1] = xx[nvar];
	++nvar;
	atoms_1.y[i__ - 1] = xx[nvar];
	++nvar;
	atoms_1.z__[i__ - 1] = xx[nvar];
    }

/*     compute and store the Hessian elements */

    hessian_(&h__[1], &hinit[1], &hstop[1], &hindex[1], &hdiag[1]);
    return 0;
} /* scan2_ */



/*     ######################################################### */
/*     ##                                                     ## */
/*     ##  subroutine modesrch  --  normal mode local search  ## */
/*     ##                                                     ## */
/*     ######################################################### */


/* Subroutine */ int modesrch_(integer *nmap, doublereal *emap, doublereal *
	range, integer *neigen, doublereal *grdmin)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int mapcheck_(integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *), eigenrot_(doublereal *, 
	    doublereal *);
    static integer i__, k;
    static doublereal step[1000], eigen[1000], vects[1000000]	/* was [1000][
	    1000] */;
    extern /* Subroutine */ int makeref_(integer *), climber_(integer *, 
	    doublereal *, doublereal *, doublereal *);
    static integer nsearch;
    extern /* Subroutine */ int makeint_(integer *);
    static doublereal minimum;


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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




/*     store the current coordinates as the reference set */

    /* Parameter adjustments */
    --emap;

    /* Function Body */
    makeref_(&c__1);

/*     convert to internal coordinates and find torsional modes */

    makeint_(&c__0);
    eigenrot_(eigen, vects);

/*     search both directions along each torsional eigenvector */

    nsearch = 0;
    i__1 = *neigen;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = omega_1.nomega;
	for (k = 1; k <= i__2; ++k) {
	    step[k - 1] = vects_ref(k, omega_1.nomega - i__ + 1);
	}
	++nsearch;
	climber_(&nsearch, &minimum, step, grdmin);
	mapcheck_(nmap, &emap[1], range, &minimum, grdmin);
	i__2 = omega_1.nomega;
	for (k = 1; k <= i__2; ++k) {
	    step[k - 1] = -vects_ref(k, omega_1.nomega - i__ + 1);
	}
	++nsearch;
	climber_(&nsearch, &minimum, step, grdmin);
	mapcheck_(nmap, &emap[1], range, &minimum, grdmin);
    }
    return 0;
} /* modesrch_ */

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
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

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
    static doublereal vnorm, matrix[500500];
    extern /* Subroutine */ int hessrot_(char *, doublereal *, ftnlen);


#define hrot_ref(a_1,a_2) hrot[(a_2)*1000 + a_1 - 1001]
#define vects_ref(a_1,a_2) vects[(a_2)*1000 + a_1]



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

/*     normalize the torsional Hessian eigenvectors */

    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vnorm = 0.;
	i__2 = omega_1.nomega;
	for (j = 1; j <= i__2; ++j) {
/* Computing 2nd power */
	    d__1 = vects_ref(j, i__);
	    vnorm += d__1 * d__1;
	}
	vnorm = sqrt(vnorm);
	i__2 = omega_1.nomega;
	for (j = 1; j <= i__2; ++j) {
	    vects_ref(j, i__) = vects_ref(j, i__) / vnorm;
	}
    }
    return 0;
} /* eigenrot_ */

#undef vects_ref
#undef hrot_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine climber  --  explore single search direction  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/* Subroutine */ int climber_(integer *nsearch, doublereal *minimum, 
	doublereal *step, doublereal *grdmin)
{
    /* Format strings */
    static char fmt_10[] = "(4x,\002Search Direction\002,i4,38x,\002<<<<<"
	    "<\002)";
    static char fmt_20[] = "(4x,\002Search Direction\002,i4,38x,\002>>>>>"
	    ">\002)";
    static char fmt_30[] = "(4x,\002Search Direction\002,i4,11x,\002Step\002"
	    ",i7,6x,f20.8)";
    static char fmt_40[] = "(4x,\002Search Direction\002,i4,11x,\002Step\002"
	    ",i7,6x,f18.6)";
    static char fmt_50[] = "(4x,\002Search Direction\002,i4,11x,\002Step\002"
	    ",i7,6x,f16.4)";
    static char fmt_60[] = "(4x,\002Search Direction\002,i4,38x,\002-----"
	    "-\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int localmin_(doublereal *, doublereal *);
    static integer i__;
    static doublereal big;
    static logical done;
    static doublereal size, estep[501];
    static integer kstep, nstep;
    extern /* Subroutine */ int getref_(integer *);
    extern doublereal energy_(void);
    extern /* Subroutine */ int makeint_(integer *), makexyz_(void);

    /* Fortran I/O blocks */
    static cilist io___69 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___71 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___72 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___73 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_60, 0 };




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




/*     convert current reference coordinates to a Z-matrix */

    /* Parameter adjustments */
    --step;

    /* Function Body */
    getref_(&c__1);
    makeint_(&c__0);

/*     set the maximum number of steps and the step size */

    done = FALSE_;
    big = 1e5;
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
		localmin_(minimum, grdmin);
		if (*minimum <= -big) {
		    *minimum = big;
		    io___69.ciunit = iounit_1.iout;
		    s_wsfe(&io___69);
		    do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(integer)
			    );
		    e_wsfe();
		} else if (*minimum >= big) {
		    *minimum = big;
		    io___70.ciunit = iounit_1.iout;
		    s_wsfe(&io___70);
		    do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(integer)
			    );
		    e_wsfe();
		} else {
		    if (inform_1.digits >= 8) {
			io___71.ciunit = iounit_1.iout;
			s_wsfe(&io___71);
			do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(
				integer));
			i__1 = kstep - 1;
			do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&(*minimum), (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    } else if (inform_1.digits >= 6) {
			io___72.ciunit = iounit_1.iout;
			s_wsfe(&io___72);
			do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(
				integer));
			i__1 = kstep - 1;
			do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&(*minimum), (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    } else {
			io___73.ciunit = iounit_1.iout;
			s_wsfe(&io___73);
			do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(
				integer));
			i__1 = kstep - 1;
			do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&(*minimum), (ftnlen)sizeof(
				doublereal));
			e_wsfe();
		    }
		}
	    }
	}
	if (kstep >= nstep && ! done) {
	    done = TRUE_;
	    io___74.ciunit = iounit_1.iout;
	    s_wsfe(&io___74);
	    do_fio(&c__1, (char *)&(*nsearch), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	++kstep;
    }
    return 0;
} /* climber_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine localmin  --  optimize local search candidate  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "localmin" is used during normal mode local search to */
/*     perform a Cartesian coordinate energy minimization */


/* Subroutine */ int localmin_(doublereal *minimum, doublereal *grdmin)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static integer i__, j;
    static doublereal xx[75000], big;
    static char mode[6];
    extern /* Subroutine */ int tncg_(char *, char *, integer *, doublereal *,
	     doublereal *, doublereal *, D_fp, S_fp, U_fp, ftnlen, ftnlen);
    static integer nvar;
    static doublereal grms;
    extern doublereal scan1_(doublereal *, doublereal *);
    extern /* Subroutine */ int scan2_(char *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, ftnlen);
    static doublereal gnorm;
    static char method[6];
    static doublereal derivs[75000]	/* was [3][25000] */;
    static logical oldverb;
    extern /* Subroutine */ int optsave_();


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




/*     initialize optimization parameters and output level */

    inform_1.iwrite = 0;
    inform_1.iprint = 0;
    oldverb = inform_1.verbose;
    inform_1.verbose = FALSE_;
    big = 1e5;

/*     translate the coordinates of each atom */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	xx[nvar - 1] = atoms_1.x[i__ - 1];
	++nvar;
	xx[nvar - 1] = atoms_1.y[i__ - 1];
	++nvar;
	xx[nvar - 1] = atoms_1.z__[i__ - 1];
    }

/*     use optimization to reach the nearest local minimum */

    s_copy(mode, "AUTO", (ftnlen)6, (ftnlen)4);
    s_copy(method, "AUTO", (ftnlen)6, (ftnlen)4);
    tncg_(mode, method, &nvar, xx, minimum, grdmin, (D_fp)scan1_, (S_fp)
	    scan2_, (U_fp)optsave_, (ftnlen)6, (ftnlen)6);
/*     call lbfgs (nvar,xx,minimum,grdmin,scan1,optsave) */
    inform_1.verbose = oldverb;

/*     untranslate the final coordinates for each atom */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++nvar;
	atoms_1.x[i__ - 1] = xx[nvar - 1];
	++nvar;
	atoms_1.y[i__ - 1] = xx[nvar - 1];
	++nvar;
	atoms_1.z__[i__ - 1] = xx[nvar - 1];
    }

/*     independently check the gradient convergence criterion */

    gradient_(minimum, derivs);
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
    grms = gnorm / sqrt((doublereal) atoms_1.n);
    if (grms > *grdmin) {
	*minimum = big;
    }
    return 0;
} /* localmin_ */

#undef derivs_ref


/* Main program alias */ int scan_ () { MAIN__ (); return 0; }
