/* embed.f -- translated by f2c (version 20050501).
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
    doublereal bnd[1000000]	/* was [1000][1000] */, vdwrad[25000], vdwmax,
	     compact, pathmax;
    logical use_invert__, use_anneal__;
} disgeo_;

#define disgeo_1 disgeo_

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
    doublereal fctmin, hguess;
    integer maxiter, nextiter;
} minima_;

#define minima_1 minima_

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
    doublereal xpfix[25000], ypfix[25000], zpfix[25000], pfix[50000]	/* 
	    was [2][25000] */, dfix[75000]	/* was [3][25000] */, afix[
	    75000]	/* was [3][25000] */, tfix[75000]	/* was [3][
	    25000] */, gfix[75000]	/* was [3][25000] */, chir[75000]	
	    /* was [3][25000] */, depth, width, rwall;
    integer npfix, ipfix[25000], kpfix[75000]	/* was [3][25000] */, ndfix, 
	    idfix[50000]	/* was [2][25000] */, nafix, iafix[75000]	
	    /* was [3][25000] */, ntfix, itfix[100000]	/* was [4][25000] */, 
	    ngfix, igfix[50000]	/* was [2][25000] */, nchir, ichir[100000]	
	    /* was [4][25000] */;
    logical use_basin__, use_wall__;
} kgeoms_;

#define kgeoms_1 kgeoms_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

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
    integer nlight, kbx[25000], kby[25000], kbz[25000], kex[25000], key[25000]
	    , kez[25000], locx[200000], locy[200000], locz[200000], rgx[
	    200000], rgy[200000], rgz[200000];
} light_;

#define light_1 light_

struct {
    doublereal ak[75000], anat[75000], afld[75000];
    integer nangle, iang[300000]	/* was [4][75000] */;
} angle_;

#define angle_1 angle_

struct {
    doublereal bk[50000], bl[50000];
    integer nbond, ibnd[100000]	/* was [2][50000] */;
} bond_;

#define bond_1 bond_

struct {
    doublereal tors1[400000]	/* was [4][100000] */, tors2[400000]	/* 
	    was [4][100000] */, tors3[400000]	/* was [4][100000] */, tors4[
	    400000]	/* was [4][100000] */, tors5[400000]	/* was [4][
	    100000] */, tors6[400000]	/* was [4][100000] */;
    integer ntors, itors[400000]	/* was [4][100000] */;
} tors_;

#define tors_1 tors_

/* Table of constant values */

static integer c__1 = 1;
static integer c__4 = 4;
static integer c__5 = 5;
static doublereal c_b128 = .25;
static integer c__1000 = 1000;
static integer c__0 = 0;
static doublereal c_b272 = 3.5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine embed  --  structures via distance geometry  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "embed" is a distance geometry routine patterned after the */
/*     ideas of Gordon Crippen, Irwin Kuntz and Tim Havel; it takes */
/*     as input a set of upper and lower bounds on the interpoint */
/*     distances, chirality restraints and torsional restraints, */
/*     and attempts to generate a set of coordinates that satisfy */
/*     the input bounds and restraints */

/*     literature references: */

/*     G. M. Crippen and T. F. Havel, "Distance Geometry and Molecular */
/*     Conformation", Research Studies Press, Letchworth U.K., 1988, */
/*     John Wiley and Sons, U.S. distributor */

/*     T. F. Havel, "An Evaluation of Computational Strategies for */
/*     Use in the Determination of Protein Structure from Distance */
/*     Constraints obtained by Nuclear Magnetic Resonance", Progress */
/*     in Biophysics and Molecular Biology, 56, 43-78 (1991) */


/* Subroutine */ int embed_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 EMBED  --  Warning, Using Metric Matri"
	    "x\002,\002 with\002,i4,\002 Negative Distances\002)";
    static char fmt_20[] = "(/,\002 EMBED  --  Warning, Using Poor Initia"
	    "l\002,\002 Coordinates\002)";
    static char fmt_30[] = "(/,\002 RMS Superposition for Original and\002"
	    ",\002 Enantiomer : \002,2f12.4)";
    static char fmt_40[] = "(/,\002 Time Required for Refinement :\002,10x,f"
	    "12.2,\002 seconds\002)";
    static char fmt_50[] = "(/,\002 Results of Distance Geometry Protocol "
	    ":\002,//,\002 Final Error Function Value :\002,10x,f16.4,//,\002"
	    " Distance Restraint Error :\002,12x,f16.4,/,\002 Hard Sphere Con"
	    "tact Error :\002,11x,f16.4,/,\002 Local Geometry Error :\002,16x"
	    ",f16.4,/,\002 Chirality-Planarity Error :\002,11x,f16.4,/,\002 T"
	    "orsional Restraint Error :\002,11x,f16.4)";
    static char fmt_60[] = "(/,\002 Radius of Gyration after Refinement  "
	    ":\002,6x,f16.4)";

    /* System generated locals */
    address a__1[4];
    integer i__1, i__2[4], i__3;
    doublereal d__1;
    olist o__1;
    cllist cl__1;
    inlist ioin__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_inqu(inlist *), f_open(olist *), f_clos(cllist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int fracdist_(char *, ftnlen), majorize_(
	    doublereal *);
    extern integer freeunit_(void);
    static integer maxinner, maxouter;
    extern /* Subroutine */ int rmserror_(char *, ftnlen);
    static doublereal a[75000];
    static integer i__, j;
    static doublereal v[75000], temp_stop__, dt, rg, temp_start__, evc[5000]	
	    /* was [1000][5] */, evl[5];
    static char ext[7];
    static integer igeo;
    static logical done;
    static integer nneg;
    static logical info;
    static doublereal secs, mass;
    static integer lext;
    extern /* Subroutine */ int eigen_(doublereal *, doublereal *, doublereal 
	    *, logical *);
    static doublereal local;
    static logical valid;
    static char title[120];
    static integer nstep;
    static logical exist;
    static doublereal chiral;
    extern /* Subroutine */ int refine_(char *, doublereal *, doublereal *, 
	    ftnlen);
    extern doublereal bnderr_(doublereal *), chirer_(doublereal *);
    extern /* Subroutine */ int getime_(doublereal *);
    static integer maxneg, ninner, nouter;
    static doublereal fctval, grdmin, bounds;
    extern doublereal vdwerr_(doublereal *), locerr_(doublereal *), torser_(
	    doublereal *);
    static doublereal derivs[3000]	/* was [3][1000] */, matrix[1000000]	
	    /* was [1000][1000] */;
    static char errtyp[7];
    extern /* Subroutine */ int dstmat_(doublereal *), metric_(doublereal *, 
	    integer *), coords_(doublereal *, doublereal *), impose_(integer *
	    , doublereal *, doublereal *, doublereal *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *), setime_(void), 
	    gyrate_(doublereal *), dmdump_(doublereal *), prtxyz_(integer *);
    static char geofile[120];
    static doublereal contact;
    extern /* Subroutine */ int chksize_(void), numeral_(integer *, char *, 
	    integer *, ftnlen);
    static doublereal rmsflip;
    extern /* Subroutine */ int explore_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, ftnlen);
    static doublereal rmsorig, torsion;

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___46 = { 0, 0, 0, fmt_60, 0 };



#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]
#define matrix_ref(a_1,a_2) matrix[(a_2)*1000 + a_1 - 1001]



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
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




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




/*     initialize any chirality restraints, then smooth the */
/*     bounds via triangle and inverse triangle inequalities; */
/*     currently these functions are performed by "distgeom" */

/*     call kchiral */
/*     if (verbose .and. n.le.130) then */
/*        title = 'Input Distance Bounds :' */
/*        call grafic (n,maxgeo,bnd,title) */
/*     end if */
/*     call geodesic */
/*     if (verbose .and. n.le.130)) then */
/*        title = 'Triangle Smoothed Bounds :' */
/*        call grafic (n,maxgeo,bnd,title) */
/*     end if */

/*     generate a distance matrix between the upper and */
/*     lower bounds, then convert to a metric matrix */

    maxinner = 3;
    maxouter = 3;
    maxneg = 2;
    nouter = 0;
    valid = FALSE_;
    while(! valid) {
	ninner = 0;
	done = FALSE_;
	while(! done) {
	    ++ninner;
	    dstmat_(matrix);
	    metric_(matrix, &nneg);
	    if (nneg <= maxneg || ninner == maxinner) {
		done = TRUE_;
	    }
	    disgeo_1.compact = 0.;
	}
	if (inform_1.verbose && nneg > maxneg) {
	    io___10.ciunit = iounit_1.iout;
	    s_wsfe(&io___10);
	    do_fio(&c__1, (char *)&nneg, (ftnlen)sizeof(integer));
	    e_wsfe();
	}

/*     find the principle components of metric matrix, then */
/*     generate the trial Cartesian coordinates */

	++nouter;
	eigen_(evl, evc, matrix, &valid);
	coords_(evl, evc);
	if (nouter == maxouter && ! valid) {
	    valid = TRUE_;
	    if (inform_1.verbose) {
		io___13.ciunit = iounit_1.iout;
		s_wsfe(&io___13);
		e_wsfe();
	    }
	}
    }

/*     superimpose embedded structure and enantiomer on reference */

    info = inform_1.verbose;
    inform_1.verbose = FALSE_;
    impose_(refer_1.nref, refer_1.xref, refer_1.yref, refer_1.zref, &
	    atoms_1.n, atoms_1.x, atoms_1.y, atoms_1.z__, &rmsorig);
    if (disgeo_1.use_invert__) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    atoms_1.x[i__ - 1] = -atoms_1.x[i__ - 1];
	}
	impose_(refer_1.nref, refer_1.xref, refer_1.yref, refer_1.zref, &
		atoms_1.n, atoms_1.x, atoms_1.y, atoms_1.z__, &rmsflip);
	if (rmsorig < rmsflip) {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		atoms_1.x[i__ - 1] = -atoms_1.x[i__ - 1];
	    }
	    impose_(refer_1.nref, refer_1.xref, refer_1.yref, refer_1.zref, &
		    atoms_1.n, atoms_1.x, atoms_1.y, atoms_1.z__, &rmsorig);
	}
	io___18.ciunit = iounit_1.iout;
	s_wsfe(&io___18);
	do_fio(&c__1, (char *)&rmsorig, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rmsflip, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    inform_1.verbose = info;

/*     compute an index of compaction for the embedded structure */

    chksize_();

/*     write out the unrefined embedded atomic coordinate set */

    if (inform_1.debug) {
	i__ = 0;
	exist = TRUE_;
	while(exist) {
	    ++i__;
	    lext = 3;
	    numeral_(&i__, ext, &lext, (ftnlen)7);
/* Writing concatenation */
	    i__2[0] = files_1.leng, a__1[0] = files_1.filename;
	    i__2[1] = 6, a__1[1] = "-embed";
	    i__2[2] = 1, a__1[2] = ".";
	    i__2[3] = lext, a__1[3] = ext;
	    s_cat(geofile, a__1, i__2, &c__4, (ftnlen)120);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = geofile;
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
	}
	igeo = freeunit_();
	o__1.oerr = 0;
	o__1.ounit = igeo;
	o__1.ofnmlen = 120;
	o__1.ofnm = geofile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	prtxyz_(&igeo);
	cl__1.cerr = 0;
	cl__1.cunit = igeo;
	cl__1.csta = 0;
	f_clos(&cl__1);
	s_copy(title, "after Embedding :", (ftnlen)120, (ftnlen)17);
	fracdist_(title, (ftnlen)120);
    }

/*     use majorization to improve initial embedded coordinates */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	matrix_ref(i__, i__) = 0.;
	i__3 = atoms_1.n;
	for (j = i__ + 1; j <= i__3; ++j) {
	    matrix_ref(j, i__) = matrix_ref(i__, j);
	}
    }
    majorize_(matrix);

/*     square the bounds for use during structure refinement */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__3 = atoms_1.n;
	for (j = 1; j <= i__3; ++j) {
/* Computing 2nd power */
	    d__1 = bnd_ref(j, i__);
	    bnd_ref(j, i__) = d__1 * d__1;
	}
    }

/*     minimize the error function via simulated annealing */

    if (inform_1.verbose) {
	setime_();
    }
    if (disgeo_1.use_anneal__) {
	inform_1.iprint = 0;
	if (inform_1.verbose) {
	    inform_1.iprint = 10;
	}
	grdmin = 1.;
	mass = 1e4;
	i__1 = atoms_1.n * 3;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    v[i__ - 1] = 0.;
	    a[i__ - 1] = 0.;
	}
	s_copy(errtyp, "FINAL", (ftnlen)7, (ftnlen)5);
	refine_(errtyp, &fctval, &grdmin, (ftnlen)7);
	nstep = 1000;
	dt = .04;
	temp_start__ = 200.;
	temp_stop__ = 200.;
	explore_(errtyp, &nstep, &dt, &mass, &temp_start__, &temp_stop__, v, 
		a, (ftnlen)7);
	nstep = 10000;
	dt = .2;
	temp_start__ = 200.;
	temp_stop__ = 0.;
	explore_(errtyp, &nstep, &dt, &mass, &temp_start__, &temp_stop__, v, 
		a, (ftnlen)7);
	grdmin = .01;
	refine_(errtyp, &fctval, &grdmin, (ftnlen)7);

/*     minimize the error function via nonlinear optimization */

    } else {
	inform_1.iprint = 0;
	if (inform_1.verbose) {
	    inform_1.iprint = 10;
	}
	grdmin = .01;
	s_copy(errtyp, "INITIAL", (ftnlen)7, (ftnlen)7);
	refine_(errtyp, &fctval, &grdmin, (ftnlen)7);
	s_copy(errtyp, "MIDDLE", (ftnlen)7, (ftnlen)6);
	refine_(errtyp, &fctval, &grdmin, (ftnlen)7);
	s_copy(errtyp, "FINAL", (ftnlen)7, (ftnlen)5);
	refine_(errtyp, &fctval, &grdmin, (ftnlen)7);
    }
    if (inform_1.verbose) {
	getime_(&secs);
	io___37.ciunit = iounit_1.iout;
	s_wsfe(&io___37);
	do_fio(&c__1, (char *)&secs, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     print the final error function and its components */

    bounds = bnderr_(derivs);
    contact = vdwerr_(derivs);
    local = locerr_(derivs);
    chiral = chirer_(derivs);
    torsion = torser_(derivs);
    io___44.ciunit = iounit_1.iout;
    s_wsfe(&io___44);
    do_fio(&c__1, (char *)&fctval, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&bounds, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&contact, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&local, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&chiral, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&torsion, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     take the root of the currently squared distance bounds */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__3 = atoms_1.n;
	for (j = 1; j <= i__3; ++j) {
	    bnd_ref(j, i__) = sqrt(bnd_ref(j, i__));
	}
    }

/*     print the final rms deviations and radius of gyration */

    s_copy(title, "after Refinement :", (ftnlen)120, (ftnlen)18);
    rmserror_(title, (ftnlen)120);
    gyrate_(&rg);
    io___46.ciunit = iounit_1.iout;
    s_wsfe(&io___46);
    do_fio(&c__1, (char *)&rg, (ftnlen)sizeof(doublereal));
    e_wsfe();
    if (inform_1.verbose && atoms_1.n <= 130) {
	dmdump_(matrix);
    }

/*     print the normalized fractional distance distribution */

    if (inform_1.debug) {
	s_copy(title, "after Refinement :", (ftnlen)120, (ftnlen)18);
	fracdist_(title, (ftnlen)120);
    }
    return 0;
} /* embed_ */

#undef matrix_ref
#undef bnd_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine kchiral  --  chirality restraint assignment  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "kchiral" determines the target value for each chirality */
/*     and planarity restraint as the signed volume of the */
/*     parallelpiped spanned by vectors from a common atom to */
/*     each of three other atoms */


/* Subroutine */ int kchiral_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Chirality and Planarity Constraints :"
	    "\002)";
    static char fmt_20[] = "(/,18x,\002Atom Numbers\002,12x,\002Signed Vol"
	    "ume\002,/)";
    static char fmt_30[] = "(i6,5x,4i6,5x,f12.4)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j;
    static doublereal c1, c2, c3;
    static integer ia, ib, ic, id;
    static doublereal xad, yad, zad, xbd, ybd, zbd, xcd, ycd, zcd;

    /* Fortran I/O blocks */
    static cilist io___64 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___65 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___66 = { 0, 0, 0, fmt_30, 0 };



#define chir_ref(a_1,a_2) kgeoms_1.chir[(a_2)*3 + a_1 - 4]
#define ichir_ref(a_1,a_2) kgeoms_1.ichir[(a_2)*4 + a_1 - 5]



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
/*     ##  kgeoms.i  --  parameters for the geometrical restraints  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xpfix      x-coordinate target for each restrained position */
/*     ypfix      y-coordinate target for each restrained position */
/*     zpfix      z-coordinate target for each restrained position */
/*     pfix       force constant and flat-well range for each position */
/*     dfix       force constant and target range for each distance */
/*     afix       force constant and target range for each angle */
/*     tfix       force constant and target range for each torsion */
/*     gfix       force constant and target range for each group distance */
/*     chir       force constant and target range for chiral centers */
/*     depth      depth of shallow Gaussian basin restraint */
/*     width      exponential width coefficient of Gaussian basin */
/*     rwall      radius of spherical droplet boundary restraint */
/*     npfix      number of position restraints to be applied */
/*     ipfix      atom number involved in each position restraint */
/*     kpfix      flags to use x-, y-, z-coordinate position restraints */
/*     ndfix      number of distance restraints to be applied */
/*     idfix      atom numbers defining each distance restraint */
/*     nafix      number of angle restraints to be applied */
/*     iafix      atom numbers defining each angle restraint */
/*     ntfix      number of torsional restraints to be applied */
/*     itfix      atom numbers defining each torsional restraint */
/*     ngfix      number of group distance restraints to be applied */
/*     igfix      group numbers defining each group distance restraint */
/*     nchir      number of chirality restraints to be applied */
/*     ichir      atom numbers defining each chirality restraint */
/*     use_basin  logical flag governing use of Gaussian basin */
/*     use_wall   logical flag governing use of droplet boundary */




/*     compute the signed volume of each parallelpiped; */
/*     if the defining atoms almost lie in a plane, then */
/*     set the signed volume to exactly zero */

    i__1 = kgeoms_1.nchir;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ichir_ref(1, i__);
	ib = ichir_ref(2, i__);
	ic = ichir_ref(3, i__);
	id = ichir_ref(4, i__);
	xad = atoms_1.x[ia - 1] - atoms_1.x[id - 1];
	yad = atoms_1.y[ia - 1] - atoms_1.y[id - 1];
	zad = atoms_1.z__[ia - 1] - atoms_1.z__[id - 1];
	xbd = atoms_1.x[ib - 1] - atoms_1.x[id - 1];
	ybd = atoms_1.y[ib - 1] - atoms_1.y[id - 1];
	zbd = atoms_1.z__[ib - 1] - atoms_1.z__[id - 1];
	xcd = atoms_1.x[ic - 1] - atoms_1.x[id - 1];
	ycd = atoms_1.y[ic - 1] - atoms_1.y[id - 1];
	zcd = atoms_1.z__[ic - 1] - atoms_1.z__[id - 1];
	c1 = ybd * zcd - zbd * ycd;
	c2 = ycd * zad - zcd * yad;
	c3 = yad * zbd - zad * ybd;
	chir_ref(1, i__) = .1;
	chir_ref(2, i__) = xad * c1 + xbd * c2 + xcd * c3;
	if ((d__1 = chir_ref(2, i__), abs(d__1)) < 1.) {
	    chir_ref(2, i__) = 0.;
	}
	chir_ref(3, i__) = chir_ref(2, i__);
    }

/*     print out the results for each restraint */

    if (inform_1.verbose) {
	if (kgeoms_1.nchir != 0) {
	    io___64.ciunit = iounit_1.iout;
	    s_wsfe(&io___64);
	    e_wsfe();
	    io___65.ciunit = iounit_1.iout;
	    s_wsfe(&io___65);
	    e_wsfe();
	}
	i__1 = kgeoms_1.nchir;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___66.ciunit = iounit_1.iout;
	    s_wsfe(&io___66);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    for (j = 1; j <= 4; ++j) {
		do_fio(&c__1, (char *)&ichir_ref(j, i__), (ftnlen)sizeof(
			integer));
	    }
	    do_fio(&c__1, (char *)&chir_ref(2, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* kchiral_ */

#undef ichir_ref
#undef chir_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine triangle  --  triangle inequality smoothing  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "triangle" smooths the upper and lower distance bounds via */
/*     the triangle inequality using a full-matrix variant of the */
/*     Floyd-Warshall shortest path algorithm; this routine is */
/*     usually much slower than the sparse matrix shortest path */
/*     methods in "geodesic" and "trifix", and should be used only */
/*     for comparison with answers generated by those routines */

/*     literature reference: */

/*     A. W. M. Dress and T. F. Havel, "Shortest-Path Problems and */
/*     Molecular Conformation", Discrete Applied Mathematics, 19, */
/*     129-144 (1988) */


/* Subroutine */ int triangle_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 TRIANGLE  --  Inconsistent Bounds;\002"
	    ",\002 Geometrically Impossible\002)";
    static char fmt_20[] = "(/,\002 Error at :\002,6x,2i6,3x,2f9.4)";
    static char fmt_30[] = "(/,\002 Traced to :\002,5x,2i6,3x,2f9.4,/,17x,2i"
	    "6,3x,2f9.4)";
    static char fmt_40[] = "(\002 TRIANGLE  --  Altered Lower Bound at\002,2"
	    "x,2i6,3x,f9.4,\002 -->\002,f9.4)";
    static char fmt_50[] = "(\002 TRIANGLE  --  Altered Upper Bound at\002,2"
	    "x,2i6,3x,f9.4,\002 -->\002,f9.4)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, ik1, ik2, jk1, jk2;
    static doublereal lij, lik, ljk, eps, uij, uik, ujk;
    extern /* Subroutine */ int fatal_(void);

    /* Fortran I/O blocks */
    static cilist io___82 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___83 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___85 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___86 = { 0, 0, 0, fmt_50, 0 };



#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]



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
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




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




/*     use full-matrix algorithm to smooth upper and lower bounds */

    eps = 1e-10;
    i__1 = atoms_1.n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = atoms_1.n - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ik1 = min(i__,k);
	    ik2 = max(i__,k);
	    lik = bnd_ref(ik2, ik1);
	    uik = bnd_ref(ik1, ik2);
	    i__3 = atoms_1.n;
	    for (j = i__ + 1; j <= i__3; ++j) {
		lij = bnd_ref(j, i__);
		uij = bnd_ref(i__, j);
		jk1 = min(j,k);
		jk2 = max(j,k);
		ljk = bnd_ref(jk2, jk1);
		ujk = bnd_ref(jk1, jk2);
/* Computing MAX */
		d__1 = lij, d__2 = lik - ujk, d__1 = max(d__1,d__2), d__2 = 
			ljk - uik;
		lij = max(d__1,d__2);
/* Computing MIN */
		d__1 = uij, d__2 = uik + ujk;
		uij = min(d__1,d__2);
		if (lij - uij > eps) {
		    io___82.ciunit = iounit_1.iout;
		    s_wsfe(&io___82);
		    e_wsfe();
		    io___83.ciunit = iounit_1.iout;
		    s_wsfe(&io___83);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&lij, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&uij, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		    io___84.ciunit = iounit_1.iout;
		    s_wsfe(&io___84);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&lik, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&uik, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ljk, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&ujk, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		    fatal_();
		}
		if (lij - bnd_ref(j, i__) > eps) {
		    io___85.ciunit = iounit_1.iout;
		    s_wsfe(&io___85);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&bnd_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&lij, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
		if (bnd_ref(i__, j) - uij > eps) {
		    io___86.ciunit = iounit_1.iout;
		    s_wsfe(&io___86);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&bnd_ref(i__, j), (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&uij, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
		bnd_ref(j, i__) = lij;
		bnd_ref(i__, j) = uij;
	    }
	}
    }
    return 0;
} /* triangle_ */

#undef bnd_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine geodesic  --  sparse matrix triangle smoothing  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "geodesic" smooths the upper and lower distance bounds via */
/*     the triangle inequality using a sparse matrix version of a */
/*     shortest path algorithm */

/*     literature reference: */

/*     G. M. Crippen and T. F. Havel, "Distance Geometry and Molecular */
/*     Conformation", Research Studies Press, Letchworth U.K., 1988, */
/*     John Wiley and Sons, U.S. distributor, see section 6-2 */


/* Subroutine */ int geodesic_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, key[50000], list[50000], stop[25000];
    extern /* Subroutine */ int sort3_(integer *, integer *, integer *);
    static integer nlist;
    static doublereal lower[25000], upper[25000];
    static integer start[25000];
    extern /* Subroutine */ int minpath_(integer *, doublereal *, doublereal *
	    , integer *, integer *, integer *);


#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]
#define idfix_ref(a_1,a_2) kgeoms_1.idfix[(a_2)*2 + a_1 - 3]



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
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  kgeoms.i  --  parameters for the geometrical restraints  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xpfix      x-coordinate target for each restrained position */
/*     ypfix      y-coordinate target for each restrained position */
/*     zpfix      z-coordinate target for each restrained position */
/*     pfix       force constant and flat-well range for each position */
/*     dfix       force constant and target range for each distance */
/*     afix       force constant and target range for each angle */
/*     tfix       force constant and target range for each torsion */
/*     gfix       force constant and target range for each group distance */
/*     chir       force constant and target range for chiral centers */
/*     depth      depth of shallow Gaussian basin restraint */
/*     width      exponential width coefficient of Gaussian basin */
/*     rwall      radius of spherical droplet boundary restraint */
/*     npfix      number of position restraints to be applied */
/*     ipfix      atom number involved in each position restraint */
/*     kpfix      flags to use x-, y-, z-coordinate position restraints */
/*     ndfix      number of distance restraints to be applied */
/*     idfix      atom numbers defining each distance restraint */
/*     nafix      number of angle restraints to be applied */
/*     iafix      atom numbers defining each angle restraint */
/*     ntfix      number of torsional restraints to be applied */
/*     itfix      atom numbers defining each torsional restraint */
/*     ngfix      number of group distance restraints to be applied */
/*     igfix      group numbers defining each group distance restraint */
/*     nchir      number of chirality restraints to be applied */
/*     ichir      atom numbers defining each chirality restraint */
/*     use_basin  logical flag governing use of Gaussian basin */
/*     use_wall   logical flag governing use of droplet boundary */




/*     build an indexed list of atoms in distance restraints */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	start[i__ - 1] = 0;
	stop[i__ - 1] = -1;
    }
    nlist = kgeoms_1.ndfix << 1;
    i__1 = kgeoms_1.ndfix;
    for (i__ = 1; i__ <= i__1; ++i__) {
	list[i__ - 1] = idfix_ref(1, i__);
	list[i__ + kgeoms_1.ndfix - 1] = idfix_ref(2, i__);
    }
    sort3_(&nlist, list, key);
    j = -1;
    i__1 = nlist;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = list[i__ - 1];
	if (k != j) {
	    start[k - 1] = i__;
	    j = k;
	}
    }
    j = -1;
    for (i__ = nlist; i__ >= 1; --i__) {
	k = list[i__ - 1];
	if (k != j) {
	    stop[k - 1] = i__;
	    j = k;
	}
    }
    i__1 = nlist;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = key[i__ - 1];
	if (k <= kgeoms_1.ndfix) {
	    list[i__ - 1] = idfix_ref(2, k);
	} else {
	    list[i__ - 1] = idfix_ref(1, k - kgeoms_1.ndfix);
	}
    }

/*     triangle smooth bounds via sparse shortest path method */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	minpath_(&i__, upper, lower, start, stop, list);
	i__2 = atoms_1.n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    bnd_ref(i__, j) = upper[j - 1];
/* Computing MAX */
	    d__1 = lower[j - 1], d__2 = bnd_ref(j, i__);
	    bnd_ref(j, i__) = max(d__1,d__2);
	}
    }
    return 0;
} /* geodesic_ */

#undef idfix_ref
#undef bnd_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine minpath  --  triangle smoothed bounds to atom  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "minpath" is a routine for finding the triangle smoothed upper */
/*     and lower bounds of each atom to a specified root atom using a */
/*     sparse variant of the Bellman-Ford shortest path algorithm */

/*     literature reference: */

/*     D. P. Bertsekas, "A Simple and Fast Label Correcting Algorithm */
/*     for Shortest Paths", Networks, 23, 703-709 (1993) */


/* Subroutine */ int minpath_(integer *root, doublereal *upper, doublereal *
	lower, integer *start, integer *stop, integer *list)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal big;
    static integer head, iarc[25000], narc, tail;
    static doublereal small;
    static logical enter;
    static integer queue[25000];
    static logical queued[25000];


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]



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
/*     ##  couple.i  --  near-neighbor atom connectivity lists   ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxn13   maximum number of atoms 1-3 connected to an atom */
/*     maxn14   maximum number of atoms 1-4 connected to an atom */
/*     maxn15   maximum number of atoms 1-5 connected to an atom */

/*     n12      number of atoms directly bonded to each atom */
/*     i12      atom numbers of atoms 1-2 connected to each atom */
/*     n13      number of atoms in a 1-3 relation to each atom */
/*     i13      atom numbers of atoms 1-3 connected to each atom */
/*     n14      number of atoms in a 1-4 relation to each atom */
/*     i14      atom numbers of atoms 1-4 connected to each atom */
/*     n15      number of atoms in a 1-5 relation to each atom */
/*     i15      atom numbers of atoms 1-5 connected to each atom */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




/*     initialize candidate atom queue and the path lengths */

    /* Parameter adjustments */
    --list;
    --stop;
    --start;
    --lower;
    --upper;

    /* Function Body */
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	queued[i__ - 1] = FALSE_;
	upper[i__] = 1e6;
	lower[i__] = 0.;
    }

/*     put the root atom into the queue of candidate atoms */

    head = *root;
    tail = *root;
    queue[*root - 1] = 0;
    queued[*root - 1] = TRUE_;
    upper[*root] = 0.;

/*     get the next candidate atom from head of queue */

    while(head != 0) {
	j = head;
	queued[j - 1] = FALSE_;
	head = queue[head - 1];

/*     make a list of arcs to the current candidate atom */

	narc = 0;
	i__1 = couple_1.n12[j - 1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = i12_ref(i__, j);
	    if (k != *root) {
		++narc;
		iarc[narc - 1] = k;
	    }
	}
	i__1 = couple_1.n13[j - 1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = i13_ref(i__, j);
	    if (k != *root) {
		++narc;
		iarc[narc - 1] = k;
	    }
	}
	i__1 = couple_1.n14[j - 1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = i14_ref(i__, j);
	    if (k != *root) {
		++narc;
		iarc[narc - 1] = k;
	    }
	}
	i__1 = stop[j];
	for (i__ = start[j]; i__ <= i__1; ++i__) {
	    k = list[i__];
	    if (k != *root) {
		++narc;
		iarc[narc - 1] = k;
	    }
	}

/*     check each arc for alteration of the path length bounds */

	i__1 = narc;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = iarc[i__ - 1];
	    if (k < j) {
		big = upper[j] + bnd_ref(k, j);
/* Computing MAX */
		d__1 = bnd_ref(j, k) - upper[j], d__2 = lower[j] - bnd_ref(k, 
			j);
		small = max(d__1,d__2);
	    } else {
		big = upper[j] + bnd_ref(j, k);
/* Computing MAX */
		d__1 = bnd_ref(k, j) - upper[j], d__2 = lower[j] - bnd_ref(j, 
			k);
		small = max(d__1,d__2);
	    }
	    enter = FALSE_;
	    if (upper[k] > big) {
		upper[k] = big;
		if (! queued[k - 1]) {
		    enter = TRUE_;
		}
	    }
	    if (lower[k] < small) {
		lower[k] = small;
		if (! queued[k - 1]) {
		    enter = TRUE_;
		}
	    }

/*     enter a new candidate atom at the tail of the queue */

	    if (enter) {
		queued[k - 1] = TRUE_;
		if (head == 0) {
		    head = k;
		    tail = k;
		    queue[k - 1] = 0;
		} else {
		    queue[tail - 1] = k;
		    queue[k - 1] = 0;
		    tail = k;
		}
	    }
	}
    }
    return 0;
} /* minpath_ */

#undef bnd_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine trifix  --  update triangle inequality bounds  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "trifix" rebuilds both the upper and lower distance bound */
/*     matrices following tightening of one or both of the bounds */
/*     between a specified pair of atoms, "p" and "q", using a */
/*     modification of Murchland's shortest path update algorithm */

/*     literature references: */

/*     P. A. Steenbrink, "Optimization of Transport Networks", John */
/*     Wiley and Sons, Bristol, 1974; see section 7.7 */

/*     R. Dionne, "Etude et Extension d'un Algorithme de Murchland", */
/*     Infor, 16, 132-146 (1978) */


/* Subroutine */ int trifix_(integer *p, integer *q)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, k, ip, iq, np, nq, pt[1000], qt[1000];
    static doublereal eps;
    static logical pun[1000], qun[1000];
    static doublereal pmin[1000], qmin[1000], pmax[1000], qmax[1000], ipmin, 
	    iqmin, ipmax, iqmax;


#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]



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
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




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




/*     initialize the set of nodes that may have changed bounds */

    eps = 1e-10;
    np = 0;
    nq = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pun[i__ - 1] = TRUE_;
	qun[i__ - 1] = TRUE_;
    }

/*     store the upper and lower bounds to "p" and "q" */

    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pmin[i__ - 1] = bnd_ref(*p, i__);
	pmax[i__ - 1] = bnd_ref(i__, *p);
    }
    i__1 = atoms_1.n;
    for (i__ = *p + 1; i__ <= i__1; ++i__) {
	pmin[i__ - 1] = bnd_ref(i__, *p);
	pmax[i__ - 1] = bnd_ref(*p, i__);
    }
    i__1 = *q;
    for (i__ = 1; i__ <= i__1; ++i__) {
	qmin[i__ - 1] = bnd_ref(*q, i__);
	qmax[i__ - 1] = bnd_ref(i__, *q);
    }
    i__1 = atoms_1.n;
    for (i__ = *q + 1; i__ <= i__1; ++i__) {
	qmin[i__ - 1] = bnd_ref(i__, *q);
	qmax[i__ - 1] = bnd_ref(*q, i__);
    }

/*     check for changes in the upper bounds to "p" and "q" */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ipmax = qmax[*p - 1] + qmax[i__ - 1];
	if (pmax[i__ - 1] > ipmax + eps) {
	    ++np;
	    pt[np - 1] = i__;
	    pmax[i__ - 1] = ipmax;
	    pun[i__ - 1] = FALSE_;
	}
	iqmax = pmax[*q - 1] + pmax[i__ - 1];
	if (qmax[i__ - 1] > iqmax + eps) {
	    ++nq;
	    qt[nq - 1] = i__;
	    qmax[i__ - 1] = iqmax;
	    qun[i__ - 1] = FALSE_;
	}
    }

/*     for node pairs whose bounds to "p" and "q" have changed, */
/*     make any needed changes to upper bound of the pair */

    i__1 = np;
    for (ip = 1; ip <= i__1; ++ip) {
	i__ = pt[ip - 1];
	ipmax = pmax[i__ - 1];
	i__2 = nq;
	for (iq = 1; iq <= i__2; ++iq) {
	    k = qt[iq - 1];
	    if (i__ < k) {
/* Computing MIN */
		d__1 = bnd_ref(i__, k), d__2 = ipmax + pmax[k - 1];
		bnd_ref(i__, k) = min(d__1,d__2);
	    } else {
/* Computing MIN */
		d__1 = bnd_ref(k, i__), d__2 = ipmax + pmax[k - 1];
		bnd_ref(k, i__) = min(d__1,d__2);
	    }
	}
    }

/*     check for changes in the lower bounds to "p" and "q" */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__1 = qmin[*p - 1] - qmax[i__ - 1], d__2 = qmin[i__ - 1] - qmax[*p - 
		1];
	ipmin = max(d__1,d__2);
	if (pmin[i__ - 1] < ipmin - eps) {
	    if (pun[i__ - 1]) {
		++np;
		pt[np - 1] = i__;
	    }
	    pmin[i__ - 1] = ipmin;
	}
/* Computing MAX */
	d__1 = pmin[*q - 1] - pmax[i__ - 1], d__2 = pmin[i__ - 1] - pmax[*q - 
		1];
	iqmin = max(d__1,d__2);
	if (qmin[i__ - 1] < iqmin - eps) {
	    if (qun[i__ - 1]) {
		++nq;
		qt[nq - 1] = i__;
	    }
	    qmin[i__ - 1] = iqmin;
	}
    }

/*     for node pairs whose bounds to "p" and "q" have changed, */
/*     make any needed changes to lower bound of the pair */

    i__1 = np;
    for (ip = 1; ip <= i__1; ++ip) {
	i__ = pt[ip - 1];
	ipmin = pmin[i__ - 1];
	ipmax = pmax[i__ - 1];
	i__2 = nq;
	for (iq = 1; iq <= i__2; ++iq) {
	    k = qt[iq - 1];
	    if (i__ < k) {
/* Computing MAX */
		d__1 = bnd_ref(k, i__), d__2 = ipmin - pmax[k - 1], d__1 = 
			max(d__1,d__2), d__2 = pmin[k - 1] - ipmax;
		bnd_ref(k, i__) = max(d__1,d__2);
	    } else {
/* Computing MAX */
		d__1 = bnd_ref(i__, k), d__2 = ipmin - pmax[k - 1], d__1 = 
			max(d__1,d__2), d__2 = pmin[k - 1] - ipmax;
		bnd_ref(i__, k) = max(d__1,d__2);
	    }
	}
    }

/*     update the upper and lower bounds to "p" and "q" */

    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	bnd_ref(*p, i__) = pmin[i__ - 1];
	bnd_ref(i__, *p) = pmax[i__ - 1];
    }
    i__1 = atoms_1.n;
    for (i__ = *p + 1; i__ <= i__1; ++i__) {
	bnd_ref(i__, *p) = pmin[i__ - 1];
	bnd_ref(*p, i__) = pmax[i__ - 1];
    }
    i__1 = *q;
    for (i__ = 1; i__ <= i__1; ++i__) {
	bnd_ref(*q, i__) = qmin[i__ - 1];
	bnd_ref(i__, *q) = qmax[i__ - 1];
    }
    i__1 = atoms_1.n;
    for (i__ = *q + 1; i__ <= i__1; ++i__) {
	bnd_ref(i__, *q) = qmin[i__ - 1];
	bnd_ref(*q, i__) = qmax[i__ - 1];
    }

/*     output the atoms updated and amount of work required */

/*     if (debug) then */
/*        write (iout,10)  p,q,np*nq */
/*  10    format (' TRIFIX  --  Bounds Update for Atoms',2i6, */
/*    &              ' with',i8,' Searches') */
/*     end if */
    return 0;
} /* trifix_ */

#undef bnd_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine grafic  --  schematic graphical matrix output  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "grafic" outputs the upper & lower triangles and diagonal */
/*     of a square matrix in a schematic form for visual inspection */


/* Subroutine */ int grafic_(integer *n, integer *np, doublereal *a, char *
	title, ftnlen title_len)
{
    /* Initialized data */

    static char ta[1] = " ";
    static char tb[1] = ".";
    static char tc[1] = "+";
    static char td[1] = "X";
    static char te[1] = "#";
    static char dash[1] = "-";
    static char digit[1*10] = "0" "1" "2" "3" "4" "5" "6" "7" "8" "9";

    /* Format strings */
    static char fmt_10[] = "(/,1x,130a1)";
    static char fmt_20[] = "(/,1x,a)";
    static char fmt_30[] = "(/,\002 Range of Above Diag Elements : \002,f13."
	    "4,\002 to\002,f13.4,/,\002 Range of Diagonal Elements :   \002,f"
	    "13.4,\002 to\002,f13.4,/,\002 Range of Below Diag Elements : "
	    "\002,f13.4,\002 to\002,f13.4)";
    static char fmt_40[] = "(/,\002 Symbol Magnitude Ordering :\002,14x,\002"
	    "# > X > + > . > ' '\002,/)";
    static char fmt_50[] = "(1x,130a1)";
    static char fmt_60[] = "(/,1x,130a1)";
    static char fmt_70[] = "()";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     i_dnnt(doublereal *);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j, k, m;
    static doublereal v, ca, cb, cc, cd, cw, cx, cy, cz, big, rcl, scl, tcl, 
	    amin, bmin, amax, dmin__, bmax, dmax__;
    static integer maxj, nrow, ndash;
    static char symbol[1*130];
    static integer minrow, maxrow;

    /* Fortran I/O blocks */
    static cilist io___138 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___140 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___149 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___150 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___168 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___169 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___170 = { 0, 0, 0, fmt_70, 0 };



#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]



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


    /* Parameter adjustments */
    a_dim1 = *np;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */


/*     set bounds of length of print row and write the header */

    minrow = 54;
    maxrow = 130;
/* Computing MIN */
    i__1 = max(*n,minrow);
    ndash = min(i__1,maxrow);
    io___138.ciunit = iounit_1.iout;
    s_wsfe(&io___138);
    i__1 = ndash;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, dash, (ftnlen)1);
    }
    e_wsfe();
    io___140.ciunit = iounit_1.iout;
    s_wsfe(&io___140);
    do_fio(&c__1, title, trimtext_(title, (ftnlen)120));
    e_wsfe();

/*     find the maximum and minimum elements of the matrix */

    big = 1e6;
    dmax__ = -big;
    dmin__ = big;
    amax = -big;
    amin = big;
    bmax = -big;
    bmin = big;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (a_ref(i__, i__) > dmax__) {
	    dmax__ = a_ref(i__, i__);
	}
	if (a_ref(i__, i__) < dmin__) {
	    dmin__ = a_ref(i__, i__);
	}
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    if (a_ref(j, i__) > amax) {
		amax = a_ref(j, i__);
	    }
	    if (a_ref(j, i__) < amin) {
		amin = a_ref(j, i__);
	    }
	    if (a_ref(i__, j) > bmax) {
		bmax = a_ref(i__, j);
	    }
	    if (a_ref(i__, j) < bmin) {
		bmin = a_ref(i__, j);
	    }
	}
    }
    io___149.ciunit = iounit_1.iout;
    s_wsfe(&io___149);
    do_fio(&c__1, (char *)&amin, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&amax, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&dmin__, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&dmax__, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&bmin, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&bmax, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     now, print out the graphical representation */

    io___150.ciunit = iounit_1.iout;
    s_wsfe(&io___150);
    e_wsfe();
    rcl = (bmax - bmin) / 5.;
    scl = (amax - amin) / 5.;
    tcl = (dmax__ - dmin__) / 9.;
    if (rcl == 0.) {
	rcl = 1.;
    }
    if (scl == 0.) {
	scl = 1.;
    }
    if (tcl == 0.) {
	tcl = 1.;
    }
    ca = amin + scl;
    cb = ca + scl;
    cc = cb + scl;
    cd = cc + scl;
    cw = bmin + rcl;
    cx = cw + rcl;
    cy = cx + rcl;
    cz = cy + rcl;
    i__1 = *n;
    i__2 = maxrow;
    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	maxj = j + maxrow - 1;
	if (maxj > *n) {
	    maxj = *n;
	}
	nrow = maxj - j + 1;
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = maxj;
	    for (k = j; k <= i__4; ++k) {
		m = k - j + 1;
		if (k < i__) {
		    v = (d__1 = a_ref(i__, k), abs(d__1));
		    if (v <= cw) {
			*(unsigned char *)&symbol[m - 1] = *(unsigned char *)&
				ta[0];
		    } else if (v <= cx) {
			*(unsigned char *)&symbol[m - 1] = *(unsigned char *)&
				tb[0];
		    } else if (v <= cy) {
			*(unsigned char *)&symbol[m - 1] = *(unsigned char *)&
				tc[0];
		    } else if (v <= cz) {
			*(unsigned char *)&symbol[m - 1] = *(unsigned char *)&
				td[0];
		    } else {
			*(unsigned char *)&symbol[m - 1] = *(unsigned char *)&
				te[0];
		    }
		} else if (k == i__) {
		    d__1 = (a_ref(i__, i__) - dmin__) / tcl;
		    *(unsigned char *)&symbol[m - 1] = *(unsigned char *)&
			    digit[i_dnnt(&d__1)];
		} else if (k > i__) {
		    v = (d__1 = a_ref(i__, k), abs(d__1));
		    if (v <= ca) {
			*(unsigned char *)&symbol[m - 1] = *(unsigned char *)&
				ta[0];
		    } else if (v <= cb) {
			*(unsigned char *)&symbol[m - 1] = *(unsigned char *)&
				tb[0];
		    } else if (v <= cc) {
			*(unsigned char *)&symbol[m - 1] = *(unsigned char *)&
				tc[0];
		    } else if (v <= cd) {
			*(unsigned char *)&symbol[m - 1] = *(unsigned char *)&
				td[0];
		    } else {
			*(unsigned char *)&symbol[m - 1] = *(unsigned char *)&
				te[0];
		    }
		}
	    }
	    io___168.ciunit = iounit_1.iout;
	    s_wsfe(&io___168);
	    i__4 = nrow;
	    for (k = 1; k <= i__4; ++k) {
		do_fio(&c__1, symbol + (k - 1), (ftnlen)1);
	    }
	    e_wsfe();
	}
	io___169.ciunit = iounit_1.iout;
	s_wsfe(&io___169);
	i__3 = ndash;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    do_fio(&c__1, dash, (ftnlen)1);
	}
	e_wsfe();
	if (maxj < *n) {
	    io___170.ciunit = iounit_1.iout;
	    s_wsfe(&io___170);
	    e_wsfe();
	}
    }
    return 0;
} /* grafic_ */

#undef a_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine dstmat  --  choose values for distance matrix  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "dstmat" selects a distance matrix containing values between */
/*     the previously smoothed upper and lower bounds; the distance */
/*     values are chosen from uniform distributions, in a triangle */
/*     correlated fashion, or using random partial metrization */


/* Subroutine */ int dstmat_(doublereal *dmx)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* Format strings */
    static char fmt_30[] = "(/,\002 Distance Matrix via Uniform Random\002"
	    ",\002 Fractions without Metrization :\002)";
    static char fmt_40[] = "(/,\002 Distance Matrix Generated via Norma"
	    "l\002,\002 Fractions without Metrization :\002)";
    static char fmt_50[] = "(/,\002 Distance Matrix Generated via Triangl"
	    "e\002,\002 Correlated Fractions :\002)";
    static char fmt_60[] = "(/,\002 Distance Matrix Generated via\002,i4,"
	    "\002-Atom\002,\002 Partial Metrization :\002)";
    static char fmt_70[] = "(/,\002 Distance Matrix Generated via Randomize"
	    "d\002,\002 Atom-Based Metrization :\002)";
    static char fmt_80[] = "(/,\002 Distance Matrix Generated via\002,i4,"
	    "\002-Atom\002,\002 Partial Metrization :\002)";
    static char fmt_90[] = "(/,\002 Distance Matrix Generated via\002,f6.2"
	    ",\002%\002,\002 Random Pairwise Metrization :\002)";
    static char fmt_100[] = "(/,\002 Distance Matrix Generated via Randomi"
	    "zed\002,\002 Pairwise Metrization :\002)";
    static char fmt_110[] = "(/,\002 Trial Distances Selected at Random fro"
	    "m\002,\002 Uniform Distribution\002)";
    static char fmt_120[] = "(/,\002 Trial Distance Beta Distribution :\002,"
	    "4x,f5.2,\002 +/-\002,f5.2,3x,\002Alpha-Beta\002,2f6.2)";
    static char fmt_130[] = "(/,\002 Average Bound Gap after Partial Metriza"
	    "tion :\002,3x,f12.4)";
    static char fmt_150[] = "(/,\002 Average Bound Gap after Partial Metriza"
	    "tion :\002,3x,f12.4)";
    static char fmt_160[] = "(/,\002 Average Bound Gap after Partial Metriza"
	    "tion :\002,3x,f12.4)";
    static char fmt_170[] = "(/,\002 Time Required for Distance Matrix :\002"
	    ",5x,f12.2,\002 seconds\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    double d_sign(doublereal *, doublereal *);
    integer i_dnnt(doublereal *);

    /* Local variables */
    static doublereal fraction;
    static integer nmetrize, i__, j, k, m;
    static doublereal gap;
    static integer mik, mjk, nik, njk;
    static doublereal eps, beta, mean, secs, corr, swap;
    static integer list[499500], next;
    extern /* Subroutine */ int sort2_(integer *, doublereal *, integer *);
    static doublereal alpha, delta, denom;
    static integer index, npair;
    static doublereal value[499500];
    static integer npart;
    static doublereal stdev;
    extern /* Subroutine */ int getime_(doublereal *);
    static char record[120];
    extern doublereal random_(void);
    static logical update;
    static char method[8], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), setime_(void), 
	    trifix_(integer *, integer *);
    extern doublereal invbeta_(doublereal *, doublereal *, doublereal *);
    static doublereal percent;
    extern /* Subroutine */ int getnumb_(char *, integer *, integer *, ftnlen)
	    ;
    static logical uniform;
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___184 = { 1, string, 1, 0, 120, 1 };
    static icilist io___185 = { 1, string, 1, 0, 120, 1 };
    static cilist io___188 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___189 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___190 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___191 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___192 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___193 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___194 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___195 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___196 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___197 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___213 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___216 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___219 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___221 = { 0, 0, 0, fmt_170, 0 };



#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]
#define dmx_ref(a_1,a_2) dmx[(a_2)*1000 + a_1]
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
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




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


    /* Parameter adjustments */
    dmx -= 1001;

    /* Function Body */


/*     initialize the method for distance element selection */

    if (first) {
	first = FALSE_;
	s_copy(method, "PAIRWISE", (ftnlen)8, (ftnlen)8);
	uniform = FALSE_;
	update = TRUE_;
	npart = 0;
	percent = 0.;
	mean = 0.;
	disgeo_1.compact = 0.;
	beta = 4.;

/*     search each line of the keyword file for options */

	i__1 = keys_1.nkey;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    next = 1;
	    s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	    gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	    upcase_(keyword, (ftnlen)20);

/*     get a distance selection method and extent of metrization */

	    if (s_cmp(keyword, "TRIAL-DISTANCE ", (ftnlen)15, (ftnlen)15) == 
		    0) {
		gettext_(record, method, &next, (ftnlen)120, (ftnlen)8);
		upcase_(method, (ftnlen)8);
		if (s_cmp(method, "HAVEL", (ftnlen)8, (ftnlen)5) == 0) {
		    getnumb_(record, &npart, &next, (ftnlen)120);
		} else if (s_cmp(method, "PARTIAL", (ftnlen)8, (ftnlen)7) == 
			0) {
		    getnumb_(record, &npart, &next, (ftnlen)120);
		} else if (s_cmp(method, "PAIRWISE", (ftnlen)8, (ftnlen)8) == 
			0) {
		    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (
			    next - 1));
		    i__2 = s_rsli(&io___184);
		    if (i__2 != 0) {
			goto L10;
		    }
		    i__2 = do_lio(&c__5, &c__1, (char *)&percent, (ftnlen)
			    sizeof(doublereal));
		    if (i__2 != 0) {
			goto L10;
		    }
		    i__2 = e_rsli();
		    if (i__2 != 0) {
			goto L10;
		    }
L10:
		    ;
		}

/*     get a choice of initial mean for the trial distribution */

	    } else if (s_cmp(keyword, "TRIAL-DISTRIBUTION ", (ftnlen)19, (
		    ftnlen)19) == 0) {
		s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next 
			- 1));
		i__2 = s_rsli(&io___185);
		if (i__2 != 0) {
		    goto L20;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&mean, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L20;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L20;
		}
L20:
		update = FALSE_;
	    }
	}

/*     set extent of partial metrization during distance selection */

	if (s_cmp(method, "HAVEL", (ftnlen)8, (ftnlen)5) == 0) {
	    if (npart <= 0 || npart >= atoms_1.n - 1) {
		npart = atoms_1.n;
	    }
	} else if (s_cmp(method, "PARTIAL", (ftnlen)8, (ftnlen)7) == 0) {
	    if (npart <= 0 || npart >= atoms_1.n - 1) {
		npart = 4;
	    }
	} else if (s_cmp(method, "PAIRWISE", (ftnlen)8, (ftnlen)8) == 0) {
	    if (percent <= 0. || percent >= 100.) {
/* Computing MIN */
		d__1 = 100., d__2 = 2e3 / (doublereal) atoms_1.n;
		percent = min(d__1,d__2);
	    }
	}

/*     set the initial distribution for selection of trial distances */

	if (s_cmp(method, "CLASSIC", (ftnlen)8, (ftnlen)7) == 0) {
	    uniform = TRUE_;
	}
	if (s_cmp(method, "TRICOR", (ftnlen)8, (ftnlen)6) == 0) {
	    uniform = TRUE_;
	}
	if (s_cmp(method, "HAVEL", (ftnlen)8, (ftnlen)5) == 0) {
	    uniform = TRUE_;
	}
	if (uniform) {
	    update = FALSE_;
	}
	if (update) {
/*           mean = 2.35d0 / log(pathmax) */
	    mean = 1.65 / pow_dd(&disgeo_1.pathmax, &c_b128);
/*           mean = 1.30d0 / (pathmax)**0.20d0 */
	}
	alpha = beta * mean / (1. - mean);
	stdev = sqrt(alpha * beta / (alpha + beta + 1.)) / (alpha + beta);
    }

/*     write out the final choice for distance matrix generation */

    if (inform_1.verbose) {
	setime_();
	if (s_cmp(method, "CLASSIC", (ftnlen)8, (ftnlen)7) == 0) {
	    io___188.ciunit = iounit_1.iout;
	    s_wsfe(&io___188);
	    e_wsfe();
	} else if (s_cmp(method, "RANDOM", (ftnlen)8, (ftnlen)6) == 0) {
	    io___189.ciunit = iounit_1.iout;
	    s_wsfe(&io___189);
	    e_wsfe();
	} else if (s_cmp(method, "TRICOR", (ftnlen)8, (ftnlen)6) == 0) {
	    io___190.ciunit = iounit_1.iout;
	    s_wsfe(&io___190);
	    e_wsfe();
	} else if (s_cmp(method, "HAVEL", (ftnlen)8, (ftnlen)5) == 0 && npart 
		< atoms_1.n) {
	    io___191.ciunit = iounit_1.iout;
	    s_wsfe(&io___191);
	    do_fio(&c__1, (char *)&npart, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else if (s_cmp(method, "HAVEL", (ftnlen)8, (ftnlen)5) == 0) {
	    io___192.ciunit = iounit_1.iout;
	    s_wsfe(&io___192);
	    e_wsfe();
	} else if (s_cmp(method, "PARTIAL", (ftnlen)8, (ftnlen)7) == 0) {
	    io___193.ciunit = iounit_1.iout;
	    s_wsfe(&io___193);
	    do_fio(&c__1, (char *)&npart, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else if (s_cmp(method, "PAIRWISE", (ftnlen)8, (ftnlen)8) == 0 && 
		percent < 100.) {
	    io___194.ciunit = iounit_1.iout;
	    s_wsfe(&io___194);
	    do_fio(&c__1, (char *)&percent, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___195.ciunit = iounit_1.iout;
	    s_wsfe(&io___195);
	    e_wsfe();
	}
    }

/*     adjust the distribution for selection of trial distances */

    if (uniform) {
	io___196.ciunit = iounit_1.iout;
	s_wsfe(&io___196);
	e_wsfe();
    } else {
	if (update) {
	    d__1 = sqrt((abs(disgeo_1.compact)));
	    alpha -= d_sign(&d__1, &disgeo_1.compact) * .2;
	    mean = alpha / (alpha + beta);
	    stdev = sqrt(alpha * beta / (alpha + beta + 1.)) / (alpha + beta);
	}
	io___197.ciunit = iounit_1.iout;
	s_wsfe(&io___197);
	do_fio(&c__1, (char *)&mean, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&stdev, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&alpha, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&beta, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     uniform or Gaussian distributed distances without metrization */

    if (s_cmp(method, "CLASSIC", (ftnlen)8, (ftnlen)7) == 0 || s_cmp(method, 
	    "RANDOM", (ftnlen)8, (ftnlen)6) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dmx_ref(i__, i__) = 0.;
	}
	i__1 = atoms_1.n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = atoms_1.n;
	    for (j = i__ + 1; j <= i__2; ++j) {
		fraction = random_();
		if (s_cmp(method, "RANDOM", (ftnlen)8, (ftnlen)6) == 0) {
		    fraction = invbeta_(&alpha, &beta, &fraction);
		}
		delta = bnd_ref(i__, j) - bnd_ref(j, i__);
		dmx_ref(j, i__) = bnd_ref(j, i__) + delta * fraction;
		dmx_ref(i__, j) = dmx_ref(j, i__);
	    }
	}

/*     Crippen's triangle correlated distance selection */

    } else if (s_cmp(method, "TRICOR", (ftnlen)8, (ftnlen)6) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dmx_ref(i__, i__) = 0.;
	}
	i__1 = atoms_1.n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = atoms_1.n;
	    for (j = i__ + 1; j <= i__2; ++j) {
		dmx_ref(j, i__) = random_();
		dmx_ref(i__, j) = dmx_ref(j, i__);
	    }
	}
	i__1 = atoms_1.n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = atoms_1.n;
	    for (j = i__ + 1; j <= i__2; ++j) {
		denom = 0.;
		dmx_ref(i__, j) = 0.;
		i__3 = atoms_1.n;
		for (k = 1; k <= i__3; ++k) {
		    if (k != i__) {
			mik = max(i__,k);
			mjk = max(j,k);
			nik = min(i__,k);
			njk = min(j,k);
			if (k == j) {
			    dmx_ref(i__, j) = dmx_ref(i__, j) + dmx_ref(j, 
				    i__);
			    denom += 1.;
			} else if (bnd_ref(njk, mjk) <= bnd_ref(nik, mik) * 
				.2) {
			    if (i__ > k) {
				corr = dmx_ref(i__, k) * .9;
			    }
			    if (k > i__) {
				corr = dmx_ref(k, i__) * .9;
			    }
			    dmx_ref(i__, j) = dmx_ref(i__, j) + corr;
			    denom += .9;
			} else if (bnd_ref(nik, mik) <= bnd_ref(njk, mjk) * 
				.2) {
			    if (j > k) {
				corr = dmx_ref(j, k) * .9;
			    }
			    if (k > j) {
				corr = dmx_ref(k, j) * .9;
			    }
			    dmx_ref(i__, j) = dmx_ref(i__, j) + corr;
			    denom += .9;
			} else if (bnd_ref(mik, nik) >= bnd_ref(njk, mjk) * 
				.9) {
			    if (j > k) {
				corr = (1. - dmx_ref(j, k)) * .5;
			    }
			    if (k > j) {
				corr = (1. - dmx_ref(k, j)) * .5;
			    }
			    dmx_ref(i__, j) = dmx_ref(i__, j) + corr;
			    denom += .5;
			} else if (bnd_ref(mjk, njk) >= bnd_ref(nik, mik) * 
				.9) {
			    if (i__ > k) {
				corr = (1. - dmx_ref(i__, k)) * .5;
			    }
			    if (k > i__) {
				corr = (1. - dmx_ref(k, i__)) * .5;
			    }
			    dmx_ref(i__, j) = dmx_ref(i__, j) + corr;
			    denom += .5;
			}
		    }
		}
		dmx_ref(i__, j) = dmx_ref(i__, j) / denom;
	    }
	}
	i__1 = atoms_1.n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = atoms_1.n;
	    for (j = i__ + 1; j <= i__2; ++j) {
		delta = bnd_ref(i__, j) - bnd_ref(j, i__);
		dmx_ref(i__, j) = bnd_ref(j, i__) + delta * dmx_ref(i__, j);
		dmx_ref(j, i__) = dmx_ref(i__, j);
	    }
	}

/*     Havel/XPLOR atom-based metrization over various distributions */

    } else if (s_cmp(method, "HAVEL", (ftnlen)8, (ftnlen)5) == 0 || s_cmp(
	    method, "PARTIAL", (ftnlen)8, (ftnlen)7) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		dmx_ref(j, i__) = bnd_ref(j, i__);
	    }
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    value[i__ - 1] = random_();
	}
	sort2_(&atoms_1.n, value, list);
	gap = 0.;
	i__1 = atoms_1.n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = list[i__ - 1];
	    i__2 = atoms_1.n;
	    for (j = i__ + 1; j <= i__2; ++j) {
		m = list[j - 1];
		fraction = random_();
		if (s_cmp(method, "PARTIAL", (ftnlen)8, (ftnlen)7) == 0) {
		    fraction = invbeta_(&alpha, &beta, &fraction);
		}
		delta = (d__1 = bnd_ref(k, m) - bnd_ref(m, k), abs(d__1));
		if (k < m) {
		    bnd_ref(k, m) = bnd_ref(m, k) + delta * fraction;
		    bnd_ref(m, k) = bnd_ref(k, m);
		} else {
		    bnd_ref(k, m) = bnd_ref(k, m) + delta * fraction;
		    bnd_ref(m, k) = bnd_ref(k, m);
		}
		if (i__ <= npart) {
		    trifix_(&k, &m);
		}
		if (i__ > npart) {
		    gap += delta;
		}
	    }
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		swap = dmx_ref(j, i__);
		dmx_ref(j, i__) = bnd_ref(j, i__);
		bnd_ref(j, i__) = swap;
	    }
	}
	if (inform_1.verbose && npart < atoms_1.n - 1) {
	    io___213.ciunit = iounit_1.iout;
	    s_wsfe(&io___213);
	    d__1 = gap / (doublereal) ((atoms_1.n - npart) * (atoms_1.n - 
		    npart - 1) / 2);
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}

/*     use partial randomized pairwise distance-based metrization */

    } else if (s_cmp(method, "PAIRWISE", (ftnlen)8, (ftnlen)8) == 0 && 
	    percent <= 10.) {
	npair = atoms_1.n * (atoms_1.n - 1) / 2;
	d__1 = percent * .01 * (doublereal) npair;
	nmetrize = i_dnnt(&d__1);
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		dmx_ref(j, i__) = bnd_ref(j, i__);
	    }
	}
	i__1 = nmetrize;
	for (i__ = 1; i__ <= i__1; ++i__) {
L140:
	    k = (integer) ((doublereal) atoms_1.n * random_()) + 1;
	    m = (integer) ((doublereal) atoms_1.n * random_()) + 1;
	    if (bnd_ref(k, m) == bnd_ref(m, k)) {
		goto L140;
	    }
	    if (k > m) {
		swap = (doublereal) k;
		k = m;
		m = (integer) swap;
	    }
	    fraction = random_();
	    fraction = invbeta_(&alpha, &beta, &fraction);
	    delta = bnd_ref(k, m) - bnd_ref(m, k);
	    bnd_ref(k, m) = bnd_ref(m, k) + delta * fraction;
	    bnd_ref(m, k) = bnd_ref(k, m);
	    trifix_(&k, &m);
	}
	gap = 0.;
	i__1 = atoms_1.n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = atoms_1.n;
	    for (j = i__; j <= i__2; ++j) {
		delta = bnd_ref(i__, j) - bnd_ref(j, i__);
		if (delta != 0.) {
		    gap += delta;
		    fraction = random_();
		    fraction = invbeta_(&alpha, &beta, &fraction);
		    bnd_ref(i__, j) = bnd_ref(j, i__) + delta * fraction;
		    bnd_ref(j, i__) = bnd_ref(i__, j);
		}
	    }
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		swap = dmx_ref(j, i__);
		dmx_ref(j, i__) = bnd_ref(j, i__);
		bnd_ref(j, i__) = swap;
	    }
	}
	if (inform_1.verbose && nmetrize < npair) {
	    io___216.ciunit = iounit_1.iout;
	    s_wsfe(&io___216);
	    d__1 = gap / (doublereal) (npair - nmetrize);
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}

/*     use randomized pairwise distance-based metrization */

    } else if (s_cmp(method, "PAIRWISE", (ftnlen)8, (ftnlen)8) == 0) {
	npair = atoms_1.n * (atoms_1.n - 1) / 2;
	d__1 = percent * .01 * (doublereal) npair;
	nmetrize = i_dnnt(&d__1);
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		dmx_ref(j, i__) = bnd_ref(j, i__);
	    }
	}
	i__1 = npair;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    value[i__ - 1] = random_();
	}
	sort2_(&npair, value, list);
	eps = 1e-10;
	gap = 0.;
	i__1 = npair;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    index = list[i__ - 1];
	    k = (integer) (((doublereal) ((atoms_1.n << 1) + 1) - sqrt((
		    doublereal) ((atoms_1.n << 2) * (atoms_1.n - 1) - (index 
		    << 3) + 9))) * .5 + eps);
	    m = atoms_1.n * (1 - k) + k * (k + 1) / 2 + index;
	    fraction = random_();
	    fraction = invbeta_(&alpha, &beta, &fraction);
	    delta = bnd_ref(k, m) - bnd_ref(m, k);
	    bnd_ref(k, m) = bnd_ref(m, k) + delta * fraction;
	    bnd_ref(m, k) = bnd_ref(k, m);
	    if (i__ <= nmetrize) {
		trifix_(&k, &m);
	    }
	    if (i__ > nmetrize) {
		gap += delta;
	    }
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		swap = dmx_ref(j, i__);
		dmx_ref(j, i__) = bnd_ref(j, i__);
		bnd_ref(j, i__) = swap;
	    }
	}
	if (inform_1.verbose && nmetrize < npair) {
	    io___219.ciunit = iounit_1.iout;
	    s_wsfe(&io___219);
	    d__1 = gap / (doublereal) (npair - nmetrize);
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     get the time required for distance matrix generation */

    if (inform_1.verbose) {
	getime_(&secs);
	io___221.ciunit = iounit_1.iout;
	s_wsfe(&io___221);
	do_fio(&c__1, (char *)&secs, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* dstmat_ */

#undef keyline_ref
#undef dmx_ref
#undef bnd_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine metric  --  computation of the metric matrix  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "metric" takes as input the trial distance matrix and computes */
/*     the metric matrix of all possible dot products between the atomic */
/*     vectors and the center of mass using the law of cosines and the */
/*     following formula for the distances to the center of mass: */

/*        dcm(i)**2 = (1/n) * sum(j=1,n)(dist(i,j)**2) */
/*                          - (1/n**2) * sum(j<k)(dist(j,k)**2) */

/*     upon output, the metric matrix is stored in the lower triangle */
/*     plus diagonal of the input trial distance matrix, the upper */
/*     triangle of the input matrix is unchanged */

/*     literature reference: */

/*     G. M. Crippen and T. F. Havel, "Stable Calculation of Coordinates */
/*     from Distance Information", Acta Cryst., A34, 282-284 (1978) */


/* Subroutine */ int metric_(doublereal *gmx, integer *nneg)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Radius of Gyration before Embedding :"
	    "\002,7x,f16.4)";
    static char fmt_20[] = "(/,\002 Atomic Distances to the Center of Mass "
	    ":\002,/)";
    static char fmt_30[] = "(6f13.4)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j;
    static doublereal rg, dcm[1000], dsq[1000], sum, total;

    /* Fortran I/O blocks */
    static cilist io___226 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___230 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___231 = { 0, 0, 0, fmt_30, 0 };



#define gmx_ref(a_1,a_2) gmx[(a_2)*1000 + a_1]



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
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




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




/*     square and sum trial distances to get radius of gyration */

    /* Parameter adjustments */
    gmx -= 1001;

    /* Function Body */
    total = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = atoms_1.n;
	for (j = i__; j <= i__2; ++j) {
/* Computing 2nd power */
	    d__1 = gmx_ref(j, i__);
	    gmx_ref(j, i__) = d__1 * d__1;
	    total += gmx_ref(j, i__);
	}
    }
/* Computing 2nd power */
    i__1 = atoms_1.n;
    total /= (doublereal) (i__1 * i__1);
    rg = sqrt(total);
    if (inform_1.verbose) {
	io___226.ciunit = iounit_1.iout;
	s_wsfe(&io___226);
	do_fio(&c__1, (char *)&rg, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     sum squared distances from each atom; the center */
/*     of mass is derived using the formula shown above */

    *nneg = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    sum += gmx_ref(i__, j);
	}
	i__2 = atoms_1.n;
	for (j = i__; j <= i__2; ++j) {
	    sum += gmx_ref(j, i__);
	}
	dsq[i__ - 1] = sum / (doublereal) atoms_1.n - total;
	dcm[i__ - 1] = sqrt((d__1 = dsq[i__ - 1], abs(d__1)));
	if (dsq[i__ - 1] < 0.) {
	    ++(*nneg);
	    dcm[i__ - 1] = -dcm[i__ - 1];
	}
    }
    if (inform_1.verbose && atoms_1.n <= 130) {
	io___230.ciunit = iounit_1.iout;
	s_wsfe(&io___230);
	e_wsfe();
	io___231.ciunit = iounit_1.iout;
	s_wsfe(&io___231);
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&dcm[i__ - 1], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }

/*     calculate the metric matrix using the law of cosines, and */
/*     place into the lower triangle of the input distance matrix */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = atoms_1.n;
	for (j = i__; j <= i__2; ++j) {
	    gmx_ref(j, i__) = (dsq[i__ - 1] + dsq[j - 1] - gmx_ref(j, i__)) * 
		    .5;
	}
    }
    return 0;
} /* metric_ */

#undef gmx_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine eigen  --  largest eigenvalues of metric metrix  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "eigen" uses the power method to compute the largest eigenvalues */
/*     and eigenvectors of the metric matrix, "valid" is set true if the */
/*     first three eigenvalues are positive */


/* Subroutine */ int eigen_(doublereal *evl, doublereal *evc, doublereal *gmx,
	 logical *valid)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Eigenvalues from Metric Matrix :\002,/)";
    static char fmt_20[] = "(5f15.4)";
    static char fmt_30[] = "(/,\002 Eigenvectors from Metric Matrix :\002,/)";
    static char fmt_40[] = "(5f15.4)";
    static char fmt_50[] = "(/,\002 Time Required for Eigenvalues :\002,9x,f"
	    "12.2,\002 seconds\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j;
    static doublereal secs, work[1000];
    static integer neigen;
    extern /* Subroutine */ int getime_(doublereal *), setime_(void), 
	    deflate_(integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___235 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___236 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___237 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___238 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___241 = { 0, 0, 0, fmt_50, 0 };



#define evc_ref(a_1,a_2) evc[(a_2)*1000 + a_1]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     initialize number of eigenvalues and convergence criteria */

    /* Parameter adjustments */
    gmx -= 1001;
    evc -= 1001;
    --evl;

    /* Function Body */
    if (inform_1.verbose) {
	setime_();
    }
    neigen = 3;

/*     compute largest eigenvalues via power method with deflation */

    deflate_(&atoms_1.n, &c__1000, &neigen, &gmx[1001], &evl[1], &evc[1001], 
	    work);

/*     check to see if the first three eigenvalues are positive */

    *valid = TRUE_;
    for (i__ = 1; i__ <= 3; ++i__) {
	if (evl[i__] < 0.) {
	    *valid = FALSE_;
	}
    }

/*     print out the eigenvalues and their eigenvectors */

    if (inform_1.verbose) {
	io___235.ciunit = iounit_1.iout;
	s_wsfe(&io___235);
	e_wsfe();
	io___236.ciunit = iounit_1.iout;
	s_wsfe(&io___236);
	i__1 = neigen;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&evl[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    if (inform_1.debug) {
	io___237.ciunit = iounit_1.iout;
	s_wsfe(&io___237);
	e_wsfe();
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___238.ciunit = iounit_1.iout;
	    s_wsfe(&io___238);
	    i__2 = neigen;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&evc_ref(i__, j), (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	}
    }

/*     get the time required for partial matrix diagonalization */

    if (inform_1.verbose) {
	getime_(&secs);
	io___241.ciunit = iounit_1.iout;
	s_wsfe(&io___241);
	do_fio(&c__1, (char *)&secs, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* eigen_ */

#undef evc_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine coords  --  converts eigenvalues to coordinates  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "coords" converts the three principal eigenvalues/vectors from */
/*     the metric matrix into atomic coordinates, and calls a routine */
/*     to compute the rms deviation from the bounds */


/* Subroutine */ int coords_(doublereal *evl, doublereal *evc)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Radius of Gyration after Embedding :\002"
	    ",8x,f16.4)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int rmserror_(char *, ftnlen);
    static integer i__, j;
    static doublereal rg;
    static char title[120];
    static integer neigen;
    extern /* Subroutine */ int gyrate_(doublereal *);

    /* Fortran I/O blocks */
    static cilist io___247 = { 0, 0, 0, fmt_10, 0 };



#define evc_ref(a_1,a_2) evc[(a_2)*1000 + a_1]



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
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




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




/*     compute coordinates from the largest eigenvalues and vectors */

    /* Parameter adjustments */
    evc -= 1001;
    --evl;

    /* Function Body */
    neigen = 3;
    i__1 = neigen;
    for (j = 1; j <= i__1; ++j) {
	evl[j] = sqrt((d__1 = evl[j], abs(d__1)));
    }
    i__1 = neigen;
    for (j = 1; j <= i__1; ++j) {
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    evc_ref(i__, j) = evl[j] * evc_ref(i__, j);
	}
    }

/*     transfer the final coordinates back to atomic vectors */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.x[i__ - 1] = evc_ref(i__, 1);
	atoms_1.y[i__ - 1] = evc_ref(i__, 2);
	atoms_1.z__[i__ - 1] = evc_ref(i__, 3);
    }

/*     find the rms bounds deviations and radius of gyration */

    if (inform_1.verbose) {
	s_copy(title, "after Embedding :", (ftnlen)120, (ftnlen)17);
	rmserror_(title, (ftnlen)120);
	gyrate_(&rg);
	io___247.ciunit = iounit_1.iout;
	s_wsfe(&io___247);
	do_fio(&c__1, (char *)&rg, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* coords_ */

#undef evc_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine chksize  --  estimate compaction or expansion  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "chksize" computes a measure of overall global structural */
/*     expansion or compaction from the number of excess upper */
/*     or lower bounds matrix violations */


/* Subroutine */ int chksize_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Index of Structure Expansion/Compactio"
	    "n :\002,7x,f12.4)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, k;
    static doublereal xi, yi, zi;
    static integer skip[25000], npair;
    static doublereal blosq;
    static integer nskip;
    static doublereal bupsq, dstsq;
    static integer nlarge, nsmall;

    /* Fortran I/O blocks */
    static cilist io___261 = { 0, 0, 0, fmt_10, 0 };



#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]



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
/*     ##  couple.i  --  near-neighbor atom connectivity lists   ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxn13   maximum number of atoms 1-3 connected to an atom */
/*     maxn14   maximum number of atoms 1-4 connected to an atom */
/*     maxn15   maximum number of atoms 1-5 connected to an atom */

/*     n12      number of atoms directly bonded to each atom */
/*     i12      atom numbers of atoms 1-2 connected to each atom */
/*     n13      number of atoms in a 1-3 relation to each atom */
/*     i13      atom numbers of atoms 1-3 connected to each atom */
/*     n14      number of atoms in a 1-4 relation to each atom */
/*     i14      atom numbers of atoms 1-4 connected to each atom */
/*     n15      number of atoms in a 1-5 relation to each atom */
/*     i15      atom numbers of atoms 1-5 connected to each atom */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




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




/*     zero out the list of atoms locally connected to each atom */

    nskip = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	skip[i__ - 1] = 0;
    }

/*     initialize counters, total pair number, and cutoff distance */

    nlarge = 0;
    nsmall = 0;
    npair = atoms_1.n * (atoms_1.n - 1) / 2;

/*     count the number of excess upper or lower bound violations */

    i__1 = atoms_1.n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
/*        do k = 1, n12(i) */
/*           skip(i12(k,i)) = i */
/*        end do */
/*        do k = 1, n13(i) */
/*           skip(i13(k,i)) = i */
/*        end do */
/*        do k = 1, n14(i) */
/*           skip(i14(k,i)) = i */
/*        end do */
	i__2 = atoms_1.n;
	for (k = i__ + 1; k <= i__2; ++k) {
	    if (skip[k - 1] == i__) {
		++nskip;
	    } else {
/* Computing 2nd power */
		d__1 = atoms_1.x[k - 1] - xi;
/* Computing 2nd power */
		d__2 = atoms_1.y[k - 1] - yi;
/* Computing 2nd power */
		d__3 = atoms_1.z__[k - 1] - zi;
		dstsq = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
/* Computing 2nd power */
		d__1 = bnd_ref(i__, k);
		bupsq = d__1 * d__1;
/* Computing 2nd power */
		d__1 = bnd_ref(k, i__);
		blosq = d__1 * d__1;
		if (dstsq > bupsq) {
		    ++nlarge;
		} else if (blosq > dstsq) {
		    ++nsmall;
		}
	    }
	}
    }

/*     set the value for the overall index of compaction */

    disgeo_1.compact = (doublereal) (nlarge - nsmall) * 100. / (doublereal) (
	    npair - nskip);
    if (inform_1.verbose) {
	io___261.ciunit = iounit_1.iout;
	s_wsfe(&io___261);
	do_fio(&c__1, (char *)&disgeo_1.compact, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* chksize_ */

#undef bnd_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine majorize  --  Guttman transform majorization  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "majorize" refines the projected coordinates by attempting to */
/*     minimize the least square residual between the trial distance */
/*     matrix and the distances computed from the coordinates */

/*     literature reference: */

/*     T. F. Havel, "An Evaluation of Computational Strategies for */
/*     Use in the Determination of Protein Structure from Distance */
/*     Constraints obtained by Nuclear Magnetic Resonance", Progress */
/*     in Biophysics and Molecular Biology, 56, 43-78 (1991) */


/* Subroutine */ int majorize_(doublereal *dmx)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Majorization to Trial Distances using"
	    "\002,\002 Constant Weights :\002,//,4x,\002Iteration\002,6x,\002"
	    "RMS Error\002,5x,\002Ave % Error\002,/)";
    static char fmt_20[] = "(5x,i5,2f16.4)";
    static char fmt_30[] = "(5x,i5,2f16.4)";
    static char fmt_40[] = "(/,\002 Radius of Gyration after Majorization "
	    ":\002,5x,f16.4)";
    static char fmt_50[] = "(/,\002 Time Required for Majorization :\002,8x,"
	    "f12.2,\002 seconds\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int rmserror_(char *, ftnlen);
    static doublereal b[25000];
    static integer i__, k;
    static doublereal rg, xi, yi, zi, xx[25000], yy[25000], zz[25000], dn1, 
	    dn2, secs;
    static integer iter;
    static doublereal dist, pairs;
    static integer niter;
    static char title[120];
    static doublereal error;
    extern /* Subroutine */ int getime_(doublereal *);
    static integer period;
    static doublereal target;
    extern /* Subroutine */ int setime_(void), gyrate_(doublereal *);
    static doublereal rmserr, average;

    /* Fortran I/O blocks */
    static cilist io___278 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___279 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___284 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___287 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___289 = { 0, 0, 0, fmt_50, 0 };



#define dmx_ref(a_1,a_2) dmx[(a_2)*1000 + a_1]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     set number of iterations and some other needed values */

    /* Parameter adjustments */
    dmx -= 1001;

    /* Function Body */
    if (inform_1.verbose) {
	setime_();
    }
    niter = 20;
    period = 5;
    pairs = (doublereal) (atoms_1.n * (atoms_1.n - 1) / 2);
    dn1 = (doublereal) (atoms_1.n - 1);
    dn2 = (doublereal) (atoms_1.n * atoms_1.n);

/*     find the average and rms error from trial distances */

    iter = 0;
    rmserr = 0.;
    average = 0.;
    i__1 = atoms_1.n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	i__2 = atoms_1.n;
	for (k = i__ + 1; k <= i__2; ++k) {
	    target = dmx_ref(k, i__);
/* Computing 2nd power */
	    d__1 = atoms_1.x[k - 1] - xi;
/* Computing 2nd power */
	    d__2 = atoms_1.y[k - 1] - yi;
/* Computing 2nd power */
	    d__3 = atoms_1.z__[k - 1] - zi;
	    dist = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	    error = dist - target;
/* Computing 2nd power */
	    d__1 = error;
	    rmserr += d__1 * d__1;
	    average += error / target;
	}
    }
    rmserr = sqrt(rmserr / pairs);
    average = average * 100. / pairs;

/*     write a header with the initial error values */

    if (inform_1.verbose) {
	io___278.ciunit = iounit_1.iout;
	s_wsfe(&io___278);
	e_wsfe();
	io___279.ciunit = iounit_1.iout;
	s_wsfe(&io___279);
	do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&rmserr, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&average, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     initialize the transformed coordinates for each atom */

    i__1 = niter;
    for (iter = 1; iter <= i__1; ++iter) {
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    b[i__ - 1] = 0.;
	    xx[i__ - 1] = 0.;
	    yy[i__ - 1] = 0.;
	    zz[i__ - 1] = 0.;

/*     form a single row of the B matrix assuming unity weights */

	    i__3 = i__ - 1;
	    for (k = 1; k <= i__3; ++k) {
/* Computing 2nd power */
		d__1 = atoms_1.x[k - 1] - xi;
/* Computing 2nd power */
		d__2 = atoms_1.y[k - 1] - yi;
/* Computing 2nd power */
		d__3 = atoms_1.z__[k - 1] - zi;
		dist = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
		b[k - 1] = -dmx_ref(k, i__) / dist;
		b[i__ - 1] -= b[k - 1];
	    }
	    i__3 = atoms_1.n;
	    for (k = i__ + 1; k <= i__3; ++k) {
/* Computing 2nd power */
		d__1 = atoms_1.x[k - 1] - xi;
/* Computing 2nd power */
		d__2 = atoms_1.y[k - 1] - yi;
/* Computing 2nd power */
		d__3 = atoms_1.z__[k - 1] - zi;
		dist = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
		b[k - 1] = -dmx_ref(k, i__) / dist;
		b[i__ - 1] -= b[k - 1];
	    }

/*     multiply the row of the B matrix by the atomic coordinates */

	    i__3 = atoms_1.n;
	    for (k = 1; k <= i__3; ++k) {
		xx[i__ - 1] += b[k - 1] * atoms_1.x[k - 1];
		yy[i__ - 1] += b[k - 1] * atoms_1.y[k - 1];
		zz[i__ - 1] += b[k - 1] * atoms_1.z__[k - 1];
	    }
	}

/*     move the intermediate values into the coordinate arrays */

	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    atoms_1.x[i__ - 1] = xx[i__ - 1];
	    atoms_1.y[i__ - 1] = yy[i__ - 1];
	    atoms_1.z__[i__ - 1] = zz[i__ - 1];
	}

/*     multiply the inverse weight matrix S+ by the coordinates */

	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xx[i__ - 1] = dn1 / dn2 * atoms_1.x[i__ - 1];
	    yy[i__ - 1] = dn1 / dn2 * atoms_1.y[i__ - 1];
	    zz[i__ - 1] = dn1 / dn2 * atoms_1.z__[i__ - 1];
	    i__3 = i__ - 1;
	    for (k = 1; k <= i__3; ++k) {
		xx[i__ - 1] -= atoms_1.x[k - 1] / dn2;
		yy[i__ - 1] -= atoms_1.y[k - 1] / dn2;
		zz[i__ - 1] -= atoms_1.z__[k - 1] / dn2;
	    }
	    i__3 = atoms_1.n;
	    for (k = i__ + 1; k <= i__3; ++k) {
		xx[i__ - 1] -= atoms_1.x[k - 1] / dn2;
		yy[i__ - 1] -= atoms_1.y[k - 1] / dn2;
		zz[i__ - 1] -= atoms_1.z__[k - 1] / dn2;
	    }
	}

/*     copy the new coordinates into their permanent arrays */

	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    atoms_1.x[i__ - 1] = xx[i__ - 1];
	    atoms_1.y[i__ - 1] = yy[i__ - 1];
	    atoms_1.z__[i__ - 1] = zz[i__ - 1];
	}

/*     find the average and rms error from trial distances */

	rmserr = 0.;
	average = 0.;
	i__2 = atoms_1.n - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    i__3 = atoms_1.n;
	    for (k = i__ + 1; k <= i__3; ++k) {
		target = dmx_ref(k, i__);
/* Computing 2nd power */
		d__1 = atoms_1.x[k - 1] - xi;
/* Computing 2nd power */
		d__2 = atoms_1.y[k - 1] - yi;
/* Computing 2nd power */
		d__3 = atoms_1.z__[k - 1] - zi;
		dist = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
		error = dist - target;
/* Computing 2nd power */
		d__1 = error;
		rmserr += d__1 * d__1;
		average += error / target;
	    }
	}
	rmserr = sqrt(rmserr / pairs);
	average = average * 100. / pairs;
	if (inform_1.verbose && iter % period == 0) {
	    io___284.ciunit = iounit_1.iout;
	    s_wsfe(&io___284);
	    do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&rmserr, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&average, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     find the rms bounds deviations and radius of gyration */

    if (inform_1.verbose) {
	s_copy(title, "after Majorization :", (ftnlen)120, (ftnlen)20);
	rmserror_(title, (ftnlen)120);
	gyrate_(&rg);
	io___287.ciunit = iounit_1.iout;
	s_wsfe(&io___287);
	do_fio(&c__1, (char *)&rg, (ftnlen)sizeof(doublereal));
	e_wsfe();

/*     get the time required for the majorization procedure */

	getime_(&secs);
	io___289.ciunit = iounit_1.iout;
	s_wsfe(&io___289);
	do_fio(&c__1, (char *)&secs, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* majorize_ */

#undef dmx_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine refine  --  minimization of initial embedding  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "refine" performs minimization of the atomic coordinates */
/*     of an initial crude embedded distance geometry structure versus */
/*     the bound, chirality, planarity and torsional error functions */


/* Subroutine */ int refine_(char *mode, doublereal *fctval, doublereal *
	grdmin, ftnlen mode_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal xx[75000];
    static integer nvar;
    extern /* Subroutine */ int lbfgs_(integer *, doublereal *, doublereal *, 
	    doublereal *, D_fp, U_fp);
    extern doublereal miderr_(), toterr_(), initerr_();
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
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




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




/*     translate the atomic coordinates to optimization variables */

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

/*     set values of parameters needed for optimization */

    s_copy(output_1.coordtype, "NONE", (ftnlen)9, (ftnlen)4);
    output_1.cyclesave = TRUE_;
/*     grdmin = 0.01d0 */
    minima_1.maxiter = nvar << 1;
    inform_1.iwrite = 0;
/*     iprint = 0 */
/*     if (verbose)  iprint = 10 */

/*     minimize initially only on the local geometry and torsions, */
/*     then on local geometry and chirality, torsions, and finally */
/*     minimize on all distance bounds, chirality and torsions */

    if (s_cmp(mode, "INITIAL", (ftnlen)7, (ftnlen)7) == 0) {
	lbfgs_(&nvar, xx, fctval, grdmin, (D_fp)initerr_, (U_fp)optsave_);
    } else if (s_cmp(mode, "MIDDLE", (ftnlen)7, (ftnlen)6) == 0) {
	lbfgs_(&nvar, xx, fctval, grdmin, (D_fp)miderr_, (U_fp)optsave_);
    } else if (s_cmp(mode, "FINAL", (ftnlen)7, (ftnlen)5) == 0) {
	lbfgs_(&nvar, xx, fctval, grdmin, (D_fp)toterr_, (U_fp)optsave_);
    }

/*     translate optimization variables back to coordinates */

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
    return 0;
} /* refine_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine explore  --  simulated annealing refinement  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "explore" uses simulated annealing on an initial crude */
/*     embedded distance geoemtry structure to refine versus the */
/*     bound, chirality, planarity and torsional error functions */


/* Subroutine */ int explore_(char *mode, integer *nstep, doublereal *dt, 
	doublereal *mass, doublereal *temp_start__, doublereal *temp_stop__, 
	doublereal *v, doublereal *a, ftnlen mode_len)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Molecular Dynamics Simulated Annealing R"
	    "efinement :\002)";
    static char fmt_20[] = "(/,\002 Steps:\002,i6,3x,\002Time/Step:\002,f6"
	    ".3,\002 ps\002,3x,\002LogMass:\002,f5.2,3x,\002Temp:\002,f6.1"
	    ",\002 to\002,f6.1)";
    static char fmt_30[] = "(/,\002 MD Step    E Total   E Potential   E Kin"
	    "etic\002,\002     Temp    MaxMove   RMS Move\002,/)";
    static char fmt_40[] = "(i6,2f13.4,f12.4,f11.2)";
    static char fmt_50[] = "(i6,2f13.4,f12.4,f11.2,2f10.4)";
    static char fmt_60[] = "(i6,2f13.4,f12.4,f11.2,2f10.4)";
    static char fmt_70[] = "(/,\002 EXPLORE  --  Simulated Annealing Unstabl"
	    "e;\002,\002 Switching to Minimization\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    double log(doublereal);
    integer do_fio(integer *, char *, ftnlen), s_cmp(char *, char *, ftnlen, 
	    ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static doublereal tau_stop__, g[75000];
    static integer i__;
    static doublereal tau_start__, xx[75000], dt2, dt_2__, dt2_2__, xbig, 
	    temp;
    static integer nvar;
    static doublereal xrms, scale, ratio, total;
    static integer istep;
    static doublereal error, prior, xmove[75000], change;
    static integer period;
    static doublereal target;
    extern doublereal miderr_(doublereal *, doublereal *), toterr_(doublereal 
	    *, doublereal *);
    static doublereal kinetic;
    extern doublereal sigmoid_(doublereal *, doublereal *), initerr_(
	    doublereal *, doublereal *);
    static doublereal tautemp;

    /* Fortran I/O blocks */
    static cilist io___303 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___304 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___311 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___312 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___320 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___322 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___323 = { 0, 0, 0, fmt_70, 0 };




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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  units.i  --  physical constants and unit conversions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     literature reference: */

/*     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended */
/*     Values of the Fundamental Physical Constants: 2006", Reviews of */
/*     Modern Physics, 80, 633-730 (2008) */

/*     The "2006 CODATA Recommended Values" are also available from */
/*     the NIST Reference on Constants, Units, and Uncertainty site */
/*     at http://physics.nist.gov/cuu/index.html */

/*     Most values below are derived from the 2006 CODATA reference */

/*     The conversion from calorie to Joule is the definition of the */
/*     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992) */

/*     The "coulomb" energy conversion factor is found by dimensional */
/*     analysis of Coulomb's Law, ie, by dividing the square of the */
/*     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is */
/*     the permittivity of vacuum (the "electric constant"); note that */
/*     eps0 is typically given in F/m, equivalent to C**2/(J-m) */

/*     The approximate value used for the Debye, 3.33564 x 10-30 C-m, */
/*     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997) */

/*     The value of "prescon" is based on definition of 1 atmosphere */
/*     as 101325 Pa set by the 10th Conference Generale des Poids et */
/*     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3 */

/*     avogadro    Avogadro's number (N) in particles/mole */
/*     lightspd    speed of light in vacuum (c) in cm/ps */
/*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K */
/*     gasconst    ideal gas constant (R) in kcal/mole/K */
/*     emass       mass of an electron in atomic mass units */
/*     joule       conversion from calories to joules */
/*     convert     conversion from kcal to g*Ang**2/ps**2 */
/*     bohr        conversion from Bohrs to Angstroms */
/*     hartree     conversion from Hartree to kcal/mole */
/*     evolt       conversion from Hartree to electron-volts */
/*     efreq       conversion from Hartree to cm-1 */
/*     coulomb     conversion from electron**2/Ang to kcal/mole */
/*     debye       conversion from electron-Ang to Debyes */
/*     prescon     conversion from kcal/mole/Ang**3 to Atm */




/*     set values of the basic simulated annealing parameters */

/*     nstep = 5000 */
/*     dt = 0.1d0 */
/*     temp_start = 200.0d0 */
/*     temp_stop = 0.0d0 */
/*     mass = 1000.0d0 */

/*     translate the atomic coordinates to annealing variables */

    /* Parameter adjustments */
    --a;
    --v;

    /* Function Body */
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

/*     initialize the velocities, accelerations and other parameters */

    dt2 = *dt * *dt;
    dt_2__ = *dt / 2.;
    dt2_2__ = dt2 / 2.;
    period = 100;
    tau_start__ = *dt * 100.;
    tau_stop__ = *dt * 10.;
    tautemp = tau_start__;

/*     print a header for the simulated annealing protocol */

    io___303.ciunit = iounit_1.iout;
    s_wsfe(&io___303);
    e_wsfe();
    io___304.ciunit = iounit_1.iout;
    s_wsfe(&io___304);
    do_fio(&c__1, (char *)&(*nstep), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*dt), (ftnlen)sizeof(doublereal));
    d__1 = log(*mass) / 2.302585092994045684;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*temp_start__), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*temp_stop__), (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     get the total error and temperature at start of dynamics */

    if (s_cmp(mode, "INITIAL", (ftnlen)7, (ftnlen)7) == 0) {
	error = initerr_(xx, g);
    } else if (s_cmp(mode, "MIDDLE", (ftnlen)7, (ftnlen)6) == 0) {
	error = miderr_(xx, g);
    } else if (s_cmp(mode, "FINAL", (ftnlen)7, (ftnlen)5) == 0) {
	error = toterr_(xx, g);
    }
    kinetic = 0.;
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = v[i__];
	kinetic += *mass * (d__1 * d__1);
    }
    kinetic = kinetic * .5 / 418.4;
    temp = kinetic * 2. / ((doublereal) nvar * .0019872066);
    total = error + kinetic;
    prior = total;
    if (inform_1.verbose) {
	io___311.ciunit = iounit_1.iout;
	s_wsfe(&io___311);
	e_wsfe();
	io___312.ciunit = iounit_1.iout;
	s_wsfe(&io___312);
	do_fio(&c__1, (char *)&c__0, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&total, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&error, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&kinetic, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&temp, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     find new positions and half-step velocities via Verlet */

    i__1 = *nstep;
    for (istep = 1; istep <= i__1; ++istep) {
	xbig = 0.;
	xrms = 0.;
	i__2 = nvar;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xmove[i__ - 1] = v[i__] * *dt + a[i__] * dt2_2__;
	    xx[i__ - 1] += xmove[i__ - 1];
	    v[i__] += a[i__] * dt_2__;
	    if ((d__1 = xmove[i__ - 1], abs(d__1)) > xbig) {
		xbig = (d__2 = xmove[i__ - 1], abs(d__2));
	    }
/* Computing 2nd power */
	    d__1 = xmove[i__ - 1];
	    xrms += d__1 * d__1;
	}
	xrms = sqrt(xrms / (doublereal) nvar);

/*     get the error function value and gradient */

	if (s_cmp(mode, "INITIAL", (ftnlen)7, (ftnlen)7) == 0) {
	    error = initerr_(xx, g);
	} else if (s_cmp(mode, "MIDDLE", (ftnlen)7, (ftnlen)6) == 0) {
	    error = miderr_(xx, g);
	} else if (s_cmp(mode, "FINAL", (ftnlen)7, (ftnlen)5) == 0) {
	    error = toterr_(xx, g);
	}

/*     use Newton's second law to get the next accelerations; */
/*     find the full-step velocities using the Verlet recursion */

	i__2 = nvar;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    a[i__] = g[i__ - 1] * -418.4 / *mass;
	    v[i__] += a[i__] * dt_2__;
	}

/*     find the total kinetic energy and system temperature */

	kinetic = 0.;
	i__2 = nvar;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    d__1 = v[i__];
	    kinetic += *mass * (d__1 * d__1);
	}
	kinetic = kinetic * .5 / 418.4;
	temp = kinetic * 2. / ((doublereal) nvar * .0019872066);
	if (temp == 0.) {
	    temp = .1;
	}

/*     set target temperature and coupling via a sigmoidal cooling */

	ratio = (doublereal) istep / (doublereal) (*nstep);
	ratio = sigmoid_(&c_b272, &ratio);
	target = *temp_start__ * (1. - ratio) + *temp_stop__ * ratio;
	tautemp = tau_start__ * (1. - ratio) + tau_stop__ * ratio;

/*     couple to external temperature bath via velocity scaling */

	scale = sqrt(*dt / tautemp * (target / temp - 1.) + 1.);
	i__2 = nvar;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    v[i__] = scale * v[i__];
	}

/*     write results for the current annealing step */

	total = error + kinetic;
	if (inform_1.verbose && istep % period == 0) {
	    io___320.ciunit = iounit_1.iout;
	    s_wsfe(&io___320);
	    do_fio(&c__1, (char *)&istep, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&total, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&error, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&kinetic, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&temp, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&xbig, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&xrms, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}

/*     check the energy change for instability in the dynamics */

	change = total - prior;
	if (change > (doublereal) atoms_1.n) {
	    i__2 = nvar;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		xx[i__ - 1] -= xmove[i__ - 1];
	    }
	    if (inform_1.verbose && istep % period != 0) {
		io___322.ciunit = iounit_1.iout;
		s_wsfe(&io___322);
		do_fio(&c__1, (char *)&istep, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&total, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&error, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&kinetic, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&temp, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&xbig, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&xrms, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    io___323.ciunit = iounit_1.iout;
	    s_wsfe(&io___323);
	    e_wsfe();
	    goto L80;
	}
    }

/*     translate annealing variables back to coordinates */

L80:
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
    return 0;
} /* explore_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine fracdist  --  fractional distance distribution  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "fracdist" computes a normalized distribution of the pairwise */
/*     fractional distances between the smoothed upper and lower bounds */

/*     literature reference: */

/*     C. M. Oshiro, J. Thomason and I. D. Kuntz, "Effects of Limited */
/*     Input Distance Constraints Upon the Distance Geometry Algorithm", */
/*     Biopolymers, 31, 1049-1064 (1991) */


/* Subroutine */ int fracdist_(char *title, ftnlen title_len)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Fractional Distance Distribution \002,a,"
	    "/)";
    static char fmt_20[] = "(8x,f8.4,8x,f8.4,8x,f8.4)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);
    integer i_dnnt(doublereal *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal fraction;
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j, k;
    static doublereal xi, yi, zi;
    static integer bin[141], sum, bin2[141], leng;
    static doublereal dist, size, range, fdist[141], fdist2[141];

    /* Fortran I/O blocks */
    static cilist io___340 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___341 = { 0, 0, 0, fmt_20, 0 };



#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]



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
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




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




/*     set the bin size and zero out the individual bins */

    size = .01;
    for (i__ = -20; i__ <= 120; ++i__) {
	bin[i__ + 20] = 0;
	bin2[i__ + 20] = 0;
    }

/*     get distribution of fractional distances between bounds */

    i__1 = atoms_1.n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	i__2 = atoms_1.n;
	for (j = i__ + 1; j <= i__2; ++j) {
/* Computing 2nd power */
	    d__1 = atoms_1.x[j - 1] - xi;
/* Computing 2nd power */
	    d__2 = atoms_1.y[j - 1] - yi;
/* Computing 2nd power */
	    d__3 = atoms_1.z__[j - 1] - zi;
	    dist = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	    range = bnd_ref(i__, j) - bnd_ref(j, i__);
	    if (range >= 1.) {
		fraction = (dist - bnd_ref(j, i__)) / range;
		d__1 = fraction / size;
		k = i_dnnt(&d__1);
/* Computing MAX */
		i__3 = -20, i__4 = min(120,k);
		k = max(i__3,i__4);
		++bin[k + 20];
		if (range >= disgeo_1.pathmax * .8) {
		    ++bin2[k + 20];
		}
	    }
	}
    }

/*     normalize the fractional distance frequency distribution */

    sum = 0;
    for (i__ = -20; i__ <= 120; ++i__) {
	sum += bin[i__ + 20];
    }
    for (i__ = -20; i__ <= 120; ++i__) {
	fdist[i__ + 20] = (doublereal) bin[i__ + 20] / (size * (doublereal) 
		sum);
    }
    sum = 0;
    for (i__ = -20; i__ <= 120; ++i__) {
	sum += bin2[i__ + 20];
    }
    for (i__ = -20; i__ <= 120; ++i__) {
	fdist2[i__ + 20] = (doublereal) bin2[i__ + 20] / (size * (doublereal) 
		sum);
    }

/*     print the normalized fractional distance probability */

    leng = trimtext_(title, (ftnlen)120);
    io___340.ciunit = iounit_1.iout;
    s_wsfe(&io___340);
    do_fio(&c__1, title, leng);
    e_wsfe();
    for (i__ = -20; i__ <= 120; ++i__) {
	io___341.ciunit = iounit_1.iout;
	s_wsfe(&io___341);
	d__1 = size * (doublereal) i__;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&fdist[i__ + 20], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&fdist2[i__ + 20], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* fracdist_ */

#undef bnd_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine rmserror  --  rms bound and restraint error  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "rmserror" computes the maximum absolute deviation and the */
/*     rms deviation from the distance bounds, and the number and */
/*     rms value of the distance restraint violations */


/* Subroutine */ int rmserror_(char *title, ftnlen title_len)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Fit to Bounds \002,a)";
    static char fmt_20[] = "(/,\002 Num Upper Bound Violations :\002,4x,i11"
	    ",\002  of \002,i12,/,\002 Num Lower Bound Violations :\002,4x,i1"
	    "1,\002  of \002,i12,/,\002 Max Upper Bound Violation :\002,4x,f1"
	    "2.4,\002  at \002,2i6,/,\002 Max Lower Bound Violation :\002,4x,"
	    "f12.4,\002  at \002,2i6,/,\002 RMS Deviation from Bounds :\002,4"
	    "x,f12.4)";
    static char fmt_30[] = "(/,\002 Num Upper Restraint Violations :\002,i"
	    "11,\002  of \002,i12,/,\002 Num Lower Restraint Violations :\002"
	    ",i11,\002  of \002,i12,/,\002 Max Upper Restraint Violation :"
	    "\002,f12.4,\002  at \002,2i6,/,\002 Max Lower Restraint Violatio"
	    "n :\002,f12.4,\002  at \002,2i6,/,\002 RMS Restraint Dist Violat"
	    "ion : \002,f12.4)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j, k, ihi, jhi, ilo, jlo;
    static doublereal rms;
    static integer leng;
    static doublereal dist;
    static integer npair;
    static doublereal himax, hierr, lomax, loerr;
    static integer nhierr, nloerr;

    /* Fortran I/O blocks */
    static cilist io___358 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___359 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___361 = { 0, 0, 0, fmt_30, 0 };



#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]
#define dfix_ref(a_1,a_2) kgeoms_1.dfix[(a_2)*3 + a_1 - 4]
#define idfix_ref(a_1,a_2) kgeoms_1.idfix[(a_2)*2 + a_1 - 3]



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
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




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
/*     ##  kgeoms.i  --  parameters for the geometrical restraints  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xpfix      x-coordinate target for each restrained position */
/*     ypfix      y-coordinate target for each restrained position */
/*     zpfix      z-coordinate target for each restrained position */
/*     pfix       force constant and flat-well range for each position */
/*     dfix       force constant and target range for each distance */
/*     afix       force constant and target range for each angle */
/*     tfix       force constant and target range for each torsion */
/*     gfix       force constant and target range for each group distance */
/*     chir       force constant and target range for chiral centers */
/*     depth      depth of shallow Gaussian basin restraint */
/*     width      exponential width coefficient of Gaussian basin */
/*     rwall      radius of spherical droplet boundary restraint */
/*     npfix      number of position restraints to be applied */
/*     ipfix      atom number involved in each position restraint */
/*     kpfix      flags to use x-, y-, z-coordinate position restraints */
/*     ndfix      number of distance restraints to be applied */
/*     idfix      atom numbers defining each distance restraint */
/*     nafix      number of angle restraints to be applied */
/*     iafix      atom numbers defining each angle restraint */
/*     ntfix      number of torsional restraints to be applied */
/*     itfix      atom numbers defining each torsional restraint */
/*     ngfix      number of group distance restraints to be applied */
/*     igfix      group numbers defining each group distance restraint */
/*     nchir      number of chirality restraints to be applied */
/*     ichir      atom numbers defining each chirality restraint */
/*     use_basin  logical flag governing use of Gaussian basin */
/*     use_wall   logical flag governing use of droplet boundary */




/*     search all atom pairs for maximal bounds deviations */

    npair = atoms_1.n * (atoms_1.n - 1) / 2;
    nloerr = 0;
    nhierr = 0;
    ilo = 0;
    jlo = 0;
    ihi = 0;
    jhi = 0;
    rms = 0.;
    lomax = 0.;
    himax = 0.;
    i__1 = atoms_1.n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = atoms_1.n;
	for (j = i__ + 1; j <= i__2; ++j) {
/* Computing 2nd power */
	    d__1 = atoms_1.x[i__ - 1] - atoms_1.x[j - 1];
/* Computing 2nd power */
	    d__2 = atoms_1.y[i__ - 1] - atoms_1.y[j - 1];
/* Computing 2nd power */
	    d__3 = atoms_1.z__[i__ - 1] - atoms_1.z__[j - 1];
	    dist = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	    dist = sqrt(dist);
	    hierr = dist - bnd_ref(i__, j);
	    if (hierr > 0.) {
		++nhierr;
/* Computing 2nd power */
		d__1 = hierr;
		rms += d__1 * d__1;
		if (hierr > himax) {
		    himax = hierr;
		    ihi = i__;
		    jhi = j;
		}
	    }
	    loerr = bnd_ref(j, i__) - dist;
	    if (loerr > 0.) {
		++nloerr;
/* Computing 2nd power */
		d__1 = loerr;
		rms += d__1 * d__1;
		if (loerr > lomax) {
		    lomax = loerr;
		    ilo = i__;
		    jlo = j;
		}
	    }
	}
    }
    rms = sqrt(rms / (doublereal) (atoms_1.n * (atoms_1.n - 1) / 2));

/*     print the maximal and rms bound deviations */

    leng = trimtext_(title, (ftnlen)120);
    io___358.ciunit = iounit_1.iout;
    s_wsfe(&io___358);
    do_fio(&c__1, title, leng);
    e_wsfe();
    io___359.ciunit = iounit_1.iout;
    s_wsfe(&io___359);
    do_fio(&c__1, (char *)&nhierr, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&npair, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nloerr, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&npair, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&himax, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ihi, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&jhi, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&lomax, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ilo, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&jlo, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&rms, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     search the list of distance restraints for violations */

    if (kgeoms_1.ndfix > 0) {
	nloerr = 0;
	nhierr = 0;
	ilo = 0;
	jlo = 0;
	ihi = 0;
	jhi = 0;
	rms = 0.;
	himax = 0.;
	lomax = 0.;
	i__1 = kgeoms_1.ndfix;
	for (k = 1; k <= i__1; ++k) {
	    i__ = idfix_ref(1, k);
	    j = idfix_ref(2, k);
/* Computing 2nd power */
	    d__1 = atoms_1.x[i__ - 1] - atoms_1.x[j - 1];
/* Computing 2nd power */
	    d__2 = atoms_1.y[i__ - 1] - atoms_1.y[j - 1];
/* Computing 2nd power */
	    d__3 = atoms_1.z__[i__ - 1] - atoms_1.z__[j - 1];
	    dist = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	    dist = sqrt(dist);
	    if (dist < dfix_ref(2, k)) {
		++nloerr;
		loerr = dfix_ref(2, k) - dist;
/* Computing 2nd power */
		d__1 = loerr;
		rms += d__1 * d__1;
		if (loerr > lomax) {
		    lomax = loerr;
		    ilo = i__;
		    jlo = j;
		}
	    } else if (dist > dfix_ref(3, k)) {
		++nhierr;
		hierr = dist - dfix_ref(3, k);
/* Computing 2nd power */
		d__1 = hierr;
		rms += d__1 * d__1;
		if (hierr > himax) {
		    himax = hierr;
		    ihi = i__;
		    jhi = j;
		}
	    }
	}
	rms = sqrt(rms / (doublereal) kgeoms_1.ndfix);

/*     print total number and rms value of restraint violations */

	io___361.ciunit = iounit_1.iout;
	s_wsfe(&io___361);
	do_fio(&c__1, (char *)&nhierr, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&kgeoms_1.ndfix, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nloerr, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&kgeoms_1.ndfix, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&himax, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ihi, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&jhi, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&lomax, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ilo, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&jlo, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&rms, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* rmserror_ */

#undef idfix_ref
#undef dfix_ref
#undef bnd_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine dmdump  --  final distance and error matrix  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "dmdump" puts the distance matrix of the final structure */
/*     into the upper half of a matrix, the distance of each atom */
/*     to the centroid on the diagonal, and the individual terms */
/*     of the bounds errors into the lower half of the matrix */


/* Subroutine */ int dmdump_(doublereal *dmd)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j;
    static doublereal sum, dist, rgsq, dist2;
    static char title[120];
    extern /* Subroutine */ int grafic_(integer *, integer *, doublereal *, 
	    char *, ftnlen);


#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]
#define dmd_ref(a_1,a_2) dmd[(a_2)*1000 + a_1]



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
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




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




/*     store the final distance matrix and bound violations */

    /* Parameter adjustments */
    dmd -= 1001;

    /* Function Body */
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dmd_ref(i__, i__) = 0.;
    }
    sum = 0.;
    i__1 = atoms_1.n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = atoms_1.n;
	for (j = i__ + 1; j <= i__2; ++j) {
/* Computing 2nd power */
	    d__1 = atoms_1.x[i__ - 1] - atoms_1.x[j - 1];
/* Computing 2nd power */
	    d__2 = atoms_1.y[i__ - 1] - atoms_1.y[j - 1];
/* Computing 2nd power */
	    d__3 = atoms_1.z__[i__ - 1] - atoms_1.z__[j - 1];
	    dist2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	    sum += dist2;
	    dmd_ref(i__, i__) = dmd_ref(i__, i__) + dist2;
	    dmd_ref(j, j) = dmd_ref(j, j) + dist2;
	    dist = sqrt(dist2);
	    dmd_ref(i__, j) = dist;
	    if (dist > bnd_ref(i__, j)) {
		dmd_ref(j, i__) = dist - bnd_ref(i__, j);
	    } else if (dist < bnd_ref(j, i__)) {
		dmd_ref(j, i__) = bnd_ref(j, i__) - dist;
	    } else {
		dmd_ref(j, i__) = 0.;
	    }
	}
    }

/*     put the distance to the centroid on the diagonal */

/* Computing 2nd power */
    i__1 = atoms_1.n;
    rgsq = sum / (doublereal) (i__1 * i__1);
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dmd_ref(i__, i__) = sqrt(dmd_ref(i__, i__) / (doublereal) atoms_1.n - 
		rgsq);
    }

/*     write out the interatomic distance and error matrices */

    s_copy(title, "Final Dist Matrix Above; DCM on Diag; Error Below :", (
	    ftnlen)120, (ftnlen)51);
    grafic_(&atoms_1.n, &c__1000, &dmd[1001], title, (ftnlen)120);
    return 0;
} /* dmdump_ */

#undef dmd_ref
#undef bnd_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  function initerr  --  initial error function and gradient  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "initerr" is the initial error function and derivatives for */
/*     a distance geometry embedding; it includes components from */
/*     the local geometry and torsional restraint errors */


doublereal initerr_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, j, nvar;
    static doublereal local;
    extern doublereal locerr_(doublereal *);
    static doublereal derivs[3000]	/* was [3][1000] */;
    extern doublereal torser_(doublereal *);
    static doublereal torsion;


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

/*     zero out the values of the atomic gradient components */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    derivs_ref(j, i__) = 0.;
	}
    }

/*     compute the local goemetry and the torsional */
/*     components of the error function and its gradient */

    local = locerr_(derivs);
    torsion = torser_(derivs);
    ret_val = local + torsion;

/*     store the atomic gradients as the optimization gradient */

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
} /* initerr_ */

#undef derivs_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  function miderr  --  second error function and gradient  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "miderr" is the secondary error function and derivatives */
/*     for a distance geometry embedding; it includes components */
/*     from the distance bounds, local geometry, chirality and */
/*     torsional restraint errors */


doublereal miderr_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, j, nvar;
    static doublereal local, chiral;
    extern doublereal bnderr_(doublereal *), chirer_(doublereal *), locerr_(
	    doublereal *);
    static doublereal bounds, derivs[3000]	/* was [3][1000] */;
    extern doublereal torser_(doublereal *);
    static doublereal torsion;


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

/*     zero out the values of the atomic gradient components */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    derivs_ref(j, i__) = 0.;
	}
    }

/*     compute the local geometry, chirality and torsional */
/*     components of the error function and its gradient */

    bounds = bnderr_(derivs);
    local = locerr_(derivs);
    chiral = chirer_(derivs);
    torsion = torser_(derivs);
    ret_val = bounds + local + chiral + torsion;

/*     store the atomic gradients as the optimization gradient */

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
} /* miderr_ */

#undef derivs_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function toterr  --  total error function and gradient  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "toterr" is the error function and derivatives for a distance */
/*     geometry embedding; it includes components from the distance */
/*     bounds, hard sphere contacts, local geometry, chirality and */
/*     torsional restraint errors */


doublereal toterr_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, j, nvar;
    static doublereal local, chiral;
    extern doublereal bnderr_(doublereal *), chirer_(doublereal *), locerr_(
	    doublereal *);
    static doublereal bounds, derivs[3000]	/* was [3][1000] */;
    extern doublereal vdwerr_(doublereal *), torser_(doublereal *);
    static doublereal contact, torsion;


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

/*     zero out the values of the atomic gradient components */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    derivs_ref(j, i__) = 0.;
	}
    }

/*     compute the distance bound, vdw, chirality and torsional */
/*     components of the error function and its gradient */

    bounds = bnderr_(derivs);
    contact = vdwerr_(derivs);
    local = locerr_(derivs);
    chiral = chirer_(derivs);
    torsion = torser_(derivs);
    ret_val = bounds + contact + local + chiral + torsion;

/*     store the atomic gradients as the optimization gradient */

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
} /* toterr_ */

#undef derivs_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function bnderr  --  computes total distance bound error  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "bnderr" is the distance bound error function and derivatives; */
/*     this version implements the original and Havel's normalized */
/*     lower bound penalty, the normalized version is preferred when */
/*     lower bounds are small (as with NMR NOE restraints), the */
/*     original penalty is needed if large lower bounds are present */


doublereal bnderr_(doublereal *derivs)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal dx, dy, dz, gx, gy, gz, gap, term, chain, scale, weigh, 
	    blosq, error, bupsq, cutsq, dstsq, buffer;


#define dfix_ref(a_1,a_2) kgeoms_1.dfix[(a_2)*3 + a_1 - 4]
#define idfix_ref(a_1,a_2) kgeoms_1.idfix[(a_2)*2 + a_1 - 3]
#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1]



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
/*     ##  kgeoms.i  --  parameters for the geometrical restraints  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xpfix      x-coordinate target for each restrained position */
/*     ypfix      y-coordinate target for each restrained position */
/*     zpfix      z-coordinate target for each restrained position */
/*     pfix       force constant and flat-well range for each position */
/*     dfix       force constant and target range for each distance */
/*     afix       force constant and target range for each angle */
/*     tfix       force constant and target range for each torsion */
/*     gfix       force constant and target range for each group distance */
/*     chir       force constant and target range for chiral centers */
/*     depth      depth of shallow Gaussian basin restraint */
/*     width      exponential width coefficient of Gaussian basin */
/*     rwall      radius of spherical droplet boundary restraint */
/*     npfix      number of position restraints to be applied */
/*     ipfix      atom number involved in each position restraint */
/*     kpfix      flags to use x-, y-, z-coordinate position restraints */
/*     ndfix      number of distance restraints to be applied */
/*     idfix      atom numbers defining each distance restraint */
/*     nafix      number of angle restraints to be applied */
/*     iafix      atom numbers defining each angle restraint */
/*     ntfix      number of torsional restraints to be applied */
/*     itfix      atom numbers defining each torsional restraint */
/*     ngfix      number of group distance restraints to be applied */
/*     igfix      group numbers defining each group distance restraint */
/*     nchir      number of chirality restraints to be applied */
/*     ichir      atom numbers defining each chirality restraint */
/*     use_basin  logical flag governing use of Gaussian basin */
/*     use_wall   logical flag governing use of droplet boundary */




/*     zero out the distance bounds error function */

    /* Parameter adjustments */
    derivs -= 4;

    /* Function Body */
    ret_val = 0.;
    scale = 10.;
    cutsq = 40.;

/*     calculate the pairwise distances between atoms */

    i__1 = kgeoms_1.ndfix;
    for (k = 1; k <= i__1; ++k) {
	i__ = idfix_ref(1, k);
	j = idfix_ref(2, k);
	dx = atoms_1.x[i__ - 1] - atoms_1.x[j - 1];
	dy = atoms_1.y[i__ - 1] - atoms_1.y[j - 1];
	dz = atoms_1.z__[i__ - 1] - atoms_1.z__[j - 1];

/*     calculate squared actual distance and bound distances; */
/*     use of a small "buffer" cleans up the final error count */

	dstsq = dx * dx + dy * dy + dz * dz;
	gap = dfix_ref(3, k) - dfix_ref(2, k);
	buffer = min(1.,gap) * .05;
/* Computing 2nd power */
	d__1 = dfix_ref(2, k) + buffer;
	blosq = d__1 * d__1;
/* Computing 2nd power */
	d__1 = dfix_ref(3, k) - buffer;
	bupsq = d__1 * d__1;

/*     error and derivatives for upper bound violation */

	if (dstsq > bupsq) {
	    weigh = scale * dfix_ref(1, k);
	    term = (dstsq - bupsq) / bupsq;
	    chain = weigh * 4. * term / bupsq;
/* Computing 2nd power */
	    d__1 = term;
	    error = weigh * (d__1 * d__1);
	    gx = dx * chain;
	    gy = dy * chain;
	    gz = dz * chain;
	    ret_val += error;
	    derivs_ref(1, i__) = derivs_ref(1, i__) + gx;
	    derivs_ref(2, i__) = derivs_ref(2, i__) + gy;
	    derivs_ref(3, i__) = derivs_ref(3, i__) + gz;
	    derivs_ref(1, j) = derivs_ref(1, j) - gx;
	    derivs_ref(2, j) = derivs_ref(2, j) - gy;
	    derivs_ref(3, j) = derivs_ref(3, j) - gz;

/*     error and derivatives for lower bound violation */

	} else if (dstsq < blosq) {
	    weigh = scale * dfix_ref(1, k);
	    if (blosq > cutsq) {
		term = (blosq - dstsq) / dstsq;
/* Computing 2nd power */
		d__1 = dstsq;
		chain = weigh * -4. * term * (blosq / (d__1 * d__1));
	    } else {
		term = (blosq - dstsq) / (blosq + dstsq);
/* Computing 2nd power */
		d__1 = blosq + dstsq;
		chain = weigh * -8. * term * (blosq / (d__1 * d__1));
	    }
/* Computing 2nd power */
	    d__1 = term;
	    error = weigh * (d__1 * d__1);
	    gx = dx * chain;
	    gy = dy * chain;
	    gz = dz * chain;
	    ret_val += error;
	    derivs_ref(1, i__) = derivs_ref(1, i__) + gx;
	    derivs_ref(2, i__) = derivs_ref(2, i__) + gy;
	    derivs_ref(3, i__) = derivs_ref(3, i__) + gz;
	    derivs_ref(1, j) = derivs_ref(1, j) - gx;
	    derivs_ref(2, j) = derivs_ref(2, j) - gy;
	    derivs_ref(3, j) = derivs_ref(3, j) - gz;
	}
    }
    return ret_val;
} /* bnderr_ */

#undef derivs_ref
#undef idfix_ref
#undef dfix_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  function vdwerr  --  computes van der Waals bound error  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "vdwerr" is the hard sphere van der Waals bound error function */
/*     and derivatives that penalizes close nonbonded contacts, */
/*     pairwise neighbors are generated via the method of lights */


doublereal vdwerr_(doublereal *derivs)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal dx, dy, dz, gx, gy, xi, yi, zi, gz;
    static integer kgy, kgz;
    static doublereal radi;
    static integer skip[1000];
    static doublereal term, chain, scale, radsq, blosq, error, dstsq, xsort[
	    200000], ysort[200000], zsort[200000];
    extern /* Subroutine */ int lights_(doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *);


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]
#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1]



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
/*     ##  couple.i  --  near-neighbor atom connectivity lists   ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxn13   maximum number of atoms 1-3 connected to an atom */
/*     maxn14   maximum number of atoms 1-4 connected to an atom */
/*     maxn15   maximum number of atoms 1-5 connected to an atom */

/*     n12      number of atoms directly bonded to each atom */
/*     i12      atom numbers of atoms 1-2 connected to each atom */
/*     n13      number of atoms in a 1-3 relation to each atom */
/*     i13      atom numbers of atoms 1-3 connected to each atom */
/*     n14      number of atoms in a 1-4 relation to each atom */
/*     i14      atom numbers of atoms 1-4 connected to each atom */
/*     n15      number of atoms in a 1-5 relation to each atom */
/*     i15      atom numbers of atoms 1-5 connected to each atom */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  light.i  --  indices for method of lights pair neighbors  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     nlight  total number of sites for method of lights calculation */
/*     kbx     low index of neighbors of each site in the x-sorted list */
/*     kby     low index of neighbors of each site in the y-sorted list */
/*     kbz     low index of neighbors of each site in the z-sorted list */
/*     kex     high index of neighbors of each site in the x-sorted list */
/*     key     high index of neighbors of each site in the y-sorted list */
/*     kez     high index of neighbors of each site in the z-sorted list */
/*     locx    pointer from x-sorted list into original interaction list */
/*     locy    pointer from y-sorted list into original interaction list */
/*     locz    pointer from z-sorted list into original interaction list */
/*     rgx     pointer from original interaction list into x-sorted list */
/*     rgy     pointer from original interaction list into y-sorted list */
/*     rgz     pointer from original interaction list into z-sorted list */




/*     zero out the distance van der Waals error function */

    /* Parameter adjustments */
    derivs -= 4;

    /* Function Body */
    ret_val = 0.;
    scale = 1.;

/*     transfer coordinates and zero out atoms to be skipped */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xsort[i__ - 1] = atoms_1.x[i__ - 1];
	ysort[i__ - 1] = atoms_1.y[i__ - 1];
	zsort[i__ - 1] = atoms_1.z__[i__ - 1];
	skip[i__ - 1] = 0;
    }

/*     use the method of lights to generate neighbors */

    lights_(&disgeo_1.vdwmax, &atoms_1.n, xsort, ysort, zsort);

/*     now, loop over all atoms computing the interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	radi = disgeo_1.vdwrad[i__ - 1];
	xi = xsort[light_1.rgx[i__ - 1] - 1];
	yi = ysort[light_1.rgy[i__ - 1] - 1];
	zi = zsort[light_1.rgz[i__ - 1] - 1];
	i__2 = couple_1.n12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    skip[i12_ref(j, i__) - 1] = i__;
	}
	i__2 = couple_1.n13[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    skip[i13_ref(j, i__) - 1] = i__;
	}
	i__2 = couple_1.n14[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    skip[i14_ref(j, i__) - 1] = i__;
	}
	i__2 = light_1.kex[i__ - 1];
	for (j = light_1.kbx[i__ - 1] + 1; j <= i__2; ++j) {
	    k = light_1.locx[j - 1];
	    if (skip[k - 1] == i__) {
		goto L10;
	    }
	    kgy = light_1.rgy[k - 1];
	    if (kgy < light_1.kby[i__ - 1] || kgy > light_1.key[i__ - 1]) {
		goto L10;
	    }
	    kgz = light_1.rgz[k - 1];
	    if (kgz < light_1.kbz[i__ - 1] || kgz > light_1.kez[i__ - 1]) {
		goto L10;
	    }

/*     calculate squared distances and bounds */

	    dx = xi - xsort[j - 1];
	    dy = yi - ysort[kgy - 1];
	    dz = zi - zsort[kgz - 1];
	    dstsq = dx * dx + dy * dy + dz * dz;
/* Computing 2nd power */
	    d__1 = radi + disgeo_1.vdwrad[k - 1];
	    radsq = d__1 * d__1;
/* Computing MIN */
	    d__1 = bnd_ref(k, i__), d__2 = bnd_ref(i__, k), d__1 = min(d__1,
		    d__2);
	    blosq = min(d__1,radsq);

/*     error and derivatives for lower bound violation */

	    if (dstsq < blosq) {
		term = (blosq - dstsq) / (blosq + dstsq);
/* Computing 2nd power */
		d__1 = blosq + dstsq;
		chain = scale * -8. * term * (blosq / (d__1 * d__1));
/* Computing 2nd power */
		d__1 = term;
		error = scale * (d__1 * d__1);
		gx = dx * chain;
		gy = dy * chain;
		gz = dz * chain;
		ret_val += error;
		derivs_ref(1, i__) = derivs_ref(1, i__) + gx;
		derivs_ref(2, i__) = derivs_ref(2, i__) + gy;
		derivs_ref(3, i__) = derivs_ref(3, i__) + gz;
		derivs_ref(1, k) = derivs_ref(1, k) - gx;
		derivs_ref(2, k) = derivs_ref(2, k) - gy;
		derivs_ref(3, k) = derivs_ref(3, k) - gz;
	    }
L10:
	    ;
	}
    }
    return ret_val;
} /* vdwerr_ */

#undef derivs_ref
#undef bnd_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function locerr  --  computes local geometry error value  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "locerr" is the local geometry error function and derivatives */
/*     including the 1-2, 1-3 and 1-4 distance bound restraints */


doublereal locerr_(doublereal *derivs)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__, ia, ib, ic, id;
    static doublereal dx, dy, dz, gx, gy, gz, term, chain, scale, blosq, 
	    error, bupsq, dstsq;


#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define itors_ref(a_1,a_2) tors_1.itors[(a_2)*4 + a_1 - 5]
#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1]



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
/*     ##  angle.i  --  bond angles within the current structure  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     ak       harmonic angle force constant (kcal/mole/rad**2) */
/*     anat     ideal bond angle or phase shift angle (degrees) */
/*     afld     periodicity for Fourier bond angle term */
/*     nangle   total number of bond angles in the system */
/*     iang     numbers of the atoms in each bond angle */




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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  bond.i  --  covalent bonds in the current structure  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     bk      bond stretch force constants (kcal/mole/Ang**2) */
/*     bl      ideal bond length values in Angstroms */
/*     nbond   total number of bond stretches in the system */
/*     ibnd    numbers of the atoms in each bond stretch */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  tors.i  --  torsional angles within the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     tors1   1-fold amplitude and phase for each torsional angle */
/*     tors2   2-fold amplitude and phase for each torsional angle */
/*     tors3   3-fold amplitude and phase for each torsional angle */
/*     tors4   4-fold amplitude and phase for each torsional angle */
/*     tors5   5-fold amplitude and phase for each torsional angle */
/*     tors6   6-fold amplitude and phase for each torsional angle */
/*     ntors   total number of torsional angles in the system */
/*     itors   numbers of the atoms in each torsional angle */




/*     zero out the local geometry error function */

    /* Parameter adjustments */
    derivs -= 4;

    /* Function Body */
    ret_val = 0.;
    scale = 10.;

/*     calculate the bounds error for bond length distances */

    i__1 = bond_1.nbond;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	i__2 = ibnd_ref(1, i__), i__3 = ibnd_ref(2, i__);
	ia = min(i__2,i__3);
/* Computing MAX */
	i__2 = ibnd_ref(1, i__), i__3 = ibnd_ref(2, i__);
	ib = max(i__2,i__3);
	dx = atoms_1.x[ia - 1] - atoms_1.x[ib - 1];
	dy = atoms_1.y[ia - 1] - atoms_1.y[ib - 1];
	dz = atoms_1.z__[ia - 1] - atoms_1.z__[ib - 1];
	dstsq = dx * dx + dy * dy + dz * dz;
	bupsq = bnd_ref(ia, ib);
	blosq = bnd_ref(ib, ia);
	if (dstsq > bupsq) {
	    term = (dstsq - bupsq) / bupsq;
	    chain = scale * 4. * term / bupsq;
/* Computing 2nd power */
	    d__1 = term;
	    error = scale * (d__1 * d__1);
	    gx = dx * chain;
	    gy = dy * chain;
	    gz = dz * chain;
	    ret_val += error;
	    derivs_ref(1, ia) = derivs_ref(1, ia) + gx;
	    derivs_ref(2, ia) = derivs_ref(2, ia) + gy;
	    derivs_ref(3, ia) = derivs_ref(3, ia) + gz;
	    derivs_ref(1, ib) = derivs_ref(1, ib) - gx;
	    derivs_ref(2, ib) = derivs_ref(2, ib) - gy;
	    derivs_ref(3, ib) = derivs_ref(3, ib) - gz;
	} else if (dstsq < blosq) {
	    term = (blosq - dstsq) / (blosq + dstsq);
/* Computing 2nd power */
	    d__1 = blosq + dstsq;
	    chain = scale * -8. * term * (blosq / (d__1 * d__1));
/* Computing 2nd power */
	    d__1 = term;
	    error = scale * (d__1 * d__1);
	    gx = dx * chain;
	    gy = dy * chain;
	    gz = dz * chain;
	    ret_val += error;
	    derivs_ref(1, ia) = derivs_ref(1, ia) + gx;
	    derivs_ref(2, ia) = derivs_ref(2, ia) + gy;
	    derivs_ref(3, ia) = derivs_ref(3, ia) + gz;
	    derivs_ref(1, ib) = derivs_ref(1, ib) - gx;
	    derivs_ref(2, ib) = derivs_ref(2, ib) - gy;
	    derivs_ref(3, ib) = derivs_ref(3, ib) - gz;
	}
    }

/*     calculate the bounds error for the bond angle distances */

    i__1 = angle_1.nangle;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	i__2 = iang_ref(1, i__), i__3 = iang_ref(3, i__);
	ia = min(i__2,i__3);
/* Computing MAX */
	i__2 = iang_ref(1, i__), i__3 = iang_ref(3, i__);
	ic = max(i__2,i__3);
	dx = atoms_1.x[ia - 1] - atoms_1.x[ic - 1];
	dy = atoms_1.y[ia - 1] - atoms_1.y[ic - 1];
	dz = atoms_1.z__[ia - 1] - atoms_1.z__[ic - 1];
	dstsq = dx * dx + dy * dy + dz * dz;
	bupsq = bnd_ref(ia, ic);
	blosq = bnd_ref(ic, ia);
	if (dstsq > bupsq) {
	    term = (dstsq - bupsq) / bupsq;
	    chain = scale * 4. * term / bupsq;
/* Computing 2nd power */
	    d__1 = term;
	    error = scale * (d__1 * d__1);
	    gx = dx * chain;
	    gy = dy * chain;
	    gz = dz * chain;
	    ret_val += error;
	    derivs_ref(1, ia) = derivs_ref(1, ia) + gx;
	    derivs_ref(2, ia) = derivs_ref(2, ia) + gy;
	    derivs_ref(3, ia) = derivs_ref(3, ia) + gz;
	    derivs_ref(1, ic) = derivs_ref(1, ic) - gx;
	    derivs_ref(2, ic) = derivs_ref(2, ic) - gy;
	    derivs_ref(3, ic) = derivs_ref(3, ic) - gz;
	} else if (dstsq < blosq) {
	    term = (blosq - dstsq) / (blosq + dstsq);
/* Computing 2nd power */
	    d__1 = blosq + dstsq;
	    chain = scale * -8. * term * (blosq / (d__1 * d__1));
/* Computing 2nd power */
	    d__1 = term;
	    error = scale * (d__1 * d__1);
	    gx = dx * chain;
	    gy = dy * chain;
	    gz = dz * chain;
	    ret_val += error;
	    derivs_ref(1, ia) = derivs_ref(1, ia) + gx;
	    derivs_ref(2, ia) = derivs_ref(2, ia) + gy;
	    derivs_ref(3, ia) = derivs_ref(3, ia) + gz;
	    derivs_ref(1, ic) = derivs_ref(1, ic) - gx;
	    derivs_ref(2, ic) = derivs_ref(2, ic) - gy;
	    derivs_ref(3, ic) = derivs_ref(3, ic) - gz;
	}
    }

/*     calculate the bounds error for the torsion angle distances */

    i__1 = tors_1.ntors;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	i__2 = itors_ref(1, i__), i__3 = itors_ref(4, i__);
	ia = min(i__2,i__3);
/* Computing MAX */
	i__2 = itors_ref(1, i__), i__3 = itors_ref(4, i__);
	id = max(i__2,i__3);
	dx = atoms_1.x[ia - 1] - atoms_1.x[id - 1];
	dy = atoms_1.y[ia - 1] - atoms_1.y[id - 1];
	dz = atoms_1.z__[ia - 1] - atoms_1.z__[id - 1];
	dstsq = dx * dx + dy * dy + dz * dz;
	bupsq = bnd_ref(ia, id);
	blosq = bnd_ref(id, ia);
	if (dstsq > bupsq) {
	    term = (dstsq - bupsq) / bupsq;
	    chain = scale * 4. * term / bupsq;
/* Computing 2nd power */
	    d__1 = term;
	    error = scale * (d__1 * d__1);
	    gx = dx * chain;
	    gy = dy * chain;
	    gz = dz * chain;
	    ret_val += error;
	    derivs_ref(1, ia) = derivs_ref(1, ia) + gx;
	    derivs_ref(2, ia) = derivs_ref(2, ia) + gy;
	    derivs_ref(3, ia) = derivs_ref(3, ia) + gz;
	    derivs_ref(1, id) = derivs_ref(1, id) - gx;
	    derivs_ref(2, id) = derivs_ref(2, id) - gy;
	    derivs_ref(3, id) = derivs_ref(3, id) - gz;
	} else if (dstsq < blosq) {
	    term = (blosq - dstsq) / (blosq + dstsq);
/* Computing 2nd power */
	    d__1 = blosq + dstsq;
	    chain = scale * -8. * term * (blosq / (d__1 * d__1));
/* Computing 2nd power */
	    d__1 = term;
	    error = scale * (d__1 * d__1);
	    gx = dx * chain;
	    gy = dy * chain;
	    gz = dz * chain;
	    ret_val += error;
	    derivs_ref(1, ia) = derivs_ref(1, ia) + gx;
	    derivs_ref(2, ia) = derivs_ref(2, ia) + gy;
	    derivs_ref(3, ia) = derivs_ref(3, ia) + gz;
	    derivs_ref(1, id) = derivs_ref(1, id) - gx;
	    derivs_ref(2, id) = derivs_ref(2, id) - gy;
	    derivs_ref(3, id) = derivs_ref(3, id) - gz;
	}
    }
    return ret_val;
} /* locerr_ */

#undef derivs_ref
#undef itors_ref
#undef iang_ref
#undef ibnd_ref
#undef bnd_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function chirer  --  computes chirality error function  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "chirer" computes the chirality error and its derivatives */
/*     with respect to atomic Cartesian coordinates as a sum the */
/*     squares of deviations of chiral volumes from target values */


doublereal chirer_(doublereal *derivs)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__;
    static doublereal c1, c2, c3;
    static integer ia, ib, ic, id;
    static doublereal dt, xad, yad, zad, xbd, ybd, zbd, xcd, ycd, zcd, vol, 
	    dedt, scale, error, dedxia, dedyia, dedzia, dedxib, dedyib, 
	    dedzib, dedxic, dedyic, dedzic, dedxid, dedyid, dedzid;


#define chir_ref(a_1,a_2) kgeoms_1.chir[(a_2)*3 + a_1 - 4]
#define ichir_ref(a_1,a_2) kgeoms_1.ichir[(a_2)*4 + a_1 - 5]
#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1]



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
/*     ##  kgeoms.i  --  parameters for the geometrical restraints  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xpfix      x-coordinate target for each restrained position */
/*     ypfix      y-coordinate target for each restrained position */
/*     zpfix      z-coordinate target for each restrained position */
/*     pfix       force constant and flat-well range for each position */
/*     dfix       force constant and target range for each distance */
/*     afix       force constant and target range for each angle */
/*     tfix       force constant and target range for each torsion */
/*     gfix       force constant and target range for each group distance */
/*     chir       force constant and target range for chiral centers */
/*     depth      depth of shallow Gaussian basin restraint */
/*     width      exponential width coefficient of Gaussian basin */
/*     rwall      radius of spherical droplet boundary restraint */
/*     npfix      number of position restraints to be applied */
/*     ipfix      atom number involved in each position restraint */
/*     kpfix      flags to use x-, y-, z-coordinate position restraints */
/*     ndfix      number of distance restraints to be applied */
/*     idfix      atom numbers defining each distance restraint */
/*     nafix      number of angle restraints to be applied */
/*     iafix      atom numbers defining each angle restraint */
/*     ntfix      number of torsional restraints to be applied */
/*     itfix      atom numbers defining each torsional restraint */
/*     ngfix      number of group distance restraints to be applied */
/*     igfix      group numbers defining each group distance restraint */
/*     nchir      number of chirality restraints to be applied */
/*     ichir      atom numbers defining each chirality restraint */
/*     use_basin  logical flag governing use of Gaussian basin */
/*     use_wall   logical flag governing use of droplet boundary */




/*     zero the chirality restraint error function */

    /* Parameter adjustments */
    derivs -= 4;

    /* Function Body */
    ret_val = 0.;
    scale = .1;

/*     find signed volume value and compute chirality error */

    i__1 = kgeoms_1.nchir;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ichir_ref(1, i__);
	ib = ichir_ref(2, i__);
	ic = ichir_ref(3, i__);
	id = ichir_ref(4, i__);
	xad = atoms_1.x[ia - 1] - atoms_1.x[id - 1];
	yad = atoms_1.y[ia - 1] - atoms_1.y[id - 1];
	zad = atoms_1.z__[ia - 1] - atoms_1.z__[id - 1];
	xbd = atoms_1.x[ib - 1] - atoms_1.x[id - 1];
	ybd = atoms_1.y[ib - 1] - atoms_1.y[id - 1];
	zbd = atoms_1.z__[ib - 1] - atoms_1.z__[id - 1];
	xcd = atoms_1.x[ic - 1] - atoms_1.x[id - 1];
	ycd = atoms_1.y[ic - 1] - atoms_1.y[id - 1];
	zcd = atoms_1.z__[ic - 1] - atoms_1.z__[id - 1];
	c1 = ybd * zcd - zbd * ycd;
	c2 = ycd * zad - zcd * yad;
	c3 = yad * zbd - zad * ybd;
	vol = xad * c1 + xbd * c2 + xcd * c3;
	dt = vol - chir_ref(2, i__);
/* Computing 2nd power */
	d__1 = dt;
	error = scale * (d__1 * d__1);
	dedt = scale * 2. * dt;

/*     compute derivative components for this interaction */

	dedxia = dedt * (ybd * zcd - zbd * ycd);
	dedyia = dedt * (zbd * xcd - xbd * zcd);
	dedzia = dedt * (xbd * ycd - ybd * xcd);
	dedxib = dedt * (zad * ycd - yad * zcd);
	dedyib = dedt * (xad * zcd - zad * xcd);
	dedzib = dedt * (yad * xcd - xad * ycd);
	dedxic = dedt * (yad * zbd - zad * ybd);
	dedyic = dedt * (zad * xbd - xad * zbd);
	dedzic = dedt * (xad * ybd - yad * xbd);
	dedxid = -dedxia - dedxib - dedxic;
	dedyid = -dedyia - dedyib - dedyic;
	dedzid = -dedzia - dedzib - dedzic;

/*     increment the chirality restraint error and derivatives */

	ret_val += error;
	derivs_ref(1, ia) = derivs_ref(1, ia) + dedxia;
	derivs_ref(2, ia) = derivs_ref(2, ia) + dedyia;
	derivs_ref(3, ia) = derivs_ref(3, ia) + dedzia;
	derivs_ref(1, ib) = derivs_ref(1, ib) + dedxib;
	derivs_ref(2, ib) = derivs_ref(2, ib) + dedyib;
	derivs_ref(3, ib) = derivs_ref(3, ib) + dedzib;
	derivs_ref(1, ic) = derivs_ref(1, ic) + dedxic;
	derivs_ref(2, ic) = derivs_ref(2, ic) + dedyic;
	derivs_ref(3, ic) = derivs_ref(3, ic) + dedzic;
	derivs_ref(1, id) = derivs_ref(1, id) + dedxid;
	derivs_ref(2, id) = derivs_ref(2, id) + dedyid;
	derivs_ref(3, id) = derivs_ref(3, id) + dedzid;
    }
    return ret_val;
} /* chirer_ */

#undef derivs_ref
#undef ichir_ref
#undef chir_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function torser  --  computes torsional error function  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "torser" computes the torsional error function and its first */
/*     derivatives with respect to the atomic Cartesian coordinates */
/*     based on the deviation of specified torsional angles from */
/*     desired values, the contained bond angles are also restrained */
/*     to avoid a numerical instability */


doublereal torser_(doublereal *derivs)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal t1, t2;
    static integer ia, ib, ic, id;
    static doublereal dt, rp, xp, yp, zp, xt, yt, zt, xu, yu, zu, tf1, tf2, 
	    rt2, ru2, rba, rcb, rdc, xba, yba, zba, xcb, ycb, zcb, xdc, xia, 
	    yia, zia, xib, yib, zib, xic, yic, zic, xid, yid, zid, dot, ydc, 
	    zdc, xca, yca, zca, xdb, ydb, zdb, xtu, ytu, ztu;
    static integer iref;
    static doublereal rrba, rrcb, rrdc, sine, rrca, rrdb, xria, yria, zria, 
	    xrib, yrib, zrib, xric, yric, zric, xrid, yrid, zrid, rtru, deddt,
	     angle, scale, bamax, bamin, cbmax, cbmin, force, dcmax, dcmin, 
	    terma, termb, termc, termd, error, camax, camin, dbmax, dbmin, 
	    dedxt, dedyt, dedzt, dedxu, dedyu, dedzu, bndfac, angfac;
    static logical bonded;
    static doublereal dedphi, dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, 
	    dedxic, dedyic, dedzic, angmin, angmax, dedxid, dedyid, cosine, 
	    dedzid, target, cosmax, cosmin;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]
#define xref_ref(a_1,a_2) refer_1.xref[(a_2)*25000 + a_1 - 25001]
#define yref_ref(a_1,a_2) refer_1.yref[(a_2)*25000 + a_1 - 25001]
#define zref_ref(a_1,a_2) refer_1.zref[(a_2)*25000 + a_1 - 25001]
#define tfix_ref(a_1,a_2) kgeoms_1.tfix[(a_2)*3 + a_1 - 4]
#define itfix_ref(a_1,a_2) kgeoms_1.itfix[(a_2)*4 + a_1 - 5]
#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1]



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
/*     ##  couple.i  --  near-neighbor atom connectivity lists   ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxn13   maximum number of atoms 1-3 connected to an atom */
/*     maxn14   maximum number of atoms 1-4 connected to an atom */
/*     maxn15   maximum number of atoms 1-5 connected to an atom */

/*     n12      number of atoms directly bonded to each atom */
/*     i12      atom numbers of atoms 1-2 connected to each atom */
/*     n13      number of atoms in a 1-3 relation to each atom */
/*     i13      atom numbers of atoms 1-3 connected to each atom */
/*     n14      number of atoms in a 1-4 relation to each atom */
/*     i14      atom numbers of atoms 1-4 connected to each atom */
/*     n15      number of atoms in a 1-5 relation to each atom */
/*     i15      atom numbers of atoms 1-5 connected to each atom */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  kgeoms.i  --  parameters for the geometrical restraints  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xpfix      x-coordinate target for each restrained position */
/*     ypfix      y-coordinate target for each restrained position */
/*     zpfix      z-coordinate target for each restrained position */
/*     pfix       force constant and flat-well range for each position */
/*     dfix       force constant and target range for each distance */
/*     afix       force constant and target range for each angle */
/*     tfix       force constant and target range for each torsion */
/*     gfix       force constant and target range for each group distance */
/*     chir       force constant and target range for chiral centers */
/*     depth      depth of shallow Gaussian basin restraint */
/*     width      exponential width coefficient of Gaussian basin */
/*     rwall      radius of spherical droplet boundary restraint */
/*     npfix      number of position restraints to be applied */
/*     ipfix      atom number involved in each position restraint */
/*     kpfix      flags to use x-, y-, z-coordinate position restraints */
/*     ndfix      number of distance restraints to be applied */
/*     idfix      atom numbers defining each distance restraint */
/*     nafix      number of angle restraints to be applied */
/*     iafix      atom numbers defining each angle restraint */
/*     ntfix      number of torsional restraints to be applied */
/*     itfix      atom numbers defining each torsional restraint */
/*     ngfix      number of group distance restraints to be applied */
/*     igfix      group numbers defining each group distance restraint */
/*     nchir      number of chirality restraints to be applied */
/*     ichir      atom numbers defining each chirality restraint */
/*     use_basin  logical flag governing use of Gaussian basin */
/*     use_wall   logical flag governing use of droplet boundary */




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




/*     zero the torsional restraint error function */

    /* Parameter adjustments */
    derivs -= 4;

    /* Function Body */
    ret_val = 0.;
    scale = .01;

/*     compute error value and derivs for torsional restraints */

    i__1 = kgeoms_1.ntfix;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = itfix_ref(1, i__);
	ib = itfix_ref(2, i__);
	ic = itfix_ref(3, i__);
	id = itfix_ref(4, i__);
	xia = atoms_1.x[ia - 1];
	yia = atoms_1.y[ia - 1];
	zia = atoms_1.z__[ia - 1];
	xib = atoms_1.x[ib - 1];
	yib = atoms_1.y[ib - 1];
	zib = atoms_1.z__[ib - 1];
	xic = atoms_1.x[ic - 1];
	yic = atoms_1.y[ic - 1];
	zic = atoms_1.z__[ic - 1];
	xid = atoms_1.x[id - 1];
	yid = atoms_1.y[id - 1];
	zid = atoms_1.z__[id - 1];
	xba = xib - xia;
	yba = yib - yia;
	zba = zib - zia;
	xcb = xic - xib;
	ycb = yic - yib;
	zcb = zic - zib;
	xdc = xid - xic;
	ydc = yid - yic;
	zdc = zid - zic;

/*     find the actual distances between the four atoms */

	rba = sqrt(xba * xba + yba * yba + zba * zba);
	rcb = sqrt(xcb * xcb + ycb * ycb + zcb * zcb);
	rdc = sqrt(xdc * xdc + ydc * ydc + zdc * zdc);

/*     see if the torsion involves four directly bonded atoms */

	k = 0;
	i__2 = couple_1.n12[ia - 1];
	for (j = 1; j <= i__2; ++j) {
	    if (i12_ref(j, ia) == ib) {
		++k;
	    }
	}
	i__2 = couple_1.n12[ib - 1];
	for (j = 1; j <= i__2; ++j) {
	    if (i12_ref(j, ib) == ic) {
		++k;
	    }
	}
	i__2 = couple_1.n12[ic - 1];
	for (j = 1; j <= i__2; ++j) {
	    if (i12_ref(j, ic) == id) {
		++k;
	    }
	}
	if (k == 3) {
	    bonded = TRUE_;
	} else {
	    bonded = FALSE_;
	}

/*     get maximum and minimum distances from distance matrix */

	if (bonded) {
	    bamax = sqrt(bnd_ref(min(ib,ia), max(ib,ia)));
	    bamin = sqrt(bnd_ref(max(ib,ia), min(ib,ia)));
	    cbmax = sqrt(bnd_ref(min(ic,ib), max(ic,ib)));
	    cbmin = sqrt(bnd_ref(max(ic,ib), min(ic,ib)));
	    dcmax = sqrt(bnd_ref(min(id,ic), max(id,ic)));
	    dcmin = sqrt(bnd_ref(max(id,ic), min(id,ic)));
	    camax = sqrt(bnd_ref(min(ic,ia), max(ic,ia)));
	    camin = sqrt(bnd_ref(max(ic,ia), min(ic,ia)));
	    dbmax = sqrt(bnd_ref(min(id,ib), max(id,ib)));
	    dbmin = sqrt(bnd_ref(max(id,ib), min(id,ib)));

/*     get maximum and minimum distances from input structure */

	} else {
	    iref = 1;
	    xria = xref_ref(ia, iref);
	    yria = yref_ref(ia, iref);
	    zria = zref_ref(ia, iref);
	    xrib = xref_ref(ib, iref);
	    yrib = yref_ref(ib, iref);
	    zrib = zref_ref(ib, iref);
	    xric = xref_ref(ic, iref);
	    yric = yref_ref(ic, iref);
	    zric = zref_ref(ic, iref);
	    xrid = xref_ref(id, iref);
	    yrid = yref_ref(id, iref);
	    zrid = zref_ref(id, iref);
/* Computing 2nd power */
	    d__1 = xrib - xria;
/* Computing 2nd power */
	    d__2 = yrib - yria;
/* Computing 2nd power */
	    d__3 = zrib - zria;
	    rrba = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
	    d__1 = xric - xrib;
/* Computing 2nd power */
	    d__2 = yric - yrib;
/* Computing 2nd power */
	    d__3 = zric - zrib;
	    rrcb = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
	    d__1 = xrid - xric;
/* Computing 2nd power */
	    d__2 = yrid - yric;
/* Computing 2nd power */
	    d__3 = zrid - zric;
	    rrdc = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
	    d__1 = xric - xria;
/* Computing 2nd power */
	    d__2 = yric - yria;
/* Computing 2nd power */
	    d__3 = zric - zria;
	    rrca = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
	    d__1 = xrid - xrib;
/* Computing 2nd power */
	    d__2 = yrid - yrib;
/* Computing 2nd power */
	    d__3 = zrid - zrib;
	    rrdb = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	    bndfac = .05;
	    angfac = .05;
	    bamax = (bndfac + 1.) * rrba;
	    bamin = (1. - bndfac) * rrba;
	    cbmax = (bndfac + 1.) * rrcb;
	    cbmin = (1. - bndfac) * rrcb;
	    dcmax = (bndfac + 1.) * rrdc;
	    dcmin = (1. - bndfac) * rrdc;
	    camax = (angfac + 1.) * rrca;
	    camin = (1. - angfac) * rrca;
	    dbmax = (angfac + 1.) * rrdb;
	    dbmin = (1. - angfac) * rrdb;
	}

/*     compute the ia-ib-ic bond angle and any error */

	dot = xba * xcb + yba * ycb + zba * zcb;
	cosine = -dot / (rba * rcb);
/* Computing MIN */
	d__1 = 1., d__2 = max(-1.,cosine);
	cosine = min(d__1,d__2);
	angle = acos(cosine) * 57.29577951308232088;
/* Computing 2nd power */
	d__1 = bamin;
/* Computing 2nd power */
	d__2 = cbmin;
/* Computing 2nd power */
	d__3 = camax;
	cosmax = (d__1 * d__1 + d__2 * d__2 - d__3 * d__3) / (bamin * 2. * 
		cbmin);
/* Computing MIN */
	d__1 = 1., d__2 = max(-1.,cosmax);
	cosmax = min(d__1,d__2);
	angmax = acos(cosmax) * 57.29577951308232088;
/* Computing 2nd power */
	d__1 = bamax;
/* Computing 2nd power */
	d__2 = cbmax;
/* Computing 2nd power */
	d__3 = camin;
	cosmin = (d__1 * d__1 + d__2 * d__2 - d__3 * d__3) / (bamax * 2. * 
		cbmax);
/* Computing MIN */
	d__1 = 1., d__2 = max(-1.,cosmin);
	cosmin = min(d__1,d__2);
	angmin = acos(cosmin) * 57.29577951308232088;
	if (angle > angmax) {
	    dt = angle - angmax;
	} else if (angle < angmin) {
	    dt = angle - angmin;
	} else {
	    dt = 0.;
	}
/* Computing 2nd power */
	d__1 = dt;
	error = scale * (d__1 * d__1);
	deddt = scale * 114.59155902616465 * dt;

/*     compute derivative components for this interaction */

	xp = zcb * yba - ycb * zba;
	yp = xcb * zba - zcb * xba;
	zp = ycb * xba - xcb * yba;
	rp = sqrt(xp * xp + yp * yp + zp * zp);
	if (rp != 0.) {
	    terma = -deddt / (rba * rba * rp);
	    termc = deddt / (rcb * rcb * rp);
	    dedxia = terma * (zba * yp - yba * zp);
	    dedyia = terma * (xba * zp - zba * xp);
	    dedzia = terma * (yba * xp - xba * yp);
	    dedxic = termc * (ycb * zp - zcb * yp);
	    dedyic = termc * (zcb * xp - xcb * zp);
	    dedzic = termc * (xcb * yp - ycb * xp);
	    dedxib = -dedxia - dedxic;
	    dedyib = -dedyia - dedyic;
	    dedzib = -dedzia - dedzic;

/*     increment the bond angle restraint error and derivatives */

	    ret_val += error;
	    derivs_ref(1, ia) = derivs_ref(1, ia) + dedxia;
	    derivs_ref(2, ia) = derivs_ref(2, ia) + dedyia;
	    derivs_ref(3, ia) = derivs_ref(3, ia) + dedzia;
	    derivs_ref(1, ib) = derivs_ref(1, ib) + dedxib;
	    derivs_ref(2, ib) = derivs_ref(2, ib) + dedyib;
	    derivs_ref(3, ib) = derivs_ref(3, ib) + dedzib;
	    derivs_ref(1, ic) = derivs_ref(1, ic) + dedxic;
	    derivs_ref(2, ic) = derivs_ref(2, ic) + dedyic;
	    derivs_ref(3, ic) = derivs_ref(3, ic) + dedzic;
	}

/*     compute the ib-ic-id bond angle and any error */

	dot = xdc * xcb + ydc * ycb + zdc * zcb;
	cosine = -dot / (rdc * rcb);
/* Computing MIN */
	d__1 = 1., d__2 = max(-1.,cosine);
	cosine = min(d__1,d__2);
	angle = acos(cosine) * 57.29577951308232088;
/* Computing 2nd power */
	d__1 = dcmin;
/* Computing 2nd power */
	d__2 = cbmin;
/* Computing 2nd power */
	d__3 = dbmax;
	cosmax = (d__1 * d__1 + d__2 * d__2 - d__3 * d__3) / (dcmin * 2. * 
		cbmin);
/* Computing MIN */
	d__1 = 1., d__2 = max(-1.,cosmax);
	cosmax = min(d__1,d__2);
	angmax = acos(cosmax) * 57.29577951308232088;
/* Computing 2nd power */
	d__1 = dcmax;
/* Computing 2nd power */
	d__2 = cbmax;
/* Computing 2nd power */
	d__3 = dbmin;
	cosmin = (d__1 * d__1 + d__2 * d__2 - d__3 * d__3) / (dcmax * 2. * 
		cbmax);
/* Computing MIN */
	d__1 = 1., d__2 = max(-1.,cosmin);
	cosmax = min(d__1,d__2);
	angmin = acos(cosmin) * 57.29577951308232088;
	if (angle > angmax) {
	    dt = angle - angmax;
	} else if (angle < angmin) {
	    dt = angle - angmin;
	} else {
	    dt = 0.;
	}
/* Computing 2nd power */
	d__1 = dt;
	error = scale * (d__1 * d__1);
	deddt = scale * 114.59155902616465 * dt;

/*     compute derivative components for this interaction */

	xp = zdc * ycb - ydc * zcb;
	yp = xdc * zcb - zdc * xcb;
	zp = ydc * xcb - xdc * ycb;
	rp = sqrt(xp * xp + yp * yp + zp * zp);
	if (rp != 0.) {
	    termb = -deddt / (rcb * rcb * rp);
	    termd = deddt / (rdc * rdc * rp);
	    dedxib = termb * (zcb * yp - ycb * zp);
	    dedyib = termb * (xcb * zp - zcb * xp);
	    dedzib = termb * (ycb * xp - xcb * yp);
	    dedxid = termd * (ydc * zp - zdc * yp);
	    dedyid = termd * (zdc * xp - xdc * zp);
	    dedzid = termd * (xdc * yp - ydc * xp);
	    dedxic = -dedxib - dedxid;
	    dedyic = -dedyib - dedyid;
	    dedzic = -dedzib - dedzid;

/*     increment the bond angle restraint error and derivatives */

	    ret_val += error;
	    derivs_ref(1, ib) = derivs_ref(1, ib) + dedxib;
	    derivs_ref(2, ib) = derivs_ref(2, ib) + dedyib;
	    derivs_ref(3, ib) = derivs_ref(3, ib) + dedzib;
	    derivs_ref(1, ic) = derivs_ref(1, ic) + dedxic;
	    derivs_ref(2, ic) = derivs_ref(2, ic) + dedyic;
	    derivs_ref(3, ic) = derivs_ref(3, ic) + dedzic;
	    derivs_ref(1, id) = derivs_ref(1, id) + dedxid;
	    derivs_ref(2, id) = derivs_ref(2, id) + dedyid;
	    derivs_ref(3, id) = derivs_ref(3, id) + dedzid;
	}

/*     compute the value of the ia-ib-ic-id torsional angle */

	xt = yba * zcb - ycb * zba;
	yt = zba * xcb - zcb * xba;
	zt = xba * ycb - xcb * yba;
	xu = ycb * zdc - ydc * zcb;
	yu = zcb * xdc - zdc * xcb;
	zu = xcb * ydc - xdc * ycb;
	xtu = yt * zu - yu * zt;
	ytu = zt * xu - zu * xt;
	ztu = xt * yu - xu * yt;
	rt2 = xt * xt + yt * yt + zt * zt;
	ru2 = xu * xu + yu * yu + zu * zu;
	rtru = sqrt(rt2 * ru2);
	if (rtru != 0.) {
	    rcb = sqrt(xcb * xcb + ycb * ycb + zcb * zcb);
	    cosine = (xt * xu + yt * yu + zt * zu) / rtru;
	    sine = (xcb * xtu + ycb * ytu + zcb * ztu) / (rcb * rtru);
/* Computing MIN */
	    d__1 = 1., d__2 = max(-1.,cosine);
	    cosine = min(d__1,d__2);
	    angle = acos(cosine) * 57.29577951308232088;
	    if (sine < 0.) {
		angle = -angle;
	    }

/*     calculate the torsional restraint error for this angle */

	    force = tfix_ref(1, i__);
	    tf1 = tfix_ref(2, i__);
	    tf2 = tfix_ref(3, i__);
	    if (angle > tf1 && angle < tf2) {
		target = angle;
	    } else if (angle > tf1 && tf1 > tf2) {
		target = angle;
	    } else if (angle < tf2 && tf1 > tf2) {
		target = angle;
	    } else {
		t1 = angle - tf1;
		t2 = angle - tf2;
		if (t1 > 180.) {
		    t1 += -360.;
		} else if (t1 < -180.) {
		    t1 += 360.;
		}
		if (t2 > 180.) {
		    t2 += -360.;
		} else if (t2 < -180.) {
		    t2 += 360.;
		}
		if (abs(t1) < abs(t2)) {
		    target = tf1;
		} else {
		    target = tf2;
		}
	    }
	    dt = angle - target;
	    if (dt > 180.) {
		dt += -360.;
	    } else if (dt < -180.) {
		dt += 360.;
	    }
/* Computing 2nd power */
	    d__1 = dt;
	    error = scale * force * (d__1 * d__1);
	    dedphi = scale * 114.59155902616465 * force * dt;

/*     chain rule terms for first derivative components */

	    xca = xic - xia;
	    yca = yic - yia;
	    zca = zic - zia;
	    xdb = xid - xib;
	    ydb = yid - yib;
	    zdb = zid - zib;
	    dedxt = dedphi * (yt * zcb - ycb * zt) / (rt2 * rcb);
	    dedyt = dedphi * (zt * xcb - zcb * xt) / (rt2 * rcb);
	    dedzt = dedphi * (xt * ycb - xcb * yt) / (rt2 * rcb);
	    dedxu = -dedphi * (yu * zcb - ycb * zu) / (ru2 * rcb);
	    dedyu = -dedphi * (zu * xcb - zcb * xu) / (ru2 * rcb);
	    dedzu = -dedphi * (xu * ycb - xcb * yu) / (ru2 * rcb);

/*     compute first derivative components for torsion angle */

	    dedxia = zcb * dedyt - ycb * dedzt;
	    dedyia = xcb * dedzt - zcb * dedxt;
	    dedzia = ycb * dedxt - xcb * dedyt;
	    dedxib = yca * dedzt - zca * dedyt + zdc * dedyu - ydc * dedzu;
	    dedyib = zca * dedxt - xca * dedzt + xdc * dedzu - zdc * dedxu;
	    dedzib = xca * dedyt - yca * dedxt + ydc * dedxu - xdc * dedyu;
	    dedxic = zba * dedyt - yba * dedzt + ydb * dedzu - zdb * dedyu;
	    dedyic = xba * dedzt - zba * dedxt + zdb * dedxu - xdb * dedzu;
	    dedzic = yba * dedxt - xba * dedyt + xdb * dedyu - ydb * dedxu;
	    dedxid = zcb * dedyu - ycb * dedzu;
	    dedyid = xcb * dedzu - zcb * dedxu;
	    dedzid = ycb * dedxu - xcb * dedyu;

/*     increment the torsional restraint error and derivatives */

	    ret_val += error;
	    derivs_ref(1, ia) = derivs_ref(1, ia) + dedxia;
	    derivs_ref(2, ia) = derivs_ref(2, ia) + dedyia;
	    derivs_ref(3, ia) = derivs_ref(3, ia) + dedzia;
	    derivs_ref(1, ib) = derivs_ref(1, ib) + dedxib;
	    derivs_ref(2, ib) = derivs_ref(2, ib) + dedyib;
	    derivs_ref(3, ib) = derivs_ref(3, ib) + dedzib;
	    derivs_ref(1, ic) = derivs_ref(1, ic) + dedxic;
	    derivs_ref(2, ic) = derivs_ref(2, ic) + dedyic;
	    derivs_ref(3, ic) = derivs_ref(3, ic) + dedzic;
	    derivs_ref(1, id) = derivs_ref(1, id) + dedxid;
	    derivs_ref(2, id) = derivs_ref(2, id) + dedyid;
	    derivs_ref(3, id) = derivs_ref(3, id) + dedzid;
	}
    }
    return ret_val;
} /* torser_ */

#undef derivs_ref
#undef itfix_ref
#undef tfix_ref
#undef zref_ref
#undef yref_ref
#undef xref_ref
#undef bnd_ref
#undef i12_ref


