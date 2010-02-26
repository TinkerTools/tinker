/* spacefill.f -- translated by f2c (version 20050501).
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
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

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
    doublereal rad[5000], eps[5000], rad4[5000], eps4[5000], reduct[5000];
} kvdws_;

#define kvdws_1 kvdws_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  program spacefill  --  surface area and volume of model  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "spacefill" computes the surface area and volume of */
/*     a structure; the van der Waals, accessible-excluded, */
/*     and contact-reentrant definitions are available */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Three Types of Area and Volume can be Co"
	    "mputed :\002,//,4x,\002(1) Van der Waals Area and Volume\002,/,4"
	    "x,\002(2) Accessible Area and Excluded Volume\002,/,4x,\002(3) C"
	    "ontact-Reentrant Area and Volume\002)";
    static char fmt_30[] = "(/,\002 Enter the Number of your Choice [1] : "
	    " \002,$)";
    static char fmt_40[] = "(i10)";
    static char fmt_60[] = "(/,\002 Enter a Value for the Probe Radius\002"
	    ",\002 [1.4 Ang] :  \002,$)";
    static char fmt_70[] = "(f20.0)";
    static char fmt_80[] = "(/,\002 Include the Hydrogen Atoms in Computat"
	    "ion\002,\002 [N] :  \002,$)";
    static char fmt_90[] = "(a120)";
    static char fmt_100[] = "(/,\002 Area and Volume for Archive Structure "
	    ":\002,5x,i8)";
    static char fmt_110[] = "(/,\002 Van der Waals Surface Area and Volume "
	    ":\002)";
    static char fmt_120[] = "(/,\002 Accessible Surface Area and Excluded Vo"
	    "lume :\002)";
    static char fmt_130[] = "(/,\002 Contact-Reentrant Surface Area and Volu"
	    "me :\002)";
    static char fmt_140[] = "(/,\002 Total Area :\002,f20.3,\002 Square Angs"
	    "troms\002)";
    static char fmt_150[] = "(\002 Total Volume :\002,f18.3,\002 Cubic Angst"
	    "roms\002)";

    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *), f_rew(alist *), f_clos(cllist *);

    /* Local variables */
    extern integer freeunit_(void);
    extern /* Subroutine */ int connolly_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer i__;
    static doublereal area;
    static integer mode;
    extern /* Subroutine */ int kvdw_(void);
    static integer next, ixyz;
    extern /* Subroutine */ int field_(void), final_(void);
    static integer frame;
    static doublereal probe;
    extern /* Subroutine */ int katom_(void);
    static doublereal value;
    static logical exist, query;
    extern /* Subroutine */ int active_(void);
    static char record[120];
    extern doublereal random_(void);
    static doublereal radius[25000], volume;
    static char answer[1], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), suffix_(char *, char 
	    *, ftnlen, ftnlen), getxyz_(void);
    static doublereal exclude;
    extern /* Subroutine */ int initial_(void), nextarg_(char *, logical *, 
	    ftnlen), gettext_(char *, char *, integer *, ftnlen, ftnlen), 
	    version_(char *, char *, ftnlen, ftnlen), readxyz_(integer *);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static icilist io___5 = { 1, string, 1, 0, 120, 1 };
    static cilist io___7 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_40, 0 };
    static icilist io___12 = { 1, string, 1, 0, 120, 1 };
    static cilist io___13 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_150, 0 };




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
/*     ##  atmtyp.i  --  atomic properties for each current atom  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     mass      atomic weight for each atom in the system */
/*     tag       integer atom labels from input coordinates file */
/*     class     atom class number for each atom in the system */
/*     atomic    atomic number for each atom in the system */
/*     valence   valence number for each atom in the system */
/*     name      atom name for each atom in the system */
/*     story     descriptive type for each atom in system */




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

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kvdws.i  --  forcefield parameters for van der Waals terms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     rad      van der Waals radius parameter for each atom type */
/*     eps      van der Waals well depth parameter for each atom type */
/*     rad4     van der Waals radius parameter in 1-4 interactions */
/*     eps4     van der Waals well depth parameter in 1-4 interactions */
/*     reduct   van der Waals reduction factor for each atom type */




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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     get the Cartesian coordinates for the system */

    initial_();
    getxyz_();

/*     determine the atoms to be used in computation; */
/*     radii can be changed via the keyword mechanism */

    field_();
    active_();
    katom_();
    kvdw_();

/*     initialize random numbers and turn on extra printing */

    value = random_();
    inform_1.verbose = TRUE_;
    inform_1.debug = TRUE_;

/*     select either vdw, excluded or molecular volume and area */

    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___5);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&mode, (ftnlen)sizeof(integer));
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
	io___8.ciunit = iounit_1.iout;
	s_wsfe(&io___8);
	e_wsfe();
	io___9.ciunit = iounit_1.input;
	s_rsfe(&io___9);
	do_fio(&c__1, (char *)&mode, (ftnlen)sizeof(integer));
	e_rsfe();
    }
    if (mode != 2 && mode != 3) {
	mode = 1;
    }

/*     set the excluded/accessible and contact/reentrant probes */

    value = 0.;
    probe = 0.;
    exclude = 0.;
    if (mode == 2 || mode == 3) {
	query = TRUE_;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___12);
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&value, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L50;
	    }
	    query = FALSE_;
	}
L50:
	if (query) {
	    io___13.ciunit = iounit_1.iout;
	    s_wsfe(&io___13);
	    e_wsfe();
	    io___14.ciunit = iounit_1.input;
	    s_rsfe(&io___14);
	    do_fio(&c__1, (char *)&value, (ftnlen)sizeof(doublereal));
	    e_rsfe();
	}
	if (value == 0.) {
	    value = 1.4;
	}
	if (mode == 2) {
	    exclude = value;
	}
	if (mode == 3) {
	    probe = value;
	}
    }

/*     decide whether to include hydrogens in the calculation */

    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___16.ciunit = iounit_1.iout;
	s_wsfe(&io___16);
	e_wsfe();
	io___17.ciunit = iounit_1.input;
	s_rsfe(&io___17);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer != 'Y') {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (atmtyp_1.atomic[i__ - 1] == 1) {
		usage_1.use[i__ - 1] = FALSE_;
	    }
	}
    }

/*     set each atomic radius to the Lennard-Jones sigma value */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    radius[i__ - 1] = kvdws_1.rad[atmtyp_1.class__[i__ - 1] - 1] / 
		    1.122462048309372981;
	} else {
	    radius[i__ - 1] = 0.;
	}
    }

/*     reopen the coordinates file and read the first structure */

    frame = 0;
    ixyz = freeunit_();
    s_copy(xyzfile, files_1.filename, (ftnlen)120, (ftnlen)120);
    suffix_(xyzfile, "xyz", (ftnlen)120, (ftnlen)3);
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
    al__1.aerr = 0;
    al__1.aunit = ixyz;
    f_rew(&al__1);
    readxyz_(&ixyz);

/*     get area and volume for successive coordinate structures */

    while(! inform_1.abort) {
	++frame;
	if (frame > 1) {
	    io___25.ciunit = iounit_1.iout;
	    s_wsfe(&io___25);
	    do_fio(&c__1, (char *)&frame, (ftnlen)sizeof(integer));
	    e_wsfe();
	}

/*     use the Connolly routines to find the area and volume */

	connolly_(&volume, &area, radius, &probe, &exclude);

/*     print out the values of the total surface area and volume */

	if (mode == 1) {
	    io___28.ciunit = iounit_1.iout;
	    s_wsfe(&io___28);
	    e_wsfe();
	} else if (mode == 2) {
	    io___29.ciunit = iounit_1.iout;
	    s_wsfe(&io___29);
	    e_wsfe();
	} else if (mode == 3) {
	    io___30.ciunit = iounit_1.iout;
	    s_wsfe(&io___30);
	    e_wsfe();
	}
	io___31.ciunit = iounit_1.iout;
	s_wsfe(&io___31);
	do_fio(&c__1, (char *)&area, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___32.ciunit = iounit_1.iout;
	s_wsfe(&io___32);
	do_fio(&c__1, (char *)&volume, (ftnlen)sizeof(doublereal));
	e_wsfe();

/*     attempt to read next structure from the coordinate file */

	readxyz_(&ixyz);
    }

/*     perform any final tasks before program exit */

    cl__1.cerr = 0;
    cl__1.cunit = ixyz;
    cl__1.csta = 0;
    f_clos(&cl__1);
    inform_1.debug = FALSE_;
    final_();
    return 0;
} /* MAIN__ */

/* Main program alias */ int spacefill_ () { MAIN__ (); return 0; }
