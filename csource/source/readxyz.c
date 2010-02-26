/* readxyz.f -- translated by f2c (version 20050501).
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
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

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
    integer ltitle;
    char title[120];
} titles_;

#define titles_1 titles_

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__25000 = 25000;
static integer c__5 = 5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine readxyz  --  input of Cartesian coordinates  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "readxyz" gets a set of Cartesian coordinates from */
/*     an external disk file */


/* Subroutine */ int readxyz_(integer *ixyz)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 READXYZ  --  Unable to Find the Cartes"
	    "ian\002,\002 Coordinates File\002)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(/,\002 READXYZ  --  The Maximum of\002,i8,\002 "
	    "Atoms\002,\002 has been Exceeded\002)";
    static char fmt_40[] = "(a120)";
    static char fmt_70[] = "(/,\002 READXYZ  --  Error in Coordinate File at"
	    " Atom\002,i6)";
    static char fmt_90[] = "(/,\002 READXYZ  --  Atom Labels not Sequential"
	    ",\002,\002 Attempting to Renumber\002)";
    static char fmt_100[] = "(/,\002 READXYZ  --  Check Connection of Ato"
	    "m\002,i6,\002 to Atom\002,i6)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3, i__4;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_rew(alist *), s_wsfe(cilist *), e_wsfe(void), 
	    s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void),
	     s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_clos(cllist *);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen), nexttext_(char *, ftnlen);
    static integer i__, j, k, m, last, size, list[25000], next;
    static logical quit;
    extern /* Subroutine */ int sort_(integer *, integer *), fatal_(void);
    static logical clash;
    static integer first;
    static logical exist, opened;
    static char record[120], string[120];
    extern /* Subroutine */ int chkxyz_(logical *);
    static logical reorder;
    extern /* Subroutine */ int getword_(char *, char *, integer *, ftnlen, 
	    ftnlen), gettext_(char *, char *, integer *, ftnlen, ftnlen), 
	    version_(char *, char *, ftnlen, ftnlen);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___7 = { 1, 0, 1, fmt_20, 0 };
    static icilist io___12 = { 1, string, 1, 0, 120, 1 };
    static cilist io___15 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___17 = { 1, 0, 1, fmt_40, 0 };
    static icilist io___18 = { 1, record, 1, 0, 120, 1 };
    static icilist io___19 = { 1, string, 1, 0, 120, 1 };
    static cilist io___20 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_100, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  titles.i  --  title for the current molecular system  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     ltitle   length in characters of the nonblank title string */
/*     title    title used to describe the current structure */




/*     initialize the total number of atoms in the system */

    atoms_1.n = 0;

/*     open the input file if it has not already been done */

    ioin__1.inerr = 0;
    ioin__1.inunit = *ixyz;
    ioin__1.infile = 0;
    ioin__1.inex = 0;
    ioin__1.inopen = &opened;
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
    if (! opened) {
/* Writing concatenation */
	i__1[0] = files_1.leng, a__1[0] = files_1.filename;
	i__1[1] = 4, a__1[1] = ".xyz";
	s_cat(xyzfile, a__1, i__1, &c__2, (ftnlen)120);
	version_(xyzfile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = xyzfile;
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
	    o__1.oerr = 0;
	    o__1.ounit = *ixyz;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = xyzfile;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	    al__1.aerr = 0;
	    al__1.aunit = *ixyz;
	    f_rew(&al__1);
	} else {
	    io___4.ciunit = iounit_1.iout;
	    s_wsfe(&io___4);
	    e_wsfe();
	    fatal_();
	}
    }

/*     read first line and return if already at end of file */

    quit = FALSE_;
    inform_1.abort = TRUE_;
    size = 0;
    while(size == 0) {
	io___7.ciunit = *ixyz;
	i__2 = s_rsfe(&io___7);
	if (i__2 != 0) {
	    goto L60;
	}
	i__2 = do_fio(&c__1, record, (ftnlen)120);
	if (i__2 != 0) {
	    goto L60;
	}
	i__2 = e_rsfe();
	if (i__2 != 0) {
	    goto L60;
	}
	size = trimtext_(record, (ftnlen)120);
    }
    inform_1.abort = FALSE_;
    quit = TRUE_;

/*     parse the title line to get the number of atoms */

    i__ = 0;
    next = 1;
    gettext_(record, string, &next, (ftnlen)120, (ftnlen)120);
    i__2 = s_rsli(&io___12);
    if (i__2 != 0) {
	goto L60;
    }
    i__2 = do_lio(&c__3, &c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
    if (i__2 != 0) {
	goto L60;
    }
    i__2 = e_rsli();
    if (i__2 != 0) {
	goto L60;
    }

/*     extract the title and determine its length */

    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
    first = nexttext_(string, (ftnlen)120);
    last = trimtext_(string, (ftnlen)120);
    if (last == 0) {
	s_copy(titles_1.title, " ", (ftnlen)120, (ftnlen)1);
	titles_1.ltitle = 0;
    } else {
	s_copy(titles_1.title, string + (first - 1), (ftnlen)120, last - (
		first - 1));
	titles_1.ltitle = trimtext_(titles_1.title, (ftnlen)120);
    }

/*     check for too many total atoms in the file */

    if (atoms_1.n > 25000) {
	io___15.ciunit = iounit_1.iout;
	s_wsfe(&io___15);
	do_fio(&c__1, (char *)&c__25000, (ftnlen)sizeof(integer));
	e_wsfe();
	fatal_();
    }

/*     initialize coordinates and connectivities for each atom */

    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	atoms_1.x[i__ - 1] = 0.;
	atoms_1.y[i__ - 1] = 0.;
	atoms_1.z__[i__ - 1] = 0.;
	atoms_1.type__[i__ - 1] = 0;
	for (j = 1; j <= 8; ++j) {
	    i12_ref(j, i__) = 0;
	}
    }

/*     read the coordinates and connectivities for each atom */

    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	next = 1;
	size = 0;
	while(size == 0) {
	    io___17.ciunit = *ixyz;
	    i__3 = s_rsfe(&io___17);
	    if (i__3 != 0) {
		goto L60;
	    }
	    i__3 = do_fio(&c__1, record, (ftnlen)120);
	    if (i__3 != 0) {
		goto L60;
	    }
	    i__3 = e_rsfe();
	    if (i__3 != 0) {
		goto L60;
	    }
	    size = trimtext_(record, (ftnlen)120);
	}
	i__3 = s_rsli(&io___18);
	if (i__3 != 0) {
	    goto L60;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&atmtyp_1.tag[i__ - 1], (ftnlen)
		sizeof(integer));
	if (i__3 != 0) {
	    goto L60;
	}
	i__3 = e_rsli();
	if (i__3 != 0) {
	    goto L60;
	}
	getword_(record, name___ref(0, i__), &next, (ftnlen)120, (ftnlen)3);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	i__3 = s_rsli(&io___19);
	if (i__3 != 0) {
	    goto L50;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&atoms_1.x[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__3 != 0) {
	    goto L50;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&atoms_1.y[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__3 != 0) {
	    goto L50;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&atoms_1.z__[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__3 != 0) {
	    goto L50;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&atoms_1.type__[i__ - 1], (ftnlen)
		sizeof(integer));
	if (i__3 != 0) {
	    goto L50;
	}
	for (j = 1; j <= 8; ++j) {
	    i__3 = do_lio(&c__3, &c__1, (char *)&i12_ref(j, i__), (ftnlen)
		    sizeof(integer));
	    if (i__3 != 0) {
		goto L50;
	    }
	}
	i__3 = e_rsli();
	if (i__3 != 0) {
	    goto L50;
	}
L50:
	;
    }
    quit = FALSE_;
L60:
    if (! opened) {
	cl__1.cerr = 0;
	cl__1.cunit = *ixyz;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     an error occurred in reading the coordinate file */

    if (quit) {
	io___20.ciunit = iounit_1.iout;
	s_wsfe(&io___20);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
	fatal_();
    }

/*     for each atom, count and sort its attached atoms */

    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	couple_1.n12[i__ - 1] = 0;
	for (j = 8; j >= 1; --j) {
	    if (i12_ref(j, i__) != 0) {
		couple_1.n12[i__ - 1] = j;
		goto L80;
	    }
	}
L80:
	sort_(&couple_1.n12[i__ - 1], &i12_ref(1, i__));
    }

/*     check for scrambled atom order and attempt to renumber */

    reorder = FALSE_;
    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	list[atmtyp_1.tag[i__ - 1] - 1] = i__;
	if (atmtyp_1.tag[i__ - 1] != i__) {
	    reorder = TRUE_;
	}
    }
    if (reorder) {
	io___23.ciunit = iounit_1.iout;
	s_wsfe(&io___23);
	e_wsfe();
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    atmtyp_1.tag[i__ - 1] = i__;
	    i__3 = couple_1.n12[i__ - 1];
	    for (j = 1; j <= i__3; ++j) {
		i12_ref(j, i__) = list[i12_ref(j, i__) - 1];
	    }
	    sort_(&couple_1.n12[i__ - 1], &i12_ref(1, i__));
	}
    }

/*     check for atom pairs with identical coordinates */

    clash = FALSE_;
    chkxyz_(&clash);

/*     make sure that all connectivities are bidirectional */

    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__3 = couple_1.n12[i__ - 1];
	for (j = 1; j <= i__3; ++j) {
	    k = i12_ref(j, i__);
	    i__4 = couple_1.n12[k - 1];
	    for (m = 1; m <= i__4; ++m) {
		if (i12_ref(m, k) == i__) {
		    goto L110;
		}
	    }
	    io___27.ciunit = iounit_1.iout;
	    s_wsfe(&io___27);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    e_wsfe();
	    fatal_();
L110:
	    ;
	}
    }
    return 0;
} /* readxyz_ */

#undef name___ref
#undef i12_ref


