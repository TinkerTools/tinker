/* readint.f -- translated by f2c (version 20050501).
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
    integer ltitle;
    char title[120];
} titles_;

#define titles_1 titles_

struct {
    integer nadd, iadd[50000]	/* was [2][25000] */, ndel, idel[50000]	/* 
	    was [2][25000] */;
} zclose_;

#define zclose_1 zclose_

struct {
    doublereal zbond[25000], zang[25000], ztors[25000];
    integer iz[100000]	/* was [4][25000] */;
} zcoord_;

#define zcoord_1 zcoord_

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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine readint  --  input of internal coordinates  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "readint" gets a set of Z-matrix internal coordinates */
/*     from an external file */


/* Subroutine */ int readint_(integer *izmt)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 READINT  --  Unable to Find the Interna"
	    "l\002,\002 Coordinates File\002)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(\002 READINT  --  The Maximum of\002,i8,\002 At"
	    "oms\002,\002 has been Exceeded\002)";
    static char fmt_40[] = "(a120)";
    static char fmt_70[] = "(\002 READZ  --  Error in Z-Matrix File at Ato"
	    "m\002,i6)";
    static char fmt_80[] = "()";
    static char fmt_90[] = "(a120)";
    static char fmt_110[] = "(a120)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3;
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
    static integer i__, j, last, size, next;
    static logical quit;
    extern /* Subroutine */ int fatal_(void);
    static integer first;
    static logical exist, opened;
    static char record[120], string[120], intfile[120];
    extern /* Subroutine */ int getword_(char *, char *, integer *, ftnlen, 
	    ftnlen), gettext_(char *, char *, integer *, ftnlen, ftnlen), 
	    version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___7 = { 1, 0, 1, fmt_20, 0 };
    static icilist io___12 = { 1, string, 1, 0, 120, 1 };
    static cilist io___15 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___17 = { 1, 0, 1, fmt_40, 0 };
    static icilist io___18 = { 1, record, 1, 0, 120, 1 };
    static icilist io___19 = { 1, string, 1, 0, 120, 1 };
    static cilist io___20 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___21 = { 1, 0, 1, fmt_80, 0 };
    static cilist io___22 = { 1, 0, 1, fmt_90, 0 };
    static icilist io___23 = { 1, record, 1, 0, 120, 1 };
    static cilist io___24 = { 1, 0, 1, fmt_110, 0 };
    static icilist io___25 = { 1, record, 1, 0, 120, 1 };



#define iz_ref(a_1,a_2) zcoord_1.iz[(a_2)*4 + a_1 - 5]
#define iadd_ref(a_1,a_2) zclose_1.iadd[(a_2)*2 + a_1 - 3]
#define idel_ref(a_1,a_2) zclose_1.idel[(a_2)*2 + a_1 - 3]
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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  zclose.i  --  ring openings and closures for Z-matrix  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     nadd   number of added bonds between Z-matrix atoms */
/*     iadd   numbers of the atom pairs defining added bonds */
/*     ndel   number of bonds between Z-matrix bonds to delete */
/*     idel   numbers of the atom pairs defining deleted bonds */




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




/*     initialize the total number of atoms in the system */

    atoms_1.n = 0;

/*     open the input file if it has not already been done */

    ioin__1.inerr = 0;
    ioin__1.inunit = *izmt;
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
	i__1[1] = 4, a__1[1] = ".int";
	s_cat(intfile, a__1, i__1, &c__2, (ftnlen)120);
	o__1.oerr = 0;
	o__1.ounit = *izmt;
	o__1.ofnmlen = 120;
	o__1.ofnm = intfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = *izmt;
	f_rew(&al__1);
	version_(intfile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = intfile;
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
	    o__1.ounit = *izmt;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = intfile;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	    al__1.aerr = 0;
	    al__1.aunit = *izmt;
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
	io___7.ciunit = *izmt;
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
	zcoord_1.zbond[i__ - 1] = 0.;
	zcoord_1.zang[i__ - 1] = 0.;
	zcoord_1.ztors[i__ - 1] = 0.;
	atoms_1.type__[i__ - 1] = 0;
	for (j = 1; j <= 4; ++j) {
	    iz_ref(j, i__) = 0;
	}
    }

/*     read the coordinates and connectivities for each atom */

    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	next = 1;
	size = 0;
	while(size == 0) {
	    io___17.ciunit = *izmt;
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
	i__3 = do_lio(&c__3, &c__1, (char *)&atoms_1.type__[i__ - 1], (ftnlen)
		sizeof(integer));
	if (i__3 != 0) {
	    goto L50;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&iz_ref(1, i__), (ftnlen)sizeof(
		integer));
	if (i__3 != 0) {
	    goto L50;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&zcoord_1.zbond[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__3 != 0) {
	    goto L50;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&iz_ref(2, i__), (ftnlen)sizeof(
		integer));
	if (i__3 != 0) {
	    goto L50;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&zcoord_1.zang[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__3 != 0) {
	    goto L50;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&iz_ref(3, i__), (ftnlen)sizeof(
		integer));
	if (i__3 != 0) {
	    goto L50;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&zcoord_1.ztors[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__3 != 0) {
	    goto L50;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&iz_ref(4, i__), (ftnlen)sizeof(
		integer));
	if (i__3 != 0) {
	    goto L50;
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
	cl__1.cunit = *izmt;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     an error occurred in reading the Z-matrix coordinates */

    if (quit) {
	io___20.ciunit = iounit_1.iout;
	s_wsfe(&io___20);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
	fatal_();
    }

/*     read in any additional bonds to be added or deleted */

    zclose_1.nadd = 0;
    zclose_1.ndel = 0;
    io___21.ciunit = *izmt;
    i__2 = s_rsfe(&io___21);
    if (i__2 != 0) {
	goto L120;
    }
    i__2 = e_rsfe();
    if (i__2 != 0) {
	goto L120;
    }
    for (i__ = 1; i__ <= 25000; ++i__) {
	io___22.ciunit = *izmt;
	i__2 = s_rsfe(&io___22);
	if (i__2 != 0) {
	    goto L120;
	}
	i__2 = do_fio(&c__1, record, (ftnlen)120);
	if (i__2 != 0) {
	    goto L120;
	}
	i__2 = e_rsfe();
	if (i__2 != 0) {
	    goto L120;
	}
	i__2 = s_rsli(&io___23);
	if (i__2 != 0) {
	    goto L100;
	}
	for (j = 1; j <= 2; ++j) {
	    i__2 = do_lio(&c__3, &c__1, (char *)&iadd_ref(j, i__), (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L100;
	    }
	}
	i__2 = e_rsli();
	if (i__2 != 0) {
	    goto L100;
	}
	zclose_1.nadd = i__;
    }
L100:
    for (i__ = 1; i__ <= 25000; ++i__) {
	io___24.ciunit = *izmt;
	i__2 = s_rsfe(&io___24);
	if (i__2 != 0) {
	    goto L120;
	}
	i__2 = do_fio(&c__1, record, (ftnlen)120);
	if (i__2 != 0) {
	    goto L120;
	}
	i__2 = e_rsfe();
	if (i__2 != 0) {
	    goto L120;
	}
	i__2 = s_rsli(&io___25);
	if (i__2 != 0) {
	    goto L120;
	}
	for (j = 1; j <= 2; ++j) {
	    i__2 = do_lio(&c__3, &c__1, (char *)&idel_ref(j, i__), (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L120;
	    }
	}
	i__2 = e_rsli();
	if (i__2 != 0) {
	    goto L120;
	}
	zclose_1.ndel = i__;
    }
L120:
    return 0;
} /* readint_ */

#undef name___ref
#undef idel_ref
#undef iadd_ref
#undef iz_ref


