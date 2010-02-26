/* superpose.f -- translated by f2c (version 20050501).
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
    doublereal wfit[25000];
    integer nfit, ifit[50000]	/* was [2][25000] */;
} align_;

#define align_1 align_

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

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  program superpose  --  optimal coordinate superposition  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "superpose" takes pairs of structures and superimposes them */
/*     in the optimal least squares sense; it will attempt to match */
/*     all atom pairs or only those specified by the user */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Two Options are Available :  (1) Fit ato"
	    "ms\002,\002 \"M\" through \"N\" from structure 1\002,/,\002 to t"
	    "he corresponding atoms of structure 2.\002,\002 Enter \"1,M,N\" "
	    "to use this option.\002,/,\002 If \"N\" is omitted, the fit uses"
	    " atoms 1\002,\002 through \"M\". If both \"M\" and \"N\" are\002"
	    ",/,\002 omitted, the fit uses all atoms; or (2)\002,\002 Individ"
	    "ual entry of atom range pairs\002,/,\002 to be used in the fitti"
	    "ng procedure.\002)";
    static char fmt_30[] = "(/,\002 Enter an Option (either 1,M,N or 2\002"
	    ",\002 [<CR>=1,0,0]) :  \002,$)";
    static char fmt_40[] = "(a120)";
    static char fmt_60[] = "(/,\002 SUPERPOSE  --  The Molecules contai"
	    "n\002,\002 Different Numbers of Atoms\002)";
    static char fmt_70[] = "(/,\002 Include Hydrogen Atoms in the Fitting"
	    "\002,\002 [Y] :  \002,$)";
    static char fmt_80[] = "(a120)";
    static char fmt_90[] = "(/,\002 On successive lines below, enter atom"
	    "\002,\002 pairs or pairs of atom ranges to use\002,/,\002 during"
	    " fitting. Entering \"4,7\" will fit\002,\002 atom 4 of structure"
	    " 1 to atom 7 of\002,/,\002 structure 2, while the entry \"4,7,9,"
	    "12\"\002,\002 will match atoms 4 through 7 from\002,/,\002 struc"
	    "ture 1 with atoms 9 through 12 of\002,\002 structure 2. Hit <RET"
	    "> to end entry\002,/,\002 of the list of pairs.\002)";
    static char fmt_100[] = "(/,\002 Enter a Pair of Atoms or Ranges :  \002"
	    ",$)";
    static char fmt_110[] = "(a120)";
    static char fmt_140[] = "(/,\002 Use Mass- or Unit-Weighted Coordinate"
	    "s\002,\002 (M or [U]) :  \002,$)";
    static char fmt_150[] = "(a120)";
    static char fmt_160[] = "(/,\002 Write Best-Fit Coordinates of 2nd Molec"
	    "ule\002,\002 [N] :  \002,$)";
    static char fmt_170[] = "(a120)";
    static char fmt_190[] = "(/,\002 Cutoff Value for Listing RMS Deviation"
	    "s\002,\002 [0.0] :  \002,$)";
    static char fmt_200[] = "(f20.0)";
    static char fmt_210[] = "(/,\002 Structure File 1 :  \002,a)";
    static char fmt_220[] = "(/,\002 Structure File 2 :  \002,a)";
    static char fmt_230[] = "(/,\002 File 1 Frame :\002,i6,13x,\002File 2 Fr"
	    "ame :\002,i6)";
    static char fmt_240[] = "(/,\002 Summary of Results from Structural\002"
	    ",\002 Superposition :\002)";
    static char fmt_250[] = "(/,\002 Root Mean Square Distance :\002,11x,f15"
	    ".6,2x,2i7)";
    static char fmt_260[] = "(/,\002   Atom in the\002,9x,\002Atom in the"
	    "\002,12x,\002Distance\002,10x,\002Weight\002/,\002 First Structu"
	    "re\002,5x,\002Second Structure\002,8x,\002Separated\002,10x,\002"
	    "in Fit\002/)";
    static char fmt_270[] = "(5x,i5,\002-\002,a3,11x,i5,\002-\002,a3,7x,f13."
	    "6,4x,f12.4)";
    static char fmt_280[] = "(/,\002 Root Mean Square Distance :\002,11x,f15"
	    ".6)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2];
    doublereal d__1, d__2, d__3;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void), f_open(olist *), 
	    f_rew(alist *), s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_clos(cllist *);

    /* Local variables */
    extern integer freeunit_(void);
    static doublereal rmsvalue;
    extern integer trimtext_(char *, ftnlen);
    static integer i__, i1, i2, n1, n2;
    static doublereal x1[25000], x2[25000], y1[25000], y2[25000], z1[25000], 
	    z2[25000];
    static logical same;
    static doublereal dist;
    static logical skip;
    static integer next, stop;
    static char file1[120], name1[3*25000], name2[3*25000];
    static integer ixyz;
    static char file2[120];
    static integer leng1, leng2;
    static doublereal mass1[25000], mass2[25000];
    extern /* Subroutine */ int field_(void);
    static integer delta;
    extern /* Subroutine */ int final_(void);
    static integer range[4];
    extern /* Subroutine */ int katom_(void);
    static logical exist;
    static integer start;
    static logical query;
    static integer ifile1, ifile2, frame1, frame2;
    static logical header;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static doublereal cutoff;
    static char answer[1];
    static integer option;
    static char string[120];
    extern /* Subroutine */ int suffix_(char *, char *, ftnlen, ftnlen), 
	    impose_(integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    , getxyz_(void);
    static integer atomic1[25000], atomic2[25000];
    extern /* Subroutine */ int prtxyz_(integer *), initial_(void), nextarg_(
	    char *, logical *, ftnlen), gettext_(char *, char *, integer *, 
	    ftnlen, ftnlen), version_(char *, char *, ftnlen, ftnlen), 
	    readxyz_(integer *);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static icilist io___21 = { 1, string, 1, 0, 120, 1 };
    static icilist io___22 = { 1, string, 1, 0, 120, 1 };
    static icilist io___23 = { 1, string, 1, 0, 120, 1 };
    static cilist io___24 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_40, 0 };
    static icilist io___28 = { 1, record, 1, 0, 120, 1 };
    static cilist io___29 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_110, 0 };
    static icilist io___38 = { 1, record, 1, 0, 120, 1 };
    static cilist io___40 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_170, 0 };
    static icilist io___45 = { 1, string, 1, 0, 120, 1 };
    static cilist io___46 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___48 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___49 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___64 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___69 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___71 = { 0, 0, 0, fmt_280, 0 };



#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define ifit_ref(a_1,a_2) align_1.ifit[(a_2)*2 + a_1 - 3]
#define name1_ref(a_0,a_1) &name1[(a_1)*3 + a_0 - 3]
#define name2_ref(a_0,a_1) &name2[(a_1)*3 + a_0 - 3]



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  align.i  --  information for superposition of structures  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     wfit    weights assigned to atom pairs during superposition */
/*     nfit    number of atoms to use in superimposing two structures */
/*     ifit    atom numbers of pairs of atoms to be superimposed */




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




/*     get atom names and masses for the first structure type */

    initial_();
    getxyz_();
    field_();
    katom_();
    s_copy(file1, files_1.filename, (ftnlen)120, (ftnlen)120);
    leng1 = trimtext_(file1, (ftnlen)120);
    n1 = atoms_1.n;
    i__1 = n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(name1_ref(0, i__), name___ref(0, i__), (ftnlen)3, (ftnlen)3);
	atomic1[i__ - 1] = atmtyp_1.atomic[i__ - 1];
	mass1[i__ - 1] = atmtyp_1.mass[i__ - 1];
    }

/*     get atom names and masses for the second structure type */

    getxyz_();
    field_();
    katom_();
    s_copy(file2, files_1.filename, (ftnlen)120, (ftnlen)120);
    leng2 = trimtext_(file2, (ftnlen)120);
    n2 = atoms_1.n;
    i__1 = n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(name2_ref(0, i__), name___ref(0, i__), (ftnlen)3, (ftnlen)3);
	atomic2[i__ - 1] = atmtyp_1.atomic[i__ - 1];
	mass2[i__ - 1] = atmtyp_1.mass[i__ - 1];
    }

/*     get atom pairs to be superimposed from command line */

    option = 0;
    start = 0;
    stop = 0;
    *(unsigned char *)answer = ' ';
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	query = FALSE_;
	i__1 = s_rsli(&io___21);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&option, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (option == 1) {
	    nextarg_(string, &exist, (ftnlen)120);
	    if (exist) {
		*(unsigned char *)answer = *(unsigned char *)string;
		i__1 = s_rsli(&io___22);
		if (i__1 != 0) {
		    goto L10;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&start, (ftnlen)sizeof(
			integer));
		if (i__1 != 0) {
		    goto L10;
		}
		i__1 = e_rsli();
		if (i__1 != 0) {
		    goto L10;
		}
		*(unsigned char *)answer = ' ';
		nextarg_(string, &exist, (ftnlen)120);
		if (exist) {
		    *(unsigned char *)answer = *(unsigned char *)string;
		    i__1 = s_rsli(&io___23);
		    if (i__1 != 0) {
			goto L10;
		    }
		    i__1 = do_lio(&c__3, &c__1, (char *)&stop, (ftnlen)sizeof(
			    integer));
		    if (i__1 != 0) {
			goto L10;
		    }
		    i__1 = e_rsli();
		    if (i__1 != 0) {
			goto L10;
		    }
		    *(unsigned char *)answer = ' ';
		}
	    }
	}
    }
L10:

/*     ask the user which pairs of atoms are to be superimposed */

    if (query) {
	io___24.ciunit = iounit_1.iout;
	s_wsfe(&io___24);
	e_wsfe();
	io___25.ciunit = iounit_1.iout;
	s_wsfe(&io___25);
	e_wsfe();
	io___26.ciunit = iounit_1.input;
	s_rsfe(&io___26);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___28);
	if (i__1 != 0) {
	    goto L50;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&option, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L50;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&start, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L50;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&stop, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L50;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L50;
	}
L50:
	if (option < 1 || option > 2) {
	    option = 1;
	    start = 0;
	    stop = 0;
	}
    }

/*     warning if structures have different numbers of atoms */

    if (option == 1) {
	if (n1 != n2 && start == 0) {
	    io___29.ciunit = iounit_1.iout;
	    s_wsfe(&io___29);
	    e_wsfe();
	}
    }

/*     setup automatic superposition with option to omit hydrogens */

    if (option == 1) {
	if (*(unsigned char *)answer == ' ') {
	    nextarg_(answer, &exist, (ftnlen)1);
	} else {
	    exist = TRUE_;
	}
	if (! exist) {
	    io___30.ciunit = iounit_1.iout;
	    s_wsfe(&io___30);
	    e_wsfe();
	    io___31.ciunit = iounit_1.input;
	    s_rsfe(&io___31);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    next = 1;
	    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
	}
	upcase_(answer, (ftnlen)1);
	if (start == 0 && stop == 0) {
	    start = 1;
	    stop = min(n1,n2);
	} else if (start != 0 && stop == 0) {
/* Computing MIN */
	    i__1 = min(n1,n2);
	    stop = min(i__1,start);
	    start = 1;
	} else if (start != 0 && stop != 0) {
	    start = max(1,start);
/* Computing MIN */
	    i__1 = min(n1,n2);
	    stop = min(i__1,stop);
	}
	align_1.nfit = 0;
	i__1 = stop;
	for (i__ = start; i__ <= i__1; ++i__) {
	    skip = FALSE_;
	    if (*(unsigned char *)answer == 'N') {
		if (atomic1[i__ - 1] <= 1 || atomic2[i__ - 1] <= 1) {
		    skip = TRUE_;
		}
	    }
	    if (! skip) {
		++align_1.nfit;
		ifit_ref(1, align_1.nfit) = i__;
		ifit_ref(2, align_1.nfit) = i__;
	    }
	}
    }

/*     manual input of the pairs of atom ranges to superimpose */

    if (option == 2) {
	io___34.ciunit = iounit_1.iout;
	s_wsfe(&io___34);
	e_wsfe();
	align_1.nfit = 0;
	while(TRUE_) {
	    for (i__ = 1; i__ <= 4; ++i__) {
		range[i__ - 1] = 0;
	    }
	    io___36.ciunit = iounit_1.iout;
	    s_wsfe(&io___36);
	    e_wsfe();
	    io___37.ciunit = iounit_1.input;
	    s_rsfe(&io___37);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    i__1 = s_rsli(&io___38);
	    if (i__1 != 0) {
		goto L120;
	    }
	    for (i__ = 1; i__ <= 4; ++i__) {
		i__1 = do_lio(&c__3, &c__1, (char *)&range[i__ - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L120;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L120;
	    }
L120:
	    if (range[0] == 0) {
		goto L130;
	    } else if (range[1] == 0) {
		++align_1.nfit;
		ifit_ref(1, align_1.nfit) = range[0];
		ifit_ref(2, align_1.nfit) = range[0];
	    } else if (range[2] == 0) {
		++align_1.nfit;
		ifit_ref(1, align_1.nfit) = range[0];
		ifit_ref(2, align_1.nfit) = range[1];
	    } else {
		delta = range[2] - range[0];
		i__1 = range[1];
		for (i__ = range[0]; i__ <= i__1; ++i__) {
		    ++align_1.nfit;
		    ifit_ref(1, align_1.nfit) = i__;
		    ifit_ref(2, align_1.nfit) = i__ + delta;
		}
	    }
	}
L130:
	;
    }

/*     decide on the weighting to use for the coordinates */

    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___40.ciunit = iounit_1.iout;
	s_wsfe(&io___40);
	e_wsfe();
	io___41.ciunit = iounit_1.input;
	s_rsfe(&io___41);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'M') {
	i__1 = align_1.nfit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    align_1.wfit[i__ - 1] = (mass1[ifit_ref(1, i__) - 1] + mass2[
		    ifit_ref(2, i__) - 1]) * .5;
	}
    } else {
	i__1 = align_1.nfit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    align_1.wfit[i__ - 1] = 1.;
	}
    }

/*     decide whether to write the best fit set of coordinates */

    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___42.ciunit = iounit_1.iout;
	s_wsfe(&io___42);
	e_wsfe();
	io___43.ciunit = iounit_1.input;
	s_rsfe(&io___43);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);

/*     chose cutoff value for output of atom pair deviations */

    cutoff = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___45);
	if (i__1 != 0) {
	    goto L180;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&cutoff, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L180;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L180;
	}
    }
L180:
    if (cutoff < 0.) {
	io___46.ciunit = iounit_1.iout;
	s_wsfe(&io___46);
	e_wsfe();
	io___47.ciunit = iounit_1.input;
	s_rsfe(&io___47);
	do_fio(&c__1, (char *)&cutoff, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }

/*     information about structures to be superimposed */

    io___48.ciunit = iounit_1.iout;
    s_wsfe(&io___48);
    do_fio(&c__1, file1, leng1);
    e_wsfe();
    io___49.ciunit = iounit_1.iout;
    s_wsfe(&io___49);
    do_fio(&c__1, file2, leng2);
    e_wsfe();

/*     reopen the coordinate files with structures to superimpose */

    ifile1 = freeunit_();
    suffix_(file1, "xyz", (ftnlen)120, (ftnlen)3);
    version_(file1, "old", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = ifile1;
    o__1.ofnmlen = 120;
    o__1.ofnm = file1;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    al__1.aerr = 0;
    al__1.aunit = ifile1;
    f_rew(&al__1);
    suffix_(file2, "xyz", (ftnlen)120, (ftnlen)3);
    version_(file2, "old", (ftnlen)120, (ftnlen)3);
    if (s_cmp(file1, file2, (ftnlen)120, (ftnlen)120) == 0) {
	same = TRUE_;
	ifile2 = ifile1;
    } else {
	same = FALSE_;
	ifile2 = freeunit_();
	o__1.oerr = 0;
	o__1.ounit = ifile2;
	o__1.ofnmlen = 120;
	o__1.ofnm = file2;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = ifile2;
	f_rew(&al__1);
    }

/*     read initial structure set from the first coordinate file */

    frame1 = 1;
    readxyz_(&ifile1);
    n1 = atoms_1.n;
    i__1 = n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x1[i__ - 1] = atoms_1.x[i__ - 1];
	y1[i__ - 1] = atoms_1.y[i__ - 1];
	z1[i__ - 1] = atoms_1.z__[i__ - 1];
    }

/*     read initial structure set from the second coordinate file */

    frame2 = 1;
    if (same) {
	frame2 = 2;
    }
    readxyz_(&ifile2);
    n2 = atoms_1.n;
    i__1 = n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x2[i__ - 1] = atoms_1.x[i__ - 1];
	y2[i__ - 1] = atoms_1.y[i__ - 1];
	z2[i__ - 1] = atoms_1.z__[i__ - 1];
    }
    if (inform_1.abort) {
	inform_1.abort = FALSE_;
	frame2 = 1;
	n2 = n1;
	i__1 = n2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x2[i__ - 1] = x1[i__ - 1];
	    y2[i__ - 1] = y1[i__ - 1];
	    z2[i__ - 1] = z1[i__ - 1];
	}
    }

/*     perform the superposition of a structure pair */

    while(! inform_1.abort) {
	io___61.ciunit = iounit_1.iout;
	s_wsfe(&io___61);
	do_fio(&c__1, (char *)&frame1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&frame2, (ftnlen)sizeof(integer));
	e_wsfe();
	io___62.ciunit = iounit_1.iout;
	s_wsfe(&io___62);
	e_wsfe();
	inform_1.verbose = TRUE_;
	impose_(&n1, x1, y1, z1, &n2, x2, y2, z2, &rmsvalue);
	io___64.ciunit = iounit_1.iout;
	s_wsfe(&io___64);
	do_fio(&c__1, (char *)&rmsvalue, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&frame1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&frame2, (ftnlen)sizeof(integer));
	e_wsfe();

/*     write out the results of the superposition */

	header = TRUE_;
	i__1 = align_1.nfit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i1 = ifit_ref(1, i__);
	    i2 = ifit_ref(2, i__);
/* Computing 2nd power */
	    d__1 = x1[i1 - 1] - x2[i2 - 1];
/* Computing 2nd power */
	    d__2 = y1[i1 - 1] - y2[i2 - 1];
/* Computing 2nd power */
	    d__3 = z1[i1 - 1] - z2[i2 - 1];
	    dist = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	    if (dist >= cutoff) {
		if (header) {
		    header = FALSE_;
		    io___69.ciunit = iounit_1.iout;
		    s_wsfe(&io___69);
		    e_wsfe();
		}
		io___70.ciunit = iounit_1.iout;
		s_wsfe(&io___70);
		do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(integer));
		do_fio(&c__1, name1_ref(0, i1), (ftnlen)3);
		do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(integer));
		do_fio(&c__1, name2_ref(0, i2), (ftnlen)3);
		do_fio(&c__1, (char *)&dist, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&align_1.wfit[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
	if (! header) {
	    io___71.ciunit = iounit_1.iout;
	    s_wsfe(&io___71);
	    do_fio(&c__1, (char *)&rmsvalue, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}

/*     create output file for superimposed second structure */

	if (*(unsigned char *)answer == 'Y') {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		atoms_1.x[i__ - 1] = x2[i__ - 1];
		atoms_1.y[i__ - 1] = y2[i__ - 1];
		atoms_1.z__[i__ - 1] = z2[i__ - 1];
	    }
	    ixyz = freeunit_();
/* Writing concatenation */
	    i__2[0] = files_1.leng, a__1[0] = file2;
	    i__2[1] = 4, a__1[1] = ".xyz";
	    s_cat(xyzfile, a__1, i__2, &c__2, (ftnlen)120);
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

/*     attempt to get next structure pair from coordinate files */

	++frame2;
	readxyz_(&ifile2);
	n2 = atoms_1.n;
	i__1 = n2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x2[i__ - 1] = atoms_1.x[i__ - 1];
	    y2[i__ - 1] = atoms_1.y[i__ - 1];
	    z2[i__ - 1] = atoms_1.z__[i__ - 1];
	}
	if (inform_1.abort) {
	    inform_1.abort = FALSE_;
	    if (same) {
		al__1.aerr = 0;
		al__1.aunit = ifile1;
		f_rew(&al__1);
		i__1 = frame1;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    readxyz_(&ifile1);
		}
	    }
	    ++frame1;
	    readxyz_(&ifile1);
	    n1 = atoms_1.n;
	    i__1 = n1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		x1[i__ - 1] = atoms_1.x[i__ - 1];
		y1[i__ - 1] = atoms_1.y[i__ - 1];
		z1[i__ - 1] = atoms_1.z__[i__ - 1];
	    }
	    if (! inform_1.abort) {
		frame2 = frame1 + 1;
		if (! same) {
		    frame2 = 1;
		    al__1.aerr = 0;
		    al__1.aunit = ifile2;
		    f_rew(&al__1);
		}
		readxyz_(&ifile2);
		n2 = atoms_1.n;
		i__1 = n2;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    x2[i__ - 1] = atoms_1.x[i__ - 1];
		    y2[i__ - 1] = atoms_1.y[i__ - 1];
		    z2[i__ - 1] = atoms_1.z__[i__ - 1];
		}
	    }
	}
    }

/*     perform any final tasks before program exit */

    cl__1.cerr = 0;
    cl__1.cunit = ifile1;
    cl__1.csta = 0;
    f_clos(&cl__1);
    if (! same) {
	cl__1.cerr = 0;
	cl__1.cunit = ifile2;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    final_();
    return 0;
} /* MAIN__ */

#undef name2_ref
#undef name1_ref
#undef ifit_ref
#undef name___ref


/* Main program alias */ int superpose_ () { MAIN__ (); return 0; }
