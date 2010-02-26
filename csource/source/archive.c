/* archive.f -- translated by f2c (version 20050501).
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
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    integer ltitle;
    char title[120];
} titles_;

#define titles_1 titles_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  program archive  --  create or extract from an archive  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "archive" is a utility program for coordinate files which */
/*     concatenates multiple coordinate sets into a single archive */
/*     file, or extracts individual coordinate sets from an archive */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter Name of the Coordinate Archive Fil"
	    "e :  \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(/,\002 Create (C), Extract (E) from or Trim (T"
	    ")\002,\002 an Archive [C] :  \002,$)";
    static char fmt_40[] = "(a120)";
    static char fmt_60[] = "(/,\002 Numbers of First & Last File and Step"
	    "\002,\002 Increment :  \002,$)";
    static char fmt_70[] = "(a120)";
    static char fmt_90[] = "(/,\002 Enter Name of the Coordinate Archiv"
	    "e\002,\002 File :  \002,$)";
    static char fmt_100[] = "(a120)";
    static char fmt_110[] = "(/,\002 Numbers of the Atoms to be Removed : "
	    " \002,$)";
    static char fmt_120[] = "(a120)";
    static char fmt_150[] = "(/,\002 Numbers of First & Last File and Ste"
	    "p\002,\002 [<CR>=Exit] :  \002,$)";
    static char fmt_160[] = "(a120)";
    static char fmt_200[] = "(/,\002 Numbers of First & Last File and Ste"
	    "p\002,\002 [<CR>=Exit] :  \002,$)";
    static char fmt_210[] = "(a120)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2[3], i__3, i__4, i__5;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), f_open(olist *), s_rsli(
	    icilist *), do_lio(integer *, integer *, char *, ftnlen), e_rsli(
	    void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_inqu(inlist *), f_rew(alist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int basefile_(char *, ftnlen);
    static char basename[120];
    extern integer freeunit_(void);
    static integer i__, k;
    static char ext[7];
    static integer now, iarc;
    static char mode[7];
    static integer step, lext, list[25000], next, stop, ixyz, leng1, leng2, 
	    lengb;
    extern /* Subroutine */ int final_(void);
    static integer nlist;
    static logical exist;
    static integer start;
    static logical query;
    extern /* Subroutine */ int active_(void);
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char answer[1], string[120];
    extern /* Subroutine */ int suffix_(char *, char *, ftnlen, ftnlen), 
	    prtarc_(integer *);
    static char arcfile[120];
    extern /* Subroutine */ int initial_(void), numeral_(integer *, char *, 
	    integer *, ftnlen), nextarg_(char *, logical *, ftnlen), gettext_(
	    char *, char *, integer *, ftnlen, ftnlen), version_(char *, char 
	    *, ftnlen, ftnlen), readxyz_(integer *);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_40, 0 };
    static icilist io___19 = { 1, string, 1, 0, 120, 1 };
    static icilist io___20 = { 1, string, 1, 0, 120, 1 };
    static icilist io___21 = { 1, string, 1, 0, 120, 1 };
    static cilist io___22 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_70, 0 };
    static icilist io___24 = { 1, record, 1, 0, 120, 1 };
    static cilist io___31 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_120, 0 };
    static icilist io___37 = { 1, record, 1, 0, 120, 1 };
    static icilist io___41 = { 1, string, 1, 0, 120, 1 };
    static icilist io___42 = { 1, string, 1, 0, 120, 1 };
    static icilist io___43 = { 1, string, 1, 0, 120, 1 };
    static cilist io___44 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_160, 0 };
    static icilist io___46 = { 1, record, 1, 0, 120, 1 };
    static icilist io___47 = { 1, string, 1, 0, 120, 1 };
    static icilist io___48 = { 1, string, 1, 0, 120, 1 };
    static icilist io___49 = { 1, string, 1, 0, 120, 1 };
    static cilist io___50 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_210, 0 };
    static icilist io___52 = { 1, record, 1, 0, 120, 1 };




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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     get the name to use for the coordinate archive file */

    initial_();
    nextarg_(arcfile, &exist, (ftnlen)120);
    if (! exist) {
	io___3.ciunit = iounit_1.iout;
	s_wsfe(&io___3);
	e_wsfe();
	io___4.ciunit = iounit_1.input;
	s_rsfe(&io___4);
	do_fio(&c__1, arcfile, (ftnlen)120);
	e_rsfe();
    }

/*     decide whether to create or extract from an archive file */

    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___6.ciunit = iounit_1.iout;
	s_wsfe(&io___6);
	e_wsfe();
	io___7.ciunit = iounit_1.input;
	s_rsfe(&io___7);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'E') {
	s_copy(mode, "EXTRACT", (ftnlen)7, (ftnlen)7);
    } else if (*(unsigned char *)answer == 'T') {
	s_copy(mode, "TRIM", (ftnlen)7, (ftnlen)4);
    } else {
	s_copy(mode, "CREATE", (ftnlen)7, (ftnlen)6);
    }

/*     open the archive file to be processed */

    iarc = freeunit_();
    basefile_(arcfile, (ftnlen)120);
    s_copy(basename, arcfile, (ftnlen)120, (ftnlen)120);
    lengb = files_1.leng;
    suffix_(arcfile, "arc", (ftnlen)120, (ftnlen)3);
    if (s_cmp(mode, "CREATE", (ftnlen)7, (ftnlen)6) == 0) {
	version_(arcfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = iarc;
	o__1.ofnmlen = 120;
	o__1.ofnm = arcfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    } else {
	version_(arcfile, "old", (ftnlen)120, (ftnlen)3);
    }

/*     concatenate individual files into a single archive file */

    if (s_cmp(mode, "CREATE", (ftnlen)7, (ftnlen)6) == 0) {
	start = 0;
	stop = 0;
	step = 0;
	query = TRUE_;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___19);
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&start, (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L50;
	    }
	    query = FALSE_;
	}
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___20);
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&stop, (ftnlen)sizeof(integer)
		    );
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L50;
	    }
	}
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___21);
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&step, (ftnlen)sizeof(integer)
		    );
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L50;
	    }
	}
L50:
	if (query) {
	    io___22.ciunit = iounit_1.iout;
	    s_wsfe(&io___22);
	    e_wsfe();
	    io___23.ciunit = iounit_1.input;
	    s_rsfe(&io___23);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    i__1 = s_rsli(&io___24);
	    if (i__1 != 0) {
		goto L80;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&start, (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L80;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&stop, (ftnlen)sizeof(integer)
		    );
	    if (i__1 != 0) {
		goto L80;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&step, (ftnlen)sizeof(integer)
		    );
	    if (i__1 != 0) {
		goto L80;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L80;
	    }
L80:
	    ;
	}
	if (stop == 0) {
	    stop = start;
	}
	if (step == 0) {
	    step = 1;
	}

/*     cycle over the user specified coordinate files */

	i__ = start;
	while(i__ >= start && i__ <= stop) {
	    ixyz = freeunit_();
	    lext = 3;
	    numeral_(&i__, ext, &lext, (ftnlen)7);
/* Writing concatenation */
	    i__2[0] = lengb, a__1[0] = basename;
	    i__2[1] = 1, a__1[1] = ".";
	    i__2[2] = lext, a__1[2] = ext;
	    s_cat(xyzfile, a__1, i__2, &c__3, (ftnlen)120);
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
	    if (! exist && i__ < 100) {
		lext = 2;
		numeral_(&i__, ext, &lext, (ftnlen)7);
/* Writing concatenation */
		i__2[0] = lengb, a__1[0] = basename;
		i__2[1] = 1, a__1[1] = ".";
		i__2[2] = lext, a__1[2] = ext;
		s_cat(xyzfile, a__1, i__2, &c__3, (ftnlen)120);
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
	    }
	    if (! exist && i__ < 10) {
		lext = 1;
		numeral_(&i__, ext, &lext, (ftnlen)7);
/* Writing concatenation */
		i__2[0] = lengb, a__1[0] = basename;
		i__2[1] = 1, a__1[1] = ".";
		i__2[2] = lext, a__1[2] = ext;
		s_cat(xyzfile, a__1, i__2, &c__3, (ftnlen)120);
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
	    }
	    if (exist) {
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
		cl__1.cerr = 0;
		cl__1.cunit = ixyz;
		cl__1.csta = 0;
		f_clos(&cl__1);
		if (i__ == start) {
		    usage_1.nuse = atoms_1.n;
		    i__1 = usage_1.nuse;
		    for (k = 1; k <= i__1; ++k) {
			usage_1.use[k - 1] = TRUE_;
		    }
		}
		prtarc_(&iarc);
	    }
	    i__ += step;
	}
    }

/*     extract individual files from a concatenated archive file */

    if (s_cmp(mode, "EXTRACT", (ftnlen)7, (ftnlen)7) == 0 || s_cmp(mode, 
	    "TRIM", (ftnlen)7, (ftnlen)4) == 0) {
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = arcfile;
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
	while(! exist) {
	    io___31.ciunit = iounit_1.iout;
	    s_wsfe(&io___31);
	    e_wsfe();
	    io___32.ciunit = iounit_1.input;
	    s_rsfe(&io___32);
	    do_fio(&c__1, arcfile, (ftnlen)120);
	    e_rsfe();
	    basefile_(arcfile, (ftnlen)120);
	    s_copy(basename, arcfile, (ftnlen)120, (ftnlen)120);
	    lengb = files_1.leng;
	    suffix_(arcfile, "arc", (ftnlen)120, (ftnlen)3);
	    version_(arcfile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = arcfile;
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
	o__1.oerr = 0;
	o__1.ounit = iarc;
	o__1.ofnmlen = 120;
	o__1.ofnm = arcfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = iarc;
	f_rew(&al__1);
	readxyz_(&iarc);
	al__1.aerr = 0;
	al__1.aunit = iarc;
	f_rew(&al__1);

/*     decide whether atoms are to be removed from each frame */

	active_();
	if (s_cmp(mode, "EXTRACT", (ftnlen)7, (ftnlen)7) == 0) {
	    usage_1.nuse = atoms_1.n;
	    i__1 = usage_1.nuse;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		usage_1.use[i__ - 1] = TRUE_;
	    }
	} else if (s_cmp(mode, "TRIM", (ftnlen)7, (ftnlen)4) == 0 && 
		usage_1.nuse == atoms_1.n) {
	    nlist = 0;
	    for (i__ = 1; i__ <= 25000; ++i__) {
		list[i__ - 1] = 0;
	    }
	    io___35.ciunit = iounit_1.iout;
	    s_wsfe(&io___35);
	    e_wsfe();
	    io___36.ciunit = iounit_1.input;
	    s_rsfe(&io___36);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    i__1 = s_rsli(&io___37);
	    if (i__1 != 0) {
		goto L130;
	    }
	    for (i__ = 1; i__ <= 25000; ++i__) {
		i__1 = do_lio(&c__3, &c__1, (char *)&list[i__ - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L130;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L130;
	    }
L130:
	    i__ = 1;
	    while(list[i__ - 1] != 0) {
/* Computing MAX */
/* Computing MIN */
		i__4 = atoms_1.n, i__5 = list[i__ - 1];
		i__1 = -atoms_1.n, i__3 = min(i__4,i__5);
		list[i__ - 1] = max(i__1,i__3);
		if (list[i__ - 1] > 0) {
		    k = list[i__ - 1];
		    if (usage_1.use[k - 1]) {
			usage_1.use[k - 1] = FALSE_;
			--usage_1.nuse;
		    }
		    ++i__;
		} else {
/* Computing MAX */
/* Computing MIN */
		    i__4 = atoms_1.n, i__5 = list[i__];
		    i__1 = -atoms_1.n, i__3 = min(i__4,i__5);
		    list[i__] = max(i__1,i__3);
		    i__4 = (i__3 = list[i__], abs(i__3));
		    for (k = (i__1 = list[i__ - 1], abs(i__1)); k <= i__4; 
			    ++k) {
			if (usage_1.use[k - 1]) {
			    usage_1.use[k - 1] = FALSE_;
			    --usage_1.nuse;
			}
		    }
		    i__ += 2;
		}
	    }
	}

/*     get the initial and final coordinate frames to extract */

	now = 1;
	leng1 = 1;
	leng2 = files_1.leng;
	i__4 = files_1.leng;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    if (*(unsigned char *)&files_1.filename[i__ - 1] == '/') {
		leng1 = i__ + 1;
	    }
	    if (*(unsigned char *)&files_1.filename[i__ - 1] == ']') {
		leng1 = i__ + 1;
	    }
	    if (*(unsigned char *)&files_1.filename[i__ - 1] == ':') {
		leng1 = i__ + 1;
	    }
	}
	i__4 = leng1;
	for (i__ = files_1.leng; i__ >= i__4; --i__) {
	    if (*(unsigned char *)&files_1.filename[i__ - 1] == '.') {
		leng2 = i__ - 1;
	    }
	}
	files_1.leng = leng2 - leng1 + 1;
	s_copy(files_1.filename, files_1.filename + (leng1 - 1), files_1.leng,
		 leng2 - (leng1 - 1));
	start = 0;
	stop = 0;
	step = 0;
	query = TRUE_;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__4 = s_rsli(&io___41);
	    if (i__4 != 0) {
		goto L140;
	    }
	    i__4 = do_lio(&c__3, &c__1, (char *)&start, (ftnlen)sizeof(
		    integer));
	    if (i__4 != 0) {
		goto L140;
	    }
	    i__4 = e_rsli();
	    if (i__4 != 0) {
		goto L140;
	    }
	    query = FALSE_;
	}
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__4 = s_rsli(&io___42);
	    if (i__4 != 0) {
		goto L140;
	    }
	    i__4 = do_lio(&c__3, &c__1, (char *)&stop, (ftnlen)sizeof(integer)
		    );
	    if (i__4 != 0) {
		goto L140;
	    }
	    i__4 = e_rsli();
	    if (i__4 != 0) {
		goto L140;
	    }
	}
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__4 = s_rsli(&io___43);
	    if (i__4 != 0) {
		goto L140;
	    }
	    i__4 = do_lio(&c__3, &c__1, (char *)&step, (ftnlen)sizeof(integer)
		    );
	    if (i__4 != 0) {
		goto L140;
	    }
	    i__4 = e_rsli();
	    if (i__4 != 0) {
		goto L140;
	    }
	}
L140:
	if (query) {
	    io___44.ciunit = iounit_1.iout;
	    s_wsfe(&io___44);
	    e_wsfe();
	    io___45.ciunit = iounit_1.input;
	    s_rsfe(&io___45);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    i__4 = s_rsli(&io___46);
	    if (i__4 != 0) {
		goto L170;
	    }
	    i__4 = do_lio(&c__3, &c__1, (char *)&start, (ftnlen)sizeof(
		    integer));
	    if (i__4 != 0) {
		goto L170;
	    }
	    i__4 = do_lio(&c__3, &c__1, (char *)&stop, (ftnlen)sizeof(integer)
		    );
	    if (i__4 != 0) {
		goto L170;
	    }
	    i__4 = do_lio(&c__3, &c__1, (char *)&step, (ftnlen)sizeof(integer)
		    );
	    if (i__4 != 0) {
		goto L170;
	    }
	    i__4 = e_rsli();
	    if (i__4 != 0) {
		goto L170;
	    }
L170:
	    ;
	}
	if (stop == 0) {
	    stop = start;
	}
	if (step == 0) {
	    step = 1;
	}

/*     loop over the individual coordinate files to be extracted */

	while(start != 0) {
	    if (start <= now) {
		now = 1;
		al__1.aerr = 0;
		al__1.aunit = iarc;
		f_rew(&al__1);
	    }
	    i__4 = start - now;
	    for (k = 1; k <= i__4; ++k) {
		readxyz_(&iarc);
	    }
	    i__ = start;
	    if (s_cmp(mode, "EXTRACT", (ftnlen)7, (ftnlen)7) == 0) {
		while(i__ >= start && i__ <= stop) {
		    lext = 3;
		    numeral_(&i__, ext, &lext, (ftnlen)7);
		    readxyz_(&iarc);
		    if (inform_1.abort) {
			goto L180;
		    }
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
		    prtarc_(&ixyz);
		    cl__1.cerr = 0;
		    cl__1.cunit = ixyz;
		    cl__1.csta = 0;
		    f_clos(&cl__1);
		    i__ += step;
		    i__4 = step - 1;
		    for (k = 1; k <= i__4; ++k) {
			readxyz_(&iarc);
		    }
		}
	    } else if (s_cmp(mode, "TRIM", (ftnlen)7, (ftnlen)4) == 0) {
		ixyz = freeunit_();
		s_copy(xyzfile, basename, (ftnlen)120, (ftnlen)120);
		suffix_(xyzfile, "arc", (ftnlen)120, (ftnlen)3);
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
		while(i__ >= start && i__ <= stop) {
		    readxyz_(&iarc);
		    if (inform_1.abort) {
			goto L180;
		    }
		    prtarc_(&ixyz);
		    i__ += step;
		    i__4 = step - 1;
		    for (k = 1; k <= i__4; ++k) {
			readxyz_(&iarc);
		    }
		}
		cl__1.cerr = 0;
		cl__1.cunit = ixyz;
		cl__1.csta = 0;
		f_clos(&cl__1);
	    }
L180:
	    now = stop;
	    start = 0;
	    stop = 0;
	    step = 0;
	    query = TRUE_;
	    nextarg_(string, &exist, (ftnlen)120);
	    if (exist) {
		i__4 = s_rsli(&io___47);
		if (i__4 != 0) {
		    goto L190;
		}
		i__4 = do_lio(&c__3, &c__1, (char *)&start, (ftnlen)sizeof(
			integer));
		if (i__4 != 0) {
		    goto L190;
		}
		i__4 = e_rsli();
		if (i__4 != 0) {
		    goto L190;
		}
		query = FALSE_;
	    }
	    nextarg_(string, &exist, (ftnlen)120);
	    if (exist) {
		i__4 = s_rsli(&io___48);
		if (i__4 != 0) {
		    goto L190;
		}
		i__4 = do_lio(&c__3, &c__1, (char *)&stop, (ftnlen)sizeof(
			integer));
		if (i__4 != 0) {
		    goto L190;
		}
		i__4 = e_rsli();
		if (i__4 != 0) {
		    goto L190;
		}
	    }
	    nextarg_(string, &exist, (ftnlen)120);
	    if (exist) {
		i__4 = s_rsli(&io___49);
		if (i__4 != 0) {
		    goto L190;
		}
		i__4 = do_lio(&c__3, &c__1, (char *)&step, (ftnlen)sizeof(
			integer));
		if (i__4 != 0) {
		    goto L190;
		}
		i__4 = e_rsli();
		if (i__4 != 0) {
		    goto L190;
		}
	    }
L190:
	    if (query) {
		io___50.ciunit = iounit_1.iout;
		s_wsfe(&io___50);
		e_wsfe();
		io___51.ciunit = iounit_1.input;
		s_rsfe(&io___51);
		do_fio(&c__1, record, (ftnlen)120);
		e_rsfe();
		i__4 = s_rsli(&io___52);
		if (i__4 != 0) {
		    goto L220;
		}
		i__4 = do_lio(&c__3, &c__1, (char *)&start, (ftnlen)sizeof(
			integer));
		if (i__4 != 0) {
		    goto L220;
		}
		i__4 = do_lio(&c__3, &c__1, (char *)&stop, (ftnlen)sizeof(
			integer));
		if (i__4 != 0) {
		    goto L220;
		}
		i__4 = do_lio(&c__3, &c__1, (char *)&step, (ftnlen)sizeof(
			integer));
		if (i__4 != 0) {
		    goto L220;
		}
		i__4 = e_rsli();
		if (i__4 != 0) {
		    goto L220;
		}
L220:
		;
	    }
	    if (stop == 0) {
		stop = start;
	    }
	    if (step == 0) {
		step = 1;
	    }
	}
    }

/*     perform any final tasks before program exit */

    cl__1.cerr = 0;
    cl__1.cunit = iarc;
    cl__1.csta = 0;
    f_clos(&cl__1);
    final_();
    return 0;
} /* MAIN__ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine prtarc  --  output of a TINKER archive file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "prtarc" writes out a set of Cartesian coordinates for */
/*     all active atoms in the TINKER XYZ archive format */


/* Subroutine */ int prtarc_(integer *iarc)
{
    /* Format strings */
    static char fmt_10[] = "(i6)";
    static char fmt_20[] = "(i6,2x,a)";
    static char fmt_30[] = "(i6,2x,a3,3f12.6,9i6)";
    static char fmt_40[] = "(i6,2x,a3,3f14.8,9i6)";
    static char fmt_50[] = "(i6,2x,a3,3f16.10,9i6)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3;
    olist o__1;
    cllist cl__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void), f_clos(cllist *);

    /* Local variables */
    static integer i__, k;
    static logical opened;
    static char arcfile[120];
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___55 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___56 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___58 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___60 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_50, 0 };



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     open output unit if not already done */

    ioin__1.inerr = 0;
    ioin__1.inunit = *iarc;
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
	i__1[1] = 4, a__1[1] = ".arc";
	s_cat(arcfile, a__1, i__1, &c__2, (ftnlen)120);
	version_(arcfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = *iarc;
	o__1.ofnmlen = 120;
	o__1.ofnm = arcfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }

/*     write out the number of atoms and the title */

    if (titles_1.ltitle == 0) {
	io___55.ciunit = *iarc;
	s_wsfe(&io___55);
	do_fio(&c__1, (char *)&usage_1.nuse, (ftnlen)sizeof(integer));
	e_wsfe();
    } else {
	io___56.ciunit = *iarc;
	s_wsfe(&io___56);
	do_fio(&c__1, (char *)&usage_1.nuse, (ftnlen)sizeof(integer));
	do_fio(&c__1, titles_1.title, titles_1.ltitle);
	e_wsfe();
    }

/*     finally, write the coordinates for each atom */

    if (inform_1.digits <= 6) {
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (usage_1.use[i__ - 1]) {
		io___58.ciunit = *iarc;
		s_wsfe(&io___58);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
		do_fio(&c__1, (char *)&atoms_1.x[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&atoms_1.y[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&atoms_1.z__[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&atoms_1.type__[i__ - 1], (ftnlen)
			sizeof(integer));
		i__3 = couple_1.n12[i__ - 1];
		for (k = 1; k <= i__3; ++k) {
		    do_fio(&c__1, (char *)&i12_ref(k, i__), (ftnlen)sizeof(
			    integer));
		}
		e_wsfe();
	    }
	}
    } else if (inform_1.digits <= 8) {
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (usage_1.use[i__ - 1]) {
		io___60.ciunit = *iarc;
		s_wsfe(&io___60);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
		do_fio(&c__1, (char *)&atoms_1.x[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&atoms_1.y[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&atoms_1.z__[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&atoms_1.type__[i__ - 1], (ftnlen)
			sizeof(integer));
		i__3 = couple_1.n12[i__ - 1];
		for (k = 1; k <= i__3; ++k) {
		    do_fio(&c__1, (char *)&i12_ref(k, i__), (ftnlen)sizeof(
			    integer));
		}
		e_wsfe();
	    }
	}
    } else {
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (usage_1.use[i__ - 1]) {
		io___61.ciunit = *iarc;
		s_wsfe(&io___61);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
		do_fio(&c__1, (char *)&atoms_1.x[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&atoms_1.y[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&atoms_1.z__[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&atoms_1.type__[i__ - 1], (ftnlen)
			sizeof(integer));
		i__3 = couple_1.n12[i__ - 1];
		for (k = 1; k <= i__3; ++k) {
		    do_fio(&c__1, (char *)&i12_ref(k, i__), (ftnlen)sizeof(
			    integer));
		}
		e_wsfe();
	    }
	}
    }
    if (! opened) {
	cl__1.cerr = 0;
	cl__1.cunit = *iarc;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    return 0;
} /* prtarc_ */

#undef name___ref
#undef i12_ref


/* Main program alias */ int archive_ () { MAIN__ (); return 0; }
