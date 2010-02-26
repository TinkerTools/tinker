/* readmol2.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine readmol2  --  input of a Sybyl MOL2 file  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "readmol2" gets a set of Sybyl MOL2 coordinates from */
/*     an external disk file */


/* Subroutine */ int readmol2_(integer *isyb)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 READMOL2  --  Unable to Find the TRIPO"
	    "S\002,\002 Sybyl MOL2 File\002)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(a120)";
    static char fmt_40[] = "(a120)";
    static char fmt_60[] = "(/,\002 READMOL2  --  The Maximum of\002,i8,\002"
	    " Atoms\002,\002 has been Exceeded\002)";
    static char fmt_70[] = "(a120)";
    static char fmt_80[] = "(a120)";
    static char fmt_100[] = "(a120)";
    static char fmt_110[] = "(a120)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_rew(alist *), s_wsfe(cilist *), e_wsfe(void), 
	    s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void),
	     s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_clos(cllist *);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j, k, m;
    static char sybylfile[120];
    static integer ia, ib, next;
    extern /* Subroutine */ int sort_(integer *, integer *), fatal_(void);
    static integer nbond;
    static logical exist, opened;
    static char atmnam[10], record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static integer number;
    static char string[120];
    extern /* Subroutine */ int getword_(char *, char *, integer *, ftnlen, 
	    ftnlen), gettext_(char *, char *, integer *, ftnlen, ftnlen), 
	    version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_40, 0 };
    static icilist io___12 = { 0, record, 0, 0, 120, 1 };
    static cilist io___14 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_80, 0 };
    static icilist io___18 = { 0, record, 0, 0, 120, 1 };
    static icilist io___21 = { 0, string, 0, 0, 120, 1 };
    static cilist io___24 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_110, 0 };
    static icilist io___26 = { 0, record, 0, 0, 120, 1 };



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




/*     open the input file if it has not already been done */

    ioin__1.inerr = 0;
    ioin__1.inunit = *isyb;
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
	i__1[1] = 5, a__1[1] = ".mol2";
	s_cat(sybylfile, a__1, i__1, &c__2, (ftnlen)120);
	version_(sybylfile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = sybylfile;
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
	    o__1.ounit = *isyb;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = sybylfile;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	    al__1.aerr = 0;
	    al__1.aunit = *isyb;
	    f_rew(&al__1);
	} else {
	    io___4.ciunit = iounit_1.iout;
	    s_wsfe(&io___4);
	    e_wsfe();
	    fatal_();
	}
    }

/*     get title line and get the number of atoms and bonds */

    for (i__ = 1; i__ <= 1000000; ++i__) {
	io___6.ciunit = *isyb;
	s_rsfe(&io___6);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, string, &next, (ftnlen)120, (ftnlen)120);
	upcase_(string, (ftnlen)120);
	if (s_cmp(string, "@<TRIPOS>MOLECULE", (ftnlen)120, (ftnlen)17) == 0) 
		{
	    io___10.ciunit = *isyb;
	    s_rsfe(&io___10);
	    do_fio(&c__1, titles_1.title, (ftnlen)120);
	    e_rsfe();
	    titles_1.ltitle = trimtext_(titles_1.title, (ftnlen)120);
	    io___11.ciunit = *isyb;
	    s_rsfe(&io___11);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    s_rsli(&io___12);
	    do_lio(&c__3, &c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&nbond, (ftnlen)sizeof(integer));
	    e_rsli();
	    goto L50;
	}
    }
L50:

/*     check for too many total atoms in the file */

    if (atoms_1.n > 25000) {
	io___14.ciunit = iounit_1.iout;
	s_wsfe(&io___14);
	do_fio(&c__1, (char *)&c__25000, (ftnlen)sizeof(integer));
	e_wsfe();
	fatal_();
    }

/*     read the atom names and coordinates */

    for (i__ = 1; i__ <= 1000000; ++i__) {
	io___15.ciunit = *isyb;
	s_rsfe(&io___15);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, string, &next, (ftnlen)120, (ftnlen)120);
	upcase_(string, (ftnlen)120);
	if (s_cmp(string, "@<TRIPOS>ATOM", (ftnlen)120, (ftnlen)13) == 0) {
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		io___17.ciunit = *isyb;
		s_rsfe(&io___17);
		do_fio(&c__1, record, (ftnlen)120);
		e_rsfe();
		s_rsli(&io___18);
		do_lio(&c__3, &c__1, (char *)&number, (ftnlen)sizeof(integer))
			;
		e_rsli();
		next = 1;
		getword_(record, atmnam, &next, (ftnlen)120, (ftnlen)10);
		s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next 
			- 1));
		s_rsli(&io___21);
		do_lio(&c__5, &c__1, (char *)&atoms_1.x[j - 1], (ftnlen)
			sizeof(doublereal));
		do_lio(&c__5, &c__1, (char *)&atoms_1.y[j - 1], (ftnlen)
			sizeof(doublereal));
		do_lio(&c__5, &c__1, (char *)&atoms_1.z__[j - 1], (ftnlen)
			sizeof(doublereal));
		e_rsli();
		getword_(record, atmnam, &next, (ftnlen)120, (ftnlen)10);
		s_copy(name___ref(0, j), atmnam, (ftnlen)3, (ftnlen)3);
		for (k = 1; k <= 3; ++k) {
		    if (*(unsigned char *)&atmnam[k - 1] == '.') {
			for (m = k; m <= 3; ++m) {
			    *(unsigned char *)name___ref(m - 1, j) = ' ';
			}
		    }
		}
		atoms_1.type__[j - 1] = 0;
	    }
	    goto L90;
	}
    }
L90:

/*     read the bond list to get attached atom lists */

    for (i__ = 1; i__ <= 1000000; ++i__) {
	io___24.ciunit = *isyb;
	s_rsfe(&io___24);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, string, &next, (ftnlen)120, (ftnlen)120);
	upcase_(string, (ftnlen)120);
	if (s_cmp(string, "@<TRIPOS>BOND", (ftnlen)120, (ftnlen)13) == 0) {
	    i__2 = nbond;
	    for (j = 1; j <= i__2; ++j) {
		io___25.ciunit = *isyb;
		s_rsfe(&io___25);
		do_fio(&c__1, record, (ftnlen)120);
		e_rsfe();
		s_rsli(&io___26);
		do_lio(&c__3, &c__1, (char *)&number, (ftnlen)sizeof(integer))
			;
		do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
		e_rsli();
		++couple_1.n12[ia - 1];
		i12_ref(couple_1.n12[ia - 1], ia) = ib;
		++couple_1.n12[ib - 1];
		i12_ref(couple_1.n12[ib - 1], ib) = ia;
	    }
	    goto L120;
	}
    }
L120:

/*     for each atom, sort its list of attached atoms */

    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	sort_(&couple_1.n12[i__ - 1], &i12_ref(1, i__));
    }
    if (! opened) {
	cl__1.cerr = 0;
	cl__1.cunit = *isyb;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    return 0;
} /* readmol2_ */

#undef name___ref
#undef i12_ref


