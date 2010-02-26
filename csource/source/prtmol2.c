/* prtmol2.f -- translated by f2c (version 20050501).
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
    doublereal bk[50000], bl[50000];
    integer nbond, ibnd[100000]	/* was [2][50000] */;
} bond_;

#define bond_1 bond_

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



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  program prtmol2  --  output of Sybyl MOL2 coordinate file  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "prtmol2" writes out a set of coordinates in Sybyl MOL2 */
/*     format to an external disk file */


/* Subroutine */ int prtmol2_(integer *isyb)
{
    /* Initialized data */

    static char digit[1*10] = "0" "1" "2" "3" "4" "5" "6" "7" "8" "9";

    /* Format strings */
    static char fmt_10[] = "(\002@<TRIPOS>MOLECULE\002)";
    static char fmt_20[] = "(\002****\002)";
    static char fmt_30[] = "(a)";
    static char fmt_40[] = "(3i8)";
    static char fmt_50[] = "(\002SMALL\002)";
    static char fmt_60[] = "(\002NO_CHARGES\002)";
    static char fmt_70[] = "(/,\002@<TRIPOS>ATOM\002)";
    static char fmt_80[] = "(i8,3x,a7,2x,3f12.6,3x,a4)";
    static char fmt_90[] = "(/,\002@<TRIPOS>BOND\002)";
    static char fmt_100[] = "(4i8)";
    static char fmt_110[] = "(/,\002@<TRIPOS>SUBSTRUCTURE\002)";
    static char fmt_120[] = "(i8,12x,a4,i8)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3;
    olist o__1;
    cllist cl__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), s_wsfe(cilist *), e_wsfe(void), do_fio(integer *,
	     char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), f_clos(cllist *);

    /* Local variables */
    static integer thousand, i__, j, k;
    static char sybylfile[120];
    static integer ones, tens;
    static logical opened;
    static char atmnam[7], number[4], atmtyp[4];
    static integer hundred;
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_120, 0 };



#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
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




/*     open output unit if not already done */

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
	version_(sybylfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = *isyb;
	o__1.ofnmlen = 120;
	o__1.ofnm = sybylfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }

/*     write the molecule record type indicator */

    io___4.ciunit = *isyb;
    s_wsfe(&io___4);
    e_wsfe();
    if (titles_1.ltitle == 0) {
	io___5.ciunit = *isyb;
	s_wsfe(&io___5);
	e_wsfe();
    } else {
	io___6.ciunit = *isyb;
	s_wsfe(&io___6);
	do_fio(&c__1, titles_1.title, titles_1.ltitle);
	e_wsfe();
    }
    io___7.ciunit = *isyb;
    s_wsfe(&io___7);
    do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&bond_1.nbond, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&c__1, (ftnlen)sizeof(integer));
    e_wsfe();
    io___8.ciunit = *isyb;
    s_wsfe(&io___8);
    e_wsfe();
    io___9.ciunit = *isyb;
    s_wsfe(&io___9);
    e_wsfe();

/*     write the atom record type indicator */

    io___10.ciunit = *isyb;
    s_wsfe(&io___10);
    e_wsfe();
    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {

/*     set the Sybyl atom_name for the atom */

	thousand = i__ / 1000;
	hundred = (i__ - thousand * 1000) / 100;
	tens = (i__ - thousand * 1000 - hundred * 100) / 10;
	ones = i__ - thousand * 1000 - hundred * 100 - tens * 10;
	*(unsigned char *)number = *(unsigned char *)&digit[thousand];
	*(unsigned char *)&number[1] = *(unsigned char *)&digit[hundred];
	*(unsigned char *)&number[2] = *(unsigned char *)&digit[tens];
	*(unsigned char *)&number[3] = *(unsigned char *)&digit[ones];
	if (*(unsigned char *)number == '0') {
	    *(unsigned char *)number = ' ';
	}
	if (*(unsigned char *)&number[1] == '0' && *(unsigned char *)number ==
		 ' ') {
	    *(unsigned char *)&number[1] = ' ';
	}
	if (*(unsigned char *)&number[2] == '0' && *(unsigned char *)&number[
		1] == ' ') {
	    *(unsigned char *)&number[2] = ' ';
	}
/* Writing concatenation */
	i__1[0] = 3, a__1[0] = name___ref(0, i__);
	i__1[1] = 4, a__1[1] = number;
	s_cat(atmnam, a__1, i__1, &c__2, (ftnlen)7);
	for (j = 1; j <= 6; ++j) {
	    while(*(unsigned char *)&atmnam[j - 1] == ' ') {
		for (k = j; k <= 6; ++k) {
		    i__3 = k;
		    s_copy(atmnam + (k - 1), atmnam + i__3, (ftnlen)1, k + 1 
			    - i__3);
		}
		*(unsigned char *)&atmnam[6] = '*';
	    }
	}
	for (j = 1; j <= 7; ++j) {
	    if (*(unsigned char *)&atmnam[j - 1] == '*') {
		*(unsigned char *)&atmnam[j - 1] = ' ';
	    }
	}

/*     set the Sybyl atom_type for the atom */

/* Writing concatenation */
	i__1[0] = 3, a__1[0] = name___ref(0, i__);
	i__1[1] = 1, a__1[1] = " ";
	s_cat(atmtyp, a__1, i__1, &c__2, (ftnlen)4);
	if (s_cmp(atmtyp, "C  ", (ftnlen)4, (ftnlen)3) == 0) {
	    if (couple_1.n12[i__ - 1] == 4) {
		s_copy(atmtyp, "C.3 ", (ftnlen)4, (ftnlen)4);
	    }
	    if (couple_1.n12[i__ - 1] == 3) {
		s_copy(atmtyp, "C.2 ", (ftnlen)4, (ftnlen)4);
	    }
	    if (couple_1.n12[i__ - 1] == 2) {
		s_copy(atmtyp, "C.1 ", (ftnlen)4, (ftnlen)4);
	    }
	} else if (s_cmp(atmtyp, "N  ", (ftnlen)4, (ftnlen)3) == 0) {
	    if (couple_1.n12[i__ - 1] >= 3) {
		s_copy(atmtyp, "N.3 ", (ftnlen)4, (ftnlen)4);
	    }
	    if (couple_1.n12[i__ - 1] == 2) {
		s_copy(atmtyp, "N.2 ", (ftnlen)4, (ftnlen)4);
	    }
	    if (couple_1.n12[i__ - 1] == 1) {
		s_copy(atmtyp, "N.1 ", (ftnlen)4, (ftnlen)4);
	    }
	} else if (s_cmp(atmtyp, "N+ ", (ftnlen)4, (ftnlen)3) == 0) {
	    s_copy(atmtyp, "N.4 ", (ftnlen)4, (ftnlen)4);
	} else if (s_cmp(atmtyp, "O  ", (ftnlen)4, (ftnlen)3) == 0) {
	    if (couple_1.n12[i__ - 1] >= 2) {
		s_copy(atmtyp, "O.3 ", (ftnlen)4, (ftnlen)4);
	    }
	    if (couple_1.n12[i__ - 1] <= 1) {
		s_copy(atmtyp, "O.2 ", (ftnlen)4, (ftnlen)4);
	    }
	} else if (s_cmp(atmtyp, "O- ", (ftnlen)4, (ftnlen)3) == 0) {
	    s_copy(atmtyp, "O.2 ", (ftnlen)4, (ftnlen)4);
	} else if (s_cmp(atmtyp, "S  ", (ftnlen)4, (ftnlen)3) == 0) {
	    if (couple_1.n12[i__ - 1] >= 2) {
		s_copy(atmtyp, "S.3 ", (ftnlen)4, (ftnlen)4);
	    }
	    if (couple_1.n12[i__ - 1] <= 1) {
		s_copy(atmtyp, "S.2 ", (ftnlen)4, (ftnlen)4);
	    }
	} else if (s_cmp(atmtyp, "P  ", (ftnlen)4, (ftnlen)3) == 0) {
	    s_copy(atmtyp, "P.3 ", (ftnlen)4, (ftnlen)4);
	} else if (s_cmp(atmtyp, "Lp ", (ftnlen)4, (ftnlen)3) == 0) {
	    s_copy(atmtyp, "LP  ", (ftnlen)4, (ftnlen)4);
	}
	io___21.ciunit = *isyb;
	s_wsfe(&io___21);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, atmnam, (ftnlen)7);
	do_fio(&c__1, (char *)&atoms_1.x[i__ - 1], (ftnlen)sizeof(doublereal))
		;
	do_fio(&c__1, (char *)&atoms_1.y[i__ - 1], (ftnlen)sizeof(doublereal))
		;
	do_fio(&c__1, (char *)&atoms_1.z__[i__ - 1], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, atmtyp, (ftnlen)4);
	e_wsfe();
    }

/*     write the bond record type indicator */

    io___22.ciunit = *isyb;
    s_wsfe(&io___22);
    e_wsfe();
    i__2 = bond_1.nbond;
    for (i__ = 1; i__ <= i__2; ++i__) {
	io___23.ciunit = *isyb;
	s_wsfe(&io___23);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	for (j = 1; j <= 2; ++j) {
	    do_fio(&c__1, (char *)&ibnd_ref(j, i__), (ftnlen)sizeof(integer));
	}
	do_fio(&c__1, (char *)&c__1, (ftnlen)sizeof(integer));
	e_wsfe();
    }

/*     write the substructure record type indicator */

    io___24.ciunit = *isyb;
    s_wsfe(&io___24);
    e_wsfe();
    io___25.ciunit = *isyb;
    s_wsfe(&io___25);
    do_fio(&c__1, (char *)&c__1, (ftnlen)sizeof(integer));
    do_fio(&c__1, "****", (ftnlen)4);
    do_fio(&c__1, (char *)&c__1, (ftnlen)sizeof(integer));
    e_wsfe();
    if (! opened) {
	cl__1.cerr = 0;
	cl__1.cunit = *isyb;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    return 0;
} /* prtmol2_ */

#undef name___ref
#undef ibnd_ref


