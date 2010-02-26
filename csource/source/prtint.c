/* prtint.f -- translated by f2c (version 20050501).
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



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine prtint  --  output of internal coordinates  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "prtint" writes out a set of Z-matrix internal */
/*     coordinates to an external disk file */


/* Subroutine */ int prtint_(integer *izmt)
{
    /* Format strings */
    static char fmt_10[] = "(i6)";
    static char fmt_20[] = "(i6,2x,a)";
    static char fmt_30[] = "(i6,2x,a3,i6)";
    static char fmt_40[] = "(i6,2x,a3,2i6,f10.5)";
    static char fmt_50[] = "(i6,2x,a3,2i6,f12.7)";
    static char fmt_60[] = "(i6,2x,a3,2i6,f14.9)";
    static char fmt_70[] = "(i6,2x,a3,2i6,f10.5,i6,f10.4)";
    static char fmt_80[] = "(i6,2x,a3,2i6,f12.7,i6,f12.6)";
    static char fmt_90[] = "(i6,2x,a3,2i6,f14.9,i6,f14.8)";
    static char fmt_100[] = "(i6,2x,a3,2i6,f10.5,i6,f10.4,i6,f10.4,i6)";
    static char fmt_110[] = "(i6,2x,a3,2i6,f12.7,i6,f12.6,i6,f12.6,i6)";
    static char fmt_120[] = "(i6,2x,a3,2i6,f14.9,i6,f14.8,i6,f14.8,i6)";
    static char fmt_130[] = "()";
    static char fmt_140[] = "(2i6)";
    static char fmt_150[] = "()";
    static char fmt_160[] = "(2i6)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2;
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
    static char zmtfile[120];
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_160, 0 };



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




/*     open output unit if not already done */

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
	s_cat(zmtfile, a__1, i__1, &c__2, (ftnlen)120);
	version_(zmtfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = *izmt;
	o__1.ofnmlen = 120;
	o__1.ofnm = zmtfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }

/*     write out the number of atoms and the title */

    if (titles_1.ltitle == 0) {
	io___3.ciunit = *izmt;
	s_wsfe(&io___3);
	do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
	e_wsfe();
    } else {
	io___4.ciunit = *izmt;
	s_wsfe(&io___4);
	do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
	do_fio(&c__1, titles_1.title, titles_1.ltitle);
	e_wsfe();
    }

/*     output of first three atoms is handled separately */

    if (atoms_1.n >= 1) {
	io___5.ciunit = *izmt;
	s_wsfe(&io___5);
	do_fio(&c__1, (char *)&c__1, (ftnlen)sizeof(integer));
	do_fio(&c__1, name___ref(0, 1), (ftnlen)3);
	do_fio(&c__1, (char *)&atoms_1.type__[0], (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (atoms_1.n >= 2) {
	if (inform_1.digits <= 6) {
	    io___6.ciunit = *izmt;
	    s_wsfe(&io___6);
	    do_fio(&c__1, (char *)&c__2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, 2), (ftnlen)3);
	    do_fio(&c__1, (char *)&atoms_1.type__[1], (ftnlen)sizeof(integer))
		    ;
	    do_fio(&c__1, (char *)&iz_ref(1, 2), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zbond[1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	} else if (inform_1.digits <= 8) {
	    io___7.ciunit = *izmt;
	    s_wsfe(&io___7);
	    do_fio(&c__1, (char *)&c__2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, 2), (ftnlen)3);
	    do_fio(&c__1, (char *)&atoms_1.type__[1], (ftnlen)sizeof(integer))
		    ;
	    do_fio(&c__1, (char *)&iz_ref(1, 2), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zbond[1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	} else {
	    io___8.ciunit = *izmt;
	    s_wsfe(&io___8);
	    do_fio(&c__1, (char *)&c__2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, 2), (ftnlen)3);
	    do_fio(&c__1, (char *)&atoms_1.type__[1], (ftnlen)sizeof(integer))
		    ;
	    do_fio(&c__1, (char *)&iz_ref(1, 2), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zbond[1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }
    if (atoms_1.n >= 3) {
	if (inform_1.digits <= 6) {
	    io___9.ciunit = *izmt;
	    s_wsfe(&io___9);
	    do_fio(&c__1, (char *)&c__3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, 3), (ftnlen)3);
	    do_fio(&c__1, (char *)&atoms_1.type__[2], (ftnlen)sizeof(integer))
		    ;
	    do_fio(&c__1, (char *)&iz_ref(1, 3), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zbond[2], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iz_ref(2, 3), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zang[2], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	} else if (inform_1.digits <= 8) {
	    io___10.ciunit = *izmt;
	    s_wsfe(&io___10);
	    do_fio(&c__1, (char *)&c__3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, 3), (ftnlen)3);
	    do_fio(&c__1, (char *)&atoms_1.type__[2], (ftnlen)sizeof(integer))
		    ;
	    do_fio(&c__1, (char *)&iz_ref(1, 3), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zbond[2], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iz_ref(2, 3), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zang[2], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	} else {
	    io___11.ciunit = *izmt;
	    s_wsfe(&io___11);
	    do_fio(&c__1, (char *)&c__3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, 3), (ftnlen)3);
	    do_fio(&c__1, (char *)&atoms_1.type__[2], (ftnlen)sizeof(integer))
		    ;
	    do_fio(&c__1, (char *)&iz_ref(1, 3), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zbond[2], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iz_ref(2, 3), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zang[2], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }

/*     convert the torsional angles to lie in standard range */

    i__2 = atoms_1.n;
    for (i__ = 4; i__ <= i__2; ++i__) {
	if (iz_ref(4, i__) == 0) {
	    while(zcoord_1.ztors[i__ - 1] < -180.) {
		zcoord_1.ztors[i__ - 1] += 360.;
	    }
	    while(zcoord_1.ztors[i__ - 1] > 180.) {
		zcoord_1.ztors[i__ - 1] += -360.;
	    }
	}
    }

/*     now, output the fourth through final atoms */

    if (inform_1.digits <= 6) {
	i__2 = atoms_1.n;
	for (i__ = 4; i__ <= i__2; ++i__) {
	    io___13.ciunit = *izmt;
	    s_wsfe(&io___13);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
	    do_fio(&c__1, (char *)&atoms_1.type__[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&iz_ref(1, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zbond[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iz_ref(2, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zang[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iz_ref(3, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.ztors[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iz_ref(4, i__), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    } else if (inform_1.digits <= 8) {
	i__2 = atoms_1.n;
	for (i__ = 4; i__ <= i__2; ++i__) {
	    io___14.ciunit = *izmt;
	    s_wsfe(&io___14);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
	    do_fio(&c__1, (char *)&atoms_1.type__[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&iz_ref(1, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zbond[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iz_ref(2, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zang[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iz_ref(3, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.ztors[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iz_ref(4, i__), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    } else {
	i__2 = atoms_1.n;
	for (i__ = 4; i__ <= i__2; ++i__) {
	    io___15.ciunit = *izmt;
	    s_wsfe(&io___15);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
	    do_fio(&c__1, (char *)&atoms_1.type__[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&iz_ref(1, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zbond[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iz_ref(2, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zang[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iz_ref(3, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.ztors[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iz_ref(4, i__), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }

/*     finally, add or delete bonds as required */

    if (zclose_1.nadd != 0 || zclose_1.ndel != 0) {
	io___16.ciunit = *izmt;
	s_wsfe(&io___16);
	e_wsfe();
    }
    i__2 = zclose_1.nadd;
    for (i__ = 1; i__ <= i__2; ++i__) {
	io___17.ciunit = *izmt;
	s_wsfe(&io___17);
	for (k = 1; k <= 2; ++k) {
	    do_fio(&c__1, (char *)&iadd_ref(k, i__), (ftnlen)sizeof(integer));
	}
	e_wsfe();
    }
    if (zclose_1.ndel != 0) {
	io___19.ciunit = *izmt;
	s_wsfe(&io___19);
	e_wsfe();
    }
    i__2 = zclose_1.ndel;
    for (i__ = 1; i__ <= i__2; ++i__) {
	io___20.ciunit = *izmt;
	s_wsfe(&io___20);
	for (k = 1; k <= 2; ++k) {
	    do_fio(&c__1, (char *)&idel_ref(k, i__), (ftnlen)sizeof(integer));
	}
	e_wsfe();
    }
    if (! opened) {
	cl__1.cerr = 0;
	cl__1.cunit = *izmt;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    return 0;
} /* prtint_ */

#undef name___ref
#undef idel_ref
#undef iadd_ref
#undef iz_ref


