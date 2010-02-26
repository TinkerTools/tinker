/* katom.f -- translated by f2c (version 20050501).
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
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal weight[5000];
    integer atmcls[5000], atmnum[5000], ligand[5000];
    char symbol[15000], describe[120000];
} katoms_;

#define katoms_1 katoms_

struct {
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine katom  --  atom type parameter assignment  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "katom" assigns an atom type definitions to each atom in */
/*     the structure and processes any new or changed values */


/* Subroutine */ int katom_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Additional Atom Type Parameters :\002,//"
	    ",5x,\002Type  Class  Symbol  Description\002,15x,\002Atomic\002,"
	    "4x,\002Mass\002,3x,\002Valence\002,/)";
    static char fmt_20[] = "(2x,i6,1x,i6,5x,a3,3x,a24,i6,f11.3,i6)";
    static char fmt_30[] = "(/,\002 KATOM   --  Too many Atom Types;\002,"
	    "\002 Increase MAXTYP\002)";
    static char fmt_50[] = "(/,\002 Additional Atom Types for\002,\002 Speci"
	    "fic Atoms :\002,//,5x,\002Atom  Class  Symbol  Description\002,1"
	    "5x,\002Atomic\002,4x,\002Mass\002,3x,\002Valence\002,/)";
    static char fmt_60[] = "(2x,i6,1x,i6,5x,a3,3x,a24,i6,f11.3,i6)";
    static char fmt_80[] = "(/,\002 Undefined Atom Types or Classes :\002,"
	    "//,\002 Type\002,10x,\002Atom Number\002,5x,\002Atom Type\002,5x,"
	    "\002Atom Class\002,/)";
    static char fmt_90[] = "(\002 Atom\002,12x,i5,10x,i5,10x,i5)";
    static char fmt_100[] = "(/,\002 Atoms with an Unusual Number of Attac"
	    "hed\002,\002 Atoms :\002,//,\002 Type\002,11x,\002Atom Name\002,"
	    "6x,\002Atom Type\002,7x,\002Expected\002,4x,\002Found\002,/)";
    static char fmt_110[] = "(\002 Valence\002,7x,i5,\002-\002,a3,8x,i5,10x,"
	    "i5,5x,i5)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, k;
    extern /* Subroutine */ int getstring_(char *, char *, integer *, ftnlen, 
	    ftnlen);
    static integer lig, cls, atn;
    static doublereal wght;
    static char symb[3];
    static integer next;
    static logical header;
    static char record[120], notice[24], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), getnumb_(char *, 
	    integer *, integer *, ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___14 = { 1, string, 1, 0, 120, 1 };
    static cilist io___15 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_30, 0 };
    static icilist io___18 = { 1, string, 1, 0, 120, 1 };
    static cilist io___19 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_110, 0 };



#define describe_ref(a_0,a_1) &katoms_1.describe[(a_1)*24 + a_0 - 24]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define story_ref(a_0,a_1) &atmtyp_1.story[(a_1)*24 + a_0 - 24]
#define symbol_ref(a_0,a_1) &katoms_1.symbol[(a_1)*3 + a_0 - 3]
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
/*     ##  katoms.i  --  forcefield parameters for the atom types  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     weight     average atomic mass of each atom type */
/*     atmcls     atom class number for each of the atom types */
/*     atmnum     atomic number for each of the atom types */
/*     ligand     number of atoms to be attached to each atom type */
/*     symbol     modified atomic symbol for each atom type */
/*     describe   string identifying each of the atom types */




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




/*     process keywords containing atom type parameters */

    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "ATOM ", (ftnlen)5, (ftnlen)5) == 0) {
	    k = 0;
	    cls = 0;
	    s_copy(symb, " ", (ftnlen)3, (ftnlen)1);
	    s_copy(notice, " ", (ftnlen)24, (ftnlen)1);
	    atn = 0;
	    wght = 0.;
	    lig = 0;
	    getnumb_(record, &k, &next, (ftnlen)120);
	    getnumb_(record, &cls, &next, (ftnlen)120);
	    if (cls == 0) {
		cls = k;
	    }
	    katoms_1.atmcls[k - 1] = cls;
	    gettext_(record, symb, &next, (ftnlen)120, (ftnlen)3);
	    getstring_(record, notice, &next, (ftnlen)120, (ftnlen)24);
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___14);
	    if (i__2 != 0) {
		goto L40;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&atn, (ftnlen)sizeof(integer))
		    ;
	    if (i__2 != 0) {
		goto L40;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&wght, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L40;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&lig, (ftnlen)sizeof(integer))
		    ;
	    if (i__2 != 0) {
		goto L40;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L40;
	    }
	    if (k >= 1 && k <= 5000) {
		if (header) {
		    header = FALSE_;
		    io___15.ciunit = iounit_1.iout;
		    s_wsfe(&io___15);
		    e_wsfe();
		}
		s_copy(symbol_ref(0, k), symb, (ftnlen)3, (ftnlen)3);
		s_copy(describe_ref(0, k), notice, (ftnlen)24, (ftnlen)24);
		katoms_1.atmnum[k - 1] = atn;
		katoms_1.weight[k - 1] = wght;
		katoms_1.ligand[k - 1] = lig;
		io___16.ciunit = iounit_1.iout;
		s_wsfe(&io___16);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&cls, (ftnlen)sizeof(integer));
		do_fio(&c__1, symb, (ftnlen)3);
		do_fio(&c__1, notice, (ftnlen)24);
		do_fio(&c__1, (char *)&atn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&wght, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&lig, (ftnlen)sizeof(integer));
		e_wsfe();
	    } else if (k >= 5000) {
		io___17.ciunit = iounit_1.iout;
		s_wsfe(&io___17);
		e_wsfe();
		inform_1.abort = TRUE_;
	    }
L40:
	    ;
	}
    }

/*     transfer atom type values to individual atoms */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = atoms_1.type__[i__ - 1];
	if (k == 0) {
	    atmtyp_1.class__[i__ - 1] = 0;
	    atmtyp_1.atomic[i__ - 1] = 0;
	    atmtyp_1.mass[i__ - 1] = 0.;
	    atmtyp_1.valence[i__ - 1] = 0;
	    s_copy(story_ref(0, i__), "Undefined Atom Type     ", (ftnlen)24, 
		    (ftnlen)24);
	} else {
	    if (s_cmp(symbol_ref(0, k), "   ", (ftnlen)3, (ftnlen)3) != 0) {
		s_copy(name___ref(0, i__), symbol_ref(0, k), (ftnlen)3, (
			ftnlen)3);
	    }
	    atmtyp_1.class__[i__ - 1] = katoms_1.atmcls[k - 1];
	    atmtyp_1.atomic[i__ - 1] = katoms_1.atmnum[k - 1];
	    atmtyp_1.mass[i__ - 1] = katoms_1.weight[k - 1];
	    atmtyp_1.valence[i__ - 1] = katoms_1.ligand[k - 1];
	    s_copy(story_ref(0, i__), describe_ref(0, k), (ftnlen)24, (ftnlen)
		    24);
	}
    }

/*     process keywords containing atom types for specific atoms */

    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "ATOM ", (ftnlen)5, (ftnlen)5) == 0) {
	    k = 0;
	    s_copy(symb, " ", (ftnlen)3, (ftnlen)1);
	    s_copy(notice, " ", (ftnlen)24, (ftnlen)1);
	    atn = 0;
	    wght = 0.;
	    lig = 0;
	    getnumb_(record, &k, &next, (ftnlen)120);
	    getnumb_(record, &cls, &next, (ftnlen)120);
	    gettext_(record, symb, &next, (ftnlen)120, (ftnlen)3);
	    getstring_(record, notice, &next, (ftnlen)120, (ftnlen)24);
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___18);
	    if (i__2 != 0) {
		goto L70;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&atn, (ftnlen)sizeof(integer))
		    ;
	    if (i__2 != 0) {
		goto L70;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&wght, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L70;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&lig, (ftnlen)sizeof(integer))
		    ;
	    if (i__2 != 0) {
		goto L70;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L70;
	    }
	    if (k < 0 && k >= -atoms_1.n) {
		if (header) {
		    header = FALSE_;
		    io___19.ciunit = iounit_1.iout;
		    s_wsfe(&io___19);
		    e_wsfe();
		}
		k = -k;
		if (cls == 0) {
		    cls = k;
		}
		atmtyp_1.class__[k - 1] = cls;
		s_copy(name___ref(0, k), symb, (ftnlen)3, (ftnlen)3);
		s_copy(story_ref(0, k), notice, (ftnlen)24, (ftnlen)24);
		atmtyp_1.atomic[k - 1] = atn;
		atmtyp_1.mass[k - 1] = wght;
		atmtyp_1.valence[k - 1] = lig;
		io___20.ciunit = iounit_1.iout;
		s_wsfe(&io___20);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&cls, (ftnlen)sizeof(integer));
		do_fio(&c__1, symb, (ftnlen)3);
		do_fio(&c__1, notice, (ftnlen)24);
		do_fio(&c__1, (char *)&atn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&wght, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&lig, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
L70:
	    ;
	}
    }

/*     check for presence of undefined atom types or classes */

    header = TRUE_;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = atoms_1.type__[i__ - 1];
	cls = atmtyp_1.class__[i__ - 1];
	if (k < 1 || k > 5000 || cls < 1 || cls > 1000) {
	    inform_1.abort = TRUE_;
	    if (header) {
		header = FALSE_;
		io___21.ciunit = iounit_1.iout;
		s_wsfe(&io___21);
		e_wsfe();
	    }
	    io___22.ciunit = iounit_1.iout;
	    s_wsfe(&io___22);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&cls, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }

/*     check the number of atoms attached to each atom */

    header = TRUE_;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (couple_1.n12[i__ - 1] != atmtyp_1.valence[i__ - 1]) {
	    if (header) {
		header = FALSE_;
		io___23.ciunit = iounit_1.iout;
		s_wsfe(&io___23);
		e_wsfe();
	    }
	    io___24.ciunit = iounit_1.iout;
	    s_wsfe(&io___24);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
	    do_fio(&c__1, (char *)&atoms_1.type__[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&atmtyp_1.valence[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&couple_1.n12[i__ - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	}
    }
    return 0;
} /* katom_ */

#undef keyline_ref
#undef symbol_ref
#undef story_ref
#undef name___ref
#undef describe_ref


