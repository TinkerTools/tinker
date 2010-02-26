/* intedit.f -- translated by f2c (version 20050501).
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
    doublereal zbond[25000], zang[25000], ztors[25000];
    integer iz[100000]	/* was [4][25000] */;
} zcoord_;

#define zcoord_1 zcoord_

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__0 = 0;
static integer c__5 = 5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  program intedit  --  edit and display Z-matrix file  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "intedit" allows the user to extract information from */
/*     or alter the values within an internal coordinates file */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 INTEDIT>  \002,$)";
    static char fmt_40[] = "(a120)";
    static char fmt_50[] = "(/,\002 Z-Matrix Internal Coordinates written to"
	    " :  \002,a)";
    static char fmt_60[] = "(/,\002 The Z-Matrix was not Changed;\002,\002 N"
	    "o File was Written\002)";
    static char fmt_70[] = "()";
    static char fmt_90[] = "(/,\002 Warning; Only\002,i6,\002 Atoms are Pres"
	    "ent\002,\002 in the Z-matrix\002)";
    static char fmt_100[] = "(/,\002 Atom Number :\002,i8)";
    static char fmt_110[] = "(\002 Atom Name :\002,6x,a4)";
    static char fmt_120[] = "(\002 Atom Type :\002,5x,a20)";
    static char fmt_130[] = "(\002 Type Number :\002,i8)";
    static char fmt_140[] = "(/,\002 Atom 1 is at the Coordinate System Orig"
	    "in\002)";
    static char fmt_150[] = "(/,\002 Internal Coordinate Structural Definiti"
	    "on :\002,/)";
    static char fmt_160[] = "(1x,2i6,17x,\002Distance Value :\002,f14.4)";
    static char fmt_170[] = "(1x,3i6,11x,\002Bond Angle Value :\002,f12.4)";
    static char fmt_180[] = "(1x,4i6,5x,\002Dihedral Angle :\002,f14.4)";
    static char fmt_190[] = "(1x,3i6,11x,\002Bond Angle Value :\002,f12.4)";
    static char fmt_200[] = "(30x,\002Chirality Flag :\002,6x,i6)";
    static char fmt_210[] = "(/,\002 Inverting Chirality of Atom : \002,i6)";
    static char fmt_220[] = "(/,\002 Invalid Atom Number\002)";
    static char fmt_230[] = "(/,\002 The Current Distance is : \002,f9.4)";
    static char fmt_240[] = "(\002 That Bond is not in the Z-matrix\002)";
    static char fmt_250[] = "(/,\002 The Current Distance is : \002,f9.4)";
    static char fmt_260[] = "(/,\002 Old Atom Type is :  \002,a20)";
    static char fmt_270[] = "(\002 New Atom Type is :  \002,a20)";
    static char fmt_280[] = "(/,\002 Invalid Atom Type; Valid Types are :"
	    "\002,/)";
    static char fmt_290[] = "(1x,3(i3,1x,a20,2x))";
    static char fmt_300[] = "(/,\002 Invalid Atom Number\002)";
    static char fmt_310[] = "(/,\002 The Bond Angle Value is :  \002,f9.4)";
    static char fmt_320[] = "(\002 That Bond Angle is not in the Z-matrix"
	    "\002)";
    static char fmt_330[] = "(/,\002 The Bond Angle Value is :  \002,f9.4)";
    static char fmt_340[] = "(/,\002 The Bond Angle Value is :  \002,f9.4)";
    static char fmt_350[] = "(/,\002 The Bond Angle Value is :  \002,f9.4)";
    static char fmt_360[] = "(\002 That Bond Angle is not in the Z-matrix"
	    "\002)";
    static char fmt_370[] = "(/,\002 Invalid Atom Number\002)";
    static char fmt_380[] = "(/,\002 The Dihedral Angle Value is :  \002,f9."
	    "4)";
    static char fmt_390[] = "(\002 That Dihedral Angle is not in the Z-mat"
	    "rix\002)";
    static char fmt_400[] = "(/,\002 The Dihedral Angle Value is :  \002,f9."
	    "4)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3, i__4;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void), s_cmp(char *, char *, ftnlen, 
	    ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);

    /* Local variables */
    extern integer freeunit_(void);
    extern doublereal geometry_(integer *, integer *, integer *, integer *);
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j, k, l, m;
    static char word[4];
    static integer next, izmt;
    extern /* Subroutine */ int field_(void), final_(void);
    static integer space;
    static doublereal value;
    extern /* Subroutine */ int zhelp_(void);
    static logical error;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static integer number[4], numcol;
    extern /* Subroutine */ int getint_(void), zvalue_(char *, doublereal *, 
	    logical *, ftnlen), prtint_(integer *);
    static integer numrow;
    static logical changed;
    extern /* Subroutine */ int initial_(void);
    static char zmtfile[120];
    extern /* Subroutine */ int getword_(char *, char *, integer *, ftnlen, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen), makexyz_(void);

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_70, 0 };
    static icilist io___17 = { 1, record, 1, 0, 120, 1 };
    static cilist io___18 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___46 = { 0, 0, 0, fmt_320, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_330, 0 };
    static cilist io___48 = { 0, 0, 0, fmt_340, 0 };
    static cilist io___49 = { 0, 0, 0, fmt_350, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_360, 0 };
    static cilist io___52 = { 0, 0, 0, fmt_370, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_380, 0 };
    static cilist io___54 = { 0, 0, 0, fmt_390, 0 };
    static cilist io___55 = { 0, 0, 0, fmt_400, 0 };



#define describe_ref(a_0,a_1) &katoms_1.describe[(a_1)*24 + a_0 - 24]
#define iz_ref(a_1,a_2) zcoord_1.iz[(a_2)*4 + a_1 - 5]
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
/*     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     zbond   bond length used to define each Z-matrix atom */
/*     zang    bond angle used to define each Z-matrix atom */
/*     ztors   angle or torsion used to define Z-matrix atom */
/*     iz      defining atom numbers for each Z-matrix atom */




/*     read coordinate file and force field definition */

    initial_();
    getint_();
    field_();

/*     print out the instructions for the program */

    next = 1;
    changed = FALSE_;
    error = FALSE_;
L10:
    zhelp_();

/*     start of main loop, examine or change Z-matrix elements */

L20:
    m = 0;
    io___5.ciunit = iounit_1.iout;
    s_wsfe(&io___5);
    e_wsfe();
    io___6.ciunit = iounit_1.input;
    s_rsfe(&io___6);
    do_fio(&c__1, record, (ftnlen)120);
    e_rsfe();

/*     interpret any user entered text command */

    space = 1;
    getword_(record, word, &space, (ftnlen)120, (ftnlen)4);
    upcase_(word, (ftnlen)4);
    if (s_cmp(word, "EXIT", (ftnlen)4, (ftnlen)4) == 0) {
	if (changed) {
	    izmt = freeunit_();
/* Writing concatenation */
	    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
	    i__1[1] = 4, a__1[1] = ".int";
	    s_cat(zmtfile, a__1, i__1, &c__2, (ftnlen)120);
	    version_(zmtfile, "new", (ftnlen)120, (ftnlen)3);
	    o__1.oerr = 0;
	    o__1.ounit = izmt;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = zmtfile;
	    o__1.orl = 0;
	    o__1.osta = "new";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	    prtint_(&izmt);
	    cl__1.cerr = 0;
	    cl__1.cunit = izmt;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	    io___12.ciunit = iounit_1.iout;
	    s_wsfe(&io___12);
	    do_fio(&c__1, zmtfile, trimtext_(zmtfile, (ftnlen)120));
	    e_wsfe();
	} else {
	    io___13.ciunit = iounit_1.iout;
	    s_wsfe(&io___13);
	    e_wsfe();
	}
	goto L410;
    } else if (s_cmp(word, "QUIT", (ftnlen)4, (ftnlen)4) == 0) {
	goto L410;
    } else if (s_cmp(word, "SHOW", (ftnlen)4, (ftnlen)4) == 0) {
	io___14.ciunit = iounit_1.iout;
	s_wsfe(&io___14);
	e_wsfe();
	prtint_(&iounit_1.iout);

/*     get the number of atoms entered by the user */

    } else {
	for (i__ = 1; i__ <= 4; ++i__) {
	    number[i__ - 1] = 0;
	}
	i__2 = s_rsli(&io___17);
	if (i__2 != 0) {
	    goto L100001;
	}
	for (i__ = 1; i__ <= 4; ++i__) {
	    i__2 = do_lio(&c__3, &c__1, (char *)&number[i__ - 1], (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L100001;
	    }
	}
	i__2 = e_rsli();
L100001:
	if (i__2 < 0) {
	    goto L80;
	}
	if (i__2 > 0) {
	    goto L10;
	}
L80:
	for (i__ = 1; i__ <= 4; ++i__) {
	    if (number[i__ - 1] != 0) {
		m = i__;
	    }
	    if (number[i__ - 1] > atoms_1.n) {
		io___18.ciunit = iounit_1.iout;
		s_wsfe(&io___18);
		do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
		e_wsfe();
		goto L20;
	    }
	}
	if (m == 0) {
	    m = 1;
	    number[0] = next;
	}
    }

/*     get information about a single specified atom */

    if (m == 1) {
	i__ = number[0];
	io___19.ciunit = iounit_1.iout;
	s_wsfe(&io___19);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
	io___20.ciunit = iounit_1.iout;
	s_wsfe(&io___20);
	do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
	e_wsfe();
	io___21.ciunit = iounit_1.iout;
	s_wsfe(&io___21);
	do_fio(&c__1, describe_ref(0, atoms_1.type__[i__ - 1]), (ftnlen)24);
	e_wsfe();
	io___22.ciunit = iounit_1.iout;
	s_wsfe(&io___22);
	do_fio(&c__1, (char *)&atoms_1.type__[i__ - 1], (ftnlen)sizeof(
		integer));
	e_wsfe();
	if (i__ == 1) {
	    io___23.ciunit = iounit_1.iout;
	    s_wsfe(&io___23);
	    e_wsfe();
	} else {
	    io___24.ciunit = iounit_1.iout;
	    s_wsfe(&io___24);
	    e_wsfe();
	    io___25.ciunit = iounit_1.iout;
	    s_wsfe(&io___25);
	    do_fio(&c__1, (char *)&iz_ref(1, i__), (ftnlen)sizeof(integer));
	    i__2 = -i__;
	    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zcoord_1.zbond[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    if (i__ > 2) {
		io___26.ciunit = iounit_1.iout;
		s_wsfe(&io___26);
		do_fio(&c__1, (char *)&iz_ref(2, i__), (ftnlen)sizeof(integer)
			);
		i__2 = -iz_ref(1, i__);
		do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
		i__3 = -i__;
		do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&zcoord_1.zang[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
		if (i__ > 3) {
		    if (iz_ref(4, i__) == 0) {
			io___27.ciunit = iounit_1.iout;
			s_wsfe(&io___27);
			do_fio(&c__1, (char *)&iz_ref(3, i__), (ftnlen)sizeof(
				integer));
			i__2 = -iz_ref(2, i__);
			do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
			i__3 = -iz_ref(1, i__);
			do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
			i__4 = -i__;
			do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&zcoord_1.ztors[i__ - 1], (
				ftnlen)sizeof(doublereal));
			e_wsfe();
		    } else {
			io___28.ciunit = iounit_1.iout;
			s_wsfe(&io___28);
			do_fio(&c__1, (char *)&iz_ref(3, i__), (ftnlen)sizeof(
				integer));
			i__2 = -iz_ref(1, i__);
			do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
			i__3 = -i__;
			do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&zcoord_1.ztors[i__ - 1], (
				ftnlen)sizeof(doublereal));
			e_wsfe();
			io___29.ciunit = iounit_1.iout;
			s_wsfe(&io___29);
			do_fio(&c__1, (char *)&iz_ref(4, i__), (ftnlen)sizeof(
				integer));
			e_wsfe();
		    }
		}
	    }
	}
	next = i__ + 1;
	if (next > atoms_1.n) {
	    next = 1;
	}

/*     chirality change for an atom was requested */

    } else if (m == 2 && number[1] < 0) {
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (iz_ref(4, i__) != 0 && iz_ref(1, i__) == number[0]) {
		changed = TRUE_;
		io___30.ciunit = iounit_1.iout;
		s_wsfe(&io___30);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		e_wsfe();
		iz_ref(4, i__) = -iz_ref(4, i__);
	    }
	}
	next = number[0];
	makexyz_();

/*     information about a specified bond or distance */

    } else if (m == 2) {
	i__ = max(number[0],number[1]);
	j = min(number[0],number[1]);
	if (min(i__,j) <= 0 && max(i__,j) > atoms_1.n) {
	    io___32.ciunit = iounit_1.iout;
	    s_wsfe(&io___32);
	    e_wsfe();
	    error = TRUE_;
	} else {
	    if (j != iz_ref(1, i__)) {
		value = geometry_(&i__, &j, &c__0, &c__0);
		io___34.ciunit = iounit_1.iout;
		s_wsfe(&io___34);
		do_fio(&c__1, (char *)&value, (ftnlen)sizeof(doublereal));
		e_wsfe();
		io___35.ciunit = iounit_1.iout;
		s_wsfe(&io___35);
		e_wsfe();
	    } else {
		io___36.ciunit = iounit_1.iout;
		s_wsfe(&io___36);
		do_fio(&c__1, (char *)&zcoord_1.zbond[i__ - 1], (ftnlen)
			sizeof(doublereal));
		e_wsfe();
		zvalue_("Bond Length", &zcoord_1.zbond[i__ - 1], &changed, (
			ftnlen)11);
		next = i__;
	    }
	}

/*     an atom type change was requested */

    } else if (m == 3 && number[1] < 0) {
	if (number[2] > 0 && number[2] <= 5000) {
	    changed = TRUE_;
	    io___37.ciunit = iounit_1.iout;
	    s_wsfe(&io___37);
	    do_fio(&c__1, describe_ref(0, atoms_1.type__[number[0] - 1]), (
		    ftnlen)24);
	    e_wsfe();
	    atoms_1.type__[number[0] - 1] = number[2];
	    io___38.ciunit = iounit_1.iout;
	    s_wsfe(&io___38);
	    do_fio(&c__1, describe_ref(0, atoms_1.type__[number[0] - 1]), (
		    ftnlen)24);
	    e_wsfe();
	} else {
	    io___39.ciunit = iounit_1.iout;
	    s_wsfe(&io___39);
	    e_wsfe();
	    numrow = 1667;
	    numcol = 2;
	    i__2 = numrow;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (i__ > numrow - 1) {
		    numcol = 1;
		}
		io___42.ciunit = iounit_1.iout;
		s_wsfe(&io___42);
		i__3 = numcol;
		for (j = 0; j <= i__3; ++j) {
		    i__4 = j * numrow + i__;
		    do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
		    do_fio(&c__1, describe_ref(0, j * numrow + i__), (ftnlen)
			    24);
		}
		e_wsfe();
	    }
	}

/*     information about a specified bond angle */

    } else if (m == 3) {
	i__ = max(number[0],number[2]);
	j = number[1];
	k = min(number[0],number[2]);
/* Computing MIN */
	i__2 = min(i__,j);
/* Computing MAX */
	i__4 = max(i__,j);
	if (min(i__2,k) <= 0 || max(i__4,k) > atoms_1.n) {
	    io___44.ciunit = iounit_1.iout;
	    s_wsfe(&io___44);
	    e_wsfe();
	    error = TRUE_;
	} else {
	    if (iz_ref(1, i__) != j) {
		value = geometry_(&i__, &j, &k, &c__0);
		io___45.ciunit = iounit_1.iout;
		s_wsfe(&io___45);
		do_fio(&c__1, (char *)&value, (ftnlen)sizeof(doublereal));
		e_wsfe();
		io___46.ciunit = iounit_1.iout;
		s_wsfe(&io___46);
		e_wsfe();
	    } else if (iz_ref(2, i__) == k) {
		io___47.ciunit = iounit_1.iout;
		s_wsfe(&io___47);
		do_fio(&c__1, (char *)&zcoord_1.zang[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
		zvalue_("Bond Angle", &zcoord_1.zang[i__ - 1], &changed, (
			ftnlen)10);
		next = i__;
	    } else if (iz_ref(3, i__) == k && iz_ref(4, i__) != 0) {
		io___48.ciunit = iounit_1.iout;
		s_wsfe(&io___48);
		do_fio(&c__1, (char *)&zcoord_1.ztors[i__ - 1], (ftnlen)
			sizeof(doublereal));
		e_wsfe();
		zvalue_("Bond Angle", &zcoord_1.ztors[i__ - 1], &changed, (
			ftnlen)10);
		next = i__;
	    } else {
		value = geometry_(&i__, &j, &k, &c__0);
		io___49.ciunit = iounit_1.iout;
		s_wsfe(&io___49);
		do_fio(&c__1, (char *)&value, (ftnlen)sizeof(doublereal));
		e_wsfe();
		io___50.ciunit = iounit_1.iout;
		s_wsfe(&io___50);
		e_wsfe();
	    }
	}

/*    information about a specified dihedral angle */

    } else if (m == 4) {
	if (number[0] > number[3]) {
	    i__ = number[0];
	    j = number[1];
	    k = number[2];
	    l = number[3];
	} else {
	    i__ = number[3];
	    j = number[2];
	    k = number[1];
	    l = number[0];
	}
/* Computing MIN */
	i__2 = min(i__,j), i__2 = min(i__2,k);
/* Computing MAX */
	i__4 = max(i__,j), i__4 = max(i__4,k);
	if (min(i__2,l) <= 0 || max(i__4,l) > atoms_1.n) {
	    io___52.ciunit = iounit_1.iout;
	    s_wsfe(&io___52);
	    e_wsfe();
	    error = TRUE_;
	} else {
	    if (iz_ref(1, i__) != j || iz_ref(2, i__) != k || iz_ref(3, i__) 
		    != l || iz_ref(4, i__) != 0) {
		value = geometry_(&i__, &j, &k, &l);
		io___53.ciunit = iounit_1.iout;
		s_wsfe(&io___53);
		do_fio(&c__1, (char *)&value, (ftnlen)sizeof(doublereal));
		e_wsfe();
		io___54.ciunit = iounit_1.iout;
		s_wsfe(&io___54);
		e_wsfe();
	    } else {
		io___55.ciunit = iounit_1.iout;
		s_wsfe(&io___55);
		do_fio(&c__1, (char *)&zcoord_1.ztors[i__ - 1], (ftnlen)
			sizeof(doublereal));
		e_wsfe();
		zvalue_("Dihedral Angle", &zcoord_1.ztors[i__ - 1], &changed, 
			(ftnlen)14);
		next = i__;
	    }
	}
    }

/*     print instructions for the program if needed */

    if (error) {
	error = FALSE_;
	zhelp_();
    }
    goto L20;

/*     perform any final tasks before program exit */

L410:
    final_();
    return 0;
} /* MAIN__ */

#undef name___ref
#undef iz_ref
#undef describe_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine zhelp  --  print Z-matrix editing instructions  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "zhelp" prints the general information and instructions */
/*     for the Z-matrix editing program */


/* Subroutine */ int zhelp_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 If a single atom number is entered, th"
	    "e\002,\002 current definition of\002,/,\002 the atom will be dis"
	    "played.\002,//,\002 If two atom numbers are entered, the outpu"
	    "t\002,\002 gives the distance\002,/,\002 between the atoms, and "
	    "asks for a new bond\002,\002 length if applicable;\002,/,\002 En"
	    "try of three atoms shows the angle, and\002,\002 entry of four a"
	    "toms\002,/,\002 will display the corresponding dihedral angle"
	    ".\002,//,\002 To change the chirality at an atom, enter\002,\002"
	    " its number and -1.\002,/,\002 To change the type of an atom, en"
	    "ter its\002,\002 number, -1, and the\002,/,\002 new atom type nu"
	    "mber.\002)";
    static char fmt_20[] = "(/,\002 A carriage return at the prompt will dis"
	    "play\002,\002 the atom last\002,/,\002 changed or the next atom "
	    "after the one just\002,\002 examined.\002,//,\002 Typing SHOW wi"
	    "ll display the contents of the\002,\002 current Z-matrix.\002,//,"
	    "\002 Entering EXIT writes a new file then stops,\002,\002 while "
	    "QUIT aborts.\002)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Fortran I/O blocks */
    static cilist io___56 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___57 = { 0, 0, 0, fmt_20, 0 };




/*     print the help and information message for Z-matrix editing */



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


    io___56.ciunit = iounit_1.iout;
    s_wsfe(&io___56);
    e_wsfe();
    io___57.ciunit = iounit_1.iout;
    s_wsfe(&io___57);
    e_wsfe();
    return 0;
} /* zhelp_ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine zvalue  --  gets user input Z-matrix value  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "zvalue" gets user supplied values for selected coordinates */
/*     as needed by the internal coordinate editing program */


/* Subroutine */ int zvalue_(char *text, doublereal *x, logical *changed, 
	ftnlen text_len)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter the New \002,a,\002 :  \002,$)";
    static char fmt_20[] = "(a120)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_rsfe(cilist *), e_rsfe(void), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static doublereal xnew;
    static char record[120];
    static integer length;
    extern /* Subroutine */ int makexyz_(void);

    /* Fortran I/O blocks */
    static cilist io___59 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___60 = { 0, 0, 0, fmt_20, 0 };
    static icilist io___63 = { 1, record, 1, 0, 120, 1 };




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




/*     ask the user for the new internal coordinate value */

    xnew = *x;
    io___59.ciunit = iounit_1.iout;
    s_wsfe(&io___59);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    io___60.ciunit = iounit_1.input;
    s_rsfe(&io___60);
    do_fio(&c__1, record, (ftnlen)120);
    e_rsfe();
    length = trimtext_(record, (ftnlen)120);
    if (length != 0) {
	i__1 = s_rsli(&io___63);
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&xnew, (ftnlen)sizeof(doublereal))
		;
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L30;
	}
L30:
	;
    }

/*     return with the altered value and recompute coordinates */

    if (xnew != *x) {
	*changed = TRUE_;
	*x = xnew;
	makexyz_();
    }
    return 0;
} /* zvalue_ */

/* Main program alias */ int intedit_ () { MAIN__ (); return 0; }
