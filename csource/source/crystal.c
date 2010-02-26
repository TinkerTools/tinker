/* crystal.f -- translated by f2c (version 20050501).
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
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

struct {
    doublereal xbox, ybox, zbox, alpha, beta, gamma, xbox2, ybox2, zbox2, 
	    box34, lvec[9]	/* was [3][3] */, recip[9]	/* was [3][3] 
	    */, volbox, beta_sin__, beta_cos__, gamma_sin__, gamma_cos__, 
	    beta_term__, gamma_term__;
    logical orthogonal, monoclinic, triclinic, octahedron;
    char spacegrp[10];
} boxes_;

#define boxes_1 boxes_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

struct {
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__25000 = 25000;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  program crystal  --  fractional coordinate manipulations  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "crystal" is a utility which converts between fractional and */
/*     Cartesian coordinates, and can generate full unit cells from */
/*     asymmetric units */


/* Main program */ int MAIN__(void)
{
    /* Initialized data */

    static char sgroup[10*52] = "P1        " "P2        " "P1(-)     " "P21 "
	    "      " "C2        " "Pc        " "Cc        " "P21/m     " "C2/"
	    "m      " "P2/c      " "P21/c     " "P21/n     " "P21/a     " 
	    "C2/c      " "P21212    " "P212121   " "C2221     " "Pca21     " 
	    "Pmn21     " "Pna21     " "Pn21a     " "Cmc21     " "Aba2      " 
	    "Fdd2      " "Pnna      " "Pmna      " "Pcca      " "Pbam      " 
	    "Pccn      " "Pbcm      " "Pbcn      " "Pbca      " "Pnma      " 
	    "Cmcm      " "Cmca      " "P41       " "P43       " "I4(-)     " 
	    "P42/n     " "I41/a     " "P41212    " "P43212    " "P4(-)21m  " 
	    "P4(-)21c  " "P4(-)m2   " "R3(-)     " "R3c       " "P63/m     " 
	    "P6(3)/mmc " "Pa3(-)    " "Fm3(-)m   " "Im3(-)m   ";

    /* Format strings */
    static char fmt_20[] = "(/,\002 The TINKER Crystal Facility can Provid"
	    "e :\002,//,4x,\002(1) Convert Fractional to Cartesian Coords\002"
	    ",/,4x,\002(2) Convert Cartesian to Fractional Coords\002,/,4x"
	    ",\002(3) Move Any Stray Molecules into Unit Cell\002,/,4x,\002(4"
	    ") Make a Unit Cell from Asymmetric Unit\002,/,4x,\002(5) Make a "
	    "Big Block from Single Unit Cell\002)";
    static char fmt_30[] = "(/,\002 Enter the Number of the Desired Choice :"
	    "  \002,$)";
    static char fmt_40[] = "(i10)";
    static char fmt_70[] = "(/,\002 Available Crystallographic Space Group"
	    "s :\002,/,/,3x,\002(1) \002,a10,3x,\002(2) \002,a10,3x,\002(3)"
	    " \002,a10,3x,\002(4) \002,a10,/,3x,\002(5) \002,a10,3x,\002(6)"
	    " \002,a10,3x,\002(7) \002,a10,3x,\002(8) \002,a10,/,3x,\002(9)"
	    " \002,a10,2x,\002(10) \002,a10,2x,\002(11) \002,a10,2x,\002(12)"
	    " \002,a10,/,2x,\002(13) \002,a10,2x,\002(14) \002,a10,2x,\002(15"
	    ") \002,a10,2x,\002(16) \002,a10,/,2x,\002(17) \002,a10,2x,\002(1"
	    "8) \002,a10,2x,\002(19) \002,a10,2x,\002(20) \002,a10,/,2x,\002("
	    "21) \002,a10,2x,\002(22) \002,a10,2x,\002(23) \002,a10,2x,\002(2"
	    "4) \002,a10,/,2x,\002(25) \002,a10,2x,\002(26) \002,a10,2x,\002("
	    "27) \002,a10,2x,\002(28) \002,a10,/,2x,\002(29) \002,a10,2x,\002"
	    "(30) \002,a10,2x,\002(31) \002,a10,2x,\002(32) \002,a10,/,2x,"
	    "\002(33) \002,a10,2x,\002(34) \002,a10,2x,\002(35) \002,a10,2x"
	    ",\002(36) \002,a10,/,2x,\002(37) \002,a10,2x,\002(38) \002,a10,2"
	    "x,\002(39) \002,a10,2x,\002(40) \002,a10,/,2x,\002(41) \002,a10,"
	    "2x,\002(42) \002,a10,2x,\002(43) \002,a10,2x,\002(44) \002,a10,/"
	    ",2x,\002(45) \002,a10,2x,\002(46) \002,a10,2x,\002(47) \002,a10,"
	    "2x,\002(48) \002,a10,/,2x,\002(49) \002,a10,2x,\002(50) \002,a10"
	    ",2x,\002(51) \002,a10,2x,\002(52) \002,a10)";
    static char fmt_80[] = "(/,\002 Enter the Number of the Desired Choice :"
	    "  \002,$)";
    static char fmt_90[] = "(i10)";
    static char fmt_110[] = "(/,\002 Enter Unit Cell Axis Lengths :  \002,$)";
    static char fmt_120[] = "(a120)";
    static char fmt_140[] = "(/,\002 Enter Unit Cell Axis Angles :   \002,$)";
    static char fmt_150[] = "(a120)";
    static char fmt_170[] = "(/,\002 Unit Cell Dimensions :      a    =\002,"
	    "f10.4,/,\002                             b    =\002,f10.4,/,\002"
	    "                             c    =\002,f10.4,/,\002            "
	    "                Alpha =\002,f10.4,/,\002                        "
	    "    Beta  =\002,f10.4,/,\002                            Gamma "
	    "=\002,f10.4)";
    static char fmt_180[] = "(/,\002 Space Group Symbol :\002,12x,a10)";
    static char fmt_190[] = "(/,\002 Enter Number of Replicates along a-, b-"
	    " and\002,\002 c-Axes [1 1 1] :   \002,$)";
    static char fmt_200[] = "(a120)";
    static char fmt_220[] = "(/,\002 CRYSTAL  --  The Maximum of\002,i8,\002"
	    " Atoms\002,\002 has been Exceeded\002)";
    static char fmt_230[] = "(/,\002 Dimensions of the\002,i3,\002 x\002,i3"
	    ",\002 x\002,i3,\002 Cell Block :\002,//,\002 New Cell Dimensions"
	    " :       a    =\002,f10.4,/,\002                             b  "
	    "  =\002,f10.4,/,\002                             c    =\002,f10."
	    "4)";
    static char fmt_240[] = "(/,\002 Attempt to Merge Fragments to Form Ful"
	    "l\002,\002 Molecules [N] :   \002,$)";
    static char fmt_250[] = "(a120)";
    static char fmt_260[] = "(/,\002 Locate Center of Unit Cell at Coordin"
	    "ate\002,\002 Origin [N] :   \002,$)";
    static char fmt_270[] = "(a120)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2];
    doublereal d__1, d__2, d__3;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void), s_cmp(char *, 
	    char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int bigblock_(integer *, integer *, integer *), 
	    molecule_(void), molmerge_(void), unitcell_(void);
    extern integer freeunit_(void);
    extern /* Subroutine */ int symmetry_(char *, ftnlen);
    static integer i__, na, nb, nc, mode, next, ixyz;
    extern /* Subroutine */ int field_(void), fatal_(void), final_(void), 
	    katom_(void);
    static logical exist, query;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), bounds_(void);
    static char answer[1], string[120];
    extern /* Subroutine */ int getxyz_(void), prtxyz_(integer *), lattice_(
	    void), initial_(void), nextarg_(char *, logical *, ftnlen), 
	    gettext_(char *, char *, integer *, ftnlen, ftnlen), version_(
	    char *, char *, ftnlen, ftnlen);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static icilist io___6 = { 1, string, 1, 0, 120, 1 };
    static cilist io___7 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___9 = { 1, 0, 1, fmt_40, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_120, 0 };
    static icilist io___17 = { 1, record, 1, 0, 120, 1 };
    static cilist io___18 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_150, 0 };
    static icilist io___20 = { 1, record, 1, 0, 120, 1 };
    static cilist io___21 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_200, 0 };
    static icilist io___28 = { 1, record, 1, 0, 120, 1 };
    static cilist io___29 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_270, 0 };



#define sgroup_ref(a_0,a_1) &sgroup[(a_1)*10 + a_0 - 10]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  bound.i  --  control of periodic boundary conditions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     polycut       cutoff distance for infinite polymer nonbonds */
/*     polycut2      square of infinite polymer nonbond cutoff */
/*     use_bounds    flag to use periodic boundary conditions */
/*     use_replica   flag to use replicates for periodic system */
/*     use_polymer   flag to mark presence of infinite polymer */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  boxes.i  --  parameters for periodic boundary conditions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     xbox        length of a-axis of periodic box in Angstroms */
/*     ybox        length of b-axis of periodic box in Angstroms */
/*     zbox        length of c-axis of periodic box in Angstroms */
/*     alpha       angle between b- and c-axes of box in degrees */
/*     beta        angle between a- and c-axes of box in degrees */
/*     gamma       angle between a- and b-axes of box in degrees */
/*     xbox2       half of the a-axis length of periodic box */
/*     ybox2       half of the b-axis length of periodic box */
/*     zbox2       half of the c-axis length of periodic box */
/*     box34       three-fourths axis length of truncated octahedron */
/*     lvec        real space lattice vectors as matrix rows */
/*     recip       reciprocal lattice vectors as matrix columns */
/*     volbox      volume in Ang**3 of the periodic box */
/*     beta_sin    sine of the beta periodic box angle */
/*     beta_cos    cosine of the beta periodic box angle */
/*     gamma_sin   sine of the gamma periodic box angle */
/*     gamma_cos   cosine of the gamma periodic box angle */
/*     beta_term   term used in generating triclinic box */
/*     gamma_term  term used in generating triclinic box */
/*     orthogonal  flag to mark periodic box as orthogonal */
/*     monoclinic  flag to mark periodic box as monoclinic */
/*     triclinic   flag to mark periodic box as triclinic */
/*     octahedron  flag to mark box as truncated octahedron */
/*     spacegrp    space group symbol for the unitcell type */




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




/*     get and read the Cartesian coordinates file */

    initial_();
    getxyz_();

/*     find out which unitcell manipulation is to be performed */

    mode = 0;
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___6);
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
	while(mode < 1 || mode > 5) {
	    mode = 0;
	    io___8.ciunit = iounit_1.iout;
	    s_wsfe(&io___8);
	    e_wsfe();
	    io___9.ciunit = iounit_1.input;
	    i__1 = s_rsfe(&io___9);
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_fio(&c__1, (char *)&mode, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L50;
	    }
L50:
	    ;
	}
    }

/*     get any cell dimensions found in the keyword list */

    unitcell_();

/*     determine the space group if it will be needed later */

    if (mode == 4) {
	for (i__ = 1; i__ <= 52; ++i__) {
	    if (s_cmp(boxes_1.spacegrp, sgroup_ref(0, i__), (ftnlen)10, (
		    ftnlen)10) == 0) {
		goto L100;
	    }
	}
L60:
	io___11.ciunit = iounit_1.iout;
	s_wsfe(&io___11);
	for (i__ = 1; i__ <= 52; ++i__) {
	    do_fio(&c__1, sgroup_ref(0, i__), (ftnlen)10);
	}
	e_wsfe();
	io___12.ciunit = iounit_1.iout;
	s_wsfe(&io___12);
	e_wsfe();
	io___13.ciunit = iounit_1.input;
	s_rsfe(&io___13);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_rsfe();
	if (i__ < 1 || i__ > 52) {
	    goto L60;
	}
	s_copy(boxes_1.spacegrp, sgroup_ref(0, i__), (ftnlen)10, (ftnlen)10);
L100:
	;
    }

/*     if not in keyfile, get the unit cell axis lengths */

    while(boxes_1.xbox == 0.) {
	io___14.ciunit = iounit_1.iout;
	s_wsfe(&io___14);
	e_wsfe();
	io___15.ciunit = iounit_1.input;
	s_rsfe(&io___15);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___17);
	if (i__1 != 0) {
	    goto L130;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&boxes_1.xbox, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L130;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&boxes_1.ybox, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L130;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&boxes_1.zbox, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L130;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L130;
	}
L130:
	if (boxes_1.ybox == 0.) {
	    boxes_1.ybox = boxes_1.xbox;
	}
	if (boxes_1.zbox == 0.) {
	    boxes_1.zbox = boxes_1.xbox;
	}
	bound_1.use_bounds__ = TRUE_;
    }

/*     if not in keyfile, get the unit cell angle values */

    while(boxes_1.alpha == 0.) {
	io___18.ciunit = iounit_1.iout;
	s_wsfe(&io___18);
	e_wsfe();
	io___19.ciunit = iounit_1.input;
	s_rsfe(&io___19);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___20);
	if (i__1 != 0) {
	    goto L160;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&boxes_1.alpha, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L160;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&boxes_1.beta, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L160;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&boxes_1.gamma, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L160;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L160;
	}
L160:
	if (boxes_1.alpha == 0.) {
	    boxes_1.alpha = 90.;
	}
	if (boxes_1.beta == 0.) {
	    boxes_1.beta = boxes_1.alpha;
	}
	if (boxes_1.gamma == 0.) {
	    boxes_1.gamma = boxes_1.alpha;
	}
	if (boxes_1.alpha == 90. && boxes_1.beta == 90. && boxes_1.gamma == 
		90.) {
	    boxes_1.orthogonal = TRUE_;
	} else if (boxes_1.alpha == 90. && boxes_1.gamma == 90.) {
	    boxes_1.monoclinic = TRUE_;
	} else {
	    boxes_1.triclinic = TRUE_;
	}
    }

/*     find constants for coordinate interconversion */

    lattice_();

/*     print out the initial cell dimensions to be used */

    io___21.ciunit = iounit_1.iout;
    s_wsfe(&io___21);
    do_fio(&c__1, (char *)&boxes_1.xbox, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.ybox, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.zbox, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.alpha, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.beta, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.gamma, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     convert Cartesian to fractional coordinates */

    if (mode != 1 && mode != 3) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    atoms_1.z__[i__ - 1] = atoms_1.z__[i__ - 1] / 
		    boxes_1.gamma_term__ / boxes_1.zbox;
	    atoms_1.y[i__ - 1] = (atoms_1.y[i__ - 1] - atoms_1.z__[i__ - 1] * 
		    boxes_1.zbox * boxes_1.beta_term__) / boxes_1.gamma_sin__ 
		    / boxes_1.ybox;
	    atoms_1.x[i__ - 1] = (atoms_1.x[i__ - 1] - atoms_1.y[i__ - 1] * 
		    boxes_1.ybox * boxes_1.gamma_cos__ - atoms_1.z__[i__ - 1] 
		    * boxes_1.zbox * boxes_1.beta_cos__) / boxes_1.xbox;
	}
    }

/*     apply the appropriate space group symmetry operators */

    if (mode == 4) {
	io___22.ciunit = iounit_1.iout;
	s_wsfe(&io___22);
	do_fio(&c__1, boxes_1.spacegrp, (ftnlen)10);
	e_wsfe();
	symmetry_(boxes_1.spacegrp, (ftnlen)10);
    }

/*     replicate the unit cell to make a block of unit cells */

    if (mode == 5) {
	na = 0;
	nb = 0;
	nc = 0;
	io___26.ciunit = iounit_1.iout;
	s_wsfe(&io___26);
	e_wsfe();
	io___27.ciunit = iounit_1.input;
	s_rsfe(&io___27);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___28);
	if (i__1 != 0) {
	    goto L210;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&na, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L210;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&nb, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L210;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&nc, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L210;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L210;
	}
L210:
	if (na == 0) {
	    na = 1;
	}
	if (nb == 0) {
	    nb = na;
	}
	if (nc == 0) {
	    nc = na;
	}
	if (na * nb * nc * atoms_1.n > 25000) {
	    io___29.ciunit = iounit_1.iout;
	    s_wsfe(&io___29);
	    do_fio(&c__1, (char *)&c__25000, (ftnlen)sizeof(integer));
	    e_wsfe();
	    fatal_();
	}
	bigblock_(&na, &nb, &nc);
	io___30.ciunit = iounit_1.iout;
	s_wsfe(&io___30);
	do_fio(&c__1, (char *)&na, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nb, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nc, (ftnlen)sizeof(integer));
	d__1 = (doublereal) na * boxes_1.xbox;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	d__2 = (doublereal) nb * boxes_1.ybox;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	d__3 = (doublereal) nc * boxes_1.zbox;
	do_fio(&c__1, (char *)&d__3, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     convert fractional to Cartesian coordinates */

    if (mode != 2 && mode != 3) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    atoms_1.x[i__ - 1] = atoms_1.x[i__ - 1] * boxes_1.xbox + 
		    atoms_1.y[i__ - 1] * boxes_1.ybox * boxes_1.gamma_cos__ + 
		    atoms_1.z__[i__ - 1] * boxes_1.zbox * boxes_1.beta_cos__;
	    atoms_1.y[i__ - 1] = atoms_1.y[i__ - 1] * boxes_1.ybox * 
		    boxes_1.gamma_sin__ + atoms_1.z__[i__ - 1] * boxes_1.zbox 
		    * boxes_1.beta_term__;
	    atoms_1.z__[i__ - 1] = atoms_1.z__[i__ - 1] * boxes_1.zbox * 
		    boxes_1.gamma_term__;
	}
    }

/*     translate any stray molecules back into the unit cell */

    if (mode == 3) {
	field_();
	katom_();
	molecule_();
	bounds_();
    }

/*     merge fragments to form complete connected molecules */

    if (mode == 4) {
	nextarg_(answer, &exist, (ftnlen)1);
	if (! exist) {
	    io___32.ciunit = iounit_1.iout;
	    s_wsfe(&io___32);
	    e_wsfe();
	    io___33.ciunit = iounit_1.input;
	    s_rsfe(&io___33);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    next = 1;
	    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
	}
	upcase_(answer, (ftnlen)1);
	if (*(unsigned char *)answer == 'Y') {
	    molmerge_();
	}
    }

/*     optionally move unit cell center to coordinate origin */

    if (mode == 4) {
	nextarg_(answer, &exist, (ftnlen)1);
	if (! exist) {
	    io___35.ciunit = iounit_1.iout;
	    s_wsfe(&io___35);
	    e_wsfe();
	    io___36.ciunit = iounit_1.input;
	    s_rsfe(&io___36);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    next = 1;
	    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
	}
	upcase_(answer, (ftnlen)1);
	if (*(unsigned char *)answer == 'Y') {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		atoms_1.z__[i__ - 1] = atoms_1.z__[i__ - 1] / 
			boxes_1.gamma_term__ / boxes_1.zbox;
		atoms_1.y[i__ - 1] = (atoms_1.y[i__ - 1] - atoms_1.z__[i__ - 
			1] * boxes_1.zbox * boxes_1.beta_term__) / 
			boxes_1.gamma_sin__ / boxes_1.ybox;
		atoms_1.x[i__ - 1] = (atoms_1.x[i__ - 1] - atoms_1.y[i__ - 1] 
			* boxes_1.ybox * boxes_1.gamma_cos__ - atoms_1.z__[
			i__ - 1] * boxes_1.zbox * boxes_1.beta_cos__) / 
			boxes_1.xbox;
	    }
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		atoms_1.x[i__ - 1] += -.5;
		atoms_1.y[i__ - 1] += -.5;
		atoms_1.z__[i__ - 1] += -.5;
	    }
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		atoms_1.x[i__ - 1] = atoms_1.x[i__ - 1] * boxes_1.xbox + 
			atoms_1.y[i__ - 1] * boxes_1.ybox * 
			boxes_1.gamma_cos__ + atoms_1.z__[i__ - 1] * 
			boxes_1.zbox * boxes_1.beta_cos__;
		atoms_1.y[i__ - 1] = atoms_1.y[i__ - 1] * boxes_1.ybox * 
			boxes_1.gamma_sin__ + atoms_1.z__[i__ - 1] * 
			boxes_1.zbox * boxes_1.beta_term__;
		atoms_1.z__[i__ - 1] = atoms_1.z__[i__ - 1] * boxes_1.zbox * 
			boxes_1.gamma_term__;
	    }
	}
    }

/*     write out the new coordinates to a file */

    ixyz = freeunit_();
    if (mode == 2) {
/* Writing concatenation */
	i__2[0] = files_1.leng, a__1[0] = files_1.filename;
	i__2[1] = 5, a__1[1] = ".frac";
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
    } else {
/* Writing concatenation */
	i__2[0] = files_1.leng, a__1[0] = files_1.filename;
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
    }
    prtxyz_(&ixyz);
    cl__1.cerr = 0;
    cl__1.cunit = ixyz;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef sgroup_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine molmerge  --  connect fragments into molecules  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "molmerge" connects fragments and removes duplicate atoms */
/*     during generation of a unit cell from an asymmetric unit */


/* Subroutine */ int molmerge_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int molecule_(void);
    static integer h__, i__, j, k, m;
    static doublereal r__;
    static integer im, in, km, kn;
    static doublereal xi, yi, zi, xr, yr, zr, eps;
    static logical join, omit[25000];
    extern /* Subroutine */ int sort_(integer *, integer *), image_(
	    doublereal *, doublereal *, doublereal *);
    static logical merge;
    extern /* Subroutine */ int delete_(integer *);


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define imol_ref(a_1,a_2) molcul_1.imol[(a_2)*2 + a_1 - 3]



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  molcul.i  --  individual molecules within current system  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     molmass   molecular weight for each molecule in the system */
/*     totmass   total weight of all the molecules in the system */
/*     nmol      total number of separate molecules in the system */
/*     kmol      contiguous list of the atoms in each molecule */
/*     imol      first and last atom of each molecule in the list */
/*     molcule   number of the molecule to which each atom belongs */




/*     parse the system to find molecules and fragments */

    molecule_();

/*     zero out the list of atoms to be deleted */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	omit[i__ - 1] = FALSE_;
    }

/*     first pass tests all pairs for duplicate atoms */

    eps = .01;
    i__1 = atoms_1.n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	i__2 = atoms_1.n;
	for (k = i__ + 1; k <= i__2; ++k) {
	    km = molcul_1.molcule[k - 1];
	    xr = atoms_1.x[k - 1] - xi;
	    yr = atoms_1.y[k - 1] - yi;
	    zr = atoms_1.z__[k - 1] - zi;
	    image_(&xr, &yr, &zr);
	    r__ = sqrt(xr * xr + yr * yr + zr * zr);
	    merge = FALSE_;
	    if (r__ < eps) {
		merge = TRUE_;
	    }

/*    translate molecular fragment to the closest image */

	    if (merge) {
		xr = xr - atoms_1.x[k - 1] + xi;
		yr = yr - atoms_1.y[k - 1] + yi;
		zr = zr - atoms_1.z__[k - 1] + zi;
		i__3 = imol_ref(2, km);
		for (j = imol_ref(1, km); j <= i__3; ++j) {
		    m = molcul_1.kmol[j - 1];
		    atoms_1.x[m - 1] += xr;
		    atoms_1.y[m - 1] += yr;
		    atoms_1.z__[m - 1] += zr;
		}

/*    connections between partially duplicated fragments */

		omit[k - 1] = TRUE_;
		i__3 = couple_1.n12[k - 1];
		for (j = 1; j <= i__3; ++j) {
		    m = i12_ref(j, k);
		    join = TRUE_;
		    i__4 = couple_1.n12[m - 1];
		    for (h__ = 1; h__ <= i__4; ++h__) {
			if (i12_ref(h__, m) == i__) {
			    join = FALSE_;
			}
		    }
		    if (join) {
			++couple_1.n12[m - 1];
			i12_ref(couple_1.n12[m - 1], m) = i__;
		    }
		    join = TRUE_;
		    i__4 = couple_1.n12[i__ - 1];
		    for (h__ = 1; h__ <= i__4; ++h__) {
			if (i12_ref(h__, i__) == m) {
			    join = FALSE_;
			}
		    }
		    if (join) {
			++couple_1.n12[i__ - 1];
			i12_ref(couple_1.n12[i__ - 1], i__) = m;
		    }
		}
	    }
	}
    }

/*     delete any duplicated atoms identical by symmetry */

    j = atoms_1.n;
    for (i__ = j; i__ >= 1; --i__) {
	if (omit[i__ - 1]) {
	    delete_(&i__);
	}
    }

/*     parse the system to find molecules and fragments */

    molecule_();

/*     second pass tests all pairs for atoms to be bonded */

    i__1 = atoms_1.n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	im = molcul_1.molcule[i__ - 1];
	in = couple_1.n12[i__ - 1];
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	i__2 = atoms_1.n;
	for (k = i__ + 1; k <= i__2; ++k) {
	    km = molcul_1.molcule[k - 1];
	    kn = couple_1.n12[k - 1];
	    xr = atoms_1.x[k - 1] - xi;
	    yr = atoms_1.y[k - 1] - yi;
	    zr = atoms_1.z__[k - 1] - zi;
	    image_(&xr, &yr, &zr);
	    r__ = sqrt(xr * xr + yr * yr + zr * zr);
	    merge = FALSE_;
	    if (im != km) {
		eps = 1.6;
		if (in > 1 && kn > 1) {
		    eps = 2.;
		}
		if (r__ < eps) {
		    merge = TRUE_;
		}
	    }

/*    translate molecular fragment to the closest image */

	    if (merge) {
		xr = xr - atoms_1.x[k - 1] + xi;
		yr = yr - atoms_1.y[k - 1] + yi;
		zr = zr - atoms_1.z__[k - 1] + zi;
		i__3 = imol_ref(2, km);
		for (j = imol_ref(1, km); j <= i__3; ++j) {
		    m = molcul_1.kmol[j - 1];
		    atoms_1.x[m - 1] += xr;
		    atoms_1.y[m - 1] += yr;
		    atoms_1.z__[m - 1] += zr;
		}

/*     connection between bonded atoms in different fragments */

		couple_1.n12[i__ - 1] = in + 1;
		i12_ref(couple_1.n12[i__ - 1], i__) = k;
		couple_1.n12[k - 1] = kn + 1;
		i12_ref(couple_1.n12[k - 1], k) = i__;
	    }
	}
    }

/*     sort the connected atom lists into ascending order */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sort_(&couple_1.n12[i__ - 1], &i12_ref(1, i__));
    }
    return 0;
} /* molmerge_ */

#undef imol_ref
#undef i12_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine cellatom  --  add new atom to the unit cell  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "cellatom" completes the addition of a symmetry related atom */
/*     to a unit cell by updating the atom type and attachment arrays */


/* Subroutine */ int cellatom_(integer *jj, integer *j)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, delta;


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




/*     attachments of replicated atom are analogous to base atom */

    delta = *jj - *j;
    couple_1.n12[*jj - 1] = couple_1.n12[*j - 1];
    i__1 = couple_1.n12[*j - 1];
    for (i__ = 1; i__ <= i__1; ++i__) {
	i12_ref(i__, *jj) = i12_ref(i__, *j) + delta;
    }
    atoms_1.type__[*jj - 1] = atoms_1.type__[*j - 1];
    s_copy(name___ref(0, *jj), name___ref(0, *j), (ftnlen)3, (ftnlen)3);
    return 0;
} /* cellatom_ */

#undef name___ref
#undef i12_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine bigblock  --  create a block of unit cells  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "bigblock" replicates the coordinates of a single unit */
/*     cell to give a larger block of repeated units */


/* Subroutine */ int bigblock_(integer *na, integer *nb, integer *nc)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int cellatom_(integer *, integer *);
    static integer i__, j, k, ii, jj, nsym;
    static doublereal trans[30000]	/* was [3][10000] */;


#define trans_ref(a_1,a_2) trans[(a_2)*3 + a_1 - 4]



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




/*     construct translation offsets for the replicated cells */

    nsym = 0;
    i__1 = *na / 2;
    for (i__ = (1 - *na) / 2; i__ <= i__1; ++i__) {
	i__2 = *nb / 2;
	for (j = (1 - *nb) / 2; j <= i__2; ++j) {
	    i__3 = *nc / 2;
	    for (k = (1 - *nc) / 2; k <= i__3; ++k) {
		++nsym;
		trans_ref(1, nsym) = (doublereal) i__;
		trans_ref(2, nsym) = (doublereal) j;
		trans_ref(3, nsym) = (doublereal) k;
	    }
	}
    }

/*     put the original cell at the top of the replica list */

    i__1 = nsym;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (trans_ref(1, i__) == 0. && trans_ref(2, i__) == 0. && trans_ref(3,
		 i__) == 0.) {
	    k = i__;
	}
    }
    for (i__ = k; i__ >= 2; --i__) {
	trans_ref(1, i__) = trans_ref(1, i__ - 1);
	trans_ref(2, i__) = trans_ref(2, i__ - 1);
	trans_ref(3, i__) = trans_ref(3, i__ - 1);
    }
    trans_ref(1, 1) = 0.;
    trans_ref(2, 1) = 0.;
    trans_ref(3, 1) = 0.;

/*     translate the original unit cell to make a block of cells */

    i__1 = nsym;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ii = (i__ - 1) * atoms_1.n;
	i__2 = atoms_1.n;
	for (j = 1; j <= i__2; ++j) {
	    jj = j + ii;
	    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + trans_ref(1, i__);
	    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + trans_ref(2, i__);
	    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + trans_ref(3, i__);
	    cellatom_(&jj, &j);
	}
    }
    atoms_1.n = nsym * atoms_1.n;
    return 0;
} /* bigblock_ */

#undef trans_ref




/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine symmetry  --  apply space group symmetry  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "symmetry" applies symmetry operators to the fractional */
/*     coordinates of the asymmetric unit in order to generate */
/*     the symmetry related atoms of the full unit cell */


/* Subroutine */ int symmetry_(char *spacegrp, ftnlen spacegrp_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int cellatom_(integer *, integer *);
    static integer i__, j, ii, jj;
    static doublereal one3, one6, fiv6, two3;
    static integer nsym;



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




/*     P1 space group  (International Tables 1) */

    if (s_cmp(spacegrp, "P1        ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 1;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    atoms_1.x[i__ - 1] = atoms_1.x[i__ - 1];
	    atoms_1.y[i__ - 1] = atoms_1.y[i__ - 1];
	    atoms_1.z__[i__ - 1] = atoms_1.z__[i__ - 1];
	}

/*     P1(-) space group  (International Tables 2) */

    } else if (s_cmp(spacegrp, "P1(-)     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 2;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P2 space group  (International Tables 3) */

    } else if (s_cmp(spacegrp, "P2        ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 2;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P21 space group  (International Tables 4) */

    } else if (s_cmp(spacegrp, "P21       ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 2;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     C2 space group  (International Tables 5) */

    } else if (s_cmp(spacegrp, "C2        ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pc space group  (International Tables 7) */

    } else if (s_cmp(spacegrp, "Pc        ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 2;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Cc space group  (International Tables 9) */

    } else if (s_cmp(spacegrp, "Cc        ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P21/m space group  (International Tables 11) */

    } else if (s_cmp(spacegrp, "P21/m     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     C2/m space group  (International Tables 12) */

    } else if (s_cmp(spacegrp, "C2/m      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P2/c space group  (International Tables 13) */

    } else if (s_cmp(spacegrp, "P2/c      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P21/c space group  (International Tables 14) */

    } else if (s_cmp(spacegrp, "P21/c     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P21/n space group  (International Tables 14) */

    } else if (s_cmp(spacegrp, "P21/n     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P21/a space group  (International Tables 14) */

    } else if (s_cmp(spacegrp, "P21/a     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     C2/c space group  (International Tables 15) */

    } else if (s_cmp(spacegrp, "C2/c      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P21212 space group  (International Tables 18) */

    } else if (s_cmp(spacegrp, "P21212    ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P212121 space group  (International Tables 19) */

    } else if (s_cmp(spacegrp, "P212121   ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     C2221 space group  (International Tables 20) */

    } else if (s_cmp(spacegrp, "C2221     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pca21 space group  (International Tables 29) */

    } else if (s_cmp(spacegrp, "Pca21     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pmn21 space group  (International Tables 31) */

    } else if (s_cmp(spacegrp, "Pmn21     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pna21 space group  (International Tables 33) */

    } else if (s_cmp(spacegrp, "Pna21     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pn21a space group  (International Tables 33) */

    } else if (s_cmp(spacegrp, "Pn21a     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Cmc21 space group  (International Tables 36) */

    } else if (s_cmp(spacegrp, "Cmc21     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Aba2 space group  (International Tables 41) */

    } else if (s_cmp(spacegrp, "Aba2      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Fdd2 space group  (International Tables 43) */

    } else if (s_cmp(spacegrp, "Fdd2      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 16;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .25;
		    atoms_1.y[jj - 1] = .25 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .25;
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = .25 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .25;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .25;
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .25;
		    atoms_1.y[jj - 1] = .75 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .75;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = .25 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .75;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .75;
		} else if (i__ == 9) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 10) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 11) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .75;
		    atoms_1.y[jj - 1] = .25 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .75;
		} else if (i__ == 12) {
		    atoms_1.x[jj - 1] = .75 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .25;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .75;
		} else if (i__ == 13) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 14) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 15) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .75;
		    atoms_1.y[jj - 1] = .75 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .25;
		} else if (i__ == 16) {
		    atoms_1.x[jj - 1] = .75 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .75;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .25;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pnna space group  (International Tables 52) */

    } else if (s_cmp(spacegrp, "Pnna      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pmna space group  (International Tables 53) */

    } else if (s_cmp(spacegrp, "Pmna      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pcca space group  (International Tables 54) */

    } else if (s_cmp(spacegrp, "Pcca      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pbam space group  (International Tables 55) */

    } else if (s_cmp(spacegrp, "Pbam      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pccn space group  (International Tables 56) */

    } else if (s_cmp(spacegrp, "Pccn      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pbcm space group  (International Tables 57) */

    } else if (s_cmp(spacegrp, "Pbcm      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pbcn space group  (International Tables 60) */

    } else if (s_cmp(spacegrp, "Pbcn      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pbca space group  (International Tables 61) */

    } else if (s_cmp(spacegrp, "Pbca      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pnma space group  (International Tables 62) */

    } else if (s_cmp(spacegrp, "Pnma      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Cmcm space group  (International Tables 63) */

    } else if (s_cmp(spacegrp, "Cmcm      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 16;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 9) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 10) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 11) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 12) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 13) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 14) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 15) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 16) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Cmca space group  (International Tables 64) */

    } else if (s_cmp(spacegrp, "Cmca      ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 16;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 9) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 10) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 11) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 12) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 13) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 14) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 15) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 16) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P41 space group  (International Tables 76) */

    } else if (s_cmp(spacegrp, "P41       ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .25;
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .75;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P43 space group  (International Tables 78) */

    } else if (s_cmp(spacegrp, "P43       ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .75;
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .25;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     I4(-) space group  (International Tables 82) */

    } else if (s_cmp(spacegrp, "I4(-)     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P42/n space group  (International Tables 86) */

    } else if (s_cmp(spacegrp, "P42/n     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     I41/a space group  (International Tables 88) */

    } else if (s_cmp(spacegrp, "I41/a     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 16;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = .75 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] + .25;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .25;
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = .25 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] + .75;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .75;
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] + .75;
		    atoms_1.y[jj - 1] = .75 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .75;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] + .25;
		    atoms_1.y[jj - 1] = .25 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .25;
		} else if (i__ == 9) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 10) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 11) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 12) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 13) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] + .25;
		    atoms_1.y[jj - 1] = .75 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = .75 - atoms_1.z__[j - 1];
		} else if (i__ == 14) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] + .75;
		    atoms_1.y[jj - 1] = .25 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = .25 - atoms_1.z__[j - 1];
		} else if (i__ == 15) {
		    atoms_1.x[jj - 1] = .25 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] + .25;
		    atoms_1.z__[jj - 1] = .25 - atoms_1.z__[j - 1];
		} else if (i__ == 16) {
		    atoms_1.x[jj - 1] = .75 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] + .75;
		    atoms_1.z__[jj - 1] = .75 - atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P41212 space group  (International Tables 92) */

    } else if (s_cmp(spacegrp, "P41212    ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .25;
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .75;
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .25 - atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .75 - atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P43212 space group  (International Tables 96) */

    } else if (s_cmp(spacegrp, "P43212    ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .75;
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .25;
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .75 - atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .25 - atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P4(-)21m space group  (International Tables 113) */

    } else if (s_cmp(spacegrp, "P4(-)21m  ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P4(-)21c space group  (International Tables 114) */

    } else if (s_cmp(spacegrp, "P4(-)21c  ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P4(-)m2 space group  (International Tables 115) */

    } else if (s_cmp(spacegrp, "P4(-)m2   ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 8;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     R3(-) space group  (International Tables 148) */

    } else if (s_cmp(spacegrp, "R3(-)     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 18;
	one3 = .33333333333333331;
	two3 = .66666666666666663;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = two3 + atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = one3 + atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = one3 + atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = two3 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = one3 + atoms_1.x[j - 1] - atoms_1.y[j 
			    - 1];
		    atoms_1.z__[jj - 1] = one3 + atoms_1.z__[j - 1];
		} else if (i__ == 9) {
		    atoms_1.x[jj - 1] = two3 + atoms_1.y[j - 1] - atoms_1.x[j 
			    - 1];
		    atoms_1.y[jj - 1] = one3 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = one3 + atoms_1.z__[j - 1];
		} else if (i__ == 10) {
		    atoms_1.x[jj - 1] = two3 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = one3 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = one3 - atoms_1.z__[j - 1];
		} else if (i__ == 11) {
		    atoms_1.x[jj - 1] = two3 + atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = one3 + atoms_1.y[j - 1] - atoms_1.x[j 
			    - 1];
		    atoms_1.z__[jj - 1] = one3 - atoms_1.z__[j - 1];
		} else if (i__ == 12) {
		    atoms_1.x[jj - 1] = two3 + atoms_1.x[j - 1] - atoms_1.y[j 
			    - 1];
		    atoms_1.y[jj - 1] = one3 + atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = one3 - atoms_1.z__[j - 1];
		} else if (i__ == 13) {
		    atoms_1.x[jj - 1] = one3 + atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = two3 + atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = two3 + atoms_1.z__[j - 1];
		} else if (i__ == 14) {
		    atoms_1.x[jj - 1] = one3 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = two3 + atoms_1.x[j - 1] - atoms_1.y[j 
			    - 1];
		    atoms_1.z__[jj - 1] = two3 + atoms_1.z__[j - 1];
		} else if (i__ == 15) {
		    atoms_1.x[jj - 1] = one3 + atoms_1.y[j - 1] - atoms_1.x[j 
			    - 1];
		    atoms_1.y[jj - 1] = two3 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = two3 + atoms_1.z__[j - 1];
		} else if (i__ == 16) {
		    atoms_1.x[jj - 1] = one3 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = two3 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = two3 - atoms_1.z__[j - 1];
		} else if (i__ == 17) {
		    atoms_1.x[jj - 1] = one3 + atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = two3 + atoms_1.y[j - 1] - atoms_1.x[j 
			    - 1];
		    atoms_1.z__[jj - 1] = two3 - atoms_1.z__[j - 1];
		} else if (i__ == 18) {
		    atoms_1.x[jj - 1] = one3 + atoms_1.x[j - 1] - atoms_1.y[j 
			    - 1];
		    atoms_1.y[jj - 1] = two3 + atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = two3 - atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     R3c space group  (International Tables 161) */

    } else if (s_cmp(spacegrp, "R3c       ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 18;
	one3 = .33333333333333331;
	two3 = .66666666666666663;
	one6 = .16666666666666666;
	fiv6 = .83333333333333337;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = two3 + atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = one3 + atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = one3 + atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = two3 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = one3 + atoms_1.x[j - 1] - atoms_1.y[j 
			    - 1];
		    atoms_1.z__[jj - 1] = one3 + atoms_1.z__[j - 1];
		} else if (i__ == 9) {
		    atoms_1.x[jj - 1] = two3 + atoms_1.y[j - 1] - atoms_1.x[j 
			    - 1];
		    atoms_1.y[jj - 1] = one3 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = one3 + atoms_1.z__[j - 1];
		} else if (i__ == 10) {
		    atoms_1.x[jj - 1] = two3 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = one3 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = fiv6 + atoms_1.z__[j - 1];
		} else if (i__ == 11) {
		    atoms_1.x[jj - 1] = two3 + atoms_1.y[j - 1] - atoms_1.x[j 
			    - 1];
		    atoms_1.y[jj - 1] = one3 + atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = fiv6 + atoms_1.z__[j - 1];
		} else if (i__ == 12) {
		    atoms_1.x[jj - 1] = two3 + atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = one3 + atoms_1.x[j - 1] - atoms_1.y[j 
			    - 1];
		    atoms_1.z__[jj - 1] = fiv6 + atoms_1.z__[j - 1];
		} else if (i__ == 13) {
		    atoms_1.x[jj - 1] = one3 + atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = two3 + atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = two3 + atoms_1.z__[j - 1];
		} else if (i__ == 14) {
		    atoms_1.x[jj - 1] = one3 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = two3 + atoms_1.x[j - 1] - atoms_1.y[j 
			    - 1];
		    atoms_1.z__[jj - 1] = two3 + atoms_1.z__[j - 1];
		} else if (i__ == 15) {
		    atoms_1.x[jj - 1] = one3 + atoms_1.y[j - 1] - atoms_1.x[j 
			    - 1];
		    atoms_1.y[jj - 1] = two3 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = two3 + atoms_1.z__[j - 1];
		} else if (i__ == 16) {
		    atoms_1.x[jj - 1] = one3 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = two3 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = one6 + atoms_1.z__[j - 1];
		} else if (i__ == 17) {
		    atoms_1.x[jj - 1] = one3 + atoms_1.y[j - 1] - atoms_1.x[j 
			    - 1];
		    atoms_1.y[jj - 1] = two3 + atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = one6 + atoms_1.z__[j - 1];
		} else if (i__ == 18) {
		    atoms_1.x[jj - 1] = one3 + atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = two3 + atoms_1.x[j - 1] - atoms_1.y[j 
			    - 1];
		    atoms_1.z__[jj - 1] = one6 + atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P63/m space group  (International Tables 176) */

    } else if (s_cmp(spacegrp, "P63/m     ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 12;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 9) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 10) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 11) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 12) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     P6(3)/mmc space group  (Intl. Tables 194, Hexagonal Close Packed) */

    } else if (s_cmp(spacegrp, "P6(3)/mmc ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 2;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Pa3(-) space group  (International Tables 205) */

    } else if (s_cmp(spacegrp, "Pa3(-)    ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 24;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 5) {
		    atoms_1.x[jj - 1] = atoms_1.z__[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.y[j - 1];
		} else if (i__ == 6) {
		    atoms_1.x[jj - 1] = atoms_1.z__[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.y[j - 1];
		} else if (i__ == 7) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.z__[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.y[j - 1] + .5;
		} else if (i__ == 8) {
		    atoms_1.x[jj - 1] = -atoms_1.z__[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.y[j - 1];
		} else if (i__ == 9) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.z__[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.x[j - 1];
		} else if (i__ == 10) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.z__[j - 1] + .5;
		    atoms_1.z__[jj - 1] = .5 - atoms_1.x[j - 1];
		} else if (i__ == 11) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.y[jj - 1] = .5 - atoms_1.z__[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.x[j - 1];
		} else if (i__ == 12) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.z__[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.x[j - 1] + .5;
		} else if (i__ == 13) {
		    atoms_1.x[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.z__[j - 1];
		} else if (i__ == 14) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.z__[j - 1];
		} else if (i__ == 15) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 16) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		} else if (i__ == 17) {
		    atoms_1.x[jj - 1] = -atoms_1.z__[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.y[j - 1];
		} else if (i__ == 18) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.z__[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.y[j - 1];
		} else if (i__ == 19) {
		    atoms_1.x[jj - 1] = atoms_1.z__[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.y[j - 1];
		} else if (i__ == 20) {
		    atoms_1.x[jj - 1] = atoms_1.z__[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.x[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.y[j - 1] + .5;
		} else if (i__ == 21) {
		    atoms_1.x[jj - 1] = -atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = -atoms_1.z__[j - 1];
		    atoms_1.z__[jj - 1] = -atoms_1.x[j - 1];
		} else if (i__ == 22) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = .5 - atoms_1.z__[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.x[j - 1] + .5;
		} else if (i__ == 23) {
		    atoms_1.x[jj - 1] = .5 - atoms_1.y[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.z__[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.x[j - 1];
		} else if (i__ == 24) {
		    atoms_1.x[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.z__[j - 1];
		    atoms_1.z__[jj - 1] = .5 - atoms_1.x[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Fm3m space group  (Intl. Tables 225, Face Centered Cubic) */

    } else if (s_cmp(spacegrp, "Fm3(-)m   ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 4;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1];
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 3) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1];
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		} else if (i__ == 4) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1];
		}
		cellatom_(&jj, &j);
	    }
	}

/*     Im3m space group  (Intl. Tables 229, Body Centered Cubic) */

    } else if (s_cmp(spacegrp, "Im3(-)m   ", (ftnlen)10, (ftnlen)10) == 0) {
	nsym = 2;
	i__1 = nsym;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ii = (i__ - 1) * atoms_1.n;
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + ii;
		if (i__ == 2) {
		    atoms_1.x[jj - 1] = atoms_1.x[j - 1] + .5;
		    atoms_1.y[jj - 1] = atoms_1.y[j - 1] + .5;
		    atoms_1.z__[jj - 1] = atoms_1.z__[j - 1] + .5;
		}
		cellatom_(&jj, &j);
	    }
	}
    }

/*     set the total number of atoms in the full unitcell */

    atoms_1.n = nsym * atoms_1.n;
    return 0;
} /* symmetry_ */

/* Main program alias */ int crystal_ () { MAIN__ (); return 0; }
