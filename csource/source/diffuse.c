/* diffuse.f -- translated by f2c (version 20050501).
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
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static integer c__800 = 800;
static integer c__2000 = 2000;



/*     ################################################################# */
/*     ##  COPYRIGHT (C)  1995  by  Yong Kong and Jay William Ponder  ## */
/*     ##                     All Rights Reserved                     ## */
/*     ################################################################# */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  program diffuse  --  find liquid self-diffusion constant  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "diffuse" finds the self-diffusion constant for a homogeneous */
/*     liquid via the Einstein relation from a set of stored molecular */
/*     dynamics frames; molecular centers of mass are unfolded and mean */
/*     squared displacements are computed versus time separation */

/*     the estimate for the self-diffusion constant in 10-5 cm**2/sec */
/*     is printed in the far right column of output and can be checked */
/*     by plotting mean square displacements as a function of the time */
/*     separation */

/*     diffusion values for very large time separation are inaccurate */
/*     due to the small amount of data; the current version requires */
/*     an orthogonal unit cell */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter Coordinate Archive File Name : "
	    " \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_40[] = "(/,\002 Numbers of First & Last Frame and Ste"
	    "p\002,\002 Increment :  \002,$)";
    static char fmt_50[] = "(a120)";
    static char fmt_80[] = "(/,\002 Enter the Time Increment in Picosecond"
	    "s\002,\002 [1.0] :  \002,$)";
    static char fmt_90[] = "(f20.0)";
    static char fmt_100[] = "(/,\002 Enter Unit Cell Axis Lengths :  \002,$)";
    static char fmt_110[] = "(a120)";
    static char fmt_130[] = "(/,\002 DIFFUSE  --  The Maximum of\002,i6,\002"
	    " Molecules\002,\002 has been Exceeded\002)";
    static char fmt_140[] = "(/,\002 Reading the Coordinates Archive File "
	    ":\002,/)";
    static char fmt_150[] = "()";
    static char fmt_170[] = "(4x,\002Processing Coordinate Frame\002,i13)";
    static char fmt_180[] = "(/,\002 DIFFUSE  --  The Maximum of\002,i6,\002"
	    " Frames has been Exceeded\002)";
    static char fmt_200[] = "(/,\002 Total Number of Coordinate Frames :\002"
	    ",i8)";
    static char fmt_210[] = "(/,\002 Mean Squared Diffusion Distance and Sel"
	    "f-Diffusion\002,\002 Constant :\002,//,5x,\002Time Step\002,5x"
	    ",\002X MSD\002,7x,\002Y MSD\002,7x,\002Z MSD\002,7x,\002R MSD"
	    "\002,4x,\002Diff Const\002,/)";
    static char fmt_220[] = "(f12.2,4f12.2,f12.4)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *)
	    , do_fio(integer *, char *, ftnlen), e_rsfe(void), f_open(olist *)
	    , f_rew(alist *), s_rsli(icilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_rsli(void), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int basefile_(char *, ftnlen), molecule_(void), 
	    unitcell_(void);
    extern integer freeunit_(void);
    static integer i__, j, k, m;
    static doublereal dx, dy, dz, xcm[1600000]	/* was [800][2000] */, ycm[
	    1600000]	/* was [800][2000] */, zcm[1600000]	/* was [800][
	    2000] */;
    static integer iarc;
    static doublereal xmid, ymid;
    static integer skip;
    static doublereal zmid, xold, yold, zold;
    static integer step;
    static doublereal xmsd[2000], ymsd[2000], zmsd[2000];
    static integer stop;
    extern /* Subroutine */ int field_(void), fatal_(void);
    static doublereal delta;
    extern /* Subroutine */ int final_(void);
    static doublereal xdiff, ydiff, zdiff, weigh;
    extern /* Subroutine */ int katom_(void);
    static integer ntime[2000];
    static logical exist;
    static integer start;
    static doublereal tstep;
    static logical query;
    static integer iframe, nframe;
    static char record[120];
    static doublereal dvalue, rvalue, xvalue, dunits, yvalue, zvalue, counts;
    static char string[120];
    extern /* Subroutine */ int suffix_(char *, char *, ftnlen, ftnlen);
    static char arcfile[120];
    extern /* Subroutine */ int initial_(void), nextarg_(char *, logical *, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen), readxyz_(
	    integer *);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_20, 0 };
    static icilist io___11 = { 1, string, 1, 0, 120, 1 };
    static icilist io___12 = { 1, string, 1, 0, 120, 1 };
    static icilist io___13 = { 1, string, 1, 0, 120, 1 };
    static cilist io___14 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___17 = { 1, record, 1, 0, 120, 1 };
    static icilist io___19 = { 1, string, 1, 0, 120, 1 };
    static cilist io___20 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_110, 0 };
    static icilist io___24 = { 1, record, 1, 0, 120, 1 };
    static cilist io___25 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___31 = { 1, 0, 1, fmt_150, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___49 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___60 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_220, 0 };



#define xcm_ref(a_1,a_2) xcm[(a_2)*800 + a_1 - 801]
#define ycm_ref(a_1,a_2) ycm[(a_2)*800 + a_1 - 801]
#define zcm_ref(a_1,a_2) zcm[(a_2)*800 + a_1 - 801]
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




/*     perform the standard initialization functions */

    initial_();

/*     try to get a filename from the command line arguments */

    nextarg_(arcfile, &exist, (ftnlen)120);
    if (exist) {
	basefile_(arcfile, (ftnlen)120);
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

/*     ask for the user specified input structure filename */

    while(! exist) {
	io___3.ciunit = iounit_1.iout;
	s_wsfe(&io___3);
	e_wsfe();
	io___4.ciunit = iounit_1.input;
	s_rsfe(&io___4);
	do_fio(&c__1, arcfile, (ftnlen)120);
	e_rsfe();
	basefile_(arcfile, (ftnlen)120);
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

/*     read the first coordinate set in the archive */

    iarc = freeunit_();
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
    readxyz_(&iarc);
    al__1.aerr = 0;
    al__1.aunit = iarc;
    f_rew(&al__1);

/*     get numbers of the coordinate frames to be processed */

    start = 0;
    stop = 0;
    step = 0;
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___11);
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&start, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L30;
	}
	query = FALSE_;
    }
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___12);
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&stop, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L30;
	}
    }
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___13);
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&step, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L30;
	}
    }
L30:
    if (query) {
	io___14.ciunit = iounit_1.iout;
	s_wsfe(&io___14);
	e_wsfe();
	io___15.ciunit = iounit_1.input;
	s_rsfe(&io___15);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___17);
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&start, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&stop, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&step, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L60;
	}
L60:
	;
    }
    if (stop == 0) {
	stop = start;
    }
    if (step == 0) {
	step = 1;
    }

/*     get the time increment between frames in picoseconds */

    tstep = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___19);
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&tstep, (ftnlen)sizeof(doublereal)
		);
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L70;
	}
    }
L70:
    if (tstep <= 0.) {
	io___20.ciunit = iounit_1.iout;
	s_wsfe(&io___20);
	e_wsfe();
	io___21.ciunit = iounit_1.input;
	s_rsfe(&io___21);
	do_fio(&c__1, (char *)&tstep, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (tstep <= 0.) {
	tstep = 1.;
    }

/*     try to find the unit cell axis lengths in the keyfile */

    unitcell_();

/*     get cell axis lengths from command line or interactively */

    while(boxes_1.xbox == 0.) {
	io___22.ciunit = iounit_1.iout;
	s_wsfe(&io___22);
	e_wsfe();
	io___23.ciunit = iounit_1.input;
	s_rsfe(&io___23);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___24);
	if (i__1 != 0) {
	    goto L120;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&boxes_1.xbox, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L120;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&boxes_1.ybox, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L120;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&boxes_1.zbox, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L120;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L120;
	}
L120:
	if (boxes_1.ybox == 0.) {
	    boxes_1.ybox = boxes_1.xbox;
	}
	if (boxes_1.zbox == 0.) {
	    boxes_1.zbox = boxes_1.xbox;
	}
    }

/*     set the half width values for the periodic box */

    boxes_1.xbox2 = boxes_1.xbox * .5;
    boxes_1.ybox2 = boxes_1.ybox * .5;
    boxes_1.zbox2 = boxes_1.zbox * .5;

/*     assign the atom parameters and count the molecules */

    field_();
    katom_();
    molecule_();

/*     check for too many iudividual molecules in the system */

    if (molcul_1.nmol > 800) {
	io___25.ciunit = iounit_1.iout;
	s_wsfe(&io___25);
	do_fio(&c__1, (char *)&c__800, (ftnlen)sizeof(integer));
	e_wsfe();
	fatal_();
    }

/*     get the archived coordinates for each frame in turn */

    io___26.ciunit = iounit_1.iout;
    s_wsfe(&io___26);
    e_wsfe();
    nframe = 0;
    iframe = start;
    skip = start;
    while(iframe >= start && iframe <= stop) {
	skip = (skip - 1) * (atoms_1.n + 1);
	i__1 = skip;
	for (j = 1; j <= i__1; ++j) {
	    io___31.ciunit = iarc;
	    i__2 = s_rsfe(&io___31);
	    if (i__2 != 0) {
		goto L160;
	    }
	    i__2 = e_rsfe();
	    if (i__2 != 0) {
		goto L160;
	    }
	}
L160:
	iframe += step;
	skip = step;
	readxyz_(&iarc);
	if (atoms_1.n == 0) {
	    goto L190;
	}
	++nframe;
	if (nframe % 100 == 0) {
	    io___32.ciunit = iounit_1.iout;
	    s_wsfe(&io___32);
	    do_fio(&c__1, (char *)&nframe, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	if (nframe > 2000) {
	    io___33.ciunit = iounit_1.iout;
	    s_wsfe(&io___33);
	    do_fio(&c__1, (char *)&c__2000, (ftnlen)sizeof(integer));
	    e_wsfe();
	    fatal_();
	}

/*     unfold each molecule to get its corrected center of mass */

	i__1 = molcul_1.nmol;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xmid = 0.;
	    ymid = 0.;
	    zmid = 0.;
	    i__2 = imol_ref(2, i__);
	    for (j = imol_ref(1, i__); j <= i__2; ++j) {
		k = molcul_1.kmol[j - 1];
		weigh = atmtyp_1.mass[k - 1];
		xmid += atoms_1.x[k - 1] * weigh;
		ymid += atoms_1.y[k - 1] * weigh;
		zmid += atoms_1.z__[k - 1] * weigh;
	    }
	    weigh = molcul_1.molmass[i__ - 1];
	    xmid /= weigh;
	    ymid /= weigh;
	    zmid /= weigh;
	    if (nframe == 1) {
		xold = xmid;
		yold = ymid;
		zold = zmid;
	    } else {
		xold = xcm_ref(i__, nframe - 1);
		yold = ycm_ref(i__, nframe - 1);
		zold = zcm_ref(i__, nframe - 1);
	    }
	    dx = xmid - xold;
	    dy = ymid - yold;
	    dz = zmid - zold;
	    while(dx > boxes_1.xbox2) {
		dx -= boxes_1.xbox;
	    }
	    while(dx < -boxes_1.xbox2) {
		dx += boxes_1.xbox;
	    }
	    while(dy > boxes_1.ybox2) {
		dy -= boxes_1.ybox;
	    }
	    while(dy < -boxes_1.ybox2) {
		dy += boxes_1.ybox;
	    }
	    while(dz > boxes_1.zbox2) {
		dz -= boxes_1.zbox;
	    }
	    while(dz < -boxes_1.zbox2) {
		dz += boxes_1.zbox;
	    }
	    xcm_ref(i__, nframe) = xold + dx;
	    ycm_ref(i__, nframe) = yold + dy;
	    zcm_ref(i__, nframe) = zold + dz;
	}
    }
L190:
    cl__1.cerr = 0;
    cl__1.cunit = iarc;
    cl__1.csta = 0;
    f_clos(&cl__1);
    io___49.ciunit = iounit_1.iout;
    s_wsfe(&io___49);
    do_fio(&c__1, (char *)&nframe, (ftnlen)sizeof(integer));
    e_wsfe();

/*     increment the squared displacements for each frame pair */

    i__1 = nframe;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ntime[i__ - 1] = 0;
	xmsd[i__ - 1] = 0.;
	ymsd[i__ - 1] = 0.;
	zmsd[i__ - 1] = 0.;
    }
    i__1 = nframe - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nframe;
	for (j = i__ + 1; j <= i__2; ++j) {
	    m = j - i__;
	    ++ntime[m - 1];
	    i__3 = molcul_1.nmol;
	    for (k = 1; k <= i__3; ++k) {
		xdiff = xcm_ref(k, j) - xcm_ref(k, i__);
		ydiff = ycm_ref(k, j) - ycm_ref(k, i__);
		zdiff = zcm_ref(k, j) - zcm_ref(k, i__);
		xmsd[m - 1] += xdiff * xdiff;
		ymsd[m - 1] += ydiff * ydiff;
		zmsd[m - 1] += zdiff * zdiff;
	    }
	}
    }

/*     get mean squared displacements and convert units; */
/*     conversion is from sq. Ang/ps to 10-5 sq. cm/sec */

    dunits = 10.;
    i__1 = nframe - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	counts = (doublereal) molcul_1.nmol * (doublereal) ntime[i__ - 1];
	xmsd[i__ - 1] *= dunits / counts;
	ymsd[i__ - 1] *= dunits / counts;
	zmsd[i__ - 1] *= dunits / counts;
    }

/*     estimate the diffusion constant via the Einstein relation */

    io___60.ciunit = iounit_1.iout;
    s_wsfe(&io___60);
    e_wsfe();
    i__1 = nframe - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	delta = tstep * (doublereal) i__;
	xvalue = xmsd[i__ - 1] / 2.;
	yvalue = ymsd[i__ - 1] / 2.;
	zvalue = zmsd[i__ - 1] / 2.;
	rvalue = (xmsd[i__ - 1] + ymsd[i__ - 1] + zmsd[i__ - 1]) / 6.;
	dvalue = rvalue / delta;
	io___67.ciunit = iounit_1.iout;
	s_wsfe(&io___67);
	do_fio(&c__1, (char *)&delta, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&xvalue, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&yvalue, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&zvalue, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rvalue, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&dvalue, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef imol_ref
#undef zcm_ref
#undef ycm_ref
#undef xcm_ref


/* Main program alias */ int diffuse_ () { MAIN__ (); return 0; }
