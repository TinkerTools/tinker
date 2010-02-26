/* prtdyn.f -- translated by f2c (version 20050501).
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
    doublereal xbox, ybox, zbox, alpha, beta, gamma, xbox2, ybox2, zbox2, 
	    box34, lvec[9]	/* was [3][3] */, recip[9]	/* was [3][3] 
	    */, volbox, beta_sin__, beta_cos__, gamma_sin__, gamma_cos__, 
	    beta_term__, gamma_term__;
    logical orthogonal, monoclinic, triclinic, octahedron;
    char spacegrp[10];
} boxes_;

#define boxes_1 boxes_

struct {
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

struct {
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    integer nfree;
    logical velsave, frcsave, uindsave;
    char integrate[10];
} mdstuf_;

#define mdstuf_1 mdstuf_

struct {
    doublereal v[75000]	/* was [3][25000] */, a[75000]	/* was [3][25000] */, 
	    aold[75000]	/* was [3][25000] */;
} moldyn_;

#define moldyn_1 moldyn_

struct {
    doublereal vcm[3000]	/* was [3][1000] */, wcm[3000]	/* was [3][
	    1000] */, lm[3000]	/* was [3][1000] */, vc[3000]	/* was [3][
	    1000] */, wc[3000]	/* was [3][1000] */;
    logical linear[1000];
} rgddyn_;

#define rgddyn_1 rgddyn_

struct {
    integer ltitle;
    char title[120];
} titles_;

#define titles_1 titles_

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine prtdyn  --  output of MD restart information  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "prtdyn" writes out the information needed to restart a */
/*     molecular dynamics trajectory to an external disk file */


/* Subroutine */ int prtdyn_(void)
{
    /* Format strings */
    static char fmt_10[] = "(\002 Number of Atoms and Title :\002)";
    static char fmt_20[] = "(i6)";
    static char fmt_30[] = "(i6,2x,a)";
    static char fmt_40[] = "(\002 Periodic Box Dimensions :\002)";
    static char fmt_50[] = "(3d26.16)";
    static char fmt_60[] = "(3d26.16)";
    static char fmt_70[] = "(\002 Current Atomic Positions :\002)";
    static char fmt_80[] = "(3d26.16)";
    static char fmt_90[] = "(\002 Current Translational Velocities :\002)";
    static char fmt_100[] = "(3d26.16)";
    static char fmt_110[] = "(\002 Current Angular Velocities :\002)";
    static char fmt_120[] = "(3d26.16)";
    static char fmt_130[] = "(\002 Current Angular Momenta :\002)";
    static char fmt_140[] = "(3d26.16)";
    static char fmt_150[] = "(\002 Current Atomic Positions :\002)";
    static char fmt_160[] = "(3d26.16)";
    static char fmt_170[] = "(\002 Current Atomic Velocities :\002)";
    static char fmt_180[] = "(3d26.16)";
    static char fmt_190[] = "(\002 Current Atomic Accelerations :\002)";
    static char fmt_200[] = "(3d26.16)";
    static char fmt_210[] = "(\002 Previous Atomic Accelerations :\002)";
    static char fmt_220[] = "(3d26.16)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_inqu(inlist *), f_open(olist *), f_rew(alist *), s_wsfe(cilist *
	    ), e_wsfe(void), do_fio(integer *, char *, ftnlen), s_cmp(char *, 
	    char *, ftnlen, ftnlen), f_clos(cllist *);

    /* Local variables */
    extern integer freeunit_(void);
    static integer i__, idyn;
    static logical exist;
    static char dynfile[120];

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_220, 0 };



#define a_ref(a_1,a_2) moldyn_1.a[(a_2)*3 + a_1 - 4]
#define v_ref(a_1,a_2) moldyn_1.v[(a_2)*3 + a_1 - 4]
#define lm_ref(a_1,a_2) rgddyn_1.lm[(a_2)*3 + a_1 - 4]
#define vcm_ref(a_1,a_2) rgddyn_1.vcm[(a_2)*3 + a_1 - 4]
#define wcm_ref(a_1,a_2) rgddyn_1.wcm[(a_2)*3 + a_1 - 4]
#define aold_ref(a_1,a_2) moldyn_1.aold[(a_2)*3 + a_1 - 4]



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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  group.i  --  partitioning of system into atom groups  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     grpmass     total mass of all the atoms in each group */
/*     wgrp        weight for each set of group-group interactions */
/*     ngrp        total number of atom groups in the system */
/*     kgrp        contiguous list of the atoms in each group */
/*     igrp        first and last atom of each group in the list */
/*     grplist     number of the group to which each atom belongs */
/*     use_group   flag to use partitioning of system into groups */
/*     use_intra   flag to include only intragroup interactions */
/*     use_inter   flag to include only intergroup interactions */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  mdstuf.i  --  control of molecular dynamics trajectory  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nfree       total number of degrees of freedom for a system */
/*     velsave     flag to save velocity vector components to a file */
/*     frcsave     flag to save force vector components to a file */
/*     uindsave    flag to save induced atomic dipoles to a file */
/*     integrate   type of molecular dynamics integration algorithm */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  moldyn.i  --  velocity and acceleration on MD trajectory  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     v       current velocity of each atom along the x,y,z-axes */
/*     a       current acceleration of each atom along x,y,z-axes */
/*     aold    previous acceleration of each atom along x,y,z-axes */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  rgddyn.i  --  velocities and momenta for rigid body MD  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     vcm     current translational velocity of each rigid body */
/*     wcm     current angular velocity of each rigid body */
/*     lm      current angular momentum of each rigid body */
/*     vc      half-step translational velocity for kinetic energy */
/*     wc      half-step angular velocity for kinetic energy */
/*     linear  logical flag to mark group as linear or nonlinear */




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




/*     update an existing restart file or open a new one */

    idyn = freeunit_();
/* Writing concatenation */
    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
    i__1[1] = 4, a__1[1] = ".dyn";
    s_cat(dynfile, a__1, i__1, &c__2, (ftnlen)120);
    ioin__1.inerr = 0;
    ioin__1.infilen = 120;
    ioin__1.infile = dynfile;
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
	o__1.ounit = idyn;
	o__1.ofnmlen = 120;
	o__1.ofnm = dynfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = idyn;
	f_rew(&al__1);
    } else {
	o__1.oerr = 0;
	o__1.ounit = idyn;
	o__1.ofnmlen = 120;
	o__1.ofnm = dynfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }

/*     save the number of atoms and the title string */

    io___4.ciunit = idyn;
    s_wsfe(&io___4);
    e_wsfe();
    if (titles_1.ltitle == 0) {
	io___5.ciunit = idyn;
	s_wsfe(&io___5);
	do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
	e_wsfe();
    } else {
	io___6.ciunit = idyn;
	s_wsfe(&io___6);
	do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
	do_fio(&c__1, titles_1.title, titles_1.ltitle);
	e_wsfe();
    }

/*     save the periodic box edge lengths and angles */

    io___7.ciunit = idyn;
    s_wsfe(&io___7);
    e_wsfe();
    io___8.ciunit = idyn;
    s_wsfe(&io___8);
    do_fio(&c__1, (char *)&boxes_1.xbox, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.ybox, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.zbox, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___9.ciunit = idyn;
    s_wsfe(&io___9);
    do_fio(&c__1, (char *)&boxes_1.alpha, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.beta, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&boxes_1.gamma, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     save rigid body positions, translational and angular velocities */

    if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 0) {
	io___10.ciunit = idyn;
	s_wsfe(&io___10);
	e_wsfe();
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    io___12.ciunit = idyn;
	    s_wsfe(&io___12);
	    do_fio(&c__1, (char *)&atoms_1.x[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&atoms_1.y[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&atoms_1.z__[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
	io___13.ciunit = idyn;
	s_wsfe(&io___13);
	e_wsfe();
	i__2 = group_1.ngrp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    io___14.ciunit = idyn;
	    s_wsfe(&io___14);
	    do_fio(&c__1, (char *)&vcm_ref(1, i__), (ftnlen)sizeof(doublereal)
		    );
	    do_fio(&c__1, (char *)&vcm_ref(2, i__), (ftnlen)sizeof(doublereal)
		    );
	    do_fio(&c__1, (char *)&vcm_ref(3, i__), (ftnlen)sizeof(doublereal)
		    );
	    e_wsfe();
	}
	io___15.ciunit = idyn;
	s_wsfe(&io___15);
	e_wsfe();
	i__2 = group_1.ngrp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    io___16.ciunit = idyn;
	    s_wsfe(&io___16);
	    do_fio(&c__1, (char *)&wcm_ref(1, i__), (ftnlen)sizeof(doublereal)
		    );
	    do_fio(&c__1, (char *)&wcm_ref(2, i__), (ftnlen)sizeof(doublereal)
		    );
	    do_fio(&c__1, (char *)&wcm_ref(3, i__), (ftnlen)sizeof(doublereal)
		    );
	    e_wsfe();
	}
	io___17.ciunit = idyn;
	s_wsfe(&io___17);
	e_wsfe();
	i__2 = group_1.ngrp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    io___18.ciunit = idyn;
	    s_wsfe(&io___18);
	    do_fio(&c__1, (char *)&lm_ref(1, i__), (ftnlen)sizeof(doublereal))
		    ;
	    do_fio(&c__1, (char *)&lm_ref(2, i__), (ftnlen)sizeof(doublereal))
		    ;
	    do_fio(&c__1, (char *)&lm_ref(3, i__), (ftnlen)sizeof(doublereal))
		    ;
	    e_wsfe();
	}

/*     save the atomic positions, velocities and accelerations */

    } else {
	io___19.ciunit = idyn;
	s_wsfe(&io___19);
	e_wsfe();
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    io___20.ciunit = idyn;
	    s_wsfe(&io___20);
	    do_fio(&c__1, (char *)&atoms_1.x[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&atoms_1.y[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&atoms_1.z__[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
	io___21.ciunit = idyn;
	s_wsfe(&io___21);
	e_wsfe();
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    io___22.ciunit = idyn;
	    s_wsfe(&io___22);
	    do_fio(&c__1, (char *)&v_ref(1, i__), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&v_ref(2, i__), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&v_ref(3, i__), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___23.ciunit = idyn;
	s_wsfe(&io___23);
	e_wsfe();
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    io___24.ciunit = idyn;
	    s_wsfe(&io___24);
	    do_fio(&c__1, (char *)&a_ref(1, i__), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&a_ref(2, i__), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&a_ref(3, i__), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	io___25.ciunit = idyn;
	s_wsfe(&io___25);
	e_wsfe();
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    io___26.ciunit = idyn;
	    s_wsfe(&io___26);
	    do_fio(&c__1, (char *)&aold_ref(1, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&aold_ref(2, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&aold_ref(3, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }

/*     close the dynamics trajectory restart file */

    cl__1.cerr = 0;
    cl__1.cunit = idyn;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* prtdyn_ */

#undef aold_ref
#undef wcm_ref
#undef vcm_ref
#undef lm_ref
#undef v_ref
#undef a_ref


