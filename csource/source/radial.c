/* radial.f -- translated by f2c (version 20050501).
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
    integer narg;
    logical listarg[21];
    char arg[2520];
} argue_;

#define argue_1 argue_

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
    doublereal vdwcut, chgcut, dplcut, mpolecut, vdwtaper, chgtaper, dpltaper,
	     mpoletaper, ewaldcut;
    logical use_ewald__, use_lights__, use_list__, use_vlist__, use_clist__, 
	    use_mlist__;
} cutoff_;

#define cutoff_1 cutoff_

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
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

struct {
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;



/*     ################################################################# */
/*     ##  COPYRIGHT (C)  1995  by  Yong Kong and Jay William Ponder  ## */
/*     ##                     All Rights Reserved                     ## */
/*     ################################################################# */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  program radial  --  compute radial distribution function  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "radial" finds the radial distribution function for a specified */
/*     pair of atom types via analysis of a set of coordinate frames */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter Coordinate Archive File Name : "
	    " \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_40[] = "(/,\002 Numbers of First & Last Frame and Ste"
	    "p\002,\002 Increment :  \002,$)";
    static char fmt_50[] = "(a120)";
    static char fmt_70[] = "(/,\002 Enter 1st & 2nd Atom Names or Type Numbe"
	    "rs :  \002,$)";
    static char fmt_80[] = "(a120)";
    static char fmt_120[] = "(/,\002 Enter Maximum Distance to Accumulate"
	    "\002,\002 [10.0 Ang] :  \002,$)";
    static char fmt_130[] = "(f20.0)";
    static char fmt_150[] = "(/,\002 Enter Width of Distance Bins [0.01 Ang]"
	    " :  \002,$)";
    static char fmt_160[] = "(f20.0)";
    static char fmt_170[] = "(/,\002 Include Intramolecular Pairs in Distrib"
	    "utions\002,\002 [N] :  \002,$)";
    static char fmt_180[] = "(a120)";
    static char fmt_190[] = "(/,\002 RADIAL  --  Too many Distance Bins; Inc"
	    "rease\002,\002 MAXBIN\002)";
    static char fmt_200[] = "(/,\002 Reading the Coordinates Archive File "
	    ":\002,/)";
    static char fmt_210[] = "()";
    static char fmt_230[] = "(4x,\002Processing Coordinate Frame\002,i13)";
    static char fmt_240[] = "(/,\002 Total Number of Coordinate Frames :\002"
	    ",i8)";
    static char fmt_250[] = "(/,\002 Pairwise Radial Distribution Function "
	    ":\002//,7x,\002First Name or Type :  \002,a6,5x,\002Second Name "
	    "or Type :  \002,a6)";
    static char fmt_260[] = "(/,5x,\002Bin\002,9x,\002Counts\002,7x,\002Dist"
	    "ance\002,7x,\002Raw g(r)\002,4x,\002Smooth g(r)\002,/)";
    static char fmt_270[] = "(i8,i15,3x,f12.4,3x,f12.4,3x,f12.4)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *)
	    , do_fio(integer *, char *, ftnlen), e_rsfe(void), f_open(olist *)
	    , f_rew(alist *), s_rsli(icilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_rsli(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int basefile_(char *, ftnlen), molecule_(void), 
	    unitcell_(void);
    extern integer freeunit_(void);
    static logical intramol;
    static integer i__, j, k;
    static doublereal gr[10000], gs[10000], dx, dy, dz, xj, yj, zj;
    static integer bin;
    static doublereal rjk;
    static integer iarc, nbin, molj, molk, skip, hist[10000];
    static doublereal rmax;
    static integer numj, numk, step, next, stop;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *), fatal_(void), final_(void);
    static char namej[3], namek[3];
    static doublereal width, pairs;
    static integer typej, typek, start;
    static logical exist, query;
    static char labelj[6], labelk[6];
    static integer iframe, nframe;
    static doublereal factor, volume, rlower, rupper, expect;
    static char answer[1], record[120], string[120];
    extern /* Subroutine */ int suffix_(char *, char *, ftnlen, ftnlen), 
	    upcase_(char *, ftnlen);
    static char arcfile[120];
    extern /* Subroutine */ int lattice_(void), initial_(void), nextarg_(char 
	    *, logical *, ftnlen), gettext_(char *, char *, integer *, ftnlen,
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
    static cilist io___20 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_80, 0 };
    static icilist io___25 = { 1, labelj, 1, 0, 6, 1 };
    static icilist io___28 = { 1, labelk, 1, 0, 6, 1 };
    static icilist io___30 = { 1, string, 1, 0, 120, 1 };
    static cilist io___31 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_130, 0 };
    static icilist io___34 = { 1, string, 1, 0, 120, 1 };
    static cilist io___35 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___52 = { 1, 0, 1, fmt_210, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___65 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___75 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___76 = { 0, 0, 0, fmt_270, 0 };



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
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  argue.i  --  command line arguments at program startup  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     maxarg    maximum number of command line arguments */

/*     narg      number of command line arguments to the program */
/*     listarg   flag to mark available command line arguments */
/*     arg       strings containing the command line arguments */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cutoff.i  --  cutoff distances for energy interactions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     vdwcut      cutoff distance for van der Waals interactions */
/*     chgcut      cutoff distance for charge-charge interactions */
/*     dplcut      cutoff distance for dipole-dipole interactions */
/*     mpolecut    cutoff distance for atomic multipole interactions */
/*     vdwtaper    distance at which van der Waals switching begins */
/*     chgtaper    distance at which charge-charge switching begins */
/*     dpltaper    distance at which dipole-dipole switching begins */
/*     mpoletaper  distance at which atomic multipole switching begins */
/*     ewaldcut    cutoff distance for direct space Ewald summation */
/*     use_ewald   logical flag governing use of Ewald summation */
/*     use_lights  logical flag governing use of method of lights */
/*     use_list    logical flag governing use of any neighbor lists */
/*     use_vlist   logical flag governing use of vdw neighbor lists */
/*     use_clist   logical flag governing use of charge neighbor lists */
/*     use_mlist   logical flag governing use of multipole neighbor lists */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  potent.i  --  usage of each potential energy component  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     use_bond    logical flag governing use of bond stretch potential */
/*     use_angle   logical flag governing use of angle bend potential */
/*     use_strbnd  logical flag governing use of stretch-bend potential */
/*     use_urey    logical flag governing use of Urey-Bradley potential */
/*     use_angang  logical flag governing use of angle-angle cross term */
/*     use_opbend  logical flag governing use of out-of-plane bend term */
/*     use_opdist  logical flag governing use of out-of-plane distance */
/*     use_improp  logical flag governing use of improper dihedral term */
/*     use_imptor  logical flag governing use of improper torsion term */
/*     use_tors    logical flag governing use of torsional potential */
/*     use_pitors  logical flag governing use of pi-orbital torsion term */
/*     use_strtor  logical flag governing use of stretch-torsion term */
/*     use_tortor  logical flag governing use of torsion-torsion term */
/*     use_vdw     logical flag governing use of vdw der Waals potential */
/*     use_charge  logical flag governing use of charge-charge potential */
/*     use_chgdpl  logical flag governing use of charge-dipole potential */
/*     use_dipole  logical flag governing use of dipole-dipole potential */
/*     use_mpole   logical flag governing use of multipole potential */
/*     use_polar   logical flag governing use of polarization term */
/*     use_rxnfld  logical flag governing use of reaction field term */
/*     use_solv    logical flag governing use of continuum solvation */
/*     use_metal   logical flag governing use of ligand field term */
/*     use_geom    logical flag governing use of geometric restraints */
/*     use_extra   logical flag governing use of extra potential term */
/*     use_born    logical flag governing use of Born radii values */
/*     use_orbit   logical flag governing use of pisystem computation */




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

/*     get the unitcell parameters and number of molecules */

    unitcell_();
    molecule_();

/*     set cutoffs small to enforce use of minimum images */

    potent_1.use_vdw__ = TRUE_;
    potent_1.use_charge__ = FALSE_;
    potent_1.use_dipole__ = FALSE_;
    potent_1.use_mpole__ = FALSE_;
    cutoff_1.use_ewald__ = FALSE_;
    cutoff_1.vdwcut = .01;
    lattice_();

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

/*     get the names of the atoms to be used in rdf computation */

    nextarg_(labelj, &exist, (ftnlen)6);
    nextarg_(labelk, &exist, (ftnlen)6);
    if (! exist) {
	io___20.ciunit = iounit_1.iout;
	s_wsfe(&io___20);
	e_wsfe();
	io___21.ciunit = iounit_1.input;
	s_rsfe(&io___21);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, labelj, &next, (ftnlen)120, (ftnlen)6);
	gettext_(record, labelk, &next, (ftnlen)120, (ftnlen)6);
    }

/*     convert the labels to either atom names or type numbers */

    s_copy(namej, "   ", (ftnlen)3, (ftnlen)3);
    typej = -1;
    i__1 = s_rsli(&io___25);
    if (i__1 != 0) {
	goto L90;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&typej, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L90;
    }
    i__1 = e_rsli();
    if (i__1 != 0) {
	goto L90;
    }
L90:
    if (typej < 0) {
	next = 1;
	gettext_(labelj, namej, &next, (ftnlen)6, (ftnlen)3);
    }
    s_copy(namek, "   ", (ftnlen)3, (ftnlen)3);
    typek = -1;
    i__1 = s_rsli(&io___28);
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&typek, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = e_rsli();
    if (i__1 != 0) {
	goto L100;
    }
L100:
    if (typek < 0) {
	next = 1;
	gettext_(labelk, namek, &next, (ftnlen)6, (ftnlen)3);
    }

/*     get maximum distance from input or minimum image convention */

    if (! bound_1.use_bounds__) {
	rmax = -1.;
	query = TRUE_;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___30);
	    if (i__1 != 0) {
		goto L110;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&rmax, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L110;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L110;
	    }
	    query = FALSE_;
	}
L110:
	if (query) {
	    io___31.ciunit = iounit_1.iout;
	    s_wsfe(&io___31);
	    e_wsfe();
	    io___32.ciunit = iounit_1.input;
	    s_rsfe(&io___32);
	    do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(doublereal));
	    e_rsfe();
	}
	if (rmax <= 0.) {
	    rmax = 10.;
	}
    } else if (boxes_1.octahedron) {
	rmax = sqrt(3.) / 4. * boxes_1.xbox;
    } else {
/* Computing MIN */
	d__1 = boxes_1.xbox2 * boxes_1.beta_sin__ * boxes_1.gamma_sin__, d__2 
		= boxes_1.ybox2 * boxes_1.gamma_sin__, d__1 = min(d__1,d__2), 
		d__2 = boxes_1.zbox2 * boxes_1.beta_sin__;
	rmax = min(d__1,d__2);
    }

/*     get the desired width of the radial distance bins */

    width = -1.;
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___34);
	if (i__1 != 0) {
	    goto L140;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&width, (ftnlen)sizeof(doublereal)
		);
	if (i__1 != 0) {
	    goto L140;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L140;
	}
	query = FALSE_;
    }
L140:
    if (query) {
	io___35.ciunit = iounit_1.iout;
	s_wsfe(&io___35);
	e_wsfe();
	io___36.ciunit = iounit_1.input;
	s_rsfe(&io___36);
	do_fio(&c__1, (char *)&width, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (width <= 0.) {
	width = .01;
    }

/*     decide whether to restrict to intermolecular atom pairs */

    intramol = FALSE_;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___39.ciunit = iounit_1.iout;
	s_wsfe(&io___39);
	e_wsfe();
	io___40.ciunit = iounit_1.input;
	s_rsfe(&io___40);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'Y') {
	intramol = TRUE_;
    }

/*     set the number of distance bins to be accumulated */

    nbin = (integer) (rmax / width);
    if (nbin > 10000) {
	io___42.ciunit = iounit_1.iout;
	s_wsfe(&io___42);
	e_wsfe();
	fatal_();
    }

/*     zero out the distance bins and distribution functions */

    i__1 = nbin;
    for (i__ = 1; i__ <= i__1; ++i__) {
	hist[i__ - 1] = 0;
	gr[i__ - 1] = 0.;
	gs[i__ - 1] = 0.;
    }

/*     get the archived coordinates for each frame in turn */

    io___47.ciunit = iounit_1.iout;
    s_wsfe(&io___47);
    e_wsfe();
    nframe = 0;
    iframe = start;
    skip = start;
    while(iframe >= start && iframe <= stop) {
	skip = (skip - 1) * (atoms_1.n + 1);
	i__1 = skip;
	for (j = 1; j <= i__1; ++j) {
	    io___52.ciunit = iarc;
	    i__2 = s_rsfe(&io___52);
	    if (i__2 != 0) {
		goto L220;
	    }
	    i__2 = e_rsfe();
	    if (i__2 != 0) {
		goto L220;
	    }
	}
L220:
	iframe += step;
	skip = step;
	readxyz_(&iarc);
	if (! inform_1.abort) {
	    ++nframe;
	    if (nframe % 100 == 0) {
		io___53.ciunit = iounit_1.iout;
		s_wsfe(&io___53);
		do_fio(&c__1, (char *)&nframe, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	    i__1 = atoms_1.n;
	    for (j = 1; j <= i__1; ++j) {
		if (s_cmp(name___ref(0, j), namej, (ftnlen)3, (ftnlen)3) == 0 
			|| atoms_1.type__[j - 1] == typej) {
		    xj = atoms_1.x[j - 1];
		    yj = atoms_1.y[j - 1];
		    zj = atoms_1.z__[j - 1];
		    molj = molcul_1.molcule[j - 1];
		    i__2 = atoms_1.n;
		    for (k = 1; k <= i__2; ++k) {
			if (s_cmp(name___ref(0, k), namek, (ftnlen)3, (ftnlen)
				3) == 0 || atoms_1.type__[k - 1] == typek) {
			    molk = molcul_1.molcule[k - 1];
			    if (intramol || molj != molk) {
				dx = atoms_1.x[k - 1] - xj;
				dy = atoms_1.y[k - 1] - yj;
				dz = atoms_1.z__[k - 1] - zj;
				image_(&dx, &dy, &dz);
				rjk = sqrt(dx * dx + dy * dy + dz * dz);
				bin = (integer) (rjk / width) + 1;
				++hist[bin - 1];
			    }
			}
		    }
		}
	    }
	}
    }

/*     ensure a valid frame is loaded and report total frames */

    if (inform_1.abort) {
	al__1.aerr = 0;
	al__1.aunit = iarc;
	f_rew(&al__1);
	readxyz_(&iarc);
    }
    cl__1.cerr = 0;
    cl__1.cunit = iarc;
    cl__1.csta = 0;
    f_clos(&cl__1);
    io___65.ciunit = iounit_1.iout;
    s_wsfe(&io___65);
    do_fio(&c__1, (char *)&nframe, (ftnlen)sizeof(integer));
    e_wsfe();

/*     count the number of occurences of each atom type */

    numj = 0;
    numk = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (s_cmp(name___ref(0, i__), namej, (ftnlen)3, (ftnlen)3) == 0 || 
		atoms_1.type__[i__ - 1] == typej) {
	    ++numj;
	}
	if (s_cmp(name___ref(0, i__), namek, (ftnlen)3, (ftnlen)3) == 0 || 
		atoms_1.type__[i__ - 1] == typek) {
	    ++numk;
	}
    }

/*     normalize the distance bins to give radial distribution */

    if (numj != 0 && numk != 0) {
	factor = (doublereal) nframe * 4.1887902047863905;
	if (bound_1.use_bounds__) {
	    pairs = (doublereal) numj * (doublereal) numk;
	    volume = boxes_1.gamma_sin__ * boxes_1.gamma_term__ * 
		    boxes_1.xbox * boxes_1.ybox * boxes_1.zbox;
	    if (boxes_1.octahedron) {
		volume *= .5;
	    }
	    factor = factor * pairs / volume;
	}
	i__1 = nbin;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    rupper = (doublereal) i__ * width;
	    rlower = rupper - width;
/* Computing 3rd power */
	    d__1 = rupper;
/* Computing 3rd power */
	    d__2 = rlower;
	    expect = factor * (d__1 * (d__1 * d__1) - d__2 * (d__2 * d__2));
	    gr[i__ - 1] = (doublereal) hist[i__ - 1] / expect;
	}
    }

/*     find the 5th degree polynomial smoothed distribution function */

    if (nbin >= 5) {
	gs[0] = (gr[0] * 69. + gr[1] * 4. - gr[2] * 6. + gr[3] * 4. - gr[4]) /
		 70.;
	gs[1] = (gr[0] * 2. + gr[1] * 27. + gr[2] * 12. - gr[3] * 8. + gr[4] *
		 2.) / 35.;
	i__1 = nbin - 2;
	for (i__ = 3; i__ <= i__1; ++i__) {
	    gs[i__ - 1] = (gr[i__ - 3] * -3. + gr[i__ - 2] * 12. + gr[i__ - 1]
		     * 17. + gr[i__] * 12. - gr[i__ + 1] * 3.) / 35.;
	}
	gs[nbin - 2] = (gr[nbin - 5] * 2. - gr[nbin - 4] * 8. + gr[nbin - 3] *
		 12. + gr[nbin - 2] * 27. + gr[nbin - 1] * 2.) / 35.;
	gs[nbin - 1] = (-gr[nbin - 5] + gr[nbin - 4] * 4. - gr[nbin - 3] * 6. 
		+ gr[nbin - 2] * 4. + gr[nbin - 1] * 69.) / 70.;
	i__1 = nbin;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__1 = 0., d__2 = gs[i__ - 1];
	    gs[i__ - 1] = max(d__1,d__2);
	}
    }

/*     output the final radial distribution function results */

    io___74.ciunit = iounit_1.iout;
    s_wsfe(&io___74);
    do_fio(&c__1, labelj, (ftnlen)6);
    do_fio(&c__1, labelk, (ftnlen)6);
    e_wsfe();
    io___75.ciunit = iounit_1.iout;
    s_wsfe(&io___75);
    e_wsfe();
    i__1 = nbin;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___76.ciunit = iounit_1.iout;
	s_wsfe(&io___76);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&hist[i__ - 1], (ftnlen)sizeof(integer));
	d__1 = ((doublereal) i__ - .5) * width;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&gr[i__ - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&gs[i__ - 1], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef name___ref


/* Main program alias */ int radial_ () { MAIN__ (); return 0; }
