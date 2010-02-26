/* potential.f -- translated by f2c (version 20050501).
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
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

struct {
    doublereal fctmin, hguess;
    integer maxiter, nextiter;
} minima_;

#define minima_1 minima_

struct {
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

struct {
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

struct {
    doublereal pgrid[3000000]	/* was [3][100000][10] */, epot[2000000]	
	    /* was [2][100000][10] */, xdpl0[10], ydpl0[10], zdpl0[10], 
	    xxqdp0[10], xyqdp0[10], xzqdp0[10], yyqdp0[10], yzqdp0[10], 
	    zzqdp0[10];
    integer nconf, npgrid[10], ipgrid[1000000]	/* was [100000][10] */, ngatm,
	     nfatm;
    logical gatm[25000], fatm[25000], use_dpl__, use_qdp__, fit_mpl__, 
	    fit_dpl__, fit_qdp__;
} potfit_;

#define potfit_1 potfit_

struct {
    doublereal pchg[25000];
    integer nion, iion[25000], jion[25000], kion[25000], chglist[25000];
} charge_;

#define charge_1 charge_

struct {
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

struct {
    integer ltitle;
    char title[120];
} titles_;

#define titles_1 titles_

struct {
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

struct {
    doublereal electric, dielec, ebuffer, c2scale, c3scale, c4scale, c5scale;
    logical neutnbr, neutcut;
} chgpot_;

#define chgpot_1 chgpot_

struct {
    doublereal bdpl[50000], sdpl[50000];
    integer ndipole, idpl[100000]	/* was [2][50000] */;
} dipole_;

#define dipole_1 dipole_

struct {
    doublereal polarity[25000], thole[25000], pdamp[25000], uind[75000]	/* 
	    was [3][25000] */, uinp[75000]	/* was [3][25000] */, uinds[
	    75000]	/* was [3][25000] */, uinps[75000]	/* was [3][
	    25000] */;
    integer npolar;
} polar_;

#define polar_1 polar_

struct {
    doublereal netchg, netdpl, netqdp[3], xdpl, ydpl, zdpl, xxqdp, xyqdp, 
	    xzqdp, yxqdp, yyqdp, yzqdp, zxqdp, zyqdp, zzqdp;
} moment_;

#define moment_1 moment_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__0 = 0;
static doublereal c_b170 = 0.;
static integer c__2 = 2;



/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  program potential  --  compute electrostatic potential  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "potential" calculates the electrostatic potential for a */
/*     molecule at a set of grid points; optionally compares to a */
/*     target potential or optimizes electrostatic parameters */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 The TINKER Electrostatic Potential Facil"
	    "ity Can :\002,//,4x,\002(1) Create an Input File for Gaussian CU"
	    "BEGEN\002,/,4x,\002(2) Get QM Potential from a Gaussian CUBE File"
	    "\002,/,4x,\002(3) Calculate the Model Potential for a System\002"
	    ",/,4x,\002(4) Compare the Model Potentials of Two Systems\002,/,"
	    "4x,\002(5) Compare a Model Potential to a Target Grid\002,/,4x"
	    ",\002(6) Fit Electrostatic Parameters to Target Grid\002)";
    static char fmt_30[] = "(/,\002 Enter the Number of the Desired Choice :"
	    "  \002,$)";
    static char fmt_40[] = "(i10)";
    static char fmt_60[] = "(/,\002 Enter the Gaussian CUBE File Name :  "
	    "\002,$)";
    static char fmt_70[] = "(a120)";
    static char fmt_80[] = "(1x,a120)";
    static char fmt_90[] = "()";
    static char fmt_100[] = "(i5)";
    static char fmt_110[] = "(i5)";
    static char fmt_120[] = "()";
    static char fmt_130[] = "(a120)";
    static char fmt_140[] = "(i8,2x,a)";
    static char fmt_150[] = "(i8,3x,3f12.6,2x,f12.4)";
    static char fmt_160[] = "(/,\002 Electrostatic Potential Written to File"
	    " :  \002,a)";
    static char fmt_170[] = "(/,\002 Conformations for Potential Analysis "
	    ":\002,i6)";
    static char fmt_220[] = "(/,\002 Enter Name of Second Coordinate File "
	    ":  \002,$)";
    static char fmt_230[] = "(a120)";
    static char fmt_240[] = "(/,\002 Enter Target Grid/Potential File Name :"
	    "  \002,$)";
    static char fmt_250[] = "(a120)";
    static char fmt_260[] = "(/,\002 Output Potential Value at Each Grid Poi"
	    "nt\002,\002 [N] :  \002,$)";
    static char fmt_270[] = "(a120)";
    static char fmt_280[] = "()";
    static char fmt_290[] = "(\002 POTENTIAL  --  Too many Grid Points;\002"
	    ",\002 Increase MAXGRID\002)";
    static char fmt_300[] = "(\002 Electrostatic Potential Grid Points :\002"
	    ",6x,i10)";
    static char fmt_310[] = "(\002 Potential Grid Points for Conformer\002,i"
	    "4,\002 :\002,2x,i10)";
    static char fmt_320[] = "(3f15.8)";
    static char fmt_330[] = "(/,\002 Gaussian CUBEGEN Input Written to File "
	    ":   \002,a)";
    static char fmt_340[] = "(/,\002 Next, run the Gaussian CUBEGEN program;"
	    " for\002,\002 example:\002,/,5x,\002cubegen 0 potential=MP2 xxx."
	    "fchk\002,\002 xxx.cube -5 h < xxx.grid\002,//,\002 See the Gauss"
	    "ian documentation for additional\002,\002 details;\002,/,\002 Af"
	    "ter CUBEGEN, rerun TINKER POTENTIAL using\002,\002 Option 2\002)";
    static char fmt_360[] = "(/,\002 Enter RMS Gradient Termination Criter"
	    "ion\002,\002 [0.001] :  \002,$)";
    static char fmt_370[] = "(f20.0)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void), f_inqu(inlist *),
	     f_open(olist *), f_rew(alist *), f_clos(cllist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int basefile_(char *, ftnlen);
    static char gridfile[120];
    static logical dotarget;
    extern integer freeunit_(void);
    extern /* Subroutine */ int potpoint_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern integer trimtext_(char *, ftnlen);
    extern /* Subroutine */ int torsions_(void);
    static integer i__, j, k;
    static doublereal x0, y0, z0, xi, yi, zi, xx[75000], xx0, xy0, xz0, yy0, 
	    yz0, zz0, pot;
    static integer mode, igrd;
    extern /* Subroutine */ int ocvm_(integer *, doublereal *, doublereal *, 
	    doublereal *, D_fp, U_fp);
    static integer nvar, ipot, next, ixyz;
    extern /* Subroutine */ int field_(void), fatal_(void);
    static logical dofit;
    extern /* Subroutine */ int bonds_(void), katom_(void);
    static integer glist[25000], flist[25000];
    extern /* Subroutine */ int rings_(void);
    static logical exist, query, docube;
    extern /* Subroutine */ int attach_(void), induce_(void);
    static logical dogrid;
    extern /* Subroutine */ int angles_(void), active_(void), getref_(integer 
	    *);
    static integer nmodel, nglist, nflist;
    static doublereal grdmin;
    static logical dopair, dofull;
    static char answer[1], record[120], string[120];
    extern /* Subroutine */ int suffix_(char *, char *, ftnlen, ftnlen), 
	    bitors_(void), kmpole_(void), kpolar_(void), getxyz_(void), 
	    upcase_(char *, ftnlen), prmvar_(integer *, doublereal *), 
	    varprm_(integer *, doublereal *, integer *, doublereal *), 
	    prtfit_(void);
    extern doublereal potfit1_();
    extern /* Subroutine */ int kcharge_(void), makeref_(integer *);
    static logical domodel;
    extern /* Subroutine */ int kdipole_(void), chkpole_(void), initial_(void)
	    , readpot_(integer *, integer *);
    static char potfile[120];
    static doublereal minimum;
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen);
    extern /* Subroutine */ int optsave_();
    extern /* Subroutine */ int cutoffs_(void), potgrid_(integer *);
    static char keyword[20];
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen), 
	    readxyz_(integer *), gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), rotpole_(void);
    static char xyzfile[120];
    extern /* Subroutine */ int potstat_(logical *, logical *, logical *, 
	    logical *);

    /* Fortran I/O blocks */
    static icilist io___13 = { 1, string, 1, 0, 120, 1 };
    static cilist io___14 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___16 = { 1, 0, 1, fmt_40, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_130, 0 };
    static icilist io___28 = { 0, record, 0, 0, 120, 1 };
    static cilist io___33 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_170, 0 };
    static icilist io___45 = { 1, string, 1, 0, 120, 1 };
    static icilist io___47 = { 1, string, 1, 0, 120, 1 };
    static icilist io___48 = { 1, string, 1, 0, 120, 1 };
    static icilist io___52 = { 1, string, 1, 0, 120, 1 };
    static cilist io___59 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___60 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___65 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___66 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___68 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___69 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_320, 0 };
    static cilist io___75 = { 0, 0, 0, fmt_330, 0 };
    static cilist io___76 = { 0, 0, 0, fmt_340, 0 };
    static icilist io___78 = { 1, string, 1, 0, 120, 1 };
    static cilist io___79 = { 0, 0, 0, fmt_360, 0 };
    static cilist io___80 = { 0, 0, 0, fmt_370, 0 };



#define epot_ref(a_1,a_2,a_3) potfit_1.epot[((a_3)*100000 + (a_2))*2 + \
a_1 - 200003]
#define pgrid_ref(a_1,a_2,a_3) potfit_1.pgrid[((a_3)*100000 + (a_2))*3 + \
a_1 - 300004]
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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  keys.i  --  contents of current keyword parameter file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nkey      number of nonblank lines in the keyword file */
/*     keyline   contents of each individual keyword file line */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  minima.i  --  general parameters for minimizations  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     fctmin    value below which function is deemed optimized */
/*     hguess    initial value for the H-matrix diagonal elements */
/*     maxiter   maximum number of iterations during optimization */
/*     nextiter  iteration number to use for the first iteration */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  output.i  --  control of coordinate output file format  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     archive    logical flag to save structures in an archive */
/*     noversion  logical flag governing use of filename versions */
/*     overwrite  logical flag to overwrite intermediate files inplace */
/*     cyclesave  logical flag to mark use of numbered cycle files */
/*     coordtype  selects Cartesian, internal, rigid body or none */




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




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  potfit.i  --  values for electrostatic potential fitting  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxpgrd   maximum dimension of electrostatic potential grid */

/*     pgrid     Cartesian coordinates of potential grid points */
/*     epot      values of electrostatic potential at grid points */
/*     xdpl0     target x-component of the molecular dipole moment */
/*     ydpl0     target y-component of the molecular dipole moment */
/*     zdpl0     target z-component of the molecular dipole moment */
/*     xxqdp0    target xx-component of the molecular quadrupole moment */
/*     xyqdp0    target xy-component of the molecular quadrupole moment */
/*     xzqdp0    target xz-component of the molecular quadrupole moment */
/*     yyqdp0    target yy-component of the molecular quadrupole moment */
/*     yzqdp0    target yz-component of the molecular quadrupole moment */
/*     zzqdp0    target zz-component of the molecular quadrupole moment */
/*     nconf     total number of conformations to be analyzed */
/*     npgrid    total number of electrostatic potential grid points */
/*     ipgrid    atom associated with each potential grid point */
/*     ngatm     total number of atoms with active potential grid points */
/*     nfatm     total number of atoms in electrostatic potential fit */
/*     gatm      flag to use potential grid points around each atom */
/*     fatm      flag to use each atom in electrostatic potential fit */
/*     use_dpl   flag to include molecular dipole in potential fit */
/*     use_qdp   flag to include molecular quadrupole in potential fit */
/*     fit_mpl   flag for atomic monopoles to vary in potential fit */
/*     fit_dpl   flag for atomic dipoles to vary in potential fit */
/*     fit_qdp   flag for atomic quadrupoles to vary in potential fit */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  units.i  --  physical constants and unit conversions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     literature reference: */

/*     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended */
/*     Values of the Fundamental Physical Constants: 2006", Reviews of */
/*     Modern Physics, 80, 633-730 (2008) */

/*     The "2006 CODATA Recommended Values" are also available from */
/*     the NIST Reference on Constants, Units, and Uncertainty site */
/*     at http://physics.nist.gov/cuu/index.html */

/*     Most values below are derived from the 2006 CODATA reference */

/*     The conversion from calorie to Joule is the definition of the */
/*     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992) */

/*     The "coulomb" energy conversion factor is found by dimensional */
/*     analysis of Coulomb's Law, ie, by dividing the square of the */
/*     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is */
/*     the permittivity of vacuum (the "electric constant"); note that */
/*     eps0 is typically given in F/m, equivalent to C**2/(J-m) */

/*     The approximate value used for the Debye, 3.33564 x 10-30 C-m, */
/*     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997) */

/*     The value of "prescon" is based on definition of 1 atmosphere */
/*     as 101325 Pa set by the 10th Conference Generale des Poids et */
/*     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3 */

/*     avogadro    Avogadro's number (N) in particles/mole */
/*     lightspd    speed of light in vacuum (c) in cm/ps */
/*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K */
/*     gasconst    ideal gas constant (R) in kcal/mole/K */
/*     emass       mass of an electron in atomic mass units */
/*     joule       conversion from calories to joules */
/*     convert     conversion from kcal to g*Ang**2/ps**2 */
/*     bohr        conversion from Bohrs to Angstroms */
/*     hartree     conversion from Hartree to kcal/mole */
/*     evolt       conversion from Hartree to electron-volts */
/*     efreq       conversion from Hartree to cm-1 */
/*     coulomb     conversion from electron**2/Ang to kcal/mole */
/*     debye       conversion from electron-Ang to Debyes */
/*     prescon     conversion from kcal/mole/Ang**3 to Atm */




/*     setup the computation and assign some default values */

    initial_();
    nmodel = 1;
    potfit_1.nconf = 0;
    dogrid = FALSE_;
    docube = FALSE_;
    domodel = FALSE_;
    dopair = FALSE_;
    dotarget = FALSE_;
    dofit = FALSE_;
    potfit_1.fit_mpl__ = TRUE_;
    potfit_1.fit_dpl__ = TRUE_;
    potfit_1.fit_qdp__ = TRUE_;

/*     initialize target molecular dipole and quadrupole values */

    potfit_1.use_dpl__ = FALSE_;
    potfit_1.use_qdp__ = FALSE_;
    for (i__ = 1; i__ <= 10; ++i__) {
	potfit_1.xdpl0[i__ - 1] = 0.;
	potfit_1.ydpl0[i__ - 1] = 0.;
	potfit_1.zdpl0[i__ - 1] = 0.;
	potfit_1.xxqdp0[i__ - 1] = 0.;
	potfit_1.xyqdp0[i__ - 1] = 0.;
	potfit_1.xzqdp0[i__ - 1] = 0.;
	potfit_1.yyqdp0[i__ - 1] = 0.;
	potfit_1.yzqdp0[i__ - 1] = 0.;
	potfit_1.zzqdp0[i__ - 1] = 0.;
    }

/*     find electrostatic potential manipulation to be performed */

    mode = 0;
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___13);
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
	io___14.ciunit = iounit_1.iout;
	s_wsfe(&io___14);
	e_wsfe();
	while(mode < 1 || mode > 6) {
	    mode = 0;
	    io___15.ciunit = iounit_1.iout;
	    s_wsfe(&io___15);
	    e_wsfe();
	    io___16.ciunit = iounit_1.input;
	    i__1 = s_rsfe(&io___16);
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
    if (mode == 1) {
	dogrid = TRUE_;
    } else if (mode == 2) {
	docube = TRUE_;
    } else if (mode == 3) {
	domodel = TRUE_;
    } else if (mode == 4) {
	nmodel = 2;
	dopair = TRUE_;
    } else if (mode == 5) {
	dotarget = TRUE_;
    } else if (mode == 6) {
	dotarget = TRUE_;
	dofit = TRUE_;
    }

/*     read the electrostatic potential from a Gaussian CUBE file */

    if (docube) {
	nextarg_(potfile, &exist, (ftnlen)120);
	if (exist) {
	    basefile_(potfile, (ftnlen)120);
	    suffix_(potfile, "pot", (ftnlen)120, (ftnlen)3);
	    version_(potfile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = potfile;
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
	while(! exist) {
	    io___18.ciunit = iounit_1.iout;
	    s_wsfe(&io___18);
	    e_wsfe();
	    io___19.ciunit = iounit_1.input;
	    s_rsfe(&io___19);
	    do_fio(&c__1, potfile, (ftnlen)120);
	    e_rsfe();
	    basefile_(potfile, (ftnlen)120);
	    suffix_(potfile, "cube", (ftnlen)120, (ftnlen)4);
	    version_(potfile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = potfile;
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
	ipot = freeunit_();
	o__1.oerr = 0;
	o__1.ounit = ipot;
	o__1.ofnmlen = 120;
	o__1.ofnm = potfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = ipot;
	f_rew(&al__1);
	io___21.ciunit = ipot;
	s_rsfe(&io___21);
	do_fio(&c__1, titles_1.title, (ftnlen)120);
	e_rsfe();
	titles_1.ltitle = trimtext_(titles_1.title, (ftnlen)120);
	io___22.ciunit = ipot;
	s_rsfe(&io___22);
	e_rsfe();
	io___23.ciunit = ipot;
	s_rsfe(&io___23);
	do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
	e_rsfe();
	io___24.ciunit = ipot;
	s_rsfe(&io___24);
	do_fio(&c__1, (char *)&potfit_1.npgrid[0], (ftnlen)sizeof(integer));
	e_rsfe();
	i__1 = atoms_1.n + 2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___25.ciunit = ipot;
	    s_rsfe(&io___25);
	    e_rsfe();
	}
	i__1 = potfit_1.npgrid[0];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___26.ciunit = ipot;
	    s_rsfe(&io___26);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    s_rsli(&io___28);
	    do_lio(&c__5, &c__1, (char *)&xi, (ftnlen)sizeof(doublereal));
	    do_lio(&c__5, &c__1, (char *)&yi, (ftnlen)sizeof(doublereal));
	    do_lio(&c__5, &c__1, (char *)&zi, (ftnlen)sizeof(doublereal));
	    do_lio(&c__5, &c__1, (char *)&pot, (ftnlen)sizeof(doublereal));
	    e_rsli();
	    pgrid_ref(1, i__, 1) = xi;
	    pgrid_ref(2, i__, 1) = yi;
	    pgrid_ref(3, i__, 1) = zi;
	    epot_ref(1, i__, 1) = pot * 627.5094688;
	}
	cl__1.cerr = 0;
	cl__1.cunit = ipot;
	cl__1.csta = 0;
	f_clos(&cl__1);
	s_copy(potfile, files_1.filename, (ftnlen)120, (ftnlen)120);
	suffix_(potfile, "pot", (ftnlen)120, (ftnlen)3);
	version_(potfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = ipot;
	o__1.ofnmlen = 120;
	o__1.ofnm = potfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = ipot;
	f_rew(&al__1);
	io___33.ciunit = ipot;
	s_wsfe(&io___33);
	do_fio(&c__1, (char *)&potfit_1.npgrid[0], (ftnlen)sizeof(integer));
	do_fio(&c__1, titles_1.title, titles_1.ltitle);
	e_wsfe();
	i__1 = potfit_1.npgrid[0];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xi = pgrid_ref(1, i__, 1);
	    yi = pgrid_ref(2, i__, 1);
	    zi = pgrid_ref(3, i__, 1);
	    pot = epot_ref(1, i__, 1);
	    io___34.ciunit = ipot;
	    s_wsfe(&io___34);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xi, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&yi, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&zi, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&pot, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	cl__1.cerr = 0;
	cl__1.cunit = ipot;
	cl__1.csta = 0;
	f_clos(&cl__1);
	io___35.ciunit = iounit_1.iout;
	s_wsfe(&io___35);
	do_fio(&c__1, potfile, trimtext_(potfile, (ftnlen)120));
	e_wsfe();
	goto L380;
    }

/*     read the first structure and get connectivity info */

    getxyz_();
    attach_();
    active_();
    bonds_();
    angles_();
    torsions_();
    bitors_();
    rings_();

/*     reopen the structure file and get all conformations */

    ixyz = freeunit_();
    s_copy(xyzfile, files_1.filename, (ftnlen)120, (ftnlen)120);
    suffix_(xyzfile, "xyz", (ftnlen)120, (ftnlen)3);
    version_(xyzfile, "old", (ftnlen)120, (ftnlen)3);
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
    while(! inform_1.abort) {
	++potfit_1.nconf;
	makeref_(&potfit_1.nconf);
	readxyz_(&ixyz);
    }
    cl__1.cerr = 0;
    cl__1.cunit = ixyz;
    cl__1.csta = 0;
    f_clos(&cl__1);
    getref_(&c__1);
    if (potfit_1.nconf > 1) {
	io___38.ciunit = iounit_1.iout;
	s_wsfe(&io___38);
	do_fio(&c__1, (char *)&potfit_1.nconf, (ftnlen)sizeof(integer));
	e_wsfe();
    }

/*     get the force field parameters and setup electrostatics */

    cutoffs_();
    field_();
    katom_();
    if (potent_1.use_charge__ || potent_1.use_chgdpl__) {
	kcharge_();
    }
    if (potent_1.use_dipole__ || potent_1.use_chgdpl__) {
	kdipole_();
    }
    if (potent_1.use_mpole__ || potent_1.use_polar__) {
	kmpole_();
    }
    if (potent_1.use_polar__) {
	kpolar_();
    }

/*     set defaults for the active grid atoms and fit atoms */

    nglist = 0;
    nflist = 0;
    potfit_1.ngatm = atoms_1.n;
    potfit_1.nfatm = atoms_1.n;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	glist[i__ - 1] = 0;
	flist[i__ - 1] = 0;
	potfit_1.gatm[i__ - 1] = TRUE_;
	potfit_1.fatm[i__ - 1] = TRUE_;
    }

/*     get control parameters and target values from keyfile */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "POTENTIAL-ATOMS ", (ftnlen)16, (ftnlen)16) == 0) {
	    i__2 = s_rsli(&io___45);
	    if (i__2 != 0) {
		goto L180;
	    }
	    i__3 = atoms_1.n;
	    for (k = nglist + 1; k <= i__3; ++k) {
		i__2 = do_lio(&c__3, &c__1, (char *)&glist[k - 1], (ftnlen)
			sizeof(integer));
		if (i__2 != 0) {
		    goto L180;
		}
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L180;
	    }
L180:
	    while(glist[nglist] != 0) {
		++nglist;
/* Computing MAX */
/* Computing MIN */
		i__4 = atoms_1.n, i__5 = glist[nglist - 1];
		i__2 = -atoms_1.n, i__3 = min(i__4,i__5);
		glist[nglist - 1] = max(i__2,i__3);
	    }
	} else if (s_cmp(keyword, "POTENTIAL-FIT ", (ftnlen)14, (ftnlen)14) ==
		 0) {
	    i__2 = s_rsli(&io___47);
	    if (i__2 != 0) {
		goto L190;
	    }
	    i__3 = atoms_1.n;
	    for (k = nflist + 1; k <= i__3; ++k) {
		i__2 = do_lio(&c__3, &c__1, (char *)&flist[k - 1], (ftnlen)
			sizeof(integer));
		if (i__2 != 0) {
		    goto L190;
		}
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L190;
	    }
L190:
	    while(flist[nflist] != 0) {
		++nflist;
/* Computing MAX */
/* Computing MIN */
		i__4 = atoms_1.n, i__5 = flist[nflist - 1];
		i__2 = -atoms_1.n, i__3 = min(i__4,i__5);
		flist[nflist - 1] = max(i__2,i__3);
	    }
	} else if (s_cmp(keyword, "FIX-MONOPOLE ", (ftnlen)13, (ftnlen)13) == 
		0) {
	    potfit_1.fit_mpl__ = FALSE_;
	} else if (s_cmp(keyword, "FIX-DIPOLE ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    potfit_1.fit_dpl__ = FALSE_;
	} else if (s_cmp(keyword, "FIX-QUADRUPOLE ", (ftnlen)15, (ftnlen)15) 
		== 0) {
	    potfit_1.fit_qdp__ = FALSE_;
	} else if (s_cmp(keyword, "TARGET-DIPOLE ", (ftnlen)14, (ftnlen)14) ==
		 0) {
	    potfit_1.use_dpl__ = TRUE_;
	    k = 1;
	    i__2 = s_rsli(&io___48);
	    if (i__2 != 0) {
		goto L200;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&x0, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L200;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&y0, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L200;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&z0, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L200;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L200;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L200;
	    }
L200:
	    potfit_1.xdpl0[k - 1] = x0;
	    potfit_1.ydpl0[k - 1] = y0;
	    potfit_1.zdpl0[k - 1] = z0;
	} else if (s_cmp(keyword, "TARGET-QUADRUPOLE ", (ftnlen)18, (ftnlen)
		18) == 0) {
	    potfit_1.use_qdp__ = TRUE_;
	    k = 1;
	    i__2 = s_rsli(&io___52);
	    if (i__2 != 0) {
		goto L210;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&xx0, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L210;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&xy0, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L210;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&xz0, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L210;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&yy0, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L210;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&yz0, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L210;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&zz0, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L210;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L210;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L210;
	    }
L210:
	    potfit_1.xxqdp0[k - 1] = xx0;
	    potfit_1.xyqdp0[k - 1] = xy0;
	    potfit_1.xzqdp0[k - 1] = xz0;
	    potfit_1.yyqdp0[k - 1] = yy0;
	    potfit_1.yzqdp0[k - 1] = yz0;
	    potfit_1.zzqdp0[k - 1] = zz0;
	}
    }

/*     set active grid atoms to only those marked for use */

    i__ = 1;
    while(glist[i__ - 1] != 0) {
	if (i__ == 1) {
	    potfit_1.ngatm = 0;
	    i__1 = atoms_1.n;
	    for (k = 1; k <= i__1; ++k) {
		potfit_1.gatm[k - 1] = FALSE_;
	    }
	}
	if (glist[i__ - 1] > 0) {
	    k = glist[i__ - 1];
	    if (! potfit_1.gatm[k - 1]) {
		potfit_1.gatm[k - 1] = TRUE_;
		++potfit_1.ngatm;
	    }
	    ++i__;
	} else {
	    i__3 = (i__2 = glist[i__], abs(i__2));
	    for (k = (i__1 = glist[i__ - 1], abs(i__1)); k <= i__3; ++k) {
		if (! potfit_1.gatm[k - 1]) {
		    potfit_1.gatm[k - 1] = TRUE_;
		    ++potfit_1.ngatm;
		}
	    }
	    i__ += 2;
	}
    }

/*     set active fitting atoms to only those marked for use */

    i__ = 1;
    while(flist[i__ - 1] != 0) {
	if (i__ == 1) {
	    potfit_1.nfatm = 0;
	    i__3 = atoms_1.n;
	    for (k = 1; k <= i__3; ++k) {
		potfit_1.fatm[k - 1] = FALSE_;
	    }
	}
	if (flist[i__ - 1] > 0) {
	    k = flist[i__ - 1];
	    if (! potfit_1.fatm[k - 1]) {
		potfit_1.fatm[k - 1] = TRUE_;
		++potfit_1.nfatm;
	    }
	    ++i__;
	} else {
	    i__2 = (i__1 = flist[i__], abs(i__1));
	    for (k = (i__3 = flist[i__ - 1], abs(i__3)); k <= i__2; ++k) {
		if (! potfit_1.fatm[k - 1]) {
		    potfit_1.fatm[k - 1] = TRUE_;
		    ++potfit_1.nfatm;
		}
	    }
	    i__ += 2;
	}
    }

/*     generate potential grid based on the molecular surface */

    if (! dotarget) {
	i__2 = potfit_1.nconf;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    getref_(&i__);
	    potgrid_(&i__);
	}
    }

/*     get name of optional second structure for comparison */

    if (dopair) {
	nextarg_(xyzfile, &exist, (ftnlen)120);
	if (exist) {
	    basefile_(xyzfile, (ftnlen)120);
	    suffix_(xyzfile, "xyz", (ftnlen)120, (ftnlen)3);
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
	while(! exist) {
	    io___59.ciunit = iounit_1.iout;
	    s_wsfe(&io___59);
	    e_wsfe();
	    io___60.ciunit = iounit_1.input;
	    s_rsfe(&io___60);
	    do_fio(&c__1, xyzfile, (ftnlen)120);
	    e_rsfe();
	    basefile_(xyzfile, (ftnlen)120);
	    suffix_(xyzfile, "xyz", (ftnlen)120, (ftnlen)3);
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
    }

/*     get optional file with grid points and target potential */

    if (dotarget) {
	nextarg_(potfile, &exist, (ftnlen)120);
	if (exist) {
	    basefile_(potfile, (ftnlen)120);
	    suffix_(potfile, "pot", (ftnlen)120, (ftnlen)3);
	    version_(potfile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = potfile;
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
	while(! exist) {
	    io___61.ciunit = iounit_1.iout;
	    s_wsfe(&io___61);
	    e_wsfe();
	    io___62.ciunit = iounit_1.input;
	    s_rsfe(&io___62);
	    do_fio(&c__1, potfile, (ftnlen)120);
	    e_rsfe();
	    basefile_(potfile, (ftnlen)120);
	    suffix_(potfile, "pot", (ftnlen)120, (ftnlen)3);
	    version_(potfile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = potfile;
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
    }

/*     decide whether to output potential at each grid point */

    dofull = FALSE_;
    if (domodel || dopair || dotarget) {
	nextarg_(answer, &exist, (ftnlen)1);
	if (! exist) {
	    io___65.ciunit = iounit_1.iout;
	    s_wsfe(&io___65);
	    e_wsfe();
	    io___66.ciunit = iounit_1.input;
	    s_rsfe(&io___66);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    next = 1;
	    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
	}
	upcase_(answer, (ftnlen)1);
	if (*(unsigned char *)answer == 'Y') {
	    dofull = TRUE_;
	}
    }

/*     read the grid points where potential will be computed */

    if (dotarget) {
	ipot = freeunit_();
	o__1.oerr = 0;
	o__1.ounit = ipot;
	o__1.ofnmlen = 120;
	o__1.ofnm = potfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = ipot;
	f_rew(&al__1);
	i__2 = potfit_1.nconf;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    readpot_(&ipot, &i__);
	}
	cl__1.cerr = 0;
	cl__1.cunit = ipot;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     output the number of potential grid points to be used */

    i__2 = potfit_1.nconf;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (i__ == 1) {
	    io___67.ciunit = iounit_1.iout;
	    s_wsfe(&io___67);
	    e_wsfe();
	}
	if (potfit_1.npgrid[i__ - 1] > 100000) {
	    io___68.ciunit = iounit_1.iout;
	    s_wsfe(&io___68);
	    e_wsfe();
	    fatal_();
	} else if (potfit_1.nconf == 1) {
	    io___69.ciunit = iounit_1.iout;
	    s_wsfe(&io___69);
	    do_fio(&c__1, (char *)&potfit_1.npgrid[0], (ftnlen)sizeof(integer)
		    );
	    e_wsfe();
	} else {
	    io___70.ciunit = iounit_1.iout;
	    s_wsfe(&io___70);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&potfit_1.npgrid[i__ - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	}
    }

/*     output grid points at which to compute QM potential */

    if (dogrid) {
	igrd = freeunit_();
	s_copy(gridfile, files_1.filename, (ftnlen)120, (ftnlen)120);
	suffix_(gridfile, "grid", (ftnlen)120, (ftnlen)4);
	version_(gridfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = igrd;
	o__1.ofnmlen = 120;
	o__1.ofnm = gridfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	i__2 = potfit_1.nconf;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = potfit_1.npgrid[j - 1];
	    for (i__ = 1; i__ <= i__3; ++i__) {
		xi = pgrid_ref(1, i__, j);
		yi = pgrid_ref(2, i__, j);
		zi = pgrid_ref(3, i__, j);
		io___74.ciunit = igrd;
		s_wsfe(&io___74);
		do_fio(&c__1, (char *)&xi, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&yi, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&zi, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
	cl__1.cerr = 0;
	cl__1.cunit = igrd;
	cl__1.csta = 0;
	f_clos(&cl__1);
	io___75.ciunit = iounit_1.iout;
	s_wsfe(&io___75);
	do_fio(&c__1, gridfile, trimtext_(gridfile, (ftnlen)120));
	e_wsfe();
	io___76.ciunit = iounit_1.iout;
	s_wsfe(&io___76);
	e_wsfe();
    }

/*     get termination criterion for fitting as RMS gradient */

    if (dofit) {
	grdmin = -1.;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__2 = s_rsli(&io___78);
	    if (i__2 != 0) {
		goto L350;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&grdmin, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L350;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L350;
	    }
	}
L350:
	if (grdmin <= 0.) {
	    io___79.ciunit = iounit_1.iout;
	    s_wsfe(&io___79);
	    e_wsfe();
	    io___80.ciunit = iounit_1.input;
	    s_rsfe(&io___80);
	    do_fio(&c__1, (char *)&grdmin, (ftnlen)sizeof(doublereal));
	    e_rsfe();
	}
	if (grdmin <= 0.) {
	    grdmin = .001;
	}
    }

/*     setup the potential computation for each input structure */

    if (! dogrid) {
	i__2 = nmodel;
	for (k = 1; k <= i__2; ++k) {
	    if (k == 2) {
		ixyz = freeunit_();
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
		i__3 = potfit_1.nconf;
		for (j = 1; j <= i__3; ++j) {
		    readxyz_(&ixyz);
		    makeref_(&j);
		}
		cl__1.cerr = 0;
		cl__1.cunit = ixyz;
		cl__1.csta = 0;
		f_clos(&cl__1);
		cutoffs_();
		field_();
		katom_();
		if (potent_1.use_charge__ || potent_1.use_chgdpl__) {
		    kcharge_();
		}
		if (potent_1.use_dipole__ || potent_1.use_chgdpl__) {
		    kdipole_();
		}
		if (potent_1.use_mpole__ || potent_1.use_polar__) {
		    kmpole_();
		}
		if (potent_1.use_polar__) {
		    kpolar_();
		}
	    }

/*     compute the value of the potential at each grid point */

	    i__3 = potfit_1.nconf;
	    for (j = 1; j <= i__3; ++j) {
		getref_(&j);
		if (potent_1.use_mpole__) {
		    chkpole_();
		    rotpole_();
		}
		if (potent_1.use_polar__) {
		    induce_();
		}
		i__1 = potfit_1.npgrid[j - 1];
		for (i__ = 1; i__ <= i__1; ++i__) {
		    xi = pgrid_ref(1, i__, j);
		    yi = pgrid_ref(2, i__, j);
		    zi = pgrid_ref(3, i__, j);
		    potpoint_(&xi, &yi, &zi, &pot);
		    epot_ref(k, i__, j) = pot;
		}
	    }
	}

/*     output overall statistics for electrostatic potential */

	potstat_(&dofull, &domodel, &dopair, &dotarget);
    }

/*     set initial electrostatic parameters and do optimization */

    if (dofit) {
	prmvar_(&nvar, xx);
	minima_1.hguess = 1e-4;
	s_copy(output_1.coordtype, "NONE", (ftnlen)9, (ftnlen)4);
	ocvm_(&nvar, xx, &minimum, &grdmin, (D_fp)potfit1_, (U_fp)optsave_);

/*     move optimized values back to parameters and output */

	varprm_(&nvar, xx, &c__0, &c_b170);
	prmvar_(&nvar, xx);

/*     evaluate the potential value at each grid point */

	i__2 = potfit_1.nconf;
	for (j = 1; j <= i__2; ++j) {
	    getref_(&j);
	    if (potent_1.use_mpole__) {
		chkpole_();
		rotpole_();
	    }
	    if (potent_1.use_polar__) {
		induce_();
	    }
	    i__3 = potfit_1.npgrid[j - 1];
	    for (i__ = 1; i__ <= i__3; ++i__) {
		xi = pgrid_ref(1, i__, j);
		yi = pgrid_ref(2, i__, j);
		zi = pgrid_ref(3, i__, j);
		potpoint_(&xi, &yi, &zi, &pot);
		epot_ref(1, i__, j) = pot;
	    }
	}

/*     output overall statistics and optimized parameters */

	potstat_(&dofull, &domodel, &dopair, &dotarget);
	prtfit_();
    }
L380:
    return 0;
} /* MAIN__ */

#undef keyline_ref
#undef pgrid_ref
#undef epot_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine readpot  --  get and assign potential grid  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "readpot" gets a set of grid points and target electrostatic */
/*     potential values from an external disk file */


/* Subroutine */ int readpot_(integer *ipot, integer *iconf)
{
    /* Format strings */
    static char fmt_10[] = "(a120)";
    static char fmt_30[] = "(a120)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void),
	     s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal r2, xi, yi, zi, big, rad[25000], dist, small;
    static char record[120];
    static integer atmnum, npoint;

    /* Fortran I/O blocks */
    static cilist io___85 = { 1, 0, 1, fmt_10, 0 };
    static icilist io___87 = { 1, record, 1, 0, 120, 1 };
    static cilist io___89 = { 1, 0, 1, fmt_30, 0 };
    static icilist io___90 = { 1, record, 1, 0, 120, 1 };



#define epot_ref(a_1,a_2,a_3) potfit_1.epot[((a_3)*100000 + (a_2))*2 + \
a_1 - 200003]
#define pgrid_ref(a_1,a_2,a_3) potfit_1.pgrid[((a_3)*100000 + (a_2))*3 + \
a_1 - 300004]
#define ipgrid_ref(a_1,a_2) potfit_1.ipgrid[(a_2)*100000 + a_1 - 100001]



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




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  potfit.i  --  values for electrostatic potential fitting  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxpgrd   maximum dimension of electrostatic potential grid */

/*     pgrid     Cartesian coordinates of potential grid points */
/*     epot      values of electrostatic potential at grid points */
/*     xdpl0     target x-component of the molecular dipole moment */
/*     ydpl0     target y-component of the molecular dipole moment */
/*     zdpl0     target z-component of the molecular dipole moment */
/*     xxqdp0    target xx-component of the molecular quadrupole moment */
/*     xyqdp0    target xy-component of the molecular quadrupole moment */
/*     xzqdp0    target xz-component of the molecular quadrupole moment */
/*     yyqdp0    target yy-component of the molecular quadrupole moment */
/*     yzqdp0    target yz-component of the molecular quadrupole moment */
/*     zzqdp0    target zz-component of the molecular quadrupole moment */
/*     nconf     total number of conformations to be analyzed */
/*     npgrid    total number of electrostatic potential grid points */
/*     ipgrid    atom associated with each potential grid point */
/*     ngatm     total number of atoms with active potential grid points */
/*     nfatm     total number of atoms in electrostatic potential fit */
/*     gatm      flag to use potential grid points around each atom */
/*     fatm      flag to use each atom in electrostatic potential fit */
/*     use_dpl   flag to include molecular dipole in potential fit */
/*     use_qdp   flag to include molecular quadrupole in potential fit */
/*     fit_mpl   flag for atomic monopoles to vary in potential fit */
/*     fit_dpl   flag for atomic dipoles to vary in potential fit */
/*     fit_qdp   flag for atomic quadrupoles to vary in potential fit */




/*     read the grid points and target potential from a file */

    npoint = 0;
    io___85.ciunit = *ipot;
    i__1 = s_rsfe(&io___85);
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = do_fio(&c__1, record, (ftnlen)120);
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = s_rsli(&io___87);
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&npoint, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = e_rsli();
    if (i__1 != 0) {
	goto L20;
    }
L20:
    i__1 = npoint;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pgrid_ref(1, i__, *iconf) = 0.;
	pgrid_ref(2, i__, *iconf) = 0.;
	pgrid_ref(3, i__, *iconf) = 0.;
	epot_ref(2, i__, *iconf) = 0.;
	io___89.ciunit = *ipot;
	i__2 = s_rsfe(&io___89);
	if (i__2 != 0) {
	    goto L40;
	}
	i__2 = do_fio(&c__1, record, (ftnlen)120);
	if (i__2 != 0) {
	    goto L40;
	}
	i__2 = e_rsfe();
	if (i__2 != 0) {
	    goto L40;
	}
	i__2 = s_rsli(&io___90);
	if (i__2 != 0) {
	    goto L40;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	if (i__2 != 0) {
	    goto L40;
	}
	for (j = 1; j <= 3; ++j) {
	    i__2 = do_lio(&c__5, &c__1, (char *)&pgrid_ref(j, i__, *iconf), (
		    ftnlen)sizeof(doublereal));
	    if (i__2 != 0) {
		goto L40;
	    }
	}
	i__2 = do_lio(&c__5, &c__1, (char *)&epot_ref(2, i__, *iconf), (
		ftnlen)sizeof(doublereal));
	if (i__2 != 0) {
	    goto L40;
	}
	i__2 = e_rsli();
	if (i__2 != 0) {
	    goto L40;
	}
L40:
	;
    }

/*     set base atomic radii from traditional Bondi values */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rad[i__ - 1] = 1.7;
	atmnum = atmtyp_1.atomic[i__ - 1];
	if (atmnum == 0) {
	    rad[i__ - 1] = 0.;
	}
	if (atmnum == 1) {
	    rad[i__ - 1] = 1.2;
	}
	if (atmnum == 2) {
	    rad[i__ - 1] = 1.4;
	}
	if (atmnum == 6) {
	    rad[i__ - 1] = 1.7;
	}
	if (atmnum == 7) {
	    rad[i__ - 1] = 1.55;
	}
	if (atmnum == 8) {
	    rad[i__ - 1] = 1.52;
	}
	if (atmnum == 9) {
	    rad[i__ - 1] = 1.47;
	}
	if (atmnum == 10) {
	    rad[i__ - 1] = 1.54;
	}
	if (atmnum == 14) {
	    rad[i__ - 1] = 2.1;
	}
	if (atmnum == 15) {
	    rad[i__ - 1] = 1.8;
	}
	if (atmnum == 16) {
	    rad[i__ - 1] = 1.8;
	}
	if (atmnum == 17) {
	    rad[i__ - 1] = 1.75;
	}
	if (atmnum == 18) {
	    rad[i__ - 1] = 1.88;
	}
	if (atmnum == 35) {
	    rad[i__ - 1] = 1.85;
	}
	if (atmnum == 36) {
	    rad[i__ - 1] = 2.02;
	}
	if (atmnum == 53) {
	    rad[i__ - 1] = 1.98;
	}
	if (atmnum == 54) {
	    rad[i__ - 1] = 2.16;
	}
    }

/*     assign each grid point to atom on molecular surface */

    big = 1e3;
    i__1 = npoint;
    for (i__ = 1; i__ <= i__1; ++i__) {
	small = big;
	xi = pgrid_ref(1, i__, *iconf);
	yi = pgrid_ref(2, i__, *iconf);
	zi = pgrid_ref(3, i__, *iconf);
	i__2 = atoms_1.n;
	for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    d__1 = xi - atoms_1.x[k - 1];
/* Computing 2nd power */
	    d__2 = yi - atoms_1.y[k - 1];
/* Computing 2nd power */
	    d__3 = zi - atoms_1.z__[k - 1];
	    r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	    dist = sqrt(r2) - rad[k - 1];
	    if (dist < small) {
		small = dist;
		ipgrid_ref(i__, *iconf) = k;
	    }
	}
    }

/*     use potential grid points only for active grid atoms */

    k = npoint;
    npoint = 0;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (potfit_1.gatm[ipgrid_ref(i__, *iconf) - 1]) {
	    ++npoint;
	    ipgrid_ref(npoint, *iconf) = ipgrid_ref(i__, *iconf);
	    pgrid_ref(1, npoint, *iconf) = pgrid_ref(1, i__, *iconf);
	    pgrid_ref(2, npoint, *iconf) = pgrid_ref(2, i__, *iconf);
	    pgrid_ref(3, npoint, *iconf) = pgrid_ref(3, i__, *iconf);
	    epot_ref(2, npoint, *iconf) = epot_ref(2, i__, *iconf);
	}
    }
    potfit_1.npgrid[*iconf - 1] = npoint;
    return 0;
} /* readpot_ */

#undef ipgrid_ref
#undef pgrid_ref
#undef epot_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine potgrid  --  generate shells of grid points  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "potgrid" generates electrostatic potential grid points in */
/*     radially distributed shells outside the molecular surface */


/* Subroutine */ int potgrid_(integer *iconf)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 POTGRID  --  Too many Surface Grid Point"
	    "s;\002,\002 Increase MAXDOT\002)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), i_dnnt(doublereal *);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal r2, xi, yi, zi, xj, yj, zj, xr, yr, zr, rad[25000], dot[
	    150000]	/* was [3][50000] */, rad2[25000];
    static integer ndot, next;
    extern /* Subroutine */ int fatal_(void);
    static doublereal round;
    static char record[120];
    static integer nshell, atmnum;
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    static integer npoint;
    extern /* Subroutine */ int sphere_(integer *, doublereal *);
    static doublereal spacing, roffset, density;
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___113 = { 1, string, 1, 0, 120, 1 };
    static icilist io___114 = { 1, string, 1, 0, 120, 1 };
    static icilist io___115 = { 1, string, 1, 0, 120, 1 };
    static cilist io___124 = { 0, 0, 0, fmt_20, 0 };



#define dot_ref(a_1,a_2) dot[(a_2)*3 + a_1 - 4]
#define pgrid_ref(a_1,a_2,a_3) potfit_1.pgrid[((a_3)*100000 + (a_2))*3 + \
a_1 - 300004]
#define ipgrid_ref(a_1,a_2) potfit_1.ipgrid[(a_2)*100000 + a_1 - 100001]
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
/*     ##  keys.i  --  contents of current keyword parameter file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nkey      number of nonblank lines in the keyword file */
/*     keyline   contents of each individual keyword file line */




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




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  potfit.i  --  values for electrostatic potential fitting  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxpgrd   maximum dimension of electrostatic potential grid */

/*     pgrid     Cartesian coordinates of potential grid points */
/*     epot      values of electrostatic potential at grid points */
/*     xdpl0     target x-component of the molecular dipole moment */
/*     ydpl0     target y-component of the molecular dipole moment */
/*     zdpl0     target z-component of the molecular dipole moment */
/*     xxqdp0    target xx-component of the molecular quadrupole moment */
/*     xyqdp0    target xy-component of the molecular quadrupole moment */
/*     xzqdp0    target xz-component of the molecular quadrupole moment */
/*     yyqdp0    target yy-component of the molecular quadrupole moment */
/*     yzqdp0    target yz-component of the molecular quadrupole moment */
/*     zzqdp0    target zz-component of the molecular quadrupole moment */
/*     nconf     total number of conformations to be analyzed */
/*     npgrid    total number of electrostatic potential grid points */
/*     ipgrid    atom associated with each potential grid point */
/*     ngatm     total number of atoms with active potential grid points */
/*     nfatm     total number of atoms in electrostatic potential fit */
/*     gatm      flag to use potential grid points around each atom */
/*     fatm      flag to use each atom in electrostatic potential fit */
/*     use_dpl   flag to include molecular dipole in potential fit */
/*     use_qdp   flag to include molecular quadrupole in potential fit */
/*     fit_mpl   flag for atomic monopoles to vary in potential fit */
/*     fit_dpl   flag for atomic dipoles to vary in potential fit */
/*     fit_qdp   flag for atomic quadrupoles to vary in potential fit */




/*     set default values for grid point generation parameters */

    npoint = 0;
    nshell = 4;
    spacing = .35;
/* Computing 2nd power */
    d__1 = spacing;
    density = 12.566370614359172 / (d__1 * d__1);
    roffset = 1.;
    round = 1e-6;

/*     check for keywords containing any altered parameters */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "POTENTIAL-SHELLS ", (ftnlen)17, (ftnlen)17) == 0) 
		{
	    i__2 = s_rsli(&io___113);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&nshell, (ftnlen)sizeof(
		    integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "POTENTIAL-SPACING ", (ftnlen)18, (ftnlen)
		18) == 0) {
	    i__2 = s_rsli(&io___114);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&spacing, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
/* Computing 2nd power */
	    d__1 = spacing;
	    density = 12.566370614359172 / (d__1 * d__1);
	} else if (s_cmp(keyword, "POTENTIAL-OFFSET ", (ftnlen)17, (ftnlen)17)
		 == 0) {
	    i__2 = s_rsli(&io___115);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&roffset, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	}
L10:
	;
    }

/*     set base atomic radii from traditional Bondi values */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rad[i__ - 1] = 1.7;
	atmnum = atmtyp_1.atomic[i__ - 1];
	if (atmnum == 0) {
	    rad[i__ - 1] = 0.;
	}
	if (atmnum == 1) {
	    rad[i__ - 1] = 1.2;
	}
	if (atmnum == 2) {
	    rad[i__ - 1] = 1.4;
	}
	if (atmnum == 6) {
	    rad[i__ - 1] = 1.7;
	}
	if (atmnum == 7) {
	    rad[i__ - 1] = 1.55;
	}
	if (atmnum == 8) {
	    rad[i__ - 1] = 1.52;
	}
	if (atmnum == 9) {
	    rad[i__ - 1] = 1.47;
	}
	if (atmnum == 10) {
	    rad[i__ - 1] = 1.54;
	}
	if (atmnum == 14) {
	    rad[i__ - 1] = 2.1;
	}
	if (atmnum == 15) {
	    rad[i__ - 1] = 1.8;
	}
	if (atmnum == 16) {
	    rad[i__ - 1] = 1.8;
	}
	if (atmnum == 17) {
	    rad[i__ - 1] = 1.75;
	}
	if (atmnum == 18) {
	    rad[i__ - 1] = 1.88;
	}
	if (atmnum == 35) {
	    rad[i__ - 1] = 1.85;
	}
	if (atmnum == 36) {
	    rad[i__ - 1] = 2.02;
	}
	if (atmnum == 53) {
	    rad[i__ - 1] = 1.98;
	}
	if (atmnum == 54) {
	    rad[i__ - 1] = 2.16;
	}
	rad[i__ - 1] += roffset;
/* Computing 2nd power */
	d__1 = rad[i__ - 1];
	rad2[i__ - 1] = d__1 * d__1;
    }

/*     find points on each of the molecular surface shells */

    i__1 = nshell;
    for (m = 1; m <= i__1; ++m) {
	if (m != 1) {
	    i__2 = atoms_1.n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		rad[i__ - 1] += spacing;
/* Computing 2nd power */
		d__1 = rad[i__ - 1];
		rad2[i__ - 1] = d__1 * d__1;
	    }
	}
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    ndot = (integer) (density * rad2[i__ - 1]);
	    if (ndot > 50000) {
		io___124.ciunit = iounit_1.iout;
		s_wsfe(&io___124);
		e_wsfe();
		fatal_();
	    }
	    sphere_(&ndot, dot);
	    i__3 = ndot;
	    for (j = 1; j <= i__3; ++j) {
		xj = xi + rad[i__ - 1] * dot_ref(1, j);
		yj = yi + rad[i__ - 1] * dot_ref(2, j);
		zj = zi + rad[i__ - 1] * dot_ref(3, j);
		d__1 = xj / round;
		xj = (doublereal) i_dnnt(&d__1) * round;
		d__1 = yj / round;
		yj = (doublereal) i_dnnt(&d__1) * round;
		d__1 = zj / round;
		zj = (doublereal) i_dnnt(&d__1) * round;
		i__4 = i__ - 1;
		for (k = 1; k <= i__4; ++k) {
		    xr = xj - atoms_1.x[k - 1];
		    yr = yj - atoms_1.y[k - 1];
		    zr = zj - atoms_1.z__[k - 1];
		    r2 = xr * xr + yr * yr + zr * zr;
		    if (r2 < rad2[k - 1]) {
			goto L30;
		    }
		}
		i__4 = atoms_1.n;
		for (k = i__ + 1; k <= i__4; ++k) {
		    xr = xj - atoms_1.x[k - 1];
		    yr = yj - atoms_1.y[k - 1];
		    zr = zj - atoms_1.z__[k - 1];
		    r2 = xr * xr + yr * yr + zr * zr;
		    if (r2 < rad2[k - 1]) {
			goto L30;
		    }
		}
		++npoint;
		ipgrid_ref(npoint, *iconf) = i__;
		pgrid_ref(1, npoint, *iconf) = xj;
		pgrid_ref(2, npoint, *iconf) = yj;
		pgrid_ref(3, npoint, *iconf) = zj;
L30:
		;
	    }
	}
    }

/*     use potential grid points only for active grid atoms */

    k = npoint;
    npoint = 0;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (potfit_1.gatm[ipgrid_ref(i__, *iconf) - 1]) {
	    ++npoint;
	    ipgrid_ref(npoint, *iconf) = ipgrid_ref(i__, *iconf);
	    pgrid_ref(1, npoint, *iconf) = pgrid_ref(1, i__, *iconf);
	    pgrid_ref(2, npoint, *iconf) = pgrid_ref(2, i__, *iconf);
	    pgrid_ref(3, npoint, *iconf) = pgrid_ref(3, i__, *iconf);
	}
    }
    potfit_1.npgrid[*iconf - 1] = npoint;
    return 0;
} /* potgrid_ */

#undef keyline_ref
#undef ipgrid_ref
#undef pgrid_ref
#undef dot_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine potpoint  --  electrostatic potential at point  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "potpoint" calculates the electrostatic potential at a grid */
/*     point "i" as the total electrostatic interaction energy of */
/*     a molecule with a positive charge located at the grid point */


/* Subroutine */ int potpoint_(doublereal *xi, doublereal *yi, doublereal *zi,
	 doublereal *pot)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer k;
    static doublereal r__;
    static integer k1, k2;
    static doublereal r2, ec, ed, ci, ei, fi, ck, em, ep;
    static integer kk;
    static doublereal xk, yk, zk, xr, yr, zr, rk2, rr1, rr3, rr5, scd, dkx, 
	    dky, dkz, scq, scu, qkx, qky, qkz, ukx, uky, ukz, rkr3, dotk, 
	    qkxx, qkxy, qkxz, qkyy, qkyz, qkzz;


#define idpl_ref(a_1,a_2) dipole_1.idpl[(a_2)*2 + a_1 - 3]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]



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
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  chgpot.i  --  specifics of charge-charge functional form  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     electric   energy factor in kcal/mole for current force field */
/*     dielec     dielectric constant for electrostatic interactions */
/*     ebuffer    electrostatic buffering constant added to distance */
/*     c2scale    factor by which 1-2 charge interactions are scaled */
/*     c3scale    factor by which 1-3 charge interactions are scaled */
/*     c4scale    factor by which 1-4 charge interactions are scaled */
/*     c5scale    factor by which 1-5 charge interactions are scaled */
/*     neutnbr    logical flag governing use of neutral group neighbors */
/*     neutcut    logical flag governing use of neutral group cutoffs */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  dipole.i  --  atom & bond dipoles for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     bdpl      magnitude of each of the dipoles (Debyes) */
/*     sdpl      position of each dipole between defining atoms */
/*     ndipole   total number of dipoles in the system */
/*     idpl      numbers of atoms that define each dipole */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  polar.i  --  polarizabilities and induced dipole moments  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     polarity  dipole polarizability for each multipole site (Ang**3) */
/*     thole     Thole polarizability damping value for each site */
/*     pdamp     value of polarizability scale factor for each site */
/*     uind      induced dipole components at each multipole site */
/*     uinp      induced dipoles in field used for energy interactions */
/*     uinds     GK or PB induced dipoles at each multipole site */
/*     uinps     induced dipoles in field used for GK or PB energy */
/*     npolar    total number of polarizable sites in the system */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  units.i  --  physical constants and unit conversions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     literature reference: */

/*     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended */
/*     Values of the Fundamental Physical Constants: 2006", Reviews of */
/*     Modern Physics, 80, 633-730 (2008) */

/*     The "2006 CODATA Recommended Values" are also available from */
/*     the NIST Reference on Constants, Units, and Uncertainty site */
/*     at http://physics.nist.gov/cuu/index.html */

/*     Most values below are derived from the 2006 CODATA reference */

/*     The conversion from calorie to Joule is the definition of the */
/*     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992) */

/*     The "coulomb" energy conversion factor is found by dimensional */
/*     analysis of Coulomb's Law, ie, by dividing the square of the */
/*     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is */
/*     the permittivity of vacuum (the "electric constant"); note that */
/*     eps0 is typically given in F/m, equivalent to C**2/(J-m) */

/*     The approximate value used for the Debye, 3.33564 x 10-30 C-m, */
/*     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997) */

/*     The value of "prescon" is based on definition of 1 atmosphere */
/*     as 101325 Pa set by the 10th Conference Generale des Poids et */
/*     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3 */

/*     avogadro    Avogadro's number (N) in particles/mole */
/*     lightspd    speed of light in vacuum (c) in cm/ps */
/*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K */
/*     gasconst    ideal gas constant (R) in kcal/mole/K */
/*     emass       mass of an electron in atomic mass units */
/*     joule       conversion from calories to joules */
/*     convert     conversion from kcal to g*Ang**2/ps**2 */
/*     bohr        conversion from Bohrs to Angstroms */
/*     hartree     conversion from Hartree to kcal/mole */
/*     evolt       conversion from Hartree to electron-volts */
/*     efreq       conversion from Hartree to cm-1 */
/*     coulomb     conversion from electron**2/Ang to kcal/mole */
/*     debye       conversion from electron-Ang to Debyes */
/*     prescon     conversion from kcal/mole/Ang**3 to Atm */




/*     zero out charge, dipole and polarizable multipole energies */

    ec = 0.;
    ed = 0.;
    em = 0.;
    ep = 0.;

/*     set charge of probe site and electrostatic constants */

    f = chgpot_1.electric / chgpot_1.dielec;
    ci = 1.;
    fi = f * ci;

/*     calculate the charge contribution to the potential */

    if (potent_1.use_charge__) {
	i__1 = charge_1.nion;
	for (k = 1; k <= i__1; ++k) {
	    kk = charge_1.iion[k - 1];
	    xr = atoms_1.x[kk - 1] - *xi;
	    yr = atoms_1.y[kk - 1] - *yi;
	    zr = atoms_1.z__[kk - 1] - *zi;
	    r2 = xr * xr + yr * yr + zr * zr;
	    r__ = sqrt(r2);
	    e = fi * charge_1.pchg[k - 1] / r__;
	    ec += e;
	}
    }

/*     calculate the bond dipole contribution to the potential */

    if (potent_1.use_dipole__) {
	i__1 = dipole_1.ndipole;
	for (k = 1; k <= i__1; ++k) {
	    k1 = idpl_ref(1, k);
	    k2 = idpl_ref(2, k);
	    xk = atoms_1.x[k2 - 1] - atoms_1.x[k1 - 1];
	    yk = atoms_1.y[k2 - 1] - atoms_1.y[k1 - 1];
	    zk = atoms_1.z__[k2 - 1] - atoms_1.z__[k1 - 1];
	    xr = atoms_1.x[k1 - 1] + xk * dipole_1.sdpl[k - 1] - *xi;
	    yr = atoms_1.y[k1 - 1] + yk * dipole_1.sdpl[k - 1] - *yi;
	    zr = atoms_1.z__[k1 - 1] + zk * dipole_1.sdpl[k - 1] - *zi;
	    r2 = xr * xr + yr * yr + zr * zr;
	    rk2 = xk * xk + yk * yk + zk * zk;
	    rkr3 = sqrt(rk2 * r2) * r2;
	    dotk = xk * xr + yk * yr + zk * zr;
	    e = fi / 4.80321 * dipole_1.bdpl[k - 1] * dotk / rkr3;
	    ed += e;
	}
    }

/*     calculate the multipole contribution to the potential */

    if (potent_1.use_mpole__ || potent_1.use_polar__) {
	i__1 = mpole_1.npole;
	for (k = 1; k <= i__1; ++k) {
	    kk = mpole_1.ipole[k - 1];
	    xr = atoms_1.x[kk - 1] - *xi;
	    yr = atoms_1.y[kk - 1] - *yi;
	    zr = atoms_1.z__[kk - 1] - *zi;
	    r2 = xr * xr + yr * yr + zr * zr;
	    r__ = sqrt(r2);
	    ck = rpole_ref(1, k);
	    dkx = rpole_ref(2, k);
	    dky = rpole_ref(3, k);
	    dkz = rpole_ref(4, k);
	    qkxx = rpole_ref(5, k);
	    qkxy = rpole_ref(6, k);
	    qkxz = rpole_ref(7, k);
	    qkyy = rpole_ref(9, k);
	    qkyz = rpole_ref(10, k);
	    qkzz = rpole_ref(13, k);
	    ukx = uind_ref(1, k);
	    uky = uind_ref(2, k);
	    ukz = uind_ref(3, k);

/*     construct some intermediate quadrupole values */

	    qkx = qkxx * xr + qkxy * yr + qkxz * zr;
	    qky = qkxy * xr + qkyy * yr + qkyz * zr;
	    qkz = qkxz * xr + qkyz * yr + qkzz * zr;

/*     calculate scalar products for permanent and induced */

	    scd = dkx * xr + dky * yr + dkz * zr;
	    scq = qkx * xr + qky * yr + qkz * zr;
	    scu = ukx * xr + uky * yr + ukz * zr;

/*     compute the energy contributions for this interaction */

	    rr1 = 1. / r__;
	    rr3 = rr1 / r2;
	    rr5 = rr3 * 3. / r2;
	    e = ck * rr1 - scd * rr3 + scq * rr5;
	    ei = -scu * rr3;

/*     increment the overall multipole and polarization energies */

	    e = fi * e;
	    ei = fi * ei;
	    em += e;
	    ep += ei;
	}
    }

/*     potential is sum of all interactions with probe site */

    *pot = ec + ed + em + ep;
    return 0;
} /* potpoint_ */

#undef rpole_ref
#undef uind_ref
#undef idpl_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function potfit1  --  potential fit error and gradient  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "potfit1" is a service routine that computes the RMS error */
/*     and gradient for electrostatic parameters fit to a potential */


doublereal potfit1_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal ret_val, d__1;

    /* Builtin functions */
    integer i_dnnt(doublereal *);
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int potpoint_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal e;
    static integer i__, j, k, m;
    static doublereal e0, ec, er, et, xi, yi, zi, eps, pot;
    static integer nvar;
    static doublereal cscale;
    extern /* Subroutine */ int induce_(void);
    static doublereal tscale;
    extern /* Subroutine */ int getref_(integer *);
    static integer npoint;
    extern /* Subroutine */ int varprm_(integer *, doublereal *, integer *, 
	    doublereal *), chkpole_(void), momfull_(void), rotpole_(void);


#define epot_ref(a_1,a_2,a_3) potfit_1.epot[((a_3)*100000 + (a_2))*2 + \
a_1 - 200003]
#define pgrid_ref(a_1,a_2,a_3) potfit_1.pgrid[((a_3)*100000 + (a_2))*3 + \
a_1 - 300004]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  moment.i  --  components of electric multipole moments  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     netchg   net electric charge for the total system */
/*     netdpl   dipole moment magnitude for the total system */
/*     netqdp   diagonal quadrupole (Qxx, Qyy, Qzz) for system */
/*     xdpl     dipole vector x-component in the global frame */
/*     ydpl     dipole vector y-component in the global frame */
/*     zdpl     dipole vector z-component in the global frame */
/*     xxqdp    quadrupole tensor xx-component in global frame */
/*     xyqdp    quadrupole tensor xy-component in global frame */
/*     xzqdp    quadrupole tensor xz-component in global frame */
/*     yxqdp    quadrupole tensor yx-component in global frame */
/*     yyqdp    quadrupole tensor yy-component in global frame */
/*     yzqdp    quadrupole tensor yz-component in global frame */
/*     zxqdp    quadrupole tensor zx-component in global frame */
/*     zyqdp    quadrupole tensor zy-component in global frame */
/*     zzqdp    quadrupole tensor zz-component in global frame */




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




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  potfit.i  --  values for electrostatic potential fitting  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxpgrd   maximum dimension of electrostatic potential grid */

/*     pgrid     Cartesian coordinates of potential grid points */
/*     epot      values of electrostatic potential at grid points */
/*     xdpl0     target x-component of the molecular dipole moment */
/*     ydpl0     target y-component of the molecular dipole moment */
/*     zdpl0     target z-component of the molecular dipole moment */
/*     xxqdp0    target xx-component of the molecular quadrupole moment */
/*     xyqdp0    target xy-component of the molecular quadrupole moment */
/*     xzqdp0    target xz-component of the molecular quadrupole moment */
/*     yyqdp0    target yy-component of the molecular quadrupole moment */
/*     yzqdp0    target yz-component of the molecular quadrupole moment */
/*     zzqdp0    target zz-component of the molecular quadrupole moment */
/*     nconf     total number of conformations to be analyzed */
/*     npgrid    total number of electrostatic potential grid points */
/*     ipgrid    atom associated with each potential grid point */
/*     ngatm     total number of atoms with active potential grid points */
/*     nfatm     total number of atoms in electrostatic potential fit */
/*     gatm      flag to use potential grid points around each atom */
/*     fatm      flag to use each atom in electrostatic potential fit */
/*     use_dpl   flag to include molecular dipole in potential fit */
/*     use_qdp   flag to include molecular quadrupole in potential fit */
/*     fit_mpl   flag for atomic monopoles to vary in potential fit */
/*     fit_dpl   flag for atomic dipoles to vary in potential fit */
/*     fit_qdp   flag for atomic quadrupoles to vary in potential fit */




/*     initialize scaling factors for error and gradient */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    npoint = 0;
    i__1 = potfit_1.nconf;
    for (j = 1; j <= i__1; ++j) {
	npoint += potfit_1.npgrid[j - 1];
    }
    cscale = 1e8 / (doublereal) potfit_1.nconf;
    tscale = 1e4 / (doublereal) potfit_1.nconf;
    eps = 1e-6;

/*     copy optimization values to electrostatic parameters */

    varprm_(&nvar, &xx[1], &c__0, &c_b170);

/*     find total error by cycling over all conformations */

    er = 0.;
    ec = 0.;
    et = 0.;
    i__1 = potfit_1.nconf;
    for (j = 1; j <= i__1; ++j) {
	getref_(&j);
	if (potent_1.use_mpole__) {
	    chkpole_();
	    rotpole_();
	}
	if (potent_1.use_polar__) {
	    induce_();
	}

/*     get the RMS potential error summed over grid points */

	i__2 = potfit_1.npgrid[j - 1];
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xi = pgrid_ref(1, i__, j);
	    yi = pgrid_ref(2, i__, j);
	    zi = pgrid_ref(3, i__, j);
	    potpoint_(&xi, &yi, &zi, &pot);
	    epot_ref(1, i__, j) = pot;
/* Computing 2nd power */
	    d__1 = epot_ref(1, i__, j) - epot_ref(2, i__, j);
	    er += d__1 * d__1;
	}

/*     get deviation from integral net molecular charge */

	momfull_();
/* Computing 2nd power */
	d__1 = moment_1.netchg - (doublereal) i_dnnt(&moment_1.netchg);
	ec += cscale * (d__1 * d__1);

/*     get deviation from dipole and quadrupole targets */

	if (potfit_1.use_dpl__) {
/* Computing 2nd power */
	    d__1 = moment_1.xdpl - potfit_1.xdpl0[j - 1];
	    et += tscale * (d__1 * d__1);
/* Computing 2nd power */
	    d__1 = moment_1.ydpl - potfit_1.ydpl0[j - 1];
	    et += tscale * (d__1 * d__1);
/* Computing 2nd power */
	    d__1 = moment_1.zdpl - potfit_1.zdpl0[j - 1];
	    et += tscale * (d__1 * d__1);
	}
	if (potfit_1.use_qdp__) {
/* Computing 2nd power */
	    d__1 = moment_1.xxqdp - potfit_1.xxqdp0[j - 1];
	    et += tscale * (d__1 * d__1);
/* Computing 2nd power */
	    d__1 = moment_1.xyqdp - potfit_1.xyqdp0[j - 1];
	    et += tscale * (d__1 * d__1);
/* Computing 2nd power */
	    d__1 = moment_1.xzqdp - potfit_1.xzqdp0[j - 1];
	    et += tscale * (d__1 * d__1);
/* Computing 2nd power */
	    d__1 = moment_1.yyqdp - potfit_1.yyqdp0[j - 1];
	    et += tscale * (d__1 * d__1);
/* Computing 2nd power */
	    d__1 = moment_1.yzqdp - potfit_1.yzqdp0[j - 1];
	    et += tscale * (d__1 * d__1);
/* Computing 2nd power */
	    d__1 = moment_1.zzqdp - potfit_1.zzqdp0[j - 1];
	    et += tscale * (d__1 * d__1);
	}
    }
    er = sqrt(er / (doublereal) npoint);
    ret_val = er + ec + et;

/*     compute numerical gradient for electrostatic parameters */

    m = nvar;
    i__1 = m;
    for (k = 1; k <= i__1; ++k) {
	d__1 = eps * -.5;
	varprm_(&nvar, &xx[1], &k, &d__1);
	er = 0.;
	ec = 0.;
	et = 0.;
	i__2 = potfit_1.nconf;
	for (j = 1; j <= i__2; ++j) {
	    getref_(&j);
	    if (potent_1.use_mpole__) {
		chkpole_();
		rotpole_();
	    }
	    if (potent_1.use_polar__) {
		induce_();
	    }
	    i__3 = potfit_1.npgrid[j - 1];
	    for (i__ = 1; i__ <= i__3; ++i__) {
		xi = pgrid_ref(1, i__, j);
		yi = pgrid_ref(2, i__, j);
		zi = pgrid_ref(3, i__, j);
		potpoint_(&xi, &yi, &zi, &pot);
		epot_ref(1, i__, j) = pot;
/* Computing 2nd power */
		d__1 = epot_ref(1, i__, j) - epot_ref(2, i__, j);
		er += d__1 * d__1;
	    }
	    momfull_();
/* Computing 2nd power */
	    d__1 = moment_1.netchg - (doublereal) i_dnnt(&moment_1.netchg);
	    ec += cscale * (d__1 * d__1);
	    if (potfit_1.use_dpl__) {
/* Computing 2nd power */
		d__1 = moment_1.xdpl - potfit_1.xdpl0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.ydpl - potfit_1.ydpl0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.zdpl - potfit_1.zdpl0[j - 1];
		et += tscale * (d__1 * d__1);
	    }
	    if (potfit_1.use_qdp__) {
/* Computing 2nd power */
		d__1 = moment_1.xxqdp - potfit_1.xxqdp0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.xyqdp - potfit_1.xyqdp0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.xzqdp - potfit_1.xzqdp0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.yyqdp - potfit_1.yyqdp0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.yzqdp - potfit_1.yzqdp0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.zzqdp - potfit_1.zzqdp0[j - 1];
		et += tscale * (d__1 * d__1);
	    }
	}
	er = sqrt(er / (doublereal) npoint);
	e0 = er + ec + et;
	d__1 = eps * .5;
	varprm_(&nvar, &xx[1], &k, &d__1);
	er = 0.;
	ec = 0.;
	et = 0.;
	i__2 = potfit_1.nconf;
	for (j = 1; j <= i__2; ++j) {
	    getref_(&j);
	    if (potent_1.use_mpole__) {
		chkpole_();
		rotpole_();
	    }
	    if (potent_1.use_polar__) {
		induce_();
	    }
	    i__3 = potfit_1.npgrid[j - 1];
	    for (i__ = 1; i__ <= i__3; ++i__) {
		xi = pgrid_ref(1, i__, j);
		yi = pgrid_ref(2, i__, j);
		zi = pgrid_ref(3, i__, j);
		potpoint_(&xi, &yi, &zi, &pot);
		epot_ref(1, i__, j) = pot;
/* Computing 2nd power */
		d__1 = epot_ref(1, i__, j) - epot_ref(2, i__, j);
		er += d__1 * d__1;
	    }
	    momfull_();
/* Computing 2nd power */
	    d__1 = moment_1.netchg - (doublereal) i_dnnt(&moment_1.netchg);
	    ec += cscale * (d__1 * d__1);
	    if (potfit_1.use_dpl__) {
/* Computing 2nd power */
		d__1 = moment_1.xdpl - potfit_1.xdpl0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.ydpl - potfit_1.ydpl0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.zdpl - potfit_1.zdpl0[j - 1];
		et += tscale * (d__1 * d__1);
	    }
	    if (potfit_1.use_qdp__) {
/* Computing 2nd power */
		d__1 = moment_1.xxqdp - potfit_1.xxqdp0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.xyqdp - potfit_1.xyqdp0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.xzqdp - potfit_1.xzqdp0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.yyqdp - potfit_1.yyqdp0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.yzqdp - potfit_1.yzqdp0[j - 1];
		et += tscale * (d__1 * d__1);
/* Computing 2nd power */
		d__1 = moment_1.zzqdp - potfit_1.zzqdp0[j - 1];
		et += tscale * (d__1 * d__1);
	    }
	}
	er = sqrt(er / (doublereal) npoint);
	e = er + ec + et;
	g[k] = (e - e0) / eps;
    }
    return ret_val;
} /* potfit1_ */

#undef pgrid_ref
#undef epot_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine prmvar  --  electrostatics to optimization  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "prmvar" determines the optimization values from the */
/*     corresponding electrostatic potential energy parameters */


/* Subroutine */ int prmvar_(integer *nvar, doublereal *xx)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Atomic Parameters Included in Potential "
	    "Fitting :\002,//,3x,\002Atom\002,10x,\002Atom Name\002,9x,\002At"
	    "om Type\002,9x,\002Parameters\002,/)";
    static char fmt_20[] = "(i6,15x,a3,10x,i6,13x,a)";
    static char fmt_30[] = "(i6,15x,a3,10x,i6,13x,a)";
    static char fmt_40[] = "(/,\002 Potential Fitting of Electrostatic Param"
	    "eters :\002,//,1x,\002Parameter\002,6x,\002Atom Type\002,9x,\002"
	    "Category\002,12x,\002Value\002,8x,\002Fixed\002,/)";
    static char fmt_50[] = "(i6,7x,i8,13x,a8,5x,f12.5)";
    static char fmt_60[] = "(4x,\002--\002,7x,i8,13x,a8,5x,f12.5,10x,\002"
	    "X\002)";
    static char fmt_70[] = "(i6,7x,i8,13x,a8,5x,f12.5)";
    static char fmt_80[] = "(4x,\002--\002,7x,i8,13x,a8,5x,f12.5,10x,\002"
	    "X\002)";
    static char fmt_90[] = "(i6,7x,i8,13x,a8,5x,f12.5)";
    static char fmt_100[] = "(4x,\002--\002,7x,i8,13x,a8,5x,f12.5,10x,\002"
	    "X\002)";
    static char fmt_110[] = "(i6,7x,i8,13x,a8,5x,f12.5)";
    static char fmt_120[] = "(4x,\002--\002,7x,i8,13x,a8,5x,f12.5,10x,\002"
	    "X\002)";
    static char fmt_130[] = "(i6,7x,i8,13x,a8,5x,f12.5)";
    static char fmt_140[] = "(4x,\002--\002,7x,i8,13x,a8,5x,f12.5,10x,\002"
	    "X\002)";
    static char fmt_150[] = "(i6,7x,i8,13x,a8,5x,f12.5)";
    static char fmt_160[] = "(4x,\002--\002,7x,i8,13x,a8,5x,f12.5,10x,\002"
	    "X\002)";
    static char fmt_170[] = "(i6,7x,i8,13x,a8,5x,f12.5)";
    static char fmt_180[] = "(4x,\002--\002,7x,i8,13x,a8,5x,f12.5,10x,\002"
	    "X\002)";
    static char fmt_190[] = "(i6,7x,i8,13x,a8,5x,f12.5)";
    static char fmt_200[] = "(4x,\002--\002,7x,i8,13x,a8,5x,f12.5,10x,\002"
	    "X\002)";
    static char fmt_210[] = "(i6,7x,i8,13x,a8,5x,f12.5)";
    static char fmt_220[] = "(4x,\002--\002,7x,i8,13x,a8,5x,f12.5,10x,\002"
	    "X\002)";
    static char fmt_230[] = "(i6,7x,i8,13x,a8,5x,f12.5)";
    static char fmt_240[] = "(4x,\002--\002,7x,i8,13x,a8,5x,f12.5,10x,\002"
	    "X\002)";
    static char fmt_250[] = "(4x,\002--\002,7x,i8,13x,a8,5x,f12.5)";
    static char fmt_260[] = "(4x,\002--\002,7x,i8,13x,a8,5x,f12.5,10x,\002"
	    "X\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, k, ii, kk, it, kt;
    static logical done;
    static doublereal dterm, qterm;
    static char prmtyp[17];

    /* Fortran I/O blocks */
    static cilist io___201 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___206 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___207 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___208 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___213 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___214 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___215 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___216 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___217 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___218 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___219 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___220 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___221 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___222 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___223 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___224 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___225 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___226 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___227 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___228 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___229 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___230 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___231 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___232 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___233 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___234 = { 0, 0, 0, fmt_260, 0 };



#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define pole_ref(a_1,a_2) mpole_1.pole[(a_2)*13 + a_1 - 14]



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
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  potfit.i  --  values for electrostatic potential fitting  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxpgrd   maximum dimension of electrostatic potential grid */

/*     pgrid     Cartesian coordinates of potential grid points */
/*     epot      values of electrostatic potential at grid points */
/*     xdpl0     target x-component of the molecular dipole moment */
/*     ydpl0     target y-component of the molecular dipole moment */
/*     zdpl0     target z-component of the molecular dipole moment */
/*     xxqdp0    target xx-component of the molecular quadrupole moment */
/*     xyqdp0    target xy-component of the molecular quadrupole moment */
/*     xzqdp0    target xz-component of the molecular quadrupole moment */
/*     yyqdp0    target yy-component of the molecular quadrupole moment */
/*     yzqdp0    target yz-component of the molecular quadrupole moment */
/*     zzqdp0    target zz-component of the molecular quadrupole moment */
/*     nconf     total number of conformations to be analyzed */
/*     npgrid    total number of electrostatic potential grid points */
/*     ipgrid    atom associated with each potential grid point */
/*     ngatm     total number of atoms with active potential grid points */
/*     nfatm     total number of atoms in electrostatic potential fit */
/*     gatm      flag to use potential grid points around each atom */
/*     fatm      flag to use each atom in electrostatic potential fit */
/*     use_dpl   flag to include molecular dipole in potential fit */
/*     use_qdp   flag to include molecular quadrupole in potential fit */
/*     fit_mpl   flag for atomic monopoles to vary in potential fit */
/*     fit_dpl   flag for atomic dipoles to vary in potential fit */
/*     fit_qdp   flag for atomic quadrupoles to vary in potential fit */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  units.i  --  physical constants and unit conversions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     literature reference: */

/*     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended */
/*     Values of the Fundamental Physical Constants: 2006", Reviews of */
/*     Modern Physics, 80, 633-730 (2008) */

/*     The "2006 CODATA Recommended Values" are also available from */
/*     the NIST Reference on Constants, Units, and Uncertainty site */
/*     at http://physics.nist.gov/cuu/index.html */

/*     Most values below are derived from the 2006 CODATA reference */

/*     The conversion from calorie to Joule is the definition of the */
/*     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992) */

/*     The "coulomb" energy conversion factor is found by dimensional */
/*     analysis of Coulomb's Law, ie, by dividing the square of the */
/*     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is */
/*     the permittivity of vacuum (the "electric constant"); note that */
/*     eps0 is typically given in F/m, equivalent to C**2/(J-m) */

/*     The approximate value used for the Debye, 3.33564 x 10-30 C-m, */
/*     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997) */

/*     The value of "prescon" is based on definition of 1 atmosphere */
/*     as 101325 Pa set by the 10th Conference Generale des Poids et */
/*     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3 */

/*     avogadro    Avogadro's number (N) in particles/mole */
/*     lightspd    speed of light in vacuum (c) in cm/ps */
/*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K */
/*     gasconst    ideal gas constant (R) in kcal/mole/K */
/*     emass       mass of an electron in atomic mass units */
/*     joule       conversion from calories to joules */
/*     convert     conversion from kcal to g*Ang**2/ps**2 */
/*     bohr        conversion from Bohrs to Angstroms */
/*     hartree     conversion from Hartree to kcal/mole */
/*     evolt       conversion from Hartree to electron-volts */
/*     efreq       conversion from Hartree to cm-1 */
/*     coulomb     conversion from electron**2/Ang to kcal/mole */
/*     debye       conversion from electron-Ang to Debyes */
/*     prescon     conversion from kcal/mole/Ang**3 to Atm */




/*     zero out the total number of optimization parameters */

    /* Parameter adjustments */
    --xx;

    /* Function Body */
    *nvar = 0;
    dterm = 1.8897261328856432;
    qterm = 10.713194571932783;

/*     list active atoms when not all are used in optimization */

    if (potfit_1.nfatm != atoms_1.n) {
	io___201.ciunit = iounit_1.iout;
	s_wsfe(&io___201);
	e_wsfe();
	i__1 = charge_1.nion;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = charge_1.iion[i__ - 1];
	    if (potfit_1.fatm[ii - 1]) {
		it = atoms_1.type__[ii - 1];
		s_copy(prmtyp, "Partial Charge", (ftnlen)17, (ftnlen)14);
		io___206.ciunit = iounit_1.iout;
		s_wsfe(&io___206);
		do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, prmtyp, (ftnlen)17);
		e_wsfe();
	    }
	}
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    if (potfit_1.fatm[ii - 1]) {
		it = atoms_1.type__[ii - 1];
		s_copy(prmtyp, "Atomic Multipoles", (ftnlen)17, (ftnlen)17);
		io___207.ciunit = iounit_1.iout;
		s_wsfe(&io___207);
		do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, prmtyp, (ftnlen)17);
		e_wsfe();
	    }
	}
    }

/*     print header information for electrostatic parameters */

    io___208.ciunit = iounit_1.iout;
    s_wsfe(&io___208);
    e_wsfe();

/*     get optimization parameters from partial charge values */

    i__1 = charge_1.nion;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = TRUE_;
	ii = charge_1.iion[i__ - 1];
	if (potfit_1.fatm[ii - 1]) {
	    done = FALSE_;
	}
	if (! done) {
	    it = atoms_1.type__[ii - 1];
	    i__2 = i__ - 1;
	    for (k = 1; k <= i__2; ++k) {
		kk = charge_1.iion[k - 1];
		kt = atoms_1.type__[kk - 1];
		if (kt == it && potfit_1.fatm[kk - 1]) {
		    done = TRUE_;
		}
	    }
	}
	if (! done) {
	    if (charge_1.pchg[i__ - 1] != 0.) {
		++(*nvar);
		xx[*nvar] = charge_1.pchg[i__ - 1];
		io___213.ciunit = iounit_1.iout;
		s_wsfe(&io___213);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Charge  ", (ftnlen)8);
		do_fio(&c__1, (char *)&charge_1.pchg[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___214.ciunit = iounit_1.iout;
		s_wsfe(&io___214);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Charge  ", (ftnlen)8);
		do_fio(&c__1, (char *)&charge_1.pchg[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     get optimization parameters from atomic multipole values */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = TRUE_;
	ii = mpole_1.ipole[i__ - 1];
	if (potfit_1.fatm[ii - 1]) {
	    done = FALSE_;
	}
	if (! done) {
	    it = atoms_1.type__[ii - 1];
	    i__2 = i__ - 1;
	    for (k = 1; k <= i__2; ++k) {
		kk = mpole_1.ipole[k - 1];
		kt = atoms_1.type__[kk - 1];
		if (kt == it && potfit_1.fatm[kk - 1]) {
		    done = TRUE_;
		}
	    }
	}
	if (! done) {
	    if (potfit_1.fit_mpl__ && pole_ref(1, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = pole_ref(1, i__);
		io___215.ciunit = iounit_1.iout;
		s_wsfe(&io___215);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Monopole", (ftnlen)8);
		do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___216.ciunit = iounit_1.iout;
		s_wsfe(&io___216);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Monopole", (ftnlen)8);
		do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	    if (potfit_1.fit_dpl__ && pole_ref(2, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = pole_ref(2, i__);
		io___217.ciunit = iounit_1.iout;
		s_wsfe(&io___217);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "X-Dipole", (ftnlen)8);
		d__1 = dterm * pole_ref(2, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___218.ciunit = iounit_1.iout;
		s_wsfe(&io___218);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "X-Dipole", (ftnlen)8);
		d__1 = dterm * pole_ref(2, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (potfit_1.fit_dpl__ && pole_ref(3, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = pole_ref(3, i__);
		io___219.ciunit = iounit_1.iout;
		s_wsfe(&io___219);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Y-Dipole", (ftnlen)8);
		d__1 = dterm * pole_ref(3, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___220.ciunit = iounit_1.iout;
		s_wsfe(&io___220);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Y-Dipole", (ftnlen)8);
		d__1 = dterm * pole_ref(3, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (potfit_1.fit_dpl__ && pole_ref(4, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = pole_ref(4, i__);
		io___221.ciunit = iounit_1.iout;
		s_wsfe(&io___221);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Z-Dipole", (ftnlen)8);
		d__1 = dterm * pole_ref(4, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___222.ciunit = iounit_1.iout;
		s_wsfe(&io___222);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Z-Dipole", (ftnlen)8);
		d__1 = dterm * pole_ref(4, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (potfit_1.fit_qdp__ && pole_ref(5, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = pole_ref(5, i__);
		io___223.ciunit = iounit_1.iout;
		s_wsfe(&io___223);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "XX-Quad ", (ftnlen)8);
		d__1 = qterm * pole_ref(5, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___224.ciunit = iounit_1.iout;
		s_wsfe(&io___224);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "XX-Quad ", (ftnlen)8);
		d__1 = qterm * pole_ref(5, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (potfit_1.fit_qdp__ && pole_ref(6, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = pole_ref(6, i__);
		io___225.ciunit = iounit_1.iout;
		s_wsfe(&io___225);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "XY-Quad ", (ftnlen)8);
		d__1 = qterm * pole_ref(6, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___226.ciunit = iounit_1.iout;
		s_wsfe(&io___226);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "XY-Quad ", (ftnlen)8);
		d__1 = qterm * pole_ref(6, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (potfit_1.fit_qdp__ && pole_ref(7, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = pole_ref(7, i__);
		io___227.ciunit = iounit_1.iout;
		s_wsfe(&io___227);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "XZ-Quad ", (ftnlen)8);
		d__1 = qterm * pole_ref(7, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___228.ciunit = iounit_1.iout;
		s_wsfe(&io___228);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "XZ-Quad ", (ftnlen)8);
		d__1 = qterm * pole_ref(7, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (potfit_1.fit_qdp__ && pole_ref(9, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = pole_ref(9, i__);
		io___229.ciunit = iounit_1.iout;
		s_wsfe(&io___229);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "YY-Quad ", (ftnlen)8);
		d__1 = qterm * pole_ref(9, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___230.ciunit = iounit_1.iout;
		s_wsfe(&io___230);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "YY-Quad ", (ftnlen)8);
		d__1 = qterm * pole_ref(9, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (potfit_1.fit_qdp__ && pole_ref(10, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = pole_ref(10, i__);
		io___231.ciunit = iounit_1.iout;
		s_wsfe(&io___231);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "YZ-Quad ", (ftnlen)8);
		d__1 = qterm * pole_ref(10, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___232.ciunit = iounit_1.iout;
		s_wsfe(&io___232);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "YZ-Quad ", (ftnlen)8);
		d__1 = qterm * pole_ref(10, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (potfit_1.fit_qdp__ && pole_ref(13, i__) != 0.) {
		io___233.ciunit = iounit_1.iout;
		s_wsfe(&io___233);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "ZZ-Quad ", (ftnlen)8);
		d__1 = qterm * pole_ref(13, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___234.ciunit = iounit_1.iout;
		s_wsfe(&io___234);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, "ZZ-Quad ", (ftnlen)8);
		d__1 = qterm * pole_ref(13, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
    }
    return 0;
} /* prmvar_ */

#undef pole_ref
#undef name___ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine varprm  --  optimization to electrostatics  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "varprm" copies the current optimization values into the */
/*     corresponding electrostatic potential energy parameters */


/* Subroutine */ int varprm_(integer *nvar, doublereal *xx, integer *ivar, 
	doublereal *eps)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, ii, kk, it, kt;
    static logical done;
    extern /* Subroutine */ int chkpole_(void), rotpole_(void);


#define pole_ref(a_1,a_2) mpole_1.pole[(a_2)*13 + a_1 - 14]



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
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




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




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  potfit.i  --  values for electrostatic potential fitting  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxpgrd   maximum dimension of electrostatic potential grid */

/*     pgrid     Cartesian coordinates of potential grid points */
/*     epot      values of electrostatic potential at grid points */
/*     xdpl0     target x-component of the molecular dipole moment */
/*     ydpl0     target y-component of the molecular dipole moment */
/*     zdpl0     target z-component of the molecular dipole moment */
/*     xxqdp0    target xx-component of the molecular quadrupole moment */
/*     xyqdp0    target xy-component of the molecular quadrupole moment */
/*     xzqdp0    target xz-component of the molecular quadrupole moment */
/*     yyqdp0    target yy-component of the molecular quadrupole moment */
/*     yzqdp0    target yz-component of the molecular quadrupole moment */
/*     zzqdp0    target zz-component of the molecular quadrupole moment */
/*     nconf     total number of conformations to be analyzed */
/*     npgrid    total number of electrostatic potential grid points */
/*     ipgrid    atom associated with each potential grid point */
/*     ngatm     total number of atoms with active potential grid points */
/*     nfatm     total number of atoms in electrostatic potential fit */
/*     gatm      flag to use potential grid points around each atom */
/*     fatm      flag to use each atom in electrostatic potential fit */
/*     use_dpl   flag to include molecular dipole in potential fit */
/*     use_qdp   flag to include molecular quadrupole in potential fit */
/*     fit_mpl   flag for atomic monopoles to vary in potential fit */
/*     fit_dpl   flag for atomic dipoles to vary in potential fit */
/*     fit_qdp   flag for atomic quadrupoles to vary in potential fit */




/*     zero out the total number of optimization parameters */

    /* Parameter adjustments */
    --xx;

    /* Function Body */
    *nvar = 0;

/*     translate optimization values back to partial charges */

    i__1 = charge_1.nion;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = TRUE_;
	ii = charge_1.iion[i__ - 1];
	if (potfit_1.fatm[ii - 1]) {
	    done = FALSE_;
	}
	if (! done) {
	    it = atoms_1.type__[ii - 1];
	    i__2 = i__ - 1;
	    for (k = 1; k <= i__2; ++k) {
		kk = charge_1.iion[k - 1];
		kt = atoms_1.type__[kk - 1];
		if (kt == it && potfit_1.fatm[kk - 1]) {
		    done = TRUE_;
		}
	    }
	}
	if (! done) {
	    if (charge_1.pchg[i__ - 1] != 0.) {
		++(*nvar);
		charge_1.pchg[i__ - 1] = xx[*nvar];
		if (*ivar == *nvar) {
		    charge_1.pchg[i__ - 1] += *eps;
		}
		i__2 = charge_1.nion;
		for (k = i__ + 1; k <= i__2; ++k) {
		    kk = charge_1.iion[k - 1];
		    kt = atoms_1.type__[kk - 1];
		    if (kt == it && potfit_1.fatm[kk - 1]) {
			charge_1.pchg[k - 1] = charge_1.pchg[i__ - 1];
		    }
		}
	    }
	}
    }

/*     translate optimization values back to atomic multipoles */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = TRUE_;
	ii = mpole_1.ipole[i__ - 1];
	if (potfit_1.fatm[ii - 1]) {
	    done = FALSE_;
	}
	if (! done) {
	    it = atoms_1.type__[ii - 1];
	    i__2 = i__ - 1;
	    for (k = 1; k <= i__2; ++k) {
		kk = mpole_1.ipole[k - 1];
		kt = atoms_1.type__[kk - 1];
		if (kt == it && potfit_1.fatm[kk - 1]) {
		    done = TRUE_;
		}
	    }
	}
	if (! done) {
	    if (potfit_1.fit_mpl__ && pole_ref(1, i__) != 0.) {
		++(*nvar);
		pole_ref(1, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    pole_ref(1, i__) = pole_ref(1, i__) + *eps;
		}
	    }
	    if (potfit_1.fit_dpl__ && pole_ref(2, i__) != 0.) {
		++(*nvar);
		pole_ref(2, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    pole_ref(2, i__) = pole_ref(2, i__) + *eps;
		}
	    }
	    if (potfit_1.fit_dpl__ && pole_ref(3, i__) != 0.) {
		++(*nvar);
		pole_ref(3, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    pole_ref(3, i__) = pole_ref(3, i__) + *eps;
		}
	    }
	    if (potfit_1.fit_dpl__ && pole_ref(4, i__) != 0.) {
		++(*nvar);
		pole_ref(4, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    pole_ref(4, i__) = pole_ref(4, i__) + *eps;
		}
	    }
	    if (potfit_1.fit_qdp__ && pole_ref(5, i__) != 0.) {
		++(*nvar);
		pole_ref(5, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    pole_ref(5, i__) = pole_ref(5, i__) + *eps;
		}
	    }
	    if (potfit_1.fit_qdp__ && pole_ref(6, i__) != 0.) {
		++(*nvar);
		pole_ref(6, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    pole_ref(6, i__) = pole_ref(6, i__) + *eps;
		}
		pole_ref(8, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    pole_ref(8, i__) = pole_ref(8, i__) + *eps;
		}
	    }
	    if (potfit_1.fit_qdp__ && pole_ref(7, i__) != 0.) {
		++(*nvar);
		pole_ref(7, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    pole_ref(7, i__) = pole_ref(7, i__) + *eps;
		}
		pole_ref(11, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    pole_ref(11, i__) = pole_ref(11, i__) + *eps;
		}
	    }
	    if (potfit_1.fit_qdp__ && pole_ref(9, i__) != 0.) {
		++(*nvar);
		pole_ref(9, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    pole_ref(9, i__) = pole_ref(9, i__) + *eps;
		}
	    }
	    if (potfit_1.fit_qdp__ && pole_ref(10, i__) != 0.) {
		++(*nvar);
		pole_ref(10, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    pole_ref(10, i__) = pole_ref(10, i__) + *eps;
		}
		pole_ref(12, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    pole_ref(12, i__) = pole_ref(12, i__) + *eps;
		}
	    }
	    if (potfit_1.fit_qdp__ && pole_ref(13, i__) != 0.) {
		pole_ref(13, i__) = -pole_ref(5, i__) - pole_ref(9, i__);
	    }
	    i__2 = mpole_1.npole;
	    for (k = i__ + 1; k <= i__2; ++k) {
		kk = mpole_1.ipole[k - 1];
		kt = atoms_1.type__[kk - 1];
		if (kt == it && potfit_1.fatm[kk - 1]) {
		    for (j = 1; j <= 13; ++j) {
			pole_ref(j, k) = pole_ref(j, i__);
		    }
		}
	    }
	}
    }

/*     check chiral multipoles and rotate into global frame */

    if (potent_1.use_mpole__) {
	chkpole_();
	rotpole_();
    }
    return 0;
} /* varprm_ */

#undef pole_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine potstat  --  electrostatic potential statistics  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "potstat" computes and prints statistics for the electrostatic */
/*     potential over a set of grid points */


/* Subroutine */ int potstat_(logical *dofull, logical *domodel, logical *
	dopair, logical *dotarget)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Electrostatic Potential at Each Gri"
	    "d\002,\002 Point :\002,/,8x,\002(Kcal/mole per unit charge)\002)";
    static char fmt_20[] = "(/,\002 Electrostatic Potential at Grid Point"
	    "s\002,\002 for Conformer\002,i4,\002 :\002,/,12x,\002(Kcal/mole "
	    "per unit charge)\002)";
    static char fmt_30[] = "(/,3x,\002Point\002,15x,\002XYZ-Coordinates\002,"
	    "15x,\002Potential\002,5x,\002Target\002,/)";
    static char fmt_40[] = "(/,3x,\002Point\002,15x,\002XYZ-Coordinates\002,"
	    "13x,\002Potential 1\002,3x,\002Potential 2\002,/)";
    static char fmt_50[] = "(/,3x,\002Point\002,15x,\002XYZ-Coordinates\002,"
	    "14x,\002Potential\002,/)";
    static char fmt_60[] = "(i8,2x,a)";
    static char fmt_70[] = "(i8,3x,3f12.6,2x,2f12.4)";
    static char fmt_80[] = "(i8,3x,3f12.6,2x,f12.4)";
    static char fmt_90[] = "(i8,3x,3f12.6,2x,f12.4)";
    static char fmt_100[] = "(/,\002 Electrostatic Potential Written to File"
	    " :  \002,a)";
    static char fmt_110[] = "(/,\002 Average Electrostatic Potential over At"
	    "oms :\002,/,6x,\002(Kcal/mole per unit charge)\002)";
    static char fmt_120[] = "(/,4x,\002Atom\002,9x,\002Points\002,8x,\002Pot"
	    "ential\002,8x,\002Target\002,8x,\002RMS Diff\002,/)";
    static char fmt_130[] = "(/,4x,\002Atom\002,9x,\002Points\002,7x,\002Pot"
	    "ential 1\002,4x,\002Potential 2\002,6x,\002RMS Diff\002,/)";
    static char fmt_140[] = "(/,4x,\002Atom\002,9x,\002Points\002,8x,\002Pot"
	    "ential\002,/)";
    static char fmt_150[] = "(i8,6x,i8,5x,f12.4,3x,f12.4,3x,f12.4)";
    static char fmt_160[] = "(i8,6x,i8,5x,f12.4)";
    static char fmt_170[] = "(/,\002 Electrostatic Potential over all Grid P"
	    "oints :\002,//,\002 Average Magnitude for Potential 1 :\002,6x,f"
	    "12.4)";
    static char fmt_180[] = "(/,\002 Electrostatic Potential over all Grid P"
	    "oints :\002,//,\002 Average Magnitude for Potential :\002,8x,f12"
	    ".4)";
    static char fmt_190[] = "(\002 Average Magnitude for Target :\002,11x,f1"
	    "2.4,/,\002 Average Signed Potential Difference :\002,4x,f12.4,/"
	    ",\002 Average Unsigned Potential Difference :\002,2x,f12.4,/,"
	    "\002 Root Mean Square Potential Difference :\002,2x,f12.4)";
    static char fmt_200[] = "(\002 Average Magnitude for Potential 2 :\002,6"
	    "x,f12.4,/,\002 Average Signed Potential Difference :\002,4x,f12."
	    "4,/,\002 Average Unsigned Potential Difference :\002,2x,f12.4,/"
	    ",\002 Root Mean Square Potential Difference :\002,2x,f12.4)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3;
    doublereal d__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), s_wsfe(cilist *), e_wsfe(void), do_fio(integer *,
	     char *, ftnlen), f_clos(cllist *);
    double sqrt(doublereal);

    /* Local variables */
    extern integer freeunit_(void), trimtext_(char *, ftnlen);
    static integer i__, j, k;
    static doublereal xi, yi, zi;
    static integer natm[25000];
    static doublereal tave, uave, rmsa[25000], rmsd;
    static integer ipot;
    static doublereal pave1, pave2, patm1[25000], patm2[25000];
    static integer npoint;
    static char potfile[120];
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___246 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___247 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___248 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___249 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___250 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___251 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___256 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___257 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___258 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___259 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___260 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___261 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___262 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___263 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___269 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___270 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___277 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___278 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___279 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___280 = { 0, 0, 0, fmt_200, 0 };



#define epot_ref(a_1,a_2,a_3) potfit_1.epot[((a_3)*100000 + (a_2))*2 + \
a_1 - 200003]
#define pgrid_ref(a_1,a_2,a_3) potfit_1.pgrid[((a_3)*100000 + (a_2))*3 + \
a_1 - 300004]
#define ipgrid_ref(a_1,a_2) potfit_1.ipgrid[(a_2)*100000 + a_1 - 100001]



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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  potfit.i  --  values for electrostatic potential fitting  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxpgrd   maximum dimension of electrostatic potential grid */

/*     pgrid     Cartesian coordinates of potential grid points */
/*     epot      values of electrostatic potential at grid points */
/*     xdpl0     target x-component of the molecular dipole moment */
/*     ydpl0     target y-component of the molecular dipole moment */
/*     zdpl0     target z-component of the molecular dipole moment */
/*     xxqdp0    target xx-component of the molecular quadrupole moment */
/*     xyqdp0    target xy-component of the molecular quadrupole moment */
/*     xzqdp0    target xz-component of the molecular quadrupole moment */
/*     yyqdp0    target yy-component of the molecular quadrupole moment */
/*     yzqdp0    target yz-component of the molecular quadrupole moment */
/*     zzqdp0    target zz-component of the molecular quadrupole moment */
/*     nconf     total number of conformations to be analyzed */
/*     npgrid    total number of electrostatic potential grid points */
/*     ipgrid    atom associated with each potential grid point */
/*     ngatm     total number of atoms with active potential grid points */
/*     nfatm     total number of atoms in electrostatic potential fit */
/*     gatm      flag to use potential grid points around each atom */
/*     fatm      flag to use each atom in electrostatic potential fit */
/*     use_dpl   flag to include molecular dipole in potential fit */
/*     use_qdp   flag to include molecular quadrupole in potential fit */
/*     fit_mpl   flag for atomic monopoles to vary in potential fit */
/*     fit_dpl   flag for atomic dipoles to vary in potential fit */
/*     fit_qdp   flag for atomic quadrupoles to vary in potential fit */




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




/*     output potential values for each model at each point */

    if (*dofull) {
	if (*domodel) {
	    ipot = freeunit_();
/* Writing concatenation */
	    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
	    i__1[1] = 4, a__1[1] = ".pot";
	    s_cat(potfile, a__1, i__1, &c__2, (ftnlen)120);
	    version_(potfile, "new", (ftnlen)120, (ftnlen)3);
	    o__1.oerr = 0;
	    o__1.ounit = ipot;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = potfile;
	    o__1.orl = 0;
	    o__1.osta = "new";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	}
	i__2 = potfit_1.nconf;
	for (j = 1; j <= i__2; ++j) {
	    if (potfit_1.nconf == 1) {
		io___246.ciunit = iounit_1.iout;
		s_wsfe(&io___246);
		e_wsfe();
	    } else {
		io___247.ciunit = iounit_1.iout;
		s_wsfe(&io___247);
		do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	    if (*dotarget) {
		io___248.ciunit = iounit_1.iout;
		s_wsfe(&io___248);
		e_wsfe();
	    } else if (*dopair) {
		io___249.ciunit = iounit_1.iout;
		s_wsfe(&io___249);
		e_wsfe();
	    } else if (*domodel) {
		io___250.ciunit = iounit_1.iout;
		s_wsfe(&io___250);
		e_wsfe();
		io___251.ciunit = ipot;
		s_wsfe(&io___251);
		do_fio(&c__1, (char *)&potfit_1.npgrid[j - 1], (ftnlen)sizeof(
			integer));
		do_fio(&c__1, titles_1.title, titles_1.ltitle);
		e_wsfe();
	    }
	    i__3 = potfit_1.npgrid[j - 1];
	    for (i__ = 1; i__ <= i__3; ++i__) {
		xi = pgrid_ref(1, i__, j);
		yi = pgrid_ref(2, i__, j);
		zi = pgrid_ref(3, i__, j);
		if (*dotarget || *dopair) {
		    io___256.ciunit = iounit_1.iout;
		    s_wsfe(&io___256);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&xi, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&yi, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&zi, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&epot_ref(1, i__, j), (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&epot_ref(2, i__, j), (ftnlen)
			    sizeof(doublereal));
		    e_wsfe();
		} else if (*domodel) {
		    io___257.ciunit = iounit_1.iout;
		    s_wsfe(&io___257);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&xi, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&yi, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&zi, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&epot_ref(1, i__, j), (ftnlen)
			    sizeof(doublereal));
		    e_wsfe();
		    io___258.ciunit = ipot;
		    s_wsfe(&io___258);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&xi, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&yi, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&zi, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&epot_ref(1, i__, j), (ftnlen)
			    sizeof(doublereal));
		    e_wsfe();
		}
	    }
	}
	if (*domodel) {
	    cl__1.cerr = 0;
	    cl__1.cunit = ipot;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	    io___259.ciunit = iounit_1.iout;
	    s_wsfe(&io___259);
	    do_fio(&c__1, potfile, trimtext_(potfile, (ftnlen)120));
	    e_wsfe();
	}
    }

/*     find average electrostatic potential around each atom */

    io___260.ciunit = iounit_1.iout;
    s_wsfe(&io___260);
    e_wsfe();
    if (*dotarget) {
	io___261.ciunit = iounit_1.iout;
	s_wsfe(&io___261);
	e_wsfe();
    } else if (*dopair) {
	io___262.ciunit = iounit_1.iout;
	s_wsfe(&io___262);
	e_wsfe();
    } else if (*domodel) {
	io___263.ciunit = iounit_1.iout;
	s_wsfe(&io___263);
	e_wsfe();
    }
    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	natm[i__ - 1] = 0;
	patm1[i__ - 1] = 0.;
	patm2[i__ - 1] = 0.;
	rmsa[i__ - 1] = 0.;
    }
    i__2 = potfit_1.nconf;
    for (j = 1; j <= i__2; ++j) {
	i__3 = potfit_1.npgrid[j - 1];
	for (i__ = 1; i__ <= i__3; ++i__) {
	    k = ipgrid_ref(i__, j);
	    ++natm[k - 1];
	    patm1[k - 1] += epot_ref(1, i__, j);
	    patm2[k - 1] += epot_ref(2, i__, j);
/* Computing 2nd power */
	    d__1 = epot_ref(1, i__, j) - epot_ref(2, i__, j);
	    rmsa[k - 1] += d__1 * d__1;
	}
    }
    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (natm[i__ - 1] != 0) {
	    patm1[i__ - 1] /= (doublereal) natm[i__ - 1];
	    patm2[i__ - 1] /= (doublereal) natm[i__ - 1];
	    rmsa[i__ - 1] = sqrt(rmsa[i__ - 1] / (doublereal) natm[i__ - 1]);
	}
	if (potfit_1.gatm[i__ - 1]) {
	    if (*dotarget || *dopair) {
		io___269.ciunit = iounit_1.iout;
		s_wsfe(&io___269);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&natm[i__ - 1], (ftnlen)sizeof(integer))
			;
		do_fio(&c__1, (char *)&patm1[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&patm2[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&rmsa[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else if (*domodel) {
		io___270.ciunit = iounit_1.iout;
		s_wsfe(&io___270);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&natm[i__ - 1], (ftnlen)sizeof(integer))
			;
		do_fio(&c__1, (char *)&patm1[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     overall averages for the sets of electrostatic potentials */

    npoint = 0;
    pave1 = 0.;
    pave2 = 0.;
    tave = 0.;
    uave = 0.;
    rmsd = 0.;
    i__2 = potfit_1.nconf;
    for (j = 1; j <= i__2; ++j) {
	npoint += potfit_1.npgrid[j - 1];
	i__3 = potfit_1.npgrid[j - 1];
	for (i__ = 1; i__ <= i__3; ++i__) {
	    pave1 += (d__1 = epot_ref(1, i__, j), abs(d__1));
	    pave2 += (d__1 = epot_ref(2, i__, j), abs(d__1));
	    tave = tave + epot_ref(1, i__, j) - epot_ref(2, i__, j);
	    uave += (d__1 = epot_ref(1, i__, j) - epot_ref(2, i__, j), abs(
		    d__1));
/* Computing 2nd power */
	    d__1 = epot_ref(1, i__, j) - epot_ref(2, i__, j);
	    rmsd += d__1 * d__1;
	}
    }
    pave1 /= (doublereal) npoint;
    pave2 /= (doublereal) npoint;
    tave /= (doublereal) npoint;
    uave /= (doublereal) npoint;
    rmsd = sqrt(rmsd / (doublereal) npoint);
    if (*dopair) {
	io___277.ciunit = iounit_1.iout;
	s_wsfe(&io___277);
	do_fio(&c__1, (char *)&pave1, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else {
	io___278.ciunit = iounit_1.iout;
	s_wsfe(&io___278);
	do_fio(&c__1, (char *)&pave1, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (*dotarget) {
	io___279.ciunit = iounit_1.iout;
	s_wsfe(&io___279);
	do_fio(&c__1, (char *)&pave2, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&tave, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&uave, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rmsd, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else if (*dopair) {
	io___280.ciunit = iounit_1.iout;
	s_wsfe(&io___280);
	do_fio(&c__1, (char *)&pave2, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&tave, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&uave, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rmsd, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* potstat_ */

#undef ipgrid_ref
#undef pgrid_ref
#undef epot_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine prtfit  --  create file with optimal parameters  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "prtfit" makes a key file containing results from fitting a */
/*     charge or multipole model to an electrostatic potential grid */


/* Subroutine */ int prtfit_(void)
{
    /* Format strings */
    static char fmt_10[] = "(a)";
    static char fmt_20[] = "(/,\002#\002,/,\002# Results of Electrostatic Po"
	    "tential Fitting\002,/,\002#\002)";
    static char fmt_30[] = "()";
    static char fmt_40[] = "(\002charge\002,4x,i5,10x,f11.4)";
    static char fmt_50[] = "()";
    static char fmt_60[] = "(\002multipole\002,1x,3i5,11x,f11.5)";
    static char fmt_70[] = "(\002multipole\002,1x,4i5,6x,f11.5)";
    static char fmt_80[] = "(\002multipole\002,1x,3i5,11x,f11.5)";
    static char fmt_90[] = "(\002multipole\002,1x,4i5,6x,f11.5)";
    static char fmt_100[] = "(\002multipole\002,1x,4i5,6x,f11.5)";
    static char fmt_110[] = "(\002multipole\002,1x,4i5,6x,f11.5)";
    static char fmt_120[] = "(36x,3f11.5)";
    static char fmt_130[] = "(36x,f11.5)";
    static char fmt_140[] = "(36x,2f11.5)";
    static char fmt_150[] = "(36x,3f11.5)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2], i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer i_dnnt(doublereal *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_cmp(char *, char *, ftnlen, ftnlen), f_clos(cllist *);

    /* Local variables */
    extern integer freeunit_(void), trimtext_(char *, ftnlen);
    static integer i__, j, k, ii, kk, it, kt, ix, iy, iz;
    static doublereal big, eps, sum;
    static integer ikey, size;
    static doublereal dterm, qterm;
    static logical header;
    static char record[120], keyfile[120];
    static logical doprint;
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___293 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___294 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___301 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___302 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___303 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___307 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___308 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___309 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___310 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___311 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___312 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___313 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___314 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___315 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___316 = { 0, 0, 0, fmt_150, 0 };



#define pole_ref(a_1,a_2) mpole_1.pole[(a_2)*13 + a_1 - 14]
#define polaxe_ref(a_0,a_1) &mpole_1.polaxe[(a_1)*8 + a_0 - 8]
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
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  keys.i  --  contents of current keyword parameter file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nkey      number of nonblank lines in the keyword file */
/*     keyline   contents of each individual keyword file line */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  potfit.i  --  values for electrostatic potential fitting  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxpgrd   maximum dimension of electrostatic potential grid */

/*     pgrid     Cartesian coordinates of potential grid points */
/*     epot      values of electrostatic potential at grid points */
/*     xdpl0     target x-component of the molecular dipole moment */
/*     ydpl0     target y-component of the molecular dipole moment */
/*     zdpl0     target z-component of the molecular dipole moment */
/*     xxqdp0    target xx-component of the molecular quadrupole moment */
/*     xyqdp0    target xy-component of the molecular quadrupole moment */
/*     xzqdp0    target xz-component of the molecular quadrupole moment */
/*     yyqdp0    target yy-component of the molecular quadrupole moment */
/*     yzqdp0    target yz-component of the molecular quadrupole moment */
/*     zzqdp0    target zz-component of the molecular quadrupole moment */
/*     nconf     total number of conformations to be analyzed */
/*     npgrid    total number of electrostatic potential grid points */
/*     ipgrid    atom associated with each potential grid point */
/*     ngatm     total number of atoms with active potential grid points */
/*     nfatm     total number of atoms in electrostatic potential fit */
/*     gatm      flag to use potential grid points around each atom */
/*     fatm      flag to use each atom in electrostatic potential fit */
/*     use_dpl   flag to include molecular dipole in potential fit */
/*     use_qdp   flag to include molecular quadrupole in potential fit */
/*     fit_mpl   flag for atomic monopoles to vary in potential fit */
/*     fit_dpl   flag for atomic dipoles to vary in potential fit */
/*     fit_qdp   flag for atomic quadrupoles to vary in potential fit */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  units.i  --  physical constants and unit conversions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     literature reference: */

/*     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended */
/*     Values of the Fundamental Physical Constants: 2006", Reviews of */
/*     Modern Physics, 80, 633-730 (2008) */

/*     The "2006 CODATA Recommended Values" are also available from */
/*     the NIST Reference on Constants, Units, and Uncertainty site */
/*     at http://physics.nist.gov/cuu/index.html */

/*     Most values below are derived from the 2006 CODATA reference */

/*     The conversion from calorie to Joule is the definition of the */
/*     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992) */

/*     The "coulomb" energy conversion factor is found by dimensional */
/*     analysis of Coulomb's Law, ie, by dividing the square of the */
/*     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is */
/*     the permittivity of vacuum (the "electric constant"); note that */
/*     eps0 is typically given in F/m, equivalent to C**2/(J-m) */

/*     The approximate value used for the Debye, 3.33564 x 10-30 C-m, */
/*     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997) */

/*     The value of "prescon" is based on definition of 1 atmosphere */
/*     as 101325 Pa set by the 10th Conference Generale des Poids et */
/*     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3 */

/*     avogadro    Avogadro's number (N) in particles/mole */
/*     lightspd    speed of light in vacuum (c) in cm/ps */
/*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K */
/*     gasconst    ideal gas constant (R) in kcal/mole/K */
/*     emass       mass of an electron in atomic mass units */
/*     joule       conversion from calories to joules */
/*     convert     conversion from kcal to g*Ang**2/ps**2 */
/*     bohr        conversion from Bohrs to Angstroms */
/*     hartree     conversion from Hartree to kcal/mole */
/*     evolt       conversion from Hartree to electron-volts */
/*     efreq       conversion from Hartree to cm-1 */
/*     coulomb     conversion from electron**2/Ang to kcal/mole */
/*     debye       conversion from electron-Ang to Debyes */
/*     prescon     conversion from kcal/mole/Ang**3 to Atm */




/*     convert dipole and quadrupole moments to atomic units */

    dterm = 1.8897261328856432;
    qterm = 10.713194571932783;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 2; j <= 4; ++j) {
	    pole_ref(j, i__) = dterm * pole_ref(j, i__);
	}
	for (j = 5; j <= 13; ++j) {
	    pole_ref(j, i__) = qterm * pole_ref(j, i__);
	}
    }

/*     regularize the multipole moments to desired precision */

    eps = 1e-5;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 13; ++j) {
	    d__1 = pole_ref(j, i__) / eps;
	    pole_ref(j, i__) = (doublereal) i_dnnt(&d__1) * eps;
	}
    }

/*     maintain traceless quadrupole at each multipole site */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = pole_ref(5, i__) + pole_ref(9, i__) + pole_ref(13, i__);
/* Computing MAX */
	d__4 = (d__1 = pole_ref(5, i__), abs(d__1)), d__5 = (d__2 = pole_ref(
		9, i__), abs(d__2)), d__4 = max(d__4,d__5), d__5 = (d__3 = 
		pole_ref(13, i__), abs(d__3));
	big = max(d__4,d__5);
	k = 0;
	if (big == (d__1 = pole_ref(5, i__), abs(d__1))) {
	    k = 5;
	}
	if (big == (d__1 = pole_ref(9, i__), abs(d__1))) {
	    k = 9;
	}
	if (big == (d__1 = pole_ref(13, i__), abs(d__1))) {
	    k = 13;
	}
	if (k != 0) {
	    pole_ref(k, i__) = pole_ref(k, i__) - sum;
	}
    }

/*     open a new keyfile to contain the optimized parameters */

    ikey = freeunit_();
/* Writing concatenation */
    i__2[0] = files_1.leng, a__1[0] = files_1.filename;
    i__2[1] = 4, a__1[1] = ".key";
    s_cat(keyfile, a__1, i__2, &c__2, (ftnlen)120);
    version_(keyfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = ikey;
    o__1.ofnmlen = 120;
    o__1.ofnm = keyfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);

/*     copy the contents of any previously existing keyfile */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	size = trimtext_(record, (ftnlen)120);
	io___293.ciunit = ikey;
	s_wsfe(&io___293);
	do_fio(&c__1, record, size);
	e_wsfe();
    }

/*     print a header for the fitted multipole parameters */

    io___294.ciunit = ikey;
    s_wsfe(&io___294);
    e_wsfe();

/*     output the optimized charge values to the keyfile */

    header = TRUE_;
    i__1 = charge_1.nion;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = charge_1.iion[i__ - 1];
	it = atoms_1.type__[ii - 1];
	doprint = FALSE_;
	if (potfit_1.fatm[ii - 1]) {
	    doprint = TRUE_;
	    i__3 = i__ - 1;
	    for (k = 1; k <= i__3; ++k) {
		kk = charge_1.iion[k - 1];
		kt = atoms_1.type__[kk - 1];
		if (potfit_1.fatm[kk - 1] && it == kt) {
		    doprint = FALSE_;
		}
	    }
	}
	if (doprint) {
	    if (header) {
		header = FALSE_;
		io___301.ciunit = ikey;
		s_wsfe(&io___301);
		e_wsfe();
	    }
	    io___302.ciunit = ikey;
	    s_wsfe(&io___302);
	    do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&charge_1.pchg[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }

/*     output the optimized multipole values to the keyfile */

    header = TRUE_;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	it = atoms_1.type__[ii - 1];
	doprint = FALSE_;
	if (potfit_1.fatm[ii - 1]) {
	    doprint = TRUE_;
	    i__3 = i__ - 1;
	    for (k = 1; k <= i__3; ++k) {
		kk = mpole_1.ipole[k - 1];
		kt = atoms_1.type__[kk - 1];
		if (potfit_1.fatm[kk - 1] && it == kt) {
		    doprint = FALSE_;
		}
	    }
	}
	if (doprint) {
	    if (header) {
		header = FALSE_;
		io___303.ciunit = ikey;
		s_wsfe(&io___303);
		e_wsfe();
	    }
	    iz = mpole_1.zaxis[i__ - 1];
	    ix = mpole_1.xaxis[i__ - 1];
	    iy = mpole_1.yaxis[i__ - 1];
	    if (iz != 0) {
		iz = atoms_1.type__[iz - 1];
	    }
	    if (ix != 0) {
		ix = atoms_1.type__[ix - 1];
	    }
	    if (iy != 0) {
		iy = atoms_1.type__[iy - 1];
	    }
	    if (s_cmp(polaxe_ref(0, i__), "Z-then-X", (ftnlen)8, (ftnlen)8) ==
		     0) {
		if (mpole_1.yaxis[i__ - 1] == 0) {
		    io___307.ciunit = ikey;
		    s_wsfe(&io___307);
		    do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&iz, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ix, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
			    doublereal));
		    e_wsfe();
		} else {
		    io___308.ciunit = ikey;
		    s_wsfe(&io___308);
		    do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&iz, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ix, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&iy, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
			    doublereal));
		    e_wsfe();
		}
	    } else if (s_cmp(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (
		    ftnlen)8) == 0) {
		if (mpole_1.yaxis[i__ - 1] == 0) {
		    io___309.ciunit = ikey;
		    s_wsfe(&io___309);
		    do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&iz, (ftnlen)sizeof(integer));
		    i__3 = -ix;
		    do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
			    doublereal));
		    e_wsfe();
		} else {
		    io___310.ciunit = ikey;
		    s_wsfe(&io___310);
		    do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&iz, (ftnlen)sizeof(integer));
		    i__3 = -ix;
		    do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&iy, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
			    doublereal));
		    e_wsfe();
		}
	    } else if (s_cmp(polaxe_ref(0, i__), "Z-Bisect", (ftnlen)8, (
		    ftnlen)8) == 0) {
		io___311.ciunit = ikey;
		s_wsfe(&io___311);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iz, (ftnlen)sizeof(integer));
		i__3 = -ix;
		do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		i__4 = -iy;
		do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else if (s_cmp(polaxe_ref(0, i__), "3-Fold", (ftnlen)8, (ftnlen)
		    6) == 0) {
		io___312.ciunit = ikey;
		s_wsfe(&io___312);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		i__3 = -iz;
		do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		i__4 = -ix;
		do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
		i__5 = -iy;
		do_fio(&c__1, (char *)&i__5, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	    io___313.ciunit = ikey;
	    s_wsfe(&io___313);
	    do_fio(&c__1, (char *)&pole_ref(2, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(3, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(4, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___314.ciunit = ikey;
	    s_wsfe(&io___314);
	    do_fio(&c__1, (char *)&pole_ref(5, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___315.ciunit = ikey;
	    s_wsfe(&io___315);
	    do_fio(&c__1, (char *)&pole_ref(8, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(9, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___316.ciunit = ikey;
	    s_wsfe(&io___316);
	    do_fio(&c__1, (char *)&pole_ref(11, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(12, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(13, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }
    cl__1.cerr = 0;
    cl__1.cunit = ikey;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* prtfit_ */

#undef keyline_ref
#undef polaxe_ref
#undef pole_ref


/* Main program alias */ int potential_ () { MAIN__ (); return 0; }
