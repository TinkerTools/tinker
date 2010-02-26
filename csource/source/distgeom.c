/* distgeom.f -- translated by f2c (version 20050501).
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
    doublereal ak[75000], anat[75000], afld[75000];
    integer nangle, iang[300000]	/* was [4][75000] */;
} angle_;

#define angle_1 angle_

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
    doublereal bnd[1000000]	/* was [1000][1000] */, vdwrad[25000], vdwmax,
	     compact, pathmax;
    logical use_invert__, use_anneal__;
} disgeo_;

#define disgeo_1 disgeo_

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
    doublereal xpfix[25000], ypfix[25000], zpfix[25000], pfix[50000]	/* 
	    was [2][25000] */, dfix[75000]	/* was [3][25000] */, afix[
	    75000]	/* was [3][25000] */, tfix[75000]	/* was [3][
	    25000] */, gfix[75000]	/* was [3][25000] */, chir[75000]	
	    /* was [3][25000] */, depth, width, rwall;
    integer npfix, ipfix[25000], kpfix[75000]	/* was [3][25000] */, ndfix, 
	    idfix[50000]	/* was [2][25000] */, nafix, iafix[75000]	
	    /* was [3][25000] */, ntfix, itfix[100000]	/* was [4][25000] */, 
	    ngfix, igfix[50000]	/* was [2][25000] */, nchir, ichir[100000]	
	    /* was [4][25000] */;
    logical use_basin__, use_wall__;
} kgeoms_;

#define kgeoms_1 kgeoms_

struct {
    doublereal xref[250000]	/* was [25000][10] */, yref[250000]	/* 
	    was [25000][10] */, zref[250000]	/* was [25000][10] */;
    integer nref[10], reftyp[250000]	/* was [25000][10] */, n12ref[250000]	
	    /* was [25000][10] */, i12ref[2000000]	/* was [8][25000][10] 
	    */, refleng[10], refltitle[10];
    char refnam[750000]	/* was [25000][10] */, reffile[1200], reftitle[1200];
} refer_;

#define refer_1 refer_

struct {
    doublereal tors1[400000]	/* was [4][100000] */, tors2[400000]	/* 
	    was [4][100000] */, tors3[400000]	/* was [4][100000] */, tors4[
	    400000]	/* was [4][100000] */, tors5[400000]	/* was [4][
	    100000] */, tors6[400000]	/* was [4][100000] */;
    integer ntors, itors[400000]	/* was [4][100000] */;
} tors_;

#define tors_1 tors_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__1000 = 1000;
static integer c__5 = 5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  program distgeom  --  produce distance geometry structures  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "distgeom" uses a metric matrix distance geometry procedure to */
/*     generate structures with interpoint distances that lie within */
/*     specified bounds, with chiral centers that maintain chirality, */
/*     and with torsional angles restrained to desired values; the */
/*     user also has the ability to interactively inspect and alter */
/*     the triangle smoothed bounds matrix prior to embedding */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 DISTGEOM  --  Too many Distance Geometry"
	    " Atoms;\002,\002 Increase MAXGEO\002)";
    static char fmt_30[] = "(/,\002 Number of Distance Geometry Structure"
	    "s\002,\002 Desired [1] :  \002,$)";
    static char fmt_40[] = "(i10)";
    static char fmt_50[] = "(/,\002 Impose Chirality Constraints on Tetrahed"
	    "ral\002,\002 Atoms [Y] :  \002,$)";
    static char fmt_60[] = "(a120)";
    static char fmt_70[] = "(/,\002 Use \"Floating\" Chirality for -XH2- and"
	    " -XH3\002,\002 Groups [N] :  \002,$)";
    static char fmt_80[] = "(a120)";
    static char fmt_90[] = "(/,\002 Impose Planarity and/or Chirality of Tri"
	    "gonal\002,\002 Atoms [Y] :  \002,$)";
    static char fmt_100[] = "(a120)";
    static char fmt_110[] = "(/,\002 Impose Torsional Planarity on Adjacent "
	    "Trigonal\002,\002 Atoms [Y] :  \002,$)";
    static char fmt_120[] = "(a120)";
    static char fmt_130[] = "(/,\002 Do You Wish to Examine or Alter the Bou"
	    "nds\002,\002 Matrix [N] :  \002,$)";
    static char fmt_140[] = "(a120)";
    static char fmt_150[] = "(/,\002 Select the Enantiomer Closest to the In"
	    "put\002,\002 Structure [N] :  \002,$)";
    static char fmt_160[] = "(a120)";
    static char fmt_170[] = "(/,\002 Refinement via Minimization or Anneal"
	    "ing\002,\002 [M or A, <CR>=A] :  \002,$)";
    static char fmt_180[] = "(a120)";
    static char fmt_190[] = "(/,\002 Interatomic Distance Bound Restraints "
	    ":\002,//,12x,\002Atom Numbers\002,7x,\002LowerBound\002,4x,\002U"
	    "pperBound\002,7x,\002Weight\002,/)";
    static char fmt_200[] = "(i6,5x,2i6,3x,2f14.4)";
    static char fmt_210[] = "(i6,5x,2i6,3x,3f14.4)";
    static char fmt_220[] = "(/,\002 Intramolecular Torsional Angle Restrain"
	    "ts :\002,//,18x,\002Atom Numbers\002,16x,\002Torsion Range\002,9"
	    "x,\002Weight\002,/)";
    static char fmt_230[] = "(i6,5x,4i6,3x,3f12.4)";
    static char fmt_250[] = "(/,\002 Bounds Smoothing via Triangle and Inver"
	    "se\002,\002 Triangle Inequality :\002)";
    static char fmt_260[] = "(/,\002 Time Required for Bounds Smoothing :"
	    "\002,4x,f12.2,\002 seconds\002)";
    static char fmt_270[] = "(/,\002 Enter an Atom Pair to Display Bounds"
	    "\002,\002 [<CR> When Done] :  \002,$)";
    static char fmt_280[] = "(a120)";
    static char fmt_290[] = "(/,\002 Lower Bound :\002,f8.3,8x,\002Upper Bou"
	    "nd :\002,f8.3)";
    static char fmt_310[] = "(/,\002 Enter New Bounds or <CR> to Leave Uncha"
	    "nged :  \002,$)";
    static char fmt_320[] = "(a120)";
    static char fmt_340[] = "(/,\002 Largest Upper Bound Distance :\002,15x,"
	    "f15.4)";
    static char fmt_350[] = "(/,\002 DISTGEOM  --  Atom\002,i6,\002 has no D"
	    "istance\002,\002 Constraints\002)";
    static char fmt_370[] = "(/,\002 Generation via Distance Geometry of Str"
	    "ucture\002,i5)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2, i__3, i__4[3];
    doublereal d__1, d__2, d__3;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_cmp(char *, char *, ftnlen, 
	    ftnlen), s_rsli(icilist *), do_lio(integer *, integer *, char *, 
	    ftnlen), e_rsli(void), s_rsfe(cilist *), do_fio(integer *, char *,
	     ftnlen), e_rsfe(void);
    double sqrt(doublereal), acos(doublereal), cos(doublereal);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int geodesic_(void);
    extern integer freeunit_(void);
    static doublereal rmsvalue, uppermax;
    extern /* Subroutine */ int torsions_(void);
    static integer i__, j, k, b1, b2, r1, r2, r3, r4;
    static doublereal t1, t2;
    static integer ia, ja, kb, ib, ic, id;
    static doublereal rab, rac, rbc, qab, qbc, qcd, qac, qbd;
    static char ext[7];
    static doublereal big1, big2, radi;
    static integer igeo, ngeo;
    static logical done;
    static doublereal secs;
    static logical info;
    static integer swap, lext, next;
    static logical quit;
    extern /* Subroutine */ int embed_(void);
    static doublereal angle;
    extern /* Subroutine */ int fatal_(void), final_(void);
    static doublereal weigh;
    static char title[120];
    extern /* Subroutine */ int bonds_(void), kgeom_(void);
    static logical exist, query;
    static doublereal hbond1, hbond2, bndfac, angfac;
    static logical header;
    static doublereal cosabc;
    extern /* Subroutine */ int grafic_(integer *, integer *, doublereal *, 
	    char *, ftnlen);
    static doublereal cosbcd, sinabc, sinbcd;
    extern /* Subroutine */ int attach_(void);
    static doublereal bndmin, bndmax;
    extern /* Subroutine */ int active_(void);
    static char record[120];
    extern /* Subroutine */ int angles_(void);
    static integer nhydro;
    static doublereal tormin, tormax, cosmin, cosmax;
    static char answer[1], letter[1], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), setime_(void), 
	    getime_(doublereal *), trifix_(integer *, integer *), impose_(
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), getxyz_(
	    void), prtxyz_(integer *);
    static char geofile[120];
    extern /* Subroutine */ int makeref_(integer *), kchiral_(void), initial_(
	    void), numeral_(integer *, char *, integer *, ftnlen), nextarg_(
	    char *, logical *, ftnlen), gettext_(char *, char *, integer *, 
	    ftnlen, ftnlen), version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };
    static icilist io___10 = { 1, string, 1, 0, 120, 1 };
    static cilist io___11 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___49 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___78 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___80 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___82 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___83 = { 0, 0, 0, fmt_280, 0 };
    static icilist io___84 = { 1, record, 1, 0, 120, 1 };
    static cilist io___88 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___89 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___90 = { 0, 0, 0, fmt_320, 0 };
    static icilist io___91 = { 1, record, 1, 0, 120, 1 };
    static cilist io___92 = { 0, 0, 0, fmt_340, 0 };
    static cilist io___94 = { 0, 0, 0, fmt_350, 0 };
    static cilist io___95 = { 0, 0, 0, fmt_370, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define bnd_ref(a_1,a_2) disgeo_1.bnd[(a_2)*1000 + a_1 - 1001]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define dfix_ref(a_1,a_2) kgeoms_1.dfix[(a_2)*3 + a_1 - 4]
#define tfix_ref(a_1,a_2) kgeoms_1.tfix[(a_2)*3 + a_1 - 4]
#define ichir_ref(a_1,a_2) kgeoms_1.ichir[(a_2)*4 + a_1 - 5]
#define idfix_ref(a_1,a_2) kgeoms_1.idfix[(a_2)*2 + a_1 - 3]
#define itfix_ref(a_1,a_2) kgeoms_1.itfix[(a_2)*4 + a_1 - 5]
#define itors_ref(a_1,a_2) tors_1.itors[(a_2)*4 + a_1 - 5]



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
/*     ##  angle.i  --  bond angles within the current structure  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     ak       harmonic angle force constant (kcal/mole/rad**2) */
/*     anat     ideal bond angle or phase shift angle (degrees) */
/*     afld     periodicity for Fourier bond angle term */
/*     nangle   total number of bond angles in the system */
/*     iang     numbers of the atoms in each bond angle */




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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  disgeo.i  --  distance geometry bounds and parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     bnd         distance geometry upper and lower bounds matrix */
/*     vdwrad      hard sphere radii for distance geometry atoms */
/*     vdwmax      maximum value of hard sphere sum for an atom pair */
/*     compact     index of local distance compaction on embedding */
/*     pathmax     maximum value of upper bound after smoothing */
/*     use_invert  flag to use enantiomer closest to input structure */
/*     use_anneal  flag to use simulated annealing refinement */




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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  kgeoms.i  --  parameters for the geometrical restraints  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xpfix      x-coordinate target for each restrained position */
/*     ypfix      y-coordinate target for each restrained position */
/*     zpfix      z-coordinate target for each restrained position */
/*     pfix       force constant and flat-well range for each position */
/*     dfix       force constant and target range for each distance */
/*     afix       force constant and target range for each angle */
/*     tfix       force constant and target range for each torsion */
/*     gfix       force constant and target range for each group distance */
/*     chir       force constant and target range for chiral centers */
/*     depth      depth of shallow Gaussian basin restraint */
/*     width      exponential width coefficient of Gaussian basin */
/*     rwall      radius of spherical droplet boundary restraint */
/*     npfix      number of position restraints to be applied */
/*     ipfix      atom number involved in each position restraint */
/*     kpfix      flags to use x-, y-, z-coordinate position restraints */
/*     ndfix      number of distance restraints to be applied */
/*     idfix      atom numbers defining each distance restraint */
/*     nafix      number of angle restraints to be applied */
/*     iafix      atom numbers defining each angle restraint */
/*     ntfix      number of torsional restraints to be applied */
/*     itfix      atom numbers defining each torsional restraint */
/*     ngfix      number of group distance restraints to be applied */
/*     igfix      group numbers defining each group distance restraint */
/*     nchir      number of chirality restraints to be applied */
/*     ichir      atom numbers defining each chirality restraint */
/*     use_basin  logical flag governing use of Gaussian basin */
/*     use_wall   logical flag governing use of droplet boundary */




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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  refer.i  --  storage of reference atomic coordinate set  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xref        reference x-coordinates for atoms in each system */
/*     yref        reference y-coordinates for atoms in each system */
/*     zref        reference z-coordinates for atoms in each system */
/*     nref        total number of atoms in each reference system */
/*     reftyp      atom types of the atoms in each reference system */
/*     n12ref      number of atoms bonded to each reference atom */
/*     i12ref      atom numbers of atoms 1-2 connected to each atom */
/*     refleng     length in characters of each reference filename */
/*     refltitle   length in characters of each reference title line */
/*     refnam      atom names of the atoms in each reference system */
/*     reffile     base filename for each reference system */
/*     reftitle    title used to describe each reference system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  tors.i  --  torsional angles within the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     tors1   1-fold amplitude and phase for each torsional angle */
/*     tors2   2-fold amplitude and phase for each torsional angle */
/*     tors3   3-fold amplitude and phase for each torsional angle */
/*     tors4   4-fold amplitude and phase for each torsional angle */
/*     tors5   5-fold amplitude and phase for each torsional angle */
/*     tors6   6-fold amplitude and phase for each torsional angle */
/*     ntors   total number of torsional angles in the system */
/*     itors   numbers of the atoms in each torsional angle */




/*     get the input structure file for the embedding */

    initial_();
    getxyz_();

/*     quit if there are too many atoms for distance geometry */

    if (atoms_1.n > 1000) {
	io___1.ciunit = iounit_1.iout;
	s_wsfe(&io___1);
	e_wsfe();
	fatal_();
    }

/*     set the lists of attached atoms and local interactions */

    attach_();
    active_();
    bonds_();
    angles_();
    torsions_();

/*     get distance bound and torsional angle restraints */

    kgeom_();

/*     store the input structure for later comparison */

    makeref_(&c__1);

/*     assign approximate radii to each of the atoms */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)letter = *(unsigned char *)name___ref(0, i__);
	if (s_cmp(name___ref(0, i__), "CH ", (ftnlen)3, (ftnlen)3) == 0) {
	    disgeo_1.vdwrad[i__ - 1] = 1.5;
	} else if (s_cmp(name___ref(0, i__), "CH2", (ftnlen)3, (ftnlen)3) == 
		0) {
	    disgeo_1.vdwrad[i__ - 1] = 1.6;
	} else if (s_cmp(name___ref(0, i__), "CH3", (ftnlen)3, (ftnlen)3) == 
		0) {
	    disgeo_1.vdwrad[i__ - 1] = 1.7;
	} else if (*(unsigned char *)letter == 'H') {
	    disgeo_1.vdwrad[i__ - 1] = .95;
	} else if (*(unsigned char *)letter == 'C') {
	    disgeo_1.vdwrad[i__ - 1] = 1.45;
	} else if (*(unsigned char *)letter == 'N') {
	    disgeo_1.vdwrad[i__ - 1] = 1.35;
	} else if (*(unsigned char *)letter == 'O') {
	    disgeo_1.vdwrad[i__ - 1] = 1.35;
	} else if (*(unsigned char *)letter == 'P') {
	    disgeo_1.vdwrad[i__ - 1] = 1.8;
	} else if (*(unsigned char *)letter == 'S') {
	    disgeo_1.vdwrad[i__ - 1] = 1.8;
	} else {
	    disgeo_1.vdwrad[i__ - 1] = .5;
	}
    }

/*     find maximum value of vdw radii sum for an atom pair */

    big1 = 0.;
    big2 = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	radi = disgeo_1.vdwrad[i__ - 1];
	if (radi > big1) {
	    big2 = big1;
	    big1 = radi;
	} else if (radi > big2) {
	    big2 = radi;
	}
    }
    disgeo_1.vdwmax = big1 + big2;

/*     set number of distance geometry structures to generate */

    ngeo = -1;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___10);
	if (i__1 != 0) {
	    goto L20;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ngeo, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L20;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L20;
	}
    }
L20:
    if (ngeo <= 0) {
	io___11.ciunit = iounit_1.iout;
	s_wsfe(&io___11);
	e_wsfe();
	io___12.ciunit = iounit_1.input;
	s_rsfe(&io___12);
	do_fio(&c__1, (char *)&ngeo, (ftnlen)sizeof(integer));
	e_rsfe();
	if (ngeo <= 0) {
	    ngeo = 1;
	}
    }

/*     enforce the original chirality of tetravalent atoms */

    kgeoms_1.nchir = 0;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___14.ciunit = iounit_1.iout;
	s_wsfe(&io___14);
	e_wsfe();
	io___15.ciunit = iounit_1.input;
	s_rsfe(&io___15);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer != 'N') {
	nextarg_(answer, &exist, (ftnlen)1);
	if (! exist) {
	    io___18.ciunit = iounit_1.iout;
	    s_wsfe(&io___18);
	    e_wsfe();
	    io___19.ciunit = iounit_1.input;
	    s_rsfe(&io___19);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    next = 1;
	    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
	}
	upcase_(answer, (ftnlen)1);
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (couple_1.n12[i__ - 1] == 4) {
		nhydro = 0;
		if (*(unsigned char *)answer == 'Y') {
		    for (j = 1; j <= 4; ++j) {
			*(unsigned char *)letter = *(unsigned char *)
				name___ref(0, i12_ref(j, i__));
			if (*(unsigned char *)letter == 'H') {
			    ++nhydro;
			}
		    }
		}
		if (nhydro < 2) {
		    ++kgeoms_1.nchir;
		    for (j = 1; j <= 4; ++j) {
			ichir_ref(j, kgeoms_1.nchir) = i12_ref(j, i__);
		    }
		}
	    }
	}
    }

/*     enforce the planarity or chirality of trigonal centers */

    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___22.ciunit = iounit_1.iout;
	s_wsfe(&io___22);
	e_wsfe();
	io___23.ciunit = iounit_1.input;
	s_rsfe(&io___23);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer != 'N') {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (couple_1.n12[i__ - 1] == 3) {
		++kgeoms_1.nchir;
		for (j = 1; j <= 3; ++j) {
		    ichir_ref(j, kgeoms_1.nchir) = i12_ref(j, i__);
		}
		ichir_ref(4, kgeoms_1.nchir) = i__;
	    }
	}
    }

/*     enforce torsional planarity on adjacent trigonal sites */

    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___24.ciunit = iounit_1.iout;
	s_wsfe(&io___24);
	e_wsfe();
	io___25.ciunit = iounit_1.input;
	s_rsfe(&io___25);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer != 'N') {
	i__1 = bond_1.nbond;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = ibnd_ref(1, i__);
	    ib = ibnd_ref(2, i__);
	    if (couple_1.n12[ia - 1] == 3 && couple_1.n12[ib - 1] == 3) {
		i__2 = couple_1.n12[ia - 1];
		for (j = 1; j <= i__2; ++j) {
		    ja = i12_ref(j, ia);
		    i__3 = couple_1.n12[ib - 1];
		    for (k = 1; k <= i__3; ++k) {
			kb = i12_ref(k, ib);
			if (ja != ib && kb != ia) {
			    ++kgeoms_1.nchir;
			    ichir_ref(1, kgeoms_1.nchir) = ja;
			    ichir_ref(2, kgeoms_1.nchir) = kb;
			    ichir_ref(3, kgeoms_1.nchir) = ia;
			    ichir_ref(4, kgeoms_1.nchir) = ib;
			}
		    }
		}
	    }
	}
    }

/*     optionally inspect and alter the smoothed bounds matrix */

    query = FALSE_;
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
	query = TRUE_;
    }

/*     choose the global enantiomer nearest to the original */

    disgeo_1.use_invert__ = FALSE_;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___34.ciunit = iounit_1.iout;
	s_wsfe(&io___34);
	e_wsfe();
	io___35.ciunit = iounit_1.input;
	s_rsfe(&io___35);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'Y') {
	disgeo_1.use_invert__ = TRUE_;
    }

/*     set the type of refinement to be used after embedding */

    disgeo_1.use_anneal__ = TRUE_;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___36.ciunit = iounit_1.iout;
	s_wsfe(&io___36);
	e_wsfe();
	io___37.ciunit = iounit_1.input;
	s_rsfe(&io___37);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'M') {
	disgeo_1.use_anneal__ = FALSE_;
    }

/*     initialize chirality and planarity restraint values */

    kchiral_();

/*     change the default distance restraint force constant */

    i__1 = kgeoms_1.ndfix;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (dfix_ref(1, i__) == 100.) {
	    dfix_ref(1, i__) = 1.;
	}
    }

/*     print lists of the interatomic distance restraints */

    if (inform_1.verbose) {
	header = TRUE_;
	i__1 = kgeoms_1.ndfix;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = idfix_ref(1, i__);
	    ib = idfix_ref(2, i__);
	    weigh = dfix_ref(1, i__);
	    bndmin = dfix_ref(2, i__);
	    bndmax = dfix_ref(3, i__);
	    if (header) {
		header = FALSE_;
		io___42.ciunit = iounit_1.iout;
		s_wsfe(&io___42);
		e_wsfe();
	    }
	    if (weigh == 1.) {
		io___43.ciunit = iounit_1.iout;
		s_wsfe(&io___43);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&bndmin, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bndmax, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___44.ciunit = iounit_1.iout;
		s_wsfe(&io___44);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&bndmin, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bndmax, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&weigh, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}

/*     print lists of the torsional angle restraints */

	header = TRUE_;
	i__1 = kgeoms_1.ntfix;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = itfix_ref(1, i__);
	    ib = itfix_ref(2, i__);
	    ic = itfix_ref(3, i__);
	    id = itfix_ref(4, i__);
	    weigh = tfix_ref(1, i__);
	    tormin = tfix_ref(2, i__);
	    tormax = tfix_ref(3, i__);
	    if (header) {
		header = FALSE_;
		io___49.ciunit = iounit_1.iout;
		s_wsfe(&io___49);
		e_wsfe();
	    }
	    io___50.ciunit = iounit_1.iout;
	    s_wsfe(&io___50);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&tormin, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&tormax, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&weigh, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     initialize the upper and lower bounds matrix */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	bnd_ref(i__, i__) = 0.;
    }
    uppermax = 1e6;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    bnd_ref(j, i__) = uppermax;
	}
    }
    i__1 = atoms_1.n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	radi = disgeo_1.vdwrad[i__ - 1];
	i__2 = atoms_1.n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    bnd_ref(j, i__) = radi + disgeo_1.vdwrad[j - 1];
	}
    }

/*     set the upper and lower bounds for 1-2 distances */

    bndfac = 0.;
/*     bndfac = 0.01d0 */
    i__1 = bond_1.nbond;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ibnd_ref(1, i__);
	ib = ibnd_ref(2, i__);
/* Computing 2nd power */
	d__1 = atoms_1.x[ia - 1] - atoms_1.x[ib - 1];
/* Computing 2nd power */
	d__2 = atoms_1.y[ia - 1] - atoms_1.y[ib - 1];
/* Computing 2nd power */
	d__3 = atoms_1.z__[ia - 1] - atoms_1.z__[ib - 1];
	rab = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	bndmin = (1. - bndfac) * rab;
	bndmax = (bndfac + 1.) * rab;
	if (ia > ib) {
	    bnd_ref(ia, ib) = bndmin;
	    bnd_ref(ib, ia) = bndmax;
	} else {
	    bnd_ref(ia, ib) = bndmax;
	    bnd_ref(ib, ia) = bndmin;
	}
    }

/*     set the upper and lower bounds for 1-3 distances */

    angfac = 0.;
/*     angfac = 0.01d0 */
    i__1 = angle_1.nangle;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
/* Computing 2nd power */
	d__1 = atoms_1.x[ia - 1] - atoms_1.x[ib - 1];
/* Computing 2nd power */
	d__2 = atoms_1.y[ia - 1] - atoms_1.y[ib - 1];
/* Computing 2nd power */
	d__3 = atoms_1.z__[ia - 1] - atoms_1.z__[ib - 1];
	rab = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
	d__1 = atoms_1.x[ib - 1] - atoms_1.x[ic - 1];
/* Computing 2nd power */
	d__2 = atoms_1.y[ib - 1] - atoms_1.y[ic - 1];
/* Computing 2nd power */
	d__3 = atoms_1.z__[ib - 1] - atoms_1.z__[ic - 1];
	rbc = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
	d__1 = atoms_1.x[ia - 1] - atoms_1.x[ic - 1];
/* Computing 2nd power */
	d__2 = atoms_1.y[ia - 1] - atoms_1.y[ic - 1];
/* Computing 2nd power */
	d__3 = atoms_1.z__[ia - 1] - atoms_1.z__[ic - 1];
	rac = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
	d__1 = rab;
/* Computing 2nd power */
	d__2 = rbc;
/* Computing 2nd power */
	d__3 = rac;
	angle = acos((d__1 * d__1 + d__2 * d__2 - d__3 * d__3) / (rab * 2. * 
		rbc));
	cosmin = cos(angle * (1. - angfac));
/* Computing MIN */
	d__1 = 3.141592653589793238, d__2 = angle * (angfac + 1.);
	cosmax = cos((min(d__1,d__2)));
/* Computing MIN */
	d__1 = bnd_ref(ia, ib), d__2 = bnd_ref(ib, ia);
	qab = min(d__1,d__2);
/* Computing MIN */
	d__1 = bnd_ref(ic, ib), d__2 = bnd_ref(ib, ic);
	qbc = min(d__1,d__2);
/* Computing 2nd power */
	d__1 = qab;
/* Computing 2nd power */
	d__2 = qbc;
	bndmin = d__1 * d__1 + d__2 * d__2 - qab * 2. * qbc * cosmin;
	bndmin = sqrt((max(0.,bndmin)));
/* Computing MAX */
	d__1 = bnd_ref(ia, ib), d__2 = bnd_ref(ib, ia);
	qab = max(d__1,d__2);
/* Computing MAX */
	d__1 = bnd_ref(ic, ib), d__2 = bnd_ref(ib, ic);
	qbc = max(d__1,d__2);
/* Computing 2nd power */
	d__1 = qab;
/* Computing 2nd power */
	d__2 = qbc;
	bndmax = d__1 * d__1 + d__2 * d__2 - qab * 2. * qbc * cosmax;
	bndmax = sqrt((max(0.,bndmax)));
	if (ia > ic) {
	    bnd_ref(ia, ic) = bndmin;
	    bnd_ref(ic, ia) = bndmax;
	} else {
	    bnd_ref(ia, ic) = bndmax;
	    bnd_ref(ic, ia) = bndmin;
	}
    }

/*     set the upper and lower bounds for 1-4 distances */

    i__1 = tors_1.ntors;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = itors_ref(1, i__);
	ib = itors_ref(2, i__);
	ic = itors_ref(3, i__);
	id = itors_ref(4, i__);
	cosmin = 1.;
	cosmax = -1.;
	i__2 = kgeoms_1.ntfix;
	for (j = 1; j <= i__2; ++j) {
	    r1 = itfix_ref(1, j);
	    r2 = itfix_ref(2, j);
	    r3 = itfix_ref(3, j);
	    r4 = itfix_ref(4, j);
	    if (ia == r1 && ib == r2 && ic == r3 && id == r4 || ia == r4 && 
		    ib == r3 && ic == r2 && id == r1) {
		t1 = tfix_ref(2, j) / 57.29577951308232088;
		t2 = tfix_ref(3, j) / 57.29577951308232088;
		if (t2 >= 0. && t1 <= 0.) {
		    cosmin = 1.;
/* Computing MIN */
		    d__1 = cos(t1), d__2 = cos(t2);
		    cosmax = min(d__1,d__2);
		} else if (t1 >= 0. && t2 <= 0.) {
/* Computing MAX */
		    d__1 = cos(t1), d__2 = cos(t2);
		    cosmin = max(d__1,d__2);
		    cosmax = -1.;
		} else if (t1 >= 0. && t2 >= t1) {
		    cosmin = cos(t1);
		    cosmax = cos(t2);
		} else if (t2 <= 0. && t1 <= t2) {
		    cosmin = cos(t2);
		    cosmax = cos(t1);
		}
		goto L240;
	    }
	}
L240:
/* Computing MIN */
	d__1 = bnd_ref(ia, ib), d__2 = bnd_ref(ib, ia);
	qab = min(d__1,d__2);
/* Computing MIN */
	d__1 = bnd_ref(ib, ic), d__2 = bnd_ref(ic, ib);
	qbc = min(d__1,d__2);
/* Computing MIN */
	d__1 = bnd_ref(ic, id), d__2 = bnd_ref(id, ic);
	qcd = min(d__1,d__2);
/* Computing MIN */
	d__1 = bnd_ref(ia, ic), d__2 = bnd_ref(ic, ia);
	qac = min(d__1,d__2);
/* Computing MIN */
	d__1 = bnd_ref(ib, id), d__2 = bnd_ref(id, ib);
	qbd = min(d__1,d__2);
/* Computing 2nd power */
	d__1 = qab;
/* Computing 2nd power */
	d__2 = qbc;
/* Computing 2nd power */
	d__3 = qac;
	cosabc = (d__1 * d__1 + d__2 * d__2 - d__3 * d__3) / (qab * 2. * qbc);
/* Computing MAX */
/* Computing 2nd power */
	d__3 = cosabc;
	d__1 = 0., d__2 = 1. - d__3 * d__3;
	sinabc = sqrt((max(d__1,d__2)));
/* Computing 2nd power */
	d__1 = qbc;
/* Computing 2nd power */
	d__2 = qcd;
/* Computing 2nd power */
	d__3 = qbd;
	cosbcd = (d__1 * d__1 + d__2 * d__2 - d__3 * d__3) / (qbc * 2. * qcd);
/* Computing MAX */
/* Computing 2nd power */
	d__3 = cosbcd;
	d__1 = 0., d__2 = 1. - d__3 * d__3;
	sinbcd = sqrt((max(d__1,d__2)));
/* Computing 2nd power */
	d__1 = qab;
/* Computing 2nd power */
	d__2 = qbc;
/* Computing 2nd power */
	d__3 = qcd;
	bndmin = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + qab * 2. * qcd * 
		cosabc * cosbcd - qab * 2. * qcd * sinabc * sinbcd * cosmin - 
		qbc * 2. * (qab * cosabc + qcd * cosbcd);
	bndmin = sqrt((max(0.,bndmin)));
/* Computing MAX */
	d__1 = bnd_ref(ia, ib), d__2 = bnd_ref(ib, ia);
	qab = max(d__1,d__2);
/* Computing MAX */
	d__1 = bnd_ref(ib, ic), d__2 = bnd_ref(ic, ib);
	qbc = max(d__1,d__2);
/* Computing MAX */
	d__1 = bnd_ref(ic, id), d__2 = bnd_ref(id, ic);
	qcd = max(d__1,d__2);
/* Computing MAX */
	d__1 = bnd_ref(ia, ic), d__2 = bnd_ref(ic, ia);
	qac = max(d__1,d__2);
/* Computing MAX */
	d__1 = bnd_ref(ib, id), d__2 = bnd_ref(id, ib);
	qbd = max(d__1,d__2);
/* Computing 2nd power */
	d__1 = qab;
/* Computing 2nd power */
	d__2 = qbc;
/* Computing 2nd power */
	d__3 = qac;
	cosabc = (d__1 * d__1 + d__2 * d__2 - d__3 * d__3) / (qab * 2. * qbc);
/* Computing MAX */
/* Computing 2nd power */
	d__3 = cosabc;
	d__1 = 0., d__2 = 1. - d__3 * d__3;
	sinabc = sqrt((max(d__1,d__2)));
/* Computing 2nd power */
	d__1 = qbc;
/* Computing 2nd power */
	d__2 = qcd;
/* Computing 2nd power */
	d__3 = qbd;
	cosbcd = (d__1 * d__1 + d__2 * d__2 - d__3 * d__3) / (qbc * 2. * qcd);
/* Computing MAX */
/* Computing 2nd power */
	d__3 = cosbcd;
	d__1 = 0., d__2 = 1. - d__3 * d__3;
	sinbcd = sqrt((max(d__1,d__2)));
/* Computing 2nd power */
	d__1 = qab;
/* Computing 2nd power */
	d__2 = qbc;
/* Computing 2nd power */
	d__3 = qcd;
	bndmax = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + qab * 2. * qcd * 
		cosabc * cosbcd - qab * 2. * qcd * sinabc * sinbcd * cosmax - 
		qbc * 2. * (qab * cosabc + qcd * cosbcd);
	bndmax = sqrt((max(0.,bndmax)));
	if (ia > id) {
	    bnd_ref(ia, id) = bndmin;
	    bnd_ref(id, ia) = bndmax;
	} else {
	    bnd_ref(ia, id) = bndmax;
	    bnd_ref(id, ia) = bndmin;
	}
    }

/*     convert distance restraints into bounds matrix elements */

    i__1 = kgeoms_1.ndfix;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = idfix_ref(1, i__);
	ib = idfix_ref(2, i__);
	bndmin = dfix_ref(2, i__);
	bndmax = dfix_ref(3, i__);
	if (ia > ib) {
	    bnd_ref(ia, ib) = bndmin;
	    bnd_ref(ib, ia) = bndmax;
	} else {
	    bnd_ref(ia, ib) = bndmax;
	    bnd_ref(ib, ia) = bndmin;
	}
    }

/*     modify lower bounds to allow hydrogen bond formation */

    hbond1 = 1.7;
    hbond2 = 2.55;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)letter = *(unsigned char *)name___ref(0, i__);
	if (*(unsigned char *)letter == 'N' || *(unsigned char *)letter == 
		'O') {
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		*(unsigned char *)letter = *(unsigned char *)name___ref(0, j);
		if (*(unsigned char *)letter == 'H') {
		    k = i12_ref(1, j);
		    *(unsigned char *)letter = *(unsigned char *)name___ref(0,
			     k);
		    if (*(unsigned char *)letter == 'N' || *(unsigned char *)
			    letter == 'O') {
			if (j > i__) {
/* Computing MIN */
			    d__1 = hbond1, d__2 = bnd_ref(j, i__);
			    bnd_ref(j, i__) = min(d__1,d__2);
			} else {
/* Computing MIN */
			    d__1 = hbond1, d__2 = bnd_ref(i__, j);
			    bnd_ref(i__, j) = min(d__1,d__2);
			}
			if (k > i__) {
/* Computing MIN */
			    d__1 = hbond2, d__2 = bnd_ref(k, i__);
			    bnd_ref(k, i__) = min(d__1,d__2);
			} else {
/* Computing MIN */
			    d__1 = hbond2, d__2 = bnd_ref(i__, k);
			    bnd_ref(i__, k) = min(d__1,d__2);
			}
		    }
		}
	    }
	}
    }

/*     use the triangle inequalities to smooth the bounds */

    if (inform_1.verbose && atoms_1.n <= 130) {
	s_copy(title, "Input Distance Bounds :", (ftnlen)120, (ftnlen)23);
	grafic_(&atoms_1.n, &c__1000, disgeo_1.bnd, title, (ftnlen)120);
    }
    io___78.ciunit = iounit_1.iout;
    s_wsfe(&io___78);
    e_wsfe();
    if (inform_1.verbose) {
	setime_();
    }
    geodesic_();
/*     call triangle */
    if (inform_1.verbose) {
	getime_(&secs);
	io___80.ciunit = iounit_1.iout;
	s_wsfe(&io___80);
	do_fio(&c__1, (char *)&secs, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     allow interactive alteration of the bounds matrix */

    done = FALSE_;
    while(query && ! done) {
	done = TRUE_;
	io___82.ciunit = iounit_1.iout;
	s_wsfe(&io___82);
	e_wsfe();
	io___83.ciunit = iounit_1.input;
	s_rsfe(&io___83);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___84);
	if (i__1 != 0) {
	    goto L330;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&b1, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L330;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&b2, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L330;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L330;
	}
	done = FALSE_;
	if (b1 < 1 || b2 > atoms_1.n || b1 == b2) {
	    goto L330;
	}
	if (b1 > b2) {
	    swap = b1;
	    b1 = b2;
	    b2 = swap;
	}
	io___88.ciunit = iounit_1.iout;
	s_wsfe(&io___88);
	do_fio(&c__1, (char *)&bnd_ref(b2, b1), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&bnd_ref(b1, b2), (ftnlen)sizeof(doublereal));
	e_wsfe();
L300:
	io___89.ciunit = iounit_1.iout;
	s_wsfe(&io___89);
	e_wsfe();
	io___90.ciunit = iounit_1.input;
	s_rsfe(&io___90);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___91);
	if (i__1 != 0) {
	    goto L330;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&bndmin, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L330;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&bndmax, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L330;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L330;
	}
	if (bndmin > bndmax) {
	    goto L300;
	}
	bnd_ref(b2, b1) = bndmin;
	bnd_ref(b1, b2) = bndmax;
	trifix_(&b1, &b2);
L330:
	;
    }

/*     display the smoothed upper and lower bounds matrix */

    if (inform_1.verbose && atoms_1.n <= 130) {
	s_copy(title, "Triangle Smoothed Bounds :", (ftnlen)120, (ftnlen)26);
	grafic_(&atoms_1.n, &c__1000, disgeo_1.bnd, title, (ftnlen)120);
    }

/*     find the largest value of an upper bound between atoms */

    disgeo_1.pathmax = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    if (disgeo_1.pathmax < bnd_ref(j, i__)) {
		disgeo_1.pathmax = bnd_ref(j, i__);
	    }
	}
    }
    io___92.ciunit = iounit_1.iout;
    s_wsfe(&io___92);
    do_fio(&c__1, (char *)&disgeo_1.pathmax, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     check for any atoms that have no distance restraints */

    quit = FALSE_;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    if (bnd_ref(j, i__) != uppermax) {
		goto L360;
	    }
	}
	i__2 = atoms_1.n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    if (bnd_ref(i__, j) != uppermax) {
		goto L360;
	    }
	}
	quit = TRUE_;
	io___94.ciunit = iounit_1.iout;
	s_wsfe(&io___94);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
L360:
	;
    }
    if (quit) {
	fatal_();
    }

/*     generate the desired number of distance geometry structures */

    i__1 = ngeo;
    for (j = 1; j <= i__1; ++j) {
	io___95.ciunit = iounit_1.iout;
	s_wsfe(&io___95);
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	e_wsfe();
	embed_();

/*     superpose the distance geometry structure on input structure */

	info = inform_1.verbose;
	inform_1.verbose = FALSE_;
	impose_(&atoms_1.n, refer_1.xref, refer_1.yref, refer_1.zref, &
		atoms_1.n, atoms_1.x, atoms_1.y, atoms_1.z__, &rmsvalue);
	inform_1.verbose = info;

/*     write out the final optimized distance geometry structure */

	lext = 3;
	numeral_(&j, ext, &lext, (ftnlen)7);
/* Writing concatenation */
	i__4[0] = files_1.leng, a__1[0] = files_1.filename;
	i__4[1] = 1, a__1[1] = ".";
	i__4[2] = lext, a__1[2] = ext;
	s_cat(geofile, a__1, i__4, &c__3, (ftnlen)120);
	version_(geofile, "new", (ftnlen)120, (ftnlen)3);
	igeo = freeunit_();
	o__1.oerr = 0;
	o__1.ounit = igeo;
	o__1.ofnmlen = 120;
	o__1.ofnm = geofile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	prtxyz_(&igeo);
	cl__1.cerr = 0;
	cl__1.cunit = igeo;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef itors_ref
#undef itfix_ref
#undef idfix_ref
#undef ichir_ref
#undef tfix_ref
#undef dfix_ref
#undef name___ref
#undef iang_ref
#undef ibnd_ref
#undef bnd_ref
#undef i12_ref


/* Main program alias */ int distgeom_ () { MAIN__ (); return 0; }
