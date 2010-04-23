/* valence.f -- translated by f2c (version 20050501).
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
    doublereal stpmin, stpmax, cappa, slpmax, angmax;
    integer intmax;
} linmin_;

#define linmin_1 linmin_

struct {
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

struct {
    doublereal poleps, polsor, p2scale, p3scale, p4scale, p5scale, d1scale, 
	    d2scale, d3scale, d4scale, u1scale, u2scale, u3scale, u4scale;
    char poltyp[6];
} polpot_;

#define polpot_1 polpot_

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
    doublereal gx[25000], gy[25000], gz[25000], gforce[75000]	/* was [3][
	    25000] */, gh[1000000], gfreq[1000];
    integer ngatom;
} qmstuf_;

#define qmstuf_1 qmstuf_

struct {
    logical fit_bond__, fit_angle__, fit_strbnd__, fit_urey__, fit_opbend__, 
	    fit_tors__, fit_struct__, fit_force__;
} valfit_;

#define valfit_1 valfit_

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
    doublereal bk[50000], bl[50000];
    integer nbond, ibnd[100000]	/* was [2][50000] */;
} bond_;

#define bond_1 bond_

struct {
    doublereal acon[2000], acon5[500], acon4[500], acon3[500], aconf[500], 
	    ang[6000]	/* was [3][2000] */, ang5[1500]	/* was [3][500] */, 
	    ang4[1500]	/* was [3][500] */, ang3[1500]	/* was [3][500] */, 
	    angf[1000]	/* was [2][500] */;
    char ka[24000], ka5[6000], ka4[6000], ka3[6000], kaf[6000];
} kangs_;

#define kangs_1 kangs_

struct {
    doublereal bcon[2000], blen[2000], bcon5[500], blen5[500], bcon4[500], 
	    blen4[500], bcon3[500], blen3[500], dlen[500];
    char kb[16000], kb5[4000], kb4[4000], kb3[4000], kel[6000];
} kbonds_;

#define kbonds_1 kbonds_

struct {
    doublereal opbn[500];
    char kopb[8000];
} kopbnd_;

#define kopbnd_1 kopbnd_

struct {
    doublereal stbn[4000]	/* was [2][2000] */;
    char ksb[24000];
} kstbnd_;

#define kstbnd_1 kstbnd_

struct {
    doublereal t1[4000]	/* was [2][2000] */, t2[4000]	/* was [2][2000] */, 
	    t3[4000]	/* was [2][2000] */, t4[4000]	/* was [2][2000] */, 
	    t5[4000]	/* was [2][2000] */, t6[4000]	/* was [2][2000] */, 
	    t15[1000]	/* was [2][500] */, t25[1000]	/* was [2][500] */, 
	    t35[1000]	/* was [2][500] */, t45[1000]	/* was [2][500] */, 
	    t55[1000]	/* was [2][500] */, t65[1000]	/* was [2][500] */, 
	    t14[1000]	/* was [2][500] */, t24[1000]	/* was [2][500] */, 
	    t34[1000]	/* was [2][500] */, t44[1000]	/* was [2][500] */, 
	    t54[1000]	/* was [2][500] */, t64[1000]	/* was [2][500] */;
    char kt[32000], kt5[8000], kt4[8000];
} ktorsn_;

#define ktorsn_1 ktorsn_

struct {
    doublereal ucon[2000], dst13[2000];
    char ku[24000];
} kurybr_;

#define kurybr_1 kurybr_

struct {
    doublereal rad[5000], eps[5000], rad4[5000], eps4[5000], reduct[5000];
} kvdws_;

#define kvdws_1 kvdws_

struct {
    doublereal opbk[75000];
    integer nopbend, iopb[75000];
} opbend_;

#define opbend_1 opbend_

struct {
    doublereal sbk[150000]	/* was [2][75000] */;
    integer nstrbnd, isb[225000]	/* was [3][75000] */;
} strbnd_;

#define strbnd_1 strbnd_

struct {
    doublereal tors1[400000]	/* was [4][100000] */, tors2[400000]	/* 
	    was [4][100000] */, tors3[400000]	/* was [4][100000] */, tors4[
	    400000]	/* was [4][100000] */, tors5[400000]	/* was [4][
	    100000] */, tors6[400000]	/* was [4][100000] */;
    integer ntors, itors[400000]	/* was [4][100000] */;
} tors_;

#define tors_1 tors_

struct {
    doublereal uk[75000], ul[75000];
    integer nurey, iury[225000]	/* was [3][75000] */;
} urey_;

#define urey_1 urey_

struct {
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
    logical use_vcorr__;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal cbnd, qbnd, bndunit;
    char bndtyp[8];
} bndpot_;

#define bndpot_1 bndpot_

struct {
    doublereal angunit, stbnunit, aaunit, opbunit, opdunit, cang, qang, pang, 
	    sang, copb, qopb, popb, sopb, copd, qopd, popd, sopd;
    char angtyp[600000], opbtyp[8];
} angpot_;

#define angpot_1 angpot_

struct {
    doublereal cury, qury, ureyunit;
} urypot_;

#define urypot_1 urypot_

struct {
    doublereal idihunit, itorunit, torsunit, ptorunit, storunit, ttorunit;
} torpot_;

#define torpot_1 torpot_

struct {
    doublereal hesscut;
} hescut_;

#define hescut_1 hescut_

struct {
    doublereal fctmin, hguess;
    integer maxiter, nextiter;
} minima_;

#define minima_1 minima_

struct {
    doublereal scale[75000];
    logical set_scale__;
} scales_;

#define scales_1 scales_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__0 = 0;
static doublereal c_b36 = 0.;
static integer c__2 = 2;
static integer c__4 = 4;
static integer c__1000 = 1000;
static doublereal c_b215 = 1.;



/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2009 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  program valence  --  derive valence force field parameters  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "valence" refines force field parameters for valence terms based */
/*     on a quantum mechanical optimized structure and frequencies */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 The TINKER Valence Parameter Facility Ca"
	    "n :\002,//,4x,\002(1) Set Initial Values for Valence Parameter"
	    "s\002,/,4x,\002(2) Compare QM and MM Vibrational Frequencies\002"
	    ",/,4x,\002(3) Force Fit of Parameters to QM Results\002,/,4x,"
	    "\002(4) Structure Fit of Parameters to QM Results\002)";
    static char fmt_30[] = "(/,\002 Enter the Number of the Desired Choice :"
	    "  \002,$)";
    static char fmt_40[] = "(i10)";
    static char fmt_70[] = "(/,\002 Enter RMS Gradient Termination Criterio"
	    "n\002,\002 [0.01] :  \002,$)";
    static char fmt_80[] = "(f20.0)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    static logical dotarget;
    extern /* Subroutine */ int valguess_(void), torsions_(void);
    static integer i__;
    static doublereal xx[75000];
    static integer mode;
    extern /* Subroutine */ int ocvm_(integer *, doublereal *, doublereal *, 
	    doublereal *, D_fp, U_fp);
    static integer nvar, next;
    extern /* Subroutine */ int field_(void);
    static logical dofit;
    extern /* Subroutine */ int bonds_(void), katom_(void);
    static doublereal value;
    static logical exist, query;
    extern /* Subroutine */ int attach_(void), angles_(void);
    static char record[120];
    static integer length;
    static doublereal grdmin;
    extern doublereal valrms_(integer *);
    static char string[120];
    extern /* Subroutine */ int getkey_(void), upcase_(char *, ftnlen), 
	    getxyz_(void), prmvar_(integer *, doublereal *), varprm_(integer *
	    , doublereal *, integer *, doublereal *), prtval_(void);
    extern doublereal valfit1_();
    extern /* Subroutine */ int readgau_(void), initial_(void);
    static logical doguess;
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen);
    static doublereal minimum;
    extern /* Subroutine */ int optsave_();
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static icilist io___8 = { 1, string, 1, 0, 120, 1 };
    static cilist io___9 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___11 = { 1, 0, 1, fmt_40, 0 };
    static icilist io___22 = { 1, string, 1, 0, 120, 1 };
    static cilist io___23 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_80, 0 };



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  linmin.i  --  parameters for line search minimization  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     stpmin   minimum step length in current line search direction */
/*     stpmax   maximum step length in current line search direction */
/*     cappa    stringency of line search (0=tight < cappa < 1=loose) */
/*     slpmax   projected gradient above which stepsize is reduced */
/*     angmax   maximum angle between search direction and -gradient */
/*     intmax   maximum number of interpolations during line search */




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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  polpot.i  --  specifics of polarization functional form  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     poleps    induced dipole convergence criterion (rms Debyes/atom) */
/*     polsor    induced dipole SOR convergence acceleration factor */
/*     p2scale   field 1-2 scale factor for energy evaluations */
/*     p3scale   field 1-3 scale factor for energy evaluations */
/*     p4scale   field 1-4 scale factor for energy evaluations */
/*     p5scale   field 1-5 scale factor for energy evaluations */
/*     d1scale   field intra-group scale factor for direct induced */
/*     d2scale   field 1-2 group scale factor for direct induced */
/*     d3scale   field 1-3 group scale factor for direct induced */
/*     d4scale   field 1-4 group scale factor for direct induced */
/*     u1scale   field intra-group scale factor for mutual induced */
/*     u2scale   field 1-2 group scale factor for mutual induced */
/*     u3scale   field 1-3 group scale factor for mutual induced */
/*     u4scale   field 1-4 group scale factor for mutual induced */
/*     poltyp    type of polarization potential (direct or mutual) */




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




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2009 by Chuanjie Wu and Jay William Ponder  ## */
/*     ##                    All Rights Reserved                     ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  qmstuf.i  --  quantum data from Gaussian 03 calculation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     gx       x-coordinate of each atom in the QM data file */
/*     gy       y-coordinate of each atom in the QM data file */
/*     gz       z-coordinate of each atom in the QM data file */
/*     gforce   force components on each atom from QM data */
/*     gh       Hessian maxtrix elements from QM data */
/*     gfreq    calculated vibrational frequencies from QM data */
/*     ngatom   number of atoms in the QM data file */




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  valfit.i  --  values for valence term parameter fitting  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     fit_bond    logical flag to fit bond stretch parameters */
/*     fit_angle   logical flag to fit angle bend parameters */
/*     fit_strbnd  logical flag to fit stretch-bend parameters */
/*     fit_urey    logical flag to fit Urey-Bradley parameters */
/*     fit_opbend  logical flag to fit out-of-plane bend parameters */
/*     fit_tors    logical flag to fit torsional parameters */
/*     fit_struct  logical flag to structure-fit valence parameters */
/*     fit_force   logical flag to force-fit valence parameters */




/*     initialization of the various modes of operation */

    initial_();
    valfit_1.fit_bond__ = TRUE_;
    valfit_1.fit_angle__ = TRUE_;
    valfit_1.fit_strbnd__ = FALSE_;
    valfit_1.fit_urey__ = FALSE_;
    valfit_1.fit_opbend__ = FALSE_;
    valfit_1.fit_tors__ = FALSE_;
    valfit_1.fit_force__ = FALSE_;
    valfit_1.fit_struct__ = FALSE_;
    doguess = FALSE_;
    dotarget = FALSE_;
    dofit = FALSE_;

/*     find out which valence term protocol is to be performed */

    mode = 0;
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___8);
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
	io___9.ciunit = iounit_1.iout;
	s_wsfe(&io___9);
	e_wsfe();
	while(mode < 1 || mode > 4) {
	    mode = 0;
	    io___10.ciunit = iounit_1.iout;
	    s_wsfe(&io___10);
	    e_wsfe();
	    io___11.ciunit = iounit_1.input;
	    i__1 = s_rsfe(&io___11);
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
	doguess = TRUE_;
    } else if (mode == 2) {
	dotarget = TRUE_;
    } else if (mode == 3) {
	dotarget = TRUE_;
	dofit = TRUE_;
	valfit_1.fit_force__ = TRUE_;
    } else if (mode == 4) {
	dotarget = TRUE_;
	dofit = TRUE_;
	valfit_1.fit_struct__ = TRUE_;
    }

/*     read the Cartesian coordinates and connectivity info */

    getxyz_();
    s_copy(xyzfile, files_1.filename, (ftnlen)120, (ftnlen)120);
    length = files_1.leng;

/*     read structure and vibrational data from Gaussian output */

    readgau_();
    s_copy(files_1.filename, xyzfile, (ftnlen)120, (ftnlen)120);
    files_1.leng = length;
    getkey_();

/*     assign estimated values to the valence parameters */

    if (doguess) {
	attach_();
	bonds_();
	angles_();
	torsions_();
	field_();
	katom_();
	valguess_();
    } else {
	mechanic_();
    }

/*     get control parameters and target values from keyfile */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "FIT-BOND ", (ftnlen)9, (ftnlen)9) == 0) {
	    valfit_1.fit_bond__ = TRUE_;
	} else if (s_cmp(keyword, "FIX-BOND ", (ftnlen)9, (ftnlen)9) == 0) {
	    valfit_1.fit_bond__ = FALSE_;
	} else if (s_cmp(keyword, "FIT-ANGLE ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    valfit_1.fit_angle__ = TRUE_;
	} else if (s_cmp(keyword, "FIX-ANGLE ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    valfit_1.fit_angle__ = FALSE_;
	} else if (s_cmp(keyword, "FIT-STRBND ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    valfit_1.fit_strbnd__ = TRUE_;
	} else if (s_cmp(keyword, "FIX-STRBND ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    valfit_1.fit_strbnd__ = FALSE_;
	} else if (s_cmp(keyword, "FIT-UREY ", (ftnlen)9, (ftnlen)9) == 0) {
	    valfit_1.fit_urey__ = TRUE_;
	} else if (s_cmp(keyword, "FIX-UREY ", (ftnlen)9, (ftnlen)9) == 0) {
	    valfit_1.fit_urey__ = FALSE_;
	} else if (s_cmp(keyword, "FIT-OPBEND ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    valfit_1.fit_opbend__ = TRUE_;
	} else if (s_cmp(keyword, "FIX-OPBEND ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    valfit_1.fit_opbend__ = FALSE_;
	} else if (s_cmp(keyword, "FIT-TORSION ", (ftnlen)12, (ftnlen)12) == 
		0) {
	    valfit_1.fit_tors__ = TRUE_;
	} else if (s_cmp(keyword, "FIX-TORSION ", (ftnlen)12, (ftnlen)12) == 
		0) {
	    valfit_1.fit_tors__ = FALSE_;
	}
    }

/*     try to increase robustness of polarization calculations */

    if (dofit && potent_1.use_polar__) {
	linmin_1.stpmax = 1.;
	polpot_1.polsor = .55;
    }

/*     comparison of QM and TINKER structure and frequencies */

    if (dotarget) {
	if (! dofit) {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		atoms_1.x[i__ - 1] = qmstuf_1.gx[i__ - 1];
		atoms_1.y[i__ - 1] = qmstuf_1.gy[i__ - 1];
		atoms_1.z__[i__ - 1] = qmstuf_1.gz[i__ - 1];
	    }
	    value = valrms_(&c__1);

/*     optimize the valence term force field parameters */

	} else {
	    prmvar_(&nvar, xx);
	    value = valrms_(&c__1);
	    grdmin = -1.;
	    nextarg_(string, &exist, (ftnlen)120);
	    if (exist) {
		i__1 = s_rsli(&io___22);
		if (i__1 != 0) {
		    goto L60;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&grdmin, (ftnlen)sizeof(
			doublereal));
		if (i__1 != 0) {
		    goto L60;
		}
		i__1 = e_rsli();
		if (i__1 != 0) {
		    goto L60;
		}
	    }
L60:
	    if (grdmin <= 0.) {
		io___23.ciunit = iounit_1.iout;
		s_wsfe(&io___23);
		e_wsfe();
		io___24.ciunit = iounit_1.input;
		s_rsfe(&io___24);
		do_fio(&c__1, (char *)&grdmin, (ftnlen)sizeof(doublereal));
		e_rsfe();
	    }
	    if (grdmin <= 0.) {
		grdmin = .01;
	    }
	    s_copy(output_1.coordtype, "NONE", (ftnlen)9, (ftnlen)4);
	    ocvm_(&nvar, xx, &minimum, &grdmin, (D_fp)valfit1_, (U_fp)
		    optsave_);
	    varprm_(&nvar, xx, &c__0, &c_b36);
	    prmvar_(&nvar, xx);
	    value = valrms_(&c__1);
	    prtval_();
	}
    }
    return 0;
} /* MAIN__ */

#undef keyline_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine valguess  --  estimate valence parameter values  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "valguess" sets approximate valence parameter values based on */
/*     quantum mechanical structure and frequency data */


/* Subroutine */ int valguess_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 VALENCE  --  The Number of Atoms is No"
	    "t\002,\002 Consistent\002)";
    static char fmt_20[] = "(/,\002 Estimated van der Waals Parameters :\002"
	    ",/)";
    static char fmt_30[] = "(\002 vdw\002,7x,i5,10x,f10.3,f11.4)";
    static char fmt_40[] = "(\002 vdw\002,7x,i5,10x,f10.3,f11.4,f9.2)";
    static char fmt_50[] = "(/,\002 Estimated Bond Stretching Parameters "
	    ":\002,/)";
    static char fmt_60[] = "(\002 bond\002,6x,2i5,5x,f10.1,f11.4)";
    static char fmt_70[] = "(/,\002 Estimated Angle Bending Parameters :\002"
	    ",/)";
    static char fmt_80[] = "(\002 angle\002,5x,3i5,f10.2,f11.2)";
    static char fmt_90[] = "(/,\002 Estimated Stretch-Bend Parameters :\002,"
	    "/)";
    static char fmt_100[] = "(\002 strbnd\002,4x,3i5,f10.2,f11.2)";
    static char fmt_110[] = "(/,\002 Estimated Urey-Bradley Parameters :\002"
	    ",/)";
    static char fmt_120[] = "(\002 ureybrad\002,2x,3i5,f10.1,f11.4)";
    static char fmt_130[] = "(/,\002 Estimated Out-of-Plane Parameters :\002"
	    ",/)";
    static char fmt_140[] = "(\002 opbend\002,4x,4i5,6x,f10.2)";
    static char fmt_150[] = "(/,\002 Estimated Torsional Parameters :\002/)";
    static char fmt_160[] = "(\002 torsion\002,3x,4i5,3x,f8.3,\002 0.0 1\002"
	    ",f8.3,\002 180.0 2\002,f8.3,\002 0.0 3\002)";

    /* System generated locals */
    address a__1[2], a__2[3], a__3[4];
    integer i__1, i__2, i__3[2], i__4[3], i__5[4], i__6;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_cmp(char *, char *, ftnlen, 
	    ftnlen), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen),
	     s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal), acos(doublereal);

    /* Local variables */
    extern doublereal bndguess_(integer *, integer *), angguess_(integer *, 
	    integer *, integer *), opbguess_(integer *, integer *, integer *, 
	    integer *);
    extern /* Subroutine */ int vdwguess_(integer *, doublereal *, doublereal 
	    *, doublereal *), torguess_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *);
    extern doublereal uryguess_(integer *, integer *, integer *);
    static integer i__, j, k, ia, ib, ic, id, nb, na;
    static char pa[4], pb[4], pc[4];
    static integer nv, nt;
    static char pd[4];
    static integer iia, iib;
    static doublereal xab;
    static integer ita, itb, itc, itd, iva, ivb, ivc, nsb;
    static doublereal yab, zab, xcb;
    static integer nop;
    static doublereal ycb, zcb, xac, yac, zac, dot;
    static char ptb[8], pta[12], ptt[16];
    static doublereal rab2, rcb2;
    static integer isba, isbb, iita, iitb;
    static logical done;
    static integer size, vnum[5000];
    extern /* Subroutine */ int fatal_(void);
    static doublereal cosine;
    extern integer number_(char *, ftnlen);
    static integer nequiv[75000];
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    , sbguess_(integer *, integer *, integer *, doublereal *, 
	    doublereal *);

    /* Fortran I/O blocks */
    static cilist io___26 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___49 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___64 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___69 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___75 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___87 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___88 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___90 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___91 = { 0, 0, 0, fmt_160, 0 };



#define t1_ref(a_1,a_2) ktorsn_1.t1[(a_2)*2 + a_1 - 3]
#define t2_ref(a_1,a_2) ktorsn_1.t2[(a_2)*2 + a_1 - 3]
#define t3_ref(a_1,a_2) ktorsn_1.t3[(a_2)*2 + a_1 - 3]
#define ka_ref(a_0,a_1) &kangs_1.ka[(a_1)*12 + a_0 - 12]
#define kb_ref(a_0,a_1) &kbonds_1.kb[(a_1)*8 + a_0 - 8]
#define kt_ref(a_0,a_1) &ktorsn_1.kt[(a_1)*16 + a_0 - 16]
#define ku_ref(a_0,a_1) &kurybr_1.ku[(a_1)*12 + a_0 - 12]
#define ang_ref(a_1,a_2) kangs_1.ang[(a_2)*3 + a_1 - 4]
#define ksb_ref(a_0,a_1) &kstbnd_1.ksb[(a_1)*12 + a_0 - 12]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define kopb_ref(a_0,a_1) &kopbnd_1.kopb[(a_1)*16 + a_0 - 16]
#define stbn_ref(a_1,a_2) kstbnd_1.stbn[(a_2)*2 + a_1 - 3]
#define iury_ref(a_1,a_2) urey_1.iury[(a_2)*3 + a_1 - 4]
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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kangs.i  --  forcefield parameters for bond angle bending  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxna    maximum number of harmonic angle bend parameter entries */
/*     maxna5   maximum number of 5-membered ring angle bend entries */
/*     maxna4   maximum number of 4-membered ring angle bend entries */
/*     maxna3   maximum number of 3-membered ring angle bend entries */
/*     maxnaf   maximum number of Fourier angle bend parameter entries */

/*     acon     force constant parameters for harmonic angle bends */
/*     acon5    force constant parameters for 5-ring angle bends */
/*     acon4    force constant parameters for 4-ring angle bends */
/*     acon3    force constant parameters for 3-ring angle bends */
/*     aconf    force constant parameters for Fourier angle bends */
/*     ang      bond angle parameters for harmonic angle bends */
/*     ang5     bond angle parameters for 5-ring angle bends */
/*     ang4     bond angle parameters for 4-ring angle bends */
/*     ang3     bond angle parameters for 3-ring angle bends */
/*     angf     phase shift angle and periodicity for Fourier bends */
/*     ka       string of atom classes for harmonic angle bends */
/*     ka5      string of atom classes for 5-ring angle bends */
/*     ka4      string of atom classes for 4-ring angle bends */
/*     ka3      string of atom classes for 3-ring angle bends */
/*     kaf      string of atom classes for Fourier angle bends */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  kbonds.i  --  forcefield parameters for bond stretching  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxnb   maximum number of bond stretch parameter entries */
/*     maxnb5  maximum number of 5-membered ring bond stretch entries */
/*     maxnb4  maximum number of 4-membered ring bond stretch entries */
/*     maxnb3  maximum number of 3-membered ring bond stretch entries */
/*     maxnel  maximum number of electronegativity bond corrections */

/*     bcon    force constant parameters for harmonic bond stretch */
/*     blen    bond length parameters for harmonic bond stretch */
/*     bcon5   force constant parameters for 5-ring bond stretch */
/*     blen5   bond length parameters for 5-ring bond stretch */
/*     bcon4   force constant parameters for 4-ring bond stretch */
/*     blen4   bond length parameters for 4-ring bond stretch */
/*     bcon3   force constant parameters for 3-ring bond stretch */
/*     blen3   bond length parameters for 3-ring bond stretch */
/*     dlen    electronegativity bond length correction parameters */
/*     kb      string of atom classes for harmonic bond stretch */
/*     kb5     string of atom classes for 5-ring bond stretch */
/*     kb4     string of atom classes for 4-ring bond stretch */
/*     kb3     string of atom classes for 3-ring bond stretch */
/*     kel     string of atom classes for electronegativity corrections */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kopbnd.i  --  forcefield parameters for out-of-plane bend  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxnopb   maximum number of out-of-plane bending entries */

/*     opbn      force constant parameters for out-of-plane bending */
/*     kopb      string of atom classes for out-of-plane bending */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  kstbnd.i  --  forcefield parameters for stretch-bend  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxnsb   maximum number of stretch-bend parameter entries */

/*     stbn     force constant parameters for stretch-bend terms */
/*     ksb      string of atom classes for stretch-bend terms */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  ktorsn.i  --  forcefield parameters for torsional angles  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxnt    maximum number of torsional angle parameter entries */
/*     maxnt5   maximum number of 5-membered ring torsion entries */
/*     maxnt4   maximum number of 4-membered ring torsion entries */

/*     t1       torsional parameters for standard 1-fold rotation */
/*     t2       torsional parameters for standard 2-fold rotation */
/*     t3       torsional parameters for standard 3-fold rotation */
/*     t4       torsional parameters for standard 4-fold rotation */
/*     t5       torsional parameters for standard 5-fold rotation */
/*     t6       torsional parameters for standard 6-fold rotation */
/*     t15      torsional parameters for 1-fold rotation in 5-ring */
/*     t25      torsional parameters for 2-fold rotation in 5-ring */
/*     t35      torsional parameters for 3-fold rotation in 5-ring */
/*     t45      torsional parameters for 4-fold rotation in 5-ring */
/*     t55      torsional parameters for 5-fold rotation in 5-ring */
/*     t65      torsional parameters for 6-fold rotation in 5-ring */
/*     t14      torsional parameters for 1-fold rotation in 4-ring */
/*     t24      torsional parameters for 2-fold rotation in 4-ring */
/*     t34      torsional parameters for 3-fold rotation in 4-ring */
/*     t44      torsional parameters for 4-fold rotation in 4-ring */
/*     t54      torsional parameters for 5-fold rotation in 4-ring */
/*     t64      torsional parameters for 6-fold rotation in 4-ring */
/*     kt       string of atom classes for torsional angles */
/*     kt5      string of atom classes for 5-ring torsions */
/*     kt4      string of atom classes for 4-ring torsions */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kurybr.i  --  forcefield parameters for Urey-Bradley terms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     maxnu   maximum number of Urey-Bradley parameter entries */

/*     ucon    force constant parameters for Urey-Bradley terms */
/*     dst13   ideal 1-3 distance parameters for Urey-Bradley terms */
/*     ku      string of atom classes for Urey-Bradley terms */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kvdws.i  --  forcefield parameters for van der Waals terms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     rad      van der Waals radius parameter for each atom type */
/*     eps      van der Waals well depth parameter for each atom type */
/*     rad4     van der Waals radius parameter in 1-4 interactions */
/*     eps4     van der Waals well depth parameter in 1-4 interactions */
/*     reduct   van der Waals reduction factor for each atom type */




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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  opbend.i  --  out-of-plane bends in the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     opbk      force constant values for out-of-plane bending */
/*     nopbend   total number of out-of-plane bends in the system */
/*     iopb      bond angle numbers used in out-of-plane bending */




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2009 by Chuanjie Wu and Jay William Ponder  ## */
/*     ##                    All Rights Reserved                     ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  qmstuf.i  --  quantum data from Gaussian 03 calculation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     gx       x-coordinate of each atom in the QM data file */
/*     gy       y-coordinate of each atom in the QM data file */
/*     gz       z-coordinate of each atom in the QM data file */
/*     gforce   force components on each atom from QM data */
/*     gh       Hessian maxtrix elements from QM data */
/*     gfreq    calculated vibrational frequencies from QM data */
/*     ngatom   number of atoms in the QM data file */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  strbnd.i  --  stretch-bends in the current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     sbk       force constants for stretch-bend terms */
/*     nstrbnd   total number of stretch-bend interactions */
/*     isb       angle and bond numbers used in stretch-bend */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  urey.i  --  Urey-Bradley interactions in the structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     uk      Urey-Bradley force constants (kcal/mole/Ang**2) */
/*     ul      ideal 1-3 distance values in Angstroms */
/*     nurey   total number of Urey-Bradley terms in the system */
/*     iury    numbers of the atoms in each Urey-Bradley interaction */




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  valfit.i  --  values for valence term parameter fitting  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     fit_bond    logical flag to fit bond stretch parameters */
/*     fit_angle   logical flag to fit angle bend parameters */
/*     fit_strbnd  logical flag to fit stretch-bend parameters */
/*     fit_urey    logical flag to fit Urey-Bradley parameters */
/*     fit_opbend  logical flag to fit out-of-plane bend parameters */
/*     fit_tors    logical flag to fit torsional parameters */
/*     fit_struct  logical flag to structure-fit valence parameters */
/*     fit_force   logical flag to force-fit valence parameters */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  vdwpot.i  --  specifics of van der Waals functional form  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     abuck       value of "A" constant in Buckingham vdw potential */
/*     bbuck       value of "B" constant in Buckingham vdw potential */
/*     cbuck       value of "C" constant in Buckingham vdw potential */
/*     ghal        value of "gamma" in buffered 14-7 vdw potential */
/*     dhal        value of "delta" in buffered 14-7 vdw potential */
/*     v2scale     factor by which 1-2 vdw interactions are scaled */
/*     v3scale     factor by which 1-3 vdw interactions are scaled */
/*     v4scale     factor by which 1-4 vdw interactions are scaled */
/*     v5scale     factor by which 1-5 vdw interactions are scaled */
/*     igauss      coefficients of Gaussian fit to vdw potential */
/*     ngauss      number of Gaussians used in fit to vdw potential */
/*     use_vcorr   flag to use long range vdw der Waals correction */
/*     vdwindex    indexing mode (atom type or class) for vdw parameters */
/*     vdwtyp      type of van der Waals potential energy function */
/*     radtyp      type of parameter (sigma or R-min) for atomic size */
/*     radsiz      atomic size provided as radius or diameter */
/*     radrule     combining rule for atomic size parameters */
/*     epsrule     combining rule for vdw well depth parameters */
/*     gausstyp    type of Gaussian fit to van der Waals potential */




/*     check the number of atoms in QM output and TINKER xyz file */

    if (atoms_1.n != qmstuf_1.ngatom) {
	io___26.ciunit = iounit_1.iout;
	s_wsfe(&io___26);
	e_wsfe();
	fatal_();
    }

/*     assign initial values to van der Waals parameters */

    nv = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ita = atmtyp_1.class__[i__ - 1];
	if (s_cmp(vdwpot_1.vdwindex, "TYPE", (ftnlen)5, (ftnlen)4) == 0) {
	    ita = atoms_1.type__[i__ - 1];
	}
	done = FALSE_;
	if (i__ > 1) {
	    i__2 = nv;
	    for (j = 1; j <= i__2; ++j) {
		if (ita == vnum[j - 1]) {
		    done = TRUE_;
		}
	    }
	}
	if (! done) {
	    ++nv;
	    vnum[nv - 1] = ita;
	    vdwguess_(&i__, &kvdws_1.rad[ita - 1], &kvdws_1.eps[ita - 1], &
		    kvdws_1.reduct[ita - 1]);
	}
    }

/*     print the initial van der Waals parameter values */

    if (nv > 0) {
	io___33.ciunit = iounit_1.iout;
	s_wsfe(&io___33);
	e_wsfe();
    }
    i__1 = nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = vnum[i__ - 1];
	if (kvdws_1.reduct[ia - 1] == 0.) {
	    io___35.ciunit = iounit_1.iout;
	    s_wsfe(&io___35);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kvdws_1.rad[ia - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kvdws_1.eps[ia - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	} else {
	    io___36.ciunit = iounit_1.iout;
	    s_wsfe(&io___36);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kvdws_1.rad[ia - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kvdws_1.eps[ia - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kvdws_1.reduct[ia - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }

/*     find and store the unique bond stretches in the system */

    nb = 0;
    i__1 = bond_1.nbond;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ibnd_ref(1, i__);
	ib = ibnd_ref(2, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	if (ita <= itb) {
/* Writing concatenation */
	    i__3[0] = 4, a__1[0] = pa;
	    i__3[1] = 4, a__1[1] = pb;
	    s_cat(ptb, a__1, i__3, &c__2, (ftnlen)8);
	} else {
/* Writing concatenation */
	    i__3[0] = 4, a__1[0] = pb;
	    i__3[1] = 4, a__1[1] = pa;
	    s_cat(ptb, a__1, i__3, &c__2, (ftnlen)8);
	}
	done = FALSE_;
	i__2 = nb;
	for (j = 1; j <= i__2; ++j) {
	    if (s_cmp(ptb, kb_ref(0, j), (ftnlen)8, (ftnlen)8) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    ++nb;
	    s_copy(kb_ref(0, nb), ptb, (ftnlen)8, (ftnlen)8);
	}
    }

/*     assign initial values to bond stretch parameters */

    k = 0;
    i__1 = bond_1.nbond;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ibnd_ref(1, i__);
	ib = ibnd_ref(2, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	if (ita <= itb) {
/* Writing concatenation */
	    i__3[0] = 4, a__1[0] = pa;
	    i__3[1] = 4, a__1[1] = pb;
	    s_cat(ptb, a__1, i__3, &c__2, (ftnlen)8);
	} else {
/* Writing concatenation */
	    i__3[0] = 4, a__1[0] = pb;
	    i__3[1] = 4, a__1[1] = pa;
	    s_cat(ptb, a__1, i__3, &c__2, (ftnlen)8);
	}
	xab = qmstuf_1.gx[ia - 1] - qmstuf_1.gx[ib - 1];
	yab = qmstuf_1.gy[ia - 1] - qmstuf_1.gy[ib - 1];
	zab = qmstuf_1.gz[ia - 1] - qmstuf_1.gz[ib - 1];
	bond_1.bl[i__ - 1] = sqrt(xab * xab + yab * yab + zab * zab);
	done = FALSE_;
	i__2 = k;
	for (j = 1; j <= i__2; ++j) {
	    if (s_cmp(ptb, kb_ref(0, j), (ftnlen)8, (ftnlen)8) == 0) {
		done = TRUE_;
		kbonds_1.blen[j - 1] += bond_1.bl[i__ - 1];
		++nequiv[j - 1];
	    }
	}
	if (! done) {
	    ++k;
	    kbonds_1.bcon[k - 1] = bndguess_(&ia, &ib);
	    kbonds_1.blen[k - 1] = bond_1.bl[i__ - 1];
	    nequiv[k - 1] = 1;
	}
    }

/*     print the initial bond stretch parameter values */

    if (nb > 0) {
	io___49.ciunit = iounit_1.iout;
	s_wsfe(&io___49);
	e_wsfe();
    }
    i__1 = nb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	kbonds_1.blen[i__ - 1] /= (doublereal) nequiv[i__ - 1];
	s_copy(ptb, kb_ref(0, i__), (ftnlen)8, (ftnlen)8);
	ia = number_(ptb, (ftnlen)4);
	ib = number_(ptb + 4, (ftnlen)4);
	io___50.ciunit = iounit_1.iout;
	s_wsfe(&io___50);
	do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&kbonds_1.bcon[i__ - 1], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&kbonds_1.blen[i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_wsfe();
    }

/*     find and store the unique angle bends in the system */

    na = 0;
    i__1 = angle_1.nangle;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pa;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pc;
	    s_cat(pta, a__2, i__4, &c__3, (ftnlen)12);
	} else {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pc;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pa;
	    s_cat(pta, a__2, i__4, &c__3, (ftnlen)12);
	}
	done = FALSE_;
	i__2 = na;
	for (j = 1; j <= i__2; ++j) {
	    if (s_cmp(pta, ka_ref(0, j), (ftnlen)12, (ftnlen)12) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    ++na;
	    s_copy(ka_ref(0, na), pta, (ftnlen)12, (ftnlen)12);
	}
    }

/*     assign initial values to angle bend parameters */

    k = 0;
    i__1 = angle_1.nangle;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pa;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pc;
	    s_cat(pta, a__2, i__4, &c__3, (ftnlen)12);
	} else {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pc;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pa;
	    s_cat(pta, a__2, i__4, &c__3, (ftnlen)12);
	}
	xab = qmstuf_1.gx[ia - 1] - qmstuf_1.gx[ib - 1];
	yab = qmstuf_1.gy[ia - 1] - qmstuf_1.gy[ib - 1];
	zab = qmstuf_1.gz[ia - 1] - qmstuf_1.gz[ib - 1];
	xcb = qmstuf_1.gx[ic - 1] - qmstuf_1.gx[ib - 1];
	ycb = qmstuf_1.gy[ic - 1] - qmstuf_1.gy[ib - 1];
	zcb = qmstuf_1.gz[ic - 1] - qmstuf_1.gz[ib - 1];
	rab2 = xab * xab + yab * yab + zab * zab;
	rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
	if (rab2 != 0. && rcb2 != 0.) {
	    dot = xab * xcb + yab * ycb + zab * zcb;
	    cosine = dot / sqrt(rab2 * rcb2);
/* Computing MIN */
	    d__1 = 1., d__2 = max(-1.,cosine);
	    cosine = min(d__1,d__2);
	    angle_1.anat[i__ - 1] = acos(cosine) * 57.29577951308232088;
	}
	done = FALSE_;
	i__2 = k;
	for (j = 1; j <= i__2; ++j) {
	    if (s_cmp(pta, ka_ref(0, j), (ftnlen)12, (ftnlen)12) == 0) {
		done = TRUE_;
		ang_ref(1, j) = ang_ref(1, j) + angle_1.anat[i__ - 1];
		++nequiv[j - 1];
	    }
	}
	if (! done) {
	    ++k;
	    kangs_1.acon[k - 1] = angguess_(&ia, &ib, &ic);
	    ang_ref(1, k) = angle_1.anat[i__ - 1];
	    nequiv[k - 1] = 1;
	}
    }

/*     print the initial angle bend parameter values */

    if (na > 0) {
	io___63.ciunit = iounit_1.iout;
	s_wsfe(&io___63);
	e_wsfe();
    }
    i__1 = na;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ang_ref(1, i__) = ang_ref(1, i__) / (doublereal) nequiv[i__ - 1];
	s_copy(pta, ka_ref(0, i__), (ftnlen)12, (ftnlen)12);
	ia = number_(pta, (ftnlen)4);
	ib = number_(pta + 4, (ftnlen)4);
	ic = number_(pta + 8, (ftnlen)4);
	io___64.ciunit = iounit_1.iout;
	s_wsfe(&io___64);
	do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&kangs_1.acon[i__ - 1], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&ang_ref(1, i__), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     assign initial values to stretch-bend parameters */

    nsb = 0;
    i__1 = angle_1.nangle;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	iva = atmtyp_1.valence[ia - 1];
	ivb = atmtyp_1.valence[ib - 1];
	ivc = atmtyp_1.valence[ic - 1];
	if (iva > 1 || ivc > 1) {
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    numeral_(&itc, pc, &size, (ftnlen)4);
	    if (ita <= itc) {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pa;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pc;
		s_cat(pta, a__2, i__4, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pc;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pa;
		s_cat(pta, a__2, i__4, &c__3, (ftnlen)12);
	    }
	    done = FALSE_;
	    i__2 = nsb;
	    for (j = 1; j <= i__2; ++j) {
		if (s_cmp(pta, ksb_ref(0, j), (ftnlen)12, (ftnlen)12) == 0) {
		    done = TRUE_;
		}
	    }
	    if (! done) {
		++nsb;
		s_copy(ksb_ref(0, nsb), pta, (ftnlen)12, (ftnlen)12);
		if (ita <= itc) {
		    sbguess_(&ia, &ib, &ic, &stbn_ref(1, nsb), &stbn_ref(2, 
			    nsb));
		} else {
		    sbguess_(&ic, &ib, &ia, &stbn_ref(1, nsb), &stbn_ref(2, 
			    nsb));
		}
	    }
	}
    }

/*     print the initial stretch-bend parameter values */

    if (nsb > 0) {
	io___69.ciunit = iounit_1.iout;
	s_wsfe(&io___69);
	e_wsfe();
    }
    i__1 = nsb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(pta, ksb_ref(0, i__), (ftnlen)12, (ftnlen)12);
	ia = number_(pta, (ftnlen)4);
	ib = number_(pta + 4, (ftnlen)4);
	ic = number_(pta + 8, (ftnlen)4);
	io___70.ciunit = iounit_1.iout;
	s_wsfe(&io___70);
	do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&stbn_ref(1, i__), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&stbn_ref(2, i__), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     assign initial values to Urey-Bradley parameters */

    k = 0;
    i__1 = urey_1.nurey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iury_ref(1, i__);
	ib = iury_ref(2, i__);
	ic = iury_ref(3, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pa;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pc;
	    s_cat(pta, a__2, i__4, &c__3, (ftnlen)12);
	} else {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pc;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pa;
	    s_cat(pta, a__2, i__4, &c__3, (ftnlen)12);
	}
	xac = qmstuf_1.gx[ia - 1] - qmstuf_1.gx[ic - 1];
	yac = qmstuf_1.gy[ia - 1] - qmstuf_1.gy[ic - 1];
	zac = qmstuf_1.gz[ia - 1] - qmstuf_1.gz[ic - 1];
	urey_1.ul[i__ - 1] = sqrt(xac * xac + yac * yac + zac * zac);
	done = FALSE_;
	i__2 = k;
	for (j = 1; j <= i__2; ++j) {
	    if (s_cmp(pta, ku_ref(0, j), (ftnlen)12, (ftnlen)12) == 0) {
		done = TRUE_;
		kurybr_1.dst13[j - 1] += urey_1.ul[i__ - 1];
		++nequiv[j - 1];
	    }
	}
	if (! done) {
	    ++k;
	    kurybr_1.ucon[k - 1] = uryguess_(&ia, &ib, &ic);
	    kurybr_1.dst13[k - 1] = urey_1.ul[i__ - 1];
	    nequiv[k - 1] = 1;
	}
    }

/*     print the initial Urey-Bradley parameter values */

    if (urey_1.nurey > 0) {
	io___74.ciunit = iounit_1.iout;
	s_wsfe(&io___74);
	e_wsfe();
    }
    i__1 = nsb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(pta, ku_ref(0, i__), (ftnlen)12, (ftnlen)12);
	ia = number_(pta, (ftnlen)4);
	ib = number_(pta + 4, (ftnlen)4);
	ic = number_(pta + 8, (ftnlen)4);
	io___75.ciunit = iounit_1.iout;
	s_wsfe(&io___75);
	do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&kurybr_1.ucon[i__ - 1], (ftnlen)sizeof(
		doublereal));
	do_fio(&c__1, (char *)&kurybr_1.dst13[i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_wsfe();
    }

/*     assign initial values to out-of-plane bend parameters */

    nop = 0;
    i__1 = bond_1.nbond;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ibnd_ref(1, i__);
	ib = ibnd_ref(2, i__);
	ic = 0;
	id = 0;
	iva = atmtyp_1.valence[ia - 1];
	ivb = atmtyp_1.valence[ib - 1];
	if (iva == 3 || ivb == 3) {
	    ita = atmtyp_1.class__[ia - 1];
	    itb = atmtyp_1.class__[ib - 1];
	    itc = 0;
	    itd = 0;
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    numeral_(&itc, pc, &size, (ftnlen)4);
	    numeral_(&itd, pd, &size, (ftnlen)4);
	    if (iva == 3) {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pb;
		i__5[1] = 4, a__3[1] = pa;
		i__5[2] = 4, a__3[2] = pc;
		i__5[3] = 4, a__3[3] = pd;
		s_cat(ptt, a__3, i__5, &c__4, (ftnlen)16);
		isba = ia;
		isbb = ib;
	    } else {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pa;
		i__5[1] = 4, a__3[1] = pb;
		i__5[2] = 4, a__3[2] = pc;
		i__5[3] = 4, a__3[3] = pd;
		s_cat(ptt, a__3, i__5, &c__4, (ftnlen)16);
		isba = ib;
		isbb = ia;
	    }
	    if (atmtyp_1.atomic[isba - 1] == 6) {
		done = FALSE_;
		i__2 = nop;
		for (j = 1; j <= i__2; ++j) {
		    if (s_cmp(ptt, kopb_ref(0, j), (ftnlen)16, (ftnlen)16) == 
			    0) {
			done = TRUE_;
		    }
		}
		if (! done) {
		    ++nop;
		    s_copy(kopb_ref(0, nop), ptt, (ftnlen)16, (ftnlen)16);
		    kopbnd_1.opbn[nop - 1] = opbguess_(&isba, &isbb, &ic, &id)
			    ;
		    i__2 = bond_1.nbond;
		    for (j = i__ + 1; j <= i__2; ++j) {
			iia = ibnd_ref(1, j);
			iib = ibnd_ref(2, j);
			if (iia == isba || iib == isba) {
			    iita = atmtyp_1.class__[iia - 1];
			    iitb = atmtyp_1.class__[iib - 1];
			    size = 4;
			    numeral_(&iita, pa, &size, (ftnlen)4);
			    numeral_(&iitb, pb, &size, (ftnlen)4);
			    if (iia == isba) {
/* Writing concatenation */
				i__5[0] = 4, a__3[0] = pb;
				i__5[1] = 4, a__3[1] = pa;
				i__5[2] = 4, a__3[2] = pc;
				i__5[3] = 4, a__3[3] = pd;
				s_cat(ptt, a__3, i__5, &c__4, (ftnlen)16);
			    } else if (iib == isba) {
/* Writing concatenation */
				i__5[0] = 4, a__3[0] = pa;
				i__5[1] = 4, a__3[1] = pb;
				i__5[2] = 4, a__3[2] = pc;
				i__5[3] = 4, a__3[3] = pd;
				s_cat(ptt, a__3, i__5, &c__4, (ftnlen)16);
			    }
			    done = FALSE_;
			    i__6 = nop;
			    for (k = 1; k <= i__6; ++k) {
				if (s_cmp(ptt, kopb_ref(0, k), (ftnlen)16, (
					ftnlen)16) == 0) {
				    done = TRUE_;
				}
			    }
			    if (! done) {
				++nop;
				s_copy(kopb_ref(0, nop), ptt, (ftnlen)16, (
					ftnlen)16);
				if (iia == isba) {
				    kopbnd_1.opbn[nop - 1] = opbguess_(&iia, &
					    iib, &ic, &id);
				} else if (iib == isba) {
				    kopbnd_1.opbn[nop - 1] = opbguess_(&iib, &
					    iia, &ic, &id);
				}
			    }
			}
		    }
		}
	    } else if (atmtyp_1.atomic[isba - 1] == 7) {
		if (atmtyp_1.valence[isbb - 1] == 3 && atmtyp_1.atomic[isbb - 
			1] == 6) {
		    ++nop;
		    s_copy(kopb_ref(0, nop), ptt, (ftnlen)16, (ftnlen)16);
		    kopbnd_1.opbn[nop - 1] = opbguess_(&isba, &isbb, &ic, &id)
			    ;
		    i__2 = bond_1.nbond;
		    for (j = 1; j <= i__2; ++j) {
			if (j != i__ && (ibnd_ref(1, j) == isba || ibnd_ref(2,
				 j) == isba)) {
			    if (ibnd_ref(1, j) == isba) {
				iia = ibnd_ref(2, j);
			    } else {
				iia = ibnd_ref(1, j);
			    }
			    size = 4;
			    numeral_(&atmtyp_1.class__[isba - 1], pa, &size, (
				    ftnlen)4);
			    numeral_(&atmtyp_1.class__[iia - 1], pb, &size, (
				    ftnlen)4);
/* Writing concatenation */
			    i__5[0] = 4, a__3[0] = pb;
			    i__5[1] = 4, a__3[1] = pa;
			    i__5[2] = 4, a__3[2] = pc;
			    i__5[3] = 4, a__3[3] = pd;
			    s_cat(ptt, a__3, i__5, &c__4, (ftnlen)16);
			    done = FALSE_;
			    i__6 = nop;
			    for (k = 1; k <= i__6; ++k) {
				if (s_cmp(ptt, ksb_ref(0, k), (ftnlen)16, (
					ftnlen)12) == 0) {
				    done = TRUE_;
				}
			    }
			    if (! done) {
				++nop;
				s_copy(kopb_ref(0, nop), ptt, (ftnlen)16, (
					ftnlen)16);
				kopbnd_1.opbn[nop - 1] = opbguess_(&isba, &
					iia, &ic, &id);
			    }
			}
		    }
		}
	    }
	}
    }

/*     print the initial out-of-plane bend parameter values */

    if (nop > 0) {
	io___87.ciunit = iounit_1.iout;
	s_wsfe(&io___87);
	e_wsfe();
    }
    i__1 = nop;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(ptt, kopb_ref(0, i__), (ftnlen)16, (ftnlen)16);
	ia = number_(ptt, (ftnlen)4);
	ib = number_(ptt + 4, (ftnlen)4);
	ic = number_(ptt + 8, (ftnlen)4);
	id = number_(ptt + 12, (ftnlen)4);
	io___88.ciunit = iounit_1.iout;
	s_wsfe(&io___88);
	do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&kopbnd_1.opbn[i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_wsfe();
    }

/*     assign initial values to torsional parameters */

    nt = 0;
    i__1 = tors_1.ntors;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = itors_ref(1, i__);
	ib = itors_ref(2, i__);
	ic = itors_ref(3, i__);
	id = itors_ref(4, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	itd = atmtyp_1.class__[id - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	numeral_(&itd, pd, &size, (ftnlen)4);
	if (itb <= itc) {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pa;
	    i__5[1] = 4, a__3[1] = pb;
	    i__5[2] = 4, a__3[2] = pc;
	    i__5[3] = 4, a__3[3] = pd;
	    s_cat(ptt, a__3, i__5, &c__4, (ftnlen)16);
	} else {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pd;
	    i__5[1] = 4, a__3[1] = pc;
	    i__5[2] = 4, a__3[2] = pb;
	    i__5[3] = 4, a__3[3] = pa;
	    s_cat(ptt, a__3, i__5, &c__4, (ftnlen)16);
	}
	done = FALSE_;
	i__2 = nt;
	for (j = 1; j <= i__2; ++j) {
	    if (s_cmp(ptt, kt_ref(0, j), (ftnlen)16, (ftnlen)16) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    ++nt;
	    s_copy(kt_ref(0, nt), ptt, (ftnlen)16, (ftnlen)16);
	    torguess_(&ia, &ib, &ic, &id, &t1_ref(1, nt), &t2_ref(1, nt), &
		    t3_ref(1, nt));
	}
    }

/*     print the initial torsional parameter values */

    if (nt > 0) {
	io___90.ciunit = iounit_1.iout;
	s_wsfe(&io___90);
	e_wsfe();
    }
    i__1 = nt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(ptt, kt_ref(0, i__), (ftnlen)16, (ftnlen)16);
	ia = number_(ptt, (ftnlen)4);
	ib = number_(ptt + 4, (ftnlen)4);
	ic = number_(ptt + 8, (ftnlen)4);
	id = number_(ptt + 12, (ftnlen)4);
	io___91.ciunit = iounit_1.iout;
	s_wsfe(&io___91);
	do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&t1_ref(1, i__), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&t2_ref(1, i__), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&t3_ref(1, i__), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* valguess_ */

#undef itors_ref
#undef iury_ref
#undef stbn_ref
#undef kopb_ref
#undef iang_ref
#undef ibnd_ref
#undef ksb_ref
#undef ang_ref
#undef ku_ref
#undef kt_ref
#undef kb_ref
#undef ka_ref
#undef t3_ref
#undef t2_ref
#undef t1_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine vdwguess  --  estimate van der Waals parameters  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "vdwguess" sets initial VDW parameters based on atom type */
/*     and connected atoms */


/* Subroutine */ int vdwguess_(integer *ia, doublereal *rad, doublereal *eps, 
	doublereal *reduce)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k, ita, itb, iva, ivb;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]



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
/*     ##  vdwpot.i  --  specifics of van der Waals functional form  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     abuck       value of "A" constant in Buckingham vdw potential */
/*     bbuck       value of "B" constant in Buckingham vdw potential */
/*     cbuck       value of "C" constant in Buckingham vdw potential */
/*     ghal        value of "gamma" in buffered 14-7 vdw potential */
/*     dhal        value of "delta" in buffered 14-7 vdw potential */
/*     v2scale     factor by which 1-2 vdw interactions are scaled */
/*     v3scale     factor by which 1-3 vdw interactions are scaled */
/*     v4scale     factor by which 1-4 vdw interactions are scaled */
/*     v5scale     factor by which 1-5 vdw interactions are scaled */
/*     igauss      coefficients of Gaussian fit to vdw potential */
/*     ngauss      number of Gaussians used in fit to vdw potential */
/*     use_vcorr   flag to use long range vdw der Waals correction */
/*     vdwindex    indexing mode (atom type or class) for vdw parameters */
/*     vdwtyp      type of van der Waals potential energy function */
/*     radtyp      type of parameter (sigma or R-min) for atomic size */
/*     radsiz      atomic size provided as radius or diameter */
/*     radrule     combining rule for atomic size parameters */
/*     epsrule     combining rule for vdw well depth parameters */
/*     gausstyp    type of Gaussian fit to van der Waals potential */




/*     set default value for radius, well depth and reduction factor */

    *rad = 1.;
    *eps = .1;
    *reduce = 0.;

/*     get atomic number and valence for the atom and its neighbor */

    ita = atmtyp_1.atomic[*ia - 1];
    iva = atmtyp_1.valence[*ia - 1];
    itb = 0;
    ivb = 0;
    i__1 = couple_1.n12[*ia - 1];
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = i12_ref(i__, *ia);
	k = atmtyp_1.atomic[j - 1];
	if (k > itb) {
	    itb = k;
	    ivb = atmtyp_1.valence[j - 1];
	}
    }

/*     assign specific values based on atom type and connectivity */

    if (ita == 1) {
	if (itb == 6) {
	    if (ivb == 3) {
		*rad = 2.98;
		*eps = .026;
		*reduce = .92;
	    } else if (ivb == 4) {
		*rad = 2.78;
		*eps = .026;
		*reduce = .91;
	    } else {
		*rad = 2.78;
		*eps = .026;
		*reduce = .91;
	    }
	} else if (itb == 7) {
	    *rad = 2.7;
	    *eps = .02;
	    *reduce = .91;
	} else if (itb == 8) {
	    *rad = 2.655;
	    *eps = .0135;
	    *reduce = .91;
	} else if (itb == 16) {
	    *rad = 3.;
	    *eps = .0265;
	    *reduce = .98;
	} else {
	    *rad = 2.98;
	    *eps = .026;
	    *reduce = .92;
	}
    } else if (ita == 6) {
	if (iva == 3) {
	    *rad = 3.8;
	    *eps = .089;
	} else if (iva == 4) {
	    *rad = 3.82;
	    *eps = .101;
	} else {
	    *rad = 3.82;
	    *eps = .101;
	}
    } else if (ita == 7) {
	if (iva == 3) {
	    *rad = 3.71;
	    *eps = .105;
	} else if (iva == 2) {
	    *rad = 3.71;
	    *eps = .11;
	} else {
	    *rad = 3.71;
	    *eps = .105;
	}
    } else if (ita == 8) {
	if (iva == 1) {
	    if (itb == 6) {
		*rad = 3.3;
		*eps = .112;
	    } else if (itb == 7) {
		*rad = 3.3;
		*eps = .112;
	    } else if (itb == 15) {
		*rad = 3.36;
		*eps = .112;
	    } else if (itb == 16) {
		*rad = 3.51;
		*eps = .112;
	    } else {
		*rad = 3.3;
		*eps = .112;
	    }
	} else if (iva == 2) {
	    if (itb == 15) {
		*rad = 3.405;
		*eps = .112;
	    } else {
		*rad = 3.405;
		*eps = .11;
	    }
	} else {
	    *rad = 3.405;
	    *eps = .11;
	}
    } else if (ita == 9) {
	if (iva == 0) {
	    *rad = 3.4;
	    *eps = .25;
	} else if (iva == 1) {
	    *rad = 3.22;
	    *eps = .12;
	} else {
	    *rad = 3.22;
	    *eps = .12;
	}
    } else if (ita == 11) {
	*rad = 3.02;
	*eps = .26;
    } else if (ita == 12) {
	*rad = 2.55;
	*eps = .85;
    } else if (ita == 15) {
	*rad = 4.45;
	*eps = .39;
    } else if (ita == 16) {
	if (iva == 2) {
	    *rad = 3.91;
	    *eps = .385;
	} else if (iva == 3) {
	    *rad = 3.91;
	    *eps = .385;
	} else if (iva == 4) {
	    *rad = 3.91;
	    *eps = .385;
	} else {
	    *rad = 3.91;
	    *eps = .385;
	}
    } else if (ita == 17) {
	if (iva == 0) {
	    *rad = 4.13;
	    *eps = .34;
	} else if (iva == 1) {
	    *rad = 4.13;
	    *eps = .34;
	} else {
	    *rad = 4.13;
	    *eps = .34;
	}
    } else if (ita == 19) {
	*rad = 3.71;
	*eps = .35;
    } else if (ita == 20) {
	*rad = 3.15;
	*eps = 1.6;
    } else if (ita == 35) {
	if (iva == 0) {
	    *rad = 4.38;
	    *eps = .43;
	} else if (iva == 1) {
	    *rad = 4.38;
	    *eps = .43;
	} else {
	    *rad = 4.38;
	    *eps = .43;
	}
    } else if (ita == 53) {
	if (iva == 0) {
	    *rad = 4.66;
	    *eps = .52;
	} else if (iva == 1) {
	    *rad = 4.66;
	    *eps = .52;
	} else {
	    *rad = 4.66;
	    *eps = .52;
	}
    }

/*     scale the vdw parameters to the desired units */

    if (s_cmp(vdwpot_1.radsiz, "RADIUS", (ftnlen)8, (ftnlen)6) == 0) {
	*rad *= .5;
    }
    if (s_cmp(vdwpot_1.radsiz, "SIGMA", (ftnlen)8, (ftnlen)5) == 0) {
	*rad /= 1.122462048309372981;
    }
    return 0;
} /* vdwguess_ */

#undef i12_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  function bndguess  --  estimate bond stretch parameters  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "bndguess" sets approximate bond strech force constants based */
/*     on atom type and connected atoms */


doublereal bndguess_(integer *ia, integer *ib)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer ita, itb, iva, ivb, tmp;



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  bndpot.i  --  specifics of bond stretch functional forms  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     cbnd      cubic coefficient in bond stretch potential */
/*     qbnd      quartic coefficient in bond stretch potential */
/*     bndunit   convert bond stretch energy to kcal/mole */
/*     bndtyp    type of bond stretch potential energy function */




/*     get the atomic number and valence of each atom */

    ita = atmtyp_1.atomic[*ia - 1];
    itb = atmtyp_1.atomic[*ib - 1];
    iva = atmtyp_1.valence[*ia - 1];
    ivb = atmtyp_1.valence[*ib - 1];

/*     reverse the atom order based on atomic number */

    if (ita > itb) {
	tmp = ita;
	ita = itb;
	itb = tmp;
	tmp = iva;
	iva = ivb;
	ivb = tmp;
    }

/*     assign estimated bond stretch force constants */

    if (ita == 1) {
	if (itb == 6) {
	    if (ivb == 3) {
		ret_val = 410.;
	    } else if (ivb == 4) {
		ret_val = 400.;
	    } else {
		ret_val = 400.;
	    }
	} else if (itb == 7) {
	    ret_val = 520.;
	} else if (itb == 8) {
	    ret_val = 560.;
	} else if (itb == 9) {
	    ret_val = 500.;
	} else if (itb == 14) {
	    ret_val = 200.;
	} else if (itb == 15) {
	    ret_val = 230.;
	} else if (itb == 16) {
	    ret_val = 260.;
	} else {
	    ret_val = 300.;
	}
    } else if (ita == 6) {
	if (itb == 6) {
	    if (iva == 3 && ivb == 3) {
		ret_val = 680.;
	    } else if (iva == 4 || ivb == 4) {
		ret_val = 385.;
	    } else {
		ret_val = 350.;
	    }
	} else if (itb == 7) {
	    if (iva == 3 && ivb == 2) {
		ret_val = 435.;
	    } else if (iva == 3 && ivb == 3) {
		ret_val = 250.;
	    } else if (iva == 4) {
		ret_val = 400.;
	    } else {
		ret_val = 450.;
	    }
	} else if (itb == 8) {
	    if (ivb == 1) {
		ret_val = 680.;
	    } else if (ivb == 2) {
		ret_val = 465.;
	    } else {
		ret_val = 465.;
	    }
	} else if (itb == 9) {
	    ret_val = 350.;
	} else if (itb == 14) {
	    ret_val = 350.;
	} else if (itb == 15) {
	    ret_val = 350.;
	} else if (itb == 16) {
	    ret_val = 216.;
	} else if (itb == 17) {
	    ret_val = 350.;
	} else {
	    ret_val = 450.;
	}
    } else if (ita == 7) {
	if (itb == 7) {
	    if (iva == 1) {
		ret_val = 1613.;
	    } else if (iva == 2 && ivb == 2) {
		ret_val = 950.;
	    } else {
		ret_val = 850.;
	    }
	} else if (itb == 8) {
	    if (ivb == 1) {
		ret_val = 900.;
	    } else {
		ret_val = 750.;
	    }
	} else if (itb == 14) {
	    ret_val = 450.;
	} else if (itb == 15) {
	    ret_val = 500.;
	} else if (itb == 16) {
	    ret_val = 550.;
	} else {
	    ret_val = 600.;
	}
    } else if (ita == 8) {
	if (itb == 8) {
	    ret_val = 750.;
	} else if (itb == 14) {
	    ret_val = 500.;
	} else if (itb == 15) {
	    if (iva == 2) {
		ret_val = 450.;
	    } else if (iva == 1) {
		ret_val = 775.;
	    } else {
		ret_val = 450.;
	    }
	} else if (itb == 16) {
	    ret_val = 606.;
	} else if (itb == 17) {
	    ret_val = 500.;
	} else {
	    ret_val = 600.;
	}
    } else if (ita == 14) {
	if (itb == 14) {
	    ret_val = 400.;
	} else if (itb == 15) {
	    ret_val = 450.;
	} else if (itb == 16) {
	    ret_val = 500.;
	} else if (itb == 17) {
	    ret_val = 650.;
	} else {
	    ret_val = 450.;
	}
    } else if (ita == 16) {
	if (itb == 16) {
	    ret_val = 188.;
	} else {
	    ret_val = 250.;
	}
    } else if (ita == 17) {
	ret_val = 300.;
    } else {
	ret_val = 350.;
    }

/*     scale the force constant to the desired units */

    ret_val /= bndpot_1.bndunit;
    return ret_val;
} /* bndguess_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function angguess  --  estimate angle bending parameters  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "angguess" sets approximate angle bend force constants based */
/*     on atom type and connected atoms */


doublereal angguess_(integer *ia, integer *ib, integer *ic)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer ita, itb, itc, iva, ivb, ivc, tmp;



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  angpot.i  --  specifics of angle bend functional forms  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     angunit    convert angle bending energy to kcal/mole */
/*     stbnunit   convert stretch-bend energy to kcal/mole */
/*     aaunit     convert angle-angle energy to kcal/mole */
/*     opbunit    convert out-of-plane bend energy to kcal/mole */
/*     opdunit    convert out-of-plane distance energy to kcal/mole */
/*     cang       cubic coefficient in angle bending potential */
/*     qang       quartic coefficient in angle bending potential */
/*     pang       quintic coefficient in angle bending potential */
/*     sang       sextic coefficient in angle bending potential */
/*     copb       cubic coefficient in out-of-plane bend potential */
/*     qopb       quartic coefficient in out-of-plane bend potential */
/*     popb       quintic coefficient in out-of-plane bend potential */
/*     sopb       sextic coefficient in out-of-plane bend potential */
/*     copd       cubic coefficient in out-of-plane distance potential */
/*     qopd       quartic coefficient in out-of-plane distance potential */
/*     popd       quintic coefficient in out-of-plane distance potential */
/*     sopd       sextic coefficient in out-of-plane distance potential */
/*     angtyp     type of angle bending function for each bond angle */
/*     opbtyp     type of out-of-plane bend potential energy function */




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




/*     get the atomic number and valence of each atom */

    ita = atmtyp_1.atomic[*ia - 1];
    itb = atmtyp_1.atomic[*ib - 1];
    itc = atmtyp_1.atomic[*ic - 1];
    iva = atmtyp_1.valence[*ia - 1];
    ivb = atmtyp_1.valence[*ib - 1];
    ivc = atmtyp_1.valence[*ic - 1];

/*     resort ja,jb,jc based on the atomic orders */

    if (ita > itc) {
	tmp = ita;
	ita = itc;
	itc = tmp;
	tmp = iva;
	iva = ivc;
	ivc = tmp;
    }

/*     assign estimated angle bend force constants */

    if (itb == 6) {
	if (ita == 1) {
	    if (ivb == 4) {
		if (itc == 1) {
		    ret_val = 34.5;
		} else if (itc == 6) {
		    ret_val = 38.;
		} else if (itc == 7) {
		    ret_val = 50.6;
		} else if (itc == 8) {
		    ret_val = 51.5;
		} else if (itc == 9) {
		    ret_val = 50.;
		} else {
		    ret_val = 35.;
		}
	    } else if (ivb == 3) {
		ret_val = 32.;
	    } else {
		ret_val = 32.;
	    }
	} else if (ita == 6) {
	    if (ivb == 4) {
		if (itc == 6) {
		    ret_val = 60.;
		} else if (itc == 7) {
		    ret_val = 80.;
		} else if (itc == 8) {
		    ret_val = 88.;
		} else if (itc == 9) {
		    ret_val = 89.;
		} else if (itc == 14) {
		    ret_val = 65.;
		} else if (itc == 15) {
		    ret_val = 60.;
		} else if (itc == 16) {
		    ret_val = 53.2;
		} else if (itc == 17) {
		    ret_val = 55.;
		} else {
		    ret_val = 50.;
		}
	    } else if (ivb == 3) {
		ret_val = 60.;
	    } else {
		ret_val = 60.;
	    }
	} else if (ita == 8) {
	    if (ivb == 4) {
		if (itc == 8) {
		    ret_val = 65.;
		} else if (itc == 9) {
		    ret_val = 65.;
		} else if (itc == 15) {
		    ret_val = 60.;
		} else if (itc == 16) {
		    ret_val = 65.;
		} else {
		    ret_val = 65.;
		}
	    } else if (ivb == 3) {
		ret_val = 50.;
	    } else {
		ret_val = 60.;
	    }
	} else {
	    ret_val = 60.;
	}
    } else if (itb == 8) {
	if (ita == 1) {
	    if (itc == 1) {
		ret_val = 34.05;
	    } else if (itc == 6) {
		ret_val = 65.;
	    } else {
		ret_val = 60.;
	    }
	} else if (ita == 6) {
	    if (itc == 6) {
		ret_val = 88.5;
	    } else if (itc == 8) {
		if (iva == 1 || ivc == 1) {
		    ret_val = 122.3;
		} else {
		    ret_val = 85.;
		}
	    } else if (itc == 15) {
		ret_val = 80.3;
	    } else {
		ret_val = 80.;
	    }
	} else {
	    ret_val = 80.;
	}
    } else if (itb == 15) {
	if (ita == 1) {
	    ret_val = 30.;
	} else if (ita == 6) {
	    if (itc == 6) {
		ret_val = 75.;
	    } else if (itc == 8) {
		ret_val = 80.;
	    } else {
		ret_val = 75.;
	    }
	} else if (ita == 8) {
	    if (itc == 8) {
		if (iva == 1 && ivc == 1) {
		    ret_val = 89.88;
		} else if (iva == 1 || ivc == 1) {
		    ret_val = 75.86;
		} else {
		    ret_val = 65.58;
		}
	    } else {
		ret_val = 70.;
	    }
	} else {
	    ret_val = 75.;
	}
    } else if (itb == 16) {
	if (ita == 1) {
	    ret_val = 30.;
	} else if (ita == 6) {
	    if (itc == 16) {
		ret_val = 72.;
	    } else {
		ret_val = 80.;
	    }
	} else if (ita == 8) {
	    if (itc == 8) {
		if (iva == 1 && ivc == 1) {
		    ret_val = 168.;
		} else if (iva == 1 || ivc == 1) {
		    ret_val = 85.;
		} else {
		    ret_val = 80.;
		}
	    } else if (itc == 16) {
		ret_val = 75.;
	    } else {
		ret_val = 75.;
	    }
	} else {
	    ret_val = 75.;
	}
    } else if (ita == 1) {
	ret_val = 35.;
    } else {
	ret_val = 65.;
    }

/*     scale the force constant to the desired units */

    ret_val /= angpot_1.angunit * 3282.8063500117441;
    return ret_val;
} /* angguess_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine sbguess  --  estimate stretch-bend parameters  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "sbguess" sets approximate stretch-bend force constants based */
/*     on atom type and connected atoms */


/* Subroutine */ int sbguess_(integer *ia, integer *ib, integer *ic, 
	doublereal *sb1, doublereal *sb2)
{
    static integer ita, itb, itc, iva, ivb, ivc;



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  angpot.i  --  specifics of angle bend functional forms  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     angunit    convert angle bending energy to kcal/mole */
/*     stbnunit   convert stretch-bend energy to kcal/mole */
/*     aaunit     convert angle-angle energy to kcal/mole */
/*     opbunit    convert out-of-plane bend energy to kcal/mole */
/*     opdunit    convert out-of-plane distance energy to kcal/mole */
/*     cang       cubic coefficient in angle bending potential */
/*     qang       quartic coefficient in angle bending potential */
/*     pang       quintic coefficient in angle bending potential */
/*     sang       sextic coefficient in angle bending potential */
/*     copb       cubic coefficient in out-of-plane bend potential */
/*     qopb       quartic coefficient in out-of-plane bend potential */
/*     popb       quintic coefficient in out-of-plane bend potential */
/*     sopb       sextic coefficient in out-of-plane bend potential */
/*     copd       cubic coefficient in out-of-plane distance potential */
/*     qopd       quartic coefficient in out-of-plane distance potential */
/*     popd       quintic coefficient in out-of-plane distance potential */
/*     sopd       sextic coefficient in out-of-plane distance potential */
/*     angtyp     type of angle bending function for each bond angle */
/*     opbtyp     type of out-of-plane bend potential energy function */




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




/*     get the atomic number and valence of each atom */

    ita = atmtyp_1.atomic[*ia - 1];
    itb = atmtyp_1.atomic[*ib - 1];
    itc = atmtyp_1.atomic[*ic - 1];
    iva = atmtyp_1.valence[*ia - 1];
    ivb = atmtyp_1.valence[*ib - 1];
    ivc = atmtyp_1.valence[*ic - 1];

/*     set initial stretch-bend parameters */

    if (ita == 1 && itc == 1) {
	*sb1 = 0.;
	*sb2 = 0.;
    } else if (itb == 6) {
	if (ita == 1) {
	    *sb1 = 11.5;
	    *sb2 = 18.7;
	} else if (itc == 1) {
	    *sb1 = 18.7;
	    *sb2 = 11.5;
	} else {
	    *sb1 = 18.7;
	    *sb2 = 18.7;
	}
    } else if (itb == 6) {
	if (ita == 1) {
	    *sb1 = 4.5;
	    *sb2 = 12.95;
	} else if (itc == 1) {
	    *sb1 = 12.95;
	    *sb2 = 4.5;
	} else {
	    *sb1 = 14.4;
	    *sb2 = 14.4;
	}
    } else if (itb == 7) {
	if (ivb >= 3) {
	    if (ita == 1) {
		*sb1 = 4.3;
		*sb2 = 7.2;
	    } else if (itc == 1) {
		*sb1 = 7.2;
		*sb2 = 4.3;
	    } else {
		*sb1 = 7.2;
		*sb2 = 7.2;
	    }
	} else {
	    if (ita == 1) {
		*sb1 = 4.3;
		*sb2 = 14.4;
	    } else if (itc == 1) {
		*sb1 = 14.4;
		*sb2 = 4.3;
	    } else {
		*sb1 = 14.4;
		*sb2 = 14.4;
	    }
	}
    } else if (itb == 14) {
	if (ita == 1) {
	    *sb1 = 8.6;
	    *sb2 = 14.4;
	} else if (itc == 1) {
	    *sb1 = 14.4;
	    *sb2 = 8.6;
	} else {
	    *sb1 = 14.4;
	    *sb2 = 14.4;
	}
    } else if (itb == 15) {
	if (ivb == 4) {
	    if (ita == 1) {
		*sb1 = 14.4;
		*sb2 = 14.4;
	    } else if (itc == 1) {
		*sb1 = 14.4;
		*sb2 = 14.4;
	    } else {
		*sb1 = 14.4;
		*sb2 = 14.4;
	    }
	} else {
	    if (ita == 1) {
		*sb1 = 8.6;
		*sb2 = 8.6;
	    } else if (itc == 1) {
		*sb1 = 8.6;
		*sb2 = 8.6;
	    } else {
		*sb1 = 8.6;
		*sb2 = 8.6;
	    }
	}
    } else if (itb == 16) {
	if (ita == 1) {
	    *sb1 = 1.45;
	    *sb2 = -5.75;
	} else if (itc == 1) {
	    *sb1 = -5.75;
	    *sb2 = 1.45;
	} else {
	    *sb1 = -5.75;
	    *sb2 = -5.75;
	}
    } else if (ita == 1 && itc > 1) {
	*sb1 = -4.5;
	*sb2 = 38.;
    } else if (ita > 1 && itc == 1) {
	*sb1 = 38.;
	*sb2 = -4.5;
    } else if (ita > 1 && itc > 1) {
	*sb1 = 38.;
	*sb2 = 38.;
    } else {
	*sb1 = 38.;
	*sb2 = 38.;
    }

/*     scale the force constant to the desired units */

    *sb1 /= angpot_1.stbnunit * 57.29577951308232088;
    *sb2 /= angpot_1.stbnunit * 57.29577951308232088;
    return 0;
} /* sbguess_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  function uryguess  --  estimate Urey-Bradley parameters  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "uryguess" sets approximate Urey-Bradley force constants */
/*     based on atom type and connected atoms */


doublereal uryguess_(integer *ia, integer *ib, integer *ic)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer ita, itb, itc, iva, ivb, ivc;



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
/*     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  urypot.i  --  specifics of Urey-Bradley functional form  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     cury       cubic coefficient in Urey-Bradley potential */
/*     qury       quartic coefficient in Urey-Bradley potential */
/*     ureyunit   convert Urey-Bradley energy to kcal/mole */




/*     get the atomic number and valence of each atom */

    ita = atmtyp_1.atomic[*ia - 1];
    itb = atmtyp_1.atomic[*ib - 1];
    itc = atmtyp_1.atomic[*ic - 1];
    iva = atmtyp_1.valence[*ia - 1];
    ivb = atmtyp_1.valence[*ib - 1];
    ivc = atmtyp_1.valence[*ic - 1];

/*     assign estimated out-of-plane parameter values */

    ret_val = 10.;

/*     scale the force constant to the desired units */

    ret_val /= urypot_1.ureyunit;
    return ret_val;
} /* uryguess_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function opbguess  --  estimate out-of-plane bend values  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "opbguess" sets approximate out-of-plane bend force constants */
/*     based on atom type and connected atoms */


doublereal opbguess_(integer *ia, integer *ib, integer *ic, integer *id)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer ita, itb, iva, ivb;



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  angpot.i  --  specifics of angle bend functional forms  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     angunit    convert angle bending energy to kcal/mole */
/*     stbnunit   convert stretch-bend energy to kcal/mole */
/*     aaunit     convert angle-angle energy to kcal/mole */
/*     opbunit    convert out-of-plane bend energy to kcal/mole */
/*     opdunit    convert out-of-plane distance energy to kcal/mole */
/*     cang       cubic coefficient in angle bending potential */
/*     qang       quartic coefficient in angle bending potential */
/*     pang       quintic coefficient in angle bending potential */
/*     sang       sextic coefficient in angle bending potential */
/*     copb       cubic coefficient in out-of-plane bend potential */
/*     qopb       quartic coefficient in out-of-plane bend potential */
/*     popb       quintic coefficient in out-of-plane bend potential */
/*     sopb       sextic coefficient in out-of-plane bend potential */
/*     copd       cubic coefficient in out-of-plane distance potential */
/*     qopd       quartic coefficient in out-of-plane distance potential */
/*     popd       quintic coefficient in out-of-plane distance potential */
/*     sopd       sextic coefficient in out-of-plane distance potential */
/*     angtyp     type of angle bending function for each bond angle */
/*     opbtyp     type of out-of-plane bend potential energy function */




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




/*     get the atomic number and valence of each atom */

    ita = atmtyp_1.atomic[*ia - 1];
    itb = atmtyp_1.atomic[*ib - 1];
    iva = atmtyp_1.valence[*ia - 1];
    ivb = atmtyp_1.valence[*ib - 1];

/*     assign estimated out-of-plane parameter values */

    ret_val = 14.4;

/*     scale the force constant to the desired units */

    ret_val /= angpot_1.opbunit * 3282.8063500117441;
    return ret_val;
} /* opbguess_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine torguess  --  estimate torsional parameters  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "torguess" set approximate torsion amplitude parameters based */
/*     on atom type and connected atoms */


/* Subroutine */ int torguess_(integer *ia, integer *ib, integer *ic, integer 
	*id, doublereal *tf1, doublereal *tf2, doublereal *tf3)
{
    static integer ita, itb, itc, itd, iva, ivb, ivc, ivd, tmp;



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  torpot.i  --  specifics of torsional functional forms  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     idihunit  convert improper dihedral energy to kcal/mole */
/*     itorunit  convert improper torsion amplitudes to kcal/mole */
/*     torsunit  convert torsional parameter amplitudes to kcal/mole */
/*     ptorunit  convert pi-orbital torsion energy to kcal/mole */
/*     storunit  convert stretch-torsion energy to kcal/mole */
/*     ttorunit  convert stretch-torsion energy to kcal/mole */




/*     get the atomic number and valence of each atom */

    ita = atmtyp_1.atomic[*ia - 1];
    itb = atmtyp_1.atomic[*ib - 1];
    itc = atmtyp_1.atomic[*ic - 1];
    itd = atmtyp_1.atomic[*id - 1];
    iva = atmtyp_1.valence[*ia - 1];
    ivb = atmtyp_1.valence[*ib - 1];
    ivc = atmtyp_1.valence[*ic - 1];
    ivd = atmtyp_1.valence[*id - 1];

/*     reorder the atoms based on the atomic numbers */

    if (itb > itc || itb == itc && ita > itd) {
	tmp = itb;
	itb = itc;
	itc = tmp;
	tmp = ivb;
	ivb = ivc;
	ivc = tmp;
	tmp = ita;
	ita = itd;
	itd = tmp;
	tmp = iva;
	iva = ivd;
	ivd = tmp;
    }

/*     assign estimated torsional parameter values */

    *tf1 = 0.;
    *tf2 = 0.;
    *tf3 = 0.;
    if (itb == 6 && itc == 6) {
	if (ita == 6 && itd == 6) {
	    if (ivb == 3 && ivc == 3) {
		if (iva == 3 && ivd == 3) {
		    *tf1 = -.335;
		    *tf2 = 2.;
		    *tf3 = 0.;
		} else if (iva == 3 && ivd == 4) {
		    *tf1 = -.305;
		    *tf2 = 2.105;
		    *tf3 = 0.;
		} else if (iva == 4 && ivd == 4) {
		    *tf1 = 0.;
		    *tf2 = 4.;
		    *tf3 = 0.;
		}
	    } else if (ivb == 3 && ivc == 4) {
		*tf1 = -.4;
		*tf2 = -.05;
		*tf3 = -.275;
	    } else if (ivb == 4 && ivc == 4) {
		*tf1 = .09;
		*tf2 = .085;
		*tf3 = .26;
	    }
	} else if (ita == 1 && itd == 1) {
	    if (ivb == 3 && ivc == 3) {
		*tf1 = 0.;
		*tf2 = 2.035;
		*tf3 = 0.;
	    } else {
		*tf1 = 0.;
		*tf2 = 0.;
		*tf3 = .15;
	    }
	} else if (ita == 1 && itd == 6) {
	    if (ivb == 4 && ivc == 4) {
		*tf1 = 0.;
		*tf2 = 0.;
		*tf3 = .17;
	    } else if (ivb == 3 && ivc == 3 && ivd == 3) {
		*tf1 = 0.;
		*tf2 = 3.05;
		*tf3 = 0.;
	    } else if (ivb == 3 && ivc == 3 && ivd == 4) {
		*tf1 = 0.;
		*tf2 = 3.05;
		*tf3 = 0.;
	    } else if (ivb == 4 && ivc == 3) {
		*tf1 = 0.;
		*tf2 = 0.;
		*tf3 = -.045;
	    }
	} else if (ita == 1 && itd == 7) {
	    if (ivb == 3 && ivc == 3) {
		*tf1 = -1.575;
		*tf2 = 1.5;
		*tf3 = 0.;
	    } else {
		*tf1 = 0.;
		*tf2 = 0.;
		*tf3 = .25;
	    }
	} else if (ita == 1 && itd == 8) {
	    *tf1 = 0.;
	    *tf2 = 0.;
	    *tf3 = .15;
	} else if (ita == 6 && itd == 8) {
	    if (ivb == 3 && ivc == 3) {
		*tf1 = 0.;
		*tf2 = 2.235;
		*tf3 = 0.;
	    } else {
		*tf1 = -.575;
		*tf2 = 0.;
		*tf3 = .64;
	    }
	} else if (ita == 8 && itd == 8) {
	    *tf1 = 1.11;
	    *tf2 = -.69;
	    *tf3 = -.59;
	} else if (ivb == 3 && ivc == 3) {
	    *tf1 = 0.;
	    *tf2 = 1.25;
	    *tf3 = 0.;
	} else {
	    *tf1 = 0.;
	    *tf2 = 0.;
	    *tf3 = .15;
	}
    } else if (itb == 6 && itc == 8) {
	if (ita == 1 && itd == 1) {
	    *tf1 = 0.;
	    *tf2 = 0.;
	    *tf3 = .135;
	} else if (ita == 1 && itd == 6) {
	    if (ivc == 3 && ivd == 3) {
		*tf1 = 0.;
		*tf2 = 2.235;
		*tf3 = 0.;
	    } else {
		*tf1 = 0.;
		*tf2 = 0.;
		*tf3 = .355;
	    }
	} else if (ita == 1) {
	    *tf1 = 0.;
	    *tf2 = 0.;
	    *tf3 = .375;
	} else if (ita == 6 && itd == 1 && ivb == 4) {
	    *tf1 = -.885;
	    *tf2 = .115;
	    *tf3 = .38;
	} else if (ita == 6 && itd == 6) {
	    *tf1 = 1.;
	    *tf2 = -.75;
	    *tf3 = .445;
	} else if (ita == 6 && itd == 1 && ivb == 3) {
	    *tf1 = 0.;
	    *tf2 = 1.175;
	    *tf3 = 0.;
	} else if (ivb == 3) {
	    *tf1 = 0.;
	    *tf2 = 1.25;
	    *tf3 = 0.;
	} else if (ivb == 4) {
	    *tf1 = 1.;
	    *tf2 = -.75;
	    *tf3 = .445;
	}
    } else if (itb == 6 && itc == 15) {
	*tf1 = 0.;
	*tf2 = 1.25;
	*tf3 = .25;
    } else if (itb == 6 && itc == 16) {
	*tf1 = 0.;
	*tf2 = 0.;
	*tf3 = .25;
    } else if (itb == 8 && itc == 15) {
	*tf1 = -1.;
	*tf2 = -.84;
	*tf3 = -.4;
    } else if (itb == 8 && itc == 16) {
	*tf1 = -.75;
	*tf2 = -1.;
	*tf3 = -.4;
    } else {
	*tf1 = 0.;
	*tf2 = .5;
	*tf3 = .25;
    }

/*     scale the amplitude values to the desired units */

    *tf1 /= torpot_1.torsunit;
    *tf2 /= torpot_1.torsunit;
    *tf3 /= torpot_1.torsunit;
    return 0;
} /* torguess_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  function valrms  --  compute structure & vibration RMSD  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "valrms" evaluates a valence parameter goodness-of-fit error */
/*     function based on comparison of forces, frequencies, bond */
/*     lengths and angles to QM results */


doublereal valrms_(integer *prtflg)
{
    /* Initialized data */

    static char axis[1*3] = "X" "Y" "Z";

    /* Format strings */
    static char fmt_10[] = "(/,\002 Comparison of Bond Lengths :\002,//,6x"
	    ",\002Bond\002,8x,\002Atoms\002,19x,\002QM Bond\002,6x,\002MM Bond"
	    "\002,8x,\002Delta\002,/)";
    static char fmt_20[] = "(4x,i5,4x,2i5,13x,3f13.4)";
    static char fmt_30[] = "(/,4x,\002Average Unsigned Difference :\002,30x,"
	    "f12.4,/,4x,\002Root Mean Square Deviation :\002,31x,f12.4)";
    static char fmt_40[] = "(/,\002 Comparison of Bond Angles :\002,//,5x"
	    ",\002Angle\002,10x,\002Atoms\002,16x,\002QM Angle\002,5x,\002MM "
	    "Angle\002,8x,\002Delta\002,/)";
    static char fmt_50[] = "(4x,i5,4x,3i5,8x,2f13.2,f13.4)";
    static char fmt_60[] = "(/,4x,\002Average Unsigned Difference :\002,30x,"
	    "f12.4,/,4x,\002Root Mean Square Deviation :\002,31x,f12.4)";
    static char fmt_70[] = "(/,\002 Comparison of Torsion Angles :\002,//,"
	    "4x,\002Torsion\002,12x,\002Atoms\002,13x,\002QM Angle\002,5x,"
	    "\002MM Angle\002,8x,\002Delta\002,/)";
    static char fmt_80[] = "(4x,i5,4x,4i5,3x,2f13.2,f13.4)";
    static char fmt_90[] = "(/,4x,\002Average Unsigned Difference :\002,30x,"
	    "f12.4,/,4x,\002Root Mean Square Deviation :\002,31x,f12.4)";
    static char fmt_100[] = "(/,\002 Comparison of Gradient Components :\002"
	    ",//,7x,\002Atom\002,14x,\002QM Grad\002,8x,\002MM Grad\002,10x"
	    ",\002Delta\002,/)";
    static char fmt_110[] = "(4x,i5,1x,a1,8x,f13.4,2x,f13.4,2x,f13.4)";
    static char fmt_120[] = "(/,4x,\002Average Unsigned Difference :\002,17x"
	    ",f12.4,/,4x,\002Root Mean Square Deviation :\002,18x,f12.4)";
    static char fmt_130[] = "(/,\002 Comparison of Hessian Elements :\002,//"
	    ",7x,\002Atom\002,14x,\002QM Hess\002,8x,\002MM Hess\002,10x,\002"
	    "Delta\002,/)";
    static char fmt_140[] = "(4x,i5,1x,a1,8x,f13.2,2x,f13.2,2x,f13.4)";
    static char fmt_150[] = "(/,4x,\002Average Unsigned Difference :\002,17x"
	    ",f12.4,/,4x,\002Root Mean Square Deviation :\002,18x,f12.4)";
    static char fmt_160[] = "(/,\002 VIBRATE  --  Too many Atoms in the Mole"
	    "cule\002)";
    static char fmt_170[] = "(/,\002 Comparison of Vibrational Frequencies "
	    ":\002,//,6x,\002Mode\002,15x,\002QM Freq\002,8x,\002MM Freq\002,"
	    "10x,\002Delta\002,/)";
    static char fmt_180[] = "(4x,i5,10x,f13.2,2x,f13.2,2x,f13.4)";
    static char fmt_190[] = "(/,4x,\002Average Unsigned Difference :\002,17x"
	    ",f12.4,/,4x,\002Root Mean Square Deviation :\002,18x,f12.4)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), e_wsfe(void);
    double sqrt(doublereal);
    integer do_fio(integer *, char *, ftnlen);
    double acos(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    extern doublereal minimiz1_();
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static doublereal a[1000], b[1000], h__[1000000];
    static integer i__, j, k, m;
    static doublereal p[1000], w[1000];
    static integer m1, m2, ia, ib, ic, id;
    static doublereal ta[1000], tb[1000], xt, yt, zt, xu, yu, zu, xx[75000], 
	    ty[1000], rt2, ru2, rcb, xab, yab, zab, xba, yba, zba, xcb, ycb, 
	    zcb, xdc, ydc, zdc, xtu, ytu, ztu, rab2, rcb2, afac, bfac, ffac, 
	    gfac, hfac, rabc, bave, aave, bond, tfac, gave, have, sine, tave, 
	    fave, arms, brms, fcut;
    static integer nvar;
    static doublereal grms, hrms, frms, trms, rtru, mass2[25000], hdiag[75000]
	    	/* was [3][25000] */;
    extern /* Subroutine */ int diagq_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal gbond, angle, delta, eigen[1000];
    extern /* Subroutine */ int lbfgs_(integer *, doublereal *, doublereal *, 
	    doublereal *, D_fp, U_fp), fatal_(void);
    static integer ihess, nfreq, hinit[75000]	/* was [3][25000] */;
    static doublereal vects[1000000]	/* was [1000][1000] */;
    static integer hstop[75000]	/* was [3][25000] */;
    static doublereal gangle, factor;
    static integer oldprt, oldwrt, hindex[1000000];
    static doublereal grdmin, cosine, energy, matrix[500500], derivs[75000]	
	    /* was [3][25000] */;
    extern /* Subroutine */ int hessian_(doublereal *, integer *, integer *, 
	    integer *, doublereal *);
    static integer olditer, oldstep;
    static doublereal minimum;
    extern /* Subroutine */ int optsave_();

    /* Fortran I/O blocks */
    static cilist io___148 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___157 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___158 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___161 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___172 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___173 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___176 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___198 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___199 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___204 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___206 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___207 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___215 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___218 = { 0, 6, 0, fmt_140, 0 };
    static cilist io___221 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___223 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___240 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___241 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___242 = { 0, 0, 0, fmt_190, 0 };



#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define hdiag_ref(a_1,a_2) hdiag[(a_2)*3 + a_1 - 4]
#define hinit_ref(a_1,a_2) hinit[(a_2)*3 + a_1 - 4]
#define hstop_ref(a_1,a_2) hstop[(a_2)*3 + a_1 - 4]
#define itors_ref(a_1,a_2) tors_1.itors[(a_2)*4 + a_1 - 5]
#define gforce_ref(a_1,a_2) qmstuf_1.gforce[(a_2)*3 + a_1 - 4]
#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1 - 4]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  hescut.i  --  cutoff value for Hessian matrix elements  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     hesscut   magnitude of smallest allowed Hessian element */




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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kangs.i  --  forcefield parameters for bond angle bending  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxna    maximum number of harmonic angle bend parameter entries */
/*     maxna5   maximum number of 5-membered ring angle bend entries */
/*     maxna4   maximum number of 4-membered ring angle bend entries */
/*     maxna3   maximum number of 3-membered ring angle bend entries */
/*     maxnaf   maximum number of Fourier angle bend parameter entries */

/*     acon     force constant parameters for harmonic angle bends */
/*     acon5    force constant parameters for 5-ring angle bends */
/*     acon4    force constant parameters for 4-ring angle bends */
/*     acon3    force constant parameters for 3-ring angle bends */
/*     aconf    force constant parameters for Fourier angle bends */
/*     ang      bond angle parameters for harmonic angle bends */
/*     ang5     bond angle parameters for 5-ring angle bends */
/*     ang4     bond angle parameters for 4-ring angle bends */
/*     ang3     bond angle parameters for 3-ring angle bends */
/*     angf     phase shift angle and periodicity for Fourier bends */
/*     ka       string of atom classes for harmonic angle bends */
/*     ka5      string of atom classes for 5-ring angle bends */
/*     ka4      string of atom classes for 4-ring angle bends */
/*     ka3      string of atom classes for 3-ring angle bends */
/*     kaf      string of atom classes for Fourier angle bends */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  kbonds.i  --  forcefield parameters for bond stretching  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxnb   maximum number of bond stretch parameter entries */
/*     maxnb5  maximum number of 5-membered ring bond stretch entries */
/*     maxnb4  maximum number of 4-membered ring bond stretch entries */
/*     maxnb3  maximum number of 3-membered ring bond stretch entries */
/*     maxnel  maximum number of electronegativity bond corrections */

/*     bcon    force constant parameters for harmonic bond stretch */
/*     blen    bond length parameters for harmonic bond stretch */
/*     bcon5   force constant parameters for 5-ring bond stretch */
/*     blen5   bond length parameters for 5-ring bond stretch */
/*     bcon4   force constant parameters for 4-ring bond stretch */
/*     blen4   bond length parameters for 4-ring bond stretch */
/*     bcon3   force constant parameters for 3-ring bond stretch */
/*     blen3   bond length parameters for 3-ring bond stretch */
/*     dlen    electronegativity bond length correction parameters */
/*     kb      string of atom classes for harmonic bond stretch */
/*     kb5     string of atom classes for 5-ring bond stretch */
/*     kb4     string of atom classes for 4-ring bond stretch */
/*     kb3     string of atom classes for 3-ring bond stretch */
/*     kel     string of atom classes for electronegativity corrections */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kopbnd.i  --  forcefield parameters for out-of-plane bend  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxnopb   maximum number of out-of-plane bending entries */

/*     opbn      force constant parameters for out-of-plane bending */
/*     kopb      string of atom classes for out-of-plane bending */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  kstbnd.i  --  forcefield parameters for stretch-bend  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxnsb   maximum number of stretch-bend parameter entries */

/*     stbn     force constant parameters for stretch-bend terms */
/*     ksb      string of atom classes for stretch-bend terms */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  ktorsn.i  --  forcefield parameters for torsional angles  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxnt    maximum number of torsional angle parameter entries */
/*     maxnt5   maximum number of 5-membered ring torsion entries */
/*     maxnt4   maximum number of 4-membered ring torsion entries */

/*     t1       torsional parameters for standard 1-fold rotation */
/*     t2       torsional parameters for standard 2-fold rotation */
/*     t3       torsional parameters for standard 3-fold rotation */
/*     t4       torsional parameters for standard 4-fold rotation */
/*     t5       torsional parameters for standard 5-fold rotation */
/*     t6       torsional parameters for standard 6-fold rotation */
/*     t15      torsional parameters for 1-fold rotation in 5-ring */
/*     t25      torsional parameters for 2-fold rotation in 5-ring */
/*     t35      torsional parameters for 3-fold rotation in 5-ring */
/*     t45      torsional parameters for 4-fold rotation in 5-ring */
/*     t55      torsional parameters for 5-fold rotation in 5-ring */
/*     t65      torsional parameters for 6-fold rotation in 5-ring */
/*     t14      torsional parameters for 1-fold rotation in 4-ring */
/*     t24      torsional parameters for 2-fold rotation in 4-ring */
/*     t34      torsional parameters for 3-fold rotation in 4-ring */
/*     t44      torsional parameters for 4-fold rotation in 4-ring */
/*     t54      torsional parameters for 5-fold rotation in 4-ring */
/*     t64      torsional parameters for 6-fold rotation in 4-ring */
/*     kt       string of atom classes for torsional angles */
/*     kt5      string of atom classes for 5-ring torsions */
/*     kt4      string of atom classes for 4-ring torsions */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kvdws.i  --  forcefield parameters for van der Waals terms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     rad      van der Waals radius parameter for each atom type */
/*     eps      van der Waals well depth parameter for each atom type */
/*     rad4     van der Waals radius parameter in 1-4 interactions */
/*     eps4     van der Waals well depth parameter in 1-4 interactions */
/*     reduct   van der Waals reduction factor for each atom type */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  linmin.i  --  parameters for line search minimization  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     stpmin   minimum step length in current line search direction */
/*     stpmax   maximum step length in current line search direction */
/*     cappa    stringency of line search (0=tight < cappa < 1=loose) */
/*     slpmax   projected gradient above which stepsize is reduced */
/*     angmax   maximum angle between search direction and -gradient */
/*     intmax   maximum number of interpolations during line search */




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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  opbend.i  --  out-of-plane bends in the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     opbk      force constant values for out-of-plane bending */
/*     nopbend   total number of out-of-plane bends in the system */
/*     iopb      bond angle numbers used in out-of-plane bending */




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




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2009 by Chuanjie Wu and Jay William Ponder  ## */
/*     ##                    All Rights Reserved                     ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  qmstuf.i  --  quantum data from Gaussian 03 calculation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     gx       x-coordinate of each atom in the QM data file */
/*     gy       y-coordinate of each atom in the QM data file */
/*     gz       z-coordinate of each atom in the QM data file */
/*     gforce   force components on each atom from QM data */
/*     gh       Hessian maxtrix elements from QM data */
/*     gfreq    calculated vibrational frequencies from QM data */
/*     ngatom   number of atoms in the QM data file */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  scales.i  --  parameter scale factors for optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     scale      multiplicative factor for each optimization parameter */
/*     set_scale  logical flag to show if scale factors have been set */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  strbnd.i  --  stretch-bends in the current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     sbk       force constants for stretch-bend terms */
/*     nstrbnd   total number of stretch-bend interactions */
/*     isb       angle and bond numbers used in stretch-bend */




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




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  valfit.i  --  values for valence term parameter fitting  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     fit_bond    logical flag to fit bond stretch parameters */
/*     fit_angle   logical flag to fit angle bend parameters */
/*     fit_strbnd  logical flag to fit stretch-bend parameters */
/*     fit_urey    logical flag to fit Urey-Bradley parameters */
/*     fit_opbend  logical flag to fit out-of-plane bend parameters */
/*     fit_tors    logical flag to fit torsional parameters */
/*     fit_struct  logical flag to structure-fit valence parameters */
/*     fit_force   logical flag to force-fit valence parameters */




/*     scale the coordinates of each active atom; use the */
/*     square root of median eigenvalue of typical Hessian */

    if (valfit_1.fit_struct__) {
	scales_1.set_scale__ = TRUE_;
	nvar = 0;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ++nvar;
	    scales_1.scale[nvar - 1] = 12.;
	    xx[nvar - 1] = qmstuf_1.gx[i__ - 1] * scales_1.scale[nvar - 1];
	    ++nvar;
	    scales_1.scale[nvar - 1] = 12.;
	    xx[nvar - 1] = qmstuf_1.gy[i__ - 1] * scales_1.scale[nvar - 1];
	    ++nvar;
	    scales_1.scale[nvar - 1] = 12.;
	    xx[nvar - 1] = qmstuf_1.gz[i__ - 1] * scales_1.scale[nvar - 1];
	}

/*     make the call to the optimization routine */

	oldstep = (integer) linmin_1.stpmax;
	olditer = minima_1.maxiter;
	oldprt = inform_1.iprint;
	oldwrt = inform_1.iwrite;
	linmin_1.stpmax = 0.;
	minima_1.maxiter = 0;
	inform_1.iprint = 0;
	inform_1.iwrite = 0;
	grdmin = 1e-4;
	s_copy(output_1.coordtype, "CARTESIAN", (ftnlen)9, (ftnlen)9);
	lbfgs_(&nvar, xx, &minimum, &grdmin, (D_fp)minimiz1_, (U_fp)optsave_);
	s_copy(output_1.coordtype, "NONE", (ftnlen)9, (ftnlen)4);
	linmin_1.stpmax = (doublereal) oldstep;
	minima_1.maxiter = olditer;
	inform_1.iprint = oldprt;
	inform_1.iwrite = oldwrt;

/*     unscale the final coordinates for active atoms */

	nvar = 0;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ++nvar;
	    atoms_1.x[i__ - 1] = xx[nvar - 1] / scales_1.scale[nvar - 1];
	    scales_1.scale[nvar - 1] = 1.;
	    ++nvar;
	    atoms_1.y[i__ - 1] = xx[nvar - 1] / scales_1.scale[nvar - 1];
	    scales_1.scale[nvar - 1] = 1.;
	    ++nvar;
	    atoms_1.z__[i__ - 1] = xx[nvar - 1] / scales_1.scale[nvar - 1];
	    scales_1.scale[nvar - 1] = 1.;
	}
    }

/*     compute the RMS between QM and TINKER bond lengths */

    bave = 0.;
    brms = 0.;
    if (valfit_1.fit_struct__) {
	if (*prtflg == 1 && bond_1.nbond != 0) {
	    io___148.ciunit = iounit_1.iout;
	    s_wsfe(&io___148);
	    e_wsfe();
	}
	i__1 = bond_1.nbond;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = ibnd_ref(1, i__);
	    ib = ibnd_ref(2, i__);
	    xab = atoms_1.x[ia - 1] - atoms_1.x[ib - 1];
	    yab = atoms_1.y[ia - 1] - atoms_1.y[ib - 1];
	    zab = atoms_1.z__[ia - 1] - atoms_1.z__[ib - 1];
	    bond = sqrt(xab * xab + yab * yab + zab * zab);
	    xab = qmstuf_1.gx[ia - 1] - qmstuf_1.gx[ib - 1];
	    yab = qmstuf_1.gy[ia - 1] - qmstuf_1.gy[ib - 1];
	    zab = qmstuf_1.gz[ia - 1] - qmstuf_1.gz[ib - 1];
	    gbond = sqrt(xab * xab + yab * yab + zab * zab);
	    delta = bond - gbond;
	    bave += abs(delta);
	    brms += delta * delta;
	    if (*prtflg == 1) {
		io___157.ciunit = iounit_1.iout;
		s_wsfe(&io___157);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&gbond, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bond, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&delta, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
	if (bond_1.nbond != 0) {
	    bave /= (doublereal) bond_1.nbond;
	}
	if (bond_1.nbond != 0) {
	    brms = sqrt(brms / (doublereal) bond_1.nbond);
	}
	if (*prtflg == 1 && bond_1.nbond != 0) {
	    io___158.ciunit = iounit_1.iout;
	    s_wsfe(&io___158);
	    do_fio(&c__1, (char *)&bave, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&brms, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     compute the RMS between QM and TINKER bond angles */

    aave = 0.;
    arms = 0.;
    if (valfit_1.fit_struct__) {
	if (*prtflg == 1 && angle_1.nangle != 0) {
	    io___161.ciunit = iounit_1.iout;
	    s_wsfe(&io___161);
	    e_wsfe();
	}
	i__1 = angle_1.nangle;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = iang_ref(1, i__);
	    ib = iang_ref(2, i__);
	    ic = iang_ref(3, i__);
	    xab = atoms_1.x[ia - 1] - atoms_1.x[ib - 1];
	    yab = atoms_1.y[ia - 1] - atoms_1.y[ib - 1];
	    zab = atoms_1.z__[ia - 1] - atoms_1.z__[ib - 1];
	    xcb = atoms_1.x[ic - 1] - atoms_1.x[ib - 1];
	    ycb = atoms_1.y[ic - 1] - atoms_1.y[ib - 1];
	    zcb = atoms_1.z__[ic - 1] - atoms_1.z__[ib - 1];
	    rab2 = xab * xab + yab * yab + zab * zab;
	    rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
	    rabc = sqrt(rab2 * rcb2);
	    if (rabc != 0.) {
		cosine = (xab * xcb + yab * ycb + zab * zcb) / rabc;
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine);
		cosine = min(d__1,d__2);
		angle = acos(cosine) * 57.29577951308232088;
	    }
	    xab = qmstuf_1.gx[ia - 1] - qmstuf_1.gx[ib - 1];
	    yab = qmstuf_1.gy[ia - 1] - qmstuf_1.gy[ib - 1];
	    zab = qmstuf_1.gz[ia - 1] - qmstuf_1.gz[ib - 1];
	    xcb = qmstuf_1.gx[ic - 1] - qmstuf_1.gx[ib - 1];
	    ycb = qmstuf_1.gy[ic - 1] - qmstuf_1.gy[ib - 1];
	    zcb = qmstuf_1.gz[ic - 1] - qmstuf_1.gz[ib - 1];
	    rab2 = xab * xab + yab * yab + zab * zab;
	    rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
	    rabc = sqrt(rab2 * rcb2);
	    if (rabc != 0.) {
		cosine = (xab * xcb + yab * ycb + zab * zcb) / rabc;
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine);
		cosine = min(d__1,d__2);
		gangle = acos(cosine) * 57.29577951308232088;
	    }
	    delta = angle - gangle;
	    aave += abs(delta);
	    arms += delta * delta;
	    if (*prtflg == 1) {
		io___172.ciunit = iounit_1.iout;
		s_wsfe(&io___172);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&gangle, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&angle, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&delta, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
	if (angle_1.nangle != 0) {
	    aave /= (doublereal) angle_1.nangle;
	}
	if (angle_1.nangle != 0) {
	    arms = sqrt(arms / (doublereal) angle_1.nangle);
	}
	if (*prtflg == 1 && angle_1.nangle != 0) {
	    io___173.ciunit = iounit_1.iout;
	    s_wsfe(&io___173);
	    do_fio(&c__1, (char *)&aave, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&arms, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     compute the RMS between QM and TINKER torsion angles */

    tave = 0.;
    trms = 0.;
    if (valfit_1.fit_struct__) {
	if (*prtflg == 1 && tors_1.ntors != 0) {
	    io___176.ciunit = iounit_1.iout;
	    s_wsfe(&io___176);
	    e_wsfe();
	}
	i__1 = tors_1.ntors;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = itors_ref(1, i__);
	    ib = itors_ref(2, i__);
	    ic = itors_ref(3, i__);
	    id = itors_ref(4, i__);
	    xba = atoms_1.x[ib - 1] - atoms_1.x[ia - 1];
	    yba = atoms_1.y[ib - 1] - atoms_1.y[ia - 1];
	    zba = atoms_1.z__[ib - 1] - atoms_1.z__[ia - 1];
	    xcb = atoms_1.x[ic - 1] - atoms_1.x[ib - 1];
	    ycb = atoms_1.y[ic - 1] - atoms_1.y[ib - 1];
	    zcb = atoms_1.z__[ic - 1] - atoms_1.z__[ib - 1];
	    xdc = atoms_1.x[id - 1] - atoms_1.x[ic - 1];
	    ydc = atoms_1.y[id - 1] - atoms_1.y[ic - 1];
	    zdc = atoms_1.z__[id - 1] - atoms_1.z__[ic - 1];
	    xt = yba * zcb - ycb * zba;
	    yt = zba * xcb - zcb * xba;
	    zt = xba * ycb - xcb * yba;
	    xu = ycb * zdc - ydc * zcb;
	    yu = zcb * xdc - zdc * xcb;
	    zu = xcb * ydc - xdc * ycb;
	    xtu = yt * zu - yu * zt;
	    ytu = zt * xu - zu * xt;
	    ztu = xt * yu - xu * yt;
	    rt2 = xt * xt + yt * yt + zt * zt;
	    ru2 = xu * xu + yu * yu + zu * zu;
	    rtru = sqrt(rt2 * ru2);
	    if (rtru != 0.) {
		rcb = sqrt(xcb * xcb + ycb * ycb + zcb * zcb);
		cosine = (xt * xu + yt * yu + zt * zu) / rtru;
		sine = (xcb * xtu + ycb * ytu + zcb * ztu) / (rcb * rtru);
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine);
		cosine = min(d__1,d__2);
		angle = acos(cosine) * 57.29577951308232088;
		if (sine < 0.) {
		    angle = -angle;
		}
	    }
	    xba = qmstuf_1.gx[ib - 1] - qmstuf_1.gx[ia - 1];
	    yba = qmstuf_1.gy[ib - 1] - qmstuf_1.gy[ia - 1];
	    zba = qmstuf_1.gz[ib - 1] - qmstuf_1.gz[ia - 1];
	    xcb = qmstuf_1.gx[ic - 1] - qmstuf_1.gx[ib - 1];
	    ycb = qmstuf_1.gy[ic - 1] - qmstuf_1.gy[ib - 1];
	    zcb = qmstuf_1.gz[ic - 1] - qmstuf_1.gz[ib - 1];
	    xdc = qmstuf_1.gx[id - 1] - qmstuf_1.gx[ic - 1];
	    ydc = qmstuf_1.gy[id - 1] - qmstuf_1.gy[ic - 1];
	    zdc = qmstuf_1.gz[id - 1] - qmstuf_1.gz[ic - 1];
	    xt = yba * zcb - ycb * zba;
	    yt = zba * xcb - zcb * xba;
	    zt = xba * ycb - xcb * yba;
	    xu = ycb * zdc - ydc * zcb;
	    yu = zcb * xdc - zdc * xcb;
	    zu = xcb * ydc - xdc * ycb;
	    xtu = yt * zu - yu * zt;
	    ytu = zt * xu - zu * xt;
	    ztu = xt * yu - xu * yt;
	    rt2 = xt * xt + yt * yt + zt * zt;
	    ru2 = xu * xu + yu * yu + zu * zu;
	    rtru = sqrt(rt2 * ru2);
	    if (rtru != 0.) {
		rcb = sqrt(xcb * xcb + ycb * ycb + zcb * zcb);
		cosine = (xt * xu + yt * yu + zt * zu) / rtru;
		sine = (xcb * xtu + ycb * ytu + zcb * ztu) / (rcb * rtru);
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine);
		cosine = min(d__1,d__2);
		gangle = acos(cosine) * 57.29577951308232088;
		if (sine < 0.) {
		    gangle = -gangle;
		}
	    }
	    delta = angle - gangle;
	    if (delta > 180.) {
		delta += -360.;
	    }
	    if (delta < -180.) {
		delta += 360.;
	    }
	    tave += abs(delta);
	    trms += delta * delta;
	    if (*prtflg == 1) {
		io___198.ciunit = iounit_1.iout;
		s_wsfe(&io___198);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&gangle, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&angle, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&delta, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
	if (tors_1.ntors != 0) {
	    tave /= (doublereal) tors_1.ntors;
	}
	if (tors_1.ntors != 0) {
	    trms = sqrt(trms / (doublereal) tors_1.ntors);
	}
	if (*prtflg == 1 && tors_1.ntors != 0) {
	    io___199.ciunit = iounit_1.iout;
	    s_wsfe(&io___199);
	    do_fio(&c__1, (char *)&tave, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&trms, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     compute the RMS between QM and TINKER gradient components */

    gave = 0.;
    grms = 0.;
    if (valfit_1.fit_force__) {
	gradient_(&energy, derivs);
	if (*prtflg == 1) {
	    io___204.ciunit = iounit_1.iout;
	    s_wsfe(&io___204);
	    e_wsfe();
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		delta = gforce_ref(j, i__) - derivs_ref(j, i__);
		gave += abs(delta);
		grms += delta * delta;
		if (*prtflg == 1) {
		    io___206.ciunit = iounit_1.iout;
		    s_wsfe(&io___206);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, axis + (j - 1), (ftnlen)1);
		    do_fio(&c__1, (char *)&gforce_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&derivs_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&delta, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	}
	gave /= (doublereal) (atoms_1.n * 3);
	grms = sqrt(grms / (doublereal) (atoms_1.n * 3));
	if (*prtflg == 1) {
	    io___207.ciunit = iounit_1.iout;
	    s_wsfe(&io___207);
	    do_fio(&c__1, (char *)&gave, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&grms, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     calculate the full Hessian matrix of second derivatives */

    hescut_1.hesscut = 0.;
    hessian_(h__, hinit, hstop, hindex, hdiag);

/*     compute the RMS between QM and TINKER Hessian elements */

    have = 0.;
    hrms = 0.;
    if (valfit_1.fit_force__) {
	if (*prtflg == 1) {
	    io___215.ciunit = iounit_1.iout;
	    s_wsfe(&io___215);
	    e_wsfe();
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		m1 = (i__ - 1) * 3 + j;
		m = m1 * (m1 + 1) / 2;
		delta = qmstuf_1.gh[m - 1] - hdiag_ref(j, i__);
		have += abs(delta);
		hrms += delta * delta;
		if (*prtflg == 1) {
		    s_wsfe(&io___218);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, axis + (j - 1), (ftnlen)1);
		    do_fio(&c__1, (char *)&qmstuf_1.gh[m - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&hdiag_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&delta, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
		m1 = (i__ - 1) * 3 + j;
		m2 = m1;
		i__2 = hstop_ref(j, i__);
		for (k = hinit_ref(j, i__); k <= i__2; ++k) {
		    ++m2;
		    m = m1 + m2 * (m2 - 1) / 2;
		    delta = qmstuf_1.gh[m - 1] - h__[k - 1];
		    have += abs(delta);
		    hrms += delta * delta;
		}
	    }
	}
	have /= (doublereal) ((atoms_1.n * 9 * atoms_1.n + atoms_1.n * 3) / 2)
		;
	hrms = sqrt(hrms / (doublereal) ((atoms_1.n * 9 * atoms_1.n + 
		atoms_1.n * 3) / 2));
	if (*prtflg == 1) {
	    io___221.ciunit = iounit_1.iout;
	    s_wsfe(&io___221);
	    do_fio(&c__1, (char *)&have, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&hrms, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     set atomic mass roots needed for vibrational analysis */

    nfreq = atoms_1.n * 3;
    if (nfreq > 1000 || nfreq * nfreq > 1000000) {
	io___223.ciunit = iounit_1.iout;
	s_wsfe(&io___223);
	e_wsfe();
	fatal_();
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mass2[i__ - 1] = sqrt(atmtyp_1.mass[i__ - 1]);
    }

/*     store upper triangle of the mass-weighted Hessian matrix */

    ihess = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    ++ihess;
	    matrix[ihess - 1] = hdiag_ref(j, i__) / atmtyp_1.mass[i__ - 1];
	    i__2 = hstop_ref(j, i__);
	    for (k = hinit_ref(j, i__); k <= i__2; ++k) {
		m = (hindex[k - 1] + 2) / 3;
		++ihess;
		matrix[ihess - 1] = h__[k - 1] / (mass2[i__ - 1] * mass2[m - 
			1]);
	    }
	}
    }

/*     diagonalize to get vibrational frequencies and normal modes */

    diagq_(&nfreq, &c__1000, &nfreq, matrix, eigen, vects, a, b, p, w, ta, tb,
	     ty);
    factor = sqrt(418.4) / .1883651567308853;
    i__1 = nfreq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	eigen[i__ - 1] = factor * d_sign(&c_b215, &eigen[i__ - 1]) * sqrt((
		d__1 = eigen[i__ - 1], abs(d__1)));
    }

/*     compute the RMS between QM and TINKER vibrational frequencies */

    fcut = 800.;
    if (valfit_1.fit_tors__) {
	fcut = 200.;
    }
    fave = 0.;
    frms = 0.;
    if (*prtflg == 1) {
	io___240.ciunit = iounit_1.iout;
	s_wsfe(&io___240);
	e_wsfe();
    }
    k = 0;
    for (i__ = nfreq; i__ >= 7; --i__) {
	if (qmstuf_1.gfreq[i__ - 7] > fcut) {
	    ++k;
	    delta = eigen[i__ - 1] - qmstuf_1.gfreq[i__ - 7];
	    fave += abs(delta);
	    frms += delta * delta;
	    if (*prtflg == 1) {
		io___241.ciunit = iounit_1.iout;
		s_wsfe(&io___241);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&qmstuf_1.gfreq[i__ - 7], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&eigen[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&delta, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
    }
    fave /= (doublereal) k;
    frms = sqrt(frms / (doublereal) k);
    if (*prtflg == 1) {
	io___242.ciunit = iounit_1.iout;
	s_wsfe(&io___242);
	do_fio(&c__1, (char *)&fave, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&frms, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     sum weighted RMS values to get overall error function */

    bfac = 100.;
    afac = 10.;
    tfac = 1.;
    gfac = 10.;
    hfac = .1;
    ffac = .1;
    ret_val = bfac * brms + afac * arms + tfac * trms + gfac * grms + hfac * 
	    hrms + ffac * frms;
    return ret_val;
} /* valrms_ */

#undef derivs_ref
#undef gforce_ref
#undef itors_ref
#undef hstop_ref
#undef hinit_ref
#undef hdiag_ref
#undef iang_ref
#undef ibnd_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  function minimiz1  --  energy and gradient for minimize  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "minimiz1" is a service routine that computes the energy and */
/*     gradient for a low storage BFGS optimization in Cartesian */
/*     coordinate space */


doublereal minimiz1_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static logical analytic;
    static doublereal e;
    static integer i__;
    static doublereal eps;
    static integer nvar;
    extern doublereal energy_(void);
    static doublereal derivs[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int numgrad_(D_fp, doublereal *, doublereal *);


#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1 - 4]



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
/*     ##  scales.i  --  parameter scale factors for optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     scale      multiplicative factor for each optimization parameter */
/*     set_scale  logical flag to show if scale factors have been set */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     use either analytical or numerical gradients */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    analytic = TRUE_;
    eps = 1e-5;

/*     translate optimization parameters to atomic coordinates */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    atoms_1.x[i__ - 1] = xx[nvar] / scales_1.scale[nvar - 1];
	    ++nvar;
	    atoms_1.y[i__ - 1] = xx[nvar] / scales_1.scale[nvar - 1];
	    ++nvar;
	    atoms_1.z__[i__ - 1] = xx[nvar] / scales_1.scale[nvar - 1];
	}
    }

/*     compute and store the energy and gradient */

    if (analytic) {
	gradient_(&e, derivs);
    } else {
	e = energy_();
	numgrad_((D_fp)energy_, derivs, &eps);
    }
    ret_val = e;

/*     store Cartesian gradient as optimization gradient */

    nvar = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++nvar;
	    g[nvar] = derivs_ref(1, i__) / scales_1.scale[nvar - 1];
	    ++nvar;
	    g[nvar] = derivs_ref(2, i__) / scales_1.scale[nvar - 1];
	    ++nvar;
	    g[nvar] = derivs_ref(3, i__) / scales_1.scale[nvar - 1];
	}
    }
    return ret_val;
} /* minimiz1_ */

#undef derivs_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine prmvar  --  valence terms to optimization  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "prmvar" determines the optimization values from the */
/*     corresponding valence potential energy parameters */


/* Subroutine */ int prmvar_(integer *nvar, doublereal *xx)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Valence Parameters Used in Structure Fit"
	    "ting :\002)";
    static char fmt_20[] = "(/,\002 Valence Parameters Used in Force Fitting"
	    " :\002)";
    static char fmt_30[] = "(/,\002 Parameter\002,10x,\002Atom Classes\002,1"
	    "0x,\002Category\002,12x,\002Value\002,5x,\002Fixed\002,/)";
    static char fmt_40[] = "(i6,5x,2i6,19x,a10,3x,f12.4)";
    static char fmt_50[] = "(i6,5x,2i6,19x,a11,2x,f12.4)";
    static char fmt_60[] = "(4x,\002--\002,5x,2i6,19x,a10,3x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_70[] = "(4x,\002--\002,5x,2i6,19x,a11,2x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_80[] = "(i6,5x,3i6,13x,a11,2x,f12.4)";
    static char fmt_90[] = "(i6,5x,3i6,13x,a11,2x,f12.4)";
    static char fmt_100[] = "(4x,\002--\002,5x,3i6,13x,a11,2x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_110[] = "(4x,\002--\002,5x,3i6,13x,a11,2x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_120[] = "(i6,5x,3i6,13x,a8,5x,f12.4)";
    static char fmt_130[] = "(i6,5x,3i6,13x,a8,5x,f12.4)";
    static char fmt_140[] = "(4x,\002--\002,5x,3i6,13x,a8,5x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_150[] = "(4x,\002--\002,5x,3i6,13x,a8,5x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_160[] = "(i6,5x,3i6,13x,a10,3x,f12.4)";
    static char fmt_170[] = "(i6,5x,3i6,13x,a9,4x,f12.4)";
    static char fmt_180[] = "(4x,\002--\002,5x,3i6,13x,a10,3x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_190[] = "(4x,\002--\002,5x,3i6,13x,a9,4x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_200[] = "(i6,5x,4i6,7x,a8,5x,f12.4)";
    static char fmt_210[] = "(4x,\002--\002,5x,4i6,7x,a8,5x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_220[] = "(i6,5x,4i6,7x,a9,4x,f12.4)";
    static char fmt_230[] = "(4x,\002--\002,5x,4i6,7x,a9,4x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_240[] = "(i6,5x,4i6,7x,a9,4x,f12.4)";
    static char fmt_250[] = "(4x,\002--\002,5x,4i6,7x,a9,4x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_260[] = "(i6,5x,4i6,7x,a9,4x,f12.4)";
    static char fmt_270[] = "(4x,\002--\002,5x,4i6,7x,a9,4x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_280[] = "(i6,5x,4i6,7x,a9,4x,f12.4)";
    static char fmt_290[] = "(4x,\002--\002,5x,4i6,7x,a9,4x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_300[] = "(i6,5x,4i6,7x,a9,4x,f12.4)";
    static char fmt_310[] = "(4x,\002--\002,5x,4i6,7x,a9,4x,f12.4,7x,\002"
	    "X\002)";
    static char fmt_320[] = "(i6,5x,4i6,7x,a9,4x,f12.4)";
    static char fmt_330[] = "(4x,\002--\002,5x,4i6,7x,a9,4x,f12.4,7x,\002"
	    "X\002)";

    /* System generated locals */
    address a__1[2], a__2[3], a__3[4];
    integer i__1, i__2[2], i__3, i__4[3], i__5[4], i__6;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), do_fio(integer *, char *, 
	    ftnlen);

    /* Local variables */
    static integer i__, k, ia, ib, ic, id, ka, ii, kb, kc, kd, kk;
    static char pa[4], pb[4], pc[4], pd[4];
    static integer ita, itb, itc, itd, kta, ktb, ktc, ktd;
    static logical done;
    static char pita[12], pitb[8], pkta[12], pktb[8];
    static integer size;
    static char pitt[16], pktt[16];
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    ;

    /* Fortran I/O blocks */
    static cilist io___255 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___256 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___257 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___274 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___275 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___276 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___277 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___285 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___286 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___287 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___288 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___291 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___292 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___293 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___294 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___295 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___296 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___297 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___298 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___306 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___307 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___308 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___309 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___310 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___311 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___312 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___313 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___314 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___315 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___316 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___317 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___318 = { 0, 0, 0, fmt_320, 0 };
    static cilist io___319 = { 0, 0, 0, fmt_330, 0 };



#define isb_ref(a_1,a_2) strbnd_1.isb[(a_2)*3 + a_1 - 4]
#define sbk_ref(a_1,a_2) strbnd_1.sbk[(a_2)*2 + a_1 - 3]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define iury_ref(a_1,a_2) urey_1.iury[(a_2)*3 + a_1 - 4]
#define tors1_ref(a_1,a_2) tors_1.tors1[(a_2)*4 + a_1 - 5]
#define tors2_ref(a_1,a_2) tors_1.tors2[(a_2)*4 + a_1 - 5]
#define tors3_ref(a_1,a_2) tors_1.tors3[(a_2)*4 + a_1 - 5]
#define tors4_ref(a_1,a_2) tors_1.tors4[(a_2)*4 + a_1 - 5]
#define tors5_ref(a_1,a_2) tors_1.tors5[(a_2)*4 + a_1 - 5]
#define tors6_ref(a_1,a_2) tors_1.tors6[(a_2)*4 + a_1 - 5]
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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  opbend.i  --  out-of-plane bends in the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     opbk      force constant values for out-of-plane bending */
/*     nopbend   total number of out-of-plane bends in the system */
/*     iopb      bond angle numbers used in out-of-plane bending */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  strbnd.i  --  stretch-bends in the current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     sbk       force constants for stretch-bend terms */
/*     nstrbnd   total number of stretch-bend interactions */
/*     isb       angle and bond numbers used in stretch-bend */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  urey.i  --  Urey-Bradley interactions in the structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     uk      Urey-Bradley force constants (kcal/mole/Ang**2) */
/*     ul      ideal 1-3 distance values in Angstroms */
/*     nurey   total number of Urey-Bradley terms in the system */
/*     iury    numbers of the atoms in each Urey-Bradley interaction */




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  valfit.i  --  values for valence term parameter fitting  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     fit_bond    logical flag to fit bond stretch parameters */
/*     fit_angle   logical flag to fit angle bend parameters */
/*     fit_strbnd  logical flag to fit stretch-bend parameters */
/*     fit_urey    logical flag to fit Urey-Bradley parameters */
/*     fit_opbend  logical flag to fit out-of-plane bend parameters */
/*     fit_tors    logical flag to fit torsional parameters */
/*     fit_struct  logical flag to structure-fit valence parameters */
/*     fit_force   logical flag to force-fit valence parameters */




/*     zero out the total number of optimization parameters */

    /* Parameter adjustments */
    --xx;

    /* Function Body */
    *nvar = 0;

/*     print a header for the parameters used in fitting */

    if (valfit_1.fit_struct__) {
	io___255.ciunit = iounit_1.iout;
	s_wsfe(&io___255);
	e_wsfe();
    } else if (valfit_1.fit_force__) {
	io___256.ciunit = iounit_1.iout;
	s_wsfe(&io___256);
	e_wsfe();
    }
    io___257.ciunit = iounit_1.iout;
    s_wsfe(&io___257);
    e_wsfe();

/*     find bond stretch force constants and target lengths */

    i__1 = bond_1.nbond;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	ia = ibnd_ref(1, i__);
	ib = ibnd_ref(2, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	if (ita <= itb) {
/* Writing concatenation */
	    i__2[0] = 4, a__1[0] = pa;
	    i__2[1] = 4, a__1[1] = pb;
	    s_cat(pitb, a__1, i__2, &c__2, (ftnlen)8);
	} else {
/* Writing concatenation */
	    i__2[0] = 4, a__1[0] = pb;
	    i__2[1] = 4, a__1[1] = pa;
	    s_cat(pitb, a__1, i__2, &c__2, (ftnlen)8);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    ka = ibnd_ref(1, k);
	    kb = ibnd_ref(2, k);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    if (kta <= ktb) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(pktb, a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(pktb, a__1, i__2, &c__2, (ftnlen)8);
	    }
	    if (s_cmp(pktb, pitb, (ftnlen)8, (ftnlen)8) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_bond__ && bond_1.bk[i__ - 1] != 0.) {
		++(*nvar);
		xx[*nvar] = bond_1.bk[i__ - 1];
		io___274.ciunit = iounit_1.iout;
		s_wsfe(&io___274);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Bond Force", (ftnlen)10);
		do_fio(&c__1, (char *)&bond_1.bk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
		++(*nvar);
		xx[*nvar] = bond_1.bl[i__ - 1];
		xx[*nvar] *= 100.;
		io___275.ciunit = iounit_1.iout;
		s_wsfe(&io___275);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Bond Length", (ftnlen)11);
		do_fio(&c__1, (char *)&bond_1.bl[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___276.ciunit = iounit_1.iout;
		s_wsfe(&io___276);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Bond Force", (ftnlen)10);
		do_fio(&c__1, (char *)&bond_1.bk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
		io___277.ciunit = iounit_1.iout;
		s_wsfe(&io___277);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Bond Length", (ftnlen)11);
		do_fio(&c__1, (char *)&bond_1.bl[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     find angle bend force constants and target angles */

    i__1 = angle_1.nangle;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pa;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pc;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	} else {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pc;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pa;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    ka = iang_ref(1, k);
	    kb = iang_ref(2, k);
	    kc = iang_ref(3, k);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    if (kta <= ktc) {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pa;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pc;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pc;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pa;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    }
	    if (s_cmp(pkta, pita, (ftnlen)12, (ftnlen)12) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_angle__ && angle_1.ak[i__ - 1] != 0.) {
		++(*nvar);
		xx[*nvar] = angle_1.ak[i__ - 1];
		io___285.ciunit = iounit_1.iout;
		s_wsfe(&io___285);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Angle Force", (ftnlen)11);
		do_fio(&c__1, (char *)&angle_1.ak[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
		++(*nvar);
		xx[*nvar] = angle_1.anat[i__ - 1];
		io___286.ciunit = iounit_1.iout;
		s_wsfe(&io___286);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Angle Value", (ftnlen)11);
		do_fio(&c__1, (char *)&angle_1.anat[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___287.ciunit = iounit_1.iout;
		s_wsfe(&io___287);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Angle Force", (ftnlen)11);
		do_fio(&c__1, (char *)&angle_1.ak[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
		io___288.ciunit = iounit_1.iout;
		s_wsfe(&io___288);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Angle Value", (ftnlen)11);
		do_fio(&c__1, (char *)&angle_1.anat[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     find stretch-bend force constant parameter values */

    i__1 = strbnd_1.nstrbnd;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	ii = isb_ref(1, i__);
	ia = iang_ref(1, ii);
	ib = iang_ref(2, ii);
	ic = iang_ref(3, ii);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pa;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pc;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	} else {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pc;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pa;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    kk = isb_ref(1, k);
	    ka = iang_ref(1, kk);
	    kb = iang_ref(2, kk);
	    kc = iang_ref(3, kk);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    if (kta <= ktc) {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pa;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pc;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pc;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pa;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    }
	    if (s_cmp(pkta, pita, (ftnlen)12, (ftnlen)12) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_strbnd__ && sbk_ref(1, i__) != 0. && sbk_ref(2, 
		    i__) != 0.) {
		++(*nvar);
		xx[*nvar] = sbk_ref(1, i__);
		io___291.ciunit = iounit_1.iout;
		s_wsfe(&io___291);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, "StrBnd-1", (ftnlen)8);
		do_fio(&c__1, (char *)&sbk_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
		++(*nvar);
		xx[*nvar] = sbk_ref(2, i__);
		io___292.ciunit = iounit_1.iout;
		s_wsfe(&io___292);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, "StrBnd-2", (ftnlen)8);
		do_fio(&c__1, (char *)&sbk_ref(2, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___293.ciunit = iounit_1.iout;
		s_wsfe(&io___293);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, "StrBnd-1", (ftnlen)8);
		do_fio(&c__1, (char *)&sbk_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
		io___294.ciunit = iounit_1.iout;
		s_wsfe(&io___294);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, "StrBnd-2", (ftnlen)8);
		do_fio(&c__1, (char *)&sbk_ref(2, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     find Urey-Bradley force constant parameter values */

    i__1 = urey_1.nurey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	ia = iury_ref(1, i__);
	ib = iury_ref(2, i__);
	ic = iury_ref(3, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pa;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pc;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	} else {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pc;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pa;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    ka = iury_ref(1, k);
	    kb = iury_ref(2, k);
	    kc = iury_ref(3, k);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    if (kta <= ktc) {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pa;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pc;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pc;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pa;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    }
	    if (s_cmp(pkta, pita, (ftnlen)12, (ftnlen)12) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_urey__ && urey_1.uk[i__ - 1] != 0.) {
		++(*nvar);
		xx[*nvar] = urey_1.uk[i__ - 1];
		io___295.ciunit = iounit_1.iout;
		s_wsfe(&io___295);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Urey Force", (ftnlen)10);
		do_fio(&c__1, (char *)&urey_1.uk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
		++(*nvar);
		xx[*nvar] = urey_1.ul[i__ - 1];
		io___296.ciunit = iounit_1.iout;
		s_wsfe(&io___296);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Urey Dist", (ftnlen)9);
		do_fio(&c__1, (char *)&urey_1.ul[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___297.ciunit = iounit_1.iout;
		s_wsfe(&io___297);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Urey Force", (ftnlen)10);
		do_fio(&c__1, (char *)&urey_1.uk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
		io___298.ciunit = iounit_1.iout;
		s_wsfe(&io___298);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Urey Dist", (ftnlen)9);
		do_fio(&c__1, (char *)&urey_1.ul[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     find out-of-plane bend force constant parameter values */

    i__1 = opbend_1.nopbend;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	ii = opbend_1.iopb[i__ - 1];
	ia = iang_ref(1, ii);
	ib = iang_ref(2, ii);
	ic = iang_ref(3, ii);
	id = iang_ref(4, ii);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	itd = atmtyp_1.class__[id - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	numeral_(&itd, pd, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pd;
	    i__5[1] = 4, a__3[1] = pb;
	    i__5[2] = 4, a__3[2] = pa;
	    i__5[3] = 4, a__3[3] = pc;
	    s_cat(pitt, a__3, i__5, &c__4, (ftnlen)16);
	} else {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pd;
	    i__5[1] = 4, a__3[1] = pb;
	    i__5[2] = 4, a__3[2] = pc;
	    i__5[3] = 4, a__3[3] = pa;
	    s_cat(pitt, a__3, i__5, &c__4, (ftnlen)16);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    kk = opbend_1.iopb[k - 1];
	    ka = iang_ref(1, kk);
	    kb = iang_ref(2, kk);
	    kc = iang_ref(3, kk);
	    kd = iang_ref(4, kk);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    ktd = atmtyp_1.class__[kd - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    numeral_(&ktd, pd, &size, (ftnlen)4);
	    if (kta <= ktc) {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pd;
		i__5[1] = 4, a__3[1] = pb;
		i__5[2] = 4, a__3[2] = pa;
		i__5[3] = 4, a__3[3] = pc;
		s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
	    } else {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pd;
		i__5[1] = 4, a__3[1] = pb;
		i__5[2] = 4, a__3[2] = pc;
		i__5[3] = 4, a__3[3] = pa;
		s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
	    }
	    if (s_cmp(pktt, pitt, (ftnlen)16, (ftnlen)16) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_opbend__ && opbend_1.opbk[i__ - 1] != 0.) {
		++(*nvar);
		xx[*nvar] = opbend_1.opbk[i__ - 1];
		io___306.ciunit = iounit_1.iout;
		s_wsfe(&io___306);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		i__3 = min(ita,itc);
		do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		i__6 = max(ita,itc);
		do_fio(&c__1, (char *)&i__6, (ftnlen)sizeof(integer));
		do_fio(&c__1, "O-P-Bend", (ftnlen)8);
		do_fio(&c__1, (char *)&opbend_1.opbk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___307.ciunit = iounit_1.iout;
		s_wsfe(&io___307);
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		i__3 = min(ita,itc);
		do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		i__6 = max(ita,itc);
		do_fio(&c__1, (char *)&i__6, (ftnlen)sizeof(integer));
		do_fio(&c__1, "O-P-Bend", (ftnlen)8);
		do_fio(&c__1, (char *)&opbend_1.opbk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     find torsional angle amplitude parameter values */

    i__1 = tors_1.ntors;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	ia = itors_ref(1, i__);
	ib = itors_ref(2, i__);
	ic = itors_ref(3, i__);
	id = itors_ref(4, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	itd = atmtyp_1.class__[id - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	numeral_(&itd, pd, &size, (ftnlen)4);
	if (itb < itc || itb == itc && ita <= itd) {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pa;
	    i__5[1] = 4, a__3[1] = pb;
	    i__5[2] = 4, a__3[2] = pc;
	    i__5[3] = 4, a__3[3] = pd;
	    s_cat(pitt, a__3, i__5, &c__4, (ftnlen)16);
	} else {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pd;
	    i__5[1] = 4, a__3[1] = pc;
	    i__5[2] = 4, a__3[2] = pb;
	    i__5[3] = 4, a__3[3] = pa;
	    s_cat(pitt, a__3, i__5, &c__4, (ftnlen)16);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    ka = itors_ref(1, k);
	    kb = itors_ref(2, k);
	    kc = itors_ref(3, k);
	    kd = itors_ref(4, k);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    ktd = atmtyp_1.class__[kd - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    numeral_(&ktd, pd, &size, (ftnlen)4);
	    if (ktb < ktc || ktb == ktc && kta <= ktd) {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pa;
		i__5[1] = 4, a__3[1] = pb;
		i__5[2] = 4, a__3[2] = pc;
		i__5[3] = 4, a__3[3] = pd;
		s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
	    } else {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pd;
		i__5[1] = 4, a__3[1] = pc;
		i__5[2] = 4, a__3[2] = pb;
		i__5[3] = 4, a__3[3] = pa;
		s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
	    }
	    if (s_cmp(pktt, pitt, (ftnlen)16, (ftnlen)16) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_tors__ && tors1_ref(1, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = tors1_ref(1, i__);
		io___308.ciunit = iounit_1.iout;
		s_wsfe(&io___308);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Torsion-1", (ftnlen)9);
		do_fio(&c__1, (char *)&tors1_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___309.ciunit = iounit_1.iout;
		s_wsfe(&io___309);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Torsion-1", (ftnlen)9);
		do_fio(&c__1, (char *)&tors1_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	    if (valfit_1.fit_tors__ && tors2_ref(1, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = tors2_ref(1, i__);
		io___310.ciunit = iounit_1.iout;
		s_wsfe(&io___310);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Torsion-2", (ftnlen)9);
		do_fio(&c__1, (char *)&tors2_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___311.ciunit = iounit_1.iout;
		s_wsfe(&io___311);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Torsion-2", (ftnlen)9);
		do_fio(&c__1, (char *)&tors2_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	    if (valfit_1.fit_tors__ && tors3_ref(1, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = tors3_ref(1, i__);
		io___312.ciunit = iounit_1.iout;
		s_wsfe(&io___312);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Torsion-3", (ftnlen)9);
		do_fio(&c__1, (char *)&tors3_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___313.ciunit = iounit_1.iout;
		s_wsfe(&io___313);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Torsion-3", (ftnlen)9);
		do_fio(&c__1, (char *)&tors3_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	    if (valfit_1.fit_tors__ && tors4_ref(1, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = tors4_ref(1, i__);
		io___314.ciunit = iounit_1.iout;
		s_wsfe(&io___314);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Torsion-4", (ftnlen)9);
		do_fio(&c__1, (char *)&tors4_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else if (tors4_ref(1, i__) != 0.) {
		io___315.ciunit = iounit_1.iout;
		s_wsfe(&io___315);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Torsion-4", (ftnlen)9);
		do_fio(&c__1, (char *)&tors4_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	    if (valfit_1.fit_tors__ && tors5_ref(1, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = tors5_ref(1, i__);
		io___316.ciunit = iounit_1.iout;
		s_wsfe(&io___316);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Torsion-5", (ftnlen)9);
		do_fio(&c__1, (char *)&tors5_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else if (tors5_ref(1, i__) != 0.) {
		io___317.ciunit = iounit_1.iout;
		s_wsfe(&io___317);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Torsion-5", (ftnlen)9);
		do_fio(&c__1, (char *)&tors5_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	    if (valfit_1.fit_tors__ && tors6_ref(1, i__) != 0.) {
		++(*nvar);
		xx[*nvar] = tors6_ref(1, i__);
		io___318.ciunit = iounit_1.iout;
		s_wsfe(&io___318);
		do_fio(&c__1, (char *)&(*nvar), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Torsion-6", (ftnlen)9);
		do_fio(&c__1, (char *)&tors6_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else if (tors6_ref(1, i__) != 0.) {
		io___319.ciunit = iounit_1.iout;
		s_wsfe(&io___319);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, "Torsion-6", (ftnlen)9);
		do_fio(&c__1, (char *)&tors6_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }
    return 0;
} /* prmvar_ */

#undef itors_ref
#undef tors6_ref
#undef tors5_ref
#undef tors4_ref
#undef tors3_ref
#undef tors2_ref
#undef tors1_ref
#undef iury_ref
#undef iang_ref
#undef ibnd_ref
#undef sbk_ref
#undef isb_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine varprm  --  optimization to valence terms  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "varprm" copies the current optimization values into the */
/*     corresponding valence potential energy parameters */


/* Subroutine */ int varprm_(integer *nvar, doublereal *xx, integer *ivar, 
	doublereal *eps)
{
    /* System generated locals */
    address a__1[2], a__2[3], a__3[4];
    integer i__1, i__2[2], i__3, i__4[3], i__5[4];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, k, ia, ib, ic, id, ka, ii, kb, kc, kd, kk;
    static char pa[4], pb[4], pc[4], pd[4];
    static integer ita, itb, itc, itd, kta, ktb, ktc, ktd;
    static logical done;
    static char pita[12], pitb[8], pkta[12], pktb[8];
    static integer size;
    static char pitt[16], pktt[16];
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    ;


#define isb_ref(a_1,a_2) strbnd_1.isb[(a_2)*3 + a_1 - 4]
#define sbk_ref(a_1,a_2) strbnd_1.sbk[(a_2)*2 + a_1 - 3]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define iury_ref(a_1,a_2) urey_1.iury[(a_2)*3 + a_1 - 4]
#define tors1_ref(a_1,a_2) tors_1.tors1[(a_2)*4 + a_1 - 5]
#define tors2_ref(a_1,a_2) tors_1.tors2[(a_2)*4 + a_1 - 5]
#define tors3_ref(a_1,a_2) tors_1.tors3[(a_2)*4 + a_1 - 5]
#define tors4_ref(a_1,a_2) tors_1.tors4[(a_2)*4 + a_1 - 5]
#define tors5_ref(a_1,a_2) tors_1.tors5[(a_2)*4 + a_1 - 5]
#define tors6_ref(a_1,a_2) tors_1.tors6[(a_2)*4 + a_1 - 5]
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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  opbend.i  --  out-of-plane bends in the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     opbk      force constant values for out-of-plane bending */
/*     nopbend   total number of out-of-plane bends in the system */
/*     iopb      bond angle numbers used in out-of-plane bending */




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
/*     ##  strbnd.i  --  stretch-bends in the current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     sbk       force constants for stretch-bend terms */
/*     nstrbnd   total number of stretch-bend interactions */
/*     isb       angle and bond numbers used in stretch-bend */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  urey.i  --  Urey-Bradley interactions in the structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     uk      Urey-Bradley force constants (kcal/mole/Ang**2) */
/*     ul      ideal 1-3 distance values in Angstroms */
/*     nurey   total number of Urey-Bradley terms in the system */
/*     iury    numbers of the atoms in each Urey-Bradley interaction */




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  valfit.i  --  values for valence term parameter fitting  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     fit_bond    logical flag to fit bond stretch parameters */
/*     fit_angle   logical flag to fit angle bend parameters */
/*     fit_strbnd  logical flag to fit stretch-bend parameters */
/*     fit_urey    logical flag to fit Urey-Bradley parameters */
/*     fit_opbend  logical flag to fit out-of-plane bend parameters */
/*     fit_tors    logical flag to fit torsional parameters */
/*     fit_struct  logical flag to structure-fit valence parameters */
/*     fit_force   logical flag to force-fit valence parameters */




/*     zero out the total number of optimization parameters */

    /* Parameter adjustments */
    --xx;

    /* Function Body */
    *nvar = 0;

/*     translate optimization values to bond stretch parameters */

    i__1 = bond_1.nbond;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	ia = ibnd_ref(1, i__);
	ib = ibnd_ref(2, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	if (ita <= itb) {
/* Writing concatenation */
	    i__2[0] = 4, a__1[0] = pa;
	    i__2[1] = 4, a__1[1] = pb;
	    s_cat(pitb, a__1, i__2, &c__2, (ftnlen)8);
	} else {
/* Writing concatenation */
	    i__2[0] = 4, a__1[0] = pb;
	    i__2[1] = 4, a__1[1] = pa;
	    s_cat(pitb, a__1, i__2, &c__2, (ftnlen)8);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    ka = ibnd_ref(1, k);
	    kb = ibnd_ref(2, k);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    if (kta <= ktb) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(pktb, a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(pktb, a__1, i__2, &c__2, (ftnlen)8);
	    }
	    if (s_cmp(pktb, pitb, (ftnlen)8, (ftnlen)8) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_bond__ && bond_1.bk[i__ - 1] != 0.) {
		++(*nvar);
		bond_1.bk[i__ - 1] = xx[*nvar];
		if (*ivar == *nvar) {
		    bond_1.bk[i__ - 1] += *eps;
		}
		++(*nvar);
		bond_1.bl[i__ - 1] = xx[*nvar];
		if (*ivar == *nvar) {
		    bond_1.bl[i__ - 1] += *eps;
		}
		bond_1.bl[i__ - 1] *= .01;
		i__3 = bond_1.nbond;
		for (k = i__ + 1; k <= i__3; ++k) {
		    ka = ibnd_ref(1, k);
		    kb = ibnd_ref(2, k);
		    kta = atmtyp_1.class__[ka - 1];
		    ktb = atmtyp_1.class__[kb - 1];
		    size = 4;
		    numeral_(&kta, pa, &size, (ftnlen)4);
		    numeral_(&ktb, pb, &size, (ftnlen)4);
		    if (kta <= ktb) {
/* Writing concatenation */
			i__2[0] = 4, a__1[0] = pa;
			i__2[1] = 4, a__1[1] = pb;
			s_cat(pktb, a__1, i__2, &c__2, (ftnlen)8);
		    } else {
/* Writing concatenation */
			i__2[0] = 4, a__1[0] = pb;
			i__2[1] = 4, a__1[1] = pa;
			s_cat(pktb, a__1, i__2, &c__2, (ftnlen)8);
		    }
		    if (s_cmp(pktb, pitb, (ftnlen)8, (ftnlen)8) == 0) {
			bond_1.bk[k - 1] = bond_1.bk[i__ - 1];
			bond_1.bl[k - 1] = bond_1.bl[i__ - 1];
		    }
		}
	    }
	}
    }

/*     translate optimization values to angle bend parameters */

    i__1 = angle_1.nangle;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pa;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pc;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	} else {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pc;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pa;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    ka = iang_ref(1, k);
	    kb = iang_ref(2, k);
	    kc = iang_ref(3, k);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    if (kta <= ktc) {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pa;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pc;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pc;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pa;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    }
	    if (s_cmp(pkta, pita, (ftnlen)12, (ftnlen)12) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_angle__ && angle_1.ak[i__ - 1] != 0.) {
		++(*nvar);
		angle_1.ak[i__ - 1] = xx[*nvar];
		if (*ivar == *nvar) {
		    angle_1.ak[i__ - 1] += *eps;
		}
		++(*nvar);
		angle_1.anat[i__ - 1] = xx[*nvar];
		if (*ivar == *nvar) {
		    angle_1.anat[i__ - 1] += *eps;
		}
		i__3 = angle_1.nangle;
		for (k = i__ + 1; k <= i__3; ++k) {
		    ka = iang_ref(1, k);
		    kb = iang_ref(2, k);
		    kc = iang_ref(3, k);
		    kta = atmtyp_1.class__[ka - 1];
		    ktb = atmtyp_1.class__[kb - 1];
		    ktc = atmtyp_1.class__[kc - 1];
		    size = 4;
		    numeral_(&kta, pa, &size, (ftnlen)4);
		    numeral_(&ktb, pb, &size, (ftnlen)4);
		    numeral_(&ktc, pc, &size, (ftnlen)4);
		    if (kta <= ktc) {
/* Writing concatenation */
			i__4[0] = 4, a__2[0] = pa;
			i__4[1] = 4, a__2[1] = pb;
			i__4[2] = 4, a__2[2] = pc;
			s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
		    } else {
/* Writing concatenation */
			i__4[0] = 4, a__2[0] = pc;
			i__4[1] = 4, a__2[1] = pb;
			i__4[2] = 4, a__2[2] = pa;
			s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
		    }
		    if (s_cmp(pkta, pita, (ftnlen)12, (ftnlen)12) == 0) {
			angle_1.ak[k - 1] = angle_1.ak[i__ - 1];
			angle_1.anat[k - 1] = angle_1.anat[i__ - 1];
		    }
		}
	    }
	}
    }

/*     translate optimization values to stretch-bend parameters */

    i__1 = strbnd_1.nstrbnd;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	ii = isb_ref(1, i__);
	ia = iang_ref(1, ii);
	ib = iang_ref(2, ii);
	ic = iang_ref(3, ii);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pa;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pc;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	} else {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pc;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pa;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    kk = isb_ref(1, k);
	    ka = iang_ref(1, kk);
	    kb = iang_ref(2, kk);
	    kc = iang_ref(3, kk);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    if (kta <= ktc) {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pa;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pc;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pc;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pa;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    }
	    if (s_cmp(pkta, pita, (ftnlen)12, (ftnlen)12) == 0) {
		done = TRUE_;
	    }
	}
	if (! done && valfit_1.fit_strbnd__) {
	    if (sbk_ref(1, i__) != 0.) {
		++(*nvar);
		sbk_ref(1, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    sbk_ref(1, i__) = sbk_ref(1, i__) + *eps;
		}
	    }
	    if (sbk_ref(2, i__) != 0.) {
		++(*nvar);
		sbk_ref(2, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    sbk_ref(2, i__) = sbk_ref(2, i__) + *eps;
		}
	    }
	    i__3 = strbnd_1.nstrbnd;
	    for (k = i__ + 1; k <= i__3; ++k) {
		kk = isb_ref(1, k);
		ka = iang_ref(1, kk);
		kb = iang_ref(2, kk);
		kc = iang_ref(3, kk);
		kta = atmtyp_1.class__[ka - 1];
		ktb = atmtyp_1.class__[kb - 1];
		ktc = atmtyp_1.class__[kc - 1];
		size = 4;
		numeral_(&kta, pa, &size, (ftnlen)4);
		numeral_(&ktb, pb, &size, (ftnlen)4);
		numeral_(&ktc, pc, &size, (ftnlen)4);
		if (kta <= ktc) {
/* Writing concatenation */
		    i__4[0] = 4, a__2[0] = pa;
		    i__4[1] = 4, a__2[1] = pb;
		    i__4[2] = 4, a__2[2] = pc;
		    s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
		} else {
/* Writing concatenation */
		    i__4[0] = 4, a__2[0] = pc;
		    i__4[1] = 4, a__2[1] = pb;
		    i__4[2] = 4, a__2[2] = pa;
		    s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
		}
		if (s_cmp(pkta, pita, (ftnlen)12, (ftnlen)12) == 0) {
		    if (kta == ita) {
			sbk_ref(1, k) = sbk_ref(1, i__);
			sbk_ref(2, k) = sbk_ref(2, i__);
		    } else {
			sbk_ref(2, k) = sbk_ref(1, i__);
			sbk_ref(1, k) = sbk_ref(2, i__);
		    }
		}
	    }
	}
    }

/*     translate optimization values to Urey-Bradley parameters */

    i__1 = urey_1.nurey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	ia = iury_ref(1, i__);
	ib = iury_ref(2, i__);
	ic = iury_ref(3, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pa;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pc;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	} else {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pc;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pa;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    ka = iury_ref(1, k);
	    kb = iury_ref(2, k);
	    kc = iury_ref(3, k);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    if (kta <= ktc) {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pa;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pc;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pc;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pa;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    }
	    if (s_cmp(pkta, pita, (ftnlen)12, (ftnlen)12) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_urey__ && urey_1.uk[i__ - 1] != 0.) {
		++(*nvar);
		urey_1.uk[i__ - 1] = xx[*nvar];
		if (*ivar == *nvar) {
		    urey_1.uk[i__ - 1] += *eps;
		}
		++(*nvar);
		urey_1.ul[i__ - 1] = xx[*nvar];
		if (*ivar == *nvar) {
		    urey_1.ul[i__ - 1] += *eps;
		}
		i__3 = urey_1.nurey;
		for (k = i__ + 1; k <= i__3; ++k) {
		    ka = iury_ref(1, k);
		    kb = iury_ref(2, k);
		    kc = iury_ref(3, k);
		    kta = atmtyp_1.class__[ka - 1];
		    ktb = atmtyp_1.class__[kb - 1];
		    ktc = atmtyp_1.class__[kc - 1];
		    size = 4;
		    numeral_(&kta, pa, &size, (ftnlen)4);
		    numeral_(&ktb, pb, &size, (ftnlen)4);
		    numeral_(&ktc, pc, &size, (ftnlen)4);
		    if (kta <= ktc) {
/* Writing concatenation */
			i__4[0] = 4, a__2[0] = pa;
			i__4[1] = 4, a__2[1] = pb;
			i__4[2] = 4, a__2[2] = pc;
			s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
		    } else {
/* Writing concatenation */
			i__4[0] = 4, a__2[0] = pc;
			i__4[1] = 4, a__2[1] = pb;
			i__4[2] = 4, a__2[2] = pa;
			s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
		    }
		    if (s_cmp(pkta, pita, (ftnlen)12, (ftnlen)12) == 0) {
			urey_1.uk[k - 1] = urey_1.uk[i__ - 1];
			urey_1.ul[k - 1] = urey_1.ul[i__ - 1];
		    }
		}
	    }
	}
    }

/*     translate optimization values to out-of-plane bend parameters */

    i__1 = opbend_1.nopbend;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	ii = opbend_1.iopb[i__ - 1];
	ia = iang_ref(1, ii);
	ib = iang_ref(2, ii);
	ic = iang_ref(3, ii);
	id = iang_ref(4, ii);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	itd = atmtyp_1.class__[id - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	numeral_(&itd, pd, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pd;
	    i__5[1] = 4, a__3[1] = pb;
	    i__5[2] = 4, a__3[2] = pa;
	    i__5[3] = 4, a__3[3] = pc;
	    s_cat(pitt, a__3, i__5, &c__4, (ftnlen)16);
	} else {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pd;
	    i__5[1] = 4, a__3[1] = pb;
	    i__5[2] = 4, a__3[2] = pc;
	    i__5[3] = 4, a__3[3] = pa;
	    s_cat(pitt, a__3, i__5, &c__4, (ftnlen)16);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    kk = opbend_1.iopb[k - 1];
	    ka = iang_ref(1, kk);
	    kb = iang_ref(2, kk);
	    kc = iang_ref(3, kk);
	    kd = iang_ref(4, kk);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    ktd = atmtyp_1.class__[kd - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    numeral_(&ktd, pd, &size, (ftnlen)4);
	    if (kta <= ktc) {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pd;
		i__5[1] = 4, a__3[1] = pb;
		i__5[2] = 4, a__3[2] = pa;
		i__5[3] = 4, a__3[3] = pc;
		s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
	    } else {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pd;
		i__5[1] = 4, a__3[1] = pb;
		i__5[2] = 4, a__3[2] = pc;
		i__5[3] = 4, a__3[3] = pa;
		s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
	    }
	    if (s_cmp(pktt, pitt, (ftnlen)16, (ftnlen)16) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_opbend__ && opbend_1.opbk[i__ - 1] != 0.) {
		++(*nvar);
		opbend_1.opbk[i__ - 1] = xx[*nvar];
		if (*ivar == *nvar) {
		    opbend_1.opbk[i__ - 1] += *eps;
		}
		i__3 = opbend_1.nopbend;
		for (k = i__ + 1; k <= i__3; ++k) {
		    kk = opbend_1.iopb[k - 1];
		    ka = iang_ref(1, kk);
		    kb = iang_ref(2, kk);
		    kc = iang_ref(3, kk);
		    kd = iang_ref(4, kk);
		    kta = atmtyp_1.class__[ka - 1];
		    ktb = atmtyp_1.class__[kb - 1];
		    ktc = atmtyp_1.class__[kc - 1];
		    ktd = atmtyp_1.class__[kd - 1];
		    size = 4;
		    numeral_(&kta, pa, &size, (ftnlen)4);
		    numeral_(&ktb, pb, &size, (ftnlen)4);
		    numeral_(&ktc, pc, &size, (ftnlen)4);
		    numeral_(&ktd, pd, &size, (ftnlen)4);
		    if (kta <= ktc) {
/* Writing concatenation */
			i__5[0] = 4, a__3[0] = pd;
			i__5[1] = 4, a__3[1] = pb;
			i__5[2] = 4, a__3[2] = pa;
			i__5[3] = 4, a__3[3] = pc;
			s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
		    } else {
/* Writing concatenation */
			i__5[0] = 4, a__3[0] = pd;
			i__5[1] = 4, a__3[1] = pb;
			i__5[2] = 4, a__3[2] = pc;
			i__5[3] = 4, a__3[3] = pa;
			s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
		    }
		    if (s_cmp(pktt, pitt, (ftnlen)16, (ftnlen)16) == 0) {
			opbend_1.opbk[k - 1] = opbend_1.opbk[i__ - 1];
		    }
		}
	    }
	}
    }

/*     translate optimization values to torsional parameters */

    i__1 = tors_1.ntors;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	ia = itors_ref(1, i__);
	ib = itors_ref(2, i__);
	ic = itors_ref(3, i__);
	id = itors_ref(4, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	itd = atmtyp_1.class__[id - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	numeral_(&itd, pd, &size, (ftnlen)4);
	if (itb < itc || itb == itc && ita <= itd) {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pa;
	    i__5[1] = 4, a__3[1] = pb;
	    i__5[2] = 4, a__3[2] = pc;
	    i__5[3] = 4, a__3[3] = pd;
	    s_cat(pitt, a__3, i__5, &c__4, (ftnlen)16);
	} else {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pd;
	    i__5[1] = 4, a__3[1] = pc;
	    i__5[2] = 4, a__3[2] = pb;
	    i__5[3] = 4, a__3[3] = pa;
	    s_cat(pitt, a__3, i__5, &c__4, (ftnlen)16);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    ka = itors_ref(1, k);
	    kb = itors_ref(2, k);
	    kc = itors_ref(3, k);
	    kd = itors_ref(4, k);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    ktd = atmtyp_1.class__[kd - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    numeral_(&ktd, pd, &size, (ftnlen)4);
	    if (ktb < ktc || ktb == ktc && kta <= ktd) {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pa;
		i__5[1] = 4, a__3[1] = pb;
		i__5[2] = 4, a__3[2] = pc;
		i__5[3] = 4, a__3[3] = pd;
		s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
	    } else {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pd;
		i__5[1] = 4, a__3[1] = pc;
		i__5[2] = 4, a__3[2] = pb;
		i__5[3] = 4, a__3[3] = pa;
		s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
	    }
	    if (s_cmp(pktt, pitt, (ftnlen)16, (ftnlen)16) == 0) {
		done = TRUE_;
	    }
	}
	if (! done && valfit_1.fit_tors__) {
	    if (tors1_ref(1, i__) != 0.) {
		++(*nvar);
		tors1_ref(1, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    tors1_ref(1, i__) = tors1_ref(1, i__) + *eps;
		}
	    }
	    if (tors2_ref(1, i__) != 0.) {
		++(*nvar);
		tors2_ref(1, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    tors2_ref(1, i__) = tors2_ref(1, i__) + *eps;
		}
	    }
	    if (tors3_ref(1, i__) != 0.) {
		++(*nvar);
		tors3_ref(1, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    tors3_ref(1, i__) = tors3_ref(1, i__) + *eps;
		}
	    }
	    if (tors4_ref(1, i__) != 0.) {
		++(*nvar);
		tors4_ref(1, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    tors4_ref(1, i__) = tors4_ref(1, i__) + *eps;
		}
	    }
	    if (tors5_ref(1, i__) != 0.) {
		++(*nvar);
		tors5_ref(1, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    tors5_ref(1, i__) = tors5_ref(1, i__) + *eps;
		}
	    }
	    if (tors6_ref(1, i__) != 0.) {
		++(*nvar);
		tors6_ref(1, i__) = xx[*nvar];
		if (*ivar == *nvar) {
		    tors6_ref(1, i__) = tors6_ref(1, i__) + *eps;
		}
	    }
	    i__3 = tors_1.ntors;
	    for (k = i__ + 1; k <= i__3; ++k) {
		ka = itors_ref(1, k);
		kb = itors_ref(2, k);
		kc = itors_ref(3, k);
		kd = itors_ref(4, k);
		kta = atmtyp_1.class__[ka - 1];
		ktb = atmtyp_1.class__[kb - 1];
		ktc = atmtyp_1.class__[kc - 1];
		ktd = atmtyp_1.class__[kd - 1];
		size = 4;
		numeral_(&kta, pa, &size, (ftnlen)4);
		numeral_(&ktb, pb, &size, (ftnlen)4);
		numeral_(&ktc, pc, &size, (ftnlen)4);
		numeral_(&ktd, pd, &size, (ftnlen)4);
		if (ktb < ktc || ktb == ktc && kta <= ktd) {
/* Writing concatenation */
		    i__5[0] = 4, a__3[0] = pa;
		    i__5[1] = 4, a__3[1] = pb;
		    i__5[2] = 4, a__3[2] = pc;
		    i__5[3] = 4, a__3[3] = pd;
		    s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
		} else {
/* Writing concatenation */
		    i__5[0] = 4, a__3[0] = pd;
		    i__5[1] = 4, a__3[1] = pc;
		    i__5[2] = 4, a__3[2] = pb;
		    i__5[3] = 4, a__3[3] = pa;
		    s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
		}
		if (s_cmp(pktt, pitt, (ftnlen)16, (ftnlen)16) == 0) {
		    tors1_ref(1, k) = tors1_ref(1, i__);
		    tors2_ref(1, k) = tors2_ref(1, i__);
		    tors3_ref(1, k) = tors3_ref(1, i__);
		    tors4_ref(1, k) = tors4_ref(1, i__);
		    tors5_ref(1, k) = tors5_ref(1, i__);
		    tors6_ref(1, k) = tors6_ref(1, i__);
		}
	    }
	}
    }
    return 0;
} /* varprm_ */

#undef itors_ref
#undef tors6_ref
#undef tors5_ref
#undef tors4_ref
#undef tors3_ref
#undef tors2_ref
#undef tors1_ref
#undef iury_ref
#undef iang_ref
#undef ibnd_ref
#undef sbk_ref
#undef isb_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  function valfit1  --  valence fit error and gradient  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "valfit1" is a service routine that computes the RMS error */
/*     and gradient for valence parameters fit to QM results */


doublereal valfit1_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal e;
    static integer i__, k;
    static doublereal e0, eps[75000];
    static integer nvar;
    static doublereal delta;
    extern doublereal valrms_(integer *);
    extern /* Subroutine */ int varprm_(integer *, doublereal *, integer *, 
	    doublereal *);



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  valfit.i  --  values for valence term parameter fitting  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     fit_bond    logical flag to fit bond stretch parameters */
/*     fit_angle   logical flag to fit angle bend parameters */
/*     fit_strbnd  logical flag to fit stretch-bend parameters */
/*     fit_urey    logical flag to fit Urey-Bradley parameters */
/*     fit_opbend  logical flag to fit out-of-plane bend parameters */
/*     fit_tors    logical flag to fit torsional parameters */
/*     fit_struct  logical flag to structure-fit valence parameters */
/*     fit_force   logical flag to force-fit valence parameters */




/*     copy optimization values to valence parameters */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    varprm_(&nvar, &xx[1], &c__0, &c_b36);

/*     set the numerical gradient step size for each parameter */

    delta = 1e-7;
    i__1 = nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	eps[i__ - 1] = delta * xx[i__];
    }

/*     get the RMS of frequencies */

    ret_val = valrms_(&c__0);

/*     compute numerical gradient for valence parameters */

    k = nvar;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = eps[i__ - 1] * -.5;
	varprm_(&nvar, &xx[1], &i__, &d__1);
	e0 = valrms_(&c__0);
	d__1 = eps[i__ - 1] * .5;
	varprm_(&nvar, &xx[1], &i__, &d__1);
	e = valrms_(&c__0);
	g[i__] = (e - e0) / eps[i__ - 1];
    }
    return ret_val;
} /* valfit1_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##   subroutine prtval  --  print fitted valence parameters  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "prtval" writes the final fitted valence parameters to the */
/*     standard output and appends the values to a key file */


/* Subroutine */ int prtval_(void)
{
    /* Format strings */
    static char fmt_10[] = "(a)";
    static char fmt_20[] = "(/,\002#\002,/,\002# Results of Valence Paramete"
	    "r Fitting\002,/,\002#\002,/)";
    static char fmt_30[] = "(\002bond\002,6x,2i5,5x,f11.2,f11.4)";
    static char fmt_40[] = "(\002angle\002,5x,3i5,f11.2,f11.2)";
    static char fmt_50[] = "(\002strbnd\002,4x,3i5,2f11.3)";
    static char fmt_60[] = "(\002ureybrad\002,2x,3i5,f11.3,f11.4)";
    static char fmt_70[] = "(\002opbend\002,4x,4i5,6x,f11.2)";
    static char fmt_80[] = "(\002torsion\002,3x,4i5,3x,f8.3,\002 0.0 1\002,f"
	    "8.3,\002 180.0 2\002,f8.3,\002 0.0 3\002)";

    /* System generated locals */
    address a__1[2], a__2[3], a__3[4];
    integer i__1[2], i__2, i__3, i__4[3], i__5[4], i__6;
    olist o__1;

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern integer freeunit_(void), trimtext_(char *, ftnlen);
    static integer i__, k, ia, ib, ic, id, ka, ii, kb, kc, kd, kk;
    static char pa[4], pb[4], pc[4], pd[4];
    static integer ita, itb, itc, itd, kta, ktb, ktc, ktd;
    static logical done;
    static char pita[12], pitb[8], pkta[12], pktb[8];
    static integer ikey, size;
    static char pitt[16], pktt[16], record[120], keyfile[120];
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    , version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___364 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___365 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___380 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___388 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___391 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___392 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___400 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___401 = { 0, 0, 0, fmt_80, 0 };



#define isb_ref(a_1,a_2) strbnd_1.isb[(a_2)*3 + a_1 - 4]
#define sbk_ref(a_1,a_2) strbnd_1.sbk[(a_2)*2 + a_1 - 3]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define iury_ref(a_1,a_2) urey_1.iury[(a_2)*3 + a_1 - 4]
#define tors1_ref(a_1,a_2) tors_1.tors1[(a_2)*4 + a_1 - 5]
#define tors2_ref(a_1,a_2) tors_1.tors2[(a_2)*4 + a_1 - 5]
#define tors3_ref(a_1,a_2) tors_1.tors3[(a_2)*4 + a_1 - 5]
#define itors_ref(a_1,a_2) tors_1.itors[(a_2)*4 + a_1 - 5]
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
/*     ##  keys.i  --  contents of current keyword parameter file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nkey      number of nonblank lines in the keyword file */
/*     keyline   contents of each individual keyword file line */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  opbend.i  --  out-of-plane bends in the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     opbk      force constant values for out-of-plane bending */
/*     nopbend   total number of out-of-plane bends in the system */
/*     iopb      bond angle numbers used in out-of-plane bending */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  strbnd.i  --  stretch-bends in the current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     sbk       force constants for stretch-bend terms */
/*     nstrbnd   total number of stretch-bend interactions */
/*     isb       angle and bond numbers used in stretch-bend */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  urey.i  --  Urey-Bradley interactions in the structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     uk      Urey-Bradley force constants (kcal/mole/Ang**2) */
/*     ul      ideal 1-3 distance values in Angstroms */
/*     nurey   total number of Urey-Bradley terms in the system */
/*     iury    numbers of the atoms in each Urey-Bradley interaction */




/*     ############################################################## */
/*     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  valfit.i  --  values for valence term parameter fitting  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     fit_bond    logical flag to fit bond stretch parameters */
/*     fit_angle   logical flag to fit angle bend parameters */
/*     fit_strbnd  logical flag to fit stretch-bend parameters */
/*     fit_urey    logical flag to fit Urey-Bradley parameters */
/*     fit_opbend  logical flag to fit out-of-plane bend parameters */
/*     fit_tors    logical flag to fit torsional parameters */
/*     fit_struct  logical flag to structure-fit valence parameters */
/*     fit_force   logical flag to force-fit valence parameters */




/*     output some definitions and parameters to a keyfile */

    ikey = freeunit_();
/* Writing concatenation */
    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
    i__1[1] = 4, a__1[1] = ".key";
    s_cat(keyfile, a__1, i__1, &c__2, (ftnlen)120);
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

    i__2 = keys_1.nkey;
    for (i__ = 1; i__ <= i__2; ++i__) {
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	size = trimtext_(record, (ftnlen)120);
	io___364.ciunit = ikey;
	s_wsfe(&io___364);
	do_fio(&c__1, record, size);
	e_wsfe();
    }

/*     print a header for the fitted valence parameters */

    if (valfit_1.fit_bond__ || valfit_1.fit_angle__ || valfit_1.fit_tors__ || 
	    valfit_1.fit_strbnd__ || valfit_1.fit_opbend__) {
	io___365.ciunit = ikey;
	s_wsfe(&io___365);
	e_wsfe();
    }

/*     output any fitted bond stretch parameter values */

    i__2 = bond_1.nbond;
    for (i__ = 1; i__ <= i__2; ++i__) {
	done = FALSE_;
	ia = ibnd_ref(1, i__);
	ib = ibnd_ref(2, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	if (ita <= itb) {
/* Writing concatenation */
	    i__1[0] = 4, a__1[0] = pa;
	    i__1[1] = 4, a__1[1] = pb;
	    s_cat(pitb, a__1, i__1, &c__2, (ftnlen)8);
	} else {
/* Writing concatenation */
	    i__1[0] = 4, a__1[0] = pb;
	    i__1[1] = 4, a__1[1] = pa;
	    s_cat(pitb, a__1, i__1, &c__2, (ftnlen)8);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    ka = ibnd_ref(1, k);
	    kb = ibnd_ref(2, k);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    if (kta <= ktb) {
/* Writing concatenation */
		i__1[0] = 4, a__1[0] = pa;
		i__1[1] = 4, a__1[1] = pb;
		s_cat(pktb, a__1, i__1, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__1[0] = 4, a__1[0] = pb;
		i__1[1] = 4, a__1[1] = pa;
		s_cat(pktb, a__1, i__1, &c__2, (ftnlen)8);
	    }
	    if (s_cmp(pktb, pitb, (ftnlen)8, (ftnlen)8) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_bond__) {
		io___380.ciunit = ikey;
		s_wsfe(&io___380);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&bond_1.bk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&bond_1.bl[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     output any fitted angle bend parameter values */

    i__2 = angle_1.nangle;
    for (i__ = 1; i__ <= i__2; ++i__) {
	done = FALSE_;
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pa;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pc;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	} else {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pc;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pa;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    ka = iang_ref(1, k);
	    kb = iang_ref(2, k);
	    kc = iang_ref(3, k);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    if (kta <= ktc) {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pa;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pc;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pc;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pa;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    }
	    if (s_cmp(pkta, pita, (ftnlen)12, (ftnlen)12) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_angle__) {
		io___388.ciunit = ikey;
		s_wsfe(&io___388);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&angle_1.ak[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&angle_1.anat[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     output any fitted stretch-bend parameter values */

    i__2 = strbnd_1.nstrbnd;
    for (i__ = 1; i__ <= i__2; ++i__) {
	done = FALSE_;
	ii = isb_ref(1, i__);
	ia = iang_ref(1, ii);
	ib = iang_ref(2, ii);
	ic = iang_ref(3, ii);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pa;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pc;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	} else {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pc;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pa;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    kk = isb_ref(1, k);
	    ka = iang_ref(1, kk);
	    kb = iang_ref(2, kk);
	    kc = iang_ref(3, kk);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    if (kta <= ktc) {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pa;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pc;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pc;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pa;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    }
	    if (s_cmp(pkta, pita, (ftnlen)12, (ftnlen)12) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_strbnd__) {
		io___391.ciunit = ikey;
		s_wsfe(&io___391);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&sbk_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&sbk_ref(2, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     output any fitted Urey-Bradley parameter values */

    i__2 = urey_1.nurey;
    for (i__ = 1; i__ <= i__2; ++i__) {
	done = FALSE_;
	ia = iury_ref(1, i__);
	ib = iury_ref(2, i__);
	ic = iury_ref(3, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pa;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pc;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	} else {
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = pc;
	    i__4[1] = 4, a__2[1] = pb;
	    i__4[2] = 4, a__2[2] = pa;
	    s_cat(pita, a__2, i__4, &c__3, (ftnlen)12);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    ka = iury_ref(1, k);
	    kb = iury_ref(2, k);
	    kc = iury_ref(3, k);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    if (kta <= ktc) {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pa;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pc;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__4[0] = 4, a__2[0] = pc;
		i__4[1] = 4, a__2[1] = pb;
		i__4[2] = 4, a__2[2] = pa;
		s_cat(pkta, a__2, i__4, &c__3, (ftnlen)12);
	    }
	    if (s_cmp(pkta, pita, (ftnlen)12, (ftnlen)12) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_urey__) {
		io___392.ciunit = ikey;
		s_wsfe(&io___392);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&urey_1.uk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&urey_1.ul[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     output any fitted out-of-plane bend parameter values */

    i__2 = opbend_1.nopbend;
    for (i__ = 1; i__ <= i__2; ++i__) {
	done = FALSE_;
	ii = opbend_1.iopb[i__ - 1];
	ia = iang_ref(1, ii);
	ib = iang_ref(2, ii);
	ic = iang_ref(3, ii);
	id = iang_ref(4, ii);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	itd = atmtyp_1.class__[id - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	numeral_(&itd, pd, &size, (ftnlen)4);
	if (ita <= itc) {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pd;
	    i__5[1] = 4, a__3[1] = pb;
	    i__5[2] = 4, a__3[2] = pa;
	    i__5[3] = 4, a__3[3] = pc;
	    s_cat(pitt, a__3, i__5, &c__4, (ftnlen)16);
	} else {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pd;
	    i__5[1] = 4, a__3[1] = pb;
	    i__5[2] = 4, a__3[2] = pc;
	    i__5[3] = 4, a__3[3] = pa;
	    s_cat(pitt, a__3, i__5, &c__4, (ftnlen)16);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    kk = opbend_1.iopb[k - 1];
	    ka = iang_ref(1, kk);
	    kb = iang_ref(2, kk);
	    kc = iang_ref(3, kk);
	    kd = iang_ref(4, kk);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    ktd = atmtyp_1.class__[kd - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    numeral_(&ktd, pd, &size, (ftnlen)4);
	    if (kta <= ktc) {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pd;
		i__5[1] = 4, a__3[1] = pb;
		i__5[2] = 4, a__3[2] = pa;
		i__5[3] = 4, a__3[3] = pc;
		s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
	    } else {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pd;
		i__5[1] = 4, a__3[1] = pb;
		i__5[2] = 4, a__3[2] = pc;
		i__5[3] = 4, a__3[3] = pa;
		s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
	    }
	    if (s_cmp(pktt, pitt, (ftnlen)16, (ftnlen)16) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_opbend__) {
		io___400.ciunit = ikey;
		s_wsfe(&io___400);
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		i__3 = min(ita,itc);
		do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		i__6 = max(ita,itc);
		do_fio(&c__1, (char *)&i__6, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&opbend_1.opbk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     output any fitted torsional parameter values */

    i__2 = tors_1.ntors;
    for (i__ = 1; i__ <= i__2; ++i__) {
	done = FALSE_;
	ia = itors_ref(1, i__);
	ib = itors_ref(2, i__);
	ic = itors_ref(3, i__);
	id = itors_ref(4, i__);
	ita = atmtyp_1.class__[ia - 1];
	itb = atmtyp_1.class__[ib - 1];
	itc = atmtyp_1.class__[ic - 1];
	itd = atmtyp_1.class__[id - 1];
	size = 4;
	numeral_(&ita, pa, &size, (ftnlen)4);
	numeral_(&itb, pb, &size, (ftnlen)4);
	numeral_(&itc, pc, &size, (ftnlen)4);
	numeral_(&itd, pd, &size, (ftnlen)4);
	if (itb < itc || itb == itc && ita <= itd) {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pa;
	    i__5[1] = 4, a__3[1] = pb;
	    i__5[2] = 4, a__3[2] = pc;
	    i__5[3] = 4, a__3[3] = pd;
	    s_cat(pitt, a__3, i__5, &c__4, (ftnlen)16);
	} else {
/* Writing concatenation */
	    i__5[0] = 4, a__3[0] = pd;
	    i__5[1] = 4, a__3[1] = pc;
	    i__5[2] = 4, a__3[2] = pb;
	    i__5[3] = 4, a__3[3] = pa;
	    s_cat(pitt, a__3, i__5, &c__4, (ftnlen)16);
	}
	i__3 = i__ - 1;
	for (k = 1; k <= i__3; ++k) {
	    ka = itors_ref(1, k);
	    kb = itors_ref(2, k);
	    kc = itors_ref(3, k);
	    kd = itors_ref(4, k);
	    kta = atmtyp_1.class__[ka - 1];
	    ktb = atmtyp_1.class__[kb - 1];
	    ktc = atmtyp_1.class__[kc - 1];
	    ktd = atmtyp_1.class__[kd - 1];
	    size = 4;
	    numeral_(&kta, pa, &size, (ftnlen)4);
	    numeral_(&ktb, pb, &size, (ftnlen)4);
	    numeral_(&ktc, pc, &size, (ftnlen)4);
	    numeral_(&ktd, pd, &size, (ftnlen)4);
	    if (ktb < ktc || ktb == ktc && kta <= ktd) {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pa;
		i__5[1] = 4, a__3[1] = pb;
		i__5[2] = 4, a__3[2] = pc;
		i__5[3] = 4, a__3[3] = pd;
		s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
	    } else {
/* Writing concatenation */
		i__5[0] = 4, a__3[0] = pd;
		i__5[1] = 4, a__3[1] = pc;
		i__5[2] = 4, a__3[2] = pb;
		i__5[3] = 4, a__3[3] = pa;
		s_cat(pktt, a__3, i__5, &c__4, (ftnlen)16);
	    }
	    if (s_cmp(pktt, pitt, (ftnlen)16, (ftnlen)16) == 0) {
		done = TRUE_;
	    }
	}
	if (! done) {
	    if (valfit_1.fit_tors__) {
		io___401.ciunit = ikey;
		s_wsfe(&io___401);
		do_fio(&c__1, (char *)&ita, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itb, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itc, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&itd, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&tors1_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&tors2_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&tors3_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }
    return 0;
} /* prtval_ */

#undef keyline_ref
#undef itors_ref
#undef tors3_ref
#undef tors2_ref
#undef tors1_ref
#undef iury_ref
#undef iang_ref
#undef ibnd_ref
#undef sbk_ref
#undef isb_ref


/* Main program alias */ int valence_ () { MAIN__ (); return 0; }
