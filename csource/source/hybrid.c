/* hybrid.f -- translated by f2c (version 20050501).
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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal lambda, vlambda, clambda, dlambda, mlambda, plambda;
    integer nmut, imut[25000], type0[25000], type1[25000], class0[25000], 
	    class1[25000];
    logical mut[25000];
} mutant_;

#define mutant_1 mutant_

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
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    doublereal weight[5000];
    integer atmcls[5000], atmnum[5000], ligand[5000];
    char symbol[15000], describe[120000];
} katoms_;

#define katoms_1 katoms_

struct {
    doublereal bk[50000], bl[50000];
    integer nbond, ibnd[100000]	/* was [2][50000] */;
} bond_;

#define bond_1 bond_

struct {
    doublereal bcon[2000], blen[2000], bcon5[500], blen5[500], bcon4[500], 
	    blen4[500], bcon3[500], blen3[500], dlen[500];
    char kb[16000], kb5[4000], kb4[4000], kb3[4000], kel[6000];
} kbonds_;

#define kbonds_1 kbonds_

struct {
    doublereal ak[75000], anat[75000], afld[75000];
    integer nangle, iang[300000]	/* was [4][75000] */;
} angle_;

#define angle_1 angle_

struct {
    doublereal acon[2000], acon5[500], acon4[500], acon3[500], aconf[500], 
	    ang[6000]	/* was [3][2000] */, ang5[1500]	/* was [3][500] */, 
	    ang4[1500]	/* was [3][500] */, ang3[1500]	/* was [3][500] */, 
	    angf[1000]	/* was [2][500] */;
    char ka[24000], ka5[6000], ka4[6000], ka3[6000], kaf[6000];
} kangs_;

#define kangs_1 kangs_

struct {
    integer bndlist[200000]	/* was [8][25000] */, anglist[700000]	/* 
	    was [28][25000] */;
} atmlst_;

#define atmlst_1 atmlst_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal stbn[4000]	/* was [2][2000] */;
    char ksb[24000];
} kstbnd_;

#define kstbnd_1 kstbnd_

struct {
    doublereal sbk[150000]	/* was [2][75000] */;
    integer nstrbnd, isb[225000]	/* was [3][75000] */;
} strbnd_;

#define strbnd_1 strbnd_

struct {
    doublereal itors1[400000]	/* was [4][100000] */, itors2[400000]	/* 
	    was [4][100000] */, itors3[400000]	/* was [4][100000] */;
    integer nitors, iitors[400000]	/* was [4][100000] */;
} imptor_;

#define imptor_1 imptor_

struct {
    doublereal ti1[1000]	/* was [2][500] */, ti2[1000]	/* was [2][
	    500] */, ti3[1000]	/* was [2][500] */;
    char kti[8000];
} kitors_;

#define kitors_1 kitors_

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
    doublereal tors1[400000]	/* was [4][100000] */, tors2[400000]	/* 
	    was [4][100000] */, tors3[400000]	/* was [4][100000] */, tors4[
	    400000]	/* was [4][100000] */, tors5[400000]	/* was [4][
	    100000] */, tors6[400000]	/* was [4][100000] */;
    integer ntors, itors[400000]	/* was [4][100000] */;
} tors_;

#define tors_1 tors_

struct {
    doublereal btcon[1500]	/* was [3][500] */;
    char kbt[8000];
} ksttor_;

#define ksttor_1 ksttor_

struct {
    doublereal kst[300000]	/* was [3][100000] */;
    integer nstrtor, ist[200000]	/* was [2][100000] */;
} strtor_;

#define strtor_1 strtor_

struct {
    doublereal rad[5000], eps[5000], rad4[5000], eps4[5000], reduct[5000];
} kvdws_;

#define kvdws_1 kvdws_

struct {
    doublereal radmin[1000000]	/* was [1000][1000] */, epsilon[1000000]	
	    /* was [1000][1000] */, radmin4[1000000]	/* was [1000][1000] */
	    , epsilon4[1000000]	/* was [1000][1000] */, radhbnd[1000000]	
	    /* was [1000][1000] */, epshbnd[1000000]	/* was [1000][1000] */
	    , kred[25000];
    integer ired[25000], nvdw, ivdw[25000], jvdw[25000], nvt, ivt[25000], jvt[
	    25000];
} vdw_;

#define vdw_1 vdw_

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
    doublereal pchg[25000];
    integer nion, iion[25000], jion[25000], kion[25000], chglist[25000];
} charge_;

#define charge_1 charge_

struct {
    doublereal chg[5000];
} kchrge_;

#define kchrge_1 kchrge_

struct {
    doublereal bdpl[50000], sdpl[50000];
    integer ndipole, idpl[100000]	/* was [2][50000] */;
} dipole_;

#define dipole_1 dipole_

struct {
    doublereal dpl[1000], dpl5[500], dpl4[500], dpl3[500], pos[1000], pos5[
	    500], pos4[500], pos3[500];
    char kd[8000], kd5[4000], kd4[4000], kd3[4000];
} kdipol_;

#define kdipol_1 kdipol_

struct {
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

struct {
    doublereal polr[5000], athl[5000];
    integer pgrp[40000]	/* was [8][5000] */;
} kpolr_;

#define kpolr_1 kpolr_

struct {
    doublereal polarity[25000], thole[25000], pdamp[25000], uind[75000]	/* 
	    was [3][25000] */, uinp[75000]	/* was [3][25000] */, uinds[
	    75000]	/* was [3][25000] */, uinps[75000]	/* was [3][
	    25000] */;
    integer npolar;
} polar_;

#define polar_1 polar_

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;



/*     ############################################################### */
/*     ##  COPYRIGHT (C) 1991 by Shawn Huston & Jay William Ponder  ## */
/*     ##                    All Rights Reserved                    ## */
/*     ############################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine hybrid  --  set chimeric force field parameters  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "hybrid" constructs the hybrid hamiltonian for a specified */
/*     initial state, final state and mutation parameter "lambda" */


/* Subroutine */ int hybrid_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Lambda Coupling Parameter for FEP :\002,"
	    "f12.3)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int hvdw_(void), hbond_(void), hatom_(void), 
	    htors_(void), hangle_(void), hmpole_(void), hpolar_(void), 
	    hcharge_(void), hdipole_(void), hstrbnd_(void), himptor_(void), 
	    hstrtor_(void);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };




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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     set the potential energy parameters for hybrid atoms */



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */


    if (mutant_1.nmut != 0) {
	io___1.ciunit = iounit_1.iout;
	s_wsfe(&io___1);
	do_fio(&c__1, (char *)&mutant_1.lambda, (ftnlen)sizeof(doublereal));
	e_wsfe();
	hatom_();
	hbond_();
	hangle_();
	hstrbnd_();
	himptor_();
	htors_();
	hstrtor_();
	hvdw_();
	hcharge_();
	hdipole_();
	hmpole_();
	hpolar_();
    }
    return 0;
} /* hybrid_ */



/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine hatom  --  assign hybrid atom parameters  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "hatom" assigns a new atom type to each hybrid site */


/* Subroutine */ int hatom_(void)
{
    /* Format strings */
    static char fmt_20[] = "(\002 HATOM  --  Too many Sites as Hybrid Atoms"
	    ";\002,\002 Increase MAXTYP\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, k, it, it0, it1, ntype;

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_20, 0 };



#define describe_ref(a_0,a_1) &katoms_1.describe[(a_1)*24 + a_0 - 24]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define story_ref(a_0,a_1) &atmtyp_1.story[(a_1)*24 + a_0 - 24]
#define symbol_ref(a_0,a_1) &katoms_1.symbol[(a_1)*3 + a_0 - 3]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




/*     find the total number of atom types currently used; */
/*     exclude the "HYB" types so that they can be reused */

    for (i__ = 1; i__ <= 5000; ++i__) {
	if (s_cmp(symbol_ref(0, i__), "   ", (ftnlen)3, (ftnlen)3) == 0 || 
		s_cmp(symbol_ref(0, i__), "HYB", (ftnlen)3, (ftnlen)3) == 0) {
	    ntype = i__ - 1;
	    goto L10;
	}
    }
L10:

/*     stop if there are too many atom types required */

    if (5000 < ntype + mutant_1.nmut) {
	inform_1.abort = TRUE_;
	io___4.ciunit = iounit_1.iout;
	s_wsfe(&io___4);
	e_wsfe();
    }

/*     create a new atom type for each of the hybrid atoms */

    i__1 = mutant_1.nmut;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = mutant_1.imut[i__ - 1];
	it = ntype + i__;
	it0 = mutant_1.type0[i__ - 1];
	it1 = mutant_1.type1[i__ - 1];
	s_copy(symbol_ref(0, it), "HYB", (ftnlen)3, (ftnlen)3);
	katoms_1.atmnum[it - 1] = 0;
	katoms_1.weight[it - 1] = mutant_1.lambda * katoms_1.weight[it1 - 1] 
		+ (1. - mutant_1.lambda) * katoms_1.weight[it0 - 1];
	katoms_1.ligand[it - 1] = 0;
	s_copy(describe_ref(0, it), "Hybrid Atom Type        ", (ftnlen)24, (
		ftnlen)24);
	atoms_1.type__[k - 1] = it;
	s_copy(name___ref(0, k), symbol_ref(0, it), (ftnlen)3, (ftnlen)3);
	atmtyp_1.atomic[k - 1] = katoms_1.atmnum[it - 1];
	atmtyp_1.mass[k - 1] = katoms_1.weight[it - 1];
	atmtyp_1.valence[k - 1] = katoms_1.ligand[it - 1];
	s_copy(story_ref(0, k), describe_ref(0, it), (ftnlen)24, (ftnlen)24);
    }
    return 0;
} /* hatom_ */

#undef symbol_ref
#undef story_ref
#undef name___ref
#undef describe_ref




/*     ######################################################### */
/*     ##                                                     ## */
/*     ##  subroutine hbond  --  find hybrid bond parameters  ## */
/*     ##                                                     ## */
/*     ######################################################### */


/*     "hbond" constructs hybrid bond stretch parameters given */
/*     an initial state, final state and "lambda" value */


/* Subroutine */ int hbond_(void)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 Hybrid Bond Stretching Parameters :\002,"
	    "//,6x,\002Atom Numbers\002,9x,\002KS\002,7x,\002Length\002,/)";
    static char fmt_40[] = "(6x,2i5,f14.3,f12.4)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, ia, ib;
    static char pa[4], pb[4], pt[8];
    static doublereal bk0, bk1, bl0, bl1;
    static integer ita, itb, size;
    static logical header;
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    ;

    /* Fortran I/O blocks */
    static cilist io___25 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_40, 0 };



#define kb_ref(a_0,a_1) &kbonds_1.kb[(a_1)*8 + a_0 - 8]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




/*     assign the hybrid parameters for individual bonds */

    header = TRUE_;
    i__1 = bond_1.nbond;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ibnd_ref(1, i__);
	ib = ibnd_ref(2, i__);
	if (mutant_1.mut[ia - 1] || mutant_1.mut[ib - 1]) {
	    ita = atmtyp_1.class__[ia - 1];
	    itb = atmtyp_1.class__[ib - 1];

/*     find the bond parameters for the initial state */

	    i__2 = mutant_1.nmut;
	    for (j = 1; j <= i__2; ++j) {
		k = mutant_1.imut[j - 1];
		if (k == ia) {
		    ita = mutant_1.class0[j - 1];
		}
		if (k == ib) {
		    itb = mutant_1.class0[j - 1];
		}
	    }
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    if (ita <= itb) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pb;
		i__3[1] = 4, a__1[1] = pa;
		s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
	    }
	    bk0 = 0.;
	    bl0 = 0.;
	    for (j = 1; j <= 2000; ++j) {
		if (s_cmp(kb_ref(0, j), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    bk0 = kbonds_1.bcon[j - 1];
		    bl0 = kbonds_1.blen[j - 1];
		    goto L10;
		}
	    }
L10:

/*     find the bond parameters for the final state */

	    i__2 = mutant_1.nmut;
	    for (j = 1; j <= i__2; ++j) {
		k = mutant_1.imut[j - 1];
		if (k == ia) {
		    ita = mutant_1.class1[j - 1];
		}
		if (k == ib) {
		    itb = mutant_1.class1[j - 1];
		}
	    }
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    if (ita <= itb) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pb;
		i__3[1] = 4, a__1[1] = pa;
		s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
	    }
	    bk1 = 0.;
	    bl1 = 0.;
	    for (j = 1; j <= 2000; ++j) {
		if (s_cmp(kb_ref(0, j), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    bk1 = kbonds_1.bcon[j - 1];
		    bl1 = kbonds_1.blen[j - 1];
		    goto L20;
		}
	    }
L20:

/*     form the hybrid parameters for the current bond */

	    if (bl0 == 0.) {
		bl0 = bl1;
	    }
	    if (bl1 == 0.) {
		bl1 = bl0;
	    }
	    bond_1.bk[i__ - 1] = mutant_1.lambda * bk1 + (1. - 
		    mutant_1.lambda) * bk0;
	    bond_1.bl[i__ - 1] = mutant_1.lambda * bl1 + (1. - 
		    mutant_1.lambda) * bl0;
	    if (inform_1.verbose) {
		if (header) {
		    header = FALSE_;
		    io___25.ciunit = iounit_1.iout;
		    s_wsfe(&io___25);
		    e_wsfe();
		}
		io___26.ciunit = iounit_1.iout;
		s_wsfe(&io___26);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&bond_1.bk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&bond_1.bl[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }
    return 0;
} /* hbond_ */

#undef ibnd_ref
#undef kb_ref




/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine hangle  --  find hybrid angle parameters  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "hangle" constructs hybrid angle bending parameters given */
/*     an initial state, final state and "lambda" value */


/* Subroutine */ int hangle_(void)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 Hybrid Angle Bending Parameters :\002,//"
	    ",6x,\002Atom Numbers\002,9x,\002KB\002,8x,\002Angle\002,/)";
    static char fmt_40[] = "(3x,3i5,2f12.3)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2, i__3[3];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, ia, ib, ic;
    static char pa[4], pb[4], pc[4], pt[12];
    static doublereal ak0, ak1;
    static integer ita, itb, itc, size;
    static doublereal anat0, anat1;
    static logical header;
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    ;

    /* Fortran I/O blocks */
    static cilist io___46 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_40, 0 };



#define ka_ref(a_0,a_1) &kangs_1.ka[(a_1)*12 + a_0 - 12]
#define ang_ref(a_1,a_2) kangs_1.ang[(a_2)*3 + a_1 - 4]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




/*     assign the hybrid parameters for individual angles */

    header = TRUE_;
    i__1 = angle_1.nangle;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	if (mutant_1.mut[ia - 1] || mutant_1.mut[ib - 1] || mutant_1.mut[ic - 
		1]) {
	    ita = atmtyp_1.class__[ia - 1];
	    itb = atmtyp_1.class__[ib - 1];
	    itc = atmtyp_1.class__[ic - 1];

/*     find the angle parameters for the initial state */

	    i__2 = mutant_1.nmut;
	    for (j = 1; j <= i__2; ++j) {
		k = mutant_1.imut[j - 1];
		if (k == ia) {
		    ita = mutant_1.class0[j - 1];
		}
		if (k == ib) {
		    itb = mutant_1.class0[j - 1];
		}
		if (k == ic) {
		    itc = mutant_1.class0[j - 1];
		}
	    }
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    numeral_(&itc, pc, &size, (ftnlen)4);
	    if (ita <= itc) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pc;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pa;
		s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	    }
	    ak0 = 0.;
	    anat0 = 0.;
	    for (j = 1; j <= 2000; ++j) {
		if (s_cmp(ka_ref(0, j), pt, (ftnlen)12, (ftnlen)12) == 0) {
		    ak0 = kangs_1.acon[j - 1];
		    anat0 = ang_ref(1, j);
		    goto L10;
		}
	    }
L10:

/*     find the angle parameters for the final state */

	    i__2 = mutant_1.nmut;
	    for (j = 1; j <= i__2; ++j) {
		k = mutant_1.imut[j - 1];
		if (k == ia) {
		    ita = mutant_1.class1[j - 1];
		}
		if (k == ib) {
		    itb = mutant_1.class1[j - 1];
		}
		if (k == ic) {
		    itc = mutant_1.class1[j - 1];
		}
	    }
	    size = 4;
	    numeral_(&ita, pa, &c__3, (ftnlen)4);
	    numeral_(&itb, pb, &c__3, (ftnlen)4);
	    numeral_(&itc, pc, &c__3, (ftnlen)4);
	    if (ita <= itc) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pc;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pa;
		s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	    }
	    ak1 = 0.;
	    anat1 = 0.;
	    for (j = 1; j <= 2000; ++j) {
		if (s_cmp(ka_ref(0, j), pt, (ftnlen)12, (ftnlen)12) == 0) {
		    ak1 = kangs_1.acon[j - 1];
		    anat1 = ang_ref(1, j);
		    goto L20;
		}
	    }
L20:

/*     form the hybrid parameters for the current angle */

	    if (anat0 == 0.) {
		anat0 = anat1;
	    }
	    if (anat1 == 0.) {
		anat1 = anat0;
	    }
	    angle_1.ak[i__ - 1] = mutant_1.lambda * ak1 + (1. - 
		    mutant_1.lambda) * ak0;
	    angle_1.anat[i__ - 1] = mutant_1.lambda * anat1 + (1. - 
		    mutant_1.lambda) * anat0;
	    if (inform_1.verbose) {
		if (header) {
		    header = FALSE_;
		    io___46.ciunit = iounit_1.iout;
		    s_wsfe(&io___46);
		    e_wsfe();
		}
		io___47.ciunit = iounit_1.iout;
		s_wsfe(&io___47);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&angle_1.ak[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&angle_1.anat[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }
    return 0;
} /* hangle_ */

#undef iang_ref
#undef ang_ref
#undef ka_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine hstrbnd  --  hybrid stretch-bend parameters  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "hstrbnd" constructs hybrid stretch-bend parameters given */
/*     an initial state, final state and "lambda" value */


/* Subroutine */ int hstrbnd_(void)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 Hybrid Stretch-Bend Parameters :\002,//,"
	    "6x,\002Atom Numbers\002,8x,\002KSB 1\002,7x,\002KSB 2\002,/)";
    static char fmt_40[] = "(3x,3i5,2f12.3)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2, i__3[3];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, ia, ib, ic;
    static char pa[4], pb[4], pc[4], pt[12];
    static integer nba, nbc, ita, itb, itc;
    static doublereal sbk0[2], sbk1[2];
    static logical used;
    static integer size;
    static logical header;
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    ;

    /* Fortran I/O blocks */
    static cilist io___68 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___69 = { 0, 0, 0, fmt_40, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define isb_ref(a_1,a_2) strbnd_1.isb[(a_2)*3 + a_1 - 4]
#define ksb_ref(a_0,a_1) &kstbnd_1.ksb[(a_1)*12 + a_0 - 12]
#define sbk_ref(a_1,a_2) strbnd_1.sbk[(a_2)*2 + a_1 - 3]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define stbn_ref(a_1,a_2) kstbnd_1.stbn[(a_2)*2 + a_1 - 3]
#define bndlist_ref(a_1,a_2) atmlst_1.bndlist[(a_2)*8 + a_1 - 9]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  atmlst.i  --  local geometry terms involving each atom  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     bndlist   list of the bond numbers involving each atom */
/*     anglist   list of the angle numbers centered on each atom */




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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




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




/*     assign hybrid parameters for the stretch-bend sites */

    header = TRUE_;
    i__1 = angle_1.nangle;
    for (i__ = 1; i__ <= i__1; ++i__) {
	used = FALSE_;
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	if (mutant_1.mut[ia - 1] || mutant_1.mut[ib - 1] || mutant_1.mut[ic - 
		1]) {
	    ita = atmtyp_1.class__[ia - 1];
	    itb = atmtyp_1.class__[ib - 1];
	    itc = atmtyp_1.class__[ic - 1];
	    i__2 = couple_1.n12[ib - 1];
	    for (j = 1; j <= i__2; ++j) {
		if (i12_ref(j, ib) == ia) {
		    nba = bndlist_ref(j, ib);
		}
		if (i12_ref(j, ib) == ic) {
		    nbc = bndlist_ref(j, ib);
		}
	    }

/*     find the stretch-bend parameters for the initial state */

	    i__2 = mutant_1.nmut;
	    for (j = 1; j <= i__2; ++j) {
		k = mutant_1.imut[j - 1];
		if (k == ia) {
		    ita = mutant_1.class0[j - 1];
		}
		if (k == ib) {
		    itb = mutant_1.class0[j - 1];
		}
		if (k == ic) {
		    itc = mutant_1.class0[j - 1];
		}
	    }
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    numeral_(&itc, pc, &size, (ftnlen)4);
	    if (ita <= itc) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pc;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pa;
		s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	    }
	    sbk0[0] = 0.;
	    sbk0[1] = 0.;
	    for (j = 1; j <= 2000; ++j) {
		if (s_cmp(ksb_ref(0, j), pt, (ftnlen)12, (ftnlen)12) == 0) {
		    used = TRUE_;
		    if (ita <= itc) {
			sbk0[0] = stbn_ref(1, j);
			sbk0[1] = stbn_ref(2, j);
		    } else {
			sbk0[0] = stbn_ref(2, j);
			sbk0[1] = stbn_ref(1, j);
		    }
		    goto L10;
		}
	    }
L10:

/*     find the stretch-bend parameters for the final state */

	    i__2 = mutant_1.nmut;
	    for (j = 1; j <= i__2; ++j) {
		k = mutant_1.imut[j - 1];
		if (k == ia) {
		    ita = mutant_1.class1[j - 1];
		}
		if (k == ib) {
		    itb = mutant_1.class1[j - 1];
		}
		if (k == ic) {
		    itc = mutant_1.class1[j - 1];
		}
	    }
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    numeral_(&itc, pc, &size, (ftnlen)4);
	    if (ita <= itc) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pc;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pa;
		s_cat(pt, a__1, i__3, &c__3, (ftnlen)12);
	    }
	    sbk1[0] = 0.;
	    sbk1[1] = 0.;
	    for (j = 1; j <= 2000; ++j) {
		if (s_cmp(ksb_ref(0, j), pt, (ftnlen)12, (ftnlen)12) == 0) {
		    used = TRUE_;
		    if (ita <= itc) {
			sbk1[0] = stbn_ref(1, j);
			sbk1[1] = stbn_ref(2, j);
		    } else {
			sbk1[0] = stbn_ref(2, j);
			sbk1[1] = stbn_ref(1, j);
		    }
		    goto L20;
		}
	    }
L20:

/*     form hybrid parameters for the current stretch-bend */

	    if (used) {
		++strbnd_1.nstrbnd;
		k = strbnd_1.nstrbnd;
		isb_ref(1, k) = i__;
		isb_ref(2, k) = nba;
		isb_ref(3, k) = nbc;
		sbk_ref(1, k) = mutant_1.lambda * sbk1[0] + (1. - 
			mutant_1.lambda) * sbk0[0];
		sbk_ref(2, k) = mutant_1.lambda * sbk1[1] + (1. - 
			mutant_1.lambda) * sbk0[1];
		if (inform_1.verbose) {
		    if (header) {
			header = FALSE_;
			io___68.ciunit = iounit_1.iout;
			s_wsfe(&io___68);
			e_wsfe();
		    }
		    io___69.ciunit = iounit_1.iout;
		    s_wsfe(&io___69);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&sbk_ref(1, i__), (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&sbk_ref(2, i__), (ftnlen)sizeof(
			    doublereal));
		    e_wsfe();
		}
	    }
	}
    }
    return 0;
} /* hstrbnd_ */

#undef bndlist_ref
#undef stbn_ref
#undef iang_ref
#undef sbk_ref
#undef ksb_ref
#undef isb_ref
#undef i12_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine himptor  --  find hybrid improper torsions  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "himptor" constructs hybrid improper torsional parameters */
/*     given an initial state, final state and "lambda" value */

/*     note this version does not handle multiple parameters at */
/*     a single trigonal site */


/* Subroutine */ int himptor_(void)
{
    /* Format strings */
    static char fmt_40[] = "(/,\002 Hybrid Improper Torsional\002,\002 Param"
	    "eters :\002,//,6x,\002Atom Numbers\002,16x,\002KIT1\002,13x,\002"
	    "KIT2\002,13x,\002KIT3\002,/)";
    static char fmt_50[] = "(1x,4i5,4x,3(f10.4,f7.1))";

    /* System generated locals */
    address a__1[4];
    integer i__1, i__2, i__3[4];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double cos(doublereal), sin(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, ia, ib, ic, id;
    static char pa[4], pb[4], pc[4], pd[4], pt[16*6], pt0[16];
    static doublereal s1_0__, s2_0__, s3_0__, v1_0__, v2_0__, v3_0__, v1_1__, 
	    v2_1__, v3_1__, s1_1__;
    static integer ita, itb, itc, itd;
    static doublereal s2_1__, s3_1__;
    static integer nti;
    static logical used;
    static integer size;
    static doublereal symm, angle;
    static char blank[16], zeros[4];
    static logical header;
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    ;

    /* Fortran I/O blocks */
    static cilist io___107 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___108 = { 0, 0, 0, fmt_50, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define pt_ref(a_0,a_1) &pt[(a_1)*16 + a_0 - 16]
#define ti1_ref(a_1,a_2) kitors_1.ti1[(a_2)*2 + a_1 - 3]
#define ti2_ref(a_1,a_2) kitors_1.ti2[(a_2)*2 + a_1 - 3]
#define ti3_ref(a_1,a_2) kitors_1.ti3[(a_2)*2 + a_1 - 3]
#define kti_ref(a_0,a_1) &kitors_1.kti[(a_1)*16 + a_0 - 16]
#define itors1_ref(a_1,a_2) imptor_1.itors1[(a_2)*4 + a_1 - 5]
#define itors2_ref(a_1,a_2) imptor_1.itors2[(a_2)*4 + a_1 - 5]
#define itors3_ref(a_1,a_2) imptor_1.itors3[(a_2)*4 + a_1 - 5]
#define iitors_ref(a_1,a_2) imptor_1.iitors[(a_2)*4 + a_1 - 5]



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  imptor.i  --  improper torsions in the current structure  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     itors1   1-fold amplitude and phase for each improper torsion */
/*     itors2   2-fold amplitude and phase for each improper torsion */
/*     itors3   3-fold amplitude and phase for each improper torsion */
/*     nitors   total number of improper torsional angles in the system */
/*     iitors   numbers of the atoms in each improper torsional angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kitors.i  --  forcefield parameters for improper torsions  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxnti   maximum number of improper torsion parameter entries */

/*     ti1      torsional parameters for improper 1-fold rotation */
/*     ti2      torsional parameters for improper 2-fold rotation */
/*     ti3      torsional parameters for improper 3-fold rotation */
/*     kti      string of atom classes for improper torsional parameters */




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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




/*     construct hybrid improper torsion parameters */

    s_copy(blank, "                ", (ftnlen)16, (ftnlen)16);
    s_copy(zeros, "0000", (ftnlen)4, (ftnlen)4);
    header = TRUE_;

/*     determine the total number of forcefield parameters */

    nti = 500;
    for (i__ = 500; i__ >= 1; --i__) {
	if (s_cmp(kti_ref(0, i__), blank, (ftnlen)16, (ftnlen)16) == 0) {
	    nti = i__ - 1;
	}
    }

/*     construct hybrid improper torsion parameters */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (couple_1.n12[i__ - 1] == 3) {
	    used = FALSE_;
	    ia = i12_ref(1, i__);
	    ib = i12_ref(2, i__);
	    ic = i__;
	    id = i12_ref(3, i__);
	    if (mutant_1.mut[ia - 1] || mutant_1.mut[ib - 1] || mutant_1.mut[
		    ic - 1] || mutant_1.mut[id - 1]) {
		ita = atmtyp_1.class__[ia - 1];
		itb = atmtyp_1.class__[ib - 1];
		itc = atmtyp_1.class__[ic - 1];
		itd = atmtyp_1.class__[id - 1];

/*     find improper torsion parameters for the initial state */

		i__2 = mutant_1.nmut;
		for (j = 1; j <= i__2; ++j) {
		    k = mutant_1.imut[j - 1];
		    if (k == ia) {
			ita = mutant_1.class0[j - 1];
		    }
		    if (k == ib) {
			itb = mutant_1.class0[j - 1];
		    }
		    if (k == ic) {
			itc = mutant_1.class0[j - 1];
		    }
		    if (k == id) {
			itd = mutant_1.class0[j - 1];
		    }
		}
		size = 4;
		numeral_(&ita, pa, &size, (ftnlen)4);
		numeral_(&itb, pb, &size, (ftnlen)4);
		numeral_(&itc, pc, &size, (ftnlen)4);
		numeral_(&itd, pd, &size, (ftnlen)4);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt_ref(0, 1), a__1, i__3, &c__4, (ftnlen)16);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pb;
		i__3[1] = 4, a__1[1] = pa;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt_ref(0, 2), a__1, i__3, &c__4, (ftnlen)16);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pd;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pb;
		s_cat(pt_ref(0, 3), a__1, i__3, &c__4, (ftnlen)16);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pa;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pb;
		s_cat(pt_ref(0, 4), a__1, i__3, &c__4, (ftnlen)16);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pb;
		i__3[1] = 4, a__1[1] = pd;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt_ref(0, 5), a__1, i__3, &c__4, (ftnlen)16);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt_ref(0, 6), a__1, i__3, &c__4, (ftnlen)16);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = zeros;
		i__3[1] = 4, a__1[1] = zeros;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = zeros;
		s_cat(pt0, a__1, i__3, &c__4, (ftnlen)16);
		symm = 1.;
		if (s_cmp(pa, pb, (ftnlen)4, (ftnlen)4) == 0 || s_cmp(pa, pd, 
			(ftnlen)4, (ftnlen)4) == 0 || s_cmp(pb, pd, (ftnlen)4,
			 (ftnlen)4) == 0) {
		    symm = 2.;
		}
		if (s_cmp(pa, pb, (ftnlen)4, (ftnlen)4) == 0 && s_cmp(pa, pd, 
			(ftnlen)4, (ftnlen)4) == 0 && s_cmp(pb, pd, (ftnlen)4,
			 (ftnlen)4) == 0) {
		    symm = 6.;
		}
		v1_0__ = 0.;
		s1_0__ = 0.;
		v2_0__ = 0.;
		s2_0__ = 0.;
		v3_0__ = 0.;
		s3_0__ = 0.;
		i__2 = nti;
		for (j = 1; j <= i__2; ++j) {
		    if (s_cmp(kti_ref(8, j), pc, (ftnlen)4, (ftnlen)4) == 0) {
			for (k = 1; k <= 6; ++k) {
			    if (s_cmp(kti_ref(0, j), pt_ref(0, k), (ftnlen)16,
				     (ftnlen)16) == 0) {
				used = TRUE_;
				v1_0__ = ti1_ref(1, j) / symm;
				s1_0__ = ti1_ref(2, j);
				v2_0__ = ti2_ref(1, j) / symm;
				s2_0__ = ti2_ref(2, j);
				v3_0__ = ti3_ref(1, j) / symm;
				s3_0__ = ti3_ref(2, j);
				goto L10;
			    }
			}
		    }
		}
		i__2 = nti;
		for (j = 1; j <= i__2; ++j) {
		    if (s_cmp(kti_ref(0, j), pt0, (ftnlen)16, (ftnlen)16) == 
			    0) {
			used = TRUE_;
			v1_0__ = ti1_ref(1, j) / symm;
			s1_0__ = ti1_ref(2, j);
			v2_0__ = ti2_ref(1, j) / symm;
			s2_0__ = ti2_ref(2, j);
			v3_0__ = ti3_ref(1, j) / symm;
			s3_0__ = ti3_ref(2, j);
			goto L10;
		    }
		}
L10:

/*     find improper torsion parameters for the final state */

		i__2 = mutant_1.nmut;
		for (j = 1; j <= i__2; ++j) {
		    k = mutant_1.imut[j - 1];
		    if (k == ia) {
			ita = mutant_1.class1[j - 1];
		    }
		    if (k == ib) {
			itb = mutant_1.class1[j - 1];
		    }
		    if (k == ic) {
			itc = mutant_1.class1[j - 1];
		    }
		    if (k == id) {
			itd = mutant_1.class1[j - 1];
		    }
		}
		size = 4;
		numeral_(&ita, pa, &size, (ftnlen)4);
		numeral_(&itb, pb, &size, (ftnlen)4);
		numeral_(&itc, pc, &size, (ftnlen)4);
		numeral_(&itd, pd, &size, (ftnlen)4);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt_ref(0, 1), a__1, i__3, &c__4, (ftnlen)16);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pb;
		i__3[1] = 4, a__1[1] = pa;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt_ref(0, 2), a__1, i__3, &c__4, (ftnlen)16);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pd;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pb;
		s_cat(pt_ref(0, 3), a__1, i__3, &c__4, (ftnlen)16);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pa;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pb;
		s_cat(pt_ref(0, 4), a__1, i__3, &c__4, (ftnlen)16);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pb;
		i__3[1] = 4, a__1[1] = pd;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt_ref(0, 5), a__1, i__3, &c__4, (ftnlen)16);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt_ref(0, 6), a__1, i__3, &c__4, (ftnlen)16);
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = zeros;
		i__3[1] = 4, a__1[1] = zeros;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = zeros;
		s_cat(pt0, a__1, i__3, &c__4, (ftnlen)16);
		symm = 1.;
		if (s_cmp(pa, pb, (ftnlen)4, (ftnlen)4) == 0 || s_cmp(pa, pd, 
			(ftnlen)4, (ftnlen)4) == 0 || s_cmp(pb, pd, (ftnlen)4,
			 (ftnlen)4) == 0) {
		    symm = 2.;
		}
		if (s_cmp(pa, pb, (ftnlen)4, (ftnlen)4) == 0 && s_cmp(pa, pd, 
			(ftnlen)4, (ftnlen)4) == 0 && s_cmp(pb, pd, (ftnlen)4,
			 (ftnlen)4) == 0) {
		    symm = 6.;
		}
		v1_1__ = 0.;
		s1_1__ = 0.;
		v2_1__ = 0.;
		s2_1__ = 0.;
		v3_1__ = 0.;
		s3_1__ = 0.;
		i__2 = nti;
		for (j = 1; j <= i__2; ++j) {
		    if (s_cmp(kti_ref(8, j), pc, (ftnlen)4, (ftnlen)4) == 0) {
			for (k = 1; k <= 6; ++k) {
			    if (s_cmp(kti_ref(0, j), pt_ref(0, k), (ftnlen)16,
				     (ftnlen)16) == 0) {
				used = TRUE_;
				v1_1__ = ti1_ref(1, j) / symm;
				s1_1__ = ti1_ref(2, j);
				v2_1__ = ti2_ref(1, j) / symm;
				s2_1__ = ti2_ref(2, j);
				v3_1__ = ti3_ref(1, j) / symm;
				s3_1__ = ti3_ref(2, j);
				goto L20;
			    }
			}
		    }
		}
		i__2 = nti;
		for (j = 1; j <= i__2; ++j) {
		    if (s_cmp(kti_ref(0, j), pt0, (ftnlen)16, (ftnlen)16) == 
			    0) {
			used = TRUE_;
			v1_1__ = ti1_ref(1, j) / symm;
			s1_1__ = ti1_ref(2, j);
			v2_1__ = ti2_ref(1, j) / symm;
			s2_1__ = ti2_ref(2, j);
			v3_1__ = ti3_ref(1, j) / symm;
			s3_1__ = ti3_ref(2, j);
			goto L20;
		    }
		}
L20:

/*     form hybrid parameters for the current improper torsion */

		if (used) {
		    i__2 = imptor_1.nitors;
		    for (j = 1; j <= i__2; ++j) {
			if (iitors_ref(3, j) == ic) {
			    k = j;
			    goto L30;
			}
		    }
		    ++imptor_1.nitors;
		    k = imptor_1.nitors;
		    iitors_ref(1, k) = ia;
		    iitors_ref(2, k) = ib;
		    iitors_ref(3, k) = ic;
		    iitors_ref(4, k) = id;
L30:
		    if (s1_0__ == 0.) {
			s1_0__ = s1_1__;
		    }
		    if (s2_0__ == 0.) {
			s2_0__ = s2_1__;
		    }
		    if (s3_0__ == 0.) {
			s3_0__ = s3_1__;
		    }
		    if (s1_1__ == 0.) {
			s1_1__ = s1_0__;
		    }
		    if (s2_1__ == 0.) {
			s2_1__ = s2_0__;
		    }
		    if (s3_1__ == 0.) {
			s3_1__ = s3_0__;
		    }
		    itors1_ref(1, k) = mutant_1.lambda * v1_1__ + (1. - 
			    mutant_1.lambda) * v1_0__;
		    itors1_ref(2, k) = mutant_1.lambda * s1_1__ + (1. - 
			    mutant_1.lambda) * s1_0__;
		    angle = itors1_ref(2, k) / 57.29577951308232088;
		    itors1_ref(3, k) = cos(angle);
		    itors1_ref(4, k) = sin(angle);
		    itors2_ref(1, k) = mutant_1.lambda * v2_1__ + (1. - 
			    mutant_1.lambda) * v2_0__;
		    itors2_ref(2, k) = mutant_1.lambda * s2_1__ + (1. - 
			    mutant_1.lambda) * s2_0__;
		    angle = itors2_ref(2, k) / 57.29577951308232088;
		    itors2_ref(3, k) = cos(angle);
		    itors2_ref(4, k) = sin(angle);
		    itors3_ref(1, k) = mutant_1.lambda * v3_1__ + (1. - 
			    mutant_1.lambda) * v3_0__;
		    itors3_ref(2, k) = mutant_1.lambda * s3_1__ + (1. - 
			    mutant_1.lambda) * s3_0__;
		    angle = itors3_ref(2, k) / 57.29577951308232088;
		    itors3_ref(3, k) = cos(angle);
		    itors3_ref(4, k) = sin(angle);
		    if (inform_1.verbose) {
			if (header) {
			    header = FALSE_;
			    io___107.ciunit = iounit_1.iout;
			    s_wsfe(&io___107);
			    e_wsfe();
			}
			ia = iitors_ref(1, i__);
			ib = iitors_ref(2, i__);
			ic = iitors_ref(3, i__);
			id = iitors_ref(4, i__);
			io___108.ciunit = iounit_1.iout;
			s_wsfe(&io___108);
			do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&itors1_ref(1, k), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&itors1_ref(2, k), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&itors2_ref(1, k), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&itors2_ref(2, k), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&itors3_ref(1, k), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&itors3_ref(2, k), (ftnlen)
				sizeof(doublereal));
			e_wsfe();
		    }
		}
	    }
	}
    }
    return 0;
} /* himptor_ */

#undef iitors_ref
#undef itors3_ref
#undef itors2_ref
#undef itors1_ref
#undef kti_ref
#undef ti3_ref
#undef ti2_ref
#undef ti1_ref
#undef pt_ref
#undef i12_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine htors  --  find hybrid torsion parameters  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "htors" constructs hybrid torsional parameters for a given */
/*     initial state, final state and "lambda" value */


/* Subroutine */ int htors_(void)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 Hybrid Torsional Parameters :\002,//,5x"
	    ",\002Atom Numbers\002,6x,\002KT1\002,7x,\002KT2\002,7x,\002KT"
	    "3\002,7x,\002KT4\002,7x,\002KT5\002,7x,\002KT6\002,/)";
    static char fmt_40[] = "(1x,4i4,1x,6(f6.2,i4))";

    /* System generated locals */
    address a__1[4], a__2[3];
    integer i__1, i__2, i__3[4], i__4[3], i__5, i__6, i__7, i__8, i__9;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double cos(doublereal), sin(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen),
	     i_dnnt(doublereal *);

    /* Local variables */
    static integer i__, j, k, ia, ib, ic, id;
    static char pa[4], pb[4], pc[4], pd[4], pt[16], pt0[16];
    static doublereal s1_0__, s2_0__, s3_0__, v1_0__, v2_0__, v3_0__, v4_0__, 
	    v5_0__, v6_0__, s4_0__;
    static integer ita, itb, itc, itd;
    static doublereal s5_0__, s6_0__, v1_1__, v2_1__, v3_1__, v4_1__, v5_1__, 
	    v6_1__, s1_1__, s2_1__, s3_1__, s4_1__, s5_1__, s6_1__;
    static integer size;
    static doublereal angle;
    static char zeros[4];
    static logical header;
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    ;

    /* Fortran I/O blocks */
    static cilist io___154 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___155 = { 0, 0, 0, fmt_40, 0 };



#define t1_ref(a_1,a_2) ktorsn_1.t1[(a_2)*2 + a_1 - 3]
#define t2_ref(a_1,a_2) ktorsn_1.t2[(a_2)*2 + a_1 - 3]
#define t3_ref(a_1,a_2) ktorsn_1.t3[(a_2)*2 + a_1 - 3]
#define t4_ref(a_1,a_2) ktorsn_1.t4[(a_2)*2 + a_1 - 3]
#define t5_ref(a_1,a_2) ktorsn_1.t5[(a_2)*2 + a_1 - 3]
#define t6_ref(a_1,a_2) ktorsn_1.t6[(a_2)*2 + a_1 - 3]
#define kt_ref(a_0,a_1) &ktorsn_1.kt[(a_1)*16 + a_0 - 16]
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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




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




/*     construct hybrid torsional parameters */

    s_copy(zeros, "0000", (ftnlen)4, (ftnlen)4);
    header = TRUE_;
    i__1 = tors_1.ntors;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = itors_ref(1, i__);
	ib = itors_ref(2, i__);
	ic = itors_ref(3, i__);
	id = itors_ref(4, i__);
	if (mutant_1.mut[ia - 1] || mutant_1.mut[ib - 1] || mutant_1.mut[ic - 
		1] || mutant_1.mut[id - 1]) {
	    ita = atmtyp_1.class__[ia - 1];
	    itb = atmtyp_1.class__[ib - 1];
	    itc = atmtyp_1.class__[ic - 1];
	    itd = atmtyp_1.class__[id - 1];

/*     find the torsion parameters for the initial state */

	    i__2 = mutant_1.nmut;
	    for (j = 1; j <= i__2; ++j) {
		k = mutant_1.imut[j - 1];
		if (k == ia) {
		    ita = mutant_1.class0[j - 1];
		}
		if (k == ib) {
		    itb = mutant_1.class0[j - 1];
		}
		if (k == ic) {
		    itc = mutant_1.class0[j - 1];
		}
		if (k == id) {
		    itd = mutant_1.class0[j - 1];
		}
	    }
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    numeral_(&itc, pc, &size, (ftnlen)4);
	    numeral_(&itd, pd, &size, (ftnlen)4);
	    if (itb < itc) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (itc < itb) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (ita <= itd) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (itd < ita) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    }
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = zeros;
	    i__4[1] = 8, a__2[1] = pt + 4;
	    i__4[2] = 4, a__2[2] = zeros;
	    s_cat(pt0, a__2, i__4, &c__3, (ftnlen)16);
	    v1_0__ = 0.;
	    s1_0__ = 0.;
	    v2_0__ = 0.;
	    s2_0__ = 0.;
	    v3_0__ = 0.;
	    s3_0__ = 0.;
	    v4_0__ = 0.;
	    s4_0__ = 0.;
	    v5_0__ = 0.;
	    s5_0__ = 0.;
	    v6_0__ = 0.;
	    s6_0__ = 0.;
	    for (j = 1; j <= 2000; ++j) {
		if (s_cmp(kt_ref(0, j), pt, (ftnlen)16, (ftnlen)16) == 0) {
		    v1_0__ = t1_ref(1, j);
		    s1_0__ = t1_ref(2, j);
		    v2_0__ = t2_ref(1, j);
		    s2_0__ = t2_ref(2, j);
		    v3_0__ = t3_ref(1, j);
		    s3_0__ = t3_ref(2, j);
		    v4_0__ = t4_ref(1, j);
		    s4_0__ = t4_ref(2, j);
		    v5_0__ = t5_ref(1, j);
		    s5_0__ = t5_ref(2, j);
		    v6_0__ = t6_ref(1, j);
		    s6_0__ = t6_ref(2, j);
		    goto L10;
		}
	    }
	    for (j = 1; j <= 2000; ++j) {
		if (s_cmp(kt_ref(0, j), pt0, (ftnlen)16, (ftnlen)16) == 0) {
		    v1_0__ = t1_ref(1, j);
		    s1_0__ = t1_ref(2, j);
		    v2_0__ = t2_ref(1, j);
		    s2_0__ = t2_ref(2, j);
		    v3_0__ = t3_ref(1, j);
		    s3_0__ = t3_ref(2, j);
		    v4_0__ = t4_ref(1, j);
		    s4_0__ = t4_ref(2, j);
		    v5_0__ = t5_ref(1, j);
		    s5_0__ = t5_ref(2, j);
		    v6_0__ = t6_ref(1, j);
		    s6_0__ = t6_ref(2, j);
		    goto L10;
		}
	    }
L10:

/*     find the torsion parameters for the final state */

	    i__2 = mutant_1.nmut;
	    for (j = 1; j <= i__2; ++j) {
		k = mutant_1.imut[j - 1];
		if (k == ia) {
		    ita = mutant_1.class1[j - 1];
		}
		if (k == ib) {
		    itb = mutant_1.class1[j - 1];
		}
		if (k == ic) {
		    itc = mutant_1.class1[j - 1];
		}
		if (k == id) {
		    itd = mutant_1.class1[j - 1];
		}
	    }
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    numeral_(&itc, pc, &size, (ftnlen)4);
	    numeral_(&itd, pd, &size, (ftnlen)4);
	    if (itb < itc) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (itc < itb) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (ita <= itd) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (itd < ita) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    }
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = zeros;
	    i__4[1] = 8, a__2[1] = pt + 4;
	    i__4[2] = 4, a__2[2] = zeros;
	    s_cat(pt0, a__2, i__4, &c__3, (ftnlen)16);
	    v1_1__ = 0.;
	    s1_1__ = 0.;
	    v2_1__ = 0.;
	    s2_1__ = 0.;
	    v3_1__ = 0.;
	    s3_1__ = 0.;
	    v4_1__ = 0.;
	    s4_1__ = 0.;
	    v5_1__ = 0.;
	    s5_1__ = 0.;
	    v6_1__ = 0.;
	    s6_1__ = 0.;
	    for (j = 1; j <= 2000; ++j) {
		if (s_cmp(kt_ref(0, j), pt, (ftnlen)16, (ftnlen)16) == 0) {
		    v1_1__ = t1_ref(1, j);
		    s1_1__ = t1_ref(2, j);
		    v2_1__ = t2_ref(1, j);
		    s2_1__ = t2_ref(2, j);
		    v3_1__ = t3_ref(1, j);
		    s3_1__ = t3_ref(2, j);
		    v4_1__ = t4_ref(1, j);
		    s4_1__ = t4_ref(2, j);
		    v5_1__ = t5_ref(1, j);
		    s5_1__ = t5_ref(2, j);
		    v6_1__ = t6_ref(1, j);
		    s6_1__ = t6_ref(2, j);
		    goto L20;
		}
	    }
	    for (j = 1; j <= 2000; ++j) {
		if (s_cmp(kt_ref(0, j), pt0, (ftnlen)16, (ftnlen)16) == 0) {
		    v1_1__ = t1_ref(1, j);
		    s1_1__ = t1_ref(2, j);
		    v2_1__ = t2_ref(1, j);
		    s2_1__ = t2_ref(2, j);
		    v3_1__ = t3_ref(1, j);
		    s3_1__ = t3_ref(2, j);
		    v4_1__ = t4_ref(1, j);
		    s4_1__ = t4_ref(2, j);
		    v5_1__ = t5_ref(1, j);
		    s5_1__ = t5_ref(2, j);
		    v6_1__ = t6_ref(1, j);
		    s6_1__ = t6_ref(2, j);
		    goto L20;
		}
	    }
L20:

/*     form the hybrid parameters for the current torsion */

	    if (s1_0__ == 0.) {
		s1_0__ = s1_1__;
	    }
	    if (s2_0__ == 0.) {
		s2_0__ = s2_1__;
	    }
	    if (s3_0__ == 0.) {
		s3_0__ = s3_1__;
	    }
	    if (s4_0__ == 0.) {
		s4_0__ = s4_1__;
	    }
	    if (s5_0__ == 0.) {
		s5_0__ = s5_1__;
	    }
	    if (s6_0__ == 0.) {
		s6_0__ = s6_1__;
	    }
	    if (s1_1__ == 0.) {
		s1_1__ = s1_0__;
	    }
	    if (s2_1__ == 0.) {
		s2_1__ = s2_0__;
	    }
	    if (s3_1__ == 0.) {
		s3_1__ = s3_0__;
	    }
	    if (s4_1__ == 0.) {
		s4_1__ = s4_0__;
	    }
	    if (s5_1__ == 0.) {
		s5_1__ = s5_0__;
	    }
	    if (s6_1__ == 0.) {
		s6_1__ = s6_0__;
	    }
	    tors1_ref(1, i__) = mutant_1.lambda * v1_1__ + (1. - 
		    mutant_1.lambda) * v1_0__;
	    tors1_ref(2, i__) = mutant_1.lambda * s1_1__ + (1. - 
		    mutant_1.lambda) * s1_0__;
	    angle = tors1_ref(2, i__) / 57.29577951308232088;
	    tors1_ref(3, i__) = cos(angle);
	    tors1_ref(4, i__) = sin(angle);
	    tors2_ref(1, i__) = mutant_1.lambda * v2_1__ + (1. - 
		    mutant_1.lambda) * v2_0__;
	    tors2_ref(2, i__) = mutant_1.lambda * s2_1__ + (1. - 
		    mutant_1.lambda) * s2_0__;
	    angle = tors2_ref(2, i__) / 57.29577951308232088;
	    tors2_ref(3, i__) = cos(angle);
	    tors2_ref(4, i__) = sin(angle);
	    tors3_ref(1, i__) = mutant_1.lambda * v3_1__ + (1. - 
		    mutant_1.lambda) * v3_0__;
	    tors3_ref(2, i__) = mutant_1.lambda * s3_1__ + (1. - 
		    mutant_1.lambda) * s3_0__;
	    angle = tors3_ref(2, i__) / 57.29577951308232088;
	    tors3_ref(3, i__) = cos(angle);
	    tors3_ref(4, i__) = sin(angle);
	    tors4_ref(1, i__) = mutant_1.lambda * v4_1__ + (1. - 
		    mutant_1.lambda) * v4_0__;
	    tors4_ref(2, i__) = mutant_1.lambda * s4_1__ + (1. - 
		    mutant_1.lambda) * s4_0__;
	    angle = tors4_ref(2, i__) / 57.29577951308232088;
	    tors4_ref(3, i__) = cos(angle);
	    tors4_ref(4, i__) = sin(angle);
	    tors5_ref(1, i__) = mutant_1.lambda * v5_1__ + (1. - 
		    mutant_1.lambda) * v5_0__;
	    tors5_ref(2, i__) = mutant_1.lambda * s5_1__ + (1. - 
		    mutant_1.lambda) * s5_0__;
	    angle = tors5_ref(2, i__) / 57.29577951308232088;
	    tors5_ref(3, i__) = cos(angle);
	    tors5_ref(4, i__) = sin(angle);
	    tors6_ref(1, i__) = mutant_1.lambda * v6_1__ + (1. - 
		    mutant_1.lambda) * v6_0__;
	    tors6_ref(2, i__) = mutant_1.lambda * s6_1__ + (1. - 
		    mutant_1.lambda) * s6_0__;
	    angle = tors6_ref(2, i__) / 57.29577951308232088;
	    tors6_ref(3, i__) = cos(angle);
	    tors6_ref(4, i__) = sin(angle);
	    if (inform_1.verbose) {
		if (header) {
		    header = FALSE_;
		    io___154.ciunit = iounit_1.iout;
		    s_wsfe(&io___154);
		    e_wsfe();
		}
		io___155.ciunit = iounit_1.iout;
		s_wsfe(&io___155);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&tors1_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		i__2 = i_dnnt(&tors1_ref(2, i__));
		do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&tors2_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		i__5 = i_dnnt(&tors2_ref(2, i__));
		do_fio(&c__1, (char *)&i__5, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&tors3_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		i__6 = i_dnnt(&tors3_ref(2, i__));
		do_fio(&c__1, (char *)&i__6, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&tors4_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		i__7 = i_dnnt(&tors4_ref(2, i__));
		do_fio(&c__1, (char *)&i__7, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&tors5_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		i__8 = i_dnnt(&tors5_ref(2, i__));
		do_fio(&c__1, (char *)&i__8, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&tors6_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		i__9 = i_dnnt(&tors6_ref(2, i__));
		do_fio(&c__1, (char *)&i__9, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	}
    }
    return 0;
} /* htors_ */

#undef itors_ref
#undef tors6_ref
#undef tors5_ref
#undef tors4_ref
#undef tors3_ref
#undef tors2_ref
#undef tors1_ref
#undef kt_ref
#undef t6_ref
#undef t5_ref
#undef t4_ref
#undef t3_ref
#undef t2_ref
#undef t1_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine hstrtor  --  hybrid stretch-torsion terms  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "hstrtor" constructs hybrid stretch-torsion parameters */
/*     given an initial state, final state and "lambda" value */


/* Subroutine */ int hstrtor_(void)
{
    /* Format strings */
    static char fmt_40[] = "(/,\002 Hybrid Stretch-Torsion Parameters :\002,"
	    "//,6x,\002Atom Numbers\002,13x,\002KST1\002,8x,\002KST2\002,8x"
	    ",\002KST3\002,/)";
    static char fmt_50[] = "(3x,4i5,3f12.3)";

    /* System generated locals */
    address a__1[4], a__2[3];
    integer i__1, i__2, i__3[4], i__4[3];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, ia, ib, ic, id;
    static char pa[4], pb[4], pc[4], pd[4], pt[16], pt0[16];
    static integer ita, itb, itc, itd;
    static doublereal kst0[3], kst1[3];
    static integer size;
    static char zeros[4];
    static logical header;
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    ;

    /* Fortran I/O blocks */
    static cilist io___178 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___179 = { 0, 0, 0, fmt_50, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define kbt_ref(a_0,a_1) &ksttor_1.kbt[(a_1)*16 + a_0 - 16]
#define ist_ref(a_1,a_2) strtor_1.ist[(a_2)*2 + a_1 - 3]
#define kst_ref(a_1,a_2) strtor_1.kst[(a_2)*3 + a_1 - 4]
#define btcon_ref(a_1,a_2) ksttor_1.btcon[(a_2)*3 + a_1 - 4]
#define itors_ref(a_1,a_2) tors_1.itors[(a_2)*4 + a_1 - 5]
#define bndlist_ref(a_1,a_2) atmlst_1.bndlist[(a_2)*8 + a_1 - 9]



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
/*     ##  atmlst.i  --  local geometry terms involving each atom  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     bndlist   list of the bond numbers involving each atom */
/*     anglist   list of the angle numbers centered on each atom */




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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  ksttor.i  --  forcefield parameters for stretch-torsions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxnbt   maximum number of stretch-torsion parameter entries */

/*     btcon    force constant parameters for stretch-torsion */
/*     kbt      string of atom classes for stretch-torsion terms */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  strtor.i  --  stretch-torsions in the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     kst       1-, 2- and 3-fold stretch-torsion force constants */
/*     nstrtor   total number of stretch-torsion interactions */
/*     ist       torsion and bond numbers used in stretch-torsion */




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




/*     assign hybrid parameters for the stretch-torsion sites */

    s_copy(zeros, "0000", (ftnlen)4, (ftnlen)4);
    header = TRUE_;
    i__1 = tors_1.ntors;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = itors_ref(1, i__);
	ib = itors_ref(2, i__);
	ic = itors_ref(3, i__);
	id = itors_ref(4, i__);
	if (mutant_1.mut[ia - 1] || mutant_1.mut[ib - 1] || mutant_1.mut[ic - 
		1] || mutant_1.mut[id - 1]) {
	    ita = atmtyp_1.class__[ia - 1];
	    itb = atmtyp_1.class__[ib - 1];
	    itc = atmtyp_1.class__[ic - 1];
	    itd = atmtyp_1.class__[id - 1];

/*     find the stretch-torsion parameters for the initial state */

	    i__2 = mutant_1.nmut;
	    for (j = 1; j <= i__2; ++j) {
		k = mutant_1.imut[j - 1];
		if (k == ia) {
		    ita = mutant_1.class0[j - 1];
		}
		if (k == ib) {
		    itb = mutant_1.class0[j - 1];
		}
		if (k == ic) {
		    itc = mutant_1.class0[j - 1];
		}
		if (k == id) {
		    itd = mutant_1.class0[j - 1];
		}
	    }
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    numeral_(&itc, pc, &size, (ftnlen)4);
	    numeral_(&itd, pd, &size, (ftnlen)4);
	    if (itb < itc) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (itc < itb) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (ita <= itd) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (itd < ita) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    }
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = zeros;
	    i__4[1] = 8, a__2[1] = pt + 4;
	    i__4[2] = 4, a__2[2] = zeros;
	    s_cat(pt0, a__2, i__4, &c__3, (ftnlen)16);
	    for (k = 1; k <= 3; ++k) {
		kst0[k - 1] = 0.;
	    }
	    for (j = 1; j <= 500; ++j) {
		if (s_cmp(kbt_ref(0, j), pt, (ftnlen)16, (ftnlen)16) == 0) {
		    for (k = 1; k <= 3; ++k) {
			kst0[k - 1] = btcon_ref(k, j);
		    }
		    goto L10;
		}
	    }
	    for (j = 1; j <= 500; ++j) {
		if (s_cmp(kbt_ref(0, j), pt0, (ftnlen)16, (ftnlen)16) == 0) {
		    for (k = 1; k <= 3; ++k) {
			kst0[k - 1] = btcon_ref(k, j);
		    }
		    goto L10;
		}
	    }
L10:

/*     find the stretch-torsion parameters for the final state */

	    i__2 = mutant_1.nmut;
	    for (j = 1; j <= i__2; ++j) {
		k = mutant_1.imut[j - 1];
		if (k == ia) {
		    ita = mutant_1.class0[j - 1];
		}
		if (k == ib) {
		    itb = mutant_1.class0[j - 1];
		}
		if (k == ic) {
		    itc = mutant_1.class0[j - 1];
		}
		if (k == id) {
		    itd = mutant_1.class0[j - 1];
		}
	    }
	    size = 4;
	    numeral_(&ita, pa, &size, (ftnlen)4);
	    numeral_(&itb, pb, &size, (ftnlen)4);
	    numeral_(&itc, pc, &size, (ftnlen)4);
	    numeral_(&itd, pd, &size, (ftnlen)4);
	    if (itb < itc) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (itc < itb) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (ita <= itd) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		i__3[2] = 4, a__1[2] = pc;
		i__3[3] = 4, a__1[3] = pd;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    } else if (itd < ita) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pd;
		i__3[1] = 4, a__1[1] = pc;
		i__3[2] = 4, a__1[2] = pb;
		i__3[3] = 4, a__1[3] = pa;
		s_cat(pt, a__1, i__3, &c__4, (ftnlen)16);
	    }
/* Writing concatenation */
	    i__4[0] = 4, a__2[0] = zeros;
	    i__4[1] = 8, a__2[1] = pt + 4;
	    i__4[2] = 4, a__2[2] = zeros;
	    s_cat(pt0, a__2, i__4, &c__3, (ftnlen)16);
	    for (k = 1; k <= 3; ++k) {
		kst1[k - 1] = 0.;
	    }
	    for (j = 1; j <= 500; ++j) {
		if (s_cmp(kbt_ref(0, j), pt, (ftnlen)16, (ftnlen)16) == 0) {
		    for (k = 1; k <= 3; ++k) {
			kst1[k - 1] = btcon_ref(k, j);
		    }
		    goto L20;
		}
	    }
	    for (j = 1; j <= 500; ++j) {
		if (s_cmp(kbt_ref(0, j), pt0, (ftnlen)16, (ftnlen)16) == 0) {
		    for (k = 1; k <= 3; ++k) {
			kst1[k - 1] = btcon_ref(k, j);
		    }
		    goto L20;
		}
	    }
L20:

/*     form hybrid parameters for the current stretch-torsion */

	    for (j = 1; j <= 3; ++j) {
		kst_ref(j, i__) = mutant_1.lambda * kst1[j - 1] + (1. - 
			mutant_1.lambda) * kst0[j - 1];
	    }
	    if (kst_ref(1, i__) == 0. && kst_ref(2, i__) == 0. && kst_ref(3, 
		    i__) == 0.) {
		if (ist_ref(1, i__) != 0) {
		    --strtor_1.nstrtor;
		    ist_ref(1, i__) = 0;
		}
	    } else {
		if (ist_ref(1, i__) != i__) {
		    ++strtor_1.nstrtor;
		    ist_ref(1, i__) = i__;
		    i__2 = couple_1.n12[ib - 1];
		    for (j = 1; j <= i__2; ++j) {
			if (i12_ref(j, ib) == ic) {
			    ist_ref(2, i__) = bndlist_ref(j, ib);
			    goto L30;
			}
		    }
L30:
		    ;
		}
		if (inform_1.verbose) {
		    if (header) {
			header = FALSE_;
			io___178.ciunit = iounit_1.iout;
			s_wsfe(&io___178);
			e_wsfe();
		    }
		    io___179.ciunit = iounit_1.iout;
		    s_wsfe(&io___179);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		    for (j = 1; j <= 3; ++j) {
			do_fio(&c__1, (char *)&kst_ref(j, i__), (ftnlen)
				sizeof(doublereal));
		    }
		    e_wsfe();
		}
	    }
	}
    }
    return 0;
} /* hstrtor_ */

#undef bndlist_ref
#undef itors_ref
#undef btcon_ref
#undef kst_ref
#undef ist_ref
#undef kbt_ref
#undef i12_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine hvdw  --  hybrid van der Waals parameters  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "hvdw" constructs hybrid van der Waals  parameters given */
/*     an initial state, final state and "lambda" value */


/* Subroutine */ int hvdw_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Hybrid van der Waals Parameters :\002,//"
	    ",7x,\002Atom Number\002,5x,\002Size\002,6x,\002Epsilon\002,/)";
    static char fmt_20[] = "(6x,i8,f14.4,f12.4)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k;
    static doublereal ep, rd;
    static integer it, it0, it1;
    static doublereal srad[5000], seps[5000], srad4[5000], seps4[5000];
    static logical header;
    static doublereal radius;

    /* Fortran I/O blocks */
    static cilist io___193 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___195 = { 0, 0, 0, fmt_20, 0 };



#define epsilon4_ref(a_1,a_2) vdw_1.epsilon4[(a_2)*1000 + a_1 - 1001]
#define radmin_ref(a_1,a_2) vdw_1.radmin[(a_2)*1000 + a_1 - 1001]
#define radmin4_ref(a_1,a_2) vdw_1.radmin4[(a_2)*1000 + a_1 - 1001]
#define epsilon_ref(a_1,a_2) vdw_1.epsilon[(a_2)*1000 + a_1 - 1001]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  vdw.i  --  van der Waals parameters for current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     radmin     minimum energy distance for each atom class pair */
/*     epsilon    well depth parameter for each atom class pair */
/*     radmin4    minimum energy distance for 1-4 interaction pairs */
/*     epsilon4   well depth parameter for 1-4 interaction pairs */
/*     radhbnd    minimum energy distance for hydrogen bonding pairs */
/*     epshbnd    well depth parameter for hydrogen bonding pairs */
/*     kred       value of reduction factor parameter for each atom */
/*     ired       attached atom from which reduction factor is applied */
/*     nvdw       total number van der Waals active sites in the system */
/*     ivdw       number of the atom for each van der Waals active site */
/*     jvdw       type or class index into vdw parameters for each atom */
/*     nvt        number of distinct vdw types/classes in the system */
/*     ivt        type/class index for each distinct vdw type or class */
/*     jvt        frequency of each vdw type or class in the system */




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




/*     assign the hybrid van der Waals parameters */

    i__1 = mutant_1.nmut;
    for (j = 1; j <= i__1; ++j) {
	i__ = mutant_1.imut[j - 1];
	if (s_cmp(vdwpot_1.vdwindex, "TYPE", (ftnlen)5, (ftnlen)4) == 0) {
	    it = atoms_1.type__[i__ - 1];
	    it0 = mutant_1.type0[j - 1];
	    it1 = mutant_1.type1[j - 1];
	} else {
	    it = atmtyp_1.class__[i__ - 1];
	    it0 = mutant_1.class0[j - 1];
	    it1 = mutant_1.class1[j - 1];
	}
	mutant_1.vlambda = mutant_1.lambda;
	kvdws_1.rad[it - 1] = mutant_1.vlambda * kvdws_1.rad[it1 - 1] + (1. - 
		mutant_1.vlambda) * kvdws_1.rad[it0 - 1];
	kvdws_1.eps[it - 1] = mutant_1.vlambda * kvdws_1.eps[it1 - 1] + (1. - 
		mutant_1.vlambda) * kvdws_1.eps[it0 - 1];
    }

/*     get the square roots of the vdw radii and well depths */

    for (i__ = 1; i__ <= 1000; ++i__) {
	srad[i__ - 1] = sqrt(kvdws_1.rad[i__ - 1]);
	seps[i__ - 1] = sqrt(kvdws_1.eps[i__ - 1]);
    }

/*     use combination rules to set pairwise vdw radii sums */

    i__1 = mutant_1.nmut;
    for (j = 1; j <= i__1; ++j) {
	i__ = mutant_1.imut[j - 1];
	it = atmtyp_1.class__[i__ - 1];
	if (s_cmp(vdwpot_1.vdwindex, "TYPE", (ftnlen)5, (ftnlen)4) == 0) {
	    it = atoms_1.type__[i__ - 1];
	}
	for (k = 1; k <= 1000; ++k) {
	    if (kvdws_1.rad[it - 1] == 0. && kvdws_1.rad[k - 1] == 0.) {
		rd = 0.;
	    } else if (s_cmp(vdwpot_1.radrule, "ARITHMETIC", (ftnlen)10, (
		    ftnlen)10) == 0) {
		rd = kvdws_1.rad[it - 1] + kvdws_1.rad[k - 1];
	    } else if (s_cmp(vdwpot_1.radrule, "GEOMETRIC", (ftnlen)9, (
		    ftnlen)9) == 0) {
		rd = srad[it - 1] * srad[k - 1] * 2.;
	    } else if (s_cmp(vdwpot_1.radrule, "CUBIC-MEAN", (ftnlen)10, (
		    ftnlen)10) == 0) {
/* Computing 3rd power */
		d__1 = kvdws_1.rad[it - 1];
/* Computing 3rd power */
		d__2 = kvdws_1.rad[k - 1];
/* Computing 2nd power */
		d__3 = kvdws_1.rad[it - 1];
/* Computing 2nd power */
		d__4 = kvdws_1.rad[k - 1];
		rd = (d__1 * (d__1 * d__1) + d__2 * (d__2 * d__2)) * 2. / (
			d__3 * d__3 + d__4 * d__4);
	    } else {
		rd = kvdws_1.rad[it - 1] + kvdws_1.rad[k - 1];
	    }
	    radmin_ref(it, k) = rd;
	    radmin_ref(k, it) = rd;
	}
    }

/*     use combination rules to set pairwise well depths */

    i__1 = mutant_1.nmut;
    for (j = 1; j <= i__1; ++j) {
	i__ = mutant_1.imut[j - 1];
	it = atmtyp_1.class__[i__ - 1];
	if (s_cmp(vdwpot_1.vdwindex, "TYPE", (ftnlen)5, (ftnlen)4) == 0) {
	    it = atoms_1.type__[i__ - 1];
	}
	for (k = 1; k <= 1000; ++k) {
	    if (kvdws_1.eps[it - 1] == 0. && kvdws_1.eps[k - 1] == 0.) {
		ep = 0.;
	    } else if (s_cmp(vdwpot_1.epsrule, "ARITHMETIC", (ftnlen)10, (
		    ftnlen)10) == 0) {
		ep = (kvdws_1.eps[it - 1] + kvdws_1.eps[k - 1]) * .5;
	    } else if (s_cmp(vdwpot_1.epsrule, "GEOMETRIC", (ftnlen)9, (
		    ftnlen)9) == 0) {
		ep = seps[it - 1] * seps[k - 1];
	    } else if (s_cmp(vdwpot_1.epsrule, "HARMONIC", (ftnlen)8, (ftnlen)
		    8) == 0) {
		ep = kvdws_1.eps[it - 1] * kvdws_1.eps[k - 1] * 2. / (
			kvdws_1.eps[it - 1] + kvdws_1.eps[k - 1]);
	    } else if (s_cmp(vdwpot_1.epsrule, "HHG", (ftnlen)3, (ftnlen)3) ==
		     0) {
/* Computing 2nd power */
		d__1 = seps[it - 1] + seps[k - 1];
		ep = kvdws_1.eps[it - 1] * kvdws_1.eps[k - 1] * 4. / (d__1 * 
			d__1);
	    } else {
		ep = seps[it - 1] * seps[k - 1];
	    }
	    epsilon_ref(it, k) = ep;
	    epsilon_ref(k, it) = ep;
	}
    }

/*     use combination rules to set pairwise 1-4 vdw radii sums */

    i__1 = mutant_1.nmut;
    for (j = 1; j <= i__1; ++j) {
	i__ = mutant_1.imut[j - 1];
	it = atmtyp_1.class__[i__ - 1];
	if (s_cmp(vdwpot_1.vdwindex, "TYPE", (ftnlen)5, (ftnlen)4) == 0) {
	    it = atoms_1.type__[i__ - 1];
	}
	for (k = 1; k <= 1000; ++k) {
	    if (kvdws_1.rad4[it - 1] == 0. && kvdws_1.rad4[k - 1] == 0.) {
		rd = 0.;
	    } else if (s_cmp(vdwpot_1.radrule, "ARITHMETIC", (ftnlen)10, (
		    ftnlen)10) == 0) {
		rd = kvdws_1.rad4[it - 1] + kvdws_1.rad4[k - 1];
	    } else if (s_cmp(vdwpot_1.radrule, "GEOMETRIC", (ftnlen)9, (
		    ftnlen)9) == 0) {
		rd = srad4[it - 1] * srad4[k - 1] * 2.;
	    } else if (s_cmp(vdwpot_1.radrule, "CUBIC-MEAN", (ftnlen)10, (
		    ftnlen)10) == 0) {
/* Computing 3rd power */
		d__1 = kvdws_1.rad4[it - 1];
/* Computing 3rd power */
		d__2 = kvdws_1.rad4[k - 1];
/* Computing 2nd power */
		d__3 = kvdws_1.rad4[it - 1];
/* Computing 2nd power */
		d__4 = kvdws_1.rad4[k - 1];
		rd = (d__1 * (d__1 * d__1) + d__2 * (d__2 * d__2)) * 2. / (
			d__3 * d__3 + d__4 * d__4);
	    } else {
		rd = kvdws_1.rad4[it - 1] + kvdws_1.rad4[k - 1];
	    }
	    radmin4_ref(it, k) = rd;
	    radmin4_ref(k, it) = rd;
	}
    }

/*     use combination rules to set pairwise 1-4 well depths */

    i__1 = mutant_1.nmut;
    for (j = 1; j <= i__1; ++j) {
	i__ = mutant_1.imut[j - 1];
	it = atmtyp_1.class__[i__ - 1];
	if (s_cmp(vdwpot_1.vdwindex, "TYPE", (ftnlen)5, (ftnlen)4) == 0) {
	    it = atoms_1.type__[i__ - 1];
	}
	for (k = 1; k <= 1000; ++k) {
	    if (kvdws_1.eps4[it - 1] == 0. && kvdws_1.eps4[k - 1] == 0.) {
		ep = 0.;
	    } else if (s_cmp(vdwpot_1.epsrule, "ARITHMETIC", (ftnlen)10, (
		    ftnlen)10) == 0) {
		ep = (kvdws_1.eps4[it - 1] + kvdws_1.eps4[k - 1]) * .5;
	    } else if (s_cmp(vdwpot_1.epsrule, "GEOMETRIC", (ftnlen)9, (
		    ftnlen)9) == 0) {
		ep = seps4[it - 1] * seps4[k - 1];
	    } else if (s_cmp(vdwpot_1.epsrule, "HARMONIC", (ftnlen)8, (ftnlen)
		    8) == 0) {
		ep = kvdws_1.eps4[it - 1] * kvdws_1.eps4[k - 1] * 2. / (
			kvdws_1.eps4[it - 1] + kvdws_1.eps4[k - 1]);
	    } else if (s_cmp(vdwpot_1.epsrule, "HHG", (ftnlen)3, (ftnlen)3) ==
		     0) {
/* Computing 2nd power */
		d__1 = seps4[it - 1] + seps4[k - 1];
		ep = kvdws_1.eps4[it - 1] * kvdws_1.eps4[k - 1] * 4. / (d__1 *
			 d__1);
	    } else {
		ep = seps4[it - 1] * seps4[k - 1];
	    }
	    epsilon4_ref(it, k) = ep;
	    epsilon4_ref(k, it) = ep;
	}
    }

/*     print the van der Waals parameters for hybrid atoms */

    header = TRUE_;
    i__1 = mutant_1.nmut;
    for (j = 1; j <= i__1; ++j) {
	i__ = mutant_1.imut[j - 1];
	it = atmtyp_1.class__[i__ - 1];
	if (s_cmp(vdwpot_1.vdwindex, "TYPE", (ftnlen)5, (ftnlen)4) == 0) {
	    it = atoms_1.type__[i__ - 1];
	}
	if (inform_1.verbose) {
	    if (header) {
		header = FALSE_;
		io___193.ciunit = iounit_1.iout;
		s_wsfe(&io___193);
		e_wsfe();
	    }
	    radius = kvdws_1.rad[it - 1];
	    if (s_cmp(vdwpot_1.radsiz, "DIAMETER", (ftnlen)8, (ftnlen)8) == 0)
		     {
		radius *= 2.;
	    }
	    if (s_cmp(vdwpot_1.radtyp, "SIGMA", (ftnlen)5, (ftnlen)5) == 0) {
		radius /= 1.122462048309372981;
	    }
	    io___195.ciunit = iounit_1.iout;
	    s_wsfe(&io___195);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&radius, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&kvdws_1.eps[it - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* hvdw_ */

#undef epsilon_ref
#undef radmin4_ref
#undef radmin_ref
#undef epsilon4_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine hcharge  --  find hybrid charge parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "hcharge" constructs hybrid charge interaction parameters */
/*     given an initial state, final state and "lambda" value */


/* Subroutine */ int hcharge_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Hybrid Atomic Partial Charge Parameter"
	    "s :\002,//,7x,\002Atom Number\002,7x,\002Charge\002,/)";
    static char fmt_30[] = "(6x,i8,5x,f12.3)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, it, it0, it1;
    static doublereal chg0, chg1, hchg;
    static logical used, header;

    /* Fortran I/O blocks */
    static cilist io___207 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___208 = { 0, 0, 0, fmt_30, 0 };




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
/*     ##  kchrge.i  --  forcefield parameters for partial charges  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     chg   partial charge parameters for each atom type */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




/*     assign the hybrid parameters for atomic charges */

    header = TRUE_;
    i__1 = mutant_1.nmut;
    for (j = 1; j <= i__1; ++j) {
	used = FALSE_;
	i__ = mutant_1.imut[j - 1];
	it = atoms_1.type__[i__ - 1];
	it0 = mutant_1.type0[j - 1];
	it1 = mutant_1.type1[j - 1];
	chg0 = 0.;
	chg1 = 0.;
	if (it0 != 0) {
	    chg0 = kchrge_1.chg[it0 - 1];
	}
	if (it1 != 0) {
	    chg1 = kchrge_1.chg[it1 - 1];
	}
	mutant_1.clambda = mutant_1.lambda;
	hchg = mutant_1.clambda * chg1 + (1. - mutant_1.clambda) * chg0;
	i__2 = charge_1.nion;
	for (k = 1; k <= i__2; ++k) {
	    if (charge_1.iion[k - 1] == i__) {
		used = TRUE_;
		charge_1.pchg[k - 1] = hchg;
		goto L10;
	    }
	}
	if (chg0 != 0. || chg1 != 0.) {
	    used = TRUE_;
	    ++charge_1.nion;
	    charge_1.iion[charge_1.nion - 1] = i__;
	    charge_1.kion[charge_1.nion - 1] = i__;
	    charge_1.pchg[charge_1.nion - 1] = hchg;
	}
L10:
	if (inform_1.verbose && used) {
	    if (header) {
		header = FALSE_;
		io___207.ciunit = iounit_1.iout;
		s_wsfe(&io___207);
		e_wsfe();
	    }
	    io___208.ciunit = iounit_1.iout;
	    s_wsfe(&io___208);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&hchg, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* hcharge_ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine hdipole  --  find hybrid dipole parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "hdipole" constructs hybrid dipole interaction parameters */
/*     given an initial state, final state and "lambda" value */


/* Subroutine */ int hdipole_(void)
{
    /* Format strings */
    static char fmt_40[] = "(/,\002 Hybrid Bond Dipole Moment Parameters "
	    ":\002,//,6x,\002Atom Numbers\002,7x,\002Moment\002,7x,\002Positi"
	    "on\002,/)";
    static char fmt_50[] = "(6x,2i5,2f15.3)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen), s_cat(char *,
	     char **, integer *, integer *, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, ia, ib;
    static char pa[4], pb[4], pt[8];
    static integer ita, itb;
    static doublereal dpl0, dpl1, pos0, pos1, hdpl;
    static logical used;
    static doublereal hpos;
    static integer size;
    static char blank[8];
    static logical header;
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    ;

    /* Fortran I/O blocks */
    static cilist io___229 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___230 = { 0, 0, 0, fmt_50, 0 };



#define kd_ref(a_0,a_1) &kdipol_1.kd[(a_1)*8 + a_0 - 8]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define idpl_ref(a_1,a_2) dipole_1.idpl[(a_2)*2 + a_1 - 3]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  kdipol.i  --  forcefield parameters for bond dipoles  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxnd    maximum number of bond dipole parameter entries */
/*     maxnd5   maximum number of 5-membered ring dipole entries */
/*     maxnd4   maximum number of 4-membered ring dipole entries */
/*     maxnd3   maximum number of 3-membered ring dipole entries */

/*     dpl      dipole moment parameters for bond dipoles */
/*     dpl5     dipole moment parameters for 5-ring dipoles */
/*     dpl4     dipole moment parameters for 4-ring dipoles */
/*     dpl3     dipole moment parameters for 3-ring dipoles */
/*     pos      dipole position parameters for bond dipoles */
/*     pos5     dipole position parameters for 5-ring dipoles */
/*     pos4     dipole position parameters for 4-ring dipoles */
/*     pos3     dipole position parameters for 3-ring dipoles */
/*     kd       string of atom classes for bond dipoles */
/*     kd5      string of atom classes for 5-ring dipoles */
/*     kd4      string of atom classes for 4-ring dipoles */
/*     kd3      string of atom classes for 3-ring dipoles */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




/*     assign the hybrid parameters for bond dipoles */

    s_copy(blank, "      ", (ftnlen)8, (ftnlen)6);
    header = TRUE_;
    i__1 = bond_1.nbond;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ibnd_ref(1, i__);
	ib = ibnd_ref(2, i__);
	if (mutant_1.mut[ia - 1] || mutant_1.mut[ib - 1]) {
	    ita = atoms_1.type__[ia - 1];
	    itb = atoms_1.type__[ib - 1];

/*     find the dipole parameters for the initial state */

	    i__2 = mutant_1.nmut;
	    for (j = 1; j <= i__2; ++j) {
		k = mutant_1.imut[j - 1];
		if (k == ia) {
		    ita = mutant_1.type0[j - 1];
		}
		if (k == ib) {
		    itb = mutant_1.type0[j - 1];
		}
	    }
	    dpl0 = 0.;
	    pos0 = .5;
	    if (ita != 0 && itb != 0) {
		size = 4;
		numeral_(&ita, pa, &size, (ftnlen)4);
		numeral_(&itb, pb, &size, (ftnlen)4);
		if (ita <= itb) {
/* Writing concatenation */
		    i__3[0] = 4, a__1[0] = pa;
		    i__3[1] = 4, a__1[1] = pb;
		    s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
		} else {
/* Writing concatenation */
		    i__3[0] = 4, a__1[0] = pb;
		    i__3[1] = 4, a__1[1] = pa;
		    s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
		}
		for (j = 1; j <= 1000; ++j) {
		    if (s_cmp(kd_ref(0, j), blank, (ftnlen)8, (ftnlen)8) == 0)
			     {
			goto L10;
		    }
		    if (s_cmp(kd_ref(0, j), pt, (ftnlen)8, (ftnlen)8) == 0) {
			if (ita <= itb) {
			    dpl0 = dipole_1.bdpl[j - 1];
			    pos0 = dipole_1.sdpl[j - 1];
			} else {
			    dpl0 = -dipole_1.bdpl[j - 1];
			    pos0 = 1. - dipole_1.sdpl[j - 1];
			}
		    }
		}
L10:
		;
	    }

/*     find the dipole parameters for the final state */

	    i__2 = mutant_1.nmut;
	    for (j = 1; j <= i__2; ++j) {
		k = mutant_1.imut[j - 1];
		if (k == ia) {
		    ita = mutant_1.type1[j - 1];
		}
		if (k == ib) {
		    itb = mutant_1.type1[j - 1];
		}
	    }
	    dpl1 = 0.;
	    pos1 = .5;
	    if (ita != 0 && itb != 0) {
		size = 4;
		numeral_(&ita, pa, &size, (ftnlen)4);
		numeral_(&itb, pb, &size, (ftnlen)4);
		if (ita <= itb) {
/* Writing concatenation */
		    i__3[0] = 4, a__1[0] = pa;
		    i__3[1] = 4, a__1[1] = pb;
		    s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
		} else {
/* Writing concatenation */
		    i__3[0] = 4, a__1[0] = pb;
		    i__3[1] = 4, a__1[1] = pa;
		    s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
		}
		for (j = 1; j <= 1000; ++j) {
		    if (s_cmp(kd_ref(0, j), blank, (ftnlen)8, (ftnlen)8) == 0)
			     {
			goto L20;
		    }
		    if (s_cmp(kd_ref(0, j), pt, (ftnlen)8, (ftnlen)8) == 0) {
			if (ita <= itb) {
			    dpl1 = dipole_1.bdpl[j - 1];
			    pos1 = dipole_1.sdpl[j - 1];
			} else {
			    dpl1 = -dipole_1.bdpl[j - 1];
			    pos1 = 1. - dipole_1.sdpl[j - 1];
			}
		    }
		}
L20:
		;
	    }

/*     form the hybrid parameters for the current dipole */

	    mutant_1.dlambda = mutant_1.lambda;
	    hdpl = mutant_1.dlambda * dpl1 + (1. - mutant_1.dlambda) * dpl0;
	    hpos = mutant_1.dlambda * pos1 + (1. - mutant_1.dlambda) * pos0;
	    used = FALSE_;
	    i__2 = dipole_1.ndipole;
	    for (j = 1; j <= i__2; ++j) {
		if (idpl_ref(1, j) == ia && idpl_ref(2, j) == ib || idpl_ref(
			1, j) == ib && idpl_ref(2, j) == ia) {
		    idpl_ref(1, j) = ia;
		    idpl_ref(2, j) = ib;
		    dipole_1.bdpl[j - 1] = hdpl;
		    dipole_1.sdpl[j - 1] = hpos;
		    used = TRUE_;
		    goto L30;
		}
	    }
	    if (hdpl != 0.) {
		++dipole_1.ndipole;
		idpl_ref(1, dipole_1.ndipole) = ia;
		idpl_ref(2, dipole_1.ndipole) = ib;
		dipole_1.bdpl[dipole_1.ndipole - 1] = hdpl;
		dipole_1.sdpl[dipole_1.ndipole - 1] = hpos;
		used = TRUE_;
	    }
L30:
	    if (inform_1.verbose && used) {
		if (header) {
		    header = FALSE_;
		    io___229.ciunit = iounit_1.iout;
		    s_wsfe(&io___229);
		    e_wsfe();
		}
		io___230.ciunit = iounit_1.iout;
		s_wsfe(&io___230);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&hdpl, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&hpos, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
    }
    return 0;
} /* hdipole_ */

#undef idpl_ref
#undef ibnd_ref
#undef kd_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine hmpole  --  find hybrid multipole parameters  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "hmpole" constructs hybrid multipole interaction parameters */
/*     given an initial state, final state and "lambda" value */


/* Subroutine */ int hmpole_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Hybrid Atomic Multipole Parameters :\002"
	    ",//,7x,\002Atom Number\002,7x,\002Charge\002,/)";
    static char fmt_30[] = "(6x,i8,5x,f12.3)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, m, it, it0, it1;
    static doublereal mpl0[13], mpl1[13], hmpl[13];
    static logical used, header;

    /* Fortran I/O blocks */
    static cilist io___243 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___244 = { 0, 0, 0, fmt_30, 0 };



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
/*     ##  kchrge.i  --  forcefield parameters for partial charges  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     chg   partial charge parameters for each atom type */




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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




/*     assign the hybrid parameters for atomic multipoles */

    header = TRUE_;
    i__1 = mutant_1.nmut;
    for (j = 1; j <= i__1; ++j) {
	used = FALSE_;
	i__ = mutant_1.imut[j - 1];
	it = atoms_1.type__[i__ - 1];
	it0 = mutant_1.type0[j - 1];
	it1 = mutant_1.type1[j - 1];
	for (m = 1; m <= 13; ++m) {
	    mpl0[m - 1] = 0.;
	    mpl1[m - 1] = 0.;
	}
	if (it0 != 0) {
	    for (m = 1; m <= 13; ++m) {
		mpl0[m - 1] = pole_ref(m, it0);
	    }
	}
	if (it1 != 0) {
	    for (m = 1; m <= 13; ++m) {
		mpl1[m - 1] = pole_ref(m, it1);
	    }
	}
	mutant_1.mlambda = mutant_1.lambda;
	for (k = 1; k <= 13; ++k) {
	    hmpl[k - 1] = mutant_1.mlambda * mpl1[k - 1] + (1. - 
		    mutant_1.mlambda) * mpl0[k - 1];
	}
	i__2 = mpole_1.npole;
	for (k = 1; k <= i__2; ++k) {
	    if (mpole_1.ipole[k - 1] == i__) {
		used = TRUE_;
		for (m = 1; m <= 13; ++m) {
		    pole_ref(m, k) = hmpl[m - 1];
		}
		goto L10;
	    }
	}
	if (mpl0[0] != 0. || mpl1[0] != 0.) {
	    used = TRUE_;
	    ++mpole_1.npole;
	    mpole_1.ipole[mpole_1.npole - 1] = i__;
	    for (m = 1; m <= 13; ++m) {
		pole_ref(m, mpole_1.npole) = hmpl[m - 1];
	    }
	}
L10:
	if (inform_1.verbose && used) {
	    if (header) {
		header = FALSE_;
		io___243.ciunit = iounit_1.iout;
		s_wsfe(&io___243);
		e_wsfe();
	    }
	    io___244.ciunit = iounit_1.iout;
	    s_wsfe(&io___244);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&hmpl[0], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* hmpole_ */

#undef pole_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine hpolar  --  find hybrid polarization parameters  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "hpolar" constructs hybrid atomic polarizability parameters */
/*     given an initial state, final state and "lambda" value */


/* Subroutine */ int hpolar_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Hybrid Dipole Polarizability Parameter"
	    "s :\002,//,7x,\002Atom Number\002,7x,\002Alpha\002,8x,\002Dam"
	    "p\002,/)";
    static char fmt_30[] = "(6x,i8,5x,2f12.3)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, it, it0, it1;
    static logical used;
    static doublereal athl0, athl1, polr0, polr1, hathl, hpolr;
    static logical header;

    /* Fortran I/O blocks */
    static cilist io___259 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___260 = { 0, 0, 0, fmt_30, 0 };




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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  kpolr.i  --  forcefield parameters for polarizability  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     polr   dipole polarizability parameters for each atom type */
/*     athl   Thole polarizability damping value for each atom type */
/*     pgrp   connected types in polarization group of each atom type */




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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




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




/*     assign the hybrid parameters for atomic polarizabilities */

    header = TRUE_;
    i__1 = mutant_1.nmut;
    for (j = 1; j <= i__1; ++j) {
	used = FALSE_;
	i__ = mutant_1.imut[j - 1];
	it = atoms_1.type__[i__ - 1];
	it0 = mutant_1.type0[j - 1];
	it1 = mutant_1.type1[j - 1];
	polr0 = 0.;
	athl0 = 0.;
	polr1 = 0.;
	athl1 = 0.;
	if (it0 != 0) {
	    polr0 = kpolr_1.polr[it0 - 1];
	    athl0 = kpolr_1.athl[it0 - 1];
	}
	if (it1 != 0) {
	    polr1 = kpolr_1.polr[it1 - 1];
	    athl1 = kpolr_1.athl[it1 - 1];
	}
	mutant_1.plambda = mutant_1.lambda;
	hpolr = mutant_1.plambda * polr1 + (1. - mutant_1.plambda) * polr0;
	hathl = mutant_1.plambda * athl1 + (1. - mutant_1.plambda) * athl0;
	i__2 = mpole_1.npole;
	for (k = 1; k <= i__2; ++k) {
	    if (mpole_1.ipole[k - 1] == i__) {
		used = TRUE_;
		polar_1.polarity[k - 1] = hpolr;
		polar_1.thole[k - 1] = hathl;
		goto L10;
	    }
	}
	if (polr0 != 0. || polr1 != 0.) {
	    used = TRUE_;
	    ++mpole_1.npole;
	    mpole_1.ipole[mpole_1.npole - 1] = i__;
	    ++polar_1.npolar;
	    polar_1.polarity[mpole_1.npole - 1] = hpolr;
	    polar_1.thole[mpole_1.npole - 1] = hathl;
	}
L10:
	if (inform_1.verbose && used) {
	    if (header) {
		header = FALSE_;
		io___259.ciunit = iounit_1.iout;
		s_wsfe(&io___259);
		e_wsfe();
	    }
	    io___260.ciunit = iounit_1.iout;
	    s_wsfe(&io___260);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&hpolr, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&hathl, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* hpolar_ */

