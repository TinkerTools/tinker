/* alchemy.f -- translated by f2c (version 20050501).
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
    doublereal aesum[25000], aeb[25000], aea[25000], aeba[25000], aeub[25000],
	     aeaa[25000], aeopb[25000], aeopd[25000], aeid[25000], aeit[25000]
	    , aet[25000], aept[25000], aebt[25000], aett[25000], aev[25000], 
	    aec[25000], aecd[25000], aed[25000], aem[25000], aep[25000], aer[
	    25000], aes[25000], aelf[25000], aeg[25000], aex[25000];
} analyz_;

#define analyz_1 analyz_

struct {
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

struct {
    doublereal esum, eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, 
	    ebt, ett, ev, ec, ecd, ed, em, ep, er, es, elf, eg, ex;
} energi_;

#define energi_1 energi_

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
    doublereal weight[5000];
    integer atmcls[5000], atmnum[5000], ligand[5000];
    char symbol[15000], describe[120000];
} katoms_;

#define katoms_1 katoms_

struct {
    doublereal lambda, vlambda, clambda, dlambda, mlambda, plambda;
    integer nmut, imut[25000], type0[25000], type1[25000], class0[25000], 
	    class1[25000];
    logical mut[25000];
} mutant_;

#define mutant_1 mutant_

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
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;



/*     ############################################################### */
/*     ##  COPYRIGHT (C) 1991 by Shawn Huston & Jay William Ponder  ## */
/*     ##                    All Rights Reserved                    ## */
/*     ############################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  program alchemy  --  perform free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "alchemy" computes the free energy difference corresponding */
/*     to a small perturbation by Boltzmann weighting the potential */
/*     energy difference over a number of sample states; current */
/*     version (incorrectly) considers the charge energy to be */
/*     intermolecular in finding the perturbation energies */

/*     variables and parameters : */

/*     nlamb   number of lambda values for energy computation */
/*     delta   step size for the perturbation in lambda */
/*     nrun    number of steps over which to calculate averages */

/*     deplus       energy change for postive delta lambda */
/*     deminus      energy change for negative delta lambda */
/*     sdep,sdem    accumulated energy changes */
/*     adep,adem    (running) average energy changes */
/*     adapb,adamb  average over block averages of free energy changes */
/*     bdep,bdem    (block) accumulated energy changes */
/*     badep,badem  (block) average energy changes */
/*     vdep,vdem    SD of energy changes from variance of subaverages */
/*     fdep,fdem    fluctuations of perturbation energies */
/*     se_ap,se_am  standard error of free energy changes */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Numbers of First and Last File to Analyz"
	    "e :  \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_40[] = "(/,\002 Enter the Lambda Increment for FEP : "
	    " \002,$)";
    static char fmt_50[] = "(f7.4)";
    static char fmt_60[] = "(/,\002 Enter the System Temperature [300 K] : "
	    " \002,$)";
    static char fmt_70[] = "(f20.0)";
    static char fmt_80[] = "(/,\002 Enter Number of Blocks for Sub-Averages "
	    "[1] :  \002,$)";
    static char fmt_90[] = "(i10)";
    static char fmt_100[] = "(/,\002 Consider only Intermolecular Perturbati"
	    "on\002,\002 Energy [N] :  \002,$)";
    static char fmt_110[] = "(a120)";
    static char fmt_120[] = "(/,\002 Calculation will Involve Full Perturbat"
	    "ion\002,\002 Energy \002)";
    static char fmt_130[] = "(/,\002 Calculation will Consider Only Intermol"
	    "ecular\002,\002 Interactions \002)";
    static char fmt_140[] = "(/,2x,\002Step\002,6x,\002EB\002,8x,\002EA\002,"
	    "7x,\002EIT\002,8x,\002ET\002,8x,\002EV\002,8x,\002EC\002,/)";
    static char fmt_150[] = "(i6,6f10.4)";
    static char fmt_170[] = "(/,2x,\002Step\002,8x,\002E0\002,10x,\002EP\002"
	    ",10x,\002EM\002,10x,\002DEP\002,9x,\002DEM\002,/)";
    static char fmt_180[] = "(i6,5f12.4)";
    static char fmt_190[] = "(/,2x,\002Block\002,6x,\002NStep\002,7x,\002BAD"
	    "EP\002,7x,\002BADEM\002,8x,\002BDAP\002,8x,\002BDAM\002,/)";
    static char fmt_200[] = "(i6,5x,i6,1x,4f12.4)";
    static char fmt_210[] = "(/,\002 Running Averages over\002,i5,\002 Step"
	    "s\002,\002 with Std Error from\002,i4,\002 Blocks :\002)";
    static char fmt_220[] = "(/,\002 Free Energy :\002)";
    static char fmt_230[] = "(/,\002 DA(+) =\002,f12.4,\002 with Std Erro"
	    "r\002,f10.4)";
    static char fmt_240[] = "(\002 DA(-) =\002,f12.4,\002 with Std Error\002"
	    ",f10.4)";
    static char fmt_250[] = "(/,\002 Potential Energy :\002)";
    static char fmt_260[] = "(/,\002 DE(+) =\002,f12.4,\002 with Fluct\002,f"
	    "10.4,\002 and Std Error\002,f10.4)";
    static char fmt_270[] = "(\002 DE(-) =\002,f12.4,\002 with Fluct\002,f10"
	    ".4,\002 and Std Error\002,f10.4)";
    static char fmt_280[] = "(/,\002 Component Energies :\002,/)";
    static char fmt_290[] = "(\002 BOND  +/- :\002,f12.4,5x,f12.4)";
    static char fmt_300[] = "(\002 ANGLE +/- :\002,f12.4,5x,f12.4)";
    static char fmt_310[] = "(\002 IMPT  +/- :\002,f12.4,5x,f12.4)";
    static char fmt_320[] = "(\002 TORS  +/- :\002,f12.4,5x,f12.4)";
    static char fmt_330[] = "(\002 VDW   +/- :\002,f12.4,5x,f12.4)";
    static char fmt_340[] = "(\002 CHG   +/- :\002,f12.4,5x,f12.4)";

    /* System generated locals */
    address a__1[3];
    integer i__1, i__2[3], i__3;
    doublereal d__1, d__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);
    double exp(doublereal), log(doublereal), sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    static integer modblock;
    extern integer freeunit_(void);
    static integer i__, j, k;
    static doublereal v, e0, ca[15000]	/* was [3][5000] */, da, cb[15000]	
	    /* was [3][5000] */, cc[15000]	/* was [3][5000] */, ct[15000]
	    	/* was [3][5000] */, cv[15000]	/* was [3][5000] */, rt, eb0, 
	    ea0, ec0, et0, ev0, bda, abm, dap, neg, dam, sbm, sap, sbp, scp, 
	    sam, scm, abp, aap, atp, avp, acp, pos, aam, stm, atm, stp, svm, 
	    svp, avm, acm, nrg[15000]	/* was [3][5000] */, cit[15000]	/* 
	    was [3][5000] */;
    static char ext[7];
    static doublereal lam0, eit0, bdam[100], adem, bdem, adep, bdep, bneg, 
	    fdep, fdem, bdap[100], lamm, lamp, sdep, sdem, sneg, vdap, vdam, 
	    bpos, vdep, temp, vdem, aitp, aitm;
    static integer lext;
    static doublereal sitm;
    static integer next, nrun;
    static doublereal spos;
    static integer stop;
    static doublereal sitp, a2dem, adem2, a2dep, adep2;
    static integer ixyz;
    static doublereal s2dem, sdem2, s2dep, sdep2, adamb, adapb, badem[100], 
	    badep[100], eaneg, ebneg, ecneg;
    static integer ilamb;
    static doublereal se_am__;
    static integer nlamb;
    static doublereal delta, sdapb, sdamb, se_ep__, se_em__, se_ap__;
    extern /* Subroutine */ int final_(void);
    static doublereal etneg, evneg, ebpos, eapos, ecpos;
    extern /* Subroutine */ int hatom_(void);
    static integer istep;
    static doublereal eplus;
    static integer nstep;
    static doublereal etpos, evpos;
    static integer start, nblock;
    static logical dogeom;
    static doublereal eitneg;
    static char record[120];
    extern doublereal energy_(void);
    static integer modrun;
    static doublereal eminus, deplus, eitpos;
    static char answer[1];
    extern /* Subroutine */ int upcase_(char *, ftnlen), hybrid_(void), 
	    getxyz_(void), initial_(void), numeral_(integer *, char *, 
	    integer *, ftnlen);
    static doublereal deminus;
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen), readxyz_(
	    integer *);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_20, 0 };
    static icilist io___6 = { 1, record, 1, 0, 120, 1 };
    static cilist io___9 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___105 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___106 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___114 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___115 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___144 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___145 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___146 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___147 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___148 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___149 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___150 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___151 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___152 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___153 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___154 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___155 = { 0, 0, 0, fmt_320, 0 };
    static cilist io___156 = { 0, 0, 0, fmt_330, 0 };
    static cilist io___157 = { 0, 0, 0, fmt_340, 0 };



#define ca_ref(a_1,a_2) ca[(a_2)*3 + a_1 - 4]
#define cb_ref(a_1,a_2) cb[(a_2)*3 + a_1 - 4]
#define cc_ref(a_1,a_2) cc[(a_2)*3 + a_1 - 4]
#define ct_ref(a_1,a_2) ct[(a_2)*3 + a_1 - 4]
#define cv_ref(a_1,a_2) cv[(a_2)*3 + a_1 - 4]
#define nrg_ref(a_1,a_2) nrg[(a_2)*3 + a_1 - 4]
#define cit_ref(a_1,a_2) cit[(a_2)*3 + a_1 - 4]



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
/*     ##  analyz.i  --  energy components partitioned over atoms  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     aesum   total potential energy partitioned over atoms */
/*     aeb     bond stretch energy partitioned over atoms */
/*     aea     angle bend energy partitioned over atoms */
/*     aeba    stretch-bend energy partitioned over atoms */
/*     aeub    Urey-Bradley energy partitioned over atoms */
/*     aeaa    angle-angle energy partitioned over atoms */
/*     aeopb   out-of-plane bend energy partitioned over atoms */
/*     aeopd   out-of-plane distance energy partitioned over atoms */
/*     aeid    improper dihedral energy partitioned over atoms */
/*     aeit    improper torsion energy partitioned over atoms */
/*     aet     torsional energy partitioned over atoms */
/*     aept    pi-orbital torsion energy partitioned over atoms */
/*     aebt    stretch-torsion energy partitioned over atoms */
/*     aett    torsion-torsion energy partitioned over atoms */
/*     aev     van der Waals energy partitioned over atoms */
/*     aec     charge-charge energy partitioned over atoms */
/*     aecd    charge-dipole energy partitioned over atoms */
/*     aed     dipole-dipole energy partitioned over atoms */
/*     aem     multipole energy partitioned over atoms */
/*     aep     polarization energy partitioned over atoms */
/*     aer     reaction field energy partitioned over atoms */
/*     aes     solvation energy partitioned over atoms */
/*     aelf    metal ligand field energy partitioned over atoms */
/*     aeg     geometric restraint energy partitioned over atoms */
/*     aex     extra energy term partitioned over atoms */




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
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




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




/*     set up the structure and mechanics calculation */

    initial_();
    getxyz_();
    mechanic_();

/*     get the numbers of the files to be used */

    start = 0;
    stop = 0;
    io___3.ciunit = iounit_1.iout;
    s_wsfe(&io___3);
    e_wsfe();
    io___4.ciunit = iounit_1.input;
    s_rsfe(&io___4);
    do_fio(&c__1, record, (ftnlen)120);
    e_rsfe();
    i__1 = s_rsli(&io___6);
    if (i__1 != 0) {
	goto L30;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&start, (ftnlen)sizeof(integer));
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
L30:
    if (start == 0) {
	start = 1;
    }
    if (stop == 0) {
	stop = start;
    }
    nrun = stop - start + 1;

/*     obtain the lambda values to be calculated */

    delta = 0.;
    io___9.ciunit = iounit_1.iout;
    s_wsfe(&io___9);
    e_wsfe();
    io___10.ciunit = iounit_1.input;
    s_rsfe(&io___10);
    do_fio(&c__1, (char *)&delta, (ftnlen)sizeof(doublereal));
    e_rsfe();
    nlamb = 3;
    lam0 = mutant_1.lambda;
/* Computing MIN */
    d__1 = 1., d__2 = mutant_1.lambda + delta;
    lamp = min(d__1,d__2);
/* Computing MAX */
    d__1 = 0., d__2 = mutant_1.lambda - delta;
    lamm = max(d__1,d__2);

/*     obtain the target temperature value */

    temp = 0.;
    io___16.ciunit = iounit_1.iout;
    s_wsfe(&io___16);
    e_wsfe();
    io___17.ciunit = iounit_1.input;
    s_rsfe(&io___17);
    do_fio(&c__1, (char *)&temp, (ftnlen)sizeof(doublereal));
    e_rsfe();
    if (temp == 0.) {
	temp = 300.;
    }
    rt = temp * .0019872066;

/*     set number of steps for running averages and block averages */

    nblock = 0;
    io___20.ciunit = iounit_1.iout;
    s_wsfe(&io___20);
    e_wsfe();
    io___21.ciunit = iounit_1.input;
    s_rsfe(&io___21);
    do_fio(&c__1, (char *)&nblock, (ftnlen)sizeof(integer));
    e_rsfe();
    if (nblock == 0) {
	nblock = 1;
    }
    nblock = nrun / nblock;

/*     decide whether to include the intramolecular energies */

    dogeom = TRUE_;
    io___23.ciunit = iounit_1.iout;
    s_wsfe(&io___23);
    e_wsfe();
    io___24.ciunit = iounit_1.input;
    s_rsfe(&io___24);
    do_fio(&c__1, record, (ftnlen)120);
    e_rsfe();
    next = 1;
    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'Y') {
	dogeom = FALSE_;
    }
    if (dogeom) {
	io___27.ciunit = iounit_1.iout;
	s_wsfe(&io___27);
	e_wsfe();
    } else {
	io___28.ciunit = iounit_1.iout;
	s_wsfe(&io___28);
	e_wsfe();
    }

/*     cycle over the coordinate files once per lambda value */

    i__1 = nlamb;
    for (ilamb = 1; ilamb <= i__1; ++ilamb) {
	i__ = start;
	istep = 0;
	if (ilamb == 2) {
	    mutant_1.lambda = lamp;
	}
	if (ilamb == 3) {
	    mutant_1.lambda = lamm;
	}
	hybrid_();

/*     read in the next molecular dynamics coordinate frame */

	while(i__ >= start && i__ <= stop) {
	    ++istep;
	    lext = 3;
	    numeral_(&i__, ext, &lext, (ftnlen)7);
/* Writing concatenation */
	    i__2[0] = files_1.leng, a__1[0] = files_1.filename;
	    i__2[1] = 1, a__1[1] = ".";
	    i__2[2] = lext, a__1[2] = ext;
	    s_cat(xyzfile, a__1, i__2, &c__3, (ftnlen)120);
	    version_(xyzfile, "old", (ftnlen)120, (ftnlen)3);
	    ixyz = freeunit_();
	    o__1.oerr = 1;
	    o__1.ounit = ixyz;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = xyzfile;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    i__3 = f_open(&o__1);
	    if (i__3 != 0) {
		goto L160;
	    }
	    readxyz_(&ixyz);
	    cl__1.cerr = 0;
	    cl__1.cunit = ixyz;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	    hatom_();

/*     select interactions for perturbation energy calculation */

	    i__3 = atoms_1.n;
	    for (j = 1; j <= i__3; ++j) {
		usage_1.use[j - 1] = FALSE_;
	    }
	    i__3 = mutant_1.nmut;
	    for (j = 1; j <= i__3; ++j) {
		usage_1.use[mutant_1.imut[j - 1] - 1] = TRUE_;
	    }
	    if (! dogeom) {
		potent_1.use_bond__ = FALSE_;
		potent_1.use_angle__ = FALSE_;
		potent_1.use_strbnd__ = FALSE_;
		potent_1.use_urey__ = FALSE_;
		potent_1.use_imptor__ = FALSE_;
		potent_1.use_tors__ = FALSE_;
		potent_1.use_strtor__ = FALSE_;
	    }

/*     compute and store energy components for the current lambda */

	    nrg_ref(ilamb, istep) = energy_();
	    cb_ref(ilamb, istep) = energi_1.eb;
	    ca_ref(ilamb, istep) = energi_1.ea;
	    cit_ref(ilamb, istep) = energi_1.eit;
	    ct_ref(ilamb, istep) = energi_1.et;
	    cv_ref(ilamb, istep) = energi_1.ev;
	    cc_ref(ilamb, istep) = energi_1.ec;
	    if (inform_1.verbose) {
		if (istep == 1) {
		    io___44.ciunit = iounit_1.iout;
		    s_wsfe(&io___44);
		    e_wsfe();
		}
		io___45.ciunit = iounit_1.iout;
		s_wsfe(&io___45);
		do_fio(&c__1, (char *)&istep, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&energi_1.eb, (ftnlen)sizeof(doublereal)
			);
		do_fio(&c__1, (char *)&energi_1.ea, (ftnlen)sizeof(doublereal)
			);
		do_fio(&c__1, (char *)&energi_1.eit, (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&energi_1.et, (ftnlen)sizeof(doublereal)
			);
		do_fio(&c__1, (char *)&energi_1.ev, (ftnlen)sizeof(doublereal)
			);
		do_fio(&c__1, (char *)&energi_1.ec, (ftnlen)sizeof(doublereal)
			);
		e_wsfe();
	    }
L160:
	    ++i__;
	}
    }
    nstep = istep;

/*     get free energy change by averaging over all frames */

    i__1 = nstep;
    for (istep = 1; istep <= i__1; ++istep) {
	e0 = nrg_ref(1, istep);
	eplus = nrg_ref(2, istep);
	eminus = nrg_ref(3, istep);
	ev0 = cv_ref(1, istep);
	evpos = cv_ref(2, istep);
	evneg = cv_ref(3, istep);
	ec0 = cc_ref(1, istep);
	ecpos = cc_ref(2, istep);
	ecneg = cc_ref(3, istep);
	if (dogeom) {
	    eb0 = cb_ref(1, istep);
	    ebpos = cb_ref(2, istep);
	    ebneg = cb_ref(3, istep);
	    ea0 = ca_ref(1, istep);
	    eapos = ca_ref(2, istep);
	    eaneg = ca_ref(3, istep);
	    eit0 = cit_ref(1, istep);
	    eitpos = cit_ref(2, istep);
	    eitneg = cit_ref(3, istep);
	    et0 = ct_ref(1, istep);
	    etpos = ct_ref(2, istep);
	    etneg = ct_ref(3, istep);
	}
	modrun = (istep - 1) % nrun;
	modblock = (istep - 1) % nblock;

/*     zero out summation variables for new running average */

	if (modrun == 0) {
	    sdep = 0.;
	    s2dep = 0.;
	    sdep2 = 0.;
	    sdem = 0.;
	    s2dem = 0.;
	    sdem2 = 0.;
	    spos = 0.;
	    sneg = 0.;
	    sdapb = 0.;
	    sdamb = 0.;
	    sbp = 0.;
	    sbm = 0.;
	    sap = 0.;
	    sam = 0.;
	    sitp = 0.;
	    sitm = 0.;
	    stp = 0.;
	    stm = 0.;
	    svp = 0.;
	    svm = 0.;
	    scp = 0.;
	    scm = 0.;
	    fdep = 0.;
	    fdem = 0.;
	    vdep = 0.;
	    vdem = 0.;
	    vdap = 0.;
	    vdam = 0.;
	    k = 0;
	}

/*     zero out summation variables for new block average */

	if (modblock == 0) {
	    bdep = 0.;
	    bdem = 0.;
	    bpos = 0.;
	    bneg = 0.;
	}
	modrun = istep % nrun;
	modblock = istep % nblock;

/*     accumulate statistics */

	deplus = eplus - e0;
	deminus = eminus - e0;
	if (inform_1.verbose) {
	    if (modblock == 1 || nblock == 1) {
		io___105.ciunit = iounit_1.iout;
		s_wsfe(&io___105);
		e_wsfe();
	    }
	    io___106.ciunit = iounit_1.iout;
	    s_wsfe(&io___106);
	    do_fio(&c__1, (char *)&istep, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&e0, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&eplus, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&eminus, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&deplus, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&deminus, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	pos = exp(-deplus / rt);
	neg = exp(-deminus / rt);
	bdep += deplus;
	bdem += deminus;
	bpos += pos;
	bneg += neg;
	sdep += deplus;
	sdem += deminus;
	s2dep += deplus * deplus;
	s2dem += deminus * deminus;
	spos += pos;
	sneg += neg;
	svp = svp + evpos - ev0;
	svm = svm + evneg - ev0;
	scp = scp + ecpos - ec0;
	scm = scm + ecneg - ec0;
	if (dogeom) {
	    sbp = sbp + ebpos - eb0;
	    sbm = sbm + ebneg - eb0;
	    sap = sap + eapos - ea0;
	    sam = sam + eaneg - ea0;
	    sitp = sitp + eitpos - eit0;
	    sitm = sitm + eitneg - eit0;
	    stp = stp + etpos - et0;
	    stm = stm + etneg - et0;
	}

/*     calculate block averages */

	if (modblock == 0) {
	    ++k;
	    badep[k - 1] = bdep / (doublereal) nblock;
	    badem[k - 1] = bdem / (doublereal) nblock;
	    bda = bpos / (doublereal) nblock;
	    bda = -rt * log(bda);
	    bdap[k - 1] += bda;
	    bda = bneg / (doublereal) nblock;
	    bda = -rt * log(bda);
	    bdam[k - 1] += bda;
	    sdapb += bdap[k - 1];
	    sdamb += bdam[k - 1];
	    if (inform_1.verbose || k == 1) {
		io___114.ciunit = iounit_1.iout;
		s_wsfe(&io___114);
		e_wsfe();
	    }
	    io___115.ciunit = iounit_1.iout;
	    s_wsfe(&io___115);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&istep, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&badep[k - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&badem[k - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&bdap[k - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&bdam[k - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}

/*     calculate running averages for potential energy */

	if (modrun == 0) {
	    adep = sdep / (doublereal) nrun;
	    adem = sdem / (doublereal) nrun;
	    a2dep = s2dep / (doublereal) nrun;
	    a2dem = s2dem / (doublereal) nrun;
	    adep2 = adep * adep;
	    adem2 = adem * adem;
	    fdep = sqrt(a2dep - adep2);
	    fdem = sqrt(a2dem - adem2);
	    i__3 = nrun / nblock;
	    for (k = 1; k <= i__3; ++k) {
/* Computing 2nd power */
		d__1 = badep[k - 1] - adep;
		v = d__1 * d__1;
		vdep += v;
/* Computing 2nd power */
		d__1 = badem[k - 1] - adem;
		v = d__1 * d__1;
		vdem += v;
	    }
	    vdep /= (doublereal) (nrun / nblock);
	    se_ep__ = sqrt(vdep / (doublereal) (nrun / nblock));
	    vdem /= (doublereal) (nrun / nblock);
	    se_em__ = sqrt(vdem / (doublereal) (nrun / nblock));

/*     calculate running averages for free energy */

	    da = spos / (doublereal) nrun;
	    da = -rt * log(da);
	    dap = da;
	    da = sneg / (doublereal) nrun;
	    da = -rt * log(da);
	    dam = da;
	    adapb = sdapb / (doublereal) (nrun / nblock);
	    adamb = sdamb / (doublereal) (nrun / nblock);
	    i__3 = nrun / nblock;
	    for (k = 1; k <= i__3; ++k) {
/* Computing 2nd power */
		d__1 = bdap[k - 1] - adapb;
		v = d__1 * d__1;
		vdap += v;
/* Computing 2nd power */
		d__1 = bdam[k - 1] - adamb;
		v = d__1 * d__1;
		vdam += v;
	    }
	    vdap /= (doublereal) (nrun / nblock);
	    se_ap__ = sqrt(vdap / (doublereal) (nrun / nblock));
	    vdam /= (doublereal) (nrun / nblock);
	    se_am__ = sqrt(vdam / (doublereal) (nrun / nblock));

/*     calculate running averages for energy components */

	    avp = svp / (doublereal) nrun;
	    avm = svm / (doublereal) nrun;
	    acp = scp / (doublereal) nrun;
	    acm = scm / (doublereal) nrun;
	    if (dogeom) {
		abp = sbp / (doublereal) nrun;
		abm = sbm / (doublereal) nrun;
		aap = sap / (doublereal) nrun;
		aam = sam / (doublereal) nrun;
		aitp = sitp / (doublereal) nrun;
		aitm = sitm / (doublereal) nrun;
		atp = stp / (doublereal) nrun;
		atm = stm / (doublereal) nrun;
	    }
	    sdep = 0.;
	    sdem = 0.;

/*     write information about running averages and block averages */

	    io___144.ciunit = iounit_1.iout;
	    s_wsfe(&io___144);
	    do_fio(&c__1, (char *)&nstep, (ftnlen)sizeof(integer));
	    i__3 = nrun / nblock;
	    do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___145.ciunit = iounit_1.iout;
	    s_wsfe(&io___145);
	    e_wsfe();
	    io___146.ciunit = iounit_1.iout;
	    s_wsfe(&io___146);
	    do_fio(&c__1, (char *)&dap, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&se_ap__, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    io___147.ciunit = iounit_1.iout;
	    s_wsfe(&io___147);
	    do_fio(&c__1, (char *)&dam, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&se_am__, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    io___148.ciunit = iounit_1.iout;
	    s_wsfe(&io___148);
	    e_wsfe();
	    io___149.ciunit = iounit_1.iout;
	    s_wsfe(&io___149);
	    do_fio(&c__1, (char *)&adep, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&fdep, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&se_ep__, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    io___150.ciunit = iounit_1.iout;
	    s_wsfe(&io___150);
	    do_fio(&c__1, (char *)&adem, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&fdem, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&se_em__, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    io___151.ciunit = iounit_1.iout;
	    s_wsfe(&io___151);
	    e_wsfe();
	    if (dogeom) {
		io___152.ciunit = iounit_1.iout;
		s_wsfe(&io___152);
		do_fio(&c__1, (char *)&abp, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&abm, (ftnlen)sizeof(doublereal));
		e_wsfe();
		io___153.ciunit = iounit_1.iout;
		s_wsfe(&io___153);
		do_fio(&c__1, (char *)&aap, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&aam, (ftnlen)sizeof(doublereal));
		e_wsfe();
		io___154.ciunit = iounit_1.iout;
		s_wsfe(&io___154);
		do_fio(&c__1, (char *)&aitp, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&aitm, (ftnlen)sizeof(doublereal));
		e_wsfe();
		io___155.ciunit = iounit_1.iout;
		s_wsfe(&io___155);
		do_fio(&c__1, (char *)&atp, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&atm, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    io___156.ciunit = iounit_1.iout;
	    s_wsfe(&io___156);
	    do_fio(&c__1, (char *)&avp, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&avm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    io___157.ciunit = iounit_1.iout;
	    s_wsfe(&io___157);
	    do_fio(&c__1, (char *)&acp, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&acm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef cit_ref
#undef nrg_ref
#undef cv_ref
#undef ct_ref
#undef cc_ref
#undef cb_ref
#undef ca_ref


/* Main program alias */ int alchemy_ () { MAIN__ (); return 0; }
