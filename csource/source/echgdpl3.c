/* echgdpl3.f -- translated by f2c (version 20050501).
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
    integer neb, nea, neba, neub, neaa, neopb, neopd, neid, neit, net, nept, 
	    nebt, nett, nev, nec, necd, ned, nem, nep, new__, ner, nes, nelf, 
	    neg, nex;
} action_;

#define action_1 action_

struct {
    doublereal aesum[25000], aeb[25000], aea[25000], aeba[25000], aeub[25000],
	     aeaa[25000], aeopb[25000], aeopd[25000], aeid[25000], aeit[25000]
	    , aet[25000], aept[25000], aebt[25000], aett[25000], aev[25000], 
	    aec[25000], aecd[25000], aed[25000], aem[25000], aep[25000], aer[
	    25000], aes[25000], aelf[25000], aeg[25000], aex[25000];
} analyz_;

#define analyz_1 analyz_

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
    doublereal xcell, ycell, zcell, xcell2, ycell2, zcell2;
    integer ncell, icell[30000]	/* was [3][10000] */;
} cell_;

#define cell_1 cell_

struct {
    doublereal pchg[25000];
    integer nion, iion[25000], jion[25000], kion[25000], chglist[25000];
} charge_;

#define charge_1 charge_

struct {
    doublereal electric, dielec, ebuffer, c2scale, c3scale, c4scale, c5scale;
    logical neutnbr, neutcut;
} chgpot_;

#define chgpot_1 chgpot_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal bdpl[50000], sdpl[50000];
    integer ndipole, idpl[100000]	/* was [2][50000] */;
} dipole_;

#define dipole_1 dipole_

struct {
    doublereal esum, eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, 
	    ebt, ett, ev, ec, ecd, ed, em, ep, er, es, elf, eg, ex;
} energi_;

#define energi_1 energi_

struct {
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    doublereal einter;
} inter_;

#define inter_1 inter_

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
    doublereal off, off2, cut, cut2, c0, c1, c2, c3, c4, c5, f0, f1, f2, f3, 
	    f4, f5, f6, f7;
} shunt_;

#define shunt_1 shunt_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine echgdpl3  --  charge-dipole energy & analysis  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "echgdpl3" calculates the charge-dipole interaction energy; */
/*     also partitions the energy among the atoms */


/* Subroutine */ int echgdpl3_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Individual Charge-Dipole\002,\002 Intera"
	    "ctions :\002,//,\002 Type\002,10x,\002Charge\002,11x,\002Dipol"
	    "e\002,17x,\002Distance\002,6x,\002Energy\002,/)";
    static char fmt_20[] = "(\002 Chg-Dpl\002,4x,i5,\002-\002,a3,\002 / \002"
	    ",i5,\002-\002,a3,1x,i5,\002-\002,a3,10x,f10.4,f12.4)";
    static char fmt_30[] = "(/,\002 Individual Charge-Dipole\002,\002 Intera"
	    "ctions :\002,//,\002 Type\002,10x,\002Charge\002,11x,\002Dipol"
	    "e\002,17x,\002Distance\002,6x,\002Energy\002,/)";
    static char fmt_40[] = "(\002 Chg-Dpl\002,4x,i5,\002-\002,a3,\002 / \002"
	    ",i5,\002-\002,a3,1x,i5,\002-\002,a3,2x,\002(XTAL)\002,2x,f10.4,f"
	    "12.4)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__;
    static integer i1, k1, k2;
    static doublereal r2, r3, r4, r5, fi, xi, yi, zi, xk, yk, zk, xr, yr, zr, 
	    rk2, fik, rkr3;
    static integer skip[25000];
    static doublereal dotk, fgrp;
    static logical huge__;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal taper;
    static logical header;
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *), switch_(char *, ftnlen), groups_(
	    logical *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *);
    static logical proceed;

    /* Fortran I/O blocks */
    static cilist io___33 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_40, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
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
/*     ##  action.i  --  total number of each energy term computed  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     neb     number of bond stretch energy terms computed */
/*     nea     number of angle bend energy terms computed */
/*     neba    number of stretch-bend energy terms computed */
/*     neub    number of Urey-Bradley energy terms computed */
/*     neaa    number of angle-angle energy terms computed */
/*     neopb   number of out-of-plane bend energy terms computed */
/*     neopd   number of out-of-plane distance energy terms computed */
/*     neid    number of improper dihedral energy terms computed */
/*     neit    number of improper torsion energy terms computed */
/*     net     number of torsional energy terms computed */
/*     nept    number of pi-orbital torsion energy terms computed */
/*     nebt    number of stretch-torsion energy terms computed */
/*     nett    number of torsion-torsion energy terms computed */
/*     nev     number of van der Waals energy terms computed */
/*     nec     number of charge-charge energy terms computed */
/*     necd    number of charge-dipole energy terms computed */
/*     ned     number of dipole-dipole energy terms computed */
/*     nem     number of multipole energy terms computed */
/*     nep     number of polarization energy terms computed */
/*     new     number of Ewald summation energy terms computed */
/*     ner     number of reaction field energy terms computed */
/*     nes     number of solvation energy terms computed */
/*     nelf    number of metal ligand field energy terms computed */
/*     neg     number of geometric restraint energy terms computed */
/*     nex     number of extra energy terms computed */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cell.i  --  periodic boundaries using replicated cells  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xcell    length of the a-axis of the complete replicated cell */
/*     ycell    length of the b-axis of the complete replicated cell */
/*     zcell    length of the c-axis of the complete replicated cell */
/*     xcell2   half the length of the a-axis of the replicated cell */
/*     ycell2   half the length of the b-axis of the replicated cell */
/*     zcell2   half the length of the c-axis of the replicated cell */
/*     ncell    total number of cell replicates for periodic boundaries */
/*     icell    offset along axes for each replicate periodic cell */




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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  inter.i  --  sum of intermolecular energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     einter   total intermolecular potential energy */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  shunt.i  --  polynomial switching function coefficients  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     off    distance at which the potential energy goes to zero */
/*     off2   square of distance at which the potential goes to zero */
/*     cut    distance at which switching of the potential begins */
/*     cut2   square of distance at which the switching begins */
/*     c0     zeroth order coefficient of multiplicative switch */
/*     c1     first order coefficient of multiplicative switch */
/*     c2     second order coefficient of multiplicative switch */
/*     c3     third order coefficient of multiplicative switch */
/*     c4     fourth order coefficient of multiplicative switch */
/*     c5     fifth order coefficient of multiplicative switch */
/*     f0     zeroth order coefficient of additive switch function */
/*     f1     first order coefficient of additive switch function */
/*     f2     second order coefficient of additive switch function */
/*     f3     third order coefficient of additive switch function */
/*     f4     fourth order coefficient of additive switch function */
/*     f5     fifth order coefficient of additive switch function */
/*     f6     sixth order coefficient of additive switch function */
/*     f7     seventh order coefficient of additive switch function */




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




/*     zero out the overall charge-dipole interaction energy */
/*     and partitioning; set up constants for the calculation */

    action_1.necd = 0;
    energi_1.ecd = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	analyz_1.aecd[i__ - 1] = 0.;
    }
    if (dipole_1.ndipole == 0 || charge_1.nion == 0) {
	return 0;
    }
    header = TRUE_;

/*     zero out the list of atoms to be skipped */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	skip[i__ - 1] = 0;
    }

/*     set conversion factor and switching function coefficients */

    f = chgpot_1.electric / (chgpot_1.dielec * 4.80321);
    switch_("CHGDPL", (ftnlen)6);

/*     get the total energy by looping over each charge-dipole pair */

    i__1 = charge_1.nion;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = charge_1.iion[i__ - 1];
	skip[i1 - 1] = i1;
	i__2 = couple_1.n12[i1 - 1];
	for (k = 1; k <= i__2; ++k) {
	    skip[i12_ref(k, i1) - 1] = i1;
	}
	xi = atoms_1.x[i1 - 1];
	yi = atoms_1.y[i1 - 1];
	zi = atoms_1.z__[i1 - 1];
	fi = f * charge_1.pchg[i__ - 1];

/*     decide whether to compute the current interaction */

	i__2 = dipole_1.ndipole;
	for (k = 1; k <= i__2; ++k) {
	    k1 = idpl_ref(1, k);
	    k2 = idpl_ref(2, k);
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i1, &k1, &k2, &c__0, &c__0, &c__0);
	    }
	    if (proceed) {
		proceed = usage_1.use[i1 - 1] || usage_1.use[k1 - 1] || 
			usage_1.use[k2 - 1];
	    }
	    if (proceed) {
		proceed = skip[k1 - 1] != i1 && skip[k2 - 1] != i1;
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xk = atoms_1.x[k2 - 1] - atoms_1.x[k1 - 1];
		yk = atoms_1.y[k2 - 1] - atoms_1.y[k1 - 1];
		zk = atoms_1.z__[k2 - 1] - atoms_1.z__[k1 - 1];
		xr = atoms_1.x[k1 - 1] + xk * dipole_1.sdpl[k - 1] - xi;
		yr = atoms_1.y[k1 - 1] + yk * dipole_1.sdpl[k - 1] - yi;
		zr = atoms_1.z__[k1 - 1] + zk * dipole_1.sdpl[k - 1] - zi;
		image_(&xr, &yr, &zr);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    fik = fi * dipole_1.bdpl[k - 1];
		    rk2 = xk * xk + yk * yk + zk * zk;
		    rkr3 = sqrt(rk2 * r2) * r2;
		    dotk = xk * xr + yk * yr + zk * zr;
		    e = fik * dotk / rkr3;

/*     use energy switching if near the cutoff distance */

		    if (r2 > shunt_1.cut2) {
			r__ = sqrt(r2);
			r3 = r2 * r__;
			r4 = r2 * r2;
			r5 = r2 * r3;
			taper = shunt_1.c5 * r5 + shunt_1.c4 * r4 + 
				shunt_1.c3 * r3 + shunt_1.c2 * r2 + 
				shunt_1.c1 * r__ + shunt_1.c0;
			e *= taper;
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
		    }

/*     increment the overall charge-dipole energy component */

		    ++action_1.necd;
		    energi_1.ecd += e;
		    analyz_1.aecd[i1 - 1] += e * .5;
		    analyz_1.aecd[k1 - 1] += e * .25;
		    analyz_1.aecd[k2 - 1] += e * .25;

/*     increment the total intermolecular energy */

		    if (molcul_1.molcule[i1 - 1] != molcul_1.molcule[k1 - 1]) 
			    {
			inter_1.einter += e;
		    }

/*     print a message if the energy of this interaction is large */

		    huge__ = abs(e) > 25.;
		    if (inform_1.debug || inform_1.verbose && huge__) {
			if (header) {
			    header = FALSE_;
			    io___33.ciunit = iounit_1.iout;
			    s_wsfe(&io___33);
			    e_wsfe();
			}
			io___34.ciunit = iounit_1.iout;
			s_wsfe(&io___34);
			do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(integer));
			do_fio(&c__1, name___ref(0, i1), (ftnlen)3);
			do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
			do_fio(&c__1, name___ref(0, k1), (ftnlen)3);
			do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
			do_fio(&c__1, name___ref(0, k2), (ftnlen)3);
			d__1 = sqrt(r2);
			do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
			e_wsfe();
		    }
		}
	    }
	}
    }

/*     for periodic boundary conditions with large cutoffs */
/*     neighbors must be found by the replicates method */

    if (! bound_1.use_replica__) {
	return 0;
    }

/*     calculate interaction energy with other unit cells */

    i__1 = charge_1.nion;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = charge_1.iion[i__ - 1];
	skip[i1 - 1] = i1;
	i__2 = couple_1.n12[i1 - 1];
	for (k = 1; k <= i__2; ++k) {
	    skip[i12_ref(k, i1) - 1] = i1;
	}
	xi = atoms_1.x[i1 - 1];
	yi = atoms_1.y[i1 - 1];
	zi = atoms_1.z__[i1 - 1];
	fi = f * charge_1.pchg[i__ - 1];

/*     decide whether to compute the current interaction */

	i__2 = dipole_1.ndipole;
	for (k = 1; k <= i__2; ++k) {
	    k1 = idpl_ref(1, k);
	    k2 = idpl_ref(2, k);
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i1, &k1, &k2, &c__0, &c__0, &c__0);
	    }
	    if (proceed) {
		proceed = usage_1.use[i1 - 1] || usage_1.use[k1 - 1] || 
			usage_1.use[k2 - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		i__3 = cell_1.ncell;
		for (j = 1; j <= i__3; ++j) {
		    xk = atoms_1.x[k2 - 1] - atoms_1.x[k1 - 1];
		    yk = atoms_1.y[k2 - 1] - atoms_1.y[k1 - 1];
		    zk = atoms_1.z__[k2 - 1] - atoms_1.z__[k1 - 1];
		    xr = atoms_1.x[k1 - 1] + xk * dipole_1.sdpl[k - 1] - xi;
		    yr = atoms_1.y[k1 - 1] + yk * dipole_1.sdpl[k - 1] - yi;
		    zr = atoms_1.z__[k1 - 1] + zk * dipole_1.sdpl[k - 1] - zi;
		    imager_(&xr, &yr, &zr, &j);
		    r2 = xr * xr + yr * yr + zr * zr;
		    if (r2 <= shunt_1.off2) {
			fik = fi * dipole_1.bdpl[k - 1];
			if (bound_1.use_polymer__) {
			    if (r2 < bound_1.polycut2) {
				if (skip[k1 - 1] == i1 || skip[k2 - 1] != i1) 
					{
				    fik = 0.;
				}
			    }
			}
			rk2 = xk * xk + yk * yk + zk * zk;
			rkr3 = sqrt(rk2 * r2) * r2;
			dotk = xk * xr + yk * yr + zk * zr;
			e = fik * dotk / rkr3;

/*     use energy switching if near the cutoff distance */

			if (r2 > shunt_1.cut2) {
			    r__ = sqrt(r2);
			    r3 = r2 * r__;
			    r4 = r2 * r2;
			    r5 = r2 * r3;
			    taper = shunt_1.c5 * r5 + shunt_1.c4 * r4 + 
				    shunt_1.c3 * r3 + shunt_1.c2 * r2 + 
				    shunt_1.c1 * r__ + shunt_1.c0;
			    e *= taper;
			}

/*     scale the interaction based on its group membership */

			if (group_1.use_group__) {
			    e *= fgrp;
			}

/*     increment the overall charge-dipole energy component */

			if (e != 0.) {
			    ++action_1.necd;
			}
			energi_1.ecd += e;
			analyz_1.aecd[i1 - 1] += e * .5;
			analyz_1.aecd[k1 - 1] += e * .25;
			analyz_1.aecd[k2 - 1] += e * .25;

/*     increment the total intermolecular energy */

			inter_1.einter += e;

/*     print a message if the energy of this interaction is large */

			huge__ = abs(e) > 25.;
			if (inform_1.debug && e != 0. || inform_1.verbose && 
				huge__) {
			    if (header) {
				header = FALSE_;
				io___36.ciunit = iounit_1.iout;
				s_wsfe(&io___36);
				e_wsfe();
			    }
			    io___37.ciunit = iounit_1.iout;
			    s_wsfe(&io___37);
			    do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(integer)
				    );
			    do_fio(&c__1, name___ref(0, i1), (ftnlen)3);
			    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer)
				    );
			    do_fio(&c__1, name___ref(0, k1), (ftnlen)3);
			    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer)
				    );
			    do_fio(&c__1, name___ref(0, k2), (ftnlen)3);
			    d__1 = sqrt(r2);
			    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(
				    doublereal));
			    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(
				    doublereal));
			    e_wsfe();
			}
		    }
		}
	    }
	}
    }
    return 0;
} /* echgdpl3_ */

#undef idpl_ref
#undef name___ref
#undef i12_ref


