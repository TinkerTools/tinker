/* egeom3.f -- translated by f2c (version 20050501).
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
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b119 = 0.;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine egeom3  --  restraint energy terms & analysis  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "egeom3" calculates the energy due to restraints on positions, */
/*     distances, angles and torsions as well as Gaussian basin and */
/*     droplet restraints; also partitions energy among the atoms */


/* Subroutine */ int egeom3_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Individual Atomic Position Restrain"
	    "t\002,\002 Terms :\002,//,\002 Type\002,8x,\002Atom Name\002,11x,"
	    "\002Target\002,\002 Position\002,7x,\002Distance\002,6x,\002Ener"
	    "gy\002,/)";
    static char fmt_20[] = "(\002 Position\002,3x,i5,\002-\002,a3,2x,f10.4,5"
	    "x,\002----\002,6x,\002----\002,1x,f10.4,f12.4)";
    static char fmt_30[] = "(\002 Position\002,3x,i5,\002-\002,a3,7x,\002-"
	    "---\002,1x,f10.4,5x,\002----\002,1x,f10.4,f12.4)";
    static char fmt_40[] = "(\002 Position\002,3x,i5,\002-\002,a3,7x,\002-"
	    "---\002,6x,\002----\002,1x,2f10.4,f12.4)";
    static char fmt_50[] = "(\002 Position\002,3x,i5,\002-\002,a3,2x,4f10.4,"
	    "f12.4)";
    static char fmt_60[] = "(/,\002 Individual Interatomic Distance\002,\002"
	    " Restraint Terms :\002,//,\002 Type\002,13x,\002Atom Names\002,1"
	    "4x,\002Ideal Range\002,4x,\002Actual\002,6x,\002Energy\002,/)";
    static char fmt_70[] = "(\002 Distance\002,3x,i5,\002-\002,a3,1x,i5"
	    ",\002-\002,a3,6x,2f8.2,f10.4,f12.4)";
    static char fmt_80[] = "(/,\002 Individual Interatomic Angle\002,\002 Re"
	    "straint Terms :\002,//,\002 Type\002,12x,\002Atom Numbers\002,13"
	    "x,\002Ideal\002,\002 Range\002,4x,\002Actual\002,6x,\002Energ"
	    "y\002,/)";
    static char fmt_90[] = "(\002 Angle\002,8x,3i5,8x,2f8.2,f10.4,f12.4)";
    static char fmt_100[] = "(/,\002 Individual Torsional Angle Restraint  "
	    " \002,\002 Terms :\002,//,\002 Type\002,12x,\002Atom Numbers\002"
	    ",13x,\002Ideal\002,\002 Range\002,4x,\002Actual\002,6x,\002Energy"
	    "\002,/)";
    static char fmt_110[] = "(\002 Torsion\002,4x,4i5,5x,2f8.2,f10.4,f12.4)";
    static char fmt_120[] = "(/,\002 Individual Intergroup Distance\002,\002"
	    " Restraint Terms :\002,//,\002 Type\002,12x,\002Group Numbers"
	    "\002,12x,\002Ideal Range\002,4x,\002Actual\002,6x,\002Energy\002"
	    ",/)";
    static char fmt_130[] = "(\002 Distance\002,8x,i5,1x,i5,9x,2f8.2,f10.4,f"
	    "12.4)";
    static char fmt_140[] = "(/,\002 Individual Chirality Restraint Terms "
	    ":\002,//,\002 Type\002,12x,\002Atom Numbers\002,13x,\002Ideal"
	    "\002,\002 Range\002,4x,\002Actual\002,6x,\002Energy\002,/)";
    static char fmt_150[] = "(\002 Chiral\002,5x,4i5,5x,2f8.2,f10.4,f12.4)";
    static char fmt_160[] = "(/,\002 Individual Gaussian Basin\002,\002 Rest"
	    "raint Terms :\002,//,\002 Type\002,13x,\002Atom Names\002,20x"
	    ",\002Ideal\002,4x,\002Actual\002,6x,\002Energy\002,/)";
    static char fmt_170[] = "(\002 Distance\002,3x,i5,\002-\002,a3,1x,i5,"
	    "\002-\002,a3,12x,2f10.4,f12.4)";
    static char fmt_180[] = "(/,\002 Individual Spherical Boundary\002,\002 "
	    "Restraint Terms :\002,//,\002 Type\002,13x,\002Atom Name\002,28x,"
	    "\002Distance\002,6x,\002Energy\002,/)";
    static char fmt_190[] = "(\002 Wall\002,12x,i5,\002-\002,a3,27x,f10.4,f1"
	    "2.4)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    double acos(doublereal), exp(doublereal);

    /* Local variables */
    static logical intermol;
    static doublereal a, b, e;
    static integer i__, j, k;
    static doublereal r__, c1, c2, c3, r2, t1, t2, r6;
    static integer ia, ib, ic, id;
    static doublereal r12, dt, ri, xi, yi, zi, xr, yr, zr, xt, yt, zt, xu, yu,
	     zu, af1, af2, cf1, df1, df2, cf2, gf1, gf2, dt2, tf1, tf2, rt2, 
	    ru2, dot, rcb, xia, yia, zia, xib, yib, zib, xic, yic, zic, xid, 
	    yid, zid, xab, yab, zab, xba, xtu, ytu, ztu, yba, zba, xcb, rab2, 
	    ycb, rcb2, zcb, xdc, ydc, zdc, xad, yad, zad, xbd, ybd, zbd, xcd, 
	    ycd, zcd, xcm, ycm, zcm, vol, fgrp, sine, size, term;
    static logical huge__;
    static doublereal rtru;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal angle, force, weigh;
    static logical header;
    static doublereal weigha, weighb, buffer, cosine, target;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;

    /* Fortran I/O blocks */
    static cilist io___15 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___80 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___81 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___93 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___94 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___110 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___111 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___117 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___118 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___125 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___126 = { 0, 0, 0, fmt_190, 0 };



#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define chir_ref(a_1,a_2) kgeoms_1.chir[(a_2)*3 + a_1 - 4]
#define dfix_ref(a_1,a_2) kgeoms_1.dfix[(a_2)*3 + a_1 - 4]
#define afix_ref(a_1,a_2) kgeoms_1.afix[(a_2)*3 + a_1 - 4]
#define gfix_ref(a_1,a_2) kgeoms_1.gfix[(a_2)*3 + a_1 - 4]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]
#define pfix_ref(a_1,a_2) kgeoms_1.pfix[(a_2)*2 + a_1 - 3]
#define tfix_ref(a_1,a_2) kgeoms_1.tfix[(a_2)*3 + a_1 - 4]
#define ichir_ref(a_1,a_2) kgeoms_1.ichir[(a_2)*4 + a_1 - 5]
#define iafix_ref(a_1,a_2) kgeoms_1.iafix[(a_2)*3 + a_1 - 4]
#define idfix_ref(a_1,a_2) kgeoms_1.idfix[(a_2)*2 + a_1 - 3]
#define igfix_ref(a_1,a_2) kgeoms_1.igfix[(a_2)*2 + a_1 - 3]
#define kpfix_ref(a_1,a_2) kgeoms_1.kpfix[(a_2)*3 + a_1 - 4]
#define itfix_ref(a_1,a_2) kgeoms_1.itfix[(a_2)*4 + a_1 - 5]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     zero out the restraint energy and partitioning terms */

    action_1.neg = 0;
    energi_1.eg = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	analyz_1.aeg[i__ - 1] = 0.;
    }

/*     compute the energy for position restraint terms */

    header = TRUE_;
    i__1 = kgeoms_1.npfix;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = kgeoms_1.ipfix[i__ - 1];
	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &c__0, &c__0, &c__0, &c__0, &c__0);
	}
	if (proceed) {
	    proceed = usage_1.use[ia - 1];
	}
	if (proceed) {
	    xr = 0.;
	    yr = 0.;
	    zr = 0.;
	    if (kpfix_ref(1, i__) != 0) {
		xr = atoms_1.x[ia - 1] - kgeoms_1.xpfix[i__ - 1];
	    }
	    if (kpfix_ref(2, i__) != 0) {
		yr = atoms_1.y[ia - 1] - kgeoms_1.ypfix[i__ - 1];
	    }
	    if (kpfix_ref(3, i__) != 0) {
		zr = atoms_1.z__[ia - 1] - kgeoms_1.zpfix[i__ - 1];
	    }
	    r__ = sqrt(xr * xr + yr * yr + zr * zr);
	    force = pfix_ref(1, i__);
/* Computing MAX */
	    d__1 = 0., d__2 = r__ - pfix_ref(2, i__);
	    dt = max(d__1,d__2);
	    dt2 = dt * dt;
	    e = force * dt2;
	    if (group_1.use_group__) {
		e *= fgrp;
	    }
	    ++action_1.neg;
	    energi_1.eg += e;
	    analyz_1.aeg[ia - 1] += e;
	    huge__ = e > 10.;
	    if (inform_1.debug || inform_1.verbose && huge__) {
		if (header) {
		    header = FALSE_;
		    io___15.ciunit = iounit_1.iout;
		    s_wsfe(&io___15);
		    e_wsfe();
		}
		dt = sqrt(dt2);
		if (kpfix_ref(2, i__) == 0 && kpfix_ref(3, i__) == 0) {
		    io___16.ciunit = iounit_1.iout;
		    s_wsfe(&io___16);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, name___ref(0, ia), (ftnlen)3);
		    do_fio(&c__1, (char *)&kgeoms_1.xpfix[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&dt, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else if (kpfix_ref(1, i__) == 0 && kpfix_ref(3, i__) == 0) {
		    io___17.ciunit = iounit_1.iout;
		    s_wsfe(&io___17);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, name___ref(0, ia), (ftnlen)3);
		    do_fio(&c__1, (char *)&kgeoms_1.ypfix[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&dt, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else if (kpfix_ref(1, i__) == 0 && kpfix_ref(2, i__) == 0) {
		    io___18.ciunit = iounit_1.iout;
		    s_wsfe(&io___18);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, name___ref(0, ia), (ftnlen)3);
		    do_fio(&c__1, (char *)&kgeoms_1.zpfix[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&dt, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else {
		    io___19.ciunit = iounit_1.iout;
		    s_wsfe(&io___19);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, name___ref(0, ia), (ftnlen)3);
		    do_fio(&c__1, (char *)&kgeoms_1.xpfix[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&kgeoms_1.ypfix[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&kgeoms_1.zpfix[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&dt, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	}
    }

/*     compute the energy for distance restraint terms */

    header = TRUE_;
    i__1 = kgeoms_1.ndfix;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = idfix_ref(1, i__);
	ib = idfix_ref(2, i__);
	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &c__0, &c__0, &c__0, &c__0);
	}
	if (proceed) {
	    proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1];
	}
	if (proceed) {
	    xr = atoms_1.x[ia - 1] - atoms_1.x[ib - 1];
	    yr = atoms_1.y[ia - 1] - atoms_1.y[ib - 1];
	    zr = atoms_1.z__[ia - 1] - atoms_1.z__[ib - 1];
	    intermol = molcul_1.molcule[ia - 1] != molcul_1.molcule[ib - 1];
	    if (bound_1.use_bounds__ && intermol) {
		image_(&xr, &yr, &zr);
	    }
	    r__ = sqrt(xr * xr + yr * yr + zr * zr);
	    force = dfix_ref(1, i__);
	    df1 = dfix_ref(2, i__);
	    df2 = dfix_ref(3, i__);
	    target = r__;
	    if (r__ < df1) {
		target = df1;
	    }
	    if (r__ > df2) {
		target = df2;
	    }
	    dt = r__ - target;
	    dt2 = dt * dt;
	    e = force * dt2;
	    if (group_1.use_group__) {
		e *= fgrp;
	    }
	    ++action_1.neg;
	    energi_1.eg += e;
	    analyz_1.aeg[ia - 1] += e * .5;
	    analyz_1.aeg[ib - 1] += e * .5;
	    if (intermol) {
		inter_1.einter += e;
	    }
	    huge__ = e > 10.;
	    if (inform_1.debug || inform_1.verbose && huge__) {
		if (header) {
		    header = FALSE_;
		    io___25.ciunit = iounit_1.iout;
		    s_wsfe(&io___25);
		    e_wsfe();
		}
		io___26.ciunit = iounit_1.iout;
		s_wsfe(&io___26);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, ia), (ftnlen)3);
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, ib), (ftnlen)3);
		do_fio(&c__1, (char *)&df1, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&df2, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&r__, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
    }

/*     compute the energy for angle restraint terms */

    header = TRUE_;
    i__1 = kgeoms_1.nafix;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iafix_ref(1, i__);
	ib = iafix_ref(2, i__);
	ic = iafix_ref(3, i__);
	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &c__0, &c__0, &c__0);
	}
	if (proceed) {
	    proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
		    usage_1.use[ic - 1];
	}
	if (proceed) {
	    xia = atoms_1.x[ia - 1];
	    yia = atoms_1.y[ia - 1];
	    zia = atoms_1.z__[ia - 1];
	    xib = atoms_1.x[ib - 1];
	    yib = atoms_1.y[ib - 1];
	    zib = atoms_1.z__[ib - 1];
	    xic = atoms_1.x[ic - 1];
	    yic = atoms_1.y[ic - 1];
	    zic = atoms_1.z__[ic - 1];
	    xab = xia - xib;
	    yab = yia - yib;
	    zab = zia - zib;
	    xcb = xic - xib;
	    ycb = yic - yib;
	    zcb = zic - zib;
	    rab2 = xab * xab + yab * yab + zab * zab;
	    rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
	    if (rab2 != 0. && rcb2 != 0.) {
		dot = xab * xcb + yab * ycb + zab * zcb;
		cosine = dot / sqrt(rab2 * rcb2);
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine);
		cosine = min(d__1,d__2);
		angle = acos(cosine) * 57.29577951308232088;
		force = afix_ref(1, i__);
		af1 = afix_ref(2, i__);
		af2 = afix_ref(3, i__);
		target = angle;
		if (angle < af1) {
		    target = af1;
		}
		if (angle > af2) {
		    target = af2;
		}
		dt = angle - target;
		dt2 = dt * dt;
		e = force * dt2;
		if (group_1.use_group__) {
		    e *= fgrp;
		}
		++action_1.neg;
		energi_1.eg += e;
		analyz_1.aeg[ib - 1] += e * .5;
		if (intermol) {
		    inter_1.einter += e;
		}
		huge__ = e > 10.;
		if (inform_1.debug || inform_1.verbose && huge__) {
		    if (header) {
			header = FALSE_;
			io___50.ciunit = iounit_1.iout;
			s_wsfe(&io___50);
			e_wsfe();
		    }
		    io___51.ciunit = iounit_1.iout;
		    s_wsfe(&io___51);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&af1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&af2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&angle, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	}
    }

/*     compute the energy for torsional restraint terms */

    header = TRUE_;
    i__1 = kgeoms_1.ntfix;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = itfix_ref(1, i__);
	ib = itfix_ref(2, i__);
	ic = itfix_ref(3, i__);
	id = itfix_ref(4, i__);
	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &c__0, &c__0);
	}
	if (proceed) {
	    proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
		    usage_1.use[ic - 1] || usage_1.use[id - 1];
	}
	if (proceed) {
	    xia = atoms_1.x[ia - 1];
	    yia = atoms_1.y[ia - 1];
	    zia = atoms_1.z__[ia - 1];
	    xib = atoms_1.x[ib - 1];
	    yib = atoms_1.y[ib - 1];
	    zib = atoms_1.z__[ib - 1];
	    xic = atoms_1.x[ic - 1];
	    yic = atoms_1.y[ic - 1];
	    zic = atoms_1.z__[ic - 1];
	    xid = atoms_1.x[id - 1];
	    yid = atoms_1.y[id - 1];
	    zid = atoms_1.z__[id - 1];
	    xba = xib - xia;
	    yba = yib - yia;
	    zba = zib - zia;
	    xcb = xic - xib;
	    ycb = yic - yib;
	    zcb = zic - zib;
	    xdc = xid - xic;
	    ydc = yid - yic;
	    zdc = zid - zic;
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
		force = tfix_ref(1, i__);
		tf1 = tfix_ref(2, i__);
		tf2 = tfix_ref(3, i__);
		if (angle > tf1 && angle < tf2) {
		    target = angle;
		} else if (angle > tf1 && tf1 > tf2) {
		    target = angle;
		} else if (angle < tf2 && tf1 > tf2) {
		    target = angle;
		} else {
		    t1 = angle - tf1;
		    t2 = angle - tf2;
		    if (t1 > 180.) {
			t1 += -360.;
		    } else if (t1 < -180.) {
			t1 += 360.;
		    }
		    if (t2 > 180.) {
			t2 += -360.;
		    } else if (t2 < -180.) {
			t2 += 360.;
		    }
		    if (abs(t1) < abs(t2)) {
			target = tf1;
		    } else {
			target = tf2;
		    }
		}
		dt = angle - target;
		if (dt > 180.) {
		    dt += -360.;
		} else if (dt < -180.) {
		    dt += 360.;
		}
		dt2 = dt * dt;
		e = force * dt2;
		if (group_1.use_group__) {
		    e *= fgrp;
		}
		++action_1.neg;
		energi_1.eg += e;
		analyz_1.aeg[ib - 1] += e * .5;
		analyz_1.aeg[ic - 1] += e * .5;
		if (molcul_1.molcule[ia - 1] != molcul_1.molcule[ib - 1] || 
			molcul_1.molcule[ia - 1] != molcul_1.molcule[ic - 1] 
			|| molcul_1.molcule[ia - 1] != molcul_1.molcule[id - 
			1]) {
		    inter_1.einter += e;
		}
		huge__ = e > 10.;
		if (inform_1.debug || inform_1.verbose && huge__) {
		    if (header) {
			header = FALSE_;
			io___80.ciunit = iounit_1.iout;
			s_wsfe(&io___80);
			e_wsfe();
		    }
		    io___81.ciunit = iounit_1.iout;
		    s_wsfe(&io___81);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&tf1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&tf2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&angle, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	}
    }

/*     compute the energy for group distance restraint terms */

    header = TRUE_;
    i__1 = kgeoms_1.ngfix;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = igfix_ref(1, i__);
	ib = igfix_ref(2, i__);
	xcm = 0.;
	ycm = 0.;
	zcm = 0.;
	i__2 = igrp_ref(2, ia);
	for (j = igrp_ref(1, ia); j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    weigh = atmtyp_1.mass[k - 1];
	    xcm += atoms_1.x[k - 1] * weigh;
	    ycm += atoms_1.y[k - 1] * weigh;
	    zcm += atoms_1.z__[k - 1] * weigh;
	}
/* Computing MAX */
	d__1 = 1., d__2 = group_1.grpmass[ia - 1];
	weigha = max(d__1,d__2);
	xr = xcm / weigha;
	yr = ycm / weigha;
	zr = zcm / weigha;
	xcm = 0.;
	ycm = 0.;
	zcm = 0.;
	i__2 = igrp_ref(2, ib);
	for (j = igrp_ref(1, ib); j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    weigh = atmtyp_1.mass[k - 1];
	    xcm += atoms_1.x[k - 1] * weigh;
	    ycm += atoms_1.y[k - 1] * weigh;
	    zcm += atoms_1.z__[k - 1] * weigh;
	}
/* Computing MAX */
	d__1 = 1., d__2 = group_1.grpmass[ib - 1];
	weighb = max(d__1,d__2);
	xr -= xcm / weighb;
	yr -= ycm / weighb;
	zr -= zcm / weighb;
	intermol = molcul_1.molcule[group_1.kgrp[igrp_ref(1, ia) - 1] - 1] != 
		molcul_1.molcule[group_1.kgrp[igrp_ref(1, ib) - 1] - 1];
	if (bound_1.use_bounds__ && intermol) {
	    image_(&xr, &yr, &zr);
	}
	r__ = sqrt(xr * xr + yr * yr + zr * zr);
	force = gfix_ref(1, i__);
	gf1 = gfix_ref(2, i__);
	gf2 = gfix_ref(3, i__);
	target = r__;
	if (r__ < gf1) {
	    target = gf1;
	}
	if (r__ > gf2) {
	    target = gf2;
	}
	dt = r__ - target;
	dt2 = dt * dt;
	e = force * dt2;
	++action_1.neg;
	energi_1.eg += e;
	size = (doublereal) (igrp_ref(2, ia) - igrp_ref(1, ia) + 1);
	i__2 = igrp_ref(2, ia);
	for (j = igrp_ref(1, ia); j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    analyz_1.aeg[k - 1] += e * .5 / size;
	}
	size = (doublereal) (igrp_ref(2, ib) - igrp_ref(1, ib) + 1);
	i__2 = igrp_ref(2, ib);
	for (j = igrp_ref(1, ib); j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    analyz_1.aeg[k - 1] += e * .5 / size;
	}
	if (intermol) {
	    inter_1.einter += e;
	}
	huge__ = e > 10.;
	if (inform_1.debug || inform_1.verbose && huge__) {
	    if (header) {
		header = FALSE_;
		io___93.ciunit = iounit_1.iout;
		s_wsfe(&io___93);
		e_wsfe();
	    }
	    io___94.ciunit = iounit_1.iout;
	    s_wsfe(&io___94);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&gf1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gf2, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&r__, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     compute the energy for chirality restraint terms */

    header = TRUE_;
    i__1 = kgeoms_1.nchir;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ichir_ref(1, i__);
	ib = ichir_ref(2, i__);
	ic = ichir_ref(3, i__);
	id = ichir_ref(4, i__);
	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &c__0, &c__0);
	}
	if (proceed) {
	    proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
		    usage_1.use[ic - 1] || usage_1.use[id - 1];
	}
	if (proceed) {
	    xad = atoms_1.x[ia - 1] - atoms_1.x[id - 1];
	    yad = atoms_1.y[ia - 1] - atoms_1.y[id - 1];
	    zad = atoms_1.z__[ia - 1] - atoms_1.z__[id - 1];
	    xbd = atoms_1.x[ib - 1] - atoms_1.x[id - 1];
	    ybd = atoms_1.y[ib - 1] - atoms_1.y[id - 1];
	    zbd = atoms_1.z__[ib - 1] - atoms_1.z__[id - 1];
	    xcd = atoms_1.x[ic - 1] - atoms_1.x[id - 1];
	    ycd = atoms_1.y[ic - 1] - atoms_1.y[id - 1];
	    zcd = atoms_1.z__[ic - 1] - atoms_1.z__[id - 1];
	    c1 = ybd * zcd - zbd * ycd;
	    c2 = ycd * zad - zcd * yad;
	    c3 = yad * zbd - zad * ybd;
	    vol = xad * c1 + xbd * c2 + xcd * c3;
	    force = chir_ref(1, i__);
	    cf1 = chir_ref(2, i__);
	    cf2 = chir_ref(3, i__);
	    target = vol;
	    if (vol < min(cf1,cf2)) {
		target = min(cf1,cf2);
	    }
	    if (vol > max(cf1,cf2)) {
		target = max(cf1,cf2);
	    }
	    dt = vol - target;
	    dt2 = dt * dt;
	    e = force * dt2;
	    if (group_1.use_group__) {
		e *= fgrp;
	    }
	    ++action_1.neg;
	    energi_1.eg += e;
	    analyz_1.aeg[ia - 1] += e * .25;
	    analyz_1.aeg[ib - 1] += e * .25;
	    analyz_1.aeg[ic - 1] += e * .25;
	    analyz_1.aeg[id - 1] += e * .25;
	    huge__ = e > 10.;
	    if (inform_1.debug || inform_1.verbose && huge__) {
		if (header) {
		    header = FALSE_;
		    io___110.ciunit = iounit_1.iout;
		    s_wsfe(&io___110);
		    e_wsfe();
		}
		io___111.ciunit = iounit_1.iout;
		s_wsfe(&io___111);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&cf1, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&cf2, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&vol, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
    }

/*     compute the energy for a Gaussian basin restraint */

    if (kgeoms_1.use_basin__) {
	header = TRUE_;
	i__1 = atoms_1.n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    i__2 = atoms_1.n;
	    for (k = i__ + 1; k <= i__2; ++k) {
		proceed = TRUE_;
		if (group_1.use_group__) {
		    groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &
			    c__0);
		}
		if (proceed) {
		    proceed = usage_1.use[i__ - 1] || usage_1.use[k - 1];
		}
		if (proceed) {
		    xr = xi - atoms_1.x[k - 1];
		    yr = yi - atoms_1.y[k - 1];
		    zr = zi - atoms_1.z__[k - 1];
		    r2 = xr * xr + yr * yr + zr * zr;
		    term = -kgeoms_1.width * r2;
		    e = 0.;
		    if (term > -50.) {
			e = kgeoms_1.depth * exp(term);
		    }
		    e -= kgeoms_1.depth;
		    if (group_1.use_group__) {
			e *= fgrp;
		    }
		    ++action_1.neg;
		    energi_1.eg += e;
		    analyz_1.aeg[i__ - 1] += e * .5;
		    analyz_1.aeg[k - 1] += e * .5;
		    huge__ = e > 10.;
		    if (inform_1.debug || inform_1.verbose && huge__) {
			if (header) {
			    header = FALSE_;
			    io___117.ciunit = iounit_1.iout;
			    s_wsfe(&io___117);
			    e_wsfe();
			}
			r__ = sqrt(r2);
			io___118.ciunit = iounit_1.iout;
			s_wsfe(&io___118);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			do_fio(&c__1, name___ref(0, k), (ftnlen)3);
			do_fio(&c__1, (char *)&c_b119, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&r__, (ftnlen)sizeof(doublereal)
				);
			do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
			e_wsfe();
		    }
		}
	    }
	}
    }

/*     compute the energy for a spherical droplet restraint */

    if (kgeoms_1.use_wall__) {
	header = TRUE_;
	buffer = 2.5;
	a = 2048.;
	b = 64.;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &c__0, &c__0, &c__0, &c__0, &
			c__0);
	    }
	    if (proceed) {
		proceed = usage_1.use[i__ - 1];
	    }
	    if (proceed) {
		xi = atoms_1.x[i__ - 1];
		yi = atoms_1.y[i__ - 1];
		zi = atoms_1.z__[i__ - 1];
/* Computing 2nd power */
		d__1 = xi;
/* Computing 2nd power */
		d__2 = yi;
/* Computing 2nd power */
		d__3 = zi;
		ri = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
		r__ = kgeoms_1.rwall + buffer - ri;
		r2 = r__ * r__;
		r6 = r2 * r2 * r2;
		r12 = r6 * r6;
		e = a / r12 - b / r6;
		if (group_1.use_group__) {
		    e *= fgrp;
		}
		++action_1.neg;
		energi_1.eg += e;
		analyz_1.aeg[i__ - 1] += e;
		huge__ = e > 10.;
		if (inform_1.debug || inform_1.verbose && huge__) {
		    if (header) {
			header = FALSE_;
			io___125.ciunit = iounit_1.iout;
			s_wsfe(&io___125);
			e_wsfe();
		    }
		    io___126.ciunit = iounit_1.iout;
		    s_wsfe(&io___126);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
		    do_fio(&c__1, (char *)&ri, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	}
    }
    return 0;
} /* egeom3_ */

#undef itfix_ref
#undef kpfix_ref
#undef igfix_ref
#undef idfix_ref
#undef iafix_ref
#undef ichir_ref
#undef tfix_ref
#undef pfix_ref
#undef igrp_ref
#undef gfix_ref
#undef afix_ref
#undef dfix_ref
#undef chir_ref
#undef name___ref


