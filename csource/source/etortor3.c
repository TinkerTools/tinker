/* etortor3.f -- translated by f2c (version 20050501).
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
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

struct {
    integer nbitor, ibitor[500000]	/* was [5][100000] */;
} bitor_;

#define bitor_1 bitor_

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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal ttx[3000]	/* was [30][100] */, tty[3000]	/* was [30][
	    100] */, tbf[90000]	/* was [900][100] */, tbx[90000]	/* 
	    was [900][100] */, tby[90000]	/* was [900][100] */, tbxy[
	    90000]	/* was [900][100] */;
    integer tnx[100], tny[100];
    char ktt[2000];
} ktrtor_;

#define ktrtor_1 ktrtor_

struct {
    doublereal idihunit, itorunit, torsunit, ptorunit, storunit, ttorunit;
} torpot_;

#define torpot_1 torpot_

struct {
    integer ntortor, itt[300000]	/* was [3][100000] */;
} tortor_;

#define tortor_1 tortor_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;



/*     ############################################################# */
/*     ##  COPYRIGHT (C) 2003 by Pengyu Ren & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine etortor3  --  torsion-torsion energy & analysis  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "etortor3" calculates the torsion-torsion potential energy; */
/*     also partitions the energy terms among the atoms */


/* Subroutine */ int etortor3_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Individual Torsion-Torsion\002,\002 Inte"
	    "ractions :\002,//,\002 Type\002,17x,\002Atom Numbers\002,13x,"
	    "\002Angle1\002,4x,\002Angle2\002,6x,\002Energy\002,/)";
    static char fmt_20[] = "(\002 TorTor\002,4x,5(1x,i5),2x,2f10.4,f12.4)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal e;
    static integer i__, k, ia, ib, ic, id, ie;
    static doublereal xh, yh;
    static integer nt;
    static doublereal xt, yt, zt, xu, yu, zu, xv, yv, zv, ft1[4], ft2[4], x1l,
	     y1l, rt2, ru2, rv2, x1u, y1u;
    static integer nlo, nhi, xlo, ylo;
    static doublereal xia, yia, zia, xib, yib, zib, xic, xtu, ytu, ztu, xuv, 
	    yuv, zuv, yic, zic, xid, yid, zid, xie, yie, zie, xba, yba, zba, 
	    xdc, ydc, zdc, xcb, ycb, zcb, xed, yed, zed, ftt[4], ft12[4];
    static integer pos1, pos2;
    static logical huge__;
    static doublereal fgrp, sign, rtru, rurv;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal angle1, angle2, value1, value2;
    static logical header;
    extern /* Subroutine */ int bcuint_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal cosine1, cosine2;
    static logical proceed;
    extern /* Subroutine */ int chkttor_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static integer itortor;

    /* Fortran I/O blocks */
    static cilist io___85 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___86 = { 0, 0, 0, fmt_20, 0 };



#define tbf_ref(a_1,a_2) ktrtor_1.tbf[(a_2)*900 + a_1 - 901]
#define tbx_ref(a_1,a_2) ktrtor_1.tbx[(a_2)*900 + a_1 - 901]
#define tby_ref(a_1,a_2) ktrtor_1.tby[(a_2)*900 + a_1 - 901]
#define itt_ref(a_1,a_2) tortor_1.itt[(a_2)*3 + a_1 - 4]
#define ttx_ref(a_1,a_2) ktrtor_1.ttx[(a_2)*30 + a_1 - 31]
#define tty_ref(a_1,a_2) ktrtor_1.tty[(a_2)*30 + a_1 - 31]
#define tbxy_ref(a_1,a_2) ktrtor_1.tbxy[(a_2)*900 + a_1 - 901]
#define ibitor_ref(a_1,a_2) bitor_1.ibitor[(a_2)*5 + a_1 - 6]



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
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  bitor.i  --  bitorsions within the current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     nbitor  total number of bitorsions in the system */
/*     ibitor  numbers of the atoms in each bitorsion */




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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  ktrtor.i  --  forcefield parameters for torsion-torsions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxntt    maximum number of torsion-torsion parameter entries */
/*     maxtgrd   maximum dimension of torsion-torsion spline grid */
/*     maxtgrd2  maximum number of torsion-torsion spline grid points */

/*     ttx       angle values for first torsion of spline grid */
/*     tty       angle values for second torsion of spline grid */
/*     tbf       function values at points on spline grid */
/*     tbx       gradient over first torsion of spline grid */
/*     tby       gradient over second torsion of spline grid */
/*     tbxy      Hessian cross components over spline grid */
/*     tnx       number of columns in torsion-torsion spline grid */
/*     tny       number of rows in torsion-torsion spline grid */
/*     ktt       string of torsion-torsion atom classes */




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
/*     ##  torpot.i  --  specifics of torsional functional forms  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     idihunit  convert improper dihedral energy to kcal/mole */
/*     itorunit  convert improper torsion amplitudes to kcal/mole */
/*     torsunit  convert torsional parameter amplitudes to kcal/mole */
/*     ptorunit  convert pi-orbital torsion energy to kcal/mole */
/*     storunit  convert stretch-torsion energy to kcal/mole */
/*     ttorunit  convert stretch-torsion energy to kcal/mole */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  tortor.i  --  torsion-torsions in the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     ntortor   total number of torsion-torsion interactions */
/*     itt       atoms and parameter indices for torsion-torsion */




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




/*     zero out the torsion-torsion energy and partitioning terms */

    action_1.nett = 0;
    energi_1.ett = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	analyz_1.aett[i__ - 1] = 0.;
    }
    header = TRUE_;

/*     calculate the torsion-torsion interaction energy term */

    i__1 = tortor_1.ntortor;
    for (itortor = 1; itortor <= i__1; ++itortor) {
	i__ = itt_ref(1, itortor);
	k = itt_ref(2, itortor);
	if (itt_ref(3, itortor) == 1) {
	    ia = ibitor_ref(1, i__);
	    ib = ibitor_ref(2, i__);
	    ic = ibitor_ref(3, i__);
	    id = ibitor_ref(4, i__);
	    ie = ibitor_ref(5, i__);
	} else {
	    ia = ibitor_ref(5, i__);
	    ib = ibitor_ref(4, i__);
	    ic = ibitor_ref(3, i__);
	    id = ibitor_ref(2, i__);
	    ie = ibitor_ref(1, i__);
	}

/*     decide whether to compute the current interaction */

	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &ie, &c__0);
	}
	if (proceed) {
	    proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
		    usage_1.use[ic - 1] || usage_1.use[id - 1] || usage_1.use[
		    ie - 1];
	}

/*     compute the value of the torsional angles */

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
	    xie = atoms_1.x[ie - 1];
	    yie = atoms_1.y[ie - 1];
	    zie = atoms_1.z__[ie - 1];
	    xba = xib - xia;
	    yba = yib - yia;
	    zba = zib - zia;
	    xcb = xic - xib;
	    ycb = yic - yib;
	    zcb = zic - zib;
	    xdc = xid - xic;
	    ydc = yid - yic;
	    zdc = zid - zic;
	    xed = xie - xid;
	    yed = yie - yid;
	    zed = zie - zid;
	    if (bound_1.use_polymer__) {
		image_(&xba, &yba, &zba);
		image_(&xcb, &ycb, &zcb);
		image_(&xdc, &ydc, &zdc);
		image_(&xed, &yed, &zed);
	    }
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
	    xv = ydc * zed - yed * zdc;
	    yv = zdc * xed - zed * xdc;
	    zv = xdc * yed - xed * ydc;
	    xuv = yu * zv - yv * zu;
	    yuv = zu * xv - zv * xu;
	    zuv = xu * yv - xv * yu;
	    rv2 = xv * xv + yv * yv + zv * zv;
	    rurv = sqrt(ru2 * rv2);
	    if (rtru != 0. && rurv != 0.) {
		cosine1 = (xt * xu + yt * yu + zt * zu) / rtru;
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine1);
		cosine1 = min(d__1,d__2);
		angle1 = acos(cosine1) * 57.29577951308232088;
		sign = xba * xu + yba * yu + zba * zu;
		if (sign < 0.) {
		    angle1 = -angle1;
		}
		value1 = angle1;
		cosine2 = (xu * xv + yu * yv + zu * zv) / rurv;
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine2);
		cosine2 = min(d__1,d__2);
		angle2 = acos(cosine2) * 57.29577951308232088;
		sign = xcb * xv + ycb * yv + zcb * zv;
		if (sign < 0.) {
		    angle2 = -angle2;
		}
		value2 = angle2;

/*     check for inverted chirality at the central atom */

		chkttor_(&ib, &ic, &id, &sign, &value1, &value2);

/*     use bicubic interpolation to compute spline values */

		nlo = 1;
		nhi = ktrtor_1.tnx[k - 1];
		while(nhi - nlo > 1) {
		    nt = (nhi + nlo) / 2;
		    if (ttx_ref(nt, k) > value1) {
			nhi = nt;
		    } else {
			nlo = nt;
		    }
		}
		xlo = nlo;
		nlo = 1;
		nhi = ktrtor_1.tny[k - 1];
		while(nhi - nlo > 1) {
		    nt = (nhi + nlo) / 2;
		    if (tty_ref(nt, k) > value2) {
			nhi = nt;
		    } else {
			nlo = nt;
		    }
		}
		ylo = nlo;
		xh = ttx_ref(xlo + 1, k) - ttx_ref(xlo, k);
		yh = tty_ref(ylo + 1, k) - tty_ref(ylo, k);
		x1l = ttx_ref(xlo, k);
		x1u = ttx_ref(xlo + 1, k);
		y1l = tty_ref(ylo, k);
		y1u = tty_ref(ylo + 1, k);
		pos2 = ylo * ktrtor_1.tnx[k - 1] + xlo;
		pos1 = pos2 - ktrtor_1.tnx[k - 1];
		ftt[0] = tbf_ref(pos1, k);
		ftt[1] = tbf_ref(pos1 + 1, k);
		ftt[2] = tbf_ref(pos2 + 1, k);
		ftt[3] = tbf_ref(pos2, k);
		ft1[0] = tbx_ref(pos1, k);
		ft1[1] = tbx_ref(pos1 + 1, k);
		ft1[2] = tbx_ref(pos2 + 1, k);
		ft1[3] = tbx_ref(pos2, k);
		ft2[0] = tby_ref(pos1, k);
		ft2[1] = tby_ref(pos1 + 1, k);
		ft2[2] = tby_ref(pos2 + 1, k);
		ft2[3] = tby_ref(pos2, k);
		ft12[0] = tbxy_ref(pos1, k);
		ft12[1] = tbxy_ref(pos1 + 1, k);
		ft12[2] = tbxy_ref(pos2 + 1, k);
		ft12[3] = tbxy_ref(pos2, k);
		bcuint_(ftt, ft1, ft2, ft12, &x1l, &x1u, &y1l, &y1u, &value1, 
			&value2, &e);
		e = torpot_1.ttorunit * e;

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    e *= fgrp;
		}

/*     increment the total torsion-torsion energy */

		++action_1.nett;
		energi_1.ett += e;
		analyz_1.aett[ib - 1] += e / 3.;
		analyz_1.aett[ic - 1] += e / 3.;
		analyz_1.aett[id - 1] += e / 3.;

/*     print a message if the energy of this interaction is large */

		huge__ = e > 3.;
		if (inform_1.debug || inform_1.verbose && huge__) {
		    if (header) {
			header = FALSE_;
			io___85.ciunit = iounit_1.iout;
			s_wsfe(&io___85);
			e_wsfe();
		    }
		    io___86.ciunit = iounit_1.iout;
		    s_wsfe(&io___86);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ie, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&angle1, (ftnlen)sizeof(doublereal))
			    ;
		    do_fio(&c__1, (char *)&angle2, (ftnlen)sizeof(doublereal))
			    ;
		    do_fio(&c__1, (char *)&e, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	}
    }
    return 0;
} /* etortor3_ */

#undef ibitor_ref
#undef tbxy_ref
#undef tty_ref
#undef ttx_ref
#undef itt_ref
#undef tby_ref
#undef tbx_ref
#undef tbf_ref


