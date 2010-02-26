/* etortor2.f -- translated by f2c (version 20050501).
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
    integer nbitor, ibitor[500000]	/* was [5][100000] */;
} bitor_;

#define bitor_1 bitor_

struct {
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

struct {
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    doublereal hessx[75000]	/* was [3][25000] */, hessy[75000]	/* 
	    was [3][25000] */, hessz[75000]	/* was [3][25000] */;
} hessn_;

#define hessn_1 hessn_

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

/* Table of constant values */

static integer c__0 = 0;



/*     ############################################################# */
/*     ##  COPYRIGHT (C) 2003 by Pengyu Ren & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine etortor2  --  atomwise torsion-torsion Hessian  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "etortor2" calculates the torsion-torsion potential energy */
/*     second derivatives with respect to Cartesian coordinates */


/* Subroutine */ int etortor2_(integer *i__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal);

    /* Local variables */
    static doublereal e;
    static integer j, k, ia, ib, ic, id, ie;
    static doublereal xh, yh;
    static integer nt;
    static doublereal xt, yt, zt, xu, yu, zu, xv, yv, zv, ft1[4], ft2[4], x1l,
	     y1l, rt2, ru2, rv2, x1u, y1u, rcb, rdc, xba, yba, zba;
    static integer nhi;
    static doublereal xcb, ycb, xia, yia, zia, xib, yib, zib, xic;
    static integer nlo;
    static doublereal yic, zic, xid, yid, zid, xie;
    static integer xlo, ylo;
    static doublereal yie, zie, zcb, xdc, ydc, zdc, xca, yca, zca, xdb, ydb, 
	    zdb, xed, yed, zed, xec, yec, zec, xtu, ytu, ztu, xuv, yuv, zuv, 
	    ftt[4], ft12[4];
    static integer pos1, pos2;
    static doublereal fgrp, sign, rtru, rurv, rcbt2, rcbu2, rdcu2, rdcv2, 
	    xycb2, xzcb2, yzcb2, xydc2, xzdc2, yzdc2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal rcbxt, rcbyt, rcbzt, rcbxu, rcbyu, rcbzu, rdcxv, rdcyv, 
	    rdczv, rdcxu, rdcyu, rdczu, angle1, angle2, da1dxt, da1dyt, 
	    da1dzt, da1dxu, da1dyu, da1dzu, da2dxv, da2dyv, value1, value2, 
	    da2dzv, da2dxu, da2dyu, da2dzu, d2eda1a1, d2eda2a2, dedang1, 
	    dedang2, d2eda1a2, da1dxia, da1dyia, da1dzia, da1dxib, da1dyib, 
	    da1dzib, da1dxic, da1dyic, da1dzic, da1dxid, da1dyid, da1dzid, 
	    da2dxie, da2dyie, da2dzie, da2dxib, da2dyib, da2dzib, cosine1, 
	    cosine2, da2dxic, da2dyic, da2dzic, da2dxid, da2dyid, da2dzid;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), bcuint2_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static logical proceed;
    static doublereal dxiaxia, dxiayia, dyiayia, dziazia, dxibxib, dyibyib, 
	    dzibzib, dxicxic, dyicyic, dziczic, dxidxid, dyidyid, dzidzid, 
	    dxiazia, dyiazia, dxibyib, dxibzib, dyibzib, dxicyic, dxiczic, 
	    dyiczic, dxidyid, dxidzid, dyidzid, dxiaxib, dxiayib, dxiazib, 
	    dyiaxib, da1dxibt, da1dyibt;
    static integer itortor;
    static doublereal da1dzibt, da1dxibu, da1dyibu, da1dzibu, da1dxict, 
	    da1dyict, da1dzict, da1dxicu, da1dyicu, da1dzicu, da2dxidu, 
	    da2dyidu, da2dzidu, da2dxidv, da2dyidv, da2dzidv, da2dxicu, 
	    da2dyicu, da2dzicu, da2dxicv, da2dyicv, da2dzicv, dyiayib, 
	    dyiazib, dziaxib, dziayib, dziazib, dxiaxic, dxiayic, dxiazic, 
	    dyiaxic, dyiayic, dyiazic, dziaxic, dziayic, dziazic, dxiaxid, 
	    dxiayid, dxiazid, dyiaxid, dyiayid, dyiazid, dziaxid, dziayid, 
	    dziazid, dxibxic, dxibyic, dxibzic, dyibxic, dyibyic, dyibzic, 
	    dzibxic, dzibyic, dzibzic, dxibxid, dxibyid, dxibzid, dyibxid, 
	    dyibyid, dyibzid, dzibxid, dzibyid, dzibzid, dxicxid, dxicyid, 
	    dxiczid, dyicxid, dyicyid, dyiczid, dzicxid, dzicyid, dziczid, 
	    dxibxib2, dyibyib2, dzibzib2, dxicxic2, dyicyic2, dziczic2, 
	    dxidxid2, dyidyid2, dzidzid2, dxiexie2, dyieyie2, dziezie2, 
	    dxibyib2, dxibzib2, dyibzib2, dxicyic2, dxiczic2, dyiczic2, 
	    dxidyid2, dxidzid2, dyidzid2, dxieyie2, dxiezie2, dyiezie2, 
	    dxibxic2, dxibyic2, dxibzic2, dyibxic2, dyibyic2, dyibzic2, 
	    dzibxic2, dzibyic2, dzibzic2, dxibxid2, dxibyid2, dxibzid2, 
	    dyibxid2, dyibyid2, dyibzid2, dzibxid2, dzibyid2, dzibzid2, 
	    dxibxie2, dxibyie2, dxibzie2, dyibxie2, dyibyie2, dyibzie2, 
	    dzibxie2, dzibyie2, dzibzie2, dxicxid2, dxicyid2, dxiczid2, 
	    dyicxid2, dyicyid2, dyiczid2, dzicxid2, dzicyid2, dziczid2, 
	    dxicxie2, dxicyie2, dxiczie2, dyicxie2, dyicyie2, dyiczie2, 
	    dzicxie2, dzicyie2, dziczie2, dxidxie2, dxidyie2, dxidzie2, 
	    dyidxie2, dyidyie2, dyidzie2, dzidxie2, dzidyie2, dzidzie2;
    extern /* Subroutine */ int chkttor_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);


#define tbf_ref(a_1,a_2) ktrtor_1.tbf[(a_2)*900 + a_1 - 901]
#define tbx_ref(a_1,a_2) ktrtor_1.tbx[(a_2)*900 + a_1 - 901]
#define tby_ref(a_1,a_2) ktrtor_1.tby[(a_2)*900 + a_1 - 901]
#define itt_ref(a_1,a_2) tortor_1.itt[(a_2)*3 + a_1 - 4]
#define ttx_ref(a_1,a_2) ktrtor_1.ttx[(a_2)*30 + a_1 - 31]
#define tty_ref(a_1,a_2) ktrtor_1.tty[(a_2)*30 + a_1 - 31]
#define tbxy_ref(a_1,a_2) ktrtor_1.tbxy[(a_2)*900 + a_1 - 901]
#define hessx_ref(a_1,a_2) hessn_1.hessx[(a_2)*3 + a_1 - 4]
#define hessy_ref(a_1,a_2) hessn_1.hessy[(a_2)*3 + a_1 - 4]
#define hessz_ref(a_1,a_2) hessn_1.hessz[(a_2)*3 + a_1 - 4]
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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  hessn.i  --  Cartesian Hessian elements for a single atom  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     hessx   Hessian elements for x-component of current atom */
/*     hessy   Hessian elements for y-component of current atom */
/*     hessz   Hessian elements for z-component of current atom */




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




/*     compute the Hessian elements of the torsion-torsions */

    i__1 = tortor_1.ntortor;
    for (itortor = 1; itortor <= i__1; ++itortor) {
	j = itt_ref(1, itortor);
	k = itt_ref(2, itortor);
	if (itt_ref(3, itortor) == 1) {
	    ia = ibitor_ref(1, j);
	    ib = ibitor_ref(2, j);
	    ic = ibitor_ref(3, j);
	    id = ibitor_ref(4, j);
	    ie = ibitor_ref(5, j);
	} else {
	    ia = ibitor_ref(5, j);
	    ib = ibitor_ref(4, j);
	    ic = ibitor_ref(3, j);
	    id = ibitor_ref(2, j);
	    ie = ibitor_ref(1, j);
	}

/*     decide whether to compute the current interaction */

	proceed = *i__ == ia || *i__ == ib || *i__ == ic || *i__ == id || *
		i__ == ie;
	if (proceed && group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &ie, &c__0);
	}

/*     compute the values of the torsional angles */

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
		rcb = sqrt(xcb * xcb + ycb * ycb + zcb * zcb);
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
		rdc = sqrt(xdc * xdc + ydc * ydc + zdc * zdc);
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
		bcuint2_(ftt, ft1, ft2, ft12, &x1l, &x1u, &y1l, &y1u, &value1,
			 &value2, &e, &dedang1, &dedang2, &d2eda1a2, &
			d2eda1a1, &d2eda2a2);
		dedang1 = sign * torpot_1.ttorunit * 57.29577951308232088 * 
			dedang1;
		dedang2 = sign * torpot_1.ttorunit * 57.29577951308232088 * 
			dedang2;
		d2eda1a1 = torpot_1.ttorunit * 3282.8063500117441 * d2eda1a1;
		d2eda2a2 = torpot_1.ttorunit * 3282.8063500117441 * d2eda2a2;
		d2eda1a2 = torpot_1.ttorunit * 3282.8063500117441 * d2eda1a2;

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    dedang1 *= fgrp;
		    dedang2 *= fgrp;
		    d2eda1a1 *= fgrp;
		    d2eda2a2 *= fgrp;
		    d2eda1a2 *= fgrp;
		}

/*     abbreviations for first derivative chain rule terms */

		xca = xic - xia;
		yca = yic - yia;
		zca = zic - zia;
		xdb = xid - xib;
		ydb = yid - yib;
		zdb = zid - zib;
		if (bound_1.use_polymer__) {
		    image_(&xca, &yca, &zca);
		    image_(&xdb, &ydb, &zdb);
		}
		da1dxt = (yt * zcb - ycb * zt) / (rt2 * rcb);
		da1dyt = (zt * xcb - zcb * xt) / (rt2 * rcb);
		da1dzt = (xt * ycb - xcb * yt) / (rt2 * rcb);
		da1dxu = -(yu * zcb - ycb * zu) / (ru2 * rcb);
		da1dyu = -(zu * xcb - zcb * xu) / (ru2 * rcb);
		da1dzu = -(xu * ycb - xcb * yu) / (ru2 * rcb);

/*     abbreviations for second derivative chain rule terms */

		xycb2 = xcb * xcb + ycb * ycb;
		xzcb2 = xcb * xcb + zcb * zcb;
		yzcb2 = ycb * ycb + zcb * zcb;
		rcbxt = rcb * -2. * da1dxt;
		rcbyt = rcb * -2. * da1dyt;
		rcbzt = rcb * -2. * da1dzt;
		rcbt2 = rcb * rt2;
		rcbxu = rcb * 2. * da1dxu;
		rcbyu = rcb * 2. * da1dyu;
		rcbzu = rcb * 2. * da1dzu;
		rcbu2 = rcb * ru2;
		da1dxibt = yca * da1dzt - zca * da1dyt;
		da1dxibu = zdc * da1dyu - ydc * da1dzu;
		da1dyibt = zca * da1dxt - xca * da1dzt;
		da1dyibu = xdc * da1dzu - zdc * da1dxu;
		da1dzibt = xca * da1dyt - yca * da1dxt;
		da1dzibu = ydc * da1dxu - xdc * da1dyu;
		da1dxict = zba * da1dyt - yba * da1dzt;
		da1dxicu = ydb * da1dzu - zdb * da1dyu;
		da1dyict = xba * da1dzt - zba * da1dxt;
		da1dyicu = zdb * da1dxu - xdb * da1dzu;
		da1dzict = yba * da1dxt - xba * da1dyt;
		da1dzicu = xdb * da1dyu - ydb * da1dxu;

/*     chain rule terms for first derivative components */

		da1dxia = zcb * da1dyt - ycb * da1dzt;
		da1dyia = xcb * da1dzt - zcb * da1dxt;
		da1dzia = ycb * da1dxt - xcb * da1dyt;
		da1dxib = da1dxibt + da1dxibu;
		da1dyib = da1dyibt + da1dyibu;
		da1dzib = da1dzibt + da1dzibu;
		da1dxic = da1dxict + da1dxicu;
		da1dyic = da1dyict + da1dyicu;
		da1dzic = da1dzict + da1dzicu;
		da1dxid = zcb * da1dyu - ycb * da1dzu;
		da1dyid = xcb * da1dzu - zcb * da1dxu;
		da1dzid = ycb * da1dxu - xcb * da1dyu;

/*     chain rule terms for second derivative components */

		dxiaxia = rcbxt * da1dxia;
		dxiayia = rcbxt * da1dyia - zcb * rcb / rt2;
		dxiazia = rcbxt * da1dzia + ycb * rcb / rt2;
		dxiaxic = rcbxt * da1dxict + xcb * xt / rcbt2;
		dxiayic = rcbxt * da1dyict - da1dzt - (xba * zcb * xcb + zba *
			 yzcb2) / rcbt2;
		dxiazic = rcbxt * da1dzict + da1dyt + (xba * ycb * xcb + yba *
			 yzcb2) / rcbt2;
		dxiaxid = 0.;
		dxiayid = 0.;
		dxiazid = 0.;
		dyiayia = rcbyt * da1dyia;
		dyiazia = rcbyt * da1dzia - xcb * rcb / rt2;
		dyiaxib = rcbyt * da1dxibt - da1dzt - (yca * zcb * ycb + zca *
			 xzcb2) / rcbt2;
		dyiaxic = rcbyt * da1dxict + da1dzt + (yba * zcb * ycb + zba *
			 xzcb2) / rcbt2;
		dyiayic = rcbyt * da1dyict + ycb * yt / rcbt2;
		dyiazic = rcbyt * da1dzict - da1dxt - (yba * xcb * ycb + xba *
			 xzcb2) / rcbt2;
		dyiaxid = 0.;
		dyiayid = 0.;
		dyiazid = 0.;
		dziazia = rcbzt * da1dzia;
		dziaxib = rcbzt * da1dxibt + da1dyt + (zca * ycb * zcb + yca *
			 xycb2) / rcbt2;
		dziayib = rcbzt * da1dyibt - da1dxt - (zca * xcb * zcb + xca *
			 xycb2) / rcbt2;
		dziaxic = rcbzt * da1dxict - da1dyt - (zba * ycb * zcb + yba *
			 xycb2) / rcbt2;
		dziayic = rcbzt * da1dyict + da1dxt + (zba * xcb * zcb + xba *
			 xycb2) / rcbt2;
		dziazic = rcbzt * da1dzict + zcb * zt / rcbt2;
		dziaxid = 0.;
		dziayid = 0.;
		dziazid = 0.;
		dxibxic = -xcb * da1dxib / (rcb * rcb) - (yca * (zba * xcb + 
			yt) - zca * (yba * xcb - zt)) / rcbt2 - (yt * zba - 
			yba * zt) * 2. * da1dxibt / rt2 - (zdc * (ydb * xcb + 
			zu) - ydc * (zdb * xcb - yu)) / rcbu2 + (yu * zdb - 
			ydb * zu) * 2. * da1dxibu / ru2;
		dxibyic = -ycb * da1dxib / (rcb * rcb) + da1dzt + da1dzu - (
			yca * (zba * ycb - xt) + zca * (xba * xcb + zcb * zba)
			) / rcbt2 - (zt * xba - zba * xt) * 2. * da1dxibt / 
			rt2 + (zdc * (xdb * xcb + zcb * zdb) + ydc * (zdb * 
			ycb + xu)) / rcbu2 + (zu * xdb - zdb * xu) * 2. * 
			da1dxibu / ru2;
		dxibxid = rcbxu * da1dxibu + xcb * xu / rcbu2;
		dxibyid = rcbyu * da1dxibu - da1dzu - (ydc * zcb * ycb + zdc *
			 xzcb2) / rcbu2;
		dxibzid = rcbzu * da1dxibu + da1dyu + (zdc * ycb * zcb + ydc *
			 xycb2) / rcbu2;
		dyibzib = ycb * da1dzib / (rcb * rcb) - (xca * (xca * xcb + 
			zcb * zca) + yca * (ycb * xca + zt)) / rcbt2 - (xt * 
			zca - xca * zt) * 2. * da1dzibt / rt2 + (ydc * (xdc * 
			ycb - zu) + xdc * (xdc * xcb + zcb * zdc)) / rcbu2 + (
			xu * zdc - xdc * zu) * 2. * da1dzibu / ru2;
		dyibxic = -xcb * da1dyib / (rcb * rcb) - da1dzt - da1dzu + (
			xca * (zba * xcb + yt) + zca * (zba * zcb + ycb * yba)
			) / rcbt2 - (yt * zba - yba * zt) * 2. * da1dyibt / 
			rt2 - (zdc * (zdb * zcb + ycb * ydb) + xdc * (zdb * 
			xcb - yu)) / rcbu2 + (yu * zdb - ydb * zu) * 2. * 
			da1dyibu / ru2;
		dyibyic = -ycb * da1dyib / (rcb * rcb) - (zca * (xba * ycb + 
			zt) - xca * (zba * ycb - xt)) / rcbt2 - (zt * xba - 
			zba * xt) * 2. * da1dyibt / rt2 - (xdc * (zdb * ycb + 
			xu) - zdc * (xdb * ycb - zu)) / rcbu2 + (zu * xdb - 
			zdb * xu) * 2. * da1dyibu / ru2;
		dyibxid = rcbxu * da1dyibu + da1dzu + (xdc * zcb * xcb + zdc *
			 yzcb2) / rcbu2;
		dyibyid = rcbyu * da1dyibu + ycb * yu / rcbu2;
		dyibzid = rcbzu * da1dyibu - da1dxu - (zdc * xcb * zcb + xdc *
			 xycb2) / rcbu2;
		dzibxic = -xcb * da1dzib / (rcb * rcb) + da1dyt + da1dyu - (
			xca * (yba * xcb - zt) + yca * (zba * zcb + ycb * yba)
			) / rcbt2 - (yt * zba - yba * zt) * 2. * da1dzibt / 
			rt2 + (ydc * (zdb * zcb + ycb * ydb) + xdc * (ydb * 
			xcb + zu)) / rcbu2 + (yu * zdb - ydb * zu) * 2. * 
			da1dzibu / ru2;
		dzibzic = -zcb * da1dzib / (rcb * rcb) - (xca * (yba * zcb + 
			xt) - yca * (xba * zcb - yt)) / rcbt2 - (xt * yba - 
			xba * yt) * 2. * da1dzibt / rt2 - (ydc * (xdb * zcb + 
			yu) - xdc * (ydb * zcb - xu)) / rcbu2 + (xu * ydb - 
			xdb * yu) * 2. * da1dzibu / ru2;
		dzibxid = rcbxu * da1dzibu - da1dyu - (xdc * ycb * xcb + ydc *
			 yzcb2) / rcbu2;
		dzibyid = rcbyu * da1dzibu + da1dxu + (ydc * xcb * ycb + xdc *
			 xzcb2) / rcbu2;
		dzibzid = rcbzu * da1dzibu + zcb * zu / rcbu2;
		dxicxid = rcbxu * da1dxicu - xcb * (zdb * ycb - ydb * zcb) / 
			rcbu2;
		dxicyid = rcbyu * da1dxicu + da1dzu + (ydb * zcb * ycb + zdb *
			 xzcb2) / rcbu2;
		dxiczid = rcbzu * da1dxicu - da1dyu - (zdb * ycb * zcb + ydb *
			 xycb2) / rcbu2;
		dyicxid = rcbxu * da1dyicu - da1dzu - (xdb * zcb * xcb + zdb *
			 yzcb2) / rcbu2;
		dyicyid = rcbyu * da1dyicu - ycb * (xdb * zcb - zdb * xcb) / 
			rcbu2;
		dyiczid = rcbzu * da1dyicu + da1dxu + (zdb * xcb * zcb + xdb *
			 xycb2) / rcbu2;
		dzicxid = rcbxu * da1dzicu + da1dyu + (xdb * ycb * xcb + ydb *
			 yzcb2) / rcbu2;
		dzicyid = rcbyu * da1dzicu - da1dxu - (ydb * xcb * ycb + xdb *
			 xzcb2) / rcbu2;
		dziczid = rcbzu * da1dzicu - zcb * (ydb * xcb - xdb * ycb) / 
			rcbu2;
		dxidxid = rcbxu * da1dxid;
		dxidyid = rcbxu * da1dyid + zcb * rcb / ru2;
		dxidzid = rcbxu * da1dzid - ycb * rcb / ru2;
		dyidyid = rcbyu * da1dyid;
		dyidzid = rcbyu * da1dzid + xcb * rcb / ru2;
		dzidzid = rcbzu * da1dzid;

/*     get some second derivative chain rule terms by difference */

		dxiaxib = -dxiaxia - dxiaxic - dxiaxid;
		dxiayib = -dxiayia - dxiayic - dxiayid;
		dxiazib = -dxiazia - dxiazic - dxiazid;
		dyiayib = -dyiayia - dyiayic - dyiayid;
		dyiazib = -dyiazia - dyiazic - dyiazid;
		dziazib = -dziazia - dziazic - dziazid;
		dxibxib = -dxiaxib - dxibxic - dxibxid;
		dxibyib = -dyiaxib - dxibyic - dxibyid;
		dxibzib = -dxiazib - dzibxic - dzibxid;
		dxibzic = -dziaxib - dxibzib - dxibzid;
		dyibyib = -dyiayib - dyibyic - dyibyid;
		dyibzic = -dziayib - dyibzib - dyibzid;
		dzibzib = -dziazib - dzibzic - dzibzid;
		dzibyic = -dyiazib - dyibzib - dzibyid;
		dxicxic = -dxiaxic - dxibxic - dxicxid;
		dxicyic = -dyiaxic - dyibxic - dxicyid;
		dxiczic = -dziaxic - dzibxic - dxiczid;
		dyicyic = -dyiayic - dyibyic - dyicyid;
		dyiczic = -dziayic - dzibyic - dyiczid;
		dziczic = -dziazic - dzibzic - dziczid;

/*     increment diagonal and off-diagonal Hessian elements */

		if (*i__ == ia) {
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedang1 * dxiaxia + 
			    d2eda1a1 * da1dxia * da1dxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedang1 * dxiayia + 
			    d2eda1a1 * da1dxia * da1dyia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedang1 * dxiazia + 
			    d2eda1a1 * da1dxia * da1dzia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedang1 * dxiayia + 
			    d2eda1a1 * da1dxia * da1dyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedang1 * dyiayia + 
			    d2eda1a1 * da1dyia * da1dyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedang1 * dyiazia + 
			    d2eda1a1 * da1dyia * da1dzia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedang1 * dxiazia + 
			    d2eda1a1 * da1dxia * da1dzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedang1 * dyiazia + 
			    d2eda1a1 * da1dyia * da1dzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedang1 * dziazia + 
			    d2eda1a1 * da1dzia * da1dzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedang1 * dxiaxib + 
			    d2eda1a1 * da1dxia * da1dxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedang1 * dyiaxib + 
			    d2eda1a1 * da1dyia * da1dxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedang1 * dziaxib + 
			    d2eda1a1 * da1dzia * da1dxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedang1 * dxiayib + 
			    d2eda1a1 * da1dxia * da1dyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedang1 * dyiayib + 
			    d2eda1a1 * da1dyia * da1dyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedang1 * dziayib + 
			    d2eda1a1 * da1dzia * da1dyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedang1 * dxiazib + 
			    d2eda1a1 * da1dxia * da1dzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedang1 * dyiazib + 
			    d2eda1a1 * da1dyia * da1dzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedang1 * dziazib + 
			    d2eda1a1 * da1dzia * da1dzib;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedang1 * dxiaxic + 
			    d2eda1a1 * da1dxia * da1dxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedang1 * dyiaxic + 
			    d2eda1a1 * da1dyia * da1dxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedang1 * dziaxic + 
			    d2eda1a1 * da1dzia * da1dxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedang1 * dxiayic + 
			    d2eda1a1 * da1dxia * da1dyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedang1 * dyiayic + 
			    d2eda1a1 * da1dyia * da1dyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedang1 * dziayic + 
			    d2eda1a1 * da1dzia * da1dyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedang1 * dxiazic + 
			    d2eda1a1 * da1dxia * da1dzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedang1 * dyiazic + 
			    d2eda1a1 * da1dyia * da1dzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedang1 * dziazic + 
			    d2eda1a1 * da1dzia * da1dzic;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedang1 * dxiaxid + 
			    d2eda1a1 * da1dxia * da1dxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedang1 * dyiaxid + 
			    d2eda1a1 * da1dyia * da1dxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedang1 * dziaxid + 
			    d2eda1a1 * da1dzia * da1dxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedang1 * dxiayid + 
			    d2eda1a1 * da1dxia * da1dyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedang1 * dyiayid + 
			    d2eda1a1 * da1dyia * da1dyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedang1 * dziayid + 
			    d2eda1a1 * da1dzia * da1dyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedang1 * dxiazid + 
			    d2eda1a1 * da1dxia * da1dzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedang1 * dyiazid + 
			    d2eda1a1 * da1dyia * da1dzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedang1 * dziazid + 
			    d2eda1a1 * da1dzia * da1dzid;
		} else if (*i__ == ib) {
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedang1 * dxibxib + 
			    d2eda1a1 * da1dxib * da1dxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedang1 * dxibyib + 
			    d2eda1a1 * da1dxib * da1dyib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedang1 * dxibzib + 
			    d2eda1a1 * da1dxib * da1dzib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedang1 * dxibyib + 
			    d2eda1a1 * da1dxib * da1dyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedang1 * dyibyib + 
			    d2eda1a1 * da1dyib * da1dyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedang1 * dyibzib + 
			    d2eda1a1 * da1dyib * da1dzib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedang1 * dxibzib + 
			    d2eda1a1 * da1dxib * da1dzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedang1 * dyibzib + 
			    d2eda1a1 * da1dyib * da1dzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedang1 * dzibzib + 
			    d2eda1a1 * da1dzib * da1dzib;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedang1 * dxiaxib + 
			    d2eda1a1 * da1dxib * da1dxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedang1 * dxiayib + 
			    d2eda1a1 * da1dyib * da1dxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedang1 * dxiazib + 
			    d2eda1a1 * da1dzib * da1dxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedang1 * dyiaxib + 
			    d2eda1a1 * da1dxib * da1dyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedang1 * dyiayib + 
			    d2eda1a1 * da1dyib * da1dyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedang1 * dyiazib + 
			    d2eda1a1 * da1dzib * da1dyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedang1 * dziaxib + 
			    d2eda1a1 * da1dxib * da1dzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedang1 * dziayib + 
			    d2eda1a1 * da1dyib * da1dzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedang1 * dziazib + 
			    d2eda1a1 * da1dzib * da1dzia;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedang1 * dxibxic + 
			    d2eda1a1 * da1dxib * da1dxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedang1 * dyibxic + 
			    d2eda1a1 * da1dyib * da1dxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedang1 * dzibxic + 
			    d2eda1a1 * da1dzib * da1dxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedang1 * dxibyic + 
			    d2eda1a1 * da1dxib * da1dyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedang1 * dyibyic + 
			    d2eda1a1 * da1dyib * da1dyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedang1 * dzibyic + 
			    d2eda1a1 * da1dzib * da1dyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedang1 * dxibzic + 
			    d2eda1a1 * da1dxib * da1dzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedang1 * dyibzic + 
			    d2eda1a1 * da1dyib * da1dzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedang1 * dzibzic + 
			    d2eda1a1 * da1dzib * da1dzic;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedang1 * dxibxid + 
			    d2eda1a1 * da1dxib * da1dxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedang1 * dyibxid + 
			    d2eda1a1 * da1dyib * da1dxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedang1 * dzibxid + 
			    d2eda1a1 * da1dzib * da1dxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedang1 * dxibyid + 
			    d2eda1a1 * da1dxib * da1dyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedang1 * dyibyid + 
			    d2eda1a1 * da1dyib * da1dyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedang1 * dzibyid + 
			    d2eda1a1 * da1dzib * da1dyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedang1 * dxibzid + 
			    d2eda1a1 * da1dxib * da1dzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedang1 * dyibzid + 
			    d2eda1a1 * da1dyib * da1dzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedang1 * dzibzid + 
			    d2eda1a1 * da1dzib * da1dzid;
		} else if (*i__ == ic) {
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedang1 * dxicxic + 
			    d2eda1a1 * da1dxic * da1dxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedang1 * dxicyic + 
			    d2eda1a1 * da1dxic * da1dyic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedang1 * dxiczic + 
			    d2eda1a1 * da1dxic * da1dzic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedang1 * dxicyic + 
			    d2eda1a1 * da1dxic * da1dyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedang1 * dyicyic + 
			    d2eda1a1 * da1dyic * da1dyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedang1 * dyiczic + 
			    d2eda1a1 * da1dyic * da1dzic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedang1 * dxiczic + 
			    d2eda1a1 * da1dxic * da1dzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedang1 * dyiczic + 
			    d2eda1a1 * da1dyic * da1dzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedang1 * dziczic + 
			    d2eda1a1 * da1dzic * da1dzic;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedang1 * dxiaxic + 
			    d2eda1a1 * da1dxic * da1dxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedang1 * dxiayic + 
			    d2eda1a1 * da1dyic * da1dxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedang1 * dxiazic + 
			    d2eda1a1 * da1dzic * da1dxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedang1 * dyiaxic + 
			    d2eda1a1 * da1dxic * da1dyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedang1 * dyiayic + 
			    d2eda1a1 * da1dyic * da1dyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedang1 * dyiazic + 
			    d2eda1a1 * da1dzic * da1dyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedang1 * dziaxic + 
			    d2eda1a1 * da1dxic * da1dzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedang1 * dziayic + 
			    d2eda1a1 * da1dyic * da1dzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedang1 * dziazic + 
			    d2eda1a1 * da1dzic * da1dzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedang1 * dxibxic + 
			    d2eda1a1 * da1dxic * da1dxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedang1 * dxibyic + 
			    d2eda1a1 * da1dyic * da1dxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedang1 * dxibzic + 
			    d2eda1a1 * da1dzic * da1dxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedang1 * dyibxic + 
			    d2eda1a1 * da1dxic * da1dyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedang1 * dyibyic + 
			    d2eda1a1 * da1dyic * da1dyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedang1 * dyibzic + 
			    d2eda1a1 * da1dzic * da1dyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedang1 * dzibxic + 
			    d2eda1a1 * da1dxic * da1dzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedang1 * dzibyic + 
			    d2eda1a1 * da1dyic * da1dzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedang1 * dzibzic + 
			    d2eda1a1 * da1dzic * da1dzib;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedang1 * dxicxid + 
			    d2eda1a1 * da1dxic * da1dxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedang1 * dyicxid + 
			    d2eda1a1 * da1dyic * da1dxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedang1 * dzicxid + 
			    d2eda1a1 * da1dzic * da1dxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedang1 * dxicyid + 
			    d2eda1a1 * da1dxic * da1dyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedang1 * dyicyid + 
			    d2eda1a1 * da1dyic * da1dyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedang1 * dzicyid + 
			    d2eda1a1 * da1dzic * da1dyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedang1 * dxiczid + 
			    d2eda1a1 * da1dxic * da1dzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedang1 * dyiczid + 
			    d2eda1a1 * da1dyic * da1dzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedang1 * dziczid + 
			    d2eda1a1 * da1dzic * da1dzid;
		} else if (*i__ == id) {
		    hessx_ref(1, id) = hessx_ref(1, id) + dedang1 * dxidxid + 
			    d2eda1a1 * da1dxid * da1dxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedang1 * dxidyid + 
			    d2eda1a1 * da1dxid * da1dyid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedang1 * dxidzid + 
			    d2eda1a1 * da1dxid * da1dzid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedang1 * dxidyid + 
			    d2eda1a1 * da1dxid * da1dyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedang1 * dyidyid + 
			    d2eda1a1 * da1dyid * da1dyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedang1 * dyidzid + 
			    d2eda1a1 * da1dyid * da1dzid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedang1 * dxidzid + 
			    d2eda1a1 * da1dxid * da1dzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedang1 * dyidzid + 
			    d2eda1a1 * da1dyid * da1dzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedang1 * dzidzid + 
			    d2eda1a1 * da1dzid * da1dzid;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedang1 * dxiaxid + 
			    d2eda1a1 * da1dxid * da1dxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedang1 * dxiayid + 
			    d2eda1a1 * da1dyid * da1dxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedang1 * dxiazid + 
			    d2eda1a1 * da1dzid * da1dxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedang1 * dyiaxid + 
			    d2eda1a1 * da1dxid * da1dyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedang1 * dyiayid + 
			    d2eda1a1 * da1dyid * da1dyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedang1 * dyiazid + 
			    d2eda1a1 * da1dzid * da1dyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedang1 * dziaxid + 
			    d2eda1a1 * da1dxid * da1dzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedang1 * dziayid + 
			    d2eda1a1 * da1dyid * da1dzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedang1 * dziazid + 
			    d2eda1a1 * da1dzid * da1dzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedang1 * dxibxid + 
			    d2eda1a1 * da1dxid * da1dxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedang1 * dxibyid + 
			    d2eda1a1 * da1dyid * da1dxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedang1 * dxibzid + 
			    d2eda1a1 * da1dzid * da1dxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedang1 * dyibxid + 
			    d2eda1a1 * da1dxid * da1dyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedang1 * dyibyid + 
			    d2eda1a1 * da1dyid * da1dyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedang1 * dyibzid + 
			    d2eda1a1 * da1dzid * da1dyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedang1 * dzibxid + 
			    d2eda1a1 * da1dxid * da1dzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedang1 * dzibyid + 
			    d2eda1a1 * da1dyid * da1dzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedang1 * dzibzid + 
			    d2eda1a1 * da1dzid * da1dzib;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedang1 * dxicxid + 
			    d2eda1a1 * da1dxid * da1dxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedang1 * dxicyid + 
			    d2eda1a1 * da1dyid * da1dxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedang1 * dxiczid + 
			    d2eda1a1 * da1dzid * da1dxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedang1 * dyicxid + 
			    d2eda1a1 * da1dxid * da1dyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedang1 * dyicyid + 
			    d2eda1a1 * da1dyid * da1dyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedang1 * dyiczid + 
			    d2eda1a1 * da1dzid * da1dyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedang1 * dzicxid + 
			    d2eda1a1 * da1dxid * da1dzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedang1 * dzicyid + 
			    d2eda1a1 * da1dyid * da1dzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedang1 * dziczid + 
			    d2eda1a1 * da1dzid * da1dzic;
		}

/*     abbreviations for first derivative chain rule terms */

		xec = xie - xic;
		yec = yie - yic;
		zec = zie - zic;
		if (bound_1.use_polymer__) {
		    image_(&xec, &yec, &zec);
		}
		da2dxu = (yu * zdc - ydc * zu) / (ru2 * rdc);
		da2dyu = (zu * xdc - zdc * xu) / (ru2 * rdc);
		da2dzu = (xu * ydc - xdc * yu) / (ru2 * rdc);
		da2dxv = -(yv * zdc - ydc * zv) / (rv2 * rdc);
		da2dyv = -(zv * xdc - zdc * xv) / (rv2 * rdc);
		da2dzv = -(xv * ydc - xdc * yv) / (rv2 * rdc);

/*     abbreviations for second derivative chain rule terms */

		xydc2 = xdc * xdc + ydc * ydc;
		xzdc2 = xdc * xdc + zdc * zdc;
		yzdc2 = ydc * ydc + zdc * zdc;
		rdcxu = rdc * -2. * da2dxu;
		rdcyu = rdc * -2. * da2dyu;
		rdczu = rdc * -2. * da2dzu;
		rdcu2 = rdc * ru2;
		rdcxv = rdc * 2. * da2dxv;
		rdcyv = rdc * 2. * da2dyv;
		rdczv = rdc * 2. * da2dzv;
		rdcv2 = rdc * rv2;
		da2dxicu = ydb * da2dzu - zdb * da2dyu;
		da2dxicv = zed * da2dyv - yed * da2dzv;
		da2dyicu = zdb * da2dxu - xdb * da2dzu;
		da2dyicv = xed * da2dzv - zed * da2dxv;
		da2dzicu = xdb * da2dyu - ydb * da2dxu;
		da2dzicv = yed * da2dxv - xed * da2dyv;
		da2dxidu = zcb * da2dyu - ycb * da2dzu;
		da2dxidv = yec * da2dzv - zec * da2dyv;
		da2dyidu = xcb * da2dzu - zcb * da2dxu;
		da2dyidv = zec * da2dxv - xec * da2dzv;
		da2dzidu = ycb * da2dxu - xcb * da2dyu;
		da2dzidv = xec * da2dyv - yec * da2dxv;

/*     chain rule terms for first derivative components */

		da2dxib = zdc * da2dyu - ydc * da2dzu;
		da2dyib = xdc * da2dzu - zdc * da2dxu;
		da2dzib = ydc * da2dxu - xdc * da2dyu;
		da2dxic = da2dxicu + da2dxicv;
		da2dyic = da2dyicu + da2dyicv;
		da2dzic = da2dzicu + da2dzicv;
		da2dxid = da2dxidu + da2dxidv;
		da2dyid = da2dyidu + da2dyidv;
		da2dzid = da2dzidu + da2dzidv;
		da2dxie = zdc * da2dyv - ydc * da2dzv;
		da2dyie = xdc * da2dzv - zdc * da2dxv;
		da2dzie = ydc * da2dxv - xdc * da2dyv;

/*     chain rule terms for second derivative components */

		dxibxib2 = rdcxu * da2dxib;
		dxibyib2 = rdcxu * da2dyib - zdc * rdc / ru2;
		dxibzib2 = rdcxu * da2dzib + ydc * rdc / ru2;
		dxibxid2 = rdcxu * da2dxidu + xdc * xu / rdcu2;
		dxibyid2 = rdcxu * da2dyidu - da2dzu - (xcb * zdc * xdc + zcb 
			* yzdc2) / rdcu2;
		dxibzid2 = rdcxu * da2dzidu + da2dyu + (xcb * ydc * xdc + ycb 
			* yzdc2) / rdcu2;
		dxibxie2 = 0.;
		dxibyie2 = 0.;
		dxibzie2 = 0.;
		dyibyib2 = rdcyu * da2dyib;
		dyibzib2 = rdcyu * da2dzib - xdc * rdc / ru2;
		dyibxic2 = rdcyu * da2dxicu - da2dzu - (ydb * zdc * ydc + zdb 
			* xzdc2) / rdcu2;
		dyibxid2 = rdcyu * da2dxidu + da2dzu + (ycb * zdc * ydc + zcb 
			* xzdc2) / rdcu2;
		dyibyid2 = rdcyu * da2dyidu + ydc * yu / rdcu2;
		dyibzid2 = rdcyu * da2dzidu - da2dxu - (ycb * xdc * ydc + xcb 
			* xzdc2) / rdcu2;
		dyibxie2 = 0.;
		dyibyie2 = 0.;
		dyibzie2 = 0.;
		dzibzib2 = rdczu * da2dzib;
		dzibxic2 = rdczu * da2dxicu + da2dyu + (zdb * ydc * zdc + ydb 
			* xydc2) / rdcu2;
		dzibyic2 = rdczu * da2dyicu - da2dxu - (zdb * xdc * zdc + xdb 
			* xydc2) / rdcu2;
		dzibxid2 = rdczu * da2dxidu - da2dyu - (zcb * ydc * zdc + ycb 
			* xydc2) / rdcu2;
		dzibyid2 = rdczu * da2dyidu + da2dxu + (zcb * xdc * zdc + xcb 
			* xydc2) / rdcu2;
		dzibzid2 = rdczu * da2dzidu + zdc * zu / rdcu2;
		dzibxie2 = 0.;
		dzibyie2 = 0.;
		dzibzie2 = 0.;
		dxicxid2 = -xdc * da2dxic / (rdc * rdc) - (ydb * (zcb * xdc + 
			yu) - zdb * (ycb * xdc - zu)) / rdcu2 - (yu * zcb - 
			ycb * zu) * 2. * da2dxicu / ru2 - (zed * (yec * xdc + 
			zv) - yed * (zec * xdc - yv)) / rdcv2 + (yv * zec - 
			yec * zv) * 2. * da2dxicv / rv2;
		dxicyid2 = -ydc * da2dxic / (rdc * rdc) + da2dzu + da2dzv - (
			ydb * (zcb * ydc - xu) + zdb * (xcb * xdc + zdc * zcb)
			) / rdcu2 - (zu * xcb - zcb * xu) * 2. * da2dxicu / 
			ru2 + (zed * (xec * xdc + zdc * zec) + yed * (zec * 
			ydc + xv)) / rdcv2 + (zv * xec - zec * xv) * 2. * 
			da2dxicv / rv2;
		dxicxie2 = rdcxv * da2dxicv + xdc * xv / rdcv2;
		dxicyie2 = rdcyv * da2dxicv - da2dzv - (yed * zdc * ydc + zed 
			* xzdc2) / rdcv2;
		dxiczie2 = rdczv * da2dxicv + da2dyv + (zed * ydc * zdc + yed 
			* xydc2) / rdcv2;
		dyiczic2 = ydc * da2dzic / (rdc * rdc) - (xdb * (xdb * xdc + 
			zdc * zdb) + ydb * (ydc * xdb + zu)) / rdcu2 - (xu * 
			zdb - xdb * zu) * 2. * da2dzicu / ru2 + (yed * (xed * 
			ydc - zv) + xed * (xed * xdc + zdc * zed)) / rdcv2 + (
			xv * zed - xed * zv) * 2. * da2dzicv / rv2;
		dyicxid2 = -xdc * da2dyic / (rdc * rdc) - da2dzu - da2dzv + (
			xdb * (zcb * xdc + yu) + zdb * (zcb * zdc + ydc * ycb)
			) / rdcu2 - (yu * zcb - ycb * zu) * 2. * da2dyicu / 
			ru2 - (zed * (zec * zdc + ydc * yec) + xed * (zec * 
			xdc - yv)) / rdcv2 + (yv * zec - yec * zv) * 2. * 
			da2dyicv / rv2;
		dyicyid2 = -ydc * da2dyic / (rdc * rdc) - (zdb * (xcb * ydc + 
			zu) - xdb * (zcb * ydc - xu)) / rdcu2 - (zu * xcb - 
			zcb * xu) * 2. * da2dyicu / ru2 - (xed * (zec * ydc + 
			xv) - zed * (xec * ydc - zv)) / rdcv2 + (zv * xec - 
			zec * xv) * 2. * da2dyicv / rv2;
		dyicxie2 = rdcxv * da2dyicv + da2dzv + (xed * zdc * xdc + zed 
			* yzdc2) / rdcv2;
		dyicyie2 = rdcyv * da2dyicv + ydc * yv / rdcv2;
		dyiczie2 = rdczv * da2dyicv - da2dxv - (zed * xdc * zdc + xed 
			* xydc2) / rdcv2;
		dzicxid2 = -xdc * da2dzic / (rdc * rdc) + da2dyu + da2dyv - (
			xdb * (ycb * xdc - zu) + ydb * (zcb * zdc + ydc * ycb)
			) / rdcu2 - (yu * zcb - ycb * zu) * 2. * da2dzicu / 
			ru2 + (yed * (zec * zdc + ydc * yec) + xed * (yec * 
			xdc + zv)) / rdcv2 + (yv * zec - yec * zv) * 2. * 
			da2dzicv / rv2;
		dziczid2 = -zdc * da2dzic / (rdc * rdc) - (xdb * (ycb * zdc + 
			xu) - ydb * (xcb * zdc - yu)) / rdcu2 - (xu * ycb - 
			xcb * yu) * 2. * da2dzicu / ru2 - (yed * (xec * zdc + 
			yv) - xed * (yec * zdc - xv)) / rdcv2 + (xv * yec - 
			xec * yv) * 2. * da2dzicv / rv2;
		dzicxie2 = rdcxv * da2dzicv - da2dyv - (xed * ydc * xdc + yed 
			* yzdc2) / rdcv2;
		dzicyie2 = rdcyv * da2dzicv + da2dxv + (yed * xdc * ydc + xed 
			* xzdc2) / rdcv2;
		dziczie2 = rdczv * da2dzicv + zdc * zv / rdcv2;
		dxidxie2 = rdcxv * da2dxidv - xdc * (zec * ydc - yec * zdc) / 
			rdcv2;
		dxidyie2 = rdcyv * da2dxidv + da2dzv + (yec * zdc * ydc + zec 
			* xzdc2) / rdcv2;
		dxidzie2 = rdczv * da2dxidv - da2dyv - (zec * ydc * zdc + yec 
			* xydc2) / rdcv2;
		dyidxie2 = rdcxv * da2dyidv - da2dzv - (xec * zdc * xdc + zec 
			* yzdc2) / rdcv2;
		dyidyie2 = rdcyv * da2dyidv - ydc * (xec * zdc - zec * xdc) / 
			rdcv2;
		dyidzie2 = rdczv * da2dyidv + da2dxv + (zec * xdc * zdc + xec 
			* xydc2) / rdcv2;
		dzidxie2 = rdcxv * da2dzidv + da2dyv + (xec * ydc * xdc + yec 
			* yzdc2) / rdcv2;
		dzidyie2 = rdcyv * da2dzidv - da2dxv - (yec * xdc * ydc + xec 
			* xzdc2) / rdcv2;
		dzidzie2 = rdczv * da2dzidv - zdc * (yec * xdc - xec * ydc) / 
			rdcv2;
		dxiexie2 = rdcxv * da2dxie;
		dxieyie2 = rdcxv * da2dyie + zdc * rdc / rv2;
		dxiezie2 = rdcxv * da2dzie - ydc * rdc / rv2;
		dyieyie2 = rdcyv * da2dyie;
		dyiezie2 = rdcyv * da2dzie + xdc * rdc / rv2;
		dziezie2 = rdczv * da2dzie;

/*     get some second derivative chain rule terms by difference */

		dxibxic2 = -dxibxib2 - dxibxid2 - dxibxie2;
		dxibyic2 = -dxibyib2 - dxibyid2 - dxibyie2;
		dxibzic2 = -dxibzib2 - dxibzid2 - dxibzie2;
		dyibyic2 = -dyibyib2 - dyibyid2 - dyibyie2;
		dyibzic2 = -dyibzib2 - dyibzid2 - dyibzie2;
		dzibzic2 = -dzibzib2 - dzibzid2 - dzibzie2;
		dxicxic2 = -dxibxic2 - dxicxid2 - dxicxie2;
		dxicyic2 = -dyibxic2 - dxicyid2 - dxicyie2;
		dxiczic2 = -dxibzic2 - dzicxid2 - dzicxie2;
		dxiczid2 = -dzibxic2 - dxiczic2 - dxiczie2;
		dyicyic2 = -dyibyic2 - dyicyid2 - dyicyie2;
		dyiczid2 = -dzibyic2 - dyiczic2 - dyiczie2;
		dziczic2 = -dzibzic2 - dziczid2 - dziczie2;
		dzicyid2 = -dyibzic2 - dyiczic2 - dzicyie2;
		dxidxid2 = -dxibxid2 - dxicxid2 - dxidxie2;
		dxidyid2 = -dyibxid2 - dyicxid2 - dxidyie2;
		dxidzid2 = -dzibxid2 - dzicxid2 - dxidzie2;
		dyidyid2 = -dyibyid2 - dyicyid2 - dyidyie2;
		dyidzid2 = -dzibyid2 - dzicyid2 - dyidzie2;
		dzidzid2 = -dzibzid2 - dziczid2 - dzidzie2;

/*     increment diagonal and off-diagonal Hessian elements */

		if (*i__ == ia) {
		    hessx_ref(1, ib) = hessx_ref(1, ib) + d2eda1a2 * da1dxia *
			     da2dxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + d2eda1a2 * da1dyia *
			     da2dxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + d2eda1a2 * da1dzia *
			     da2dxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + d2eda1a2 * da1dxia *
			     da2dyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + d2eda1a2 * da1dyia *
			     da2dyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + d2eda1a2 * da1dzia *
			     da2dyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + d2eda1a2 * da1dxia *
			     da2dzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + d2eda1a2 * da1dyia *
			     da2dzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + d2eda1a2 * da1dzia *
			     da2dzib;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + d2eda1a2 * da1dxia *
			     da2dxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + d2eda1a2 * da1dyia *
			     da2dxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + d2eda1a2 * da1dzia *
			     da2dxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + d2eda1a2 * da1dxia *
			     da2dyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + d2eda1a2 * da1dyia *
			     da2dyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + d2eda1a2 * da1dzia *
			     da2dyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + d2eda1a2 * da1dxia *
			     da2dzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + d2eda1a2 * da1dyia *
			     da2dzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + d2eda1a2 * da1dzia *
			     da2dzic;
		    hessx_ref(1, id) = hessx_ref(1, id) + d2eda1a2 * da1dxia *
			     da2dxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + d2eda1a2 * da1dyia *
			     da2dxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + d2eda1a2 * da1dzia *
			     da2dxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + d2eda1a2 * da1dxia *
			     da2dyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + d2eda1a2 * da1dyia *
			     da2dyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + d2eda1a2 * da1dzia *
			     da2dyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + d2eda1a2 * da1dxia *
			     da2dzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + d2eda1a2 * da1dyia *
			     da2dzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + d2eda1a2 * da1dzia *
			     da2dzid;
		    hessx_ref(1, ie) = hessx_ref(1, ie) + d2eda1a2 * da1dxia *
			     da2dxie;
		    hessy_ref(1, ie) = hessy_ref(1, ie) + d2eda1a2 * da1dyia *
			     da2dxie;
		    hessz_ref(1, ie) = hessz_ref(1, ie) + d2eda1a2 * da1dzia *
			     da2dxie;
		    hessx_ref(2, ie) = hessx_ref(2, ie) + d2eda1a2 * da1dxia *
			     da2dyie;
		    hessy_ref(2, ie) = hessy_ref(2, ie) + d2eda1a2 * da1dyia *
			     da2dyie;
		    hessz_ref(2, ie) = hessz_ref(2, ie) + d2eda1a2 * da1dzia *
			     da2dyie;
		    hessx_ref(3, ie) = hessx_ref(3, ie) + d2eda1a2 * da1dxia *
			     da2dzie;
		    hessy_ref(3, ie) = hessy_ref(3, ie) + d2eda1a2 * da1dyia *
			     da2dzie;
		    hessz_ref(3, ie) = hessz_ref(3, ie) + d2eda1a2 * da1dzia *
			     da2dzie;
		} else if (*i__ == ib) {
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedang2 * dxibxib2 
			    + d2eda2a2 * da2dxib * da2dxib + d2eda1a2 * 2. * 
			    da1dxib * da2dxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedang2 * dxibyib2 
			    + d2eda2a2 * da2dxib * da2dyib + d2eda1a2 * 
			    da1dxib * da2dyib + d2eda1a2 * da2dxib * da1dyib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedang2 * dxibzib2 
			    + d2eda2a2 * da2dxib * da2dzib + d2eda1a2 * 
			    da1dxib * da2dzib + d2eda1a2 * da2dxib * da1dzib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedang2 * dxibyib2 
			    + d2eda2a2 * da2dxib * da2dyib + d2eda1a2 * 
			    da1dxib * da2dyib + d2eda1a2 * da2dxib * da1dyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedang2 * dyibyib2 
			    + d2eda2a2 * da2dyib * da2dyib + d2eda1a2 * 
			    da1dyib * da2dyib + d2eda1a2 * da2dyib * da1dyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedang2 * dyibzib2 
			    + d2eda2a2 * da2dyib * da2dzib + d2eda1a2 * 
			    da1dyib * da2dzib + d2eda1a2 * da2dyib * da1dzib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedang2 * dxibzib2 
			    + d2eda2a2 * da2dxib * da2dzib + d2eda1a2 * 
			    da1dxib * da2dzib + d2eda1a2 * da2dxib * da1dzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedang2 * dyibzib2 
			    + d2eda2a2 * da2dyib * da2dzib + d2eda1a2 * 
			    da1dyib * da2dzib + d2eda1a2 * da2dyib * da1dzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedang2 * dzibzib2 
			    + d2eda2a2 * da2dzib * da2dzib + d2eda1a2 * 
			    da1dzib * da2dzib + d2eda1a2 * da2dzib * da1dzib;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + d2eda1a2 * da2dxib *
			     da1dxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + d2eda1a2 * da2dyib *
			     da1dxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + d2eda1a2 * da2dzib *
			     da1dxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + d2eda1a2 * da2dxib *
			     da1dyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + d2eda1a2 * da2dyib *
			     da1dyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + d2eda1a2 * da2dzib *
			     da1dyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + d2eda1a2 * da2dxib *
			     da1dzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + d2eda1a2 * da2dyib *
			     da1dzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + d2eda1a2 * da2dzib *
			     da1dzia;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedang2 * dxibxic2 
			    + d2eda2a2 * da2dxib * da2dxic + d2eda1a2 * 
			    da1dxib * da2dxic + d2eda1a2 * da2dxib * da1dxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedang2 * dyibxic2 
			    + d2eda2a2 * da2dyib * da2dxic + d2eda1a2 * 
			    da1dyib * da2dxic + d2eda1a2 * da2dyib * da1dxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedang2 * dzibxic2 
			    + d2eda2a2 * da2dzib * da2dxic + d2eda1a2 * 
			    da1dzib * da2dxic + d2eda1a2 * da2dzib * da1dxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedang2 * dxibyic2 
			    + d2eda2a2 * da2dxib * da2dyic + d2eda1a2 * 
			    da1dxib * da2dyic + d2eda1a2 * da2dxib * da1dyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedang2 * dyibyic2 
			    + d2eda2a2 * da2dyib * da2dyic + d2eda1a2 * 
			    da1dyib * da2dyic + d2eda1a2 * da2dyib * da1dyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedang2 * dzibyic2 
			    + d2eda2a2 * da2dzib * da2dyic + d2eda1a2 * 
			    da1dzib * da2dyic + d2eda1a2 * da2dzib * da1dyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedang2 * dxibzic2 
			    + d2eda2a2 * da2dxib * da2dzic + d2eda1a2 * 
			    da1dxib * da2dzic + d2eda1a2 * da2dxib * da1dzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedang2 * dyibzic2 
			    + d2eda2a2 * da2dyib * da2dzic + d2eda1a2 * 
			    da1dyib * da2dzic + d2eda1a2 * da2dyib * da1dzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedang2 * dzibzic2 
			    + d2eda2a2 * da2dzib * da2dzic + d2eda1a2 * 
			    da1dzib * da2dzic + d2eda1a2 * da2dzib * da1dzic;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedang2 * dxibxid2 
			    + d2eda2a2 * da2dxib * da2dxid + d2eda1a2 * 
			    da1dxib * da2dxid + d2eda1a2 * da2dxib * da1dxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedang2 * dyibxid2 
			    + d2eda2a2 * da2dyib * da2dxid + d2eda1a2 * 
			    da1dyib * da2dxid + d2eda1a2 * da2dyib * da1dxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedang2 * dzibxid2 
			    + d2eda2a2 * da2dzib * da2dxid + d2eda1a2 * 
			    da1dzib * da2dxid + d2eda1a2 * da2dzib * da1dxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedang2 * dxibyid2 
			    + d2eda2a2 * da2dxib * da2dyid + d2eda1a2 * 
			    da1dxib * da2dyid + d2eda1a2 * da2dxib * da1dyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedang2 * dyibyid2 
			    + d2eda2a2 * da2dyib * da2dyid + d2eda1a2 * 
			    da1dyib * da2dyid + d2eda1a2 * da2dyib * da1dyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedang2 * dzibyid2 
			    + d2eda2a2 * da2dzib * da2dyid + d2eda1a2 * 
			    da1dzib * da2dyid + d2eda1a2 * da2dzib * da1dyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedang2 * dxibzid2 
			    + d2eda2a2 * da2dxib * da2dzid + d2eda1a2 * 
			    da1dxib * da2dzid + d2eda1a2 * da2dxib * da1dzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedang2 * dyibzid2 
			    + d2eda2a2 * da2dyib * da2dzid + d2eda1a2 * 
			    da1dyib * da2dzid + d2eda1a2 * da2dyib * da1dzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedang2 * dzibzid2 
			    + d2eda2a2 * da2dzib * da2dzid + d2eda1a2 * 
			    da1dzib * da2dzid + d2eda1a2 * da2dzib * da1dzid;
		    hessx_ref(1, ie) = hessx_ref(1, ie) + dedang2 * dxibxie2 
			    + d2eda2a2 * da2dxib * da2dxie + d2eda1a2 * 
			    da1dxib * da2dxie;
		    hessy_ref(1, ie) = hessy_ref(1, ie) + dedang2 * dyibxie2 
			    + d2eda2a2 * da2dyib * da2dxie + d2eda1a2 * 
			    da1dyib * da2dxie;
		    hessz_ref(1, ie) = hessz_ref(1, ie) + dedang2 * dzibxie2 
			    + d2eda2a2 * da2dzib * da2dxie + d2eda1a2 * 
			    da1dzib * da2dxie;
		    hessx_ref(2, ie) = hessx_ref(2, ie) + dedang2 * dxibyie2 
			    + d2eda2a2 * da2dxib * da2dyie + d2eda1a2 * 
			    da1dxib * da2dyie;
		    hessy_ref(2, ie) = hessy_ref(2, ie) + dedang2 * dyibyie2 
			    + d2eda2a2 * da2dyib * da2dyie + d2eda1a2 * 
			    da1dyib * da2dyie;
		    hessz_ref(2, ie) = hessz_ref(2, ie) + dedang2 * dzibyie2 
			    + d2eda2a2 * da2dzib * da2dyie + d2eda1a2 * 
			    da1dzib * da2dyie;
		    hessx_ref(3, ie) = hessx_ref(3, ie) + dedang2 * dxibzie2 
			    + d2eda2a2 * da2dxib * da2dzie + d2eda1a2 * 
			    da1dxib * da2dzie;
		    hessy_ref(3, ie) = hessy_ref(3, ie) + dedang2 * dyibzie2 
			    + d2eda2a2 * da2dyib * da2dzie + d2eda1a2 * 
			    da1dyib * da2dzie;
		    hessz_ref(3, ie) = hessz_ref(3, ie) + dedang2 * dzibzie2 
			    + d2eda2a2 * da2dzib * da2dzie + d2eda1a2 * 
			    da1dzib * da2dzie;
		} else if (*i__ == ic) {
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedang2 * dxicxic2 
			    + d2eda2a2 * da2dxic * da2dxic + d2eda1a2 * 
			    da1dxic * da2dxic + d2eda1a2 * da2dxic * da1dxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedang2 * dxicyic2 
			    + d2eda2a2 * da2dxic * da2dyic + d2eda1a2 * 
			    da1dxic * da2dyic + d2eda1a2 * da2dxic * da1dyic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedang2 * dxiczic2 
			    + d2eda2a2 * da2dxic * da2dzic + d2eda1a2 * 
			    da1dxic * da2dzic + d2eda1a2 * da2dxic * da1dzic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedang2 * dxicyic2 
			    + d2eda2a2 * da2dxic * da2dyic + d2eda1a2 * 
			    da1dxic * da2dyic + d2eda1a2 * da2dxic * da1dyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedang2 * dyicyic2 
			    + d2eda2a2 * da2dyic * da2dyic + d2eda1a2 * 
			    da1dyic * da2dyic + d2eda1a2 * da2dyic * da1dyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedang2 * dyiczic2 
			    + d2eda2a2 * da2dyic * da2dzic + d2eda1a2 * 
			    da1dyic * da2dzic + d2eda1a2 * da2dyic * da1dzic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedang2 * dxiczic2 
			    + d2eda2a2 * da2dxic * da2dzic + d2eda1a2 * 
			    da1dxic * da2dzic + d2eda1a2 * da2dxic * da1dzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedang2 * dyiczic2 
			    + d2eda2a2 * da2dyic * da2dzic + d2eda1a2 * 
			    da1dyic * da2dzic + d2eda1a2 * da2dyic * da1dzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedang2 * dziczic2 
			    + d2eda2a2 * da2dzic * da2dzic + d2eda1a2 * 
			    da1dzic * da2dzic + d2eda1a2 * da2dzic * da1dzic;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + d2eda1a2 * da2dxic *
			     da1dxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + d2eda1a2 * da2dyic *
			     da1dxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + d2eda1a2 * da2dzic *
			     da1dxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + d2eda1a2 * da2dxic *
			     da1dyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + d2eda1a2 * da2dyic *
			     da1dyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + d2eda1a2 * da2dzic *
			     da1dyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + d2eda1a2 * da2dxic *
			     da1dzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + d2eda1a2 * da2dyic *
			     da1dzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + d2eda1a2 * da2dzic *
			     da1dzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedang2 * dxibxic2 
			    + d2eda2a2 * da2dxic * da2dxib + d2eda1a2 * 
			    da1dxic * da2dxib + d2eda1a2 * da2dxic * da1dxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedang2 * dxibyic2 
			    + d2eda2a2 * da2dyic * da2dxib + d2eda1a2 * 
			    da1dyic * da2dxib + d2eda1a2 * da2dyic * da1dxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedang2 * dxibzic2 
			    + d2eda2a2 * da2dzic * da2dxib + d2eda1a2 * 
			    da1dzic * da2dxib + d2eda1a2 * da2dzic * da1dxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedang2 * dyibxic2 
			    + d2eda2a2 * da2dxic * da2dyib + d2eda1a2 * 
			    da1dxic * da2dyib + d2eda1a2 * da2dxic * da1dyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedang2 * dyibyic2 
			    + d2eda2a2 * da2dyic * da2dyib + d2eda1a2 * 
			    da1dyic * da2dyib + d2eda1a2 * da2dyic * da1dyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedang2 * dyibzic2 
			    + d2eda2a2 * da2dzic * da2dyib + d2eda1a2 * 
			    da1dzic * da2dyib + d2eda1a2 * da2dzic * da1dyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedang2 * dzibxic2 
			    + d2eda2a2 * da2dxic * da2dzib + d2eda1a2 * 
			    da1dxic * da2dzib + d2eda1a2 * da2dxic * da1dzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedang2 * dzibyic2 
			    + d2eda2a2 * da2dyic * da2dzib + d2eda1a2 * 
			    da1dyic * da2dzib + d2eda1a2 * da2dyic * da1dzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedang2 * dzibzic2 
			    + d2eda2a2 * da2dzic * da2dzib + d2eda1a2 * 
			    da1dzic * da2dzib + d2eda1a2 * da2dzic * da1dzib;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedang2 * dxicxid2 
			    + d2eda2a2 * da2dxic * da2dxid + d2eda1a2 * 
			    da1dxic * da2dxid + d2eda1a2 * da2dxic * da1dxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedang2 * dyicxid2 
			    + d2eda2a2 * da2dyic * da2dxid + d2eda1a2 * 
			    da1dyic * da2dxid + d2eda1a2 * da2dyic * da1dxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedang2 * dzicxid2 
			    + d2eda2a2 * da2dzic * da2dxid + d2eda1a2 * 
			    da1dzic * da2dxid + d2eda1a2 * da2dzic * da1dxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedang2 * dxicyid2 
			    + d2eda2a2 * da2dxic * da2dyid + d2eda1a2 * 
			    da1dxic * da2dyid + d2eda1a2 * da2dxic * da1dyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedang2 * dyicyid2 
			    + d2eda2a2 * da2dyic * da2dyid + d2eda1a2 * 
			    da1dyic * da2dyid + d2eda1a2 * da2dyic * da1dyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedang2 * dzicyid2 
			    + d2eda2a2 * da2dzic * da2dyid + d2eda1a2 * 
			    da1dzic * da2dyid + d2eda1a2 * da2dzic * da1dyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedang2 * dxiczid2 
			    + d2eda2a2 * da2dxic * da2dzid + d2eda1a2 * 
			    da1dxic * da2dzid + d2eda1a2 * da2dxic * da1dzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedang2 * dyiczid2 
			    + d2eda2a2 * da2dyic * da2dzid + d2eda1a2 * 
			    da1dyic * da2dzid + d2eda1a2 * da2dyic * da1dzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedang2 * dziczid2 
			    + d2eda2a2 * da2dzic * da2dzid + d2eda1a2 * 
			    da1dzic * da2dzid + d2eda1a2 * da2dzic * da1dzid;
		    hessx_ref(1, ie) = hessx_ref(1, ie) + dedang2 * dxicxie2 
			    + d2eda2a2 * da2dxic * da2dxie + d2eda1a2 * 
			    da1dxic * da2dxie;
		    hessy_ref(1, ie) = hessy_ref(1, ie) + dedang2 * dyicxie2 
			    + d2eda2a2 * da2dyic * da2dxie + d2eda1a2 * 
			    da1dyic * da2dxie;
		    hessz_ref(1, ie) = hessz_ref(1, ie) + dedang2 * dzicxie2 
			    + d2eda2a2 * da2dzic * da2dxie + d2eda1a2 * 
			    da1dzic * da2dxie;
		    hessx_ref(2, ie) = hessx_ref(2, ie) + dedang2 * dxicyie2 
			    + d2eda2a2 * da2dxic * da2dyie + d2eda1a2 * 
			    da1dxic * da2dyie;
		    hessy_ref(2, ie) = hessy_ref(2, ie) + dedang2 * dyicyie2 
			    + d2eda2a2 * da2dyic * da2dyie + d2eda1a2 * 
			    da1dyic * da2dyie;
		    hessz_ref(2, ie) = hessz_ref(2, ie) + dedang2 * dzicyie2 
			    + d2eda2a2 * da2dzic * da2dyie + d2eda1a2 * 
			    da1dzic * da2dyie;
		    hessx_ref(3, ie) = hessx_ref(3, ie) + dedang2 * dxiczie2 
			    + d2eda2a2 * da2dxic * da2dzie + d2eda1a2 * 
			    da1dxic * da2dzie;
		    hessy_ref(3, ie) = hessy_ref(3, ie) + dedang2 * dyiczie2 
			    + d2eda2a2 * da2dyic * da2dzie + d2eda1a2 * 
			    da1dyic * da2dzie;
		    hessz_ref(3, ie) = hessz_ref(3, ie) + dedang2 * dziczie2 
			    + d2eda2a2 * da2dzic * da2dzie + d2eda1a2 * 
			    da1dzic * da2dzie;
		} else if (*i__ == id) {
		    hessx_ref(1, id) = hessx_ref(1, id) + dedang2 * dxidxid2 
			    + d2eda2a2 * da2dxid * da2dxid + d2eda1a2 * 
			    da1dxid * da2dxid + d2eda1a2 * da2dxid * da1dxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedang2 * dxidyid2 
			    + d2eda2a2 * da2dxid * da2dyid + d2eda1a2 * 
			    da1dxid * da2dyid + d2eda1a2 * da2dxid * da1dyid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedang2 * dxidzid2 
			    + d2eda2a2 * da2dxid * da2dzid + d2eda1a2 * 
			    da1dxid * da2dzid + d2eda1a2 * da2dxid * da1dzid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedang2 * dxidyid2 
			    + d2eda2a2 * da2dxid * da2dyid + d2eda1a2 * 
			    da1dxid * da2dyid + d2eda1a2 * da2dxid * da1dyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedang2 * dyidyid2 
			    + d2eda2a2 * da2dyid * da2dyid + d2eda1a2 * 
			    da1dyid * da2dyid + d2eda1a2 * da2dyid * da1dyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedang2 * dyidzid2 
			    + d2eda2a2 * da2dyid * da2dzid + d2eda1a2 * 
			    da1dyid * da2dzid + d2eda1a2 * da2dyid * da1dzid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedang2 * dxidzid2 
			    + d2eda2a2 * da2dxid * da2dzid + d2eda1a2 * 
			    da1dxid * da2dzid + d2eda1a2 * da2dxid * da1dzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedang2 * dyidzid2 
			    + d2eda2a2 * da2dyid * da2dzid + d2eda1a2 * 
			    da1dyid * da2dzid + d2eda1a2 * da2dyid * da1dzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedang2 * dzidzid2 
			    + d2eda2a2 * da2dzid * da2dzid + d2eda1a2 * 
			    da1dzid * da2dzid + d2eda1a2 * da2dzid * da1dzid;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + d2eda1a2 * da2dxid *
			     da1dxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + d2eda1a2 * da2dyid *
			     da1dxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + d2eda1a2 * da2dzid *
			     da1dxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + d2eda1a2 * da2dxid *
			     da1dyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + d2eda1a2 * da2dyid *
			     da1dyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + d2eda1a2 * da2dzid *
			     da1dyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + d2eda1a2 * da2dxid *
			     da1dzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + d2eda1a2 * da2dyid *
			     da1dzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + d2eda1a2 * da2dzid *
			     da1dzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedang2 * dxibxid2 
			    + d2eda2a2 * da2dxid * da2dxib + d2eda1a2 * 
			    da1dxid * da2dxib + d2eda1a2 * da2dxid * da1dxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedang2 * dxibyid2 
			    + d2eda2a2 * da2dyid * da2dxib + d2eda1a2 * 
			    da1dyid * da2dxib + d2eda1a2 * da2dyid * da1dxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedang2 * dxibzid2 
			    + d2eda2a2 * da2dzid * da2dxib + d2eda1a2 * 
			    da1dzid * da2dxib + d2eda1a2 * da2dzid * da1dxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedang2 * dyibxid2 
			    + d2eda2a2 * da2dxid * da2dyib + d2eda1a2 * 
			    da1dxid * da2dyib + d2eda1a2 * da2dxid * da1dyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedang2 * dyibyid2 
			    + d2eda2a2 * da2dyid * da2dyib + d2eda1a2 * 
			    da1dyid * da2dyib + d2eda1a2 * da2dyid * da1dyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedang2 * dyibzid2 
			    + d2eda2a2 * da2dzid * da2dyib + d2eda1a2 * 
			    da1dzid * da2dyib + d2eda1a2 * da2dzid * da1dyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedang2 * dzibxid2 
			    + d2eda2a2 * da2dxid * da2dzib + d2eda1a2 * 
			    da1dxid * da2dzib + d2eda1a2 * da2dxid * da1dzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedang2 * dzibyid2 
			    + d2eda2a2 * da2dyid * da2dzib + d2eda1a2 * 
			    da1dyid * da2dzib + d2eda1a2 * da2dyid * da1dzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedang2 * dzibzid2 
			    + d2eda2a2 * da2dzid * da2dzib + d2eda1a2 * 
			    da1dzid * da2dzib + d2eda1a2 * da2dzid * da1dzib;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedang2 * dxicxid2 
			    + d2eda2a2 * da2dxid * da2dxic + d2eda1a2 * 
			    da1dxid * da2dxic + d2eda1a2 * da2dxid * da1dxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedang2 * dxicyid2 
			    + d2eda2a2 * da2dyid * da2dxic + d2eda1a2 * 
			    da1dyid * da2dxic + d2eda1a2 * da2dyid * da1dxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedang2 * dxiczid2 
			    + d2eda2a2 * da2dzid * da2dxic + d2eda1a2 * 
			    da1dzid * da2dxic + d2eda1a2 * da2dzid * da1dxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedang2 * dyicxid2 
			    + d2eda2a2 * da2dxid * da2dyic + d2eda1a2 * 
			    da1dxid * da2dyic + d2eda1a2 * da2dxid * da1dyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedang2 * dyicyid2 
			    + d2eda2a2 * da2dyid * da2dyic + d2eda1a2 * 
			    da1dyid * da2dyic + d2eda1a2 * da2dyid * da1dyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedang2 * dyiczid2 
			    + d2eda2a2 * da2dzid * da2dyic + d2eda1a2 * 
			    da1dzid * da2dyic + d2eda1a2 * da2dzid * da1dyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedang2 * dzicxid2 
			    + d2eda2a2 * da2dxid * da2dzic + d2eda1a2 * 
			    da1dxid * da2dzic + d2eda1a2 * da2dxid * da1dzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedang2 * dzicyid2 
			    + d2eda2a2 * da2dyid * da2dzic + d2eda1a2 * 
			    da1dyid * da2dzic + d2eda1a2 * da2dyid * da1dzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedang2 * dziczid2 
			    + d2eda2a2 * da2dzid * da2dzic + d2eda1a2 * 
			    da1dzid * da2dzic + d2eda1a2 * da2dzid * da1dzic;
		    hessx_ref(1, ie) = hessx_ref(1, ie) + dedang2 * dxidxie2 
			    + d2eda2a2 * da2dxid * da2dxie + d2eda1a2 * 
			    da1dxid * da2dxie;
		    hessy_ref(1, ie) = hessy_ref(1, ie) + dedang2 * dyidxie2 
			    + d2eda2a2 * da2dyid * da2dxie + d2eda1a2 * 
			    da1dyid * da2dxie;
		    hessz_ref(1, ie) = hessz_ref(1, ie) + dedang2 * dzidxie2 
			    + d2eda2a2 * da2dzid * da2dxie + d2eda1a2 * 
			    da1dzid * da2dxie;
		    hessx_ref(2, ie) = hessx_ref(2, ie) + dedang2 * dxidyie2 
			    + d2eda2a2 * da2dxid * da2dyie + d2eda1a2 * 
			    da1dxid * da2dyie;
		    hessy_ref(2, ie) = hessy_ref(2, ie) + dedang2 * dyidyie2 
			    + d2eda2a2 * da2dyid * da2dyie + d2eda1a2 * 
			    da1dyid * da2dyie;
		    hessz_ref(2, ie) = hessz_ref(2, ie) + dedang2 * dzidyie2 
			    + d2eda2a2 * da2dzid * da2dyie + d2eda1a2 * 
			    da1dzid * da2dyie;
		    hessx_ref(3, ie) = hessx_ref(3, ie) + dedang2 * dxidzie2 
			    + d2eda2a2 * da2dxid * da2dzie + d2eda1a2 * 
			    da1dxid * da2dzie;
		    hessy_ref(3, ie) = hessy_ref(3, ie) + dedang2 * dyidzie2 
			    + d2eda2a2 * da2dyid * da2dzie + d2eda1a2 * 
			    da1dyid * da2dzie;
		    hessz_ref(3, ie) = hessz_ref(3, ie) + dedang2 * dzidzie2 
			    + d2eda2a2 * da2dzid * da2dzie + d2eda1a2 * 
			    da1dzid * da2dzie;
		} else if (*i__ == ie) {
		    hessx_ref(1, ie) = hessx_ref(1, ie) + dedang2 * dxiexie2 
			    + d2eda2a2 * da2dxie * da2dxie;
		    hessy_ref(1, ie) = hessy_ref(1, ie) + dedang2 * dxieyie2 
			    + d2eda2a2 * da2dxie * da2dyie;
		    hessz_ref(1, ie) = hessz_ref(1, ie) + dedang2 * dxiezie2 
			    + d2eda2a2 * da2dxie * da2dzie;
		    hessx_ref(2, ie) = hessx_ref(2, ie) + dedang2 * dxieyie2 
			    + d2eda2a2 * da2dxie * da2dyie;
		    hessy_ref(2, ie) = hessy_ref(2, ie) + dedang2 * dyieyie2 
			    + d2eda2a2 * da2dyie * da2dyie;
		    hessz_ref(2, ie) = hessz_ref(2, ie) + dedang2 * dyiezie2 
			    + d2eda2a2 * da2dyie * da2dzie;
		    hessx_ref(3, ie) = hessx_ref(3, ie) + dedang2 * dxiezie2 
			    + d2eda2a2 * da2dxie * da2dzie;
		    hessy_ref(3, ie) = hessy_ref(3, ie) + dedang2 * dyiezie2 
			    + d2eda2a2 * da2dyie * da2dzie;
		    hessz_ref(3, ie) = hessz_ref(3, ie) + dedang2 * dziezie2 
			    + d2eda2a2 * da2dzie * da2dzie;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + d2eda1a2 * da2dxie *
			     da1dxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + d2eda1a2 * da2dyie *
			     da1dxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + d2eda1a2 * da2dzie *
			     da1dxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + d2eda1a2 * da2dxie *
			     da1dyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + d2eda1a2 * da2dyie *
			     da1dyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + d2eda1a2 * da2dzie *
			     da1dyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + d2eda1a2 * da2dxie *
			     da1dzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + d2eda1a2 * da2dyie *
			     da1dzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + d2eda1a2 * da2dzie *
			     da1dzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedang2 * dxibxie2 
			    + d2eda2a2 * da2dxie * da2dxib + d2eda1a2 * 
			    da2dxie * da1dxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedang2 * dxibyie2 
			    + d2eda2a2 * da2dyie * da2dxib + d2eda1a2 * 
			    da2dyie * da1dxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedang2 * dxibzie2 
			    + d2eda2a2 * da2dzie * da2dxib + d2eda1a2 * 
			    da2dzie * da1dxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedang2 * dyibxie2 
			    + d2eda2a2 * da2dxie * da2dyib + d2eda1a2 * 
			    da2dxie * da1dyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedang2 * dyibyie2 
			    + d2eda2a2 * da2dyie * da2dyib + d2eda1a2 * 
			    da2dyie * da1dyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedang2 * dyibzie2 
			    + d2eda2a2 * da2dzie * da2dyib + d2eda1a2 * 
			    da2dzie * da1dyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedang2 * dzibxie2 
			    + d2eda2a2 * da2dxie * da2dzib + d2eda1a2 * 
			    da2dxie * da1dzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedang2 * dzibyie2 
			    + d2eda2a2 * da2dyie * da2dzib + d2eda1a2 * 
			    da2dyie * da1dzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedang2 * dzibzie2 
			    + d2eda2a2 * da2dzie * da2dzib + d2eda1a2 * 
			    da2dzie * da1dzib;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedang2 * dxicxie2 
			    + d2eda2a2 * da2dxie * da2dxic + d2eda1a2 * 
			    da2dxie * da1dxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedang2 * dxicyie2 
			    + d2eda2a2 * da2dyie * da2dxic + d2eda1a2 * 
			    da2dyie * da1dxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedang2 * dxiczie2 
			    + d2eda2a2 * da2dzie * da2dxic + d2eda1a2 * 
			    da2dzie * da1dxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedang2 * dyicxie2 
			    + d2eda2a2 * da2dxie * da2dyic + d2eda1a2 * 
			    da2dxie * da1dyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedang2 * dyicyie2 
			    + d2eda2a2 * da2dyie * da2dyic + d2eda1a2 * 
			    da2dyie * da1dyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedang2 * dyiczie2 
			    + d2eda2a2 * da2dzie * da2dyic + d2eda1a2 * 
			    da2dzie * da1dyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedang2 * dzicxie2 
			    + d2eda2a2 * da2dxie * da2dzic + d2eda1a2 * 
			    da2dxie * da1dzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedang2 * dzicyie2 
			    + d2eda2a2 * da2dyie * da2dzic + d2eda1a2 * 
			    da2dyie * da1dzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedang2 * dziczie2 
			    + d2eda2a2 * da2dzie * da2dzic + d2eda1a2 * 
			    da2dzie * da1dzic;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedang2 * dxidxie2 
			    + d2eda2a2 * da2dxid * da2dxie + d2eda1a2 * 
			    da1dxid * da2dxie;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedang2 * dxidyie2 
			    + d2eda2a2 * da2dxid * da2dyie + d2eda1a2 * 
			    da1dxid * da2dyie;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedang2 * dxidzie2 
			    + d2eda2a2 * da2dxid * da2dzie + d2eda1a2 * 
			    da1dxid * da2dzie;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedang2 * dyidxie2 
			    + d2eda2a2 * da2dyid * da2dxie + d2eda1a2 * 
			    da1dyid * da2dxie;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedang2 * dyidyie2 
			    + d2eda2a2 * da2dyid * da2dyie + d2eda1a2 * 
			    da1dyid * da2dyie;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedang2 * dyidzie2 
			    + d2eda2a2 * da2dyid * da2dzie + d2eda1a2 * 
			    da1dyid * da2dzie;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedang2 * dzidxie2 
			    + d2eda2a2 * da2dzid * da2dxie + d2eda1a2 * 
			    da1dzid * da2dxie;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedang2 * dzidyie2 
			    + d2eda2a2 * da2dzid * da2dyie + d2eda1a2 * 
			    da1dzid * da2dyie;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedang2 * dzidzie2 
			    + d2eda2a2 * da2dzid * da2dzie + d2eda1a2 * 
			    da1dzid * da2dzie;
		}
	    }
	}
    }
    return 0;
} /* etortor2_ */

#undef ibitor_ref
#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef tbxy_ref
#undef tty_ref
#undef ttx_ref
#undef itt_ref
#undef tby_ref
#undef tbx_ref
#undef tbf_ref


