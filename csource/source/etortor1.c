/* etortor1.f -- translated by f2c (version 20050501).
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
    doublereal desum[75000]	/* was [3][25000] */, deb[75000]	/* 
	    was [3][25000] */, dea[75000]	/* was [3][25000] */, deba[
	    75000]	/* was [3][25000] */, deub[75000]	/* was [3][
	    25000] */, deaa[75000]	/* was [3][25000] */, deopb[75000]	
	    /* was [3][25000] */, deopd[75000]	/* was [3][25000] */, deid[
	    75000]	/* was [3][25000] */, deit[75000]	/* was [3][
	    25000] */, det[75000]	/* was [3][25000] */, dept[75000]	
	    /* was [3][25000] */, debt[75000]	/* was [3][25000] */, dett[
	    75000]	/* was [3][25000] */, dev[75000]	/* was [3][
	    25000] */, dec[75000]	/* was [3][25000] */, decd[75000]	
	    /* was [3][25000] */, ded[75000]	/* was [3][25000] */, dem[
	    75000]	/* was [3][25000] */, dep[75000]	/* was [3][
	    25000] */, der[75000]	/* was [3][25000] */, des[75000]	
	    /* was [3][25000] */, delf[75000]	/* was [3][25000] */, deg[
	    75000]	/* was [3][25000] */, dex[75000]	/* was [3][
	    25000] */;
} deriv_;

#define deriv_1 deriv_

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

struct {
    doublereal vir[9]	/* was [3][3] */;
} virial_;

#define virial_1 virial_

/* Table of constant values */

static integer c__0 = 0;



/*     ############################################################# */
/*     ##  COPYRIGHT (C) 2003 by Pengyu Ren & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine etortor1  --  torsion-torsion energy & derivs  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "etortor1" calculates the torsion-torsion energy and first */
/*     derivatives with respect to Cartesian coordinates */


/* Subroutine */ int etortor1_(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal);

    /* Local variables */
    static doublereal e;
    static integer i__, k, ia, ib, ic, id, ie;
    static doublereal xh, yh;
    static integer nt;
    static doublereal xt, yt, zt, xu, yu, zu, xv, yv, zv, ft1[4], ft2[4], x1l,
	     y1l, rt2, ru2, rv2, x1u, y1u, rcb, rdc;
    static integer nhi;
    static doublereal xia, yia, zia, xib;
    static integer nlo;
    static doublereal yib, zib;
    static integer xlo, ylo;
    static doublereal xic, yic, zic, xid, yid, zid, xie, yie, xtu, ytu, ztu, 
	    xuv, yuv, zuv, zie, xba, yba, zba, xdc, ydc, zdc, xcb, ycb, zcb, 
	    xed, yed, zed, xca, yca, zca, xdb, ydb, zdb, xec, yec, zec, vxx, 
	    vyy, vzz, vyx;
    static integer pos1, pos2;
    static doublereal vzx, vzy, ftt[4], ft12[4], vxx2, vyx2, vyy2, vzz2, vzx2,
	     vzy2, fgrp, sign, rtru, rurv;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal dedxt, dedyt, dedzt, dedxu, dedyu, dedzu, angle1, 
	    angle2, dedxu2, dedyu2, value1, value2, dedzu2, dedxv2, dedyv2, 
	    dedzv2, dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, dedxic, 
	    dedyic, dedzic, dedxid, dedyid, dedzid, dedang1, dedang2;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal dedxib2, dedyib2, dedzib2, dedxic2, dedyic2, dedzic2, 
	    dedxid2, dedyid2, dedzid2, dedxie2, dedyie2, dedzie2, cosine1, 
	    cosine2;
    extern /* Subroutine */ int bcuint1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static logical proceed;
    extern /* Subroutine */ int chkttor_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static integer itortor;


#define tbf_ref(a_1,a_2) ktrtor_1.tbf[(a_2)*900 + a_1 - 901]
#define tbx_ref(a_1,a_2) ktrtor_1.tbx[(a_2)*900 + a_1 - 901]
#define tby_ref(a_1,a_2) ktrtor_1.tby[(a_2)*900 + a_1 - 901]
#define itt_ref(a_1,a_2) tortor_1.itt[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define ttx_ref(a_1,a_2) ktrtor_1.ttx[(a_2)*30 + a_1 - 31]
#define tty_ref(a_1,a_2) ktrtor_1.tty[(a_2)*30 + a_1 - 31]
#define dett_ref(a_1,a_2) deriv_1.dett[(a_2)*3 + a_1 - 4]
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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  deriv.i  --  Cartesian coordinate derivative components  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     desum   total energy Cartesian coordinate derivatives */
/*     deb     bond stretch Cartesian coordinate derivatives */
/*     dea     angle bend Cartesian coordinate derivatives */
/*     deba    stretch-bend Cartesian coordinate derivatives */
/*     deub    Urey-Bradley Cartesian coordinate derivatives */
/*     deaa    angle-angle Cartesian coordinate derivatives */
/*     deopb   out-of-plane bend Cartesian coordinate derivatives */
/*     deopd   out-of-plane distance Cartesian coordinate derivatives */
/*     deid    improper dihedral Cartesian coordinate derivatives */
/*     deit    improper torsion Cartesian coordinate derivatives */
/*     det     torsional Cartesian coordinate derivatives */
/*     dept    pi-orbital torsion Cartesian coordinate derivatives */
/*     debt    stretch-torsion Cartesian coordinate derivatives */
/*     dett    torsion-torsion Cartesian coordinate derivatives */
/*     dev     van der Waals Cartesian coordinate derivatives */
/*     dec     charge-charge Cartesian coordinate derivatives */
/*     decd    charge-dipole Cartesian coordinate derivatives */
/*     ded     dipole-dipole Cartesian coordinate derivatives */
/*     dem     multipole Cartesian coordinate derivatives */
/*     dep     polarization Cartesian coordinate derivatives */
/*     der     reaction field Cartesian coordinate derivatives */
/*     des     solvation Cartesian coordinate derivatives */
/*     delf    metal ligand field Cartesian coordinate derivatives */
/*     deg     geometric restraint Cartesian coordinate derivatives */
/*     dex     extra energy term Cartesian coordinate derivatives */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  virial.i  --  components of internal virial tensor  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     vir    total internal virial Cartesian tensor components */




/*     zero out the torsion-torsion energy and first derivatives */

    energi_1.ett = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dett_ref(1, i__) = 0.;
	dett_ref(2, i__) = 0.;
	dett_ref(3, i__) = 0.;
    }

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
		bcuint1_(ftt, ft1, ft2, ft12, &x1l, &x1u, &y1l, &y1u, &value1,
			 &value2, &e, &dedang1, &dedang2);
		e = torpot_1.ttorunit * e;
		dedang1 = sign * torpot_1.ttorunit * 57.29577951308232088 * 
			dedang1;
		dedang2 = sign * torpot_1.ttorunit * 57.29577951308232088 * 
			dedang2;

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    e *= fgrp;
		    dedang1 *= fgrp;
		    dedang2 *= fgrp;
		}

/*     chain rule terms for first angle derivative components */

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
		dedxt = dedang1 * (yt * zcb - ycb * zt) / (rt2 * rcb);
		dedyt = dedang1 * (zt * xcb - zcb * xt) / (rt2 * rcb);
		dedzt = dedang1 * (xt * ycb - xcb * yt) / (rt2 * rcb);
		dedxu = -dedang1 * (yu * zcb - ycb * zu) / (ru2 * rcb);
		dedyu = -dedang1 * (zu * xcb - zcb * xu) / (ru2 * rcb);
		dedzu = -dedang1 * (xu * ycb - xcb * yu) / (ru2 * rcb);

/*     compute first derivative components for first angle */

		dedxia = zcb * dedyt - ycb * dedzt;
		dedyia = xcb * dedzt - zcb * dedxt;
		dedzia = ycb * dedxt - xcb * dedyt;
		dedxib = yca * dedzt - zca * dedyt + zdc * dedyu - ydc * 
			dedzu;
		dedyib = zca * dedxt - xca * dedzt + xdc * dedzu - zdc * 
			dedxu;
		dedzib = xca * dedyt - yca * dedxt + ydc * dedxu - xdc * 
			dedyu;
		dedxic = zba * dedyt - yba * dedzt + ydb * dedzu - zdb * 
			dedyu;
		dedyic = xba * dedzt - zba * dedxt + zdb * dedxu - xdb * 
			dedzu;
		dedzic = yba * dedxt - xba * dedyt + xdb * dedyu - ydb * 
			dedxu;
		dedxid = zcb * dedyu - ycb * dedzu;
		dedyid = xcb * dedzu - zcb * dedxu;
		dedzid = ycb * dedxu - xcb * dedyu;

/*     chain rule terms for second angle derivative components */

		xec = xie - xic;
		yec = yie - yic;
		zec = zie - zic;
		if (bound_1.use_polymer__) {
		    image_(&xdb, &ydb, &zdb);
		    image_(&xec, &yec, &zec);
		}
		dedxu2 = dedang2 * (yu * zdc - ydc * zu) / (ru2 * rdc);
		dedyu2 = dedang2 * (zu * xdc - zdc * xu) / (ru2 * rdc);
		dedzu2 = dedang2 * (xu * ydc - xdc * yu) / (ru2 * rdc);
		dedxv2 = -dedang2 * (yv * zdc - ydc * zv) / (rv2 * rdc);
		dedyv2 = -dedang2 * (zv * xdc - zdc * xv) / (rv2 * rdc);
		dedzv2 = -dedang2 * (xv * ydc - xdc * yv) / (rv2 * rdc);

/*     compute first derivative components for second angle */

		dedxib2 = zdc * dedyu2 - ydc * dedzu2;
		dedyib2 = xdc * dedzu2 - zdc * dedxu2;
		dedzib2 = ydc * dedxu2 - xdc * dedyu2;
		dedxic2 = ydb * dedzu2 - zdb * dedyu2 + zed * dedyv2 - yed * 
			dedzv2;
		dedyic2 = zdb * dedxu2 - xdb * dedzu2 + xed * dedzv2 - zed * 
			dedxv2;
		dedzic2 = xdb * dedyu2 - ydb * dedxu2 + yed * dedxv2 - xed * 
			dedyv2;
		dedxid2 = zcb * dedyu2 - ycb * dedzu2 + yec * dedzv2 - zec * 
			dedyv2;
		dedyid2 = xcb * dedzu2 - zcb * dedxu2 + zec * dedxv2 - xec * 
			dedzv2;
		dedzid2 = ycb * dedxu2 - xcb * dedyu2 + xec * dedyv2 - yec * 
			dedxv2;
		dedxie2 = zdc * dedyv2 - ydc * dedzv2;
		dedyie2 = xdc * dedzv2 - zdc * dedxv2;
		dedzie2 = ydc * dedxv2 - xdc * dedyv2;

/*     increment the torsion-torsion energy and gradient */

		energi_1.ett += e;
		dett_ref(1, ia) = dett_ref(1, ia) + dedxia;
		dett_ref(2, ia) = dett_ref(2, ia) + dedyia;
		dett_ref(3, ia) = dett_ref(3, ia) + dedzia;
		dett_ref(1, ib) = dett_ref(1, ib) + dedxib + dedxib2;
		dett_ref(2, ib) = dett_ref(2, ib) + dedyib + dedyib2;
		dett_ref(3, ib) = dett_ref(3, ib) + dedzib + dedzib2;
		dett_ref(1, ic) = dett_ref(1, ic) + dedxic + dedxic2;
		dett_ref(2, ic) = dett_ref(2, ic) + dedyic + dedyic2;
		dett_ref(3, ic) = dett_ref(3, ic) + dedzic + dedzic2;
		dett_ref(1, id) = dett_ref(1, id) + dedxid + dedxid2;
		dett_ref(2, id) = dett_ref(2, id) + dedyid + dedyid2;
		dett_ref(3, id) = dett_ref(3, id) + dedzid + dedzid2;
		dett_ref(1, ie) = dett_ref(1, ie) + dedxie2;
		dett_ref(2, ie) = dett_ref(2, ie) + dedyie2;
		dett_ref(3, ie) = dett_ref(3, ie) + dedzie2;

/*     increment the internal virial tensor components */

		vxx = xcb * (dedxic + dedxid) - xba * dedxia + xdc * dedxid;
		vyx = ycb * (dedxic + dedxid) - yba * dedxia + ydc * dedxid;
		vzx = zcb * (dedxic + dedxid) - zba * dedxia + zdc * dedxid;
		vyy = ycb * (dedyic + dedyid) - yba * dedyia + ydc * dedyid;
		vzy = zcb * (dedyic + dedyid) - zba * dedyia + zdc * dedyid;
		vzz = zcb * (dedzic + dedzid) - zba * dedzia + zdc * dedzid;
		vxx2 = xdc * (dedxid2 + dedxie2) - xcb * dedxib2 + xed * 
			dedxie2;
		vyx2 = ydc * (dedxid2 + dedxie2) - ycb * dedxib2 + yed * 
			dedxie2;
		vzx2 = zdc * (dedxid2 + dedxie2) - zcb * dedxib2 + zed * 
			dedxie2;
		vyy2 = ydc * (dedyid2 + dedyie2) - ycb * dedyib2 + yed * 
			dedyie2;
		vzy2 = zdc * (dedyid2 + dedyie2) - zcb * dedyib2 + zed * 
			dedyie2;
		vzz2 = zdc * (dedzid2 + dedzie2) - zcb * dedzib2 + zed * 
			dedzie2;
		vir_ref(1, 1) = vir_ref(1, 1) + vxx + vxx2;
		vir_ref(2, 1) = vir_ref(2, 1) + vyx + vyx2;
		vir_ref(3, 1) = vir_ref(3, 1) + vzx + vzx2;
		vir_ref(1, 2) = vir_ref(1, 2) + vyx + vyx2;
		vir_ref(2, 2) = vir_ref(2, 2) + vyy + vyy2;
		vir_ref(3, 2) = vir_ref(3, 2) + vzy + vzy2;
		vir_ref(1, 3) = vir_ref(1, 3) + vzx + vzx2;
		vir_ref(2, 3) = vir_ref(2, 3) + vzy + vzy2;
		vir_ref(3, 3) = vir_ref(3, 3) + vzz + vzz2;
	    }
	}
    }
    return 0;
} /* etortor1_ */

#undef ibitor_ref
#undef tbxy_ref
#undef dett_ref
#undef tty_ref
#undef ttx_ref
#undef vir_ref
#undef itt_ref
#undef tby_ref
#undef tbx_ref
#undef tbf_ref


