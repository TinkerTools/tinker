/* etors2.f -- translated by f2c (version 20050501).
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
    doublereal m2[25000], deform, difft, diffv, diffc;
    logical use_smooth__, use_dem__, use_gda__, use_tophat__, use_stophat__;
} warp_;

#define warp_1 warp_

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
    doublereal idihunit, itorunit, torsunit, ptorunit, storunit, ttorunit;
} torpot_;

#define torpot_1 torpot_

struct {
    doublereal tors1[400000]	/* was [4][100000] */, tors2[400000]	/* 
	    was [4][100000] */, tors3[400000]	/* was [4][100000] */, tors4[
	    400000]	/* was [4][100000] */, tors5[400000]	/* was [4][
	    100000] */, tors6[400000]	/* was [4][100000] */;
    integer ntors, itors[400000]	/* was [4][100000] */;
} tors_;

#define tors_1 tors_

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine etors2  --  atom-by-atom torsional Hessian  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "etors2" calculates the second derivatives of the torsional */
/*     energy for a single atom */


/* Subroutine */ int etors2_(integer *i__)
{
    extern /* Subroutine */ int etors2a_(integer *), etors2b_(integer *);



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
/*     ##  warp.i  --  parameters for potential surface smoothing  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     m2           second moment of the GDA gaussian for each atom */
/*     deform       value of the smoothing deformation parameter */
/*     difft        diffusion coefficient for torsional potential */
/*     diffv        diffusion coefficient for van der Waals potential */
/*     diffc        diffusion coefficient for charge-charge potential */
/*     use_smooth   flag to use a potential energy smoothing method */
/*     use_dem      flag to use diffusion equation method potential */
/*     use_gda      flag to use gaussian density annealing potential */
/*     use_tophat   flag to use analytical tophat smoothed potential */
/*     use_stophat  flag to use shifted tophat smoothed potential */




/*     choose standard or potential energy smoothing version */

    if (warp_1.use_smooth__) {
	etors2b_(i__);
    } else {
	etors2a_(i__);
    }
    return 0;
} /* etors2_ */



/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  subroutine etors2a  --  standard torsional Hessian  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     "etors2a" calculates the second derivatives of the torsional */
/*     energy for a single atom using a standard sum of Fourier terms */


/* Subroutine */ int etors2a_(integer *i__)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal dxiazid, dyiaxid, dyiayid, dyiazid, dziaxid, dziayid, 
	    dziazid, dxibxic, dxibyic, dxibzic, dyibxic, dyibyic, dyibzic, 
	    dzibxic, dzibyic, dzibzic, dxibxid, dxibyid, dxibzid, dyibxid, 
	    dyibyid, dyibzid, dzibxid, dzibyid, dzibzid, dxicxid, dxicyid, 
	    dxiczid, dyicxid, dyicyid, dyiczid, dzicxid, dzicyid, dziczid, 
	    dphidxia, dphidyia, dphidzia, dphidxib, dphidyib, dphidzib, 
	    dphidxic, dphidyic, dphidzic, dphidxid, dphidyid, dphidzid, c1, 
	    c2, c3, c4, c5, c6, dphidxibt, dphidyibt, dphidzibt, dphidxibu, 
	    dphidyibu, dphidzibu, s1, s2, s3, v1, v2, v3, v4, v5, v6, s4, s5, 
	    s6, dphidxict, dphidyict, dphidzict, dphidxicu, dphidyicu, 
	    dphidzicu;
    static integer ia, ib, ic, id;
    static doublereal xt, yt, zt, xu, yu, zu, rt2, ru2, rcb, xba, yba, zba, 
	    xcb, ycb, zcb, xdc, xia, yia, zia, xib, yib, zib, xic, yic, zic, 
	    xid, yid, zid, ydc, zdc, xca, yca, zca, xdb, ydb, zdb, xtu, ytu, 
	    ztu, fgrp, sine, rtru, dphi1, dphi2, dphi3, dphi4, dphi5, dphi6, 
	    rcbt2, rcbu2, sine2, sine3, sine4, sine5, sine6, xycb2, xzcb2, 
	    yzcb2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal d2phi1, d2phi2, d2phi3, d2phi4, d2phi5, d2phi6, rcbxt, 
	    rcbyt, rcbzt, rcbxu, rcbyu, rcbzu;
    static integer ktors;
    static doublereal dedphi, cosine;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal cosine2, cosine3, cosine4, cosine5, cosine6, d2edphi2;
    static logical proceed;
    static doublereal dxiaxia, dxiayia, dyiayia, dxibxib, dziazia, dyibyib, 
	    dzibzib, dxicxic, dyicyic, dziczic, dxidxid, dyidyid, dzidzid, 
	    dphidxt, dphidyt, dphidzt, dphidxu, dphidyu, dphidzu, dxiazia, 
	    dyiazia, dxibyib, dxibzib, dyibzib, dxicyic, dxiczic, dyiczic, 
	    dxidyid, dxidzid, dyidzid, dxiaxib, dxiayib, dxiazib, dyiaxib, 
	    dyiayib, dyiazib, dziaxib, dziayib, dziazib, dxiaxic, dxiayic, 
	    dxiazic, dyiaxic, dyiayic, dyiazic, dziaxic, dziayic, dziazic, 
	    dxiaxid, dxiayid;


#define tors1_ref(a_1,a_2) tors_1.tors1[(a_2)*4 + a_1 - 5]
#define tors2_ref(a_1,a_2) tors_1.tors2[(a_2)*4 + a_1 - 5]
#define tors3_ref(a_1,a_2) tors_1.tors3[(a_2)*4 + a_1 - 5]
#define tors4_ref(a_1,a_2) tors_1.tors4[(a_2)*4 + a_1 - 5]
#define tors5_ref(a_1,a_2) tors_1.tors5[(a_2)*4 + a_1 - 5]
#define tors6_ref(a_1,a_2) tors_1.tors6[(a_2)*4 + a_1 - 5]
#define hessx_ref(a_1,a_2) hessn_1.hessx[(a_2)*3 + a_1 - 4]
#define hessy_ref(a_1,a_2) hessn_1.hessy[(a_2)*3 + a_1 - 4]
#define hessz_ref(a_1,a_2) hessn_1.hessz[(a_2)*3 + a_1 - 4]
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




/*     calculate the torsional angle energy term */

    i__1 = tors_1.ntors;
    for (ktors = 1; ktors <= i__1; ++ktors) {
	ia = itors_ref(1, ktors);
	ib = itors_ref(2, ktors);
	ic = itors_ref(3, ktors);
	id = itors_ref(4, ktors);

/*     decide whether to compute the current interaction */

	proceed = *i__ == ia || *i__ == ib || *i__ == ic || *i__ == id;
	if (proceed && group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &c__0, &c__0);
	}

/*     compute the value of the torsional angle */

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
	    if (bound_1.use_polymer__) {
		image_(&xba, &yba, &zba);
		image_(&xcb, &ycb, &zcb);
		image_(&xdc, &ydc, &zdc);
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
	    if (rtru != 0.) {
		rcb = sqrt(xcb * xcb + ycb * ycb + zcb * zcb);
		cosine = (xt * xu + yt * yu + zt * zu) / rtru;
		sine = (xcb * xtu + ycb * ytu + zcb * ztu) / (rcb * rtru);

/*     set the torsional parameters for this angle */

		v1 = tors1_ref(1, ktors);
		c1 = tors1_ref(3, ktors);
		s1 = tors1_ref(4, ktors);
		v2 = tors2_ref(1, ktors);
		c2 = tors2_ref(3, ktors);
		s2 = tors2_ref(4, ktors);
		v3 = tors3_ref(1, ktors);
		c3 = tors3_ref(3, ktors);
		s3 = tors3_ref(4, ktors);
		v4 = tors4_ref(1, ktors);
		c4 = tors4_ref(3, ktors);
		s4 = tors4_ref(4, ktors);
		v5 = tors5_ref(1, ktors);
		c5 = tors5_ref(3, ktors);
		s5 = tors5_ref(4, ktors);
		v6 = tors6_ref(1, ktors);
		c6 = tors6_ref(3, ktors);
		s6 = tors6_ref(4, ktors);

/*     compute the multiple angle trigonometry and the phase terms */

		cosine2 = cosine * cosine - sine * sine;
		sine2 = cosine * 2. * sine;
		cosine3 = cosine * cosine2 - sine * sine2;
		sine3 = cosine * sine2 + sine * cosine2;
		cosine4 = cosine * cosine3 - sine * sine3;
		sine4 = cosine * sine3 + sine * cosine3;
		cosine5 = cosine * cosine4 - sine * sine4;
		sine5 = cosine * sine4 + sine * cosine4;
		cosine6 = cosine * cosine5 - sine * sine5;
		sine6 = cosine * sine5 + sine * cosine5;
		dphi1 = cosine * s1 - sine * c1;
		dphi2 = (cosine2 * s2 - sine2 * c2) * 2.;
		dphi3 = (cosine3 * s3 - sine3 * c3) * 3.;
		dphi4 = (cosine4 * s4 - sine4 * c4) * 4.;
		dphi5 = (cosine5 * s5 - sine5 * c5) * 5.;
		dphi6 = (cosine6 * s6 - sine6 * c6) * 6.;
		d2phi1 = -(cosine * c1 + sine * s1);
		d2phi2 = (cosine2 * c2 + sine2 * s2) * -4.;
		d2phi3 = (cosine3 * c3 + sine3 * s3) * -9.;
		d2phi4 = (cosine4 * c4 + sine4 * s4) * -16.;
		d2phi5 = (cosine5 * c5 + sine5 * s5) * -25.;
		d2phi6 = (cosine6 * c6 + sine6 * s6) * -36.;

/*     calculate the torsional master chain rule terms */

		dedphi = torpot_1.torsunit * (v1 * dphi1 + v2 * dphi2 + v3 * 
			dphi3 + v4 * dphi4 + v5 * dphi5 + v6 * dphi6);
		d2edphi2 = torpot_1.torsunit * (v1 * d2phi1 + v2 * d2phi2 + 
			v3 * d2phi3 + v4 * d2phi4 + v5 * d2phi5 + v6 * d2phi6)
			;

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    dedphi *= fgrp;
		    d2edphi2 *= fgrp;
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
		dphidxt = (yt * zcb - ycb * zt) / (rt2 * rcb);
		dphidyt = (zt * xcb - zcb * xt) / (rt2 * rcb);
		dphidzt = (xt * ycb - xcb * yt) / (rt2 * rcb);
		dphidxu = -(yu * zcb - ycb * zu) / (ru2 * rcb);
		dphidyu = -(zu * xcb - zcb * xu) / (ru2 * rcb);
		dphidzu = -(xu * ycb - xcb * yu) / (ru2 * rcb);

/*     abbreviations for second derivative chain rule terms */

		xycb2 = xcb * xcb + ycb * ycb;
		xzcb2 = xcb * xcb + zcb * zcb;
		yzcb2 = ycb * ycb + zcb * zcb;
		rcbxt = rcb * -2. * dphidxt;
		rcbyt = rcb * -2. * dphidyt;
		rcbzt = rcb * -2. * dphidzt;
		rcbt2 = rcb * rt2;
		rcbxu = rcb * 2. * dphidxu;
		rcbyu = rcb * 2. * dphidyu;
		rcbzu = rcb * 2. * dphidzu;
		rcbu2 = rcb * ru2;
		dphidxibt = yca * dphidzt - zca * dphidyt;
		dphidxibu = zdc * dphidyu - ydc * dphidzu;
		dphidyibt = zca * dphidxt - xca * dphidzt;
		dphidyibu = xdc * dphidzu - zdc * dphidxu;
		dphidzibt = xca * dphidyt - yca * dphidxt;
		dphidzibu = ydc * dphidxu - xdc * dphidyu;
		dphidxict = zba * dphidyt - yba * dphidzt;
		dphidxicu = ydb * dphidzu - zdb * dphidyu;
		dphidyict = xba * dphidzt - zba * dphidxt;
		dphidyicu = zdb * dphidxu - xdb * dphidzu;
		dphidzict = yba * dphidxt - xba * dphidyt;
		dphidzicu = xdb * dphidyu - ydb * dphidxu;

/*     chain rule terms for first derivative components */

		dphidxia = zcb * dphidyt - ycb * dphidzt;
		dphidyia = xcb * dphidzt - zcb * dphidxt;
		dphidzia = ycb * dphidxt - xcb * dphidyt;
		dphidxib = dphidxibt + dphidxibu;
		dphidyib = dphidyibt + dphidyibu;
		dphidzib = dphidzibt + dphidzibu;
		dphidxic = dphidxict + dphidxicu;
		dphidyic = dphidyict + dphidyicu;
		dphidzic = dphidzict + dphidzicu;
		dphidxid = zcb * dphidyu - ycb * dphidzu;
		dphidyid = xcb * dphidzu - zcb * dphidxu;
		dphidzid = ycb * dphidxu - xcb * dphidyu;

/*     chain rule terms for second derivative components */

		dxiaxia = rcbxt * dphidxia;
		dxiayia = rcbxt * dphidyia - zcb * rcb / rt2;
		dxiazia = rcbxt * dphidzia + ycb * rcb / rt2;
		dxiaxic = rcbxt * dphidxict + xcb * xt / rcbt2;
		dxiayic = rcbxt * dphidyict - dphidzt - (xba * zcb * xcb + 
			zba * yzcb2) / rcbt2;
		dxiazic = rcbxt * dphidzict + dphidyt + (xba * ycb * xcb + 
			yba * yzcb2) / rcbt2;
		dxiaxid = 0.;
		dxiayid = 0.;
		dxiazid = 0.;
		dyiayia = rcbyt * dphidyia;
		dyiazia = rcbyt * dphidzia - xcb * rcb / rt2;
		dyiaxib = rcbyt * dphidxibt - dphidzt - (yca * zcb * ycb + 
			zca * xzcb2) / rcbt2;
		dyiaxic = rcbyt * dphidxict + dphidzt + (yba * zcb * ycb + 
			zba * xzcb2) / rcbt2;
		dyiayic = rcbyt * dphidyict + ycb * yt / rcbt2;
		dyiazic = rcbyt * dphidzict - dphidxt - (yba * xcb * ycb + 
			xba * xzcb2) / rcbt2;
		dyiaxid = 0.;
		dyiayid = 0.;
		dyiazid = 0.;
		dziazia = rcbzt * dphidzia;
		dziaxib = rcbzt * dphidxibt + dphidyt + (zca * ycb * zcb + 
			yca * xycb2) / rcbt2;
		dziayib = rcbzt * dphidyibt - dphidxt - (zca * xcb * zcb + 
			xca * xycb2) / rcbt2;
		dziaxic = rcbzt * dphidxict - dphidyt - (zba * ycb * zcb + 
			yba * xycb2) / rcbt2;
		dziayic = rcbzt * dphidyict + dphidxt + (zba * xcb * zcb + 
			xba * xycb2) / rcbt2;
		dziazic = rcbzt * dphidzict + zcb * zt / rcbt2;
		dziaxid = 0.;
		dziayid = 0.;
		dziazid = 0.;
		dxibxic = -xcb * dphidxib / (rcb * rcb) - (yca * (zba * xcb + 
			yt) - zca * (yba * xcb - zt)) / rcbt2 - (yt * zba - 
			yba * zt) * 2. * dphidxibt / rt2 - (zdc * (ydb * xcb 
			+ zu) - ydc * (zdb * xcb - yu)) / rcbu2 + (yu * zdb - 
			ydb * zu) * 2. * dphidxibu / ru2;
		dxibyic = -ycb * dphidxib / (rcb * rcb) + dphidzt + dphidzu - 
			(yca * (zba * ycb - xt) + zca * (xba * xcb + zcb * 
			zba)) / rcbt2 - (zt * xba - zba * xt) * 2. * 
			dphidxibt / rt2 + (zdc * (xdb * xcb + zcb * zdb) + 
			ydc * (zdb * ycb + xu)) / rcbu2 + (zu * xdb - zdb * 
			xu) * 2. * dphidxibu / ru2;
		dxibxid = rcbxu * dphidxibu + xcb * xu / rcbu2;
		dxibyid = rcbyu * dphidxibu - dphidzu - (ydc * zcb * ycb + 
			zdc * xzcb2) / rcbu2;
		dxibzid = rcbzu * dphidxibu + dphidyu + (zdc * ycb * zcb + 
			ydc * xycb2) / rcbu2;
		dyibzib = ycb * dphidzib / (rcb * rcb) - (xca * (xca * xcb + 
			zcb * zca) + yca * (ycb * xca + zt)) / rcbt2 - (xt * 
			zca - xca * zt) * 2. * dphidzibt / rt2 + (ydc * (xdc *
			 ycb - zu) + xdc * (xdc * xcb + zcb * zdc)) / rcbu2 + 
			(xu * zdc - xdc * zu) * 2. * dphidzibu / ru2;
		dyibxic = -xcb * dphidyib / (rcb * rcb) - dphidzt - dphidzu + 
			(xca * (zba * xcb + yt) + zca * (zba * zcb + ycb * 
			yba)) / rcbt2 - (yt * zba - yba * zt) * 2. * 
			dphidyibt / rt2 - (zdc * (zdb * zcb + ycb * ydb) + 
			xdc * (zdb * xcb - yu)) / rcbu2 + (yu * zdb - ydb * 
			zu) * 2. * dphidyibu / ru2;
		dyibyic = -ycb * dphidyib / (rcb * rcb) - (zca * (xba * ycb + 
			zt) - xca * (zba * ycb - xt)) / rcbt2 - (zt * xba - 
			zba * xt) * 2. * dphidyibt / rt2 - (xdc * (zdb * ycb 
			+ xu) - zdc * (xdb * ycb - zu)) / rcbu2 + (zu * xdb - 
			zdb * xu) * 2. * dphidyibu / ru2;
		dyibxid = rcbxu * dphidyibu + dphidzu + (xdc * zcb * xcb + 
			zdc * yzcb2) / rcbu2;
		dyibyid = rcbyu * dphidyibu + ycb * yu / rcbu2;
		dyibzid = rcbzu * dphidyibu - dphidxu - (zdc * xcb * zcb + 
			xdc * xycb2) / rcbu2;
		dzibxic = -xcb * dphidzib / (rcb * rcb) + dphidyt + dphidyu - 
			(xca * (yba * xcb - zt) + yca * (zba * zcb + ycb * 
			yba)) / rcbt2 - (yt * zba - yba * zt) * 2. * 
			dphidzibt / rt2 + (ydc * (zdb * zcb + ycb * ydb) + 
			xdc * (ydb * xcb + zu)) / rcbu2 + (yu * zdb - ydb * 
			zu) * 2. * dphidzibu / ru2;
		dzibzic = -zcb * dphidzib / (rcb * rcb) - (xca * (yba * zcb + 
			xt) - yca * (xba * zcb - yt)) / rcbt2 - (xt * yba - 
			xba * yt) * 2. * dphidzibt / rt2 - (ydc * (xdb * zcb 
			+ yu) - xdc * (ydb * zcb - xu)) / rcbu2 + (xu * ydb - 
			xdb * yu) * 2. * dphidzibu / ru2;
		dzibxid = rcbxu * dphidzibu - dphidyu - (xdc * ycb * xcb + 
			ydc * yzcb2) / rcbu2;
		dzibyid = rcbyu * dphidzibu + dphidxu + (ydc * xcb * ycb + 
			xdc * xzcb2) / rcbu2;
		dzibzid = rcbzu * dphidzibu + zcb * zu / rcbu2;
		dxicxid = rcbxu * dphidxicu - xcb * (zdb * ycb - ydb * zcb) / 
			rcbu2;
		dxicyid = rcbyu * dphidxicu + dphidzu + (ydb * zcb * ycb + 
			zdb * xzcb2) / rcbu2;
		dxiczid = rcbzu * dphidxicu - dphidyu - (zdb * ycb * zcb + 
			ydb * xycb2) / rcbu2;
		dyicxid = rcbxu * dphidyicu - dphidzu - (xdb * zcb * xcb + 
			zdb * yzcb2) / rcbu2;
		dyicyid = rcbyu * dphidyicu - ycb * (xdb * zcb - zdb * xcb) / 
			rcbu2;
		dyiczid = rcbzu * dphidyicu + dphidxu + (zdb * xcb * zcb + 
			xdb * xycb2) / rcbu2;
		dzicxid = rcbxu * dphidzicu + dphidyu + (xdb * ycb * xcb + 
			ydb * yzcb2) / rcbu2;
		dzicyid = rcbyu * dphidzicu - dphidxu - (ydb * xcb * ycb + 
			xdb * xzcb2) / rcbu2;
		dziczid = rcbzu * dphidzicu - zcb * (ydb * xcb - xdb * ycb) / 
			rcbu2;
		dxidxid = rcbxu * dphidxid;
		dxidyid = rcbxu * dphidyid + zcb * rcb / ru2;
		dxidzid = rcbxu * dphidzid - ycb * rcb / ru2;
		dyidyid = rcbyu * dphidyid;
		dyidzid = rcbyu * dphidzid + xcb * rcb / ru2;
		dzidzid = rcbzu * dphidzid;

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
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedphi * dxiaxia + 
			    d2edphi2 * dphidxia * dphidxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedphi * dxiayia + 
			    d2edphi2 * dphidxia * dphidyia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedphi * dxiazia + 
			    d2edphi2 * dphidxia * dphidzia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedphi * dxiayia + 
			    d2edphi2 * dphidxia * dphidyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedphi * dyiayia + 
			    d2edphi2 * dphidyia * dphidyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedphi * dyiazia + 
			    d2edphi2 * dphidyia * dphidzia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedphi * dxiazia + 
			    d2edphi2 * dphidxia * dphidzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedphi * dyiazia + 
			    d2edphi2 * dphidyia * dphidzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedphi * dziazia + 
			    d2edphi2 * dphidzia * dphidzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedphi * dxiaxib + 
			    d2edphi2 * dphidxia * dphidxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedphi * dyiaxib + 
			    d2edphi2 * dphidyia * dphidxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedphi * dziaxib + 
			    d2edphi2 * dphidzia * dphidxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedphi * dxiayib + 
			    d2edphi2 * dphidxia * dphidyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedphi * dyiayib + 
			    d2edphi2 * dphidyia * dphidyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedphi * dziayib + 
			    d2edphi2 * dphidzia * dphidyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedphi * dxiazib + 
			    d2edphi2 * dphidxia * dphidzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedphi * dyiazib + 
			    d2edphi2 * dphidyia * dphidzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedphi * dziazib + 
			    d2edphi2 * dphidzia * dphidzib;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedphi * dxiaxic + 
			    d2edphi2 * dphidxia * dphidxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedphi * dyiaxic + 
			    d2edphi2 * dphidyia * dphidxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedphi * dziaxic + 
			    d2edphi2 * dphidzia * dphidxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedphi * dxiayic + 
			    d2edphi2 * dphidxia * dphidyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedphi * dyiayic + 
			    d2edphi2 * dphidyia * dphidyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedphi * dziayic + 
			    d2edphi2 * dphidzia * dphidyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedphi * dxiazic + 
			    d2edphi2 * dphidxia * dphidzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedphi * dyiazic + 
			    d2edphi2 * dphidyia * dphidzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedphi * dziazic + 
			    d2edphi2 * dphidzia * dphidzic;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedphi * dxiaxid + 
			    d2edphi2 * dphidxia * dphidxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedphi * dyiaxid + 
			    d2edphi2 * dphidyia * dphidxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedphi * dziaxid + 
			    d2edphi2 * dphidzia * dphidxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedphi * dxiayid + 
			    d2edphi2 * dphidxia * dphidyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedphi * dyiayid + 
			    d2edphi2 * dphidyia * dphidyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedphi * dziayid + 
			    d2edphi2 * dphidzia * dphidyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedphi * dxiazid + 
			    d2edphi2 * dphidxia * dphidzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedphi * dyiazid + 
			    d2edphi2 * dphidyia * dphidzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedphi * dziazid + 
			    d2edphi2 * dphidzia * dphidzid;
		} else if (*i__ == ib) {
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedphi * dxibxib + 
			    d2edphi2 * dphidxib * dphidxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedphi * dxibyib + 
			    d2edphi2 * dphidxib * dphidyib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedphi * dxibzib + 
			    d2edphi2 * dphidxib * dphidzib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedphi * dxibyib + 
			    d2edphi2 * dphidxib * dphidyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedphi * dyibyib + 
			    d2edphi2 * dphidyib * dphidyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedphi * dyibzib + 
			    d2edphi2 * dphidyib * dphidzib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedphi * dxibzib + 
			    d2edphi2 * dphidxib * dphidzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedphi * dyibzib + 
			    d2edphi2 * dphidyib * dphidzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedphi * dzibzib + 
			    d2edphi2 * dphidzib * dphidzib;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedphi * dxiaxib + 
			    d2edphi2 * dphidxib * dphidxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedphi * dxiayib + 
			    d2edphi2 * dphidyib * dphidxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedphi * dxiazib + 
			    d2edphi2 * dphidzib * dphidxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedphi * dyiaxib + 
			    d2edphi2 * dphidxib * dphidyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedphi * dyiayib + 
			    d2edphi2 * dphidyib * dphidyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedphi * dyiazib + 
			    d2edphi2 * dphidzib * dphidyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedphi * dziaxib + 
			    d2edphi2 * dphidxib * dphidzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedphi * dziayib + 
			    d2edphi2 * dphidyib * dphidzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedphi * dziazib + 
			    d2edphi2 * dphidzib * dphidzia;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedphi * dxibxic + 
			    d2edphi2 * dphidxib * dphidxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedphi * dyibxic + 
			    d2edphi2 * dphidyib * dphidxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedphi * dzibxic + 
			    d2edphi2 * dphidzib * dphidxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedphi * dxibyic + 
			    d2edphi2 * dphidxib * dphidyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedphi * dyibyic + 
			    d2edphi2 * dphidyib * dphidyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedphi * dzibyic + 
			    d2edphi2 * dphidzib * dphidyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedphi * dxibzic + 
			    d2edphi2 * dphidxib * dphidzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedphi * dyibzic + 
			    d2edphi2 * dphidyib * dphidzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedphi * dzibzic + 
			    d2edphi2 * dphidzib * dphidzic;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedphi * dxibxid + 
			    d2edphi2 * dphidxib * dphidxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedphi * dyibxid + 
			    d2edphi2 * dphidyib * dphidxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedphi * dzibxid + 
			    d2edphi2 * dphidzib * dphidxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedphi * dxibyid + 
			    d2edphi2 * dphidxib * dphidyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedphi * dyibyid + 
			    d2edphi2 * dphidyib * dphidyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedphi * dzibyid + 
			    d2edphi2 * dphidzib * dphidyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedphi * dxibzid + 
			    d2edphi2 * dphidxib * dphidzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedphi * dyibzid + 
			    d2edphi2 * dphidyib * dphidzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedphi * dzibzid + 
			    d2edphi2 * dphidzib * dphidzid;
		} else if (*i__ == ic) {
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedphi * dxicxic + 
			    d2edphi2 * dphidxic * dphidxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedphi * dxicyic + 
			    d2edphi2 * dphidxic * dphidyic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedphi * dxiczic + 
			    d2edphi2 * dphidxic * dphidzic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedphi * dxicyic + 
			    d2edphi2 * dphidxic * dphidyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedphi * dyicyic + 
			    d2edphi2 * dphidyic * dphidyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedphi * dyiczic + 
			    d2edphi2 * dphidyic * dphidzic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedphi * dxiczic + 
			    d2edphi2 * dphidxic * dphidzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedphi * dyiczic + 
			    d2edphi2 * dphidyic * dphidzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedphi * dziczic + 
			    d2edphi2 * dphidzic * dphidzic;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedphi * dxiaxic + 
			    d2edphi2 * dphidxic * dphidxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedphi * dxiayic + 
			    d2edphi2 * dphidyic * dphidxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedphi * dxiazic + 
			    d2edphi2 * dphidzic * dphidxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedphi * dyiaxic + 
			    d2edphi2 * dphidxic * dphidyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedphi * dyiayic + 
			    d2edphi2 * dphidyic * dphidyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedphi * dyiazic + 
			    d2edphi2 * dphidzic * dphidyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedphi * dziaxic + 
			    d2edphi2 * dphidxic * dphidzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedphi * dziayic + 
			    d2edphi2 * dphidyic * dphidzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedphi * dziazic + 
			    d2edphi2 * dphidzic * dphidzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedphi * dxibxic + 
			    d2edphi2 * dphidxic * dphidxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedphi * dxibyic + 
			    d2edphi2 * dphidyic * dphidxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedphi * dxibzic + 
			    d2edphi2 * dphidzic * dphidxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedphi * dyibxic + 
			    d2edphi2 * dphidxic * dphidyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedphi * dyibyic + 
			    d2edphi2 * dphidyic * dphidyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedphi * dyibzic + 
			    d2edphi2 * dphidzic * dphidyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedphi * dzibxic + 
			    d2edphi2 * dphidxic * dphidzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedphi * dzibyic + 
			    d2edphi2 * dphidyic * dphidzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedphi * dzibzic + 
			    d2edphi2 * dphidzic * dphidzib;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedphi * dxicxid + 
			    d2edphi2 * dphidxic * dphidxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedphi * dyicxid + 
			    d2edphi2 * dphidyic * dphidxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedphi * dzicxid + 
			    d2edphi2 * dphidzic * dphidxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedphi * dxicyid + 
			    d2edphi2 * dphidxic * dphidyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedphi * dyicyid + 
			    d2edphi2 * dphidyic * dphidyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedphi * dzicyid + 
			    d2edphi2 * dphidzic * dphidyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedphi * dxiczid + 
			    d2edphi2 * dphidxic * dphidzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedphi * dyiczid + 
			    d2edphi2 * dphidyic * dphidzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedphi * dziczid + 
			    d2edphi2 * dphidzic * dphidzid;
		} else if (*i__ == id) {
		    hessx_ref(1, id) = hessx_ref(1, id) + dedphi * dxidxid + 
			    d2edphi2 * dphidxid * dphidxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedphi * dxidyid + 
			    d2edphi2 * dphidxid * dphidyid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedphi * dxidzid + 
			    d2edphi2 * dphidxid * dphidzid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedphi * dxidyid + 
			    d2edphi2 * dphidxid * dphidyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedphi * dyidyid + 
			    d2edphi2 * dphidyid * dphidyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedphi * dyidzid + 
			    d2edphi2 * dphidyid * dphidzid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedphi * dxidzid + 
			    d2edphi2 * dphidxid * dphidzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedphi * dyidzid + 
			    d2edphi2 * dphidyid * dphidzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedphi * dzidzid + 
			    d2edphi2 * dphidzid * dphidzid;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedphi * dxiaxid + 
			    d2edphi2 * dphidxid * dphidxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedphi * dxiayid + 
			    d2edphi2 * dphidyid * dphidxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedphi * dxiazid + 
			    d2edphi2 * dphidzid * dphidxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedphi * dyiaxid + 
			    d2edphi2 * dphidxid * dphidyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedphi * dyiayid + 
			    d2edphi2 * dphidyid * dphidyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedphi * dyiazid + 
			    d2edphi2 * dphidzid * dphidyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedphi * dziaxid + 
			    d2edphi2 * dphidxid * dphidzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedphi * dziayid + 
			    d2edphi2 * dphidyid * dphidzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedphi * dziazid + 
			    d2edphi2 * dphidzid * dphidzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedphi * dxibxid + 
			    d2edphi2 * dphidxid * dphidxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedphi * dxibyid + 
			    d2edphi2 * dphidyid * dphidxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedphi * dxibzid + 
			    d2edphi2 * dphidzid * dphidxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedphi * dyibxid + 
			    d2edphi2 * dphidxid * dphidyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedphi * dyibyid + 
			    d2edphi2 * dphidyid * dphidyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedphi * dyibzid + 
			    d2edphi2 * dphidzid * dphidyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedphi * dzibxid + 
			    d2edphi2 * dphidxid * dphidzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedphi * dzibyid + 
			    d2edphi2 * dphidyid * dphidzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedphi * dzibzid + 
			    d2edphi2 * dphidzid * dphidzib;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedphi * dxicxid + 
			    d2edphi2 * dphidxid * dphidxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedphi * dxicyid + 
			    d2edphi2 * dphidyid * dphidxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedphi * dxiczid + 
			    d2edphi2 * dphidzid * dphidxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedphi * dyicxid + 
			    d2edphi2 * dphidxid * dphidyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedphi * dyicyid + 
			    d2edphi2 * dphidyid * dphidyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedphi * dyiczid + 
			    d2edphi2 * dphidzid * dphidyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedphi * dzicxid + 
			    d2edphi2 * dphidxid * dphidzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedphi * dzicyid + 
			    d2edphi2 * dphidyid * dphidzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedphi * dziczid + 
			    d2edphi2 * dphidzid * dphidzic;
		}
	    }
	}
    }
    return 0;
} /* etors2a_ */

#undef itors_ref
#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef tors6_ref
#undef tors5_ref
#undef tors4_ref
#undef tors3_ref
#undef tors2_ref
#undef tors1_ref




/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  subroutine etors2b  --  smoothed torsional Hessian  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     "etors2b" calculates the second derivatives of the torsional */
/*     energy for a single atom for use with potential energy */
/*     smoothing methods */


/* Subroutine */ int etors2b_(integer *i__)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double exp(doublereal), sin(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal dyiaxid, dyiayid, dyiazid, dziaxid, dziayid, dziazid, 
	    dxibxic, dxibyic, dxibzic, dyibxic, dyibyic, dyibzic, dzibxic, 
	    dzibyic, dzibzic, dxibxid, dxibyid, dxibzid, dyibxid, dyibyid, 
	    dyibzid, dzibxid, dzibyid, dzibzid, dxicxid, dxicyid, dxiczid, 
	    dyicxid, dyicyid, dyiczid, dzicxid, dzicyid, dziczid, dphidxia, 
	    dphidyia, dphidzia, dphidxib, dphidyib, dphidzib, dphidxic, 
	    dphidyic, dphidzic, dphidxid, dphidyid, dphidzid, c1, c2, c3, c4, 
	    c5, c6, dphidxibt, dphidyibt, dphidzibt, dphidxibu, dphidyibu, s1,
	     s2, s3, v1, v2, v3, v4, v5, v6, s4, s5, s6, dphidzibu, dphidxict,
	     dphidyict, dphidzict, dphidxicu, dphidyicu, dphidzicu;
    static integer ia, ib, ic, id;
    static doublereal xt, yt, zt, xu, yu, zu, rt2, ru2, rcb, xba, yba, zba, 
	    xcb, ycb, zcb, xdc, xia, yia, zia, xib, yib, zib, xic, yic, zic, 
	    xid, yid, zid, ydc, zdc, xca, yca, zca, xdb, ydb, zdb, xtu, ytu, 
	    ztu, fgrp, sine, rtru, damp1, damp2, damp3, damp4, damp5, damp6, 
	    dphi1, dphi2, dphi3, dphi4, dphi5, dphi6, rcbt2, rcbu2, sine2, 
	    sine3, sine4, sine5, sine6, xycb2, xzcb2, yzcb2, d2phi1, d2phi2, 
	    d2phi3, d2phi4, d2phi5, d2phi6, width, rcbxt, rcbyt, rcbzt, rcbxu,
	     rcbyu, rcbzu, wterm;
    static integer ktors;
    static doublereal dedphi, cosine;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal cosine2, cosine3, cosine4, cosine5, cosine6, d2edphi2;
    static logical proceed;
    static doublereal dxiaxia, dxiayia, dyiayia, dxibxib, dziazia, dyibyib, 
	    dzibzib, dxicxic, dyicyic, dziczic, dxidxid, dyidyid, dzidzid, 
	    dphidxt, dphidyt, dphidzt, dphidxu, dphidyu, dphidzu, dxiazia, 
	    dyiazia, dxibyib, dxibzib, dyibzib, dxicyic, dxiczic, dyiczic, 
	    dxidyid, dxidzid, dyidzid, dxiaxib, dxiayib, dxiazib, dyiaxib, 
	    dyiayib, dyiazib, dziaxib, dziayib, dziazib, dxiaxic, dxiayic, 
	    dxiazic, dyiaxic, dyiayic, dyiazic, dziaxic, dziayic, dziazic, 
	    dxiaxid, dxiayid, dxiazid;


#define tors1_ref(a_1,a_2) tors_1.tors1[(a_2)*4 + a_1 - 5]
#define tors2_ref(a_1,a_2) tors_1.tors2[(a_2)*4 + a_1 - 5]
#define tors3_ref(a_1,a_2) tors_1.tors3[(a_2)*4 + a_1 - 5]
#define tors4_ref(a_1,a_2) tors_1.tors4[(a_2)*4 + a_1 - 5]
#define tors5_ref(a_1,a_2) tors_1.tors5[(a_2)*4 + a_1 - 5]
#define tors6_ref(a_1,a_2) tors_1.tors6[(a_2)*4 + a_1 - 5]
#define hessx_ref(a_1,a_2) hessn_1.hessx[(a_2)*3 + a_1 - 4]
#define hessy_ref(a_1,a_2) hessn_1.hessy[(a_2)*3 + a_1 - 4]
#define hessz_ref(a_1,a_2) hessn_1.hessz[(a_2)*3 + a_1 - 4]
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
/*     ##  warp.i  --  parameters for potential surface smoothing  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     m2           second moment of the GDA gaussian for each atom */
/*     deform       value of the smoothing deformation parameter */
/*     difft        diffusion coefficient for torsional potential */
/*     diffv        diffusion coefficient for van der Waals potential */
/*     diffc        diffusion coefficient for charge-charge potential */
/*     use_smooth   flag to use a potential energy smoothing method */
/*     use_dem      flag to use diffusion equation method potential */
/*     use_gda      flag to use gaussian density annealing potential */
/*     use_tophat   flag to use analytical tophat smoothed potential */
/*     use_stophat  flag to use shifted tophat smoothed potential */




/*     set the extent of smoothing to be performed */

    width = warp_1.difft * warp_1.deform;
    if (width <= 0.) {
	damp1 = 1.;
	damp2 = 1.;
	damp3 = 1.;
	damp4 = 1.;
	damp5 = 1.;
	damp6 = 1.;
    } else if (warp_1.use_dem__) {
	damp1 = exp(-width);
	damp2 = exp(width * -4.);
	damp3 = exp(width * -9.);
	damp4 = exp(width * -16.);
	damp5 = exp(width * -25.);
	damp6 = exp(width * -36.);
    } else if (warp_1.use_gda__) {
	wterm = warp_1.difft / 12.;
    } else if (warp_1.use_tophat__ || warp_1.use_stophat__) {
	damp1 = 0.;
	damp2 = 0.;
	damp3 = 0.;
	damp4 = 0.;
	damp5 = 0.;
	damp6 = 0.;
	if (width < 3.141592653589793238) {
	    damp1 = sin(width) / width;
	}
	wterm = width * 2.;
	if (wterm < 3.141592653589793238) {
	    damp2 = sin(wterm) / wterm;
	}
	wterm = width * 3.;
	if (wterm < 3.141592653589793238) {
	    damp3 = sin(wterm) / wterm;
	}
	wterm = width * 4.;
	if (wterm < 3.141592653589793238) {
	    damp4 = sin(wterm) / wterm;
	}
	wterm = width * 5.;
	if (wterm < 3.141592653589793238) {
	    damp5 = sin(wterm) / wterm;
	}
	wterm = width * 6.;
	if (wterm < 3.141592653589793238) {
	    damp6 = sin(wterm) / wterm;
	}
    }

/*     calculate the torsional angle energy term */

    i__1 = tors_1.ntors;
    for (ktors = 1; ktors <= i__1; ++ktors) {
	ia = itors_ref(1, ktors);
	ib = itors_ref(2, ktors);
	ic = itors_ref(3, ktors);
	id = itors_ref(4, ktors);

/*     decide whether to compute the current interaction */

	proceed = *i__ == ia || *i__ == ib || *i__ == ic || *i__ == id;
	if (proceed && group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &c__0, &c__0);
	}

/*     compute the value of the torsional angle */

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

/*     set the torsional parameters for this angle */

		v1 = tors1_ref(1, ktors);
		c1 = tors1_ref(3, ktors);
		s1 = tors1_ref(4, ktors);
		v2 = tors2_ref(1, ktors);
		c2 = tors2_ref(3, ktors);
		s2 = tors2_ref(4, ktors);
		v3 = tors3_ref(1, ktors);
		c3 = tors3_ref(3, ktors);
		s3 = tors3_ref(4, ktors);
		v4 = tors4_ref(1, ktors);
		c4 = tors4_ref(3, ktors);
		s4 = tors4_ref(4, ktors);
		v5 = tors5_ref(1, ktors);
		c5 = tors5_ref(3, ktors);
		s5 = tors5_ref(4, ktors);
		v6 = tors6_ref(1, ktors);
		c6 = tors6_ref(3, ktors);
		s6 = tors6_ref(4, ktors);

/*     compute the multiple angle trigonometry and the phase terms */

		cosine2 = cosine * cosine - sine * sine;
		sine2 = cosine * 2. * sine;
		cosine3 = cosine * cosine2 - sine * sine2;
		sine3 = cosine * sine2 + sine * cosine2;
		cosine4 = cosine * cosine3 - sine * sine3;
		sine4 = cosine * sine3 + sine * cosine3;
		cosine5 = cosine * cosine4 - sine * sine4;
		sine5 = cosine * sine4 + sine * cosine4;
		cosine6 = cosine * cosine5 - sine * sine5;
		sine6 = cosine * sine5 + sine * cosine5;
		dphi1 = cosine * s1 - sine * c1;
		dphi2 = (cosine2 * s2 - sine2 * c2) * 2.;
		dphi3 = (cosine3 * s3 - sine3 * c3) * 3.;
		dphi4 = (cosine4 * s4 - sine4 * c4) * 4.;
		dphi5 = (cosine5 * s5 - sine5 * c5) * 5.;
		dphi6 = (cosine6 * s6 - sine6 * c6) * 6.;
		d2phi1 = -(cosine * c1 + sine * s1);
		d2phi2 = (cosine2 * c2 + sine2 * s2) * -4.;
		d2phi3 = (cosine3 * c3 + sine3 * s3) * -9.;
		d2phi4 = (cosine4 * c4 + sine4 * s4) * -16.;
		d2phi5 = (cosine5 * c5 + sine5 * s5) * -25.;
		d2phi6 = (cosine6 * c6 + sine6 * s6) * -36.;

/*     transform the potential function via smoothing */

		if (warp_1.use_gda__) {
		    width = wterm * (warp_1.m2[ia - 1] + warp_1.m2[ib - 1] + 
			    warp_1.m2[ic - 1] + warp_1.m2[id - 1]);
		    damp1 = exp(-width);
		    damp2 = exp(width * -4.);
		    damp3 = exp(width * -9.);
		    damp4 = exp(width * -16.);
		    damp5 = exp(width * -25.);
		    damp6 = exp(width * -36.);
		}
		dphi1 *= damp1;
		dphi2 *= damp2;
		dphi3 *= damp3;
		dphi4 *= damp4;
		dphi5 *= damp5;
		dphi6 *= damp6;
		d2phi1 *= damp1;
		d2phi2 *= damp2;
		d2phi3 *= damp3;
		d2phi4 *= damp4;
		d2phi5 *= damp5;
		d2phi6 *= damp6;

/*     calculate the torsional master chain rule terms */

		dedphi = torpot_1.torsunit * (v1 * dphi1 + v2 * dphi2 + v3 * 
			dphi3 + v4 * dphi4 + v5 * dphi5 + v6 * dphi6);
		d2edphi2 = torpot_1.torsunit * (v1 * d2phi1 + v2 * d2phi2 + 
			v3 * d2phi3 + v4 * d2phi4 + v5 * d2phi5 + v6 * d2phi6)
			;

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    dedphi *= fgrp;
		    d2edphi2 *= fgrp;
		}

/*     abbreviations for first derivative chain rule terms */

		xca = xic - xia;
		yca = yic - yia;
		zca = zic - zia;
		xdb = xid - xib;
		ydb = yid - yib;
		zdb = zid - zib;
		dphidxt = (yt * zcb - ycb * zt) / (rt2 * rcb);
		dphidyt = (zt * xcb - zcb * xt) / (rt2 * rcb);
		dphidzt = (xt * ycb - xcb * yt) / (rt2 * rcb);
		dphidxu = -(yu * zcb - ycb * zu) / (ru2 * rcb);
		dphidyu = -(zu * xcb - zcb * xu) / (ru2 * rcb);
		dphidzu = -(xu * ycb - xcb * yu) / (ru2 * rcb);

/*     abbreviations for second derivative chain rule terms */

		xycb2 = xcb * xcb + ycb * ycb;
		xzcb2 = xcb * xcb + zcb * zcb;
		yzcb2 = ycb * ycb + zcb * zcb;
		rcbxt = rcb * -2. * dphidxt;
		rcbyt = rcb * -2. * dphidyt;
		rcbzt = rcb * -2. * dphidzt;
		rcbt2 = rcb * rt2;
		rcbxu = rcb * 2. * dphidxu;
		rcbyu = rcb * 2. * dphidyu;
		rcbzu = rcb * 2. * dphidzu;
		rcbu2 = rcb * ru2;
		dphidxibt = yca * dphidzt - zca * dphidyt;
		dphidxibu = zdc * dphidyu - ydc * dphidzu;
		dphidyibt = zca * dphidxt - xca * dphidzt;
		dphidyibu = xdc * dphidzu - zdc * dphidxu;
		dphidzibt = xca * dphidyt - yca * dphidxt;
		dphidzibu = ydc * dphidxu - xdc * dphidyu;
		dphidxict = zba * dphidyt - yba * dphidzt;
		dphidxicu = ydb * dphidzu - zdb * dphidyu;
		dphidyict = xba * dphidzt - zba * dphidxt;
		dphidyicu = zdb * dphidxu - xdb * dphidzu;
		dphidzict = yba * dphidxt - xba * dphidyt;
		dphidzicu = xdb * dphidyu - ydb * dphidxu;

/*     chain rule terms for first derivative components */

		dphidxia = zcb * dphidyt - ycb * dphidzt;
		dphidyia = xcb * dphidzt - zcb * dphidxt;
		dphidzia = ycb * dphidxt - xcb * dphidyt;
		dphidxib = dphidxibt + dphidxibu;
		dphidyib = dphidyibt + dphidyibu;
		dphidzib = dphidzibt + dphidzibu;
		dphidxic = dphidxict + dphidxicu;
		dphidyic = dphidyict + dphidyicu;
		dphidzic = dphidzict + dphidzicu;
		dphidxid = zcb * dphidyu - ycb * dphidzu;
		dphidyid = xcb * dphidzu - zcb * dphidxu;
		dphidzid = ycb * dphidxu - xcb * dphidyu;

/*     chain rule terms for second derivative components */

		dxiaxia = rcbxt * dphidxia;
		dxiayia = rcbxt * dphidyia - zcb * rcb / rt2;
		dxiazia = rcbxt * dphidzia + ycb * rcb / rt2;
		dxiaxic = rcbxt * dphidxict + xcb * xt / rcbt2;
		dxiayic = rcbxt * dphidyict - dphidzt - (xba * zcb * xcb + 
			zba * yzcb2) / rcbt2;
		dxiazic = rcbxt * dphidzict + dphidyt + (xba * ycb * xcb + 
			yba * yzcb2) / rcbt2;
		dxiaxid = 0.;
		dxiayid = 0.;
		dxiazid = 0.;
		dyiayia = rcbyt * dphidyia;
		dyiazia = rcbyt * dphidzia - xcb * rcb / rt2;
		dyiaxib = rcbyt * dphidxibt - dphidzt - (yca * zcb * ycb + 
			zca * xzcb2) / rcbt2;
		dyiaxic = rcbyt * dphidxict + dphidzt + (yba * zcb * ycb + 
			zba * xzcb2) / rcbt2;
		dyiayic = rcbyt * dphidyict + ycb * yt / rcbt2;
		dyiazic = rcbyt * dphidzict - dphidxt - (yba * xcb * ycb + 
			xba * xzcb2) / rcbt2;
		dyiaxid = 0.;
		dyiayid = 0.;
		dyiazid = 0.;
		dziazia = rcbzt * dphidzia;
		dziaxib = rcbzt * dphidxibt + dphidyt + (zca * ycb * zcb + 
			yca * xycb2) / rcbt2;
		dziayib = rcbzt * dphidyibt - dphidxt - (zca * xcb * zcb + 
			xca * xycb2) / rcbt2;
		dziaxic = rcbzt * dphidxict - dphidyt - (zba * ycb * zcb + 
			yba * xycb2) / rcbt2;
		dziayic = rcbzt * dphidyict + dphidxt + (zba * xcb * zcb + 
			xba * xycb2) / rcbt2;
		dziazic = rcbzt * dphidzict + zcb * zt / rcbt2;
		dziaxid = 0.;
		dziayid = 0.;
		dziazid = 0.;
		dxibxic = -xcb * dphidxib / (rcb * rcb) - (yca * (zba * xcb + 
			yt) - zca * (yba * xcb - zt)) / rcbt2 - (yt * zba - 
			yba * zt) * 2. * dphidxibt / rt2 - (zdc * (ydb * xcb 
			+ zu) - ydc * (zdb * xcb - yu)) / rcbu2 + (yu * zdb - 
			ydb * zu) * 2. * dphidxibu / ru2;
		dxibyic = -ycb * dphidxib / (rcb * rcb) + dphidzt + dphidzu - 
			(yca * (zba * ycb - xt) + zca * (xba * xcb + zcb * 
			zba)) / rcbt2 - (zt * xba - zba * xt) * 2. * 
			dphidxibt / rt2 + (zdc * (xdb * xcb + zcb * zdb) + 
			ydc * (zdb * ycb + xu)) / rcbu2 + (zu * xdb - zdb * 
			xu) * 2. * dphidxibu / ru2;
		dxibxid = rcbxu * dphidxibu + xcb * xu / rcbu2;
		dxibyid = rcbyu * dphidxibu - dphidzu - (ydc * zcb * ycb + 
			zdc * xzcb2) / rcbu2;
		dxibzid = rcbzu * dphidxibu + dphidyu + (zdc * ycb * zcb + 
			ydc * xycb2) / rcbu2;
		dyibzib = ycb * dphidzib / (rcb * rcb) - (xca * (xca * xcb + 
			zcb * zca) + yca * (ycb * xca + zt)) / rcbt2 - (xt * 
			zca - xca * zt) * 2. * dphidzibt / rt2 + (ydc * (xdc *
			 ycb - zu) + xdc * (xdc * xcb + zcb * zdc)) / rcbu2 + 
			(xu * zdc - xdc * zu) * 2. * dphidzibu / ru2;
		dyibxic = -xcb * dphidyib / (rcb * rcb) - dphidzt - dphidzu + 
			(xca * (zba * xcb + yt) + zca * (zba * zcb + ycb * 
			yba)) / rcbt2 - (yt * zba - yba * zt) * 2. * 
			dphidyibt / rt2 - (zdc * (zdb * zcb + ycb * ydb) + 
			xdc * (zdb * xcb - yu)) / rcbu2 + (yu * zdb - ydb * 
			zu) * 2. * dphidyibu / ru2;
		dyibyic = -ycb * dphidyib / (rcb * rcb) - (zca * (xba * ycb + 
			zt) - xca * (zba * ycb - xt)) / rcbt2 - (zt * xba - 
			zba * xt) * 2. * dphidyibt / rt2 - (xdc * (zdb * ycb 
			+ xu) - zdc * (xdb * ycb - zu)) / rcbu2 + (zu * xdb - 
			zdb * xu) * 2. * dphidyibu / ru2;
		dyibxid = rcbxu * dphidyibu + dphidzu + (xdc * zcb * xcb + 
			zdc * yzcb2) / rcbu2;
		dyibyid = rcbyu * dphidyibu + ycb * yu / rcbu2;
		dyibzid = rcbzu * dphidyibu - dphidxu - (zdc * xcb * zcb + 
			xdc * xycb2) / rcbu2;
		dzibxic = -xcb * dphidzib / (rcb * rcb) + dphidyt + dphidyu - 
			(xca * (yba * xcb - zt) + yca * (zba * zcb + ycb * 
			yba)) / rcbt2 - (yt * zba - yba * zt) * 2. * 
			dphidzibt / rt2 + (ydc * (zdb * zcb + ycb * ydb) + 
			xdc * (ydb * xcb + zu)) / rcbu2 + (yu * zdb - ydb * 
			zu) * 2. * dphidzibu / ru2;
		dzibzic = -zcb * dphidzib / (rcb * rcb) - (xca * (yba * zcb + 
			xt) - yca * (xba * zcb - yt)) / rcbt2 - (xt * yba - 
			xba * yt) * 2. * dphidzibt / rt2 - (ydc * (xdb * zcb 
			+ yu) - xdc * (ydb * zcb - xu)) / rcbu2 + (xu * ydb - 
			xdb * yu) * 2. * dphidzibu / ru2;
		dzibxid = rcbxu * dphidzibu - dphidyu - (xdc * ycb * xcb + 
			ydc * yzcb2) / rcbu2;
		dzibyid = rcbyu * dphidzibu + dphidxu + (ydc * xcb * ycb + 
			xdc * xzcb2) / rcbu2;
		dzibzid = rcbzu * dphidzibu + zcb * zu / rcbu2;
		dxicxid = rcbxu * dphidxicu - xcb * (zdb * ycb - ydb * zcb) / 
			rcbu2;
		dxicyid = rcbyu * dphidxicu + dphidzu + (ydb * zcb * ycb + 
			zdb * xzcb2) / rcbu2;
		dxiczid = rcbzu * dphidxicu - dphidyu - (zdb * ycb * zcb + 
			ydb * xycb2) / rcbu2;
		dyicxid = rcbxu * dphidyicu - dphidzu - (xdb * zcb * xcb + 
			zdb * yzcb2) / rcbu2;
		dyicyid = rcbyu * dphidyicu - ycb * (xdb * zcb - zdb * xcb) / 
			rcbu2;
		dyiczid = rcbzu * dphidyicu + dphidxu + (zdb * xcb * zcb + 
			xdb * xycb2) / rcbu2;
		dzicxid = rcbxu * dphidzicu + dphidyu + (xdb * ycb * xcb + 
			ydb * yzcb2) / rcbu2;
		dzicyid = rcbyu * dphidzicu - dphidxu - (ydb * xcb * ycb + 
			xdb * xzcb2) / rcbu2;
		dziczid = rcbzu * dphidzicu - zcb * (ydb * xcb - xdb * ycb) / 
			rcbu2;
		dxidxid = rcbxu * dphidxid;
		dxidyid = rcbxu * dphidyid + zcb * rcb / ru2;
		dxidzid = rcbxu * dphidzid - ycb * rcb / ru2;
		dyidyid = rcbyu * dphidyid;
		dyidzid = rcbyu * dphidzid + xcb * rcb / ru2;
		dzidzid = rcbzu * dphidzid;

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
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedphi * dxiaxia + 
			    d2edphi2 * dphidxia * dphidxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedphi * dxiayia + 
			    d2edphi2 * dphidxia * dphidyia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedphi * dxiazia + 
			    d2edphi2 * dphidxia * dphidzia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedphi * dxiayia + 
			    d2edphi2 * dphidxia * dphidyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedphi * dyiayia + 
			    d2edphi2 * dphidyia * dphidyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedphi * dyiazia + 
			    d2edphi2 * dphidyia * dphidzia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedphi * dxiazia + 
			    d2edphi2 * dphidxia * dphidzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedphi * dyiazia + 
			    d2edphi2 * dphidyia * dphidzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedphi * dziazia + 
			    d2edphi2 * dphidzia * dphidzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedphi * dxiaxib + 
			    d2edphi2 * dphidxia * dphidxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedphi * dyiaxib + 
			    d2edphi2 * dphidyia * dphidxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedphi * dziaxib + 
			    d2edphi2 * dphidzia * dphidxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedphi * dxiayib + 
			    d2edphi2 * dphidxia * dphidyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedphi * dyiayib + 
			    d2edphi2 * dphidyia * dphidyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedphi * dziayib + 
			    d2edphi2 * dphidzia * dphidyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedphi * dxiazib + 
			    d2edphi2 * dphidxia * dphidzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedphi * dyiazib + 
			    d2edphi2 * dphidyia * dphidzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedphi * dziazib + 
			    d2edphi2 * dphidzia * dphidzib;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedphi * dxiaxic + 
			    d2edphi2 * dphidxia * dphidxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedphi * dyiaxic + 
			    d2edphi2 * dphidyia * dphidxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedphi * dziaxic + 
			    d2edphi2 * dphidzia * dphidxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedphi * dxiayic + 
			    d2edphi2 * dphidxia * dphidyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedphi * dyiayic + 
			    d2edphi2 * dphidyia * dphidyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedphi * dziayic + 
			    d2edphi2 * dphidzia * dphidyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedphi * dxiazic + 
			    d2edphi2 * dphidxia * dphidzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedphi * dyiazic + 
			    d2edphi2 * dphidyia * dphidzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedphi * dziazic + 
			    d2edphi2 * dphidzia * dphidzic;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedphi * dxiaxid + 
			    d2edphi2 * dphidxia * dphidxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedphi * dyiaxid + 
			    d2edphi2 * dphidyia * dphidxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedphi * dziaxid + 
			    d2edphi2 * dphidzia * dphidxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedphi * dxiayid + 
			    d2edphi2 * dphidxia * dphidyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedphi * dyiayid + 
			    d2edphi2 * dphidyia * dphidyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedphi * dziayid + 
			    d2edphi2 * dphidzia * dphidyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedphi * dxiazid + 
			    d2edphi2 * dphidxia * dphidzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedphi * dyiazid + 
			    d2edphi2 * dphidyia * dphidzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedphi * dziazid + 
			    d2edphi2 * dphidzia * dphidzid;
		} else if (*i__ == ib) {
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedphi * dxibxib + 
			    d2edphi2 * dphidxib * dphidxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedphi * dxibyib + 
			    d2edphi2 * dphidxib * dphidyib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedphi * dxibzib + 
			    d2edphi2 * dphidxib * dphidzib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedphi * dxibyib + 
			    d2edphi2 * dphidxib * dphidyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedphi * dyibyib + 
			    d2edphi2 * dphidyib * dphidyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedphi * dyibzib + 
			    d2edphi2 * dphidyib * dphidzib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedphi * dxibzib + 
			    d2edphi2 * dphidxib * dphidzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedphi * dyibzib + 
			    d2edphi2 * dphidyib * dphidzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedphi * dzibzib + 
			    d2edphi2 * dphidzib * dphidzib;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedphi * dxiaxib + 
			    d2edphi2 * dphidxib * dphidxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedphi * dxiayib + 
			    d2edphi2 * dphidyib * dphidxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedphi * dxiazib + 
			    d2edphi2 * dphidzib * dphidxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedphi * dyiaxib + 
			    d2edphi2 * dphidxib * dphidyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedphi * dyiayib + 
			    d2edphi2 * dphidyib * dphidyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedphi * dyiazib + 
			    d2edphi2 * dphidzib * dphidyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedphi * dziaxib + 
			    d2edphi2 * dphidxib * dphidzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedphi * dziayib + 
			    d2edphi2 * dphidyib * dphidzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedphi * dziazib + 
			    d2edphi2 * dphidzib * dphidzia;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedphi * dxibxic + 
			    d2edphi2 * dphidxib * dphidxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedphi * dyibxic + 
			    d2edphi2 * dphidyib * dphidxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedphi * dzibxic + 
			    d2edphi2 * dphidzib * dphidxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedphi * dxibyic + 
			    d2edphi2 * dphidxib * dphidyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedphi * dyibyic + 
			    d2edphi2 * dphidyib * dphidyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedphi * dzibyic + 
			    d2edphi2 * dphidzib * dphidyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedphi * dxibzic + 
			    d2edphi2 * dphidxib * dphidzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedphi * dyibzic + 
			    d2edphi2 * dphidyib * dphidzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedphi * dzibzic + 
			    d2edphi2 * dphidzib * dphidzic;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedphi * dxibxid + 
			    d2edphi2 * dphidxib * dphidxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedphi * dyibxid + 
			    d2edphi2 * dphidyib * dphidxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedphi * dzibxid + 
			    d2edphi2 * dphidzib * dphidxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedphi * dxibyid + 
			    d2edphi2 * dphidxib * dphidyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedphi * dyibyid + 
			    d2edphi2 * dphidyib * dphidyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedphi * dzibyid + 
			    d2edphi2 * dphidzib * dphidyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedphi * dxibzid + 
			    d2edphi2 * dphidxib * dphidzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedphi * dyibzid + 
			    d2edphi2 * dphidyib * dphidzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedphi * dzibzid + 
			    d2edphi2 * dphidzib * dphidzid;
		} else if (*i__ == ic) {
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedphi * dxicxic + 
			    d2edphi2 * dphidxic * dphidxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedphi * dxicyic + 
			    d2edphi2 * dphidxic * dphidyic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedphi * dxiczic + 
			    d2edphi2 * dphidxic * dphidzic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedphi * dxicyic + 
			    d2edphi2 * dphidxic * dphidyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedphi * dyicyic + 
			    d2edphi2 * dphidyic * dphidyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedphi * dyiczic + 
			    d2edphi2 * dphidyic * dphidzic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedphi * dxiczic + 
			    d2edphi2 * dphidxic * dphidzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedphi * dyiczic + 
			    d2edphi2 * dphidyic * dphidzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedphi * dziczic + 
			    d2edphi2 * dphidzic * dphidzic;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedphi * dxiaxic + 
			    d2edphi2 * dphidxic * dphidxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedphi * dxiayic + 
			    d2edphi2 * dphidyic * dphidxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedphi * dxiazic + 
			    d2edphi2 * dphidzic * dphidxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedphi * dyiaxic + 
			    d2edphi2 * dphidxic * dphidyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedphi * dyiayic + 
			    d2edphi2 * dphidyic * dphidyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedphi * dyiazic + 
			    d2edphi2 * dphidzic * dphidyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedphi * dziaxic + 
			    d2edphi2 * dphidxic * dphidzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedphi * dziayic + 
			    d2edphi2 * dphidyic * dphidzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedphi * dziazic + 
			    d2edphi2 * dphidzic * dphidzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedphi * dxibxic + 
			    d2edphi2 * dphidxic * dphidxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedphi * dxibyic + 
			    d2edphi2 * dphidyic * dphidxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedphi * dxibzic + 
			    d2edphi2 * dphidzic * dphidxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedphi * dyibxic + 
			    d2edphi2 * dphidxic * dphidyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedphi * dyibyic + 
			    d2edphi2 * dphidyic * dphidyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedphi * dyibzic + 
			    d2edphi2 * dphidzic * dphidyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedphi * dzibxic + 
			    d2edphi2 * dphidxic * dphidzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedphi * dzibyic + 
			    d2edphi2 * dphidyic * dphidzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedphi * dzibzic + 
			    d2edphi2 * dphidzic * dphidzib;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedphi * dxicxid + 
			    d2edphi2 * dphidxic * dphidxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedphi * dyicxid + 
			    d2edphi2 * dphidyic * dphidxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedphi * dzicxid + 
			    d2edphi2 * dphidzic * dphidxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedphi * dxicyid + 
			    d2edphi2 * dphidxic * dphidyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedphi * dyicyid + 
			    d2edphi2 * dphidyic * dphidyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedphi * dzicyid + 
			    d2edphi2 * dphidzic * dphidyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedphi * dxiczid + 
			    d2edphi2 * dphidxic * dphidzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedphi * dyiczid + 
			    d2edphi2 * dphidyic * dphidzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedphi * dziczid + 
			    d2edphi2 * dphidzic * dphidzid;
		} else if (*i__ == id) {
		    hessx_ref(1, id) = hessx_ref(1, id) + dedphi * dxidxid + 
			    d2edphi2 * dphidxid * dphidxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedphi * dxidyid + 
			    d2edphi2 * dphidxid * dphidyid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedphi * dxidzid + 
			    d2edphi2 * dphidxid * dphidzid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedphi * dxidyid + 
			    d2edphi2 * dphidxid * dphidyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedphi * dyidyid + 
			    d2edphi2 * dphidyid * dphidyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedphi * dyidzid + 
			    d2edphi2 * dphidyid * dphidzid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedphi * dxidzid + 
			    d2edphi2 * dphidxid * dphidzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedphi * dyidzid + 
			    d2edphi2 * dphidyid * dphidzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedphi * dzidzid + 
			    d2edphi2 * dphidzid * dphidzid;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedphi * dxiaxid + 
			    d2edphi2 * dphidxid * dphidxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedphi * dxiayid + 
			    d2edphi2 * dphidyid * dphidxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedphi * dxiazid + 
			    d2edphi2 * dphidzid * dphidxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedphi * dyiaxid + 
			    d2edphi2 * dphidxid * dphidyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedphi * dyiayid + 
			    d2edphi2 * dphidyid * dphidyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedphi * dyiazid + 
			    d2edphi2 * dphidzid * dphidyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedphi * dziaxid + 
			    d2edphi2 * dphidxid * dphidzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedphi * dziayid + 
			    d2edphi2 * dphidyid * dphidzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedphi * dziazid + 
			    d2edphi2 * dphidzid * dphidzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedphi * dxibxid + 
			    d2edphi2 * dphidxid * dphidxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedphi * dxibyid + 
			    d2edphi2 * dphidyid * dphidxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedphi * dxibzid + 
			    d2edphi2 * dphidzid * dphidxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedphi * dyibxid + 
			    d2edphi2 * dphidxid * dphidyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedphi * dyibyid + 
			    d2edphi2 * dphidyid * dphidyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedphi * dyibzid + 
			    d2edphi2 * dphidzid * dphidyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedphi * dzibxid + 
			    d2edphi2 * dphidxid * dphidzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedphi * dzibyid + 
			    d2edphi2 * dphidyid * dphidzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedphi * dzibzid + 
			    d2edphi2 * dphidzid * dphidzib;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedphi * dxicxid + 
			    d2edphi2 * dphidxid * dphidxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedphi * dxicyid + 
			    d2edphi2 * dphidyid * dphidxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedphi * dxiczid + 
			    d2edphi2 * dphidzid * dphidxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedphi * dyicxid + 
			    d2edphi2 * dphidxid * dphidyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedphi * dyicyid + 
			    d2edphi2 * dphidyid * dphidyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedphi * dyiczid + 
			    d2edphi2 * dphidzid * dphidyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedphi * dzicxid + 
			    d2edphi2 * dphidxid * dphidzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedphi * dzicyid + 
			    d2edphi2 * dphidyid * dphidzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedphi * dziczid + 
			    d2edphi2 * dphidzid * dphidzic;
		}
	    }
	}
    }
    return 0;
} /* etors2b_ */

#undef itors_ref
#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef tors6_ref
#undef tors5_ref
#undef tors4_ref
#undef tors3_ref
#undef tors2_ref
#undef tors1_ref


