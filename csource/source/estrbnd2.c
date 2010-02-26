/* estrbnd2.f -- translated by f2c (version 20050501).
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
    doublereal angunit, stbnunit, aaunit, opbunit, opdunit, cang, qang, pang, 
	    sang, copb, qopb, popb, sopb, copd, qopd, popd, sopd;
    char angtyp[600000], opbtyp[8];
} angpot_;

#define angpot_1 angpot_

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
    doublereal sbk[150000]	/* was [2][75000] */;
    integer nstrbnd, isb[225000]	/* was [3][75000] */;
} strbnd_;

#define strbnd_1 strbnd_

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine estrbnd2  --  stretch-bend Hessian; analytical  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "estrbnd2" calculates the stretch-bend potential energy */
/*     second derivatives with respect to Cartesian coordinates */


/* Subroutine */ int estrbnd2_(integer *iatom)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal);

    /* Local variables */
    static doublereal drxiaxia, dtxiaxia, dtxiayia, dtxiazia, dtxibxib, 
	    dtxibyib, dtxibzib, dtxicxic, dtxicyic, dtxiczic, dtyiayia, 
	    dtyiazia, dtziazia, dtyibyib, dtyibzib, dtzibzib, dtyicyic, 
	    dtyiczic, dtziczic, dtxibxia, dtxibyia, dtxibzia, dtyibxia, 
	    dtyibyia, dtyibzia, dtzibxia, dtzibyia, dtzibzia, dtxibxic, 
	    dtxibyic, dtxibzic, dtyibxic, dtyibyic, dtyibzic, dtzibxic, 
	    dtzibyic, dtzibzic, dtxiaxic, dtxiayic, dtxiazic, dtyiaxic, 
	    dtyiayic, dtyiazic, dtziaxic;
    static integer i__, j, k;
    static doublereal dtziayic, dtziazic, drxiayia, drxiazia, drxibxib, 
	    drxibyib, drxibzib, drxicxic, drxicyic, drxiczic, dryiayia, 
	    dryiazia, drziazia, dryibyib, dryibzib, drzibzib, dryicyic, 
	    dryiczic, drziczic;
    static integer ia, ib, ic;
    static doublereal dr, dt, rp, xp, yp, zp, dr1, dr2, rp2, rab, rcb, xab, 
	    yab, zab, xcb, ycb, xia, yia, zia, xib, yib, dot, zib, xic, yic, 
	    zic, zcb, rab2, rcb2, xabp, yabp, xrab, yrab, fgrp, zrab, term, 
	    xrcb, yrcb, zrcb, zabp, xcbp, ycbp, zcbp, term1, term2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal angle, force1, force2, cosine;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal ddrdxia, ddrdyia, ddtdxia, ddtdyia, ddtdzia, ddtdxib, 
	    ddtdyib, ddtdzib, ddtdxic, ddtdyic, ddtdzic, ddrdzia, ddrdxib, 
	    ddrdyib, ddrdzib, ddrdxic, ddrdyic, ddrdzic;
    static logical proceed;
    static integer istrbnd;


#define isb_ref(a_1,a_2) strbnd_1.isb[(a_2)*3 + a_1 - 4]
#define sbk_ref(a_1,a_2) strbnd_1.sbk[(a_2)*2 + a_1 - 3]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define hessx_ref(a_1,a_2) hessn_1.hessx[(a_2)*3 + a_1 - 4]
#define hessy_ref(a_1,a_2) hessn_1.hessy[(a_2)*3 + a_1 - 4]
#define hessz_ref(a_1,a_2) hessn_1.hessz[(a_2)*3 + a_1 - 4]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  strbnd.i  --  stretch-bends in the current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     sbk       force constants for stretch-bend terms */
/*     nstrbnd   total number of stretch-bend interactions */
/*     isb       angle and bond numbers used in stretch-bend */




/*     compute the Hessian elements of the stretch-bends */

    i__1 = strbnd_1.nstrbnd;
    for (istrbnd = 1; istrbnd <= i__1; ++istrbnd) {
	i__ = isb_ref(1, istrbnd);
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	force1 = sbk_ref(1, istrbnd);
	force2 = sbk_ref(2, istrbnd);

/*     decide whether to compute the current interaction */

	proceed = *iatom == ia || *iatom == ib || *iatom == ic;
	if (proceed && group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &c__0, &c__0, &c__0);
	}

/*     get the coordinates of the atoms in the angle */

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

/*     compute the value of the bond angle */

	    xab = xia - xib;
	    yab = yia - yib;
	    zab = zia - zib;
	    xcb = xic - xib;
	    ycb = yic - yib;
	    zcb = zic - zib;
	    if (bound_1.use_polymer__) {
		image_(&xab, &yab, &zab);
		image_(&xcb, &ycb, &zcb);
	    }
	    rab2 = xab * xab + yab * yab + zab * zab;
	    rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
	    if (rab2 != 0. && rcb2 != 0.) {
		rab = sqrt(rab2);
		rcb = sqrt(rcb2);
		xp = ycb * zab - zcb * yab;
		yp = zcb * xab - xcb * zab;
		zp = xcb * yab - ycb * xab;
		rp = sqrt(xp * xp + yp * yp + zp * zp);
		rp = max(rp,.001);
		dot = xab * xcb + yab * ycb + zab * zcb;
		cosine = dot / (rab * rcb);
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine);
		cosine = min(d__1,d__2);
		angle = acos(cosine) * 57.29577951308232088;

/*     first derivatives of angle with respect to coordinates */

		dt = angle - angle_1.anat[i__ - 1];
		term1 = -57.29577951308232088 / (rab2 * rp);
		term2 = 57.29577951308232088 / (rcb2 * rp);
		ddtdxia = term1 * (yab * zp - zab * yp);
		ddtdyia = term1 * (zab * xp - xab * zp);
		ddtdzia = term1 * (xab * yp - yab * xp);
		ddtdxic = term2 * (ycb * zp - zcb * yp);
		ddtdyic = term2 * (zcb * xp - xcb * zp);
		ddtdzic = term2 * (xcb * yp - ycb * xp);
		ddtdxib = -ddtdxia - ddtdxic;
		ddtdyib = -ddtdyia - ddtdyic;
		ddtdzib = -ddtdzia - ddtdzic;

/*     abbreviations used in defining chain rule terms */

		xrab = xab * 2. / rab2;
		yrab = yab * 2. / rab2;
		zrab = zab * 2. / rab2;
		xrcb = xcb * 2. / rcb2;
		yrcb = ycb * 2. / rcb2;
		zrcb = zcb * 2. / rcb2;
		rp2 = 1. / (rp * rp);
		xabp = (yab * zp - zab * yp) * rp2;
		yabp = (zab * xp - xab * zp) * rp2;
		zabp = (xab * yp - yab * xp) * rp2;
		xcbp = (ycb * zp - zcb * yp) * rp2;
		ycbp = (zcb * xp - xcb * zp) * rp2;
		zcbp = (xcb * yp - ycb * xp) * rp2;

/*     second derivatives of angle with respect to coordinates */

		dtxiaxia = term1 * (xab * xcb - dot) + ddtdxia * (xcbp - xrab)
			;
		dtxiayia = term1 * (zp + yab * xcb) + ddtdxia * (ycbp - yrab);
		dtxiazia = term1 * (zab * xcb - yp) + ddtdxia * (zcbp - zrab);
		dtyiayia = term1 * (yab * ycb - dot) + ddtdyia * (ycbp - yrab)
			;
		dtyiazia = term1 * (xp + zab * ycb) + ddtdyia * (zcbp - zrab);
		dtziazia = term1 * (zab * zcb - dot) + ddtdzia * (zcbp - zrab)
			;
		dtxicxic = term2 * (dot - xab * xcb) - ddtdxic * (xabp + xrcb)
			;
		dtxicyic = term2 * (zp - ycb * xab) - ddtdxic * (yabp + yrcb);
		dtxiczic = -term2 * (yp + zcb * xab) - ddtdxic * (zabp + zrcb)
			;
		dtyicyic = term2 * (dot - yab * ycb) - ddtdyic * (yabp + yrcb)
			;
		dtyiczic = term2 * (xp - zcb * yab) - ddtdyic * (zabp + zrcb);
		dtziczic = term2 * (dot - zab * zcb) - ddtdzic * (zabp + zrcb)
			;
		dtxiaxic = term1 * (yab * yab + zab * zab) - ddtdxia * xabp;
		dtxiayic = -term1 * xab * yab - ddtdxia * yabp;
		dtxiazic = -term1 * xab * zab - ddtdxia * zabp;
		dtyiaxic = -term1 * xab * yab - ddtdyia * xabp;
		dtyiayic = term1 * (xab * xab + zab * zab) - ddtdyia * yabp;
		dtyiazic = -term1 * yab * zab - ddtdyia * zabp;
		dtziaxic = -term1 * xab * zab - ddtdzia * xabp;
		dtziayic = -term1 * yab * zab - ddtdzia * yabp;
		dtziazic = term1 * (xab * xab + yab * yab) - ddtdzia * zabp;

/*     more angle deviation derivatives resulting from symmetry */

		dtxibxia = -dtxiaxia - dtxiaxic;
		dtxibyia = -dtxiayia - dtyiaxic;
		dtxibzia = -dtxiazia - dtziaxic;
		dtyibxia = -dtxiayia - dtxiayic;
		dtyibyia = -dtyiayia - dtyiayic;
		dtyibzia = -dtyiazia - dtziayic;
		dtzibxia = -dtxiazia - dtxiazic;
		dtzibyia = -dtyiazia - dtyiazic;
		dtzibzia = -dtziazia - dtziazic;
		dtxibxic = -dtxicxic - dtxiaxic;
		dtxibyic = -dtxicyic - dtxiayic;
		dtxibzic = -dtxiczic - dtxiazic;
		dtyibxic = -dtxicyic - dtyiaxic;
		dtyibyic = -dtyicyic - dtyiayic;
		dtyibzic = -dtyiczic - dtyiazic;
		dtzibxic = -dtxiczic - dtziaxic;
		dtzibyic = -dtyiczic - dtziayic;
		dtzibzic = -dtziczic - dtziazic;
		dtxibxib = -dtxibxia - dtxibxic;
		dtxibyib = -dtxibyia - dtxibyic;
		dtxibzib = -dtxibzia - dtxibzic;
		dtyibyib = -dtyibyia - dtyibyic;
		dtyibzib = -dtyibzia - dtyibzic;
		dtzibzib = -dtzibzia - dtzibzic;

/*     compute the values of the bond length deviations */

		j = isb_ref(2, istrbnd);
		k = isb_ref(3, istrbnd);
		term = angpot_1.stbnunit * force1;
		dr1 = term * (rab - bond_1.bl[j - 1]);
		term1 = term / rab;
		term = angpot_1.stbnunit * force2;
		dr2 = term * (rcb - bond_1.bl[k - 1]);
		term2 = term / rcb;
		dr = dr1 + dr2;

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    dr *= fgrp;
		    dr1 *= fgrp;
		    dr2 *= fgrp;
		    term1 *= fgrp;
		    term2 *= fgrp;
		}

/*     first derivatives of bond length with respect to coordinates */

		ddrdxia = term1 * xab;
		ddrdyia = term1 * yab;
		ddrdzia = term1 * zab;
		ddrdxic = term2 * xcb;
		ddrdyic = term2 * ycb;
		ddrdzic = term2 * zcb;
		ddrdxib = -ddrdxia - ddrdxic;
		ddrdyib = -ddrdyia - ddrdyic;
		ddrdzib = -ddrdzia - ddrdzic;

/*     abbreviations used in defining chain rule terms */

		xab /= rab;
		yab /= rab;
		zab /= rab;
		xcb /= rcb;
		ycb /= rcb;
		zcb /= rcb;

/*     second derivatives of bond length with respect to coordinates */

		drxiaxia = term1 * (1. - xab * xab);
		drxiayia = -term1 * xab * yab;
		drxiazia = -term1 * xab * zab;
		dryiayia = term1 * (1. - yab * yab);
		dryiazia = -term1 * yab * zab;
		drziazia = term1 * (1. - zab * zab);
		drxicxic = term2 * (1. - xcb * xcb);
		drxicyic = -term2 * xcb * ycb;
		drxiczic = -term2 * xcb * zcb;
		dryicyic = term2 * (1. - ycb * ycb);
		dryiczic = -term2 * ycb * zcb;
		drziczic = term2 * (1. - zcb * zcb);
		drxibxib = drxiaxia + drxicxic;
		drxibyib = drxiayia + drxicyic;
		drxibzib = drxiazia + drxiczic;
		dryibyib = dryiayia + dryicyic;
		dryibzib = dryiazia + dryiczic;
		drzibzib = drziazia + drziczic;

/*     increment diagonal and non-diagonal Hessian elements */

		if (ia == *iatom) {
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dt * drxiaxia + dr *
			     dtxiaxia + ddtdxia * 2. * ddrdxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dt * drxiayia + dr *
			     dtxiayia + ddtdxia * ddrdyia + ddtdyia * ddrdxia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dt * drxiazia + dr *
			     dtxiazia + ddtdxia * ddrdzia + ddtdzia * ddrdxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dt * drxiayia + dr *
			     dtxiayia + ddtdyia * ddrdxia + ddtdxia * ddrdyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dt * dryiayia + dr *
			     dtyiayia + ddtdyia * 2. * ddrdyia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dt * dryiazia + dr *
			     dtyiazia + ddtdyia * ddrdzia + ddtdzia * ddrdyia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dt * drxiazia + dr *
			     dtxiazia + ddtdzia * ddrdxia + ddtdxia * ddrdzia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dt * dryiazia + dr *
			     dtyiazia + ddtdzia * ddrdyia + ddtdyia * ddrdzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dt * drziazia + dr *
			     dtziazia + ddtdzia * 2. * ddrdzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) - dt * drxiaxia + dr *
			     dtxibxia + ddtdxia * ddrdxib + ddtdxib * ddrdxia;
		    hessx_ref(2, ib) = hessx_ref(2, ib) - dt * drxiayia + dr *
			     dtxibyia + ddtdxia * ddrdyib + ddtdyib * ddrdxia;
		    hessx_ref(3, ib) = hessx_ref(3, ib) - dt * drxiazia + dr *
			     dtxibzia + ddtdxia * ddrdzib + ddtdzib * ddrdxia;
		    hessy_ref(1, ib) = hessy_ref(1, ib) - dt * drxiayia + dr *
			     dtyibxia + ddtdyia * ddrdxib + ddtdxib * ddrdyia;
		    hessy_ref(2, ib) = hessy_ref(2, ib) - dt * dryiayia + dr *
			     dtyibyia + ddtdyia * ddrdyib + ddtdyib * ddrdyia;
		    hessy_ref(3, ib) = hessy_ref(3, ib) - dt * dryiazia + dr *
			     dtyibzia + ddtdyia * ddrdzib + ddtdzib * ddrdyia;
		    hessz_ref(1, ib) = hessz_ref(1, ib) - dt * drxiazia + dr *
			     dtzibxia + ddtdzia * ddrdxib + ddtdxib * ddrdzia;
		    hessz_ref(2, ib) = hessz_ref(2, ib) - dt * dryiazia + dr *
			     dtzibyia + ddtdzia * ddrdyib + ddtdyib * ddrdzia;
		    hessz_ref(3, ib) = hessz_ref(3, ib) - dt * drziazia + dr *
			     dtzibzia + ddtdzia * ddrdzib + ddtdzib * ddrdzia;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dr * dtxiaxic + 
			    ddtdxia * ddrdxic + ddtdxic * ddrdxia;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dr * dtxiayic + 
			    ddtdxia * ddrdyic + ddtdyic * ddrdxia;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dr * dtxiazic + 
			    ddtdxia * ddrdzic + ddtdzic * ddrdxia;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dr * dtyiaxic + 
			    ddtdyia * ddrdxic + ddtdxic * ddrdyia;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dr * dtyiayic + 
			    ddtdyia * ddrdyic + ddtdyic * ddrdyia;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dr * dtyiazic + 
			    ddtdyia * ddrdzic + ddtdzic * ddrdyia;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dr * dtziaxic + 
			    ddtdzia * ddrdxic + ddtdxic * ddrdzia;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dr * dtziayic + 
			    ddtdzia * ddrdyic + ddtdyic * ddrdzia;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dr * dtziazic + 
			    ddtdzia * ddrdzic + ddtdzic * ddrdzia;
		} else if (ib == *iatom) {
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dt * drxibxib + dr *
			     dtxibxib + ddtdxib * 2. * ddrdxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dt * drxibyib + dr *
			     dtxibyib + ddtdxib * ddrdyib + ddtdyib * ddrdxib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dt * drxibzib + dr *
			     dtxibzib + ddtdxib * ddrdzib + ddtdzib * ddrdxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dt * drxibyib + dr *
			     dtxibyib + ddtdyib * ddrdxib + ddtdxib * ddrdyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dt * dryibyib + dr *
			     dtyibyib + ddtdyib * 2. * ddrdyib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dt * dryibzib + dr *
			     dtyibzib + ddtdyib * ddrdzib + ddtdzib * ddrdyib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dt * drxibzib + dr *
			     dtxibzib + ddtdzib * ddrdxib + ddtdxib * ddrdzib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dt * dryibzib + dr *
			     dtyibzib + ddtdzib * ddrdyib + ddtdyib * ddrdzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dt * drzibzib + dr *
			     dtzibzib + ddtdzib * 2. * ddrdzib;
		    hessx_ref(1, ia) = hessx_ref(1, ia) - dt * drxiaxia + dr *
			     dtxibxia + ddtdxib * ddrdxia + ddtdxia * ddrdxib;
		    hessx_ref(2, ia) = hessx_ref(2, ia) - dt * drxiayia + dr *
			     dtxibyia + ddtdxib * ddrdyia + ddtdyia * ddrdxib;
		    hessx_ref(3, ia) = hessx_ref(3, ia) - dt * drxiazia + dr *
			     dtxibzia + ddtdxib * ddrdzia + ddtdzia * ddrdxib;
		    hessy_ref(1, ia) = hessy_ref(1, ia) - dt * drxiayia + dr *
			     dtyibxia + ddtdyib * ddrdxia + ddtdxia * ddrdyib;
		    hessy_ref(2, ia) = hessy_ref(2, ia) - dt * dryiayia + dr *
			     dtyibyia + ddtdyib * ddrdyia + ddtdyia * ddrdyib;
		    hessy_ref(3, ia) = hessy_ref(3, ia) - dt * dryiazia + dr *
			     dtyibzia + ddtdyib * ddrdzia + ddtdzia * ddrdyib;
		    hessz_ref(1, ia) = hessz_ref(1, ia) - dt * drxiazia + dr *
			     dtzibxia + ddtdzib * ddrdxia + ddtdxia * ddrdzib;
		    hessz_ref(2, ia) = hessz_ref(2, ia) - dt * dryiazia + dr *
			     dtzibyia + ddtdzib * ddrdyia + ddtdyia * ddrdzib;
		    hessz_ref(3, ia) = hessz_ref(3, ia) - dt * drziazia + dr *
			     dtzibzia + ddtdzib * ddrdzia + ddtdzia * ddrdzib;
		    hessx_ref(1, ic) = hessx_ref(1, ic) - dt * drxicxic + dr *
			     dtxibxic + ddtdxib * ddrdxic + ddtdxic * ddrdxib;
		    hessx_ref(2, ic) = hessx_ref(2, ic) - dt * drxicyic + dr *
			     dtxibyic + ddtdxib * ddrdyic + ddtdyic * ddrdxib;
		    hessx_ref(3, ic) = hessx_ref(3, ic) - dt * drxiczic + dr *
			     dtxibzic + ddtdxib * ddrdzic + ddtdzic * ddrdxib;
		    hessy_ref(1, ic) = hessy_ref(1, ic) - dt * drxicyic + dr *
			     dtyibxic + ddtdyib * ddrdxic + ddtdxic * ddrdyib;
		    hessy_ref(2, ic) = hessy_ref(2, ic) - dt * dryicyic + dr *
			     dtyibyic + ddtdyib * ddrdyic + ddtdyic * ddrdyib;
		    hessy_ref(3, ic) = hessy_ref(3, ic) - dt * dryiczic + dr *
			     dtyibzic + ddtdyib * ddrdzic + ddtdzic * ddrdyib;
		    hessz_ref(1, ic) = hessz_ref(1, ic) - dt * drxiczic + dr *
			     dtzibxic + ddtdzib * ddrdxic + ddtdxic * ddrdzib;
		    hessz_ref(2, ic) = hessz_ref(2, ic) - dt * dryiczic + dr *
			     dtzibyic + ddtdzib * ddrdyic + ddtdyic * ddrdzib;
		    hessz_ref(3, ic) = hessz_ref(3, ic) - dt * drziczic + dr *
			     dtzibzic + ddtdzib * ddrdzic + ddtdzic * ddrdzib;
		} else if (ic == *iatom) {
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dt * drxicxic + dr *
			     dtxicxic + ddtdxic * 2. * ddrdxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dt * drxicyic + dr *
			     dtxicyic + ddtdxic * ddrdyic + ddtdyic * ddrdxic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dt * drxiczic + dr *
			     dtxiczic + ddtdxic * ddrdzic + ddtdzic * ddrdxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dt * drxicyic + dr *
			     dtxicyic + ddtdyic * ddrdxic + ddtdxic * ddrdyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dt * dryicyic + dr *
			     dtyicyic + ddtdyic * 2. * ddrdyic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dt * dryiczic + dr *
			     dtyiczic + ddtdyic * ddrdzic + ddtdzic * ddrdyic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dt * drxiczic + dr *
			     dtxiczic + ddtdzic * ddrdxic + ddtdxic * ddrdzic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dt * dryiczic + dr *
			     dtyiczic + ddtdzic * ddrdyic + ddtdyic * ddrdzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dt * drziczic + dr *
			     dtziczic + ddtdzic * 2. * ddrdzic;
		    hessx_ref(1, ib) = hessx_ref(1, ib) - dt * drxicxic + dr *
			     dtxibxic + ddtdxic * ddrdxib + ddtdxib * ddrdxic;
		    hessx_ref(2, ib) = hessx_ref(2, ib) - dt * drxicyic + dr *
			     dtxibyic + ddtdxic * ddrdyib + ddtdyib * ddrdxic;
		    hessx_ref(3, ib) = hessx_ref(3, ib) - dt * drxiczic + dr *
			     dtxibzic + ddtdxic * ddrdzib + ddtdzib * ddrdxic;
		    hessy_ref(1, ib) = hessy_ref(1, ib) - dt * drxicyic + dr *
			     dtyibxic + ddtdyic * ddrdxib + ddtdxib * ddrdyic;
		    hessy_ref(2, ib) = hessy_ref(2, ib) - dt * dryicyic + dr *
			     dtyibyic + ddtdyic * ddrdyib + ddtdyib * ddrdyic;
		    hessy_ref(3, ib) = hessy_ref(3, ib) - dt * dryiczic + dr *
			     dtyibzic + ddtdyic * ddrdzib + ddtdzib * ddrdyic;
		    hessz_ref(1, ib) = hessz_ref(1, ib) - dt * drxiczic + dr *
			     dtzibxic + ddtdzic * ddrdxib + ddtdxib * ddrdzic;
		    hessz_ref(2, ib) = hessz_ref(2, ib) - dt * dryiczic + dr *
			     dtzibyic + ddtdzic * ddrdyib + ddtdyib * ddrdzic;
		    hessz_ref(3, ib) = hessz_ref(3, ib) - dt * drziczic + dr *
			     dtzibzic + ddtdzic * ddrdzib + ddtdzib * ddrdzic;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dr * dtxiaxic + 
			    ddtdxic * ddrdxia + ddtdxia * ddrdxic;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dr * dtyiaxic + 
			    ddtdxic * ddrdyia + ddtdyia * ddrdxic;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dr * dtziaxic + 
			    ddtdxic * ddrdzia + ddtdzia * ddrdxic;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dr * dtxiayic + 
			    ddtdyic * ddrdxia + ddtdxia * ddrdyic;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dr * dtyiayic + 
			    ddtdyic * ddrdyia + ddtdyia * ddrdyic;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dr * dtziayic + 
			    ddtdyic * ddrdzia + ddtdzia * ddrdyic;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dr * dtxiazic + 
			    ddtdzic * ddrdxia + ddtdxia * ddrdzic;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dr * dtyiazic + 
			    ddtdzic * ddrdyia + ddtdyia * ddrdzic;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dr * dtziazic + 
			    ddtdzic * ddrdzia + ddtdzia * ddrdzic;
		}
	    }
	}
    }
    return 0;
} /* estrbnd2_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef iang_ref
#undef sbk_ref
#undef isb_ref


