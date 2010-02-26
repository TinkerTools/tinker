/* eopdist2.f -- translated by f2c (version 20050501).
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
    doublereal opdk[25000];
    integer nopdist, iopd[100000]	/* was [4][25000] */;
} opdist_;

#define opdist_1 opdist_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine eopdist2  --  atomwise out-plane dist Hessian  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "eopdist2" calculates second derivatives of the out-of-plane */
/*     distance energy for a single atom via the central atom height */


/* Subroutine */ int eopdist2_(integer *i__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal dyiaxid, dyiayid, dyiazid, dziaxid, dziayid, dziazid, 
	    dxibxic, dxibyic, dxibzic, dyibxic, dyibyic, dyibzic, dzibxic, 
	    dzibyic, dzibzic, dxibxid, dxibyid, dxibzid, dyibxid, dyibyid, 
	    dyibzid, dzibxid, dzibyid, dzibzid, dxicxid, dxicyid, dxiczid, 
	    dyicxid, dyicyid, dyiczid, dzicxid, dzicyid, dziczid, e;
    static integer ia, ib, ic, id;
    static doublereal dt, xt, yt, zt, dt2, dt3, dt4, rt2, xad, yad, zad, xbd, 
	    ybd, xia, yia, zia, xib, yib, dot, zib, xic, yic, zic, xid, yid, 
	    zid, zbd, xcd, ycd, zcd, xtd, ytd, ztd, xyt, xzt, yxt, yzt, zxt, 
	    zyt, drt2, fgrp, ddrt, term, dotx, doty, dotz;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal deddt, force, termd, termx, termy, termz;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;
    static doublereal dxiaxia, dxiayia, dyiayia, dxibxib, dziazia, dyibyib, 
	    dzibzib, dxicxic, dyicyic, dziczic, dxidxid, dyidyid, dzidzid, 
	    dxiazia, dyiazia, dxibyib, dxibzib, dyibzib;
    static integer kopdist;
    static doublereal dxicyic, dxiczic, dyiczic, dxidyid, dxidzid, dyidzid, 
	    dxiaxib, dxiayib, dxiazib, dyiaxib, dyiayib, dyiazib, dziaxib, 
	    dziayib, dziazib, dxiaxic, dxiayic, dxiazic, dyiaxic, dyiayic, 
	    dyiazic, dziaxic, dziayic, dziazic, dxiaxid, dxiayid, dxiazid;


#define iopd_ref(a_1,a_2) opdist_1.iopd[(a_2)*4 + a_1 - 5]
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
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  opdist.i  --  out-of-plane distances in current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     opdk      force constant values for out-of-plane distance */
/*     nopdist   total number of out-of-plane distances in the system */
/*     iopb      numbers of the atoms in each out-of-plane distance */




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




/*     compute Hessian elements for the out-of-plane distances */

    i__1 = opdist_1.nopdist;
    for (kopdist = 1; kopdist <= i__1; ++kopdist) {
	ia = iopd_ref(1, kopdist);
	ib = iopd_ref(2, kopdist);
	ic = iopd_ref(3, kopdist);
	id = iopd_ref(4, kopdist);
	force = opdist_1.opdk[kopdist - 1];

/*     decide whether to compute the current interaction */

	proceed = *i__ == ia || *i__ == ib || *i__ == ic || *i__ == id;
	if (proceed && group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &c__0, &c__0);
	}

/*     get the coordinates of the defining atoms */

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

/*     compute the out-of-plane distance for central atom */

	    xad = xia - xid;
	    yad = yia - yid;
	    zad = zia - zid;
	    xbd = xib - xid;
	    ybd = yib - yid;
	    zbd = zib - zid;
	    xcd = xic - xid;
	    ycd = yic - yid;
	    zcd = zic - zid;
	    if (bound_1.use_polymer__) {
		image_(&xad, &yad, &zad);
		image_(&xbd, &ybd, &zbd);
		image_(&xcd, &ycd, &zcd);
	    }
	    xt = ybd * zcd - zbd * ycd;
	    yt = zbd * xcd - xbd * zcd;
	    zt = xbd * ycd - ybd * xcd;
	    rt2 = xt * xt + yt * yt + zt * zt;
	    dot = xt * xad + yt * yad + zt * zad;
	    drt2 = dot / rt2;
	    dt2 = dot * drt2;
	    deddt = angpot_1.opdunit * force * (angpot_1.copd * 3. * dt + 2. 
		    + angpot_1.qopd * 4. * dt2 + angpot_1.popd * 5. * dt3 + 
		    angpot_1.sopd * 6. * dt4);

/*     scale the interaction based on its group membership */

	    if (group_1.use_group__) {
		e *= fgrp;
		deddt *= fgrp;
	    }

/*     abbreviations for second derivative chain rule terms */

	    term = deddt / rt2;
	    termx = term * xt;
	    termy = term * yt;
	    termz = term * zt;
	    termd = term * dot;
	    xtd = xad - xt * 2. * drt2;
	    ytd = yad - yt * 2. * drt2;
	    ztd = zad - zt * 2. * drt2;
	    xyt = xcd * ytd - ycd * xtd;
	    xzt = xbd * ztd - zbd * xtd;
	    yxt = ybd * xtd - xbd * ytd;
	    yzt = ycd * ztd - zcd * ytd;
	    zxt = zcd * xtd - xcd * ztd;
	    zyt = zbd * ytd - ybd * ztd;
	    ddrt = dot * drt2;
	    dotx = dot * xad;
	    doty = dot * yad;
	    dotz = dot * zad;

/*     chain rule terms for second derivative components */

	    dxiaxia = termx * xt;
	    dxiayia = termx * yt;
	    dxiazia = termx * zt;
	    dxiaxib = termx * yzt;
	    dxiayib = termx * zxt + termd * zcd;
	    dxiazib = termx * xyt - termd * ycd;
	    dxiaxic = termx * zyt;
	    dxiayic = termx * xzt - termd * zbd;
	    dxiazic = termx * yxt + termd * ybd;
	    dyiayia = termy * yt;
	    dyiazia = termy * zt;
	    dyiaxib = termy * yzt - termd * zcd;
	    dyiayib = termy * zxt;
	    dyiazib = termy * xyt + termd * xcd;
	    dyiaxic = termy * zyt + termd * zbd;
	    dyiayic = termy * xzt;
	    dyiazic = termy * yxt - termd * xbd;
	    dziazia = termz * zt;
	    dziaxib = termz * yzt + termd * ycd;
	    dziayib = termz * zxt - termd * xcd;
	    dziazib = termz * xyt;
	    dziaxic = termz * zyt - termd * ybd;
	    dziayic = termz * xzt + termd * xbd;
	    dziazic = termz * yxt;
	    dxibxib = term * (yzt * yzt - ddrt * (ycd * ycd + zcd * zcd));
	    dxibyib = term * (yzt * zxt + ddrt * xcd * ycd);
	    dxibzib = term * (yzt * xyt + ddrt * xcd * zcd);
	    dxibxic = term * (yzt * zyt + ddrt * (ybd * ycd + zbd * zcd));
	    dxibyic = term * (xzt * yzt - ddrt * (xbd * ycd + zt) + dotz);
	    dxibzic = term * (yxt * yzt - ddrt * (xbd * zcd - yt) - doty);
	    dyibyib = term * (zxt * zxt - ddrt * (xcd * xcd + zcd * zcd));
	    dyibzib = term * (zxt * xyt + ddrt * ycd * zcd);
	    dyibxic = term * (zyt * zxt - ddrt * (ybd * xcd - zt) - dotz);
	    dyibyic = term * (zxt * xzt + ddrt * (xbd * xcd + zbd * zcd));
	    dyibzic = term * (yxt * zxt - ddrt * (ybd * zcd + xt) + dotx);
	    dzibzib = term * (xyt * xyt - ddrt * (xcd * xcd + ycd * ycd));
	    dzibxic = term * (zyt * xyt - ddrt * (zbd * xcd + yt) + doty);
	    dzibyic = term * (xzt * xyt - ddrt * (zbd * ycd - xt) - dotx);
	    dzibzic = term * (xyt * yxt + ddrt * (xbd * xcd + ybd * ycd));
	    dxicxic = term * (zyt * zyt - ddrt * (ybd * ybd + zbd * zbd));
	    dxicyic = term * (zyt * xzt + ddrt * xbd * ybd);
	    dxiczic = term * (zyt * yxt + ddrt * xbd * zbd);
	    dyicyic = term * (xzt * xzt - ddrt * (xbd * xbd + zbd * zbd));
	    dyiczic = term * (xzt * yxt + ddrt * ybd * zbd);
	    dziczic = term * (yxt * yxt - ddrt * (xbd * xbd + ybd * ybd));

/*     get some second derivative chain rule terms by difference */

	    dxiaxid = -dxiaxia - dxiaxib - dxiaxic;
	    dxiayid = -dxiayia - dxiayib - dxiayic;
	    dxiazid = -dxiazia - dxiazib - dxiazic;
	    dyiaxid = -dxiayia - dyiaxib - dyiaxic;
	    dyiayid = -dyiayia - dyiayib - dyiayic;
	    dyiazid = -dyiazia - dyiazib - dyiazic;
	    dziaxid = -dxiazia - dziaxib - dziaxic;
	    dziayid = -dyiazia - dziayib - dziayic;
	    dziazid = -dziazia - dziazib - dziazic;
	    dxibxid = -dxiaxib - dxibxib - dxibxic;
	    dxibyid = -dyiaxib - dxibyib - dxibyic;
	    dxibzid = -dziaxib - dxibzib - dxibzic;
	    dyibxid = -dxiayib - dxibyib - dyibxic;
	    dyibyid = -dyiayib - dyibyib - dyibyic;
	    dyibzid = -dziayib - dyibzib - dyibzic;
	    dzibxid = -dxiazib - dxibzib - dzibxic;
	    dzibyid = -dyiazib - dyibzib - dzibyic;
	    dzibzid = -dziazib - dzibzib - dzibzic;
	    dxicxid = -dxiaxic - dxibxic - dxicxic;
	    dxicyid = -dyiaxic - dyibxic - dxicyic;
	    dxiczid = -dziaxic - dzibxic - dxiczic;
	    dyicxid = -dxiayic - dxibyic - dxicyic;
	    dyicyid = -dyiayic - dyibyic - dyicyic;
	    dyiczid = -dziayic - dzibyic - dyiczic;
	    dzicxid = -dxiazic - dxibzic - dxiczic;
	    dzicyid = -dyiazic - dyibzic - dyiczic;
	    dziczid = -dziazic - dzibzic - dziczic;
	    dxidxid = -dxiaxid - dxibxid - dxicxid;
	    dxidyid = -dxiayid - dxibyid - dxicyid;
	    dxidzid = -dxiazid - dxibzid - dxiczid;
	    dyidyid = -dyiayid - dyibyid - dyicyid;
	    dyidzid = -dyiazid - dyibzid - dyiczid;
	    dzidzid = -dziazid - dzibzid - dziczid;

/*     increment diagonal and off-diagonal Hessian elements */

	    if (*i__ == ia) {
		hessx_ref(1, ia) = hessx_ref(1, ia) + dxiaxia;
		hessy_ref(1, ia) = hessy_ref(1, ia) + dxiayia;
		hessz_ref(1, ia) = hessz_ref(1, ia) + dxiazia;
		hessx_ref(2, ia) = hessx_ref(2, ia) + dxiayia;
		hessy_ref(2, ia) = hessy_ref(2, ia) + dyiayia;
		hessz_ref(2, ia) = hessz_ref(2, ia) + dyiazia;
		hessx_ref(3, ia) = hessx_ref(3, ia) + dxiazia;
		hessy_ref(3, ia) = hessy_ref(3, ia) + dyiazia;
		hessz_ref(3, ia) = hessz_ref(3, ia) + dziazia;
		hessx_ref(1, ib) = hessx_ref(1, ib) + dxiaxib;
		hessy_ref(1, ib) = hessy_ref(1, ib) + dyiaxib;
		hessz_ref(1, ib) = hessz_ref(1, ib) + dziaxib;
		hessx_ref(2, ib) = hessx_ref(2, ib) + dxiayib;
		hessy_ref(2, ib) = hessy_ref(2, ib) + dyiayib;
		hessz_ref(2, ib) = hessz_ref(2, ib) + dziayib;
		hessx_ref(3, ib) = hessx_ref(3, ib) + dxiazib;
		hessy_ref(3, ib) = hessy_ref(3, ib) + dyiazib;
		hessz_ref(3, ib) = hessz_ref(3, ib) + dziazib;
		hessx_ref(1, ic) = hessx_ref(1, ic) + dxiaxic;
		hessy_ref(1, ic) = hessy_ref(1, ic) + dyiaxic;
		hessz_ref(1, ic) = hessz_ref(1, ic) + dziaxic;
		hessx_ref(2, ic) = hessx_ref(2, ic) + dxiayic;
		hessy_ref(2, ic) = hessy_ref(2, ic) + dyiayic;
		hessz_ref(2, ic) = hessz_ref(2, ic) + dziayic;
		hessx_ref(3, ic) = hessx_ref(3, ic) + dxiazic;
		hessy_ref(3, ic) = hessy_ref(3, ic) + dyiazic;
		hessz_ref(3, ic) = hessz_ref(3, ic) + dziazic;
		hessx_ref(1, id) = hessx_ref(1, id) + dxiaxid;
		hessy_ref(1, id) = hessy_ref(1, id) + dyiaxid;
		hessz_ref(1, id) = hessz_ref(1, id) + dziaxid;
		hessx_ref(2, id) = hessx_ref(2, id) + dxiayid;
		hessy_ref(2, id) = hessy_ref(2, id) + dyiayid;
		hessz_ref(2, id) = hessz_ref(2, id) + dziayid;
		hessx_ref(3, id) = hessx_ref(3, id) + dxiazid;
		hessy_ref(3, id) = hessy_ref(3, id) + dyiazid;
		hessz_ref(3, id) = hessz_ref(3, id) + dziazid;
	    } else if (*i__ == ib) {
		hessx_ref(1, ib) = hessx_ref(1, ib) + dxibxib;
		hessy_ref(1, ib) = hessy_ref(1, ib) + dxibyib;
		hessz_ref(1, ib) = hessz_ref(1, ib) + dxibzib;
		hessx_ref(2, ib) = hessx_ref(2, ib) + dxibyib;
		hessy_ref(2, ib) = hessy_ref(2, ib) + dyibyib;
		hessz_ref(2, ib) = hessz_ref(2, ib) + dyibzib;
		hessx_ref(3, ib) = hessx_ref(3, ib) + dxibzib;
		hessy_ref(3, ib) = hessy_ref(3, ib) + dyibzib;
		hessz_ref(3, ib) = hessz_ref(3, ib) + dzibzib;
		hessx_ref(1, ia) = hessx_ref(1, ia) + dxiaxib;
		hessy_ref(1, ia) = hessy_ref(1, ia) + dxiayib;
		hessz_ref(1, ia) = hessz_ref(1, ia) + dxiazib;
		hessx_ref(2, ia) = hessx_ref(2, ia) + dyiaxib;
		hessy_ref(2, ia) = hessy_ref(2, ia) + dyiayib;
		hessz_ref(2, ia) = hessz_ref(2, ia) + dyiazib;
		hessx_ref(3, ia) = hessx_ref(3, ia) + dziaxib;
		hessy_ref(3, ia) = hessy_ref(3, ia) + dziayib;
		hessz_ref(3, ia) = hessz_ref(3, ia) + dziazib;
		hessx_ref(1, ic) = hessx_ref(1, ic) + dxibxic;
		hessy_ref(1, ic) = hessy_ref(1, ic) + dyibxic;
		hessz_ref(1, ic) = hessz_ref(1, ic) + dzibxic;
		hessx_ref(2, ic) = hessx_ref(2, ic) + dxibyic;
		hessy_ref(2, ic) = hessy_ref(2, ic) + dyibyic;
		hessz_ref(2, ic) = hessz_ref(2, ic) + dzibyic;
		hessx_ref(3, ic) = hessx_ref(3, ic) + dxibzic;
		hessy_ref(3, ic) = hessy_ref(3, ic) + dyibzic;
		hessz_ref(3, ic) = hessz_ref(3, ic) + dzibzic;
		hessx_ref(1, id) = hessx_ref(1, id) + dxibxid;
		hessy_ref(1, id) = hessy_ref(1, id) + dyibxid;
		hessz_ref(1, id) = hessz_ref(1, id) + dzibxid;
		hessx_ref(2, id) = hessx_ref(2, id) + dxibyid;
		hessy_ref(2, id) = hessy_ref(2, id) + dyibyid;
		hessz_ref(2, id) = hessz_ref(2, id) + dzibyid;
		hessx_ref(3, id) = hessx_ref(3, id) + dxibzid;
		hessy_ref(3, id) = hessy_ref(3, id) + dyibzid;
		hessz_ref(3, id) = hessz_ref(3, id) + dzibzid;
	    } else if (*i__ == ic) {
		hessx_ref(1, ic) = hessx_ref(1, ic) + dxicxic;
		hessy_ref(1, ic) = hessy_ref(1, ic) + dxicyic;
		hessz_ref(1, ic) = hessz_ref(1, ic) + dxiczic;
		hessx_ref(2, ic) = hessx_ref(2, ic) + dxicyic;
		hessy_ref(2, ic) = hessy_ref(2, ic) + dyicyic;
		hessz_ref(2, ic) = hessz_ref(2, ic) + dyiczic;
		hessx_ref(3, ic) = hessx_ref(3, ic) + dxiczic;
		hessy_ref(3, ic) = hessy_ref(3, ic) + dyiczic;
		hessz_ref(3, ic) = hessz_ref(3, ic) + dziczic;
		hessx_ref(1, ia) = hessx_ref(1, ia) + dxiaxic;
		hessy_ref(1, ia) = hessy_ref(1, ia) + dxiayic;
		hessz_ref(1, ia) = hessz_ref(1, ia) + dxiazic;
		hessx_ref(2, ia) = hessx_ref(2, ia) + dyiaxic;
		hessy_ref(2, ia) = hessy_ref(2, ia) + dyiayic;
		hessz_ref(2, ia) = hessz_ref(2, ia) + dyiazic;
		hessx_ref(3, ia) = hessx_ref(3, ia) + dziaxic;
		hessy_ref(3, ia) = hessy_ref(3, ia) + dziayic;
		hessz_ref(3, ia) = hessz_ref(3, ia) + dziazic;
		hessx_ref(1, ib) = hessx_ref(1, ib) + dxibxic;
		hessy_ref(1, ib) = hessy_ref(1, ib) + dxibyic;
		hessz_ref(1, ib) = hessz_ref(1, ib) + dxibzic;
		hessx_ref(2, ib) = hessx_ref(2, ib) + dyibxic;
		hessy_ref(2, ib) = hessy_ref(2, ib) + dyibyic;
		hessz_ref(2, ib) = hessz_ref(2, ib) + dyibzic;
		hessx_ref(3, ib) = hessx_ref(3, ib) + dzibxic;
		hessy_ref(3, ib) = hessy_ref(3, ib) + dzibyic;
		hessz_ref(3, ib) = hessz_ref(3, ib) + dzibzic;
		hessx_ref(1, id) = hessx_ref(1, id) + dxicxid;
		hessy_ref(1, id) = hessy_ref(1, id) + dyicxid;
		hessz_ref(1, id) = hessz_ref(1, id) + dzicxid;
		hessx_ref(2, id) = hessx_ref(2, id) + dxicyid;
		hessy_ref(2, id) = hessy_ref(2, id) + dyicyid;
		hessz_ref(2, id) = hessz_ref(2, id) + dzicyid;
		hessx_ref(3, id) = hessx_ref(3, id) + dxiczid;
		hessy_ref(3, id) = hessy_ref(3, id) + dyiczid;
		hessz_ref(3, id) = hessz_ref(3, id) + dziczid;
	    } else if (*i__ == id) {
		hessx_ref(1, id) = hessx_ref(1, id) + dxidxid;
		hessy_ref(1, id) = hessy_ref(1, id) + dxidyid;
		hessz_ref(1, id) = hessz_ref(1, id) + dxidzid;
		hessx_ref(2, id) = hessx_ref(2, id) + dxidyid;
		hessy_ref(2, id) = hessy_ref(2, id) + dyidyid;
		hessz_ref(2, id) = hessz_ref(2, id) + dyidzid;
		hessx_ref(3, id) = hessx_ref(3, id) + dxidzid;
		hessy_ref(3, id) = hessy_ref(3, id) + dyidzid;
		hessz_ref(3, id) = hessz_ref(3, id) + dzidzid;
		hessx_ref(1, ia) = hessx_ref(1, ia) + dxiaxid;
		hessy_ref(1, ia) = hessy_ref(1, ia) + dxiayid;
		hessz_ref(1, ia) = hessz_ref(1, ia) + dxiazid;
		hessx_ref(2, ia) = hessx_ref(2, ia) + dyiaxid;
		hessy_ref(2, ia) = hessy_ref(2, ia) + dyiayid;
		hessz_ref(2, ia) = hessz_ref(2, ia) + dyiazid;
		hessx_ref(3, ia) = hessx_ref(3, ia) + dziaxid;
		hessy_ref(3, ia) = hessy_ref(3, ia) + dziayid;
		hessz_ref(3, ia) = hessz_ref(3, ia) + dziazid;
		hessx_ref(1, ib) = hessx_ref(1, ib) + dxibxid;
		hessy_ref(1, ib) = hessy_ref(1, ib) + dxibyid;
		hessz_ref(1, ib) = hessz_ref(1, ib) + dxibzid;
		hessx_ref(2, ib) = hessx_ref(2, ib) + dyibxid;
		hessy_ref(2, ib) = hessy_ref(2, ib) + dyibyid;
		hessz_ref(2, ib) = hessz_ref(2, ib) + dyibzid;
		hessx_ref(3, ib) = hessx_ref(3, ib) + dzibxid;
		hessy_ref(3, ib) = hessy_ref(3, ib) + dzibyid;
		hessz_ref(3, ib) = hessz_ref(3, ib) + dzibzid;
		hessx_ref(1, ic) = hessx_ref(1, ic) + dxicxid;
		hessy_ref(1, ic) = hessy_ref(1, ic) + dxicyid;
		hessz_ref(1, ic) = hessz_ref(1, ic) + dxiczid;
		hessx_ref(2, ic) = hessx_ref(2, ic) + dyicxid;
		hessy_ref(2, ic) = hessy_ref(2, ic) + dyicyid;
		hessz_ref(2, ic) = hessz_ref(2, ic) + dyiczid;
		hessx_ref(3, ic) = hessx_ref(3, ic) + dzicxid;
		hessy_ref(3, ic) = hessy_ref(3, ic) + dzicyid;
		hessz_ref(3, ic) = hessz_ref(3, ic) + dziczid;
	    }
	}
    }
    return 0;
} /* eopdist2_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef iopd_ref


