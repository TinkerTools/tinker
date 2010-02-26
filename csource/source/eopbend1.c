/* eopbend1.f -- translated by f2c (version 20050501).
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
    doublereal opbk[75000];
    integer nopbend, iopb[75000];
} opbend_;

#define opbend_1 opbend_

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
static doublereal c_b6 = 1.;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine eopbend1  --  out-of-plane energy and derivs  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "eopbend1" computes the out-of-plane bend potential energy and */
/*     first derivatives at trigonal centers via a Wilson-Decius-Cross */
/*     or Allinger angle */


/* Subroutine */ int eopbend1_(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal), acos(doublereal), d_sign(doublereal *, 
	    doublereal *);

    /* Local variables */
    static doublereal e;
    static integer i__, ia, ib, ic, id;
    static doublereal cc, ee, dt, dt2, dt3, dt4, xia, yia, zia, xib, dot, yib,
	     zib, xic, yic, zic, xid, yid, zid, xab, yab, zab, xcb, ycb, zcb, 
	    xdb, ydb, zdb, xad, yad, zad, xcd, ycd, zcd, vxx, rab2, vyy, rad2,
	     bkk2, rdb2, rcd2, rcb2, vzz, vyx, vzx, vzy, fgrp, term;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal deddt, angle, force, dedxia, dedyia, dedcos, dedzia, 
	    dedxib, dedyib, dedzib, dedxic, dedyic, dedzic, dedxid, dedyid, 
	    dedzid, cosine;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal dccdxia, dccdyia, dccdzia, dccdxic, dccdyic, dccdzic, 
	    dccdxid, dccdyid, dccdzid, deedxia, deedyia, deedzia, deedxic, 
	    deedyic, deedzic, deedxid, deedyid;
    static integer iopbend;
    static doublereal deedzid;
    static logical proceed;


#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define deopb_ref(a_1,a_2) deriv_1.deopb[(a_2)*3 + a_1 - 4]



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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  opbend.i  --  out-of-plane bends in the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     opbk      force constant values for out-of-plane bending */
/*     nopbend   total number of out-of-plane bends in the system */
/*     iopb      bond angle numbers used in out-of-plane bending */




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




/*     zero out out-of-plane energy and first derivatives */

    energi_1.eopb = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	deopb_ref(1, i__) = 0.;
	deopb_ref(2, i__) = 0.;
	deopb_ref(3, i__) = 0.;
    }

/*     calculate the out-of-plane bending energy and derivatives */

    i__1 = opbend_1.nopbend;
    for (iopbend = 1; iopbend <= i__1; ++iopbend) {
	i__ = opbend_1.iopb[iopbend - 1];
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	id = iang_ref(4, i__);
	force = opbend_1.opbk[iopbend - 1];

/*     decide whether to compute the current interaction */

	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &c__0, &c__0);
	}
	if (proceed) {
	    proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
		    usage_1.use[ic - 1] || usage_1.use[id - 1];
	}

/*     get the coordinates of the atoms at trigonal center */

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

/*     compute the out-of-plane bending angle */

	    xab = xia - xib;
	    yab = yia - yib;
	    zab = zia - zib;
	    xcb = xic - xib;
	    ycb = yic - yib;
	    zcb = zic - zib;
	    xdb = xid - xib;
	    ydb = yid - yib;
	    zdb = zid - zib;
	    xad = xia - xid;
	    yad = yia - yid;
	    zad = zia - zid;
	    xcd = xic - xid;
	    ycd = yic - yid;
	    zcd = zic - zid;
	    if (bound_1.use_polymer__) {
		image_(&xab, &yab, &zab);
		image_(&xcb, &ycb, &zcb);
		image_(&xdb, &ydb, &zdb);
		image_(&xad, &yad, &zad);
		image_(&xcd, &ycd, &zcd);
	    }

/*     W-D-C angle between A-B-C plane and B-D vector for D-B<AC */

	    if (s_cmp(angpot_1.opbtyp, "W-D-C", (ftnlen)8, (ftnlen)5) == 0) {
		rab2 = xab * xab + yab * yab + zab * zab;
		rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
		dot = xab * xcb + yab * ycb + zab * zcb;
		cc = rab2 * rcb2 - dot * dot;

/*     Allinger angle between A-C-D plane and D-B vector for D-B<AC */

	    } else if (s_cmp(angpot_1.opbtyp, "ALLINGER", (ftnlen)8, (ftnlen)
		    8) == 0) {
		rad2 = xad * xad + yad * yad + zad * zad;
		rcd2 = xcd * xcd + ycd * ycd + zcd * zcd;
		dot = xad * xcd + yad * ycd + zad * zcd;
		cc = rad2 * rcd2 - dot * dot;
	    }

/*     find the out-of-plane angle bending energy */

	    ee = xdb * (yab * zcb - zab * ycb) + ydb * (zab * xcb - xab * zcb)
		     + zdb * (xab * ycb - yab * xcb);
	    rdb2 = xdb * xdb + ydb * ydb + zdb * zdb;
	    if (rdb2 != 0. && cc != 0.) {
		bkk2 = rdb2 - ee * ee / cc;
		cosine = sqrt(bkk2 / rdb2);
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine);
		cosine = min(d__1,d__2);
		angle = acos(cosine) * 57.29577951308232088;
		dt = angle;
		dt2 = dt * dt;
		dt3 = dt2 * dt;
		dt4 = dt2 * dt2;
		e = angpot_1.opbunit * force * dt2 * (angpot_1.copb * dt + 1. 
			+ angpot_1.qopb * dt2 + angpot_1.popb * dt3 + 
			angpot_1.sopb * dt4);
		deddt = angpot_1.opbunit * force * dt * 57.29577951308232088 *
			 (angpot_1.copb * 3. * dt + 2. + angpot_1.qopb * 4. * 
			dt2 + angpot_1.popb * 5. * dt3 + angpot_1.sopb * 6. * 
			dt4);
		dedcos = -deddt * d_sign(&c_b6, &ee) / sqrt(cc * bkk2);

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    e *= fgrp;
		    dedcos *= fgrp;
		}

/*     chain rule terms for first derivative components */

		if (s_cmp(angpot_1.opbtyp, "W-D-C", (ftnlen)8, (ftnlen)5) == 
			0) {
		    term = ee / cc;
		    dccdxia = (xab * rcb2 - xcb * dot) * term;
		    dccdyia = (yab * rcb2 - ycb * dot) * term;
		    dccdzia = (zab * rcb2 - zcb * dot) * term;
		    dccdxic = (xcb * rab2 - xab * dot) * term;
		    dccdyic = (ycb * rab2 - yab * dot) * term;
		    dccdzic = (zcb * rab2 - zab * dot) * term;
		    dccdxid = 0.;
		    dccdyid = 0.;
		    dccdzid = 0.;
		} else if (s_cmp(angpot_1.opbtyp, "ALLINGER", (ftnlen)8, (
			ftnlen)8) == 0) {
		    term = ee / cc;
		    dccdxia = (xad * rcd2 - xcd * dot) * term;
		    dccdyia = (yad * rcd2 - ycd * dot) * term;
		    dccdzia = (zad * rcd2 - zcd * dot) * term;
		    dccdxic = (xcd * rad2 - xad * dot) * term;
		    dccdyic = (ycd * rad2 - yad * dot) * term;
		    dccdzic = (zcd * rad2 - zad * dot) * term;
		    dccdxid = -dccdxia - dccdxic;
		    dccdyid = -dccdyia - dccdyic;
		    dccdzid = -dccdzia - dccdzic;
		}
		term = ee / rdb2;
		deedxia = ydb * zcb - zdb * ycb;
		deedyia = zdb * xcb - xdb * zcb;
		deedzia = xdb * ycb - ydb * xcb;
		deedxic = yab * zdb - zab * ydb;
		deedyic = zab * xdb - xab * zdb;
		deedzic = xab * ydb - yab * xdb;
		deedxid = ycb * zab - zcb * yab + xdb * term;
		deedyid = zcb * xab - xcb * zab + ydb * term;
		deedzid = xcb * yab - ycb * xab + zdb * term;

/*     compute first derivative components for this angle */

		dedxia = dedcos * (dccdxia + deedxia);
		dedyia = dedcos * (dccdyia + deedyia);
		dedzia = dedcos * (dccdzia + deedzia);
		dedxic = dedcos * (dccdxic + deedxic);
		dedyic = dedcos * (dccdyic + deedyic);
		dedzic = dedcos * (dccdzic + deedzic);
		dedxid = dedcos * (dccdxid + deedxid);
		dedyid = dedcos * (dccdyid + deedyid);
		dedzid = dedcos * (dccdzid + deedzid);
		dedxib = -dedxia - dedxic - dedxid;
		dedyib = -dedyia - dedyic - dedyid;
		dedzib = -dedzia - dedzic - dedzid;

/*     increment the out-of-plane bending energy and gradient */

		energi_1.eopb += e;
		deopb_ref(1, ia) = deopb_ref(1, ia) + dedxia;
		deopb_ref(2, ia) = deopb_ref(2, ia) + dedyia;
		deopb_ref(3, ia) = deopb_ref(3, ia) + dedzia;
		deopb_ref(1, ib) = deopb_ref(1, ib) + dedxib;
		deopb_ref(2, ib) = deopb_ref(2, ib) + dedyib;
		deopb_ref(3, ib) = deopb_ref(3, ib) + dedzib;
		deopb_ref(1, ic) = deopb_ref(1, ic) + dedxic;
		deopb_ref(2, ic) = deopb_ref(2, ic) + dedyic;
		deopb_ref(3, ic) = deopb_ref(3, ic) + dedzic;
		deopb_ref(1, id) = deopb_ref(1, id) + dedxid;
		deopb_ref(2, id) = deopb_ref(2, id) + dedyid;
		deopb_ref(3, id) = deopb_ref(3, id) + dedzid;

/*     increment the internal virial tensor components */

		vxx = xab * dedxia + xcb * dedxic + xdb * dedxid;
		vyx = yab * dedxia + ycb * dedxic + ydb * dedxid;
		vzx = zab * dedxia + zcb * dedxic + zdb * dedxid;
		vyy = yab * dedyia + ycb * dedyic + ydb * dedyid;
		vzy = zab * dedyia + zcb * dedyic + zdb * dedyid;
		vzz = zab * dedzia + zcb * dedzic + zdb * dedzid;
		vir_ref(1, 1) = vir_ref(1, 1) + vxx;
		vir_ref(2, 1) = vir_ref(2, 1) + vyx;
		vir_ref(3, 1) = vir_ref(3, 1) + vzx;
		vir_ref(1, 2) = vir_ref(1, 2) + vyx;
		vir_ref(2, 2) = vir_ref(2, 2) + vyy;
		vir_ref(3, 2) = vir_ref(3, 2) + vzy;
		vir_ref(1, 3) = vir_ref(1, 3) + vzx;
		vir_ref(2, 3) = vir_ref(2, 3) + vzy;
		vir_ref(3, 3) = vir_ref(3, 3) + vzz;
	    }
	}
    }
    return 0;
} /* eopbend1_ */

#undef deopb_ref
#undef iang_ref
#undef vir_ref


