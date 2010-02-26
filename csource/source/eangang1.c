/* eangang1.f -- translated by f2c (version 20050501).
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
    doublereal kaa[100000];
    integer nangang, iaa[200000]	/* was [2][100000] */;
} angang_;

#define angang_1 angang_

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



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine eangang1  --  angle-angle energy & derivatives  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "eangang1" calculates the angle-angle potential energy and */
/*     first derivatives with respect to Cartesian coordinates */


/* Subroutine */ int eangang1_(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal);

    /* Local variables */
    static doublereal e;
    static integer i__, k, ia, ib, ic, id, ie;
    static doublereal rp, rq, xp, yp, zp, xq, yq, zq, dt1, dt2, xia, yia, zia,
	     xib, dot, yib, zib, xic, yic, zic, xid, yid, zid, xie, yie, zie, 
	    xab, yab, zab, xcb, ycb, zcb, xdb, ydb, zdb, xeb, yeb, zeb, vxx, 
	    rab2, vyy, rcb2, rdb2, reb2, vzz, vyx, vzx, vzy, fgrp;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal angle, terma, termc, termd, terme, deddt1, deddt2, 
	    dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, dedxic, dedyic, 
	    dedzic, dedxid, dedyid, dedzid, dedxie, cosine, dedyie, dedzie;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static integer iangang;
    static logical proceed;


#define iaa_ref(a_1,a_2) angang_1.iaa[(a_2)*2 + a_1 - 3]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define deaa_ref(a_1,a_2) deriv_1.deaa[(a_2)*3 + a_1 - 4]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  angang.i  --  angle-angle terms in current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     kaa       force constant for angle-angle cross terms */
/*     nangang   total number of angle-angle interactions */
/*     iaa       angle numbers used in each angle-angle term */




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




/*     zero out the angle-angle energy and first derivatives */

    energi_1.eaa = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	deaa_ref(1, i__) = 0.;
	deaa_ref(2, i__) = 0.;
	deaa_ref(3, i__) = 0.;
    }

/*     find the energy of each angle-angle interaction */

    i__1 = angang_1.nangang;
    for (iangang = 1; iangang <= i__1; ++iangang) {
	i__ = iaa_ref(1, iangang);
	k = iaa_ref(2, iangang);
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	id = iang_ref(1, k);
	ie = iang_ref(3, k);

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
	    xid = atoms_1.x[id - 1];
	    yid = atoms_1.y[id - 1];
	    zid = atoms_1.z__[id - 1];
	    xie = atoms_1.x[ie - 1];
	    yie = atoms_1.y[ie - 1];
	    zie = atoms_1.z__[ie - 1];

/*     compute the values of the two bond angles */

	    xab = xia - xib;
	    yab = yia - yib;
	    zab = zia - zib;
	    xcb = xic - xib;
	    ycb = yic - yib;
	    zcb = zic - zib;
	    xdb = xid - xib;
	    ydb = yid - yib;
	    zdb = zid - zib;
	    xeb = xie - xib;
	    yeb = yie - yib;
	    zeb = zie - zib;
	    if (bound_1.use_polymer__) {
		image_(&xab, &yab, &zab);
		image_(&xcb, &ycb, &zcb);
		image_(&xdb, &ydb, &zdb);
		image_(&xeb, &yeb, &zeb);
	    }
	    rab2 = xab * xab + yab * yab + zab * zab;
	    rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
	    rdb2 = xdb * xdb + ydb * ydb + zdb * zdb;
	    reb2 = xeb * xeb + yeb * yeb + zeb * zeb;
	    xp = ycb * zab - zcb * yab;
	    yp = zcb * xab - xcb * zab;
	    zp = xcb * yab - ycb * xab;
	    xq = yeb * zdb - zeb * ydb;
	    yq = zeb * xdb - xeb * zdb;
	    zq = xeb * ydb - yeb * xdb;
	    rp = sqrt(xp * xp + yp * yp + zp * zp);
	    rq = sqrt(xq * xq + yq * yq + zq * zq);
	    if (rp * rq != 0.) {
		dot = xab * xcb + yab * ycb + zab * zcb;
		cosine = dot / sqrt(rab2 * rcb2);
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine);
		cosine = min(d__1,d__2);
		angle = acos(cosine) * 57.29577951308232088;
		dt1 = angle - angle_1.anat[i__ - 1];
		dot = xdb * xeb + ydb * yeb + zdb * zeb;
		cosine = dot / sqrt(rdb2 * reb2);
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine);
		cosine = min(d__1,d__2);
		angle = acos(cosine) * 57.29577951308232088;
		dt2 = angle - angle_1.anat[k - 1];

/*     get the energy and master chain rule terms for derivatives */

		e = angpot_1.aaunit * angang_1.kaa[iangang - 1] * dt1 * dt2;
		deddt1 = e * 57.29577951308232088 / dt1;
		deddt2 = e * 57.29577951308232088 / dt2;

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    e *= fgrp;
		    deddt1 *= fgrp;
		    deddt2 *= fgrp;
		}

/*     find chain rule terms for the first bond angle deviation */

		terma = -deddt1 / (rab2 * rp);
		termc = deddt1 / (rcb2 * rp);
		dedxia = terma * (yab * zp - zab * yp);
		dedyia = terma * (zab * xp - xab * zp);
		dedzia = terma * (xab * yp - yab * xp);
		dedxic = termc * (ycb * zp - zcb * yp);
		dedyic = termc * (zcb * xp - xcb * zp);
		dedzic = termc * (xcb * yp - ycb * xp);

/*     find chain rule terms for the second bond angle deviation */

		termd = -deddt2 / (rdb2 * rq);
		terme = deddt2 / (reb2 * rq);
		dedxid = termd * (ydb * zq - zdb * yq);
		dedyid = termd * (zdb * xq - xdb * zq);
		dedzid = termd * (xdb * yq - ydb * xq);
		dedxie = terme * (yeb * zq - zeb * yq);
		dedyie = terme * (zeb * xq - xeb * zq);
		dedzie = terme * (xeb * yq - yeb * xq);

/*     get the central atom derivative terms by difference */

		dedxib = -dedxia - dedxic - dedxid - dedxie;
		dedyib = -dedyia - dedyic - dedyid - dedyie;
		dedzib = -dedzia - dedzic - dedzid - dedzie;

/*     increment the total angle-angle energy and derivatives */

		energi_1.eaa += e;
		deaa_ref(1, ia) = deaa_ref(1, ia) + dedxia;
		deaa_ref(2, ia) = deaa_ref(2, ia) + dedyia;
		deaa_ref(3, ia) = deaa_ref(3, ia) + dedzia;
		deaa_ref(1, ib) = deaa_ref(1, ib) + dedxib;
		deaa_ref(2, ib) = deaa_ref(2, ib) + dedyib;
		deaa_ref(3, ib) = deaa_ref(3, ib) + dedzib;
		deaa_ref(1, ic) = deaa_ref(1, ic) + dedxic;
		deaa_ref(2, ic) = deaa_ref(2, ic) + dedyic;
		deaa_ref(3, ic) = deaa_ref(3, ic) + dedzic;
		deaa_ref(1, id) = deaa_ref(1, id) + dedxid;
		deaa_ref(2, id) = deaa_ref(2, id) + dedyid;
		deaa_ref(3, id) = deaa_ref(3, id) + dedzid;
		deaa_ref(1, ie) = deaa_ref(1, ie) + dedxie;
		deaa_ref(2, ie) = deaa_ref(2, ie) + dedyie;
		deaa_ref(3, ie) = deaa_ref(3, ie) + dedzie;

/*     increment the internal virial tensor components */

		vxx = xab * dedxia + xcb * dedxic + xdb * dedxid + xeb * 
			dedxie;
		vyx = yab * dedxia + ycb * dedxic + ydb * dedxid + yeb * 
			dedxie;
		vzx = zab * dedxia + zcb * dedxic + zdb * dedxid + zeb * 
			dedxie;
		vyy = yab * dedyia + ycb * dedyic + ydb * dedyid + yeb * 
			dedyie;
		vzy = zab * dedyia + zcb * dedyic + zdb * dedyid + zeb * 
			dedyie;
		vzz = zab * dedzia + zcb * dedzic + zdb * dedzid + zeb * 
			dedzie;
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
} /* eangang1_ */

#undef iang_ref
#undef deaa_ref
#undef vir_ref
#undef iaa_ref


