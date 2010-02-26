/* eopdist1.f -- translated by f2c (version 20050501).
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
    doublereal opdk[25000];
    integer nopdist, iopd[100000]	/* was [4][25000] */;
} opdist_;

#define opdist_1 opdist_

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
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine eopdist1  --  out-of-plane dist energy & derivs  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "eopdist1" computes the out-of-plane distance potential */
/*     energy and first derivatives at trigonal centers via */
/*     the central atom height */


/* Subroutine */ int eopdist1_(void)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e;
    static integer i__, ia, ib, ic, id;
    static doublereal dt, xt, yt, zt, dt2, dt3, dt4, rt2, xia, yia, zia, xib, 
	    dot, yib, zib, xic, yic, zic, xid, yid, zid, xad, yad, zad, xbd, 
	    ybd, zbd, xcd, ycd, zcd, xtd, ytd, ztd, vxx, vyx, vyy, vzx, vzz, 
	    vzy, drt2, fgrp;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal deddt, force, dedxia, dedyia, dedzia, dedxib, dedyib, 
	    dedzib, dedxic, dedyic, dedzic, dedxid, dedyid, dedzid;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;


#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define iopd_ref(a_1,a_2) opdist_1.iopd[(a_2)*4 + a_1 - 5]
#define deopd_ref(a_1,a_2) deriv_1.deopd[(a_2)*3 + a_1 - 4]



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

    energi_1.eopd = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	deopd_ref(1, i__) = 0.;
	deopd_ref(2, i__) = 0.;
	deopd_ref(3, i__) = 0.;
    }

/*     calculate the out-of-plane distance energy and derivatives */

    i__1 = opdist_1.nopdist;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iopd_ref(1, i__);
	ib = iopd_ref(2, i__);
	ic = iopd_ref(3, i__);
	id = iopd_ref(4, i__);
	force = opdist_1.opdk[i__ - 1];

/*     decide whether to compute the current interaction */

	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &c__0, &c__0);
	}
	if (proceed) {
	    proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
		    usage_1.use[ic - 1] || usage_1.use[id - 1];
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
	    dt = sqrt(dt2);
	    dt3 = dt2 * dt;
	    dt4 = dt2 * dt2;

/*     find the out-of-plane energy and master chain rule terms */

	    e = angpot_1.opdunit * force * dt2 * (angpot_1.copd * dt + 1. + 
		    angpot_1.qopd * dt2 + angpot_1.popd * dt3 + angpot_1.sopd 
		    * dt4);
	    deddt = angpot_1.opdunit * force * drt2 * (angpot_1.copd * 3. * 
		    dt + 2. + angpot_1.qopd * 4. * dt2 + angpot_1.popd * 5. * 
		    dt3 + angpot_1.sopd * 6. * dt4);

/*     scale the interaction based on its group membership */

	    if (group_1.use_group__) {
		e *= fgrp;
		deddt *= fgrp;
	    }

/*     chain rule terms for first derivative components */

	    xtd = xad - xt * drt2;
	    ytd = yad - yt * drt2;
	    ztd = zad - zt * drt2;

/*     compute derivative components for this interaction */

	    dedxia = deddt * xt;
	    dedyia = deddt * yt;
	    dedzia = deddt * zt;
	    dedxib = deddt * (ycd * ztd - zcd * ytd);
	    dedyib = deddt * (zcd * xtd - xcd * ztd);
	    dedzib = deddt * (xcd * ytd - ycd * xtd);
	    dedxic = deddt * (zbd * ytd - ybd * ztd);
	    dedyic = deddt * (xbd * ztd - zbd * xtd);
	    dedzic = deddt * (ybd * xtd - xbd * ytd);

/*     get some derivative components by difference */

	    dedxid = -dedxia - dedxib - dedxic;
	    dedyid = -dedyia - dedyib - dedyic;
	    dedzid = -dedzia - dedzib - dedzic;

/*     increment the out-of-plane distance energy and gradient */

	    energi_1.eopd += e;
	    deopd_ref(1, ia) = deopd_ref(1, ia) + dedxia;
	    deopd_ref(2, ia) = deopd_ref(2, ia) + dedyia;
	    deopd_ref(3, ia) = deopd_ref(3, ia) + dedzia;
	    deopd_ref(1, ib) = deopd_ref(1, ib) + dedxib;
	    deopd_ref(2, ib) = deopd_ref(2, ib) + dedyib;
	    deopd_ref(3, ib) = deopd_ref(3, ib) + dedzib;
	    deopd_ref(1, ic) = deopd_ref(1, ic) + dedxic;
	    deopd_ref(2, ic) = deopd_ref(2, ic) + dedyic;
	    deopd_ref(3, ic) = deopd_ref(3, ic) + dedzic;
	    deopd_ref(1, id) = deopd_ref(1, id) + dedxid;
	    deopd_ref(2, id) = deopd_ref(2, id) + dedyid;
	    deopd_ref(3, id) = deopd_ref(3, id) + dedzid;

/*     increment the internal virial tensor components */

	    vxx = xad * dedxia + xbd * dedxib + xcd * dedxic;
	    vyx = yad * dedxia + ybd * dedxib + ycd * dedxic;
	    vzx = zad * dedxia + zbd * dedxib + zcd * dedxic;
	    vyy = yad * dedyia + ybd * dedyib + ycd * dedyic;
	    vzy = zad * dedyia + zbd * dedyib + zcd * dedyic;
	    vzz = zad * dedzia + zbd * dedzib + zcd * dedzic;
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
    return 0;
} /* eopdist1_ */

#undef deopd_ref
#undef iopd_ref
#undef vir_ref


