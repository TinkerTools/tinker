/* eangle.f -- translated by f2c (version 20050501).
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

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine eangle  --  angle bending potential energy  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "eangle" calculates the angle bending potential energy; */
/*     projected in-plane angles at trigonal centers, special */
/*     linear or Fourier angle bending terms are optionally used */


/* Subroutine */ int eangle_(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal), acos(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal e;
    static integer i__, ia, ib, ic, id;
    static doublereal dt, xt, yt, zt, dt2, dt3, dt4, rt2, xab, yab, zab, xcb, 
	    ycb, xia, yia, zia, xib, dot, yib, zib, xic, yic, zic, xid, yid, 
	    zid, zcb, xad, yad, zad, xbd, ybd, zbd, xcd, ycd, zcd, xip, yip, 
	    zip, xap, yap, zap, xcp, ycp, zcp, rab2, rcb2, rap2, rcp2, fold, 
	    fgrp, ideal;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal angle, delta, force, factor, cosine;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;


#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define angtyp_ref(a_0,a_1) &angpot_1.angtyp[(a_1)*8 + a_0 - 8]



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




/*     zero out the angle bending energy component */

    energi_1.ea = 0.;

/*     calculate the bond angle bending energy term */

    i__1 = angle_1.nangle;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	id = iang_ref(4, i__);
	ideal = angle_1.anat[i__ - 1];
	force = angle_1.ak[i__ - 1];

/*     decide whether to compute the current interaction */

	proceed = TRUE_;
	if (s_cmp(angtyp_ref(0, i__), "IN-PLANE", (ftnlen)8, (ftnlen)8) == 0) 
		{
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &c__0, &c__0);
	    }
	    if (proceed) {
		proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
			usage_1.use[ic - 1] || usage_1.use[id - 1];
	    }
	} else {
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &ia, &ib, &ic, &c__0, &c__0, &c__0);
	    }
	    if (proceed) {
		proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
			usage_1.use[ic - 1];
	    }
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

/*     compute the bond angle bending energy */

	    if (s_cmp(angtyp_ref(0, i__), "IN-PLANE", (ftnlen)8, (ftnlen)8) !=
		     0) {
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
		    dot = xab * xcb + yab * ycb + zab * zcb;
		    cosine = dot / sqrt(rab2 * rcb2);
/* Computing MIN */
		    d__1 = 1., d__2 = max(-1.,cosine);
		    cosine = min(d__1,d__2);
		    angle = acos(cosine) * 57.29577951308232088;
		    if (s_cmp(angtyp_ref(0, i__), "HARMONIC", (ftnlen)8, (
			    ftnlen)8) == 0) {
			dt = angle - ideal;
			dt2 = dt * dt;
			dt3 = dt2 * dt;
			dt4 = dt2 * dt2;
			e = angpot_1.angunit * force * dt2 * (angpot_1.cang * 
				dt + 1. + angpot_1.qang * dt2 + angpot_1.pang 
				* dt3 + angpot_1.sang * dt4);
		    } else if (s_cmp(angtyp_ref(0, i__), "LINEAR", (ftnlen)8, 
			    (ftnlen)6) == 0) {
			factor = angpot_1.angunit * 2. * 3282.8063500117441;
			e = factor * force * (cosine + 1.);
		    } else if (s_cmp(angtyp_ref(0, i__), "FOURIER", (ftnlen)8,
			     (ftnlen)7) == 0) {
			fold = angle_1.afld[i__ - 1];
/* Computing 2nd power */
			d__1 = 57.29577951308232088 / fold;
			factor = angpot_1.angunit * 2. * (d__1 * d__1);
			cosine = cos((fold * angle - ideal) / 
				57.29577951308232088);
			e = factor * force * (cosine + 1.);
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
		    }

/*     increment the total bond angle bending energy */

		    energi_1.ea += e;
		}

/*     compute the projected in-plane angle bend energy */

	    } else {
		xid = atoms_1.x[id - 1];
		yid = atoms_1.y[id - 1];
		zid = atoms_1.z__[id - 1];
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
		xt = yad * zcd - zad * ycd;
		yt = zad * xcd - xad * zcd;
		zt = xad * ycd - yad * xcd;
		rt2 = xt * xt + yt * yt + zt * zt;
		delta = -(xt * xbd + yt * ybd + zt * zbd) / rt2;
		xip = xib + xt * delta;
		yip = yib + yt * delta;
		zip = zib + zt * delta;
		xap = xia - xip;
		yap = yia - yip;
		zap = zia - zip;
		xcp = xic - xip;
		ycp = yic - yip;
		zcp = zic - zip;
		if (bound_1.use_polymer__) {
		    image_(&xap, &yap, &zap);
		    image_(&xcp, &ycp, &zcp);
		}
		rap2 = xap * xap + yap * yap + zap * zap;
		rcp2 = xcp * xcp + ycp * ycp + zcp * zcp;
		if (rap2 != 0. && rcp2 != 0.) {
		    dot = xap * xcp + yap * ycp + zap * zcp;
		    cosine = dot / sqrt(rap2 * rcp2);
/* Computing MIN */
		    d__1 = 1., d__2 = max(-1.,cosine);
		    cosine = min(d__1,d__2);
		    angle = acos(cosine) * 57.29577951308232088;
		    dt = angle - ideal;
		    dt2 = dt * dt;
		    dt3 = dt2 * dt;
		    dt4 = dt2 * dt2;
		    e = angpot_1.angunit * force * dt2 * (angpot_1.cang * dt 
			    + 1. + angpot_1.qang * dt2 + angpot_1.pang * dt3 
			    + angpot_1.sang * dt4);

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
		    }

/*     increment the total bond angle bending energy */

		    energi_1.ea += e;
		}
	    }
	}
    }
    return 0;
} /* eangle_ */

#undef angtyp_ref
#undef iang_ref


