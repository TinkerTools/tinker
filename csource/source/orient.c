/* orient.f -- translated by f2c (version 20050501).
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
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    doublereal xrb[25000], yrb[25000], zrb[25000], rbc[6000]	/* was [6][
	    1000] */;
    logical use_rigid__;
} rigid_;

#define rigid_1 rigid_

struct {
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

/* Table of constant values */

static integer c__3 = 3;
static doublereal c_b7 = 3.141592653589793238;



/*     ############################################################## */
/*     ##  COPYRIGHT (C) 1997 by Rohit Pappu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine orient  --  rigid body reference coordinates  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "orient" computes a set of reference Cartesian coordinates */
/*     in standard orientation for each rigid body atom group */


/* Subroutine */ int orient_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    extern /* Subroutine */ int xyzrigid_(void);
    static doublereal a[9]	/* was [3][3] */;
    static integer i__, j, k;
    static doublereal phi, xcm, ycm, zcm, psi, cphi, cpsi;
    static integer init;
    static doublereal sphi, spsi;
    static integer stop;
    static doublereal theta, xterm, yterm, zterm, ctheta, stheta;


#define a_ref(a_1,a_2) a[(a_2)*3 + a_1 - 4]
#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]



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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  rigid.i  --  rigid body coordinates for atom groups  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     xrb         rigid body reference x-coordinate for each atom */
/*     yrb         rigid body reference y-coordinate for each atom */
/*     zrb         rigid body reference z-coordinate for each atom */
/*     rbc         current rigid body coordinates for each group */
/*     use_rigid   flag to mark use of rigid body coordinate system */




/*     use current coordinates as default reference coordinates */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rigid_1.xrb[i__ - 1] = atoms_1.x[i__ - 1];
	rigid_1.yrb[i__ - 1] = atoms_1.y[i__ - 1];
	rigid_1.zrb[i__ - 1] = atoms_1.z__[i__ - 1];
    }

/*     compute the rigid body coordinates for each atom group */

    xyzrigid_();

/*     get the center of mass and Euler angles for each group */

    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xcm = rbc_ref(1, i__);
	ycm = rbc_ref(2, i__);
	zcm = rbc_ref(3, i__);
	phi = rbc_ref(4, i__);
	theta = rbc_ref(5, i__);
	psi = rbc_ref(6, i__);
	cphi = cos(phi);
	sphi = sin(phi);
	ctheta = cos(theta);
	stheta = sin(theta);
	cpsi = cos(psi);
	spsi = sin(psi);

/*     construct the rotation matrix from Euler angle values */

	a_ref(1, 1) = ctheta * cphi;
	a_ref(2, 1) = spsi * stheta * cphi - cpsi * sphi;
	a_ref(3, 1) = cpsi * stheta * cphi + spsi * sphi;
	a_ref(1, 2) = ctheta * sphi;
	a_ref(2, 2) = spsi * stheta * sphi + cpsi * cphi;
	a_ref(3, 2) = cpsi * stheta * sphi - spsi * cphi;
	a_ref(1, 3) = -stheta;
	a_ref(2, 3) = ctheta * spsi;
	a_ref(3, 3) = ctheta * cpsi;

/*     translate and rotate each atom group into inertial frame */

	init = igrp_ref(1, i__);
	stop = igrp_ref(2, i__);
	i__2 = stop;
	for (j = init; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    xterm = atoms_1.x[k - 1] - xcm;
	    yterm = atoms_1.y[k - 1] - ycm;
	    zterm = atoms_1.z__[k - 1] - zcm;
	    rigid_1.xrb[k - 1] = a_ref(1, 1) * xterm + a_ref(1, 2) * yterm + 
		    a_ref(1, 3) * zterm;
	    rigid_1.yrb[k - 1] = a_ref(2, 1) * xterm + a_ref(2, 2) * yterm + 
		    a_ref(2, 3) * zterm;
	    rigid_1.zrb[k - 1] = a_ref(3, 1) * xterm + a_ref(3, 2) * yterm + 
		    a_ref(3, 3) * zterm;
	}
    }
    return 0;
} /* orient_ */

#undef igrp_ref
#undef rbc_ref
#undef a_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine xyzrigid  --  determine rigid body coordinates  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "xyzrigid" computes the center of mass and Euler angle rigid */
/*     body coordinates for each atom group in the system */

/*     literature reference: */

/*     Herbert Goldstein, "Classical Mechanics, 2nd Edition", */
/*     Addison-Wesley, Reading, MA, 1980; see the Euler angle */
/*     xyz convention in Appendix B */


/* Subroutine */ int xyzrigid_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int roteuler_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal a[9]	/* was [3][3] */;
    static integer i__, j, k, m;
    static doublereal xx, xy, xz, yy, yz, zz, vec[9]	/* was [3][3] */, phi,
	     dot, xcm, ycm, zcm, psi;
    static integer init, stop;
    static doublereal work1[3], work2[3], weigh, theta, total, xterm, yterm, 
	    zterm;
    extern /* Subroutine */ int jacobi_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal moment[3], tensor[9]	/* was [3][3] */;


#define a_ref(a_1,a_2) a[(a_2)*3 + a_1 - 4]
#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]
#define vec_ref(a_1,a_2) vec[(a_2)*3 + a_1 - 4]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]
#define tensor_ref(a_1,a_2) tensor[(a_2)*3 + a_1 - 4]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  atmtyp.i  --  atomic properties for each current atom  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     mass      atomic weight for each atom in the system */
/*     tag       integer atom labels from input coordinates file */
/*     class     atom class number for each atom in the system */
/*     atomic    atomic number for each atom in the system */
/*     valence   valence number for each atom in the system */
/*     name      atom name for each atom in the system */
/*     story     descriptive type for each atom in system */




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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  rigid.i  --  rigid body coordinates for atom groups  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     xrb         rigid body reference x-coordinate for each atom */
/*     yrb         rigid body reference y-coordinate for each atom */
/*     zrb         rigid body reference z-coordinate for each atom */
/*     rbc         current rigid body coordinates for each group */
/*     use_rigid   flag to mark use of rigid body coordinate system */




/*     get the first and last atom in the current group */

    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	init = igrp_ref(1, i__);
	stop = igrp_ref(2, i__);

/*     compute the position of the group center of mass */

	total = 0.;
	xcm = 0.;
	ycm = 0.;
	zcm = 0.;
	i__2 = stop;
	for (j = init; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    weigh = atmtyp_1.mass[k - 1];
	    total += weigh;
	    xcm += atoms_1.x[k - 1] * weigh;
	    ycm += atoms_1.y[k - 1] * weigh;
	    zcm += atoms_1.z__[k - 1] * weigh;
	}
	if (total != 0.) {
	    xcm /= total;
	    ycm /= total;
	    zcm /= total;
	}

/*     compute and then diagonalize the inertia tensor */

	xx = 0.;
	xy = 0.;
	xz = 0.;
	yy = 0.;
	yz = 0.;
	zz = 0.;
	i__2 = stop;
	for (j = init; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    weigh = atmtyp_1.mass[k - 1];
	    xterm = atoms_1.x[k - 1] - xcm;
	    yterm = atoms_1.y[k - 1] - ycm;
	    zterm = atoms_1.z__[k - 1] - zcm;
	    xx += xterm * xterm * weigh;
	    xy += xterm * yterm * weigh;
	    xz += xterm * zterm * weigh;
	    yy += yterm * yterm * weigh;
	    yz += yterm * zterm * weigh;
	    zz += zterm * zterm * weigh;
	}
	tensor_ref(1, 1) = yy + zz;
	tensor_ref(2, 1) = -xy;
	tensor_ref(3, 1) = -xz;
	tensor_ref(1, 2) = -xy;
	tensor_ref(2, 2) = xx + zz;
	tensor_ref(3, 2) = -yz;
	tensor_ref(1, 3) = -xz;
	tensor_ref(2, 3) = -yz;
	tensor_ref(3, 3) = xx + yy;
	jacobi_(&c__3, &c__3, tensor, moment, vec, work1, work2);

/*     select the direction for each principle moment axis */

	for (m = 1; m <= 2; ++m) {
	    i__2 = stop;
	    for (j = init; j <= i__2; ++j) {
		k = group_1.kgrp[j - 1];
		xterm = vec_ref(1, m) * (atoms_1.x[k - 1] - xcm);
		yterm = vec_ref(2, m) * (atoms_1.y[k - 1] - ycm);
		zterm = vec_ref(3, m) * (atoms_1.z__[k - 1] - zcm);
		dot = xterm + yterm + zterm;
		if (dot < 0.) {
		    vec_ref(1, m) = -vec_ref(1, m);
		    vec_ref(2, m) = -vec_ref(2, m);
		    vec_ref(3, m) = -vec_ref(3, m);
		}
		if (dot != 0.) {
		    goto L10;
		}
	    }
L10:
	    ;
	}

/*     moment axes must give a right-handed coordinate system */

	xterm = vec_ref(1, 1) * (vec_ref(2, 2) * vec_ref(3, 3) - vec_ref(2, 3)
		 * vec_ref(3, 2));
	yterm = vec_ref(2, 1) * (vec_ref(1, 3) * vec_ref(3, 2) - vec_ref(1, 2)
		 * vec_ref(3, 3));
	zterm = vec_ref(3, 1) * (vec_ref(1, 2) * vec_ref(2, 3) - vec_ref(1, 3)
		 * vec_ref(2, 2));
	dot = xterm + yterm + zterm;
	if (dot < 0.) {
	    for (j = 1; j <= 3; ++j) {
		vec_ref(j, 3) = -vec_ref(j, 3);
	    }
	}

/*     principal moment axes form rows of Euler rotation matrix */

	for (k = 1; k <= 3; ++k) {
	    for (j = 1; j <= 3; ++j) {
		a_ref(k, j) = vec_ref(j, k);
	    }
	}

/*     compute Euler angles consistent with the rotation matrix */

	roteuler_(a, &phi, &theta, &psi);

/*     set the rigid body coordinates for each atom group */

	rbc_ref(1, i__) = xcm;
	rbc_ref(2, i__) = ycm;
	rbc_ref(3, i__) = zcm;
	rbc_ref(4, i__) = phi;
	rbc_ref(5, i__) = theta;
	rbc_ref(6, i__) = psi;
    }
    return 0;
} /* xyzrigid_ */

#undef tensor_ref
#undef igrp_ref
#undef vec_ref
#undef rbc_ref
#undef a_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine roteuler  --  rotation matrix to Euler angles   ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "roteuler" computes a set of Euler angle values consistent */
/*     with an input rotation matrix */


/* Subroutine */ int roteuler_(doublereal *a, doublereal *phi, doublereal *
	theta, doublereal *psi)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double asin(doublereal), cos(doublereal), acos(doublereal), atan(
	    doublereal), sin(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal b[3];
    static integer i__;
    static doublereal eps, cphi;
    static logical flip[3];
    static doublereal cpsi, sphi, spsi, ctheta, stheta;


#define a_ref(a_1,a_2) a[(a_2)*3 + a_1]



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




/*     set the tolerance for Euler angles and rotation elements */

    /* Parameter adjustments */
    a -= 4;

    /* Function Body */
    eps = 1e-8;

/*     get a trial value of theta from a single rotation element */

/* Computing MIN */
/* Computing MAX */
    d__3 = -1., d__4 = -a_ref(1, 3);
    d__1 = 1., d__2 = max(d__3,d__4);
    *theta = asin((min(d__1,d__2)));
    ctheta = cos(*theta);
    stheta = -a_ref(1, 3);

/*     set the phi/psi difference when theta is either 90 or -90 */

    if (abs(ctheta) <= eps) {
	*phi = 0.;
	if ((d__1 = a_ref(3, 1), abs(d__1)) < eps) {
/* Computing MIN */
/* Computing MAX */
	    d__3 = -1., d__4 = -a_ref(2, 1) / a_ref(1, 3);
	    d__1 = 1., d__2 = max(d__3,d__4);
	    *psi = asin((min(d__1,d__2)));
	} else if ((d__1 = a_ref(2, 1), abs(d__1)) < eps) {
/* Computing MIN */
/* Computing MAX */
	    d__3 = -1., d__4 = -a_ref(3, 1) / a_ref(1, 3);
	    d__1 = 1., d__2 = max(d__3,d__4);
	    *psi = acos((min(d__1,d__2)));
	} else {
	    *psi = atan(a_ref(2, 1) / a_ref(3, 1));
	}

/*     set the phi and psi values for all other theta values */

    } else {
	if ((d__1 = a_ref(1, 1), abs(d__1)) < eps) {
/* Computing MIN */
/* Computing MAX */
	    d__3 = -1., d__4 = a_ref(1, 2) / ctheta;
	    d__1 = 1., d__2 = max(d__3,d__4);
	    *phi = asin((min(d__1,d__2)));
	} else if ((d__1 = a_ref(1, 2), abs(d__1)) < eps) {
/* Computing MIN */
/* Computing MAX */
	    d__3 = -1., d__4 = a_ref(1, 1) / ctheta;
	    d__1 = 1., d__2 = max(d__3,d__4);
	    *phi = acos((min(d__1,d__2)));
	} else {
	    *phi = atan(a_ref(1, 2) / a_ref(1, 1));
	}
	if ((d__1 = a_ref(3, 3), abs(d__1)) < eps) {
/* Computing MIN */
/* Computing MAX */
	    d__3 = -1., d__4 = a_ref(2, 3) / ctheta;
	    d__1 = 1., d__2 = max(d__3,d__4);
	    *psi = asin((min(d__1,d__2)));
	} else if ((d__1 = a_ref(2, 3), abs(d__1)) < eps) {
/* Computing MIN */
/* Computing MAX */
	    d__3 = -1., d__4 = a_ref(3, 3) / ctheta;
	    d__1 = 1., d__2 = max(d__3,d__4);
	    *psi = acos((min(d__1,d__2)));
	} else {
	    *psi = atan(a_ref(2, 3) / a_ref(3, 3));
	}
    }

/*     find sine and cosine of the trial phi and psi values */

    cphi = cos(*phi);
    sphi = sin(*phi);
    cpsi = cos(*psi);
    spsi = sin(*psi);

/*     reconstruct the diagonal of the rotation matrix */

    b[0] = ctheta * cphi;
    b[1] = spsi * stheta * sphi + cpsi * cphi;
    b[2] = ctheta * cpsi;

/*     compare the correct matrix diagonal to rebuilt diagonal */

    for (i__ = 1; i__ <= 3; ++i__) {
	flip[i__ - 1] = FALSE_;
	if ((d__1 = a_ref(i__, i__) - b[i__ - 1], abs(d__1)) > eps) {
	    flip[i__ - 1] = TRUE_;
	}
    }

/*     alter Euler angles to get correct rotation matrix values */

    if (flip[0] && flip[1]) {
	*phi -= d_sign(&c_b7, phi);
    }
    if (flip[0] && flip[2]) {
	*theta = -(*theta) + d_sign(&c_b7, theta);
    }
    if (flip[1] && flip[2]) {
	*psi -= d_sign(&c_b7, psi);
    }
    return 0;
} /* roteuler_ */

#undef a_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine rigidxyz  --  rigid body to Cartesian coords  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "rigidxyz" computes Cartesian coordinates for a rigid body */
/*     group via rotation and translation of reference coordinates */

/*     literature reference: */

/*     Herbert Goldstein, "Classical Mechanics, 2nd Edition", */
/*     Addison-Wesley, Reading, MA, 1980; see the Euler angle */
/*     xyz convention in Appendix B */


/* Subroutine */ int rigidxyz_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal a[9]	/* was [3][3] */;
    static integer i__, j, k;
    static doublereal phi, xcm, ycm, zcm, psi, cphi, cpsi;
    static integer init;
    static doublereal sphi, spsi;
    static integer stop;
    static doublereal theta, xterm, yterm, zterm, ctheta, stheta;


#define a_ref(a_1,a_2) a[(a_2)*3 + a_1 - 4]
#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]



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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  rigid.i  --  rigid body coordinates for atom groups  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     xrb         rigid body reference x-coordinate for each atom */
/*     yrb         rigid body reference y-coordinate for each atom */
/*     zrb         rigid body reference z-coordinate for each atom */
/*     rbc         current rigid body coordinates for each group */
/*     use_rigid   flag to mark use of rigid body coordinate system */




/*     get the center of mass and Euler angles for each group */

    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xcm = rbc_ref(1, i__);
	ycm = rbc_ref(2, i__);
	zcm = rbc_ref(3, i__);
	phi = rbc_ref(4, i__);
	theta = rbc_ref(5, i__);
	psi = rbc_ref(6, i__);
	cphi = cos(phi);
	sphi = sin(phi);
	ctheta = cos(theta);
	stheta = sin(theta);
	cpsi = cos(psi);
	spsi = sin(psi);

/*     construct the rotation matrix from Euler angle values */

	a_ref(1, 1) = ctheta * cphi;
	a_ref(2, 1) = spsi * stheta * cphi - cpsi * sphi;
	a_ref(3, 1) = cpsi * stheta * cphi + spsi * sphi;
	a_ref(1, 2) = ctheta * sphi;
	a_ref(2, 2) = spsi * stheta * sphi + cpsi * cphi;
	a_ref(3, 2) = cpsi * stheta * sphi - spsi * cphi;
	a_ref(1, 3) = -stheta;
	a_ref(2, 3) = ctheta * spsi;
	a_ref(3, 3) = ctheta * cpsi;

/*     rotate and translate reference coordinates into global frame */

	init = igrp_ref(1, i__);
	stop = igrp_ref(2, i__);
	i__2 = stop;
	for (j = init; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    xterm = rigid_1.xrb[k - 1];
	    yterm = rigid_1.yrb[k - 1];
	    zterm = rigid_1.zrb[k - 1];
	    atoms_1.x[k - 1] = a_ref(1, 1) * xterm + a_ref(2, 1) * yterm + 
		    a_ref(3, 1) * zterm + xcm;
	    atoms_1.y[k - 1] = a_ref(1, 2) * xterm + a_ref(2, 2) * yterm + 
		    a_ref(3, 2) * zterm + ycm;
	    atoms_1.z__[k - 1] = a_ref(1, 3) * xterm + a_ref(2, 3) * yterm + 
		    a_ref(3, 3) * zterm + zcm;
	}
    }
    return 0;
} /* rigidxyz_ */

#undef igrp_ref
#undef rbc_ref
#undef a_ref


