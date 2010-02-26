/* epitors2.f -- translated by f2c (version 20050501).
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
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

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
    doublereal kpit[100000];
    integer npitors, ipit[600000]	/* was [6][100000] */;
} pitors_;

#define pitors_1 pitors_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

struct {
    doublereal idihunit, itorunit, torsunit, ptorunit, storunit, ttorunit;
} torpot_;

#define torpot_1 torpot_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine epitors2  --  out-of-plane bend Hessian; numer  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "epitors2" calculates the second derivatives of the pi-orbital */
/*     torsion energy for a single atom using finite difference methods */


/* Subroutine */ int epitors2_(integer *i__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static logical twosided;
    static integer j;
    extern /* Subroutine */ int epitors2a_(integer *);
    static doublereal d0[75000]	/* was [3][25000] */;
    static integer ia, ib, ic, id, ie, ig;
    static doublereal old, eps, fgrp, term;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;
    static integer ipitors;


#define d0_ref(a_1,a_2) d0[(a_2)*3 + a_1 - 4]
#define dept_ref(a_1,a_2) deriv_1.dept[(a_2)*3 + a_1 - 4]
#define ipit_ref(a_1,a_2) pitors_1.ipit[(a_2)*6 + a_1 - 7]
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

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  pitors.i  --  pi-orbital torsions in the current structure  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     kpit      2-fold pi-orbital torsional force constants */
/*     npitors   total number of pi-orbital torsional interactions */
/*     ipit      numbers of the atoms in each pi-orbital torsion */




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




/*     set stepsize for derivatives and default group weight */

    eps = 1e-5;
    fgrp = 1.;
    twosided = FALSE_;
    if (atoms_1.n <= 50) {
	twosided = TRUE_;
    }

/*     compute numerical out-of-plane Hessian for current atom */

    i__1 = pitors_1.npitors;
    for (ipitors = 1; ipitors <= i__1; ++ipitors) {
	ia = ipit_ref(1, ipitors);
	ib = ipit_ref(2, ipitors);
	ic = ipit_ref(3, ipitors);
	id = ipit_ref(4, ipitors);
	ie = ipit_ref(5, ipitors);
	ig = ipit_ref(6, ipitors);

/*     decide whether to compute the current interaction */

	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &ie, &ig);
	}
	if (proceed) {
	    proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
		    usage_1.use[ic - 1] || usage_1.use[id - 1] || usage_1.use[
		    ie - 1] || usage_1.use[ig - 1];
	}
	if (proceed) {
	    term = fgrp / eps;

/*     find first derivatives for the base structure */

	    if (! twosided) {
		epitors2a_(&ipitors);
		for (j = 1; j <= 3; ++j) {
		    d0_ref(j, ia) = dept_ref(j, ia);
		    d0_ref(j, ib) = dept_ref(j, ib);
		    d0_ref(j, ic) = dept_ref(j, ic);
		    d0_ref(j, id) = dept_ref(j, id);
		    d0_ref(j, ie) = dept_ref(j, ie);
		    d0_ref(j, ig) = dept_ref(j, ig);
		}
	    }

/*     find numerical x-components via perturbed structures */

	    old = atoms_1.x[*i__ - 1];
	    if (twosided) {
		atoms_1.x[*i__ - 1] -= eps * .5;
		epitors2a_(&ipitors);
		for (j = 1; j <= 3; ++j) {
		    d0_ref(j, ia) = dept_ref(j, ia);
		    d0_ref(j, ib) = dept_ref(j, ib);
		    d0_ref(j, ic) = dept_ref(j, ic);
		    d0_ref(j, id) = dept_ref(j, id);
		    d0_ref(j, ie) = dept_ref(j, ie);
		    d0_ref(j, ig) = dept_ref(j, ig);
		}
	    }
	    atoms_1.x[*i__ - 1] += eps;
	    epitors2a_(&ipitors);
	    atoms_1.x[*i__ - 1] = old;
	    for (j = 1; j <= 3; ++j) {
		hessx_ref(j, ia) = hessx_ref(j, ia) + term * (dept_ref(j, ia) 
			- d0_ref(j, ia));
		hessx_ref(j, ib) = hessx_ref(j, ib) + term * (dept_ref(j, ib) 
			- d0_ref(j, ib));
		hessx_ref(j, ic) = hessx_ref(j, ic) + term * (dept_ref(j, ic) 
			- d0_ref(j, ic));
		hessx_ref(j, id) = hessx_ref(j, id) + term * (dept_ref(j, id) 
			- d0_ref(j, id));
		hessx_ref(j, ie) = hessx_ref(j, ie) + term * (dept_ref(j, ie) 
			- d0_ref(j, ie));
		hessx_ref(j, ig) = hessx_ref(j, ig) + term * (dept_ref(j, ig) 
			- d0_ref(j, ig));
	    }

/*     find numerical y-components via perturbed structures */

	    old = atoms_1.y[*i__ - 1];
	    if (twosided) {
		atoms_1.y[*i__ - 1] -= eps * .5;
		epitors2a_(&ipitors);
		for (j = 1; j <= 3; ++j) {
		    d0_ref(j, ia) = dept_ref(j, ia);
		    d0_ref(j, ib) = dept_ref(j, ib);
		    d0_ref(j, ic) = dept_ref(j, ic);
		    d0_ref(j, id) = dept_ref(j, id);
		    d0_ref(j, ie) = dept_ref(j, ie);
		    d0_ref(j, ig) = dept_ref(j, ig);
		}
	    }
	    atoms_1.y[*i__ - 1] += eps;
	    epitors2a_(&ipitors);
	    atoms_1.y[*i__ - 1] = old;
	    for (j = 1; j <= 3; ++j) {
		hessy_ref(j, ia) = hessy_ref(j, ia) + term * (dept_ref(j, ia) 
			- d0_ref(j, ia));
		hessy_ref(j, ib) = hessy_ref(j, ib) + term * (dept_ref(j, ib) 
			- d0_ref(j, ib));
		hessy_ref(j, ic) = hessy_ref(j, ic) + term * (dept_ref(j, ic) 
			- d0_ref(j, ic));
		hessy_ref(j, id) = hessy_ref(j, id) + term * (dept_ref(j, id) 
			- d0_ref(j, id));
		hessy_ref(j, ie) = hessy_ref(j, ie) + term * (dept_ref(j, ie) 
			- d0_ref(j, ie));
		hessy_ref(j, ig) = hessy_ref(j, ig) + term * (dept_ref(j, ig) 
			- d0_ref(j, ig));
	    }

/*     find numerical z-components via perturbed structures */

	    old = atoms_1.z__[*i__ - 1];
	    if (twosided) {
		atoms_1.z__[*i__ - 1] -= eps * .5;
		epitors2a_(&ipitors);
		for (j = 1; j <= 3; ++j) {
		    d0_ref(j, ia) = dept_ref(j, ia);
		    d0_ref(j, ib) = dept_ref(j, ib);
		    d0_ref(j, ic) = dept_ref(j, ic);
		    d0_ref(j, id) = dept_ref(j, id);
		    d0_ref(j, ie) = dept_ref(j, ie);
		    d0_ref(j, ig) = dept_ref(j, ig);
		}
	    }
	    atoms_1.z__[*i__ - 1] += eps;
	    epitors2a_(&ipitors);
	    atoms_1.z__[*i__ - 1] = old;
	    for (j = 1; j <= 3; ++j) {
		hessz_ref(j, ia) = hessz_ref(j, ia) + term * (dept_ref(j, ia) 
			- d0_ref(j, ia));
		hessz_ref(j, ib) = hessz_ref(j, ib) + term * (dept_ref(j, ib) 
			- d0_ref(j, ib));
		hessz_ref(j, ic) = hessz_ref(j, ic) + term * (dept_ref(j, ic) 
			- d0_ref(j, ic));
		hessz_ref(j, id) = hessz_ref(j, id) + term * (dept_ref(j, id) 
			- d0_ref(j, id));
		hessz_ref(j, ie) = hessz_ref(j, ie) + term * (dept_ref(j, ie) 
			- d0_ref(j, ie));
		hessz_ref(j, ig) = hessz_ref(j, ig) + term * (dept_ref(j, ig) 
			- d0_ref(j, ig));
	    }
	}
    }
    return 0;
} /* epitors2_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef ipit_ref
#undef dept_ref
#undef d0_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine epitors2a  --  pi-orbital torsion derivatives  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "epitors2a" calculates the pi-orbital torsion first derivatives; */
/*     used in computation of finite difference second derivatives */


/* Subroutine */ int epitors2a_(integer *i__)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c2, s2, v2;
    static integer ia, ib, ic, id, ie, ig;
    static doublereal xt, yt, zt, xu, yu, zu, rt2, ru2, rdc, xad, yad, xia, 
	    yia, zia, xib, yib, zib, xic, yic, zic, xid, yid, zid, xie, yie, 
	    zie, xig, yig, zig, xip, yip, zip, xiq, yiq, ziq, zad, xbd, ybd, 
	    zbd, xec, yec, zec, xtu, ytu, ztu, xgc, ygc, zgc, xcp, ycp, zcp, 
	    xdc, ydc, zdc, xqd, yqd, zqd, xdp, ydp, zdp, xqc, yqc, zqc, sine, 
	    rtru, dphi2, sine2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal dedxt, dedyt, dedzt, dedxu, dedyu, dedzu, dedphi, 
	    dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, dedxic, dedyic, 
	    dedzic, dedxid, dedyid, dedzid, dedxie, dedyie, dedzie, cosine, 
	    dedxig, dedyig, dedzig, dedxip, dedyip, dedzip, dedxiq, dedyiq, 
	    dedziq, cosine2;


#define dept_ref(a_1,a_2) deriv_1.dept[(a_2)*3 + a_1 - 4]
#define ipit_ref(a_1,a_2) pitors_1.ipit[(a_2)*6 + a_1 - 7]



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
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  pitors.i  --  pi-orbital torsions in the current structure  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     kpit      2-fold pi-orbital torsional force constants */
/*     npitors   total number of pi-orbital torsional interactions */
/*     ipit      numbers of the atoms in each pi-orbital torsion */




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




/*     set the atom numbers for this pi-orbital torsion */

    ia = ipit_ref(1, *i__);
    ib = ipit_ref(2, *i__);
    ic = ipit_ref(3, *i__);
    id = ipit_ref(4, *i__);
    ie = ipit_ref(5, *i__);
    ig = ipit_ref(6, *i__);

/*     compute the value of the pi-orbital torsion angle */

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
    xig = atoms_1.x[ig - 1];
    yig = atoms_1.y[ig - 1];
    zig = atoms_1.z__[ig - 1];
    xad = xia - xid;
    yad = yia - yid;
    zad = zia - zid;
    xbd = xib - xid;
    ybd = yib - yid;
    zbd = zib - zid;
    xec = xie - xic;
    yec = yie - yic;
    zec = zie - zic;
    xgc = xig - xic;
    ygc = yig - yic;
    zgc = zig - zic;
    if (bound_1.use_polymer__) {
	image_(&xad, &yad, &zad);
	image_(&xbd, &ybd, &zbd);
	image_(&xec, &yec, &zec);
	image_(&xgc, &ygc, &zgc);
    }
    xip = yad * zbd - ybd * zad + xic;
    yip = zad * xbd - zbd * xad + yic;
    zip = xad * ybd - xbd * yad + zic;
    xiq = yec * zgc - ygc * zec + xid;
    yiq = zec * xgc - zgc * xec + yid;
    ziq = xec * ygc - xgc * yec + zid;
    xcp = xic - xip;
    ycp = yic - yip;
    zcp = zic - zip;
    xdc = xid - xic;
    ydc = yid - yic;
    zdc = zid - zic;
    xqd = xiq - xid;
    yqd = yiq - yid;
    zqd = ziq - zid;
    if (bound_1.use_polymer__) {
	image_(&xcp, &ycp, &zcp);
	image_(&xdc, &ydc, &zdc);
	image_(&xqd, &yqd, &zqd);
    }
    xt = ycp * zdc - ydc * zcp;
    yt = zcp * xdc - zdc * xcp;
    zt = xcp * ydc - xdc * ycp;
    xu = ydc * zqd - yqd * zdc;
    yu = zdc * xqd - zqd * xdc;
    zu = xdc * yqd - xqd * ydc;
    xtu = yt * zu - yu * zt;
    ytu = zt * xu - zu * xt;
    ztu = xt * yu - xu * yt;
    rt2 = xt * xt + yt * yt + zt * zt;
    ru2 = xu * xu + yu * yu + zu * zu;
    rtru = sqrt(rt2 * ru2);
    if (rtru != 0.) {
	rdc = sqrt(xdc * xdc + ydc * ydc + zdc * zdc);
	cosine = (xt * xu + yt * yu + zt * zu) / rtru;
	sine = (xdc * xtu + ydc * ytu + zdc * ztu) / (rdc * rtru);

/*     set the pi-orbital torsion parameters for this angle */

	v2 = pitors_1.kpit[*i__ - 1];
	c2 = -1.;
	s2 = 0.;

/*     compute the multiple angle trigonometry and the phase terms */

	cosine2 = cosine * cosine - sine * sine;
	sine2 = cosine * 2. * sine;
	dphi2 = (cosine2 * s2 - sine2 * c2) * 2.;

/*     calculate pi-orbital torsion energy and master chain rule term */

	dedphi = torpot_1.ptorunit * v2 * dphi2;

/*     chain rule terms for first derivative components */

	xdp = xid - xip;
	ydp = yid - yip;
	zdp = zid - zip;
	xqc = xiq - xic;
	yqc = yiq - yic;
	zqc = ziq - zic;
	dedxt = dedphi * (yt * zdc - ydc * zt) / (rt2 * rdc);
	dedyt = dedphi * (zt * xdc - zdc * xt) / (rt2 * rdc);
	dedzt = dedphi * (xt * ydc - xdc * yt) / (rt2 * rdc);
	dedxu = -dedphi * (yu * zdc - ydc * zu) / (ru2 * rdc);
	dedyu = -dedphi * (zu * xdc - zdc * xu) / (ru2 * rdc);
	dedzu = -dedphi * (xu * ydc - xdc * yu) / (ru2 * rdc);

/*     compute first derivative components for pi-orbital angle */

	dedxip = zdc * dedyt - ydc * dedzt;
	dedyip = xdc * dedzt - zdc * dedxt;
	dedzip = ydc * dedxt - xdc * dedyt;
	dedxic = ydp * dedzt - zdp * dedyt + zqd * dedyu - yqd * dedzu;
	dedyic = zdp * dedxt - xdp * dedzt + xqd * dedzu - zqd * dedxu;
	dedzic = xdp * dedyt - ydp * dedxt + yqd * dedxu - xqd * dedyu;
	dedxid = zcp * dedyt - ycp * dedzt + yqc * dedzu - zqc * dedyu;
	dedyid = xcp * dedzt - zcp * dedxt + zqc * dedxu - xqc * dedzu;
	dedzid = ycp * dedxt - xcp * dedyt + xqc * dedyu - yqc * dedxu;
	dedxiq = zdc * dedyu - ydc * dedzu;
	dedyiq = xdc * dedzu - zdc * dedxu;
	dedziq = ydc * dedxu - xdc * dedyu;

/*     compute first derivative components for individual atoms */

	dedxia = ybd * dedzip - zbd * dedyip;
	dedyia = zbd * dedxip - xbd * dedzip;
	dedzia = xbd * dedyip - ybd * dedxip;
	dedxib = zad * dedyip - yad * dedzip;
	dedyib = xad * dedzip - zad * dedxip;
	dedzib = yad * dedxip - xad * dedyip;
	dedxie = ygc * dedziq - zgc * dedyiq;
	dedyie = zgc * dedxiq - xgc * dedziq;
	dedzie = xgc * dedyiq - ygc * dedxiq;
	dedxig = zec * dedyiq - yec * dedziq;
	dedyig = xec * dedziq - zec * dedxiq;
	dedzig = yec * dedxiq - xec * dedyiq;
	dedxic = dedxic + dedxip - dedxie - dedxig;
	dedyic = dedyic + dedyip - dedyie - dedyig;
	dedzic = dedzic + dedzip - dedzie - dedzig;
	dedxid = dedxid + dedxiq - dedxia - dedxib;
	dedyid = dedyid + dedyiq - dedyia - dedyib;
	dedzid = dedzid + dedziq - dedzia - dedzib;

/*     increment the total pi-orbital torsion energy and gradient */

	dept_ref(1, ia) = dedxia;
	dept_ref(2, ia) = dedyia;
	dept_ref(3, ia) = dedzia;
	dept_ref(1, ib) = dedxib;
	dept_ref(2, ib) = dedyib;
	dept_ref(3, ib) = dedzib;
	dept_ref(1, ic) = dedxic;
	dept_ref(2, ic) = dedyic;
	dept_ref(3, ic) = dedzic;
	dept_ref(1, id) = dedxid;
	dept_ref(2, id) = dedyid;
	dept_ref(3, id) = dedzid;
	dept_ref(1, ie) = dedxie;
	dept_ref(2, ie) = dedyie;
	dept_ref(3, ie) = dedzie;
	dept_ref(1, ig) = dedxig;
	dept_ref(2, ig) = dedyig;
	dept_ref(3, ig) = dedzig;
    }
    return 0;
} /* epitors2a_ */

#undef ipit_ref
#undef dept_ref


