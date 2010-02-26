/* rgdstep.f -- translated by f2c (version 20050501).
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
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal vcm[3000]	/* was [3][1000] */, wcm[3000]	/* was [3][
	    1000] */, lm[3000]	/* was [3][1000] */, vc[3000]	/* was [3][
	    1000] */, wc[3000]	/* was [3][1000] */;
    logical linear[1000];
} rgddyn_;

#define rgddyn_1 rgddyn_

struct {
    doublereal vir[9]	/* was [3][3] */;
} virial_;

#define virial_1 virial_

/* Table of constant values */

static integer c__3 = 3;



/*     ########################################################### */
/*     ##                 COPYRIGHT (C) 2001 by                 ## */
/*     ##  Andrey Kutepov, Marina A. Vorobieva & Jay W. Ponder  ## */
/*     ##                  All Rights Reserved                  ## */
/*     ########################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine rgdstep  --  rigid body molecular dynamics step  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "rgdstep" performs a single molecular dynamics time step */
/*     for a rigid body calculation */

/*     literature reference: */

/*     W. Smith, "Hail Euler and Farewell: Rotational Motion in the */
/*     Laboratory Frame", CCP5 Newsletter, February 2005 */

/*     based on an original algorithm developed by Andrey Kutapov */
/*     and Marina A. Vorobieva, VNIITF, Russian Federal Nuclear */
/*     Center, Chelyabinsk, Russia, February 2001 */


/* Subroutine */ int rgdstep_(integer *istep, doublereal *dt)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 RGDSTEP  --  Angular Momentum Converge"
	    "nce\002,\002 Failure\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *), 
	    cholesky_(integer *, doublereal *, doublereal *), pressure_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer i__, j, k;
    static doublereal x2, y2, z2, fc[3], rc[3], tc[3], fx, fy, fz, xm[25000], 
	    ym[25000], zm[25000], xp[25000], yp[25000], xr, yr, zr, zp[25000],
	     dfi[3], eps, vcp[3], wcp[3], ekin[9]	/* was [3][3] */;
    static integer iter;
    static doublereal epot, temp, pres;
    static integer size;
    static doublereal etot, arot[9]	/* was [3][3] */;
    static integer stop;
    extern /* Subroutine */ int fatal_(void);
    static doublereal delta, weigh, rcold[3], inert[6], eksum;
    static integer start;
    static doublereal dfiold[3];
    extern /* Subroutine */ int mdsave_(integer *, doublereal *, doublereal *)
	    ;
    static doublereal derivs[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int rotrgd_(doublereal *, doublereal *), mdstat_(
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *), prterr_(void);
    static doublereal stress[9]	/* was [3][3] */;
    extern /* Subroutine */ int temper2_(doublereal *, doublereal *, 
	    doublereal *), kinetic_(doublereal *, doublereal *), linbody_(
	    integer *, doublereal *, doublereal *);
    static integer maxiter;

    /* Fortran I/O blocks */
    static cilist io___39 = { 0, 0, 0, fmt_10, 0 };



#define lm_ref(a_1,a_2) rgddyn_1.lm[(a_2)*3 + a_1 - 4]
#define vc_ref(a_1,a_2) rgddyn_1.vc[(a_2)*3 + a_1 - 4]
#define wc_ref(a_1,a_2) rgddyn_1.wc[(a_2)*3 + a_1 - 4]
#define vcm_ref(a_1,a_2) rgddyn_1.vcm[(a_2)*3 + a_1 - 4]
#define wcm_ref(a_1,a_2) rgddyn_1.wcm[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]
#define arot_ref(a_1,a_2) arot[(a_2)*3 + a_1 - 4]
#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1 - 4]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  rgddyn.i  --  velocities and momenta for rigid body MD  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     vcm     current translational velocity of each rigid body */
/*     wcm     current angular velocity of each rigid body */
/*     lm      current angular momentum of each rigid body */
/*     vc      half-step translational velocity for kinetic energy */
/*     wc      half-step angular velocity for kinetic energy */
/*     linear  logical flag to mark group as linear or nonlinear */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  units.i  --  physical constants and unit conversions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     literature reference: */

/*     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended */
/*     Values of the Fundamental Physical Constants: 2006", Reviews of */
/*     Modern Physics, 80, 633-730 (2008) */

/*     The "2006 CODATA Recommended Values" are also available from */
/*     the NIST Reference on Constants, Units, and Uncertainty site */
/*     at http://physics.nist.gov/cuu/index.html */

/*     Most values below are derived from the 2006 CODATA reference */

/*     The conversion from calorie to Joule is the definition of the */
/*     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992) */

/*     The "coulomb" energy conversion factor is found by dimensional */
/*     analysis of Coulomb's Law, ie, by dividing the square of the */
/*     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is */
/*     the permittivity of vacuum (the "electric constant"); note that */
/*     eps0 is typically given in F/m, equivalent to C**2/(J-m) */

/*     The approximate value used for the Debye, 3.33564 x 10-30 C-m, */
/*     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997) */

/*     The value of "prescon" is based on definition of 1 atmosphere */
/*     as 101325 Pa set by the 10th Conference Generale des Poids et */
/*     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3 */

/*     avogadro    Avogadro's number (N) in particles/mole */
/*     lightspd    speed of light in vacuum (c) in cm/ps */
/*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K */
/*     gasconst    ideal gas constant (R) in kcal/mole/K */
/*     emass       mass of an electron in atomic mass units */
/*     joule       conversion from calories to joules */
/*     convert     conversion from kcal to g*Ang**2/ps**2 */
/*     bohr        conversion from Bohrs to Angstroms */
/*     hartree     conversion from Hartree to kcal/mole */
/*     evolt       conversion from Hartree to electron-volts */
/*     efreq       conversion from Hartree to cm-1 */
/*     coulomb     conversion from electron**2/Ang to kcal/mole */
/*     debye       conversion from electron-Ang to Debyes */
/*     prescon     conversion from kcal/mole/Ang**3 to Atm */




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




/*     set iteration limit and tolerance for angular momenta */

    maxiter = 15;
    eps = 1e-12;

/*     get the energy and atomic forces prior to the step */

    gradient_(&epot, derivs);

/*     perform the integration step for each rigid body */

    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	start = igrp_ref(1, i__);
	stop = igrp_ref(2, i__);
	size = stop - start + 1;
	for (j = 1; j <= 3; ++j) {
	    rc[j - 1] = 0.;
	}
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    weigh = atmtyp_1.mass[k - 1];
	    rc[0] += atoms_1.x[k - 1] * weigh;
	    rc[1] += atoms_1.y[k - 1] * weigh;
	    rc[2] += atoms_1.z__[k - 1] * weigh;
	}
	for (j = 1; j <= 3; ++j) {
	    rc[j - 1] /= group_1.grpmass[i__ - 1];
	}

/*     find center of mass offsets only for first step */

	if (*istep == 1) {
	    i__2 = stop;
	    for (j = start; j <= i__2; ++j) {
		k = group_1.kgrp[j - 1];
		xm[k - 1] = atoms_1.x[k - 1] - rc[0];
		ym[k - 1] = atoms_1.y[k - 1] - rc[1];
		zm[k - 1] = atoms_1.z__[k - 1] - rc[2];
	    }
	}

/*     compute the force and torque components for rigid body */

	for (j = 1; j <= 3; ++j) {
	    fc[j - 1] = 0.;
	    tc[j - 1] = 0.;
	}
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    xr = atoms_1.x[k - 1] - rc[0];
	    yr = atoms_1.y[k - 1] - rc[1];
	    zr = atoms_1.z__[k - 1] - rc[2];
	    fx = derivs_ref(1, k) * -418.4;
	    fy = derivs_ref(2, k) * -418.4;
	    fz = derivs_ref(3, k) * -418.4;
	    fc[0] += fx;
	    fc[1] += fy;
	    fc[2] += fz;
	    tc[0] = tc[0] + yr * fz - zr * fy;
	    tc[1] = tc[1] + zr * fx - xr * fz;
	    tc[2] = tc[2] + xr * fy - yr * fx;
	}

/*     update the translational velocity of the center of mass */

	for (j = 1; j <= 3; ++j) {
	    vcp[j - 1] = vcm_ref(j, i__) + *dt * fc[j - 1] / group_1.grpmass[
		    i__ - 1];
	    vc_ref(j, i__) = (vcm_ref(j, i__) + vcp[j - 1]) * .5;
	    vcm_ref(j, i__) = vcp[j - 1];
	}

/*     update the coordinates of the group center of mass */

	for (j = 1; j <= 3; ++j) {
	    rcold[j - 1] = rc[j - 1];
	    rc[j - 1] += *dt * vcp[j - 1];
	}

/*     single atom groups are treated as a separate case */

	if (size == 1) {
	    k = group_1.kgrp[igrp_ref(1, i__) - 1];
	    atoms_1.x[k - 1] = rc[0];
	    atoms_1.y[k - 1] = rc[1];
	    atoms_1.z__[k - 1] = rc[2];
	    for (j = 1; j <= 3; ++j) {
		wcm_ref(j, i__) = 0.;
		lm_ref(j, i__) = 0.;
	    }

/*     get impulse moment in fixed space coordinate system */

	} else {
	    for (j = 1; j <= 3; ++j) {
		lm_ref(j, i__) = lm_ref(j, i__) + *dt * tc[j - 1];
		dfi[j - 1] = *dt * wcm_ref(j, i__);
		dfiold[j - 1] = dfi[j - 1];
	    }

/*     use iterative scheme to converge the angular momenta */

	    iter = 0;
	    delta = 1.;
	    while(delta > eps && iter < maxiter) {
		++iter;
		rotrgd_(dfi, arot);

/*     calculate the inertia tensor from rotated coordinates */

		for (j = 1; j <= 6; ++j) {
		    inert[j - 1] = 0.;
		}
		i__2 = stop;
		for (j = start; j <= i__2; ++j) {
		    k = group_1.kgrp[j - 1];
		    xr = arot_ref(1, 1) * xm[k - 1] + arot_ref(1, 2) * ym[k - 
			    1] + arot_ref(1, 3) * zm[k - 1];
		    yr = arot_ref(2, 1) * xm[k - 1] + arot_ref(2, 2) * ym[k - 
			    1] + arot_ref(2, 3) * zm[k - 1];
		    zr = arot_ref(3, 1) * xm[k - 1] + arot_ref(3, 2) * ym[k - 
			    1] + arot_ref(3, 3) * zm[k - 1];
		    x2 = xr * xr;
		    y2 = yr * yr;
		    z2 = zr * zr;
		    weigh = atmtyp_1.mass[k - 1];
		    inert[0] += weigh * (y2 + z2);
		    inert[1] -= weigh * xr * yr;
		    inert[2] -= weigh * xr * zr;
		    inert[3] += weigh * (x2 + z2);
		    inert[4] -= weigh * yr * zr;
		    inert[5] += weigh * (x2 + y2);
		    xp[k - 1] = xr;
		    yp[k - 1] = yr;
		    zp[k - 1] = zr;
		}

/*     compute the angular velocity from the relation L=Iw */

		for (j = 1; j <= 3; ++j) {
		    wcp[j - 1] = lm_ref(j, i__);
		}
		if (rgddyn_1.linear[i__ - 1]) {
		    linbody_(&i__, inert, wcp);
		} else {
		    cholesky_(&c__3, inert, wcp);
		}
		delta = 0.;
		for (j = 1; j <= 3; ++j) {
		    dfi[j - 1] = *dt * .5 * (wcm_ref(j, i__) + wcp[j - 1]);
		    delta += (d__1 = dfi[j - 1] - dfiold[j - 1], abs(d__1));
		    dfiold[j - 1] = dfi[j - 1];
		}
	    }

/*     check to make sure the angular momenta converged */

	    if (delta > eps) {
		io___39.ciunit = iounit_1.iout;
		s_wsfe(&io___39);
		e_wsfe();
		prterr_();
		fatal_();
	    }

/*     set the final angular velocity and atomic coordinates */

	    for (j = 1; j <= 3; ++j) {
		dfi[j - 1] = *dt * wcp[j - 1];
	    }
	    rotrgd_(dfi, arot);
	    i__2 = stop;
	    for (j = start; j <= i__2; ++j) {
		k = group_1.kgrp[j - 1];
		xr = atoms_1.x[k - 1] - rcold[0];
		yr = atoms_1.y[k - 1] - rcold[1];
		zr = atoms_1.z__[k - 1] - rcold[2];
		atoms_1.x[k - 1] = arot_ref(1, 1) * xr + arot_ref(1, 2) * yr 
			+ arot_ref(1, 3) * zr + rc[0];
		atoms_1.y[k - 1] = arot_ref(2, 1) * xr + arot_ref(2, 2) * yr 
			+ arot_ref(2, 3) * zr + rc[1];
		atoms_1.z__[k - 1] = arot_ref(3, 1) * xr + arot_ref(3, 2) * 
			yr + arot_ref(3, 3) * zr + rc[2];
	    }
	    for (j = 1; j <= 3; ++j) {
		wc_ref(j, i__) = (wcm_ref(j, i__) + wcp[j - 1]) * .5;
		wcm_ref(j, i__) = wcp[j - 1];
	    }
	}
    }

/*     update the distance to center of mass for each atom */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xm[i__ - 1] = xp[i__ - 1];
	ym[i__ - 1] = yp[i__ - 1];
	zm[i__ - 1] = zp[i__ - 1];
    }

/*     make center of mass correction to virial for rigid body */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vir_ref(1, 1) = vir_ref(1, 1) - xm[i__ - 1] * derivs_ref(1, i__);
	vir_ref(2, 1) = vir_ref(2, 1) - ym[i__ - 1] * derivs_ref(1, i__);
	vir_ref(3, 1) = vir_ref(3, 1) - zm[i__ - 1] * derivs_ref(1, i__);
	vir_ref(1, 2) = vir_ref(1, 2) - xm[i__ - 1] * derivs_ref(2, i__);
	vir_ref(2, 2) = vir_ref(2, 2) - ym[i__ - 1] * derivs_ref(2, i__);
	vir_ref(3, 2) = vir_ref(3, 2) - zm[i__ - 1] * derivs_ref(2, i__);
	vir_ref(1, 3) = vir_ref(1, 3) - xm[i__ - 1] * derivs_ref(3, i__);
	vir_ref(2, 3) = vir_ref(2, 3) - ym[i__ - 1] * derivs_ref(3, i__);
	vir_ref(3, 3) = vir_ref(3, 3) - zm[i__ - 1] * derivs_ref(3, i__);
    }

/*     accumulate the kinetic energy and its outer product */

    kinetic_(&eksum, ekin);

/*     compute and control the temperature and pressure */

    temper2_(dt, &eksum, &temp);
    pressure_(dt, &epot, ekin, &temp, &pres, stress);

/*     system energy is sum of kinetic and potential energies */

    etot = eksum + epot;

/*     compute statistics and save trajectory for this step */

    mdstat_(istep, dt, &etot, &epot, &eksum, &temp, &pres);
    mdsave_(istep, dt, &epot);
    return 0;
} /* rgdstep_ */

#undef derivs_ref
#undef arot_ref
#undef igrp_ref
#undef vir_ref
#undef wcm_ref
#undef vcm_ref
#undef wc_ref
#undef vc_ref
#undef lm_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine rotrgd  --  rigid dynamics rotation matrix  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "rotrgd" finds the rotation matrix for a rigid body due */
/*     to a single step of dynamics */


/* Subroutine */ int rotrgd_(doublereal *dfi, doublereal *arot)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal x, y, z__, xc, yc, zc, xs, ys, zs, sine, anorm, cosine, 
	    coterm;


#define arot_ref(a_1,a_2) arot[(a_2)*3 + a_1]



/*     construct rotation matrix from angular distance */

    /* Parameter adjustments */
    arot -= 4;
    --dfi;

    /* Function Body */
/* Computing 2nd power */
    d__1 = dfi[1];
/* Computing 2nd power */
    d__2 = dfi[2];
/* Computing 2nd power */
    d__3 = dfi[3];
    anorm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    cosine = cos(anorm);
    sine = sin(anorm);
    coterm = 1. - cosine;
    if (anorm <= 0.) {
	anorm = 1.;
    }
    x = dfi[1] / anorm;
    y = dfi[2] / anorm;
    z__ = dfi[3] / anorm;
    xc = x * coterm;
    yc = y * coterm;
    zc = z__ * coterm;
    xs = x * sine;
    ys = y * sine;
    zs = z__ * sine;
    arot_ref(1, 1) = xc * x + cosine;
    arot_ref(2, 1) = xc * y + zs;
    arot_ref(3, 1) = xc * z__ - ys;
    arot_ref(1, 2) = yc * x - zs;
    arot_ref(2, 2) = yc * y + cosine;
    arot_ref(3, 2) = yc * z__ + xs;
    arot_ref(1, 3) = zc * x + ys;
    arot_ref(2, 3) = zc * y - xs;
    arot_ref(3, 3) = zc * z__ + cosine;
    return 0;
} /* rotrgd_ */

#undef arot_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine linbody  --  angular velocity of linear body  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "linbody" finds the angular velocity of a linear rigid body */
/*     given the inertia tensor and angular momentum */


/* Subroutine */ int linbody_(integer *i__, doublereal *inert, doublereal *
	wcp)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal b1, b2, r1[3], r2[3], r3[3], w1, w2, a11, a12, a22, 
	    rmin, rmol[3], rinv;


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




/*     construct a normalized vector along the molecular axis */

    /* Parameter adjustments */
    --wcp;
    --inert;

    /* Function Body */
    j = group_1.kgrp[igrp_ref(1, *i__) - 1];
    k = group_1.kgrp[igrp_ref(2, *i__) - 1];
    rmol[0] = atoms_1.x[k - 1] - atoms_1.x[j - 1];
    rmol[1] = atoms_1.y[k - 1] - atoms_1.y[j - 1];
    rmol[2] = atoms_1.z__[k - 1] - atoms_1.z__[j - 1];
/* Computing 2nd power */
    d__1 = rmol[0];
/* Computing 2nd power */
    d__2 = rmol[1];
/* Computing 2nd power */
    d__3 = rmol[2];
    rinv = 1. / sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    for (j = 1; j <= 3; ++j) {
	rmol[j - 1] *= rinv;
    }

/*     find two orthogonal vectors to complete coordinate frame */

    k = 1;
    rmin = abs(rmol[0]);
    for (j = 2; j <= 3; ++j) {
	if ((d__1 = rmol[j - 1], abs(d__1)) < rmin) {
	    k = j;
	    rmin = (d__1 = rmol[j - 1], abs(d__1));
	}
    }
    for (j = 1; j <= 3; ++j) {
	r1[j - 1] = -rmol[k - 1] * rmol[j - 1];
    }
    r1[k - 1] += 1.;
/* Computing 2nd power */
    d__1 = r1[0];
/* Computing 2nd power */
    d__2 = r1[1];
/* Computing 2nd power */
    d__3 = r1[2];
    rinv = 1. / sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    for (j = 1; j <= 3; ++j) {
	r1[j - 1] *= rinv;
    }
    r2[0] = r1[1] * rmol[2] - r1[2] * rmol[1];
    r2[1] = r1[2] * rmol[0] - r1[0] * rmol[2];
    r2[2] = r1[0] * rmol[1] - r1[1] * rmol[0];

/*     solve the 2-by-2 linear system for angular velocity */

    r3[0] = inert[1] * r1[0] + inert[2] * r1[1] + inert[3] * r1[2];
    r3[1] = inert[2] * r1[0] + inert[4] * r1[1] + inert[5] * r1[2];
    r3[2] = inert[3] * r1[0] + inert[5] * r1[1] + inert[6] * r1[2];
    a11 = r1[0] * r3[0] + r1[1] * r3[1] + r1[2] * r3[2];
    r3[0] = inert[1] * r2[0] + inert[2] * r2[1] + inert[3] * r2[2];
    r3[1] = inert[2] * r2[0] + inert[4] * r2[1] + inert[5] * r2[2];
    r3[2] = inert[3] * r2[0] + inert[5] * r2[1] + inert[6] * r2[2];
    a12 = r1[0] * r3[0] + r1[1] * r3[1] + r1[2] * r3[2];
    a22 = r2[0] * r3[0] + r2[1] * r3[1] + r2[2] * r3[2];
    b1 = r1[0] * wcp[1] + r1[1] * wcp[2] + r1[2] * wcp[3];
    b2 = r2[0] * wcp[1] + r2[1] * wcp[2] + r2[2] * wcp[3];
    w1 = (a12 * b2 - a22 * b1) / (a12 * a12 - a11 * a22);
    w2 = (b2 - a12 * w1) / a22;
    for (j = 1; j <= 3; ++j) {
	wcp[j] = w1 * r1[j - 1] + w2 * r2[j - 1];
    }
    return 0;
} /* linbody_ */

#undef igrp_ref


