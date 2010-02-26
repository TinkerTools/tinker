/* mdrest.f -- translated by f2c (version 20050501).
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
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nfree;
    logical velsave, frcsave, uindsave;
    char integrate[10];
} mdstuf_;

#define mdstuf_1 mdstuf_

struct {
    doublereal v[75000]	/* was [3][25000] */, a[75000]	/* was [3][25000] */, 
	    aold[75000]	/* was [3][25000] */;
} moldyn_;

#define moldyn_1 moldyn_

struct {
    doublereal vcm[3000]	/* was [3][1000] */, wcm[3000]	/* was [3][
	    1000] */, lm[3000]	/* was [3][1000] */, vc[3000]	/* was [3][
	    1000] */, wc[3000]	/* was [3][1000] */;
    logical linear[1000];
} rgddyn_;

#define rgddyn_1 rgddyn_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine mdrest  --  stop system translation & rotation  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "mdrest" finds and removes any translational or rotational */
/*     kinetic energy of the overall system center of mass */


/* Subroutine */ int mdrest_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 System Linear Velocity :  \002,3d12.2,/"
	    ",\002 Translational Kinetic Energy :\002,10x,f12.4,\002 Kcal/mole"
	    "\002)";
    static char fmt_20[] = "(/,\002 System Angular Velocity : \002,3d12.2,/"
	    ",\002 Rotational Kinetic Energy :\002,13x,f12.4,\002 Kcal/mol"
	    "e\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), do_fio(
	    integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k;
    static doublereal xx, xy, yy, xz, zz, yz, eps, xcm[1000], ycm[1000], zcm[
	    1000], mang[3], vang[3], xdel, ydel, zdel, erot, vtot[3], xtot, 
	    ytot, ztot, weigh, etrans;
    extern /* Subroutine */ int invert_(integer *, integer *, doublereal *);
    static doublereal tensor[9]	/* was [3][3] */, totmass;

    /* Fortran I/O blocks */
    static cilist io___28 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_20, 0 };



#define v_ref(a_1,a_2) moldyn_1.v[(a_2)*3 + a_1 - 4]
#define vcm_ref(a_1,a_2) rgddyn_1.vcm[(a_2)*3 + a_1 - 4]
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
/*     ##  inform.i  --  control values for I/O and program flow  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     digits    decimal places output for energy and coordinates */
/*     iprint    steps between status printing (0=no printing) */
/*     iwrite    steps between coordinate dumps (0=no dumps) */
/*     isend     steps between socket communication (0=no sockets) */
/*     verbose   logical flag to turn on extra information */
/*     debug     logical flag to turn on full debug printing */
/*     holdup    logical flag to wait for carriage return on exit */
/*     abort     logical flag to stop execution at next chance */




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
/*     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  mdstuf.i  --  control of molecular dynamics trajectory  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nfree       total number of degrees of freedom for a system */
/*     velsave     flag to save velocity vector components to a file */
/*     frcsave     flag to save force vector components to a file */
/*     uindsave    flag to save induced atomic dipoles to a file */
/*     integrate   type of molecular dynamics integration algorithm */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  moldyn.i  --  velocity and acceleration on MD trajectory  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     v       current velocity of each atom along the x,y,z-axes */
/*     a       current acceleration of each atom along x,y,z-axes */
/*     aold    previous acceleration of each atom along x,y,z-axes */




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




/*     zero out the total mass and overall linear velocity */

    totmass = 0.;
    for (j = 1; j <= 3; ++j) {
	vtot[j - 1] = 0.;
    }

/*     compute linear velocity of the system center of mass */

    if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 0) {
	i__1 = group_1.ngrp;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    weigh = group_1.grpmass[i__ - 1];
	    totmass += weigh;
	    for (j = 1; j <= 3; ++j) {
		vtot[j - 1] += vcm_ref(j, i__) * weigh;
	    }
	}
    } else {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    weigh = atmtyp_1.mass[i__ - 1];
	    totmass += weigh;
	    for (j = 1; j <= 3; ++j) {
		vtot[j - 1] += v_ref(j, i__) * weigh;
	    }
	}
    }

/*     compute translational kinetic energy of overall system */

    etrans = 0.;
    for (j = 1; j <= 3; ++j) {
	vtot[j - 1] /= totmass;
/* Computing 2nd power */
	d__1 = vtot[j - 1];
	etrans += d__1 * d__1;
    }
    etrans = etrans * .5 * totmass / 418.4;

/*     find the center of mass coordinates of the overall system */

    if (! bound_1.use_bounds__) {
	xtot = 0.;
	ytot = 0.;
	ztot = 0.;
	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		xcm[i__ - 1] = 0.;
		ycm[i__ - 1] = 0.;
		zcm[i__ - 1] = 0.;
		i__2 = igrp_ref(2, i__);
		for (j = igrp_ref(1, i__); j <= i__2; ++j) {
		    k = group_1.kgrp[j - 1];
		    weigh = atmtyp_1.mass[k - 1];
		    xcm[i__ - 1] += atoms_1.x[k - 1] * weigh;
		    ycm[i__ - 1] += atoms_1.y[k - 1] * weigh;
		    zcm[i__ - 1] += atoms_1.z__[k - 1] * weigh;
		}
		xtot += xcm[i__ - 1];
		ytot += ycm[i__ - 1];
		ztot += zcm[i__ - 1];
/* Computing MAX */
		d__1 = 1., d__2 = group_1.grpmass[i__ - 1];
		weigh = max(d__1,d__2);
		xcm[i__ - 1] /= weigh;
		ycm[i__ - 1] /= weigh;
		zcm[i__ - 1] /= weigh;
	    }
	} else {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		weigh = atmtyp_1.mass[i__ - 1];
		xtot += atoms_1.x[i__ - 1] * weigh;
		ytot += atoms_1.y[i__ - 1] * weigh;
		ztot += atoms_1.z__[i__ - 1] * weigh;
	    }
	}
	xtot /= totmass;
	ytot /= totmass;
	ztot /= totmass;

/*     compute the angular momentum of the overall system */

	for (j = 1; j <= 3; ++j) {
	    mang[j - 1] = 0.;
	}
	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		weigh = group_1.grpmass[i__ - 1];
		mang[0] += (ycm[i__ - 1] * vcm_ref(3, i__) - zcm[i__ - 1] * 
			vcm_ref(2, i__)) * weigh;
		mang[1] += (zcm[i__ - 1] * vcm_ref(1, i__) - xcm[i__ - 1] * 
			vcm_ref(3, i__)) * weigh;
		mang[2] += (xcm[i__ - 1] * vcm_ref(2, i__) - ycm[i__ - 1] * 
			vcm_ref(1, i__)) * weigh;
	    }
	} else {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		weigh = atmtyp_1.mass[i__ - 1];
		mang[0] += (atoms_1.y[i__ - 1] * v_ref(3, i__) - atoms_1.z__[
			i__ - 1] * v_ref(2, i__)) * weigh;
		mang[1] += (atoms_1.z__[i__ - 1] * v_ref(1, i__) - atoms_1.x[
			i__ - 1] * v_ref(3, i__)) * weigh;
		mang[2] += (atoms_1.x[i__ - 1] * v_ref(2, i__) - atoms_1.y[
			i__ - 1] * v_ref(1, i__)) * weigh;
	    }
	}
	mang[0] -= (ytot * vtot[2] - ztot * vtot[1]) * totmass;
	mang[1] -= (ztot * vtot[0] - xtot * vtot[2]) * totmass;
	mang[2] -= (xtot * vtot[1] - ytot * vtot[0]) * totmass;

/*     calculate the moment of inertia tensor */

	xx = 0.;
	xy = 0.;
	xz = 0.;
	yy = 0.;
	yz = 0.;
	zz = 0.;
	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		weigh = group_1.grpmass[i__ - 1];
		xdel = xcm[i__ - 1] - xtot;
		ydel = ycm[i__ - 1] - ytot;
		zdel = zcm[i__ - 1] - ztot;
		xx += xdel * xdel * weigh;
		xy += xdel * ydel * weigh;
		xz += xdel * zdel * weigh;
		yy += ydel * ydel * weigh;
		yz += ydel * zdel * weigh;
		zz += zdel * zdel * weigh;
	    }
	} else {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		weigh = atmtyp_1.mass[i__ - 1];
		xdel = atoms_1.x[i__ - 1] - xtot;
		ydel = atoms_1.y[i__ - 1] - ytot;
		zdel = atoms_1.z__[i__ - 1] - ztot;
		xx += xdel * xdel * weigh;
		xy += xdel * ydel * weigh;
		xz += xdel * zdel * weigh;
		yy += ydel * ydel * weigh;
		yz += ydel * zdel * weigh;
		zz += zdel * zdel * weigh;
	    }
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

/*     fix to avoid singularity for one- or two-body systems */

	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    if (group_1.ngrp <= 2) {
		eps = 1e-6;
		tensor_ref(1, 1) = tensor_ref(1, 1) + eps;
		tensor_ref(2, 2) = tensor_ref(2, 2) + eps;
		tensor_ref(3, 3) = tensor_ref(3, 3) + eps;
	    }
	} else {
	    if (atoms_1.n <= 2) {
		eps = 1e-6;
		tensor_ref(1, 1) = tensor_ref(1, 1) + eps;
		tensor_ref(2, 2) = tensor_ref(2, 2) + eps;
		tensor_ref(3, 3) = tensor_ref(3, 3) + eps;
	    }
	}

/*     diagonalize the moment of inertia tensor */

	invert_(&c__3, &c__3, tensor);

/*     compute angular velocity and rotational kinetic energy */

	erot = 0.;
	for (i__ = 1; i__ <= 3; ++i__) {
	    vang[i__ - 1] = 0.;
	    for (j = 1; j <= 3; ++j) {
		vang[i__ - 1] += tensor_ref(i__, j) * mang[j - 1];
	    }
	    erot += vang[i__ - 1] * mang[i__ - 1];
	}
	erot = erot * .5 / 418.4;
    }

/*     eliminate any translation of the overall system */

    if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 0) {
	i__1 = group_1.ngrp;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		vcm_ref(j, i__) = vcm_ref(j, i__) - vtot[j - 1];
	    }
	}
    } else {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		v_ref(j, i__) = v_ref(j, i__) - vtot[j - 1];
	    }
	}
    }

/*     print the translational velocity of the overall system */

    if (inform_1.verbose) {
	io___28.ciunit = iounit_1.iout;
	s_wsfe(&io___28);
	for (i__ = 1; i__ <= 3; ++i__) {
	    do_fio(&c__1, (char *)&vtot[i__ - 1], (ftnlen)sizeof(doublereal));
	}
	do_fio(&c__1, (char *)&etrans, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     eliminate any rotation about the system center of mass */

    if (! bound_1.use_bounds__) {
	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		xdel = xcm[i__ - 1] - xtot;
		ydel = ycm[i__ - 1] - ytot;
		zdel = zcm[i__ - 1] - ztot;
		vcm_ref(1, i__) = vcm_ref(1, i__) - vang[1] * zdel + vang[2] *
			 ydel;
		vcm_ref(2, i__) = vcm_ref(2, i__) - vang[2] * xdel + vang[0] *
			 zdel;
		vcm_ref(3, i__) = vcm_ref(3, i__) - vang[0] * ydel + vang[1] *
			 xdel;
	    }
	} else {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		xdel = atoms_1.x[i__ - 1] - xtot;
		ydel = atoms_1.y[i__ - 1] - ytot;
		zdel = atoms_1.z__[i__ - 1] - ztot;
		v_ref(1, i__) = v_ref(1, i__) - vang[1] * zdel + vang[2] * 
			ydel;
		v_ref(2, i__) = v_ref(2, i__) - vang[2] * xdel + vang[0] * 
			zdel;
		v_ref(3, i__) = v_ref(3, i__) - vang[0] * ydel + vang[1] * 
			xdel;
	    }
	}

/*     print the angular velocity of the overall system */

	if (inform_1.verbose) {
	    io___29.ciunit = iounit_1.iout;
	    s_wsfe(&io___29);
	    for (i__ = 1; i__ <= 3; ++i__) {
		do_fio(&c__1, (char *)&vang[i__ - 1], (ftnlen)sizeof(
			doublereal));
	    }
	    do_fio(&c__1, (char *)&erot, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* mdrest_ */

#undef tensor_ref
#undef igrp_ref
#undef vcm_ref
#undef v_ref


