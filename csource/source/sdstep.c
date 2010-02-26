/* sdstep.f -- translated by f2c (version 20050501).
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
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

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
    doublereal krat[25000];
    integer nrat, nratx, irat[50000]	/* was [2][25000] */, iratx[25000], 
	    kratx[25000];
    logical ratimage[25000], use_rattle__;
} shake_;

#define shake_1 shake_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal vir[9]	/* was [3][3] */;
} virial_;

#define virial_1 virial_

struct {
    doublereal kelvin0, kelvin, atmsph, tautemp, taupres, compress, collide, 
	    xnh[2], vnh[2], qnh[2], gnh[2], volmove;
    integer voltrial;
    logical isothermal, isobaric, anisotrop;
    char thermostat[11], barostat[10], volscale[9];
} bath_;

#define bath_1 bath_

struct {
    doublereal friction, fgamma[25000];
    logical use_sdarea__;
} stodyn_;

#define stodyn_1 stodyn_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal rad[5000], eps[5000], rad4[5000], eps4[5000], reduct[5000];
} kvdws_;

#define kvdws_1 kvdws_



/*     ############################################################## */
/*     ##  COPYRIGHT (C) 1998 by Rohit Pappu & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                    ## */
/*     ############################################################## */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine sdstep  --  Verlet stochastic dynamics step  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "sdstep" performs a single stochastic dynamics time step */
/*     via a velocity Verlet integration algorithm */

/*     literature references: */

/*     M. P. Allen, "Brownian Dynamics Simulation of a Chemical */
/*     Reaction in Solution", Molecular Physics, 40, 1073-1087 (1980) */

/*     F. Guarnieri and W. C. Still, "A Rapidly Convergent Simulation */
/*     Method: Mixed Monte Carlo/Stochastic Dynamics", Journal of */
/*     Computational Chemistry, 15, 1302-1310 (1994) */


/* Subroutine */ int sdstep_(integer *istep, doublereal *dt)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *), 
	    pressure_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer i__, j;
    static doublereal vxx, vyx, vyy, vzx, vzz, vzy, ekin[9]	/* was [3][3] 
	    */, temp, term, epot, pres, etot, xold[25000], yold[25000], zold[
	    25000], afric[25000], pfric[25000], prand[75000]	/* was [3][
	    25000] */, vfric[25000], vrand[75000]	/* was [3][25000] */, 
	    eksum;
    extern /* Subroutine */ int mdsave_(integer *, doublereal *, doublereal *)
	    ;
    static doublereal derivs[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int sdterm_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *), rattle_(
	    doublereal *, doublereal *, doublereal *, doublereal *), mdstat_(
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static doublereal stress[9]	/* was [3][3] */;
    extern /* Subroutine */ int rattle2_(doublereal *), kinetic_(doublereal *,
	     doublereal *);


#define a_ref(a_1,a_2) moldyn_1.a[(a_2)*3 + a_1 - 4]
#define v_ref(a_1,a_2) moldyn_1.v[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define prand_ref(a_1,a_2) prand[(a_2)*3 + a_1 - 4]
#define vrand_ref(a_1,a_2) vrand[(a_2)*3 + a_1 - 4]
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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  shake.i  --  definition of Shake/Rattle constraints  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     krat         ideal distance value for rattle constraint */
/*     nrat         number of rattle distance constraints to apply */
/*     nratx        number of atom group spatial constraints to apply */
/*     irat         atom numbers of atoms in a rattle constraint */
/*     iratx        group number of group in a spatial constraint */
/*     kratx        spatial constraint type (1=plane, 2=line, 3=point) */
/*     ratimage     flag to use minimum image for rattle constraint */
/*     use_rattle   logical flag to set use of rattle contraints */




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




/*     get frictional and random terms for position and velocity */

    sdterm_(istep, dt, pfric, vfric, afric, prand, vrand);

/*     store the current atom positions, then find new atom */
/*     positions and half-step velocities via Verlet recursion */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    xold[i__ - 1] = atoms_1.x[i__ - 1];
	    yold[i__ - 1] = atoms_1.y[i__ - 1];
	    zold[i__ - 1] = atoms_1.z__[i__ - 1];
	    atoms_1.x[i__ - 1] = atoms_1.x[i__ - 1] + v_ref(1, i__) * vfric[
		    i__ - 1] + a_ref(1, i__) * afric[i__ - 1] + prand_ref(1, 
		    i__);
	    atoms_1.y[i__ - 1] = atoms_1.y[i__ - 1] + v_ref(2, i__) * vfric[
		    i__ - 1] + a_ref(2, i__) * afric[i__ - 1] + prand_ref(2, 
		    i__);
	    atoms_1.z__[i__ - 1] = atoms_1.z__[i__ - 1] + v_ref(3, i__) * 
		    vfric[i__ - 1] + a_ref(3, i__) * afric[i__ - 1] + 
		    prand_ref(3, i__);
	    for (j = 1; j <= 3; ++j) {
		v_ref(j, i__) = v_ref(j, i__) * pfric[i__ - 1] + a_ref(j, i__)
			 * .5 * vfric[i__ - 1];
	    }
	}
    }

/*     get constraint-corrected positions and half-step velocities */

    if (shake_1.use_rattle__) {
	rattle_(dt, xold, yold, zold);
    }

/*     get the potential energy and atomic forces */

    gradient_(&epot, derivs);

/*     use Newton's second law to get the next accelerations; */
/*     find the full-step velocities using the Verlet recursion */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    for (j = 1; j <= 3; ++j) {
		a_ref(j, i__) = derivs_ref(j, i__) * -418.4 / atmtyp_1.mass[
			i__ - 1];
		v_ref(j, i__) = v_ref(j, i__) + a_ref(j, i__) * .5 * vfric[
			i__ - 1] + vrand_ref(j, i__);
	    }
	}
    }

/*     correct internal virial to account for frictional forces */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    term = vfric[i__ - 1] / *dt - 1.;
	    vxx = term * atoms_1.x[i__ - 1] * derivs_ref(1, i__);
	    vyx = term * .5 * (atoms_1.y[i__ - 1] * derivs_ref(1, i__) + 
		    atoms_1.x[i__ - 1] * derivs_ref(2, i__));
	    vzx = term * .5 * (atoms_1.z__[i__ - 1] * derivs_ref(1, i__) + 
		    atoms_1.x[i__ - 1] * derivs_ref(3, i__));
	    vyy = term * atoms_1.y[i__ - 1] * derivs_ref(2, i__);
	    vzy = term * .5 * (atoms_1.z__[i__ - 1] * derivs_ref(2, i__) + 
		    atoms_1.y[i__ - 1] * derivs_ref(3, i__));
	    vzz = term * atoms_1.z__[i__ - 1] * derivs_ref(3, i__);
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

/*     find the constraint-corrected full-step velocities */

    if (shake_1.use_rattle__) {
	rattle2_(dt);
    }

/*     accumulate the kinetic energy and its outer product */

    kinetic_(&eksum, ekin);

/*     compute and control the temperature and pressure */

    temp = eksum * 2. / ((doublereal) mdstuf_1.nfree * .0019872066);
    pressure_(dt, &epot, ekin, &temp, &pres, stress);

/*     system energy is sum of kinetic and potential energies */

    etot = eksum + epot;

/*     compute statistics and save trajectory for this step */

    mdstat_(istep, dt, &etot, &epot, &eksum, &temp, &pres);
    mdsave_(istep, dt, &epot);
    return 0;
} /* sdstep_ */

#undef derivs_ref
#undef vrand_ref
#undef prand_ref
#undef vir_ref
#undef v_ref
#undef a_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine sdterm  --  frictional and random SD terms  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "sdterm" gets frictional and random force terms needed to */
/*     update positions and velocities via stochastic dynamics */


/* Subroutine */ int sdterm_(integer *istep, doublereal *dt, doublereal *
	pfric, doublereal *vfric, doublereal *afric, doublereal *prand, 
	doublereal *vrand)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal gdt, rho, ktm, gdt2, gdt3, gdt4, gdt5, gdt6, gdt7, gdt8,
	     gdt9, egdt, rhoc, psig, vsig, pterm, pnorm, vterm, vnorm;
    extern /* Subroutine */ int sdarea_(integer *);
    extern doublereal normal_(void);


#define prand_ref(a_1,a_2) prand[(a_2)*3 + a_1]
#define vrand_ref(a_1,a_2) vrand[(a_2)*3 + a_1]



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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  bath.i  --  temperature and pressure control parameters  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxnose     maximum length of the Nose-Hoover chain */

/*     kelvin0     target value for the system temperature (K) */
/*     kelvin      variable target temperature for thermostat (K) */
/*     atmsph      target value for the system pressure (atm) */
/*     tautemp     time constant for Berendsen thermostat (psec) */
/*     taupres     time constant for Berendsen barostat (psec) */
/*     compress    isothermal compressibility of medium (atm-1) */
/*     collide     collision frequency for Andersen thermostat */
/*     xnh         position of each chained Nose-Hoover thermostat */
/*     vnh         velocity of each chained Nose-Hoover thermostat */
/*     qnh         mass for each chained Nose-Hoover thermostat */
/*     gnh         coupling between chained Nose-Hoover thermostats */
/*     volmove     maximum volume move for Monte Carlo barostat (Ang**3) */
/*     voltrial    mean number of steps between Monte Carlo moves */
/*     isothermal  logical flag governing use of temperature control */
/*     isobaric    logical flag governing use of pressure control */
/*     anisotrop   logical flag governing use of anisotropic pressure */
/*     thermostat  choice of temperature control method to be used */
/*     barostat    choice of pressure control method to be used */
/*     volscale    choice of scaling method for Monte Carlo barostat */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  stodyn.i  --  frictional coefficients for SD trajectory  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     friction    global frictional coefficient for exposed particle */
/*     fgamma      atomic frictional coefficients for each atom */
/*     use_sdarea  logical flag to use surface area friction scaling */




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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     set the atomic friction coefficients to the global value */

    /* Parameter adjustments */
    vrand -= 4;
    prand -= 4;
    --afric;
    --vfric;
    --pfric;

    /* Function Body */
    if (*istep == 1) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    stodyn_1.fgamma[i__ - 1] = stodyn_1.friction;
	}
    }

/*     set the value of the friction coefficient for each atom */

    if (stodyn_1.use_sdarea__) {
	sdarea_(istep);
    }

/*     get the frictional and random terms for stochastic dynamics */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    gdt = stodyn_1.fgamma[i__ - 1] * *dt;

/*     stochastic dynamics reduces to simple MD for zero friction */

	    if (gdt <= 0.) {
		pfric[i__] = 1.;
		vfric[i__] = *dt;
		afric[i__] = *dt * .5 * *dt;
		for (j = 1; j <= 3; ++j) {
		    prand_ref(j, i__) = 0.;
		    vrand_ref(j, i__) = 0.;
		}

/*     analytical expressions when friction coefficient is large */

	    } else {
		if (gdt >= .05) {
		    egdt = exp(-gdt);
		    pfric[i__] = egdt;
		    vfric[i__] = (1. - egdt) / stodyn_1.fgamma[i__ - 1];
		    afric[i__] = (*dt - vfric[i__]) / stodyn_1.fgamma[i__ - 1]
			    ;
		    pterm = gdt * 2. - 3. + (4. - egdt) * egdt;
/* Computing 2nd power */
		    d__1 = egdt;
		    vterm = 1. - d__1 * d__1;
/* Computing 2nd power */
		    d__1 = 1. - egdt;
		    rho = d__1 * d__1 / sqrt(pterm * vterm);

/*     use series expansions when friction coefficient is small */

		} else {
		    gdt2 = gdt * gdt;
		    gdt3 = gdt * gdt2;
		    gdt4 = gdt2 * gdt2;
		    gdt5 = gdt2 * gdt3;
		    gdt6 = gdt3 * gdt3;
		    gdt7 = gdt3 * gdt4;
		    gdt8 = gdt4 * gdt4;
		    gdt9 = gdt4 * gdt5;
/* Computing 2nd power */
		    d__1 = stodyn_1.fgamma[i__ - 1];
		    afric[i__] = (gdt2 / 2. - gdt3 / 6. + gdt4 / 24. - gdt5 / 
			    120. + gdt6 / 720. - gdt7 / 5040. + gdt8 / 40320. 
			    - gdt9 / 362880.) / (d__1 * d__1);
		    vfric[i__] = *dt - stodyn_1.fgamma[i__ - 1] * afric[i__];
		    pfric[i__] = 1. - stodyn_1.fgamma[i__ - 1] * vfric[i__];
		    pterm = gdt3 * 2. / 3. - gdt4 / 2. + gdt5 * 7. / 30. - 
			    gdt6 / 12. + gdt7 * 31. / 1260. - gdt8 / 160. + 
			    gdt9 * 127. / 90720.;
		    vterm = gdt * 2. - gdt2 * 2. + gdt3 * 4. / 3. - gdt4 * 2. 
			    / 3. + gdt5 * 4. / 15. - gdt6 * 4. / 45. + gdt7 * 
			    8. / 315. - gdt8 * 2. / 315. + gdt9 * 4. / 2835.;
		    rho = sqrt(3.) * (.5 - gdt / 16. - gdt2 * 17. / 1280. + 
			    gdt3 * 17. / 6144. + gdt4 * 40967. / 34406400. - 
			    gdt5 * 57203. / 275251200. - gdt6 * 1429487. / 
			    13212057600. + gdt7 * 1877509. / 105696460800.);
		}

/*     compute random terms to thermostat the nonzero friction case */

		ktm = bath_1.kelvin * .83144725 / atmtyp_1.mass[i__ - 1];
		psig = sqrt(ktm * pterm) / stodyn_1.fgamma[i__ - 1];
		vsig = sqrt(ktm * vterm);
		rhoc = sqrt(1. - rho * rho);
		for (j = 1; j <= 3; ++j) {
		    pnorm = normal_();
		    vnorm = normal_();
		    prand_ref(j, i__) = psig * pnorm;
		    vrand_ref(j, i__) = vsig * (rho * pnorm + rhoc * vnorm);
		}
	    }
	}
    }
    return 0;
} /* sdterm_ */

#undef vrand_ref
#undef prand_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine sdarea  --  scale SD friction coefficients  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "sdarea" optionally scales the atomic friction coefficient */
/*     of each atom based on its accessible surface area */

/*     literature reference: */

/*     S. Yun-Yi, W. Lu and W. F. van Gunsteren, "On the Approximation */
/*     of Solvent Effects on the Conformation and Dynamics of */
/*     Cyclosporin A by Stochastic Dynamics Simulation Techniques", */
/*     Molecular Simulation, 1, 369-383 (1988) */


/* Subroutine */ int sdarea_(integer *istep)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ int surfatom_(integer *, doublereal *, doublereal 
	    *);
    static integer i__;
    static doublereal area, probe, ratio, radius[25000];
    static integer resurf, modstep;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]



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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  couple.i  --  near-neighbor atom connectivity lists   ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxn13   maximum number of atoms 1-3 connected to an atom */
/*     maxn14   maximum number of atoms 1-4 connected to an atom */
/*     maxn15   maximum number of atoms 1-5 connected to an atom */

/*     n12      number of atoms directly bonded to each atom */
/*     i12      atom numbers of atoms 1-2 connected to each atom */
/*     n13      number of atoms in a 1-3 relation to each atom */
/*     i13      atom numbers of atoms 1-3 connected to each atom */
/*     n14      number of atoms in a 1-4 relation to each atom */
/*     i14      atom numbers of atoms 1-4 connected to each atom */
/*     n15      number of atoms in a 1-5 relation to each atom */
/*     i15      atom numbers of atoms 1-5 connected to each atom */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kvdws.i  --  forcefield parameters for van der Waals terms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     rad      van der Waals radius parameter for each atom type */
/*     eps      van der Waals well depth parameter for each atom type */
/*     rad4     van der Waals radius parameter in 1-4 interactions */
/*     eps4     van der Waals well depth parameter in 1-4 interactions */
/*     reduct   van der Waals reduction factor for each atom type */




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
/*     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  stodyn.i  --  frictional coefficients for SD trajectory  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     friction    global frictional coefficient for exposed particle */
/*     fgamma      atomic frictional coefficients for each atom */
/*     use_sdarea  logical flag to use surface area friction scaling */




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




/*     determine new friction coefficients every few SD steps */

    resurf = 100;
    modstep = *istep % resurf;
    if (modstep != 1) {
	return 0;
    }

/*     set the atomic radii to estimates of sigma values */

    probe = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	radius[i__ - 1] = kvdws_1.rad[atmtyp_1.class__[i__ - 1] - 1] / 
		1.122462048309372981;
	if (radius[i__ - 1] != 0.) {
	    radius[i__ - 1] += probe;
	}
    }

/*     scale atomic friction coefficients by accessible area */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    if (radius[i__ - 1] != 0.) {
		surfatom_(&i__, &area, radius);
/* Computing 2nd power */
		d__1 = radius[i__ - 1];
		ratio = area / (d__1 * d__1 * 12.566370614359172);
		stodyn_1.fgamma[i__ - 1] = ratio * stodyn_1.friction;
	    }
	}
    }

/*     monovalent atoms with zero radius get attached atom value */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    if (radius[i__ - 1] == 0. && couple_1.n12[i__ - 1] == 1) {
		stodyn_1.fgamma[i__ - 1] = stodyn_1.fgamma[i12_ref(1, i__) - 
			1];
	    }
	}
    }
    return 0;
} /* sdarea_ */

#undef i12_ref


