/* temper.f -- translated by f2c (version 20050501).
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
    doublereal kelvin0, kelvin, atmsph, tautemp, taupres, compress, collide, 
	    xnh[2], vnh[2], qnh[2], gnh[2], volmove;
    integer voltrial;
    logical isothermal, isobaric, anisotrop;
    char thermostat[11], barostat[10], volscale[9];
} bath_;

#define bath_1 bath_

struct {
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

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

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

/* Table of constant values */

static doublereal c_b11 = .66666666666666663;



/*     ############################################################# */
/*     ##  COPYRIGHT (C) 2003 by Alan Grossfield & Jay W. Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine temper  --  thermostat applied at half step  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "temper" applies a velocity correction at the half time step */
/*     as needed for the Nose-Hoover extended system thermostat */

/*     literature references: */

/*     D. Frenkel and B. Smit, "Understanding Molecular Simulation, */
/*     2nd Edition", Academic Press, San Diego, CA, 2002; see Appendix */
/*     E.2 for implementation details */

/*     G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein, */
/*     "Explicit Reversible Integrators for Extended Systems Dynamics", */
/*     Molecular Physics, 87, 1117-1157 (1996) */


/* Subroutine */ int temper_(doublereal *dt)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double exp(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal dt2, dt4, dt8, ekt, ekin[9]	/* was [3][3] */, 
	    scale, eksum;
    extern /* Subroutine */ int kinetic_(doublereal *, doublereal *);


#define v_ref(a_1,a_2) moldyn_1.v[(a_2)*3 + a_1 - 4]
#define vcm_ref(a_1,a_2) rgddyn_1.vcm[(a_2)*3 + a_1 - 4]
#define wcm_ref(a_1,a_2) rgddyn_1.wcm[(a_2)*3 + a_1 - 4]



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




/*     make half-step velocity correction for Nose-Hoover system */

    if (s_cmp(bath_1.thermostat, "NOSE-HOOVER", (ftnlen)11, (ftnlen)11) == 0) 
	    {
	kinetic_(&eksum, ekin);
	dt2 = *dt / 2.;
	dt4 = *dt / 4.;
	dt8 = *dt / 8.;
	ekt = bath_1.kelvin * .0019872066;
	bath_1.gnh[1] = (bath_1.qnh[0] * bath_1.vnh[0] * bath_1.vnh[0] - ekt) 
		/ bath_1.qnh[1];
	bath_1.vnh[1] += bath_1.gnh[1] * dt4;
	bath_1.vnh[0] *= exp(-bath_1.vnh[1] * dt8);
	bath_1.gnh[0] = (eksum * 2. - (doublereal) mdstuf_1.nfree * ekt) / 
		bath_1.qnh[0];
	bath_1.vnh[0] += bath_1.gnh[0] * dt4;
	bath_1.vnh[0] *= exp(-bath_1.vnh[1] * dt8);
	bath_1.xnh[0] += bath_1.vnh[0] * dt2;
	bath_1.xnh[1] += bath_1.vnh[1] * dt2;
	scale = exp(-bath_1.vnh[0] * dt2);
	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    vcm_ref(j, i__) = scale * vcm_ref(j, i__);
		    wcm_ref(j, i__) = scale * wcm_ref(j, i__);
		}
	    }
	} else {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (usage_1.use[i__ - 1]) {
		    for (j = 1; j <= 3; ++j) {
			v_ref(j, i__) = scale * v_ref(j, i__);
		    }
		}
	    }
	}
	eksum = eksum * scale * scale;
	bath_1.vnh[0] *= exp(-bath_1.vnh[1] * dt8);
	bath_1.gnh[0] = (eksum * 2. - (doublereal) mdstuf_1.nfree * ekt) / 
		bath_1.qnh[0];
	bath_1.vnh[0] += bath_1.gnh[0] * dt4;
	bath_1.vnh[0] *= exp(-bath_1.vnh[1] * dt8);
	bath_1.gnh[1] = (bath_1.qnh[0] * bath_1.vnh[0] * bath_1.vnh[0] - ekt) 
		/ bath_1.qnh[1];
	bath_1.vnh[1] += bath_1.gnh[1] * dt4;
    }
    return 0;
} /* temper_ */

#undef wcm_ref
#undef vcm_ref
#undef v_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine temper2  --  thermostat applied at full step  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "temper2" computes the instantaneous temperature and applies a */
/*     thermostat via Berendsen or Bussi-Parrinello velocity scaling, */
/*     Andersen stochastic collisions or Nose-Hoover extended system */

/*     literature references: */

/*     H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren, */
/*     A. DiNola and J. R. Hauk, "Molecular Dynamics with Coupling */
/*     to an External Bath", Journal of Chemical Physics, 81, */
/*     3684-3690 (1984) */

/*     G. Bussi and M. Parrinello, "Stochastic Thermostats: Comparison */
/*     of Local and Global Schemes", Computer Physics Communications, */
/*     179, 26-29 (2008) */

/*     H. C. Andersen, "Molecular Dynamics Simulations at Constant */
/*     Pressure and/or Temperature", Journal of Chemical Physics, */
/*     72, 2384-2393 (1980) */


/* Subroutine */ int temper2_(doublereal *dt, doublereal *eksum, doublereal *
	temp)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal), exp(doublereal), pow_dd(doublereal *, doublereal 
	    *);

    /* Local variables */
    static doublereal c__, d__;
    static integer i__, j;
    static doublereal r__, s, si, kt, dt2, dt4, dt8, ekt, rate, scale, speed, 
	    trial;
    extern doublereal random_(void), normal_(void);


#define v_ref(a_1,a_2) moldyn_1.v[(a_2)*3 + a_1 - 4]
#define vcm_ref(a_1,a_2) rgddyn_1.vcm[(a_2)*3 + a_1 - 4]
#define wcm_ref(a_1,a_2) rgddyn_1.wcm[(a_2)*3 + a_1 - 4]



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




/*     get the instantaneous temperature from the kinetic energy */

    *temp = *eksum * 2. / ((doublereal) mdstuf_1.nfree * .0019872066);
    if (! bath_1.isothermal) {
	return 0;
    }

/*     couple to external temperature bath via Berendsen scaling */

    if (s_cmp(bath_1.thermostat, "BERENDSEN", (ftnlen)11, (ftnlen)9) == 0) {
	if (*temp == 0.) {
	    *temp = .1;
	}
	scale = sqrt(*dt / bath_1.tautemp * (bath_1.kelvin / *temp - 1.) + 1.)
		;
	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    vcm_ref(j, i__) = scale * vcm_ref(j, i__);
		    wcm_ref(j, i__) = scale * wcm_ref(j, i__);
		}
	    }
	} else {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (usage_1.use[i__ - 1]) {
		    for (j = 1; j <= 3; ++j) {
			v_ref(j, i__) = scale * v_ref(j, i__);
		    }
		}
	    }
	}

/*     couple to external temperature bath via Bussi scaling */

    } else if (s_cmp(bath_1.thermostat, "BUSSI", (ftnlen)11, (ftnlen)5) == 0) 
	    {
	if (*temp == 0.) {
	    *temp = .1;
	}
	c__ = exp(-(*dt) / bath_1.tautemp);
	d__ = (1. - c__) * (bath_1.kelvin / *temp) / (doublereal) 
		mdstuf_1.nfree;
	r__ = normal_();
	s = 0.;
	i__1 = mdstuf_1.nfree - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    si = normal_();
	    s += si * si;
	}
	scale = c__ + (s + r__ * r__) * d__ + r__ * 2. * sqrt(c__ * d__);
	scale = sqrt(scale);
	if (r__ + sqrt(c__ / d__) < 0.) {
	    scale = -scale;
	}
	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    vcm_ref(j, i__) = scale * vcm_ref(j, i__);
		    wcm_ref(j, i__) = scale * wcm_ref(j, i__);
		}
	    }
	} else {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (usage_1.use[i__ - 1]) {
		    for (j = 1; j <= 3; ++j) {
			v_ref(j, i__) = scale * v_ref(j, i__);
		    }
		}
	    }
	}

/*     select random velocities via Andersen stochastic collisions */

    } else if (s_cmp(bath_1.thermostat, "ANDERSEN", (ftnlen)11, (ftnlen)8) == 
	    0) {
	kt = bath_1.kelvin * .83144725;
	rate = *dt * 1e3 * bath_1.collide;
	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    d__1 = (doublereal) group_1.ngrp;
	    rate /= pow_dd(&d__1, &c_b11);
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		trial = random_();
		if (trial < rate) {
		    speed = sqrt(kt / group_1.grpmass[i__ - 1]);
		    for (j = 1; j <= 3; ++j) {
			vcm_ref(j, i__) = speed * normal_();
		    }
		}
	    }
	} else {
	    d__1 = (doublereal) usage_1.nuse;
	    rate /= pow_dd(&d__1, &c_b11);
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (usage_1.use[i__ - 1]) {
		    trial = random_();
		    if (trial < rate) {
			speed = sqrt(kt / atmtyp_1.mass[i__ - 1]);
			for (j = 1; j <= 3; ++j) {
			    v_ref(j, i__) = speed * normal_();
			}
		    }
		}
	    }
	}

/*     make full-step velocity correction for Nose-Hoover system */

    } else if (s_cmp(bath_1.thermostat, "NOSE-HOOVER", (ftnlen)11, (ftnlen)11)
	     == 0) {
	dt2 = *dt / 2.;
	dt4 = *dt / 4.;
	dt8 = *dt / 8.;
	ekt = bath_1.kelvin * .0019872066;
	bath_1.gnh[1] = (bath_1.qnh[0] * bath_1.vnh[0] * bath_1.vnh[0] - ekt) 
		/ bath_1.qnh[1];
	bath_1.vnh[1] += bath_1.gnh[1] * dt4;
	bath_1.vnh[0] *= exp(-bath_1.vnh[1] * dt8);
	bath_1.gnh[0] = (*eksum * 2. - (doublereal) mdstuf_1.nfree * ekt) / 
		bath_1.qnh[0];
	bath_1.vnh[0] += bath_1.gnh[0] * dt4;
	bath_1.vnh[0] *= exp(-bath_1.vnh[1] * dt8);
	bath_1.xnh[0] += bath_1.vnh[0] * dt2;
	bath_1.xnh[1] += bath_1.vnh[1] * dt2;
	scale = exp(-bath_1.vnh[0] * dt2);
	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    vcm_ref(j, i__) = scale * vcm_ref(j, i__);
		    wcm_ref(j, i__) = scale * wcm_ref(j, i__);
		}
	    }
	} else {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (usage_1.use[i__ - 1]) {
		    for (j = 1; j <= 3; ++j) {
			v_ref(j, i__) = scale * v_ref(j, i__);
		    }
		}
	    }
	}
	*eksum = *eksum * scale * scale;
	bath_1.vnh[0] *= exp(-bath_1.vnh[1] * dt8);
	bath_1.gnh[0] = (*eksum * 2. - (doublereal) mdstuf_1.nfree * ekt) / 
		bath_1.qnh[0];
	bath_1.vnh[0] += bath_1.gnh[0] * dt4;
	bath_1.vnh[0] *= exp(-bath_1.vnh[1] * dt8);
	bath_1.gnh[1] = (bath_1.qnh[0] * bath_1.vnh[0] * bath_1.vnh[0] - ekt) 
		/ bath_1.qnh[1];
	bath_1.vnh[1] += bath_1.gnh[1] * dt4;
    }
    return 0;
} /* temper2_ */

#undef wcm_ref
#undef vcm_ref
#undef v_ref


