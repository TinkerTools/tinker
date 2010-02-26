/* beeman.f -- translated by f2c (version 20050501).
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



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine beeman  --  Beeman molecular dynamics step  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "beeman" performs a single molecular dynamics time step */
/*     by means of a Beeman multistep recursion formula; the */
/*     actual coefficients are Brooks' "Better Beeman" values */

/*     literature references: */

/*     D. Beeman, "Some Multistep Methods for Use in Molecular */
/*     Dynamics Calculations", Journal of Computational Physics, */
/*     20, 130-139 (1976) */

/*     B. R. Brooks, "Algorithms for Molecular Dynamics at Constant */
/*     Temperature and Pressure", DCRT Report, NIH, April 1988 */


/* Subroutine */ int beeman_(integer *istep, doublereal *dt)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *), 
	    pressure_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer i__, j;
    static doublereal dt_8__, dt2_8__, ekin[9]	/* was [3][3] */, temp, epot, 
	    xold[25000], pres, etot, yold[25000], zold[25000], eksum, xterm, 
	    yterm, zterm;
    extern /* Subroutine */ int mdsave_(integer *, doublereal *, doublereal *)
	    ;
    static doublereal derivs[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int temper_(doublereal *), rattle_(doublereal *, 
	    doublereal *, doublereal *, doublereal *), mdstat_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal stress[9]	/* was [3][3] */;
    extern /* Subroutine */ int rattle2_(doublereal *), temper2_(doublereal *,
	     doublereal *, doublereal *), kinetic_(doublereal *, doublereal *)
	    ;


#define a_ref(a_1,a_2) moldyn_1.a[(a_2)*3 + a_1 - 4]
#define v_ref(a_1,a_2) moldyn_1.v[(a_2)*3 + a_1 - 4]
#define aold_ref(a_1,a_2) moldyn_1.aold[(a_2)*3 + a_1 - 4]
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




/*     set some time values for the dynamics integration */

    dt_8__ = *dt * .125;
    dt2_8__ = *dt * dt_8__;

/*     make half-step temperature and pressure corrections */

    temper_(dt);

/*     store the current atom positions, then find new atom */
/*     positions and half-step velocities via Beeman recursion */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    xold[i__ - 1] = atoms_1.x[i__ - 1];
	    yold[i__ - 1] = atoms_1.y[i__ - 1];
	    zold[i__ - 1] = atoms_1.z__[i__ - 1];
	    xterm = a_ref(1, i__) * 5. - aold_ref(1, i__);
	    yterm = a_ref(2, i__) * 5. - aold_ref(2, i__);
	    zterm = a_ref(3, i__) * 5. - aold_ref(3, i__);
	    atoms_1.x[i__ - 1] = atoms_1.x[i__ - 1] + v_ref(1, i__) * *dt + 
		    xterm * dt2_8__;
	    atoms_1.y[i__ - 1] = atoms_1.y[i__ - 1] + v_ref(2, i__) * *dt + 
		    yterm * dt2_8__;
	    atoms_1.z__[i__ - 1] = atoms_1.z__[i__ - 1] + v_ref(3, i__) * *dt 
		    + zterm * dt2_8__;
	    v_ref(1, i__) = v_ref(1, i__) + xterm * dt_8__;
	    v_ref(2, i__) = v_ref(2, i__) + yterm * dt_8__;
	    v_ref(3, i__) = v_ref(3, i__) + zterm * dt_8__;
	}
    }

/*     get constraint-corrected positions and half-step velocities */

    if (shake_1.use_rattle__) {
	rattle_(dt, xold, yold, zold);
    }

/*     get the potential energy and atomic forces */

    gradient_(&epot, derivs);

/*     use Newton's second law to get the next accelerations; */
/*     find the full-step velocities using the Beeman recursion */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    for (j = 1; j <= 3; ++j) {
		aold_ref(j, i__) = a_ref(j, i__);
		a_ref(j, i__) = derivs_ref(j, i__) * -418.4 / atmtyp_1.mass[
			i__ - 1];
		v_ref(j, i__) = v_ref(j, i__) + (a_ref(j, i__) * 3. + 
			aold_ref(j, i__)) * dt_8__;
	    }
	}
    }

/*     find the constraint-corrected full-step velocities */

    if (shake_1.use_rattle__) {
	rattle2_(dt);
    }

/*     accumulate the kinetic energy and its outer product */

    kinetic_(&eksum, ekin);

/*     make full-step temperature and pressure corrections */

    temper2_(dt, &eksum, &temp);
    pressure_(dt, &epot, ekin, &temp, &pres, stress);

/*     system energy is sum of kinetic and potential energies */

    etot = eksum + epot;

/*     compute statistics and save trajectory for this step */

    mdstat_(istep, dt, &etot, &epot, &eksum, &temp, &pres);
    mdsave_(istep, dt, &epot);
    return 0;
} /* beeman_ */

#undef derivs_ref
#undef aold_ref
#undef v_ref
#undef a_ref


