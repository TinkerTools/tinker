/* pressure.f -- translated by f2c (version 20050501).
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
    doublereal kelvin0, kelvin, atmsph, tautemp, taupres, compress, collide, 
	    xnh[2], vnh[2], qnh[2], gnh[2], volmove;
    integer voltrial;
    logical isothermal, isobaric, anisotrop;
    char thermostat[11], barostat[10], volscale[9];
} bath_;

#define bath_1 bath_

struct {
    doublereal xbox, ybox, zbox, alpha, beta, gamma, xbox2, ybox2, zbox2, 
	    box34, lvec[9]	/* was [3][3] */, recip[9]	/* was [3][3] 
	    */, volbox, beta_sin__, beta_cos__, gamma_sin__, gamma_cos__, 
	    beta_term__, gamma_term__;
    logical orthogonal, monoclinic, triclinic, octahedron;
    char spacegrp[10];
} boxes_;

#define boxes_1 boxes_

struct {
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

struct {
    doublereal vir[9]	/* was [3][3] */;
} virial_;

#define virial_1 virial_

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
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine pressure  --  constant pressure via barostat  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "pressure" uses the internal virial to find the pressure */
/*     in a periodic box and maintains a constant desired pressure */
/*     via a barostat method */


/* Subroutine */ int pressure_(doublereal *dt, doublereal *epot, doublereal *
	ekin, doublereal *temp, doublereal *pres, doublereal *stress)
{
    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int pscale_(doublereal *, doublereal *, 
	    doublereal *);
    static doublereal factor;
    extern /* Subroutine */ int pmonte_(doublereal *, doublereal *);


#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define ekin_ref(a_1,a_2) ekin[(a_2)*3 + a_1]
#define stress_ref(a_1,a_2) stress[(a_2)*3 + a_1]



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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  boxes.i  --  parameters for periodic boundary conditions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     xbox        length of a-axis of periodic box in Angstroms */
/*     ybox        length of b-axis of periodic box in Angstroms */
/*     zbox        length of c-axis of periodic box in Angstroms */
/*     alpha       angle between b- and c-axes of box in degrees */
/*     beta        angle between a- and c-axes of box in degrees */
/*     gamma       angle between a- and b-axes of box in degrees */
/*     xbox2       half of the a-axis length of periodic box */
/*     ybox2       half of the b-axis length of periodic box */
/*     zbox2       half of the c-axis length of periodic box */
/*     box34       three-fourths axis length of truncated octahedron */
/*     lvec        real space lattice vectors as matrix rows */
/*     recip       reciprocal lattice vectors as matrix columns */
/*     volbox      volume in Ang**3 of the periodic box */
/*     beta_sin    sine of the beta periodic box angle */
/*     beta_cos    cosine of the beta periodic box angle */
/*     gamma_sin   sine of the gamma periodic box angle */
/*     gamma_cos   cosine of the gamma periodic box angle */
/*     beta_term   term used in generating triclinic box */
/*     gamma_term  term used in generating triclinic box */
/*     orthogonal  flag to mark periodic box as orthogonal */
/*     monoclinic  flag to mark periodic box as monoclinic */
/*     triclinic   flag to mark periodic box as triclinic */
/*     octahedron  flag to mark box as truncated octahedron */
/*     spacegrp    space group symbol for the unitcell type */




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




/*     only necessary if periodic boundaries are in use */

    /* Parameter adjustments */
    stress -= 4;
    ekin -= 4;

    /* Function Body */
    if (! bound_1.use_bounds__) {
	return 0;
    }

/*     calculate the stress tensor for anisotropic systems */

    factor = 68568.4112 / boxes_1.volbox;
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    stress_ref(j, i__) = factor * (ekin_ref(j, i__) * 2. - vir_ref(j, 
		    i__));
	}
    }

/*     set isotropic pressure to the average of tensor diagonal */

    *pres = (stress_ref(1, 1) + stress_ref(2, 2) + stress_ref(3, 3)) / 3.;

/*     use either the Berendsen or Monte Carlo barostat method */

    if (bath_1.isobaric) {
	if (s_cmp(bath_1.barostat, "BERENDSEN", (ftnlen)10, (ftnlen)9) == 0) {
	    pscale_(dt, pres, &stress[4]);
	} else if (s_cmp(bath_1.barostat, "MONTECARLO", (ftnlen)10, (ftnlen)
		10) == 0) {
	    pmonte_(epot, temp);
	}
    }
    return 0;
} /* pressure_ */

#undef stress_ref
#undef ekin_ref
#undef vir_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine pscale  --  Berendsen barostat via scaling  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "pscale" implements a Berendsen barostat by scaling the */
/*     coordinates and box dimensions via coupling to an external */
/*     constant pressure bath */

/*     literature references: */

/*     H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren, */
/*     A. DiNola and J. R. Hauk, "Molecular Dynamics with Coupling */
/*     to an External Bath", Journal of Chemical Physics, 81, */
/*     3684-3690 (1984) */

/*     S. E. Feller, Y. Zhang, R. W. Pastor, B. R. Brooks, "Constant */
/*     Pressure Molecular Dynamics Simulation: The Langevin Piston */
/*     Method", Journal of Chemical Physics, 103, 4613-4621 (1995) */

/*     code for anisotropic pressure coupling was provided by Guido */
/*     Raos, Dipartimento di Chimica, Politecnico di Milano, Italy */


/* Subroutine */ int pscale_(doublereal *dt, doublereal *pres, doublereal *
	stress)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal), acos(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal xcm, ycm, zcm, hbox[9]	/* was [3][3] */, temp[9]	
	    /* was [3][3] */;
    static integer stop;
    static doublereal scale, weigh, third;
    static integer start;
    static doublereal xmove, ymove, zmove, ascale[9]	/* was [3][3] */, 
	    cosine;
    extern /* Subroutine */ int lattice_(void);


#define hbox_ref(a_1,a_2) hbox[(a_2)*3 + a_1 - 4]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]
#define temp_ref(a_1,a_2) temp[(a_2)*3 + a_1 - 4]
#define ascale_ref(a_1,a_2) ascale[(a_2)*3 + a_1 - 4]
#define stress_ref(a_1,a_2) stress[(a_2)*3 + a_1]



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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  boxes.i  --  parameters for periodic boundary conditions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     xbox        length of a-axis of periodic box in Angstroms */
/*     ybox        length of b-axis of periodic box in Angstroms */
/*     zbox        length of c-axis of periodic box in Angstroms */
/*     alpha       angle between b- and c-axes of box in degrees */
/*     beta        angle between a- and c-axes of box in degrees */
/*     gamma       angle between a- and b-axes of box in degrees */
/*     xbox2       half of the a-axis length of periodic box */
/*     ybox2       half of the b-axis length of periodic box */
/*     zbox2       half of the c-axis length of periodic box */
/*     box34       three-fourths axis length of truncated octahedron */
/*     lvec        real space lattice vectors as matrix rows */
/*     recip       reciprocal lattice vectors as matrix columns */
/*     volbox      volume in Ang**3 of the periodic box */
/*     beta_sin    sine of the beta periodic box angle */
/*     beta_cos    cosine of the beta periodic box angle */
/*     gamma_sin   sine of the gamma periodic box angle */
/*     gamma_cos   cosine of the gamma periodic box angle */
/*     beta_term   term used in generating triclinic box */
/*     gamma_term  term used in generating triclinic box */
/*     orthogonal  flag to mark periodic box as orthogonal */
/*     monoclinic  flag to mark periodic box as monoclinic */
/*     triclinic   flag to mark periodic box as triclinic */
/*     octahedron  flag to mark box as truncated octahedron */
/*     spacegrp    space group symbol for the unitcell type */




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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     find the isotropic scale factor for constant pressure */

    /* Parameter adjustments */
    stress -= 4;

    /* Function Body */
    if (! bath_1.anisotrop) {
	scale = 1.;
	third = .33333333333333331;
	d__1 = *dt * bath_1.compress / bath_1.taupres * (*pres - 
		bath_1.atmsph) + 1.;
	scale = pow_dd(&d__1, &third);

/*     modify the current periodic box dimension values */

	boxes_1.xbox = scale * boxes_1.xbox;
	boxes_1.ybox = scale * boxes_1.ybox;
	boxes_1.zbox = scale * boxes_1.zbox;

/*     propagate the new box dimensions to other lattice values */

	lattice_();

/*     couple to pressure bath via atom scaling in Cartesian space */

	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) != 
		0) {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (usage_1.use[i__ - 1]) {
		    atoms_1.x[i__ - 1] = scale * atoms_1.x[i__ - 1];
		    atoms_1.y[i__ - 1] = scale * atoms_1.y[i__ - 1];
		    atoms_1.z__[i__ - 1] = scale * atoms_1.z__[i__ - 1];
		}
	    }

/*     couple to pressure bath via center of mass of rigid bodies */

	} else {
	    scale += -1.;
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		start = igrp_ref(1, i__);
		stop = igrp_ref(2, i__);
		xcm = 0.;
		ycm = 0.;
		zcm = 0.;
		i__2 = stop;
		for (j = start; j <= i__2; ++j) {
		    k = group_1.kgrp[j - 1];
		    weigh = atmtyp_1.mass[k - 1];
		    xcm += atoms_1.x[k - 1] * weigh;
		    ycm += atoms_1.y[k - 1] * weigh;
		    zcm += atoms_1.z__[k - 1] * weigh;
		}
		xmove = scale * xcm / group_1.grpmass[i__ - 1];
		ymove = scale * ycm / group_1.grpmass[i__ - 1];
		zmove = scale * zcm / group_1.grpmass[i__ - 1];
		i__2 = stop;
		for (j = start; j <= i__2; ++j) {
		    k = group_1.kgrp[j - 1];
		    atoms_1.x[k - 1] += xmove;
		    atoms_1.y[k - 1] += ymove;
		    atoms_1.z__[k - 1] += zmove;
		}
	    }
	}

/*     find the anisotropic scale factors for constant pressure */

    } else {
	scale = *dt * bath_1.compress / (bath_1.taupres * 3.);
	for (i__ = 1; i__ <= 3; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		if (j == i__) {
		    ascale_ref(j, i__) = scale * (stress_ref(i__, i__) - 
			    bath_1.atmsph) + 1.;
		} else {
		    ascale_ref(j, i__) = scale * stress_ref(j, i__);
		}
	    }
	}

/*     modify the current periodic box dimension values */

	temp_ref(1, 1) = boxes_1.xbox;
	temp_ref(2, 1) = 0.;
	temp_ref(3, 1) = 0.;
	temp_ref(1, 2) = boxes_1.ybox * boxes_1.gamma_cos__;
	temp_ref(2, 2) = boxes_1.ybox * boxes_1.gamma_sin__;
	temp_ref(3, 2) = 0.;
	temp_ref(1, 3) = boxes_1.zbox * boxes_1.beta_cos__;
	temp_ref(2, 3) = boxes_1.zbox * boxes_1.beta_term__;
	temp_ref(3, 3) = boxes_1.zbox * boxes_1.gamma_term__;
	for (i__ = 1; i__ <= 3; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		hbox_ref(j, i__) = 0.;
		for (k = 1; k <= 3; ++k) {
		    hbox_ref(j, i__) = hbox_ref(j, i__) + ascale_ref(j, k) * 
			    temp_ref(k, i__);
		}
	    }
	}
/* Computing 2nd power */
	d__1 = hbox_ref(1, 1);
/* Computing 2nd power */
	d__2 = hbox_ref(2, 1);
/* Computing 2nd power */
	d__3 = hbox_ref(3, 1);
	boxes_1.xbox = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
	d__1 = hbox_ref(1, 2);
/* Computing 2nd power */
	d__2 = hbox_ref(2, 2);
/* Computing 2nd power */
	d__3 = hbox_ref(3, 2);
	boxes_1.ybox = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
	d__1 = hbox_ref(1, 3);
/* Computing 2nd power */
	d__2 = hbox_ref(2, 3);
/* Computing 2nd power */
	d__3 = hbox_ref(3, 3);
	boxes_1.zbox = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	if (boxes_1.monoclinic) {
	    cosine = (hbox_ref(1, 1) * hbox_ref(1, 3) + hbox_ref(2, 1) * 
		    hbox_ref(2, 3) + hbox_ref(3, 1) * hbox_ref(3, 3)) / (
		    boxes_1.xbox * boxes_1.zbox);
	    boxes_1.beta = acos(cosine) * 57.29577951308232088;
	} else if (boxes_1.triclinic) {
	    cosine = (hbox_ref(1, 2) * hbox_ref(1, 3) + hbox_ref(2, 2) * 
		    hbox_ref(2, 3) + hbox_ref(3, 2) * hbox_ref(3, 3)) / (
		    boxes_1.ybox * boxes_1.zbox);
	    boxes_1.alpha = acos(cosine) * 57.29577951308232088;
	    cosine = (hbox_ref(1, 1) * hbox_ref(1, 3) + hbox_ref(2, 1) * 
		    hbox_ref(2, 3) + hbox_ref(3, 1) * hbox_ref(3, 3)) / (
		    boxes_1.xbox * boxes_1.zbox);
	    boxes_1.beta = acos(cosine) * 57.29577951308232088;
	    cosine = (hbox_ref(1, 1) * hbox_ref(1, 2) + hbox_ref(2, 1) * 
		    hbox_ref(2, 2) + hbox_ref(3, 1) * hbox_ref(3, 2)) / (
		    boxes_1.xbox * boxes_1.ybox);
	    boxes_1.gamma = acos(cosine) * 57.29577951308232088;
	}

/*     propagate the new box dimensions to other lattice values */

	lattice_();

/*     couple to pressure bath via atom scaling in Cartesian space */

	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) != 
		0) {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (usage_1.use[i__ - 1]) {
		    atoms_1.x[i__ - 1] = ascale_ref(1, 1) * atoms_1.x[i__ - 1]
			     + ascale_ref(1, 2) * atoms_1.y[i__ - 1] + 
			    ascale_ref(1, 3) * atoms_1.z__[i__ - 1];
		    atoms_1.y[i__ - 1] = ascale_ref(2, 1) * atoms_1.x[i__ - 1]
			     + ascale_ref(2, 2) * atoms_1.y[i__ - 1] + 
			    ascale_ref(2, 3) * atoms_1.z__[i__ - 1];
		    atoms_1.z__[i__ - 1] = ascale_ref(3, 1) * atoms_1.x[i__ - 
			    1] + ascale_ref(3, 2) * atoms_1.y[i__ - 1] + 
			    ascale_ref(3, 3) * atoms_1.z__[i__ - 1];
		}
	    }

/*     couple to pressure bath via center of mass of rigid bodies */

	} else {
	    ascale_ref(1, 1) = ascale_ref(1, 1) - 1.;
	    ascale_ref(2, 2) = ascale_ref(2, 2) - 1.;
	    ascale_ref(3, 3) = ascale_ref(3, 3) - 1.;
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		start = igrp_ref(1, i__);
		stop = igrp_ref(2, i__);
		xcm = 0.;
		ycm = 0.;
		zcm = 0.;
		i__2 = stop;
		for (j = start; j <= i__2; ++j) {
		    k = group_1.kgrp[j - 1];
		    weigh = atmtyp_1.mass[k - 1];
		    xcm += atoms_1.x[k - 1] * weigh;
		    ycm = xcm + atoms_1.y[k - 1] * weigh;
		    zcm = xcm + atoms_1.z__[k - 1] * weigh;
		}
		xcm /= group_1.grpmass[i__ - 1];
		ycm /= group_1.grpmass[i__ - 1];
		zcm /= group_1.grpmass[i__ - 1];
		xmove = ascale_ref(1, 1) * xcm + ascale_ref(1, 2) * ycm + 
			ascale_ref(1, 3) * zcm;
		ymove = ascale_ref(2, 1) * xcm + ascale_ref(2, 2) * ycm + 
			ascale_ref(2, 3) * zcm;
		zmove = ascale_ref(3, 1) * xcm + ascale_ref(3, 2) * ycm + 
			ascale_ref(3, 3) * zcm;
		i__2 = stop;
		for (j = start; j <= i__2; ++j) {
		    k = group_1.kgrp[j - 1];
		    atoms_1.x[k - 1] += xmove;
		    atoms_1.y[k - 1] += ymove;
		    atoms_1.z__[k - 1] += zmove;
		}
	    }
	}
    }
    return 0;
} /* pscale_ */

#undef stress_ref
#undef ascale_ref
#undef temp_ref
#undef igrp_ref
#undef hbox_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine pmonte  --  Monte Carlo barostat trial moves  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "pmonte" implements a Monte Carlo barostat via random trial */
/*     changes in the periodic box volume and shape */

/*     literature references: */

/*     D. Frenkel and B. Smit, "Understanding Molecular Simulation, */
/*     2nd Edition", Academic Press, San Diego, CA, 2002; Section 5.4 */

/*     original version written by Alan Grossfield, January 2004; */
/*     this version only implements isotropic volume changes */


/* Subroutine */ int pmonte_(doublereal *epot, doublereal *temp)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal de, dv, kt, xcm, ycm, zcm, lnv, diff, enew, step, vold, 
	    term, xold[25000], yold[25000];
    static integer stop;
    static doublereal zold[25000], scale, weigh, third;
    static integer start;
    static doublereal xmove, ymove, zmove;
    extern doublereal random_(void), energy_(void);
    extern /* Subroutine */ int lattice_(void);
    static doublereal xboxold, yboxold, zboxold, expterm;


#define imol_ref(a_1,a_2) molcul_1.imol[(a_2)*2 + a_1 - 3]
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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  boxes.i  --  parameters for periodic boundary conditions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     xbox        length of a-axis of periodic box in Angstroms */
/*     ybox        length of b-axis of periodic box in Angstroms */
/*     zbox        length of c-axis of periodic box in Angstroms */
/*     alpha       angle between b- and c-axes of box in degrees */
/*     beta        angle between a- and c-axes of box in degrees */
/*     gamma       angle between a- and b-axes of box in degrees */
/*     xbox2       half of the a-axis length of periodic box */
/*     ybox2       half of the b-axis length of periodic box */
/*     zbox2       half of the c-axis length of periodic box */
/*     box34       three-fourths axis length of truncated octahedron */
/*     lvec        real space lattice vectors as matrix rows */
/*     recip       reciprocal lattice vectors as matrix columns */
/*     volbox      volume in Ang**3 of the periodic box */
/*     beta_sin    sine of the beta periodic box angle */
/*     beta_cos    cosine of the beta periodic box angle */
/*     gamma_sin   sine of the gamma periodic box angle */
/*     gamma_cos   cosine of the gamma periodic box angle */
/*     beta_term   term used in generating triclinic box */
/*     gamma_term  term used in generating triclinic box */
/*     orthogonal  flag to mark periodic box as orthogonal */
/*     monoclinic  flag to mark periodic box as monoclinic */
/*     triclinic   flag to mark periodic box as triclinic */
/*     octahedron  flag to mark box as truncated octahedron */
/*     spacegrp    space group symbol for the unitcell type */




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
/*     ##  molcul.i  --  individual molecules within current system  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     molmass   molecular weight for each molecule in the system */
/*     totmass   total weight of all the molecules in the system */
/*     nmol      total number of separate molecules in the system */
/*     kmol      contiguous list of the atoms in each molecule */
/*     imol      first and last atom of each molecule in the list */
/*     molcule   number of the molecule to which each atom belongs */




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




/*     make volume move, save old box size and coordinates */

    if (random_() < 1. / (doublereal) bath_1.voltrial) {
	xboxold = boxes_1.xbox;
	yboxold = boxes_1.ybox;
	zboxold = boxes_1.zbox;
	vold = boxes_1.volbox;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xold[i__ - 1] = atoms_1.x[i__ - 1];
	    yold[i__ - 1] = atoms_1.y[i__ - 1];
	    zold[i__ - 1] = atoms_1.z__[i__ - 1];
	}

/*     get scale factor that reflects the chosen volume change */

	step = bath_1.volmove * (random_() * 2. - 1.);
	boxes_1.volbox += step;
	third = .33333333333333331;
	d__1 = boxes_1.volbox / vold;
	scale = pow_dd(&d__1, &third);
	if (boxes_1.monoclinic) {
	    term = 1. / boxes_1.beta_sin__;
	    scale = (scale - 1.) * pow_dd(&term, &third) + 1.;
	} else if (boxes_1.triclinic) {
	    term = 1. / (boxes_1.gamma_sin__ * boxes_1.gamma_term__);
	    scale = (scale - 1.) * pow_dd(&term, &third) + 1.;
	} else if (boxes_1.octahedron) {
	    term = 2.;
	    scale = (scale - 1.) * pow_dd(&term, &third) + 1.;
	}

/*     set the new box dimensions and other lattice values */

	boxes_1.xbox *= scale;
	boxes_1.ybox *= scale;
	boxes_1.zbox *= scale;
	lattice_();

/*     scale the coordinates by groups, molecules or atoms */

	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    diff = scale - 1.;
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		xcm = 0.;
		ycm = 0.;
		zcm = 0.;
		start = igrp_ref(1, i__);
		stop = igrp_ref(2, i__);
		i__2 = stop;
		for (j = start; j <= i__2; ++j) {
		    k = group_1.kgrp[j - 1];
		    weigh = atmtyp_1.mass[k - 1];
		    xcm += atoms_1.x[k - 1] * weigh;
		    ycm += atoms_1.y[k - 1] * weigh;
		    zcm += atoms_1.z__[k - 1] * weigh;
		}
		xmove = diff * xcm / group_1.grpmass[i__ - 1];
		ymove = diff * ycm / group_1.grpmass[i__ - 1];
		zmove = diff * zcm / group_1.grpmass[i__ - 1];
		i__2 = stop;
		for (j = start; j <= i__2; ++j) {
		    k = group_1.kgrp[j - 1];
		    atoms_1.x[k - 1] += xmove;
		    atoms_1.y[k - 1] += ymove;
		    atoms_1.z__[k - 1] += zmove;
		}
	    }
	} else if (s_cmp(bath_1.volscale, "MOLECULAR", (ftnlen)9, (ftnlen)9) 
		== 0) {
	    diff = scale - 1.;
	    i__1 = molcul_1.nmol;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		xcm = 0.;
		ycm = 0.;
		zcm = 0.;
		start = imol_ref(1, i__);
		stop = imol_ref(2, i__);
		i__2 = stop;
		for (j = start; j <= i__2; ++j) {
		    k = molcul_1.kmol[j - 1];
		    weigh = atmtyp_1.mass[k - 1];
		    xcm += atoms_1.x[k - 1] * weigh;
		    ycm += atoms_1.y[k - 1] * weigh;
		    zcm += atoms_1.z__[k - 1] * weigh;
		}
		xmove = diff * xcm / molcul_1.molmass[i__ - 1];
		ymove = diff * ycm / molcul_1.molmass[i__ - 1];
		zmove = diff * zcm / molcul_1.molmass[i__ - 1];
		i__2 = stop;
		for (j = start; j <= i__2; ++j) {
		    k = molcul_1.kmol[j - 1];
		    if (usage_1.use[k - 1]) {
			atoms_1.x[k - 1] += xmove;
			atoms_1.y[k - 1] += ymove;
			atoms_1.z__[k - 1] += zmove;
		    }
		}
	    }
	} else {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (usage_1.use[i__ - 1]) {
		    atoms_1.x[i__ - 1] = scale * atoms_1.x[i__ - 1];
		    atoms_1.y[i__ - 1] = scale * atoms_1.y[i__ - 1];
		    atoms_1.z__[i__ - 1] = scale * atoms_1.z__[i__ - 1];
		}
	    }
	}

/*     find the energy change and PV term for the trial move */

	enew = energy_();
	de = enew - *epot;
	dv = bath_1.atmsph * (boxes_1.volbox - vold) / 68568.4112;

/*     set the entropy of mixing term based on system type */

	if (bath_1.isothermal) {
	    kt = bath_1.kelvin0 * .0019872066;
	} else {
	    kt = *temp * .0019872066;
	}
	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    lnv = (doublereal) group_1.ngrp * kt * log(boxes_1.volbox / vold);
	} else if (s_cmp(bath_1.volscale, "MOLECULAR", (ftnlen)9, (ftnlen)9) 
		== 0) {
	    lnv = (doublereal) molcul_1.nmol * kt * log(boxes_1.volbox / vold)
		    ;
	} else {
	    lnv = (doublereal) usage_1.nuse * kt * log(boxes_1.volbox / vold);
	}

/*     acceptance ratio from energy change, PV and mixing term */

	term = -(de + dv - lnv) / kt;
	expterm = exp(term);

/*     reject the step, restore old box size and coordinates */

	if (random_() > expterm) {
	    boxes_1.xbox = xboxold;
	    boxes_1.ybox = yboxold;
	    boxes_1.zbox = zboxold;
	    lattice_();
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		atoms_1.x[i__ - 1] = xold[i__ - 1];
		atoms_1.y[i__ - 1] = yold[i__ - 1];
		atoms_1.z__[i__ - 1] = zold[i__ - 1];
	    }
	}
    }
    return 0;
} /* pmonte_ */

#undef igrp_ref
#undef imol_ref


