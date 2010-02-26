/* mdstat.f -- translated by f2c (version 20050501).
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
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

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
    doublereal vdwcut, chgcut, dplcut, mpolecut, vdwtaper, chgtaper, dpltaper,
	     mpoletaper, ewaldcut;
    logical use_ewald__, use_lights__, use_list__, use_vlist__, use_clist__, 
	    use_mlist__;
} cutoff_;

#define cutoff_1 cutoff_

struct {
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    doublereal einter;
} inter_;

#define inter_1 inter_

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
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal m2[25000], deform, difft, diffv, diffc;
    logical use_smooth__, use_dem__, use_gda__, use_tophat__, use_stophat__;
} warp_;

#define warp_1 warp_

/* Table of constant values */

static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine mdstat  --  compute averages over a trajectory  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "mdstat" is called at each molecular dynamics time step to */
/*     form statistics on various average values and fluctuations, */
/*     and to periodically save the state of the trajectory */


/* Subroutine */ int mdstat_(integer *istep, doublereal *dt, doublereal *etot,
	 doublereal *epot, doublereal *ekin, doublereal *temp, doublereal *
	pres)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002    MD Step      E Total   E Potential"
	    "\002,\002   E Kinetic       Temp       Pres\002,/)";
    static char fmt_20[] = "(/,\002    MD Step      E Total   E Potential"
	    "\002,\002   E Kinetic       Temp\002,/)";
    static char fmt_30[] = "(i10,2f14.4,f12.4,2f11.2)";
    static char fmt_40[] = "(i10,2f14.4,f12.4,f11.2)";
    static char fmt_50[] = "(/,\002 Average Values for the last\002,i6,\002 "
	    "out of\002,i9,\002 Dynamics Steps\002)";
    static char fmt_60[] = "(/,\002 Simulation Time\002,5x,f15.4,\002 Picose"
	    "cond\002)";
    static char fmt_70[] = "(\002 Total Energy\002,8x,f15.4,\002 Kcal/mol"
	    "e\002,3x,\002(+/-\002,f9.4,\002)\002)";
    static char fmt_80[] = "(\002 Potential Energy\002,4x,f15.4,\002 Kcal/mo"
	    "le\002,3x,\002(+/-\002,f9.4,\002)\002)";
    static char fmt_90[] = "(\002 Kinetic Energy\002,6x,f15.4,\002 Kcal/mol"
	    "e\002,3x,\002(+/-\002,f9.4,\002)\002)";
    static char fmt_100[] = "(\002 Intermolecular\002,6x,f15.4,\002 Kcal/m"
	    "ole\002,3x,\002(+/-\002,f9.4,\002)\002)";
    static char fmt_110[] = "(\002 Temperature\002,9x,f15.2,\002 Kelvin\002,"
	    "6x,\002(+/-\002,f9.2,\002)\002)";
    static char fmt_120[] = "(\002 Pressure\002,12x,f15.2,\002 Atmosphere"
	    "\002,2x,\002(+/-\002,f9.2,\002)\002)";
    static char fmt_130[] = "(\002 Density\002,13x,f15.4,\002 Grams/cc\002,4"
	    "x,\002(+/-\002,f9.4,\002)\002)";
    static char fmt_140[] = "(\002 Deformation\002,9x,f15.3,\002 Sqr Angs"
	    "\002)";

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static doublereal ekin_ave__, dens_ave__, eint_ave__, temp_ave__, 
	    epot_ave__, pres_ave__, etot_ave__, kinfluct, eint_sum__, 
	    ekin_sum__, dens_sum__, intfluct, temp_sum__, epot_sum__, 
	    etot_sum__, potfluct, pres_sum__, ekin2_ave__, dens2_ave__, 
	    eint2_ave__, temp2_ave__, epot2_ave__, pres2_ave__, etot2_ave__, 
	    ekin2_sum__, dens2_sum__, kinfluct2, eint2_sum__, intfluct2, 
	    epot2_sum__, temp2_sum__, pres2_sum__, etot2_sum__, potfluct2, 
	    fluctuate, fluctuate2, pico, dens, dfluct, tfluct, pfluct, 
	    dfluct2, pfluct2, tfluct2;
    static integer modstep;

    /* Fortran I/O blocks */
    static cilist io___16 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___52 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___58 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___59 = { 0, 0, 0, fmt_140, 0 };




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cutoff.i  --  cutoff distances for energy interactions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     vdwcut      cutoff distance for van der Waals interactions */
/*     chgcut      cutoff distance for charge-charge interactions */
/*     dplcut      cutoff distance for dipole-dipole interactions */
/*     mpolecut    cutoff distance for atomic multipole interactions */
/*     vdwtaper    distance at which van der Waals switching begins */
/*     chgtaper    distance at which charge-charge switching begins */
/*     dpltaper    distance at which dipole-dipole switching begins */
/*     mpoletaper  distance at which atomic multipole switching begins */
/*     ewaldcut    cutoff distance for direct space Ewald summation */
/*     use_ewald   logical flag governing use of Ewald summation */
/*     use_lights  logical flag governing use of method of lights */
/*     use_list    logical flag governing use of any neighbor lists */
/*     use_vlist   logical flag governing use of vdw neighbor lists */
/*     use_clist   logical flag governing use of charge neighbor lists */
/*     use_mlist   logical flag governing use of multipole neighbor lists */




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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  inter.i  --  sum of intermolecular energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     einter   total intermolecular potential energy */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  warp.i  --  parameters for potential surface smoothing  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     m2           second moment of the GDA gaussian for each atom */
/*     deform       value of the smoothing deformation parameter */
/*     difft        diffusion coefficient for torsional potential */
/*     diffv        diffusion coefficient for van der Waals potential */
/*     diffc        diffusion coefficient for charge-charge potential */
/*     use_smooth   flag to use a potential energy smoothing method */
/*     use_dem      flag to use diffusion equation method potential */
/*     use_gda      flag to use gaussian density annealing potential */
/*     use_tophat   flag to use analytical tophat smoothed potential */
/*     use_stophat  flag to use shifted tophat smoothed potential */




/*     set number of steps for block averages of properties */

    modstep = *istep % inform_1.iprint;

/*     zero out summation variables for new averaging period */

    if (modstep == 1 || inform_1.iprint == 1) {
	etot_sum__ = 0.;
	etot2_sum__ = 0.;
	epot_sum__ = 0.;
	epot2_sum__ = 0.;
	ekin_sum__ = 0.;
	ekin2_sum__ = 0.;
	eint_sum__ = 0.;
	eint2_sum__ = 0.;
	temp_sum__ = 0.;
	temp2_sum__ = 0.;
	pres_sum__ = 0.;
	pres2_sum__ = 0.;
	dens_sum__ = 0.;
	dens2_sum__ = 0.;
    }

/*     print energy, temperature and pressure for current step */

    if (inform_1.verbose) {
	if (modstep == 1) {
	    if (bound_1.use_bounds__ && s_cmp(mdstuf_1.integrate, "STOCHASTIC"
		    , (ftnlen)10, (ftnlen)10) != 0) {
		io___16.ciunit = iounit_1.iout;
		s_wsfe(&io___16);
		e_wsfe();
	    } else {
		io___17.ciunit = iounit_1.iout;
		s_wsfe(&io___17);
		e_wsfe();
	    }
	}
	if (bound_1.use_bounds__ && s_cmp(mdstuf_1.integrate, "STOCHASTIC", (
		ftnlen)10, (ftnlen)10) != 0) {
	    io___18.ciunit = iounit_1.iout;
	    s_wsfe(&io___18);
	    do_fio(&c__1, (char *)&(*istep), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*etot), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*epot), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ekin), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*temp), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*pres), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___19.ciunit = iounit_1.iout;
	    s_wsfe(&io___19);
	    do_fio(&c__1, (char *)&(*istep), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*etot), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*epot), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ekin), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*temp), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     print header for the averages over a group of recent steps */

    if (modstep == 0) {
	pico = (doublereal) (*istep) * *dt;
	io___21.ciunit = iounit_1.iout;
	s_wsfe(&io___21);
	do_fio(&c__1, (char *)&inform_1.iprint, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*istep), (ftnlen)sizeof(integer));
	e_wsfe();
	io___22.ciunit = iounit_1.iout;
	s_wsfe(&io___22);
	do_fio(&c__1, (char *)&pico, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     compute total energy and fluctuation for recent steps */

    etot_sum__ += *etot;
/* Computing 2nd power */
    d__1 = *etot;
    etot2_sum__ += d__1 * d__1;
    if (modstep == 0) {
	etot_ave__ = etot_sum__ / (doublereal) inform_1.iprint;
	etot2_ave__ = etot2_sum__ / (doublereal) inform_1.iprint;
/* Computing 2nd power */
	d__1 = etot_ave__;
	fluctuate2 = etot2_ave__ - d__1 * d__1;
	if (fluctuate2 > 0.) {
	    fluctuate = sqrt(fluctuate2);
	} else {
	    fluctuate = 0.;
	}
	io___27.ciunit = iounit_1.iout;
	s_wsfe(&io___27);
	do_fio(&c__1, (char *)&etot_ave__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&fluctuate, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     compute average potential energy and its fluctuation */

    epot_sum__ += *epot;
/* Computing 2nd power */
    d__1 = *epot;
    epot2_sum__ += d__1 * d__1;
    if (modstep == 0) {
	epot_ave__ = epot_sum__ / (doublereal) inform_1.iprint;
	epot2_ave__ = epot2_sum__ / (doublereal) inform_1.iprint;
/* Computing 2nd power */
	d__1 = epot_ave__;
	potfluct2 = epot2_ave__ - d__1 * d__1;
	if (potfluct2 > 0.) {
	    potfluct = sqrt(potfluct2);
	} else {
	    potfluct = 0.;
	}
	io___32.ciunit = iounit_1.iout;
	s_wsfe(&io___32);
	do_fio(&c__1, (char *)&epot_ave__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&potfluct, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     compute average kinetic energy and its fluctuation */

    ekin_sum__ += *ekin;
/* Computing 2nd power */
    d__1 = *ekin;
    ekin2_sum__ += d__1 * d__1;
    if (modstep == 0) {
	ekin_ave__ = ekin_sum__ / (doublereal) inform_1.iprint;
	ekin2_ave__ = ekin2_sum__ / (doublereal) inform_1.iprint;
/* Computing 2nd power */
	d__1 = ekin_ave__;
	kinfluct2 = ekin2_ave__ - d__1 * d__1;
	if (kinfluct2 > 0.) {
	    kinfluct = sqrt(kinfluct2);
	} else {
	    kinfluct = 0.;
	}
	io___37.ciunit = iounit_1.iout;
	s_wsfe(&io___37);
	do_fio(&c__1, (char *)&ekin_ave__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&kinfluct, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     compute average intermolecular energy and its fluctuation */

    if (molcul_1.nmol != 1 && molcul_1.nmol != atoms_1.n && ! 
	    cutoff_1.use_ewald__) {
	eint_sum__ += inter_1.einter;
/* Computing 2nd power */
	d__1 = inter_1.einter;
	eint2_sum__ += d__1 * d__1;
	if (modstep == 0) {
	    eint_ave__ = eint_sum__ / (doublereal) inform_1.iprint;
	    eint2_ave__ = eint2_sum__ / (doublereal) inform_1.iprint;
/* Computing 2nd power */
	    d__1 = eint_ave__;
	    intfluct2 = eint2_ave__ - d__1 * d__1;
	    if (intfluct2 > 0.) {
		intfluct = sqrt(intfluct2);
	    } else {
		intfluct = 0.;
	    }
	    io___42.ciunit = iounit_1.iout;
	    s_wsfe(&io___42);
	    do_fio(&c__1, (char *)&eint_ave__, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&intfluct, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     compute the average temperature and its fluctuation */

    temp_sum__ += *temp;
/* Computing 2nd power */
    d__1 = *temp;
    temp2_sum__ += d__1 * d__1;
    if (modstep == 0) {
	temp_ave__ = temp_sum__ / (doublereal) inform_1.iprint;
	temp2_ave__ = temp2_sum__ / (doublereal) inform_1.iprint;
/* Computing 2nd power */
	d__1 = temp_ave__;
	tfluct2 = temp2_ave__ - d__1 * d__1;
	if (tfluct2 > 0.) {
	    tfluct = sqrt(tfluct2);
	} else {
	    tfluct = 0.;
	}
	io___47.ciunit = iounit_1.iout;
	s_wsfe(&io___47);
	do_fio(&c__1, (char *)&temp_ave__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&tfluct, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     compute the average pressure and its fluctuation */

    if (bound_1.use_bounds__) {
	pres_sum__ += *pres;
/* Computing 2nd power */
	d__1 = *pres;
	pres2_sum__ += d__1 * d__1;
	if (modstep == 0) {
	    pres_ave__ = pres_sum__ / (doublereal) inform_1.iprint;
	    pres2_ave__ = pres2_sum__ / (doublereal) inform_1.iprint;
/* Computing 2nd power */
	    d__1 = pres_ave__;
	    pfluct2 = pres2_ave__ - d__1 * d__1;
	    if (pfluct2 > 0.) {
		pfluct = sqrt(pfluct2);
	    } else {
		pfluct = 0.;
	    }
	    io___52.ciunit = iounit_1.iout;
	    s_wsfe(&io___52);
	    do_fio(&c__1, (char *)&pres_ave__, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&pfluct, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}

/*     compute the average density and its fluctuation */

	dens = 1e24 / boxes_1.volbox * (molcul_1.totmass / 6.02214179e23);
	dens_sum__ += dens;
/* Computing 2nd power */
	d__1 = dens;
	dens2_sum__ += d__1 * d__1;
	if (modstep == 0) {
	    dens_ave__ = dens_sum__ / (doublereal) inform_1.iprint;
	    dens2_ave__ = dens2_sum__ / (doublereal) inform_1.iprint;
/* Computing 2nd power */
	    d__1 = dens_ave__;
	    dfluct2 = dens2_ave__ - d__1 * d__1;
	    if (dfluct2 > 0.) {
		dfluct = sqrt(dfluct2);
	    } else {
		dfluct = 0.;
	    }
	    io___58.ciunit = iounit_1.iout;
	    s_wsfe(&io___58);
	    do_fio(&c__1, (char *)&dens_ave__, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&dfluct, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     note deformation value for potential energy smoothing */

    if (warp_1.use_smooth__) {
	if (modstep == 0) {
	    io___59.ciunit = iounit_1.iout;
	    s_wsfe(&io___59);
	    do_fio(&c__1, (char *)&warp_1.deform, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* mdstat_ */

