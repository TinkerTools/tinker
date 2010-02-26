/* dynamic.f -- translated by f2c (version 20050501).
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
    doublereal bk[50000], bl[50000];
    integer nbond, ibnd[100000]	/* was [2][50000] */;
} bond_;

#define bond_1 bond_

struct {
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

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
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

struct {
    doublereal rsolv[25000], asolv[25000], rborn[25000], drb[25000], drbp[
	    25000], drobc[25000], doffset, p1, p2, p3, p4, p5, gpol[25000], 
	    shct[25000], aobc[25000], bobc[25000], gobc[25000], vsolv[25000], 
	    wace[1000000]	/* was [1000][1000] */, s2ace[1000000]	/* 
	    was [1000][1000] */, uace[1000000]	/* was [1000][1000] */;
    char solvtyp[8], borntyp[8];
} solute_;

#define solute_1 solute_

struct {
    doublereal friction, fgamma[25000];
    logical use_sdarea__;
} stodyn_;

#define stodyn_1 stodyn_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  program dynamic  --  run molecular or stochastic dynamics  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "dynamic" computes a molecular dynamics trajectory in any of */
/*     several statistical mechanical ensembles with optional periodic */
/*     boundaries and optional coupling to temperature and pressure baths */
/*     alternatively a stochastic dynamics trajectory can be generated */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Enter the Number of Dynamics Steps to b"
	    "e\002,\002 Taken :  \002,$)";
    static char fmt_30[] = "(i10)";
    static char fmt_50[] = "(/,\002 Enter the Time Step Length in Femtosecon"
	    "ds\002,\002 [1.0] :  \002,$)";
    static char fmt_60[] = "(f20.0)";
    static char fmt_90[] = "(/,\002 Enter Time between Dumps in Picosecond"
	    "s\002,\002 [0.1] :  \002,$)";
    static char fmt_100[] = "(f20.0)";
    static char fmt_130[] = "(/,\002 Available Statistical Mechanical Ensemb"
	    "les :\002,//,4x,\002(1) Microcanonical (NVE)\002,/,4x,\002(2) Ca"
	    "nonical (NVT)\002,/,4x,\002(3) Isoenthalpic-Isobaric (NPH)\002,/"
	    ",4x,\002(4) Isothermal-Isobaric (NPT)\002,//,\002 Enter the Numb"
	    "er of the Desired Choice\002,\002 [1] :  \002,$)";
    static char fmt_140[] = "(i10)";
    static char fmt_170[] = "(/,\002 Enter the Desired Temperature in Degr"
	    "ees\002,\002 K [298] :  \002,$)";
    static char fmt_180[] = "(f20.0)";
    static char fmt_210[] = "(/,\002 Enter the Desired Pressure in Atm\002"
	    ",\002 [1.0] :  \002,$)";
    static char fmt_220[] = "(f20.0)";
    static char fmt_250[] = "(/,\002 Available Simulation Control Modes :"
	    "\002,//,4x,\002(1) Constant Total Energy Value (E)\002,/,4x,\002"
	    "(2) Constant Temperature via Thermostat (T)\002,//,\002 Enter th"
	    "e Number of the Desired Choice\002,\002 [1] :  \002,$)";
    static char fmt_260[] = "(i10)";
    static char fmt_290[] = "(/,\002 Enter the Desired Temperature in Degr"
	    "ees\002,\002 K [298] :  \002,$)";
    static char fmt_300[] = "(f20.0)";
    static char fmt_320[] = "(/,\002 Molecular Dynamics Trajectory via\002"
	    ",\002 Velocity Verlet Algorithm\002)";
    static char fmt_330[] = "(/,\002 Stochastic Dynamics Trajectory via\002"
	    ",\002 Velocity Verlet Algorithm\002)";
    static char fmt_340[] = "(/,\002 Molecular Dynamics Trajectory via\002"
	    ",\002 Rigid Body Algorithm\002)";
    static char fmt_350[] = "(/,\002 Molecular Dynamics Trajectory via\002"
	    ",\002 Modified Beeman Algorithm\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void), i_dnnt(
	    doublereal *), s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    static doublereal dt;
    static integer mode;
    extern /* Subroutine */ int final_(void);
    static integer istep, nstep;
    static logical exist, query;
    extern /* Subroutine */ int beeman_(integer *, doublereal *), mdinit_(
	    void);
    static doublereal dtdump;
    static char string[120];
    extern /* Subroutine */ int verlet_(integer *, doublereal *), sdstep_(
	    integer *, doublereal *), mdrest_(void), getxyz_(void), initial_(
	    void), shakeup_(void), nextarg_(char *, logical *, ftnlen);
    static integer modstep;
    extern /* Subroutine */ int rgdstep_(integer *, doublereal *);

    /* Fortran I/O blocks */
    static icilist io___4 = { 1, string, 1, 0, 120, 1 };
    static cilist io___6 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_30, 0 };
    static icilist io___9 = { 1, string, 1, 0, 120, 1 };
    static cilist io___10 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___11 = { 1, 0, 0, fmt_60, 0 };
    static icilist io___13 = { 1, string, 1, 0, 120, 1 };
    static cilist io___14 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___15 = { 1, 0, 0, fmt_100, 0 };
    static icilist io___17 = { 1, string, 1, 0, 120, 1 };
    static cilist io___18 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___19 = { 1, 0, 0, fmt_140, 0 };
    static icilist io___20 = { 1, string, 1, 0, 120, 1 };
    static cilist io___21 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___22 = { 1, 0, 0, fmt_180, 0 };
    static icilist io___23 = { 1, string, 1, 0, 120, 1 };
    static cilist io___24 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___25 = { 1, 0, 0, fmt_220, 0 };
    static icilist io___26 = { 1, string, 1, 0, 120, 1 };
    static cilist io___27 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___28 = { 1, 0, 0, fmt_260, 0 };
    static icilist io___29 = { 1, string, 1, 0, 120, 1 };
    static cilist io___30 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___31 = { 1, 0, 0, fmt_300, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_320, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_330, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_340, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_350, 0 };




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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  bond.i  --  covalent bonds in the current structure  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     bk      bond stretch force constants (kcal/mole/Ang**2) */
/*     bl      ideal bond length values in Angstroms */
/*     nbond   total number of bond stretches in the system */
/*     ibnd    numbers of the atoms in each bond stretch */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  potent.i  --  usage of each potential energy component  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     use_bond    logical flag governing use of bond stretch potential */
/*     use_angle   logical flag governing use of angle bend potential */
/*     use_strbnd  logical flag governing use of stretch-bend potential */
/*     use_urey    logical flag governing use of Urey-Bradley potential */
/*     use_angang  logical flag governing use of angle-angle cross term */
/*     use_opbend  logical flag governing use of out-of-plane bend term */
/*     use_opdist  logical flag governing use of out-of-plane distance */
/*     use_improp  logical flag governing use of improper dihedral term */
/*     use_imptor  logical flag governing use of improper torsion term */
/*     use_tors    logical flag governing use of torsional potential */
/*     use_pitors  logical flag governing use of pi-orbital torsion term */
/*     use_strtor  logical flag governing use of stretch-torsion term */
/*     use_tortor  logical flag governing use of torsion-torsion term */
/*     use_vdw     logical flag governing use of vdw der Waals potential */
/*     use_charge  logical flag governing use of charge-charge potential */
/*     use_chgdpl  logical flag governing use of charge-dipole potential */
/*     use_dipole  logical flag governing use of dipole-dipole potential */
/*     use_mpole   logical flag governing use of multipole potential */
/*     use_polar   logical flag governing use of polarization term */
/*     use_rxnfld  logical flag governing use of reaction field term */
/*     use_solv    logical flag governing use of continuum solvation */
/*     use_metal   logical flag governing use of ligand field term */
/*     use_geom    logical flag governing use of geometric restraints */
/*     use_extra   logical flag governing use of extra potential term */
/*     use_born    logical flag governing use of Born radii values */
/*     use_orbit   logical flag governing use of pisystem computation */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  solute.i  --  parameters for continuum solvation models  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     rsolv     atomic radius of each atom for continuum solvation */
/*     asolv     atomic surface area solvation parameters */
/*     rborn     Born radius of each atom for GB/SA solvation */
/*     drb       solvation derivatives with respect to Born radii */
/*     drbp      GK polarization derivatives with respect to Born radii */
/*     drobc     chain rule term for Onufriev-Bashford-Case radii */
/*     doffset   dielectric offset to continuum solvation atomic radii */
/*     p1        single-atom scale factor for analytical Still radii */
/*     p2        1-2 interaction scale factor for analytical Still radii */
/*     p3        1-3 interaction scale factor for analytical Still radii */
/*     p4        nonbonded scale factor for analytical Still radii */
/*     p5        soft cutoff parameter for analytical Still radii */
/*     gpol      polarization self-energy values for each atom */
/*     shct      overlap scale factors for Hawkins-Cramer-Truhlar radii */
/*     aobc      alpha values for Onufriev-Bashford-Case radii */
/*     bobc      beta values for Onufriev-Bashford-Case radii */
/*     gobc      gamma values for Onufriev-Bashford-Case radii */
/*     vsolv     atomic volume of each atom for use with ACE */
/*     wace      "omega" values for atom class pairs for use with ACE */
/*     s2ace     "sigma^2" values for atom class pairs for use with ACE */
/*     uace      "mu" values for atom class pairs for use with ACE */
/*     solvtyp   type of continuum solvation energy model in use */
/*     borntyp   method to be used for the Born radius computation */




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




/*     set up the structure and molecular mechanics calculation */

    initial_();
    getxyz_();
    mechanic_();

/*     initialize the temperature, pressure and coupling baths */

    bath_1.kelvin = 0.;
    bath_1.atmsph = 0.;
    bath_1.isothermal = FALSE_;
    bath_1.isobaric = FALSE_;

/*     initialize the simulation length as number of time steps */

    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___4);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&nstep, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	query = FALSE_;
    }
L10:
    if (query) {
	io___6.ciunit = iounit_1.iout;
	s_wsfe(&io___6);
	e_wsfe();
	io___7.ciunit = iounit_1.input;
	s_rsfe(&io___7);
	do_fio(&c__1, (char *)&nstep, (ftnlen)sizeof(integer));
	e_rsfe();
    }

/*     get the length of the dynamics time step in picoseconds */

    dt = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___9);
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&dt, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L40;
	}
    }
L40:
    while(dt < 0.) {
	io___10.ciunit = iounit_1.iout;
	s_wsfe(&io___10);
	e_wsfe();
	io___11.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___11);
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = do_fio(&c__1, (char *)&dt, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L70;
	}
	if (dt <= 0.) {
	    dt = 1.;
	}
L70:
	;
    }
    dt *= .001;

/*     set bounds on the Berendsen bath coupling parameters */

    bath_1.tautemp = max(bath_1.tautemp,dt);
    bath_1.taupres = max(bath_1.taupres,dt);

/*     set the time between trajectory snapshot coordinate dumps */

    dtdump = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___13);
	if (i__1 != 0) {
	    goto L80;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&dtdump, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L80;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L80;
	}
    }
L80:
    while(dtdump < 0.) {
	io___14.ciunit = iounit_1.iout;
	s_wsfe(&io___14);
	e_wsfe();
	io___15.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___15);
	if (i__1 != 0) {
	    goto L110;
	}
	i__1 = do_fio(&c__1, (char *)&dtdump, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L110;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L110;
	}
	if (dtdump <= 0.) {
	    dtdump = .1;
	}
L110:
	;
    }
    d__1 = dtdump / dt;
    inform_1.iwrite = i_dnnt(&d__1);

/*     get choice of statistical ensemble for periodic system */

    if (bound_1.use_bounds__) {
	mode = -1;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___17);
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&mode, (ftnlen)sizeof(integer)
		    );
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L120;
	    }
	}
L120:
	while(mode < 1 || mode > 4) {
	    io___18.ciunit = iounit_1.iout;
	    s_wsfe(&io___18);
	    e_wsfe();
	    io___19.ciunit = iounit_1.input;
	    i__1 = s_rsfe(&io___19);
	    if (i__1 != 0) {
		goto L150;
	    }
	    i__1 = do_fio(&c__1, (char *)&mode, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L150;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L150;
	    }
	    if (mode <= 0) {
		mode = 1;
	    }
L150:
	    ;
	}
	if (mode == 2 || mode == 4) {
	    bath_1.isothermal = TRUE_;
	    bath_1.kelvin = -1.;
	    nextarg_(string, &exist, (ftnlen)120);
	    if (exist) {
		i__1 = s_rsli(&io___20);
		if (i__1 != 0) {
		    goto L160;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&bath_1.kelvin, (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L160;
		}
		i__1 = e_rsli();
		if (i__1 != 0) {
		    goto L160;
		}
	    }
L160:
	    while(bath_1.kelvin < 0.) {
		io___21.ciunit = iounit_1.iout;
		s_wsfe(&io___21);
		e_wsfe();
		io___22.ciunit = iounit_1.input;
		i__1 = s_rsfe(&io___22);
		if (i__1 != 0) {
		    goto L190;
		}
		i__1 = do_fio(&c__1, (char *)&bath_1.kelvin, (ftnlen)sizeof(
			doublereal));
		if (i__1 != 0) {
		    goto L190;
		}
		i__1 = e_rsfe();
		if (i__1 != 0) {
		    goto L190;
		}
		if (bath_1.kelvin <= 0.) {
		    bath_1.kelvin = 298.;
		}
L190:
		;
	    }
	    bath_1.kelvin0 = bath_1.kelvin;
	}
	if (mode == 3 || mode == 4) {
	    bath_1.isobaric = TRUE_;
	    bath_1.atmsph = -1.;
	    nextarg_(string, &exist, (ftnlen)120);
	    if (exist) {
		i__1 = s_rsli(&io___23);
		if (i__1 != 0) {
		    goto L200;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&bath_1.atmsph, (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L200;
		}
		i__1 = e_rsli();
		if (i__1 != 0) {
		    goto L200;
		}
	    }
L200:
	    while(bath_1.atmsph < 0.) {
		io___24.ciunit = iounit_1.iout;
		s_wsfe(&io___24);
		e_wsfe();
		io___25.ciunit = iounit_1.input;
		i__1 = s_rsfe(&io___25);
		if (i__1 != 0) {
		    goto L230;
		}
		i__1 = do_fio(&c__1, (char *)&bath_1.atmsph, (ftnlen)sizeof(
			doublereal));
		if (i__1 != 0) {
		    goto L230;
		}
		i__1 = e_rsfe();
		if (i__1 != 0) {
		    goto L230;
		}
		if (bath_1.atmsph <= 0.) {
		    bath_1.atmsph = 1.;
		}
L230:
		;
	    }
	}
    }

/*     use constant energy or temperature for nonperiodic system */

    if (! bound_1.use_bounds__) {
	mode = -1;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___26);
	    if (i__1 != 0) {
		goto L240;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&mode, (ftnlen)sizeof(integer)
		    );
	    if (i__1 != 0) {
		goto L240;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L240;
	    }
	}
L240:
	while(mode < 1 || mode > 2) {
	    io___27.ciunit = iounit_1.iout;
	    s_wsfe(&io___27);
	    e_wsfe();
	    io___28.ciunit = iounit_1.input;
	    i__1 = s_rsfe(&io___28);
	    if (i__1 != 0) {
		goto L270;
	    }
	    i__1 = do_fio(&c__1, (char *)&mode, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L270;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L270;
	    }
	    if (mode <= 0) {
		mode = 1;
	    }
L270:
	    ;
	}
	if (mode == 2) {
	    bath_1.isothermal = TRUE_;
	    bath_1.kelvin = -1.;
	    nextarg_(string, &exist, (ftnlen)120);
	    if (exist) {
		i__1 = s_rsli(&io___29);
		if (i__1 != 0) {
		    goto L280;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&bath_1.kelvin, (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L280;
		}
		i__1 = e_rsli();
		if (i__1 != 0) {
		    goto L280;
		}
	    }
L280:
	    while(bath_1.kelvin < 0.) {
		io___30.ciunit = iounit_1.iout;
		s_wsfe(&io___30);
		e_wsfe();
		io___31.ciunit = iounit_1.input;
		i__1 = s_rsfe(&io___31);
		if (i__1 != 0) {
		    goto L310;
		}
		i__1 = do_fio(&c__1, (char *)&bath_1.kelvin, (ftnlen)sizeof(
			doublereal));
		if (i__1 != 0) {
		    goto L310;
		}
		i__1 = e_rsfe();
		if (i__1 != 0) {
		    goto L310;
		}
		if (bath_1.kelvin <= 0.) {
		    bath_1.kelvin = 298.;
		}
L310:
		;
	    }
	}
    }

/*     initialize any rattle constraints and setup dynamics */

    shakeup_();
    mdinit_();

/*     print out a header line for the dynamics computation */

    if (s_cmp(mdstuf_1.integrate, "VERLET", (ftnlen)10, (ftnlen)6) == 0) {
	io___32.ciunit = iounit_1.iout;
	s_wsfe(&io___32);
	e_wsfe();
    } else if (s_cmp(mdstuf_1.integrate, "STOCHASTIC", (ftnlen)10, (ftnlen)10)
	     == 0) {
	io___33.ciunit = iounit_1.iout;
	s_wsfe(&io___33);
	e_wsfe();
    } else if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) 
	    == 0) {
	io___34.ciunit = iounit_1.iout;
	s_wsfe(&io___34);
	e_wsfe();
    } else {
	io___35.ciunit = iounit_1.iout;
	s_wsfe(&io___35);
	e_wsfe();
    }

/*     integrate equations of motion to take a time step */

    i__1 = nstep;
    for (istep = 1; istep <= i__1; ++istep) {
	if (s_cmp(mdstuf_1.integrate, "VERLET", (ftnlen)10, (ftnlen)6) == 0) {
	    verlet_(&istep, &dt);
	} else if (s_cmp(mdstuf_1.integrate, "STOCHASTIC", (ftnlen)10, (
		ftnlen)10) == 0) {
	    sdstep_(&istep, &dt);
	} else if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)
		9) == 0) {
	    rgdstep_(&istep, &dt);
	} else {
	    beeman_(&istep, &dt);
	}

/*     remove center of mass translation and rotation if needed */

	modstep = istep % inform_1.iprint;
	if (modstep == 0 && usage_1.nuse == atoms_1.n) {
	    if (bath_1.isothermal && s_cmp(mdstuf_1.integrate, "STOCHASTIC", (
		    ftnlen)10, (ftnlen)10) != 0 && s_cmp(bath_1.thermostat, 
		    "ANDERSEN", (ftnlen)11, (ftnlen)8) != 0) {
		mdrest_();
	    }
	}
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

/* Main program alias */ int dynamic_ () { MAIN__ (); return 0; }
