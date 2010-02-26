/* anneal.f -- translated by f2c (version 20050501).
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

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__3 = 3;
static doublereal c_b71 = 10.;
static doublereal c_b103 = 3.5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  program anneal  --  molecular dynamics simulated annealing  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "anneal" performs a simulated annealing protocol by means of */
/*     variable temperature molecular dynamics using either linear, */
/*     exponential or sigmoidal cooling schedules */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Enter the Initial and Final Temperatures"
	    " in\002,\002 Degrees K [1000,0] :  \002,$)";
    static char fmt_30[] = "(a120)";
    static char fmt_60[] = "(/,\002 Enter the Number of Equilibration Steps "
	    "[0] :  \002,$)";
    static char fmt_70[] = "(i10)";
    static char fmt_100[] = "(/,\002 Enter the Number of Cooling Protocol St"
	    "eps\002,\002 [2000] :  \002,$)";
    static char fmt_110[] = "(i10)";
    static char fmt_130[] = "(/,\002 Use Linear, Sigmoidal or Exponential Co"
	    "oling\002,\002 Protocol ([L], S or E) :  \002,$)";
    static char fmt_140[] = "(a120)";
    static char fmt_160[] = "(/,\002 Enter the Time Step Length in Femtoseco"
	    "nds\002,\002 [1.0] :  \002,$)";
    static char fmt_170[] = "(f20.0)";
    static char fmt_200[] = "(/,\002 Enter Time between Dumps in Picosecond"
	    "s\002,\002 [0.1] :  \002,$)";
    static char fmt_210[] = "(f20.0)";
    static char fmt_240[] = "(/,\002 Increase Atomic Weights by a Factor o"
	    "f\002,\002 10^x [x=0.0] :  \002,$)";
    static char fmt_250[] = "(f20.0)";
    static char fmt_280[] = "(/,\002 Enter Final Desired Deformation Paramet"
	    "er\002,\002 [0.0] :  \002,$)";
    static char fmt_290[] = "(f20.0)";
    static char fmt_310[] = "(/,\002 Simulated Annealing Equilibration Phas"
	    "e\002)";
    static char fmt_320[] = "(/,\002 Steps:\002,i6,3x,\002Time/Step:\002,f6."
	    "3,\002 ps\002,3x,\002LogMass:\002,f5.2,3x,\002Temp:\002,f7.1,"
	    "\002 to\002,f7.1)";
    static char fmt_330[] = "(/,\002 Simulated Annealing Cooling Protocol"
	    "\002)";
    static char fmt_340[] = "(/,\002 Steps:\002,i6,3x,\002Time/Step:\002,f6."
	    "3,\002 ps\002,3x,\002LogMass:\002,f5.2,3x,\002Temp:\002,f7.1,"
	    "\002 to\002,f7.1)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_dnnt(doublereal *);
    double pow_dd(doublereal *, doublereal *);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double exp(doublereal);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    static integer i__;
    static doublereal dt, hot, cold;
    static integer next;
    extern /* Subroutine */ int final_(void);
    static doublereal sharp, ratio, tight, loose;
    static integer istep, nstep;
    static logical exist;
    static doublereal fuzzy;
    extern /* Subroutine */ int beeman_(integer *, doublereal *);
    static doublereal factor;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static integer nequil;
    static doublereal dtdump;
    static char answer[1], string[120];
    extern /* Subroutine */ int mdinit_(void), verlet_(integer *, doublereal *
	    ), sdstep_(integer *, doublereal *), mdrest_(void), getxyz_(void);
    extern doublereal sigmoid_(doublereal *, doublereal *);
    extern /* Subroutine */ int initial_(void), shakeup_(void);
    static doublereal logmass;
    static integer modstep;
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), rgdstep_(
	    integer *, doublereal *), gettext_(char *, char *, integer *, 
	    ftnlen, ftnlen);
    static char cooltyp[8];

    /* Fortran I/O blocks */
    static icilist io___5 = { 1, string, 1, 0, 120, 1 };
    static icilist io___6 = { 1, string, 1, 0, 120, 1 };
    static cilist io___7 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_30, 0 };
    static icilist io___10 = { 1, record, 1, 0, 120, 1 };
    static icilist io___12 = { 1, string, 1, 0, 120, 1 };
    static cilist io___13 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___14 = { 1, 0, 0, fmt_70, 0 };
    static icilist io___16 = { 1, string, 1, 0, 120, 1 };
    static cilist io___17 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___18 = { 1, 0, 0, fmt_110, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_140, 0 };
    static icilist io___25 = { 1, string, 1, 0, 120, 1 };
    static cilist io___26 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___27 = { 1, 0, 0, fmt_170, 0 };
    static icilist io___29 = { 1, string, 1, 0, 120, 1 };
    static cilist io___30 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___31 = { 1, 0, 0, fmt_210, 0 };
    static icilist io___33 = { 1, string, 1, 0, 120, 1 };
    static cilist io___34 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___35 = { 1, 0, 0, fmt_250, 0 };
    static icilist io___39 = { 1, string, 1, 0, 120, 1 };
    static cilist io___40 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___41 = { 1, 0, 0, fmt_290, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___46 = { 0, 0, 0, fmt_320, 0 };
    static cilist io___49 = { 0, 0, 0, fmt_330, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_340, 0 };




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




/*     set up the structure and mechanics calculation */

    initial_();
    getxyz_();
    mechanic_();

/*     get choice of statistical ensemble for periodic system */

    hot = -1.;
    cold = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___5);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&hot, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    }
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___6);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&cold, (ftnlen)sizeof(doublereal))
		;
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    }
L10:
    while(hot < 0. || cold < 0.) {
	hot = -1.;
	cold = -1.;
	io___7.ciunit = iounit_1.iout;
	s_wsfe(&io___7);
	e_wsfe();
	io___8.ciunit = iounit_1.input;
	s_rsfe(&io___8);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___10);
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&hot, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&cold, (ftnlen)sizeof(doublereal))
		;
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L40;
	}
L40:
	if (hot <= 0.) {
	    hot = 1e3;
	}
	if (cold <= 0.) {
	    cold = 0.;
	}
    }

/*     set the number of steps of initial equilibration */

    nequil = -1;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___12);
	if (i__1 != 0) {
	    goto L50;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&nequil, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L50;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L50;
	}
    }
L50:
    while(nequil < 0) {
	io___13.ciunit = iounit_1.iout;
	s_wsfe(&io___13);
	e_wsfe();
	io___14.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___14);
	if (i__1 != 0) {
	    goto L80;
	}
	i__1 = do_fio(&c__1, (char *)&nequil, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L80;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L80;
	}
	if (nequil <= 0) {
	    nequil = 0;
	}
L80:
	;
    }

/*     set the number of dynamics steps for the cooling protocol */

    nstep = -1;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___16);
	if (i__1 != 0) {
	    goto L90;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&nstep, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L90;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L90;
	}
    }
L90:
    while(nstep < 0) {
	io___17.ciunit = iounit_1.iout;
	s_wsfe(&io___17);
	e_wsfe();
	io___18.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___18);
	if (i__1 != 0) {
	    goto L120;
	}
	i__1 = do_fio(&c__1, (char *)&nstep, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L120;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L120;
	}
	if (nstep <= 0) {
	    nstep = 2000;
	}
L120:
	;
    }

/*     decide which annealing cooling protocol to use */

    s_copy(cooltyp, "LINEAR", (ftnlen)8, (ftnlen)6);
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___21.ciunit = iounit_1.iout;
	s_wsfe(&io___21);
	e_wsfe();
	io___22.ciunit = iounit_1.input;
	s_rsfe(&io___22);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'S') {
	s_copy(cooltyp, "SIGMOID", (ftnlen)8, (ftnlen)7);
    }
    if (*(unsigned char *)answer == 'E') {
	s_copy(cooltyp, "EXPONENT", (ftnlen)8, (ftnlen)8);
    }

/*     get the length of the dynamics time step in picoseconds */

    dt = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___25);
	if (i__1 != 0) {
	    goto L150;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&dt, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L150;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L150;
	}
    }
L150:
    while(dt < 0.) {
	io___26.ciunit = iounit_1.iout;
	s_wsfe(&io___26);
	e_wsfe();
	io___27.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___27);
	if (i__1 != 0) {
	    goto L180;
	}
	i__1 = do_fio(&c__1, (char *)&dt, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L180;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L180;
	}
	if (dt <= 0.) {
	    dt = 1.;
	}
L180:
	;
    }
    dt *= .001;

/*     set the time between trajectory snapshot coordinate dumps */

    dtdump = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___29);
	if (i__1 != 0) {
	    goto L190;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&dtdump, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L190;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L190;
	}
    }
L190:
    while(dtdump < 0.) {
	io___30.ciunit = iounit_1.iout;
	s_wsfe(&io___30);
	e_wsfe();
	io___31.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___31);
	if (i__1 != 0) {
	    goto L220;
	}
	i__1 = do_fio(&c__1, (char *)&dtdump, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L220;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L220;
	}
	if (dtdump <= 0.) {
	    dtdump = .1;
	}
L220:
	;
    }
    d__1 = dtdump / dt;
    inform_1.iwrite = i_dnnt(&d__1);

/*     get factor by which atomic weights are to be increased */

    logmass = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___33);
	if (i__1 != 0) {
	    goto L230;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&logmass, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L230;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L230;
	}
    }
L230:
    while(logmass < 0.) {
	io___34.ciunit = iounit_1.iout;
	s_wsfe(&io___34);
	e_wsfe();
	io___35.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___35);
	if (i__1 != 0) {
	    goto L260;
	}
	i__1 = do_fio(&c__1, (char *)&logmass, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L260;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L260;
	}
	if (logmass <= 0.) {
	    logmass = 0.;
	}
L260:
	;
    }
    if (logmass > 0.) {
	factor = pow_dd(&c_b71, &logmass);
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    atmtyp_1.mass[i__ - 1] *= factor;
	}
    }

/*     rate of deformation change for potential surface smoothing */

    if (warp_1.use_smooth__) {
	sharp = -1.;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___39);
	    if (i__1 != 0) {
		goto L270;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&sharp, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L270;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L270;
	    }
	}
L270:
	while(sharp < 0.) {
	    io___40.ciunit = iounit_1.iout;
	    s_wsfe(&io___40);
	    e_wsfe();
	    io___41.ciunit = iounit_1.input;
	    i__1 = s_rsfe(&io___41);
	    if (i__1 != 0) {
		goto L300;
	    }
	    i__1 = do_fio(&c__1, (char *)&sharp, (ftnlen)sizeof(doublereal));
	    if (i__1 != 0) {
		goto L300;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L300;
	    }
	    if (sharp <= 0.) {
		sharp = 0.;
	    }
L300:
	    ;
	}
	fuzzy = warp_1.deform - sharp;
	if (fuzzy <= 0.) {
	    fuzzy = 0.;
	}
    }

/*     set values for temperature, pressure and coupling baths */

    bath_1.isothermal = TRUE_;
    bath_1.isobaric = FALSE_;
    loose = dt * 100.;
    tight = dt * 10.;
    bath_1.kelvin = hot;
    bath_1.kelvin0 = bath_1.kelvin;
    bath_1.tautemp = loose;

/*     initialize any rattle constraints and setup dynamics */

    shakeup_();
    mdinit_();

/*     print out a header lines for the equilibration phase */

    if (nequil != 0) {
	io___45.ciunit = iounit_1.iout;
	s_wsfe(&io___45);
	e_wsfe();
	io___46.ciunit = iounit_1.iout;
	s_wsfe(&io___46);
	do_fio(&c__1, (char *)&nequil, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&dt, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&logmass, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&hot, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&hot, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     take the dynamics steps for the equilibration phase */

    i__1 = nequil;
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
	modstep = istep % inform_1.iprint;
	if (usage_1.nuse == atoms_1.n && modstep == 0) {
	    mdrest_();
	}
    }

/*     start the cooling phase from the end of equilibration phase */

    if (nequil != 0) {
	mdinit_();
    }

/*     print out a header lines for the cooling protocol */

    io___49.ciunit = iounit_1.iout;
    s_wsfe(&io___49);
    e_wsfe();
    io___50.ciunit = iounit_1.iout;
    s_wsfe(&io___50);
    do_fio(&c__1, (char *)&nstep, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&dt, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&logmass, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&hot, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&cold, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     set target temperature using the desired cooling protocol */

    i__1 = nstep;
    for (istep = 1; istep <= i__1; ++istep) {
	ratio = (doublereal) istep / (doublereal) nstep;
	if (s_cmp(cooltyp, "SIGMOID", (ftnlen)8, (ftnlen)7) == 0) {
	    ratio = sigmoid_(&c_b103, &ratio);
	} else if (s_cmp(cooltyp, "EXPONENT", (ftnlen)8, (ftnlen)8) == 0) {
	    ratio = 1. - exp(ratio * -5.);
	}
	bath_1.kelvin = hot * (1. - ratio) + cold * ratio;
	bath_1.kelvin0 = bath_1.kelvin;
	bath_1.tautemp = loose * (1. - ratio) + tight * ratio;

/*     set the deformation value if potential smoothing is used */

	if (warp_1.use_smooth__) {
/* Computing 3rd power */
	    d__1 = 1. - (doublereal) istep / (doublereal) nstep;
	    ratio = d__1 * (d__1 * d__1);
	    warp_1.deform = sharp + ratio * fuzzy;
	}

/*     integrate equations of motion to take a time step */

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
	    if (s_cmp(mdstuf_1.integrate, "STOCHASTIC", (ftnlen)10, (ftnlen)
		    10) != 0 && s_cmp(bath_1.thermostat, "ANDERSEN", (ftnlen)
		    11, (ftnlen)8) != 0) {
		mdrest_();
	    }
	}
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

/* Main program alias */ int anneal_ () { MAIN__ (); return 0; }
