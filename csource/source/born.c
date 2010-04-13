/* born.f -- translated by f2c (version 20050501).
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
    doublereal electric, dielec, ebuffer, c2scale, c3scale, c4scale, c5scale;
    logical neutnbr, neutcut;
} chgpot_;

#define chgpot_1 chgpot_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

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
    doublereal pbe, apbe[25000], pbr[25000], pbep[75000]	/* was [3][
	    25000] */, pbfp[75000]	/* was [3][25000] */, pbtp[75000]	
	    /* was [3][25000] */, pbeuind[75000]	/* was [3][25000] */, 
	    pbeuinp[75000]	/* was [3][25000] */, grid[3], gcent[3], 
	    cgrid[3], cgcent[3], fgrid[3], fgcent[3], ionr[10], ionc[10], 
	    pdie, sdie, srad, swin, sdens, smin;
    integer ionn, dime[3], ionq[10];
    char pbtyp[20], pbsoln[20], bcfl[20], srfm[20], chgm[20];
} pb_;

#define pb_1 pb_

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
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

struct {
    doublereal vir[9]	/* was [3][3] */;
} virial_;

#define virial_1 virial_

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b21 = 3.141592653589793238;
static doublereal c_b22 = 6.;
static doublereal c_b23 = 3.;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine born  --  Born radii for continuum solvation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "born" computes the Born radius of each atom for use with */
/*     the various continuum solvation models */

/*     literature references: */

/*     W. C. Still, A. Tempczyk, R. C. Hawley and T. Hendrickson, */
/*     "A Semianalytical Treatment of Solvation for Molecular */
/*     Mechanics and Dynamics", J. Amer. Chem. Soc., 112, 6127-6129 */
/*     (1990)  ("Onion" Method; see supplimentary material) */

/*     D. Qiu, P. S. Shenkin, F. P. Hollinger and W. C. Still, "The */
/*     GB/SA Continuum Model for Solvation. A Fast Analytical Method */
/*     for the Calculation of Approximate Radii", J. Phys. Chem. A, */
/*     101, 3005-3014 (1997)  (Analytical Still Method) */

/*     G. D. Hawkins, C. J. Cramer and D. G. Truhlar, "Parametrized */
/*     Models of Aqueous Free Energies of Solvation Based on Pairwise */
/*     Descreening of Solute Atomic Charges from a Dielectric Medium", */
/*     J. Phys. Chem., 100, 19824-19839 (1996)  (HCT Method) */

/*     A. Onufriev, D. Bashford and D. A. Case, "Exploring Protein */
/*     Native States and Large-Scale Conformational Changes with a */
/*     Modified Generalized Born Model", PROTEINS, 55, 383-394 (2004) */
/*     (OBC Method) */

/*     T. Grycuk, "Deficiency of the Coulomb-field Approximation */
/*     in the Generalized Born Model: An Improved Formula for Born */
/*     Radii Evaluation", J. Chem. Phys., 119, 4817-4826 (2003) */
/*     (Grycuk Method) */

/*     M. Schaefer, C. Bartels and M. Karplus, "Solution Conformations */
/*     and Thermodynamics of Structured Peptides: Molecular Dynamics */
/*     Simulation with an Implicit Solvation Model", J. Mol. Biol., */
/*     284, 835-848 (1998)  (ACE Method) */


/* Subroutine */ int born_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Born Radii for Individual Atoms :\002,/)";
    static char fmt_20[] = "(1x,5(i7,f8.3))";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double cos(doublereal), sqrt(doublereal), log(doublereal), tanh(
	    doublereal), pow_dd(doublereal *, doublereal *), exp(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal fraction;
    extern /* Subroutine */ int surfatom_(integer *, doublereal *, doublereal 
	    *);
    static integer i__, j, k;
    static doublereal r__, t, b0, l2, l4, r2, r3, r4, u2, u4, ri;
    static integer it;
    static doublereal rk;
    static integer kt;
    static doublereal sk, xi, yi, zi, lr, ur, xr, yr, zr;
    extern /* Subroutine */ int apbsempole_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal sk2, l4r, u4r, ccf, lik, gpi, pi43, uik, pos[75000]	
	    /* was [3][25000] */, rmu, sum, lik2, uik2, pip5, sum2, sum3, 
	    area, beta;
    static logical done;
    static doublereal roff[25000];
    static integer skip[25000];
    static doublereal rold, term, rvdw, tsum, p5inv, gamma, alpha, gself, 
	    theta, shell, inner, ratio, third, total, tinit, outer, tchain, 
	    pbpole[325000]	/* was [13][25000] */, bornmax, expterm;

    /* Fortran I/O blocks */
    static cilist io___69 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_20, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define pos_ref(a_1,a_2) pos[(a_2)*3 + a_1 - 4]
#define uace_ref(a_1,a_2) solute_1.uace[(a_2)*1000 + a_1 - 1001]
#define wace_ref(a_1,a_2) solute_1.wace[(a_2)*1000 + a_1 - 1001]
#define s2ace_ref(a_1,a_2) solute_1.s2ace[(a_2)*1000 + a_1 - 1001]
#define pbpole_ref(a_1,a_2) pbpole[(a_2)*13 + a_1 - 14]



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
/*     ##  chgpot.i  --  specifics of charge-charge functional form  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     electric   energy factor in kcal/mole for current force field */
/*     dielec     dielectric constant for electrostatic interactions */
/*     ebuffer    electrostatic buffering constant added to distance */
/*     c2scale    factor by which 1-2 charge interactions are scaled */
/*     c3scale    factor by which 1-3 charge interactions are scaled */
/*     c4scale    factor by which 1-4 charge interactions are scaled */
/*     c5scale    factor by which 1-5 charge interactions are scaled */
/*     neutnbr    logical flag governing use of neutral group neighbors */
/*     neutcut    logical flag governing use of neutral group cutoffs */




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




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  pb.i  --  parameters for Poisson-Boltzmann solvation  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     pbe      Poisson-Boltzman permanent multipole solvation energy */
/*     apbe     Poisson-Boltzman permanent multipole energy over atoms */
/*     pbr      Poisson-Boltzman cavity radii for atom types */
/*     pbep     Poisson-Boltzman energies on permanent multipoles */
/*     pbfp     Poisson-Boltzman forces on permanent multipoles */
/*     pbtp     Poisson-Boltzman torques on permanent multipoles */
/*     pbeuind  Poisson-Boltzman field due to induced dipoles */
/*     pbeuinp  Poisson-Boltzman field due to non-local induced dipoles */

/*     APBS configuration parameters (see APBS documentation for more details) */
/*     In the column on the right are possible values for each variable, with */
/*     the default values given in brackets. Note that only a subset of APBS */
/*     options are supported and/or are appropriate for use with AMOEBA. */

/*     pbtyp                                   lpbe */

/*     At some point AMOEBA with the non-linear PBE could be supported, but */
/*     there is only have theory for energies (no gradients). */

/*     pbsoln                                  mg-auto, [mg-manual] */

/*     Currently there is only limited support for focusing calculations, */
/*     which is a powerful feature of APBS. The current requirement is */
/*     that energies and forces must all be calculated using the finest */
/*     solution. */

/*     bcfl     boundary conditions            zero, sdh, [mdh] */
/*     chgm     multipole discretization       spl4 */

/*     other charge discretization methods are not appropriate for AMOEBA */

/*     srfm     surface method                 mol, smol, [spl4] */

/*     spl4 is required for forces calculations, although mol is useful for */
/*     comparison with generalized Kirkwood */

/*     dime     number of grid points          [65, 65, 65] */
/*     grid     grid spacing (mg-manual)       fxn of "dime" */
/*     cgrid    coarse grid spacing            fxn of "dime" */
/*     fgrid    fine grid spacing              cgrid / 2 */

/*     stable results require grid spacing to be fine enough to keep */
/*     multipoles inside the dielectric boundary (2.5 * grid < PBR) */

/*     gcent    grid center  (mg-manual)       center of mass */
/*     cgcent   coarse grid center             center of mass */
/*     fgcent   fine grid center               center of mass */
/*     pdie     solute/homogeneous dieletric   [1.0] */
/*     sdie     solvent dieletric              [78.3] */
/*     ionn     number of ion species          [0] */
/*     ionc     ion concentration (M)          [0.0] */
/*     ionq     ion charge (electrons)         [1.0] */
/*     ionr     ion radius (A)                 [2.0] */
/*     srad     solvent probe radius (A)       [1.4] */
/*     swin     surface spline window width    [0.3] */
/*     sdens    density of surface points      [10.0] */

/*     additional parameter to facilitate default grid setup */

/*     smin     minimum distance between an    [10.0] */
/*              atomic center and the grid */
/*              boundary (A) */




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




/*     set offset modified radii and OBC chain rule factor */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	roff[i__ - 1] = solute_1.rsolv[i__ - 1] - solute_1.doffset;
	solute_1.drobc[i__ - 1] = 1.;
    }

/*     get the Born radii via the numerical "Onion" method */

    if (s_cmp(solute_1.borntyp, "ONION", (ftnlen)8, (ftnlen)5) == 0) {
	tinit = .1;
	ratio = 1.5;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    t = tinit;
	    rold = roff[i__ - 1];
	    total = 0.;
	    done = FALSE_;
	    while(! done) {
		roff[i__ - 1] += t * .5;
		surfatom_(&i__, &area, roff);
/* Computing 2nd power */
		d__1 = roff[i__ - 1];
		fraction = area / (d__1 * d__1 * 12.566370614359172);
		if (fraction < .99) {
		    inner = roff[i__ - 1] - t * .5;
		    outer = inner + t;
		    shell = 1. / inner - 1. / outer;
		    total += fraction * shell;
		    roff[i__ - 1] += t * .5;
		    t = ratio * t;
		} else {
		    inner = roff[i__ - 1] - t * .5;
		    total += 1. / inner;
		    done = TRUE_;
		}
	    }
	    solute_1.rborn[i__ - 1] = 1. / total;
	    roff[i__ - 1] = rold;
	}

/*     get the Born radii via the analytical Still method; */
/*     note this code only loops over the variable parts */

    } else if (s_cmp(solute_1.borntyp, "STILL", (ftnlen)8, (ftnlen)5) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    skip[i__ - 1] = 0;
	}
	p5inv = 1. / solute_1.p5;
	pip5 = solute_1.p5 * 3.141592653589793238;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    gpi = solute_1.gpol[i__ - 1];
	    skip[i__ - 1] = i__;
	    i__2 = couple_1.n12[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		skip[i12_ref(j, i__) - 1] = i__;
	    }
	    i__2 = couple_1.n13[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		skip[i13_ref(j, i__) - 1] = i__;
	    }
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		if (skip[k - 1] != i__) {
		    xr = atoms_1.x[k - 1] - xi;
		    yr = atoms_1.y[k - 1] - yi;
		    zr = atoms_1.z__[k - 1] - zi;
/* Computing 2nd power */
		    d__1 = xr;
/* Computing 2nd power */
		    d__2 = yr;
/* Computing 2nd power */
		    d__3 = zr;
		    r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
		    r4 = r2 * r2;
		    rvdw = solute_1.rsolv[i__ - 1] + solute_1.rsolv[k - 1];
		    ratio = r2 / (rvdw * rvdw);
		    if (ratio > p5inv) {
			ccf = 1.;
		    } else {
			theta = ratio * pip5;
			term = (1. - cos(theta)) * .5;
			ccf = term * term;
		    }
		    gpi += solute_1.p4 * ccf * solute_1.vsolv[k - 1] / r4;
		}
	    }
	    solute_1.rborn[i__ - 1] = chgpot_1.electric * -.5 / gpi;
	}

/*     get the Born radii via the Hawkins-Cramer-Truhlar method */

    } else if (s_cmp(solute_1.borntyp, "HCT", (ftnlen)8, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    ri = roff[i__ - 1];
	    sum = 1. / ri;
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		if (i__ != k) {
		    xr = atoms_1.x[k - 1] - xi;
		    yr = atoms_1.y[k - 1] - yi;
		    zr = atoms_1.z__[k - 1] - zi;
/* Computing 2nd power */
		    d__1 = xr;
/* Computing 2nd power */
		    d__2 = yr;
/* Computing 2nd power */
		    d__3 = zr;
		    r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
		    r__ = sqrt(r2);
		    rk = roff[k - 1];
		    sk = rk * solute_1.shct[k - 1];
		    sk2 = sk * sk;
		    if (ri < r__ + sk) {
/* Computing MAX */
			d__2 = ri, d__3 = (d__1 = r__ - sk, abs(d__1));
			lik = 1. / max(d__2,d__3);
			uik = 1. / (r__ + sk);
			lik2 = lik * lik;
			uik2 = uik * uik;
			term = lik - uik + r__ * .25 * (uik2 - lik2) + .5 / 
				r__ * log(uik / lik) + sk2 * .25 / r__ * (
				lik2 - uik2);
			if (ri < sk - r__) {
			    term += (1. / ri - lik) * 2.;
			}
			sum -= term * .5;
		    }
		}
	    }
	    solute_1.rborn[i__ - 1] = 1. / sum;
	}

/*     get the Born radii via the Onufriev-Bashford-Case method */

    } else if (s_cmp(solute_1.borntyp, "OBC", (ftnlen)8, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    ri = roff[i__ - 1];
	    sum = 0.;
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		if (i__ != k) {
		    xr = atoms_1.x[k - 1] - xi;
		    yr = atoms_1.y[k - 1] - yi;
		    zr = atoms_1.z__[k - 1] - zi;
/* Computing 2nd power */
		    d__1 = xr;
/* Computing 2nd power */
		    d__2 = yr;
/* Computing 2nd power */
		    d__3 = zr;
		    r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
		    r__ = sqrt(r2);
		    rk = roff[k - 1];
		    sk = rk * solute_1.shct[k - 1];
		    sk2 = sk * sk;
		    if (ri < r__ + sk) {
/* Computing MAX */
			d__2 = ri, d__3 = (d__1 = r__ - sk, abs(d__1));
			lik = 1. / max(d__2,d__3);
			uik = 1. / (r__ + sk);
			lik2 = lik * lik;
			uik2 = uik * uik;
			term = lik - uik + r__ * .25 * (uik2 - lik2) + .5 / 
				r__ * log(uik / lik) + sk2 * .25 / r__ * (
				lik2 - uik2);
			if (ri < sk - r__) {
			    term += (1. / ri - lik) * 2.;
			}
			sum += term * .5;
		    }
		}
	    }
	    alpha = solute_1.aobc[i__ - 1];
	    beta = solute_1.bobc[i__ - 1];
	    gamma = solute_1.gobc[i__ - 1];
	    sum = ri * sum;
	    sum2 = sum * sum;
	    sum3 = sum * sum2;
	    tsum = tanh(alpha * sum - beta * sum2 + gamma * sum3);
	    solute_1.rborn[i__ - 1] = 1. / ri - tsum / solute_1.rsolv[i__ - 1]
		    ;
	    solute_1.rborn[i__ - 1] = 1. / solute_1.rborn[i__ - 1];
	    tchain = ri * (alpha - beta * 2. * sum + gamma * 3. * sum2);
	    solute_1.drobc[i__ - 1] = (1. - tsum * tsum) * tchain / 
		    solute_1.rsolv[i__ - 1];
	}

/*     get the Born radii via Grycuk's modified HCT method */

    } else if (s_cmp(solute_1.borntyp, "GRYCUK", (ftnlen)8, (ftnlen)6) == 0) {
	third = .33333333333333331;
	pi43 = third * 4. * 3.141592653589793238;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    solute_1.rborn[i__ - 1] = 0.;
	    ri = solute_1.rsolv[i__ - 1];
	    if (ri > 0.) {
		xi = atoms_1.x[i__ - 1];
		yi = atoms_1.y[i__ - 1];
		zi = atoms_1.z__[i__ - 1];
/* Computing 3rd power */
		d__1 = ri;
		sum = pi43 / (d__1 * (d__1 * d__1));
		i__2 = atoms_1.n;
		for (k = 1; k <= i__2; ++k) {
		    rk = solute_1.rsolv[k - 1];
		    if (i__ != k && rk > 0.) {
			xr = atoms_1.x[k - 1] - xi;
			yr = atoms_1.y[k - 1] - yi;
			zr = atoms_1.z__[k - 1] - zi;
/* Computing 2nd power */
			d__1 = xr;
/* Computing 2nd power */
			d__2 = yr;
/* Computing 2nd power */
			d__3 = zr;
			r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
			r__ = sqrt(r2);
			sk = rk * solute_1.shct[k - 1];
			sk2 = sk * sk;
			if (ri + r__ < sk) {
			    lik = ri;
			    uik = sk - r__;
/* Computing 3rd power */
			    d__1 = uik;
/* Computing 3rd power */
			    d__2 = lik;
			    sum += pi43 * (1. / (d__1 * (d__1 * d__1)) - 1. / 
				    (d__2 * (d__2 * d__2)));
			}
			uik = r__ + sk;
			if (ri + r__ < sk) {
			    lik = sk - r__;
			} else if (r__ < ri + sk) {
			    lik = ri;
			} else {
			    lik = r__ - sk;
			}
			l2 = lik * lik;
			l4 = l2 * l2;
			lr = lik * r__;
			l4r = l4 * r__;
			u2 = uik * uik;
			u4 = u2 * u2;
			ur = uik * r__;
			u4r = u4 * r__;
			term = ((r2 - sk2) * 3. + u2 * 6. - ur * 8.) / u4r - (
				(r2 - sk2) * 3. + l2 * 6. - lr * 8.) / l4r;
			sum -= term * 3.141592653589793238 / 12.;
		    }
		}
		d__1 = sum / pi43;
		solute_1.rborn[i__ - 1] = pow_dd(&d__1, &third);
		if (solute_1.rborn[i__ - 1] <= 0.) {
		    solute_1.rborn[i__ - 1] = 1e-4;
		}
		solute_1.rborn[i__ - 1] = 1. / solute_1.rborn[i__ - 1];
	    }
	}

/*     get the Born radii via analytical continuum electrostatics */

    } else if (s_cmp(solute_1.borntyp, "ACE", (ftnlen)8, (ftnlen)3) == 0) {
	third = .33333333333333331;
	b0 = 0.;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    b0 += solute_1.vsolv[i__ - 1];
	}
	d__1 = b0 * .75 / 3.141592653589793238;
	b0 = pow_dd(&d__1, &third);
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    ri = solute_1.rsolv[i__ - 1];
	    it = atmtyp_1.class__[i__ - 1];
	    gself = 1. / ri + wace_ref(it, it) * 2.;
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		if (k != i__) {
		    xr = atoms_1.x[k - 1] - xi;
		    yr = atoms_1.y[k - 1] - yi;
		    zr = atoms_1.z__[k - 1] - zi;
		    kt = atmtyp_1.class__[k - 1];
/* Computing 2nd power */
		    d__1 = xr;
/* Computing 2nd power */
		    d__2 = yr;
/* Computing 2nd power */
		    d__3 = zr;
		    r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
		    r3 = r2 * sqrt(r2);
		    r4 = r2 * r2;
		    expterm = wace_ref(it, kt) * exp(-r2 / s2ace_ref(it, kt));
/* Computing 4th power */
		    d__1 = uace_ref(it, kt), d__1 *= d__1;
		    rmu = r4 + d__1 * d__1;
/* Computing 4th power */
		    d__1 = r3 / rmu, d__1 *= d__1;
		    term = solute_1.vsolv[k - 1] / 25.132741228718345 * (d__1 
			    * d__1);
		    gself -= (expterm + term) * 2.;
		}
	    }
	    if (gself >= .5 / b0) {
		solute_1.rborn[i__ - 1] = 1. / gself;
	    } else {
		solute_1.rborn[i__ - 1] = b0 * 2. * (b0 * gself + 1.);
	    }
	}

/*     get the "perfect" Born radii via Poisson-Boltzmann */

    } else if (s_cmp(solute_1.borntyp, "PERFECT", (ftnlen)8, (ftnlen)7) == 0) 
	    {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    pos_ref(1, i__) = atoms_1.x[i__ - 1];
	    pos_ref(2, i__) = atoms_1.y[i__ - 1];
	    pos_ref(3, i__) = atoms_1.z__[i__ - 1];
	    for (j = 1; j <= 13; ++j) {
		pbpole_ref(j, i__) = 0.;
	    }
	}
	term = chgpot_1.electric * -.5 * (1. - 1. / pb_1.sdie);
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    pbpole_ref(1, i__) = 1.;
	    apbsempole_(&atoms_1.n, pos, solute_1.rsolv, pbpole, &pb_1.pbe, 
		    pb_1.apbe, pb_1.pbep, pb_1.pbfp, pb_1.pbtp);
	    pbpole_ref(1, i__) = 0.;
	    solute_1.rborn[i__ - 1] = term / pb_1.pbe;
	}
    }

/*     make sure the final values are in a reasonable range */

    bornmax = 500.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (solute_1.rborn[i__ - 1] < 0. || solute_1.rborn[i__ - 1] > bornmax)
		 {
	    solute_1.rborn[i__ - 1] = bornmax;
	}
    }

/*     write out the final Born radius value for each atom */

    if (inform_1.debug) {
	io___69.ciunit = iounit_1.iout;
	s_wsfe(&io___69);
	e_wsfe();
	k = 1;
	while(k <= atoms_1.n) {
	    io___70.ciunit = iounit_1.iout;
	    s_wsfe(&io___70);
/* Computing MIN */
	    i__2 = k + 4;
	    i__1 = min(i__2,atoms_1.n);
	    for (i__ = k; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&solute_1.rborn[i__ - 1], (ftnlen)
			sizeof(doublereal));
	    }
	    e_wsfe();
	    k += 5;
	}
    }
    return 0;
} /* born_ */

#undef pbpole_ref
#undef s2ace_ref
#undef wace_ref
#undef uace_ref
#undef pos_ref
#undef i13_ref
#undef i12_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine born1  --  Born radii chain rule derivatives  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "born1" computes derivatives of the Born radii with respect */
/*     to atomic coordinates and increments total energy derivatives */
/*     and virial components for potentials involving Born radii */


/* Subroutine */ int born1_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal), cos(doublereal), sin(doublereal), log(doublereal)
	    , pow_dd(doublereal *, doublereal *), exp(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal r__, r2, r3, r4, r6, t1, t2, t3, de, ri;
    static integer it;
    static doublereal rk;
    static integer kt;
    static doublereal sk, xi, yi, zi, vk, vi, xr, yr, zr, de1, de2, rb2, sk2, 
	    ws2, ccf, gpi, lik, rbi, pi43, dbr, uik, vxx, vyx, vyy, vzx, vzz, 
	    vzy, rbi2, lik2, lik3, s2ik, uik2, uik3, uik4, pip5, dccf, dlik, 
	    dedx, dedy, dedz, duik, roff[25000], cosq;
    static integer skip[25000];
    static doublereal term, sinq, p5inv, theta, third, ratio, rusum, factor, 
	    expterm;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define des_ref(a_1,a_2) deriv_1.des[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define uace_ref(a_1,a_2) solute_1.uace[(a_2)*1000 + a_1 - 1001]
#define wace_ref(a_1,a_2) solute_1.wace[(a_2)*1000 + a_1 - 1001]
#define s2ace_ref(a_1,a_2) solute_1.s2ace[(a_2)*1000 + a_1 - 1001]



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
/*     ##  chgpot.i  --  specifics of charge-charge functional form  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     electric   energy factor in kcal/mole for current force field */
/*     dielec     dielectric constant for electrostatic interactions */
/*     ebuffer    electrostatic buffering constant added to distance */
/*     c2scale    factor by which 1-2 charge interactions are scaled */
/*     c3scale    factor by which 1-3 charge interactions are scaled */
/*     c4scale    factor by which 1-4 charge interactions are scaled */
/*     c5scale    factor by which 1-5 charge interactions are scaled */
/*     neutnbr    logical flag governing use of neutral group neighbors */
/*     neutcut    logical flag governing use of neutral group cutoffs */




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

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  virial.i  --  components of internal virial tensor  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     vir    total internal virial Cartesian tensor components */




/*     compute atomic radii modified by the dielectric offset */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	roff[i__ - 1] = solute_1.rsolv[i__ - 1] - solute_1.doffset;
    }

/*     get Born radius chain rule components for the Still method */

    if (s_cmp(solute_1.borntyp, "STILL", (ftnlen)8, (ftnlen)5) == 0) {
	p5inv = 1. / solute_1.p5;
	pip5 = solute_1.p5 * 3.141592653589793238;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    skip[i__ - 1] = 0;
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    skip[i__ - 1] = i__;
	    i__2 = couple_1.n12[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		skip[i12_ref(j, i__) - 1] = i__;
	    }
	    i__2 = couple_1.n13[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		skip[i13_ref(j, i__) - 1] = i__;
	    }
/* Computing 2nd power */
	    d__1 = solute_1.rborn[i__ - 1];
	    gpi = d__1 * d__1 * 2. / chgpot_1.electric;
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		if (skip[k - 1] != i__) {
		    xr = atoms_1.x[k - 1] - xi;
		    yr = atoms_1.y[k - 1] - yi;
		    zr = atoms_1.z__[k - 1] - zi;
		    vk = solute_1.vsolv[k - 1];
/* Computing 2nd power */
		    d__1 = xr;
/* Computing 2nd power */
		    d__2 = yr;
/* Computing 2nd power */
		    d__3 = zr;
		    r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
		    r__ = sqrt(r2);
		    r6 = r2 * r2 * r2;
/* Computing 2nd power */
		    d__1 = solute_1.rsolv[i__ - 1] + solute_1.rsolv[k - 1];
		    ratio = r2 / (d__1 * d__1);
		    if (ratio > p5inv) {
			ccf = 1.;
			dccf = 0.;
		    } else {
			theta = ratio * pip5;
			cosq = cos(theta);
			term = (1. - cosq) * .5;
			ccf = term * term;
			sinq = sin(theta);
			dccf = term * 2. * sinq * pip5 * ratio;
		    }
		    de = solute_1.drb[i__ - 1] * solute_1.p4 * gpi * vk * (
			    ccf * 4. - dccf) / r6;

/*     increment the overall continuum solvation derivatives */

		    dedx = de * xr;
		    dedy = de * yr;
		    dedz = de * zr;
		    des_ref(1, i__) = des_ref(1, i__) + dedx;
		    des_ref(2, i__) = des_ref(2, i__) + dedy;
		    des_ref(3, i__) = des_ref(3, i__) + dedz;
		    des_ref(1, k) = des_ref(1, k) - dedx;
		    des_ref(2, k) = des_ref(2, k) - dedy;
		    des_ref(3, k) = des_ref(3, k) - dedz;

/*     increment the internal virial tensor components */

		    vxx = xr * dedx;
		    vyx = yr * dedx;
		    vzx = zr * dedx;
		    vyy = yr * dedy;
		    vzy = zr * dedy;
		    vzz = zr * dedz;
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
	}

/*     get Born radius chain rule components for the HCT method */

    } else if (s_cmp(solute_1.borntyp, "HCT", (ftnlen)8, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    ri = roff[i__ - 1];
	    rb2 = solute_1.rborn[i__ - 1] * solute_1.rborn[i__ - 1];
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		if (k != i__) {
		    xr = atoms_1.x[k - 1] - xi;
		    yr = atoms_1.y[k - 1] - yi;
		    zr = atoms_1.z__[k - 1] - zi;
		    rk = roff[k - 1];
		    sk = rk * solute_1.shct[k - 1];
		    sk2 = sk * sk;
/* Computing 2nd power */
		    d__1 = xr;
/* Computing 2nd power */
		    d__2 = yr;
/* Computing 2nd power */
		    d__3 = zr;
		    r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
		    r__ = sqrt(r2);
		    if (ri < r__ + sk) {
/* Computing MAX */
			d__2 = ri, d__3 = (d__1 = r__ - sk, abs(d__1));
			lik = 1. / max(d__2,d__3);
			uik = 1. / (r__ + sk);
			lik2 = lik * lik;
			uik2 = uik * uik;
			lik3 = lik * lik2;
			uik3 = uik * uik2;
			dlik = 1.;
			if (ri >= r__ - sk) {
			    dlik = 0.;
			}
			duik = 1.;
			t1 = lik2 * .5 + sk2 * .25 * lik3 / r__ - (lik / r__ 
				+ lik3 * r__) * .25;
			t2 = uik2 * -.5 - sk2 * .25 * uik3 / r__ + (uik / r__ 
				+ uik3 * r__) * .25;
			t3 = (sk2 / r2 + 1.) * .125 * (lik2 - uik2) + log(uik 
				/ lik) * .25 / r2;
			de = solute_1.drb[i__ - 1] * rb2 * (dlik * t1 + duik *
				 t2 + t3) / r__;

/*     increment the overall continuum solvation derivatives */

			dedx = de * xr;
			dedy = de * yr;
			dedz = de * zr;
			des_ref(1, i__) = des_ref(1, i__) + dedx;
			des_ref(2, i__) = des_ref(2, i__) + dedy;
			des_ref(3, i__) = des_ref(3, i__) + dedz;
			des_ref(1, k) = des_ref(1, k) - dedx;
			des_ref(2, k) = des_ref(2, k) - dedy;
			des_ref(3, k) = des_ref(3, k) - dedz;

/*     increment the internal virial tensor components */

			vxx = xr * dedx;
			vyx = yr * dedx;
			vzx = zr * dedx;
			vyy = yr * dedy;
			vzy = zr * dedy;
			vzz = zr * dedz;
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
	    }
	}

/*     get Born radius chain rule components for the OBC method */

    } else if (s_cmp(solute_1.borntyp, "OBC", (ftnlen)8, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    ri = roff[i__ - 1];
	    rb2 = solute_1.rborn[i__ - 1] * solute_1.rborn[i__ - 1] * 
		    solute_1.drobc[i__ - 1];
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		if (k != i__) {
		    xr = atoms_1.x[k - 1] - xi;
		    yr = atoms_1.y[k - 1] - yi;
		    zr = atoms_1.z__[k - 1] - zi;
		    rk = roff[k - 1];
		    sk = rk * solute_1.shct[k - 1];
		    sk2 = sk * sk;
/* Computing 2nd power */
		    d__1 = xr;
/* Computing 2nd power */
		    d__2 = yr;
/* Computing 2nd power */
		    d__3 = zr;
		    r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
		    r__ = sqrt(r2);
		    if (ri < r__ + sk) {
/* Computing MAX */
			d__2 = ri, d__3 = (d__1 = r__ - sk, abs(d__1));
			lik = 1. / max(d__2,d__3);
			uik = 1. / (r__ + sk);
			lik2 = lik * lik;
			uik2 = uik * uik;
			lik3 = lik * lik2;
			uik3 = uik * uik2;
			dlik = 1.;
			if (ri >= r__ - sk) {
			    dlik = 0.;
			}
			duik = 1.;
			t1 = lik2 * .5 + sk2 * .25 * lik3 / r__ - (lik / r__ 
				+ lik3 * r__) * .25;
			t2 = uik2 * -.5 - sk2 * .25 * uik3 / r__ + (uik / r__ 
				+ uik3 * r__) * .25;
			t3 = (sk2 / r2 + 1.) * .125 * (lik2 - uik2) + log(uik 
				/ lik) * .25 / r2;
			de = solute_1.drb[i__ - 1] * rb2 * (dlik * t1 + duik *
				 t2 + t3) / r__;

/*     increment the overall permanent solvation derivatives */

			dedx = de * xr;
			dedy = de * yr;
			dedz = de * zr;
			des_ref(1, i__) = des_ref(1, i__) + dedx;
			des_ref(2, i__) = des_ref(2, i__) + dedy;
			des_ref(3, i__) = des_ref(3, i__) + dedz;
			des_ref(1, k) = des_ref(1, k) - dedx;
			des_ref(2, k) = des_ref(2, k) - dedy;
			des_ref(3, k) = des_ref(3, k) - dedz;

/*     increment the internal virial tensor components */

			vxx = xr * dedx;
			vyx = yr * dedx;
			vzx = zr * dedx;
			vyy = yr * dedy;
			vzy = zr * dedy;
			vzz = zr * dedz;
			vir_ref(1, 1) = vir_ref(1, 1) + vxx;
			vir_ref(2, 1) = vir_ref(2, 1) + vyx;
			vir_ref(3, 1) = vir_ref(3, 1) + vzx;
			vir_ref(1, 2) = vir_ref(1, 2) + vyx;
			vir_ref(2, 2) = vir_ref(2, 2) + vyy;
			vir_ref(3, 2) = vir_ref(3, 2) + vzy;
			vir_ref(1, 3) = vir_ref(1, 3) + vzx;
			vir_ref(2, 3) = vir_ref(2, 3) + vzy;
			vir_ref(3, 3) = vir_ref(3, 3) + vzz;

/*     increment the polarization solvation derivatives */

			if (potent_1.use_mpole__ || potent_1.use_polar__) {
			    de = solute_1.drbp[i__ - 1] * rb2 * (dlik * t1 + 
				    duik * t2 + t3) / r__;
			    dedx = de * xr;
			    dedy = de * yr;
			    dedz = de * zr;
			    des_ref(1, i__) = des_ref(1, i__) + dedx;
			    des_ref(2, i__) = des_ref(2, i__) + dedy;
			    des_ref(3, i__) = des_ref(3, i__) + dedz;
			    des_ref(1, k) = des_ref(1, k) - dedx;
			    des_ref(2, k) = des_ref(2, k) - dedy;
			    des_ref(3, k) = des_ref(3, k) - dedz;

/*     increment the internal virial tensor components */

			    vxx = xr * dedx;
			    vyx = yr * dedx;
			    vzx = zr * dedx;
			    vyy = yr * dedy;
			    vzy = zr * dedy;
			    vzz = zr * dedz;
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
		}
	    }
	}

/*     get Born radius chain rule components for Grycuk's HCT method */

    } else if (s_cmp(solute_1.borntyp, "GRYCUK", (ftnlen)8, (ftnlen)6) == 0) {
	third = .33333333333333331;
	pi43 = third * 4. * 3.141592653589793238;
	d__1 = third * 2.;
	factor = -pow_dd(&c_b21, &third) * pow_dd(&c_b22, &d__1) / 9.;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ri = solute_1.rsolv[i__ - 1];
	    if (ri > 0.) {
		xi = atoms_1.x[i__ - 1];
		yi = atoms_1.y[i__ - 1];
		zi = atoms_1.z__[i__ - 1];
		term = pi43 / pow_dd(&solute_1.rborn[i__ - 1], &c_b23);
		d__1 = third * 4.;
		term = factor / pow_dd(&term, &d__1);
		i__2 = atoms_1.n;
		for (k = 1; k <= i__2; ++k) {
		    rk = solute_1.rsolv[k - 1];
		    if (k != i__ && rk > 0.) {
			xr = atoms_1.x[k - 1] - xi;
			yr = atoms_1.y[k - 1] - yi;
			zr = atoms_1.z__[k - 1] - zi;
			sk = rk * solute_1.shct[k - 1];
			sk2 = sk * sk;
/* Computing 2nd power */
			d__1 = xr;
/* Computing 2nd power */
			d__2 = yr;
/* Computing 2nd power */
			d__3 = zr;
			r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
			r__ = sqrt(r2);
			de = 0.;
			if (ri + r__ < sk) {
			    uik = sk - r__;
/* Computing 4th power */
			    d__1 = uik, d__1 *= d__1;
			    de = -12.566370614359172 / (d__1 * d__1);
			}
			if (ri + r__ < sk) {
			    lik = sk - r__;
/* Computing 4th power */
			    d__1 = lik, d__1 *= d__1;
			    de += (sk2 - sk * 4. * r__ + r2 * 17.) * 
				    .78539816339744828 / (r2 * (d__1 * d__1));
			} else if (r__ < ri + sk) {
			    lik = ri;
/* Computing 4th power */
			    d__1 = lik, d__1 *= d__1;
			    de += (ri * 2. * ri - sk2 - r2) * 
				    .78539816339744828 / (r2 * (d__1 * d__1));
			} else {
			    lik = r__ - sk;
/* Computing 4th power */
			    d__1 = lik, d__1 *= d__1;
			    de += (sk2 - sk * 4. * r__ + r2) * 
				    .78539816339744828 / (r2 * (d__1 * d__1));
			}
			uik = r__ + sk;
/* Computing 4th power */
			d__1 = uik, d__1 *= d__1;
			de -= (sk2 + sk * 4. * r__ + r2) * .78539816339744828 
				/ (r2 * (d__1 * d__1));
			dbr = term * de / r__;
			de = dbr * solute_1.drb[i__ - 1];

/*     increment the overall permanent solvation derivatives */

			dedx = de * xr;
			dedy = de * yr;
			dedz = de * zr;
			des_ref(1, i__) = des_ref(1, i__) + dedx;
			des_ref(2, i__) = des_ref(2, i__) + dedy;
			des_ref(3, i__) = des_ref(3, i__) + dedz;
			des_ref(1, k) = des_ref(1, k) - dedx;
			des_ref(2, k) = des_ref(2, k) - dedy;
			des_ref(3, k) = des_ref(3, k) - dedz;

/*     increment the internal virial tensor components */

			vxx = xr * dedx;
			vyx = yr * dedx;
			vzx = zr * dedx;
			vyy = yr * dedy;
			vzy = zr * dedy;
			vzz = zr * dedz;
			vir_ref(1, 1) = vir_ref(1, 1) + vxx;
			vir_ref(2, 1) = vir_ref(2, 1) + vyx;
			vir_ref(3, 1) = vir_ref(3, 1) + vzx;
			vir_ref(1, 2) = vir_ref(1, 2) + vyx;
			vir_ref(2, 2) = vir_ref(2, 2) + vyy;
			vir_ref(3, 2) = vir_ref(3, 2) + vzy;
			vir_ref(1, 3) = vir_ref(1, 3) + vzx;
			vir_ref(2, 3) = vir_ref(2, 3) + vzy;
			vir_ref(3, 3) = vir_ref(3, 3) + vzz;

/*     increment the polarization solvation derivatives */

			if (potent_1.use_mpole__ || potent_1.use_polar__) {
			    de = dbr * solute_1.drbp[i__ - 1];
			    dedx = de * xr;
			    dedy = de * yr;
			    dedz = de * zr;
			    des_ref(1, i__) = des_ref(1, i__) + dedx;
			    des_ref(2, i__) = des_ref(2, i__) + dedy;
			    des_ref(3, i__) = des_ref(3, i__) + dedz;
			    des_ref(1, k) = des_ref(1, k) - dedx;
			    des_ref(2, k) = des_ref(2, k) - dedy;
			    des_ref(3, k) = des_ref(3, k) - dedz;

/*     increment the internal virial tensor components */

			    vxx = xr * dedx;
			    vyx = yr * dedx;
			    vzx = zr * dedx;
			    vyy = yr * dedy;
			    vzy = zr * dedy;
			    vzz = zr * dedz;
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
		}
	    }
	}

/*     get Born radius chain rule components for the ACE method */

    } else if (s_cmp(solute_1.borntyp, "ACE", (ftnlen)8, (ftnlen)3) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    it = atmtyp_1.class__[i__ - 1];
	    vi = solute_1.vsolv[i__ - 1];
	    rbi = solute_1.rborn[i__ - 1];
	    rbi2 = rbi * rbi;
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		if (k != i__) {
		    xr = xi - atoms_1.x[k - 1];
		    yr = yi - atoms_1.y[k - 1];
		    zr = zi - atoms_1.z__[k - 1];
		    kt = atmtyp_1.class__[k - 1];
		    vk = solute_1.vsolv[k - 1];
		    s2ik = 1. / s2ace_ref(it, kt);
		    ws2 = wace_ref(it, kt) * s2ik;
/* Computing 4th power */
		    d__1 = uace_ref(it, kt), d__1 *= d__1;
		    uik4 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = xr;
/* Computing 2nd power */
		    d__2 = yr;
/* Computing 2nd power */
		    d__3 = zr;
		    r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
		    r__ = sqrt(r2);
		    r3 = r2 * r__;
		    r4 = r2 * r2;
		    r6 = r2 * r4;
		    rusum = r4 + uik4;
		    ratio = r3 / rusum;
		    expterm = exp(-r2 * s2ik);
		    de1 = r__ * -4. * ws2 * expterm;
/* Computing 2nd power */
		    d__1 = rusum;
		    de2 = r2 * 3. / rusum - r6 * 4. / (d__1 * d__1);
/* Computing 3rd power */
		    d__1 = ratio;
		    de = solute_1.drb[i__ - 1] * rbi2 * (de1 + vk * (d__1 * (
			    d__1 * d__1)) * de2 / 3.141592653589793238) / r__;

/*     increment the overall continuum solvation derivatives */

		    dedx = de * xr;
		    dedy = de * yr;
		    dedz = de * zr;
		    des_ref(1, i__) = des_ref(1, i__) + dedx;
		    des_ref(2, i__) = des_ref(2, i__) + dedy;
		    des_ref(3, i__) = des_ref(3, i__) + dedz;
		    des_ref(1, k) = des_ref(1, k) - dedx;
		    des_ref(2, k) = des_ref(2, k) - dedy;
		    des_ref(3, k) = des_ref(3, k) - dedz;

/*     increment the internal virial tensor components */

		    vxx = xr * dedx;
		    vyx = yr * dedx;
		    vzx = zr * dedx;
		    vyy = yr * dedy;
		    vzy = zr * dedy;
		    vzz = zr * dedz;
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
	}
    }
    return 0;
} /* born1_ */

#undef s2ace_ref
#undef wace_ref
#undef uace_ref
#undef vir_ref
#undef des_ref
#undef i13_ref
#undef i12_ref


