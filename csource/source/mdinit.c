/* mdinit.f -- translated by f2c (version 20050501).
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
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

struct {
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

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
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

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
    doublereal xrb[25000], yrb[25000], zrb[25000], rbc[6000]	/* was [6][
	    1000] */;
    logical use_rigid__;
} rigid_;

#define rigid_1 rigid_

struct {
    doublereal krat[25000];
    integer nrat, nratx, irat[50000]	/* was [2][25000] */, iratx[25000], 
	    kratx[25000];
    logical ratimage[25000], use_rattle__;
} shake_;

#define shake_1 shake_

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
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine mdinit  --  initialize a dynamics trajectory  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "mdinit" initializes the velocities and accelerations */
/*     for a molecular dynamics trajectory, including restarts */


/* Subroutine */ int mdinit_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 MDINIT  --  Warning, Mass of Group\002,i"
	    "6,\002 Set to 1.0 for Dynamics\002)";
    static char fmt_30[] = "(/,\002 MDINIT  --  Warning, Mass of Atom\002,"
	    "i6,\002 Set to 1.0 for Dynamics\002)";
    static char fmt_40[] = "(/,\002 Number of Degrees of Freedom for Dynamic"
	    "s :\002,i10)";
    static char fmt_50[] = "(/,\002 MDINIT  --  No Degrees of Freedom for Dy"
	    "namics\002)";

    /* System generated locals */
    address a__1[2], a__2[3];
    integer i__1, i__2, i__3[2], i__4[3];
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_inqu(inlist *), f_open(olist *), f_rew(alist *), f_clos(cllist *
	    );

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    extern integer freeunit_(void);
    static doublereal e;
    static integer i__, j;
    static doublereal vec[3];
    static char ext[7];
    static integer idyn, size, next, lext;
    extern /* Subroutine */ int fatal_(void);
    static doublereal speed;
    static logical exist;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static doublereal derivs[75000]	/* was [3][25000] */;
    static char string[120];
    extern /* Subroutine */ int ranvec_(doublereal *), mdrest_(void), 
	    readdyn_(integer *), lattice_(void);
    static char dynfile[120];
    extern /* Subroutine */ int grpline_(void), numeral_(integer *, char *, 
	    integer *, ftnlen);
    extern doublereal maxwell_(doublereal *, doublereal *);
    extern /* Subroutine */ int getword_(char *, char *, integer *, ftnlen, 
	    ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___6 = { 1, string, 1, 0, 120, 1 };
    static icilist io___7 = { 1, string, 1, 0, 120, 1 };
    static icilist io___8 = { 1, string, 1, 0, 120, 1 };
    static icilist io___9 = { 1, string, 1, 0, 120, 1 };
    static icilist io___10 = { 1, string, 1, 0, 120, 1 };
    static icilist io___11 = { 1, string, 1, 0, 120, 1 };
    static icilist io___13 = { 1, string, 1, 0, 120, 1 };
    static icilist io___14 = { 1, string, 1, 0, 120, 1 };
    static icilist io___15 = { 1, string, 1, 0, 120, 1 };
    static icilist io___16 = { 1, string, 1, 0, 120, 1 };
    static cilist io___17 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_50, 0 };



#define a_ref(a_1,a_2) moldyn_1.a[(a_2)*3 + a_1 - 4]
#define v_ref(a_1,a_2) moldyn_1.v[(a_2)*3 + a_1 - 4]
#define lm_ref(a_1,a_2) rgddyn_1.lm[(a_2)*3 + a_1 - 4]
#define vcm_ref(a_1,a_2) rgddyn_1.vcm[(a_2)*3 + a_1 - 4]
#define wcm_ref(a_1,a_2) rgddyn_1.wcm[(a_2)*3 + a_1 - 4]
#define aold_ref(a_1,a_2) moldyn_1.aold[(a_2)*3 + a_1 - 4]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]
#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1 - 4]
#define keyline_ref(a_0,a_1) &keys_1.keyline[(a_1)*120 + a_0 - 120]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  files.i  --  name and number of current structure files  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     nprior     number of previously existing cycle files */
/*     ldir       length in characters of the directory name */
/*     leng       length in characters of the base filename */
/*     filename   base filename used by default for all files */
/*     outfile    output filename used for intermediate results */




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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  keys.i  --  contents of current keyword parameter file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nkey      number of nonblank lines in the keyword file */
/*     keyline   contents of each individual keyword file line */




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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  rigid.i  --  rigid body coordinates for atom groups  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     xrb         rigid body reference x-coordinate for each atom */
/*     yrb         rigid body reference y-coordinate for each atom */
/*     zrb         rigid body reference z-coordinate for each atom */
/*     rbc         current rigid body coordinates for each group */
/*     use_rigid   flag to mark use of rigid body coordinate system */




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




/*     set default parameters for the dynamics trajectory */

    s_copy(mdstuf_1.integrate, "BEEMAN", (ftnlen)10, (ftnlen)6);
    mdstuf_1.nfree = 0;
    mdstuf_1.velsave = FALSE_;
    mdstuf_1.frcsave = FALSE_;
    mdstuf_1.uindsave = FALSE_;
    stodyn_1.friction = 91.;
    stodyn_1.use_sdarea__ = FALSE_;
    inform_1.iprint = 100;

/*     set default values for temperature and pressure control */

    s_copy(bath_1.thermostat, "BUSSI", (ftnlen)11, (ftnlen)5);
    bath_1.tautemp = .2;
    bath_1.collide = .1;
    for (i__ = 1; i__ <= 2; ++i__) {
	bath_1.xnh[i__ - 1] = 0.;
	bath_1.vnh[i__ - 1] = 0.;
	bath_1.qnh[i__ - 1] = 0.;
	bath_1.gnh[i__ - 1] = 0.;
    }
    s_copy(bath_1.barostat, "BERENDSEN", (ftnlen)10, (ftnlen)9);
    bath_1.anisotrop = FALSE_;
    bath_1.taupres = 2.;
    bath_1.compress = 4.6e-5;
    bath_1.voltrial = 20;
    bath_1.volmove = 100.;
    s_copy(bath_1.volscale, "ATOMIC", (ftnlen)9, (ftnlen)6);

/*     check for keywords containing any altered parameters */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "INTEGRATE ", (ftnlen)10, (ftnlen)10) == 0) {
	    getword_(record, mdstuf_1.integrate, &next, (ftnlen)120, (ftnlen)
		    10);
	    upcase_(mdstuf_1.integrate, (ftnlen)10);
	} else if (s_cmp(keyword, "DEGREES-FREEDOM ", (ftnlen)16, (ftnlen)16) 
		== 0) {
	    i__2 = s_rsli(&io___6);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&mdstuf_1.nfree, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "SAVE-VELOCITY ", (ftnlen)14, (ftnlen)14) ==
		 0) {
	    mdstuf_1.velsave = TRUE_;
	} else if (s_cmp(keyword, "SAVE-FORCE ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    mdstuf_1.frcsave = TRUE_;
	} else if (s_cmp(keyword, "SAVE-INDUCED ", (ftnlen)13, (ftnlen)13) == 
		0) {
	    mdstuf_1.uindsave = TRUE_;
	} else if (s_cmp(keyword, "FRICTION ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___7);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&stodyn_1.friction, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "FRICTION-SCALING ", (ftnlen)17, (ftnlen)17)
		 == 0) {
	    stodyn_1.use_sdarea__ = TRUE_;
	} else if (s_cmp(keyword, "PRINTOUT ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___8);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&inform_1.iprint, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "THERMOSTAT ", (ftnlen)11, (ftnlen)11) == 0)
		 {
	    getword_(record, bath_1.thermostat, &next, (ftnlen)120, (ftnlen)
		    11);
	    upcase_(bath_1.thermostat, (ftnlen)11);
	} else if (s_cmp(keyword, "TAU-TEMPERATURE ", (ftnlen)16, (ftnlen)16) 
		== 0) {
	    i__2 = s_rsli(&io___9);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&bath_1.tautemp, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "COLLISION ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    i__2 = s_rsli(&io___10);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&bath_1.collide, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "NOSE-MASS ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    i__2 = s_rsli(&io___11);
	    if (i__2 != 0) {
		goto L10;
	    }
	    for (j = 1; j <= 2; ++j) {
		i__2 = do_lio(&c__5, &c__1, (char *)&bath_1.qnh[j - 1], (
			ftnlen)sizeof(doublereal));
		if (i__2 != 0) {
		    goto L10;
		}
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "BAROSTAT ", (ftnlen)9, (ftnlen)9) == 0) {
	    getword_(record, bath_1.barostat, &next, (ftnlen)120, (ftnlen)10);
	    upcase_(bath_1.barostat, (ftnlen)10);
	} else if (s_cmp(keyword, "ANISO-PRESSURE ", (ftnlen)15, (ftnlen)15) 
		== 0) {
	    bath_1.anisotrop = TRUE_;
	} else if (s_cmp(keyword, "TAU-PRESSURE ", (ftnlen)13, (ftnlen)13) == 
		0) {
	    i__2 = s_rsli(&io___13);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&bath_1.taupres, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "COMPRESS ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___14);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&bath_1.compress, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "VOLUME-TRIAL ", (ftnlen)13, (ftnlen)13) == 
		0) {
	    i__2 = s_rsli(&io___15);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&bath_1.voltrial, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "VOLUME-MOVE ", (ftnlen)12, (ftnlen)12) == 
		0) {
	    i__2 = s_rsli(&io___16);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&bath_1.volmove, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "VOLUME-SCALE ", (ftnlen)13, (ftnlen)13) == 
		0) {
	    getword_(record, bath_1.volscale, &next, (ftnlen)120, (ftnlen)9);
	    upcase_(bath_1.volscale, (ftnlen)9);
	}
L10:
	;
    }

/*     make sure all atoms or groups have a nonzero mass */

    if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 0) {
	i__1 = group_1.ngrp;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (group_1.grpmass[i__ - 1] <= 0.) {
		group_1.grpmass[i__ - 1] = 1.;
		if (igrp_ref(1, i__) <= igrp_ref(2, i__)) {
		    molcul_1.totmass += 1.;
		    io___17.ciunit = iounit_1.iout;
		    s_wsfe(&io___17);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    e_wsfe();
		}
	    }
	}
    } else {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (usage_1.use[i__ - 1] && atmtyp_1.mass[i__ - 1] <= 0.) {
		atmtyp_1.mass[i__ - 1] = 1.;
		molcul_1.totmass += 1.;
		io___18.ciunit = iounit_1.iout;
		s_wsfe(&io___18);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	}
    }

/*     enforce use of velocity Verlet with Andersen thermostat */

    if (s_cmp(bath_1.thermostat, "ANDERSEN", (ftnlen)11, (ftnlen)8) == 0) {
	if (s_cmp(mdstuf_1.integrate, "BEEMAN", (ftnlen)10, (ftnlen)6) == 0) {
	    s_copy(mdstuf_1.integrate, "VERLET", (ftnlen)10, (ftnlen)6);
	}
    }

/*     set masses for the thermostats in a Nose-Hoover chain */

    if (s_cmp(bath_1.thermostat, "NOSE-HOOVER", (ftnlen)11, (ftnlen)11) == 0) 
	    {
	if (bath_1.qnh[0] == 0.) {
	    bath_1.qnh[0] = .1;
	}
	for (j = 2; j <= 2; ++j) {
	    if (bath_1.qnh[j - 1] == 0.) {
		bath_1.qnh[j - 1] = bath_1.qnh[j - 2];
	    }
	}
    }

/*     set the number of degrees of freedom for the system */

    if (mdstuf_1.nfree == 0) {
	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    grpline_();
	    mdstuf_1.nfree = group_1.ngrp * 6;
	    i__1 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		size = igrp_ref(2, i__) - igrp_ref(1, i__) + 1;
		if (size == 1) {
		    mdstuf_1.nfree += -3;
		}
		if (rgddyn_1.linear[i__ - 1]) {
		    --mdstuf_1.nfree;
		}
	    }
	} else {
	    mdstuf_1.nfree = usage_1.nuse * 3;
	}
	if (shake_1.use_rattle__) {
	    mdstuf_1.nfree -= shake_1.nrat;
	    i__1 = shake_1.nratx;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		mdstuf_1.nfree -= shake_1.kratx[i__ - 1];
	    }
	}
	if (bath_1.isothermal && s_cmp(mdstuf_1.integrate, "STOCHASTIC", (
		ftnlen)10, (ftnlen)10) != 0 && s_cmp(bath_1.thermostat, "AND"
		"ERSEN", (ftnlen)11, (ftnlen)8) != 0) {
	    if (bound_1.use_bounds__) {
		mdstuf_1.nfree += -3;
	    } else {
		mdstuf_1.nfree += -6;
	    }
	}
    }

/*     check for a nonzero number of degrees of freedom */

    if (mdstuf_1.nfree < 0) {
	mdstuf_1.nfree = 0;
    }
    if (inform_1.debug) {
	io___20.ciunit = iounit_1.iout;
	s_wsfe(&io___20);
	do_fio(&c__1, (char *)&mdstuf_1.nfree, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (mdstuf_1.nfree == 0) {
	io___21.ciunit = iounit_1.iout;
	s_wsfe(&io___21);
	e_wsfe();
	fatal_();
    }

/*     try to restart using prior velocities and accelerations */

/* Writing concatenation */
    i__3[0] = files_1.leng, a__1[0] = files_1.filename;
    i__3[1] = 4, a__1[1] = ".dyn";
    s_cat(dynfile, a__1, i__3, &c__2, (ftnlen)120);
    version_(dynfile, "old", (ftnlen)120, (ftnlen)3);
    ioin__1.inerr = 0;
    ioin__1.infilen = 120;
    ioin__1.infile = dynfile;
    ioin__1.inex = &exist;
    ioin__1.inopen = 0;
    ioin__1.innum = 0;
    ioin__1.innamed = 0;
    ioin__1.inname = 0;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    if (exist) {
	idyn = freeunit_();
	o__1.oerr = 0;
	o__1.ounit = idyn;
	o__1.ofnmlen = 120;
	o__1.ofnm = dynfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = idyn;
	f_rew(&al__1);
	readdyn_(&idyn);
	cl__1.cerr = 0;
	cl__1.cunit = idyn;
	cl__1.csta = 0;
	f_clos(&cl__1);
	lattice_();

/*     set translational velocities for rigid body dynamics */

    } else if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) 
	    == 0) {
	gradient_(&e, derivs);
	i__1 = group_1.ngrp;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    speed = maxwell_(&group_1.grpmass[i__ - 1], &bath_1.kelvin);
	    ranvec_(vec);
	    for (j = 1; j <= 3; ++j) {
		vcm_ref(j, i__) = speed * vec[j - 1];
		wcm_ref(j, i__) = 0.;
		lm_ref(j, i__) = 0.;
	    }
	}
	if (usage_1.nuse == atoms_1.n) {
	    mdrest_();
	}

/*     set velocities and accelerations for Cartesian dynamics */

    } else {
	gradient_(&e, derivs);
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (usage_1.use[i__ - 1]) {
		speed = maxwell_(&atmtyp_1.mass[i__ - 1], &bath_1.kelvin);
		ranvec_(vec);
		for (j = 1; j <= 3; ++j) {
		    v_ref(j, i__) = speed * vec[j - 1];
		    a_ref(j, i__) = derivs_ref(j, i__) * -418.4 / 
			    atmtyp_1.mass[i__ - 1];
		    aold_ref(j, i__) = a_ref(j, i__);
		}
	    } else {
		for (j = 1; j <= 3; ++j) {
		    v_ref(j, i__) = 0.;
		    a_ref(j, i__) = 0.;
		    aold_ref(j, i__) = 0.;
		}
	    }
	}
	if (usage_1.nuse == atoms_1.n) {
	    mdrest_();
	}
    }

/*     check for any prior dynamics coordinate sets */

    i__ = 0;
    exist = TRUE_;
    while(exist) {
	++i__;
	lext = 3;
	numeral_(&i__, ext, &lext, (ftnlen)7);
/* Writing concatenation */
	i__4[0] = files_1.leng, a__2[0] = files_1.filename;
	i__4[1] = 1, a__2[1] = ".";
	i__4[2] = lext, a__2[2] = ext;
	s_cat(dynfile, a__2, i__4, &c__3, (ftnlen)120);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = dynfile;
	ioin__1.inex = &exist;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! exist && i__ < 100) {
	    lext = 2;
	    numeral_(&i__, ext, &lext, (ftnlen)7);
/* Writing concatenation */
	    i__4[0] = files_1.leng, a__2[0] = files_1.filename;
	    i__4[1] = 1, a__2[1] = ".";
	    i__4[2] = lext, a__2[2] = ext;
	    s_cat(dynfile, a__2, i__4, &c__3, (ftnlen)120);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = dynfile;
	    ioin__1.inex = &exist;
	    ioin__1.inopen = 0;
	    ioin__1.innum = 0;
	    ioin__1.innamed = 0;
	    ioin__1.inname = 0;
	    ioin__1.inacc = 0;
	    ioin__1.inseq = 0;
	    ioin__1.indir = 0;
	    ioin__1.infmt = 0;
	    ioin__1.inform = 0;
	    ioin__1.inunf = 0;
	    ioin__1.inrecl = 0;
	    ioin__1.innrec = 0;
	    ioin__1.inblank = 0;
	    f_inqu(&ioin__1);
	}
	if (! exist && i__ < 10) {
	    lext = 1;
	    numeral_(&i__, ext, &lext, (ftnlen)7);
/* Writing concatenation */
	    i__4[0] = files_1.leng, a__2[0] = files_1.filename;
	    i__4[1] = 1, a__2[1] = ".";
	    i__4[2] = lext, a__2[2] = ext;
	    s_cat(dynfile, a__2, i__4, &c__3, (ftnlen)120);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = dynfile;
	    ioin__1.inex = &exist;
	    ioin__1.inopen = 0;
	    ioin__1.innum = 0;
	    ioin__1.innamed = 0;
	    ioin__1.inname = 0;
	    ioin__1.inacc = 0;
	    ioin__1.inseq = 0;
	    ioin__1.indir = 0;
	    ioin__1.infmt = 0;
	    ioin__1.inform = 0;
	    ioin__1.inunf = 0;
	    ioin__1.inrecl = 0;
	    ioin__1.innrec = 0;
	    ioin__1.inblank = 0;
	    f_inqu(&ioin__1);
	}
    }
    files_1.nprior = i__ - 1;
    return 0;
} /* mdinit_ */

#undef keyline_ref
#undef derivs_ref
#undef igrp_ref
#undef aold_ref
#undef wcm_ref
#undef vcm_ref
#undef lm_ref
#undef v_ref
#undef a_ref


