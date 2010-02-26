/* initial.f -- translated by f2c (version 20050501).
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
    doublereal wfit[25000];
    integer nfit, ifit[50000]	/* was [2][25000] */;
} align_;

#define align_1 align_

struct {
    integer narg;
    logical listarg[21];
    char arg[2520];
} argue_;

#define argue_1 argue_

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
    doublereal xcell, ycell, zcell, xcell2, ycell2, zcell2;
    integer ncell, icell[30000]	/* was [3][10000] */;
} cell_;

#define cell_1 cell_

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
    doublereal stpmin, stpmax, cappa, slpmax, angmax;
    integer intmax;
} linmin_;

#define linmin_1 linmin_

struct {
    doublereal fctmin, hguess;
    integer maxiter, nextiter;
} minima_;

#define minima_1 minima_

struct {
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

struct {
    doublereal lambda, vlambda, clambda, dlambda, mlambda, plambda;
    integer nmut, imut[25000], type0[25000], type1[25000], class0[25000], 
	    class1[25000];
    logical mut[25000];
} mutant_;

#define mutant_1 mutant_

struct {
    doublereal lbuffer, lbuf2, vbuf2, cbuf2, mbuf2;
    integer nvlst[25000], vlst[45000000]	/* was [1800][25000] */, 
	    nelst[25000], elst[30000000]	/* was [1200][25000] */;
    logical dovlst, doclst, domlst;
} neigh_;

#define neigh_1 neigh_

struct {
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

struct {
    integer nprm;
    char prmline[3000000];
} params_;

#define params_1 params_

struct {
    doublereal xpdb[25000], ypdb[25000], zpdb[25000];
    integer npdb, resnum[25000], npdb12[25000], ipdb12[200000]	/* was [8][
	    25000] */, pdblist[25000];
    char pdbtyp[150000], atmnam[100000], resnam[75000], chntyp[20], altsym[1],
	     instyp[20];
} pdb_;

#define pdb_1 pdb_

struct {
    doublereal tiny, small, huge__;
} precis_;

#define precis_1 precis_

struct {
    doublereal xrb[25000], yrb[25000], zrb[25000], rbc[6000]	/* was [6][
	    1000] */;
    logical use_rigid__;
} rigid_;

#define rigid_1 rigid_

struct {
    doublereal scale[75000];
    logical set_scale__;
} scales_;

#define scales_1 scales_

struct {
    integer nseq, nchain, ichain[20000]	/* was [2][10000] */, seqtyp[10000];
    char seq[30000], chnnam[10000];
} sequen_;

#define sequen_1 sequen_

struct {
    integer runtyp, cstep;
    doublereal cdt, cenergy, cdx[25000], cdy[25000], cdz[25000];
    logical use_socket__, skt_init__, skt_close__;
} socket_;

#define socket_1 socket_

struct {
    doublereal m2[25000], deform, difft, diffv, diffc;
    logical use_smooth__, use_dem__, use_gda__, use_tophat__, use_stophat__;
} warp_;

#define warp_1 warp_

struct {
    integer nadd, iadd[50000]	/* was [2][25000] */, ndel, idel[50000]	/* 
	    was [2][25000] */;
} zclose_;

#define zclose_1 zclose_

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine initial  --  initial values and program setup  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "initial" sets up original values for some parameters */
/*     and variables that might not otherwise get initialized */


/* Subroutine */ int initial_(void)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int promo_(void), command_(void);
    extern doublereal precise_(integer *);
    extern /* Subroutine */ int initres_(void);



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  align.i  --  information for superposition of structures  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     wfit    weights assigned to atom pairs during superposition */
/*     nfit    number of atoms to use in superimposing two structures */
/*     ifit    atom numbers of pairs of atoms to be superimposed */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  argue.i  --  command line arguments at program startup  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     maxarg    maximum number of command line arguments */

/*     narg      number of command line arguments to the program */
/*     listarg   flag to mark available command line arguments */
/*     arg       strings containing the command line arguments */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cell.i  --  periodic boundaries using replicated cells  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xcell    length of the a-axis of the complete replicated cell */
/*     ycell    length of the b-axis of the complete replicated cell */
/*     zcell    length of the c-axis of the complete replicated cell */
/*     xcell2   half the length of the a-axis of the replicated cell */
/*     ycell2   half the length of the b-axis of the replicated cell */
/*     zcell2   half the length of the c-axis of the replicated cell */
/*     ncell    total number of cell replicates for periodic boundaries */
/*     icell    offset along axes for each replicate periodic cell */




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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  linmin.i  --  parameters for line search minimization  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     stpmin   minimum step length in current line search direction */
/*     stpmax   maximum step length in current line search direction */
/*     cappa    stringency of line search (0=tight < cappa < 1=loose) */
/*     slpmax   projected gradient above which stepsize is reduced */
/*     angmax   maximum angle between search direction and -gradient */
/*     intmax   maximum number of interpolations during line search */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  minima.i  --  general parameters for minimizations  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     fctmin    value below which function is deemed optimized */
/*     hguess    initial value for the H-matrix diagonal elements */
/*     maxiter   maximum number of iterations during optimization */
/*     nextiter  iteration number to use for the first iteration */




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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mutant.i  --  parameters for free energy perturbation  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     lambda    generic weighting of the initial and final states */
/*     vlambda   weighting of initial and final states for vdw */
/*     clambda   weighting of initial and final states for charges */
/*     dlambda   weighting of initial and final states for dipoles */
/*     mlambda   weighting of initial and final states for multipoles */
/*     plambda   weighting of initial and final states for polarization */
/*     nmut      number of atoms mutated from initial to final state */
/*     imut      atomic sites differing in initial and final state */
/*     type0     atom type of each atom in the initial state system */
/*     type1     atom type of each atom in the final state system */
/*     class0    atom class of each atom in the initial state system */
/*     class1    atom class of each atom in the final state system */
/*     mut       true if an atom is to be mutated, false otherwise */




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  neigh.i  --  pairwise neighbor list indices and storage  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     lbuffer     width of the neighbor list buffer region */
/*     lbuf2       motion squared needed to trigger list rebuild */
/*     vbuf2       square of vdw cutoff plus neighbor list buffer */
/*     cbuf2       square of charge cutoff plus neighbor list buffer */
/*     mbuf2       square of multipole cutoff plus neighbor list buffer */
/*     nvlst       number of sites in list for each vdw site */
/*     vlst        site numbers in neighbor list of each vdw site */
/*     nelst       number of sites in list for each electrostatic site */
/*     elst        site numbers in list of each electrostatic site */
/*     dovlst      logical flag to rebuild vdw neighbor list */
/*     doclst      logical flag to rebuild charge neighbor list */
/*     domlst      logical flag to rebuild multipole neighbor list */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  output.i  --  control of coordinate output file format  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     archive    logical flag to save structures in an archive */
/*     noversion  logical flag governing use of filename versions */
/*     overwrite  logical flag to overwrite intermediate files inplace */
/*     cyclesave  logical flag to mark use of numbered cycle files */
/*     coordtype  selects Cartesian, internal, rigid body or none */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  params.i  --  contents of force field parameter file  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     nprm      number of nonblank lines in the parameter file */
/*     prmline   contents of each individual parameter file line */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  pdb.i  --  definition of a Protein Data Bank structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xpdb      x-coordinate of each atom stored in PDB format */
/*     ypdb      y-coordinate of each atom stored in PDB format */
/*     zpdb      z-coordinate of each atom stored in PDB format */
/*     npdb      number of atoms stored in Protein Data Bank format */
/*     resnum    number of the residue to which each atom belongs */
/*     npdb12    number of atoms directly bonded to each CONECT atom */
/*     ipdb12    atom numbers of atoms connected to each CONECT atom */
/*     pdblist   list of the Protein Data Bank atom number of each atom */
/*     pdbtyp    Protein Data Bank record type assigned to each atom */
/*     atmnam    Protein Data Bank atom name assigned to each atom */
/*     resnam    Protein Data Bank residue name assigned to each atom */
/*     chntyp    string with PDB chain identifiers to be included */
/*     altsym    string with PDB alternate locations to be included */
/*     instyp    string with PDB insertion records to be included */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  precis.i  --  values of machine precision tolerances  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     tiny    the smallest positive floating point value */
/*     small   the smallest relative floating point spacing */
/*     huge    the largest relative floating point spacing */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  scales.i  --  parameter scale factors for optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     scale      multiplicative factor for each optimization parameter */
/*     set_scale  logical flag to show if scale factors have been set */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  sequen.i  --  sequence information for a biopolymer  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nseq     total number of residues in biopolymer sequences */
/*     nchain   number of separate biopolymer sequence chains */
/*     ichain   first and last residue in each biopolymer chain */
/*     seqtyp   residue type for each residue in the sequence */
/*     seq      three-letter code for each residue in the sequence */
/*     chnnam   one-letter identifier for each sequence chain */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2002  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  socket.i  --  control parameters for socket communication  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     runtyp      calculation type for passing socket information */
/*     cstep       current optimization or dynamics step number */
/*     cdt         current dynamics cumulative simulation time */
/*     cenergy     current potential energy from simulation */
/*     cdx         current gradient components along the x-axis */
/*     cdy         current gradient components along the y-axis */
/*     cdz         current gradient components along the z-axis */
/*     use_socket  logical flag governing use of external sockets */
/*     skt_init    logical flag to indicate socket initialization */
/*     skt_close   logical flag to indicate socket shutdown */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  zclose.i  --  ring openings and closures for Z-matrix  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     nadd   number of added bonds between Z-matrix atoms */
/*     iadd   numbers of the atom pairs defining added bonds */
/*     ndel   number of bonds between Z-matrix bonds to delete */
/*     idel   numbers of the atom pairs defining deleted bonds */




/*     number of atoms used in superposition */

    align_1.nfit = 0;

/*     number of command line arguments */

    argue_1.narg = 0;

/*     number of atoms in the system */

    atoms_1.n = 0;

/*     flags for temperature and pressure baths */

    bath_1.isothermal = FALSE_;
    bath_1.isobaric = FALSE_;

/*     flags for periodic boundaries */

    bound_1.use_bounds__ = FALSE_;
    bound_1.use_replica__ = FALSE_;
    bound_1.use_polymer__ = FALSE_;

/*     number of unit cell replicates */

    cell_1.ncell = 0;

/*     flag for use of atom groups */

    group_1.use_group__ = FALSE_;

/*     highest numbered previous cycle file */

    files_1.nprior = 0;

/*     information levels within the program */

    inform_1.verbose = FALSE_;
    inform_1.debug = FALSE_;
    inform_1.abort = FALSE_;

/*     default input/output unit numbers */

    iounit_1.input = 5;
    iounit_1.iout = 6;

/*     number of lines in the keyfile */

    keys_1.nkey = 0;

/*     default parameters used by line search */

    linmin_1.stpmin = 0.;
    linmin_1.stpmax = 0.;
    linmin_1.cappa = 0.;
    linmin_1.slpmax = 0.;
    linmin_1.angmax = 0.;
    linmin_1.intmax = 0;

/*     default parameters used by optimizations */

    minima_1.fctmin = 0.;
    minima_1.maxiter = 0;
    minima_1.nextiter = 0;
    inform_1.iprint = -1;
    inform_1.iwrite = -1;

/*     number of molecules in the system */

    molcul_1.nmol = 0;

/*     number of mutated atoms in the system */

    mutant_1.nmut = 0;

/*     flags for rebuilding of neighbor lists */

    neigh_1.dovlst = TRUE_;
    neigh_1.doclst = TRUE_;
    neigh_1.domlst = TRUE_;

/*     type of coordinates file */

    s_copy(output_1.coordtype, "NONE", (ftnlen)9, (ftnlen)4);

/*     number of lines in the parameter file */

    params_1.nprm = 0;

/*     number of atoms in Protein Data Bank format */

    pdb_1.npdb = 0;

/*     flag for use of rigid bodies */

    rigid_1.use_rigid__ = FALSE_;

/*     flag to show setting of optimization scale factors */

    scales_1.set_scale__ = FALSE_;

/*     number of residues and chains in biopolymer sequence */

    sequen_1.nseq = 0;
    sequen_1.nchain = 0;

/*     flags for external Java socket communication */

    socket_1.skt_init__ = FALSE_;
    socket_1.use_socket__ = FALSE_;

/*     flags for potential energy smoothing */

    warp_1.use_smooth__ = FALSE_;
    warp_1.use_dem__ = FALSE_;
    warp_1.use_gda__ = FALSE_;
    warp_1.use_tophat__ = FALSE_;
    warp_1.use_stophat__ = FALSE_;

/*     number of bonds added or deleted from Z-matrix */

    zclose_1.nadd = 0;
    zclose_1.ndel = 0;

/*     display program info and copyright notice */

    promo_();

/*     names of biopolymer residue types */

    initres_();

/*     determine a set of machine precision values */

    precis_1.tiny = precise_(&c__1);
    precis_1.small = precise_(&c__2);
    precis_1.huge__ = precise_(&c__3);

/*     get any command line arguments to the program */

    command_();
    return 0;
} /* initial_ */

