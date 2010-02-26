/* mdsave.f -- translated by f2c (version 20050501).
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
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

struct {
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

struct {
    doublereal polarity[25000], thole[25000], pdamp[25000], uind[75000]	/* 
	    was [3][25000] */, uinp[75000]	/* was [3][25000] */, uinds[
	    75000]	/* was [3][25000] */, uinps[75000]	/* was [3][
	    25000] */;
    integer npolar;
} polar_;

#define polar_1 polar_

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
    doublereal vcm[3000]	/* was [3][1000] */, wcm[3000]	/* was [3][
	    1000] */, lm[3000]	/* was [3][1000] */, vc[3000]	/* was [3][
	    1000] */, wc[3000]	/* was [3][1000] */;
    logical linear[1000];
} rgddyn_;

#define rgddyn_1 rgddyn_

struct {
    integer runtyp, cstep;
    doublereal cdt, cenergy, cdx[25000], cdy[25000], cdz[25000];
    logical use_socket__, skt_init__, skt_close__;
} socket_;

#define socket_1 socket_

struct {
    integer ltitle;
    char title[120];
} titles_;

#define titles_1 titles_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine mdsave  --  save trajectory and restart files  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "mdsave" writes molecular dynamics trajectory snapshots and */
/*     auxiliary files with velocity, force or induced dipole data; */
/*     also checks for user requested termination of a simulation */


/* Subroutine */ int mdsave_(integer *istep, doublereal *dt, doublereal *epot)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Instantaneous Values for Frame saved a"
	    "t\002,i10,\002 Dynamics Steps\002)";
    static char fmt_20[] = "(/,\002 Current Time\002,8x,f15.4,\002 Picosec"
	    "ond\002)";
    static char fmt_30[] = "(\002 Current Potential\002,3x,f15.4,\002 Kcal/m"
	    "ole\002)";
    static char fmt_40[] = "(\002 Lattice Lengths\002,6x,3f14.6)";
    static char fmt_50[] = "(\002 Lattice Angles\002,7x,3f14.6)";
    static char fmt_60[] = "(\002 Frame Number\002,13x,i10)";
    static char fmt_70[] = "(\002 Coordinate File\002,12x,a)";
    static char fmt_80[] = "(i6,2x,a)";
    static char fmt_90[] = "(i6,3x,d13.6,3x,d13.6,3x,d13.6)";
    static char fmt_100[] = "(i6,3x,d13.6,3x,d13.6,3x,d13.6)";
    static char fmt_110[] = "(i6,2x,a)";
    static char fmt_120[] = "(i6,2x,a3,3x,d13.6,3x,d13.6,3x,d13.6)";
    static char fmt_130[] = "(\002 Velocity File\002,15x,a)";
    static char fmt_140[] = "(i6,2x,a)";
    static char fmt_150[] = "(i6,2x,a3,3x,d13.6,3x,d13.6,3x,d13.6)";
    static char fmt_160[] = "(\002 Force Vector File\002,11x,a)";
    static char fmt_170[] = "(i6,2x,a)";
    static char fmt_180[] = "(i6,2x,a3,3f12.6)";
    static char fmt_190[] = "(\002 Induced Dipole File\002,10x,a)";
    static char fmt_200[] = "(/,\002 MDSAVE  --  Dynamics Calculation Endin"
	    "g\002,\002 due to User Request\002)";
    static char fmt_210[] = "()";

    /* System generated locals */
    address a__1[3], a__2[4], a__3[2];
    integer i__1[3], i__2[4], i__3, i__4[2];
    doublereal d__1;
    olist o__1;
    cllist cl__1;
    inlist ioin__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_inqu(inlist *), f_open(olist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_clos(cllist *), s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern integer freeunit_(void), trimtext_(char *, ftnlen);
    static integer i__, j, k;
    static doublereal wt;
    static char ext[7];
    static integer iend, iind, ifrc;
    static doublereal pico;
    static integer ivel, lext, ixyz;
    extern /* Subroutine */ int fatal_(void);
    static integer idump;
    static logical exist;
    extern /* Subroutine */ int sktdyn_(integer *, doublereal *, doublereal *)
	    , suffix_(char *, char *, ftnlen, ftnlen), prtdyn_(void), prtxyz_(
	    integer *);
    static char endfile[120], frcfile[120], indfile[120], velfile[120];
    extern /* Subroutine */ int openend_(integer *, char *, ftnlen), numeral_(
	    integer *, char *, integer *, ftnlen);
    static integer moddump;
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_210, 0 };



#define a_ref(a_1,a_2) moldyn_1.a[(a_2)*3 + a_1 - 4]
#define v_ref(a_1,a_2) moldyn_1.v[(a_2)*3 + a_1 - 4]
#define vcm_ref(a_1,a_2) rgddyn_1.vcm[(a_2)*3 + a_1 - 4]
#define wcm_ref(a_1,a_2) rgddyn_1.wcm[(a_2)*3 + a_1 - 4]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  polar.i  --  polarizabilities and induced dipole moments  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     polarity  dipole polarizability for each multipole site (Ang**3) */
/*     thole     Thole polarizability damping value for each site */
/*     pdamp     value of polarizability scale factor for each site */
/*     uind      induced dipole components at each multipole site */
/*     uinp      induced dipoles in field used for energy interactions */
/*     uinds     GK or PB induced dipoles at each multipole site */
/*     uinps     induced dipoles in field used for GK or PB energy */
/*     npolar    total number of polarizable sites in the system */




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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  titles.i  --  title for the current molecular system  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     ltitle   length in characters of the nonblank title string */
/*     title    title used to describe the current structure */




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




/*     send data via external socket communication if desired */

    if (! socket_1.skt_init__ || socket_1.use_socket__) {
	sktdyn_(istep, dt, epot);
    }

/*     check number of steps between trajectory file dumps */

    moddump = *istep % inform_1.iwrite;
    if (moddump != 0) {
	return 0;
    }

/*     get the sequence number of the current trajectory frame */

    idump = files_1.nprior + *istep / inform_1.iwrite;
    lext = 3;
    numeral_(&idump, ext, &lext, (ftnlen)7);

/*     print header for the instantaneous values at current step */

    pico = (doublereal) (*istep) * *dt;
    io___6.ciunit = iounit_1.iout;
    s_wsfe(&io___6);
    do_fio(&c__1, (char *)&(*istep), (ftnlen)sizeof(integer));
    e_wsfe();
    io___7.ciunit = iounit_1.iout;
    s_wsfe(&io___7);
    do_fio(&c__1, (char *)&pico, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___8.ciunit = iounit_1.iout;
    s_wsfe(&io___8);
    do_fio(&c__1, (char *)&(*epot), (ftnlen)sizeof(doublereal));
    e_wsfe();
    if (bound_1.use_bounds__) {
	io___9.ciunit = iounit_1.iout;
	s_wsfe(&io___9);
	do_fio(&c__1, (char *)&boxes_1.xbox, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&boxes_1.ybox, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&boxes_1.zbox, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___10.ciunit = iounit_1.iout;
	s_wsfe(&io___10);
	do_fio(&c__1, (char *)&boxes_1.alpha, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&boxes_1.beta, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&boxes_1.gamma, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    io___11.ciunit = iounit_1.iout;
    s_wsfe(&io___11);
    do_fio(&c__1, (char *)&idump, (ftnlen)sizeof(integer));
    e_wsfe();

/*     update the information needed to restart the trajectory */

    prtdyn_();

/*     save coordinates to an archive or numbered structure file */

    ixyz = freeunit_();
    if (output_1.archive) {
	s_copy(xyzfile, files_1.filename, (ftnlen)120, files_1.leng);
	suffix_(xyzfile, "arc", (ftnlen)120, (ftnlen)3);
	version_(xyzfile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = xyzfile;
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
	    openend_(&ixyz, xyzfile, (ftnlen)120);
	} else {
	    o__1.oerr = 0;
	    o__1.ounit = ixyz;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = xyzfile;
	    o__1.orl = 0;
	    o__1.osta = "new";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	}
    } else {
/* Writing concatenation */
	i__1[0] = files_1.leng, a__1[0] = files_1.filename;
	i__1[1] = 1, a__1[1] = ".";
	i__1[2] = lext, a__1[2] = ext;
	s_cat(xyzfile, a__1, i__1, &c__3, (ftnlen)120);
	version_(xyzfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = ixyz;
	o__1.ofnmlen = 120;
	o__1.ofnm = xyzfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }
    prtxyz_(&ixyz);
    cl__1.cerr = 0;
    cl__1.cunit = ixyz;
    cl__1.csta = 0;
    f_clos(&cl__1);
    io___15.ciunit = iounit_1.iout;
    s_wsfe(&io___15);
    do_fio(&c__1, xyzfile, trimtext_(xyzfile, (ftnlen)120));
    e_wsfe();

/*     save the velocity vector components at the current step */

    if (mdstuf_1.velsave) {
	ivel = freeunit_();
	if (output_1.archive) {
	    s_copy(velfile, files_1.filename, (ftnlen)120, files_1.leng);
	    suffix_(velfile, "vel", (ftnlen)120, (ftnlen)3);
	    version_(velfile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = velfile;
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
		openend_(&ivel, velfile, (ftnlen)120);
	    } else {
		o__1.oerr = 0;
		o__1.ounit = ivel;
		o__1.ofnmlen = 120;
		o__1.ofnm = velfile;
		o__1.orl = 0;
		o__1.osta = "new";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
	    }
	} else {
/* Writing concatenation */
	    i__2[0] = files_1.leng, a__2[0] = files_1.filename;
	    i__2[1] = 1, a__2[1] = ".";
	    i__2[2] = lext, a__2[2] = ext;
	    i__2[3] = 1, a__2[3] = "v";
	    s_cat(velfile, a__2, i__2, &c__4, (ftnlen)120);
	    version_(velfile, "new", (ftnlen)120, (ftnlen)3);
	    o__1.oerr = 0;
	    o__1.ounit = ivel;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = velfile;
	    o__1.orl = 0;
	    o__1.osta = "new";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	}
	if (s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10, (ftnlen)9) == 
		0) {
	    io___18.ciunit = ivel;
	    s_wsfe(&io___18);
	    do_fio(&c__1, (char *)&group_1.ngrp, (ftnlen)sizeof(integer));
	    do_fio(&c__1, titles_1.title, titles_1.ltitle);
	    e_wsfe();
	    i__3 = group_1.ngrp;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		io___20.ciunit = ivel;
		s_wsfe(&io___20);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&vcm_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
		io___22.ciunit = ivel;
		s_wsfe(&io___22);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&wcm_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	} else {
	    io___23.ciunit = ivel;
	    s_wsfe(&io___23);
	    do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
	    do_fio(&c__1, titles_1.title, titles_1.ltitle);
	    e_wsfe();
	    i__3 = atoms_1.n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		io___24.ciunit = ivel;
		s_wsfe(&io___24);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&v_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	}
	cl__1.cerr = 0;
	cl__1.cunit = ivel;
	cl__1.csta = 0;
	f_clos(&cl__1);
	io___25.ciunit = iounit_1.iout;
	s_wsfe(&io___25);
	do_fio(&c__1, velfile, trimtext_(velfile, (ftnlen)120));
	e_wsfe();
    }

/*     save the force vector components for the current step */

    if (mdstuf_1.frcsave && s_cmp(mdstuf_1.integrate, "RIGIDBODY", (ftnlen)10,
	     (ftnlen)9) != 0) {
	ifrc = freeunit_();
	if (output_1.archive) {
	    s_copy(frcfile, files_1.filename, (ftnlen)120, files_1.leng);
	    suffix_(frcfile, "frc", (ftnlen)120, (ftnlen)3);
	    version_(frcfile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = frcfile;
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
		openend_(&ifrc, frcfile, (ftnlen)120);
	    } else {
		o__1.oerr = 0;
		o__1.ounit = ifrc;
		o__1.ofnmlen = 120;
		o__1.ofnm = frcfile;
		o__1.orl = 0;
		o__1.osta = "new";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
	    }
	} else {
/* Writing concatenation */
	    i__2[0] = files_1.leng, a__2[0] = files_1.filename;
	    i__2[1] = 1, a__2[1] = ".";
	    i__2[2] = lext, a__2[2] = ext;
	    i__2[3] = 1, a__2[3] = "f";
	    s_cat(frcfile, a__2, i__2, &c__4, (ftnlen)120);
	    version_(frcfile, "new", (ftnlen)120, (ftnlen)3);
	    o__1.oerr = 0;
	    o__1.ounit = ifrc;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = frcfile;
	    o__1.orl = 0;
	    o__1.osta = "new";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	}
	io___28.ciunit = ifrc;
	s_wsfe(&io___28);
	do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
	do_fio(&c__1, titles_1.title, titles_1.ltitle);
	e_wsfe();
	i__3 = atoms_1.n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    wt = atmtyp_1.mass[i__ - 1] / 418.4;
	    io___30.ciunit = ifrc;
	    s_wsfe(&io___30);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
	    for (j = 1; j <= 3; ++j) {
		d__1 = wt * a_ref(j, i__);
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();
	}
	cl__1.cerr = 0;
	cl__1.cunit = ifrc;
	cl__1.csta = 0;
	f_clos(&cl__1);
	io___31.ciunit = iounit_1.iout;
	s_wsfe(&io___31);
	do_fio(&c__1, frcfile, trimtext_(frcfile, (ftnlen)120));
	e_wsfe();
    }

/*     save the current induced dipole moment at each site */

    if (mdstuf_1.uindsave && potent_1.use_polar__) {
	iind = freeunit_();
	if (output_1.archive) {
	    s_copy(indfile, files_1.filename, (ftnlen)120, files_1.leng);
	    suffix_(indfile, "uind", (ftnlen)120, (ftnlen)4);
	    version_(indfile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = indfile;
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
		openend_(&iind, indfile, (ftnlen)120);
	    } else {
		o__1.oerr = 0;
		o__1.ounit = iind;
		o__1.ofnmlen = 120;
		o__1.ofnm = indfile;
		o__1.orl = 0;
		o__1.osta = "new";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
	    }
	} else {
/* Writing concatenation */
	    i__2[0] = files_1.leng, a__2[0] = files_1.filename;
	    i__2[1] = 1, a__2[1] = ".";
	    i__2[2] = lext, a__2[2] = ext;
	    i__2[3] = 1, a__2[3] = "u";
	    s_cat(indfile, a__2, i__2, &c__4, (ftnlen)120);
	    version_(indfile, "new", (ftnlen)120, (ftnlen)3);
	    o__1.oerr = 0;
	    o__1.ounit = iind;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = indfile;
	    o__1.orl = 0;
	    o__1.osta = "new";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	}
	io___34.ciunit = iind;
	s_wsfe(&io___34);
	do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
	do_fio(&c__1, titles_1.title, titles_1.ltitle);
	e_wsfe();
	i__3 = mpole_1.npole;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    if (polar_1.polarity[i__ - 1] != 0.) {
		k = mpole_1.ipole[i__ - 1];
		io___36.ciunit = iind;
		s_wsfe(&io___36);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, k), (ftnlen)3);
		for (j = 1; j <= 3; ++j) {
		    d__1 = uind_ref(j, i__) * 4.80321;
		    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		}
		e_wsfe();
	    }
	}
	cl__1.cerr = 0;
	cl__1.cunit = iind;
	cl__1.csta = 0;
	f_clos(&cl__1);
	io___37.ciunit = iounit_1.iout;
	s_wsfe(&io___37);
	do_fio(&c__1, indfile, trimtext_(indfile, (ftnlen)120));
	e_wsfe();
    }

/*     test for requested termination of the dynamics calculation */

    s_copy(endfile, "tinker.end", (ftnlen)120, (ftnlen)10);
    ioin__1.inerr = 0;
    ioin__1.infilen = 120;
    ioin__1.infile = endfile;
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
    if (! exist) {
/* Writing concatenation */
	i__4[0] = files_1.leng, a__3[0] = files_1.filename;
	i__4[1] = 4, a__3[1] = ".end";
	s_cat(endfile, a__3, i__4, &c__2, (ftnlen)120);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = endfile;
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
	    iend = freeunit_();
	    o__1.oerr = 0;
	    o__1.ounit = iend;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = endfile;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	    cl__1.cerr = 0;
	    cl__1.cunit = iend;
	    cl__1.csta = "delete";
	    f_clos(&cl__1);
	}
    }
    if (exist) {
	io___40.ciunit = iounit_1.iout;
	s_wsfe(&io___40);
	e_wsfe();
	fatal_();
    }

/*     skip an extra line to keep the output formating neat */

    moddump = *istep % inform_1.iprint;
    if (inform_1.verbose && moddump != 0) {
	io___41.ciunit = iounit_1.iout;
	s_wsfe(&io___41);
	e_wsfe();
    }
    return 0;
} /* mdsave_ */

#undef uind_ref
#undef name___ref
#undef wcm_ref
#undef vcm_ref
#undef v_ref
#undef a_ref


