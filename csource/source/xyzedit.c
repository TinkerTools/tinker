/* xyzedit.f -- translated by f2c (version 20050501).
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
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal vdwcut, chgcut, dplcut, mpolecut, vdwtaper, chgtaper, dpltaper,
	     mpoletaper, ewaldcut;
    logical use_ewald__, use_lights__, use_list__, use_vlist__, use_clist__, 
	    use_mlist__;
} cutoff_;

#define cutoff_1 cutoff_

struct {
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer ltitle;
    char title[120];
} titles_;

#define titles_1 titles_

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

struct {
    doublereal xref[250000]	/* was [25000][10] */, yref[250000]	/* 
	    was [25000][10] */, zref[250000]	/* was [25000][10] */;
    integer nref[10], reftyp[250000]	/* was [25000][10] */, n12ref[250000]	
	    /* was [25000][10] */, i12ref[2000000]	/* was [8][25000][10] 
	    */, refleng[10], refltitle[10];
    char refnam[750000]	/* was [25000][10] */, reffile[1200], reftitle[1200];
} refer_;

#define refer_1 refer_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  program xyzedit  --  editing of Cartesian coordinates  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "xyzedit" provides for modification and manipulation */
/*     of the contents of a Cartesian coordinates file */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 The XYZ Coordinate Editing Facility can "
	    "Provide :\002,//,4x,\002(1) Offset the Numbers of the Current At"
	    "oms\002,/,4x,\002(2) Deletion of Individual Specified Atoms\002,"
	    "/,4x,\002(3) Deletion of Specified Types of Atoms\002,/,4x,\002("
	    "4) Deletion of Atoms outside Cutoff Range\002,/,4x,\002(5) Inser"
	    "tion of Individual Specified Atoms\002,/,4x,\002(6) Replace Old "
	    "Atom Type with a New Type\002,/,4x,\002(7) Assign Connectivities"
	    " based on Distance\002,/,4x,\002(8) Convert Units from Bohrs to "
	    "Angstroms\002,/,4x,\002(9) Invert thru Origin to give Mirror Ima"
	    "ge\002,/,3x,\002(10) Translate All Atoms by an X,Y,Z-Vector\002,"
	    "/,3x,\002(11) Translate Center of Mass to the Origin\002,/,3x"
	    ",\002(12) Translate a Specified Atom to the Origin\002,/,3x,\002"
	    "(13) Translate and Rotate to Inertial Frame\002,/,3x,\002(14) Mo"
	    "ve to Specified Rigid Body Coordinates\002,/,3x,\002(15) Move St"
	    "ray Molecules into Periodic Box\002,/,3x,\002(16) Create and Fil"
	    "l a Periodic Boundary Box\002,/,3x,\002(17) Soak Current Molecul"
	    "e in Box of Solvent\002,/,3x,\002(18) Append a Second XYZ file t"
	    "o Current One\002)";
    static char fmt_30[] = "(/,\002 Number of the Desired Choice [<CR>=Exit]"
	    " :  \002,$)";
    static char fmt_40[] = "(i10)";
    static char fmt_70[] = "(/,\002 Offset used to Renumber the Current Atom"
	    "s :  \002,$)";
    static char fmt_80[] = "(i10)";
    static char fmt_90[] = "(/,\002 Numbers of the Atoms to be Removed : "
	    " \002,$)";
    static char fmt_100[] = "(a120)";
    static char fmt_120[] = "(/,\002 Atom Types to be Removed :  \002,$)";
    static char fmt_130[] = "(a120)";
    static char fmt_170[] = "(/,\002 Numbers of the Atoms to be Inserted : "
	    " \002,$)";
    static char fmt_180[] = "(a120)";
    static char fmt_210[] = "(/,\002 Numbers of the Old and New Atom Types :"
	    "  \002,$)";
    static char fmt_220[] = "(a120)";
    static char fmt_230[] = "(/,\002 Enter Translation Vector Components : "
	    " \002,$)";
    static char fmt_240[] = "(a120)";
    static char fmt_260[] = "(/,\002 Number of the Atom to Move to the Origi"
	    "n :  \002,$)";
    static char fmt_270[] = "(i10)";
    static char fmt_280[] = "(/,\002 Enter Rigid Body Coordinates :  \002,$)";
    static char fmt_290[] = "(a120)";
    static char fmt_310[] = "(/,\002 Enter Number of Molecules in Box :  "
	    "\002,$)";
    static char fmt_320[] = "(i10)";
    static char fmt_330[] = "(/,\002 Enter Periodic Box Dimensions (X,Y,Z) :"
	    "  \002,$)";
    static char fmt_340[] = "(a120)";
    static char fmt_360[] = "(i6)";
    static char fmt_370[] = "(i6,2x,a)";
    static char fmt_380[] = "(i6,2x,a3,3f12.6,5i6)";
    static char fmt_390[] = "(/,\002 New Coordinates written to :  \002,a)";
    static char fmt_400[] = "(i6)";
    static char fmt_410[] = "(i6,2x,a)";
    static char fmt_420[] = "(i6,2x,a3,3f12.6,5i6)";
    static char fmt_430[] = "(/,\002 New Coordinates written to File :  \002"
	    ",a)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3, i__4[2], i__5, i__6;
    doublereal d__1, d__2, d__3;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    double sqrt(doublereal), cos(doublereal), sin(doublereal);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int molecule_(void), unitcell_(void);
    extern integer freeunit_(void);
    static doublereal a[9]	/* was [3][3] */;
    static integer i__, j, k, nmolecule;
    static doublereal ri;
    static integer it;
    static doublereal xi, yi, zi, xr, yr, zr, dij, rad[25000], phi, rij, xcm, 
	    ycm, zcm, psi, cut2, cphi;
    static integer mode, keep[25000];
    extern /* Subroutine */ int soak_(void);
    static doublereal cpsi, sphi, xran, yran;
    static integer list[25000];
    static doublereal zran, norm, spsi, xbox, ybox, zbox;
    extern /* Subroutine */ int sort_(integer *, integer *);
    static integer ixyz;
    static doublereal dist2;
    extern /* Subroutine */ int sort4_(integer *, integer *), field_(void), 
	    final_(void), merge_(integer *);
    static integer nmode;
    static doublereal weigh, theta;
    extern /* Subroutine */ int katom_(void);
    static integer natom;
    static doublereal xorig;
    static integer nlist;
    static doublereal yorig, zorig;
    static logical write;
    extern /* Subroutine */ int delete_(integer *);
    static doublereal ctheta;
    static char record[120];
    extern doublereal random_(void);
    extern /* Subroutine */ int active_(void);
    static integer offset, atmnum, origin;
    static doublereal stheta;
    extern /* Subroutine */ int insert_(integer *), bounds_(void), getxyz_(
	    void), prtxyz_(integer *), makeref_(integer *), lattice_(void), 
	    initial_(void), inertia_(integer *), cutoffs_(void);
    static integer oldtype;
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen);
    static char xyzfile[120];
    static integer newtype;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___7 = { 1, 0, 1, fmt_40, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___9 = { 1, 0, 0, fmt_80, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_100, 0 };
    static icilist io___16 = { 1, record, 1, 0, 120, 1 };
    static cilist io___18 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_130, 0 };
    static icilist io___20 = { 1, record, 1, 0, 120, 1 };
    static cilist io___29 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_180, 0 };
    static icilist io___31 = { 1, record, 1, 0, 120, 1 };
    static cilist io___34 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_220, 0 };
    static icilist io___36 = { 1, record, 1, 0, 120, 1 };
    static cilist io___45 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___46 = { 0, 0, 0, fmt_240, 0 };
    static icilist io___47 = { 1, record, 1, 0, 120, 1 };
    static cilist io___53 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___54 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_290, 0 };
    static icilist io___64 = { 1, record, 1, 0, 120, 1 };
    static cilist io___72 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___73 = { 0, 0, 0, fmt_320, 0 };
    static cilist io___78 = { 0, 0, 0, fmt_330, 0 };
    static cilist io___79 = { 0, 0, 0, fmt_340, 0 };
    static icilist io___80 = { 1, record, 1, 0, 120, 1 };
    static cilist io___83 = { 0, 0, 0, fmt_360, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_370, 0 };
    static cilist io___89 = { 0, 0, 0, fmt_380, 0 };
    static cilist io___90 = { 0, 0, 0, fmt_390, 0 };
    static cilist io___91 = { 0, 0, 0, fmt_400, 0 };
    static cilist io___92 = { 0, 0, 0, fmt_410, 0 };
    static cilist io___93 = { 0, 0, 0, fmt_420, 0 };
    static cilist io___94 = { 0, 0, 0, fmt_430, 0 };



#define a_ref(a_1,a_2) a[(a_2)*3 + a_1 - 4]
#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]



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




/*     initialize various constants and the write flag */

    initial_();
    offset = 0;
    write = FALSE_;

/*     read in the coordinates and force field definition */

    getxyz_();
    field_();
    katom_();

/*     present a list of possible coordinate modifications */

    io___3.ciunit = iounit_1.iout;
    s_wsfe(&io___3);
    e_wsfe();

/*     get the desired type of coordinate file modification */

L20:
    nmode = 18;
    mode = -1;
    while(mode < 0 || mode > nmode) {
	mode = 0;
	io___6.ciunit = iounit_1.iout;
	s_wsfe(&io___6);
	e_wsfe();
	io___7.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___7);
	if (i__1 != 0) {
	    goto L100001;
	}
	i__1 = do_fio(&c__1, (char *)&mode, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L100001;
	}
	i__1 = e_rsfe();
L100001:
	if (i__1 < 0) {
	    goto L50;
	}
	if (i__1 > 0) {
	    goto L20;
	}
L50:
	;
    }

/*     get the offset value to be used in atom renumbering */

    if (mode == 1) {
L60:
	io___8.ciunit = iounit_1.iout;
	s_wsfe(&io___8);
	e_wsfe();
	io___9.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___9);
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_fio(&c__1, (char *)&offset, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L60;
	}
	write = TRUE_;
    }

/*     remove a specified list of individual atoms */

    if (mode == 2) {
	nlist = 0;
	for (i__ = 1; i__ <= 25000; ++i__) {
	    list[i__ - 1] = 0;
	}
	io___13.ciunit = iounit_1.iout;
	s_wsfe(&io___13);
	e_wsfe();
	io___14.ciunit = iounit_1.input;
	s_rsfe(&io___14);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___16);
	if (i__1 != 0) {
	    goto L110;
	}
	for (i__ = 1; i__ <= 25000; ++i__) {
	    i__1 = do_lio(&c__3, &c__1, (char *)&list[i__ - 1], (ftnlen)
		    sizeof(integer));
	    if (i__1 != 0) {
		goto L110;
	    }
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L110;
	}
L110:
	while(list[nlist] != 0) {
	    ++nlist;
	}
	i__1 = nlist;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (list[i__ - 1] > atoms_1.n) {
		list[i__ - 1] = atoms_1.n;
	    }
	    if (list[i__ - 1] < -atoms_1.n) {
		list[i__ - 1] = -atoms_1.n;
	    }
	}
	sort4_(&nlist, list);
	for (i__ = nlist; i__ >= 1; --i__) {
	    if (i__ > 1) {
		if (list[i__ - 2] < 0) {
		    i__3 = (i__2 = list[i__ - 2], abs(i__2));
		    for (j = (i__1 = list[i__ - 1], abs(i__1)); j >= i__3; 
			    --j) {
			delete_(&j);
		    }
		} else if (list[i__ - 1] > 0) {
		    delete_(&list[i__ - 1]);
		}
	    } else if (list[i__ - 1] > 0) {
		delete_(&list[i__ - 1]);
	    }
	}
	write = TRUE_;
	goto L20;
    }

/*     remove all atoms with any of a specified list of atom types */

    if (mode == 3) {
	nlist = 0;
	for (i__ = 1; i__ <= 25000; ++i__) {
	    list[i__ - 1] = 0;
	}
	io___18.ciunit = iounit_1.iout;
	s_wsfe(&io___18);
	e_wsfe();
	io___19.ciunit = iounit_1.input;
	s_rsfe(&io___19);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__3 = s_rsli(&io___20);
	if (i__3 != 0) {
	    goto L140;
	}
	for (i__ = 1; i__ <= 25000; ++i__) {
	    i__3 = do_lio(&c__3, &c__1, (char *)&list[i__ - 1], (ftnlen)
		    sizeof(integer));
	    if (i__3 != 0) {
		goto L140;
	    }
	}
	i__3 = e_rsli();
	if (i__3 != 0) {
	    goto L140;
	}
L140:
	while(list[nlist] != 0) {
	    ++nlist;
	}
	natom = atoms_1.n;
	for (i__ = natom; i__ >= 1; --i__) {
	    it = atoms_1.type__[i__ - 1];
	    i__3 = nlist;
	    for (j = 1; j <= i__3; ++j) {
		if (list[j - 1] == it) {
		    delete_(&i__);
		    goto L150;
		}
	    }
L150:
	    ;
	}
	write = TRUE_;
	goto L20;
    }

/*     remove all atoms that are inactive and lie outside all cutoffs */

    if (mode == 4) {
	active_();
	cutoffs_();
	nlist = 0;
	i__3 = atoms_1.n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    keep[i__ - 1] = 0;
	}
	cut2 = 0.;
	if (cutoff_1.vdwcut <= 1e3) {
/* Computing MAX */
/* Computing 2nd power */
	    d__2 = cutoff_1.vdwcut;
	    d__1 = d__2 * d__2;
	    cut2 = max(d__1,cut2);
	}
	if (cutoff_1.chgcut <= 1e3) {
/* Computing MAX */
/* Computing 2nd power */
	    d__2 = cutoff_1.chgcut;
	    d__1 = d__2 * d__2;
	    cut2 = max(d__1,cut2);
	}
	if (cutoff_1.dplcut <= 1e3) {
/* Computing MAX */
/* Computing 2nd power */
	    d__2 = cutoff_1.dplcut;
	    d__1 = d__2 * d__2;
	    cut2 = max(d__1,cut2);
	}
	if (cutoff_1.mpolecut <= 1e3) {
/* Computing MAX */
/* Computing 2nd power */
	    d__2 = cutoff_1.mpolecut;
	    d__1 = d__2 * d__2;
	    cut2 = max(d__1,cut2);
	}
	if (cut2 == 0.) {
	    cut2 = 1e16;
	}
	i__3 = atoms_1.n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    if (! usage_1.use[i__ - 1]) {
		i__1 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__1; ++j) {
		    keep[i12_ref(j, i__) - 1] = i__;
		}
		i__1 = couple_1.n13[i__ - 1];
		for (j = 1; j <= i__1; ++j) {
		    keep[i13_ref(j, i__) - 1] = i__;
		}
		i__1 = couple_1.n14[i__ - 1];
		for (j = 1; j <= i__1; ++j) {
		    keep[i14_ref(j, i__) - 1] = i__;
		}
		xi = atoms_1.x[i__ - 1];
		yi = atoms_1.y[i__ - 1];
		zi = atoms_1.z__[i__ - 1];
		i__1 = atoms_1.n;
		for (j = 1; j <= i__1; ++j) {
		    if (usage_1.use[j - 1]) {
			if (keep[j - 1] == i__) {
			    goto L160;
			}
/* Computing 2nd power */
			d__1 = atoms_1.x[j - 1] - xi;
/* Computing 2nd power */
			d__2 = atoms_1.y[j - 1] - yi;
/* Computing 2nd power */
			d__3 = atoms_1.z__[j - 1] - zi;
			dist2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
			if (dist2 <= cut2) {
			    goto L160;
			}
		    }
		}
		++nlist;
		list[nlist - 1] = i__;
L160:
		;
	    }
	}
	for (i__ = nlist; i__ >= 1; --i__) {
	    delete_(&list[i__ - 1]);
	}
	write = TRUE_;
	goto L20;
    }

/*     insert a specified list of individual atoms */

    if (mode == 5) {
	nlist = 0;
	for (i__ = 1; i__ <= 25000; ++i__) {
	    list[i__ - 1] = 0;
	}
	io___29.ciunit = iounit_1.iout;
	s_wsfe(&io___29);
	e_wsfe();
	io___30.ciunit = iounit_1.input;
	s_rsfe(&io___30);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__3 = s_rsli(&io___31);
	if (i__3 != 0) {
	    goto L190;
	}
	for (i__ = 1; i__ <= 25000; ++i__) {
	    i__3 = do_lio(&c__3, &c__1, (char *)&list[i__ - 1], (ftnlen)
		    sizeof(integer));
	    if (i__3 != 0) {
		goto L190;
	    }
	}
	i__3 = e_rsli();
	if (i__3 != 0) {
	    goto L190;
	}
L190:
	while(list[nlist] != 0) {
	    ++nlist;
	}
	sort4_(&nlist, list);
	for (i__ = nlist; i__ >= 1; --i__) {
	    if (i__ > 1) {
		if (list[i__ - 2] < 0) {
		    i__2 = (i__1 = list[i__ - 1], abs(i__1));
		    for (j = (i__3 = list[i__ - 2], abs(i__3)); j <= i__2; 
			    ++j) {
			insert_(&j);
		    }
		} else if (list[i__ - 1] > 0) {
		    insert_(&list[i__ - 1]);
		}
	    } else if (list[i__ - 1] > 0) {
		insert_(&list[i__ - 1]);
	    }
	}
	write = TRUE_;
	goto L20;
    }

/*     get an old atom type and new atom type for replacement */

    if (mode == 6) {
L200:
	oldtype = 0;
	newtype = 0;
	io___34.ciunit = iounit_1.iout;
	s_wsfe(&io___34);
	e_wsfe();
	io___35.ciunit = iounit_1.input;
	s_rsfe(&io___35);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__2 = s_rsli(&io___36);
	if (i__2 != 0) {
	    goto L200;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&oldtype, (ftnlen)sizeof(integer))
		;
	if (i__2 != 0) {
	    goto L200;
	}
	i__2 = do_lio(&c__3, &c__1, (char *)&newtype, (ftnlen)sizeof(integer))
		;
	if (i__2 != 0) {
	    goto L200;
	}
	i__2 = e_rsli();
	if (i__2 != 0) {
	    goto L200;
	}
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (atoms_1.type__[i__ - 1] == oldtype) {
		atoms_1.type__[i__ - 1] = newtype;
	    }
	}
	katom_();
	write = TRUE_;
	goto L20;
    }

/*     assign atom connectivities based on interatomic distances */

    if (mode == 7) {
	unitcell_();
	lattice_();
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    rad[i__ - 1] = .77;
	    atmnum = atmtyp_1.atomic[i__ - 1];
	    if (atmnum == 0) {
		rad[i__ - 1] = 0.;
	    }
	    if (atmnum == 1) {
		rad[i__ - 1] = .37;
	    }
	    if (atmnum == 2) {
		rad[i__ - 1] = .32;
	    }
	    if (atmnum == 6) {
		rad[i__ - 1] = .77;
	    }
	    if (atmnum == 7) {
		rad[i__ - 1] = .75;
	    }
	    if (atmnum == 8) {
		rad[i__ - 1] = .73;
	    }
	    if (atmnum == 9) {
		rad[i__ - 1] = .71;
	    }
	    if (atmnum == 10) {
		rad[i__ - 1] = .69;
	    }
	    if (atmnum == 14) {
		rad[i__ - 1] = 1.11;
	    }
	    if (atmnum == 15) {
		rad[i__ - 1] = 1.06;
	    }
	    if (atmnum == 16) {
		rad[i__ - 1] = 1.02;
	    }
	    if (atmnum == 17) {
		rad[i__ - 1] = .99;
	    }
	    if (atmnum == 18) {
		rad[i__ - 1] = .97;
	    }
	    if (atmnum == 35) {
		rad[i__ - 1] = 1.14;
	    }
	    if (atmnum == 36) {
		rad[i__ - 1] = 1.1;
	    }
	    if (atmnum == 53) {
		rad[i__ - 1] = 1.33;
	    }
	    if (atmnum == 54) {
		rad[i__ - 1] = 1.3;
	    }
	    rad[i__ - 1] *= 1.1;
	}
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    couple_1.n12[i__ - 1] = 0;
	}
	i__2 = atoms_1.n - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    ri = rad[i__ - 1];
	    i__3 = atoms_1.n;
	    for (j = i__ + 1; j <= i__3; ++j) {
		xr = atoms_1.x[j - 1] - xi;
		yr = atoms_1.y[j - 1] - yi;
		zr = atoms_1.z__[j - 1] - zi;
		rij = ri + rad[j - 1];
		dij = sqrt(xr * xr + yr * yr + zr * zr);
		if (dij < rij) {
		    ++couple_1.n12[i__ - 1];
		    i12_ref(couple_1.n12[i__ - 1], i__) = j;
		    ++couple_1.n12[j - 1];
		    i12_ref(couple_1.n12[j - 1], j) = i__;
		}
	    }
	}
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sort_(&couple_1.n12[i__ - 1], &i12_ref(1, i__));
	}
	write = TRUE_;
	goto L20;
    }

/*     convert the coordinate units from Bohrs to Angstroms */

    if (mode == 8) {
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    atoms_1.x[i__ - 1] *= .52917720859;
	    atoms_1.y[i__ - 1] *= .52917720859;
	    atoms_1.z__[i__ - 1] *= .52917720859;
	}
	write = TRUE_;
	goto L20;
    }

/*     get mirror image by inverting coordinates through origin */

    if (mode == 9) {
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    atoms_1.x[i__ - 1] = -atoms_1.x[i__ - 1];
	    atoms_1.y[i__ - 1] = -atoms_1.y[i__ - 1];
	    atoms_1.z__[i__ - 1] = -atoms_1.z__[i__ - 1];
	}
	write = TRUE_;
	goto L20;
    }

/*     translate the entire system by a specified x,y,z-vector */

    if (mode == 10) {
	xr = 0.;
	yr = 0.;
	zr = 0.;
	io___45.ciunit = iounit_1.iout;
	s_wsfe(&io___45);
	e_wsfe();
	io___46.ciunit = iounit_1.input;
	s_rsfe(&io___46);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__2 = s_rsli(&io___47);
	if (i__2 != 0) {
	    goto L250;
	}
	i__2 = do_lio(&c__5, &c__1, (char *)&xr, (ftnlen)sizeof(doublereal));
	if (i__2 != 0) {
	    goto L250;
	}
	i__2 = do_lio(&c__5, &c__1, (char *)&yr, (ftnlen)sizeof(doublereal));
	if (i__2 != 0) {
	    goto L250;
	}
	i__2 = do_lio(&c__5, &c__1, (char *)&zr, (ftnlen)sizeof(doublereal));
	if (i__2 != 0) {
	    goto L250;
	}
	i__2 = e_rsli();
	if (i__2 != 0) {
	    goto L250;
	}
L250:
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    atoms_1.x[i__ - 1] += xr;
	    atoms_1.y[i__ - 1] += yr;
	    atoms_1.z__[i__ - 1] += zr;
	}
	write = TRUE_;
	goto L20;
    }

/*     translate the center of mass to the coordinate origin */

    if (mode == 11) {
	xcm = 0.;
	ycm = 0.;
	zcm = 0.;
	norm = 0.;
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    weigh = atmtyp_1.mass[i__ - 1];
	    xcm += atoms_1.x[i__ - 1] * weigh;
	    ycm += atoms_1.y[i__ - 1] * weigh;
	    zcm += atoms_1.z__[i__ - 1] * weigh;
	    norm += weigh;
	}
	xcm /= norm;
	ycm /= norm;
	zcm /= norm;
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    atoms_1.x[i__ - 1] -= xcm;
	    atoms_1.y[i__ - 1] -= ycm;
	    atoms_1.z__[i__ - 1] -= zcm;
	}
	write = TRUE_;
	goto L20;
    }

/*     translate to place a specified atom at the origin */

    if (mode == 12) {
	io___53.ciunit = iounit_1.iout;
	s_wsfe(&io___53);
	e_wsfe();
	io___54.ciunit = iounit_1.input;
	s_rsfe(&io___54);
	do_fio(&c__1, (char *)&origin, (ftnlen)sizeof(integer));
	e_rsfe();
	xorig = atoms_1.x[origin - 1];
	yorig = atoms_1.y[origin - 1];
	zorig = atoms_1.z__[origin - 1];
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    atoms_1.x[i__ - 1] -= xorig;
	    atoms_1.y[i__ - 1] -= yorig;
	    atoms_1.z__[i__ - 1] -= zorig;
	}
	write = TRUE_;
	goto L20;
    }

/*     translate and rotate into standard orientation */

    if (mode == 13) {
	inertia_(&c__2);
	write = TRUE_;
	goto L20;
    }

/*     translate and rotate to specified rigid body coordinates */

    if (mode == 14) {
	xcm = 0.;
	ycm = 0.;
	zcm = 0.;
	phi = 0.;
	theta = 0.;
	psi = 0.;
	io___62.ciunit = iounit_1.iout;
	s_wsfe(&io___62);
	e_wsfe();
	io___63.ciunit = iounit_1.input;
	s_rsfe(&io___63);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__2 = s_rsli(&io___64);
	if (i__2 != 0) {
	    goto L300;
	}
	i__2 = do_lio(&c__5, &c__1, (char *)&xcm, (ftnlen)sizeof(doublereal));
	if (i__2 != 0) {
	    goto L300;
	}
	i__2 = do_lio(&c__5, &c__1, (char *)&ycm, (ftnlen)sizeof(doublereal));
	if (i__2 != 0) {
	    goto L300;
	}
	i__2 = do_lio(&c__5, &c__1, (char *)&zcm, (ftnlen)sizeof(doublereal));
	if (i__2 != 0) {
	    goto L300;
	}
	i__2 = do_lio(&c__5, &c__1, (char *)&phi, (ftnlen)sizeof(doublereal));
	if (i__2 != 0) {
	    goto L300;
	}
	i__2 = do_lio(&c__5, &c__1, (char *)&theta, (ftnlen)sizeof(doublereal)
		);
	if (i__2 != 0) {
	    goto L300;
	}
	i__2 = do_lio(&c__5, &c__1, (char *)&psi, (ftnlen)sizeof(doublereal));
	if (i__2 != 0) {
	    goto L300;
	}
	i__2 = e_rsli();
	if (i__2 != 0) {
	    goto L300;
	}
L300:
	inertia_(&c__2);
	phi /= 57.29577951308232088;
	theta /= 57.29577951308232088;
	psi /= 57.29577951308232088;
	cphi = cos(phi);
	sphi = sin(phi);
	ctheta = cos(theta);
	stheta = sin(theta);
	cpsi = cos(psi);
	spsi = sin(psi);
	a_ref(1, 1) = ctheta * cphi;
	a_ref(2, 1) = spsi * stheta * cphi - cpsi * sphi;
	a_ref(3, 1) = cpsi * stheta * cphi + spsi * sphi;
	a_ref(1, 2) = ctheta * sphi;
	a_ref(2, 2) = spsi * stheta * sphi + cpsi * cphi;
	a_ref(3, 2) = cpsi * stheta * sphi - spsi * cphi;
	a_ref(1, 3) = -stheta;
	a_ref(2, 3) = ctheta * spsi;
	a_ref(3, 3) = ctheta * cpsi;
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xorig = atoms_1.x[i__ - 1];
	    yorig = atoms_1.y[i__ - 1];
	    zorig = atoms_1.z__[i__ - 1];
	    atoms_1.x[i__ - 1] = a_ref(1, 1) * xorig + a_ref(2, 1) * yorig + 
		    a_ref(3, 1) * zorig + xcm;
	    atoms_1.y[i__ - 1] = a_ref(1, 2) * xorig + a_ref(2, 2) * yorig + 
		    a_ref(3, 2) * zorig + ycm;
	    atoms_1.z__[i__ - 1] = a_ref(1, 3) * xorig + a_ref(2, 3) * yorig 
		    + a_ref(3, 3) * zorig + zcm;
	}
	write = TRUE_;
	goto L20;
    }

/*     move stray molecules back into original periodic box */

    if (mode == 15) {
	unitcell_();
	if (bound_1.use_bounds__) {
	    lattice_();
	    molecule_();
	    bounds_();
	    write = TRUE_;
	}
	goto L20;
    }

/*     create a random box full of the current coordinates file */

    if (mode == 16) {
	io___72.ciunit = iounit_1.iout;
	s_wsfe(&io___72);
	e_wsfe();
	io___73.ciunit = iounit_1.input;
	s_rsfe(&io___73);
	do_fio(&c__1, (char *)&nmolecule, (ftnlen)sizeof(integer));
	e_rsfe();
	xbox = 0.;
	ybox = 0.;
	zbox = 0.;
	while(xbox == 0.) {
	    io___78.ciunit = iounit_1.iout;
	    s_wsfe(&io___78);
	    e_wsfe();
	    io___79.ciunit = iounit_1.input;
	    s_rsfe(&io___79);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    i__2 = s_rsli(&io___80);
	    if (i__2 != 0) {
		goto L350;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&xbox, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L350;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&ybox, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L350;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&zbox, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L350;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L350;
	    }
L350:
	    if (ybox == 0.) {
		ybox = xbox;
	    }
	    if (zbox == 0.) {
		zbox = xbox;
	    }
	}
	ixyz = freeunit_();
/* Writing concatenation */
	i__4[0] = files_1.leng, a__1[0] = files_1.filename;
	i__4[1] = 4, a__1[1] = ".xyz";
	s_cat(xyzfile, a__1, i__4, &c__2, (ftnlen)120);
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
	if (titles_1.ltitle == 0) {
	    io___83.ciunit = ixyz;
	    s_wsfe(&io___83);
	    i__2 = atoms_1.n * nmolecule;
	    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___84.ciunit = ixyz;
	    s_wsfe(&io___84);
	    i__2 = atoms_1.n * nmolecule;
	    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, titles_1.title, titles_1.ltitle);
	    e_wsfe();
	}
	i__2 = nmolecule;
	for (k = 1; k <= i__2; ++k) {
	    offset = (k - 1) * atoms_1.n;
	    xran = xbox * random_();
	    yran = ybox * random_();
	    zran = zbox * random_();
	    i__3 = atoms_1.n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		io___89.ciunit = ixyz;
		s_wsfe(&io___89);
		i__1 = i__ + offset;
		do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
		d__1 = atoms_1.x[i__ - 1] + xran;
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		d__2 = atoms_1.y[i__ - 1] + yran;
		do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
		d__3 = atoms_1.z__[i__ - 1] + zran;
		do_fio(&c__1, (char *)&d__3, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&atoms_1.type__[i__ - 1], (ftnlen)
			sizeof(integer));
		i__5 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__5; ++j) {
		    i__6 = i12_ref(j, i__) + offset;
		    do_fio(&c__1, (char *)&i__6, (ftnlen)sizeof(integer));
		}
		e_wsfe();
	    }
	}
	cl__1.cerr = 0;
	cl__1.cunit = ixyz;
	cl__1.csta = 0;
	f_clos(&cl__1);
	io___90.ciunit = iounit_1.iout;
	s_wsfe(&io___90);
	do_fio(&c__1, xyzfile, (ftnlen)120);
	e_wsfe();
	write = FALSE_;
    }

/*     solvate the current system by insertion into a solvent box */

    if (mode == 17) {
	soak_();
	write = TRUE_;
    }

/*     append a second file to the current coordinates file */

    if (mode == 18) {
	makeref_(&c__1);
	getxyz_();
	merge_(&c__1);
	write = TRUE_;
	goto L20;
    }

/*     write out the new coordinates file in its new format */

    if (write) {
	ixyz = freeunit_();
/* Writing concatenation */
	i__4[0] = files_1.leng, a__1[0] = files_1.filename;
	i__4[1] = 4, a__1[1] = ".xyz";
	s_cat(xyzfile, a__1, i__4, &c__2, (ftnlen)120);
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
	if (offset == 0) {
	    prtxyz_(&ixyz);
	} else {
	    if (titles_1.ltitle == 0) {
		io___91.ciunit = ixyz;
		s_wsfe(&io___91);
		do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
		e_wsfe();
	    } else {
		io___92.ciunit = ixyz;
		s_wsfe(&io___92);
		do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
		do_fio(&c__1, titles_1.title, titles_1.ltitle);
		e_wsfe();
	    }
	    i__2 = atoms_1.n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		io___93.ciunit = ixyz;
		s_wsfe(&io___93);
		i__3 = i__ + offset;
		do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
		do_fio(&c__1, (char *)&atoms_1.x[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&atoms_1.y[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&atoms_1.z__[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&atoms_1.type__[i__ - 1], (ftnlen)
			sizeof(integer));
		i__1 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__1; ++j) {
		    i__6 = i12_ref(j, i__) + offset;
		    do_fio(&c__1, (char *)&i__6, (ftnlen)sizeof(integer));
		}
		e_wsfe();
	    }
	}
	cl__1.cerr = 0;
	cl__1.cunit = ixyz;
	cl__1.csta = 0;
	f_clos(&cl__1);
	io___94.ciunit = iounit_1.iout;
	s_wsfe(&io___94);
	do_fio(&c__1, xyzfile, (ftnlen)120);
	e_wsfe();
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef name___ref
#undef i14_ref
#undef i13_ref
#undef i12_ref
#undef a_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine soak  --  place a solute into a solvent box  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "soak" takes a currently defined solute system and places */
/*     it into a solvent box, with removal of any solvent molecules */
/*     that overlap the solute */


/* Subroutine */ int soak_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Enter Name of Solvent Box Coordinates "
	    ":  \002,$)";
    static char fmt_30[] = "(a120)";

    /* System generated locals */
    integer i__1, i__2;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void), f_open(olist *), f_rew(alist *), 
	    f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int molecule_(void), unitcell_(void);
    extern integer freeunit_(void);
    static char solvfile[120];
    static integer i__, k;
    static doublereal xi, yi, zi, xr, yr, zr, rik2;
    static integer ntot;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *), merge_(integer *);
    static doublereal close;
    static integer isolv;
    static doublereal close2;
    extern /* Subroutine */ int delete_(integer *);
    static logical remove[25000];
    extern /* Subroutine */ int suffix_(char *, char *, ftnlen, ftnlen), 
	    makeref_(integer *), lattice_(void), version_(char *, char *, 
	    ftnlen, ftnlen), readxyz_(integer *);

    /* Fortran I/O blocks */
    static cilist io___95 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___96 = { 0, 0, 0, fmt_30, 0 };




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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  refer.i  --  storage of reference atomic coordinate set  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xref        reference x-coordinates for atoms in each system */
/*     yref        reference y-coordinates for atoms in each system */
/*     zref        reference z-coordinates for atoms in each system */
/*     nref        total number of atoms in each reference system */
/*     reftyp      atom types of the atoms in each reference system */
/*     n12ref      number of atoms bonded to each reference atom */
/*     i12ref      atom numbers of atoms 1-2 connected to each atom */
/*     refleng     length in characters of each reference filename */
/*     refltitle   length in characters of each reference title line */
/*     refnam      atom names of the atoms in each reference system */
/*     reffile     base filename for each reference system */
/*     reftitle    title used to describe each reference system */




/*     make a copy of the solute coordinates and connectivities */

    makeref_(&c__1);

/*     read the coordinates for the solvent box */

L10:
    io___95.ciunit = iounit_1.iout;
    s_wsfe(&io___95);
    e_wsfe();
    io___96.ciunit = iounit_1.input;
    s_rsfe(&io___96);
    do_fio(&c__1, solvfile, (ftnlen)120);
    e_rsfe();
    suffix_(solvfile, "xyz", (ftnlen)120, (ftnlen)3);
    version_(solvfile, "old", (ftnlen)120, (ftnlen)3);
    isolv = freeunit_();
    o__1.oerr = 1;
    o__1.ounit = isolv;
    o__1.ofnmlen = 120;
    o__1.ofnm = solvfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L10;
    }
    al__1.aerr = 0;
    al__1.aunit = isolv;
    f_rew(&al__1);
    readxyz_(&isolv);
    cl__1.cerr = 0;
    cl__1.cunit = isolv;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     combine solute and solvent into a single coordinate set */

    merge_(&c__1);

/*     count number of molecules and set lattice parameters */

    molecule_();
    unitcell_();
    lattice_();

/*     initialize the list of solvent molecules to be deleted */

    i__1 = molcul_1.nmol;
    for (i__ = 1; i__ <= i__1; ++i__) {
	remove[i__ - 1] = FALSE_;
    }

/*     search for close contacts between solute and solvent */

    close = 1.5;
    close2 = close * close;
    i__1 = refer_1.nref[0];
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	i__2 = atoms_1.n;
	for (k = refer_1.nref[0] + 1; k <= i__2; ++k) {
	    xr = atoms_1.x[k - 1] - xi;
	    yr = atoms_1.y[k - 1] - yi;
	    zr = atoms_1.z__[k - 1] - zi;
	    image_(&xr, &yr, &zr);
	    rik2 = xr * xr + yr * yr + zr * zr;
	    if (rik2 < close2) {
		remove[molcul_1.molcule[k - 1] - 1] = TRUE_;
	    }
	}
    }

/*     remove solvent molecules that are too close to the solute */

    ntot = atoms_1.n;
    i__1 = refer_1.nref[0] + 1;
    for (i__ = ntot; i__ >= i__1; --i__) {
	if (remove[molcul_1.molcule[i__ - 1] - 1]) {
	    delete_(&i__);
	}
    }
    return 0;
} /* soak_ */

/* Main program alias */ int xyzedit_ () { MAIN__ (); return 0; }
