/* poledit.f -- translated by f2c (version 20050501).
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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

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
    doublereal mp[25000], dpx[25000], dpy[25000], dpz[25000], q20[25000], 
	    q21c[25000], q21s[25000], q22c[25000], q22s[25000];
} dma_;

#define dma_1 dma_

struct {
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

struct {
    doublereal polarity[25000], thole[25000], pdamp[25000], uind[75000]	/* 
	    was [3][25000] */, uinp[75000]	/* was [3][25000] */, uinds[
	    75000]	/* was [3][25000] */, uinps[75000]	/* was [3][
	    25000] */;
    integer npolar;
} polar_;

#define polar_1 polar_

struct {
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

struct {
    doublereal polr[5000], athl[5000];
    integer pgrp[40000]	/* was [8][5000] */;
} kpolr_;

#define kpolr_1 kpolr_

struct {
    integer np11[25000], ip11[2500000]	/* was [100][25000] */, np12[25000], 
	    ip12[1250000]	/* was [50][25000] */, np13[25000], ip13[
	    1250000]	/* was [50][25000] */, np14[25000], ip14[1250000]	
	    /* was [50][25000] */;
} polgrp_;

#define polgrp_1 polgrp_

struct {
    doublereal poleps, polsor, p2scale, p3scale, p4scale, p5scale, d1scale, 
	    d2scale, d3scale, d4scale, u1scale, u2scale, u3scale, u4scale;
    char poltyp[6];
} polpot_;

#define polpot_1 polpot_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__2 = 2;



/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2000 by P. Bagossi, P. Ren & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  program poledit  --  manipulate atomic multipole values  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "poledit" provides for modification and manipulation of */
/*     the atomic multipole electrostatic models used in TINKER */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 The Multipole Editing Facility can Provi"
	    "de :\002,//,4x,\002(1) Multipole Parameters from GDMA Output\002"
	    ",/,4x,\002(2) Alter Local Coordinate Frame Definitions\002,/,4x"
	    ",\002(3) Removal of Intramolecular Polarization\002)";
    static char fmt_30[] = "(/,\002 Enter the Number of the Desired Choice :"
	    "  \002,$)";
    static char fmt_40[] = "(i10)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);

    /* Local variables */
    extern /* Subroutine */ int readgdma_(void), fixframe_(void), setframe_(
	    void), rotframe_(void), fixpolar_(void), intrapol_(void), 
	    setpolar_(void), prtpolar_(void), molsetup_(void);
    static integer mode;
    extern /* Subroutine */ int field_(void), katom_(void);
    static logical exist, query;
    extern /* Subroutine */ int attach_(void), kmpole_(void), kpolar_(void);
    static char string[120];
    extern /* Subroutine */ int getxyz_(void), initial_(void), nextarg_(char *
	    , logical *, ftnlen), initprm_(void);

    /* Fortran I/O blocks */
    static icilist io___5 = { 1, string, 1, 0, 120, 1 };
    static cilist io___6 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___8 = { 1, 0, 1, fmt_40, 0 };




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




/*     get the desired type of coordinate file modification */

    initial_();
    mode = 0;
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___5);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&mode, (ftnlen)sizeof(integer));
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
	while(mode < 1 || mode > 3) {
	    mode = 0;
	    io___7.ciunit = iounit_1.iout;
	    s_wsfe(&io___7);
	    e_wsfe();
	    io___8.ciunit = iounit_1.input;
	    i__1 = s_rsfe(&io___8);
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_fio(&c__1, (char *)&mode, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L50;
	    }
L50:
	    ;
	}
    }

/*     perform the desired multipole manipulation operation */

    if (mode == 1) {
	potent_1.use_mpole__ = TRUE_;
	potent_1.use_polar__ = TRUE_;
	readgdma_();
	initprm_();
	molsetup_();
	setframe_();
	rotframe_();
	setpolar_();
	intrapol_();
	fixpolar_();
	prtpolar_();
    } else if (mode == 2) {
	getxyz_();
	attach_();
	field_();
	katom_();
	kmpole_();
	kpolar_();
	fixframe_();
	prtpolar_();
    } else if (mode == 3) {
	getxyz_();
	attach_();
	field_();
	katom_();
	kmpole_();
	kpolar_();
	intrapol_();
	fixpolar_();
	prtpolar_();
    }
    return 0;
} /* MAIN__ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine readgdma  --  get information from GDMA output  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "readgdma" takes the DMA output in spherical harmonics from */
/*     the GDMA program and converts to Cartesian multipoles in */
/*     the global coordinate frame */

/*     this version is compatible with the formatted output from */
/*     GDMA 2.2.04 released by Anthony Stone in Fall 2008 */


/* Subroutine */ int readgdma_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter GDMA Output File Name :  \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(a120)";
    static char fmt_40[] = "()";
    static char fmt_60[] = "(a120)";
    static char fmt_80[] = "(a120)";
    static char fmt_90[] = "(a120)";
    static char fmt_120[] = "(/,\002 Global Frame Cartesian Multipole Moment"
	    "s :\002)";
    static char fmt_130[] = "(/,\002 Site:\002,i8,9x,\002Name:\002,3x,a3,7x"
	    ",\002Atomic Number:\002,i8)";
    static char fmt_140[] = "(/,\002 Coordinates:\002,5x,3f15.6)";
    static char fmt_150[] = "(/,\002 Charge:\002,10x,f15.5)";
    static char fmt_160[] = "(\002 Dipole:\002,10x,3f15.5)";
    static char fmt_170[] = "(\002 Quadrupole:\002,6x,f15.5)";
    static char fmt_180[] = "(18x,2f15.5)";
    static char fmt_190[] = "(18x,3f15.5)";

    /* System generated locals */
    integer i__1, i__2;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *)
	    , do_fio(integer *, char *, ftnlen), e_rsfe(void), f_open(olist *)
	    , f_rew(alist *), s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(
	    icilist *), do_lio(integer *, integer *, char *, ftnlen), e_rsli(
	    void);
    double sqrt(doublereal);
    integer f_clos(cllist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int basefile_(char *, ftnlen);
    static logical use_bohr__;
    extern integer freeunit_(void);
    static integer i__, j, k, idma;
    static logical done;
    static doublereal term;
    static integer next;
    static logical exist;
    extern /* Subroutine */ int match1_(integer *, char *, ftnlen);
    static char atmnam[3], record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), suffix_(char *, char 
	    *, ftnlen, ftnlen);
    static char dmafile[120];
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), gettext_(
	    char *, char *, integer *, ftnlen, ftnlen), version_(char *, char 
	    *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___12 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___15 = { 1, 0, 1, fmt_30, 0 };
    static icilist io___18 = { 0, record+14, 0, 0, 10, 1 };
    static icilist io___19 = { 0, record+29, 0, 0, 10, 1 };
    static icilist io___20 = { 0, record+44, 0, 0, 10, 1 };
    static cilist io___21 = { 1, 0, 1, fmt_40, 0 };
    static cilist io___24 = { 1, 0, 1, fmt_60, 0 };
    static cilist io___26 = { 1, 0, 1, fmt_80, 0 };
    static icilist io___28 = { 1, record+16, 1, 0, 104, 1 };
    static cilist io___30 = { 1, 0, 1, fmt_90, 0 };
    static icilist io___31 = { 1, record, 1, 0, 120, 1 };
    static icilist io___33 = { 1, atmnam, 1, 0, 3, 1 };
    static cilist io___34 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_190, 0 };



#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]



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
/*     ##  COPYRIGHT (C)  2005  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  dma.i  --  distributed multipole analysis components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     mp        atomic monopole charge values from DMA */
/*     dpx       atomic dipole moment x-component from DMA */
/*     dpy       atomic dipole moment y-component from DMA */
/*     dpz       atomic dipole moment z-component from DMA */
/*     q20       atomic Q20 quadrupole component from DMA (zz) */
/*     q21c      atomic Q21c quadrupole component from DMA (xz) */
/*     q21s      atomic Q21s quadrupole component from DMA (yz) */
/*     q22c      atomic Q22c quadrupole component from DMA (xx-yy) */
/*     q22s      atomic Q22s quadrupole component from DMA (xy) */




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




/*     zero out the atomic coordinates and DMA values */

    for (i__ = 1; i__ <= 25000; ++i__) {
	atoms_1.x[i__ - 1] = 0.;
	atoms_1.y[i__ - 1] = 0.;
	atoms_1.z__[i__ - 1] = 0.;
	dma_1.mp[i__ - 1] = 0.;
	dma_1.dpx[i__ - 1] = 0.;
	dma_1.dpy[i__ - 1] = 0.;
	dma_1.dpz[i__ - 1] = 0.;
	dma_1.q20[i__ - 1] = 0.;
	dma_1.q21c[i__ - 1] = 0.;
	dma_1.q21s[i__ - 1] = 0.;
	dma_1.q22c[i__ - 1] = 0.;
	dma_1.q22s[i__ - 1] = 0.;
    }

/*     try to get a filename from the command line arguments */

    nextarg_(dmafile, &exist, (ftnlen)120);
    if (exist) {
	basefile_(dmafile, (ftnlen)120);
	suffix_(dmafile, "dma", (ftnlen)120, (ftnlen)3);
	version_(dmafile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = dmafile;
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

/*     ask for the user specified GDMA output filename */

    while(! exist) {
	io___12.ciunit = iounit_1.iout;
	s_wsfe(&io___12);
	e_wsfe();
	io___13.ciunit = iounit_1.input;
	s_rsfe(&io___13);
	do_fio(&c__1, dmafile, (ftnlen)120);
	e_rsfe();
	basefile_(dmafile, (ftnlen)120);
	suffix_(dmafile, "dma", (ftnlen)120, (ftnlen)3);
	version_(dmafile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = dmafile;
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

/*     first open and then read the GDMA output file */

    idma = freeunit_();
    o__1.oerr = 0;
    o__1.ounit = idma;
    o__1.ofnmlen = 120;
    o__1.ofnm = dmafile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);

/*     get coordinates and multipoles from GDMA output file */

    i__ = 0;
    al__1.aerr = 0;
    al__1.aunit = idma;
    f_rew(&al__1);
    while(TRUE_) {
	io___15.ciunit = idma;
	i__1 = s_rsfe(&io___15);
	if (i__1 != 0) {
	    goto L50;
	}
	i__1 = do_fio(&c__1, record, (ftnlen)120);
	if (i__1 != 0) {
	    goto L50;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L50;
	}
	if (i__ != 0) {
	    match1_(&i__, record, (ftnlen)120);
	}
	if (s_cmp(record + 11, "x =", (ftnlen)3, (ftnlen)3) == 0) {
	    ++i__;
	    next = 1;
	    gettext_(record, name___ref(0, i__), &next, (ftnlen)120, (ftnlen)
		    3);
	    s_rsli(&io___18);
	    do_lio(&c__5, &c__1, (char *)&atoms_1.x[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_rsli();
	    s_rsli(&io___19);
	    do_lio(&c__5, &c__1, (char *)&atoms_1.y[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_rsli();
	    s_rsli(&io___20);
	    do_lio(&c__5, &c__1, (char *)&atoms_1.z__[i__ - 1], (ftnlen)
		    sizeof(doublereal));
	    e_rsli();
	    io___21.ciunit = idma;
	    i__1 = s_rsfe(&io___21);
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L50;
	    }
	} else if (s_cmp(record, "Total multipoles", (ftnlen)16, (ftnlen)16) 
		== 0) {
	    goto L50;
	}
    }
L50:
    atoms_1.n = i__;

/*     convert quadrupole from spherical harmonic to Cartesian */

    term = sqrt(.75);
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rpole_ref(1, i__) = dma_1.mp[i__ - 1];
	rpole_ref(2, i__) = dma_1.dpx[i__ - 1];
	rpole_ref(3, i__) = dma_1.dpy[i__ - 1];
	rpole_ref(4, i__) = dma_1.dpz[i__ - 1];
	rpole_ref(5, i__) = dma_1.q20[i__ - 1] * -.5 + term * dma_1.q22c[i__ 
		- 1];
	rpole_ref(6, i__) = term * dma_1.q22s[i__ - 1];
	rpole_ref(7, i__) = term * dma_1.q21c[i__ - 1];
	rpole_ref(8, i__) = rpole_ref(6, i__);
	rpole_ref(9, i__) = dma_1.q20[i__ - 1] * -.5 - term * dma_1.q22c[i__ 
		- 1];
	rpole_ref(10, i__) = term * dma_1.q21s[i__ - 1];
	rpole_ref(11, i__) = rpole_ref(7, i__);
	rpole_ref(12, i__) = rpole_ref(10, i__);
	rpole_ref(13, i__) = dma_1.q20[i__ - 1];
    }

/*     check for GDMA coordinate values in atomic units */

    use_bohr__ = FALSE_;
    al__1.aerr = 0;
    al__1.aunit = idma;
    f_rew(&al__1);
    while(TRUE_) {
	io___24.ciunit = idma;
	i__1 = s_rsfe(&io___24);
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = do_fio(&c__1, record, (ftnlen)120);
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L70;
	}
	if (s_cmp(record, "Positions and radii in bohr", (ftnlen)27, (ftnlen)
		27) == 0) {
	    use_bohr__ = TRUE_;
	    goto L70;
	}
    }
L70:

/*     convert coordinates from Bohrs to Angstroms if needed */

    if (use_bohr__) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    atoms_1.x[i__ - 1] *= .52917720859;
	    atoms_1.y[i__ - 1] *= .52917720859;
	    atoms_1.z__[i__ - 1] *= .52917720859;
	}
    }

/*     find atomic numbers in verbose GDMA output if available */

    done = FALSE_;
    al__1.aerr = 0;
    al__1.aunit = idma;
    f_rew(&al__1);
    while(TRUE_) {
	io___26.ciunit = idma;
	i__1 = s_rsfe(&io___26);
	if (i__1 != 0) {
	    goto L100;
	}
	i__1 = do_fio(&c__1, record, (ftnlen)120);
	if (i__1 != 0) {
	    goto L100;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L100;
	}
	if (s_cmp(record, "Nuclear charges:", (ftnlen)16, (ftnlen)16) == 0) {
	    k = min(atoms_1.n,20);
	    i__1 = s_rsli(&io___28);
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__2 = k;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__1 = do_lio(&c__3, &c__1, (char *)&atmtyp_1.atomic[i__ - 1],
			 (ftnlen)sizeof(integer));
		if (i__1 != 0) {
		    goto L100;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L100;
	    }
	    while(k != atoms_1.n) {
		j = k + 1;
/* Computing MIN */
		i__1 = atoms_1.n, i__2 = k + 20;
		k = min(i__1,i__2);
		io___30.ciunit = idma;
		i__1 = s_rsfe(&io___30);
		if (i__1 != 0) {
		    goto L100;
		}
		i__1 = do_fio(&c__1, record, (ftnlen)120);
		if (i__1 != 0) {
		    goto L100;
		}
		i__1 = e_rsfe();
		if (i__1 != 0) {
		    goto L100;
		}
		i__1 = s_rsli(&io___31);
		if (i__1 != 0) {
		    goto L100;
		}
		i__2 = k;
		for (i__ = j; i__ <= i__2; ++i__) {
		    i__1 = do_lio(&c__3, &c__1, (char *)&atmtyp_1.atomic[i__ 
			    - 1], (ftnlen)sizeof(integer));
		    if (i__1 != 0) {
			goto L100;
		    }
		}
		i__1 = e_rsli();
		if (i__1 != 0) {
		    goto L100;
		}
	    }
	    done = TRUE_;
	}
    }
L100:
    cl__1.cerr = 0;
    cl__1.cunit = idma;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     attempt to get atomic numbers from GDMA atom names */

    if (! done) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    atmtyp_1.atomic[i__ - 1] = 0;
	    s_copy(atmnam, name___ref(0, i__), (ftnlen)3, (ftnlen)3);
	    upcase_(atmnam, (ftnlen)3);
	    if (s_cmp(atmnam, "SI", (ftnlen)2, (ftnlen)2) == 0) {
		atmtyp_1.atomic[i__ - 1] = 14;
	    } else if (s_cmp(atmnam, "CL", (ftnlen)2, (ftnlen)2) == 0) {
		atmtyp_1.atomic[i__ - 1] = 17;
	    } else if (s_cmp(atmnam, "BR", (ftnlen)2, (ftnlen)2) == 0) {
		atmtyp_1.atomic[i__ - 1] = 35;
	    } else if (*(unsigned char *)atmnam == 'H') {
		atmtyp_1.atomic[i__ - 1] = 1;
	    } else if (*(unsigned char *)atmnam == 'B') {
		atmtyp_1.atomic[i__ - 1] = 5;
	    } else if (*(unsigned char *)atmnam == 'C') {
		atmtyp_1.atomic[i__ - 1] = 6;
	    } else if (*(unsigned char *)atmnam == 'N') {
		atmtyp_1.atomic[i__ - 1] = 7;
	    } else if (*(unsigned char *)atmnam == 'O') {
		atmtyp_1.atomic[i__ - 1] = 8;
	    } else if (*(unsigned char *)atmnam == 'F') {
		atmtyp_1.atomic[i__ - 1] = 9;
	    } else if (*(unsigned char *)atmnam == 'P') {
		atmtyp_1.atomic[i__ - 1] = 15;
	    } else if (*(unsigned char *)atmnam == 'S') {
		atmtyp_1.atomic[i__ - 1] = 16;
	    } else if (*(unsigned char *)atmnam == 'I') {
		atmtyp_1.atomic[i__ - 1] = 53;
	    } else {
		i__2 = s_rsli(&io___33);
		if (i__2 != 0) {
		    goto L110;
		}
		i__2 = do_lio(&c__3, &c__1, (char *)&atmtyp_1.atomic[i__ - 1],
			 (ftnlen)sizeof(integer));
		if (i__2 != 0) {
		    goto L110;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L110;
		}
L110:
		;
	    }
	}
    }

/*     print the global frame Cartesian atomic multipoles */

    io___34.ciunit = iounit_1.iout;
    s_wsfe(&io___34);
    e_wsfe();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___35.ciunit = iounit_1.iout;
	s_wsfe(&io___35);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
	do_fio(&c__1, (char *)&atmtyp_1.atomic[i__ - 1], (ftnlen)sizeof(
		integer));
	e_wsfe();
	io___36.ciunit = iounit_1.iout;
	s_wsfe(&io___36);
	do_fio(&c__1, (char *)&atoms_1.x[i__ - 1], (ftnlen)sizeof(doublereal))
		;
	do_fio(&c__1, (char *)&atoms_1.y[i__ - 1], (ftnlen)sizeof(doublereal))
		;
	do_fio(&c__1, (char *)&atoms_1.z__[i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_wsfe();
	io___37.ciunit = iounit_1.iout;
	s_wsfe(&io___37);
	do_fio(&c__1, (char *)&rpole_ref(1, i__), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___38.ciunit = iounit_1.iout;
	s_wsfe(&io___38);
	do_fio(&c__1, (char *)&rpole_ref(2, i__), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rpole_ref(3, i__), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rpole_ref(4, i__), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___39.ciunit = iounit_1.iout;
	s_wsfe(&io___39);
	do_fio(&c__1, (char *)&rpole_ref(5, i__), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___40.ciunit = iounit_1.iout;
	s_wsfe(&io___40);
	do_fio(&c__1, (char *)&rpole_ref(8, i__), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&rpole_ref(9, i__), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___41.ciunit = iounit_1.iout;
	s_wsfe(&io___41);
	do_fio(&c__1, (char *)&rpole_ref(11, i__), (ftnlen)sizeof(doublereal))
		;
	do_fio(&c__1, (char *)&rpole_ref(12, i__), (ftnlen)sizeof(doublereal))
		;
	do_fio(&c__1, (char *)&rpole_ref(13, i__), (ftnlen)sizeof(doublereal))
		;
	e_wsfe();
    }

/*     convert the dipole and quadrupole moments to Angstroms, */
/*     quadrupole divided by 3 for use as traceless values */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (k = 2; k <= 4; ++k) {
	    rpole_ref(k, i__) = rpole_ref(k, i__) * .52917720859;
	}
	for (k = 5; k <= 13; ++k) {
	    rpole_ref(k, i__) = rpole_ref(k, i__) * .28002851809110429 / 3.;
	}
    }
    return 0;
} /* readgdma_ */

#undef rpole_ref
#undef name___ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine match1  --  match first value from GDMA output  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "match1" finds and stores the first multipole component found */
/*     on a line of output from Stone's GDMA program */


/* Subroutine */ int match1_(integer *i__, char *record, ftnlen record_len)
{
    /* System generated locals */
    icilist ici__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);

    /* Local variables */
    extern /* Subroutine */ int match2_(integer *, char *, ftnlen);



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
/*     ##  COPYRIGHT (C)  2005  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  dma.i  --  distributed multipole analysis components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     mp        atomic monopole charge values from DMA */
/*     dpx       atomic dipole moment x-component from DMA */
/*     dpy       atomic dipole moment y-component from DMA */
/*     dpz       atomic dipole moment z-component from DMA */
/*     q20       atomic Q20 quadrupole component from DMA (zz) */
/*     q21c      atomic Q21c quadrupole component from DMA (xz) */
/*     q21s      atomic Q21s quadrupole component from DMA (yz) */
/*     q22c      atomic Q22c quadrupole component from DMA (xx-yy) */
/*     q22s      atomic Q22s quadrupole component from DMA (xy) */




/*     store first multipole component on line of GDMA output */

    if (s_cmp(record + 19, "Q00 ", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 25;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.mp[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
    } else if (s_cmp(record + 19, "Q10 ", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 25;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.dpz[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match2_(i__, record, (ftnlen)120);
    } else if (s_cmp(record + 19, "Q11c", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 25;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.dpx[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match2_(i__, record, (ftnlen)120);
    } else if (s_cmp(record + 19, "Q11s", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 25;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.dpy[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match2_(i__, record, (ftnlen)120);
    } else if (s_cmp(record + 19, "Q20 ", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 25;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.q20[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match2_(i__, record, (ftnlen)120);
    } else if (s_cmp(record + 19, "Q21c", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 25;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.q21c[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match2_(i__, record, (ftnlen)120);
    } else if (s_cmp(record + 19, "Q21s", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 25;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.q21s[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match2_(i__, record, (ftnlen)120);
    } else if (s_cmp(record + 19, "Q22c", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 25;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.q22c[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match2_(i__, record, (ftnlen)120);
    } else if (s_cmp(record + 19, "Q22s", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 25;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.q22s[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match2_(i__, record, (ftnlen)120);
    }
    return 0;
} /* match1_ */



/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine match2  --  match second value from GDMA output  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "match2" finds and stores the second multipole component found */
/*     on a line of output from Stone's GDMA program */


/* Subroutine */ int match2_(integer *i__, char *record, ftnlen record_len)
{
    /* System generated locals */
    icilist ici__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);

    /* Local variables */
    extern /* Subroutine */ int match3_(integer *, char *, ftnlen);



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
/*     ##  COPYRIGHT (C)  2005  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  dma.i  --  distributed multipole analysis components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     mp        atomic monopole charge values from DMA */
/*     dpx       atomic dipole moment x-component from DMA */
/*     dpy       atomic dipole moment y-component from DMA */
/*     dpz       atomic dipole moment z-component from DMA */
/*     q20       atomic Q20 quadrupole component from DMA (zz) */
/*     q21c      atomic Q21c quadrupole component from DMA (xz) */
/*     q21s      atomic Q21s quadrupole component from DMA (yz) */
/*     q22c      atomic Q22c quadrupole component from DMA (xx-yy) */
/*     q22s      atomic Q22s quadrupole component from DMA (xy) */




/*     store second multipole component on line of GDMA output */

    if (s_cmp(record + 38, "Q11c", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 44;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.dpx[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match3_(i__, record, (ftnlen)120);
    } else if (s_cmp(record + 38, "Q11s", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 44;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.dpy[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match3_(i__, record, (ftnlen)120);
    } else if (s_cmp(record + 38, "Q21c", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 44;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.q21c[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match3_(i__, record, (ftnlen)120);
    } else if (s_cmp(record + 38, "Q21s", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 44;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.q21s[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match3_(i__, record, (ftnlen)120);
    } else if (s_cmp(record + 38, "Q22c", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 44;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.q22c[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match3_(i__, record, (ftnlen)120);
    } else if (s_cmp(record + 38, "Q22s", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 44;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.q22s[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
	match3_(i__, record, (ftnlen)120);
    }
    return 0;
} /* match2_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine match3  --  match third value from GDMA output  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "match3" finds and stores the third multipole component found */
/*     on a line of output from Stone's GDMA program */


/* Subroutine */ int match3_(integer *i__, char *record, ftnlen record_len)
{
    /* System generated locals */
    icilist ici__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);



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
/*     ##  COPYRIGHT (C)  2005  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  dma.i  --  distributed multipole analysis components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     mp        atomic monopole charge values from DMA */
/*     dpx       atomic dipole moment x-component from DMA */
/*     dpy       atomic dipole moment y-component from DMA */
/*     dpz       atomic dipole moment z-component from DMA */
/*     q20       atomic Q20 quadrupole component from DMA (zz) */
/*     q21c      atomic Q21c quadrupole component from DMA (xz) */
/*     q21s      atomic Q21s quadrupole component from DMA (yz) */
/*     q22c      atomic Q22c quadrupole component from DMA (xx-yy) */
/*     q22s      atomic Q22s quadrupole component from DMA (xy) */




/*     store third multipole component on line of GDMA output */

    if (s_cmp(record + 57, "Q11s", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 63;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.dpy[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
    } else if (s_cmp(record + 57, "Q21s", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 63;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.q21s[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
    } else if (s_cmp(record + 57, "Q22c", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 63;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.q22c[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
    } else if (s_cmp(record + 57, "Q22s", (ftnlen)4, (ftnlen)4) == 0) {
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 11;
	ici__1.iciunit = record + 63;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&dma_1.q22s[*i__ - 1], (ftnlen)sizeof(
		doublereal));
	e_rsli();
    }
    return 0;
} /* match3_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine molsetup  --  set molecule for polarization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "molsetup" generates trial parameters needed to perform */
/*     polarizable multipole calculations on a structure read */
/*     from a GDMA output file */


/* Subroutine */ int molsetup_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern integer freeunit_(void);
    static integer i__, j;
    static doublereal ri, xi, yi, zi, xr, yr, zr, dij, rad[25000], rij;
    static integer size;
    extern /* Subroutine */ int sort_(integer *, integer *);
    static integer ixyz;
    static doublereal sixth;
    static integer atmnum;
    extern /* Subroutine */ int prtxyz_(integer *);


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define story_ref(a_0,a_1) &atmtyp_1.story[(a_1)*24 + a_0 - 24]
#define polaxe_ref(a_0,a_1) &mpole_1.polaxe[(a_1)*8 + a_0 - 8]



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




/*     set base atomic radii from covalent radius values */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
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

/*     assign atom connectivities based on interatomic distances */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	couple_1.n12[i__ - 1] = 0;
    }
    i__1 = atoms_1.n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	ri = rad[i__ - 1];
	i__2 = atoms_1.n;
	for (j = i__ + 1; j <= i__2; ++j) {
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
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sort_(&couple_1.n12[i__ - 1], &i12_ref(1, i__));
    }

/*     assign unique atom types and set the valence values */

    size = min(20,files_1.leng);
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atoms_1.type__[i__ - 1] = i__;
	atmtyp_1.valence[i__ - 1] = couple_1.n12[i__ - 1];
	s_copy(story_ref(0, i__), files_1.filename, (ftnlen)24, size);
    }

/*     create a file with coordinates and connectivities */

    ixyz = freeunit_();
    prtxyz_(&ixyz);

/*     assign atomic mass and polarizability by atomic number */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atmnum = atmtyp_1.atomic[i__ - 1];
	atmtyp_1.mass[i__ - 1] = 1.;
	polar_1.polarity[i__ - 1] = 0.;
	polar_1.thole[i__ - 1] = .39;
	if (atmnum == 1) {
	    atmtyp_1.mass[i__ - 1] = 1.008;
	    polar_1.polarity[i__ - 1] = .496;
	} else if (atmnum == 5) {
	    atmtyp_1.mass[i__ - 1] = 10.81;
	    polar_1.polarity[i__ - 1] = 1.6;
	} else if (atmnum == 6) {
	    atmtyp_1.mass[i__ - 1] = 12.011;
	    polar_1.polarity[i__ - 1] = 1.334;
	} else if (atmnum == 7) {
	    atmtyp_1.mass[i__ - 1] = 14.007;
	    polar_1.polarity[i__ - 1] = 1.073;
	} else if (atmnum == 8) {
	    atmtyp_1.mass[i__ - 1] = 15.999;
	    polar_1.polarity[i__ - 1] = .837;
	} else if (atmnum == 9) {
	    atmtyp_1.mass[i__ - 1] = 18.998;
	} else if (atmnum == 14) {
	    atmtyp_1.mass[i__ - 1] = 28.086;
	} else if (atmnum == 15) {
	    atmtyp_1.mass[i__ - 1] = 30.974;
	    polar_1.polarity[i__ - 1] = 1.828;
	} else if (atmnum == 16) {
	    atmtyp_1.mass[i__ - 1] = 32.066;
	    polar_1.polarity[i__ - 1] = 3.3;
	} else if (atmnum == 17) {
	    atmtyp_1.mass[i__ - 1] = 35.453;
	    polar_1.polarity[i__ - 1] = 4.;
	} else if (atmnum == 35) {
	    atmtyp_1.mass[i__ - 1] = 79.904;
	    polar_1.polarity[i__ - 1] = 5.65;
	} else if (atmnum == 53) {
	    atmtyp_1.mass[i__ - 1] = 126.904;
	    polar_1.polarity[i__ - 1] = 7.25;
	}
    }

/*     set atomic multipole and polarizability scaling values */

    mpole_1.npole = atoms_1.n;
    polar_1.npolar = atoms_1.n;
    sixth = .16666666666666666;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mpole_1.ipole[i__ - 1] = i__;
	mpole_1.pollist[i__ - 1] = i__;
	mpole_1.polsiz[i__ - 1] = 13;
	s_copy(polaxe_ref(0, i__), "Z-then-X", (ftnlen)8, (ftnlen)8);
	if (polar_1.thole[i__ - 1] == 0.) {
	    polar_1.pdamp[i__ - 1] = 0.;
	} else {
	    polar_1.pdamp[i__ - 1] = pow_dd(&polar_1.polarity[i__ - 1], &
		    sixth);
	}
    }
    return 0;
} /* molsetup_ */

#undef polaxe_ref
#undef story_ref
#undef i12_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine setframe  --  define local coordinate frames  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "setframe" assigns a local coordinate frame at each atomic */
/*     multipole site using high priority connected atoms along axes */


/* Subroutine */ int setframe_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Local Frame Definition for Multipole Sit"
	    "es :\002)";
    static char fmt_20[] = "(/,5x,\002Atom\002,5x,\002Name\002,6x,\002Axis T"
	    "ype\002,5x,\002Z Axis\002,2x,\002X Axis\002,2x,\002Y Axis\002,/)";
    static char fmt_30[] = "(i8,6x,a3,7x,a8,2x,3i8)";
    static char fmt_40[] = "(/,\002 Enter Altered Local Frame Definitio"
	    "n\002,\002 [<CR>=Exit] :  \002,$)";
    static char fmt_50[] = "(a120)";
    static char fmt_70[] = "(/,\002 Local Frame Definition for Multipole Sit"
	    "es :\002)";
    static char fmt_80[] = "(/,5x,\002Atom\002,5x,\002Name\002,6x,\002Axis T"
	    "ype\002,5x,\002Z Axis\002,2x,\002X Axis\002,2x,\002Y Axis\002,/)";
    static char fmt_90[] = "(i8,6x,a3,7x,a8,2x,3i8)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen),
	     s_rsfe(cilist *), e_rsfe(void), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);

    /* Local variables */
    extern integer priority_(integer *, integer *, integer *);
    static integer i__, j, k, m, ia, ib, kb, ic, id, kab, kac, kad, kbc, kbd, 
	    kcd, big;
    static logical alter, query;
    static char record[120];
    extern doublereal random_(void);

    /* Fortran I/O blocks */
    static cilist io___73 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___75 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___78 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___79 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___81 = { 1, record, 1, 0, 120, 1 };
    static cilist io___82 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___83 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_90, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define polaxe_ref(a_0,a_1) &mpole_1.polaxe[(a_1)*8 + a_0 - 8]



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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




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




/*     automatically assign local frame for an isolated atom */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = couple_1.n12[i__ - 1];
	if (j == 0) {
	    mpole_1.zaxis[i__ - 1] = 0;
	    mpole_1.xaxis[i__ - 1] = 0;
	    mpole_1.yaxis[i__ - 1] = 0;

/*     automatically assign local frame for a monovalent atom */

	} else if (j == 1) {
	    ia = i12_ref(1, i__);
	    mpole_1.zaxis[i__ - 1] = ia;
	    if (couple_1.n12[ia - 1] == 1) {
		mpole_1.xaxis[i__ - 1] = 0;
	    } else {
		m = 0;
		i__2 = couple_1.n12[ia - 1];
		for (k = 1; k <= i__2; ++k) {
		    ib = i12_ref(k, ia);
		    kb = atmtyp_1.atomic[ib - 1];
		    if (kb > m && ib != i__) {
			mpole_1.xaxis[i__ - 1] = ib;
			m = kb;
		    }
		}
	    }
	    mpole_1.yaxis[i__ - 1] = 0;

/*     automatically assign local frame for a divalent atom */

	} else if (j == 2) {
	    ia = i12_ref(1, i__);
	    ib = i12_ref(2, i__);
	    kab = priority_(&i__, &ia, &ib);
	    if (kab == ia) {
		mpole_1.zaxis[i__ - 1] = ia;
		mpole_1.xaxis[i__ - 1] = ib;
	    } else if (kab == ib) {
		mpole_1.zaxis[i__ - 1] = ib;
		mpole_1.xaxis[i__ - 1] = ia;
	    } else {
		mpole_1.zaxis[i__ - 1] = ia;
		mpole_1.xaxis[i__ - 1] = ib;
		s_copy(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8);
	    }
	    mpole_1.yaxis[i__ - 1] = 0;

/*     automatically assign local frame for a trivalent atom */

	} else if (j == 3) {
	    ia = i12_ref(1, i__);
	    ib = i12_ref(2, i__);
	    ic = i12_ref(3, i__);
	    kab = priority_(&i__, &ia, &ib);
	    kac = priority_(&i__, &ia, &ic);
	    kbc = priority_(&i__, &ib, &ic);
	    if (kab == 0 && kac == 0) {
		mpole_1.zaxis[i__ - 1] = ia;
		mpole_1.xaxis[i__ - 1] = ib;
	    } else if (kab == ia && kac == ia) {
		mpole_1.zaxis[i__ - 1] = ia;
		mpole_1.xaxis[i__ - 1] = ib;
		if (kbc == ic) {
		    mpole_1.xaxis[i__ - 1] = ic;
		}
	    } else if (kab == ib && kbc == ib) {
		mpole_1.zaxis[i__ - 1] = ib;
		mpole_1.xaxis[i__ - 1] = ia;
		if (kac == ic) {
		    mpole_1.xaxis[i__ - 1] = ic;
		}
	    } else if (kac == ic && kbc == ic) {
		mpole_1.zaxis[i__ - 1] = ic;
		mpole_1.xaxis[i__ - 1] = ia;
		if (kab == ib) {
		    mpole_1.xaxis[i__ - 1] = ib;
		}
	    } else if (kab == 0) {
		mpole_1.zaxis[i__ - 1] = ia;
		mpole_1.xaxis[i__ - 1] = ib;
		s_copy(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8);
	    } else if (kac == 0) {
		mpole_1.zaxis[i__ - 1] = ia;
		mpole_1.xaxis[i__ - 1] = ic;
		s_copy(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8);
	    } else if (kbc == 0) {
		mpole_1.zaxis[i__ - 1] = ib;
		mpole_1.xaxis[i__ - 1] = ic;
		s_copy(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8);
	    }

/*     automatically assign local frame for a tetravalent atom */

	} else if (j == 4) {
	    ia = i12_ref(1, i__);
	    ib = i12_ref(2, i__);
	    ic = i12_ref(3, i__);
	    id = i12_ref(4, i__);
	    kab = priority_(&i__, &ia, &ib);
	    kac = priority_(&i__, &ia, &ic);
	    kad = priority_(&i__, &ia, &id);
	    kbc = priority_(&i__, &ib, &ic);
	    kbd = priority_(&i__, &ib, &id);
	    kcd = priority_(&i__, &ic, &id);
	    if (kab == 0 && kac == 0 && kad == 0) {
		mpole_1.zaxis[i__ - 1] = ia;
		mpole_1.xaxis[i__ - 1] = ib;
	    } else if (kab == ia && kac == ia && kad == ia) {
		mpole_1.zaxis[i__ - 1] = ia;
		mpole_1.xaxis[i__ - 1] = ib;
		if (kbc == ic && kcd == ic) {
		    mpole_1.xaxis[i__ - 1] = ic;
		}
		if (kbd == id && kcd == id) {
		    mpole_1.xaxis[i__ - 1] = id;
		}
		if (kbc == ic && kcd == 0) {
		    mpole_1.xaxis[i__ - 1] = ic;
		}
	    } else if (kab == ib && kbc == ib && kbd == ib) {
		mpole_1.zaxis[i__ - 1] = ib;
		mpole_1.xaxis[i__ - 1] = ia;
		if (kac == ic && kcd == ic) {
		    mpole_1.xaxis[i__ - 1] = ic;
		}
		if (kad == id && kcd == id) {
		    mpole_1.xaxis[i__ - 1] = id;
		}
		if (kac == ic && kcd == 0) {
		    mpole_1.xaxis[i__ - 1] = ic;
		}
	    } else if (kac == ic && kbc == ic && kcd == ic) {
		mpole_1.zaxis[i__ - 1] = ic;
		mpole_1.xaxis[i__ - 1] = ia;
		if (kab == ib && kbd == ib) {
		    mpole_1.xaxis[i__ - 1] = ib;
		}
		if (kad == id && kbd == id) {
		    mpole_1.xaxis[i__ - 1] = id;
		}
		if (kab == ib && kbd == 0) {
		    mpole_1.xaxis[i__ - 1] = ib;
		}
	    } else if (kad == id && kbd == id && kcd == id) {
		mpole_1.zaxis[i__ - 1] = id;
		mpole_1.xaxis[i__ - 1] = ia;
		if (kab == ib && kbc == ib) {
		    mpole_1.xaxis[i__ - 1] = ib;
		}
		if (kac == ic && kbc == ic) {
		    mpole_1.xaxis[i__ - 1] = ic;
		}
		if (kab == ib && kbc == 0) {
		    mpole_1.xaxis[i__ - 1] = ib;
		}
	    } else if (kab == 0 && kac == ia && kad == ia) {
		mpole_1.zaxis[i__ - 1] = ia;
		mpole_1.xaxis[i__ - 1] = ib;
		s_copy(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8);
	    } else if (kac == 0 && kab == ia && kad == ia) {
		mpole_1.zaxis[i__ - 1] = ia;
		mpole_1.xaxis[i__ - 1] = ic;
		s_copy(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8);
	    } else if (kad == 0 && kab == ia && kac == ia) {
		mpole_1.zaxis[i__ - 1] = ia;
		mpole_1.xaxis[i__ - 1] = id;
		s_copy(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8);
	    } else if (kbc == 0 && kab == ib && kbd == ib) {
		mpole_1.zaxis[i__ - 1] = ib;
		mpole_1.xaxis[i__ - 1] = ic;
		s_copy(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8);
	    } else if (kbd == 0 && kab == ib && kbc == ib) {
		mpole_1.zaxis[i__ - 1] = ib;
		mpole_1.xaxis[i__ - 1] = id;
		s_copy(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8);
	    } else if (kcd == 0 && kac == ic && kbc == ic) {
		mpole_1.zaxis[i__ - 1] = ic;
		mpole_1.xaxis[i__ - 1] = id;
		s_copy(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8);
	    }
	}
    }

/*     list the local frame definition for each multipole site */

    io___73.ciunit = iounit_1.iout;
    s_wsfe(&io___73);
    e_wsfe();
    io___74.ciunit = iounit_1.iout;
    s_wsfe(&io___74);
    e_wsfe();
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___75.ciunit = iounit_1.iout;
	s_wsfe(&io___75);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
	do_fio(&c__1, polaxe_ref(0, i__), (ftnlen)8);
	do_fio(&c__1, (char *)&mpole_1.zaxis[i__ - 1], (ftnlen)sizeof(integer)
		);
	do_fio(&c__1, (char *)&mpole_1.xaxis[i__ - 1], (ftnlen)sizeof(integer)
		);
	do_fio(&c__1, (char *)&mpole_1.yaxis[i__ - 1], (ftnlen)sizeof(integer)
		);
	e_wsfe();
    }

/*     allow the user to manually alter local coordinate frames */

    query = TRUE_;
    alter = FALSE_;
    while(query) {
	i__ = 0;
	ia = 0;
	ib = 0;
	ic = 0;
	io___78.ciunit = iounit_1.iout;
	s_wsfe(&io___78);
	e_wsfe();
	io___79.ciunit = iounit_1.input;
	s_rsfe(&io___79);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___81);
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L60;
	}
L60:
	if (i__ == 0) {
	    query = FALSE_;
	} else {
	    alter = TRUE_;
	    if (ia > 0 && ib > 0) {
		s_copy(polaxe_ref(0, i__), "Z-then-X", (ftnlen)8, (ftnlen)8);
	    }
	    if (ia < 0 || ib < 0) {
		s_copy(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8);
	    }
	    if (ib < 0 && ic < 0) {
		s_copy(polaxe_ref(0, i__), "Z-Bisect", (ftnlen)8, (ftnlen)8);
	    }
/* Computing MAX */
	    i__1 = max(ia,ib);
	    if (max(i__1,ic) < 0) {
		s_copy(polaxe_ref(0, i__), "3-Fold", (ftnlen)8, (ftnlen)6);
	    }
	    mpole_1.zaxis[i__ - 1] = abs(ia);
	    mpole_1.xaxis[i__ - 1] = abs(ib);
	    mpole_1.yaxis[i__ - 1] = abs(ic);
	}
    }

/*     repeat local frame list if definitions were altered */

    if (alter) {
	io___82.ciunit = iounit_1.iout;
	s_wsfe(&io___82);
	e_wsfe();
	io___83.ciunit = iounit_1.iout;
	s_wsfe(&io___83);
	e_wsfe();
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___84.ciunit = iounit_1.iout;
	    s_wsfe(&io___84);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
	    do_fio(&c__1, polaxe_ref(0, i__), (ftnlen)8);
	    do_fio(&c__1, (char *)&mpole_1.zaxis[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&mpole_1.xaxis[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&mpole_1.yaxis[i__ - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	}
    }

/*     create random dummy atoms to define local frames as needed */

    big = 0;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (mpole_1.zaxis[i__ - 1] == 0) {
	    mpole_1.zaxis[i__ - 1] = atoms_1.n + 1;
	    mpole_1.xaxis[i__ - 1] = atoms_1.n + 2;
	    big = atoms_1.n + 2;
	} else if (mpole_1.xaxis[i__ - 1] == 0) {
	    mpole_1.xaxis[i__ - 1] = atoms_1.n + 1;
/* Computing MAX */
	    i__2 = big, i__3 = atoms_1.n + 1;
	    big = max(i__2,i__3);
	}
    }
    if (big > atoms_1.n) {
	i__1 = big;
	for (i__ = atoms_1.n + 1; i__ <= i__1; ++i__) {
	    atoms_1.x[i__ - 1] = random_();
	    atoms_1.y[i__ - 1] = random_();
	    atoms_1.z__[i__ - 1] = random_();
	}
    }
    return 0;
} /* setframe_ */

#undef polaxe_ref
#undef name___ref
#undef i12_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function priority  --  atom priority for axis assignment  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "priority" decides which of two connected atoms should be */
/*     preferred during construction of a local coordinate frame */
/*     and returns that atom number; if the two atoms have equal */
/*     priority then a zero is returned */


integer priority_(integer *i__, integer *ia, integer *ib)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer k, m, ka, kb;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]



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




/*     get priority based on atomic number and connected atoms */

    ka = atmtyp_1.atomic[*ia - 1];
    kb = atmtyp_1.atomic[*ib - 1];
    if (ka > kb) {
	ret_val = *ia;
    } else if (kb > ka) {
	ret_val = *ib;
    } else {
	ka = 0;
	i__1 = couple_1.n12[*ia - 1];
	for (k = 1; k <= i__1; ++k) {
	    m = i12_ref(k, *ia);
	    if (*i__ != m) {
		m = atmtyp_1.atomic[m - 1];
		if (m > ka) {
		    ka = m;
		}
	    }
	}
	kb = 0;
	i__1 = couple_1.n12[*ib - 1];
	for (k = 1; k <= i__1; ++k) {
	    m = i12_ref(k, *ib);
	    if (*i__ != m) {
		m = atmtyp_1.atomic[m - 1];
		if (m > kb) {
		    kb = m;
		}
	    }
	}
	if (couple_1.n12[*ia - 1] < couple_1.n12[*ib - 1]) {
	    ret_val = *ia;
	} else if (couple_1.n12[*ib - 1] < couple_1.n12[*ia - 1]) {
	    ret_val = *ib;
	} else if (ka > kb) {
	    ret_val = *ia;
	} else if (kb > ka) {
	    ret_val = *ib;
	} else {
	    ret_val = 0;
	}
    }
    return ret_val;
} /* priority_ */

#undef i12_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine rotframe  --  convert multipoles to local frame  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "rotframe" takes the global multipole moments and rotates them */
/*     into the local coordinate frame defined at each atomic site */


/* Subroutine */ int rotframe_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Local Frame Cartesian Multipole Moment"
	    "s :\002)";
    static char fmt_20[] = "(/,\002 Atom:\002,i8,9x,\002Name:\002,3x,a3,7x"
	    ",\002Atomic Number:\002,i8)";
    static char fmt_30[] = "(/,\002 No Atomic Multipole Moments for this S"
	    "ite\002)";
    static char fmt_40[] = "(/,\002 Atom:\002,i8,9x,\002Name:\002,3x,a3,7x"
	    ",\002Atomic Number:\002,i8)";
    static char fmt_50[] = "(/,\002 Local Frame:\002,12x,a8,6x,3i8)";
    static char fmt_60[] = "(/,\002 Charge:\002,10x,f15.5)";
    static char fmt_70[] = "(\002 Dipole:\002,10x,3f15.5)";
    static char fmt_80[] = "(\002 Quadrupole:\002,6x,f15.5)";
    static char fmt_90[] = "(18x,2f15.5)";
    static char fmt_100[] = "(18x,3f15.5)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal a[9]	/* was [3][3] */;
    static integer i__, j, ii, ixaxe, iyaxe, izaxe;
    extern /* Subroutine */ int rotmat_(integer *, doublereal *), invert_(
	    integer *, integer *, doublereal *), chkpole_(void), rotsite_(
	    integer *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___93 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___95 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___96 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___100 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___101 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___102 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___103 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___104 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___105 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___106 = { 0, 0, 0, fmt_100, 0 };



#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define pole_ref(a_1,a_2) mpole_1.pole[(a_2)*13 + a_1 - 14]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define polaxe_ref(a_0,a_1) &mpole_1.polaxe[(a_1)*8 + a_0 - 8]



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




/*     store the global multipoles in the local frame array */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 13; ++j) {
	    pole_ref(j, i__) = rpole_ref(j, i__);
	}
    }

/*     rotate the multipoles from global frame to local frame */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rotmat_(&i__, a);
	invert_(&c__3, &c__3, a);
	rotsite_(&i__, a);
    }

/*     copy the rotated multipoles back to local frame array */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 13; ++j) {
	    pole_ref(j, i__) = rpole_ref(j, i__);
	}
    }

/*     check the sign of multipole components at chiral sites */

    chkpole_();

/*     convert dipole and quadrupole moments back to atomic units */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rpole_ref(1, i__) = pole_ref(1, i__);
	for (j = 2; j <= 4; ++j) {
	    rpole_ref(j, i__) = pole_ref(j, i__) / .52917720859;
	}
	for (j = 5; j <= 13; ++j) {
	    rpole_ref(j, i__) = pole_ref(j, i__) * 3. / .28002851809110429;
	}
    }

/*     print the local frame Cartesian atomic multipoles */

    io___93.ciunit = iounit_1.iout;
    s_wsfe(&io___93);
    e_wsfe();
    i__1 = mpole_1.npole;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = mpole_1.pollist[ii - 1];
	if (i__ == 0) {
	    io___95.ciunit = iounit_1.iout;
	    s_wsfe(&io___95);
	    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
	    do_fio(&c__1, (char *)&atmtyp_1.atomic[ii - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	    io___96.ciunit = iounit_1.iout;
	    s_wsfe(&io___96);
	    e_wsfe();
	} else {
	    izaxe = mpole_1.zaxis[i__ - 1];
	    ixaxe = mpole_1.xaxis[i__ - 1];
	    iyaxe = mpole_1.yaxis[i__ - 1];
	    if (izaxe > atoms_1.n) {
		izaxe = 0;
	    }
	    if (ixaxe > atoms_1.n) {
		ixaxe = 0;
	    }
	    if (iyaxe < 0) {
		iyaxe = -iyaxe;
	    }
	    io___100.ciunit = iounit_1.iout;
	    s_wsfe(&io___100);
	    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
	    do_fio(&c__1, (char *)&atmtyp_1.atomic[ii - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	    io___101.ciunit = iounit_1.iout;
	    s_wsfe(&io___101);
	    do_fio(&c__1, polaxe_ref(0, i__), (ftnlen)8);
	    do_fio(&c__1, (char *)&izaxe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ixaxe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iyaxe, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___102.ciunit = iounit_1.iout;
	    s_wsfe(&io___102);
	    do_fio(&c__1, (char *)&rpole_ref(1, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___103.ciunit = iounit_1.iout;
	    s_wsfe(&io___103);
	    do_fio(&c__1, (char *)&rpole_ref(2, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&rpole_ref(3, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&rpole_ref(4, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___104.ciunit = iounit_1.iout;
	    s_wsfe(&io___104);
	    do_fio(&c__1, (char *)&rpole_ref(5, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___105.ciunit = iounit_1.iout;
	    s_wsfe(&io___105);
	    do_fio(&c__1, (char *)&rpole_ref(8, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&rpole_ref(9, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___106.ciunit = iounit_1.iout;
	    s_wsfe(&io___106);
	    do_fio(&c__1, (char *)&rpole_ref(11, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&rpole_ref(12, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&rpole_ref(13, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* rotframe_ */

#undef polaxe_ref
#undef rpole_ref
#undef pole_ref
#undef name___ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine fixframe  --  alter the local frame definition  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "fixframe" is a service routine that alters the local frame */
/*     definition for specified atoms */


/* Subroutine */ int fixframe_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Local Frame Definition for Multipole Sit"
	    "es :\002)";
    static char fmt_20[] = "(/,5x,\002Atom\002,5x,\002Name\002,6x,\002Axis T"
	    "ype\002,5x,\002Z Axis\002,2x,\002X Axis\002,2x,\002Y Axis\002,/)";
    static char fmt_30[] = "(i8,6x,a3,10x,\002--\002,11x,\002--\002,6x,"
	    "\002--\002,6x,\002--\002)";
    static char fmt_40[] = "(i8,6x,a3,7x,a8,2x,3i8)";
    static char fmt_50[] = "(/,\002 Enter Altered Local Frame Definitio"
	    "n\002,\002 [<CR>=Exit] :  \002,$)";
    static char fmt_60[] = "(a120)";
    static char fmt_80[] = "(/,\002 Local Frame Definition for Multipole Sit"
	    "es :\002)";
    static char fmt_90[] = "(/,5x,\002Atom\002,5x,\002Name\002,6x,\002Axis T"
	    "ype\002,5x,\002Z Axis\002,2x,\002X Axis\002,2x,\002Y Axis\002,/)";
    static char fmt_100[] = "(i8,6x,a3,10x,\002--\002,11x,\002--\002,6x,\002"
	    "--\002,6x,\002--\002)";
    static char fmt_110[] = "(i8,6x,a3,7x,a8,2x,3i8)";
    static char fmt_130[] = "(/,\002 Multipoles With Altered Local Frame Def"
	    "inition :\002)";
    static char fmt_140[] = "(/,\002 Atom:\002,i8,9x,\002Name:\002,3x,a3,7x"
	    ",\002Atomic Number:\002,i8)";
    static char fmt_150[] = "(/,\002 No Atomic Multipole Moments for this Si"
	    "te\002)";
    static char fmt_160[] = "(/,\002 Atom:\002,i8,9x,\002Name:\002,3x,a3,7x"
	    ",\002Atomic Number:\002,i8)";
    static char fmt_170[] = "(/,\002 Local Frame:\002,12x,a8,6x,3i8)";
    static char fmt_180[] = "(/,\002 Charge:\002,10x,f15.5)";
    static char fmt_190[] = "(\002 Dipole:\002,10x,3f15.5)";
    static char fmt_200[] = "(\002 Quadrupole:\002,6x,f15.5)";
    static char fmt_210[] = "(18x,2f15.5)";
    static char fmt_220[] = "(18x,3f15.5)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen),
	     s_rsfe(cilist *), e_rsfe(void), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_dnnt(doublereal *);

    /* Local variables */
    static doublereal a[9]	/* was [3][3] */;
    static integer i__, j, k, ia, ib, ic, ii;
    static doublereal ci, cj, big, eps, sum;
    static integer ixaxe, iyaxe, izaxe;
    static logical alter, query;
    static char record[120];
    extern /* Subroutine */ int rotmat_(integer *, doublereal *), invert_(
	    integer *, integer *, doublereal *), chkpole_(void), rotpole_(
	    void), rotsite_(integer *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___107 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___108 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___111 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___115 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___121 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___122 = { 0, 0, 0, fmt_60, 0 };
    static icilist io___124 = { 1, record, 1, 0, 120, 1 };
    static cilist io___125 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___126 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___127 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___128 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___137 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___138 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___139 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___140 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___141 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___142 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___143 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___144 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___145 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___146 = { 0, 0, 0, fmt_220, 0 };



#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define pole_ref(a_1,a_2) mpole_1.pole[(a_2)*13 + a_1 - 14]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define polaxe_ref(a_0,a_1) &mpole_1.polaxe[(a_1)*8 + a_0 - 8]



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
/*     ##  kpolr.i  --  forcefield parameters for polarizability  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     polr   dipole polarizability parameters for each atom type */
/*     athl   Thole polarizability damping value for each atom type */
/*     pgrp   connected types in polarization group of each atom type */




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




/*     rotate the multipole components into the global frame */

    rotpole_();

/*     list the local frame definition for each multipole site */

    io___107.ciunit = iounit_1.iout;
    s_wsfe(&io___107);
    e_wsfe();
    io___108.ciunit = iounit_1.iout;
    s_wsfe(&io___108);
    e_wsfe();
    i__1 = atoms_1.n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = mpole_1.pollist[ii - 1];
	if (i__ == 0) {
	    io___111.ciunit = iounit_1.iout;
	    s_wsfe(&io___111);
	    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
	    e_wsfe();
	} else {
	    izaxe = mpole_1.zaxis[i__ - 1];
	    ixaxe = mpole_1.xaxis[i__ - 1];
	    iyaxe = mpole_1.yaxis[i__ - 1];
	    if (izaxe > atoms_1.n) {
		izaxe = 0;
	    }
	    if (ixaxe > atoms_1.n) {
		ixaxe = 0;
	    }
	    if (iyaxe < 0) {
		iyaxe = -iyaxe;
	    }
	    io___115.ciunit = iounit_1.iout;
	    s_wsfe(&io___115);
	    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
	    do_fio(&c__1, polaxe_ref(0, i__), (ftnlen)8);
	    do_fio(&c__1, (char *)&izaxe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ixaxe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iyaxe, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }

/*     allow the user to manually alter local coordinate frames */

    query = TRUE_;
    alter = FALSE_;
    while(query) {
	ii = 0;
	ia = 0;
	ib = 0;
	ic = 0;
	io___121.ciunit = iounit_1.iout;
	s_wsfe(&io___121);
	e_wsfe();
	io___122.ciunit = iounit_1.input;
	s_rsfe(&io___122);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___124);
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ii, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L70;
	}
L70:
	i__ = mpole_1.pollist[ii - 1];
	if (ii == 0) {
	    query = FALSE_;
	} else if (i__ != 0) {
	    alter = TRUE_;
	    if (ia > 0 && ib > 0) {
		s_copy(polaxe_ref(0, i__), "Z-then-X", (ftnlen)8, (ftnlen)8);
	    }
	    if (ia < 0 || ib < 0) {
		s_copy(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8);
	    }
	    if (ib < 0 && ic < 0) {
		s_copy(polaxe_ref(0, i__), "Z-Bisect", (ftnlen)8, (ftnlen)8);
	    }
/* Computing MAX */
	    i__1 = max(ia,ib);
	    if (max(i__1,ic) < 0) {
		s_copy(polaxe_ref(0, i__), "3-Fold", (ftnlen)8, (ftnlen)6);
	    }
	    mpole_1.zaxis[i__ - 1] = abs(ia);
	    mpole_1.xaxis[i__ - 1] = abs(ib);
	    mpole_1.yaxis[i__ - 1] = abs(ic);
	}
    }

/*     repeat local frame list if definitions were altered */

    if (alter) {
	io___125.ciunit = iounit_1.iout;
	s_wsfe(&io___125);
	e_wsfe();
	io___126.ciunit = iounit_1.iout;
	s_wsfe(&io___126);
	e_wsfe();
	i__1 = mpole_1.npole;
	for (ii = 1; ii <= i__1; ++ii) {
	    i__ = mpole_1.pollist[ii - 1];
	    if (i__ == 0) {
		io___127.ciunit = iounit_1.iout;
		s_wsfe(&io___127);
		do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
		e_wsfe();
	    } else {
		izaxe = mpole_1.zaxis[i__ - 1];
		ixaxe = mpole_1.xaxis[i__ - 1];
		iyaxe = mpole_1.yaxis[i__ - 1];
		if (izaxe > atoms_1.n) {
		    izaxe = 0;
		}
		if (ixaxe > atoms_1.n) {
		    ixaxe = 0;
		}
		if (iyaxe < 0) {
		    iyaxe = -iyaxe;
		}
		io___128.ciunit = iounit_1.iout;
		s_wsfe(&io___128);
		do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
		do_fio(&c__1, polaxe_ref(0, i__), (ftnlen)8);
		do_fio(&c__1, (char *)&izaxe, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ixaxe, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iyaxe, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	}
    }

/*     store the global multipoles in the local frame array */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 13; ++j) {
	    pole_ref(j, i__) = rpole_ref(j, i__);
	}
    }

/*     rotate the multipoles from global frame to local frame */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rotmat_(&i__, a);
	invert_(&c__3, &c__3, a);
	rotsite_(&i__, a);
    }

/*     copy the rotated multipoles back to local frame array */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 13; ++j) {
	    pole_ref(j, i__) = rpole_ref(j, i__);
	}
    }

/*     check the sign of multipole components at chiral sites */

    chkpole_();

/*     convert dipole and quadrupole moments back to atomic units */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pole_ref(1, i__) = pole_ref(1, i__);
	for (j = 2; j <= 4; ++j) {
	    pole_ref(j, i__) = pole_ref(j, i__) / .52917720859;
	}
	for (j = 5; j <= 13; ++j) {
	    pole_ref(j, i__) = pole_ref(j, i__) * 3. / .28002851809110429;
	}
    }

/*     regularize the multipole moments to desired precision */

    eps = 1e-5;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 13; ++j) {
	    d__1 = pole_ref(j, i__) / eps;
	    pole_ref(j, i__) = (doublereal) i_dnnt(&d__1) * eps;
	}
    }

/*     maintain integer net charge for the whole system */

    k = 0;
    big = 0.;
    sum = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum += pole_ref(1, i__);
	ci = (d__1 = pole_ref(1, i__), abs(d__1));
	if (ci > big) {
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		cj = (d__1 = pole_ref(1, j), abs(d__1));
		if (i__ != j && ci == cj) {
		    goto L120;
		}
	    }
	    k = i__;
	    big = ci;
L120:
	    ;
	}
    }
    sum -= (doublereal) i_dnnt(&sum);
    if (k != 0) {
	pole_ref(1, k) = pole_ref(1, k) - sum;
    }

/*     maintain traceless quadrupole at each multipole site */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = pole_ref(5, i__) + pole_ref(9, i__) + pole_ref(13, i__);
/* Computing MAX */
	d__4 = (d__1 = pole_ref(5, i__), abs(d__1)), d__5 = (d__2 = pole_ref(
		9, i__), abs(d__2)), d__4 = max(d__4,d__5), d__5 = (d__3 = 
		pole_ref(13, i__), abs(d__3));
	big = max(d__4,d__5);
	k = 0;
	if (big == (d__1 = pole_ref(5, i__), abs(d__1))) {
	    k = 5;
	}
	if (big == (d__1 = pole_ref(9, i__), abs(d__1))) {
	    k = 9;
	}
	if (big == (d__1 = pole_ref(13, i__), abs(d__1))) {
	    k = 13;
	}
	if (k != 0) {
	    pole_ref(k, i__) = pole_ref(k, i__) - sum;
	}
    }

/*     print the altered local frame atomic multipole values */

    io___137.ciunit = iounit_1.iout;
    s_wsfe(&io___137);
    e_wsfe();
    i__1 = atoms_1.n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = mpole_1.pollist[ii - 1];
	if (i__ == 0) {
	    io___138.ciunit = iounit_1.iout;
	    s_wsfe(&io___138);
	    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
	    do_fio(&c__1, (char *)&atmtyp_1.atomic[ii - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	    io___139.ciunit = iounit_1.iout;
	    s_wsfe(&io___139);
	    e_wsfe();
	} else {
	    izaxe = mpole_1.zaxis[i__ - 1];
	    ixaxe = mpole_1.xaxis[i__ - 1];
	    iyaxe = mpole_1.yaxis[i__ - 1];
	    if (izaxe > atoms_1.n) {
		izaxe = 0;
	    }
	    if (ixaxe > atoms_1.n) {
		ixaxe = 0;
	    }
	    if (iyaxe < 0) {
		iyaxe = -iyaxe;
	    }
	    io___140.ciunit = iounit_1.iout;
	    s_wsfe(&io___140);
	    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
	    do_fio(&c__1, (char *)&atmtyp_1.atomic[ii - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	    io___141.ciunit = iounit_1.iout;
	    s_wsfe(&io___141);
	    do_fio(&c__1, polaxe_ref(0, i__), (ftnlen)8);
	    do_fio(&c__1, (char *)&izaxe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ixaxe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iyaxe, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___142.ciunit = iounit_1.iout;
	    s_wsfe(&io___142);
	    do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___143.ciunit = iounit_1.iout;
	    s_wsfe(&io___143);
	    do_fio(&c__1, (char *)&pole_ref(2, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(3, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(4, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___144.ciunit = iounit_1.iout;
	    s_wsfe(&io___144);
	    do_fio(&c__1, (char *)&pole_ref(5, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___145.ciunit = iounit_1.iout;
	    s_wsfe(&io___145);
	    do_fio(&c__1, (char *)&pole_ref(8, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(9, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___146.ciunit = iounit_1.iout;
	    s_wsfe(&io___146);
	    do_fio(&c__1, (char *)&pole_ref(11, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(12, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(13, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* fixframe_ */

#undef polaxe_ref
#undef rpole_ref
#undef pole_ref
#undef name___ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine setpolar  --  define polarization and groups  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "setpolar" assigns atomic polarizability and Thole damping */
/*     parameters and allows user alteration of these values */


/* Subroutine */ int setpolar_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Atomic Polarizabilities for Multipole Si"
	    "tes :\002)";
    static char fmt_20[] = "(/,5x,\002Atom\002,5x,\002Name\002,7x,\002Polari"
	    "ze\002,10x,\002Thole\002,/)";
    static char fmt_30[] = "(i8,6x,a3,12x,\002--\002,13x,\002--\002)";
    static char fmt_40[] = "(i8,6x,a3,4x,f12.4,3x,f12.4)";
    static char fmt_50[] = "(/,\002 Enter Atom Number & Polarizability Val"
	    "ues\002,\002 [<CR>=Exit] :  \002,$)";
    static char fmt_60[] = "(a120)";
    static char fmt_80[] = "(/,\002 Atomic Polarizabilities for Multipole Si"
	    "tes :\002)";
    static char fmt_90[] = "(/,5x,\002Atom\002,5x,\002Name\002,7x,\002Polari"
	    "ze\002,10x,\002Thole\002,/)";
    static char fmt_100[] = "(i8,6x,a3,12x,\002--\002,13x,\002--\002)";
    static char fmt_110[] = "(i8,6x,a3,4x,f12.4,3x,f12.4)";
    static char fmt_120[] = "(/,\002 The default is to place all Atoms into "
	    "one\002,\002 Polarization Group;\002,/,\002 This can be altered "
	    "by entering a series of\002,\002 Bonded Atom Pairs\002,/,\002 th"
	    "at separate the Molecule into distinct\002,\002 Polarization Gro"
	    "ups\002)";
    static char fmt_130[] = "(/,\002 Enter a Bond between Polarization Gro"
	    "ups\002,\002 [<CR>=Exit] :  \002,$)";
    static char fmt_140[] = "(a120)";
    static char fmt_160[] = "(/,\002 Polarization Groups for Multipole Sites"
	    " :\002)";
    static char fmt_170[] = "(/,5x,\002Atom\002,5x,\002Name\002,7x,\002Polar"
	    "ization Group\002,\002 Definition\002,/)";
    static char fmt_180[] = "(i8,6x,a3,8x,20i6)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen),
	     s_rsfe(cilist *), e_rsfe(void), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);

    /* Local variables */
    extern /* Subroutine */ int polargrp_(void);
    static integer i__, j, k, ia, ib, ii;
    static doublereal thl, pol;
    static logical alter, query;
    static char record[120];

    /* Fortran I/O blocks */
    static cilist io___147 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___148 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___151 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___152 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___157 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___158 = { 0, 0, 0, fmt_60, 0 };
    static icilist io___160 = { 1, record, 1, 0, 120, 1 };
    static cilist io___161 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___162 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___163 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___164 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___165 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___169 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___170 = { 0, 0, 0, fmt_140, 0 };
    static icilist io___171 = { 1, record, 1, 0, 120, 1 };
    static cilist io___172 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___173 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___175 = { 0, 0, 0, fmt_180, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define pgrp_ref(a_1,a_2) kpolr_1.pgrp[(a_2)*8 + a_1 - 9]



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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  kpolr.i  --  forcefield parameters for polarizability  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     polr   dipole polarizability parameters for each atom type */
/*     athl   Thole polarizability damping value for each atom type */
/*     pgrp   connected types in polarization group of each atom type */




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
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  polgrp.i  --  polarizable site group connectivity lists   ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxp11   maximum number of atoms in a polarization group */
/*     maxp12   maximum number of atoms in groups 1-2 to an atom */
/*     maxp13   maximum number of atoms in groups 1-3 to an atom */
/*     maxp14   maximum number of atoms in groups 1-4 to an atom */

/*     np11     number of atoms in polarization group of each atom */
/*     ip11     atom numbers of atoms in same group as each atom */
/*     np12     number of atoms in groups 1-2 to each atom */
/*     ip12     atom numbers of atoms in groups 1-2 to each atom */
/*     np13     number of atoms in groups 1-3 to each atom */
/*     ip13     atom numbers of atoms in groups 1-3 to each atom */
/*     np14     number of atoms in groups 1-4 to each atom */
/*     ip14     atom numbers of atoms in groups 1-4 to each atom */




/*     list the polariability values for each multipole site */

    io___147.ciunit = iounit_1.iout;
    s_wsfe(&io___147);
    e_wsfe();
    io___148.ciunit = iounit_1.iout;
    s_wsfe(&io___148);
    e_wsfe();
    i__1 = atoms_1.n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = mpole_1.pollist[ii - 1];
	if (i__ == 0) {
	    io___151.ciunit = iounit_1.iout;
	    s_wsfe(&io___151);
	    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
	    e_wsfe();
	} else {
	    io___152.ciunit = iounit_1.iout;
	    s_wsfe(&io___152);
	    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
	    do_fio(&c__1, (char *)&polar_1.polarity[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&polar_1.thole[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }

/*     allow the user to manually alter polarizability values */

    query = TRUE_;
    alter = FALSE_;
    while(query) {
	ii = 0;
	pol = 0.;
	thl = .39;
	io___157.ciunit = iounit_1.iout;
	s_wsfe(&io___157);
	e_wsfe();
	io___158.ciunit = iounit_1.input;
	s_rsfe(&io___158);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___160);
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ii, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&pol, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&thl, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L70;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L70;
	}
L70:
	i__ = mpole_1.pollist[ii - 1];
	if (ii == 0) {
	    query = FALSE_;
	} else if (i__ != 0) {
	    alter = TRUE_;
	    polar_1.polarity[i__ - 1] = pol;
	    polar_1.thole[i__ - 1] = thl;
	}
    }

/*     repeat polarizability values if parameters were altered */

    if (alter) {
	io___161.ciunit = iounit_1.iout;
	s_wsfe(&io___161);
	e_wsfe();
	io___162.ciunit = iounit_1.iout;
	s_wsfe(&io___162);
	e_wsfe();
	i__1 = atoms_1.n;
	for (ii = 1; ii <= i__1; ++ii) {
	    i__ = mpole_1.pollist[ii - 1];
	    if (i__ == 0) {
		io___163.ciunit = iounit_1.iout;
		s_wsfe(&io___163);
		do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
		e_wsfe();
	    } else {
		io___164.ciunit = iounit_1.iout;
		s_wsfe(&io___164);
		do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
		do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
		do_fio(&c__1, (char *)&polar_1.polarity[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&polar_1.thole[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     use bonded atoms as initial guess at polarization groups */

    io___165.ciunit = iounit_1.iout;
    s_wsfe(&io___165);
    e_wsfe();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = couple_1.n12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    pgrp_ref(j, i__) = i12_ref(j, i__);
	}
    }

/*     get the bonds that separate the polarization groups */

    query = TRUE_;
    while(query) {
	ia = 0;
	ib = 0;
	io___169.ciunit = iounit_1.iout;
	s_wsfe(&io___169);
	e_wsfe();
	io___170.ciunit = iounit_1.input;
	s_rsfe(&io___170);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	i__1 = s_rsli(&io___171);
	if (i__1 != 0) {
	    goto L150;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L150;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L150;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L150;
	}
L150:
	if (ia == 0 || ib == 0) {
	    query = FALSE_;
	} else {
	    i__1 = couple_1.n12[ia - 1];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (pgrp_ref(i__, ia) == ib) {
		    i__2 = couple_1.n12[ia - 1];
		    for (j = i__ + 1; j <= i__2; ++j) {
			pgrp_ref(j - 1, ia) = pgrp_ref(j, ia);
		    }
		    pgrp_ref(couple_1.n12[ia - 1], ia) = 0;
		}
	    }
	    i__1 = couple_1.n12[ib - 1];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (pgrp_ref(i__, ib) == ia) {
		    i__2 = couple_1.n12[ib - 1];
		    for (j = i__ + 1; j <= i__2; ++j) {
			pgrp_ref(j - 1, ib) = pgrp_ref(j, ib);
		    }
		    pgrp_ref(couple_1.n12[ib - 1], ib) = 0;
		}
	    }
	}
    }
    polargrp_();

/*     list the polarization group for each multipole site */

    io___172.ciunit = iounit_1.iout;
    s_wsfe(&io___172);
    e_wsfe();
    io___173.ciunit = iounit_1.iout;
    s_wsfe(&io___173);
    e_wsfe();
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = 0;
	for (j = 1; j <= 8; ++j) {
	    if (pgrp_ref(j, i__) != 0) {
		k = j;
	    }
	}
	io___175.ciunit = iounit_1.iout;
	s_wsfe(&io___175);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
	i__2 = k;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&pgrp_ref(j, i__), (ftnlen)sizeof(integer));
	}
	e_wsfe();
    }
    return 0;
} /* setpolar_ */

#undef pgrp_ref
#undef name___ref
#undef i12_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine intrapol  --  alter multipoles for polarization  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "intrapol" finds an output set of TINKER multipole parameters */
/*     which when used with an intergroup polarization model will */
/*     give the same electrostatic potential around the molecule as */
/*     the input set of multipole parameters with all atoms in one */
/*     polarization group */

/*     for example, the input parameters could be from a distributed */
/*     multipole analysis of a molecular wavefunction and the output */
/*     will be the parameter values that achieve the same potential */
/*     in the presence of intergroup (intramolecular) polarization */


/* Subroutine */ int intrapol_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Multipoles after Removal of Intergrou"
	    "p\002,\002 Polarization :\002)";
    static char fmt_20[] = "(/,\002 Atom:\002,i8,9x,\002Name:\002,3x,a3,7x"
	    ",\002Atomic Number:\002,i8)";
    static char fmt_30[] = "(/,\002 No Atomic Multipole Moments for this S"
	    "ite\002)";
    static char fmt_40[] = "(/,\002 Atom:\002,i8,9x,\002Name:\002,3x,a3,7x"
	    ",\002Atomic Number:\002,i8)";
    static char fmt_50[] = "(/,\002 Local Frame:\002,12x,a8,6x,3i8)";
    static char fmt_60[] = "(/,\002 Charge:\002,10x,f15.5)";
    static char fmt_70[] = "(\002 Dipole:\002,10x,3f15.5)";
    static char fmt_80[] = "(\002 Quadrupole:\002,6x,f15.5)";
    static char fmt_90[] = "(18x,2f15.5)";
    static char fmt_100[] = "(18x,3f15.5)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal a[9]	/* was [3][3] */;
    static integer i__, j;
    extern /* Subroutine */ int intrapol1_(void);
    static integer ii, ixaxe, iyaxe, izaxe;
    extern /* Subroutine */ int rotmat_(integer *, doublereal *), invert_(
	    integer *, integer *, doublereal *), rotpole_(void), rotsite_(
	    integer *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___179 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___181 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___182 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___186 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___187 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___188 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___189 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___190 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___191 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___192 = { 0, 0, 0, fmt_100, 0 };



#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define pole_ref(a_1,a_2) mpole_1.pole[(a_2)*13 + a_1 - 14]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define polaxe_ref(a_0,a_1) &mpole_1.polaxe[(a_1)*8 + a_0 - 8]



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




/*     rotate the multipole components into the global frame */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rotpole_();
    }

/*     compute induced dipoles to be removed from QM multipoles */

    intrapol1_();

/*     remove induced dipole from global frame multipoles */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rpole_ref(2, i__) = rpole_ref(2, i__) - uind_ref(1, i__);
	rpole_ref(3, i__) = rpole_ref(3, i__) - uind_ref(2, i__);
	rpole_ref(4, i__) = rpole_ref(4, i__) - uind_ref(3, i__);
	for (j = 1; j <= 13; ++j) {
	    pole_ref(j, i__) = rpole_ref(j, i__);
	}
    }

/*     rotate the multipoles from global frame to local frame */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rotmat_(&i__, a);
	invert_(&c__3, &c__3, a);
	rotsite_(&i__, a);
    }

/*     copy the rotated multipoles back to local frame array */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 13; ++j) {
	    pole_ref(j, i__) = rpole_ref(j, i__);
	}
    }

/*     convert dipole and quadrupole moments back to atomic units */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rpole_ref(1, i__) = pole_ref(1, i__);
	for (j = 2; j <= 4; ++j) {
	    rpole_ref(j, i__) = pole_ref(j, i__) / .52917720859;
	}
	for (j = 5; j <= 13; ++j) {
	    rpole_ref(j, i__) = pole_ref(j, i__) * 3. / .28002851809110429;
	}
    }

/*     print multipoles with intergroup polarization removed */

    io___179.ciunit = iounit_1.iout;
    s_wsfe(&io___179);
    e_wsfe();
    i__1 = atoms_1.n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = mpole_1.pollist[ii - 1];
	if (i__ == 0) {
	    io___181.ciunit = iounit_1.iout;
	    s_wsfe(&io___181);
	    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
	    do_fio(&c__1, (char *)&atmtyp_1.atomic[ii - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	    io___182.ciunit = iounit_1.iout;
	    s_wsfe(&io___182);
	    e_wsfe();
	} else {
	    izaxe = mpole_1.zaxis[i__ - 1];
	    ixaxe = mpole_1.xaxis[i__ - 1];
	    iyaxe = mpole_1.yaxis[i__ - 1];
	    if (izaxe > atoms_1.n) {
		izaxe = 0;
	    }
	    if (ixaxe > atoms_1.n) {
		ixaxe = 0;
	    }
	    if (iyaxe < 0) {
		iyaxe = -iyaxe;
	    }
	    io___186.ciunit = iounit_1.iout;
	    s_wsfe(&io___186);
	    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
	    do_fio(&c__1, (char *)&atmtyp_1.atomic[ii - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	    io___187.ciunit = iounit_1.iout;
	    s_wsfe(&io___187);
	    do_fio(&c__1, polaxe_ref(0, i__), (ftnlen)8);
	    do_fio(&c__1, (char *)&izaxe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ixaxe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iyaxe, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___188.ciunit = iounit_1.iout;
	    s_wsfe(&io___188);
	    do_fio(&c__1, (char *)&rpole_ref(1, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___189.ciunit = iounit_1.iout;
	    s_wsfe(&io___189);
	    do_fio(&c__1, (char *)&rpole_ref(2, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&rpole_ref(3, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&rpole_ref(4, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___190.ciunit = iounit_1.iout;
	    s_wsfe(&io___190);
	    do_fio(&c__1, (char *)&rpole_ref(5, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___191.ciunit = iounit_1.iout;
	    s_wsfe(&io___191);
	    do_fio(&c__1, (char *)&rpole_ref(8, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&rpole_ref(9, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___192.ciunit = iounit_1.iout;
	    s_wsfe(&io___192);
	    do_fio(&c__1, (char *)&rpole_ref(11, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&rpole_ref(12, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&rpole_ref(13, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* intrapol_ */

#undef polaxe_ref
#undef rpole_ref
#undef uind_ref
#undef pole_ref
#undef name___ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine intrapol1  --  induced dipoles via double loop  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "intrapol1" is a service routine that computes induced dipole */
/*     moments for use during removal of intergroup polarization */


/* Subroutine */ int intrapol1_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Determination of Intergroup Induced\002"
	    ",\002 Dipoles :\002,//,4x,\002Iter\002,8x,\002RMS Change (Debyes)"
	    "\002,/)";
    static char fmt_20[] = "(i8,7x,f16.10)";
    static char fmt_30[] = "(/,\002 INTRAPOL1  --  Warning, Induced Dipole"
	    "s\002,\002 are not Converged\002)";
    static char fmt_40[] = "(/,\002 Intergroup Induced Dipoles to be Remove"
	    "d\002,\002 (Debyes) :\002)";
    static char fmt_50[] = "(/,4x,\002Atom\002,14x,\002X\002,11x,\002Y\002,1"
	    "1x,\002Z\002,9x,\002Total\002/)";
    static char fmt_60[] = "(i8,5x,4f12.4)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k;
    static doublereal r__, r2, ci, ck, fi[3], fk[3];
    static integer ii, kk, ix, iz, kx, kz;
    static doublereal xr, yr, zr, rr3, rr5, rr7, pdi, dir, dix, diy, diz, dkx,
	     dky, dkz, dkr, qir, pti, qkr, eps, uir, ukr, qix, qiy, qiz, uix, 
	    uiy, uiz, ukx, uky, ukz, qkx, qky, qkz, damp;
    static logical done;
    static integer iter;
    static doublereal udir[75000]	/* was [3][25000] */, uold[75000]	
	    /* was [3][25000] */, norm, qixx, qixy, qixz, qiyy, qiyz, qizz, 
	    qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, field[75000]	/* was [3][
	    25000] */;
    extern /* Subroutine */ int fatal_(void);
    static doublereal scale3, scale5, scale7, pgamma, pscale[25000], epsold;
    extern /* Subroutine */ int prterr_(void);
    static integer maxiter;

    /* Fortran I/O blocks */
    static cilist io___266 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___267 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___268 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___269 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___270 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___272 = { 0, 0, 0, fmt_60, 0 };



#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define ip12_ref(a_1,a_2) polgrp_1.ip12[(a_2)*50 + a_1 - 51]
#define ip13_ref(a_1,a_2) polgrp_1.ip13[(a_2)*50 + a_1 - 51]
#define ip14_ref(a_1,a_2) polgrp_1.ip14[(a_2)*50 + a_1 - 51]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define udir_ref(a_1,a_2) udir[(a_2)*3 + a_1 - 4]
#define uold_ref(a_1,a_2) uold[(a_2)*3 + a_1 - 4]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]



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
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  polgrp.i  --  polarizable site group connectivity lists   ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxp11   maximum number of atoms in a polarization group */
/*     maxp12   maximum number of atoms in groups 1-2 to an atom */
/*     maxp13   maximum number of atoms in groups 1-3 to an atom */
/*     maxp14   maximum number of atoms in groups 1-4 to an atom */

/*     np11     number of atoms in polarization group of each atom */
/*     ip11     atom numbers of atoms in same group as each atom */
/*     np12     number of atoms in groups 1-2 to each atom */
/*     ip12     atom numbers of atoms in groups 1-2 to each atom */
/*     np13     number of atoms in groups 1-3 to each atom */
/*     ip13     atom numbers of atoms in groups 1-3 to each atom */
/*     np14     number of atoms in groups 1-4 to each atom */
/*     ip14     atom numbers of atoms in groups 1-4 to each atom */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  polpot.i  --  specifics of polarization functional form  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     poleps    induced dipole convergence criterion (rms Debyes/atom) */
/*     polsor    induced dipole SOR convergence acceleration factor */
/*     p2scale   field 1-2 scale factor for energy evaluations */
/*     p3scale   field 1-3 scale factor for energy evaluations */
/*     p4scale   field 1-4 scale factor for energy evaluations */
/*     p5scale   field 1-5 scale factor for energy evaluations */
/*     d1scale   field intra-group scale factor for direct induced */
/*     d2scale   field 1-2 group scale factor for direct induced */
/*     d3scale   field 1-3 group scale factor for direct induced */
/*     d4scale   field 1-4 group scale factor for direct induced */
/*     u1scale   field intra-group scale factor for mutual induced */
/*     u2scale   field 1-2 group scale factor for mutual induced */
/*     u3scale   field 1-3 group scale factor for mutual induced */
/*     u4scale   field 1-4 group scale factor for mutual induced */
/*     poltyp    type of polarization potential (direct or mutual) */




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




/*     zero out the induced dipole and the field at each site */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    uind_ref(j, i__) = 0.;
	    field_ref(j, i__) = 0.;
	}
    }

/*     compute the direct induced dipole moment at each atom */

    i__1 = mpole_1.npole - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	iz = mpole_1.zaxis[ii - 1];
	ix = mpole_1.xaxis[ii - 1];
	pdi = polar_1.pdamp[i__ - 1];
	pti = polar_1.thole[i__ - 1];
	ci = rpole_ref(1, i__);
	dix = rpole_ref(2, i__);
	diy = rpole_ref(3, i__);
	diz = rpole_ref(4, i__);
	qixx = rpole_ref(5, i__);
	qixy = rpole_ref(6, i__);
	qixz = rpole_ref(7, i__);
	qiyy = rpole_ref(9, i__);
	qiyz = rpole_ref(10, i__);
	qizz = rpole_ref(13, i__);
	i__2 = mpole_1.npole;
	for (j = i__ + 1; j <= i__2; ++j) {
	    pscale[mpole_1.ipole[j - 1] - 1] = 1.;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[ip11_ref(j, ii) - 1] = polpot_1.d1scale;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[ip12_ref(j, ii) - 1] = polpot_1.d2scale;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[ip13_ref(j, ii) - 1] = polpot_1.d3scale;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[ip14_ref(j, ii) - 1] = polpot_1.d4scale;
	}
	i__2 = mpole_1.npole;
	for (k = i__ + 1; k <= i__2; ++k) {
	    kk = mpole_1.ipole[k - 1];
	    kz = mpole_1.zaxis[kk - 1];
	    kx = mpole_1.xaxis[kk - 1];
	    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
	    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
	    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
	    r2 = xr * xr + yr * yr + zr * zr;
	    r__ = sqrt(r2);
	    ck = rpole_ref(1, k);
	    dkx = rpole_ref(2, k);
	    dky = rpole_ref(3, k);
	    dkz = rpole_ref(4, k);
	    qkxx = rpole_ref(5, k);
	    qkxy = rpole_ref(6, k);
	    qkxz = rpole_ref(7, k);
	    qkyy = rpole_ref(9, k);
	    qkyz = rpole_ref(10, k);
	    qkzz = rpole_ref(13, k);
	    scale3 = pscale[kk - 1];
	    scale5 = pscale[kk - 1];
	    scale7 = pscale[kk - 1];
	    damp = pdi * polar_1.pdamp[k - 1];
	    if (damp != 0.) {
/* Computing MIN */
		d__1 = pti, d__2 = polar_1.thole[k - 1];
		pgamma = min(d__1,d__2);
/* Computing 3rd power */
		d__1 = r__ / damp;
		damp = -pgamma * (d__1 * (d__1 * d__1));
		if (damp > -50.) {
		    scale3 *= 1. - exp(damp);
		    scale5 *= 1. - (1. - damp) * exp(damp);
/* Computing 2nd power */
		    d__1 = damp;
		    scale7 *= 1. - (1. - damp + d__1 * d__1 * .6) * exp(damp);
		}
	    }
	    rr3 = scale3 / (r__ * r2);
	    rr5 = scale5 * 3. / (r__ * r2 * r2);
	    rr7 = scale7 * 15. / (r__ * r2 * r2 * r2);
	    dir = dix * xr + diy * yr + diz * zr;
	    qix = qixx * xr + qixy * yr + qixz * zr;
	    qiy = qixy * xr + qiyy * yr + qiyz * zr;
	    qiz = qixz * xr + qiyz * yr + qizz * zr;
	    qir = qix * xr + qiy * yr + qiz * zr;
	    dkr = dkx * xr + dky * yr + dkz * zr;
	    qkx = qkxx * xr + qkxy * yr + qkxz * zr;
	    qky = qkxy * xr + qkyy * yr + qkyz * zr;
	    qkz = qkxz * xr + qkyz * yr + qkzz * zr;
	    qkr = qkx * xr + qky * yr + qkz * zr;
	    fi[0] = -xr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dkx + 
		    rr5 * 2. * qkx;
	    fi[1] = -yr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dky + 
		    rr5 * 2. * qky;
	    fi[2] = -zr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dkz + 
		    rr5 * 2. * qkz;
	    fk[0] = xr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * dix - rr5 
		    * 2. * qix;
	    fk[1] = yr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * diy - rr5 
		    * 2. * qiy;
	    fk[2] = zr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * diz - rr5 
		    * 2. * qiz;
	    for (j = 1; j <= 3; ++j) {
		field_ref(j, i__) = field_ref(j, i__) + fi[j - 1];
		field_ref(j, k) = field_ref(j, k) + fk[j - 1];
	    }
	}
    }

/*     compute induced dipoles as polarizability times field */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    uind_ref(j, i__) = polar_1.polarity[i__ - 1] * field_ref(j, i__);
	}
    }

/*     for direct-only models set mutual scale factors to zero */

    if (s_cmp(polpot_1.poltyp, "DIRECT", (ftnlen)6, (ftnlen)6) == 0) {
	polpot_1.u1scale = 0.;
	polpot_1.u2scale = 0.;
	polpot_1.u3scale = 0.;
	polpot_1.u4scale = 0.;
    }

/*     set tolerances for computation of mutual induced dipoles */

    done = FALSE_;
    maxiter = 500;
    iter = 0;
    eps = 1.;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    udir_ref(j, i__) = uind_ref(j, i__);
	}
    }

/*     compute mutual induced dipole moments by an iterative method */

    while(! done) {
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		field_ref(j, i__) = 0.;
	    }
	}
	i__1 = mpole_1.npole - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    iz = mpole_1.zaxis[ii - 1];
	    ix = mpole_1.xaxis[ii - 1];
	    pdi = polar_1.pdamp[i__ - 1];
	    pti = polar_1.thole[i__ - 1];
	    uix = uind_ref(1, i__);
	    uiy = uind_ref(2, i__);
	    uiz = uind_ref(3, i__);
	    i__2 = mpole_1.npole;
	    for (j = i__ + 1; j <= i__2; ++j) {
		pscale[mpole_1.ipole[j - 1] - 1] = 0.;
	    }
	    i__2 = polgrp_1.np11[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[ip11_ref(j, ii) - 1] = polpot_1.u1scale - 
			polpot_1.d1scale;
	    }
	    i__2 = polgrp_1.np12[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[ip12_ref(j, ii) - 1] = polpot_1.u2scale - 
			polpot_1.d2scale;
	    }
	    i__2 = polgrp_1.np13[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[ip13_ref(j, ii) - 1] = polpot_1.u3scale - 
			polpot_1.d3scale;
	    }
	    i__2 = polgrp_1.np14[ii - 1];
	    for (j = 1; j <= i__2; ++j) {
		pscale[ip14_ref(j, ii) - 1] = polpot_1.u4scale - 
			polpot_1.d4scale;
	    }
	    i__2 = mpole_1.npole;
	    for (k = i__ + 1; k <= i__2; ++k) {
		kk = mpole_1.ipole[k - 1];
		kz = mpole_1.zaxis[kk - 1];
		kx = mpole_1.xaxis[kk - 1];
		xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
		yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
		zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
		r2 = xr * xr + yr * yr + zr * zr;
		r__ = sqrt(r2);
		ukx = uind_ref(1, k);
		uky = uind_ref(2, k);
		ukz = uind_ref(3, k);
		scale3 = pscale[kk - 1];
		scale5 = pscale[kk - 1];
		damp = pdi * polar_1.pdamp[k - 1];
		if (damp != 0.) {
/* Computing MIN */
		    d__1 = pti, d__2 = polar_1.thole[k - 1];
		    pgamma = min(d__1,d__2);
/* Computing 3rd power */
		    d__1 = r__ / damp;
		    damp = -pgamma * (d__1 * (d__1 * d__1));
		    if (damp > -50.) {
			scale3 *= 1. - exp(damp);
			scale5 *= 1. - (1. - damp) * exp(damp);
		    }
		}
		rr3 = scale3 / (r__ * r2);
		rr5 = scale5 * 3. / (r__ * r2 * r2);
		uir = xr * uix + yr * uiy + zr * uiz;
		ukr = xr * ukx + yr * uky + zr * ukz;
		fi[0] = -rr3 * ukx + rr5 * ukr * xr;
		fi[1] = -rr3 * uky + rr5 * ukr * yr;
		fi[2] = -rr3 * ukz + rr5 * ukr * zr;
		fk[0] = -rr3 * uix + rr5 * uir * xr;
		fk[1] = -rr3 * uiy + rr5 * uir * yr;
		fk[2] = -rr3 * uiz + rr5 * uir * zr;
		for (j = 1; j <= 3; ++j) {
		    field_ref(j, i__) = field_ref(j, i__) + fi[j - 1];
		    field_ref(j, k) = field_ref(j, k) + fk[j - 1];
		}
	    }
	}

/*     check to see if the mutual induced dipoles have converged */

	++iter;
	epsold = eps;
	eps = 0.;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		uold_ref(j, i__) = uind_ref(j, i__);
		uind_ref(j, i__) = udir_ref(j, i__) + polar_1.polarity[i__ - 
			1] * field_ref(j, i__);
		uind_ref(j, i__) = uold_ref(j, i__) + polpot_1.polsor * (
			uind_ref(j, i__) - uold_ref(j, i__));
/* Computing 2nd power */
		d__1 = uind_ref(j, i__) - uold_ref(j, i__);
		eps += d__1 * d__1;
	    }
	}
	eps = sqrt(eps / (doublereal) polar_1.npolar) * 4.80321;
	if (iter == 1) {
	    io___266.ciunit = iounit_1.iout;
	    s_wsfe(&io___266);
	    e_wsfe();
	}
	io___267.ciunit = iounit_1.iout;
	s_wsfe(&io___267);
	do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
	e_wsfe();
	if (eps < polpot_1.poleps) {
	    done = TRUE_;
	}
	if (eps > epsold) {
	    done = TRUE_;
	}
	if (iter >= maxiter) {
	    done = TRUE_;
	}
    }

/*     terminate the calculation if dipoles failed to converge */

    if (eps > polpot_1.poleps) {
	io___268.ciunit = iounit_1.iout;
	s_wsfe(&io___268);
	e_wsfe();
	prterr_();
	fatal_();
    }

/*     print out a list of the final induced dipole moments */

    io___269.ciunit = iounit_1.iout;
    s_wsfe(&io___269);
    e_wsfe();
    io___270.ciunit = iounit_1.iout;
    s_wsfe(&io___270);
    e_wsfe();
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = mpole_1.ipole[i__ - 1];
/* Computing 2nd power */
	d__1 = uind_ref(1, i__);
/* Computing 2nd power */
	d__2 = uind_ref(2, i__);
/* Computing 2nd power */
	d__3 = uind_ref(3, i__);
	norm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	io___272.ciunit = iounit_1.iout;
	s_wsfe(&io___272);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	for (j = 1; j <= 3; ++j) {
	    d__1 = uind_ref(j, i__) * 4.80321;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	}
	d__2 = norm * 4.80321;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* intrapol1_ */

#undef rpole_ref
#undef field_ref
#undef uold_ref
#undef udir_ref
#undef uind_ref
#undef ip14_ref
#undef ip13_ref
#undef ip12_ref
#undef ip11_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine fixpolar  --  postprocess multipole moments  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "fixpolar" averages multipoles over equivalent sites, sets */
/*     symmetric components at achiral atoms to zero, and maintains */
/*     an integer net charge and traceless quadrupoles */


/* Subroutine */ int fixpolar_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Average the Multipole Moments of Equival"
	    "ent\002,\002 Atoms [N] :  \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_40[] = "(/,\002 Remove Multipole Components Zeroed by"
	    "\002,\002 Symmetry [N] :  \002,$)";
    static char fmt_50[] = "(a120)";
    static char fmt_70[] = "(/,\002 Final Multipole Moments for the AMOEBA F"
	    "orce\002,\002 Field :\002)";
    static char fmt_80[] = "(/,\002 Atom:\002,i8,9x,\002Name:\002,3x,a3,7x"
	    ",\002Atomic Number:\002,i8)";
    static char fmt_90[] = "(/,\002 No Atomic Multipole Moments for this S"
	    "ite\002)";
    static char fmt_100[] = "(/,\002 Atom:\002,i8,9x,\002Name:\002,3x,a3,7x"
	    ",\002Atomic Number:\002,i8)";
    static char fmt_110[] = "(/,\002 Local Frame:\002,12x,a8,6x,3i8)";
    static char fmt_120[] = "(/,\002 Charge:\002,10x,f15.5)";
    static char fmt_130[] = "(\002 Dipole:\002,10x,3f15.5)";
    static char fmt_140[] = "(\002 Quadrupole:\002,6x,f15.5)";
    static char fmt_150[] = "(18x,2f15.5)";
    static char fmt_160[] = "(18x,3f15.5)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void), s_cmp(char *, char *, ftnlen, 
	    ftnlen), i_dnnt(doublereal *);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal ci;
    static integer ii, kk;
    static doublereal cj;
    static integer nx;
    static doublereal big, eps, sum, pave[13];
    static integer list[8], next, ixaxe, iyaxe, izaxe, nlist;
    static logical yzero;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char answer[1];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static cilist io___273 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___274 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___288 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___289 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___296 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___297 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___298 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___302 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___303 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___304 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___305 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___306 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___307 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___308 = { 0, 0, 0, fmt_160, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define pole_ref(a_1,a_2) mpole_1.pole[(a_2)*13 + a_1 - 14]
#define polaxe_ref(a_0,a_1) &mpole_1.polaxe[(a_1)*8 + a_0 - 8]



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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




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




/*     optionally average multipoles for equivalent atoms */

    io___273.ciunit = iounit_1.iout;
    s_wsfe(&io___273);
    e_wsfe();
    io___274.ciunit = iounit_1.input;
    s_rsfe(&io___274);
    do_fio(&c__1, record, (ftnlen)120);
    e_rsfe();
    next = 1;
    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    upcase_(answer, (ftnlen)1);

/*     perform averaging for equivalent monovalent atoms */

    if (*(unsigned char *)answer == 'Y') {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    nlist = 0;
	    i__2 = couple_1.n12[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		k = i12_ref(j, i__);
		if (couple_1.n12[k - 1] == 1) {
		    i__3 = nlist;
		    for (m = 1; m <= i__3; ++m) {
			if (list[m - 1] == atmtyp_1.atomic[k - 1]) {
			    goto L30;
			}
		    }
		    ++nlist;
		    list[nlist - 1] = atmtyp_1.atomic[k - 1];
L30:
		    ;
		}
	    }
	    i__2 = nlist;
	    for (ii = 1; ii <= i__2; ++ii) {
		kk = list[ii - 1];
		for (j = 1; j <= 13; ++j) {
		    pave[j - 1] = 0.;
		}
		nx = 0;
		i__3 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__3; ++j) {
		    k = i12_ref(j, i__);
		    if (mpole_1.pollist[k - 1] != 0) {
			if (atmtyp_1.atomic[k - 1] == kk && couple_1.n12[k - 
				1] == 1) {
			    ++nx;
			    for (m = 1; m <= 13; ++m) {
				pave[m - 1] += pole_ref(m, k);
			    }
			}
		    }
		}
		if (nx >= 2) {
		    for (j = 1; j <= 13; ++j) {
			pave[j - 1] /= (doublereal) nx;
		    }
		    i__3 = couple_1.n12[i__ - 1];
		    for (j = 1; j <= i__3; ++j) {
			k = i12_ref(j, i__);
			if (mpole_1.pollist[k - 1] != 0) {
			    if (atmtyp_1.atomic[k - 1] == kk && couple_1.n12[
				    k - 1] == 1) {
				for (m = 1; m <= 13; ++m) {
				    pole_ref(m, k) = pave[m - 1];
				}
			    }
			}
		    }
		}
	    }
	}
    }

/*     optionally set symmetric multipole components to zero */

    io___288.ciunit = iounit_1.iout;
    s_wsfe(&io___288);
    e_wsfe();
    io___289.ciunit = iounit_1.input;
    s_rsfe(&io___289);
    do_fio(&c__1, record, (ftnlen)120);
    e_rsfe();
    next = 1;
    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    upcase_(answer, (ftnlen)1);

/*     remove multipole components that are zero by symmetry */

    if (*(unsigned char *)answer == 'Y') {
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    yzero = FALSE_;
	    if (mpole_1.yaxis[i__ - 1] == 0) {
		yzero = TRUE_;
	    }
	    if (s_cmp(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8) ==
		     0) {
		yzero = TRUE_;
	    }
	    if (s_cmp(polaxe_ref(0, i__), "Z-Bisect", (ftnlen)8, (ftnlen)8) ==
		     0) {
		yzero = TRUE_;
	    }
	    if (mpole_1.zaxis[i__ - 1] == 0 || mpole_1.zaxis[i__ - 1] > 
		    atoms_1.n) {
		pole_ref(13, i__) = 0.;
	    }
	    if (mpole_1.xaxis[i__ - 1] == 0 || mpole_1.xaxis[i__ - 1] > 
		    atoms_1.n) {
		pole_ref(2, i__) = 0.;
		pole_ref(5, i__) = pole_ref(13, i__) * -.5;
		pole_ref(7, i__) = 0.;
		pole_ref(9, i__) = pole_ref(5, i__);
		pole_ref(11, i__) = 0.;
	    }
	    if (yzero) {
		pole_ref(3, i__) = 0.;
		pole_ref(6, i__) = 0.;
		pole_ref(8, i__) = 0.;
		pole_ref(10, i__) = 0.;
		pole_ref(12, i__) = 0.;
	    }
	}
    }

/*     convert dipole and quadrupole moments back to atomic units */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pole_ref(1, i__) = pole_ref(1, i__);
	for (j = 2; j <= 4; ++j) {
	    pole_ref(j, i__) = pole_ref(j, i__) / .52917720859;
	}
	for (j = 5; j <= 13; ++j) {
	    pole_ref(j, i__) = pole_ref(j, i__) * 3. / .28002851809110429;
	}
    }

/*     regularize the multipole moments to desired precision */

    eps = 1e-5;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 13; ++j) {
	    d__1 = pole_ref(j, i__) / eps;
	    pole_ref(j, i__) = (doublereal) i_dnnt(&d__1) * eps;
	}
    }

/*     maintain integer net charge for the whole system */

    k = 0;
    big = 0.;
    sum = 0.;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum += pole_ref(1, i__);
	ci = (d__1 = pole_ref(1, i__), abs(d__1));
	if (ci > big) {
	    i__2 = atoms_1.n;
	    for (j = 1; j <= i__2; ++j) {
		cj = (d__1 = pole_ref(1, j), abs(d__1));
		if (i__ != j && ci == cj) {
		    goto L60;
		}
	    }
	    k = i__;
	    big = ci;
L60:
	    ;
	}
    }
    sum -= (doublereal) i_dnnt(&sum);
    if (k != 0) {
	pole_ref(1, k) = pole_ref(1, k) - sum;
    }

/*     maintain traceless quadrupole at each multipole site */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = pole_ref(5, i__) + pole_ref(9, i__) + pole_ref(13, i__);
/* Computing MAX */
	d__4 = (d__1 = pole_ref(5, i__), abs(d__1)), d__5 = (d__2 = pole_ref(
		9, i__), abs(d__2)), d__4 = max(d__4,d__5), d__5 = (d__3 = 
		pole_ref(13, i__), abs(d__3));
	big = max(d__4,d__5);
	k = 0;
	if (big == (d__1 = pole_ref(5, i__), abs(d__1))) {
	    k = 5;
	}
	if (big == (d__1 = pole_ref(9, i__), abs(d__1))) {
	    k = 9;
	}
	if (big == (d__1 = pole_ref(13, i__), abs(d__1))) {
	    k = 13;
	}
	if (k != 0) {
	    pole_ref(k, i__) = pole_ref(k, i__) - sum;
	}
    }

/*     print the final post-processed multipoles for AMOEBA */

    io___296.ciunit = iounit_1.iout;
    s_wsfe(&io___296);
    e_wsfe();
    i__1 = atoms_1.n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = mpole_1.pollist[ii - 1];
	if (i__ == 0) {
	    io___297.ciunit = iounit_1.iout;
	    s_wsfe(&io___297);
	    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
	    do_fio(&c__1, (char *)&atmtyp_1.atomic[ii - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	    io___298.ciunit = iounit_1.iout;
	    s_wsfe(&io___298);
	    e_wsfe();
	} else {
	    izaxe = mpole_1.zaxis[i__ - 1];
	    ixaxe = mpole_1.xaxis[i__ - 1];
	    iyaxe = mpole_1.yaxis[i__ - 1];
	    if (izaxe > atoms_1.n) {
		izaxe = 0;
	    }
	    if (ixaxe > atoms_1.n) {
		ixaxe = 0;
	    }
	    if (iyaxe < 0) {
		iyaxe = -iyaxe;
	    }
	    io___302.ciunit = iounit_1.iout;
	    s_wsfe(&io___302);
	    do_fio(&c__1, (char *)&ii, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, ii), (ftnlen)3);
	    do_fio(&c__1, (char *)&atmtyp_1.atomic[ii - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	    io___303.ciunit = iounit_1.iout;
	    s_wsfe(&io___303);
	    do_fio(&c__1, polaxe_ref(0, i__), (ftnlen)8);
	    do_fio(&c__1, (char *)&izaxe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ixaxe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iyaxe, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___304.ciunit = iounit_1.iout;
	    s_wsfe(&io___304);
	    do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___305.ciunit = iounit_1.iout;
	    s_wsfe(&io___305);
	    do_fio(&c__1, (char *)&pole_ref(2, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(3, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(4, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___306.ciunit = iounit_1.iout;
	    s_wsfe(&io___306);
	    do_fio(&c__1, (char *)&pole_ref(5, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___307.ciunit = iounit_1.iout;
	    s_wsfe(&io___307);
	    do_fio(&c__1, (char *)&pole_ref(8, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(9, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	    io___308.ciunit = iounit_1.iout;
	    s_wsfe(&io___308);
	    do_fio(&c__1, (char *)&pole_ref(11, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(12, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pole_ref(13, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* fixpolar_ */

#undef polaxe_ref
#undef pole_ref
#undef name___ref
#undef i12_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine prtpolar  --  create file with final multipoles  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "prtpolar" makes a key file containing results from distributed */
/*     multipole analysis or removal of intramolecular polarization */


/* Subroutine */ int prtpolar_(void)
{
    /* Format strings */
    static char fmt_10[] = "(a)";
    static char fmt_20[] = "()";
    static char fmt_30[] = "(\002atom\002,6x,2i5,4x,a3,4x,\002\"\002,a20"
	    ",\002\"\002,i10,f10.3,i5)";
    static char fmt_40[] = "()";
    static char fmt_50[] = "(\002multipole\002,1x,3i5,11x,f11.5)";
    static char fmt_60[] = "(\002multipole\002,1x,4i5,6x,f11.5)";
    static char fmt_70[] = "(\002multipole\002,1x,3i5,11x,f11.5)";
    static char fmt_80[] = "(\002multipole\002,1x,4i5,6x,f11.5)";
    static char fmt_90[] = "(\002multipole\002,1x,4i5,6x,f11.5)";
    static char fmt_100[] = "(\002multipole\002,1x,4i5,6x,f11.5)";
    static char fmt_110[] = "(36x,3f11.5)";
    static char fmt_120[] = "(36x,f11.5)";
    static char fmt_130[] = "(36x,2f11.5)";
    static char fmt_140[] = "(36x,3f11.5)";
    static char fmt_150[] = "()";
    static char fmt_160[] = "(\002polarize\002,2x,i5,9x,2f11.3,3x,20i5)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3, i__4, i__5;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_cmp(char *, char *, ftnlen, ftnlen), f_clos(cllist *);

    /* Local variables */
    extern integer freeunit_(void), trimtext_(char *, ftnlen);
    static integer i__, j, k, it, ikey, size, ixaxe, iyaxe, izaxe;
    static char record[120];
    static logical dofull;
    static char keyfile[120];
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___314 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___315 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___317 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___318 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___323 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___324 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___325 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___326 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___327 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___328 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___329 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___330 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___331 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___332 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___333 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___336 = { 0, 0, 0, fmt_160, 0 };



#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define pole_ref(a_1,a_2) mpole_1.pole[(a_2)*13 + a_1 - 14]
#define pgrp_ref(a_1,a_2) kpolr_1.pgrp[(a_2)*8 + a_1 - 9]
#define story_ref(a_0,a_1) &atmtyp_1.story[(a_1)*24 + a_0 - 24]
#define polaxe_ref(a_0,a_1) &mpole_1.polaxe[(a_1)*8 + a_0 - 8]
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
/*     ##  kpolr.i  --  forcefield parameters for polarizability  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     polr   dipole polarizability parameters for each atom type */
/*     athl   Thole polarizability damping value for each atom type */
/*     pgrp   connected types in polarization group of each atom type */




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




/*     output some definitions and parameters to a keyfile */

    ikey = freeunit_();
/* Writing concatenation */
    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
    i__1[1] = 4, a__1[1] = ".key";
    s_cat(keyfile, a__1, i__1, &c__2, (ftnlen)120);
    version_(keyfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = ikey;
    o__1.ofnmlen = 120;
    o__1.ofnm = keyfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);

/*     copy the contents of any previously existing keyfile */

    i__2 = keys_1.nkey;
    for (i__ = 1; i__ <= i__2; ++i__) {
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	size = trimtext_(record, (ftnlen)120);
	io___314.ciunit = ikey;
	s_wsfe(&io___314);
	do_fio(&c__1, record, size);
	e_wsfe();
    }
    if (keys_1.nkey != 0) {
	io___315.ciunit = ikey;
	s_wsfe(&io___315);
	e_wsfe();
    }

/*     output the atom definitions to the keyfile as appropriate */

    dofull = TRUE_;
    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (atoms_1.type__[i__ - 1] != i__) {
	    dofull = FALSE_;
	}
    }
    if (dofull) {
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    io___317.ciunit = ikey;
	    s_wsfe(&io___317);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
	    do_fio(&c__1, story_ref(0, i__), (ftnlen)24);
	    do_fio(&c__1, (char *)&atmtyp_1.atomic[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&atmtyp_1.mass[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&atmtyp_1.valence[i__ - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	}
	if (atoms_1.n != 0) {
	    io___318.ciunit = ikey;
	    s_wsfe(&io___318);
	    e_wsfe();
	}
    }

/*     output the local frame multipole values to the keyfile */

    i__2 = mpole_1.npole;
    for (i__ = 1; i__ <= i__2; ++i__) {
	it = mpole_1.ipole[i__ - 1];
	if (! dofull) {
	    it = -it;
	}
	izaxe = mpole_1.zaxis[i__ - 1];
	ixaxe = mpole_1.xaxis[i__ - 1];
	iyaxe = mpole_1.yaxis[i__ - 1];
	if (izaxe > atoms_1.n) {
	    izaxe = 0;
	}
	if (ixaxe > atoms_1.n) {
	    ixaxe = 0;
	}
	if (iyaxe < 0) {
	    iyaxe = -iyaxe;
	}
	if (s_cmp(polaxe_ref(0, i__), "Z-then-X", (ftnlen)8, (ftnlen)8) == 0) 
		{
	    if (mpole_1.yaxis[i__ - 1] == 0) {
		io___323.ciunit = ikey;
		s_wsfe(&io___323);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&izaxe, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ixaxe, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___324.ciunit = ikey;
		s_wsfe(&io___324);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&izaxe, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ixaxe, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iyaxe, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(polaxe_ref(0, i__), "Bisector", (ftnlen)8, (ftnlen)8)
		 == 0) {
	    if (mpole_1.yaxis[i__ - 1] == 0) {
		io___325.ciunit = ikey;
		s_wsfe(&io___325);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		i__3 = -izaxe;
		do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		i__4 = -ixaxe;
		do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___326.ciunit = ikey;
		s_wsfe(&io___326);
		do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
		i__3 = -izaxe;
		do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		i__4 = -ixaxe;
		do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iyaxe, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(polaxe_ref(0, i__), "Z-Bisect", (ftnlen)8, (ftnlen)8)
		 == 0) {
	    io___327.ciunit = ikey;
	    s_wsfe(&io___327);
	    do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&izaxe, (ftnlen)sizeof(integer));
	    i__3 = -ixaxe;
	    do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
	    i__4 = -iyaxe;
	    do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	} else if (s_cmp(polaxe_ref(0, i__), "3-Fold", (ftnlen)8, (ftnlen)6) 
		== 0) {
	    io___328.ciunit = ikey;
	    s_wsfe(&io___328);
	    do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
	    i__3 = -izaxe;
	    do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
	    i__4 = -ixaxe;
	    do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
	    i__5 = -iyaxe;
	    do_fio(&c__1, (char *)&i__5, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&pole_ref(1, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
	io___329.ciunit = ikey;
	s_wsfe(&io___329);
	do_fio(&c__1, (char *)&pole_ref(2, i__), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&pole_ref(3, i__), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&pole_ref(4, i__), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___330.ciunit = ikey;
	s_wsfe(&io___330);
	do_fio(&c__1, (char *)&pole_ref(5, i__), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___331.ciunit = ikey;
	s_wsfe(&io___331);
	do_fio(&c__1, (char *)&pole_ref(8, i__), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&pole_ref(9, i__), (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___332.ciunit = ikey;
	s_wsfe(&io___332);
	do_fio(&c__1, (char *)&pole_ref(11, i__), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&pole_ref(12, i__), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&pole_ref(13, i__), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     output the polarizability parameters to the keyfile */

    if (dofull) {
	if (atoms_1.n != 0) {
	    io___333.ciunit = ikey;
	    s_wsfe(&io___333);
	    e_wsfe();
	}
	i__2 = mpole_1.npole;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    k = 0;
	    for (j = 1; j <= 8; ++j) {
		if (pgrp_ref(j, i__) != 0) {
		    k = j;
		}
	    }
	    io___336.ciunit = ikey;
	    s_wsfe(&io___336);
	    do_fio(&c__1, (char *)&mpole_1.ipole[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&polar_1.polarity[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&polar_1.thole[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    i__3 = k;
	    for (j = 1; j <= i__3; ++j) {
		do_fio(&c__1, (char *)&pgrp_ref(j, i__), (ftnlen)sizeof(
			integer));
	    }
	    e_wsfe();
	}
    }
    cl__1.cerr = 0;
    cl__1.cunit = ikey;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* prtpolar_ */

#undef keyline_ref
#undef polaxe_ref
#undef story_ref
#undef pgrp_ref
#undef pole_ref
#undef name___ref


/* Main program alias */ int poledit_ () { MAIN__ (); return 0; }
