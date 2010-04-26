/* polarize.f -- translated by f2c (version 20050501).
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
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

struct {
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

struct {
    doublereal polarity[25000], thole[25000], pdamp[25000], uind[75000]	/* 
	    was [3][25000] */, uinp[75000]	/* was [3][25000] */, uinds[
	    75000]	/* was [3][25000] */, uinps[75000]	/* was [3][
	    25000] */;
    integer npolar;
} polar_;

#define polar_1 polar_

struct {
    doublereal poleps, polsor, p2scale, p3scale, p4scale, p5scale, d1scale, 
	    d2scale, d3scale, d4scale, u1scale, u2scale, u3scale, u4scale;
    char poltyp[6];
} polpot_;

#define polpot_1 polpot_

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
    integer np11[25000], ip11[2500000]	/* was [100][25000] */, np12[25000], 
	    ip12[1250000]	/* was [50][25000] */, np13[25000], ip13[
	    1250000]	/* was [50][25000] */, np14[25000], ip14[1250000]	
	    /* was [50][25000] */;
} polgrp_;

#define polgrp_1 polgrp_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;



/*     ############################################################# */
/*     ##  COPYRIGHT (C) 2001 by Pengyu Ren & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  program polarize  --  compute the molecular polarizability  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "polarize" computes the molecular polarizability by applying */
/*     an external field along each axis followed by diagonalization */
/*     of the resulting polarizability tensor */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 POLARIZE  --  Dipole Polarizability\002"
	    ",\002 is Not in Use\002)";
    static char fmt_20[] = "(/,\002 Additive Molecular Polarizability :\002,"
	    "f15.4)";
    static char fmt_30[] = "(/,\002 Additive Total Polarizability :\002,4x,f"
	    "15.4)";
    static char fmt_40[] = "(/,\002 Molecular Polarizability Tensor :\002,/)";
    static char fmt_50[] = "(/,\002 Total Polarizability Tensor:\002,/)";
    static char fmt_60[] = "(15x,3f12.4,/,15x,3f12.4,/,15x,3f12.4)";
    static char fmt_70[] = "(/,\002 Polarizability Tensor Eigenvalues :\002,"
	    "/)";
    static char fmt_80[] = "(15x,3f12.4)";
    static char fmt_90[] = "(/,\002 Interactive Molecular Polarizability "
	    ":\002,f12.4)";
    static char fmt_100[] = "(/,\002 Interactive Total Polarizability :\002,"
	    "4x,f12.4)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int molecule_(void);
    static integer i__;
    static doublereal addu, umol[3], work1[3], work2[3];
    extern /* Subroutine */ int field_(void);
    static doublereal alpha[9]	/* was [3][3] */;
    extern /* Subroutine */ int fatal_(void), jacobi_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal dalpha[3], malpha, valpha[9]	/* was [3][3] */;
    extern /* Subroutine */ int kpolar_(void), getxyz_(void);
    static doublereal exfield[3];
    extern /* Subroutine */ int initial_(void), moluind_(doublereal *, 
	    doublereal *);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_100, 0 };



#define alpha_ref(a_1,a_2) alpha[(a_2)*3 + a_1 - 4]



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




/*     get the coordinates and required force field parameters */

    initial_();
    getxyz_();
    field_();
    molecule_();
    kpolar_();

/*     sum atomic polarizabilities to get additive molecular value */

    if (! potent_1.use_polar__) {
	io___1.ciunit = iounit_1.iout;
	s_wsfe(&io___1);
	e_wsfe();
	fatal_();
    }
    addu = 0.;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	addu = polar_1.polarity[i__ - 1] + addu;
    }
    if (molcul_1.nmol == 1) {
	io___4.ciunit = iounit_1.iout;
	s_wsfe(&io___4);
	do_fio(&c__1, (char *)&addu, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else {
	io___5.ciunit = iounit_1.iout;
	s_wsfe(&io___5);
	do_fio(&c__1, (char *)&addu, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     compute each column of the polarizability tensor */

    for (i__ = 1; i__ <= 3; ++i__) {
	exfield[i__ - 1] = 0.;
    }
    exfield[0] = .001;
    moluind_(exfield, umol);
    alpha_ref(1, 1) = umol[0] / exfield[0];
    alpha_ref(2, 1) = umol[1] / exfield[0];
    alpha_ref(3, 1) = umol[2] / exfield[0];
    for (i__ = 1; i__ <= 3; ++i__) {
	exfield[i__ - 1] = 0.;
    }
    exfield[1] = .001;
    moluind_(exfield, umol);
    alpha_ref(1, 2) = umol[0] / exfield[1];
    alpha_ref(2, 2) = umol[1] / exfield[1];
    alpha_ref(3, 2) = umol[2] / exfield[1];
    for (i__ = 1; i__ <= 3; ++i__) {
	exfield[i__ - 1] = 0.;
    }
    exfield[2] = .001;
    moluind_(exfield, umol);
    alpha_ref(1, 3) = umol[0] / exfield[2];
    alpha_ref(2, 3) = umol[1] / exfield[2];
    alpha_ref(3, 3) = umol[2] / exfield[2];

/*     print out the full polarizability tensor */

    if (molcul_1.nmol == 1) {
	io___9.ciunit = iounit_1.iout;
	s_wsfe(&io___9);
	e_wsfe();
    } else {
	io___10.ciunit = iounit_1.iout;
	s_wsfe(&io___10);
	e_wsfe();
    }
    io___11.ciunit = iounit_1.iout;
    s_wsfe(&io___11);
    do_fio(&c__1, (char *)&alpha_ref(1, 1), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&alpha_ref(1, 2), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&alpha_ref(1, 3), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&alpha_ref(2, 1), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&alpha_ref(2, 2), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&alpha_ref(2, 3), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&alpha_ref(3, 1), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&alpha_ref(3, 2), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&alpha_ref(3, 3), (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     diagonalize the tensor and get molecular polarizability */

    jacobi_(&c__3, &c__3, alpha, dalpha, valpha, work1, work2);
    io___16.ciunit = iounit_1.iout;
    s_wsfe(&io___16);
    e_wsfe();
    io___17.ciunit = iounit_1.iout;
    s_wsfe(&io___17);
    do_fio(&c__1, (char *)&dalpha[0], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&dalpha[1], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&dalpha[2], (ftnlen)sizeof(doublereal));
    e_wsfe();
    malpha = (dalpha[0] + dalpha[1] + dalpha[2]) / 3.;
    if (molcul_1.nmol == 1) {
	io___19.ciunit = iounit_1.iout;
	s_wsfe(&io___19);
	do_fio(&c__1, (char *)&malpha, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else {
	io___20.ciunit = iounit_1.iout;
	s_wsfe(&io___20);
	do_fio(&c__1, (char *)&malpha, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* MAIN__ */

#undef alpha_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine moluind  --  molecular induced dipole in field  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "moluind" computes the molecular induced dipole components */
/*     in the presence of an external electric field */


/* Subroutine */ int moluind_(doublereal *exfield, doublereal *umol)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Determination of Induced Dipole\002,\002"
	    " Moments :\002,//,4x,\002Iter\002,8x,\002RMS Change (Debyes)\002"
	    ",/)";
    static char fmt_20[] = "(i8,7x,f16.10)";
    static char fmt_30[] = "(/,\002 MOLUIND  --  Warning, Induced Dipoles"
	    "\002,\002 are not Converged\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j;
    static doublereal eps;
    static logical done;
    static integer iter;
    static doublereal udir[75000]	/* was [3][25000] */, uold[75000]	
	    /* was [3][25000] */, field[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int ufield_(doublereal *);
    static doublereal epsold;
    static integer maxiter;

    /* Fortran I/O blocks */
    static cilist io___31 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_30, 0 };



#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define udir_ref(a_1,a_2) udir[(a_2)*3 + a_1 - 4]
#define uold_ref(a_1,a_2) uold[(a_2)*3 + a_1 - 4]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1 - 4]



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




/*     set induced dipoles to polarizability times external field */

    /* Parameter adjustments */
    --umol;
    --exfield;

    /* Function Body */
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    uind_ref(j, i__) = polar_1.polarity[i__ - 1] * exfield[j];
	}
    }

/*     compute direct induced dipole moments from direct field */

    if (s_cmp(polpot_1.poltyp, "DIRECT", (ftnlen)6, (ftnlen)6) == 0) {
	ufield_(field);
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		uind_ref(j, i__) = uind_ref(j, i__) + polar_1.polarity[i__ - 
			1] * field_ref(j, i__);
	    }
	}
    }

/*     compute mutual induced dipole moments by an iterative method */

    if (s_cmp(polpot_1.poltyp, "MUTUAL", (ftnlen)6, (ftnlen)6) == 0) {
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

/*     check to see if the mutual induced dipoles have converged */

	while(! done) {
	    ufield_(field);
	    ++iter;
	    epsold = eps;
	    eps = 0.;
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    uold_ref(j, i__) = uind_ref(j, i__);
		    uind_ref(j, i__) = udir_ref(j, i__) + polar_1.polarity[
			    i__ - 1] * field_ref(j, i__);
		    uind_ref(j, i__) = uold_ref(j, i__) + polpot_1.polsor * (
			    uind_ref(j, i__) - uold_ref(j, i__));
/* Computing 2nd power */
		    d__1 = uind_ref(j, i__) - uold_ref(j, i__);
		    eps += d__1 * d__1;
		}
	    }
	    eps = sqrt(eps / (doublereal) polar_1.npolar) * 4.80321;
	    if (inform_1.debug) {
		if (iter == 1) {
		    io___31.ciunit = iounit_1.iout;
		    s_wsfe(&io___31);
		    e_wsfe();
		}
		io___32.ciunit = iounit_1.iout;
		s_wsfe(&io___32);
		do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
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

/*     print a warning if induced dipoles failed to converge */

	if (eps > polpot_1.poleps) {
	    io___33.ciunit = iounit_1.iout;
	    s_wsfe(&io___33);
	    e_wsfe();
	}
    }

/*     sum up the total molecular induced dipole components */

    for (j = 1; j <= 3; ++j) {
	umol[j] = 0.;
    }
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	umol[1] += uind_ref(1, i__);
	umol[2] += uind_ref(2, i__);
	umol[3] += uind_ref(3, i__);
    }
    return 0;
} /* moluind_ */

#undef field_ref
#undef uold_ref
#undef udir_ref
#undef uind_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine ufield  --  electric field from induced dipoles  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "ufield" finds the field at each polarizable site due to the */
/*     induced dipoles at the other sites using Thole's method to */
/*     damp the field at close range */


/* Subroutine */ int ufield_(doublereal *field)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal r__, r2, fi[3], fk[3];
    static integer ii, kk;
    static doublereal xr, yr, zr, rr3, rr5, pdi, pti, uir, ukr, uix, uiy, uiz,
	     ukx, uky, ukz, damp, scale3, scale5, pgamma, pscale[25000];


#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define ip12_ref(a_1,a_2) polgrp_1.ip12[(a_2)*50 + a_1 - 51]
#define ip13_ref(a_1,a_2) polgrp_1.ip13[(a_2)*50 + a_1 - 51]
#define ip14_ref(a_1,a_2) polgrp_1.ip14[(a_2)*50 + a_1 - 51]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define field_ref(a_1,a_2) field[(a_2)*3 + a_1]



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




/*     zero out the value of the electric field at each site */

    /* Parameter adjustments */
    field -= 4;

    /* Function Body */
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    field_ref(j, i__) = 0.;
	}
    }

/*     loop over pairs of sites incrementing the electric field */

    i__1 = mpole_1.npole - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	pdi = polar_1.pdamp[i__ - 1];
	pti = polar_1.thole[i__ - 1];
	uix = uind_ref(1, i__);
	uiy = uind_ref(2, i__);
	uiz = uind_ref(3, i__);
	i__2 = mpole_1.npole;
	for (j = i__ + 1; j <= i__2; ++j) {
	    pscale[mpole_1.ipole[j - 1] - 1] = 1.;
	}
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[ip11_ref(j, ii) - 1] = polpot_1.u1scale;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[ip12_ref(j, ii) - 1] = polpot_1.u2scale;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[ip13_ref(j, ii) - 1] = polpot_1.u3scale;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[ip14_ref(j, ii) - 1] = polpot_1.u4scale;
	}
	i__2 = mpole_1.npole;
	for (k = i__ + 1; k <= i__2; ++k) {
	    kk = mpole_1.ipole[k - 1];
	    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
	    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
	    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
	    r2 = xr * xr + yr * yr + zr * zr;
	    ukx = uind_ref(1, k);
	    uky = uind_ref(2, k);
	    ukz = uind_ref(3, k);
	    r__ = sqrt(r2);
	    rr3 = 1. / (r__ * r2);
	    rr5 = rr3 * 3. / r2;
	    uir = xr * uix + yr * uiy + zr * uiz;
	    ukr = xr * ukx + yr * uky + zr * ukz;

/*     adjust the field to account for polarization damping */

	    scale3 = 1.;
	    scale5 = 1.;
	    damp = pdi * polar_1.pdamp[k - 1];
	    if (damp != 0.) {
/* Computing MIN */
		d__1 = pti, d__2 = polar_1.thole[k - 1];
		pgamma = min(d__1,d__2);
/* Computing 3rd power */
		d__1 = r__ / damp;
		damp = -pgamma * (d__1 * (d__1 * d__1));
		if (damp > -50.) {
		    scale3 = 1. - exp(damp);
		    scale5 = 1. - (1. - damp) * exp(damp);
		}
	    }
	    fi[0] = -rr3 * ukx * scale3 + rr5 * ukr * xr * scale5;
	    fi[1] = -rr3 * uky * scale3 + rr5 * ukr * yr * scale5;
	    fi[2] = -rr3 * ukz * scale3 + rr5 * ukr * zr * scale5;
	    fk[0] = -rr3 * uix * scale3 + rr5 * uir * xr * scale5;
	    fk[1] = -rr3 * uiy * scale3 + rr5 * uir * yr * scale5;
	    fk[2] = -rr3 * uiz * scale3 + rr5 * uir * zr * scale5;
	    for (j = 1; j <= 3; ++j) {
		fi[j - 1] *= pscale[kk - 1];
		fk[j - 1] *= pscale[kk - 1];
	    }

/*     increment the field at each site due to this interaction */

	    for (j = 1; j <= 3; ++j) {
		field_ref(j, i__) = field_ref(j, i__) + fi[j - 1];
		field_ref(j, k) = field_ref(j, k) + fk[j - 1];
	    }
	}
    }
    return 0;
} /* ufield_ */

#undef field_ref
#undef uind_ref
#undef ip14_ref
#undef ip13_ref
#undef ip12_ref
#undef ip11_ref


/* Main program alias */ int polarize_ () { MAIN__ (); return 0; }
