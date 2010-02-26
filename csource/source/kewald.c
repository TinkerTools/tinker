/* kewald.f -- translated by f2c (version 20050501).
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
    doublereal aewald;
    char boundary[7];
} ewald_;

#define ewald_1 ewald_

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
    doublereal bsmod1[100], bsmod2[100], bsmod3[100], table[1200]	/* 
	    was [400][3] */, qgrid[2000000]	/* was [2][100][100][100] */, 
	    qfac[1000000]	/* was [100][100][100] */, thetai1[1000000]	
	    /* was [4][10][25000] */, thetai2[1000000]	/* was [4][10][25000] 
	    */, thetai3[1000000]	/* was [4][10][25000] */;
    integer nfft1, nfft2, nfft3, bsorder, iprime[45]	/* was [15][3] */, 
	    igrid[75000]	/* was [3][25000] */;
} pme_;

#define pme_1 pme_

struct {
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__3 = 3;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine kewald  --  Ewald sum parameter assignment  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "kewald" assigns particle mesh Ewald parameters and options */
/*     for a periodic system */


/* Subroutine */ int kewald_(void)
{
    /* Initialized data */

    static integer multi[54] = { 2,4,6,8,10,12,16,18,20,24,30,32,36,40,48,50,
	    54,60,64,72,80,90,96,100,108,120,128,144,150,160,162,180,192,200,
	    216,240,250,256,270,288,300,320,324,360,384,400,432,450,480,486,
	    500,512,540,576 };

    /* Format strings */
    static char fmt_30[] = "(/,\002 KEWALD  --  B-Spline Order Too Large;"
	    "\002,\002 Increase MAXORDER\002)";
    static char fmt_40[] = "(/,\002 KEWALD  --  FFT Charge Grid Too Large"
	    ";\002,\002 Increase MAXFFT\002)";
    static char fmt_50[] = "(/,\002 KEWALD  --  Warning, Small Charge Gri"
	    "d\002,\002 may give Poor Accuracy\002)";
    static char fmt_60[] = "(/,\002 Smooth Particle Mesh Ewald Parameters "
	    ":\002,//,4x,\002Ewald Coefficient\002,6x,\002Charge Grid\002,"
	    "\002 Dimensions\002,6x,\002B-Spline Order\002,//,8x,f8.4,11x,3i6"
	    ",12x,i6)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int ewaldcof_(doublereal *, doublereal *), 
	    fftsetup_(void);
    static integer i__, k;
    static doublereal dens, rmax;
    static integer next, ifft1, ifft2, ifft3;
    extern /* Subroutine */ int fatal_(void);
    static doublereal delta;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), moduli_(void);
    static char string[120];
    extern /* Subroutine */ int extent_(doublereal *), lattice_(void);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___13 = { 1, string, 0, 0, 120, 1 };
    static icilist io___14 = { 1, string, 1, 0, 120, 1 };
    static icilist io___15 = { 1, string, 0, 0, 120, 1 };
    static cilist io___17 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_60, 0 };



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
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  ewald.i  --  parameters and options for Ewald summation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     aewald     Ewald convergence coefficient value (Ang-1) */
/*     boundary   Ewald boundary condition; none, tinfoil or vacuum */




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
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  pme.i  --  values for particle mesh Ewald summation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxorder   maximum order of the B-spline approximation */
/*     maxprime   maximum number of prime factors of FFT dimension */
/*     maxtable   maximum size of the FFT table array */

/*     bsmod1     B-spline moduli along the a-axis direction */
/*     bsmod2     B-spline moduli along the b-axis direction */
/*     bsmod3     B-spline moduli along the c-axis direction */
/*     table      intermediate array used by the FFT calculation */
/*     qgrid      values on the particle mesh Ewald charge grid */
/*     qfac       prefactors for particle mesh Ewald charge grid */
/*     thetai1    B-spline coefficients along the a-axis */
/*     thetai2    B-spline coefficients along the b-axis */
/*     thetai3    B-spline coefficients along the c-axis */
/*     nfft1      number of grid points along the a-axis direction */
/*     nfft2      number of grid points along the b-axis direction */
/*     nfft3      number of grid points along the c-axis direction */
/*     bsorder    order of the PME B-spline approximation */
/*     iprime     prime factorization of each FFT dimension */
/*     igrid      initial Ewald charge grid values for B-spline */



/*     Ewald grid needs even numbers with factors of 2, 3 and 5 */



/*     default boundary treatment, B-spline order and grid density */

    s_copy(ewald_1.boundary, "TINFOIL", (ftnlen)7, (ftnlen)7);
    pme_1.bsorder = 5;
    dens = 1.2;

/*     estimate an optimal value for the Ewald coefficient */

    ewaldcof_(&ewald_1.aewald, &cutoff_1.ewaldcut);

/*     set the system extent for nonperiodic Ewald summation */

    if (! bound_1.use_bounds__) {
	extent_(&rmax);
	boxes_1.xbox = (rmax + cutoff_1.ewaldcut) * 2.;
	boxes_1.ybox = boxes_1.xbox;
	boxes_1.zbox = boxes_1.xbox;
	boxes_1.alpha = 90.;
	boxes_1.beta = 90.;
	boxes_1.gamma = 90.;
	boxes_1.orthogonal = TRUE_;
	lattice_();
	s_copy(ewald_1.boundary, "NONE", (ftnlen)7, (ftnlen)4);
	dens = .75;
    }

/*     get default grid size from periodic system dimensions */

    delta = 1e-8;
    ifft1 = (integer) (boxes_1.xbox * dens - delta) + 1;
    ifft2 = (integer) (boxes_1.ybox * dens - delta) + 1;
    ifft3 = (integer) (boxes_1.zbox * dens - delta) + 1;

/*     search keywords for Ewald summation commands */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "EWALD-ALPHA ", (ftnlen)12, (ftnlen)12) == 0) {
	    i__2 = s_rsli(&io___13);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&ewald_1.aewald, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	} else if (s_cmp(keyword, "EWALD-BOUNDARY ", (ftnlen)15, (ftnlen)15) 
		== 0) {
	    s_copy(ewald_1.boundary, "VACUUM", (ftnlen)7, (ftnlen)6);
	} else if (s_cmp(keyword, "PME-GRID ", (ftnlen)9, (ftnlen)9) == 0) {
	    ifft1 = 0;
	    ifft2 = 0;
	    ifft3 = 0;
	    i__2 = s_rsli(&io___14);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ifft1, (ftnlen)sizeof(
		    integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ifft2, (ftnlen)sizeof(
		    integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ifft3, (ftnlen)sizeof(
		    integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
L10:
	    if (ifft2 == 0) {
		ifft2 = ifft1;
	    }
	    if (ifft3 == 0) {
		ifft3 = ifft1;
	    }
	} else if (s_cmp(keyword, "PME-ORDER ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    i__2 = s_rsli(&io___15);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&pme_1.bsorder, (ftnlen)
		    sizeof(integer));
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
	}
L20:
	;
    }

/*     grid size must be even, with prime factors of 2, 3 and 5 */

    pme_1.nfft1 = 100;
    pme_1.nfft2 = 100;
    pme_1.nfft3 = 100;
    for (i__ = 54; i__ >= 1; --i__) {
	k = multi[i__ - 1];
	if (k <= 100) {
	    if (k >= ifft1) {
		pme_1.nfft1 = k;
	    }
	    if (k >= ifft2) {
		pme_1.nfft2 = k;
	    }
	    if (k >= ifft3) {
		pme_1.nfft3 = k;
	    }
	}
    }

/*     check the B-spline order and charge grid dimension */

    if (pme_1.bsorder > 10) {
	io___17.ciunit = iounit_1.iout;
	s_wsfe(&io___17);
	e_wsfe();
	fatal_();
    }
/* Computing MAX */
    i__1 = max(pme_1.nfft1,pme_1.nfft2);
    if (max(i__1,pme_1.nfft3) > 100) {
	io___18.ciunit = iounit_1.iout;
	s_wsfe(&io___18);
	e_wsfe();
	fatal_();
    } else if (pme_1.nfft1 < ifft1 || pme_1.nfft2 < ifft2 || pme_1.nfft3 < 
	    ifft3) {
	io___19.ciunit = iounit_1.iout;
	s_wsfe(&io___19);
	e_wsfe();
    }

/*     initialize the PME arrays that can be precomputed */

    moduli_();
    fftsetup_();

/*     print a message listing some of the Ewald parameters */

    if (inform_1.verbose) {
	io___20.ciunit = iounit_1.iout;
	s_wsfe(&io___20);
	do_fio(&c__1, (char *)&ewald_1.aewald, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&pme_1.nfft1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&pme_1.nfft2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&pme_1.nfft3, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&pme_1.bsorder, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    return 0;
} /* kewald_ */

#undef keyline_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine ewaldcof  --  estimation of Ewald coefficient  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "ewaldcof" finds a value of the Ewald coefficient such that */
/*     all terms beyond the specified cutoff distance will have an */
/*     value less than a specified tolerance */


/* Subroutine */ int ewaldcof_(doublereal *alpha, doublereal *cutoff)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k;
    static doublereal x, y, eps, xhi, xlo;
    extern doublereal erfc_(doublereal *);
    static doublereal ratio;



/*     set the tolerance value; use of 1.0d-8 results in large */
/*     Ewald coefficients that ensure continuity in the gradient */

    eps = 1e-8;

/*     get approximate value from cutoff and tolerance */

    ratio = eps + 1.;
    x = .5;
    i__ = 0;
    while(ratio >= eps) {
	++i__;
	x *= 2.;
	y = x * *cutoff;
	ratio = erfc_(&y) / *cutoff;
    }

/*     use a binary search to refine the coefficient */

    k = i__ + 60;
    xlo = 0.;
    xhi = x;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x = (xlo + xhi) / 2.;
	y = x * *cutoff;
	ratio = erfc_(&y) / *cutoff;
	if (ratio >= eps) {
	    xlo = x;
	} else {
	    xhi = x;
	}
    }
    *alpha = x;
    return 0;
} /* ewaldcof_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine extent  --  find maximum interatomic distance  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "extent" finds the largest interatomic distance in a system */


/* Subroutine */ int extent_(doublereal *rmax)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, k;
    static doublereal r2, xi, yi, zi, xk, yk, zk;



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




/*     search all atom pairs to find the largest distance */

    *rmax = 0.;
    i__1 = atoms_1.n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	i__2 = atoms_1.n;
	for (k = i__ + 1; k <= i__2; ++k) {
	    xk = atoms_1.x[k - 1];
	    yk = atoms_1.y[k - 1];
	    zk = atoms_1.z__[k - 1];
/* Computing 2nd power */
	    d__1 = xi - xk;
/* Computing 2nd power */
	    d__2 = yi - yk;
/* Computing 2nd power */
	    d__3 = zi - zk;
	    r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	    *rmax = max(r2,*rmax);
	}
    }
    *rmax = sqrt(*rmax);
    return 0;
} /* extent_ */



/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine moduli  --  store the inverse DFT moduli  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "moduli" sets the moduli of the inverse discrete Fourier */
/*     transform of the B-splines */


/* Subroutine */ int moduli_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal x, array[10];
    extern /* Subroutine */ int dftmod_(doublereal *, doublereal *, integer *,
	     integer *), bspline_(doublereal *, integer *, doublereal *);
    static doublereal bsarray[100];



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
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  pme.i  --  values for particle mesh Ewald summation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxorder   maximum order of the B-spline approximation */
/*     maxprime   maximum number of prime factors of FFT dimension */
/*     maxtable   maximum size of the FFT table array */

/*     bsmod1     B-spline moduli along the a-axis direction */
/*     bsmod2     B-spline moduli along the b-axis direction */
/*     bsmod3     B-spline moduli along the c-axis direction */
/*     table      intermediate array used by the FFT calculation */
/*     qgrid      values on the particle mesh Ewald charge grid */
/*     qfac       prefactors for particle mesh Ewald charge grid */
/*     thetai1    B-spline coefficients along the a-axis */
/*     thetai2    B-spline coefficients along the b-axis */
/*     thetai3    B-spline coefficients along the c-axis */
/*     nfft1      number of grid points along the a-axis direction */
/*     nfft2      number of grid points along the b-axis direction */
/*     nfft3      number of grid points along the c-axis direction */
/*     bsorder    order of the PME B-spline approximation */
/*     iprime     prime factorization of each FFT dimension */
/*     igrid      initial Ewald charge grid values for B-spline */




/*     compute and load the moduli values */

    x = 0.;
    bspline_(&x, &pme_1.bsorder, array);
    for (i__ = 1; i__ <= 100; ++i__) {
	bsarray[i__ - 1] = 0.;
    }
    i__1 = pme_1.bsorder + 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	bsarray[i__ - 1] = array[i__ - 2];
    }
    dftmod_(pme_1.bsmod1, bsarray, &pme_1.nfft1, &pme_1.bsorder);
    dftmod_(pme_1.bsmod2, bsarray, &pme_1.nfft2, &pme_1.bsorder);
    dftmod_(pme_1.bsmod3, bsarray, &pme_1.nfft3, &pme_1.bsorder);
    return 0;
} /* moduli_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine dftmod  --  discrete Fourier transform modulus  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "dftmod" computes the modulus of the discrete Fourier transform */
/*     of "bsarray" and stores it in "bsmod" */


/* Subroutine */ int dftmod_(doublereal *bsmod, doublereal *bsarray, integer *
	nfft, integer *order)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), pow_di(doublereal *, integer *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal arg, eps, sum1, sum2, zeta;
    static integer jcut, order2;
    static doublereal factor;



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




/*     get the modulus of the discrete Fourier transform */

    /* Parameter adjustments */
    --bsarray;
    --bsmod;

    /* Function Body */
    factor = 6.2831853071795862 / (doublereal) (*nfft);
    i__1 = *nfft;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum1 = 0.;
	sum2 = 0.;
	i__2 = *nfft;
	for (j = 1; j <= i__2; ++j) {
	    arg = factor * (doublereal) ((i__ - 1) * (j - 1));
	    sum1 += bsarray[j] * cos(arg);
	    sum2 += bsarray[j] * sin(arg);
	}
/* Computing 2nd power */
	d__1 = sum1;
/* Computing 2nd power */
	d__2 = sum2;
	bsmod[i__] = d__1 * d__1 + d__2 * d__2;
    }

/*     fix for exponential Euler spline interpolation failure */

    eps = 1e-7;
    if (bsmod[1] < eps) {
	bsmod[1] = bsmod[2] * .5;
    }
    i__1 = *nfft - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (bsmod[i__] < eps) {
	    bsmod[i__] = (bsmod[i__ - 1] + bsmod[i__ + 1]) * .5;
	}
    }
    if (bsmod[*nfft] < eps) {
	bsmod[*nfft] = bsmod[*nfft - 1] * .5;
    }

/*     compute and apply the optimal zeta coefficient */

    jcut = 50;
    order2 = *order << 1;
    i__1 = *nfft;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = i__ - 1;
	if (i__ > *nfft / 2) {
	    k -= *nfft;
	}
	if (k == 0) {
	    zeta = 1.;
	} else {
	    sum1 = 1.;
	    sum2 = 1.;
	    factor = (doublereal) k * 3.141592653589793238 / (doublereal) (*
		    nfft);
	    i__2 = jcut;
	    for (j = 1; j <= i__2; ++j) {
		arg = factor / (factor + (doublereal) j * 
			3.141592653589793238);
		sum1 += pow_di(&arg, order);
		sum2 += pow_di(&arg, &order2);
	    }
	    i__2 = jcut;
	    for (j = 1; j <= i__2; ++j) {
		arg = factor / (factor - (doublereal) j * 
			3.141592653589793238);
		sum1 += pow_di(&arg, order);
		sum2 += pow_di(&arg, &order2);
	    }
	    zeta = sum2 / sum1;
	}
/* Computing 2nd power */
	d__1 = zeta;
	bsmod[i__] *= d__1 * d__1;
    }
    return 0;
} /* dftmod_ */

