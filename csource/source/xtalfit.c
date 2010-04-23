/* xtalfit.f -- translated by f2c (version 20050501).
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
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

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
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

struct {
    doublereal e0_lattice__, moment_0__;
    integer nxtal, nvary, ivary[50], vary[100]	/* was [2][50] */, iresid[100]
	    ;
    char rsdtyp[2000], vartyp[1000];
} xtal1_;

#define xtal1_1 xtal1_

struct {
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

struct {
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

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
    doublereal pchg[25000];
    integer nion, iion[25000], jion[25000], kion[25000], chglist[25000];
} charge_;

#define charge_1 charge_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal bdpl[50000], sdpl[50000];
    integer ndipole, idpl[100000]	/* was [2][50000] */;
} dipole_;

#define dipole_1 dipole_

struct {
    doublereal xfrac[25000], yfrac[25000], zfrac[25000];
} fracs_;

#define fracs_1 fracs_

struct {
    doublereal rad[5000], eps[5000], rad4[5000], eps4[5000], reduct[5000];
} kvdws_;

#define kvdws_1 kvdws_

struct {
    doublereal radmin[1000000]	/* was [1000][1000] */, epsilon[1000000]	
	    /* was [1000][1000] */, radmin4[1000000]	/* was [1000][1000] */
	    , epsilon4[1000000]	/* was [1000][1000] */, radhbnd[1000000]	
	    /* was [1000][1000] */, epshbnd[1000000]	/* was [1000][1000] */
	    , kred[25000];
    integer ired[25000], nvdw, ivdw[25000], jvdw[25000], nvt, ivt[25000], jvt[
	    25000];
} vdw_;

#define vdw_1 vdw_

struct {
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__100 = 100;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  program xtalfit  --  fit parameters to crystal structures  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "xtalfit" computes an optimized set of potential energy */
/*     parameters for user specified van der Waals and electrostatic */
/*     interactions by fitting to crystal structure, lattice energy */
/*     and monomer dipole moment data */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 The Following Parameters can be Fit fo"
	    "r\002,\002 each Atom Type :\002,//,4x,\002(1) Van der Waals Atom"
	    "ic Radius\002,/,4x,\002(2) Van der Waals Well Depth\002,/,4x,"
	    "\002(3) Hydrogen Atom Reduction Factor\002,/,4x,\002(4) Atomic P"
	    "artial Charge\002,/,4x,\002(5) Bond Dipole Moment Magnitude\002,"
	    "/,4x,\002(6) Bond Dipole Moment Position\002)";
    static char fmt_30[] = "(/,\002 Enter Parameter Type then Atom Clas"
	    "s\002,\002 or Type(s) :  \002,$)";
    static char fmt_40[] = "(a120)";
    static char fmt_70[] = "(/,\002 Enter RMS Gradient Termination Criterio"
	    "n\002,\002 [0.1] :  \002,$)";
    static char fmt_80[] = "(f20.0)";
    static char fmt_100[] = "(/,\002 Enter Number of Crystals to be Used [1]"
	    " :  \002,$)";
    static char fmt_110[] = "(i10)";
    static char fmt_130[] = "(/,\002 Enter Lattice Energy Value [<CR> to omi"
	    "t]\002,\002 :  \002,$)";
    static char fmt_140[] = "(f20.0)";
    static char fmt_160[] = "(/,\002 Enter Dipole Moment Value [<CR> to om"
	    "it]\002,\002 :   \002,$)";
    static char fmt_170[] = "(f20.0)";
    static char fmt_180[] = "(/,\002 File Name of Crystal Structure\002,i4"
	    ",\002 :  \002,5x,a35,/,\002 Number of Molecules per Unit Cell "
	    ":\002,i13)";
    static char fmt_190[] = "(\002 Value of Crystal Lattice Energy :  \002,f"
	    "13.2)";
    static char fmt_200[] = "(\002 Value of Molecular Dipole Moment : \002,f"
	    "13.2)";
    static char fmt_210[] = "(/,\002 Initial Values of the Parameters :\002,"
	    "/)";
    static char fmt_220[] = "(3x,\002(\002,i2,\002)\002,2x,a20,\002Atom Cl"
	    "ass\002,i5,4x,f12.4)";
    static char fmt_230[] = "(3x,\002(\002,i2,\002)\002,2x,a20,\002Atom Ty"
	    "pe \002,i5,4x,f12.4)";
    static char fmt_240[] = "(3x,\002(\002,i2,\002)\002,2x,a20,\002Bond Ty"
	    "pe \002,2i5,f12.4)";
    static char fmt_250[] = "(/,\002 Final Values of the Parameters and Scal"
	    "ed\002,\002 Derivatives :\002,/)";
    static char fmt_260[] = "(3x,\002(\002,i2,\002)\002,2x,a20,\002Atom Cl"
	    "ass\002,i5,4x,2f12.4)";
    static char fmt_270[] = "(3x,\002(\002,i2,\002)\002,2x,a20,\002Atom Ty"
	    "pe \002,i5,4x,2f12.4)";
    static char fmt_280[] = "(3x,\002(\002,i2,\002)\002,2x,a20,\002Bond Ty"
	    "pe \002,2i5,2f12.4)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsli(icilist *), do_lio(integer 
	    *, integer *, char *, ftnlen), e_rsli(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal jacobian[5000]	/* was [100][50] */;
    extern /* Subroutine */ int mechanic_(void);
    static doublereal f[100], g[50];
    static integer i__;
    static doublereal xx[50], xhi[50], xlo[50];
    static integer atom1, atom2;
    static char label[20*6];
    extern /* Subroutine */ int final_(void);
    static integer ixtal;
    static logical exist, query;
    static char record[120];
    static integer nresid;
    static doublereal grdmin;
    extern /* Subroutine */ int potoff_(void);
    static char string[120];
    extern /* Subroutine */ int square_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, U_fp, U_fp), getxyz_(void);
    static integer prmtyp;
    extern /* Subroutine */ int initial_(void), nextarg_(char *, logical *, 
	    ftnlen);
    extern /* Subroutine */ int xtalerr_();
    extern /* Subroutine */ int xtalprm_(char *, integer *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int xtalwrt_();

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, fmt_10, 0 };
    static icilist io___9 = { 1, string, 1, 0, 120, 1 };
    static icilist io___10 = { 1, string, 1, 0, 120, 1 };
    static icilist io___11 = { 1, string, 1, 0, 120, 1 };
    static cilist io___12 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_40, 0 };
    static icilist io___15 = { 1, record, 1, 0, 120, 1 };
    static icilist io___17 = { 1, string, 1, 0, 120, 1 };
    static cilist io___18 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_80, 0 };
    static icilist io___20 = { 1, string, 1, 0, 120, 1 };
    static cilist io___21 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_110, 0 };
    static icilist io___24 = { 1, string, 1, 0, 120, 1 };
    static cilist io___25 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_140, 0 };
    static icilist io___27 = { 1, string, 1, 0, 120, 1 };
    static cilist io___28 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___46 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___48 = { 0, 0, 0, fmt_280, 0 };



#define vary_ref(a_1,a_2) xtal1_1.vary[(a_2)*2 + a_1 - 3]
#define label_ref(a_0,a_1) &label[(a_1)*20 + a_0 - 20]
#define rsdtyp_ref(a_0,a_1) &xtal1_1.rsdtyp[(a_1)*20 + a_0 - 20]
#define vartyp_ref(a_0,a_1) &xtal1_1.vartyp[(a_1)*20 + a_0 - 20]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  xtals.i  --  crystal structures for parameter fitting  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     e0_lattice   ideal lattice energy for the current crystal */
/*     moment_0     ideal dipole moment for monomer from crystal */
/*     nxtal        number of crystal structures to be stored */
/*     nvary        number of potential parameters to optimize */
/*     ivary        index for the types of potential parameters */
/*     vary         atom numbers involved in potential parameters */
/*     iresid       crystal structure to which each residual refers */
/*     rsdtyp       experimental variable for each of the residuals */
/*     vartyp       type of potential parameter to be optimized */




/*     initialize some variables and print the header information */

    initial_();
    xtal1_1.nvary = 0;
    nresid = 0;
    io___2.ciunit = iounit_1.iout;
    s_wsfe(&io___2);
    e_wsfe();

/*     get the types of potential parameters to be optimized */

    query = TRUE_;
    while(query) {
	prmtyp = -1;
	atom1 = 0;
	atom2 = 0;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___9);
	    if (i__1 != 0) {
		goto L20;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&prmtyp, (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L20;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L20;
	    }
	}
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___10);
	    if (i__1 != 0) {
		goto L20;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&atom1, (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L20;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L20;
	    }
	}
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___11);
	    if (i__1 != 0) {
		goto L20;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&atom2, (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L20;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L20;
	    }
	}
L20:
	if (prmtyp != 0) {
	    prmtyp = 0;
	    io___12.ciunit = iounit_1.iout;
	    s_wsfe(&io___12);
	    e_wsfe();
	    io___13.ciunit = iounit_1.input;
	    s_rsfe(&io___13);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    i__1 = s_rsli(&io___15);
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&prmtyp, (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&atom1, (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&atom2, (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L50;
	    }
L50:
	    ;
	}
	if (prmtyp == 0) {
	    query = FALSE_;
	} else {
	    query = TRUE_;
	    ++xtal1_1.nvary;
	    xtal1_1.ivary[xtal1_1.nvary - 1] = prmtyp;
	    if (prmtyp < 5) {
		vary_ref(1, xtal1_1.nvary) = atom1;
	    } else {
		vary_ref(1, xtal1_1.nvary) = min(atom1,atom2);
		vary_ref(2, xtal1_1.nvary) = max(atom1,atom2);
	    }
	}
    }

/*     get the termination criterion as RMS gradient over parameters */

    grdmin = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___17);
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&grdmin, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L60;
	}
    }
L60:
    if (grdmin <= 0.) {
	io___18.ciunit = iounit_1.iout;
	s_wsfe(&io___18);
	e_wsfe();
	io___19.ciunit = iounit_1.input;
	s_rsfe(&io___19);
	do_fio(&c__1, (char *)&grdmin, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (grdmin <= 0.) {
	grdmin = .1;
    }

/*     get number of crystal structures to include in optimization */

    xtal1_1.nxtal = 0;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___20);
	if (i__1 != 0) {
	    goto L90;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&xtal1_1.nxtal, (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L90;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L90;
	}
    }
L90:
    if (xtal1_1.nxtal <= 0) {
	io___21.ciunit = iounit_1.iout;
	s_wsfe(&io___21);
	e_wsfe();
	io___22.ciunit = iounit_1.input;
	s_rsfe(&io___22);
	do_fio(&c__1, (char *)&xtal1_1.nxtal, (ftnlen)sizeof(integer));
	e_rsfe();
    }
    if (xtal1_1.nxtal == 0) {
	xtal1_1.nxtal = 1;
    }

/*     get the structural data for each crystal in turn */

    i__1 = xtal1_1.nxtal;
    for (ixtal = 1; ixtal <= i__1; ++ixtal) {
	getxyz_();
	mechanic_();

/*     get an ideal value for the crystal lattice energy */

	query = TRUE_;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__2 = s_rsli(&io___24);
	    if (i__2 != 0) {
		goto L120;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&xtal1_1.e0_lattice__, (
		    ftnlen)sizeof(doublereal));
	    if (i__2 != 0) {
		goto L120;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L120;
	    }
	    query = FALSE_;
	}
L120:
	if (query) {
	    io___25.ciunit = iounit_1.iout;
	    s_wsfe(&io___25);
	    e_wsfe();
	    io___26.ciunit = iounit_1.input;
	    s_rsfe(&io___26);
	    do_fio(&c__1, (char *)&xtal1_1.e0_lattice__, (ftnlen)sizeof(
		    doublereal));
	    e_rsfe();
	}
	if (xtal1_1.e0_lattice__ > 0.) {
	    xtal1_1.e0_lattice__ = -xtal1_1.e0_lattice__;
	}

/*     get an ideal value for the monomer dipole moment */

	query = TRUE_;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__2 = s_rsli(&io___27);
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&xtal1_1.moment_0__, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L150;
	    }
	    query = FALSE_;
	}
L150:
	if (query) {
	    io___28.ciunit = iounit_1.iout;
	    s_wsfe(&io___28);
	    e_wsfe();
	    io___29.ciunit = iounit_1.input;
	    s_rsfe(&io___29);
	    do_fio(&c__1, (char *)&xtal1_1.moment_0__, (ftnlen)sizeof(
		    doublereal));
	    e_rsfe();
	}

/*     set the types of residuals for use in optimization */

	for (i__ = 1; i__ <= 6; ++i__) {
	    xtal1_1.iresid[nresid + i__ - 1] = ixtal;
	}
	s_copy(rsdtyp_ref(0, nresid + 1), "Force a-Axis", (ftnlen)20, (ftnlen)
		12);
	s_copy(rsdtyp_ref(0, nresid + 2), "Force b-Axis", (ftnlen)20, (ftnlen)
		12);
	s_copy(rsdtyp_ref(0, nresid + 3), "Force c-Axis", (ftnlen)20, (ftnlen)
		12);
	s_copy(rsdtyp_ref(0, nresid + 4), "Force Alpha", (ftnlen)20, (ftnlen)
		11);
	s_copy(rsdtyp_ref(0, nresid + 5), "Force Beta", (ftnlen)20, (ftnlen)
		10);
	s_copy(rsdtyp_ref(0, nresid + 6), "Force Gamma", (ftnlen)20, (ftnlen)
		11);
	nresid += 6;

/*     print molecules per unit cell, lattice energy and dipole */

	io___31.ciunit = iounit_1.iout;
	s_wsfe(&io___31);
	do_fio(&c__1, (char *)&ixtal, (ftnlen)sizeof(integer));
	do_fio(&c__1, files_1.filename, (ftnlen)35);
	do_fio(&c__1, (char *)&molcul_1.nmol, (ftnlen)sizeof(integer));
	e_wsfe();
	if (xtal1_1.e0_lattice__ != 0.) {
	    ++nresid;
	    xtal1_1.iresid[nresid - 1] = ixtal;
	    s_copy(rsdtyp_ref(0, nresid), "Lattice Energy", (ftnlen)20, (
		    ftnlen)14);
	    io___32.ciunit = iounit_1.iout;
	    s_wsfe(&io___32);
	    do_fio(&c__1, (char *)&xtal1_1.e0_lattice__, (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
	if (xtal1_1.moment_0__ != 0.) {
	    ++nresid;
	    xtal1_1.iresid[nresid - 1] = ixtal;
	    s_copy(rsdtyp_ref(0, nresid), "Dipole Moment", (ftnlen)20, (
		    ftnlen)13);
	    io___33.ciunit = iounit_1.iout;
	    s_wsfe(&io___33);
	    do_fio(&c__1, (char *)&xtal1_1.moment_0__, (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}

/*     set the initial values of the parameters */

	xtalprm_("STORE", &ixtal, xx, (ftnlen)5);
    }

/*     turn off all local interactions and extra terms */

    potoff_();
    potent_1.use_vdw__ = TRUE_;
    potent_1.use_charge__ = TRUE_;
    potent_1.use_chgdpl__ = TRUE_;
    potent_1.use_dipole__ = TRUE_;

/*     types of variables for use in optimization */

    s_copy(label_ref(0, 1), "Atomic Radius", (ftnlen)20, (ftnlen)13);
    s_copy(label_ref(0, 2), "Well Depth", (ftnlen)20, (ftnlen)10);
    s_copy(label_ref(0, 3), "H Reduction", (ftnlen)20, (ftnlen)11);
    s_copy(label_ref(0, 4), "Partial Charge", (ftnlen)20, (ftnlen)14);
    s_copy(label_ref(0, 5), "Dipole Magnitude", (ftnlen)20, (ftnlen)16);
    s_copy(label_ref(0, 6), "Dipole Position", (ftnlen)20, (ftnlen)15);
    i__1 = xtal1_1.nvary;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(vartyp_ref(0, i__), label_ref(0, xtal1_1.ivary[i__ - 1]), (
		ftnlen)20, (ftnlen)20);
    }

/*     print the initial parameter values */

    io___36.ciunit = iounit_1.iout;
    s_wsfe(&io___36);
    e_wsfe();
    i__1 = xtal1_1.nvary;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (xtal1_1.ivary[i__ - 1] <= 3) {
	    io___37.ciunit = iounit_1.iout;
	    s_wsfe(&io___37);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, vartyp_ref(0, i__), (ftnlen)20);
	    do_fio(&c__1, (char *)&vary_ref(1, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xx[i__ - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (xtal1_1.ivary[i__ - 1] == 4 || xtal1_1.ivary[i__ - 1] == 5)
		 {
	    io___38.ciunit = iounit_1.iout;
	    s_wsfe(&io___38);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, vartyp_ref(0, i__), (ftnlen)20);
	    do_fio(&c__1, (char *)&vary_ref(1, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xx[i__ - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (xtal1_1.ivary[i__ - 1] == 6) {
	    io___39.ciunit = iounit_1.iout;
	    s_wsfe(&io___39);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, vartyp_ref(0, i__), (ftnlen)20);
	    do_fio(&c__1, (char *)&vary_ref(1, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&vary_ref(2, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xx[i__ - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     set upper and lower bounds based on the parameter type */

    i__1 = xtal1_1.nvary;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (xtal1_1.ivary[i__ - 1] == 1) {
	    xlo[i__ - 1] = xx[i__ - 1] * .9;
	    xhi[i__ - 1] = xx[i__ - 1] * 1.1;
	} else if (xtal1_1.ivary[i__ - 1] == 2) {
	    xlo[i__ - 1] = xx[i__ - 1] * .8;
	    xhi[i__ - 1] = xx[i__ - 1] * 1.2;
	} else if (xtal1_1.ivary[i__ - 1] == 3) {
	    xlo[i__ - 1] = xx[i__ - 1] * .8;
	    xhi[i__ - 1] = xx[i__ - 1] * 1.2;
	    if (xhi[i__ - 1] > 1.) {
		xhi[i__ - 1] = 1.;
	    }
	} else if (xtal1_1.ivary[i__ - 1] == 4) {
/* Computing MIN */
	    d__2 = xx[i__ - 1] - (d__1 = xx[i__ - 1], abs(d__1)) * .2, d__3 = 
		    xx[i__ - 1] - .2;
	    xlo[i__ - 1] = min(d__2,d__3);
/* Computing MAX */
	    d__2 = xx[i__ - 1] + (d__1 = xx[i__ - 1], abs(d__1)) * .2, d__3 = 
		    xx[i__ - 1] + .2;
	    xhi[i__ - 1] = max(d__2,d__3);
	} else if (xtal1_1.ivary[i__ - 1] == 5) {
/* Computing MIN */
	    d__2 = xx[i__ - 1] - (d__1 = xx[i__ - 1], abs(d__1)) * .2, d__3 = 
		    xx[i__ - 1] - .2;
	    xlo[i__ - 1] = min(d__2,d__3);
/* Computing MAX */
	    d__2 = xx[i__ - 1] + (d__1 = xx[i__ - 1], abs(d__1)) * .2, d__3 = 
		    xx[i__ - 1] + .2;
	    xhi[i__ - 1] = max(d__2,d__3);
	} else if (xtal1_1.ivary[i__ - 1] == 6) {
	    xlo[i__ - 1] = xx[i__ - 1] * .8;
	    xhi[i__ - 1] = xx[i__ - 1] * 1.2;
	    if (xlo[i__ - 1] < .2) {
		xhi[i__ - 1] = .2;
	    }
	    if (xhi[i__ - 1] > .8) {
		xhi[i__ - 1] = .8;
	    }
	}
    }

/*     use nonlinear least squares to refine the parameters */

    square_(&nresid, &xtal1_1.nvary, xlo, xhi, xx, f, g, jacobian, &c__100, &
	    grdmin, (U_fp)xtalerr_, (U_fp)xtalwrt_);

/*     print the final parameter values */

    io___45.ciunit = iounit_1.iout;
    s_wsfe(&io___45);
    e_wsfe();
    i__1 = xtal1_1.nvary;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (xtal1_1.ivary[i__ - 1] <= 3) {
	    io___46.ciunit = iounit_1.iout;
	    s_wsfe(&io___46);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, vartyp_ref(0, i__), (ftnlen)20);
	    do_fio(&c__1, (char *)&vary_ref(1, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xx[i__ - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g[i__ - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (xtal1_1.ivary[i__ - 1] == 4 || xtal1_1.ivary[i__ - 1] == 5)
		 {
	    io___47.ciunit = iounit_1.iout;
	    s_wsfe(&io___47);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, vartyp_ref(0, i__), (ftnlen)20);
	    do_fio(&c__1, (char *)&vary_ref(1, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xx[i__ - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g[i__ - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (xtal1_1.ivary[i__ - 1] == 6) {
	    io___48.ciunit = iounit_1.iout;
	    s_wsfe(&io___48);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, vartyp_ref(0, i__), (ftnlen)20);
	    do_fio(&c__1, (char *)&vary_ref(1, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&vary_ref(2, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xx[i__ - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&g[i__ - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef vartyp_ref
#undef rsdtyp_ref
#undef label_ref
#undef vary_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine xtalprm  --  energy/optimization conversion  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "xtalprm" stores or retrieves a crystal structure; used */
/*     to make a previously stored structure the currently active */
/*     structure, or to store a structure for later use; only */
/*     provides for the intermolecular energy terms */


/* Subroutine */ int xtalprm_(char *mode, integer *ixtal, doublereal *xx, 
	ftnlen mode_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int molecule_(void);
    static integer ndipoles[10], i__, j, k;
    static doublereal moment_0s__[10];
    static integer ns[10];
    static doublereal xs[250000]	/* was [25000][10] */, ys[250000]	
	    /* was [25000][10] */, zs[250000]	/* was [25000][10] */, 
	    e0_lattices__[10];
    static integer i12s[1000000]	/* was [4][25000][10] */, i13s[
	    3000000]	/* was [12][25000][10] */, i14s[9000000]	/* 
	    was [36][25000][10] */, n12s[250000]	/* was [25000][10] */,
	     n13s[250000]	/* was [25000][10] */, n14s[250000]	/* 
	    was [25000][10] */;
    static doublereal xmid, ymid, zmid;
    static integer atom1, atom2;
    static doublereal betas[10];
    static integer ireds[250000]	/* was [25000][10] */, idpls[1000000]	
	    /* was [2][50000][10] */;
    static doublereal pchgs[250000]	/* was [25000][10] */, bdpls[500000]	
	    /* was [50000][10] */, kreds[250000]	/* was [25000][10] */;
    static char names[3*25000*10];
    static integer iions[250000]	/* was [25000][10] */;
    static doublereal sdpls[500000]	/* was [50000][10] */;
    static integer nions[10], ivdws[250000]	/* was [25000][10] */, nvdws[
	    10];
    static doublereal xboxs[10];
    static integer types[250000]	/* was [25000][10] */;
    static doublereal yboxs[10], zboxs[10], gammas[10], alphas[10], xfracs[
	    250000]	/* was [25000][10] */, masses[250000]	/* was [25000]
	    [10] */, yfracs[250000]	/* was [25000][10] */, zfracs[250000]	
	    /* was [25000][10] */;
    extern /* Subroutine */ int bounds_(void);
    static integer prmtyp;
    extern /* Subroutine */ int lattice_(void);
    static integer classes[250000]	/* was [25000][10] */;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define xs_ref(a_1,a_2) xs[(a_2)*25000 + a_1 - 25001]
#define ys_ref(a_1,a_2) ys[(a_2)*25000 + a_1 - 25001]
#define zs_ref(a_1,a_2) zs[(a_2)*25000 + a_1 - 25001]
#define i12s_ref(a_1,a_2,a_3) i12s[((a_3)*25000 + (a_2))*4 + a_1 - 100005]
#define i13s_ref(a_1,a_2,a_3) i13s[((a_3)*25000 + (a_2))*12 + a_1 - \
300013]
#define i14s_ref(a_1,a_2,a_3) i14s[((a_3)*25000 + (a_2))*36 + a_1 - \
900037]
#define n12s_ref(a_1,a_2) n12s[(a_2)*25000 + a_1 - 25001]
#define n13s_ref(a_1,a_2) n13s[(a_2)*25000 + a_1 - 25001]
#define n14s_ref(a_1,a_2) n14s[(a_2)*25000 + a_1 - 25001]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define idpl_ref(a_1,a_2) dipole_1.idpl[(a_2)*2 + a_1 - 3]
#define imol_ref(a_1,a_2) molcul_1.imol[(a_2)*2 + a_1 - 3]
#define vary_ref(a_1,a_2) xtal1_1.vary[(a_2)*2 + a_1 - 3]
#define ireds_ref(a_1,a_2) ireds[(a_2)*25000 + a_1 - 25001]
#define idpls_ref(a_1,a_2,a_3) idpls[((a_3)*50000 + (a_2))*2 + a_1 - \
100003]
#define pchgs_ref(a_1,a_2) pchgs[(a_2)*25000 + a_1 - 25001]
#define bdpls_ref(a_1,a_2) bdpls[(a_2)*50000 + a_1 - 50001]
#define kreds_ref(a_1,a_2) kreds[(a_2)*25000 + a_1 - 25001]
#define names_ref(a_0,a_1,a_2) &names[((a_2)*25000 + a_1)*3 + a_0 - 75003]
#define iions_ref(a_1,a_2) iions[(a_2)*25000 + a_1 - 25001]
#define sdpls_ref(a_1,a_2) sdpls[(a_2)*50000 + a_1 - 50001]
#define ivdws_ref(a_1,a_2) ivdws[(a_2)*25000 + a_1 - 25001]
#define types_ref(a_1,a_2) types[(a_2)*25000 + a_1 - 25001]
#define radmin_ref(a_1,a_2) vdw_1.radmin[(a_2)*1000 + a_1 - 1001]
#define xfracs_ref(a_1,a_2) xfracs[(a_2)*25000 + a_1 - 25001]
#define masses_ref(a_1,a_2) masses[(a_2)*25000 + a_1 - 25001]
#define yfracs_ref(a_1,a_2) yfracs[(a_2)*25000 + a_1 - 25001]
#define zfracs_ref(a_1,a_2) zfracs[(a_2)*25000 + a_1 - 25001]
#define classes_ref(a_1,a_2) classes[(a_2)*25000 + a_1 - 25001]
#define epsilon_ref(a_1,a_2) vdw_1.epsilon[(a_2)*1000 + a_1 - 1001]



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
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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
/*     ##  dipole.i  --  atom & bond dipoles for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     bdpl      magnitude of each of the dipoles (Debyes) */
/*     sdpl      position of each dipole between defining atoms */
/*     ndipole   total number of dipoles in the system */
/*     idpl      numbers of atoms that define each dipole */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  fracs.i  --  atom distances to molecular center of mass  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xfrac   fractional coordinate along a-axis of center of mass */
/*     yfrac   fractional coordinate along b-axis of center of mass */
/*     zfrac   fractional coordinate along c-axis of center of mass */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kvdws.i  --  forcefield parameters for van der Waals terms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     rad      van der Waals radius parameter for each atom type */
/*     eps      van der Waals well depth parameter for each atom type */
/*     rad4     van der Waals radius parameter in 1-4 interactions */
/*     eps4     van der Waals well depth parameter in 1-4 interactions */
/*     reduct   van der Waals reduction factor for each atom type */




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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  vdw.i  --  van der Waals parameters for current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     radmin     minimum energy distance for each atom class pair */
/*     epsilon    well depth parameter for each atom class pair */
/*     radmin4    minimum energy distance for 1-4 interaction pairs */
/*     epsilon4   well depth parameter for 1-4 interaction pairs */
/*     radhbnd    minimum energy distance for hydrogen bonding pairs */
/*     epshbnd    well depth parameter for hydrogen bonding pairs */
/*     kred       value of reduction factor parameter for each atom */
/*     ired       attached atom from which reduction factor is applied */
/*     nvdw       total number van der Waals active sites in the system */
/*     ivdw       number of the atom for each van der Waals active site */
/*     jvdw       type or class index into vdw parameters for each atom */
/*     nvt        number of distinct vdw types/classes in the system */
/*     ivt        type/class index for each distinct vdw type or class */
/*     jvt        frequency of each vdw type or class in the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  xtals.i  --  crystal structures for parameter fitting  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     e0_lattice   ideal lattice energy for the current crystal */
/*     moment_0     ideal dipole moment for monomer from crystal */
/*     nxtal        number of crystal structures to be stored */
/*     nvary        number of potential parameters to optimize */
/*     ivary        index for the types of potential parameters */
/*     vary         atom numbers involved in potential parameters */
/*     iresid       crystal structure to which each residual refers */
/*     rsdtyp       experimental variable for each of the residuals */
/*     vartyp       type of potential parameter to be optimized */




/*     number of atoms, atomic coordinates and other info */

    /* Parameter adjustments */
    --xx;

    /* Function Body */
    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
	ns[*ixtal - 1] = atoms_1.n;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xs_ref(i__, *ixtal) = atoms_1.x[i__ - 1];
	    ys_ref(i__, *ixtal) = atoms_1.y[i__ - 1];
	    zs_ref(i__, *ixtal) = atoms_1.z__[i__ - 1];
	    types_ref(i__, *ixtal) = atoms_1.type__[i__ - 1];
	    classes_ref(i__, *ixtal) = atmtyp_1.class__[i__ - 1];
	    s_copy(names_ref(0, i__, *ixtal), name___ref(0, i__), (ftnlen)3, (
		    ftnlen)3);
	    masses_ref(i__, *ixtal) = atmtyp_1.mass[i__ - 1];
	}
    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
	atoms_1.n = ns[*ixtal - 1];
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    atoms_1.x[i__ - 1] = xs_ref(i__, *ixtal);
	    atoms_1.y[i__ - 1] = ys_ref(i__, *ixtal);
	    atoms_1.z__[i__ - 1] = zs_ref(i__, *ixtal);
	    atoms_1.type__[i__ - 1] = types_ref(i__, *ixtal);
	    atmtyp_1.class__[i__ - 1] = classes_ref(i__, *ixtal);
	    s_copy(name___ref(0, i__), names_ref(0, i__, *ixtal), (ftnlen)3, (
		    ftnlen)3);
	    atmtyp_1.mass[i__ - 1] = masses_ref(i__, *ixtal);
	}
    }

/*     lists of the attached atoms and neighbors */

    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    n12s_ref(i__, *ixtal) = couple_1.n12[i__ - 1];
	    i__2 = couple_1.n12[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		i12s_ref(j, i__, *ixtal) = i12_ref(j, i__);
	    }
	    n13s_ref(i__, *ixtal) = couple_1.n13[i__ - 1];
	    i__2 = couple_1.n13[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		i13s_ref(j, i__, *ixtal) = i13_ref(j, i__);
	    }
	    n14s_ref(i__, *ixtal) = couple_1.n14[i__ - 1];
	    i__2 = couple_1.n14[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		i14s_ref(j, i__, *ixtal) = i14_ref(j, i__);
	    }
	}
    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    couple_1.n12[i__ - 1] = n12s_ref(i__, *ixtal);
	    i__2 = couple_1.n12[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		i12_ref(j, i__) = i12s_ref(j, i__, *ixtal);
	    }
	    couple_1.n13[i__ - 1] = n13s_ref(i__, *ixtal);
	    i__2 = couple_1.n13[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		i13_ref(j, i__) = i13s_ref(j, i__, *ixtal);
	    }
	    couple_1.n14[i__ - 1] = n14s_ref(i__, *ixtal);
	    i__2 = couple_1.n14[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		i14_ref(j, i__) = i14s_ref(j, i__, *ixtal);
	    }
	}
    }

/*     lattice type and unit cell parameters */

    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
	xboxs[*ixtal - 1] = boxes_1.xbox;
	yboxs[*ixtal - 1] = boxes_1.ybox;
	zboxs[*ixtal - 1] = boxes_1.zbox;
	alphas[*ixtal - 1] = boxes_1.alpha;
	betas[*ixtal - 1] = boxes_1.beta;
	gammas[*ixtal - 1] = boxes_1.gamma;
    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
	boxes_1.xbox = xboxs[*ixtal - 1];
	boxes_1.ybox = yboxs[*ixtal - 1];
	boxes_1.zbox = zboxs[*ixtal - 1];
	boxes_1.alpha = alphas[*ixtal - 1];
	boxes_1.beta = betas[*ixtal - 1];
	boxes_1.gamma = gammas[*ixtal - 1];
    }

/*     find number of molecules and atoms in molecules; */
/*     type of crystal lattice; enforce periodic bounds */

    if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
	molecule_();
	lattice_();
	bounds_();
    }

/*     fractional coordinates of molecular center of mass */

    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
	i__1 = molcul_1.nmol;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xmid = 0.;
	    ymid = 0.;
	    zmid = 0.;
	    i__2 = imol_ref(2, i__);
	    for (j = imol_ref(1, i__); j <= i__2; ++j) {
		k = molcul_1.kmol[j - 1];
		xmid += atoms_1.x[k - 1] * atmtyp_1.mass[k - 1];
		ymid += atoms_1.y[k - 1] * atmtyp_1.mass[k - 1];
		zmid += atoms_1.z__[k - 1] * atmtyp_1.mass[k - 1];
	    }
	    zmid /= boxes_1.gamma_term__;
	    ymid = (ymid - zmid * boxes_1.beta_term__) / boxes_1.gamma_sin__;
	    xmid = xmid - ymid * boxes_1.gamma_cos__ - zmid * 
		    boxes_1.beta_cos__;
	    xfracs_ref(i__, *ixtal) = xmid / (boxes_1.xbox * molcul_1.molmass[
		    i__ - 1]);
	    yfracs_ref(i__, *ixtal) = ymid / (boxes_1.ybox * molcul_1.molmass[
		    i__ - 1]);
	    zfracs_ref(i__, *ixtal) = zmid / (boxes_1.zbox * molcul_1.molmass[
		    i__ - 1]);
	}
    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
	i__1 = molcul_1.nmol;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    fracs_1.xfrac[i__ - 1] = xfracs_ref(i__, *ixtal);
	    fracs_1.yfrac[i__ - 1] = yfracs_ref(i__, *ixtal);
	    fracs_1.zfrac[i__ - 1] = zfracs_ref(i__, *ixtal);
	}
    }

/*     number and types of van der Waals sites */

    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
	nvdws[*ixtal - 1] = vdw_1.nvdw;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ivdws_ref(i__, *ixtal) = vdw_1.ivdw[i__ - 1];
	    ireds_ref(i__, *ixtal) = vdw_1.ired[i__ - 1];
	    kreds_ref(i__, *ixtal) = vdw_1.kred[i__ - 1];
	}
    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
	vdw_1.nvdw = nvdws[*ixtal - 1];
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    vdw_1.ivdw[i__ - 1] = ivdws_ref(i__, *ixtal);
	    vdw_1.ired[i__ - 1] = ireds_ref(i__, *ixtal);
	    vdw_1.kred[i__ - 1] = kreds_ref(i__, *ixtal);
	}
    }

/*     number and types of atomic partial charges */

    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
	nions[*ixtal - 1] = charge_1.nion;
	i__1 = charge_1.nion;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iions_ref(i__, *ixtal) = charge_1.iion[i__ - 1];
	    pchgs_ref(i__, *ixtal) = charge_1.pchg[i__ - 1];
	}
    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
	charge_1.nion = nions[*ixtal - 1];
	i__1 = charge_1.nion;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    charge_1.iion[i__ - 1] = iions_ref(i__, *ixtal);
	    charge_1.pchg[i__ - 1] = pchgs_ref(i__, *ixtal);
	}
    }

/*     number and types of bond dipole moments */

    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
	ndipoles[*ixtal - 1] = dipole_1.ndipole;
	i__1 = dipole_1.ndipole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    idpls_ref(1, i__, *ixtal) = idpl_ref(1, i__);
	    idpls_ref(2, i__, *ixtal) = idpl_ref(2, i__);
	    bdpls_ref(i__, *ixtal) = dipole_1.bdpl[i__ - 1];
	    sdpls_ref(i__, *ixtal) = dipole_1.sdpl[i__ - 1];
	}
    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
	dipole_1.ndipole = ndipoles[*ixtal - 1];
	i__1 = dipole_1.ndipole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    idpl_ref(1, i__) = idpls_ref(1, i__, *ixtal);
	    idpl_ref(2, i__) = idpls_ref(2, i__, *ixtal);
	    dipole_1.bdpl[i__ - 1] = bdpls_ref(i__, *ixtal);
	    dipole_1.sdpl[i__ - 1] = sdpls_ref(i__, *ixtal);
	}
    }

/*     values of ideal lattice energy and dipole moment */

    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
	e0_lattices__[*ixtal - 1] = xtal1_1.e0_lattice__;
	moment_0s__[*ixtal - 1] = xtal1_1.moment_0__;
    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
	xtal1_1.e0_lattice__ = e0_lattices__[*ixtal - 1];
	xtal1_1.moment_0__ = moment_0s__[*ixtal - 1];
    }

/*     store or reset values of the optimization variables */

    i__1 = xtal1_1.nvary;
    for (j = 1; j <= i__1; ++j) {
	prmtyp = xtal1_1.ivary[j - 1];
	atom1 = vary_ref(1, j);
	if (prmtyp == 1) {
	    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
		xx[j] = kvdws_1.rad[atom1 - 1];
	    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
		kvdws_1.rad[atom1 - 1] = xx[j];
		for (i__ = 1; i__ <= 1000; ++i__) {
		    radmin_ref(i__, atom1) = kvdws_1.rad[i__ - 1] + 
			    kvdws_1.rad[atom1 - 1];
		    radmin_ref(atom1, i__) = radmin_ref(i__, atom1);
		}
	    }
	} else if (prmtyp == 2) {
	    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
		xx[j] = kvdws_1.eps[atom1 - 1];
	    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
		kvdws_1.eps[atom1 - 1] = (d__1 = xx[j], abs(d__1));
		for (i__ = 1; i__ <= 1000; ++i__) {
		    epsilon_ref(i__, atom1) = sqrt(kvdws_1.eps[i__ - 1] * 
			    kvdws_1.eps[atom1 - 1]);
		    epsilon_ref(atom1, i__) = epsilon_ref(i__, atom1);
		}
	    }
	} else if (prmtyp == 3) {
	    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
		i__2 = atoms_1.n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    if (atmtyp_1.class__[i__ - 1] == atom1) {
			xx[j] = vdw_1.kred[i__ - 1];
			goto L10;
		    }
		}
	    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
		i__2 = atoms_1.n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    if (atmtyp_1.class__[i__ - 1] == atom1) {
			vdw_1.kred[i__ - 1] = xx[j];
		    }
		}
	    }
	} else if (prmtyp == 4) {
	    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
		i__2 = charge_1.nion;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    if (atoms_1.type__[charge_1.iion[i__ - 1] - 1] == atom1) {
			xx[j] = charge_1.pchg[i__ - 1];
			goto L10;
		    }
		}
	    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
		i__2 = atoms_1.n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    if (atoms_1.type__[charge_1.iion[i__ - 1] - 1] == atom1) {
			charge_1.pchg[i__ - 1] = xx[j];
		    }
		}
	    }
	} else if (prmtyp == 5) {
	    atom2 = vary_ref(2, j);
	    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
		i__2 = dipole_1.ndipole;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    if (atoms_1.type__[idpl_ref(1, i__) - 1] == atom1 && 
			    atoms_1.type__[idpl_ref(2, i__) - 1] == atom2) {
			xx[j] = dipole_1.bdpl[i__ - 1];
			goto L10;
		    }
		}
	    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
		i__2 = dipole_1.ndipole;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    if (atoms_1.type__[idpl_ref(1, i__) - 1] == atom1 && 
			    atoms_1.type__[idpl_ref(2, i__) - 1] == atom2) {
			dipole_1.bdpl[i__ - 1] = xx[j];
		    }
		}
	    }
	} else if (prmtyp == 6) {
	    atom2 = vary_ref(2, j);
	    if (s_cmp(mode, "STORE", (ftnlen)5, (ftnlen)5) == 0) {
		i__2 = dipole_1.ndipole;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    if (atoms_1.type__[idpl_ref(1, i__) - 1] == atom1 && 
			    atoms_1.type__[idpl_ref(2, i__) - 1] == atom2) {
			xx[j] = dipole_1.sdpl[i__ - 1];
			goto L10;
		    }
		}
	    } else if (s_cmp(mode, "RESET", (ftnlen)5, (ftnlen)5) == 0) {
		i__2 = dipole_1.ndipole;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    if (atoms_1.type__[idpl_ref(1, i__) - 1] == atom1 && 
			    atoms_1.type__[idpl_ref(2, i__) - 1] == atom2) {
			dipole_1.sdpl[i__ - 1] = xx[j];
		    }
		}
	    }
	}
L10:
	;
    }
    return 0;
} /* xtalprm_ */

#undef epsilon_ref
#undef classes_ref
#undef zfracs_ref
#undef yfracs_ref
#undef masses_ref
#undef xfracs_ref
#undef radmin_ref
#undef types_ref
#undef ivdws_ref
#undef sdpls_ref
#undef iions_ref
#undef names_ref
#undef kreds_ref
#undef bdpls_ref
#undef pchgs_ref
#undef idpls_ref
#undef ireds_ref
#undef vary_ref
#undef imol_ref
#undef idpl_ref
#undef name___ref
#undef n14s_ref
#undef n13s_ref
#undef n12s_ref
#undef i14s_ref
#undef i13s_ref
#undef i12s_ref
#undef zs_ref
#undef ys_ref
#undef xs_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  subroutine xtalerr  --  error function for xtalfit  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     "xtalerr" computes an error function value derived from */
/*     derivatives with respect to lattice parameters, lattice */
/*     energy and monomer dipole moments */


/* Subroutine */ int xtalerr_(integer *nresid, integer *nvaried, doublereal *
	xx, doublereal *resid)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer nion_old__, nvdw_old__;
    extern /* Subroutine */ int xtalmove_(void);
    static integer i__;
    static doublereal e_lattice__;
    static integer i1, i2;
    static doublereal e_monomer__, xi, yi, zi, fik, eps;
    static integer ndipole_old__;
    static doublereal temp;
    static integer n_old__;
    static doublereal x_dpm__, y_dpm__, z_dpm__;
    static integer ixtal;
    static doublereal e_beta__, g_beta__, e_xtal__;
    extern doublereal energy_(void);
    static doublereal e_xbox__, moment, e_ybox__, e_zbox__, g_xbox__, 
	    g_ybox__, g_zbox__, e_gamma__, g_gamma__, e_alpha__, g_alpha__;
    extern /* Subroutine */ int xtalprm_(char *, integer *, doublereal *, 
	    ftnlen);


#define idpl_ref(a_1,a_2) dipole_1.idpl[(a_2)*2 + a_1 - 3]



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
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  dipole.i  --  atom & bond dipoles for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     bdpl      magnitude of each of the dipoles (Debyes) */
/*     sdpl      position of each dipole between defining atoms */
/*     ndipole   total number of dipoles in the system */
/*     idpl      numbers of atoms that define each dipole */




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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  vdw.i  --  van der Waals parameters for current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     radmin     minimum energy distance for each atom class pair */
/*     epsilon    well depth parameter for each atom class pair */
/*     radmin4    minimum energy distance for 1-4 interaction pairs */
/*     epsilon4   well depth parameter for 1-4 interaction pairs */
/*     radhbnd    minimum energy distance for hydrogen bonding pairs */
/*     epshbnd    well depth parameter for hydrogen bonding pairs */
/*     kred       value of reduction factor parameter for each atom */
/*     ired       attached atom from which reduction factor is applied */
/*     nvdw       total number van der Waals active sites in the system */
/*     ivdw       number of the atom for each van der Waals active site */
/*     jvdw       type or class index into vdw parameters for each atom */
/*     nvt        number of distinct vdw types/classes in the system */
/*     ivt        type/class index for each distinct vdw type or class */
/*     jvt        frequency of each vdw type or class in the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  xtals.i  --  crystal structures for parameter fitting  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     e0_lattice   ideal lattice energy for the current crystal */
/*     moment_0     ideal dipole moment for monomer from crystal */
/*     nxtal        number of crystal structures to be stored */
/*     nvary        number of potential parameters to optimize */
/*     ivary        index for the types of potential parameters */
/*     vary         atom numbers involved in potential parameters */
/*     iresid       crystal structure to which each residual refers */
/*     rsdtyp       experimental variable for each of the residuals */
/*     vartyp       type of potential parameter to be optimized */




/*     zero out the number of residual functions */

    /* Parameter adjustments */
    --resid;
    --xx;

    /* Function Body */
    *nresid = 0;

/*     set the values of the potential energy parameters */
/*     and get the energy of the base structure */

    i__1 = xtal1_1.nxtal;
    for (ixtal = 1; ixtal <= i__1; ++ixtal) {
	xtalprm_("RESET", &ixtal, &xx[1], (ftnlen)5);
	e_xtal__ = energy_();

/*     get energy derivative with respect to lattice a-axis */

	eps = 1e-5;
	temp = boxes_1.xbox;
	boxes_1.xbox += eps;
	xtalmove_();
	e_xbox__ = energy_();
	boxes_1.xbox = temp;
	g_xbox__ = (e_xbox__ - e_xtal__) / eps;

/*     get energy derivative with respect to lattice b-axis */

	temp = boxes_1.ybox;
	boxes_1.ybox += eps;
	xtalmove_();
	e_ybox__ = energy_();
	boxes_1.ybox = temp;
	g_ybox__ = (e_ybox__ - e_xtal__) / eps;

/*     get energy derivative with respect to lattice c-axis */

	temp = boxes_1.zbox;
	boxes_1.zbox += eps;
	xtalmove_();
	e_zbox__ = energy_();
	boxes_1.zbox = temp;
	g_zbox__ = (e_zbox__ - e_xtal__) / eps;

/*     get energy derivative with respect to lattice alpha */

	temp = boxes_1.alpha;
	boxes_1.alpha += eps * 57.29577951308232088;
	xtalmove_();
	e_alpha__ = energy_();
	boxes_1.alpha = temp;
	g_alpha__ = (e_alpha__ - e_xtal__) / eps;

/*     get energy derivative with respect to lattice beta */

	temp = boxes_1.beta;
	boxes_1.beta += eps * 57.29577951308232088;
	xtalmove_();
	e_beta__ = energy_();
	boxes_1.beta = temp;
	g_beta__ = (e_beta__ - e_xtal__) / eps;

/*     get energy derivative with respect to lattice gamma */

	temp = boxes_1.gamma;
	boxes_1.gamma += eps * 57.29577951308232088;
	xtalmove_();
	e_gamma__ = energy_();
	boxes_1.gamma = temp;
	g_gamma__ = (e_gamma__ - e_xtal__) / eps;
	xtalmove_();

/*     setup to compute properties of monomer; assumes that */
/*     molecules are contiguous in the coordinates file */

	bound_1.use_bounds__ = FALSE_;
	bound_1.use_replica__ = FALSE_;
	n_old__ = atoms_1.n;
	nvdw_old__ = vdw_1.nvdw;
	nion_old__ = charge_1.nion;
	ndipole_old__ = dipole_1.ndipole;
	atoms_1.n /= molcul_1.nmol;
	vdw_1.nvdw /= molcul_1.nmol;
	charge_1.nion /= molcul_1.nmol;
	dipole_1.ndipole /= molcul_1.nmol;

/*     compute the crystal lattice energy */

	e_monomer__ = energy_();
	e_lattice__ = (e_xtal__ - molcul_1.nmol * e_monomer__) / (doublereal) 
		molcul_1.nmol;

/*     compute the dipole moment of the monomer */

	x_dpm__ = 0.;
	y_dpm__ = 0.;
	z_dpm__ = 0.;
	i__2 = dipole_1.ndipole;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i1 = idpl_ref(1, i__);
	    i2 = idpl_ref(2, i__);
	    xi = atoms_1.x[i2 - 1] - atoms_1.x[i1 - 1];
	    yi = atoms_1.y[i2 - 1] - atoms_1.y[i1 - 1];
	    zi = atoms_1.z__[i2 - 1] - atoms_1.z__[i1 - 1];
	    fik = dipole_1.bdpl[i__ - 1] / sqrt(xi * xi + yi * yi + zi * zi);
	    x_dpm__ -= xi * fik;
	    y_dpm__ -= yi * fik;
	    z_dpm__ -= zi * fik;
	}
/* Computing 2nd power */
	d__1 = x_dpm__;
/* Computing 2nd power */
	d__2 = y_dpm__;
/* Computing 2nd power */
	d__3 = z_dpm__;
	moment = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);

/*     convert back from monomer to full crystal */

	atoms_1.n = n_old__;
	vdw_1.nvdw = nvdw_old__;
	charge_1.nion = nion_old__;
	dipole_1.ndipole = ndipole_old__;
	bound_1.use_bounds__ = TRUE_;
	bound_1.use_replica__ = TRUE_;

/*     compute the residual vector as a collection */
/*     of the forces on the crystal lattice, plus any */
/*     lattice energy and dipole moment deviations */

	++(*nresid);
	resid[*nresid] = g_xbox__;
	++(*nresid);
	resid[*nresid] = g_ybox__;
	++(*nresid);
	resid[*nresid] = g_zbox__;
	++(*nresid);
	resid[*nresid] = g_alpha__;
	++(*nresid);
	resid[*nresid] = g_beta__;
	++(*nresid);
	resid[*nresid] = g_gamma__;
	if (xtal1_1.e0_lattice__ != 0.) {
	    ++(*nresid);
	    resid[*nresid] = e_lattice__ - xtal1_1.e0_lattice__;
	}
	if (xtal1_1.moment_0__ != 0.) {
	    ++(*nresid);
	    resid[*nresid] = (moment - xtal1_1.moment_0__) * 2.;
	}
    }
    return 0;
} /* xtalerr_ */

#undef idpl_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine xtalmove  --  translation of rigid molecules  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "xtalmove" converts fractional to Cartesian coordinates for */
/*     rigid molecules during fitting of force field parameters to */
/*     crystal structure data */


/* Subroutine */ int xtalmove_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal xmid;
    static integer init;
    static doublereal ymid, zmid, xoff[25000], yoff[25000], zoff[25000];
    static integer stop;
    static doublereal weigh;
    extern /* Subroutine */ int lattice_(void);


#define imol_ref(a_1,a_2) molcul_1.imol[(a_2)*2 + a_1 - 3]



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
/*     ##  fracs.i  --  atom distances to molecular center of mass  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xfrac   fractional coordinate along a-axis of center of mass */
/*     yfrac   fractional coordinate along b-axis of center of mass */
/*     zfrac   fractional coordinate along c-axis of center of mass */




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




/*     get values for fractional coordinate interconversion */

    lattice_();

/*     locate the center of mass of each molecule */

    i__1 = molcul_1.nmol;
    for (i__ = 1; i__ <= i__1; ++i__) {
	init = imol_ref(1, i__);
	stop = imol_ref(2, i__);
	xmid = 0.;
	ymid = 0.;
	zmid = 0.;
	i__2 = stop;
	for (j = init; j <= i__2; ++j) {
	    k = molcul_1.kmol[j - 1];
	    weigh = atmtyp_1.mass[k - 1];
	    xmid += atoms_1.x[k - 1] * weigh;
	    ymid += atoms_1.y[k - 1] * weigh;
	    zmid += atoms_1.z__[k - 1] * weigh;
	}
	weigh = molcul_1.molmass[i__ - 1];
	xmid /= weigh;
	ymid /= weigh;
	zmid /= weigh;

/*     save atomic coordinates relative to center of mass */

	i__2 = stop;
	for (j = init; j <= i__2; ++j) {
	    k = molcul_1.kmol[j - 1];
	    xoff[k - 1] = atoms_1.x[k - 1] - xmid;
	    yoff[k - 1] = atoms_1.y[k - 1] - ymid;
	    zoff[k - 1] = atoms_1.z__[k - 1] - zmid;
	}

/*     convert fractional center of mass to Cartesian coordinates */

	xmid = fracs_1.xfrac[i__ - 1] * boxes_1.xbox + fracs_1.yfrac[i__ - 1] 
		* boxes_1.ybox * boxes_1.gamma_cos__ + fracs_1.zfrac[i__ - 1] 
		* boxes_1.zbox * boxes_1.beta_cos__;
	ymid = fracs_1.yfrac[i__ - 1] * boxes_1.ybox * boxes_1.gamma_sin__ + 
		fracs_1.zfrac[i__ - 1] * boxes_1.zbox * boxes_1.beta_term__;
	zmid = fracs_1.zfrac[i__ - 1] * boxes_1.zbox * boxes_1.gamma_term__;

/*     translate coordinates via offset from center of mass */

	i__2 = stop;
	for (j = init; j <= i__2; ++j) {
	    k = molcul_1.kmol[j - 1];
	    atoms_1.x[k - 1] = xoff[k - 1] + xmid;
	    atoms_1.y[k - 1] = yoff[k - 1] + ymid;
	    atoms_1.z__[k - 1] = zoff[k - 1] + zmid;
	}
    }
    return 0;
} /* xtalmove_ */

#undef imol_ref




/*     ######################################################## */
/*     ##                                                    ## */
/*     ##  subroutine xtalwrt  --  write current parameters  ## */
/*     ##                                                    ## */
/*     ######################################################## */


/*     "xtalwrt" is a utility that prints intermediate results */
/*     during fitting of force field parameters to crystal data */


/* Subroutine */ int xtalwrt_(integer *niter, doublereal *xx, doublereal *gs, 
	integer *nresid, doublereal *f)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Parameters and Scaled Derivatives a"
	    "t\002,\002 Iteration\002,i4,\002 :\002,/)";
    static char fmt_20[] = "(3x,\002(\002,i2,\002)\002,2x,a20,\002Atom Clas"
	    "s\002,i5,4x,2f12.4)";
    static char fmt_30[] = "(3x,\002(\002,i2,\002)\002,2x,a20,\002Atom Type"
	    " \002,i5,4x,2f12.4)";
    static char fmt_40[] = "(3x,\002(\002,i2,\002)\002,2x,a20,\002Bond Type"
	    " \002,2i5,2f12.4)";
    static char fmt_50[] = "(/,\002 Residual Error Function Values at Iterat"
	    "ion\002,i4,\002 :\002,/)";
    static char fmt_60[] = "(3x,\002(\002,i2,\002)\002,2x,a20,2x,\002Crysta"
	    "l\002,i4,4x,f12.4)";
    static char fmt_70[] = "()";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__;

    /* Fortran I/O blocks */
    static cilist io___139 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___141 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___142 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___143 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___144 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___145 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___146 = { 0, 0, 0, fmt_70, 0 };



#define vary_ref(a_1,a_2) xtal1_1.vary[(a_2)*2 + a_1 - 3]
#define rsdtyp_ref(a_0,a_1) &xtal1_1.rsdtyp[(a_1)*20 + a_0 - 20]
#define vartyp_ref(a_0,a_1) &xtal1_1.vartyp[(a_1)*20 + a_0 - 20]



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
/*     ##  xtals.i  --  crystal structures for parameter fitting  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     e0_lattice   ideal lattice energy for the current crystal */
/*     moment_0     ideal dipole moment for monomer from crystal */
/*     nxtal        number of crystal structures to be stored */
/*     nvary        number of potential parameters to optimize */
/*     ivary        index for the types of potential parameters */
/*     vary         atom numbers involved in potential parameters */
/*     iresid       crystal structure to which each residual refers */
/*     rsdtyp       experimental variable for each of the residuals */
/*     vartyp       type of potential parameter to be optimized */




/*     write the values of parameters and scaled derivatives */

    /* Parameter adjustments */
    --f;
    --gs;
    --xx;

    /* Function Body */
    io___139.ciunit = iounit_1.iout;
    s_wsfe(&io___139);
    do_fio(&c__1, (char *)&(*niter), (ftnlen)sizeof(integer));
    e_wsfe();
    i__1 = xtal1_1.nvary;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (xtal1_1.ivary[i__ - 1] <= 3) {
	    io___141.ciunit = iounit_1.iout;
	    s_wsfe(&io___141);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, vartyp_ref(0, i__), (ftnlen)20);
	    do_fio(&c__1, (char *)&vary_ref(1, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xx[i__], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gs[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (xtal1_1.ivary[i__ - 1] == 4 || xtal1_1.ivary[i__ - 1] == 5)
		 {
	    io___142.ciunit = iounit_1.iout;
	    s_wsfe(&io___142);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, vartyp_ref(0, i__), (ftnlen)20);
	    do_fio(&c__1, (char *)&vary_ref(1, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xx[i__], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gs[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (xtal1_1.ivary[i__ - 1] == 6) {
	    io___143.ciunit = iounit_1.iout;
	    s_wsfe(&io___143);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, vartyp_ref(0, i__), (ftnlen)20);
	    do_fio(&c__1, (char *)&vary_ref(1, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&vary_ref(2, i__), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xx[i__], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&gs[i__], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     write the values of the residual functions */

    io___144.ciunit = iounit_1.iout;
    s_wsfe(&io___144);
    do_fio(&c__1, (char *)&(*niter), (ftnlen)sizeof(integer));
    e_wsfe();
    i__1 = *nresid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___145.ciunit = iounit_1.iout;
	s_wsfe(&io___145);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, rsdtyp_ref(0, i__), (ftnlen)20);
	do_fio(&c__1, (char *)&xtal1_1.iresid[i__ - 1], (ftnlen)sizeof(
		integer));
	do_fio(&c__1, (char *)&f[i__], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    io___146.ciunit = iounit_1.iout;
    s_wsfe(&io___146);
    e_wsfe();
    return 0;
} /* xtalwrt_ */

#undef vartyp_ref
#undef rsdtyp_ref
#undef vary_ref


/* Main program alias */ int xtalfit_ () { MAIN__ (); return 0; }
