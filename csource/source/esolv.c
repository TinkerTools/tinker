/* esolv.f -- translated by f2c (version 20050501).
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
    doublereal esum, eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, 
	    ebt, ett, ev, ec, ecd, ed, em, ep, er, es, elf, eg, ex;
} energi_;

#define energi_1 energi_

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
    doublereal m2[25000], deform, difft, diffv, diffc;
    logical use_smooth__, use_dem__, use_gda__, use_tophat__, use_stophat__;
} warp_;

#define warp_1 warp_

struct {
    doublereal pchg[25000];
    integer nion, iion[25000], jion[25000], kion[25000], chglist[25000];
} charge_;

#define charge_1 charge_

struct {
    doublereal electric, dielec, ebuffer, c2scale, c3scale, c4scale, c5scale;
    logical neutnbr, neutcut;
} chgpot_;

#define chgpot_1 chgpot_

struct {
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    doublereal off, off2, cut, cut2, c0, c1, c2, c3, c4, c5, f0, f1, f2, f3, 
	    f4, f5, f6, f7;
} shunt_;

#define shunt_1 shunt_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal gkr[5000], gkc;
} gk_;

#define gk_1 gk_

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
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

struct {
    doublereal rad[5000], eps[5000], rad4[5000], eps4[5000], reduct[5000];
} kvdws_;

#define kvdws_1 kvdws_

struct {
    doublereal solvprs, surften, spcut, spoff, stcut, stoff, rcav[25000], 
	    rdisp[25000], cdisp[25000];
} npolar_;

#define npolar_1 npolar_

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
    doublereal radmin[1000000]	/* was [1000][1000] */, epsilon[1000000]	
	    /* was [1000][1000] */, radmin4[1000000]	/* was [1000][1000] */
	    , epsilon4[1000000]	/* was [1000][1000] */, radhbnd[1000000]	
	    /* was [1000][1000] */, epshbnd[1000000]	/* was [1000][1000] */
	    , kred[25000];
    integer ired[25000], nvdw, ivdw[25000], jvdw[25000], nvt, ivt[25000], jvt[
	    25000];
} vdw_;

#define vdw_1 vdw_

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################################ */
/*     ##         COPYRIGHT (C)  1993  by  Jay William Ponder        ## */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ######################################################## */
/*     ##                                                    ## */
/*     ##  subroutine esolv  --  continuum solvation energy  ## */
/*     ##                                                    ## */
/*     ######################################################## */


/*     "esolv" calculates the continuum solvation energy for */
/*     surface area, generalized Born, generalized Kirkwood */
/*     and Poisson-Boltzmann solvation models */


/* Subroutine */ int esolv_(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal e;
    static integer i__;
    static doublereal ai, rb, ri, aes[25000];
    extern /* Subroutine */ int egk_(void), epb_(void), enp_(doublereal *, 
	    doublereal *);
    static doublereal ecav, term;
    extern /* Subroutine */ int egb0a_(void), egb0b_(void);
    static doublereal edisp, probe;
    extern /* Subroutine */ int surface_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);



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
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




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




/*     zero out the continuum solvation energy */

    energi_1.es = 0.;

/*     set a value for the solvent molecule probe radius */

    probe = 1.4;

/*     total solvation energy for surface area only models */

    if (s_cmp(solute_1.solvtyp, "ASP", (ftnlen)8, (ftnlen)3) == 0 || s_cmp(
	    solute_1.solvtyp, "SASA", (ftnlen)8, (ftnlen)4) == 0) {
	surface_(&energi_1.es, aes, solute_1.rsolv, solute_1.asolv, &probe);

/*     nonpolar energy for Onion GB method via exact area */

    } else if (s_cmp(solute_1.solvtyp, "ONION", (ftnlen)8, (ftnlen)5) == 0) {
	surface_(&energi_1.es, aes, solute_1.rsolv, solute_1.asolv, &probe);

/*     nonpolar energy as cavity formation plus dispersion */

    } else if (s_cmp(solute_1.solvtyp, "GK", (ftnlen)8, (ftnlen)2) == 0 || 
	    s_cmp(solute_1.solvtyp, "PB", (ftnlen)8, (ftnlen)2) == 0) {
	enp_(&ecav, &edisp);
	energi_1.es = ecav + edisp;

/*     nonpolar energy for GB via ACE surface area approximation */

    } else {
	term = 12.566370614359172;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ai = solute_1.asolv[i__ - 1];
	    ri = solute_1.rsolv[i__ - 1];
	    rb = solute_1.rborn[i__ - 1];
	    if (rb != 0.) {
/* Computing 2nd power */
		d__1 = ri + probe;
/* Computing 6th power */
		d__2 = ri / rb, d__2 *= d__2;
		e = ai * term * (d__1 * d__1) * (d__2 * (d__2 * d__2));
		energi_1.es += e;
	    }
	}
    }

/*     get polarization energy term for the solvation methods */

    if (s_cmp(solute_1.solvtyp, "GK", (ftnlen)8, (ftnlen)2) == 0) {
	egk_();
    } else if (s_cmp(solute_1.solvtyp, "PB", (ftnlen)8, (ftnlen)2) == 0) {
	epb_();
    } else if (potent_1.use_born__) {
	if (warp_1.use_smooth__) {
	    egb0b_();
	} else {
	    egb0a_();
	}
    }
    return 0;
} /* esolv_ */



/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine egb0a  --  generalized Born polarization energy  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "egb0a" calculates the generalized Born polarization energy */
/*     for the GB/SA solvation models */


/* Subroutine */ int egb0a_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, k;
    static doublereal r__, r2, r3, r4, r5, r6, r7, fi;
    static integer ii, kk;
    static doublereal xi, yi, zi, xr, yr, zr, rb2, rm2, fgb, fik, fgm, fgrp;
    static logical usei;
    static doublereal taper, shift, trans, dwater;
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *);
    static logical proceed;



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
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  shunt.i  --  polynomial switching function coefficients  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     off    distance at which the potential energy goes to zero */
/*     off2   square of distance at which the potential goes to zero */
/*     cut    distance at which switching of the potential begins */
/*     cut2   square of distance at which the switching begins */
/*     c0     zeroth order coefficient of multiplicative switch */
/*     c1     first order coefficient of multiplicative switch */
/*     c2     second order coefficient of multiplicative switch */
/*     c3     third order coefficient of multiplicative switch */
/*     c4     fourth order coefficient of multiplicative switch */
/*     c5     fifth order coefficient of multiplicative switch */
/*     f0     zeroth order coefficient of additive switch function */
/*     f1     first order coefficient of additive switch function */
/*     f2     second order coefficient of additive switch function */
/*     f3     third order coefficient of additive switch function */
/*     f4     fourth order coefficient of additive switch function */
/*     f5     fifth order coefficient of additive switch function */
/*     f6     sixth order coefficient of additive switch function */
/*     f7     seventh order coefficient of additive switch function */




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




/*     set the solvent dielectric and energy conversion factor */

    if (charge_1.nion == 0) {
	return 0;
    }
    dwater = 78.3;
    f = -chgpot_1.electric * (1. - 1. / dwater);

/*     set cutoff distances and switching function coefficients */

    switch_("CHARGE", (ftnlen)6);

/*     calculate GB electrostatic polarization energy term */

    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = charge_1.iion[ii - 1];
	usei = usage_1.use[i__ - 1];
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	fi = f * charge_1.pchg[ii - 1];

/*     decide whether to compute the current interaction */

	i__2 = charge_1.nion;
	for (kk = ii; kk <= i__2; ++kk) {
	    k = charge_1.iion[kk - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xr = xi - atoms_1.x[k - 1];
		yr = yi - atoms_1.y[k - 1];
		zr = zi - atoms_1.z__[k - 1];
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    fik = fi * charge_1.pchg[kk - 1];
		    rb2 = solute_1.rborn[i__ - 1] * solute_1.rborn[k - 1];
		    fgb = sqrt(r2 + rb2 * exp(r2 * -.25 / rb2));
		    e = fik / fgb;

/*     use shifted energy switching if near the cutoff distance */

/* Computing 2nd power */
		    d__1 = (shunt_1.off + shunt_1.cut) * .5;
		    rm2 = d__1 * d__1;
		    fgm = sqrt(rm2 + rb2 * exp(rm2 * -.25 / rb2));
		    shift = fik / fgm;
		    e -= shift;
		    if (r2 > shunt_1.cut2) {
			r__ = sqrt(r2);
			r3 = r2 * r__;
			r4 = r2 * r2;
			r5 = r2 * r3;
			r6 = r3 * r3;
			r7 = r3 * r4;
			taper = shunt_1.c5 * r5 + shunt_1.c4 * r4 + 
				shunt_1.c3 * r3 + shunt_1.c2 * r2 + 
				shunt_1.c1 * r__ + shunt_1.c0;
			trans = fik * (shunt_1.f7 * r7 + shunt_1.f6 * r6 + 
				shunt_1.f5 * r5 + shunt_1.f4 * r4 + 
				shunt_1.f3 * r3 + shunt_1.f2 * r2 + 
				shunt_1.f1 * r__ + shunt_1.f0);
			e = e * taper + trans;
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
		    }

/*     increment the overall GB solvation energy component */

		    if (i__ == k) {
			e *= .5;
		    }
		    energi_1.es += e;
		}
	    }
	}
    }
    return 0;
} /* egb0a_ */



/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine egb0b  --  GB polarization energy for smoothing  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "egb0b" calculates the generalized Born polarization energy */
/*     for the GB/SA solvation models for use with potential smoothing */
/*     methods via analogy to the smoothing of Coulomb's law */


/* Subroutine */ int egb0b_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, k;
    static doublereal r2, fi;
    static integer ii, kk;
    static doublereal xi, yi, zi, xr, yr, zr, rb2, fgb, fik;
    extern doublereal erf_(doublereal *);
    static doublereal fgrp;
    static logical usei;
    static doublereal width, sterm, dwater;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;



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
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




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




/*     set the solvent dielectric and energy conversion factor */

    if (charge_1.nion == 0) {
	return 0;
    }
    dwater = 78.3;
    f = -chgpot_1.electric * (1. - 1. / dwater);

/*     set the extent of smoothing to be performed */

    sterm = .5 / sqrt(warp_1.diffc);

/*     calculate GB electrostatic polarization energy term */

    i__1 = charge_1.nion;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = charge_1.iion[ii - 1];
	usei = usage_1.use[i__ - 1];
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	fi = f * charge_1.pchg[ii - 1];

/*     decide whether to compute the current interaction */

	i__2 = charge_1.nion;
	for (kk = ii; kk <= i__2; ++kk) {
	    k = charge_1.iion[kk - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xr = xi - atoms_1.x[k - 1];
		yr = yi - atoms_1.y[k - 1];
		zr = zi - atoms_1.z__[k - 1];
		r2 = xr * xr + yr * yr + zr * zr;
		fik = fi * charge_1.pchg[kk - 1];
		rb2 = solute_1.rborn[i__ - 1] * solute_1.rborn[k - 1];
		fgb = sqrt(r2 + rb2 * exp(r2 * -.25 / rb2));
		e = fik / fgb;

/*     use a smoothable GB analogous to Coulomb's law solution */

		if (warp_1.deform > 0.) {
		    width = warp_1.deform + rb2 * .15 * exp(rb2 * -.006 / 
			    warp_1.deform);
		    width = sterm / sqrt(width);
		    d__1 = width * fgb;
		    e *= erf_(&d__1);
		}

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    e *= fgrp;
		}

/*     increment the overall GB solvation energy component */

		if (i__ == k) {
		    e *= .5;
		}
		energi_1.es += e;
	    }
	}
    }
    return 0;
} /* egb0b_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine egk  --  generalized Kirkwood solvation model  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "egk" calculates the generalized Kirkwood electrostatic */
/*     solvation free energy for the GK/NP implicit solvation model */


/* Subroutine */ int egk_(void)
{
    extern /* Subroutine */ int egk0a_(void), ediff_(void);



/*     compute the generalized Kirkwood electrostatic energy */

    egk0a_();

/*     correct the solvation energy for vacuum to polarized state */

    ediff_();
    return 0;
} /* egk_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine egk0a  --  find generalized Kirkwood energy  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "egk0a" calculates the electrostatic portion of the continuum */
/*     solvation energy via the generalized Kirkwood model */


/* Subroutine */ int egk0a_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal a[15]	/* was [5][3] */, e;
    static integer i__, k;
    static doublereal r2, expcdexpc, fc, fd, ci, ei, ck, gf, gc[10];
    static integer ii, kk;
    static doublereal fq, xi, yi, zi, xr, yr, zr, gf2, gf3, gf5, gf7, rb2, 
	    gf9, xr2, yr2, zr2, rbi, rbk, dxi, dyi, dzi, dxk, dyk, dzk, ewi, 
	    ewk, gux[10], guy[10], uxi, uyi, uzi, uxk, uyk, uzk, guz[10], 
	    ewii, fgrp, expc, ewki;
    static logical usei;
    static doublereal esym, gqxx[10], qxxi, qxyi, qxzi, qyyi, qyzi, qzzi, 
	    qxxk, qxyk, qxzk, qyyk, qyzk, qzzk, gqxy[10], gqxz[10], gqyy[10], 
	    gqyz[10], gqzz[10], expc1, dexpc, esymi, dwater;
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *);
    static logical proceed;
    extern /* Subroutine */ int chkpole_(void);
    static doublereal expterm;
    extern /* Subroutine */ int rotpole_(void);


#define a_ref(a_1,a_2) a[(a_2)*5 + a_1 - 0]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define uinds_ref(a_1,a_2) polar_1.uinds[(a_2)*3 + a_1 - 4]



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
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  gk.i  --  parameters for generalized Kirkwood solvation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     gkr      generalized Kirkwood cavity radii for atom types */
/*     gkc      tuning parameter exponent in the f(GB) function */




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
/*     ##  shunt.i  --  polynomial switching function coefficients  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     off    distance at which the potential energy goes to zero */
/*     off2   square of distance at which the potential goes to zero */
/*     cut    distance at which switching of the potential begins */
/*     cut2   square of distance at which the switching begins */
/*     c0     zeroth order coefficient of multiplicative switch */
/*     c1     first order coefficient of multiplicative switch */
/*     c2     second order coefficient of multiplicative switch */
/*     c3     third order coefficient of multiplicative switch */
/*     c4     fourth order coefficient of multiplicative switch */
/*     c5     fifth order coefficient of multiplicative switch */
/*     f0     zeroth order coefficient of additive switch function */
/*     f1     first order coefficient of additive switch function */
/*     f2     second order coefficient of additive switch function */
/*     f3     third order coefficient of additive switch function */
/*     f4     fourth order coefficient of additive switch function */
/*     f5     fifth order coefficient of additive switch function */
/*     f6     sixth order coefficient of additive switch function */
/*     f7     seventh order coefficient of additive switch function */




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




/*     set the bulk dielectric constant to the water value */

    if (mpole_1.npole == 0) {
	return 0;
    }
    dwater = 78.3;
    fc = chgpot_1.electric * 1. * (1. - dwater) / (dwater * 1. + 0.);
    fd = chgpot_1.electric * 2. * (1. - dwater) / (dwater * 2. + 1.);
    fq = chgpot_1.electric * 3. * (1. - dwater) / (dwater * 3. + 2.);

/*     set cutoff distances and switching function coefficients */

    switch_("MPOLE", (ftnlen)5);

/*     setup the multipoles for solvation only calculations */

    if (! potent_1.use_mpole__ && ! potent_1.use_polar__) {
	chkpole_();
	rotpole_();
    }

/*     calculate GK electrostatic solvation free energy */

    i__1 = mpole_1.npole;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = mpole_1.ipole[ii - 1];
	usei = usage_1.use[i__ - 1];
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	rbi = solute_1.rborn[i__ - 1];
	ci = rpole_ref(1, ii);
	uxi = rpole_ref(2, ii);
	uyi = rpole_ref(3, ii);
	uzi = rpole_ref(4, ii);
	qxxi = rpole_ref(5, ii);
	qxyi = rpole_ref(6, ii);
	qxzi = rpole_ref(7, ii);
	qyyi = rpole_ref(9, ii);
	qyzi = rpole_ref(10, ii);
	qzzi = rpole_ref(13, ii);
	dxi = uinds_ref(1, ii);
	dyi = uinds_ref(2, ii);
	dzi = uinds_ref(3, ii);

/*     decide whether to compute the current interaction */

	i__2 = mpole_1.npole;
	for (kk = ii; kk <= i__2; ++kk) {
	    k = mpole_1.ipole[kk - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xr = atoms_1.x[k - 1] - xi;
		yr = atoms_1.y[k - 1] - yi;
		zr = atoms_1.z__[k - 1] - zi;
		xr2 = xr * xr;
		yr2 = yr * yr;
		zr2 = zr * zr;
		r2 = xr2 + yr2 + zr2;
		if (r2 <= shunt_1.off2) {
		    rbk = solute_1.rborn[k - 1];
		    ck = rpole_ref(1, kk);
		    uxk = rpole_ref(2, kk);
		    uyk = rpole_ref(3, kk);
		    uzk = rpole_ref(4, kk);
		    qxxk = rpole_ref(5, kk);
		    qxyk = rpole_ref(6, kk);
		    qxzk = rpole_ref(7, kk);
		    qyyk = rpole_ref(9, kk);
		    qyzk = rpole_ref(10, kk);
		    qzzk = rpole_ref(13, kk);
		    dxk = uinds_ref(1, kk);
		    dyk = uinds_ref(2, kk);
		    dzk = uinds_ref(3, kk);
		    rb2 = rbi * rbk;
		    expterm = exp(-r2 / (gk_1.gkc * rb2));
		    expc = expterm / gk_1.gkc;
		    dexpc = -2. / (gk_1.gkc * rbi * rbk);
		    gf2 = 1. / (r2 + rb2 * expterm);
		    gf = sqrt(gf2);
		    gf3 = gf2 * gf;
		    gf5 = gf3 * gf2;
		    gf7 = gf5 * gf2;
		    gf9 = gf7 * gf2;

/*     reaction potential auxiliary terms */

		    a_ref(0, 0) = gf;
		    a_ref(1, 0) = -gf3;
		    a_ref(2, 0) = gf5 * 3.;
		    a_ref(3, 0) = gf7 * -15.;
		    a_ref(4, 0) = gf9 * 105.;

/*     reaction potential gradient auxiliary terms */

		    expc1 = 1. - expc;
		    a_ref(0, 1) = expc1 * a_ref(1, 0);
		    a_ref(1, 1) = expc1 * a_ref(2, 0);
		    a_ref(2, 1) = expc1 * a_ref(3, 0);
		    a_ref(3, 1) = expc1 * a_ref(4, 0);

/*     second reaction potential gradient auxiliary terms */

		    expcdexpc = -expc * dexpc;
		    a_ref(0, 2) = expc1 * a_ref(1, 1) + expcdexpc * a_ref(1, 
			    0);
		    a_ref(1, 2) = expc1 * a_ref(2, 1) + expcdexpc * a_ref(2, 
			    0);
		    a_ref(2, 2) = expc1 * a_ref(3, 1) + expcdexpc * a_ref(3, 
			    0);

/*     multiply the auxillary terms by their dieletric functions */

		    a_ref(0, 0) = fc * a_ref(0, 0);
		    a_ref(0, 1) = fc * a_ref(0, 1);
		    a_ref(0, 2) = fc * a_ref(0, 2);
		    a_ref(1, 0) = fd * a_ref(1, 0);
		    a_ref(1, 1) = fd * a_ref(1, 1);
		    a_ref(1, 2) = fd * a_ref(1, 2);
		    a_ref(2, 0) = fq * a_ref(2, 0);
		    a_ref(2, 1) = fq * a_ref(2, 1);
		    a_ref(2, 2) = fq * a_ref(2, 2);

/*     unweighted reaction potential tensor */

		    gc[0] = a_ref(0, 0);
		    gux[0] = xr * a_ref(1, 0);
		    guy[0] = yr * a_ref(1, 0);
		    guz[0] = zr * a_ref(1, 0);
		    gqxx[0] = xr2 * a_ref(2, 0);
		    gqyy[0] = yr2 * a_ref(2, 0);
		    gqzz[0] = zr2 * a_ref(2, 0);
		    gqxy[0] = xr * yr * a_ref(2, 0);
		    gqxz[0] = xr * zr * a_ref(2, 0);
		    gqyz[0] = yr * zr * a_ref(2, 0);

/*     unweighted reaction potential gradient tensor */

		    gc[1] = xr * a_ref(0, 1);
		    gc[2] = yr * a_ref(0, 1);
		    gc[3] = zr * a_ref(0, 1);
		    gux[1] = a_ref(1, 0) + xr2 * a_ref(1, 1);
		    gux[2] = xr * yr * a_ref(1, 1);
		    gux[3] = xr * zr * a_ref(1, 1);
		    guy[1] = gux[2];
		    guy[2] = a_ref(1, 0) + yr2 * a_ref(1, 1);
		    guy[3] = yr * zr * a_ref(1, 1);
		    guz[1] = gux[3];
		    guz[2] = guy[3];
		    guz[3] = a_ref(1, 0) + zr2 * a_ref(1, 1);
		    gqxx[1] = xr * (a_ref(2, 0) * 2. + xr2 * a_ref(2, 1));
		    gqxx[2] = yr * xr2 * a_ref(2, 1);
		    gqxx[3] = zr * xr2 * a_ref(2, 1);
		    gqyy[1] = xr * yr2 * a_ref(2, 1);
		    gqyy[2] = yr * (a_ref(2, 0) * 2. + yr2 * a_ref(2, 1));
		    gqyy[3] = zr * yr2 * a_ref(2, 1);
		    gqzz[1] = xr * zr2 * a_ref(2, 1);
		    gqzz[2] = yr * zr2 * a_ref(2, 1);
		    gqzz[3] = zr * (a_ref(2, 0) * 2. + zr2 * a_ref(2, 1));
		    gqxy[1] = yr * (a_ref(2, 0) + xr2 * a_ref(2, 1));
		    gqxy[2] = xr * (a_ref(2, 0) + yr2 * a_ref(2, 1));
		    gqxy[3] = zr * xr * yr * a_ref(2, 1);
		    gqxz[1] = zr * (a_ref(2, 0) + xr2 * a_ref(2, 1));
		    gqxz[2] = gqxy[3];
		    gqxz[3] = xr * (a_ref(2, 0) + zr2 * a_ref(2, 1));
		    gqyz[1] = gqxy[3];
		    gqyz[2] = zr * (a_ref(2, 0) + yr2 * a_ref(2, 1));
		    gqyz[3] = yr * (a_ref(2, 0) + zr2 * a_ref(2, 1));

/*     unweighted second reaction potential gradient tensor */

		    gc[4] = a_ref(0, 1) + xr2 * a_ref(0, 2);
		    gc[5] = xr * yr * a_ref(0, 2);
		    gc[6] = xr * zr * a_ref(0, 2);
		    gc[7] = a_ref(0, 1) + yr2 * a_ref(0, 2);
		    gc[8] = yr * zr * a_ref(0, 2);
		    gc[9] = a_ref(0, 1) + zr2 * a_ref(0, 2);
		    gux[4] = xr * (a_ref(1, 1) + a_ref(1, 1) * 2. + xr2 * 
			    a_ref(1, 2));
		    gux[5] = yr * (a_ref(1, 1) + xr2 * a_ref(1, 2));
		    gux[6] = zr * (a_ref(1, 1) + xr2 * a_ref(1, 2));
		    gux[7] = xr * (a_ref(1, 1) + yr2 * a_ref(1, 2));
		    gux[8] = zr * xr * yr * a_ref(1, 2);
		    gux[9] = xr * (a_ref(1, 1) + zr2 * a_ref(1, 2));
		    guy[4] = yr * (a_ref(1, 1) + xr2 * a_ref(1, 2));
		    guy[5] = xr * (a_ref(1, 1) + yr2 * a_ref(1, 2));
		    guy[6] = gux[8];
		    guy[7] = yr * (a_ref(1, 1) + a_ref(1, 1) * 2. + yr2 * 
			    a_ref(1, 2));
		    guy[8] = zr * (a_ref(1, 1) + yr2 * a_ref(1, 2));
		    guy[9] = yr * (a_ref(1, 1) + zr2 * a_ref(1, 2));
		    guz[4] = zr * (a_ref(1, 1) + xr2 * a_ref(1, 2));
		    guz[5] = gux[8];
		    guz[6] = xr * (a_ref(1, 1) + zr2 * a_ref(1, 2));
		    guz[7] = zr * (a_ref(1, 1) + yr2 * a_ref(1, 2));
		    guz[8] = yr * (a_ref(1, 1) + zr2 * a_ref(1, 2));
		    guz[9] = zr * (a_ref(1, 1) + a_ref(1, 1) * 2. + zr2 * 
			    a_ref(1, 2));
		    gqxx[4] = a_ref(2, 0) * 2. + xr2 * (a_ref(2, 1) * 5. + 
			    xr2 * a_ref(2, 2));
		    gqxx[5] = yr * xr * (a_ref(2, 1) * 2. + xr2 * a_ref(2, 2))
			    ;
		    gqxx[6] = zr * xr * (a_ref(2, 1) * 2. + xr2 * a_ref(2, 2))
			    ;
		    gqxx[7] = xr2 * (a_ref(2, 1) + yr2 * a_ref(2, 2));
		    gqxx[8] = zr * yr * xr2 * a_ref(2, 2);
		    gqxx[9] = xr2 * (a_ref(2, 1) + zr2 * a_ref(2, 2));
		    gqyy[4] = yr2 * (a_ref(2, 1) + xr2 * a_ref(2, 2));
		    gqyy[5] = xr * yr * (a_ref(2, 1) * 2. + yr2 * a_ref(2, 2))
			    ;
		    gqyy[6] = xr * zr * yr2 * a_ref(2, 2);
		    gqyy[7] = a_ref(2, 0) * 2. + yr2 * (a_ref(2, 1) * 5. + 
			    yr2 * a_ref(2, 2));
		    gqyy[8] = yr * zr * (a_ref(2, 1) * 2. + yr2 * a_ref(2, 2))
			    ;
		    gqyy[9] = yr2 * (a_ref(2, 1) + zr2 * a_ref(2, 2));
		    gqzz[4] = zr2 * (a_ref(2, 1) + xr2 * a_ref(2, 2));
		    gqzz[5] = xr * yr * zr2 * a_ref(2, 2);
		    gqzz[6] = xr * zr * (a_ref(2, 1) * 2. + zr2 * a_ref(2, 2))
			    ;
		    gqzz[7] = zr2 * (a_ref(2, 1) + yr2 * a_ref(2, 2));
		    gqzz[8] = yr * zr * (a_ref(2, 1) * 2. + zr2 * a_ref(2, 2))
			    ;
		    gqzz[9] = a_ref(2, 0) * 2. + zr2 * (a_ref(2, 1) * 5. + 
			    zr2 * a_ref(2, 2));
		    gqxy[4] = xr * yr * (a_ref(2, 1) * 3. + xr2 * a_ref(2, 2))
			    ;
		    gqxy[5] = a_ref(2, 0) + (xr2 + yr2) * a_ref(2, 1) + xr2 * 
			    yr2 * a_ref(2, 2);
		    gqxy[6] = zr * yr * (a_ref(2, 1) + xr2 * a_ref(2, 2));
		    gqxy[7] = xr * yr * (a_ref(2, 1) * 3. + yr2 * a_ref(2, 2))
			    ;
		    gqxy[8] = zr * xr * (a_ref(2, 1) + yr2 * a_ref(2, 2));
		    gqxy[9] = xr * yr * (a_ref(2, 1) + zr2 * a_ref(2, 2));
		    gqxz[4] = xr * zr * (a_ref(2, 1) * 3. + xr2 * a_ref(2, 2))
			    ;
		    gqxz[5] = yr * zr * (a_ref(2, 1) + xr2 * a_ref(2, 2));
		    gqxz[6] = a_ref(2, 0) + (xr2 + zr2) * a_ref(2, 1) + xr2 * 
			    zr2 * a_ref(2, 2);
		    gqxz[7] = xr * zr * (a_ref(2, 1) + yr2 * a_ref(2, 2));
		    gqxz[8] = xr * yr * (a_ref(2, 1) + zr2 * a_ref(2, 2));
		    gqxz[9] = xr * zr * (a_ref(2, 1) * 3. + zr2 * a_ref(2, 2))
			    ;
		    gqyz[4] = zr * yr * (a_ref(2, 1) + xr2 * a_ref(2, 2));
		    gqyz[5] = xr * zr * (a_ref(2, 1) + yr2 * a_ref(2, 2));
		    gqyz[6] = xr * yr * (a_ref(2, 1) + zr2 * a_ref(2, 2));
		    gqyz[7] = yr * zr * (a_ref(2, 1) * 3. + yr2 * a_ref(2, 2))
			    ;
		    gqyz[8] = a_ref(2, 0) + (yr2 + zr2) * a_ref(2, 1) + yr2 * 
			    zr2 * a_ref(2, 2);
		    gqyz[9] = yr * zr * (a_ref(2, 1) * 3. + zr2 * a_ref(2, 2))
			    ;

/*     electrostatic solvation free energy of the permanent multipoles */
/*     in their own GK reaction potential */

		    esym = ci * ck * gc[0] - uxi * (uxk * gux[1] + uyk * guy[
			    1] + uzk * guz[1]) - uyi * (uxk * gux[2] + uyk * 
			    guy[2] + uzk * guz[2]) - uzi * (uxk * gux[3] + 
			    uyk * guy[3] + uzk * guz[3]);
		    ewi = ci * (uxk * gc[1] + uyk * gc[2] + uzk * gc[3]) - ck 
			    * (uxi * gux[0] + uyi * guy[0] + uzi * guz[0]) + 
			    ci * (qxxk * gc[4] + qyyk * gc[7] + qzzk * gc[9] 
			    + (qxyk * gc[5] + qxzk * gc[6] + qyzk * gc[8]) * 
			    2.) + ck * (qxxi * gqxx[0] + qyyi * gqyy[0] + 
			    qzzi * gqzz[0] + (qxyi * gqxy[0] + qxzi * gqxz[0] 
			    + qyzi * gqyz[0]) * 2.) - uxi * (qxxk * gux[4] + 
			    qyyk * gux[7] + qzzk * gux[9] + (qxyk * gux[5] + 
			    qxzk * gux[6] + qyzk * gux[8]) * 2.) - uyi * (
			    qxxk * guy[4] + qyyk * guy[7] + qzzk * guy[9] + (
			    qxyk * guy[5] + qxzk * guy[6] + qyzk * guy[8]) * 
			    2.) - uzi * (qxxk * guz[4] + qyyk * guz[7] + qzzk 
			    * guz[9] + (qxyk * guz[5] + qxzk * guz[6] + qyzk *
			     guz[8]) * 2.) + uxk * (qxxi * gqxx[1] + qyyi * 
			    gqyy[1] + qzzi * gqzz[1] + (qxyi * gqxy[1] + qxzi 
			    * gqxz[1] + qyzi * gqyz[1]) * 2.) + uyk * (qxxi * 
			    gqxx[2] + qyyi * gqyy[2] + qzzi * gqzz[2] + (qxyi 
			    * gqxy[2] + qxzi * gqxz[2] + qyzi * gqyz[2]) * 2.)
			     + uzk * (qxxi * gqxx[3] + qyyi * gqyy[3] + qzzi *
			     gqzz[3] + (qxyi * gqxy[3] + qxzi * gqxz[3] + 
			    qyzi * gqyz[3]) * 2.) + qxxi * (qxxk * gqxx[4] + 
			    qyyk * gqxx[7] + qzzk * gqxx[9] + (qxyk * gqxx[5] 
			    + qxzk * gqxx[6] + qyzk * gqxx[8]) * 2.) + qyyi * 
			    (qxxk * gqyy[4] + qyyk * gqyy[7] + qzzk * gqyy[9] 
			    + (qxyk * gqyy[5] + qxzk * gqyy[6] + qyzk * gqyy[
			    8]) * 2.) + qzzi * (qxxk * gqzz[4] + qyyk * gqzz[
			    7] + qzzk * gqzz[9] + (qxyk * gqzz[5] + qxzk * 
			    gqzz[6] + qyzk * gqzz[8]) * 2.) + (qxyi * (qxxk * 
			    gqxy[4] + qyyk * gqxy[7] + qzzk * gqxy[9] + (qxyk 
			    * gqxy[5] + qxzk * gqxy[6] + qyzk * gqxy[8]) * 2.)
			     + qxzi * (qxxk * gqxz[4] + qyyk * gqxz[7] + qzzk 
			    * gqxz[9] + (qxyk * gqxz[5] + qxzk * gqxz[6] + 
			    qyzk * gqxz[8]) * 2.) + qyzi * (qxxk * gqyz[4] + 
			    qyyk * gqyz[7] + qzzk * gqyz[9] + (qxyk * gqyz[5] 
			    + qxzk * gqyz[6] + qyzk * gqyz[8]) * 2.)) * 2.;
		    ewk = ci * (uxk * gux[0] + uyk * guy[0] + uzk * guz[0]) - 
			    ck * (uxi * gc[1] + uyi * gc[2] + uzi * gc[3]) + 
			    ci * (qxxk * gqxx[0] + qyyk * gqyy[0] + qzzk * 
			    gqzz[0] + (qxyk * gqxy[0] + qxzk * gqxz[0] + qyzk 
			    * gqyz[0]) * 2.) + ck * (qxxi * gc[4] + qyyi * gc[
			    7] + qzzi * gc[9] + (qxyi * gc[5] + qxzi * gc[6] 
			    + qyzi * gc[8]) * 2.) - uxi * (qxxk * gqxx[1] + 
			    qyyk * gqyy[1] + qzzk * gqzz[1] + (qxyk * gqxy[1] 
			    + qxzk * gqxz[1] + qyzk * gqyz[1]) * 2.) - uyi * (
			    qxxk * gqxx[2] + qyyk * gqyy[2] + qzzk * gqzz[2] 
			    + (qxyk * gqxy[2] + qxzk * gqxz[2] + qyzk * gqyz[
			    2]) * 2.) - uzi * (qxxk * gqxx[3] + qyyk * gqyy[3]
			     + qzzk * gqzz[3] + (qxyk * gqxy[3] + qxzk * gqxz[
			    3] + qyzk * gqyz[3]) * 2.) + uxk * (qxxi * gux[4] 
			    + qyyi * gux[7] + qzzi * gux[9] + (qxyi * gux[5] 
			    + qxzi * gux[6] + qyzi * gux[8]) * 2.) + uyk * (
			    qxxi * guy[4] + qyyi * guy[7] + qzzi * guy[9] + (
			    qxyi * guy[5] + qxzi * guy[6] + qyzi * guy[8]) * 
			    2.) + uzk * (qxxi * guz[4] + qyyi * guz[7] + qzzi 
			    * guz[9] + (qxyi * guz[5] + qxzi * guz[6] + qyzi *
			     guz[8]) * 2.) + qxxi * (qxxk * gqxx[4] + qyyk * 
			    gqyy[4] + qzzk * gqzz[4] + (qxyk * gqxy[4] + qxzk 
			    * gqxz[4] + qyzk * gqyz[4]) * 2.) + qyyi * (qxxk *
			     gqxx[7] + qyyk * gqyy[7] + qzzk * gqzz[7] + (
			    qxyk * gqxy[7] + qxzk * gqxz[7] + qyzk * gqyz[7]) 
			    * 2.) + qzzi * (qxxk * gqxx[9] + qyyk * gqyy[9] + 
			    qzzk * gqzz[9] + (qxyk * gqxy[9] + qxzk * gqxz[9] 
			    + qyzk * gqyz[9]) * 2.) + (qxyi * (qxxk * gqxx[5] 
			    + qyyk * gqyy[5] + qzzk * gqzz[5] + (qxyk * gqxy[
			    5] + qxzk * gqxz[5] + qyzk * gqyz[5]) * 2.) + 
			    qxzi * (qxxk * gqxx[6] + qyyk * gqyy[6] + qzzk * 
			    gqzz[6] + (qxyk * gqxy[6] + qxzk * gqxz[6] + qyzk 
			    * gqyz[6]) * 2.) + qyzi * (qxxk * gqxx[8] + qyyk *
			     gqyy[8] + qzzk * gqzz[8] + (qxyk * gqxy[8] + 
			    qxzk * gqxz[8] + qyzk * gqyz[8]) * 2.)) * 2.;

/*     electrostatic solvation free energy of the permenant multipoles */
/*     in the GK reaction potential of the induced dipoles */

		    esymi = -uxi * (dxk * gux[1] + dyk * guy[1] + dzk * guz[1]
			    ) - uyi * (dxk * gux[2] + dyk * guy[2] + dzk * 
			    guz[2]) - uzi * (dxk * gux[3] + dyk * guy[3] + 
			    dzk * guz[3]) - uxk * (dxi * gux[1] + dyi * guy[1]
			     + dzi * guz[1]) - uyk * (dxi * gux[2] + dyi * 
			    guy[2] + dzi * guz[2]) - uzk * (dxi * gux[3] + 
			    dyi * guy[3] + dzi * guz[3]);
		    ewii = ci * (dxk * gc[1] + dyk * gc[2] + dzk * gc[3]) - 
			    ck * (dxi * gux[0] + dyi * guy[0] + dzi * guz[0]) 
			    - dxi * (qxxk * gux[4] + qyyk * gux[7] + qzzk * 
			    gux[9] + (qxyk * gux[5] + qxzk * gux[6] + qyzk * 
			    gux[8]) * 2.) - dyi * (qxxk * guy[4] + qyyk * guy[
			    7] + qzzk * guy[9] + (qxyk * guy[5] + qxzk * guy[
			    6] + qyzk * guy[8]) * 2.) - dzi * (qxxk * guz[4] 
			    + qyyk * guz[7] + qzzk * guz[9] + (qxyk * guz[5] 
			    + qxzk * guz[6] + qyzk * guz[8]) * 2.) + dxk * (
			    qxxi * gqxx[1] + qyyi * gqyy[1] + qzzi * gqzz[1] 
			    + (qxyi * gqxy[1] + qxzi * gqxz[1] + qyzi * gqyz[
			    1]) * 2.) + dyk * (qxxi * gqxx[2] + qyyi * gqyy[2]
			     + qzzi * gqzz[2] + (qxyi * gqxy[2] + qxzi * gqxz[
			    2] + qyzi * gqyz[2]) * 2.) + dzk * (qxxi * gqxx[3]
			     + qyyi * gqyy[3] + qzzi * gqzz[3] + (qxyi * gqxy[
			    3] + qxzi * gqxz[3] + qyzi * gqyz[3]) * 2.);
		    ewki = ci * (dxk * gux[0] + dyk * guy[0] + dzk * guz[0]) 
			    - ck * (dxi * gc[1] + dyi * gc[2] + dzi * gc[3]) 
			    - dxi * (qxxk * gqxx[1] + qyyk * gqyy[1] + qzzk * 
			    gqzz[1] + (qxyk * gqxy[1] + qxzk * gqxz[1] + qyzk 
			    * gqyz[1]) * 2.) - dyi * (qxxk * gqxx[2] + qyyk * 
			    gqyy[2] + qzzk * gqzz[2] + (qxyk * gqxy[2] + qxzk 
			    * gqxz[2] + qyzk * gqyz[2]) * 2.) - dzi * (qxxk * 
			    gqxx[3] + qyyk * gqyy[3] + qzzk * gqzz[3] + (qxyk 
			    * gqxy[3] + qxzk * gqxz[3] + qyzk * gqyz[3]) * 2.)
			     + dxk * (qxxi * gux[4] + qyyi * gux[7] + qzzi * 
			    gux[9] + (qxyi * gux[5] + qxzi * gux[6] + qyzi * 
			    gux[8]) * 2.) + dyk * (qxxi * guy[4] + qyyi * guy[
			    7] + qzzi * guy[9] + (qxyi * guy[5] + qxzi * guy[
			    6] + qyzi * guy[8]) * 2.) + dzk * (qxxi * guz[4] 
			    + qyyi * guz[7] + qzzi * guz[9] + (qxyi * guz[5] 
			    + qxzi * guz[6] + qyzi * guz[8]) * 2.);

/*     total permanent and induced energies for this interaction */

		    e = esym + (ewi + ewk) * .5;
		    ei = (esymi + (ewii + ewki) * .5) * .5;

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
			ei *= fgrp;
		    }

/*     increment the total GK electrostatic solvation energy */

		    if (i__ == k) {
			e *= .5;
			ei *= .5;
		    }
		    energi_1.es = energi_1.es + e + ei;
		}
	    }
	}
    }
    return 0;
} /* egk0a_ */

#undef uinds_ref
#undef rpole_ref
#undef a_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine ediff  --  correction for vacuum to SCRF energy  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "ediff" calculates the energy of polarizing the vacuum induced */
/*     dipoles to their SCRF polarized values */


/* Subroutine */ int ediff_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal f;
    static integer i__, j, k;
    static doublereal r__, r2;
    static integer ii, kk;
    static doublereal ei, ci, ck;
    static integer ix;
    static doublereal sc[6];
    static integer iz, kx, kz;
    static doublereal xr, yr, zr, rr1, rr3, rr5, rr7, pdi, dix, diy, diz, pti,
	     dkx, dky, dkz, qix, qiy, qiz, uix, uiy, uiz, ukx, uky, ukz, qkx, 
	    qky, qkz, sci[8], gli[3], damp, fikp, fgrp;
    static logical usei, usek;
    static doublereal qixx, qixy, qixz, qiyy, qiyz, qizz, qkxx, qkxy, qkxz, 
	    qkyy, qkyz, qkzz;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal scale3, scale5, scale7, pgamma, pscale[25000];
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *);
    static logical proceed;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define uinds_ref(a_1,a_2) polar_1.uinds[(a_2)*3 + a_1 - 4]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  shunt.i  --  polynomial switching function coefficients  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     off    distance at which the potential energy goes to zero */
/*     off2   square of distance at which the potential goes to zero */
/*     cut    distance at which switching of the potential begins */
/*     cut2   square of distance at which the switching begins */
/*     c0     zeroth order coefficient of multiplicative switch */
/*     c1     first order coefficient of multiplicative switch */
/*     c2     second order coefficient of multiplicative switch */
/*     c3     third order coefficient of multiplicative switch */
/*     c4     fourth order coefficient of multiplicative switch */
/*     c5     fifth order coefficient of multiplicative switch */
/*     f0     zeroth order coefficient of additive switch function */
/*     f1     first order coefficient of additive switch function */
/*     f2     second order coefficient of additive switch function */
/*     f3     third order coefficient of additive switch function */
/*     f4     fourth order coefficient of additive switch function */
/*     f5     fifth order coefficient of additive switch function */
/*     f6     sixth order coefficient of additive switch function */
/*     f7     seventh order coefficient of additive switch function */




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




/*     set conversion factor, cutoff and scaling coefficients */

    if (mpole_1.npole == 0) {
	return 0;
    }
    f = chgpot_1.electric / chgpot_1.dielec;
    switch_("MPOLE", (ftnlen)5);

/*     set array needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pscale[i__ - 1] = 1.;
    }

/*     calculate the multipole interaction energy term */

    i__1 = mpole_1.npole - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	iz = mpole_1.zaxis[i__ - 1];
	ix = mpole_1.xaxis[i__ - 1];
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
	uix = uinds_ref(1, i__) - uind_ref(1, i__);
	uiy = uinds_ref(2, i__) - uind_ref(2, i__);
	uiz = uinds_ref(3, i__) - uind_ref(3, i__);
	usei = usage_1.use[ii - 1] || usage_1.use[iz - 1] || usage_1.use[ix - 
		1];
	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i14_ref(j, ii) - 1] = polpot_1.p4scale;
	    i__3 = polgrp_1.np11[ii - 1];
	    for (k = 1; k <= i__3; ++k) {
		if (i14_ref(j, ii) == ip11_ref(k, ii)) {
		    pscale[i14_ref(j, ii) - 1] = pscale[i14_ref(j, ii) - 1] * 
			    .5;
		}
	    }
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	}

/*     decide whether to compute the current interaction */

	i__2 = mpole_1.npole;
	for (k = i__ + 1; k <= i__2; ++k) {
	    kk = mpole_1.ipole[k - 1];
	    kz = mpole_1.zaxis[k - 1];
	    kx = mpole_1.xaxis[k - 1];
	    usek = usage_1.use[kk - 1] || usage_1.use[kz - 1] || usage_1.use[
		    kx - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &ii, &kk, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (! group_1.use_intra__) {
		proceed = TRUE_;
	    }
	    if (proceed) {
		proceed = usei || usek;
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
		yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
		zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
		image_(&xr, &yr, &zr);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
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
		    ukx = uinds_ref(1, k) - uind_ref(1, k);
		    uky = uinds_ref(2, k) - uind_ref(2, k);
		    ukz = uinds_ref(3, k) - uind_ref(3, k);

/*     construct some intermediate quadrupole values */

		    qix = qixx * xr + qixy * yr + qixz * zr;
		    qiy = qixy * xr + qiyy * yr + qiyz * zr;
		    qiz = qixz * xr + qiyz * yr + qizz * zr;
		    qkx = qkxx * xr + qkxy * yr + qkxz * zr;
		    qky = qkxy * xr + qkyy * yr + qkyz * zr;
		    qkz = qkxz * xr + qkyz * yr + qkzz * zr;

/*     calculate the scalar products for permanent multipoles */

		    sc[2] = dix * xr + diy * yr + diz * zr;
		    sc[3] = dkx * xr + dky * yr + dkz * zr;
		    sc[4] = qix * xr + qiy * yr + qiz * zr;
		    sc[5] = qkx * xr + qky * yr + qkz * zr;

/*     calculate the scalar products for polarization components */

		    sci[1] = uix * dkx + dix * ukx + uiy * dky + diy * uky + 
			    uiz * dkz + diz * ukz;
		    sci[2] = uix * xr + uiy * yr + uiz * zr;
		    sci[3] = ukx * xr + uky * yr + ukz * zr;
		    sci[6] = qix * ukx + qiy * uky + qiz * ukz;
		    sci[7] = qkx * uix + qky * uiy + qkz * uiz;

/*     calculate the gl functions for polarization components */

		    gli[0] = ck * sci[2] - ci * sci[3] + sci[1];
		    gli[1] = (sci[6] - sci[7]) * 2. - sci[2] * sc[3] - sc[2] *
			     sci[3];
		    gli[2] = sci[2] * sc[5] - sci[3] * sc[4];

/*     compute the energy contributions for this interaction */

		    rr1 = 1. / r__;
		    rr3 = rr1 / r2;
		    rr5 = rr3 * 3. / r2;
		    rr7 = rr5 * 5. / r2;
		    scale3 = 1.;
		    scale5 = 1.;
		    scale7 = 1.;
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
/* Computing 2nd power */
			    d__1 = damp;
			    scale7 = 1. - (1. - damp + d__1 * d__1 * .6) * 
				    exp(damp);
			}
		    }
		    ei = gli[0] * rr3 * scale3 + gli[1] * rr5 * scale5 + gli[
			    2] * rr7 * scale7;

/*     make the adjustment for scaled interactions */

		    fikp = f * pscale[kk - 1];
		    ei = fikp * .5 * ei;

/*     scale the interaction based on its group membership; */
/*     polarization cannot be group scaled as it is not pairwise */

		    if (group_1.use_group__) {
			ei *= fgrp;
		    }

/*     increment the total GK electrostatic solvation energy */

		    energi_1.es += ei;
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i12_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i13_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i14_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    pscale[i15_ref(j, ii) - 1] = 1.;
	}
    }
    return 0;
} /* ediff_ */

#undef uinds_ref
#undef rpole_ref
#undef uind_ref
#undef ip11_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine epb  --  Poisson-Boltzmann solvation model  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "epb" calculates the continuum solvation energy via the */
/*     Poisson-Boltzmann plus nonpolar implicit solvation */


/* Subroutine */ int epb_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int pbempole_(void);
    static doublereal e;
    static integer i__, ii;
    extern /* Subroutine */ int ediff_(void);


#define pbep_ref(a_1,a_2) pb_1.pbep[(a_2)*3 + a_1 - 4]
#define uinds_ref(a_1,a_2) polar_1.uinds[(a_2)*3 + a_1 - 4]



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
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




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




/*     compute the electrostatic energy via Poisson-Boltzmann */

    if (potent_1.use_polar__) {
	e = 0.;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    e = e + uinds_ref(1, i__) * pbep_ref(1, ii) + uinds_ref(2, i__) * 
		    pbep_ref(2, ii) + uinds_ref(3, i__) * pbep_ref(3, ii);
	}
	e = chgpot_1.electric * -.5 * e;
	pb_1.pbe += e;
    } else {
	pbempole_();
    }

/*     increment solvation energy by Poisson-Boltzmann results */

    energi_1.es += pb_1.pbe;

/*     correct the solvation energy for vacuum to polarized state */

    ediff_();
    return 0;
} /* epb_ */

#undef uinds_ref
#undef pbep_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine pbempole  --  permanent multipole PB energy  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "pbempole" calculates the permanent multipole PB energy, */
/*     field, forces and torques */


/* Subroutine */ int pbempole_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, ii;
    extern /* Subroutine */ int apbsempole_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal pos[75000]	/* was [3][25000] */, pbpole[325000]	
	    /* was [13][25000] */;


#define pos_ref(a_1,a_2) pos[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
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




/*     initialization of coordinates and multipoles for APBS */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pos_ref(1, i__) = atoms_1.x[i__ - 1];
	pos_ref(2, i__) = atoms_1.y[i__ - 1];
	pos_ref(3, i__) = atoms_1.z__[i__ - 1];
	for (j = 1; j <= 13; ++j) {
	    pbpole_ref(j, i__) = 0.;
	}
    }

/*     copy the permanent moments into an array for APBS */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	pbpole_ref(1, ii) = rpole_ref(1, i__);
	for (j = 2; j <= 4; ++j) {
	    pbpole_ref(j, ii) = rpole_ref(j, i__);
	}
	for (j = 5; j <= 13; ++j) {
	    pbpole_ref(j, ii) = rpole_ref(j, i__) * 3.;
	}
    }

/*     numerical solution of the Poisson-Boltzmann equation */

    apbsempole_(&atoms_1.n, pos, solute_1.rsolv, pbpole, &pb_1.pbe, pb_1.apbe,
	     pb_1.pbep, pb_1.pbfp, pb_1.pbtp);
    return 0;
} /* pbempole_ */

#undef pbpole_ref
#undef rpole_ref
#undef pos_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine enp  --  cavity/dispersion nonpolar solvation  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "enp" calculates the nonpolar continuum solvation energy */
/*     as a sum of cavity and dispersion terms */


/* Subroutine */ int enp_(doublereal *ecav, doublereal *edisp)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int ewca_(doublereal *);
    static doublereal reff, evol, reff2, reff3, reff4, reff5, taper, esurf, 
	    aesurf[25000];
    extern /* Subroutine */ int switch_(char *, ftnlen), volume_(doublereal *,
	     doublereal *, doublereal *), surface_(doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *);
    static doublereal exclude;



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
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




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




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  npolar.i  --  nonpolar cavity & dispersion parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     epso      water oxygen eps for implicit dispersion term */
/*     epsh      water hydrogen eps for implicit dispersion term */
/*     rmino     water oxygen Rmin for implicit dispersion term */
/*     rminh     water hydrogen Rmin for implicit dispersion term */
/*     awater    water number density at standard temp & pressure */
/*     slevy     enthalpy-to-free energy scale factor for dispersion */

/*     solvprs   limiting microscopic solvent pressure value */
/*     surften   limiting macroscopic surface tension value */
/*     spcut     starting radius for solvent pressure tapering */
/*     spoff     cutoff radius for solvent pressure tapering */
/*     stcut     starting radius for surface tension tapering */
/*     stoff     cutoff radius for surface tension tapering */
/*     rcav      atomic radius of each atom for cavitation energy */
/*     rdisp     atomic radius of each atom for dispersion energy */
/*     cdisp     maximum dispersion energy for each atom */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  shunt.i  --  polynomial switching function coefficients  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     off    distance at which the potential energy goes to zero */
/*     off2   square of distance at which the potential goes to zero */
/*     cut    distance at which switching of the potential begins */
/*     cut2   square of distance at which the switching begins */
/*     c0     zeroth order coefficient of multiplicative switch */
/*     c1     first order coefficient of multiplicative switch */
/*     c2     second order coefficient of multiplicative switch */
/*     c3     third order coefficient of multiplicative switch */
/*     c4     fourth order coefficient of multiplicative switch */
/*     c5     fifth order coefficient of multiplicative switch */
/*     f0     zeroth order coefficient of additive switch function */
/*     f1     first order coefficient of additive switch function */
/*     f2     second order coefficient of additive switch function */
/*     f3     third order coefficient of additive switch function */
/*     f4     fourth order coefficient of additive switch function */
/*     f5     fifth order coefficient of additive switch function */
/*     f6     sixth order coefficient of additive switch function */
/*     f7     seventh order coefficient of additive switch function */




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




/*     zero out the nonpolar continuum solvation energy terms */

    *ecav = 0.;
    *edisp = 0.;

/*     compute SASA and effective radius needed for cavity term */

    exclude = 1.4;
    surface_(&esurf, aesurf, npolar_1.rcav, solute_1.asolv, &exclude);
    reff = sqrt(esurf / (npolar_1.surften * 3.141592653589793238)) * .5;
    reff2 = reff * reff;
    reff3 = reff2 * reff;
    reff4 = reff3 * reff;
    reff5 = reff4 * reff;

/*     compute solvent excluded volume for needed for small solutes */

    if (reff < npolar_1.spoff) {
	volume_(&evol, npolar_1.rcav, &exclude);
	evol *= npolar_1.solvprs;
    }

/*     find cavity energy from only the solvent excluded volume */

    if (reff <= npolar_1.spcut) {
	*ecav = evol;

/*     find cavity energy from only a tapered volume term */

    } else if (reff > npolar_1.spcut && reff <= npolar_1.stoff) {
	switch_("GKV", (ftnlen)3);
	taper = shunt_1.c5 * reff5 + shunt_1.c4 * reff4 + shunt_1.c3 * reff3 
		+ shunt_1.c2 * reff2 + shunt_1.c1 * reff + shunt_1.c0;
	*ecav = taper * evol;

/*     find cavity energy using both volume and SASA terms */

    } else if (reff > npolar_1.stoff && reff <= npolar_1.spoff) {
	switch_("GKV", (ftnlen)3);
	taper = shunt_1.c5 * reff5 + shunt_1.c4 * reff4 + shunt_1.c3 * reff3 
		+ shunt_1.c2 * reff2 + shunt_1.c1 * reff + shunt_1.c0;
	*ecav = taper * evol;
	switch_("GKSA", (ftnlen)4);
	taper = shunt_1.c5 * reff5 + shunt_1.c4 * reff4 + shunt_1.c3 * reff3 
		+ shunt_1.c2 * reff2 + shunt_1.c1 * reff + shunt_1.c0;
	taper = 1. - taper;
	*ecav += taper * esurf;

/*     find cavity energy from only a tapered SASA term */

    } else if (reff > npolar_1.spoff && reff <= npolar_1.stcut) {
	switch_("GKSA", (ftnlen)4);
	taper = shunt_1.c5 * reff5 + shunt_1.c4 * reff4 + shunt_1.c3 * reff3 
		+ shunt_1.c2 * reff2 + shunt_1.c1 * reff + shunt_1.c0;
	taper = 1. - taper;
	*ecav = taper * esurf;

/*     find cavity energy from only a SASA-based term */

    } else {
	*ecav = esurf;
    }

/*     find the Weeks-Chandler-Andersen dispersion energy */

    ewca_(edisp);
/*     call ewcax (edisp) */
    return 0;
} /* enp_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine ewca  --  WCA dispersion energy for solvation  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "ewca" find the Weeks-Chandler-Andersen dispersion energy */
/*     of a solute using an HCT-like method */


/* Subroutine */ int ewca_(doublereal *edisp)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e;
    static integer i__, k;
    static doublereal r__, r2, ah, ao, ri, rk, sk, xi, yi, zi, xr, yr, zr, 
	    sk2, lik, uik, sum, lik2, lik3, lik4, lik5, uik2, uik3, uik4, 
	    uik5, iwca, lik10, lik11, lik12, uik10, uik11, irep, epsi, term, 
	    rmax, uik12, shctd, idisp, emixh, rmini, emixo, rmixh, rmixo, 
	    rmixh7, rmixo7, offset;



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  npolar.i  --  nonpolar cavity & dispersion parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     epso      water oxygen eps for implicit dispersion term */
/*     epsh      water hydrogen eps for implicit dispersion term */
/*     rmino     water oxygen Rmin for implicit dispersion term */
/*     rminh     water hydrogen Rmin for implicit dispersion term */
/*     awater    water number density at standard temp & pressure */
/*     slevy     enthalpy-to-free energy scale factor for dispersion */

/*     solvprs   limiting microscopic solvent pressure value */
/*     surften   limiting macroscopic surface tension value */
/*     spcut     starting radius for solvent pressure tapering */
/*     spoff     cutoff radius for solvent pressure tapering */
/*     stcut     starting radius for surface tension tapering */
/*     stoff     cutoff radius for surface tension tapering */
/*     rcav      atomic radius of each atom for cavitation energy */
/*     rdisp     atomic radius of each atom for dispersion energy */
/*     cdisp     maximum dispersion energy for each atom */




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
/*     nvt        number of distinct van der Waals types in the system */
/*     ivt        number of each distinct vdw type/class in the system */
/*     jvt        frequency of each vdw type or class in the system */




/*     zero out the Weeks-Chandler-Andersen dispersion energy */

    *edisp = 0.;

/*     set overlap scale factor for HCT descreening method */

    shctd = .81;
    offset = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	npolar_1.rdisp[i__ - 1] = kvdws_1.rad[atmtyp_1.class__[i__ - 1] - 1] 
		+ offset;
    }

/*     find the Weeks-Chandler-Andersen dispersion energy */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	epsi = kvdws_1.eps[atmtyp_1.class__[i__ - 1] - 1];
	rmini = kvdws_1.rad[atmtyp_1.class__[i__ - 1] - 1];
/* Computing 2nd power */
	d__1 = sqrt(.11) + sqrt(epsi);
	emixo = epsi * .44 / (d__1 * d__1);
/* Computing 3rd power */
	d__1 = rmini;
/* Computing 2nd power */
	d__2 = rmini;
	rmixo = (d__1 * (d__1 * d__1) + 4.9347068906249989) * 2. / (d__2 * 
		d__2 + 2.8985062499999996);
/* Computing 7th power */
	d__1 = rmixo, d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
	rmixo7 = d__2 * (d__1 * d__1);
	ao = emixo * rmixo7;
/* Computing 2nd power */
	d__1 = sqrt(.0135) + sqrt(epsi);
	emixh = epsi * .053999999999999999 / (d__1 * d__1);
/* Computing 3rd power */
	d__1 = rmini;
/* Computing 2nd power */
	d__2 = rmini;
	rmixh = (d__1 * (d__1 * d__1) + 2.3393951718749992) * 2. / (d__2 * 
		d__2 + 1.7622562499999996);
/* Computing 7th power */
	d__1 = rmixh, d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
	rmixh7 = d__2 * (d__1 * d__1);
	ah = emixh * rmixh7;
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	ri = npolar_1.rdisp[i__ - 1];

/*     remove  contribution due to solvent displaced by solute atoms */

	sum = 0.;
	i__2 = atoms_1.n;
	for (k = 1; k <= i__2; ++k) {
	    if (i__ != k) {
		xr = atoms_1.x[k - 1] - xi;
		yr = atoms_1.y[k - 1] - yi;
		zr = atoms_1.z__[k - 1] - zi;
		r2 = xr * xr + yr * yr + zr * zr;
		r__ = sqrt(r2);
		rk = npolar_1.rdisp[k - 1];
/*              sk = rk * shct(k) */
		sk = rk * shctd;
		sk2 = sk * sk;
		if (ri < r__ + sk) {
/* Computing MAX */
		    d__1 = ri, d__2 = r__ - sk;
		    rmax = max(d__1,d__2);
		    lik = rmax;
		    lik2 = lik * lik;
		    lik3 = lik2 * lik;
		    lik4 = lik3 * lik;
		    if (lik < rmixo) {
/* Computing MIN */
			d__1 = r__ + sk;
			uik = min(d__1,rmixo);
			uik2 = uik * uik;
			uik3 = uik2 * uik;
			uik4 = uik3 * uik;
			term = 12.566370614359172 / (r__ * 48.) * ((lik4 - 
				uik4) * 3. - r__ * 8. * (lik3 - uik3) + (r2 - 
				sk2) * 6. * (lik2 - uik2));
			iwca = -emixo * term;
			sum += iwca;
		    }
		    if (lik < rmixh) {
/* Computing MIN */
			d__1 = r__ + sk;
			uik = min(d__1,rmixh);
			uik2 = uik * uik;
			uik3 = uik2 * uik;
			uik4 = uik3 * uik;
			term = 12.566370614359172 / (r__ * 48.) * ((lik4 - 
				uik4) * 3. - r__ * 8. * (lik3 - uik3) + (r2 - 
				sk2) * 6. * (lik2 - uik2));
			iwca = emixh * -2. * term;
			sum += iwca;
		    }
		    uik = r__ + sk;
		    uik2 = uik * uik;
		    uik3 = uik2 * uik;
		    uik4 = uik3 * uik;
		    uik5 = uik4 * uik;
		    uik10 = uik5 * uik5;
		    uik11 = uik10 * uik;
		    uik12 = uik11 * uik;
		    if (uik > rmixo) {
			lik = max(rmax,rmixo);
			lik2 = lik * lik;
			lik3 = lik2 * lik;
			lik4 = lik3 * lik;
			lik5 = lik4 * lik;
			lik10 = lik5 * lik5;
			lik11 = lik10 * lik;
			lik12 = lik11 * lik;
			term = 12.566370614359172 / (r__ * 120. * lik5 * uik5)
				 * (uik * 15. * lik * r__ * (uik4 - lik4) - 
				uik2 * 10. * lik2 * (uik3 - lik3) + (sk2 - r2)
				 * 6. * (uik5 - lik5));
			idisp = ao * -2. * term;
			term = 12.566370614359172 / (r__ * 2640. * lik12 * 
				uik12) * (uik * 120. * lik * r__ * (uik11 - 
				lik11) - uik2 * 66. * lik2 * (uik10 - lik10) 
				+ (sk2 - r2) * 55. * (uik12 - lik12));
			irep = ao * rmixo7 * term;
			sum = sum + irep + idisp;
		    }
		    if (uik > rmixh) {
			lik = max(rmax,rmixh);
			lik2 = lik * lik;
			lik3 = lik2 * lik;
			lik4 = lik3 * lik;
			lik5 = lik4 * lik;
			lik10 = lik5 * lik5;
			lik11 = lik10 * lik;
			lik12 = lik11 * lik;
			term = 12.566370614359172 / (r__ * 120. * lik5 * uik5)
				 * (uik * 15. * lik * r__ * (uik4 - lik4) - 
				uik2 * 10. * lik2 * (uik3 - lik3) + (sk2 - r2)
				 * 6. * (uik5 - lik5));
			idisp = ah * -4. * term;
			term = 12.566370614359172 / (r__ * 2640. * lik12 * 
				uik12) * (uik * 120. * lik * r__ * (uik11 - 
				lik11) - uik2 * 66. * lik2 * (uik10 - lik10) 
				+ (sk2 - r2) * 55. * (uik12 - lik12));
			irep = ah * 2. * rmixh7 * term;
			sum = sum + irep + idisp;
		    }
		}
	    }
	}

/*     increment the overall dispersion energy component */

	e = npolar_1.cdisp[i__ - 1] - sum * .033427999999999999;
	*edisp += e;
    }
    return 0;
} /* ewca_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine ewcax  --  alternative WCA dispersion energy  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "ewcax" finds the Weeks-Chandler-Anderson dispersion energy */
/*     of a solute using a numerical "onion shell" method */


/* Subroutine */ int ewcax_(doublereal *edisp)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal fraction;
    extern /* Subroutine */ int surfatom_(integer *, doublereal *, doublereal 
	    *);
    static doublereal e;
    static integer i__, j, k;
    static doublereal t, her7, oer7, area, her14;
    static logical done;
    static doublereal oer14, roff[25000], epsi, rmax, shell, inner, ratio, 
	    rinit, tinit, rmini, epsoi, epshi, outer, rmult, offset, rminoi, 
	    rminhi, rswitch;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  npolar.i  --  nonpolar cavity & dispersion parameters  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     epso      water oxygen eps for implicit dispersion term */
/*     epsh      water hydrogen eps for implicit dispersion term */
/*     rmino     water oxygen Rmin for implicit dispersion term */
/*     rminh     water hydrogen Rmin for implicit dispersion term */
/*     awater    water number density at standard temp & pressure */
/*     slevy     enthalpy-to-free energy scale factor for dispersion */

/*     solvprs   limiting microscopic solvent pressure value */
/*     surften   limiting macroscopic surface tension value */
/*     spcut     starting radius for solvent pressure tapering */
/*     spoff     cutoff radius for solvent pressure tapering */
/*     stcut     starting radius for surface tension tapering */
/*     stoff     cutoff radius for surface tension tapering */
/*     rcav      atomic radius of each atom for cavitation energy */
/*     rdisp     atomic radius of each atom for dispersion energy */
/*     cdisp     maximum dispersion energy for each atom */




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
/*     nvt        number of distinct van der Waals types in the system */
/*     ivt        number of each distinct vdw type/class in the system */
/*     jvt        frequency of each vdw type or class in the system */




/*     zero out the Weeks-Chandler-Andersen dispersion energy */

    *edisp = 0.;

/*     set parameters for high accuracy numerical shells */

/*     tinit = 0.2d0 */
/*     rinit = 1.0d0 */
/*     rmult = 1.5d0 */
/*     rswitch = 7.0d0 */
/*     rmax = 12.0d0 */

/*     set parameters for medium accuracy numerical shells */

    tinit = 1.;
    rinit = 1.;
    rmult = 2.;
    rswitch = 5.;
    rmax = 9.;

/*     set parameters for low accuracy numerical shells */

/*     tinit = 1.0d0 */
/*     rinit = 1.0d0 */
/*     rmult = 2.0d0 */
/*     rswitch = 4.0d0 */
/*     rmax = 7.0d0 */

/*     set parameters for atomic radii and probe radii */

    offset = .55;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	npolar_1.rdisp[i__ - 1] = kvdws_1.rad[atmtyp_1.class__[i__ - 1] - 1] 
		+ .27;
	roff[i__ - 1] = npolar_1.rdisp[i__ - 1] + offset;
    }

/*     compute the dispersion energy for each atom in the system */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	epsi = kvdws_1.eps[atmtyp_1.class__[i__ - 1] - 1];
	rmini = kvdws_1.rad[atmtyp_1.class__[i__ - 1] - 1];
/* Computing 2nd power */
	d__1 = sqrt(.11) + sqrt(epsi);
	epsoi = epsi * .44 / (d__1 * d__1);
/* Computing 3rd power */
	d__1 = rmini;
/* Computing 2nd power */
	d__2 = rmini;
	rminoi = (d__1 * (d__1 * d__1) + 4.9347068906249989) * 2. / (d__2 * 
		d__2 + 2.8985062499999996);
/* Computing 2nd power */
	d__1 = sqrt(.0135) + sqrt(epsi);
	epshi = epsi * .053999999999999999 / (d__1 * d__1);
/* Computing 3rd power */
	d__1 = rmini;
/* Computing 2nd power */
	d__2 = rmini;
	rminhi = (d__1 * (d__1 * d__1) + 2.3393951718749992) * 2. / (d__2 * 
		d__2 + 1.7622562499999996);
/* Computing 7th power */
	d__1 = rminhi, d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
	her7 = epshi * (d__2 * (d__1 * d__1));
/* Computing 7th power */
	d__1 = rminoi, d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
	oer7 = epsoi * (d__2 * (d__1 * d__1));
/* Computing 14th power */
	d__1 = rminhi, d__1 *= d__1, d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
	her14 = epshi * (d__2 * (d__1 * d__1));
/* Computing 14th power */
	d__1 = rminoi, d__1 *= d__1, d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
	oer14 = epsoi * (d__2 * (d__1 * d__1));

/*     alter radii values for atoms attached to current atom */

	roff[i__ - 1] = npolar_1.rdisp[i__ - 1];
	i__2 = couple_1.n12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    k = i12_ref(j, i__);
	    roff[k - 1] = npolar_1.rdisp[k - 1];
	}
	i__2 = couple_1.n13[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    k = i13_ref(j, i__);
	    roff[k - 1] = npolar_1.rdisp[k - 1];
	}
	i__2 = couple_1.n14[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    k = i14_ref(j, i__);
	    roff[k - 1] = npolar_1.rdisp[k - 1];
	}
	i__2 = couple_1.n15[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    k = i15_ref(j, i__);
	    roff[k - 1] = npolar_1.rdisp[k - 1];
	}

/*     get the dispersion energy via a series of "onion" shells */

	t = tinit;
	ratio = rinit;
	e = 0.;
	done = FALSE_;
	while(! done) {
	    inner = roff[i__ - 1];
	    outer = inner + t;
	    roff[i__ - 1] = (inner + outer) * .5;
	    surfatom_(&i__, &area, roff);
/* Computing 2nd power */
	    d__1 = roff[i__ - 1];
	    fraction = area / (d__1 * d__1 * 12.566370614359172);
	    if (outer < rminoi) {
/* Computing 3rd power */
		d__1 = outer;
/* Computing 3rd power */
		d__2 = inner;
		shell = (d__1 * (d__1 * d__1) - d__2 * (d__2 * d__2)) / 3.;
		e -= epsoi * fraction * shell;
	    } else if (inner > rminoi) {
/* Computing 4th power */
		d__1 = inner, d__1 *= d__1;
/* Computing 4th power */
		d__2 = outer, d__2 *= d__2;
		shell = (1. / (d__1 * d__1) - 1. / (d__2 * d__2)) / 4.;
		e -= oer7 * 2. * fraction * shell;
/* Computing 11th power */
		d__1 = inner, d__2 = d__1, d__1 *= d__1, d__2 *= d__1, d__1 *=
			 d__1;
/* Computing 11th power */
		d__3 = outer, d__4 = d__3, d__3 *= d__3, d__4 *= d__3, d__3 *=
			 d__3;
		shell = (1. / (d__2 * (d__1 * d__1)) - 1. / (d__4 * (d__3 * 
			d__3))) / 11.;
		e += oer14 * fraction * shell;
	    } else {
/* Computing 3rd power */
		d__1 = rminoi;
/* Computing 3rd power */
		d__2 = inner;
		shell = (d__1 * (d__1 * d__1) - d__2 * (d__2 * d__2)) / 3.;
		e -= epsoi * fraction * shell;
/* Computing 4th power */
		d__1 = rminoi, d__1 *= d__1;
/* Computing 4th power */
		d__2 = outer, d__2 *= d__2;
		shell = (1. / (d__1 * d__1) - 1. / (d__2 * d__2)) / 4.;
		e -= oer7 * 2. * fraction * shell;
/* Computing 11th power */
		d__1 = rminoi, d__2 = d__1, d__1 *= d__1, d__2 *= d__1, d__1 
			*= d__1;
/* Computing 11th power */
		d__3 = outer, d__4 = d__3, d__3 *= d__3, d__4 *= d__3, d__3 *=
			 d__3;
		shell = (1. / (d__2 * (d__1 * d__1)) - 1. / (d__4 * (d__3 * 
			d__3))) / 11.;
		e += oer14 * fraction * shell;
	    }
	    if (outer < rminhi) {
/* Computing 3rd power */
		d__1 = outer;
/* Computing 3rd power */
		d__2 = inner;
		shell = (d__1 * (d__1 * d__1) - d__2 * (d__2 * d__2)) / 3.;
		e -= epshi * 2. * fraction * shell;
	    } else if (inner > rminhi) {
/* Computing 4th power */
		d__1 = inner, d__1 *= d__1;
/* Computing 4th power */
		d__2 = outer, d__2 *= d__2;
		shell = (1. / (d__1 * d__1) - 1. / (d__2 * d__2)) / 4.;
		e -= her7 * 4. * fraction * shell;
/* Computing 11th power */
		d__1 = inner, d__2 = d__1, d__1 *= d__1, d__2 *= d__1, d__1 *=
			 d__1;
/* Computing 11th power */
		d__3 = outer, d__4 = d__3, d__3 *= d__3, d__4 *= d__3, d__3 *=
			 d__3;
		shell = (1. / (d__2 * (d__1 * d__1)) - 1. / (d__4 * (d__3 * 
			d__3))) / 11.;
		e += her14 * 2. * fraction * shell;
	    } else {
/* Computing 3rd power */
		d__1 = rminhi;
/* Computing 3rd power */
		d__2 = inner;
		shell = (d__1 * (d__1 * d__1) - d__2 * (d__2 * d__2)) / 3.;
		e -= epshi * 2. * fraction * shell;
/* Computing 4th power */
		d__1 = rminhi, d__1 *= d__1;
/* Computing 4th power */
		d__2 = outer, d__2 *= d__2;
		shell = (1. / (d__1 * d__1) - 1. / (d__2 * d__2)) / 4.;
		e -= her7 * 4. * fraction * shell;
/* Computing 11th power */
		d__1 = rminhi, d__2 = d__1, d__1 *= d__1, d__2 *= d__1, d__1 
			*= d__1;
/* Computing 11th power */
		d__3 = outer, d__4 = d__3, d__3 *= d__3, d__4 *= d__3, d__3 *=
			 d__3;
		shell = (1. / (d__2 * (d__1 * d__1)) - 1. / (d__4 * (d__3 * 
			d__3))) / 11.;
		e += her14 * 2. * fraction * shell;
	    }
	    if (outer > rmax) {
		done = TRUE_;
	    }
	    if (fraction > .99 && outer > rminoi) {
		done = TRUE_;
	    }
	    if (done) {
/* Computing 4th power */
		d__1 = outer, d__1 *= d__1;
		e -= oer7 * 2. * fraction / (d__1 * d__1 * 4.);
/* Computing 11th power */
		d__1 = outer, d__2 = d__1, d__1 *= d__1, d__2 *= d__1, d__1 *=
			 d__1;
		e += oer14 * fraction / (d__2 * (d__1 * d__1) * 11.);
/* Computing 4th power */
		d__1 = outer, d__1 *= d__1;
		e -= her7 * 4. * fraction / (d__1 * d__1 * 4.);
/* Computing 11th power */
		d__1 = outer, d__2 = d__1, d__1 *= d__1, d__2 *= d__1, d__1 *=
			 d__1;
		e += her14 * 2. * fraction / (d__2 * (d__1 * d__1) * 11.);
	    }
	    roff[i__ - 1] += t * .5;
	    if (outer > rswitch) {
		ratio = rmult;
	    }
	    t = ratio * t;
	}

/*     increment the overall WCA dispersion energy component */

	e *= .42006863689679841;
	*edisp += e;

/*     reset the radii values for atoms attached to current atom */

	roff[i__ - 1] = npolar_1.rdisp[i__ - 1] + offset;
	i__2 = couple_1.n12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    k = i12_ref(j, i__);
	    roff[k - 1] = npolar_1.rdisp[k - 1] + offset;
	}
	i__2 = couple_1.n13[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    k = i13_ref(j, i__);
	    roff[k - 1] = npolar_1.rdisp[k - 1] + offset;
	}
	i__2 = couple_1.n14[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    k = i14_ref(j, i__);
	    roff[k - 1] = npolar_1.rdisp[k - 1] + offset;
	}
	i__2 = couple_1.n15[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    k = i15_ref(j, i__);
	    roff[k - 1] = npolar_1.rdisp[k - 1] + offset;
	}
    }
    return 0;
} /* ewcax_ */

#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref


