/* esolv1.f -- translated by f2c (version 20050501).
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
    doublereal einter;
} inter_;

#define inter_1 inter_

struct {
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

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
    doublereal vir[9]	/* was [3][3] */;
} virial_;

#define virial_1 virial_

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
    doublereal poleps, polsor, p2scale, p3scale, p4scale, p5scale, d1scale, 
	    d2scale, d3scale, d4scale, u1scale, u2scale, u3scale, u4scale;
    char poltyp[6];
} polpot_;

#define polpot_1 polpot_

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
    doublereal m2scale, m3scale, m4scale, m5scale;
} mplpot_;

#define mplpot_1 mplpot_

struct {
    integer np11[25000], ip11[2500000]	/* was [100][25000] */, np12[25000], 
	    ip12[1250000]	/* was [50][25000] */, np13[25000], ip13[
	    1250000]	/* was [50][25000] */, np14[25000], ip14[1250000]	
	    /* was [50][25000] */;
} polgrp_;

#define polgrp_1 polgrp_

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



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine esolv1  --  solvation energy and derivatives  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "esolv1" calculates the continuum solvation energy and */
/*     first derivatives with respect to Cartesian coordinates */
/*     for surface area, generalized Born, generalized Kirkwood */
/*     and Poisson-Boltzmann solvation models */


/* Subroutine */ int esolv1_(void)
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
    extern /* Subroutine */ int egk1_(void), epb1_(void), enp1_(doublereal *, 
	    doublereal *);
    static doublereal ecav, term;
    extern /* Subroutine */ int egb1a_(void), egb1b_(void), born1_(void);
    static doublereal edisp, probe;
    extern /* Subroutine */ int surface1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


#define des_ref(a_1,a_2) deriv_1.des[(a_2)*3 + a_1 - 4]



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




/*     zero out the continuum solvation energy and derivatives */

    energi_1.es = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	solute_1.drb[i__ - 1] = 0.;
	solute_1.drbp[i__ - 1] = 0.;
	des_ref(1, i__) = 0.;
	des_ref(2, i__) = 0.;
	des_ref(3, i__) = 0.;
    }

/*     set a value for the solvent molecule probe radius */

    probe = 1.4;

/*     solvation energy and derivs for surface area only models */

    if (s_cmp(solute_1.solvtyp, "ASP", (ftnlen)8, (ftnlen)3) == 0 || s_cmp(
	    solute_1.solvtyp, "SASA", (ftnlen)8, (ftnlen)4) == 0) {
	surface1_(&energi_1.es, aes, deriv_1.des, solute_1.rsolv, 
		solute_1.asolv, &probe);

/*     nonpolar energy and derivs for Onion method via exact area */

    } else if (s_cmp(solute_1.solvtyp, "ONION", (ftnlen)8, (ftnlen)5) == 0) {
	surface1_(&energi_1.es, aes, deriv_1.des, solute_1.rsolv, 
		solute_1.asolv, &probe);

/*     nonpolar energy and derivs as cavity plus dispersion */

    } else if (s_cmp(solute_1.solvtyp, "GK", (ftnlen)8, (ftnlen)2) == 0 || 
	    s_cmp(solute_1.solvtyp, "PB", (ftnlen)8, (ftnlen)2) == 0) {
	enp1_(&ecav, &edisp);
	energi_1.es = ecav + edisp;

/*     nonpolar energy and derivs via ACE area approximation */

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
		solute_1.drb[i__ - 1] -= e * 6. / rb;
	    }
	}
    }

/*     get polarization energy and derivatives for solvation methods */

    if (s_cmp(solute_1.solvtyp, "GK", (ftnlen)8, (ftnlen)2) == 0) {
	egk1_();
    } else if (s_cmp(solute_1.solvtyp, "PB", (ftnlen)8, (ftnlen)2) == 0) {
	epb1_();
    } else if (potent_1.use_born__) {
	if (warp_1.use_smooth__) {
	    egb1b_();
	} else {
	    egb1a_();
	}
	born1_();
    }
    return 0;
} /* esolv1_ */

#undef des_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine egb1a  --  generalized Born energy and derivs  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "egb1a" calculates the generalized Born electrostatic energy */
/*     and first derivatives of the GB/SA solvation models */

/*     note application of distance cutoff scaling directly to */
/*     the Born radii chain rule term "derb" is an approximation */


/* Subroutine */ int egb1a_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, k;
    static doublereal r__, r2, r3, r4, r5, r6, r7, de, fi;
    static integer ii, kk;
    static doublereal xi, yi, zi, xr, yr, zr, rb2, rm2, fgb, fik, rbi, rbk, 
	    fgb2, vxx, vyx, vyy, vzx, vzz, vzy, derb, drbi, drbk, dedx, dedy, 
	    dedz, fgrp;
    static logical usei;
    static doublereal taper, shift, trans, dtaper, dwater, dtrans;
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *);
    static logical proceed;
    static doublereal expterm;


#define des_ref(a_1,a_2) deriv_1.des[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  inter.i  --  sum of intermolecular energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     einter   total intermolecular potential energy */




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
	rbi = solute_1.rborn[i__ - 1];

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
		    r__ = sqrt(r2);
		    rbk = solute_1.rborn[k - 1];
		    fik = fi * charge_1.pchg[kk - 1];
		    rb2 = rbi * rbk;
		    expterm = exp(r2 * -.25 / rb2);
		    fgb2 = r2 + rb2 * expterm;
		    fgb = sqrt(fgb2);
		    e = fik / fgb;
		    de = -e * (r__ - r__ * .25 * expterm) / fgb2;
		    derb = -e * expterm * (r2 * .125 / rb2 + .5) / fgb2;

/*     use energy switching if near the cutoff distance */

/* Computing 2nd power */
		    d__1 = (shunt_1.off + shunt_1.cut) * .5;
		    rm2 = d__1 * d__1;
		    shift = fik / sqrt(rm2 + rb2 * exp(rm2 * -.25 / rb2));
		    e -= shift;
		    if (r2 > shunt_1.cut2) {
			r3 = r2 * r__;
			r4 = r2 * r2;
			r5 = r2 * r3;
			r6 = r3 * r3;
			r7 = r3 * r4;
			taper = shunt_1.c5 * r5 + shunt_1.c4 * r4 + 
				shunt_1.c3 * r3 + shunt_1.c2 * r2 + 
				shunt_1.c1 * r__ + shunt_1.c0;
			dtaper = shunt_1.c5 * 5. * r4 + shunt_1.c4 * 4. * r3 
				+ shunt_1.c3 * 3. * r2 + shunt_1.c2 * 2. * 
				r__ + shunt_1.c1;
			trans = fik * (shunt_1.f7 * r7 + shunt_1.f6 * r6 + 
				shunt_1.f5 * r5 + shunt_1.f4 * r4 + 
				shunt_1.f3 * r3 + shunt_1.f2 * r2 + 
				shunt_1.f1 * r__ + shunt_1.f0);
			dtrans = fik * (shunt_1.f7 * 7. * r6 + shunt_1.f6 * 
				6. * r5 + shunt_1.f5 * 5. * r4 + shunt_1.f4 * 
				4. * r3 + shunt_1.f3 * 3. * r2 + shunt_1.f2 * 
				2. * r__ + shunt_1.f1);
			derb *= taper;
			de = e * dtaper + de * taper + dtrans;
			e = e * taper + trans;
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
			de *= fgrp;
			derb *= fgrp;
		    }

/*     increment the overall energy and derivative expressions */

		    if (i__ == k) {
			e *= .5;
			energi_1.es += e;
			drbi = derb * rbk;
			solute_1.drb[i__ - 1] += drbi;
		    } else {
			energi_1.es += e;
			de /= r__;
			dedx = de * xr;
			dedy = de * yr;
			dedz = de * zr;
			des_ref(1, i__) = des_ref(1, i__) + dedx;
			des_ref(2, i__) = des_ref(2, i__) + dedy;
			des_ref(3, i__) = des_ref(3, i__) + dedz;
			des_ref(1, k) = des_ref(1, k) - dedx;
			des_ref(2, k) = des_ref(2, k) - dedy;
			des_ref(3, k) = des_ref(3, k) - dedz;
			drbi = derb * rbk;
			drbk = derb * rbi;
			solute_1.drb[i__ - 1] += drbi;
			solute_1.drb[k - 1] += drbk;

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

/*     increment the total intermolecular energy */

		    if (molcul_1.molcule[i__ - 1] != molcul_1.molcule[k - 1]) 
			    {
			inter_1.einter += e;
		    }
		}
	    }
	}
    }
    return 0;
} /* egb1a_ */

#undef vir_ref
#undef des_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine egb1b  --  GB energy and derivs for smoothing  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "egb1b" calculates the generalized Born energy and first */
/*     derivatives of the GB/SA solvation models for use with */
/*     potential smoothing methods */


/* Subroutine */ int egb1b_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, k;
    static doublereal r__, r2, de, fi;
    static integer ii, kk;
    static doublereal xi, yi, zi, xr, yr, zr, rb2, fgb, fik, rbi, rbk;
    extern doublereal erf_(doublereal *);
    static doublereal fgb2, vxx, vyx, vyy, vzx, vzz, vzy, derb, drbi, drbk, 
	    dedx, dedy, dedz, fgrp, term;
    static logical usei;
    static doublereal bterm, width, rterm, sterm, wterm, dwater;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;
    static doublereal erfterm, expterm;


#define des_ref(a_1,a_2) deriv_1.des[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  inter.i  --  sum of intermolecular energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     einter   total intermolecular potential energy */




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

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  virial.i  --  components of internal virial tensor  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     vir    total internal virial Cartesian tensor components */




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
	rbi = solute_1.rborn[i__ - 1];

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
		r__ = sqrt(r2);
		rbk = solute_1.rborn[k - 1];
		fik = fi * charge_1.pchg[kk - 1];
		rb2 = rbi * rbk;
		expterm = exp(r2 * -.25 / rb2);
		fgb2 = r2 + rb2 * expterm;
		fgb = sqrt(fgb2);
		e = fik / fgb;
		de = -e * (r__ - r__ * .25 * expterm) / fgb2;
		derb = -e * expterm * (r2 * .125 / rb2 + .5) / fgb2;

/*     use a smoothable GB analogous to the Coulomb solution */

		if (warp_1.deform > 0.) {
		    wterm = exp(rb2 * -.006 / warp_1.deform);
		    width = sterm / sqrt(warp_1.deform + rb2 * .15 * wterm);
		    d__1 = width * fgb;
		    erfterm = erf_(&d__1);
/* Computing 2nd power */
		    d__1 = width * fgb;
		    term = width * exp(-(d__1 * d__1)) / 1.772453850905516027;
		    rterm = term * (r__ * 2. - r__ * .5 * expterm) / fgb;
/* Computing 2nd power */
		    d__1 = width / sterm;
		    bterm = term * (expterm * (r2 * .25 / rb2 + 1.) / fgb - 
			    fgb * (d__1 * d__1) * wterm * (.15 - rb2 * 9e-4 / 
			    warp_1.deform));
		    derb = derb * erfterm + e * bterm;
		    de = de * erfterm + e * rterm;
		    e *= erfterm;
		}

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    e *= fgrp;
		    de *= fgrp;
		    derb *= fgrp;
		}

/*     increment the overall energy and derivative expressions */

		if (i__ == k) {
		    e *= .5;
		    energi_1.es += e;
		    drbi = derb * rbk;
		    solute_1.drb[i__ - 1] += drbi;
		} else {
		    energi_1.es += e;
		    de /= r__;
		    dedx = de * xr;
		    dedy = de * yr;
		    dedz = de * zr;
		    des_ref(1, i__) = des_ref(1, i__) + dedx;
		    des_ref(2, i__) = des_ref(2, i__) + dedy;
		    des_ref(3, i__) = des_ref(3, i__) + dedz;
		    des_ref(1, k) = des_ref(1, k) - dedx;
		    des_ref(2, k) = des_ref(2, k) - dedy;
		    des_ref(3, k) = des_ref(3, k) - dedz;
		    drbi = derb * rbk;
		    drbk = derb * rbi;
		    solute_1.drb[i__ - 1] += drbi;
		    solute_1.drb[k - 1] += drbk;

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

/*     increment the total intermolecular energy */

		if (molcul_1.molcule[i__ - 1] != molcul_1.molcule[k - 1]) {
		    inter_1.einter += e;
		}
	    }
	}
    }
    return 0;
} /* egb1b_ */

#undef vir_ref
#undef des_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine egk1  --  generalized Kirkwood energy & derivs  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "egk1" calculates the continuum solvation energy and derivatives */
/*     via the generalized Kirkwood plus nonpolar implicit solvation */


/* Subroutine */ int egk1_(void)
{
    extern /* Subroutine */ int egk1a_(void), born1_(void), ediff1_(void);



/*     compute the generalized Kirkwood energy and gradient */



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


    egk1a_();
    born1_();

/*     correct energy and derivatives for vacuum to polarized state */

    ediff1_();
    return 0;
} /* egk1_ */



/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine egk1a  --  find GK energy and derivatives  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "egk1a" calculates the electrostatic portion of the continuum */
/*     solvation energy and derivatives via the generalized Kirkwood */
/*     model */


/* Subroutine */ int egk1a_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal a[24]	/* was [6][4] */, b[15]	/* was [5][3] */, e;
    static integer i__, k;
    static doublereal r2, expcdexpc, fc, fd, ci, ei, ck, gf, gc[30];
    static integer ii, kk;
    static doublereal fq, xi, yi, zi, xr, yr, zr, gf2, gf3, gf5, gf7, rb2, 
	    gf9, xr2, yr2, zr2, gf11, fid[3], dxi, dyi, dzi, dxk, dyk, dzk, 
	    rbi, rbk, pxi, pyi, pzi, pxk, uxi, uyi, uzi, uxk, uyk, uzk, pyk, 
	    pzk, sxi, syi, szi, sxk, syk, szk, ewi, ewk, vxx, vyx, vyy, vzx, 
	    vzz, vzy, fkd[3], trq[75000]	/* was [3][25000] */, gux[30],
	     guy[30], guz[30], fidg[9]	/* was [3][3] */, fkdg[9]	/* 
	    was [3][3] */, drbi, drbk, dedx, dedy, dedz, dpbi, dpbk, fgrp, 
	    dpdx, dpdy, dpdz, expc, qxxi, qxyi, qxzi, qyyi, qyzi, qzzi, qxxk, 
	    qxyk, qxzk, qyyk, qyzk, qzzk, esym, ewii, ewki, trqi[75000]	/* 
	    was [3][25000] */, gqxx[30], gqxy[30], gqxz[30], gqyy[30], gqyz[
	    30], gqzz[30], expc1;
    static logical usei;
    static doublereal dgfdr, dexpc, expcr, esymi, duvdr, dewidr, dewkdr, 
	    dewidx, dwater, dexpcr, dewkdx, dewidy, dewkdy, dewidz, dewkdz, 
	    dsumdr, dsymdr, dpwidx, dpwkdx, dpwidy, dpwkdy, dpwidz, dpwkdz, 
	    dwipdr, dwkpdr;
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *), torque2_(doublereal *, doublereal *);
    static logical proceed;
    extern /* Subroutine */ int chkpole_(void);
    static doublereal desymdr, desymdx, desymdy, desymdz, expterm, dpsymdx, 
	    dpsymdy, dpsymdz;
    extern /* Subroutine */ int rotpole_(void);


#define a_ref(a_1,a_2) a[(a_2)*6 + a_1 - 0]
#define b_ref(a_1,a_2) b[(a_2)*5 + a_1 - 0]
#define des_ref(a_1,a_2) deriv_1.des[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define trq_ref(a_1,a_2) trq[(a_2)*3 + a_1 - 4]
#define fidg_ref(a_1,a_2) fidg[(a_2)*3 + a_1 - 4]
#define fkdg_ref(a_1,a_2) fkdg[(a_2)*3 + a_1 - 4]
#define trqi_ref(a_1,a_2) trqi[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define uinds_ref(a_1,a_2) polar_1.uinds[(a_2)*3 + a_1 - 4]
#define uinps_ref(a_1,a_2) polar_1.uinps[(a_2)*3 + a_1 - 4]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  inter.i  --  sum of intermolecular energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     einter   total intermolecular potential energy */




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

/*     zero out local accumulation arrays for derivatives */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (k = 1; k <= 3; ++k) {
	    trq_ref(k, i__) = 0.;
	    trqi_ref(k, i__) = 0.;
	}
    }

/*     calculate GK electrostatic solvation free energy */

    i__1 = mpole_1.npole;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = mpole_1.ipole[ii - 1];
	usei = usage_1.use[i__ - 1];
	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
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
	pxi = uinps_ref(1, ii);
	pyi = uinps_ref(2, ii);
	pzi = uinps_ref(3, ii);
	sxi = dxi + pxi;
	syi = dyi + pyi;
	szi = dzi + pzi;
	rbi = solute_1.rborn[i__ - 1];

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
		    pxk = uinps_ref(1, kk);
		    pyk = uinps_ref(2, kk);
		    pzk = uinps_ref(3, kk);
		    sxk = dxk + pxk;
		    syk = dyk + pyk;
		    szk = dzk + pzk;
		    rb2 = rbi * rbk;
		    expterm = exp(-r2 / (gk_1.gkc * rb2));
		    expc = expterm / gk_1.gkc;
		    expcr = r2 * expterm / (gk_1.gkc * gk_1.gkc * rb2 * rb2);
		    dexpc = -2. / (gk_1.gkc * rb2);
		    dexpcr = 2. / (gk_1.gkc * rb2 * rb2);
		    dgfdr = expterm * .5 * (r2 / (rb2 * gk_1.gkc) + 1.);
		    gf2 = 1. / (r2 + rb2 * expterm);
		    gf = sqrt(gf2);
		    gf3 = gf2 * gf;
		    gf5 = gf3 * gf2;
		    gf7 = gf5 * gf2;
		    gf9 = gf7 * gf2;
		    gf11 = gf9 * gf2;

/*     reaction potential auxiliary terms */

		    a_ref(0, 0) = gf;
		    a_ref(1, 0) = -gf3;
		    a_ref(2, 0) = gf5 * 3.;
		    a_ref(3, 0) = gf7 * -15.;
		    a_ref(4, 0) = gf9 * 105.;
		    a_ref(5, 0) = gf11 * -945.;

/*     Born radii derivatives of reaction potential auxiliary terms */

		    b_ref(0, 0) = dgfdr * a_ref(1, 0);
		    b_ref(1, 0) = dgfdr * a_ref(2, 0);
		    b_ref(2, 0) = dgfdr * a_ref(3, 0);
		    b_ref(3, 0) = dgfdr * a_ref(4, 0);
		    b_ref(4, 0) = dgfdr * a_ref(5, 0);

/*     reaction potential gradient auxiliary terms */

		    expc1 = 1. - expc;
		    a_ref(0, 1) = expc1 * a_ref(1, 0);
		    a_ref(1, 1) = expc1 * a_ref(2, 0);
		    a_ref(2, 1) = expc1 * a_ref(3, 0);
		    a_ref(3, 1) = expc1 * a_ref(4, 0);
		    a_ref(4, 1) = expc1 * a_ref(5, 0);

/*     Born radii derivs of reaction potential gradient auxiliary terms */

		    b_ref(0, 1) = b_ref(1, 0) - expcr * a_ref(1, 0) - expc * 
			    b_ref(1, 0);
		    b_ref(1, 1) = b_ref(2, 0) - expcr * a_ref(2, 0) - expc * 
			    b_ref(2, 0);
		    b_ref(2, 1) = b_ref(3, 0) - expcr * a_ref(3, 0) - expc * 
			    b_ref(3, 0);
		    b_ref(3, 1) = b_ref(4, 0) - expcr * a_ref(4, 0) - expc * 
			    b_ref(4, 0);

/*     2nd reaction potential gradient auxiliary terms */

		    expcdexpc = -expc * dexpc;
		    a_ref(0, 2) = expc1 * a_ref(1, 1) + expcdexpc * a_ref(1, 
			    0);
		    a_ref(1, 2) = expc1 * a_ref(2, 1) + expcdexpc * a_ref(2, 
			    0);
		    a_ref(2, 2) = expc1 * a_ref(3, 1) + expcdexpc * a_ref(3, 
			    0);
		    a_ref(3, 2) = expc1 * a_ref(4, 1) + expcdexpc * a_ref(4, 
			    0);

/*     Born radii derivatives of the 2nd reaction potential */
/*     gradient auxiliary terms */

		    b_ref(0, 2) = b_ref(1, 1) - (expcr * (a_ref(1, 1) + dexpc 
			    * a_ref(1, 0)) + expc * (b_ref(1, 1) + dexpcr * 
			    a_ref(1, 0) + dexpc * b_ref(1, 0)));
		    b_ref(1, 2) = b_ref(2, 1) - (expcr * (a_ref(2, 1) + dexpc 
			    * a_ref(2, 0)) + expc * (b_ref(2, 1) + dexpcr * 
			    a_ref(2, 0) + dexpc * b_ref(2, 0)));
		    b_ref(2, 2) = b_ref(3, 1) - (expcr * (a_ref(3, 1) + dexpc 
			    * a_ref(3, 0)) + expc * (b_ref(3, 1) + dexpcr * 
			    a_ref(3, 0) + dexpc * b_ref(3, 0)));

/*     3rd reaction potential gradient auxiliary terms */

		    expcdexpc *= 2.;
		    a_ref(0, 3) = expc1 * a_ref(1, 2) + expcdexpc * a_ref(1, 
			    1);
		    a_ref(1, 3) = expc1 * a_ref(2, 2) + expcdexpc * a_ref(2, 
			    1);
		    a_ref(2, 3) = expc1 * a_ref(3, 2) + expcdexpc * a_ref(3, 
			    1);
/* Computing 2nd power */
		    d__1 = dexpc;
		    expcdexpc = -expc * (d__1 * d__1);
		    a_ref(0, 3) = a_ref(0, 3) + expcdexpc * a_ref(1, 0);
		    a_ref(1, 3) = a_ref(1, 3) + expcdexpc * a_ref(2, 0);
		    a_ref(2, 3) = a_ref(2, 3) + expcdexpc * a_ref(3, 0);

/*     multiply the auxillary terms by their dieletric functions */

		    a_ref(0, 0) = fc * a_ref(0, 0);
		    a_ref(0, 1) = fc * a_ref(0, 1);
		    a_ref(0, 2) = fc * a_ref(0, 2);
		    a_ref(0, 3) = fc * a_ref(0, 3);
		    b_ref(0, 0) = fc * b_ref(0, 0);
		    b_ref(0, 1) = fc * b_ref(0, 1);
		    b_ref(0, 2) = fc * b_ref(0, 2);
		    a_ref(1, 0) = fd * a_ref(1, 0);
		    a_ref(1, 1) = fd * a_ref(1, 1);
		    a_ref(1, 2) = fd * a_ref(1, 2);
		    a_ref(1, 3) = fd * a_ref(1, 3);
		    b_ref(1, 0) = fd * b_ref(1, 0);
		    b_ref(1, 1) = fd * b_ref(1, 1);
		    b_ref(1, 2) = fd * b_ref(1, 2);
		    a_ref(2, 0) = fq * a_ref(2, 0);
		    a_ref(2, 1) = fq * a_ref(2, 1);
		    a_ref(2, 2) = fq * a_ref(2, 2);
		    a_ref(2, 3) = fq * a_ref(2, 3);
		    b_ref(2, 0) = fq * b_ref(2, 0);
		    b_ref(2, 1) = fq * b_ref(2, 1);
		    b_ref(2, 2) = fq * b_ref(2, 2);

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

/*     Born radii derivs of unweighted reaction potential tensor */

		    gc[20] = b_ref(0, 0);
		    gux[20] = xr * b_ref(1, 0);
		    guy[20] = yr * b_ref(1, 0);
		    guz[20] = zr * b_ref(1, 0);
		    gqxx[20] = xr2 * b_ref(2, 0);
		    gqyy[20] = yr2 * b_ref(2, 0);
		    gqzz[20] = zr2 * b_ref(2, 0);
		    gqxy[20] = xr * yr * b_ref(2, 0);
		    gqxz[20] = xr * zr * b_ref(2, 0);
		    gqyz[20] = yr * zr * b_ref(2, 0);

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

/*     Born derivs of the unweighted reaction potential gradient tensor */

		    gc[21] = xr * b_ref(0, 1);
		    gc[22] = yr * b_ref(0, 1);
		    gc[23] = zr * b_ref(0, 1);
		    gux[21] = b_ref(1, 0) + xr2 * b_ref(1, 1);
		    gux[22] = xr * yr * b_ref(1, 1);
		    gux[23] = xr * zr * b_ref(1, 1);
		    guy[21] = gux[22];
		    guy[22] = b_ref(1, 0) + yr2 * b_ref(1, 1);
		    guy[23] = yr * zr * b_ref(1, 1);
		    guz[21] = gux[23];
		    guz[22] = guy[23];
		    guz[23] = b_ref(1, 0) + zr2 * b_ref(1, 1);
		    gqxx[21] = xr * (b_ref(2, 0) * 2. + xr2 * b_ref(2, 1));
		    gqxx[22] = yr * xr2 * b_ref(2, 1);
		    gqxx[23] = zr * xr2 * b_ref(2, 1);
		    gqyy[21] = xr * yr2 * b_ref(2, 1);
		    gqyy[22] = yr * (b_ref(2, 0) * 2. + yr2 * b_ref(2, 1));
		    gqyy[23] = zr * yr2 * b_ref(2, 1);
		    gqzz[21] = xr * zr2 * b_ref(2, 1);
		    gqzz[22] = yr * zr2 * b_ref(2, 1);
		    gqzz[23] = zr * (b_ref(2, 0) * 2. + zr2 * b_ref(2, 1));
		    gqxy[21] = yr * (b_ref(2, 0) + xr2 * b_ref(2, 1));
		    gqxy[22] = xr * (b_ref(2, 0) + yr2 * b_ref(2, 1));
		    gqxy[23] = zr * xr * yr * b_ref(2, 1);
		    gqxz[21] = zr * (b_ref(2, 0) + xr2 * b_ref(2, 1));
		    gqxz[22] = gqxy[23];
		    gqxz[23] = xr * (b_ref(2, 0) + zr2 * b_ref(2, 1));
		    gqyz[21] = gqxy[23];
		    gqyz[22] = zr * (b_ref(2, 0) + yr2 * b_ref(2, 1));
		    gqyz[23] = yr * (b_ref(2, 0) + zr2 * b_ref(2, 1));

/*     unweighted 2nd reaction potential gradient tensor */

		    gc[4] = a_ref(0, 1) + xr2 * a_ref(0, 2);
		    gc[5] = xr * yr * a_ref(0, 2);
		    gc[6] = xr * zr * a_ref(0, 2);
		    gc[7] = a_ref(0, 1) + yr2 * a_ref(0, 2);
		    gc[8] = yr * zr * a_ref(0, 2);
		    gc[9] = a_ref(0, 1) + zr2 * a_ref(0, 2);
		    gux[4] = xr * (a_ref(1, 1) * 3. + xr2 * a_ref(1, 2));
		    gux[5] = yr * (a_ref(1, 1) + xr2 * a_ref(1, 2));
		    gux[6] = zr * (a_ref(1, 1) + xr2 * a_ref(1, 2));
		    gux[7] = xr * (a_ref(1, 1) + yr2 * a_ref(1, 2));
		    gux[8] = zr * xr * yr * a_ref(1, 2);
		    gux[9] = xr * (a_ref(1, 1) + zr2 * a_ref(1, 2));
		    guy[4] = yr * (a_ref(1, 1) + xr2 * a_ref(1, 2));
		    guy[5] = xr * (a_ref(1, 1) + yr2 * a_ref(1, 2));
		    guy[6] = gux[8];
		    guy[7] = yr * (a_ref(1, 1) * 3. + yr2 * a_ref(1, 2));
		    guy[8] = zr * (a_ref(1, 1) + yr2 * a_ref(1, 2));
		    guy[9] = yr * (a_ref(1, 1) + zr2 * a_ref(1, 2));
		    guz[4] = zr * (a_ref(1, 1) + xr2 * a_ref(1, 2));
		    guz[5] = gux[8];
		    guz[6] = xr * (a_ref(1, 1) + zr2 * a_ref(1, 2));
		    guz[7] = zr * (a_ref(1, 1) + yr2 * a_ref(1, 2));
		    guz[8] = yr * (a_ref(1, 1) + zr2 * a_ref(1, 2));
		    guz[9] = zr * (a_ref(1, 1) * 3. + zr2 * a_ref(1, 2));
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

/*     Born radii derivatives of the unweighted 2nd reaction */
/*     potential gradient tensor */

		    gc[24] = b_ref(0, 1) + xr2 * b_ref(0, 2);
		    gc[25] = xr * yr * b_ref(0, 2);
		    gc[26] = xr * zr * b_ref(0, 2);
		    gc[27] = b_ref(0, 1) + yr2 * b_ref(0, 2);
		    gc[28] = yr * zr * b_ref(0, 2);
		    gc[29] = b_ref(0, 1) + zr2 * b_ref(0, 2);
		    gux[24] = xr * (b_ref(1, 1) * 3. + xr2 * b_ref(1, 2));
		    gux[25] = yr * (b_ref(1, 1) + xr2 * b_ref(1, 2));
		    gux[26] = zr * (b_ref(1, 1) + xr2 * b_ref(1, 2));
		    gux[27] = xr * (b_ref(1, 1) + yr2 * b_ref(1, 2));
		    gux[28] = zr * xr * yr * b_ref(1, 2);
		    gux[29] = xr * (b_ref(1, 1) + zr2 * b_ref(1, 2));
		    guy[24] = yr * (b_ref(1, 1) + xr2 * b_ref(1, 2));
		    guy[25] = xr * (b_ref(1, 1) + yr2 * b_ref(1, 2));
		    guy[26] = gux[28];
		    guy[27] = yr * (b_ref(1, 1) * 3. + yr2 * b_ref(1, 2));
		    guy[28] = zr * (b_ref(1, 1) + yr2 * b_ref(1, 2));
		    guy[29] = yr * (b_ref(1, 1) + zr2 * b_ref(1, 2));
		    guz[24] = zr * (b_ref(1, 1) + xr2 * b_ref(1, 2));
		    guz[25] = gux[28];
		    guz[26] = xr * (b_ref(1, 1) + zr2 * b_ref(1, 2));
		    guz[27] = zr * (b_ref(1, 1) + yr2 * b_ref(1, 2));
		    guz[28] = yr * (b_ref(1, 1) + zr2 * b_ref(1, 2));
		    guz[29] = zr * (b_ref(1, 1) * 3. + zr2 * b_ref(1, 2));
		    gqxx[24] = b_ref(2, 0) * 2. + xr2 * (b_ref(2, 1) * 5. + 
			    xr2 * b_ref(2, 2));
		    gqxx[25] = yr * xr * (b_ref(2, 1) * 2. + xr2 * b_ref(2, 2)
			    );
		    gqxx[26] = zr * xr * (b_ref(2, 1) * 2. + xr2 * b_ref(2, 2)
			    );
		    gqxx[27] = xr2 * (b_ref(2, 1) + yr2 * b_ref(2, 2));
		    gqxx[28] = zr * yr * xr2 * b_ref(2, 2);
		    gqxx[29] = xr2 * (b_ref(2, 1) + zr2 * b_ref(2, 2));
		    gqyy[24] = yr2 * (b_ref(2, 1) + xr2 * b_ref(2, 2));
		    gqyy[25] = xr * yr * (b_ref(2, 1) * 2. + yr2 * b_ref(2, 2)
			    );
		    gqyy[26] = xr * zr * yr2 * b_ref(2, 2);
		    gqyy[27] = b_ref(2, 0) * 2. + yr2 * (b_ref(2, 1) * 5. + 
			    yr2 * b_ref(2, 2));
		    gqyy[28] = yr * zr * (b_ref(2, 1) * 2. + yr2 * b_ref(2, 2)
			    );
		    gqyy[29] = yr2 * (b_ref(2, 1) + zr2 * b_ref(2, 2));
		    gqzz[24] = zr2 * (b_ref(2, 1) + xr2 * b_ref(2, 2));
		    gqzz[25] = xr * yr * zr2 * b_ref(2, 2);
		    gqzz[26] = xr * zr * (b_ref(2, 1) * 2. + zr2 * b_ref(2, 2)
			    );
		    gqzz[27] = zr2 * (b_ref(2, 1) + yr2 * b_ref(2, 2));
		    gqzz[28] = yr * zr * (b_ref(2, 1) * 2. + zr2 * b_ref(2, 2)
			    );
		    gqzz[29] = b_ref(2, 0) * 2. + zr2 * (b_ref(2, 1) * 5. + 
			    zr2 * b_ref(2, 2));
		    gqxy[24] = xr * yr * (b_ref(2, 1) * 3. + xr2 * b_ref(2, 2)
			    );
		    gqxy[25] = b_ref(2, 0) + (xr2 + yr2) * b_ref(2, 1) + xr2 *
			     yr2 * b_ref(2, 2);
		    gqxy[26] = zr * yr * (b_ref(2, 1) + xr2 * b_ref(2, 2));
		    gqxy[27] = xr * yr * (b_ref(2, 1) * 3. + yr2 * b_ref(2, 2)
			    );
		    gqxy[28] = zr * xr * (b_ref(2, 1) + yr2 * b_ref(2, 2));
		    gqxy[29] = xr * yr * (b_ref(2, 1) + zr2 * b_ref(2, 2));
		    gqxz[24] = xr * zr * (b_ref(2, 1) * 3. + xr2 * b_ref(2, 2)
			    );
		    gqxz[25] = yr * zr * (b_ref(2, 1) + xr2 * b_ref(2, 2));
		    gqxz[26] = b_ref(2, 0) + (xr2 + zr2) * b_ref(2, 1) + xr2 *
			     zr2 * b_ref(2, 2);
		    gqxz[27] = xr * zr * (b_ref(2, 1) + yr2 * b_ref(2, 2));
		    gqxz[28] = xr * yr * (b_ref(2, 1) + zr2 * b_ref(2, 2));
		    gqxz[29] = xr * zr * (b_ref(2, 1) * 3. + zr2 * b_ref(2, 2)
			    );
		    gqyz[24] = zr * yr * (b_ref(2, 1) + xr2 * b_ref(2, 2));
		    gqyz[25] = xr * zr * (b_ref(2, 1) + yr2 * b_ref(2, 2));
		    gqyz[26] = xr * yr * (b_ref(2, 1) + zr2 * b_ref(2, 2));
		    gqyz[27] = yr * zr * (b_ref(2, 1) * 3. + yr2 * b_ref(2, 2)
			    );
		    gqyz[28] = b_ref(2, 0) + (yr2 + zr2) * b_ref(2, 1) + yr2 *
			     zr2 * b_ref(2, 2);
		    gqyz[29] = yr * zr * (b_ref(2, 1) * 3. + zr2 * b_ref(2, 2)
			    );

/*     unweighted 3rd reaction potential gradient tensor */

		    gc[10] = xr * (a_ref(0, 2) * 3. + xr2 * a_ref(0, 3));
		    gc[11] = yr * (a_ref(0, 2) + xr2 * a_ref(0, 3));
		    gc[12] = zr * (a_ref(0, 2) + xr2 * a_ref(0, 3));
		    gc[13] = xr * (a_ref(0, 2) + yr2 * a_ref(0, 3));
		    gc[14] = xr * yr * zr * a_ref(0, 3);
		    gc[15] = xr * (a_ref(0, 2) + zr2 * a_ref(0, 3));
		    gc[16] = yr * (a_ref(0, 2) * 3. + yr2 * a_ref(0, 3));
		    gc[17] = zr * (a_ref(0, 2) + yr2 * a_ref(0, 3));
		    gc[18] = yr * (a_ref(0, 2) + zr2 * a_ref(0, 3));
		    gc[19] = zr * (a_ref(0, 2) * 3. + zr2 * a_ref(0, 3));
		    gux[10] = a_ref(1, 1) * 3. + xr2 * (a_ref(1, 2) * 6. + 
			    xr2 * a_ref(1, 3));
		    gux[11] = xr * yr * (a_ref(1, 2) * 3. + xr2 * a_ref(1, 3))
			    ;
		    gux[12] = xr * zr * (a_ref(1, 2) * 3. + xr2 * a_ref(1, 3))
			    ;
		    gux[13] = a_ref(1, 1) + (xr2 + yr2) * a_ref(1, 2) + xr2 * 
			    yr2 * a_ref(1, 3);
		    gux[14] = yr * zr * (a_ref(1, 2) + xr2 * a_ref(1, 3));
		    gux[15] = a_ref(1, 1) + (xr2 + zr2) * a_ref(1, 2) + xr2 * 
			    zr2 * a_ref(1, 3);
		    gux[16] = xr * yr * (a_ref(1, 2) * 3. + yr2 * a_ref(1, 3))
			    ;
		    gux[17] = xr * zr * (a_ref(1, 2) + yr2 * a_ref(1, 3));
		    gux[18] = xr * yr * (a_ref(1, 2) + zr2 * a_ref(1, 3));
		    gux[19] = xr * zr * (a_ref(1, 2) * 3. + zr2 * a_ref(1, 3))
			    ;
		    guy[10] = gux[11];
		    guy[11] = gux[13];
		    guy[12] = gux[14];
		    guy[13] = gux[16];
		    guy[14] = gux[17];
		    guy[15] = gux[18];
		    guy[16] = a_ref(1, 1) * 3. + yr2 * (a_ref(1, 2) * 6. + 
			    yr2 * a_ref(1, 3));
		    guy[17] = yr * zr * (a_ref(1, 2) * 3. + yr2 * a_ref(1, 3))
			    ;
		    guy[18] = a_ref(1, 1) + (yr2 + zr2) * a_ref(1, 2) + yr2 * 
			    zr2 * a_ref(1, 3);
		    guy[19] = yr * zr * (a_ref(1, 2) * 3. + zr2 * a_ref(1, 3))
			    ;
		    guz[10] = gux[12];
		    guz[11] = gux[14];
		    guz[12] = gux[15];
		    guz[13] = gux[17];
		    guz[14] = gux[18];
		    guz[15] = gux[19];
		    guz[16] = guy[17];
		    guz[17] = guy[18];
		    guz[18] = guy[19];
		    guz[19] = a_ref(1, 1) * 3. + zr2 * (a_ref(1, 2) * 6. + 
			    zr2 * a_ref(1, 3));
		    gqxx[10] = xr * (a_ref(2, 1) * 12. + xr2 * (a_ref(2, 2) * 
			    9. + xr2 * a_ref(2, 3)));
		    gqxx[11] = yr * (a_ref(2, 1) * 2. + xr2 * (a_ref(2, 2) * 
			    5. + xr2 * a_ref(2, 3)));
		    gqxx[12] = zr * (a_ref(2, 1) * 2. + xr2 * (a_ref(2, 2) * 
			    5. + xr2 * a_ref(2, 3)));
		    gqxx[13] = xr * (a_ref(2, 1) * 2. + yr2 * 2. * a_ref(2, 2)
			     + xr2 * (a_ref(2, 2) + yr2 * a_ref(2, 3)));
		    gqxx[14] = xr * yr * zr * (a_ref(2, 2) * 2. + xr2 * a_ref(
			    2, 3));
		    gqxx[15] = xr * (a_ref(2, 1) * 2. + zr2 * 2. * a_ref(2, 2)
			     + xr2 * (a_ref(2, 2) + zr2 * a_ref(2, 3)));
		    gqxx[16] = yr * xr2 * (a_ref(2, 2) * 3. + yr2 * a_ref(2, 
			    3));
		    gqxx[17] = zr * xr2 * (a_ref(2, 2) + yr2 * a_ref(2, 3));
		    gqxx[18] = yr * xr2 * (a_ref(2, 2) + zr2 * a_ref(2, 3));
		    gqxx[19] = zr * xr2 * (a_ref(2, 2) * 3. + zr2 * a_ref(2, 
			    3));
		    gqxy[10] = yr * (a_ref(2, 1) * 3. + xr2 * (a_ref(2, 2) * 
			    6. + xr2 * a_ref(2, 3)));
		    gqxy[11] = xr * ((a_ref(2, 1) + yr2 * a_ref(2, 2)) * 3. + 
			    xr2 * (a_ref(2, 2) + yr2 * a_ref(2, 3)));
		    gqxy[12] = xr * yr * zr * (a_ref(2, 2) * 3. + xr2 * a_ref(
			    2, 3));
		    gqxy[13] = yr * ((a_ref(2, 1) + xr2 * a_ref(2, 2)) * 3. + 
			    yr2 * (a_ref(2, 2) + xr2 * a_ref(2, 3)));
		    gqxy[14] = zr * (a_ref(2, 1) + (yr2 + xr2) * a_ref(2, 2) 
			    + yr2 * xr2 * a_ref(2, 3));
		    gqxy[15] = yr * (a_ref(2, 1) + (xr2 + zr2) * a_ref(2, 2) 
			    + xr2 * zr2 * a_ref(2, 3));
		    gqxy[16] = xr * ((a_ref(2, 1) + yr2 * a_ref(2, 2)) * 3. + 
			    yr2 * (a_ref(2, 2) * 3. + yr2 * a_ref(2, 3)));
		    gqxy[17] = xr * yr * zr * (a_ref(2, 2) * 3. + yr2 * a_ref(
			    2, 3));
		    gqxy[18] = xr * (a_ref(2, 1) + (yr2 + zr2) * a_ref(2, 2) 
			    + yr2 * zr2 * a_ref(2, 3));
		    gqxy[19] = xr * yr * zr * (a_ref(2, 2) * 3. + zr2 * a_ref(
			    2, 3));
		    gqxz[10] = zr * (a_ref(2, 1) * 3. + xr2 * (a_ref(2, 2) * 
			    6. + xr2 * a_ref(2, 3)));
		    gqxz[11] = xr * yr * zr * (a_ref(2, 2) * 3. + xr2 * a_ref(
			    2, 3));
		    gqxz[12] = xr * ((a_ref(2, 1) + zr2 * a_ref(2, 2)) * 3. + 
			    xr2 * (a_ref(2, 2) + zr2 * a_ref(2, 3)));
		    gqxz[13] = zr * (a_ref(2, 1) + (xr2 + yr2) * a_ref(2, 2) 
			    + xr2 * yr2 * a_ref(2, 3));
		    gqxz[14] = yr * (a_ref(2, 1) + (xr2 + zr2) * a_ref(2, 2) 
			    + zr2 * xr2 * a_ref(2, 3));
		    gqxz[15] = zr * ((a_ref(2, 1) + xr2 * a_ref(2, 2)) * 3. + 
			    zr2 * (a_ref(2, 2) + xr2 * a_ref(2, 3)));
		    gqxz[16] = xr * yr * zr * (a_ref(2, 2) * 3. + yr2 * a_ref(
			    2, 3));
		    gqxz[17] = xr * (a_ref(2, 1) + (zr2 + yr2) * a_ref(2, 2) 
			    + zr2 * yr2 * a_ref(2, 3));
		    gqxz[18] = xr * yr * zr * (a_ref(2, 2) * 3. + zr2 * a_ref(
			    2, 3));
		    gqxz[19] = xr * (a_ref(2, 1) * 3. + zr2 * (a_ref(2, 2) * 
			    6. + zr2 * a_ref(2, 3)));
		    gqyy[10] = xr * yr2 * (a_ref(2, 2) * 3. + xr2 * a_ref(2, 
			    3));
		    gqyy[11] = yr * (a_ref(2, 1) * 2. + xr2 * 2. * a_ref(2, 2)
			     + yr2 * (a_ref(2, 2) + xr2 * a_ref(2, 3)));
		    gqyy[12] = zr * yr2 * (a_ref(2, 2) + xr2 * a_ref(2, 3));
		    gqyy[13] = xr * (a_ref(2, 1) * 2. + yr2 * (a_ref(2, 2) * 
			    5. + yr2 * a_ref(2, 3)));
		    gqyy[14] = xr * yr * zr * (a_ref(2, 2) * 2. + yr2 * a_ref(
			    2, 3));
		    gqyy[15] = xr * yr2 * (a_ref(2, 2) + zr2 * a_ref(2, 3));
		    gqyy[16] = yr * (a_ref(2, 1) * 12. + yr2 * (a_ref(2, 2) * 
			    9. + yr2 * a_ref(2, 3)));
		    gqyy[17] = zr * (a_ref(2, 1) * 2. + yr2 * (a_ref(2, 2) * 
			    5. + yr2 * a_ref(2, 3)));
		    gqyy[18] = yr * (a_ref(2, 1) * 2. + zr2 * 2. * a_ref(2, 2)
			     + yr2 * (a_ref(2, 2) + zr2 * a_ref(2, 3)));
		    gqyy[19] = zr * yr2 * (a_ref(2, 2) * 3. + zr2 * a_ref(2, 
			    3));
		    gqyz[10] = xr * yr * zr * (a_ref(2, 2) * 3. + xr2 * a_ref(
			    2, 3));
		    gqyz[11] = zr * (a_ref(2, 1) + (xr2 + yr2) * a_ref(2, 2) 
			    + xr2 * yr2 * a_ref(2, 3));
		    gqyz[12] = yr * (a_ref(2, 1) + (xr2 + zr2) * a_ref(2, 2) 
			    + xr2 * zr2 * a_ref(2, 3));
		    gqyz[13] = xr * yr * zr * (a_ref(2, 2) * 3. + yr2 * a_ref(
			    2, 3));
		    gqyz[14] = xr * (a_ref(2, 1) + (yr2 + zr2) * a_ref(2, 2) 
			    + yr2 * zr2 * a_ref(2, 3));
		    gqyz[15] = xr * yr * zr * (a_ref(2, 2) * 3. + zr2 * a_ref(
			    2, 3));
		    gqyz[16] = zr * (a_ref(2, 1) * 3. + yr2 * (a_ref(2, 2) * 
			    6. + yr2 * a_ref(2, 3)));
		    gqyz[17] = yr * ((a_ref(2, 1) + zr2 * a_ref(2, 2)) * 3. + 
			    yr2 * (a_ref(2, 2) + zr2 * a_ref(2, 3)));
		    gqyz[18] = zr * ((a_ref(2, 1) + yr2 * a_ref(2, 2)) * 3. + 
			    zr2 * (a_ref(2, 2) + yr2 * a_ref(2, 3)));
		    gqyz[19] = yr * (a_ref(2, 1) * 3. + zr2 * (a_ref(2, 2) * 
			    6. + zr2 * a_ref(2, 3)));
		    gqzz[10] = xr * zr2 * (a_ref(2, 2) * 3. + xr2 * a_ref(2, 
			    3));
		    gqzz[11] = yr * (zr2 * a_ref(2, 2) + xr2 * (zr2 * a_ref(2,
			     3)));
		    gqzz[12] = zr * (a_ref(2, 1) * 2. + xr2 * 2. * a_ref(2, 2)
			     + zr2 * (a_ref(2, 2) + xr2 * a_ref(2, 3)));
		    gqzz[13] = xr * zr2 * (a_ref(2, 2) + yr2 * a_ref(2, 3));
		    gqzz[14] = xr * yr * zr * (a_ref(2, 2) * 2. + zr2 * a_ref(
			    2, 3));
		    gqzz[15] = xr * (a_ref(2, 1) * 2. + zr2 * (a_ref(2, 2) * 
			    5. + zr2 * a_ref(2, 3)));
		    gqzz[16] = yr * zr2 * (a_ref(2, 2) * 3. + yr2 * a_ref(2, 
			    3));
		    gqzz[17] = zr * (a_ref(2, 1) * 2. + yr2 * 2. * a_ref(2, 2)
			     + zr2 * (a_ref(2, 2) + yr2 * a_ref(2, 3)));
		    gqzz[18] = yr * (a_ref(2, 1) * 2. + zr2 * (a_ref(2, 2) * 
			    5. + zr2 * a_ref(2, 3)));
		    gqzz[19] = zr * (a_ref(2, 1) * 12. + zr2 * (a_ref(2, 2) * 
			    9. + zr2 * a_ref(2, 3)));

/*     electrostatic solvation energy of the permanent multipoles */
/*     in their own GK reaction potential */

		    esym = ci * ck * gc[0] - (uxi * (uxk * gux[1] + uyk * guy[
			    1] + uzk * guz[1]) + uyi * (uxk * gux[2] + uyk * 
			    guy[2] + uzk * guz[2]) + uzi * (uxk * gux[3] + 
			    uyk * guy[3] + uzk * guz[3]));
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

		    desymdx = ci * ck * gc[1] - (uxi * (uxk * gux[4] + uyk * 
			    guy[4] + uzk * guz[4]) + uyi * (uxk * gux[5] + 
			    uyk * guy[5] + uzk * guz[5]) + uzi * (uxk * gux[6]
			     + uyk * guy[6] + uzk * guz[6]));
		    dewidx = ci * (uxk * gc[4] + uyk * gc[5] + uzk * gc[6]) - 
			    ck * (uxi * gux[1] + uyi * guy[1] + uzi * guz[1]) 
			    + ci * (qxxk * gc[10] + qyyk * gc[13] + qzzk * gc[
			    15] + (qxyk * gc[11] + qxzk * gc[12] + qyzk * gc[
			    14]) * 2.) + ck * (qxxi * gqxx[1] + qyyi * gqyy[1]
			     + qzzi * gqzz[1] + (qxyi * gqxy[1] + qxzi * gqxz[
			    1] + qyzi * gqyz[1]) * 2.) - uxi * (qxxk * gux[10]
			     + qyyk * gux[13] + qzzk * gux[15] + (qxyk * gux[
			    11] + qxzk * gux[12] + qyzk * gux[14]) * 2.) - 
			    uyi * (qxxk * guy[10] + qyyk * guy[13] + qzzk * 
			    guy[15] + (qxyk * guy[11] + qxzk * guy[12] + qyzk 
			    * guy[14]) * 2.) - uzi * (qxxk * guz[10] + qyyk * 
			    guz[13] + qzzk * guz[15] + (qxyk * guz[11] + qxzk 
			    * guz[12] + qyzk * guz[14]) * 2.) + uxk * (qxxi * 
			    gqxx[4] + qyyi * gqyy[4] + qzzi * gqzz[4] + (qxyi 
			    * gqxy[4] + qxzi * gqxz[4] + qyzi * gqyz[4]) * 2.)
			     + uyk * (qxxi * gqxx[5] + qyyi * gqyy[5] + qzzi *
			     gqzz[5] + (qxyi * gqxy[5] + qxzi * gqxz[5] + 
			    qyzi * gqyz[5]) * 2.) + uzk * (qxxi * gqxx[6] + 
			    qyyi * gqyy[6] + qzzi * gqzz[6] + (qxyi * gqxy[6] 
			    + qxzi * gqxz[6] + qyzi * gqyz[6]) * 2.) + qxxi * 
			    (qxxk * gqxx[10] + qyyk * gqxx[13] + qzzk * gqxx[
			    15] + (qxyk * gqxx[11] + qxzk * gqxx[12] + qyzk * 
			    gqxx[14]) * 2.) + qyyi * (qxxk * gqyy[10] + qyyk *
			     gqyy[13] + qzzk * gqyy[15] + (qxyk * gqyy[11] + 
			    qxzk * gqyy[12] + qyzk * gqyy[14]) * 2.) + qzzi * 
			    (qxxk * gqzz[10] + qyyk * gqzz[13] + qzzk * gqzz[
			    15] + (qxyk * gqzz[11] + qxzk * gqzz[12] + qyzk * 
			    gqzz[14]) * 2.) + (qxyi * (qxxk * gqxy[10] + qyyk 
			    * gqxy[13] + qzzk * gqxy[15] + (qxyk * gqxy[11] + 
			    qxzk * gqxy[12] + qyzk * gqxy[14]) * 2.) + qxzi * 
			    (qxxk * gqxz[10] + qyyk * gqxz[13] + qzzk * gqxz[
			    15] + (qxyk * gqxz[11] + qxzk * gqxz[12] + qyzk * 
			    gqxz[14]) * 2.) + qyzi * (qxxk * gqyz[10] + qyyk *
			     gqyz[13] + qzzk * gqyz[15] + (qxyk * gqyz[11] + 
			    qxzk * gqyz[12] + qyzk * gqyz[14]) * 2.)) * 2.;
		    dewkdx = ci * (uxk * gux[1] + uyk * guy[1] + uzk * guz[1])
			     - ck * (uxi * gc[4] + uyi * gc[5] + uzi * gc[6]) 
			    + ci * (qxxk * gqxx[1] + qyyk * gqyy[1] + qzzk * 
			    gqzz[1] + (qxyk * gqxy[1] + qxzk * gqxz[1] + qyzk 
			    * gqyz[1]) * 2.) + ck * (qxxi * gc[10] + qyyi * 
			    gc[13] + qzzi * gc[15] + (qxyi * gc[11] + qxzi * 
			    gc[12] + qyzi * gc[14]) * 2.) - uxi * (qxxk * 
			    gqxx[4] + qyyk * gqyy[4] + qzzk * gqzz[4] + (qxyk 
			    * gqxy[4] + qxzk * gqxz[4] + qyzk * gqyz[4]) * 2.)
			     - uyi * (qxxk * gqxx[5] + qyyk * gqyy[5] + qzzk *
			     gqzz[5] + (qxyk * gqxy[5] + qxzk * gqxz[5] + 
			    qyzk * gqyz[5]) * 2.) - uzi * (qxxk * gqxx[6] + 
			    qyyk * gqyy[6] + qzzk * gqzz[6] + (qxyk * gqxy[6] 
			    + qxzk * gqxz[6] + qyzk * gqyz[6]) * 2.) + uxk * (
			    qxxi * gux[10] + qyyi * gux[13] + qzzi * gux[15] 
			    + (qxyi * gux[11] + qxzi * gux[12] + qyzi * gux[
			    14]) * 2.) + uyk * (qxxi * guy[10] + qyyi * guy[
			    13] + qzzi * guy[15] + (qxyi * guy[11] + qxzi * 
			    guy[12] + qyzi * guy[14]) * 2.) + uzk * (qxxi * 
			    guz[10] + qyyi * guz[13] + qzzi * guz[15] + (qxyi 
			    * guz[11] + qxzi * guz[12] + qyzi * guz[14]) * 2.)
			     + qxxi * (qxxk * gqxx[10] + qyyk * gqyy[10] + 
			    qzzk * gqzz[10] + (qxyk * gqxy[10] + qxzk * gqxz[
			    10] + qyzk * gqyz[10]) * 2.) + qyyi * (qxxk * 
			    gqxx[13] + qyyk * gqyy[13] + qzzk * gqzz[13] + (
			    qxyk * gqxy[13] + qxzk * gqxz[13] + qyzk * gqyz[
			    13]) * 2.) + qzzi * (qxxk * gqxx[15] + qyyk * 
			    gqyy[15] + qzzk * gqzz[15] + (qxyk * gqxy[15] + 
			    qxzk * gqxz[15] + qyzk * gqyz[15]) * 2.) + (qxyi *
			     (qxxk * gqxx[11] + qyyk * gqyy[11] + qzzk * gqzz[
			    11] + (qxyk * gqxy[11] + qxzk * gqxz[11] + qyzk * 
			    gqyz[11]) * 2.) + qxzi * (qxxk * gqxx[12] + qyyk *
			     gqyy[12] + qzzk * gqzz[12] + (qxyk * gqxy[12] + 
			    qxzk * gqxz[12] + qyzk * gqyz[12]) * 2.) + qyzi * 
			    (qxxk * gqxx[14] + qyyk * gqyy[14] + qzzk * gqzz[
			    14] + (qxyk * gqxy[14] + qxzk * gqxz[14] + qyzk * 
			    gqyz[14]) * 2.)) * 2.;
		    dedx = desymdx + (dewidx + dewkdx) * .5;

		    desymdy = ci * ck * gc[2] - (uxi * (uxk * gux[5] + uyk * 
			    guy[5] + uzk * guz[5]) + uyi * (uxk * gux[7] + 
			    uyk * guy[7] + uzk * guz[7]) + uzi * (uxk * gux[8]
			     + uyk * guy[8] + uzk * guz[8]));
		    dewidy = ci * (uxk * gc[5] + uyk * gc[7] + uzk * gc[8]) - 
			    ck * (uxi * gux[2] + uyi * guy[2] + uzi * guz[2]) 
			    + ci * (qxxk * gc[11] + qyyk * gc[16] + qzzk * gc[
			    18] + (qxyk * gc[13] + qxzk * gc[14] + qyzk * gc[
			    17]) * 2.) + ck * (qxxi * gqxx[2] + qyyi * gqyy[2]
			     + qzzi * gqzz[2] + (qxyi * gqxy[2] + qxzi * gqxz[
			    2] + qyzi * gqyz[2]) * 2.) - uxi * (qxxk * gux[11]
			     + qyyk * gux[16] + qzzk * gux[18] + (qxyk * gux[
			    13] + qxzk * gux[14] + qyzk * gux[17]) * 2.) - 
			    uyi * (qxxk * guy[11] + qyyk * guy[16] + qzzk * 
			    guy[18] + (qxyk * guy[13] + qxzk * guy[14] + qyzk 
			    * guy[17]) * 2.) - uzi * (qxxk * guz[11] + qyyk * 
			    guz[16] + qzzk * guz[18] + (qxyk * guz[13] + qxzk 
			    * guz[14] + qyzk * guz[17]) * 2.) + uxk * (qxxi * 
			    gqxx[5] + qyyi * gqyy[5] + qzzi * gqzz[5] + (qxyi 
			    * gqxy[5] + qxzi * gqxz[5] + qyzi * gqyz[5]) * 2.)
			     + uyk * (qxxi * gqxx[7] + qyyi * gqyy[7] + qzzi *
			     gqzz[7] + (qxyi * gqxy[7] + qxzi * gqxz[7] + 
			    qyzi * gqyz[7]) * 2.) + uzk * (qxxi * gqxx[8] + 
			    qyyi * gqyy[8] + qzzi * gqzz[8] + (qxyi * gqxy[8] 
			    + qxzi * gqxz[8] + qyzi * gqyz[8]) * 2.) + qxxi * 
			    (qxxk * gqxx[11] + qyyk * gqxx[16] + qzzk * gqxx[
			    18] + (qxyk * gqxx[13] + qxzk * gqxx[14] + qyzk * 
			    gqxx[17]) * 2.) + qyyi * (qxxk * gqyy[11] + qyyk *
			     gqyy[16] + qzzk * gqyy[18] + (qxyk * gqyy[13] + 
			    qxzk * gqyy[14] + qyzk * gqyy[17]) * 2.) + qzzi * 
			    (qxxk * gqzz[11] + qyyk * gqzz[16] + qzzk * gqzz[
			    18] + (qxyk * gqzz[13] + qxzk * gqzz[14] + qyzk * 
			    gqzz[17]) * 2.) + (qxyi * (qxxk * gqxy[11] + qyyk 
			    * gqxy[16] + qzzk * gqxy[18] + (qxyk * gqxy[13] + 
			    qxzk * gqxy[14] + qyzk * gqxy[17]) * 2.) + qxzi * 
			    (qxxk * gqxz[11] + qyyk * gqxz[16] + qzzk * gqxz[
			    18] + (qxyk * gqxz[13] + qxzk * gqxz[14] + qyzk * 
			    gqxz[17]) * 2.) + qyzi * (qxxk * gqyz[11] + qyyk *
			     gqyz[16] + qzzk * gqyz[18] + (qxyk * gqyz[13] + 
			    qxzk * gqyz[14] + qyzk * gqyz[17]) * 2.)) * 2.;
		    dewkdy = ci * (uxk * gux[2] + uyk * guy[2] + uzk * guz[2])
			     - ck * (uxi * gc[5] + uyi * gc[7] + uzi * gc[8]) 
			    + ci * (qxxk * gqxx[2] + qyyk * gqyy[2] + qzzk * 
			    gqzz[2] + (qxyk * gqxy[2] + qxzk * gqxz[2] + qyzk 
			    * gqyz[2]) * 2.) + ck * (qxxi * gc[11] + qyyi * 
			    gc[16] + qzzi * gc[18] + (qxyi * gc[13] + qxzi * 
			    gc[14] + qyzi * gc[17]) * 2.) - uxi * (qxxk * 
			    gqxx[5] + qyyk * gqyy[5] + qzzk * gqzz[5] + (qxyk 
			    * gqxy[5] + qxzk * gqxz[5] + qyzk * gqyz[5]) * 2.)
			     - uyi * (qxxk * gqxx[7] + qyyk * gqyy[7] + qzzk *
			     gqzz[7] + (qxyk * gqxy[7] + qxzk * gqxz[7] + 
			    qyzk * gqyz[7]) * 2.) - uzi * (qxxk * gqxx[8] + 
			    qyyk * gqyy[8] + qzzk * gqzz[8] + (qxyk * gqxy[8] 
			    + qxzk * gqxz[8] + qyzk * gqyz[8]) * 2.) + uxk * (
			    qxxi * gux[11] + qyyi * gux[16] + qzzi * gux[18] 
			    + (qxyi * gux[13] + qxzi * gux[14] + qyzi * gux[
			    17]) * 2.) + uyk * (qxxi * guy[11] + qyyi * guy[
			    16] + qzzi * guy[18] + (qxyi * guy[13] + qxzi * 
			    guy[14] + qyzi * guy[17]) * 2.) + uzk * (qxxi * 
			    guz[11] + qyyi * guz[16] + qzzi * guz[18] + (qxyi 
			    * guz[13] + qxzi * guz[14] + qyzi * guz[17]) * 2.)
			     + qxxi * (qxxk * gqxx[11] + qyyk * gqyy[11] + 
			    qzzk * gqzz[11] + (qxyk * gqxy[11] + qxzk * gqxz[
			    11] + qyzk * gqyz[11]) * 2.) + qyyi * (qxxk * 
			    gqxx[16] + qyyk * gqyy[16] + qzzk * gqzz[16] + (
			    qxyk * gqxy[16] + qxzk * gqxz[16] + qyzk * gqyz[
			    16]) * 2.) + qzzi * (qxxk * gqxx[18] + qyyk * 
			    gqyy[18] + qzzk * gqzz[18] + (qxyk * gqxy[18] + 
			    qxzk * gqxz[18] + qyzk * gqyz[18]) * 2.) + (qxyi *
			     (qxxk * gqxx[13] + qyyk * gqyy[13] + qzzk * gqzz[
			    13] + (qxyk * gqxy[13] + qxzk * gqxz[13] + qyzk * 
			    gqyz[13]) * 2.) + qxzi * (qxxk * gqxx[14] + qyyk *
			     gqyy[14] + qzzk * gqzz[14] + (qxyk * gqxy[14] + 
			    qxzk * gqxz[14] + qyzk * gqyz[14]) * 2.) + qyzi * 
			    (qxxk * gqxx[17] + qyyk * gqyy[17] + qzzk * gqzz[
			    17] + (qxyk * gqxy[17] + qxzk * gqxz[17] + qyzk * 
			    gqyz[17]) * 2.)) * 2.;
		    dedy = desymdy + (dewidy + dewkdy) * .5;

		    desymdz = ci * ck * gc[3] - (uxi * (uxk * gux[6] + uyk * 
			    guy[6] + uzk * guz[6]) + uyi * (uxk * gux[8] + 
			    uyk * guy[8] + uzk * guz[8]) + uzi * (uxk * gux[9]
			     + uyk * guy[9] + uzk * guz[9]));
		    dewidz = ci * (uxk * gc[6] + uyk * gc[8] + uzk * gc[9]) - 
			    ck * (uxi * gux[3] + uyi * guy[3] + uzi * guz[3]) 
			    + ci * (qxxk * gc[12] + qyyk * gc[17] + qzzk * gc[
			    19] + (qxyk * gc[14] + qxzk * gc[15] + qyzk * gc[
			    18]) * 2.) + ck * (qxxi * gqxx[3] + qyyi * gqyy[3]
			     + qzzi * gqzz[3] + (qxyi * gqxy[3] + qxzi * gqxz[
			    3] + qyzi * gqyz[3]) * 2.) - uxi * (qxxk * gux[12]
			     + qyyk * gux[17] + qzzk * gux[19] + (qxyk * gux[
			    14] + qxzk * gux[15] + qyzk * gux[18]) * 2.) - 
			    uyi * (qxxk * guy[12] + qyyk * guy[17] + qzzk * 
			    guy[19] + (qxyk * guy[14] + qxzk * guy[15] + qyzk 
			    * guy[18]) * 2.) - uzi * (qxxk * guz[12] + qyyk * 
			    guz[17] + qzzk * guz[19] + (qxyk * guz[14] + qxzk 
			    * guz[15] + qyzk * guz[18]) * 2.) + uxk * (qxxi * 
			    gqxx[6] + qyyi * gqyy[6] + qzzi * gqzz[6] + (qxyi 
			    * gqxy[6] + qxzi * gqxz[6] + qyzi * gqyz[6]) * 2.)
			     + uyk * (qxxi * gqxx[8] + qyyi * gqyy[8] + qzzi *
			     gqzz[8] + (qxyi * gqxy[8] + qxzi * gqxz[8] + 
			    qyzi * gqyz[8]) * 2.) + uzk * (qxxi * gqxx[9] + 
			    qyyi * gqyy[9] + qzzi * gqzz[9] + (qxyi * gqxy[9] 
			    + qxzi * gqxz[9] + qyzi * gqyz[9]) * 2.) + qxxi * 
			    (qxxk * gqxx[12] + qyyk * gqxx[17] + qzzk * gqxx[
			    19] + (qxyk * gqxx[14] + qxzk * gqxx[15] + qyzk * 
			    gqxx[18]) * 2.) + qyyi * (qxxk * gqyy[12] + qyyk *
			     gqyy[17] + qzzk * gqyy[19] + (qxyk * gqyy[14] + 
			    qxzk * gqyy[15] + qyzk * gqyy[18]) * 2.) + qzzi * 
			    (qxxk * gqzz[12] + qyyk * gqzz[17] + qzzk * gqzz[
			    19] + (qxyk * gqzz[14] + qxzk * gqzz[15] + qyzk * 
			    gqzz[18]) * 2.) + (qxyi * (qxxk * gqxy[12] + qyyk 
			    * gqxy[17] + qzzk * gqxy[19] + (qxyk * gqxy[14] + 
			    qxzk * gqxy[15] + qyzk * gqxy[18]) * 2.) + qxzi * 
			    (qxxk * gqxz[12] + qyyk * gqxz[17] + qzzk * gqxz[
			    19] + (qxyk * gqxz[14] + qxzk * gqxz[15] + qyzk * 
			    gqxz[18]) * 2.) + qyzi * (qxxk * gqyz[12] + qyyk *
			     gqyz[17] + qzzk * gqyz[19] + (qxyk * gqyz[14] + 
			    qxzk * gqyz[15] + qyzk * gqyz[18]) * 2.)) * 2.;
		    dewkdz = ci * (uxk * gux[3] + uyk * guy[3] + uzk * guz[3])
			     - ck * (uxi * gc[6] + uyi * gc[8] + uzi * gc[9]) 
			    + ci * (qxxk * gqxx[3] + qyyk * gqyy[3] + qzzk * 
			    gqzz[3] + (qxyk * gqxy[3] + qxzk * gqxz[3] + qyzk 
			    * gqyz[3]) * 2.) + ck * (qxxi * gc[12] + qyyi * 
			    gc[17] + qzzi * gc[19] + (qxyi * gc[14] + qxzi * 
			    gc[15] + qyzi * gc[18]) * 2.) - uxi * (qxxk * 
			    gqxx[6] + qyyk * gqyy[6] + qzzk * gqzz[6] + (qxyk 
			    * gqxy[6] + qxzk * gqxz[6] + qyzk * gqyz[6]) * 2.)
			     - uyi * (qxxk * gqxx[8] + qyyk * gqyy[8] + qzzk *
			     gqzz[8] + (qxyk * gqxy[8] + qxzk * gqxz[8] + 
			    qyzk * gqyz[8]) * 2.) - uzi * (qxxk * gqxx[9] + 
			    qyyk * gqyy[9] + qzzk * gqzz[9] + (qxyk * gqxy[9] 
			    + qxzk * gqxz[9] + qyzk * gqyz[9]) * 2.) + uxk * (
			    qxxi * gux[12] + qyyi * gux[17] + qzzi * gux[19] 
			    + (qxyi * gux[14] + qxzi * gux[15] + qyzi * gux[
			    18]) * 2.) + uyk * (qxxi * guy[12] + qyyi * guy[
			    17] + qzzi * guy[19] + (qxyi * guy[14] + qxzi * 
			    guy[15] + qyzi * guy[18]) * 2.) + uzk * (qxxi * 
			    guz[12] + qyyi * guz[17] + qzzi * guz[19] + (qxyi 
			    * guz[14] + qxzi * guz[15] + qyzi * guz[18]) * 2.)
			     + qxxi * (qxxk * gqxx[12] + qyyk * gqyy[12] + 
			    qzzk * gqzz[12] + (qxyk * gqxy[12] + qxzk * gqxz[
			    12] + qyzk * gqyz[12]) * 2.) + qyyi * (qxxk * 
			    gqxx[17] + qyyk * gqyy[17] + qzzk * gqzz[17] + (
			    qxyk * gqxy[17] + qxzk * gqxz[17] + qyzk * gqyz[
			    17]) * 2.) + qzzi * (qxxk * gqxx[19] + qyyk * 
			    gqyy[19] + qzzk * gqzz[19] + (qxyk * gqxy[19] + 
			    qxzk * gqxz[19] + qyzk * gqyz[19]) * 2.) + (qxyi *
			     (qxxk * gqxx[14] + qyyk * gqyy[14] + qzzk * gqzz[
			    14] + (qxyk * gqxy[14] + qxzk * gqxz[14] + qyzk * 
			    gqyz[14]) * 2.) + qxzi * (qxxk * gqxx[15] + qyyk *
			     gqyy[15] + qzzk * gqzz[15] + (qxyk * gqxy[15] + 
			    qxzk * gqxz[15] + qyzk * gqyz[15]) * 2.) + qyzi * 
			    (qxxk * gqxx[18] + qyyk * gqyy[18] + qzzk * gqzz[
			    18] + (qxyk * gqxy[18] + qxzk * gqxz[18] + qyzk * 
			    gqyz[18]) * 2.)) * 2.;
		    dedz = desymdz + (dewidz + dewkdz) * .5;

		    desymdr = ci * ck * gc[20] - (uxi * (uxk * gux[21] + uyk *
			     guy[21] + uzk * guz[21]) + uyi * (uxk * gux[22] 
			    + uyk * guy[22] + uzk * guz[22]) + uzi * (uxk * 
			    gux[23] + uyk * guy[23] + uzk * guz[23]));
		    dewidr = ci * (uxk * gc[21] + uyk * gc[22] + uzk * gc[23])
			     - ck * (uxi * gux[20] + uyi * guy[20] + uzi * 
			    guz[20]) + ci * (qxxk * gc[24] + qyyk * gc[27] + 
			    qzzk * gc[29] + (qxyk * gc[25] + qxzk * gc[26] + 
			    qyzk * gc[28]) * 2.) + ck * (qxxi * gqxx[20] + 
			    qyyi * gqyy[20] + qzzi * gqzz[20] + (qxyi * gqxy[
			    20] + qxzi * gqxz[20] + qyzi * gqyz[20]) * 2.) - 
			    uxi * (qxxk * gux[24] + qyyk * gux[27] + qzzk * 
			    gux[29] + (qxyk * gux[25] + qxzk * gux[26] + qyzk 
			    * gux[28]) * 2.) - uyi * (qxxk * guy[24] + qyyk * 
			    guy[27] + qzzk * guy[29] + (qxyk * guy[25] + qxzk 
			    * guy[26] + qyzk * guy[28]) * 2.) - uzi * (qxxk * 
			    guz[24] + qyyk * guz[27] + qzzk * guz[29] + (qxyk 
			    * guz[25] + qxzk * guz[26] + qyzk * guz[28]) * 2.)
			     + uxk * (qxxi * gqxx[21] + qyyi * gqyy[21] + 
			    qzzi * gqzz[21] + (qxyi * gqxy[21] + qxzi * gqxz[
			    21] + qyzi * gqyz[21]) * 2.) + uyk * (qxxi * gqxx[
			    22] + qyyi * gqyy[22] + qzzi * gqzz[22] + (qxyi * 
			    gqxy[22] + qxzi * gqxz[22] + qyzi * gqyz[22]) * 
			    2.) + uzk * (qxxi * gqxx[23] + qyyi * gqyy[23] + 
			    qzzi * gqzz[23] + (qxyi * gqxy[23] + qxzi * gqxz[
			    23] + qyzi * gqyz[23]) * 2.) + qxxi * (qxxk * 
			    gqxx[24] + qyyk * gqxx[27] + qzzk * gqxx[29] + (
			    qxyk * gqxx[25] + qxzk * gqxx[26] + qyzk * gqxx[
			    28]) * 2.) + qyyi * (qxxk * gqyy[24] + qyyk * 
			    gqyy[27] + qzzk * gqyy[29] + (qxyk * gqyy[25] + 
			    qxzk * gqyy[26] + qyzk * gqyy[28]) * 2.) + qzzi * 
			    (qxxk * gqzz[24] + qyyk * gqzz[27] + qzzk * gqzz[
			    29] + (qxyk * gqzz[25] + qxzk * gqzz[26] + qyzk * 
			    gqzz[28]) * 2.) + (qxyi * (qxxk * gqxy[24] + qyyk 
			    * gqxy[27] + qzzk * gqxy[29] + (qxyk * gqxy[25] + 
			    qxzk * gqxy[26] + qyzk * gqxy[28]) * 2.) + qxzi * 
			    (qxxk * gqxz[24] + qyyk * gqxz[27] + qzzk * gqxz[
			    29] + (qxyk * gqxz[25] + qxzk * gqxz[26] + qyzk * 
			    gqxz[28]) * 2.) + qyzi * (qxxk * gqyz[24] + qyyk *
			     gqyz[27] + qzzk * gqyz[29] + (qxyk * gqyz[25] + 
			    qxzk * gqyz[26] + qyzk * gqyz[28]) * 2.)) * 2.;
		    dewkdr = ci * (uxk * gux[20] + uyk * guy[20] + uzk * guz[
			    20]) - ck * (uxi * gc[21] + uyi * gc[22] + uzi * 
			    gc[23]) + ci * (qxxk * gqxx[20] + qyyk * gqyy[20] 
			    + qzzk * gqzz[20] + (qxyk * gqxy[20] + qxzk * 
			    gqxz[20] + qyzk * gqyz[20]) * 2.) + ck * (qxxi * 
			    gc[24] + qyyi * gc[27] + qzzi * gc[29] + (qxyi * 
			    gc[25] + qxzi * gc[26] + qyzi * gc[28]) * 2.) - 
			    uxi * (qxxk * gqxx[21] + qyyk * gqyy[21] + qzzk * 
			    gqzz[21] + (qxyk * gqxy[21] + qxzk * gqxz[21] + 
			    qyzk * gqyz[21]) * 2.) - uyi * (qxxk * gqxx[22] + 
			    qyyk * gqyy[22] + qzzk * gqzz[22] + (qxyk * gqxy[
			    22] + qxzk * gqxz[22] + qyzk * gqyz[22]) * 2.) - 
			    uzi * (qxxk * gqxx[23] + qyyk * gqyy[23] + qzzk * 
			    gqzz[23] + (qxyk * gqxy[23] + qxzk * gqxz[23] + 
			    qyzk * gqyz[23]) * 2.) + uxk * (qxxi * gux[24] + 
			    qyyi * gux[27] + qzzi * gux[29] + (qxyi * gux[25] 
			    + qxzi * gux[26] + qyzi * gux[28]) * 2.) + uyk * (
			    qxxi * guy[24] + qyyi * guy[27] + qzzi * guy[29] 
			    + (qxyi * guy[25] + qxzi * guy[26] + qyzi * guy[
			    28]) * 2.) + uzk * (qxxi * guz[24] + qyyi * guz[
			    27] + qzzi * guz[29] + (qxyi * guz[25] + qxzi * 
			    guz[26] + qyzi * guz[28]) * 2.) + qxxi * (qxxk * 
			    gqxx[24] + qyyk * gqyy[24] + qzzk * gqzz[24] + (
			    qxyk * gqxy[24] + qxzk * gqxz[24] + qyzk * gqyz[
			    24]) * 2.) + qyyi * (qxxk * gqxx[27] + qyyk * 
			    gqyy[27] + qzzk * gqzz[27] + (qxyk * gqxy[27] + 
			    qxzk * gqxz[27] + qyzk * gqyz[27]) * 2.) + qzzi * 
			    (qxxk * gqxx[29] + qyyk * gqyy[29] + qzzk * gqzz[
			    29] + (qxyk * gqxy[29] + qxzk * gqxz[29] + qyzk * 
			    gqyz[29]) * 2.) + (qxyi * (qxxk * gqxx[25] + qyyk 
			    * gqyy[25] + qzzk * gqzz[25] + (qxyk * gqxy[25] + 
			    qxzk * gqxz[25] + qyzk * gqyz[25]) * 2.) + qxzi * 
			    (qxxk * gqxx[26] + qyyk * gqyy[26] + qzzk * gqzz[
			    26] + (qxyk * gqxy[26] + qxzk * gqxz[26] + qyzk * 
			    gqyz[26]) * 2.) + qyzi * (qxxk * gqxx[28] + qyyk *
			     gqyy[28] + qzzk * gqzz[28] + (qxyk * gqxy[28] + 
			    qxzk * gqxz[28] + qyzk * gqyz[28]) * 2.)) * 2.;
		    dsumdr = desymdr + (dewidr + dewkdr) * .5;
		    drbi = rbk * dsumdr;
		    drbk = rbi * dsumdr;

/*     torque on permanent dipoles due to permanent reaction field */

		    if (i__ != k) {
			fid[0] = uxk * gux[1] + uyk * gux[2] + uzk * gux[3] + 
				(ck * gux[0] + qxxk * gux[4] + qyyk * gux[7] 
				+ qzzk * gux[9] + (qxyk * gux[5] + qxzk * gux[
				6] + qyzk * gux[8]) * 2. + ck * gc[1] + qxxk *
				 gqxx[1] + qyyk * gqyy[1] + qzzk * gqzz[1] + (
				qxyk * gqxy[1] + qxzk * gqxz[1] + qyzk * gqyz[
				1]) * 2.) * .5;
			fid[1] = uxk * guy[1] + uyk * guy[2] + uzk * guy[3] + 
				(ck * guy[0] + qxxk * guy[4] + qyyk * guy[7] 
				+ qzzk * guy[9] + (qxyk * guy[5] + qxzk * guy[
				6] + qyzk * guy[8]) * 2. + ck * gc[2] + qxxk *
				 gqxx[2] + qyyk * gqyy[2] + qzzk * gqzz[2] + (
				qxyk * gqxy[2] + qxzk * gqxz[2] + qyzk * gqyz[
				2]) * 2.) * .5;
			fid[2] = uxk * guz[1] + uyk * guz[2] + uzk * guz[3] + 
				(ck * guz[0] + qxxk * guz[4] + qyyk * guz[7] 
				+ qzzk * guz[9] + (qxyk * guz[5] + qxzk * guz[
				6] + qyzk * guz[8]) * 2. + ck * gc[3] + qxxk *
				 gqxx[3] + qyyk * gqyy[3] + qzzk * gqzz[3] + (
				qxyk * gqxy[3] + qxzk * gqxz[3] + qyzk * gqyz[
				3]) * 2.) * .5;
			fkd[0] = uxi * gux[1] + uyi * gux[2] + uzi * gux[3] - 
				(ci * gux[0] + qxxi * gux[4] + qyyi * gux[7] 
				+ qzzi * gux[9] + (qxyi * gux[5] + qxzi * gux[
				6] + qyzi * gux[8]) * 2. + ci * gc[1] + qxxi *
				 gqxx[1] + qyyi * gqyy[1] + qzzi * gqzz[1] + (
				qxyi * gqxy[1] + qxzi * gqxz[1] + qyzi * gqyz[
				1]) * 2.) * .5;
			fkd[1] = uxi * guy[1] + uyi * guy[2] + uzi * guy[3] - 
				(ci * guy[0] + qxxi * guy[4] + qyyi * guy[7] 
				+ qzzi * guy[9] + (qxyi * guy[5] + qxzi * guy[
				6] + qyzi * guy[8]) * 2. + ci * gc[2] + qxxi *
				 gqxx[2] + qyyi * gqyy[2] + qzzi * gqzz[2] + (
				qxyi * gqxy[2] + qxzi * gqxz[2] + qyzi * gqyz[
				2]) * 2.) * .5;
			fkd[2] = uxi * guz[1] + uyi * guz[2] + uzi * guz[3] - 
				(ci * guz[0] + qxxi * guz[4] + qyyi * guz[7] 
				+ qzzi * guz[9] + (qxyi * guz[5] + qxzi * guz[
				6] + qyzi * guz[8]) * 2. + ci * gc[3] + qxxi *
				 gqxx[3] + qyyi * gqyy[3] + qzzi * gqzz[3] + (
				qxyi * gqxy[3] + qxzi * gqxz[3] + qyzi * gqyz[
				3]) * 2.) * .5;
			trq_ref(1, i__) = trq_ref(1, i__) + uyi * fid[2] - 
				uzi * fid[1];
			trq_ref(2, i__) = trq_ref(2, i__) + uzi * fid[0] - 
				uxi * fid[2];
			trq_ref(3, i__) = trq_ref(3, i__) + uxi * fid[1] - 
				uyi * fid[0];
			trq_ref(1, k) = trq_ref(1, k) + uyk * fkd[2] - uzk * 
				fkd[1];
			trq_ref(2, k) = trq_ref(2, k) + uzk * fkd[0] - uxk * 
				fkd[2];
			trq_ref(3, k) = trq_ref(3, k) + uxk * fkd[1] - uyk * 
				fkd[0];

/*     torque on quadrupoles due to permanent reaction field gradient */

			fidg_ref(1, 1) = (ck * gqxx[0] + uxk * gqxx[1] + uyk *
				 gqxx[2] + uzk * gqxx[3] + qxxk * gqxx[4] + 
				qyyk * gqxx[7] + qzzk * gqxx[9] + (qxyk * 
				gqxx[5] + qxzk * gqxx[6] + qyzk * gqxx[8]) * 
				2. + ck * gc[4] + uxk * gux[4] + uyk * guy[4] 
				+ uzk * guz[4] + qxxk * gqxx[4] + qyyk * gqyy[
				4] + qzzk * gqzz[4] + (qxyk * gqxy[4] + qxzk *
				 gqxz[4] + qyzk * gqyz[4]) * 2.) * -.5;
			fidg_ref(1, 2) = (ck * gqxy[0] + uxk * gqxy[1] + uyk *
				 gqxy[2] + uzk * gqxy[3] + qxxk * gqxy[4] + 
				qyyk * gqxy[7] + qzzk * gqxy[9] + (qxyk * 
				gqxy[5] + qxzk * gqxy[6] + qyzk * gqxy[8]) * 
				2. + ck * gc[5] + uxk * gux[5] + uyk * guy[5] 
				+ uzk * guz[5] + qxxk * gqxx[5] + qyyk * gqyy[
				5] + qzzk * gqzz[5] + (qxyk * gqxy[5] + qxzk *
				 gqxz[5] + qyzk * gqyz[5]) * 2.) * -.5;
			fidg_ref(1, 3) = (ck * gqxz[0] + uxk * gqxz[1] + uyk *
				 gqxz[2] + uzk * gqxz[3] + qxxk * gqxz[4] + 
				qyyk * gqxz[7] + qzzk * gqxz[9] + (qxyk * 
				gqxz[5] + qxzk * gqxz[6] + qyzk * gqxz[8]) * 
				2. + ck * gc[6] + uxk * gux[6] + uyk * guy[6] 
				+ uzk * guz[6] + qxxk * gqxx[6] + qyyk * gqyy[
				6] + qzzk * gqzz[6] + (qxyk * gqxy[6] + qxzk *
				 gqxz[6] + qyzk * gqyz[6]) * 2.) * -.5;
			fidg_ref(2, 2) = (ck * gqyy[0] + uxk * gqyy[1] + uyk *
				 gqyy[2] + uzk * gqyy[3] + qxxk * gqyy[4] + 
				qyyk * gqyy[7] + qzzk * gqyy[9] + (qxyk * 
				gqyy[5] + qxzk * gqyy[6] + qyzk * gqyy[8]) * 
				2. + ck * gc[7] + uxk * gux[7] + uyk * guy[7] 
				+ uzk * guz[7] + qxxk * gqxx[7] + qyyk * gqyy[
				7] + qzzk * gqzz[7] + (qxyk * gqxy[7] + qxzk *
				 gqxz[7] + qyzk * gqyz[7]) * 2.) * -.5;
			fidg_ref(2, 3) = (ck * gqyz[0] + uxk * gqyz[1] + uyk *
				 gqyz[2] + uzk * gqyz[3] + qxxk * gqyz[4] + 
				qyyk * gqyz[7] + qzzk * gqyz[9] + (qxyk * 
				gqyz[5] + qxzk * gqyz[6] + qyzk * gqyz[8]) * 
				2. + ck * gc[8] + uxk * gux[8] + uyk * guy[8] 
				+ uzk * guz[8] + qxxk * gqxx[8] + qyyk * gqyy[
				8] + qzzk * gqzz[8] + (qxyk * gqxy[8] + qxzk *
				 gqxz[8] + qyzk * gqyz[8]) * 2.) * -.5;
			fidg_ref(3, 3) = (ck * gqzz[0] + uxk * gqzz[1] + uyk *
				 gqzz[2] + uzk * gqzz[3] + qxxk * gqzz[4] + 
				qyyk * gqzz[7] + qzzk * gqzz[9] + (qxyk * 
				gqzz[5] + qxzk * gqzz[6] + qyzk * gqzz[8]) * 
				2. + ck * gc[9] + uxk * gux[9] + uyk * guy[9] 
				+ uzk * guz[9] + qxxk * gqxx[9] + qyyk * gqyy[
				9] + qzzk * gqzz[9] + (qxyk * gqxy[9] + qxzk *
				 gqxz[9] + qyzk * gqyz[9]) * 2.) * -.5;
			fidg_ref(2, 1) = fidg_ref(1, 2);
			fidg_ref(3, 1) = fidg_ref(1, 3);
			fidg_ref(3, 2) = fidg_ref(2, 3);
			fkdg_ref(1, 1) = (ci * gqxx[0] - uxi * gqxx[1] - uyi *
				 gqxx[2] - uzi * gqxx[3] + qxxi * gqxx[4] + 
				qyyi * gqxx[7] + qzzi * gqxx[9] + (qxyi * 
				gqxx[5] + qxzi * gqxx[6] + qyzi * gqxx[8]) * 
				2. + ci * gc[4] - uxi * gux[4] - uyi * guy[4] 
				- uzi * guz[4] + qxxi * gqxx[4] + qyyi * gqyy[
				4] + qzzi * gqzz[4] + (qxyi * gqxy[4] + qxzi *
				 gqxz[4] + qyzi * gqyz[4]) * 2.) * -.5;
			fkdg_ref(1, 2) = (ci * gqxy[0] - uxi * gqxy[1] - uyi *
				 gqxy[2] - uzi * gqxy[3] + qxxi * gqxy[4] + 
				qyyi * gqxy[7] + qzzi * gqxy[9] + (qxyi * 
				gqxy[5] + qxzi * gqxy[6] + qyzi * gqxy[8]) * 
				2. + ci * gc[5] - uxi * gux[5] - uyi * guy[5] 
				- uzi * guz[5] + qxxi * gqxx[5] + qyyi * gqyy[
				5] + qzzi * gqzz[5] + (qxyi * gqxy[5] + qxzi *
				 gqxz[5] + qyzi * gqyz[5]) * 2.) * -.5;
			fkdg_ref(1, 3) = (ci * gqxz[0] - uxi * gqxz[1] - uyi *
				 gqxz[2] - uzi * gqxz[3] + qxxi * gqxz[4] + 
				qyyi * gqxz[7] + qzzi * gqxz[9] + (qxyi * 
				gqxz[5] + qxzi * gqxz[6] + qyzi * gqxz[8]) * 
				2. + ci * gc[6] - uxi * gux[6] - uyi * guy[6] 
				- uzi * guz[6] + qxxi * gqxx[6] + qyyi * gqyy[
				6] + qzzi * gqzz[6] + (qxyi * gqxy[6] + qxzi *
				 gqxz[6] + qyzi * gqyz[6]) * 2.) * -.5;
			fkdg_ref(2, 2) = (ci * gqyy[0] - uxi * gqyy[1] - uyi *
				 gqyy[2] - uzi * gqyy[3] + qxxi * gqyy[4] + 
				qyyi * gqyy[7] + qzzi * gqyy[9] + (qxyi * 
				gqyy[5] + qxzi * gqyy[6] + qyzi * gqyy[8]) * 
				2. + ci * gc[7] - uxi * gux[7] - uyi * guy[7] 
				- uzi * guz[7] + qxxi * gqxx[7] + qyyi * gqyy[
				7] + qzzi * gqzz[7] + (qxyi * gqxy[7] + qxzi *
				 gqxz[7] + qyzi * gqyz[7]) * 2.) * -.5;
			fkdg_ref(2, 3) = (ci * gqyz[0] - uxi * gqyz[1] - uyi *
				 gqyz[2] - uzi * gqyz[3] + qxxi * gqyz[4] + 
				qyyi * gqyz[7] + qzzi * gqyz[9] + (qxyi * 
				gqyz[5] + qxzi * gqyz[6] + qyzi * gqyz[8]) * 
				2. + ci * gc[8] - uxi * gux[8] - uyi * guy[8] 
				- uzi * guz[8] + qxxi * gqxx[8] + qyyi * gqyy[
				8] + qzzi * gqzz[8] + (qxyi * gqxy[8] + qxzi *
				 gqxz[8] + qyzi * gqyz[8]) * 2.) * -.5;
			fkdg_ref(3, 3) = (ci * gqzz[0] - uxi * gqzz[1] - uyi *
				 gqzz[2] - uzi * gqzz[3] + qxxi * gqzz[4] + 
				qyyi * gqzz[7] + qzzi * gqzz[9] + (qxyi * 
				gqzz[5] + qxzi * gqzz[6] + qyzi * gqzz[8]) * 
				2. + ci * gc[9] - uxi * gux[9] - uyi * guy[9] 
				- uzi * guz[9] + qxxi * gqxx[9] + qyyi * gqyy[
				9] + qzzi * gqzz[9] + (qxyi * gqxy[9] + qxzi *
				 gqxz[9] + qyzi * gqyz[9]) * 2.) * -.5;
			fkdg_ref(2, 1) = fkdg_ref(1, 2);
			fkdg_ref(3, 1) = fkdg_ref(1, 3);
			fkdg_ref(3, 2) = fkdg_ref(2, 3);
			trq_ref(1, i__) = trq_ref(1, i__) + (qxyi * fidg_ref(
				1, 3) + qyyi * fidg_ref(2, 3) + qyzi * 
				fidg_ref(3, 3) - qxzi * fidg_ref(1, 2) - qyzi 
				* fidg_ref(2, 2) - qzzi * fidg_ref(3, 2)) * 
				2.;
			trq_ref(2, i__) = trq_ref(2, i__) + (qxzi * fidg_ref(
				1, 1) + qyzi * fidg_ref(2, 1) + qzzi * 
				fidg_ref(3, 1) - qxxi * fidg_ref(1, 3) - qxyi 
				* fidg_ref(2, 3) - qxzi * fidg_ref(3, 3)) * 
				2.;
			trq_ref(3, i__) = trq_ref(3, i__) + (qxxi * fidg_ref(
				1, 2) + qxyi * fidg_ref(2, 2) + qxzi * 
				fidg_ref(3, 2) - qxyi * fidg_ref(1, 1) - qyyi 
				* fidg_ref(2, 1) - qyzi * fidg_ref(3, 1)) * 
				2.;
			trq_ref(1, k) = trq_ref(1, k) + (qxyk * fkdg_ref(1, 3)
				 + qyyk * fkdg_ref(2, 3) + qyzk * fkdg_ref(3, 
				3) - qxzk * fkdg_ref(1, 2) - qyzk * fkdg_ref(
				2, 2) - qzzk * fkdg_ref(3, 2)) * 2.;
			trq_ref(2, k) = trq_ref(2, k) + (qxzk * fkdg_ref(1, 1)
				 + qyzk * fkdg_ref(2, 1) + qzzk * fkdg_ref(3, 
				1) - qxxk * fkdg_ref(1, 3) - qxyk * fkdg_ref(
				2, 3) - qxzk * fkdg_ref(3, 3)) * 2.;
			trq_ref(3, k) = trq_ref(3, k) + (qxxk * fkdg_ref(1, 2)
				 + qxyk * fkdg_ref(2, 2) + qxzk * fkdg_ref(3, 
				2) - qxyk * fkdg_ref(1, 1) - qyyk * fkdg_ref(
				2, 1) - qyzk * fkdg_ref(3, 1)) * 2.;
		    }

/*     electrostatic solvation energy of the permanent multipoles in */
/*     the GK reaction potential of the induced dipoles */

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

/*     electrostatic solvation free energy gradient of the permanent */
/*     multipoles in the reaction potential of the induced dipoles */

		    dpsymdx = -uxi * (sxk * gux[4] + syk * guy[4] + szk * guz[
			    4]) - uyi * (sxk * gux[5] + syk * guy[5] + szk * 
			    guz[5]) - uzi * (sxk * gux[6] + syk * guy[6] + 
			    szk * guz[6]) - uxk * (sxi * gux[4] + syi * guy[4]
			     + szi * guz[4]) - uyk * (sxi * gux[5] + syi * 
			    guy[5] + szi * guz[5]) - uzk * (sxi * gux[6] + 
			    syi * guy[6] + szi * guz[6]);
		    dpwidx = ci * (sxk * gc[4] + syk * gc[5] + szk * gc[6]) - 
			    ck * (sxi * gux[1] + syi * guy[1] + szi * guz[1]) 
			    - sxi * (qxxk * gux[10] + qyyk * gux[13] + qzzk * 
			    gux[15] + (qxyk * gux[11] + qxzk * gux[12] + qyzk 
			    * gux[14]) * 2.) - syi * (qxxk * guy[10] + qyyk * 
			    guy[13] + qzzk * guy[15] + (qxyk * guy[11] + qxzk 
			    * guy[12] + qyzk * guy[14]) * 2.) - szi * (qxxk * 
			    guz[10] + qyyk * guz[13] + qzzk * guz[15] + (qxyk 
			    * guz[11] + qxzk * guz[12] + qyzk * guz[14]) * 2.)
			     + sxk * (qxxi * gqxx[4] + qyyi * gqyy[4] + qzzi *
			     gqzz[4] + (qxyi * gqxy[4] + qxzi * gqxz[4] + 
			    qyzi * gqyz[4]) * 2.) + syk * (qxxi * gqxx[5] + 
			    qyyi * gqyy[5] + qzzi * gqzz[5] + (qxyi * gqxy[5] 
			    + qxzi * gqxz[5] + qyzi * gqyz[5]) * 2.) + szk * (
			    qxxi * gqxx[6] + qyyi * gqyy[6] + qzzi * gqzz[6] 
			    + (qxyi * gqxy[6] + qxzi * gqxz[6] + qyzi * gqyz[
			    6]) * 2.);
		    dpwkdx = ci * (sxk * gux[1] + syk * guy[1] + szk * guz[1])
			     - ck * (sxi * gc[4] + syi * gc[5] + szi * gc[6]) 
			    - sxi * (qxxk * gqxx[4] + qyyk * gqyy[4] + qzzk * 
			    gqzz[4] + (qxyk * gqxy[4] + qxzk * gqxz[4] + qyzk 
			    * gqyz[4]) * 2.) - syi * (qxxk * gqxx[5] + qyyk * 
			    gqyy[5] + qzzk * gqzz[5] + (qxyk * gqxy[5] + qxzk 
			    * gqxz[5] + qyzk * gqyz[5]) * 2.) - szi * (qxxk * 
			    gqxx[6] + qyyk * gqyy[6] + qzzk * gqzz[6] + (qxyk 
			    * gqxy[6] + qxzk * gqxz[6] + qyzk * gqyz[6]) * 2.)
			     + sxk * (qxxi * gux[10] + qyyi * gux[13] + qzzi *
			     gux[15] + (qxyi * gux[11] + qxzi * gux[12] + 
			    qyzi * gux[14]) * 2.) + syk * (qxxi * guy[10] + 
			    qyyi * guy[13] + qzzi * guy[15] + (qxyi * guy[11] 
			    + qxzi * guy[12] + qyzi * guy[14]) * 2.) + szk * (
			    qxxi * guz[10] + qyyi * guz[13] + qzzi * guz[15] 
			    + (qxyi * guz[11] + qxzi * guz[12] + qyzi * guz[
			    14]) * 2.);
		    dpdx = (dpsymdx + (dpwidx + dpwkdx) * .5) * .5;
		    dpsymdy = -uxi * (sxk * gux[5] + syk * guy[5] + szk * guz[
			    5]) - uyi * (sxk * gux[7] + syk * guy[7] + szk * 
			    guz[7]) - uzi * (sxk * gux[8] + syk * guy[8] + 
			    szk * guz[8]) - uxk * (sxi * gux[5] + syi * guy[5]
			     + szi * guz[5]) - uyk * (sxi * gux[7] + syi * 
			    guy[7] + szi * guz[7]) - uzk * (sxi * gux[8] + 
			    syi * guy[8] + szi * guz[8]);
		    dpwidy = ci * (sxk * gc[5] + syk * gc[7] + szk * gc[8]) - 
			    ck * (sxi * gux[2] + syi * guy[2] + szi * guz[2]) 
			    - sxi * (qxxk * gux[11] + qyyk * gux[16] + qzzk * 
			    gux[18] + (qxyk * gux[13] + qxzk * gux[14] + qyzk 
			    * gux[17]) * 2.) - syi * (qxxk * guy[11] + qyyk * 
			    guy[16] + qzzk * guy[18] + (qxyk * guy[13] + qxzk 
			    * guy[14] + qyzk * guy[17]) * 2.) - szi * (qxxk * 
			    guz[11] + qyyk * guz[16] + qzzk * guz[18] + (qxyk 
			    * guz[13] + qxzk * guz[14] + qyzk * guz[17]) * 2.)
			     + sxk * (qxxi * gqxx[5] + qyyi * gqyy[5] + qzzi *
			     gqzz[5] + (qxyi * gqxy[5] + qxzi * gqxz[5] + 
			    qyzi * gqyz[5]) * 2.) + syk * (qxxi * gqxx[7] + 
			    qyyi * gqyy[7] + qzzi * gqzz[7] + (qxyi * gqxy[7] 
			    + qxzi * gqxz[7] + qyzi * gqyz[7]) * 2.) + szk * (
			    qxxi * gqxx[8] + qyyi * gqyy[8] + qzzi * gqzz[8] 
			    + (qxyi * gqxy[8] + qxzi * gqxz[8] + qyzi * gqyz[
			    8]) * 2.);
		    dpwkdy = ci * (sxk * gux[2] + syk * guy[2] + szk * guz[2])
			     - ck * (sxi * gc[5] + syi * gc[7] + szi * gc[8]) 
			    - sxi * (qxxk * gqxx[5] + qyyk * gqyy[5] + qzzk * 
			    gqzz[5] + (qxyk * gqxy[5] + qxzk * gqxz[5] + qyzk 
			    * gqyz[5]) * 2.) - syi * (qxxk * gqxx[7] + qyyk * 
			    gqyy[7] + qzzk * gqzz[7] + (qxyk * gqxy[7] + qxzk 
			    * gqxz[7] + qyzk * gqyz[7]) * 2.) - szi * (qxxk * 
			    gqxx[8] + qyyk * gqyy[8] + qzzk * gqzz[8] + (qxyk 
			    * gqxy[8] + qxzk * gqxz[8] + qyzk * gqyz[8]) * 2.)
			     + sxk * (qxxi * gux[11] + qyyi * gux[16] + qzzi *
			     gux[18] + (qxyi * gux[13] + qxzi * gux[14] + 
			    qyzi * gux[17]) * 2.) + syk * (qxxi * guy[11] + 
			    qyyi * guy[16] + qzzi * guy[18] + (qxyi * guy[13] 
			    + qxzi * guy[14] + qyzi * guy[17]) * 2.) + szk * (
			    qxxi * guz[11] + qyyi * guz[16] + qzzi * guz[18] 
			    + (qxyi * guz[13] + qxzi * guz[14] + qyzi * guz[
			    17]) * 2.);
		    dpdy = (dpsymdy + (dpwidy + dpwkdy) * .5) * .5;
		    dpsymdz = -uxi * (sxk * gux[6] + syk * guy[6] + szk * guz[
			    6]) - uyi * (sxk * gux[8] + syk * guy[8] + szk * 
			    guz[8]) - uzi * (sxk * gux[9] + syk * guy[9] + 
			    szk * guz[9]) - uxk * (sxi * gux[6] + syi * guy[6]
			     + szi * guz[6]) - uyk * (sxi * gux[8] + syi * 
			    guy[8] + szi * guz[8]) - uzk * (sxi * gux[9] + 
			    syi * guy[9] + szi * guz[9]);
		    dpwidz = ci * (sxk * gc[6] + syk * gc[8] + szk * gc[9]) - 
			    ck * (sxi * gux[3] + syi * guy[3] + szi * guz[3]) 
			    - sxi * (qxxk * gux[12] + qyyk * gux[17] + qzzk * 
			    gux[19] + (qxyk * gux[14] + qxzk * gux[15] + qyzk 
			    * gux[18]) * 2.) - syi * (qxxk * guy[12] + qyyk * 
			    guy[17] + qzzk * guy[19] + (qxyk * guy[14] + qxzk 
			    * guy[15] + qyzk * guy[18]) * 2.) - szi * (qxxk * 
			    guz[12] + qyyk * guz[17] + qzzk * guz[19] + (qxyk 
			    * guz[14] + qxzk * guz[15] + qyzk * guz[18]) * 2.)
			     + sxk * (qxxi * gqxx[6] + qyyi * gqyy[6] + qzzi *
			     gqzz[6] + (qxyi * gqxy[6] + qxzi * gqxz[6] + 
			    qyzi * gqyz[6]) * 2.) + syk * (qxxi * gqxx[8] + 
			    qyyi * gqyy[8] + qzzi * gqzz[8] + (qxyi * gqxy[8] 
			    + qxzi * gqxz[8] + qyzi * gqyz[8]) * 2.) + szk * (
			    qxxi * gqxx[9] + qyyi * gqyy[9] + qzzi * gqzz[9] 
			    + (qxyi * gqxy[9] + qxzi * gqxz[9] + qyzi * gqyz[
			    9]) * 2.);
		    dpwkdz = ci * (sxk * gux[3] + syk * guy[3] + szk * guz[3])
			     - ck * (sxi * gc[6] + syi * gc[8] + szi * gc[9]) 
			    - sxi * (qxxk * gqxx[6] + qyyk * gqyy[6] + qzzk * 
			    gqzz[6] + (qxyk * gqxy[6] + qxzk * gqxz[6] + qyzk 
			    * gqyz[6]) * 2.) - syi * (qxxk * gqxx[8] + qyyk * 
			    gqyy[8] + qzzk * gqzz[8] + (qxyk * gqxy[8] + qxzk 
			    * gqxz[8] + qyzk * gqyz[8]) * 2.) - szi * (qxxk * 
			    gqxx[9] + qyyk * gqyy[9] + qzzk * gqzz[9] + (qxyk 
			    * gqxy[9] + qxzk * gqxz[9] + qyzk * gqyz[9]) * 2.)
			     + sxk * (qxxi * gux[12] + qyyi * gux[17] + qzzi *
			     gux[19] + (qxyi * gux[14] + qxzi * gux[15] + 
			    qyzi * gux[18]) * 2.) + syk * (qxxi * guy[12] + 
			    qyyi * guy[17] + qzzi * guy[19] + (qxyi * guy[14] 
			    + qxzi * guy[15] + qyzi * guy[18]) * 2.) + szk * (
			    qxxi * guz[12] + qyyi * guz[17] + qzzi * guz[19] 
			    + (qxyi * guz[14] + qxzi * guz[15] + qyzi * guz[
			    18]) * 2.);
		    dpdz = (dpsymdz + (dpwidz + dpwkdz) * .5) * .5;

/*     effective radii chain rule terms for the */
/*     electrostatic solvation free energy gradient of the permanent */
/*     multipoles in the reaction potential of the induced dipoles */

		    dsymdr = -uxi * (sxk * gux[21] + syk * guy[21] + szk * 
			    guz[21]) - uyi * (sxk * gux[22] + syk * guy[22] + 
			    szk * guz[22]) - uzi * (sxk * gux[23] + syk * guy[
			    23] + szk * guz[23]) - uxk * (sxi * gux[21] + syi 
			    * guy[21] + szi * guz[21]) - uyk * (sxi * gux[22] 
			    + syi * guy[22] + szi * guz[22]) - uzk * (sxi * 
			    gux[23] + syi * guy[23] + szi * guz[23]);
		    dwipdr = ci * (sxk * gc[21] + syk * gc[22] + szk * gc[23])
			     - ck * (sxi * gux[20] + syi * guy[20] + szi * 
			    guz[20]) - sxi * (qxxk * gux[24] + qyyk * gux[27] 
			    + qzzk * gux[29] + (qxyk * gux[25] + qxzk * gux[
			    26] + qyzk * gux[28]) * 2.) - syi * (qxxk * guy[
			    24] + qyyk * guy[27] + qzzk * guy[29] + (qxyk * 
			    guy[25] + qxzk * guy[26] + qyzk * guy[28]) * 2.) 
			    - szi * (qxxk * guz[24] + qyyk * guz[27] + qzzk * 
			    guz[29] + (qxyk * guz[25] + qxzk * guz[26] + qyzk 
			    * guz[28]) * 2.) + sxk * (qxxi * gqxx[21] + qyyi *
			     gqyy[21] + qzzi * gqzz[21] + (qxyi * gqxy[21] + 
			    qxzi * gqxz[21] + qyzi * gqyz[21]) * 2.) + syk * (
			    qxxi * gqxx[22] + qyyi * gqyy[22] + qzzi * gqzz[
			    22] + (qxyi * gqxy[22] + qxzi * gqxz[22] + qyzi * 
			    gqyz[22]) * 2.) + szk * (qxxi * gqxx[23] + qyyi * 
			    gqyy[23] + qzzi * gqzz[23] + (qxyi * gqxy[23] + 
			    qxzi * gqxz[23] + qyzi * gqyz[23]) * 2.);
		    dwkpdr = ci * (sxk * gux[20] + syk * guy[20] + szk * guz[
			    20]) - ck * (sxi * gc[21] + syi * gc[22] + szi * 
			    gc[23]) - sxi * (qxxk * gqxx[21] + qyyk * gqyy[21]
			     + qzzk * gqzz[21] + (qxyk * gqxy[21] + qxzk * 
			    gqxz[21] + qyzk * gqyz[21]) * 2.) - syi * (qxxk * 
			    gqxx[22] + qyyk * gqyy[22] + qzzk * gqzz[22] + (
			    qxyk * gqxy[22] + qxzk * gqxz[22] + qyzk * gqyz[
			    22]) * 2.) - szi * (qxxk * gqxx[23] + qyyk * gqyy[
			    23] + qzzk * gqzz[23] + (qxyk * gqxy[23] + qxzk * 
			    gqxz[23] + qyzk * gqyz[23]) * 2.) + sxk * (qxxi * 
			    gux[24] + qyyi * gux[27] + qzzi * gux[29] + (qxyi 
			    * gux[25] + qxzi * gux[26] + qyzi * gux[28]) * 2.)
			     + syk * (qxxi * guy[24] + qyyi * guy[27] + qzzi *
			     guy[29] + (qxyi * guy[25] + qxzi * guy[26] + 
			    qyzi * guy[28]) * 2.) + szk * (qxxi * guz[24] + 
			    qyyi * guz[27] + qzzi * guz[29] + (qxyi * guz[25] 
			    + qxzi * guz[26] + qyzi * guz[28]) * 2.);
		    dsumdr = dsymdr + (dwipdr + dwkpdr) * .5;
		    dpbi = rbk * .5 * dsumdr;
		    dpbk = rbi * .5 * dsumdr;

/*     mutual polarization electrostatic solvation free energy gradient */

		    if (s_cmp(polpot_1.poltyp, "MUTUAL", (ftnlen)6, (ftnlen)6)
			     == 0) {
			dpdx -= (dxi * (pxk * gux[4] + pyk * gux[5] + pzk * 
				gux[6]) + dyi * (pxk * guy[4] + pyk * guy[5] 
				+ pzk * guy[6]) + dzi * (pxk * guz[4] + pyk * 
				guz[5] + pzk * guz[6]) + dxk * (pxi * gux[4] 
				+ pyi * gux[5] + pzi * gux[6]) + dyk * (pxi * 
				guy[4] + pyi * guy[5] + pzi * guy[6]) + dzk * 
				(pxi * guz[4] + pyi * guz[5] + pzi * guz[6])) 
				* .5;
			dpdy -= (dxi * (pxk * gux[5] + pyk * gux[7] + pzk * 
				gux[8]) + dyi * (pxk * guy[5] + pyk * guy[7] 
				+ pzk * guy[8]) + dzi * (pxk * guz[5] + pyk * 
				guz[7] + pzk * guz[8]) + dxk * (pxi * gux[5] 
				+ pyi * gux[7] + pzi * gux[8]) + dyk * (pxi * 
				guy[5] + pyi * guy[7] + pzi * guy[8]) + dzk * 
				(pxi * guz[5] + pyi * guz[7] + pzi * guz[8])) 
				* .5;
			dpdz -= (dxi * (pxk * gux[6] + pyk * gux[8] + pzk * 
				gux[9]) + dyi * (pxk * guy[6] + pyk * guy[8] 
				+ pzk * guy[9]) + dzi * (pxk * guz[6] + pyk * 
				guz[8] + pzk * guz[9]) + dxk * (pxi * gux[6] 
				+ pyi * gux[8] + pzi * gux[9]) + dyk * (pxi * 
				guy[6] + pyi * guy[8] + pzi * guy[9]) + dzk * 
				(pxi * guz[6] + pyi * guz[8] + pzi * guz[9])) 
				* .5;
			duvdr = dxi * (pxk * gux[21] + pyk * gux[22] + pzk * 
				gux[23]) + dyi * (pxk * guy[21] + pyk * guy[
				22] + pzk * guy[23]) + dzi * (pxk * guz[21] + 
				pyk * guz[22] + pzk * guz[23]) + dxk * (pxi * 
				gux[21] + pyi * gux[22] + pzi * gux[23]) + 
				dyk * (pxi * guy[21] + pyi * guy[22] + pzi * 
				guy[23]) + dzk * (pxi * guz[21] + pyi * guz[
				22] + pzi * guz[23]);
			dpbi -= rbk * .5 * duvdr;
			dpbk -= rbi * .5 * duvdr;
		    }

/*     torque due to induced reaction field on permanent dipoles */

		    fid[0] = (sxk * gux[1] + syk * guy[1] + szk * guz[1]) * 
			    .5;
		    fid[1] = (sxk * gux[2] + syk * guy[2] + szk * guz[2]) * 
			    .5;
		    fid[2] = (sxk * gux[3] + syk * guy[3] + szk * guz[3]) * 
			    .5;
		    fkd[0] = (sxi * gux[1] + syi * guy[1] + szi * guz[1]) * 
			    .5;
		    fkd[1] = (sxi * gux[2] + syi * guy[2] + szi * guz[2]) * 
			    .5;
		    fkd[2] = (sxi * gux[3] + syi * guy[3] + szi * guz[3]) * 
			    .5;
		    if (i__ == k) {
			fid[0] *= .5;
			fid[1] *= .5;
			fid[2] *= .5;
			fkd[0] *= .5;
			fkd[1] *= .5;
			fkd[2] *= .5;
		    }
		    trqi_ref(1, i__) = trqi_ref(1, i__) + uyi * fid[2] - uzi *
			     fid[1];
		    trqi_ref(2, i__) = trqi_ref(2, i__) + uzi * fid[0] - uxi *
			     fid[2];
		    trqi_ref(3, i__) = trqi_ref(3, i__) + uxi * fid[1] - uyi *
			     fid[0];
		    trqi_ref(1, k) = trqi_ref(1, k) + uyk * fkd[2] - uzk * 
			    fkd[1];
		    trqi_ref(2, k) = trqi_ref(2, k) + uzk * fkd[0] - uxk * 
			    fkd[2];
		    trqi_ref(3, k) = trqi_ref(3, k) + uxk * fkd[1] - uyk * 
			    fkd[0];

/*     torque due to induced reaction field gradient on quadrupoles */

		    fidg_ref(1, 1) = (sxk * gqxx[1] + syk * gqxx[2] + szk * 
			    gqxx[3] + (sxk * gux[4] + syk * guy[4] + szk * 
			    guz[4])) * -.25;
		    fidg_ref(1, 2) = (sxk * gqxy[1] + syk * gqxy[2] + szk * 
			    gqxy[3] + (sxk * gux[5] + syk * guy[5] + szk * 
			    guz[5])) * -.25;
		    fidg_ref(1, 3) = (sxk * gqxz[1] + syk * gqxz[2] + szk * 
			    gqxz[3] + (sxk * gux[6] + syk * guy[6] + szk * 
			    guz[6])) * -.25;
		    fidg_ref(2, 2) = (sxk * gqyy[1] + syk * gqyy[2] + szk * 
			    gqyy[3] + (sxk * gux[7] + syk * guy[7] + szk * 
			    guz[7])) * -.25;
		    fidg_ref(2, 3) = (sxk * gqyz[1] + syk * gqyz[2] + szk * 
			    gqyz[3] + (sxk * gux[8] + syk * guy[8] + szk * 
			    guz[8])) * -.25;
		    fidg_ref(3, 3) = (sxk * gqzz[1] + syk * gqzz[2] + szk * 
			    gqzz[3] + (sxk * gux[9] + syk * guy[9] + szk * 
			    guz[9])) * -.25;
		    fidg_ref(2, 1) = fidg_ref(1, 2);
		    fidg_ref(3, 1) = fidg_ref(1, 3);
		    fidg_ref(3, 2) = fidg_ref(2, 3);
		    fkdg_ref(1, 1) = (sxi * gqxx[1] + syi * gqxx[2] + szi * 
			    gqxx[3] + (sxi * gux[4] + syi * guy[4] + szi * 
			    guz[4])) * .25;
		    fkdg_ref(1, 2) = (sxi * gqxy[1] + syi * gqxy[2] + szi * 
			    gqxy[3] + (sxi * gux[5] + syi * guy[5] + szi * 
			    guz[5])) * .25;
		    fkdg_ref(1, 3) = (sxi * gqxz[1] + syi * gqxz[2] + szi * 
			    gqxz[3] + (sxi * gux[6] + syi * guy[6] + szi * 
			    guz[6])) * .25;
		    fkdg_ref(2, 2) = (sxi * gqyy[1] + syi * gqyy[2] + szi * 
			    gqyy[3] + (sxi * gux[7] + syi * guy[7] + szi * 
			    guz[7])) * .25;
		    fkdg_ref(2, 3) = (sxi * gqyz[1] + syi * gqyz[2] + szi * 
			    gqyz[3] + (sxi * gux[8] + syi * guy[8] + szi * 
			    guz[8])) * .25;
		    fkdg_ref(3, 3) = (sxi * gqzz[1] + syi * gqzz[2] + szi * 
			    gqzz[3] + (sxi * gux[9] + syi * guy[9] + szi * 
			    guz[9])) * .25;
		    fkdg_ref(2, 1) = fkdg_ref(1, 2);
		    fkdg_ref(3, 1) = fkdg_ref(1, 3);
		    fkdg_ref(3, 2) = fkdg_ref(2, 3);
		    if (i__ == k) {
			fidg_ref(1, 1) = fidg_ref(1, 1) * .5;
			fidg_ref(1, 2) = fidg_ref(1, 2) * .5;
			fidg_ref(1, 3) = fidg_ref(1, 3) * .5;
			fidg_ref(2, 1) = fidg_ref(2, 1) * .5;
			fidg_ref(2, 2) = fidg_ref(2, 2) * .5;
			fidg_ref(2, 3) = fidg_ref(2, 3) * .5;
			fidg_ref(3, 1) = fidg_ref(3, 1) * .5;
			fidg_ref(3, 2) = fidg_ref(3, 2) * .5;
			fidg_ref(3, 3) = fidg_ref(3, 3) * .5;
			fkdg_ref(1, 1) = fkdg_ref(1, 1) * .5;
			fkdg_ref(1, 2) = fkdg_ref(1, 2) * .5;
			fkdg_ref(1, 3) = fkdg_ref(1, 3) * .5;
			fkdg_ref(2, 1) = fkdg_ref(2, 1) * .5;
			fkdg_ref(2, 2) = fkdg_ref(2, 2) * .5;
			fkdg_ref(2, 3) = fkdg_ref(2, 3) * .5;
			fkdg_ref(3, 1) = fkdg_ref(3, 1) * .5;
			fkdg_ref(3, 2) = fkdg_ref(3, 2) * .5;
			fkdg_ref(3, 3) = fkdg_ref(3, 3) * .5;
		    }
		    trqi_ref(1, i__) = trqi_ref(1, i__) + (qxyi * fidg_ref(1, 
			    3) + qyyi * fidg_ref(2, 3) + qyzi * fidg_ref(3, 3)
			     - qxzi * fidg_ref(1, 2) - qyzi * fidg_ref(2, 2) 
			    - qzzi * fidg_ref(3, 2)) * 2.;
		    trqi_ref(2, i__) = trqi_ref(2, i__) + (qxzi * fidg_ref(1, 
			    1) + qyzi * fidg_ref(2, 1) + qzzi * fidg_ref(3, 1)
			     - qxxi * fidg_ref(1, 3) - qxyi * fidg_ref(2, 3) 
			    - qxzi * fidg_ref(3, 3)) * 2.;
		    trqi_ref(3, i__) = trqi_ref(3, i__) + (qxxi * fidg_ref(1, 
			    2) + qxyi * fidg_ref(2, 2) + qxzi * fidg_ref(3, 2)
			     - qxyi * fidg_ref(1, 1) - qyyi * fidg_ref(2, 1) 
			    - qyzi * fidg_ref(3, 1)) * 2.;
		    trqi_ref(1, k) = trqi_ref(1, k) + (qxyk * fkdg_ref(1, 3) 
			    + qyyk * fkdg_ref(2, 3) + qyzk * fkdg_ref(3, 3) - 
			    qxzk * fkdg_ref(1, 2) - qyzk * fkdg_ref(2, 2) - 
			    qzzk * fkdg_ref(3, 2)) * 2.;
		    trqi_ref(2, k) = trqi_ref(2, k) + (qxzk * fkdg_ref(1, 1) 
			    + qyzk * fkdg_ref(2, 1) + qzzk * fkdg_ref(3, 1) - 
			    qxxk * fkdg_ref(1, 3) - qxyk * fkdg_ref(2, 3) - 
			    qxzk * fkdg_ref(3, 3)) * 2.;
		    trqi_ref(3, k) = trqi_ref(3, k) + (qxxk * fkdg_ref(1, 2) 
			    + qxyk * fkdg_ref(2, 2) + qxzk * fkdg_ref(3, 2) - 
			    qxyk * fkdg_ref(1, 1) - qyyk * fkdg_ref(2, 1) - 
			    qyzk * fkdg_ref(3, 1)) * 2.;

/*     total permanent and induced energies for this interaction */

		    e = esym + (ewi + ewk) * .5;
		    ei = (esymi + (ewii + ewki) * .5) * .5;

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
			dedx *= fgrp;
			dedy *= fgrp;
			dedz *= fgrp;
			drbi *= fgrp;
			drbk *= fgrp;
			ei *= fgrp;
			dpdx *= fgrp;
			dpdy *= fgrp;
			dpdz *= fgrp;
			dpbi *= fgrp;
			dpbk *= fgrp;
		    }

/*     increment the overall energy and derivative expressions */

		    if (i__ == k) {
			e *= .5;
			ei *= .5;
			energi_1.es = energi_1.es + e + ei;
			solute_1.drb[i__ - 1] += drbi;
			solute_1.drbp[i__ - 1] += dpbi;
		    } else {
			energi_1.es = energi_1.es + e + ei;
			des_ref(1, i__) = des_ref(1, i__) - dedx - dpdx;
			des_ref(2, i__) = des_ref(2, i__) - dedy - dpdy;
			des_ref(3, i__) = des_ref(3, i__) - dedz - dpdz;
			des_ref(1, k) = des_ref(1, k) + dedx + dpdx;
			des_ref(2, k) = des_ref(2, k) + dedy + dpdy;
			des_ref(3, k) = des_ref(3, k) + dedz + dpdz;
			solute_1.drb[i__ - 1] += drbi;
			solute_1.drb[k - 1] += drbk;
			solute_1.drbp[i__ - 1] += dpbi;
			solute_1.drbp[k - 1] += dpbk;

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

/*     increment the total intermolecular energy */

		    if (molcul_1.molcule[i__ - 1] != molcul_1.molcule[k - 1]) 
			    {
			inter_1.einter = inter_1.einter + e + ei;
		    }
		}
	    }
	}
    }
    torque2_(trq, deriv_1.des);
    torque2_(trqi, deriv_1.des);
    return 0;
} /* egk1a_ */

#undef uinps_ref
#undef uinds_ref
#undef rpole_ref
#undef trqi_ref
#undef fkdg_ref
#undef fidg_ref
#undef trq_ref
#undef vir_ref
#undef des_ref
#undef b_ref
#undef a_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine ediff1  --  correct vacuum to SCRF derivatives  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "ediff1" calculates the energy and derivatives of polarizing */
/*     the vacuum induced dipoles to their SCRF polarized values */


/* Subroutine */ int ediff1_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal f;
    static integer i__, j, k;
    static doublereal r__, r2;
    static integer ii, kk;
    static doublereal ei, ci, di[3];
    static integer ix;
    static doublereal qi[9];
    static integer iz, kx, kz;
    static doublereal ck, dk[3], qk[9], sc[10], xr, yr, zr, rr1, rr3, rr5, 
	    rr7, rr9, gfd, gfi[6], pdi, pti, qir[3], qkr[3], gli[7], sci[8], 
	    gti[6], dsc3, dsc5, dsc7, dsc9, psc3, psc5, psc7, psc9, damp, 
	    fdir[3], qidk[3], qkdi[3], glip[7], fgrp, qiuk[3], qkui[3], dixr[
	    3], dkxr[3], scip[8], trqi[75000]	/* was [3][25000] */;
    static logical usei, usek;
    static doublereal ddsc3[3], ddsc5[3], ddsc7[3], ftm1i[75000]	/* 
	    was [3][25000] */, ftm2i[3], ttm1i[75000]	/* was [3][25000] */, 
	    ttm2i[3], ttm3i[3];
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal dixdk[3], dixuk[3], dkxui[3], qiqkr[3], qkqir[3], qiukp[
	    3], qkuip[3], qixqk[3], rxqir[3], rxqkr[3], scale3, scale5, 
	    scale7, scale9, dscale[25000], pgamma, pscale[25000], uscale[
	    25000], findmp[3], fridmp[3], dixukp[3], dkxuip[3], dixqkr[3], 
	    scale3i, scale5i, scale7i, uixqkr[3], ukxqir[3], rxqiuk[3], 
	    rxqkui[3], dkxqir[3], rxqikr[3], rxqkir[3], rxqidk[3], rxqkdi[3];
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *), torque2_(doublereal *, doublereal *);
    static logical proceed;
    static doublereal qkrxqir[3], uixqkrp[3], ukxqirp[3], rxqiukp[3], rxqkuip[
	    3];


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define des_ref(a_1,a_2) deriv_1.des[(a_2)*3 + a_1 - 4]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define ip12_ref(a_1,a_2) polgrp_1.ip12[(a_2)*50 + a_1 - 51]
#define ip13_ref(a_1,a_2) polgrp_1.ip13[(a_2)*50 + a_1 - 51]
#define ip14_ref(a_1,a_2) polgrp_1.ip14[(a_2)*50 + a_1 - 51]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define uinp_ref(a_1,a_2) polar_1.uinp[(a_2)*3 + a_1 - 4]
#define trqi_ref(a_1,a_2) trqi[(a_2)*3 + a_1 - 4]
#define ftm1i_ref(a_1,a_2) ftm1i[(a_2)*3 + a_1 - 4]
#define ttm1i_ref(a_1,a_2) ttm1i[(a_2)*3 + a_1 - 4]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]
#define uinds_ref(a_1,a_2) polar_1.uinds[(a_2)*3 + a_1 - 4]
#define uinps_ref(a_1,a_2) polar_1.uinps[(a_2)*3 + a_1 - 4]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  inter.i  --  sum of intermolecular energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     einter   total intermolecular potential energy */




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
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mplpot.i  --  specifics of atomic multipole functions  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     m2scale   factor by which 1-2 multipole interactions are scaled */
/*     m3scale   factor by which 1-3 multipole interactions are scaled */
/*     m4scale   factor by which 1-4 multipole interactions are scaled */
/*     m5scale   factor by which 1-5 multipole interactions are scaled */




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

/*     set arrays needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pscale[i__ - 1] = 1.;
	dscale[i__ - 1] = 1.;
	uscale[i__ - 1] = 1.;
    }

/*     zero out local accumulation arrays for derivatives */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	trqi_ref(1, i__) = 0.;
	trqi_ref(2, i__) = 0.;
	trqi_ref(3, i__) = 0.;
    }

/*     initialise temporary force and torque accumulators */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    ftm1i_ref(j, i__) = 0.;
	    ttm1i_ref(j, i__) = 0.;
	}
    }

/*     set scale factors for permanent multipole and induced terms */

    i__1 = mpole_1.npole - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	iz = mpole_1.zaxis[i__ - 1];
	ix = mpole_1.xaxis[i__ - 1];
	pdi = polar_1.pdamp[i__ - 1];
	pti = polar_1.thole[i__ - 1];
	ci = rpole_ref(1, i__);
	di[0] = rpole_ref(2, i__);
	di[1] = rpole_ref(3, i__);
	di[2] = rpole_ref(4, i__);
	qi[0] = rpole_ref(5, i__);
	qi[1] = rpole_ref(6, i__);
	qi[2] = rpole_ref(7, i__);
	qi[3] = rpole_ref(8, i__);
	qi[4] = rpole_ref(9, i__);
	qi[5] = rpole_ref(10, i__);
	qi[6] = rpole_ref(11, i__);
	qi[7] = rpole_ref(12, i__);
	qi[8] = rpole_ref(13, i__);
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
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = polpot_1.d1scale;
	    uscale[ip11_ref(j, ii) - 1] = polpot_1.u1scale;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = polpot_1.d2scale;
	    uscale[ip12_ref(j, ii) - 1] = polpot_1.u2scale;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = polpot_1.d3scale;
	    uscale[ip13_ref(j, ii) - 1] = polpot_1.u3scale;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = polpot_1.d4scale;
	    uscale[ip14_ref(j, ii) - 1] = polpot_1.u4scale;
	}
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
	    if (! proceed) {
		goto L10;
	    }
	    ck = rpole_ref(1, k);
	    dk[0] = rpole_ref(2, k);
	    dk[1] = rpole_ref(3, k);
	    dk[2] = rpole_ref(4, k);
	    qk[0] = rpole_ref(5, k);
	    qk[1] = rpole_ref(6, k);
	    qk[2] = rpole_ref(7, k);
	    qk[3] = rpole_ref(8, k);
	    qk[4] = rpole_ref(9, k);
	    qk[5] = rpole_ref(10, k);
	    qk[6] = rpole_ref(11, k);
	    qk[7] = rpole_ref(12, k);
	    qk[8] = rpole_ref(13, k);
	    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
	    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
	    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
	    image_(&xr, &yr, &zr);
	    r2 = xr * xr + yr * yr + zr * zr;
	    if (r2 <= shunt_1.off2) {
		r__ = sqrt(r2);
		rr1 = 1. / r__;
		rr3 = rr1 / r2;
		rr5 = rr3 * 3. / r2;
		rr7 = rr5 * 5. / r2;
		rr9 = rr7 * 7. / r2;
		scale3 = 1.;
		scale5 = 1.;
		scale7 = 1.;
		scale9 = 1.;
		for (j = 1; j <= 3; ++j) {
		    ddsc3[j - 1] = 0.;
		    ddsc5[j - 1] = 0.;
		    ddsc7[j - 1] = 0.;
		}

/*     apply Thole polarization damping to scale factors */

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
			scale7 = 1. - (1. - damp + d__1 * d__1 * .6) * exp(
				damp);
/* Computing 2nd power */
			d__1 = damp;
/* Computing 3rd power */
			d__2 = damp;
			scale9 = 1. - (1. - damp + (d__1 * d__1 * 18. - d__2 *
				 (d__2 * d__2) * 9.) / 35.) * exp(damp);
			ddsc3[0] = damp * -3. * exp(damp) * xr / r2;
			ddsc3[1] = damp * -3. * exp(damp) * yr / r2;
			ddsc3[2] = damp * -3. * exp(damp) * zr / r2;
			ddsc5[0] = -damp * ddsc3[0];
			ddsc5[1] = -damp * ddsc3[1];
			ddsc5[2] = -damp * ddsc3[2];
			ddsc7[0] = (-.2 - damp * .6) * ddsc5[0];
			ddsc7[1] = (-.2 - damp * .6) * ddsc5[1];
			ddsc7[2] = (-.2 - damp * .6) * ddsc5[2];
		    }
		}
		scale3i = scale3 * uscale[kk - 1];
		scale5i = scale5 * uscale[kk - 1];
		scale7i = scale7 * uscale[kk - 1];
		dsc3 = scale3 * dscale[kk - 1];
		dsc5 = scale5 * dscale[kk - 1];
		dsc7 = scale7 * dscale[kk - 1];
		dsc9 = scale9 * dscale[kk - 1];
		psc3 = scale3 * pscale[kk - 1];
		psc5 = scale5 * pscale[kk - 1];
		psc7 = scale7 * pscale[kk - 1];
		psc9 = scale9 * pscale[kk - 1];

/*     construct auxiliary vectors for permanent terms */

		dixdk[0] = di[1] * dk[2] - di[2] * dk[1];
		dixdk[1] = di[2] * dk[0] - di[0] * dk[2];
		dixdk[2] = di[0] * dk[1] - di[1] * dk[0];
		dixr[0] = di[1] * zr - di[2] * yr;
		dixr[1] = di[2] * xr - di[0] * zr;
		dixr[2] = di[0] * yr - di[1] * xr;
		dkxr[0] = dk[1] * zr - dk[2] * yr;
		dkxr[1] = dk[2] * xr - dk[0] * zr;
		dkxr[2] = dk[0] * yr - dk[1] * xr;
		qir[0] = qi[0] * xr + qi[3] * yr + qi[6] * zr;
		qir[1] = qi[1] * xr + qi[4] * yr + qi[7] * zr;
		qir[2] = qi[2] * xr + qi[5] * yr + qi[8] * zr;
		qkr[0] = qk[0] * xr + qk[3] * yr + qk[6] * zr;
		qkr[1] = qk[1] * xr + qk[4] * yr + qk[7] * zr;
		qkr[2] = qk[2] * xr + qk[5] * yr + qk[8] * zr;
		qiqkr[0] = qi[0] * qkr[0] + qi[3] * qkr[1] + qi[6] * qkr[2];
		qiqkr[1] = qi[1] * qkr[0] + qi[4] * qkr[1] + qi[7] * qkr[2];
		qiqkr[2] = qi[2] * qkr[0] + qi[5] * qkr[1] + qi[8] * qkr[2];
		qkqir[0] = qk[0] * qir[0] + qk[3] * qir[1] + qk[6] * qir[2];
		qkqir[1] = qk[1] * qir[0] + qk[4] * qir[1] + qk[7] * qir[2];
		qkqir[2] = qk[2] * qir[0] + qk[5] * qir[1] + qk[8] * qir[2];
		qixqk[0] = qi[1] * qk[2] + qi[4] * qk[5] + qi[7] * qk[8] - qi[
			2] * qk[1] - qi[5] * qk[4] - qi[8] * qk[7];
		qixqk[1] = qi[2] * qk[0] + qi[5] * qk[3] + qi[8] * qk[6] - qi[
			0] * qk[2] - qi[3] * qk[5] - qi[6] * qk[8];
		qixqk[2] = qi[0] * qk[1] + qi[3] * qk[4] + qi[6] * qk[7] - qi[
			1] * qk[0] - qi[4] * qk[3] - qi[7] * qk[6];
		rxqir[0] = yr * qir[2] - zr * qir[1];
		rxqir[1] = zr * qir[0] - xr * qir[2];
		rxqir[2] = xr * qir[1] - yr * qir[0];
		rxqkr[0] = yr * qkr[2] - zr * qkr[1];
		rxqkr[1] = zr * qkr[0] - xr * qkr[2];
		rxqkr[2] = xr * qkr[1] - yr * qkr[0];
		rxqikr[0] = yr * qiqkr[2] - zr * qiqkr[1];
		rxqikr[1] = zr * qiqkr[0] - xr * qiqkr[2];
		rxqikr[2] = xr * qiqkr[1] - yr * qiqkr[0];
		rxqkir[0] = yr * qkqir[2] - zr * qkqir[1];
		rxqkir[1] = zr * qkqir[0] - xr * qkqir[2];
		rxqkir[2] = xr * qkqir[1] - yr * qkqir[0];
		qkrxqir[0] = qkr[1] * qir[2] - qkr[2] * qir[1];
		qkrxqir[1] = qkr[2] * qir[0] - qkr[0] * qir[2];
		qkrxqir[2] = qkr[0] * qir[1] - qkr[1] * qir[0];
		qidk[0] = qi[0] * dk[0] + qi[3] * dk[1] + qi[6] * dk[2];
		qidk[1] = qi[1] * dk[0] + qi[4] * dk[1] + qi[7] * dk[2];
		qidk[2] = qi[2] * dk[0] + qi[5] * dk[1] + qi[8] * dk[2];
		qkdi[0] = qk[0] * di[0] + qk[3] * di[1] + qk[6] * di[2];
		qkdi[1] = qk[1] * di[0] + qk[4] * di[1] + qk[7] * di[2];
		qkdi[2] = qk[2] * di[0] + qk[5] * di[1] + qk[8] * di[2];
		dixqkr[0] = di[1] * qkr[2] - di[2] * qkr[1];
		dixqkr[1] = di[2] * qkr[0] - di[0] * qkr[2];
		dixqkr[2] = di[0] * qkr[1] - di[1] * qkr[0];
		dkxqir[0] = dk[1] * qir[2] - dk[2] * qir[1];
		dkxqir[1] = dk[2] * qir[0] - dk[0] * qir[2];
		dkxqir[2] = dk[0] * qir[1] - dk[1] * qir[0];
		rxqidk[0] = yr * qidk[2] - zr * qidk[1];
		rxqidk[1] = zr * qidk[0] - xr * qidk[2];
		rxqidk[2] = xr * qidk[1] - yr * qidk[0];
		rxqkdi[0] = yr * qkdi[2] - zr * qkdi[1];
		rxqkdi[1] = zr * qkdi[0] - xr * qkdi[2];
		rxqkdi[2] = xr * qkdi[1] - yr * qkdi[0];

/*     get intermediate variables for permanent energy terms */

		sc[2] = di[0] * xr + di[1] * yr + di[2] * zr;
		sc[3] = dk[0] * xr + dk[1] * yr + dk[2] * zr;
		sc[4] = qir[0] * xr + qir[1] * yr + qir[2] * zr;
		sc[5] = qkr[0] * xr + qkr[1] * yr + qkr[2] * zr;

/*     construct auxiliary vectors for induced terms */

		dixuk[0] = di[1] * uinds_ref(3, k) - di[2] * uinds_ref(2, k);
		dixuk[1] = di[2] * uinds_ref(1, k) - di[0] * uinds_ref(3, k);
		dixuk[2] = di[0] * uinds_ref(2, k) - di[1] * uinds_ref(1, k);
		dkxui[0] = dk[1] * uinds_ref(3, i__) - dk[2] * uinds_ref(2, 
			i__);
		dkxui[1] = dk[2] * uinds_ref(1, i__) - dk[0] * uinds_ref(3, 
			i__);
		dkxui[2] = dk[0] * uinds_ref(2, i__) - dk[1] * uinds_ref(1, 
			i__);
		dixukp[0] = di[1] * uinps_ref(3, k) - di[2] * uinps_ref(2, k);
		dixukp[1] = di[2] * uinps_ref(1, k) - di[0] * uinps_ref(3, k);
		dixukp[2] = di[0] * uinps_ref(2, k) - di[1] * uinps_ref(1, k);
		dkxuip[0] = dk[1] * uinps_ref(3, i__) - dk[2] * uinps_ref(2, 
			i__);
		dkxuip[1] = dk[2] * uinps_ref(1, i__) - dk[0] * uinps_ref(3, 
			i__);
		dkxuip[2] = dk[0] * uinps_ref(2, i__) - dk[1] * uinps_ref(1, 
			i__);
		qiuk[0] = qi[0] * uinds_ref(1, k) + qi[3] * uinds_ref(2, k) + 
			qi[6] * uinds_ref(3, k);
		qiuk[1] = qi[1] * uinds_ref(1, k) + qi[4] * uinds_ref(2, k) + 
			qi[7] * uinds_ref(3, k);
		qiuk[2] = qi[2] * uinds_ref(1, k) + qi[5] * uinds_ref(2, k) + 
			qi[8] * uinds_ref(3, k);
		qkui[0] = qk[0] * uinds_ref(1, i__) + qk[3] * uinds_ref(2, 
			i__) + qk[6] * uinds_ref(3, i__);
		qkui[1] = qk[1] * uinds_ref(1, i__) + qk[4] * uinds_ref(2, 
			i__) + qk[7] * uinds_ref(3, i__);
		qkui[2] = qk[2] * uinds_ref(1, i__) + qk[5] * uinds_ref(2, 
			i__) + qk[8] * uinds_ref(3, i__);
		qiukp[0] = qi[0] * uinps_ref(1, k) + qi[3] * uinps_ref(2, k) 
			+ qi[6] * uinps_ref(3, k);
		qiukp[1] = qi[1] * uinps_ref(1, k) + qi[4] * uinps_ref(2, k) 
			+ qi[7] * uinps_ref(3, k);
		qiukp[2] = qi[2] * uinps_ref(1, k) + qi[5] * uinps_ref(2, k) 
			+ qi[8] * uinps_ref(3, k);
		qkuip[0] = qk[0] * uinps_ref(1, i__) + qk[3] * uinps_ref(2, 
			i__) + qk[6] * uinps_ref(3, i__);
		qkuip[1] = qk[1] * uinps_ref(1, i__) + qk[4] * uinps_ref(2, 
			i__) + qk[7] * uinps_ref(3, i__);
		qkuip[2] = qk[2] * uinps_ref(1, i__) + qk[5] * uinps_ref(2, 
			i__) + qk[8] * uinps_ref(3, i__);
		uixqkr[0] = uinds_ref(2, i__) * qkr[2] - uinds_ref(3, i__) * 
			qkr[1];
		uixqkr[1] = uinds_ref(3, i__) * qkr[0] - uinds_ref(1, i__) * 
			qkr[2];
		uixqkr[2] = uinds_ref(1, i__) * qkr[1] - uinds_ref(2, i__) * 
			qkr[0];
		ukxqir[0] = uinds_ref(2, k) * qir[2] - uinds_ref(3, k) * qir[
			1];
		ukxqir[1] = uinds_ref(3, k) * qir[0] - uinds_ref(1, k) * qir[
			2];
		ukxqir[2] = uinds_ref(1, k) * qir[1] - uinds_ref(2, k) * qir[
			0];
		uixqkrp[0] = uinps_ref(2, i__) * qkr[2] - uinps_ref(3, i__) * 
			qkr[1];
		uixqkrp[1] = uinps_ref(3, i__) * qkr[0] - uinps_ref(1, i__) * 
			qkr[2];
		uixqkrp[2] = uinps_ref(1, i__) * qkr[1] - uinps_ref(2, i__) * 
			qkr[0];
		ukxqirp[0] = uinps_ref(2, k) * qir[2] - uinps_ref(3, k) * qir[
			1];
		ukxqirp[1] = uinps_ref(3, k) * qir[0] - uinps_ref(1, k) * qir[
			2];
		ukxqirp[2] = uinps_ref(1, k) * qir[1] - uinps_ref(2, k) * qir[
			0];
		rxqiuk[0] = yr * qiuk[2] - zr * qiuk[1];
		rxqiuk[1] = zr * qiuk[0] - xr * qiuk[2];
		rxqiuk[2] = xr * qiuk[1] - yr * qiuk[0];
		rxqkui[0] = yr * qkui[2] - zr * qkui[1];
		rxqkui[1] = zr * qkui[0] - xr * qkui[2];
		rxqkui[2] = xr * qkui[1] - yr * qkui[0];
		rxqiukp[0] = yr * qiukp[2] - zr * qiukp[1];
		rxqiukp[1] = zr * qiukp[0] - xr * qiukp[2];
		rxqiukp[2] = xr * qiukp[1] - yr * qiukp[0];
		rxqkuip[0] = yr * qkuip[2] - zr * qkuip[1];
		rxqkuip[1] = zr * qkuip[0] - xr * qkuip[2];
		rxqkuip[2] = xr * qkuip[1] - yr * qkuip[0];

/*     get intermediate variables for induction energy terms */

		sci[0] = uinds_ref(1, i__) * dk[0] + uinds_ref(2, i__) * dk[1]
			 + uinds_ref(3, i__) * dk[2] + di[0] * uinds_ref(1, k)
			 + di[1] * uinds_ref(2, k) + di[2] * uinds_ref(3, k);
		sci[1] = uinds_ref(1, i__) * uinds_ref(1, k) + uinds_ref(2, 
			i__) * uinds_ref(2, k) + uinds_ref(3, i__) * 
			uinds_ref(3, k);
		sci[2] = uinds_ref(1, i__) * xr + uinds_ref(2, i__) * yr + 
			uinds_ref(3, i__) * zr;
		sci[3] = uinds_ref(1, k) * xr + uinds_ref(2, k) * yr + 
			uinds_ref(3, k) * zr;
		sci[6] = qir[0] * uinds_ref(1, k) + qir[1] * uinds_ref(2, k) 
			+ qir[2] * uinds_ref(3, k);
		sci[7] = qkr[0] * uinds_ref(1, i__) + qkr[1] * uinds_ref(2, 
			i__) + qkr[2] * uinds_ref(3, i__);
		scip[0] = uinps_ref(1, i__) * dk[0] + uinps_ref(2, i__) * dk[
			1] + uinps_ref(3, i__) * dk[2] + di[0] * uinps_ref(1, 
			k) + di[1] * uinps_ref(2, k) + di[2] * uinps_ref(3, k)
			;
		scip[1] = uinds_ref(1, i__) * uinps_ref(1, k) + uinds_ref(2, 
			i__) * uinps_ref(2, k) + uinds_ref(3, i__) * 
			uinps_ref(3, k) + uinps_ref(1, i__) * uinds_ref(1, k) 
			+ uinps_ref(2, i__) * uinds_ref(2, k) + uinps_ref(3, 
			i__) * uinds_ref(3, k);
		scip[2] = uinps_ref(1, i__) * xr + uinps_ref(2, i__) * yr + 
			uinps_ref(3, i__) * zr;
		scip[3] = uinps_ref(1, k) * xr + uinps_ref(2, k) * yr + 
			uinps_ref(3, k) * zr;
		scip[6] = qir[0] * uinps_ref(1, k) + qir[1] * uinps_ref(2, k) 
			+ qir[2] * uinps_ref(3, k);
		scip[7] = qkr[0] * uinps_ref(1, i__) + qkr[1] * uinps_ref(2, 
			i__) + qkr[2] * uinps_ref(3, i__);

/*     calculate the gl functions for potential energy */

		gli[0] = ck * sci[2] - ci * sci[3];
		gli[1] = -sc[2] * sci[3] - sci[2] * sc[3];
		gli[2] = sci[2] * sc[5] - sci[3] * sc[4];
		gli[5] = sci[0];
		gli[6] = (sci[6] - sci[7]) * 2.;
		glip[0] = ck * scip[2] - ci * scip[3];
		glip[1] = -sc[2] * scip[3] - scip[2] * sc[3];
		glip[2] = scip[2] * sc[5] - scip[3] * sc[4];
		glip[5] = scip[0];
		glip[6] = (scip[6] - scip[7]) * 2.;

/*     get the permanent multipole and induced energies */

		ei = (rr3 * (gli[0] + gli[5]) * psc3 + rr5 * (gli[1] + gli[6])
			 * psc5 + rr7 * gli[2] * psc7) * .5;
		ei = f * ei;
		energi_1.es += ei;

/*     intermediate variables for the induced-permanent terms */

		gfi[0] = rr5 * .5 * ((gli[0] + gli[5]) * psc3 + (glip[0] + 
			glip[5]) * dsc3 + scip[1] * scale3i) + rr7 * .5 * ((
			gli[6] + gli[1]) * psc5 + (glip[6] + glip[1]) * dsc5 
			- (sci[2] * scip[3] + scip[2] * sci[3]) * scale5i) + 
			rr9 * .5 * (gli[2] * psc7 + glip[2] * dsc7);
		gfi[1] = -rr3 * ck + rr5 * sc[3] - rr7 * sc[5];
		gfi[2] = rr3 * ci + rr5 * sc[2] + rr7 * sc[4];
		gfi[3] = rr5 * 2.;
		gfi[4] = rr7 * (sci[3] * psc7 + scip[3] * dsc7);
		gfi[5] = -rr7 * (sci[2] * psc7 + scip[2] * dsc7);

/*     get the induced force */

		ftm2i[0] = gfi[0] * xr + (-rr3 * ck * (uinds_ref(1, i__) * 
			psc3 + uinps_ref(1, i__) * dsc3) + rr5 * sc[3] * (
			uinds_ref(1, i__) * psc5 + uinps_ref(1, i__) * dsc5) 
			- rr7 * sc[5] * (uinds_ref(1, i__) * psc7 + uinps_ref(
			1, i__) * dsc7)) * .5 + (rr3 * ci * (uinds_ref(1, k) *
			 psc3 + uinps_ref(1, k) * dsc3) + rr5 * sc[2] * (
			uinds_ref(1, k) * psc5 + uinps_ref(1, k) * dsc5) + 
			rr7 * sc[4] * (uinds_ref(1, k) * psc7 + uinps_ref(1, 
			k) * dsc7)) * .5 + rr5 * scale5i * (sci[3] * 
			uinps_ref(1, i__) + scip[3] * uinds_ref(1, i__) + sci[
			2] * uinps_ref(1, k) + scip[2] * uinds_ref(1, k)) * 
			.5 + (sci[3] * psc5 + scip[3] * dsc5) * .5 * rr5 * di[
			0] + (sci[2] * psc5 + scip[2] * dsc5) * .5 * rr5 * dk[
			0] + gfi[3] * .5 * ((qkui[0] - qiuk[0]) * psc5 + (
			qkuip[0] - qiukp[0]) * dsc5) + gfi[4] * qir[0] + gfi[
			5] * qkr[0];
		ftm2i[1] = gfi[0] * yr + (-rr3 * ck * (uinds_ref(2, i__) * 
			psc3 + uinps_ref(2, i__) * dsc3) + rr5 * sc[3] * (
			uinds_ref(2, i__) * psc5 + uinps_ref(2, i__) * dsc5) 
			- rr7 * sc[5] * (uinds_ref(2, i__) * psc7 + uinps_ref(
			2, i__) * dsc7)) * .5 + (rr3 * ci * (uinds_ref(2, k) *
			 psc3 + uinps_ref(2, k) * dsc3) + rr5 * sc[2] * (
			uinds_ref(2, k) * psc5 + uinps_ref(2, k) * dsc5) + 
			rr7 * sc[4] * (uinds_ref(2, k) * psc7 + uinps_ref(2, 
			k) * dsc7)) * .5 + rr5 * scale5i * (sci[3] * 
			uinps_ref(2, i__) + scip[3] * uinds_ref(2, i__) + sci[
			2] * uinps_ref(2, k) + scip[2] * uinds_ref(2, k)) * 
			.5 + (sci[3] * psc5 + scip[3] * dsc5) * .5 * rr5 * di[
			1] + (sci[2] * psc5 + scip[2] * dsc5) * .5 * rr5 * dk[
			1] + gfi[3] * .5 * ((qkui[1] - qiuk[1]) * psc5 + (
			qkuip[1] - qiukp[1]) * dsc5) + gfi[4] * qir[1] + gfi[
			5] * qkr[1];
		ftm2i[2] = gfi[0] * zr + (-rr3 * ck * (uinds_ref(3, i__) * 
			psc3 + uinps_ref(3, i__) * dsc3) + rr5 * sc[3] * (
			uinds_ref(3, i__) * psc5 + uinps_ref(3, i__) * dsc5) 
			- rr7 * sc[5] * (uinds_ref(3, i__) * psc7 + uinps_ref(
			3, i__) * dsc7)) * .5 + (rr3 * ci * (uinds_ref(3, k) *
			 psc3 + uinps_ref(3, k) * dsc3) + rr5 * sc[2] * (
			uinds_ref(3, k) * psc5 + uinps_ref(3, k) * dsc5) + 
			rr7 * sc[4] * (uinds_ref(3, k) * psc7 + uinps_ref(3, 
			k) * dsc7)) * .5 + rr5 * scale5i * (sci[3] * 
			uinps_ref(3, i__) + scip[3] * uinds_ref(3, i__) + sci[
			2] * uinps_ref(3, k) + scip[2] * uinds_ref(3, k)) * 
			.5 + (sci[3] * psc5 + scip[3] * dsc5) * .5 * rr5 * di[
			2] + (sci[2] * psc5 + scip[2] * dsc5) * .5 * rr5 * dk[
			2] + gfi[3] * .5 * ((qkui[2] - qiuk[2]) * psc5 + (
			qkuip[2] - qiukp[2]) * dsc5) + gfi[4] * qir[2] + gfi[
			5] * qkr[2];

/*     intermediate values needed for partially excluded interactions */

		fridmp[0] = (rr3 * ((gli[0] + gli[5]) * pscale[kk - 1] + (
			glip[0] + glip[5]) * dscale[kk - 1]) * ddsc3[0] + rr5 
			* ((gli[1] + gli[6]) * pscale[kk - 1] + (glip[1] + 
			glip[6]) * dscale[kk - 1]) * ddsc5[0] + rr7 * (gli[2] 
			* pscale[kk - 1] + glip[2] * dscale[kk - 1]) * ddsc7[
			0]) * .5;
		fridmp[1] = (rr3 * ((gli[0] + gli[5]) * pscale[kk - 1] + (
			glip[0] + glip[5]) * dscale[kk - 1]) * ddsc3[1] + rr5 
			* ((gli[1] + gli[6]) * pscale[kk - 1] + (glip[1] + 
			glip[6]) * dscale[kk - 1]) * ddsc5[1] + rr7 * (gli[2] 
			* pscale[kk - 1] + glip[2] * dscale[kk - 1]) * ddsc7[
			1]) * .5;
		fridmp[2] = (rr3 * ((gli[0] + gli[5]) * pscale[kk - 1] + (
			glip[0] + glip[5]) * dscale[kk - 1]) * ddsc3[2] + rr5 
			* ((gli[1] + gli[6]) * pscale[kk - 1] + (glip[1] + 
			glip[6]) * dscale[kk - 1]) * ddsc5[2] + rr7 * (gli[2] 
			* pscale[kk - 1] + glip[2] * dscale[kk - 1]) * ddsc7[
			2]) * .5;

/*     get the induced-induced derivative terms */

		findmp[0] = uscale[kk - 1] * .5 * (scip[1] * rr3 * ddsc3[0] - 
			rr5 * ddsc5[0] * (sci[2] * scip[3] + scip[2] * sci[3])
			);
		findmp[1] = uscale[kk - 1] * .5 * (scip[1] * rr3 * ddsc3[1] - 
			rr5 * ddsc5[1] * (sci[2] * scip[3] + scip[2] * sci[3])
			);
		findmp[2] = uscale[kk - 1] * .5 * (scip[1] * rr3 * ddsc3[2] - 
			rr5 * ddsc5[2] * (sci[2] * scip[3] + scip[2] * sci[3])
			);

/*     handle of scaling for partially excluded interactions */

		ftm2i[0] = ftm2i[0] - fridmp[0] - findmp[0];
		ftm2i[1] = ftm2i[1] - fridmp[1] - findmp[1];
		ftm2i[2] = ftm2i[2] - fridmp[2] - findmp[2];

/*     correction to convert mutual to direct polarization force */

		if (s_cmp(polpot_1.poltyp, "DIRECT", (ftnlen)6, (ftnlen)6) == 
			0) {
		    gfd = (rr5 * scip[1] * scale3i - rr7 * (scip[2] * sci[3] 
			    + sci[2] * scip[3]) * scale5i) * .5;
		    fdir[0] = gfd * xr + rr5 * .5 * scale5i * (sci[3] * 
			    uinps_ref(1, i__) + scip[3] * uinds_ref(1, i__) + 
			    sci[2] * uinps_ref(1, k) + scip[2] * uinds_ref(1, 
			    k));
		    fdir[1] = gfd * yr + rr5 * .5 * scale5i * (sci[3] * 
			    uinps_ref(2, i__) + scip[3] * uinds_ref(2, i__) + 
			    sci[2] * uinps_ref(2, k) + scip[2] * uinds_ref(2, 
			    k));
		    fdir[2] = gfd * zr + rr5 * .5 * scale5i * (sci[3] * 
			    uinps_ref(3, i__) + scip[3] * uinds_ref(3, i__) + 
			    sci[2] * uinps_ref(3, k) + scip[2] * uinds_ref(3, 
			    k));
		    ftm2i[0] = ftm2i[0] - fdir[0] + findmp[0];
		    ftm2i[1] = ftm2i[1] - fdir[1] + findmp[1];
		    ftm2i[2] = ftm2i[2] - fdir[2] + findmp[2];
		}

/*     now perform the torque calculation */
/*     intermediate terms for torque between multipoles i and k */

		gti[1] = (sci[3] * psc5 + scip[3] * dsc5) * .5 * rr5;
		gti[2] = (sci[2] * psc5 + scip[2] * dsc5) * .5 * rr5;
		gti[3] = gfi[3];
		gti[4] = gfi[4];
		gti[5] = gfi[5];

/*     calculate the induced torque components */

		ttm2i[0] = -rr3 * (dixuk[0] * psc3 + dixukp[0] * dsc3) * .5 + 
			gti[1] * dixr[0] + gti[3] * ((ukxqir[0] + rxqiuk[0]) *
			 psc5 + (ukxqirp[0] + rxqiukp[0]) * dsc5) * .5 - gti[
			4] * rxqir[0];
		ttm2i[1] = -rr3 * (dixuk[1] * psc3 + dixukp[1] * dsc3) * .5 + 
			gti[1] * dixr[1] + gti[3] * ((ukxqir[1] + rxqiuk[1]) *
			 psc5 + (ukxqirp[1] + rxqiukp[1]) * dsc5) * .5 - gti[
			4] * rxqir[1];
		ttm2i[2] = -rr3 * (dixuk[2] * psc3 + dixukp[2] * dsc3) * .5 + 
			gti[1] * dixr[2] + gti[3] * ((ukxqir[2] + rxqiuk[2]) *
			 psc5 + (ukxqirp[2] + rxqiukp[2]) * dsc5) * .5 - gti[
			4] * rxqir[2];
		ttm3i[0] = -rr3 * (dkxui[0] * psc3 + dkxuip[0] * dsc3) * .5 + 
			gti[2] * dkxr[0] - gti[3] * ((uixqkr[0] + rxqkui[0]) *
			 psc5 + (uixqkrp[0] + rxqkuip[0]) * dsc5) * .5 - gti[
			5] * rxqkr[0];
		ttm3i[1] = -rr3 * (dkxui[1] * psc3 + dkxuip[1] * dsc3) * .5 + 
			gti[2] * dkxr[1] - gti[3] * ((uixqkr[1] + rxqkui[1]) *
			 psc5 + (uixqkrp[1] + rxqkuip[1]) * dsc5) * .5 - gti[
			5] * rxqkr[1];
		ttm3i[2] = -rr3 * (dkxui[2] * psc3 + dkxuip[2] * dsc3) * .5 + 
			gti[2] * dkxr[2] - gti[3] * ((uixqkr[2] + rxqkui[2]) *
			 psc5 + (uixqkrp[2] + rxqkuip[2]) * dsc5) * .5 - gti[
			5] * rxqkr[2];

/*     update force and torque on site k */

		ftm1i_ref(1, k) = ftm1i_ref(1, k) + ftm2i[0];
		ftm1i_ref(2, k) = ftm1i_ref(2, k) + ftm2i[1];
		ftm1i_ref(3, k) = ftm1i_ref(3, k) + ftm2i[2];
		ttm1i_ref(1, k) = ttm1i_ref(1, k) + ttm3i[0];
		ttm1i_ref(2, k) = ttm1i_ref(2, k) + ttm3i[1];
		ttm1i_ref(3, k) = ttm1i_ref(3, k) + ttm3i[2];

/*     update force and torque on site i */

		ftm1i_ref(1, i__) = ftm1i_ref(1, i__) - ftm2i[0];
		ftm1i_ref(2, i__) = ftm1i_ref(2, i__) - ftm2i[1];
		ftm1i_ref(3, i__) = ftm1i_ref(3, i__) - ftm2i[2];
		ttm1i_ref(1, i__) = ttm1i_ref(1, i__) + ttm2i[0];
		ttm1i_ref(2, i__) = ttm1i_ref(2, i__) + ttm2i[1];
		ttm1i_ref(3, i__) = ttm1i_ref(3, i__) + ttm2i[2];

/*     construct auxiliary vectors for induced terms */

		dixuk[0] = di[1] * uind_ref(3, k) - di[2] * uind_ref(2, k);
		dixuk[1] = di[2] * uind_ref(1, k) - di[0] * uind_ref(3, k);
		dixuk[2] = di[0] * uind_ref(2, k) - di[1] * uind_ref(1, k);
		dkxui[0] = dk[1] * uind_ref(3, i__) - dk[2] * uind_ref(2, i__)
			;
		dkxui[1] = dk[2] * uind_ref(1, i__) - dk[0] * uind_ref(3, i__)
			;
		dkxui[2] = dk[0] * uind_ref(2, i__) - dk[1] * uind_ref(1, i__)
			;
		dixukp[0] = di[1] * uinp_ref(3, k) - di[2] * uinp_ref(2, k);
		dixukp[1] = di[2] * uinp_ref(1, k) - di[0] * uinp_ref(3, k);
		dixukp[2] = di[0] * uinp_ref(2, k) - di[1] * uinp_ref(1, k);
		dkxuip[0] = dk[1] * uinp_ref(3, i__) - dk[2] * uinp_ref(2, 
			i__);
		dkxuip[1] = dk[2] * uinp_ref(1, i__) - dk[0] * uinp_ref(3, 
			i__);
		dkxuip[2] = dk[0] * uinp_ref(2, i__) - dk[1] * uinp_ref(1, 
			i__);
		qiuk[0] = qi[0] * uind_ref(1, k) + qi[3] * uind_ref(2, k) + 
			qi[6] * uind_ref(3, k);
		qiuk[1] = qi[1] * uind_ref(1, k) + qi[4] * uind_ref(2, k) + 
			qi[7] * uind_ref(3, k);
		qiuk[2] = qi[2] * uind_ref(1, k) + qi[5] * uind_ref(2, k) + 
			qi[8] * uind_ref(3, k);
		qkui[0] = qk[0] * uind_ref(1, i__) + qk[3] * uind_ref(2, i__) 
			+ qk[6] * uind_ref(3, i__);
		qkui[1] = qk[1] * uind_ref(1, i__) + qk[4] * uind_ref(2, i__) 
			+ qk[7] * uind_ref(3, i__);
		qkui[2] = qk[2] * uind_ref(1, i__) + qk[5] * uind_ref(2, i__) 
			+ qk[8] * uind_ref(3, i__);
		qiukp[0] = qi[0] * uinp_ref(1, k) + qi[3] * uinp_ref(2, k) + 
			qi[6] * uinp_ref(3, k);
		qiukp[1] = qi[1] * uinp_ref(1, k) + qi[4] * uinp_ref(2, k) + 
			qi[7] * uinp_ref(3, k);
		qiukp[2] = qi[2] * uinp_ref(1, k) + qi[5] * uinp_ref(2, k) + 
			qi[8] * uinp_ref(3, k);
		qkuip[0] = qk[0] * uinp_ref(1, i__) + qk[3] * uinp_ref(2, i__)
			 + qk[6] * uinp_ref(3, i__);
		qkuip[1] = qk[1] * uinp_ref(1, i__) + qk[4] * uinp_ref(2, i__)
			 + qk[7] * uinp_ref(3, i__);
		qkuip[2] = qk[2] * uinp_ref(1, i__) + qk[5] * uinp_ref(2, i__)
			 + qk[8] * uinp_ref(3, i__);
		uixqkr[0] = uind_ref(2, i__) * qkr[2] - uind_ref(3, i__) * 
			qkr[1];
		uixqkr[1] = uind_ref(3, i__) * qkr[0] - uind_ref(1, i__) * 
			qkr[2];
		uixqkr[2] = uind_ref(1, i__) * qkr[1] - uind_ref(2, i__) * 
			qkr[0];
		ukxqir[0] = uind_ref(2, k) * qir[2] - uind_ref(3, k) * qir[1];
		ukxqir[1] = uind_ref(3, k) * qir[0] - uind_ref(1, k) * qir[2];
		ukxqir[2] = uind_ref(1, k) * qir[1] - uind_ref(2, k) * qir[0];
		uixqkrp[0] = uinp_ref(2, i__) * qkr[2] - uinp_ref(3, i__) * 
			qkr[1];
		uixqkrp[1] = uinp_ref(3, i__) * qkr[0] - uinp_ref(1, i__) * 
			qkr[2];
		uixqkrp[2] = uinp_ref(1, i__) * qkr[1] - uinp_ref(2, i__) * 
			qkr[0];
		ukxqirp[0] = uinp_ref(2, k) * qir[2] - uinp_ref(3, k) * qir[1]
			;
		ukxqirp[1] = uinp_ref(3, k) * qir[0] - uinp_ref(1, k) * qir[2]
			;
		ukxqirp[2] = uinp_ref(1, k) * qir[1] - uinp_ref(2, k) * qir[0]
			;
		rxqiuk[0] = yr * qiuk[2] - zr * qiuk[1];
		rxqiuk[1] = zr * qiuk[0] - xr * qiuk[2];
		rxqiuk[2] = xr * qiuk[1] - yr * qiuk[0];
		rxqkui[0] = yr * qkui[2] - zr * qkui[1];
		rxqkui[1] = zr * qkui[0] - xr * qkui[2];
		rxqkui[2] = xr * qkui[1] - yr * qkui[0];
		rxqiukp[0] = yr * qiukp[2] - zr * qiukp[1];
		rxqiukp[1] = zr * qiukp[0] - xr * qiukp[2];
		rxqiukp[2] = xr * qiukp[1] - yr * qiukp[0];
		rxqkuip[0] = yr * qkuip[2] - zr * qkuip[1];
		rxqkuip[1] = zr * qkuip[0] - xr * qkuip[2];
		rxqkuip[2] = xr * qkuip[1] - yr * qkuip[0];

/*     get intermediate variables for induction energy terms */

		sci[0] = uind_ref(1, i__) * dk[0] + uind_ref(2, i__) * dk[1] 
			+ uind_ref(3, i__) * dk[2] + di[0] * uind_ref(1, k) + 
			di[1] * uind_ref(2, k) + di[2] * uind_ref(3, k);
		sci[1] = uind_ref(1, i__) * uind_ref(1, k) + uind_ref(2, i__) 
			* uind_ref(2, k) + uind_ref(3, i__) * uind_ref(3, k);
		sci[2] = uind_ref(1, i__) * xr + uind_ref(2, i__) * yr + 
			uind_ref(3, i__) * zr;
		sci[3] = uind_ref(1, k) * xr + uind_ref(2, k) * yr + uind_ref(
			3, k) * zr;
		sci[6] = qir[0] * uind_ref(1, k) + qir[1] * uind_ref(2, k) + 
			qir[2] * uind_ref(3, k);
		sci[7] = qkr[0] * uind_ref(1, i__) + qkr[1] * uind_ref(2, i__)
			 + qkr[2] * uind_ref(3, i__);
		scip[0] = uinp_ref(1, i__) * dk[0] + uinp_ref(2, i__) * dk[1] 
			+ uinp_ref(3, i__) * dk[2] + di[0] * uinp_ref(1, k) + 
			di[1] * uinp_ref(2, k) + di[2] * uinp_ref(3, k);
		scip[1] = uind_ref(1, i__) * uinp_ref(1, k) + uind_ref(2, i__)
			 * uinp_ref(2, k) + uind_ref(3, i__) * uinp_ref(3, k) 
			+ uinp_ref(1, i__) * uind_ref(1, k) + uinp_ref(2, i__)
			 * uind_ref(2, k) + uinp_ref(3, i__) * uind_ref(3, k);
		scip[2] = uinp_ref(1, i__) * xr + uinp_ref(2, i__) * yr + 
			uinp_ref(3, i__) * zr;
		scip[3] = uinp_ref(1, k) * xr + uinp_ref(2, k) * yr + 
			uinp_ref(3, k) * zr;
		scip[6] = qir[0] * uinp_ref(1, k) + qir[1] * uinp_ref(2, k) + 
			qir[2] * uinp_ref(3, k);
		scip[7] = qkr[0] * uinp_ref(1, i__) + qkr[1] * uinp_ref(2, 
			i__) + qkr[2] * uinp_ref(3, i__);

/*     calculate the gl functions for potential energy */

		gli[0] = ck * sci[2] - ci * sci[3];
		gli[1] = -sc[2] * sci[3] - sci[2] * sc[3];
		gli[2] = sci[2] * sc[5] - sci[3] * sc[4];
		gli[5] = sci[0];
		gli[6] = (sci[6] - sci[7]) * 2.;
		glip[0] = ck * scip[2] - ci * scip[3];
		glip[1] = -sc[2] * scip[3] - scip[2] * sc[3];
		glip[2] = scip[2] * sc[5] - scip[3] * sc[4];
		glip[5] = scip[0];
		glip[6] = (scip[6] - scip[7]) * 2.;

/*     get the permanent multipole and induced energies */

		ei = (rr3 * (gli[0] + gli[5]) * psc3 + rr5 * (gli[1] + gli[6])
			 * psc5 + rr7 * gli[2] * psc7) * -.5;
		ei = f * ei;
		energi_1.es += ei;

/*     intermediate variables for the induced-permanent terms */

		gfi[0] = rr5 * .5 * ((gli[0] + gli[5]) * psc3 + (glip[0] + 
			glip[5]) * dsc3 + scip[1] * scale3i) + rr7 * .5 * ((
			gli[6] + gli[1]) * psc5 + (glip[6] + glip[1]) * dsc5 
			- (sci[2] * scip[3] + scip[2] * sci[3]) * scale5i) + 
			rr9 * .5 * (gli[2] * psc7 + glip[2] * dsc7);
		gfi[1] = -rr3 * ck + rr5 * sc[3] - rr7 * sc[5];
		gfi[2] = rr3 * ci + rr5 * sc[2] + rr7 * sc[4];
		gfi[3] = rr5 * 2.;
		gfi[4] = rr7 * (sci[3] * psc7 + scip[3] * dsc7);
		gfi[5] = -rr7 * (sci[2] * psc7 + scip[2] * dsc7);

/*     get the induced force */

		ftm2i[0] = gfi[0] * xr + (-rr3 * ck * (uind_ref(1, i__) * 
			psc3 + uinp_ref(1, i__) * dsc3) + rr5 * sc[3] * (
			uind_ref(1, i__) * psc5 + uinp_ref(1, i__) * dsc5) - 
			rr7 * sc[5] * (uind_ref(1, i__) * psc7 + uinp_ref(1, 
			i__) * dsc7)) * .5 + (rr3 * ci * (uind_ref(1, k) * 
			psc3 + uinp_ref(1, k) * dsc3) + rr5 * sc[2] * (
			uind_ref(1, k) * psc5 + uinp_ref(1, k) * dsc5) + rr7 *
			 sc[4] * (uind_ref(1, k) * psc7 + uinp_ref(1, k) * 
			dsc7)) * .5 + rr5 * scale5i * (sci[3] * uinp_ref(1, 
			i__) + scip[3] * uind_ref(1, i__) + sci[2] * uinp_ref(
			1, k) + scip[2] * uind_ref(1, k)) * .5 + (sci[3] * 
			psc5 + scip[3] * dsc5) * .5 * rr5 * di[0] + (sci[2] * 
			psc5 + scip[2] * dsc5) * .5 * rr5 * dk[0] + gfi[3] * 
			.5 * ((qkui[0] - qiuk[0]) * psc5 + (qkuip[0] - qiukp[
			0]) * dsc5) + gfi[4] * qir[0] + gfi[5] * qkr[0];
		ftm2i[1] = gfi[0] * yr + (-rr3 * ck * (uind_ref(2, i__) * 
			psc3 + uinp_ref(2, i__) * dsc3) + rr5 * sc[3] * (
			uind_ref(2, i__) * psc5 + uinp_ref(2, i__) * dsc5) - 
			rr7 * sc[5] * (uind_ref(2, i__) * psc7 + uinp_ref(2, 
			i__) * dsc7)) * .5 + (rr3 * ci * (uind_ref(2, k) * 
			psc3 + uinp_ref(2, k) * dsc3) + rr5 * sc[2] * (
			uind_ref(2, k) * psc5 + uinp_ref(2, k) * dsc5) + rr7 *
			 sc[4] * (uind_ref(2, k) * psc7 + uinp_ref(2, k) * 
			dsc7)) * .5 + rr5 * scale5i * (sci[3] * uinp_ref(2, 
			i__) + scip[3] * uind_ref(2, i__) + sci[2] * uinp_ref(
			2, k) + scip[2] * uind_ref(2, k)) * .5 + (sci[3] * 
			psc5 + scip[3] * dsc5) * .5 * rr5 * di[1] + (sci[2] * 
			psc5 + scip[2] * dsc5) * .5 * rr5 * dk[1] + gfi[3] * 
			.5 * ((qkui[1] - qiuk[1]) * psc5 + (qkuip[1] - qiukp[
			1]) * dsc5) + gfi[4] * qir[1] + gfi[5] * qkr[1];
		ftm2i[2] = gfi[0] * zr + (-rr3 * ck * (uind_ref(3, i__) * 
			psc3 + uinp_ref(3, i__) * dsc3) + rr5 * sc[3] * (
			uind_ref(3, i__) * psc5 + uinp_ref(3, i__) * dsc5) - 
			rr7 * sc[5] * (uind_ref(3, i__) * psc7 + uinp_ref(3, 
			i__) * dsc7)) * .5 + (rr3 * ci * (uind_ref(3, k) * 
			psc3 + uinp_ref(3, k) * dsc3) + rr5 * sc[2] * (
			uind_ref(3, k) * psc5 + uinp_ref(3, k) * dsc5) + rr7 *
			 sc[4] * (uind_ref(3, k) * psc7 + uinp_ref(3, k) * 
			dsc7)) * .5 + rr5 * scale5i * (sci[3] * uinp_ref(3, 
			i__) + scip[3] * uind_ref(3, i__) + sci[2] * uinp_ref(
			3, k) + scip[2] * uind_ref(3, k)) * .5 + (sci[3] * 
			psc5 + scip[3] * dsc5) * .5 * rr5 * di[2] + (sci[2] * 
			psc5 + scip[2] * dsc5) * .5 * rr5 * dk[2] + gfi[3] * 
			.5 * ((qkui[2] - qiuk[2]) * psc5 + (qkuip[2] - qiukp[
			2]) * dsc5) + gfi[4] * qir[2] + gfi[5] * qkr[2];

/*     intermediate values needed for partially excluded interactions */

		fridmp[0] = (rr3 * ((gli[0] + gli[5]) * pscale[kk - 1] + (
			glip[0] + glip[5]) * dscale[kk - 1]) * ddsc3[0] + rr5 
			* ((gli[1] + gli[6]) * pscale[kk - 1] + (glip[1] + 
			glip[6]) * dscale[kk - 1]) * ddsc5[0] + rr7 * (gli[2] 
			* pscale[kk - 1] + glip[2] * dscale[kk - 1]) * ddsc7[
			0]) * .5;
		fridmp[1] = (rr3 * ((gli[0] + gli[5]) * pscale[kk - 1] + (
			glip[0] + glip[5]) * dscale[kk - 1]) * ddsc3[1] + rr5 
			* ((gli[1] + gli[6]) * pscale[kk - 1] + (glip[1] + 
			glip[6]) * dscale[kk - 1]) * ddsc5[1] + rr7 * (gli[2] 
			* pscale[kk - 1] + glip[2] * dscale[kk - 1]) * ddsc7[
			1]) * .5;
		fridmp[2] = (rr3 * ((gli[0] + gli[5]) * pscale[kk - 1] + (
			glip[0] + glip[5]) * dscale[kk - 1]) * ddsc3[2] + rr5 
			* ((gli[1] + gli[6]) * pscale[kk - 1] + (glip[1] + 
			glip[6]) * dscale[kk - 1]) * ddsc5[2] + rr7 * (gli[2] 
			* pscale[kk - 1] + glip[2] * dscale[kk - 1]) * ddsc7[
			2]) * .5;

/*     get the induced-induced derivative terms */

		findmp[0] = uscale[kk - 1] * .5 * (scip[1] * rr3 * ddsc3[0] - 
			rr5 * ddsc5[0] * (sci[2] * scip[3] + scip[2] * sci[3])
			);
		findmp[1] = uscale[kk - 1] * .5 * (scip[1] * rr3 * ddsc3[1] - 
			rr5 * ddsc5[1] * (sci[2] * scip[3] + scip[2] * sci[3])
			);
		findmp[2] = uscale[kk - 1] * .5 * (scip[1] * rr3 * ddsc3[2] - 
			rr5 * ddsc5[2] * (sci[2] * scip[3] + scip[2] * sci[3])
			);

/*     handle of scaling for partially excluded interactions */

		ftm2i[0] = ftm2i[0] - fridmp[0] - findmp[0];
		ftm2i[1] = ftm2i[1] - fridmp[1] - findmp[1];
		ftm2i[2] = ftm2i[2] - fridmp[2] - findmp[2];

/*     correction to convert mutual to direct polarization force */

		if (s_cmp(polpot_1.poltyp, "DIRECT", (ftnlen)6, (ftnlen)6) == 
			0) {
		    gfd = (rr5 * scip[1] * scale3i - rr7 * (scip[2] * sci[3] 
			    + sci[2] * scip[3]) * scale5i) * .5;
		    fdir[0] = gfd * xr + rr5 * .5 * scale5i * (sci[3] * 
			    uinp_ref(1, i__) + scip[3] * uind_ref(1, i__) + 
			    sci[2] * uinp_ref(1, k) + scip[2] * uind_ref(1, k)
			    );
		    fdir[1] = gfd * yr + rr5 * .5 * scale5i * (sci[3] * 
			    uinp_ref(2, i__) + scip[3] * uind_ref(2, i__) + 
			    sci[2] * uinp_ref(2, k) + scip[2] * uind_ref(2, k)
			    );
		    fdir[2] = gfd * zr + rr5 * .5 * scale5i * (sci[3] * 
			    uinp_ref(3, i__) + scip[3] * uind_ref(3, i__) + 
			    sci[2] * uinp_ref(3, k) + scip[2] * uind_ref(3, k)
			    );
		    ftm2i[0] = ftm2i[0] - fdir[0] + findmp[0];
		    ftm2i[1] = ftm2i[1] - fdir[1] + findmp[1];
		    ftm2i[2] = ftm2i[2] - fdir[2] + findmp[2];
		}

/*     now perform the torque calculation */
/*     intermediate terms for torque between multipoles i and k */

		gti[1] = (sci[3] * psc5 + scip[3] * dsc5) * .5 * rr5;
		gti[2] = (sci[2] * psc5 + scip[2] * dsc5) * .5 * rr5;
		gti[3] = gfi[3];
		gti[4] = gfi[4];
		gti[5] = gfi[5];

/*     calculate the induced torque components */

		ttm2i[0] = -rr3 * (dixuk[0] * psc3 + dixukp[0] * dsc3) * .5 + 
			gti[1] * dixr[0] + gti[3] * ((ukxqir[0] + rxqiuk[0]) *
			 psc5 + (ukxqirp[0] + rxqiukp[0]) * dsc5) * .5 - gti[
			4] * rxqir[0];
		ttm2i[1] = -rr3 * (dixuk[1] * psc3 + dixukp[1] * dsc3) * .5 + 
			gti[1] * dixr[1] + gti[3] * ((ukxqir[1] + rxqiuk[1]) *
			 psc5 + (ukxqirp[1] + rxqiukp[1]) * dsc5) * .5 - gti[
			4] * rxqir[1];
		ttm2i[2] = -rr3 * (dixuk[2] * psc3 + dixukp[2] * dsc3) * .5 + 
			gti[1] * dixr[2] + gti[3] * ((ukxqir[2] + rxqiuk[2]) *
			 psc5 + (ukxqirp[2] + rxqiukp[2]) * dsc5) * .5 - gti[
			4] * rxqir[2];
		ttm3i[0] = -rr3 * (dkxui[0] * psc3 + dkxuip[0] * dsc3) * .5 + 
			gti[2] * dkxr[0] - gti[3] * ((uixqkr[0] + rxqkui[0]) *
			 psc5 + (uixqkrp[0] + rxqkuip[0]) * dsc5) * .5 - gti[
			5] * rxqkr[0];
		ttm3i[1] = -rr3 * (dkxui[1] * psc3 + dkxuip[1] * dsc3) * .5 + 
			gti[2] * dkxr[1] - gti[3] * ((uixqkr[1] + rxqkui[1]) *
			 psc5 + (uixqkrp[1] + rxqkuip[1]) * dsc5) * .5 - gti[
			5] * rxqkr[1];
		ttm3i[2] = -rr3 * (dkxui[2] * psc3 + dkxuip[2] * dsc3) * .5 + 
			gti[2] * dkxr[2] - gti[3] * ((uixqkr[2] + rxqkui[2]) *
			 psc5 + (uixqkrp[2] + rxqkuip[2]) * dsc5) * .5 - gti[
			5] * rxqkr[2];

/*     update force and torque on site k */

		ftm1i_ref(1, k) = ftm1i_ref(1, k) - ftm2i[0];
		ftm1i_ref(2, k) = ftm1i_ref(2, k) - ftm2i[1];
		ftm1i_ref(3, k) = ftm1i_ref(3, k) - ftm2i[2];
		ttm1i_ref(1, k) = ttm1i_ref(1, k) - ttm3i[0];
		ttm1i_ref(2, k) = ttm1i_ref(2, k) - ttm3i[1];
		ttm1i_ref(3, k) = ttm1i_ref(3, k) - ttm3i[2];

/*     update force and torque on site i */

		ftm1i_ref(1, i__) = ftm1i_ref(1, i__) + ftm2i[0];
		ftm1i_ref(2, i__) = ftm1i_ref(2, i__) + ftm2i[1];
		ftm1i_ref(3, i__) = ftm1i_ref(3, i__) + ftm2i[2];
		ttm1i_ref(1, i__) = ttm1i_ref(1, i__) - ttm2i[0];
		ttm1i_ref(2, i__) = ttm1i_ref(2, i__) - ttm2i[1];
		ttm1i_ref(3, i__) = ttm1i_ref(3, i__) - ttm2i[2];
	    }
L10:
	    ;
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
	i__2 = polgrp_1.np11[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip11_ref(j, ii) - 1] = 1.;
	    uscale[ip11_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip12_ref(j, ii) - 1] = 1.;
	    uscale[ip12_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip13_ref(j, ii) - 1] = 1.;
	    uscale[ip13_ref(j, ii) - 1] = 1.;
	}
	i__2 = polgrp_1.np14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    dscale[ip14_ref(j, ii) - 1] = 1.;
	    uscale[ip14_ref(j, ii) - 1] = 1.;
	}
    }

/*     increment the total forces and torques */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
	des_ref(1, ii) = des_ref(1, ii) - f * ftm1i_ref(1, i__);
	des_ref(2, ii) = des_ref(2, ii) - f * ftm1i_ref(2, i__);
	des_ref(3, ii) = des_ref(3, ii) - f * ftm1i_ref(3, i__);
	trqi_ref(1, ii) = trqi_ref(1, ii) + f * ttm1i_ref(1, i__);
	trqi_ref(2, ii) = trqi_ref(2, ii) + f * ttm1i_ref(2, i__);
	trqi_ref(3, ii) = trqi_ref(3, ii) + f * ttm1i_ref(3, i__);
    }
    torque2_(trqi, deriv_1.des);
    return 0;
} /* ediff1_ */

#undef uinps_ref
#undef uinds_ref
#undef rpole_ref
#undef ttm1i_ref
#undef ftm1i_ref
#undef trqi_ref
#undef uinp_ref
#undef uind_ref
#undef ip14_ref
#undef ip13_ref
#undef ip12_ref
#undef ip11_ref
#undef des_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine epb1  --  Poisson-Boltzmann energy and derivs  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "epb1" calculates the continuum solvation energy and derivatives */
/*     via the Poisson-Boltzmann plus nonpolar implicit solvation */


/* Subroutine */ int epb1_(void)
{
    extern /* Subroutine */ int epb1a_(void), ediff1_(void);



/*     compute the energy and gradients via Poisson-Boltzmann */

    epb1a_();

/*     correct energy and derivatives for vacuum to polarized state */

    ediff1_();
    return 0;
} /* epb1_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine epb1a  --  PB solvation energy and derivatives  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "epb1a" calculates the solvation energy and gradients for the */
/*     PB/NP solvation model */


/* Subroutine */ int epb1a_(void)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int pbempole_(void), pbdirectpolforce_(doublereal 
	    *, doublereal *, doublereal *, doublereal *);
    static integer i__, j;
    extern /* Subroutine */ int pbmutualpolforce_(doublereal *, doublereal *, 
	    doublereal *);
    static integer ii;
    extern /* Subroutine */ int apbsinduce_(doublereal *, doublereal *);
    static doublereal sum;
    extern /* Subroutine */ int apbsnlinduce_(doublereal *, doublereal *);
    static doublereal detor[75000]	/* was [3][25000] */, polgrd[75000]	
	    /* was [3][25000] */;
    extern /* Subroutine */ int torque2_(doublereal *, doublereal *);
    static doublereal directf[75000]	/* was [3][25000] */, indpole[75000]	
	    /* was [3][25000] */, directt[75000]	/* was [3][25000] */, 
	    inppole[75000]	/* was [3][25000] */, mutualf[75000]	/* 
	    was [3][25000] */;


#define des_ref(a_1,a_2) deriv_1.des[(a_2)*3 + a_1 - 4]
#define pbep_ref(a_1,a_2) pb_1.pbep[(a_2)*3 + a_1 - 4]
#define pbfp_ref(a_1,a_2) pb_1.pbfp[(a_2)*3 + a_1 - 4]
#define detor_ref(a_1,a_2) detor[(a_2)*3 + a_1 - 4]
#define uinds_ref(a_1,a_2) polar_1.uinds[(a_2)*3 + a_1 - 4]
#define uinps_ref(a_1,a_2) polar_1.uinps[(a_2)*3 + a_1 - 4]
#define polgrd_ref(a_1,a_2) polgrd[(a_2)*3 + a_1 - 4]
#define directf_ref(a_1,a_2) directf[(a_2)*3 + a_1 - 4]
#define pbeuind_ref(a_1,a_2) pb_1.pbeuind[(a_2)*3 + a_1 - 4]
#define indpole_ref(a_1,a_2) indpole[(a_2)*3 + a_1 - 4]
#define directt_ref(a_1,a_2) directt[(a_2)*3 + a_1 - 4]
#define pbeuinp_ref(a_1,a_2) pb_1.pbeuinp[(a_2)*3 + a_1 - 4]
#define inppole_ref(a_1,a_2) inppole[(a_2)*3 + a_1 - 4]
#define mutualf_ref(a_1,a_2) mutualf[(a_2)*3 + a_1 - 4]



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




/*     induced dipole continuum energy via their */
/*     interaction with the permanent multipoles */

    if (potent_1.use_polar__) {
	sum = 0.;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    sum = sum + uinds_ref(1, i__) * pbep_ref(1, ii) + uinds_ref(2, 
		    i__) * pbep_ref(2, ii) + uinds_ref(3, i__) * pbep_ref(3, 
		    ii);
	}
	sum = chgpot_1.electric * -.5 * sum;
	pb_1.pbe += sum;

/*     initialize induced dipole continuum energy gradients */

	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		indpole_ref(j, i__) = 0.;
		inppole_ref(j, i__) = 0.;
		directf_ref(j, i__) = 0.;
		directt_ref(j, i__) = 0.;
		mutualf_ref(j, i__) = 0.;
		polgrd_ref(j, i__) = 0.;
	    }
	}

/*     copy induced electrostatics into atom-based arrays */

	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    for (j = 1; j <= 3; ++j) {
		indpole_ref(j, ii) = uinds_ref(j, i__);
		inppole_ref(j, ii) = uinps_ref(j, i__);
	    }
	}

/*     for direct polarization, the reaction field due to the */
/*     induced dipoles still needs to be computed because */
/*     the mutual portion of "apbsinduce" was not called */

	if (s_cmp(polpot_1.poltyp, "DIRECT", (ftnlen)6, (ftnlen)6) == 0) {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (j = 1; j <= 3; ++j) {
		    pbeuind_ref(j, i__) = 0.;
		    pbeuinp_ref(j, i__) = 0.;
		}
	    }
	    apbsinduce_(indpole, pb_1.pbeuind);
	    apbsnlinduce_(inppole, pb_1.pbeuinp);
	}

/*     compute direct induced dipole continuum solvation energy */
/*     gradients using potentials saved during the SCRF convergence */

	pbdirectpolforce_(indpole, inppole, directf, directt);

/*     convert torques due to induced dipole reaction field acting */
/*     on permanent multipoles into forces on adjacent atoms */

	torque2_(directt, polgrd);
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    polgrd_ref(1, i__) = polgrd_ref(1, i__) - directf_ref(1, i__);
	    polgrd_ref(2, i__) = polgrd_ref(2, i__) - directf_ref(2, i__);
	    polgrd_ref(3, i__) = polgrd_ref(3, i__) - directf_ref(3, i__);
	}

/*     compute mutual induced dipole solvation energy gradients */

	if (s_cmp(polpot_1.poltyp, "MUTUAL", (ftnlen)6, (ftnlen)6) == 0) {
	    pbmutualpolforce_(indpole, inppole, mutualf);
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		polgrd_ref(1, i__) = polgrd_ref(1, i__) - mutualf_ref(1, i__);
		polgrd_ref(2, i__) = polgrd_ref(2, i__) - mutualf_ref(2, i__);
		polgrd_ref(3, i__) = polgrd_ref(3, i__) - mutualf_ref(3, i__);
	    }
	}

/*     add induced dipole continuum solvation energy gradients */
/*     to overall polarization energy gradients */

	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    des_ref(1, i__) = des_ref(1, i__) + polgrd_ref(1, i__);
	    des_ref(2, i__) = des_ref(2, i__) + polgrd_ref(2, i__);
	    des_ref(3, i__) = des_ref(3, i__) + polgrd_ref(3, i__);
	}

/*     if polarization is off, get the permanent reaction field */

    } else {
	pbempole_();
    }

/*     increment solvation energy by Poisson-Boltzmann results */

    energi_1.es += pb_1.pbe;

/*     convert torques on permanent moments due to their own reaction */
/*     field into forces on adjacent atoms */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    detor_ref(j, i__) = 0.;
	}
    }
    torque2_(pb_1.pbtp, detor);

/*     add permanent reaction field forces to the torque results */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	des_ref(1, i__) = des_ref(1, i__) - pbfp_ref(1, i__) + detor_ref(1, 
		i__);
	des_ref(2, i__) = des_ref(2, i__) - pbfp_ref(2, i__) + detor_ref(2, 
		i__);
	des_ref(3, i__) = des_ref(3, i__) - pbfp_ref(3, i__) + detor_ref(3, 
		i__);
    }
    return 0;
} /* epb1a_ */

#undef mutualf_ref
#undef inppole_ref
#undef pbeuinp_ref
#undef directt_ref
#undef indpole_ref
#undef pbeuind_ref
#undef directf_ref
#undef polgrd_ref
#undef uinps_ref
#undef uinds_ref
#undef detor_ref
#undef pbfp_ref
#undef pbep_ref
#undef des_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine enp1  --  cavity/dispersion energy and derivs  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "enp1" calculates the nonpolar continuum solvation energy */
/*     and derivatives as a sum of cavity and dispersion terms */


/* Subroutine */ int enp1_(doublereal *ecav, doublereal *edisp)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal dtapersa;
    static integer i__;
    static doublereal reff, evol, dvol[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int ewca1_(doublereal *);
    static doublereal reff2, reff3, reff4, reff5, dreff, esurf, dsurf[75000]	
	    /* was [3][25000] */, taperv, aesurf[25000];
    extern /* Subroutine */ int switch_(char *, ftnlen), volume_(doublereal *,
	     doublereal *, doublereal *), volume1_(doublereal *, doublereal *,
	     doublereal *);
    static doublereal exclude, tapersa, dtaperv;
    extern /* Subroutine */ int surface1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


#define des_ref(a_1,a_2) deriv_1.des[(a_2)*3 + a_1 - 4]
#define dvol_ref(a_1,a_2) dvol[(a_2)*3 + a_1 - 4]
#define dsurf_ref(a_1,a_2) dsurf[(a_2)*3 + a_1 - 4]



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




/*     zero out the nonpolar solvation energy and first derivatives */

    *ecav = 0.;
    *edisp = 0.;
/*     compute SASA and effective radius needed for cavity term */

    exclude = 1.4;
    surface1_(&esurf, aesurf, dsurf, npolar_1.rcav, solute_1.asolv, &exclude);
    reff = sqrt(esurf / (npolar_1.surften * 3.141592653589793238)) * .5;
    dreff = reff / (esurf * 2.);
    reff2 = reff * reff;
    reff3 = reff2 * reff;
    reff4 = reff3 * reff;
    reff5 = reff4 * reff;

/*     compute solvent excluded volume for needed for small solutes */

    if (reff < npolar_1.spoff) {
	volume_(&evol, npolar_1.rcav, &exclude);
	evol *= npolar_1.solvprs;
	volume1_(npolar_1.rcav, &exclude, dvol);
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dvol_ref(1, i__) = dvol_ref(1, i__) * npolar_1.solvprs;
	    dvol_ref(2, i__) = dvol_ref(2, i__) * npolar_1.solvprs;
	    dvol_ref(3, i__) = dvol_ref(3, i__) * npolar_1.solvprs;
	}
    }

/*     find cavity energy from only the solvent excluded volume */

    if (reff <= npolar_1.spcut) {
	*ecav = evol;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    des_ref(1, i__) = des_ref(1, i__) + dvol_ref(1, i__);
	    des_ref(2, i__) = des_ref(2, i__) + dvol_ref(2, i__);
	    des_ref(3, i__) = des_ref(3, i__) + dvol_ref(3, i__);
	}

/*     find cavity energy from only a tapered volume term */

    } else if (reff > npolar_1.spcut && reff <= npolar_1.stoff) {
	switch_("GKV", (ftnlen)3);
	taperv = shunt_1.c5 * reff5 + shunt_1.c4 * reff4 + shunt_1.c3 * reff3 
		+ shunt_1.c2 * reff2 + shunt_1.c1 * reff + shunt_1.c0;
	dtaperv = (shunt_1.c5 * 5. * reff4 + shunt_1.c4 * 4. * reff3 + 
		shunt_1.c3 * 3. * reff2 + shunt_1.c2 * 2. * reff + shunt_1.c1)
		 * dreff;
	*ecav = evol * taperv;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    des_ref(1, i__) = des_ref(1, i__) + taperv * dvol_ref(1, i__) + 
		    evol * dtaperv * dsurf_ref(1, i__);
	    des_ref(2, i__) = des_ref(2, i__) + taperv * dvol_ref(2, i__) + 
		    evol * dtaperv * dsurf_ref(2, i__);
	    des_ref(3, i__) = des_ref(3, i__) + taperv * dvol_ref(3, i__) + 
		    evol * dtaperv * dsurf_ref(3, i__);
	}

/*     find cavity energy using both volume and SASA terms */

    } else if (reff > npolar_1.stoff && reff <= npolar_1.spoff) {
	switch_("GKV", (ftnlen)3);
	taperv = shunt_1.c5 * reff5 + shunt_1.c4 * reff4 + shunt_1.c3 * reff3 
		+ shunt_1.c2 * reff2 + shunt_1.c1 * reff + shunt_1.c0;
	dtaperv = (shunt_1.c5 * 5. * reff4 + shunt_1.c4 * 4. * reff3 + 
		shunt_1.c3 * 3. * reff2 + shunt_1.c2 * 2. * reff + shunt_1.c1)
		 * dreff;
	switch_("GKSA", (ftnlen)4);
	tapersa = shunt_1.c5 * reff5 + shunt_1.c4 * reff4 + shunt_1.c3 * 
		reff3 + shunt_1.c2 * reff2 + shunt_1.c1 * reff + shunt_1.c0;
	tapersa = 1. - tapersa;
	dtapersa = (shunt_1.c5 * 5. * reff4 + shunt_1.c4 * 4. * reff3 + 
		shunt_1.c3 * 3. * reff2 + shunt_1.c2 * 2. * reff + shunt_1.c1)
		 * dreff;
	dtapersa = -dtapersa;
	*ecav = evol * taperv;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    des_ref(1, i__) = des_ref(1, i__) + taperv * dvol_ref(1, i__) + 
		    evol * dtaperv * dsurf_ref(1, i__);
	    des_ref(2, i__) = des_ref(2, i__) + taperv * dvol_ref(2, i__) + 
		    evol * dtaperv * dsurf_ref(2, i__);
	    des_ref(3, i__) = des_ref(3, i__) + taperv * dvol_ref(3, i__) + 
		    evol * dtaperv * dsurf_ref(3, i__);
	}
	*ecav += tapersa * esurf;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    des_ref(1, i__) = des_ref(1, i__) + (tapersa + esurf * dtapersa) *
		     dsurf_ref(1, i__);
	    des_ref(2, i__) = des_ref(2, i__) + (tapersa + esurf * dtapersa) *
		     dsurf_ref(2, i__);
	    des_ref(3, i__) = des_ref(3, i__) + (tapersa + esurf * dtapersa) *
		     dsurf_ref(3, i__);
	}

/*     find cavity energy from only a tapered SASA term */

    } else if (reff > npolar_1.spoff && reff <= npolar_1.stcut) {
	switch_("GKSA", (ftnlen)4);
	tapersa = shunt_1.c5 * reff5 + shunt_1.c4 * reff4 + shunt_1.c3 * 
		reff3 + shunt_1.c2 * reff2 + shunt_1.c1 * reff + shunt_1.c0;
	tapersa = 1. - tapersa;
	dtapersa = (shunt_1.c5 * 5. * reff4 + shunt_1.c4 * 4. * reff3 + 
		shunt_1.c3 * 3. * reff2 + shunt_1.c2 * 2. * reff + shunt_1.c1)
		 * dreff;
	dtapersa = -dtapersa;
	*ecav = tapersa * esurf;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    des_ref(1, i__) = des_ref(1, i__) + (tapersa + esurf * dtapersa) *
		     dsurf_ref(1, i__);
	    des_ref(2, i__) = des_ref(2, i__) + (tapersa + esurf * dtapersa) *
		     dsurf_ref(2, i__);
	    des_ref(3, i__) = des_ref(3, i__) + (tapersa + esurf * dtapersa) *
		     dsurf_ref(3, i__);
	}

/*     find cavity energy from only a SASA-based term */

    } else {
	*ecav = esurf;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    des_ref(1, i__) = des_ref(1, i__) + dsurf_ref(1, i__);
	    des_ref(2, i__) = des_ref(2, i__) + dsurf_ref(2, i__);
	    des_ref(3, i__) = des_ref(3, i__) + dsurf_ref(3, i__);
	}
    }

/*     find the continuum dispersion solvation energy */

    ewca1_(edisp);
    return 0;
} /* enp1_ */

#undef dsurf_ref
#undef dvol_ref
#undef des_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine ewca1  --  find WCA dispersion energy & derivs  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "ewca1" finds the Weeks-Chandler-Anderson dispersion energy */
/*     and derivatives of a solute */


/* Subroutine */ int ewca1_(doublereal *edisp)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e;
    static integer i__, k;
    static doublereal r__, r2, r3, ah, de, ao, dl, du, ri, rk, sk, xi, yi, zi,
	     xr, yr, zr, sk2, lik, uik, sum, lik2, lik3, lik4, lik5, lik6, 
	    uik2, uik3, uik4, uik5, uik6, lik10, iwca, lik11, lik12, lik13, 
	    uik10, uik11, irep, epsi, term, rmax, uik12, uik13, dedx, dedy, 
	    dedz, shctd, idisp, emixh, rmini, emixo, rmixh, rmixo, rmixh7, 
	    rmixo7, offset;


#define des_ref(a_1,a_2) deriv_1.des[(a_2)*3 + a_1 - 4]



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
/*     nvt        number of distinct vdw types/classes in the system */
/*     ivt        type/class index for each distinct vdw type or class */
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

/*     find the WCA dispersion energy and gradient components */

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
	ri = npolar_1.rdisp[i__ - 1];

/*     remove the contribution due to solvent displaced by solute atoms */

	xi = atoms_1.x[i__ - 1];
	yi = atoms_1.y[i__ - 1];
	zi = atoms_1.z__[i__ - 1];
	sum = 0.;
	i__2 = atoms_1.n;
	for (k = 1; k <= i__2; ++k) {
	    if (i__ != k) {
		xr = xi - atoms_1.x[k - 1];
		yr = yi - atoms_1.y[k - 1];
		zr = zi - atoms_1.z__[k - 1];
		r2 = xr * xr + yr * yr + zr * zr;
		r__ = sqrt(r2);
		r3 = r__ * r2;
		rk = npolar_1.rdisp[k - 1];
/*              sk = rk * shct(k) */
		sk = rk * shctd;
		sk2 = sk * sk;
		if (ri < r__ + sk) {
		    de = 0.;
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
			if (ri > r__ - sk) {
			    dl = -lik2 + r2 * 2. + sk2 * 2.;
			    dl *= lik2;
			} else {
			    dl = -lik3 + lik2 * 4. * r__ - lik * 6. * r2 + 
				    lik * 2. * sk2 + r3 * 4. - r__ * 4. * sk2;
			    dl *= lik;
			}
			if (r__ + sk > rmixo) {
			    du = -uik2 + r2 * 2. + sk2 * 2.;
			    du = -du * uik2;
			} else {
			    du = -uik3 + uik2 * 4. * r__ - uik * 6. * r2 + 
				    uik * 2. * sk2 + r3 * 4. - r__ * 4. * sk2;
			    du = -du * uik;
			}
			iwca = -emixo * term;
			de -= emixo * 3.141592653589793238 * (dl + du) / (r2 *
				 4.);
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
			if (ri > r__ - sk) {
			    dl = -lik2 + r2 * 2. + sk2 * 2.;
			    dl *= lik2;
			} else {
			    dl = -lik3 + lik2 * 4. * r__ - lik * 6. * r2 + 
				    lik * 2. * sk2 + r3 * 4. - r__ * 4. * sk2;
			    dl *= lik;
			}
			if (r__ + sk > rmixh) {
			    du = -uik2 + r2 * 2. + sk2 * 2.;
			    du = -du * uik2;
			} else {
			    du = -uik3 + uik2 * 4. * r__ - uik * 6. * r2 + 
				    uik * 2. * sk2 + r3 * 4. - r__ * 4. * sk2;
			    du = -du * uik;
			}
			iwca = emixh * -2. * term;
			de -= emixh * 2. * 3.141592653589793238 * (dl + du) / 
				(r2 * 4.);
			sum += iwca;
		    }
		    uik = r__ + sk;
		    uik2 = uik * uik;
		    uik3 = uik2 * uik;
		    uik4 = uik3 * uik;
		    uik5 = uik4 * uik;
		    uik6 = uik5 * uik;
		    uik10 = uik5 * uik5;
		    uik11 = uik10 * uik;
		    uik12 = uik11 * uik;
		    uik13 = uik12 * uik;
		    if (uik > rmixo) {
			lik = max(rmax,rmixo);
			lik2 = lik * lik;
			lik3 = lik2 * lik;
			lik4 = lik3 * lik;
			lik5 = lik4 * lik;
			lik6 = lik5 * lik;
			lik10 = lik5 * lik5;
			lik11 = lik10 * lik;
			lik12 = lik11 * lik;
			lik13 = lik12 * lik;
			term = 12.566370614359172 / (r__ * 120. * lik5 * uik5)
				 * (uik * 15. * lik * r__ * (uik4 - lik4) - 
				uik2 * 10. * lik2 * (uik3 - lik3) + (sk2 - r2)
				 * 6. * (uik5 - lik5));
			if (ri > r__ - sk || rmax < rmixo) {
			    dl = lik2 * -5. + r2 * 3. + sk2 * 3.;
			    dl = -dl / lik5;
			} else {
			    dl = lik3 * 5. - lik * 33. * r2 - lik * 3. * sk2 
				    + (lik2 * r__ + r3 - r__ * sk2) * 15.;
			    dl /= lik6;
			}
			du = uik3 * 5. - uik * 33. * r2 - uik * 3. * sk2 + (
				uik2 * r__ + r3 - r__ * sk2) * 15.;
			du = -du / uik6;
			idisp = ao * -2. * term;
			de -= ao * 2. * 3.141592653589793238 * (dl + du) / (
				r2 * 15.);
			term = 12.566370614359172 / (r__ * 2640. * lik12 * 
				uik12) * (uik * 120. * lik * r__ * (uik11 - 
				lik11) - uik2 * 66. * lik2 * (uik10 - lik10) 
				+ (sk2 - r2) * 55. * (uik12 - lik12));
			if (ri > r__ - sk || rmax < rmixo) {
			    dl = lik2 * -6. + r2 * 5. + sk2 * 5.;
			    dl = -dl / lik12;
			} else {
			    dl = lik3 * 6. - lik * 125. * r2 - lik * 5. * sk2 
				    + (lik2 * r__ + r3 - r__ * sk2) * 60.;
			    dl /= lik13;
			}
			du = uik3 * 6. - uik * 125. * r2 - uik * 5. * sk2 + (
				uik2 * r__ + r3 - r__ * sk2) * 60.;
			du = -du / uik13;
			irep = ao * rmixo7 * term;
			de += ao * rmixo7 * 3.141592653589793238 * (dl + du) /
				 (r2 * 60.);
			sum = sum + irep + idisp;
		    }
		    if (uik > rmixh) {
			lik = max(rmax,rmixh);
			lik2 = lik * lik;
			lik3 = lik2 * lik;
			lik4 = lik3 * lik;
			lik5 = lik4 * lik;
			lik6 = lik5 * lik;
			lik10 = lik5 * lik5;
			lik11 = lik10 * lik;
			lik12 = lik11 * lik;
			lik13 = lik12 * lik;
			term = 12.566370614359172 / (r__ * 120. * lik5 * uik5)
				 * (uik * 15. * lik * r__ * (uik4 - lik4) - 
				uik2 * 10. * lik2 * (uik3 - lik3) + (sk2 - r2)
				 * 6. * (uik5 - lik5));
			if (ri > r__ - sk || rmax < rmixh) {
			    dl = lik2 * -5. + r2 * 3. + sk2 * 3.;
			    dl = -dl / lik5;
			} else {
			    dl = lik3 * 5. - lik * 33. * r2 - lik * 3. * sk2 
				    + (lik2 * r__ + r3 - r__ * sk2) * 15.;
			    dl /= lik6;
			}
			du = uik3 * 5. - uik * 33. * r2 - uik * 3. * sk2 + (
				uik2 * r__ + r3 - r__ * sk2) * 15.;
			du = -du / uik6;
			idisp = ah * -4. * term;
			de -= ah * 4. * 3.141592653589793238 * (dl + du) / (
				r2 * 15.);
			term = 12.566370614359172 / (r__ * 2640. * lik12 * 
				uik12) * (uik * 120. * lik * r__ * (uik11 - 
				lik11) - uik2 * 66. * lik2 * (uik10 - lik10) 
				+ (sk2 - r2) * 55. * (uik12 - lik12));
			if (ri > r__ - sk || rmax < rmixh) {
			    dl = lik2 * -6. + r2 * 5. + sk2 * 5.;
			    dl = -dl / lik12;
			} else {
			    dl = lik3 * 6. - lik * 125. * r2 - lik * 5. * sk2 
				    + (lik2 * r__ + r3 - r__ * sk2) * 60.;
			    dl /= lik13;
			}
			du = uik3 * 6. - uik * 125. * r2 - uik * 5. * sk2 + (
				uik2 * r__ + r3 - r__ * sk2) * 60.;
			du = -du / uik13;
			irep = ah * 2. * rmixh7 * term;
			de += ah * rmixh7 * 3.141592653589793238 * (dl + du) /
				 (r2 * 30.);
			sum = sum + irep + idisp;
		    }

/*     increment the individual dispersion gradient components */

		    de = -de / r__ * 1. * .033428;
		    dedx = de * xr;
		    dedy = de * yr;
		    dedz = de * zr;
		    des_ref(1, i__) = des_ref(1, i__) + dedx;
		    des_ref(2, i__) = des_ref(2, i__) + dedy;
		    des_ref(3, i__) = des_ref(3, i__) + dedz;
		    des_ref(1, k) = des_ref(1, k) - dedx;
		    des_ref(2, k) = des_ref(2, k) - dedy;
		    des_ref(3, k) = des_ref(3, k) - dedz;
		}
	    }
	}

/*     increment the overall dispersion energy component */

	e = npolar_1.cdisp[i__ - 1] - sum * .033427999999999999;
	*edisp += e;
    }
    return 0;
} /* ewca1_ */

#undef des_ref


