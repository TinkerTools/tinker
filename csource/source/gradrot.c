/* gradrot.f -- translated by f2c (version 20050501).
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
    doublereal tesum[1000], teb[1000], tea[1000], teba[1000], teub[1000], 
	    teaa[1000], teopb[1000], teopd[1000], teid[1000], teit[1000], tet[
	    1000], tept[1000], tebt[1000], tett[1000], tev[1000], tec[1000], 
	    tecd[1000], ted[1000], tem[1000], tep[1000], ter[1000], tes[1000],
	     telf[1000], teg[1000], tex[1000];
} domega_;

#define domega_1 domega_

struct {
    doublereal dihed[1000];
    integer nomega, iomega[2000]	/* was [2][1000] */, zline[1000];
} omega_;

#define omega_1 omega_

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
    integer nrot, rot[25000];
    logical use_short__;
} rotate_;

#define rotate_1 rotate_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine gradrot  --  energy and torsional derivs  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "gradrot" calls subroutines to calculate the potential */
/*     energy and its torsional first derivatives */


/* Subroutine */ int gradrot_(doublereal *energy, doublereal *derivs)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static doublereal g[75000]	/* was [3][25000] */;
    static integer i__, j, k, base;
    static doublereal norm, xatom, yatom, zatom, xdist, ydist, zdist, xterm, 
	    yterm, zterm;
    static integer partner;
    extern /* Subroutine */ int rotlist_(integer *, integer *);


#define dea_ref(a_1,a_2) deriv_1.dea[(a_2)*3 + a_1 - 4]
#define deb_ref(a_1,a_2) deriv_1.deb[(a_2)*3 + a_1 - 4]
#define dec_ref(a_1,a_2) deriv_1.dec[(a_2)*3 + a_1 - 4]
#define ded_ref(a_1,a_2) deriv_1.ded[(a_2)*3 + a_1 - 4]
#define deg_ref(a_1,a_2) deriv_1.deg[(a_2)*3 + a_1 - 4]
#define dem_ref(a_1,a_2) deriv_1.dem[(a_2)*3 + a_1 - 4]
#define dep_ref(a_1,a_2) deriv_1.dep[(a_2)*3 + a_1 - 4]
#define der_ref(a_1,a_2) deriv_1.der[(a_2)*3 + a_1 - 4]
#define des_ref(a_1,a_2) deriv_1.des[(a_2)*3 + a_1 - 4]
#define det_ref(a_1,a_2) deriv_1.det[(a_2)*3 + a_1 - 4]
#define dev_ref(a_1,a_2) deriv_1.dev[(a_2)*3 + a_1 - 4]
#define dex_ref(a_1,a_2) deriv_1.dex[(a_2)*3 + a_1 - 4]
#define deaa_ref(a_1,a_2) deriv_1.deaa[(a_2)*3 + a_1 - 4]
#define deba_ref(a_1,a_2) deriv_1.deba[(a_2)*3 + a_1 - 4]
#define decd_ref(a_1,a_2) deriv_1.decd[(a_2)*3 + a_1 - 4]
#define deid_ref(a_1,a_2) deriv_1.deid[(a_2)*3 + a_1 - 4]
#define delf_ref(a_1,a_2) deriv_1.delf[(a_2)*3 + a_1 - 4]
#define debt_ref(a_1,a_2) deriv_1.debt[(a_2)*3 + a_1 - 4]
#define deub_ref(a_1,a_2) deriv_1.deub[(a_2)*3 + a_1 - 4]
#define deit_ref(a_1,a_2) deriv_1.deit[(a_2)*3 + a_1 - 4]
#define dept_ref(a_1,a_2) deriv_1.dept[(a_2)*3 + a_1 - 4]
#define dett_ref(a_1,a_2) deriv_1.dett[(a_2)*3 + a_1 - 4]
#define deopb_ref(a_1,a_2) deriv_1.deopb[(a_2)*3 + a_1 - 4]
#define deopd_ref(a_1,a_2) deriv_1.deopd[(a_2)*3 + a_1 - 4]
#define iomega_ref(a_1,a_2) omega_1.iomega[(a_2)*2 + a_1 - 3]



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

/*     ######################################################### */
/*     ##                                                     ## */
/*     ##  domega.i  --  derivative components over torsions  ## */
/*     ##                                                     ## */
/*     ######################################################### */


/*     tesum   total energy derivatives over torsions */
/*     teb     bond stretch derivatives over torsions */
/*     tea     angle bend derivatives over torsions */
/*     teba    stretch-bend derivatives over torsions */
/*     teub    Urey-Bradley derivatives over torsions */
/*     teaa    angle-angle derivatives over torsions */
/*     teopb   out-of-plane bend derivatives over torsions */
/*     teopd   out-of-plane distance derivatives over torsions */
/*     teid    improper dihedral derivatives over torsions */
/*     teit    improper torsion derivatives over torsions */
/*     tet     torsional derivatives over torsions */
/*     tept    pi-orbital torsion derivatives over torsions */
/*     tebt    stretch-torsion derivatives over torsions */
/*     tett    torsion-torsion derivatives over torsions */
/*     tev     van der Waals derivatives over torsions */
/*     tec     charge-charge derivatives over torsions */
/*     tecd    charge-dipole derivatives over torsions */
/*     ted     dipole-dipole derivatives over torsions */
/*     tem     atomic multipole derivatives over torsions */
/*     tep     polarization derivatives over torsions */
/*     ter     reaction field derivatives over torsions */
/*     tes     solvation derivatives over torsions */
/*     telf    metal ligand field derivatives over torsions */
/*     teg     geometric restraint derivatives over torsions */
/*     tex     extra energy term derivatives over torsions */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  rotate.i  --  molecule partitions for rotation of a bond  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     nrot        total number of atoms moving when bond rotates */
/*     rot         atom numbers of atoms moving when bond rotates */
/*     use_short   logical flag governing use of shortest atom list */




/*     zero out individual components of torsional derivatives */

    /* Parameter adjustments */
    --derivs;

    /* Function Body */
    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	domega_1.teb[i__ - 1] = 0.;
	domega_1.tea[i__ - 1] = 0.;
	domega_1.teba[i__ - 1] = 0.;
	domega_1.teub[i__ - 1] = 0.;
	domega_1.teaa[i__ - 1] = 0.;
	domega_1.teopb[i__ - 1] = 0.;
	domega_1.teopd[i__ - 1] = 0.;
	domega_1.teid[i__ - 1] = 0.;
	domega_1.teit[i__ - 1] = 0.;
	domega_1.tet[i__ - 1] = 0.;
	domega_1.tept[i__ - 1] = 0.;
	domega_1.tebt[i__ - 1] = 0.;
	domega_1.tett[i__ - 1] = 0.;
	domega_1.tev[i__ - 1] = 0.;
	domega_1.tec[i__ - 1] = 0.;
	domega_1.tecd[i__ - 1] = 0.;
	domega_1.ted[i__ - 1] = 0.;
	domega_1.tem[i__ - 1] = 0.;
	domega_1.tep[i__ - 1] = 0.;
	domega_1.ter[i__ - 1] = 0.;
	domega_1.tes[i__ - 1] = 0.;
	domega_1.telf[i__ - 1] = 0.;
	domega_1.teg[i__ - 1] = 0.;
	domega_1.tex[i__ - 1] = 0.;
    }

/*     calculate the energy and Cartesian first derivatives */

    gradient_(energy, g);

/*     transform Cartesian derivatives to torsional space */

    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	base = iomega_ref(1, i__);
	partner = iomega_ref(2, i__);
	rotlist_(&base, &partner);
	xdist = atoms_1.x[base - 1] - atoms_1.x[partner - 1];
	ydist = atoms_1.y[base - 1] - atoms_1.y[partner - 1];
	zdist = atoms_1.z__[base - 1] - atoms_1.z__[partner - 1];
/* Computing 2nd power */
	d__1 = xdist;
/* Computing 2nd power */
	d__2 = ydist;
/* Computing 2nd power */
	d__3 = zdist;
	norm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	xdist /= norm;
	ydist /= norm;
	zdist /= norm;
	i__2 = rotate_1.nrot;
	for (j = 1; j <= i__2; ++j) {
	    k = rotate_1.rot[j - 1];
	    xatom = atoms_1.x[k - 1] - atoms_1.x[base - 1];
	    yatom = atoms_1.y[k - 1] - atoms_1.y[base - 1];
	    zatom = atoms_1.z__[k - 1] - atoms_1.z__[base - 1];
	    xterm = ydist * zatom - zdist * yatom;
	    yterm = zdist * xatom - xdist * zatom;
	    zterm = xdist * yatom - ydist * xatom;
	    domega_1.teb[i__ - 1] = domega_1.teb[i__ - 1] + deb_ref(1, k) * 
		    xterm + deb_ref(2, k) * yterm + deb_ref(3, k) * zterm;
	    domega_1.tea[i__ - 1] = domega_1.tea[i__ - 1] + dea_ref(1, k) * 
		    xterm + dea_ref(2, k) * yterm + dea_ref(3, k) * zterm;
	    domega_1.teba[i__ - 1] = domega_1.teba[i__ - 1] + deba_ref(1, k) *
		     xterm + deba_ref(2, k) * yterm + deba_ref(3, k) * zterm;
	    domega_1.teub[i__ - 1] = domega_1.teub[i__ - 1] + deub_ref(1, k) *
		     xterm + deub_ref(2, k) * yterm + deub_ref(3, k) * zterm;
	    domega_1.teaa[i__ - 1] = domega_1.teaa[i__ - 1] + deaa_ref(1, k) *
		     xterm + deaa_ref(2, k) * yterm + deaa_ref(3, k) * zterm;
	    domega_1.teopb[i__ - 1] = domega_1.teopb[i__ - 1] + deopb_ref(1, 
		    k) * xterm + deopb_ref(2, k) * yterm + deopb_ref(3, k) * 
		    zterm;
	    domega_1.teopd[i__ - 1] = domega_1.teopd[i__ - 1] + deopd_ref(1, 
		    k) * xterm + deopd_ref(2, k) * yterm + deopd_ref(3, k) * 
		    zterm;
	    domega_1.teid[i__ - 1] = domega_1.teid[i__ - 1] + deid_ref(1, k) *
		     xterm + deid_ref(2, k) * yterm + deid_ref(3, k) * zterm;
	    domega_1.teit[i__ - 1] = domega_1.teit[i__ - 1] + deit_ref(1, k) *
		     xterm + deit_ref(2, k) * yterm + deit_ref(3, k) * zterm;
	    domega_1.tet[i__ - 1] = domega_1.tet[i__ - 1] + det_ref(1, k) * 
		    xterm + det_ref(2, k) * yterm + det_ref(3, k) * zterm;
	    domega_1.tept[i__ - 1] = domega_1.tept[i__ - 1] + dept_ref(1, k) *
		     xterm + dept_ref(2, k) * yterm + dept_ref(3, k) * zterm;
	    domega_1.tebt[i__ - 1] = domega_1.tebt[i__ - 1] + debt_ref(1, k) *
		     xterm + debt_ref(2, k) * yterm + debt_ref(3, k) * zterm;
	    domega_1.tett[i__ - 1] = domega_1.tett[i__ - 1] + dett_ref(1, k) *
		     xterm + dett_ref(2, k) * yterm + dett_ref(3, k) * zterm;
	    domega_1.tev[i__ - 1] = domega_1.tev[i__ - 1] + dev_ref(1, k) * 
		    xterm + dev_ref(2, k) * yterm + dev_ref(3, k) * zterm;
	    domega_1.tec[i__ - 1] = domega_1.tec[i__ - 1] + dec_ref(1, k) * 
		    xterm + dec_ref(2, k) * yterm + dec_ref(3, k) * zterm;
	    domega_1.tecd[i__ - 1] = domega_1.tecd[i__ - 1] + decd_ref(1, k) *
		     xterm + decd_ref(2, k) * yterm + decd_ref(3, k) * zterm;
	    domega_1.ted[i__ - 1] = domega_1.ted[i__ - 1] + ded_ref(1, k) * 
		    xterm + ded_ref(2, k) * yterm + ded_ref(3, k) * zterm;
	    domega_1.tem[i__ - 1] = domega_1.tem[i__ - 1] + dem_ref(1, k) * 
		    xterm + dem_ref(2, k) * yterm + dem_ref(3, k) * zterm;
	    domega_1.tep[i__ - 1] = domega_1.tep[i__ - 1] + dep_ref(1, k) * 
		    xterm + dep_ref(2, k) * yterm + dep_ref(3, k) * zterm;
	    domega_1.ter[i__ - 1] = domega_1.ter[i__ - 1] + der_ref(1, k) * 
		    xterm + der_ref(2, k) * yterm + der_ref(3, k) * zterm;
	    domega_1.tes[i__ - 1] = domega_1.tes[i__ - 1] + des_ref(1, k) * 
		    xterm + des_ref(2, k) * yterm + des_ref(3, k) * zterm;
	    domega_1.telf[i__ - 1] = domega_1.telf[i__ - 1] + delf_ref(1, k) *
		     xterm + delf_ref(2, k) * yterm + delf_ref(3, k) * zterm;
	    domega_1.teg[i__ - 1] = domega_1.teg[i__ - 1] + deg_ref(1, k) * 
		    xterm + deg_ref(2, k) * yterm + deg_ref(3, k) * zterm;
	    domega_1.tex[i__ - 1] = domega_1.tex[i__ - 1] + dex_ref(1, k) * 
		    xterm + dex_ref(2, k) * yterm + dex_ref(3, k) * zterm;
	}
    }

/*     sum up to give the total torsional first derivatives */

    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	domega_1.tesum[i__ - 1] = domega_1.teb[i__ - 1] + domega_1.tea[i__ - 
		1] + domega_1.teba[i__ - 1] + domega_1.teub[i__ - 1] + 
		domega_1.teaa[i__ - 1] + domega_1.teopb[i__ - 1] + 
		domega_1.teopd[i__ - 1] + domega_1.teid[i__ - 1] + 
		domega_1.teit[i__ - 1] + domega_1.tet[i__ - 1] + 
		domega_1.tept[i__ - 1] + domega_1.tebt[i__ - 1] + 
		domega_1.tett[i__ - 1] + domega_1.tev[i__ - 1] + domega_1.tec[
		i__ - 1] + domega_1.tecd[i__ - 1] + domega_1.ted[i__ - 1] + 
		domega_1.tem[i__ - 1] + domega_1.tep[i__ - 1] + domega_1.ter[
		i__ - 1] + domega_1.tes[i__ - 1] + domega_1.telf[i__ - 1] + 
		domega_1.teg[i__ - 1] + domega_1.tex[i__ - 1];
	derivs[i__] = domega_1.tesum[i__ - 1];
    }
    return 0;
} /* gradrot_ */

#undef iomega_ref
#undef deopd_ref
#undef deopb_ref
#undef dett_ref
#undef dept_ref
#undef deit_ref
#undef deub_ref
#undef debt_ref
#undef delf_ref
#undef deid_ref
#undef decd_ref
#undef deba_ref
#undef deaa_ref
#undef dex_ref
#undef dev_ref
#undef det_ref
#undef des_ref
#undef der_ref
#undef dep_ref
#undef dem_ref
#undef deg_ref
#undef ded_ref
#undef dec_ref
#undef deb_ref
#undef dea_ref


