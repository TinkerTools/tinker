/* erxnfld.f -- translated by f2c (version 20050501).
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
    doublereal electric, dielec, ebuffer, c2scale, c3scale, c4scale, c5scale;
    logical neutnbr, neutcut;
} chgpot_;

#define chgpot_1 chgpot_

struct {
    doublereal esum, eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, 
	    ebt, ett, ev, ec, ecd, ed, em, ep, er, es, elf, eg, ex;
} energi_;

#define energi_1 energi_

struct {
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

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
    doublereal b1[520]	/* was [40][13] */, b2[520]	/* was [40][13] */;
    integer ijk[216]	/* was [6][6][6] */;
} rxnfld_;

#define rxnfld_1 rxnfld_

struct {
    doublereal rfsize, rfbulkd;
    integer rfterms;
} rxnpot_;

#define rxnpot_1 rxnpot_

/* Table of constant values */

static integer c__3 = 3;



/*     ############################################################ */
/*     ##  COPYRIGHT (C) 1996 by Yong Kong & Jay William Ponder  ## */
/*     ##                  All Rights Reserved                   ## */
/*     ############################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine erxnfld  --  reaction field potential energy  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "erxnfld" calculates the macroscopic reaction field energy */
/*     arising from a set of atomic multipoles */

/*     literature reference: */

/*     Y. Kong and J. W. Ponder, "Reaction Field Methods for Off-Center */
/*     Multipoles", Journal of Chemical Physics, 107, 481-492 (1997) */


/* Subroutine */ int erxnfld_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal r2;
    static integer ii, kk, ix, iz, kx, kz;
    static doublereal xr, yr, zr, eik, rpi[13], rpk[13];
    static logical usei, usek;
    extern /* Subroutine */ int erfik_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *), switch_(
	    char *, ftnlen), ijkpts_(void), chkpole_(void), rotpole_(void);


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




/*     zero out the macroscopic reaction field energy */

    energi_1.er = 0.;

/*     set the switching function coefficients */

    switch_("MPOLE", (ftnlen)5);

/*     check the sign of multipole components at chiral sites */

    chkpole_();

/*     rotate the multipole components into the global frame */

    rotpole_();

/*     compute the indices used in reaction field calculations */

    ijkpts_();

/*     calculate the reaction field interaction energy term */

    i__1 = mpole_1.npole;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = mpole_1.ipole[ii - 1];
	iz = mpole_1.zaxis[ii - 1];
	ix = mpole_1.xaxis[ii - 1];
	usei = usage_1.use[i__ - 1] || usage_1.use[iz - 1] || usage_1.use[ix 
		- 1];
	i__2 = mpole_1.polsiz[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    rpi[j - 1] = rpole_ref(j, ii);
	}
	i__2 = mpole_1.npole;
	for (kk = ii; kk <= i__2; ++kk) {
	    k = mpole_1.ipole[kk - 1];
	    kz = mpole_1.zaxis[kk - 1];
	    kx = mpole_1.xaxis[kk - 1];
	    usek = usage_1.use[k - 1] || usage_1.use[kz - 1] || usage_1.use[
		    kx - 1];
	    if (usei || usek) {
		xr = atoms_1.x[k - 1] - atoms_1.x[i__ - 1];
		yr = atoms_1.y[k - 1] - atoms_1.y[i__ - 1];
		zr = atoms_1.z__[k - 1] - atoms_1.z__[i__ - 1];
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    i__3 = mpole_1.polsiz[kk - 1];
		    for (j = 1; j <= i__3; ++j) {
			rpk[j - 1] = rpole_ref(j, kk);
		    }
		    erfik_(&ii, &kk, &i__, &k, rpi, rpk, &eik);
		    energi_1.er += eik;
		}
	    }
	}
    }
    return 0;
} /* erxnfld_ */

#undef rpole_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine erfik  --  reaction field energy of site pair   ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "erfik" compute the reaction field energy due to a single pair */
/*     of atomic multipoles */


/* Subroutine */ int erfik_(integer *ii, integer *kk, integer *i__, integer *
	k, doublereal *rpi, doublereal *rpk, doublereal *eik)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static doublereal d__;
    static integer m;
    static doublereal d2;
    static integer n1, n2, fi, fj, nn;
    static doublereal xi, yi, zi, xk, yk, zk, ri2, rk2, xi2, yi2, zi2, xk2, 
	    yk2, zk2;
    extern doublereal d1d2_(integer *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    static integer fii;
    static doublereal m2t2[13];
    static integer p_e1__, p_e2__, p_s1__, p_s2__;
    static doublereal term, size;
    static integer isiz, ksiz;
    static doublereal size2, rklij[169]	/* was [13][13] */, ratio;
    static integer ind1_x__[13], ind2_x__[13], ind1_y__[13], ind2_y__[13], 
	    ind1_z__[13], ind2_z__[13];
    static doublereal factor;
    extern /* Subroutine */ int rfindex_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);


#define b1_ref(a_1,a_2) rxnfld_1.b1[(a_2)*40 + a_1 - 41]
#define b2_ref(a_1,a_2) rxnfld_1.b2[(a_2)*40 + a_1 - 41]
#define rklij_ref(a_1,a_2) rklij[(a_2)*13 + a_1 - 14]



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
/*     ##  rxnfld.i  --  reaction field matrix elements and indices  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     b1    first reaction field matrix element array */
/*     b2    second reaction field matrix element array */
/*     ijk   indices into the reaction field element arrays */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  rxnpot.i  --  specifics of reaction field functional form  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     rfsize    radius of reaction field sphere centered at origin */
/*     rfbulkd   bulk dielectric constant of reaction field continuum */
/*     rfterms   number of terms to use in reaction field summation */




/*     get numbers of the atoms */

    /* Parameter adjustments */
    --rpk;
    --rpi;

    /* Function Body */
    isiz = mpole_1.polsiz[*ii - 1];
    ksiz = mpole_1.polsiz[*kk - 1];
    xi = atoms_1.x[*i__ - 1];
    yi = atoms_1.y[*i__ - 1];
    zi = atoms_1.z__[*i__ - 1];
    xk = atoms_1.x[*k - 1];
    yk = atoms_1.y[*k - 1];
    zk = atoms_1.z__[*k - 1];
    d__ = xi * xk + yi * yk + zi * zk;
    ri2 = xi * xi + yi * yi + zi * zi;
    rk2 = xk * xk + yk * yk + zk * zk;

/*     set highest order of multipoles at each site (M=0, D=1, Q=2) */

    *eik = 0.;
    n1 = 2;
    n2 = 2;
    nn = rxnpot_1.rfterms;
    ratio = rxnpot_1.rfbulkd / chgpot_1.dielec;
    factor = chgpot_1.electric * (1. - ratio);
    if (*i__ == *k) {
	factor *= .5;
    }
    size = 1. / rxnpot_1.rfsize;
    size2 = size * size;

/*     get the values of the indices */

    i__1 = n1 + 1;
    m = (pow_ii(&c__3, &i__1) - 1) / 2;
    rfindex_(&n1, &m, ind1_x__, ind1_y__, ind1_z__, &p_s1__, &p_e1__);
    i__1 = n2 + 1;
    m = (pow_ii(&c__3, &i__1) - 1) / 2;
    rfindex_(&n2, &m, ind2_x__, ind2_y__, ind2_z__, &p_s2__, &p_e2__);

/*     initialize the stored matrix element arrays */

    i__1 = p_e1__;
    for (fi = 1; fi <= i__1; ++fi) {
	i__2 = p_e2__;
	for (fj = 1; fj <= i__2; ++fj) {
	    b1_ref(fi, fj) = 0.;
	    b2_ref(fi, fj) = 0.;
	}
    }

/*     explicit formula for the 0th summation term */

    if (nn >= 0) {
	*eik = size * rpi[1] * rpk[1] / ratio;
	size *= size2;
    }

/*     explicit formula for the 1st summation term */

    if (nn >= 1) {
	b2_ref(1, 1) = d__;
	b2_ref(1, 2) = xi;
	b2_ref(1, 3) = yi;
	b2_ref(1, 4) = zi;
	b2_ref(2, 1) = xk;
	b2_ref(3, 1) = yk;
	b2_ref(4, 1) = zk;
	b2_ref(2, 2) = 1.;
	b2_ref(3, 3) = 1.;
	b2_ref(4, 4) = 1.;
	for (fi = 1; fi <= 4; ++fi) {
	    m2t2[fi - 1] = 0.;
	    for (fj = 1; fj <= 4; ++fj) {
		m2t2[fi - 1] += b2_ref(fi, fj) * rpk[fj];
	    }
	}
	term = 0.;
	for (fi = 1; fi <= 4; ++fi) {
	    term += rpi[fi] * m2t2[fi - 1];
	}
	term = size * 2. * term / (ratio * 2. + 1.);
	*eik += term;
	size *= size2;
    }

/*     explicit formula for the 2nd summation term */

    if (nn >= 2) {
	b2_ref(1, 1) = (d__ * 3. * d__ - ri2 * rk2) * .5;
	b2_ref(1, 2) = xi * 3. * d__ - xk * ri2;
	b2_ref(1, 3) = yi * 3. * d__ - yk * ri2;
	b2_ref(1, 4) = zi * 3. * d__ - zk * ri2;
	b2_ref(1, 5) = xi * 3. * xi - ri2;
	b2_ref(1, 6) = xi * 3. * yi;
	b2_ref(1, 7) = xi * 3. * zi;
	b2_ref(1, 8) = b2_ref(1, 6);
	b2_ref(1, 9) = yi * 3. * yi - ri2;
	b2_ref(1, 10) = yi * 3. * zi;
	b2_ref(1, 11) = b2_ref(1, 7);
	b2_ref(1, 12) = b2_ref(1, 10);
	b2_ref(1, 13) = zi * 3. * zi - ri2;
	b2_ref(2, 1) = xk * 3. * d__ - xi * rk2;
	b2_ref(2, 2) = d__ * 3. + xi * xk;
	b2_ref(2, 3) = xk * 3. * yi - xi * 2. * yk;
	b2_ref(2, 4) = zi * 3. * xk - xi * 2. * zk;
	b2_ref(2, 5) = xi * 4.;
	b2_ref(2, 6) = yi * 3.;
	b2_ref(2, 7) = zi * 3.;
	b2_ref(2, 8) = b2_ref(2, 6);
	b2_ref(2, 9) = xi * -2.;
	b2_ref(2, 11) = b2_ref(2, 7);
	b2_ref(2, 13) = b2_ref(2, 9);
	b2_ref(3, 1) = yk * 3. * d__ - yi * rk2;
	b2_ref(3, 2) = yk * 3. * xi - yi * 2. * xk;
	b2_ref(3, 3) = d__ * 3. + yi * yk;
	b2_ref(3, 4) = yk * 3. * zi - yi * 2. * zk;
	b2_ref(3, 5) = yi * -2.;
	b2_ref(3, 6) = xi * 3.;
	b2_ref(3, 8) = b2_ref(3, 6);
	b2_ref(3, 9) = yi * 4.;
	b2_ref(3, 10) = zi * 3.;
	b2_ref(3, 12) = b2_ref(3, 10);
	b2_ref(3, 13) = b2_ref(3, 5);
	b2_ref(4, 1) = zk * 3. * d__ - zi * rk2;
	b2_ref(4, 2) = zk * 3. * xi - zi * 2. * xk;
	b2_ref(4, 3) = zk * 3. * yi - zi * 2. * yk;
	b2_ref(4, 4) = d__ * 3. + zi * zk;
	b2_ref(4, 5) = zi * -2.;
	b2_ref(4, 7) = xi * 3.;
	b2_ref(4, 9) = b2_ref(4, 5);
	b2_ref(4, 10) = yi * 3.;
	b2_ref(4, 11) = b2_ref(4, 7);
	b2_ref(4, 12) = b2_ref(4, 10);
	b2_ref(4, 13) = zi * 4.;
	b2_ref(5, 1) = xk * 3. * xk - rk2;
	b2_ref(5, 2) = xk * 4.;
	b2_ref(5, 3) = yk * -2.;
	b2_ref(5, 4) = zk * -2.;
	b2_ref(5, 5) = 4.;
	b2_ref(5, 9) = -2.;
	b2_ref(5, 13) = -2.;
	b2_ref(6, 1) = xk * 3. * yk;
	b2_ref(6, 2) = yk * 3.;
	b2_ref(6, 3) = xk * 3.;
	b2_ref(6, 6) = 3.;
	b2_ref(6, 8) = 3.;
	b2_ref(7, 1) = xk * 3. * zk;
	b2_ref(7, 2) = zk * 3.;
	b2_ref(7, 4) = xk * 3.;
	b2_ref(7, 7) = 3.;
	b2_ref(7, 11) = 3.;
	b2_ref(8, 1) = b2_ref(6, 1);
	b2_ref(8, 2) = b2_ref(6, 2);
	b2_ref(8, 3) = b2_ref(6, 3);
	b2_ref(8, 6) = 3.;
	b2_ref(8, 8) = 3.;
	b2_ref(9, 1) = yk * 3. * yk - rk2;
	b2_ref(9, 2) = xk * -2.;
	b2_ref(9, 3) = yk * 4.;
	b2_ref(9, 4) = zk * -2.;
	b2_ref(9, 5) = -2.;
	b2_ref(9, 9) = 4.;
	b2_ref(9, 13) = -2.;
	b2_ref(10, 1) = yk * 3. * zk;
	b2_ref(10, 3) = zk * 3.;
	b2_ref(10, 4) = yk * 3.;
	b2_ref(10, 10) = 3.;
	b2_ref(10, 12) = 3.;
	b2_ref(11, 1) = b2_ref(7, 1);
	b2_ref(11, 2) = b2_ref(7, 2);
	b2_ref(11, 4) = b2_ref(7, 4);
	b2_ref(11, 7) = 3.;
	b2_ref(11, 11) = 3.;
	b2_ref(12, 1) = b2_ref(10, 1);
	b2_ref(12, 3) = b2_ref(10, 3);
	b2_ref(12, 4) = b2_ref(10, 4);
	b2_ref(12, 10) = 3.;
	b2_ref(12, 12) = 3.;
	b2_ref(13, 1) = zk * 3. * zk - rk2;
	b2_ref(13, 2) = xk * -2.;
	b2_ref(13, 3) = yk * -2.;
	b2_ref(13, 4) = zk * 4.;
	b2_ref(13, 5) = -2.;
	b2_ref(13, 9) = -2.;
	b2_ref(13, 13) = 4.;
	i__1 = isiz;
	for (fi = 1; fi <= i__1; ++fi) {
	    m2t2[fi - 1] = 0.;
	    i__2 = ksiz;
	    for (fj = 1; fj <= i__2; ++fj) {
		m2t2[fi - 1] += b2_ref(fi, fj) * rpk[fj];
	    }
	}
	term = 0.;
	i__1 = isiz;
	for (fi = 1; fi <= i__1; ++fi) {
	    term += rpi[fi] * m2t2[fi - 1];
	}
	term = size * 3. * term / (ratio * 3. + 2.);
	*eik += term;
	size *= size2;
    }

/*     explicit formula for the 3rd summation term */

    if (nn >= 3) {
	d2 = d__ * d__;
	xi2 = xi * xi;
	yi2 = yi * yi;
	zi2 = zi * zi;
	xk2 = xk * xk;
	yk2 = yk * yk;
	zk2 = zk * zk;
	b1_ref(1, 1) = d__ * (d2 * 2.5 - ri2 * 1.5 * rk2);
	b1_ref(1, 2) = d2 * 7.5 * xi - xk * 3. * ri2 * d__ - xi * 1.5 * ri2 * 
		rk2;
	b1_ref(1, 3) = d2 * 7.5 * yi - yk * 3. * ri2 * d__ - yi * 1.5 * ri2 * 
		rk2;
	b1_ref(1, 4) = d2 * 7.5 * zi - zk * 3. * ri2 * d__ - zi * 1.5 * ri2 * 
		rk2;
	b1_ref(1, 5) = d__ * 15. * xi2 - ri2 * 3. * (d__ + xi * 2. * xk);
	b1_ref(1, 6) = xi * 15. * yi * d__ - ri2 * 3. * (xi * yk + xk * yi);
	b1_ref(1, 7) = xi * 15. * zi * d__ - ri2 * 3. * (xi * zk + xk * zi);
	b1_ref(1, 8) = b1_ref(1, 6);
	b1_ref(1, 9) = d__ * 15. * yi2 - ri2 * 3. * (d__ + yi * 2. * yk);
	b1_ref(1, 10) = yi * 15. * zi * d__ - ri2 * 3. * (yi * zk + yk * zi);
	b1_ref(1, 11) = b1_ref(1, 7);
	b1_ref(1, 12) = b1_ref(1, 10);
	b1_ref(1, 13) = d__ * 15. * zi2 - ri2 * 3. * (d__ + zi * 2. * zk);
	b1_ref(2, 1) = d2 * 7.5 * xk - xi * 3. * rk2 * d__ - xk * 1.5 * ri2 * 
		rk2;
	b1_ref(2, 2) = d2 * 7.5 + xi * 9. * xk * d__ - xi2 * 3. * rk2 - xk2 * 
		3. * ri2 - ri2 * 1.5 * rk2;
	b1_ref(2, 3) = ((xk * 5. * yi - xi * 2. * yk) * d__ - xi * yi * rk2 - 
		xk * yk * ri2) * 3.;
	b1_ref(2, 4) = ((xk * 5. * zi - xi * 2. * zk) * d__ - xi * zi * rk2 - 
		xk * zk * ri2) * 3.;
	b1_ref(2, 5) = xi * 24. * yi * yk + xi * 24. * zi * zk + xi2 * 18. * 
		xk - xk * 9. * yi2 - xk * 9. * zi2;
	b1_ref(2, 6) = (yi * 8. * xk * xi - xi2 * 3. * yk + yi2 * 4. * yk - 
		yk * zi2 + yi * 5. * zi * zk) * 3.;
	b1_ref(2, 7) = zi * 15. * yi * yk + zi2 * 12. * zk - xi2 * 9. * zk - 
		zk * 3. * yi2 + zi * 24. * xk * xi;
	b1_ref(2, 8) = b1_ref(2, 6);
	b1_ref(2, 9) = xi2 * -9. * xk + xk * 12. * yi2 - xk * 3. * zi2 - xi * 
		18. * yi * yk - xi * 6. * zi * zk;
	b1_ref(2, 10) = zi * 15. * xk * yi - zi * 6. * xi * yk - yi * 6. * xi 
		* zk;
	b1_ref(2, 11) = b1_ref(2, 7);
	b1_ref(2, 12) = b1_ref(2, 10);
	b1_ref(2, 13) = xi * -6. * yi * yk - xi2 * 9. * xk - xk * 3. * yi2 + 
		xk * 12. * zi2 - xi * 18. * zi * zk;
	b1_ref(3, 1) = d2 * 7.5 * yk - yi * 3. * rk2 * d__ - yk * 1.5 * ri2 * 
		rk2;
	b1_ref(3, 2) = ((xi * 5. * yk - xk * 2. * yi) * d__ - xi * yi * rk2 - 
		xk * yk * ri2) * 3.;
	b1_ref(3, 3) = d2 * 7.5 + yi * 9. * yk * d__ - yi2 * 3. * rk2 - yk2 * 
		3. * ri2 - ri2 * 1.5 * rk2;
	b1_ref(3, 4) = ((yk * 5. * zi - yi * 2. * zk) * d__ - yi * zi * rk2 - 
		yk * zk * ri2) * 3.;
	b1_ref(3, 5) = yi2 * -9. * yk - yi * 6. * zi * zk - yi * 18. * xk * 
		xi + xi2 * 12. * yk - yk * 3. * zi2;
	b1_ref(3, 6) = xi2 * 12. * xk + xi * 15. * zi * zk - xk * 9. * yi2 - 
		xk * 3. * zi2 + xi * 24. * yi * yk;
	b1_ref(3, 7) = zi * 15. * xi * yk - yi * 6. * xi * zk - zi * 6. * xk *
		 yi;
	b1_ref(3, 8) = b1_ref(3, 6);
	b1_ref(3, 9) = xi2 * -9. * yk + yi2 * 18. * yk - yk * 9. * zi2 + yi * 
		24. * xk * xi + yi * 24. * zi * zk;
	b1_ref(3, 10) = zi * 24. * yi * yk - xi2 * 3. * zk - zk * 9. * yi2 + 
		zi2 * 12. * zk + zi * 15. * xk * xi;
	b1_ref(3, 11) = b1_ref(3, 7);
	b1_ref(3, 12) = b1_ref(3, 10);
	b1_ref(3, 13) = xi2 * -3. * yk - yi2 * 9. * yk + yk * 12. * zi2 - yi *
		 18. * zi * zk - yi * 6. * xk * xi;
	b1_ref(4, 1) = d2 * 7.5 * zk - zi * 3. * rk2 * d__ - zk * 1.5 * ri2 * 
		rk2;
	b1_ref(4, 2) = ((xi * 5. * zk - xk * 2. * zi) * d__ - xi * zi * rk2 - 
		xk * zk * ri2) * 3.;
	b1_ref(4, 3) = ((yi * 5. * zk - yk * 2. * zi) * d__ - yi * zi * rk2 - 
		yk * zk * ri2) * 3.;
	b1_ref(4, 4) = d2 * 7.5 + zi * 9. * zk * d__ - zi2 * 3. * rk2 - zk2 * 
		3. * ri2 - ri2 * 1.5 * rk2;
	b1_ref(4, 5) = xi2 * 12. * zk - zk * 3. * yi2 - zi2 * 9. * zk - zi * 
		18. * xk * xi - zi * 6. * yi * yk;
	b1_ref(4, 6) = yi * 15. * xi * zk - zi * 6. * xi * yk - zi * 6. * xk *
		 yi;
	b1_ref(4, 7) = xi * 24. * zi * zk + xi2 * 12. * xk - xk * 3. * yi2 - 
		xk * 9. * zi2 + xi * 15. * yi * yk;
	b1_ref(4, 8) = b1_ref(4, 6);
	b1_ref(4, 9) = zi * -6. * xk * xi - zi2 * 9. * zk - xi2 * 3. * zk + 
		zk * 12. * yi2 - zi * 18. * yi * yk;
	b1_ref(4, 10) = yi * 15. * xk * xi + yi2 * 12. * yk - yk * 9. * zi2 + 
		yi * 24. * zi * zk - xi2 * 3. * yk;
	b1_ref(4, 11) = b1_ref(4, 7);
	b1_ref(4, 12) = b1_ref(4, 10);
	b1_ref(4, 13) = zi * 24. * xk * xi + zi2 * 18. * zk - xi2 * 9. * zk - 
		zk * 9. * yi2 + zi * 24. * yi * yk;
	b1_ref(5, 1) = d__ * 15. * xk2 - rk2 * 3. * (d__ + xi * 2. * xk);
	b1_ref(5, 2) = xi * 18. * xk2 + xk * 24. * yi * yk + xk * 24. * zi * 
		zk - xi * 9. * yk2 - xi * 9. * zk2;
	b1_ref(5, 3) = yi * 12. * xk2 - yk2 * 9. * yi - yi * 3. * zk2 - xk * 
		18. * xi * yk - yk * 6. * zi * zk;
	b1_ref(5, 4) = zk2 * -9. * zi - zk * 6. * yi * yk - xk * 18. * xi * 
		zk + zi * 12. * xk2 - zi * 3. * yk2;
	b1_ref(5, 5) = zi * 24. * zk + yi * 24. * yk + xi * 36. * xk;
	b1_ref(5, 6) = xi * -18. * yk + yi * 24. * xk;
	b1_ref(5, 7) = xi * -18. * zk + zi * 24. * xk;
	b1_ref(5, 8) = b1_ref(5, 6);
	b1_ref(5, 9) = zi * -6. * zk - yi * 18. * yk - xi * 18. * xk;
	b1_ref(5, 10) = (yi * zk + zi * yk) * -6.;
	b1_ref(5, 11) = b1_ref(5, 7);
	b1_ref(5, 12) = b1_ref(5, 10);
	b1_ref(5, 13) = yi * -6. * yk - xi * 18. * xk - zi * 18. * zk;
	b1_ref(6, 1) = xk * 15. * yk * d__ - rk2 * 3. * (xi * yk + xk * yi);
	b1_ref(6, 2) = yi * -9. * xk2 + yk2 * 12. * yi - yi * 3. * zk2 + xk * 
		24. * xi * yk + yk * 15. * zi * zk;
	b1_ref(6, 3) = xi * 12. * xk2 + xk * 15. * zi * zk - xi * 9. * yk2 - 
		xi * 3. * zk2 + xk * 24. * yi * yk;
	b1_ref(6, 4) = xk * -6. * yi * zk - yk * 6. * xi * zk + zi * 15. * xk 
		* yk;
	b1_ref(6, 5) = yi * -18. * xk + xi * 24. * yk;
	b1_ref(6, 6) = yi * 24. * yk + xi * 24. * xk + zi * 15. * zk;
	b1_ref(6, 7) = yi * -6. * zk + zi * 15. * yk;
	b1_ref(6, 8) = b1_ref(6, 6);
	b1_ref(6, 9) = xi * -18. * yk + yi * 24. * xk;
	b1_ref(6, 10) = xi * -6. * zk + zi * 15. * xk;
	b1_ref(6, 11) = b1_ref(6, 7);
	b1_ref(6, 12) = b1_ref(6, 10);
	b1_ref(6, 13) = yi * -6. * xk - xi * 6. * yk;
	b1_ref(7, 1) = xk * 15. * zk * d__ - rk2 * 3. * (xi * zk + xk * zi);
	b1_ref(7, 2) = zk * 15. * yi * yk + zk2 * 12. * zi - zi * 9. * xk2 - 
		zi * 3. * yk2 + xk * 24. * xi * zk;
	b1_ref(7, 3) = zi * -6. * xk * yk - yk * 6. * xi * zk + xk * 15. * yi 
		* zk;
	b1_ref(7, 4) = xi * 12. * xk2 - xi * 3. * yk2 - xi * 9. * zk2 + xk * 
		15. * yi * yk + xk * 24. * zi * zk;
	b1_ref(7, 5) = zi * -18. * xk + xi * 24. * zk;
	b1_ref(7, 6) = zi * -6. * yk + yi * 15. * zk;
	b1_ref(7, 7) = xi * 24. * xk + zi * 24. * zk + yi * 15. * yk;
	b1_ref(7, 8) = b1_ref(7, 6);
	b1_ref(7, 9) = zi * -6. * xk - xi * 6. * zk;
	b1_ref(7, 10) = xi * -6. * yk + yi * 15. * xk;
	b1_ref(7, 11) = b1_ref(7, 7);
	b1_ref(7, 12) = b1_ref(7, 10);
	b1_ref(7, 13) = xi * -18. * zk + zi * 24. * xk;
	b1_ref(9, 1) = d__ * 15. * yk2 - rk2 * 3. * (d__ + yi * 2. * yk);
	b1_ref(9, 2) = xi * -9. * xk2 + xi * 12. * yk2 - xi * 3. * zk2 - xk * 
		18. * yi * yk - xk * 6. * zi * zk;
	b1_ref(9, 3) = yi * -9. * xk2 + yk2 * 18. * yi - yi * 9. * zk2 + yk * 
		24. * zi * zk + xk * 24. * xi * yk;
	b1_ref(9, 4) = zi * 12. * yk2 - zk * 18. * yi * yk - zi * 3. * xk2 - 
		zk2 * 9. * zi - xk * 6. * xi * zk;
	b1_ref(9, 5) = xi * -18. * xk - zi * 6. * zk - yi * 18. * yk;
	b1_ref(9, 6) = yi * -18. * xk + xi * 24. * yk;
	b1_ref(9, 7) = zi * -6. * xk - xi * 6. * zk;
	b1_ref(9, 8) = b1_ref(9, 6);
	b1_ref(9, 9) = xi * 24. * xk + zi * 24. * zk + yi * 36. * yk;
	b1_ref(9, 10) = yi * -18. * zk + zi * 24. * yk;
	b1_ref(9, 11) = b1_ref(9, 7);
	b1_ref(9, 12) = b1_ref(9, 10);
	b1_ref(9, 13) = yi * -18. * yk - xi * 6. * xk - zi * 18. * zk;
	b1_ref(10, 1) = yk * 15. * zk * d__ - rk2 * 3. * (yi * zk + yk * zi);
	b1_ref(10, 2) = zi * -6. * xk * yk - xk * 6. * yi * zk + yk * 15. * 
		xi * zk;
	b1_ref(10, 3) = zk2 * 12. * zi + xk * 15. * xi * zk - zi * 3. * xk2 - 
		zi * 9. * yk2 + zk * 24. * yi * yk;
	b1_ref(10, 4) = xk * 15. * xi * yk + yk2 * 12. * yi - yi * 3. * xk2 - 
		yi * 9. * zk2 + yk * 24. * zi * zk;
	b1_ref(10, 5) = yi * -6. * zk - zi * 6. * yk;
	b1_ref(10, 6) = zi * -6. * xk + xi * 15. * zk;
	b1_ref(10, 7) = yi * -6. * xk + xi * 15. * yk;
	b1_ref(10, 8) = b1_ref(10, 6);
	b1_ref(10, 9) = yi * 24. * zk - zi * 18. * yk;
	b1_ref(10, 10) = xi * 15. * xk + zi * 24. * zk + yi * 24. * yk;
	b1_ref(10, 11) = b1_ref(10, 7);
	b1_ref(10, 12) = b1_ref(10, 10);
	b1_ref(10, 13) = yi * -18. * zk + zi * 24. * yk;
	b1_ref(13, 1) = d__ * 15. * zk2 - rk2 * 3. * (d__ + zi * 2. * zk);
	b1_ref(13, 2) = xi * 12. * zk2 - xk * 18. * zi * zk - xi * 9. * xk2 - 
		xi * 3. * yk2 - xk * 6. * yi * yk;
	b1_ref(13, 3) = yi * 12. * zk2 - yi * 3. * xk2 - yk2 * 9. * yi - yk * 
		18. * zi * zk - xk * 6. * xi * yk;
	b1_ref(13, 4) = zi * -9. * xk2 - zi * 9. * yk2 + zk2 * 18. * zi + xk *
		 24. * xi * zk + zk * 24. * yi * yk;
	b1_ref(13, 5) = yi * -6. * yk - zi * 18. * zk - xi * 18. * xk;
	b1_ref(13, 6) = yi * -6. * xk - xi * 6. * yk;
	b1_ref(13, 7) = xi * 24. * zk - zi * 18. * xk;
	b1_ref(13, 8) = b1_ref(13, 6);
	b1_ref(13, 9) = yi * -18. * yk - xi * 6. * xk - zi * 18. * zk;
	b1_ref(13, 10) = yi * 24. * zk - zi * 18. * yk;
	b1_ref(13, 11) = b1_ref(13, 7);
	b1_ref(13, 12) = b1_ref(13, 10);
	b1_ref(13, 13) = zi * 36. * zk + xi * 24. * xk + yi * 24. * yk;
	i__1 = isiz;
	for (fi = 1; fi <= i__1; ++fi) {
	    b1_ref(8, fi) = b1_ref(6, fi);
	    b1_ref(11, fi) = b1_ref(7, fi);
	    b1_ref(12, fi) = b1_ref(10, fi);
	}
	i__1 = isiz;
	for (fi = 1; fi <= i__1; ++fi) {
	    m2t2[fi - 1] = 0.;
	    i__2 = ksiz;
	    for (fj = 1; fj <= i__2; ++fj) {
		m2t2[fi - 1] += b1_ref(fi, fj) * rpk[fj];
	    }
	}
	term = 0.;
	i__1 = isiz;
	for (fi = 1; fi <= i__1; ++fi) {
	    term += rpi[fi] * m2t2[fi - 1];
	}
	term = size * 4. * term / (ratio * 4. + 3.);
	*eik += term;
	size *= size2;
    }

/*     recursive formulation of 4th through nth summation terms */

    i__1 = nn;
    for (fii = 4; fii <= i__1; ++fii) {
	i__2 = p_e1__;
	for (fi = 1; fi <= i__2; ++fi) {
	    if (fi == 8) {
		i__3 = p_e2__;
		for (fj = 1; fj <= i__3; ++fj) {
		    rklij_ref(fi, fj) = rklij_ref(6, fj);
		}
	    } else if (fi == 11) {
		i__3 = p_e2__;
		for (fj = 1; fj <= i__3; ++fj) {
		    rklij_ref(fi, fj) = rklij_ref(7, fj);
		}
	    } else if (fi == 12) {
		i__3 = p_e2__;
		for (fj = 1; fj <= i__3; ++fj) {
		    rklij_ref(fi, fj) = rklij_ref(10, fj);
		}
	    } else {
		i__3 = p_e2__;
		for (fj = 1; fj <= i__3; ++fj) {
		    if (fj == 8) {
			rklij_ref(fi, fj) = rklij_ref(fi, 6);
		    } else if (fj == 11) {
			rklij_ref(fi, fj) = rklij_ref(fi, 7);
		    } else if (fj == 12) {
			rklij_ref(fi, fj) = rklij_ref(fi, 10);
		    } else {
			rklij_ref(fi, fj) = d1d2_(&fii, &xi, &yi, &zi, &xk, &
				yk, &zk, &d__, &ri2, &rk2, &ind1_x__[fi - 1], 
				&ind1_y__[fi - 1], &ind1_z__[fi - 1], &
				ind2_x__[fj - 1], &ind2_y__[fj - 1], &
				ind2_z__[fj - 1]);
		    }
		}
	    }
	}

/*     update storage of the last two sets of matrix elements */

	i__2 = p_e1__;
	for (fi = 1; fi <= i__2; ++fi) {
	    i__3 = p_e2__;
	    for (fj = 1; fj <= i__3; ++fj) {
		b2_ref(fj, fi) = b1_ref(fj, fi);
		b1_ref(fj, fi) = rklij_ref(fj, fi);
	    }
	}

/*     compute interaction energy between the two multipole sites */

	i__2 = isiz;
	for (fi = 1; fi <= i__2; ++fi) {
	    m2t2[fi - 1] = 0.;
	    i__3 = ksiz;
	    for (fj = 1; fj <= i__3; ++fj) {
		m2t2[fi - 1] += rklij_ref(fi, fj) * rpk[fj];
	    }
	}
	term = 0.;
	i__2 = isiz;
	for (fi = 1; fi <= i__2; ++fi) {
	    term += rpi[fi] * m2t2[fi - 1];
	}
	term = term * size * (doublereal) (fii + 1) / ((doublereal) (fii + 1) 
		* ratio + (doublereal) fii);
	*eik += term;
	size *= size2;
    }
    *eik = factor * *eik;
    return 0;
} /* erfik_ */

#undef rklij_ref
#undef b2_ref
#undef b1_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine rfindex  --  reaction field indices for sites  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "rfindex" finds indices for each multipole site for use */
/*     in computing reaction field energetics */


/* Subroutine */ int rfindex_(integer *n, integer *m, integer *ind_x__, 
	integer *ind_y__, integer *ind_z__, integer *p_s__, integer *p_e__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k;



    /* Parameter adjustments */
    --ind_z__;
    --ind_y__;
    --ind_x__;

    /* Function Body */
    *p_s__ = 1;
    *p_e__ = 1;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ind_x__[i__] = 0;
	ind_y__[i__] = 0;
	ind_z__[i__] = 0;
    }
    k = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *p_e__;
	for (j = *p_s__; j <= i__2; ++j) {
	    ++k;
	    ind_x__[k] = ind_x__[j] + 1;
	    ind_y__[k] = ind_y__[j];
	    ind_z__[k] = ind_z__[j];
	}
	i__2 = *p_e__;
	for (j = *p_s__; j <= i__2; ++j) {
	    ++k;
	    ind_x__[k] = ind_x__[j];
	    ind_y__[k] = ind_y__[j] + 1;
	    ind_z__[k] = ind_z__[j];
	}
	i__2 = *p_e__;
	for (j = *p_s__; j <= i__2; ++j) {
	    ++k;
	    ind_x__[k] = ind_x__[j];
	    ind_y__[k] = ind_y__[j];
	    ind_z__[k] = ind_z__[j] + 1;
	}
	*p_s__ = *p_e__ + 1;
	*p_e__ = k;
    }
    return 0;
} /* rfindex_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine ijkpts  --  indices into "b1" and "b2" arrays  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "ijkpts" stores a set of indices used during calculation */
/*     of macroscopic reaction field energetics */

/* Subroutine */ int ijkpts_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer i__, j, k;


#define ijk_ref(a_1,a_2,a_3) rxnfld_1.ijk[((a_3)*6 + (a_2))*6 + a_1 - 0]



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  rxnfld.i  --  reaction field matrix elements and indices  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     b1    first reaction field matrix element array */
/*     b2    second reaction field matrix element array */
/*     ijk   indices into the reaction field element arrays */




/*     find and store indices into the "b1" and "b2" arrays */

    for (i__ = 0; i__ <= 5; ++i__) {
	for (j = 0; j <= 5; ++j) {
	    for (k = 0; k <= 5; ++k) {
		i__1 = i__ + j + k;
		i__2 = j + k;
		ijk_ref(i__, j, k) = (pow_ii(&c__3, &i__1) + pow_ii(&c__3, &
			i__2) + pow_ii(&c__3, &k) - 1) / 2;
	    }
	}
    }
    return 0;
} /* ijkpts_ */

#undef ijk_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function d1d2  --  recursive summation element utility  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "d1d2" is a utility function used in computation of the */
/*     reaction field recursive summation elements */


doublereal d1d2_(integer *n, doublereal *x1, doublereal *y1, doublereal *z1, 
	doublereal *x2, doublereal *y2, doublereal *z2, doublereal *d__, 
	doublereal *r1sq, doublereal *r2sq, integer *i__, integer *j, integer 
	*k, integer *s, integer *t, integer *u)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal f, g;
    static integer is, it, iu, js, jt, ju, ks, kt, ku;


#define b1_ref(a_1,a_2) rxnfld_1.b1[(a_2)*40 + a_1 - 41]
#define b2_ref(a_1,a_2) rxnfld_1.b2[(a_2)*40 + a_1 - 41]
#define ijk_ref(a_1,a_2,a_3) rxnfld_1.ijk[((a_3)*6 + (a_2))*6 + a_1 - 0]



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  rxnfld.i  --  reaction field matrix elements and indices  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     b1    first reaction field matrix element array */
/*     b2    second reaction field matrix element array */
/*     ijk   indices into the reaction field element arrays */




    if (*n < *i__ + *j + *k || *n < *s + *t + *u) {
	ret_val = 0.;
	return ret_val;
    }
    is = *i__ * *s;
    it = *i__ * *t;
    iu = *i__ * *u;
    js = *j * *s;
    jt = *j * *t;
    ju = *j * *u;
    ks = *k * *s;
    kt = *k * *t;
    ku = *k * *u;
    f = *d__ * b1_ref(ijk_ref(*i__, *j, *k), ijk_ref(*s, *t, *u));
    g = *r1sq * *r2sq * b2_ref(ijk_ref(*i__, *j, *k), ijk_ref(*s, *t, *u));
    if (*i__ != 0) {
	f += *i__ * *x2 * b1_ref(ijk_ref(*i__ - 1, *j, *k), ijk_ref(*s, *t, *
		u));
	g += *i__ * 2. * *x1 * *r2sq * b2_ref(ijk_ref(*i__ - 1, *j, *k), 
		ijk_ref(*s, *t, *u));
	if (*i__ != 1) {
	    g += *i__ * (*i__ - 1) * *r2sq * b2_ref(ijk_ref(*i__ - 2, *j, *k),
		     ijk_ref(*s, *t, *u));
	}
    }
    if (*j != 0) {
	f += *j * *y2 * b1_ref(ijk_ref(*i__, *j - 1, *k), ijk_ref(*s, *t, *u))
		;
	g += *j * 2. * *y1 * *r2sq * b2_ref(ijk_ref(*i__, *j - 1, *k), 
		ijk_ref(*s, *t, *u));
	if (*j != 1) {
	    g += *j * (*j - 1) * *r2sq * b2_ref(ijk_ref(*i__, *j - 2, *k), 
		    ijk_ref(*s, *t, *u));
	}
    }
    if (*k != 0) {
	f += *k * *z2 * b1_ref(ijk_ref(*i__, *j, *k - 1), ijk_ref(*s, *t, *u))
		;
	g += *k * 2. * *z1 * *r2sq * b2_ref(ijk_ref(*i__, *j, *k - 1), 
		ijk_ref(*s, *t, *u));
	if (*k != 1) {
	    g += *k * (*k - 1) * *r2sq * b2_ref(ijk_ref(*i__, *j, *k - 2), 
		    ijk_ref(*s, *t, *u));
	}
    }
    if (*s != 0) {
	f += *s * *x1 * b1_ref(ijk_ref(*i__, *j, *k), ijk_ref(*s - 1, *t, *u))
		;
	g += *s * 2. * *x2 * *r1sq * b2_ref(ijk_ref(*i__, *j, *k), ijk_ref(*s 
		- 1, *t, *u));
	if (*s != 1) {
	    g += *s * (*s - 1) * *r1sq * b2_ref(ijk_ref(*i__, *j, *k), 
		    ijk_ref(*s - 2, *t, *u));
	}
    }
    if (*t != 0) {
	f += *t * *y1 * b1_ref(ijk_ref(*i__, *j, *k), ijk_ref(*s, *t - 1, *u))
		;
	g += *t * 2. * *y2 * *r1sq * b2_ref(ijk_ref(*i__, *j, *k), ijk_ref(*s,
		 *t - 1, *u));
	if (*t != 1) {
	    g += *t * (*t - 1) * *r1sq * b2_ref(ijk_ref(*i__, *j, *k), 
		    ijk_ref(*s, *t - 2, *u));
	}
    }
    if (*u != 0) {
	f += *u * *z1 * b1_ref(ijk_ref(*i__, *j, *k), ijk_ref(*s, *t, *u - 1))
		;
	g += *u * 2. * *z2 * *r1sq * b2_ref(ijk_ref(*i__, *j, *k), ijk_ref(*s,
		 *t, *u - 1));
	if (*u != 1) {
	    g += *u * (*u - 1) * *r1sq * b2_ref(ijk_ref(*i__, *j, *k), 
		    ijk_ref(*s, *t, *u - 2));
	}
    }
    if (is != 0) {
	f += is * b1_ref(ijk_ref(*i__ - 1, *j, *k), ijk_ref(*s - 1, *t, *u));
	g += is * 4. * *x1 * *x2 * b2_ref(ijk_ref(*i__ - 1, *j, *k), ijk_ref(*
		s - 1, *t, *u));
	if (*i__ != 1) {
	    g += (*i__ - 1) * 2. * is * *x2 * b2_ref(ijk_ref(*i__ - 2, *j, *k)
		    , ijk_ref(*s - 1, *t, *u));
	    if (*s != 1) {
		g += (*i__ - 1) * (*s - 1) * is * b2_ref(ijk_ref(*i__ - 2, *j,
			 *k), ijk_ref(*s - 2, *t, *u));
	    }
	}
	if (*s != 1) {
	    g += (*s - 1) * 2. * is * *x1 * b2_ref(ijk_ref(*i__ - 1, *j, *k), 
		    ijk_ref(*s - 2, *t, *u));
	}
    }
    if (jt != 0) {
	f += jt * b1_ref(ijk_ref(*i__, *j - 1, *k), ijk_ref(*s, *t - 1, *u));
	g += jt * 4. * *y1 * *y2 * b2_ref(ijk_ref(*i__, *j - 1, *k), ijk_ref(*
		s, *t - 1, *u));
	if (*j != 1) {
	    g += (*j - 1) * 2. * jt * *y2 * b2_ref(ijk_ref(*i__, *j - 2, *k), 
		    ijk_ref(*s, *t - 1, *u));
	    if (*t != 1) {
		g += (*j - 1) * (*t - 1) * jt * b2_ref(ijk_ref(*i__, *j - 2, *
			k), ijk_ref(*s, *t - 2, *u));
	    }
	}
	if (*t != 1) {
	    g += (*t - 1) * 2. * jt * *y1 * b2_ref(ijk_ref(*i__, *j - 1, *k), 
		    ijk_ref(*s, *t - 2, *u));
	}
    }
    if (ku != 0) {
	f += ku * b1_ref(ijk_ref(*i__, *j, *k - 1), ijk_ref(*s, *t, *u - 1));
	g += ku * 4. * *z1 * *z2 * b2_ref(ijk_ref(*i__, *j, *k - 1), ijk_ref(*
		s, *t, *u - 1));
	if (*k != 1) {
	    g += (*k - 1) * 2. * ku * *z2 * b2_ref(ijk_ref(*i__, *j, *k - 2), 
		    ijk_ref(*s, *t, *u - 1));
	    if (*u != 1) {
		g += (*k - 1) * (*u - 1) * ku * b2_ref(ijk_ref(*i__, *j, *k - 
			2), ijk_ref(*s, *t, *u - 2));
	    }
	}
	if (*u != 1) {
	    g += (*u - 1) * 2. * ku * *z1 * b2_ref(ijk_ref(*i__, *j, *k - 1), 
		    ijk_ref(*s, *t, *u - 2));
	}
    }
    if (it != 0) {
	g += it * 4. * *x1 * *y2 * b2_ref(ijk_ref(*i__ - 1, *j, *k), ijk_ref(*
		s, *t - 1, *u));
	if (*i__ != 1) {
	    g += (*i__ - 1) * 2. * it * *y2 * b2_ref(ijk_ref(*i__ - 2, *j, *k)
		    , ijk_ref(*s, *t - 1, *u));
	    if (*t != 1) {
		g += (*i__ - 1) * (*t - 1) * it * b2_ref(ijk_ref(*i__ - 2, *j,
			 *k), ijk_ref(*s, *t - 2, *u));
	    }
	}
	if (*t != 1) {
	    g += (*t - 1) * 2. * it * *x1 * b2_ref(ijk_ref(*i__ - 1, *j, *k), 
		    ijk_ref(*s, *t - 2, *u));
	}
    }
    if (iu != 0) {
	g += iu * 4. * *x1 * *z2 * b2_ref(ijk_ref(*i__ - 1, *j, *k), ijk_ref(*
		s, *t, *u - 1));
	if (*i__ != 1) {
	    g += (*i__ - 1) * 2. * iu * *z2 * b2_ref(ijk_ref(*i__ - 2, *j, *k)
		    , ijk_ref(*s, *t, *u - 1));
	    if (*u != 1) {
		g += (*i__ - 1) * (*u - 1) * iu * b2_ref(ijk_ref(*i__ - 2, *j,
			 *k), ijk_ref(*s, *t, *u - 2));
	    }
	}
	if (*u != 1) {
	    g += (*u - 1) * 2. * iu * *x1 * b2_ref(ijk_ref(*i__ - 1, *j, *k), 
		    ijk_ref(*s, *t, *u - 2));
	}
    }
    if (js != 0) {
	g += js * 4. * *y1 * *x2 * b2_ref(ijk_ref(*i__, *j - 1, *k), ijk_ref(*
		s - 1, *t, *u));
	if (*j != 1) {
	    g += (*j - 1) * 2. * js * *x2 * b2_ref(ijk_ref(*i__, *j - 2, *k), 
		    ijk_ref(*s - 1, *t, *u));
	    if (*s != 1) {
		g += (*j - 1) * (*s - 1) * js * b2_ref(ijk_ref(*i__, *j - 2, *
			k), ijk_ref(*s - 2, *t, *u));
	    }
	}
	if (*s != 1) {
	    g += (*s - 1) * 2. * js * *y1 * b2_ref(ijk_ref(*i__, *j - 1, *k), 
		    ijk_ref(*s - 2, *t, *u));
	}
    }
    if (ju != 0) {
	g += ju * 4. * *y1 * *z2 * b2_ref(ijk_ref(*i__, *j - 1, *k), ijk_ref(*
		s, *t, *u - 1));
	if (*j != 1) {
	    g += (*j - 1) * 2. * ju * *z2 * b2_ref(ijk_ref(*i__, *j - 2, *k), 
		    ijk_ref(*s, *t, *u - 1));
	    if (*u != 1) {
		g += (*j - 1) * (*u - 1) * ju * b2_ref(ijk_ref(*i__, *j - 2, *
			k), ijk_ref(*s, *t, *u - 2));
	    }
	}
	if (*u != 1) {
	    g += (*u - 1) * 2. * ju * *y1 * b2_ref(ijk_ref(*i__, *j - 1, *k), 
		    ijk_ref(*s, *t, *u - 2));
	}
    }
    if (ks != 0) {
	g += ks * 4. * *z1 * *x2 * b2_ref(ijk_ref(*i__, *j, *k - 1), ijk_ref(*
		s - 1, *t, *u));
	if (*k != 1) {
	    g += (*k - 1) * 2. * ks * *x2 * b2_ref(ijk_ref(*i__, *j, *k - 2), 
		    ijk_ref(*s - 1, *t, *u));
	    if (*s != 1) {
		g += (*k - 1) * (*s - 1) * ks * b2_ref(ijk_ref(*i__, *j, *k - 
			2), ijk_ref(*s - 2, *t, *u));
	    }
	}
	if (*s != 1) {
	    g += (*s - 1) * 2. * ks * *z1 * b2_ref(ijk_ref(*i__, *j, *k - 1), 
		    ijk_ref(*s - 2, *t, *u));
	}
    }
    if (kt != 0) {
	g += kt * 4. * *z1 * *y2 * b2_ref(ijk_ref(*i__, *j, *k - 1), ijk_ref(*
		s, *t - 1, *u));
	if (*k != 1) {
	    g += (*k - 1) * 2. * kt * *y2 * b2_ref(ijk_ref(*i__, *j, *k - 2), 
		    ijk_ref(*s, *t - 1, *u));
	    if (*t != 1) {
		g += (*k - 1) * (*t - 1) * kt * b2_ref(ijk_ref(*i__, *j, *k - 
			2), ijk_ref(*s, *t - 2, *u));
	    }
	}
	if (*t != 1) {
	    g += (*t - 1) * 2. * kt * *z1 * b2_ref(ijk_ref(*i__, *j, *k - 1), 
		    ijk_ref(*s, *t - 2, *u));
	}
    }
    f = (doublereal) ((*n << 1) - 1) * f;
    g = (doublereal) (*n - 1) * g;
    ret_val = (f - g) / (doublereal) (*n);
    return ret_val;
} /* d1d2_ */

#undef ijk_ref
#undef b2_ref
#undef b1_ref


