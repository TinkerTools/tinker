/* switch.f -- translated by f2c (version 20050501).
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
    doublereal vdwcut, chgcut, dplcut, mpolecut, vdwtaper, chgtaper, dpltaper,
	     mpoletaper, ewaldcut;
    logical use_ewald__, use_lights__, use_list__, use_vlist__, use_clist__, 
	    use_mlist__;
} cutoff_;

#define cutoff_1 cutoff_

struct {
    doublereal solvprs, surften, spcut, spoff, stcut, stoff, rcav[25000], 
	    rdisp[25000], cdisp[25000];
} npolar_;

#define npolar_1 npolar_

struct {
    doublereal off, off2, cut, cut2, c0, c1, c2, c3, c4, c5, f0, f1, f2, f3, 
	    f4, f5, f6, f7;
} shunt_;

#define shunt_1 shunt_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine switch  --  get switching function coefficients  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "switch" sets the coeffcients used by the fifth and seventh */
/*     order polynomial switching functions for spherical cutoffs */


/* Subroutine */ int switch_(char *mode, ftnlen mode_len)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static doublereal off3, off4, off5, off6, off7, cut3, cut4, cut5, cut6, 
	    cut7, term, denom;
    extern /* Subroutine */ int replica_(doublereal *);



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




/*     get the switching window for the current potential type */

    if (s_cmp(mode, "VDW", (ftnlen)3, (ftnlen)3) == 0) {
	shunt_1.off = cutoff_1.vdwcut;
	shunt_1.cut = cutoff_1.vdwtaper;
    } else if (s_cmp(mode, "CHARGE", (ftnlen)6, (ftnlen)6) == 0) {
	shunt_1.off = cutoff_1.chgcut;
	shunt_1.cut = cutoff_1.chgtaper;
    } else if (s_cmp(mode, "CHGDPL", (ftnlen)6, (ftnlen)6) == 0) {
	shunt_1.off = sqrt(cutoff_1.chgcut * cutoff_1.dplcut);
	shunt_1.cut = sqrt(cutoff_1.chgtaper * cutoff_1.dpltaper);
    } else if (s_cmp(mode, "DIPOLE", (ftnlen)6, (ftnlen)6) == 0) {
	shunt_1.off = cutoff_1.dplcut;
	shunt_1.cut = cutoff_1.dpltaper;
    } else if (s_cmp(mode, "MPOLE", (ftnlen)5, (ftnlen)5) == 0) {
	shunt_1.off = cutoff_1.mpolecut;
	shunt_1.cut = cutoff_1.mpoletaper;
    } else if (s_cmp(mode, "EWALD", (ftnlen)5, (ftnlen)5) == 0) {
	shunt_1.off = cutoff_1.ewaldcut;
	shunt_1.cut = cutoff_1.ewaldcut;
    } else if (s_cmp(mode, "GKV", (ftnlen)3, (ftnlen)3) == 0) {
	shunt_1.off = npolar_1.spoff;
	shunt_1.cut = npolar_1.spcut;
    } else if (s_cmp(mode, "GKSA", (ftnlen)4, (ftnlen)4) == 0) {
	shunt_1.off = npolar_1.stcut;
	shunt_1.cut = npolar_1.stoff;
    } else {
/* Computing MIN */
	d__1 = min(cutoff_1.vdwcut,cutoff_1.chgcut), d__1 = min(d__1,
		cutoff_1.dplcut);
	shunt_1.off = min(d__1,cutoff_1.mpolecut);
/* Computing MIN */
	d__1 = min(cutoff_1.vdwtaper,cutoff_1.chgtaper), d__1 = min(d__1,
		cutoff_1.dpltaper);
	shunt_1.cut = min(d__1,cutoff_1.mpoletaper);
    }

/*     test for replicate periodic boundaries at this cutoff */

    replica_(&shunt_1.off);

/*     set switching coefficients to zero for truncation cutoffs */

    shunt_1.c0 = 0.;
    shunt_1.c1 = 0.;
    shunt_1.c2 = 0.;
    shunt_1.c3 = 0.;
    shunt_1.c4 = 0.;
    shunt_1.c5 = 0.;
    shunt_1.f0 = 0.;
    shunt_1.f1 = 0.;
    shunt_1.f2 = 0.;
    shunt_1.f3 = 0.;
    shunt_1.f4 = 0.;
    shunt_1.f5 = 0.;
    shunt_1.f6 = 0.;
    shunt_1.f7 = 0.;

/*     store the powers of the switching window cutoffs */

    shunt_1.off2 = shunt_1.off * shunt_1.off;
    off3 = shunt_1.off2 * shunt_1.off;
    off4 = shunt_1.off2 * shunt_1.off2;
    off5 = shunt_1.off2 * off3;
    off6 = off3 * off3;
    off7 = off3 * off4;
    shunt_1.cut2 = shunt_1.cut * shunt_1.cut;
    cut3 = shunt_1.cut2 * shunt_1.cut;
    cut4 = shunt_1.cut2 * shunt_1.cut2;
    cut5 = shunt_1.cut2 * cut3;
    cut6 = cut3 * cut3;
    cut6 = cut3 * cut3;
    cut7 = cut3 * cut4;

/*     get 5th degree multiplicative switching function coefficients */

    if (shunt_1.cut < shunt_1.off) {
/* Computing 5th power */
	d__1 = shunt_1.off - shunt_1.cut, d__2 = d__1, d__1 *= d__1;
	denom = d__2 * (d__1 * d__1);
	shunt_1.c0 = shunt_1.off * shunt_1.off2 * (shunt_1.off2 - shunt_1.off 
		* 5. * shunt_1.cut + shunt_1.cut2 * 10.) / denom;
	shunt_1.c1 = shunt_1.off2 * -30. * shunt_1.cut2 / denom;
	shunt_1.c2 = (shunt_1.off2 * shunt_1.cut + shunt_1.off * shunt_1.cut2)
		 * 30. / denom;
	shunt_1.c3 = (shunt_1.off2 + shunt_1.off * 4. * shunt_1.cut + 
		shunt_1.cut2) * -10. / denom;
	shunt_1.c4 = (shunt_1.off + shunt_1.cut) * 15. / denom;
	shunt_1.c5 = -6. / denom;
    }

/*     get 7th degree additive switching function coefficients */

    if (shunt_1.cut < shunt_1.off && s_cmp(mode, "CHARGE", (ftnlen)6, (ftnlen)
	    6) == 0) {
	term = shunt_1.cut * 9.3 * shunt_1.off / (shunt_1.off - shunt_1.cut);
	denom = cut7 - cut6 * 7. * shunt_1.off + cut5 * 21. * shunt_1.off2 - 
		cut4 * 35. * off3 + cut3 * 35. * off4 - shunt_1.cut2 * 21. * 
		off5 + shunt_1.cut * 7. * off6 - off7;
	denom = term * denom;
	shunt_1.f0 = cut3 * off3 * (shunt_1.cut * -39. + shunt_1.off * 64.) / 
		denom;
	shunt_1.f1 = shunt_1.cut2 * shunt_1.off2 * (shunt_1.cut2 * 117. - 
		shunt_1.cut * 100. * shunt_1.off - shunt_1.off2 * 192.) / 
		denom;
	shunt_1.f2 = shunt_1.cut * shunt_1.off * (cut3 * -117. - shunt_1.cut2 
		* 84. * shunt_1.off + shunt_1.cut * 534. * shunt_1.off2 + 
		off3 * 192.) / denom;
	shunt_1.f3 = (cut4 * 39. + cut3 * 212. * shunt_1.off - shunt_1.cut2 * 
		450. * shunt_1.off2 - shunt_1.cut * 612. * off3 - off4 * 64.) 
		/ denom;
	shunt_1.f4 = (cut3 * -92. + shunt_1.cut2 * 66. * shunt_1.off + 
		shunt_1.cut * 684. * shunt_1.off2 + off3 * 217.) / denom;
	shunt_1.f5 = (shunt_1.cut2 * 42. - shunt_1.cut * 300. * shunt_1.off - 
		shunt_1.off2 * 267.) / denom;
	shunt_1.f6 = (shunt_1.cut * 36. + shunt_1.off * 139.) / denom;
	shunt_1.f7 = -25. / denom;
    }
    return 0;
} /* switch_ */

