/* etors1.f -- translated by f2c (version 20050501).
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
    doublereal m2[25000], deform, difft, diffv, diffc;
    logical use_smooth__, use_dem__, use_gda__, use_tophat__, use_stophat__;
} warp_;

#define warp_1 warp_

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
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    doublereal idihunit, itorunit, torsunit, ptorunit, storunit, ttorunit;
} torpot_;

#define torpot_1 torpot_

struct {
    doublereal tors1[400000]	/* was [4][100000] */, tors2[400000]	/* 
	    was [4][100000] */, tors3[400000]	/* was [4][100000] */, tors4[
	    400000]	/* was [4][100000] */, tors5[400000]	/* was [4][
	    100000] */, tors6[400000]	/* was [4][100000] */;
    integer ntors, itors[400000]	/* was [4][100000] */;
} tors_;

#define tors_1 tors_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal vir[9]	/* was [3][3] */;
} virial_;

#define virial_1 virial_

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine etors1  --  torsional energy & derivatives  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "etors1" calculates the torsional potential energy and first */
/*     derivatives with respect to Cartesian coordinates */


/* Subroutine */ int etors1_(void)
{
    extern /* Subroutine */ int etors1a_(void), etors1b_(void);



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




/*     choose standard or potential energy smoothing version */



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


    if (warp_1.use_smooth__) {
	etors1b_();
    } else {
	etors1a_();
    }
    return 0;
} /* etors1_ */



/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine etors1a  --  standard torsional energy & derivs  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "etors1a" calculates the torsional potential energy and first */
/*     derivatives with respect to Cartesian coordinates using a */
/*     standard sum of Fourier terms */


/* Subroutine */ int etors1a_(void)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e;
    static integer i__;
    static doublereal c1, c2, c3, c4, c5, c6, s1, s2, s3, v1, v2, v3, v4, v5, 
	    v6, s4, s5, s6;
    static integer ia, ib, ic, id;
    static doublereal xt, yt, zt, xu, yu, zu, rt2, ru2, rcb, xia, yia, zia, 
	    xib, yib, zib, xic, yic, zic, xid, yid, zid, xba, yba, zba, xdc, 
	    ydc, zdc, xcb, ycb, zcb, xca, yca, zca, xdb, ydb, zdb, xtu, ytu, 
	    ztu, vxx, vyx, vyy, vzx, vzz, vzy, phi1, phi2, phi3, phi4, phi5, 
	    phi6, fgrp, sine, rtru, dphi1, dphi2, dphi3, dphi4, dphi5, dphi6, 
	    sine2, sine3, sine4, sine5, sine6;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal dedxt, dedyt, dedzt, dedxu, dedyu, dedzu, dedphi, 
	    dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, dedxic, dedyic, 
	    dedzic, dedxid, dedyid, dedzid, cosine;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal cosine2, cosine3, cosine4, cosine5, cosine6;
    static logical proceed;


#define det_ref(a_1,a_2) deriv_1.det[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define tors1_ref(a_1,a_2) tors_1.tors1[(a_2)*4 + a_1 - 5]
#define tors2_ref(a_1,a_2) tors_1.tors2[(a_2)*4 + a_1 - 5]
#define tors3_ref(a_1,a_2) tors_1.tors3[(a_2)*4 + a_1 - 5]
#define tors4_ref(a_1,a_2) tors_1.tors4[(a_2)*4 + a_1 - 5]
#define tors5_ref(a_1,a_2) tors_1.tors5[(a_2)*4 + a_1 - 5]
#define tors6_ref(a_1,a_2) tors_1.tors6[(a_2)*4 + a_1 - 5]
#define itors_ref(a_1,a_2) tors_1.itors[(a_2)*4 + a_1 - 5]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  torpot.i  --  specifics of torsional functional forms  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     idihunit  convert improper dihedral energy to kcal/mole */
/*     itorunit  convert improper torsion amplitudes to kcal/mole */
/*     torsunit  convert torsional parameter amplitudes to kcal/mole */
/*     ptorunit  convert pi-orbital torsion energy to kcal/mole */
/*     storunit  convert stretch-torsion energy to kcal/mole */
/*     ttorunit  convert stretch-torsion energy to kcal/mole */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  tors.i  --  torsional angles within the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     tors1   1-fold amplitude and phase for each torsional angle */
/*     tors2   2-fold amplitude and phase for each torsional angle */
/*     tors3   3-fold amplitude and phase for each torsional angle */
/*     tors4   4-fold amplitude and phase for each torsional angle */
/*     tors5   5-fold amplitude and phase for each torsional angle */
/*     tors6   6-fold amplitude and phase for each torsional angle */
/*     ntors   total number of torsional angles in the system */
/*     itors   numbers of the atoms in each torsional angle */




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




/*     zero out the torsional energy and first derivatives */

    energi_1.et = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	det_ref(1, i__) = 0.;
	det_ref(2, i__) = 0.;
	det_ref(3, i__) = 0.;
    }

/*     calculate the torsional angle energy and first derivatives */

    i__1 = tors_1.ntors;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = itors_ref(1, i__);
	ib = itors_ref(2, i__);
	ic = itors_ref(3, i__);
	id = itors_ref(4, i__);

/*     decide whether to compute the current interaction */

	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &c__0, &c__0);
	}
	if (proceed) {
	    proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
		    usage_1.use[ic - 1] || usage_1.use[id - 1];
	}

/*     compute the value of the torsional angle */

	if (proceed) {
	    xia = atoms_1.x[ia - 1];
	    yia = atoms_1.y[ia - 1];
	    zia = atoms_1.z__[ia - 1];
	    xib = atoms_1.x[ib - 1];
	    yib = atoms_1.y[ib - 1];
	    zib = atoms_1.z__[ib - 1];
	    xic = atoms_1.x[ic - 1];
	    yic = atoms_1.y[ic - 1];
	    zic = atoms_1.z__[ic - 1];
	    xid = atoms_1.x[id - 1];
	    yid = atoms_1.y[id - 1];
	    zid = atoms_1.z__[id - 1];
	    xba = xib - xia;
	    yba = yib - yia;
	    zba = zib - zia;
	    xcb = xic - xib;
	    ycb = yic - yib;
	    zcb = zic - zib;
	    xdc = xid - xic;
	    ydc = yid - yic;
	    zdc = zid - zic;
	    if (bound_1.use_polymer__) {
		image_(&xba, &yba, &zba);
		image_(&xcb, &ycb, &zcb);
		image_(&xdc, &ydc, &zdc);
	    }
	    xt = yba * zcb - ycb * zba;
	    yt = zba * xcb - zcb * xba;
	    zt = xba * ycb - xcb * yba;
	    xu = ycb * zdc - ydc * zcb;
	    yu = zcb * xdc - zdc * xcb;
	    zu = xcb * ydc - xdc * ycb;
	    xtu = yt * zu - yu * zt;
	    ytu = zt * xu - zu * xt;
	    ztu = xt * yu - xu * yt;
	    rt2 = xt * xt + yt * yt + zt * zt;
	    ru2 = xu * xu + yu * yu + zu * zu;
	    rtru = sqrt(rt2 * ru2);
	    if (rtru != 0.) {
		rcb = sqrt(xcb * xcb + ycb * ycb + zcb * zcb);
		cosine = (xt * xu + yt * yu + zt * zu) / rtru;
		sine = (xcb * xtu + ycb * ytu + zcb * ztu) / (rcb * rtru);

/*     set the torsional parameters for this angle */

		v1 = tors1_ref(1, i__);
		c1 = tors1_ref(3, i__);
		s1 = tors1_ref(4, i__);
		v2 = tors2_ref(1, i__);
		c2 = tors2_ref(3, i__);
		s2 = tors2_ref(4, i__);
		v3 = tors3_ref(1, i__);
		c3 = tors3_ref(3, i__);
		s3 = tors3_ref(4, i__);
		v4 = tors4_ref(1, i__);
		c4 = tors4_ref(3, i__);
		s4 = tors4_ref(4, i__);
		v5 = tors5_ref(1, i__);
		c5 = tors5_ref(3, i__);
		s5 = tors5_ref(4, i__);
		v6 = tors6_ref(1, i__);
		c6 = tors6_ref(3, i__);
		s6 = tors6_ref(4, i__);

/*     compute the multiple angle trigonometry and the phase terms */

		cosine2 = cosine * cosine - sine * sine;
		sine2 = cosine * 2. * sine;
		cosine3 = cosine * cosine2 - sine * sine2;
		sine3 = cosine * sine2 + sine * cosine2;
		cosine4 = cosine * cosine3 - sine * sine3;
		sine4 = cosine * sine3 + sine * cosine3;
		cosine5 = cosine * cosine4 - sine * sine4;
		sine5 = cosine * sine4 + sine * cosine4;
		cosine6 = cosine * cosine5 - sine * sine5;
		sine6 = cosine * sine5 + sine * cosine5;
		phi1 = cosine * c1 + sine * s1 + 1.;
		phi2 = cosine2 * c2 + sine2 * s2 + 1.;
		phi3 = cosine3 * c3 + sine3 * s3 + 1.;
		phi4 = cosine4 * c4 + sine4 * s4 + 1.;
		phi5 = cosine5 * c5 + sine5 * s5 + 1.;
		phi6 = cosine6 * c6 + sine6 * s6 + 1.;
		dphi1 = cosine * s1 - sine * c1;
		dphi2 = (cosine2 * s2 - sine2 * c2) * 2.;
		dphi3 = (cosine3 * s3 - sine3 * c3) * 3.;
		dphi4 = (cosine4 * s4 - sine4 * c4) * 4.;
		dphi5 = (cosine5 * s5 - sine5 * c5) * 5.;
		dphi6 = (cosine6 * s6 - sine6 * c6) * 6.;

/*     calculate torsional energy and master chain rule term */

		e = torpot_1.torsunit * (v1 * phi1 + v2 * phi2 + v3 * phi3 + 
			v4 * phi4 + v5 * phi5 + v6 * phi6);
		dedphi = torpot_1.torsunit * (v1 * dphi1 + v2 * dphi2 + v3 * 
			dphi3 + v4 * dphi4 + v5 * dphi5 + v6 * dphi6);

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    e *= fgrp;
		    dedphi *= fgrp;
		}

/*     chain rule terms for first derivative components */

		xca = xic - xia;
		yca = yic - yia;
		zca = zic - zia;
		xdb = xid - xib;
		ydb = yid - yib;
		zdb = zid - zib;
		if (bound_1.use_polymer__) {
		    image_(&xca, &yca, &zca);
		    image_(&xdb, &ydb, &zdb);
		}
		dedxt = dedphi * (yt * zcb - ycb * zt) / (rt2 * rcb);
		dedyt = dedphi * (zt * xcb - zcb * xt) / (rt2 * rcb);
		dedzt = dedphi * (xt * ycb - xcb * yt) / (rt2 * rcb);
		dedxu = -dedphi * (yu * zcb - ycb * zu) / (ru2 * rcb);
		dedyu = -dedphi * (zu * xcb - zcb * xu) / (ru2 * rcb);
		dedzu = -dedphi * (xu * ycb - xcb * yu) / (ru2 * rcb);

/*     compute first derivative components for this angle */

		dedxia = zcb * dedyt - ycb * dedzt;
		dedyia = xcb * dedzt - zcb * dedxt;
		dedzia = ycb * dedxt - xcb * dedyt;
		dedxib = yca * dedzt - zca * dedyt + zdc * dedyu - ydc * 
			dedzu;
		dedyib = zca * dedxt - xca * dedzt + xdc * dedzu - zdc * 
			dedxu;
		dedzib = xca * dedyt - yca * dedxt + ydc * dedxu - xdc * 
			dedyu;
		dedxic = zba * dedyt - yba * dedzt + ydb * dedzu - zdb * 
			dedyu;
		dedyic = xba * dedzt - zba * dedxt + zdb * dedxu - xdb * 
			dedzu;
		dedzic = yba * dedxt - xba * dedyt + xdb * dedyu - ydb * 
			dedxu;
		dedxid = zcb * dedyu - ycb * dedzu;
		dedyid = xcb * dedzu - zcb * dedxu;
		dedzid = ycb * dedxu - xcb * dedyu;

/*     increment the total torsional angle energy and gradient */

		energi_1.et += e;
		det_ref(1, ia) = det_ref(1, ia) + dedxia;
		det_ref(2, ia) = det_ref(2, ia) + dedyia;
		det_ref(3, ia) = det_ref(3, ia) + dedzia;
		det_ref(1, ib) = det_ref(1, ib) + dedxib;
		det_ref(2, ib) = det_ref(2, ib) + dedyib;
		det_ref(3, ib) = det_ref(3, ib) + dedzib;
		det_ref(1, ic) = det_ref(1, ic) + dedxic;
		det_ref(2, ic) = det_ref(2, ic) + dedyic;
		det_ref(3, ic) = det_ref(3, ic) + dedzic;
		det_ref(1, id) = det_ref(1, id) + dedxid;
		det_ref(2, id) = det_ref(2, id) + dedyid;
		det_ref(3, id) = det_ref(3, id) + dedzid;

/*     increment the internal virial tensor components */

		vxx = xcb * (dedxic + dedxid) - xba * dedxia + xdc * dedxid;
		vyx = ycb * (dedxic + dedxid) - yba * dedxia + ydc * dedxid;
		vzx = zcb * (dedxic + dedxid) - zba * dedxia + zdc * dedxid;
		vyy = ycb * (dedyic + dedyid) - yba * dedyia + ydc * dedyid;
		vzy = zcb * (dedyic + dedyid) - zba * dedyia + zdc * dedyid;
		vzz = zcb * (dedzic + dedzid) - zba * dedzia + zdc * dedzid;
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
	}
    }
    return 0;
} /* etors1a_ */

#undef itors_ref
#undef tors6_ref
#undef tors5_ref
#undef tors4_ref
#undef tors3_ref
#undef tors2_ref
#undef tors1_ref
#undef vir_ref
#undef det_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine etors1b  --  smoothed torsional energy & derivs  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "etors1b" calculates the torsional potential energy and first */
/*     derivatives with respect to Cartesian coordinates for use with */
/*     potential energy smoothing methods */


/* Subroutine */ int etors1b_(void)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double exp(doublereal), sin(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal e;
    static integer i__;
    static doublereal c1, c2, c3, c4, c5, c6, s1, s2, s3, v1, v2, v3, v4, v5, 
	    v6, s4, s5, s6;
    static integer ia, ib, ic, id;
    static doublereal xt, yt, zt, xu, yu, zu, rt2, ru2, rcb, xia, yia, zia, 
	    xib, yib, zib, xic, yic, zic, xid, yid, zid, xba, yba, zba, xdc, 
	    ydc, zdc, xcb, ycb, zcb, xca, yca, zca, xdb, ydb, zdb, xtu, ytu, 
	    ztu, vxx, vyx, vyy, vzx, vzz, vzy, phi1, phi2, phi3, phi4, phi5, 
	    phi6, fgrp, sine, rtru, damp1, damp2, damp3, damp4, damp5, damp6, 
	    dphi1, dphi2, dphi3, dphi4, dphi5, dphi6, sine2, sine3, sine4, 
	    sine5, sine6, dedxt, dedyt, dedzt, dedxu, dedyu, width, dedzu, 
	    wterm, dedphi, dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, 
	    dedxic, dedyic, dedzic, dedxid, dedyid, dedzid, cosine;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal cosine2, cosine3, cosine4, cosine5, cosine6;
    static logical proceed;


#define det_ref(a_1,a_2) deriv_1.det[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define tors1_ref(a_1,a_2) tors_1.tors1[(a_2)*4 + a_1 - 5]
#define tors2_ref(a_1,a_2) tors_1.tors2[(a_2)*4 + a_1 - 5]
#define tors3_ref(a_1,a_2) tors_1.tors3[(a_2)*4 + a_1 - 5]
#define tors4_ref(a_1,a_2) tors_1.tors4[(a_2)*4 + a_1 - 5]
#define tors5_ref(a_1,a_2) tors_1.tors5[(a_2)*4 + a_1 - 5]
#define tors6_ref(a_1,a_2) tors_1.tors6[(a_2)*4 + a_1 - 5]
#define itors_ref(a_1,a_2) tors_1.itors[(a_2)*4 + a_1 - 5]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  torpot.i  --  specifics of torsional functional forms  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     idihunit  convert improper dihedral energy to kcal/mole */
/*     itorunit  convert improper torsion amplitudes to kcal/mole */
/*     torsunit  convert torsional parameter amplitudes to kcal/mole */
/*     ptorunit  convert pi-orbital torsion energy to kcal/mole */
/*     storunit  convert stretch-torsion energy to kcal/mole */
/*     ttorunit  convert stretch-torsion energy to kcal/mole */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  tors.i  --  torsional angles within the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     tors1   1-fold amplitude and phase for each torsional angle */
/*     tors2   2-fold amplitude and phase for each torsional angle */
/*     tors3   3-fold amplitude and phase for each torsional angle */
/*     tors4   4-fold amplitude and phase for each torsional angle */
/*     tors5   5-fold amplitude and phase for each torsional angle */
/*     tors6   6-fold amplitude and phase for each torsional angle */
/*     ntors   total number of torsional angles in the system */
/*     itors   numbers of the atoms in each torsional angle */




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




/*     zero out the torsional energy and first derivatives */

    energi_1.et = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	det_ref(1, i__) = 0.;
	det_ref(2, i__) = 0.;
	det_ref(3, i__) = 0.;
    }

/*     set the extent of smoothing to be performed */

    width = warp_1.difft * warp_1.deform;
    if (width <= 0.) {
	damp1 = 1.;
	damp2 = 1.;
	damp3 = 1.;
	damp4 = 1.;
	damp5 = 1.;
	damp6 = 1.;
    } else if (warp_1.use_dem__) {
	damp1 = exp(-width);
	damp2 = exp(width * -4.);
	damp3 = exp(width * -9.);
	damp4 = exp(width * -16.);
	damp5 = exp(width * -25.);
	damp6 = exp(width * -36.);
    } else if (warp_1.use_gda__) {
	wterm = warp_1.difft / 12.;
    } else if (warp_1.use_tophat__ || warp_1.use_stophat__) {
	damp1 = 0.;
	damp2 = 0.;
	damp3 = 0.;
	damp4 = 0.;
	damp5 = 0.;
	damp6 = 0.;
	if (width < 3.141592653589793238) {
	    damp1 = sin(width) / width;
	}
	wterm = width * 2.;
	if (wterm < 3.141592653589793238) {
	    damp2 = sin(wterm) / wterm;
	}
	wterm = width * 3.;
	if (wterm < 3.141592653589793238) {
	    damp3 = sin(wterm) / wterm;
	}
	wterm = width * 4.;
	if (wterm < 3.141592653589793238) {
	    damp4 = sin(wterm) / wterm;
	}
	wterm = width * 5.;
	if (wterm < 3.141592653589793238) {
	    damp5 = sin(wterm) / wterm;
	}
	wterm = width * 6.;
	if (wterm < 3.141592653589793238) {
	    damp6 = sin(wterm) / wterm;
	}
    }

/*     calculate the torsional angle energy and first derivatives */

    i__1 = tors_1.ntors;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = itors_ref(1, i__);
	ib = itors_ref(2, i__);
	ic = itors_ref(3, i__);
	id = itors_ref(4, i__);

/*     decide whether to compute the current interaction */

	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &c__0, &c__0);
	}
	if (proceed) {
	    proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
		    usage_1.use[ic - 1] || usage_1.use[id - 1];
	}

/*     compute the value of the torsional angle */

	if (proceed) {
	    xia = atoms_1.x[ia - 1];
	    yia = atoms_1.y[ia - 1];
	    zia = atoms_1.z__[ia - 1];
	    xib = atoms_1.x[ib - 1];
	    yib = atoms_1.y[ib - 1];
	    zib = atoms_1.z__[ib - 1];
	    xic = atoms_1.x[ic - 1];
	    yic = atoms_1.y[ic - 1];
	    zic = atoms_1.z__[ic - 1];
	    xid = atoms_1.x[id - 1];
	    yid = atoms_1.y[id - 1];
	    zid = atoms_1.z__[id - 1];
	    xba = xib - xia;
	    yba = yib - yia;
	    zba = zib - zia;
	    xcb = xic - xib;
	    ycb = yic - yib;
	    zcb = zic - zib;
	    xdc = xid - xic;
	    ydc = yid - yic;
	    zdc = zid - zic;
	    xt = yba * zcb - ycb * zba;
	    yt = zba * xcb - zcb * xba;
	    zt = xba * ycb - xcb * yba;
	    xu = ycb * zdc - ydc * zcb;
	    yu = zcb * xdc - zdc * xcb;
	    zu = xcb * ydc - xdc * ycb;
	    xtu = yt * zu - yu * zt;
	    ytu = zt * xu - zu * xt;
	    ztu = xt * yu - xu * yt;
	    rt2 = xt * xt + yt * yt + zt * zt;
	    ru2 = xu * xu + yu * yu + zu * zu;
	    rtru = sqrt(rt2 * ru2);
	    if (rtru != 0.) {
		rcb = sqrt(xcb * xcb + ycb * ycb + zcb * zcb);
		cosine = (xt * xu + yt * yu + zt * zu) / rtru;
		sine = (xcb * xtu + ycb * ytu + zcb * ztu) / (rcb * rtru);

/*     set the torsional parameters for this angle */

		v1 = tors1_ref(1, i__);
		c1 = tors1_ref(3, i__);
		s1 = tors1_ref(4, i__);
		v2 = tors2_ref(1, i__);
		c2 = tors2_ref(3, i__);
		s2 = tors2_ref(4, i__);
		v3 = tors3_ref(1, i__);
		c3 = tors3_ref(3, i__);
		s3 = tors3_ref(4, i__);
		v4 = tors4_ref(1, i__);
		c4 = tors4_ref(3, i__);
		s4 = tors4_ref(4, i__);
		v5 = tors5_ref(1, i__);
		c5 = tors5_ref(3, i__);
		s5 = tors5_ref(4, i__);
		v6 = tors6_ref(1, i__);
		c6 = tors6_ref(3, i__);
		s6 = tors6_ref(4, i__);

/*     compute the multiple angle trigonometry and the phase terms */

		cosine2 = cosine * cosine - sine * sine;
		sine2 = cosine * 2. * sine;
		cosine3 = cosine * cosine2 - sine * sine2;
		sine3 = cosine * sine2 + sine * cosine2;
		cosine4 = cosine * cosine3 - sine * sine3;
		sine4 = cosine * sine3 + sine * cosine3;
		cosine5 = cosine * cosine4 - sine * sine4;
		sine5 = cosine * sine4 + sine * cosine4;
		cosine6 = cosine * cosine5 - sine * sine5;
		sine6 = cosine * sine5 + sine * cosine5;
		phi1 = cosine * c1 + sine * s1 + 1.;
		phi2 = cosine2 * c2 + sine2 * s2 + 1.;
		phi3 = cosine3 * c3 + sine3 * s3 + 1.;
		phi4 = cosine4 * c4 + sine4 * s4 + 1.;
		phi5 = cosine5 * c5 + sine5 * s5 + 1.;
		phi6 = cosine6 * c6 + sine6 * s6 + 1.;
		dphi1 = cosine * s1 - sine * c1;
		dphi2 = (cosine2 * s2 - sine2 * c2) * 2.;
		dphi3 = (cosine3 * s3 - sine3 * c3) * 3.;
		dphi4 = (cosine4 * s4 - sine4 * c4) * 4.;
		dphi5 = (cosine5 * s5 - sine5 * c5) * 5.;
		dphi6 = (cosine6 * s6 - sine6 * c6) * 6.;

/*     transform the potential function via smoothing */

		if (warp_1.use_gda__) {
		    width = wterm * (warp_1.m2[ia - 1] + warp_1.m2[ib - 1] + 
			    warp_1.m2[ic - 1] + warp_1.m2[id - 1]);
		    damp1 = exp(-width);
		    damp2 = exp(width * -4.);
		    damp3 = exp(width * -9.);
		    damp4 = exp(width * -16.);
		    damp5 = exp(width * -25.);
		    damp6 = exp(width * -36.);
		}
		phi1 *= damp1;
		phi2 *= damp2;
		phi3 *= damp3;
		phi4 *= damp4;
		phi5 *= damp5;
		phi6 *= damp6;
		dphi1 *= damp1;
		dphi2 *= damp2;
		dphi3 *= damp3;
		dphi4 *= damp4;
		dphi5 *= damp5;
		dphi6 *= damp6;

/*     calculate torsional energy and master chain rule term */

		e = torpot_1.torsunit * (v1 * phi1 + v2 * phi2 + v3 * phi3 + 
			v4 * phi4 + v5 * phi5 + v6 * phi6);
		dedphi = torpot_1.torsunit * (v1 * dphi1 + v2 * dphi2 + v3 * 
			dphi3 + v4 * dphi4 + v5 * dphi5 + v6 * dphi6);

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    e *= fgrp;
		    dedphi *= fgrp;
		}

/*     chain rule terms for first derivative components */

		xca = xic - xia;
		yca = yic - yia;
		zca = zic - zia;
		xdb = xid - xib;
		ydb = yid - yib;
		zdb = zid - zib;
		dedxt = dedphi * (yt * zcb - ycb * zt) / (rt2 * rcb);
		dedyt = dedphi * (zt * xcb - zcb * xt) / (rt2 * rcb);
		dedzt = dedphi * (xt * ycb - xcb * yt) / (rt2 * rcb);
		dedxu = -dedphi * (yu * zcb - ycb * zu) / (ru2 * rcb);
		dedyu = -dedphi * (zu * xcb - zcb * xu) / (ru2 * rcb);
		dedzu = -dedphi * (xu * ycb - xcb * yu) / (ru2 * rcb);

/*     compute first derivative components for this angle */

		dedxia = zcb * dedyt - ycb * dedzt;
		dedyia = xcb * dedzt - zcb * dedxt;
		dedzia = ycb * dedxt - xcb * dedyt;
		dedxib = yca * dedzt - zca * dedyt + zdc * dedyu - ydc * 
			dedzu;
		dedyib = zca * dedxt - xca * dedzt + xdc * dedzu - zdc * 
			dedxu;
		dedzib = xca * dedyt - yca * dedxt + ydc * dedxu - xdc * 
			dedyu;
		dedxic = zba * dedyt - yba * dedzt + ydb * dedzu - zdb * 
			dedyu;
		dedyic = xba * dedzt - zba * dedxt + zdb * dedxu - xdb * 
			dedzu;
		dedzic = yba * dedxt - xba * dedyt + xdb * dedyu - ydb * 
			dedxu;
		dedxid = zcb * dedyu - ycb * dedzu;
		dedyid = xcb * dedzu - zcb * dedxu;
		dedzid = ycb * dedxu - xcb * dedyu;

/*     increment the total torsional angle energy and gradient */

		energi_1.et += e;
		det_ref(1, ia) = det_ref(1, ia) + dedxia;
		det_ref(2, ia) = det_ref(2, ia) + dedyia;
		det_ref(3, ia) = det_ref(3, ia) + dedzia;
		det_ref(1, ib) = det_ref(1, ib) + dedxib;
		det_ref(2, ib) = det_ref(2, ib) + dedyib;
		det_ref(3, ib) = det_ref(3, ib) + dedzib;
		det_ref(1, ic) = det_ref(1, ic) + dedxic;
		det_ref(2, ic) = det_ref(2, ic) + dedyic;
		det_ref(3, ic) = det_ref(3, ic) + dedzic;
		det_ref(1, id) = det_ref(1, id) + dedxid;
		det_ref(2, id) = det_ref(2, id) + dedyid;
		det_ref(3, id) = det_ref(3, id) + dedzid;

/*     increment the internal virial tensor components */

		vxx = xcb * (dedxic + dedxid) - xba * dedxia + xdc * dedxid;
		vyx = ycb * (dedxic + dedxid) - yba * dedxia + ydc * dedxid;
		vzx = zcb * (dedxic + dedxid) - zba * dedxia + zdc * dedxid;
		vyy = ycb * (dedyic + dedyid) - yba * dedyia + ydc * dedyid;
		vzy = zcb * (dedyic + dedyid) - zba * dedyia + zdc * dedyid;
		vzz = zcb * (dedzic + dedzid) - zba * dedzia + zdc * dedzid;
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
	}
    }
    return 0;
} /* etors1b_ */

#undef itors_ref
#undef tors6_ref
#undef tors5_ref
#undef tors4_ref
#undef tors3_ref
#undef tors2_ref
#undef tors1_ref
#undef vir_ref
#undef det_ref


