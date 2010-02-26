/* egeom2.f -- translated by f2c (version 20050501).
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
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    doublereal hessx[75000]	/* was [3][25000] */, hessy[75000]	/* 
	    was [3][25000] */, hessz[75000]	/* was [3][25000] */;
} hessn_;

#define hessn_1 hessn_

struct {
    doublereal xpfix[25000], ypfix[25000], zpfix[25000], pfix[50000]	/* 
	    was [2][25000] */, dfix[75000]	/* was [3][25000] */, afix[
	    75000]	/* was [3][25000] */, tfix[75000]	/* was [3][
	    25000] */, gfix[75000]	/* was [3][25000] */, chir[75000]	
	    /* was [3][25000] */, depth, width, rwall;
    integer npfix, ipfix[25000], kpfix[75000]	/* was [3][25000] */, ndfix, 
	    idfix[50000]	/* was [2][25000] */, nafix, iafix[75000]	
	    /* was [3][25000] */, ntfix, itfix[100000]	/* was [4][25000] */, 
	    ngfix, igfix[50000]	/* was [2][25000] */, nchir, ichir[100000]	
	    /* was [4][25000] */;
    logical use_basin__, use_wall__;
} kgeoms_;

#define kgeoms_1 kgeoms_

struct {
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine egeom2  --  atom-by-atom restraint Hessian  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "egeom2" calculates second derivatives of restraints */
/*     on positions, distances, angles and torsions as well */
/*     as Gaussian basin and spherical droplet restraints */

/*     note that the Hessian is discontinuous when an upper and */
/*     lower bound range is used instead of a single distance */


/* Subroutine */ int egeom2_(integer *i__)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal a, b;
    static integer j, k, m;
    static doublereal r__, c1, c2, c3, t1, t2, r2, r6, de;
    static integer ia, ib, ic, id;
    static doublereal r12, dt, ri, xi, rp, yi, zi, xp, yp, xr, yr, zr, zp, xt,
	     yt, zt, xu, yu, zu, af1, af2, cf1, d2e[9]	/* was [3][3] */, df1,
	     df2, gf1, gf2, cf2, dt2, tf1, tf2, ri2, rp2, rt2, ru2, rcb, xab, 
	    xia, yia, zia, xib, yib, dot, zib, xic, yic, zic, xid, yid, zid, 
	    yab, zab, xba, yba, zba, xcb, ycb, zcb, xdc, ydc, zdc, xca, yca, 
	    zca, xdb, ydb, zdb, xad, yad, zad, xbd, ybd, zbd, xcd, ycd, zcd, 
	    rab2, rcb2, xpo, ypo, zpo, xtu, ytu, ztu, xcm, ycm, zcm, vol;
    static integer kang;
    static doublereal dedr, xabp, fgrp, sine, term, xrab;
    static integer kpos;
    static doublereal yrab, zrab, xrcb, yrcb, zrcb, yabp, zabp, xcbp, ycbp, 
	    zcbp, rtru, rcbt2, rcbu2, xycb2, xzcb2, yzcb2, d2edr2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal deddt, angle, force;
    static integer kchir;
    static doublereal terma, termc;
    static integer kdist;
    static doublereal weigh, ratio, rcbxt, rcbyt, termx;
    static integer ktors;
    static doublereal termy, termz, rcbzt, rcbxu, rcbyu, rcbzu, d2eddt2, 
	    dedphi, weigha, weighb, buffer;
    static logical linear;
    static doublereal cosine, target;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal d2edphi2, ddtdxia, ddtdyia, ddtdzia, ddtdxib, ddtdyib, 
	    ddtdzib, ddtdxic, ddtdyic, ddtdzic, dxiaxia, dyiayia, dziazia, 
	    dxibxib, dyibyib, dzibzib, dxicxic, dyicyic, dphidxt, dphidyt, 
	    dphidzt, dphidxu, dphidyu, dphidzu, dziczic, dxidxid, dyidyid, 
	    dzidzid, dxiayia, dxiazia, dyiazia, dxibyib, dxibzib, dyibzib, 
	    dxicyic, dxiczic, dyiczic, dxidyid, dxidzid, dyidzid, dxiaxib, 
	    dxiayib, dxiazib, dyiaxib, dyiayib, dyiazib, dziaxib, dziayib, 
	    dziazib, dxiaxic, dxiayic, dxiazic, dyiaxic, dyiayic, dyiazic, 
	    dziaxic, dziayic, dziazic, dxiaxid, dxiayid, dxiazid, dyiaxid, 
	    dyiayid, dyiazid, dziaxid, dziayid, dziazid, dxibxia, dxibyia, 
	    dxibzia, dyibxia, dyibyia, dyibzia, dzibxia, dzibyia, dzibzia, 
	    dxibxic, dxibyic, dxibzic, dyibxic, dyibyic, dyibzic, dzibxic, 
	    dzibyic, dzibzic, dxibxid, dxibyid, dxibzid, dyibxid, dyibyid, 
	    dyibzid, dzibxid, dphidxia, dphidyia, dphidzia, dphidxib, 
	    dphidyib, dphidzib, dphidxic, dphidyic, dphidzic, dphidxid, 
	    dphidyid, dphidzid, dzibyid, dzibzid, dxicxid, dxicyid, dxiczid, 
	    dyicxid, dyicyid, dyiczid, dzicxid, dzicyid, dziczid, ddtdxid, 
	    ddtdyid, ddtdzid, expterm;
    static logical proceed, intermol;
    static doublereal dphidxibt, dphidyibt, dphidzibt, dphidxibu, dphidyibu, 
	    dphidzibu, dphidxict, dphidyict, dphidzict, dphidxicu, dphidyicu, 
	    dphidzicu;


#define d2e_ref(a_1,a_2) d2e[(a_2)*3 + a_1 - 4]
#define chir_ref(a_1,a_2) kgeoms_1.chir[(a_2)*3 + a_1 - 4]
#define afix_ref(a_1,a_2) kgeoms_1.afix[(a_2)*3 + a_1 - 4]
#define dfix_ref(a_1,a_2) kgeoms_1.dfix[(a_2)*3 + a_1 - 4]
#define gfix_ref(a_1,a_2) kgeoms_1.gfix[(a_2)*3 + a_1 - 4]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]
#define pfix_ref(a_1,a_2) kgeoms_1.pfix[(a_2)*2 + a_1 - 3]
#define tfix_ref(a_1,a_2) kgeoms_1.tfix[(a_2)*3 + a_1 - 4]
#define ichir_ref(a_1,a_2) kgeoms_1.ichir[(a_2)*4 + a_1 - 5]
#define iafix_ref(a_1,a_2) kgeoms_1.iafix[(a_2)*3 + a_1 - 4]
#define idfix_ref(a_1,a_2) kgeoms_1.idfix[(a_2)*2 + a_1 - 3]
#define igfix_ref(a_1,a_2) kgeoms_1.igfix[(a_2)*2 + a_1 - 3]
#define kpfix_ref(a_1,a_2) kgeoms_1.kpfix[(a_2)*3 + a_1 - 4]
#define itfix_ref(a_1,a_2) kgeoms_1.itfix[(a_2)*4 + a_1 - 5]
#define hessx_ref(a_1,a_2) hessn_1.hessx[(a_2)*3 + a_1 - 4]
#define hessy_ref(a_1,a_2) hessn_1.hessy[(a_2)*3 + a_1 - 4]
#define hessz_ref(a_1,a_2) hessn_1.hessz[(a_2)*3 + a_1 - 4]



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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  hessn.i  --  Cartesian Hessian elements for a single atom  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     hessx   Hessian elements for x-component of current atom */
/*     hessy   Hessian elements for y-component of current atom */
/*     hessz   Hessian elements for z-component of current atom */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  kgeoms.i  --  parameters for the geometrical restraints  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xpfix      x-coordinate target for each restrained position */
/*     ypfix      y-coordinate target for each restrained position */
/*     zpfix      z-coordinate target for each restrained position */
/*     pfix       force constant and flat-well range for each position */
/*     dfix       force constant and target range for each distance */
/*     afix       force constant and target range for each angle */
/*     tfix       force constant and target range for each torsion */
/*     gfix       force constant and target range for each group distance */
/*     chir       force constant and target range for chiral centers */
/*     depth      depth of shallow Gaussian basin restraint */
/*     width      exponential width coefficient of Gaussian basin */
/*     rwall      radius of spherical droplet boundary restraint */
/*     npfix      number of position restraints to be applied */
/*     ipfix      atom number involved in each position restraint */
/*     kpfix      flags to use x-, y-, z-coordinate position restraints */
/*     ndfix      number of distance restraints to be applied */
/*     idfix      atom numbers defining each distance restraint */
/*     nafix      number of angle restraints to be applied */
/*     iafix      atom numbers defining each angle restraint */
/*     ntfix      number of torsional restraints to be applied */
/*     itfix      atom numbers defining each torsional restraint */
/*     ngfix      number of group distance restraints to be applied */
/*     igfix      group numbers defining each group distance restraint */
/*     nchir      number of chirality restraints to be applied */
/*     ichir      atom numbers defining each chirality restraint */
/*     use_basin  logical flag governing use of Gaussian basin */
/*     use_wall   logical flag governing use of droplet boundary */




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




/*     compute the Hessian elements for position restraints */

    i__1 = kgeoms_1.npfix;
    for (kpos = 1; kpos <= i__1; ++kpos) {
	ia = kgeoms_1.ipfix[kpos - 1];
	proceed = *i__ == ia;
	if (proceed && group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &c__0, &c__0, &c__0, &c__0, &c__0);
	}
	if (proceed) {
	    xr = 0.;
	    yr = 0.;
	    zr = 0.;
	    if (kpfix_ref(1, *i__) != 0) {
		xr = atoms_1.x[ia - 1] - kgeoms_1.xpfix[kpos - 1];
	    }
	    if (kpfix_ref(2, *i__) != 0) {
		yr = atoms_1.y[ia - 1] - kgeoms_1.ypfix[kpos - 1];
	    }
	    if (kpfix_ref(3, *i__) != 0) {
		zr = atoms_1.z__[ia - 1] - kgeoms_1.zpfix[kpos - 1];
	    }
	    r2 = xr * xr + yr * yr + zr * zr;
	    r__ = sqrt(r2);
	    force = pfix_ref(1, kpos);
/* Computing MAX */
	    d__1 = 0., d__2 = r__ - pfix_ref(2, kpos);
	    dt = max(d__1,d__2);
	    dt2 = dt * dt;
	    deddt = force * 2.;

/*     scale the interaction based on its group membership */

	    if (group_1.use_group__) {
		deddt *= fgrp;
	    }

/*     set the chain rule terms for the Hessian elements */

	    if (r__ == 0.) {
		de = deddt;
		term = 0.;
	    } else {
		de = deddt * dt / r__;
		term = (deddt - de) / r2;
	    }
	    termx = term * xr;
	    termy = term * yr;
	    termz = term * zr;
	    d2e_ref(1, 1) = termx * xr + de;
	    d2e_ref(1, 2) = termx * yr;
	    d2e_ref(1, 3) = termx * zr;
	    d2e_ref(2, 1) = d2e_ref(1, 2);
	    d2e_ref(2, 2) = termy * yr + de;
	    d2e_ref(2, 3) = termy * zr;
	    d2e_ref(3, 1) = d2e_ref(1, 3);
	    d2e_ref(3, 2) = d2e_ref(2, 3);
	    d2e_ref(3, 3) = termz * zr + de;

/*     increment diagonal and non-diagonal Hessian elements */

	    for (j = 1; j <= 3; ++j) {
		hessx_ref(j, ia) = hessx_ref(j, ia) + d2e_ref(1, j);
		hessy_ref(j, ia) = hessy_ref(j, ia) + d2e_ref(2, j);
		hessz_ref(j, ia) = hessz_ref(j, ia) + d2e_ref(3, j);
	    }
	}
    }

/*     compute the Hessian elements for distance restraints */

    i__1 = kgeoms_1.ndfix;
    for (kdist = 1; kdist <= i__1; ++kdist) {
	ia = idfix_ref(1, kdist);
	ib = idfix_ref(2, kdist);
	proceed = *i__ == ia || *i__ == ib;
	if (proceed && group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &c__0, &c__0, &c__0, &c__0);
	}
	if (proceed) {
	    if (*i__ == ib) {
		ib = ia;
		ia = *i__;
	    }
	    xr = atoms_1.x[ia - 1] - atoms_1.x[ib - 1];
	    yr = atoms_1.y[ia - 1] - atoms_1.y[ib - 1];
	    zr = atoms_1.z__[ia - 1] - atoms_1.z__[ib - 1];
	    intermol = molcul_1.molcule[ia - 1] != molcul_1.molcule[ib - 1];
	    if (bound_1.use_bounds__ && intermol) {
		image_(&xr, &yr, &zr);
	    }
	    r2 = xr * xr + yr * yr + zr * zr;
	    r__ = sqrt(r2);
	    force = dfix_ref(1, kdist);
	    df1 = dfix_ref(2, kdist);
	    df2 = dfix_ref(3, kdist);
	    target = r__;
	    if (r__ < df1) {
		target = df1;
	    }
	    if (r__ > df2) {
		target = df2;
	    }
	    dt = r__ - target;
	    deddt = force * 2.;

/*     scale the interaction based on its group membership */

	    if (group_1.use_group__) {
		deddt *= fgrp;
	    }

/*     set the chain rule terms for the Hessian elements */

	    if (r__ == 0.) {
		r__ = 1e-4;
		r2 = r__ * r__;
	    }
	    de = deddt * dt / r__;
	    term = (deddt - de) / r2;
	    termx = term * xr;
	    termy = term * yr;
	    termz = term * zr;
	    d2e_ref(1, 1) = termx * xr + de;
	    d2e_ref(1, 2) = termx * yr;
	    d2e_ref(1, 3) = termx * zr;
	    d2e_ref(2, 1) = d2e_ref(1, 2);
	    d2e_ref(2, 2) = termy * yr + de;
	    d2e_ref(2, 3) = termy * zr;
	    d2e_ref(3, 1) = d2e_ref(1, 3);
	    d2e_ref(3, 2) = d2e_ref(2, 3);
	    d2e_ref(3, 3) = termz * zr + de;

/*     increment diagonal and non-diagonal Hessian elements */

	    for (j = 1; j <= 3; ++j) {
		hessx_ref(j, ia) = hessx_ref(j, ia) + d2e_ref(1, j);
		hessy_ref(j, ia) = hessy_ref(j, ia) + d2e_ref(2, j);
		hessz_ref(j, ia) = hessz_ref(j, ia) + d2e_ref(3, j);
		hessx_ref(j, ib) = hessx_ref(j, ib) - d2e_ref(1, j);
		hessy_ref(j, ib) = hessy_ref(j, ib) - d2e_ref(2, j);
		hessz_ref(j, ib) = hessz_ref(j, ib) - d2e_ref(3, j);
	    }
	}
    }

/*     compute the Hessian elements for angle restraints */

    i__1 = kgeoms_1.nafix;
    for (kang = 1; kang <= i__1; ++kang) {
	ia = iafix_ref(1, kang);
	ib = iafix_ref(2, kang);
	ic = iafix_ref(3, kang);
	proceed = *i__ == ia || *i__ == ib || *i__ == ic;
	if (proceed && group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &c__0, &c__0, &c__0);
	}
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
	    xab = xia - xib;
	    yab = yia - yib;
	    zab = zia - zib;
	    xcb = xic - xib;
	    ycb = yic - yib;
	    zcb = zic - zib;
	    rab2 = xab * xab + yab * yab + zab * zab;
	    rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
	    if (rab2 != 0. && rcb2 != 0.) {
		xp = ycb * zab - zcb * yab;
		yp = zcb * xab - xcb * zab;
		zp = xcb * yab - ycb * xab;
		rp = sqrt(xp * xp + yp * yp + zp * zp);
		dot = xab * xcb + yab * ycb + zab * zcb;
		cosine = dot / sqrt(rab2 * rcb2);
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine);
		cosine = min(d__1,d__2);
		angle = acos(cosine) * 57.29577951308232088;
		force = afix_ref(1, kang);
		af1 = afix_ref(2, kang);
		af2 = afix_ref(3, kang);
		target = angle;
		if (angle < af1) {
		    target = af1;
		}
		if (angle > af2) {
		    target = af2;
		}
		dt = angle - target;
		dt2 = dt * dt;
		deddt = force * 2. * dt * 57.29577951308232088;
		d2eddt2 = force * 2. * 3282.8063500117441;

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    deddt *= fgrp;
		    d2eddt2 *= fgrp;
		}

/*     construct an orthogonal direction for linear angles */

		linear = FALSE_;
		if (rp < 1e-6) {
		    linear = TRUE_;
		    if (xab != 0. && yab != 0.) {
			xp = -yab;
			yp = xab;
			zp = 0.;
		    } else if (xab == 0. && yab == 0.) {
			xp = 1.;
			yp = 0.;
			zp = 0.;
		    } else if (xab != 0. && yab == 0.) {
			xp = 0.;
			yp = 1.;
			zp = 0.;
		    } else if (xab == 0. && yab != 0.) {
			xp = 1.;
			yp = 0.;
			zp = 0.;
		    }
		    rp = sqrt(xp * xp + yp * yp + zp * zp);
		}

/*     first derivatives of bond angle with respect to coordinates */

L10:
		terma = -1. / (rab2 * rp);
		termc = 1. / (rcb2 * rp);
		ddtdxia = terma * (yab * zp - zab * yp);
		ddtdyia = terma * (zab * xp - xab * zp);
		ddtdzia = terma * (xab * yp - yab * xp);
		ddtdxic = termc * (ycb * zp - zcb * yp);
		ddtdyic = termc * (zcb * xp - xcb * zp);
		ddtdzic = termc * (xcb * yp - ycb * xp);
		ddtdxib = -ddtdxia - ddtdxic;
		ddtdyib = -ddtdyia - ddtdyic;
		ddtdzib = -ddtdzia - ddtdzic;

/*     abbreviations used in defining chain rule terms */

		xrab = xab * 2. / rab2;
		yrab = yab * 2. / rab2;
		zrab = zab * 2. / rab2;
		xrcb = xcb * 2. / rcb2;
		yrcb = ycb * 2. / rcb2;
		zrcb = zcb * 2. / rcb2;
		rp2 = 1. / (rp * rp);
		xabp = (yab * zp - zab * yp) * rp2;
		yabp = (zab * xp - xab * zp) * rp2;
		zabp = (xab * yp - yab * xp) * rp2;
		xcbp = (ycb * zp - zcb * yp) * rp2;
		ycbp = (zcb * xp - xcb * zp) * rp2;
		zcbp = (xcb * yp - ycb * xp) * rp2;

/*     chain rule terms for second derivative components */

		dxiaxia = terma * (xab * xcb - dot) + ddtdxia * (xcbp - xrab);
		dxiayia = terma * (zp + yab * xcb) + ddtdxia * (ycbp - yrab);
		dxiazia = terma * (zab * xcb - yp) + ddtdxia * (zcbp - zrab);
		dyiayia = terma * (yab * ycb - dot) + ddtdyia * (ycbp - yrab);
		dyiazia = terma * (xp + zab * ycb) + ddtdyia * (zcbp - zrab);
		dziazia = terma * (zab * zcb - dot) + ddtdzia * (zcbp - zrab);
		dxicxic = termc * (dot - xab * xcb) - ddtdxic * (xabp + xrcb);
		dxicyic = termc * (zp - ycb * xab) - ddtdxic * (yabp + yrcb);
		dxiczic = -termc * (yp + zcb * xab) - ddtdxic * (zabp + zrcb);
		dyicyic = termc * (dot - yab * ycb) - ddtdyic * (yabp + yrcb);
		dyiczic = termc * (xp - zcb * yab) - ddtdyic * (zabp + zrcb);
		dziczic = termc * (dot - zab * zcb) - ddtdzic * (zabp + zrcb);
		dxiaxic = terma * (yab * yab + zab * zab) - ddtdxia * xabp;
		dxiayic = -terma * xab * yab - ddtdxia * yabp;
		dxiazic = -terma * xab * zab - ddtdxia * zabp;
		dyiaxic = -terma * xab * yab - ddtdyia * xabp;
		dyiayic = terma * (xab * xab + zab * zab) - ddtdyia * yabp;
		dyiazic = -terma * yab * zab - ddtdyia * zabp;
		dziaxic = -terma * xab * zab - ddtdzia * xabp;
		dziayic = -terma * yab * zab - ddtdzia * yabp;
		dziazic = terma * (xab * xab + yab * yab) - ddtdzia * zabp;

/*     get some second derivative chain rule terms by difference */

		dxibxia = -dxiaxia - dxiaxic;
		dxibyia = -dxiayia - dyiaxic;
		dxibzia = -dxiazia - dziaxic;
		dyibxia = -dxiayia - dxiayic;
		dyibyia = -dyiayia - dyiayic;
		dyibzia = -dyiazia - dziayic;
		dzibxia = -dxiazia - dxiazic;
		dzibyia = -dyiazia - dyiazic;
		dzibzia = -dziazia - dziazic;
		dxibxic = -dxicxic - dxiaxic;
		dxibyic = -dxicyic - dxiayic;
		dxibzic = -dxiczic - dxiazic;
		dyibxic = -dxicyic - dyiaxic;
		dyibyic = -dyicyic - dyiayic;
		dyibzic = -dyiczic - dyiazic;
		dzibxic = -dxiczic - dziaxic;
		dzibyic = -dyiczic - dziayic;
		dzibzic = -dziczic - dziazic;
		dxibxib = -dxibxia - dxibxic;
		dxibyib = -dxibyia - dxibyic;
		dxibzib = -dxibzia - dxibzic;
		dyibyib = -dyibyia - dyibyic;
		dyibzib = -dyibzia - dyibzic;
		dzibzib = -dzibzia - dzibzic;

/*     increment diagonal and off-diagonal Hessian elements */

		if (ia == *i__) {
		    hessx_ref(1, ia) = hessx_ref(1, ia) + deddt * dxiaxia + 
			    d2eddt2 * ddtdxia * ddtdxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + deddt * dxiayia + 
			    d2eddt2 * ddtdxia * ddtdyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + deddt * dxiazia + 
			    d2eddt2 * ddtdxia * ddtdzia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + deddt * dxiayia + 
			    d2eddt2 * ddtdyia * ddtdxia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + deddt * dyiayia + 
			    d2eddt2 * ddtdyia * ddtdyia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + deddt * dyiazia + 
			    d2eddt2 * ddtdyia * ddtdzia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + deddt * dxiazia + 
			    d2eddt2 * ddtdzia * ddtdxia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + deddt * dyiazia + 
			    d2eddt2 * ddtdzia * ddtdyia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + deddt * dziazia + 
			    d2eddt2 * ddtdzia * ddtdzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + deddt * dxibxia + 
			    d2eddt2 * ddtdxia * ddtdxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + deddt * dyibxia + 
			    d2eddt2 * ddtdxia * ddtdyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + deddt * dzibxia + 
			    d2eddt2 * ddtdxia * ddtdzib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + deddt * dxibyia + 
			    d2eddt2 * ddtdyia * ddtdxib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + deddt * dyibyia + 
			    d2eddt2 * ddtdyia * ddtdyib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + deddt * dzibyia + 
			    d2eddt2 * ddtdyia * ddtdzib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + deddt * dxibzia + 
			    d2eddt2 * ddtdzia * ddtdxib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + deddt * dyibzia + 
			    d2eddt2 * ddtdzia * ddtdyib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + deddt * dzibzia + 
			    d2eddt2 * ddtdzia * ddtdzib;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + deddt * dxiaxic + 
			    d2eddt2 * ddtdxia * ddtdxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + deddt * dxiayic + 
			    d2eddt2 * ddtdxia * ddtdyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + deddt * dxiazic + 
			    d2eddt2 * ddtdxia * ddtdzic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + deddt * dyiaxic + 
			    d2eddt2 * ddtdyia * ddtdxic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + deddt * dyiayic + 
			    d2eddt2 * ddtdyia * ddtdyic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + deddt * dyiazic + 
			    d2eddt2 * ddtdyia * ddtdzic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + deddt * dziaxic + 
			    d2eddt2 * ddtdzia * ddtdxic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + deddt * dziayic + 
			    d2eddt2 * ddtdzia * ddtdyic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + deddt * dziazic + 
			    d2eddt2 * ddtdzia * ddtdzic;
		} else if (ib == *i__) {
		    hessx_ref(1, ib) = hessx_ref(1, ib) + deddt * dxibxib + 
			    d2eddt2 * ddtdxib * ddtdxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + deddt * dxibyib + 
			    d2eddt2 * ddtdxib * ddtdyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + deddt * dxibzib + 
			    d2eddt2 * ddtdxib * ddtdzib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + deddt * dxibyib + 
			    d2eddt2 * ddtdyib * ddtdxib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + deddt * dyibyib + 
			    d2eddt2 * ddtdyib * ddtdyib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + deddt * dyibzib + 
			    d2eddt2 * ddtdyib * ddtdzib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + deddt * dxibzib + 
			    d2eddt2 * ddtdzib * ddtdxib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + deddt * dyibzib + 
			    d2eddt2 * ddtdzib * ddtdyib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + deddt * dzibzib + 
			    d2eddt2 * ddtdzib * ddtdzib;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + deddt * dxibxia + 
			    d2eddt2 * ddtdxib * ddtdxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + deddt * dxibyia + 
			    d2eddt2 * ddtdxib * ddtdyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + deddt * dxibzia + 
			    d2eddt2 * ddtdxib * ddtdzia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + deddt * dyibxia + 
			    d2eddt2 * ddtdyib * ddtdxia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + deddt * dyibyia + 
			    d2eddt2 * ddtdyib * ddtdyia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + deddt * dyibzia + 
			    d2eddt2 * ddtdyib * ddtdzia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + deddt * dzibxia + 
			    d2eddt2 * ddtdzib * ddtdxia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + deddt * dzibyia + 
			    d2eddt2 * ddtdzib * ddtdyia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + deddt * dzibzia + 
			    d2eddt2 * ddtdzib * ddtdzia;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + deddt * dxibxic + 
			    d2eddt2 * ddtdxib * ddtdxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + deddt * dxibyic + 
			    d2eddt2 * ddtdxib * ddtdyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + deddt * dxibzic + 
			    d2eddt2 * ddtdxib * ddtdzic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + deddt * dyibxic + 
			    d2eddt2 * ddtdyib * ddtdxic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + deddt * dyibyic + 
			    d2eddt2 * ddtdyib * ddtdyic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + deddt * dyibzic + 
			    d2eddt2 * ddtdyib * ddtdzic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + deddt * dzibxic + 
			    d2eddt2 * ddtdzib * ddtdxic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + deddt * dzibyic + 
			    d2eddt2 * ddtdzib * ddtdyic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + deddt * dzibzic + 
			    d2eddt2 * ddtdzib * ddtdzic;
		} else if (ic == *i__) {
		    hessx_ref(1, ic) = hessx_ref(1, ic) + deddt * dxicxic + 
			    d2eddt2 * ddtdxic * ddtdxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + deddt * dxicyic + 
			    d2eddt2 * ddtdxic * ddtdyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + deddt * dxiczic + 
			    d2eddt2 * ddtdxic * ddtdzic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + deddt * dxicyic + 
			    d2eddt2 * ddtdyic * ddtdxic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + deddt * dyicyic + 
			    d2eddt2 * ddtdyic * ddtdyic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + deddt * dyiczic + 
			    d2eddt2 * ddtdyic * ddtdzic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + deddt * dxiczic + 
			    d2eddt2 * ddtdzic * ddtdxic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + deddt * dyiczic + 
			    d2eddt2 * ddtdzic * ddtdyic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + deddt * dziczic + 
			    d2eddt2 * ddtdzic * ddtdzic;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + deddt * dxibxic + 
			    d2eddt2 * ddtdxic * ddtdxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + deddt * dyibxic + 
			    d2eddt2 * ddtdxic * ddtdyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + deddt * dzibxic + 
			    d2eddt2 * ddtdxic * ddtdzib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + deddt * dxibyic + 
			    d2eddt2 * ddtdyic * ddtdxib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + deddt * dyibyic + 
			    d2eddt2 * ddtdyic * ddtdyib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + deddt * dzibyic + 
			    d2eddt2 * ddtdyic * ddtdzib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + deddt * dxibzic + 
			    d2eddt2 * ddtdzic * ddtdxib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + deddt * dyibzic + 
			    d2eddt2 * ddtdzic * ddtdyib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + deddt * dzibzic + 
			    d2eddt2 * ddtdzic * ddtdzib;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + deddt * dxiaxic + 
			    d2eddt2 * ddtdxic * ddtdxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + deddt * dyiaxic + 
			    d2eddt2 * ddtdxic * ddtdyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + deddt * dziaxic + 
			    d2eddt2 * ddtdxic * ddtdzia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + deddt * dxiayic + 
			    d2eddt2 * ddtdyic * ddtdxia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + deddt * dyiayic + 
			    d2eddt2 * ddtdyic * ddtdyia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + deddt * dziayic + 
			    d2eddt2 * ddtdyic * ddtdzia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + deddt * dxiazic + 
			    d2eddt2 * ddtdzic * ddtdxia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + deddt * dyiazic + 
			    d2eddt2 * ddtdzic * ddtdyia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + deddt * dziazic + 
			    d2eddt2 * ddtdzic * ddtdzia;
		}

/*     construct a second orthogonal direction for linear angles */

		if (linear) {
		    linear = FALSE_;
		    xpo = xp;
		    ypo = yp;
		    zpo = zp;
		    xp = ypo * zab - zpo * yab;
		    yp = zpo * xab - xpo * zab;
		    zp = xpo * yab - ypo * xab;
		    rp = sqrt(xp * xp + yp * yp + zp * zp);
		    goto L10;
		}
	    }
	}
    }

/*     compute the Hessian elements for torsion restraints */

    i__1 = kgeoms_1.ntfix;
    for (ktors = 1; ktors <= i__1; ++ktors) {
	ia = itfix_ref(1, ktors);
	ib = itfix_ref(2, ktors);
	ic = itfix_ref(3, ktors);
	id = itfix_ref(4, ktors);
	proceed = *i__ == ia || *i__ == ib || *i__ == ic || *i__ == id;
	if (proceed && group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &c__0, &c__0);
	}
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
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine);
		cosine = min(d__1,d__2);
		angle = acos(cosine) * 57.29577951308232088;
		if (sine < 0.) {
		    angle = -angle;
		}

/*     calculate the pseudoenergy master chain rule terms */

		force = tfix_ref(1, ktors);
		tf1 = tfix_ref(2, ktors);
		tf2 = tfix_ref(3, ktors);
		if (angle > tf1 && angle < tf2) {
		    target = angle;
		} else if (angle > tf1 && tf1 > tf2) {
		    target = angle;
		} else if (angle < tf2 && tf1 > tf2) {
		    target = angle;
		} else {
		    t1 = angle - tf1;
		    t2 = angle - tf2;
		    if (t1 > 180.) {
			t1 += -360.;
		    } else if (t1 < -180.) {
			t1 += 360.;
		    }
		    if (t2 > 180.) {
			t2 += -360.;
		    } else if (t2 < -180.) {
			t2 += 360.;
		    }
		    if (abs(t1) < abs(t2)) {
			target = tf1;
		    } else {
			target = tf2;
		    }
		}
		dt = angle - target;
		if (dt > 180.) {
		    dt += -360.;
		} else if (dt < -180.) {
		    dt += 360.;
		}
		dedphi = force * 114.59155902616465 * dt;
		d2edphi2 = force * 6565.6127000234883;

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    dedphi *= fgrp;
		    d2edphi2 *= fgrp;
		}

/*     abbreviations for first derivative chain rule terms */

		xca = xic - xia;
		yca = yic - yia;
		zca = zic - zia;
		xdb = xid - xib;
		ydb = yid - yib;
		zdb = zid - zib;
		dphidxt = (yt * zcb - ycb * zt) / (rt2 * rcb);
		dphidyt = (zt * xcb - zcb * xt) / (rt2 * rcb);
		dphidzt = (xt * ycb - xcb * yt) / (rt2 * rcb);
		dphidxu = -(yu * zcb - ycb * zu) / (ru2 * rcb);
		dphidyu = -(zu * xcb - zcb * xu) / (ru2 * rcb);
		dphidzu = -(xu * ycb - xcb * yu) / (ru2 * rcb);

/*     abbreviations for second derivative chain rule terms */

		xycb2 = xcb * xcb + ycb * ycb;
		xzcb2 = xcb * xcb + zcb * zcb;
		yzcb2 = ycb * ycb + zcb * zcb;
		rcbxt = rcb * -2. * dphidxt;
		rcbyt = rcb * -2. * dphidyt;
		rcbzt = rcb * -2. * dphidzt;
		rcbt2 = rcb * rt2;
		rcbxu = rcb * 2. * dphidxu;
		rcbyu = rcb * 2. * dphidyu;
		rcbzu = rcb * 2. * dphidzu;
		rcbu2 = rcb * ru2;
		dphidxibt = yca * dphidzt - zca * dphidyt;
		dphidxibu = zdc * dphidyu - ydc * dphidzu;
		dphidyibt = zca * dphidxt - xca * dphidzt;
		dphidyibu = xdc * dphidzu - zdc * dphidxu;
		dphidzibt = xca * dphidyt - yca * dphidxt;
		dphidzibu = ydc * dphidxu - xdc * dphidyu;
		dphidxict = zba * dphidyt - yba * dphidzt;
		dphidxicu = ydb * dphidzu - zdb * dphidyu;
		dphidyict = xba * dphidzt - zba * dphidxt;
		dphidyicu = zdb * dphidxu - xdb * dphidzu;
		dphidzict = yba * dphidxt - xba * dphidyt;
		dphidzicu = xdb * dphidyu - ydb * dphidxu;

/*     chain rule terms for first derivative components */

		dphidxia = zcb * dphidyt - ycb * dphidzt;
		dphidyia = xcb * dphidzt - zcb * dphidxt;
		dphidzia = ycb * dphidxt - xcb * dphidyt;
		dphidxib = dphidxibt + dphidxibu;
		dphidyib = dphidyibt + dphidyibu;
		dphidzib = dphidzibt + dphidzibu;
		dphidxic = dphidxict + dphidxicu;
		dphidyic = dphidyict + dphidyicu;
		dphidzic = dphidzict + dphidzicu;
		dphidxid = zcb * dphidyu - ycb * dphidzu;
		dphidyid = xcb * dphidzu - zcb * dphidxu;
		dphidzid = ycb * dphidxu - xcb * dphidyu;

/*     chain rule terms for second derivative components */

		dxiaxia = rcbxt * dphidxia;
		dxiayia = rcbxt * dphidyia - zcb * rcb / rt2;
		dxiazia = rcbxt * dphidzia + ycb * rcb / rt2;
		dxiaxic = rcbxt * dphidxict + xcb * xt / rcbt2;
		dxiayic = rcbxt * dphidyict - dphidzt - (xba * zcb * xcb + 
			zba * yzcb2) / rcbt2;
		dxiazic = rcbxt * dphidzict + dphidyt + (xba * ycb * xcb + 
			yba * yzcb2) / rcbt2;
		dxiaxid = 0.;
		dxiayid = 0.;
		dxiazid = 0.;
		dyiayia = rcbyt * dphidyia;
		dyiazia = rcbyt * dphidzia - xcb * rcb / rt2;
		dyiaxib = rcbyt * dphidxibt - dphidzt - (yca * zcb * ycb + 
			zca * xzcb2) / rcbt2;
		dyiaxic = rcbyt * dphidxict + dphidzt + (yba * zcb * ycb + 
			zba * xzcb2) / rcbt2;
		dyiayic = rcbyt * dphidyict + ycb * yt / rcbt2;
		dyiazic = rcbyt * dphidzict - dphidxt - (yba * xcb * ycb + 
			xba * xzcb2) / rcbt2;
		dyiaxid = 0.;
		dyiayid = 0.;
		dyiazid = 0.;
		dziazia = rcbzt * dphidzia;
		dziaxib = rcbzt * dphidxibt + dphidyt + (zca * ycb * zcb + 
			yca * xycb2) / rcbt2;
		dziayib = rcbzt * dphidyibt - dphidxt - (zca * xcb * zcb + 
			xca * xycb2) / rcbt2;
		dziaxic = rcbzt * dphidxict - dphidyt - (zba * ycb * zcb + 
			yba * xycb2) / rcbt2;
		dziayic = rcbzt * dphidyict + dphidxt + (zba * xcb * zcb + 
			xba * xycb2) / rcbt2;
		dziazic = rcbzt * dphidzict + zcb * zt / rcbt2;
		dziaxid = 0.;
		dziayid = 0.;
		dziazid = 0.;
		dxibxic = -xcb * dphidxib / (rcb * rcb) - (yca * (zba * xcb + 
			yt) - zca * (yba * xcb - zt)) / rcbt2 - (yt * zba - 
			yba * zt) * 2. * dphidxibt / rt2 - (zdc * (ydb * xcb 
			+ zu) - ydc * (zdb * xcb - yu)) / rcbu2 + (yu * zdb - 
			ydb * zu) * 2. * dphidxibu / ru2;
		dxibyic = -ycb * dphidxib / (rcb * rcb) + dphidzt + dphidzu - 
			(yca * (zba * ycb - xt) + zca * (xba * xcb + zcb * 
			zba)) / rcbt2 - (zt * xba - zba * xt) * 2. * 
			dphidxibt / rt2 + (zdc * (xdb * xcb + zcb * zdb) + 
			ydc * (zdb * ycb + xu)) / rcbu2 + (zu * xdb - zdb * 
			xu) * 2. * dphidxibu / ru2;
		dxibxid = rcbxu * dphidxibu + xcb * xu / rcbu2;
		dxibyid = rcbyu * dphidxibu - dphidzu - (ydc * zcb * ycb + 
			zdc * xzcb2) / rcbu2;
		dxibzid = rcbzu * dphidxibu + dphidyu + (zdc * ycb * zcb + 
			ydc * xycb2) / rcbu2;
		dyibzib = ycb * dphidzib / (rcb * rcb) - (xca * (xca * xcb + 
			zcb * zca) + yca * (ycb * xca + zt)) / rcbt2 - (xt * 
			zca - xca * zt) * 2. * dphidzibt / rt2 + (ydc * (xdc *
			 ycb - zu) + xdc * (xdc * xcb + zcb * zdc)) / rcbu2 + 
			(xu * zdc - xdc * zu) * 2. * dphidzibu / ru2;
		dyibxic = -xcb * dphidyib / (rcb * rcb) - dphidzt - dphidzu + 
			(xca * (zba * xcb + yt) + zca * (zba * zcb + ycb * 
			yba)) / rcbt2 - (yt * zba - yba * zt) * 2. * 
			dphidyibt / rt2 - (zdc * (zdb * zcb + ycb * ydb) + 
			xdc * (zdb * xcb - yu)) / rcbu2 + (yu * zdb - ydb * 
			zu) * 2. * dphidyibu / ru2;
		dyibyic = -ycb * dphidyib / (rcb * rcb) - (zca * (xba * ycb + 
			zt) - xca * (zba * ycb - xt)) / rcbt2 - (zt * xba - 
			zba * xt) * 2. * dphidyibt / rt2 - (xdc * (zdb * ycb 
			+ xu) - zdc * (xdb * ycb - zu)) / rcbu2 + (zu * xdb - 
			zdb * xu) * 2. * dphidyibu / ru2;
		dyibxid = rcbxu * dphidyibu + dphidzu + (xdc * zcb * xcb + 
			zdc * yzcb2) / rcbu2;
		dyibyid = rcbyu * dphidyibu + ycb * yu / rcbu2;
		dyibzid = rcbzu * dphidyibu - dphidxu - (zdc * xcb * zcb + 
			xdc * xycb2) / rcbu2;
		dzibxic = -xcb * dphidzib / (rcb * rcb) + dphidyt + dphidyu - 
			(xca * (yba * xcb - zt) + yca * (zba * zcb + ycb * 
			yba)) / rcbt2 - (yt * zba - yba * zt) * 2. * 
			dphidzibt / rt2 + (ydc * (zdb * zcb + ycb * ydb) + 
			xdc * (ydb * xcb + zu)) / rcbu2 + (yu * zdb - ydb * 
			zu) * 2. * dphidzibu / ru2;
		dzibzic = -zcb * dphidzib / (rcb * rcb) - (xca * (yba * zcb + 
			xt) - yca * (xba * zcb - yt)) / rcbt2 - (xt * yba - 
			xba * yt) * 2. * dphidzibt / rt2 - (ydc * (xdb * zcb 
			+ yu) - xdc * (ydb * zcb - xu)) / rcbu2 + (xu * ydb - 
			xdb * yu) * 2. * dphidzibu / ru2;
		dzibxid = rcbxu * dphidzibu - dphidyu - (xdc * ycb * xcb + 
			ydc * yzcb2) / rcbu2;
		dzibyid = rcbyu * dphidzibu + dphidxu + (ydc * xcb * ycb + 
			xdc * xzcb2) / rcbu2;
		dzibzid = rcbzu * dphidzibu + zcb * zu / rcbu2;
		dxicxid = rcbxu * dphidxicu - xcb * (zdb * ycb - ydb * zcb) / 
			rcbu2;
		dxicyid = rcbyu * dphidxicu + dphidzu + (ydb * zcb * ycb + 
			zdb * xzcb2) / rcbu2;
		dxiczid = rcbzu * dphidxicu - dphidyu - (zdb * ycb * zcb + 
			ydb * xycb2) / rcbu2;
		dyicxid = rcbxu * dphidyicu - dphidzu - (xdb * zcb * xcb + 
			zdb * yzcb2) / rcbu2;
		dyicyid = rcbyu * dphidyicu - ycb * (xdb * zcb - zdb * xcb) / 
			rcbu2;
		dyiczid = rcbzu * dphidyicu + dphidxu + (zdb * xcb * zcb + 
			xdb * xycb2) / rcbu2;
		dzicxid = rcbxu * dphidzicu + dphidyu + (xdb * ycb * xcb + 
			ydb * yzcb2) / rcbu2;
		dzicyid = rcbyu * dphidzicu - dphidxu - (ydb * xcb * ycb + 
			xdb * xzcb2) / rcbu2;
		dziczid = rcbzu * dphidzicu - zcb * (ydb * xcb - xdb * ycb) / 
			rcbu2;
		dxidxid = rcbxu * dphidxid;
		dxidyid = rcbxu * dphidyid + zcb * rcb / ru2;
		dxidzid = rcbxu * dphidzid - ycb * rcb / ru2;
		dyidyid = rcbyu * dphidyid;
		dyidzid = rcbyu * dphidzid + xcb * rcb / ru2;
		dzidzid = rcbzu * dphidzid;

/*     get some second derivative chain rule terms by difference */

		dxiaxib = -dxiaxia - dxiaxic - dxiaxid;
		dxiayib = -dxiayia - dxiayic - dxiayid;
		dxiazib = -dxiazia - dxiazic - dxiazid;
		dyiayib = -dyiayia - dyiayic - dyiayid;
		dyiazib = -dyiazia - dyiazic - dyiazid;
		dziazib = -dziazia - dziazic - dziazid;
		dxibxib = -dxiaxib - dxibxic - dxibxid;
		dxibyib = -dyiaxib - dxibyic - dxibyid;
		dxibzib = -dxiazib - dzibxic - dzibxid;
		dxibzic = -dziaxib - dxibzib - dxibzid;
		dyibyib = -dyiayib - dyibyic - dyibyid;
		dyibzic = -dziayib - dyibzib - dyibzid;
		dzibzib = -dziazib - dzibzic - dzibzid;
		dzibyic = -dyiazib - dyibzib - dzibyid;
		dxicxic = -dxiaxic - dxibxic - dxicxid;
		dxicyic = -dyiaxic - dyibxic - dxicyid;
		dxiczic = -dziaxic - dzibxic - dxiczid;
		dyicyic = -dyiayic - dyibyic - dyicyid;
		dyiczic = -dziayic - dzibyic - dyiczid;
		dziczic = -dziazic - dzibzic - dziczid;

/*     increment diagonal and off-diagonal Hessian elements */

		if (*i__ == ia) {
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedphi * dxiaxia + 
			    d2edphi2 * dphidxia * dphidxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedphi * dxiayia + 
			    d2edphi2 * dphidxia * dphidyia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedphi * dxiazia + 
			    d2edphi2 * dphidxia * dphidzia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedphi * dxiayia + 
			    d2edphi2 * dphidxia * dphidyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedphi * dyiayia + 
			    d2edphi2 * dphidyia * dphidyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedphi * dyiazia + 
			    d2edphi2 * dphidyia * dphidzia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedphi * dxiazia + 
			    d2edphi2 * dphidxia * dphidzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedphi * dyiazia + 
			    d2edphi2 * dphidyia * dphidzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedphi * dziazia + 
			    d2edphi2 * dphidzia * dphidzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedphi * dxiaxib + 
			    d2edphi2 * dphidxia * dphidxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedphi * dyiaxib + 
			    d2edphi2 * dphidyia * dphidxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedphi * dziaxib + 
			    d2edphi2 * dphidzia * dphidxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedphi * dxiayib + 
			    d2edphi2 * dphidxia * dphidyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedphi * dyiayib + 
			    d2edphi2 * dphidyia * dphidyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedphi * dziayib + 
			    d2edphi2 * dphidzia * dphidyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedphi * dxiazib + 
			    d2edphi2 * dphidxia * dphidzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedphi * dyiazib + 
			    d2edphi2 * dphidyia * dphidzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedphi * dziazib + 
			    d2edphi2 * dphidzia * dphidzib;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedphi * dxiaxic + 
			    d2edphi2 * dphidxia * dphidxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedphi * dyiaxic + 
			    d2edphi2 * dphidyia * dphidxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedphi * dziaxic + 
			    d2edphi2 * dphidzia * dphidxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedphi * dxiayic + 
			    d2edphi2 * dphidxia * dphidyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedphi * dyiayic + 
			    d2edphi2 * dphidyia * dphidyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedphi * dziayic + 
			    d2edphi2 * dphidzia * dphidyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedphi * dxiazic + 
			    d2edphi2 * dphidxia * dphidzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedphi * dyiazic + 
			    d2edphi2 * dphidyia * dphidzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedphi * dziazic + 
			    d2edphi2 * dphidzia * dphidzic;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedphi * dxiaxid + 
			    d2edphi2 * dphidxia * dphidxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedphi * dyiaxid + 
			    d2edphi2 * dphidyia * dphidxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedphi * dziaxid + 
			    d2edphi2 * dphidzia * dphidxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedphi * dxiayid + 
			    d2edphi2 * dphidxia * dphidyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedphi * dyiayid + 
			    d2edphi2 * dphidyia * dphidyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedphi * dziayid + 
			    d2edphi2 * dphidzia * dphidyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedphi * dxiazid + 
			    d2edphi2 * dphidxia * dphidzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedphi * dyiazid + 
			    d2edphi2 * dphidyia * dphidzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedphi * dziazid + 
			    d2edphi2 * dphidzia * dphidzid;
		} else if (*i__ == ib) {
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedphi * dxibxib + 
			    d2edphi2 * dphidxib * dphidxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedphi * dxibyib + 
			    d2edphi2 * dphidxib * dphidyib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedphi * dxibzib + 
			    d2edphi2 * dphidxib * dphidzib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedphi * dxibyib + 
			    d2edphi2 * dphidxib * dphidyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedphi * dyibyib + 
			    d2edphi2 * dphidyib * dphidyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedphi * dyibzib + 
			    d2edphi2 * dphidyib * dphidzib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedphi * dxibzib + 
			    d2edphi2 * dphidxib * dphidzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedphi * dyibzib + 
			    d2edphi2 * dphidyib * dphidzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedphi * dzibzib + 
			    d2edphi2 * dphidzib * dphidzib;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedphi * dxiaxib + 
			    d2edphi2 * dphidxib * dphidxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedphi * dxiayib + 
			    d2edphi2 * dphidyib * dphidxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedphi * dxiazib + 
			    d2edphi2 * dphidzib * dphidxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedphi * dyiaxib + 
			    d2edphi2 * dphidxib * dphidyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedphi * dyiayib + 
			    d2edphi2 * dphidyib * dphidyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedphi * dyiazib + 
			    d2edphi2 * dphidzib * dphidyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedphi * dziaxib + 
			    d2edphi2 * dphidxib * dphidzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedphi * dziayib + 
			    d2edphi2 * dphidyib * dphidzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedphi * dziazib + 
			    d2edphi2 * dphidzib * dphidzia;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedphi * dxibxic + 
			    d2edphi2 * dphidxib * dphidxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedphi * dyibxic + 
			    d2edphi2 * dphidyib * dphidxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedphi * dzibxic + 
			    d2edphi2 * dphidzib * dphidxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedphi * dxibyic + 
			    d2edphi2 * dphidxib * dphidyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedphi * dyibyic + 
			    d2edphi2 * dphidyib * dphidyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedphi * dzibyic + 
			    d2edphi2 * dphidzib * dphidyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedphi * dxibzic + 
			    d2edphi2 * dphidxib * dphidzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedphi * dyibzic + 
			    d2edphi2 * dphidyib * dphidzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedphi * dzibzic + 
			    d2edphi2 * dphidzib * dphidzic;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedphi * dxibxid + 
			    d2edphi2 * dphidxib * dphidxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedphi * dyibxid + 
			    d2edphi2 * dphidyib * dphidxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedphi * dzibxid + 
			    d2edphi2 * dphidzib * dphidxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedphi * dxibyid + 
			    d2edphi2 * dphidxib * dphidyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedphi * dyibyid + 
			    d2edphi2 * dphidyib * dphidyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedphi * dzibyid + 
			    d2edphi2 * dphidzib * dphidyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedphi * dxibzid + 
			    d2edphi2 * dphidxib * dphidzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedphi * dyibzid + 
			    d2edphi2 * dphidyib * dphidzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedphi * dzibzid + 
			    d2edphi2 * dphidzib * dphidzid;
		} else if (*i__ == ic) {
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedphi * dxicxic + 
			    d2edphi2 * dphidxic * dphidxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedphi * dxicyic + 
			    d2edphi2 * dphidxic * dphidyic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedphi * dxiczic + 
			    d2edphi2 * dphidxic * dphidzic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedphi * dxicyic + 
			    d2edphi2 * dphidxic * dphidyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedphi * dyicyic + 
			    d2edphi2 * dphidyic * dphidyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedphi * dyiczic + 
			    d2edphi2 * dphidyic * dphidzic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedphi * dxiczic + 
			    d2edphi2 * dphidxic * dphidzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedphi * dyiczic + 
			    d2edphi2 * dphidyic * dphidzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedphi * dziczic + 
			    d2edphi2 * dphidzic * dphidzic;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedphi * dxiaxic + 
			    d2edphi2 * dphidxic * dphidxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedphi * dxiayic + 
			    d2edphi2 * dphidyic * dphidxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedphi * dxiazic + 
			    d2edphi2 * dphidzic * dphidxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedphi * dyiaxic + 
			    d2edphi2 * dphidxic * dphidyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedphi * dyiayic + 
			    d2edphi2 * dphidyic * dphidyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedphi * dyiazic + 
			    d2edphi2 * dphidzic * dphidyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedphi * dziaxic + 
			    d2edphi2 * dphidxic * dphidzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedphi * dziayic + 
			    d2edphi2 * dphidyic * dphidzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedphi * dziazic + 
			    d2edphi2 * dphidzic * dphidzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedphi * dxibxic + 
			    d2edphi2 * dphidxic * dphidxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedphi * dxibyic + 
			    d2edphi2 * dphidyic * dphidxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedphi * dxibzic + 
			    d2edphi2 * dphidzic * dphidxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedphi * dyibxic + 
			    d2edphi2 * dphidxic * dphidyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedphi * dyibyic + 
			    d2edphi2 * dphidyic * dphidyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedphi * dyibzic + 
			    d2edphi2 * dphidzic * dphidyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedphi * dzibxic + 
			    d2edphi2 * dphidxic * dphidzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedphi * dzibyic + 
			    d2edphi2 * dphidyic * dphidzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedphi * dzibzic + 
			    d2edphi2 * dphidzic * dphidzib;
		    hessx_ref(1, id) = hessx_ref(1, id) + dedphi * dxicxid + 
			    d2edphi2 * dphidxic * dphidxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedphi * dyicxid + 
			    d2edphi2 * dphidyic * dphidxid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedphi * dzicxid + 
			    d2edphi2 * dphidzic * dphidxid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedphi * dxicyid + 
			    d2edphi2 * dphidxic * dphidyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedphi * dyicyid + 
			    d2edphi2 * dphidyic * dphidyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedphi * dzicyid + 
			    d2edphi2 * dphidzic * dphidyid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedphi * dxiczid + 
			    d2edphi2 * dphidxic * dphidzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedphi * dyiczid + 
			    d2edphi2 * dphidyic * dphidzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedphi * dziczid + 
			    d2edphi2 * dphidzic * dphidzid;
		} else if (*i__ == id) {
		    hessx_ref(1, id) = hessx_ref(1, id) + dedphi * dxidxid + 
			    d2edphi2 * dphidxid * dphidxid;
		    hessy_ref(1, id) = hessy_ref(1, id) + dedphi * dxidyid + 
			    d2edphi2 * dphidxid * dphidyid;
		    hessz_ref(1, id) = hessz_ref(1, id) + dedphi * dxidzid + 
			    d2edphi2 * dphidxid * dphidzid;
		    hessx_ref(2, id) = hessx_ref(2, id) + dedphi * dxidyid + 
			    d2edphi2 * dphidxid * dphidyid;
		    hessy_ref(2, id) = hessy_ref(2, id) + dedphi * dyidyid + 
			    d2edphi2 * dphidyid * dphidyid;
		    hessz_ref(2, id) = hessz_ref(2, id) + dedphi * dyidzid + 
			    d2edphi2 * dphidyid * dphidzid;
		    hessx_ref(3, id) = hessx_ref(3, id) + dedphi * dxidzid + 
			    d2edphi2 * dphidxid * dphidzid;
		    hessy_ref(3, id) = hessy_ref(3, id) + dedphi * dyidzid + 
			    d2edphi2 * dphidyid * dphidzid;
		    hessz_ref(3, id) = hessz_ref(3, id) + dedphi * dzidzid + 
			    d2edphi2 * dphidzid * dphidzid;
		    hessx_ref(1, ia) = hessx_ref(1, ia) + dedphi * dxiaxid + 
			    d2edphi2 * dphidxid * dphidxia;
		    hessy_ref(1, ia) = hessy_ref(1, ia) + dedphi * dxiayid + 
			    d2edphi2 * dphidyid * dphidxia;
		    hessz_ref(1, ia) = hessz_ref(1, ia) + dedphi * dxiazid + 
			    d2edphi2 * dphidzid * dphidxia;
		    hessx_ref(2, ia) = hessx_ref(2, ia) + dedphi * dyiaxid + 
			    d2edphi2 * dphidxid * dphidyia;
		    hessy_ref(2, ia) = hessy_ref(2, ia) + dedphi * dyiayid + 
			    d2edphi2 * dphidyid * dphidyia;
		    hessz_ref(2, ia) = hessz_ref(2, ia) + dedphi * dyiazid + 
			    d2edphi2 * dphidzid * dphidyia;
		    hessx_ref(3, ia) = hessx_ref(3, ia) + dedphi * dziaxid + 
			    d2edphi2 * dphidxid * dphidzia;
		    hessy_ref(3, ia) = hessy_ref(3, ia) + dedphi * dziayid + 
			    d2edphi2 * dphidyid * dphidzia;
		    hessz_ref(3, ia) = hessz_ref(3, ia) + dedphi * dziazid + 
			    d2edphi2 * dphidzid * dphidzia;
		    hessx_ref(1, ib) = hessx_ref(1, ib) + dedphi * dxibxid + 
			    d2edphi2 * dphidxid * dphidxib;
		    hessy_ref(1, ib) = hessy_ref(1, ib) + dedphi * dxibyid + 
			    d2edphi2 * dphidyid * dphidxib;
		    hessz_ref(1, ib) = hessz_ref(1, ib) + dedphi * dxibzid + 
			    d2edphi2 * dphidzid * dphidxib;
		    hessx_ref(2, ib) = hessx_ref(2, ib) + dedphi * dyibxid + 
			    d2edphi2 * dphidxid * dphidyib;
		    hessy_ref(2, ib) = hessy_ref(2, ib) + dedphi * dyibyid + 
			    d2edphi2 * dphidyid * dphidyib;
		    hessz_ref(2, ib) = hessz_ref(2, ib) + dedphi * dyibzid + 
			    d2edphi2 * dphidzid * dphidyib;
		    hessx_ref(3, ib) = hessx_ref(3, ib) + dedphi * dzibxid + 
			    d2edphi2 * dphidxid * dphidzib;
		    hessy_ref(3, ib) = hessy_ref(3, ib) + dedphi * dzibyid + 
			    d2edphi2 * dphidyid * dphidzib;
		    hessz_ref(3, ib) = hessz_ref(3, ib) + dedphi * dzibzid + 
			    d2edphi2 * dphidzid * dphidzib;
		    hessx_ref(1, ic) = hessx_ref(1, ic) + dedphi * dxicxid + 
			    d2edphi2 * dphidxid * dphidxic;
		    hessy_ref(1, ic) = hessy_ref(1, ic) + dedphi * dxicyid + 
			    d2edphi2 * dphidyid * dphidxic;
		    hessz_ref(1, ic) = hessz_ref(1, ic) + dedphi * dxiczid + 
			    d2edphi2 * dphidzid * dphidxic;
		    hessx_ref(2, ic) = hessx_ref(2, ic) + dedphi * dyicxid + 
			    d2edphi2 * dphidxid * dphidyic;
		    hessy_ref(2, ic) = hessy_ref(2, ic) + dedphi * dyicyid + 
			    d2edphi2 * dphidyid * dphidyic;
		    hessz_ref(2, ic) = hessz_ref(2, ic) + dedphi * dyiczid + 
			    d2edphi2 * dphidzid * dphidyic;
		    hessx_ref(3, ic) = hessx_ref(3, ic) + dedphi * dzicxid + 
			    d2edphi2 * dphidxid * dphidzic;
		    hessy_ref(3, ic) = hessy_ref(3, ic) + dedphi * dzicyid + 
			    d2edphi2 * dphidyid * dphidzic;
		    hessz_ref(3, ic) = hessz_ref(3, ic) + dedphi * dziczid + 
			    d2edphi2 * dphidzid * dphidzic;
		}
	    }
	}
    }

/*     compute the Hessian elements for group distance restraints */

    i__1 = kgeoms_1.ngfix;
    for (kdist = 1; kdist <= i__1; ++kdist) {
	ia = igfix_ref(1, kdist);
	ib = igfix_ref(2, kdist);
	proceed = group_1.grplist[*i__ - 1] == ia || group_1.grplist[*i__ - 1]
		 == ib;
	if (proceed) {
	    if (group_1.grplist[*i__ - 1] == ib) {
		ib = ia;
		ia = group_1.grplist[*i__ - 1];
	    }
	    xcm = 0.;
	    ycm = 0.;
	    zcm = 0.;
	    i__2 = igrp_ref(2, ia);
	    for (j = igrp_ref(1, ia); j <= i__2; ++j) {
		k = group_1.kgrp[j - 1];
		weigh = atmtyp_1.mass[k - 1];
		xcm += atoms_1.x[k - 1] * weigh;
		ycm += atoms_1.y[k - 1] * weigh;
		zcm += atoms_1.z__[k - 1] * weigh;
	    }
/* Computing MAX */
	    d__1 = 1., d__2 = group_1.grpmass[ia - 1];
	    weigha = max(d__1,d__2);
	    xr = xcm / weigha;
	    yr = ycm / weigha;
	    zr = zcm / weigha;
	    xcm = 0.;
	    ycm = 0.;
	    zcm = 0.;
	    i__2 = igrp_ref(2, ib);
	    for (j = igrp_ref(1, ib); j <= i__2; ++j) {
		k = group_1.kgrp[j - 1];
		weigh = atmtyp_1.mass[k - 1];
		xcm += atoms_1.x[k - 1] * weigh;
		ycm += atoms_1.y[k - 1] * weigh;
		zcm += atoms_1.z__[k - 1] * weigh;
	    }
/* Computing MAX */
	    d__1 = 1., d__2 = group_1.grpmass[ib - 1];
	    weighb = max(d__1,d__2);
	    xr -= xcm / weighb;
	    yr -= ycm / weighb;
	    zr -= zcm / weighb;
	    intermol = molcul_1.molcule[group_1.kgrp[igrp_ref(1, ia) - 1] - 1]
		     != molcul_1.molcule[group_1.kgrp[igrp_ref(1, ib) - 1] - 
		    1];
	    if (bound_1.use_bounds__ && intermol) {
		image_(&xr, &yr, &zr);
	    }
	    r2 = xr * xr + yr * yr + zr * zr;
	    r__ = sqrt(r2);
	    force = gfix_ref(1, kdist);
	    gf1 = gfix_ref(2, kdist);
	    gf2 = gfix_ref(3, kdist);
	    target = r__;
	    if (r__ < gf1) {
		target = gf1;
	    }
	    if (r__ > gf2) {
		target = gf2;
	    }
	    dt = r__ - target;
	    deddt = force * 2.;

/*     set the chain rule terms for the Hessian elements */

	    if (r__ == 0.) {
		r__ = 1e-4;
		r2 = r__ * r__;
	    }
	    de = deddt * dt / r__;
	    term = (deddt - de) / r2;
	    termx = term * xr;
	    termy = term * yr;
	    termz = term * zr;
	    d2e_ref(1, 1) = termx * xr + de;
	    d2e_ref(1, 2) = termx * yr;
	    d2e_ref(1, 3) = termx * zr;
	    d2e_ref(2, 1) = d2e_ref(1, 2);
	    d2e_ref(2, 2) = termy * yr + de;
	    d2e_ref(2, 3) = termy * zr;
	    d2e_ref(3, 1) = d2e_ref(1, 3);
	    d2e_ref(3, 2) = d2e_ref(2, 3);
	    d2e_ref(3, 3) = termz * zr + de;

/*     increment diagonal and non-diagonal Hessian elements */

	    i__2 = igrp_ref(2, ia);
	    for (k = igrp_ref(1, ia); k <= i__2; ++k) {
		m = group_1.kgrp[k - 1];
		ratio = atmtyp_1.mass[*i__ - 1] * atmtyp_1.mass[m - 1] / (
			weigha * weigha);
		for (j = 1; j <= 3; ++j) {
		    hessx_ref(j, m) = hessx_ref(j, m) + d2e_ref(1, j) * ratio;
		    hessy_ref(j, m) = hessy_ref(j, m) + d2e_ref(2, j) * ratio;
		    hessz_ref(j, m) = hessz_ref(j, m) + d2e_ref(3, j) * ratio;
		}
	    }
	    i__2 = igrp_ref(2, ib);
	    for (k = igrp_ref(1, ib); k <= i__2; ++k) {
		m = group_1.kgrp[k - 1];
		ratio = atmtyp_1.mass[*i__ - 1] * atmtyp_1.mass[m - 1] / (
			weigha * weighb);
		for (j = 1; j <= 3; ++j) {
		    hessx_ref(j, m) = hessx_ref(j, m) - d2e_ref(1, j) * ratio;
		    hessy_ref(j, m) = hessy_ref(j, m) - d2e_ref(2, j) * ratio;
		    hessz_ref(j, m) = hessz_ref(j, m) - d2e_ref(3, j) * ratio;
		}
	    }
	}
    }

/*     compute the Hessian elements for chirality restraints */

    i__1 = kgeoms_1.nchir;
    for (kchir = 1; kchir <= i__1; ++kchir) {
	ia = ichir_ref(1, kchir);
	ib = ichir_ref(2, kchir);
	ic = ichir_ref(3, kchir);
	id = ichir_ref(4, kchir);
	proceed = *i__ == ia || *i__ == ib || *i__ == ic || *i__ == id;
	if (proceed && group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &c__0, &c__0);
	}
	if (proceed) {
	    xad = atoms_1.x[ia - 1] - atoms_1.x[id - 1];
	    yad = atoms_1.y[ia - 1] - atoms_1.y[id - 1];
	    zad = atoms_1.z__[ia - 1] - atoms_1.z__[id - 1];
	    xbd = atoms_1.x[ib - 1] - atoms_1.x[id - 1];
	    ybd = atoms_1.y[ib - 1] - atoms_1.y[id - 1];
	    zbd = atoms_1.z__[ib - 1] - atoms_1.z__[id - 1];
	    xcd = atoms_1.x[ic - 1] - atoms_1.x[id - 1];
	    ycd = atoms_1.y[ic - 1] - atoms_1.y[id - 1];
	    zcd = atoms_1.z__[ic - 1] - atoms_1.z__[id - 1];
	    c1 = ybd * zcd - zbd * ycd;
	    c2 = ycd * zad - zcd * yad;
	    c3 = yad * zbd - zad * ybd;
	    vol = xad * c1 + xbd * c2 + xcd * c3;
	    force = chir_ref(1, kchir);
	    cf1 = chir_ref(2, kchir);
	    cf2 = chir_ref(3, kchir);
	    target = vol;
	    if (vol < min(cf1,cf2)) {
		target = min(cf1,cf2);
	    }
	    if (vol > max(cf1,cf2)) {
		target = max(cf1,cf2);
	    }
	    dt = vol - target;
	    dt2 = dt * dt;
	    deddt = force * 2. * dt;
	    d2eddt2 = force * 2.;

/*     scale the interaction based on its group membership */

	    if (group_1.use_group__) {
		deddt *= fgrp;
		d2eddt2 *= fgrp;
	    }

/*     chain rule terms for first derivative components */

	    term = sqrt(d2eddt2);
	    ddtdxia = term * (ybd * zcd - zbd * ycd);
	    ddtdyia = term * (zbd * xcd - xbd * zcd);
	    ddtdzia = term * (xbd * ycd - ybd * xcd);
	    ddtdxib = term * (zad * ycd - yad * zcd);
	    ddtdyib = term * (xad * zcd - zad * xcd);
	    ddtdzib = term * (yad * xcd - xad * ycd);
	    ddtdxic = term * (yad * zbd - zad * ybd);
	    ddtdyic = term * (zad * xbd - xad * zbd);
	    ddtdzic = term * (xad * ybd - yad * xbd);
	    ddtdxid = -ddtdxia - ddtdxib - ddtdxic;
	    ddtdyid = -ddtdyia - ddtdyib - ddtdyic;
	    ddtdzid = -ddtdzia - ddtdzib - ddtdzic;

/*     chain rule terms for second derivative components (*deddt) */

	    dyiaxib = -deddt * zcd;
	    dziaxib = deddt * ycd;
	    dxiayib = deddt * zcd;
	    dziayib = -deddt * xcd;
	    dxiazib = -deddt * ycd;
	    dyiazib = deddt * xcd;
	    dyiaxic = deddt * zbd;
	    dziaxic = -deddt * ybd;
	    dxiayic = -deddt * zbd;
	    dziayic = deddt * xbd;
	    dxiazic = deddt * ybd;
	    dyiazic = -deddt * xbd;
	    dyibxic = -deddt * zad;
	    dzibxic = deddt * yad;
	    dxibyic = deddt * zad;
	    dzibyic = -deddt * xad;
	    dxibzic = -deddt * yad;
	    dyibzic = deddt * xad;
	    dyiaxid = -dyiaxib - dyiaxic;
	    dziaxid = -dziaxib - dziaxic;
	    dxiayid = -dxiayib - dxiayic;
	    dziayid = -dziayib - dziayic;
	    dxiazid = -dxiazib - dxiazic;
	    dyiazid = -dyiazib - dyiazic;
	    dyibxid = -dxiayib - dyibxic;
	    dzibxid = -dxiazib - dzibxic;
	    dxibyid = -dyiaxib - dxibyic;
	    dzibyid = -dyiazib - dzibyic;
	    dxibzid = -dziaxib - dxibzic;
	    dyibzid = -dziayib - dyibzic;
	    dyicxid = -dxiayic - dxibyic;
	    dzicxid = -dxiazic - dxibzic;
	    dxicyid = -dyiaxic - dyibxic;
	    dzicyid = -dyiazic - dyibzic;
	    dxiczid = -dziaxic - dzibxic;
	    dyiczid = -dziayic - dzibyic;

/*     increment diagonal and off-diagonal Hessian elements */

	    if (*i__ == ia) {
		hessx_ref(1, ia) = hessx_ref(1, ia) + ddtdxia * ddtdxia;
		hessy_ref(1, ia) = hessy_ref(1, ia) + ddtdxia * ddtdyia;
		hessz_ref(1, ia) = hessz_ref(1, ia) + ddtdxia * ddtdzia;
		hessx_ref(2, ia) = hessx_ref(2, ia) + ddtdxia * ddtdyia;
		hessy_ref(2, ia) = hessy_ref(2, ia) + ddtdyia * ddtdyia;
		hessz_ref(2, ia) = hessz_ref(2, ia) + ddtdyia * ddtdzia;
		hessx_ref(3, ia) = hessx_ref(3, ia) + ddtdxia * ddtdzia;
		hessy_ref(3, ia) = hessy_ref(3, ia) + ddtdyia * ddtdzia;
		hessz_ref(3, ia) = hessz_ref(3, ia) + ddtdzia * ddtdzia;
		hessx_ref(1, ib) = hessx_ref(1, ib) + ddtdxia * ddtdxib;
		hessy_ref(1, ib) = hessy_ref(1, ib) + ddtdyia * ddtdxib + 
			dyiaxib;
		hessz_ref(1, ib) = hessz_ref(1, ib) + ddtdzia * ddtdxib + 
			dziaxib;
		hessx_ref(2, ib) = hessx_ref(2, ib) + ddtdxia * ddtdyib + 
			dxiayib;
		hessy_ref(2, ib) = hessy_ref(2, ib) + ddtdyia * ddtdyib;
		hessz_ref(2, ib) = hessz_ref(2, ib) + ddtdzia * ddtdyib + 
			dziayib;
		hessx_ref(3, ib) = hessx_ref(3, ib) + ddtdxia * ddtdzib + 
			dxiazib;
		hessy_ref(3, ib) = hessy_ref(3, ib) + ddtdyia * ddtdzib + 
			dyiazib;
		hessz_ref(3, ib) = hessz_ref(3, ib) + ddtdzia * ddtdzib;
		hessx_ref(1, ic) = hessx_ref(1, ic) + ddtdxia * ddtdxic;
		hessy_ref(1, ic) = hessy_ref(1, ic) + ddtdyia * ddtdxic + 
			dyiaxic;
		hessz_ref(1, ic) = hessz_ref(1, ic) + ddtdzia * ddtdxic + 
			dziaxic;
		hessx_ref(2, ic) = hessx_ref(2, ic) + ddtdxia * ddtdyic + 
			dxiayic;
		hessy_ref(2, ic) = hessy_ref(2, ic) + ddtdyia * ddtdyic;
		hessz_ref(2, ic) = hessz_ref(2, ic) + ddtdzia * ddtdyic + 
			dziayic;
		hessx_ref(3, ic) = hessx_ref(3, ic) + ddtdxia * ddtdzic + 
			dxiazic;
		hessy_ref(3, ic) = hessy_ref(3, ic) + ddtdyia * ddtdzic + 
			dyiazic;
		hessz_ref(3, ic) = hessz_ref(3, ic) + ddtdzia * ddtdzic;
		hessx_ref(1, id) = hessx_ref(1, id) + ddtdxia * ddtdxid;
		hessy_ref(1, id) = hessy_ref(1, id) + ddtdyia * ddtdxid + 
			dyiaxid;
		hessz_ref(1, id) = hessz_ref(1, id) + ddtdzia * ddtdxid + 
			dziaxid;
		hessx_ref(2, id) = hessx_ref(2, id) + ddtdxia * ddtdyid + 
			dxiayid;
		hessy_ref(2, id) = hessy_ref(2, id) + ddtdyia * ddtdyid;
		hessz_ref(2, id) = hessz_ref(2, id) + ddtdzia * ddtdyid + 
			dziayid;
		hessx_ref(3, id) = hessx_ref(3, id) + ddtdxia * ddtdzid + 
			dxiazid;
		hessy_ref(3, id) = hessy_ref(3, id) + ddtdyia * ddtdzid + 
			dyiazid;
		hessz_ref(3, id) = hessz_ref(3, id) + ddtdzia * ddtdzid;
	    } else if (*i__ == ib) {
		hessx_ref(1, ib) = hessx_ref(1, ib) + ddtdxib * ddtdxib;
		hessy_ref(1, ib) = hessy_ref(1, ib) + ddtdxib * ddtdyib;
		hessz_ref(1, ib) = hessz_ref(1, ib) + ddtdxib * ddtdzib;
		hessx_ref(2, ib) = hessx_ref(2, ib) + ddtdxib * ddtdyib;
		hessy_ref(2, ib) = hessy_ref(2, ib) + ddtdyib * ddtdyib;
		hessz_ref(2, ib) = hessz_ref(2, ib) + ddtdyib * ddtdzib;
		hessx_ref(3, ib) = hessx_ref(3, ib) + ddtdxib * ddtdzib;
		hessy_ref(3, ib) = hessy_ref(3, ib) + ddtdyib * ddtdzib;
		hessz_ref(3, ib) = hessz_ref(3, ib) + ddtdzib * ddtdzib;
		hessx_ref(1, ia) = hessx_ref(1, ia) + ddtdxib * ddtdxia;
		hessy_ref(1, ia) = hessy_ref(1, ia) + ddtdyib * ddtdxia + 
			dxiayib;
		hessz_ref(1, ia) = hessz_ref(1, ia) + ddtdzib * ddtdxia + 
			dxiazib;
		hessx_ref(2, ia) = hessx_ref(2, ia) + ddtdxib * ddtdyia + 
			dyiaxib;
		hessy_ref(2, ia) = hessy_ref(2, ia) + ddtdyib * ddtdyia;
		hessz_ref(2, ia) = hessz_ref(2, ia) + ddtdzib * ddtdyia + 
			dyiazib;
		hessx_ref(3, ia) = hessx_ref(3, ia) + ddtdxib * ddtdzia + 
			dziaxib;
		hessy_ref(3, ia) = hessy_ref(3, ia) + ddtdyib * ddtdzia + 
			dziayib;
		hessz_ref(3, ia) = hessz_ref(3, ia) + ddtdzib * ddtdzia;
		hessx_ref(1, ic) = hessx_ref(1, ic) + ddtdxib * ddtdxic;
		hessy_ref(1, ic) = hessy_ref(1, ic) + ddtdyib * ddtdxic + 
			dyibxic;
		hessz_ref(1, ic) = hessz_ref(1, ic) + ddtdzib * ddtdxic + 
			dzibxic;
		hessx_ref(2, ic) = hessx_ref(2, ic) + ddtdxib * ddtdyic + 
			dxibyic;
		hessy_ref(2, ic) = hessy_ref(2, ic) + ddtdyib * ddtdyic;
		hessz_ref(2, ic) = hessz_ref(2, ic) + ddtdzib * ddtdyic + 
			dzibyic;
		hessx_ref(3, ic) = hessx_ref(3, ic) + ddtdxib * ddtdzic + 
			dxibzic;
		hessy_ref(3, ic) = hessy_ref(3, ic) + ddtdyib * ddtdzic + 
			dyibzic;
		hessz_ref(3, ic) = hessz_ref(3, ic) + ddtdzib * ddtdzic;
		hessx_ref(1, id) = hessx_ref(1, id) + ddtdxib * ddtdxid;
		hessy_ref(1, id) = hessy_ref(1, id) + ddtdyib * ddtdxid + 
			dyibxid;
		hessz_ref(1, id) = hessz_ref(1, id) + ddtdzib * ddtdxid + 
			dzibxid;
		hessx_ref(2, id) = hessx_ref(2, id) + ddtdxib * ddtdyid + 
			dxibyid;
		hessy_ref(2, id) = hessy_ref(2, id) + ddtdyib * ddtdyid;
		hessz_ref(2, id) = hessz_ref(2, id) + ddtdzib * ddtdyid + 
			dzibyid;
		hessx_ref(3, id) = hessx_ref(3, id) + ddtdxib * ddtdzid + 
			dxibzid;
		hessy_ref(3, id) = hessy_ref(3, id) + ddtdyib * ddtdzid + 
			dyibzid;
		hessz_ref(3, id) = hessz_ref(3, id) + ddtdzib * ddtdzid;
	    } else if (*i__ == ic) {
		hessx_ref(1, ic) = hessx_ref(1, ic) + ddtdxic * ddtdxic;
		hessy_ref(1, ic) = hessy_ref(1, ic) + ddtdxic * ddtdyic;
		hessz_ref(1, ic) = hessz_ref(1, ic) + ddtdxic * ddtdzic;
		hessx_ref(2, ic) = hessx_ref(2, ic) + ddtdxic * ddtdyic;
		hessy_ref(2, ic) = hessy_ref(2, ic) + ddtdyic * ddtdyic;
		hessz_ref(2, ic) = hessz_ref(2, ic) + ddtdyic * ddtdzic;
		hessx_ref(3, ic) = hessx_ref(3, ic) + ddtdxic * ddtdzic;
		hessy_ref(3, ic) = hessy_ref(3, ic) + ddtdyic * ddtdzic;
		hessz_ref(3, ic) = hessz_ref(3, ic) + ddtdzic * ddtdzic;
		hessx_ref(1, ia) = hessx_ref(1, ia) + ddtdxic * ddtdxia;
		hessy_ref(1, ia) = hessy_ref(1, ia) + ddtdyic * ddtdxia + 
			dxiayic;
		hessz_ref(1, ia) = hessz_ref(1, ia) + ddtdzic * ddtdxia + 
			dxiazic;
		hessx_ref(2, ia) = hessx_ref(2, ia) + ddtdxic * ddtdyia + 
			dyiaxic;
		hessy_ref(2, ia) = hessy_ref(2, ia) + ddtdyic * ddtdyia;
		hessz_ref(2, ia) = hessz_ref(2, ia) + ddtdzic * ddtdyia + 
			dyiazic;
		hessx_ref(3, ia) = hessx_ref(3, ia) + ddtdxic * ddtdzia + 
			dziaxic;
		hessy_ref(3, ia) = hessy_ref(3, ia) + ddtdyic * ddtdzia + 
			dziayic;
		hessz_ref(3, ia) = hessz_ref(3, ia) + ddtdzic * ddtdzia;
		hessx_ref(1, ib) = hessx_ref(1, ib) + ddtdxic * ddtdxib;
		hessy_ref(1, ib) = hessy_ref(1, ib) + ddtdyic * ddtdxib + 
			dxibyic;
		hessz_ref(1, ib) = hessz_ref(1, ib) + ddtdzic * ddtdxib + 
			dxibzic;
		hessx_ref(2, ib) = hessx_ref(2, ib) + ddtdxic * ddtdyib + 
			dyibxic;
		hessy_ref(2, ib) = hessy_ref(2, ib) + ddtdyic * ddtdyib;
		hessz_ref(2, ib) = hessz_ref(2, ib) + ddtdzic * ddtdyib + 
			dyibzic;
		hessx_ref(3, ib) = hessx_ref(3, ib) + ddtdxic * ddtdzib + 
			dzibxic;
		hessy_ref(3, ib) = hessy_ref(3, ib) + ddtdyic * ddtdzib + 
			dzibyic;
		hessz_ref(3, ib) = hessz_ref(3, ib) + ddtdzic * ddtdzib;
		hessx_ref(1, id) = hessx_ref(1, id) + ddtdxic * ddtdxid;
		hessy_ref(1, id) = hessy_ref(1, id) + ddtdyic * ddtdxid + 
			dyicxid;
		hessz_ref(1, id) = hessz_ref(1, id) + ddtdzic * ddtdxid + 
			dzicxid;
		hessx_ref(2, id) = hessx_ref(2, id) + ddtdxic * ddtdyid + 
			dxicyid;
		hessy_ref(2, id) = hessy_ref(2, id) + ddtdyic * ddtdyid;
		hessz_ref(2, id) = hessz_ref(2, id) + ddtdzic * ddtdyid + 
			dzicyid;
		hessx_ref(3, id) = hessx_ref(3, id) + ddtdxic * ddtdzid + 
			dxiczid;
		hessy_ref(3, id) = hessy_ref(3, id) + ddtdyic * ddtdzid + 
			dyiczid;
		hessz_ref(3, id) = hessz_ref(3, id) + ddtdzic * ddtdzid;
	    } else if (*i__ == id) {
		hessx_ref(1, id) = hessx_ref(1, id) + ddtdxid * ddtdxid;
		hessy_ref(1, id) = hessy_ref(1, id) + ddtdxid * ddtdyid;
		hessz_ref(1, id) = hessz_ref(1, id) + ddtdxid * ddtdzid;
		hessx_ref(2, id) = hessx_ref(2, id) + ddtdxid * ddtdyid;
		hessy_ref(2, id) = hessy_ref(2, id) + ddtdyid * ddtdyid;
		hessz_ref(2, id) = hessz_ref(2, id) + ddtdyid * ddtdzid;
		hessx_ref(3, id) = hessx_ref(3, id) + ddtdxid * ddtdzid;
		hessy_ref(3, id) = hessy_ref(3, id) + ddtdyid * ddtdzid;
		hessz_ref(3, id) = hessz_ref(3, id) + ddtdzid * ddtdzid;
		hessx_ref(1, ia) = hessx_ref(1, ia) + ddtdxid * ddtdxia;
		hessy_ref(1, ia) = hessy_ref(1, ia) + ddtdyid * ddtdxia + 
			dxiayid;
		hessz_ref(1, ia) = hessz_ref(1, ia) + ddtdzid * ddtdxia + 
			dxiazid;
		hessx_ref(2, ia) = hessx_ref(2, ia) + ddtdxid * ddtdyia + 
			dyiaxid;
		hessy_ref(2, ia) = hessy_ref(2, ia) + ddtdyid * ddtdyia;
		hessz_ref(2, ia) = hessz_ref(2, ia) + ddtdzid * ddtdyia + 
			dyiazid;
		hessx_ref(3, ia) = hessx_ref(3, ia) + ddtdxid * ddtdzia + 
			dziaxid;
		hessy_ref(3, ia) = hessy_ref(3, ia) + ddtdyid * ddtdzia + 
			dziayid;
		hessz_ref(3, ia) = hessz_ref(3, ia) + ddtdzid * ddtdzia;
		hessx_ref(1, ib) = hessx_ref(1, ib) + ddtdxid * ddtdxib;
		hessy_ref(1, ib) = hessy_ref(1, ib) + ddtdyid * ddtdxib + 
			dxibyid;
		hessz_ref(1, ib) = hessz_ref(1, ib) + ddtdzid * ddtdxib + 
			dxibzid;
		hessx_ref(2, ib) = hessx_ref(2, ib) + ddtdxid * ddtdyib + 
			dyibxid;
		hessy_ref(2, ib) = hessy_ref(2, ib) + ddtdyid * ddtdyib;
		hessz_ref(2, ib) = hessz_ref(2, ib) + ddtdzid * ddtdyib + 
			dyibzid;
		hessx_ref(3, ib) = hessx_ref(3, ib) + ddtdxid * ddtdzib + 
			dzibxid;
		hessy_ref(3, ib) = hessy_ref(3, ib) + ddtdyid * ddtdzib + 
			dzibyid;
		hessz_ref(3, ib) = hessz_ref(3, ib) + ddtdzid * ddtdzib;
		hessx_ref(1, ic) = hessx_ref(1, ic) + ddtdxid * ddtdxic;
		hessy_ref(1, ic) = hessy_ref(1, ic) + ddtdyid * ddtdxic + 
			dxicyid;
		hessz_ref(1, ic) = hessz_ref(1, ic) + ddtdzid * ddtdxic + 
			dxiczid;
		hessx_ref(2, ic) = hessx_ref(2, ic) + ddtdxid * ddtdyic + 
			dyicxid;
		hessy_ref(2, ic) = hessy_ref(2, ic) + ddtdyid * ddtdyic;
		hessz_ref(2, ic) = hessz_ref(2, ic) + ddtdzid * ddtdyic + 
			dyiczid;
		hessx_ref(3, ic) = hessx_ref(3, ic) + ddtdxid * ddtdzic + 
			dzicxid;
		hessy_ref(3, ic) = hessy_ref(3, ic) + ddtdyid * ddtdzic + 
			dzicyid;
		hessz_ref(3, ic) = hessz_ref(3, ic) + ddtdzid * ddtdzic;
	    }
	}
    }

/*     compute Hessian elements for a Gaussian basin restraint */

    if (kgeoms_1.use_basin__) {
	xi = atoms_1.x[*i__ - 1];
	yi = atoms_1.y[*i__ - 1];
	zi = atoms_1.z__[*i__ - 1];
	i__1 = atoms_1.n;
	for (k = 1; k <= i__1; ++k) {
	    proceed = k != *i__;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, i__, &k, &c__0, &c__0, &c__0, &c__0);
	    }
	    if (proceed) {
		xr = xi - atoms_1.x[k - 1];
		yr = yi - atoms_1.y[k - 1];
		zr = zi - atoms_1.z__[k - 1];
		r2 = xr * xr + yr * yr + zr * zr;
		term = -kgeoms_1.width * r2;
		expterm = 0.;
		if (term > -50.) {
		    expterm = kgeoms_1.depth * kgeoms_1.width * exp(term);
		}
		dedr = expterm * -2.;
		d2edr2 = (term * -4. - 2.) * expterm;

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    dedr *= fgrp;
		    d2edr2 *= fgrp;
		}

/*     set the chain rule terms for the Hessian elements */

		if (r2 == 0.) {
		    term = 0.;
		} else {
		    term = (d2edr2 - dedr) / r2;
		}
		termx = term * xr;
		termy = term * yr;
		termz = term * zr;
		d2e_ref(1, 1) = termx * xr + dedr;
		d2e_ref(1, 2) = termx * yr;
		d2e_ref(1, 3) = termx * zr;
		d2e_ref(2, 1) = d2e_ref(1, 2);
		d2e_ref(2, 2) = termy * yr + dedr;
		d2e_ref(2, 3) = termy * zr;
		d2e_ref(3, 1) = d2e_ref(1, 3);
		d2e_ref(3, 2) = d2e_ref(2, 3);
		d2e_ref(3, 3) = termz * zr + dedr;

/*     increment diagonal and non-diagonal Hessian elements */

		for (j = 1; j <= 3; ++j) {
		    hessx_ref(j, *i__) = hessx_ref(j, *i__) + d2e_ref(1, j);
		    hessy_ref(j, *i__) = hessy_ref(j, *i__) + d2e_ref(2, j);
		    hessz_ref(j, *i__) = hessz_ref(j, *i__) + d2e_ref(3, j);
		    hessx_ref(j, k) = hessx_ref(j, k) - d2e_ref(1, j);
		    hessy_ref(j, k) = hessy_ref(j, k) - d2e_ref(2, j);
		    hessz_ref(j, k) = hessz_ref(j, k) - d2e_ref(3, j);
		}
	    }
	}
    }

/*     compute Hessian elements for a spherical droplet restraint */

    if (kgeoms_1.use_wall__) {
	buffer = 2.5;
	a = 2048.;
	b = 64.;
	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, i__, &c__0, &c__0, &c__0, &c__0, &c__0);
	}
	if (proceed) {
	    xi = atoms_1.x[*i__ - 1];
	    yi = atoms_1.y[*i__ - 1];
	    zi = atoms_1.z__[*i__ - 1];
/* Computing 2nd power */
	    d__1 = xi;
/* Computing 2nd power */
	    d__2 = yi;
/* Computing 2nd power */
	    d__3 = zi;
	    ri2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	    ri = sqrt(ri2);
	    r__ = kgeoms_1.rwall + buffer - ri;
	    r2 = r__ * r__;
	    r6 = r2 * r2 * r2;
	    r12 = r6 * r6;
	    if (ri == 0.) {
		ri = 1.;
		ri2 = 1.;
	    }
	    dedr = (a * 12. / r12 - b * 6. / r6) / (r__ * ri);
	    d2edr2 = (a * 156. / r12 - b * 42. / r6) / (r2 * ri2);

/*     scale the interaction based on its group membership */

	    if (group_1.use_group__) {
		dedr *= fgrp;
		d2edr2 *= fgrp;
	    }

/*     set the chain rule terms for the Hessian elements */

	    d2edr2 -= dedr / ri2;
	    termx = d2edr2 * xi;
	    termy = d2edr2 * yi;
	    termz = d2edr2 * zi;
	    d2e_ref(1, 1) = termx * xi + dedr;
	    d2e_ref(1, 2) = termx * yi;
	    d2e_ref(1, 3) = termx * zi;
	    d2e_ref(2, 1) = d2e_ref(1, 2);
	    d2e_ref(2, 2) = termy * yi + dedr;
	    d2e_ref(2, 3) = termy * zi;
	    d2e_ref(3, 1) = d2e_ref(1, 3);
	    d2e_ref(3, 2) = d2e_ref(2, 3);
	    d2e_ref(3, 3) = termz * zi + dedr;

/*     increment diagonal and non-diagonal Hessian elements */

	    for (j = 1; j <= 3; ++j) {
		hessx_ref(j, *i__) = hessx_ref(j, *i__) + d2e_ref(1, j);
		hessy_ref(j, *i__) = hessy_ref(j, *i__) + d2e_ref(2, j);
		hessz_ref(j, *i__) = hessz_ref(j, *i__) + d2e_ref(3, j);
	    }
	}
    }
    return 0;
} /* egeom2_ */

#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef itfix_ref
#undef kpfix_ref
#undef igfix_ref
#undef idfix_ref
#undef iafix_ref
#undef ichir_ref
#undef tfix_ref
#undef pfix_ref
#undef igrp_ref
#undef gfix_ref
#undef dfix_ref
#undef afix_ref
#undef chir_ref
#undef d2e_ref


