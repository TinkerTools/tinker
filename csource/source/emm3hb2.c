/* emm3hb2.f -- translated by f2c (version 20050501).
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
    integer bndlist[200000]	/* was [8][25000] */, anglist[700000]	/* 
	    was [28][25000] */;
} atmlst_;

#define atmlst_1 atmlst_

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
    doublereal bk[50000], bl[50000];
    integer nbond, ibnd[100000]	/* was [2][50000] */;
} bond_;

#define bond_1 bond_

struct {
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

struct {
    doublereal xcell, ycell, zcell, xcell2, ycell2, zcell2;
    integer ncell, icell[30000]	/* was [3][10000] */;
} cell_;

#define cell_1 cell_

struct {
    doublereal electric, dielec, ebuffer, c2scale, c3scale, c4scale, c5scale;
    logical neutnbr, neutcut;
} chgpot_;

#define chgpot_1 chgpot_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

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
    doublereal off, off2, cut, cut2, c0, c1, c2, c3, c4, c5, f0, f1, f2, f3, 
	    f4, f5, f6, f7;
} shunt_;

#define shunt_1 shunt_

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
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine emm3hb2  --  atomwise MM3 vdw & hbond Hessian  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "emm3hb2" calculates the MM3 exp-6 van der Waals and directional */
/*     charge transfer hydrogen bonding second derivatives for a single */
/*     atom at a time */

/*     note this version only partially incorporates the directional */
/*     hydrogen bonding term into the Hessian calculation */

/*     literature references: */

/*     J.-H. Lii and N. L. Allinger, "Directional Hydrogen Bonding in */
/*     the MM3 Force Field. I", Journal of Physical Organic Chemistry, */
/*     7, 591-609 (1994) */

/*     J.-H. Lii and N. L. Allinger, "Directional Hydrogen Bonding in */
/*     the MM3 Force Field. II", Journal of Computational Chemistry, */
/*     19, 1001-1016 (1998) */


/* Subroutine */ int emm3hb2_(integer *iatom, doublereal *xred, doublereal *
	yred, doublereal *zred)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal expmerge, e;
    static integer i__, j, k;
    static doublereal p, p2, p6, de;
    static integer ia, ib, ii, kk, ic;
    static doublereal p12;
    static integer it, iv, kt, kv;
    static doublereal xi, yi, zi, rp, rv, xp, xr, yr, zr, yp, zp, d2e, rab, 
	    xab, yab, zab, xcb, ycb, xia, yia;
    static integer iv14[25000];
    static doublereal zia, rik, dot, eps, xib, yib, zib, xic, yic, zic, zcb, 
	    rab2, rcb2, rik2, rik3, rik4, rik5, redi, redk, fgrp, sine, term[
	    9]	/* was [3][3] */;
    static integer list[5];
    static doublereal redi2, d2edx, d2edy, d2edz, ideal;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static integer jcell;
    static doublereal redik, rediv, redkv, taper, fterm, ratio;
    static integer nlist;
    static doublereal rediv2, fcbuck;
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal vscale[25000], rediiv, redivk, redikv, dtaper, cosine, 
	    expcut;
    extern /* Subroutine */ int switch_(char *, ftnlen);
    static doublereal rvterm;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal d2taper, expmin2, expcut2, rvterm2;
    static logical proceed;
    static doublereal redivkv, expterm;


#define epsilon4_ref(a_1,a_2) vdw_1.epsilon4[(a_2)*1000 + a_1 - 1001]
#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define term_ref(a_1,a_2) term[(a_2)*3 + a_1 - 4]
#define hessx_ref(a_1,a_2) hessn_1.hessx[(a_2)*3 + a_1 - 4]
#define hessy_ref(a_1,a_2) hessn_1.hessy[(a_2)*3 + a_1 - 4]
#define hessz_ref(a_1,a_2) hessn_1.hessz[(a_2)*3 + a_1 - 4]
#define radmin_ref(a_1,a_2) vdw_1.radmin[(a_2)*1000 + a_1 - 1001]
#define radmin4_ref(a_1,a_2) vdw_1.radmin4[(a_2)*1000 + a_1 - 1001]
#define radhbnd_ref(a_1,a_2) vdw_1.radhbnd[(a_2)*1000 + a_1 - 1001]
#define epshbnd_ref(a_1,a_2) vdw_1.epshbnd[(a_2)*1000 + a_1 - 1001]
#define bndlist_ref(a_1,a_2) atmlst_1.bndlist[(a_2)*8 + a_1 - 9]
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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  atmlst.i  --  local geometry terms involving each atom  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     bndlist   list of the bond numbers involving each atom */
/*     anglist   list of the angle numbers centered on each atom */




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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  bond.i  --  covalent bonds in the current structure  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     bk      bond stretch force constants (kcal/mole/Ang**2) */
/*     bl      ideal bond length values in Angstroms */
/*     nbond   total number of bond stretches in the system */
/*     ibnd    numbers of the atoms in each bond stretch */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cell.i  --  periodic boundaries using replicated cells  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xcell    length of the a-axis of the complete replicated cell */
/*     ycell    length of the b-axis of the complete replicated cell */
/*     zcell    length of the c-axis of the complete replicated cell */
/*     xcell2   half the length of the a-axis of the replicated cell */
/*     ycell2   half the length of the b-axis of the replicated cell */
/*     zcell2   half the length of the c-axis of the replicated cell */
/*     ncell    total number of cell replicates for periodic boundaries */
/*     icell    offset along axes for each replicate periodic cell */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  vdwpot.i  --  specifics of van der Waals functional form  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     abuck      value of "A" constant in Buckingham vdw potential */
/*     bbuck      value of "B" constant in Buckingham vdw potential */
/*     cbuck      value of "C" constant in Buckingham vdw potential */
/*     ghal       value of "gamma" in buffered 14-7 vdw potential */
/*     dhal       value of "delta" in buffered 14-7 vdw potential */
/*     v2scale    factor by which 1-2 vdw interactions are scaled */
/*     v3scale    factor by which 1-3 vdw interactions are scaled */
/*     v4scale    factor by which 1-4 vdw interactions are scaled */
/*     v5scale    factor by which 1-5 vdw interactions are scaled */
/*     igauss     coefficients of Gaussian fit to vdw potential */
/*     ngauss     number of Gaussians used in fit to vdw potential */
/*     vdwindex   indexing mode (atom type or class) for vdw parameters */
/*     vdwtyp     type of van der Waals potential energy function */
/*     radtyp     type of parameter (sigma or R-min) for atomic size */
/*     radsiz     atomic size provided as radius or diameter */
/*     radrule    combining rule for atomic size parameters */
/*     epsrule    combining rule for vdw well depth parameters */
/*     gausstyp   type of Gaussian fit to van der Waals potential */




/*     set arrays needed to scale connected atom interactions */

    /* Parameter adjustments */
    --zred;
    --yred;
    --xred;

    /* Function Body */
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vscale[i__ - 1] = 1.;
	iv14[i__ - 1] = 0;
    }

/*     set the coefficients for the switching function */

    switch_("VDW", (ftnlen)3);

/*     special cutoffs for very short and very long range terms */

    expmin2 = .01;
    expcut = 2.;
    expcut2 = expcut * expcut;
/* Computing 6th power */
    d__1 = expcut, d__1 *= d__1;
/* Computing 12th power */
    d__2 = expcut, d__2 *= d__2, d__2 *= d__2;
    expmerge = (vdwpot_1.abuck * exp(-vdwpot_1.bbuck / expcut) - 
	    vdwpot_1.cbuck * (d__1 * (d__1 * d__1))) / (d__2 * (d__2 * d__2));

/*     check to see if the atom of interest is a vdw site */

    nlist = 0;
    i__1 = vdw_1.nvdw;
    for (k = 1; k <= i__1; ++k) {
	if (vdw_1.ivdw[k - 1] == *iatom) {
	    ++nlist;
	    list[nlist - 1] = *iatom;
	    goto L10;
	}
    }
    return 0;
L10:

/*     determine the atoms involved via reduction factors */

    nlist = 1;
    list[nlist - 1] = *iatom;
    i__1 = couple_1.n12[*iatom - 1];
    for (k = 1; k <= i__1; ++k) {
	i__ = i12_ref(k, *iatom);
	if (vdw_1.ired[i__ - 1] == *iatom) {
	    ++nlist;
	    list[nlist - 1] = i__;
	}
    }

/*     find van der Waals Hessian elements for involved atoms */

    i__1 = nlist;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = list[ii - 1];
	iv = vdw_1.ired[i__ - 1];
	redi = vdw_1.kred[i__ - 1];
	if (i__ != iv) {
	    rediv = 1. - redi;
	    redi2 = redi * redi;
	    rediv2 = rediv * rediv;
	    rediiv = redi * rediv;
	}
	it = vdw_1.jvdw[i__ - 1];
	xi = xred[i__];
	yi = yred[i__];
	zi = zred[i__];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i12_ref(j, i__) - 1] = vdwpot_1.v2scale;
	}
	i__2 = couple_1.n13[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i13_ref(j, i__) - 1] = vdwpot_1.v3scale;
	}
	i__2 = couple_1.n14[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i14_ref(j, i__) - 1] = vdwpot_1.v4scale;
	    iv14[i14_ref(j, i__) - 1] = i__;
	}
	i__2 = couple_1.n15[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i15_ref(j, i__) - 1] = vdwpot_1.v5scale;
	}

/*     decide whether to compute the current interaction */

	i__2 = vdw_1.nvdw;
	for (kk = 1; kk <= i__2; ++kk) {
	    k = vdw_1.ivdw[kk - 1];
	    kv = vdw_1.ired[k - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = k != i__;
	    }

/*     compute the Hessian elements for this interaction */

	    if (proceed) {
		kt = vdw_1.jvdw[k - 1];
		xr = xi - xred[k];
		yr = yi - yred[k];
		zr = zi - zred[k];
		image_(&xr, &yr, &zr);
		rik2 = xr * xr + yr * yr + zr * zr;

/*     check for an interaction distance less than the cutoff */

		if (rik2 <= shunt_1.off2) {
		    fterm = 1.;
		    rv = radmin_ref(kt, it);
		    eps = epsilon_ref(kt, it);
		    if (iv14[k - 1] == i__) {
			rv = radmin4_ref(kt, it);
			eps = epsilon4_ref(kt, it);
		    } else if (radhbnd_ref(kt, it) != 0.) {
			rv = radhbnd_ref(kt, it);
			eps = epshbnd_ref(kt, it) / chgpot_1.dielec;
			if (atmtyp_1.atomic[i__ - 1] == 1) {
			    ia = i__;
			    ib = i12_ref(1, i__);
			    ic = k;
			} else {
			    ia = k;
			    ib = i12_ref(1, k);
			    ic = i__;
			}
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
			image_(&xcb, &ycb, &zcb);
			xp = ycb * zab - zcb * yab;
			yp = zcb * xab - xcb * zab;
			zp = xcb * yab - ycb * xab;
			rp = sqrt(xp * xp + yp * yp + zp * zp);
			rab2 = xab * xab + yab * yab + zab * zab;
			rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
			rcb2 = max(1e-4,rcb2);
			dot = xab * xcb + yab * ycb + zab * zcb;
			cosine = dot / sqrt(rab2 * rcb2);
/* Computing 2nd power */
			d__2 = cosine;
			sine = sqrt((d__1 = 1. - d__2 * d__2, abs(d__1)));
			rab = sqrt(rab2);
			ideal = bond_1.bl[bndlist_ref(1, ia) - 1];
			ratio = rab / ideal;
			fterm = cosine * ratio;
		    }
		    eps *= vscale[k - 1];
		    p2 = rv * rv / rik2;
		    p6 = p2 * p2 * p2;
		    rik = sqrt(rik2);
		    if (p2 <= expmin2) {
			e = 0.;
			de = 0.;
			d2e = 0.;
		    } else if (p2 <= expcut2) {
			p = sqrt(p2);
			rvterm = -vdwpot_1.bbuck / rv;
			rvterm2 = rvterm * rvterm;
			expterm = vdwpot_1.abuck * exp(-vdwpot_1.bbuck / p);
			fcbuck = fterm * vdwpot_1.cbuck * p6;
			e = eps * (expterm - fcbuck);
			de = eps * (rvterm * expterm + fcbuck * 6. / rik);
			d2e = eps * (rvterm2 * expterm - fcbuck * 42. / rik2);
		    } else {
			p12 = p6 * p6;
			e = expmerge * eps * p12;
			de = e * -12. / rik;
			d2e = e * 156. / rik2;
		    }

/*     use energy switching if near the cutoff distance */

		    if (rik2 > shunt_1.cut2) {
			rik3 = rik2 * rik;
			rik4 = rik2 * rik2;
			rik5 = rik2 * rik3;
			taper = shunt_1.c5 * rik5 + shunt_1.c4 * rik4 + 
				shunt_1.c3 * rik3 + shunt_1.c2 * rik2 + 
				shunt_1.c1 * rik + shunt_1.c0;
			dtaper = shunt_1.c5 * 5. * rik4 + shunt_1.c4 * 4. * 
				rik3 + shunt_1.c3 * 3. * rik2 + shunt_1.c2 * 
				2. * rik + shunt_1.c1;
			d2taper = shunt_1.c5 * 20. * rik3 + shunt_1.c4 * 12. *
				 rik2 + shunt_1.c3 * 6. * rik + shunt_1.c2 * 
				2.;
			d2e = e * d2taper + de * 2. * dtaper + d2e * taper;
			de = e * dtaper + de * taper;
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			de *= fgrp;
			d2e *= fgrp;
		    }

/*     get chain rule terms for van der Waals Hessian elements */

		    de /= rik;
		    d2e = (d2e - de) / rik2;
		    d2edx = d2e * xr;
		    d2edy = d2e * yr;
		    d2edz = d2e * zr;
		    term_ref(1, 1) = d2edx * xr + de;
		    term_ref(1, 2) = d2edx * yr;
		    term_ref(1, 3) = d2edx * zr;
		    term_ref(2, 1) = term_ref(1, 2);
		    term_ref(2, 2) = d2edy * yr + de;
		    term_ref(2, 3) = d2edy * zr;
		    term_ref(3, 1) = term_ref(1, 3);
		    term_ref(3, 2) = term_ref(2, 3);
		    term_ref(3, 3) = d2edz * zr + de;

/*     increment diagonal and non-diagonal Hessian elements */

		    if (i__ == *iatom) {
			if (i__ == iv && k == kv) {
			    for (j = 1; j <= 3; ++j) {
				hessx_ref(j, i__) = hessx_ref(j, i__) + 
					term_ref(1, j);
				hessy_ref(j, i__) = hessy_ref(j, i__) + 
					term_ref(2, j);
				hessz_ref(j, i__) = hessz_ref(j, i__) + 
					term_ref(3, j);
				hessx_ref(j, k) = hessx_ref(j, k) - term_ref(
					1, j);
				hessy_ref(j, k) = hessy_ref(j, k) - term_ref(
					2, j);
				hessz_ref(j, k) = hessz_ref(j, k) - term_ref(
					3, j);
			    }
			} else if (k == kv) {
			    for (j = 1; j <= 3; ++j) {
				hessx_ref(j, i__) = hessx_ref(j, i__) + 
					term_ref(1, j) * redi2;
				hessy_ref(j, i__) = hessy_ref(j, i__) + 
					term_ref(2, j) * redi2;
				hessz_ref(j, i__) = hessz_ref(j, i__) + 
					term_ref(3, j) * redi2;
				hessx_ref(j, k) = hessx_ref(j, k) - term_ref(
					1, j) * redi;
				hessy_ref(j, k) = hessy_ref(j, k) - term_ref(
					2, j) * redi;
				hessz_ref(j, k) = hessz_ref(j, k) - term_ref(
					3, j) * redi;
				hessx_ref(j, iv) = hessx_ref(j, iv) + 
					term_ref(1, j) * rediiv;
				hessy_ref(j, iv) = hessy_ref(j, iv) + 
					term_ref(2, j) * rediiv;
				hessz_ref(j, iv) = hessz_ref(j, iv) + 
					term_ref(3, j) * rediiv;
			    }
			} else if (i__ == iv) {
			    redk = vdw_1.kred[k - 1];
			    redkv = 1. - redk;
			    for (j = 1; j <= 3; ++j) {
				hessx_ref(j, i__) = hessx_ref(j, i__) + 
					term_ref(1, j);
				hessy_ref(j, i__) = hessy_ref(j, i__) + 
					term_ref(2, j);
				hessz_ref(j, i__) = hessz_ref(j, i__) + 
					term_ref(3, j);
				hessx_ref(j, k) = hessx_ref(j, k) - term_ref(
					1, j) * redk;
				hessy_ref(j, k) = hessy_ref(j, k) - term_ref(
					2, j) * redk;
				hessz_ref(j, k) = hessz_ref(j, k) - term_ref(
					3, j) * redk;
				hessx_ref(j, kv) = hessx_ref(j, kv) - 
					term_ref(1, j) * redkv;
				hessy_ref(j, kv) = hessy_ref(j, kv) - 
					term_ref(2, j) * redkv;
				hessz_ref(j, kv) = hessz_ref(j, kv) - 
					term_ref(3, j) * redkv;
			    }
			} else {
			    redk = vdw_1.kred[k - 1];
			    redkv = 1. - redk;
			    redik = redi * redk;
			    redikv = redi * redkv;
			    for (j = 1; j <= 3; ++j) {
				hessx_ref(j, i__) = hessx_ref(j, i__) + 
					term_ref(1, j) * redi2;
				hessy_ref(j, i__) = hessy_ref(j, i__) + 
					term_ref(2, j) * redi2;
				hessz_ref(j, i__) = hessz_ref(j, i__) + 
					term_ref(3, j) * redi2;
				hessx_ref(j, k) = hessx_ref(j, k) - term_ref(
					1, j) * redik;
				hessy_ref(j, k) = hessy_ref(j, k) - term_ref(
					2, j) * redik;
				hessz_ref(j, k) = hessz_ref(j, k) - term_ref(
					3, j) * redik;
				hessx_ref(j, iv) = hessx_ref(j, iv) + 
					term_ref(1, j) * rediiv;
				hessy_ref(j, iv) = hessy_ref(j, iv) + 
					term_ref(2, j) * rediiv;
				hessz_ref(j, iv) = hessz_ref(j, iv) + 
					term_ref(3, j) * rediiv;
				hessx_ref(j, kv) = hessx_ref(j, kv) - 
					term_ref(1, j) * redikv;
				hessy_ref(j, kv) = hessy_ref(j, kv) - 
					term_ref(2, j) * redikv;
				hessz_ref(j, kv) = hessz_ref(j, kv) - 
					term_ref(3, j) * redikv;
			    }
			}
		    } else if (iv == *iatom) {
			if (k == kv) {
			    for (j = 1; j <= 3; ++j) {
				hessx_ref(j, i__) = hessx_ref(j, i__) + 
					term_ref(1, j) * rediiv;
				hessy_ref(j, i__) = hessy_ref(j, i__) + 
					term_ref(2, j) * rediiv;
				hessz_ref(j, i__) = hessz_ref(j, i__) + 
					term_ref(3, j) * rediiv;
				hessx_ref(j, k) = hessx_ref(j, k) - term_ref(
					1, j) * rediv;
				hessy_ref(j, k) = hessy_ref(j, k) - term_ref(
					2, j) * rediv;
				hessz_ref(j, k) = hessz_ref(j, k) - term_ref(
					3, j) * rediv;
				hessx_ref(j, iv) = hessx_ref(j, iv) + 
					term_ref(1, j) * rediv2;
				hessy_ref(j, iv) = hessy_ref(j, iv) + 
					term_ref(2, j) * rediv2;
				hessz_ref(j, iv) = hessz_ref(j, iv) + 
					term_ref(3, j) * rediv2;
			    }
			} else {
			    redk = vdw_1.kred[k - 1];
			    redkv = 1. - redk;
			    redivk = rediv * redk;
			    redivkv = rediv * redkv;
			    for (j = 1; j <= 3; ++j) {
				hessx_ref(j, i__) = hessx_ref(j, i__) + 
					term_ref(1, j) * rediiv;
				hessy_ref(j, i__) = hessy_ref(j, i__) + 
					term_ref(2, j) * rediiv;
				hessz_ref(j, i__) = hessz_ref(j, i__) + 
					term_ref(3, j) * rediiv;
				hessx_ref(j, k) = hessx_ref(j, k) - term_ref(
					1, j) * redivk;
				hessy_ref(j, k) = hessy_ref(j, k) - term_ref(
					2, j) * redivk;
				hessz_ref(j, k) = hessz_ref(j, k) - term_ref(
					3, j) * redivk;
				hessx_ref(j, iv) = hessx_ref(j, iv) + 
					term_ref(1, j) * rediv2;
				hessy_ref(j, iv) = hessy_ref(j, iv) + 
					term_ref(2, j) * rediv2;
				hessz_ref(j, iv) = hessz_ref(j, iv) + 
					term_ref(3, j) * rediv2;
				hessx_ref(j, kv) = hessx_ref(j, kv) - 
					term_ref(1, j) * redivkv;
				hessy_ref(j, kv) = hessy_ref(j, kv) - 
					term_ref(2, j) * redivkv;
				hessz_ref(j, kv) = hessz_ref(j, kv) - 
					term_ref(3, j) * redivkv;
			    }
			}
		    }
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i12_ref(j, i__) - 1] = 1.;
	}
	i__2 = couple_1.n13[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i13_ref(j, i__) - 1] = 1.;
	}
	i__2 = couple_1.n14[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i14_ref(j, i__) - 1] = 1.;
	}
	i__2 = couple_1.n15[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i15_ref(j, i__) - 1] = 1.;
	}
    }

/*     for periodic boundary conditions with large cutoffs */
/*     neighbors must be found by the replicates method */

    if (! bound_1.use_replica__) {
	return 0;
    }

/*     calculate interaction energy with other unit cells */

    i__1 = nlist;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = list[ii - 1];
	iv = vdw_1.ired[i__ - 1];
	redi = vdw_1.kred[i__ - 1];
	if (i__ != iv) {
	    rediv = 1. - redi;
	    redi2 = redi * redi;
	    rediv2 = rediv * rediv;
	    rediiv = redi * rediv;
	}
	it = vdw_1.jvdw[i__ - 1];
	xi = xred[i__];
	yi = yred[i__];
	zi = zred[i__];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i12_ref(j, i__) - 1] = vdwpot_1.v2scale;
	}
	i__2 = couple_1.n13[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i13_ref(j, i__) - 1] = vdwpot_1.v3scale;
	}
	i__2 = couple_1.n14[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i14_ref(j, i__) - 1] = vdwpot_1.v4scale;
	    iv14[i14_ref(j, i__) - 1] = i__;
	}
	i__2 = couple_1.n15[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i15_ref(j, i__) - 1] = vdwpot_1.v5scale;
	}

/*     decide whether to compute the current interaction */

	i__2 = vdw_1.nvdw;
	for (kk = 1; kk <= i__2; ++kk) {
	    k = vdw_1.ivdw[kk - 1];
	    kv = vdw_1.ired[k - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }

/*     compute the Hessian elements for this interaction */

	    if (proceed) {
		kt = vdw_1.jvdw[k - 1];
		i__3 = cell_1.ncell;
		for (jcell = 1; jcell <= i__3; ++jcell) {
		    xr = xi - xred[k];
		    yr = yi - yred[k];
		    zr = zi - zred[k];
		    imager_(&xr, &yr, &zr, &jcell);
		    rik2 = xr * xr + yr * yr + zr * zr;

/*     check for an interaction distance less than the cutoff */

		    if (rik2 <= shunt_1.off2) {
			fterm = 1.;
			rv = radmin_ref(kt, it);
			eps = epsilon_ref(kt, it);
			if (radhbnd_ref(kt, it) != 0.) {
			    rv = radhbnd_ref(kt, it);
			    eps = epshbnd_ref(kt, it) / chgpot_1.dielec;
			    if (atmtyp_1.atomic[i__ - 1] == 1) {
				ia = i__;
				ib = i12_ref(1, i__);
				ic = k;
			    } else {
				ia = k;
				ib = i12_ref(1, k);
				ic = i__;
			    }
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
			    imager_(&xcb, &ycb, &zcb, &jcell);
			    xp = ycb * zab - zcb * yab;
			    yp = zcb * xab - xcb * zab;
			    zp = xcb * yab - ycb * xab;
			    rp = sqrt(xp * xp + yp * yp + zp * zp);
			    rab2 = xab * xab + yab * yab + zab * zab;
			    rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
			    rcb2 = max(1e-4,rcb2);
			    dot = xab * xcb + yab * ycb + zab * zcb;
			    cosine = dot / sqrt(rab2 * rcb2);
/* Computing 2nd power */
			    d__2 = cosine;
			    sine = sqrt((d__1 = 1. - d__2 * d__2, abs(d__1)));
			    rab = sqrt(rab2);
			    ideal = bond_1.bl[bndlist_ref(1, ia) - 1];
			    ratio = rab / ideal;
			    fterm = cosine * ratio;
			}
			if (bound_1.use_polymer__) {
			    if (rik2 <= bound_1.polycut2) {
				if (iv14[k - 1] == i__) {
				    fterm = 1.;
				    rv = radmin4_ref(kt, it);
				    eps = epsilon4_ref(kt, it);
				}
				eps *= vscale[k - 1];
			    }
			}
			p2 = rv * rv / rik2;
			p6 = p2 * p2 * p2;
			rik = sqrt(rik2);
			if (p2 <= expmin2) {
			    e = 0.;
			    de = 0.;
			    d2e = 0.;
			} else if (p2 <= expcut2) {
			    p = sqrt(p2);
			    rvterm = -vdwpot_1.bbuck / rv;
			    rvterm2 = rvterm * rvterm;
			    expterm = vdwpot_1.abuck * exp(-vdwpot_1.bbuck / 
				    p);
			    fcbuck = fterm * vdwpot_1.cbuck * p6;
			    e = eps * (expterm - fcbuck);
			    de = eps * (rvterm * expterm + fcbuck * 6. / rik);
			    d2e = eps * (rvterm2 * expterm - fcbuck * 42. / 
				    rik2);
			} else {
			    p12 = p6 * p6;
			    e = expmerge * eps * p12;
			    de = e * -12. / rik;
			    d2e = e * 156. / rik2;
			}

/*     use energy switching if near the cutoff distance */

			if (rik2 > shunt_1.cut2) {
			    rik3 = rik2 * rik;
			    rik4 = rik2 * rik2;
			    rik5 = rik2 * rik3;
			    taper = shunt_1.c5 * rik5 + shunt_1.c4 * rik4 + 
				    shunt_1.c3 * rik3 + shunt_1.c2 * rik2 + 
				    shunt_1.c1 * rik + shunt_1.c0;
			    dtaper = shunt_1.c5 * 5. * rik4 + shunt_1.c4 * 4. 
				    * rik3 + shunt_1.c3 * 3. * rik2 + 
				    shunt_1.c2 * 2. * rik + shunt_1.c1;
			    d2taper = shunt_1.c5 * 20. * rik3 + shunt_1.c4 * 
				    12. * rik2 + shunt_1.c3 * 6. * rik + 
				    shunt_1.c2 * 2.;
			    d2e = e * d2taper + de * 2. * dtaper + d2e * 
				    taper;
			    de = e * dtaper + de * taper;
			}

/*     scale the interaction based on its group membership */

			if (group_1.use_group__) {
			    de *= fgrp;
			    d2e *= fgrp;
			}

/*     get chain rule terms for van der Waals Hessian elements */

			de /= rik;
			d2e = (d2e - de) / rik2;
			d2edx = d2e * xr;
			d2edy = d2e * yr;
			d2edz = d2e * zr;
			term_ref(1, 1) = d2edx * xr + de;
			term_ref(1, 2) = d2edx * yr;
			term_ref(1, 3) = d2edx * zr;
			term_ref(2, 1) = term_ref(1, 2);
			term_ref(2, 2) = d2edy * yr + de;
			term_ref(2, 3) = d2edy * zr;
			term_ref(3, 1) = term_ref(1, 3);
			term_ref(3, 2) = term_ref(2, 3);
			term_ref(3, 3) = d2edz * zr + de;

/*     increment diagonal and non-diagonal Hessian elements */

			if (i__ == *iatom) {
			    if (i__ == iv && k == kv) {
				for (j = 1; j <= 3; ++j) {
				    hessx_ref(j, i__) = hessx_ref(j, i__) + 
					    term_ref(1, j);
				    hessy_ref(j, i__) = hessy_ref(j, i__) + 
					    term_ref(2, j);
				    hessz_ref(j, i__) = hessz_ref(j, i__) + 
					    term_ref(3, j);
				    hessx_ref(j, k) = hessx_ref(j, k) - 
					    term_ref(1, j);
				    hessy_ref(j, k) = hessy_ref(j, k) - 
					    term_ref(2, j);
				    hessz_ref(j, k) = hessz_ref(j, k) - 
					    term_ref(3, j);
				}
			    } else if (k == kv) {
				for (j = 1; j <= 3; ++j) {
				    hessx_ref(j, i__) = hessx_ref(j, i__) + 
					    term_ref(1, j) * redi2;
				    hessy_ref(j, i__) = hessy_ref(j, i__) + 
					    term_ref(2, j) * redi2;
				    hessz_ref(j, i__) = hessz_ref(j, i__) + 
					    term_ref(3, j) * redi2;
				    hessx_ref(j, k) = hessx_ref(j, k) - 
					    term_ref(1, j) * redi;
				    hessy_ref(j, k) = hessy_ref(j, k) - 
					    term_ref(2, j) * redi;
				    hessz_ref(j, k) = hessz_ref(j, k) - 
					    term_ref(3, j) * redi;
				    hessx_ref(j, iv) = hessx_ref(j, iv) + 
					    term_ref(1, j) * rediiv;
				    hessy_ref(j, iv) = hessy_ref(j, iv) + 
					    term_ref(2, j) * rediiv;
				    hessz_ref(j, iv) = hessz_ref(j, iv) + 
					    term_ref(3, j) * rediiv;
				}
			    } else if (i__ == iv) {
				redk = vdw_1.kred[k - 1];
				redkv = 1. - redk;
				for (j = 1; j <= 3; ++j) {
				    hessx_ref(j, i__) = hessx_ref(j, i__) + 
					    term_ref(1, j);
				    hessy_ref(j, i__) = hessy_ref(j, i__) + 
					    term_ref(2, j);
				    hessz_ref(j, i__) = hessz_ref(j, i__) + 
					    term_ref(3, j);
				    hessx_ref(j, k) = hessx_ref(j, k) - 
					    term_ref(1, j) * redk;
				    hessy_ref(j, k) = hessy_ref(j, k) - 
					    term_ref(2, j) * redk;
				    hessz_ref(j, k) = hessz_ref(j, k) - 
					    term_ref(3, j) * redk;
				    hessx_ref(j, kv) = hessx_ref(j, kv) - 
					    term_ref(1, j) * redkv;
				    hessy_ref(j, kv) = hessy_ref(j, kv) - 
					    term_ref(2, j) * redkv;
				    hessz_ref(j, kv) = hessz_ref(j, kv) - 
					    term_ref(3, j) * redkv;
				}
			    } else {
				redk = vdw_1.kred[k - 1];
				redkv = 1. - redk;
				redik = redi * redk;
				redikv = redi * redkv;
				for (j = 1; j <= 3; ++j) {
				    hessx_ref(j, i__) = hessx_ref(j, i__) + 
					    term_ref(1, j) * redi2;
				    hessy_ref(j, i__) = hessy_ref(j, i__) + 
					    term_ref(2, j) * redi2;
				    hessz_ref(j, i__) = hessz_ref(j, i__) + 
					    term_ref(3, j) * redi2;
				    hessx_ref(j, k) = hessx_ref(j, k) - 
					    term_ref(1, j) * redik;
				    hessy_ref(j, k) = hessy_ref(j, k) - 
					    term_ref(2, j) * redik;
				    hessz_ref(j, k) = hessz_ref(j, k) - 
					    term_ref(3, j) * redik;
				    hessx_ref(j, iv) = hessx_ref(j, iv) + 
					    term_ref(1, j) * rediiv;
				    hessy_ref(j, iv) = hessy_ref(j, iv) + 
					    term_ref(2, j) * rediiv;
				    hessz_ref(j, iv) = hessz_ref(j, iv) + 
					    term_ref(3, j) * rediiv;
				    hessx_ref(j, kv) = hessx_ref(j, kv) - 
					    term_ref(1, j) * redikv;
				    hessy_ref(j, kv) = hessy_ref(j, kv) - 
					    term_ref(2, j) * redikv;
				    hessz_ref(j, kv) = hessz_ref(j, kv) - 
					    term_ref(3, j) * redikv;
				}
			    }
			} else if (iv == *iatom) {
			    if (k == kv) {
				for (j = 1; j <= 3; ++j) {
				    hessx_ref(j, i__) = hessx_ref(j, i__) + 
					    term_ref(1, j) * rediiv;
				    hessy_ref(j, i__) = hessy_ref(j, i__) + 
					    term_ref(2, j) * rediiv;
				    hessz_ref(j, i__) = hessz_ref(j, i__) + 
					    term_ref(3, j) * rediiv;
				    hessx_ref(j, k) = hessx_ref(j, k) - 
					    term_ref(1, j) * rediv;
				    hessy_ref(j, k) = hessy_ref(j, k) - 
					    term_ref(2, j) * rediv;
				    hessz_ref(j, k) = hessz_ref(j, k) - 
					    term_ref(3, j) * rediv;
				    hessx_ref(j, iv) = hessx_ref(j, iv) + 
					    term_ref(1, j) * rediv2;
				    hessy_ref(j, iv) = hessy_ref(j, iv) + 
					    term_ref(2, j) * rediv2;
				    hessz_ref(j, iv) = hessz_ref(j, iv) + 
					    term_ref(3, j) * rediv2;
				}
			    } else {
				redk = vdw_1.kred[k - 1];
				redkv = 1. - redk;
				redivk = rediv * redk;
				redivkv = rediv * redkv;
				for (j = 1; j <= 3; ++j) {
				    hessx_ref(j, i__) = hessx_ref(j, i__) + 
					    term_ref(1, j) * rediiv;
				    hessy_ref(j, i__) = hessy_ref(j, i__) + 
					    term_ref(2, j) * rediiv;
				    hessz_ref(j, i__) = hessz_ref(j, i__) + 
					    term_ref(3, j) * rediiv;
				    hessx_ref(j, k) = hessx_ref(j, k) - 
					    term_ref(1, j) * redivk;
				    hessy_ref(j, k) = hessy_ref(j, k) - 
					    term_ref(2, j) * redivk;
				    hessz_ref(j, k) = hessz_ref(j, k) - 
					    term_ref(3, j) * redivk;
				    hessx_ref(j, iv) = hessx_ref(j, iv) + 
					    term_ref(1, j) * rediv2;
				    hessy_ref(j, iv) = hessy_ref(j, iv) + 
					    term_ref(2, j) * rediv2;
				    hessz_ref(j, iv) = hessz_ref(j, iv) + 
					    term_ref(3, j) * rediv2;
				    hessx_ref(j, kv) = hessx_ref(j, kv) - 
					    term_ref(1, j) * redivkv;
				    hessy_ref(j, kv) = hessy_ref(j, kv) - 
					    term_ref(2, j) * redivkv;
				    hessz_ref(j, kv) = hessz_ref(j, kv) - 
					    term_ref(3, j) * redivkv;
				}
			    }
			}
		    }
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i12_ref(j, i__) - 1] = 1.;
	}
	i__2 = couple_1.n13[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i13_ref(j, i__) - 1] = 1.;
	}
	i__2 = couple_1.n14[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i14_ref(j, i__) - 1] = 1.;
	}
	i__2 = couple_1.n15[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    vscale[i15_ref(j, i__) - 1] = 1.;
	}
    }
    return 0;
} /* emm3hb2_ */

#undef epsilon_ref
#undef bndlist_ref
#undef epshbnd_ref
#undef radhbnd_ref
#undef radmin4_ref
#undef radmin_ref
#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef term_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref
#undef epsilon4_ref


