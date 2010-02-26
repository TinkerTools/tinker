/* egauss.f -- translated by f2c (version 20050501).
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
    doublereal m2[25000], deform, difft, diffv, diffc;
    logical use_smooth__, use_dem__, use_gda__, use_tophat__, use_stophat__;
} warp_;

#define warp_1 warp_

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
    doublereal xcell, ycell, zcell, xcell2, ycell2, zcell2;
    integer ncell, icell[30000]	/* was [3][10000] */;
} cell_;

#define cell_1 cell_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

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
    doublereal radmin[1000000]	/* was [1000][1000] */, epsilon[1000000]	
	    /* was [1000][1000] */, radmin4[1000000]	/* was [1000][1000] */
	    , epsilon4[1000000]	/* was [1000][1000] */, radhbnd[1000000]	
	    /* was [1000][1000] */, epshbnd[1000000]	/* was [1000][1000] */
	    , kred[25000];
    integer ired[25000], nvdw, ivdw[25000], jvdw[25000];
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
    integer nlight, kbx[25000], kby[25000], kbz[25000], kex[25000], key[25000]
	    , kez[25000], locx[200000], locy[200000], locz[200000], rgx[
	    200000], rgy[200000], rgz[200000];
} light_;

#define light_1 light_

struct {
    doublereal lbuffer, lbuf2, vbuf2, cbuf2, mbuf2;
    integer nvlst[25000], vlst[45000000]	/* was [1800][25000] */, 
	    nelst[25000], elst[30000000]	/* was [1200][25000] */;
    logical dovlst, doclst, domlst;
} neigh_;

#define neigh_1 neigh_

/* Table of constant values */

static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine egauss  --  Gaussian van der Waals energy  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "egauss" calculates the Gaussian expansion van der Waals energy */


/* Subroutine */ int egauss_(void)
{
    extern /* Subroutine */ int egauss0a_(void), egauss0b_(void), egauss0c_(
	    void), egauss0d_(void);



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




/*     choose the method for summing over pairwise interactions */



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
	egauss0d_();
    } else if (cutoff_1.use_vlist__) {
	egauss0c_();
    } else if (cutoff_1.use_lights__) {
	egauss0b_();
    } else {
	egauss0a_();
    }
    return 0;
} /* egauss_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine egauss0a  --  double loop Gaussian vdw energy  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "egauss0a" calculates the Gaussian expansion van der Waals */
/*     energy using a pairwise double loop */


/* Subroutine */ int egauss0a_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal a[10], b[10], e;
    static integer i__, j, k, m, ii, kk, it, iv, kt, kv;
    static doublereal xi, yi, zi, xr, yr, zr;
    static integer iv14[25000];
    static doublereal rdn, rik, eps, rad2, rik2, rik3, rik4, rik5, fgrp, xred[
	    25000], yred[25000], zred[25000];
    static logical usei;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal taper;
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal vscale[25000];
    extern /* Subroutine */ int switch_(char *, ftnlen);
    static doublereal expcut;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;
    static doublereal expterm;


#define epsilon4_ref(a_1,a_2) vdw_1.epsilon4[(a_2)*1000 + a_1 - 1001]
#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define radmin_ref(a_1,a_2) vdw_1.radmin[(a_2)*1000 + a_1 - 1001]
#define igauss_ref(a_1,a_2) vdwpot_1.igauss[(a_2)*2 + a_1 - 3]
#define radmin4_ref(a_1,a_2) vdw_1.radmin4[(a_2)*1000 + a_1 - 1001]
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




/*     zero out the van der Waals energy contribution */

    energi_1.ev = 0.;

/*     set arrays needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vscale[i__ - 1] = 1.;
	iv14[i__ - 1] = 0;
    }

/*     set cutoff distances and switching function coefficients */

    switch_("VDW", (ftnlen)3);
    expcut = -50.;

/*     apply any reduction factor to the atomic coordinates */

    i__1 = vdw_1.nvdw;
    for (k = 1; k <= i__1; ++k) {
	i__ = vdw_1.ivdw[k - 1];
	iv = vdw_1.ired[i__ - 1];
	rdn = vdw_1.kred[i__ - 1];
	xred[i__ - 1] = rdn * (atoms_1.x[i__ - 1] - atoms_1.x[iv - 1]) + 
		atoms_1.x[iv - 1];
	yred[i__ - 1] = rdn * (atoms_1.y[i__ - 1] - atoms_1.y[iv - 1]) + 
		atoms_1.y[iv - 1];
	zred[i__ - 1] = rdn * (atoms_1.z__[i__ - 1] - atoms_1.z__[iv - 1]) + 
		atoms_1.z__[iv - 1];
    }

/*     find the van der Waals energy via double loop search */

    i__1 = vdw_1.nvdw - 1;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = vdw_1.ivdw[ii - 1];
	iv = vdw_1.ired[i__ - 1];
	it = vdw_1.jvdw[i__ - 1];
	xi = xred[i__ - 1];
	yi = yred[i__ - 1];
	zi = zred[i__ - 1];
	usei = usage_1.use[i__ - 1] || usage_1.use[iv - 1];

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
	for (kk = ii + 1; kk <= i__2; ++kk) {
	    k = vdw_1.ivdw[kk - 1];
	    kv = vdw_1.ired[k - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1] || usage_1.use[kv - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		kt = vdw_1.jvdw[k - 1];
		xr = xi - xred[k - 1];
		yr = yi - yred[k - 1];
		zr = zi - zred[k - 1];
		image_(&xr, &yr, &zr);
		rik2 = xr * xr + yr * yr + zr * zr;

/*     check for an interaction distance less than the cutoff */

		if (rik2 <= shunt_1.off2) {
/* Computing 2nd power */
		    d__1 = radmin_ref(kt, it);
		    rad2 = d__1 * d__1;
		    eps = epsilon_ref(kt, it);
		    if (iv14[k - 1] == i__) {
/* Computing 2nd power */
			d__1 = radmin4_ref(kt, it);
			rad2 = d__1 * d__1;
			eps = epsilon4_ref(kt, it);
		    }
		    eps *= vscale[k - 1];
		    i__3 = vdwpot_1.ngauss;
		    for (j = 1; j <= i__3; ++j) {
			a[j - 1] = igauss_ref(1, j) * eps;
			b[j - 1] = igauss_ref(2, j) / rad2;
		    }
		    e = 0.;
		    i__3 = vdwpot_1.ngauss;
		    for (j = 1; j <= i__3; ++j) {
			expterm = -b[j - 1] * rik2;
			if (expterm > expcut) {
			    e += a[j - 1] * exp(expterm);
			}
		    }

/*     use energy switching if near the cutoff distance */

		    if (rik2 > shunt_1.cut2) {
			rik = sqrt(rik2);
			rik3 = rik2 * rik;
			rik4 = rik2 * rik2;
			rik5 = rik2 * rik3;
			taper = shunt_1.c5 * rik5 + shunt_1.c4 * rik4 + 
				shunt_1.c3 * rik3 + shunt_1.c2 * rik2 + 
				shunt_1.c1 * rik + shunt_1.c0;
			e *= taper;
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
		    }

/*     increment the overall van der Waals energy components */

		    energi_1.ev += e;
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

    i__1 = vdw_1.nvdw;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = vdw_1.ivdw[ii - 1];
	iv = vdw_1.ired[i__ - 1];
	it = vdw_1.jvdw[i__ - 1];
	xi = xred[i__ - 1];
	yi = yred[i__ - 1];
	zi = zred[i__ - 1];
	usei = usage_1.use[i__ - 1] || usage_1.use[iv - 1];

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
	for (kk = ii; kk <= i__2; ++kk) {
	    k = vdw_1.ivdw[kk - 1];
	    kv = vdw_1.ired[k - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1] || usage_1.use[kv - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		kt = vdw_1.jvdw[k - 1];
		i__3 = cell_1.ncell;
		for (m = 1; m <= i__3; ++m) {
		    xr = xi - xred[k - 1];
		    yr = yi - yred[k - 1];
		    zr = zi - zred[k - 1];
		    imager_(&xr, &yr, &zr, &m);
		    rik2 = xr * xr + yr * yr + zr * zr;

/*     check for an interaction distance less than the cutoff */

		    if (rik2 <= shunt_1.off2) {
/* Computing 2nd power */
			d__1 = radmin_ref(kt, it);
			rad2 = d__1 * d__1;
			eps = epsilon_ref(kt, it);
			if (bound_1.use_polymer__) {
			    if (rik2 <= bound_1.polycut2) {
				if (iv14[k - 1] == i__) {
/* Computing 2nd power */
				    d__1 = radmin4_ref(kt, it);
				    rad2 = d__1 * d__1;
				    eps = epsilon4_ref(kt, it);
				}
				eps *= vscale[k - 1];
			    }
			}
			i__4 = vdwpot_1.ngauss;
			for (j = 1; j <= i__4; ++j) {
			    a[j - 1] = igauss_ref(1, j) * eps;
			    b[j - 1] = igauss_ref(2, j) / rad2;
			}
			e = 0.;
			i__4 = vdwpot_1.ngauss;
			for (j = 1; j <= i__4; ++j) {
			    expterm = -b[j - 1] * rik2;
			    if (expterm > expcut) {
				e += a[j - 1] * exp(expterm);
			    }
			}

/*     use energy switching if near the cutoff distance */

			if (rik2 > shunt_1.cut2) {
			    rik = sqrt(rik2);
			    rik3 = rik2 * rik;
			    rik4 = rik2 * rik2;
			    rik5 = rik2 * rik3;
			    taper = shunt_1.c5 * rik5 + shunt_1.c4 * rik4 + 
				    shunt_1.c3 * rik3 + shunt_1.c2 * rik2 + 
				    shunt_1.c1 * rik + shunt_1.c0;
			    e *= taper;
			}

/*     scale the interaction based on its group membership */

			if (group_1.use_group__) {
			    e *= fgrp;
			}

/*     increment the overall van der Waals energy component; */
/*     interaction of an atom with its own image counts half */

			if (i__ == k) {
			    e *= .5;
			}
			energi_1.ev += e;
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
} /* egauss0a_ */

#undef epsilon_ref
#undef radmin4_ref
#undef igauss_ref
#undef radmin_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref
#undef epsilon4_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine egauss0b  --  Gaussian vdw energy via lights  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "egauss0b" calculates the Gaussian expansion van der Waals energy */
/*     using the method of lights */


/* Subroutine */ int egauss0b_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), exp(doublereal), sqrt(
	    doublereal);

    /* Local variables */
    static doublereal a[10], b[10], e;
    static integer i__, j, k, m, ii, kk, it, iv, kt, kv;
    static doublereal xi, yi, zi, xr, yr, zr;
    static integer iv14[25000], kgy, kgz;
    static doublereal eps, rdn, rik, rad2, rik2, rik3, rik4, rik5, fgrp, xred[
	    25000], yred[25000];
    static integer stop;
    static doublereal zred[25000];
    static logical usei;
    static doublereal taper;
    static logical prime;
    static integer start;
    static doublereal xsort[200000], ysort[200000], zsort[200000], vscale[
	    25000];
    static logical repeat;
    extern /* Subroutine */ int lights_(doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *), switch_(char *, ftnlen);
    static doublereal expcut;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;
    static doublereal expterm;


#define epsilon4_ref(a_1,a_2) vdw_1.epsilon4[(a_2)*1000 + a_1 - 1001]
#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define radmin_ref(a_1,a_2) vdw_1.radmin[(a_2)*1000 + a_1 - 1001]
#define igauss_ref(a_1,a_2) vdwpot_1.igauss[(a_2)*2 + a_1 - 3]
#define radmin4_ref(a_1,a_2) vdw_1.radmin4[(a_2)*1000 + a_1 - 1001]
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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  light.i  --  indices for method of lights pair neighbors  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     nlight  total number of sites for method of lights calculation */
/*     kbx     low index of neighbors of each site in the x-sorted list */
/*     kby     low index of neighbors of each site in the y-sorted list */
/*     kbz     low index of neighbors of each site in the z-sorted list */
/*     kex     high index of neighbors of each site in the x-sorted list */
/*     key     high index of neighbors of each site in the y-sorted list */
/*     kez     high index of neighbors of each site in the z-sorted list */
/*     locx    pointer from x-sorted list into original interaction list */
/*     locy    pointer from y-sorted list into original interaction list */
/*     locz    pointer from z-sorted list into original interaction list */
/*     rgx     pointer from original interaction list into x-sorted list */
/*     rgy     pointer from original interaction list into y-sorted list */
/*     rgz     pointer from original interaction list into z-sorted list */




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




/*     zero out the van der Waals energy contribution */

    energi_1.ev = 0.;

/*     set arrays needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vscale[i__ - 1] = 1.;
	iv14[i__ - 1] = 0;
    }

/*     set cutoff distances and switching function coefficients */

    switch_("VDW", (ftnlen)3);
    expcut = -50.;

/*     apply any reduction factor to the atomic coordinates */

    i__1 = vdw_1.nvdw;
    for (j = 1; j <= i__1; ++j) {
	i__ = vdw_1.ivdw[j - 1];
	iv = vdw_1.ired[i__ - 1];
	rdn = vdw_1.kred[i__ - 1];
	xred[j - 1] = rdn * (atoms_1.x[i__ - 1] - atoms_1.x[iv - 1]) + 
		atoms_1.x[iv - 1];
	yred[j - 1] = rdn * (atoms_1.y[i__ - 1] - atoms_1.y[iv - 1]) + 
		atoms_1.y[iv - 1];
	zred[j - 1] = rdn * (atoms_1.z__[i__ - 1] - atoms_1.z__[iv - 1]) + 
		atoms_1.z__[iv - 1];
    }

/*     transfer the interaction site coordinates to sorting arrays */

    i__1 = vdw_1.nvdw;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xsort[i__ - 1] = xred[i__ - 1];
	ysort[i__ - 1] = yred[i__ - 1];
	zsort[i__ - 1] = zred[i__ - 1];
    }

/*     use the method of lights to generate neighbors */

    lights_(&shunt_1.off, &vdw_1.nvdw, xsort, ysort, zsort);

/*     loop over all atoms computing the interactions */

    i__1 = vdw_1.nvdw;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = vdw_1.ivdw[ii - 1];
	iv = vdw_1.ired[i__ - 1];
	it = vdw_1.jvdw[i__ - 1];
	xi = xsort[light_1.rgx[ii - 1] - 1];
	yi = ysort[light_1.rgy[ii - 1] - 1];
	zi = zsort[light_1.rgz[ii - 1] - 1];
	usei = usage_1.use[i__ - 1] || usage_1.use[iv - 1];

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

/*     loop over method of lights neighbors of current atom */

	if (light_1.kbx[ii - 1] <= light_1.kex[ii - 1]) {
	    repeat = FALSE_;
	    start = light_1.kbx[ii - 1] + 1;
	    stop = light_1.kex[ii - 1];
	} else {
	    repeat = TRUE_;
	    start = 1;
	    stop = light_1.kex[ii - 1];
	}
L10:
	i__2 = stop;
	for (m = start; m <= i__2; ++m) {
	    kk = light_1.locx[m - 1];
	    kgy = light_1.rgy[kk - 1];
	    if (light_1.kby[ii - 1] <= light_1.key[ii - 1]) {
		if (kgy < light_1.kby[ii - 1] || kgy > light_1.key[ii - 1]) {
		    goto L20;
		}
	    } else {
		if (kgy < light_1.kby[ii - 1] && kgy > light_1.key[ii - 1]) {
		    goto L20;
		}
	    }
	    kgz = light_1.rgz[kk - 1];
	    if (light_1.kbz[ii - 1] <= light_1.kez[ii - 1]) {
		if (kgz < light_1.kbz[ii - 1] || kgz > light_1.kez[ii - 1]) {
		    goto L20;
		}
	    } else {
		if (kgz < light_1.kbz[ii - 1] && kgz > light_1.kez[ii - 1]) {
		    goto L20;
		}
	    }
	    k = vdw_1.ivdw[kk - (kk - 1) / vdw_1.nvdw * vdw_1.nvdw - 1];
	    kv = vdw_1.ired[k - 1];
	    prime = kk <= vdw_1.nvdw;

/*     decide whether to compute the current interaction */

	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1] || usage_1.use[kv - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		kt = vdw_1.jvdw[k - 1];
		xr = xi - xsort[m - 1];
		yr = yi - ysort[kgy - 1];
		zr = zi - zsort[kgz - 1];
		if (bound_1.use_bounds__) {
		    if (abs(xr) > cell_1.xcell2) {
			xr -= d_sign(&cell_1.xcell, &xr);
		    }
		    if (abs(yr) > cell_1.ycell2) {
			yr -= d_sign(&cell_1.ycell, &yr);
		    }
		    if (abs(zr) > cell_1.zcell2) {
			zr -= d_sign(&cell_1.zcell, &zr);
		    }
		    if (boxes_1.monoclinic) {
			xr += zr * boxes_1.beta_cos__;
			zr *= boxes_1.beta_sin__;
		    } else if (boxes_1.triclinic) {
			xr = xr + yr * boxes_1.gamma_cos__ + zr * 
				boxes_1.beta_cos__;
			yr = yr * boxes_1.gamma_sin__ + zr * 
				boxes_1.beta_term__;
			zr *= boxes_1.gamma_term__;
		    }
		}
		rik2 = xr * xr + yr * yr + zr * zr;

/*     check for an interaction distance less than the cutoff */

		if (rik2 <= shunt_1.off2) {
/* Computing 2nd power */
		    d__1 = radmin_ref(kt, it);
		    rad2 = d__1 * d__1;
		    eps = epsilon_ref(kt, it);
		    if (prime) {
			if (iv14[k - 1] == i__) {
/* Computing 2nd power */
			    d__1 = radmin4_ref(kt, it);
			    rad2 = d__1 * d__1;
			    eps = epsilon4_ref(kt, it);
			}
			eps *= vscale[k - 1];
		    }
		    i__3 = vdwpot_1.ngauss;
		    for (j = 1; j <= i__3; ++j) {
			a[j - 1] = igauss_ref(1, j) * eps;
			b[j - 1] = igauss_ref(2, j) / rad2;
		    }
		    e = 0.;
		    i__3 = vdwpot_1.ngauss;
		    for (j = 1; j <= i__3; ++j) {
			expterm = -b[j - 1] * rik2;
			if (expterm > expcut) {
			    e += a[j - 1] * exp(expterm);
			}
		    }

/*     use energy switching if near the cutoff distance */

		    if (rik2 > shunt_1.cut2) {
			rik = sqrt(rik2);
			rik3 = rik2 * rik;
			rik4 = rik2 * rik2;
			rik5 = rik2 * rik3;
			taper = shunt_1.c5 * rik5 + shunt_1.c4 * rik4 + 
				shunt_1.c3 * rik3 + shunt_1.c2 * rik2 + 
				shunt_1.c1 * rik + shunt_1.c0;
			e *= taper;
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
		    }

/*     increment the overall van der Waals energy component */

		    energi_1.ev += e;
		}
	    }
L20:
	    ;
	}
	if (repeat) {
	    repeat = FALSE_;
	    start = light_1.kbx[ii - 1] + 1;
	    stop = light_1.nlight;
	    goto L10;
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
} /* egauss0b_ */

#undef epsilon_ref
#undef radmin4_ref
#undef igauss_ref
#undef radmin_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref
#undef epsilon4_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine egauss0c  --  Gaussian vdw energy via list  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "egauss0c" calculates the Gaussian expansion van der Waals */
/*     energy using a pairwise neighbor list */


/* Subroutine */ int egauss0c_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal a[10], b[10], e;
    static integer i__, j, k, ii, kk, it, iv, kt, kv;
    static doublereal xi, yi, zi, xr, yr, zr;
    static integer iv14[25000];
    static doublereal rdn, rik, eps, rad2, rik2, rik3, rik4, rik5, fgrp, xred[
	    25000], yred[25000], zred[25000];
    static logical usei;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal taper, vscale[25000], expcut;
    extern /* Subroutine */ int switch_(char *, ftnlen), groups_(logical *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *);
    static logical proceed;
    static doublereal expterm;


#define epsilon4_ref(a_1,a_2) vdw_1.epsilon4[(a_2)*1000 + a_1 - 1001]
#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define vlst_ref(a_1,a_2) neigh_1.vlst[(a_2)*1800 + a_1 - 1801]
#define radmin_ref(a_1,a_2) vdw_1.radmin[(a_2)*1000 + a_1 - 1001]
#define igauss_ref(a_1,a_2) vdwpot_1.igauss[(a_2)*2 + a_1 - 3]
#define radmin4_ref(a_1,a_2) vdw_1.radmin4[(a_2)*1000 + a_1 - 1001]
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




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  neigh.i  --  pairwise neighbor list indices and storage  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     lbuffer     width of the neighbor list buffer region */
/*     lbuf2       motion squared needed to trigger list rebuild */
/*     vbuf2       square of vdw cutoff plus neighbor list buffer */
/*     cbuf2       square of charge cutoff plus neighbor list buffer */
/*     mbuf2       square of multipole cutoff plus neighbor list buffer */
/*     nvlst       number of sites in list for each vdw site */
/*     vlst        site numbers in neighbor list of each vdw site */
/*     nelst       number of sites in list for each electrostatic site */
/*     elst        site numbers in list of each electrostatic site */
/*     dovlst      logical flag to rebuild vdw neighbor list */
/*     doclst      logical flag to rebuild charge neighbor list */
/*     domlst      logical flag to rebuild multipole neighbor list */




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




/*     zero out the van der Waals energy contribution */

    energi_1.ev = 0.;

/*     set arrays needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vscale[i__ - 1] = 1.;
	iv14[i__ - 1] = 0;
    }

/*     set cutoff distances and switching function coefficients */

    switch_("VDW", (ftnlen)3);
    expcut = -50.;

/*     apply any reduction factor to the atomic coordinates */

    i__1 = vdw_1.nvdw;
    for (k = 1; k <= i__1; ++k) {
	i__ = vdw_1.ivdw[k - 1];
	iv = vdw_1.ired[i__ - 1];
	rdn = vdw_1.kred[i__ - 1];
	xred[i__ - 1] = rdn * (atoms_1.x[i__ - 1] - atoms_1.x[iv - 1]) + 
		atoms_1.x[iv - 1];
	yred[i__ - 1] = rdn * (atoms_1.y[i__ - 1] - atoms_1.y[iv - 1]) + 
		atoms_1.y[iv - 1];
	zred[i__ - 1] = rdn * (atoms_1.z__[i__ - 1] - atoms_1.z__[iv - 1]) + 
		atoms_1.z__[iv - 1];
    }

/*     find the van der Waals energy via neighbor list search */

    i__1 = vdw_1.nvdw;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = vdw_1.ivdw[ii - 1];
	iv = vdw_1.ired[i__ - 1];
	it = vdw_1.jvdw[i__ - 1];
	xi = xred[i__ - 1];
	yi = yred[i__ - 1];
	zi = zred[i__ - 1];
	usei = usage_1.use[i__ - 1] || usage_1.use[iv - 1];

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

	i__2 = neigh_1.nvlst[ii - 1];
	for (kk = 1; kk <= i__2; ++kk) {
	    k = vdw_1.ivdw[vlst_ref(kk, ii) - 1];
	    kv = vdw_1.ired[k - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1] || usage_1.use[kv - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		kt = vdw_1.jvdw[k - 1];
		xr = xi - xred[k - 1];
		yr = yi - yred[k - 1];
		zr = zi - zred[k - 1];
		image_(&xr, &yr, &zr);
		rik2 = xr * xr + yr * yr + zr * zr;

/*     check for an interaction distance less than the cutoff */

		if (rik2 <= shunt_1.off2) {
/* Computing 2nd power */
		    d__1 = radmin_ref(kt, it);
		    rad2 = d__1 * d__1;
		    eps = epsilon_ref(kt, it);
		    if (iv14[k - 1] == i__) {
/* Computing 2nd power */
			d__1 = radmin4_ref(kt, it);
			rad2 = d__1 * d__1;
			eps = epsilon4_ref(kt, it);
		    }
		    eps *= vscale[k - 1];
		    i__3 = vdwpot_1.ngauss;
		    for (j = 1; j <= i__3; ++j) {
			a[j - 1] = igauss_ref(1, j) * eps;
			b[j - 1] = igauss_ref(2, j) / rad2;
		    }
		    e = 0.;
		    i__3 = vdwpot_1.ngauss;
		    for (j = 1; j <= i__3; ++j) {
			expterm = -b[j - 1] * rik2;
			if (expterm > expcut) {
			    e += a[j - 1] * exp(expterm);
			}
		    }

/*     use energy switching if near the cutoff distance */

		    if (rik2 > shunt_1.cut2) {
			rik = sqrt(rik2);
			rik3 = rik2 * rik;
			rik4 = rik2 * rik2;
			rik5 = rik2 * rik3;
			taper = shunt_1.c5 * rik5 + shunt_1.c4 * rik4 + 
				shunt_1.c3 * rik3 + shunt_1.c2 * rik2 + 
				shunt_1.c1 * rik + shunt_1.c0;
			e *= taper;
		    }

/*     scale the interaction based on its group membership */

		    if (group_1.use_group__) {
			e *= fgrp;
		    }

/*     increment the overall van der Waals energy components */

		    energi_1.ev += e;
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
} /* egauss0c_ */

#undef epsilon_ref
#undef radmin4_ref
#undef igauss_ref
#undef radmin_ref
#undef vlst_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref
#undef epsilon4_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine egauss0d  --  Gaussian vdw energy for smoothing  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "egauss0d" calculates the Gaussian expansion van der Waals */
/*     energy for use with potential energy smoothing */


/* Subroutine */ int egauss0d_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal expterm2, a[10], b[10], e;
    static integer i__, j, k;
    static doublereal t1, t2;
    static integer ii, kk, iv, it, kv, kt;
    static doublereal xi, yi, zi, xr, yr, zr;
    extern doublereal erf_(doublereal *);
    static integer iv14[25000];
    static doublereal rdn, rik, eps, rad2, rik2, fgrp, xred[25000], term, 
	    yred[25000], zred[25000];
    static logical usei;
    static doublereal width, broot, wterm, vscale[25000], expcut;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static logical proceed;
    static doublereal expterm;


#define epsilon4_ref(a_1,a_2) vdw_1.epsilon4[(a_2)*1000 + a_1 - 1001]
#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define radmin_ref(a_1,a_2) vdw_1.radmin[(a_2)*1000 + a_1 - 1001]
#define igauss_ref(a_1,a_2) vdwpot_1.igauss[(a_2)*2 + a_1 - 3]
#define radmin4_ref(a_1,a_2) vdw_1.radmin4[(a_2)*1000 + a_1 - 1001]
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




/*     zero out the van der Waals energy contribution */

    energi_1.ev = 0.;

/*     set arrays needed to scale connected atom interactions */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vscale[i__ - 1] = 1.;
	iv14[i__ - 1] = 0;
    }

/*     set the extent of smoothing to be performed */

    expcut = -50.;
    width = 0.;
    if (warp_1.use_dem__) {
	width = warp_1.diffv * 4. * warp_1.deform;
    } else if (warp_1.use_gda__) {
	wterm = warp_1.diffv * .66666666666666663;
    } else if (warp_1.use_tophat__) {
/* Computing MAX */
	d__1 = warp_1.diffv * warp_1.deform;
	width = max(d__1,1e-4);
    }

/*     apply any reduction factor to the atomic coordinates */

    i__1 = vdw_1.nvdw;
    for (k = 1; k <= i__1; ++k) {
	i__ = vdw_1.ivdw[k - 1];
	iv = vdw_1.ired[i__ - 1];
	rdn = vdw_1.kred[i__ - 1];
	xred[i__ - 1] = rdn * (atoms_1.x[i__ - 1] - atoms_1.x[iv - 1]) + 
		atoms_1.x[iv - 1];
	yred[i__ - 1] = rdn * (atoms_1.y[i__ - 1] - atoms_1.y[iv - 1]) + 
		atoms_1.y[iv - 1];
	zred[i__ - 1] = rdn * (atoms_1.z__[i__ - 1] - atoms_1.z__[iv - 1]) + 
		atoms_1.z__[iv - 1];
    }

/*     find the van der Waals energy via double loop search */

    i__1 = vdw_1.nvdw - 1;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = vdw_1.ivdw[ii - 1];
	iv = vdw_1.ired[i__ - 1];
	it = vdw_1.jvdw[i__ - 1];
	xi = xred[i__ - 1];
	yi = yred[i__ - 1];
	zi = zred[i__ - 1];
	usei = usage_1.use[i__ - 1] || usage_1.use[iv - 1];

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
	for (kk = ii + 1; kk <= i__2; ++kk) {
	    k = vdw_1.ivdw[kk - 1];
	    kv = vdw_1.ired[k - 1];
	    proceed = TRUE_;
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &i__, &k, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    if (proceed) {
		proceed = usei || usage_1.use[k - 1] || usage_1.use[kv - 1];
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		kt = vdw_1.jvdw[k - 1];
		xr = xi - xred[k - 1];
		yr = yi - yred[k - 1];
		zr = zi - zred[k - 1];
		rik2 = xr * xr + yr * yr + zr * zr;

/*     check for an interaction distance less than the cutoff */

/* Computing 2nd power */
		d__1 = radmin_ref(kt, it);
		rad2 = d__1 * d__1;
		eps = epsilon_ref(kt, it);
		if (iv14[k - 1] == i__) {
/* Computing 2nd power */
		    d__1 = radmin4_ref(kt, it);
		    rad2 = d__1 * d__1;
		    eps = epsilon4_ref(kt, it);
		}
		eps *= vscale[k - 1];
		i__3 = vdwpot_1.ngauss;
		for (j = 1; j <= i__3; ++j) {
		    a[j - 1] = igauss_ref(1, j) * eps;
		    b[j - 1] = igauss_ref(2, j) / rad2;
		}
		e = 0.;

/*     transform the potential function via smoothing */

		if (warp_1.use_tophat__) {
		    rik = sqrt(rik2);
		    i__3 = vdwpot_1.ngauss;
		    for (j = 1; j <= i__3; ++j) {
			broot = sqrt(b[j - 1]);
/* Computing 2nd power */
			d__1 = rik + width;
			expterm = -b[j - 1] * (d__1 * d__1);
			if (expterm > expcut) {
			    expterm = exp(expterm);
			} else {
			    expterm = 0.;
			}
/* Computing 2nd power */
			d__1 = width - rik;
			expterm2 = -b[j - 1] * (d__1 * d__1);
			if (expterm2 > expcut) {
			    expterm2 = exp(expterm2);
			} else {
			    expterm2 = 0.;
			}
			term = broot * (expterm - expterm2);
			d__1 = broot * (rik + width);
			d__2 = broot * (width - rik);
			term += b[j - 1] * 1.772453850905516027 * rik * (erf_(
				&d__1) + erf_(&d__2));
			e += term * a[j - 1] / (b[j - 1] * b[j - 1] * broot);
		    }
/* Computing 3rd power */
		    d__1 = width;
		    e = e * 3. / (rik * 8. * (d__1 * (d__1 * d__1)));
		} else {
		    if (warp_1.use_gda__) {
			width = wterm * (warp_1.m2[i__ - 1] + warp_1.m2[k - 1]
				);
		    }
		    i__3 = vdwpot_1.ngauss;
		    for (j = 1; j <= i__3; ++j) {
			t1 = b[j - 1] * width + 1.;
/* Computing 3rd power */
			d__1 = t1;
			t2 = sqrt(d__1 * (d__1 * d__1));
			expterm = -b[j - 1] * rik2 / t1;
			if (expterm > expcut) {
			    e += a[j - 1] / t2 * exp(expterm);
			}
		    }
		}

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    e *= fgrp;
		}

/*     increment the overall van der Waals energy components */

		energi_1.ev += e;
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
} /* egauss0d_ */

#undef epsilon_ref
#undef radmin4_ref
#undef igauss_ref
#undef radmin_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref
#undef epsilon4_ref


