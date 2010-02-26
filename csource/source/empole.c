/* empole.f -- translated by f2c (version 20050501).
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
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

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
    doublereal m2scale, m3scale, m4scale, m5scale;
} mplpot_;

#define mplpot_1 mplpot_

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
    doublereal lbuffer, lbuf2, vbuf2, cbuf2, mbuf2;
    integer nvlst[25000], vlst[45000000]	/* was [1800][25000] */, 
	    nelst[25000], elst[30000000]	/* was [1200][25000] */;
    logical dovlst, doclst, domlst;
} neigh_;

#define neigh_1 neigh_

struct {
    doublereal aewald;
    char boundary[7];
} ewald_;

#define ewald_1 ewald_

struct {
    doublereal bsmod1[100], bsmod2[100], bsmod3[100], table[1200]	/* 
	    was [400][3] */, qgrid[2000000]	/* was [2][100][100][100] */, 
	    qfac[1000000]	/* was [100][100][100] */, thetai1[1000000]	
	    /* was [4][10][25000] */, thetai2[1000000]	/* was [4][10][25000] 
	    */, thetai3[1000000]	/* was [4][10][25000] */;
    integer nfft1, nfft2, nfft3, bsorder, iprime[45]	/* was [15][3] */, 
	    igrid[75000]	/* was [3][25000] */;
} pme_;

#define pme_1 pme_

/* Table of constant values */

static integer c__0 = 0;



/*     ############################################################# */
/*     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine empole  --  multipole & polarization energy  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "empole" calculates the electrostatic energy due to */
/*     atomic multipole and dipole polarizability interactions */


/* Subroutine */ int empole_(void)
{
    extern /* Subroutine */ int empole0a_(void), empole0b_(void), empole0c_(
	    void), empole0d_(void);



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




/*     choose the method for summing over multipole interactions */



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


    if (cutoff_1.use_ewald__) {
	if (cutoff_1.use_mlist__) {
	    empole0d_();
	} else {
	    empole0c_();
	}
    } else {
	if (cutoff_1.use_mlist__) {
	    empole0b_();
	} else {
	    empole0a_();
	}
    }

/*     zero out potential energies which are not in use */

    if (! potent_1.use_mpole__) {
	energi_1.em = 0.;
    }
    if (! potent_1.use_polar__) {
	energi_1.ep = 0.;
    }
    return 0;
} /* empole_ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine empole0a  --  double loop multipole energy  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "empole0a" calculates the atomic multipole and dipole */
/*     polarizability interaction energy using a double loop */


/* Subroutine */ int empole0a_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__, r2;
    static integer ii, kk;
    static doublereal ei, fm;
    static integer ix;
    static doublereal fp;
    static integer iz, kx, kz;
    static doublereal ci, ck, sc[10], gl[5], xr, yr, zr, rr1, rr3, rr5, rr7, 
	    rr9, pdi, dix, diy, diz, pti, dkx, dky, dkz, qix, qiy, qiz, uix, 
	    uiy, uiz, ukx, uky, ukz, qkx, qky, qkz, sci[8], gli[3], damp, 
	    fgrp;
    static logical usei, usek;
    static doublereal qixx, qixy, qixz, qiyy, qiyz, qizz, qkxx, qkxy, qkxz, 
	    qkyy, qkyz, qkzz;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal scale3, scale5, scale7, pgamma, mscale[25000];
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal pscale[25000];
    extern /* Subroutine */ int induce_(void), switch_(char *, ftnlen), 
	    groups_(logical *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    static logical proceed;
    extern /* Subroutine */ int chkpole_(void);
    static doublereal expdamp;
    extern /* Subroutine */ int rotpole_(void);


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
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




/*     zero out the multipole and polarization energies */

    energi_1.em = 0.;
    energi_1.ep = 0.;

/*     check the sign of multipole components at chiral sites */

    chkpole_();

/*     rotate the multipole components into the global frame */

    rotpole_();

/*     compute the induced dipoles at each polarizable atom */

    induce_();

/*     set arrays needed to scale connected atom interactions */

    if (mpole_1.npole == 0) {
	return 0;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mscale[i__ - 1] = 1.;
	pscale[i__ - 1] = 1.;
    }

/*     set conversion factor, cutoff and switching coefficients */

    f = chgpot_1.electric / chgpot_1.dielec;
    switch_("MPOLE", (ftnlen)5);

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
	uix = uind_ref(1, i__);
	uiy = uind_ref(2, i__);
	uiz = uind_ref(3, i__);
	usei = usage_1.use[ii - 1] || usage_1.use[iz - 1] || usage_1.use[ix - 
		1];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = mplpot_1.m2scale;
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = mplpot_1.m3scale;
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = mplpot_1.m4scale;
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
	    mscale[i15_ref(j, ii) - 1] = mplpot_1.m5scale;
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
		    ukx = uind_ref(1, k);
		    uky = uind_ref(2, k);
		    ukz = uind_ref(3, k);

/*     construct some intermediate quadrupole values */

		    qix = qixx * xr + qixy * yr + qixz * zr;
		    qiy = qixy * xr + qiyy * yr + qiyz * zr;
		    qiz = qixz * xr + qiyz * yr + qizz * zr;
		    qkx = qkxx * xr + qkxy * yr + qkxz * zr;
		    qky = qkxy * xr + qkyy * yr + qkyz * zr;
		    qkz = qkxz * xr + qkyz * yr + qkzz * zr;

/*     calculate the scalar products for permanent multipoles */

		    sc[1] = dix * dkx + diy * dky + diz * dkz;
		    sc[2] = dix * xr + diy * yr + diz * zr;
		    sc[3] = dkx * xr + dky * yr + dkz * zr;
		    sc[4] = qix * xr + qiy * yr + qiz * zr;
		    sc[5] = qkx * xr + qky * yr + qkz * zr;
		    sc[6] = qix * dkx + qiy * dky + qiz * dkz;
		    sc[7] = qkx * dix + qky * diy + qkz * diz;
		    sc[8] = qix * qkx + qiy * qky + qiz * qkz;
		    sc[9] = (qixy * qkxy + qixz * qkxz + qiyz * qkyz) * 2. + 
			    qixx * qkxx + qiyy * qkyy + qizz * qkzz;

/*     calculate the scalar products for polarization components */

		    sci[1] = uix * dkx + dix * ukx + uiy * dky + diy * uky + 
			    uiz * dkz + diz * ukz;
		    sci[2] = uix * xr + uiy * yr + uiz * zr;
		    sci[3] = ukx * xr + uky * yr + ukz * zr;
		    sci[6] = qix * ukx + qiy * uky + qiz * ukz;
		    sci[7] = qkx * uix + qky * uiy + qkz * uiz;

/*     calculate the gl functions for permanent multipoles */

		    gl[0] = ci * ck;
		    gl[1] = ck * sc[2] - ci * sc[3] + sc[1];
		    gl[2] = ci * sc[5] + ck * sc[4] - sc[2] * sc[3] + (sc[6] 
			    - sc[7] + sc[9]) * 2.;
		    gl[3] = sc[2] * sc[5] - sc[3] * sc[4] - sc[8] * 4.;
		    gl[4] = sc[4] * sc[5];

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
		    rr9 = rr7 * 7. / r2;
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
			    expdamp = exp(damp);
			    scale3 = 1. - expdamp;
			    scale5 = 1. - (1. - damp) * expdamp;
/* Computing 2nd power */
			    d__1 = damp;
			    scale7 = 1. - (1. - damp + d__1 * d__1 * .6) * 
				    expdamp;
			}
		    }
		    e = gl[0] * rr1 + gl[1] * rr3 + gl[2] * rr5 + gl[3] * rr7 
			    + gl[4] * rr9;
		    ei = gli[0] * rr3 * scale3 + gli[1] * rr5 * scale5 + gli[
			    2] * rr7 * scale7;

/*     apply the energy adjustments for scaled interactions */

		    fm = f * mscale[kk - 1];
		    fp = f * pscale[kk - 1];
		    e = fm * e;
		    ei = fp * .5 * ei;

/*     scale the interaction based on its group membership; */
/*     polarization cannot be group scaled as it is not pairwise */

		    if (group_1.use_group__) {
			e *= fgrp;
/*                    ei = ei * fgrp */
		    }

/*     increment the overall multipole and polarization energies */

		    energi_1.em += e;
		    energi_1.ep += ei;
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = 1.;
	    pscale[i12_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = 1.;
	    pscale[i13_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = 1.;
	    pscale[i14_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i15_ref(j, ii) - 1] = 1.;
	    pscale[i15_ref(j, ii) - 1] = 1.;
	}
    }

/*     for periodic boundary conditions with large cutoffs */
/*     neighbors must be found by the replicates method */

    if (! bound_1.use_replica__) {
	return 0;
    }

/*     calculate interaction energy with other unit cells */

    i__1 = mpole_1.npole;
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
	uix = uind_ref(1, i__);
	uiy = uind_ref(2, i__);
	uiz = uind_ref(3, i__);
	usei = usage_1.use[ii - 1] || usage_1.use[iz - 1] || usage_1.use[ix - 
		1];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = mplpot_1.m2scale;
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = mplpot_1.m3scale;
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = mplpot_1.m4scale;
	    pscale[i14_ref(j, ii) - 1] = polpot_1.p4scale;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i15_ref(j, ii) - 1] = mplpot_1.m5scale;
	    pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	}

/*     decide whether to compute the current interaction */

	i__2 = mpole_1.npole;
	for (k = i__; k <= i__2; ++k) {
	    kk = mpole_1.ipole[k - 1];
	    kz = mpole_1.zaxis[k - 1];
	    kx = mpole_1.xaxis[k - 1];
	    usek = usage_1.use[kk - 1] || usage_1.use[kz - 1] || usage_1.use[
		    kx - 1];
	    if (group_1.use_group__) {
		groups_(&proceed, &fgrp, &ii, &kk, &c__0, &c__0, &c__0, &c__0)
			;
	    }
	    proceed = TRUE_;
	    if (proceed) {
		proceed = usei || usek;
	    }

/*     compute the energy contribution for this interaction */

	    if (proceed) {
		i__3 = cell_1.ncell;
		for (j = 1; j <= i__3; ++j) {
		    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
		    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
		    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
		    imager_(&xr, &yr, &zr, &j);
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
			ukx = uind_ref(1, k);
			uky = uind_ref(2, k);
			ukz = uind_ref(3, k);

/*     construct some intermediate quadrupole values */

			qix = qixx * xr + qixy * yr + qixz * zr;
			qiy = qixy * xr + qiyy * yr + qiyz * zr;
			qiz = qixz * xr + qiyz * yr + qizz * zr;
			qkx = qkxx * xr + qkxy * yr + qkxz * zr;
			qky = qkxy * xr + qkyy * yr + qkyz * zr;
			qkz = qkxz * xr + qkyz * yr + qkzz * zr;

/*     calculate the scalar products for permanent multipoles */

			sc[1] = dix * dkx + diy * dky + diz * dkz;
			sc[2] = dix * xr + diy * yr + diz * zr;
			sc[3] = dkx * xr + dky * yr + dkz * zr;
			sc[4] = qix * xr + qiy * yr + qiz * zr;
			sc[5] = qkx * xr + qky * yr + qkz * zr;
			sc[6] = qix * dkx + qiy * dky + qiz * dkz;
			sc[7] = qkx * dix + qky * diy + qkz * diz;
			sc[8] = qix * qkx + qiy * qky + qiz * qkz;
			sc[9] = (qixy * qkxy + qixz * qkxz + qiyz * qkyz) * 
				2. + qixx * qkxx + qiyy * qkyy + qizz * qkzz;

/*     calculate the scalar products for polarization components */

			sci[1] = uix * dkx + dix * ukx + uiy * dky + diy * 
				uky + uiz * dkz + diz * ukz;
			sci[2] = uix * xr + uiy * yr + uiz * zr;
			sci[3] = ukx * xr + uky * yr + ukz * zr;
			sci[6] = qix * ukx + qiy * uky + qiz * ukz;
			sci[7] = qkx * uix + qky * uiy + qkz * uiz;

/*     calculate the gl functions for permanent multipoles */

			gl[0] = ci * ck;
			gl[1] = ck * sc[2] - ci * sc[3] + sc[1];
			gl[2] = ci * sc[5] + ck * sc[4] - sc[2] * sc[3] + (sc[
				6] - sc[7] + sc[9]) * 2.;
			gl[3] = sc[2] * sc[5] - sc[3] * sc[4] - sc[8] * 4.;
			gl[4] = sc[4] * sc[5];

/*     calculate the gl functions for polarization components */

			gli[0] = ck * sci[2] - ci * sci[3] + sci[1];
			gli[1] = (sci[6] - sci[7]) * 2. - sci[2] * sc[3] - sc[
				2] * sci[3];
			gli[2] = sci[2] * sc[5] - sci[3] * sc[4];

/*     compute the energy contributions for this interaction */

			rr1 = 1. / r__;
			rr3 = rr1 / r2;
			rr5 = rr3 * 3. / r2;
			rr7 = rr5 * 5. / r2;
			rr9 = rr7 * 7. / r2;
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
				expdamp = exp(damp);
				scale3 = 1. - expdamp;
				scale5 = 1. - (1. - damp) * expdamp;
/* Computing 2nd power */
				d__1 = damp;
				scale7 = 1. - (1. - damp + d__1 * d__1 * .6) *
					 expdamp;
			    }
			}
			e = gl[0] * rr1 + gl[1] * rr3 + gl[2] * rr5 + gl[3] * 
				rr7 + gl[4] * rr9;
			ei = gli[0] * rr3 * scale3 + gli[1] * rr5 * scale5 + 
				gli[2] * rr7 * scale7;

/*     apply the energy adjustments for scaled interactions */

			fm = f;
			fp = f;
			if (bound_1.use_polymer__) {
			    if (r2 <= bound_1.polycut2) {
				fm *= mscale[kk - 1];
				fp *= pscale[kk - 1];
			    }
			}
			e = fm * e;
			ei = fp * .5 * ei;

/*     scale the interaction based on its group membership; */
/*     polarization cannot be group scaled as it is not pairwise */

			if (group_1.use_group__) {
			    e *= fgrp;
/*                       ei = ei * fgrp */
			}

/*     increment the overall multipole and polarization energies */

			if (ii == kk) {
			    e *= .5;
			    ei *= .5;
			}
			energi_1.em += e;
			energi_1.ep += ei;
		    }
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = 1.;
	    pscale[i12_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = 1.;
	    pscale[i13_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = 1.;
	    pscale[i14_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i15_ref(j, ii) - 1] = 1.;
	    pscale[i15_ref(j, ii) - 1] = 1.;
	}
    }
    return 0;
} /* empole0a_ */

#undef rpole_ref
#undef uind_ref
#undef ip11_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine empole0b  --  neighbor list multipole energy  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "empole0b" calculates the atomic multipole and dipole */
/*     polarizability interaction energy using a neighbor list */


/* Subroutine */ int empole0b_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k;
    static doublereal r__, r2;
    static integer ii, kk;
    static doublereal ei, fm;
    static integer ix;
    static doublereal fp;
    static integer iz, kx, kz;
    static doublereal ci, ck, sc[10], gl[5], xr, yr, zr, rr1, rr3, rr5, rr7, 
	    rr9;
    static integer kkk;
    static doublereal pdi, dix, diy, pti, diz, dkx, dky, dkz, qix, qiy, uix, 
	    uiy, uiz, ukx, uky, ukz, qiz, qkx, qky, qkz, sci[8], gli[3], damp,
	     fgrp;
    static logical usei, usek;
    static doublereal qixx, qixy, qixz, qiyy, qiyz, qizz, qkxx, qkxy, qkxz, 
	    qkyy, qkyz, qkzz;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal scale3, scale5, scale7, pgamma, mscale[25000], pscale[
	    25000];
    extern /* Subroutine */ int induce_(void), switch_(char *, ftnlen), 
	    groups_(logical *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    static logical proceed;
    extern /* Subroutine */ int chkpole_(void);
    static doublereal expdamp;
    extern /* Subroutine */ int rotpole_(void);


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define elst_ref(a_1,a_2) neigh_1.elst[(a_2)*1200 + a_1 - 1201]
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




/*     zero out the multipole and polarization energies */

    energi_1.em = 0.;
    energi_1.ep = 0.;

/*     check the sign of multipole components at chiral sites */

    chkpole_();

/*     rotate the multipole components into the global frame */

    rotpole_();

/*     compute the induced dipoles at each polarizable atom */

    induce_();

/*     set arrays needed to scale connected atom interactions */

    if (mpole_1.npole == 0) {
	return 0;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mscale[i__ - 1] = 1.;
	pscale[i__ - 1] = 1.;
    }

/*     set conversion factor, cutoff and switching coefficients */

    f = chgpot_1.electric / chgpot_1.dielec;
    switch_("MPOLE", (ftnlen)5);

/*     calculate the multipole interaction energy term */

    i__1 = mpole_1.npole;
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
	uix = uind_ref(1, i__);
	uiy = uind_ref(2, i__);
	uiz = uind_ref(3, i__);
	usei = usage_1.use[ii - 1] || usage_1.use[iz - 1] || usage_1.use[ix - 
		1];

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = mplpot_1.m2scale;
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = mplpot_1.m3scale;
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = mplpot_1.m4scale;
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
	    mscale[i15_ref(j, ii) - 1] = mplpot_1.m5scale;
	    pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	}

/*     decide whether to compute the current interaction */

	i__2 = neigh_1.nelst[i__ - 1];
	for (kkk = 1; kkk <= i__2; ++kkk) {
	    k = elst_ref(kkk, i__);
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
		    ukx = uind_ref(1, k);
		    uky = uind_ref(2, k);
		    ukz = uind_ref(3, k);

/*     construct some intermediate quadrupole values */

		    qix = qixx * xr + qixy * yr + qixz * zr;
		    qiy = qixy * xr + qiyy * yr + qiyz * zr;
		    qiz = qixz * xr + qiyz * yr + qizz * zr;
		    qkx = qkxx * xr + qkxy * yr + qkxz * zr;
		    qky = qkxy * xr + qkyy * yr + qkyz * zr;
		    qkz = qkxz * xr + qkyz * yr + qkzz * zr;

/*     calculate the scalar products for permanent multipoles */

		    sc[1] = dix * dkx + diy * dky + diz * dkz;
		    sc[2] = dix * xr + diy * yr + diz * zr;
		    sc[3] = dkx * xr + dky * yr + dkz * zr;
		    sc[4] = qix * xr + qiy * yr + qiz * zr;
		    sc[5] = qkx * xr + qky * yr + qkz * zr;
		    sc[6] = qix * dkx + qiy * dky + qiz * dkz;
		    sc[7] = qkx * dix + qky * diy + qkz * diz;
		    sc[8] = qix * qkx + qiy * qky + qiz * qkz;
		    sc[9] = (qixy * qkxy + qixz * qkxz + qiyz * qkyz) * 2. + 
			    qixx * qkxx + qiyy * qkyy + qizz * qkzz;

/*     calculate the scalar products for polarization components */

		    sci[1] = uix * dkx + dix * ukx + uiy * dky + diy * uky + 
			    uiz * dkz + diz * ukz;
		    sci[2] = uix * xr + uiy * yr + uiz * zr;
		    sci[3] = ukx * xr + uky * yr + ukz * zr;
		    sci[6] = qix * ukx + qiy * uky + qiz * ukz;
		    sci[7] = qkx * uix + qky * uiy + qkz * uiz;

/*     calculate the gl functions for permanent multipoles */

		    gl[0] = ci * ck;
		    gl[1] = ck * sc[2] - ci * sc[3] + sc[1];
		    gl[2] = ci * sc[5] + ck * sc[4] - sc[2] * sc[3] + (sc[6] 
			    - sc[7] + sc[9]) * 2.;
		    gl[3] = sc[2] * sc[5] - sc[3] * sc[4] - sc[8] * 4.;
		    gl[4] = sc[4] * sc[5];

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
		    rr9 = rr7 * 7. / r2;
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
			    expdamp = exp(damp);
			    scale3 = 1. - expdamp;
			    scale5 = 1. - (1. - damp) * expdamp;
/* Computing 2nd power */
			    d__1 = damp;
			    scale7 = 1. - (1. - damp + d__1 * d__1 * .6) * 
				    expdamp;
			}
		    }
		    e = gl[0] * rr1 + gl[1] * rr3 + gl[2] * rr5 + gl[3] * rr7 
			    + gl[4] * rr9;
		    ei = gli[0] * rr3 * scale3 + gli[1] * rr5 * scale5 + gli[
			    2] * rr7 * scale7;

/*     apply the energy adjustments for scaled interactions */

		    fm = f * mscale[kk - 1];
		    fp = f * pscale[kk - 1];
		    e = fm * e;
		    ei = fp * .5 * ei;

/*     scale the interaction based on its group membership; */
/*     polarization cannot be group scaled as it is not pairwise */

		    if (group_1.use_group__) {
			e *= fgrp;
/*                    ei = ei * fgrp */
		    }

/*     increment the overall multipole and polarization energies */

		    energi_1.em += e;
		    energi_1.ep += ei;
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = 1.;
	    pscale[i12_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = 1.;
	    pscale[i13_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = 1.;
	    pscale[i14_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i15_ref(j, ii) - 1] = 1.;
	    pscale[i15_ref(j, ii) - 1] = 1.;
	}
    }
    return 0;
} /* empole0b_ */

#undef rpole_ref
#undef elst_ref
#undef uind_ref
#undef ip11_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine empole0c  --  Ewald multipole energy via loop  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "empole0c" calculates the atomic multipole and dipole */
/*     polarizability interactions energy using a particle mesh */
/*     Ewald summation and double loop */


/* Subroutine */ int empole0c_(void)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal e, f;
    static integer i__;
    static doublereal ci, ei;
    static integer ii;
    static doublereal xd, yd, zd, xu, yu, zu, cii, dii, qii, dix, diy, uii, 
	    diz, uix, uiy, uiz, term, qixx, qixy, qixz, qiyy, qiyz, qizz, 
	    fterm;
    extern /* Subroutine */ int induce_(void), ereal0c_(void), emrecip_(void),
	     chkpole_(void), rotpole_(void);


#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
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
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  ewald.i  --  parameters and options for Ewald summation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     aewald     Ewald convergence coefficient value (Ang-1) */
/*     boundary   Ewald boundary condition; none, tinfoil or vacuum */




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




/*     zero out the multipole and polarization energies */

    energi_1.em = 0.;
    energi_1.ep = 0.;

/*     set the energy unit conversion factor */

    f = chgpot_1.electric / chgpot_1.dielec;

/*     check the sign of multipole components at chiral sites */

    chkpole_();

/*     rotate the multipole components into the global frame */

    rotpole_();

/*     compute the induced dipoles at each polarizable atom */

    induce_();

/*     compute the self-energy portion of the Ewald summation */

    term = ewald_1.aewald * 2. * ewald_1.aewald;
    fterm = -f * ewald_1.aewald / 1.772453850905516027;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
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
	uix = uind_ref(1, i__);
	uiy = uind_ref(2, i__);
	uiz = uind_ref(3, i__);
	cii = ci * ci;
	dii = dix * dix + diy * diy + diz * diz;
	qii = qixx * qixx + qiyy * qiyy + qizz * qizz + (qixy * qixy + qixz * 
		qixz + qiyz * qiyz) * 2.;
	uii = dix * uix + diy * uiy + diz * uiz;
	e = fterm * (cii + term * (dii / 3. + term * 2. * qii / 5.));
	ei = fterm * term * uii / 3.;
	energi_1.em += e;
	energi_1.ep += ei;
    }

/*     compute the cell dipole boundary correction term */

    if (s_cmp(ewald_1.boundary, "VACUUM", (ftnlen)7, (ftnlen)6) == 0) {
	xd = 0.;
	yd = 0.;
	zd = 0.;
	xu = 0.;
	yu = 0.;
	zu = 0.;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    dix = rpole_ref(2, i__);
	    diy = rpole_ref(3, i__);
	    diz = rpole_ref(4, i__);
	    uix = uind_ref(1, i__);
	    uiy = uind_ref(2, i__);
	    uiz = uind_ref(3, i__);
	    xd = xd + dix + rpole_ref(1, i__) * atoms_1.x[ii - 1];
	    yd = yd + diy + rpole_ref(1, i__) * atoms_1.y[ii - 1];
	    zd = zd + diz + rpole_ref(1, i__) * atoms_1.z__[ii - 1];
	    xu += uix;
	    yu += uiy;
	    zu += uiz;
	}
	term = f * .66666666666666663 * (3.141592653589793238 / 
		boxes_1.volbox);
	energi_1.em += term * (xd * xd + yd * yd + zd * zd);
	energi_1.ep += term * (xd * xu + yd * yu + zd * zu);
    }

/*     compute the reciprocal space part of the Ewald summation */

    emrecip_();

/*     compute the real space portion of the Ewald summation */

    ereal0c_();
    return 0;
} /* empole0c_ */

#undef rpole_ref
#undef uind_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine ereal0c  --  real space mpole energy via loop  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "ereal0c" evaluates the real space portion of the Ewald sum */
/*     energy due to atomic multipoles and dipole polarizability */
/*     using a double loop */

/*     literature reference: */

/*     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)", */
/*     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/) */


/* Subroutine */ int ereal0c_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k, m;
    static doublereal r__, r2;
    static integer ii, kk;
    static doublereal ei, ci, ck, sc[10], gl[5], bn[5], xr, yr, zr, rr1, rr3, 
	    rr5, rr7, rr9, pdi, dix, diy, diz, pti, dkx, dky, dkz, qix, qiy, 
	    qiz, qkx, uix, uiy, uiz, ukx, uky, ukz, qky, qkz, sci[8], gli[3], 
	    bfac;
    extern doublereal erfc_(doublereal *);
    static doublereal damp, efix, qixx, qixy, qixz, qiyy, qiyz, qizz, qkxx, 
	    qkxy, qkxz, qkyy, qkyz, qkzz, exp2a, alsq2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal eifix, scale3, scale5, scale7, alsq2n, pgamma, mscale[
	    25000];
    extern /* Subroutine */ int imager_(doublereal *, doublereal *, 
	    doublereal *, integer *);
    static doublereal ralpha, pscale[25000];
    extern /* Subroutine */ int switch_(char *, ftnlen);
    static doublereal expdamp;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
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
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  ewald.i  --  parameters and options for Ewald summation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     aewald     Ewald convergence coefficient value (Ang-1) */
/*     boundary   Ewald boundary condition; none, tinfoil or vacuum */




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




/*     set arrays needed to scale connected atom interactions */

    if (mpole_1.npole == 0) {
	return 0;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mscale[i__ - 1] = 1.;
	pscale[i__ - 1] = 1.;
    }

/*     set conversion factor, cutoff and switching coefficients */

    f = chgpot_1.electric / chgpot_1.dielec;
    switch_("EWALD", (ftnlen)5);

/*     compute the real space portion of the Ewald summation */

    i__1 = mpole_1.npole - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
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
	uix = uind_ref(1, i__);
	uiy = uind_ref(2, i__);
	uiz = uind_ref(3, i__);

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = mplpot_1.m2scale;
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = mplpot_1.m3scale;
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = mplpot_1.m4scale;
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
	    mscale[i15_ref(j, ii) - 1] = mplpot_1.m5scale;
	    pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	}
	i__2 = mpole_1.npole;
	for (k = i__ + 1; k <= i__2; ++k) {
	    kk = mpole_1.ipole[k - 1];
	    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
	    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
	    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
	    image_(&xr, &yr, &zr);
	    r2 = xr * xr + yr * yr + zr * zr;
	    if (r2 <= shunt_1.cut2) {
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
		ukx = uind_ref(1, k);
		uky = uind_ref(2, k);
		ukz = uind_ref(3, k);

/*     calculate the error function damping terms */

		ralpha = ewald_1.aewald * r__;
		bn[0] = erfc_(&ralpha) / r__;
/* Computing 2nd power */
		d__1 = ewald_1.aewald;
		alsq2 = d__1 * d__1 * 2.;
		alsq2n = 0.;
		if (ewald_1.aewald > 0.) {
		    alsq2n = 1. / (ewald_1.aewald * 1.772453850905516027);
		}
/* Computing 2nd power */
		d__1 = ralpha;
		exp2a = exp(-(d__1 * d__1));
		for (m = 1; m <= 4; ++m) {
		    bfac = (doublereal) (m + m - 1);
		    alsq2n = alsq2 * alsq2n;
		    bn[m] = (bfac * bn[m - 1] + alsq2n * exp2a) / r2;
		}

/*     construct some intermediate quadrupole values */

		qix = qixx * xr + qixy * yr + qixz * zr;
		qiy = qixy * xr + qiyy * yr + qiyz * zr;
		qiz = qixz * xr + qiyz * yr + qizz * zr;
		qkx = qkxx * xr + qkxy * yr + qkxz * zr;
		qky = qkxy * xr + qkyy * yr + qkyz * zr;
		qkz = qkxz * xr + qkyz * yr + qkzz * zr;

/*     calculate the scalar products for permanent multipoles */

		sc[1] = dix * dkx + diy * dky + diz * dkz;
		sc[2] = dix * xr + diy * yr + diz * zr;
		sc[3] = dkx * xr + dky * yr + dkz * zr;
		sc[4] = qix * xr + qiy * yr + qiz * zr;
		sc[5] = qkx * xr + qky * yr + qkz * zr;
		sc[6] = qix * dkx + qiy * dky + qiz * dkz;
		sc[7] = qkx * dix + qky * diy + qkz * diz;
		sc[8] = qix * qkx + qiy * qky + qiz * qkz;
		sc[9] = (qixy * qkxy + qixz * qkxz + qiyz * qkyz) * 2. + qixx 
			* qkxx + qiyy * qkyy + qizz * qkzz;

/*     calculate the scalar products for polarization components */

		sci[1] = uix * dkx + dix * ukx + uiy * dky + diy * uky + uiz *
			 dkz + diz * ukz;
		sci[2] = uix * xr + uiy * yr + uiz * zr;
		sci[3] = ukx * xr + uky * yr + ukz * zr;
		sci[6] = qix * ukx + qiy * uky + qiz * ukz;
		sci[7] = qkx * uix + qky * uiy + qkz * uiz;

/*     calculate the gl functions for permanent multipoles */

		gl[0] = ci * ck;
		gl[1] = ck * sc[2] - ci * sc[3] + sc[1];
		gl[2] = ci * sc[5] + ck * sc[4] - sc[2] * sc[3] + (sc[6] - sc[
			7] + sc[9]) * 2.;
		gl[3] = sc[2] * sc[5] - sc[3] * sc[4] - sc[8] * 4.;
		gl[4] = sc[4] * sc[5];

/*     calculate the gl functions for polarization components */

		gli[0] = ck * sci[2] - ci * sci[3] + sci[1];
		gli[1] = (sci[6] - sci[7]) * 2. - sci[2] * sc[3] - sc[2] * 
			sci[3];
		gli[2] = sci[2] * sc[5] - sci[3] * sc[4];

/*     compute the energy contributions for this interaction */

		e = gl[0] * bn[0] + gl[1] * bn[1] + gl[2] * bn[2] + gl[3] * 
			bn[3] + gl[4] * bn[4];
		ei = gli[0] * bn[1] + gli[1] * bn[2] + gli[2] * bn[3];

/*     full real space energies needed for scaled interactions */

		rr1 = 1. / r__;
		rr3 = rr1 / r2;
		rr5 = rr3 * 3. / r2;
		rr7 = rr5 * 5. / r2;
		rr9 = rr7 * 7. / r2;
		scale3 = pscale[kk - 1];
		scale5 = pscale[kk - 1];
		scale7 = pscale[kk - 1];
		damp = pdi * polar_1.pdamp[k - 1];
		if (damp != 0.) {
/* Computing MIN */
		    d__1 = pti, d__2 = polar_1.thole[k - 1];
		    pgamma = min(d__1,d__2);
/* Computing 3rd power */
		    d__1 = r__ / damp;
		    damp = -pgamma * (d__1 * (d__1 * d__1));
		    if (damp > -50.) {
			expdamp = exp(damp);
			scale3 *= 1. - expdamp;
			scale5 *= 1. - (1. - damp) * expdamp;
/* Computing 2nd power */
			d__1 = damp;
			scale7 *= 1. - (1. - damp + d__1 * d__1 * .6) * 
				expdamp;
		    }
		}
		efix = gl[0] * rr1 + gl[1] * rr3 + gl[2] * rr5 + gl[3] * rr7 
			+ gl[4] * rr9;
		eifix = gli[0] * rr3 * (1. - scale3) + gli[1] * rr5 * (1. - 
			scale5) + gli[2] * rr7 * (1. - scale7);

/*     apply the energy adjustments for scaled interactions */

		e -= efix * (1. - mscale[kk - 1]);
		ei -= eifix;

/*     increment the overall multipole and polarization energies */

		e = f * e;
		ei = f * .5 * ei;
		energi_1.em += e;
		energi_1.ep += ei;
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = 1.;
	    pscale[i12_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = 1.;
	    pscale[i13_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = 1.;
	    pscale[i14_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i15_ref(j, ii) - 1] = 1.;
	    pscale[i15_ref(j, ii) - 1] = 1.;
	}
    }

/*     for periodic boundary conditions with large cutoffs */
/*     neighbors must be found by the replicates method */

    if (! bound_1.use_replica__) {
	return 0;
    }

/*     calculate interaction energy with other unit cells */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
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
	uix = uind_ref(1, i__);
	uiy = uind_ref(2, i__);
	uiz = uind_ref(3, i__);

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = mplpot_1.m2scale;
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = mplpot_1.m3scale;
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = mplpot_1.m4scale;
	    pscale[i14_ref(j, ii) - 1] = polpot_1.p4scale;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i15_ref(j, ii) - 1] = mplpot_1.m5scale;
	    pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	}
	i__2 = mpole_1.npole;
	for (k = i__; k <= i__2; ++k) {
	    kk = mpole_1.ipole[k - 1];
	    i__3 = cell_1.ncell;
	    for (j = 1; j <= i__3; ++j) {
		xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
		yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
		zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
		imager_(&xr, &yr, &zr, &j);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.cut2) {
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
		    ukx = uind_ref(1, k);
		    uky = uind_ref(2, k);
		    ukz = uind_ref(3, k);

/*     calculate the error function damping terms */

		    ralpha = ewald_1.aewald * r__;
		    bn[0] = erfc_(&ralpha) / r__;
/* Computing 2nd power */
		    d__1 = ewald_1.aewald;
		    alsq2 = d__1 * d__1 * 2.;
		    alsq2n = 0.;
		    if (ewald_1.aewald > 0.) {
			alsq2n = 1. / (ewald_1.aewald * 1.772453850905516027);
		    }
/* Computing 2nd power */
		    d__1 = ralpha;
		    exp2a = exp(-(d__1 * d__1));
		    for (m = 1; m <= 4; ++m) {
			bfac = (doublereal) (m + m - 1);
			alsq2n = alsq2 * alsq2n;
			bn[m] = (bfac * bn[m - 1] + alsq2n * exp2a) / r2;
		    }

/*     construct some intermediate quadrupole values */

		    qix = qixx * xr + qixy * yr + qixz * zr;
		    qiy = qixy * xr + qiyy * yr + qiyz * zr;
		    qiz = qixz * xr + qiyz * yr + qizz * zr;
		    qkx = qkxx * xr + qkxy * yr + qkxz * zr;
		    qky = qkxy * xr + qkyy * yr + qkyz * zr;
		    qkz = qkxz * xr + qkyz * yr + qkzz * zr;

/*     calculate the scalar products for permanent multipoles */

		    sc[1] = dix * dkx + diy * dky + diz * dkz;
		    sc[2] = dix * xr + diy * yr + diz * zr;
		    sc[3] = dkx * xr + dky * yr + dkz * zr;
		    sc[4] = qix * xr + qiy * yr + qiz * zr;
		    sc[5] = qkx * xr + qky * yr + qkz * zr;
		    sc[6] = qix * dkx + qiy * dky + qiz * dkz;
		    sc[7] = qkx * dix + qky * diy + qkz * diz;
		    sc[8] = qix * qkx + qiy * qky + qiz * qkz;
		    sc[9] = (qixy * qkxy + qixz * qkxz + qiyz * qkyz) * 2. + 
			    qixx * qkxx + qiyy * qkyy + qizz * qkzz;

/*     calculate the scalar products for polarization components */

		    sci[1] = uix * dkx + dix * ukx + uiy * dky + diy * uky + 
			    uiz * dkz + diz * ukz;
		    sci[2] = uix * xr + uiy * yr + uiz * zr;
		    sci[3] = ukx * xr + uky * yr + ukz * zr;
		    sci[6] = qix * ukx + qiy * uky + qiz * ukz;
		    sci[7] = qkx * uix + qky * uiy + qkz * uiz;

/*     calculate the gl functions for permanent multipoles */

		    gl[0] = ci * ck;
		    gl[1] = ck * sc[2] - ci * sc[3] + sc[1];
		    gl[2] = ci * sc[5] + ck * sc[4] - sc[2] * sc[3] + (sc[6] 
			    - sc[7] + sc[9]) * 2.;
		    gl[3] = sc[2] * sc[5] - sc[3] * sc[4] - sc[8] * 4.;
		    gl[4] = sc[4] * sc[5];

/*     calculate the gl functions for polarization components */

		    gli[0] = ck * sci[2] - ci * sci[3] + sci[1];
		    gli[1] = (sci[6] - sci[7]) * 2. - sci[2] * sc[3] - sc[2] *
			     sci[3];
		    gli[2] = sci[2] * sc[5] - sci[3] * sc[4];

/*     compute the energy contributions for this interaction */

		    e = gl[0] * bn[0] + gl[1] * bn[1] + gl[2] * bn[2] + gl[3] 
			    * bn[3] + gl[4] * bn[4];
		    ei = gli[0] * bn[1] + gli[1] * bn[2] + gli[2] * bn[3];

/*     full real space energies needeed for scaled interactions */

		    rr1 = 1. / r__;
		    rr3 = rr1 / r2;
		    rr5 = rr3 * 3. / r2;
		    rr7 = rr5 * 5. / r2;
		    rr9 = rr7 * 7. / r2;
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
			    expdamp = exp(damp);
			    scale3 = 1. - expdamp;
			    scale5 = 1. - (1. - damp) * expdamp;
/* Computing 2nd power */
			    d__1 = damp;
			    scale7 = 1. - (1. - damp + d__1 * d__1 * .6) * 
				    expdamp;
			    if (bound_1.use_polymer__ && r2 <= 
				    bound_1.polycut2) {
				scale3 *= pscale[kk - 1];
				scale5 *= pscale[kk - 1];
				scale7 *= pscale[kk - 1];
			    }
			}
		    }
		    efix = gl[0] * rr1 + gl[1] * rr3 + gl[2] * rr5 + gl[3] * 
			    rr7 + gl[4] * rr9;
		    eifix = gli[0] * rr3 * (1. - scale3) + gli[1] * rr5 * (1. 
			    - scale5) + gli[2] * rr7 * (1. - scale7);

/*     apply the energy adjustments for scaled interactions */

		    if (bound_1.use_polymer__ && r2 <= bound_1.polycut2) {
			e -= efix * (1. - mscale[kk - 1]);
		    }
		    ei -= eifix;

/*     increment the overall multipole and polarization energies */

		    e = f * e;
		    ei = f * .5 * ei;
		    if (ii == kk) {
			e *= .5;
			ei *= .5;
		    }
		    energi_1.em += e;
		    energi_1.ep += ei;
		}
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = 1.;
	    pscale[i12_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = 1.;
	    pscale[i13_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = 1.;
	    pscale[i14_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i15_ref(j, ii) - 1] = 1.;
	    pscale[i15_ref(j, ii) - 1] = 1.;
	}
    }
    return 0;
} /* ereal0c_ */

#undef rpole_ref
#undef uind_ref
#undef ip11_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine empole0d  --  Ewald multipole energy via list  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "empole0d" calculates the atomic multipole and dipole */
/*     polarizability interaction energy using a particle mesh */
/*     Ewald summation and a neighbor list */


/* Subroutine */ int empole0d_(void)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal e, f;
    static integer i__;
    static doublereal ci, ei;
    static integer ii;
    static doublereal xd, yd, zd, xu, yu, zu, cii, dii, qii, dix, diy, uii, 
	    diz, uix, uiy, uiz, term, qixx, qixy, qixz, qiyy, qiyz, qizz, 
	    fterm;
    extern /* Subroutine */ int induce_(void), ereal0d_(void), emrecip_(void),
	     chkpole_(void), rotpole_(void);


#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
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
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  ewald.i  --  parameters and options for Ewald summation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     aewald     Ewald convergence coefficient value (Ang-1) */
/*     boundary   Ewald boundary condition; none, tinfoil or vacuum */




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




/*     zero out the multipole and polarization energies */

    energi_1.em = 0.;
    energi_1.ep = 0.;

/*     set the energy unit conversion factor */

    f = chgpot_1.electric / chgpot_1.dielec;

/*     check the sign of multipole components at chiral sites */

    chkpole_();

/*     rotate the multipole components into the global frame */

    rotpole_();

/*     compute the induced dipoles at each polarizable atom */

    induce_();

/*     compute the self-energy portion of the Ewald summation */

    term = ewald_1.aewald * 2. * ewald_1.aewald;
    fterm = -f * ewald_1.aewald / 1.772453850905516027;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
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
	uix = uind_ref(1, i__);
	uiy = uind_ref(2, i__);
	uiz = uind_ref(3, i__);
	cii = ci * ci;
	dii = dix * dix + diy * diy + diz * diz;
	qii = qixx * qixx + qiyy * qiyy + qizz * qizz + (qixy * qixy + qixz * 
		qixz + qiyz * qiyz) * 2.;
	uii = dix * uix + diy * uiy + diz * uiz;
	e = fterm * (cii + term * (dii / 3. + term * 2. * qii / 5.));
	ei = fterm * term * uii / 3.;
	energi_1.em += e;
	energi_1.ep += ei;
    }

/*     compute the cell dipole boundary correction term */

    if (s_cmp(ewald_1.boundary, "VACUUM", (ftnlen)7, (ftnlen)6) == 0) {
	xd = 0.;
	yd = 0.;
	zd = 0.;
	xu = 0.;
	yu = 0.;
	zu = 0.;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = mpole_1.ipole[i__ - 1];
	    dix = rpole_ref(2, i__);
	    diy = rpole_ref(3, i__);
	    diz = rpole_ref(4, i__);
	    uix = uind_ref(1, i__);
	    uiy = uind_ref(2, i__);
	    uiz = uind_ref(3, i__);
	    xd = xd + dix + rpole_ref(1, i__) * atoms_1.x[ii - 1];
	    yd = yd + diy + rpole_ref(1, i__) * atoms_1.y[ii - 1];
	    zd = zd + diz + rpole_ref(1, i__) * atoms_1.z__[ii - 1];
	    xu += uix;
	    yu += uiy;
	    zu += uiz;
	}
	term = f * .66666666666666663 * (3.141592653589793238 / 
		boxes_1.volbox);
	energi_1.em += term * (xd * xd + yd * yd + zd * zd);
	energi_1.ep += term * (xd * xu + yd * yu + zd * zu);
    }

/*     compute the reciprocal space part of the Ewald summation */

    emrecip_();

/*     compute the real space portion of the Ewald summation */

    ereal0d_();
    return 0;
} /* empole0d_ */

#undef rpole_ref
#undef uind_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine ereal0d  --  real space mpole energy via list  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "ereal0d" evaluates the real space portion of the Ewald sum */
/*     energy due to atomic multipoles and dipole polarizability */
/*     using a neighbor list */

/*     literature reference: */

/*     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)", */
/*     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/) */


/* Subroutine */ int ereal0d_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal e, f;
    static integer i__, j, k, m;
    static doublereal r__, r2;
    static integer ii, kk;
    static doublereal ei, ci, ck, sc[10], gl[5], bn[5], xr, yr, zr, rr1, rr3, 
	    rr5, rr7, rr9;
    static integer kkk;
    static doublereal pdi, dix, diy, pti, diz, dkx, dky, dkz, qix, qiy, qiz, 
	    uix, uiy, uiz, ukx, uky, ukz, qkx, qky, qkz, sci[8], gli[3], bfac;
    extern doublereal erfc_(doublereal *);
    static doublereal damp, efix, qixx, qixy, qixz, qiyy, qiyz, qizz, qkxx, 
	    qkxy, qkxz, qkyy, qkyz, qkzz, exp2a, alsq2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal eifix, scale3, scale5, scale7, alsq2n, pgamma, mscale[
	    25000], ralpha, pscale[25000];
    extern /* Subroutine */ int switch_(char *, ftnlen);
    static doublereal expdamp;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define i15_ref(a_1,a_2) couple_1.i15[(a_2)*216 + a_1 - 217]
#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define elst_ref(a_1,a_2) neigh_1.elst[(a_2)*1200 + a_1 - 1201]
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
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  ewald.i  --  parameters and options for Ewald summation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     aewald     Ewald convergence coefficient value (Ang-1) */
/*     boundary   Ewald boundary condition; none, tinfoil or vacuum */




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




/*     set arrays needed to scale connected atom interactions */

    if (mpole_1.npole == 0) {
	return 0;
    }
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mscale[i__ - 1] = 1.;
	pscale[i__ - 1] = 1.;
    }

/*     set conversion factor, cutoff and switching coefficients */

    f = chgpot_1.electric / chgpot_1.dielec;
    switch_("EWALD", (ftnlen)5);

/*     compute the real space portion of the Ewald summation */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = mpole_1.ipole[i__ - 1];
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
	uix = uind_ref(1, i__);
	uiy = uind_ref(2, i__);
	uiz = uind_ref(3, i__);

/*     set interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = mplpot_1.m2scale;
	    pscale[i12_ref(j, ii) - 1] = polpot_1.p2scale;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = mplpot_1.m3scale;
	    pscale[i13_ref(j, ii) - 1] = polpot_1.p3scale;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = mplpot_1.m4scale;
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
	    mscale[i15_ref(j, ii) - 1] = mplpot_1.m5scale;
	    pscale[i15_ref(j, ii) - 1] = polpot_1.p5scale;
	}
	i__2 = neigh_1.nelst[i__ - 1];
	for (kkk = 1; kkk <= i__2; ++kkk) {
	    k = elst_ref(kkk, i__);
	    kk = mpole_1.ipole[k - 1];
	    xr = atoms_1.x[kk - 1] - atoms_1.x[ii - 1];
	    yr = atoms_1.y[kk - 1] - atoms_1.y[ii - 1];
	    zr = atoms_1.z__[kk - 1] - atoms_1.z__[ii - 1];
	    image_(&xr, &yr, &zr);
	    r2 = xr * xr + yr * yr + zr * zr;
	    if (r2 <= shunt_1.cut2) {
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
		ukx = uind_ref(1, k);
		uky = uind_ref(2, k);
		ukz = uind_ref(3, k);

/*     calculate the error function damping terms */

		ralpha = ewald_1.aewald * r__;
		bn[0] = erfc_(&ralpha) / r__;
/* Computing 2nd power */
		d__1 = ewald_1.aewald;
		alsq2 = d__1 * d__1 * 2.;
		alsq2n = 0.;
		if (ewald_1.aewald > 0.) {
		    alsq2n = 1. / (ewald_1.aewald * 1.772453850905516027);
		}
/* Computing 2nd power */
		d__1 = ralpha;
		exp2a = exp(-(d__1 * d__1));
		for (m = 1; m <= 4; ++m) {
		    bfac = (doublereal) (m + m - 1);
		    alsq2n = alsq2 * alsq2n;
		    bn[m] = (bfac * bn[m - 1] + alsq2n * exp2a) / r2;
		}

/*     construct some intermediate quadrupole values */

		qix = qixx * xr + qixy * yr + qixz * zr;
		qiy = qixy * xr + qiyy * yr + qiyz * zr;
		qiz = qixz * xr + qiyz * yr + qizz * zr;
		qkx = qkxx * xr + qkxy * yr + qkxz * zr;
		qky = qkxy * xr + qkyy * yr + qkyz * zr;
		qkz = qkxz * xr + qkyz * yr + qkzz * zr;

/*     calculate the scalar products for permanent multipoles */

		sc[1] = dix * dkx + diy * dky + diz * dkz;
		sc[2] = dix * xr + diy * yr + diz * zr;
		sc[3] = dkx * xr + dky * yr + dkz * zr;
		sc[4] = qix * xr + qiy * yr + qiz * zr;
		sc[5] = qkx * xr + qky * yr + qkz * zr;
		sc[6] = qix * dkx + qiy * dky + qiz * dkz;
		sc[7] = qkx * dix + qky * diy + qkz * diz;
		sc[8] = qix * qkx + qiy * qky + qiz * qkz;
		sc[9] = (qixy * qkxy + qixz * qkxz + qiyz * qkyz) * 2. + qixx 
			* qkxx + qiyy * qkyy + qizz * qkzz;

/*     calculate the scalar products for polarization components */

		sci[1] = uix * dkx + dix * ukx + uiy * dky + diy * uky + uiz *
			 dkz + diz * ukz;
		sci[2] = uix * xr + uiy * yr + uiz * zr;
		sci[3] = ukx * xr + uky * yr + ukz * zr;
		sci[6] = qix * ukx + qiy * uky + qiz * ukz;
		sci[7] = qkx * uix + qky * uiy + qkz * uiz;

/*     calculate the gl functions for permanent multipoles */

		gl[0] = ci * ck;
		gl[1] = ck * sc[2] - ci * sc[3] + sc[1];
		gl[2] = ci * sc[5] + ck * sc[4] - sc[2] * sc[3] + (sc[6] - sc[
			7] + sc[9]) * 2.;
		gl[3] = sc[2] * sc[5] - sc[3] * sc[4] - sc[8] * 4.;
		gl[4] = sc[4] * sc[5];

/*     calculate the gl functions for polarization components */

		gli[0] = ck * sci[2] - ci * sci[3] + sci[1];
		gli[1] = (sci[6] - sci[7]) * 2. - sci[2] * sc[3] - sc[2] * 
			sci[3];
		gli[2] = sci[2] * sc[5] - sci[3] * sc[4];

/*     compute the energy contributions for this interaction */

		e = gl[0] * bn[0] + gl[1] * bn[1] + gl[2] * bn[2] + gl[3] * 
			bn[3] + gl[4] * bn[4];
		ei = gli[0] * bn[1] + gli[1] * bn[2] + gli[2] * bn[3];

/*     full real space energies needeed for scaled interactions */

		rr1 = 1. / r__;
		rr3 = rr1 / r2;
		rr5 = rr3 * 3. / r2;
		rr7 = rr5 * 5. / r2;
		rr9 = rr7 * 7. / r2;
		scale3 = pscale[kk - 1];
		scale5 = pscale[kk - 1];
		scale7 = pscale[kk - 1];
		damp = pdi * polar_1.pdamp[k - 1];
		if (damp != 0.) {
/* Computing MIN */
		    d__1 = pti, d__2 = polar_1.thole[k - 1];
		    pgamma = min(d__1,d__2);
/* Computing 3rd power */
		    d__1 = r__ / damp;
		    damp = -pgamma * (d__1 * (d__1 * d__1));
		    if (damp > -50.) {
			expdamp = exp(damp);
			scale3 *= 1. - expdamp;
			scale5 *= 1. - (1. - damp) * expdamp;
/* Computing 2nd power */
			d__1 = damp;
			scale7 *= 1. - (1. - damp + d__1 * d__1 * .6) * 
				expdamp;
		    }
		}
		efix = gl[0] * rr1 + gl[1] * rr3 + gl[2] * rr5 + gl[3] * rr7 
			+ gl[4] * rr9;
		eifix = gli[0] * rr3 * (1. - scale3) + gli[1] * rr5 * (1. - 
			scale5) + gli[2] * rr7 * (1. - scale7);

/*     apply the energy adjustments for scaled interactions */

		e -= efix * (1. - mscale[kk - 1]);
		ei -= eifix;

/*     increment the overall multipole and polarization energies */

		e = f * e;
		ei = f * .5 * ei;
		energi_1.em += e;
		energi_1.ep += ei;
	    }
	}

/*     reset interaction scaling coefficients for connected atoms */

	i__2 = couple_1.n12[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i12_ref(j, ii) - 1] = 1.;
	    pscale[i12_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n13[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i13_ref(j, ii) - 1] = 1.;
	    pscale[i13_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n14[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i14_ref(j, ii) - 1] = 1.;
	    pscale[i14_ref(j, ii) - 1] = 1.;
	}
	i__2 = couple_1.n15[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    mscale[i15_ref(j, ii) - 1] = 1.;
	    pscale[i15_ref(j, ii) - 1] = 1.;
	}
    }
    return 0;
} /* ereal0d_ */

#undef rpole_ref
#undef elst_ref
#undef uind_ref
#undef ip11_ref
#undef i15_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine emrecip  --  PME recip space multipole energy  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "emrecip" evaluates the reciprocal space portion of the particle */
/*     mesh Ewald energy due to atomic multipole interactions and */
/*     dipole polarizability */

/*     literature reference: */

/*     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate */
/*     Representation of Electrostatics in Classical Force Fields: */
/*     Efficient Implementation of Multipolar Interactions in */
/*     Biomolecular Simulations", Journal of Chemical Physics, 120, */
/*     73-87 (2004) */

/*     modifications for nonperiodic systems suggested by Tom Darden */
/*     during May 2007 */


/* Subroutine */ int emrecip_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal), cos(doublereal);

    /* Local variables */
    extern /* Subroutine */ int fftfront_(void);
    static doublereal a[9]	/* was [3][3] */, e;
    static integer i__, j, k;
    static doublereal h1, h2, h3;
    static integer k1, k2, k3, m1, m2, m3;
    static doublereal r1, r2, r3;
    extern /* Subroutine */ int grid_mpole__(doublereal *), fphi_mpole__(
	    doublereal *), cmp_to_fmp__(doublereal *, doublereal *);
    static integer nf1, nf2, nf3, nff;
    static doublereal cmp[250000]	/* was [10][25000] */, fmp[250000]	
	    /* was [10][25000] */, hsq, fphi[500000]	/* was [20][25000] */,
	     term;
    static integer ntot;
    extern /* Subroutine */ int bspline_fill__(void);
    static doublereal denom, fuind[75000]	/* was [3][25000] */, pterm;
    extern /* Subroutine */ int fftback_(void);
    static doublereal expterm, volterm;


#define a_ref(a_1,a_2) a[(a_2)*3 + a_1 - 4]
#define cmp_ref(a_1,a_2) cmp[(a_2)*10 + a_1 - 11]
#define fmp_ref(a_1,a_2) fmp[(a_2)*10 + a_1 - 11]
#define qfac_ref(a_1,a_2,a_3) pme_1.qfac[((a_3)*100 + (a_2))*100 + a_1 - \
10101]
#define fphi_ref(a_1,a_2) fphi[(a_2)*20 + a_1 - 21]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define recip_ref(a_1,a_2) boxes_1.recip[(a_2)*3 + a_1 - 4]
#define fuind_ref(a_1,a_2) fuind[(a_2)*3 + a_1 - 4]
#define qgrid_ref(a_1,a_2,a_3,a_4) pme_1.qgrid[(((a_4)*100 + (a_3))*100 + \
(a_2))*2 + a_1 - 20203]
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
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  ewald.i  --  parameters and options for Ewald summation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     aewald     Ewald convergence coefficient value (Ang-1) */
/*     boundary   Ewald boundary condition; none, tinfoil or vacuum */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  pme.i  --  values for particle mesh Ewald summation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxorder   maximum order of the B-spline approximation */
/*     maxprime   maximum number of prime factors of FFT dimension */
/*     maxtable   maximum size of the FFT table array */

/*     bsmod1     B-spline moduli along the a-axis direction */
/*     bsmod2     B-spline moduli along the b-axis direction */
/*     bsmod3     B-spline moduli along the c-axis direction */
/*     table      intermediate array used by the FFT calculation */
/*     qgrid      values on the particle mesh Ewald charge grid */
/*     qfac       prefactors for particle mesh Ewald charge grid */
/*     thetai1    B-spline coefficients along the a-axis */
/*     thetai2    B-spline coefficients along the b-axis */
/*     thetai3    B-spline coefficients along the c-axis */
/*     nfft1      number of grid points along the a-axis direction */
/*     nfft2      number of grid points along the b-axis direction */
/*     nfft3      number of grid points along the c-axis direction */
/*     bsorder    order of the PME B-spline approximation */
/*     iprime     prime factorization of each FFT dimension */
/*     igrid      initial Ewald charge grid values for B-spline */




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




/*     return if the Ewald coefficient is zero */

    if (ewald_1.aewald < 1e-6) {
	return 0;
    }

/*     copy the multipole moments into local storage areas */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cmp_ref(1, i__) = rpole_ref(1, i__);
	cmp_ref(2, i__) = rpole_ref(2, i__);
	cmp_ref(3, i__) = rpole_ref(3, i__);
	cmp_ref(4, i__) = rpole_ref(4, i__);
	cmp_ref(5, i__) = rpole_ref(5, i__);
	cmp_ref(6, i__) = rpole_ref(9, i__);
	cmp_ref(7, i__) = rpole_ref(13, i__);
	cmp_ref(8, i__) = rpole_ref(6, i__) * 2.;
	cmp_ref(9, i__) = rpole_ref(7, i__) * 2.;
	cmp_ref(10, i__) = rpole_ref(10, i__) * 2.;
    }

/*     compute the arrays of B-spline coefficients */

    if (! potent_1.use_polar__) {
	bspline_fill__();
    }

/*     convert Cartesian multipoles to fractional coordinates */

    cmp_to_fmp__(cmp, fmp);

/*     assign PME grid and perform 3-D FFT forward transform */

    grid_mpole__(fmp);
    fftfront_();

/*     make the scalar summation over reciprocal lattice */

    ntot = pme_1.nfft1 * pme_1.nfft2 * pme_1.nfft3;
/* Computing 2nd power */
    d__1 = 3.141592653589793238 / ewald_1.aewald;
    pterm = d__1 * d__1;
    volterm = boxes_1.volbox * 3.141592653589793238;
    nff = pme_1.nfft1 * pme_1.nfft2;
    nf1 = (pme_1.nfft1 + 1) / 2;
    nf2 = (pme_1.nfft2 + 1) / 2;
    nf3 = (pme_1.nfft3 + 1) / 2;
    i__1 = ntot - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k3 = i__ / nff + 1;
	j = i__ - (k3 - 1) * nff;
	k2 = j / pme_1.nfft1 + 1;
	k1 = j - (k2 - 1) * pme_1.nfft1 + 1;
	m1 = k1 - 1;
	m2 = k2 - 1;
	m3 = k3 - 1;
	if (k1 > nf1) {
	    m1 -= pme_1.nfft1;
	}
	if (k2 > nf2) {
	    m2 -= pme_1.nfft2;
	}
	if (k3 > nf3) {
	    m3 -= pme_1.nfft3;
	}
	r1 = (doublereal) m1;
	r2 = (doublereal) m2;
	r3 = (doublereal) m3;
	h1 = recip_ref(1, 1) * r1 + recip_ref(1, 2) * r2 + recip_ref(1, 3) * 
		r3;
	h2 = recip_ref(2, 1) * r1 + recip_ref(2, 2) * r2 + recip_ref(2, 3) * 
		r3;
	h3 = recip_ref(3, 1) * r1 + recip_ref(3, 2) * r2 + recip_ref(3, 3) * 
		r3;
	hsq = h1 * h1 + h2 * h2 + h3 * h3;
	term = -pterm * hsq;
	expterm = 0.;
	if (term > -50.) {
	    denom = volterm * hsq * pme_1.bsmod1[k1 - 1] * pme_1.bsmod2[k2 - 
		    1] * pme_1.bsmod3[k3 - 1];
	    expterm = exp(term) / denom;
	    if (! bound_1.use_bounds__) {
		expterm *= 1. - cos(boxes_1.xbox * 3.141592653589793238 * 
			sqrt(hsq));
	    } else if (boxes_1.octahedron) {
		if ((m1 + m2 + m3) % 2 != 0) {
		    expterm = 0.;
		}
	    }
	}
	qfac_ref(k1, k2, k3) = expterm;
    }

/*     account for the zeroth grid point for a finite system */

    qfac_ref(1, 1, 1) = 0.;
    if (! bound_1.use_bounds__) {
	expterm = 1.5707963267948966 / boxes_1.xbox;
	qfac_ref(1, 1, 1) = expterm;
    }

/*     complete the transformation of the charge grid */

    i__1 = pme_1.nfft3;
    for (k = 1; k <= i__1; ++k) {
	i__2 = pme_1.nfft2;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = pme_1.nfft1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		term = qfac_ref(i__, j, k);
		qgrid_ref(1, i__, j, k) = term * qgrid_ref(1, i__, j, k);
		qgrid_ref(2, i__, j, k) = term * qgrid_ref(2, i__, j, k);
	    }
	}
    }

/*     perform 3-D FFT backward transform and get potential */

    fftback_();
    fphi_mpole__(fphi);

/*     sum over multipoles and increment total multipole energy */

    e = 0.;
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (k = 1; k <= 10; ++k) {
	    e += fmp_ref(k, i__) * fphi_ref(k, i__);
	}
    }
    e = chgpot_1.electric * .5 * e;
    energi_1.em += e;

/*     convert Cartesian induced dipoles to fractional coordinates */

    if (potent_1.use_polar__) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    a_ref(1, i__) = (doublereal) pme_1.nfft1 * recip_ref(i__, 1);
	    a_ref(2, i__) = (doublereal) pme_1.nfft2 * recip_ref(i__, 2);
	    a_ref(3, i__) = (doublereal) pme_1.nfft3 * recip_ref(i__, 3);
	}
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (k = 1; k <= 3; ++k) {
		fuind_ref(k, i__) = a_ref(k, 1) * uind_ref(1, i__) + a_ref(k, 
			2) * uind_ref(2, i__) + a_ref(k, 3) * uind_ref(3, i__)
			;
	    }
	}

/*     sum over induced dipoles and increment total induced energy */

	e = 0.;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (k = 1; k <= 3; ++k) {
		e += fuind_ref(k, i__) * fphi_ref(k + 1, i__);
	    }
	}
	e = chgpot_1.electric * .5 * e;
	energi_1.ep += e;
    }
    return 0;
} /* emrecip_ */

#undef rpole_ref
#undef qgrid_ref
#undef fuind_ref
#undef recip_ref
#undef uind_ref
#undef fphi_ref
#undef qfac_ref
#undef fmp_ref
#undef cmp_ref
#undef a_ref


