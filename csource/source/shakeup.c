/* shakeup.f -- translated by f2c (version 20050501).
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
    doublereal ak[75000], anat[75000], afld[75000];
    integer nangle, iang[300000]	/* was [4][75000] */;
} angle_;

#define angle_1 angle_

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
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

struct {
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

struct {
    doublereal krat[25000];
    integer nrat, nratx, irat[50000]	/* was [2][25000] */, iratx[25000], 
	    kratx[25000];
    logical ratimage[25000], use_rattle__;
} shake_;

#define shake_1 shake_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    integer nring3, iring3[30000]	/* was [3][10000] */, nring4, iring4[
	    40000]	/* was [4][10000] */, nring5, iring5[50000]	/* 
	    was [5][10000] */, nring6, iring6[60000]	/* was [6][10000] */;
} ring_;

#define ring_1 ring_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine shakeup  --  setup of rattle constraints  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "shakeup" initializes any holonomic constraints for use */
/*     with the rattle algorithm during molecular dynamics */


/* Subroutine */ int shakeup_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double cos(doublereal), sqrt(doublereal);
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void);

    /* Local variables */
    extern /* Subroutine */ int chkangle_(integer *, integer *, integer *);
    static integer i__, j, k, m, ia, ib, nh, ic, ja, jb;
    static doublereal rab, rac, rbc;
    static logical done;
    static integer next, ilist;
    static char record[120];
    static doublereal cosine;
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120], rattyp[9];
    extern /* Subroutine */ int getnumb_(char *, integer *, integer *, ftnlen)
	    , getword_(char *, char *, integer *, ftnlen, ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___19 = { 1, string, 1, 0, 120, 1 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define irat_ref(a_1,a_2) shake_1.irat[(a_2)*2 + a_1 - 3]
#define bndlist_ref(a_1,a_2) atmlst_1.bndlist[(a_2)*8 + a_1 - 9]
#define anglist_ref(a_1,a_2) atmlst_1.anglist[(a_2)*28 + a_1 - 29]
#define keyline_ref(a_0,a_1) &keys_1.keyline[(a_1)*120 + a_0 - 120]



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
/*     ##  angle.i  --  bond angles within the current structure  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     ak       harmonic angle force constant (kcal/mole/rad**2) */
/*     anat     ideal bond angle or phase shift angle (degrees) */
/*     afld     periodicity for Fourier bond angle term */
/*     nangle   total number of bond angles in the system */
/*     iang     numbers of the atoms in each bond angle */




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
/*     ##  keys.i  --  contents of current keyword parameter file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nkey      number of nonblank lines in the keyword file */
/*     keyline   contents of each individual keyword file line */




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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  shake.i  --  definition of Shake/Rattle constraints  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     krat         ideal distance value for rattle constraint */
/*     nrat         number of rattle distance constraints to apply */
/*     nratx        number of atom group spatial constraints to apply */
/*     irat         atom numbers of atoms in a rattle constraint */
/*     iratx        group number of group in a spatial constraint */
/*     kratx        spatial constraint type (1=plane, 2=line, 3=point) */
/*     ratimage     flag to use minimum image for rattle constraint */
/*     use_rattle   logical flag to set use of rattle contraints */




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




/*     zero out the number of distance and spatial constraints */

    shake_1.nrat = 0;
    shake_1.nratx = 0;
    shake_1.use_rattle__ = TRUE_;

/*     process keywords containing generic constraint options */

    i__1 = keys_1.nkey;
    for (k = 1; k <= i__1; ++k) {
	next = 1;
	s_copy(record, keyline_ref(0, k), (ftnlen)120, (ftnlen)120);
	upcase_(record, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);

/*     get the distance constraint types for the rattle method */

	if (s_cmp(keyword, "RATTLE ", (ftnlen)7, (ftnlen)7) == 0) {
	    getword_(record, rattyp, &next, (ftnlen)120, (ftnlen)9);

/*     constrain all bond lengths at their ideal values */

	    if (s_cmp(rattyp, "BONDS", (ftnlen)5, (ftnlen)5) == 0) {
		i__2 = bond_1.nbond;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ia = ibnd_ref(1, i__);
		    ib = ibnd_ref(2, i__);
		    if (usage_1.use[ia - 1] || usage_1.use[ib - 1]) {
			++shake_1.nrat;
			irat_ref(1, shake_1.nrat) = ia;
			irat_ref(2, shake_1.nrat) = ib;
			shake_1.krat[shake_1.nrat - 1] = bond_1.bl[i__ - 1];
		    }
		}

/*     constrain bonds and independent angles at ideal values */

	    } else if (s_cmp(rattyp, "ANGLES", (ftnlen)6, (ftnlen)6) == 0) {
		i__2 = bond_1.nbond;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ia = ibnd_ref(1, i__);
		    ib = ibnd_ref(2, i__);
		    if (usage_1.use[ia - 1] || usage_1.use[ib - 1]) {
			++shake_1.nrat;
			irat_ref(1, shake_1.nrat) = ia;
			irat_ref(2, shake_1.nrat) = ib;
			shake_1.krat[shake_1.nrat - 1] = bond_1.bl[i__ - 1];
		    }
		}
		i__2 = atoms_1.n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    if (couple_1.n12[i__ - 1] > 1) {
			i__3 = (couple_1.n12[i__ - 1] << 1) - 3;
			for (j = 1; j <= i__3; ++j) {
			    ilist = anglist_ref(j, i__);
			    ia = iang_ref(1, ilist);
			    ib = iang_ref(2, ilist);
			    ic = iang_ref(3, ilist);
			    if (usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
				    usage_1.use[ic - 1]) {
				i__4 = couple_1.n12[ib - 1];
				for (m = 1; m <= i__4; ++m) {
				    if (i12_ref(m, ib) == ia) {
					rab = bond_1.bl[bndlist_ref(m, ib) - 
						1];
				    } else if (i12_ref(m, ib) == ic) {
					rbc = bond_1.bl[bndlist_ref(m, ib) - 
						1];
				    }
				}
				cosine = cos(angle_1.anat[ilist - 1] / 
					57.29577951308232088);
				rac = sqrt(rab * rab + rbc * rbc - rab * 2. * 
					rbc * cosine);
				++shake_1.nrat;
				irat_ref(1, shake_1.nrat) = ia;
				irat_ref(2, shake_1.nrat) = ic;
				shake_1.krat[shake_1.nrat - 1] = rac;
				chkangle_(&ia, &ib, &ic);
			    }
			}
		    }
		}

/*     fix bond length in diatomics to give a rigid molecule */

	    } else if (s_cmp(rattyp, "DIATOMIC", (ftnlen)8, (ftnlen)8) == 0) {
		i__2 = bond_1.nbond;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ia = ibnd_ref(1, i__);
		    ib = ibnd_ref(2, i__);
		    if (couple_1.n12[ia - 1] == 1 && couple_1.n12[ib - 1] == 
			    1) {
			if (usage_1.use[ia - 1] || usage_1.use[ib - 1]) {
			    ++shake_1.nrat;
			    irat_ref(1, shake_1.nrat) = ia;
			    irat_ref(2, shake_1.nrat) = ib;
			    shake_1.krat[shake_1.nrat - 1] = bond_1.bl[i__ - 
				    1];
			}
		    }
		}

/*     fix bonds and angle in triatomics to give a rigid molecule */

	    } else if (s_cmp(rattyp, "TRIATOMIC", (ftnlen)9, (ftnlen)9) == 0) 
		    {
		i__2 = angle_1.nangle;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ia = iang_ref(1, i__);
		    ib = iang_ref(2, i__);
		    ic = iang_ref(3, i__);
		    if (couple_1.n12[ia - 1] + couple_1.n12[ib - 1] + 
			    couple_1.n12[ic - 1] == 4) {
			rab = bond_1.bl[bndlist_ref(1, ia) - 1];
			rbc = bond_1.bl[bndlist_ref(1, ic) - 1];
			cosine = cos(angle_1.anat[i__ - 1] / 
				57.29577951308232088);
/* Computing 2nd power */
			d__1 = rab;
/* Computing 2nd power */
			d__2 = rbc;
			rac = sqrt(d__1 * d__1 + d__2 * d__2 - rab * 2. * rbc 
				* cosine);
			if (usage_1.use[ia - 1] || usage_1.use[ib - 1]) {
			    ++shake_1.nrat;
			    irat_ref(1, shake_1.nrat) = ia;
			    irat_ref(2, shake_1.nrat) = ib;
			    shake_1.krat[shake_1.nrat - 1] = rab;
			}
			if (usage_1.use[ib - 1] || usage_1.use[ic - 1]) {
			    ++shake_1.nrat;
			    irat_ref(1, shake_1.nrat) = ib;
			    irat_ref(2, shake_1.nrat) = ic;
			    shake_1.krat[shake_1.nrat - 1] = rbc;
			}
			if (usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
				usage_1.use[ic - 1]) {
			    ++shake_1.nrat;
			    irat_ref(1, shake_1.nrat) = ia;
			    irat_ref(2, shake_1.nrat) = ic;
			    shake_1.krat[shake_1.nrat - 1] = rac;
			}
		    }
		}

/*     fix bonds and angles of each water to give a rigid molecule */

	    } else if (s_cmp(rattyp, "WATER", (ftnlen)5, (ftnlen)5) == 0) {
		i__2 = atoms_1.n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    nh = 0;
		    if (atmtyp_1.atomic[i__ - 1] == 8) {
			i__3 = couple_1.n12[i__ - 1];
			for (j = 1; j <= i__3; ++j) {
			    if (atmtyp_1.atomic[i12_ref(j, i__) - 1] == 1) {
				++nh;
			    }
			}
		    }
		    if (nh >= 2) {
			i__3 = couple_1.n12[i__ - 1];
			for (j = 1; j <= i__3; ++j) {
			    ilist = bndlist_ref(j, i__);
			    ia = ibnd_ref(1, ilist);
			    ib = ibnd_ref(2, ilist);
			    if (usage_1.use[ia - 1] || usage_1.use[ib - 1]) {
				++shake_1.nrat;
				irat_ref(1, shake_1.nrat) = ia;
				irat_ref(2, shake_1.nrat) = ib;
				shake_1.krat[shake_1.nrat - 1] = bond_1.bl[
					ilist - 1];
			    }
			}
			i__3 = (couple_1.n12[i__ - 1] << 1) - 3;
			for (j = 1; j <= i__3; ++j) {
			    ilist = anglist_ref(j, i__);
			    ia = iang_ref(1, ilist);
			    ib = iang_ref(2, ilist);
			    ic = iang_ref(3, ilist);
			    if (usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
				    usage_1.use[ic - 1]) {
				i__4 = couple_1.n12[ib - 1];
				for (m = 1; m <= i__4; ++m) {
				    if (i12_ref(m, ib) == ia) {
					rab = bond_1.bl[bndlist_ref(m, ib) - 
						1];
				    } else if (i12_ref(m, ib) == ic) {
					rbc = bond_1.bl[bndlist_ref(m, ib) - 
						1];
				    }
				}
				cosine = cos(angle_1.anat[ilist - 1] / 
					57.29577951308232088);
				rac = sqrt(rab * rab + rbc * rbc - rab * 2. * 
					rbc * cosine);
				++shake_1.nrat;
				irat_ref(1, shake_1.nrat) = ia;
				irat_ref(2, shake_1.nrat) = ic;
				shake_1.krat[shake_1.nrat - 1] = rac;
				chkangle_(&ia, &ib, &ic);
			    }
			}
		    }
		}

/*     fix all bonds to hydrogen atoms at their ideal length */

	    } else {
		i__2 = bond_1.nbond;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ia = ibnd_ref(1, i__);
		    ib = ibnd_ref(2, i__);
		    if (atmtyp_1.atomic[ia - 1] == 1 || atmtyp_1.atomic[ib - 
			    1] == 1) {
			if (usage_1.use[ia - 1] || usage_1.use[ib - 1]) {
			    ++shake_1.nrat;
			    irat_ref(1, shake_1.nrat) = ia;
			    irat_ref(2, shake_1.nrat) = ib;
			    shake_1.krat[shake_1.nrat - 1] = bond_1.bl[i__ - 
				    1];
			}
		    }
		}
	    }
	}
    }

/*     process keywords containing specific distance constraints */

    i__1 = keys_1.nkey;
    for (k = 1; k <= i__1; ++k) {
	next = 1;
	s_copy(record, keyline_ref(0, k), (ftnlen)120, (ftnlen)120);
	upcase_(record, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	if (s_cmp(keyword, "RATTLE-DISTANCE ", (ftnlen)16, (ftnlen)16) == 0) {
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    rab = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___19);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&rab, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
L10:
	    if (rab == 0.) {
		i__2 = couple_1.n12[ia - 1];
		for (i__ = 1; i__ <= i__2; ++i__) {
		    if (i12_ref(i__, ia) == ib) {
			rab = bond_1.bl[bndlist_ref(i__, ia) - 1];
		    }
		}
	    }
	    if (rab == 0.) {
/* Computing 2nd power */
		d__1 = atoms_1.x[ia - 1] - atoms_1.x[ib - 1];
/* Computing 2nd power */
		d__2 = atoms_1.y[ia - 1] - atoms_1.y[ib - 1];
/* Computing 2nd power */
		d__3 = atoms_1.z__[ia - 1] - atoms_1.z__[ib - 1];
		rab = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	    }
	    done = FALSE_;
	    i__2 = shake_1.nrat;
	    for (j = 1; j <= i__2; ++j) {
		ja = irat_ref(1, j);
		jb = irat_ref(2, j);
		if (ia == ja && ib == jb || ia == jb && ib == ja) {
		    done = TRUE_;
		    shake_1.krat[j - 1] = rab;
		}
	    }
	    if (! done) {
		++shake_1.nrat;
		irat_ref(1, shake_1.nrat) = ia;
		irat_ref(2, shake_1.nrat) = ib;
		shake_1.krat[shake_1.nrat - 1] = rab;
	    }

/*     process keywords containing atom group spatial constraints */

	} else if (s_cmp(keyword, "RATTLE-PLANE ", (ftnlen)13, (ftnlen)13) == 
		0) {
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    ++shake_1.nratx;
	    shake_1.iratx[shake_1.nratx - 1] = ia;
	    shake_1.kratx[shake_1.nratx - 1] = 1;
	} else if (s_cmp(keyword, "RATTLE-LINE ", (ftnlen)12, (ftnlen)12) == 
		0) {
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    ++shake_1.nratx;
	    shake_1.iratx[shake_1.nratx - 1] = ia;
	    shake_1.kratx[shake_1.nratx - 1] = 2;
	} else if (s_cmp(keyword, "RATTLE-ORIGIN ", (ftnlen)14, (ftnlen)14) ==
		 0) {
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    ++shake_1.nratx;
	    shake_1.iratx[shake_1.nratx - 1] = ia;
	    shake_1.kratx[shake_1.nratx - 1] = 3;
	}
    }

/*     find and remove any duplicate distance constraints */

    i__1 = shake_1.nrat - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = irat_ref(1, i__);
	ib = irat_ref(2, i__);
	i__2 = shake_1.nrat;
	for (j = i__ + 1; j <= i__2; ++j) {
	    ja = irat_ref(1, j);
	    jb = irat_ref(2, j);
	    if (ia == ja && ib == jb || ia == jb && ib == ja) {
		shake_1.krat[j - 1] = -1.;
	    }
	}
    }
    k = shake_1.nrat;
    for (i__ = k; i__ >= 1; --i__) {
	if (shake_1.krat[i__ - 1] < 0.) {
	    i__1 = k - 1;
	    for (j = i__; j <= i__1; ++j) {
		irat_ref(1, j) = irat_ref(1, j + 1);
		irat_ref(2, j) = irat_ref(2, j + 1);
		shake_1.krat[j - 1] = shake_1.krat[j];
	    }
	    --shake_1.nrat;
	}
    }

/*     set flag to apply minimum image to intermolecular constraints */

    i__1 = shake_1.nrat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = irat_ref(1, i__);
	ib = irat_ref(2, i__);
	if (bound_1.use_bounds__ && molcul_1.molcule[ia - 1] != 
		molcul_1.molcule[ib - 1]) {
	    shake_1.ratimage[i__ - 1] = TRUE_;
	} else {
	    shake_1.ratimage[i__ - 1] = FALSE_;
	}
    }

/*     if no rattle constraints are used, turn off its use */

    if (shake_1.nrat == 0 && shake_1.nratx == 0) {
	shake_1.use_rattle__ = FALSE_;
    }
    return 0;
} /* shakeup_ */

#undef keyline_ref
#undef anglist_ref
#undef bndlist_ref
#undef irat_ref
#undef iang_ref
#undef ibnd_ref
#undef i12_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine chkangle  --  eliminate redundant constraints  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "chkangle" tests angles to be constrained for their presence */
/*     in small rings and removes constraints that are redundant */

/*     note this version correctly handles isolated small rings, */
/*     but should remove one additional redundant constraint for */
/*     each ring fusion */


/* Subroutine */ int chkangle_(integer *ia, integer *ib, integer *ic)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, id, ie, imin;
    static logical remove;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  ring.i  --  number and location of small ring structures  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     nring3   total number of 3-membered rings in the system */
/*     iring3   numbers of the atoms involved in each 3-ring */
/*     nring4   total number of 4-membered rings in the system */
/*     iring4   numbers of the atoms involved in each 4-ring */
/*     nring5   total number of 5-membered rings in the system */
/*     iring5   numbers of the atoms involved in each 5-ring */
/*     nring6   total number of 6-membered rings in the system */
/*     iring6   numbers of the atoms involved in each 6-ring */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  shake.i  --  definition of Shake/Rattle constraints  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     krat         ideal distance value for rattle constraint */
/*     nrat         number of rattle distance constraints to apply */
/*     nratx        number of atom group spatial constraints to apply */
/*     irat         atom numbers of atoms in a rattle constraint */
/*     iratx        group number of group in a spatial constraint */
/*     kratx        spatial constraint type (1=plane, 2=line, 3=point) */
/*     ratimage     flag to use minimum image for rattle constraint */
/*     use_rattle   logical flag to set use of rattle contraints */




/*     all internal angles of 3-membered rings are redundant */

    remove = FALSE_;
    if (ring_1.nring3 != 0) {
	i__1 = couple_1.n12[*ia - 1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j = i12_ref(i__, *ia);
	    if (j == *ic) {
		remove = TRUE_;
	    }
	}
    }

/*     for 4-membered rings, two internal angles are redundant */

    if (ring_1.nring4 != 0) {
	i__1 = couple_1.n12[*ia - 1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    id = i12_ref(i__, *ia);
	    if (id != *ib) {
		i__2 = couple_1.n12[id - 1];
		for (j = 1; j <= i__2; ++j) {
		    k = i12_ref(j, id);
		    if (k == *ic) {
/* Computing MIN */
			i__3 = min(*ia,*ib), i__3 = min(i__3,*ic);
			imin = min(i__3,id);
			if (*ib == imin) {
			    remove = TRUE_;
			}
			if (id == imin) {
			    remove = TRUE_;
			}
		    }
		}
	    }
	}
    }

/*     for 5-membered rings, one internal angle is redundant */

    if (ring_1.nring5 != 0) {
	i__1 = couple_1.n12[*ia - 1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    id = i12_ref(i__, *ia);
	    if (id != *ib && id != *ic) {
		i__2 = couple_1.n12[*ic - 1];
		for (j = 1; j <= i__2; ++j) {
		    ie = i12_ref(j, *ic);
		    if (ie != *ib && ie != *ia) {
			i__3 = couple_1.n12[id - 1];
			for (k = 1; k <= i__3; ++k) {
			    if (i12_ref(k, id) == ie) {
/* Computing MIN */
				i__4 = min(*ia,*ib), i__4 = min(i__4,*ic), 
					i__4 = min(i__4,id);
				imin = min(i__4,ie);
				if (*ib == imin) {
				    remove = TRUE_;
				}
			    }
			}
		    }
		}
	    }
	}
    }

/*     if the constraint if redundant, remove it from the list */

    if (remove) {
	--shake_1.nrat;
    }
    return 0;
} /* chkangle_ */

#undef i12_ref


