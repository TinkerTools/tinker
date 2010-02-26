/* orbital.f -- translated by f2c (version 20050501).
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
    doublereal bk[50000], bl[50000];
    integer nbond, ibnd[100000]	/* was [2][50000] */;
} bond_;

#define bond_1 bond_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

struct {
    integer norbit, iorbit[100], reorbit, piperp[300]	/* was [3][100] */, 
	    nbpi, ibpi[600]	/* was [3][200] */, ntpi, itpi[800]	/* 
	    was [2][400] */;
    logical listpi[25000];
} piorbs_;

#define piorbs_1 piorbs_

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
    doublereal tors1[400000]	/* was [4][100000] */, tors2[400000]	/* 
	    was [4][100000] */, tors3[400000]	/* was [4][100000] */, tors4[
	    400000]	/* was [4][100000] */, tors5[400000]	/* was [4][
	    100000] */, tors6[400000]	/* was [4][100000] */;
    integer ntors, itors[400000]	/* was [4][100000] */;
} tors_;

#define tors_1 tors_

struct {
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine orbital  --  setup for pisystem calculation  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "orbital" finds and organizes lists of atoms in a pisystem, */
/*     bonds connecting pisystem atoms and torsions whose two */
/*     central atoms are both pisystem atoms */


/* Subroutine */ int orbital_(void)
{
    /* Format strings */
    static char fmt_20[] = "(\002 ORBITAL  --  Too many Pi-Orbitals;\002,"
	    "\002 Increase MAXPI\002)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k, m, ib, ic, iorb, jorb, next;
    extern /* Subroutine */ int fatal_(void);
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int piplane_(void);
    static integer piatoms[100];
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___7 = { 1, string, 1, 0, 120, 1 };
    static cilist io___10 = { 0, 0, 0, fmt_20, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define ibpi_ref(a_1,a_2) piorbs_1.ibpi[(a_2)*3 + a_1 - 4]
#define itpi_ref(a_1,a_2) piorbs_1.itpi[(a_2)*2 + a_1 - 3]
#define itors_ref(a_1,a_2) tors_1.itors[(a_2)*4 + a_1 - 5]
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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  piorbs.i  --  conjugated system in the current structure  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     norbit    total number of pisystem orbitals in the system */
/*     iorbit    numbers of the atoms containing pisystem orbitals */
/*     reorbit   number of evaluations between orbital updates */
/*     piperp    atoms defining a normal plane to each orbital */
/*     nbpi      total number of bonds affected by the pisystem */
/*     ibpi      bond and piatom numbers for each pisystem bond */
/*     ntpi      total number of torsions affected by the pisystem */
/*     itpi      torsion and pibond numbers for each pisystem torsion */
/*     listpi    atom list indicating whether each atom has an orbital */




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




/*     set the default values for the pisystem variables */

    for (i__ = 1; i__ <= 100; ++i__) {
	piatoms[i__ - 1] = 0;
    }
    piorbs_1.reorbit = 1;

/*     check the keywords for any lists of pisystem atoms */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "PISYSTEM ", (ftnlen)9, (ftnlen)9) == 0) {
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___7);
	    if (i__2 != 0) {
		goto L10;
	    }
	    for (k = 1; k <= 100; ++k) {
		i__2 = do_lio(&c__3, &c__1, (char *)&piatoms[k - 1], (ftnlen)
			sizeof(integer));
		if (i__2 != 0) {
		    goto L10;
		}
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	}
L10:
	;
    }

/*     quit if no pisystem was found for consideration */

    if (piatoms[0] == 0) {
	potent_1.use_orbit__ = FALSE_;
	return 0;
    } else {
	potent_1.use_orbit__ = TRUE_;
    }

/*     organize and make lists of the pisystem atoms */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	piorbs_1.listpi[i__ - 1] = FALSE_;
    }
    i__ = 1;
    while(piatoms[i__ - 1] != 0) {
	if (piatoms[i__ - 1] > 0) {
	    piorbs_1.listpi[piatoms[i__ - 1] - 1] = TRUE_;
	    ++i__;
	} else {
	    i__1 = piatoms[i__];
	    for (j = -piatoms[i__ - 1]; j <= i__1; ++j) {
		piorbs_1.listpi[j - 1] = TRUE_;
	    }
	    i__ += 2;
	}
    }
    piorbs_1.norbit = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (piorbs_1.listpi[i__ - 1]) {
	    ++piorbs_1.norbit;
	    piorbs_1.iorbit[piorbs_1.norbit - 1] = i__;
	}
    }

/*     quit if the molecule contains too many piorbitals */

    if (piorbs_1.norbit > 100) {
	io___10.ciunit = iounit_1.iout;
	s_wsfe(&io___10);
	e_wsfe();
	fatal_();
    }

/*     find three atoms which define a plane */
/*     perpendicular to each orbital */

    piplane_();

/*     find and store the pisystem bonds */

    piorbs_1.nbpi = 0;
    i__1 = piorbs_1.norbit - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iorb = piorbs_1.iorbit[i__ - 1];
	i__2 = piorbs_1.norbit;
	for (j = i__; j <= i__2; ++j) {
	    jorb = piorbs_1.iorbit[j - 1];
	    i__3 = couple_1.n12[iorb - 1];
	    for (k = 1; k <= i__3; ++k) {
		if (i12_ref(k, iorb) == jorb) {
		    ++piorbs_1.nbpi;
		    i__4 = bond_1.nbond;
		    for (m = 1; m <= i__4; ++m) {
			if (iorb == ibnd_ref(1, m) && jorb == ibnd_ref(2, m)) 
				{
			    ibpi_ref(1, piorbs_1.nbpi) = m;
			    ibpi_ref(2, piorbs_1.nbpi) = i__;
			    ibpi_ref(3, piorbs_1.nbpi) = j;
			    goto L30;
			}
		    }
L30:
		    ;
		}
	    }
	}
    }

/*     find and store the pisystem torsions */

    piorbs_1.ntpi = 0;
    i__1 = tors_1.ntors;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ib = itors_ref(2, i__);
	ic = itors_ref(3, i__);
	if (piorbs_1.listpi[ib - 1] && piorbs_1.listpi[ic - 1]) {
	    i__2 = piorbs_1.nbpi;
	    for (j = 1; j <= i__2; ++j) {
		k = ibpi_ref(1, j);
		if (ib == ibnd_ref(1, k) && ic == ibnd_ref(2, k) || ib == 
			ibnd_ref(2, k) && ic == ibnd_ref(1, k)) {
		    ++piorbs_1.ntpi;
		    itpi_ref(1, piorbs_1.ntpi) = i__;
		    itpi_ref(2, piorbs_1.ntpi) = j;
		    goto L40;
		}
	    }
L40:
	    ;
	}
    }
    return 0;
} /* orbital_ */

#undef keyline_ref
#undef itors_ref
#undef itpi_ref
#undef ibpi_ref
#undef ibnd_ref
#undef i12_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine piplane  --  plane perpendicular to orbital  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "piplane" selects the three atoms which specify the plane */
/*     perpendicular to each p-orbital; the current version will */
/*     fail in certain situations, including ketenes, allenes, */
/*     and isolated or adjacent triple bonds */


/* Subroutine */ int piplane_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 PIPLANE  --  Failure to Define a\002,"
	    "\002 p-Orbital Plane for Atom\002,i6)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j, beta;
    static logical done;
    static integer iorb, gamma, alpha;
    extern /* Subroutine */ int fatal_(void);
    static integer trial, attach, atmnum;

    /* Fortran I/O blocks */
    static cilist io___26 = { 0, 0, 0, fmt_10, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define piperp_ref(a_1,a_2) piorbs_1.piperp[(a_2)*3 + a_1 - 4]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  piorbs.i  --  conjugated system in the current structure  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     norbit    total number of pisystem orbitals in the system */
/*     iorbit    numbers of the atoms containing pisystem orbitals */
/*     reorbit   number of evaluations between orbital updates */
/*     piperp    atoms defining a normal plane to each orbital */
/*     nbpi      total number of bonds affected by the pisystem */
/*     ibpi      bond and piatom numbers for each pisystem bond */
/*     ntpi      total number of torsions affected by the pisystem */
/*     itpi      torsion and pibond numbers for each pisystem torsion */
/*     listpi    atom list indicating whether each atom has an orbital */




/*     for each pisystem atom, find a set of atoms which define */
/*     the p-orbital's plane based on piatom's atomic number and */
/*     the number and type of attached atoms */

    i__1 = piorbs_1.norbit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iorb = piorbs_1.iorbit[i__ - 1];
	attach = couple_1.n12[iorb - 1];
	atmnum = atmtyp_1.atomic[iorb - 1];
	done = FALSE_;

/*     most common case of an atom bonded to three atoms */

	if (attach == 3) {
	    piperp_ref(1, i__) = i12_ref(1, iorb);
	    piperp_ref(2, i__) = i12_ref(2, iorb);
	    piperp_ref(3, i__) = i12_ref(3, iorb);
	    done = TRUE_;

/*     any non-alkyne atom bonded to exactly two atoms */

	} else if (attach == 2 && atmnum != 6) {
	    piperp_ref(1, i__) = iorb;
	    piperp_ref(2, i__) = i12_ref(1, iorb);
	    piperp_ref(3, i__) = i12_ref(2, iorb);
	    done = TRUE_;

/*     atom bonded to four different atoms (usually two lone */
/*     pairs and two "real" atoms); use the "real" atoms */

	} else if (attach == 4) {
	    piperp_ref(1, i__) = iorb;
	    i__2 = couple_1.n12[iorb - 1];
	    for (j = 1; j <= i__2; ++j) {
		trial = i12_ref(j, iorb);
		if (atmtyp_1.atomic[trial - 1] != 0) {
		    if (piperp_ref(2, i__) == 0) {
			piperp_ref(2, i__) = trial;
		    } else {
			piperp_ref(3, i__) = trial;
			done = TRUE_;
		    }
		}
	    }

/*     "carbonyl"-type oxygen atom, third atom is any atom */
/*     attached to the "carbonyl carbon"; fails for ketenes */

	} else if (attach == 1 && atmnum == 8) {
	    alpha = i12_ref(1, iorb);
	    beta = i12_ref(1, alpha);
	    if (beta == iorb) {
		beta = i12_ref(2, alpha);
	    }
	    piperp_ref(1, i__) = iorb;
	    piperp_ref(2, i__) = alpha;
	    piperp_ref(3, i__) = beta;
	    done = TRUE_;

/*     an sp nitrogen atom, third atom must be a gamma atom */

	} else if (attach == 1 && atmnum == 7) {
	    alpha = i12_ref(1, iorb);
	    i__2 = couple_1.n12[alpha - 1];
	    for (j = 1; j <= i__2; ++j) {
		trial = i12_ref(j, alpha);
		if (trial != iorb && piorbs_1.listpi[trial - 1] && 
			couple_1.n12[trial - 1] == 3) {
		    beta = trial;
		    done = TRUE_;
		}
	    }
	    gamma = i12_ref(1, beta);
	    if (gamma == alpha) {
		gamma = i12_ref(2, beta);
	    }
	    piperp_ref(1, i__) = iorb;
	    piperp_ref(2, i__) = alpha;
	    piperp_ref(3, i__) = gamma;

/*     an sp carbon atom; third atom must be an atom attached */
/*     to the non-sp piatom bonded to the original carbon */

	} else if (attach == 2 && atmnum == 6) {
	    alpha = i12_ref(1, iorb);
	    if (couple_1.n12[alpha - 1] == 2 && atmtyp_1.atomic[alpha - 1] == 
		    6 || couple_1.n12[alpha - 1] == 1 && atmtyp_1.atomic[
		    alpha - 1] == 7) {
		alpha = i12_ref(2, iorb);
	    }
	    i__2 = couple_1.n12[iorb - 1];
	    for (j = 1; j <= i__2; ++j) {
		trial = i12_ref(j, iorb);
		if (trial != iorb && trial != alpha && piorbs_1.listpi[trial 
			- 1] && couple_1.n12[trial - 1] == 3) {
		    beta = trial;
		    done = TRUE_;
		}
	    }
	    i__2 = couple_1.n12[alpha - 1];
	    for (j = 1; j <= i__2; ++j) {
		trial = i12_ref(j, alpha);
		if (trial != iorb && trial != alpha && piorbs_1.listpi[trial 
			- 1] && couple_1.n12[trial - 1] == 3) {
		    beta = trial;
		    done = TRUE_;
		}
	    }
	    gamma = i12_ref(1, beta);
	    if (gamma == iorb || gamma == alpha) {
		gamma = i12_ref(2, beta);
	    }
	    piperp_ref(1, i__) = iorb;
	    piperp_ref(2, i__) = alpha;
	    piperp_ref(3, i__) = gamma;
	}

/*     quit if the p-orbital plane remains undefined */

	if (! done) {
	    io___26.ciunit = iounit_1.iout;
	    s_wsfe(&io___26);
	    do_fio(&c__1, (char *)&iorb, (ftnlen)sizeof(integer));
	    e_wsfe();
	    fatal_();
	}
    }
    return 0;
} /* piplane_ */

#undef piperp_ref
#undef i12_ref


