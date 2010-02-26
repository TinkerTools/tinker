/* rings.f -- translated by f2c (version 20050501).
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
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

struct {
    integer nbitor, ibitor[500000]	/* was [5][100000] */;
} bitor_;

#define bitor_1 bitor_

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
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nring3, iring3[30000]	/* was [3][10000] */, nring4, iring4[
	    40000]	/* was [4][10000] */, nring5, iring5[50000]	/* 
	    was [5][10000] */, nring6, iring6[60000]	/* was [6][10000] */;
} ring_;

#define ring_1 ring_

struct {
    doublereal tors1[400000]	/* was [4][100000] */, tors2[400000]	/* 
	    was [4][100000] */, tors3[400000]	/* was [4][100000] */, tors4[
	    400000]	/* was [4][100000] */, tors5[400000]	/* was [4][
	    100000] */, tors6[400000]	/* was [4][100000] */;
    integer ntors, itors[400000]	/* was [4][100000] */;
} tors_;

#define tors_1 tors_

/* Table of constant values */

static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  subroutine rings  --  locate and store small rings  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     "rings" searches the structure for small rings and stores */
/*     their constituent atoms, and optionally reduces large rings */
/*     into their component smaller rings */

/*     note by default reducible rings are not removed since they */
/*     are needed for force field parameter assignment */


/* Subroutine */ int rings_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 RINGS  --  Too many 3-Membered Rings\002)"
	    ;
    static char fmt_30[] = "(/,\002 RINGS  --  Too many 4-Membered Rings\002)"
	    ;
    static char fmt_50[] = "(/,\002 RINGS  --  Too many 5-Membered Rings\002)"
	    ;
    static char fmt_70[] = "(/,\002 RINGS  --  Too many\002,\002 6-Membered "
	    "Rings\002)";
    static char fmt_90[] = "(/,\002 Three-Membered Rings Contained\002,\002 "
	    "in the Structure :\002,//,11x,\002Ring\002,14x,\002Atoms in Rin"
	    "g\002,/)";
    static char fmt_100[] = "(9x,i5,10x,3i6)";
    static char fmt_110[] = "(/,\002 Four-Membered Rings Contained\002,\002 "
	    "in the Structure :\002,//,11x,\002Ring\002,17x,\002Atoms in Rin"
	    "g\002,/)";
    static char fmt_120[] = "(9x,i5,10x,4i6)";
    static char fmt_130[] = "(/,\002 Five-Membered Rings Contained\002,\002 "
	    "in the Structure :\002,//,11x,\002Ring\002,20x,\002Atoms in Rin"
	    "g\002,/)";
    static char fmt_140[] = "(9x,i5,10x,5i6)";
    static char fmt_150[] = "(/,\002 Six-Membered Rings Contained\002,\002 i"
	    "n the Structure :\002,//,11x,\002Ring\002,23x,\002Atoms in Rin"
	    "g\002,/)";
    static char fmt_160[] = "(9x,i5,10x,6i6)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int torsions_(void);
    static integer i__, j, k, m, ia, ib, ic, id, ie, ig, imax, list[25000], 
	    list1, list2, list3, list4;
    extern /* Subroutine */ int fatal_(void), bonds_(void);
    static logical reduce;
    extern /* Subroutine */ int angles_(void), bitors_(void);

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_160, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define itors_ref(a_1,a_2) tors_1.itors[(a_2)*4 + a_1 - 5]
#define iring3_ref(a_1,a_2) ring_1.iring3[(a_2)*3 + a_1 - 4]
#define iring4_ref(a_1,a_2) ring_1.iring4[(a_2)*4 + a_1 - 5]
#define iring5_ref(a_1,a_2) ring_1.iring5[(a_2)*5 + a_1 - 6]
#define iring6_ref(a_1,a_2) ring_1.iring6[(a_2)*6 + a_1 - 7]
#define ibitor_ref(a_1,a_2) bitor_1.ibitor[(a_2)*5 + a_1 - 6]



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
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  bitor.i  --  bitorsions within the current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     nbitor  total number of bitorsions in the system */
/*     ibitor  numbers of the atoms in each bitorsion */




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
/*     ##  inform.i  --  control values for I/O and program flow  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     digits    decimal places output for energy and coordinates */
/*     iprint    steps between status printing (0=no printing) */
/*     iwrite    steps between coordinate dumps (0=no dumps) */
/*     isend     steps between socket communication (0=no sockets) */
/*     verbose   logical flag to turn on extra information */
/*     debug     logical flag to turn on full debug printing */
/*     holdup    logical flag to wait for carriage return on exit */
/*     abort     logical flag to stop execution at next chance */




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




/*     zero out the number of small rings in the structure */

    reduce = FALSE_;
    ring_1.nring3 = 0;
    ring_1.nring4 = 0;
    ring_1.nring5 = 0;
    ring_1.nring6 = 0;

/*     parse to find bonds, angles, torsions and bitorsions */

    if (bond_1.nbond == 0) {
	bonds_();
    }
    if (angle_1.nangle == 0) {
	angles_();
    }
    if (tors_1.ntors == 0) {
	torsions_();
    }
    if (bitor_1.nbitor == 0) {
	bitors_();
    }

/*     search for and store all of the 3-membered rings */

    i__1 = angle_1.nangle;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = iang_ref(1, i__);
	ib = iang_ref(2, i__);
	ic = iang_ref(3, i__);
	if (ib < ia && ib < ic) {
	    i__2 = couple_1.n12[ia - 1];
	    for (j = 1; j <= i__2; ++j) {
		if (i12_ref(j, ia) == ic) {
		    ++ring_1.nring3;
		    if (ring_1.nring3 > 10000) {
			io___7.ciunit = iounit_1.iout;
			s_wsfe(&io___7);
			e_wsfe();
			fatal_();
		    }
		    iring3_ref(1, ring_1.nring3) = ia;
		    iring3_ref(2, ring_1.nring3) = ib;
		    iring3_ref(3, ring_1.nring3) = ic;
		    goto L20;
		}
	    }
L20:
	    ;
	}
    }

/*     search for and store all of the 4-membered rings */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	list[i__ - 1] = 0;
    }
    i__1 = tors_1.ntors;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = itors_ref(1, i__);
	ib = itors_ref(2, i__);
	ic = itors_ref(3, i__);
	id = itors_ref(4, i__);
	if (ia < ic && id < ib) {
	    i__2 = couple_1.n12[ia - 1];
	    for (j = 1; j <= i__2; ++j) {
		if (i12_ref(j, ia) == id) {
		    ++ring_1.nring4;
		    if (ring_1.nring4 > 10000) {
			io___10.ciunit = iounit_1.iout;
			s_wsfe(&io___10);
			e_wsfe();
			fatal_();
		    }
		    iring4_ref(1, ring_1.nring4) = ia;
		    iring4_ref(2, ring_1.nring4) = ib;
		    iring4_ref(3, ring_1.nring4) = ic;
		    iring4_ref(4, ring_1.nring4) = id;

/*     remove the ring if it is reducible into smaller rings */

		    if (reduce) {
			list[ia - 1] = ring_1.nring4;
			list[ib - 1] = ring_1.nring4;
			list[ic - 1] = ring_1.nring4;
			list[id - 1] = ring_1.nring4;
			i__3 = ring_1.nring3;
			for (m = 1; m <= i__3; ++m) {
			    list1 = list[iring3_ref(1, m) - 1];
			    list2 = list[iring3_ref(2, m) - 1];
			    list3 = list[iring3_ref(3, m) - 1];
			    if (list1 == ring_1.nring4 && list2 == 
				    ring_1.nring4 && list3 == ring_1.nring4) {
				--ring_1.nring4;
				list[ia - 1] = 0;
				list[ib - 1] = 0;
				list[ic - 1] = 0;
				list[id - 1] = 0;
				goto L40;
			    }
			}
		    }
		    goto L40;
		}
	    }
L40:
	    ;
	}
    }

/*     search for and store all of the 5-membered rings */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	list[i__ - 1] = 0;
    }
    i__1 = bitor_1.nbitor;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ibitor_ref(1, i__);
	ib = ibitor_ref(2, i__);
	ic = ibitor_ref(3, i__);
	id = ibitor_ref(4, i__);
	ie = ibitor_ref(5, i__);
	if (ia < id && ie < ib && min(ia,ie) < ic) {
	    i__2 = couple_1.n12[ia - 1];
	    for (j = 1; j <= i__2; ++j) {
		if (i12_ref(j, ia) == ie) {
		    ++ring_1.nring5;
		    if (ring_1.nring5 > 10000) {
			io___16.ciunit = iounit_1.iout;
			s_wsfe(&io___16);
			e_wsfe();
			fatal_();
		    }
		    iring5_ref(1, ring_1.nring5) = ia;
		    iring5_ref(2, ring_1.nring5) = ib;
		    iring5_ref(3, ring_1.nring5) = ic;
		    iring5_ref(4, ring_1.nring5) = id;
		    iring5_ref(5, ring_1.nring5) = ie;

/*     remove the ring if it is reducible into smaller rings */

		    if (reduce) {
			list[ia - 1] = ring_1.nring5;
			list[ib - 1] = ring_1.nring5;
			list[ic - 1] = ring_1.nring5;
			list[id - 1] = ring_1.nring5;
			list[ie - 1] = ring_1.nring5;
			i__3 = ring_1.nring3;
			for (m = 1; m <= i__3; ++m) {
			    list1 = list[iring3_ref(1, m) - 1];
			    list2 = list[iring3_ref(2, m) - 1];
			    list3 = list[iring3_ref(3, m) - 1];
			    if (list1 == ring_1.nring5 && list2 == 
				    ring_1.nring5 && list3 == ring_1.nring5) {
				--ring_1.nring5;
				list[ia - 1] = 0;
				list[ib - 1] = 0;
				list[ic - 1] = 0;
				list[id - 1] = 0;
				list[ie - 1] = 0;
				goto L60;
			    }
			}
		    }
		    goto L60;
		}
	    }
L60:
	    ;
	}
    }

/*     search for and store all of the 6-membered rings */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	list[i__ - 1] = 0;
    }
    i__1 = bitor_1.nbitor;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ibitor_ref(1, i__);
	ib = ibitor_ref(2, i__);
	ic = ibitor_ref(3, i__);
	id = ibitor_ref(4, i__);
	ie = ibitor_ref(5, i__);
/* Computing MAX */
	i__2 = max(ia,ib), i__2 = max(i__2,ic), i__2 = max(i__2,id);
	imax = max(i__2,ie);
	i__2 = couple_1.n12[ia - 1];
	for (j = 1; j <= i__2; ++j) {
	    ig = i12_ref(j, ia);
	    if (ig > imax) {
		i__3 = couple_1.n12[ie - 1];
		for (k = 1; k <= i__3; ++k) {
		    if (i12_ref(k, ie) == ig) {
			++ring_1.nring6;
			if (ring_1.nring6 > 10000) {
			    io___20.ciunit = iounit_1.iout;
			    s_wsfe(&io___20);
			    e_wsfe();
			    fatal_();
			}
			iring6_ref(1, ring_1.nring6) = ia;
			iring6_ref(2, ring_1.nring6) = ib;
			iring6_ref(3, ring_1.nring6) = ic;
			iring6_ref(4, ring_1.nring6) = id;
			iring6_ref(5, ring_1.nring6) = ie;
			iring6_ref(6, ring_1.nring6) = ig;

/*     remove the ring if it is reducible into smaller rings */

			if (reduce) {
			    list[ia - 1] = ring_1.nring6;
			    list[ib - 1] = ring_1.nring6;
			    list[ic - 1] = ring_1.nring6;
			    list[id - 1] = ring_1.nring6;
			    list[ie - 1] = ring_1.nring6;
			    list[ig - 1] = ring_1.nring6;
			    i__4 = ring_1.nring3;
			    for (m = 1; m <= i__4; ++m) {
				list1 = list[iring3_ref(1, m) - 1];
				list2 = list[iring3_ref(2, m) - 1];
				list3 = list[iring3_ref(3, m) - 1];
				if (list1 == ring_1.nring6 && list2 == 
					ring_1.nring6 && list3 == 
					ring_1.nring6) {
				    --ring_1.nring6;
				    list[ia - 1] = 0;
				    list[ib - 1] = 0;
				    list[ic - 1] = 0;
				    list[id - 1] = 0;
				    list[ie - 1] = 0;
				    list[ig - 1] = 0;
				    goto L80;
				}
			    }
			    i__4 = ring_1.nring4;
			    for (m = 1; m <= i__4; ++m) {
				list1 = list[iring4_ref(1, m) - 1];
				list2 = list[iring4_ref(2, m) - 1];
				list3 = list[iring4_ref(3, m) - 1];
				list4 = list[iring4_ref(4, m) - 1];
				if (list1 == ring_1.nring6 && list2 == 
					ring_1.nring6 && list3 == 
					ring_1.nring6 && list4 == 
					ring_1.nring6) {
				    --ring_1.nring6;
				    list[ia - 1] = 0;
				    list[ib - 1] = 0;
				    list[ic - 1] = 0;
				    list[id - 1] = 0;
				    list[ie - 1] = 0;
				    list[ig - 1] = 0;
				    goto L80;
				}
			    }
			}
L80:
			;
		    }
		}
	    }
	}
    }

/*     print out lists of the small rings in the structure */

    if (inform_1.debug) {
	if (ring_1.nring3 > 0) {
	    io___22.ciunit = iounit_1.iout;
	    s_wsfe(&io___22);
	    e_wsfe();
	    i__1 = ring_1.nring3;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___23.ciunit = iounit_1.iout;
		s_wsfe(&io___23);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&iring3_ref(j, i__), (ftnlen)sizeof(
			    integer));
		}
		e_wsfe();
	    }
	}
	if (ring_1.nring4 > 0) {
	    io___24.ciunit = iounit_1.iout;
	    s_wsfe(&io___24);
	    e_wsfe();
	    i__1 = ring_1.nring4;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___25.ciunit = iounit_1.iout;
		s_wsfe(&io___25);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 4; ++j) {
		    do_fio(&c__1, (char *)&iring4_ref(j, i__), (ftnlen)sizeof(
			    integer));
		}
		e_wsfe();
	    }
	}
	if (ring_1.nring5 > 0) {
	    io___26.ciunit = iounit_1.iout;
	    s_wsfe(&io___26);
	    e_wsfe();
	    i__1 = ring_1.nring5;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___27.ciunit = iounit_1.iout;
		s_wsfe(&io___27);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 5; ++j) {
		    do_fio(&c__1, (char *)&iring5_ref(j, i__), (ftnlen)sizeof(
			    integer));
		}
		e_wsfe();
	    }
	}
	if (ring_1.nring6 > 0) {
	    io___28.ciunit = iounit_1.iout;
	    s_wsfe(&io___28);
	    e_wsfe();
	    i__1 = ring_1.nring6;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___29.ciunit = iounit_1.iout;
		s_wsfe(&io___29);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 6; ++j) {
		    do_fio(&c__1, (char *)&iring6_ref(j, i__), (ftnlen)sizeof(
			    integer));
		}
		e_wsfe();
	    }
	}
    }
    return 0;
} /* rings_ */

#undef ibitor_ref
#undef iring6_ref
#undef iring5_ref
#undef iring4_ref
#undef iring3_ref
#undef itors_ref
#undef iang_ref
#undef i12_ref


