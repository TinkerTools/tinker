/* extra3.f -- translated by f2c (version 20050501).
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
    integer neb, nea, neba, neub, neaa, neopb, neopd, neid, neit, net, nept, 
	    nebt, nett, nev, nec, necd, ned, nem, nep, new__, ner, nes, nelf, 
	    neg, nex;
} action_;

#define action_1 action_

struct {
    doublereal aesum[25000], aeb[25000], aea[25000], aeba[25000], aeub[25000],
	     aeaa[25000], aeopb[25000], aeopd[25000], aeid[25000], aeit[25000]
	    , aet[25000], aept[25000], aebt[25000], aett[25000], aev[25000], 
	    aec[25000], aecd[25000], aed[25000], aem[25000], aep[25000], aer[
	    25000], aes[25000], aelf[25000], aeg[25000], aex[25000];
} analyz_;

#define analyz_1 analyz_

struct {
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

struct {
    doublereal esum, eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, 
	    ebt, ett, ev, ec, ecd, ed, em, ep, er, es, elf, eg, ex;
} energi_;

#define energi_1 energi_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine extra3  --  user defined extra potentials  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "extra3" calculates any additional user defined potential */
/*     contribution and also partitions the energy among the atoms */


/* Subroutine */ int extra3_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;



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
/*     ##  action.i  --  total number of each energy term computed  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     neb     number of bond stretch energy terms computed */
/*     nea     number of angle bend energy terms computed */
/*     neba    number of stretch-bend energy terms computed */
/*     neub    number of Urey-Bradley energy terms computed */
/*     neaa    number of angle-angle energy terms computed */
/*     neopb   number of out-of-plane bend energy terms computed */
/*     neopd   number of out-of-plane distance energy terms computed */
/*     neid    number of improper dihedral energy terms computed */
/*     neit    number of improper torsion energy terms computed */
/*     net     number of torsional energy terms computed */
/*     nept    number of pi-orbital torsion energy terms computed */
/*     nebt    number of stretch-torsion energy terms computed */
/*     nett    number of torsion-torsion energy terms computed */
/*     nev     number of van der Waals energy terms computed */
/*     nec     number of charge-charge energy terms computed */
/*     necd    number of charge-dipole energy terms computed */
/*     ned     number of dipole-dipole energy terms computed */
/*     nem     number of multipole energy terms computed */
/*     nep     number of polarization energy terms computed */
/*     new     number of Ewald summation energy terms computed */
/*     ner     number of reaction field energy terms computed */
/*     nes     number of solvation energy terms computed */
/*     nelf    number of metal ligand field energy terms computed */
/*     neg     number of geometric restraint energy terms computed */
/*     nex     number of extra energy terms computed */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  analyz.i  --  energy components partitioned over atoms  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     aesum   total potential energy partitioned over atoms */
/*     aeb     bond stretch energy partitioned over atoms */
/*     aea     angle bend energy partitioned over atoms */
/*     aeba    stretch-bend energy partitioned over atoms */
/*     aeub    Urey-Bradley energy partitioned over atoms */
/*     aeaa    angle-angle energy partitioned over atoms */
/*     aeopb   out-of-plane bend energy partitioned over atoms */
/*     aeopd   out-of-plane distance energy partitioned over atoms */
/*     aeid    improper dihedral energy partitioned over atoms */
/*     aeit    improper torsion energy partitioned over atoms */
/*     aet     torsional energy partitioned over atoms */
/*     aept    pi-orbital torsion energy partitioned over atoms */
/*     aebt    stretch-torsion energy partitioned over atoms */
/*     aett    torsion-torsion energy partitioned over atoms */
/*     aev     van der Waals energy partitioned over atoms */
/*     aec     charge-charge energy partitioned over atoms */
/*     aecd    charge-dipole energy partitioned over atoms */
/*     aed     dipole-dipole energy partitioned over atoms */
/*     aem     multipole energy partitioned over atoms */
/*     aep     polarization energy partitioned over atoms */
/*     aer     reaction field energy partitioned over atoms */
/*     aes     solvation energy partitioned over atoms */
/*     aelf    metal ligand field energy partitioned over atoms */
/*     aeg     geometric restraint energy partitioned over atoms */
/*     aex     extra energy term partitioned over atoms */




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




/*     zero out the energy due to extra potential terms */

    action_1.nex = 0;
    energi_1.ex = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	analyz_1.aex[i__ - 1] = 0.;
    }

/*     add any user-defined extra potentials and partitioning */

/*     e = ...... */
/*     nex = nex + 1 */
/*     ex = ex + e */
/*     aex(i) = ...... */

    return 0;
} /* extra3_ */

