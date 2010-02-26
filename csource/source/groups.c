/* groups.f -- translated by f2c (version 20050501).
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
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine groups  --  group membership of set of atoms  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "groups" tests a set of atoms to see if all are members */
/*     of a single atom group or a pair of atom groups; if so, */
/*     then the correct intra- or intergroup weight is assigned */


/* Subroutine */ int groups_(logical *proceed, doublereal *weigh, integer *ia,
	 integer *ib, integer *ic, integer *id, integer *ie, integer *ig)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer iga, igb, igc, igd, ige, igg, gmin, gmax, nset;


#define wgrp_ref(a_1,a_2) group_1.wgrp[(a_2)*1001 + a_1 - 0]



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




/*     determine the number of atoms in the set to be compared */

    nset = 0;
    *weigh = 0.;
    if (*ig != 0) {
	nset = 6;
    } else if (*ie != 0) {
	nset = 5;
    } else if (*id != 0) {
	nset = 4;
    } else if (*ic != 0) {
	nset = 3;
    } else if (*ib != 0) {
	nset = 2;
    } else if (*ia != 0) {
	nset = 1;
    }

/*     check group membership for a set containing one atom */

    if (nset == 1) {
	iga = group_1.grplist[*ia - 1];
	*weigh = wgrp_ref(iga, iga);

/*     check group membership for a set containing two atoms */

    } else if (nset == 2) {
	iga = group_1.grplist[*ia - 1];
	igb = group_1.grplist[*ib - 1];
	*weigh = wgrp_ref(iga, igb);

/*     check group membership for a set containing three atoms */

    } else if (nset == 3) {
	iga = group_1.grplist[*ia - 1];
	igb = group_1.grplist[*ib - 1];
	igc = group_1.grplist[*ic - 1];
	if (iga == igb || igb == igc) {
	    *weigh = wgrp_ref(iga, igc);
	} else if (iga == igc) {
	    *weigh = wgrp_ref(iga, igb);
	}

/*     check group membership for a set containing four atoms */

    } else if (nset == 4) {
	iga = group_1.grplist[*ia - 1];
	igb = group_1.grplist[*ib - 1];
	igc = group_1.grplist[*ic - 1];
	igd = group_1.grplist[*id - 1];
/* Computing MIN */
	i__1 = min(iga,igb), i__1 = min(i__1,igc);
	gmin = min(i__1,igd);
/* Computing MAX */
	i__1 = max(iga,igb), i__1 = max(i__1,igc);
	gmax = max(i__1,igd);
	if ((iga == gmin || iga == gmax) && (igb == gmin || igb == gmax) && (
		igc == gmin || igc == gmax) && (igd == gmin || igd == gmax)) {
	    *weigh = wgrp_ref(gmin, gmax);
	}

/*     check group membership for a set containing five atoms */

    } else if (nset == 5) {
	iga = group_1.grplist[*ia - 1];
	igb = group_1.grplist[*ib - 1];
	igc = group_1.grplist[*ic - 1];
	igd = group_1.grplist[*id - 1];
	ige = group_1.grplist[*ie - 1];
/* Computing MIN */
	i__1 = min(iga,igb), i__1 = min(i__1,igc), i__1 = min(i__1,igd);
	gmin = min(i__1,ige);
/* Computing MAX */
	i__1 = max(iga,igb), i__1 = max(i__1,igc), i__1 = max(i__1,igd);
	gmax = max(i__1,ige);
	if ((iga == gmin || iga == gmax) && (igb == gmin || igb == gmax) && (
		igc == gmin || igc == gmax) && (igd == gmin || igd == gmax) &&
		 (ige == gmin || ige == gmax)) {
	    *weigh = wgrp_ref(gmin, gmax);
	}

/*     check group membership for a set containing five atoms */

    } else if (nset == 6) {
	iga = group_1.grplist[*ia - 1];
	igb = group_1.grplist[*ib - 1];
	igc = group_1.grplist[*ic - 1];
	igd = group_1.grplist[*id - 1];
	ige = group_1.grplist[*ie - 1];
	igg = group_1.grplist[*ig - 1];
/* Computing MIN */
	i__1 = min(iga,igb), i__1 = min(i__1,igc), i__1 = min(i__1,igd), i__1 
		= min(i__1,ige);
	gmin = min(i__1,igg);
/* Computing MAX */
	i__1 = max(iga,igb), i__1 = max(i__1,igc), i__1 = max(i__1,igd), i__1 
		= max(i__1,ige);
	gmax = max(i__1,igg);
	if ((iga == gmin || iga == gmax) && (igb == gmin || igb == gmax) && (
		igc == gmin || igc == gmax) && (igd == gmin || igd == gmax) &&
		 (ige == gmin || ige == gmax) && (igg == gmin || igg == gmax))
		 {
	    *weigh = wgrp_ref(gmin, gmax);
	}
    }

/*     interaction will be used if its group has nonzero weight */

    if (*weigh == 0.) {
	*proceed = FALSE_;
    } else {
	*proceed = TRUE_;
    }
    return 0;
} /* groups_ */

#undef wgrp_ref


