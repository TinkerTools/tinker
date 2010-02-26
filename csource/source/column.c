/* column.f -- translated by f2c (version 20050501).
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



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine column  --  access Hessian elements by column  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "column" takes the off-diagonal Hessian elements stored */
/*     as sparse rows and sets up indices to allow column access */


/* Subroutine */ int column_(integer *nvar, integer *hinit, integer *hstop, 
	integer *hindex, integer *cinit, integer *cstop, integer *cindex, 
	integer *cvalue)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, m;



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




/*     zero out the start and end marker for each column */

    /* Parameter adjustments */
    --cvalue;
    --cindex;
    --cstop;
    --cinit;
    --hindex;
    --hstop;
    --hinit;

    /* Function Body */
    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cinit[i__] = 0;
	cstop[i__] = 0;
    }

/*     count the number of elements in each column */

    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = hstop[i__];
	for (j = hinit[i__]; j <= i__2; ++j) {
	    k = hindex[j];
	    ++cstop[k];
	}
    }

/*     set each start marker just past last element for its column */

    cinit[1] = cstop[1] + 1;
    i__1 = *nvar;
    for (i__ = 2; i__ <= i__1; ++i__) {
	cinit[i__] = cinit[i__ - 1] + cstop[i__];
    }

/*     set column index by scanning rows in reverse order */

    for (i__ = *nvar; i__ >= 1; --i__) {
	i__1 = hstop[i__];
	for (j = hinit[i__]; j <= i__1; ++j) {
	    k = hindex[j];
	    m = cinit[k] - 1;
	    cinit[k] = m;
	    cindex[m] = i__;
	    cvalue[m] = j;
	}
    }

/*     convert from number of elements to end marker for column */

    i__1 = *nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cstop[i__] = cinit[i__] + cstop[i__] - 1;
    }
    return 0;
} /* column_ */

