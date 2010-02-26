/* initres.f -- translated by f2c (version 20050501).
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
    char amino[93], nuclz[36], amino1[31], nuclz1[12];
} resdue_;

#define resdue_1 resdue_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine initres  --  setup biopolymer residue names  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "initres" sets names for biopolymer residue types used in */
/*     PDB file conversion and automated generation of structures */


/* Subroutine */ int initres_(void)
{
    /* Initialized data */

    static char acid1[1*31] = "G" "A" "V" "L" "I" "S" "T" "C" "C" "P" "F" 
	    "Y" "W" "H" "U" "Z" "D" "N" "E" "Q" "M" "K" "R" "O" "B" "J" "f" 
	    "a" "n" "m" "X";
    static char acid3[3*31] = "GLY" "ALA" "VAL" "LEU" "ILE" "SER" "THR" "CYS" 
	    "CYX" "PRO" "PHE" "TYR" "TRP" "HIS" "HID" "HIE" "ASP" "ASN" "GLU" 
	    "GLN" "MET" "LYS" "ARG" "ORN" "AIB" "PCA" "FOR" "ACE" "NH2" "NME" 
	    "UNK";
    static char base1[1*12] = "A" "G" "C" "U" "D" "B" "I" "T" "1" "2" "3" 
	    "X";
    static char base3[3*12] = "A  " "G  " "C  " "U  " "DA " "DG " "DC " "DT " 
	    "MP " "DP " "TP " "UNK";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;


#define acid3_ref(a_0,a_1) &acid3[(a_1)*3 + a_0 - 3]
#define base3_ref(a_0,a_1) &base3[(a_1)*3 + a_0 - 3]
#define amino_ref(a_0,a_1) &resdue_1.amino[(a_1)*3 + a_0 - 3]
#define nuclz_ref(a_0,a_1) &resdue_1.nuclz[(a_1)*3 + a_0 - 3]



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
/*     ##  resdue.i  --  standard biopolymer residue abbreviations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     amino    three-letter abbreviations for amino acids types */
/*     nuclz    three-letter abbreviations for nucleic acids types */
/*     amino1   one-letter abbreviations for amino acids types */
/*     nuclz1   one-letter abbreviations for nucleic acids types */




/*     set values for the 1- and 3-letter amino acid names */

    for (i__ = 1; i__ <= 31; ++i__) {
	s_copy(amino_ref(0, i__), acid3_ref(0, i__), (ftnlen)3, (ftnlen)3);
	*(unsigned char *)&resdue_1.amino1[i__ - 1] = *(unsigned char *)&
		acid1[i__ - 1];
    }

/*     set values for the 1- and 3-letter nucleic acid names */

    for (i__ = 1; i__ <= 12; ++i__) {
	s_copy(nuclz_ref(0, i__), base3_ref(0, i__), (ftnlen)3, (ftnlen)3);
	*(unsigned char *)&resdue_1.nuclz1[i__ - 1] = *(unsigned char *)&
		base1[i__ - 1];
    }
    return 0;
} /* initres_ */

#undef nuclz_ref
#undef amino_ref
#undef base3_ref
#undef acid3_ref


