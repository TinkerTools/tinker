/* makeref.f -- translated by f2c (version 20050501).
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
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

struct {
    doublereal xref[250000]	/* was [25000][10] */, yref[250000]	/* 
	    was [25000][10] */, zref[250000]	/* was [25000][10] */;
    integer nref[10], reftyp[250000]	/* was [25000][10] */, n12ref[250000]	
	    /* was [25000][10] */, i12ref[2000000]	/* was [8][25000][10] 
	    */, refleng[10], refltitle[10];
    char refnam[750000]	/* was [25000][10] */, reffile[1200], reftitle[1200];
} refer_;

#define refer_1 refer_

struct {
    integer ltitle;
    char title[120];
} titles_;

#define titles_1 titles_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine makeref  --  copy structure to reference area  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "makeref" copies the information contained in the "xyz" file */
/*     of the current structure into corresponding reference areas */


/* Subroutine */ int makeref_(integer *iref)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j;


#define reftitle_ref(a_0,a_1) &refer_1.reftitle[(a_1)*120 + a_0 - 120]
#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define xref_ref(a_1,a_2) refer_1.xref[(a_2)*25000 + a_1 - 25001]
#define yref_ref(a_1,a_2) refer_1.yref[(a_2)*25000 + a_1 - 25001]
#define zref_ref(a_1,a_2) refer_1.zref[(a_2)*25000 + a_1 - 25001]
#define i12ref_ref(a_1,a_2,a_3) refer_1.i12ref[((a_3)*25000 + (a_2))*8 + \
a_1 - 200009]
#define n12ref_ref(a_1,a_2) refer_1.n12ref[(a_2)*25000 + a_1 - 25001]
#define refnam_ref(a_0,a_1,a_2) &refer_1.refnam[((a_2)*25000 + a_1)*3 + \
a_0 - 75003]
#define reftyp_ref(a_1,a_2) refer_1.reftyp[(a_2)*25000 + a_1 - 25001]
#define reffile_ref(a_0,a_1) &refer_1.reffile[(a_1)*120 + a_0 - 120]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  files.i  --  name and number of current structure files  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     nprior     number of previously existing cycle files */
/*     ldir       length in characters of the directory name */
/*     leng       length in characters of the base filename */
/*     filename   base filename used by default for all files */
/*     outfile    output filename used for intermediate results */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  refer.i  --  storage of reference atomic coordinate set  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xref        reference x-coordinates for atoms in each system */
/*     yref        reference y-coordinates for atoms in each system */
/*     zref        reference z-coordinates for atoms in each system */
/*     nref        total number of atoms in each reference system */
/*     reftyp      atom types of the atoms in each reference system */
/*     n12ref      number of atoms bonded to each reference atom */
/*     i12ref      atom numbers of atoms 1-2 connected to each atom */
/*     refleng     length in characters of each reference filename */
/*     refltitle   length in characters of each reference title line */
/*     refnam      atom names of the atoms in each reference system */
/*     reffile     base filename for each reference system */
/*     reftitle    title used to describe each reference system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  titles.i  --  title for the current molecular system  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     ltitle   length in characters of the nonblank title string */
/*     title    title used to describe the current structure */




/*     copy the filename and title line for the structure */

    s_copy(reffile_ref(0, *iref), files_1.filename, (ftnlen)120, (ftnlen)120);
    refer_1.refleng[*iref - 1] = files_1.leng;
    s_copy(reftitle_ref(0, *iref), titles_1.title, (ftnlen)120, (ftnlen)120);
    refer_1.refltitle[*iref - 1] = titles_1.ltitle;

/*     copy the coordinates, type and connectivity of each atom */

    refer_1.nref[*iref - 1] = atoms_1.n;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(refnam_ref(0, i__, *iref), name___ref(0, i__), (ftnlen)3, (
		ftnlen)3);
	xref_ref(i__, *iref) = atoms_1.x[i__ - 1];
	yref_ref(i__, *iref) = atoms_1.y[i__ - 1];
	zref_ref(i__, *iref) = atoms_1.z__[i__ - 1];
	reftyp_ref(i__, *iref) = atoms_1.type__[i__ - 1];
	n12ref_ref(i__, *iref) = couple_1.n12[i__ - 1];
	i__2 = couple_1.n12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    i12ref_ref(j, i__, *iref) = i12_ref(j, i__);
	}
    }
    return 0;
} /* makeref_ */

#undef reffile_ref
#undef reftyp_ref
#undef refnam_ref
#undef n12ref_ref
#undef i12ref_ref
#undef zref_ref
#undef yref_ref
#undef xref_ref
#undef name___ref
#undef i12_ref
#undef reftitle_ref


