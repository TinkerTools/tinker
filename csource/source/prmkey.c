/* prmkey.f -- translated by f2c (version 20050501).
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
    doublereal angunit, stbnunit, aaunit, opbunit, opdunit, cang, qang, pang, 
	    sang, copb, qopb, popb, sopb, copd, qopd, popd, sopd;
    char angtyp[600000], opbtyp[8];
} angpot_;

#define angpot_1 angpot_

struct {
    doublereal cbnd, qbnd, bndunit;
    char bndtyp[8];
} bndpot_;

#define bndpot_1 bndpot_

struct {
    doublereal electric, dielec, ebuffer, c2scale, c3scale, c4scale, c5scale;
    logical neutnbr, neutcut;
} chgpot_;

#define chgpot_1 chgpot_

struct {
    integer biotyp[10000];
    char forcefield[20];
} fields_;

#define fields_1 fields_

struct {
    doublereal m2scale, m3scale, m4scale, m5scale;
} mplpot_;

#define mplpot_1 mplpot_

struct {
    doublereal poleps, polsor, p2scale, p3scale, p4scale, p5scale, d1scale, 
	    d2scale, d3scale, d4scale, u1scale, u2scale, u3scale, u4scale;
    char poltyp[6];
} polpot_;

#define polpot_1 polpot_

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
    doublereal rfsize, rfbulkd;
    integer rfterms;
} rxnpot_;

#define rxnpot_1 rxnpot_

struct {
    doublereal idihunit, itorunit, torsunit, ptorunit, storunit, ttorunit;
} torpot_;

#define torpot_1 torpot_

struct {
    doublereal cury, qury, ureyunit;
} urypot_;

#define urypot_1 urypot_

struct {
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__3 = 3;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine prmkey  --  interpret force field keywords  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "prmkey" parses a text string to extract keywords related to */
/*     force field potential energy functional forms and constants */


/* Subroutine */ int prmkey_(char *text, ftnlen text_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);

    /* Local variables */
    static integer next;
    static char value[4], record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int potoff_(void), getword_(char *, char *, 
	    integer *, ftnlen, ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___6 = { 1, string, 1, 0, 120, 1 };
    static icilist io___7 = { 1, string, 1, 0, 120, 1 };
    static icilist io___8 = { 1, string, 1, 0, 120, 1 };
    static icilist io___9 = { 1, string, 1, 0, 120, 1 };
    static icilist io___10 = { 1, string, 1, 0, 120, 1 };
    static icilist io___11 = { 1, string, 1, 0, 120, 1 };
    static icilist io___12 = { 1, string, 1, 0, 120, 1 };
    static icilist io___13 = { 1, string, 1, 0, 120, 1 };
    static icilist io___14 = { 1, string, 1, 0, 120, 1 };
    static icilist io___15 = { 1, string, 1, 0, 120, 1 };
    static icilist io___16 = { 1, string, 1, 0, 120, 1 };
    static icilist io___17 = { 1, string, 1, 0, 120, 1 };
    static icilist io___18 = { 1, string, 1, 0, 120, 1 };
    static icilist io___19 = { 1, string, 1, 0, 120, 1 };
    static icilist io___20 = { 1, string, 1, 0, 120, 1 };
    static icilist io___21 = { 1, string, 1, 0, 120, 1 };
    static icilist io___22 = { 1, string, 1, 0, 120, 1 };
    static icilist io___23 = { 1, string, 1, 0, 120, 1 };
    static icilist io___24 = { 1, string, 1, 0, 120, 1 };
    static icilist io___25 = { 1, string, 1, 0, 120, 1 };
    static icilist io___26 = { 1, string, 1, 0, 120, 1 };
    static icilist io___27 = { 1, string, 1, 0, 120, 1 };
    static icilist io___28 = { 1, string, 1, 0, 120, 1 };
    static icilist io___29 = { 1, string, 1, 0, 120, 1 };
    static icilist io___30 = { 1, string, 1, 0, 120, 1 };
    static icilist io___31 = { 1, string, 1, 0, 120, 1 };
    static icilist io___32 = { 1, string, 1, 0, 120, 1 };
    static icilist io___33 = { 1, string, 1, 0, 120, 1 };
    static icilist io___34 = { 1, string, 1, 0, 120, 1 };
    static icilist io___35 = { 1, string, 1, 0, 120, 1 };
    static icilist io___36 = { 1, string, 1, 0, 120, 1 };
    static icilist io___37 = { 1, string, 1, 0, 120, 1 };
    static icilist io___38 = { 1, string, 1, 0, 120, 1 };
    static icilist io___39 = { 1, string, 1, 0, 120, 1 };
    static icilist io___40 = { 1, string, 1, 0, 120, 1 };
    static icilist io___41 = { 1, string, 1, 0, 120, 1 };
    static icilist io___42 = { 1, string, 1, 0, 120, 1 };
    static icilist io___43 = { 1, string, 1, 0, 120, 1 };
    static icilist io___44 = { 1, string, 1, 0, 120, 1 };
    static icilist io___45 = { 1, string, 1, 0, 120, 1 };
    static icilist io___46 = { 1, string, 1, 0, 120, 1 };
    static icilist io___47 = { 1, string, 1, 0, 120, 1 };
    static icilist io___48 = { 1, string, 1, 0, 120, 1 };
    static icilist io___49 = { 1, string, 1, 0, 120, 1 };
    static icilist io___50 = { 1, string, 1, 0, 120, 1 };
    static icilist io___51 = { 1, string, 1, 0, 120, 1 };
    static icilist io___52 = { 1, string, 1, 0, 120, 1 };
    static icilist io___53 = { 1, string, 1, 0, 120, 1 };
    static icilist io___54 = { 1, string, 1, 0, 120, 1 };
    static icilist io___55 = { 1, string, 1, 0, 120, 1 };
    static icilist io___56 = { 1, string, 1, 0, 120, 1 };
    static icilist io___57 = { 1, string, 1, 0, 120, 1 };
    static icilist io___58 = { 1, string, 1, 0, 120, 1 };
    static icilist io___59 = { 1, string, 1, 0, 120, 1 };
    static icilist io___60 = { 1, string, 1, 0, 120, 1 };
    static icilist io___61 = { 1, string, 1, 0, 120, 1 };
    static icilist io___62 = { 1, string, 1, 0, 120, 1 };
    static icilist io___63 = { 1, string, 1, 0, 120, 1 };
    static icilist io___64 = { 1, string, 1, 0, 120, 1 };
    static icilist io___65 = { 1, string, 1, 0, 120, 1 };
    static icilist io___66 = { 1, string, 1, 0, 120, 1 };
    static icilist io___67 = { 1, string, 1, 0, 120, 1 };
    static icilist io___68 = { 1, string, 1, 0, 120, 1 };
    static icilist io___69 = { 1, string, 1, 0, 120, 1 };




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
/*     ##  angpot.i  --  specifics of angle bend functional forms  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     angunit    convert angle bending energy to kcal/mole */
/*     stbnunit   convert stretch-bend energy to kcal/mole */
/*     aaunit     convert angle-angle energy to kcal/mole */
/*     opbunit    convert out-of-plane bend energy to kcal/mole */
/*     opdunit    convert out-of-plane distance energy to kcal/mole */
/*     cang       cubic coefficient in angle bending potential */
/*     qang       quartic coefficient in angle bending potential */
/*     pang       quintic coefficient in angle bending potential */
/*     sang       sextic coefficient in angle bending potential */
/*     copb       cubic coefficient in out-of-plane bend potential */
/*     qopb       quartic coefficient in out-of-plane bend potential */
/*     popb       quintic coefficient in out-of-plane bend potential */
/*     sopb       sextic coefficient in out-of-plane bend potential */
/*     copd       cubic coefficient in out-of-plane distance potential */
/*     qopd       quartic coefficient in out-of-plane distance potential */
/*     popd       quintic coefficient in out-of-plane distance potential */
/*     sopd       sextic coefficient in out-of-plane distance potential */
/*     angtyp     type of angle bending function for each bond angle */
/*     opbtyp     type of out-of-plane bend potential energy function */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  bndpot.i  --  specifics of bond stretch functional forms  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     cbnd      cubic coefficient in bond stretch potential */
/*     qbnd      quartic coefficient in bond stretch potential */
/*     bndunit   convert bond stretch energy to kcal/mole */
/*     bndtyp    type of bond stretch potential energy function */




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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  fields.i  --  molecular mechanics force field description  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     biotyp       force field atom type of each biopolymer type */
/*     forcefield   string used to describe the current forcefield */




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
/*     ##  rxnpot.i  --  specifics of reaction field functional form  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     rfsize    radius of reaction field sphere centered at origin */
/*     rfbulkd   bulk dielectric constant of reaction field continuum */
/*     rfterms   number of terms to use in reaction field summation */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  torpot.i  --  specifics of torsional functional forms  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     idihunit  convert improper dihedral energy to kcal/mole */
/*     itorunit  convert improper torsion amplitudes to kcal/mole */
/*     torsunit  convert torsional parameter amplitudes to kcal/mole */
/*     ptorunit  convert pi-orbital torsion energy to kcal/mole */
/*     storunit  convert stretch-torsion energy to kcal/mole */
/*     ttorunit  convert stretch-torsion energy to kcal/mole */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  urypot.i  --  specifics of Urey-Bradley functional form  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     cury       cubic coefficient in Urey-Bradley potential */
/*     qury       quartic coefficient in Urey-Bradley potential */
/*     ureyunit   convert Urey-Bradley energy to kcal/mole */




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




/*     parse the line to extract any possible keyword */

    s_copy(record, text, (ftnlen)120, (ftnlen)120);
    next = 1;
    upcase_(record, (ftnlen)120);
    gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));

/*     select the individual force field potential terms */

    if (s_cmp(keyword, "BONDTERM ", (ftnlen)9, (ftnlen)9) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_bond__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_bond__ = FALSE_;
	}
    } else if (s_cmp(keyword, "ANGLETERM ", (ftnlen)10, (ftnlen)10) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_angle__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_angle__ = FALSE_;
	}
    } else if (s_cmp(keyword, "STRBNDTERM ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_strbnd__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_strbnd__ = FALSE_;
	}
    } else if (s_cmp(keyword, "UREYTERM ", (ftnlen)9, (ftnlen)9) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_urey__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_urey__ = FALSE_;
	}
    } else if (s_cmp(keyword, "ANGANGTERM ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_angang__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_angang__ = FALSE_;
	}
    } else if (s_cmp(keyword, "OPBENDTERM ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_opbend__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_opbend__ = FALSE_;
	}
    } else if (s_cmp(keyword, "OPDISTTERM ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_opdist__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_opdist__ = FALSE_;
	}
    } else if (s_cmp(keyword, "IMPROPTERM ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_improp__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_improp__ = FALSE_;
	}
    } else if (s_cmp(keyword, "IMPTORSTERM ", (ftnlen)12, (ftnlen)12) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_imptor__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_imptor__ = FALSE_;
	}
    } else if (s_cmp(keyword, "TORSIONTERM ", (ftnlen)12, (ftnlen)12) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_tors__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_tors__ = FALSE_;
	}
    } else if (s_cmp(keyword, "PITORSTERM ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_pitors__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_pitors__ = FALSE_;
	}
    } else if (s_cmp(keyword, "STRTORTERM ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_strtor__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_strtor__ = FALSE_;
	}
    } else if (s_cmp(keyword, "TORTORTERM ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_tortor__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_tortor__ = FALSE_;
	}
    } else if (s_cmp(keyword, "VDWTERM ", (ftnlen)8, (ftnlen)8) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_vdw__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_vdw__ = FALSE_;
	}
    } else if (s_cmp(keyword, "CHARGETERM ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_charge__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_charge__ = FALSE_;
	}
    } else if (s_cmp(keyword, "CHGDPLTERM ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_chgdpl__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_chgdpl__ = FALSE_;
	}
    } else if (s_cmp(keyword, "DIPOLETERM ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_dipole__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_dipole__ = FALSE_;
	}
    } else if (s_cmp(keyword, "MPOLETERM ", (ftnlen)10, (ftnlen)10) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_mpole__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_mpole__ = FALSE_;
	}
    } else if (s_cmp(keyword, "POLARIZETERM ", (ftnlen)13, (ftnlen)13) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_polar__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_polar__ = FALSE_;
	}
    } else if (s_cmp(keyword, "RXNFIELDTERM ", (ftnlen)13, (ftnlen)13) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_rxnfld__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_rxnfld__ = FALSE_;
	}
    } else if (s_cmp(keyword, "SOLVATETERM ", (ftnlen)12, (ftnlen)12) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_solv__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_solv__ = FALSE_;
	}
    } else if (s_cmp(keyword, "METALTERM ", (ftnlen)12, (ftnlen)10) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_metal__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_metal__ = FALSE_;
	}
    } else if (s_cmp(keyword, "RESTRAINTERM ", (ftnlen)13, (ftnlen)13) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_geom__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_geom__ = FALSE_;
	}
    } else if (s_cmp(keyword, "EXTRATERM ", (ftnlen)10, (ftnlen)10) == 0) {
	getword_(record, value, &next, (ftnlen)120, (ftnlen)4);
	if (s_cmp(value, "ONLY", (ftnlen)4, (ftnlen)4) == 0) {
	    potoff_();
	}
	potent_1.use_extra__ = TRUE_;
	if (s_cmp(value, "NONE", (ftnlen)4, (ftnlen)4) == 0) {
	    potent_1.use_extra__ = FALSE_;
	}
    }

/*     select the name of the force field parameter set */

    if (s_cmp(keyword, "FORCEFIELD ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, fields_1.forcefield, &next, (ftnlen)120, (ftnlen)20);

/*     set control parameters for bond stretching potentials */

    } else if (s_cmp(keyword, "BONDTYPE ", (ftnlen)9, (ftnlen)9) == 0) {
	getword_(record, bndpot_1.bndtyp, &next, (ftnlen)120, (ftnlen)8);
    } else if (s_cmp(keyword, "BONDUNIT ", (ftnlen)9, (ftnlen)9) == 0) {
	i__1 = s_rsli(&io___6);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&bndpot_1.bndunit, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "BOND-CUBIC ", (ftnlen)11, (ftnlen)11) == 0) {
	i__1 = s_rsli(&io___7);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&bndpot_1.cbnd, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "BOND-QUARTIC ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___8);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&bndpot_1.qbnd, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}

/*     set control parameters for bond angle bending potentials */

    } else if (s_cmp(keyword, "ANGLEUNIT ", (ftnlen)10, (ftnlen)10) == 0) {
	i__1 = s_rsli(&io___9);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.angunit, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "ANGLE-CUBIC ", (ftnlen)12, (ftnlen)12) == 0) {
	i__1 = s_rsli(&io___10);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.cang, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "ANGLE-QUARTIC ", (ftnlen)14, (ftnlen)14) == 0) 
	    {
	i__1 = s_rsli(&io___11);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.qang, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "ANGLE-PENTIC ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___12);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.pang, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "ANGLE-SEXTIC ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___13);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.sang, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}

/*     set control parameters for stretch-bend potentials */

    } else if (s_cmp(keyword, "STRBNDUNIT ", (ftnlen)11, (ftnlen)11) == 0) {
	i__1 = s_rsli(&io___14);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.stbnunit, (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}

/*     set control parameters for Urey-Bradley potentials */

    } else if (s_cmp(keyword, "UREYUNIT ", (ftnlen)9, (ftnlen)9) == 0) {
	i__1 = s_rsli(&io___15);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&urypot_1.ureyunit, (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "UREY-CUBIC ", (ftnlen)11, (ftnlen)11) == 0) {
	i__1 = s_rsli(&io___16);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&urypot_1.cury, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "UREY-QUARTIC ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___17);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&urypot_1.qury, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}

/*     set control parameters for out-of-plane bend potentials */

    } else if (s_cmp(keyword, "OPBENDTYPE ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, angpot_1.opbtyp, &next, (ftnlen)120, (ftnlen)8);
    } else if (s_cmp(keyword, "OPBENDUNIT ", (ftnlen)11, (ftnlen)11) == 0) {
	i__1 = s_rsli(&io___18);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.opbunit, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "OPBEND-CUBIC ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___19);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.copb, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "OPBEND-QUARTIC ", (ftnlen)15, (ftnlen)15) == 0)
	     {
	i__1 = s_rsli(&io___20);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.qopb, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "OPBEND-PENTIC ", (ftnlen)14, (ftnlen)14) == 0) 
	    {
	i__1 = s_rsli(&io___21);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.popb, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "OPBEND-SEXTIC ", (ftnlen)14, (ftnlen)14) == 0) 
	    {
	i__1 = s_rsli(&io___22);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.sopb, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}

/*     set control parameters for out-of-plane distance potentials */

    } else if (s_cmp(keyword, "OPDISTUNIT ", (ftnlen)11, (ftnlen)11) == 0) {
	i__1 = s_rsli(&io___23);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.opdunit, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "OPDIST-CUBIC ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___24);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.copd, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "OPDIST-QUARTIC ", (ftnlen)15, (ftnlen)15) == 0)
	     {
	i__1 = s_rsli(&io___25);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.qopd, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "OPDIST-PENTIC ", (ftnlen)14, (ftnlen)14) == 0) 
	    {
	i__1 = s_rsli(&io___26);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.popd, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "OPDIST-SEXTIC ", (ftnlen)14, (ftnlen)14) == 0) 
	    {
	i__1 = s_rsli(&io___27);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.sopd, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}

/*     set control parameters for other local geometry potentials */

    } else if (s_cmp(keyword, "ANGANGUNIT ", (ftnlen)11, (ftnlen)11) == 0) {
	i__1 = s_rsli(&io___28);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&angpot_1.aaunit, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "IMPROPUNIT ", (ftnlen)11, (ftnlen)11) == 0) {
	i__1 = s_rsli(&io___29);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&torpot_1.idihunit, (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "IMPTORUNIT ", (ftnlen)11, (ftnlen)11) == 0) {
	i__1 = s_rsli(&io___30);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&torpot_1.itorunit, (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "TORSIONUNIT ", (ftnlen)12, (ftnlen)12) == 0) {
	i__1 = s_rsli(&io___31);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&torpot_1.torsunit, (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "PITORSUNIT ", (ftnlen)11, (ftnlen)11) == 0) {
	i__1 = s_rsli(&io___32);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&torpot_1.ptorunit, (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "STRTORUNIT ", (ftnlen)11, (ftnlen)11) == 0) {
	i__1 = s_rsli(&io___33);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&torpot_1.storunit, (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "TORTORUNIT ", (ftnlen)11, (ftnlen)11) == 0) {
	i__1 = s_rsli(&io___34);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&torpot_1.ttorunit, (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}

/*     set control parameters for van der Waals potentials */

    } else if (s_cmp(keyword, "VDWINDEX ", (ftnlen)9, (ftnlen)9) == 0) {
	getword_(record, vdwpot_1.vdwindex, &next, (ftnlen)120, (ftnlen)5);
    } else if (s_cmp(keyword, "VDWTYPE ", (ftnlen)8, (ftnlen)8) == 0) {
	getword_(record, vdwpot_1.vdwtyp, &next, (ftnlen)120, (ftnlen)13);
    } else if (s_cmp(keyword, "RADIUSTYPE ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, vdwpot_1.radtyp, &next, (ftnlen)120, (ftnlen)5);
    } else if (s_cmp(keyword, "RADIUSSIZE ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, vdwpot_1.radsiz, &next, (ftnlen)120, (ftnlen)8);
    } else if (s_cmp(keyword, "RADIUSRULE ", (ftnlen)11, (ftnlen)11) == 0) {
	getword_(record, vdwpot_1.radrule, &next, (ftnlen)120, (ftnlen)10);
    } else if (s_cmp(keyword, "EPSILONRULE ", (ftnlen)12, (ftnlen)12) == 0) {
	getword_(record, vdwpot_1.epsrule, &next, (ftnlen)120, (ftnlen)10);
    } else if (s_cmp(keyword, "GAUSSTYPE ", (ftnlen)14, (ftnlen)10) == 0) {
	getword_(record, vdwpot_1.gausstyp, &next, (ftnlen)120, (ftnlen)8);
    } else if (s_cmp(keyword, "A-EXPTERM ", (ftnlen)10, (ftnlen)10) == 0) {
	i__1 = s_rsli(&io___35);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&vdwpot_1.abuck, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "B-EXPTERM ", (ftnlen)10, (ftnlen)10) == 0) {
	i__1 = s_rsli(&io___36);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&vdwpot_1.bbuck, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "C-EXPTERM ", (ftnlen)10, (ftnlen)10) == 0) {
	i__1 = s_rsli(&io___37);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&vdwpot_1.cbuck, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "GAMMA-HALGREN ", (ftnlen)14, (ftnlen)14) == 0) 
	    {
	i__1 = s_rsli(&io___38);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&vdwpot_1.ghal, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "DELTA-HALGREN ", (ftnlen)14, (ftnlen)14) == 0) 
	    {
	i__1 = s_rsli(&io___39);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&vdwpot_1.dhal, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "VDW-12-SCALE ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___40);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&vdwpot_1.v2scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (vdwpot_1.v2scale > 1.) {
	    vdwpot_1.v2scale = 1. / vdwpot_1.v2scale;
	}
    } else if (s_cmp(keyword, "VDW-13-SCALE ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___41);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&vdwpot_1.v3scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (vdwpot_1.v3scale > 1.) {
	    vdwpot_1.v3scale = 1. / vdwpot_1.v3scale;
	}
    } else if (s_cmp(keyword, "VDW-14-SCALE ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___42);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&vdwpot_1.v4scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (vdwpot_1.v4scale > 1.) {
	    vdwpot_1.v4scale = 1. / vdwpot_1.v4scale;
	}
    } else if (s_cmp(keyword, "VDW-15-SCALE ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___43);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&vdwpot_1.v5scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (vdwpot_1.v5scale > 1.) {
	    vdwpot_1.v5scale = 1. / vdwpot_1.v5scale;
	}

/*     set control parameters for charge-charge potentials */

    } else if (s_cmp(keyword, "ELECTRIC ", (ftnlen)9, (ftnlen)9) == 0) {
	i__1 = s_rsli(&io___44);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&chgpot_1.electric, (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "DIELECTRIC ", (ftnlen)11, (ftnlen)11) == 0) {
	i__1 = s_rsli(&io___45);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&chgpot_1.dielec, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "CHG-BUFFER ", (ftnlen)11, (ftnlen)11) == 0) {
	i__1 = s_rsli(&io___46);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&chgpot_1.ebuffer, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "NEUTRAL-GROUPS ", (ftnlen)15, (ftnlen)15) == 0)
	     {
	chgpot_1.neutcut = TRUE_;
    } else if (s_cmp(keyword, "NEIGHBOR-GROUPS ", (ftnlen)16, (ftnlen)16) == 
	    0) {
	chgpot_1.neutnbr = TRUE_;
    } else if (s_cmp(keyword, "CHG-12-SCALE ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___47);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&chgpot_1.c2scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (chgpot_1.c2scale > 1.) {
	    chgpot_1.c2scale = 1. / chgpot_1.c2scale;
	}
    } else if (s_cmp(keyword, "CHG-13-SCALE ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___48);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&chgpot_1.c3scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (chgpot_1.c3scale > 1.) {
	    chgpot_1.c3scale = 1. / chgpot_1.c3scale;
	}
    } else if (s_cmp(keyword, "CHG-14-SCALE ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___49);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&chgpot_1.c4scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (chgpot_1.c4scale > 1.) {
	    chgpot_1.c4scale = 1. / chgpot_1.c4scale;
	}
    } else if (s_cmp(keyword, "CHG-15-SCALE ", (ftnlen)13, (ftnlen)13) == 0) {
	i__1 = s_rsli(&io___50);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&chgpot_1.c5scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (chgpot_1.c5scale > 1.) {
	    chgpot_1.c5scale = 1. / chgpot_1.c5scale;
	}

/*     set control parameters for atomic multipole potentials */

    } else if (s_cmp(keyword, "MPOLE-12-SCALE ", (ftnlen)15, (ftnlen)15) == 0)
	     {
	i__1 = s_rsli(&io___51);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&mplpot_1.m2scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (mplpot_1.m2scale > 1.) {
	    mplpot_1.m2scale = 1. / mplpot_1.m2scale;
	}
    } else if (s_cmp(keyword, "MPOLE-13-SCALE ", (ftnlen)15, (ftnlen)15) == 0)
	     {
	i__1 = s_rsli(&io___52);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&mplpot_1.m3scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (mplpot_1.m3scale > 1.) {
	    mplpot_1.m3scale = 1. / mplpot_1.m3scale;
	}
    } else if (s_cmp(keyword, "MPOLE-14-SCALE ", (ftnlen)15, (ftnlen)15) == 0)
	     {
	i__1 = s_rsli(&io___53);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&mplpot_1.m4scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (mplpot_1.m4scale > 1.) {
	    mplpot_1.m4scale = 1. / mplpot_1.m4scale;
	}
    } else if (s_cmp(keyword, "MPOLE-15-SCALE ", (ftnlen)15, (ftnlen)15) == 0)
	     {
	i__1 = s_rsli(&io___54);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&mplpot_1.m5scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (mplpot_1.m5scale > 1.) {
	    mplpot_1.m5scale = 1. / mplpot_1.m5scale;
	}

/*     set control parameters for polarization potentials */

    } else if (s_cmp(keyword, "POLARIZATION ", (ftnlen)13, (ftnlen)13) == 0) {
	getword_(record, polpot_1.poltyp, &next, (ftnlen)120, (ftnlen)6);
    } else if (s_cmp(keyword, "POLAR-EPS ", (ftnlen)10, (ftnlen)10) == 0) {
	i__1 = s_rsli(&io___55);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.poleps, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "POLAR-SOR ", (ftnlen)10, (ftnlen)10) == 0) {
	i__1 = s_rsli(&io___56);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.polsor, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    } else if (s_cmp(keyword, "POLAR-12-SCALE ", (ftnlen)15, (ftnlen)15) == 0)
	     {
	i__1 = s_rsli(&io___57);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.p2scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (polpot_1.p2scale > 1.) {
	    polpot_1.p2scale = 1. / polpot_1.p2scale;
	}
    } else if (s_cmp(keyword, "POLAR-13-SCALE ", (ftnlen)15, (ftnlen)15) == 0)
	     {
	i__1 = s_rsli(&io___58);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.p3scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (polpot_1.p3scale > 1.) {
	    polpot_1.p3scale = 1. / polpot_1.p3scale;
	}
    } else if (s_cmp(keyword, "POLAR-14-SCALE ", (ftnlen)15, (ftnlen)15) == 0)
	     {
	i__1 = s_rsli(&io___59);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.p4scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (polpot_1.p4scale > 1.) {
	    polpot_1.p4scale = 1. / polpot_1.p4scale;
	}
    } else if (s_cmp(keyword, "POLAR-15-SCALE ", (ftnlen)15, (ftnlen)15) == 0)
	     {
	i__1 = s_rsli(&io___60);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.p5scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (polpot_1.p5scale > 1.) {
	    polpot_1.p5scale = 1. / polpot_1.p5scale;
	}
    } else if (s_cmp(keyword, "DIRECT-11-SCALE ", (ftnlen)16, (ftnlen)16) == 
	    0) {
	i__1 = s_rsli(&io___61);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.d1scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (polpot_1.d1scale > 1.) {
	    polpot_1.d1scale = 1. / polpot_1.d1scale;
	}
    } else if (s_cmp(keyword, "DIRECT-12-SCALE ", (ftnlen)16, (ftnlen)16) == 
	    0) {
	i__1 = s_rsli(&io___62);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.d2scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (polpot_1.d2scale > 1.) {
	    polpot_1.d2scale = 1. / polpot_1.d2scale;
	}
    } else if (s_cmp(keyword, "DIRECT-13-SCALE ", (ftnlen)16, (ftnlen)16) == 
	    0) {
	i__1 = s_rsli(&io___63);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.d3scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (polpot_1.d3scale > 1.) {
	    polpot_1.d3scale = 1. / polpot_1.d3scale;
	}
    } else if (s_cmp(keyword, "DIRECT-14-SCALE ", (ftnlen)16, (ftnlen)16) == 
	    0) {
	i__1 = s_rsli(&io___64);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.d4scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (polpot_1.d4scale > 1.) {
	    polpot_1.d4scale = 1. / polpot_1.d4scale;
	}
    } else if (s_cmp(keyword, "MUTUAL-11-SCALE ", (ftnlen)16, (ftnlen)16) == 
	    0) {
	i__1 = s_rsli(&io___65);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.u1scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (polpot_1.u1scale > 1.) {
	    polpot_1.u1scale = 1. / polpot_1.u1scale;
	}
    } else if (s_cmp(keyword, "MUTUAL-12-SCALE ", (ftnlen)16, (ftnlen)16) == 
	    0) {
	i__1 = s_rsli(&io___66);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.u2scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (polpot_1.u2scale > 1.) {
	    polpot_1.u2scale = 1. / polpot_1.u2scale;
	}
    } else if (s_cmp(keyword, "MUTUAL-13-SCALE ", (ftnlen)16, (ftnlen)16) == 
	    0) {
	i__1 = s_rsli(&io___67);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.u3scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (polpot_1.u3scale > 1.) {
	    polpot_1.u3scale = 1. / polpot_1.u3scale;
	}
    } else if (s_cmp(keyword, "MUTUAL-14-SCALE ", (ftnlen)16, (ftnlen)16) == 
	    0) {
	i__1 = s_rsli(&io___68);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&polpot_1.u4scale, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	if (polpot_1.u4scale > 1.) {
	    polpot_1.u4scale = 1. / polpot_1.u4scale;
	}

/*     set control parameters for reaction field potentials */

    } else if (s_cmp(keyword, "REACTIONFIELD ", (ftnlen)14, (ftnlen)14) == 0) 
	    {
	i__1 = s_rsli(&io___69);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&rxnpot_1.rfsize, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&rxnpot_1.rfbulkd, (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&rxnpot_1.rfterms, (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
    }

/*     jump directly to the end if any error was detected */

L10:
    return 0;
} /* prmkey_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine potoff  --  turn off all potential functions  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "potoff" clears the forcefield definition by turning off */
/*     the use of each of the potential energy functions */


/* Subroutine */ int potoff_(void)
{


/*     turn off the use of each of the potential energy functions */



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


    potent_1.use_bond__ = FALSE_;
    potent_1.use_angle__ = FALSE_;
    potent_1.use_strbnd__ = FALSE_;
    potent_1.use_urey__ = FALSE_;
    potent_1.use_angang__ = FALSE_;
    potent_1.use_opbend__ = FALSE_;
    potent_1.use_opdist__ = FALSE_;
    potent_1.use_improp__ = FALSE_;
    potent_1.use_imptor__ = FALSE_;
    potent_1.use_tors__ = FALSE_;
    potent_1.use_pitors__ = FALSE_;
    potent_1.use_strtor__ = FALSE_;
    potent_1.use_tortor__ = FALSE_;
    potent_1.use_vdw__ = FALSE_;
    potent_1.use_charge__ = FALSE_;
    potent_1.use_chgdpl__ = FALSE_;
    potent_1.use_dipole__ = FALSE_;
    potent_1.use_mpole__ = FALSE_;
    potent_1.use_polar__ = FALSE_;
    potent_1.use_rxnfld__ = FALSE_;
    potent_1.use_solv__ = FALSE_;
    potent_1.use_metal__ = FALSE_;
    potent_1.use_geom__ = FALSE_;
    potent_1.use_extra__ = FALSE_;
    return 0;
} /* potoff_ */

