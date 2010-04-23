/* kvdw.f -- translated by f2c (version 20050501).
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
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

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
    doublereal radhb[500], epshb[500];
    char khb[4000];
} khbond_;

#define khbond_1 khbond_

struct {
    doublereal rad[5000], eps[5000], rad4[5000], eps4[5000], reduct[5000];
} kvdws_;

#define kvdws_1 kvdws_

struct {
    doublereal radpr[500], epspr[500];
    char kvpr[4000];
} kvdwpr_;

#define kvdwpr_1 kvdwpr_

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
    doublereal radmin[1000000]	/* was [1000][1000] */, epsilon[1000000]	
	    /* was [1000][1000] */, radmin4[1000000]	/* was [1000][1000] */
	    , epsilon4[1000000]	/* was [1000][1000] */, radhbnd[1000000]	
	    /* was [1000][1000] */, epshbnd[1000000]	/* was [1000][1000] */
	    , kred[25000];
    integer ired[25000], nvdw, ivdw[25000], jvdw[25000], nvt, ivt[25000], jvt[
	    25000];
} vdw_;

#define vdw_1 vdw_

struct {
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
    logical use_vcorr__;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__1000 = 1000;
static integer c__3 = 3;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine kvdw  --  van der Waals parameter assignment  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "kvdw" assigns the parameters to be used in computing the */
/*     van der Waals interactions and processes any new or changed */
/*     values for these parameters */


/* Subroutine */ int kvdw_(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Additional van der Waals Parameters :"
	    "\002,//,5x,\002Atom Class\002,10x,\002Size\002,6x,\002Epsilon"
	    "\002,5x,\002Reduction\002,/)";
    static char fmt_30[] = "(/,\002 Additional van der Waals Parameters :"
	    "\002,//,5x,\002Atom Type\002,11x,\002Size\002,6x,\002Epsilon\002"
	    ",5x,\002Reduction\002,/)";
    static char fmt_40[] = "(4x,i6,8x,2f12.4,f12.3)";
    static char fmt_50[] = "(/,\002 KVDW  --  Only Atom Classes through\002,"
	    "i4,\002 are Allowed\002)";
    static char fmt_70[] = "(/,\002 Additional 1-4 van der Waals\002,\002 Pa"
	    "rameters :\002,//,5x,\002Atom Class\002,10x,\002Size\002,6x,\002"
	    "Epsilon\002,/)";
    static char fmt_80[] = "(/,\002 Additional 1-4 van der Waals\002,\002 Pa"
	    "rameters :\002,//,5x,\002Atom Type\002,11x,\002Size\002,6x,\002E"
	    "psilon\002,/)";
    static char fmt_90[] = "(4x,i6,8x,2f12.4)";
    static char fmt_100[] = "(/,\002 KVDW  --  Only Atom Classes through\002"
	    ",i4,\002 are Allowed\002)";
    static char fmt_110[] = "(/,\002 Additional van der Waals Parameter"
	    "s\002,\002 for Specific Pairs :\002,//,5x,\002Atom Classes\002,6"
	    "x,\002Size Sum\002,4x,\002Epsilon\002,/)";
    static char fmt_120[] = "(/,\002 Additional van der Waals Parameter"
	    "s\002,\002 for Specific Pairs :\002,//,5x,\002Atom Types\002,8x"
	    ",\002Size Sum\002,4x,\002Epsilon\002,/)";
    static char fmt_130[] = "(6x,2i4,4x,2f12.4)";
    static char fmt_140[] = "(/,\002 KVDW  --  Too many Special VDW Pai"
	    "r\002,\002 Parameters\002)";
    static char fmt_160[] = "(/,\002 Additional van der Waals Hydrogen\002"
	    ",\002 Bonding Parameters :\002,//,5x,\002Atom Classes\002,6x,"
	    "\002Size Sum\002,4x,\002Epsilon\002,/)";
    static char fmt_170[] = "(/,\002 Additional van der Waals Hydrogen\002"
	    ",\002 Bonding Parameters :\002,//,5x,\002Atom Types\002,8x,\002S"
	    "ize Sum\002,4x,\002Epsilon\002,/)";
    static char fmt_180[] = "(6x,2i4,4x,2f12.4)";
    static char fmt_190[] = "(/,\002 KVDW  --  Too many Hydrogen Bonding P"
	    "air\002,\002 Parameters\002)";
    static char fmt_210[] = "(/,\002 KVDW  --  Unable to Index VDW Parameter"
	    "s;\002,\002 Increase MAXCLASS\002)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2], i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, k, ia, ib;
    static doublereal ep, rd;
    static char pa[4], pb[4];
    static integer it;
    static char pt[8];
    static doublereal rdn, srad[5000];
    static integer size;
    static doublereal seps[5000];
    static integer next;
    static doublereal srad4[5000], seps4[5000];
    static char blank[8];
    static logical header;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    extern integer number_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int getnumb_(char *, integer *, integer *, ftnlen)
	    , numeral_(integer *, char *, integer *, ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___12 = { 1, string, 1, 0, 120, 1 };
    static cilist io___13 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___17 = { 1, string, 1, 0, 120, 1 };
    static cilist io___18 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_100, 0 };
    static icilist io___24 = { 1, string, 1, 0, 120, 1 };
    static cilist io___25 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_140, 0 };
    static icilist io___33 = { 1, string, 1, 0, 120, 1 };
    static cilist io___34 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_210, 0 };



#define epsilon4_ref(a_1,a_2) vdw_1.epsilon4[(a_2)*1000 + a_1 - 1001]
#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define khb_ref(a_0,a_1) &khbond_1.khb[(a_1)*8 + a_0 - 8]
#define kvpr_ref(a_0,a_1) &kvdwpr_1.kvpr[(a_1)*8 + a_0 - 8]
#define radmin_ref(a_1,a_2) vdw_1.radmin[(a_2)*1000 + a_1 - 1001]
#define igauss_ref(a_1,a_2) vdwpot_1.igauss[(a_2)*2 + a_1 - 3]
#define radmin4_ref(a_1,a_2) vdw_1.radmin4[(a_2)*1000 + a_1 - 1001]
#define radhbnd_ref(a_1,a_2) vdw_1.radhbnd[(a_2)*1000 + a_1 - 1001]
#define epshbnd_ref(a_1,a_2) vdw_1.epshbnd[(a_2)*1000 + a_1 - 1001]
#define keyline_ref(a_0,a_1) &keys_1.keyline[(a_1)*120 + a_0 - 120]
#define epsilon_ref(a_1,a_2) vdw_1.epsilon[(a_2)*1000 + a_1 - 1001]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  keys.i  --  contents of current keyword parameter file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nkey      number of nonblank lines in the keyword file */
/*     keyline   contents of each individual keyword file line */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  khbond.i  --  forcefield parameters for H-bonding terms  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxnhb   maximum number of hydrogen bonding pair entries */

/*     radhb    radius parameter for hydrogen bonding pairs */
/*     epshb    well depth parameter for hydrogen bonding pairs */
/*     khb      string of atom types for hydrogen bonding pairs */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kvdws.i  --  forcefield parameters for van der Waals terms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     rad      van der Waals radius parameter for each atom type */
/*     eps      van der Waals well depth parameter for each atom type */
/*     rad4     van der Waals radius parameter in 1-4 interactions */
/*     eps4     van der Waals well depth parameter in 1-4 interactions */
/*     reduct   van der Waals reduction factor for each atom type */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kvdwpr.i  --  forcefield parameters for special vdw terms  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxnvp   maximum number of special van der Waals pair entries */

/*     radpr    radius parameter for special van der Waals pairs */
/*     epspr    well depth parameter for special van der Waals pairs */
/*     kvpr     string of atom classes for special van der Waals pairs */




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
/*     ##  vdw.i  --  van der Waals parameters for current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     radmin     minimum energy distance for each atom class pair */
/*     epsilon    well depth parameter for each atom class pair */
/*     radmin4    minimum energy distance for 1-4 interaction pairs */
/*     epsilon4   well depth parameter for 1-4 interaction pairs */
/*     radhbnd    minimum energy distance for hydrogen bonding pairs */
/*     epshbnd    well depth parameter for hydrogen bonding pairs */
/*     kred       value of reduction factor parameter for each atom */
/*     ired       attached atom from which reduction factor is applied */
/*     nvdw       total number van der Waals active sites in the system */
/*     ivdw       number of the atom for each van der Waals active site */
/*     jvdw       type or class index into vdw parameters for each atom */
/*     nvt        number of distinct vdw types/classes in the system */
/*     ivt        type/class index for each distinct vdw type or class */
/*     jvt        frequency of each vdw type or class in the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  vdwpot.i  --  specifics of van der Waals functional form  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     abuck       value of "A" constant in Buckingham vdw potential */
/*     bbuck       value of "B" constant in Buckingham vdw potential */
/*     cbuck       value of "C" constant in Buckingham vdw potential */
/*     ghal        value of "gamma" in buffered 14-7 vdw potential */
/*     dhal        value of "delta" in buffered 14-7 vdw potential */
/*     v2scale     factor by which 1-2 vdw interactions are scaled */
/*     v3scale     factor by which 1-3 vdw interactions are scaled */
/*     v4scale     factor by which 1-4 vdw interactions are scaled */
/*     v5scale     factor by which 1-5 vdw interactions are scaled */
/*     igauss      coefficients of Gaussian fit to vdw potential */
/*     ngauss      number of Gaussians used in fit to vdw potential */
/*     use_vcorr   flag to use long range vdw der Waals correction */
/*     vdwindex    indexing mode (atom type or class) for vdw parameters */
/*     vdwtyp      type of van der Waals potential energy function */
/*     radtyp      type of parameter (sigma or R-min) for atomic size */
/*     radsiz      atomic size provided as radius or diameter */
/*     radrule     combining rule for atomic size parameters */
/*     epsrule     combining rule for vdw well depth parameters */
/*     gausstyp    type of Gaussian fit to van der Waals potential */




/*     process keywords containing van der Waals parameters */

    s_copy(blank, "        ", (ftnlen)8, (ftnlen)8);
    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "VDW ", (ftnlen)4, (ftnlen)4) == 0) {
	    getnumb_(record, &k, &next, (ftnlen)120);
	    if (k >= 1 && k <= 1000) {
		rd = kvdws_1.rad[k - 1];
		ep = kvdws_1.eps[k - 1];
		rdn = kvdws_1.reduct[k - 1];
		s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next 
			- 1));
		i__2 = s_rsli(&io___12);
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&rd, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&ep, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&rdn, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L10;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L10;
		}
L10:
		if (header) {
		    header = FALSE_;
		    if (s_cmp(vdwpot_1.vdwindex, "CLASS", (ftnlen)5, (ftnlen)
			    5) == 0) {
			io___13.ciunit = iounit_1.iout;
			s_wsfe(&io___13);
			e_wsfe();
		    } else {
			io___14.ciunit = iounit_1.iout;
			s_wsfe(&io___14);
			e_wsfe();
		    }
		}
		kvdws_1.rad[k - 1] = rd;
		kvdws_1.eps[k - 1] = ep;
		kvdws_1.reduct[k - 1] = rdn;
		io___15.ciunit = iounit_1.iout;
		s_wsfe(&io___15);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&rd, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&ep, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&rdn, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else if (k > 1000) {
		io___16.ciunit = iounit_1.iout;
		s_wsfe(&io___16);
		do_fio(&c__1, (char *)&c__1000, (ftnlen)sizeof(integer));
		e_wsfe();
		inform_1.abort = TRUE_;
	    }
	}
    }

/*     process keywords containing 1-4 van der Waals parameters */

    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "VDW14 ", (ftnlen)6, (ftnlen)6) == 0) {
	    getnumb_(record, &k, &next, (ftnlen)120);
	    if (k >= 1 && k <= 1000) {
		rd = kvdws_1.rad4[k - 1];
		ep = kvdws_1.eps4[k - 1];
		s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next 
			- 1));
		i__2 = s_rsli(&io___17);
		if (i__2 != 0) {
		    goto L60;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&rd, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L60;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&ep, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L60;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L60;
		}
L60:
		if (header) {
		    header = FALSE_;
		    if (s_cmp(vdwpot_1.vdwindex, "CLASS", (ftnlen)5, (ftnlen)
			    5) == 0) {
			io___18.ciunit = iounit_1.iout;
			s_wsfe(&io___18);
			e_wsfe();
		    } else {
			io___19.ciunit = iounit_1.iout;
			s_wsfe(&io___19);
			e_wsfe();
		    }
		}
		kvdws_1.rad4[k - 1] = rd;
		kvdws_1.eps4[k - 1] = ep;
		io___20.ciunit = iounit_1.iout;
		s_wsfe(&io___20);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&rd, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&ep, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else if (k > 1000) {
		io___21.ciunit = iounit_1.iout;
		s_wsfe(&io___21);
		do_fio(&c__1, (char *)&c__1000, (ftnlen)sizeof(integer));
		e_wsfe();
		inform_1.abort = TRUE_;
	    }
	}
    }

/*     process keywords containing specific pair vdw parameters */

    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "VDWPR ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    rd = 0.;
	    ep = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___24);
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&rd, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&ep, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L150;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L150;
	    }
	    if (header) {
		header = FALSE_;
		if (s_cmp(vdwpot_1.vdwindex, "CLASS", (ftnlen)5, (ftnlen)5) ==
			 0) {
		    io___25.ciunit = iounit_1.iout;
		    s_wsfe(&io___25);
		    e_wsfe();
		} else {
		    io___26.ciunit = iounit_1.iout;
		    s_wsfe(&io___26);
		    e_wsfe();
		}
	    }
	    io___27.ciunit = iounit_1.iout;
	    s_wsfe(&io___27);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&rd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ep, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    size = 4;
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    if (ia <= ib) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pb;
		i__3[1] = 4, a__1[1] = pa;
		s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
	    }
	    for (k = 1; k <= 500; ++k) {
		if (s_cmp(kvpr_ref(0, k), blank, (ftnlen)8, (ftnlen)8) == 0 ||
			 s_cmp(kvpr_ref(0, k), pt, (ftnlen)8, (ftnlen)8) == 0)
			 {
		    s_copy(kvpr_ref(0, k), pt, (ftnlen)8, (ftnlen)8);
		    kvdwpr_1.radpr[k - 1] = rd;
		    kvdwpr_1.epspr[k - 1] = ep;
		    goto L150;
		}
	    }
	    io___32.ciunit = iounit_1.iout;
	    s_wsfe(&io___32);
	    e_wsfe();
	    inform_1.abort = TRUE_;
L150:
	    ;
	}
    }

/*     process keywords containing hydrogen bonding vdw parameters */

    header = TRUE_;
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "HBOND ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    rd = 0.;
	    ep = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___33);
	    if (i__2 != 0) {
		goto L200;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L200;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L200;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&rd, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L200;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&ep, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L200;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L200;
	    }
	    if (header) {
		header = FALSE_;
		if (s_cmp(vdwpot_1.vdwindex, "CLASS", (ftnlen)5, (ftnlen)5) ==
			 0) {
		    io___34.ciunit = iounit_1.iout;
		    s_wsfe(&io___34);
		    e_wsfe();
		} else {
		    io___35.ciunit = iounit_1.iout;
		    s_wsfe(&io___35);
		    e_wsfe();
		}
	    }
	    io___36.ciunit = iounit_1.iout;
	    s_wsfe(&io___36);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&rd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ep, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    size = 4;
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    if (ia <= ib) {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pa;
		i__3[1] = 4, a__1[1] = pb;
		s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__1[0] = pb;
		i__3[1] = 4, a__1[1] = pa;
		s_cat(pt, a__1, i__3, &c__2, (ftnlen)8);
	    }
	    for (k = 1; k <= 500; ++k) {
		if (s_cmp(khb_ref(0, k), blank, (ftnlen)8, (ftnlen)8) == 0 || 
			s_cmp(khb_ref(0, k), pt, (ftnlen)8, (ftnlen)8) == 0) {
		    s_copy(khb_ref(0, k), pt, (ftnlen)8, (ftnlen)8);
		    khbond_1.radhb[k - 1] = rd;
		    khbond_1.epshb[k - 1] = ep;
		    goto L200;
		}
	    }
	    io___37.ciunit = iounit_1.iout;
	    s_wsfe(&io___37);
	    e_wsfe();
	    inform_1.abort = TRUE_;
L200:
	    ;
	}
    }

/*     use atom class or type as index into vdw parameters */

    k = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vdw_1.jvdw[i__ - 1] = atmtyp_1.class__[i__ - 1];
	if (s_cmp(vdwpot_1.vdwindex, "TYPE", (ftnlen)5, (ftnlen)4) == 0) {
	    vdw_1.jvdw[i__ - 1] = atoms_1.type__[i__ - 1];
	}
/* Computing MAX */
	i__2 = k, i__4 = vdw_1.jvdw[i__ - 1];
	k = max(i__2,i__4);
    }
    if (k > 1000) {
	io___38.ciunit = iounit_1.iout;
	s_wsfe(&io___38);
	e_wsfe();
	inform_1.abort = TRUE_;
    }

/*     count the number of vdw types and their frequencies */

    vdw_1.nvt = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	it = vdw_1.jvdw[i__ - 1];
	i__2 = vdw_1.nvt;
	for (k = 1; k <= i__2; ++k) {
	    if (vdw_1.ivt[k - 1] == it) {
		++vdw_1.jvt[k - 1];
		goto L220;
	    }
	}
	++vdw_1.nvt;
	vdw_1.ivt[vdw_1.nvt - 1] = it;
	vdw_1.jvt[vdw_1.nvt - 1] = 1;
L220:
	;
    }

/*     get the vdw radii and well depths for each atom type */

    for (i__ = 1; i__ <= 5000; ++i__) {
	if (kvdws_1.rad4[i__ - 1] == 0.) {
	    kvdws_1.rad4[i__ - 1] = kvdws_1.rad[i__ - 1];
	}
	if (kvdws_1.eps4[i__ - 1] == 0.) {
	    kvdws_1.eps4[i__ - 1] = kvdws_1.eps[i__ - 1];
	}
	if (s_cmp(vdwpot_1.radtyp, "SIGMA", (ftnlen)5, (ftnlen)5) == 0) {
	    kvdws_1.rad[i__ - 1] *= 1.122462048309372981;
	    kvdws_1.rad4[i__ - 1] *= 1.122462048309372981;
	}
	if (s_cmp(vdwpot_1.radsiz, "DIAMETER", (ftnlen)8, (ftnlen)8) == 0) {
	    kvdws_1.rad[i__ - 1] *= .5;
	    kvdws_1.rad4[i__ - 1] *= .5;
	}
	srad[i__ - 1] = sqrt(kvdws_1.rad[i__ - 1]);
	kvdws_1.eps[i__ - 1] = (d__1 = kvdws_1.eps[i__ - 1], abs(d__1));
	seps[i__ - 1] = sqrt(kvdws_1.eps[i__ - 1]);
	srad4[i__ - 1] = sqrt(kvdws_1.rad4[i__ - 1]);
	kvdws_1.eps4[i__ - 1] = (d__1 = kvdws_1.eps4[i__ - 1], abs(d__1));
	seps4[i__ - 1] = sqrt(kvdws_1.eps4[i__ - 1]);
    }

/*     use combination rules to set pairwise vdw radii sums */

    for (i__ = 1; i__ <= 1000; ++i__) {
	for (k = i__; k <= 1000; ++k) {
	    if (kvdws_1.rad[i__ - 1] == 0. && kvdws_1.rad[k - 1] == 0.) {
		rd = 0.;
	    } else if (s_cmp(vdwpot_1.radrule, "ARITHMETIC", (ftnlen)10, (
		    ftnlen)10) == 0) {
		rd = kvdws_1.rad[i__ - 1] + kvdws_1.rad[k - 1];
	    } else if (s_cmp(vdwpot_1.radrule, "GEOMETRIC", (ftnlen)9, (
		    ftnlen)9) == 0) {
		rd = srad[i__ - 1] * srad[k - 1] * 2.;
	    } else if (s_cmp(vdwpot_1.radrule, "CUBIC-MEAN", (ftnlen)10, (
		    ftnlen)10) == 0) {
/* Computing 3rd power */
		d__1 = kvdws_1.rad[i__ - 1];
/* Computing 3rd power */
		d__2 = kvdws_1.rad[k - 1];
/* Computing 2nd power */
		d__3 = kvdws_1.rad[i__ - 1];
/* Computing 2nd power */
		d__4 = kvdws_1.rad[k - 1];
		rd = (d__1 * (d__1 * d__1) + d__2 * (d__2 * d__2)) * 2. / (
			d__3 * d__3 + d__4 * d__4);
	    } else {
		rd = kvdws_1.rad[i__ - 1] + kvdws_1.rad[k - 1];
	    }
	    radmin_ref(i__, k) = rd;
	    radmin_ref(k, i__) = rd;
	}
    }

/*     use combination rules to set pairwise well depths */

    for (i__ = 1; i__ <= 1000; ++i__) {
	for (k = i__; k <= 1000; ++k) {
	    if (kvdws_1.eps[i__ - 1] == 0. && kvdws_1.eps[k - 1] == 0.) {
		ep = 0.;
	    } else if (s_cmp(vdwpot_1.epsrule, "ARITHMETIC", (ftnlen)10, (
		    ftnlen)10) == 0) {
		ep = (kvdws_1.eps[i__ - 1] + kvdws_1.eps[k - 1]) * .5;
	    } else if (s_cmp(vdwpot_1.epsrule, "GEOMETRIC", (ftnlen)9, (
		    ftnlen)9) == 0) {
		ep = seps[i__ - 1] * seps[k - 1];
	    } else if (s_cmp(vdwpot_1.epsrule, "HARMONIC", (ftnlen)8, (ftnlen)
		    8) == 0) {
		ep = kvdws_1.eps[i__ - 1] * kvdws_1.eps[k - 1] * 2. / (
			kvdws_1.eps[i__ - 1] + kvdws_1.eps[k - 1]);
	    } else if (s_cmp(vdwpot_1.epsrule, "HHG", (ftnlen)3, (ftnlen)3) ==
		     0) {
/* Computing 2nd power */
		d__1 = seps[i__ - 1] + seps[k - 1];
		ep = kvdws_1.eps[i__ - 1] * kvdws_1.eps[k - 1] * 4. / (d__1 * 
			d__1);
	    } else {
		ep = seps[i__ - 1] * seps[k - 1];
	    }
	    epsilon_ref(i__, k) = ep;
	    epsilon_ref(k, i__) = ep;
	}
    }

/*     use combination rules to set pairwise 1-4 vdw radii sums */

    for (i__ = 1; i__ <= 1000; ++i__) {
	for (k = i__; k <= 1000; ++k) {
	    if (kvdws_1.rad4[i__ - 1] == 0. && kvdws_1.rad4[k - 1] == 0.) {
		rd = 0.;
	    } else if (s_cmp(vdwpot_1.radrule, "ARITHMETIC", (ftnlen)10, (
		    ftnlen)10) == 0) {
		rd = kvdws_1.rad4[i__ - 1] + kvdws_1.rad4[k - 1];
	    } else if (s_cmp(vdwpot_1.radrule, "GEOMETRIC", (ftnlen)9, (
		    ftnlen)9) == 0) {
		rd = srad4[i__ - 1] * srad4[k - 1] * 2.;
	    } else if (s_cmp(vdwpot_1.radrule, "CUBIC-MEAN", (ftnlen)10, (
		    ftnlen)10) == 0) {
/* Computing 3rd power */
		d__1 = kvdws_1.rad4[i__ - 1];
/* Computing 3rd power */
		d__2 = kvdws_1.rad4[k - 1];
/* Computing 2nd power */
		d__3 = kvdws_1.rad4[i__ - 1];
/* Computing 2nd power */
		d__4 = kvdws_1.rad4[k - 1];
		rd = (d__1 * (d__1 * d__1) + d__2 * (d__2 * d__2)) * 2. / (
			d__3 * d__3 + d__4 * d__4);
	    } else {
		rd = kvdws_1.rad4[i__ - 1] + kvdws_1.rad4[k - 1];
	    }
	    radmin4_ref(i__, k) = rd;
	    radmin4_ref(k, i__) = rd;
	}
    }

/*     use combination rules to set pairwise 1-4 well depths */

    for (i__ = 1; i__ <= 1000; ++i__) {
	for (k = i__; k <= 1000; ++k) {
	    if (kvdws_1.eps4[i__ - 1] == 0. && kvdws_1.eps4[k - 1] == 0.) {
		ep = 0.;
	    } else if (s_cmp(vdwpot_1.epsrule, "ARITHMETIC", (ftnlen)10, (
		    ftnlen)10) == 0) {
		ep = (kvdws_1.eps4[i__ - 1] + kvdws_1.eps4[k - 1]) * .5;
	    } else if (s_cmp(vdwpot_1.epsrule, "GEOMETRIC", (ftnlen)9, (
		    ftnlen)9) == 0) {
		ep = seps4[i__ - 1] * seps4[k - 1];
	    } else if (s_cmp(vdwpot_1.epsrule, "HARMONIC", (ftnlen)8, (ftnlen)
		    8) == 0) {
		ep = kvdws_1.eps4[i__ - 1] * kvdws_1.eps4[k - 1] * 2. / (
			kvdws_1.eps4[i__ - 1] + kvdws_1.eps4[k - 1]);
	    } else if (s_cmp(vdwpot_1.epsrule, "HHG", (ftnlen)3, (ftnlen)3) ==
		     0) {
/* Computing 2nd power */
		d__1 = seps4[i__ - 1] + seps4[k - 1];
		ep = kvdws_1.eps4[i__ - 1] * kvdws_1.eps4[k - 1] * 4. / (d__1 
			* d__1);
	    } else {
		ep = seps4[i__ - 1] * seps4[k - 1];
	    }
	    epsilon4_ref(i__, k) = ep;
	    epsilon4_ref(k, i__) = ep;
	}
    }

/*     vdw reduction factor information for each individual atom */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vdw_1.kred[i__ - 1] = kvdws_1.reduct[vdw_1.jvdw[i__ - 1] - 1];
	if (couple_1.n12[i__ - 1] != 1 || vdw_1.kred[i__ - 1] == 0.) {
	    vdw_1.ired[i__ - 1] = i__;
	} else {
	    vdw_1.ired[i__ - 1] = i12_ref(1, i__);
	}
    }

/*     radii and well depths for special atom class pairs */

    for (i__ = 1; i__ <= 500; ++i__) {
	if (s_cmp(kvpr_ref(0, i__), blank, (ftnlen)8, (ftnlen)8) == 0) {
	    goto L230;
	}
	ia = number_(kvpr_ref(0, i__), (ftnlen)4);
	ib = number_(kvpr_ref(4, i__), (ftnlen)4);
	if (kvdws_1.rad[ia - 1] == 0.) {
	    kvdws_1.rad[ia - 1] = .001;
	}
	if (kvdws_1.rad[ib - 1] == 0.) {
	    kvdws_1.rad[ib - 1] = .001;
	}
	if (s_cmp(vdwpot_1.radtyp, "SIGMA", (ftnlen)5, (ftnlen)5) == 0) {
	    kvdwpr_1.radpr[i__ - 1] *= 1.122462048309372981;
	}
	radmin_ref(ia, ib) = kvdwpr_1.radpr[i__ - 1];
	radmin_ref(ib, ia) = kvdwpr_1.radpr[i__ - 1];
	epsilon_ref(ia, ib) = (d__1 = kvdwpr_1.epspr[i__ - 1], abs(d__1));
	epsilon_ref(ib, ia) = (d__1 = kvdwpr_1.epspr[i__ - 1], abs(d__1));
	radmin4_ref(ia, ib) = kvdwpr_1.radpr[i__ - 1];
	radmin4_ref(ib, ia) = kvdwpr_1.radpr[i__ - 1];
	epsilon4_ref(ia, ib) = (d__1 = kvdwpr_1.epspr[i__ - 1], abs(d__1));
	epsilon4_ref(ib, ia) = (d__1 = kvdwpr_1.epspr[i__ - 1], abs(d__1));
    }
L230:

/*     radii and well depths for hydrogen bonding pairs */

    if (s_cmp(vdwpot_1.vdwtyp, "MM3-HBOND", (ftnlen)13, (ftnlen)9) == 0) {
	for (i__ = 1; i__ <= 1000; ++i__) {
	    for (k = 1; k <= 1000; ++k) {
		radhbnd_ref(k, i__) = 0.;
		epshbnd_ref(k, i__) = 0.;
	    }
	}
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(khb_ref(0, i__), blank, (ftnlen)8, (ftnlen)8) == 0) {
		goto L240;
	    }
	    ia = number_(khb_ref(0, i__), (ftnlen)4);
	    ib = number_(khb_ref(4, i__), (ftnlen)4);
	    if (kvdws_1.rad[ia - 1] == 0.) {
		kvdws_1.rad[ia - 1] = .001;
	    }
	    if (kvdws_1.rad[ib - 1] == 0.) {
		kvdws_1.rad[ib - 1] = .001;
	    }
	    if (s_cmp(vdwpot_1.radtyp, "SIGMA", (ftnlen)5, (ftnlen)5) == 0) {
		khbond_1.radhb[i__ - 1] *= 1.122462048309372981;
	    }
	    radhbnd_ref(ia, ib) = khbond_1.radhb[i__ - 1];
	    radhbnd_ref(ib, ia) = khbond_1.radhb[i__ - 1];
	    epshbnd_ref(ia, ib) = (d__1 = khbond_1.epshb[i__ - 1], abs(d__1));
	    epshbnd_ref(ib, ia) = (d__1 = khbond_1.epshb[i__ - 1], abs(d__1));
	}
L240:
	;
    }

/*     set coefficients for Gaussian fit to eps=1 and radmin=1 */

    if (s_cmp(vdwpot_1.vdwtyp, "GAUSSIAN", (ftnlen)13, (ftnlen)8) == 0) {
	if (s_cmp(vdwpot_1.gausstyp, "LJ-4", (ftnlen)8, (ftnlen)4) == 0) {
	    vdwpot_1.ngauss = 4;
	    igauss_ref(1, 1) = 846706.7;
	    igauss_ref(2, 1) = 19.483929383599527;
	    igauss_ref(1, 2) = 2713.651;
	    igauss_ref(2, 2) = 9.256482463446396;
	    igauss_ref(1, 3) = -9.699172;
	    igauss_ref(2, 3) = 2.3313232628966012;
	    igauss_ref(1, 4) = -.715442;
	    igauss_ref(2, 4) = .80587196185480869;
	} else if (s_cmp(vdwpot_1.gausstyp, "LJ-2", (ftnlen)8, (ftnlen)4) == 
		0) {
	    vdwpot_1.ngauss = 2;
	    igauss_ref(1, 1) = 14487.1;
	    igauss_ref(2, 1) = 11.404150184702447;
	    igauss_ref(1, 2) = -5.55338;
	    igauss_ref(2, 2) = 1.5438568576991818;
	} else if (s_cmp(vdwpot_1.gausstyp, "MM3-2", (ftnlen)8, (ftnlen)5) == 
		0) {
	    vdwpot_1.ngauss = 2;
	    igauss_ref(1, 1) = 2438.886;
	    igauss_ref(2, 1) = 9.342616;
	    igauss_ref(1, 2) = -6.197368;
	    igauss_ref(2, 2) = 1.564486;
	} else if (s_cmp(vdwpot_1.gausstyp, "MM2-2", (ftnlen)8, (ftnlen)5) == 
		0) {
	    vdwpot_1.ngauss = 2;
	    igauss_ref(1, 1) = 3423.562;
	    igauss_ref(2, 1) = 9.692821;
	    igauss_ref(1, 2) = -6.50376;
	    igauss_ref(2, 2) = 1.585344;
	} else if (s_cmp(vdwpot_1.gausstyp, "IN-PLACE", (ftnlen)8, (ftnlen)8) 
		== 0) {
	    vdwpot_1.ngauss = 2;
	    igauss_ref(1, 1) = 500.;
	    igauss_ref(2, 1) = 6.143;
	    igauss_ref(1, 2) = -18.831;
	    igauss_ref(2, 2) = 2.209;
	}
    }

/*     remove zero-sized atoms from the list of vdw sites */

    vdw_1.nvdw = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (kvdws_1.rad[vdw_1.jvdw[i__ - 1] - 1] != 0.) {
	    ++vdw_1.nvdw;
	    vdw_1.ivdw[vdw_1.nvdw - 1] = i__;
	}
    }

/*     turn off the van der Waals potential if it is not used */

    if (vdw_1.nvdw == 0) {
	potent_1.use_vdw__ = FALSE_;
    }
    return 0;
} /* kvdw_ */

#undef epsilon_ref
#undef keyline_ref
#undef epshbnd_ref
#undef radhbnd_ref
#undef radmin4_ref
#undef igauss_ref
#undef radmin_ref
#undef kvpr_ref
#undef khb_ref
#undef i12_ref
#undef epsilon4_ref


