/* testpair.f -- translated by f2c (version 20050501).
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
    doublereal vdwcut, chgcut, dplcut, mpolecut, vdwtaper, chgtaper, dpltaper,
	     mpoletaper, ewaldcut;
    logical use_ewald__, use_lights__, use_list__, use_vlist__, use_clist__, 
	    use_mlist__;
} cutoff_;

#define cutoff_1 cutoff_

struct {
    doublereal desum[75000]	/* was [3][25000] */, deb[75000]	/* 
	    was [3][25000] */, dea[75000]	/* was [3][25000] */, deba[
	    75000]	/* was [3][25000] */, deub[75000]	/* was [3][
	    25000] */, deaa[75000]	/* was [3][25000] */, deopb[75000]	
	    /* was [3][25000] */, deopd[75000]	/* was [3][25000] */, deid[
	    75000]	/* was [3][25000] */, deit[75000]	/* was [3][
	    25000] */, det[75000]	/* was [3][25000] */, dept[75000]	
	    /* was [3][25000] */, debt[75000]	/* was [3][25000] */, dett[
	    75000]	/* was [3][25000] */, dev[75000]	/* was [3][
	    25000] */, dec[75000]	/* was [3][25000] */, decd[75000]	
	    /* was [3][25000] */, ded[75000]	/* was [3][25000] */, dem[
	    75000]	/* was [3][25000] */, dep[75000]	/* was [3][
	    25000] */, der[75000]	/* was [3][25000] */, des[75000]	
	    /* was [3][25000] */, delf[75000]	/* was [3][25000] */, deg[
	    75000]	/* was [3][25000] */, dex[75000]	/* was [3][
	    25000] */;
} deriv_;

#define deriv_1 deriv_

struct {
    doublereal esum, eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, 
	    ebt, ett, ev, ec, ecd, ed, em, ep, er, es, elf, eg, ex;
} energi_;

#define energi_1 energi_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nlight, kbx[25000], kby[25000], kbz[25000], kex[25000], key[25000]
	    , kez[25000], locx[200000], locy[200000], locz[200000], rgx[
	    200000], rgy[200000], rgz[200000];
} light_;

#define light_1 light_

struct {
    doublereal lbuffer, lbuf2, vbuf2, cbuf2, mbuf2;
    integer nvlst[25000], vlst[45000000]	/* was [1800][25000] */, 
	    nelst[25000], elst[30000000]	/* was [1200][25000] */;
    logical dovlst, doclst, domlst;
} neigh_;

#define neigh_1 neigh_

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
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  program testpair  --  time various neighbor pair schemes  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "testpair" performs a set of timing tests to compare the */
/*     evaluation of potential energy and energy/gradient using */
/*     different methods for finding pairwise neighbors */


/* Main program */ int MAIN__(void)
{
    /* Initialized data */

    static char axis[1*3] = "X" "Y" "Z";

    /* Format strings */
    static char fmt_20[] = "(/,\002 Enter Desired Number of Repetitions [1] "
	    ":  \002,$)";
    static char fmt_30[] = "(i10)";
    static char fmt_40[] = "(/,\002 Computation Overhead :\002,11x,\002Tim"
	    "e\002,7x,i5,\002 Evaluations\002)";
    static char fmt_50[] = "(/,\002 Double Nested Loop\002,9x,f10.3)";
    static char fmt_80[] = "(\002 Method of Lights\002,11x,f10.3)";
    static char fmt_90[] = "(\002 Pair Neighbor List\002,9x,f10.3)";
    static char fmt_100[] = "(/,\002 Potential Energy Only :\002,10x,\002T"
	    "ime\002,14x,\002Evdw\002,12x,\002Eelect\002)";
    static char fmt_110[] = "(/,\002 Double Nested Loop\002,9x,f10.3,3x,f15."
	    "4,3x,f15.4)";
    static char fmt_120[] = "(\002 Method of Lights\002,11x,f10.3,3x,f15.4,3"
	    "x,f15.4)";
    static char fmt_130[] = "(\002 Pair Neighbor List\002,9x,f10.3,3x,f15.4,"
	    "3x,f15.4)";
    static char fmt_140[] = "(/,\002 Energy and Gradient :\002,12x,\002Tim"
	    "e\002,14x,\002Dvdw\002,12x,\002Delect\002)";
    static char fmt_150[] = "(/,\002 Double Nested Loop\002,9x,f10.3,3x,f15."
	    "4,3x,f15.4)";
    static char fmt_160[] = "(\002 Method of Lights\002,11x,f10.3,3x,f15.4,3"
	    "x,f15.4)";
    static char fmt_170[] = "(\002 Pair Neighbor List\002,9x,f10.3,3x,f15.4,"
	    "3x,f15.4)";
    static char fmt_180[] = "(/,\002 Comparison of Nonbond Gradients from"
	    "\002,\002 the Neighbor Methods:\002,//,11x,\002Component\002,16x,"
	    "\002Loop\002,10x,\002Lights\002,12x,\002List\002,/)";
    static char fmt_190[] = "(10x,i6,\002 (\002,a1,\002)\002,4x,3f16.4)";
    static char fmt_200[] = "(/,\002 Gradients Computed via all the Neighbor"
	    " Methods\002,\002 are Identical\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void), s_cmp(char *, 
	    char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    static integer i__, j, k, m;
    static doublereal r2, xi, yi, zi, xr, yr, zr, off;
    extern /* Subroutine */ int elj_(void);
    static integer kgy, kgz;
    static doublereal eps;
    extern /* Subroutine */ int elj1_(void);
    static doublereal off2;
    extern /* Subroutine */ int ehal_(void);
    static doublereal erms;
    static integer stop;
    static doublereal vrms;
    extern /* Subroutine */ int ehal1_(void), image_(doublereal *, doublereal 
	    *, doublereal *), fatal_(void), ebuck_(void);
    static logical match;
    extern /* Subroutine */ int final_(void);
    static integer npair;
    static doublereal dloop[75000]	/* was [3][25000] */, dlist[75000]	
	    /* was [3][25000] */;
    static integer nterm;
    static logical exist;
    static integer start;
    static logical query;
    extern /* Subroutine */ int ebuck1_(void), emm3hb_(void);
    static doublereal xsort[200000], ysort[200000], zsort[200000];
    static logical header;
    extern /* Subroutine */ int emm3hb1_(void), getime_(doublereal *);
    static integer ncalls;
    static doublereal dlight[75000]	/* was [3][25000] */;
    static logical repeat;
    extern /* Subroutine */ int empole_(void);
    static char string[120];
    extern /* Subroutine */ int setime_(void), lights_(doublereal *, integer *
	    , doublereal *, doublereal *, doublereal *), nblist_(void), 
	    egauss_(void), getxyz_(void), empole1_(void), egauss1_(void), 
	    echarge_(void);
    static doublereal elapsed;
    extern /* Subroutine */ int initial_(void), nextarg_(char *, logical *, 
	    ftnlen), echarge1_(void);

    /* Fortran I/O blocks */
    static icilist io___6 = { 1, string, 1, 0, 120, 1 };
    static cilist io___7 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___8 = { 1, 0, 0, fmt_30, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___49 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___54 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___55 = { 0, 0, 0, fmt_200, 0 };



#define dec_ref(a_1,a_2) deriv_1.dec[(a_2)*3 + a_1 - 4]
#define dem_ref(a_1,a_2) deriv_1.dem[(a_2)*3 + a_1 - 4]
#define dep_ref(a_1,a_2) deriv_1.dep[(a_2)*3 + a_1 - 4]
#define dev_ref(a_1,a_2) deriv_1.dev[(a_2)*3 + a_1 - 4]
#define dloop_ref(a_1,a_2) dloop[(a_2)*3 + a_1 - 4]
#define dlist_ref(a_1,a_2) dlist[(a_2)*3 + a_1 - 4]
#define dlight_ref(a_1,a_2) dlight[(a_2)*3 + a_1 - 4]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cutoff.i  --  cutoff distances for energy interactions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     vdwcut      cutoff distance for van der Waals interactions */
/*     chgcut      cutoff distance for charge-charge interactions */
/*     dplcut      cutoff distance for dipole-dipole interactions */
/*     mpolecut    cutoff distance for atomic multipole interactions */
/*     vdwtaper    distance at which van der Waals switching begins */
/*     chgtaper    distance at which charge-charge switching begins */
/*     dpltaper    distance at which dipole-dipole switching begins */
/*     mpoletaper  distance at which atomic multipole switching begins */
/*     ewaldcut    cutoff distance for direct space Ewald summation */
/*     use_ewald   logical flag governing use of Ewald summation */
/*     use_lights  logical flag governing use of method of lights */
/*     use_list    logical flag governing use of any neighbor lists */
/*     use_vlist   logical flag governing use of vdw neighbor lists */
/*     use_clist   logical flag governing use of charge neighbor lists */
/*     use_mlist   logical flag governing use of multipole neighbor lists */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  deriv.i  --  Cartesian coordinate derivative components  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     desum   total energy Cartesian coordinate derivatives */
/*     deb     bond stretch Cartesian coordinate derivatives */
/*     dea     angle bend Cartesian coordinate derivatives */
/*     deba    stretch-bend Cartesian coordinate derivatives */
/*     deub    Urey-Bradley Cartesian coordinate derivatives */
/*     deaa    angle-angle Cartesian coordinate derivatives */
/*     deopb   out-of-plane bend Cartesian coordinate derivatives */
/*     deopd   out-of-plane distance Cartesian coordinate derivatives */
/*     deid    improper dihedral Cartesian coordinate derivatives */
/*     deit    improper torsion Cartesian coordinate derivatives */
/*     det     torsional Cartesian coordinate derivatives */
/*     dept    pi-orbital torsion Cartesian coordinate derivatives */
/*     debt    stretch-torsion Cartesian coordinate derivatives */
/*     dett    torsion-torsion Cartesian coordinate derivatives */
/*     dev     van der Waals Cartesian coordinate derivatives */
/*     dec     charge-charge Cartesian coordinate derivatives */
/*     decd    charge-dipole Cartesian coordinate derivatives */
/*     ded     dipole-dipole Cartesian coordinate derivatives */
/*     dem     multipole Cartesian coordinate derivatives */
/*     dep     polarization Cartesian coordinate derivatives */
/*     der     reaction field Cartesian coordinate derivatives */
/*     des     solvation Cartesian coordinate derivatives */
/*     delf    metal ligand field Cartesian coordinate derivatives */
/*     deg     geometric restraint Cartesian coordinate derivatives */
/*     dex     extra energy term Cartesian coordinate derivatives */




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
/*     ##  light.i  --  indices for method of lights pair neighbors  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     nlight  total number of sites for method of lights calculation */
/*     kbx     low index of neighbors of each site in the x-sorted list */
/*     kby     low index of neighbors of each site in the y-sorted list */
/*     kbz     low index of neighbors of each site in the z-sorted list */
/*     kex     high index of neighbors of each site in the x-sorted list */
/*     key     high index of neighbors of each site in the y-sorted list */
/*     kez     high index of neighbors of each site in the z-sorted list */
/*     locx    pointer from x-sorted list into original interaction list */
/*     locy    pointer from y-sorted list into original interaction list */
/*     locz    pointer from z-sorted list into original interaction list */
/*     rgx     pointer from original interaction list into x-sorted list */
/*     rgy     pointer from original interaction list into y-sorted list */
/*     rgz     pointer from original interaction list into z-sorted list */




/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  neigh.i  --  pairwise neighbor list indices and storage  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     lbuffer     width of the neighbor list buffer region */
/*     lbuf2       motion squared needed to trigger list rebuild */
/*     vbuf2       square of vdw cutoff plus neighbor list buffer */
/*     cbuf2       square of charge cutoff plus neighbor list buffer */
/*     mbuf2       square of multipole cutoff plus neighbor list buffer */
/*     nvlst       number of sites in list for each vdw site */
/*     vlst        site numbers in neighbor list of each vdw site */
/*     nelst       number of sites in list for each electrostatic site */
/*     elst        site numbers in list of each electrostatic site */
/*     dovlst      logical flag to rebuild vdw neighbor list */
/*     doclst      logical flag to rebuild charge neighbor list */
/*     domlst      logical flag to rebuild multipole neighbor list */




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




/*     read the molecular system and setup molecular mechanics */

    initial_();
    getxyz_();
    mechanic_();

/*     get the number calculation cycles to be performed */

    ncalls = 0;
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___6);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	query = FALSE_;
    }
L10:
    if (query) {
	io___7.ciunit = iounit_1.iout;
	s_wsfe(&io___7);
	e_wsfe();
	io___8.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___8);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L10;
	}
    }
    if (ncalls == 0) {
	ncalls = 1;
    }

/*     initialize number of pairs and generic cutoff distance */

    npair = 0;
    nterm = 0;
    if (potent_1.use_vdw__) {
	++nterm;
    }
    if (potent_1.use_charge__) {
	++nterm;
    }
    if (potent_1.use_mpole__ || potent_1.use_polar__) {
	++nterm;
    }
    nterm *= ncalls;
    off = 5.;
    off2 = off * off;

/*     get the timing for setup of double nested loop */

    setime_();
    i__1 = nterm;
    for (m = 1; m <= i__1; ++m) {
	i__2 = atoms_1.n - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xi = atoms_1.x[i__ - 1];
	    yi = atoms_1.y[i__ - 1];
	    zi = atoms_1.z__[i__ - 1];
	    i__3 = atoms_1.n;
	    for (j = i__ + 1; j <= i__3; ++j) {
		xr = atoms_1.x[j - 1] - xi;
		yr = atoms_1.y[j - 1] - yi;
		zr = atoms_1.z__[j - 1] - zi;
		image_(&xr, &yr, &zr);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 < off2) {
		    ++npair;
		}
	    }
	}
    }
    getime_(&elapsed);
    io___24.ciunit = iounit_1.iout;
    s_wsfe(&io___24);
    do_fio(&c__1, (char *)&ncalls, (ftnlen)sizeof(integer));
    e_wsfe();
    io___25.ciunit = iounit_1.iout;
    s_wsfe(&io___25);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     get the timing for setup of method of lights */

    setime_();
    i__1 = nterm;
    for (m = 1; m <= i__1; ++m) {
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xsort[i__ - 1] = atoms_1.x[i__ - 1];
	    ysort[i__ - 1] = atoms_1.y[i__ - 1];
	    zsort[i__ - 1] = atoms_1.z__[i__ - 1];
	}
	lights_(&off, &atoms_1.n, xsort, ysort, zsort);
	i__2 = atoms_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xi = xsort[light_1.rgx[i__ - 1] - 1];
	    yi = ysort[light_1.rgy[i__ - 1] - 1];
	    zi = zsort[light_1.rgz[i__ - 1] - 1];
	    if (light_1.kbx[i__ - 1] <= light_1.kex[i__ - 1]) {
		repeat = FALSE_;
		start = light_1.kbx[i__ - 1] + 1;
		stop = light_1.kex[i__ - 1];
	    } else {
		repeat = TRUE_;
		start = 1;
		stop = light_1.kex[i__ - 1];
	    }
L60:
	    i__3 = stop;
	    for (j = start; j <= i__3; ++j) {
		k = light_1.locx[j - 1];
		kgy = light_1.rgy[k - 1];
		if (light_1.kby[i__ - 1] <= light_1.key[i__ - 1]) {
		    if (kgy < light_1.kby[i__ - 1] || kgy > light_1.key[i__ - 
			    1]) {
			goto L70;
		    }
		} else {
		    if (kgy < light_1.kby[i__ - 1] && kgy > light_1.key[i__ - 
			    1]) {
			goto L70;
		    }
		}
		kgz = light_1.rgz[k - 1];
		if (light_1.kbz[i__ - 1] <= light_1.kez[i__ - 1]) {
		    if (kgz < light_1.kbz[i__ - 1] || kgz > light_1.kez[i__ - 
			    1]) {
			goto L70;
		    }
		} else {
		    if (kgz < light_1.kbz[i__ - 1] && kgz > light_1.kez[i__ - 
			    1]) {
			goto L70;
		    }
		}
		xr = xi - xsort[j - 1];
		yr = yi - ysort[kgy - 1];
		zr = zi - zsort[kgz - 1];
		image_(&xr, &yr, &zr);
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 < off2) {
		    ++npair;
		}
L70:
		;
	    }
	    if (repeat) {
		repeat = FALSE_;
		start = light_1.kbx[i__ - 1] + 1;
		stop = light_1.nlight;
		goto L60;
	    }
	}
    }
    getime_(&elapsed);
    io___35.ciunit = iounit_1.iout;
    s_wsfe(&io___35);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    e_wsfe();
    if (npair < 0) {
	fatal_();
    }

/*     get the timing for setup of pair neighbor list */

    setime_();
    i__1 = ncalls;
    for (m = 1; m <= i__1; ++m) {
	neigh_1.dovlst = TRUE_;
	neigh_1.doclst = TRUE_;
	neigh_1.domlst = TRUE_;
	nblist_();
    }
    getime_(&elapsed);
    io___36.ciunit = iounit_1.iout;
    s_wsfe(&io___36);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     get the timing for energy terms via double nested loop */

    cutoff_1.use_lights__ = FALSE_;
    cutoff_1.use_vlist__ = FALSE_;
    cutoff_1.use_clist__ = FALSE_;
    cutoff_1.use_mlist__ = FALSE_;
    setime_();
    i__1 = ncalls;
    for (k = 1; k <= i__1; ++k) {
	if (potent_1.use_vdw__) {
	    if (s_cmp(vdwpot_1.vdwtyp, "LENNARD-JONES", (ftnlen)13, (ftnlen)
		    13) == 0) {
		elj_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "BUCKINGHAM", (ftnlen)13, (ftnlen)10) 
		    == 0) {
		ebuck_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "MM3-HBOND", (ftnlen)13, (ftnlen)9) == 
		    0) {
		emm3hb_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "BUFFERED-14-7", (ftnlen)13, (ftnlen)
		    13) == 0) {
		ehal_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "GAUSSIAN", (ftnlen)13, (ftnlen)8) == 
		    0) {
		egauss_();
	    }
	}
	if (potent_1.use_charge__) {
	    echarge_();
	}
	if (potent_1.use_mpole__ || potent_1.use_polar__) {
	    empole_();
	}
    }
    getime_(&elapsed);
    io___37.ciunit = iounit_1.iout;
    s_wsfe(&io___37);
    e_wsfe();
    io___38.ciunit = iounit_1.iout;
    s_wsfe(&io___38);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&energi_1.ev, (ftnlen)sizeof(doublereal));
    d__1 = energi_1.ec + energi_1.em + energi_1.ep;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     get the timing for energy terms via method of lights */

    cutoff_1.use_lights__ = TRUE_;
    cutoff_1.use_vlist__ = FALSE_;
    cutoff_1.use_clist__ = FALSE_;
    cutoff_1.use_mlist__ = FALSE_;
    setime_();
    i__1 = ncalls;
    for (k = 1; k <= i__1; ++k) {
	if (potent_1.use_vdw__) {
	    if (s_cmp(vdwpot_1.vdwtyp, "LENNARD-JONES", (ftnlen)13, (ftnlen)
		    13) == 0) {
		elj_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "BUCKINGHAM", (ftnlen)13, (ftnlen)10) 
		    == 0) {
		ebuck_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "MM3-HBOND", (ftnlen)13, (ftnlen)9) == 
		    0) {
		emm3hb_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "BUFFERED-14-7", (ftnlen)13, (ftnlen)
		    13) == 0) {
		ehal_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "GAUSSIAN", (ftnlen)13, (ftnlen)8) == 
		    0) {
		egauss_();
	    }
	}
	if (potent_1.use_charge__) {
	    echarge_();
	}
	if (potent_1.use_mpole__ || potent_1.use_polar__) {
	    empole_();
	}
    }
    getime_(&elapsed);
    io___39.ciunit = iounit_1.iout;
    s_wsfe(&io___39);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&energi_1.ev, (ftnlen)sizeof(doublereal));
    d__1 = energi_1.ec + energi_1.em + energi_1.ep;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     get the timing for energy terms via pair neighbor list */

    cutoff_1.use_lights__ = FALSE_;
    cutoff_1.use_vlist__ = TRUE_;
    cutoff_1.use_clist__ = TRUE_;
    cutoff_1.use_mlist__ = TRUE_;
    nblist_();
    setime_();
    i__1 = ncalls;
    for (k = 1; k <= i__1; ++k) {
	if (potent_1.use_vdw__) {
	    if (s_cmp(vdwpot_1.vdwtyp, "LENNARD-JONES", (ftnlen)13, (ftnlen)
		    13) == 0) {
		elj_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "BUCKINGHAM", (ftnlen)13, (ftnlen)10) 
		    == 0) {
		ebuck_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "MM3-HBOND", (ftnlen)13, (ftnlen)9) == 
		    0) {
		emm3hb_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "BUFFERED-14-7", (ftnlen)13, (ftnlen)
		    13) == 0) {
		ehal_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "GAUSSIAN", (ftnlen)13, (ftnlen)8) == 
		    0) {
		egauss_();
	    }
	}
	if (potent_1.use_charge__) {
	    echarge_();
	}
	if (potent_1.use_mpole__ || potent_1.use_polar__) {
	    empole_();
	}
    }
    getime_(&elapsed);
    io___40.ciunit = iounit_1.iout;
    s_wsfe(&io___40);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&energi_1.ev, (ftnlen)sizeof(doublereal));
    d__1 = energi_1.ec + energi_1.em + energi_1.ep;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     get the timing for gradient via double nested loop */

    cutoff_1.use_lights__ = FALSE_;
    cutoff_1.use_vlist__ = FALSE_;
    cutoff_1.use_clist__ = FALSE_;
    cutoff_1.use_mlist__ = FALSE_;
    setime_();
    i__1 = ncalls;
    for (k = 1; k <= i__1; ++k) {
	if (potent_1.use_vdw__) {
	    if (s_cmp(vdwpot_1.vdwtyp, "LENNARD-JONES", (ftnlen)13, (ftnlen)
		    13) == 0) {
		elj1_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "BUCKINGHAM", (ftnlen)13, (ftnlen)10) 
		    == 0) {
		ebuck1_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "MM3-HBOND", (ftnlen)13, (ftnlen)9) == 
		    0) {
		emm3hb1_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "BUFFERED-14-7", (ftnlen)13, (ftnlen)
		    13) == 0) {
		ehal1_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "GAUSSIAN", (ftnlen)13, (ftnlen)8) == 
		    0) {
		egauss1_();
	    }
	}
	if (potent_1.use_charge__) {
	    echarge1_();
	}
	if (potent_1.use_mpole__ || potent_1.use_polar__) {
	    empole1_();
	}
    }
    getime_(&elapsed);

/*     store the double loop gradient and get rms values */

    vrms = 0.;
    erms = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    dloop_ref(j, i__) = dev_ref(j, i__) + dec_ref(j, i__) + dem_ref(j,
		     i__) + dep_ref(j, i__);
/* Computing 2nd power */
	    d__1 = dev_ref(j, i__);
	    vrms += d__1 * d__1;
/* Computing 2nd power */
	    d__1 = dec_ref(j, i__);
/* Computing 2nd power */
	    d__2 = dem_ref(j, i__);
/* Computing 2nd power */
	    d__3 = dep_ref(j, i__);
	    erms = erms + d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	}
    }
    vrms = sqrt(vrms / (doublereal) atoms_1.n);
    erms = sqrt(erms / (doublereal) atoms_1.n);
    io___44.ciunit = iounit_1.iout;
    s_wsfe(&io___44);
    e_wsfe();
    io___45.ciunit = iounit_1.iout;
    s_wsfe(&io___45);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&vrms, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&erms, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     get the timing for gradient via method of lights */

    cutoff_1.use_lights__ = TRUE_;
    cutoff_1.use_vlist__ = FALSE_;
    cutoff_1.use_clist__ = FALSE_;
    cutoff_1.use_mlist__ = FALSE_;
    setime_();
    i__1 = ncalls;
    for (k = 1; k <= i__1; ++k) {
	if (potent_1.use_vdw__) {
	    if (s_cmp(vdwpot_1.vdwtyp, "LENNARD-JONES", (ftnlen)13, (ftnlen)
		    13) == 0) {
		elj1_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "BUCKINGHAM", (ftnlen)13, (ftnlen)10) 
		    == 0) {
		ebuck1_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "MM3-HBOND", (ftnlen)13, (ftnlen)9) == 
		    0) {
		emm3hb1_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "BUFFERED-14-7", (ftnlen)13, (ftnlen)
		    13) == 0) {
		ehal1_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "GAUSSIAN", (ftnlen)13, (ftnlen)8) == 
		    0) {
		egauss1_();
	    }
	}
	if (potent_1.use_charge__) {
	    echarge1_();
	}
	if (potent_1.use_mpole__ || potent_1.use_polar__) {
	    empole1_();
	}
    }
    getime_(&elapsed);

/*     store the method of lights gradient and get rms values */

    vrms = 0.;
    erms = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    dlight_ref(j, i__) = dev_ref(j, i__) + dec_ref(j, i__) + dem_ref(
		    j, i__) + dep_ref(j, i__);
/* Computing 2nd power */
	    d__1 = dev_ref(j, i__);
	    vrms += d__1 * d__1;
/* Computing 2nd power */
	    d__1 = dec_ref(j, i__);
/* Computing 2nd power */
	    d__2 = dem_ref(j, i__);
/* Computing 2nd power */
	    d__3 = dep_ref(j, i__);
	    erms = erms + d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	}
    }
    vrms = sqrt(vrms / (doublereal) atoms_1.n);
    erms = sqrt(erms / (doublereal) atoms_1.n);
    io___47.ciunit = iounit_1.iout;
    s_wsfe(&io___47);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&vrms, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&erms, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     get the timing for gradient via pair neighbor list */

    cutoff_1.use_lights__ = FALSE_;
    cutoff_1.use_vlist__ = TRUE_;
    cutoff_1.use_clist__ = TRUE_;
    cutoff_1.use_mlist__ = TRUE_;
    setime_();
    i__1 = ncalls;
    for (k = 1; k <= i__1; ++k) {
	if (potent_1.use_vdw__) {
	    if (s_cmp(vdwpot_1.vdwtyp, "LENNARD-JONES", (ftnlen)13, (ftnlen)
		    13) == 0) {
		elj1_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "BUCKINGHAM", (ftnlen)13, (ftnlen)10) 
		    == 0) {
		ebuck1_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "MM3-HBOND", (ftnlen)13, (ftnlen)9) == 
		    0) {
		emm3hb1_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "BUFFERED-14-7", (ftnlen)13, (ftnlen)
		    13) == 0) {
		ehal1_();
	    }
	    if (s_cmp(vdwpot_1.vdwtyp, "GAUSSIAN", (ftnlen)13, (ftnlen)8) == 
		    0) {
		egauss1_();
	    }
	}
	if (potent_1.use_charge__) {
	    echarge1_();
	}
	if (potent_1.use_mpole__ || potent_1.use_polar__) {
	    empole1_();
	}
    }
    getime_(&elapsed);

/*     get the pair neighbor list gradient rms values */

    vrms = 0.;
    erms = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    dlist_ref(j, i__) = dev_ref(j, i__) + dec_ref(j, i__) + dem_ref(j,
		     i__) + dep_ref(j, i__);
/* Computing 2nd power */
	    d__1 = dev_ref(j, i__);
	    vrms += d__1 * d__1;
/* Computing 2nd power */
	    d__1 = dec_ref(j, i__);
/* Computing 2nd power */
	    d__2 = dem_ref(j, i__);
/* Computing 2nd power */
	    d__3 = dep_ref(j, i__);
	    erms = erms + d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	}
    }
    vrms = sqrt(vrms / (doublereal) atoms_1.n);
    erms = sqrt(erms / (doublereal) atoms_1.n);
    io___49.ciunit = iounit_1.iout;
    s_wsfe(&io___49);
    do_fio(&c__1, (char *)&elapsed, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&vrms, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&erms, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*     compare the nonbond gradients from the various methods */

    eps = 1e-4;
    match = TRUE_;
    header = TRUE_;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    if ((d__1 = dlight_ref(j, i__) - dloop_ref(j, i__), abs(d__1)) > 
		    eps || (d__2 = dlist_ref(j, i__) - dloop_ref(j, i__), abs(
		    d__2)) > eps) {
		if (header) {
		    match = FALSE_;
		    header = FALSE_;
		    io___53.ciunit = iounit_1.iout;
		    s_wsfe(&io___53);
		    e_wsfe();
		}
		io___54.ciunit = iounit_1.iout;
		s_wsfe(&io___54);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, axis + (j - 1), (ftnlen)1);
		do_fio(&c__1, (char *)&dloop_ref(j, i__), (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&dlight_ref(j, i__), (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&dlist_ref(j, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     success if all the neighbor methods give the same gradient */

    if (match) {
	io___55.ciunit = iounit_1.iout;
	s_wsfe(&io___55);
	e_wsfe();
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef dlight_ref
#undef dlist_ref
#undef dloop_ref
#undef dev_ref
#undef dep_ref
#undef dem_ref
#undef dec_ref


/* Main program alias */ int testpair_ () { MAIN__ (); return 0; }
