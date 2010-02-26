/* hessian.f -- translated by f2c (version 20050501).
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
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal vdwcut, chgcut, dplcut, mpolecut, vdwtaper, chgtaper, dpltaper,
	     mpoletaper, ewaldcut;
    logical use_ewald__, use_lights__, use_list__, use_vlist__, use_clist__, 
	    use_mlist__;
} cutoff_;

#define cutoff_1 cutoff_

struct {
    doublereal hesscut;
} hescut_;

#define hescut_1 hescut_

struct {
    doublereal hessx[75000]	/* was [3][25000] */, hessy[75000]	/* 
	    was [3][25000] */, hessz[75000]	/* was [3][25000] */;
} hessn_;

#define hessn_1 hessn_

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
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

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
    doublereal xrb[25000], yrb[25000], zrb[25000], rbc[6000]	/* was [6][
	    1000] */;
    logical use_rigid__;
} rigid_;

#define rigid_1 rigid_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal radmin[1000000]	/* was [1000][1000] */, epsilon[1000000]	
	    /* was [1000][1000] */, radmin4[1000000]	/* was [1000][1000] */
	    , epsilon4[1000000]	/* was [1000][1000] */, radhbnd[1000000]	
	    /* was [1000][1000] */, epshbnd[1000000]	/* was [1000][1000] */
	    , kred[25000];
    integer ired[25000], nvdw, ivdw[25000], jvdw[25000];
} vdw_;

#define vdw_1 vdw_

struct {
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_

/* Table of constant values */

static integer c__1 = 1;
static integer c_b12 = 1000000;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine hessian  --  atom-by-atom Hessian elements  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "hessian" calls subroutines to calculate the Hessian elements */
/*     for each atom in turn with respect to Cartesian coordinates */


/* Subroutine */ int hessian_(doublereal *h__, integer *hinit, integer *hstop,
	 integer *hindex, doublereal *hdiag)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Current Atom :\002,i7,6x,\002Total Atoms"
	    " :\002,i7,/,\002 Current Required Hessian Storage : \002,i12,/"
	    ",\002 Maximum Allowed Hessian Storage :  \002,i12,/,\002 Minimum"
	    " Significant Hessian Value :\002,f12.5)";
    static char fmt_20[] = "(/,\002 HESSIAN  --  Increase MAXHESS\002,\002 a"
	    "nd/or HESSCUT\002)";
    static char fmt_30[] = "(\002 HESSIAN  --\002,i11,\002 Elements\002,f9"
	    ".2,\002 % Off-diag H\002,f8.2,\002 % Storage\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), do_fio(
	    integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int estrbnd2_(integer *), erxnfld2_(integer *), 
	    eopdist2_(integer *), eimprop2_(integer *), eimptor2_(integer *), 
	    epitors2_(integer *), etortor2_(integer *), estrtor2_(integer *);
    static integer i__, j, k, ii;
    static doublereal rdn;
    extern /* Subroutine */ int elj2_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    static logical keep[25000];
    static doublereal hmax;
    extern /* Subroutine */ int born_(void);
    static doublereal xred[25000], yred[25000], zred[25000];
    extern /* Subroutine */ int ehal2_(integer *, doublereal *, doublereal *, 
	    doublereal *), fatal_(void), piscf_(void);
    static integer nhess;
    extern /* Subroutine */ int ebond2_(integer *), ebuck2_(integer *, 
	    doublereal *, doublereal *, doublereal *), egeom2_(integer *), 
	    extra2_(integer *), esolv2_(integer *), eurey2_(integer *), 
	    etors2_(integer *), emm3hb2_(integer *, doublereal *, doublereal *
	    , doublereal *);
    static doublereal filled;
    extern /* Subroutine */ int induce_(void);
    static doublereal cutoff;
    extern /* Subroutine */ int bounds_(void), nblist_(void), eangle2_(
	    integer *), emetal2_(integer *), empole2_(integer *), egauss2_(
	    integer *, doublereal *, doublereal *, doublereal *), replica_(
	    doublereal *), chkpole_(void);
    static doublereal percent;
    extern /* Subroutine */ int echarge2_(integer *), eangang2_(integer *), 
	    rotpole_(void), echgdpl2_(integer *), eopbend2_(integer *), 
	    edipole2_(integer *);

    /* Fortran I/O blocks */
    static cilist io___13 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_30, 0 };



#define hdiag_ref(a_1,a_2) hdiag[(a_2)*3 + a_1]
#define hinit_ref(a_1,a_2) hinit[(a_2)*3 + a_1]
#define hessx_ref(a_1,a_2) hessn_1.hessx[(a_2)*3 + a_1 - 4]
#define hessy_ref(a_1,a_2) hessn_1.hessy[(a_2)*3 + a_1 - 4]
#define hessz_ref(a_1,a_2) hessn_1.hessz[(a_2)*3 + a_1 - 4]
#define hstop_ref(a_1,a_2) hstop[(a_2)*3 + a_1]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  bound.i  --  control of periodic boundary conditions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     polycut       cutoff distance for infinite polymer nonbonds */
/*     polycut2      square of infinite polymer nonbond cutoff */
/*     use_bounds    flag to use periodic boundary conditions */
/*     use_replica   flag to use replicates for periodic system */
/*     use_polymer   flag to mark presence of infinite polymer */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  hescut.i  --  cutoff value for Hessian matrix elements  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     hesscut   magnitude of smallest allowed Hessian element */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  hessn.i  --  Cartesian Hessian elements for a single atom  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     hessx   Hessian elements for x-component of current atom */
/*     hessy   Hessian elements for y-component of current atom */
/*     hessz   Hessian elements for z-component of current atom */




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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  rigid.i  --  rigid body coordinates for atom groups  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     xrb         rigid body reference x-coordinate for each atom */
/*     yrb         rigid body reference y-coordinate for each atom */
/*     zrb         rigid body reference z-coordinate for each atom */
/*     rbc         current rigid body coordinates for each group */
/*     use_rigid   flag to mark use of rigid body coordinate system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




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




/*     zero out total number of indexed Hessian elements */

    /* Parameter adjustments */
    hdiag -= 4;
    --hindex;
    hstop -= 4;
    hinit -= 4;
    --h__;

    /* Function Body */
    nhess = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    hinit_ref(j, i__) = 1;
	    hstop_ref(j, i__) = 0;
	    hdiag_ref(j, i__) = 0.;
	}
    }

/*     maintain any periodic boundary conditions */

    if (bound_1.use_bounds__ && ! rigid_1.use_rigid__) {
	bounds_();
    }

/*     update the pairwise interaction neighbor lists */

    if (cutoff_1.use_list__) {
	nblist_();
    }

/*     many implicit solvation models require Born radii */

    if (potent_1.use_born__) {
	born_();
    }

/*     alter bond and torsion constants for pisystem */

    if (potent_1.use_orbit__) {
	piscf_();
    }

/*     compute the induced dipoles at polarizable atoms */

    if (potent_1.use_polar__) {
	chkpole_();
	rotpole_();
	induce_();
    }

/*     calculate the "reduced" atomic coordinates */

    if (potent_1.use_vdw__) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ii = vdw_1.ired[i__ - 1];
	    rdn = vdw_1.kred[i__ - 1];
	    xred[i__ - 1] = rdn * (atoms_1.x[i__ - 1] - atoms_1.x[ii - 1]) + 
		    atoms_1.x[ii - 1];
	    yred[i__ - 1] = rdn * (atoms_1.y[i__ - 1] - atoms_1.y[ii - 1]) + 
		    atoms_1.y[ii - 1];
	    zred[i__ - 1] = rdn * (atoms_1.z__[i__ - 1] - atoms_1.z__[ii - 1])
		     + atoms_1.z__[ii - 1];
	}
    }

/*     zero out the Hessian elements for the current atom */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    i__2 = atoms_1.n;
	    for (k = 1; k <= i__2; ++k) {
		for (j = 1; j <= 3; ++j) {
		    hessx_ref(j, k) = 0.;
		    hessy_ref(j, k) = 0.;
		    hessz_ref(j, k) = 0.;
		}
	    }

/*     remove any previous use of the replicates method */

	    cutoff = 0.;
	    replica_(&cutoff);

/*     call the local geometry Hessian component routines */

	    if (potent_1.use_bond__) {
		ebond2_(&i__);
	    }
	    if (potent_1.use_angle__) {
		eangle2_(&i__);
	    }
	    if (potent_1.use_strbnd__) {
		estrbnd2_(&i__);
	    }
	    if (potent_1.use_urey__) {
		eurey2_(&i__);
	    }
	    if (potent_1.use_angang__) {
		eangang2_(&i__);
	    }
	    if (potent_1.use_opbend__) {
		eopbend2_(&i__);
	    }
	    if (potent_1.use_opdist__) {
		eopdist2_(&i__);
	    }
	    if (potent_1.use_improp__) {
		eimprop2_(&i__);
	    }
	    if (potent_1.use_imptor__) {
		eimptor2_(&i__);
	    }
	    if (potent_1.use_tors__) {
		etors2_(&i__);
	    }
	    if (potent_1.use_pitors__) {
		epitors2_(&i__);
	    }
	    if (potent_1.use_strtor__) {
		estrtor2_(&i__);
	    }
	    if (potent_1.use_tortor__) {
		etortor2_(&i__);
	    }

/*     call the van der Waals Hessian component routines */

	    if (potent_1.use_vdw__) {
		if (s_cmp(vdwpot_1.vdwtyp, "LENNARD-JONES", (ftnlen)13, (
			ftnlen)13) == 0) {
		    elj2_(&i__, xred, yred, zred);
		} else if (s_cmp(vdwpot_1.vdwtyp, "BUCKINGHAM", (ftnlen)13, (
			ftnlen)10) == 0) {
		    ebuck2_(&i__, xred, yred, zred);
		} else if (s_cmp(vdwpot_1.vdwtyp, "MM3-HBOND", (ftnlen)13, (
			ftnlen)9) == 0) {
		    emm3hb2_(&i__, xred, yred, zred);
		} else if (s_cmp(vdwpot_1.vdwtyp, "BUFFERED-14-7", (ftnlen)13,
			 (ftnlen)13) == 0) {
		    ehal2_(&i__, xred, yred, zred);
		} else if (s_cmp(vdwpot_1.vdwtyp, "GAUSSIAN", (ftnlen)13, (
			ftnlen)8) == 0) {
		    egauss2_(&i__, xred, yred, zred);
		}
	    }

/*     call the electrostatic Hessian component routines */

	    if (potent_1.use_charge__) {
		echarge2_(&i__);
	    }
	    if (potent_1.use_chgdpl__) {
		echgdpl2_(&i__);
	    }
	    if (potent_1.use_dipole__) {
		edipole2_(&i__);
	    }
	    if (potent_1.use_mpole__ || potent_1.use_polar__) {
		empole2_(&i__);
	    }
	    if (potent_1.use_rxnfld__) {
		erxnfld2_(&i__);
	    }

/*     call any miscellaneous Hessian component routines */

	    if (potent_1.use_solv__) {
		esolv2_(&i__);
	    }
	    if (potent_1.use_metal__) {
		emetal2_(&i__);
	    }
	    if (potent_1.use_geom__) {
		egeom2_(&i__);
	    }
	    if (potent_1.use_extra__) {
		extra2_(&i__);
	    }

/*     set the diagonal Hessian matrix elements */

	    hdiag_ref(1, i__) = hdiag_ref(1, i__) + hessx_ref(1, i__);
	    hdiag_ref(2, i__) = hdiag_ref(2, i__) + hessy_ref(2, i__);
	    hdiag_ref(3, i__) = hdiag_ref(3, i__) + hessz_ref(3, i__);

/*     search each 3x3 block to see which blocks will be kept */

	    i__2 = atoms_1.n;
	    for (k = i__ + 1; k <= i__2; ++k) {
		keep[k - 1] = FALSE_;
		if (usage_1.use[k - 1]) {
/* Computing MAX */
		    d__10 = (d__1 = hessx_ref(1, k), abs(d__1)), d__11 = (
			    d__2 = hessx_ref(2, k), abs(d__2)), d__10 = max(
			    d__10,d__11), d__11 = (d__3 = hessx_ref(3, k), 
			    abs(d__3)), d__10 = max(d__10,d__11), d__11 = (
			    d__4 = hessy_ref(1, k), abs(d__4)), d__10 = max(
			    d__10,d__11), d__11 = (d__5 = hessy_ref(2, k), 
			    abs(d__5)), d__10 = max(d__10,d__11), d__11 = (
			    d__6 = hessy_ref(3, k), abs(d__6)), d__10 = max(
			    d__10,d__11), d__11 = (d__7 = hessz_ref(1, k), 
			    abs(d__7)), d__10 = max(d__10,d__11), d__11 = (
			    d__8 = hessz_ref(2, k), abs(d__8)), d__10 = max(
			    d__10,d__11), d__11 = (d__9 = hessz_ref(3, k), 
			    abs(d__9));
		    hmax = max(d__10,d__11);
		    if (hmax >= hescut_1.hesscut) {
			keep[k - 1] = TRUE_;
		    }
		}
	    }

/*     copy selected off-diagonal Hessian elements for current */
/*     atom into an indexed master list of Hessian elements; */
/*     if any elements of 3x3 block are kept, keep them all */

	    hinit_ref(1, i__) = nhess + 1;
	    for (j = 2; j <= 3; ++j) {
		++nhess;
		hindex[nhess] = i__ * 3 + j - 3;
		h__[nhess] = hessx_ref(j, i__);
	    }
	    i__2 = atoms_1.n;
	    for (k = i__ + 1; k <= i__2; ++k) {
		if (keep[k - 1]) {
		    for (j = 1; j <= 3; ++j) {
			++nhess;
			hindex[nhess] = k * 3 + j - 3;
			h__[nhess] = hessx_ref(j, k);
		    }
		}
	    }
	    hstop_ref(1, i__) = nhess;
	    hinit_ref(2, i__) = nhess + 1;
	    ++nhess;
	    hindex[nhess] = i__ * 3;
	    h__[nhess] = hessy_ref(3, i__);
	    i__2 = atoms_1.n;
	    for (k = i__ + 1; k <= i__2; ++k) {
		if (keep[k - 1]) {
		    for (j = 1; j <= 3; ++j) {
			++nhess;
			hindex[nhess] = k * 3 + j - 3;
			h__[nhess] = hessy_ref(j, k);
		    }
		}
	    }
	    hstop_ref(2, i__) = nhess;
	    hinit_ref(3, i__) = nhess + 1;
	    i__2 = atoms_1.n;
	    for (k = i__ + 1; k <= i__2; ++k) {
		if (keep[k - 1]) {
		    for (j = 1; j <= 3; ++j) {
			++nhess;
			hindex[nhess] = k * 3 + j - 3;
			h__[nhess] = hessz_ref(j, k);
		    }
		}
	    }
	    hstop_ref(3, i__) = nhess;

/*     check for storage of too many Hessian elements */

	    if (nhess > 1000000) {
		io___13.ciunit = iounit_1.iout;
		s_wsfe(&io___13);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nhess, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&c_b12, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&hescut_1.hesscut, (ftnlen)sizeof(
			doublereal));
		e_wsfe();
		io___14.ciunit = iounit_1.iout;
		s_wsfe(&io___14);
		e_wsfe();
		fatal_();
	    }
	}
    }

/*     print message telling how much storage was finally used */

    if (inform_1.verbose) {
	percent = (doublereal) nhess * 100. / (doublereal) (atoms_1.n * 3 * (
		atoms_1.n * 3 - 1) / 2);
	filled = (doublereal) nhess * 100. / 1e6;
	io___17.ciunit = iounit_1.iout;
	s_wsfe(&io___17);
	do_fio(&c__1, (char *)&nhess, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&percent, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&filled, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* hessian_ */

#undef hstop_ref
#undef hessz_ref
#undef hessy_ref
#undef hessx_ref
#undef hinit_ref
#undef hdiag_ref


