/* mechanic.f -- translated by f2c (version 20050501).
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
    doublereal vdwcut, chgcut, dplcut, mpolecut, vdwtaper, chgtaper, dpltaper,
	     mpoletaper, ewaldcut;
    logical use_ewald__, use_lights__, use_list__, use_vlist__, use_clist__, 
	    use_mlist__;
} cutoff_;

#define cutoff_1 cutoff_

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
    logical use_vcorr__;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine mechanic  --  initialize molecular mechanics  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "mechanic" sets up needed parameters for the potential energy */
/*     calculation and reads in many of the user selectable options */


/* Subroutine */ int mechanic_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 MECHANIC  --  Some Required Potential En"
	    "ergy\002,\002 Parameters are Undefined\002)";

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void);

    /* Local variables */
    extern /* Subroutine */ int molecule_(void), unitcell_(void), torsions_(
	    void), kvdw_(void), field_(void), fatal_(void), kbond_(void), 
	    kgeom_(void), bonds_(void), katom_(void), rings_(void), ksolv_(
	    void), kurey_(void), ktors_(void), kangle_(void), attach_(void), 
	    kewald_(void), angles_(void), active_(void), kmetal_(void), 
	    kmpole_(void), kpolar_(void), korbit_(void), mutate_(void), 
	    bitors_(void), kcharge_(void), kangang_(void), kopbend_(void), 
	    lattice_(void), kdipole_(void), orbital_(void), flatten_(void), 
	    kstrbnd_(void), cutoffs_(void), kopdist_(void), cluster_(void), 
	    kimprop_(void), kimptor_(void), polymer_(void), kpitors_(void), 
	    ktortor_(void), kstrtor_(void);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };




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




/*     set the bonded connectivity lists and active atoms */



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


    attach_();
    active_();

/*     find bonds, angles, torsions, bitorsions and small rings */

    bonds_();
    angles_();
    torsions_();
    bitors_();
    rings_();

/*     find unit cell type, lattice parameters and cutoff values */

    unitcell_();
    lattice_();
    polymer_();
    cutoffs_();

/*     setup needed for potential energy smoothing methods */

    flatten_();

/*     get the force field parameters and assign atom types */

    field_();
    katom_();

/*     assign atoms to molcules and set the atom groups */

    molecule_();
    cluster_();

/*     find any pisystem atoms, bonds and torsional angles */

    orbital_();

/*     assign bond, angle and cross term potential parameters */

    if (potent_1.use_bond__ || potent_1.use_strbnd__ || potent_1.use_strtor__ 
	    || potent_1.use_vdw__ && s_cmp(vdwpot_1.vdwtyp, "MM3-HBOND", (
	    ftnlen)13, (ftnlen)9) == 0) {
	kbond_();
    }
    if (potent_1.use_angle__ || potent_1.use_strbnd__ || 
	    potent_1.use_angang__) {
	kangle_();
    }
    if (potent_1.use_strbnd__) {
	kstrbnd_();
    }
    if (potent_1.use_urey__) {
	kurey_();
    }
    if (potent_1.use_angang__) {
	kangang_();
    }

/*     assign out-of-plane deformation potential parameters */

    if (potent_1.use_angle__ || potent_1.use_opbend__) {
	kopbend_();
    }
    if (potent_1.use_angle__ || potent_1.use_opdist__) {
	kopdist_();
    }
    if (potent_1.use_improp__) {
	kimprop_();
    }
    if (potent_1.use_imptor__) {
	kimptor_();
    }

/*     assign torsion and torsion cross term potential parameters */

    if (potent_1.use_tors__ || potent_1.use_strtor__ || potent_1.use_tortor__)
	     {
	ktors_();
    }
    if (potent_1.use_pitors__) {
	kpitors_();
    }
    if (potent_1.use_strtor__) {
	kstrtor_();
    }
    if (potent_1.use_tortor__) {
	ktortor_();
    }

/*     assign van der Waals and electrostatic potential parameters */

    if (potent_1.use_vdw__ || potent_1.use_solv__) {
	kvdw_();
    }
    if (potent_1.use_charge__ || potent_1.use_chgdpl__ || potent_1.use_solv__)
	     {
	kcharge_();
    }
    if (potent_1.use_dipole__ || potent_1.use_chgdpl__) {
	kdipole_();
    }
    if (potent_1.use_mpole__ || potent_1.use_polar__ || potent_1.use_solv__ ||
	     potent_1.use_rxnfld__) {
	kmpole_();
    }
    if (potent_1.use_polar__ || potent_1.use_solv__) {
	kpolar_();
    }
    if (cutoff_1.use_ewald__) {
	kewald_();
    }

/*     assign solvation, metal, pisystem and restraint parameters */

    if (potent_1.use_solv__) {
	ksolv_();
    }
    if (potent_1.use_metal__) {
	kmetal_();
    }
    if (potent_1.use_orbit__) {
	korbit_();
    }
    if (potent_1.use_geom__) {
	kgeom_();
    }

/*     set hybrid parameter values for free energy perturbation */

    mutate_();

/*     quit if essential parameter information is missing */

    if (inform_1.abort) {
	io___1.ciunit = iounit_1.iout;
	s_wsfe(&io___1);
	e_wsfe();
	fatal_();
    }
    return 0;
} /* mechanic_ */

