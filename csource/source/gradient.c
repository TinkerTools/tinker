/* gradient.f -- translated by f2c (version 20050501).
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
    doublereal einter;
} inter_;

#define inter_1 inter_

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
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_

struct {
    doublereal vir[9]	/* was [3][3] */;
} virial_;

#define virial_1 virial_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine gradient  --  find energy & gradient components  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "gradient" calls subroutines to calculate the potential energy */
/*     and first derivatives with respect to Cartesian coordinates */


/* Subroutine */ int gradient_(doublereal *energy, doublereal *derivs)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int estrbnd1_(void), erxnfld1_(void), eopdist1_(
	    void), eimprop1_(void), eimptor1_(void), epitors1_(void), 
	    etortor1_(void), estrtor1_(void);
    static integer i__, j;
    extern /* Subroutine */ int elj1_(void), born_(void), ehal1_(void), 
	    piscf_(void), ebond1_(void), ebuck1_(void), egeom1_(void), 
	    extra1_(void), esolv1_(void), eurey1_(void), etors1_(void), 
	    emm3hb1_(void);
    static doublereal cutoff;
    extern /* Subroutine */ int bounds_(void), nblist_(void), eangle1_(void), 
	    emetal1_(void), empole1_(void), egauss1_(void), replica_(
	    doublereal *), echarge1_(void), eangang1_(void), echgdpl1_(void), 
	    eopbend1_(void), edipole1_(void);


#define dea_ref(a_1,a_2) deriv_1.dea[(a_2)*3 + a_1 - 4]
#define deb_ref(a_1,a_2) deriv_1.deb[(a_2)*3 + a_1 - 4]
#define dec_ref(a_1,a_2) deriv_1.dec[(a_2)*3 + a_1 - 4]
#define ded_ref(a_1,a_2) deriv_1.ded[(a_2)*3 + a_1 - 4]
#define deg_ref(a_1,a_2) deriv_1.deg[(a_2)*3 + a_1 - 4]
#define dem_ref(a_1,a_2) deriv_1.dem[(a_2)*3 + a_1 - 4]
#define dep_ref(a_1,a_2) deriv_1.dep[(a_2)*3 + a_1 - 4]
#define der_ref(a_1,a_2) deriv_1.der[(a_2)*3 + a_1 - 4]
#define des_ref(a_1,a_2) deriv_1.des[(a_2)*3 + a_1 - 4]
#define det_ref(a_1,a_2) deriv_1.det[(a_2)*3 + a_1 - 4]
#define dev_ref(a_1,a_2) deriv_1.dev[(a_2)*3 + a_1 - 4]
#define dex_ref(a_1,a_2) deriv_1.dex[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define deaa_ref(a_1,a_2) deriv_1.deaa[(a_2)*3 + a_1 - 4]
#define deba_ref(a_1,a_2) deriv_1.deba[(a_2)*3 + a_1 - 4]
#define decd_ref(a_1,a_2) deriv_1.decd[(a_2)*3 + a_1 - 4]
#define deid_ref(a_1,a_2) deriv_1.deid[(a_2)*3 + a_1 - 4]
#define delf_ref(a_1,a_2) deriv_1.delf[(a_2)*3 + a_1 - 4]
#define debt_ref(a_1,a_2) deriv_1.debt[(a_2)*3 + a_1 - 4]
#define deub_ref(a_1,a_2) deriv_1.deub[(a_2)*3 + a_1 - 4]
#define deit_ref(a_1,a_2) deriv_1.deit[(a_2)*3 + a_1 - 4]
#define dept_ref(a_1,a_2) deriv_1.dept[(a_2)*3 + a_1 - 4]
#define dett_ref(a_1,a_2) deriv_1.dett[(a_2)*3 + a_1 - 4]
#define deopb_ref(a_1,a_2) deriv_1.deopb[(a_2)*3 + a_1 - 4]
#define deopd_ref(a_1,a_2) deriv_1.deopd[(a_2)*3 + a_1 - 4]
#define desum_ref(a_1,a_2) deriv_1.desum[(a_2)*3 + a_1 - 4]
#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  inter.i  --  sum of intermolecular energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     einter   total intermolecular potential energy */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  virial.i  --  components of internal virial tensor  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     vir    total internal virial Cartesian tensor components */




/*     zero out each of the potential energy components */

    /* Parameter adjustments */
    derivs -= 4;

    /* Function Body */
    energi_1.eb = 0.;
    energi_1.ea = 0.;
    energi_1.eba = 0.;
    energi_1.eub = 0.;
    energi_1.eaa = 0.;
    energi_1.eopb = 0.;
    energi_1.eopd = 0.;
    energi_1.eid = 0.;
    energi_1.eit = 0.;
    energi_1.et = 0.;
    energi_1.ept = 0.;
    energi_1.ebt = 0.;
    energi_1.ett = 0.;
    energi_1.ev = 0.;
    energi_1.ec = 0.;
    energi_1.ecd = 0.;
    energi_1.ed = 0.;
    energi_1.em = 0.;
    energi_1.ep = 0.;
    energi_1.er = 0.;
    energi_1.es = 0.;
    energi_1.elf = 0.;
    energi_1.eg = 0.;
    energi_1.ex = 0.;

/*     zero out each of the first derivative components */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    deb_ref(j, i__) = 0.;
	    dea_ref(j, i__) = 0.;
	    deba_ref(j, i__) = 0.;
	    deub_ref(j, i__) = 0.;
	    deaa_ref(j, i__) = 0.;
	    deopb_ref(j, i__) = 0.;
	    deopd_ref(j, i__) = 0.;
	    deid_ref(j, i__) = 0.;
	    deit_ref(j, i__) = 0.;
	    det_ref(j, i__) = 0.;
	    dept_ref(j, i__) = 0.;
	    debt_ref(j, i__) = 0.;
	    dett_ref(j, i__) = 0.;
	    dev_ref(j, i__) = 0.;
	    dec_ref(j, i__) = 0.;
	    decd_ref(j, i__) = 0.;
	    ded_ref(j, i__) = 0.;
	    dem_ref(j, i__) = 0.;
	    dep_ref(j, i__) = 0.;
	    der_ref(j, i__) = 0.;
	    des_ref(j, i__) = 0.;
	    delf_ref(j, i__) = 0.;
	    deg_ref(j, i__) = 0.;
	    dex_ref(j, i__) = 0.;
	}
    }

/*     zero out the virial and the intermolecular energy */

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    vir_ref(j, i__) = 0.;
	}
    }
    inter_1.einter = 0.;

/*     maintain any periodic boundary conditions */

    if (bound_1.use_bounds__ && ! rigid_1.use_rigid__) {
	bounds_();
    }

/*     update the pairwise interaction neighbor lists */

    if (cutoff_1.use_list__) {
	nblist_();
    }

/*     remove any previous use of the replicates method */

    cutoff = 0.;
    replica_(&cutoff);

/*     many implicit solvation models require Born radii */

    if (potent_1.use_born__) {
	born_();
    }

/*     alter bond and torsion constants for pisystem */

    if (potent_1.use_orbit__) {
	piscf_();
    }

/*     call the local geometry energy and gradient routines */

    if (potent_1.use_bond__) {
	ebond1_();
    }
    if (potent_1.use_angle__) {
	eangle1_();
    }
    if (potent_1.use_strbnd__) {
	estrbnd1_();
    }
    if (potent_1.use_urey__) {
	eurey1_();
    }
    if (potent_1.use_angang__) {
	eangang1_();
    }
    if (potent_1.use_opbend__) {
	eopbend1_();
    }
    if (potent_1.use_opdist__) {
	eopdist1_();
    }
    if (potent_1.use_improp__) {
	eimprop1_();
    }
    if (potent_1.use_imptor__) {
	eimptor1_();
    }
    if (potent_1.use_tors__) {
	etors1_();
    }
    if (potent_1.use_pitors__) {
	epitors1_();
    }
    if (potent_1.use_strtor__) {
	estrtor1_();
    }
    if (potent_1.use_tortor__) {
	etortor1_();
    }

/*     call the van der Waals energy and gradient routines */

    if (potent_1.use_vdw__) {
	if (s_cmp(vdwpot_1.vdwtyp, "LENNARD-JONES", (ftnlen)13, (ftnlen)13) ==
		 0) {
	    elj1_();
	}
	if (s_cmp(vdwpot_1.vdwtyp, "BUCKINGHAM", (ftnlen)13, (ftnlen)10) == 0)
		 {
	    ebuck1_();
	}
	if (s_cmp(vdwpot_1.vdwtyp, "MM3-HBOND", (ftnlen)13, (ftnlen)9) == 0) {
	    emm3hb1_();
	}
	if (s_cmp(vdwpot_1.vdwtyp, "BUFFERED-14-7", (ftnlen)13, (ftnlen)13) ==
		 0) {
	    ehal1_();
	}
	if (s_cmp(vdwpot_1.vdwtyp, "GAUSSIAN", (ftnlen)13, (ftnlen)8) == 0) {
	    egauss1_();
	}
    }

/*     call the electrostatic energy and gradient routines */

    if (potent_1.use_charge__) {
	echarge1_();
    }
    if (potent_1.use_chgdpl__) {
	echgdpl1_();
    }
    if (potent_1.use_dipole__) {
	edipole1_();
    }
    if (potent_1.use_mpole__ || potent_1.use_polar__) {
	empole1_();
    }
    if (potent_1.use_rxnfld__) {
	erxnfld1_();
    }

/*     call any miscellaneous energy and gradient routines */

    if (potent_1.use_solv__) {
	esolv1_();
    }
    if (potent_1.use_metal__) {
	emetal1_();
    }
    if (potent_1.use_geom__) {
	egeom1_();
    }
    if (potent_1.use_extra__) {
	extra1_();
    }

/*     sum up to get the total energy and first derivatives */

    energi_1.esum = energi_1.eb + energi_1.ea + energi_1.eba + energi_1.eub + 
	    energi_1.eaa + energi_1.eopb + energi_1.eopd + energi_1.eid + 
	    energi_1.eit + energi_1.et + energi_1.ept + energi_1.ebt + 
	    energi_1.ett + energi_1.ev + energi_1.ec + energi_1.ecd + 
	    energi_1.ed + energi_1.em + energi_1.ep + energi_1.er + 
	    energi_1.es + energi_1.elf + energi_1.eg + energi_1.ex;
    *energy = energi_1.esum;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    desum_ref(j, i__) = deb_ref(j, i__) + dea_ref(j, i__) + deba_ref(
		    j, i__) + deub_ref(j, i__) + deaa_ref(j, i__) + deopb_ref(
		    j, i__) + deopd_ref(j, i__) + deid_ref(j, i__) + deit_ref(
		    j, i__) + det_ref(j, i__) + dept_ref(j, i__) + debt_ref(j,
		     i__) + dett_ref(j, i__) + dev_ref(j, i__) + dec_ref(j, 
		    i__) + decd_ref(j, i__) + ded_ref(j, i__) + dem_ref(j, 
		    i__) + dep_ref(j, i__) + der_ref(j, i__) + des_ref(j, i__)
		     + delf_ref(j, i__) + deg_ref(j, i__) + dex_ref(j, i__);
	    derivs_ref(j, i__) = desum_ref(j, i__);
	}
    }
    return 0;
} /* gradient_ */

#undef derivs_ref
#undef desum_ref
#undef deopd_ref
#undef deopb_ref
#undef dett_ref
#undef dept_ref
#undef deit_ref
#undef deub_ref
#undef debt_ref
#undef delf_ref
#undef deid_ref
#undef decd_ref
#undef deba_ref
#undef deaa_ref
#undef vir_ref
#undef dex_ref
#undef dev_ref
#undef det_ref
#undef des_ref
#undef der_ref
#undef dep_ref
#undef dem_ref
#undef deg_ref
#undef ded_ref
#undef dec_ref
#undef deb_ref
#undef dea_ref


