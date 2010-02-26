/* analysis.f -- translated by f2c (version 20050501).
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
    doublereal aesum[25000], aeb[25000], aea[25000], aeba[25000], aeub[25000],
	     aeaa[25000], aeopb[25000], aeopd[25000], aeid[25000], aeit[25000]
	    , aet[25000], aept[25000], aebt[25000], aett[25000], aev[25000], 
	    aec[25000], aecd[25000], aed[25000], aem[25000], aep[25000], aer[
	    25000], aes[25000], aelf[25000], aeg[25000], aex[25000];
} analyz_;

#define analyz_1 analyz_

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
    doublereal esum, eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, 
	    ebt, ett, ev, ec, ecd, ed, em, ep, er, es, elf, eg, ex;
} energi_;

#define energi_1 energi_

struct {
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

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
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
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
/*     ##  subroutine analysis  --  energy components and analysis  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "analysis" calls the series of routines needed to calculate */
/*     the potential energy and perform energy partitioning analysis */
/*     in terms of type of interaction or atom number */


/* Subroutine */ int analysis_(doublereal *energy)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int estrbnd3_(void), erxnfld3_(void), eopdist3_(
	    void), eimprop3_(void), eimptor3_(void), epitors3_(void), 
	    etortor3_(void), estrtor3_(void);
    static integer i__;
    extern /* Subroutine */ int elj3_(void), born_(void), ehal3_(void), 
	    piscf_(void), ebond3_(void), ebuck3_(void), egeom3_(void), 
	    extra3_(void), esolv3_(void), eurey3_(void), etors3_(void), 
	    emm3hb3_(void);
    static doublereal cutoff;
    extern /* Subroutine */ int bounds_(void), nblist_(void), eangle3_(void), 
	    emetal3_(void), empole3_(void), egauss3_(void), replica_(
	    doublereal *), echarge3_(void), eangang3_(void), echgdpl3_(void), 
	    eopbend3_(void), edipole3_(void);



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
/*     ##  analyz.i  --  energy components partitioned over atoms  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     aesum   total potential energy partitioned over atoms */
/*     aeb     bond stretch energy partitioned over atoms */
/*     aea     angle bend energy partitioned over atoms */
/*     aeba    stretch-bend energy partitioned over atoms */
/*     aeub    Urey-Bradley energy partitioned over atoms */
/*     aeaa    angle-angle energy partitioned over atoms */
/*     aeopb   out-of-plane bend energy partitioned over atoms */
/*     aeopd   out-of-plane distance energy partitioned over atoms */
/*     aeid    improper dihedral energy partitioned over atoms */
/*     aeit    improper torsion energy partitioned over atoms */
/*     aet     torsional energy partitioned over atoms */
/*     aept    pi-orbital torsion energy partitioned over atoms */
/*     aebt    stretch-torsion energy partitioned over atoms */
/*     aett    torsion-torsion energy partitioned over atoms */
/*     aev     van der Waals energy partitioned over atoms */
/*     aec     charge-charge energy partitioned over atoms */
/*     aecd    charge-dipole energy partitioned over atoms */
/*     aed     dipole-dipole energy partitioned over atoms */
/*     aem     multipole energy partitioned over atoms */
/*     aep     polarization energy partitioned over atoms */
/*     aer     reaction field energy partitioned over atoms */
/*     aes     solvation energy partitioned over atoms */
/*     aelf    metal ligand field energy partitioned over atoms */
/*     aeg     geometric restraint energy partitioned over atoms */
/*     aex     extra energy term partitioned over atoms */




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




/*     zero out each of the potential energy components */

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

/*     zero out energy partitioning components for each atom */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	analyz_1.aeb[i__ - 1] = 0.;
	analyz_1.aea[i__ - 1] = 0.;
	analyz_1.aeba[i__ - 1] = 0.;
	analyz_1.aeub[i__ - 1] = 0.;
	analyz_1.aeaa[i__ - 1] = 0.;
	analyz_1.aeopb[i__ - 1] = 0.;
	analyz_1.aeopd[i__ - 1] = 0.;
	analyz_1.aeid[i__ - 1] = 0.;
	analyz_1.aeit[i__ - 1] = 0.;
	analyz_1.aet[i__ - 1] = 0.;
	analyz_1.aept[i__ - 1] = 0.;
	analyz_1.aebt[i__ - 1] = 0.;
	analyz_1.aett[i__ - 1] = 0.;
	analyz_1.aev[i__ - 1] = 0.;
	analyz_1.aec[i__ - 1] = 0.;
	analyz_1.aecd[i__ - 1] = 0.;
	analyz_1.aed[i__ - 1] = 0.;
	analyz_1.aem[i__ - 1] = 0.;
	analyz_1.aep[i__ - 1] = 0.;
	analyz_1.aer[i__ - 1] = 0.;
	analyz_1.aes[i__ - 1] = 0.;
	analyz_1.aelf[i__ - 1] = 0.;
	analyz_1.aeg[i__ - 1] = 0.;
	analyz_1.aex[i__ - 1] = 0.;
    }

/*     zero out the total intermolecular energy */

    inter_1.einter = 0.;

/*     maintain any periodic boundary conditions */

    if (bound_1.use_bounds__ && ! group_1.use_group__) {
	bounds_();
    }

/*     remove any previous use of the replicates method */

    cutoff = 0.;
    replica_(&cutoff);

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

/*     call the local geometry energy component routines */

    if (potent_1.use_bond__) {
	ebond3_();
    }
    if (potent_1.use_angle__) {
	eangle3_();
    }
    if (potent_1.use_strbnd__) {
	estrbnd3_();
    }
    if (potent_1.use_urey__) {
	eurey3_();
    }
    if (potent_1.use_angang__) {
	eangang3_();
    }
    if (potent_1.use_opbend__) {
	eopbend3_();
    }
    if (potent_1.use_opdist__) {
	eopdist3_();
    }
    if (potent_1.use_improp__) {
	eimprop3_();
    }
    if (potent_1.use_imptor__) {
	eimptor3_();
    }
    if (potent_1.use_tors__) {
	etors3_();
    }
    if (potent_1.use_pitors__) {
	epitors3_();
    }
    if (potent_1.use_strtor__) {
	estrtor3_();
    }
    if (potent_1.use_tortor__) {
	etortor3_();
    }

/*     call the van der Waals energy component routines */

    if (potent_1.use_vdw__) {
	if (s_cmp(vdwpot_1.vdwtyp, "LENNARD-JONES", (ftnlen)13, (ftnlen)13) ==
		 0) {
	    elj3_();
	}
	if (s_cmp(vdwpot_1.vdwtyp, "BUCKINGHAM", (ftnlen)13, (ftnlen)10) == 0)
		 {
	    ebuck3_();
	}
	if (s_cmp(vdwpot_1.vdwtyp, "MM3-HBOND", (ftnlen)13, (ftnlen)9) == 0) {
	    emm3hb3_();
	}
	if (s_cmp(vdwpot_1.vdwtyp, "BUFFERED-14-7", (ftnlen)13, (ftnlen)13) ==
		 0) {
	    ehal3_();
	}
	if (s_cmp(vdwpot_1.vdwtyp, "GAUSSIAN", (ftnlen)13, (ftnlen)8) == 0) {
	    egauss3_();
	}
    }

/*     call the electrostatic energy component routines */

    if (potent_1.use_charge__) {
	echarge3_();
    }
    if (potent_1.use_chgdpl__) {
	echgdpl3_();
    }
    if (potent_1.use_dipole__) {
	edipole3_();
    }
    if (potent_1.use_mpole__ || potent_1.use_polar__) {
	empole3_();
    }
    if (potent_1.use_rxnfld__) {
	erxnfld3_();
    }

/*     call any miscellaneous energy component routines */

    if (potent_1.use_solv__) {
	esolv3_();
    }
    if (potent_1.use_metal__) {
	emetal3_();
    }
    if (potent_1.use_geom__) {
	egeom3_();
    }
    if (potent_1.use_extra__) {
	extra3_();
    }

/*     sum up to give the total potential energy */

    energi_1.esum = energi_1.eb + energi_1.ea + energi_1.eba + energi_1.eub + 
	    energi_1.eaa + energi_1.eopb + energi_1.eopd + energi_1.eid + 
	    energi_1.eit + energi_1.et + energi_1.ept + energi_1.ebt + 
	    energi_1.ett + energi_1.ev + energi_1.ec + energi_1.ecd + 
	    energi_1.ed + energi_1.em + energi_1.ep + energi_1.er + 
	    energi_1.es + energi_1.elf + energi_1.eg + energi_1.ex;
    *energy = energi_1.esum;

/*     sum up to give the total potential energy per atom */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	analyz_1.aesum[i__ - 1] = analyz_1.aeb[i__ - 1] + analyz_1.aea[i__ - 
		1] + analyz_1.aeba[i__ - 1] + analyz_1.aeub[i__ - 1] + 
		analyz_1.aeaa[i__ - 1] + analyz_1.aeopb[i__ - 1] + 
		analyz_1.aeopd[i__ - 1] + analyz_1.aeid[i__ - 1] + 
		analyz_1.aeit[i__ - 1] + analyz_1.aet[i__ - 1] + 
		analyz_1.aept[i__ - 1] + analyz_1.aebt[i__ - 1] + 
		analyz_1.aett[i__ - 1] + analyz_1.aev[i__ - 1] + analyz_1.aec[
		i__ - 1] + analyz_1.aecd[i__ - 1] + analyz_1.aed[i__ - 1] + 
		analyz_1.aem[i__ - 1] + analyz_1.aep[i__ - 1] + analyz_1.aer[
		i__ - 1] + analyz_1.aes[i__ - 1] + analyz_1.aelf[i__ - 1] + 
		analyz_1.aeg[i__ - 1] + analyz_1.aex[i__ - 1];
    }
    return 0;
} /* analysis_ */

