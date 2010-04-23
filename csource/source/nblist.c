/* nblist.f -- translated by f2c (version 20050501).
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
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

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
    doublereal xbox, ybox, zbox, alpha, beta, gamma, xbox2, ybox2, zbox2, 
	    box34, lvec[9]	/* was [3][3] */, recip[9]	/* was [3][3] 
	    */, volbox, beta_sin__, beta_cos__, gamma_sin__, gamma_cos__, 
	    beta_term__, gamma_term__;
    logical orthogonal, monoclinic, triclinic, octahedron;
    char spacegrp[10];
} boxes_;

#define boxes_1 boxes_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal lbuffer, lbuf2, vbuf2, cbuf2, mbuf2;
    integer nvlst[25000], vlst[45000000]	/* was [1800][25000] */, 
	    nelst[25000], elst[30000000]	/* was [1200][25000] */;
    logical dovlst, doclst, domlst;
} neigh_;

#define neigh_1 neigh_

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
    integer nlight, kbx[25000], kby[25000], kbz[25000], kex[25000], key[25000]
	    , kez[25000], locx[200000], locy[200000], locz[200000], rgx[
	    200000], rgy[200000], rgz[200000];
} light_;

#define light_1 light_

struct {
    doublereal pchg[25000];
    integer nion, iion[25000], jion[25000], kion[25000], chglist[25000];
} charge_;

#define charge_1 charge_

struct {
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_



/*     ############################################################### */
/*     ##  COPYRIGHT (C) 2006 by David Gohara & Jay William Ponder  ## */
/*     ##                    All Rights Reserved                    ## */
/*     ############################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine nblist  --  maintain pairwise neighbor lists  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "nblist" constructs and maintains nonbonded pair neighbor lists */
/*     for vdw and electrostatic interactions */

/*     the routines "vlist", "clist" and "mlist" each implement two */
/*     methods for complete list rebuilds and two methods for list */
/*     updates; only one of each pair should be uncommented */


/* Subroutine */ int nblist_(void)
{
    extern /* Subroutine */ int clist_(void), mlist_(void), vlist_(void);



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




/*     update the vdw and electrostatic neighbor lists */



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


    if (potent_1.use_vdw__ && cutoff_1.use_vlist__) {
	vlist_();
    }
    if (potent_1.use_charge__ && cutoff_1.use_clist__) {
	clist_();
    }
    if ((potent_1.use_mpole__ || potent_1.use_polar__) && 
	    cutoff_1.use_mlist__) {
	mlist_();
    }
    return 0;
} /* nblist_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine vlist  --  build van der Waals neighbor lists  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "vlist" performs an update or a complete rebuild of the */
/*     van der Waals neighbor list */


/* Subroutine */ int vlist_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 VLIST  --  Pairwise Neighbor List canno"
	    "t\002,\002 be used with Replicas\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, k;
    static doublereal r2;
    static integer ii, iv;
    static doublereal xi, yi, zi, xr, yr, zr, rdn, xred[25000], yred[25000], 
	    zred[25000], xold[25000], yold[25000], zold[25000];
    extern /* Subroutine */ int fatal_(void);
    static logical reset[25000];
    extern /* Subroutine */ int vfull_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *), imagen_(doublereal *
	    , doublereal *, doublereal *);
    static doublereal radius;
    extern /* Subroutine */ int vbuild_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *), 
	    replica_(doublereal *);

    /* Fortran I/O blocks */
    static cilist io___9 = { 0, 0, 0, fmt_10, 0 };




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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  boxes.i  --  parameters for periodic boundary conditions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     xbox        length of a-axis of periodic box in Angstroms */
/*     ybox        length of b-axis of periodic box in Angstroms */
/*     zbox        length of c-axis of periodic box in Angstroms */
/*     alpha       angle between b- and c-axes of box in degrees */
/*     beta        angle between a- and c-axes of box in degrees */
/*     gamma       angle between a- and b-axes of box in degrees */
/*     xbox2       half of the a-axis length of periodic box */
/*     ybox2       half of the b-axis length of periodic box */
/*     zbox2       half of the c-axis length of periodic box */
/*     box34       three-fourths axis length of truncated octahedron */
/*     lvec        real space lattice vectors as matrix rows */
/*     recip       reciprocal lattice vectors as matrix columns */
/*     volbox      volume in Ang**3 of the periodic box */
/*     beta_sin    sine of the beta periodic box angle */
/*     beta_cos    cosine of the beta periodic box angle */
/*     gamma_sin   sine of the gamma periodic box angle */
/*     gamma_cos   cosine of the gamma periodic box angle */
/*     beta_term   term used in generating triclinic box */
/*     gamma_term  term used in generating triclinic box */
/*     orthogonal  flag to mark periodic box as orthogonal */
/*     monoclinic  flag to mark periodic box as monoclinic */
/*     triclinic   flag to mark periodic box as triclinic */
/*     octahedron  flag to mark box as truncated octahedron */
/*     spacegrp    space group symbol for the unitcell type */




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




/*     apply reduction factors to find coordinates for each site */

    i__1 = vdw_1.nvdw;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = vdw_1.ivdw[i__ - 1];
	iv = vdw_1.ired[ii - 1];
	rdn = vdw_1.kred[ii - 1];
	xred[i__ - 1] = rdn * (atoms_1.x[ii - 1] - atoms_1.x[iv - 1]) + 
		atoms_1.x[iv - 1];
	yred[i__ - 1] = rdn * (atoms_1.y[ii - 1] - atoms_1.y[iv - 1]) + 
		atoms_1.y[iv - 1];
	zred[i__ - 1] = rdn * (atoms_1.z__[ii - 1] - atoms_1.z__[iv - 1]) + 
		atoms_1.z__[iv - 1];
    }

/*     neighbor list cannot be used with the replicates method */

    radius = sqrt(neigh_1.vbuf2);
    replica_(&radius);
    if (bound_1.use_replica__) {
	io___9.ciunit = iounit_1.iout;
	s_wsfe(&io___9);
	e_wsfe();
	fatal_();
    }

/*     perform a complete list build instead of an update */

    if (neigh_1.dovlst) {
	neigh_1.dovlst = FALSE_;
	if (boxes_1.octahedron) {
	    i__1 = vdw_1.nvdw;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		vbuild_(&i__, xred, yred, zred, xold, yold, zold);
	    }
	} else {
	    vfull_(xred, yred, zred, xold, yold, zold);
	}
	return 0;
    }

/*     update sites whose displacement exceeds half the buffer */

    i__1 = vdw_1.nvdw;
    for (i__ = 1; i__ <= i__1; ++i__) {
	reset[i__ - 1] = FALSE_;
	xi = xred[i__ - 1];
	yi = yred[i__ - 1];
	zi = zred[i__ - 1];
	xr = xi - xold[i__ - 1];
	yr = yi - yold[i__ - 1];
	zr = zi - zold[i__ - 1];
	imagen_(&xr, &yr, &zr);
	r2 = xr * xr + yr * yr + zr * zr;
	if (neigh_1.dovlst || r2 >= neigh_1.lbuf2) {
	    vbuild_(&i__, xred, yred, zred, xold, yold, zold);
	    reset[i__ - 1] = TRUE_;
	    i__2 = i__ - 1;
	    for (k = 1; k <= i__2; ++k) {
		if (! reset[k - 1]) {
		    xr = xi - xred[k - 1];
		    yr = yi - yred[k - 1];
		    zr = zi - zred[k - 1];
		    imagen_(&xr, &yr, &zr);
		    r2 = xr * xr + yr * yr + zr * zr;

/*     perform an update by rebuilding lower numbered neighbors */

		    if (r2 <= neigh_1.vbuf2) {
			vbuild_(&k, xred, yred, zred, xold, yold, zold);
			reset[k - 1] = TRUE_;
		    }
		}
	    }
	}
    }
    return 0;
} /* vlist_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine vbuild  --  make vdw pair list for one site  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "vbuild" performs a complete rebuild of the van der Waals */
/*     pair neighbor list for a single site */


/* Subroutine */ int vbuild_(integer *i__, doublereal *xred, doublereal *yred,
	 doublereal *zred, doublereal *xold, doublereal *yold, doublereal *
	zold)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 VBUILD  --  Too many Neighbors;\002,\002"
	    " Increase MAXVLST\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer k;
    static doublereal r2, xi, yi, zi, xr, yr, zr;
    extern /* Subroutine */ int fatal_(void), imagen_(doublereal *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___30 = { 0, 0, 0, fmt_10, 0 };



#define vlst_ref(a_1,a_2) neigh_1.vlst[(a_2)*1800 + a_1 - 1801]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




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




/*     zero out the list and get coordinates for the site */

    /* Parameter adjustments */
    --zold;
    --yold;
    --xold;
    --zred;
    --yred;
    --xred;

    /* Function Body */
    neigh_1.nvlst[*i__ - 1] = 0;
    xi = xred[*i__];
    yi = yred[*i__];
    zi = zred[*i__];

/*     generate all neighbors for the site being rebuilt */

    i__1 = vdw_1.nvdw;
    for (k = *i__ + 1; k <= i__1; ++k) {
	xr = xi - xred[k];
	yr = yi - yred[k];
	zr = zi - zred[k];
	imagen_(&xr, &yr, &zr);
	r2 = xr * xr + yr * yr + zr * zr;
	if (r2 <= neigh_1.vbuf2) {
	    ++neigh_1.nvlst[*i__ - 1];
	    vlst_ref(neigh_1.nvlst[*i__ - 1], *i__) = k;
	}
    }

/*     store new coordinates to reflect update of the site */

    xold[*i__] = xi;
    yold[*i__] = yi;
    zold[*i__] = zi;

/*     check to see if the neighbor list is too long */

    if (neigh_1.nvlst[*i__ - 1] >= 1800) {
	io___30.ciunit = iounit_1.iout;
	s_wsfe(&io___30);
	e_wsfe();
	fatal_();
    }
    return 0;
} /* vbuild_ */

#undef vlst_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine vfull  --  make vdw pair list for all sites  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "vfull" performs a complete rebuild of the van der Waals */
/*     pair neighbor list for all sites using the method of lights */


/* Subroutine */ int vfull_(doublereal *xred, doublereal *yred, doublereal *
	zred, doublereal *xold, doublereal *yold, doublereal *zold)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 VFULL  --  Too many Neighbors;\002,\002 "
	    "Increase MAXVLST\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k;
    static doublereal r2, xi, yi, zi, xr, yr, zr, off;
    static integer kgy, kgz, stop;
    extern /* Subroutine */ int fatal_(void);
    static integer start;
    static doublereal xsort[200000], ysort[200000], zsort[200000];
    extern /* Subroutine */ int imagen_(doublereal *, doublereal *, 
	    doublereal *);
    static logical repeat;
    extern /* Subroutine */ int lights_(doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___50 = { 0, 0, 0, fmt_30, 0 };



#define vlst_ref(a_1,a_2) neigh_1.vlst[(a_2)*1800 + a_1 - 1801]



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




/*     transfer interaction site coordinates to sorting arrays */

    /* Parameter adjustments */
    --zold;
    --yold;
    --xold;
    --zred;
    --yred;
    --xred;

    /* Function Body */
    i__1 = vdw_1.nvdw;
    for (i__ = 1; i__ <= i__1; ++i__) {
	neigh_1.nvlst[i__ - 1] = 0;
	xold[i__] = xred[i__];
	yold[i__] = yred[i__];
	zold[i__] = zred[i__];
	xsort[i__ - 1] = xred[i__];
	ysort[i__ - 1] = yred[i__];
	zsort[i__ - 1] = zred[i__];
    }

/*     use the method of lights to generate neighbors */

    off = sqrt(neigh_1.vbuf2);
    lights_(&off, &vdw_1.nvdw, xsort, ysort, zsort);

/*     loop over all atoms computing the interactions */

    i__1 = vdw_1.nvdw;
    for (i__ = 1; i__ <= i__1; ++i__) {
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
L10:
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    k = light_1.locx[j - 1];
	    kgy = light_1.rgy[k - 1];
	    if (light_1.kby[i__ - 1] <= light_1.key[i__ - 1]) {
		if (kgy < light_1.kby[i__ - 1] || kgy > light_1.key[i__ - 1]) 
			{
		    goto L20;
		}
	    } else {
		if (kgy < light_1.kby[i__ - 1] && kgy > light_1.key[i__ - 1]) 
			{
		    goto L20;
		}
	    }
	    kgz = light_1.rgz[k - 1];
	    if (light_1.kbz[i__ - 1] <= light_1.kez[i__ - 1]) {
		if (kgz < light_1.kbz[i__ - 1] || kgz > light_1.kez[i__ - 1]) 
			{
		    goto L20;
		}
	    } else {
		if (kgz < light_1.kbz[i__ - 1] && kgz > light_1.kez[i__ - 1]) 
			{
		    goto L20;
		}
	    }
	    xr = xi - xsort[j - 1];
	    yr = yi - ysort[kgy - 1];
	    zr = zi - zsort[kgz - 1];
	    imagen_(&xr, &yr, &zr);
	    r2 = xr * xr + yr * yr + zr * zr;
	    if (r2 <= neigh_1.vbuf2) {
		if (i__ < k) {
		    ++neigh_1.nvlst[i__ - 1];
		    vlst_ref(neigh_1.nvlst[i__ - 1], i__) = k;
		} else {
		    ++neigh_1.nvlst[k - 1];
		    vlst_ref(neigh_1.nvlst[k - 1], k) = i__;
		}
	    }
L20:
	    ;
	}
	if (repeat) {
	    repeat = FALSE_;
	    start = light_1.kbx[i__ - 1] + 1;
	    stop = light_1.nlight;
	    goto L10;
	}
    }

/*     check to see if the neighbor lists are too long */

    i__1 = vdw_1.nvdw;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (neigh_1.nvlst[i__ - 1] >= 1800) {
	    io___50.ciunit = iounit_1.iout;
	    s_wsfe(&io___50);
	    e_wsfe();
	    fatal_();
	}
    }
    return 0;
} /* vfull_ */

#undef vlst_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine clist  --  build partial charge neighbor lists  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "clist" performs and update or a complete rebuild of the */
/*     electrostatic neighbor lists for partial charges */


/* Subroutine */ int clist_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 CLIST  --  Pairwise Neighbor List canno"
	    "t\002,\002 be used with Replicas\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, k;
    static doublereal r2;
    static integer ic, kc;
    static doublereal xi, yi, zi, xr, yr, zr, xold[25000], yold[25000], zold[
	    25000];
    extern /* Subroutine */ int fatal_(void), cfull_(doublereal *, doublereal 
	    *, doublereal *);
    static logical reset[25000];
    extern /* Subroutine */ int imagen_(doublereal *, doublereal *, 
	    doublereal *), cbuild_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal radius;
    extern /* Subroutine */ int replica_(doublereal *);

    /* Fortran I/O blocks */
    static cilist io___52 = { 0, 0, 0, fmt_10, 0 };




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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  boxes.i  --  parameters for periodic boundary conditions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     xbox        length of a-axis of periodic box in Angstroms */
/*     ybox        length of b-axis of periodic box in Angstroms */
/*     zbox        length of c-axis of periodic box in Angstroms */
/*     alpha       angle between b- and c-axes of box in degrees */
/*     beta        angle between a- and c-axes of box in degrees */
/*     gamma       angle between a- and b-axes of box in degrees */
/*     xbox2       half of the a-axis length of periodic box */
/*     ybox2       half of the b-axis length of periodic box */
/*     zbox2       half of the c-axis length of periodic box */
/*     box34       three-fourths axis length of truncated octahedron */
/*     lvec        real space lattice vectors as matrix rows */
/*     recip       reciprocal lattice vectors as matrix columns */
/*     volbox      volume in Ang**3 of the periodic box */
/*     beta_sin    sine of the beta periodic box angle */
/*     beta_cos    cosine of the beta periodic box angle */
/*     gamma_sin   sine of the gamma periodic box angle */
/*     gamma_cos   cosine of the gamma periodic box angle */
/*     beta_term   term used in generating triclinic box */
/*     gamma_term  term used in generating triclinic box */
/*     orthogonal  flag to mark periodic box as orthogonal */
/*     monoclinic  flag to mark periodic box as monoclinic */
/*     triclinic   flag to mark periodic box as triclinic */
/*     octahedron  flag to mark box as truncated octahedron */
/*     spacegrp    space group symbol for the unitcell type */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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




/*     neighbor list cannot be used with the replicates method */

    radius = sqrt(neigh_1.cbuf2);
    replica_(&radius);
    if (bound_1.use_replica__) {
	io___52.ciunit = iounit_1.iout;
	s_wsfe(&io___52);
	e_wsfe();
	fatal_();
    }

/*     perform a complete list build instead of an update */

    if (neigh_1.doclst) {
	neigh_1.doclst = FALSE_;
	if (boxes_1.octahedron) {
	    i__1 = charge_1.nion;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		cbuild_(&i__, xold, yold, zold);
	    }
	} else {
	    cfull_(xold, yold, zold);
	}
	return 0;
    }

/*     update sites whose displacement exceeds half the buffer */

    i__1 = charge_1.nion;
    for (i__ = 1; i__ <= i__1; ++i__) {
	reset[i__ - 1] = FALSE_;
	ic = charge_1.kion[i__ - 1];
	xi = atoms_1.x[ic - 1];
	yi = atoms_1.y[ic - 1];
	zi = atoms_1.z__[ic - 1];
	xr = xi - xold[i__ - 1];
	yr = yi - yold[i__ - 1];
	zr = zi - zold[i__ - 1];
	imagen_(&xr, &yr, &zr);
	r2 = xr * xr + yr * yr + zr * zr;
	if (neigh_1.doclst || r2 >= neigh_1.lbuf2) {
	    cbuild_(&i__, xold, yold, zold);
	    reset[i__ - 1] = TRUE_;
	    i__2 = i__ - 1;
	    for (k = 1; k <= i__2; ++k) {
		if (! reset[k - 1]) {
		    kc = charge_1.kion[k - 1];
		    xr = xi - atoms_1.x[kc - 1];
		    yr = yi - atoms_1.y[kc - 1];
		    zr = zi - atoms_1.z__[kc - 1];
		    imagen_(&xr, &yr, &zr);
		    r2 = xr * xr + yr * yr + zr * zr;

/*     rebuild the list of any possible lower numbered neighbor */

		    if (r2 <= neigh_1.cbuf2) {
			cbuild_(&k, xold, yold, zold);
			reset[k - 1] = TRUE_;
		    }
		}
	    }
	}
    }
    return 0;
} /* clist_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine cbuild  --  make charge pair list for one site  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "cbuild" performs a complete rebuild of the partial charge */
/*     electrostatic neighbor list for a single site */


/* Subroutine */ int cbuild_(integer *i__, doublereal *xold, doublereal *yold,
	 doublereal *zold)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 CBUILD  --  Too many Neighbors;\002,\002"
	    " Increase MAXELST\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer k;
    static doublereal r2;
    static integer ic, kc;
    static doublereal xi, yi, zi, xr, yr, zr;
    extern /* Subroutine */ int fatal_(void), imagen_(doublereal *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___78 = { 0, 0, 0, fmt_10, 0 };



#define elst_ref(a_1,a_2) neigh_1.elst[(a_2)*1200 + a_1 - 1201]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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




/*     zero out the list and get coordinates for the site */

    /* Parameter adjustments */
    --zold;
    --yold;
    --xold;

    /* Function Body */
    neigh_1.nelst[*i__ - 1] = 0;
    ic = charge_1.kion[*i__ - 1];
    xi = atoms_1.x[ic - 1];
    yi = atoms_1.y[ic - 1];
    zi = atoms_1.z__[ic - 1];

/*     generate all neighbors for the site being rebuilt */

    i__1 = charge_1.nion;
    for (k = *i__ + 1; k <= i__1; ++k) {
	kc = charge_1.kion[k - 1];
	xr = xi - atoms_1.x[kc - 1];
	yr = yi - atoms_1.y[kc - 1];
	zr = zi - atoms_1.z__[kc - 1];
	imagen_(&xr, &yr, &zr);
	r2 = xr * xr + yr * yr + zr * zr;
	if (r2 <= neigh_1.cbuf2) {
	    ++neigh_1.nelst[*i__ - 1];
	    elst_ref(neigh_1.nelst[*i__ - 1], *i__) = k;
	}
    }

/*     store new coordinates to reflect update of the site */

    xold[*i__] = xi;
    yold[*i__] = yi;
    zold[*i__] = zi;

/*     check to see if the neighbor list is too long */

    if (neigh_1.nelst[*i__ - 1] >= 1200) {
	io___78.ciunit = iounit_1.iout;
	s_wsfe(&io___78);
	e_wsfe();
	fatal_();
    }
    return 0;
} /* cbuild_ */

#undef elst_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine cfull  --  make charge pair list for all sites  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "cfull" performs a complete rebuild of the partial charge */
/*     pair neighbor list for all sites using the method of lights */


/* Subroutine */ int cfull_(doublereal *xold, doublereal *yold, doublereal *
	zold)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 CFULL  --  Too many Neighbors;\002,\002 "
	    "Increase MAXELST\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k;
    static doublereal r2;
    static integer ic;
    static doublereal xi, yi, zi, xr, yr, zr, off;
    static integer kgy, kgz, stop;
    extern /* Subroutine */ int fatal_(void);
    static integer start;
    static doublereal xsort[200000], ysort[200000], zsort[200000];
    extern /* Subroutine */ int imagen_(doublereal *, doublereal *, 
	    doublereal *);
    static logical repeat;
    extern /* Subroutine */ int lights_(doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___99 = { 0, 0, 0, fmt_30, 0 };



#define elst_ref(a_1,a_2) neigh_1.elst[(a_2)*1200 + a_1 - 1201]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




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




/*     transfer interaction site coordinates to sorting arrays */

    /* Parameter adjustments */
    --zold;
    --yold;
    --xold;

    /* Function Body */
    i__1 = charge_1.nion;
    for (i__ = 1; i__ <= i__1; ++i__) {
	neigh_1.nelst[i__ - 1] = 0;
	ic = charge_1.kion[i__ - 1];
	xold[i__] = atoms_1.x[ic - 1];
	yold[i__] = atoms_1.y[ic - 1];
	zold[i__] = atoms_1.z__[ic - 1];
	xsort[i__ - 1] = atoms_1.x[ic - 1];
	ysort[i__ - 1] = atoms_1.y[ic - 1];
	zsort[i__ - 1] = atoms_1.z__[ic - 1];
    }

/*     use the method of lights to generate neighbors */

    off = sqrt(neigh_1.cbuf2);
    lights_(&off, &charge_1.nion, xsort, ysort, zsort);

/*     loop over all atoms computing the interactions */

    i__1 = charge_1.nion;
    for (i__ = 1; i__ <= i__1; ++i__) {
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
L10:
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    k = light_1.locx[j - 1];
	    kgy = light_1.rgy[k - 1];
	    if (light_1.kby[i__ - 1] <= light_1.key[i__ - 1]) {
		if (kgy < light_1.kby[i__ - 1] || kgy > light_1.key[i__ - 1]) 
			{
		    goto L20;
		}
	    } else {
		if (kgy < light_1.kby[i__ - 1] && kgy > light_1.key[i__ - 1]) 
			{
		    goto L20;
		}
	    }
	    kgz = light_1.rgz[k - 1];
	    if (light_1.kbz[i__ - 1] <= light_1.kez[i__ - 1]) {
		if (kgz < light_1.kbz[i__ - 1] || kgz > light_1.kez[i__ - 1]) 
			{
		    goto L20;
		}
	    } else {
		if (kgz < light_1.kbz[i__ - 1] && kgz > light_1.kez[i__ - 1]) 
			{
		    goto L20;
		}
	    }
	    xr = xi - xsort[j - 1];
	    yr = yi - ysort[kgy - 1];
	    zr = zi - zsort[kgz - 1];
	    imagen_(&xr, &yr, &zr);
	    r2 = xr * xr + yr * yr + zr * zr;
	    if (r2 <= neigh_1.cbuf2) {
		if (i__ < k) {
		    ++neigh_1.nelst[i__ - 1];
		    elst_ref(neigh_1.nelst[i__ - 1], i__) = k;
		} else {
		    ++neigh_1.nelst[k - 1];
		    elst_ref(neigh_1.nelst[k - 1], k) = i__;
		}
	    }
L20:
	    ;
	}
	if (repeat) {
	    repeat = FALSE_;
	    start = light_1.kbx[i__ - 1] + 1;
	    stop = light_1.nlight;
	    goto L10;
	}
    }

/*     check to see if the neighbor lists are too long */

    i__1 = charge_1.nion;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (neigh_1.nelst[i__ - 1] >= 1200) {
	    io___99.ciunit = iounit_1.iout;
	    s_wsfe(&io___99);
	    e_wsfe();
	    fatal_();
	}
    }
    return 0;
} /* cfull_ */

#undef elst_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine mlist  --  build atom multipole neighbor lists  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "mlist" performs and update or a complete rebuild of the */
/*     electrostatic neighbor lists for atomic multipoles */


/* Subroutine */ int mlist_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 MLIST  --  Pairwise Neighbor List canno"
	    "t\002,\002 be used with Replicas\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, k;
    static doublereal r2;
    static integer ii, kk;
    static doublereal xi, yi, zi, xr, yr, zr, xold[25000], yold[25000], zold[
	    25000];
    extern /* Subroutine */ int fatal_(void), mfull_(doublereal *, doublereal 
	    *, doublereal *);
    static logical reset[25000];
    extern /* Subroutine */ int imagen_(doublereal *, doublereal *, 
	    doublereal *), mbuild_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal radius;
    extern /* Subroutine */ int replica_(doublereal *);

    /* Fortran I/O blocks */
    static cilist io___101 = { 0, 0, 0, fmt_10, 0 };




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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  boxes.i  --  parameters for periodic boundary conditions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     xbox        length of a-axis of periodic box in Angstroms */
/*     ybox        length of b-axis of periodic box in Angstroms */
/*     zbox        length of c-axis of periodic box in Angstroms */
/*     alpha       angle between b- and c-axes of box in degrees */
/*     beta        angle between a- and c-axes of box in degrees */
/*     gamma       angle between a- and b-axes of box in degrees */
/*     xbox2       half of the a-axis length of periodic box */
/*     ybox2       half of the b-axis length of periodic box */
/*     zbox2       half of the c-axis length of periodic box */
/*     box34       three-fourths axis length of truncated octahedron */
/*     lvec        real space lattice vectors as matrix rows */
/*     recip       reciprocal lattice vectors as matrix columns */
/*     volbox      volume in Ang**3 of the periodic box */
/*     beta_sin    sine of the beta periodic box angle */
/*     beta_cos    cosine of the beta periodic box angle */
/*     gamma_sin   sine of the gamma periodic box angle */
/*     gamma_cos   cosine of the gamma periodic box angle */
/*     beta_term   term used in generating triclinic box */
/*     gamma_term  term used in generating triclinic box */
/*     orthogonal  flag to mark periodic box as orthogonal */
/*     monoclinic  flag to mark periodic box as monoclinic */
/*     triclinic   flag to mark periodic box as triclinic */
/*     octahedron  flag to mark box as truncated octahedron */
/*     spacegrp    space group symbol for the unitcell type */




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




/*     neighbor list cannot be used with the replicates method */

    radius = sqrt(neigh_1.mbuf2);
    replica_(&radius);
    if (bound_1.use_replica__) {
	io___101.ciunit = iounit_1.iout;
	s_wsfe(&io___101);
	e_wsfe();
	fatal_();
    }

/*     perform a complete list build instead of an update */

    if (neigh_1.domlst) {
	neigh_1.domlst = FALSE_;
	if (boxes_1.octahedron) {
	    i__1 = mpole_1.npole;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		mbuild_(&i__, xold, yold, zold);
	    }
	} else {
	    mfull_(xold, yold, zold);
	}
	return 0;
    }

/*     update sites whose displacement exceeds half the buffer */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	reset[i__ - 1] = FALSE_;
	ii = mpole_1.ipole[i__ - 1];
	xi = atoms_1.x[ii - 1];
	yi = atoms_1.y[ii - 1];
	zi = atoms_1.z__[ii - 1];
	xr = xi - xold[i__ - 1];
	yr = yi - yold[i__ - 1];
	zr = zi - zold[i__ - 1];
	imagen_(&xr, &yr, &zr);
	r2 = xr * xr + yr * yr + zr * zr;
	if (neigh_1.domlst || r2 >= neigh_1.lbuf2) {
	    mbuild_(&i__, xold, yold, zold);
	    reset[i__ - 1] = TRUE_;
	    i__2 = i__ - 1;
	    for (k = 1; k <= i__2; ++k) {
		if (! reset[k - 1]) {
		    kk = mpole_1.ipole[k - 1];
		    xr = xi - atoms_1.x[kk - 1];
		    yr = yi - atoms_1.y[kk - 1];
		    zr = zi - atoms_1.z__[kk - 1];
		    imagen_(&xr, &yr, &zr);
		    r2 = xr * xr + yr * yr + zr * zr;

/*     rebuild the list of any possible lower numbered neighbor */

		    if (r2 <= neigh_1.mbuf2) {
			mbuild_(&k, xold, yold, zold);
			reset[k - 1] = TRUE_;
		    }
		}
	    }
	}
    }
    return 0;
} /* mlist_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine mbuild  --  make mpole pair list for one site  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "mbuild" performs a complete rebuild of the atomic multipole */
/*     electrostatic neighbor list for a single site */


/* Subroutine */ int mbuild_(integer *i__, doublereal *xold, doublereal *yold,
	 doublereal *zold)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 MBUILD  --  Too many Neighbors;\002,\002"
	    " Increase MAXELST\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer k;
    static doublereal r2;
    static integer ii, kk;
    static doublereal xi, yi, zi, xr, yr, zr;
    extern /* Subroutine */ int fatal_(void), imagen_(doublereal *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___127 = { 0, 0, 0, fmt_10, 0 };



#define elst_ref(a_1,a_2) neigh_1.elst[(a_2)*1200 + a_1 - 1201]



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




/*     zero out the list and get coordinates for the site */

    /* Parameter adjustments */
    --zold;
    --yold;
    --xold;

    /* Function Body */
    neigh_1.nelst[*i__ - 1] = 0;
    ii = mpole_1.ipole[*i__ - 1];
    xi = atoms_1.x[ii - 1];
    yi = atoms_1.y[ii - 1];
    zi = atoms_1.z__[ii - 1];

/*     generate all neighbors for the site being rebuilt */

    i__1 = mpole_1.npole;
    for (k = *i__ + 1; k <= i__1; ++k) {
	kk = mpole_1.ipole[k - 1];
	xr = xi - atoms_1.x[kk - 1];
	yr = yi - atoms_1.y[kk - 1];
	zr = zi - atoms_1.z__[kk - 1];
	imagen_(&xr, &yr, &zr);
	r2 = xr * xr + yr * yr + zr * zr;
	if (r2 <= neigh_1.mbuf2) {
	    ++neigh_1.nelst[*i__ - 1];
	    elst_ref(neigh_1.nelst[*i__ - 1], *i__) = k;
	}
    }

/*     store new coordinates to reflect update of the site */

    xold[*i__] = xi;
    yold[*i__] = yi;
    zold[*i__] = zi;

/*     check to see if the neighbor list is too long */

    if (neigh_1.nelst[*i__ - 1] >= 1200) {
	io___127.ciunit = iounit_1.iout;
	s_wsfe(&io___127);
	e_wsfe();
	fatal_();
    }
    return 0;
} /* mbuild_ */

#undef elst_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine mfull  --  make mpole pair list for all sites  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "mfull" performs a complete rebuild of the partial charge */
/*     pair neighbor list for all sites using the method of lights */


/* Subroutine */ int mfull_(doublereal *xold, doublereal *yold, doublereal *
	zold)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 MFULL  --  Too many Neighbors;\002,\002 "
	    "Increase MAXELST\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, j, k;
    static doublereal r2;
    static integer ii;
    static doublereal xi, yi, zi, xr, yr, zr, off;
    static integer kgy, kgz, stop;
    extern /* Subroutine */ int fatal_(void);
    static integer start;
    static doublereal xsort[200000], ysort[200000], zsort[200000];
    extern /* Subroutine */ int imagen_(doublereal *, doublereal *, 
	    doublereal *);
    static logical repeat;
    extern /* Subroutine */ int lights_(doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___148 = { 0, 0, 0, fmt_30, 0 };



#define elst_ref(a_1,a_2) neigh_1.elst[(a_2)*1200 + a_1 - 1201]



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




/*     transfer interaction site coordinates to sorting arrays */

    /* Parameter adjustments */
    --zold;
    --yold;
    --xold;

    /* Function Body */
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	neigh_1.nelst[i__ - 1] = 0;
	ii = mpole_1.ipole[i__ - 1];
	xold[i__] = atoms_1.x[ii - 1];
	yold[i__] = atoms_1.y[ii - 1];
	zold[i__] = atoms_1.z__[ii - 1];
	xsort[i__ - 1] = atoms_1.x[ii - 1];
	ysort[i__ - 1] = atoms_1.y[ii - 1];
	zsort[i__ - 1] = atoms_1.z__[ii - 1];
    }

/*     use the method of lights to generate neighbors */

    off = sqrt(neigh_1.mbuf2);
    lights_(&off, &mpole_1.npole, xsort, ysort, zsort);

/*     loop over all atoms computing the interactions */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
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
L10:
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    k = light_1.locx[j - 1];
	    kgy = light_1.rgy[k - 1];
	    if (light_1.kby[i__ - 1] <= light_1.key[i__ - 1]) {
		if (kgy < light_1.kby[i__ - 1] || kgy > light_1.key[i__ - 1]) 
			{
		    goto L20;
		}
	    } else {
		if (kgy < light_1.kby[i__ - 1] && kgy > light_1.key[i__ - 1]) 
			{
		    goto L20;
		}
	    }
	    kgz = light_1.rgz[k - 1];
	    if (light_1.kbz[i__ - 1] <= light_1.kez[i__ - 1]) {
		if (kgz < light_1.kbz[i__ - 1] || kgz > light_1.kez[i__ - 1]) 
			{
		    goto L20;
		}
	    } else {
		if (kgz < light_1.kbz[i__ - 1] && kgz > light_1.kez[i__ - 1]) 
			{
		    goto L20;
		}
	    }
	    xr = xi - xsort[j - 1];
	    yr = yi - ysort[kgy - 1];
	    zr = zi - zsort[kgz - 1];
	    imagen_(&xr, &yr, &zr);
	    r2 = xr * xr + yr * yr + zr * zr;
	    if (r2 <= neigh_1.mbuf2) {
		if (i__ < k) {
		    ++neigh_1.nelst[i__ - 1];
		    elst_ref(neigh_1.nelst[i__ - 1], i__) = k;
		} else {
		    ++neigh_1.nelst[k - 1];
		    elst_ref(neigh_1.nelst[k - 1], k) = i__;
		}
	    }
L20:
	    ;
	}
	if (repeat) {
	    repeat = FALSE_;
	    start = light_1.kbx[i__ - 1] + 1;
	    stop = light_1.nlight;
	    goto L10;
	}
    }

/*     check to see if the neighbor lists are too long */

    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (neigh_1.nelst[i__ - 1] >= 1200) {
	    io___148.ciunit = iounit_1.iout;
	    s_wsfe(&io___148);
	    e_wsfe();
	    fatal_();
	}
    }
    return 0;
} /* mfull_ */

#undef elst_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine imagen  --  neighbor minimum image distance  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "imagen" takes the components of pairwise distance between */
/*     two points and converts to the components of the minimum */
/*     image distance; fast version for neighbor list generation */


/* Subroutine */ int imagen_(doublereal *xr, doublereal *yr, doublereal *zr)
{
    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  boxes.i  --  parameters for periodic boundary conditions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     xbox        length of a-axis of periodic box in Angstroms */
/*     ybox        length of b-axis of periodic box in Angstroms */
/*     zbox        length of c-axis of periodic box in Angstroms */
/*     alpha       angle between b- and c-axes of box in degrees */
/*     beta        angle between a- and c-axes of box in degrees */
/*     gamma       angle between a- and b-axes of box in degrees */
/*     xbox2       half of the a-axis length of periodic box */
/*     ybox2       half of the b-axis length of periodic box */
/*     zbox2       half of the c-axis length of periodic box */
/*     box34       three-fourths axis length of truncated octahedron */
/*     lvec        real space lattice vectors as matrix rows */
/*     recip       reciprocal lattice vectors as matrix columns */
/*     volbox      volume in Ang**3 of the periodic box */
/*     beta_sin    sine of the beta periodic box angle */
/*     beta_cos    cosine of the beta periodic box angle */
/*     gamma_sin   sine of the gamma periodic box angle */
/*     gamma_cos   cosine of the gamma periodic box angle */
/*     beta_term   term used in generating triclinic box */
/*     gamma_term  term used in generating triclinic box */
/*     orthogonal  flag to mark periodic box as orthogonal */
/*     monoclinic  flag to mark periodic box as monoclinic */
/*     triclinic   flag to mark periodic box as triclinic */
/*     octahedron  flag to mark box as truncated octahedron */
/*     spacegrp    space group symbol for the unitcell type */




/*     for orthogonal lattice, find the desired image directly */

    if (boxes_1.orthogonal) {
	*xr = abs(*xr);
	*yr = abs(*yr);
	*zr = abs(*zr);
	if (*xr > boxes_1.xbox2) {
	    *xr -= boxes_1.xbox;
	}
	if (*yr > boxes_1.ybox2) {
	    *yr -= boxes_1.ybox;
	}
	if (*zr > boxes_1.zbox2) {
	    *zr -= boxes_1.zbox;
	}

/*     for monoclinic lattice, convert "xr" and "zr" specially */

    } else if (boxes_1.monoclinic) {
	if (abs(*xr) > boxes_1.xbox2) {
	    *xr -= d_sign(&boxes_1.xbox, xr);
	}
	if (abs(*yr) > boxes_1.ybox2) {
	    *yr -= d_sign(&boxes_1.ybox, yr);
	}
	if (abs(*zr) > boxes_1.zbox2) {
	    *zr -= d_sign(&boxes_1.zbox, zr);
	}
	*xr += *zr * boxes_1.beta_cos__;
	*zr *= boxes_1.beta_sin__;

/*     for triclinic lattice, use general conversion equations */

    } else if (boxes_1.triclinic) {
	if (abs(*xr) > boxes_1.xbox2) {
	    *xr -= d_sign(&boxes_1.xbox, xr);
	}
	if (abs(*yr) > boxes_1.ybox2) {
	    *yr -= d_sign(&boxes_1.ybox, yr);
	}
	if (abs(*zr) > boxes_1.zbox2) {
	    *zr -= d_sign(&boxes_1.zbox, zr);
	}
	*xr = *xr + *yr * boxes_1.gamma_cos__ + *zr * boxes_1.beta_cos__;
	*yr = *yr * boxes_1.gamma_sin__ + *zr * boxes_1.beta_term__;
	*zr *= boxes_1.gamma_term__;

/*     for truncated octahedron, remove the corner pieces */

    } else if (boxes_1.octahedron) {
	if (abs(*xr) > boxes_1.xbox2) {
	    *xr -= d_sign(&boxes_1.xbox, xr);
	}
	if (abs(*yr) > boxes_1.ybox2) {
	    *yr -= d_sign(&boxes_1.ybox, yr);
	}
	if (abs(*zr) > boxes_1.zbox2) {
	    *zr -= d_sign(&boxes_1.zbox, zr);
	}
	if (abs(*xr) + abs(*yr) + abs(*zr) > boxes_1.box34) {
	    *xr -= d_sign(&boxes_1.xbox2, xr);
	    *yr -= d_sign(&boxes_1.ybox2, yr);
	    *zr -= d_sign(&boxes_1.zbox2, zr);
	}
    }
    return 0;
} /* imagen_ */

