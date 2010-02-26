/* sktstuff.f -- translated by f2c (version 20050501).
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
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

struct {
    doublereal polarity[25000], thole[25000], pdamp[25000], uind[75000]	/* 
	    was [3][25000] */, uinp[75000]	/* was [3][25000] */, uinds[
	    75000]	/* was [3][25000] */, uinps[75000]	/* was [3][
	    25000] */;
    integer npolar;
} polar_;

#define polar_1 polar_

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
    integer runtyp, cstep;
    doublereal cdt, cenergy, cdx[25000], cdy[25000], cdz[25000];
    logical use_socket__, skt_init__, skt_close__;
} socket_;

#define socket_1 socket_

struct {
    doublereal v[75000]	/* was [3][25000] */, a[75000]	/* was [3][25000] */, 
	    aold[75000]	/* was [3][25000] */;
} moldyn_;

#define moldyn_1 moldyn_

struct {
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

struct {
    doublereal pchg[25000];
    integer nion, iion[25000], jion[25000], kion[25000], chglist[25000];
} charge_;

#define charge_1 charge_

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
    integer biotyp[10000];
    char forcefield[20];
} fields_;

#define fields_1 fields_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_



/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2002 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine sktopt  --  send current optimization info  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "sktopt" sends the current optimization info via a socket */


/* Subroutine */ int sktopt_(integer *ncycle, doublereal *eopt)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k;
    extern /* Subroutine */ int setenergy_(doublereal *);
    static doublereal px[25000], py[25000], pz[25000];
    extern /* Subroutine */ int needupdate_(integer *), setinduced_(integer *,
	     doublereal *, doublereal *, doublereal *), setupdated_(void), 
	    getmonitor_(void);
    static integer flag__;
    extern /* Subroutine */ int setgradients_(integer *, doublereal *, 
	    doublereal *, doublereal *), setcoordinates_(integer *, 
	    doublereal *, doublereal *, doublereal *), releasemonitor_(void), 
	    sktinit_(void), setstep_(integer *);


#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]
#define desum_ref(a_1,a_2) deriv_1.desum[(a_2)*3 + a_1 - 4]



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  polar.i  --  polarizabilities and induced dipole moments  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     polarity  dipole polarizability for each multipole site (Ang**3) */
/*     thole     Thole polarizability damping value for each site */
/*     pdamp     value of polarizability scale factor for each site */
/*     uind      induced dipole components at each multipole site */
/*     uinp      induced dipoles in field used for energy interactions */
/*     uinds     GK or PB induced dipoles at each multipole site */
/*     uinps     induced dipoles in field used for GK or PB energy */
/*     npolar    total number of polarizable sites in the system */




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
/*     ##  COPYRIGHT (C)  2002  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  socket.i  --  control parameters for socket communication  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     runtyp      calculation type for passing socket information */
/*     cstep       current optimization or dynamics step number */
/*     cdt         current dynamics cumulative simulation time */
/*     cenergy     current potential energy from simulation */
/*     cdx         current gradient components along the x-axis */
/*     cdy         current gradient components along the y-axis */
/*     cdz         current gradient components along the z-axis */
/*     use_socket  logical flag governing use of external sockets */
/*     skt_init    logical flag to indicate socket initialization */
/*     skt_close   logical flag to indicate socket shutdown */




/*     check to see if the Server has been created */

    socket_1.runtyp = 2;
    if (! socket_1.skt_init__) {
	sktinit_();
    }
    if (! socket_1.use_socket__) {
	return 0;
    }

/*     save the current step number and energy */

    socket_1.cstep = *ncycle;
    socket_1.cenergy = *eopt;

/*     check to see if an update is needed */

    flag__ = 1;
    if (! socket_1.skt_close__) {
	needupdate_(&flag__);
    }
    if (flag__ == 0) {
	return 0;
    }

/*     get the monitor for the update structure */

    getmonitor_();

/*     load current optimization data */

    setcoordinates_(&atoms_1.n, atoms_1.x, atoms_1.y, atoms_1.z__);
    setstep_(ncycle);
    setenergy_(eopt);
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	socket_1.cdx[i__ - 1] = desum_ref(1, i__);
	socket_1.cdy[i__ - 1] = desum_ref(2, i__);
	socket_1.cdz[i__ - 1] = desum_ref(3, i__);
	px[i__ - 1] = 0.;
	py[i__ - 1] = 0.;
	pz[i__ - 1] = 0.;
    }
    setgradients_(&atoms_1.n, socket_1.cdx, socket_1.cdy, socket_1.cdz);
    if (potent_1.use_polar__) {
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = mpole_1.ipole[i__ - 1];
	    px[k - 1] = uind_ref(1, i__);
	    py[k - 1] = uind_ref(2, i__);
	    pz[k - 1] = uind_ref(3, i__);
	}
	setinduced_(&atoms_1.n, px, py, pz);
    }

/*     release the monitor for the system stucture */

    setupdated_();
    releasemonitor_();
    return 0;
} /* sktopt_ */

#undef desum_ref
#undef uind_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine sktdyn  --   send the current dynamics info  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "sktdyn" sends the current dynamics info via a socket */


/* Subroutine */ int sktdyn_(integer *istep, doublereal *dt, doublereal *epot)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k;
    extern /* Subroutine */ int setenergy_(doublereal *);
    static doublereal ax[25000], ay[25000], az[25000], px[25000], py[25000], 
	    pz[25000], vx[25000], vy[25000], vz[25000];
    extern /* Subroutine */ int needupdate_(integer *), setinduced_(integer *,
	     doublereal *, doublereal *, doublereal *), setupdated_(void), 
	    getmonitor_(void), setvelocity_(integer *, doublereal *, 
	    doublereal *, doublereal *);
    static integer flag__;
    static doublereal time;
    extern /* Subroutine */ int setcoordinates_(integer *, doublereal *, 
	    doublereal *, doublereal *), releasemonitor_(void), settime_(
	    doublereal *), sktinit_(void), setacceleration_(integer *, 
	    doublereal *, doublereal *, doublereal *);


#define a_ref(a_1,a_2) moldyn_1.a[(a_2)*3 + a_1 - 4]
#define v_ref(a_1,a_2) moldyn_1.v[(a_2)*3 + a_1 - 4]
#define uind_ref(a_1,a_2) polar_1.uind[(a_2)*3 + a_1 - 4]



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  moldyn.i  --  velocity and acceleration on MD trajectory  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     v       current velocity of each atom along the x,y,z-axes */
/*     a       current acceleration of each atom along x,y,z-axes */
/*     aold    previous acceleration of each atom along x,y,z-axes */




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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  polar.i  --  polarizabilities and induced dipole moments  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     polarity  dipole polarizability for each multipole site (Ang**3) */
/*     thole     Thole polarizability damping value for each site */
/*     pdamp     value of polarizability scale factor for each site */
/*     uind      induced dipole components at each multipole site */
/*     uinp      induced dipoles in field used for energy interactions */
/*     uinds     GK or PB induced dipoles at each multipole site */
/*     uinps     induced dipoles in field used for GK or PB energy */
/*     npolar    total number of polarizable sites in the system */




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
/*     ##  COPYRIGHT (C)  2002  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  socket.i  --  control parameters for socket communication  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     runtyp      calculation type for passing socket information */
/*     cstep       current optimization or dynamics step number */
/*     cdt         current dynamics cumulative simulation time */
/*     cenergy     current potential energy from simulation */
/*     cdx         current gradient components along the x-axis */
/*     cdy         current gradient components along the y-axis */
/*     cdz         current gradient components along the z-axis */
/*     use_socket  logical flag governing use of external sockets */
/*     skt_init    logical flag to indicate socket initialization */
/*     skt_close   logical flag to indicate socket shutdown */




/*     check to see if the Java objects have been created */

    socket_1.runtyp = 1;
    if (! socket_1.skt_init__) {
	sktinit_();
    }
    if (! socket_1.use_socket__) {
	return 0;
    }

/*     save the current step number, time and energy */

    socket_1.cstep = *istep;
    socket_1.cdt = *dt;
    socket_1.cenergy = *epot;

/*     check to see if we need to update the system info */

    flag__ = 1;
    if (! socket_1.skt_close__) {
	needupdate_(&flag__);
    }
    if (flag__ == 0) {
	return 0;
    }

/*     get the monitor for the update structure */

    getmonitor_();

/*     load the current dynamics information */

    setcoordinates_(&atoms_1.n, atoms_1.x, atoms_1.y, atoms_1.z__);
    time = (doublereal) (*istep) * *dt;
    settime_(&time);
    setenergy_(epot);
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vx[i__ - 1] = v_ref(1, i__);
	vy[i__ - 1] = v_ref(2, i__);
	vz[i__ - 1] = v_ref(3, i__);
	ax[i__ - 1] = a_ref(1, i__);
	ay[i__ - 1] = a_ref(2, i__);
	az[i__ - 1] = a_ref(3, i__);
	px[i__ - 1] = 0.;
	py[i__ - 1] = 0.;
	pz[i__ - 1] = 0.;
    }
    setvelocity_(&atoms_1.n, vx, vy, vz);
    setacceleration_(&atoms_1.n, ax, ay, az);
    if (potent_1.use_polar__) {
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = mpole_1.ipole[i__ - 1];
	    px[k - 1] = uind_ref(1, i__);
	    py[k - 1] = uind_ref(2, i__);
	    pz[k - 1] = uind_ref(3, i__);
	}
	setinduced_(&atoms_1.n, px, py, pz);
    }

/*     release the monitor for the update stucture */

    setupdated_();
    releasemonitor_();
    return 0;
} /* sktdyn_ */

#undef uind_ref
#undef v_ref
#undef a_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine sktinit  --  initialize socket communication  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "sktinit" sets up socket communication with the graphical */
/*     user interface by starting a Java virtual machine, initiating */
/*     a server, and loading an object with system information */


/* Subroutine */ int sktinit_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 SKTINIT  --  Unable to Start Server fo"
	    "r\002,\002 Java GUI Communication\002,/,\002 Check the LD_LIBRAR"
	    "Y_PATH and CLASSPATH\002,\002 Environment Variables\002,/)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int setconnectivity_(integer *, integer *, 
	    integer *, integer *, integer *);
    static integer i__;
    extern /* Subroutine */ int setstory_(integer *, char *, ftnlen);
    static integer b1[25000], b2[25000], b3[25000], b4[25000];
    extern /* Subroutine */ int setcharge_(integer *, doublereal *), 
	    chksocket_(integer *), createjvm_(integer *), setatomic_(integer *
	    , integer *), setkeyword_(integer *, char *, ftnlen);
    static integer flag__;
    extern /* Subroutine */ int createupdate_(integer *, integer *, integer *,
	     integer *), createserver_(integer *), createsystem_(integer *, 
	    integer *, integer *), setatomtypes_(integer *, integer *), 
	    setforcefield_(char *, ftnlen), setcoordinates_(integer *, 
	    doublereal *, doublereal *, doublereal *), setfile_(char *, 
	    ftnlen), setname_(integer *, char *, ftnlen), setmass_(integer *, 
	    doublereal *);

    /* Fortran I/O blocks */
    static cilist io___21 = { 0, 0, 0, fmt_10, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define story_ref(a_0,a_1) &atmtyp_1.story[(a_1)*24 + a_0 - 24]
#define keyline_ref(a_0,a_1) &keys_1.keyline[(a_1)*120 + a_0 - 120]



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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  fields.i  --  molecular mechanics force field description  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     biotyp       force field atom type of each biopolymer type */
/*     forcefield   string used to describe the current forcefield */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  keys.i  --  contents of current keyword parameter file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nkey      number of nonblank lines in the keyword file */
/*     keyline   contents of each individual keyword file line */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  polar.i  --  polarizabilities and induced dipole moments  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     polarity  dipole polarizability for each multipole site (Ang**3) */
/*     thole     Thole polarizability damping value for each site */
/*     pdamp     value of polarizability scale factor for each site */
/*     uind      induced dipole components at each multipole site */
/*     uinp      induced dipoles in field used for energy interactions */
/*     uinds     GK or PB induced dipoles at each multipole site */
/*     uinps     induced dipoles in field used for GK or PB energy */
/*     npolar    total number of polarizable sites in the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2002  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  socket.i  --  control parameters for socket communication  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     runtyp      calculation type for passing socket information */
/*     cstep       current optimization or dynamics step number */
/*     cdt         current dynamics cumulative simulation time */
/*     cenergy     current potential energy from simulation */
/*     cdx         current gradient components along the x-axis */
/*     cdy         current gradient components along the y-axis */
/*     cdz         current gradient components along the z-axis */
/*     use_socket  logical flag governing use of external sockets */
/*     skt_init    logical flag to indicate socket initialization */
/*     skt_close   logical flag to indicate socket shutdown */




/*     set initialization flag and test for socket usage */

    socket_1.skt_init__ = TRUE_;
    socket_1.use_socket__ = TRUE_;
    chksocket_(&flag__);
    if (flag__ == 0) {
	socket_1.use_socket__ = FALSE_;
	return 0;
    }

/*     create the Java Virtual Machine */

    createjvm_(&flag__);
    if (flag__ == 0) {
	socket_1.use_socket__ = FALSE_;
	if (inform_1.debug) {
	    io___21.ciunit = iounit_1.iout;
	    s_wsfe(&io___21);
	    e_wsfe();
	}
	return 0;
    }

/*     create the TINKER system object */

    createsystem_(&atoms_1.n, &keys_1.nkey, &flag__);
    if (flag__ == 0) {
	socket_1.use_socket__ = FALSE_;
	return 0;
    }

/*     load the coordinates and keyfile information */

    setfile_(files_1.filename, (ftnlen)120);
    setforcefield_(fields_1.forcefield, (ftnlen)20);
    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	setkeyword_(&i__, keyline_ref(0, i__), (ftnlen)120);
    }
    setcoordinates_(&atoms_1.n, atoms_1.x, atoms_1.y, atoms_1.z__);
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b1[i__ - 1] = i12_ref(1, i__);
	b2[i__ - 1] = i12_ref(2, i__);
	b3[i__ - 1] = i12_ref(3, i__);
	b4[i__ - 1] = i12_ref(4, i__);
    }
    setconnectivity_(&atoms_1.n, b1, b2, b3, b4);

/*     load atom type information for the parameter set */

    setatomtypes_(&atoms_1.n, atoms_1.type__);
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	setname_(&i__, name___ref(0, i__), (ftnlen)3);
	setstory_(&i__, story_ref(0, i__), (ftnlen)24);
    }
    setatomic_(&atoms_1.n, atmtyp_1.atomic);
    setmass_(&atoms_1.n, atmtyp_1.mass);
    setcharge_(&atoms_1.n, charge_1.pchg);

/*     create the TINKER server */

    createserver_(&flag__);
    if (flag__ == 0) {
	socket_1.use_socket__ = FALSE_;
	return 0;
    }

/*     create the update object */

    createupdate_(&atoms_1.n, &socket_1.runtyp, &polar_1.npolar, &flag__);
    if (flag__ == 0) {
	socket_1.use_socket__ = FALSE_;
	return 0;
    }
    return 0;
} /* sktinit_ */

#undef keyline_ref
#undef story_ref
#undef name___ref
#undef i12_ref




/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine sktkill  --  shutdown the server and JVM  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "sktkill" closes the server and Java virtual machine */


/* Subroutine */ int sktkill_(void)
{
    extern /* Subroutine */ int destroyjvm_(void), destroyserver_(void), 
	    sktdyn_(integer *, doublereal *, doublereal *), sktopt_(integer *,
	     doublereal *);



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




/*     check to see if there is anything to close */



/*     ################################################### */
/*     ##  COPYRIGHT (C)  2002  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  socket.i  --  control parameters for socket communication  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     runtyp      calculation type for passing socket information */
/*     cstep       current optimization or dynamics step number */
/*     cdt         current dynamics cumulative simulation time */
/*     cenergy     current potential energy from simulation */
/*     cdx         current gradient components along the x-axis */
/*     cdy         current gradient components along the y-axis */
/*     cdz         current gradient components along the z-axis */
/*     use_socket  logical flag governing use of external sockets */
/*     skt_init    logical flag to indicate socket initialization */
/*     skt_close   logical flag to indicate socket shutdown */


    if (! socket_1.use_socket__) {
	return 0;
    }
    socket_1.skt_close__ = TRUE_;

/*     load the final simulation results */

    if (socket_1.runtyp == 1) {
	sktdyn_(&socket_1.cstep, &socket_1.cdt, &socket_1.cenergy);
    }
    if (socket_1.runtyp == 2) {
	sktopt_(&socket_1.cstep, &socket_1.cenergy);
    }

/*     shutdown the TINKER server */

    destroyserver_();

/*     shutdown the Java virtual machine */

    destroyjvm_();
    return 0;
} /* sktkill_ */

