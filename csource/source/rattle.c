/* rattle.f -- translated by f2c (version 20050501).
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
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

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
    doublereal v[75000]	/* was [3][25000] */, a[75000]	/* was [3][25000] */, 
	    aold[75000]	/* was [3][25000] */;
} moldyn_;

#define moldyn_1 moldyn_

struct {
    doublereal krat[25000];
    integer nrat, nratx, irat[50000]	/* was [2][25000] */, iratx[25000], 
	    kratx[25000];
    logical ratimage[25000], use_rattle__;
} shake_;

#define shake_1 shake_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal vir[9]	/* was [3][3] */;
} virial_;

#define virial_1 virial_

/* Table of constant values */

static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine rattle  --  distance and spatial constraints  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "rattle" implements the first portion of the rattle algorithm */
/*     by correcting atomic positions and half-step velocities to */
/*     maintain interatomic distance and absolute spatial constraints */

/*     literature reference: */

/*     H. C. Andersen, "RATTLE: A Velocity Version of the SHAKE */
/*     Algorithm for Molecular Dynamics Calculations", Journal of */
/*     Computational Physics, 52, 24-34 (1983) */


/* Subroutine */ int rattle_(doublereal *dt, doublereal *xold, doublereal *
	yold, doublereal *zold)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 RATTLE  --  Warning, Distance Constrai"
	    "nts\002,\002 not Satisfied\002)";
    static char fmt_20[] = "(\002 RATTLE   --  Distance Constraints met a"
	    "t\002,i6,\002 Iterations\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, ia, ib;
    static doublereal xo, yo, zo, xr, yr, zr, xv, yv, zv, rma, rmb, dot, eps, 
	    sor;
    static integer mode;
    static logical done;
    static doublereal term;
    static integer stop;
    static doublereal dist2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *), fatal_(void);
    static doublereal delta, weigh;
    static logical moved[25000];
    static integer niter, start;
    static doublereal xterm, yterm, zterm;
    static logical update[25000];
    extern /* Subroutine */ int prterr_(void);
    static integer maxiter;

    /* Fortran I/O blocks */
    static cilist io___26 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_20, 0 };



#define v_ref(a_1,a_2) moldyn_1.v[(a_2)*3 + a_1 - 4]
#define irat_ref(a_1,a_2) shake_1.irat[(a_2)*2 + a_1 - 3]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  shake.i  --  definition of Shake/Rattle constraints  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     krat         ideal distance value for rattle constraint */
/*     nrat         number of rattle distance constraints to apply */
/*     nratx        number of atom group spatial constraints to apply */
/*     irat         atom numbers of atoms in a rattle constraint */
/*     iratx        group number of group in a spatial constraint */
/*     kratx        spatial constraint type (1=plane, 2=line, 3=point) */
/*     ratimage     flag to use minimum image for rattle constraint */
/*     use_rattle   logical flag to set use of rattle contraints */




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




/*     initialize the lists of atoms previously corrected */

    /* Parameter adjustments */
    --zold;
    --yold;
    --xold;

    /* Function Body */
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    moved[i__ - 1] = TRUE_;
	} else {
	    moved[i__ - 1] = FALSE_;
	}
	update[i__ - 1] = FALSE_;
    }

/*     set the iteration counter, termination and tolerance */

    maxiter = 100;
    sor = 1.25;
    eps = 1e-6;

/*     apply rattle to distances and half-step velocity values */

    niter = 0;
    done = FALSE_;
    while(! done && niter < maxiter) {
	++niter;
	done = TRUE_;
	i__1 = shake_1.nrat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = irat_ref(1, i__);
	    ib = irat_ref(2, i__);
	    if (moved[ia - 1] || moved[ib - 1]) {
		xr = atoms_1.x[ib - 1] - atoms_1.x[ia - 1];
		yr = atoms_1.y[ib - 1] - atoms_1.y[ia - 1];
		zr = atoms_1.z__[ib - 1] - atoms_1.z__[ia - 1];
		if (shake_1.ratimage[i__ - 1]) {
		    image_(&xr, &yr, &zr);
		}
/* Computing 2nd power */
		d__1 = xr;
/* Computing 2nd power */
		d__2 = yr;
/* Computing 2nd power */
		d__3 = zr;
		dist2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
/* Computing 2nd power */
		d__1 = shake_1.krat[i__ - 1];
		delta = d__1 * d__1 - dist2;
		if (abs(delta) > eps) {
		    done = FALSE_;
		    update[ia - 1] = TRUE_;
		    update[ib - 1] = TRUE_;
		    xo = xold[ib] - xold[ia];
		    yo = yold[ib] - yold[ia];
		    zo = zold[ib] - zold[ia];
		    if (shake_1.ratimage[i__ - 1]) {
			image_(&xo, &yo, &zo);
		    }
		    dot = xr * xo + yr * yo + zr * zo;
		    rma = 1. / atmtyp_1.mass[ia - 1];
		    rmb = 1. / atmtyp_1.mass[ib - 1];
		    term = sor * delta / ((rma + rmb) * 2. * dot);
		    xterm = xo * term;
		    yterm = yo * term;
		    zterm = zo * term;
		    atoms_1.x[ia - 1] -= xterm * rma;
		    atoms_1.y[ia - 1] -= yterm * rma;
		    atoms_1.z__[ia - 1] -= zterm * rma;
		    atoms_1.x[ib - 1] += xterm * rmb;
		    atoms_1.y[ib - 1] += yterm * rmb;
		    atoms_1.z__[ib - 1] += zterm * rmb;
		    rma /= *dt;
		    rmb /= *dt;
		    v_ref(1, ia) = v_ref(1, ia) - xterm * rma;
		    v_ref(2, ia) = v_ref(2, ia) - yterm * rma;
		    v_ref(3, ia) = v_ref(3, ia) - zterm * rma;
		    v_ref(1, ib) = v_ref(1, ib) + xterm * rmb;
		    v_ref(2, ib) = v_ref(2, ib) + yterm * rmb;
		    v_ref(3, ib) = v_ref(3, ib) + zterm * rmb;
		}
	    }
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    moved[i__ - 1] = update[i__ - 1];
	    update[i__ - 1] = FALSE_;
	}
    }

/*     write information on the number of iterations needed */

    if (niter == maxiter) {
	io___26.ciunit = iounit_1.iout;
	s_wsfe(&io___26);
	e_wsfe();
	prterr_();
	fatal_();
    } else if (inform_1.debug) {
	io___27.ciunit = iounit_1.iout;
	s_wsfe(&io___27);
	do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	e_wsfe();
    }

/*     apply group position and velocity constraints via exact reset */

    i__1 = shake_1.nratx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = shake_1.iratx[i__ - 1];
	mode = shake_1.kratx[i__ - 1];
	xr = 0.;
	yr = 0.;
	zr = 0.;
	xv = 0.;
	yv = 0.;
	zv = 0.;
	start = igrp_ref(1, ia);
	stop = igrp_ref(2, ia);
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    weigh = atmtyp_1.mass[k - 1] / group_1.grpmass[ia - 1];
	    if (mode > 2) {
		xr += atoms_1.x[k - 1] * weigh;
		xv += v_ref(1, k) * weigh;
	    }
	    if (mode > 1) {
		yr += atoms_1.y[k - 1] * weigh;
		yv += v_ref(2, k) * weigh;
	    }
	    zr += atoms_1.z__[k - 1] * weigh;
	    zv += v_ref(3, k) * weigh;
	}
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    atoms_1.x[k - 1] -= xr;
	    atoms_1.y[k - 1] -= yr;
	    atoms_1.z__[k - 1] -= zr;
	    v_ref(1, k) = v_ref(1, k) - xv;
	    v_ref(2, k) = v_ref(2, k) - yv;
	    v_ref(3, k) = v_ref(3, k) - zv;
	}
    }
    return 0;
} /* rattle_ */

#undef igrp_ref
#undef irat_ref
#undef v_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine rattle2  --  rattle atom velocity constraints  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "rattle2" implements the second portion of the rattle algorithm */
/*     by correcting the full-step velocities in order to maintain */
/*     interatomic distance constraints */


/* Subroutine */ int rattle2_(doublereal *dt)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 RATTLE2  --  Warning, Velocity Constrain"
	    "ts\002,\002 not Satisfied\002)";
    static char fmt_20[] = "(\002 RATTLE2  --  Velocity Constraints met a"
	    "t\002,i6,\002 Iterations\002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, ia, ib;
    static doublereal xr, yr, zr, xv, yv, zv, rma, rmb, dot, eps, sor, vxx, 
	    vyx, vyy, vzx, vzz, vzy;
    static integer mode;
    static logical done;
    static doublereal term;
    static integer stop;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *), fatal_(void);
    static doublereal weigh;
    static logical moved[25000];
    static integer niter, start;
    static doublereal vterm, xterm, yterm, zterm;
    static logical update[25000];
    extern /* Subroutine */ int prterr_(void);
    static integer maxiter;

    /* Fortran I/O blocks */
    static cilist io___67 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___68 = { 0, 0, 0, fmt_20, 0 };



#define v_ref(a_1,a_2) moldyn_1.v[(a_2)*3 + a_1 - 4]
#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]
#define irat_ref(a_1,a_2) shake_1.irat[(a_2)*2 + a_1 - 3]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  shake.i  --  definition of Shake/Rattle constraints  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     krat         ideal distance value for rattle constraint */
/*     nrat         number of rattle distance constraints to apply */
/*     nratx        number of atom group spatial constraints to apply */
/*     irat         atom numbers of atoms in a rattle constraint */
/*     iratx        group number of group in a spatial constraint */
/*     kratx        spatial constraint type (1=plane, 2=line, 3=point) */
/*     ratimage     flag to use minimum image for rattle constraint */
/*     use_rattle   logical flag to set use of rattle contraints */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  units.i  --  physical constants and unit conversions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     literature reference: */

/*     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended */
/*     Values of the Fundamental Physical Constants: 2006", Reviews of */
/*     Modern Physics, 80, 633-730 (2008) */

/*     The "2006 CODATA Recommended Values" are also available from */
/*     the NIST Reference on Constants, Units, and Uncertainty site */
/*     at http://physics.nist.gov/cuu/index.html */

/*     Most values below are derived from the 2006 CODATA reference */

/*     The conversion from calorie to Joule is the definition of the */
/*     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992) */

/*     The "coulomb" energy conversion factor is found by dimensional */
/*     analysis of Coulomb's Law, ie, by dividing the square of the */
/*     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is */
/*     the permittivity of vacuum (the "electric constant"); note that */
/*     eps0 is typically given in F/m, equivalent to C**2/(J-m) */

/*     The approximate value used for the Debye, 3.33564 x 10-30 C-m, */
/*     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997) */

/*     The value of "prescon" is based on definition of 1 atmosphere */
/*     as 101325 Pa set by the 10th Conference Generale des Poids et */
/*     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3 */

/*     avogadro    Avogadro's number (N) in particles/mole */
/*     lightspd    speed of light in vacuum (c) in cm/ps */
/*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K */
/*     gasconst    ideal gas constant (R) in kcal/mole/K */
/*     emass       mass of an electron in atomic mass units */
/*     joule       conversion from calories to joules */
/*     convert     conversion from kcal to g*Ang**2/ps**2 */
/*     bohr        conversion from Bohrs to Angstroms */
/*     hartree     conversion from Hartree to kcal/mole */
/*     evolt       conversion from Hartree to electron-volts */
/*     efreq       conversion from Hartree to cm-1 */
/*     coulomb     conversion from electron**2/Ang to kcal/mole */
/*     debye       conversion from electron-Ang to Debyes */
/*     prescon     conversion from kcal/mole/Ang**3 to Atm */




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

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  virial.i  --  components of internal virial tensor  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     vir    total internal virial Cartesian tensor components */




/*     initialize the lists of atoms previously corrected */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    moved[i__ - 1] = TRUE_;
	} else {
	    moved[i__ - 1] = FALSE_;
	}
	update[i__ - 1] = FALSE_;
    }

/*     set the iteration counter, termination and tolerance */

    maxiter = 100;
    niter = 0;
    done = FALSE_;
    sor = 1.25;
    eps = 1e-6 / *dt;
    vterm = 2. / (*dt * 418.4);

/*     apply the rattle algorithm to correct the velocities */

    while(! done && niter < maxiter) {
	++niter;
	done = TRUE_;
	i__1 = shake_1.nrat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = irat_ref(1, i__);
	    ib = irat_ref(2, i__);
	    if (moved[ia - 1] || moved[ib - 1]) {
		xr = atoms_1.x[ib - 1] - atoms_1.x[ia - 1];
		yr = atoms_1.y[ib - 1] - atoms_1.y[ia - 1];
		zr = atoms_1.z__[ib - 1] - atoms_1.z__[ia - 1];
		if (shake_1.ratimage[i__ - 1]) {
		    image_(&xr, &yr, &zr);
		}
		xv = v_ref(1, ib) - v_ref(1, ia);
		yv = v_ref(2, ib) - v_ref(2, ia);
		zv = v_ref(3, ib) - v_ref(3, ia);
		dot = xr * xv + yr * yv + zr * zv;
		rma = 1. / atmtyp_1.mass[ia - 1];
		rmb = 1. / atmtyp_1.mass[ib - 1];
/* Computing 2nd power */
		d__1 = shake_1.krat[i__ - 1];
		term = -dot / ((rma + rmb) * (d__1 * d__1));
		if (abs(term) > eps) {
		    done = FALSE_;
		    update[ia - 1] = TRUE_;
		    update[ib - 1] = TRUE_;
		    term = sor * term;
		    xterm = xr * term;
		    yterm = yr * term;
		    zterm = zr * term;
		    v_ref(1, ia) = v_ref(1, ia) - xterm * rma;
		    v_ref(2, ia) = v_ref(2, ia) - yterm * rma;
		    v_ref(3, ia) = v_ref(3, ia) - zterm * rma;
		    v_ref(1, ib) = v_ref(1, ib) + xterm * rmb;
		    v_ref(2, ib) = v_ref(2, ib) + yterm * rmb;
		    v_ref(3, ib) = v_ref(3, ib) + zterm * rmb;

/*     increment the internal virial tensor components */

		    vxx = xr * xterm * vterm;
		    vyx = yr * xterm * vterm;
		    vzx = zr * xterm * vterm;
		    vyy = yr * yterm * vterm;
		    vzy = zr * yterm * vterm;
		    vzz = zr * zterm * vterm;
		    vir_ref(1, 1) = vir_ref(1, 1) - vxx;
		    vir_ref(2, 1) = vir_ref(2, 1) - vyx;
		    vir_ref(3, 1) = vir_ref(3, 1) - vzx;
		    vir_ref(1, 2) = vir_ref(1, 2) - vyx;
		    vir_ref(2, 2) = vir_ref(2, 2) - vyy;
		    vir_ref(3, 2) = vir_ref(3, 2) - vzy;
		    vir_ref(1, 3) = vir_ref(1, 3) - vzx;
		    vir_ref(2, 3) = vir_ref(2, 3) - vzy;
		    vir_ref(3, 3) = vir_ref(3, 3) - vzz;
		}
	    }
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    moved[i__ - 1] = update[i__ - 1];
	    update[i__ - 1] = FALSE_;
	}
    }

/*     write information on the number of iterations needed */

    if (niter == maxiter) {
	io___67.ciunit = iounit_1.iout;
	s_wsfe(&io___67);
	e_wsfe();
	prterr_();
	fatal_();
    } else if (inform_1.debug) {
	io___68.ciunit = iounit_1.iout;
	s_wsfe(&io___68);
	do_fio(&c__1, (char *)&niter, (ftnlen)sizeof(integer));
	e_wsfe();
    }

/*     apply any atom group velocity constraints via exact reset */

    i__1 = shake_1.nratx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = shake_1.iratx[i__ - 1];
	mode = shake_1.kratx[i__ - 1];
	xv = 0.;
	yv = 0.;
	zv = 0.;
	start = igrp_ref(1, ia);
	stop = igrp_ref(2, ia);
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    weigh = atmtyp_1.mass[k - 1] / group_1.grpmass[ia - 1];
	    if (mode > 2) {
		xv += v_ref(1, k) * weigh;
	    }
	    if (mode > 1) {
		yv += v_ref(2, k) * weigh;
	    }
	    zv += v_ref(3, k) * weigh;
	}
	i__2 = stop;
	for (j = start; j <= i__2; ++j) {
	    k = group_1.kgrp[j - 1];
	    v_ref(1, k) = v_ref(1, k) - xv;
	    v_ref(2, k) = v_ref(2, k) - yv;
	    v_ref(3, k) = v_ref(3, k) - zv;
	}
    }
    return 0;
} /* rattle2_ */

#undef igrp_ref
#undef irat_ref
#undef vir_ref
#undef v_ref


