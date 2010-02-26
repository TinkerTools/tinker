/* inertia.f -- translated by f2c (version 20050501).
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
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine inertia  --  principal moments of inertia  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "inertia" computes the principal moments of inertia for the */
/*     system, and optionally translates the center of mass to the */
/*     origin and rotates the principal axes onto the global axes */

/*        mode = 1     print the moments and principal axes */
/*        mode = 2     move coordinates to standard orientation */
/*        mode = 3     perform both of the above operations */

/*     literature reference: */

/*     Herbert Goldstein, "Classical Mechanics, 2nd Edition", */
/*     Addison-Wesley, Reading, MA, 1980; see the Euler angle */
/*     xyz convention in Appendix B */


/* Subroutine */ int inertia_(integer *mode)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Moments of Inertia and Principal Axes "
	    ":\002,//,13x,\002Moments (amu Ang^2)\002,10x,\002X-, Y- and Z-Co"
	    "mponents of Axes\002)";
    static char fmt_30[] = "(3(/,11x,f16.3,10x,3f12.6))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static doublereal a[9]	/* was [3][3] */;
    static integer i__, j, k;
    static doublereal xx, xy, xz, yy, yz, zz, vec[9]	/* was [3][3] */, dot,
	     xcm, ycm, zcm;
    static logical move;
    static doublereal work1[3], work2[3], weigh, total;
    static logical print;
    static doublereal xterm, yterm, zterm;
    extern /* Subroutine */ int jacobi_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal moment[3], tensor[9]	/* was [3][3] */;

    /* Fortran I/O blocks */
    static cilist io___26 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_30, 0 };



#define a_ref(a_1,a_2) a[(a_2)*3 + a_1 - 4]
#define vec_ref(a_1,a_2) vec[(a_2)*3 + a_1 - 4]
#define tensor_ref(a_1,a_2) tensor[(a_2)*3 + a_1 - 4]



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




/*     decide upon the type of output desired */

    print = FALSE_;
    move = FALSE_;
    if (*mode == 1 || *mode == 3) {
	print = TRUE_;
    }
    if (*mode == 2 || *mode == 3) {
	move = TRUE_;
    }

/*     compute the position of the center of mass */

    total = 0.;
    xcm = 0.;
    ycm = 0.;
    zcm = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	weigh = atmtyp_1.mass[i__ - 1];
	total += weigh;
	xcm += atoms_1.x[i__ - 1] * weigh;
	ycm += atoms_1.y[i__ - 1] * weigh;
	zcm += atoms_1.z__[i__ - 1] * weigh;
    }
    xcm /= total;
    ycm /= total;
    zcm /= total;

/*     compute and then diagonalize the inertia tensor */

    xx = 0.;
    xy = 0.;
    xz = 0.;
    yy = 0.;
    yz = 0.;
    zz = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	weigh = atmtyp_1.mass[i__ - 1];
	xterm = atoms_1.x[i__ - 1] - xcm;
	yterm = atoms_1.y[i__ - 1] - ycm;
	zterm = atoms_1.z__[i__ - 1] - zcm;
	xx += xterm * xterm * weigh;
	xy += xterm * yterm * weigh;
	xz += xterm * zterm * weigh;
	yy += yterm * yterm * weigh;
	yz += yterm * zterm * weigh;
	zz += zterm * zterm * weigh;
    }
    tensor_ref(1, 1) = yy + zz;
    tensor_ref(2, 1) = -xy;
    tensor_ref(3, 1) = -xz;
    tensor_ref(1, 2) = -xy;
    tensor_ref(2, 2) = xx + zz;
    tensor_ref(3, 2) = -yz;
    tensor_ref(1, 3) = -xz;
    tensor_ref(2, 3) = -yz;
    tensor_ref(3, 3) = xx + yy;
    jacobi_(&c__3, &c__3, tensor, moment, vec, work1, work2);

/*     select the direction for each principal moment axis */

    for (i__ = 1; i__ <= 2; ++i__) {
	i__1 = atoms_1.n;
	for (j = 1; j <= i__1; ++j) {
	    xterm = vec_ref(1, i__) * (atoms_1.x[j - 1] - xcm);
	    yterm = vec_ref(2, i__) * (atoms_1.y[j - 1] - ycm);
	    zterm = vec_ref(3, i__) * (atoms_1.z__[j - 1] - zcm);
	    dot = xterm + yterm + zterm;
	    if (dot < 0.) {
		for (k = 1; k <= 3; ++k) {
		    vec_ref(k, i__) = -vec_ref(k, i__);
		}
	    }
	    if (dot != 0.) {
		goto L10;
	    }
	}
L10:
	;
    }

/*     moment axes must give a right-handed coordinate system */

    xterm = vec_ref(1, 1) * (vec_ref(2, 2) * vec_ref(3, 3) - vec_ref(2, 3) * 
	    vec_ref(3, 2));
    yterm = vec_ref(2, 1) * (vec_ref(1, 3) * vec_ref(3, 2) - vec_ref(1, 2) * 
	    vec_ref(3, 3));
    zterm = vec_ref(3, 1) * (vec_ref(1, 2) * vec_ref(2, 3) - vec_ref(1, 3) * 
	    vec_ref(2, 2));
    dot = xterm + yterm + zterm;
    if (dot < 0.) {
	for (j = 1; j <= 3; ++j) {
	    vec_ref(j, 3) = -vec_ref(j, 3);
	}
    }

/*     print the moments of inertia and the principal axes */

    if (print) {
	io___26.ciunit = iounit_1.iout;
	s_wsfe(&io___26);
	e_wsfe();
	io___27.ciunit = iounit_1.iout;
	s_wsfe(&io___27);
	for (i__ = 1; i__ <= 3; ++i__) {
	    do_fio(&c__1, (char *)&moment[i__ - 1], (ftnlen)sizeof(doublereal)
		    );
	    do_fio(&c__1, (char *)&vec_ref(1, i__), (ftnlen)sizeof(doublereal)
		    );
	    do_fio(&c__1, (char *)&vec_ref(2, i__), (ftnlen)sizeof(doublereal)
		    );
	    do_fio(&c__1, (char *)&vec_ref(3, i__), (ftnlen)sizeof(doublereal)
		    );
	}
	e_wsfe();
    }

/*     principal moment axes form rows of Euler rotation matrix */

    if (move) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		a_ref(i__, j) = vec_ref(j, i__);
	    }
	}

/*     translate to origin, then apply Euler rotation matrix */

	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xterm = atoms_1.x[i__ - 1] - xcm;
	    yterm = atoms_1.y[i__ - 1] - ycm;
	    zterm = atoms_1.z__[i__ - 1] - zcm;
	    atoms_1.x[i__ - 1] = a_ref(1, 1) * xterm + a_ref(1, 2) * yterm + 
		    a_ref(1, 3) * zterm;
	    atoms_1.y[i__ - 1] = a_ref(2, 1) * xterm + a_ref(2, 2) * yterm + 
		    a_ref(2, 3) * zterm;
	    atoms_1.z__[i__ - 1] = a_ref(3, 1) * xterm + a_ref(3, 2) * yterm 
		    + a_ref(3, 3) * zterm;
	}
    }
    return 0;
} /* inertia_ */

#undef tensor_ref
#undef vec_ref
#undef a_ref


