/* impose.f -- translated by f2c (version 20050501).
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
    doublereal wfit[25000];
    integer nfit, ifit[50000]	/* was [2][25000] */;
} align_;

#define align_1 align_

struct {
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

/* Table of constant values */

static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine impose  --  superimpose two coordinate sets  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "impose" performs the least squares best superposition */
/*     of two atomic coordinate sets via a quaternion method; */
/*     upon return, the first coordinate set is unchanged while */
/*     the second set is translated and rotated to give best fit; */
/*     the final root mean square fit is returned in "rmsvalue" */


/* Subroutine */ int impose_(integer *n1, doublereal *x1, doublereal *y1, 
	doublereal *z1, integer *n2, doublereal *x2, doublereal *y2, 
	doublereal *z2, doublereal *rmsvalue)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 IMPOSE  --  Input Coordinates\002,12x,f1"
	    "2.6)";
    static char fmt_30[] = "(\002 IMPOSE  --  After Translation\002,12x,f12."
	    "6)";
    static char fmt_40[] = "(\002 IMPOSE  --  After Rotation\002,15x,f12.6)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__;
    static doublereal xmid, ymid, zmid;
    extern /* Subroutine */ int center_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *);
    extern doublereal rmsfit_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int quatfit_(integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, doublereal 
	    *);

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_40, 0 };



#define ifit_ref(a_1,a_2) align_1.ifit[(a_2)*2 + a_1 - 3]



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  align.i  --  information for superposition of structures  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     wfit    weights assigned to atom pairs during superposition */
/*     nfit    number of atoms to use in superimposing two structures */
/*     ifit    atom numbers of pairs of atoms to be superimposed */




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




/*     superimpose the full structures if not specified */

    /* Parameter adjustments */
    --z2;
    --y2;
    --x2;
    --z1;
    --y1;
    --x1;

    /* Function Body */
    if (align_1.nfit == 0) {
	align_1.nfit = min(*n1,*n2);
	i__1 = align_1.nfit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ifit_ref(1, i__) = i__;
	    ifit_ref(2, i__) = i__;
	    align_1.wfit[i__ - 1] = 1.;
	}
    }

/*     if the weights are all zero, set them to unity */

    i__1 = align_1.nfit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (align_1.wfit[i__ - 1] != 0.) {
	    goto L10;
	}
    }
    i__1 = align_1.nfit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	align_1.wfit[i__ - 1] = 1.;
    }
L10:

/*     find the rms fit of input coordinates */

    if (inform_1.verbose) {
	*rmsvalue = rmsfit_(&x1[1], &y1[1], &z1[1], &x2[1], &y2[1], &z2[1]);
	io___2.ciunit = iounit_1.iout;
	s_wsfe(&io___2);
	do_fio(&c__1, (char *)&(*rmsvalue), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     superimpose the centroids of active atom pairs */

    center_(n1, &x1[1], &y1[1], &z1[1], n2, &x2[1], &y2[1], &z2[1], &xmid, &
	    ymid, &zmid);
    if (inform_1.verbose) {
	*rmsvalue = rmsfit_(&x1[1], &y1[1], &z1[1], &x2[1], &y2[1], &z2[1]);
	io___6.ciunit = iounit_1.iout;
	s_wsfe(&io___6);
	do_fio(&c__1, (char *)&(*rmsvalue), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     use a quaternion method to achieve the superposition */

    quatfit_(n1, &x1[1], &y1[1], &z1[1], n2, &x2[1], &y2[1], &z2[1]);
    *rmsvalue = rmsfit_(&x1[1], &y1[1], &z1[1], &x2[1], &y2[1], &z2[1]);
    if (inform_1.verbose) {
	io___7.ciunit = iounit_1.iout;
	s_wsfe(&io___7);
	do_fio(&c__1, (char *)&(*rmsvalue), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     translate both coordinate sets so as to return */
/*     the first set to its original position */

    i__1 = *n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x1[i__] += xmid;
	y1[i__] += ymid;
	z1[i__] += zmid;
    }
    i__1 = *n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x2[i__] += xmid;
	y2[i__] += ymid;
	z2[i__] += zmid;
    }
    return 0;
} /* impose_ */

#undef ifit_ref


