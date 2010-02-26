/* flatten.f -- translated by f2c (version 20050501).
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
    integer biotyp[10000];
    char forcefield[20];
} fields_;

#define fields_1 fields_

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
    integer nkey;
    char keyline[1200000];
} keys_;

#define keys_1 keys_

struct {
    doublereal m2[25000], deform, difft, diffv, diffc;
    logical use_smooth__, use_dem__, use_gda__, use_tophat__, use_stophat__;
} warp_;

#define warp_1 warp_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine flatten  --  set potential smoothing parameters  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "flatten" sets the type of smoothing method and the extent of */
/*     surface deformation for use with potential energy smoothing */


/* Subroutine */ int flatten_(void)
{
    /* Format strings */
    static char fmt_30[] = "(/,\002 Enter the Initial Mean Squared Gaussia"
	    "n\002,\002 Width [200.0] :  \002,$)";
    static char fmt_40[] = "(/,\002 Enter Length Scale for Potential Surfac"
	    "e\002,\002 Averaging [0.0] :  \002,$)";
    static char fmt_50[] = "(/,\002 Enter the Potential Surface Smoothing"
	    "\002,\002 Parameter [0.0] :  \002,$)";
    static char fmt_60[] = "(a120)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *, char 
	    *, ftnlen), e_rsfe(void);

    /* Local variables */
    static integer i__, next;
    static logical exist;
    static char stype[7];
    static logical query;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), getword_(
	    char *, char *, integer *, ftnlen, ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___8 = { 1, string, 1, 0, 120, 1 };
    static icilist io___9 = { 1, string, 1, 0, 120, 1 };
    static icilist io___10 = { 1, string, 1, 0, 120, 1 };
    static icilist io___11 = { 1, string, 1, 0, 120, 1 };
    static icilist io___13 = { 1, string, 1, 0, 120, 1 };
    static cilist io___14 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_60, 0 };
    static icilist io___18 = { 1, record, 1, 0, 120, 1 };



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
/*     ##  keys.i  --  contents of current keyword parameter file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     nkey      number of nonblank lines in the keyword file */
/*     keyline   contents of each individual keyword file line */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  warp.i  --  parameters for potential surface smoothing  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     m2           second moment of the GDA gaussian for each atom */
/*     deform       value of the smoothing deformation parameter */
/*     difft        diffusion coefficient for torsional potential */
/*     diffv        diffusion coefficient for van der Waals potential */
/*     diffc        diffusion coefficient for charge-charge potential */
/*     use_smooth   flag to use a potential energy smoothing method */
/*     use_dem      flag to use diffusion equation method potential */
/*     use_gda      flag to use gaussian density annealing potential */
/*     use_tophat   flag to use analytical tophat smoothed potential */
/*     use_stophat  flag to use shifted tophat smoothed potential */




/*     set defaults for deformation and diffusion coefficients */

    query = TRUE_;
    warp_1.deform = 0.;
    warp_1.difft = .0225;
    warp_1.diffv = 1.;
    warp_1.diffc = 1.;

/*     get any keywords related to potential energy smoothing */

    i__1 = keys_1.nkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next = 1;
	s_copy(record, keyline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "SMOOTHING ", (ftnlen)10, (ftnlen)10) == 0) {
	    warp_1.use_smooth__ = TRUE_;
	    warp_1.use_dem__ = FALSE_;
	    warp_1.use_gda__ = FALSE_;
	    warp_1.use_tophat__ = FALSE_;
	    warp_1.use_stophat__ = FALSE_;
	    getword_(record, stype, &next, (ftnlen)120, (ftnlen)7);
	    upcase_(stype, (ftnlen)7);
	    if (s_cmp(stype, "DEM", (ftnlen)7, (ftnlen)3) == 0) {
		warp_1.use_dem__ = TRUE_;
	    }
	    if (s_cmp(stype, "GDA", (ftnlen)7, (ftnlen)3) == 0) {
		warp_1.use_gda__ = TRUE_;
	    }
	    if (s_cmp(stype, "TOPHAT", (ftnlen)7, (ftnlen)6) == 0) {
		warp_1.use_tophat__ = TRUE_;
	    }
	    if (s_cmp(stype, "STOPHAT", (ftnlen)7, (ftnlen)7) == 0) {
		warp_1.use_stophat__ = TRUE_;
	    }
	} else if (s_cmp(keyword, "DEFORM ", (ftnlen)7, (ftnlen)7) == 0) {
	    i__2 = s_rsli(&io___8);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&warp_1.deform, (ftnlen)
		    sizeof(doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	    query = FALSE_;
	} else if (s_cmp(keyword, "DIFFUSE-TORSION ", (ftnlen)16, (ftnlen)16) 
		== 0) {
	    i__2 = s_rsli(&io___9);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&warp_1.difft, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "DIFFUSE-VDW ", (ftnlen)12, (ftnlen)12) == 
		0) {
	    i__2 = s_rsli(&io___10);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&warp_1.diffv, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	} else if (s_cmp(keyword, "DIFFUSE-CHARGE ", (ftnlen)15, (ftnlen)15) 
		== 0) {
	    i__2 = s_rsli(&io___11);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&warp_1.diffc, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
	}
L10:
	;
    }

/*     try to get the deformation value from the command line */

    if (warp_1.use_smooth__) {
	if (query) {
	    nextarg_(string, &exist, (ftnlen)120);
	    if (exist) {
		i__1 = s_rsli(&io___13);
		if (i__1 != 0) {
		    goto L20;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&warp_1.deform, (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L20;
		}
		i__1 = e_rsli();
		if (i__1 != 0) {
		    goto L20;
		}
		query = FALSE_;
	    }
L20:
	    ;
	}

/*     ask for the potential surface deformation to be used */

	if (query) {
	    if (warp_1.use_gda__) {
		warp_1.deform = 200.;
		io___14.ciunit = iounit_1.iout;
		s_wsfe(&io___14);
		e_wsfe();
	    } else if (warp_1.use_tophat__ || warp_1.use_stophat__) {
		warp_1.deform = 0.;
		io___15.ciunit = iounit_1.iout;
		s_wsfe(&io___15);
		e_wsfe();
	    } else {
		warp_1.deform = 0.;
		io___16.ciunit = iounit_1.iout;
		s_wsfe(&io___16);
		e_wsfe();
	    }
	    io___17.ciunit = iounit_1.input;
	    s_rsfe(&io___17);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    i__1 = s_rsli(&io___18);
	    if (i__1 != 0) {
		goto L70;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&warp_1.deform, (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L70;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L70;
	    }
L70:
	    ;
	}
    }

/*     set second moment of Gaussian on each atom for GDA methods */

    if (warp_1.use_gda__) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    warp_1.m2[i__ - 1] = warp_1.deform;
	}
    }
    return 0;
} /* flatten_ */

#undef keyline_ref


