/* active.f -- translated by f2c (version 20050501).
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
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine active  --  set the list of active atoms  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "active" sets the list of atoms that are used during */
/*     each potential energy function calculation */


/* Subroutine */ int active_(void)
{
    /* Format strings */
    static char fmt_40[] = "(/,\002 Active Site Spheres used to\002,\002 Sel"
	    "ect Active Atoms :\002,//,3x,\002Atom Center\002,11x,\002Coordin"
	    "ates\002,12x,\002Radius\002,6x,\002# Active Atoms\002)";
    static char fmt_50[] = "(2x,i8,6x,3f9.2,2x,f9.2,7x,i8)";
    static char fmt_70[] = "(/,\002 List of Active Atoms for Energy\002,\002"
	    " Calculations :\002,/)";
    static char fmt_80[] = "(3x,10i7)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfe(
	    cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, next;
    static doublereal dist2;
    static integer fixed[25000], mobile[25000], nfixed;
    static char record[120];
    static integer center;
    static doublereal radius;
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    static doublereal radius2;
    static integer nmobile, nsphere;
    static doublereal xcenter, ycenter, zcenter;
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___12 = { 1, string, 1, 0, 120, 1 };
    static icilist io___13 = { 1, string, 1, 0, 120, 1 };
    static icilist io___19 = { 1, string, 1, 0, 120, 1 };
    static icilist io___20 = { 1, string, 1, 0, 120, 1 };
    static cilist io___21 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_80, 0 };



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     set defaults for the numbers and lists of active atoms */

    usage_1.nuse = atoms_1.n;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	usage_1.use[i__ - 1] = TRUE_;
    }
    nmobile = 0;
    nfixed = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mobile[i__ - 1] = 0;
	fixed[i__ - 1] = 0;
    }
    nsphere = 0;

/*     get any keywords containing active atom parameters */

    i__1 = keys_1.nkey;
    for (j = 1; j <= i__1; ++j) {
	next = 1;
	s_copy(record, keyline_ref(0, j), (ftnlen)120, (ftnlen)120);
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));

/*     get any lists of atoms whose coordinates are active */

	if (s_cmp(keyword, "ACTIVE ", (ftnlen)7, (ftnlen)7) == 0) {
	    i__2 = s_rsli(&io___12);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__3 = atoms_1.n;
	    for (i__ = nmobile + 1; i__ <= i__3; ++i__) {
		i__2 = do_lio(&c__3, &c__1, (char *)&mobile[i__ - 1], (ftnlen)
			sizeof(integer));
		if (i__2 != 0) {
		    goto L10;
		}
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L10;
	    }
L10:
	    while(mobile[nmobile] != 0) {
		++nmobile;
/* Computing MAX */
/* Computing MIN */
		i__4 = atoms_1.n, i__5 = mobile[nmobile - 1];
		i__2 = -atoms_1.n, i__3 = min(i__4,i__5);
		mobile[nmobile - 1] = max(i__2,i__3);
	    }

/*     get any lists of atoms whose coordinates are inactive */

	} else if (s_cmp(keyword, "INACTIVE ", (ftnlen)9, (ftnlen)9) == 0) {
	    i__2 = s_rsli(&io___13);
	    if (i__2 != 0) {
		goto L20;
	    }
	    i__3 = atoms_1.n;
	    for (i__ = nfixed + 1; i__ <= i__3; ++i__) {
		i__2 = do_lio(&c__3, &c__1, (char *)&fixed[i__ - 1], (ftnlen)
			sizeof(integer));
		if (i__2 != 0) {
		    goto L20;
		}
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L20;
	    }
L20:
	    while(fixed[nfixed] != 0) {
		++nfixed;
/* Computing MAX */
/* Computing MIN */
		i__4 = atoms_1.n, i__5 = fixed[nfixed - 1];
		i__2 = -atoms_1.n, i__3 = min(i__4,i__5);
		fixed[nfixed - 1] = max(i__2,i__3);
	    }

/*     get the center and radius of the sphere of active atoms */

	} else if (s_cmp(keyword, "SPHERE ", (ftnlen)7, (ftnlen)7) == 0) {
	    center = 0;
	    xcenter = 0.;
	    ycenter = 0.;
	    zcenter = 0.;
	    radius = 0.;
	    i__2 = s_rsli(&io___19);
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&xcenter, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&ycenter, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&zcenter, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&radius, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L30;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L30;
	    }
L30:
	    if (radius == 0.) {
		i__2 = s_rsli(&io___20);
		if (i__2 != 0) {
		    goto L60;
		}
		i__2 = do_lio(&c__3, &c__1, (char *)&center, (ftnlen)sizeof(
			integer));
		if (i__2 != 0) {
		    goto L60;
		}
		i__2 = do_lio(&c__5, &c__1, (char *)&radius, (ftnlen)sizeof(
			doublereal));
		if (i__2 != 0) {
		    goto L60;
		}
		i__2 = e_rsli();
		if (i__2 != 0) {
		    goto L60;
		}
		xcenter = atoms_1.x[center - 1];
		ycenter = atoms_1.y[center - 1];
		zcenter = atoms_1.z__[center - 1];
	    }
	    ++nsphere;
	    if (nsphere == 1) {
		usage_1.nuse = 0;
		i__2 = atoms_1.n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    usage_1.use[i__ - 1] = FALSE_;
		}
		if (inform_1.verbose) {
		    io___21.ciunit = iounit_1.iout;
		    s_wsfe(&io___21);
		    e_wsfe();
		}
	    }
	    radius2 = radius * radius;
	    i__2 = atoms_1.n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (! usage_1.use[i__ - 1]) {
/* Computing 2nd power */
		    d__1 = atoms_1.x[i__ - 1] - xcenter;
/* Computing 2nd power */
		    d__2 = atoms_1.y[i__ - 1] - ycenter;
/* Computing 2nd power */
		    d__3 = atoms_1.z__[i__ - 1] - zcenter;
		    dist2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
		    if (dist2 <= radius2) {
			++usage_1.nuse;
			usage_1.use[i__ - 1] = TRUE_;
		    }
		}
	    }
	    if (inform_1.verbose) {
		io___24.ciunit = iounit_1.iout;
		s_wsfe(&io___24);
		do_fio(&c__1, (char *)&center, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&xcenter, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&ycenter, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&zcenter, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&radius, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&usage_1.nuse, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
L60:
	    ;
	}
    }

/*     set active atoms to those marked as not inactive */

    i__ = 1;
    while(fixed[i__ - 1] != 0) {
	if (fixed[i__ - 1] > 0) {
	    j = fixed[i__ - 1];
	    if (usage_1.use[j - 1]) {
		usage_1.use[fixed[i__ - 1] - 1] = FALSE_;
		--usage_1.nuse;
	    }
	    ++i__;
	} else {
	    i__3 = (i__2 = fixed[i__], abs(i__2));
	    for (j = (i__1 = fixed[i__ - 1], abs(i__1)); j <= i__3; ++j) {
		if (usage_1.use[j - 1]) {
		    usage_1.use[j - 1] = FALSE_;
		    --usage_1.nuse;
		}
	    }
	    i__ += 2;
	}
    }

/*     set active atoms to only those marked as active */

    i__ = 1;
    while(mobile[i__ - 1] != 0) {
	if (i__ == 1) {
	    usage_1.nuse = 0;
	    i__3 = atoms_1.n;
	    for (j = 1; j <= i__3; ++j) {
		usage_1.use[j - 1] = FALSE_;
	    }
	}
	if (mobile[i__ - 1] > 0) {
	    j = mobile[i__ - 1];
	    if (! usage_1.use[j - 1]) {
		usage_1.use[j - 1] = TRUE_;
		++usage_1.nuse;
	    }
	    ++i__;
	} else {
	    i__2 = (i__1 = mobile[i__], abs(i__1));
	    for (j = (i__3 = mobile[i__ - 1], abs(i__3)); j <= i__2; ++j) {
		if (! usage_1.use[j - 1]) {
		    usage_1.use[j - 1] = TRUE_;
		    ++usage_1.nuse;
		}
	    }
	    i__ += 2;
	}
    }

/*     use logical array to set the list of active atoms */

    j = 0;
    i__2 = atoms_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (usage_1.use[i__ - 1]) {
	    ++j;
	    usage_1.iuse[j - 1] = i__;
	}
    }

/*     output the final list of the active atoms */

    if (inform_1.debug && usage_1.nuse > 0 && usage_1.nuse < atoms_1.n) {
	io___25.ciunit = iounit_1.iout;
	s_wsfe(&io___25);
	e_wsfe();
	io___26.ciunit = iounit_1.iout;
	s_wsfe(&io___26);
	i__2 = usage_1.nuse;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_fio(&c__1, (char *)&usage_1.iuse[i__ - 1], (ftnlen)sizeof(
		    integer));
	}
	e_wsfe();
    }
    return 0;
} /* active_ */

#undef keyline_ref


