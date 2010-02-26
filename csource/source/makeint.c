/* makeint.f -- translated by f2c (version 20050501).
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
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nadd, iadd[50000]	/* was [2][25000] */, ndel, idel[50000]	/* 
	    was [2][25000] */;
} zclose_;

#define zclose_1 zclose_

struct {
    doublereal zbond[25000], zang[25000], ztors[25000];
    integer iz[100000]	/* was [4][25000] */;
} zcoord_;

#define zcoord_1 zcoord_

struct {
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine makeint  --  convert Cartesian to internal  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "makeint" converts Cartesian to internal coordinates where */
/*     selection of internal coordinates is controlled by "mode" */

/*        mode = 0     automatic internal coordinates */
/*        mode = 1     manual selection of coordinates */
/*        mode = 2     use existing structure as a template */
/*        mode = 3     use dihedral angles in all cases */


/* Subroutine */ int makeint_(integer *mode)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Atom Number to be Defined [\002,i5,\002]"
	    " :  \002,$)";
    static char fmt_30[] = "(i10)";
    static char fmt_40[] = "(/,\002 Already Defined that Atom; Choose Anot"
	    "her\002)";
    static char fmt_50[] = "(/,\002 MAKEINT  --  Connectivity Error\002,\002"
	    " in defining Atom\002,i6)";
    static char fmt_60[] = "(/,\002 MAKEINT  --  Connectivity Error\002,\002"
	    " in defining Atom\002,i6)";
    static char fmt_70[] = "(/,\002 Specify with Dihedral Angle or Second"
	    "\002,\002 Bond Angle (\002,a8,\002) :  \002,$)";
    static char fmt_80[] = "(a120)";
    static char fmt_90[] = "(/,\002 MAKEINT  --  Connectivity Error\002,\002"
	    " in defining Atom\002,i6)";
    static char fmt_100[] = "(/,\002 MAKEINT  --  Connectivity Error\002,"
	    "\002 in defining Atom\002,i6)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_rsfe(cilist *), e_rsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern integer adjacent_(integer *, integer *, integer *, logical *, 
	    integer *, integer *);
    extern doublereal geometry_(integer *, integer *, integer *, integer *);
    static integer i__, j, i1, i2, i3, i4, i5, iz0[25001], iz1[25000];
    static doublereal sign;
    static logical more;
    static integer next;
    extern /* Subroutine */ int fatal_(void);
    static integer trial;
    static char record[120], phrase[8];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char answer[1], default__[1];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___12 = { 1, 0, 0, fmt_30, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_100, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define iz_ref(a_1,a_2) zcoord_1.iz[(a_2)*4 + a_1 - 5]
#define iadd_ref(a_1,a_2) zclose_1.iadd[(a_2)*2 + a_1 - 3]



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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  zclose.i  --  ring openings and closures for Z-matrix  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     nadd   number of added bonds between Z-matrix atoms */
/*     iadd   numbers of the atom pairs defining added bonds */
/*     ndel   number of bonds between Z-matrix bonds to delete */
/*     idel   numbers of the atom pairs defining deleted bonds */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     zbond   bond length used to define each Z-matrix atom */
/*     zang    bond angle used to define each Z-matrix atom */
/*     ztors   angle or torsion used to define Z-matrix atom */
/*     iz      defining atom numbers for each Z-matrix atom */




/*     zero out local values used for the defining atoms */

    i1 = 0;
    i2 = 0;
    i3 = 0;
    i4 = 0;
    i5 = 0;
    iz0[0] = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iz0[i__] = 0;
	iz1[i__ - 1] = 0;
    }

/*     zero out the coordinates, defining atoms and closures */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zcoord_1.zbond[i__ - 1] = 0.;
	zcoord_1.zang[i__ - 1] = 0.;
	zcoord_1.ztors[i__ - 1] = 0.;
    }
    if (*mode != 2) {
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 4; ++j) {
		iz_ref(j, i__) = 0;
	    }
	}
	zclose_1.nadd = 0;
	zclose_1.ndel = 0;
    }

/*     first, decide which of the atoms to define next */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*mode == 1) {
	    trial = i1 + 1;
L10:
	    io___11.ciunit = iounit_1.iout;
	    s_wsfe(&io___11);
	    do_fio(&c__1, (char *)&trial, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___12.ciunit = iounit_1.input;
	    i__2 = s_rsfe(&io___12);
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L10;
	    }
	    i__2 = e_rsfe();
	    if (i__2 != 0) {
		goto L10;
	    }
	    if (i1 == 0) {
		i1 = trial;
	    }
	    if (iz0[i1] != 0) {
		io___13.ciunit = iounit_1.iout;
		s_wsfe(&io___13);
		e_wsfe();
		if (i1 == trial) {
		    ++trial;
		}
		goto L10;
	    }
	} else {
	    i1 = i__;
	}

/*     define the bond length for the current atom */

	if (i__ >= 2) {
	    if (*mode == 2) {
		i2 = iz_ref(1, i1);
	    } else {
		i2 = adjacent_(&i1, &c__0, mode, &more, iz0, iz1);
		if (i2 == 0) {
		    io___15.ciunit = iounit_1.iout;
		    s_wsfe(&io___15);
		    do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(integer));
		    e_wsfe();
		    fatal_();
		}
	    }
	    zcoord_1.zbond[i1 - 1] = geometry_(&i1, &i2, &c__0, &c__0);
	}

/*     define the bond angle for the current atom */

	if (i__ >= 3) {
	    if (*mode == 2) {
		i3 = iz_ref(2, i1);
	    } else {
		i3 = adjacent_(&i2, &i1, mode, &more, iz0, iz1);
		if (i3 == 0) {
		    io___16.ciunit = iounit_1.iout;
		    s_wsfe(&io___16);
		    do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(integer));
		    e_wsfe();
		    fatal_();
		}
	    }
	    zcoord_1.zang[i1 - 1] = geometry_(&i1, &i2, &i3, &c__0);
	}

/*     decide whether to use a dihedral or second bond angle; */
/*     then find the value of the angle */

	if (i__ >= 4) {
	    if (*mode == 3) {
		*(unsigned char *)answer = 'D';
	    } else if (*mode == 2) {
		if (iz_ref(4, i1) == 0) {
		    *(unsigned char *)answer = 'D';
		} else {
		    *(unsigned char *)answer = 'B';
		}
	    } else if (*mode == 1) {
		if (more) {
		    s_copy(phrase, "D or [B]", (ftnlen)8, (ftnlen)8);
		    *(unsigned char *)default__ = 'B';
		} else {
		    s_copy(phrase, "[D] or B", (ftnlen)8, (ftnlen)8);
		    *(unsigned char *)default__ = 'D';
		}
		io___20.ciunit = iounit_1.iout;
		s_wsfe(&io___20);
		do_fio(&c__1, phrase, (ftnlen)8);
		e_wsfe();
		io___21.ciunit = iounit_1.input;
		s_rsfe(&io___21);
		do_fio(&c__1, record, (ftnlen)120);
		e_rsfe();
		next = 1;
		gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
		upcase_(answer, (ftnlen)1);
		if (*(unsigned char *)answer != 'B' && *(unsigned char *)
			answer != 'D') {
		    *(unsigned char *)answer = *(unsigned char *)default__;
		}
	    } else if (*mode == 0) {
		if (more) {
		    *(unsigned char *)answer = 'B';
		} else {
		    *(unsigned char *)answer = 'D';
		}
	    }
	    if (*(unsigned char *)answer == 'B') {
		if (*mode == 2) {
		    i4 = iz_ref(3, i1);
		} else {
		    i4 = adjacent_(&i2, &i3, mode, &more, iz0, iz1);
		    if (i4 == 0) {
			io___24.ciunit = iounit_1.iout;
			s_wsfe(&io___24);
			do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(integer));
			e_wsfe();
			fatal_();
		    }
		}
		zcoord_1.ztors[i1 - 1] = geometry_(&i1, &i2, &i4, &c__0);
		i5 = 1;
		sign = geometry_(&i1, &i2, &i3, &i4);
		if (sign > 0.) {
		    i5 = -1;
		}
	    } else if (*(unsigned char *)answer == 'D') {
		if (*mode == 2) {
		    i4 = iz_ref(3, i1);
		} else {
		    i4 = adjacent_(&i3, &i2, mode, &more, iz0, iz1);
		    if (i4 == 0) {
			io___26.ciunit = iounit_1.iout;
			s_wsfe(&io___26);
			do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(integer));
			e_wsfe();
			fatal_();
		    }
		}
		i5 = 0;
		zcoord_1.ztors[i1 - 1] = geometry_(&i1, &i2, &i3, &i4);
	    }
	}

/*     transfer defining atoms to permanent array; */
/*     mark the current atom as finished */

	iz_ref(1, i1) = iz0[i2];
	iz_ref(2, i1) = iz0[i3];
	iz_ref(3, i1) = iz0[i4];
	iz_ref(4, i1) = i5;
	iz0[i1] = i__;
	iz1[i1 - 1] = i2;
    }

/*     add any bonds needed to make ring closures */

    zclose_1.nadd = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = couple_1.n12[i__ - 1];
	for (j = 1; j <= i__2; ++j) {
	    if (iz0[i__] < iz0[i12_ref(j, i__)] && iz1[i12_ref(j, i__) - 1] !=
		     i__) {
		++zclose_1.nadd;
		iadd_ref(1, zclose_1.nadd) = iz0[i__];
		iadd_ref(2, zclose_1.nadd) = iz0[i12_ref(j, i__)];
	    }
	}
    }
    return 0;
} /* makeint_ */

#undef iadd_ref
#undef iz_ref
#undef i12_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function adjacent  --  atom adjacent to specified atom  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "adjacent" finds an atom connected to atom "i1" other than */
/*     atom "i2"; if no such atom exists, then the closest atom */
/*     in space is returned */

/*     variables and parameters : */

/*     mode   whether "makeint" is in manual mode, automatic, etc. */
/*     more   returned true if there is more than one previously */
/*              defined atom other than "i2" which is directly */
/*              connected (adjacent) to atom "i1" */
/*     iz0    line number of the Z-matrix on which an atom is */
/*              defined, 0 if not yet defined */
/*     iz1    line number of the Z-matrix on which the atom used */
/*              defining the bond length to a given atom is defined */


integer adjacent_(integer *i1, integer *i2, integer *mode, logical *more, 
	integer *iz0, integer *iz1)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 ADJACENT  --  Atom\002,i6,\002 not Attac"
	    "hed\002,\002 to any Prior Atom\002)";
    static char fmt_20[] = "(/,\002 ADJACENT  --  Atom\002,i6,\002 not Attac"
	    "hed\002,\002 to any Prior Atom\002)";
    static char fmt_40[] = "(/,\002 ADJACENT  --  Atom\002,i6,\002 is the on"
	    "ly\002,\002 Connected Atom\002)";
    static char fmt_60[] = "(\002 Choose a Connected Atom (\002,2i6,\002) "
	    ":  \002,$)";
    static char fmt_70[] = "(\002 Choose a Connected Atom (\002,3i6,\002) "
	    ":  \002,$)";
    static char fmt_80[] = "(\002 Choose a Connected Atom (\002,4i6,\002) "
	    ":  \002,$)";
    static char fmt_90[] = "(\002 Choose a Connected Atom (\002,5i6,\002) "
	    ":  \002,$)";
    static char fmt_100[] = "(\002 Choose a Connected Atom (\002,6i6,\002) :"
	    "  \002,$)";
    static char fmt_110[] = "(\002 Choose a Connected Atom (\002,7i6,\002) :"
	    "  \002,$)";
    static char fmt_120[] = "(\002 Choose a Connected Atom (\002,8i6,\002) :"
	    "  \002,$)";
    static char fmt_130[] = "(i10)";

    /* System generated locals */
    integer ret_val, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal);
    integer s_rsfe(cilist *), e_rsfe(void);

    /* Local variables */
    static integer i__, j, k, ic[8], nc;
    static doublereal dist, short__;

    /* Fortran I/O blocks */
    static cilist io___31 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___44 = { 1, 0, 0, fmt_130, 0 };



#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define idel_ref(a_1,a_2) zclose_1.idel[(a_2)*2 + a_1 - 3]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  zclose.i  --  ring openings and closures for Z-matrix  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     nadd   number of added bonds between Z-matrix atoms */
/*     iadd   numbers of the atom pairs defining added bonds */
/*     ndel   number of bonds between Z-matrix bonds to delete */
/*     idel   numbers of the atom pairs defining deleted bonds */




/*     get a list of eligible atoms bonded to the atom of interest */

    /* Parameter adjustments */
    --iz1;

    /* Function Body */
    nc = 0;
    *more = FALSE_;
    i__1 = couple_1.n12[*i1 - 1];
    for (j = 1; j <= i__1; ++j) {
	i__ = i12_ref(j, *i1);
	if (iz0[i__] != 0 && i__ != *i2) {
	    if (*i2 == 0) {
		++nc;
		ic[nc - 1] = i__;
	    } else {
		if (iz1[i__] == *i1 || iz1[*i1] == i__) {
		    ++nc;
		    ic[nc - 1] = i__;
		}
	    }
	}
    }
    if (nc > 1) {
	*more = TRUE_;
    }

/*     if no bonded atom is eligible, use the nearest neighbor */

    if (nc == 0) {
	ret_val = 0;
	if (*mode == 1) {
	    io___31.ciunit = iounit_1.iout;
	    s_wsfe(&io___31);
	    do_fio(&c__1, (char *)&(*i1), (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    short__ = 1e6;
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (iz0[i__] != 0 && i__ != *i1 && i__ != *i2) {
/* Computing 2nd power */
		    d__1 = atoms_1.x[i__ - 1] - atoms_1.x[*i1 - 1];
/* Computing 2nd power */
		    d__2 = atoms_1.y[i__ - 1] - atoms_1.y[*i1 - 1];
/* Computing 2nd power */
		    d__3 = atoms_1.z__[i__ - 1] - atoms_1.z__[*i1 - 1];
		    dist = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
		    if (dist < short__) {
			short__ = dist;
			ret_val = i__;
		    }
		}
	    }
	    if (*i2 == 0) {
		++zclose_1.ndel;
		idel_ref(1, zclose_1.ndel) = ret_val;
		idel_ref(2, zclose_1.ndel) = *i1;
		if (inform_1.debug) {
		    io___34.ciunit = iounit_1.iout;
		    s_wsfe(&io___34);
		    do_fio(&c__1, (char *)&(*i1), (ftnlen)sizeof(integer));
		    e_wsfe();
		}
	    }
	}

/*     for automatic mode, always use the first eligible bonded atom */

    } else if (*mode == 0) {
	ret_val = ic[0];

/*     for torsion mode, use an adjacent atom bonded to undefined atoms */

    } else if (*mode == 3) {
	ret_val = ic[0];
	i__1 = nc;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = couple_1.n12[ic[k - 1] - 1];
	    for (j = 1; j <= i__2; ++j) {
		i__ = i12_ref(j, ic[k - 1]);
		if (iz0[i__] != 0 && i__ != *i1) {
		    ret_val = ic[k - 1];
		    goto L30;
		}
	    }
	}
L30:

/*     if only one directly bonded atom is eligible, then use it */

	;
    } else if (nc == 1) {
	ret_val = ic[0];
	if (*mode == 1 || inform_1.debug) {
	    io___36.ciunit = iounit_1.iout;
	    s_wsfe(&io___36);
	    do_fio(&c__1, (char *)&ic[0], (ftnlen)sizeof(integer));
	    e_wsfe();
	}

/*     ask the user which eligible bonded atom to use as adjacent */

    } else {
L50:
	if (nc == 2) {
	    io___37.ciunit = iounit_1.iout;
	    s_wsfe(&io___37);
	    i__1 = nc;
	    for (j = 1; j <= i__1; ++j) {
		do_fio(&c__1, (char *)&ic[j - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	} else if (nc == 3) {
	    io___38.ciunit = iounit_1.iout;
	    s_wsfe(&io___38);
	    i__1 = nc;
	    for (j = 1; j <= i__1; ++j) {
		do_fio(&c__1, (char *)&ic[j - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	} else if (nc == 4) {
	    io___39.ciunit = iounit_1.iout;
	    s_wsfe(&io___39);
	    i__1 = nc;
	    for (j = 1; j <= i__1; ++j) {
		do_fio(&c__1, (char *)&ic[j - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	} else if (nc == 5) {
	    io___40.ciunit = iounit_1.iout;
	    s_wsfe(&io___40);
	    i__1 = nc;
	    for (j = 1; j <= i__1; ++j) {
		do_fio(&c__1, (char *)&ic[j - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	} else if (nc == 6) {
	    io___41.ciunit = iounit_1.iout;
	    s_wsfe(&io___41);
	    i__1 = nc;
	    for (j = 1; j <= i__1; ++j) {
		do_fio(&c__1, (char *)&ic[j - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	} else if (nc == 7) {
	    io___42.ciunit = iounit_1.iout;
	    s_wsfe(&io___42);
	    i__1 = nc;
	    for (j = 1; j <= i__1; ++j) {
		do_fio(&c__1, (char *)&ic[j - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	} else if (nc == 8) {
	    io___43.ciunit = iounit_1.iout;
	    s_wsfe(&io___43);
	    i__1 = nc;
	    for (j = 1; j <= i__1; ++j) {
		do_fio(&c__1, (char *)&ic[j - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	}
	io___44.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___44);
	if (i__1 != 0) {
	    goto L50;
	}
	i__1 = do_fio(&c__1, (char *)&ret_val, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L50;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L50;
	}
	if (ret_val == 0) {
	    ret_val = ic[0];
	} else {
	    i__1 = nc;
	    for (j = 1; j <= i__1; ++j) {
		if (ic[j - 1] == ret_val) {
		    goto L140;
		}
	    }
	    goto L50;
L140:
	    ;
	}
    }
    return ret_val;
} /* adjacent_ */

#undef idel_ref
#undef i12_ref


