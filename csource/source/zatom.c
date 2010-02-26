/* zatom.f -- translated by f2c (version 20050501).
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
    integer biotyp[10000];
    char forcefield[20];
} fields_;

#define fields_1 fields_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal weight[5000];
    integer atmcls[5000], atmnum[5000], ligand[5000];
    char symbol[15000], describe[120000];
} katoms_;

#define katoms_1 katoms_

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

/* Table of constant values */

static integer c__25000 = 25000;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine zatom  --  adds a single atom to Z-matrix  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "zatom" adds an atom to the end of the current Z-matrix */
/*     and then increments the atom counter; atom type, defining */
/*     atoms and internal coordinates are passed as arguments */


/* Subroutine */ int zatom_(integer *bionum, doublereal *bond, doublereal *
	angle, doublereal *dihed, integer *iz1, integer *iz2, integer *iz3, 
	integer *iz4)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 ZATOM  --  The Maximum of\002,i8,\002 At"
	    "oms\002,\002 has been Exceeded\002)";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int fatal_(void);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };



#define iz_ref(a_1,a_2) zcoord_1.iz[(a_2)*4 + a_1 - 5]
#define iadd_ref(a_1,a_2) zclose_1.iadd[(a_2)*2 + a_1 - 3]
#define idel_ref(a_1,a_2) zclose_1.idel[(a_2)*2 + a_1 - 3]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define symbol_ref(a_0,a_1) &katoms_1.symbol[(a_1)*3 + a_0 - 3]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  katoms.i  --  forcefield parameters for the atom types  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     weight     average atomic mass of each atom type */
/*     atmcls     atom class number for each of the atom types */
/*     atmnum     atomic number for each of the atom types */
/*     ligand     number of atoms to be attached to each atom type */
/*     symbol     modified atomic symbol for each atom type */
/*     describe   string identifying each of the atom types */




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




/*     fill various arrays with information for this atom */

    if (*bionum > 0) {
	atoms_1.type__[atoms_1.n - 1] = fields_1.biotyp[*bionum - 1];
	if (atoms_1.type__[atoms_1.n - 1] != 0) {
	    s_copy(name___ref(0, atoms_1.n), symbol_ref(0, atoms_1.type__[
		    atoms_1.n - 1]), (ftnlen)3, (ftnlen)3);
	} else {
	    s_copy(name___ref(0, atoms_1.n), "   ", (ftnlen)3, (ftnlen)3);
	}
	zcoord_1.zbond[atoms_1.n - 1] = *bond;
	zcoord_1.zang[atoms_1.n - 1] = *angle;
	zcoord_1.ztors[atoms_1.n - 1] = *dihed;
	if (zcoord_1.ztors[atoms_1.n - 1] < -180.) {
	    zcoord_1.ztors[atoms_1.n - 1] += 360.;
	} else if (zcoord_1.ztors[atoms_1.n - 1] > 180.) {
	    zcoord_1.ztors[atoms_1.n - 1] += -360.;
	}
	iz_ref(1, atoms_1.n) = *iz1;
	iz_ref(2, atoms_1.n) = *iz2;
	iz_ref(3, atoms_1.n) = *iz3;
	iz_ref(4, atoms_1.n) = *iz4;

/*     increment atom counter and check for too many atoms */

	++atoms_1.n;
	if (atoms_1.n > 25000) {
	    io___1.ciunit = iounit_1.iout;
	    s_wsfe(&io___1);
	    do_fio(&c__1, (char *)&c__25000, (ftnlen)sizeof(integer));
	    e_wsfe();
	    fatal_();
	}

/*     add an extra bond to make a ring closure */

    } else if (*bionum == -1) {
	++zclose_1.nadd;
	iadd_ref(1, zclose_1.nadd) = *iz1;
	iadd_ref(2, zclose_1.nadd) = *iz2;

/*     delete an extra bond to make separate molecules */

    } else if (*bionum == -2) {
	++zclose_1.ndel;
	idel_ref(1, zclose_1.ndel) = *iz1;
	idel_ref(2, zclose_1.ndel) = *iz2;
    }
    return 0;
} /* zatom_ */

#undef symbol_ref
#undef name___ref
#undef idel_ref
#undef iadd_ref
#undef iz_ref


