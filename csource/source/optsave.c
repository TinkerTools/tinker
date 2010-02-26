/* optsave.f -- translated by f2c (version 20050501).
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
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal dihed[1000];
    integer nomega, iomega[2000]	/* was [2][1000] */, zline[1000];
} omega_;

#define omega_1 omega_

struct {
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

struct {
    doublereal scale[75000];
    logical set_scale__;
} scales_;

#define scales_1 scales_

struct {
    integer runtyp, cstep;
    doublereal cdt, cenergy, cdx[25000], cdy[25000], cdz[25000];
    logical use_socket__, skt_init__, skt_close__;
} socket_;

#define socket_1 socket_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

struct {
    doublereal zbond[25000], zang[25000], ztors[25000];
    integer iz[100000]	/* was [4][25000] */;
} zcoord_;

#define zcoord_1 zcoord_

/* Table of constant values */

static integer c__3 = 3;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine optsave  --  save optimization info and results  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "optsave" is used by the optimizers to write imtermediate */
/*     coordinates and other relevant information; also checks for */
/*     user requested termination of an optimization */


/* Subroutine */ int optsave_(integer *ncycle, doublereal *f, doublereal *xx)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 OPTSAVE  --  Optimization Calculation En"
	    "ding\002,\002 due to User Request\002)";

    /* System generated locals */
    address a__1[3], a__2[2];
    integer i__1, i__2[3], i__3[2];
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_inqu(inlist *), f_open(olist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_rew(alist *), f_clos(cllist *), s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    extern integer freeunit_(void);
    static integer i__;
    static char ext[7];
    static integer iend, nvar, iopt, lext;
    extern /* Subroutine */ int fatal_(void);
    static logical exist;
    extern /* Subroutine */ int suffix_(char *, char *, ftnlen, ftnlen), 
	    prtint_(integer *), sktopt_(integer *, doublereal *), prtxyz_(
	    integer *);
    static char endfile[120];
    extern /* Subroutine */ int openend_(integer *, char *, ftnlen);
    static char optfile[120];
    extern /* Subroutine */ int numeral_(integer *, char *, integer *, ftnlen)
	    , version_(char *, char *, ftnlen, ftnlen), makexyz_(void);

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 0, 0, fmt_10, 0 };




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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  output.i  --  control of coordinate output file format  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     archive    logical flag to save structures in an archive */
/*     noversion  logical flag governing use of filename versions */
/*     overwrite  logical flag to overwrite intermediate files inplace */
/*     cyclesave  logical flag to mark use of numbered cycle files */
/*     coordtype  selects Cartesian, internal, rigid body or none */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  scales.i  --  parameter scale factors for optimization  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     scale      multiplicative factor for each optimization parameter */
/*     set_scale  logical flag to show if scale factors have been set */




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     zbond   bond length used to define each Z-matrix atom */
/*     zang    bond angle used to define each Z-matrix atom */
/*     ztors   angle or torsion used to define Z-matrix atom */
/*     iz      defining atom numbers for each Z-matrix atom */




/*     nothing to do if coordinate type is undefined */

    /* Parameter adjustments */
    --xx;

    /* Function Body */
    if (s_cmp(output_1.coordtype, "NONE", (ftnlen)9, (ftnlen)4) == 0) {
	return 0;
    }

/*     check scaling factors for optimization parameters */

    if (! scales_1.set_scale__) {
	scales_1.set_scale__ = TRUE_;
	if (s_cmp(output_1.coordtype, "CARTESIAN", (ftnlen)9, (ftnlen)9) == 0)
		 {
	    i__1 = atoms_1.n * 3;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		scales_1.scale[i__ - 1] = 1.;
	    }
	} else if (s_cmp(output_1.coordtype, "INTERNAL", (ftnlen)9, (ftnlen)8)
		 == 0) {
	    i__1 = omega_1.nomega;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		scales_1.scale[i__ - 1] = 1.;
	    }
	}
    }

/*     transform optimization parameters back to coordinates */

    if (s_cmp(output_1.coordtype, "CARTESIAN", (ftnlen)9, (ftnlen)9) == 0) {
	nvar = 0;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (usage_1.use[i__ - 1]) {
		++nvar;
		atoms_1.x[i__ - 1] = xx[nvar] / scales_1.scale[nvar - 1];
		++nvar;
		atoms_1.y[i__ - 1] = xx[nvar] / scales_1.scale[nvar - 1];
		++nvar;
		atoms_1.z__[i__ - 1] = xx[nvar] / scales_1.scale[nvar - 1];
	    }
	}
    } else if (s_cmp(output_1.coordtype, "INTERNAL", (ftnlen)9, (ftnlen)8) == 
	    0) {
	i__1 = omega_1.nomega;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    omega_1.dihed[i__ - 1] = xx[i__] / scales_1.scale[i__ - 1];
	    zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] = omega_1.dihed[i__ - 
		    1] * 57.29577951308232088;
	}
    }

/*     get name of archive or intermediate coordinates file */

    iopt = freeunit_();
    if (output_1.cyclesave) {
	if (output_1.archive) {
	    s_copy(optfile, files_1.filename, (ftnlen)120, files_1.leng);
	    suffix_(optfile, "arc", (ftnlen)120, (ftnlen)3);
	    version_(optfile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = optfile;
	    ioin__1.inex = &exist;
	    ioin__1.inopen = 0;
	    ioin__1.innum = 0;
	    ioin__1.innamed = 0;
	    ioin__1.inname = 0;
	    ioin__1.inacc = 0;
	    ioin__1.inseq = 0;
	    ioin__1.indir = 0;
	    ioin__1.infmt = 0;
	    ioin__1.inform = 0;
	    ioin__1.inunf = 0;
	    ioin__1.inrecl = 0;
	    ioin__1.innrec = 0;
	    ioin__1.inblank = 0;
	    f_inqu(&ioin__1);
	    if (exist) {
		openend_(&iopt, optfile, (ftnlen)120);
	    } else {
		o__1.oerr = 0;
		o__1.ounit = iopt;
		o__1.ofnmlen = 120;
		o__1.ofnm = optfile;
		o__1.orl = 0;
		o__1.osta = "new";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
	    }
	} else {
	    lext = 3;
	    numeral_(ncycle, ext, &lext, (ftnlen)7);
/* Writing concatenation */
	    i__2[0] = files_1.leng, a__1[0] = files_1.filename;
	    i__2[1] = 1, a__1[1] = ".";
	    i__2[2] = lext, a__1[2] = ext;
	    s_cat(optfile, a__1, i__2, &c__3, (ftnlen)120);
	    version_(optfile, "new", (ftnlen)120, (ftnlen)3);
	    o__1.oerr = 0;
	    o__1.ounit = iopt;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = optfile;
	    o__1.orl = 0;
	    o__1.osta = "new";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	}
    } else {
	s_copy(optfile, files_1.outfile, (ftnlen)120, (ftnlen)120);
	version_(optfile, "old", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = iopt;
	o__1.ofnmlen = 120;
	o__1.ofnm = optfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = iopt;
	f_rew(&al__1);
    }

/*     update intermediate file with desired coordinate type */

    if (s_cmp(output_1.coordtype, "CARTESIAN", (ftnlen)9, (ftnlen)9) == 0) {
	prtxyz_(&iopt);
    } else if (s_cmp(output_1.coordtype, "INTERNAL", (ftnlen)9, (ftnlen)8) == 
	    0) {
	prtint_(&iopt);
    } else if (s_cmp(output_1.coordtype, "RIGIDBODY", (ftnlen)9, (ftnlen)9) ==
	     0) {
	prtxyz_(&iopt);
    }
    cl__1.cerr = 0;
    cl__1.cunit = iopt;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     send data via external socket communication if desired */

    if (! socket_1.skt_init__ || socket_1.use_socket__) {
	if (s_cmp(output_1.coordtype, "INTERNAL", (ftnlen)9, (ftnlen)8) == 0) 
		{
	    makexyz_();
	}
	sktopt_(ncycle, f);
    }

/*     test for requested termination of the optimization */

    s_copy(endfile, "tinker.end", (ftnlen)120, (ftnlen)10);
    ioin__1.inerr = 0;
    ioin__1.infilen = 120;
    ioin__1.infile = endfile;
    ioin__1.inex = &exist;
    ioin__1.inopen = 0;
    ioin__1.innum = 0;
    ioin__1.innamed = 0;
    ioin__1.inname = 0;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    if (! exist) {
/* Writing concatenation */
	i__3[0] = files_1.leng, a__2[0] = files_1.filename;
	i__3[1] = 4, a__2[1] = ".end";
	s_cat(endfile, a__2, i__3, &c__2, (ftnlen)120);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = endfile;
	ioin__1.inex = &exist;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (exist) {
	    iend = freeunit_();
	    o__1.oerr = 0;
	    o__1.ounit = iend;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = endfile;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	    cl__1.cerr = 0;
	    cl__1.cunit = iend;
	    cl__1.csta = "delete";
	    f_clos(&cl__1);
	}
    }
    if (exist) {
	io___10.ciunit = iounit_1.iout;
	s_wsfe(&io___10);
	e_wsfe();
	fatal_();
    }
    return 0;
} /* optsave_ */

