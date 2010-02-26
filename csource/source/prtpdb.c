/* prtpdb.f -- translated by f2c (version 20050501).
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
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

struct {
    doublereal xpdb[25000], ypdb[25000], zpdb[25000];
    integer npdb, resnum[25000], npdb12[25000], ipdb12[200000]	/* was [8][
	    25000] */, pdblist[25000];
    char pdbtyp[150000], atmnam[100000], resnam[75000], chntyp[20], altsym[1],
	     instyp[20];
} pdb_;

#define pdb_1 pdb_

struct {
    integer nseq, nchain, ichain[20000]	/* was [2][10000] */, seqtyp[10000];
    char seq[30000], chnnam[10000];
} sequen_;

#define sequen_1 sequen_

struct {
    integer ltitle;
    char title[120];
} titles_;

#define titles_1 titles_

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine prtpdb  --  output of Protein Data Bank file  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "prtpdb" writes out a set of Protein Data Bank coordinates */
/*     to an external disk file */


/* Subroutine */ int prtpdb_(integer *ipdb)
{
    /* Format strings */
    static char fmt_10[] = "(\002HEADER\002,/,\002COMPND\002,/,\002SOURCE"
	    "\002)";
    static char fmt_20[] = "(\002HEADER\002,4x,a,/,\002COMPND\002,/,\002SOUR"
	    "CE\002)";
    static char fmt_30[] = "(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)";
    static char fmt_40[] = "(\002CONECT\002,5i5)";
    static char fmt_50[] = "(\002END\002)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3;
    olist o__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), s_wsfe(cilist *), e_wsfe(void), do_fio(integer *,
	     char *, ftnlen), s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, k, stop;
    static char chain[1*10000];
    static integer resid[10000], start;
    static logical opened;
    static char pdbfile[120], chnname[1], resname[3];
    static integer resnumb;
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_50, 0 };



#define ipdb12_ref(a_1,a_2) pdb_1.ipdb12[(a_2)*8 + a_1 - 9]
#define ichain_ref(a_1,a_2) sequen_1.ichain[(a_2)*2 + a_1 - 3]
#define atmnam_ref(a_0,a_1) &pdb_1.atmnam[(a_1)*4 + a_0 - 4]
#define resnam_ref(a_0,a_1) &pdb_1.resnam[(a_1)*3 + a_0 - 3]
#define pdbtyp_ref(a_0,a_1) &pdb_1.pdbtyp[(a_1)*6 + a_0 - 6]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  pdb.i  --  definition of a Protein Data Bank structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xpdb      x-coordinate of each atom stored in PDB format */
/*     ypdb      y-coordinate of each atom stored in PDB format */
/*     zpdb      z-coordinate of each atom stored in PDB format */
/*     npdb      number of atoms stored in Protein Data Bank format */
/*     resnum    number of the residue to which each atom belongs */
/*     npdb12    number of atoms directly bonded to each CONECT atom */
/*     ipdb12    atom numbers of atoms connected to each CONECT atom */
/*     pdblist   list of the Protein Data Bank atom number of each atom */
/*     pdbtyp    Protein Data Bank record type assigned to each atom */
/*     atmnam    Protein Data Bank atom name assigned to each atom */
/*     resnam    Protein Data Bank residue name assigned to each atom */
/*     chntyp    string with PDB chain identifiers to be included */
/*     altsym    string with PDB alternate locations to be included */
/*     instyp    string with PDB insertion records to be included */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  sequen.i  --  sequence information for a biopolymer  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nseq     total number of residues in biopolymer sequences */
/*     nchain   number of separate biopolymer sequence chains */
/*     ichain   first and last residue in each biopolymer chain */
/*     seqtyp   residue type for each residue in the sequence */
/*     seq      three-letter code for each residue in the sequence */
/*     chnnam   one-letter identifier for each sequence chain */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  titles.i  --  title for the current molecular system  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     ltitle   length in characters of the nonblank title string */
/*     title    title used to describe the current structure */




/*     open output unit if not already done */

    ioin__1.inerr = 0;
    ioin__1.inunit = *ipdb;
    ioin__1.infile = 0;
    ioin__1.inex = 0;
    ioin__1.inopen = &opened;
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
    if (! opened) {
/* Writing concatenation */
	i__1[0] = files_1.leng, a__1[0] = files_1.filename;
	i__1[1] = 4, a__1[1] = ".pdb";
	s_cat(pdbfile, a__1, i__1, &c__2, (ftnlen)120);
	version_(pdbfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = *ipdb;
	o__1.ofnmlen = 120;
	o__1.ofnm = pdbfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }

/*     write out the header lines and the title */

    if (titles_1.ltitle == 0) {
	io___3.ciunit = *ipdb;
	s_wsfe(&io___3);
	e_wsfe();
    } else {
	io___4.ciunit = *ipdb;
	s_wsfe(&io___4);
	do_fio(&c__1, titles_1.title, titles_1.ltitle);
	e_wsfe();
    }

/*     find the chain name and chain position for each residue */

    i__2 = sequen_1.nchain;
    for (i__ = 1; i__ <= i__2; ++i__) {
	start = ichain_ref(1, i__);
	stop = ichain_ref(2, i__);
	i__3 = stop;
	for (k = start; k <= i__3; ++k) {
	    resid[k - 1] = k - start + 1;
	    *(unsigned char *)&chain[k - 1] = *(unsigned char *)&
		    sequen_1.chnnam[i__ - 1];
	}
    }

/*     change some TINKER residue names to match PDB standards */

    i__2 = pdb_1.npdb;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (s_cmp(resnam_ref(0, i__), "CYX", (ftnlen)3, (ftnlen)3) == 0) {
	    s_copy(resnam_ref(0, i__), "CYS", (ftnlen)3, (ftnlen)3);
	}
	if (s_cmp(resnam_ref(0, i__), "HID", (ftnlen)3, (ftnlen)3) == 0) {
	    s_copy(resnam_ref(0, i__), "HIS", (ftnlen)3, (ftnlen)3);
	}
	if (s_cmp(resnam_ref(0, i__), "HIE", (ftnlen)3, (ftnlen)3) == 0) {
	    s_copy(resnam_ref(0, i__), "HIS", (ftnlen)3, (ftnlen)3);
	}
	if (s_cmp(resnam_ref(0, i__), "HIP", (ftnlen)3, (ftnlen)3) == 0) {
	    s_copy(resnam_ref(0, i__), "HIS", (ftnlen)3, (ftnlen)3);
	}
    }

/*     next, write the coordinates for each PDB atom */

    i__2 = pdb_1.npdb;
    for (i__ = 1; i__ <= i__2; ++i__) {
	s_copy(resname, resnam_ref(0, i__), (ftnlen)3, (ftnlen)3);
	if (s_cmp(resname + 1, "  ", (ftnlen)2, (ftnlen)2) == 0) {
/* Writing concatenation */
	    i__1[0] = 2, a__1[0] = "  ";
	    i__1[1] = 1, a__1[1] = resname;
	    s_cat(resname, a__1, i__1, &c__2, (ftnlen)3);
	}
	if (*(unsigned char *)&resname[2] == ' ') {
/* Writing concatenation */
	    i__1[0] = 1, a__1[0] = " ";
	    i__1[1] = 2, a__1[1] = resname;
	    s_cat(resname, a__1, i__1, &c__2, (ftnlen)3);
	}
	if (s_cmp(pdbtyp_ref(0, i__), "ATOM  ", (ftnlen)6, (ftnlen)6) == 0) {
	    resnumb = resid[pdb_1.resnum[i__ - 1] - 1];
	    *(unsigned char *)chnname = *(unsigned char *)&chain[pdb_1.resnum[
		    i__ - 1] - 1];
	} else {
	    resnumb = pdb_1.resnum[i__ - 1];
	    *(unsigned char *)chnname = ' ';
	}
	io___14.ciunit = *ipdb;
	s_wsfe(&io___14);
	do_fio(&c__1, pdbtyp_ref(0, i__), (ftnlen)6);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, atmnam_ref(0, i__), (ftnlen)4);
	do_fio(&c__1, resname, (ftnlen)3);
	do_fio(&c__1, chnname, (ftnlen)1);
	do_fio(&c__1, (char *)&resnumb, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&pdb_1.xpdb[i__ - 1], (ftnlen)sizeof(doublereal)
		);
	do_fio(&c__1, (char *)&pdb_1.ypdb[i__ - 1], (ftnlen)sizeof(doublereal)
		);
	do_fio(&c__1, (char *)&pdb_1.zpdb[i__ - 1], (ftnlen)sizeof(doublereal)
		);
	e_wsfe();
    }

/*     finally, write any connectivity records for PDB atoms */

    i__2 = pdb_1.npdb;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (pdb_1.npdb12[i__ - 1] != 0) {
	    io___15.ciunit = *ipdb;
	    s_wsfe(&io___15);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    i__3 = pdb_1.npdb12[i__ - 1];
	    for (k = 1; k <= i__3; ++k) {
		do_fio(&c__1, (char *)&ipdb12_ref(k, i__), (ftnlen)sizeof(
			integer));
	    }
	    e_wsfe();
	}
    }
    io___16.ciunit = *ipdb;
    s_wsfe(&io___16);
    e_wsfe();
/*     close (unit=ipdb) */
    return 0;
} /* prtpdb_ */

#undef pdbtyp_ref
#undef resnam_ref
#undef atmnam_ref
#undef ichain_ref
#undef ipdb12_ref


