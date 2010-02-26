/* readpdb.f -- translated by f2c (version 20050501).
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
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

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

struct {
    char amino[93], nuclz[36], amino1[31], nuclz1[12];
} resdue_;

#define resdue_1 resdue_

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static integer c__10000 = 10000;
static integer c__25000 = 25000;
static integer c__3 = 3;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine readpdb  --  input of Protein Data Bank file  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "readpdb" gets a set of Protein Data Bank coordinates */
/*     from an external disk file */


/* Subroutine */ int readpdb_(integer *ipdb)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 READPDB  --  Unable to Find the Protei"
	    "n\002,\002 Data Bank File\002)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3)";
    static char fmt_40[] = "(/,\002 READPDB  --  The Maximum of\002,i6,\002 "
	    "Residues has been Exceeded\002)";
    static char fmt_60[] = "(i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3)";
    static char fmt_80[] = "(/,\002 READPDB  --  The Maximum of\002,i6,\002 "
	    "Atoms has been Exceeded\002)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_rew(alist *), s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void),
	     s_cmp(char *, char *, ftnlen, ftnlen), s_rsfi(icilist *), e_rsfi(
	    void), i_indx(char *, char *, ftnlen, ftnlen), f_clos(cllist *);

    /* Local variables */
    static char namelast[3];
    extern integer trimtext_(char *, ftnlen);
    static integer i__;
    static doublereal xx, yy, zz;
    static integer nres;
    static char chain[1];
    extern /* Subroutine */ int fatal_(void);
    static logical exist, opened;
    static char chnatm[1*25000], altloc[1];
    static integer serial;
    static char remark[6], insert[1], letter[1], record[120], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), fixpdb_(char *, char 
	    *, ftnlen, ftnlen);
    static char pdbfile[120];
    extern /* Subroutine */ int scanpdb_(integer *);
    static char atmname[4], resname[3], chnlast[1];
    static integer residue, reslast;
    static char inslast[1];
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___9 = { 1, 0, 1, fmt_20, 0 };
    static icilist io___13 = { 0, string, 0, fmt_30, 120, 1 };
    static cilist io___25 = { 0, 0, 0, fmt_40, 0 };
    static icilist io___27 = { 0, string, 0, fmt_60, 120, 1 };
    static cilist io___28 = { 0, 0, 0, fmt_80, 0 };



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




/*     open the input file if it has not already been done */

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
	version_(pdbfile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = pdbfile;
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
	    o__1.oerr = 0;
	    o__1.ounit = *ipdb;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = pdbfile;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	    al__1.aerr = 0;
	    al__1.aunit = *ipdb;
	    f_rew(&al__1);
	} else {
	    io___4.ciunit = iounit_1.iout;
	    s_wsfe(&io___4);
	    e_wsfe();
	    fatal_();
	}
    }

/*     get alternate sites, chains and insertions to be used */

    scanpdb_(ipdb);

/*     initialize title, atom and residue counters and name */

    s_copy(titles_1.title, " ", (ftnlen)120, (ftnlen)1);
    titles_1.ltitle = 0;
    pdb_1.npdb = 0;
    nres = 0;
    reslast = 10000;
    s_copy(namelast, "   ", (ftnlen)3, (ftnlen)3);
    *(unsigned char *)chnlast = ' ';

/*     process individual atoms from the Protein Data Bank file */

    while(TRUE_) {
	io___9.ciunit = *ipdb;
	i__2 = s_rsfe(&io___9);
	if (i__2 != 0) {
	    goto L90;
	}
	i__2 = do_fio(&c__1, record, (ftnlen)120);
	if (i__2 != 0) {
	    goto L90;
	}
	i__2 = e_rsfe();
	if (i__2 != 0) {
	    goto L90;
	}
	upcase_(record, (ftnlen)120);
	s_copy(remark, record, (ftnlen)6, (ftnlen)6);
	if (s_cmp(remark, "HEADER", (ftnlen)6, (ftnlen)6) == 0) {
	    s_copy(titles_1.title, record + 10, (ftnlen)120, (ftnlen)60);
	    titles_1.ltitle = trimtext_(titles_1.title, (ftnlen)120);
	} else if (s_cmp(remark, "TITLE ", (ftnlen)6, (ftnlen)6) == 0) {
	    if (titles_1.ltitle == 0) {
		s_copy(titles_1.title, record + 10, (ftnlen)120, (ftnlen)60);
		titles_1.ltitle = trimtext_(titles_1.title, (ftnlen)120);
	    }
	} else if (s_cmp(remark, "ATOM  ", (ftnlen)6, (ftnlen)6) == 0) {
	    s_copy(string, record + 6, (ftnlen)120, (ftnlen)114);
	    s_rsfi(&io___13);
	    do_fio(&c__1, (char *)&serial, (ftnlen)sizeof(integer));
	    do_fio(&c__1, atmname, (ftnlen)4);
	    do_fio(&c__1, altloc, (ftnlen)1);
	    do_fio(&c__1, resname, (ftnlen)3);
	    do_fio(&c__1, chain, (ftnlen)1);
	    do_fio(&c__1, (char *)&residue, (ftnlen)sizeof(integer));
	    do_fio(&c__1, insert, (ftnlen)1);
	    do_fio(&c__1, (char *)&xx, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&yy, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&zz, (ftnlen)sizeof(doublereal));
	    e_rsfi();
	    if (s_cmp(resname, "  ", (ftnlen)2, (ftnlen)2) == 0) {
		s_copy(resname, resname + 2, (ftnlen)3, (ftnlen)1);
	    }
	    if (*(unsigned char *)resname == ' ') {
		s_copy(resname, resname + 1, (ftnlen)3, (ftnlen)2);
	    }
	    if (*(unsigned char *)chain != ' ' && i_indx(pdb_1.chntyp, chain, 
		    (ftnlen)20, (ftnlen)1) == 0) {
		goto L50;
	    }
	    if (*(unsigned char *)altloc != ' ' && *(unsigned char *)altloc !=
		     *(unsigned char *)pdb_1.altsym) {
		goto L50;
	    }
	    if (*(unsigned char *)insert != ' ' && i_indx(pdb_1.instyp, 
		    insert, (ftnlen)20, (ftnlen)1) == 0) {
		goto L50;
	    }
	    fixpdb_(resname, atmname, (ftnlen)3, (ftnlen)4);
	    if (residue != reslast || s_cmp(resname, namelast, (ftnlen)3, (
		    ftnlen)3) != 0 || *(unsigned char *)chain != *(unsigned 
		    char *)chnlast || *(unsigned char *)insert != *(unsigned 
		    char *)inslast) {
		++nres;
		reslast = residue;
		s_copy(namelast, resname, (ftnlen)3, (ftnlen)3);
		*(unsigned char *)chnlast = *(unsigned char *)chain;
		*(unsigned char *)inslast = *(unsigned char *)insert;
		if (nres > 10000) {
		    io___25.ciunit = iounit_1.iout;
		    s_wsfe(&io___25);
		    do_fio(&c__1, (char *)&c__10000, (ftnlen)sizeof(integer));
		    e_wsfe();
		    fatal_();
		}
	    }
	    ++pdb_1.npdb;
	    pdb_1.xpdb[pdb_1.npdb - 1] = xx;
	    pdb_1.ypdb[pdb_1.npdb - 1] = yy;
	    pdb_1.zpdb[pdb_1.npdb - 1] = zz;
	    s_copy(pdbtyp_ref(0, pdb_1.npdb), remark, (ftnlen)6, (ftnlen)6);
	    s_copy(atmnam_ref(0, pdb_1.npdb), atmname, (ftnlen)4, (ftnlen)4);
	    s_copy(resnam_ref(0, pdb_1.npdb), resname, (ftnlen)3, (ftnlen)3);
	    pdb_1.resnum[pdb_1.npdb - 1] = nres;
	    *(unsigned char *)&chnatm[pdb_1.npdb - 1] = *(unsigned char *)
		    chain;
L50:
	    ;
	} else if (s_cmp(remark, "HETATM", (ftnlen)6, (ftnlen)6) == 0) {
	    s_copy(string, record + 6, (ftnlen)120, (ftnlen)114);
	    s_rsfi(&io___27);
	    do_fio(&c__1, (char *)&serial, (ftnlen)sizeof(integer));
	    do_fio(&c__1, atmname, (ftnlen)4);
	    do_fio(&c__1, altloc, (ftnlen)1);
	    do_fio(&c__1, resname, (ftnlen)3);
	    do_fio(&c__1, chain, (ftnlen)1);
	    do_fio(&c__1, (char *)&residue, (ftnlen)sizeof(integer));
	    do_fio(&c__1, insert, (ftnlen)1);
	    do_fio(&c__1, (char *)&xx, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&yy, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&zz, (ftnlen)sizeof(doublereal));
	    e_rsfi();
	    if (s_cmp(resname, "  ", (ftnlen)2, (ftnlen)2) == 0) {
		s_copy(resname, resname + 2, (ftnlen)3, (ftnlen)1);
	    }
	    if (*(unsigned char *)resname == ' ') {
		s_copy(resname, resname + 1, (ftnlen)3, (ftnlen)2);
	    }
	    if (*(unsigned char *)chain != ' ' && i_indx(pdb_1.chntyp, chain, 
		    (ftnlen)20, (ftnlen)1) == 0) {
		goto L70;
	    }
	    if (*(unsigned char *)altloc != ' ' && *(unsigned char *)altloc !=
		     *(unsigned char *)pdb_1.altsym) {
		goto L70;
	    }
	    if (*(unsigned char *)insert != ' ' && i_indx(pdb_1.instyp, 
		    insert, (ftnlen)20, (ftnlen)1) == 0) {
		goto L70;
	    }
	    fixpdb_(resname, atmname, (ftnlen)3, (ftnlen)4);
	    ++pdb_1.npdb;
	    pdb_1.xpdb[pdb_1.npdb - 1] = xx;
	    pdb_1.ypdb[pdb_1.npdb - 1] = yy;
	    pdb_1.zpdb[pdb_1.npdb - 1] = zz;
	    s_copy(pdbtyp_ref(0, pdb_1.npdb), remark, (ftnlen)6, (ftnlen)6);
	    s_copy(atmnam_ref(0, pdb_1.npdb), atmname, (ftnlen)4, (ftnlen)4);
	    s_copy(resnam_ref(0, pdb_1.npdb), resname, (ftnlen)3, (ftnlen)3);
	    pdb_1.resnum[pdb_1.npdb - 1] = 0;
	    *(unsigned char *)&chnatm[pdb_1.npdb - 1] = *(unsigned char *)
		    chain;
L70:
	    ;
	} else if (s_cmp(remark, "ENDMDL", (ftnlen)6, (ftnlen)6) == 0) {
	    goto L90;
	} else if (s_cmp(remark, "END   ", (ftnlen)6, (ftnlen)6) == 0) {
	    goto L90;
	}
	if (pdb_1.npdb > 25000) {
	    io___28.ciunit = iounit_1.iout;
	    s_wsfe(&io___28);
	    do_fio(&c__1, (char *)&c__25000, (ftnlen)sizeof(integer));
	    e_wsfe();
	    fatal_();
	}
    }
L90:

/*     set the total sequence length and chain termini information */

    if (pdb_1.npdb != 0) {
	sequen_1.nseq = pdb_1.npdb;
	sequen_1.nchain = 0;
	*(unsigned char *)chnlast = '#';
	i__2 = pdb_1.npdb;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (s_cmp(pdbtyp_ref(0, i__), "ATOM  ", (ftnlen)6, (ftnlen)6) == 
		    0) {
		*(unsigned char *)letter = *(unsigned char *)&chnatm[i__ - 1];
		if (*(unsigned char *)letter != *(unsigned char *)chnlast) {
		    ++sequen_1.nchain;
		    ichain_ref(1, sequen_1.nchain) = pdb_1.resnum[i__ - 1];
		    *(unsigned char *)&sequen_1.chnnam[sequen_1.nchain - 1] = 
			    *(unsigned char *)letter;
		    *(unsigned char *)chnlast = *(unsigned char *)letter;
		} else {
		    ichain_ref(2, sequen_1.nchain) = pdb_1.resnum[i__ - 1];
		}
	    }
	}
    }

/*     close the PDB file and quit if there are no coordinates */

    if (pdb_1.npdb == 0) {
	inform_1.abort = TRUE_;
    }
    if (! opened) {
	cl__1.cerr = 0;
	cl__1.cunit = *ipdb;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    return 0;
} /* readpdb_ */

#undef pdbtyp_ref
#undef resnam_ref
#undef atmnam_ref
#undef ichain_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine scanpdb  --  PDB chains, alternates and inserts  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "scanpdb" reads the first model in a Protein Data Bank file and */
/*     sets chains, alternate sites and insertion records to be used */


/* Subroutine */ int scanpdb_(integer *ipdb)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* Format strings */
    static char fmt_10[] = "(a120)";
    static char fmt_20[] = "(10x,a1,4x,a1,4x,a1)";
    static char fmt_40[] = "(/,\002 Enter the Chain Names to Include\002,"
	    "\002 (\002,a,\002) :  \002,$)";
    static char fmt_50[] = "(a20)";
    static char fmt_60[] = "(/,\002 Enter a Set of Alternate Atom Location"
	    "s\002,\002 from (\002,a,\002) :  \002,$)";
    static char fmt_70[] = "(a120)";
    static char fmt_80[] = "(/,\002 Enter the Insert Records to Include\002"
	    ",\002 (\002,a,\002) :  \002,$)";
    static char fmt_90[] = "(a20)";

    /* System generated locals */
    address a__1[2], a__2[3];
    integer i__1, i__2[2], i__3[3];
    alist al__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void),
	     s_cmp(char *, char *, ftnlen, ftnlen), s_rsfi(icilist *), e_rsfi(
	    void), i_indx(char *, char *, ftnlen, ftnlen), f_rew(alist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__;
    static logical done;
    static integer nalt, nins, next;
    static char text[20], chain[1], blank[20];
    static logical exist;
    static char altloc[1];
    static integer length;
    static char remark[6], record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char insert[1], string[120], alttyp[20], chnlast[1], chntemp[20], 
	    altlast[1];
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen);
    static char inslast[1], instemp[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static cilist io___40 = { 1, 0, 1, fmt_10, 0 };
    static icilist io___44 = { 0, string, 0, fmt_20, 120, 1 };
    static cilist io___52 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___56 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___57 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___59 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___60 = { 0, 0, 0, fmt_90, 0 };




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




/*     only proceed if this routine was not already used */

    if (! first) {
	return 0;
    }
    first = FALSE_;

/*     initialize chain, alternate site and insertion lists */

    sequen_1.nchain = 0;
    nalt = 0;
    nins = 0;
    *(unsigned char *)chnlast = '#';
    *(unsigned char *)altlast = '#';
    *(unsigned char *)inslast = '#';
    s_copy(blank, "                    ", (ftnlen)20, (ftnlen)20);
    s_copy(pdb_1.chntyp, "####################", (ftnlen)20, (ftnlen)20);
    *(unsigned char *)pdb_1.altsym = ' ';
    s_copy(alttyp, blank, (ftnlen)20, (ftnlen)20);
    s_copy(pdb_1.instyp, blank, (ftnlen)20, (ftnlen)20);

/*     scan for multiple chains, alternate locations and inserts */

    done = FALSE_;
    while(! done) {
	io___40.ciunit = *ipdb;
	i__1 = s_rsfe(&io___40);
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = do_fio(&c__1, record, (ftnlen)120);
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L30;
	}
	upcase_(record, (ftnlen)120);
	s_copy(remark, record, (ftnlen)6, (ftnlen)6);
	s_copy(string, record + 6, (ftnlen)120, (ftnlen)114);
	if (s_cmp(remark, "ATOM  ", (ftnlen)6, (ftnlen)6) == 0 || s_cmp(
		remark, "HETATM", (ftnlen)6, (ftnlen)6) == 0) {
	    s_rsfi(&io___44);
	    do_fio(&c__1, altloc, (ftnlen)1);
	    do_fio(&c__1, chain, (ftnlen)1);
	    do_fio(&c__1, insert, (ftnlen)1);
	    e_rsfi();
	    if (*(unsigned char *)chain != *(unsigned char *)chnlast) {
		if (i_indx(pdb_1.chntyp, chain, (ftnlen)20, (ftnlen)1) == 0) {
		    ++sequen_1.nchain;
		    *(unsigned char *)&pdb_1.chntyp[sequen_1.nchain - 1] = *(
			    unsigned char *)chain;
		    *(unsigned char *)chnlast = *(unsigned char *)chain;
		}
	    }
	    if (*(unsigned char *)altloc != *(unsigned char *)altlast) {
		if (i_indx(alttyp, altloc, (ftnlen)20, (ftnlen)1) == 0) {
		    ++nalt;
		    *(unsigned char *)&alttyp[nalt - 1] = *(unsigned char *)
			    altloc;
		    *(unsigned char *)altlast = *(unsigned char *)altloc;
		}
	    }
	    if (*(unsigned char *)insert != *(unsigned char *)inslast) {
		if (i_indx(pdb_1.instyp, insert, (ftnlen)20, (ftnlen)1) == 0) 
			{
		    ++nins;
		    *(unsigned char *)&pdb_1.instyp[nins - 1] = *(unsigned 
			    char *)insert;
		    *(unsigned char *)inslast = *(unsigned char *)insert;
		}
	    }
	} else if (s_cmp(remark, "ENDMDL", (ftnlen)6, (ftnlen)6) == 0) {
	    done = TRUE_;
	} else if (s_cmp(remark, "END   ", (ftnlen)6, (ftnlen)6) == 0) {
	    done = TRUE_;
	}
    }
L30:
    al__1.aerr = 0;
    al__1.aunit = *ipdb;
    f_rew(&al__1);

/*     find out which of the multiple chains will be used */

    if (sequen_1.nchain > 1) {
	nextarg_(chntemp, &exist, (ftnlen)20);
	if (! exist) {
	    s_copy(chntemp, blank, (ftnlen)20, (ftnlen)20);
	    if (*(unsigned char *)pdb_1.chntyp == ' ') {
		s_copy(string, "BLANK", (ftnlen)120, (ftnlen)5);
		length = 5;
	    } else {
		*(unsigned char *)string = *(unsigned char *)pdb_1.chntyp;
		length = 1;
	    }
	    i__1 = sequen_1.nchain;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		if (*(unsigned char *)&pdb_1.chntyp[i__ - 1] == ' ') {
/* Writing concatenation */
		    i__2[0] = length, a__1[0] = string;
		    i__2[1] = 6, a__1[1] = " BLANK";
		    s_cat(string, a__1, i__2, &c__2, (ftnlen)120);
		    length += 6;
		} else {
/* Writing concatenation */
		    i__3[0] = length, a__2[0] = string;
		    i__3[1] = 1, a__2[1] = " ";
		    i__3[2] = 1, a__2[2] = pdb_1.chntyp + (i__ - 1);
		    s_cat(string, a__2, i__3, &c__3, (ftnlen)120);
		    length += 2;
		}
	    }
/* Writing concatenation */
	    i__2[0] = length, a__1[0] = string;
	    i__2[1] = 6, a__1[1] = " [ALL]";
	    s_cat(string, a__1, i__2, &c__2, (ftnlen)120);
	    length += 6;
	    io___52.ciunit = iounit_1.iout;
	    s_wsfe(&io___52);
	    do_fio(&c__1, string, length);
	    e_wsfe();
	    io___53.ciunit = iounit_1.input;
	    s_rsfe(&io___53);
	    do_fio(&c__1, chntemp, (ftnlen)20);
	    e_rsfe();
	}
	upcase_(chntemp, (ftnlen)20);
	next = 1;
	gettext_(chntemp, text, &next, (ftnlen)20, (ftnlen)20);
	if (s_cmp(text, blank, (ftnlen)20, (ftnlen)20) == 0 || s_cmp(text, 
		"ALL ", (ftnlen)20, (ftnlen)4) == 0) {
	    s_copy(pdb_1.chntyp, pdb_1.chntyp, (ftnlen)20, sequen_1.nchain);
	} else {
	    sequen_1.nchain = 1;
	    s_copy(pdb_1.chntyp, chntemp, (ftnlen)20, (ftnlen)1);
	}
    }

/*     find out which of the alternate locations will be used */

    if (nalt > 0) {
	nextarg_(pdb_1.altsym, &exist, (ftnlen)1);
	if (! exist) {
/* Writing concatenation */
	    i__3[0] = 1, a__2[0] = "[";
	    i__3[1] = 1, a__2[1] = alttyp;
	    i__3[2] = 1, a__2[2] = "]";
	    s_cat(string, a__2, i__3, &c__3, (ftnlen)3);
	    length = 3;
	    i__1 = nalt;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* Writing concatenation */
		i__3[0] = length, a__2[0] = string;
		i__3[1] = 1, a__2[1] = " ";
		i__3[2] = 1, a__2[2] = alttyp + (i__ - 1);
		s_cat(string, a__2, i__3, &c__3, (ftnlen)120);
		length += 2;
	    }
	    io___56.ciunit = iounit_1.iout;
	    s_wsfe(&io___56);
	    do_fio(&c__1, string, length);
	    e_wsfe();
	    io___57.ciunit = iounit_1.input;
	    s_rsfe(&io___57);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    next = 1;
	    gettext_(record, pdb_1.altsym, &next, (ftnlen)120, (ftnlen)1);
	}
	if (*(unsigned char *)pdb_1.altsym == ' ') {
	    *(unsigned char *)pdb_1.altsym = *(unsigned char *)alttyp;
	}
	upcase_(pdb_1.altsym, (ftnlen)1);
    }

/*     find out which of the insert records will be used */

    if (nins > 0) {
	nextarg_(instemp, &exist, (ftnlen)20);
	if (! exist) {
	    s_copy(instemp, blank, (ftnlen)20, (ftnlen)20);
	    *(unsigned char *)string = *(unsigned char *)pdb_1.instyp;
	    length = 1;
	    i__1 = nins;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* Writing concatenation */
		i__3[0] = length, a__2[0] = string;
		i__3[1] = 1, a__2[1] = " ";
		i__3[2] = 1, a__2[2] = pdb_1.instyp + (i__ - 1);
		s_cat(string, a__2, i__3, &c__3, (ftnlen)120);
		length += 2;
	    }
/* Writing concatenation */
	    i__2[0] = length, a__1[0] = string;
	    i__2[1] = 11, a__1[1] = " [ALL] NONE";
	    s_cat(string, a__1, i__2, &c__2, (ftnlen)120);
	    length += 11;
	    io___59.ciunit = iounit_1.iout;
	    s_wsfe(&io___59);
	    do_fio(&c__1, string, length);
	    e_wsfe();
	    io___60.ciunit = iounit_1.input;
	    s_rsfe(&io___60);
	    do_fio(&c__1, instemp, (ftnlen)20);
	    e_rsfe();
	}
	upcase_(instemp, (ftnlen)20);
	next = 1;
	gettext_(instemp, text, &next, (ftnlen)20, (ftnlen)20);
	if (s_cmp(text, blank, (ftnlen)20, (ftnlen)20) == 0 || s_cmp(text, 
		"ALL ", (ftnlen)20, (ftnlen)4) == 0) {
	    s_copy(pdb_1.instyp, pdb_1.instyp, (ftnlen)20, nins);
	} else if (s_cmp(text, "NONE ", (ftnlen)20, (ftnlen)5) == 0) {
	    s_copy(pdb_1.instyp, blank, (ftnlen)20, (ftnlen)20);
	} else {
	    s_copy(pdb_1.instyp, instemp, (ftnlen)20, (ftnlen)20);
	}
    }
    return 0;
} /* scanpdb_ */



/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine fixpdb  --  standard PDB atom and residue names  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "fixpdb" corrects problems with PDB files by converting residue */
/*     and atom names to the standard forms used by TINKER */


/* Subroutine */ int fixpdb_(char *resname, char *atmname, ftnlen resname_len,
	 ftnlen atmname_len)
{
    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static char restype[7];


#define amino_ref(a_0,a_1) &resdue_1.amino[(a_1)*3 + a_0 - 3]
#define nuclz_ref(a_0,a_1) &resdue_1.nuclz[(a_1)*3 + a_0 - 3]



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
/*     ##  resdue.i  --  standard biopolymer residue abbreviations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     amino    three-letter abbreviations for amino acids types */
/*     nuclz    three-letter abbreviations for nucleic acids types */
/*     amino1   one-letter abbreviations for amino acids types */
/*     nuclz1   one-letter abbreviations for nucleic acids types */




/*     convert traditional 3-letter base names to PDB names */

    if (s_cmp(resname, "ADE", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "A  ", (ftnlen)3, (ftnlen)3);
    }
    if (s_cmp(resname, "GUA", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "G  ", (ftnlen)3, (ftnlen)3);
    }
    if (s_cmp(resname, "CYT", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "C  ", (ftnlen)3, (ftnlen)3);
    }
    if (s_cmp(resname, "URA", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "U  ", (ftnlen)3, (ftnlen)3);
    }
    if (s_cmp(resname, "DAD", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "DA ", (ftnlen)3, (ftnlen)3);
    }
    if (s_cmp(resname, "DGU", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "DG ", (ftnlen)3, (ftnlen)3);
    }
    if (s_cmp(resname, "DCY", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "DC ", (ftnlen)3, (ftnlen)3);
    }
    if (s_cmp(resname, "THY", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "DT ", (ftnlen)3, (ftnlen)3);
    }

/*     convert unusual names for protonated histidine residues */

    if (s_cmp(resname, "HSD", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "HID", (ftnlen)3, (ftnlen)3);
    }
    if (s_cmp(resname, "HSE", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "HIE", (ftnlen)3, (ftnlen)3);
    }
    if (s_cmp(resname, "HSP", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "HIS", (ftnlen)3, (ftnlen)3);
    }
    if (s_cmp(resname, "HSH", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "HIS", (ftnlen)3, (ftnlen)3);
    }
    if (s_cmp(resname, "HIP", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "HIS", (ftnlen)3, (ftnlen)3);
    }
    if (s_cmp(resname, "HIH", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "HIS", (ftnlen)3, (ftnlen)3);
    }

/*     convert unusual names for terminal capping residues */

    if (s_cmp(resname, "NMA", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "NME", (ftnlen)3, (ftnlen)3);
    }

/*     convert nonstandard names for water molecules */

    if (s_cmp(resname, "H2O", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "HOH", (ftnlen)3, (ftnlen)3);
    }
    if (s_cmp(resname, "WAT", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(resname, "HOH", (ftnlen)3, (ftnlen)3);
    }

/*     decide whether residue is protein or nucleic acid */

    s_copy(restype, "UNKNOWN", (ftnlen)7, (ftnlen)7);
    for (i__ = 1; i__ <= 31; ++i__) {
	if (s_cmp(resname, amino_ref(0, i__), (ftnlen)3, (ftnlen)3) == 0) {
	    s_copy(restype, "PROTEIN", (ftnlen)7, (ftnlen)7);
	}
    }
    for (i__ = 1; i__ <= 12; ++i__) {
	if (s_cmp(resname, nuclz_ref(0, i__), (ftnlen)3, (ftnlen)3) == 0) {
	    s_copy(restype, "NUCLEIC", (ftnlen)7, (ftnlen)7);
	}
    }

/*     convert any generically used unusual atom names */

    if (s_cmp(atmname, " HN ", (ftnlen)4, (ftnlen)4) == 0) {
	s_copy(atmname, " H  ", (ftnlen)4, (ftnlen)4);
    }

/*     convert unusual names in protein terminal residues */

    if (s_cmp(restype, "PROTEIN", (ftnlen)7, (ftnlen)7) == 0) {
	if (s_cmp(atmname, "1H  ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HN1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HT1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2H  ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HN2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HT2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3H  ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H3 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HN3", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H3 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HT3", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H3 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " O1 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " O  ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " OT1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " O  ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "OCT1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " O  ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " O2 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " OXT", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " OT2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " OXT", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "OCT2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " OXT", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " OT ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " OXT", (ftnlen)4, (ftnlen)4);
	}
    }

/*     convert unusual names common to many nucleotides */

    if (s_cmp(restype, "NUCLEIC", (ftnlen)7, (ftnlen)7) == 0) {
	if (s_cmp(atmname, " O1P", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " OP1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " O2P", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " OP2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " O3P", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " OP3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HOP", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HOP2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HOP", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HOP3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " C1*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " C1'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " C2*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " C2'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " C3*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " C3'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " C4*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " C4'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " C5*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " C5'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " O2*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " O2'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " O3*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " O3'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " O4*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " O4'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " O5*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " O5'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " H1*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H1'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " H2*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1H2*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2H2*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "H2''", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " H3*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H3'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " H4*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H4'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1H5*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H5'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2H5*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "H5''", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HO*", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HO2'", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " H3T", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HO3'", (ftnlen)4, (ftnlen)4);
	}
    }

/*     glycine residue  (GLY) */

    if (s_cmp(resname, "GLY", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HA ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HA2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HA1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HA3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HA ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HA3", (ftnlen)4, (ftnlen)4);
	}

/*     alanine residue  (ALA) */

    } else if (s_cmp(resname, "ALA", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}

/*     valine residue  (VAL) */

    } else if (s_cmp(resname, "VAL", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG11", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG12", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG13", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HG2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG21", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HG2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG22", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HG2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG23", (ftnlen)4, (ftnlen)4);
	}

/*     leucine residue  (LEU) */

    } else if (s_cmp(resname, "LEU", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HD1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD11", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HD1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD12", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HD1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD13", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HD2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD21", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HD2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD22", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HD2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD23", (ftnlen)4, (ftnlen)4);
	}

/*     isoleucine residue  (ILE) */

    } else if (s_cmp(resname, "ILE", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, " CD ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " CD1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG12", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HG11", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG13", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG13", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HG2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG21", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HG2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG22", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HG2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG23", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HD1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD11", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HD1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD11", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HD1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD12", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HD2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD12", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HD1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD13", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HD3", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD13", (ftnlen)4, (ftnlen)4);
	}

/*     serine residue  (SER) */

    } else if (s_cmp(resname, "SER", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, " OG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " OG ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HOG", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG ", (ftnlen)4, (ftnlen)4);
	}

/*     threonine residue  (THR) */

    } else if (s_cmp(resname, "THR", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, " OG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " OG1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " CG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " CG2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HOG", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HOG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HG2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG21", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HG2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG22", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HG2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HG23", (ftnlen)4, (ftnlen)4);
	}

/*     cysteine residue  (CYS) */

    } else if (s_cmp(resname, "CYS", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, " SG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " SG ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HSG", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG ", (ftnlen)4, (ftnlen)4);
	}

/*     proline residue  (PRO) */

    } else if (s_cmp(resname, "PRO", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HD ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HD1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HD ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HT1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HT2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H3 ", (ftnlen)4, (ftnlen)4);
	}

/*     phenylalanine residue  (PHE) */

    } else if (s_cmp(resname, "PHE", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}

/*     tyrosine residue  (TYR) */

    } else if (s_cmp(resname, "TYR", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, " HOH", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HH ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}

/*     tryptophan residue  (TRP) */

    } else if (s_cmp(resname, "TRP", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HNE", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE1", (ftnlen)4, (ftnlen)4);
	}

/*     histidine (HD and HE) residue  (HIS) */

    } else if (s_cmp(resname, "HIS", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HD ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HE ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HND", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HND1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HNE", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HNE2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE2", (ftnlen)4, (ftnlen)4);
	}

/*     histidine (HD only) residue  (HID) */

    } else if (s_cmp(resname, "HID", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HD ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HE ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HND", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HND1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD1", (ftnlen)4, (ftnlen)4);
	}

/*     histidine (HE only) residue  (HIE) */

    } else if (s_cmp(resname, "HIE", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HD ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HE ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HNE", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HNE2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE2", (ftnlen)4, (ftnlen)4);
	}

/*     aspartic acid residue  (ASP) */

    } else if (s_cmp(resname, "ASP", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}

/*     asparagine residue  (ASN) */

    } else if (s_cmp(resname, "ASN", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, " OD ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " OD1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " ND ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " ND2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HD2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD21", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HND1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD21", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HD2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD22", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HND2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HD22", (ftnlen)4, (ftnlen)4);
	}

/*     glutamic acid residue  (GLU) */

    } else if (s_cmp(resname, "GLU", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}

/*     glutamine residue  (GLN) */

    } else if (s_cmp(resname, "GLN", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, " OE ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " OE1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " NE ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " NE2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HE2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HE21", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HNE1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HE21", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HE2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HE22", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HNE2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HE22", (ftnlen)4, (ftnlen)4);
	}

/*     methionine residue  (MET) */

    } else if (s_cmp(resname, "MET", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HE ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HE ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HE ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE3", (ftnlen)4, (ftnlen)4);
	}

/*     lysine residue  (LYS) */

    } else if (s_cmp(resname, "LYS", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HD ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HD1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HD ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HE ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HE1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HE ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HZ ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HZ1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HNZ1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HZ1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HZ ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HZ2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HNZ2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HZ2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HZ ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HZ3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HNZ3", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HZ3", (ftnlen)4, (ftnlen)4);
	}

/*     arginine residue  (ARG) */

    } else if (s_cmp(resname, "ARG", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HD ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HD1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HD ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HH1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HH11", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HN11", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HH11", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HH1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HH12", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HN12", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HH12", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HH2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HH21", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HN21", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HH21", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HH2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HH22", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HN22", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HH22", (ftnlen)4, (ftnlen)4);
	}

/*     ornithine residue  (ORN) */

    } else if (s_cmp(resname, "ORN", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HD ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HD1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HD ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HD3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HE ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HNE1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE1", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HE ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HNE2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HE ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HNE3", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HE3", (ftnlen)4, (ftnlen)4);
	}

/*     methylalanine residue  (AIB) */

    } else if (s_cmp(resname, "AIB", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HB11", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HB12", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HB13", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HB2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HB21", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HB22", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HB2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, "HB23", (ftnlen)4, (ftnlen)4);
	}

/*     pyroglutamic acid residue  (PCA) */

    } else if (s_cmp(resname, "PCA", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HB1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HB ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HB3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG2", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HG1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HG ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " HG3", (ftnlen)4, (ftnlen)4);
	}

/*     N-terminal acetyl residue  (ACE) */

    } else if (s_cmp(resname, "ACE", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, " CY ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " C  ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " CAY", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " CH3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " CA ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " CH3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " OY ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " O  ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1H  ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HY1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HH31", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2H  ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HY2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HH32", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3H  ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H3 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HY3", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H3 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HH33", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H3 ", (ftnlen)4, (ftnlen)4);
	}

/*     N-terminal formyl residue  (FOR) */

    } else if (s_cmp(resname, "FOR", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, " CY ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " C  ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " OY ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " O  ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HY ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H  ", (ftnlen)4, (ftnlen)4);
	}

/*     C-terminal N-methylamide residue  (NME) */

    } else if (s_cmp(resname, "NME", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, " NT ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " N  ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " CT ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " CH3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " CAT", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " CH3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " CA ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " CH3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HNT", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H  ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1H  ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1HA ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HT1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HH31", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2H  ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2HA ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HT2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HH32", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3H  ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H3 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3HA ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H3 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HT3", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H3 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "HH33", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H3 ", (ftnlen)4, (ftnlen)4);
	}

/*     C-terminal amide residue  (NH2) */

    } else if (s_cmp(resname, "NH2", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, " NT ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " N  ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1H  ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2H  ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HT1", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, " HT2", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
	}

/*     adenosine residue  (A) */

    } else if (s_cmp(resname, "A  ", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1H6 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H61", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2H6 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H62", (ftnlen)4, (ftnlen)4);
	}

/*     guanosine residue  (G) */

    } else if (s_cmp(resname, "G  ", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1H2 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H21", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2H2 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H22", (ftnlen)4, (ftnlen)4);
	}

/*     cytidine residue  (C) */

    } else if (s_cmp(resname, "C  ", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1H4 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H41", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2H4 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H42", (ftnlen)4, (ftnlen)4);
	}

/*     deoxyadenosine residue  (DA) */

    } else if (s_cmp(resname, "DA ", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1H6 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H61", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2H6 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H62", (ftnlen)4, (ftnlen)4);
	}

/*     deoxyguanosine residue  (DG) */

    } else if (s_cmp(resname, "DG ", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1H2 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H21", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2H2 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H22", (ftnlen)4, (ftnlen)4);
	}

/*     deoxycytidine residue  (DC) */

    } else if (s_cmp(resname, "DC ", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, "1H4 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H41", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2H4 ", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H42", (ftnlen)4, (ftnlen)4);
	}

/*     deoxythymidine residue  (DT) */

    } else if (s_cmp(resname, "DT ", (ftnlen)3, (ftnlen)3) == 0) {
	if (s_cmp(atmname, " C5M", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " C7 ", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "1H5M", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H71", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "2H5M", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H72", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(atmname, "3H5M", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(atmname, " H73", (ftnlen)4, (ftnlen)4);
	}
    }
    return 0;
} /* fixpdb_ */

#undef nuclz_ref
#undef amino_ref


