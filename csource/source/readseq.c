/* readseq.f -- translated by f2c (version 20050501).
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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    char amino[93], nuclz[36], amino1[31], nuclz1[12];
} resdue_;

#define resdue_1 resdue_

struct {
    integer nseq, nchain, ichain[20000]	/* was [2][10000] */, seqtyp[10000];
    char seq[30000], chnnam[10000];
} sequen_;

#define sequen_1 sequen_

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine readseq  --  read biopolymer sequence file  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "readseq" gets a biopolymer sequence containing one or more */
/*     separate chains from an external file; all lines containing */
/*     sequence must begin with the starting sequence number, the */
/*     actual sequence is read from subsequent nonblank characters */


/* Subroutine */ int readseq_(integer *iseq)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 READSEQ  --  Unable to Find the Biopolym"
	    "er\002,\002 Sequence File\002)";
    static char fmt_20[] = "(a120)";

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
	     s_cmp(char *, char *, ftnlen, ftnlen), f_clos(cllist *);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j;
    static logical done;
    static char word[3];
    static integer next;
    extern /* Subroutine */ int fatal_(void);
    static logical exist, opened;
    static char record[120];
    static integer length, number;
    static char letter[1], seqfile[120];
    extern /* Subroutine */ int getnumb_(char *, integer *, integer *, ftnlen)
	    , getword_(char *, char *, integer *, ftnlen, ftnlen), gettext_(
	    char *, char *, integer *, ftnlen, ftnlen), version_(char *, char 
	    *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___6 = { 1, 0, 1, fmt_20, 0 };



#define seq_ref(a_0,a_1) &sequen_1.seq[(a_1)*3 + a_0 - 3]
#define amino_ref(a_0,a_1) &resdue_1.amino[(a_1)*3 + a_0 - 3]
#define nuclz_ref(a_0,a_1) &resdue_1.nuclz[(a_1)*3 + a_0 - 3]
#define ichain_ref(a_1,a_2) sequen_1.ichain[(a_2)*2 + a_1 - 3]



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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




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




/*     open the input file if it has not already been done */

    ioin__1.inerr = 0;
    ioin__1.inunit = *iseq;
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
	i__1[1] = 4, a__1[1] = ".seq";
	s_cat(seqfile, a__1, i__1, &c__2, (ftnlen)120);
	version_(seqfile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = seqfile;
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
	    o__1.ounit = *iseq;
	    o__1.ofnmlen = 120;
	    o__1.ofnm = seqfile;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
	    al__1.aerr = 0;
	    al__1.aunit = *iseq;
	    f_rew(&al__1);
	} else {
	    io___4.ciunit = iounit_1.iout;
	    s_wsfe(&io___4);
	    e_wsfe();
	    fatal_();
	}
    }

/*     zero out the number and type of residues */

    sequen_1.nseq = 0;
    sequen_1.nchain = 0;
    for (i__ = 1; i__ <= 10000; ++i__) {
	s_copy(seq_ref(0, i__), "   ", (ftnlen)3, (ftnlen)3);
    }

/*     read in the biopolymer sequence file */

    while(TRUE_) {
	io___6.ciunit = *iseq;
	i__2 = s_rsfe(&io___6);
	if (i__2 != 0) {
	    goto L30;
	}
	i__2 = do_fio(&c__1, record, (ftnlen)120);
	if (i__2 != 0) {
	    goto L30;
	}
	i__2 = e_rsfe();
	if (i__2 != 0) {
	    goto L30;
	}
	length = trimtext_(record, (ftnlen)120);
	next = 1;
	gettext_(record, letter, &next, (ftnlen)120, (ftnlen)1);
	if (*(unsigned char *)letter >= '0' && *(unsigned char *)letter <= 
		'9') {
	    next = 1;
	    *(unsigned char *)letter = ' ';
	}
	getnumb_(record, &number, &next, (ftnlen)120);
	if (number == 1) {
	    ++sequen_1.nchain;
	    ichain_ref(1, sequen_1.nchain) = sequen_1.nseq + 1;
	    *(unsigned char *)&sequen_1.chnnam[sequen_1.nchain - 1] = *(
		    unsigned char *)letter;
	}
	done = FALSE_;
	while(! done) {
	    getword_(record, word, &next, (ftnlen)120, (ftnlen)3);
	    if (s_cmp(word, "   ", (ftnlen)3, (ftnlen)3) == 0) {
		done = TRUE_;
	    } else {
		++sequen_1.nseq;
		s_copy(seq_ref(0, sequen_1.nseq), word, (ftnlen)3, (ftnlen)3);
	    }
	}
    }
L30:

/*     set the last residue in each sequence chain */

    i__2 = sequen_1.nchain - 1;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ichain_ref(2, i__) = ichain_ref(1, i__ + 1) - 1;
    }
    if (sequen_1.nchain != 0) {
	ichain_ref(2, sequen_1.nchain) = sequen_1.nseq;
    }

/*     find the residue type for each sequence element */

    i__2 = sequen_1.nseq;
    for (i__ = 1; i__ <= i__2; ++i__) {
	sequen_1.seqtyp[i__ - 1] = 0;
	for (j = 1; j <= 31; ++j) {
	    if (s_cmp(seq_ref(0, i__), amino_ref(0, j), (ftnlen)3, (ftnlen)3) 
		    == 0) {
		sequen_1.seqtyp[i__ - 1] = j;
		goto L40;
	    }
	}
	for (j = 1; j <= 12; ++j) {
	    if (s_cmp(seq_ref(0, i__), nuclz_ref(0, j), (ftnlen)3, (ftnlen)3) 
		    == 0) {
		sequen_1.seqtyp[i__ - 1] = j;
		goto L40;
	    }
	}
L40:
	;
    }
    if (! opened) {
	cl__1.cerr = 0;
	cl__1.cunit = *iseq;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    return 0;
} /* readseq_ */

#undef ichain_ref
#undef nuclz_ref
#undef amino_ref
#undef seq_ref


