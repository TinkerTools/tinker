/* prtseq.f -- translated by f2c (version 20050501).
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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine prtseq  --  output of biopolymer sequence  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "prtseq" writes out a biopolymer sequence to an external */
/*     disk file with 15 residues per line and distinct chains */
/*     separated by blank lines */


/* Subroutine */ int prtseq_(integer *iseq)
{
    /* Format strings */
    static char fmt_10[] = "()";
    static char fmt_20[] = "(3x,a1,i6,1x,15(1x,a3))";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3;
    olist o__1;
    cllist cl__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), s_wsfe(cilist *), e_wsfe(void), do_fio(integer *,
	     char *, ftnlen), f_clos(cllist *);

    /* Local variables */
    static integer i__, k, smin, smax, size, stop, start;
    static logical opened;
    static char letter[1], seqfile[120];
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_20, 0 };



#define seq_ref(a_0,a_1) &sequen_1.seq[(a_1)*3 + a_0 - 3]
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




/*     open output unit if not already done */

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
	version_(seqfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = *iseq;
	o__1.ofnmlen = 120;
	o__1.ofnm = seqfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }

/*     write out a three-letter code sequence file */

    i__2 = sequen_1.nchain;
    for (i__ = 1; i__ <= i__2; ++i__) {
	*(unsigned char *)letter = *(unsigned char *)&sequen_1.chnnam[i__ - 1]
		;
	start = ichain_ref(1, i__);
	stop = ichain_ref(2, i__);
	size = stop - start + 1;
	smax = 0;
	while(smax < size) {
	    smin = smax + 1;
	    smax += 15;
	    smax = min(smax,size);
	    if (i__ != 1 && smin == 1) {
		io___10.ciunit = *iseq;
		s_wsfe(&io___10);
		e_wsfe();
	    }
	    io___11.ciunit = *iseq;
	    s_wsfe(&io___11);
	    do_fio(&c__1, letter, (ftnlen)1);
	    do_fio(&c__1, (char *)&smin, (ftnlen)sizeof(integer));
	    i__3 = smax;
	    for (k = smin; k <= i__3; ++k) {
		do_fio(&c__1, seq_ref(0, k + start - 1), (ftnlen)3);
	    }
	    e_wsfe();
	}
    }
    if (! opened) {
	cl__1.cerr = 0;
	cl__1.cunit = *iseq;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    return 0;
} /* prtseq_ */

#undef ichain_ref
#undef seq_ref


