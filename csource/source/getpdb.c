/* getpdb.f -- translated by f2c (version 20050501).
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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

/* Table of constant values */

static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine getpdb  --  get a Protein Data Bank file  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "getpdb" asks for a Protein Data Bank file name, */
/*     then reads in the coordinates file */


/* Subroutine */ int getpdb_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter Protein Data Bank File Name :  "
	    "\002,$)";
    static char fmt_20[] = "(a120)";

    /* System generated locals */
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *)
	    , do_fio(integer *, char *, ftnlen), e_rsfe(void), f_open(olist *)
	    , f_rew(alist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int basefile_(char *, ftnlen);
    extern integer freeunit_(void);
    static integer ipdb;
    static logical exist;
    extern /* Subroutine */ int suffix_(char *, char *, ftnlen, ftnlen), 
	    readpdb_(integer *);
    static char pdbfile[120];
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), version_(
	    char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_20, 0 };




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




/*     try to get a filename from the command line arguments */

    nextarg_(pdbfile, &exist, (ftnlen)120);
    if (exist) {
	basefile_(pdbfile, (ftnlen)120);
	suffix_(pdbfile, "pdb", (ftnlen)120, (ftnlen)3);
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
    }

/*     ask for the user specified input structure filename */

    while(! exist) {
	io___3.ciunit = iounit_1.iout;
	s_wsfe(&io___3);
	e_wsfe();
	io___4.ciunit = iounit_1.input;
	s_rsfe(&io___4);
	do_fio(&c__1, pdbfile, (ftnlen)120);
	e_rsfe();
	basefile_(pdbfile, (ftnlen)120);
	suffix_(pdbfile, "pdb", (ftnlen)120, (ftnlen)3);
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
    }

/*     first open and then read the PDB coordinates file */

    ipdb = freeunit_();
    o__1.oerr = 0;
    o__1.ounit = ipdb;
    o__1.ofnmlen = 120;
    o__1.ofnm = pdbfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    al__1.aerr = 0;
    al__1.aunit = ipdb;
    f_rew(&al__1);
    readpdb_(&ipdb);
    cl__1.cerr = 0;
    cl__1.cunit = ipdb;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* getpdb_ */

