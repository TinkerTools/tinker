/* getint.f -- translated by f2c (version 20050501).
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
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

/* Table of constant values */

static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine getint  --  get internal coordinate structure  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "getint" asks for an internal coordinate file name, then reads */
/*     the internal coordinates and computes Cartesian coordinates */


/* Subroutine */ int getint_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter Internal Coordinate File Name : "
	    " \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(/,\002 GETINT  --  Internal Coordinates File"
	    "\002,\002 does not Contain Any Atoms\002)";

    /* System generated locals */
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *)
	    , do_fio(integer *, char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *), f_rew(alist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int basefile_(char *, ftnlen);
    extern integer freeunit_(void);
    static integer izmt;
    extern /* Subroutine */ int fatal_(void);
    static logical clash, exist;
    extern /* Subroutine */ int suffix_(char *, char *, ftnlen, ftnlen), 
	    chkxyz_(logical *), readint_(integer *), connect_(void);
    static char intfile[120];
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), version_(
	    char *, char *, ftnlen, ftnlen), makexyz_(void);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_30, 0 };




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
/*     ##  output.i  --  control of coordinate output file format  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     archive    logical flag to save structures in an archive */
/*     noversion  logical flag governing use of filename versions */
/*     overwrite  logical flag to overwrite intermediate files inplace */
/*     cyclesave  logical flag to mark use of numbered cycle files */
/*     coordtype  selects Cartesian, internal, rigid body or none */




/*     try to get a filename from the command line arguments */

    nextarg_(intfile, &exist, (ftnlen)120);
    if (exist) {
	basefile_(intfile, (ftnlen)120);
	suffix_(intfile, "int", (ftnlen)120, (ftnlen)3);
	version_(intfile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = intfile;
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
	do_fio(&c__1, intfile, (ftnlen)120);
	e_rsfe();
	basefile_(intfile, (ftnlen)120);
	suffix_(intfile, "int", (ftnlen)120, (ftnlen)3);
	version_(intfile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = intfile;
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

/*     first open and then read the internal coordinates file */

    s_copy(output_1.coordtype, "INTERNAL", (ftnlen)9, (ftnlen)8);
    izmt = freeunit_();
    o__1.oerr = 0;
    o__1.ounit = izmt;
    o__1.ofnmlen = 120;
    o__1.ofnm = intfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    al__1.aerr = 0;
    al__1.aunit = izmt;
    f_rew(&al__1);
    readint_(&izmt);
    cl__1.cerr = 0;
    cl__1.cunit = izmt;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     quit if the internal coordinates file contains no atoms */

    if (inform_1.abort) {
	io___6.ciunit = iounit_1.iout;
	s_wsfe(&io___6);
	e_wsfe();
	fatal_();
    }

/*     convert internal to Cartesian coordinates */

    connect_();
    makexyz_();

/*     check for atoms with identical coordinates */

    clash = FALSE_;
    chkxyz_(&clash);
    return 0;
} /* getint_ */

