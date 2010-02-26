/* xyzint.f -- translated by f2c (version 20050501).
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
    integer ltitle;
    char title[120];
} titles_;

#define titles_1 titles_

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  program xyzint  --  Cartesian to internal coordinates  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "xyzint" takes as input a Cartesian coordinates file, then */
/*     converts to and writes out an internal coordinates file */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Title :  \002,a)";
    static char fmt_20[] = "(/,\002 Template (T), Dihedrals (D), Manual (M"
	    ")\002,\002 or Automatic [A] :  \002,$)";
    static char fmt_30[] = "(a120)";
    static char fmt_40[] = "(/,\002 XYZINT  --  Warning, Template File wa"
	    "s\002,\002 not Found\002)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2];
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_rsfe(cilist *), e_rsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_inqu(inlist *), f_open(olist *), f_rew(alist *), f_clos(cllist *
	    );

    /* Local variables */
    extern integer freeunit_(void);
    static integer mode, next, izmt;
    extern /* Subroutine */ int final_(void);
    static logical exist;
    static char record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char answer[1];
    extern /* Subroutine */ int prtint_(integer *), getxyz_(void), readint_(
	    integer *), makeint_(integer *), initial_(void);
    static char intfile[120];
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), gettext_(
	    char *, char *, integer *, ftnlen, ftnlen), version_(char *, char 
	    *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_40, 0 };




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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  titles.i  --  title for the current molecular system  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     ltitle   length in characters of the nonblank title string */
/*     title    title used to describe the current structure */




/*     get and read the Cartesian coordinates file */

    initial_();
    getxyz_();
    io___1.ciunit = iounit_1.iout;
    s_wsfe(&io___1);
    do_fio(&c__1, titles_1.title, titles_1.ltitle);
    e_wsfe();

/*     set the mode for conversion to internal coordinates */

    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___4.ciunit = iounit_1.iout;
	s_wsfe(&io___4);
	e_wsfe();
	io___5.ciunit = iounit_1.input;
	s_rsfe(&io___5);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    mode = 0;
    if (*(unsigned char *)answer == 'M') {
	mode = 1;
    } else if (*(unsigned char *)answer == 'T') {
	mode = 2;
/* Writing concatenation */
	i__1[0] = files_1.leng, a__1[0] = files_1.filename;
	i__1[1] = 4, a__1[1] = ".int";
	s_cat(intfile, a__1, i__1, &c__2, (ftnlen)120);
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
	if (exist) {
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
	} else {
	    mode = 0;
	    io___11.ciunit = iounit_1.iout;
	    s_wsfe(&io___11);
	    e_wsfe();
	}
    } else if (*(unsigned char *)answer == 'D') {
	mode = 3;
    }

/*     convert from Cartesian to internal coordinates */

    makeint_(&mode);

/*     write out the internal coordinates file */

    izmt = freeunit_();
/* Writing concatenation */
    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
    i__1[1] = 4, a__1[1] = ".int";
    s_cat(intfile, a__1, i__1, &c__2, (ftnlen)120);
    version_(intfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = izmt;
    o__1.ofnmlen = 120;
    o__1.ofnm = intfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    prtint_(&izmt);
    cl__1.cerr = 0;
    cl__1.cunit = izmt;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

/* Main program alias */ int xyzint_ () { MAIN__ (); return 0; }
