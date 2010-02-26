/* sybylxyz.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  program sybylxyz  --  Sybyl MOL2 to Cartesian coordinates  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "sybylxyz" takes as input a Sybyl MOL2 coordinates file, */
/*     converts to and then writes out Cartesian coordinates */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Title :  \002,a)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2];
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern integer freeunit_(void);
    static integer ixyz;
    extern /* Subroutine */ int final_(void), getmol2_(void), prtxyz_(integer 
	    *), initial_(void), version_(char *, char *, ftnlen, ftnlen);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };




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




/*     get and read the Sybyl MOL2 format file */

    initial_();
    getmol2_();
    io___1.ciunit = iounit_1.iout;
    s_wsfe(&io___1);
    do_fio(&c__1, titles_1.title, titles_1.ltitle);
    e_wsfe();

/*     write out the Cartesian coordinates file */

    ixyz = freeunit_();
/* Writing concatenation */
    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
    i__1[1] = 4, a__1[1] = ".xyz";
    s_cat(xyzfile, a__1, i__1, &c__2, (ftnlen)120);
    version_(xyzfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = ixyz;
    o__1.ofnmlen = 120;
    o__1.ofnm = xyzfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    prtxyz_(&ixyz);
    cl__1.cerr = 0;
    cl__1.cunit = ixyz;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

/* Main program alias */ int sybylxyz_ () { MAIN__ (); return 0; }
