/* prterr.f -- translated by f2c (version 20050501).
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
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

/* Table of constant values */

static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine prterr  --  output coordinates upon error  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "prterr" writes out a set of coordinates to a disk */
/*     file prior to aborting on a serious error */


/* Subroutine */ int prterr_(void)
{
    /* System generated locals */
    address a__1[2];
    integer i__1[2];
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), s_cmp(char *, char *, ftnlen, ftnlen), f_clos(
	    cllist *);

    /* Local variables */
    extern integer freeunit_(void);
    static char errorfile[120];
    static integer ierr;
    extern /* Subroutine */ int prtint_(integer *), prtxyz_(integer *), 
	    version_(char *, char *, ftnlen, ftnlen);



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
/*     ##  output.i  --  control of coordinate output file format  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     archive    logical flag to save structures in an archive */
/*     noversion  logical flag governing use of filename versions */
/*     overwrite  logical flag to overwrite intermediate files inplace */
/*     cyclesave  logical flag to mark use of numbered cycle files */
/*     coordtype  selects Cartesian, internal, rigid body or none */




/*     write the current coordinates to a file after an error */

    ierr = freeunit_();
/* Writing concatenation */
    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
    i__1[1] = 4, a__1[1] = ".err";
    s_cat(errorfile, a__1, i__1, &c__2, (ftnlen)120);
    version_(errorfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = ierr;
    o__1.ofnmlen = 120;
    o__1.ofnm = errorfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    if (s_cmp(output_1.coordtype, "CARTESIAN", (ftnlen)9, (ftnlen)9) == 0) {
	prtxyz_(&ierr);
    } else if (s_cmp(output_1.coordtype, "INTERNAL", (ftnlen)9, (ftnlen)8) == 
	    0) {
	prtint_(&ierr);
    } else if (s_cmp(output_1.coordtype, "RIGIDBODY", (ftnlen)9, (ftnlen)9) ==
	     0) {
	prtxyz_(&ierr);
    }
    cl__1.cerr = 0;
    cl__1.cunit = ierr;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* prterr_ */

