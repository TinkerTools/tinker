/* xyzpdb.f -- translated by f2c (version 20050501).
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

/* Table of constant values */

static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  program xyzpdb  --  Cartesian to Protein Data Bank file  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "xyzpdb" takes as input a Cartesian coordinates file, */
/*     then converts to and writes out a Protein Data Bank file */


/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    address a__1[2];
    integer i__1[2];
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_rew(alist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int molecule_(void);
    extern integer freeunit_(void);
    static integer ipdb, ixyz;
    extern /* Subroutine */ int field_(void), final_(void), katom_(void), 
	    prtpdb_(integer *), suffix_(char *, char *, ftnlen, ftnlen), 
	    getxyz_(void), makepdb_(void);
    static char pdbfile[120];
    extern /* Subroutine */ int initial_(void), version_(char *, char *, 
	    ftnlen, ftnlen), readxyz_(integer *);
    static char xyzfile[120];



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




/*     get the Cartesian coordinates file for the system */

    initial_();
    getxyz_();

/*     get atomic number of each atom and count the molecules */

    field_();
    katom_();
    molecule_();

/*     open the Protein Data Bank file to be used for output */

    ipdb = freeunit_();
/* Writing concatenation */
    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
    i__1[1] = 4, a__1[1] = ".pdb";
    s_cat(pdbfile, a__1, i__1, &c__2, (ftnlen)120);
    version_(pdbfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = ipdb;
    o__1.ofnmlen = 120;
    o__1.ofnm = pdbfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);

/*     reopen the coordinates file and read the first structure */

    ixyz = freeunit_();
    s_copy(xyzfile, files_1.filename, (ftnlen)120, (ftnlen)120);
    suffix_(xyzfile, "xyz", (ftnlen)120, (ftnlen)3);
    version_(xyzfile, "old", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = ixyz;
    o__1.ofnmlen = 120;
    o__1.ofnm = xyzfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    al__1.aerr = 0;
    al__1.aunit = ixyz;
    f_rew(&al__1);
    readxyz_(&ixyz);

/*     add each successive coordinate frame to the PDB file */

    while(! inform_1.abort) {
	makepdb_();
	prtpdb_(&ipdb);
	readxyz_(&ixyz);
    }

/*     perform any final tasks before program exit */

    cl__1.cerr = 0;
    cl__1.cunit = ipdb;
    cl__1.csta = 0;
    f_clos(&cl__1);
    final_();
    return 0;
} /* MAIN__ */

/* Main program alias */ int xyzpdb_ () { MAIN__ (); return 0; }
