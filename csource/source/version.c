/* version.f -- translated by f2c (version 20050501).
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

struct {
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

/* Table of constant values */

static integer c__6 = 6;
static integer c__5 = 5;
static integer c__4 = 4;
static integer c__3 = 3;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine version  --  create version number for file  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "version" checks the name of a file about to be opened; if */
/*     if "old" status is passed, the name of the highest current */
/*     version is returned; if "new" status is passed the filename */
/*     of the next available unused version is generated */


/* Subroutine */ int version_(char *filename, char *status, ftnlen 
	filename_len, ftnlen status_len)
{
    /* Initialized data */

    static char digit[1*10] = "0" "1" "2" "3" "4" "5" "6" "7" "8" "9";

    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter File Name for Coordinate Output "
	    ":  \002,$)";
    static char fmt_20[] = "(a120)";

    /* System generated locals */
    address a__1[6], a__2[5], a__3[4], a__4[3];
    integer i__1, i__2[6], i__3[5], i__4[4], i__5[3];
    inlist ioin__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_inqu(inlist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void);

    /* Local variables */
    static integer thousand;
    extern integer trimtext_(char *, ftnlen);
    static integer i__, leng, ones, tens;
    static logical exist;
    static char oldfile[120];
    static integer hundred;
    static char newfile[120];
    extern /* Subroutine */ int lowcase_(char *, ftnlen), nextarg_(char *, 
	    logical *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_20, 0 };




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




/*     process the filename and status variables */

    lowcase_(status, (ftnlen)3);
    leng = trimtext_(filename, (ftnlen)120);

/*     check for attempt to access an unversioned file */

    if (! output_1.noversion) {
	i__1 = leng - 2;
	if (s_cmp(filename + i__1, "_1", leng - i__1, (ftnlen)2) == 0) {
	    s_copy(filename, filename, (ftnlen)120, leng - 2);
	    if (s_cmp(status, "old", (ftnlen)3, (ftnlen)3) == 0) {
		return 0;
	    }
	}
    }

/*     no change is needed if the file doesn't exist */

    exist = FALSE_;
    if (leng != 0) {
	ioin__1.inerr = 0;
	ioin__1.infilen = leng;
	ioin__1.infile = filename;
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
    if (! exist) {
	return 0;
    }

/*     set initial values for the current and next versions */

    s_copy(newfile, filename, (ftnlen)120, (ftnlen)120);
    s_copy(oldfile, filename, (ftnlen)120, (ftnlen)120);

/*     append an artificial version number to the filename; */
/*     currently handles up to 10000 versions of a file */

    if (! output_1.noversion) {
	i__ = 1;
	while(exist) {
	    ++i__;
	    s_copy(oldfile, newfile, (ftnlen)120, (ftnlen)120);
	    thousand = i__ / 1000;
	    hundred = (i__ - thousand * 1000) / 100;
	    tens = (i__ - thousand * 1000 - hundred * 100) / 10;
	    ones = i__ - thousand * 1000 - hundred * 100 - tens * 10;
	    if (thousand != 0) {
/* Writing concatenation */
		i__2[0] = leng, a__1[0] = filename;
		i__2[1] = 1, a__1[1] = "_";
		i__2[2] = 1, a__1[2] = digit + thousand;
		i__2[3] = 1, a__1[3] = digit + hundred;
		i__2[4] = 1, a__1[4] = digit + tens;
		i__2[5] = 1, a__1[5] = digit + ones;
		s_cat(newfile, a__1, i__2, &c__6, (ftnlen)120);
	    } else if (hundred != 0) {
/* Writing concatenation */
		i__3[0] = leng, a__2[0] = filename;
		i__3[1] = 1, a__2[1] = "_";
		i__3[2] = 1, a__2[2] = digit + hundred;
		i__3[3] = 1, a__2[3] = digit + tens;
		i__3[4] = 1, a__2[4] = digit + ones;
		s_cat(newfile, a__2, i__3, &c__5, (ftnlen)120);
	    } else if (tens != 0) {
/* Writing concatenation */
		i__4[0] = leng, a__3[0] = filename;
		i__4[1] = 1, a__3[1] = "_";
		i__4[2] = 1, a__3[2] = digit + tens;
		i__4[3] = 1, a__3[3] = digit + ones;
		s_cat(newfile, a__3, i__4, &c__4, (ftnlen)120);
	    } else {
/* Writing concatenation */
		i__5[0] = leng, a__4[0] = filename;
		i__5[1] = 1, a__4[1] = "_";
		i__5[2] = 1, a__4[2] = digit + ones;
		s_cat(newfile, a__4, i__5, &c__3, (ftnlen)120);
	    }
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = newfile;
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
    }

/*     set the file name based on the requested status */

    if (s_cmp(status, "old", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(filename, oldfile, (ftnlen)120, (ftnlen)120);
    } else if (s_cmp(status, "new", (ftnlen)3, (ftnlen)3) == 0) {
	s_copy(filename, newfile, (ftnlen)120, (ftnlen)120);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = filename;
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
	    nextarg_(filename, &exist, (ftnlen)120);
	    if (exist) {
		ioin__1.inerr = 0;
		ioin__1.infilen = 120;
		ioin__1.infile = filename;
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
	    } else {
		exist = TRUE_;
	    }
	    while(exist) {
		io___11.ciunit = iounit_1.iout;
		s_wsfe(&io___11);
		e_wsfe();
		io___12.ciunit = iounit_1.input;
		s_rsfe(&io___12);
		do_fio(&c__1, filename, (ftnlen)120);
		e_rsfe();
		ioin__1.inerr = 0;
		ioin__1.infilen = 120;
		ioin__1.infile = filename;
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
	}
    }
    return 0;
} /* version_ */

