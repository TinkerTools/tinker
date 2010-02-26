/* openend.f -- translated by f2c (version 20050501).
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



/*     ################################################### */
/*     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine openend  --  open a file positioned for append  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "openend" opens a file on a Fortran unit such that the position */
/*     is set to the bottom for appending to the end of the file */

/*     note this routine is system dependent since the Fortran 90 */
/*     standard is not supported by many Fortran 77 compilers; only */
/*     one of the various implementations below should be activated */
/*     by removing comment characters */


/* Subroutine */ int openend_(integer *iunit, char *name__, ftnlen name_len)
{
    /* System generated locals */
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *);



/*     standard Fortran 90, unavailable in some Fortran 77 compilers */

/*     open (unit=iunit,file=name,status='old',position='append') */

/*     common extension supported by many Fortran 77 compilers */

    o__1.oerr = 0;
    o__1.ounit = *iunit;
    o__1.ofnmlen = 120;
    o__1.ofnm = name__;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = "append";
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);

/*     some Fortran 77 compilers open files for append by default */

/*     open (unit=iunit,file=name,status='old') */

/*     manually read to the end of file, slow but always correct */

/*     open (unit=iunit,file=name,status='old') */
/*     dowhile (.true.) */
/*        read (iunit,10,err=20,end=20) */
/*  10    format () */
/*     end do */
/*  20 continue */
    return 0;
} /* openend_ */

