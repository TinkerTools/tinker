/* document.f -- translated by f2c (version 20050501).
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

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  program document  --  make documentation lists from source  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "document" generates a formatted description of all the code */
/*     modules or common blocks, an index of routines called by each */
/*     source code module, a listing of all valid keywords, a list of */
/*     include file dependencies as needed by a Unix-style Makefile, */
/*     or a formatted force field parameter set summary */

/*     note the logical variable "wiki" should be set true to make */
/*     output suitable for inclusion in the TINKER User's Guide */
/*     under MediaWiki */


/* Main program */ int MAIN__(void)
{
    /* Initialized data */

    static char fname[20*71] = "ADJACENT            " "ANGGUESS            " 
	    "ANORM               " "BETACF              " "BETAI            "
	    "   " "BMAX                " "BNDERR              " "BNDGUESS    "
	    "        " "CHIRER              " "CJKM                " "D1D2   "
	    "             " "DEPTH               " "DIST2               " 
	    "DOT                 " "ENERGY              " "ERF              "
	    "   " "ERFC                " "ERFINV              " "FREEUNIT    "
	    "        " "GAMMLN              " "GDA2                " "GEOMETR"
	    "Y            " "INITERR             " "INVBETA             " 
	    "LOCERR              " "MAXWELL             " "MCM1             "
	    "   " "MCMSTEP             " "MIDERR              " "MINIMIZ1    "
	    "        " "MINIROT1            " "MINRIGID1           " "NEWTON1"
	    "             " "NEWTROT1            " "NEXTTEXT            " 
	    "NORMAL              " "NUMBER              " "OPBGUESS         "
	    "   " "OPTIMIZ1            " "OPTIROT1            " "OPTRIGID1   "
	    "        " "PATH1               " "PAULING1            " "POTFIT1"
	    "             " "POTNRG              " "PRECISE             " 
	    "PRIORITY            " "PROPERTY            " "PSS1             "
	    "   " "PSSRGD1             " "PSSROT1             " "PTINCY      "
	    "        " "RANDOM              " "RMSFIT              " "ROTANG "
	    "             " "ROTCHECK            " "SADDLE1             " 
	    "SCAN1               " "SIGMOID             " "SNIFFER1         "
	    "   " "TORSER              " "TOTERR              " "TRANSIT     "
	    "        " "TRIMTEXT            " "TRIPLE              " "VALFIT1"
	    "             " "VALRMS              " "VDWERR              " 
	    "VECANG              " "WATSON1             " "XTALLAT1         "
	    "   ";

    /* Format strings */
    static char fmt_10[] = "(/,\002 The TINKER Document Facility can Provide"
	    " :\002,//,4x,\002(1) List of Routines from a Source File\002,/,4"
	    "x,\002(2) List of Calls made by each Routine\002,/,4x,\002(3) Li"
	    "st of Common Blocks from Source\002,/,4x,\002(4) List of the TIN"
	    "KER Option Keywords\002,/,4x,\002(5) List of Include File Depend"
	    "encies\002,/,4x,\002(6) Documentation from a Parameter File\002)";
    static char fmt_30[] = "(/,\002 Enter the Number of the Desired Choice :"
	    "  \002,$)";
    static char fmt_40[] = "(i10)";
    static char fmt_50[] = "(/,\002 Enter Name of Source Code Listing File :"
	    "  \002,$)";
    static char fmt_60[] = "(a120)";
    static char fmt_70[] = "(a120)";
    static char fmt_80[] = "(///)";
    static char fmt_90[] = "(a120)";
    static char fmt_110[] = "(\002'''\002,a,\002'''\002,/)";
    static char fmt_120[] = "(a,/)";
    static char fmt_130[] = "(a)";
    static char fmt_140[] = "(a)";
    static char fmt_150[] = "()";
    static char fmt_160[] = "(/,\002 Source Documentation Written To:  \002,"
	    "a)";
    static char fmt_170[] = "(a120)";
    static char fmt_190[] = "(a,/)";
    static char fmt_200[] = "(a)";
    static char fmt_210[] = "(5x,a)";
    static char fmt_220[] = "(/,\002 Source Documentation Written To:  \002,"
	    "a)";
    static char fmt_230[] = "(a120)";
    static char fmt_240[] = "(///)";
    static char fmt_250[] = "(a120)";
    static char fmt_270[] = "(\002'''\002,a,\002'''\002,/)";
    static char fmt_280[] = "(a,/)";
    static char fmt_290[] = "(a,/)";
    static char fmt_300[] = "(a)";
    static char fmt_310[] = "()";
    static char fmt_320[] = "(/,\002 Source Documentation Written To:  \002,"
	    "a)";
    static char fmt_330[] = "(a120)";
    static char fmt_350[] = "(a)";
    static char fmt_360[] = "(/,\002 Keyword Listing Written To:  \002,a)";
    static char fmt_370[] = "(a120)";
    static char fmt_390[] = "(/,\002 File Dependencies in Makefile Format "
	    ":\002,/)";
    static char fmt_400[] = "(a)";
    static char fmt_410[] = "(/,\002 Parameter Listing Written To:  \002,a)";

    /* System generated locals */
    address a__1[3], a__2[2];
    integer i__1, i__2[3], i__3, i__4, i__5[2];
    olist o__1;
    cllist cl__1;
    inlist ioin__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsli(icilist *), do_lio(integer 
	    *, integer *, char *, ftnlen), e_rsli(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void), f_inqu(inlist *),
	     f_open(olist *), s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen),
	     s_copy(char *, char *, ftnlen, ftnlen);
    integer f_clos(cllist *), i_indx(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern integer freeunit_(void), trimtext_(char *, ftnlen), nexttext_(char 
	    *, ftnlen);
    static integer i__, j, k;
    static char key[20*1000];
    static integer idoc, mode, leng;
    static logical done;
    static char info[120*100*1000];
    static integer link[1000], isrc, last;
    static logical wiki;
    static integer nkey, next;
    extern /* Subroutine */ int sort6_(integer *, char *, ftnlen), sort7_(
	    integer *, char *, integer *, ftnlen);
    static char field[2048];
    extern /* Subroutine */ int final_(void);
    static integer nline[1000];
    extern /* Subroutine */ int sort10_(integer *, char *, ftnlen);
    static logical exist;
    static integer nunit, start;
    static char fname1[20], fname2[20], record[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char module[20];
    extern /* Subroutine */ int getprm_(void), suffix_(char *, char *, ftnlen,
	     ftnlen);
    static char string[120];
    extern /* Subroutine */ int prtprm_(integer *);
    static char docfile[120], srcfile[120];
    extern /* Subroutine */ int initial_(void), lowcase_(char *, ftnlen), 
	    nextarg_(char *, logical *, ftnlen), getword_(char *, char *, 
	    integer *, ftnlen, ftnlen);
    static char keylast[20], keyword[20], routine[120*1000];
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen), 
	    gettext_(char *, char *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static icilist io___7 = { 1, string, 1, 0, 120, 1 };
    static cilist io___8 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___9 = { 1, 0, 0, fmt_40, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___15 = { 1, 0, 1, fmt_70, 0 };
    static cilist io___21 = { 1, 0, 1, fmt_80, 0 };
    static cilist io___24 = { 1, 0, 1, fmt_90, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___40 = { 1, 0, 1, fmt_170, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___46 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___48 = { 1, 0, 1, fmt_230, 0 };
    static cilist io___50 = { 1, 0, 1, fmt_240, 0 };
    static cilist io___51 = { 1, 0, 1, fmt_250, 0 };
    static cilist io___52 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___54 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___55 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___56 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___57 = { 0, 0, 0, fmt_320, 0 };
    static cilist io___59 = { 1, 0, 1, fmt_330, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_350, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_360, 0 };
    static cilist io___64 = { 1, 0, 1, fmt_370, 0 };
    static cilist io___65 = { 0, 0, 0, fmt_390, 0 };
    static cilist io___66 = { 0, 0, 0, fmt_400, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_410, 0 };



#define key_ref(a_0,a_1) &key[(a_1)*20 + a_0 - 20]
#define info_ref(a_0,a_1,a_2) &info[((a_2)*100 + a_1)*120 + a_0 - 12120]
#define fname_ref(a_0,a_1) &fname[(a_1)*20 + a_0 - 20]
#define routine_ref(a_0,a_1) &routine[(a_1)*120 + a_0 - 120]



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



/*     list of the Fortran functions in the TINKER package */



/*     set flag to format for TINKER User's Guide under MediaWiki */

    wiki = TRUE_;

/*     find out what documentation the user wants to generate */

    initial_();
    io___3.ciunit = iounit_1.iout;
    s_wsfe(&io___3);
    e_wsfe();
    mode = 0;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___7);
	if (i__1 != 0) {
	    goto L20;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&mode, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L20;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L20;
	}
    }
L20:
    while(mode < 1 || mode > 6) {
	io___8.ciunit = iounit_1.iout;
	s_wsfe(&io___8);
	e_wsfe();
	io___9.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___9);
	if (i__1 != 0) {
	    goto L20;
	}
	i__1 = do_fio(&c__1, (char *)&mode, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L20;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L20;
	}
    }

/*     get the filename and open the input source code file */

    if (mode != 6) {
	isrc = freeunit_();
	nextarg_(srcfile, &exist, (ftnlen)120);
	if (exist) {
	    suffix_(srcfile, "txt", (ftnlen)120, (ftnlen)3);
	    version_(srcfile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = srcfile;
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
	while(! exist) {
	    io___12.ciunit = iounit_1.iout;
	    s_wsfe(&io___12);
	    e_wsfe();
	    io___13.ciunit = iounit_1.input;
	    s_rsfe(&io___13);
	    do_fio(&c__1, srcfile, (ftnlen)120);
	    e_rsfe();
	    suffix_(srcfile, "txt", (ftnlen)120, (ftnlen)3);
	    version_(srcfile, "old", (ftnlen)120, (ftnlen)3);
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 120;
	    ioin__1.infile = srcfile;
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
	o__1.oerr = 0;
	o__1.ounit = isrc;
	o__1.ofnmlen = 120;
	o__1.ofnm = srcfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }

/*     get a list of routines and descriptions from source code */

    if (mode == 1) {
	nunit = 0;
	while(TRUE_) {
	    io___15.ciunit = isrc;
	    i__1 = s_rsfe(&io___15);
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = do_fio(&c__1, record, (ftnlen)120);
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L100;
	    }
	    if (s_cmp(record, "c     ## ", (ftnlen)9, (ftnlen)9) == 0) {
		next = 10;
		getword_(record, module, &next, (ftnlen)120, (ftnlen)20);
		lowcase_(module, (ftnlen)20);
		upcase_(module, (ftnlen)1);
		if (s_cmp(module, "Subroutine", (ftnlen)20, (ftnlen)10) == 0 
			|| s_cmp(module, "Function", (ftnlen)20, (ftnlen)8) ==
			 0 || s_cmp(module, "Program", (ftnlen)20, (ftnlen)7) 
			== 0) {
		    ++nunit;
		    getword_(record, routine_ref(0, nunit), &next, (ftnlen)
			    120, (ftnlen)120);
		    upcase_(routine_ref(0, nunit), (ftnlen)120);
		    leng = trimtext_(routine_ref(0, nunit), (ftnlen)120);
/* Writing concatenation */
		    i__2[0] = leng, a__1[0] = routine_ref(0, nunit);
		    i__2[1] = 1, a__1[1] = " ";
		    i__2[2] = 20, a__1[2] = module;
		    s_cat(routine_ref(0, nunit), a__1, i__2, &c__3, (ftnlen)
			    120);
		    io___21.ciunit = isrc;
		    i__1 = s_rsfe(&io___21);
		    if (i__1 != 0) {
			goto L100;
		    }
		    i__1 = e_rsfe();
		    if (i__1 != 0) {
			goto L100;
		    }
		    k = 0;
		    done = FALSE_;
		    while(! done) {
			io___24.ciunit = isrc;
			i__1 = s_rsfe(&io___24);
			if (i__1 != 0) {
			    goto L100;
			}
			i__1 = do_fio(&c__1, record, (ftnlen)120);
			if (i__1 != 0) {
			    goto L100;
			}
			i__1 = e_rsfe();
			if (i__1 != 0) {
			    goto L100;
			}
			leng = trimtext_(record, (ftnlen)120);
			if (leng < 7) {
			    done = TRUE_;
			} else if (*(unsigned char *)record == ' ') {
			    done = TRUE_;
			} else {
			    ++k;
			    s_copy(info_ref(0, k, nunit), record + 6, (ftnlen)
				    120, leng - 6);
			}
		    }
		    nline[nunit - 1] = k;
		}
	    }
	}
L100:
	cl__1.cerr = 0;
	cl__1.cunit = isrc;
	cl__1.csta = 0;
	f_clos(&cl__1);
	sort7_(&nunit, routine, link, (ftnlen)120);
	idoc = freeunit_();
	s_copy(docfile, "routines.doc", (ftnlen)120, (ftnlen)12);
	version_(docfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = idoc;
	o__1.ofnmlen = 120;
	o__1.ofnm = docfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	i__1 = nunit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s_copy(string, routine_ref(0, i__), (ftnlen)120, (ftnlen)120);
	    leng = trimtext_(string, (ftnlen)120);
	    if (wiki) {
		io___31.ciunit = idoc;
		s_wsfe(&io___31);
		do_fio(&c__1, string, leng);
		e_wsfe();
	    } else {
		io___32.ciunit = idoc;
		s_wsfe(&io___32);
		do_fio(&c__1, string, leng);
		e_wsfe();
	    }
	    last = 0;
	    j = link[i__ - 1];
	    i__3 = nline[j - 1];
	    for (k = 1; k <= i__3; ++k) {
		s_copy(string, info_ref(0, k, j), (ftnlen)120, (ftnlen)120);
		if (wiki) {
		    if (k == 1) {
			leng = trimtext_(string, (ftnlen)120);
			s_copy(field, string, leng, leng);
			last = leng;
		    } else {
			++last;
			*(unsigned char *)&field[last - 1] = ' ';
			leng = trimtext_(string, (ftnlen)120);
			i__4 = last;
			s_copy(field + i__4, string, last + leng - i__4, leng)
				;
			last += leng;
		    }
		} else {
		    s_copy(string, info_ref(0, k, j), (ftnlen)120, (ftnlen)
			    120);
		    leng = trimtext_(string, (ftnlen)120);
		    io___36.ciunit = idoc;
		    s_wsfe(&io___36);
		    do_fio(&c__1, string, leng);
		    e_wsfe();
		}
	    }
	    if (wiki && last != 0) {
		io___37.ciunit = idoc;
		s_wsfe(&io___37);
		do_fio(&c__1, field, last);
		e_wsfe();
	    }
	    if (nline[j - 1] != 0) {
		io___38.ciunit = idoc;
		s_wsfe(&io___38);
		e_wsfe();
	    }
	}
	cl__1.cerr = 0;
	cl__1.cunit = idoc;
	cl__1.csta = 0;
	f_clos(&cl__1);
	io___39.ciunit = iounit_1.iout;
	s_wsfe(&io___39);
	do_fio(&c__1, docfile, trimtext_(docfile, (ftnlen)120));
	e_wsfe();
    }

/*     get a list of the calls made by each source code routine */

    if (mode == 2) {
	nunit = 0;
	while(TRUE_) {
	    io___40.ciunit = isrc;
	    i__1 = s_rsfe(&io___40);
	    if (i__1 != 0) {
		goto L180;
	    }
	    i__1 = do_fio(&c__1, record, (ftnlen)120);
	    if (i__1 != 0) {
		goto L180;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L180;
	    }
	    upcase_(record, (ftnlen)120);
	    if (*(unsigned char *)record != 'C') {
		next = 1;
		getword_(record, module, &next, (ftnlen)120, (ftnlen)20);
		if (s_cmp(module, "SUBROUTINE", (ftnlen)20, (ftnlen)10) == 0 
			|| s_cmp(module, "FUNCTION", (ftnlen)20, (ftnlen)8) ==
			 0 || s_cmp(module, "PROGRAM", (ftnlen)20, (ftnlen)7) 
			== 0) {
		    ++nunit;
		    getword_(record, routine_ref(0, nunit), &next, (ftnlen)
			    120, (ftnlen)120);
		    nline[nunit - 1] = 0;
		} else {
		    next = i_indx(record, " CALL ", (ftnlen)120, (ftnlen)6);
		    if (next != 0) {
			next += 6;
			getword_(record, keyword, &next, (ftnlen)120, (ftnlen)
				20);
			++nline[nunit - 1];
			s_copy(info_ref(0, nline[nunit - 1], nunit), keyword, 
				(ftnlen)120, (ftnlen)20);
		    } else {
			for (i__ = 1; i__ <= 71; ++i__) {
			    leng = trimtext_(fname_ref(0, i__), (ftnlen)20);
/* Writing concatenation */
			    i__5[0] = leng, a__2[0] = fname_ref(0, i__);
			    i__5[1] = 1, a__2[1] = "(";
			    s_cat(fname1, a__2, i__5, &c__2, (ftnlen)20);
/* Writing concatenation */
			    i__5[0] = leng, a__2[0] = fname_ref(0, i__);
			    i__5[1] = 2, a__2[1] = " (";
			    s_cat(fname2, a__2, i__5, &c__2, (ftnlen)20);
			    if (i_indx(record, fname1, (ftnlen)120, leng + 1) 
				    != 0 || i_indx(record, fname2, (ftnlen)
				    120, leng + 2) != 0) {
				++nline[nunit - 1];
				s_copy(info_ref(0, nline[nunit - 1], nunit), 
					fname_ref(0, i__), (ftnlen)120, (
					ftnlen)20);
			    }
			}
		    }
		}
	    }
	}
L180:
	cl__1.cerr = 0;
	cl__1.cunit = isrc;
	cl__1.csta = 0;
	f_clos(&cl__1);
	sort7_(&nunit, routine, link, (ftnlen)120);
	idoc = freeunit_();
	s_copy(docfile, "calls.doc", (ftnlen)120, (ftnlen)9);
	version_(docfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = idoc;
	o__1.ofnmlen = 120;
	o__1.ofnm = docfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	i__1 = nunit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s_copy(string, routine_ref(0, i__), (ftnlen)120, (ftnlen)120);
	    leng = trimtext_(string, (ftnlen)120);
	    j = link[i__ - 1];
	    sort10_(&nline[j - 1], info_ref(0, 1, j), (ftnlen)120);
	    if (wiki) {
		s_copy(field, string, (ftnlen)2048, leng);
		i__3 = nline[j - 1];
		for (k = 1; k <= i__3; ++k) {
		    leng = trimtext_(info_ref(0, k, j), (ftnlen)120);
		    last = trimtext_(field, (ftnlen)2048);
/* Writing concatenation */
		    i__2[0] = last, a__1[0] = field;
		    i__2[1] = 4, a__1[1] = "    ";
		    i__2[2] = leng, a__1[2] = info_ref(0, k, j);
		    s_cat(field, a__1, i__2, &c__3, (ftnlen)2048);
		}
		leng = trimtext_(field, (ftnlen)2048);
		io___44.ciunit = idoc;
		s_wsfe(&io___44);
		do_fio(&c__1, field, leng);
		e_wsfe();
	    } else {
		io___45.ciunit = idoc;
		s_wsfe(&io___45);
		do_fio(&c__1, string, leng);
		e_wsfe();
		i__3 = nline[j - 1];
		for (k = 1; k <= i__3; ++k) {
		    s_copy(string, info_ref(0, k, j), (ftnlen)120, (ftnlen)
			    120);
		    leng = trimtext_(string, (ftnlen)120);
		    io___46.ciunit = idoc;
		    s_wsfe(&io___46);
		    do_fio(&c__1, string, leng);
		    e_wsfe();
		}
	    }
	}
	cl__1.cerr = 0;
	cl__1.cunit = idoc;
	cl__1.csta = 0;
	f_clos(&cl__1);
	io___47.ciunit = iounit_1.iout;
	s_wsfe(&io___47);
	do_fio(&c__1, docfile, trimtext_(docfile, (ftnlen)120));
	e_wsfe();
    }

/*     get a list of common block contents from source code */

    if (mode == 3) {
	nunit = 0;
	while(TRUE_) {
	    io___48.ciunit = isrc;
	    i__1 = s_rsfe(&io___48);
	    if (i__1 != 0) {
		goto L260;
	    }
	    i__1 = do_fio(&c__1, record, (ftnlen)120);
	    if (i__1 != 0) {
		goto L260;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L260;
	    }
	    if (s_cmp(record, "c     ## ", (ftnlen)9, (ftnlen)9) == 0) {
		next = i_indx(record, ".i  --", (ftnlen)120, (ftnlen)6);
		if (next != 0) {
		    ++nunit;
		    leng = trimtext_(record, (ftnlen)120);
		    upcase_(record + 10, next - 11);
		    s_copy(string, record + 10, (ftnlen)120, next - 11);
		    start = 20;
		    if (wiki) {
			start = trimtext_(string, (ftnlen)120) + 5;
		    }
		    i__1 = next + 7;
		    s_copy(string + (start - 1), record + i__1, 120 - (start 
			    - 1), leng - 4 - i__1);
		    s_copy(routine_ref(0, nunit), string, (ftnlen)120, (
			    ftnlen)120);
		    io___50.ciunit = isrc;
		    i__1 = s_rsfe(&io___50);
		    if (i__1 != 0) {
			goto L260;
		    }
		    i__1 = e_rsfe();
		    if (i__1 != 0) {
			goto L260;
		    }
		    k = 0;
		    done = FALSE_;
		    while(! done) {
			io___51.ciunit = isrc;
			i__1 = s_rsfe(&io___51);
			if (i__1 != 0) {
			    goto L260;
			}
			i__1 = do_fio(&c__1, record, (ftnlen)120);
			if (i__1 != 0) {
			    goto L260;
			}
			i__1 = e_rsfe();
			if (i__1 != 0) {
			    goto L260;
			}
			leng = trimtext_(record, (ftnlen)120);
			if (*(unsigned char *)record == ' ') {
			    done = TRUE_;
			} else if (leng >= 7) {
			    ++k;
			    next = 7;
			    getword_(record, string, &next, (ftnlen)120, (
				    ftnlen)120);
			    s_copy(record, record + (next - 1), (ftnlen)120, 
				    leng - (next - 1));
			    next = nexttext_(record, (ftnlen)120);
			    leng = trimtext_(record, (ftnlen)120);
			    start = 20;
			    if (wiki) {
				start = trimtext_(string, (ftnlen)120) + 5;
			    }
			    s_copy(string + (start - 1), record + (next - 1), 
				    120 - (start - 1), leng - (next - 1));
			    s_copy(info_ref(0, k, nunit), string, (ftnlen)120,
				     (ftnlen)120);
			}
		    }
		    nline[nunit - 1] = k;
		}
	    }
	}
L260:
	cl__1.cerr = 0;
	cl__1.cunit = isrc;
	cl__1.csta = 0;
	f_clos(&cl__1);
	sort7_(&nunit, routine, link, (ftnlen)120);
	idoc = freeunit_();
	s_copy(docfile, "common.doc", (ftnlen)120, (ftnlen)10);
	version_(docfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = idoc;
	o__1.ofnmlen = 120;
	o__1.ofnm = docfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	i__1 = nunit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s_copy(string, routine_ref(0, i__), (ftnlen)120, (ftnlen)120);
	    leng = trimtext_(string, (ftnlen)120);
	    if (wiki) {
		io___52.ciunit = idoc;
		s_wsfe(&io___52);
		do_fio(&c__1, string, leng);
		e_wsfe();
	    } else {
		io___53.ciunit = idoc;
		s_wsfe(&io___53);
		do_fio(&c__1, string, leng);
		e_wsfe();
	    }
	    j = link[i__ - 1];
	    if (wiki) {
		i__3 = nline[j - 1];
		for (k = 1; k <= i__3; ++k) {
		    s_copy(string, info_ref(0, k, j), (ftnlen)120, (ftnlen)
			    120);
		    leng = trimtext_(string, (ftnlen)120);
		    io___54.ciunit = idoc;
		    s_wsfe(&io___54);
		    do_fio(&c__1, string, leng);
		    e_wsfe();
		}
	    } else {
		i__3 = nline[j - 1];
		for (k = 1; k <= i__3; ++k) {
		    s_copy(string, info_ref(0, k, j), (ftnlen)120, (ftnlen)
			    120);
		    leng = trimtext_(string, (ftnlen)120);
		    io___55.ciunit = idoc;
		    s_wsfe(&io___55);
		    do_fio(&c__1, string, leng);
		    e_wsfe();
		}
	    }
	    if (nline[j - 1] != 0) {
		io___56.ciunit = idoc;
		s_wsfe(&io___56);
		e_wsfe();
	    }
	}
	cl__1.cerr = 0;
	cl__1.cunit = idoc;
	cl__1.csta = 0;
	f_clos(&cl__1);
	io___57.ciunit = iounit_1.iout;
	s_wsfe(&io___57);
	do_fio(&c__1, docfile, trimtext_(docfile, (ftnlen)120));
	e_wsfe();
    }

/*     get the keyword values from the source code listing */

    if (mode == 4) {
	nkey = 0;
	while(TRUE_) {
	    io___59.ciunit = isrc;
	    i__1 = s_rsfe(&io___59);
	    if (i__1 != 0) {
		goto L340;
	    }
	    i__1 = do_fio(&c__1, record, (ftnlen)120);
	    if (i__1 != 0) {
		goto L340;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L340;
	    }
	    next = i_indx(record, "if (keyword(", (ftnlen)120, (ftnlen)12);
	    if (next != 0) {
		next = i_indx(record, ".eq.", (ftnlen)120, (ftnlen)4);
		if (next != 0) {
		    next = i_indx(record, "'", (ftnlen)120, (ftnlen)1);
		    if (next != 0) {
			++next;
			getword_(record, keyword, &next, (ftnlen)120, (ftnlen)
				20);
			upcase_(keyword, (ftnlen)20);
			++nkey;
			s_copy(key_ref(0, nkey), keyword, (ftnlen)20, (ftnlen)
				20);
		    }
		}
	    }
	}
L340:
	cl__1.cerr = 0;
	cl__1.cunit = isrc;
	cl__1.csta = 0;
	f_clos(&cl__1);
	sort6_(&nkey, key, (ftnlen)20);
	s_copy(keylast, "                    ", (ftnlen)20, (ftnlen)20);
	idoc = freeunit_();
	s_copy(docfile, "keyword.doc", (ftnlen)120, (ftnlen)11);
	version_(docfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = idoc;
	o__1.ofnmlen = 120;
	o__1.ofnm = docfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	i__1 = nkey;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s_copy(keyword, key_ref(0, i__), (ftnlen)20, (ftnlen)20);
	    leng = trimtext_(keyword, (ftnlen)20);
	    if (s_cmp(keyword, keylast, (ftnlen)20, (ftnlen)20) != 0) {
		io___62.ciunit = idoc;
		s_wsfe(&io___62);
		do_fio(&c__1, keyword, leng);
		e_wsfe();
		s_copy(keylast, keyword, (ftnlen)20, (ftnlen)20);
	    }
	}
	cl__1.cerr = 0;
	cl__1.cunit = idoc;
	cl__1.csta = 0;
	f_clos(&cl__1);
	io___63.ciunit = iounit_1.iout;
	s_wsfe(&io___63);
	do_fio(&c__1, docfile, trimtext_(docfile, (ftnlen)120));
	e_wsfe();
    }

/*     get the included files from the source code listing */

    if (mode == 5) {
	nkey = 0;
	while(TRUE_) {
	    io___64.ciunit = isrc;
	    i__1 = s_rsfe(&io___64);
	    if (i__1 != 0) {
		goto L380;
	    }
	    i__1 = do_fio(&c__1, record, (ftnlen)120);
	    if (i__1 != 0) {
		goto L380;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L380;
	    }
	    next = 1;
	    getword_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	    if (s_cmp(keyword, "include", (ftnlen)20, (ftnlen)7) == 0) {
		gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
		s_copy(keyword, keyword + 1, (ftnlen)20, trimtext_(keyword, (
			ftnlen)20) - 2);
		++nkey;
		s_copy(key_ref(0, nkey), keyword, (ftnlen)20, (ftnlen)20);
	    }
	}
L380:
	cl__1.cerr = 0;
	cl__1.cunit = isrc;
	cl__1.csta = 0;
	f_clos(&cl__1);
	sort6_(&nkey, key, (ftnlen)20);
	s_copy(keylast, "                    ", (ftnlen)20, (ftnlen)20);
	leng = i_indx(srcfile, ".", (ftnlen)120, (ftnlen)1);
/* Writing concatenation */
	i__5[0] = leng - 1, a__2[0] = srcfile;
	i__5[1] = 3, a__2[1] = ".o:";
	s_cat(field, a__2, i__5, &c__2, (ftnlen)2048);
	i__1 = nkey;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s_copy(keyword, key_ref(0, i__), (ftnlen)20, (ftnlen)20);
	    leng = trimtext_(keyword, (ftnlen)20);
	    if (s_cmp(keyword, keylast, (ftnlen)20, (ftnlen)20) != 0) {
		last = trimtext_(field, (ftnlen)2048);
/* Writing concatenation */
		i__2[0] = last, a__1[0] = field;
		i__2[1] = 1, a__1[1] = " ";
		i__2[2] = leng, a__1[2] = keyword;
		s_cat(field, a__1, i__2, &c__3, (ftnlen)2048);
		s_copy(keylast, keyword, (ftnlen)20, (ftnlen)20);
	    }
	}
	io___65.ciunit = iounit_1.iout;
	s_wsfe(&io___65);
	e_wsfe();
	leng = trimtext_(field, (ftnlen)2048);
	io___66.ciunit = iounit_1.iout;
	s_wsfe(&io___66);
	do_fio(&c__1, field, leng);
	e_wsfe();
    }

/*     get a force field parameter file and write a listing */

    if (mode == 6) {
	getprm_();
	idoc = freeunit_();
	s_copy(docfile, "parameter.doc", (ftnlen)120, (ftnlen)13);
	version_(docfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = idoc;
	o__1.ofnmlen = 120;
	o__1.ofnm = docfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	prtprm_(&idoc);
	cl__1.cerr = 0;
	cl__1.cunit = idoc;
	cl__1.csta = 0;
	f_clos(&cl__1);
	io___67.ciunit = iounit_1.iout;
	s_wsfe(&io___67);
	do_fio(&c__1, docfile, trimtext_(docfile, (ftnlen)120));
	e_wsfe();
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef routine_ref
#undef fname_ref
#undef info_ref
#undef key_ref


/* Main program alias */ int document_ () { MAIN__ (); return 0; }
