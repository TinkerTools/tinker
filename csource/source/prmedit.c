/* prmedit.f -- translated by f2c (version 20050501).
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
    doublereal angunit, stbnunit, aaunit, opbunit, opdunit, cang, qang, pang, 
	    sang, copb, qopb, popb, sopb, copd, qopd, popd, sopd;
    char angtyp[600000], opbtyp[8];
} angpot_;

#define angpot_1 angpot_

struct {
    doublereal cbnd, qbnd, bndunit;
    char bndtyp[8];
} bndpot_;

#define bndpot_1 bndpot_

struct {
    integer nprm;
    char prmline[3000000];
} params_;

#define params_1 params_

struct {
    doublereal cury, qury, ureyunit;
} urypot_;

#define urypot_1 urypot_

struct {
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__0 = 0;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  2004  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  program prmedit  --  edit and renumber parameter files  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "prmedit" reformats an existing parameter file, and revises */
/*     type and class numbers based on the "atom" parameter ordering */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 The Parameter Editing Facility can Provi"
	    "de :\002,//,4x,\002(1) Format Individual Parameter Records\002,/"
	    ",4x,\002(2) Reorder Individual Parameter Records\002,/,4x,\002(3"
	    ") Renumber the Atom Types, and Reorder\002,/,4x,\002(4) Renumber"
	    " the Atom Classes, and Reorder\002,/,4x,\002(5) Renumber Types a"
	    "nd Classes, and Reorder\002/,4x,\002(6) Sort and Format Multipol"
	    "e Parameters\002)";
    static char fmt_30[] = "(/,\002 Enter the Number of the Desired Choice :"
	    "  \002,$)";
    static char fmt_40[] = "(i10)";
    static char fmt_60[] = "(/,\002 Reformated Parameter File Written To: "
	    " \002,a)";
    static char fmt_70[] = "(/,\002 Renumbered Parameter File Written To: "
	    " \002,a)";
    static char fmt_80[] = "(/,\002 Sorted Multipole Values Written To:  "
	    "\002,a)";

    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern integer freeunit_(void);
    extern /* Subroutine */ int prmorder_(integer *, logical *, logical *), 
	    polesort_(integer *);
    extern integer trimtext_(char *, ftnlen);
    static integer mode, iprm;
    extern /* Subroutine */ int final_(void);
    static logical exist, query;
    extern /* Subroutine */ int getprm_(void);
    static logical dotype;
    static char string[120];
    static logical doclass;
    extern /* Subroutine */ int initial_(void);
    static char prmfile[120];
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), prmform_(
	    integer *), version_(char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___5 = { 1, string, 1, 0, 120, 1 };
    static cilist io___6 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___8 = { 1, 0, 1, fmt_40, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_80, 0 };




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




/*     read and store the original force field parameter file */

    initial_();
    getprm_();

/*     get the desired type of parameter file modification */

    mode = 0;
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___5);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&mode, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L10;
	}
	query = FALSE_;
    }
L10:
    if (query) {
	io___6.ciunit = iounit_1.iout;
	s_wsfe(&io___6);
	e_wsfe();
	while(mode < 1 || mode > 6) {
	    mode = 0;
	    io___7.ciunit = iounit_1.iout;
	    s_wsfe(&io___7);
	    e_wsfe();
	    io___8.ciunit = iounit_1.input;
	    i__1 = s_rsfe(&io___8);
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_fio(&c__1, (char *)&mode, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L50;
	    }
L50:
	    ;
	}
    }

/*     set the renumbering operations to be performed */

    dotype = FALSE_;
    doclass = FALSE_;
    if (mode == 3) {
	dotype = TRUE_;
    }
    if (mode == 4) {
	doclass = TRUE_;
    }
    if (mode == 5) {
	dotype = TRUE_;
	doclass = TRUE_;
    }

/*     format records in the original parameter file */

    if (mode == 1) {
	iprm = freeunit_();
	s_copy(prmfile, "parameter.prm", (ftnlen)120, (ftnlen)13);
	version_(prmfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = iprm;
	o__1.ofnmlen = 120;
	o__1.ofnm = prmfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	prmform_(&iprm);
	io___13.ciunit = iounit_1.iout;
	s_wsfe(&io___13);
	do_fio(&c__1, prmfile, trimtext_(prmfile, (ftnlen)120));
	e_wsfe();
	cl__1.cerr = 0;
	cl__1.cunit = iprm;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     reorder and renumber the original parameter file */

    if (mode >= 2 && mode <= 5) {
	iprm = freeunit_();
	s_copy(prmfile, "parameter.prm", (ftnlen)120, (ftnlen)13);
	version_(prmfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = iprm;
	o__1.ofnmlen = 120;
	o__1.ofnm = prmfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	prmorder_(&iprm, &dotype, &doclass);
	io___14.ciunit = iounit_1.iout;
	s_wsfe(&io___14);
	do_fio(&c__1, prmfile, trimtext_(prmfile, (ftnlen)120));
	e_wsfe();
	cl__1.cerr = 0;
	cl__1.cunit = iprm;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     sort the atomic multipole parameters by atom type */

    if (mode == 6) {
	iprm = freeunit_();
	s_copy(prmfile, "multipole.prm", (ftnlen)120, (ftnlen)13);
	version_(prmfile, "new", (ftnlen)120, (ftnlen)3);
	o__1.oerr = 0;
	o__1.ounit = iprm;
	o__1.ofnmlen = 120;
	o__1.ofnm = prmfile;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	polesort_(&iprm);
	io___15.ciunit = iounit_1.iout;
	s_wsfe(&io___15);
	do_fio(&c__1, prmfile, trimtext_(prmfile, (ftnlen)120));
	e_wsfe();
	cl__1.cerr = 0;
	cl__1.cunit = iprm;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    final_();
    return 0;
} /* MAIN__ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine prmform  --  reformat individual parameters  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "prmform" formats each individual parameter record to conform */
/*     to a consistent text layout */


/* Subroutine */ int prmform_(integer *iprm)
{
    /* Format strings */
    static char fmt_20[] = "(\002atom\002,6x,2i5,4x,a3,3x,a26,1x,i5,f10.3,i5)"
	    ;
    static char fmt_30[] = "(\002atom\002,6x,i5,4x,a3,3x,a26,1x,i5,f10.3,i5)";
    static char fmt_40[] = "(a)";
    static char fmt_60[] = "(\002vdw\002,7x,i5,10x,2f11.4)";
    static char fmt_70[] = "(\002vdw\002,7x,i5,10x,2f11.4,f11.3)";
    static char fmt_90[] = "(\002vdw14\002,5x,i5,10x,2f11.4)";
    static char fmt_110[] = "(\002vdwpr\002,5x,2i5,5x,2f11.4)";
    static char fmt_130[] = "(\002hbond\002,5x,2i5,5x,2f11.4)";
    static char fmt_150[] = "(\002bond\002,6x,2i5,5x,f11.2,f11.4)";
    static char fmt_160[] = "(\002bond\002,6x,2i5,5x,f11.3,f11.4)";
    static char fmt_180[] = "(\002bond5\002,5x,2i5,5x,f11.2,f11.4)";
    static char fmt_190[] = "(\002bond5\002,5x,2i5,5x,f11.3,f11.4)";
    static char fmt_210[] = "(\002bond4\002,5x,2i5,5x,f11.2,f11.4)";
    static char fmt_220[] = "(\002bond4\002,5x,2i5,5x,f11.3,f11.4)";
    static char fmt_240[] = "(\002bond3\002,5x,2i5,5x,f11.2,f11.4)";
    static char fmt_250[] = "(\002bond3\002,5x,2i5,5x,f11.3,f11.4)";
    static char fmt_270[] = "(\002electneg\002,2x,3i5,11x,f11.4)";
    static char fmt_290[] = "(\002angle\002,5x,3i5,f11.2,f11.2)";
    static char fmt_300[] = "(\002angle\002,5x,3i5,f11.3,f11.2)";
    static char fmt_310[] = "(\002angle\002,5x,3i5,f11.2,3f11.2)";
    static char fmt_320[] = "(\002angle\002,5x,3i5,f11.3,3f11.2)";
    static char fmt_340[] = "(\002angle5\002,4x,3i5,f11.2,f11.2)";
    static char fmt_350[] = "(\002angle5\002,4x,3i5,f11.3,f11.2)";
    static char fmt_360[] = "(\002angle5\002,4x,3i5,f11.2,3f11.2)";
    static char fmt_370[] = "(\002angle5\002,4x,3i5,f11.3,3f11.2)";
    static char fmt_390[] = "(\002angle4\002,4x,3i5,f11.2,f11.2)";
    static char fmt_400[] = "(\002angle4\002,4x,3i5,f11.3,f11.2)";
    static char fmt_410[] = "(\002angle4\002,4x,3i5,f11.2,3f11.2)";
    static char fmt_420[] = "(\002angle4\002,4x,3i5,f11.3,3f11.2)";
    static char fmt_440[] = "(\002angle3\002,4x,3i5,f11.2,f11.2)";
    static char fmt_450[] = "(\002angle3\002,4x,3i5,f11.3,f11.2)";
    static char fmt_460[] = "(\002angle3\002,4x,3i5,f11.2,3f11.2)";
    static char fmt_470[] = "(\002angle3\002,4x,3i5,f11.3,3f11.2)";
    static char fmt_490[] = "(\002anglef\002,4x,3i5,f11.2,f11.2,f11.1)";
    static char fmt_500[] = "(\002anglef\002,4x,3i5,f11.3,f11.2,f11.1)";
    static char fmt_520[] = "(\002strbnd\002,4x,3i5,2f11.2)";
    static char fmt_530[] = "(\002strbnd\002,4x,3i5,2f11.3)";
    static char fmt_550[] = "(\002ureybrad\002,2x,3i5,f11.2,f11.4)";
    static char fmt_560[] = "(\002ureybrad\002,2x,3i5,f11.3,f11.4)";
    static char fmt_580[] = "(\002angang\002,4x,i5,10x,3f11.2)";
    static char fmt_590[] = "(\002angang\002,4x,i5,10x,3f11.3)";
    static char fmt_610[] = "(\002opbend\002,4x,4i5,6x,f11.2)";
    static char fmt_620[] = "(\002opbend\002,4x,4i5,6x,f11.3)";
    static char fmt_640[] = "(\002opdist\002,4x,4i5,6x,f11.2)";
    static char fmt_650[] = "(\002opdist\002,4x,4i5,6x,f11.3)";
    static char fmt_670[] = "(\002improper\002,2x,4i5,6x,2f11.2)";
    static char fmt_690[] = "(\002imptors\002,3x,4i5,6x,6(f11.3,f7.1,i3))";
    static char fmt_710[] = "(\002torsion\002,3x,4i5,3x,f8.3,f4.1,i2,f8.3,f6"
	    ".1,i2,f8.3,f4.1,i2)";
    static char fmt_720[] = "(\002torsion\002,3x,4i5,6x,2(f11.3,f7.1,i3))";
    static char fmt_730[] = "(\002torsion\002,3x,4i5,3x,6(f8.3,f6.1,i2))";
    static char fmt_750[] = "(\002torsion5\002,2x,4i5,3x,f8.3,f4.1,i2,f8.3,f"
	    "6.1,i2,f8.3,f4.1,i2)";
    static char fmt_760[] = "(\002torsion5\002,2x,4i5,6x,2(f11.3,f7.1,i3))";
    static char fmt_770[] = "(\002torsion5\002,2x,4i5,3x,6(f8.3,f6.1,i2))";
    static char fmt_790[] = "(\002torsion4\002,2x,4i5,3x,f8.3,f4.1,i2,f8.3,f"
	    "6.1,i2,f8.3,f4.1,i2)";
    static char fmt_800[] = "(\002torsion4\002,2x,4i5,6x,2(f11.3,f7.1,i3))";
    static char fmt_810[] = "(\002torsion4\002,2x,4i5,3x,6(f8.3,f6.1,i2))";
    static char fmt_830[] = "(\002pitors\002,4x,2i5,5x,f11.2)";
    static char fmt_850[] = "(\002strtors\002,3x,4i5,1x,3f11.3)";
    static char fmt_870[] = "(\002tortors\002,3x,5i5,5x,2i5)";
    static char fmt_890[] = "(f8.1,1x,f8.1,1x,f11.5)";
    static char fmt_910[] = "(\002charge\002,4x,i5,10x,f11.4)";
    static char fmt_930[] = "(\002dipole\002,4x,2i5,5x,f11.4,f11.3)";
    static char fmt_950[] = "(\002dipole5\002,3x,2i5,5x,f11.4,f11.3)";
    static char fmt_970[] = "(\002dipole4\002,3x,2i5,5x,f11.4,f11.3)";
    static char fmt_990[] = "(\002dipole3\002,3x,2i5,5x,f11.4,f11.3)";
    static char fmt_1020[] = "(\002multipole\002,1x,3i5,11x,f11.5)";
    static char fmt_1030[] = "(\002multipole\002,1x,4i5,6x,f11.5)";
    static char fmt_1050[] = "(36x,3f11.5)";
    static char fmt_1070[] = "(36x,f11.5)";
    static char fmt_1090[] = "(36x,2f11.5)";
    static char fmt_1110[] = "(36x,3f11.5)";
    static char fmt_1130[] = "(\002polarize\002,2x,i5,5x,2f11.3,3x,8i5)";
    static char fmt_1150[] = "(\002piatom\002,4x,i5,10x,f11.1,2f11.3)";
    static char fmt_1170[] = "(\002pibond\002,4x,2i5,5x,f11.3,f11.4)";
    static char fmt_1190[] = "(\002pibond5\002,3x,2i5,5x,f11.3,f11.4)";
    static char fmt_1210[] = "(\002pibond4\002,3x,2i5,5x,f11.3,f11.4)";
    static char fmt_1220[] = "(\002metal\002,5x,i5,a)";
    static char fmt_1240[] = "(\002biotype\002,3x,i5,4x,a3,5x,a30,i5)";
    static char fmt_1250[] = "()";
    static char fmt_1260[] = "(a)";

    /* System generated locals */
    address a__1[4];
    integer i__1, i__2[4], i__3;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j;
    extern /* Subroutine */ int getstring_(char *, char *, integer *, ftnlen, 
	    ftnlen);
    static doublereal bd, fc;
    static integer ia, ib, ic, id, ie;
    static doublereal an;
    static integer ig[8];
    static doublereal dl;
    static integer kg;
    static doublereal dk, cg, ep, rd, ds, dp;
    static integer ft[6];
    static doublereal vd, tf, pl, el;
    static integer kt;
    static doublereal pr, ps, pt, iz;
    static integer nx, ny;
    static doublereal rp, ss, ts, vt[6], tx, ty, st[6], aa1, ba1, ba2, aa2, 
	    aa3, an1, an2, an3, bt1, bt2, bt3, pl1, pl2, pl3;
    static integer lig, atn;
    static doublereal rdn, thl, pol;
    static char sym[3];
    static integer nxy;
    static char note[24];
    static doublereal wght;
    static integer next;
    extern /* Subroutine */ int sort_(integer *, integer *);
    static char blank[30], record[120];
    static integer length;
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int getnumb_(char *, integer *, integer *, ftnlen)
	    , getword_(char *, char *, integer *, ftnlen, ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___30 = { 1, string, 1, 0, 120, 1 };
    static cilist io___31 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_40, 0 };
    static icilist io___37 = { 1, string, 1, 0, 120, 1 };
    static cilist io___38 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_70, 0 };
    static icilist io___40 = { 1, string, 1, 0, 120, 1 };
    static cilist io___41 = { 0, 0, 0, fmt_90, 0 };
    static icilist io___42 = { 1, string, 1, 0, 120, 1 };
    static cilist io___43 = { 0, 0, 0, fmt_110, 0 };
    static icilist io___44 = { 1, string, 1, 0, 120, 1 };
    static cilist io___45 = { 0, 0, 0, fmt_130, 0 };
    static icilist io___48 = { 1, string, 1, 0, 120, 1 };
    static cilist io___49 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_160, 0 };
    static icilist io___51 = { 1, string, 1, 0, 120, 1 };
    static cilist io___52 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_190, 0 };
    static icilist io___54 = { 1, string, 1, 0, 120, 1 };
    static cilist io___55 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___56 = { 0, 0, 0, fmt_220, 0 };
    static icilist io___57 = { 1, string, 1, 0, 120, 1 };
    static cilist io___58 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___59 = { 0, 0, 0, fmt_250, 0 };
    static icilist io___62 = { 1, string, 1, 0, 120, 1 };
    static cilist io___63 = { 0, 0, 0, fmt_270, 0 };
    static icilist io___67 = { 1, string, 1, 0, 120, 1 };
    static cilist io___68 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___69 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___71 = { 0, 0, 0, fmt_320, 0 };
    static icilist io___72 = { 1, string, 1, 0, 120, 1 };
    static cilist io___73 = { 0, 0, 0, fmt_340, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_350, 0 };
    static cilist io___75 = { 0, 0, 0, fmt_360, 0 };
    static cilist io___76 = { 0, 0, 0, fmt_370, 0 };
    static icilist io___77 = { 1, string, 1, 0, 120, 1 };
    static cilist io___78 = { 0, 0, 0, fmt_390, 0 };
    static cilist io___79 = { 0, 0, 0, fmt_400, 0 };
    static cilist io___80 = { 0, 0, 0, fmt_410, 0 };
    static cilist io___81 = { 0, 0, 0, fmt_420, 0 };
    static icilist io___82 = { 1, string, 1, 0, 120, 1 };
    static cilist io___83 = { 0, 0, 0, fmt_440, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_450, 0 };
    static cilist io___85 = { 0, 0, 0, fmt_460, 0 };
    static cilist io___86 = { 0, 0, 0, fmt_470, 0 };
    static icilist io___89 = { 1, string, 1, 0, 120, 1 };
    static cilist io___90 = { 0, 0, 0, fmt_490, 0 };
    static cilist io___91 = { 0, 0, 0, fmt_500, 0 };
    static icilist io___94 = { 1, string, 1, 0, 120, 1 };
    static cilist io___95 = { 0, 0, 0, fmt_520, 0 };
    static cilist io___96 = { 0, 0, 0, fmt_530, 0 };
    static icilist io___98 = { 1, string, 1, 0, 120, 1 };
    static cilist io___99 = { 0, 0, 0, fmt_550, 0 };
    static cilist io___100 = { 0, 0, 0, fmt_560, 0 };
    static icilist io___104 = { 1, string, 1, 0, 120, 1 };
    static cilist io___105 = { 0, 0, 0, fmt_580, 0 };
    static cilist io___106 = { 0, 0, 0, fmt_590, 0 };
    static icilist io___108 = { 1, string, 1, 0, 120, 1 };
    static cilist io___109 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___110 = { 0, 0, 0, fmt_620, 0 };
    static icilist io___111 = { 1, string, 1, 0, 120, 1 };
    static cilist io___112 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___113 = { 0, 0, 0, fmt_650, 0 };
    static icilist io___116 = { 1, string, 1, 0, 120, 1 };
    static cilist io___117 = { 0, 0, 0, fmt_670, 0 };
    static icilist io___122 = { 1, string, 1, 0, 120, 1 };
    static cilist io___124 = { 0, 0, 0, fmt_690, 0 };
    static icilist io___125 = { 1, string, 1, 0, 120, 1 };
    static cilist io___126 = { 0, 0, 0, fmt_710, 0 };
    static cilist io___127 = { 0, 0, 0, fmt_720, 0 };
    static cilist io___128 = { 0, 0, 0, fmt_730, 0 };
    static icilist io___129 = { 1, string, 1, 0, 120, 1 };
    static cilist io___130 = { 0, 0, 0, fmt_750, 0 };
    static cilist io___131 = { 0, 0, 0, fmt_760, 0 };
    static cilist io___132 = { 0, 0, 0, fmt_770, 0 };
    static icilist io___133 = { 1, string, 1, 0, 120, 1 };
    static cilist io___134 = { 0, 0, 0, fmt_790, 0 };
    static cilist io___135 = { 0, 0, 0, fmt_800, 0 };
    static cilist io___136 = { 0, 0, 0, fmt_810, 0 };
    static icilist io___138 = { 1, string, 1, 0, 120, 1 };
    static cilist io___139 = { 0, 0, 0, fmt_830, 0 };
    static icilist io___143 = { 1, string, 1, 0, 120, 1 };
    static cilist io___144 = { 0, 0, 0, fmt_850, 0 };
    static icilist io___148 = { 1, string, 1, 0, 120, 1 };
    static cilist io___149 = { 0, 0, 0, fmt_870, 0 };
    static icilist io___151 = { 1, record, 1, 0, 120, 1 };
    static cilist io___155 = { 0, 0, 0, fmt_890, 0 };
    static icilist io___157 = { 1, string, 1, 0, 120, 1 };
    static cilist io___158 = { 0, 0, 0, fmt_910, 0 };
    static icilist io___161 = { 1, string, 1, 0, 120, 1 };
    static cilist io___162 = { 0, 0, 0, fmt_930, 0 };
    static icilist io___163 = { 1, string, 1, 0, 120, 1 };
    static cilist io___164 = { 0, 0, 0, fmt_950, 0 };
    static icilist io___165 = { 1, string, 1, 0, 120, 1 };
    static cilist io___166 = { 0, 0, 0, fmt_970, 0 };
    static icilist io___167 = { 1, string, 1, 0, 120, 1 };
    static cilist io___168 = { 0, 0, 0, fmt_990, 0 };
    static icilist io___170 = { 1, string, 1, 0, 120, 1 };
    static icilist io___171 = { 1, string, 1, 0, 120, 1 };
    static cilist io___172 = { 0, 0, 0, fmt_1020, 0 };
    static cilist io___173 = { 0, 0, 0, fmt_1030, 0 };
    static icilist io___174 = { 1, record, 1, 0, 120, 1 };
    static cilist io___178 = { 0, 0, 0, fmt_1050, 0 };
    static icilist io___179 = { 1, record, 1, 0, 120, 1 };
    static cilist io___180 = { 0, 0, 0, fmt_1070, 0 };
    static icilist io___181 = { 1, record, 1, 0, 120, 1 };
    static cilist io___182 = { 0, 0, 0, fmt_1090, 0 };
    static icilist io___183 = { 1, record, 1, 0, 120, 1 };
    static cilist io___184 = { 0, 0, 0, fmt_1110, 0 };
    static icilist io___188 = { 1, string, 1, 0, 120, 1 };
    static cilist io___190 = { 0, 0, 0, fmt_1130, 0 };
    static icilist io___194 = { 1, string, 1, 0, 120, 1 };
    static cilist io___195 = { 0, 0, 0, fmt_1150, 0 };
    static icilist io___198 = { 1, string, 1, 0, 120, 1 };
    static cilist io___199 = { 0, 0, 0, fmt_1170, 0 };
    static icilist io___200 = { 1, string, 1, 0, 120, 1 };
    static cilist io___201 = { 0, 0, 0, fmt_1190, 0 };
    static icilist io___202 = { 1, string, 1, 0, 120, 1 };
    static cilist io___203 = { 0, 0, 0, fmt_1210, 0 };
    static cilist io___204 = { 0, 0, 0, fmt_1220, 0 };
    static icilist io___205 = { 1, string, 1, 0, 120, 1 };
    static icilist io___206 = { 1, string, 1, 0, 120, 1 };
    static cilist io___207 = { 0, 0, 0, fmt_1240, 0 };
    static cilist io___208 = { 0, 0, 0, fmt_1250, 0 };
    static cilist io___209 = { 0, 0, 0, fmt_1260, 0 };



#define prmline_ref(a_0,a_1) &params_1.prmline[(a_1)*120 + a_0 - 120]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  angpot.i  --  specifics of angle bend functional forms  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     angunit    convert angle bending energy to kcal/mole */
/*     stbnunit   convert stretch-bend energy to kcal/mole */
/*     aaunit     convert angle-angle energy to kcal/mole */
/*     opbunit    convert out-of-plane bend energy to kcal/mole */
/*     opdunit    convert out-of-plane distance energy to kcal/mole */
/*     cang       cubic coefficient in angle bending potential */
/*     qang       quartic coefficient in angle bending potential */
/*     pang       quintic coefficient in angle bending potential */
/*     sang       sextic coefficient in angle bending potential */
/*     copb       cubic coefficient in out-of-plane bend potential */
/*     qopb       quartic coefficient in out-of-plane bend potential */
/*     popb       quintic coefficient in out-of-plane bend potential */
/*     sopb       sextic coefficient in out-of-plane bend potential */
/*     copd       cubic coefficient in out-of-plane distance potential */
/*     qopd       quartic coefficient in out-of-plane distance potential */
/*     popd       quintic coefficient in out-of-plane distance potential */
/*     sopd       sextic coefficient in out-of-plane distance potential */
/*     angtyp     type of angle bending function for each bond angle */
/*     opbtyp     type of out-of-plane bend potential energy function */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  bndpot.i  --  specifics of bond stretch functional forms  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     cbnd      cubic coefficient in bond stretch potential */
/*     qbnd      quartic coefficient in bond stretch potential */
/*     bndunit   convert bond stretch energy to kcal/mole */
/*     bndtyp    type of bond stretch potential energy function */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  math.i  --  mathematical and geometrical constants  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     radian   conversion factor from radians to degrees */
/*     pi       numerical value of the geometric constant */
/*     sqrtpi   numerical value of the square root of Pi */
/*     logten   numerical value of the natural log of ten */
/*     sqrttwo  numerical value of the square root of two */
/*     twosix   numerical value of the sixth root of two */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  params.i  --  contents of force field parameter file  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     nprm      number of nonblank lines in the parameter file */
/*     prmline   contents of each individual parameter file line */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  urypot.i  --  specifics of Urey-Bradley functional form  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     cury       cubic coefficient in Urey-Bradley potential */
/*     qury       quartic coefficient in Urey-Bradley potential */
/*     ureyunit   convert Urey-Bradley energy to kcal/mole */




/*     reformat and print the various parameters */

    i__ = 0;
    s_copy(blank, "                              ", (ftnlen)30, (ftnlen)30);
    while(i__ < params_1.nprm) {
	++i__;
	s_copy(record, prmline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	length = trimtext_(record, (ftnlen)120);
	next = 1;
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	if (s_cmp(keyword, "ATOM ", (ftnlen)5, (ftnlen)5) == 0) {
	    ia = -1;
	    ib = -1;
	    s_copy(sym, "   ", (ftnlen)3, (ftnlen)3);
	    s_copy(note, "                        ", (ftnlen)24, (ftnlen)24);
	    atn = 0;
	    wght = 0.;
	    lig = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    gettext_(record, sym, &next, (ftnlen)120, (ftnlen)3);
	    getstring_(record, note, &next, (ftnlen)120, (ftnlen)24);
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___30);
	    if (i__1 != 0) {
		goto L10;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&atn, (ftnlen)sizeof(integer))
		    ;
	    if (i__1 != 0) {
		goto L10;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&wght, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L10;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&lig, (ftnlen)sizeof(integer))
		    ;
	    if (i__1 != 0) {
		goto L10;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L10;
	    }
L10:
	    length = trimtext_(note, (ftnlen)24);
/* Writing concatenation */
	    i__2[0] = 1, a__1[0] = "\"";
	    i__2[1] = length, a__1[1] = note;
	    i__2[2] = 1, a__1[2] = "\"";
	    i__2[3] = 30, a__1[3] = blank;
	    s_cat(string, a__1, i__2, &c__4, (ftnlen)120);
	    if (ib >= 0) {
		io___31.ciunit = *iprm;
		s_wsfe(&io___31);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, sym, (ftnlen)3);
		do_fio(&c__1, string, (ftnlen)26);
		do_fio(&c__1, (char *)&atn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&wght, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&lig, (ftnlen)sizeof(integer));
		e_wsfe();
	    } else if (ia >= 0) {
		io___32.ciunit = *iprm;
		s_wsfe(&io___32);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, sym, (ftnlen)3);
		do_fio(&c__1, string, (ftnlen)26);
		do_fio(&c__1, (char *)&atn, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&wght, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&lig, (ftnlen)sizeof(integer));
		e_wsfe();
	    } else {
		io___33.ciunit = *iprm;
		s_wsfe(&io___33);
		do_fio(&c__1, record, length);
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "VDW ", (ftnlen)4, (ftnlen)4) == 0) {
	    ia = 0;
	    rd = 0.;
	    ep = 0.;
	    rdn = 0.;
	    i__1 = s_rsli(&io___37);
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&rd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ep, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&rdn, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L50;
	    }
L50:
	    if (rdn == 0.) {
		io___38.ciunit = *iprm;
		s_wsfe(&io___38);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&rd, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&ep, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___39.ciunit = *iprm;
		s_wsfe(&io___39);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&rd, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&ep, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&rdn, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "VDW14 ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    rd = 0.;
	    ep = 0.;
	    i__1 = s_rsli(&io___40);
	    if (i__1 != 0) {
		goto L80;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L80;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&rd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L80;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ep, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L80;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L80;
	    }
L80:
	    io___41.ciunit = *iprm;
	    s_wsfe(&io___41);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&rd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ep, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "VDWPR ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    rd = 0.;
	    ep = 0.;
	    i__1 = s_rsli(&io___42);
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&rd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ep, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L100;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L100;
	    }
L100:
	    io___43.ciunit = *iprm;
	    s_wsfe(&io___43);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&rd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ep, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "HBOND ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    rd = 0.;
	    ep = 0.;
	    i__1 = s_rsli(&io___44);
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&rd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ep, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L120;
	    }
L120:
	    io___45.ciunit = *iprm;
	    s_wsfe(&io___45);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&rd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ep, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "BOND ", (ftnlen)5, (ftnlen)5) == 0) {
	    ia = 0;
	    ib = 0;
	    fc = 0.;
	    bd = 0.;
	    i__1 = s_rsli(&io___48);
	    if (i__1 != 0) {
		goto L140;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L140;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L140;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L140;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&bd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L140;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L140;
	    }
L140:
	    if (bndpot_1.bndunit < 10.) {
		io___49.ciunit = *iprm;
		s_wsfe(&io___49);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bd, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___50.ciunit = *iprm;
		s_wsfe(&io___50);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bd, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "BOND5 ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    fc = 0.;
	    bd = 0.;
	    i__1 = s_rsli(&io___51);
	    if (i__1 != 0) {
		goto L170;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L170;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L170;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L170;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&bd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L170;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L170;
	    }
L170:
	    if (bndpot_1.bndunit < 10.) {
		io___52.ciunit = *iprm;
		s_wsfe(&io___52);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bd, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___53.ciunit = *iprm;
		s_wsfe(&io___53);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bd, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "BOND4 ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    fc = 0.;
	    bd = 0.;
	    i__1 = s_rsli(&io___54);
	    if (i__1 != 0) {
		goto L200;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L200;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L200;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L200;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&bd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L200;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L200;
	    }
L200:
	    if (bndpot_1.bndunit < 10.) {
		io___55.ciunit = *iprm;
		s_wsfe(&io___55);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bd, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___56.ciunit = *iprm;
		s_wsfe(&io___56);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bd, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "BOND3 ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    fc = 0.;
	    bd = 0.;
	    i__1 = s_rsli(&io___57);
	    if (i__1 != 0) {
		goto L230;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L230;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L230;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L230;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&bd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L230;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L230;
	    }
L230:
	    if (bndpot_1.bndunit < 10.) {
		io___58.ciunit = *iprm;
		s_wsfe(&io___58);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bd, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___59.ciunit = *iprm;
		s_wsfe(&io___59);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bd, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "ELECTNEG ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    dl = 0.;
	    i__1 = s_rsli(&io___62);
	    if (i__1 != 0) {
		goto L260;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L260;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L260;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L260;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&dl, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L260;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L260;
	    }
L260:
	    io___63.ciunit = *iprm;
	    s_wsfe(&io___63);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dl, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "ANGLE ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    an1 = 0.;
	    an2 = 0.;
	    an3 = 0.;
	    i__1 = s_rsli(&io___67);
	    if (i__1 != 0) {
		goto L280;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L280;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L280;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L280;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L280;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L280;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L280;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L280;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L280;
	    }
L280:
	    if (an2 == 0. && an3 == 0.) {
		if (angpot_1.angunit < .0030461741978670856) {
		    io___68.ciunit = *iprm;
		    s_wsfe(&io___68);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else {
		    io___69.ciunit = *iprm;
		    s_wsfe(&io___69);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    } else {
		if (angpot_1.angunit < .0030461741978670856) {
		    io___70.ciunit = *iprm;
		    s_wsfe(&io___70);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an3, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else {
		    io___71.ciunit = *iprm;
		    s_wsfe(&io___71);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an3, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	} else if (s_cmp(keyword, "ANGLE5 ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    an1 = 0.;
	    an2 = 0.;
	    an3 = 0.;
	    i__1 = s_rsli(&io___72);
	    if (i__1 != 0) {
		goto L330;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L330;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L330;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L330;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L330;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L330;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L330;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L330;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L330;
	    }
L330:
	    if (an2 == 0. && an3 == 0.) {
		if (angpot_1.angunit < .0030461741978670856) {
		    io___73.ciunit = *iprm;
		    s_wsfe(&io___73);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else {
		    io___74.ciunit = *iprm;
		    s_wsfe(&io___74);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    } else {
		if (angpot_1.angunit < .0030461741978670856) {
		    io___75.ciunit = *iprm;
		    s_wsfe(&io___75);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an3, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else {
		    io___76.ciunit = *iprm;
		    s_wsfe(&io___76);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an3, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	} else if (s_cmp(keyword, "ANGLE4 ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    an1 = 0.;
	    an2 = 0.;
	    an3 = 0.;
	    i__1 = s_rsli(&io___77);
	    if (i__1 != 0) {
		goto L380;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L380;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L380;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L380;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L380;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L380;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L380;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L380;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L380;
	    }
L380:
	    if (an2 == 0. && an3 == 0.) {
		if (angpot_1.angunit < .0030461741978670856) {
		    io___78.ciunit = *iprm;
		    s_wsfe(&io___78);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else {
		    io___79.ciunit = *iprm;
		    s_wsfe(&io___79);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    } else {
		if (angpot_1.angunit < .0030461741978670856) {
		    io___80.ciunit = *iprm;
		    s_wsfe(&io___80);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an3, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else {
		    io___81.ciunit = *iprm;
		    s_wsfe(&io___81);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an3, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	} else if (s_cmp(keyword, "ANGLE3 ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    an1 = 0.;
	    an2 = 0.;
	    an3 = 0.;
	    i__1 = s_rsli(&io___82);
	    if (i__1 != 0) {
		goto L430;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L430;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L430;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L430;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L430;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L430;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L430;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L430;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L430;
	    }
L430:
	    if (an2 == 0. && an3 == 0.) {
		if (angpot_1.angunit < .0030461741978670856) {
		    io___83.ciunit = *iprm;
		    s_wsfe(&io___83);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else {
		    io___84.ciunit = *iprm;
		    s_wsfe(&io___84);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    } else {
		if (angpot_1.angunit < .0030461741978670856) {
		    io___85.ciunit = *iprm;
		    s_wsfe(&io___85);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an3, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		} else {
		    io___86.ciunit = *iprm;
		    s_wsfe(&io___86);
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&an3, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    }
	} else if (s_cmp(keyword, "ANGLEF ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    an = 0.;
	    pr = 0.;
	    i__1 = s_rsli(&io___89);
	    if (i__1 != 0) {
		goto L480;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L480;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L480;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L480;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L480;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L480;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pr, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L480;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L480;
	    }
L480:
	    if (angpot_1.angunit < .0030461741978670856) {
		io___90.ciunit = *iprm;
		s_wsfe(&io___90);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&an, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&pr, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___91.ciunit = *iprm;
		s_wsfe(&io___91);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&an, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&pr, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "STRBND ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    ba1 = 0.;
	    ba2 = 0.;
	    i__1 = s_rsli(&io___94);
	    if (i__1 != 0) {
		goto L510;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L510;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L510;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L510;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ba1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L510;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ba2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L510;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L510;
	    }
L510:
	    if (angpot_1.stbnunit < .17453292519943295) {
		io___95.ciunit = *iprm;
		s_wsfe(&io___95);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ba1, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&ba2, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___96.ciunit = *iprm;
		s_wsfe(&io___96);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ba1, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&ba2, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "UREYBRAD ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    ds = 0.;
	    i__1 = s_rsli(&io___98);
	    if (i__1 != 0) {
		goto L540;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L540;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L540;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L540;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L540;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ds, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L540;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L540;
	    }
L540:
	    if (urypot_1.ureyunit < 10.) {
		io___99.ciunit = *iprm;
		s_wsfe(&io___99);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&ds, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___100.ciunit = *iprm;
		s_wsfe(&io___100);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&ds, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "ANGANG ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    aa1 = 0.;
	    aa2 = 0.;
	    aa3 = 0.;
	    i__1 = s_rsli(&io___104);
	    if (i__1 != 0) {
		goto L570;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L570;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&aa1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L570;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&aa2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L570;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&aa3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L570;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L570;
	    }
L570:
	    if (abs(angpot_1.aaunit) < .0030461741978670856) {
		io___105.ciunit = *iprm;
		s_wsfe(&io___105);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&aa1, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&aa2, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&aa3, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___106.ciunit = *iprm;
		s_wsfe(&io___106);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&aa1, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&aa2, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&aa3, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "OPBEND ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    fc = 0.;
	    i__1 = s_rsli(&io___108);
	    if (i__1 != 0) {
		goto L600;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L600;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L600;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L600;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L600;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L600;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L600;
	    }
L600:
	    if (angpot_1.opbunit < .0030461741978670856) {
		io___109.ciunit = *iprm;
		s_wsfe(&io___109);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___110.ciunit = *iprm;
		s_wsfe(&io___110);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "OPDIST ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    fc = 0.;
	    i__1 = s_rsli(&io___111);
	    if (i__1 != 0) {
		goto L630;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L630;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L630;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L630;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L630;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L630;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L630;
	    }
L630:
	    if (angpot_1.opdunit < 10.) {
		io___112.ciunit = *iprm;
		s_wsfe(&io___112);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___113.ciunit = *iprm;
		s_wsfe(&io___113);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&fc, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "IMPROPER ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    dk = 0.;
	    vd = 0.;
	    i__1 = s_rsli(&io___116);
	    if (i__1 != 0) {
		goto L660;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L660;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L660;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L660;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L660;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&dk, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L660;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&vd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L660;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L660;
	    }
L660:
	    io___117.ciunit = *iprm;
	    s_wsfe(&io___117);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dk, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&vd, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "IMPTORS ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    for (j = 1; j <= 6; ++j) {
		vt[j - 1] = 0.;
		st[j - 1] = 0.;
		ft[j - 1] = 0;
	    }
	    i__1 = s_rsli(&io___122);
	    if (i__1 != 0) {
		goto L680;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L680;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L680;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L680;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L680;
	    }
	    for (j = 1; j <= 6; ++j) {
		i__1 = do_lio(&c__5, &c__1, (char *)&vt[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L680;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&st[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L680;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&ft[j - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L680;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L680;
	    }
L680:
	    kt = 0;
	    for (j = 1; j <= 6; ++j) {
		if (ft[j - 1] != 0) {
		    kt = j;
		}
	    }
	    io___124.ciunit = *iprm;
	    s_wsfe(&io___124);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    i__1 = kt;
	    for (j = 1; j <= i__1; ++j) {
		do_fio(&c__1, (char *)&vt[j - 1], (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&st[j - 1], (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&ft[j - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	} else if (s_cmp(keyword, "TORSION ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    for (j = 1; j <= 6; ++j) {
		vt[j - 1] = 0.;
		st[j - 1] = 0.;
		ft[j - 1] = 0;
	    }
	    i__1 = s_rsli(&io___125);
	    if (i__1 != 0) {
		goto L700;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L700;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L700;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L700;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L700;
	    }
	    for (j = 1; j <= 6; ++j) {
		i__1 = do_lio(&c__5, &c__1, (char *)&vt[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L700;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&st[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L700;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&ft[j - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L700;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L700;
	    }
L700:
	    kt = 0;
	    for (j = 1; j <= 6; ++j) {
		if (ft[j - 1] != 0) {
		    kt = j;
		}
	    }
	    if (kt == 3 && st[0] == 0. && st[1] == 180. && st[2] == 0.) {
		io___126.ciunit = *iprm;
		s_wsfe(&io___126);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		i__1 = kt;
		for (j = 1; j <= i__1; ++j) {
		    do_fio(&c__1, (char *)&vt[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&st[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&ft[j - 1], (ftnlen)sizeof(integer))
			    ;
		}
		e_wsfe();
	    } else if (kt <= 2) {
		io___127.ciunit = *iprm;
		s_wsfe(&io___127);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		i__1 = kt;
		for (j = 1; j <= i__1; ++j) {
		    do_fio(&c__1, (char *)&vt[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&st[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&ft[j - 1], (ftnlen)sizeof(integer))
			    ;
		}
		e_wsfe();
	    } else {
		io___128.ciunit = *iprm;
		s_wsfe(&io___128);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		i__1 = kt;
		for (j = 1; j <= i__1; ++j) {
		    do_fio(&c__1, (char *)&vt[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&st[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&ft[j - 1], (ftnlen)sizeof(integer))
			    ;
		}
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "TORSION5 ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    for (j = 1; j <= 6; ++j) {
		vt[j - 1] = 0.;
		st[j - 1] = 0.;
		ft[j - 1] = 0;
	    }
	    i__1 = s_rsli(&io___129);
	    if (i__1 != 0) {
		goto L740;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L740;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L740;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L740;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L740;
	    }
	    for (j = 1; j <= 6; ++j) {
		i__1 = do_lio(&c__5, &c__1, (char *)&vt[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L740;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&st[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L740;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&ft[j - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L740;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L740;
	    }
L740:
	    kt = 0;
	    for (j = 1; j <= 6; ++j) {
		if (ft[j - 1] != 0) {
		    kt = j;
		}
	    }
	    if (kt == 3 && st[0] == 0. && st[1] == 180. && st[2] == 0.) {
		io___130.ciunit = *iprm;
		s_wsfe(&io___130);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		i__1 = kt;
		for (j = 1; j <= i__1; ++j) {
		    do_fio(&c__1, (char *)&vt[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&st[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&ft[j - 1], (ftnlen)sizeof(integer))
			    ;
		}
		e_wsfe();
	    } else if (kt <= 2) {
		io___131.ciunit = *iprm;
		s_wsfe(&io___131);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		i__1 = kt;
		for (j = 1; j <= i__1; ++j) {
		    do_fio(&c__1, (char *)&vt[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&st[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&ft[j - 1], (ftnlen)sizeof(integer))
			    ;
		}
		e_wsfe();
	    } else {
		io___132.ciunit = *iprm;
		s_wsfe(&io___132);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		i__1 = kt;
		for (j = 1; j <= i__1; ++j) {
		    do_fio(&c__1, (char *)&vt[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&st[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&ft[j - 1], (ftnlen)sizeof(integer))
			    ;
		}
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "TORSION4 ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    for (j = 1; j <= 6; ++j) {
		vt[j - 1] = 0.;
		st[j - 1] = 0.;
		ft[j - 1] = 0;
	    }
	    i__1 = s_rsli(&io___133);
	    if (i__1 != 0) {
		goto L780;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L780;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L780;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L780;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L780;
	    }
	    for (j = 1; j <= 6; ++j) {
		i__1 = do_lio(&c__5, &c__1, (char *)&vt[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L780;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&st[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L780;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&ft[j - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L780;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L780;
	    }
L780:
	    kt = 0;
	    for (j = 1; j <= 6; ++j) {
		if (ft[j - 1] != 0) {
		    kt = j;
		}
	    }
	    if (kt == 3 && st[0] == 0. && st[1] == 180. && st[2] == 0.) {
		io___134.ciunit = *iprm;
		s_wsfe(&io___134);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		i__1 = kt;
		for (j = 1; j <= i__1; ++j) {
		    do_fio(&c__1, (char *)&vt[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&st[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&ft[j - 1], (ftnlen)sizeof(integer))
			    ;
		}
		e_wsfe();
	    } else if (kt <= 2) {
		io___135.ciunit = *iprm;
		s_wsfe(&io___135);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		i__1 = kt;
		for (j = 1; j <= i__1; ++j) {
		    do_fio(&c__1, (char *)&vt[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&st[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&ft[j - 1], (ftnlen)sizeof(integer))
			    ;
		}
		e_wsfe();
	    } else {
		io___136.ciunit = *iprm;
		s_wsfe(&io___136);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		i__1 = kt;
		for (j = 1; j <= i__1; ++j) {
		    do_fio(&c__1, (char *)&vt[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&st[j - 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&ft[j - 1], (ftnlen)sizeof(integer))
			    ;
		}
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "PITORS ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    pt = 0.;
	    i__1 = s_rsli(&io___138);
	    if (i__1 != 0) {
		goto L820;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L820;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L820;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pt, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L820;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L820;
	    }
L820:
	    io___139.ciunit = *iprm;
	    s_wsfe(&io___139);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&pt, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "STRTORS ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    bt1 = 0.;
	    bt2 = 0.;
	    bt3 = 0.;
	    i__1 = s_rsli(&io___143);
	    if (i__1 != 0) {
		goto L840;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L840;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L840;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L840;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L840;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&bt1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L840;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&bt2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L840;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&bt3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L840;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L840;
	    }
L840:
	    io___144.ciunit = *iprm;
	    s_wsfe(&io___144);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&bt1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&bt2, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&bt3, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "TORTORS ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    ie = 0;
	    nx = 0;
	    ny = 0;
	    i__1 = s_rsli(&io___148);
	    if (i__1 != 0) {
		goto L860;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L860;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L860;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L860;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L860;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ie, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L860;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&nx, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L860;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ny, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L860;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L860;
	    }
L860:
	    io___149.ciunit = *iprm;
	    s_wsfe(&io___149);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ie, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nx, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ny, (ftnlen)sizeof(integer));
	    e_wsfe();
	    nxy = nx * ny;
	    i__1 = nxy;
	    for (j = 1; j <= i__1; ++j) {
		++i__;
		s_copy(record, prmline_ref(0, i__), (ftnlen)120, (ftnlen)120);
		i__3 = s_rsli(&io___151);
		if (i__3 != 0) {
		    goto L880;
		}
		i__3 = do_lio(&c__5, &c__1, (char *)&tx, (ftnlen)sizeof(
			doublereal));
		if (i__3 != 0) {
		    goto L880;
		}
		i__3 = do_lio(&c__5, &c__1, (char *)&ty, (ftnlen)sizeof(
			doublereal));
		if (i__3 != 0) {
		    goto L880;
		}
		i__3 = do_lio(&c__5, &c__1, (char *)&tf, (ftnlen)sizeof(
			doublereal));
		if (i__3 != 0) {
		    goto L880;
		}
		i__3 = e_rsli();
		if (i__3 != 0) {
		    goto L880;
		}
L880:
		io___155.ciunit = *iprm;
		s_wsfe(&io___155);
		do_fio(&c__1, (char *)&tx, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&ty, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&tf, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "CHARGE ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    cg = 0.;
	    i__1 = s_rsli(&io___157);
	    if (i__1 != 0) {
		goto L900;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L900;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&cg, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L900;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L900;
	    }
L900:
	    io___158.ciunit = *iprm;
	    s_wsfe(&io___158);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&cg, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "DIPOLE ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    dp = 0.;
	    ps = .5;
	    i__1 = s_rsli(&io___161);
	    if (i__1 != 0) {
		goto L920;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L920;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L920;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&dp, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L920;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ps, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L920;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L920;
	    }
L920:
	    io___162.ciunit = *iprm;
	    s_wsfe(&io___162);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dp, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ps, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "DIPOLE5 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    dp = 0.;
	    ps = .5;
	    i__1 = s_rsli(&io___163);
	    if (i__1 != 0) {
		goto L940;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L940;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L940;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&dp, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L940;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ps, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L940;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L940;
	    }
L940:
	    io___164.ciunit = *iprm;
	    s_wsfe(&io___164);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dp, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ps, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "DIPOLE4 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    dp = 0.;
	    ps = .5;
	    i__1 = s_rsli(&io___165);
	    if (i__1 != 0) {
		goto L960;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L960;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L960;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&dp, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L960;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ps, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L960;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L960;
	    }
L960:
	    io___166.ciunit = *iprm;
	    s_wsfe(&io___166);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dp, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ps, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "DIPOLE3 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    dp = 0.;
	    ps = .5;
	    i__1 = s_rsli(&io___167);
	    if (i__1 != 0) {
		goto L980;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L980;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L980;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&dp, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L980;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ps, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L980;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L980;
	    }
L980:
	    io___168.ciunit = *iprm;
	    s_wsfe(&io___168);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dp, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ps, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "MULTIPOLE ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    pl = 0.;
	    i__1 = s_rsli(&io___170);
	    if (i__1 != 0) {
		goto L1000;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1000;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1000;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1000;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1000;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1000;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1000;
	    }
	    goto L1010;
L1000:
	    id = 0;
	    i__1 = s_rsli(&io___171);
	    if (i__1 != 0) {
		goto L1010;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1010;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1010;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1010;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1010;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1010;
	    }
L1010:
	    if (id == 0) {
		io___172.ciunit = *iprm;
		s_wsfe(&io___172);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&pl, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		io___173.ciunit = *iprm;
		s_wsfe(&io___173);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&pl, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    ++i__;
	    s_copy(record, prmline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	    i__1 = s_rsli(&io___174);
	    if (i__1 != 0) {
		goto L1040;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1040;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1040;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1040;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1040;
	    }
L1040:
	    io___178.ciunit = *iprm;
	    s_wsfe(&io___178);
	    do_fio(&c__1, (char *)&pl1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&pl2, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&pl3, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    ++i__;
	    s_copy(record, prmline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	    i__1 = s_rsli(&io___179);
	    if (i__1 != 0) {
		goto L1060;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1060;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1060;
	    }
L1060:
	    io___180.ciunit = *iprm;
	    s_wsfe(&io___180);
	    do_fio(&c__1, (char *)&pl1, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    ++i__;
	    s_copy(record, prmline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	    i__1 = s_rsli(&io___181);
	    if (i__1 != 0) {
		goto L1080;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1080;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1080;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1080;
	    }
L1080:
	    io___182.ciunit = *iprm;
	    s_wsfe(&io___182);
	    do_fio(&c__1, (char *)&pl1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&pl2, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    ++i__;
	    s_copy(record, prmline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	    i__1 = s_rsli(&io___183);
	    if (i__1 != 0) {
		goto L1100;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1100;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1100;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1100;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1100;
	    }
L1100:
	    io___184.ciunit = *iprm;
	    s_wsfe(&io___184);
	    do_fio(&c__1, (char *)&pl1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&pl2, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&pl3, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "POLARIZE ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    pol = 0.;
	    thl = 0.;
	    for (j = 1; j <= 8; ++j) {
		ig[j - 1] = 0;
	    }
	    i__1 = s_rsli(&io___188);
	    if (i__1 != 0) {
		goto L1120;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1120;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pol, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1120;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&thl, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1120;
	    }
	    for (j = 1; j <= 8; ++j) {
		i__1 = do_lio(&c__3, &c__1, (char *)&ig[j - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L1120;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1120;
	    }
L1120:
	    kg = 0;
	    for (j = 1; j <= 8; ++j) {
		if (ig[j - 1] != 0) {
		    kg = j;
		}
	    }
	    sort_(&kg, ig);
	    io___190.ciunit = *iprm;
	    s_wsfe(&io___190);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&pol, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&thl, (ftnlen)sizeof(doublereal));
	    i__1 = kg;
	    for (j = 1; j <= i__1; ++j) {
		do_fio(&c__1, (char *)&ig[j - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	} else if (s_cmp(keyword, "PIATOM ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    el = 0.;
	    iz = 0.;
	    rp = 0.;
	    i__1 = s_rsli(&io___194);
	    if (i__1 != 0) {
		goto L1140;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1140;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&el, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1140;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&iz, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1140;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&rp, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1140;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1140;
	    }
L1140:
	    io___195.ciunit = *iprm;
	    s_wsfe(&io___195);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&el, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&iz, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&rp, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "PIBOND ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ss = 0.;
	    ts = 0.;
	    i__1 = s_rsli(&io___198);
	    if (i__1 != 0) {
		goto L1160;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1160;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1160;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ss, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1160;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ts, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1160;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1160;
	    }
L1160:
	    io___199.ciunit = *iprm;
	    s_wsfe(&io___199);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ss, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ts, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "PIBOND5 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ss = 0.;
	    ts = 0.;
	    i__1 = s_rsli(&io___200);
	    if (i__1 != 0) {
		goto L1180;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1180;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1180;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ss, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1180;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ts, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1180;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1180;
	    }
L1180:
	    io___201.ciunit = *iprm;
	    s_wsfe(&io___201);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ss, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ts, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "PIBOND4 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ss = 0.;
	    ts = 0.;
	    i__1 = s_rsli(&io___202);
	    if (i__1 != 0) {
		goto L1200;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1200;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1200;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ss, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1200;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ts, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L1200;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1200;
	    }
L1200:
	    io___203.ciunit = *iprm;
	    s_wsfe(&io___203);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ss, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ts, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (s_cmp(keyword, "METAL ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    io___204.ciunit = *iprm;
	    s_wsfe(&io___204);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "BIOTYPE ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    s_copy(sym, "   ", (ftnlen)3, (ftnlen)3);
	    s_copy(note, "                        ", (ftnlen)24, (ftnlen)24);
	    i__1 = s_rsli(&io___205);
	    if (i__1 != 0) {
		goto L1230;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1230;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1230;
	    }
	    getword_(record, sym, &next, (ftnlen)120, (ftnlen)3);
	    getstring_(record, note, &next, (ftnlen)120, (ftnlen)24);
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___206);
	    if (i__1 != 0) {
		goto L1230;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L1230;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L1230;
	    }
L1230:
	    length = trimtext_(note, (ftnlen)24);
/* Writing concatenation */
	    i__2[0] = 1, a__1[0] = "\"";
	    i__2[1] = length, a__1[1] = note;
	    i__2[2] = 1, a__1[2] = "\"";
	    i__2[3] = 30, a__1[3] = blank;
	    s_cat(string, a__1, i__2, &c__4, (ftnlen)120);
	    io___207.ciunit = *iprm;
	    s_wsfe(&io___207);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, sym, (ftnlen)3);
	    do_fio(&c__1, string, (ftnlen)30);
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else if (length == 0) {
	    io___208.ciunit = *iprm;
	    s_wsfe(&io___208);
	    e_wsfe();
	} else {
	    io___209.ciunit = *iprm;
	    s_wsfe(&io___209);
	    do_fio(&c__1, record, length);
	    e_wsfe();
	}
    }
    return 0;
} /* prmform_ */

#undef prmline_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine prmorder  --  reorder atom types and classes  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "prmorder" places a list of atom type or class numbers into */
/*     canonical order for potential energy parameter definitions */


/* Subroutine */ int prmorder_(integer *iprm, logical *dotype, logical *
	doclass)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter Starting Number for Atom Types [1]"
	    " :  \002,$)";
    static char fmt_20[] = "(i10)";
    static char fmt_30[] = "(/,\002 Enter Starting Number for Atom Classes ["
	    "1] :  \002,$)";
    static char fmt_40[] = "(i10)";
    static char fmt_50[] = "(\002atom\002,6x,2i5,a)";
    static char fmt_60[] = "(\002vdw\002,7x,i5,a)";
    static char fmt_70[] = "(\002vdw14\002,5x,i5,a)";
    static char fmt_80[] = "(\002vdwpr\002,5x,2i5,a)";
    static char fmt_90[] = "(\002hbond\002,5x,2i5,a)";
    static char fmt_100[] = "(\002bond\002,6x,2i5,a)";
    static char fmt_110[] = "(\002bond5\002,5x,2i5,a)";
    static char fmt_120[] = "(\002bond4\002,5x,2i5,a)";
    static char fmt_130[] = "(\002bond3\002,5x,2i5,a)";
    static char fmt_140[] = "(\002electneg\002,2x,3i5,a)";
    static char fmt_150[] = "(\002angle\002,5x,3i5,a)";
    static char fmt_160[] = "(\002angle5\002,4x,3i5,a)";
    static char fmt_170[] = "(\002angle4\002,4x,3i5,a)";
    static char fmt_180[] = "(\002angle3\002,4x,3i5,a)";
    static char fmt_190[] = "(\002anglef\002,4x,3i5,a)";
    static char fmt_200[] = "(\002strbnd\002,4x,3i5,a)";
    static char fmt_210[] = "(\002ureybrad\002,2x,3i5,a)";
    static char fmt_220[] = "(\002angang\002,4x,i5,a)";
    static char fmt_230[] = "(\002opbend\002,4x,4i5,a)";
    static char fmt_240[] = "(\002opdist\002,4x,4i5,a)";
    static char fmt_250[] = "(\002improper\002,2x,4i5,a)";
    static char fmt_260[] = "(\002imptors\002,3x,4i5,a)";
    static char fmt_270[] = "(\002torsion\002,3x,4i5,a)";
    static char fmt_280[] = "(\002torsion5\002,2x,4i5,a)";
    static char fmt_290[] = "(\002torsion4\002,2x,4i5,a)";
    static char fmt_300[] = "(\002pitors\002,4x,2i5,a)";
    static char fmt_310[] = "(\002strtors\002,3x,4i5,a)";
    static char fmt_320[] = "(\002tortors\002,3x,5i5,a)";
    static char fmt_330[] = "(\002charge\002,4x,i5,a)";
    static char fmt_340[] = "(\002dipole\002,4x,2i5,a)";
    static char fmt_350[] = "(\002dipole5\002,3x,2i5,a)";
    static char fmt_360[] = "(\002dipole4\002,3x,2i5,a)";
    static char fmt_370[] = "(\002dipole3\002,3x,2i5,a)";
    static char fmt_380[] = "(\002multipole\002,1x,3i5,a)";
    static char fmt_390[] = "(\002multipole\002,1x,4i5,a)";
    static char fmt_410[] = "(\002polarize\002,2x,i5,f20.3,f11.3,3x,8i5)";
    static char fmt_420[] = "(\002piatom\002,4x,i5,a)";
    static char fmt_430[] = "(\002pibond\002,4x,2i5,a)";
    static char fmt_440[] = "(\002pibond5\002,3x,2i5,a)";
    static char fmt_450[] = "(\002pibond4\002,3x,2i5,a)";
    static char fmt_460[] = "(\002metal\002,5x,i5,a)";
    static char fmt_480[] = "(\002biotype\002,a,i5)";
    static char fmt_490[] = "()";
    static char fmt_500[] = "(a)";

    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3, i__4, i__5, i__6[2];
    char ch__1[150];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void), s_cmp(char *, char *, ftnlen, 
	    ftnlen), i_sign(integer *, integer *), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j;
    extern /* Subroutine */ int getstring_(char *, char *, integer *, ftnlen, 
	    ftnlen);
    static integer ia, ib, ic, id, kc, ie, ig[8], kg, it, kt;
    static doublereal thl, pol;
    static integer next;
    extern /* Subroutine */ int sort_(integer *, integer *);
    static char blank[30];
    static integer itype[5001], iclass[1001], length;
    static char record[120];
    static integer offset;
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int getnumb_(char *, integer *, integer *, ftnlen)
	    , getword_(char *, char *, integer *, ftnlen, ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen), prmsort_(integer *, integer *, integer *, integer *, 
	    integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___221 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___222 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___224 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___225 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___231 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___232 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___233 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___234 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___235 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___236 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___237 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___238 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___239 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___240 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___241 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___242 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___243 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___244 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___245 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___246 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___247 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___248 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___249 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___250 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___251 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___252 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___253 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___254 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___255 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___256 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___257 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___258 = { 0, 0, 0, fmt_320, 0 };
    static cilist io___259 = { 0, 0, 0, fmt_330, 0 };
    static cilist io___260 = { 0, 0, 0, fmt_340, 0 };
    static cilist io___261 = { 0, 0, 0, fmt_350, 0 };
    static cilist io___262 = { 0, 0, 0, fmt_360, 0 };
    static cilist io___263 = { 0, 0, 0, fmt_370, 0 };
    static cilist io___264 = { 0, 0, 0, fmt_380, 0 };
    static cilist io___265 = { 0, 0, 0, fmt_390, 0 };
    static icilist io___272 = { 1, string, 1, 0, 120, 1 };
    static cilist io___273 = { 0, 0, 0, fmt_410, 0 };
    static cilist io___274 = { 0, 0, 0, fmt_420, 0 };
    static cilist io___275 = { 0, 0, 0, fmt_430, 0 };
    static cilist io___276 = { 0, 0, 0, fmt_440, 0 };
    static cilist io___277 = { 0, 0, 0, fmt_450, 0 };
    static cilist io___278 = { 0, 0, 0, fmt_460, 0 };
    static icilist io___279 = { 1, string, 1, 0, 120, 1 };
    static icilist io___280 = { 1, string, 1, 0, 120, 1 };
    static cilist io___281 = { 0, 0, 0, fmt_480, 0 };
    static cilist io___282 = { 0, 0, 0, fmt_490, 0 };
    static cilist io___283 = { 0, 0, 0, fmt_500, 0 };



#define prmline_ref(a_0,a_1) &params_1.prmline[(a_1)*120 + a_0 - 120]



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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  params.i  --  contents of force field parameter file  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     nprm      number of nonblank lines in the parameter file */
/*     prmline   contents of each individual parameter file line */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  vdwpot.i  --  specifics of van der Waals functional form  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     abuck      value of "A" constant in Buckingham vdw potential */
/*     bbuck      value of "B" constant in Buckingham vdw potential */
/*     cbuck      value of "C" constant in Buckingham vdw potential */
/*     ghal       value of "gamma" in buffered 14-7 vdw potential */
/*     dhal       value of "delta" in buffered 14-7 vdw potential */
/*     v2scale    factor by which 1-2 vdw interactions are scaled */
/*     v3scale    factor by which 1-3 vdw interactions are scaled */
/*     v4scale    factor by which 1-4 vdw interactions are scaled */
/*     v5scale    factor by which 1-5 vdw interactions are scaled */
/*     igauss     coefficients of Gaussian fit to vdw potential */
/*     ngauss     number of Gaussians used in fit to vdw potential */
/*     vdwindex   indexing mode (atom type or class) for vdw parameters */
/*     vdwtyp     type of van der Waals potential energy function */
/*     radtyp     type of parameter (sigma or R-min) for atomic size */
/*     radsiz     atomic size provided as radius or diameter */
/*     radrule    combining rule for atomic size parameters */
/*     epsrule    combining rule for vdw well depth parameters */
/*     gausstyp   type of Gaussian fit to van der Waals potential */




/*     zero out the storage for atom types and classes */

    ia = 0;
    ib = 0;
    ic = 0;
    id = 0;
    ie = 0;
    kt = 0;
    kc = 0;
    for (i__ = 0; i__ <= 5000; ++i__) {
	itype[i__] = 0;
    }
    for (i__ = 0; i__ <= 1000; ++i__) {
	iclass[i__] = 0;
    }
    s_copy(blank, "                              ", (ftnlen)30, (ftnlen)30);

/*     get the starting numbers for atom types and classes */

    if (*dotype) {
	io___221.ciunit = iounit_1.iout;
	s_wsfe(&io___221);
	e_wsfe();
	io___222.ciunit = iounit_1.input;
	s_rsfe(&io___222);
	do_fio(&c__1, (char *)&offset, (ftnlen)sizeof(integer));
	e_rsfe();
	if (offset > 0) {
	    kt = offset - 1;
	}
    }
    if (*doclass) {
	io___224.ciunit = iounit_1.iout;
	s_wsfe(&io___224);
	e_wsfe();
	io___225.ciunit = iounit_1.input;
	s_rsfe(&io___225);
	do_fio(&c__1, (char *)&offset, (ftnlen)sizeof(integer));
	e_rsfe();
	if (offset > 0) {
	    kc = offset - 1;
	}
    }

/*     count and order the atom types and atom classes */

    i__1 = params_1.nprm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(record, prmline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	next = 1;
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "ATOM ", (ftnlen)5, (ftnlen)5) == 0) {
	    it = 0;
	    ic = 0;
	    getnumb_(record, &it, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    if (ic == 0) {
		ic = it;
	    }
	    if (itype[it] == 0) {
		++kt;
		if (*dotype) {
		    itype[it] = kt;
		} else {
		    itype[it] = it;
		}
	    }
	    if (iclass[ic] == 0) {
		++kc;
		if (*doclass) {
		    iclass[ic] = kc;
		} else {
		    iclass[ic] = ic;
		}
	    }
	}
    }

/*     reorder, renumber and print the various parameters */

    i__1 = params_1.nprm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(record, prmline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	length = trimtext_(record, (ftnlen)120);
	next = 1;
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "ATOM ", (ftnlen)5, (ftnlen)5) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    if (ib == 0) {
		ib = ia;
	    }
	    ia = itype[ia];
	    ib = iclass[ib];
	    io___231.ciunit = *iprm;
	    s_wsfe(&io___231);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "VDW ", (ftnlen)4, (ftnlen)4) == 0) {
	    ia = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    if (s_cmp(vdwpot_1.vdwindex, "TYPE", (ftnlen)5, (ftnlen)4) == 0) {
		ia = itype[ia];
	    } else {
		ia = iclass[ia];
	    }
	    io___232.ciunit = *iprm;
	    s_wsfe(&io___232);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "VDW14 ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    if (s_cmp(vdwpot_1.vdwindex, "TYPE", (ftnlen)5, (ftnlen)4) == 0) {
		ia = itype[ia];
	    } else {
		ia = iclass[ia];
	    }
	    io___233.ciunit = *iprm;
	    s_wsfe(&io___233);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "VDWPR ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    if (s_cmp(vdwpot_1.vdwindex, "TYPE", (ftnlen)5, (ftnlen)4) == 0) {
		ia = itype[ia];
		ib = itype[ib];
	    } else {
		ia = iclass[ia];
		ib = iclass[ib];
	    }
	    prmsort_(&c__2, &ia, &ib, &c__0, &c__0, &c__0);
	    io___234.ciunit = *iprm;
	    s_wsfe(&io___234);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "HBOND ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    if (s_cmp(vdwpot_1.vdwindex, "TYPE", (ftnlen)5, (ftnlen)4) == 0) {
		ia = itype[ia];
		ib = itype[ib];
	    } else {
		ia = iclass[ia];
		ib = iclass[ib];
	    }
	    prmsort_(&c__2, &ia, &ib, &c__0, &c__0, &c__0);
	    io___235.ciunit = *iprm;
	    s_wsfe(&io___235);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "BOND ", (ftnlen)5, (ftnlen)5) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    prmsort_(&c__2, &ia, &ib, &c__0, &c__0, &c__0);
	    io___236.ciunit = *iprm;
	    s_wsfe(&io___236);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "BOND5 ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    prmsort_(&c__2, &ia, &ib, &c__0, &c__0, &c__0);
	    io___237.ciunit = *iprm;
	    s_wsfe(&io___237);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "BOND4 ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    prmsort_(&c__2, &ia, &ib, &c__0, &c__0, &c__0);
	    io___238.ciunit = *iprm;
	    s_wsfe(&io___238);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "BOND3 ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    prmsort_(&c__2, &ia, &ib, &c__0, &c__0, &c__0);
	    io___239.ciunit = *iprm;
	    s_wsfe(&io___239);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "ELECTNEG ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    io___240.ciunit = *iprm;
	    s_wsfe(&io___240);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "ANGLE ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    prmsort_(&c__3, &ia, &ib, &ic, &c__0, &c__0);
	    io___241.ciunit = *iprm;
	    s_wsfe(&io___241);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "ANGLE5 ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    prmsort_(&c__3, &ia, &ib, &ic, &c__0, &c__0);
	    io___242.ciunit = *iprm;
	    s_wsfe(&io___242);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "ANGLE4 ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    prmsort_(&c__3, &ia, &ib, &ic, &c__0, &c__0);
	    io___243.ciunit = *iprm;
	    s_wsfe(&io___243);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "ANGLE3 ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    prmsort_(&c__3, &ia, &ib, &ic, &c__0, &c__0);
	    io___244.ciunit = *iprm;
	    s_wsfe(&io___244);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "ANGLEF ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    prmsort_(&c__3, &ia, &ib, &ic, &c__0, &c__0);
	    io___245.ciunit = *iprm;
	    s_wsfe(&io___245);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "STRBND ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    io___246.ciunit = *iprm;
	    s_wsfe(&io___246);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "UREYBRAD ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    prmsort_(&c__3, &ia, &ib, &ic, &c__0, &c__0);
	    io___247.ciunit = *iprm;
	    s_wsfe(&io___247);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "ANGANG ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    ia = iclass[ia];
	    io___248.ciunit = *iprm;
	    s_wsfe(&io___248);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "OPBEND ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    getnumb_(record, &id, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    id = iclass[id];
	    prmsort_(&c__2, &ic, &id, &c__0, &c__0, &c__0);
	    io___249.ciunit = *iprm;
	    s_wsfe(&io___249);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "OPDIST ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    getnumb_(record, &id, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    id = iclass[id];
	    prmsort_(&c__2, &ib, &ic, &c__0, &c__0, &c__0);
	    prmsort_(&c__2, &ib, &id, &c__0, &c__0, &c__0);
	    prmsort_(&c__2, &ic, &id, &c__0, &c__0, &c__0);
	    io___250.ciunit = *iprm;
	    s_wsfe(&io___250);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "IMPROPER ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    getnumb_(record, &id, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    id = iclass[id];
	    io___251.ciunit = *iprm;
	    s_wsfe(&io___251);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "IMPTORS ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    getnumb_(record, &id, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    id = iclass[id];
	    io___252.ciunit = *iprm;
	    s_wsfe(&io___252);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "TORSION ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    getnumb_(record, &id, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    id = iclass[id];
	    prmsort_(&c__4, &ia, &ib, &ic, &id, &c__0);
	    io___253.ciunit = *iprm;
	    s_wsfe(&io___253);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "TORSION5 ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    getnumb_(record, &id, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    id = iclass[id];
	    prmsort_(&c__4, &ia, &ib, &ic, &id, &c__0);
	    io___254.ciunit = *iprm;
	    s_wsfe(&io___254);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "TORSION4 ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    getnumb_(record, &id, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    id = iclass[id];
	    prmsort_(&c__4, &ia, &ib, &ic, &id, &c__0);
	    io___255.ciunit = *iprm;
	    s_wsfe(&io___255);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "PITORS ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    prmsort_(&c__2, &ia, &ib, &c__0, &c__0, &c__0);
	    io___256.ciunit = *iprm;
	    s_wsfe(&io___256);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "STRTORS ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    getnumb_(record, &id, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    id = iclass[id];
	    prmsort_(&c__4, &ia, &ib, &ic, &id, &c__0);
	    io___257.ciunit = *iprm;
	    s_wsfe(&io___257);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "TORTORS ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    ie = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    getnumb_(record, &id, &next, (ftnlen)120);
	    getnumb_(record, &ie, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    ic = iclass[ic];
	    id = iclass[id];
	    ie = iclass[ie];
	    io___258.ciunit = *iprm;
	    s_wsfe(&io___258);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ie, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "CHARGE ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    ia = itype[ia];
	    io___259.ciunit = *iprm;
	    s_wsfe(&io___259);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "DIPOLE ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    ia = itype[ia];
	    ib = itype[ib];
	    io___260.ciunit = *iprm;
	    s_wsfe(&io___260);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "DIPOLE5 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    ia = itype[ia];
	    ib = itype[ib];
	    io___261.ciunit = *iprm;
	    s_wsfe(&io___261);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "DIPOLE4 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    ia = itype[ia];
	    ib = itype[ib];
	    io___262.ciunit = *iprm;
	    s_wsfe(&io___262);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "DIPOLE3 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    ia = itype[ia];
	    ib = itype[ib];
	    io___263.ciunit = *iprm;
	    s_wsfe(&io___263);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "MULTIPOLE ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    getnumb_(record, &id, &next, (ftnlen)120);
	    ia = itype[ia];
	    ib = i_sign(&c__1, &ib) * itype[abs(ib)];
	    ic = i_sign(&c__1, &ic) * itype[abs(ic)];
	    id = i_sign(&c__1, &id) * itype[abs(id)];
	    if (id == 0) {
		io___264.ciunit = *iprm;
		s_wsfe(&io___264);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, record + (next - 1), length - (next - 1));
		e_wsfe();
	    } else {
		io___265.ciunit = *iprm;
		s_wsfe(&io___265);
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, record + (next - 1), length - (next - 1));
		e_wsfe();
	    }
	} else if (s_cmp(keyword, "POLARIZE ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    pol = 0.;
	    thl = 0.;
	    kg = 0;
	    for (j = 1; j <= 8; ++j) {
		ig[j - 1] = 0;
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___272);
	    if (i__2 != 0) {
		goto L400;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L400;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&pol, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L400;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&thl, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L400;
	    }
	    for (j = 1; j <= 8; ++j) {
		i__2 = do_lio(&c__3, &c__1, (char *)&ig[j - 1], (ftnlen)
			sizeof(integer));
		if (i__2 != 0) {
		    goto L400;
		}
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L400;
	    }
L400:
	    ia = itype[ia];
	    for (j = 1; j <= 8; ++j) {
		if (ig[j - 1] != 0) {
		    kg = j;
		    ig[j - 1] = itype[ig[j - 1]];
		}
	    }
	    sort_(&kg, ig);
	    io___273.ciunit = *iprm;
	    s_wsfe(&io___273);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&pol, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&thl, (ftnlen)sizeof(doublereal));
	    i__2 = kg;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&ig[j - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	} else if (s_cmp(keyword, "PIATOM ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    ia = iclass[ia];
	    io___274.ciunit = *iprm;
	    s_wsfe(&io___274);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "PIBOND ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    prmsort_(&c__2, &ia, &ib, &c__0, &c__0, &c__0);
	    io___275.ciunit = *iprm;
	    s_wsfe(&io___275);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "PIBOND5 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    prmsort_(&c__2, &ia, &ib, &c__0, &c__0, &c__0);
	    io___276.ciunit = *iprm;
	    s_wsfe(&io___276);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "PIBOND4 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    ia = iclass[ia];
	    ib = iclass[ib];
	    prmsort_(&c__2, &ia, &ib, &c__0, &c__0, &c__0);
	    io___277.ciunit = *iprm;
	    s_wsfe(&io___277);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "METAL ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    ia = iclass[ia];
	    io___278.ciunit = *iprm;
	    s_wsfe(&io___278);
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, record + (next - 1), length - (next - 1));
	    e_wsfe();
	} else if (s_cmp(keyword, "BIOTYPE ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___279);
	    if (i__2 != 0) {
		goto L470;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L470;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L470;
	    }
	    getword_(record, string, &next, (ftnlen)120, (ftnlen)120);
	    getstring_(record, string, &next, (ftnlen)120, (ftnlen)120);
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__2 = s_rsli(&io___280);
	    if (i__2 != 0) {
		goto L470;
	    }
	    i__2 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__2 != 0) {
		goto L470;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L470;
	    }
L470:
	    if (ib != 0) {
		ib = itype[ib];
		if (ib == 0 && kt < 999) {
		    ib = 999;
		}
	    }
/* Computing MIN */
/* Computing MAX */
	    i__4 = 1, i__5 = 51 - next;
	    i__2 = 30, i__3 = max(i__4,i__5);
	    length = min(i__2,i__3);
	    io___281.ciunit = *iprm;
	    s_wsfe(&io___281);
/* Writing concatenation */
	    i__6[0] = next - 7, a__1[0] = record + 7;
	    i__6[1] = length, a__1[1] = blank;
	    s_cat(ch__1, a__1, i__6, &c__2, (ftnlen)150);
	    do_fio(&c__1, ch__1, next - 7 + length);
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else if (length == 0) {
	    io___282.ciunit = *iprm;
	    s_wsfe(&io___282);
	    e_wsfe();
	} else {
	    io___283.ciunit = *iprm;
	    s_wsfe(&io___283);
	    do_fio(&c__1, record, length);
	    e_wsfe();
	}
    }
    return 0;
} /* prmorder_ */

#undef prmline_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine prmsort  --  reorder atom types and classes  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "prmsort" places a list of atom type or class numbers into */
/*     canonical order for potential energy parameter definitions */


/* Subroutine */ int prmsort_(integer *index, integer *ia, integer *ib, 
	integer *ic, integer *id, integer *ie)
{
    static integer temp;



/*     put the atom type or class numbers into canonical order */

    if (*index == 2) {
	if (*ia > *ib) {
	    temp = *ia;
	    *ia = *ib;
	    *ib = temp;
	}
    } else if (*index == 3) {
	if (*ia > *ic) {
	    temp = *ia;
	    *ia = *ic;
	    *ic = temp;
	}
    } else if (*index == 4) {
	if (*ib > *ic || *ib == *ic && *ia > *id) {
	    temp = *ib;
	    *ib = *ic;
	    *ic = temp;
	    temp = *ia;
	    *ia = *id;
	    *id = temp;
	}
    } else if (*index == 5) {
	if (*ib > *id || *ib == *id && *ia > *ie) {
	    temp = *ib;
	    *ib = *id;
	    *id = temp;
	    temp = *ia;
	    *ia = *ie;
	    *ie = temp;
	}
    }
    return 0;
} /* prmsort_ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine polesort  --  sort multipoles by atom type  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "polesort" sorts a set of atomic multipole parameters based */
/*     on the atom types of centers involved */


/* Subroutine */ int polesort_(integer *iprm)
{
    /* Format strings */
    static char fmt_10[] = "(\002multipole \002,4i5,6x,f11.5)";
    static char fmt_30[] = "(\002multipole \002,3i5,11x,f11.5)";
    static char fmt_50[] = "(36x,3f11.5)";
    static char fmt_60[] = "(36x,f11.5)";
    static char fmt_70[] = "(36x,2f11.5)";
    static char fmt_80[] = "(36x,3f11.5)";

    /* System generated locals */
    address a__1[4];
    integer i__1, i__2[4], i__3;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j, n;
    static doublereal v1, v2, v3;
    static integer ia, ib, ic, id;
    static char pa[4], pb[4], pc[4], pd[4];
    static integer key[25000], line[25000], size;
    static char list[16*25000];
    static integer next;
    extern /* Subroutine */ int sort7_(integer *, char *, integer *, ftnlen);
    static char record[120];
    static integer length;
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int getnumb_(char *, integer *, integer *, ftnlen)
	    , numeral_(integer *, char *, integer *, ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static icilist io___305 = { 1, string, 1, 0, 120, 1 };
    static cilist io___307 = { 0, 0, 0, fmt_10, 0 };
    static icilist io___308 = { 1, string, 1, 0, 120, 1 };
    static cilist io___309 = { 0, 0, 0, fmt_30, 0 };
    static icilist io___310 = { 1, record, 1, 0, 120, 1 };
    static cilist io___313 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___314 = { 1, record, 1, 0, 120, 1 };
    static cilist io___315 = { 0, 0, 0, fmt_60, 0 };
    static icilist io___316 = { 1, record, 1, 0, 120, 1 };
    static cilist io___317 = { 0, 0, 0, fmt_70, 0 };
    static icilist io___318 = { 1, record, 1, 0, 120, 1 };
    static cilist io___319 = { 0, 0, 0, fmt_80, 0 };



#define list_ref(a_0,a_1) &list[(a_1)*16 + a_0 - 16]
#define prmline_ref(a_0,a_1) &params_1.prmline[(a_1)*120 + a_0 - 120]



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
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  params.i  --  contents of force field parameter file  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     nprm      number of nonblank lines in the parameter file */
/*     prmline   contents of each individual parameter file line */




/*     find and store atom types for the multipole parameters */

    n = 0;
    i__1 = params_1.nprm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(record, prmline_ref(0, i__), (ftnlen)120, (ftnlen)120);
	length = trimtext_(record, (ftnlen)120);
	next = 1;
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);
	if (s_cmp(keyword, "MULTIPOLE ", (ftnlen)10, (ftnlen)10) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &ib, &next, (ftnlen)120);
	    getnumb_(record, &ic, &next, (ftnlen)120);
	    getnumb_(record, &id, &next, (ftnlen)120);
	    ia = abs(ia);
	    ib = abs(ib);
	    ic = abs(ic);
	    id = abs(id);
	    size = 4;
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    numeral_(&id, pd, &size, (ftnlen)4);
	    ++n;
	    line[n - 1] = i__;
/* Writing concatenation */
	    i__2[0] = 4, a__1[0] = pa;
	    i__2[1] = 4, a__1[1] = pb;
	    i__2[2] = 4, a__1[2] = pc;
	    i__2[3] = 4, a__1[3] = pd;
	    s_cat(list_ref(0, n), a__1, i__2, &c__4, (ftnlen)16);
	}
    }

/*     sort the parameters based on the atom type numbers */

    sort7_(&n, list, key, (ftnlen)16);

/*     format and output the sorted multipole parameters */

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = line[key[i__ - 1] - 1];
	s_copy(record, prmline_ref(0, j), (ftnlen)120, (ftnlen)120);
	next = 1;
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	ia = 0;
	ib = 0;
	ic = 0;
	id = 0;
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	i__3 = s_rsli(&io___305);
	if (i__3 != 0) {
	    goto L20;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	if (i__3 != 0) {
	    goto L20;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	if (i__3 != 0) {
	    goto L20;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	if (i__3 != 0) {
	    goto L20;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	if (i__3 != 0) {
	    goto L20;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
	if (i__3 != 0) {
	    goto L20;
	}
	i__3 = e_rsli();
	if (i__3 != 0) {
	    goto L20;
	}
	io___307.ciunit = *iprm;
	s_wsfe(&io___307);
	do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
	e_wsfe();
	goto L40;
L20:
	i__3 = s_rsli(&io___308);
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = e_rsli();
	if (i__3 != 0) {
	    goto L90;
	}
	io___309.ciunit = *iprm;
	s_wsfe(&io___309);
	do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
	e_wsfe();
L40:
	++j;
	s_copy(record, prmline_ref(0, j), (ftnlen)120, (ftnlen)120);
	i__3 = s_rsli(&io___310);
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&v2, (ftnlen)sizeof(doublereal));
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&v3, (ftnlen)sizeof(doublereal));
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = e_rsli();
	if (i__3 != 0) {
	    goto L90;
	}
	io___313.ciunit = *iprm;
	s_wsfe(&io___313);
	do_fio(&c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&v3, (ftnlen)sizeof(doublereal));
	e_wsfe();
	++j;
	s_copy(record, prmline_ref(0, j), (ftnlen)120, (ftnlen)120);
	i__3 = s_rsli(&io___314);
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = e_rsli();
	if (i__3 != 0) {
	    goto L90;
	}
	io___315.ciunit = *iprm;
	s_wsfe(&io___315);
	do_fio(&c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
	e_wsfe();
	++j;
	s_copy(record, prmline_ref(0, j), (ftnlen)120, (ftnlen)120);
	i__3 = s_rsli(&io___316);
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&v2, (ftnlen)sizeof(doublereal));
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = e_rsli();
	if (i__3 != 0) {
	    goto L90;
	}
	io___317.ciunit = *iprm;
	s_wsfe(&io___317);
	do_fio(&c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(doublereal));
	e_wsfe();
	++j;
	s_copy(record, prmline_ref(0, j), (ftnlen)120, (ftnlen)120);
	i__3 = s_rsli(&io___318);
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&v2, (ftnlen)sizeof(doublereal));
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = do_lio(&c__5, &c__1, (char *)&v3, (ftnlen)sizeof(doublereal));
	if (i__3 != 0) {
	    goto L90;
	}
	i__3 = e_rsli();
	if (i__3 != 0) {
	    goto L90;
	}
	io___319.ciunit = *iprm;
	s_wsfe(&io___319);
	do_fio(&c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&v3, (ftnlen)sizeof(doublereal));
	e_wsfe();
L90:
	;
    }
    return 0;
} /* polesort_ */

#undef prmline_ref
#undef list_ref


/* Main program alias */ int prmedit_ () { MAIN__ (); return 0; }
