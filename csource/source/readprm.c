/* readprm.f -- translated by f2c (version 20050501).
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
    integer biotyp[10000];
    char forcefield[20];
} fields_;

#define fields_1 fields_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal anan[3000]	/* was [3][1000] */;
} kanang_;

#define kanang_1 kanang_

struct {
    doublereal acon[2000], acon5[500], acon4[500], acon3[500], aconf[500], 
	    ang[6000]	/* was [3][2000] */, ang5[1500]	/* was [3][500] */, 
	    ang4[1500]	/* was [3][500] */, ang3[1500]	/* was [3][500] */, 
	    angf[1000]	/* was [2][500] */;
    char ka[24000], ka5[6000], ka4[6000], ka3[6000], kaf[6000];
} kangs_;

#define kangs_1 kangs_

struct {
    doublereal weight[5000];
    integer atmcls[5000], atmnum[5000], ligand[5000];
    char symbol[15000], describe[120000];
} katoms_;

#define katoms_1 katoms_

struct {
    doublereal bcon[2000], blen[2000], bcon5[500], blen5[500], bcon4[500], 
	    blen4[500], bcon3[500], blen3[500], dlen[500];
    char kb[16000], kb5[4000], kb4[4000], kb3[4000], kel[6000];
} kbonds_;

#define kbonds_1 kbonds_

struct {
    doublereal chg[5000];
} kchrge_;

#define kchrge_1 kchrge_

struct {
    doublereal dpl[1000], dpl5[500], dpl4[500], dpl3[500], pos[1000], pos5[
	    500], pos4[500], pos3[500];
    char kd[8000], kd5[4000], kd4[4000], kd3[4000];
} kdipol_;

#define kdipol_1 kdipol_

struct {
    doublereal radhb[500], epshb[500];
    char khb[4000];
} khbond_;

#define khbond_1 khbond_

struct {
    doublereal dcon[500], tdi[500];
    char kdi[8000];
} kiprop_;

#define kiprop_1 kiprop_

struct {
    doublereal ti1[1000]	/* was [2][500] */, ti2[1000]	/* was [2][
	    500] */, ti3[1000]	/* was [2][500] */;
    char kti[8000];
} kitors_;

#define kitors_1 kitors_

struct {
    doublereal multip[26000]	/* was [13][2000] */;
    char mpaxis[16000], kmp[32000];
} kmulti_;

#define kmulti_1 kmulti_

struct {
    doublereal opbn[500];
    char kopb[8000];
} kopbnd_;

#define kopbnd_1 kopbnd_

struct {
    doublereal opds[500];
    char kopd[8000];
} kopdst_;

#define kopdst_1 kopdst_

struct {
    doublereal electron[1000], ionize[1000], repulse[1000], sslope[500], 
	    tslope[500], sslope5[200], tslope5[200], sslope4[200], tslope4[
	    200];
    char kpi[4000], kpi5[1600], kpi4[1600];
} korbs_;

#define korbs_1 korbs_

struct {
    doublereal ptcon[500];
    char kpt[4000];
} kpitor_;

#define kpitor_1 kpitor_

struct {
    doublereal polr[5000], athl[5000];
    integer pgrp[40000]	/* was [8][5000] */;
} kpolr_;

#define kpolr_1 kpolr_

struct {
    doublereal stbn[4000]	/* was [2][2000] */;
    char ksb[24000];
} kstbnd_;

#define kstbnd_1 kstbnd_

struct {
    doublereal btcon[1500]	/* was [3][500] */;
    char kbt[8000];
} ksttor_;

#define ksttor_1 ksttor_

struct {
    doublereal t1[4000]	/* was [2][2000] */, t2[4000]	/* was [2][2000] */, 
	    t3[4000]	/* was [2][2000] */, t4[4000]	/* was [2][2000] */, 
	    t5[4000]	/* was [2][2000] */, t6[4000]	/* was [2][2000] */, 
	    t15[1000]	/* was [2][500] */, t25[1000]	/* was [2][500] */, 
	    t35[1000]	/* was [2][500] */, t45[1000]	/* was [2][500] */, 
	    t55[1000]	/* was [2][500] */, t65[1000]	/* was [2][500] */, 
	    t14[1000]	/* was [2][500] */, t24[1000]	/* was [2][500] */, 
	    t34[1000]	/* was [2][500] */, t44[1000]	/* was [2][500] */, 
	    t54[1000]	/* was [2][500] */, t64[1000]	/* was [2][500] */;
    char kt[32000], kt5[8000], kt4[8000];
} ktorsn_;

#define ktorsn_1 ktorsn_

struct {
    doublereal ttx[3000]	/* was [30][100] */, tty[3000]	/* was [30][
	    100] */, tbf[90000]	/* was [900][100] */, tbx[90000]	/* 
	    was [900][100] */, tby[90000]	/* was [900][100] */, tbxy[
	    90000]	/* was [900][100] */;
    integer tnx[100], tny[100];
    char ktt[2000];
} ktrtor_;

#define ktrtor_1 ktrtor_

struct {
    doublereal ucon[2000], dst13[2000];
    char ku[24000];
} kurybr_;

#define kurybr_1 kurybr_

struct {
    doublereal rad[5000], eps[5000], rad4[5000], eps4[5000], reduct[5000];
} kvdws_;

#define kvdws_1 kvdws_

struct {
    doublereal radpr[500], epspr[500];
    char kvpr[4000];
} kvdwpr_;

#define kvdwpr_1 kvdwpr_

struct {
    integer nprm;
    char prmline[3000000];
} params_;

#define params_1 params_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static integer c__2 = 2;
static integer c__4 = 4;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine readprm  --  input of force field parameters  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "readprm" processes the potential energy parameter file */
/*     in order to define the default force field parameters */


/* Subroutine */ int readprm_(void)
{
    /* Format strings */
    static char fmt_10[] = "()";
    static char fmt_20[] = "()";
    static char fmt_30[] = "(a)";
    static char fmt_40[] = "(/,\002 READPRM  --  Too many Atom Types;\002"
	    ",\002 Increase MAXTYP\002)";
    static char fmt_50[] = "(/,\002 READPRM  --  Too many Atom Classes;\002"
	    ",\002 Increase MAXCLASS\002)";

    /* System generated locals */
    address a__1[2], a__2[3], a__3[4], a__4[5];
    integer i__1, i__2[2], i__3[3], i__4[4], i__5, i__6[5];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), do_fio(integer *, char *, ftnlen), s_rsli(icilist *), 
	    do_lio(integer *, integer *, char *, ftnlen), e_rsli(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int torphase_(integer *, doublereal *, doublereal 
	    *);
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j;
    extern /* Subroutine */ int getstring_(char *, char *, integer *, ftnlen, 
	    ftnlen);
    static doublereal bd, fc;
    static integer ia, ib, ic, id, ie, nb, na, nd;
    static doublereal ep, rd;
    static integer pg[8];
    static doublereal an;
    static integer nu, nt, nx, ny, ft[6];
    static doublereal pr, ds, dk, vd, cg, dp, aa1, ba1, ba2, aa2, aa3, ps, dl,
	     el, pt, iz, rp, ss;
    static integer nb5, nb4, nb3, na5, na4, na3, nd5, nd4, nd3, nt4, nt5;
    static doublereal an1, an2, an3, bt1, bt2, bt3, ts, vt[6], st[6], pl[13], 
	    tx[900], ty[900], tf[900];
    static char pa[4], pb[4], pc[4], pd[4], pe[4];
    static integer naf, nhb, ndi, lig, nel, nsb, nbt, npi, nti, nmp, cls, atn,
	     npt, nvp, ntt;
    static doublereal rdn, pol, thl;
    static char axt[8];
    static integer nxy, npi4, npi5, nopb, nopd, iprm;
    static doublereal wght;
    static integer size, next;
    extern /* Subroutine */ int sort9_(integer *, doublereal *), fatal_(void);
    static logical header;
    static char record[120];
    static integer length;
    static char string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), prmkey_(char *, 
	    ftnlen), getnumb_(char *, integer *, integer *, ftnlen), numeral_(
	    integer *, char *, integer *, ftnlen), getword_(char *, char *, 
	    integer *, ftnlen, ftnlen);
    static char keyword[20];
    extern /* Subroutine */ int gettext_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static cilist io___41 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___49 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___51 = { 1, string, 1, 0, 120, 1 };
    static icilist io___55 = { 1, string, 1, 0, 120, 1 };
    static icilist io___56 = { 1, string, 1, 0, 120, 1 };
    static icilist io___58 = { 1, string, 1, 0, 120, 1 };
    static icilist io___61 = { 1, string, 1, 0, 120, 1 };
    static icilist io___64 = { 1, string, 1, 0, 120, 1 };
    static icilist io___65 = { 1, string, 1, 0, 120, 1 };
    static icilist io___66 = { 1, string, 1, 0, 120, 1 };
    static icilist io___67 = { 1, string, 1, 0, 120, 1 };
    static icilist io___70 = { 1, string, 1, 0, 120, 1 };
    static icilist io___75 = { 1, string, 1, 0, 120, 1 };
    static icilist io___76 = { 1, string, 1, 0, 120, 1 };
    static icilist io___77 = { 1, string, 1, 0, 120, 1 };
    static icilist io___78 = { 1, string, 1, 0, 120, 1 };
    static icilist io___81 = { 1, string, 1, 0, 120, 1 };
    static icilist io___84 = { 1, string, 1, 0, 120, 1 };
    static icilist io___86 = { 1, string, 1, 0, 120, 1 };
    static icilist io___90 = { 1, string, 1, 0, 120, 1 };
    static icilist io___92 = { 1, string, 1, 0, 120, 1 };
    static icilist io___94 = { 1, string, 1, 0, 120, 1 };
    static icilist io___97 = { 1, string, 1, 0, 120, 1 };
    static icilist io___102 = { 1, string, 1, 0, 120, 1 };
    static icilist io___104 = { 1, string, 1, 0, 120, 1 };
    static icilist io___105 = { 1, string, 1, 0, 120, 1 };
    static icilist io___106 = { 1, string, 1, 0, 120, 1 };
    static icilist io___108 = { 1, string, 1, 0, 120, 1 };
    static icilist io___112 = { 1, string, 1, 0, 120, 1 };
    static icilist io___120 = { 1, string, 1, 0, 120, 1 };
    static icilist io___121 = { 1, record, 1, 0, 120, 1 };
    static icilist io___124 = { 1, string, 1, 0, 120, 1 };
    static icilist io___127 = { 1, string, 1, 0, 120, 1 };
    static icilist io___128 = { 1, string, 1, 0, 120, 1 };
    static icilist io___129 = { 1, string, 1, 0, 120, 1 };
    static icilist io___130 = { 1, string, 1, 0, 120, 1 };
    static icilist io___133 = { 1, string, 1, 0, 120, 1 };
    static icilist io___134 = { 1, string, 1, 0, 120, 1 };
    static icilist io___135 = { 1, record, 1, 0, 120, 1 };
    static icilist io___136 = { 1, record, 1, 0, 120, 1 };
    static icilist io___137 = { 1, record, 1, 0, 120, 1 };
    static icilist io___138 = { 1, record, 1, 0, 120, 1 };
    static icilist io___142 = { 1, string, 1, 0, 120, 1 };
    static icilist io___146 = { 1, string, 1, 0, 120, 1 };
    static icilist io___149 = { 1, string, 1, 0, 120, 1 };
    static icilist io___150 = { 1, string, 1, 0, 120, 1 };
    static icilist io___151 = { 1, string, 1, 0, 120, 1 };
    static icilist io___152 = { 1, string, 1, 0, 120, 1 };
    static icilist io___153 = { 1, string, 1, 0, 120, 1 };
    static icilist io___154 = { 1, string, 1, 0, 120, 1 };
    static cilist io___155 = { 0, 0, 0, fmt_40, 0 };



#define describe_ref(a_0,a_1) &katoms_1.describe[(a_1)*24 + a_0 - 24]
#define t1_ref(a_1,a_2) ktorsn_1.t1[(a_2)*2 + a_1 - 3]
#define t2_ref(a_1,a_2) ktorsn_1.t2[(a_2)*2 + a_1 - 3]
#define t3_ref(a_1,a_2) ktorsn_1.t3[(a_2)*2 + a_1 - 3]
#define t4_ref(a_1,a_2) ktorsn_1.t4[(a_2)*2 + a_1 - 3]
#define t5_ref(a_1,a_2) ktorsn_1.t5[(a_2)*2 + a_1 - 3]
#define t6_ref(a_1,a_2) ktorsn_1.t6[(a_2)*2 + a_1 - 3]
#define ka_ref(a_0,a_1) &kangs_1.ka[(a_1)*12 + a_0 - 12]
#define kb_ref(a_0,a_1) &kbonds_1.kb[(a_1)*8 + a_0 - 8]
#define kd_ref(a_0,a_1) &kdipol_1.kd[(a_1)*8 + a_0 - 8]
#define t14_ref(a_1,a_2) ktorsn_1.t14[(a_2)*2 + a_1 - 3]
#define t15_ref(a_1,a_2) ktorsn_1.t15[(a_2)*2 + a_1 - 3]
#define t25_ref(a_1,a_2) ktorsn_1.t25[(a_2)*2 + a_1 - 3]
#define t35_ref(a_1,a_2) ktorsn_1.t35[(a_2)*2 + a_1 - 3]
#define t45_ref(a_1,a_2) ktorsn_1.t45[(a_2)*2 + a_1 - 3]
#define t55_ref(a_1,a_2) ktorsn_1.t55[(a_2)*2 + a_1 - 3]
#define t65_ref(a_1,a_2) ktorsn_1.t65[(a_2)*2 + a_1 - 3]
#define t24_ref(a_1,a_2) ktorsn_1.t24[(a_2)*2 + a_1 - 3]
#define t34_ref(a_1,a_2) ktorsn_1.t34[(a_2)*2 + a_1 - 3]
#define t44_ref(a_1,a_2) ktorsn_1.t44[(a_2)*2 + a_1 - 3]
#define t54_ref(a_1,a_2) ktorsn_1.t54[(a_2)*2 + a_1 - 3]
#define t64_ref(a_1,a_2) ktorsn_1.t64[(a_2)*2 + a_1 - 3]
#define kt_ref(a_0,a_1) &ktorsn_1.kt[(a_1)*16 + a_0 - 16]
#define ku_ref(a_0,a_1) &kurybr_1.ku[(a_1)*12 + a_0 - 12]
#define ka3_ref(a_0,a_1) &kangs_1.ka3[(a_1)*12 + a_0 - 12]
#define ka4_ref(a_0,a_1) &kangs_1.ka4[(a_1)*12 + a_0 - 12]
#define ka5_ref(a_0,a_1) &kangs_1.ka5[(a_1)*12 + a_0 - 12]
#define kb5_ref(a_0,a_1) &kbonds_1.kb5[(a_1)*8 + a_0 - 8]
#define kb4_ref(a_0,a_1) &kbonds_1.kb4[(a_1)*8 + a_0 - 8]
#define kb3_ref(a_0,a_1) &kbonds_1.kb3[(a_1)*8 + a_0 - 8]
#define kd5_ref(a_0,a_1) &kdipol_1.kd5[(a_1)*8 + a_0 - 8]
#define kd4_ref(a_0,a_1) &kdipol_1.kd4[(a_1)*8 + a_0 - 8]
#define kd3_ref(a_0,a_1) &kdipol_1.kd3[(a_1)*8 + a_0 - 8]
#define ti1_ref(a_1,a_2) kitors_1.ti1[(a_2)*2 + a_1 - 3]
#define ti2_ref(a_1,a_2) kitors_1.ti2[(a_2)*2 + a_1 - 3]
#define ti3_ref(a_1,a_2) kitors_1.ti3[(a_2)*2 + a_1 - 3]
#define kt4_ref(a_0,a_1) &ktorsn_1.kt4[(a_1)*16 + a_0 - 16]
#define kt5_ref(a_0,a_1) &ktorsn_1.kt5[(a_1)*16 + a_0 - 16]
#define kaf_ref(a_0,a_1) &kangs_1.kaf[(a_1)*12 + a_0 - 12]
#define khb_ref(a_0,a_1) &khbond_1.khb[(a_1)*8 + a_0 - 8]
#define ang_ref(a_1,a_2) kangs_1.ang[(a_2)*3 + a_1 - 4]
#define kdi_ref(a_0,a_1) &kiprop_1.kdi[(a_1)*16 + a_0 - 16]
#define kel_ref(a_0,a_1) &kbonds_1.kel[(a_1)*12 + a_0 - 12]
#define tbf_ref(a_1,a_2) ktrtor_1.tbf[(a_2)*900 + a_1 - 901]
#define ksb_ref(a_0,a_1) &kstbnd_1.ksb[(a_1)*12 + a_0 - 12]
#define kbt_ref(a_0,a_1) &ksttor_1.kbt[(a_1)*16 + a_0 - 16]
#define kpi_ref(a_0,a_1) &korbs_1.kpi[(a_1)*8 + a_0 - 8]
#define kti_ref(a_0,a_1) &kitors_1.kti[(a_1)*16 + a_0 - 16]
#define kmp_ref(a_0,a_1) &kmulti_1.kmp[(a_1)*16 + a_0 - 16]
#define kpt_ref(a_0,a_1) &kpitor_1.kpt[(a_1)*8 + a_0 - 8]
#define ktt_ref(a_0,a_1) &ktrtor_1.ktt[(a_1)*20 + a_0 - 20]
#define ttx_ref(a_1,a_2) ktrtor_1.ttx[(a_2)*30 + a_1 - 31]
#define tty_ref(a_1,a_2) ktrtor_1.tty[(a_2)*30 + a_1 - 31]
#define ang3_ref(a_1,a_2) kangs_1.ang3[(a_2)*3 + a_1 - 4]
#define ang4_ref(a_1,a_2) kangs_1.ang4[(a_2)*3 + a_1 - 4]
#define ang5_ref(a_1,a_2) kangs_1.ang5[(a_2)*3 + a_1 - 4]
#define kpi4_ref(a_0,a_1) &korbs_1.kpi4[(a_1)*8 + a_0 - 8]
#define kpi5_ref(a_0,a_1) &korbs_1.kpi5[(a_1)*8 + a_0 - 8]
#define angf_ref(a_1,a_2) kangs_1.angf[(a_2)*2 + a_1 - 3]
#define anan_ref(a_1,a_2) kanang_1.anan[(a_2)*3 + a_1 - 4]
#define kopb_ref(a_0,a_1) &kopbnd_1.kopb[(a_1)*16 + a_0 - 16]
#define kopd_ref(a_0,a_1) &kopdst_1.kopd[(a_1)*16 + a_0 - 16]
#define stbn_ref(a_1,a_2) kstbnd_1.stbn[(a_2)*2 + a_1 - 3]
#define pgrp_ref(a_1,a_2) kpolr_1.pgrp[(a_2)*8 + a_1 - 9]
#define kvpr_ref(a_0,a_1) &kvdwpr_1.kvpr[(a_1)*8 + a_0 - 8]
#define btcon_ref(a_1,a_2) ksttor_1.btcon[(a_2)*3 + a_1 - 4]
#define symbol_ref(a_0,a_1) &katoms_1.symbol[(a_1)*3 + a_0 - 3]
#define multip_ref(a_1,a_2) kmulti_1.multip[(a_2)*13 + a_1 - 14]
#define mpaxis_ref(a_0,a_1) &kmulti_1.mpaxis[(a_1)*8 + a_0 - 8]
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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  fields.i  --  molecular mechanics force field description  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     biotyp       force field atom type of each biopolymer type */
/*     forcefield   string used to describe the current forcefield */




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

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kanang.i  --  forcefield parameters for angle-angle terms  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     anan   angle-angle cross term parameters for each atom class */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kangs.i  --  forcefield parameters for bond angle bending  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxna    maximum number of harmonic angle bend parameter entries */
/*     maxna5   maximum number of 5-membered ring angle bend entries */
/*     maxna4   maximum number of 4-membered ring angle bend entries */
/*     maxna3   maximum number of 3-membered ring angle bend entries */
/*     maxnaf   maximum number of Fourier angle bend parameter entries */

/*     acon     force constant parameters for harmonic angle bends */
/*     acon5    force constant parameters for 5-ring angle bends */
/*     acon4    force constant parameters for 4-ring angle bends */
/*     acon3    force constant parameters for 3-ring angle bends */
/*     aconf    force constant parameters for Fourier angle bends */
/*     ang      bond angle parameters for harmonic angle bends */
/*     ang5     bond angle parameters for 5-ring angle bends */
/*     ang4     bond angle parameters for 4-ring angle bends */
/*     ang3     bond angle parameters for 3-ring angle bends */
/*     angf     phase shift angle and periodicity for Fourier bends */
/*     ka       string of atom classes for harmonic angle bends */
/*     ka5      string of atom classes for 5-ring angle bends */
/*     ka4      string of atom classes for 4-ring angle bends */
/*     ka3      string of atom classes for 3-ring angle bends */
/*     kaf      string of atom classes for Fourier angle bends */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  katoms.i  --  forcefield parameters for the atom types  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     weight     average atomic mass of each atom type */
/*     atmcls     atom class number for each of the atom types */
/*     atmnum     atomic number for each of the atom types */
/*     ligand     number of atoms to be attached to each atom type */
/*     symbol     modified atomic symbol for each atom type */
/*     describe   string identifying each of the atom types */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  kbonds.i  --  forcefield parameters for bond stretching  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxnb   maximum number of bond stretch parameter entries */
/*     maxnb5  maximum number of 5-membered ring bond stretch entries */
/*     maxnb4  maximum number of 4-membered ring bond stretch entries */
/*     maxnb3  maximum number of 3-membered ring bond stretch entries */
/*     maxnel  maximum number of electronegativity bond corrections */

/*     bcon    force constant parameters for harmonic bond stretch */
/*     blen    bond length parameters for harmonic bond stretch */
/*     bcon5   force constant parameters for 5-ring bond stretch */
/*     blen5   bond length parameters for 5-ring bond stretch */
/*     bcon4   force constant parameters for 4-ring bond stretch */
/*     blen4   bond length parameters for 4-ring bond stretch */
/*     bcon3   force constant parameters for 3-ring bond stretch */
/*     blen3   bond length parameters for 3-ring bond stretch */
/*     dlen    electronegativity bond length correction parameters */
/*     kb      string of atom classes for harmonic bond stretch */
/*     kb5     string of atom classes for 5-ring bond stretch */
/*     kb4     string of atom classes for 4-ring bond stretch */
/*     kb3     string of atom classes for 3-ring bond stretch */
/*     kel     string of atom classes for electronegativity corrections */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  kchrge.i  --  forcefield parameters for partial charges  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     chg   partial charge parameters for each atom type */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  kdipol.i  --  forcefield parameters for bond dipoles  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxnd    maximum number of bond dipole parameter entries */
/*     maxnd5   maximum number of 5-membered ring dipole entries */
/*     maxnd4   maximum number of 4-membered ring dipole entries */
/*     maxnd3   maximum number of 3-membered ring dipole entries */

/*     dpl      dipole moment parameters for bond dipoles */
/*     dpl5     dipole moment parameters for 5-ring dipoles */
/*     dpl4     dipole moment parameters for 4-ring dipoles */
/*     dpl3     dipole moment parameters for 3-ring dipoles */
/*     pos      dipole position parameters for bond dipoles */
/*     pos5     dipole position parameters for 5-ring dipoles */
/*     pos4     dipole position parameters for 4-ring dipoles */
/*     pos3     dipole position parameters for 3-ring dipoles */
/*     kd       string of atom classes for bond dipoles */
/*     kd5      string of atom classes for 5-ring dipoles */
/*     kd4      string of atom classes for 4-ring dipoles */
/*     kd3      string of atom classes for 3-ring dipoles */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  khbond.i  --  forcefield parameters for H-bonding terms  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxnhb   maximum number of hydrogen bonding pair entries */

/*     radhb    radius parameter for hydrogen bonding pairs */
/*     epshb    well depth parameter for hydrogen bonding pairs */
/*     khb      string of atom types for hydrogen bonding pairs */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kiprop.i  --  forcefield parameters for improper dihedral  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxndi   maximum number of improper dihedral parameter entries */

/*     dcon     force constant parameters for improper dihedrals */
/*     tdi      ideal dihedral angle values for improper dihedrals */
/*     kdi      string of atom classes for improper dihedral angles */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kitors.i  --  forcefield parameters for improper torsions  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxnti   maximum number of improper torsion parameter entries */

/*     ti1      torsional parameters for improper 1-fold rotation */
/*     ti2      torsional parameters for improper 2-fold rotation */
/*     ti3      torsional parameters for improper 3-fold rotation */
/*     kti      string of atom classes for improper torsional parameters */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kmulti.i  --  forcefield parameters for atomic multipoles  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxnmp   maximum number of atomic multipole parameter entries */

/*     multip   atomic monopole, dipole and quadrupole values */
/*     mpaxis   type of local axis definition for atomic multipoles */
/*     kmp      string of atom types for atomic multipoles */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kopbnd.i  --  forcefield parameters for out-of-plane bend  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxnopb   maximum number of out-of-plane bending entries */

/*     opbn      force constant parameters for out-of-plane bending */
/*     kopb      string of atom classes for out-of-plane bending */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kopdst.i  --  forcefield parameters for out-plane distance  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     maxnopd   maximum number of out-of-plane distance entries */

/*     opds      force constant parameters for out-of-plane distance */
/*     kopd      string of atom classes for out-of-plane distance */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  korbs.i  --  forcefield parameters for pisystem orbitals  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxnpi     maximum number of pisystem bond parameter entries */
/*     maxnpi5    maximum number of 5-membered ring pibond entries */
/*     maxnpi4    maximum number of 4-membered ring pibond entries */

/*     electron   number of pi-electrons for each atom class */
/*     ionize     ionization potential for each atom class */
/*     repulse    repulsion integral value for each atom class */
/*     sslope     slope for bond stretch vs. pi-bond order */
/*     tslope     slope for 2-fold torsion vs. pi-bond order */
/*     sslope5    slope for 5-ring bond stretch vs. pi-bond order */
/*     tslope5    slope for 5-ring 2-fold torsion vs. pi-bond order */
/*     sslope4    slope for 4-ring bond stretch vs. pi-bond order */
/*     tslope4    slope for 4-ring 2-fold torsion vs. pi-bond order */
/*     kpi        string of atom classes for pisystem bonds */
/*     kpi5       string of atom classes for 5-ring pisystem bonds */
/*     kpi4       string of atom classes for 4-ring pisystem bonds */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kpitor.i  --  forcefield parameters for pi-orbit torsions  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxnpt   maximum number of pi-orbital torsion parameter entries */

/*     ptcon    force constant parameters for pi-orbital torsions */
/*     kpt      string of atom classes for pi-orbital torsion terms */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  kpolr.i  --  forcefield parameters for polarizability  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     polr   dipole polarizability parameters for each atom type */
/*     athl   Thole polarizability damping value for each atom type */
/*     pgrp   connected types in polarization group of each atom type */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  kstbnd.i  --  forcefield parameters for stretch-bend  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxnsb   maximum number of stretch-bend parameter entries */

/*     stbn     force constant parameters for stretch-bend terms */
/*     ksb      string of atom classes for stretch-bend terms */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  ksttor.i  --  forcefield parameters for stretch-torsions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxnbt   maximum number of stretch-torsion parameter entries */

/*     btcon    force constant parameters for stretch-torsion */
/*     kbt      string of atom classes for stretch-torsion terms */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  ktorsn.i  --  forcefield parameters for torsional angles  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxnt    maximum number of torsional angle parameter entries */
/*     maxnt5   maximum number of 5-membered ring torsion entries */
/*     maxnt4   maximum number of 4-membered ring torsion entries */

/*     t1       torsional parameters for standard 1-fold rotation */
/*     t2       torsional parameters for standard 2-fold rotation */
/*     t3       torsional parameters for standard 3-fold rotation */
/*     t4       torsional parameters for standard 4-fold rotation */
/*     t5       torsional parameters for standard 5-fold rotation */
/*     t6       torsional parameters for standard 6-fold rotation */
/*     t15      torsional parameters for 1-fold rotation in 5-ring */
/*     t25      torsional parameters for 2-fold rotation in 5-ring */
/*     t35      torsional parameters for 3-fold rotation in 5-ring */
/*     t45      torsional parameters for 4-fold rotation in 5-ring */
/*     t55      torsional parameters for 5-fold rotation in 5-ring */
/*     t65      torsional parameters for 6-fold rotation in 5-ring */
/*     t14      torsional parameters for 1-fold rotation in 4-ring */
/*     t24      torsional parameters for 2-fold rotation in 4-ring */
/*     t34      torsional parameters for 3-fold rotation in 4-ring */
/*     t44      torsional parameters for 4-fold rotation in 4-ring */
/*     t54      torsional parameters for 5-fold rotation in 4-ring */
/*     t64      torsional parameters for 6-fold rotation in 4-ring */
/*     kt       string of atom classes for torsional angles */
/*     kt5      string of atom classes for 5-ring torsions */
/*     kt4      string of atom classes for 4-ring torsions */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  ktrtor.i  --  forcefield parameters for torsion-torsions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxntt    maximum number of torsion-torsion parameter entries */
/*     maxtgrd   maximum dimension of torsion-torsion spline grid */
/*     maxtgrd2  maximum number of torsion-torsion spline grid points */

/*     ttx       angle values for first torsion of spline grid */
/*     tty       angle values for second torsion of spline grid */
/*     tbf       function values at points on spline grid */
/*     tbx       gradient over first torsion of spline grid */
/*     tby       gradient over second torsion of spline grid */
/*     tbxy      Hessian cross components over spline grid */
/*     tnx       number of columns in torsion-torsion spline grid */
/*     tny       number of rows in torsion-torsion spline grid */
/*     ktt       string of torsion-torsion atom classes */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kurybr.i  --  forcefield parameters for Urey-Bradley terms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     maxnu   maximum number of Urey-Bradley parameter entries */

/*     ucon    force constant parameters for Urey-Bradley terms */
/*     dst13   ideal 1-3 distance parameters for Urey-Bradley terms */
/*     ku      string of atom classes for Urey-Bradley terms */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  kvdws.i  --  forcefield parameters for van der Waals terms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     rad      van der Waals radius parameter for each atom type */
/*     eps      van der Waals well depth parameter for each atom type */
/*     rad4     van der Waals radius parameter in 1-4 interactions */
/*     eps4     van der Waals well depth parameter in 1-4 interactions */
/*     reduct   van der Waals reduction factor for each atom type */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  kvdwpr.i  --  forcefield parameters for special vdw terms  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     maxnvp   maximum number of special van der Waals pair entries */

/*     radpr    radius parameter for special van der Waals pairs */
/*     epspr    well depth parameter for special van der Waals pairs */
/*     kvpr     string of atom classes for special van der Waals pairs */




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




/*     initialize the counters for some parameter types */

    nvp = 0;
    nhb = 0;
    nb = 0;
    nb5 = 0;
    nb4 = 0;
    nb3 = 0;
    nel = 0;
    na = 0;
    na5 = 0;
    na4 = 0;
    na3 = 0;
    naf = 0;
    nsb = 0;
    nu = 0;
    nopb = 0;
    nopd = 0;
    ndi = 0;
    nti = 0;
    nt = 0;
    nt5 = 0;
    nt4 = 0;
    npt = 0;
    nbt = 0;
    ntt = 0;
    nd = 0;
    nd5 = 0;
    nd4 = 0;
    nd3 = 0;
    nmp = 0;
    npi = 0;
    npi5 = 0;
    npi4 = 0;

/*     number of characters in an atom number text string */

    size = 4;

/*     set blank line header before echoed comment lines */

    header = TRUE_;

/*     process each line of the parameter file, first */
/*     extract the keyword at the start of each line */

    iprm = 0;
    while(iprm < params_1.nprm) {
	++iprm;
	s_copy(record, prmline_ref(0, iprm), (ftnlen)120, (ftnlen)120);
	next = 1;
	gettext_(record, keyword, &next, (ftnlen)120, (ftnlen)20);
	upcase_(keyword, (ftnlen)20);

/*     check for a force field modification keyword */

	prmkey_(record, (ftnlen)120);

/*     comment line to be echoed to the output */

	if (s_cmp(keyword, "ECHO ", (ftnlen)5, (ftnlen)5) == 0) {
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    length = trimtext_(string, (ftnlen)120);
	    if (header) {
		header = FALSE_;
		io___41.ciunit = iounit_1.iout;
		s_wsfe(&io___41);
		e_wsfe();
	    }
	    if (length == 0) {
		io___42.ciunit = iounit_1.iout;
		s_wsfe(&io___42);
		e_wsfe();
	    } else {
		io___43.ciunit = iounit_1.iout;
		s_wsfe(&io___43);
		do_fio(&c__1, string, length);
		e_wsfe();
	    }

/*     atom type definitions and parameters */

	} else if (s_cmp(keyword, "ATOM ", (ftnlen)5, (ftnlen)5) == 0) {
	    ia = 0;
	    cls = 0;
	    atn = 0;
	    wght = 0.;
	    lig = 0;
	    getnumb_(record, &ia, &next, (ftnlen)120);
	    getnumb_(record, &cls, &next, (ftnlen)120);
	    if (cls == 0) {
		cls = ia;
	    }
	    katoms_1.atmcls[ia - 1] = cls;
	    if (ia >= 5000) {
		io___49.ciunit = iounit_1.iout;
		s_wsfe(&io___49);
		e_wsfe();
		fatal_();
	    } else if (cls >= 1000) {
		io___50.ciunit = iounit_1.iout;
		s_wsfe(&io___50);
		e_wsfe();
		fatal_();
	    }
	    if (ia != 0) {
		gettext_(record, symbol_ref(0, ia), &next, (ftnlen)120, (
			ftnlen)3);
		getstring_(record, describe_ref(0, ia), &next, (ftnlen)120, (
			ftnlen)24);
		s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next 
			- 1));
		i__1 = s_rsli(&io___51);
		if (i__1 != 0) {
		    goto L60;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&atn, (ftnlen)sizeof(
			integer));
		if (i__1 != 0) {
		    goto L60;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&wght, (ftnlen)sizeof(
			doublereal));
		if (i__1 != 0) {
		    goto L60;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&lig, (ftnlen)sizeof(
			integer));
		if (i__1 != 0) {
		    goto L60;
		}
		i__1 = e_rsli();
		if (i__1 != 0) {
		    goto L60;
		}
L60:
		katoms_1.atmnum[ia - 1] = atn;
		katoms_1.weight[ia - 1] = wght;
		katoms_1.ligand[ia - 1] = lig;
	    }

/*     van der Waals parameters for individual atom types */

	} else if (s_cmp(keyword, "VDW ", (ftnlen)4, (ftnlen)4) == 0) {
	    ia = 0;
	    rd = 0.;
	    ep = 0.;
	    rdn = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___55);
	    if (i__1 != 0) {
		goto L70;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L70;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&rd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L70;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ep, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L70;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&rdn, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L70;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L70;
	    }
L70:
	    if (ia != 0) {
		kvdws_1.rad[ia - 1] = rd;
		kvdws_1.eps[ia - 1] = ep;
		kvdws_1.reduct[ia - 1] = rdn;
	    }

/*     van der Waals 1-4 parameters for individual atom types */

	} else if (s_cmp(keyword, "VDW14 ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    rd = 0.;
	    ep = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___56);
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
	    if (ia != 0) {
		kvdws_1.rad4[ia - 1] = rd;
		kvdws_1.eps4[ia - 1] = ep;
	    }

/*     van der Waals parameters for specific atom pairs */

	} else if (s_cmp(keyword, "VDWPR ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    rd = 0.;
	    ep = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___58);
	    if (i__1 != 0) {
		goto L90;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L90;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L90;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&rd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L90;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ep, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L90;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L90;
	    }
L90:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++nvp;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(kvpr_ref(0, nvp), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(kvpr_ref(0, nvp), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    kvdwpr_1.radpr[nvp - 1] = rd;
	    kvdwpr_1.epspr[nvp - 1] = ep;

/*     van der Waals parameters for hydrogen bonding pairs */

	} else if (s_cmp(keyword, "HBOND ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    rd = 0.;
	    ep = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___61);
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
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++nhb;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(khb_ref(0, nhb), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(khb_ref(0, nhb), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    khbond_1.radhb[nhb - 1] = rd;
	    khbond_1.epshb[nhb - 1] = ep;

/*     bond stretching parameters */

	} else if (s_cmp(keyword, "BOND ", (ftnlen)5, (ftnlen)5) == 0) {
	    ia = 0;
	    ib = 0;
	    fc = 0.;
	    bd = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___64);
	    if (i__1 != 0) {
		goto L110;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L110;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L110;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L110;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&bd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L110;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L110;
	    }
L110:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++nb;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(kb_ref(0, nb), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(kb_ref(0, nb), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    kbonds_1.bcon[nb - 1] = fc;
	    kbonds_1.blen[nb - 1] = bd;

/*     bond stretching parameters for 5-membered rings */

	} else if (s_cmp(keyword, "BOND5 ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    fc = 0.;
	    bd = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___65);
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
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&bd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L120;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L120;
	    }
L120:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++nb5;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(kb5_ref(0, nb5), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(kb5_ref(0, nb5), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    kbonds_1.bcon5[nb5 - 1] = fc;
	    kbonds_1.blen5[nb5 - 1] = bd;

/*     bond stretching parameters for 4-membered rings */

	} else if (s_cmp(keyword, "BOND4 ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    fc = 0.;
	    bd = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___66);
	    if (i__1 != 0) {
		goto L130;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L130;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L130;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L130;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&bd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L130;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L130;
	    }
L130:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++nb4;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(kb4_ref(0, nb4), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(kb4_ref(0, nb4), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    kbonds_1.bcon4[nb4 - 1] = fc;
	    kbonds_1.blen4[nb4 - 1] = bd;

/*     bond stretching parameters for 3-membered rings */

	} else if (s_cmp(keyword, "BOND3 ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    fc = 0.;
	    bd = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___67);
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
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++nb3;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(kb3_ref(0, nb3), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(kb3_ref(0, nb3), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    kbonds_1.bcon3[nb3 - 1] = fc;
	    kbonds_1.blen3[nb3 - 1] = bd;

/*     electronegativity bond length correction parameters */

	} else if (s_cmp(keyword, "ELECTNEG ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    dl = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___70);
	    if (i__1 != 0) {
		goto L150;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L150;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L150;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L150;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&dl, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L150;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L150;
	    }
L150:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    ++nel;
	    if (ia <= ic) {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pa;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pc;
		s_cat(kel_ref(0, nel), a__2, i__3, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pc;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pa;
		s_cat(kel_ref(0, nel), a__2, i__3, &c__3, (ftnlen)12);
	    }
	    kbonds_1.dlen[nel - 1] = dl;

/*     bond angle bending parameters */

	} else if (s_cmp(keyword, "ANGLE ", (ftnlen)6, (ftnlen)6) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    an1 = 0.;
	    an2 = 0.;
	    an3 = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___75);
	    if (i__1 != 0) {
		goto L160;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L160;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L160;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L160;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L160;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L160;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L160;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L160;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L160;
	    }
L160:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    ++na;
	    if (ia <= ic) {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pa;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pc;
		s_cat(ka_ref(0, na), a__2, i__3, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pc;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pa;
		s_cat(ka_ref(0, na), a__2, i__3, &c__3, (ftnlen)12);
	    }
	    kangs_1.acon[na - 1] = fc;
	    ang_ref(1, na) = an1;
	    ang_ref(2, na) = an2;
	    ang_ref(3, na) = an3;

/*     angle bending parameters for 5-membered rings */

	} else if (s_cmp(keyword, "ANGLE5 ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    an1 = 0.;
	    an2 = 0.;
	    an3 = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___76);
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
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L170;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L170;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L170;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L170;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L170;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L170;
	    }
L170:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    ++na5;
	    if (ia <= ic) {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pa;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pc;
		s_cat(ka5_ref(0, na5), a__2, i__3, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pc;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pa;
		s_cat(ka5_ref(0, na5), a__2, i__3, &c__3, (ftnlen)12);
	    }
	    kangs_1.acon5[na5 - 1] = fc;
	    ang5_ref(1, na5) = an1;
	    ang5_ref(2, na5) = an2;
	    ang5_ref(3, na5) = an3;

/*     angle bending parameters for 4-membered rings */

	} else if (s_cmp(keyword, "ANGLE4 ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    an1 = 0.;
	    an2 = 0.;
	    an3 = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___77);
	    if (i__1 != 0) {
		goto L180;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L180;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L180;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L180;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L180;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L180;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L180;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L180;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L180;
	    }
L180:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    ++na4;
	    if (ia <= ic) {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pa;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pc;
		s_cat(ka4_ref(0, na4), a__2, i__3, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pc;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pa;
		s_cat(ka4_ref(0, na4), a__2, i__3, &c__3, (ftnlen)12);
	    }
	    kangs_1.acon4[na4 - 1] = fc;
	    ang4_ref(1, na4) = an1;
	    ang4_ref(2, na4) = an2;
	    ang4_ref(3, na4) = an3;

/*     angle bending parameters for 3-membered rings */

	} else if (s_cmp(keyword, "ANGLE3 ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    an1 = 0.;
	    an2 = 0.;
	    an3 = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___78);
	    if (i__1 != 0) {
		goto L190;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L190;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L190;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L190;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L190;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L190;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L190;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L190;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L190;
	    }
L190:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    ++na3;
	    if (ia <= ic) {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pa;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pc;
		s_cat(ka3_ref(0, na3), a__2, i__3, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pc;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pa;
		s_cat(ka3_ref(0, na3), a__2, i__3, &c__3, (ftnlen)12);
	    }
	    kangs_1.acon3[na3 - 1] = fc;
	    ang3_ref(1, na3) = an1;
	    ang3_ref(2, na3) = an2;
	    ang3_ref(3, na3) = an3;

/*     Fourier bond angle bending parameters */

	} else if (s_cmp(keyword, "ANGLEF ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    an = 0.;
	    pr = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___81);
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
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L200;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L200;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&an, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L200;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pr, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L200;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L200;
	    }
L200:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    ++naf;
	    if (ia <= ic) {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pa;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pc;
		s_cat(kaf_ref(0, naf), a__2, i__3, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pc;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pa;
		s_cat(kaf_ref(0, naf), a__2, i__3, &c__3, (ftnlen)12);
	    }
	    kangs_1.aconf[naf - 1] = fc;
	    angf_ref(1, naf) = an;
	    angf_ref(2, naf) = pr;

/*     stretch-bend parameters */

	} else if (s_cmp(keyword, "STRBND ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    ba1 = 0.;
	    ba2 = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___84);
	    if (i__1 != 0) {
		goto L210;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L210;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L210;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L210;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ba1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L210;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ba2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L210;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L210;
	    }
L210:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    ++nsb;
	    if (ia <= ic) {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pa;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pc;
		s_cat(ksb_ref(0, nsb), a__2, i__3, &c__3, (ftnlen)12);
		stbn_ref(1, nsb) = ba1;
		stbn_ref(2, nsb) = ba2;
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pc;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pa;
		s_cat(ksb_ref(0, nsb), a__2, i__3, &c__3, (ftnlen)12);
		stbn_ref(1, nsb) = ba2;
		stbn_ref(2, nsb) = ba1;
	    }

/*     Urey-Bradley parameters */

	} else if (s_cmp(keyword, "UREYBRAD ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    fc = 0.;
	    ds = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___86);
	    if (i__1 != 0) {
		goto L220;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L220;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L220;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L220;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L220;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ds, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L220;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L220;
	    }
L220:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    ++nu;
	    if (ia <= ic) {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pa;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pc;
		s_cat(ku_ref(0, nu), a__2, i__3, &c__3, (ftnlen)12);
	    } else {
/* Writing concatenation */
		i__3[0] = 4, a__2[0] = pc;
		i__3[1] = 4, a__2[1] = pb;
		i__3[2] = 4, a__2[2] = pa;
		s_cat(ku_ref(0, nu), a__2, i__3, &c__3, (ftnlen)12);
	    }
	    kurybr_1.ucon[nu - 1] = fc;
	    kurybr_1.dst13[nu - 1] = ds;

/*     angle-angle parameters */

	} else if (s_cmp(keyword, "ANGANG ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    aa1 = 0.;
	    aa2 = 0.;
	    aa3 = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___90);
	    if (i__1 != 0) {
		goto L230;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L230;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&aa1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L230;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&aa2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L230;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&aa3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L230;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L230;
	    }
L230:
	    if (ia != 0) {
		anan_ref(1, ia) = aa1;
		anan_ref(2, ia) = aa2;
		anan_ref(3, ia) = aa3;
	    }

/*     out-of-plane bend parameters */

	} else if (s_cmp(keyword, "OPBEND ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    fc = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___92);
	    if (i__1 != 0) {
		goto L240;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L240;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L240;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L240;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L240;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L240;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L240;
	    }
L240:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    numeral_(&id, pd, &size, (ftnlen)4);
	    ++nopb;
/* Writing concatenation */
	    i__4[0] = 4, a__3[0] = pa;
	    i__4[1] = 4, a__3[1] = pb;
	    i__4[2] = 4, a__3[2] = pc;
	    i__4[3] = 4, a__3[3] = pd;
	    s_cat(kopb_ref(0, nopb), a__3, i__4, &c__4, (ftnlen)16);
	    kopbnd_1.opbn[nopb - 1] = fc;

/*     out-of-plane distance parameters */

	} else if (s_cmp(keyword, "OPDIST ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    fc = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___94);
	    if (i__1 != 0) {
		goto L250;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L250;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L250;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L250;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L250;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&fc, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L250;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L250;
	    }
L250:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    numeral_(&id, pd, &size, (ftnlen)4);
	    ++nopd;
/* Writing concatenation */
	    i__4[0] = 4, a__3[0] = pa;
	    i__4[1] = 4, a__3[1] = pb;
	    i__4[2] = 4, a__3[2] = pc;
	    i__4[3] = 4, a__3[3] = pd;
	    s_cat(kopd_ref(0, nopd), a__3, i__4, &c__4, (ftnlen)16);
	    kopdst_1.opds[nopd - 1] = fc;

/*     improper dihedral parameters */

	} else if (s_cmp(keyword, "IMPROPER ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    dk = 0.;
	    vd = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___97);
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
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L260;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&dk, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L260;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&vd, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L260;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L260;
	    }
L260:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    numeral_(&id, pd, &size, (ftnlen)4);
	    ++ndi;
/* Writing concatenation */
	    i__4[0] = 4, a__3[0] = pa;
	    i__4[1] = 4, a__3[1] = pb;
	    i__4[2] = 4, a__3[2] = pc;
	    i__4[3] = 4, a__3[3] = pd;
	    s_cat(kdi_ref(0, ndi), a__3, i__4, &c__4, (ftnlen)16);
	    kiprop_1.dcon[ndi - 1] = dk;
	    kiprop_1.tdi[ndi - 1] = vd;

/*     improper torsional parameters */

	} else if (s_cmp(keyword, "IMPTORS ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    for (i__ = 1; i__ <= 6; ++i__) {
		vt[i__ - 1] = 0.;
		st[i__ - 1] = 0.;
		ft[i__ - 1] = 0;
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___102);
	    if (i__1 != 0) {
		goto L270;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L270;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L270;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L270;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L270;
	    }
	    for (j = 1; j <= 6; ++j) {
		i__1 = do_lio(&c__5, &c__1, (char *)&vt[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L270;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&st[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L270;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&ft[j - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L270;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L270;
	    }
L270:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    numeral_(&id, pd, &size, (ftnlen)4);
	    ++nti;
/* Writing concatenation */
	    i__4[0] = 4, a__3[0] = pa;
	    i__4[1] = 4, a__3[1] = pb;
	    i__4[2] = 4, a__3[2] = pc;
	    i__4[3] = 4, a__3[3] = pd;
	    s_cat(kti_ref(0, nti), a__3, i__4, &c__4, (ftnlen)16);
	    torphase_(ft, vt, st);
	    ti1_ref(1, nti) = vt[0];
	    ti1_ref(2, nti) = st[0];
	    ti2_ref(1, nti) = vt[1];
	    ti2_ref(2, nti) = st[1];
	    ti3_ref(1, nti) = vt[2];
	    ti3_ref(2, nti) = st[2];

/*     torsional parameters */

	} else if (s_cmp(keyword, "TORSION ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    for (i__ = 1; i__ <= 6; ++i__) {
		vt[i__ - 1] = 0.;
		st[i__ - 1] = 0.;
		ft[i__ - 1] = 0;
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___104);
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
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L280;
	    }
	    for (j = 1; j <= 6; ++j) {
		i__1 = do_lio(&c__5, &c__1, (char *)&vt[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L280;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&st[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L280;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&ft[j - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L280;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L280;
	    }
L280:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    numeral_(&id, pd, &size, (ftnlen)4);
	    ++nt;
	    if (ib < ic) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pa;
		i__4[1] = 4, a__3[1] = pb;
		i__4[2] = 4, a__3[2] = pc;
		i__4[3] = 4, a__3[3] = pd;
		s_cat(kt_ref(0, nt), a__3, i__4, &c__4, (ftnlen)16);
	    } else if (ic < ib) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pd;
		i__4[1] = 4, a__3[1] = pc;
		i__4[2] = 4, a__3[2] = pb;
		i__4[3] = 4, a__3[3] = pa;
		s_cat(kt_ref(0, nt), a__3, i__4, &c__4, (ftnlen)16);
	    } else if (ia <= id) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pa;
		i__4[1] = 4, a__3[1] = pb;
		i__4[2] = 4, a__3[2] = pc;
		i__4[3] = 4, a__3[3] = pd;
		s_cat(kt_ref(0, nt), a__3, i__4, &c__4, (ftnlen)16);
	    } else if (id < ia) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pd;
		i__4[1] = 4, a__3[1] = pc;
		i__4[2] = 4, a__3[2] = pb;
		i__4[3] = 4, a__3[3] = pa;
		s_cat(kt_ref(0, nt), a__3, i__4, &c__4, (ftnlen)16);
	    }
	    torphase_(ft, vt, st);
	    t1_ref(1, nt) = vt[0];
	    t1_ref(2, nt) = st[0];
	    t2_ref(1, nt) = vt[1];
	    t2_ref(2, nt) = st[1];
	    t3_ref(1, nt) = vt[2];
	    t3_ref(2, nt) = st[2];
	    t4_ref(1, nt) = vt[3];
	    t4_ref(2, nt) = st[3];
	    t5_ref(1, nt) = vt[4];
	    t5_ref(2, nt) = st[4];
	    t6_ref(1, nt) = vt[5];
	    t6_ref(2, nt) = st[5];

/*     torsional parameters for 5-membered rings */

	} else if (s_cmp(keyword, "TORSION5 ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    for (i__ = 1; i__ <= 6; ++i__) {
		vt[i__ - 1] = 0.;
		st[i__ - 1] = 0.;
		ft[i__ - 1] = 0;
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___105);
	    if (i__1 != 0) {
		goto L290;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L290;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L290;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L290;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L290;
	    }
	    for (j = 1; j <= 6; ++j) {
		i__1 = do_lio(&c__5, &c__1, (char *)&vt[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L290;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&st[j - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L290;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&ft[j - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L290;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L290;
	    }
L290:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    numeral_(&id, pd, &size, (ftnlen)4);
	    ++nt5;
	    if (ib < ic) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pa;
		i__4[1] = 4, a__3[1] = pb;
		i__4[2] = 4, a__3[2] = pc;
		i__4[3] = 4, a__3[3] = pd;
		s_cat(kt5_ref(0, nt5), a__3, i__4, &c__4, (ftnlen)16);
	    } else if (ic < ib) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pd;
		i__4[1] = 4, a__3[1] = pc;
		i__4[2] = 4, a__3[2] = pb;
		i__4[3] = 4, a__3[3] = pa;
		s_cat(kt5_ref(0, nt5), a__3, i__4, &c__4, (ftnlen)16);
	    } else if (ia <= id) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pa;
		i__4[1] = 4, a__3[1] = pb;
		i__4[2] = 4, a__3[2] = pc;
		i__4[3] = 4, a__3[3] = pd;
		s_cat(kt5_ref(0, nt5), a__3, i__4, &c__4, (ftnlen)16);
	    } else if (id < ia) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pd;
		i__4[1] = 4, a__3[1] = pc;
		i__4[2] = 4, a__3[2] = pb;
		i__4[3] = 4, a__3[3] = pa;
		s_cat(kt5_ref(0, nt5), a__3, i__4, &c__4, (ftnlen)16);
	    }
	    torphase_(ft, vt, st);
	    t15_ref(1, nt5) = vt[0];
	    t15_ref(2, nt5) = st[0];
	    t25_ref(1, nt5) = vt[1];
	    t25_ref(2, nt5) = st[1];
	    t35_ref(1, nt5) = vt[2];
	    t35_ref(2, nt5) = st[2];
	    t45_ref(1, nt5) = vt[3];
	    t45_ref(2, nt5) = st[3];
	    t55_ref(1, nt5) = vt[4];
	    t55_ref(2, nt5) = st[4];
	    t65_ref(1, nt5) = vt[5];
	    t65_ref(2, nt5) = st[5];

/*     torsional parameters for 4-membered rings */

	} else if (s_cmp(keyword, "TORSION4 ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    for (i__ = 1; i__ <= 6; ++i__) {
		vt[i__ - 1] = 0.;
		st[i__ - 1] = 0.;
		ft[i__ - 1] = 0;
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___106);
	    if (i__1 != 0) {
		goto L300;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L300;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L300;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L300;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L300;
	    }
	    for (i__ = 1; i__ <= 6; ++i__) {
		i__1 = do_lio(&c__5, &c__1, (char *)&vt[i__ - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L300;
		}
		i__1 = do_lio(&c__5, &c__1, (char *)&st[i__ - 1], (ftnlen)
			sizeof(doublereal));
		if (i__1 != 0) {
		    goto L300;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&ft[i__ - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L300;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L300;
	    }
L300:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    numeral_(&id, pd, &size, (ftnlen)4);
	    ++nt4;
	    if (ib < ic) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pa;
		i__4[1] = 4, a__3[1] = pb;
		i__4[2] = 4, a__3[2] = pc;
		i__4[3] = 4, a__3[3] = pd;
		s_cat(kt4_ref(0, nt4), a__3, i__4, &c__4, (ftnlen)16);
	    } else if (ic < ib) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pd;
		i__4[1] = 4, a__3[1] = pc;
		i__4[2] = 4, a__3[2] = pb;
		i__4[3] = 4, a__3[3] = pa;
		s_cat(kt4_ref(0, nt4), a__3, i__4, &c__4, (ftnlen)16);
	    } else if (ia <= id) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pa;
		i__4[1] = 4, a__3[1] = pb;
		i__4[2] = 4, a__3[2] = pc;
		i__4[3] = 4, a__3[3] = pd;
		s_cat(kt4_ref(0, nt4), a__3, i__4, &c__4, (ftnlen)16);
	    } else if (id < ia) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pd;
		i__4[1] = 4, a__3[1] = pc;
		i__4[2] = 4, a__3[2] = pb;
		i__4[3] = 4, a__3[3] = pa;
		s_cat(kt4_ref(0, nt4), a__3, i__4, &c__4, (ftnlen)16);
	    }
	    torphase_(ft, vt, st);
	    t14_ref(1, nt4) = vt[0];
	    t14_ref(2, nt4) = st[0];
	    t24_ref(1, nt4) = vt[1];
	    t24_ref(2, nt4) = st[1];
	    t34_ref(1, nt4) = vt[2];
	    t34_ref(2, nt4) = st[2];
	    t44_ref(1, nt4) = vt[3];
	    t44_ref(2, nt4) = st[3];
	    t54_ref(1, nt4) = vt[4];
	    t54_ref(2, nt4) = st[4];
	    t64_ref(1, nt4) = vt[5];
	    t64_ref(2, nt4) = st[5];

/*     pi-orbital torsion parameters */

	} else if (s_cmp(keyword, "PITORS ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    pt = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___108);
	    if (i__1 != 0) {
		goto L310;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L310;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L310;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pt, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L310;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L310;
	    }
L310:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++npt;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(kpt_ref(0, npt), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(kpt_ref(0, npt), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    kpitor_1.ptcon[npt - 1] = pt;

/*     stretch-torsion parameters */

	} else if (s_cmp(keyword, "STRTORS ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    bt1 = 0.;
	    bt2 = 0.;
	    bt3 = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___112);
	    if (i__1 != 0) {
		goto L320;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L320;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L320;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L320;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L320;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&bt1, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L320;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&bt2, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L320;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&bt3, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L320;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L320;
	    }
L320:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    numeral_(&id, pd, &size, (ftnlen)4);
	    ++nbt;
	    if (ib < ic) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pa;
		i__4[1] = 4, a__3[1] = pb;
		i__4[2] = 4, a__3[2] = pc;
		i__4[3] = 4, a__3[3] = pd;
		s_cat(kbt_ref(0, nbt), a__3, i__4, &c__4, (ftnlen)16);
	    } else if (ic < ib) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pd;
		i__4[1] = 4, a__3[1] = pc;
		i__4[2] = 4, a__3[2] = pb;
		i__4[3] = 4, a__3[3] = pa;
		s_cat(kbt_ref(0, nbt), a__3, i__4, &c__4, (ftnlen)16);
	    } else if (ia <= id) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pa;
		i__4[1] = 4, a__3[1] = pb;
		i__4[2] = 4, a__3[2] = pc;
		i__4[3] = 4, a__3[3] = pd;
		s_cat(kbt_ref(0, nbt), a__3, i__4, &c__4, (ftnlen)16);
	    } else if (id < ia) {
/* Writing concatenation */
		i__4[0] = 4, a__3[0] = pd;
		i__4[1] = 4, a__3[1] = pc;
		i__4[2] = 4, a__3[2] = pb;
		i__4[3] = 4, a__3[3] = pa;
		s_cat(kbt_ref(0, nbt), a__3, i__4, &c__4, (ftnlen)16);
	    }
	    btcon_ref(1, nbt) = bt1;
	    btcon_ref(2, nbt) = bt2;
	    btcon_ref(3, nbt) = bt3;

/*     torsion-torsion parameters */

	} else if (s_cmp(keyword, "TORTORS ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    ie = 0;
	    nx = 0;
	    ny = 0;
	    nxy = 0;
	    for (i__ = 1; i__ <= 900; ++i__) {
		tx[i__ - 1] = 0.;
		ty[i__ - 1] = 0.;
		tf[i__ - 1] = 0.;
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___120);
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
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L330;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ie, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L330;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&nx, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L330;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ny, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L330;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L330;
	    }
	    nxy = nx * ny;
	    i__1 = nxy;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		++iprm;
		s_copy(record, prmline_ref(0, iprm), (ftnlen)120, (ftnlen)120)
			;
		i__5 = s_rsli(&io___121);
		if (i__5 != 0) {
		    goto L330;
		}
		i__5 = do_lio(&c__5, &c__1, (char *)&tx[i__ - 1], (ftnlen)
			sizeof(doublereal));
		if (i__5 != 0) {
		    goto L330;
		}
		i__5 = do_lio(&c__5, &c__1, (char *)&ty[i__ - 1], (ftnlen)
			sizeof(doublereal));
		if (i__5 != 0) {
		    goto L330;
		}
		i__5 = do_lio(&c__5, &c__1, (char *)&tf[i__ - 1], (ftnlen)
			sizeof(doublereal));
		if (i__5 != 0) {
		    goto L330;
		}
		i__5 = e_rsli();
		if (i__5 != 0) {
		    goto L330;
		}
	    }
L330:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    numeral_(&id, pd, &size, (ftnlen)4);
	    numeral_(&ie, pe, &size, (ftnlen)4);
	    ++ntt;
/* Writing concatenation */
	    i__6[0] = 4, a__4[0] = pa;
	    i__6[1] = 4, a__4[1] = pb;
	    i__6[2] = 4, a__4[2] = pc;
	    i__6[3] = 4, a__4[3] = pd;
	    i__6[4] = 4, a__4[4] = pe;
	    s_cat(ktt_ref(0, ntt), a__4, i__6, &c__5, (ftnlen)20);
	    nx = nxy;
	    sort9_(&nx, tx);
	    ny = nxy;
	    sort9_(&ny, ty);
	    ktrtor_1.tnx[ntt - 1] = nx;
	    ktrtor_1.tny[ntt - 1] = ny;
	    i__1 = nx;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ttx_ref(i__, ntt) = tx[i__ - 1];
	    }
	    i__1 = ny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		tty_ref(i__, ntt) = ty[i__ - 1];
	    }
	    i__1 = nxy;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		tbf_ref(i__, ntt) = tf[i__ - 1];
	    }

/*     atomic partial charge parameters */

	} else if (s_cmp(keyword, "CHARGE ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    cg = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___124);
	    if (i__1 != 0) {
		goto L340;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L340;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&cg, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L340;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L340;
	    }
L340:
	    if (ia != 0) {
		kchrge_1.chg[ia - 1] = cg;
	    }

/*     bond dipole moment parameters */

	} else if (s_cmp(keyword, "DIPOLE ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    dp = 0.;
	    ps = .5;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___127);
	    if (i__1 != 0) {
		goto L350;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L350;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L350;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&dp, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L350;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ps, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L350;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L350;
	    }
L350:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++nd;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(kd_ref(0, nd), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(kd_ref(0, nd), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    kdipol_1.dpl[nd - 1] = dp;
	    kdipol_1.pos[nd - 1] = ps;

/*     bond dipole moment parameters for 5-membered rings */

	} else if (s_cmp(keyword, "DIPOLE5 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    dp = 0.;
	    ps = .5;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___128);
	    if (i__1 != 0) {
		goto L360;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L360;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L360;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&dp, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L360;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ps, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L360;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L360;
	    }
L360:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++nd5;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(kd5_ref(0, nd5), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(kd5_ref(0, nd5), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    kdipol_1.dpl5[nd5 - 1] = dp;
	    kdipol_1.pos5[nd5 - 1] = ps;

/*     bond dipole moment parameters for 4-membered rings */

	} else if (s_cmp(keyword, "DIPOLE4 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    dp = 0.;
	    ps = .5;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___129);
	    if (i__1 != 0) {
		goto L370;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L370;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L370;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&dp, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L370;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ps, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L370;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L370;
	    }
L370:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++nd4;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(kd4_ref(0, nd4), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(kd4_ref(0, nd4), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    kdipol_1.dpl4[nd4 - 1] = dp;
	    kdipol_1.pos4[nd4 - 1] = ps;

/*     bond dipole moment parameters for 3-membered rings */

	} else if (s_cmp(keyword, "DIPOLE3 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    dp = 0.;
	    ps = .5;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___130);
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
	    i__1 = do_lio(&c__5, &c__1, (char *)&dp, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L380;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ps, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L380;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L380;
	    }
L380:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++nd3;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(kd3_ref(0, nd3), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(kd3_ref(0, nd3), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    kdipol_1.dpl3[nd3 - 1] = dp;
	    kdipol_1.pos3[nd3 - 1] = ps;

/*     atomic multipole moment parameters */

	} else if (s_cmp(keyword, "MULTIPOLE ", (ftnlen)10, (ftnlen)10) == 0) 
		{
	    ia = 0;
	    ib = 0;
	    ic = 0;
	    id = 0;
	    s_copy(axt, "Z-then-X", (ftnlen)8, (ftnlen)8);
	    for (i__ = 1; i__ <= 13; ++i__) {
		pl[i__ - 1] = 0.;
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___133);
	    if (i__1 != 0) {
		goto L390;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L390;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L390;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L390;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&id, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L390;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl[0], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L390;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L390;
	    }
	    goto L400;
L390:
	    id = 0;
	    i__1 = s_rsli(&io___134);
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl[0], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L410;
	    }
L400:
	    ++iprm;
	    s_copy(record, prmline_ref(0, iprm), (ftnlen)120, (ftnlen)120);
	    i__1 = s_rsli(&io___135);
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl[1], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl[2], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl[3], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L410;
	    }
	    ++iprm;
	    s_copy(record, prmline_ref(0, iprm), (ftnlen)120, (ftnlen)120);
	    i__1 = s_rsli(&io___136);
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl[4], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L410;
	    }
	    ++iprm;
	    s_copy(record, prmline_ref(0, iprm), (ftnlen)120, (ftnlen)120);
	    i__1 = s_rsli(&io___137);
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl[7], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl[8], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L410;
	    }
	    ++iprm;
	    s_copy(record, prmline_ref(0, iprm), (ftnlen)120, (ftnlen)120);
	    i__1 = s_rsli(&io___138);
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl[10], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl[11], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pl[12], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L410;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L410;
	    }
L410:
	    if (ib < 0 || ic < 0) {
		s_copy(axt, "Bisector", (ftnlen)8, (ftnlen)8);
	    }
	    if (ic < 0 && id < 0) {
		s_copy(axt, "Z-Bisect", (ftnlen)8, (ftnlen)8);
	    }
/* Computing MAX */
	    i__1 = max(ib,ic);
	    if (max(i__1,id) < 0) {
		s_copy(axt, "3-Fold", (ftnlen)8, (ftnlen)6);
	    }
	    ib = abs(ib);
	    ic = abs(ic);
	    id = abs(id);
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    numeral_(&ic, pc, &size, (ftnlen)4);
	    numeral_(&id, pd, &size, (ftnlen)4);
	    ++nmp;
/* Writing concatenation */
	    i__4[0] = 4, a__3[0] = pa;
	    i__4[1] = 4, a__3[1] = pb;
	    i__4[2] = 4, a__3[2] = pc;
	    i__4[3] = 4, a__3[3] = pd;
	    s_cat(kmp_ref(0, nmp), a__3, i__4, &c__4, (ftnlen)16);
	    s_copy(mpaxis_ref(0, nmp), axt, (ftnlen)8, (ftnlen)8);
	    multip_ref(1, nmp) = pl[0];
	    multip_ref(2, nmp) = pl[1];
	    multip_ref(3, nmp) = pl[2];
	    multip_ref(4, nmp) = pl[3];
	    multip_ref(5, nmp) = pl[4];
	    multip_ref(6, nmp) = pl[7];
	    multip_ref(7, nmp) = pl[10];
	    multip_ref(8, nmp) = pl[7];
	    multip_ref(9, nmp) = pl[8];
	    multip_ref(10, nmp) = pl[11];
	    multip_ref(11, nmp) = pl[10];
	    multip_ref(12, nmp) = pl[11];
	    multip_ref(13, nmp) = pl[12];

/*     atomic dipole polarizability parameters */

	} else if (s_cmp(keyword, "POLARIZE ", (ftnlen)9, (ftnlen)9) == 0) {
	    ia = 0;
	    pol = 0.;
	    thl = 0.;
	    for (i__ = 1; i__ <= 8; ++i__) {
		pg[i__ - 1] = 0;
	    }
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___142);
	    if (i__1 != 0) {
		goto L420;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L420;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&pol, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L420;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&thl, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L420;
	    }
	    for (i__ = 1; i__ <= 8; ++i__) {
		i__1 = do_lio(&c__3, &c__1, (char *)&pg[i__ - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L420;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L420;
	    }
L420:
	    if (ia != 0) {
		kpolr_1.polr[ia - 1] = pol;
		kpolr_1.athl[ia - 1] = thl;
		for (i__ = 1; i__ <= 8; ++i__) {
		    pgrp_ref(i__, ia) = pg[i__ - 1];
		}
	    }

/*     conjugated pisystem atom parameters */

	} else if (s_cmp(keyword, "PIATOM ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    el = 0.;
	    iz = 0.;
	    rp = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___146);
	    if (i__1 != 0) {
		goto L430;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L430;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&el, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L430;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&iz, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L430;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&rp, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L430;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L430;
	    }
L430:
	    if (ia != 0) {
		korbs_1.electron[ia - 1] = el;
		korbs_1.ionize[ia - 1] = iz;
		korbs_1.repulse[ia - 1] = rp;
	    }

/*     conjugated pisystem bond parameters */

	} else if (s_cmp(keyword, "PIBOND ", (ftnlen)7, (ftnlen)7) == 0) {
	    ia = 0;
	    ib = 0;
	    ss = 0.;
	    ts = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___149);
	    if (i__1 != 0) {
		goto L440;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L440;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L440;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ss, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L440;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ts, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L440;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L440;
	    }
L440:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++npi;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(kpi_ref(0, npi), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(kpi_ref(0, npi), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    korbs_1.sslope[npi - 1] = ss;
	    korbs_1.tslope[npi - 1] = ts;

/*     conjugated pisystem bond parameters for 5-membered rings */

	} else if (s_cmp(keyword, "PIBOND5 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ss = 0.;
	    ts = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___150);
	    if (i__1 != 0) {
		goto L450;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L450;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L450;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ss, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L450;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ts, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L450;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L450;
	    }
L450:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++npi5;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(kpi5_ref(0, npi5), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(kpi5_ref(0, npi5), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    korbs_1.sslope5[npi5 - 1] = ss;
	    korbs_1.tslope5[npi5 - 1] = ts;

/*     conjugated pisystem bond parameters for 4-membered rings */

	} else if (s_cmp(keyword, "PIBOND4 ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    ss = 0.;
	    ts = 0.;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___151);
	    if (i__1 != 0) {
		goto L460;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L460;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L460;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ss, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L460;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&ts, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L460;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L460;
	    }
L460:
	    numeral_(&ia, pa, &size, (ftnlen)4);
	    numeral_(&ib, pb, &size, (ftnlen)4);
	    ++npi4;
	    if (ia <= ib) {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pa;
		i__2[1] = 4, a__1[1] = pb;
		s_cat(kpi4_ref(0, npi4), a__1, i__2, &c__2, (ftnlen)8);
	    } else {
/* Writing concatenation */
		i__2[0] = 4, a__1[0] = pb;
		i__2[1] = 4, a__1[1] = pa;
		s_cat(kpi4_ref(0, npi4), a__1, i__2, &c__2, (ftnlen)8);
	    }
	    korbs_1.sslope4[npi4 - 1] = ss;
	    korbs_1.tslope4[npi4 - 1] = ts;

/*     metal ligand field splitting parameters */

	} else if (s_cmp(keyword, "METAL ", (ftnlen)6, (ftnlen)6) == 0) {
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___152);
	    if (i__1 != 0) {
		goto L470;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L470;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L470;
	    }
L470:

/*     biopolymer atom type conversion definitions */

	    ;
	} else if (s_cmp(keyword, "BIOTYPE ", (ftnlen)8, (ftnlen)8) == 0) {
	    ia = 0;
	    ib = 0;
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___153);
	    if (i__1 != 0) {
		goto L480;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L480;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L480;
	    }
	    getword_(record, string, &next, (ftnlen)120, (ftnlen)120);
	    getstring_(record, string, &next, (ftnlen)120, (ftnlen)120);
	    s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1))
		    ;
	    i__1 = s_rsli(&io___154);
	    if (i__1 != 0) {
		goto L480;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    if (i__1 != 0) {
		goto L480;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L480;
	    }
L480:
	    if (ia >= 10000) {
		io___155.ciunit = iounit_1.iout;
		s_wsfe(&io___155);
		e_wsfe();
/* L490: */
		fatal_();
	    }
	    if (ia != 0) {
		fields_1.biotyp[ia - 1] = ib;
	    }
	}
    }
    return 0;
} /* readprm_ */

#undef prmline_ref
#undef mpaxis_ref
#undef multip_ref
#undef symbol_ref
#undef btcon_ref
#undef kvpr_ref
#undef pgrp_ref
#undef stbn_ref
#undef kopd_ref
#undef kopb_ref
#undef anan_ref
#undef angf_ref
#undef kpi5_ref
#undef kpi4_ref
#undef ang5_ref
#undef ang4_ref
#undef ang3_ref
#undef tty_ref
#undef ttx_ref
#undef ktt_ref
#undef kpt_ref
#undef kmp_ref
#undef kti_ref
#undef kpi_ref
#undef kbt_ref
#undef ksb_ref
#undef tbf_ref
#undef kel_ref
#undef kdi_ref
#undef ang_ref
#undef khb_ref
#undef kaf_ref
#undef kt5_ref
#undef kt4_ref
#undef ti3_ref
#undef ti2_ref
#undef ti1_ref
#undef kd3_ref
#undef kd4_ref
#undef kd5_ref
#undef kb3_ref
#undef kb4_ref
#undef kb5_ref
#undef ka5_ref
#undef ka4_ref
#undef ka3_ref
#undef ku_ref
#undef kt_ref
#undef t64_ref
#undef t54_ref
#undef t44_ref
#undef t34_ref
#undef t24_ref
#undef t65_ref
#undef t55_ref
#undef t45_ref
#undef t35_ref
#undef t25_ref
#undef t15_ref
#undef t14_ref
#undef kd_ref
#undef kb_ref
#undef ka_ref
#undef t6_ref
#undef t5_ref
#undef t4_ref
#undef t3_ref
#undef t2_ref
#undef t1_ref
#undef describe_ref


