/* prtprm.f -- translated by f2c (version 20050501).
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
    doublereal electric, dielec, ebuffer, c2scale, c3scale, c4scale, c5scale;
    logical neutnbr, neutcut;
} chgpot_;

#define chgpot_1 chgpot_

struct {
    integer biotyp[10000];
    char forcefield[20];
} fields_;

#define fields_1 fields_

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
    doublereal m2scale, m3scale, m4scale, m5scale;
} mplpot_;

#define mplpot_1 mplpot_

struct {
    doublereal poleps, polsor, p2scale, p3scale, p4scale, p5scale, d1scale, 
	    d2scale, d3scale, d4scale, u1scale, u2scale, u3scale, u4scale;
    char poltyp[6];
} polpot_;

#define polpot_1 polpot_

struct {
    doublereal cury, qury, ureyunit;
} urypot_;

#define urypot_1 urypot_

struct {
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
    logical use_vcorr__;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_

/* Table of constant values */

static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine prtprm  --  output of force field parameters  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "prtprm" writes out a formatted listing of the default */
/*     set of potential energy parameters for a force field */


/* Subroutine */ int prtprm_(integer *itxt)
{
    /* Format strings */
    static char fmt_10[] = "(//,15x,\002TINKER Force Field Parameters for"
	    " \002,a20)";
    static char fmt_20[] = "(//,15x,\002Force Field Atom Definitions\002,//,"
	    "54x,\002Atomic\002,4x,\002Atomic\002,/,5x,\002Type\002,3x,\002Cl"
	    "ass\002,3x,\002Symbol\002,3x,\002Description\002,14x,\002Numbe"
	    "r\002,4x,\002Weight\002,3x,\002Valence\002,/)";
    static char fmt_30[] = "(3x,i5,3x,i5,5x,a3,5x,a24,i5,f12.3,i7)";
    static char fmt_40[] = "(a1,//,15x,\002TINKER Force Field Parameters f"
	    "or \002,a20)";
    static char fmt_50[] = "(//,15x,\002Van der Waals Parameters\002,///,2"
	    "0x,\002Class\002,7x,\002Radius\002,6x,\002Epsilon\002,4x,\002Red"
	    "uction\002,/)";
    static char fmt_60[] = "(//,15x,\002Van der Waals Parameters\002,///,2"
	    "0x,\002Type\002,8x,\002Radius\002,6x,\002Epsilon\002,4x,\002Redu"
	    "ction\002,/)";
    static char fmt_70[] = "(10x,i5,5x,i4,2x,3f12.3)";
    static char fmt_80[] = "(//,15x,\002Van der Waals Scaling Factors\002,//"
	    "/,20x,\0021-2 Atoms\002,f17.3,/,20x,\0021-3 Atoms\002,f17.3,/,20"
	    "x,\0021-4 Atoms\002,f17.3,/,20x,\0021-5 Atoms\002,f17.3)";
    static char fmt_90[] = "(//,15x,\002Van der Waals Parameters for 1-"
	    "4\002,\002 Interactions\002,///,20x,\002Class\002,7x,\002Radiu"
	    "s\002,6x,\002Epsilon\002,/)";
    static char fmt_100[] = "(//,15x,\002Van der Waals Parameters for 1-4"
	    "\002,\002 Interactions\002,///,20x,\002Type\002,8x,\002Radius"
	    "\002,6x,\002Epsilon\002,/)";
    static char fmt_110[] = "(10x,i5,5x,i4,2x,2f12.3)";
    static char fmt_120[] = "(//,15x,\002Van der Waals Parameters for Atom P"
	    "airs\002,///,22x,\002Classes\002,7x,\002Radii Sum\002,4x,\002Eps"
	    "ilon\002,/)";
    static char fmt_130[] = "(//,15x,\002Van der Waals Parameters for Atom P"
	    "airs\002,///,23x,\002Types\002,8x,\002Radii Sum\002,4x,\002Epsil"
	    "on\002,/)";
    static char fmt_140[] = "(10x,i5,5x,i4,\002-\002,i4,2x,2f12.3)";
    static char fmt_160[] = "(//,15x,\002Hydrogen Bonding Parameters for Ato"
	    "m Pairs\002,///,22x,\002Classes\002,7x,\002Radii Sum\002,4x,\002"
	    "Epsilon\002,/)";
    static char fmt_170[] = "(//,15x,\002Hydrogen Bonding Parameters for Ato"
	    "m Pairs\002,///,23x,\002Types\002,8x,\002Radii Sum\002,4x,\002Ep"
	    "silon\002,/)";
    static char fmt_180[] = "(10x,i5,5x,i4,\002-\002,i4,2x,2f12.3)";
    static char fmt_200[] = "(a1,//,15x,\002TINKER Force Field Parameters fo"
	    "r \002,a20)";
    static char fmt_210[] = "(//,15x,\002Bond Stretching Parameters\002,///,"
	    "22x,\002Classes\002,15x,\002KS\002,7x,\002Length\002,/)";
    static char fmt_220[] = "(10x,i5,5x,i4,\002-\002,i4,6x,f12.3,f12.4)";
    static char fmt_240[] = "(//,15x,\0025-Membered Ring Stretch Parameter"
	    "s\002,///,22x,\002Classes\002,15x,\002KS\002,7x,\002Length\002,/)"
	    ;
    static char fmt_250[] = "(10x,i5,5x,i4,\002-\002,i4,6x,f12.3,f12.4)";
    static char fmt_270[] = "(//,15x,\0024-Membered Ring Stretch Parameter"
	    "s\002,///,22x,\002Classes\002,15x,\002KS\002,7x,\002Length\002,/)"
	    ;
    static char fmt_310[] = "(10x,i5,5x,i4,\002-\002,i4,6x,f12.3,f12.4)";
    static char fmt_300[] = "(//,15x,\0023-Membered Ring Stretch Parameter"
	    "s\002,///,22x,\002Classes\002,15x,\002KS\002,7x,\002Length\002,/)"
	    ;
    static char fmt_330[] = "(//,15x,\002Higher Order Stretching Constant"
	    "s\002,///,20x,\002Cubic\002,f17.3,/,20x,\002Quartic\002,f15.3)";
    static char fmt_340[] = "(//,15x,\002Electronegativity Bond Length Param"
	    "eters\002,///,25x,\002Classes\002,21x,\002dLength\002,/)";
    static char fmt_350[] = "(10x,i5,5x,i4,\002-\002,i4,\002-\002,i4,14x,f12"
	    ".4)";
    static char fmt_370[] = "(a1,//,15x,\002TINKER Force Field Parameters fo"
	    "r \002,a20)";
    static char fmt_380[] = "(//,15x,\002Angle Bending Parameters\002,///,18"
	    "x,\002Classes\002,11x,\002KB\002,6x,\002Value 1\002,5x,\002Value"
	    " 2\002,5x,\002Value 3\002,/,44x,\002(R-X-R)\002,5x,\002(R-X-H"
	    ")\002,5x,\002(H-X-H)\002,/)";
    static char fmt_390[] = "(3x,i5,5x,i4,\002-\002,i4,\002-\002,i4,2f12.3)";
    static char fmt_400[] = "(3x,i5,5x,i4,\002-\002,i4,\002-\002,i4,4f12.3)";
    static char fmt_420[] = "(//,17x,\0025-Membered Ring Bend Parameters\002"
	    ",///,18x,\002Classes\002,11x,\002KB\002,6x,\002Value 1\002,5x"
	    ",\002Value 2\002,5x,\002Value 3\002,/,44x,\002(R-X-R)\002,5x,"
	    "\002(R-X-H)\002,5x,\002(H-X-H)\002,/)";
    static char fmt_430[] = "(3x,i5,5x,i4,\002-\002,i4,\002-\002,i4,2f12.3)";
    static char fmt_440[] = "(3x,i5,5x,i4,\002-\002,i4,\002-\002,i4,4f12.3)";
    static char fmt_460[] = "(//,15x,\0024-Membered Ring Bend Parameters\002"
	    ",///,18x,\002Classes\002,11x,\002KB\002,6x,\002Value 1\002,5x"
	    ",\002Value 2\002,5x,\002Value 3\002,/,44x,\002(R-X-R)\002,5x,"
	    "\002(R-X-H)\002,5x,\002(H-X-H)\002,/)";
    static char fmt_470[] = "(3x,i5,5x,i4,\002-\002,i4,\002-\002,i4,2f12.3)";
    static char fmt_480[] = "(3x,i5,5x,i4,\002-\002,i4,\002-\002,i4,4f12.3)";
    static char fmt_500[] = "(//,15x,\0023-Membered Ring Bend Parameters\002"
	    ",///,18x,\002Classes\002,11x,\002KB\002,6x,\002Value 1\002,5x"
	    ",\002Value 2\002,5x,\002Value 3\002,/,44x,\002(R-X-R)\002,5x,"
	    "\002(R-X-H)\002,5x,\002(H-X-H)\002,/)";
    static char fmt_510[] = "(3x,i5,5x,i4,\002-\002,i4,\002-\002,i4,2f12.3)";
    static char fmt_520[] = "(3x,i5,5x,i4,\002-\002,i4,\002-\002,i4,4f12.3)";
    static char fmt_540[] = "(//,15x,\002Fourier Angle Bending Parameters"
	    "\002,///,18x,\002Classes\002,11x,\002KB\002,8x,\002Shift\002,6x"
	    ",\002Period\002,/)";
    static char fmt_550[] = "(3x,i5,5x,i4,\002-\002,i4,\002-\002,i4,3f12.3)";
    static char fmt_570[] = "(//,15x,\002Higher Order Bending Constants\002,"
	    "///,20x,\002Cubic\002,d17.3,/,20x,\002Quartic\002,d15.3,/,20x"
	    ",\002Pentic\002,d16.3,/,20x,\002Sextic\002,d16.3)";
    static char fmt_580[] = "(a1,//,15x,\002TINKER Force Field Parameters fo"
	    "r \002,a20)";
    static char fmt_590[] = "(//,15x,\002Stretch-Bend Parameters\002,///,1"
	    "8x,\002Classes\002,10x,\002KSB1\002,8x,\002KSB2\002,/)";
    static char fmt_600[] = "(3x,i5,5x,i4,\002-\002,i4,\002-\002,i4,2f12.3)";
    static char fmt_620[] = "(a1,//,15x,\002TINKER Force Field Parameters fo"
	    "r \002,a20)";
    static char fmt_630[] = "(//,15x,\002Urey-Bradley Parameters\002,///,1"
	    "8x,\002Classes\002,11x,\002KB\002,6x,\002Distance\002,/)";
    static char fmt_640[] = "(3x,i5,5x,i4,\002-\002,i4,\002-\002,i4,f12.3,f1"
	    "2.4)";
    static char fmt_660[] = "(//,15x,\002Higher Order Urey-Bradley Constant"
	    "s\002,///,20x,\002Cubic\002,f17.3,/,20x,\002Quartic\002,f15.3)";
    static char fmt_670[] = "(a1,//,15x,\002TINKER Force Field Parameters fo"
	    "r \002,a20)";
    static char fmt_680[] = "(//,15x,\002Angle-Angle Parameters\002,///,20x"
	    ",\002Class\002,9x,\002KAA 1\002,7x,\002KAA 2\002,7x,\002KAA 3"
	    "\002,/,33x,\002(R-X-R)\002,5x,\002(R-X-H)\002,5x,\002(H-X-H)\002"
	    ",/)";
    static char fmt_690[] = "(8x,i5,7x,i4,3x,3f12.3)";
    static char fmt_700[] = "(a1,//,15x,\002TINKER Force Field Parameters fo"
	    "r \002,a20)";
    static char fmt_710[] = "(//,15x,\002Out-of-Plane Bend Parameters\002,//"
	    "/,26x,\002Classes\002,11x,\002KOPB\002,/)";
    static char fmt_720[] = "(8x,i5,5x,i4,\002-\002,i4,\002-\002,i4,\002-"
	    "\002,i4,f12.3)";
    static char fmt_740[] = "(a1,//,15x,\002TINKER Force Field Parameters fo"
	    "r \002,a20)";
    static char fmt_750[] = "(//,15x,\002Out-of-Plane Distance Parameters"
	    "\002,///,26x,\002Classes\002,11x,\002KOPD\002,/)";
    static char fmt_760[] = "(8x,i5,5x,i4,\002-\002,i4,\002-\002,i4,\002-"
	    "\002,i4,f12.3)";
    static char fmt_780[] = "(a1,//,15x,\002TINKER Force Field Parameters fo"
	    "r \002,a20)";
    static char fmt_790[] = "(//,15x,\002Improper Dihedral Parameters\002,//"
	    "/,20x,\002Classes\002,12x,\002KID\002,7x,\002Target\002,/)";
    static char fmt_800[] = "(2x,i5,5x,i4,\002-\002,i4,\002-\002,i4,\002-"
	    "\002,i4,f12.3,f12.4)";
    static char fmt_820[] = "(a1,//,15x,\002TINKER Force Field Parameters fo"
	    "r \002,a20)";
    static char fmt_830[] = "(//,15x,\002Improper Torsion Parameters\002,///"
	    ",17x,\002Classes\002,15x,\002KTI Values\002,/)";
    static char fmt_840[] = "(2x,i5,2x,i4,\002-\002,i4,\002-\002,i4,\002-"
	    "\002,i4,2x,3(f8.3,f6.1,i2))";
    static char fmt_860[] = "(a1,//,15x,\002TINKER Force Field Parameters fo"
	    "r \002,a20)";
    static char fmt_870[] = "(//,15x,\002Torsional Parameters\002,///,17x"
	    ",\002Classes\002,15x,\002KT Values\002,/)";
    static char fmt_880[] = "(2x,i5,2x,i4,\002-\002,i4,\002-\002,i4,\002-"
	    "\002,i4,2x,6(f8.3,f6.1,i2))";
    static char fmt_900[] = "(//,15x,\0025-Membered Ring Torsion Parameter"
	    "s\002,///,17x,\002Classes\002,15x,\002KT Values\002,/)";
    static char fmt_910[] = "(2x,i5,2x,i4,\002-\002,i4,\002-\002,i4,\002-"
	    "\002,i4,2x,6(f8.3,f6.1,i2))";
    static char fmt_930[] = "(//,15x,\0024-Membered Ring Torsion Parameter"
	    "s\002,///,17x,\002Classes\002,15x,\002KT Values\002,/)";
    static char fmt_940[] = "(2x,i5,2x,i4,\002-\002,i4,\002-\002,i4,\002-"
	    "\002,i4,2x,6(f8.3,f6.1,i2))";
    static char fmt_960[] = "(a1,//,15x,\002TINKER Force Field Parameters fo"
	    "r \002,a20)";
    static char fmt_970[] = "(//,15x,\002Pi-Orbital Torsion Parameters\002,/"
	    "//,18x,\002Classes\002,15x,\002KPT\002,/)";
    static char fmt_980[] = "(6x,i5,5x,i4,\002-\002,i4,6x,f12.3)";
    static char fmt_1000[] = "(a1,//,15x,\002TINKER Force Field Parameters f"
	    "or \002,a20)";
    static char fmt_1010[] = "(//,15x,\002Stretch-Torsion Parameters\002,///"
	    ",17x,\002Classes\002,15x,\002KST1\002,8x,\002KST2\002,8x,\002KST3"
	    "\002,/)";
    static char fmt_1020[] = "(2x,i5,2x,i4,\002-\002,i4,\002-\002,i4,\002"
	    "-\002,i4,3x,3f12.3)";
    static char fmt_1040[] = "(a1,//,15x,\002TINKER Force Field Parameters f"
	    "or \002,a20)";
    static char fmt_1050[] = "(//,15x,\002Torsion-Torsion Parameters\002,///"
	    ",19x,\002Classes\002,18x,\002KNX\002,9x,\002KNY\002)";
    static char fmt_1060[] = "(/,2x,i5,2x,i4,\002-\002,i4,\002-\002,i4,\002"
	    "-\002,i4,\002-\002,i4,2x,2i12,/)";
    static char fmt_1070[] = "(3x,6f12.4)";
    static char fmt_1090[] = "(a1,//,15x,\002TINKER Force Field Parameters f"
	    "or \002,a20)";
    static char fmt_1100[] = "(//,15x,\002Atomic Partial Charge Parameter"
	    "s\002,///,24x,\002Type\002,9x,\002Partial Chg\002,/)";
    static char fmt_1110[] = "(12x,i5,7x,i3,6x,f12.3)";
    static char fmt_1120[] = "(//,15x,\002Atomic Partial Charge Scaling Fact"
	    "ors\002,///,20x,\0021-2 Atoms\002,f17.3,/,20x,\0021-3 Atoms\002,"
	    "f17.3,/,20x,\0021-4 Atoms\002,f17.3,/,20x,\0021-5 Atoms\002,f17."
	    "3)";
    static char fmt_1130[] = "(a1,//,15x,\002TINKER Force Field Parameters f"
	    "or \002,a20)";
    static char fmt_1140[] = "(//,15x,\002Bond Dipole Moment Parameters\002,"
	    "///,25x,\002Types\002,10x,\002Bond Dipole\002,4x,\002Position"
	    "\002,/)";
    static char fmt_1150[] = "(12x,i5,5x,i4,\002-\002,i4,6x,2f12.3)";
    static char fmt_1170[] = "(//,15x,\0025-Membered Ring Bond Dipole Parame"
	    "ters\002,///,25x,\002Types\002,10x,\002Bond Dipole\002,4x,\002Po"
	    "sition\002,/)";
    static char fmt_1180[] = "(12x,i5,5x,i4,\002-\002,i4,6x,2f12.3)";
    static char fmt_1200[] = "(//,15x,\0024-Membered Ring Bond Dipole Parame"
	    "ters\002,///,25x,\002Types\002,10x,\002Bond Dipole\002,4x,\002Po"
	    "sition\002,/)";
    static char fmt_1210[] = "(12x,i5,5x,i4,\002-\002,i4,6x,2f12.3)";
    static char fmt_1230[] = "(//,15x,\0023-Membered Ring Bond Dipole Parame"
	    "ters\002,///,25x,\002Types\002,10x,\002Bond Dipole\002,4x,\002Po"
	    "sition\002,/)";
    static char fmt_1240[] = "(12x,i5,5x,i4,\002-\002,i4,6x,2f12.3)";
    static char fmt_1260[] = "(a1,//,15x,\002TINKER Force Field Parameters f"
	    "or \002,a20)";
    static char fmt_1270[] = "(//,17x,\002Atomic Multipole Parameters\002,//"
	    "/,11x,\002Type\002,7x,\002Axis Types\002,8x,\002Frame\002,9x,"
	    "\002Multipoles (M-D-Q)\002,/)";
    static char fmt_1280[] = "(2x,i5,3x,i4,3x,i4,2x,i4,2x,i4,5x,a8,2x,f10.5,"
	    "/,48x,3f10.5,/,48x,f10.5,/,48x,2f10.5,/,48x,3f10.5)";
    static char fmt_1300[] = "(//,15x,\002Atomic Multipole Scaling Factor"
	    "s\002,///,20x,\0021-2 Atoms\002,f17.3,/,20x,\0021-3 Atoms\002,f1"
	    "7.3,/,20x,\0021-4 Atoms\002,f17.3,/,20x,\0021-5 Atoms\002,f17.3)";
    static char fmt_1310[] = "(a1,//,15x,\002TINKER Force Field Parameters f"
	    "or \002,a20)";
    static char fmt_1320[] = "(//,15x,\002Dipole Polarizability Parameter"
	    "s\002,///,23x,\002Type\002,7x,\002Alpha\002,6x,\002Damp\002,6x"
	    ",\002Group Atom Types\002,/)";
    static char fmt_1330[] = "(10x,i5,7x,i4,3x,2f10.3)";
    static char fmt_1340[] = "(10x,i5,7x,i4,3x,2f10.3,4x,6i5)";
    static char fmt_1350[] = "(//,15x,\002Direct Polarizability Scaling Fact"
	    "ors\002,///,20x,\0021-1 Groups\002,f15.3,/,20x,\0021-2 Groups"
	    "\002,f15.3,/,20x,\0021-3 Groups\002,f15.3,/,20x,\0021-4 Group"
	    "s\002,f15.3)";
    static char fmt_1360[] = "(//,15x,\002Mutual Polarizability Scaling Fact"
	    "ors\002,///,20x,\0021-1 Groups\002,f15.3,/,20x,\0021-2 Groups"
	    "\002,f15.3,/,20x,\0021-3 Groups\002,f15.3,/,20x,\0021-4 Group"
	    "s\002,f15.3)";
    static char fmt_1370[] = "(//,15x,\002Polarizability Energy Scaling Fact"
	    "ors\002,///,20x,\0021-2 Atoms\002,f16.3,/,20x,\0021-3 Atoms\002,"
	    "f16.3,/,20x,\0021-4 Atoms\002,f16.3,/,20x,\0021-5 Atoms\002,f16."
	    "3)";
    static char fmt_1380[] = "(a1,//,15x,\002TINKER Force Field Parameters f"
	    "or \002,a20)";
    static char fmt_1390[] = "(//,15x,\002Conjugated Pisystem Atom Paramet"
	    "ers\002,///,20x,\002Class\002,3x,\002Electron\002,3x,\002Ionizat"
	    "ion\002,3x,\002Repulsion\002,/)";
    static char fmt_1400[] = "(8x,i5,7x,i4,f10.1,2x,2f12.3)";
    static char fmt_1410[] = "(//,15x,\002Conjugated Pisystem Bond Paramet"
	    "ers\002,///,20x,\002Classes\002,8x,\002d Force\002,4x,\002d Leng"
	    "th\002,/)";
    static char fmt_1420[] = "(8x,i5,5x,i4,\002-\002,i4,3x,f12.3,f12.3)";
    static char fmt_1440[] = "(//,15x,\0025-Membered Ring Pisystem Bond Para"
	    "meters\002,///,20x,\002Classes\002,8x,\002d Force\002,4x,\002d L"
	    "ength\002,/)";
    static char fmt_1450[] = "(8x,i5,5x,i4,\002-\002,i4,3x,f12.3,f12.3)";
    static char fmt_1470[] = "(//,15x,\0024-Membered Ring Pisystem Bond Para"
	    "meters\002,///,20x,\002Classes\002,8x,\002d Force\002,4x,\002d L"
	    "ength\002,/)";
    static char fmt_1480[] = "(8x,i5,5x,i4,\002-\002,i4,3x,f12.3,f12.3)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), do_fio(
	    integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static char formfeed[1];
    static integer i__, j, k, k1, k2, k3, k4, k5, npg, fold[6];
    static doublereal phase[6], ampli[6];
    static logical exist;
    static char blank3[3], blank8[8], blank20[20], blank12[12], blank16[16];
    extern integer number_(char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___9 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_330, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_340, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_350, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_370, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_380, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_390, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_400, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_420, 0 };
    static cilist io___48 = { 0, 0, 0, fmt_430, 0 };
    static cilist io___49 = { 0, 0, 0, fmt_440, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_460, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_470, 0 };
    static cilist io___52 = { 0, 0, 0, fmt_480, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_500, 0 };
    static cilist io___54 = { 0, 0, 0, fmt_510, 0 };
    static cilist io___55 = { 0, 0, 0, fmt_520, 0 };
    static cilist io___56 = { 0, 0, 0, fmt_540, 0 };
    static cilist io___57 = { 0, 0, 0, fmt_550, 0 };
    static cilist io___58 = { 0, 0, 0, fmt_570, 0 };
    static cilist io___59 = { 0, 0, 0, fmt_580, 0 };
    static cilist io___60 = { 0, 0, 0, fmt_590, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___64 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___65 = { 0, 0, 0, fmt_660, 0 };
    static cilist io___66 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_680, 0 };
    static cilist io___68 = { 0, 0, 0, fmt_690, 0 };
    static cilist io___69 = { 0, 0, 0, fmt_700, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_710, 0 };
    static cilist io___72 = { 0, 0, 0, fmt_720, 0 };
    static cilist io___73 = { 0, 0, 0, fmt_740, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_750, 0 };
    static cilist io___75 = { 0, 0, 0, fmt_760, 0 };
    static cilist io___76 = { 0, 0, 0, fmt_780, 0 };
    static cilist io___77 = { 0, 0, 0, fmt_790, 0 };
    static cilist io___78 = { 0, 0, 0, fmt_800, 0 };
    static cilist io___79 = { 0, 0, 0, fmt_820, 0 };
    static cilist io___80 = { 0, 0, 0, fmt_830, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_840, 0 };
    static cilist io___85 = { 0, 0, 0, fmt_860, 0 };
    static cilist io___86 = { 0, 0, 0, fmt_870, 0 };
    static cilist io___87 = { 0, 0, 0, fmt_880, 0 };
    static cilist io___88 = { 0, 0, 0, fmt_900, 0 };
    static cilist io___89 = { 0, 0, 0, fmt_910, 0 };
    static cilist io___90 = { 0, 0, 0, fmt_930, 0 };
    static cilist io___91 = { 0, 0, 0, fmt_940, 0 };
    static cilist io___92 = { 0, 0, 0, fmt_960, 0 };
    static cilist io___93 = { 0, 0, 0, fmt_970, 0 };
    static cilist io___94 = { 0, 0, 0, fmt_980, 0 };
    static cilist io___95 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___96 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___97 = { 0, 0, 0, fmt_1020, 0 };
    static cilist io___98 = { 0, 0, 0, fmt_1040, 0 };
    static cilist io___99 = { 0, 0, 0, fmt_1050, 0 };
    static cilist io___101 = { 0, 0, 0, fmt_1060, 0 };
    static cilist io___102 = { 0, 0, 0, fmt_1070, 0 };
    static cilist io___103 = { 0, 0, 0, fmt_1090, 0 };
    static cilist io___104 = { 0, 0, 0, fmt_1100, 0 };
    static cilist io___105 = { 0, 0, 0, fmt_1110, 0 };
    static cilist io___106 = { 0, 0, 0, fmt_1120, 0 };
    static cilist io___107 = { 0, 0, 0, fmt_1130, 0 };
    static cilist io___108 = { 0, 0, 0, fmt_1140, 0 };
    static cilist io___109 = { 0, 0, 0, fmt_1150, 0 };
    static cilist io___110 = { 0, 0, 0, fmt_1170, 0 };
    static cilist io___111 = { 0, 0, 0, fmt_1180, 0 };
    static cilist io___112 = { 0, 0, 0, fmt_1200, 0 };
    static cilist io___113 = { 0, 0, 0, fmt_1210, 0 };
    static cilist io___114 = { 0, 0, 0, fmt_1230, 0 };
    static cilist io___115 = { 0, 0, 0, fmt_1240, 0 };
    static cilist io___116 = { 0, 0, 0, fmt_1260, 0 };
    static cilist io___117 = { 0, 0, 0, fmt_1270, 0 };
    static cilist io___118 = { 0, 0, 0, fmt_1280, 0 };
    static cilist io___119 = { 0, 0, 0, fmt_1300, 0 };
    static cilist io___120 = { 0, 0, 0, fmt_1310, 0 };
    static cilist io___121 = { 0, 0, 0, fmt_1320, 0 };
    static cilist io___123 = { 0, 0, 0, fmt_1330, 0 };
    static cilist io___124 = { 0, 0, 0, fmt_1340, 0 };
    static cilist io___125 = { 0, 0, 0, fmt_1350, 0 };
    static cilist io___126 = { 0, 0, 0, fmt_1360, 0 };
    static cilist io___127 = { 0, 0, 0, fmt_1370, 0 };
    static cilist io___128 = { 0, 0, 0, fmt_1380, 0 };
    static cilist io___129 = { 0, 0, 0, fmt_1390, 0 };
    static cilist io___130 = { 0, 0, 0, fmt_1400, 0 };
    static cilist io___131 = { 0, 0, 0, fmt_1410, 0 };
    static cilist io___132 = { 0, 0, 0, fmt_1420, 0 };
    static cilist io___133 = { 0, 0, 0, fmt_1440, 0 };
    static cilist io___134 = { 0, 0, 0, fmt_1450, 0 };
    static cilist io___135 = { 0, 0, 0, fmt_1470, 0 };
    static cilist io___136 = { 0, 0, 0, fmt_1480, 0 };



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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  chgpot.i  --  specifics of charge-charge functional form  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     electric   energy factor in kcal/mole for current force field */
/*     dielec     dielectric constant for electrostatic interactions */
/*     ebuffer    electrostatic buffering constant added to distance */
/*     c2scale    factor by which 1-2 charge interactions are scaled */
/*     c3scale    factor by which 1-3 charge interactions are scaled */
/*     c4scale    factor by which 1-4 charge interactions are scaled */
/*     c5scale    factor by which 1-5 charge interactions are scaled */
/*     neutnbr    logical flag governing use of neutral group neighbors */
/*     neutcut    logical flag governing use of neutral group cutoffs */




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
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  mplpot.i  --  specifics of atomic multipole functions  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     m2scale   factor by which 1-2 multipole interactions are scaled */
/*     m3scale   factor by which 1-3 multipole interactions are scaled */
/*     m4scale   factor by which 1-4 multipole interactions are scaled */
/*     m5scale   factor by which 1-5 multipole interactions are scaled */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  polpot.i  --  specifics of polarization functional form  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     poleps    induced dipole convergence criterion (rms Debyes/atom) */
/*     polsor    induced dipole SOR convergence acceleration factor */
/*     p2scale   field 1-2 scale factor for energy evaluations */
/*     p3scale   field 1-3 scale factor for energy evaluations */
/*     p4scale   field 1-4 scale factor for energy evaluations */
/*     p5scale   field 1-5 scale factor for energy evaluations */
/*     d1scale   field intra-group scale factor for direct induced */
/*     d2scale   field 1-2 group scale factor for direct induced */
/*     d3scale   field 1-3 group scale factor for direct induced */
/*     d4scale   field 1-4 group scale factor for direct induced */
/*     u1scale   field intra-group scale factor for mutual induced */
/*     u2scale   field 1-2 group scale factor for mutual induced */
/*     u3scale   field 1-3 group scale factor for mutual induced */
/*     u4scale   field 1-4 group scale factor for mutual induced */
/*     poltyp    type of polarization potential (direct or mutual) */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  vdwpot.i  --  specifics of van der Waals functional form  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     abuck       value of "A" constant in Buckingham vdw potential */
/*     bbuck       value of "B" constant in Buckingham vdw potential */
/*     cbuck       value of "C" constant in Buckingham vdw potential */
/*     ghal        value of "gamma" in buffered 14-7 vdw potential */
/*     dhal        value of "delta" in buffered 14-7 vdw potential */
/*     v2scale     factor by which 1-2 vdw interactions are scaled */
/*     v3scale     factor by which 1-3 vdw interactions are scaled */
/*     v4scale     factor by which 1-4 vdw interactions are scaled */
/*     v5scale     factor by which 1-5 vdw interactions are scaled */
/*     igauss      coefficients of Gaussian fit to vdw potential */
/*     ngauss      number of Gaussians used in fit to vdw potential */
/*     use_vcorr   flag to use long range vdw der Waals correction */
/*     vdwindex    indexing mode (atom type or class) for vdw parameters */
/*     vdwtyp      type of van der Waals potential energy function */
/*     radtyp      type of parameter (sigma or R-min) for atomic size */
/*     radsiz      atomic size provided as radius or diameter */
/*     radrule     combining rule for atomic size parameters */
/*     epsrule     combining rule for vdw well depth parameters */
/*     gausstyp    type of Gaussian fit to van der Waals potential */




/*     define blank character strings of various lengths */

    s_copy(blank3, "   ", (ftnlen)3, (ftnlen)3);
    s_copy(blank8, "        ", (ftnlen)8, (ftnlen)8);
    s_copy(blank12, "            ", (ftnlen)12, (ftnlen)12);
    s_copy(blank16, "                ", (ftnlen)16, (ftnlen)16);
    s_copy(blank20, "                    ", (ftnlen)20, (ftnlen)20);

/*     set the string value of the formfeed character (Ctrl-L) */

    *(unsigned char *)formfeed = '\f';

/*     force field atom type definitions */

    exist = FALSE_;
    for (i__ = 1; i__ <= 5000; ++i__) {
	if (s_cmp(symbol_ref(0, i__), blank3, (ftnlen)3, (ftnlen)3) != 0) {
	    exist = TRUE_;
	}
    }
    if (exist) {
	io___9.ciunit = *itxt;
	s_wsfe(&io___9);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___10.ciunit = *itxt;
	s_wsfe(&io___10);
	e_wsfe();
	for (i__ = 1; i__ <= 5000; ++i__) {
	    if (s_cmp(symbol_ref(0, i__), blank3, (ftnlen)3, (ftnlen)3) != 0) 
		    {
		io___11.ciunit = *itxt;
		s_wsfe(&io___11);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&katoms_1.atmcls[i__ - 1], (ftnlen)
			sizeof(integer));
		do_fio(&c__1, symbol_ref(0, i__), (ftnlen)3);
		do_fio(&c__1, describe_ref(0, i__), (ftnlen)24);
		do_fio(&c__1, (char *)&katoms_1.atmnum[i__ - 1], (ftnlen)
			sizeof(integer));
		do_fio(&c__1, (char *)&katoms_1.weight[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&katoms_1.ligand[i__ - 1], (ftnlen)
			sizeof(integer));
		e_wsfe();
	    }
	}
    }

/*     van der Waals parameters for atom types */

    exist = FALSE_;
    for (i__ = 1; i__ <= 5000; ++i__) {
	if (kvdws_1.rad[i__ - 1] != 0.) {
	    exist = TRUE_;
	}
    }
    if (exist) {
	io___12.ciunit = *itxt;
	s_wsfe(&io___12);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	if (s_cmp(vdwpot_1.vdwindex, "CLASS", (ftnlen)5, (ftnlen)5) == 0) {
	    io___13.ciunit = *itxt;
	    s_wsfe(&io___13);
	    e_wsfe();
	} else {
	    io___14.ciunit = *itxt;
	    s_wsfe(&io___14);
	    e_wsfe();
	}
	k = 0;
	for (i__ = 1; i__ <= 5000; ++i__) {
	    if (kvdws_1.rad[i__ - 1] != 0.) {
		++k;
		io___16.ciunit = *itxt;
		s_wsfe(&io___16);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kvdws_1.rad[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&kvdws_1.eps[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&kvdws_1.reduct[i__ - 1], (ftnlen)
			sizeof(doublereal));
		e_wsfe();
	    }
	}

/*     van der Waals scaling parameters */

	io___17.ciunit = *itxt;
	s_wsfe(&io___17);
	do_fio(&c__1, (char *)&vdwpot_1.v2scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&vdwpot_1.v3scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&vdwpot_1.v4scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&vdwpot_1.v5scale, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     van der Waals 1-4 parameters for atom types */

    exist = FALSE_;
    for (i__ = 1; i__ <= 5000; ++i__) {
	if (kvdws_1.rad4[i__ - 1] != 0.) {
	    exist = TRUE_;
	}
    }
    if (exist) {
	if (s_cmp(vdwpot_1.vdwindex, "CLASS", (ftnlen)5, (ftnlen)5) == 0) {
	    io___18.ciunit = *itxt;
	    s_wsfe(&io___18);
	    e_wsfe();
	} else {
	    io___19.ciunit = *itxt;
	    s_wsfe(&io___19);
	    e_wsfe();
	}
	k = 0;
	for (i__ = 1; i__ <= 5000; ++i__) {
	    if (kvdws_1.rad4[i__ - 1] != 0.) {
		++k;
		io___20.ciunit = *itxt;
		s_wsfe(&io___20);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kvdws_1.rad4[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&kvdws_1.eps4[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     van der Waals parameters for specific atom pairs */

    if (s_cmp(kvpr_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	if (s_cmp(vdwpot_1.vdwindex, "CLASS", (ftnlen)5, (ftnlen)5) == 0) {
	    io___21.ciunit = *itxt;
	    s_wsfe(&io___21);
	    e_wsfe();
	} else {
	    io___22.ciunit = *itxt;
	    s_wsfe(&io___22);
	    e_wsfe();
	}
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kvpr_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L150;
	    }
	    k1 = number_(kvpr_ref(0, i__), (ftnlen)4);
	    k2 = number_(kvpr_ref(4, i__), (ftnlen)4);
	    io___25.ciunit = *itxt;
	    s_wsfe(&io___25);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kvdwpr_1.radpr[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kvdwpr_1.epspr[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L150:
	;
    }

/*     hydrogen bonding parameters for specific atom pairs */

    if (s_cmp(khb_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	if (s_cmp(vdwpot_1.vdwindex, "CLASS", (ftnlen)5, (ftnlen)5) == 0) {
	    io___26.ciunit = *itxt;
	    s_wsfe(&io___26);
	    e_wsfe();
	} else {
	    io___27.ciunit = *itxt;
	    s_wsfe(&io___27);
	    e_wsfe();
	}
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(khb_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L190;
	    }
	    k1 = number_(khb_ref(0, i__), (ftnlen)4);
	    k2 = number_(khb_ref(4, i__), (ftnlen)4);
	    io___28.ciunit = *itxt;
	    s_wsfe(&io___28);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&khbond_1.radhb[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&khbond_1.epshb[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L190:
	;
    }

/*     bond stretching parameters */

    if (s_cmp(kb_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	io___29.ciunit = *itxt;
	s_wsfe(&io___29);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___30.ciunit = *itxt;
	s_wsfe(&io___30);
	e_wsfe();
	for (i__ = 1; i__ <= 2000; ++i__) {
	    if (s_cmp(kb_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L230;
	    }
	    k1 = number_(kb_ref(0, i__), (ftnlen)4);
	    k2 = number_(kb_ref(4, i__), (ftnlen)4);
	    io___31.ciunit = *itxt;
	    s_wsfe(&io___31);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kbonds_1.bcon[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kbonds_1.blen[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L230:
	;
    }

/*     bond stretching parameters for 5-membered rings */

    if (s_cmp(kb5_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	io___32.ciunit = *itxt;
	s_wsfe(&io___32);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kb5_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L260;
	    }
	    k1 = number_(kb5_ref(0, i__), (ftnlen)4);
	    k2 = number_(kb5_ref(4, i__), (ftnlen)4);
	    io___33.ciunit = *itxt;
	    s_wsfe(&io___33);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kbonds_1.bcon5[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kbonds_1.blen5[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L260:
	;
    }

/*     bond stretching parameters for 4-membered rings */

    if (s_cmp(kb4_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	io___34.ciunit = *itxt;
	s_wsfe(&io___34);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kb4_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L290;
	    }
	    k1 = number_(kb4_ref(0, i__), (ftnlen)4);
	    k2 = number_(kb4_ref(4, i__), (ftnlen)4);
	    io___35.ciunit = *itxt;
	    s_wsfe(&io___35);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kbonds_1.bcon4[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kbonds_1.blen4[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
/* L280: */
	}
L290:
	;
    }

/*     bond stretching parameters for 3-membered rings */

    if (s_cmp(kb3_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	io___36.ciunit = *itxt;
	s_wsfe(&io___36);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kb3_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L320;
	    }
	    k1 = number_(kb3_ref(0, i__), (ftnlen)4);
	    k2 = number_(kb3_ref(4, i__), (ftnlen)4);
	    io___37.ciunit = *itxt;
	    s_wsfe(&io___37);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kbonds_1.bcon3[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kbonds_1.blen3[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L320:
	;
    }

/*     cubic and quartic bond stretching parameters */

    if (bndpot_1.cbnd != 0. || bndpot_1.qbnd != 0.) {
	io___38.ciunit = *itxt;
	s_wsfe(&io___38);
	do_fio(&c__1, (char *)&bndpot_1.cbnd, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&bndpot_1.qbnd, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     electronegativity bond length correction parameters */

    if (s_cmp(kel_ref(0, 1), blank12, (ftnlen)12, (ftnlen)12) != 0) {
	io___39.ciunit = *itxt;
	s_wsfe(&io___39);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kel_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12) == 0) 
		    {
		goto L360;
	    }
	    k1 = number_(kel_ref(0, i__), (ftnlen)4);
	    k2 = number_(kel_ref(4, i__), (ftnlen)4);
	    k3 = number_(kel_ref(8, i__), (ftnlen)4);
	    io___41.ciunit = *itxt;
	    s_wsfe(&io___41);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kbonds_1.dlen[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L360:
	;
    }

/*     bond angle bending parameters */

    if (s_cmp(ka_ref(0, 1), blank12, (ftnlen)12, (ftnlen)12) != 0) {
	io___42.ciunit = *itxt;
	s_wsfe(&io___42);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___43.ciunit = *itxt;
	s_wsfe(&io___43);
	e_wsfe();
	for (i__ = 1; i__ <= 2000; ++i__) {
	    if (s_cmp(ka_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12) == 0) {
		goto L410;
	    }
	    k1 = number_(ka_ref(0, i__), (ftnlen)4);
	    k2 = number_(ka_ref(4, i__), (ftnlen)4);
	    k3 = number_(ka_ref(8, i__), (ftnlen)4);
	    if (ang_ref(2, i__) == 0. && ang_ref(3, i__) == 0.) {
		io___44.ciunit = *itxt;
		s_wsfe(&io___44);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kangs_1.acon[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&ang_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___45.ciunit = *itxt;
		s_wsfe(&io___45);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kangs_1.acon[i__ - 1], (ftnlen)sizeof(
			doublereal));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&ang_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	}
L410:
	;
    }

/*     bond angle bending parameters for 5-membered rings */

    if (s_cmp(ka5_ref(0, 1), blank12, (ftnlen)12, (ftnlen)12) != 0) {
	io___47.ciunit = *itxt;
	s_wsfe(&io___47);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(ka5_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12) == 0) 
		    {
		goto L450;
	    }
	    k1 = number_(ka5_ref(0, i__), (ftnlen)4);
	    k2 = number_(ka5_ref(4, i__), (ftnlen)4);
	    k3 = number_(ka5_ref(8, i__), (ftnlen)4);
	    if (ang5_ref(2, i__) == 0. && ang5_ref(3, i__) == 0.) {
		io___48.ciunit = *itxt;
		s_wsfe(&io___48);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kangs_1.acon5[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&ang5_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___49.ciunit = *itxt;
		s_wsfe(&io___49);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kangs_1.acon5[i__ - 1], (ftnlen)sizeof(
			doublereal));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&ang5_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	}
L450:
	;
    }

/*     bond angle bending parameters for 4-membered rings */

    if (s_cmp(ka4_ref(0, 1), blank12, (ftnlen)12, (ftnlen)12) != 0) {
	io___50.ciunit = *itxt;
	s_wsfe(&io___50);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(ka4_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12) == 0) 
		    {
		goto L490;
	    }
	    k1 = number_(ka4_ref(0, i__), (ftnlen)4);
	    k2 = number_(ka4_ref(4, i__), (ftnlen)4);
	    k3 = number_(ka4_ref(8, i__), (ftnlen)4);
	    if (ang4_ref(2, i__) == 0. && ang4_ref(3, i__) == 0.) {
		io___51.ciunit = *itxt;
		s_wsfe(&io___51);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kangs_1.acon4[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&ang4_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___52.ciunit = *itxt;
		s_wsfe(&io___52);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kangs_1.acon4[i__ - 1], (ftnlen)sizeof(
			doublereal));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&ang4_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	}
L490:
	;
    }

/*     bond angle bending parameters for 3-membered rings */

    if (s_cmp(ka3_ref(0, 1), blank12, (ftnlen)12, (ftnlen)12) != 0) {
	io___53.ciunit = *itxt;
	s_wsfe(&io___53);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(ka3_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12) == 0) 
		    {
		goto L530;
	    }
	    k1 = number_(ka3_ref(0, i__), (ftnlen)4);
	    k2 = number_(ka3_ref(4, i__), (ftnlen)4);
	    k3 = number_(ka3_ref(8, i__), (ftnlen)4);
	    if (ang3_ref(2, i__) == 0. && ang3_ref(3, i__) == 0.) {
		io___54.ciunit = *itxt;
		s_wsfe(&io___54);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kangs_1.acon3[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&ang3_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___55.ciunit = *itxt;
		s_wsfe(&io___55);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kangs_1.acon3[i__ - 1], (ftnlen)sizeof(
			doublereal));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&ang3_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	}
L530:
	;
    }

/*     Fourier bond angle bending parameters */

    if (s_cmp(kaf_ref(0, 1), blank12, (ftnlen)12, (ftnlen)12) != 0) {
	io___56.ciunit = *itxt;
	s_wsfe(&io___56);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kaf_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12) == 0) 
		    {
		goto L560;
	    }
	    k1 = number_(kaf_ref(0, i__), (ftnlen)4);
	    k2 = number_(kaf_ref(4, i__), (ftnlen)4);
	    k3 = number_(kaf_ref(8, i__), (ftnlen)4);
	    io___57.ciunit = *itxt;
	    s_wsfe(&io___57);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kangs_1.aconf[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    for (j = 1; j <= 2; ++j) {
		do_fio(&c__1, (char *)&angf_ref(j, i__), (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	}
L560:
	;
    }

/*     cubic through sextic bond angle bending parameters */

    if (angpot_1.cang != 0. || angpot_1.qang != 0. || angpot_1.pang != 0. || 
	    angpot_1.sang != 0.) {
	io___58.ciunit = *itxt;
	s_wsfe(&io___58);
	do_fio(&c__1, (char *)&angpot_1.cang, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&angpot_1.qang, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&angpot_1.pang, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&angpot_1.sang, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     stretch-bend parameters */

    if (s_cmp(ksb_ref(0, 1), blank12, (ftnlen)12, (ftnlen)12) != 0) {
	io___59.ciunit = *itxt;
	s_wsfe(&io___59);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___60.ciunit = *itxt;
	s_wsfe(&io___60);
	e_wsfe();
	for (i__ = 1; i__ <= 2000; ++i__) {
	    if (s_cmp(ksb_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12) == 0) 
		    {
		goto L610;
	    }
	    k1 = number_(ksb_ref(0, i__), (ftnlen)4);
	    k2 = number_(ksb_ref(4, i__), (ftnlen)4);
	    k3 = number_(ksb_ref(8, i__), (ftnlen)4);
	    io___61.ciunit = *itxt;
	    s_wsfe(&io___61);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&stbn_ref(1, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&stbn_ref(2, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L610:
	;
    }

/*     Urey-Bradley parameters */

    if (s_cmp(ku_ref(0, 1), blank12, (ftnlen)12, (ftnlen)12) != 0) {
	io___62.ciunit = *itxt;
	s_wsfe(&io___62);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___63.ciunit = *itxt;
	s_wsfe(&io___63);
	e_wsfe();
	for (i__ = 1; i__ <= 2000; ++i__) {
	    if (s_cmp(ku_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12) == 0) {
		goto L650;
	    }
	    k1 = number_(ku_ref(0, i__), (ftnlen)4);
	    k2 = number_(ku_ref(4, i__), (ftnlen)4);
	    k3 = number_(ku_ref(8, i__), (ftnlen)4);
	    io___64.ciunit = *itxt;
	    s_wsfe(&io___64);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kurybr_1.ucon[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kurybr_1.dst13[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L650:
	;
    }

/*     cubic and quartic Urey-Bradley parameters */

    if (urypot_1.cury != 0. || urypot_1.qury != 0.) {
	io___65.ciunit = *itxt;
	s_wsfe(&io___65);
	do_fio(&c__1, (char *)&urypot_1.cury, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&urypot_1.qury, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     angle-angle parameters */

    exist = FALSE_;
    for (i__ = 1; i__ <= 1000; ++i__) {
	for (k = 1; k <= 3; ++k) {
	    if (anan_ref(k, i__) != 0.) {
		exist = TRUE_;
	    }
	}
    }
    if (exist) {
	io___66.ciunit = *itxt;
	s_wsfe(&io___66);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___67.ciunit = *itxt;
	s_wsfe(&io___67);
	e_wsfe();
	k = 0;
	for (i__ = 1; i__ <= 1000; ++i__) {
	    if (anan_ref(1, i__) != 0. || anan_ref(2, i__) != 0. || anan_ref(
		    3, i__) != 0.) {
		++k;
		io___68.ciunit = *itxt;
		s_wsfe(&io___68);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&anan_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	}
    }

/*     out-of-plane bending parameters */

    if (s_cmp(kopb_ref(0, 1), blank16, (ftnlen)16, (ftnlen)16) != 0) {
	io___69.ciunit = *itxt;
	s_wsfe(&io___69);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___70.ciunit = *itxt;
	s_wsfe(&io___70);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kopb_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16) == 0)
		     {
		goto L730;
	    }
	    k1 = number_(kopb_ref(0, i__), (ftnlen)4);
	    k2 = number_(kopb_ref(4, i__), (ftnlen)4);
	    k3 = number_(kopb_ref(8, i__), (ftnlen)4);
	    k4 = number_(kopb_ref(12, i__), (ftnlen)4);
	    io___72.ciunit = *itxt;
	    s_wsfe(&io___72);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k4, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kopbnd_1.opbn[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L730:
	;
    }

/*     out-of-plane distance parameters */

    if (s_cmp(kopd_ref(0, 1), blank16, (ftnlen)16, (ftnlen)16) != 0) {
	io___73.ciunit = *itxt;
	s_wsfe(&io___73);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___74.ciunit = *itxt;
	s_wsfe(&io___74);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kopd_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16) == 0)
		     {
		goto L770;
	    }
	    k1 = number_(kopd_ref(0, i__), (ftnlen)4);
	    k2 = number_(kopd_ref(4, i__), (ftnlen)4);
	    k3 = number_(kopd_ref(8, i__), (ftnlen)4);
	    k4 = number_(kopd_ref(12, i__), (ftnlen)4);
	    io___75.ciunit = *itxt;
	    s_wsfe(&io___75);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k4, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kopdst_1.opds[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L770:
	;
    }

/*     improper dihedral parameters */

    if (s_cmp(kdi_ref(0, 1), blank16, (ftnlen)16, (ftnlen)16) != 0) {
	io___76.ciunit = *itxt;
	s_wsfe(&io___76);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___77.ciunit = *itxt;
	s_wsfe(&io___77);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kdi_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16) == 0) 
		    {
		goto L810;
	    }
	    k1 = number_(kdi_ref(0, i__), (ftnlen)4);
	    k2 = number_(kdi_ref(4, i__), (ftnlen)4);
	    k3 = number_(kdi_ref(8, i__), (ftnlen)4);
	    k4 = number_(kdi_ref(12, i__), (ftnlen)4);
	    io___78.ciunit = *itxt;
	    s_wsfe(&io___78);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k4, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kiprop_1.dcon[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kiprop_1.tdi[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L810:
	;
    }

/*     improper torsional parameters */

    if (s_cmp(kti_ref(0, 1), blank16, (ftnlen)16, (ftnlen)16) != 0) {
	io___79.ciunit = *itxt;
	s_wsfe(&io___79);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___80.ciunit = *itxt;
	s_wsfe(&io___80);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kti_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16) == 0) 
		    {
		goto L850;
	    }
	    k1 = number_(kti_ref(0, i__), (ftnlen)4);
	    k2 = number_(kti_ref(4, i__), (ftnlen)4);
	    k3 = number_(kti_ref(8, i__), (ftnlen)4);
	    k4 = number_(kti_ref(12, i__), (ftnlen)4);
	    j = 0;
	    if (ti1_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 1;
		ampli[j - 1] = ti1_ref(1, i__);
		phase[j - 1] = ti1_ref(2, i__);
	    }
	    if (ti2_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 2;
		ampli[j - 1] = ti2_ref(1, i__);
		phase[j - 1] = ti2_ref(2, i__);
	    }
	    if (ti3_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 3;
		ampli[j - 1] = ti3_ref(1, i__);
		phase[j - 1] = ti3_ref(2, i__);
	    }
	    io___84.ciunit = *itxt;
	    s_wsfe(&io___84);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k4, (ftnlen)sizeof(integer));
	    i__1 = j;
	    for (k = 1; k <= i__1; ++k) {
		do_fio(&c__1, (char *)&ampli[k - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&phase[k - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&fold[k - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	}
L850:
	;
    }

/*     torsional angle parameters */

    if (s_cmp(kt_ref(0, 1), blank16, (ftnlen)16, (ftnlen)16) != 0) {
	io___85.ciunit = *itxt;
	s_wsfe(&io___85);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___86.ciunit = *itxt;
	s_wsfe(&io___86);
	e_wsfe();
	for (i__ = 1; i__ <= 2000; ++i__) {
	    if (s_cmp(kt_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16) == 0) {
		goto L890;
	    }
	    k1 = number_(kt_ref(0, i__), (ftnlen)4);
	    k2 = number_(kt_ref(4, i__), (ftnlen)4);
	    k3 = number_(kt_ref(8, i__), (ftnlen)4);
	    k4 = number_(kt_ref(12, i__), (ftnlen)4);
	    j = 0;
	    if (t1_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 1;
		ampli[j - 1] = t1_ref(1, i__);
		phase[j - 1] = t1_ref(2, i__);
	    }
	    if (t2_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 2;
		ampli[j - 1] = t2_ref(1, i__);
		phase[j - 1] = t2_ref(2, i__);
	    }
	    if (t3_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 3;
		ampli[j - 1] = t3_ref(1, i__);
		phase[j - 1] = t3_ref(2, i__);
	    }
	    if (t4_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 4;
		ampli[j - 1] = t4_ref(1, i__);
		phase[j - 1] = t4_ref(2, i__);
	    }
	    if (t5_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 5;
		ampli[j - 1] = t5_ref(1, i__);
		phase[j - 1] = t5_ref(2, i__);
	    }
	    if (t6_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 6;
		ampli[j - 1] = t6_ref(1, i__);
		phase[j - 1] = t6_ref(2, i__);
	    }
	    io___87.ciunit = *itxt;
	    s_wsfe(&io___87);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k4, (ftnlen)sizeof(integer));
	    i__1 = j;
	    for (k = 1; k <= i__1; ++k) {
		do_fio(&c__1, (char *)&ampli[k - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&phase[k - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&fold[k - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	}
L890:
	;
    }

/*     torsional angle parameters for 5-membered rings */

    if (s_cmp(kt5_ref(0, 1), blank16, (ftnlen)16, (ftnlen)16) != 0) {
	io___88.ciunit = *itxt;
	s_wsfe(&io___88);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kt5_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16) == 0) 
		    {
		goto L920;
	    }
	    k1 = number_(kt5_ref(0, i__), (ftnlen)4);
	    k2 = number_(kt5_ref(4, i__), (ftnlen)4);
	    k3 = number_(kt5_ref(8, i__), (ftnlen)4);
	    k4 = number_(kt5_ref(12, i__), (ftnlen)4);
	    j = 0;
	    if (t15_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 1;
		ampli[j - 1] = t15_ref(1, i__);
		phase[j - 1] = t15_ref(2, i__);
	    }
	    if (t25_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 2;
		ampli[j - 1] = t25_ref(1, i__);
		phase[j - 1] = t25_ref(2, i__);
	    }
	    if (t35_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 3;
		ampli[j - 1] = t35_ref(1, i__);
		phase[j - 1] = t35_ref(2, i__);
	    }
	    if (t45_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 4;
		ampli[j - 1] = t45_ref(1, i__);
		phase[j - 1] = t45_ref(2, i__);
	    }
	    if (t55_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 5;
		ampli[j - 1] = t55_ref(1, i__);
		phase[j - 1] = t55_ref(2, i__);
	    }
	    if (t65_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 6;
		ampli[j - 1] = t65_ref(1, i__);
		phase[j - 1] = t65_ref(2, i__);
	    }
	    io___89.ciunit = *itxt;
	    s_wsfe(&io___89);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k4, (ftnlen)sizeof(integer));
	    i__1 = j;
	    for (k = 1; k <= i__1; ++k) {
		do_fio(&c__1, (char *)&ampli[k - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&phase[k - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&fold[k - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	}
L920:
	;
    }

/*     torsional angle parameters for 4-membered rings */

    if (s_cmp(kt4_ref(0, 1), blank16, (ftnlen)16, (ftnlen)16) != 0) {
	io___90.ciunit = *itxt;
	s_wsfe(&io___90);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kt4_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16) == 0) 
		    {
		goto L950;
	    }
	    k1 = number_(kt4_ref(0, i__), (ftnlen)4);
	    k2 = number_(kt4_ref(4, i__), (ftnlen)4);
	    k3 = number_(kt4_ref(8, i__), (ftnlen)4);
	    k4 = number_(kt4_ref(12, i__), (ftnlen)4);
	    j = 0;
	    if (t14_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 1;
		ampli[j - 1] = t14_ref(1, i__);
		phase[j - 1] = t14_ref(2, i__);
	    }
	    if (t24_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 2;
		ampli[j - 1] = t24_ref(1, i__);
		phase[j - 1] = t24_ref(2, i__);
	    }
	    if (t34_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 3;
		ampli[j - 1] = t34_ref(1, i__);
		phase[j - 1] = t34_ref(2, i__);
	    }
	    if (t44_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 4;
		ampli[j - 1] = t44_ref(1, i__);
		phase[j - 1] = t44_ref(2, i__);
	    }
	    if (t54_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 5;
		ampli[j - 1] = t54_ref(1, i__);
		phase[j - 1] = t54_ref(2, i__);
	    }
	    if (t64_ref(1, i__) != 0.) {
		++j;
		fold[j - 1] = 6;
		ampli[j - 1] = t64_ref(1, i__);
		phase[j - 1] = t64_ref(2, i__);
	    }
	    io___91.ciunit = *itxt;
	    s_wsfe(&io___91);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k4, (ftnlen)sizeof(integer));
	    i__1 = j;
	    for (k = 1; k <= i__1; ++k) {
		do_fio(&c__1, (char *)&ampli[k - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&phase[k - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&fold[k - 1], (ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	}
L950:
	;
    }

/*     pi-orbital torsion parameters */

    if (s_cmp(kpt_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	io___92.ciunit = *itxt;
	s_wsfe(&io___92);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___93.ciunit = *itxt;
	s_wsfe(&io___93);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kpt_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L990;
	    }
	    k1 = number_(kpt_ref(0, i__), (ftnlen)4);
	    k2 = number_(kpt_ref(4, i__), (ftnlen)4);
	    io___94.ciunit = *itxt;
	    s_wsfe(&io___94);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kpitor_1.ptcon[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L990:
	;
    }

/*     stretch-torsion parameters */

    if (s_cmp(kbt_ref(0, 1), blank16, (ftnlen)16, (ftnlen)16) != 0) {
	io___95.ciunit = *itxt;
	s_wsfe(&io___95);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___96.ciunit = *itxt;
	s_wsfe(&io___96);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kbt_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16) == 0) 
		    {
		goto L1030;
	    }
	    k1 = number_(kbt_ref(0, i__), (ftnlen)4);
	    k2 = number_(kbt_ref(4, i__), (ftnlen)4);
	    k3 = number_(kbt_ref(8, i__), (ftnlen)4);
	    k4 = number_(kbt_ref(12, i__), (ftnlen)4);
	    io___97.ciunit = *itxt;
	    s_wsfe(&io___97);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k4, (ftnlen)sizeof(integer));
	    for (j = 1; j <= 3; ++j) {
		do_fio(&c__1, (char *)&btcon_ref(j, i__), (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	}
L1030:
	;
    }

/*     torsion-torsion parameters */

    if (s_cmp(ktt_ref(0, 1), blank20, (ftnlen)20, (ftnlen)20) != 0) {
	io___98.ciunit = *itxt;
	s_wsfe(&io___98);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___99.ciunit = *itxt;
	s_wsfe(&io___99);
	e_wsfe();
	for (i__ = 1; i__ <= 100; ++i__) {
	    if (s_cmp(ktt_ref(0, i__), blank20, (ftnlen)20, (ftnlen)20) == 0) 
		    {
		goto L1080;
	    }
	    k1 = number_(ktt_ref(0, i__), (ftnlen)4);
	    k2 = number_(ktt_ref(4, i__), (ftnlen)4);
	    k3 = number_(ktt_ref(8, i__), (ftnlen)4);
	    k4 = number_(ktt_ref(12, i__), (ftnlen)4);
	    k5 = number_(ktt_ref(16, i__), (ftnlen)4);
	    io___101.ciunit = *itxt;
	    s_wsfe(&io___101);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k4, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k5, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ktrtor_1.tnx[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&ktrtor_1.tny[i__ - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	    k = ktrtor_1.tnx[i__ - 1] * ktrtor_1.tny[i__ - 1];
	    io___102.ciunit = *itxt;
	    s_wsfe(&io___102);
	    i__1 = k;
	    for (j = 1; j <= i__1; ++j) {
		do_fio(&c__1, (char *)&tbf_ref(j, i__), (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	}
L1080:
	;
    }

/*     atomic partial charge parameters */

    exist = FALSE_;
    for (i__ = 1; i__ <= 5000; ++i__) {
	if (kchrge_1.chg[i__ - 1] != 0.) {
	    exist = TRUE_;
	}
    }
    if (exist) {
	io___103.ciunit = *itxt;
	s_wsfe(&io___103);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___104.ciunit = *itxt;
	s_wsfe(&io___104);
	e_wsfe();
	k = 0;
	for (i__ = 1; i__ <= 5000; ++i__) {
	    if (kchrge_1.chg[i__ - 1] != 0.) {
		++k;
		io___105.ciunit = *itxt;
		s_wsfe(&io___105);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kchrge_1.chg[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}

/*     atomic partial charge scaling parameters */

	io___106.ciunit = *itxt;
	s_wsfe(&io___106);
	do_fio(&c__1, (char *)&chgpot_1.c2scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&chgpot_1.c3scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&chgpot_1.c4scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&chgpot_1.c5scale, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     bond dipole moment parameters */

    if (s_cmp(kd_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	io___107.ciunit = *itxt;
	s_wsfe(&io___107);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___108.ciunit = *itxt;
	s_wsfe(&io___108);
	e_wsfe();
	for (i__ = 1; i__ <= 1000; ++i__) {
	    if (s_cmp(kd_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L1160;
	    }
	    k1 = number_(kd_ref(0, i__), (ftnlen)4);
	    k2 = number_(kd_ref(4, i__), (ftnlen)4);
	    io___109.ciunit = *itxt;
	    s_wsfe(&io___109);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kdipol_1.dpl[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kdipol_1.pos[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L1160:
	;
    }

/*     bond dipole moment parameters for 5-membered rings */

    if (s_cmp(kd5_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	io___110.ciunit = *itxt;
	s_wsfe(&io___110);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kd5_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L1190;
	    }
	    k1 = number_(kd5_ref(0, i__), (ftnlen)4);
	    k2 = number_(kd5_ref(4, i__), (ftnlen)4);
	    io___111.ciunit = *itxt;
	    s_wsfe(&io___111);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kdipol_1.dpl5[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kdipol_1.pos5[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L1190:
	;
    }

/*     bond dipole moment parameters for 4-membered rings */

    if (s_cmp(kd4_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	io___112.ciunit = *itxt;
	s_wsfe(&io___112);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kd4_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L1220;
	    }
	    k1 = number_(kd4_ref(0, i__), (ftnlen)4);
	    k2 = number_(kd4_ref(4, i__), (ftnlen)4);
	    io___113.ciunit = *itxt;
	    s_wsfe(&io___113);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kdipol_1.dpl4[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kdipol_1.pos4[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L1220:
	;
    }

/*     bond dipole moment parameters for 3-membered rings */

    if (s_cmp(kd3_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	io___114.ciunit = *itxt;
	s_wsfe(&io___114);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kd3_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L1250;
	    }
	    k1 = number_(kd3_ref(0, i__), (ftnlen)4);
	    k2 = number_(kd3_ref(4, i__), (ftnlen)4);
	    io___115.ciunit = *itxt;
	    s_wsfe(&io___115);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&kdipol_1.dpl3[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&kdipol_1.pos3[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L1250:
	;
    }

/*     atomic multipole electrostatic parameters */

    if (s_cmp(kmp_ref(0, 1), blank16, (ftnlen)16, (ftnlen)16) != 0) {
	io___116.ciunit = *itxt;
	s_wsfe(&io___116);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___117.ciunit = *itxt;
	s_wsfe(&io___117);
	e_wsfe();
	for (i__ = 1; i__ <= 2000; ++i__) {
	    if (s_cmp(kmp_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16) == 0) 
		    {
		goto L1290;
	    }
	    k1 = number_(kmp_ref(0, i__), (ftnlen)4);
	    k2 = number_(kmp_ref(4, i__), (ftnlen)4);
	    k3 = number_(kmp_ref(8, i__), (ftnlen)4);
	    k4 = number_(kmp_ref(12, i__), (ftnlen)4);
	    io___118.ciunit = *itxt;
	    s_wsfe(&io___118);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k4, (ftnlen)sizeof(integer));
	    do_fio(&c__1, mpaxis_ref(0, i__), (ftnlen)8);
	    do_fio(&c__1, (char *)&multip_ref(1, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&multip_ref(2, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&multip_ref(3, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&multip_ref(4, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&multip_ref(5, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&multip_ref(8, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&multip_ref(9, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&multip_ref(11, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&multip_ref(12, i__), (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&multip_ref(13, i__), (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L1290:

/*     atomic multipole scaling parameters */

	io___119.ciunit = *itxt;
	s_wsfe(&io___119);
	do_fio(&c__1, (char *)&mplpot_1.m2scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&mplpot_1.m3scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&mplpot_1.m4scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&mplpot_1.m5scale, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     atomic dipole polarizability parameters */

    exist = FALSE_;
    for (i__ = 1; i__ <= 5000; ++i__) {
	if (kpolr_1.polr[i__ - 1] != 0.) {
	    exist = TRUE_;
	}
    }
    if (exist) {
	io___120.ciunit = *itxt;
	s_wsfe(&io___120);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___121.ciunit = *itxt;
	s_wsfe(&io___121);
	e_wsfe();
	k = 0;
	for (i__ = 1; i__ <= 5000; ++i__) {
	    if (kpolr_1.polr[i__ - 1] != 0.) {
		++k;
		npg = 0;
		for (j = 1; j <= 8; ++j) {
		    if (pgrp_ref(j, i__) != 0) {
			++npg;
		    }
		}
		if (npg == 0) {
		    io___123.ciunit = *itxt;
		    s_wsfe(&io___123);
		    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&kpolr_1.polr[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&kpolr_1.athl[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    e_wsfe();
		} else {
		    io___124.ciunit = *itxt;
		    s_wsfe(&io___124);
		    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&kpolr_1.polr[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&kpolr_1.athl[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    i__1 = npg;
		    for (j = 1; j <= i__1; ++j) {
			do_fio(&c__1, (char *)&pgrp_ref(j, i__), (ftnlen)
				sizeof(integer));
		    }
		    e_wsfe();
		}
	    }
	}

/*     dipole polarizability scaling parameters */

	io___125.ciunit = *itxt;
	s_wsfe(&io___125);
	do_fio(&c__1, (char *)&polpot_1.d1scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&polpot_1.d2scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&polpot_1.d3scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&polpot_1.d4scale, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___126.ciunit = *itxt;
	s_wsfe(&io___126);
	do_fio(&c__1, (char *)&polpot_1.u1scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&polpot_1.u2scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&polpot_1.u3scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&polpot_1.u4scale, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___127.ciunit = *itxt;
	s_wsfe(&io___127);
	do_fio(&c__1, (char *)&polpot_1.p2scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&polpot_1.p3scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&polpot_1.p4scale, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&polpot_1.p5scale, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     conjugated pisystem atom parameters */

    exist = FALSE_;
    for (i__ = 1; i__ <= 1000; ++i__) {
	if (korbs_1.ionize[i__ - 1] != 0.) {
	    exist = TRUE_;
	}
    }
    if (exist) {
	io___128.ciunit = *itxt;
	s_wsfe(&io___128);
	do_fio(&c__1, formfeed, (ftnlen)1);
	do_fio(&c__1, fields_1.forcefield, (ftnlen)20);
	e_wsfe();
	io___129.ciunit = *itxt;
	s_wsfe(&io___129);
	e_wsfe();
	k = 0;
	for (i__ = 1; i__ <= 1000; ++i__) {
	    if (korbs_1.ionize[i__ - 1] != 0.) {
		++k;
		io___130.ciunit = *itxt;
		s_wsfe(&io___130);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&korbs_1.electron[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&korbs_1.ionize[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&korbs_1.repulse[i__ - 1], (ftnlen)
			sizeof(doublereal));
		e_wsfe();
	    }
	}
    }

/*     conjugated pisystem bond parameters */

    if (s_cmp(kpi_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	io___131.ciunit = *itxt;
	s_wsfe(&io___131);
	e_wsfe();
	for (i__ = 1; i__ <= 500; ++i__) {
	    if (s_cmp(kpi_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L1430;
	    }
	    k1 = number_(kpi_ref(0, i__), (ftnlen)4);
	    k2 = number_(kpi_ref(4, i__), (ftnlen)4);
	    io___132.ciunit = *itxt;
	    s_wsfe(&io___132);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&korbs_1.sslope[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&korbs_1.tslope[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L1430:
	;
    }

/*     conjugated pisystem bond parameters for 5-membered rings */

    if (s_cmp(kpi5_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	io___133.ciunit = *itxt;
	s_wsfe(&io___133);
	e_wsfe();
	for (i__ = 1; i__ <= 200; ++i__) {
	    if (s_cmp(kpi5_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L1460;
	    }
	    k1 = number_(kpi5_ref(0, i__), (ftnlen)4);
	    k2 = number_(kpi5_ref(4, i__), (ftnlen)4);
	    io___134.ciunit = *itxt;
	    s_wsfe(&io___134);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&korbs_1.sslope5[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&korbs_1.tslope5[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L1460:
	;
    }

/*     conjugated pisystem bond parameters for 4-membered rings */

    if (s_cmp(kpi4_ref(0, 1), blank8, (ftnlen)8, (ftnlen)8) != 0) {
	io___135.ciunit = *itxt;
	s_wsfe(&io___135);
	e_wsfe();
	for (i__ = 1; i__ <= 200; ++i__) {
	    if (s_cmp(kpi4_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8) == 0) {
		goto L1490;
	    }
	    k1 = number_(kpi4_ref(0, i__), (ftnlen)4);
	    k2 = number_(kpi4_ref(4, i__), (ftnlen)4);
	    io___136.ciunit = *itxt;
	    s_wsfe(&io___136);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&k2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&korbs_1.sslope4[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&korbs_1.tslope4[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
L1490:
	;
    }
    return 0;
} /* prtprm_ */

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


