/* testgrad.f -- translated by f2c (version 20050501).
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
    doublereal x[25000], y[25000], z__[25000];
    integer n, type__[25000];
} atoms_;

#define atoms_1 atoms_

struct {
    doublereal desum[75000]	/* was [3][25000] */, deb[75000]	/* 
	    was [3][25000] */, dea[75000]	/* was [3][25000] */, deba[
	    75000]	/* was [3][25000] */, deub[75000]	/* was [3][
	    25000] */, deaa[75000]	/* was [3][25000] */, deopb[75000]	
	    /* was [3][25000] */, deopd[75000]	/* was [3][25000] */, deid[
	    75000]	/* was [3][25000] */, deit[75000]	/* was [3][
	    25000] */, det[75000]	/* was [3][25000] */, dept[75000]	
	    /* was [3][25000] */, debt[75000]	/* was [3][25000] */, dett[
	    75000]	/* was [3][25000] */, dev[75000]	/* was [3][
	    25000] */, dec[75000]	/* was [3][25000] */, decd[75000]	
	    /* was [3][25000] */, ded[75000]	/* was [3][25000] */, dem[
	    75000]	/* was [3][25000] */, dep[75000]	/* was [3][
	    25000] */, der[75000]	/* was [3][25000] */, des[75000]	
	    /* was [3][25000] */, delf[75000]	/* was [3][25000] */, deg[
	    75000]	/* was [3][25000] */, dex[75000]	/* was [3][
	    25000] */;
} deriv_;

#define deriv_1 deriv_

struct {
    doublereal esum, eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, 
	    ebt, ett, ev, ec, ecd, ed, em, ep, er, es, elf, eg, ex;
} energi_;

#define energi_1 energi_

struct {
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    doublereal einter;
} inter_;

#define inter_1 inter_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal rsolv[25000], asolv[25000], rborn[25000], drb[25000], drbp[
	    25000], drobc[25000], doffset, p1, p2, p3, p4, p5, gpol[25000], 
	    shct[25000], aobc[25000], bobc[25000], gobc[25000], vsolv[25000], 
	    wace[1000000]	/* was [1000][1000] */, s2ace[1000000]	/* 
	    was [1000][1000] */, uace[1000000]	/* was [1000][1000] */;
    char solvtyp[8], borntyp[8];
} solute_;

#define solute_1 solute_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__1 = 1;
static integer c__5 = 5;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  program testgrad  --  derivative test; Cartesian version  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "testgrad" computes and compares the analytical and numerical */
/*     gradient vectors of the potential energy function with respect */
/*     to Cartesian coordinates */


/* Main program */ int MAIN__(void)
{
    /* Initialized data */

    static char axis[1*3] = "X" "Y" "Z";

    /* Format strings */
    static char fmt_10[] = "(/,\002 Compute the Analytical Gradient Vector ["
	    "Y] :  \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(/,\002 Compute the Numerical Gradient Vector [Y"
	    "] :   \002,$)";
    static char fmt_40[] = "(a120)";
    static char fmt_60[] = "(/,\002 Enter a Numerical Stepsize [\002,d7.1"
	    ",\002 Ang] :  \002,$)";
    static char fmt_70[] = "(f20.0)";
    static char fmt_80[] = "(/,\002 Output Breakdown by Gradient Componen"
	    "t\002,\002 [N] :  \002,$)";
    static char fmt_90[] = "(a120)";
    static char fmt_100[] = "(/,\002 Total Potential Energy :\002,8x,f20.8"
	    ",\002 Kcal/mole\002)";
    static char fmt_110[] = "(/,\002 Total Potential Energy :\002,8x,f18.6"
	    ",\002 Kcal/mole\002)";
    static char fmt_120[] = "(/,\002 Total Potential Energy :\002,8x,f16.4"
	    ",\002 Kcal/mole\002)";
    static char fmt_130[] = "(/,\002 Potential Energy Breakdown by Individ"
	    "ual\002,\002 Components :\002)";
    static char fmt_140[] = "(/,\002  Energy\002,7x,\002EB\002,14x,\002EA"
	    "\002,14x,\002EBA\002,13x,\002EUB\002,/,\002  Terms\002,8x,\002EAA"
	    "\002,13x,\002EOPB\002,12x,\002EOPD\002,12x,\002EID\002,/,15x,"
	    "\002EIT\002,13x,\002ET\002,14x,\002EPT\002,13x,\002EBT\002,/,15x,"
	    "\002ETT\002,13x,\002EV\002,14x,\002EC\002,14x,\002ECD\002,/,15x"
	    ",\002ED\002,14x,\002EM\002,14x,\002EP\002,14x,\002ER\002,/,15x"
	    ",\002ES\002,14x,\002ELF\002,13x,\002EG\002,14x,\002EX\002)";
    static char fmt_150[] = "(/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16."
	    "8,/,6x,4f16.8,/,6x,4f16.8)";
    static char fmt_160[] = "(/,\002  Energy\002,6x,\002EB\002,12x,\002EA"
	    "\002,12x,\002EBA\002,11x,\002EUB\002,11x,\002EAA\002,/,\002  Ter"
	    "ms\002,7x,\002EOPB\002,10x,\002EOPD\002,10x,\002EID\002,11x,\002"
	    "EIT\002,11x,\002ET\002,/,14x,\002EPT\002,11x,\002EBT\002,11x,"
	    "\002ETT\002,11x,\002EV\002,12x,\002EC\002,/,14x,\002ECD\002,11x"
	    ",\002ED\002,12x,\002EM\002,12x,\002EP\002,12x,\002ER\002,/,14x"
	    ",\002ES\002,12x,\002ELF\002,11x,\002EG\002,12x,\002EX\002)";
    static char fmt_170[] = "(/,6x,5f14.6,/,6x,5f14.6,/,6x,5f14.6,/,6x,5f14."
	    "6,/,6x,4f14.6)";
    static char fmt_180[] = "(/,\002  Energy\002,6x,\002EB\002,10x,\002EA"
	    "\002,10x,\002EBA\002,9x,\002EUB\002,9x,\002EAA\002,9x,\002EOP"
	    "B\002,/,\002  Terms\002,7x,\002EOPD\002,8x,\002EID\002,9x,\002EIT"
	    "\002,9x,\002ET\002,10x,\002EPT\002,9x,\002EBT\002,/,14x,\002ET"
	    "T\002,9x,\002EV\002,10x,\002EC\002,10x,\002ECD\002,9x,\002ED\002"
	    ",10x,\002EM\002,/,14x,\002EP\002,10x,\002ER\002,10x,\002ES\002,1"
	    "0x,\002ELF\002,9x,\002EG\002,10x,\002EX\002)";
    static char fmt_190[] = "(/,6x,6f12.4,/,6x,6f12.4,/,6x,6f12.4,/,6x,6f12."
	    "4)";
    static char fmt_200[] = "(/,\002 Cartesian Gradient Breakdown by Individ"
	    "ual\002,\002 Components :\002)";
    static char fmt_210[] = "(/,2x,\002Atom\002,9x,\002d EB\002,12x,\002d E"
	    "A\002,12x,\002d EBA\002,11x,\002d EUB\002,/,2x,\002Axis\002,9x"
	    ",\002d EAA\002,11x,\002d EOPB\002,10x,\002d EOPD\002,10x,\002d E"
	    "ID\002,/,2x,\002Type\002,9x,\002d EIT\002,11x,\002d ET\002,12x"
	    ",\002d EPT\002,11x,\002d EBT\002,/,15x,\002d ETT\002,11x,\002d EV"
	    "\002,12x,\002d EC\002,12x,\002d ECD\002,/,15x,\002d ED\002,12x"
	    ",\002d EM\002,12x,\002d EP\002,12x,\002d ER\002,/,15x,\002d E"
	    "S\002,12x,\002d ELF\002,11x,\002d EG\002,12x,\002d EX\002)";
    static char fmt_220[] = "(/,2x,\002Atom\002,8x,\002d EB\002,10x,\002d E"
	    "A\002,10x,\002d EBA\002,9x,\002d EUB\002,9x,\002d EAA\002,/,2x"
	    ",\002Axis\002,8x,\002d EOPB\002,8x,\002d EOPD\002,8x,\002d EI"
	    "D\002,9x,\002d EIT\002,9x,\002d ET\002,/,2x,\002Type\002,8x,\002"
	    "d EPT\002,9x,\002d EBT\002,9x,\002d ETT\002,9x,\002d EV\002,10x"
	    ",\002d EC\002,/,14x,\002d ECD\002,9x,\002d ED\002,10x,\002d E"
	    "M\002,10x,\002d EP\002,10x,\002d ER\002,/,14x,\002d ES\002,10x"
	    ",\002d ELF\002,9x,\002d EG\002,10x,\002d EX\002)";
    static char fmt_230[] = "(/,2x,\002Atom\002,6x,\002d EB\002,8x,\002d E"
	    "A\002,8x,\002d EBA\002,7x,\002d EUB\002,7x,\002d EAA\002,7x,\002"
	    "d EOPB\002,/,2x,\002Axis\002,6x,\002d EOPD\002,6x,\002d EID\002,"
	    "7x,\002d EIT\002,7x,\002d ET\002,8x,\002d EPT\002,7x,\002d EB"
	    "T\002,/,2x,\002Type\002,6x,\002d ETT\002,7x,\002d EV\002,8x,\002"
	    "d EC\002,8x,\002d ECD\002,7x,\002d ED\002,8x,\002d EM\002,/,12x"
	    ",\002d EP\002,8x,\002d ER\002,8x,\002d ES\002,8x,\002d ELF\002,7"
	    "x,\002d EG\002,8x,\002d EX\002)";
    static char fmt_240[] = "(/,i6,4f16.8,/,5x,a1,4f16.8,/,\002 Anlyt\002,4f"
	    "16.8,/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8)";
    static char fmt_250[] = "(/,i6,5f14.6,/,5x,a1,5f14.6,/,\002 Anlyt\002,5f"
	    "14.6,/,6x,5f14.6,/,6x,4f14.6)";
    static char fmt_260[] = "(/,i6,6f12.4,/,5x,a1,6f12.4,/,\002 Anlyt\002,6f"
	    "12.4,/,6x,6f12.4)";
    static char fmt_270[] = "(/,i6,4f16.8,/,5x,a1,4f16.8,/,\002 Numer\002,4f"
	    "16.8,/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8)";
    static char fmt_280[] = "(/,i6,5f14.6,/,5x,a1,5f14.6,/,\002 Numer\002,5f"
	    "14.6,/,6x,5f14.6,/,6x,4f14.6)";
    static char fmt_290[] = "(/,i6,6f12.4,/,5x,a1,6f12.4,/,\002 Numer\002,6f"
	    "12.4,/,6x,6f12.4)";
    static char fmt_300[] = "(/,\002 Cartesian Gradient Breakdown over Indiv"
	    "idual Atoms :\002)";
    static char fmt_310[] = "(/,2x,\002Type\002,4x,\002Atom\002,10x,\002dE"
	    "/dX\002,11x,\002dE/dY\002,11x,\002dE/dZ\002,11x,\002Norm\002,/)";
    static char fmt_320[] = "(/,2x,\002Type\002,6x,\002Atom\002,11x,\002dE"
	    "/dX\002,9x,\002dE/dY\002,9x,\002dE/dZ\002,11x,\002Norm\002,/)";
    static char fmt_330[] = "(/,2x,\002Type\002,6x,\002Atom\002,14x,\002dE"
	    "/dX\002,7x,\002dE/dY\002,7x,\002dE/dZ\002,10x,\002Norm\002,/)";
    static char fmt_340[] = "(\002 Anlyt\002,i8,1x,3f16.8,f16.8)";
    static char fmt_350[] = "(\002 Anlyt\002,2x,i8,3x,3f14.6,2x,f14.6)";
    static char fmt_360[] = "(\002 Anlyt\002,2x,i8,7x,3f12.4,2x,f12.4)";
    static char fmt_370[] = "(\002 Numer\002,i8,1x,3f16.8,f16.8)";
    static char fmt_380[] = "(\002 Numer\002,2x,i8,3x,3f14.6,2x,f14.6)";
    static char fmt_390[] = "(\002 Numer\002,2x,i8,7x,3f12.4,2x,f12.4)";
    static char fmt_400[] = "(/,\002 Total Gradient Norm and RMS Gradient pe"
	    "r Atom :\002)";
    static char fmt_410[] = "(/,\002 Anlyt\002,6x,\002Total Gradient Norm Va"
	    "lue\002,6x,f20.8)";
    static char fmt_420[] = "(/,\002 Anlyt\002,6x,\002Total Gradient Norm Va"
	    "lue\002,6x,f18.6)";
    static char fmt_430[] = "(/,\002 Anlyt\002,6x,\002Total Gradient Norm Va"
	    "lue\002,6x,f16.4)";
    static char fmt_440[] = "(\002 Numer\002,6x,\002Total Gradient Norm Va"
	    "lue\002,6x,f20.8)";
    static char fmt_450[] = "(\002 Numer\002,6x,\002Total Gradient Norm Va"
	    "lue\002,6x,f18.6)";
    static char fmt_460[] = "(\002 Numer\002,6x,\002Total Gradient Norm Va"
	    "lue\002,6x,f16.4)";
    static char fmt_470[] = "(/,\002 Anlyt\002,6x,\002RMS Gradient over All "
	    "Atoms\002,4x,f20.8)";
    static char fmt_480[] = "(/,\002 Anlyt\002,6x,\002RMS Gradient over All "
	    "Atoms\002,4x,f18.6)";
    static char fmt_490[] = "(/,\002 Anlyt\002,6x,\002RMS Gradient over All "
	    "Atoms\002,4x,f16.4)";
    static char fmt_500[] = "(\002 Numer\002,6x,\002RMS Gradient over All At"
	    "oms\002,4x,f20.8)";
    static char fmt_510[] = "(\002 Numer\002,6x,\002RMS Gradient over All At"
	    "oms\002,4x,f18.6)";
    static char fmt_520[] = "(\002 Numer\002,6x,\002RMS Gradient over All At"
	    "oms\002,4x,f16.4)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void), gradient_(doublereal *, 
	    doublereal *);
    static logical doanalyt;
    static doublereal ntotnorm, f;
    static integer i__, j;
    static doublereal f0, ea0, eb0, ec0, ed0, eg0, em0, ep0, er0, es0, et0, 
	    ev0, ex0, old, eps, rms, eaa0, eba0, ecd0, eid0, elf0, ebt0, eub0,
	     eit0, eps0, ept0, ett0, ndea[75000]	/* was [3][25000] */, 
	    ndeb[75000]	/* was [3][25000] */, ndec[75000]	/* was [3][
	    25000] */, nded[75000]	/* was [3][25000] */, ndeg[75000]	
	    /* was [3][25000] */, ndem[75000]	/* was [3][25000] */, ndep[
	    75000]	/* was [3][25000] */, ndet[75000]	/* was [3][
	    25000] */, nder[75000]	/* was [3][25000] */, ndev[75000]	
	    /* was [3][25000] */, ndes[75000]	/* was [3][25000] */, ndex[
	    75000]	/* was [3][25000] */, etot;
    static integer next;
    static doublereal nrms, eopb0, eopd0, ndeaa[75000]	/* was [3][25000] */, 
	    ndeba[75000]	/* was [3][25000] */, ndecd[75000]	/* 
	    was [3][25000] */, ndeid[75000]	/* was [3][25000] */, ndelf[
	    75000]	/* was [3][25000] */, ndebt[75000]	/* was [3][
	    25000] */, ndeub[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int final_(void);
    static doublereal ndeit[75000]	/* was [3][25000] */, ndept[75000]	
	    /* was [3][25000] */, ndett[75000]	/* was [3][25000] */, detot[
	    75000]	/* was [3][25000] */;
    static logical exist, query;
    static doublereal ndeopb[75000]	/* was [3][25000] */, ndeopd[75000]	
	    /* was [3][25000] */;
    static char record[120];
    static doublereal denorm[25000];
    extern doublereal energy_(void);
    static doublereal ndetot[75000]	/* was [3][25000] */;
    static logical dofull;
    static char answer[1], string[120];
    extern /* Subroutine */ int upcase_(char *, ftnlen), getxyz_(void), 
	    initial_(void);
    static doublereal ndenorm[25000];
    static logical donumer;
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), gettext_(
	    char *, char *, integer *, ftnlen, ftnlen);
    static doublereal totnorm;

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_40, 0 };
    static icilist io___16 = { 1, string, 1, 0, 120, 1 };
    static cilist io___17 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___18 = { 1, 0, 0, fmt_70, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___92 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___93 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___94 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___95 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___96 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___97 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___98 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___99 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___100 = { 0, 0, 0, fmt_320, 0 };
    static cilist io___101 = { 0, 0, 0, fmt_330, 0 };
    static cilist io___105 = { 0, 0, 0, fmt_340, 0 };
    static cilist io___106 = { 0, 0, 0, fmt_350, 0 };
    static cilist io___107 = { 0, 0, 0, fmt_360, 0 };
    static cilist io___109 = { 0, 0, 0, fmt_370, 0 };
    static cilist io___110 = { 0, 0, 0, fmt_380, 0 };
    static cilist io___111 = { 0, 0, 0, fmt_390, 0 };
    static cilist io___112 = { 0, 0, 0, fmt_400, 0 };
    static cilist io___113 = { 0, 0, 0, fmt_410, 0 };
    static cilist io___114 = { 0, 0, 0, fmt_420, 0 };
    static cilist io___115 = { 0, 0, 0, fmt_430, 0 };
    static cilist io___116 = { 0, 0, 0, fmt_440, 0 };
    static cilist io___117 = { 0, 0, 0, fmt_450, 0 };
    static cilist io___118 = { 0, 0, 0, fmt_460, 0 };
    static cilist io___120 = { 0, 0, 0, fmt_470, 0 };
    static cilist io___121 = { 0, 0, 0, fmt_480, 0 };
    static cilist io___122 = { 0, 0, 0, fmt_490, 0 };
    static cilist io___124 = { 0, 0, 0, fmt_500, 0 };
    static cilist io___125 = { 0, 0, 0, fmt_510, 0 };
    static cilist io___126 = { 0, 0, 0, fmt_520, 0 };



#define dea_ref(a_1,a_2) deriv_1.dea[(a_2)*3 + a_1 - 4]
#define deb_ref(a_1,a_2) deriv_1.deb[(a_2)*3 + a_1 - 4]
#define dec_ref(a_1,a_2) deriv_1.dec[(a_2)*3 + a_1 - 4]
#define ded_ref(a_1,a_2) deriv_1.ded[(a_2)*3 + a_1 - 4]
#define deg_ref(a_1,a_2) deriv_1.deg[(a_2)*3 + a_1 - 4]
#define dem_ref(a_1,a_2) deriv_1.dem[(a_2)*3 + a_1 - 4]
#define dep_ref(a_1,a_2) deriv_1.dep[(a_2)*3 + a_1 - 4]
#define der_ref(a_1,a_2) deriv_1.der[(a_2)*3 + a_1 - 4]
#define des_ref(a_1,a_2) deriv_1.des[(a_2)*3 + a_1 - 4]
#define det_ref(a_1,a_2) deriv_1.det[(a_2)*3 + a_1 - 4]
#define dev_ref(a_1,a_2) deriv_1.dev[(a_2)*3 + a_1 - 4]
#define dex_ref(a_1,a_2) deriv_1.dex[(a_2)*3 + a_1 - 4]
#define deaa_ref(a_1,a_2) deriv_1.deaa[(a_2)*3 + a_1 - 4]
#define deba_ref(a_1,a_2) deriv_1.deba[(a_2)*3 + a_1 - 4]
#define decd_ref(a_1,a_2) deriv_1.decd[(a_2)*3 + a_1 - 4]
#define deid_ref(a_1,a_2) deriv_1.deid[(a_2)*3 + a_1 - 4]
#define ndea_ref(a_1,a_2) ndea[(a_2)*3 + a_1 - 4]
#define ndeb_ref(a_1,a_2) ndeb[(a_2)*3 + a_1 - 4]
#define ndec_ref(a_1,a_2) ndec[(a_2)*3 + a_1 - 4]
#define delf_ref(a_1,a_2) deriv_1.delf[(a_2)*3 + a_1 - 4]
#define nded_ref(a_1,a_2) nded[(a_2)*3 + a_1 - 4]
#define debt_ref(a_1,a_2) deriv_1.debt[(a_2)*3 + a_1 - 4]
#define deub_ref(a_1,a_2) deriv_1.deub[(a_2)*3 + a_1 - 4]
#define ndeg_ref(a_1,a_2) ndeg[(a_2)*3 + a_1 - 4]
#define ndem_ref(a_1,a_2) ndem[(a_2)*3 + a_1 - 4]
#define deit_ref(a_1,a_2) deriv_1.deit[(a_2)*3 + a_1 - 4]
#define ndep_ref(a_1,a_2) ndep[(a_2)*3 + a_1 - 4]
#define ndet_ref(a_1,a_2) ndet[(a_2)*3 + a_1 - 4]
#define nder_ref(a_1,a_2) nder[(a_2)*3 + a_1 - 4]
#define dept_ref(a_1,a_2) deriv_1.dept[(a_2)*3 + a_1 - 4]
#define ndev_ref(a_1,a_2) ndev[(a_2)*3 + a_1 - 4]
#define ndes_ref(a_1,a_2) ndes[(a_2)*3 + a_1 - 4]
#define ndex_ref(a_1,a_2) ndex[(a_2)*3 + a_1 - 4]
#define dett_ref(a_1,a_2) deriv_1.dett[(a_2)*3 + a_1 - 4]
#define ndeaa_ref(a_1,a_2) ndeaa[(a_2)*3 + a_1 - 4]
#define ndeba_ref(a_1,a_2) ndeba[(a_2)*3 + a_1 - 4]
#define ndecd_ref(a_1,a_2) ndecd[(a_2)*3 + a_1 - 4]
#define ndeid_ref(a_1,a_2) ndeid[(a_2)*3 + a_1 - 4]
#define ndelf_ref(a_1,a_2) ndelf[(a_2)*3 + a_1 - 4]
#define deopb_ref(a_1,a_2) deriv_1.deopb[(a_2)*3 + a_1 - 4]
#define deopd_ref(a_1,a_2) deriv_1.deopd[(a_2)*3 + a_1 - 4]
#define ndebt_ref(a_1,a_2) ndebt[(a_2)*3 + a_1 - 4]
#define ndeub_ref(a_1,a_2) ndeub[(a_2)*3 + a_1 - 4]
#define ndeit_ref(a_1,a_2) ndeit[(a_2)*3 + a_1 - 4]
#define ndept_ref(a_1,a_2) ndept[(a_2)*3 + a_1 - 4]
#define ndett_ref(a_1,a_2) ndett[(a_2)*3 + a_1 - 4]
#define detot_ref(a_1,a_2) detot[(a_2)*3 + a_1 - 4]
#define ndeopb_ref(a_1,a_2) ndeopb[(a_2)*3 + a_1 - 4]
#define ndeopd_ref(a_1,a_2) ndeopd[(a_2)*3 + a_1 - 4]
#define ndetot_ref(a_1,a_2) ndetot[(a_2)*3 + a_1 - 4]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  atoms.i  --  number, position and type of current atoms  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     x       current x-coordinate for each atom in the system */
/*     y       current y-coordinate for each atom in the system */
/*     z       current z-coordinate for each atom in the system */
/*     n       total number of atoms in the current system */
/*     type    atom type number for each atom in the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  deriv.i  --  Cartesian coordinate derivative components  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     desum   total energy Cartesian coordinate derivatives */
/*     deb     bond stretch Cartesian coordinate derivatives */
/*     dea     angle bend Cartesian coordinate derivatives */
/*     deba    stretch-bend Cartesian coordinate derivatives */
/*     deub    Urey-Bradley Cartesian coordinate derivatives */
/*     deaa    angle-angle Cartesian coordinate derivatives */
/*     deopb   out-of-plane bend Cartesian coordinate derivatives */
/*     deopd   out-of-plane distance Cartesian coordinate derivatives */
/*     deid    improper dihedral Cartesian coordinate derivatives */
/*     deit    improper torsion Cartesian coordinate derivatives */
/*     det     torsional Cartesian coordinate derivatives */
/*     dept    pi-orbital torsion Cartesian coordinate derivatives */
/*     debt    stretch-torsion Cartesian coordinate derivatives */
/*     dett    torsion-torsion Cartesian coordinate derivatives */
/*     dev     van der Waals Cartesian coordinate derivatives */
/*     dec     charge-charge Cartesian coordinate derivatives */
/*     decd    charge-dipole Cartesian coordinate derivatives */
/*     ded     dipole-dipole Cartesian coordinate derivatives */
/*     dem     multipole Cartesian coordinate derivatives */
/*     dep     polarization Cartesian coordinate derivatives */
/*     der     reaction field Cartesian coordinate derivatives */
/*     des     solvation Cartesian coordinate derivatives */
/*     delf    metal ligand field Cartesian coordinate derivatives */
/*     deg     geometric restraint Cartesian coordinate derivatives */
/*     dex     extra energy term Cartesian coordinate derivatives */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  energi.i  --  individual potential energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     esum   total potential energy of the system */
/*     eb     bond stretch potential energy of the system */
/*     ea     angle bend potential energy of the system */
/*     eba    stretch-bend potential energy of the system */
/*     eub    Urey-Bradley potential energy of the system */
/*     eaa    angle-angle potential energy of the system */
/*     eopb   out-of-plane bend potential energy of the system */
/*     eopd   out-of-plane distance potential energy of the system */
/*     eid    improper dihedral potential energy of the system */
/*     eit    improper torsion potential energy of the system */
/*     et     torsional potential energy of the system */
/*     ept    pi-orbital torsion potential energy of the system */
/*     ebt    stretch-torsion potential energy of the system */
/*     ett    torsion-torsion potential energy of the system */
/*     ev     van der Waals potential energy of the system */
/*     ec     charge-charge potential energy of the system */
/*     ecd    charge-dipole potential energy of the system */
/*     ed     dipole-dipole potential energy of the system */
/*     em     atomic multipole potential energy of the system */
/*     ep     polarization potential energy of the system */
/*     er     reaction field potential energy of the system */
/*     es     solvation potential energy of the system */
/*     elf    metal ligand field potential energy of the system */
/*     eg     geometric restraint potential energy of the system */
/*     ex     extra term potential energy of the system */




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  inter.i  --  sum of intermolecular energy components  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     einter   total intermolecular potential energy */




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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  solute.i  --  parameters for continuum solvation models  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     rsolv     atomic radius of each atom for continuum solvation */
/*     asolv     atomic surface area solvation parameters */
/*     rborn     Born radius of each atom for GB/SA solvation */
/*     drb       solvation derivatives with respect to Born radii */
/*     drbp      GK polarization derivatives with respect to Born radii */
/*     drobc     chain rule term for Onufriev-Bashford-Case radii */
/*     doffset   dielectric offset to continuum solvation atomic radii */
/*     p1        single-atom scale factor for analytical Still radii */
/*     p2        1-2 interaction scale factor for analytical Still radii */
/*     p3        1-3 interaction scale factor for analytical Still radii */
/*     p4        nonbonded scale factor for analytical Still radii */
/*     p5        soft cutoff parameter for analytical Still radii */
/*     gpol      polarization self-energy values for each atom */
/*     shct      overlap scale factors for Hawkins-Cramer-Truhlar radii */
/*     aobc      alpha values for Onufriev-Bashford-Case radii */
/*     bobc      beta values for Onufriev-Bashford-Case radii */
/*     gobc      gamma values for Onufriev-Bashford-Case radii */
/*     vsolv     atomic volume of each atom for use with ACE */
/*     wace      "omega" values for atom class pairs for use with ACE */
/*     s2ace     "sigma^2" values for atom class pairs for use with ACE */
/*     uace      "mu" values for atom class pairs for use with ACE */
/*     solvtyp   type of continuum solvation energy model in use */
/*     borntyp   method to be used for the Born radius computation */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     set up the structure and mechanics calculation */

    initial_();
    getxyz_();
    mechanic_();

/*     decide whether to do an analytical gradient calculation */

    doanalyt = TRUE_;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___5.ciunit = iounit_1.iout;
	s_wsfe(&io___5);
	e_wsfe();
	io___6.ciunit = iounit_1.input;
	s_rsfe(&io___6);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'N') {
	doanalyt = FALSE_;
    }

/*     decide whether to do a numerical gradient calculation */

    donumer = TRUE_;
    nextarg_(answer, &exist, (ftnlen)1);
    if (! exist) {
	io___10.ciunit = iounit_1.iout;
	s_wsfe(&io___10);
	e_wsfe();
	io___11.ciunit = iounit_1.input;
	s_rsfe(&io___11);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    }
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'N') {
	donumer = FALSE_;
    }

/*     get the stepsize for numerical gradient calculation */

    if (donumer) {
	eps = -1.;
	eps0 = 1e-5;
	query = TRUE_;
	nextarg_(string, &exist, (ftnlen)120);
	if (exist) {
	    i__1 = s_rsli(&io___16);
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&eps, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L50;
	    }
	    query = FALSE_;
	}
L50:
	if (query) {
	    io___17.ciunit = iounit_1.iout;
	    s_wsfe(&io___17);
	    do_fio(&c__1, (char *)&eps0, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    io___18.ciunit = iounit_1.input;
	    i__1 = s_rsfe(&io___18);
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = do_fio(&c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
	    if (i__1 != 0) {
		goto L50;
	    }
	    i__1 = e_rsfe();
	    if (i__1 != 0) {
		goto L50;
	    }
	}
	if (eps <= 0.) {
	    eps = eps0;
	}
    }

/*     decide whether to output results by gradient component */

    dofull = TRUE_;
    if (atoms_1.n > 100) {
	dofull = FALSE_;
	nextarg_(answer, &exist, (ftnlen)1);
	if (! exist) {
	    io___20.ciunit = iounit_1.iout;
	    s_wsfe(&io___20);
	    e_wsfe();
	    io___21.ciunit = iounit_1.input;
	    s_rsfe(&io___21);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    next = 1;
	    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
	}
	upcase_(answer, (ftnlen)1);
	if (*(unsigned char *)answer == 'Y') {
	    dofull = TRUE_;
	}
    }

/*     compute the analytical gradient components */

    if (doanalyt) {
	gradient_(&etot, detot);
    }

/*     print the total potential energy of the system */

    if (doanalyt) {
	if (inform_1.digits >= 8) {
	    io___24.ciunit = iounit_1.iout;
	    s_wsfe(&io___24);
	    do_fio(&c__1, (char *)&etot, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (inform_1.digits >= 6) {
	    io___25.ciunit = iounit_1.iout;
	    s_wsfe(&io___25);
	    do_fio(&c__1, (char *)&etot, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___26.ciunit = iounit_1.iout;
	    s_wsfe(&io___26);
	    do_fio(&c__1, (char *)&etot, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}

/*     print the energy breakdown over individual components */

	io___27.ciunit = iounit_1.iout;
	s_wsfe(&io___27);
	e_wsfe();
	if (inform_1.digits >= 8) {
	    io___28.ciunit = iounit_1.iout;
	    s_wsfe(&io___28);
	    e_wsfe();
	    io___29.ciunit = iounit_1.iout;
	    s_wsfe(&io___29);
	    do_fio(&c__1, (char *)&energi_1.eb, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ea, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eba, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eub, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eaa, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eopb, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eopd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eid, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eit, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.et, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ept, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ebt, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ett, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ev, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ec, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ecd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ed, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.em, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ep, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.er, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.es, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.elf, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eg, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ex, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (inform_1.digits >= 6) {
	    io___30.ciunit = iounit_1.iout;
	    s_wsfe(&io___30);
	    e_wsfe();
	    io___31.ciunit = iounit_1.iout;
	    s_wsfe(&io___31);
	    do_fio(&c__1, (char *)&energi_1.eb, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ea, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eba, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eub, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eaa, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eopb, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eopd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eid, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eit, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.et, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ept, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ebt, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ett, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ev, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ec, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ecd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ed, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.em, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ep, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.er, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.es, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.elf, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eg, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ex, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___32.ciunit = iounit_1.iout;
	    s_wsfe(&io___32);
	    e_wsfe();
	    io___33.ciunit = iounit_1.iout;
	    s_wsfe(&io___33);
	    do_fio(&c__1, (char *)&energi_1.eb, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ea, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eba, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eub, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eaa, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eopb, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eopd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eid, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eit, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.et, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ept, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ebt, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ett, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ev, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ec, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ecd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ed, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.em, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ep, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.er, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.es, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.elf, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.eg, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&energi_1.ex, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     print a header for the gradients of individual potentials */

    if (dofull) {
	io___34.ciunit = iounit_1.iout;
	s_wsfe(&io___34);
	e_wsfe();
	if (inform_1.digits >= 8) {
	    io___35.ciunit = iounit_1.iout;
	    s_wsfe(&io___35);
	    e_wsfe();
	} else if (inform_1.digits >= 6) {
	    io___36.ciunit = iounit_1.iout;
	    s_wsfe(&io___36);
	    e_wsfe();
	} else {
	    io___37.ciunit = iounit_1.iout;
	    s_wsfe(&io___37);
	    e_wsfe();
	}
    }

/*     get the Cartesian component two-sided numerical gradients */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (donumer && usage_1.use[i__ - 1]) {
	    for (j = 1; j <= 3; ++j) {
		if (j == 1) {
		    old = atoms_1.x[i__ - 1];
		    atoms_1.x[i__ - 1] -= eps * .5;
		} else if (j == 2) {
		    old = atoms_1.y[i__ - 1];
		    atoms_1.y[i__ - 1] -= eps * .5;
		} else if (j == 3) {
		    old = atoms_1.z__[i__ - 1];
		    atoms_1.z__[i__ - 1] -= eps * .5;
		}
		f0 = energy_();
		eb0 = energi_1.eb;
		ea0 = energi_1.ea;
		eba0 = energi_1.eba;
		eub0 = energi_1.eub;
		eaa0 = energi_1.eaa;
		eopb0 = energi_1.eopb;
		eopd0 = energi_1.eopd;
		eid0 = energi_1.eid;
		eit0 = energi_1.eit;
		et0 = energi_1.et;
		ept0 = energi_1.ept;
		ebt0 = energi_1.ebt;
		ett0 = energi_1.ett;
		ev0 = energi_1.ev;
		ec0 = energi_1.ec;
		ecd0 = energi_1.ecd;
		ed0 = energi_1.ed;
		em0 = energi_1.em;
		ep0 = energi_1.ep;
		er0 = energi_1.er;
		es0 = energi_1.es;
		elf0 = energi_1.elf;
		eg0 = energi_1.eg;
		ex0 = energi_1.ex;
		if (j == 1) {
		    atoms_1.x[i__ - 1] += eps;
		} else if (j == 2) {
		    atoms_1.y[i__ - 1] += eps;
		} else if (j == 3) {
		    atoms_1.z__[i__ - 1] += eps;
		}
		f = energy_();
		if (j == 1) {
		    atoms_1.x[i__ - 1] = old;
		} else if (j == 2) {
		    atoms_1.y[i__ - 1] = old;
		} else if (j == 3) {
		    atoms_1.z__[i__ - 1] = old;
		}
		ndeb_ref(j, i__) = (energi_1.eb - eb0) / eps;
		ndea_ref(j, i__) = (energi_1.ea - ea0) / eps;
		ndeba_ref(j, i__) = (energi_1.eba - eba0) / eps;
		ndeub_ref(j, i__) = (energi_1.eub - eub0) / eps;
		ndeaa_ref(j, i__) = (energi_1.eaa - eaa0) / eps;
		ndeopb_ref(j, i__) = (energi_1.eopb - eopb0) / eps;
		ndeopd_ref(j, i__) = (energi_1.eopd - eopd0) / eps;
		ndeid_ref(j, i__) = (energi_1.eid - eid0) / eps;
		ndeit_ref(j, i__) = (energi_1.eit - eit0) / eps;
		ndet_ref(j, i__) = (energi_1.et - et0) / eps;
		ndept_ref(j, i__) = (energi_1.ept - ept0) / eps;
		ndebt_ref(j, i__) = (energi_1.ebt - ebt0) / eps;
		ndett_ref(j, i__) = (energi_1.ett - ett0) / eps;
		ndev_ref(j, i__) = (energi_1.ev - ev0) / eps;
		ndec_ref(j, i__) = (energi_1.ec - ec0) / eps;
		ndecd_ref(j, i__) = (energi_1.ecd - ecd0) / eps;
		nded_ref(j, i__) = (energi_1.ed - ed0) / eps;
		ndem_ref(j, i__) = (energi_1.em - em0) / eps;
		ndep_ref(j, i__) = (energi_1.ep - ep0) / eps;
		nder_ref(j, i__) = (energi_1.er - er0) / eps;
		ndes_ref(j, i__) = (energi_1.es - es0) / eps;
		ndelf_ref(j, i__) = (energi_1.elf - elf0) / eps;
		ndeg_ref(j, i__) = (energi_1.eg - eg0) / eps;
		ndex_ref(j, i__) = (energi_1.ex - ex0) / eps;
		ndetot_ref(j, i__) = (f - f0) / eps;
	    }
	}

/*     print analytical gradients of each energy term for each atom */

	if (dofull && usage_1.use[i__ - 1]) {
	    for (j = 1; j <= 3; ++j) {
		if (doanalyt) {
		    if (inform_1.digits >= 8) {
			io___92.ciunit = iounit_1.iout;
			s_wsfe(&io___92);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&deb_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dea_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deba_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deub_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			do_fio(&c__1, (char *)&deaa_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deopb_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deopd_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deid_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deit_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&det_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dept_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&debt_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dett_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dev_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dec_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&decd_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ded_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dem_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dep_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&der_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&des_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&delf_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deg_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dex_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			e_wsfe();
		    } else if (inform_1.digits >= 6) {
			io___93.ciunit = iounit_1.iout;
			s_wsfe(&io___93);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&deb_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dea_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deba_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deub_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deaa_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			do_fio(&c__1, (char *)&deopb_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deopd_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deid_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deit_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&det_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dept_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&debt_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dett_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dev_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dec_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&decd_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ded_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dem_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dep_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&der_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&des_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&delf_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deg_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dex_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			e_wsfe();
		    } else {
			io___94.ciunit = iounit_1.iout;
			s_wsfe(&io___94);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&deb_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dea_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deba_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deub_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deaa_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deopb_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			do_fio(&c__1, (char *)&deopd_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deid_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deit_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&det_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dept_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&debt_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dett_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dev_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dec_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&decd_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ded_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dem_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dep_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&der_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&des_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&delf_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&deg_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&dex_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			e_wsfe();
		    }
		}

/*     print numerical gradients of each energy term for each atom */

		if (donumer) {
		    if (inform_1.digits >= 8) {
			io___95.ciunit = iounit_1.iout;
			s_wsfe(&io___95);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&ndeb_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndea_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeba_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeub_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			do_fio(&c__1, (char *)&ndeaa_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeopb_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeopd_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeid_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeit_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndet_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndept_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndebt_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndett_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndev_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndec_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndecd_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&nded_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndem_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndep_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&nder_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndes_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndelf_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeg_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndex_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			e_wsfe();
		    } else if (inform_1.digits >= 6) {
			io___96.ciunit = iounit_1.iout;
			s_wsfe(&io___96);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&ndeb_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndea_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeba_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeub_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeaa_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			do_fio(&c__1, (char *)&ndeopb_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeopd_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeid_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeit_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndet_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndept_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndebt_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndett_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndev_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndec_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndecd_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&nded_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndem_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndep_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&nder_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndes_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndelf_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeg_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndex_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			e_wsfe();
		    } else {
			io___97.ciunit = iounit_1.iout;
			s_wsfe(&io___97);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&ndeb_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndea_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeba_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeub_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeaa_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeopb_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, axis + (j - 1), (ftnlen)1);
			do_fio(&c__1, (char *)&ndeopd_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeid_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeit_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndet_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndept_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndebt_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndett_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndev_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndec_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndecd_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&nded_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndem_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndep_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&nder_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndes_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndelf_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndeg_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&ndex_ref(j, i__), (ftnlen)
				sizeof(doublereal));
			e_wsfe();
		    }
		}
	    }
	}
    }

/*     print the total gradient components for each atom */

    io___98.ciunit = iounit_1.iout;
    s_wsfe(&io___98);
    e_wsfe();
    if (inform_1.digits >= 8) {
	io___99.ciunit = iounit_1.iout;
	s_wsfe(&io___99);
	e_wsfe();
    } else if (inform_1.digits >= 6) {
	io___100.ciunit = iounit_1.iout;
	s_wsfe(&io___100);
	e_wsfe();
    } else {
	io___101.ciunit = iounit_1.iout;
	s_wsfe(&io___101);
	e_wsfe();
    }
    totnorm = 0.;
    ntotnorm = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (doanalyt && usage_1.use[i__ - 1]) {
/* Computing 2nd power */
	    d__1 = detot_ref(1, i__);
/* Computing 2nd power */
	    d__2 = detot_ref(2, i__);
/* Computing 2nd power */
	    d__3 = detot_ref(3, i__);
	    denorm[i__ - 1] = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	    totnorm += denorm[i__ - 1];
	    denorm[i__ - 1] = sqrt(denorm[i__ - 1]);
	    if (inform_1.digits >= 8) {
		io___105.ciunit = iounit_1.iout;
		s_wsfe(&io___105);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&detot_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		do_fio(&c__1, (char *)&denorm[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else if (inform_1.digits >= 6) {
		io___106.ciunit = iounit_1.iout;
		s_wsfe(&io___106);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&detot_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		do_fio(&c__1, (char *)&denorm[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___107.ciunit = iounit_1.iout;
		s_wsfe(&io___107);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&detot_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		do_fio(&c__1, (char *)&denorm[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
	if (donumer && usage_1.use[i__ - 1]) {
/* Computing 2nd power */
	    d__1 = ndetot_ref(1, i__);
/* Computing 2nd power */
	    d__2 = ndetot_ref(2, i__);
/* Computing 2nd power */
	    d__3 = ndetot_ref(3, i__);
	    ndenorm[i__ - 1] = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	    ntotnorm += ndenorm[i__ - 1];
	    ndenorm[i__ - 1] = sqrt(ndenorm[i__ - 1]);
	    if (inform_1.digits >= 8) {
		io___109.ciunit = iounit_1.iout;
		s_wsfe(&io___109);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&ndetot_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		do_fio(&c__1, (char *)&ndenorm[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else if (inform_1.digits >= 6) {
		io___110.ciunit = iounit_1.iout;
		s_wsfe(&io___110);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&ndetot_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		do_fio(&c__1, (char *)&ndenorm[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___111.ciunit = iounit_1.iout;
		s_wsfe(&io___111);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		for (j = 1; j <= 3; ++j) {
		    do_fio(&c__1, (char *)&ndetot_ref(j, i__), (ftnlen)sizeof(
			    doublereal));
		}
		do_fio(&c__1, (char *)&ndenorm[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     print the total norm for the analytical gradient */

    io___112.ciunit = iounit_1.iout;
    s_wsfe(&io___112);
    e_wsfe();
    if (doanalyt) {
	totnorm = sqrt(totnorm);
	if (inform_1.digits >= 8) {
	    io___113.ciunit = iounit_1.iout;
	    s_wsfe(&io___113);
	    do_fio(&c__1, (char *)&totnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (inform_1.digits >= 6) {
	    io___114.ciunit = iounit_1.iout;
	    s_wsfe(&io___114);
	    do_fio(&c__1, (char *)&totnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___115.ciunit = iounit_1.iout;
	    s_wsfe(&io___115);
	    do_fio(&c__1, (char *)&totnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     print the total norm for the numerical gradient */

    if (donumer) {
	ntotnorm = sqrt(ntotnorm);
	if (inform_1.digits >= 8) {
	    io___116.ciunit = iounit_1.iout;
	    s_wsfe(&io___116);
	    do_fio(&c__1, (char *)&ntotnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (inform_1.digits >= 6) {
	    io___117.ciunit = iounit_1.iout;
	    s_wsfe(&io___117);
	    do_fio(&c__1, (char *)&ntotnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___118.ciunit = iounit_1.iout;
	    s_wsfe(&io___118);
	    do_fio(&c__1, (char *)&ntotnorm, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     print the rms per atom norm for the analytical gradient */

    if (doanalyt) {
	rms = totnorm / sqrt((doublereal) usage_1.nuse);
	if (inform_1.digits >= 8) {
	    io___120.ciunit = iounit_1.iout;
	    s_wsfe(&io___120);
	    do_fio(&c__1, (char *)&rms, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (inform_1.digits >= 6) {
	    io___121.ciunit = iounit_1.iout;
	    s_wsfe(&io___121);
	    do_fio(&c__1, (char *)&rms, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___122.ciunit = iounit_1.iout;
	    s_wsfe(&io___122);
	    do_fio(&c__1, (char *)&rms, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     print the rms per atom norm for the numerical gradient */

    if (donumer) {
	nrms = ntotnorm / sqrt((doublereal) usage_1.nuse);
	if (inform_1.digits >= 8) {
	    io___124.ciunit = iounit_1.iout;
	    s_wsfe(&io___124);
	    do_fio(&c__1, (char *)&nrms, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (inform_1.digits >= 6) {
	    io___125.ciunit = iounit_1.iout;
	    s_wsfe(&io___125);
	    do_fio(&c__1, (char *)&nrms, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___126.ciunit = iounit_1.iout;
	    s_wsfe(&io___126);
	    do_fio(&c__1, (char *)&nrms, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef ndetot_ref
#undef ndeopd_ref
#undef ndeopb_ref
#undef detot_ref
#undef ndett_ref
#undef ndept_ref
#undef ndeit_ref
#undef ndeub_ref
#undef ndebt_ref
#undef deopd_ref
#undef deopb_ref
#undef ndelf_ref
#undef ndeid_ref
#undef ndecd_ref
#undef ndeba_ref
#undef ndeaa_ref
#undef dett_ref
#undef ndex_ref
#undef ndes_ref
#undef ndev_ref
#undef dept_ref
#undef nder_ref
#undef ndet_ref
#undef ndep_ref
#undef deit_ref
#undef ndem_ref
#undef ndeg_ref
#undef deub_ref
#undef debt_ref
#undef nded_ref
#undef delf_ref
#undef ndec_ref
#undef ndeb_ref
#undef ndea_ref
#undef deid_ref
#undef decd_ref
#undef deba_ref
#undef deaa_ref
#undef dex_ref
#undef dev_ref
#undef det_ref
#undef des_ref
#undef der_ref
#undef dep_ref
#undef dem_ref
#undef deg_ref
#undef ded_ref
#undef dec_ref
#undef deb_ref
#undef dea_ref


/* Main program alias */ int testgrad_ () { MAIN__ (); return 0; }
