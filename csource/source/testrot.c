/* testrot.f -- translated by f2c (version 20050501).
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
    doublereal tesum[1000], teb[1000], tea[1000], teba[1000], teub[1000], 
	    teaa[1000], teopb[1000], teopd[1000], teid[1000], teit[1000], tet[
	    1000], tept[1000], tebt[1000], tett[1000], tev[1000], tec[1000], 
	    tecd[1000], ted[1000], tem[1000], tep[1000], ter[1000], tes[1000],
	     telf[1000], teg[1000], tex[1000];
} domega_;

#define domega_1 domega_

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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    doublereal dihed[1000];
    integer nomega, iomega[2000]	/* was [2][1000] */, zline[1000];
} omega_;

#define omega_1 omega_

struct {
    doublereal zbond[25000], zang[25000], ztors[25000];
    integer iz[100000]	/* was [4][25000] */;
} zcoord_;

#define zcoord_1 zcoord_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  program testrot  --  derivative test; torsional version  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "testrot" computes and compares the analytical and numerical */
/*     gradient vectors of the potential energy function with respect */
/*     to rotatable torsional angles */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_20[] = "(/,\002 Enter Numerical Derivative Stepsize ["
	    "\002,d7.1,\002 Deg] :  \002,$)";
    static char fmt_30[] = "(f20.0)";
    static char fmt_40[] = "(/,\002 Total Potential Energy :\002,8x,f20.8"
	    ",\002 Kcal/mole\002)";
    static char fmt_50[] = "(/,\002 Total Potential Energy :\002,8x,f18.6"
	    ",\002 Kcal/mole\002)";
    static char fmt_60[] = "(/,\002 Total Potential Energy :\002,8x,f16.4"
	    ",\002 Kcal/mole\002)";
    static char fmt_70[] = "(/,\002 Potential Energy Breakdown by Individua"
	    "l\002,\002 Components :\002)";
    static char fmt_80[] = "(/,\002  Energy\002,7x,\002EB\002,14x,\002EA\002"
	    ",14x,\002EBA\002,13x,\002EUB\002,/,\002  Terms\002,8x,\002EAA"
	    "\002,13x,\002EOPB\002,12x,\002EOPD\002,12x,\002EID\002,/,15x,"
	    "\002EIT\002,13x,\002ET\002,14x,\002EPT\002,13x,\002EBT\002,/,15x,"
	    "\002ETT\002,13x,\002EV\002,14x,\002EC\002,14x,\002ECD\002,/,15x"
	    ",\002ED\002,14x,\002EM\002,14x,\002EP\002,14x,\002ER\002,/,15x"
	    ",\002ES\002,14x,\002ELF\002,13x,\002EG\002,14x,\002EX\002)";
    static char fmt_90[] = "(/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8"
	    ",/,6x,4f16.8,/,6x,4f16.8)";
    static char fmt_100[] = "(/,\002  Energy\002,6x,\002EB\002,12x,\002EA"
	    "\002,12x,\002EBA\002,11x,\002EUB\002,11x,\002EAA\002,/,\002  Ter"
	    "ms\002,7x,\002EOPB\002,10x,\002EOPD\002,10x,\002EID\002,11x,\002"
	    "EIT\002,11x,\002ET\002,/,14x,\002EPT\002,11x,\002EBT\002,11x,"
	    "\002ETT\002,11x,\002EV\002,12x,\002EC\002,/,14x,\002ECD\002,11x"
	    ",\002ED\002,12x,\002EM\002,12x,\002EP\002,12x,\002ER\002,/,14x"
	    ",\002ES\002,12x,\002ELF\002,11x,\002EG\002,12x,\002EX\002)";
    static char fmt_110[] = "(/,6x,5f14.6,/,6x,5f14.6,/,6x,5f14.6,/,6x,5f14."
	    "6,/,6x,4f14.6)";
    static char fmt_120[] = "(/,\002  Energy\002,6x,\002EB\002,10x,\002EA"
	    "\002,10x,\002EBA\002,9x,\002EUB\002,9x,\002EAA\002,9x,\002EOP"
	    "B\002,/,\002  Terms\002,7x,\002EOPD\002,8x,\002EID\002,9x,\002EIT"
	    "\002,9x,\002ET\002,10x,\002EPT\002,9x,\002EBT\002,/,14x,\002ET"
	    "T\002,9x,\002EV\002,10x,\002EC\002,10x,\002ECD\002,9x,\002ED\002"
	    ",10x,\002EM\002,/,14x,\002EP\002,10x,\002ER\002,10x,\002ES\002,1"
	    "0x,\002ELF\002,9x,\002EG\002,10x,\002EX\002)";
    static char fmt_130[] = "(/,6x,6f12.4,/,6x,6f12.4,/,6x,6f12.4,/,6x,6f12."
	    "4)";
    static char fmt_140[] = "(/,\002 Torsional Gradient Breakdown by Individ"
	    "ual\002,\002 Components :\002)";
    static char fmt_150[] = "(/,2x,\002Atom\002,9x,\002d EB\002,12x,\002d E"
	    "A\002,12x,\002d EBA\002,11x,\002d EUB\002,/,2x,\002Axis\002,9x"
	    ",\002d EAA\002,11x,\002d EOPB\002,10x,\002d EOPD\002,10x,\002d E"
	    "ID\002,/,2x,\002Type\002,9x,\002d EIT\002,11x,\002d ET\002,12x"
	    ",\002d EPT\002,11x,\002d EBT\002,/,15x,\002d ETT\002,11x,\002d EV"
	    "\002,12x,\002d EC\002,12x,\002d ECD\002,/,15x,\002d ED\002,12x"
	    ",\002d EM\002,12x,\002d EP\002,12x,\002d ER\002,/,15x,\002d E"
	    "S\002,12x,\002d ELF\002,11x,\002d EG\002,12x,\002d EX\002)";
    static char fmt_160[] = "(/,2x,\002Atom\002,8x,\002d EB\002,10x,\002d E"
	    "A\002,10x,\002d EBA\002,9x,\002d EUB\002,9x,\002d EAA\002,/,2x"
	    ",\002Axis\002,8x,\002d EOPB\002,8x,\002d EOPD\002,8x,\002d EI"
	    "D\002,9x,\002d EIT\002,9x,\002d ET\002,/,2x,\002Type\002,8x,\002"
	    "d EPT\002,9x,\002d EBT\002,9x,\002d ETT\002,9x,\002d EV\002,10x"
	    ",\002d EC\002,/,14x,\002d ECD\002,9x,\002d ED\002,10x,\002d E"
	    "M\002,10x,\002d EP\002,10x,\002d ER\002,/,14x,\002d ES\002,10x"
	    ",\002d ELF\002,9x,\002d EG\002,10x,\002d EX\002)";
    static char fmt_170[] = "(/,2x,\002Atom\002,6x,\002d EB\002,8x,\002d E"
	    "A\002,8x,\002d EBA\002,7x,\002d EUB\002,7x,\002d EAA\002,7x,\002"
	    "d EOPB\002,/,2x,\002Axis\002,6x,\002d EOPD\002,6x,\002d EID\002,"
	    "7x,\002d EIT\002,7x,\002d ET\002,8x,\002d EPT\002,7x,\002d EB"
	    "T\002,/,2x,\002Type\002,6x,\002d ETT\002,7x,\002d EV\002,8x,\002"
	    "d EC\002,8x,\002d ECD\002,7x,\002d ED\002,8x,\002d EM\002,/,12x"
	    ",\002d EP\002,8x,\002d ER\002,8x,\002d ES\002,8x,\002d ELF\002,7"
	    "x,\002d EG\002,8x,\002d EX\002)";
    static char fmt_180[] = "(/,i6,4f16.8,/,i6,4f16.8,/,\002 Anlyt\002,4f16."
	    "8,/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8)";
    static char fmt_190[] = "(/,i6,5f14.6,/,i6,5f14.6,/,\002 Anlyt\002,5f14."
	    "6,/,6x,5f14.6,/,6x,4f14.6)";
    static char fmt_200[] = "(/,i6,6f12.4,/,i6,6f12.4,/,\002 Anlyt\002,6f12."
	    "4,/,6x,6f12.4)";
    static char fmt_210[] = "(/,i6,4f16.8,/,i6,4f16.8,/,\002 Numer\002,4f16."
	    "8,/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8)";
    static char fmt_220[] = "(/,i6,5f14.6,/,i6,5f14.6,/,\002 Numer\002,5f14."
	    "6,/,6x,5f14.6,/,6x,4f14.6)";
    static char fmt_230[] = "(/,i6,6f12.4,/,i6,6f12.4,/,\002 Numer\002,6f12."
	    "4,/,6x,6f12.4)";
    static char fmt_240[] = "(//,5x,\002Torsion\002,19x,\002Anlyt Deriv\002,"
	    "9x,\002Numer Deriv\002,/)";
    static char fmt_250[] = "(//,5x,\002Torsion\002,18x,\002Anlyt Deriv\002,"
	    "7x,\002Numer Deriv\002,/)";
    static char fmt_260[] = "(//,5x,\002Torsion\002,17x,\002Anlyt Deriv\002,"
	    "5x,\002Numer Deriv\002,/)";
    static char fmt_270[] = "(1x,i5,\002-\002,i5,10x,2f20.8)";
    static char fmt_280[] = "(1x,i5,\002-\002,i5,10x,2f18.6)";
    static char fmt_290[] = "(1x,i5,\002-\002,i5,10x,2f16.4)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12, d__13, d__14, d__15, d__16, d__17, d__18, d__19, 
	    d__20, d__21, d__22, d__23, d__24;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_wsfe(cilist *), do_fio(integer *, char *, ftnlen),
	     e_wsfe(void), s_rsfe(cilist *), e_rsfe(void);

    /* Local variables */
    extern /* Subroutine */ int mechanic_(void);
    static doublereal e;
    static integer i__;
    static doublereal e0, ea0, eb0, ec0, ed0, eg0, em0, ep0, er0, es0, et0, 
	    ev0, ex0, eps, eaa0, eba0, ecd0, eid0, elf0, ebt0, eub0, eit0, 
	    ept0, ett0, etot, eopb0, eopd0, delta;
    extern /* Subroutine */ int final_(void);
    static logical exist, query;
    static doublereal delta0;
    extern doublereal energy_(void);
    static doublereal derivs[1000], nderiv[1000];
    static char string[120];
    extern /* Subroutine */ int getint_(void), initial_(void), gradrot_(
	    doublereal *, doublereal *), nextarg_(char *, logical *, ftnlen), 
	    initrot_(void), makexyz_(void);

    /* Fortran I/O blocks */
    static icilist io___6 = { 1, string, 1, 0, 120, 1 };
    static cilist io___7 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___8 = { 1, 0, 0, fmt_30, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___54 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___55 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___56 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___57 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___58 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___59 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___60 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___64 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___65 = { 0, 0, 0, fmt_290, 0 };



#define iomega_ref(a_1,a_2) omega_1.iomega[(a_2)*2 + a_1 - 3]



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

/*     ######################################################### */
/*     ##                                                     ## */
/*     ##  domega.i  --  derivative components over torsions  ## */
/*     ##                                                     ## */
/*     ######################################################### */


/*     tesum   total energy derivatives over torsions */
/*     teb     bond stretch derivatives over torsions */
/*     tea     angle bend derivatives over torsions */
/*     teba    stretch-bend derivatives over torsions */
/*     teub    Urey-Bradley derivatives over torsions */
/*     teaa    angle-angle derivatives over torsions */
/*     teopb   out-of-plane bend derivatives over torsions */
/*     teopd   out-of-plane distance derivatives over torsions */
/*     teid    improper dihedral derivatives over torsions */
/*     teit    improper torsion derivatives over torsions */
/*     tet     torsional derivatives over torsions */
/*     tept    pi-orbital torsion derivatives over torsions */
/*     tebt    stretch-torsion derivatives over torsions */
/*     tett    torsion-torsion derivatives over torsions */
/*     tev     van der Waals derivatives over torsions */
/*     tec     charge-charge derivatives over torsions */
/*     tecd    charge-dipole derivatives over torsions */
/*     ted     dipole-dipole derivatives over torsions */
/*     tem     atomic multipole derivatives over torsions */
/*     tep     polarization derivatives over torsions */
/*     ter     reaction field derivatives over torsions */
/*     tes     solvation derivatives over torsions */
/*     telf    metal ligand field derivatives over torsions */
/*     teg     geometric restraint derivatives over torsions */
/*     tex     extra energy term derivatives over torsions */




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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  omega.i  --  dihedrals for torsional space computations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     dihed    current value in radians of each dihedral angle */
/*     nomega   number of dihedral angles allowed to rotate */
/*     iomega   numbers of two atoms defining rotation axis */
/*     zline    line number in Z-matrix of each dihedral angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  zcoord.i  --  Z-matrix internal coordinate definitions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     zbond   bond length used to define each Z-matrix atom */
/*     zang    bond angle used to define each Z-matrix atom */
/*     ztors   angle or torsion used to define Z-matrix atom */
/*     iz      defining atom numbers for each Z-matrix atom */




/*     set up the molecular mechanics calculation */

    initial_();
    getint_();
    mechanic_();
    initrot_();

/*     get the stepsize for numerical gradient calculation */

    delta = -1.;
    delta0 = .001;
    query = TRUE_;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___6);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&delta, (ftnlen)sizeof(doublereal)
		);
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
	io___7.ciunit = iounit_1.iout;
	s_wsfe(&io___7);
	do_fio(&c__1, (char *)&delta0, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___8.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___8);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = do_fio(&c__1, (char *)&delta, (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L10;
	}
    }
    if (delta <= 0.) {
	delta = delta0;
    }
    eps = -delta / 57.29577951308232088;

/*     make the call to get analytical torsional derivatives */

    gradrot_(&etot, derivs);

/*     print the total potential energy of the system */

    if (inform_1.digits >= 8) {
	io___12.ciunit = iounit_1.iout;
	s_wsfe(&io___12);
	do_fio(&c__1, (char *)&etot, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else if (inform_1.digits >= 6) {
	io___13.ciunit = iounit_1.iout;
	s_wsfe(&io___13);
	do_fio(&c__1, (char *)&etot, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else {
	io___14.ciunit = iounit_1.iout;
	s_wsfe(&io___14);
	do_fio(&c__1, (char *)&etot, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     print the energy breakdown over individual components */

    io___15.ciunit = iounit_1.iout;
    s_wsfe(&io___15);
    e_wsfe();
    if (inform_1.digits >= 8) {
	io___16.ciunit = iounit_1.iout;
	s_wsfe(&io___16);
	e_wsfe();
	io___17.ciunit = iounit_1.iout;
	s_wsfe(&io___17);
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
	io___18.ciunit = iounit_1.iout;
	s_wsfe(&io___18);
	e_wsfe();
	io___19.ciunit = iounit_1.iout;
	s_wsfe(&io___19);
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
	io___20.ciunit = iounit_1.iout;
	s_wsfe(&io___20);
	e_wsfe();
	io___21.ciunit = iounit_1.iout;
	s_wsfe(&io___21);
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

/*     print a header for the gradients of individual potentials */

    io___22.ciunit = iounit_1.iout;
    s_wsfe(&io___22);
    e_wsfe();
    if (inform_1.digits >= 8) {
	io___23.ciunit = iounit_1.iout;
	s_wsfe(&io___23);
	e_wsfe();
    } else if (inform_1.digits >= 6) {
	io___24.ciunit = iounit_1.iout;
	s_wsfe(&io___24);
	e_wsfe();
    } else {
	io___25.ciunit = iounit_1.iout;
	s_wsfe(&io___25);
	e_wsfe();
    }

/*     get numerical derivatives for each of the rotatable torsions */

    i__1 = omega_1.nomega;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] += delta / 2.;
	makexyz_();
	e0 = energy_();
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
	zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] -= delta;
	makexyz_();
	e = energy_();
	zcoord_1.ztors[omega_1.zline[i__ - 1] - 1] += delta / 2.;
	nderiv[i__ - 1] = (e - e0) / eps;

/*     print analytical gradients of each energy term for each atom */

	if (inform_1.digits >= 8) {
	    io___54.ciunit = iounit_1.iout;
	    s_wsfe(&io___54);
	    do_fio(&c__1, (char *)&iomega_ref(2, i__), (ftnlen)sizeof(integer)
		    );
	    do_fio(&c__1, (char *)&domega_1.teb[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tea[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teba[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teub[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iomega_ref(1, i__), (ftnlen)sizeof(integer)
		    );
	    do_fio(&c__1, (char *)&domega_1.teaa[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teopb[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teopd[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teid[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teit[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tet[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tept[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tebt[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tett[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tev[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tec[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tecd[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.ted[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tem[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tep[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.ter[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tes[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.telf[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teg[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tex[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	} else if (inform_1.digits >= 6) {
	    io___55.ciunit = iounit_1.iout;
	    s_wsfe(&io___55);
	    do_fio(&c__1, (char *)&iomega_ref(2, i__), (ftnlen)sizeof(integer)
		    );
	    do_fio(&c__1, (char *)&domega_1.teb[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tea[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teba[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teub[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teaa[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iomega_ref(1, i__), (ftnlen)sizeof(integer)
		    );
	    do_fio(&c__1, (char *)&domega_1.teopb[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teopd[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teid[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teit[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tet[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tept[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tebt[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tett[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tev[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tec[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tecd[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.ted[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tem[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tep[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.ter[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tes[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.telf[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teg[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tex[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	} else {
	    io___56.ciunit = iounit_1.iout;
	    s_wsfe(&io___56);
	    do_fio(&c__1, (char *)&iomega_ref(2, i__), (ftnlen)sizeof(integer)
		    );
	    do_fio(&c__1, (char *)&domega_1.teb[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tea[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teba[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teub[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teaa[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teopb[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&iomega_ref(1, i__), (ftnlen)sizeof(integer)
		    );
	    do_fio(&c__1, (char *)&domega_1.teopd[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teid[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teit[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tet[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tept[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tebt[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tett[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tev[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tec[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tecd[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.ted[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tem[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tep[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.ter[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tes[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.telf[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.teg[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&domega_1.tex[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}

/*     print numerical gradients of each energy term for each atom */

	if (inform_1.digits >= 8) {
	    io___57.ciunit = iounit_1.iout;
	    s_wsfe(&io___57);
	    do_fio(&c__1, (char *)&iomega_ref(2, i__), (ftnlen)sizeof(integer)
		    );
	    d__1 = (energi_1.eb - eb0) / eps;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    d__2 = (energi_1.ea - ea0) / eps;
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    d__3 = (energi_1.eba - eba0) / eps;
	    do_fio(&c__1, (char *)&d__3, (ftnlen)sizeof(doublereal));
	    d__4 = (energi_1.eub - eub0) / eps;
	    do_fio(&c__1, (char *)&d__4, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&iomega_ref(1, i__), (ftnlen)sizeof(integer)
		    );
	    d__5 = (energi_1.eaa - eaa0) / eps;
	    do_fio(&c__1, (char *)&d__5, (ftnlen)sizeof(doublereal));
	    d__6 = (energi_1.eopb - eopb0) / eps;
	    do_fio(&c__1, (char *)&d__6, (ftnlen)sizeof(doublereal));
	    d__7 = (energi_1.eopd - eopd0) / eps;
	    do_fio(&c__1, (char *)&d__7, (ftnlen)sizeof(doublereal));
	    d__8 = (energi_1.eid - eid0) / eps;
	    do_fio(&c__1, (char *)&d__8, (ftnlen)sizeof(doublereal));
	    d__9 = (energi_1.eit - eit0) / eps;
	    do_fio(&c__1, (char *)&d__9, (ftnlen)sizeof(doublereal));
	    d__10 = (energi_1.et - et0) / eps;
	    do_fio(&c__1, (char *)&d__10, (ftnlen)sizeof(doublereal));
	    d__11 = (energi_1.ept - ept0) / eps;
	    do_fio(&c__1, (char *)&d__11, (ftnlen)sizeof(doublereal));
	    d__12 = (energi_1.ebt - ebt0) / eps;
	    do_fio(&c__1, (char *)&d__12, (ftnlen)sizeof(doublereal));
	    d__13 = (energi_1.ett - ett0) / eps;
	    do_fio(&c__1, (char *)&d__13, (ftnlen)sizeof(doublereal));
	    d__14 = (energi_1.ev - ev0) / eps;
	    do_fio(&c__1, (char *)&d__14, (ftnlen)sizeof(doublereal));
	    d__15 = (energi_1.ec - ec0) / eps;
	    do_fio(&c__1, (char *)&d__15, (ftnlen)sizeof(doublereal));
	    d__16 = (energi_1.ecd - ecd0) / eps;
	    do_fio(&c__1, (char *)&d__16, (ftnlen)sizeof(doublereal));
	    d__17 = (energi_1.ed - ed0) / eps;
	    do_fio(&c__1, (char *)&d__17, (ftnlen)sizeof(doublereal));
	    d__18 = (energi_1.em - em0) / eps;
	    do_fio(&c__1, (char *)&d__18, (ftnlen)sizeof(doublereal));
	    d__19 = (energi_1.ep - ep0) / eps;
	    do_fio(&c__1, (char *)&d__19, (ftnlen)sizeof(doublereal));
	    d__20 = (energi_1.er - er0) / eps;
	    do_fio(&c__1, (char *)&d__20, (ftnlen)sizeof(doublereal));
	    d__21 = (energi_1.es - es0) / eps;
	    do_fio(&c__1, (char *)&d__21, (ftnlen)sizeof(doublereal));
	    d__22 = (energi_1.elf - elf0) / eps;
	    do_fio(&c__1, (char *)&d__22, (ftnlen)sizeof(doublereal));
	    d__23 = (energi_1.eg - eg0) / eps;
	    do_fio(&c__1, (char *)&d__23, (ftnlen)sizeof(doublereal));
	    d__24 = (energi_1.ex - ex0) / eps;
	    do_fio(&c__1, (char *)&d__24, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (inform_1.digits >= 6) {
	    io___58.ciunit = iounit_1.iout;
	    s_wsfe(&io___58);
	    do_fio(&c__1, (char *)&iomega_ref(2, i__), (ftnlen)sizeof(integer)
		    );
	    d__1 = (energi_1.eb - eb0) / eps;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    d__2 = (energi_1.ea - ea0) / eps;
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    d__3 = (energi_1.eba - eba0) / eps;
	    do_fio(&c__1, (char *)&d__3, (ftnlen)sizeof(doublereal));
	    d__4 = (energi_1.eub - eub0) / eps;
	    do_fio(&c__1, (char *)&d__4, (ftnlen)sizeof(doublereal));
	    d__5 = (energi_1.eaa - eaa0) / eps;
	    do_fio(&c__1, (char *)&d__5, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&iomega_ref(1, i__), (ftnlen)sizeof(integer)
		    );
	    d__6 = (energi_1.eopb - eopb0) / eps;
	    do_fio(&c__1, (char *)&d__6, (ftnlen)sizeof(doublereal));
	    d__7 = (energi_1.eopd - eopd0) / eps;
	    do_fio(&c__1, (char *)&d__7, (ftnlen)sizeof(doublereal));
	    d__8 = (energi_1.eid - eid0) / eps;
	    do_fio(&c__1, (char *)&d__8, (ftnlen)sizeof(doublereal));
	    d__9 = (energi_1.eit - eit0) / eps;
	    do_fio(&c__1, (char *)&d__9, (ftnlen)sizeof(doublereal));
	    d__10 = (energi_1.et - et0) / eps;
	    do_fio(&c__1, (char *)&d__10, (ftnlen)sizeof(doublereal));
	    d__11 = (energi_1.ept - ept0) / eps;
	    do_fio(&c__1, (char *)&d__11, (ftnlen)sizeof(doublereal));
	    d__12 = (energi_1.ebt - ebt0) / eps;
	    do_fio(&c__1, (char *)&d__12, (ftnlen)sizeof(doublereal));
	    d__13 = (energi_1.ett - ett0) / eps;
	    do_fio(&c__1, (char *)&d__13, (ftnlen)sizeof(doublereal));
	    d__14 = (energi_1.ev - ev0) / eps;
	    do_fio(&c__1, (char *)&d__14, (ftnlen)sizeof(doublereal));
	    d__15 = (energi_1.ec - ec0) / eps;
	    do_fio(&c__1, (char *)&d__15, (ftnlen)sizeof(doublereal));
	    d__16 = (energi_1.ecd - ecd0) / eps;
	    do_fio(&c__1, (char *)&d__16, (ftnlen)sizeof(doublereal));
	    d__17 = (energi_1.ed - ed0) / eps;
	    do_fio(&c__1, (char *)&d__17, (ftnlen)sizeof(doublereal));
	    d__18 = (energi_1.em - em0) / eps;
	    do_fio(&c__1, (char *)&d__18, (ftnlen)sizeof(doublereal));
	    d__19 = (energi_1.ep - ep0) / eps;
	    do_fio(&c__1, (char *)&d__19, (ftnlen)sizeof(doublereal));
	    d__20 = (energi_1.er - er0) / eps;
	    do_fio(&c__1, (char *)&d__20, (ftnlen)sizeof(doublereal));
	    d__21 = (energi_1.es - es0) / eps;
	    do_fio(&c__1, (char *)&d__21, (ftnlen)sizeof(doublereal));
	    d__22 = (energi_1.elf - elf0) / eps;
	    do_fio(&c__1, (char *)&d__22, (ftnlen)sizeof(doublereal));
	    d__23 = (energi_1.eg - eg0) / eps;
	    do_fio(&c__1, (char *)&d__23, (ftnlen)sizeof(doublereal));
	    d__24 = (energi_1.ex - ex0) / eps;
	    do_fio(&c__1, (char *)&d__24, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___59.ciunit = iounit_1.iout;
	    s_wsfe(&io___59);
	    do_fio(&c__1, (char *)&iomega_ref(2, i__), (ftnlen)sizeof(integer)
		    );
	    d__1 = (energi_1.eb - eb0) / eps;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    d__2 = (energi_1.ea - ea0) / eps;
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    d__3 = (energi_1.eba - eba0) / eps;
	    do_fio(&c__1, (char *)&d__3, (ftnlen)sizeof(doublereal));
	    d__4 = (energi_1.eub - eub0) / eps;
	    do_fio(&c__1, (char *)&d__4, (ftnlen)sizeof(doublereal));
	    d__5 = (energi_1.eaa - eaa0) / eps;
	    do_fio(&c__1, (char *)&d__5, (ftnlen)sizeof(doublereal));
	    d__6 = (energi_1.eopb - eopb0) / eps;
	    do_fio(&c__1, (char *)&d__6, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&iomega_ref(1, i__), (ftnlen)sizeof(integer)
		    );
	    d__7 = (energi_1.eopd - eopd0) / eps;
	    do_fio(&c__1, (char *)&d__7, (ftnlen)sizeof(doublereal));
	    d__8 = (energi_1.eid - eid0) / eps;
	    do_fio(&c__1, (char *)&d__8, (ftnlen)sizeof(doublereal));
	    d__9 = (energi_1.eit - eit0) / eps;
	    do_fio(&c__1, (char *)&d__9, (ftnlen)sizeof(doublereal));
	    d__10 = (energi_1.et - et0) / eps;
	    do_fio(&c__1, (char *)&d__10, (ftnlen)sizeof(doublereal));
	    d__11 = (energi_1.ept - ept0) / eps;
	    do_fio(&c__1, (char *)&d__11, (ftnlen)sizeof(doublereal));
	    d__12 = (energi_1.ebt - ebt0) / eps;
	    do_fio(&c__1, (char *)&d__12, (ftnlen)sizeof(doublereal));
	    d__13 = (energi_1.ett - ett0) / eps;
	    do_fio(&c__1, (char *)&d__13, (ftnlen)sizeof(doublereal));
	    d__14 = (energi_1.ev - ev0) / eps;
	    do_fio(&c__1, (char *)&d__14, (ftnlen)sizeof(doublereal));
	    d__15 = (energi_1.ec - ec0) / eps;
	    do_fio(&c__1, (char *)&d__15, (ftnlen)sizeof(doublereal));
	    d__16 = (energi_1.ecd - ecd0) / eps;
	    do_fio(&c__1, (char *)&d__16, (ftnlen)sizeof(doublereal));
	    d__17 = (energi_1.ed - ed0) / eps;
	    do_fio(&c__1, (char *)&d__17, (ftnlen)sizeof(doublereal));
	    d__18 = (energi_1.em - em0) / eps;
	    do_fio(&c__1, (char *)&d__18, (ftnlen)sizeof(doublereal));
	    d__19 = (energi_1.ep - ep0) / eps;
	    do_fio(&c__1, (char *)&d__19, (ftnlen)sizeof(doublereal));
	    d__20 = (energi_1.er - er0) / eps;
	    do_fio(&c__1, (char *)&d__20, (ftnlen)sizeof(doublereal));
	    d__21 = (energi_1.es - es0) / eps;
	    do_fio(&c__1, (char *)&d__21, (ftnlen)sizeof(doublereal));
	    d__22 = (energi_1.elf - elf0) / eps;
	    do_fio(&c__1, (char *)&d__22, (ftnlen)sizeof(doublereal));
	    d__23 = (energi_1.eg - eg0) / eps;
	    do_fio(&c__1, (char *)&d__23, (ftnlen)sizeof(doublereal));
	    d__24 = (energi_1.ex - ex0) / eps;
	    do_fio(&c__1, (char *)&d__24, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     print a header for the analytical vs. numerical comparison */

    if (inform_1.digits >= 8) {
	io___60.ciunit = iounit_1.iout;
	s_wsfe(&io___60);
	e_wsfe();
    } else if (inform_1.digits >= 6) {
	io___61.ciunit = iounit_1.iout;
	s_wsfe(&io___61);
	e_wsfe();
    } else {
	io___62.ciunit = iounit_1.iout;
	s_wsfe(&io___62);
	e_wsfe();
    }

/*     print comparison of analytical and numerical derivatives */

    if (inform_1.digits >= 8) {
	i__1 = omega_1.nomega;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___63.ciunit = iounit_1.iout;
	    s_wsfe(&io___63);
	    do_fio(&c__1, (char *)&iomega_ref(2, i__), (ftnlen)sizeof(integer)
		    );
	    do_fio(&c__1, (char *)&iomega_ref(1, i__), (ftnlen)sizeof(integer)
		    );
	    do_fio(&c__1, (char *)&derivs[i__ - 1], (ftnlen)sizeof(doublereal)
		    );
	    do_fio(&c__1, (char *)&nderiv[i__ - 1], (ftnlen)sizeof(doublereal)
		    );
	    e_wsfe();
	}
    } else if (inform_1.digits >= 6) {
	i__1 = omega_1.nomega;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___64.ciunit = iounit_1.iout;
	    s_wsfe(&io___64);
	    do_fio(&c__1, (char *)&iomega_ref(2, i__), (ftnlen)sizeof(integer)
		    );
	    do_fio(&c__1, (char *)&iomega_ref(1, i__), (ftnlen)sizeof(integer)
		    );
	    do_fio(&c__1, (char *)&derivs[i__ - 1], (ftnlen)sizeof(doublereal)
		    );
	    do_fio(&c__1, (char *)&nderiv[i__ - 1], (ftnlen)sizeof(doublereal)
		    );
	    e_wsfe();
	}
    } else {
	i__1 = omega_1.nomega;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___65.ciunit = iounit_1.iout;
	    s_wsfe(&io___65);
	    do_fio(&c__1, (char *)&iomega_ref(2, i__), (ftnlen)sizeof(integer)
		    );
	    do_fio(&c__1, (char *)&iomega_ref(1, i__), (ftnlen)sizeof(integer)
		    );
	    do_fio(&c__1, (char *)&derivs[i__ - 1], (ftnlen)sizeof(doublereal)
		    );
	    do_fio(&c__1, (char *)&nderiv[i__ - 1], (ftnlen)sizeof(doublereal)
		    );
	    e_wsfe();
	}
    }

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */

#undef iomega_ref


/* Main program alias */ int testrot_ () { MAIN__ (); return 0; }
