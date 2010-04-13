/* analyze.f -- translated by f2c (version 20050501).
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
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

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
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

struct {
    doublereal xbox, ybox, zbox, alpha, beta, gamma, xbox2, ybox2, zbox2, 
	    box34, lvec[9]	/* was [3][3] */, recip[9]	/* was [3][3] 
	    */, volbox, beta_sin__, beta_cos__, gamma_sin__, gamma_cos__, 
	    beta_term__, gamma_term__;
    logical orthogonal, monoclinic, triclinic, octahedron;
    char spacegrp[10];
} boxes_;

#define boxes_1 boxes_

struct {
    doublereal vdwcut, chgcut, dplcut, mpolecut, vdwtaper, chgtaper, dpltaper,
	     mpoletaper, ewaldcut;
    logical use_ewald__, use_lights__, use_list__, use_vlist__, use_clist__, 
	    use_mlist__;
} cutoff_;

#define cutoff_1 cutoff_

struct {
    doublereal aewald;
    char boundary[7];
} ewald_;

#define ewald_1 ewald_

struct {
    integer biotyp[10000];
    char forcefield[20];
} fields_;

#define fields_1 fields_

struct {
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

struct {
    doublereal bsmod1[100], bsmod2[100], bsmod3[100], table[1200]	/* 
	    was [400][3] */, qgrid[2000000]	/* was [2][100][100][100] */, 
	    qfac[1000000]	/* was [100][100][100] */, thetai1[1000000]	
	    /* was [4][10][25000] */, thetai2[1000000]	/* was [4][10][25000] 
	    */, thetai3[1000000]	/* was [4][10][25000] */;
    integer nfft1, nfft2, nfft3, bsorder, iprime[45]	/* was [15][3] */, 
	    igrid[75000]	/* was [3][25000] */;
} pme_;

#define pme_1 pme_

struct {
    logical use_bond__, use_angle__, use_strbnd__, use_urey__, use_angang__, 
	    use_opbend__, use_opdist__, use_improp__, use_imptor__, 
	    use_tors__, use_pitors__, use_strtor__, use_tortor__, use_vdw__, 
	    use_charge__, use_chgdpl__, use_dipole__, use_mpole__, 
	    use_polar__, use_rxnfld__, use_solv__, use_metal__, use_geom__, 
	    use_extra__, use_born__, use_orbit__;
} potent_;

#define potent_1 potent_

struct {
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_

struct {
    doublereal kaa[100000];
    integer nangang, iaa[200000]	/* was [2][100000] */;
} angang_;

#define angang_1 angang_

struct {
    doublereal ak[75000], anat[75000], afld[75000];
    integer nangle, iang[300000]	/* was [4][75000] */;
} angle_;

#define angle_1 angle_

struct {
    doublereal angunit, stbnunit, aaunit, opbunit, opdunit, cang, qang, pang, 
	    sang, copb, qopb, popb, sopb, copd, qopd, popd, sopd;
    char angtyp[600000], opbtyp[8];
} angpot_;

#define angpot_1 angpot_

struct {
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

struct {
    integer nbitor, ibitor[500000]	/* was [5][100000] */;
} bitor_;

#define bitor_1 bitor_

struct {
    doublereal bk[50000], bl[50000];
    integer nbond, ibnd[100000]	/* was [2][50000] */;
} bond_;

#define bond_1 bond_

struct {
    doublereal pchg[25000];
    integer nion, iion[25000], jion[25000], kion[25000], chglist[25000];
} charge_;

#define charge_1 charge_

struct {
    doublereal bdpl[50000], sdpl[50000];
    integer ndipole, idpl[100000]	/* was [2][50000] */;
} dipole_;

#define dipole_1 dipole_

struct {
    doublereal kprop[100000], vprop[100000];
    integer niprop, iiprop[400000]	/* was [4][100000] */;
} improp_;

#define improp_1 improp_

struct {
    doublereal itors1[400000]	/* was [4][100000] */, itors2[400000]	/* 
	    was [4][100000] */, itors3[400000]	/* was [4][100000] */;
    integer nitors, iitors[400000]	/* was [4][100000] */;
} imptor_;

#define imptor_1 imptor_

struct {
    doublereal electron[1000], ionize[1000], repulse[1000], sslope[500], 
	    tslope[500], sslope5[200], tslope5[200], sslope4[200], tslope4[
	    200];
    char kpi[4000], kpi5[1600], kpi4[1600];
} korbs_;

#define korbs_1 korbs_

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
    doublereal rad[5000], eps[5000], rad4[5000], eps4[5000], reduct[5000];
} kvdws_;

#define kvdws_1 kvdws_

struct {
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

struct {
    doublereal opbk[75000];
    integer nopbend, iopb[75000];
} opbend_;

#define opbend_1 opbend_

struct {
    doublereal opdk[25000];
    integer nopdist, iopd[100000]	/* was [4][25000] */;
} opdist_;

#define opdist_1 opdist_

struct {
    integer norbit, iorbit[100], reorbit, piperp[300]	/* was [3][100] */, 
	    nbpi, ibpi[600]	/* was [3][200] */, ntpi, itpi[800]	/* 
	    was [2][400] */;
    logical listpi[25000];
} piorbs_;

#define piorbs_1 piorbs_

struct {
    doublereal bkpi[200], blpi[200], kslope[200], lslope[200], torsp2[400];
} pistuf_;

#define pistuf_1 pistuf_

struct {
    doublereal kpit[100000];
    integer npitors, ipit[600000]	/* was [6][100000] */;
} pitors_;

#define pitors_1 pitors_

struct {
    doublereal polarity[25000], thole[25000], pdamp[25000], uind[75000]	/* 
	    was [3][25000] */, uinp[75000]	/* was [3][25000] */, uinds[
	    75000]	/* was [3][25000] */, uinps[75000]	/* was [3][
	    25000] */;
    integer npolar;
} polar_;

#define polar_1 polar_

struct {
    integer np11[25000], ip11[2500000]	/* was [100][25000] */, np12[25000], 
	    ip12[1250000]	/* was [50][25000] */, np13[25000], ip13[
	    1250000]	/* was [50][25000] */, np14[25000], ip14[1250000]	
	    /* was [50][25000] */;
} polgrp_;

#define polgrp_1 polgrp_

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
    doublereal sbk[150000]	/* was [2][75000] */;
    integer nstrbnd, isb[225000]	/* was [3][75000] */;
} strbnd_;

#define strbnd_1 strbnd_

struct {
    doublereal kst[300000]	/* was [3][100000] */;
    integer nstrtor, ist[200000]	/* was [2][100000] */;
} strtor_;

#define strtor_1 strtor_

struct {
    doublereal tors1[400000]	/* was [4][100000] */, tors2[400000]	/* 
	    was [4][100000] */, tors3[400000]	/* was [4][100000] */, tors4[
	    400000]	/* was [4][100000] */, tors5[400000]	/* was [4][
	    100000] */, tors6[400000]	/* was [4][100000] */;
    integer ntors, itors[400000]	/* was [4][100000] */;
} tors_;

#define tors_1 tors_

struct {
    integer ntortor, itt[300000]	/* was [3][100000] */;
} tortor_;

#define tortor_1 tortor_

struct {
    doublereal uk[75000], ul[75000];
    integer nurey, iury[225000]	/* was [3][75000] */;
} urey_;

#define urey_1 urey_

struct {
    doublereal radmin[1000000]	/* was [1000][1000] */, epsilon[1000000]	
	    /* was [1000][1000] */, radmin4[1000000]	/* was [1000][1000] */
	    , epsilon4[1000000]	/* was [1000][1000] */, radhbnd[1000000]	
	    /* was [1000][1000] */, epshbnd[1000000]	/* was [1000][1000] */
	    , kred[25000];
    integer ired[25000], nvdw, ivdw[25000], jvdw[25000], nvt, ivt[25000], jvt[
	    25000];
} vdw_;

#define vdw_1 vdw_

struct {
    doublereal einter;
} inter_;

#define inter_1 inter_

struct {
    integer neb, nea, neba, neub, neaa, neopb, neopd, neid, neit, net, nept, 
	    nebt, nett, nev, nec, necd, ned, nem, nep, new__, ner, nes, nelf, 
	    neg, nex;
} action_;

#define action_1 action_

struct {
    doublereal esum, eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, 
	    ebt, ett, ev, ec, ecd, ed, em, ep, er, es, elf, eg, ex;
} energi_;

#define energi_1 energi_

struct {
    doublereal electric, dielec, ebuffer, c2scale, c3scale, c4scale, c5scale;
    logical neutnbr, neutcut;
} chgpot_;

#define chgpot_1 chgpot_

struct {
    doublereal netchg, netdpl, netqdp[3], xdpl, ydpl, zdpl, xxqdp, xyqdp, 
	    xzqdp, yxqdp, yyqdp, yzqdp, zxqdp, zyqdp, zzqdp;
} moment_;

#define moment_1 moment_

struct {
    doublereal vir[9]	/* was [3][3] */;
} virial_;

#define virial_1 virial_

struct {
    doublereal aesum[25000], aeb[25000], aea[25000], aeba[25000], aeub[25000],
	     aeaa[25000], aeopb[25000], aeopd[25000], aeid[25000], aeit[25000]
	    , aet[25000], aept[25000], aebt[25000], aett[25000], aev[25000], 
	    aec[25000], aecd[25000], aed[25000], aem[25000], aep[25000], aer[
	    25000], aes[25000], aelf[25000], aeg[25000], aex[25000];
} analyz_;

#define analyz_1 analyz_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  program analyze  --  energy partitioning and analysis  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "analyze" computes and displays the total potential energy; */
/*     options are provided to display system and force field info, */
/*     partition the energy by atom or by potential function type, */
/*     show force field parameters by atom; output the large energy */
/*     interactions and find electrostatic and inertial properties */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 The TINKER Analysis Facility can Provide"
	    " :\002,//,\002 General System and Force Field Information [G]"
	    "\002,/,\002 Force Field Parameters for Interactions [P]\002,/"
	    ",\002 Total Potential Energy and its Components [E]\002,/,\002 E"
	    "nergy Breakdown over Each of the Atoms [A]\002,/,\002 List of th"
	    "e Large Individual Interactions [L]\002,/,\002 Details for All I"
	    "ndividual Interactions [D]\002,/,\002 Electrostatic, Inertial & "
	    "Virial Properties [M]\002)";
    static char fmt_30[] = "(/,\002 Enter the Desired Analysis Types\002,"
	    "\002 [G,P,E,A,L,D,M] :  \002,$)";
    static char fmt_40[] = "(a120)";
    static char fmt_60[] = "(/,\002 List Atoms for which Output is Desire"
	    "d\002,\002 [ALL] :  \002/,\002    >  \002,$)";
    static char fmt_70[] = "(a120)";
    static char fmt_90[] = "(/,\002 Analysis for Archive Structure :\002,8x,"
	    "i8)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *), f_rew(alist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int systyze_(void), mechanic_(void);
    static logical dodetail, doenergy;
    extern integer freeunit_(void);
    extern /* Subroutine */ int paramyze_(logical *);
    static logical dosystem;
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j, list[20], ixyz;
    extern /* Subroutine */ int final_(void);
    static integer frame;
    static logical exist, active[25000];
    static char record[120];
    static logical doatom;
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char letter[1], string[120];
    extern /* Subroutine */ int suffix_(char *, char *, ftnlen, ftnlen), 
	    getxyz_(void), analyz4_(void), analyz6_(void), analyz8_(void);
    static logical dolarge, doparam;
    extern /* Subroutine */ int initial_(void), nextarg_(char *, logical *, 
	    ftnlen), enrgyze_(void), version_(char *, char *, ftnlen, ftnlen);
    static logical doprops;
    extern /* Subroutine */ int readxyz_(integer *), atomyze_(logical *);
    static char xyzfile[120];
    extern /* Subroutine */ int propyze_(void);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___5 = { 1, 0, 0, fmt_40, 0 };
    static icilist io___16 = { 1, string, 1, 0, 120, 1 };
    static cilist io___17 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_70, 0 };
    static icilist io___20 = { 1, record, 1, 0, 120, 1 };
    static cilist io___26 = { 0, 0, 0, fmt_90, 0 };




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




/*     set up the structure and mechanics calculation */

    initial_();
    getxyz_();
    mechanic_();

/*     get the desired types of analysis to be performed */

    nextarg_(string, &exist, (ftnlen)120);
    if (! exist) {
	io___3.ciunit = iounit_1.iout;
	s_wsfe(&io___3);
	e_wsfe();
L20:
	io___4.ciunit = iounit_1.iout;
	s_wsfe(&io___4);
	e_wsfe();
	io___5.ciunit = iounit_1.input;
	i__1 = s_rsfe(&io___5);
	if (i__1 != 0) {
	    goto L20;
	}
	i__1 = do_fio(&c__1, string, (ftnlen)120);
	if (i__1 != 0) {
	    goto L20;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L20;
	}
    }

/*     set option control flags based desired analysis types */

    dosystem = FALSE_;
    doparam = FALSE_;
    doenergy = FALSE_;
    doatom = FALSE_;
    dolarge = FALSE_;
    dodetail = FALSE_;
    doprops = FALSE_;
    upcase_(string, (ftnlen)120);
    i__1 = trimtext_(string, (ftnlen)120);
    for (i__ = 1; i__ <= i__1; ++i__) {
	*(unsigned char *)letter = *(unsigned char *)&string[i__ - 1];
	if (*(unsigned char *)letter == 'G') {
	    dosystem = TRUE_;
	}
	if (*(unsigned char *)letter == 'P') {
	    doparam = TRUE_;
	}
	if (*(unsigned char *)letter == 'E') {
	    doenergy = TRUE_;
	}
	if (*(unsigned char *)letter == 'A') {
	    doatom = TRUE_;
	}
	if (*(unsigned char *)letter == 'L') {
	    dolarge = TRUE_;
	}
	if (*(unsigned char *)letter == 'D') {
	    dodetail = TRUE_;
	}
	if (*(unsigned char *)letter == 'M') {
	    doprops = TRUE_;
	}
    }

/*     get the list of atoms for which output is desired */

    if (doatom || doparam) {
	for (i__ = 1; i__ <= 20; ++i__) {
	    list[i__ - 1] = 0;
	}
	if (exist) {
	    for (i__ = 1; i__ <= 20; ++i__) {
		nextarg_(string, &exist, (ftnlen)120);
		if (! exist) {
		    goto L50;
		}
		i__1 = s_rsli(&io___16);
		if (i__1 != 0) {
		    goto L50;
		}
		i__1 = do_lio(&c__3, &c__1, (char *)&list[i__ - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L50;
		}
		i__1 = e_rsli();
		if (i__1 != 0) {
		    goto L50;
		}
	    }
L50:
	    ;
	} else {
	    io___17.ciunit = iounit_1.iout;
	    s_wsfe(&io___17);
	    e_wsfe();
	    io___18.ciunit = iounit_1.input;
	    s_rsfe(&io___18);
	    do_fio(&c__1, record, (ftnlen)120);
	    e_rsfe();
	    i__1 = s_rsli(&io___20);
	    if (i__1 != 0) {
		goto L80;
	    }
	    for (i__ = 1; i__ <= 20; ++i__) {
		i__1 = do_lio(&c__3, &c__1, (char *)&list[i__ - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		    goto L80;
		}
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L80;
	    }
L80:
	    ;
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    active[i__ - 1] = TRUE_;
	}
	i__ = 1;
	while(list[i__ - 1] != 0) {
	    if (i__ == 1) {
		i__1 = atoms_1.n;
		for (j = 1; j <= i__1; ++j) {
		    active[j - 1] = FALSE_;
		}
	    }
	    if (list[i__ - 1] > 0) {
		active[list[i__ - 1] - 1] = TRUE_;
		++i__;
	    } else {
		i__3 = (i__2 = list[i__], abs(i__2));
		for (j = (i__1 = list[i__ - 1], abs(i__1)); j <= i__3; ++j) {
		    active[j - 1] = TRUE_;
		}
		i__ += 2;
	    }
	}
    }

/*     setup to write out the large individual energy terms */

    if (dolarge) {
	inform_1.verbose = TRUE_;
    } else {
	inform_1.verbose = FALSE_;
    }

/*     setup to write out all of the individual energy terms */

    if (dodetail) {
	doenergy = TRUE_;
	inform_1.debug = TRUE_;
    } else {
	inform_1.debug = FALSE_;
    }

/*     reopen the coordinates file and read the first structure */

    frame = 0;
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

/*     get info on the molecular system and force field */

    if (dosystem) {
	systyze_();
    }

/*     get parameters used for molecular mechanics potentials */

    if (doparam) {
	paramyze_(active);
    }

/*     perform analysis for each successive coordinate structure */

    while(! inform_1.abort) {
	++frame;
	if (frame > 1) {
	    io___26.ciunit = iounit_1.iout;
	    s_wsfe(&io___26);
	    do_fio(&c__1, (char *)&frame, (ftnlen)sizeof(integer));
	    e_wsfe();
	}

/*     make the call to compute the potential energy */

	if (doenergy || doatom || dolarge) {
	    enrgyze_();
	}

/*     energy partitioning by potential energy components */

	if (doenergy) {
	    if (inform_1.digits >= 8) {
		analyz8_();
	    } else if (inform_1.digits >= 6) {
		analyz6_();
	    } else {
		analyz4_();
	    }
	}

/*     get various electrostatic and inertial properties */

	if (doprops) {
	    inform_1.debug = FALSE_;
	    propyze_();
	    if (dodetail) {
		inform_1.debug = TRUE_;
	    }
	}

/*     energy partitioning over the individual atoms */

	if (doatom) {
	    atomyze_(active);
	}

/*     attempt to read next structure from the coordinate file */

	readxyz_(&ixyz);
    }

/*     perform any final tasks before program exit */

    cl__1.cerr = 0;
    cl__1.cunit = ixyz;
    cl__1.csta = 0;
    f_clos(&cl__1);
    if (dodetail) {
	inform_1.debug = FALSE_;
    }
    final_();
    return 0;
} /* MAIN__ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine systyze  --  system & force field information  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "systyze" is an auxiliary routine for the analyze program */
/*     that prints general information about the molecular system */
/*     the the force field model */


/* Subroutine */ int systyze_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Overall System Contents :\002,//,\002 Nu"
	    "mber of Atoms\002,25x,i8,/,\002 Number of Molecules\002,21x,i8,/,"
	    "\002 Total System Mass\002,19x,f12.4)";
    static char fmt_20[] = "(\002 System Density\002,22x,f12.4)";
    static char fmt_30[] = "(/,\002 Periodic Boundary Conditions :\002,//"
	    ",\002 a-Axis Length\002,23x,f12.4,/,\002 b-Axis Length\002,23x,f"
	    "12.4,/,\002 c-Axis Length\002,23x,f12.4,/,\002 Alpha Angle\002,2"
	    "5x,f12.4,/,\002 Beta Angle\002,26x,f12.4,/,\002 Gamma Angle\002,"
	    "25x,f12.4,/,\002 Cell Volume\002,25x,f12.4,/,\002 Lattice Typ"
	    "e\002,16x,a20)";
    static char fmt_40[] = "(\002 Space Group\002,17x,a20)";
    static char fmt_50[] = "(/,\002 Force Field Name :\002,10x,a20)";
    static char fmt_60[] = "(/,\002 VDW Function\002,16x,a20,/,\002 Size Des"
	    "criptor\002,13x,a20,/,\002 Size Unit Type\002,14x,a20,/,\002 Siz"
	    "e Combining Rule\002,9x,a20,/,\002 Well Depth Rule\002,13x,a20)";
    static char fmt_70[] = "(\002 VDW Cutoff\002,26x,f12.4)";
    static char fmt_80[] = "()";
    static char fmt_90[] = "(\002 Electrostatics\002,14x,a20)";
    static char fmt_100[] = "(\002 Electrostatics\002,14x,a20)";
    static char fmt_110[] = "(\002 Electrostatics\002,14x,a20)";
    static char fmt_120[] = "(\002 Electrostatics\002,14x,a20)";
    static char fmt_130[] = "(/,\002 Particle Mesh Ewald :\002,//,\002 Ewald"
	    " Coefficient\002,19x,f12.4,/,\002 Real-Space Cutoff\002,19x,f12."
	    "4,/,\002 Grid Dimensions\002,21x,3i4,/,\002 B-Spline Order\002,2"
	    "6x,i8,/,\002 Boundary Condition\002,10x,a20)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal dens;
    static char label[20*5], value[20];
    extern /* Subroutine */ int justify_(char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___27 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_130, 0 };



#define label_ref(a_0,a_1) &label[(a_1)*20 + a_0 - 20]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  bound.i  --  control of periodic boundary conditions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     polycut       cutoff distance for infinite polymer nonbonds */
/*     polycut2      square of infinite polymer nonbond cutoff */
/*     use_bounds    flag to use periodic boundary conditions */
/*     use_replica   flag to use replicates for periodic system */
/*     use_polymer   flag to mark presence of infinite polymer */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  boxes.i  --  parameters for periodic boundary conditions  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     xbox        length of a-axis of periodic box in Angstroms */
/*     ybox        length of b-axis of periodic box in Angstroms */
/*     zbox        length of c-axis of periodic box in Angstroms */
/*     alpha       angle between b- and c-axes of box in degrees */
/*     beta        angle between a- and c-axes of box in degrees */
/*     gamma       angle between a- and b-axes of box in degrees */
/*     xbox2       half of the a-axis length of periodic box */
/*     ybox2       half of the b-axis length of periodic box */
/*     zbox2       half of the c-axis length of periodic box */
/*     box34       three-fourths axis length of truncated octahedron */
/*     lvec        real space lattice vectors as matrix rows */
/*     recip       reciprocal lattice vectors as matrix columns */
/*     volbox      volume in Ang**3 of the periodic box */
/*     beta_sin    sine of the beta periodic box angle */
/*     beta_cos    cosine of the beta periodic box angle */
/*     gamma_sin   sine of the gamma periodic box angle */
/*     gamma_cos   cosine of the gamma periodic box angle */
/*     beta_term   term used in generating triclinic box */
/*     gamma_term  term used in generating triclinic box */
/*     orthogonal  flag to mark periodic box as orthogonal */
/*     monoclinic  flag to mark periodic box as monoclinic */
/*     triclinic   flag to mark periodic box as triclinic */
/*     octahedron  flag to mark box as truncated octahedron */
/*     spacegrp    space group symbol for the unitcell type */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cutoff.i  --  cutoff distances for energy interactions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     vdwcut      cutoff distance for van der Waals interactions */
/*     chgcut      cutoff distance for charge-charge interactions */
/*     dplcut      cutoff distance for dipole-dipole interactions */
/*     mpolecut    cutoff distance for atomic multipole interactions */
/*     vdwtaper    distance at which van der Waals switching begins */
/*     chgtaper    distance at which charge-charge switching begins */
/*     dpltaper    distance at which dipole-dipole switching begins */
/*     mpoletaper  distance at which atomic multipole switching begins */
/*     ewaldcut    cutoff distance for direct space Ewald summation */
/*     use_ewald   logical flag governing use of Ewald summation */
/*     use_lights  logical flag governing use of method of lights */
/*     use_list    logical flag governing use of any neighbor lists */
/*     use_vlist   logical flag governing use of vdw neighbor lists */
/*     use_clist   logical flag governing use of charge neighbor lists */
/*     use_mlist   logical flag governing use of multipole neighbor lists */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  ewald.i  --  parameters and options for Ewald summation  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     aewald     Ewald convergence coefficient value (Ang-1) */
/*     boundary   Ewald boundary condition; none, tinfoil or vacuum */




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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  molcul.i  --  individual molecules within current system  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     molmass   molecular weight for each molecule in the system */
/*     totmass   total weight of all the molecules in the system */
/*     nmol      total number of separate molecules in the system */
/*     kmol      contiguous list of the atoms in each molecule */
/*     imol      first and last atom of each molecule in the list */
/*     molcule   number of the molecule to which each atom belongs */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  pme.i  --  values for particle mesh Ewald summation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxorder   maximum order of the B-spline approximation */
/*     maxprime   maximum number of prime factors of FFT dimension */
/*     maxtable   maximum size of the FFT table array */

/*     bsmod1     B-spline moduli along the a-axis direction */
/*     bsmod2     B-spline moduli along the b-axis direction */
/*     bsmod3     B-spline moduli along the c-axis direction */
/*     table      intermediate array used by the FFT calculation */
/*     qgrid      values on the particle mesh Ewald charge grid */
/*     qfac       prefactors for particle mesh Ewald charge grid */
/*     thetai1    B-spline coefficients along the a-axis */
/*     thetai2    B-spline coefficients along the b-axis */
/*     thetai3    B-spline coefficients along the c-axis */
/*     nfft1      number of grid points along the a-axis direction */
/*     nfft2      number of grid points along the b-axis direction */
/*     nfft3      number of grid points along the c-axis direction */
/*     bsorder    order of the PME B-spline approximation */
/*     iprime     prime factorization of each FFT dimension */
/*     igrid      initial Ewald charge grid values for B-spline */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  potent.i  --  usage of each potential energy component  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     use_bond    logical flag governing use of bond stretch potential */
/*     use_angle   logical flag governing use of angle bend potential */
/*     use_strbnd  logical flag governing use of stretch-bend potential */
/*     use_urey    logical flag governing use of Urey-Bradley potential */
/*     use_angang  logical flag governing use of angle-angle cross term */
/*     use_opbend  logical flag governing use of out-of-plane bend term */
/*     use_opdist  logical flag governing use of out-of-plane distance */
/*     use_improp  logical flag governing use of improper dihedral term */
/*     use_imptor  logical flag governing use of improper torsion term */
/*     use_tors    logical flag governing use of torsional potential */
/*     use_pitors  logical flag governing use of pi-orbital torsion term */
/*     use_strtor  logical flag governing use of stretch-torsion term */
/*     use_tortor  logical flag governing use of torsion-torsion term */
/*     use_vdw     logical flag governing use of vdw der Waals potential */
/*     use_charge  logical flag governing use of charge-charge potential */
/*     use_chgdpl  logical flag governing use of charge-dipole potential */
/*     use_dipole  logical flag governing use of dipole-dipole potential */
/*     use_mpole   logical flag governing use of multipole potential */
/*     use_polar   logical flag governing use of polarization term */
/*     use_rxnfld  logical flag governing use of reaction field term */
/*     use_solv    logical flag governing use of continuum solvation */
/*     use_metal   logical flag governing use of ligand field term */
/*     use_geom    logical flag governing use of geometric restraints */
/*     use_extra   logical flag governing use of extra potential term */
/*     use_born    logical flag governing use of Born radii values */
/*     use_orbit   logical flag governing use of pisystem computation */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  units.i  --  physical constants and unit conversions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     literature reference: */

/*     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended */
/*     Values of the Fundamental Physical Constants: 2006", Reviews of */
/*     Modern Physics, 80, 633-730 (2008) */

/*     The "2006 CODATA Recommended Values" are also available from */
/*     the NIST Reference on Constants, Units, and Uncertainty site */
/*     at http://physics.nist.gov/cuu/index.html */

/*     Most values below are derived from the 2006 CODATA reference */

/*     The conversion from calorie to Joule is the definition of the */
/*     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992) */

/*     The "coulomb" energy conversion factor is found by dimensional */
/*     analysis of Coulomb's Law, ie, by dividing the square of the */
/*     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is */
/*     the permittivity of vacuum (the "electric constant"); note that */
/*     eps0 is typically given in F/m, equivalent to C**2/(J-m) */

/*     The approximate value used for the Debye, 3.33564 x 10-30 C-m, */
/*     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997) */

/*     The value of "prescon" is based on definition of 1 atmosphere */
/*     as 101325 Pa set by the 10th Conference Generale des Poids et */
/*     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3 */

/*     avogadro    Avogadro's number (N) in particles/mole */
/*     lightspd    speed of light in vacuum (c) in cm/ps */
/*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K */
/*     gasconst    ideal gas constant (R) in kcal/mole/K */
/*     emass       mass of an electron in atomic mass units */
/*     joule       conversion from calories to joules */
/*     convert     conversion from kcal to g*Ang**2/ps**2 */
/*     bohr        conversion from Bohrs to Angstroms */
/*     hartree     conversion from Hartree to kcal/mole */
/*     evolt       conversion from Hartree to electron-volts */
/*     efreq       conversion from Hartree to cm-1 */
/*     coulomb     conversion from electron**2/Ang to kcal/mole */
/*     debye       conversion from electron-Ang to Debyes */
/*     prescon     conversion from kcal/mole/Ang**3 to Atm */




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




/*     info on number of atoms, molecules and mass */

    if (atoms_1.n != 0) {
	io___27.ciunit = iounit_1.iout;
	s_wsfe(&io___27);
	do_fio(&c__1, (char *)&atoms_1.n, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&molcul_1.nmol, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&molcul_1.totmass, (ftnlen)sizeof(doublereal));
	e_wsfe();
	if (bound_1.use_bounds__) {
	    dens = 1e24 / boxes_1.volbox * (molcul_1.totmass / 6.02214179e23);
	    io___29.ciunit = iounit_1.iout;
	    s_wsfe(&io___29);
	    do_fio(&c__1, (char *)&dens, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     periodic box dimensions and crystal lattice type */

    if (bound_1.use_bounds__) {
	s_copy(value, "ORTHOGONAL", (ftnlen)20, (ftnlen)10);
	if (boxes_1.monoclinic) {
	    s_copy(value, "MONOCLINIC", (ftnlen)20, (ftnlen)10);
	}
	if (boxes_1.triclinic) {
	    s_copy(value, "TRICLINIC", (ftnlen)20, (ftnlen)9);
	}
	if (boxes_1.octahedron) {
	    s_copy(value, "TRUNCATED OCTAHEDRON", (ftnlen)20, (ftnlen)20);
	}
	justify_(value, (ftnlen)20);
	io___31.ciunit = iounit_1.iout;
	s_wsfe(&io___31);
	do_fio(&c__1, (char *)&boxes_1.xbox, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&boxes_1.ybox, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&boxes_1.zbox, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&boxes_1.alpha, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&boxes_1.beta, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&boxes_1.gamma, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&boxes_1.volbox, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, value, (ftnlen)20);
	e_wsfe();
	if (s_cmp(boxes_1.spacegrp, "          ", (ftnlen)10, (ftnlen)10) != 
		0) {
	    s_copy(value, boxes_1.spacegrp, (ftnlen)20, (ftnlen)10);
	    justify_(value, (ftnlen)20);
	    io___32.ciunit = iounit_1.iout;
	    s_wsfe(&io___32);
	    do_fio(&c__1, value, (ftnlen)20);
	    e_wsfe();
	}
    }

/*     info on force field potential energy function */

    s_copy(value, fields_1.forcefield, (ftnlen)20, (ftnlen)20);
    justify_(value, (ftnlen)20);
    io___33.ciunit = iounit_1.iout;
    s_wsfe(&io___33);
    do_fio(&c__1, value, (ftnlen)20);
    e_wsfe();

/*     details of vdw potential energy functional form */

    if (potent_1.use_vdw__) {
	s_copy(label_ref(0, 1), vdwpot_1.vdwtyp, (ftnlen)20, (ftnlen)13);
	s_copy(label_ref(0, 2), vdwpot_1.radtyp, (ftnlen)20, (ftnlen)5);
	s_copy(label_ref(0, 3), vdwpot_1.radsiz, (ftnlen)20, (ftnlen)8);
	s_copy(label_ref(0, 4), vdwpot_1.radrule, (ftnlen)20, (ftnlen)10);
	s_copy(label_ref(0, 5), vdwpot_1.epsrule, (ftnlen)20, (ftnlen)10);
	for (i__ = 1; i__ <= 5; ++i__) {
	    justify_(label_ref(0, i__), (ftnlen)20);
	}
	io___36.ciunit = iounit_1.iout;
	s_wsfe(&io___36);
	for (i__ = 1; i__ <= 5; ++i__) {
	    do_fio(&c__1, label_ref(0, i__), (ftnlen)20);
	}
	e_wsfe();
	if (cutoff_1.vdwcut <= 1e3) {
	    io___37.ciunit = iounit_1.iout;
	    s_wsfe(&io___37);
	    do_fio(&c__1, (char *)&cutoff_1.vdwcut, (ftnlen)sizeof(doublereal)
		    );
	    e_wsfe();
	}
    }

/*     details of electrostatic energy functional form */

    if (potent_1.use_charge__ || potent_1.use_dipole__ || 
	    potent_1.use_mpole__ || potent_1.use_polar__) {
	io___38.ciunit = iounit_1.iout;
	s_wsfe(&io___38);
	e_wsfe();
    }
    if (potent_1.use_charge__) {
	s_copy(value, "PARTIAL CHARGE", (ftnlen)20, (ftnlen)14);
	justify_(value, (ftnlen)20);
	io___39.ciunit = iounit_1.iout;
	s_wsfe(&io___39);
	do_fio(&c__1, value, (ftnlen)20);
	e_wsfe();
    }
    if (potent_1.use_dipole__) {
	s_copy(value, "BOND DIPOLE", (ftnlen)20, (ftnlen)11);
	justify_(value, (ftnlen)20);
	io___40.ciunit = iounit_1.iout;
	s_wsfe(&io___40);
	do_fio(&c__1, value, (ftnlen)20);
	e_wsfe();
    }
    if (potent_1.use_mpole__) {
	s_copy(value, "ATOMIC MULTIPOLE", (ftnlen)20, (ftnlen)16);
	justify_(value, (ftnlen)20);
	io___41.ciunit = iounit_1.iout;
	s_wsfe(&io___41);
	do_fio(&c__1, value, (ftnlen)20);
	e_wsfe();
    }
    if (potent_1.use_polar__) {
	s_copy(value, "INDUCED DIPOLE", (ftnlen)20, (ftnlen)14);
	justify_(value, (ftnlen)20);
	io___42.ciunit = iounit_1.iout;
	s_wsfe(&io___42);
	do_fio(&c__1, value, (ftnlen)20);
	e_wsfe();
    }

/*     details of particle mesh Ewald calculation */

    if (cutoff_1.use_ewald__) {
	s_copy(value, ewald_1.boundary, (ftnlen)20, (ftnlen)7);
	justify_(value, (ftnlen)20);
	io___43.ciunit = iounit_1.iout;
	s_wsfe(&io___43);
	do_fio(&c__1, (char *)&ewald_1.aewald, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&cutoff_1.ewaldcut, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&pme_1.nfft1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&pme_1.nfft2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&pme_1.nfft3, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&pme_1.bsorder, (ftnlen)sizeof(integer));
	do_fio(&c__1, value, (ftnlen)20);
	e_wsfe();
    }
    return 0;
} /* systyze_ */

#undef label_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine paramyze  --  force field parameter analysis  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "paramyze" prints the force field parameters used in the */
/*     computation of each of the potential energy terms */


/* Subroutine */ int paramyze_(logical *active)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Interactions and Sites :\002,/)";
    static char fmt_20[] = "(\002 Bond Stretches\002,19x,i15)";
    static char fmt_30[] = "(\002 Angle Bends\002,22x,i15)";
    static char fmt_40[] = "(\002 Stretch-Bends\002,20x,i15)";
    static char fmt_50[] = "(\002 Urey-Bradley\002,21x,i15)";
    static char fmt_70[] = "(\002 Out-of-Plane Bends\002,15x,i15)";
    static char fmt_80[] = "(\002 Out-of-Plane Distances\002,11x,i15)";
    static char fmt_90[] = "(\002 Improper Dihedrals\002,15x,i15)";
    static char fmt_100[] = "(\002 Improper Torsions\002,16x,i15)";
    static char fmt_110[] = "(\002 Torsional Angles\002,17x,i15)";
    static char fmt_120[] = "(\002 Pi-Orbital Torsions\002,14x,i15)";
    static char fmt_130[] = "(\002 Stretch-Torsions\002,17x,i15)";
    static char fmt_140[] = "(\002 Torsion-Torsions\002,17x,i15)";
    static char fmt_150[] = "(\002 Van der Waals Sites\002,14x,i15)";
    static char fmt_160[] = "(\002 Atomic Partial Charges\002,11x,i15)";
    static char fmt_170[] = "(\002 Bond Dipole Moments\002,14x,i15)";
    static char fmt_180[] = "(\002 Atomic Multipoles\002,16x,i15)";
    static char fmt_190[] = "(\002 Polarizable Sites\002,16x,i15)";
    static char fmt_200[] = "(\002 Pisystem Atoms\002,19x,i15)";
    static char fmt_210[] = "(\002 Conjugated Pi-Bonds\002,14x,i15)";
    static char fmt_220[] = "(/,\002 Atom Type Definition Parameters :\002,/"
	    "/,3x,\002Atom\002,2x,\002Symbol\002,2x,\002Type\002,2x,\002Clas"
	    "s\002,2x,\002Atomic\002,3x,\002Mass\002,2x,\002Valence\002,2x"
	    ",\002Description\002,/)";
    static char fmt_230[] = "(i6,5x,a3,2i7,i6,f10.3,i5,5x,a24)";
    static char fmt_240[] = "(/,\002 Van der Waals Parameters :\002,//,10x"
	    ",\002Atom Number\002,7x,\002Size\002,3x,\002Epsilon\002,3x,\002S"
	    "ize 1-4\002,3x,\002Eps 1-4\002,3x,\002Reduction\002,/)";
    static char fmt_250[] = "(i6,3x,i6,7x,2f10.4)";
    static char fmt_260[] = "(i6,3x,i6,7x,2f10.4,22x,f10.4)";
    static char fmt_270[] = "(i6,3x,i6,7x,2f10.4,1x,2f10.4)";
    static char fmt_280[] = "(i6,3x,i6,7x,2f10.4,1x,2f10.4,1x,f10.4)";
    static char fmt_290[] = "(/,\002 Bond Stretching Parameters :\002,//,1"
	    "0x,\002Atom Numbers\002,25x,\002KS\002,7x,\002Bond\002,/)";
    static char fmt_300[] = "(i6,3x,2i6,19x,f10.3,f10.4)";
    static char fmt_310[] = "(/,\002 Angle Bending Parameters :\002,//,13x"
	    ",\002Atom Numbers\002,22x,\002KB\002,6x,\002Angle\002,3x,\002Fold"
	    "\002,4x,\002Type\002,/)";
    static char fmt_320[] = "(i6,3x,3i6,13x,2f10.3)";
    static char fmt_330[] = "(i6,3x,3i6,13x,2f10.3,9x,\002In-Plane\002)";
    static char fmt_340[] = "(i6,3x,3i6,13x,2f10.3,9x,\002Linear\002)";
    static char fmt_350[] = "(i6,3x,3i6,13x,2f10.3,f7.1,2x,\002Fourier\002)";
    static char fmt_360[] = "(/,\002 Stretch-Bend Parameters :\002,//,13x"
	    ",\002Atom Numbers\002,8x,\002KSB 1\002,5x,\002KSB 2\002,6x,\002A"
	    "ngle\002,3x,\002Bond 1\002,3x,\002Bond 2\002,/)";
    static char fmt_370[] = "(i6,3x,3i6,1x,2f10.3,2x,f9.3,2f9.4)";
    static char fmt_380[] = "(/,\002 Urey-Bradley Parameters :\002,//,13x"
	    ",\002Atom Numbers\002,21x,\002KUB\002,4x,\002Distance\002,/)";
    static char fmt_390[] = "(i6,3x,3i6,13x,f10.3,f10.4)";
    static char fmt_400[] = "(/,\002 Out-of-Plane Bend Parameters :\002,//,1"
	    "7x,\002Atom Numbers\002,19x,\002KOPB\002,/)";
    static char fmt_410[] = "(i6,3x,4i6,9x,f10.3)";
    static char fmt_420[] = "(/,\002 Out-of-Plane Distance Parameters :\002,"
	    "//,17x,\002Atom Numbers\002,19x,\002KOPD\002,/)";
    static char fmt_430[] = "(i6,3x,4i6,9x,f10.3)";
    static char fmt_440[] = "(/,\002 Improper Dihedral Parameters :\002,//,1"
	    "7x,\002Atom Numbers\002,19x,\002KID\002,4x,\002Dihedral\002,/)";
    static char fmt_450[] = "(i6,3x,4i6,9x,2f10.4)";
    static char fmt_460[] = "(/,\002 Improper Torsion Parameters :\002,//,17"
	    "x,\002Atom Numbers\002,11x,\002Amplitude, Phase and Periodicit"
	    "y\002,/)";
    static char fmt_470[] = "(i6,3x,4i6)";
    static char fmt_480[] = "(i6,3x,4i6,10x,f10.3,f8.1,i4)";
    static char fmt_490[] = "(i6,3x,4i6,2x,2(f10.3,f6.1,i4))";
    static char fmt_500[] = "(i6,3x,4i6,4x,3(f8.3,i4,\002/\002,i1))";
    static char fmt_510[] = "(/,\002 Torsional Angle Parameters :\002,//,1"
	    "7x,\002Atom Numbers\002,11x,\002Amplitude, Phase and Periodicit"
	    "y\002,/)";
    static char fmt_520[] = "(i6,3x,4i6)";
    static char fmt_530[] = "(i6,3x,4i6,4x,6(f8.3,i4,\002/\002,i1))";
    static char fmt_540[] = "(/,\002 Pi-Orbital Torsion Parameters :\002,//,"
	    "10x,\002Atom Numbers\002,19x,\002Amplitude\002,/)";
    static char fmt_550[] = "(i6,3x,2i6,19x,f10.4)";
    static char fmt_560[] = "(/,\002 Stretch-Torsion Parameters :\002,//,1"
	    "7x,\002Atom Numbers\002,10x,\002Length\002,5x,\002Torsion Term"
	    "s\002,/)";
    static char fmt_570[] = "(i6,3x,4i6,2x,f10.4,1x,3(f8.3,i4,\002/\002,i1))";
    static char fmt_580[] = "(/,\002 Torsion-Torsion Parameters :\002,//,2"
	    "0x,\002Atom Numbers\002,18x,\002Spline Grid\002,/)";
    static char fmt_590[] = "(i6,3x,5i6,10x,2i6)";
    static char fmt_600[] = "(/,\002 Atomic Partial Charge Parameters :\002,"
	    "/,45x,\002Neighbor\002,3x,\002Cutoff\002,/,10x,\002Atom Numbe"
	    "r\002,13x,\002Charge\002,7x,\002Site\002,6x,\002Site\002,/)";
    static char fmt_610[] = "(i6,3x,i6,15x,f10.4)";
    static char fmt_620[] = "(i6,3x,i6,15x,f10.4,5x,i6,4x,i6)";
    static char fmt_630[] = "(/,\002 Bond Dipole Moment Parameters :\002,//,"
	    "10x,\002Atom Numbers\002,22x,\002Dipole\002,3x,\002Position\002,"
	    "/)";
    static char fmt_640[] = "(i6,3x,2i6,19x,f10.4,f10.3)";
    static char fmt_650[] = "(/,\002 Atomic Multipole Parameters :\002,//,11"
	    "x,\002Atom\002,3x,\002Z-Axis\002,1x,\002X-Axis\002,1x,\002Y-Axi"
	    "s\002,2x,\002Frame\002,11x,\002Multipole Moments\002,/)";
    static char fmt_660[] = "(i6,3x,i6,1x,i7,17x,a8,2x,f9.5,/,50x,3f9.5,/,50"
	    "x,f9.5,/,50x,2f9.5,/,50x,3f9.5)";
    static char fmt_670[] = "(i6,3x,i6,1x,2i7,10x,a8,2x,f9.5,/,50x,3f9.5,/,5"
	    "0x,f9.5,/,50x,2f9.5,/,50x,3f9.5)";
    static char fmt_680[] = "(i6,3x,i6,1x,3i7,3x,a8,2x,f9.5,/,50x,3f9.5,/,50"
	    "x,f9.5,/,50x,2f9.5,/,50x,3f9.5)";
    static char fmt_690[] = "(/,\002 Dipole Polarizability Parameters :\002,"
	    "//,10x,\002Atom Number\002,5x,\002Alpha\002,5x,\002Damp\002,6x"
	    ",\002Polarization Group\002,/)";
    static char fmt_700[] = "(i6,3x,i6,6x,f10.4,f9.3,3x,20i6)";
    static char fmt_710[] = "(/,\002 Empirical Solvation Parameters :\002,//"
	    ",10x,\002Atom Number\002,13x,\002Radius\002,3x,\002ASP Value\002"
	    ",/)";
    static char fmt_720[] = "(i6,3x,i6,15x,2f10.4)";
    static char fmt_730[] = "(/,\002 Conjugated Pi-Atom Parameters :\002,//,"
	    "10x,\002Atom Number\002,14x,\002Nelect\002,6x,\002Ionize\002,4x"
	    ",\002Repulsion\002,/)";
    static char fmt_740[] = "(i6,3x,i6,17x,f8.1,3x,f10.4,2x,f10.4)";
    static char fmt_750[] = "(/,\002 Conjugated Pi-Bond Parameters :\002,//,"
	    "10x,\002Atom Numbers\002,21x,\002K Slope\002,3x,\002L Slope\002,"
	    "/)";
    static char fmt_760[] = "(i6,3x,2i6,19x,2f10.4)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen),
	     s_cmp(char *, char *, ftnlen, ftnlen), i_dnnt(doublereal *);

    /* Local variables */
    static integer i__, j, k, ia, ib, ic, id, ie, ig;
    static doublereal bla, blc, mpl[13];
    static integer fold[6];
    static doublereal radj, rad4j, ampli[6], phase[6];
    static integer ixaxe, iyaxe, izaxe;
    static logical header;

    /* Fortran I/O blocks */
    static cilist io___44 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___46 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___48 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___49 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___52 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___54 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___55 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___56 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___57 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___58 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___59 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___60 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___61 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___62 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___64 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___68 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___73 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___76 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___77 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___80 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___81 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___83 = { 0, 0, 0, fmt_310, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_320, 0 };
    static cilist io___85 = { 0, 0, 0, fmt_330, 0 };
    static cilist io___86 = { 0, 0, 0, fmt_340, 0 };
    static cilist io___87 = { 0, 0, 0, fmt_350, 0 };
    static cilist io___88 = { 0, 0, 0, fmt_360, 0 };
    static cilist io___91 = { 0, 0, 0, fmt_370, 0 };
    static cilist io___92 = { 0, 0, 0, fmt_380, 0 };
    static cilist io___93 = { 0, 0, 0, fmt_390, 0 };
    static cilist io___95 = { 0, 0, 0, fmt_400, 0 };
    static cilist io___96 = { 0, 0, 0, fmt_410, 0 };
    static cilist io___97 = { 0, 0, 0, fmt_420, 0 };
    static cilist io___98 = { 0, 0, 0, fmt_430, 0 };
    static cilist io___99 = { 0, 0, 0, fmt_440, 0 };
    static cilist io___100 = { 0, 0, 0, fmt_450, 0 };
    static cilist io___101 = { 0, 0, 0, fmt_460, 0 };
    static cilist io___105 = { 0, 0, 0, fmt_470, 0 };
    static cilist io___106 = { 0, 0, 0, fmt_480, 0 };
    static cilist io___107 = { 0, 0, 0, fmt_490, 0 };
    static cilist io___108 = { 0, 0, 0, fmt_500, 0 };
    static cilist io___109 = { 0, 0, 0, fmt_510, 0 };
    static cilist io___110 = { 0, 0, 0, fmt_520, 0 };
    static cilist io___111 = { 0, 0, 0, fmt_530, 0 };
    static cilist io___114 = { 0, 0, 0, fmt_540, 0 };
    static cilist io___115 = { 0, 0, 0, fmt_550, 0 };
    static cilist io___116 = { 0, 0, 0, fmt_560, 0 };
    static cilist io___117 = { 0, 0, 0, fmt_570, 0 };
    static cilist io___118 = { 0, 0, 0, fmt_580, 0 };
    static cilist io___119 = { 0, 0, 0, fmt_590, 0 };
    static cilist io___120 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___121 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___122 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___123 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___124 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___125 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___130 = { 0, 0, 0, fmt_660, 0 };
    static cilist io___131 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___132 = { 0, 0, 0, fmt_680, 0 };
    static cilist io___133 = { 0, 0, 0, fmt_690, 0 };
    static cilist io___134 = { 0, 0, 0, fmt_700, 0 };
    static cilist io___135 = { 0, 0, 0, fmt_710, 0 };
    static cilist io___136 = { 0, 0, 0, fmt_720, 0 };
    static cilist io___137 = { 0, 0, 0, fmt_730, 0 };
    static cilist io___138 = { 0, 0, 0, fmt_740, 0 };
    static cilist io___139 = { 0, 0, 0, fmt_750, 0 };
    static cilist io___140 = { 0, 0, 0, fmt_760, 0 };



#define ip11_ref(a_1,a_2) polgrp_1.ip11[(a_2)*100 + a_1 - 101]
#define isb_ref(a_1,a_2) strbnd_1.isb[(a_2)*3 + a_1 - 4]
#define sbk_ref(a_1,a_2) strbnd_1.sbk[(a_2)*2 + a_1 - 3]
#define ist_ref(a_1,a_2) strtor_1.ist[(a_2)*2 + a_1 - 3]
#define itt_ref(a_1,a_2) tortor_1.itt[(a_2)*3 + a_1 - 4]
#define kst_ref(a_1,a_2) strtor_1.kst[(a_2)*3 + a_1 - 4]
#define ibnd_ref(a_1,a_2) bond_1.ibnd[(a_2)*2 + a_1 - 3]
#define iang_ref(a_1,a_2) angle_1.iang[(a_2)*4 + a_1 - 5]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define idpl_ref(a_1,a_2) dipole_1.idpl[(a_2)*2 + a_1 - 3]
#define iopd_ref(a_1,a_2) opdist_1.iopd[(a_2)*4 + a_1 - 5]
#define pole_ref(a_1,a_2) mpole_1.pole[(a_2)*13 + a_1 - 14]
#define ibpi_ref(a_1,a_2) piorbs_1.ibpi[(a_2)*3 + a_1 - 4]
#define ipit_ref(a_1,a_2) pitors_1.ipit[(a_2)*6 + a_1 - 7]
#define iury_ref(a_1,a_2) urey_1.iury[(a_2)*3 + a_1 - 4]
#define tors1_ref(a_1,a_2) tors_1.tors1[(a_2)*4 + a_1 - 5]
#define tors2_ref(a_1,a_2) tors_1.tors2[(a_2)*4 + a_1 - 5]
#define tors3_ref(a_1,a_2) tors_1.tors3[(a_2)*4 + a_1 - 5]
#define tors4_ref(a_1,a_2) tors_1.tors4[(a_2)*4 + a_1 - 5]
#define tors5_ref(a_1,a_2) tors_1.tors5[(a_2)*4 + a_1 - 5]
#define tors6_ref(a_1,a_2) tors_1.tors6[(a_2)*4 + a_1 - 5]
#define itors_ref(a_1,a_2) tors_1.itors[(a_2)*4 + a_1 - 5]
#define story_ref(a_0,a_1) &atmtyp_1.story[(a_1)*24 + a_0 - 24]
#define itors1_ref(a_1,a_2) imptor_1.itors1[(a_2)*4 + a_1 - 5]
#define itors2_ref(a_1,a_2) imptor_1.itors2[(a_2)*4 + a_1 - 5]
#define itors3_ref(a_1,a_2) imptor_1.itors3[(a_2)*4 + a_1 - 5]
#define ibitor_ref(a_1,a_2) bitor_1.ibitor[(a_2)*5 + a_1 - 6]
#define angtyp_ref(a_0,a_1) &angpot_1.angtyp[(a_1)*8 + a_0 - 8]
#define iiprop_ref(a_1,a_2) improp_1.iiprop[(a_2)*4 + a_1 - 5]
#define iitors_ref(a_1,a_2) imptor_1.iitors[(a_2)*4 + a_1 - 5]
#define polaxe_ref(a_0,a_1) &mpole_1.polaxe[(a_1)*8 + a_0 - 8]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  angang.i  --  angle-angle terms in current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     kaa       force constant for angle-angle cross terms */
/*     nangang   total number of angle-angle interactions */
/*     iaa       angle numbers used in each angle-angle term */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  angle.i  --  bond angles within the current structure  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     ak       harmonic angle force constant (kcal/mole/rad**2) */
/*     anat     ideal bond angle or phase shift angle (degrees) */
/*     afld     periodicity for Fourier bond angle term */
/*     nangle   total number of bond angles in the system */
/*     iang     numbers of the atoms in each bond angle */




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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  atmtyp.i  --  atomic properties for each current atom  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     mass      atomic weight for each atom in the system */
/*     tag       integer atom labels from input coordinates file */
/*     class     atom class number for each atom in the system */
/*     atomic    atomic number for each atom in the system */
/*     valence   valence number for each atom in the system */
/*     name      atom name for each atom in the system */
/*     story     descriptive type for each atom in system */




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
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  bitor.i  --  bitorsions within the current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     nbitor  total number of bitorsions in the system */
/*     ibitor  numbers of the atoms in each bitorsion */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  bond.i  --  covalent bonds in the current structure  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     bk      bond stretch force constants (kcal/mole/Ang**2) */
/*     bl      ideal bond length values in Angstroms */
/*     nbond   total number of bond stretches in the system */
/*     ibnd    numbers of the atoms in each bond stretch */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  charge.i  --  partial charges for the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     pchg      magnitude of the partial charges (e-) */
/*     nion      total number of partial charges in system */
/*     iion      number of the atom site for each partial charge */
/*     jion      neighbor generation site for each partial charge */
/*     kion      cutoff switching site for each partial charge */
/*     chglist   partial charge site for each atom (0=no charge) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  dipole.i  --  atom & bond dipoles for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     bdpl      magnitude of each of the dipoles (Debyes) */
/*     sdpl      position of each dipole between defining atoms */
/*     ndipole   total number of dipoles in the system */
/*     idpl      numbers of atoms that define each dipole */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  improp.i  --  improper dihedrals in the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     kprop    force constant values for improper dihedral angles */
/*     vprop    ideal improper dihedral angle value in degrees */
/*     niprop   total number of improper dihedral angles in the system */
/*     iiprop   numbers of the atoms in each improper dihedral angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  imptor.i  --  improper torsions in the current structure  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     itors1   1-fold amplitude and phase for each improper torsion */
/*     itors2   2-fold amplitude and phase for each improper torsion */
/*     itors3   3-fold amplitude and phase for each improper torsion */
/*     nitors   total number of improper torsional angles in the system */
/*     iitors   numbers of the atoms in each improper torsional angle */




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
/*     ##  mpole.i  --  multipole components for current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     maxpole   max components (monopole=1,dipole=4,quadrupole=13) */

/*     pole      multipole values for each site in the local frame */
/*     rpole     multipoles rotated to the global coordinate system */
/*     npole     total number of multipole sites in the system */
/*     ipole     number of the atom for each multipole site */
/*     polsiz    number of multipole components at each atom */
/*     pollist   multipole site for each atom (0=no multipole) */
/*     zaxis     number of the z-axis defining atom for each site */
/*     xaxis     number of the x-axis defining atom for each site */
/*     yaxis     number of the y-axis defining atom for each site */
/*     polaxe    local axis type for each multipole site */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  opbend.i  --  out-of-plane bends in the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     opbk      force constant values for out-of-plane bending */
/*     nopbend   total number of out-of-plane bends in the system */
/*     iopb      bond angle numbers used in out-of-plane bending */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  opdist.i  --  out-of-plane distances in current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     opdk      force constant values for out-of-plane distance */
/*     nopdist   total number of out-of-plane distances in the system */
/*     iopb      numbers of the atoms in each out-of-plane distance */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  piorbs.i  --  conjugated system in the current structure  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     norbit    total number of pisystem orbitals in the system */
/*     iorbit    numbers of the atoms containing pisystem orbitals */
/*     reorbit   number of evaluations between orbital updates */
/*     piperp    atoms defining a normal plane to each orbital */
/*     nbpi      total number of bonds affected by the pisystem */
/*     ibpi      bond and piatom numbers for each pisystem bond */
/*     ntpi      total number of torsions affected by the pisystem */
/*     itpi      torsion and pibond numbers for each pisystem torsion */
/*     listpi    atom list indicating whether each atom has an orbital */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  pistuf.i  --  bonds and torsions in the current pisystem  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     bkpi     bond stretch force constants for pi-bond order of 1.0 */
/*     blpi     ideal bond length values for a pi-bond order of 1.0 */
/*     kslope   rate of force constant decrease with bond order decrease */
/*     lslope   rate of bond length increase with a bond order decrease */
/*     torsp2   2-fold torsional energy barrier for pi-bond order of 1.0 */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  pitors.i  --  pi-orbital torsions in the current structure  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     kpit      2-fold pi-orbital torsional force constants */
/*     npitors   total number of pi-orbital torsional interactions */
/*     ipit      numbers of the atoms in each pi-orbital torsion */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  polar.i  --  polarizabilities and induced dipole moments  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     polarity  dipole polarizability for each multipole site (Ang**3) */
/*     thole     Thole polarizability damping value for each site */
/*     pdamp     value of polarizability scale factor for each site */
/*     uind      induced dipole components at each multipole site */
/*     uinp      induced dipoles in field used for energy interactions */
/*     uinds     GK or PB induced dipoles at each multipole site */
/*     uinps     induced dipoles in field used for GK or PB energy */
/*     npolar    total number of polarizable sites in the system */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  polgrp.i  --  polarizable site group connectivity lists   ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     maxp11   maximum number of atoms in a polarization group */
/*     maxp12   maximum number of atoms in groups 1-2 to an atom */
/*     maxp13   maximum number of atoms in groups 1-3 to an atom */
/*     maxp14   maximum number of atoms in groups 1-4 to an atom */

/*     np11     number of atoms in polarization group of each atom */
/*     ip11     atom numbers of atoms in same group as each atom */
/*     np12     number of atoms in groups 1-2 to each atom */
/*     ip12     atom numbers of atoms in groups 1-2 to each atom */
/*     np13     number of atoms in groups 1-3 to each atom */
/*     ip13     atom numbers of atoms in groups 1-3 to each atom */
/*     np14     number of atoms in groups 1-4 to each atom */
/*     ip14     atom numbers of atoms in groups 1-4 to each atom */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  potent.i  --  usage of each potential energy component  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     use_bond    logical flag governing use of bond stretch potential */
/*     use_angle   logical flag governing use of angle bend potential */
/*     use_strbnd  logical flag governing use of stretch-bend potential */
/*     use_urey    logical flag governing use of Urey-Bradley potential */
/*     use_angang  logical flag governing use of angle-angle cross term */
/*     use_opbend  logical flag governing use of out-of-plane bend term */
/*     use_opdist  logical flag governing use of out-of-plane distance */
/*     use_improp  logical flag governing use of improper dihedral term */
/*     use_imptor  logical flag governing use of improper torsion term */
/*     use_tors    logical flag governing use of torsional potential */
/*     use_pitors  logical flag governing use of pi-orbital torsion term */
/*     use_strtor  logical flag governing use of stretch-torsion term */
/*     use_tortor  logical flag governing use of torsion-torsion term */
/*     use_vdw     logical flag governing use of vdw der Waals potential */
/*     use_charge  logical flag governing use of charge-charge potential */
/*     use_chgdpl  logical flag governing use of charge-dipole potential */
/*     use_dipole  logical flag governing use of dipole-dipole potential */
/*     use_mpole   logical flag governing use of multipole potential */
/*     use_polar   logical flag governing use of polarization term */
/*     use_rxnfld  logical flag governing use of reaction field term */
/*     use_solv    logical flag governing use of continuum solvation */
/*     use_metal   logical flag governing use of ligand field term */
/*     use_geom    logical flag governing use of geometric restraints */
/*     use_extra   logical flag governing use of extra potential term */
/*     use_born    logical flag governing use of Born radii values */
/*     use_orbit   logical flag governing use of pisystem computation */




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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  strbnd.i  --  stretch-bends in the current structure  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     sbk       force constants for stretch-bend terms */
/*     nstrbnd   total number of stretch-bend interactions */
/*     isb       angle and bond numbers used in stretch-bend */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  strtor.i  --  stretch-torsions in the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     kst       1-, 2- and 3-fold stretch-torsion force constants */
/*     nstrtor   total number of stretch-torsion interactions */
/*     ist       torsion and bond numbers used in stretch-torsion */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  tors.i  --  torsional angles within the current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     tors1   1-fold amplitude and phase for each torsional angle */
/*     tors2   2-fold amplitude and phase for each torsional angle */
/*     tors3   3-fold amplitude and phase for each torsional angle */
/*     tors4   4-fold amplitude and phase for each torsional angle */
/*     tors5   5-fold amplitude and phase for each torsional angle */
/*     tors6   6-fold amplitude and phase for each torsional angle */
/*     ntors   total number of torsional angles in the system */
/*     itors   numbers of the atoms in each torsional angle */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  tortor.i  --  torsion-torsions in the current structure  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     ntortor   total number of torsion-torsion interactions */
/*     itt       atoms and parameter indices for torsion-torsion */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  units.i  --  physical constants and unit conversions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     literature reference: */

/*     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended */
/*     Values of the Fundamental Physical Constants: 2006", Reviews of */
/*     Modern Physics, 80, 633-730 (2008) */

/*     The "2006 CODATA Recommended Values" are also available from */
/*     the NIST Reference on Constants, Units, and Uncertainty site */
/*     at http://physics.nist.gov/cuu/index.html */

/*     Most values below are derived from the 2006 CODATA reference */

/*     The conversion from calorie to Joule is the definition of the */
/*     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992) */

/*     The "coulomb" energy conversion factor is found by dimensional */
/*     analysis of Coulomb's Law, ie, by dividing the square of the */
/*     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is */
/*     the permittivity of vacuum (the "electric constant"); note that */
/*     eps0 is typically given in F/m, equivalent to C**2/(J-m) */

/*     The approximate value used for the Debye, 3.33564 x 10-30 C-m, */
/*     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997) */

/*     The value of "prescon" is based on definition of 1 atmosphere */
/*     as 101325 Pa set by the 10th Conference Generale des Poids et */
/*     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3 */

/*     avogadro    Avogadro's number (N) in particles/mole */
/*     lightspd    speed of light in vacuum (c) in cm/ps */
/*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K */
/*     gasconst    ideal gas constant (R) in kcal/mole/K */
/*     emass       mass of an electron in atomic mass units */
/*     joule       conversion from calories to joules */
/*     convert     conversion from kcal to g*Ang**2/ps**2 */
/*     bohr        conversion from Bohrs to Angstroms */
/*     hartree     conversion from Hartree to kcal/mole */
/*     evolt       conversion from Hartree to electron-volts */
/*     efreq       conversion from Hartree to cm-1 */
/*     coulomb     conversion from electron**2/Ang to kcal/mole */
/*     debye       conversion from electron-Ang to Debyes */
/*     prescon     conversion from kcal/mole/Ang**3 to Atm */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  urey.i  --  Urey-Bradley interactions in the structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     uk      Urey-Bradley force constants (kcal/mole/Ang**2) */
/*     ul      ideal 1-3 distance values in Angstroms */
/*     nurey   total number of Urey-Bradley terms in the system */
/*     iury    numbers of the atoms in each Urey-Bradley interaction */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  vdw.i  --  van der Waals parameters for current structure  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     radmin     minimum energy distance for each atom class pair */
/*     epsilon    well depth parameter for each atom class pair */
/*     radmin4    minimum energy distance for 1-4 interaction pairs */
/*     epsilon4   well depth parameter for 1-4 interaction pairs */
/*     radhbnd    minimum energy distance for hydrogen bonding pairs */
/*     epshbnd    well depth parameter for hydrogen bonding pairs */
/*     kred       value of reduction factor parameter for each atom */
/*     ired       attached atom from which reduction factor is applied */
/*     nvdw       total number van der Waals active sites in the system */
/*     ivdw       number of the atom for each van der Waals active site */
/*     jvdw       type or class index into vdw parameters for each atom */
/*     nvt        number of distinct van der Waals types in the system */
/*     ivt        number of each distinct vdw type/class in the system */
/*     jvt        frequency of each vdw type or class in the system */




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




/*     number of each type of interaction and site */

    /* Parameter adjustments */
    --active;

    /* Function Body */
    if (atoms_1.n != 0) {
	io___44.ciunit = iounit_1.iout;
	s_wsfe(&io___44);
	e_wsfe();
    }
    if (bond_1.nbond != 0) {
	io___45.ciunit = iounit_1.iout;
	s_wsfe(&io___45);
	do_fio(&c__1, (char *)&bond_1.nbond, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (angle_1.nangle != 0) {
	io___46.ciunit = iounit_1.iout;
	s_wsfe(&io___46);
	do_fio(&c__1, (char *)&angle_1.nangle, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (strbnd_1.nstrbnd != 0) {
	io___47.ciunit = iounit_1.iout;
	s_wsfe(&io___47);
	do_fio(&c__1, (char *)&strbnd_1.nstrbnd, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (urey_1.nurey != 0) {
	io___48.ciunit = iounit_1.iout;
	s_wsfe(&io___48);
	do_fio(&c__1, (char *)&urey_1.nurey, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (angang_1.nangang != 0) {
	io___49.ciunit = iounit_1.iout;
	s_wsfe(&io___49);
	do_fio(&c__1, (char *)&angang_1.nangang, (ftnlen)sizeof(integer));
	e_wsfe();
/* L60: */
    }
    if (opbend_1.nopbend != 0) {
	io___50.ciunit = iounit_1.iout;
	s_wsfe(&io___50);
	do_fio(&c__1, (char *)&opbend_1.nopbend, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (opdist_1.nopdist != 0) {
	io___51.ciunit = iounit_1.iout;
	s_wsfe(&io___51);
	do_fio(&c__1, (char *)&opdist_1.nopdist, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (improp_1.niprop != 0) {
	io___52.ciunit = iounit_1.iout;
	s_wsfe(&io___52);
	do_fio(&c__1, (char *)&improp_1.niprop, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (imptor_1.nitors != 0) {
	io___53.ciunit = iounit_1.iout;
	s_wsfe(&io___53);
	do_fio(&c__1, (char *)&imptor_1.nitors, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (tors_1.ntors != 0) {
	io___54.ciunit = iounit_1.iout;
	s_wsfe(&io___54);
	do_fio(&c__1, (char *)&tors_1.ntors, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (pitors_1.npitors != 0) {
	io___55.ciunit = iounit_1.iout;
	s_wsfe(&io___55);
	do_fio(&c__1, (char *)&pitors_1.npitors, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (strtor_1.nstrtor != 0) {
	io___56.ciunit = iounit_1.iout;
	s_wsfe(&io___56);
	do_fio(&c__1, (char *)&strtor_1.nstrtor, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (tortor_1.ntortor != 0) {
	io___57.ciunit = iounit_1.iout;
	s_wsfe(&io___57);
	do_fio(&c__1, (char *)&tortor_1.ntortor, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (vdw_1.nvdw != 0) {
	io___58.ciunit = iounit_1.iout;
	s_wsfe(&io___58);
	do_fio(&c__1, (char *)&vdw_1.nvdw, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (charge_1.nion != 0) {
	io___59.ciunit = iounit_1.iout;
	s_wsfe(&io___59);
	do_fio(&c__1, (char *)&charge_1.nion, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (dipole_1.ndipole != 0) {
	io___60.ciunit = iounit_1.iout;
	s_wsfe(&io___60);
	do_fio(&c__1, (char *)&dipole_1.ndipole, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (mpole_1.npole != 0) {
	io___61.ciunit = iounit_1.iout;
	s_wsfe(&io___61);
	do_fio(&c__1, (char *)&mpole_1.npole, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (polar_1.npolar != 0) {
	io___62.ciunit = iounit_1.iout;
	s_wsfe(&io___62);
	do_fio(&c__1, (char *)&polar_1.npolar, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (piorbs_1.norbit != 0) {
	io___63.ciunit = iounit_1.iout;
	s_wsfe(&io___63);
	do_fio(&c__1, (char *)&piorbs_1.norbit, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (piorbs_1.nbpi != 0) {
	io___64.ciunit = iounit_1.iout;
	s_wsfe(&io___64);
	do_fio(&c__1, (char *)&piorbs_1.nbpi, (ftnlen)sizeof(integer));
	e_wsfe();
    }

/*     parameters used for molecular mechanics atom types */

    header = TRUE_;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (active[i__]) {
	    if (header) {
		header = FALSE_;
		io___67.ciunit = iounit_1.iout;
		s_wsfe(&io___67);
		e_wsfe();
	    }
	    io___68.ciunit = iounit_1.iout;
	    s_wsfe(&io___68);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
	    do_fio(&c__1, (char *)&atoms_1.type__[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&atmtyp_1.class__[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&atmtyp_1.atomic[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&atmtyp_1.mass[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&atmtyp_1.valence[i__ - 1], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, story_ref(0, i__), (ftnlen)24);
	    e_wsfe();
	}
    }

/*     parameters used for van der Waals interactions */

    if (potent_1.use_vdw__) {
	header = TRUE_;
	k = 0;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (active[i__]) {
		++k;
		if (header) {
		    header = FALSE_;
		    io___70.ciunit = iounit_1.iout;
		    s_wsfe(&io___70);
		    e_wsfe();
		}
		j = atmtyp_1.class__[i__ - 1];
		if (kvdws_1.rad[j - 1] == kvdws_1.rad4[j - 1] && kvdws_1.eps[
			j - 1] == kvdws_1.eps4[j - 1]) {
		    radj = kvdws_1.rad[j - 1];
		    if (s_cmp(vdwpot_1.radsiz, "DIAMETER", (ftnlen)8, (ftnlen)
			    8) == 0) {
			radj *= 2.;
		    }
		    if (s_cmp(vdwpot_1.radtyp, "SIGMA", (ftnlen)5, (ftnlen)5) 
			    == 0) {
			radj /= 1.122462048309372981;
		    }
		    if (kvdws_1.reduct[j - 1] == 0.) {
			io___73.ciunit = iounit_1.iout;
			s_wsfe(&io___73);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&radj, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&kvdws_1.eps[j - 1], (ftnlen)
				sizeof(doublereal));
			e_wsfe();
		    } else {
			io___74.ciunit = iounit_1.iout;
			s_wsfe(&io___74);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&radj, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&kvdws_1.eps[j - 1], (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&kvdws_1.reduct[j - 1], (ftnlen)
				sizeof(doublereal));
			e_wsfe();
		    }
		} else {
		    radj = kvdws_1.rad[j - 1];
		    rad4j = kvdws_1.rad4[j - 1];
		    if (s_cmp(vdwpot_1.radsiz, "DIAMETER", (ftnlen)8, (ftnlen)
			    8) == 0) {
			radj *= 2.;
			rad4j *= 2.;
		    }
		    if (s_cmp(vdwpot_1.radtyp, "SIGMA", (ftnlen)5, (ftnlen)5) 
			    == 0) {
			radj /= 1.122462048309372981;
			rad4j /= 1.122462048309372981;
		    }
		    if (kvdws_1.reduct[j - 1] == 0.) {
			io___76.ciunit = iounit_1.iout;
			s_wsfe(&io___76);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&radj, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&kvdws_1.eps[j - 1], (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&rad4j, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&kvdws_1.eps4[j - 1], (ftnlen)
				sizeof(doublereal));
			e_wsfe();
		    } else {
			io___77.ciunit = iounit_1.iout;
			s_wsfe(&io___77);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&radj, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&kvdws_1.eps[j - 1], (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&rad4j, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&kvdws_1.eps4[j - 1], (ftnlen)
				sizeof(doublereal));
			do_fio(&c__1, (char *)&kvdws_1.reduct[j - 1], (ftnlen)
				sizeof(doublereal));
			e_wsfe();
		    }
		}
	    }
	}
    }

/*     parameters used for bond stretching interactions */

    if (potent_1.use_bond__) {
	header = TRUE_;
	i__1 = bond_1.nbond;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = ibnd_ref(1, i__);
	    ib = ibnd_ref(2, i__);
	    if (active[ia] || active[ib]) {
		if (header) {
		    header = FALSE_;
		    io___80.ciunit = iounit_1.iout;
		    s_wsfe(&io___80);
		    e_wsfe();
		}
		io___81.ciunit = iounit_1.iout;
		s_wsfe(&io___81);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&bond_1.bk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&bond_1.bl[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     parameters used for angle bending interactions */

    if (potent_1.use_angle__) {
	header = TRUE_;
	i__1 = angle_1.nangle;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = iang_ref(1, i__);
	    ib = iang_ref(2, i__);
	    ic = iang_ref(3, i__);
	    if (active[ia] || active[ib] || active[ic]) {
		if (header) {
		    header = FALSE_;
		    io___83.ciunit = iounit_1.iout;
		    s_wsfe(&io___83);
		    e_wsfe();
		}
		if (s_cmp(angtyp_ref(0, i__), "HARMONIC", (ftnlen)8, (ftnlen)
			8) == 0) {
		    io___84.ciunit = iounit_1.iout;
		    s_wsfe(&io___84);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&angle_1.ak[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&angle_1.anat[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    e_wsfe();
		} else if (s_cmp(angtyp_ref(0, i__), "IN-PLANE", (ftnlen)8, (
			ftnlen)8) == 0) {
		    io___85.ciunit = iounit_1.iout;
		    s_wsfe(&io___85);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&angle_1.ak[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&angle_1.anat[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    e_wsfe();
		} else if (s_cmp(angtyp_ref(0, i__), "IN-PLANE", (ftnlen)8, (
			ftnlen)8) == 0) {
		    io___86.ciunit = iounit_1.iout;
		    s_wsfe(&io___86);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&angle_1.ak[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&angle_1.anat[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    e_wsfe();
		} else if (s_cmp(angtyp_ref(0, i__), "FOURIER ", (ftnlen)8, (
			ftnlen)8) == 0) {
		    io___87.ciunit = iounit_1.iout;
		    s_wsfe(&io___87);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&angle_1.ak[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&angle_1.anat[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&angle_1.afld[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    e_wsfe();
		}
	    }
	}
    }

/*     parameters used for stretch-bend interactions */

    if (potent_1.use_strbnd__) {
	header = TRUE_;
	i__1 = strbnd_1.nstrbnd;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = isb_ref(1, i__);
	    ia = iang_ref(1, k);
	    ib = iang_ref(2, k);
	    ic = iang_ref(3, k);
	    if (active[ia] || active[ib] || active[ic]) {
		if (header) {
		    header = FALSE_;
		    io___88.ciunit = iounit_1.iout;
		    s_wsfe(&io___88);
		    e_wsfe();
		}
		bla = 0.;
		blc = 0.;
		if (isb_ref(2, i__) != 0) {
		    bla = bond_1.bl[isb_ref(2, i__) - 1];
		}
		if (isb_ref(3, i__) != 0) {
		    blc = bond_1.bl[isb_ref(3, i__) - 1];
		}
		io___91.ciunit = iounit_1.iout;
		s_wsfe(&io___91);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&sbk_ref(1, i__), (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&sbk_ref(2, i__), (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&angle_1.anat[k - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&bla, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&blc, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
    }

/*     parameters used for Urey-Bradley interactions */

    if (potent_1.use_urey__) {
	header = TRUE_;
	i__1 = urey_1.nurey;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = iury_ref(1, i__);
	    ib = iury_ref(2, i__);
	    ic = iury_ref(3, i__);
	    if (active[ia] || active[ic]) {
		if (header) {
		    header = FALSE_;
		    io___92.ciunit = iounit_1.iout;
		    s_wsfe(&io___92);
		    e_wsfe();
		}
		io___93.ciunit = iounit_1.iout;
		s_wsfe(&io___93);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&urey_1.uk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&urey_1.ul[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     parameters used for out-of-plane bend interactions */

    if (potent_1.use_opbend__) {
	header = TRUE_;
	i__1 = opbend_1.nopbend;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = opbend_1.iopb[i__ - 1];
	    ia = iang_ref(1, k);
	    ib = iang_ref(2, k);
	    ic = iang_ref(3, k);
	    id = iang_ref(4, k);
	    if (active[ia] || active[ib] || active[ic] || active[id]) {
		if (header) {
		    header = FALSE_;
		    io___95.ciunit = iounit_1.iout;
		    s_wsfe(&io___95);
		    e_wsfe();
		}
		io___96.ciunit = iounit_1.iout;
		s_wsfe(&io___96);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&opbend_1.opbk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     parameters used for out-of-plane distance interactions */

    if (potent_1.use_opdist__) {
	header = TRUE_;
	i__1 = opdist_1.nopdist;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = iopd_ref(1, i__);
	    ib = iopd_ref(2, i__);
	    ic = iopd_ref(3, i__);
	    id = iopd_ref(4, i__);
	    if (active[ia] || active[ib] || active[ic] || active[id]) {
		if (header) {
		    header = FALSE_;
		    io___97.ciunit = iounit_1.iout;
		    s_wsfe(&io___97);
		    e_wsfe();
		}
		io___98.ciunit = iounit_1.iout;
		s_wsfe(&io___98);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&opdist_1.opdk[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     parameters used for improper dihedral interactions */

    if (potent_1.use_improp__) {
	header = TRUE_;
	i__1 = improp_1.niprop;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = iiprop_ref(1, i__);
	    ib = iiprop_ref(2, i__);
	    ic = iiprop_ref(3, i__);
	    id = iiprop_ref(4, i__);
	    if (active[ia] || active[ib] || active[ic] || active[id]) {
		if (header) {
		    header = FALSE_;
		    io___99.ciunit = iounit_1.iout;
		    s_wsfe(&io___99);
		    e_wsfe();
		}
		io___100.ciunit = iounit_1.iout;
		s_wsfe(&io___100);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&improp_1.kprop[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&improp_1.vprop[i__ - 1], (ftnlen)
			sizeof(doublereal));
		e_wsfe();
	    }
	}
    }

/*     parameters used for improper torsion interactions */

    if (potent_1.use_imptor__) {
	header = TRUE_;
	i__1 = imptor_1.nitors;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = iitors_ref(1, i__);
	    ib = iitors_ref(2, i__);
	    ic = iitors_ref(3, i__);
	    id = iitors_ref(4, i__);
	    if (active[ia] || active[ib] || active[ic] || active[id]) {
		if (header) {
		    header = FALSE_;
		    io___101.ciunit = iounit_1.iout;
		    s_wsfe(&io___101);
		    e_wsfe();
		}
		j = 0;
		if (itors1_ref(1, i__) != 0.) {
		    ++j;
		    fold[j - 1] = 1;
		    ampli[j - 1] = itors1_ref(1, i__);
		    phase[j - 1] = itors1_ref(2, i__);
		}
		if (itors2_ref(1, i__) != 0.) {
		    ++j;
		    fold[j - 1] = 2;
		    ampli[j - 1] = itors2_ref(1, i__);
		    phase[j - 1] = itors2_ref(2, i__);
		}
		if (itors3_ref(1, i__) != 0.) {
		    ++j;
		    fold[j - 1] = 3;
		    ampli[j - 1] = itors3_ref(1, i__);
		    phase[j - 1] = itors3_ref(2, i__);
		}
		if (j == 0) {
		    io___105.ciunit = iounit_1.iout;
		    s_wsfe(&io___105);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		    e_wsfe();
		} else if (j == 1) {
		    io___106.ciunit = iounit_1.iout;
		    s_wsfe(&io___106);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ampli[0], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&phase[0], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&fold[0], (ftnlen)sizeof(integer));
		    e_wsfe();
		} else if (j == 2) {
		    io___107.ciunit = iounit_1.iout;
		    s_wsfe(&io___107);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		    i__2 = j;
		    for (k = 1; k <= i__2; ++k) {
			do_fio(&c__1, (char *)&ampli[k - 1], (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&phase[k - 1], (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&fold[k - 1], (ftnlen)sizeof(
				integer));
		    }
		    e_wsfe();
		} else {
		    io___108.ciunit = iounit_1.iout;
		    s_wsfe(&io___108);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		    i__2 = j;
		    for (k = 1; k <= i__2; ++k) {
			do_fio(&c__1, (char *)&ampli[k - 1], (ftnlen)sizeof(
				doublereal));
			i__3 = i_dnnt(&phase[k - 1]);
			do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&fold[k - 1], (ftnlen)sizeof(
				integer));
		    }
		    e_wsfe();
		}
	    }
	}
    }

/*     parameters used for torsional interactions */

    if (potent_1.use_tors__) {
	header = TRUE_;
	i__1 = tors_1.ntors;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = itors_ref(1, i__);
	    ib = itors_ref(2, i__);
	    ic = itors_ref(3, i__);
	    id = itors_ref(4, i__);
	    if (active[ia] || active[ib] || active[ic] || active[id]) {
		if (header) {
		    header = FALSE_;
		    io___109.ciunit = iounit_1.iout;
		    s_wsfe(&io___109);
		    e_wsfe();
		}
		j = 0;
		if (tors1_ref(1, i__) != 0.) {
		    ++j;
		    fold[j - 1] = 1;
		    ampli[j - 1] = tors1_ref(1, i__);
		    phase[j - 1] = tors1_ref(2, i__);
		}
		if (tors2_ref(1, i__) != 0.) {
		    ++j;
		    fold[j - 1] = 2;
		    ampli[j - 1] = tors2_ref(1, i__);
		    phase[j - 1] = tors2_ref(2, i__);
		}
		if (tors3_ref(1, i__) != 0.) {
		    ++j;
		    fold[j - 1] = 3;
		    ampli[j - 1] = tors3_ref(1, i__);
		    phase[j - 1] = tors3_ref(2, i__);
		}
		if (tors4_ref(1, i__) != 0.) {
		    ++j;
		    fold[j - 1] = 4;
		    ampli[j - 1] = tors4_ref(1, i__);
		    phase[j - 1] = tors4_ref(2, i__);
		}
		if (tors5_ref(1, i__) != 0.) {
		    ++j;
		    fold[j - 1] = 5;
		    ampli[j - 1] = tors5_ref(1, i__);
		    phase[j - 1] = tors5_ref(2, i__);
		}
		if (tors6_ref(1, i__) != 0.) {
		    ++j;
		    fold[j - 1] = 6;
		    ampli[j - 1] = tors6_ref(1, i__);
		    phase[j - 1] = tors6_ref(2, i__);
		}
		if (j == 0) {
		    io___110.ciunit = iounit_1.iout;
		    s_wsfe(&io___110);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		    e_wsfe();
		} else {
		    io___111.ciunit = iounit_1.iout;
		    s_wsfe(&io___111);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		    i__3 = j;
		    for (k = 1; k <= i__3; ++k) {
			do_fio(&c__1, (char *)&ampli[k - 1], (ftnlen)sizeof(
				doublereal));
			i__2 = i_dnnt(&phase[k - 1]);
			do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&fold[k - 1], (ftnlen)sizeof(
				integer));
		    }
		    e_wsfe();
		}
	    }
	}
    }

/*     parameters used for pi-orbital torsion interactions */

    if (potent_1.use_pitors__) {
	header = TRUE_;
	i__1 = pitors_1.npitors;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = ipit_ref(1, i__);
	    ib = ipit_ref(2, i__);
	    ic = ipit_ref(3, i__);
	    id = ipit_ref(4, i__);
	    ie = ipit_ref(5, i__);
	    ig = ipit_ref(6, i__);
	    if (active[ia] || active[ib] || active[ic] || active[id] || 
		    active[ie] || active[ig]) {
		if (header) {
		    header = FALSE_;
		    io___114.ciunit = iounit_1.iout;
		    s_wsfe(&io___114);
		    e_wsfe();
		}
		io___115.ciunit = iounit_1.iout;
		s_wsfe(&io___115);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&pitors_1.kpit[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     parameters used for stretch-torsion interactions */

    if (potent_1.use_strtor__) {
	header = TRUE_;
	i__1 = strtor_1.nstrtor;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = ist_ref(1, i__);
	    ia = itors_ref(1, k);
	    ib = itors_ref(2, k);
	    ic = itors_ref(3, k);
	    id = itors_ref(4, k);
	    if (active[ia] || active[ib] || active[ic] || active[id]) {
		if (header) {
		    header = FALSE_;
		    io___116.ciunit = iounit_1.iout;
		    s_wsfe(&io___116);
		    e_wsfe();
		}
		j = 0;
		if (kst_ref(1, i__) != 0.) {
		    ++j;
		    fold[j - 1] = 1;
		    ampli[j - 1] = kst_ref(1, i__);
		    phase[j - 1] = tors1_ref(2, k);
		}
		if (kst_ref(2, i__) != 0.) {
		    ++j;
		    fold[j - 1] = 2;
		    ampli[j - 1] = kst_ref(2, i__);
		    phase[j - 1] = tors2_ref(2, k);
		}
		if (kst_ref(3, i__) != 0.) {
		    ++j;
		    fold[j - 1] = 3;
		    ampli[j - 1] = kst_ref(3, i__);
		    phase[j - 1] = tors3_ref(2, k);
		}
		io___117.ciunit = iounit_1.iout;
		s_wsfe(&io___117);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&bond_1.bl[ist_ref(2, i__) - 1], (
			ftnlen)sizeof(doublereal));
		i__2 = j;
		for (k = 1; k <= i__2; ++k) {
		    do_fio(&c__1, (char *)&ampli[k - 1], (ftnlen)sizeof(
			    doublereal));
		    i__3 = i_dnnt(&phase[k - 1]);
		    do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&fold[k - 1], (ftnlen)sizeof(
			    integer));
		}
		e_wsfe();
	    }
	}
    }

/*     parameters used for torsion-torsion interactions */

    if (potent_1.use_tortor__) {
	header = TRUE_;
	i__1 = tortor_1.ntortor;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = itt_ref(1, i__);
	    ia = ibitor_ref(1, k);
	    ib = ibitor_ref(2, k);
	    ic = ibitor_ref(3, k);
	    id = ibitor_ref(4, k);
	    ie = ibitor_ref(5, k);
	    if (active[ia] || active[ib] || active[ic] || active[id] || 
		    active[ie]) {
		if (header) {
		    header = FALSE_;
		    io___118.ciunit = iounit_1.iout;
		    s_wsfe(&io___118);
		    e_wsfe();
		}
		j = itt_ref(2, i__);
		io___119.ciunit = iounit_1.iout;
		s_wsfe(&io___119);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&id, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ie, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ktrtor_1.tnx[j - 1], (ftnlen)sizeof(
			integer));
		do_fio(&c__1, (char *)&ktrtor_1.tny[j - 1], (ftnlen)sizeof(
			integer));
		e_wsfe();
	    }
	}
    }

/*     parameters used for atomic partial charges */

    if (potent_1.use_charge__ || potent_1.use_chgdpl__) {
	header = TRUE_;
	i__1 = charge_1.nion;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = charge_1.iion[i__ - 1];
	    ib = charge_1.jion[i__ - 1];
	    ic = charge_1.kion[i__ - 1];
	    if (active[ia] || active[ic]) {
		if (header) {
		    header = FALSE_;
		    io___120.ciunit = iounit_1.iout;
		    s_wsfe(&io___120);
		    e_wsfe();
		}
		if (ia == ib && ia == ic) {
		    io___121.ciunit = iounit_1.iout;
		    s_wsfe(&io___121);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&charge_1.pchg[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    e_wsfe();
		} else {
		    io___122.ciunit = iounit_1.iout;
		    s_wsfe(&io___122);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&charge_1.pchg[i__ - 1], (ftnlen)
			    sizeof(doublereal));
		    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		    e_wsfe();
		}
	    }
	}
    }

/*     parameters used for bond dipole moments */

    if (potent_1.use_dipole__ || potent_1.use_chgdpl__) {
	header = TRUE_;
	i__1 = dipole_1.ndipole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = idpl_ref(1, i__);
	    ib = idpl_ref(2, i__);
	    if (active[ia] || active[ib]) {
		if (header) {
		    header = FALSE_;
		    io___123.ciunit = iounit_1.iout;
		    s_wsfe(&io___123);
		    e_wsfe();
		}
		io___124.ciunit = iounit_1.iout;
		s_wsfe(&io___124);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&dipole_1.bdpl[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&dipole_1.sdpl[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }

/*     parameters used for atomic multipole moments */

    if (potent_1.use_mpole__) {
	header = TRUE_;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = mpole_1.ipole[i__ - 1];
	    if (active[ia]) {
		if (header) {
		    header = FALSE_;
		    io___125.ciunit = iounit_1.iout;
		    s_wsfe(&io___125);
		    e_wsfe();
		}
		izaxe = mpole_1.zaxis[i__ - 1];
		ixaxe = mpole_1.xaxis[i__ - 1];
		iyaxe = mpole_1.yaxis[i__ - 1];
		if (izaxe > atoms_1.n) {
		    izaxe = 0;
		}
		if (ixaxe > atoms_1.n) {
		    ixaxe = 0;
		}
		if (iyaxe < 0) {
		    iyaxe = -iyaxe;
		}
		mpl[0] = pole_ref(1, i__);
		for (j = 2; j <= 4; ++j) {
		    mpl[j - 1] = pole_ref(j, i__) / .52917720859;
		}
		for (j = 5; j <= 13; ++j) {
		    mpl[j - 1] = pole_ref(j, i__) * 3. / .28002851809110429;
		}
		if (ixaxe == 0) {
		    io___130.ciunit = iounit_1.iout;
		    s_wsfe(&io___130);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&izaxe, (ftnlen)sizeof(integer));
		    do_fio(&c__1, polaxe_ref(0, i__), (ftnlen)8);
		    for (j = 1; j <= 5; ++j) {
			do_fio(&c__1, (char *)&mpl[j - 1], (ftnlen)sizeof(
				doublereal));
		    }
		    do_fio(&c__1, (char *)&mpl[7], (ftnlen)sizeof(doublereal))
			    ;
		    do_fio(&c__1, (char *)&mpl[8], (ftnlen)sizeof(doublereal))
			    ;
		    for (j = 11; j <= 13; ++j) {
			do_fio(&c__1, (char *)&mpl[j - 1], (ftnlen)sizeof(
				doublereal));
		    }
		    e_wsfe();
		} else if (iyaxe == 0) {
		    io___131.ciunit = iounit_1.iout;
		    s_wsfe(&io___131);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&izaxe, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ixaxe, (ftnlen)sizeof(integer));
		    do_fio(&c__1, polaxe_ref(0, i__), (ftnlen)8);
		    for (j = 1; j <= 5; ++j) {
			do_fio(&c__1, (char *)&mpl[j - 1], (ftnlen)sizeof(
				doublereal));
		    }
		    do_fio(&c__1, (char *)&mpl[7], (ftnlen)sizeof(doublereal))
			    ;
		    do_fio(&c__1, (char *)&mpl[8], (ftnlen)sizeof(doublereal))
			    ;
		    for (j = 11; j <= 13; ++j) {
			do_fio(&c__1, (char *)&mpl[j - 1], (ftnlen)sizeof(
				doublereal));
		    }
		    e_wsfe();
		} else {
		    io___132.ciunit = iounit_1.iout;
		    s_wsfe(&io___132);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&izaxe, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&ixaxe, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&iyaxe, (ftnlen)sizeof(integer));
		    do_fio(&c__1, polaxe_ref(0, i__), (ftnlen)8);
		    for (j = 1; j <= 5; ++j) {
			do_fio(&c__1, (char *)&mpl[j - 1], (ftnlen)sizeof(
				doublereal));
		    }
		    do_fio(&c__1, (char *)&mpl[7], (ftnlen)sizeof(doublereal))
			    ;
		    do_fio(&c__1, (char *)&mpl[8], (ftnlen)sizeof(doublereal))
			    ;
		    for (j = 11; j <= 13; ++j) {
			do_fio(&c__1, (char *)&mpl[j - 1], (ftnlen)sizeof(
				doublereal));
		    }
		    e_wsfe();
		}
	    }
	}
    }

/*     parameters used for dipole polarizability */

    if (potent_1.use_polar__) {
	header = TRUE_;
	i__1 = mpole_1.npole;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = mpole_1.ipole[i__ - 1];
	    if (active[ia]) {
		if (header) {
		    header = FALSE_;
		    io___133.ciunit = iounit_1.iout;
		    s_wsfe(&io___133);
		    e_wsfe();
		}
		io___134.ciunit = iounit_1.iout;
		s_wsfe(&io___134);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&polar_1.polarity[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&polar_1.thole[i__ - 1], (ftnlen)sizeof(
			doublereal));
		i__3 = polgrp_1.np11[ia - 1];
		for (j = 1; j <= i__3; ++j) {
		    do_fio(&c__1, (char *)&ip11_ref(j, ia), (ftnlen)sizeof(
			    integer));
		}
		e_wsfe();
	    }
	}
    }

/*     parameters used for empirical solvation */

    if (potent_1.use_solv__) {
	header = TRUE_;
	k = 0;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (active[i__]) {
		++k;
		if (header) {
		    header = FALSE_;
		    io___135.ciunit = iounit_1.iout;
		    s_wsfe(&io___135);
		    e_wsfe();
		}
		io___136.ciunit = iounit_1.iout;
		s_wsfe(&io___136);
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&solute_1.rsolv[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&solute_1.asolv[i__ - 1], (ftnlen)
			sizeof(doublereal));
		e_wsfe();
	    }
	}
    }

/*     parameters used for conjugated pisystem atoms */

    if (potent_1.use_orbit__) {
	header = TRUE_;
	i__1 = piorbs_1.norbit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = piorbs_1.iorbit[i__ - 1];
	    j = atmtyp_1.class__[ia - 1];
	    if (header) {
		header = FALSE_;
		io___137.ciunit = iounit_1.iout;
		s_wsfe(&io___137);
		e_wsfe();
	    }
	    io___138.ciunit = iounit_1.iout;
	    s_wsfe(&io___138);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&korbs_1.electron[j - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&korbs_1.ionize[j - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&korbs_1.repulse[j - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }

/*     parameters used for conjugated pibond interactions */

    if (potent_1.use_orbit__) {
	header = TRUE_;
	i__1 = piorbs_1.nbpi;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = ibpi_ref(2, i__);
	    ib = ibpi_ref(3, i__);
	    if (header) {
		header = FALSE_;
		io___139.ciunit = iounit_1.iout;
		s_wsfe(&io___139);
		e_wsfe();
	    }
	    io___140.ciunit = iounit_1.iout;
	    s_wsfe(&io___140);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&pistuf_1.kslope[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&pistuf_1.lslope[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	    e_wsfe();
	}
    }
    return 0;
} /* paramyze_ */

#undef polaxe_ref
#undef iitors_ref
#undef iiprop_ref
#undef angtyp_ref
#undef ibitor_ref
#undef itors3_ref
#undef itors2_ref
#undef itors1_ref
#undef story_ref
#undef itors_ref
#undef tors6_ref
#undef tors5_ref
#undef tors4_ref
#undef tors3_ref
#undef tors2_ref
#undef tors1_ref
#undef iury_ref
#undef ipit_ref
#undef ibpi_ref
#undef pole_ref
#undef iopd_ref
#undef idpl_ref
#undef name___ref
#undef iang_ref
#undef ibnd_ref
#undef kst_ref
#undef itt_ref
#undef ist_ref
#undef sbk_ref
#undef isb_ref
#undef ip11_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine enrgyze  --  compute & report energy analysis  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "enrgyze" is an auxiliary routine for the analyze program */
/*     that performs the energy analysis and prints the total and */
/*     intermolecular energies */


/* Subroutine */ int enrgyze_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Total Potential Energy :\002,4x,f20.8"
	    ",\002 Kcal/mole\002)";
    static char fmt_20[] = "(/,\002 Total Potential Energy :\002,4x,d20.8"
	    ",\002 Kcal/mole\002)";
    static char fmt_30[] = "(/,\002 Total Potential Energy :\002,6x,f18.6"
	    ",\002 Kcal/mole\002)";
    static char fmt_40[] = "(/,\002 Total Potential Energy :\002,6x,d18.6"
	    ",\002 Kcal/mole\002)";
    static char fmt_50[] = "(/,\002 Total Potential Energy :\002,8x,f16.4"
	    ",\002 Kcal/mole\002)";
    static char fmt_60[] = "(/,\002 Total Potential Energy :\002,8x,d16.4"
	    ",\002 Kcal/mole\002)";
    static char fmt_70[] = "(/,\002 Intermolecular Energy :\002,5x,f20.8,"
	    "\002 Kcal/mole\002)";
    static char fmt_80[] = "(/,\002 Intermolecular Energy :\002,5x,d20.8,"
	    "\002 Kcal/mole\002)";
    static char fmt_90[] = "(/,\002 Intermolecular Energy :\002,7x,f18.6,"
	    "\002 Kcal/mole\002)";
    static char fmt_100[] = "(/,\002 Intermolecular Energy :\002,7x,d18.6"
	    ",\002 Kcal/mole\002)";
    static char fmt_110[] = "(/,\002 Intermolecular Energy :\002,9x,f16.4"
	    ",\002 Kcal/mole\002)";
    static char fmt_120[] = "(/,\002 Intermolecular Energy :\002,9x,d16.4"
	    ",\002 Kcal/mole\002)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int analysis_(doublereal *);
    static doublereal energy;

    /* Fortran I/O blocks */
    static cilist io___142 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___143 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___144 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___145 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___146 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___147 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___148 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___149 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___150 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___151 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___152 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___153 = { 0, 0, 0, fmt_120, 0 };




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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cutoff.i  --  cutoff distances for energy interactions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     vdwcut      cutoff distance for van der Waals interactions */
/*     chgcut      cutoff distance for charge-charge interactions */
/*     dplcut      cutoff distance for dipole-dipole interactions */
/*     mpolecut    cutoff distance for atomic multipole interactions */
/*     vdwtaper    distance at which van der Waals switching begins */
/*     chgtaper    distance at which charge-charge switching begins */
/*     dpltaper    distance at which dipole-dipole switching begins */
/*     mpoletaper  distance at which atomic multipole switching begins */
/*     ewaldcut    cutoff distance for direct space Ewald summation */
/*     use_ewald   logical flag governing use of Ewald summation */
/*     use_lights  logical flag governing use of method of lights */
/*     use_list    logical flag governing use of any neighbor lists */
/*     use_vlist   logical flag governing use of vdw neighbor lists */
/*     use_clist   logical flag governing use of charge neighbor lists */
/*     use_mlist   logical flag governing use of multipole neighbor lists */




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

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  molcul.i  --  individual molecules within current system  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     molmass   molecular weight for each molecule in the system */
/*     totmass   total weight of all the molecules in the system */
/*     nmol      total number of separate molecules in the system */
/*     kmol      contiguous list of the atoms in each molecule */
/*     imol      first and last atom of each molecule in the list */
/*     molcule   number of the molecule to which each atom belongs */




/*     perform the energy analysis by atom and component */

    analysis_(&energy);

/*     print out the total potential energy of the system */

    if (inform_1.digits >= 8) {
	if (abs(energy) < 1e10) {
	    io___142.ciunit = iounit_1.iout;
	    s_wsfe(&io___142);
	    do_fio(&c__1, (char *)&energy, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___143.ciunit = iounit_1.iout;
	    s_wsfe(&io___143);
	    do_fio(&c__1, (char *)&energy, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    } else if (inform_1.digits >= 6) {
	if (abs(energy) < 1e10) {
	    io___144.ciunit = iounit_1.iout;
	    s_wsfe(&io___144);
	    do_fio(&c__1, (char *)&energy, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___145.ciunit = iounit_1.iout;
	    s_wsfe(&io___145);
	    do_fio(&c__1, (char *)&energy, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    } else {
	if (abs(energy) < 1e10) {
	    io___146.ciunit = iounit_1.iout;
	    s_wsfe(&io___146);
	    do_fio(&c__1, (char *)&energy, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___147.ciunit = iounit_1.iout;
	    s_wsfe(&io___147);
	    do_fio(&c__1, (char *)&energy, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     intermolecular energy for systems with multiple molecules */

    if (molcul_1.nmol > 1 && molcul_1.nmol < atoms_1.n && ! 
	    cutoff_1.use_ewald__) {
	if (inform_1.digits >= 8) {
	    if (abs(inter_1.einter) < 1e10) {
		io___148.ciunit = iounit_1.iout;
		s_wsfe(&io___148);
		do_fio(&c__1, (char *)&inter_1.einter, (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___149.ciunit = iounit_1.iout;
		s_wsfe(&io___149);
		do_fio(&c__1, (char *)&inter_1.einter, (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	} else if (inform_1.digits >= 6) {
	    if (abs(inter_1.einter) < 1e10) {
		io___150.ciunit = iounit_1.iout;
		s_wsfe(&io___150);
		do_fio(&c__1, (char *)&inter_1.einter, (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___151.ciunit = iounit_1.iout;
		s_wsfe(&io___151);
		do_fio(&c__1, (char *)&inter_1.einter, (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	} else {
	    if (abs(inter_1.einter) < 1e10) {
		io___152.ciunit = iounit_1.iout;
		s_wsfe(&io___152);
		do_fio(&c__1, (char *)&inter_1.einter, (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    } else {
		io___153.ciunit = iounit_1.iout;
		s_wsfe(&io___153);
		do_fio(&c__1, (char *)&inter_1.einter, (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }
    return 0;
} /* enrgyze_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine analyz4  --  low precision energy components  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "analyz4" prints the energy to 4 decimal places and number */
/*     of interactions for each component of the potential energy */


/* Subroutine */ int analyz4_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Energy Component Breakdown :\002,11x,"
	    "\002Kcal/mole\002,6x,\002Interactions\002/)";
    static char fmt_20[] = "(\002 Bond Stretching\002,17x,f16.4,i15)";
    static char fmt_30[] = "(\002 Angle Bending\002,19x,f16.4,i15)";
    static char fmt_40[] = "(\002 Stretch-Bend\002,20x,f16.4,i15)";
    static char fmt_50[] = "(\002 Urey-Bradley\002,20x,f16.4,i15)";
    static char fmt_60[] = "(\002 Angle-Angle\002,21x,f16.4,i15)";
    static char fmt_70[] = "(\002 Out-of-Plane Bend\002,15x,f16.4,i15)";
    static char fmt_80[] = "(\002 Out-of-Plane Distance\002,11x,f16.4,i15)";
    static char fmt_90[] = "(\002 Improper Dihedral\002,15x,f16.4,i15)";
    static char fmt_100[] = "(\002 Improper Torsion\002,16x,f16.4,i15)";
    static char fmt_110[] = "(\002 Torsional Angle\002,17x,f16.4,i15)";
    static char fmt_120[] = "(\002 Pi-Orbital Torsion\002,14x,f16.4,i15)";
    static char fmt_130[] = "(\002 Stretch-Torsion\002,17x,f16.4,i15)";
    static char fmt_140[] = "(\002 Torsion-Torsion\002,17x,f16.4,i15)";
    static char fmt_150[] = "(\002 Van der Waals\002,19x,f16.4,i15)";
    static char fmt_160[] = "(\002 Van der Waals\002,19x,d16.4,i15)";
    static char fmt_170[] = "(\002 Charge-Charge\002,19x,f16.4,i15)";
    static char fmt_180[] = "(\002 Charge-Charge\002,19x,d16.4,i15)";
    static char fmt_190[] = "(\002 Charge-Dipole\002,19x,f16.4,i15)";
    static char fmt_200[] = "(\002 Charge-Dipole\002,19x,d16.4,i15)";
    static char fmt_210[] = "(\002 Dipole-Dipole\002,19x,f16.4,i15)";
    static char fmt_220[] = "(\002 Dipole-Dipole\002,19x,d16.4,i15)";
    static char fmt_230[] = "(\002 Atomic Multipoles\002,15x,f16.4,i15)";
    static char fmt_240[] = "(\002 Atomic Multipoles\002,15x,d16.4,i15)";
    static char fmt_250[] = "(\002 Polarization\002,20x,f16.4,i15)";
    static char fmt_260[] = "(\002 Polarization\002,20x,d16.4,i15)";
    static char fmt_270[] = "(\002 Reaction Field\002,18x,f16.4,i15)";
    static char fmt_280[] = "(\002 Continuum Solvation\002,13x,f16.4,i15)";
    static char fmt_290[] = "(\002 Metal Ligand Field\002,14x,f16.4,i15)";
    static char fmt_300[] = "(\002 Geometric Restraints\002,12x,f16.4,i15)";
    static char fmt_310[] = "(\002 Extra Energy Terms\002,14x,f16.4,i15)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___154 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___155 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___156 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___157 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___158 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___159 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___160 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___161 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___162 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___163 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___164 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___165 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___166 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___167 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___168 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___169 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___170 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___171 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___172 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___173 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___174 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___175 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___176 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___177 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___178 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___179 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___180 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___181 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___182 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___183 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___184 = { 0, 0, 0, fmt_310, 0 };




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  action.i  --  total number of each energy term computed  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     neb     number of bond stretch energy terms computed */
/*     nea     number of angle bend energy terms computed */
/*     neba    number of stretch-bend energy terms computed */
/*     neub    number of Urey-Bradley energy terms computed */
/*     neaa    number of angle-angle energy terms computed */
/*     neopb   number of out-of-plane bend energy terms computed */
/*     neopd   number of out-of-plane distance energy terms computed */
/*     neid    number of improper dihedral energy terms computed */
/*     neit    number of improper torsion energy terms computed */
/*     net     number of torsional energy terms computed */
/*     nept    number of pi-orbital torsion energy terms computed */
/*     nebt    number of stretch-torsion energy terms computed */
/*     nett    number of torsion-torsion energy terms computed */
/*     nev     number of van der Waals energy terms computed */
/*     nec     number of charge-charge energy terms computed */
/*     necd    number of charge-dipole energy terms computed */
/*     ned     number of dipole-dipole energy terms computed */
/*     nem     number of multipole energy terms computed */
/*     nep     number of polarization energy terms computed */
/*     new     number of Ewald summation energy terms computed */
/*     ner     number of reaction field energy terms computed */
/*     nes     number of solvation energy terms computed */
/*     nelf    number of metal ligand field energy terms computed */
/*     neg     number of geometric restraint energy terms computed */
/*     nex     number of extra energy terms computed */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cutoff.i  --  cutoff distances for energy interactions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     vdwcut      cutoff distance for van der Waals interactions */
/*     chgcut      cutoff distance for charge-charge interactions */
/*     dplcut      cutoff distance for dipole-dipole interactions */
/*     mpolecut    cutoff distance for atomic multipole interactions */
/*     vdwtaper    distance at which van der Waals switching begins */
/*     chgtaper    distance at which charge-charge switching begins */
/*     dpltaper    distance at which dipole-dipole switching begins */
/*     mpoletaper  distance at which atomic multipole switching begins */
/*     ewaldcut    cutoff distance for direct space Ewald summation */
/*     use_ewald   logical flag governing use of Ewald summation */
/*     use_lights  logical flag governing use of method of lights */
/*     use_list    logical flag governing use of any neighbor lists */
/*     use_vlist   logical flag governing use of vdw neighbor lists */
/*     use_clist   logical flag governing use of charge neighbor lists */
/*     use_mlist   logical flag governing use of multipole neighbor lists */




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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     write out each energy component to the desired precision */



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  potent.i  --  usage of each potential energy component  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     use_bond    logical flag governing use of bond stretch potential */
/*     use_angle   logical flag governing use of angle bend potential */
/*     use_strbnd  logical flag governing use of stretch-bend potential */
/*     use_urey    logical flag governing use of Urey-Bradley potential */
/*     use_angang  logical flag governing use of angle-angle cross term */
/*     use_opbend  logical flag governing use of out-of-plane bend term */
/*     use_opdist  logical flag governing use of out-of-plane distance */
/*     use_improp  logical flag governing use of improper dihedral term */
/*     use_imptor  logical flag governing use of improper torsion term */
/*     use_tors    logical flag governing use of torsional potential */
/*     use_pitors  logical flag governing use of pi-orbital torsion term */
/*     use_strtor  logical flag governing use of stretch-torsion term */
/*     use_tortor  logical flag governing use of torsion-torsion term */
/*     use_vdw     logical flag governing use of vdw der Waals potential */
/*     use_charge  logical flag governing use of charge-charge potential */
/*     use_chgdpl  logical flag governing use of charge-dipole potential */
/*     use_dipole  logical flag governing use of dipole-dipole potential */
/*     use_mpole   logical flag governing use of multipole potential */
/*     use_polar   logical flag governing use of polarization term */
/*     use_rxnfld  logical flag governing use of reaction field term */
/*     use_solv    logical flag governing use of continuum solvation */
/*     use_metal   logical flag governing use of ligand field term */
/*     use_geom    logical flag governing use of geometric restraints */
/*     use_extra   logical flag governing use of extra potential term */
/*     use_born    logical flag governing use of Born radii values */
/*     use_orbit   logical flag governing use of pisystem computation */


    io___154.ciunit = iounit_1.iout;
    s_wsfe(&io___154);
    e_wsfe();
    if (potent_1.use_bond__ && action_1.neb != 0) {
	io___155.ciunit = iounit_1.iout;
	s_wsfe(&io___155);
	do_fio(&c__1, (char *)&energi_1.eb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neb, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_angle__ && action_1.nea != 0) {
	io___156.ciunit = iounit_1.iout;
	s_wsfe(&io___156);
	do_fio(&c__1, (char *)&energi_1.ea, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nea, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_strbnd__ && action_1.neba != 0) {
	io___157.ciunit = iounit_1.iout;
	s_wsfe(&io___157);
	do_fio(&c__1, (char *)&energi_1.eba, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neba, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_urey__ && action_1.neub != 0) {
	io___158.ciunit = iounit_1.iout;
	s_wsfe(&io___158);
	do_fio(&c__1, (char *)&energi_1.eub, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neub, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_angang__ && action_1.neaa != 0) {
	io___159.ciunit = iounit_1.iout;
	s_wsfe(&io___159);
	do_fio(&c__1, (char *)&energi_1.eaa, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neaa, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_opbend__ && action_1.neopb != 0) {
	io___160.ciunit = iounit_1.iout;
	s_wsfe(&io___160);
	do_fio(&c__1, (char *)&energi_1.eopb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neopb, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_opdist__ && action_1.neopd != 0) {
	io___161.ciunit = iounit_1.iout;
	s_wsfe(&io___161);
	do_fio(&c__1, (char *)&energi_1.eopd, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neopd, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_improp__ && action_1.neid != 0) {
	io___162.ciunit = iounit_1.iout;
	s_wsfe(&io___162);
	do_fio(&c__1, (char *)&energi_1.eid, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neid, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_imptor__ && action_1.neit != 0) {
	io___163.ciunit = iounit_1.iout;
	s_wsfe(&io___163);
	do_fio(&c__1, (char *)&energi_1.eit, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neit, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_tors__ && action_1.net != 0) {
	io___164.ciunit = iounit_1.iout;
	s_wsfe(&io___164);
	do_fio(&c__1, (char *)&energi_1.et, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.net, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_pitors__ && action_1.nept != 0) {
	io___165.ciunit = iounit_1.iout;
	s_wsfe(&io___165);
	do_fio(&c__1, (char *)&energi_1.ept, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nept, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_strtor__ && action_1.nebt != 0) {
	io___166.ciunit = iounit_1.iout;
	s_wsfe(&io___166);
	do_fio(&c__1, (char *)&energi_1.ebt, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nebt, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_tortor__ && action_1.nett != 0) {
	io___167.ciunit = iounit_1.iout;
	s_wsfe(&io___167);
	do_fio(&c__1, (char *)&energi_1.ett, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nett, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_vdw__ && action_1.nev != 0) {
	if (abs(energi_1.ev) < 1e10) {
	    io___168.ciunit = iounit_1.iout;
	    s_wsfe(&io___168);
	    do_fio(&c__1, (char *)&energi_1.ev, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nev, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___169.ciunit = iounit_1.iout;
	    s_wsfe(&io___169);
	    do_fio(&c__1, (char *)&energi_1.ev, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nev, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_charge__ && action_1.nec != 0) {
	if (abs(energi_1.ec) < 1e10) {
	    io___170.ciunit = iounit_1.iout;
	    s_wsfe(&io___170);
	    do_fio(&c__1, (char *)&energi_1.ec, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nec, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___171.ciunit = iounit_1.iout;
	    s_wsfe(&io___171);
	    do_fio(&c__1, (char *)&energi_1.ec, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nec, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_chgdpl__ && action_1.necd != 0) {
	if (abs(energi_1.ecd) < 1e10) {
	    io___172.ciunit = iounit_1.iout;
	    s_wsfe(&io___172);
	    do_fio(&c__1, (char *)&energi_1.ecd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.necd, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___173.ciunit = iounit_1.iout;
	    s_wsfe(&io___173);
	    do_fio(&c__1, (char *)&energi_1.ecd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.necd, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_dipole__ && action_1.ned != 0) {
	if (abs(energi_1.ed) < 1e10) {
	    io___174.ciunit = iounit_1.iout;
	    s_wsfe(&io___174);
	    do_fio(&c__1, (char *)&energi_1.ed, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.ned, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___175.ciunit = iounit_1.iout;
	    s_wsfe(&io___175);
	    do_fio(&c__1, (char *)&energi_1.ed, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.ned, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_mpole__ && action_1.nem != 0) {
	if (abs(energi_1.em) < 1e10) {
	    io___176.ciunit = iounit_1.iout;
	    s_wsfe(&io___176);
	    do_fio(&c__1, (char *)&energi_1.em, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nem, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___177.ciunit = iounit_1.iout;
	    s_wsfe(&io___177);
	    do_fio(&c__1, (char *)&energi_1.em, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nem, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_polar__ && action_1.nep != 0) {
	if (abs(energi_1.ep) < 1e10) {
	    io___178.ciunit = iounit_1.iout;
	    s_wsfe(&io___178);
	    do_fio(&c__1, (char *)&energi_1.ep, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nep, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___179.ciunit = iounit_1.iout;
	    s_wsfe(&io___179);
	    do_fio(&c__1, (char *)&energi_1.ep, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nep, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_rxnfld__ && action_1.ner != 0) {
	io___180.ciunit = iounit_1.iout;
	s_wsfe(&io___180);
	do_fio(&c__1, (char *)&energi_1.er, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.ner, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_solv__ && action_1.nes != 0) {
	io___181.ciunit = iounit_1.iout;
	s_wsfe(&io___181);
	do_fio(&c__1, (char *)&energi_1.es, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nes, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_metal__ && action_1.nelf != 0) {
	io___182.ciunit = iounit_1.iout;
	s_wsfe(&io___182);
	do_fio(&c__1, (char *)&energi_1.elf, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nelf, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_geom__ && action_1.neg != 0) {
	io___183.ciunit = iounit_1.iout;
	s_wsfe(&io___183);
	do_fio(&c__1, (char *)&energi_1.eg, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neg, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_extra__ && action_1.nex != 0) {
	io___184.ciunit = iounit_1.iout;
	s_wsfe(&io___184);
	do_fio(&c__1, (char *)&energi_1.ex, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nex, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    return 0;
} /* analyz4_ */



/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine analyz6  --  medium precision energy components  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "analyz6" prints the energy to 6 decimal places and number */
/*     of interactions for each component of the potential energy */


/* Subroutine */ int analyz6_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Energy Component Breakdown :\002,11x,"
	    "\002Kcal/mole\002,6x,\002Interactions\002/)";
    static char fmt_20[] = "(\002 Bond Stretching\002,15x,f18.6,i15)";
    static char fmt_30[] = "(\002 Angle Bending\002,17x,f18.6,i15)";
    static char fmt_40[] = "(\002 Stretch-Bend\002,18x,f18.6,i15)";
    static char fmt_50[] = "(\002 Urey-Bradley\002,18x,f18.6,i15)";
    static char fmt_60[] = "(\002 Angle-Angle\002,19x,f18.6,i15)";
    static char fmt_70[] = "(\002 Out-of-Plane Bend\002,13x,f18.6,i15)";
    static char fmt_80[] = "(\002 Out-of-Plane Distance\002,9x,f18.6,i15)";
    static char fmt_90[] = "(\002 Improper Dihedral\002,13x,f18.6,i15)";
    static char fmt_100[] = "(\002 Improper Torsion\002,14x,f18.6,i15)";
    static char fmt_110[] = "(\002 Torsional Angle\002,15x,f18.6,i15)";
    static char fmt_120[] = "(\002 Pi-Orbital Torsion\002,12x,f18.6,i15)";
    static char fmt_130[] = "(\002 Stretch-Torsion\002,15x,f18.6,i15)";
    static char fmt_140[] = "(\002 Torsion-Torsion\002,15x,f18.6,i15)";
    static char fmt_150[] = "(\002 Van der Waals\002,17x,f18.6,i15)";
    static char fmt_160[] = "(\002 Van der Waals\002,17x,d18.6,i15)";
    static char fmt_170[] = "(\002 Charge-Charge\002,17x,f18.6,i15)";
    static char fmt_180[] = "(\002 Charge-Charge\002,17x,d18.6,i15)";
    static char fmt_190[] = "(\002 Charge-Dipole\002,17x,f18.6,i15)";
    static char fmt_200[] = "(\002 Charge-Dipole\002,17x,d18.6,i15)";
    static char fmt_210[] = "(\002 Dipole-Dipole\002,17x,f18.6,i15)";
    static char fmt_220[] = "(\002 Dipole-Dipole\002,17x,d18.6,i15)";
    static char fmt_230[] = "(\002 Atomic Multipoles\002,13x,f18.6,i15)";
    static char fmt_240[] = "(\002 Atomic Multipoles\002,13x,d18.6,i15)";
    static char fmt_250[] = "(\002 Polarization\002,18x,f18.6,i15)";
    static char fmt_260[] = "(\002 Polarization\002,18x,d18.6,i15)";
    static char fmt_270[] = "(\002 Reaction Field\002,16x,f18.6,i15)";
    static char fmt_280[] = "(\002 Continuum Solvation\002,11x,f18.6,i15)";
    static char fmt_290[] = "(\002 Metal Ligand Field\002,12x,f18.6,i15)";
    static char fmt_300[] = "(\002 Geometric Restraints\002,10x,f18.6,i15)";
    static char fmt_310[] = "(\002 Extra Energy Terms\002,12x,f18.6,i15)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___185 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___186 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___187 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___188 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___189 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___190 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___191 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___192 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___193 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___194 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___195 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___196 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___197 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___198 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___199 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___200 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___201 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___202 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___203 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___204 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___205 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___206 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___207 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___208 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___209 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___210 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___211 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___212 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___213 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___214 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___215 = { 0, 0, 0, fmt_310, 0 };




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  action.i  --  total number of each energy term computed  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     neb     number of bond stretch energy terms computed */
/*     nea     number of angle bend energy terms computed */
/*     neba    number of stretch-bend energy terms computed */
/*     neub    number of Urey-Bradley energy terms computed */
/*     neaa    number of angle-angle energy terms computed */
/*     neopb   number of out-of-plane bend energy terms computed */
/*     neopd   number of out-of-plane distance energy terms computed */
/*     neid    number of improper dihedral energy terms computed */
/*     neit    number of improper torsion energy terms computed */
/*     net     number of torsional energy terms computed */
/*     nept    number of pi-orbital torsion energy terms computed */
/*     nebt    number of stretch-torsion energy terms computed */
/*     nett    number of torsion-torsion energy terms computed */
/*     nev     number of van der Waals energy terms computed */
/*     nec     number of charge-charge energy terms computed */
/*     necd    number of charge-dipole energy terms computed */
/*     ned     number of dipole-dipole energy terms computed */
/*     nem     number of multipole energy terms computed */
/*     nep     number of polarization energy terms computed */
/*     new     number of Ewald summation energy terms computed */
/*     ner     number of reaction field energy terms computed */
/*     nes     number of solvation energy terms computed */
/*     nelf    number of metal ligand field energy terms computed */
/*     neg     number of geometric restraint energy terms computed */
/*     nex     number of extra energy terms computed */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cutoff.i  --  cutoff distances for energy interactions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     vdwcut      cutoff distance for van der Waals interactions */
/*     chgcut      cutoff distance for charge-charge interactions */
/*     dplcut      cutoff distance for dipole-dipole interactions */
/*     mpolecut    cutoff distance for atomic multipole interactions */
/*     vdwtaper    distance at which van der Waals switching begins */
/*     chgtaper    distance at which charge-charge switching begins */
/*     dpltaper    distance at which dipole-dipole switching begins */
/*     mpoletaper  distance at which atomic multipole switching begins */
/*     ewaldcut    cutoff distance for direct space Ewald summation */
/*     use_ewald   logical flag governing use of Ewald summation */
/*     use_lights  logical flag governing use of method of lights */
/*     use_list    logical flag governing use of any neighbor lists */
/*     use_vlist   logical flag governing use of vdw neighbor lists */
/*     use_clist   logical flag governing use of charge neighbor lists */
/*     use_mlist   logical flag governing use of multipole neighbor lists */




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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     write out each energy component to the desired precision */



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  potent.i  --  usage of each potential energy component  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     use_bond    logical flag governing use of bond stretch potential */
/*     use_angle   logical flag governing use of angle bend potential */
/*     use_strbnd  logical flag governing use of stretch-bend potential */
/*     use_urey    logical flag governing use of Urey-Bradley potential */
/*     use_angang  logical flag governing use of angle-angle cross term */
/*     use_opbend  logical flag governing use of out-of-plane bend term */
/*     use_opdist  logical flag governing use of out-of-plane distance */
/*     use_improp  logical flag governing use of improper dihedral term */
/*     use_imptor  logical flag governing use of improper torsion term */
/*     use_tors    logical flag governing use of torsional potential */
/*     use_pitors  logical flag governing use of pi-orbital torsion term */
/*     use_strtor  logical flag governing use of stretch-torsion term */
/*     use_tortor  logical flag governing use of torsion-torsion term */
/*     use_vdw     logical flag governing use of vdw der Waals potential */
/*     use_charge  logical flag governing use of charge-charge potential */
/*     use_chgdpl  logical flag governing use of charge-dipole potential */
/*     use_dipole  logical flag governing use of dipole-dipole potential */
/*     use_mpole   logical flag governing use of multipole potential */
/*     use_polar   logical flag governing use of polarization term */
/*     use_rxnfld  logical flag governing use of reaction field term */
/*     use_solv    logical flag governing use of continuum solvation */
/*     use_metal   logical flag governing use of ligand field term */
/*     use_geom    logical flag governing use of geometric restraints */
/*     use_extra   logical flag governing use of extra potential term */
/*     use_born    logical flag governing use of Born radii values */
/*     use_orbit   logical flag governing use of pisystem computation */


    io___185.ciunit = iounit_1.iout;
    s_wsfe(&io___185);
    e_wsfe();
    if (potent_1.use_bond__ && action_1.neb != 0) {
	io___186.ciunit = iounit_1.iout;
	s_wsfe(&io___186);
	do_fio(&c__1, (char *)&energi_1.eb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neb, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_angle__ && action_1.nea != 0) {
	io___187.ciunit = iounit_1.iout;
	s_wsfe(&io___187);
	do_fio(&c__1, (char *)&energi_1.ea, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nea, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_strbnd__ && action_1.neba != 0) {
	io___188.ciunit = iounit_1.iout;
	s_wsfe(&io___188);
	do_fio(&c__1, (char *)&energi_1.eba, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neba, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_urey__ && action_1.neub != 0) {
	io___189.ciunit = iounit_1.iout;
	s_wsfe(&io___189);
	do_fio(&c__1, (char *)&energi_1.eub, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neub, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_angang__ && action_1.neaa != 0) {
	io___190.ciunit = iounit_1.iout;
	s_wsfe(&io___190);
	do_fio(&c__1, (char *)&energi_1.eaa, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neaa, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_opbend__ && action_1.neopb != 0) {
	io___191.ciunit = iounit_1.iout;
	s_wsfe(&io___191);
	do_fio(&c__1, (char *)&energi_1.eopb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neopb, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_opdist__ && action_1.neopd != 0) {
	io___192.ciunit = iounit_1.iout;
	s_wsfe(&io___192);
	do_fio(&c__1, (char *)&energi_1.eopd, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neopd, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_improp__ && action_1.neid != 0) {
	io___193.ciunit = iounit_1.iout;
	s_wsfe(&io___193);
	do_fio(&c__1, (char *)&energi_1.eid, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neid, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_imptor__ && action_1.neit != 0) {
	io___194.ciunit = iounit_1.iout;
	s_wsfe(&io___194);
	do_fio(&c__1, (char *)&energi_1.eit, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neit, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_tors__ && action_1.net != 0) {
	io___195.ciunit = iounit_1.iout;
	s_wsfe(&io___195);
	do_fio(&c__1, (char *)&energi_1.et, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.net, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_pitors__ && action_1.nept != 0) {
	io___196.ciunit = iounit_1.iout;
	s_wsfe(&io___196);
	do_fio(&c__1, (char *)&energi_1.ept, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nept, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_strtor__ && action_1.nebt != 0) {
	io___197.ciunit = iounit_1.iout;
	s_wsfe(&io___197);
	do_fio(&c__1, (char *)&energi_1.ebt, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nebt, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_tortor__ && action_1.nett != 0) {
	io___198.ciunit = iounit_1.iout;
	s_wsfe(&io___198);
	do_fio(&c__1, (char *)&energi_1.ett, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nett, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_vdw__ && action_1.nev != 0) {
	if (abs(energi_1.ev) < 1e10) {
	    io___199.ciunit = iounit_1.iout;
	    s_wsfe(&io___199);
	    do_fio(&c__1, (char *)&energi_1.ev, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nev, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___200.ciunit = iounit_1.iout;
	    s_wsfe(&io___200);
	    do_fio(&c__1, (char *)&energi_1.ev, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nev, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_charge__ && action_1.nec != 0) {
	if (abs(energi_1.ec) < 1e10) {
	    io___201.ciunit = iounit_1.iout;
	    s_wsfe(&io___201);
	    do_fio(&c__1, (char *)&energi_1.ec, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nec, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___202.ciunit = iounit_1.iout;
	    s_wsfe(&io___202);
	    do_fio(&c__1, (char *)&energi_1.ec, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nec, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_chgdpl__ && action_1.necd != 0) {
	if (abs(energi_1.ecd) < 1e10) {
	    io___203.ciunit = iounit_1.iout;
	    s_wsfe(&io___203);
	    do_fio(&c__1, (char *)&energi_1.ecd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.necd, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___204.ciunit = iounit_1.iout;
	    s_wsfe(&io___204);
	    do_fio(&c__1, (char *)&energi_1.ecd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.necd, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_dipole__ && action_1.ned != 0) {
	if (abs(energi_1.ed) < 1e10) {
	    io___205.ciunit = iounit_1.iout;
	    s_wsfe(&io___205);
	    do_fio(&c__1, (char *)&energi_1.ed, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.ned, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___206.ciunit = iounit_1.iout;
	    s_wsfe(&io___206);
	    do_fio(&c__1, (char *)&energi_1.ed, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.ned, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_mpole__ && action_1.nem != 0) {
	if (abs(energi_1.em) < 1e10) {
	    io___207.ciunit = iounit_1.iout;
	    s_wsfe(&io___207);
	    do_fio(&c__1, (char *)&energi_1.em, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nem, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___208.ciunit = iounit_1.iout;
	    s_wsfe(&io___208);
	    do_fio(&c__1, (char *)&energi_1.em, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nem, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_polar__ && action_1.nep != 0) {
	if (abs(energi_1.ep) < 1e10) {
	    io___209.ciunit = iounit_1.iout;
	    s_wsfe(&io___209);
	    do_fio(&c__1, (char *)&energi_1.ep, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nep, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___210.ciunit = iounit_1.iout;
	    s_wsfe(&io___210);
	    do_fio(&c__1, (char *)&energi_1.ep, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nep, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_rxnfld__ && action_1.ner != 0) {
	io___211.ciunit = iounit_1.iout;
	s_wsfe(&io___211);
	do_fio(&c__1, (char *)&energi_1.er, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.ner, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_solv__ && action_1.nes != 0) {
	io___212.ciunit = iounit_1.iout;
	s_wsfe(&io___212);
	do_fio(&c__1, (char *)&energi_1.es, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nes, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_metal__ && action_1.nelf != 0) {
	io___213.ciunit = iounit_1.iout;
	s_wsfe(&io___213);
	do_fio(&c__1, (char *)&energi_1.elf, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nelf, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_geom__ && action_1.neg != 0) {
	io___214.ciunit = iounit_1.iout;
	s_wsfe(&io___214);
	do_fio(&c__1, (char *)&energi_1.eg, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neg, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_extra__ && action_1.nex != 0) {
	io___215.ciunit = iounit_1.iout;
	s_wsfe(&io___215);
	do_fio(&c__1, (char *)&energi_1.ex, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nex, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    return 0;
} /* analyz6_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine analyz8  --  high precision energy components  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "analyz8" prints the energy to 8 decimal places and number */
/*     of interactions for each component of the potential energy */


/* Subroutine */ int analyz8_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Energy Component Breakdown :\002,11x,"
	    "\002Kcal/mole\002,6x,\002Interactions\002/)";
    static char fmt_20[] = "(\002 Bond Stretching\002,13x,f20.8,i15)";
    static char fmt_30[] = "(\002 Angle Bending\002,15x,f20.8,i15)";
    static char fmt_40[] = "(\002 Stretch-Bend\002,16x,f20.8,i15)";
    static char fmt_50[] = "(\002 Urey-Bradley\002,16x,f20.8,i15)";
    static char fmt_60[] = "(\002 Angle-Angle\002,17x,f20.8,i15)";
    static char fmt_70[] = "(\002 Out-of-Plane Bend\002,11x,f20.8,i15)";
    static char fmt_80[] = "(\002 Out-of-Plane Distance\002,7x,f20.8,i15)";
    static char fmt_90[] = "(\002 Improper Dihedral\002,11x,f20.8,i15)";
    static char fmt_100[] = "(\002 Improper Torsion\002,12x,f20.8,i15)";
    static char fmt_110[] = "(\002 Torsional Angle\002,13x,f20.8,i15)";
    static char fmt_120[] = "(\002 Pi-Orbital Torsion\002,10x,f20.8,i15)";
    static char fmt_130[] = "(\002 Stretch-Torsion\002,13x,f20.8,i15)";
    static char fmt_140[] = "(\002 Torsion-Torsion\002,13x,f20.8,i15)";
    static char fmt_150[] = "(\002 Van der Waals\002,15x,f20.8,i15)";
    static char fmt_160[] = "(\002 Van der Waals\002,15x,d20.8,i15)";
    static char fmt_170[] = "(\002 Charge-Charge\002,15x,f20.8,i15)";
    static char fmt_180[] = "(\002 Charge-Charge\002,15x,d20.8,i15)";
    static char fmt_190[] = "(\002 Charge-Dipole\002,15x,f20.8,i15)";
    static char fmt_200[] = "(\002 Charge-Dipole\002,15x,d20.8,i15)";
    static char fmt_210[] = "(\002 Dipole-Dipole\002,15x,f20.8,i15)";
    static char fmt_220[] = "(\002 Dipole-Dipole\002,15x,d20.8,i15)";
    static char fmt_230[] = "(\002 Atomic Multipoles\002,11x,f20.8,i15)";
    static char fmt_240[] = "(\002 Atomic Multipoles\002,11x,d20.8,i15)";
    static char fmt_250[] = "(\002 Polarization\002,16x,f20.8,i15)";
    static char fmt_260[] = "(\002 Polarization\002,16x,d20.8,i15)";
    static char fmt_270[] = "(\002 Reaction Field\002,14x,f20.8,i15)";
    static char fmt_280[] = "(\002 Continuum Solvation\002,9x,f20.8,i15)";
    static char fmt_290[] = "(\002 Metal Ligand Field\002,10x,f20.8,i15)";
    static char fmt_300[] = "(\002 Geometric Restraints\002,8x,f20.8,i15)";
    static char fmt_310[] = "(\002 Extra Energy Terms\002,10x,f20.8,i15)";

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___216 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___217 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___218 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___219 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___220 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___221 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___222 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___223 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___224 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___225 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___226 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___227 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___228 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___229 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___230 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___231 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___232 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___233 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___234 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___235 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___236 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___237 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___238 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___239 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___240 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___241 = { 0, 0, 0, fmt_260, 0 };
    static cilist io___242 = { 0, 0, 0, fmt_270, 0 };
    static cilist io___243 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___244 = { 0, 0, 0, fmt_290, 0 };
    static cilist io___245 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___246 = { 0, 0, 0, fmt_310, 0 };




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  action.i  --  total number of each energy term computed  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     neb     number of bond stretch energy terms computed */
/*     nea     number of angle bend energy terms computed */
/*     neba    number of stretch-bend energy terms computed */
/*     neub    number of Urey-Bradley energy terms computed */
/*     neaa    number of angle-angle energy terms computed */
/*     neopb   number of out-of-plane bend energy terms computed */
/*     neopd   number of out-of-plane distance energy terms computed */
/*     neid    number of improper dihedral energy terms computed */
/*     neit    number of improper torsion energy terms computed */
/*     net     number of torsional energy terms computed */
/*     nept    number of pi-orbital torsion energy terms computed */
/*     nebt    number of stretch-torsion energy terms computed */
/*     nett    number of torsion-torsion energy terms computed */
/*     nev     number of van der Waals energy terms computed */
/*     nec     number of charge-charge energy terms computed */
/*     necd    number of charge-dipole energy terms computed */
/*     ned     number of dipole-dipole energy terms computed */
/*     nem     number of multipole energy terms computed */
/*     nep     number of polarization energy terms computed */
/*     new     number of Ewald summation energy terms computed */
/*     ner     number of reaction field energy terms computed */
/*     nes     number of solvation energy terms computed */
/*     nelf    number of metal ligand field energy terms computed */
/*     neg     number of geometric restraint energy terms computed */
/*     nex     number of extra energy terms computed */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  cutoff.i  --  cutoff distances for energy interactions  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     vdwcut      cutoff distance for van der Waals interactions */
/*     chgcut      cutoff distance for charge-charge interactions */
/*     dplcut      cutoff distance for dipole-dipole interactions */
/*     mpolecut    cutoff distance for atomic multipole interactions */
/*     vdwtaper    distance at which van der Waals switching begins */
/*     chgtaper    distance at which charge-charge switching begins */
/*     dpltaper    distance at which dipole-dipole switching begins */
/*     mpoletaper  distance at which atomic multipole switching begins */
/*     ewaldcut    cutoff distance for direct space Ewald summation */
/*     use_ewald   logical flag governing use of Ewald summation */
/*     use_lights  logical flag governing use of method of lights */
/*     use_list    logical flag governing use of any neighbor lists */
/*     use_vlist   logical flag governing use of vdw neighbor lists */
/*     use_clist   logical flag governing use of charge neighbor lists */
/*     use_mlist   logical flag governing use of multipole neighbor lists */




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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     write out each energy component to the desired precision */



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  potent.i  --  usage of each potential energy component  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     use_bond    logical flag governing use of bond stretch potential */
/*     use_angle   logical flag governing use of angle bend potential */
/*     use_strbnd  logical flag governing use of stretch-bend potential */
/*     use_urey    logical flag governing use of Urey-Bradley potential */
/*     use_angang  logical flag governing use of angle-angle cross term */
/*     use_opbend  logical flag governing use of out-of-plane bend term */
/*     use_opdist  logical flag governing use of out-of-plane distance */
/*     use_improp  logical flag governing use of improper dihedral term */
/*     use_imptor  logical flag governing use of improper torsion term */
/*     use_tors    logical flag governing use of torsional potential */
/*     use_pitors  logical flag governing use of pi-orbital torsion term */
/*     use_strtor  logical flag governing use of stretch-torsion term */
/*     use_tortor  logical flag governing use of torsion-torsion term */
/*     use_vdw     logical flag governing use of vdw der Waals potential */
/*     use_charge  logical flag governing use of charge-charge potential */
/*     use_chgdpl  logical flag governing use of charge-dipole potential */
/*     use_dipole  logical flag governing use of dipole-dipole potential */
/*     use_mpole   logical flag governing use of multipole potential */
/*     use_polar   logical flag governing use of polarization term */
/*     use_rxnfld  logical flag governing use of reaction field term */
/*     use_solv    logical flag governing use of continuum solvation */
/*     use_metal   logical flag governing use of ligand field term */
/*     use_geom    logical flag governing use of geometric restraints */
/*     use_extra   logical flag governing use of extra potential term */
/*     use_born    logical flag governing use of Born radii values */
/*     use_orbit   logical flag governing use of pisystem computation */


    io___216.ciunit = iounit_1.iout;
    s_wsfe(&io___216);
    e_wsfe();
    if (potent_1.use_bond__ && action_1.neb != 0) {
	io___217.ciunit = iounit_1.iout;
	s_wsfe(&io___217);
	do_fio(&c__1, (char *)&energi_1.eb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neb, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_angle__ && action_1.nea != 0) {
	io___218.ciunit = iounit_1.iout;
	s_wsfe(&io___218);
	do_fio(&c__1, (char *)&energi_1.ea, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nea, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_strbnd__ && action_1.neba != 0) {
	io___219.ciunit = iounit_1.iout;
	s_wsfe(&io___219);
	do_fio(&c__1, (char *)&energi_1.eba, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neba, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_urey__ && action_1.neub != 0) {
	io___220.ciunit = iounit_1.iout;
	s_wsfe(&io___220);
	do_fio(&c__1, (char *)&energi_1.eub, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neub, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_angang__ && action_1.neaa != 0) {
	io___221.ciunit = iounit_1.iout;
	s_wsfe(&io___221);
	do_fio(&c__1, (char *)&energi_1.eaa, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neaa, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_opbend__ && action_1.neopb != 0) {
	io___222.ciunit = iounit_1.iout;
	s_wsfe(&io___222);
	do_fio(&c__1, (char *)&energi_1.eopb, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neopb, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_opdist__ && action_1.neopd != 0) {
	io___223.ciunit = iounit_1.iout;
	s_wsfe(&io___223);
	do_fio(&c__1, (char *)&energi_1.eopd, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neopd, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_improp__ && action_1.neid != 0) {
	io___224.ciunit = iounit_1.iout;
	s_wsfe(&io___224);
	do_fio(&c__1, (char *)&energi_1.eid, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neid, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_imptor__ && action_1.neit != 0) {
	io___225.ciunit = iounit_1.iout;
	s_wsfe(&io___225);
	do_fio(&c__1, (char *)&energi_1.eit, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neit, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_tors__ && action_1.net != 0) {
	io___226.ciunit = iounit_1.iout;
	s_wsfe(&io___226);
	do_fio(&c__1, (char *)&energi_1.et, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.net, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_pitors__ && action_1.nept != 0) {
	io___227.ciunit = iounit_1.iout;
	s_wsfe(&io___227);
	do_fio(&c__1, (char *)&energi_1.ept, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nept, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_strtor__ && action_1.nebt != 0) {
	io___228.ciunit = iounit_1.iout;
	s_wsfe(&io___228);
	do_fio(&c__1, (char *)&energi_1.ebt, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nebt, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_tortor__ && action_1.nett != 0) {
	io___229.ciunit = iounit_1.iout;
	s_wsfe(&io___229);
	do_fio(&c__1, (char *)&energi_1.ett, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nett, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_vdw__ && action_1.nev != 0) {
	if (abs(energi_1.ev) < 1e10) {
	    io___230.ciunit = iounit_1.iout;
	    s_wsfe(&io___230);
	    do_fio(&c__1, (char *)&energi_1.ev, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nev, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___231.ciunit = iounit_1.iout;
	    s_wsfe(&io___231);
	    do_fio(&c__1, (char *)&energi_1.ev, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nev, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_charge__ && action_1.nec != 0) {
	if (abs(energi_1.ec) < 1e10) {
	    io___232.ciunit = iounit_1.iout;
	    s_wsfe(&io___232);
	    do_fio(&c__1, (char *)&energi_1.ec, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nec, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___233.ciunit = iounit_1.iout;
	    s_wsfe(&io___233);
	    do_fio(&c__1, (char *)&energi_1.ec, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nec, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_chgdpl__ && action_1.necd != 0) {
	if (abs(energi_1.ecd) < 1e10) {
	    io___234.ciunit = iounit_1.iout;
	    s_wsfe(&io___234);
	    do_fio(&c__1, (char *)&energi_1.ecd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.necd, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___235.ciunit = iounit_1.iout;
	    s_wsfe(&io___235);
	    do_fio(&c__1, (char *)&energi_1.ecd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.necd, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_dipole__ && action_1.ned != 0) {
	if (abs(energi_1.ed) < 1e10) {
	    io___236.ciunit = iounit_1.iout;
	    s_wsfe(&io___236);
	    do_fio(&c__1, (char *)&energi_1.ed, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.ned, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___237.ciunit = iounit_1.iout;
	    s_wsfe(&io___237);
	    do_fio(&c__1, (char *)&energi_1.ed, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.ned, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_mpole__ && action_1.nem != 0) {
	if (abs(energi_1.em) < 1e10) {
	    io___238.ciunit = iounit_1.iout;
	    s_wsfe(&io___238);
	    do_fio(&c__1, (char *)&energi_1.em, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nem, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___239.ciunit = iounit_1.iout;
	    s_wsfe(&io___239);
	    do_fio(&c__1, (char *)&energi_1.em, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nem, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_polar__ && action_1.nep != 0) {
	if (abs(energi_1.ep) < 1e10) {
	    io___240.ciunit = iounit_1.iout;
	    s_wsfe(&io___240);
	    do_fio(&c__1, (char *)&energi_1.ep, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nep, (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    io___241.ciunit = iounit_1.iout;
	    s_wsfe(&io___241);
	    do_fio(&c__1, (char *)&energi_1.ep, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&action_1.nep, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (potent_1.use_rxnfld__ && action_1.ner != 0) {
	io___242.ciunit = iounit_1.iout;
	s_wsfe(&io___242);
	do_fio(&c__1, (char *)&energi_1.er, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.ner, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_solv__ && action_1.nes != 0) {
	io___243.ciunit = iounit_1.iout;
	s_wsfe(&io___243);
	do_fio(&c__1, (char *)&energi_1.es, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nes, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_metal__ && action_1.nelf != 0) {
	io___244.ciunit = iounit_1.iout;
	s_wsfe(&io___244);
	do_fio(&c__1, (char *)&energi_1.elf, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nelf, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_geom__ && action_1.neg != 0) {
	io___245.ciunit = iounit_1.iout;
	s_wsfe(&io___245);
	do_fio(&c__1, (char *)&energi_1.eg, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.neg, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (potent_1.use_extra__ && action_1.nex != 0) {
	io___246.ciunit = iounit_1.iout;
	s_wsfe(&io___246);
	do_fio(&c__1, (char *)&energi_1.ex, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&action_1.nex, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    return 0;
} /* analyz8_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine propyze  --  electrostatic & inertial analysis  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "propyze" finds and prints the total charge, dipole moment */
/*     components, radius of gyration and moments of inertia */


/* Subroutine */ int propyze_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Total Electric Charge :\002,13x,f12.5"
	    ",\002 Electrons\002)";
    static char fmt_20[] = "(/,\002 Dipole Moment Magnitude :\002,11x,f12.3"
	    ",\002 Debyes\002,//,\002 Dipole X,Y,Z-Components :\002,11x,3f12."
	    "3)";
    static char fmt_30[] = "(/,\002 Quadrupole Moment Tensor :\002,10x,3f12."
	    "3,/,6x,\002(Buckinghams)\002,18x,3f12.3,/,37x,3f12.3)";
    static char fmt_40[] = "(/,\002 Principal Axes Quadrupole :\002,9x,3f12."
	    "3)";
    static char fmt_50[] = "(/,\002 Dielectric Constant :\002,15x,f12.3)";
    static char fmt_60[] = "(\002 Effective Total Charge :\002,12x,f12.5,"
	    "\002 Electrons\002)";
    static char fmt_70[] = "(\002 Effective Dipole Moment :\002,11x,f12.3"
	    ",\002 Debyes\002)";
    static char fmt_80[] = "(/,\002 Radius of Gyration :\002,16x,f12.3,\002 "
	    "Angstroms\002)";
    static char fmt_90[] = "(/,\002 Internal Virial Tensor :\002,12x,3f12.3,"
	    "/,37x,3f12.3,/,37x,3f12.3)";

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int gradient_(doublereal *, doublereal *);
    static integer i__;
    static doublereal rg, energy, derivs[75000]	/* was [3][25000] */;
    extern /* Subroutine */ int gyrate_(doublereal *), inertia_(integer *), 
	    moments_(void);

    /* Fortran I/O blocks */
    static cilist io___247 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___248 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___249 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___250 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___251 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___252 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___253 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___255 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___258 = { 0, 0, 0, fmt_90, 0 };



#define vir_ref(a_1,a_2) virial_1.vir[(a_2)*3 + a_1 - 4]



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
/*     ##  moment.i  --  components of electric multipole moments  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     netchg   net electric charge for the total system */
/*     netdpl   dipole moment magnitude for the total system */
/*     netqdp   diagonal quadrupole (Qxx, Qyy, Qzz) for system */
/*     xdpl     dipole vector x-component in the global frame */
/*     ydpl     dipole vector y-component in the global frame */
/*     zdpl     dipole vector z-component in the global frame */
/*     xxqdp    quadrupole tensor xx-component in global frame */
/*     xyqdp    quadrupole tensor xy-component in global frame */
/*     xzqdp    quadrupole tensor xz-component in global frame */
/*     yxqdp    quadrupole tensor yx-component in global frame */
/*     yyqdp    quadrupole tensor yy-component in global frame */
/*     yzqdp    quadrupole tensor yz-component in global frame */
/*     zxqdp    quadrupole tensor zx-component in global frame */
/*     zyqdp    quadrupole tensor zy-component in global frame */
/*     zzqdp    quadrupole tensor zz-component in global frame */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  virial.i  --  components of internal virial tensor  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     vir    total internal virial Cartesian tensor components */




/*     get the total charge, dipole and quadrupole moments */

    moments_();
    io___247.ciunit = iounit_1.iout;
    s_wsfe(&io___247);
    do_fio(&c__1, (char *)&moment_1.netchg, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___248.ciunit = iounit_1.iout;
    s_wsfe(&io___248);
    do_fio(&c__1, (char *)&moment_1.netdpl, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&moment_1.xdpl, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&moment_1.ydpl, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&moment_1.zdpl, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___249.ciunit = iounit_1.iout;
    s_wsfe(&io___249);
    do_fio(&c__1, (char *)&moment_1.xxqdp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&moment_1.xyqdp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&moment_1.xzqdp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&moment_1.yxqdp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&moment_1.yyqdp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&moment_1.yzqdp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&moment_1.zxqdp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&moment_1.zyqdp, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&moment_1.zzqdp, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___250.ciunit = iounit_1.iout;
    s_wsfe(&io___250);
    do_fio(&c__1, (char *)&moment_1.netqdp[0], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&moment_1.netqdp[1], (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&moment_1.netqdp[2], (ftnlen)sizeof(doublereal));
    e_wsfe();
    if (chgpot_1.dielec != 1.) {
	io___251.ciunit = iounit_1.iout;
	s_wsfe(&io___251);
	do_fio(&c__1, (char *)&chgpot_1.dielec, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___252.ciunit = iounit_1.iout;
	s_wsfe(&io___252);
	d__1 = moment_1.netchg / sqrt(chgpot_1.dielec);
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___253.ciunit = iounit_1.iout;
	s_wsfe(&io___253);
	d__1 = moment_1.netdpl / sqrt(chgpot_1.dielec);
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*     get the radius of gyration and moments of inertia */

    gyrate_(&rg);
    io___255.ciunit = iounit_1.iout;
    s_wsfe(&io___255);
    do_fio(&c__1, (char *)&rg, (ftnlen)sizeof(doublereal));
    e_wsfe();
    inertia_(&c__1);

/*     get the internal virial tensor via gradient calculation */

    gradient_(&energy, derivs);
    io___258.ciunit = iounit_1.iout;
    s_wsfe(&io___258);
    for (i__ = 1; i__ <= 3; ++i__) {
	do_fio(&c__1, (char *)&vir_ref(1, i__), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&vir_ref(2, i__), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&vir_ref(3, i__), (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    return 0;
} /* propyze_ */

#undef vir_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine atomyze  --  individual atom energy analysis  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "atomyze" prints the potential energy components broken */
/*     down by atom and to a choice of precision */


/* Subroutine */ int atomyze_(logical *active)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Potential Energy Breakdown over Atoms "
	    ":\002)";
    static char fmt_20[] = "(/,\002  Atom\002,9x,\002EB\002,14x,\002EA\002,1"
	    "4x,\002EBA\002,13x,\002EUB\002,/,15x,\002EAA\002,13x,\002EOPB"
	    "\002,12x,\002EOPD\002,12x,\002EID\002,/,15x,\002EIT\002,13x,\002"
	    "ET\002,14x,\002EPT\002,13x,\002EBT\002,/,15x,\002ETT\002,13x,"
	    "\002EV\002,14x,\002EC\002,14x,\002ECD\002,/,15x,\002ED\002,14x"
	    ",\002EM\002,14x,\002EP\002,14x,\002ER\002,/,15x,\002ES\002,14x"
	    ",\002ELF\002,13x,\002EG\002,14x,\002EX\002)";
    static char fmt_30[] = "(/,i6,4f16.8,/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8"
	    ",/,6x,4f16.8,/,6x,4f16.8)";
    static char fmt_40[] = "(/,\002  Atom\002,8x,\002EB\002,12x,\002EA\002,1"
	    "2x,\002EBA\002,11x,\002EUB\002,11x,\002EAA\002,/,14x,\002EOPB"
	    "\002,10x,\002EOPD\002,10x,\002EID\002,11x,\002EIT\002,11x,\002E"
	    "T\002,/,14x,\002EPT\002,11x,\002EBT\002,11x,\002ETT\002,11x,\002"
	    "EV\002,12x,\002EC\002,/,14x,\002ECD\002,11x,\002ED\002,12x,\002EM"
	    "\002,12x,\002EP\002,12x,\002ER\002,/,14x,\002ES\002,12x,\002EL"
	    "F\002,11x,\002EG\002,12x,\002EX\002)";
    static char fmt_50[] = "(/,i6,5f14.6,/,6x,5f14.6,/,6x,5f14.6,/,6x,5f14.6"
	    ",/,6x,4f14.6)";
    static char fmt_60[] = "(/,\002  Atom\002,8x,\002EB\002,10x,\002EA\002,1"
	    "0x,\002EBA\002,9x,\002EUB\002,9x,\002EAA\002,9x,\002EOPB\002,/,1"
	    "4x,\002EOPD\002,8x,\002EID\002,9x,\002EIT\002,9x,\002ET\002,10x"
	    ",\002EPT\002,9x,\002EBT\002,/,14x,\002ETT\002,9x,\002EV\002,10x"
	    ",\002EC\002,10x,\002ECD\002,9x,\002ED\002,10x,\002EM\002,/,14x"
	    ",\002EP\002,10x,\002ER\002,10x,\002ES\002,10x,\002ELF\002,9x,"
	    "\002EG\002,10x,\002EX\002)";
    static char fmt_70[] = "(/,i6,6f12.4,/,6x,6f12.4,/,6x,6f12.4,/,6x,6f12.4)"
	    ;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__;

    /* Fortran I/O blocks */
    static cilist io___260 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___261 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___263 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___264 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___265 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___266 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___267 = { 0, 0, 0, fmt_70, 0 };




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
/*     ##  analyz.i  --  energy components partitioned over atoms  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     aesum   total potential energy partitioned over atoms */
/*     aeb     bond stretch energy partitioned over atoms */
/*     aea     angle bend energy partitioned over atoms */
/*     aeba    stretch-bend energy partitioned over atoms */
/*     aeub    Urey-Bradley energy partitioned over atoms */
/*     aeaa    angle-angle energy partitioned over atoms */
/*     aeopb   out-of-plane bend energy partitioned over atoms */
/*     aeopd   out-of-plane distance energy partitioned over atoms */
/*     aeid    improper dihedral energy partitioned over atoms */
/*     aeit    improper torsion energy partitioned over atoms */
/*     aet     torsional energy partitioned over atoms */
/*     aept    pi-orbital torsion energy partitioned over atoms */
/*     aebt    stretch-torsion energy partitioned over atoms */
/*     aett    torsion-torsion energy partitioned over atoms */
/*     aev     van der Waals energy partitioned over atoms */
/*     aec     charge-charge energy partitioned over atoms */
/*     aecd    charge-dipole energy partitioned over atoms */
/*     aed     dipole-dipole energy partitioned over atoms */
/*     aem     multipole energy partitioned over atoms */
/*     aep     polarization energy partitioned over atoms */
/*     aer     reaction field energy partitioned over atoms */
/*     aes     solvation energy partitioned over atoms */
/*     aelf    metal ligand field energy partitioned over atoms */
/*     aeg     geometric restraint energy partitioned over atoms */
/*     aex     extra energy term partitioned over atoms */




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




/*     energy partitioning over the individual atoms */

    /* Parameter adjustments */
    --active;

    /* Function Body */
    io___260.ciunit = iounit_1.iout;
    s_wsfe(&io___260);
    e_wsfe();
    if (inform_1.digits >= 8) {
	io___261.ciunit = iounit_1.iout;
	s_wsfe(&io___261);
	e_wsfe();
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (active[i__]) {
		io___263.ciunit = iounit_1.iout;
		s_wsfe(&io___263);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&analyz_1.aeb[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aea[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeba[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeub[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeaa[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeopb[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeopd[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeid[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeit[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aet[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aept[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aebt[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aett[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aev[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aec[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aecd[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aed[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aem[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aep[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aer[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aes[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aelf[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeg[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aex[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    } else if (inform_1.digits >= 6) {
	io___264.ciunit = iounit_1.iout;
	s_wsfe(&io___264);
	e_wsfe();
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (active[i__]) {
		io___265.ciunit = iounit_1.iout;
		s_wsfe(&io___265);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&analyz_1.aeb[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aea[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeba[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeub[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeaa[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeopb[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeopd[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeid[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeit[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aet[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aept[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aebt[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aett[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aev[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aec[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aecd[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aed[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aem[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aep[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aer[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aes[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aelf[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeg[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aex[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    } else {
	io___266.ciunit = iounit_1.iout;
	s_wsfe(&io___266);
	e_wsfe();
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (active[i__]) {
		io___267.ciunit = iounit_1.iout;
		s_wsfe(&io___267);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&analyz_1.aeb[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aea[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeba[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeub[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeaa[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeopb[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeopd[i__ - 1], (ftnlen)
			sizeof(doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeid[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeit[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aet[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aept[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aebt[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aett[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aev[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aec[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aecd[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aed[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aem[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aep[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aer[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aes[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aelf[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aeg[i__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&analyz_1.aex[i__ - 1], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
	    }
	}
    }
    return 0;
} /* atomyze_ */

/* Main program alias */ int analyze_ () { MAIN__ (); return 0; }
