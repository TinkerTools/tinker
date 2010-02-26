/* initprm.f -- translated by f2c (version 20050501).
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
    doublereal rfsize, rfbulkd;
    integer rfterms;
} rxnpot_;

#define rxnpot_1 rxnpot_

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
    doublereal cury, qury, ureyunit;
} urypot_;

#define urypot_1 urypot_

struct {
    doublereal idihunit, itorunit, torsunit, ptorunit, storunit, ttorunit;
} torpot_;

#define torpot_1 torpot_

struct {
    doublereal abuck, bbuck, cbuck, ghal, dhal, v2scale, v3scale, v4scale, 
	    v5scale, igauss[20]	/* was [2][10] */;
    integer ngauss;
    char vdwindex[5], vdwtyp[13], radtyp[5], radsiz[8], radrule[10], epsrule[
	    10], gausstyp[8];
} vdwpot_;

#define vdwpot_1 vdwpot_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine initprm  --  initialize force field parameters  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "initprm" completely initializes a force field by setting all */
/*     parameters to zero and using defaults for control values */


/* Subroutine */ int initprm_(void)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j;
    static char blank3[3], blank8[8], blank20[20], blank12[12], blank24[24], 
	    blank16[16];


#define describe_ref(a_0,a_1) &katoms_1.describe[(a_1)*24 + a_0 - 24]
#define ka_ref(a_0,a_1) &kangs_1.ka[(a_1)*12 + a_0 - 12]
#define kb_ref(a_0,a_1) &kbonds_1.kb[(a_1)*8 + a_0 - 8]
#define kd_ref(a_0,a_1) &kdipol_1.kd[(a_1)*8 + a_0 - 8]
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
#define kt4_ref(a_0,a_1) &ktorsn_1.kt4[(a_1)*16 + a_0 - 16]
#define kt5_ref(a_0,a_1) &ktorsn_1.kt5[(a_1)*16 + a_0 - 16]
#define kaf_ref(a_0,a_1) &kangs_1.kaf[(a_1)*12 + a_0 - 12]
#define khb_ref(a_0,a_1) &khbond_1.khb[(a_1)*8 + a_0 - 8]
#define kdi_ref(a_0,a_1) &kiprop_1.kdi[(a_1)*16 + a_0 - 16]
#define kel_ref(a_0,a_1) &kbonds_1.kel[(a_1)*12 + a_0 - 12]
#define ksb_ref(a_0,a_1) &kstbnd_1.ksb[(a_1)*12 + a_0 - 12]
#define kbt_ref(a_0,a_1) &ksttor_1.kbt[(a_1)*16 + a_0 - 16]
#define kpi_ref(a_0,a_1) &korbs_1.kpi[(a_1)*8 + a_0 - 8]
#define kti_ref(a_0,a_1) &kitors_1.kti[(a_1)*16 + a_0 - 16]
#define kmp_ref(a_0,a_1) &kmulti_1.kmp[(a_1)*16 + a_0 - 16]
#define kpt_ref(a_0,a_1) &kpitor_1.kpt[(a_1)*8 + a_0 - 8]
#define ktt_ref(a_0,a_1) &ktrtor_1.ktt[(a_1)*20 + a_0 - 20]
#define kpi4_ref(a_0,a_1) &korbs_1.kpi4[(a_1)*8 + a_0 - 8]
#define kpi5_ref(a_0,a_1) &korbs_1.kpi5[(a_1)*8 + a_0 - 8]
#define anan_ref(a_1,a_2) kanang_1.anan[(a_2)*3 + a_1 - 4]
#define kopb_ref(a_0,a_1) &kopbnd_1.kopb[(a_1)*16 + a_0 - 16]
#define kopd_ref(a_0,a_1) &kopdst_1.kopd[(a_1)*16 + a_0 - 16]
#define stbn_ref(a_1,a_2) kstbnd_1.stbn[(a_2)*2 + a_1 - 3]
#define pgrp_ref(a_1,a_2) kpolr_1.pgrp[(a_2)*8 + a_1 - 9]
#define kvpr_ref(a_0,a_1) &kvdwpr_1.kvpr[(a_1)*8 + a_0 - 8]
#define symbol_ref(a_0,a_1) &katoms_1.symbol[(a_1)*3 + a_0 - 3]



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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  rxnpot.i  --  specifics of reaction field functional form  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     rfsize    radius of reaction field sphere centered at origin */
/*     rfbulkd   bulk dielectric constant of reaction field continuum */
/*     rfterms   number of terms to use in reaction field summation */




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

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  torpot.i  --  specifics of torsional functional forms  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     idihunit  convert improper dihedral energy to kcal/mole */
/*     itorunit  convert improper torsion amplitudes to kcal/mole */
/*     torsunit  convert torsional parameter amplitudes to kcal/mole */
/*     ptorunit  convert pi-orbital torsion energy to kcal/mole */
/*     storunit  convert stretch-torsion energy to kcal/mole */
/*     ttorunit  convert stretch-torsion energy to kcal/mole */




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




/*     define blank character strings of various lengths */

    s_copy(blank3, "   ", (ftnlen)3, (ftnlen)3);
    s_copy(blank8, "        ", (ftnlen)8, (ftnlen)8);
    s_copy(blank12, "            ", (ftnlen)12, (ftnlen)12);
    s_copy(blank16, "                ", (ftnlen)16, (ftnlen)16);
    s_copy(blank20, "                    ", (ftnlen)20, (ftnlen)20);
    s_copy(blank24, "                        ", (ftnlen)24, (ftnlen)24);

/*     initialize strings of parameter atom types and classes */

    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kvpr_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(khb_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 2000; ++i__) {
	s_copy(kb_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kb5_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kb4_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kb3_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kel_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12);
    }
    for (i__ = 1; i__ <= 2000; ++i__) {
	s_copy(ka_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(ka5_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(ka4_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(ka3_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kaf_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12);
    }
    for (i__ = 1; i__ <= 2000; ++i__) {
	s_copy(ksb_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12);
    }
    for (i__ = 1; i__ <= 2000; ++i__) {
	s_copy(ku_ref(0, i__), blank12, (ftnlen)12, (ftnlen)12);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kopb_ref(0, i__), blank8, (ftnlen)16, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kopd_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kdi_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kti_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16);
    }
    for (i__ = 1; i__ <= 2000; ++i__) {
	s_copy(kt_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kt5_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kt4_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kpt_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kbt_ref(0, i__), blank16, (ftnlen)16, (ftnlen)16);
    }
    for (i__ = 1; i__ <= 100; ++i__) {
	s_copy(ktt_ref(0, i__), blank20, (ftnlen)20, (ftnlen)20);
    }
    for (i__ = 1; i__ <= 1000; ++i__) {
	s_copy(kd_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kd5_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kd4_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kd3_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 2000; ++i__) {
	s_copy(kmp_ref(0, i__), blank12, (ftnlen)16, (ftnlen)12);
    }
    for (i__ = 1; i__ <= 500; ++i__) {
	s_copy(kpi_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 200; ++i__) {
	s_copy(kpi5_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }
    for (i__ = 1; i__ <= 200; ++i__) {
	s_copy(kpi4_ref(0, i__), blank8, (ftnlen)8, (ftnlen)8);
    }

/*     initialize some of the force field parameters */

    s_copy(fields_1.forcefield, blank20, (ftnlen)20, (ftnlen)20);
    for (i__ = 1; i__ <= 5000; ++i__) {
	s_copy(symbol_ref(0, i__), blank3, (ftnlen)3, (ftnlen)3);
	katoms_1.atmcls[i__ - 1] = 0;
	katoms_1.atmnum[i__ - 1] = 0;
	katoms_1.weight[i__ - 1] = 0.;
	katoms_1.ligand[i__ - 1] = 0;
	s_copy(describe_ref(0, i__), blank24, (ftnlen)24, (ftnlen)24);
	kvdws_1.rad[i__ - 1] = 0.;
	kvdws_1.eps[i__ - 1] = 0.;
	kvdws_1.rad4[i__ - 1] = 0.;
	kvdws_1.eps4[i__ - 1] = 0.;
	kvdws_1.reduct[i__ - 1] = 0.;
	kchrge_1.chg[i__ - 1] = 0.;
	kpolr_1.polr[i__ - 1] = 0.;
	kpolr_1.athl[i__ - 1] = 0.;
	for (j = 1; j <= 8; ++j) {
	    pgrp_ref(j, i__) = 0;
	}
    }
    for (i__ = 1; i__ <= 1000; ++i__) {
	for (j = 1; j <= 2; ++j) {
	    stbn_ref(j, i__) = 0.;
	}
	for (j = 1; j <= 3; ++j) {
	    anan_ref(j, i__) = 0.;
	}
	korbs_1.electron[i__ - 1] = 0.;
	korbs_1.ionize[i__ - 1] = 0.;
	korbs_1.repulse[i__ - 1] = 0.;
    }
    for (i__ = 1; i__ <= 10000; ++i__) {
	fields_1.biotyp[i__ - 1] = 0;
    }

/*     set default control parameters for local geometry terms */

    s_copy(bndpot_1.bndtyp, "HARMONIC", (ftnlen)8, (ftnlen)8);
    bndpot_1.bndunit = 1.;
    bndpot_1.cbnd = 0.;
    bndpot_1.qbnd = 0.;
    angpot_1.angunit = 3.0461741978670857e-4;
    angpot_1.cang = 0.;
    angpot_1.qang = 0.;
    angpot_1.pang = 0.;
    angpot_1.sang = 0.;
    angpot_1.stbnunit = .017453292519943295;
    urypot_1.ureyunit = 1.;
    urypot_1.cury = 0.;
    urypot_1.qury = 0.;
    angpot_1.aaunit = 3.0461741978670857e-4;
    s_copy(angpot_1.opbtyp, "W-D-C", (ftnlen)8, (ftnlen)5);
    angpot_1.opbunit = 3.0461741978670857e-4;
    angpot_1.copb = 0.;
    angpot_1.qopb = 0.;
    angpot_1.popb = 0.;
    angpot_1.sopb = 0.;
    angpot_1.opdunit = 1.;
    angpot_1.copd = 0.;
    angpot_1.qopd = 0.;
    angpot_1.popd = 0.;
    angpot_1.sopd = 0.;
    torpot_1.idihunit = 1.;
    torpot_1.itorunit = 1.;
    torpot_1.torsunit = 1.;
    torpot_1.ptorunit = 1.;
    torpot_1.storunit = 1.;
    torpot_1.ttorunit = 1.;

/*     set default control parameters for van der Waals terms */

    s_copy(vdwpot_1.vdwindex, "CLASS", (ftnlen)5, (ftnlen)5);
    s_copy(vdwpot_1.vdwtyp, "LENNARD-JONES", (ftnlen)13, (ftnlen)13);
    s_copy(vdwpot_1.radrule, "ARITHMETIC", (ftnlen)10, (ftnlen)10);
    s_copy(vdwpot_1.radtyp, "R-MIN", (ftnlen)5, (ftnlen)5);
    s_copy(vdwpot_1.radsiz, "RADIUS", (ftnlen)8, (ftnlen)6);
    s_copy(vdwpot_1.epsrule, "GEOMETRIC", (ftnlen)10, (ftnlen)9);
    s_copy(vdwpot_1.gausstyp, "NONE", (ftnlen)8, (ftnlen)4);
    vdwpot_1.ngauss = 0;
    vdwpot_1.abuck = 0.;
    vdwpot_1.bbuck = 0.;
    vdwpot_1.cbuck = 0.;
    vdwpot_1.ghal = .12;
    vdwpot_1.dhal = .07;
    vdwpot_1.v2scale = 0.;
    vdwpot_1.v3scale = 0.;
    vdwpot_1.v4scale = 1.;
    vdwpot_1.v5scale = 1.;

/*     set default control parameters for charge-charge terms */

    chgpot_1.electric = 332.063709;
    chgpot_1.dielec = 1.;
    chgpot_1.ebuffer = 0.;
    chgpot_1.c2scale = 0.;
    chgpot_1.c3scale = 0.;
    chgpot_1.c4scale = 1.;
    chgpot_1.c5scale = 1.;
    chgpot_1.neutnbr = FALSE_;
    chgpot_1.neutcut = FALSE_;

/*     set default control parameters for polarizable multipoles */

    mplpot_1.m2scale = 0.;
    mplpot_1.m3scale = 0.;
    mplpot_1.m4scale = 1.;
    mplpot_1.m5scale = 1.;
    polpot_1.p2scale = 0.;
    polpot_1.p3scale = 0.;
    polpot_1.p4scale = 1.;
    polpot_1.p5scale = 1.;

/*     set default control parameters for induced dipoles */

    s_copy(polpot_1.poltyp, "MUTUAL", (ftnlen)6, (ftnlen)6);
    polpot_1.poleps = 1e-6;
    polpot_1.polsor = .7;
    polpot_1.d1scale = 0.;
    polpot_1.d2scale = 1.;
    polpot_1.d3scale = 1.;
    polpot_1.d4scale = 1.;
    polpot_1.u1scale = 1.;
    polpot_1.u2scale = 1.;
    polpot_1.u3scale = 1.;
    polpot_1.u4scale = 1.;

/*     set default control parameters for continuum solvation */

    s_copy(solute_1.solvtyp, blank8, (ftnlen)8, (ftnlen)8);
    s_copy(solute_1.borntyp, blank8, (ftnlen)8, (ftnlen)8);

/*     set default control parameters for reaction field */

    rxnpot_1.rfsize = 1e6;
    rxnpot_1.rfbulkd = 80.;
    rxnpot_1.rfterms = 1;
    return 0;
} /* initprm_ */

#undef symbol_ref
#undef kvpr_ref
#undef pgrp_ref
#undef stbn_ref
#undef kopd_ref
#undef kopb_ref
#undef anan_ref
#undef kpi5_ref
#undef kpi4_ref
#undef ktt_ref
#undef kpt_ref
#undef kmp_ref
#undef kti_ref
#undef kpi_ref
#undef kbt_ref
#undef ksb_ref
#undef kel_ref
#undef kdi_ref
#undef khb_ref
#undef kaf_ref
#undef kt5_ref
#undef kt4_ref
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
#undef kd_ref
#undef kb_ref
#undef ka_ref
#undef describe_ref


