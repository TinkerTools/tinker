/* protein.f -- translated by f2c (version 20050501).
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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer nseq, nchain, ichain[20000]	/* was [2][10000] */, seqtyp[10000];
    char seq[30000], chnnam[10000];
} sequen_;

#define sequen_1 sequen_

struct {
    integer ltitle;
    char title[120];
} titles_;

#define titles_1 titles_

struct {
    doublereal phi[10000], psi[10000], omega[10000], chi[40000]	/* was [4][
	    10000] */;
    integer chiral[10000], disulf[10000];
} phipsi_;

#define phipsi_1 phipsi_

struct {
    char amino[93], nuclz[36], amino1[31], nuclz1[12];
} resdue_;

#define resdue_1 resdue_

struct {
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

struct {
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

struct {
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    doublereal weight[5000];
    integer atmcls[5000], atmnum[5000], ligand[5000];
    char symbol[15000], describe[120000];
} katoms_;

#define katoms_1 katoms_

struct {
    doublereal xpfix[25000], ypfix[25000], zpfix[25000], pfix[50000]	/* 
	    was [2][25000] */, dfix[75000]	/* was [3][25000] */, afix[
	    75000]	/* was [3][25000] */, tfix[75000]	/* was [3][
	    25000] */, gfix[75000]	/* was [3][25000] */, chir[75000]	
	    /* was [3][25000] */, depth, width, rwall;
    integer npfix, ipfix[25000], kpfix[75000]	/* was [3][25000] */, ndfix, 
	    idfix[50000]	/* was [2][25000] */, nafix, iafix[75000]	
	    /* was [3][25000] */, ntfix, itfix[100000]	/* was [4][25000] */, 
	    ngfix, igfix[50000]	/* was [2][25000] */, nchir, ichir[100000]	
	    /* was [4][25000] */;
    logical use_basin__, use_wall__;
} kgeoms_;

#define kgeoms_1 kgeoms_

struct {
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

struct {
    logical archive, noversion, overwrite, cyclesave;
    char coordtype[9];
} output_;

#define output_1 output_

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
    doublereal xrb[25000], yrb[25000], zrb[25000], rbc[6000]	/* was [6][
	    1000] */;
    logical use_rigid__;
} rigid_;

#define rigid_1 rigid_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__5 = 5;
static integer c__3 = 3;
static doublereal c_b84 = 0.;
static integer c__0 = 0;
static doublereal c_b91 = 1.46;
static doublereal c_b97 = 1.51;
static doublereal c_b98 = 110.7;
static doublereal c_b102 = 30.;
static doublereal c_b103 = 150.;
static doublereal c_b104 = 180.;
static integer c_n2 = -2;
static doublereal c_b120 = 1.22;
static doublereal c_b121 = 122.5;
static doublereal c_b123 = 1.02;
static doublereal c_b124 = 121.;
static doublereal c_b126 = 1.11;
static doublereal c_b127 = 109.5;
static doublereal c_b128 = 107.9;
static doublereal c_b143 = 1.12;
static doublereal c_b144 = 120.;
static doublereal c_b209 = 109.4;
static integer c_n1 = -1;
static doublereal c_b223 = 1.5;
static doublereal c_b230 = 111.6;
static doublereal c_b261 = 1.25;
static doublereal c_b262 = 117.;
static doublereal c_b273 = -120.;
static doublereal c_b332 = -60.;
static doublereal c_b393 = 108.;
static doublereal c_b402 = 1.34;
static doublereal c_b403 = 112.7;
static doublereal c_b449 = 119.;
static doublereal c_b465 = 118.;
static integer c__355 = 355;
static integer c__506 = 506;
static integer c__6 = 6;
static integer c__13 = 13;
static doublereal c_b517 = 1.54;
static doublereal c_b519 = 107.8;
static integer c__14 = 14;
static integer c__21 = 21;
static integer c__23 = 23;
static integer c__25 = 25;
static integer c__22 = 22;
static integer c__24 = 24;
static integer c__26 = 26;
static integer c__33 = 33;
static integer c__35 = 35;
static integer c__37 = 37;
static integer c__39 = 39;
static integer c__34 = 34;
static integer c__36 = 36;
static integer c__38 = 38;
static integer c__40 = 40;
static integer c__47 = 47;
static integer c__49 = 49;
static integer c__51 = 51;
static integer c__53 = 53;
static integer c__48 = 48;
static integer c__50 = 50;
static integer c__52 = 52;
static doublereal c_b681 = 110.;
static doublereal c_b687 = 109.;
static integer c__54 = 54;
static integer c__61 = 61;
static integer c__63 = 63;
static doublereal c_b715 = 1.41;
static doublereal c_b716 = 107.5;
static integer c__62 = 62;
static doublereal c_b721 = 106.7;
static integer c__64 = 64;
static doublereal c_b729 = .94;
static doublereal c_b730 = 106.9;
static integer c__71 = 71;
static integer c__73 = 73;
static integer c__75 = 75;
static doublereal c_b744 = 107.7;
static integer c__72 = 72;
static integer c__74 = 74;
static integer c__76 = 76;
static integer c__83 = 83;
static integer c__85 = 85;
static doublereal c_b776 = 1.82;
static integer c__84 = 84;
static doublereal c_b782 = 112.;
static integer c__86 = 86;
static doublereal c_b791 = 96.;
static integer c__93 = 93;
static integer c__95 = 95;
static integer c__94 = 94;
static integer c__101 = 101;
static doublereal c_b839 = 107.;
static integer c__103 = 103;
static integer c__410 = 410;
static integer c__105 = 105;
static integer c__102 = 102;
static integer c__104 = 104;
static integer c__411 = 411;
static integer c__106 = 106;
static integer c__113 = 113;
static integer c__115 = 115;
static integer c__116 = 116;
static doublereal c_b909 = 1.39;
static integer c__118 = 118;
static integer c__120 = 120;
static integer c__114 = 114;
static integer c__117 = 117;
static doublereal c_b949 = 1.1;
static integer c__119 = 119;
static integer c__121 = 121;
static integer c__128 = 128;
static integer c__130 = 130;
static integer c__131 = 131;
static integer c__133 = 133;
static integer c__135 = 135;
static integer c__136 = 136;
static doublereal c_b1013 = 1.36;
static integer c__129 = 129;
static integer c__132 = 132;
static integer c__134 = 134;
static integer c__137 = 137;
static doublereal c_b1048 = .97;
static integer c__144 = 144;
static integer c__146 = 146;
static integer c__147 = 147;
static doublereal c_b1062 = 1.35;
static doublereal c_b1063 = 126.;
static integer c__149 = 149;
static integer c__150 = 150;
static integer c__152 = 152;
static integer c__153 = 153;
static integer c__155 = 155;
static integer c__157 = 157;
static integer c__159 = 159;
static integer c__145 = 145;
static integer c__148 = 148;
static integer c__151 = 151;
static doublereal c_b1128 = 1.05;
static integer c__154 = 154;
static integer c__156 = 156;
static integer c__158 = 158;
static integer c__160 = 160;
static integer c__167 = 167;
static integer c__169 = 169;
static integer c__170 = 170;
static integer c__172 = 172;
static integer c__174 = 174;
static integer c__176 = 176;
static integer c__168 = 168;
static integer c__171 = 171;
static integer c__173 = 173;
static integer c__175 = 175;
static integer c__177 = 177;
static integer c__184 = 184;
static integer c__186 = 186;
static integer c__187 = 187;
static integer c__189 = 189;
static integer c__191 = 191;
static integer c__193 = 193;
static integer c__185 = 185;
static integer c__188 = 188;
static integer c__190 = 190;
static integer c__192 = 192;
static integer c__200 = 200;
static integer c__202 = 202;
static integer c__203 = 203;
static integer c__204 = 204;
static integer c__206 = 206;
static integer c__208 = 208;
static integer c__201 = 201;
static integer c__205 = 205;
static integer c__207 = 207;
static integer c__209 = 209;
static integer c__216 = 216;
static integer c__218 = 218;
static integer c__219 = 219;
static integer c__217 = 217;
static integer c__226 = 226;
static integer c__228 = 228;
static integer c__229 = 229;
static integer c__230 = 230;
static doublereal c_b1378 = 124.;
static integer c__227 = 227;
static integer c__231 = 231;
static integer c__238 = 238;
static integer c__240 = 240;
static integer c__242 = 242;
static integer c__243 = 243;
static integer c__239 = 239;
static integer c__241 = 241;
static integer c__250 = 250;
static integer c__252 = 252;
static integer c__254 = 254;
static integer c__255 = 255;
static integer c__256 = 256;
static integer c__251 = 251;
static integer c__253 = 253;
static integer c__257 = 257;
static integer c__264 = 264;
static integer c__266 = 266;
static integer c__268 = 268;
static integer c__269 = 269;
static doublereal c_b1509 = 96.3;
static integer c__265 = 265;
static integer c__267 = 267;
static integer c__270 = 270;
static integer c__277 = 277;
static integer c__279 = 279;
static integer c__281 = 281;
static integer c__283 = 283;
static integer c__285 = 285;
static integer c__278 = 278;
static integer c__280 = 280;
static integer c__282 = 282;
static integer c__284 = 284;
static doublereal c_b1600 = 108.8;
static integer c__286 = 286;
static integer c__293 = 293;
static integer c__295 = 295;
static integer c__297 = 297;
static integer c__299 = 299;
static doublereal c_b1636 = 1.45;
static integer c__301 = 301;
static integer c__302 = 302;
static integer c__294 = 294;
static integer c__296 = 296;
static integer c__298 = 298;
static integer c__300 = 300;
static integer c__303 = 303;
static integer c__310 = 310;
static integer c__312 = 312;
static integer c__314 = 314;
static integer c__316 = 316;
static integer c__311 = 311;
static integer c__313 = 313;
static integer c__315 = 315;
static integer c__317 = 317;
static integer c__323 = 323;
static integer c__324 = 324;
static integer c__331 = 331;
static integer c__333 = 333;
static integer c__335 = 335;
static integer c__336 = 336;
static integer c__332 = 332;
static integer c__334 = 334;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  program protein  --  build a polypeptide from sequence  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "protein" builds the internal and Cartesian coordinates */
/*     of a polypeptide from amino acid sequence and torsional */
/*     angle values for the peptide backbone and side chains */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter Name to be used for Output Files :"
	    "  \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(/,\002 Enter Title :  \002,$)";
    static char fmt_40[] = "(a120)";

    /* System generated locals */
    address a__1[2];
    integer i__1[2];
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int basefile_(char *, ftnlen), prochain_(void), 
	    molecule_(void);
    extern integer freeunit_(void), trimtext_(char *, ftnlen);
    static integer i__, mode, iseq, izmt, ixyz;
    extern /* Subroutine */ int field_(void), final_(void);
    static logical clash;
    static integer natom;
    static logical exist;
    extern /* Subroutine */ int delete_(integer *), attach_(void), getkey_(
	    void), getseq_(void), prtseq_(integer *), chkxyz_(logical *), 
	    prtint_(integer *), prtxyz_(integer *);
    static char seqfile[120];
    extern /* Subroutine */ int initial_(void);
    static char intfile[120];
    extern /* Subroutine */ int connect_(void), inertia_(integer *), makeint_(
	    integer *), pauling_(void), nextarg_(char *, logical *, ftnlen), 
	    version_(char *, char *, ftnlen, ftnlen), makexyz_(void);
    static char xyzfile[120];

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___3 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_40, 0 };




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
/*     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     iout    Fortran I/O unit for main output (default=6) */
/*     input   Fortran I/O unit for main input (default=5) */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  sequen.i  --  sequence information for a biopolymer  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nseq     total number of residues in biopolymer sequences */
/*     nchain   number of separate biopolymer sequence chains */
/*     ichain   first and last residue in each biopolymer chain */
/*     seqtyp   residue type for each residue in the sequence */
/*     seq      three-letter code for each residue in the sequence */
/*     chnnam   one-letter identifier for each sequence chain */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  titles.i  --  title for the current molecular system  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     ltitle   length in characters of the nonblank title string */
/*     title    title used to describe the current structure */




/*     get the name to use for the output structure files */

    initial_();
    nextarg_(files_1.filename, &exist, (ftnlen)120);
    if (! exist) {
	io___2.ciunit = iounit_1.iout;
	s_wsfe(&io___2);
	e_wsfe();
	io___3.ciunit = iounit_1.input;
	s_rsfe(&io___3);
	do_fio(&c__1, files_1.filename, (ftnlen)120);
	e_rsfe();
    }
    basefile_(files_1.filename, (ftnlen)120);

/*     get the title line for the output files */

    io___4.ciunit = iounit_1.iout;
    s_wsfe(&io___4);
    e_wsfe();
    io___5.ciunit = iounit_1.input;
    s_rsfe(&io___5);
    do_fio(&c__1, titles_1.title, (ftnlen)120);
    e_rsfe();
    titles_1.ltitle = trimtext_(titles_1.title, (ftnlen)120);

/*     read the keyfile and force field parameter file */

    getkey_();
    field_();

/*     get the sequence, build a Z-matrix, convert to Cartesians */

    getseq_();
    prochain_();
    connect_();
    attach_();
    molecule_();
    makexyz_();

/*     perform a packing calculation for multiple chains */

    if (sequen_1.nchain > 1) {
	pauling_();
	inertia_(&c__2);
    }

/*     remove any dummy atoms from Cartesian coordinates */

    natom = atoms_1.n;
    for (i__ = natom; i__ >= 1; --i__) {
	if (atoms_1.type__[i__ - 1] == 0) {
	    delete_(&i__);
	}
    }

/*     convert to internal and Cartesian coordinates */

    mode = 0;
    makeint_(&mode);
    makexyz_();

/*     check for atom pairs with identical coordinates */

    clash = FALSE_;
    chkxyz_(&clash);

/*     write out a amino acid sequence file */

    iseq = freeunit_();
/* Writing concatenation */
    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
    i__1[1] = 4, a__1[1] = ".seq";
    s_cat(seqfile, a__1, i__1, &c__2, (ftnlen)120);
    version_(seqfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = iseq;
    o__1.ofnmlen = 120;
    o__1.ofnm = seqfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    prtseq_(&iseq);
    cl__1.cerr = 0;
    cl__1.cunit = iseq;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     write out an internal coordinates file */

    izmt = freeunit_();
/* Writing concatenation */
    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
    i__1[1] = 4, a__1[1] = ".int";
    s_cat(intfile, a__1, i__1, &c__2, (ftnlen)120);
    version_(intfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = izmt;
    o__1.ofnmlen = 120;
    o__1.ofnm = intfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    prtint_(&izmt);
    cl__1.cerr = 0;
    cl__1.cunit = izmt;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     write out a Cartesian coordinates file */

    ixyz = freeunit_();
/* Writing concatenation */
    i__1[0] = files_1.leng, a__1[0] = files_1.filename;
    i__1[1] = 4, a__1[1] = ".xyz";
    s_cat(xyzfile, a__1, i__1, &c__2, (ftnlen)120);
    version_(xyzfile, "new", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = ixyz;
    o__1.ofnmlen = 120;
    o__1.ofnm = xyzfile;
    o__1.orl = 0;
    o__1.osta = "new";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    prtxyz_(&ixyz);
    cl__1.cerr = 0;
    cl__1.cunit = ixyz;
    cl__1.csta = 0;
    f_clos(&cl__1);

/*     perform any final tasks before program exit */

    final_();
    return 0;
} /* MAIN__ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine getseq  --  amino acid sequence and angles  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "getseq" asks the user for the amino acid sequence */
/*     and torsional angle values needed to define a peptide */


/* Subroutine */ int getseq_(void)
{
    /* Initialized data */

    static char ucase[1*26] = "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" 
	    "L" "M" "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X" "Y" "Z";

    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter One Residue per Line, Free Forma"
	    "t: \002,\002 1 or 3 Letter Code, then\002,/,\002 Phi/Psi/Omega ("
	    "3F), Chi Angles (4F), then\002,\002 Disulfide Partner if a\002,/,"
	    "\002 CYS (I), and D/L Chirality as desired (A1)\002,//,\002 The "
	    "allowed N-Cap Residues are Acetyl=ACE or\002,\002 Formyl=FOR, an"
	    "d the\002,/,\002 possible C-Cap Residues are N-MethylAmide=NM"
	    "E\002,\002 or Amide=NH2\002,//,\002 Use Residue=MOL to Begin a N"
	    "ew Strand,\002,\002 Residue=<CR> to End Entry\002)";
    static char fmt_20[] = "(/,\002 Enter Residue\002,i4,\002 :  \002,$)";
    static char fmt_30[] = "(a120)";
    static char fmt_50[] = "(/,\002 GETSEQ  --  Amino Acid Type \002,a3,\002"
	    " is Not Supported\002)";
    static char fmt_60[] = "(\002 GETSEQ  --  Error in Disulfide Bond\002"
	    ",\002 at Residue\002,i5)";
    static char fmt_70[] = "(\002 GETSEQ  --  Error in Disulfide Bond\002"
	    ",\002 at Residue\002,i5,\002 or\002,i5)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen),
	     s_rsfe(cilist *), e_rsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j;
    static char name__[3];
    static logical done;
    static char chir[1*10000];
    static integer next;
    static char record[120];
    static integer length;
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char string[120];
    extern /* Subroutine */ int getword_(char *, char *, integer *, ftnlen, 
	    ftnlen);

    /* Fortran I/O blocks */
    static cilist io___17 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_30, 0 };
    static icilist io___29 = { 1, string, 1, 0, 120, 1 };
    static cilist io___30 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_70, 0 };



#define chi_ref(a_1,a_2) phipsi_1.chi[(a_2)*4 + a_1 - 5]
#define seq_ref(a_0,a_1) &sequen_1.seq[(a_1)*3 + a_0 - 3]
#define amino_ref(a_0,a_1) &resdue_1.amino[(a_1)*3 + a_0 - 3]
#define ichain_ref(a_1,a_2) sequen_1.ichain[(a_2)*2 + a_1 - 3]



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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  phipsi.i  --  phi-psi-omega-chi angles for a protein  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     phi      value of the phi angle for each amino acid residue */
/*     psi      value of the psi angle for each amino acid residue */
/*     omega    value of the omega angle for each amino acid residue */
/*     chi      values of the chi angles for each amino acid residue */
/*     chiral   chirality of each amino acid residue (1=L, -1=D) */
/*     disulf   residue joined to each residue via a disulfide link */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  resdue.i  --  standard biopolymer residue abbreviations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     amino    three-letter abbreviations for amino acids types */
/*     nuclz    three-letter abbreviations for nucleic acids types */
/*     amino1   one-letter abbreviations for amino acids types */
/*     nuclz1   one-letter abbreviations for nucleic acids types */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  sequen.i  --  sequence information for a biopolymer  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nseq     total number of residues in biopolymer sequences */
/*     nchain   number of separate biopolymer sequence chains */
/*     ichain   first and last residue in each biopolymer chain */
/*     seqtyp   residue type for each residue in the sequence */
/*     seq      three-letter code for each residue in the sequence */
/*     chnnam   one-letter identifier for each sequence chain */




/*     provide a header to explain the method of sequence input */

    io___17.ciunit = iounit_1.iout;
    s_wsfe(&io___17);
    e_wsfe();

/*     initially, assume that only a single strand is present */

    sequen_1.nchain = 1;
    ichain_ref(1, 1) = 1;
    *(unsigned char *)&sequen_1.chnnam[0] = ' ';

/*     get the amino acid sequence data and dihedral angle values */

    i__ = 0;
    done = FALSE_;
    while(! done) {
	++i__;
	phipsi_1.phi[i__ - 1] = 0.;
	phipsi_1.psi[i__ - 1] = 0.;
	phipsi_1.omega[i__ - 1] = 0.;
	for (j = 1; j <= 4; ++j) {
	    chi_ref(j, i__) = 0.;
	}
	phipsi_1.disulf[i__ - 1] = 0;
	*(unsigned char *)&chir[i__ - 1] = ' ';
	io___22.ciunit = iounit_1.iout;
	s_wsfe(&io___22);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
	io___23.ciunit = iounit_1.input;
	s_rsfe(&io___23);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	upcase_(record, (ftnlen)120);
	next = 1;
	getword_(record, name__, &next, (ftnlen)120, (ftnlen)3);
	length = trimtext_(name__, (ftnlen)3);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	i__1 = s_rsli(&io___29);
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&phipsi_1.phi[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&phipsi_1.psi[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&phipsi_1.omega[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L40;
	}
	for (j = 1; j <= 4; ++j) {
	    i__1 = do_lio(&c__5, &c__1, (char *)&chi_ref(j, i__), (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L40;
	    }
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&phipsi_1.disulf[i__ - 1], (
		ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L40;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L40;
	}
L40:
	getword_(record, chir + (i__ - 1), &next, (ftnlen)120, (ftnlen)1);

/*     handle special names used for certain amino acids */

	if (s_cmp(name__, "CYH", (ftnlen)3, (ftnlen)3) == 0) {
	    s_copy(name__, "CYS", (ftnlen)3, (ftnlen)3);
	}
	if (s_cmp(name__, "CSS", (ftnlen)3, (ftnlen)3) == 0) {
	    s_copy(name__, "CYX", (ftnlen)3, (ftnlen)3);
	}
	if (s_cmp(name__, "HIP", (ftnlen)3, (ftnlen)3) == 0) {
	    s_copy(name__, "HIS", (ftnlen)3, (ftnlen)3);
	}

/*     disulfide bridged residues are cystine instead of cysteine */

	if (*(unsigned char *)name__ == 'C' && phipsi_1.disulf[i__ - 1] != 0) 
		{
	    length = 3;
	    s_copy(name__, "CYX", (ftnlen)3, (ftnlen)3);
	}

/*     process and store the current amino acid residue type */

	if (s_cmp(name__, "MOL", (ftnlen)3, (ftnlen)3) == 0) {
	    --i__;
	    ichain_ref(2, sequen_1.nchain) = i__;
	    ++sequen_1.nchain;
	    ichain_ref(1, sequen_1.nchain) = i__ + 1;
	} else {
	    if (s_cmp(name__, "   ", (ftnlen)3, (ftnlen)3) == 0) {
		done = TRUE_;
		sequen_1.nseq = i__ - 1;
		ichain_ref(2, sequen_1.nchain) = sequen_1.nseq;
	    } else {
		s_copy(seq_ref(0, i__), amino_ref(0, 31), (ftnlen)3, (ftnlen)
			3);
		sequen_1.seqtyp[i__ - 1] = 0;
		if (length == 1) {
		    for (j = 1; j <= 31; ++j) {
			if (*(unsigned char *)name__ == *(unsigned char *)&
				resdue_1.amino1[j - 1]) {
			    s_copy(seq_ref(0, i__), amino_ref(0, j), (ftnlen)
				    3, (ftnlen)3);
			    sequen_1.seqtyp[i__ - 1] = j;
			}
		    }
		} else if (length == 3) {
		    for (j = 1; j <= 31; ++j) {
			if (s_cmp(name__, amino_ref(0, j), (ftnlen)3, (ftnlen)
				3) == 0) {
			    s_copy(seq_ref(0, i__), amino_ref(0, j), (ftnlen)
				    3, (ftnlen)3);
			    sequen_1.seqtyp[i__ - 1] = j;
			}
		    }
		}
		if (sequen_1.seqtyp[i__ - 1] == 0) {
		    --i__;
		    io___30.ciunit = iounit_1.iout;
		    s_wsfe(&io___30);
		    do_fio(&c__1, name__, (ftnlen)3);
		    e_wsfe();
		}
	    }
	}
    }

/*     set chain identifiers if multiple chains are present */

    if (sequen_1.nchain > 1) {
	i__1 = sequen_1.nchain;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    *(unsigned char *)&sequen_1.chnnam[i__ - 1] = *(unsigned char *)&
		    ucase[i__ - 1];
	}
    }

/*     set default values for the phi-psi-omega-chi angles; */
/*     use extended values if no phi-psi values were given */

    i__1 = sequen_1.nseq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (phipsi_1.phi[i__ - 1] == 0. && phipsi_1.psi[i__ - 1] == 0.) {
	    phipsi_1.phi[i__ - 1] = -135.;
	    phipsi_1.psi[i__ - 1] = 135.;
	}
	if (phipsi_1.omega[i__ - 1] == 0.) {
	    phipsi_1.omega[i__ - 1] = 180.;
	}
	if (chi_ref(1, i__) == 0.) {
	    for (j = 1; j <= 4; ++j) {
		chi_ref(j, i__) = 180.;
		if (s_cmp(seq_ref(0, i__), "PRO", (ftnlen)3, (ftnlen)3) == 0) 
			{
		    chi_ref(j, i__) = 0.;
		}
		if (s_cmp(seq_ref(0, i__), "PCA", (ftnlen)3, (ftnlen)3) == 0) 
			{
		    chi_ref(j, i__) = 0.;
		}
	    }
	    if (s_cmp(seq_ref(0, i__), "PHE", (ftnlen)3, (ftnlen)3) == 0) {
		chi_ref(2, i__) = 90.;
	    }
	    if (s_cmp(seq_ref(0, i__), "TYR", (ftnlen)3, (ftnlen)3) == 0) {
		chi_ref(2, i__) = 90.;
	    }
	    if (s_cmp(seq_ref(0, i__), "TRP", (ftnlen)3, (ftnlen)3) == 0) {
		chi_ref(2, i__) = 90.;
	    }
	    if (s_cmp(seq_ref(0, i__), "HIS", (ftnlen)3, (ftnlen)3) == 0) {
		chi_ref(2, i__) = 90.;
	    }
	    if (s_cmp(seq_ref(0, i__), "HID", (ftnlen)3, (ftnlen)3) == 0) {
		chi_ref(2, i__) = 90.;
	    }
	    if (s_cmp(seq_ref(0, i__), "HIE", (ftnlen)3, (ftnlen)3) == 0) {
		chi_ref(2, i__) = 90.;
	    }
	}

/*     check for the presence of any disulfide bonds */

	if (phipsi_1.disulf[i__ - 1] != 0) {
	    if (s_cmp(seq_ref(0, i__), "CYX", (ftnlen)3, (ftnlen)3) != 0) {
		io___31.ciunit = iounit_1.iout;
		s_wsfe(&io___31);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	    if (i__ < phipsi_1.disulf[i__ - 1] && phipsi_1.disulf[
		    phipsi_1.disulf[i__ - 1] - 1] != i__) {
		io___32.ciunit = iounit_1.iout;
		s_wsfe(&io___32);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&phipsi_1.disulf[i__ - 1], (ftnlen)
			sizeof(integer));
		e_wsfe();
	    }
	}

/*     check the D/L chirality of the residues */

	if (*(unsigned char *)&chir[i__ - 1] == 'D') {
	    phipsi_1.chiral[i__ - 1] = -1;
	} else {
	    phipsi_1.chiral[i__ - 1] = 1;
	}
    }
    return 0;
} /* getseq_ */

#undef ichain_ref
#undef amino_ref
#undef seq_ref
#undef chi_ref




/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine prochain  --  build polypeptide backbone  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "prochain" builds up the internal coordinates for an amino */
/*     acid sequence from the phi, psi, omega and chi values */


/* Subroutine */ int prochain_(void)
{
    /* Initialized data */

    static integer ntyp[31] = { 1,7,15,27,41,55,65,77,87,96,107,122,138,161,
	    178,194,210,220,232,244,258,271,287,304,318,325,0,0,0,0,1 };
    static integer catyp[31] = { 2,8,16,28,42,56,66,78,88,97,108,123,139,162,
	    179,195,211,221,233,245,259,272,288,305,319,326,0,0,0,0,2 };
    static integer ctyp[31] = { 3,9,17,29,43,57,67,79,89,98,109,124,140,163,
	    180,196,212,222,234,246,260,273,289,306,320,327,0,0,0,0,3 };
    static integer hntyp[31] = { 4,10,18,30,44,58,68,80,90,0,110,125,141,164,
	    181,197,213,223,235,247,261,274,290,307,321,328,0,0,0,0,4 };
    static integer otyp[31] = { 5,11,19,31,45,59,69,81,91,99,111,126,142,165,
	    182,198,214,224,236,248,262,275,291,308,322,329,0,0,0,0,5 };
    static integer hatyp[31] = { 6,12,20,32,46,60,70,82,92,100,112,127,143,
	    166,183,199,215,225,237,249,263,276,292,309,0,330,0,0,0,0,6 };
    static integer nntyp[31] = { 350,356,362,368,374,380,386,392,398,404,412,
	    418,424,430,436,442,448,454,460,466,472,478,484,490,496,325,0,0,0,
	    0,350 };
    static integer cantyp[31] = { 351,357,363,369,375,381,387,393,399,405,413,
	    419,425,431,437,443,449,455,461,467,473,479,485,491,497,326,0,340,
	    0,0,351 };
    static integer cntyp[31] = { 352,358,364,370,376,382,388,394,400,406,414,
	    420,426,432,438,444,450,456,462,468,474,480,486,492,498,327,337,
	    342,0,0,352 };
    static integer hnntyp[31] = { 353,359,365,371,377,383,389,395,401,407,415,
	    421,427,433,439,445,451,457,463,469,475,481,487,493,499,328,0,0,0,
	    0,353 };
    static integer ontyp[31] = { 354,360,366,372,378,384,390,396,402,408,416,
	    422,428,434,440,446,452,458,464,470,476,482,488,494,500,329,339,
	    343,0,0,354 };
    static integer hantyp[31] = { 355,361,367,373,379,385,391,397,403,409,417,
	    423,429,435,441,447,453,459,465,471,477,483,489,495,0,330,338,341,
	    0,0,355 };
    static integer nctyp[31] = { 501,507,513,519,525,531,537,543,549,555,560,
	    566,572,578,584,590,596,602,608,614,620,626,632,638,644,0,0,0,344,
	    346,501 };
    static integer cactyp[31] = { 502,508,514,520,526,532,538,544,550,556,561,
	    567,573,579,585,591,597,603,609,615,621,627,633,639,645,0,0,0,0,
	    348,502 };
    static integer cctyp[31] = { 503,509,515,521,527,533,539,545,551,557,562,
	    568,574,580,586,592,598,604,610,616,622,628,634,640,646,0,0,0,0,0,
	    503 };
    static integer hnctyp[31] = { 504,510,516,522,528,534,540,546,552,0,563,
	    569,575,581,587,593,599,605,611,617,623,629,635,641,647,0,0,0,345,
	    347,504 };
    static integer octyp[31] = { 505,511,517,523,529,535,541,547,553,558,564,
	    570,576,582,588,594,600,606,612,618,624,630,636,642,648,0,0,0,0,0,
	    505 };
    static integer hactyp[31] = { 506,512,518,524,530,536,542,548,554,559,565,
	    571,577,583,589,595,601,607,613,619,625,631,637,643,0,0,0,0,0,349,
	    506 };

    /* Format strings */
    static char fmt_10[] = "(/,\002 Cyclize the Polypeptide Chain [N] :  "
	    "\002,$)";
    static char fmt_20[] = "(a120)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, k, m, ci[10000], ni[10000], cai[10000], next;
    extern /* Subroutine */ int zatom_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *);
    static logical cyclic;
    static char record[120];
    static logical single;
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static char answer[1], resname[3];
    extern /* Subroutine */ int proside_(char *, integer *, integer *, 
	    integer *, integer *, ftnlen), gettext_(char *, char *, integer *,
	     ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___52 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_20, 0 };



#define amino_ref(a_0,a_1) &resdue_1.amino[(a_1)*3 + a_0 - 3]
#define ichain_ref(a_1,a_2) sequen_1.ichain[(a_2)*2 + a_1 - 3]



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

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  phipsi.i  --  phi-psi-omega-chi angles for a protein  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     phi      value of the phi angle for each amino acid residue */
/*     psi      value of the psi angle for each amino acid residue */
/*     omega    value of the omega angle for each amino acid residue */
/*     chi      values of the chi angles for each amino acid residue */
/*     chiral   chirality of each amino acid residue (1=L, -1=D) */
/*     disulf   residue joined to each residue via a disulfide link */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  resdue.i  --  standard biopolymer residue abbreviations  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     amino    three-letter abbreviations for amino acids types */
/*     nuclz    three-letter abbreviations for nucleic acids types */
/*     amino1   one-letter abbreviations for amino acids types */
/*     nuclz1   one-letter abbreviations for nucleic acids types */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  sequen.i  --  sequence information for a biopolymer  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nseq     total number of residues in biopolymer sequences */
/*     nchain   number of separate biopolymer sequence chains */
/*     ichain   first and last residue in each biopolymer chain */
/*     seqtyp   residue type for each residue in the sequence */
/*     seq      three-letter code for each residue in the sequence */
/*     chnnam   one-letter identifier for each sequence chain */



/*     biopolymer atom types for amino acid backbone atoms */


/*     biopolymer atom types for N-terminal backbone atoms */


/*     biopolymer atom types for C-terminal backbone atoms */



/*     determine whether the peptide chain is cyclic */

    cyclic = FALSE_;
    io___52.ciunit = iounit_1.iout;
    s_wsfe(&io___52);
    e_wsfe();
    io___53.ciunit = iounit_1.input;
    s_rsfe(&io___53);
    do_fio(&c__1, record, (ftnlen)120);
    e_rsfe();
    next = 1;
    gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    upcase_(answer, (ftnlen)1);
    if (*(unsigned char *)answer == 'Y') {
	cyclic = TRUE_;
    }

/*     initialize the atom counter to the first atom */

    atoms_1.n = 1;

/*     set atom counter and the first residue number and type */

    i__1 = sequen_1.nchain;
    for (m = 1; m <= i__1; ++m) {
	single = FALSE_;
	if (ichain_ref(1, m) == ichain_ref(2, m)) {
	    single = TRUE_;
	}
	i__ = ichain_ref(1, m);
	k = sequen_1.seqtyp[i__ - 1];
	s_copy(resname, amino_ref(0, k), (ftnlen)3, (ftnlen)3);

/*     build the first residue for a cyclic peptide */

	if (cyclic) {
	    if (m == 1) {
		ni[i__ - 1] = atoms_1.n;
		zatom_(&ntyp[k - 1], &c_b84, &c_b84, &c_b84, &c__0, &c__0, &
			c__0, &c__0);
		cai[i__ - 1] = atoms_1.n;
		zatom_(&catyp[k - 1], &c_b91, &c_b84, &c_b84, &ni[i__ - 1], &
			c__0, &c__0, &c__0);
		ci[i__ - 1] = atoms_1.n;
		zatom_(&ctyp[k - 1], &c_b97, &c_b98, &c_b84, &cai[i__ - 1], &
			ni[i__ - 1], &c__0, &c__0);
	    } else {
		ni[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 1;
		i__3 = atoms_1.n - 2;
		i__4 = atoms_1.n - 3;
		zatom_(&ntyp[k - 1], &c_b102, &c_b103, &c_b104, &i__2, &i__3, 
			&i__4, &c__0);
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 1;
		zatom_(&c_n2, &c_b84, &c_b84, &c_b84, &i__2, &i__3, &c__0, &
			c__0);
		cai[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 3;
		zatom_(&catyp[k - 1], &c_b91, &c_b103, &c_b104, &ni[i__ - 1], 
			&i__2, &i__3, &c__0);
		ci[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 3;
		zatom_(&ctyp[k - 1], &c_b97, &c_b98, &c_b104, &cai[i__ - 1], &
			ni[i__ - 1], &i__2, &c__0);
	    }
	    d__1 = phipsi_1.psi[i__ - 1] - 180.;
	    zatom_(&otyp[k - 1], &c_b120, &c_b121, &d__1, &ci[i__ - 1], &cai[
		    i__ - 1], &ni[i__ - 1], &c__0);
	    d__1 = phipsi_1.phi[i__ - 1] - 180.;
	    zatom_(&hntyp[k - 1], &c_b123, &c_b124, &d__1, &ni[i__ - 1], &cai[
		    i__ - 1], &ci[i__ - 1], &c__0);
	    i__2 = -phipsi_1.chiral[i__ - 1];
	    zatom_(&hatyp[k - 1], &c_b126, &c_b127, &c_b128, &cai[i__ - 1], &
		    ni[i__ - 1], &ci[i__ - 1], &i__2);
	    proside_(resname, &i__, &cai[i__ - 1], &ni[i__ - 1], &ci[i__ - 1],
		     (ftnlen)3);

/*     build the first residue as an N-terminal formyl group */

	} else if (s_cmp(resname, "FOR", (ftnlen)3, (ftnlen)3) == 0) {
	    if (m == 1) {
		ci[i__ - 1] = atoms_1.n;
		zatom_(&cntyp[k - 1], &c_b84, &c_b84, &c_b84, &c__0, &c__0, &
			c__0, &c__0);
		ni[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 1;
		zatom_(&ontyp[k - 1], &c_b120, &c_b84, &c_b84, &i__2, &c__0, &
			c__0, &c__0);
		cai[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 1;
		zatom_(&hantyp[k - 1], &c_b143, &c_b144, &c_b84, &i__2, &i__3,
			 &c__0, &c__0);
	    } else {
		ci[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 1;
		i__3 = atoms_1.n - 2;
		i__4 = atoms_1.n - 3;
		zatom_(&cntyp[k - 1], &c_b102, &c_b103, &c_b104, &i__2, &i__3,
			 &i__4, &c__0);
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 1;
		zatom_(&c_n2, &c_b84, &c_b84, &c_b84, &i__2, &i__3, &c__0, &
			c__0);
		ni[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 1;
		i__3 = atoms_1.n - 2;
		i__4 = atoms_1.n - 3;
		zatom_(&ontyp[k - 1], &c_b120, &c_b103, &c_b104, &i__2, &i__3,
			 &i__4, &c__0);
		cai[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 1;
		i__4 = atoms_1.n - 3;
		zatom_(&hantyp[k - 1], &c_b143, &c_b144, &c_b84, &i__2, &i__3,
			 &i__4, &c__0);
	    }
	    phipsi_1.psi[i__ - 1] = 180.;

/*     build the first residue as an N-terminal acetyl group */

	} else if (s_cmp(resname, "ACE", (ftnlen)3, (ftnlen)3) == 0) {
	    if (m == 1) {
		cai[i__ - 1] = atoms_1.n;
		zatom_(&cantyp[k - 1], &c_b84, &c_b84, &c_b84, &c__0, &c__0, &
			c__0, &c__0);
		ci[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 1;
		zatom_(&cntyp[k - 1], &c_b97, &c_b84, &c_b84, &i__2, &c__0, &
			c__0, &c__0);
		i__2 = atoms_1.n - 1;
		i__3 = atoms_1.n - 2;
		zatom_(&ontyp[k - 1], &c_b120, &c_b121, &c_b84, &i__2, &i__3, 
			&c__0, &c__0);
	    } else {
		cai[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 1;
		i__3 = atoms_1.n - 2;
		i__4 = atoms_1.n - 3;
		zatom_(&cantyp[k - 1], &c_b102, &c_b103, &c_b104, &i__2, &
			i__3, &i__4, &c__0);
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 1;
		zatom_(&c_n2, &c_b84, &c_b84, &c_b84, &i__2, &i__3, &c__0, &
			c__0);
		ci[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 1;
		i__3 = atoms_1.n - 2;
		i__4 = atoms_1.n - 3;
		zatom_(&cntyp[k - 1], &c_b97, &c_b103, &c_b104, &i__2, &i__3, 
			&i__4, &c__0);
		i__2 = atoms_1.n - 1;
		i__3 = atoms_1.n - 2;
		i__4 = atoms_1.n - 3;
		zatom_(&ontyp[k - 1], &c_b120, &c_b121, &c_b84, &i__2, &i__3, 
			&i__4, &c__0);
	    }
	    ni[i__ - 1] = atoms_1.n;
	    i__2 = atoms_1.n - 3;
	    i__3 = atoms_1.n - 2;
	    i__4 = atoms_1.n - 1;
	    zatom_(&hantyp[k - 1], &c_b126, &c_b128, &c_b84, &i__2, &i__3, &
		    i__4, &c__0);
	    i__2 = atoms_1.n - 4;
	    i__3 = atoms_1.n - 3;
	    i__4 = atoms_1.n - 1;
	    zatom_(&hantyp[k - 1], &c_b126, &c_b128, &c_b209, &i__2, &i__3, &
		    i__4, &c__1);
	    i__2 = atoms_1.n - 5;
	    i__3 = atoms_1.n - 4;
	    i__4 = atoms_1.n - 2;
	    zatom_(&hantyp[k - 1], &c_b126, &c_b128, &c_b209, &i__2, &i__3, &
		    i__4, &c_n1);
	    phipsi_1.psi[i__ - 1] = 180.;

/*     build the first residue as a proline */

	} else if (s_cmp(resname, "PRO", (ftnlen)3, (ftnlen)3) == 0) {
	    if (m == 1) {
		ni[i__ - 1] = atoms_1.n;
		zatom_(&nntyp[k - 1], &c_b84, &c_b84, &c_b84, &c__0, &c__0, &
			c__0, &c__0);
		cai[i__ - 1] = atoms_1.n;
		zatom_(&cantyp[k - 1], &c_b223, &c_b84, &c_b84, &ni[i__ - 1], 
			&c__0, &c__0, &c__0);
		ci[i__ - 1] = atoms_1.n;
		if (single) {
		    zatom_(&cctyp[k - 1], &c_b97, &c_b230, &c_b84, &cai[i__ - 
			    1], &ni[i__ - 1], &c__0, &c__0);
		} else {
		    zatom_(&cntyp[k - 1], &c_b97, &c_b98, &c_b84, &cai[i__ - 
			    1], &ni[i__ - 1], &c__0, &c__0);
		}
	    } else {
		ni[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 1;
		i__3 = atoms_1.n - 2;
		i__4 = atoms_1.n - 3;
		zatom_(&nntyp[k - 1], &c_b102, &c_b103, &c_b104, &i__2, &i__3,
			 &i__4, &c__0);
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 1;
		zatom_(&c_n2, &c_b84, &c_b84, &c_b84, &i__2, &i__3, &c__0, &
			c__0);
		cai[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 3;
		zatom_(&cantyp[k - 1], &c_b223, &c_b103, &c_b104, &ni[i__ - 1]
			, &i__2, &i__3, &c__0);
		ci[i__ - 1] = atoms_1.n;
		if (single) {
		    i__2 = atoms_1.n - 3;
		    zatom_(&cctyp[k - 1], &c_b97, &c_b230, &c_b104, &cai[i__ 
			    - 1], &ni[i__ - 1], &i__2, &c__0);
		} else {
		    i__2 = atoms_1.n - 3;
		    zatom_(&cntyp[k - 1], &c_b97, &c_b98, &c_b104, &cai[i__ - 
			    1], &ni[i__ - 1], &i__2, &c__0);
		}
	    }
	    if (single) {
		d__1 = phipsi_1.psi[0] - 180.;
		zatom_(&octyp[k - 1], &c_b261, &c_b262, &d__1, &ci[i__ - 1], &
			cai[i__ - 1], &ni[i__ - 1], &c__0);
	    } else {
		d__1 = phipsi_1.psi[0] - 180.;
		zatom_(&ontyp[k - 1], &c_b120, &c_b121, &d__1, &ci[i__ - 1], &
			cai[i__ - 1], &ni[i__ - 1], &c__0);
	    }
	    zatom_(&hnntyp[k - 1], &c_b123, &c_b127, &c_b84, &ni[i__ - 1], &
		    cai[i__ - 1], &ci[i__ - 1], &c__0);
	    zatom_(&hnntyp[k - 1], &c_b123, &c_b127, &c_b273, &ni[i__ - 1], &
		    cai[i__ - 1], &ci[i__ - 1], &c__0);
	    i__2 = -phipsi_1.chiral[i__ - 1];
	    zatom_(&hantyp[k - 1], &c_b126, &c_b127, &c_b128, &cai[i__ - 1], &
		    ni[i__ - 1], &ci[i__ - 1], &i__2);
	    proside_(resname, &i__, &cai[i__ - 1], &ni[i__ - 1], &ci[i__ - 1],
		     (ftnlen)3);

/*     build the first residue as a pyroglutamic acid */

	} else if (s_cmp(resname, "PCA", (ftnlen)3, (ftnlen)3) == 0) {
	    if (m == 1) {
		ni[i__ - 1] = atoms_1.n;
		zatom_(&nntyp[k - 1], &c_b84, &c_b84, &c_b84, &c__0, &c__0, &
			c__0, &c__0);
		cai[i__ - 1] = atoms_1.n;
		zatom_(&cantyp[k - 1], &c_b223, &c_b84, &c_b84, &ni[i__ - 1], 
			&c__0, &c__0, &c__0);
		ci[i__ - 1] = atoms_1.n;
		if (single) {
		    zatom_(&cctyp[k - 1], &c_b97, &c_b230, &c_b84, &cai[i__ - 
			    1], &ni[i__ - 1], &c__0, &c__0);
		} else {
		    zatom_(&cntyp[k - 1], &c_b97, &c_b98, &c_b84, &cai[i__ - 
			    1], &ni[i__ - 1], &c__0, &c__0);
		}
	    } else {
		ni[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 1;
		i__3 = atoms_1.n - 2;
		i__4 = atoms_1.n - 3;
		zatom_(&nntyp[k - 1], &c_b102, &c_b103, &c_b104, &i__2, &i__3,
			 &i__4, &c__0);
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 1;
		zatom_(&c_n2, &c_b84, &c_b84, &c_b84, &i__2, &i__3, &c__0, &
			c__0);
		cai[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 3;
		zatom_(&cantyp[k - 1], &c_b223, &c_b103, &c_b104, &ni[i__ - 1]
			, &i__2, &i__3, &c__0);
		ci[i__ - 1] = atoms_1.n;
		if (single) {
		    i__2 = atoms_1.n - 3;
		    zatom_(&cctyp[k - 1], &c_b97, &c_b230, &c_b104, &cai[i__ 
			    - 1], &ni[i__ - 1], &i__2, &c__0);
		} else {
		    i__2 = atoms_1.n - 3;
		    zatom_(&cntyp[k - 1], &c_b97, &c_b98, &c_b104, &cai[i__ - 
			    1], &ni[i__ - 1], &i__2, &c__0);
		}
	    }
	    if (single) {
		d__1 = phipsi_1.psi[0] - 180.;
		zatom_(&octyp[k - 1], &c_b261, &c_b262, &d__1, &ci[i__ - 1], &
			cai[i__ - 1], &ni[i__ - 1], &c__0);
	    } else {
		d__1 = phipsi_1.psi[0] - 180.;
		zatom_(&ontyp[k - 1], &c_b120, &c_b121, &d__1, &ci[i__ - 1], &
			cai[i__ - 1], &ni[i__ - 1], &c__0);
	    }
	    zatom_(&hnntyp[k - 1], &c_b123, &c_b127, &c_b332, &ni[i__ - 1], &
		    cai[i__ - 1], &ci[i__ - 1], &c__0);
	    i__2 = -phipsi_1.chiral[i__ - 1];
	    zatom_(&hantyp[k - 1], &c_b126, &c_b127, &c_b128, &cai[i__ - 1], &
		    ni[i__ - 1], &ci[i__ - 1], &i__2);
	    proside_(resname, &i__, &cai[i__ - 1], &ni[i__ - 1], &ci[i__ - 1],
		     (ftnlen)3);

/*     build the first residue for all other standard amino acids */

	} else {
	    if (m == 1) {
		ni[i__ - 1] = atoms_1.n;
		zatom_(&nntyp[k - 1], &c_b84, &c_b84, &c_b84, &c__0, &c__0, &
			c__0, &c__0);
		cai[i__ - 1] = atoms_1.n;
		zatom_(&cantyp[k - 1], &c_b223, &c_b84, &c_b84, &ni[i__ - 1], 
			&c__0, &c__0, &c__0);
		ci[i__ - 1] = atoms_1.n;
		if (single) {
		    zatom_(&cctyp[k - 1], &c_b97, &c_b230, &c_b84, &cai[i__ - 
			    1], &ni[i__ - 1], &c__0, &c__0);
		} else {
		    zatom_(&cntyp[k - 1], &c_b97, &c_b98, &c_b84, &cai[i__ - 
			    1], &ni[i__ - 1], &c__0, &c__0);
		}
	    } else {
		ni[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 1;
		i__3 = atoms_1.n - 2;
		i__4 = atoms_1.n - 3;
		zatom_(&nntyp[k - 1], &c_b102, &c_b103, &c_b104, &i__2, &i__3,
			 &i__4, &c__0);
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 1;
		zatom_(&c_n2, &c_b84, &c_b84, &c_b84, &i__2, &i__3, &c__0, &
			c__0);
		cai[i__ - 1] = atoms_1.n;
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 3;
		zatom_(&cantyp[k - 1], &c_b223, &c_b103, &c_b104, &ni[i__ - 1]
			, &i__2, &i__3, &c__0);
		ci[i__ - 1] = atoms_1.n;
		if (single) {
		    i__2 = atoms_1.n - 3;
		    zatom_(&cctyp[k - 1], &c_b97, &c_b230, &c_b104, &cai[i__ 
			    - 1], &ni[i__ - 1], &i__2, &c__0);
		} else {
		    i__2 = atoms_1.n - 3;
		    zatom_(&cntyp[k - 1], &c_b97, &c_b98, &c_b104, &cai[i__ - 
			    1], &ni[i__ - 1], &i__2, &c__0);
		}
	    }
	    if (single) {
		d__1 = phipsi_1.psi[0] - 180.;
		zatom_(&octyp[k - 1], &c_b261, &c_b262, &d__1, &ci[i__ - 1], &
			cai[i__ - 1], &ni[i__ - 1], &c__0);
	    } else {
		d__1 = phipsi_1.psi[0] - 180.;
		zatom_(&ontyp[k - 1], &c_b120, &c_b121, &d__1, &ci[i__ - 1], &
			cai[i__ - 1], &ni[i__ - 1], &c__0);
	    }
	    zatom_(&hnntyp[k - 1], &c_b123, &c_b127, &phipsi_1.phi[i__ - 1], &
		    ni[i__ - 1], &cai[i__ - 1], &ci[i__ - 1], &c__0);
	    i__2 = atoms_1.n - 1;
	    zatom_(&hnntyp[k - 1], &c_b123, &c_b127, &c_b393, &ni[i__ - 1], &
		    cai[i__ - 1], &i__2, &c__1);
	    i__2 = atoms_1.n - 2;
	    zatom_(&hnntyp[k - 1], &c_b123, &c_b127, &c_b393, &ni[i__ - 1], &
		    cai[i__ - 1], &i__2, &c_n1);
	    i__2 = -phipsi_1.chiral[i__ - 1];
	    zatom_(&hantyp[k - 1], &c_b126, &c_b127, &c_b128, &cai[i__ - 1], &
		    ni[i__ - 1], &ci[i__ - 1], &i__2);
	    proside_(resname, &i__, &cai[i__ - 1], &ni[i__ - 1], &ci[i__ - 1],
		     (ftnlen)3);
	}

/*     build atoms for residues in the middle of the chain */

	i__2 = ichain_ref(2, m) - 1;
	for (i__ = ichain_ref(1, m) + 1; i__ <= i__2; ++i__) {
	    k = sequen_1.seqtyp[i__ - 1];
	    s_copy(resname, amino_ref(0, k), (ftnlen)3, (ftnlen)3);
	    ni[i__ - 1] = atoms_1.n;
	    zatom_(&ntyp[k - 1], &c_b402, &c_b403, &phipsi_1.psi[i__ - 2], &
		    ci[i__ - 2], &cai[i__ - 2], &ni[i__ - 2], &c__0);
	    cai[i__ - 1] = atoms_1.n;
	    zatom_(&catyp[k - 1], &c_b91, &c_b124, &phipsi_1.omega[i__ - 2], &
		    ni[i__ - 1], &ci[i__ - 2], &cai[i__ - 2], &c__0);
	    ci[i__ - 1] = atoms_1.n;
	    zatom_(&ctyp[k - 1], &c_b97, &c_b230, &phipsi_1.phi[i__ - 1], &
		    cai[i__ - 1], &ni[i__ - 1], &ci[i__ - 2], &c__0);
	    d__1 = phipsi_1.psi[i__ - 1] - 180.;
	    zatom_(&otyp[k - 1], &c_b120, &c_b121, &d__1, &ci[i__ - 1], &cai[
		    i__ - 1], &ni[i__ - 1], &c__0);
	    d__1 = phipsi_1.phi[i__ - 1] - 180.;
	    zatom_(&hntyp[k - 1], &c_b123, &c_b124, &d__1, &ni[i__ - 1], &cai[
		    i__ - 1], &ci[i__ - 1], &c__0);
	    i__3 = -phipsi_1.chiral[i__ - 1];
	    zatom_(&hatyp[k - 1], &c_b126, &c_b127, &c_b128, &cai[i__ - 1], &
		    ni[i__ - 1], &ci[i__ - 1], &i__3);
	    proside_(resname, &i__, &cai[i__ - 1], &ni[i__ - 1], &ci[i__ - 1],
		     (ftnlen)3);
	}

/*     set the number and type of the last residue */

	i__ = ichain_ref(2, m);
	k = sequen_1.seqtyp[i__ - 1];
	s_copy(resname, amino_ref(0, k), (ftnlen)3, (ftnlen)3);

/*     build the last residue for a cyclic peptide */

	if (cyclic) {
	    ni[i__ - 1] = atoms_1.n;
	    zatom_(&ntyp[k - 1], &c_b402, &c_b403, &phipsi_1.psi[i__ - 2], &
		    ci[i__ - 2], &cai[i__ - 2], &ni[i__ - 2], &c__0);
	    cai[i__ - 1] = atoms_1.n;
	    zatom_(&catyp[k - 1], &c_b91, &c_b124, &phipsi_1.omega[i__ - 2], &
		    ni[i__ - 1], &ci[i__ - 2], &cai[i__ - 2], &c__0);
	    ci[i__ - 1] = atoms_1.n;
	    zatom_(&ctyp[k - 1], &c_b97, &c_b230, &phipsi_1.phi[i__ - 1], &
		    cai[i__ - 1], &ni[i__ - 1], &ci[i__ - 2], &c__0);
	    zatom_(&c_n1, &c_b84, &c_b84, &c_b84, ni, &ci[i__ - 1], &c__0, &
		    c__0);
	    d__1 = phipsi_1.psi[i__ - 1] - 180.;
	    zatom_(&otyp[k - 1], &c_b120, &c_b121, &d__1, &ci[i__ - 1], &cai[
		    i__ - 1], &ni[i__ - 1], &c__0);
	    d__1 = phipsi_1.phi[i__ - 1] - 180.;
	    zatom_(&hntyp[k - 1], &c_b123, &c_b124, &d__1, &ni[i__ - 1], &cai[
		    i__ - 1], &ci[i__ - 1], &c__0);
	    i__2 = -phipsi_1.chiral[i__ - 1];
	    zatom_(&hatyp[k - 1], &c_b126, &c_b127, &c_b128, &cai[i__ - 1], &
		    ni[i__ - 1], &ci[i__ - 1], &i__2);
	    proside_(resname, &i__, &cai[i__ - 1], &ni[i__ - 1], &ci[i__ - 1],
		     (ftnlen)3);

/*     build the last residue as a C-terminal amide */

	} else if (s_cmp(resname, "NH2", (ftnlen)3, (ftnlen)3) == 0) {
	    zatom_(&nctyp[k - 1], &c_b402, &c_b403, &phipsi_1.psi[i__ - 2], &
		    ci[i__ - 2], &cai[i__ - 2], &ni[i__ - 2], &c__0);
	    i__2 = atoms_1.n - 1;
	    zatom_(&hnctyp[k - 1], &c_b123, &c_b449, &c_b84, &i__2, &ci[i__ - 
		    2], &cai[i__ - 2], &c__0);
	    i__2 = atoms_1.n - 2;
	    zatom_(&hnctyp[k - 1], &c_b123, &c_b449, &c_b104, &i__2, &ci[i__ 
		    - 2], &cai[i__ - 2], &c__0);

/*     build the last residue as a C-terminal N-methylamide */

	} else if (s_cmp(resname, "NME", (ftnlen)3, (ftnlen)3) == 0) {
	    zatom_(&nctyp[k - 1], &c_b402, &c_b403, &phipsi_1.psi[i__ - 2], &
		    ci[i__ - 2], &cai[i__ - 2], &ni[i__ - 2], &c__0);
	    i__2 = atoms_1.n - 1;
	    zatom_(&cactyp[k - 1], &c_b91, &c_b124, &c_b104, &i__2, &ci[i__ - 
		    2], &cai[i__ - 2], &c__0);
	    i__2 = atoms_1.n - 2;
	    i__3 = atoms_1.n - 1;
	    zatom_(&hnctyp[k - 1], &c_b123, &c_b465, &c_b124, &i__2, &ci[i__ 
		    - 2], &i__3, &c__1);
	    i__2 = atoms_1.n - 2;
	    i__3 = atoms_1.n - 3;
	    zatom_(&hactyp[k - 1], &c_b126, &c_b127, &c_b104, &i__2, &i__3, &
		    ci[i__ - 2], &c__0);
	    i__2 = atoms_1.n - 3;
	    i__3 = atoms_1.n - 4;
	    i__4 = atoms_1.n - 1;
	    zatom_(&hactyp[k - 1], &c_b126, &c_b127, &c_b127, &i__2, &i__3, &
		    i__4, &c__1);
	    i__2 = atoms_1.n - 4;
	    i__3 = atoms_1.n - 5;
	    i__4 = atoms_1.n - 2;
	    zatom_(&hactyp[k - 1], &c_b126, &c_b127, &c_b127, &i__2, &i__3, &
		    i__4, &c_n1);

/*     build the last residue for all other standard amino acids */

	} else {
	    if (! single) {
		ni[i__ - 1] = atoms_1.n;
		zatom_(&nctyp[k - 1], &c_b402, &c_b403, &phipsi_1.psi[i__ - 2]
			, &ci[i__ - 2], &cai[i__ - 2], &ni[i__ - 2], &c__0);
		cai[i__ - 1] = atoms_1.n;
		zatom_(&cactyp[k - 1], &c_b91, &c_b124, &phipsi_1.omega[i__ - 
			2], &ni[i__ - 1], &ci[i__ - 2], &cai[i__ - 2], &c__0);
		ci[i__ - 1] = atoms_1.n;
		zatom_(&cctyp[k - 1], &c_b97, &c_b230, &phipsi_1.phi[i__ - 1],
			 &cai[i__ - 1], &ni[i__ - 1], &ci[i__ - 2], &c__0);
		d__1 = phipsi_1.psi[i__ - 1] - 180.;
		zatom_(&octyp[k - 1], &c_b261, &c_b262, &d__1, &ci[i__ - 1], &
			cai[i__ - 1], &ni[i__ - 1], &c__0);
		d__1 = phipsi_1.phi[i__ - 1] - 180.;
		zatom_(&hnctyp[k - 1], &c_b123, &c_b124, &d__1, &ni[i__ - 1], 
			&cai[i__ - 1], &ci[i__ - 1], &c__0);
		i__2 = -phipsi_1.chiral[i__ - 1];
		zatom_(&hactyp[k - 1], &c_b126, &c_b127, &c_b128, &cai[i__ - 
			1], &ni[i__ - 1], &ci[i__ - 1], &i__2);
		proside_(resname, &i__, &cai[i__ - 1], &ni[i__ - 1], &ci[i__ 
			- 1], (ftnlen)3);
	    }
	    zatom_(&octyp[k - 1], &c_b261, &c_b262, &phipsi_1.psi[i__ - 1], &
		    ci[i__ - 1], &cai[i__ - 1], &ni[i__ - 1], &c__0);
	}
    }

/*     finally, set the total number of atoms */

    --atoms_1.n;
    return 0;
} /* prochain_ */

#undef ichain_ref
#undef amino_ref




/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine proside  --  build amino acid side chain  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "proside" builds the side chain for a single amino acid */
/*     residue in terms of internal coordinates */

/*     resname   3-letter name of current amino acid residue */
/*     i         number of the current amino acid residue */
/*     cai       atom number of alpha carbon in residue i */
/*     ni        atom number of amide nitrogen in residue i */
/*     ci        atom number of carbonyl carbon in residue i */


/* Subroutine */ int proside_(char *resname, integer *i__, integer *cai, 
	integer *ni, integer *ci, ftnlen resname_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int zatom_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *);


#define chi_ref(a_1,a_2) phipsi_1.chi[(a_2)*4 + a_1 - 5]



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
/*     ##  phipsi.i  --  phi-psi-omega-chi angles for a protein  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     phi      value of the phi angle for each amino acid residue */
/*     psi      value of the psi angle for each amino acid residue */
/*     omega    value of the omega angle for each amino acid residue */
/*     chi      values of the chi angles for each amino acid residue */
/*     chiral   chirality of each amino acid residue (1=L, -1=D) */
/*     disulf   residue joined to each residue via a disulfide link */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  sequen.i  --  sequence information for a biopolymer  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nseq     total number of residues in biopolymer sequences */
/*     nchain   number of separate biopolymer sequence chains */
/*     ichain   first and last residue in each biopolymer chain */
/*     seqtyp   residue type for each residue in the sequence */
/*     seq      three-letter code for each residue in the sequence */
/*     chnnam   one-letter identifier for each sequence chain */




/*     glycine residue  (GLY) */

    if (s_cmp(resname, "GLY", (ftnlen)3, (ftnlen)3) == 0) {
	if (*i__ == 1) {
	    zatom_(&c__355, &c_b126, &c_b127, &c_b128, cai, ni, ci, &
		    phipsi_1.chiral[*i__ - 1]);
	} else if (*i__ == sequen_1.nseq) {
	    zatom_(&c__506, &c_b126, &c_b127, &c_b128, cai, ni, ci, &
		    phipsi_1.chiral[*i__ - 1]);
	} else {
	    zatom_(&c__6, &c_b126, &c_b127, &c_b128, cai, ni, ci, &
		    phipsi_1.chiral[*i__ - 1]);
	}

/*     alanine residue  (ALA) */

    } else if (s_cmp(resname, "ALA", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__13, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__14, &c_b126, &c_b209, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	zatom_(&c__14, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 2;
	zatom_(&c__14, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);

/*     valine residue  (VAL) */

    } else if (s_cmp(resname, "VAL", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__21, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__23, &c_b517, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	zatom_(&c__25, &c_b517, &c_b127, &c_b127, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 2;
	zatom_(&c__22, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	zatom_(&c__24, &c_b126, &c_b209, &c_b104, &i__1, &i__2, cai, &c__0);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 1;
	zatom_(&c__24, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 2;
	zatom_(&c__24, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c_n1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	zatom_(&c__26, &c_b126, &c_b209, &c_b104, &i__1, &i__2, cai, &c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 1;
	zatom_(&c__26, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 2;
	zatom_(&c__26, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c_n1);

/*     leucine residue  (LEU) */

    } else if (s_cmp(resname, "LEU", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__33, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__35, &c_b517, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__37, &c_b517, &c_b127, &chi_ref(2, *i__), &i__1, &i__2, cai,
		 &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__39, &c_b517, &c_b127, &c_b209, &i__1, &i__2, &i__3, &c_n1);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	zatom_(&c__34, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	zatom_(&c__34, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 4;
	zatom_(&c__36, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	zatom_(&c__38, &c_b126, &c_b209, &c_b104, &i__1, &i__2, &i__3, &c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 1;
	zatom_(&c__38, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 2;
	zatom_(&c__38, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c_n1);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 10;
	zatom_(&c__40, &c_b126, &c_b209, &c_b104, &i__1, &i__2, &i__3, &c__0);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 1;
	zatom_(&c__40, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 2;
	zatom_(&c__40, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c_n1);

/*     isoleucine residue  (ILE) */

    } else if (s_cmp(resname, "ILE", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__47, &c_b517, &c_b127, &c_b127, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__49, &c_b517, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	zatom_(&c__51, &c_b517, &c_b127, &c_b127, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	zatom_(&c__53, &c_b517, &c_b127, &chi_ref(2, *i__), &i__1, &i__2, cai,
		 &c__0);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	zatom_(&c__48, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 2;
	zatom_(&c__50, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 3;
	zatom_(&c__50, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c_n1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 6;
	zatom_(&c__52, &c_b126, &c_b681, &c_b104, &i__1, &i__2, &i__3, &c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 1;
	zatom_(&c__52, &c_b126, &c_b681, &c_b687, &i__1, &i__2, &i__3, &c__1);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 2;
	zatom_(&c__52, &c_b126, &c_b681, &c_b687, &i__1, &i__2, &i__3, &c_n1);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 10;
	zatom_(&c__54, &c_b126, &c_b681, &c_b104, &i__1, &i__2, &i__3, &c__0);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 1;
	zatom_(&c__54, &c_b126, &c_b681, &c_b687, &i__1, &i__2, &i__3, &c__1);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 2;
	zatom_(&c__54, &c_b126, &c_b681, &c_b687, &i__1, &i__2, &i__3, &c_n1);

/*     serine residue  (SER) */

    } else if (s_cmp(resname, "SER", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__61, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__63, &c_b715, &c_b716, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	zatom_(&c__62, &c_b126, &c_b209, &c_b721, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 2;
	zatom_(&c__62, &c_b126, &c_b209, &c_b721, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	zatom_(&c__64, &c_b729, &c_b730, &chi_ref(2, *i__), &i__1, &i__2, cai,
		 &c__0);

/*     threonine residue  (THR) */

    } else if (s_cmp(resname, "THR", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__71, &c_b517, &c_b127, &c_b127, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__73, &c_b715, &c_b716, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	zatom_(&c__75, &c_b517, &c_b127, &c_b744, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 2;
	zatom_(&c__72, &c_b126, &c_b209, &c_b721, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	zatom_(&c__74, &c_b729, &c_b730, &chi_ref(2, *i__), &i__1, &i__2, cai,
		 &c__0);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 5;
	zatom_(&c__76, &c_b126, &c_b681, &c_b104, &i__1, &i__2, cai, &c__0);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 1;
	zatom_(&c__76, &c_b126, &c_b681, &c_b687, &i__1, &i__2, &i__3, &c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 2;
	zatom_(&c__76, &c_b126, &c_b681, &c_b687, &i__1, &i__2, &i__3, &c_n1);

/*     cysteine residue  (CYS) */

    } else if (s_cmp(resname, "CYS", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__83, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__85, &c_b776, &c_b687, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	zatom_(&c__84, &c_b126, &c_b209, &c_b782, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 2;
	zatom_(&c__84, &c_b126, &c_b209, &c_b782, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	zatom_(&c__86, &c_b402, &c_b791, &chi_ref(2, *i__), &i__1, &i__2, cai,
		 &c__0);

/*     cystine residue  (CYX) */

    } else if (s_cmp(resname, "CYX", (ftnlen)3, (ftnlen)3) == 0) {
	if (phipsi_1.disulf[*i__ - 1] > *i__) {
	    zatom_(&c__93, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		    phipsi_1.chiral[*i__ - 1]);
	    i__1 = atoms_1.n - 1;
	    zatom_(&c__95, &c_b776, &c_b687, &chi_ref(1, *i__), &i__1, cai, 
		    ni, &c__0);
	    i__1 = atoms_1.n - 2;
	    i__2 = atoms_1.n - 1;
	    zatom_(&c__94, &c_b126, &c_b209, &c_b782, &i__1, cai, &i__2, &
		    c__1);
	    i__1 = atoms_1.n - 3;
	    i__2 = atoms_1.n - 2;
	    zatom_(&c__94, &c_b126, &c_b209, &c_b782, &i__1, cai, &i__2, &
		    c_n1);
	    phipsi_1.disulf[*i__ - 1] = atoms_1.n - 3;
	} else if (phipsi_1.disulf[*i__ - 1] < *i__) {
	    zatom_(&c__93, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		    phipsi_1.chiral[*i__ - 1]);
	    i__1 = atoms_1.n - 1;
	    zatom_(&c__95, &c_b776, &c_b687, &chi_ref(1, *i__), &i__1, cai, 
		    ni, &c__0);
	    i__1 = atoms_1.n - 2;
	    i__2 = atoms_1.n - 1;
	    zatom_(&c__94, &c_b126, &c_b209, &c_b782, &i__1, cai, &i__2, &
		    c__1);
	    i__1 = atoms_1.n - 3;
	    i__2 = atoms_1.n - 2;
	    zatom_(&c__94, &c_b126, &c_b209, &c_b782, &i__1, cai, &i__2, &
		    c_n1);
	    i__1 = atoms_1.n - 3;
	    zatom_(&c_n1, &c_b84, &c_b84, &c_b84, &phipsi_1.disulf[
		    phipsi_1.disulf[*i__ - 1] - 1], &i__1, &c__0, &c__0);
	}

/*     proline residue  (PRO) */

    } else if (s_cmp(resname, "PRO", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__101, &c_b517, &c_b839, &c_b127, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__103, &c_b517, &c_b839, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	if (*i__ == 1) {
	    i__1 = atoms_1.n - 1;
	    i__2 = atoms_1.n - 2;
	    zatom_(&c__410, &c_b517, &c_b839, &chi_ref(2, *i__), &i__1, &i__2,
		     cai, &c__0);
	} else {
	    i__1 = atoms_1.n - 1;
	    i__2 = atoms_1.n - 2;
	    zatom_(&c__105, &c_b517, &c_b839, &chi_ref(2, *i__), &i__1, &i__2,
		     cai, &c__0);
	}
	i__1 = atoms_1.n - 1;
	zatom_(&c_n1, &c_b84, &c_b84, &c_b84, ni, &i__1, &c__0, &c__0);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 2;
	zatom_(&c__102, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	zatom_(&c__102, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 3;
	zatom_(&c__104, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 4;
	zatom_(&c__104, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c_n1)
		;
	if (*i__ == 1) {
	    i__1 = atoms_1.n - 5;
	    i__2 = atoms_1.n - 6;
	    zatom_(&c__411, &c_b126, &c_b209, &c_b209, &i__1, &i__2, ni, &
		    c__1);
	    i__1 = atoms_1.n - 6;
	    i__2 = atoms_1.n - 7;
	    zatom_(&c__411, &c_b126, &c_b209, &c_b209, &i__1, &i__2, ni, &
		    c_n1);
	} else {
	    i__1 = atoms_1.n - 5;
	    i__2 = atoms_1.n - 6;
	    zatom_(&c__106, &c_b126, &c_b209, &c_b209, &i__1, &i__2, ni, &
		    c__1);
	    i__1 = atoms_1.n - 6;
	    i__2 = atoms_1.n - 7;
	    zatom_(&c__106, &c_b126, &c_b209, &c_b209, &i__1, &i__2, ni, &
		    c_n1);
	}

/*     phenylalanine residue  (PHE) */

    } else if (s_cmp(resname, "PHE", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__113, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__115, &c_b223, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__116, &c_b909, &c_b144, &chi_ref(2, *i__), &i__1, &i__2, 
		cai, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__116, &c_b909, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__118, &c_b909, &c_b144, &c_b104, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 5;
	zatom_(&c__118, &c_b909, &c_b144, &c_b104, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 5;
	zatom_(&c__120, &c_b909, &c_b144, &c_b84, &i__1, &i__2, &i__3, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	zatom_(&c_n1, &c_b84, &c_b84, &c_b84, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 6;
	zatom_(&c__114, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 7;
	zatom_(&c__114, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 5;
	zatom_(&c__117, &c_b949, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 5;
	zatom_(&c__117, &c_b949, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 5;
	zatom_(&c__119, &c_b949, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 6;
	zatom_(&c__119, &c_b949, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 8;
	zatom_(&c__121, &c_b949, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;

/*     tyrosine residue  (TYR) */

    } else if (s_cmp(resname, "TYR", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__128, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__130, &c_b223, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__131, &c_b909, &c_b144, &chi_ref(2, *i__), &i__1, &i__2, 
		cai, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__131, &c_b909, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__133, &c_b909, &c_b144, &c_b104, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 5;
	zatom_(&c__133, &c_b909, &c_b144, &c_b104, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 5;
	zatom_(&c__135, &c_b909, &c_b144, &c_b84, &i__1, &i__2, &i__3, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	zatom_(&c_n1, &c_b84, &c_b84, &c_b84, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__136, &c_b1013, &c_b144, &c_b144, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 7;
	zatom_(&c__129, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 8;
	zatom_(&c__129, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 6;
	zatom_(&c__132, &c_b949, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 6;
	zatom_(&c__132, &c_b949, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 6;
	zatom_(&c__134, &c_b949, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 7;
	zatom_(&c__134, &c_b949, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	zatom_(&c__137, &c_b1048, &c_b393, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;

/*     tryptophan residue  (TRP) */

    } else if (s_cmp(resname, "TRP", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__144, &c_b517, &c_b127, &c_b127, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__146, &c_b223, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__147, &c_b1062, &c_b1063, &chi_ref(2, *i__), &i__1, &i__2, 
		cai, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__149, &c_b1062, &c_b1063, &c_b393, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__150, &c_b1062, &c_b393, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__152, &c_b1062, &c_b393, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 1;
	zatom_(&c_n1, &c_b84, &c_b84, &c_b84, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 1;
	i__3 = atoms_1.n - 2;
	zatom_(&c__153, &c_b1062, &c_b144, &c_b104, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 1;
	zatom_(&c__155, &c_b1062, &c_b144, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 3;
	zatom_(&c__157, &c_b1062, &c_b144, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 6;
	zatom_(&c__159, &c_b1062, &c_b144, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	zatom_(&c_n1, &c_b84, &c_b84, &c_b84, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 9;
	zatom_(&c__145, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 11;
	i__2 = atoms_1.n - 10;
	zatom_(&c__145, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 8;
	zatom_(&c__148, &c_b949, &c_b1063, &c_b1063, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 8;
	zatom_(&c__151, &c_b1128, &c_b1063, &c_b1063, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 6;
	zatom_(&c__154, &c_b949, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 6;
	zatom_(&c__156, &c_b949, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 7;
	zatom_(&c__158, &c_b949, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 9;
	zatom_(&c__160, &c_b949, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;

/*     histidine (HD and HE) residue  (HIS) */

    } else if (s_cmp(resname, "HIS", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__167, &c_b517, &c_b127, &c_b127, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__169, &c_b223, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__170, &c_b1062, &c_b1063, &chi_ref(2, *i__), &i__1, &i__2, 
		cai, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__172, &c_b1062, &c_b1063, &c_b393, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__174, &c_b1062, &c_b393, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 3;
	zatom_(&c__176, &c_b1062, &c_b393, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	zatom_(&c_n1, &c_b84, &c_b84, &c_b84, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	zatom_(&c__168, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 6;
	zatom_(&c__168, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	zatom_(&c__171, &c_b123, &c_b1063, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 4;
	zatom_(&c__173, &c_b949, &c_b1063, &c_b1063, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 5;
	zatom_(&c__175, &c_b949, &c_b1063, &c_b1063, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 7;
	zatom_(&c__177, &c_b123, &c_b1063, &c_b1063, &i__1, &i__2, &i__3, &
		c__1);

/*     histidine (HD only) residue  (HID) */

    } else if (s_cmp(resname, "HID", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__184, &c_b517, &c_b127, &c_b127, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__186, &c_b223, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__187, &c_b1062, &c_b1063, &chi_ref(2, *i__), &i__1, &i__2, 
		cai, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__189, &c_b1062, &c_b1063, &c_b393, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__191, &c_b1062, &c_b393, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 3;
	zatom_(&c__193, &c_b1062, &c_b393, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	zatom_(&c_n1, &c_b84, &c_b84, &c_b84, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	zatom_(&c__185, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 6;
	zatom_(&c__185, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	zatom_(&c__188, &c_b123, &c_b1063, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 4;
	zatom_(&c__190, &c_b949, &c_b1063, &c_b1063, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 5;
	zatom_(&c__192, &c_b949, &c_b1063, &c_b1063, &i__1, &i__2, &i__3, &
		c__1);

/*     histidine (HE only) residue  (HIE) */

    } else if (s_cmp(resname, "HIE", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__200, &c_b517, &c_b127, &c_b127, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__202, &c_b223, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__203, &c_b1062, &c_b1063, &chi_ref(2, *i__), &i__1, &i__2, 
		cai, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__204, &c_b1062, &c_b1063, &c_b393, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__206, &c_b1062, &c_b393, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 3;
	zatom_(&c__208, &c_b1062, &c_b393, &c_b84, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	zatom_(&c_n1, &c_b84, &c_b84, &c_b84, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	zatom_(&c__201, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 6;
	zatom_(&c__201, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 3;
	zatom_(&c__205, &c_b949, &c_b1063, &c_b1063, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 4;
	zatom_(&c__207, &c_b949, &c_b1063, &c_b1063, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 6;
	zatom_(&c__209, &c_b123, &c_b1063, &c_b1063, &i__1, &i__2, &i__3, &
		c__1);

/*     aspartic acid residue  (ASP) */

    } else if (s_cmp(resname, "ASP", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__216, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__218, &c_b97, &c_b519, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__219, &c_b261, &c_b262, &chi_ref(2, *i__), &i__1, &i__2, 
		cai, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__219, &c_b261, &c_b262, &c_b1063, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	zatom_(&c__217, &c_b126, &c_b209, &c_b128, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	zatom_(&c__217, &c_b126, &c_b209, &c_b128, &i__1, cai, &i__2, &c_n1);

/*     asparagine residue  (ASN) */

    } else if (s_cmp(resname, "ASN", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__226, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__228, &c_b97, &c_b519, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__229, &c_b120, &c_b121, &chi_ref(2, *i__), &i__1, &i__2, 
		cai, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__230, &c_b402, &c_b403, &c_b1378, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	zatom_(&c__227, &c_b126, &c_b209, &c_b128, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	zatom_(&c__227, &c_b126, &c_b209, &c_b128, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 6;
	zatom_(&c__231, &c_b123, &c_b449, &c_b84, &i__1, &i__2, &i__3, &c__0);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 1;
	zatom_(&c__231, &c_b123, &c_b449, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;

/*     glutamic acid residue  (GLU) */

    } else if (s_cmp(resname, "GLU", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__238, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__240, &c_b517, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__242, &c_b97, &c_b519, &chi_ref(2, *i__), &i__1, &i__2, cai,
		 &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__243, &c_b261, &c_b262, &chi_ref(3, *i__), &i__1, &i__2, &
		i__3, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__243, &c_b261, &c_b262, &c_b1063, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	zatom_(&c__239, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	zatom_(&c__239, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 5;
	zatom_(&c__241, &c_b126, &c_b209, &c_b128, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 6;
	zatom_(&c__241, &c_b126, &c_b209, &c_b128, &i__1, &i__2, &i__3, &c_n1)
		;

/*     glutamine residue  (GLN) */

    } else if (s_cmp(resname, "GLN", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__250, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__252, &c_b517, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__254, &c_b97, &c_b519, &chi_ref(2, *i__), &i__1, &i__2, cai,
		 &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__255, &c_b120, &c_b121, &chi_ref(3, *i__), &i__1, &i__2, &
		i__3, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__256, &c_b402, &c_b403, &c_b1378, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	zatom_(&c__251, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	zatom_(&c__251, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 5;
	zatom_(&c__253, &c_b126, &c_b209, &c_b128, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 6;
	zatom_(&c__253, &c_b126, &c_b209, &c_b128, &i__1, &i__2, &i__3, &c_n1)
		;
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	zatom_(&c__257, &c_b123, &c_b449, &c_b84, &i__1, &i__2, &i__3, &c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 1;
	zatom_(&c__257, &c_b123, &c_b449, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;

/*     methionine residue  (MET) */

    } else if (s_cmp(resname, "MET", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__264, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__266, &c_b517, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__268, &c_b776, &c_b687, &chi_ref(2, *i__), &i__1, &i__2, 
		cai, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__269, &c_b776, &c_b1509, &chi_ref(3, *i__), &i__1, &i__2, &
		i__3, &c__0);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	zatom_(&c__265, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	zatom_(&c__265, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 4;
	zatom_(&c__267, &c_b126, &c_b209, &c_b782, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 5;
	zatom_(&c__267, &c_b126, &c_b209, &c_b782, &i__1, &i__2, &i__3, &c_n1)
		;
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	zatom_(&c__270, &c_b126, &c_b782, &c_b104, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 1;
	zatom_(&c__270, &c_b126, &c_b782, &c_b209, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 2;
	zatom_(&c__270, &c_b126, &c_b782, &c_b209, &i__1, &i__2, &i__3, &c_n1)
		;

/*     lysine residue  (LYS) */

    } else if (s_cmp(resname, "LYS", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__277, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__279, &c_b517, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__281, &c_b517, &c_b127, &chi_ref(2, *i__), &i__1, &i__2, 
		cai, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__283, &c_b517, &c_b127, &chi_ref(3, *i__), &i__1, &i__2, &
		i__3, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__285, &c_b223, &c_b127, &chi_ref(4, *i__), &i__1, &i__2, &
		i__3, &c__0);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	zatom_(&c__278, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	zatom_(&c__278, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 5;
	zatom_(&c__280, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 6;
	zatom_(&c__280, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c_n1)
		;
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 6;
	zatom_(&c__282, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 7;
	zatom_(&c__282, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c_n1)
		;
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 7;
	zatom_(&c__284, &c_b126, &c_b209, &c_b1600, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 8;
	zatom_(&c__284, &c_b126, &c_b209, &c_b1600, &i__1, &i__2, &i__3, &
		c_n1);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 11;
	zatom_(&c__286, &c_b123, &c_b127, &c_b104, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 1;
	zatom_(&c__286, &c_b123, &c_b127, &c_b127, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 11;
	i__2 = atoms_1.n - 12;
	i__3 = atoms_1.n - 2;
	zatom_(&c__286, &c_b123, &c_b127, &c_b127, &i__1, &i__2, &i__3, &c_n1)
		;

/*     arginine residue  (ARG) */

    } else if (s_cmp(resname, "ARG", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__293, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__295, &c_b517, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__297, &c_b517, &c_b127, &chi_ref(2, *i__), &i__1, &i__2, 
		cai, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__299, &c_b1636, &c_b127, &chi_ref(3, *i__), &i__1, &i__2, &
		i__3, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__301, &c_b1062, &c_b144, &chi_ref(4, *i__), &i__1, &i__2, &
		i__3, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__302, &c_b1062, &c_b144, &c_b104, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 1;
	zatom_(&c__302, &c_b1062, &c_b144, &c_b144, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 6;
	zatom_(&c__294, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 7;
	zatom_(&c__294, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 7;
	zatom_(&c__296, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 8;
	zatom_(&c__296, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c_n1)
		;
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 8;
	zatom_(&c__298, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 9;
	zatom_(&c__298, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c_n1)
		;
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 9;
	zatom_(&c__300, &c_b123, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 11;
	zatom_(&c__303, &c_b123, &c_b144, &c_b104, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 1;
	zatom_(&c__303, &c_b123, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 12;
	i__3 = atoms_1.n - 13;
	zatom_(&c__303, &c_b123, &c_b144, &c_b104, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 11;
	i__2 = atoms_1.n - 13;
	i__3 = atoms_1.n - 1;
	zatom_(&c__303, &c_b123, &c_b144, &c_b144, &i__1, &i__2, &i__3, &c__1)
		;

/*     ornithine residue  (ORN) */

    } else if (s_cmp(resname, "ORN", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__310, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__312, &c_b517, &c_b127, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__314, &c_b517, &c_b127, &chi_ref(2, *i__), &i__1, &i__2, 
		cai, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__316, &c_b223, &c_b127, &chi_ref(3, *i__), &i__1, &i__2, &
		i__3, &c__0);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	zatom_(&c__311, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	zatom_(&c__311, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 4;
	zatom_(&c__313, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 5;
	zatom_(&c__313, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c_n1)
		;
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 5;
	zatom_(&c__315, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 6;
	zatom_(&c__315, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c_n1)
		;
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	zatom_(&c__317, &c_b123, &c_b127, &c_b104, &i__1, &i__2, &i__3, &c__0)
		;
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 1;
	zatom_(&c__317, &c_b123, &c_b127, &c_b127, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 2;
	zatom_(&c__317, &c_b123, &c_b127, &c_b127, &i__1, &i__2, &i__3, &c_n1)
		;

/*     methylalanine residue  (AIB) */

    } else if (s_cmp(resname, "AIB", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = -phipsi_1.chiral[*i__ - 1];
	zatom_(&c__323, &c_b517, &c_b127, &c_b519, cai, ni, ci, &i__1);
	zatom_(&c__323, &c_b517, &c_b127, &c_b519, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 2;
	zatom_(&c__324, &c_b126, &c_b209, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 1;
	zatom_(&c__324, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 2;
	zatom_(&c__324, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 4;
	zatom_(&c__324, &c_b126, &c_b209, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 1;
	zatom_(&c__324, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 2;
	zatom_(&c__324, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);

/*     pyroglutamic acid residue  (PCA) */

    } else if (s_cmp(resname, "PCA", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__331, &c_b517, &c_b839, &c_b127, cai, ni, ci, &
		phipsi_1.chiral[*i__ - 1]);
	i__1 = atoms_1.n - 1;
	zatom_(&c__333, &c_b517, &c_b839, &chi_ref(1, *i__), &i__1, cai, ni, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__335, &c_b517, &c_b839, &chi_ref(2, *i__), &i__1, &i__2, 
		cai, &c__0);
	i__1 = atoms_1.n - 1;
	zatom_(&c_n1, &c_b84, &c_b84, &c_b84, ni, &i__1, &c__0, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__336, &c_b120, &c_b1063, &c_b1063, &i__1, ni, &i__2, &c__1);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	zatom_(&c__332, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	zatom_(&c__332, &c_b126, &c_b209, &c_b209, &i__1, cai, &i__2, &c_n1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 4;
	zatom_(&c__334, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c__1)
		;
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 5;
	zatom_(&c__334, &c_b126, &c_b209, &c_b209, &i__1, &i__2, &i__3, &c_n1)
		;

/*     unknown residue  (UNK) */

    } else if (s_cmp(resname, "UNK", (ftnlen)3, (ftnlen)3) == 0) {
	if (*i__ == 1) {
	    zatom_(&c__355, &c_b126, &c_b127, &c_b128, cai, ni, ci, &
		    phipsi_1.chiral[*i__ - 1]);
	} else if (*i__ == sequen_1.nseq) {
	    zatom_(&c__506, &c_b126, &c_b127, &c_b128, cai, ni, ci, &
		    phipsi_1.chiral[*i__ - 1]);
	} else {
	    zatom_(&c__6, &c_b126, &c_b127, &c_b128, cai, ni, ci, &
		    phipsi_1.chiral[*i__ - 1]);
	}
    }
    return 0;
} /* proside_ */

#undef chi_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine pauling  --  pack multiple polypeptide chains  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "pauling" uses a rigid body optimization to approximately */
/*     pack multiple polypeptide chains */


/* Subroutine */ int pauling_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int rigidxyz_(void);
    static integer i__, j, k;
    static doublereal xx[1000];
    static integer nvar;
    extern /* Subroutine */ int ocvm_(integer *, doublereal *, doublereal *, 
	    doublereal *, D_fp, U_fp);
    static doublereal grdmin;
    extern /* Subroutine */ int potoff_(void), orient_(void);
    static doublereal minimum;
    extern /* Subroutine */ int optsave_();
    extern doublereal pauling1_();


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]
#define gfix_ref(a_1,a_2) kgeoms_1.gfix[(a_2)*3 + a_1 - 4]
#define imol_ref(a_1,a_2) molcul_1.imol[(a_2)*2 + a_1 - 3]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]
#define pfix_ref(a_1,a_2) kgeoms_1.pfix[(a_2)*2 + a_1 - 3]
#define wgrp_ref(a_1,a_2) group_1.wgrp[(a_2)*1001 + a_1 - 0]
#define igfix_ref(a_1,a_2) kgeoms_1.igfix[(a_2)*2 + a_1 - 3]
#define kpfix_ref(a_1,a_2) kgeoms_1.kpfix[(a_2)*3 + a_1 - 4]



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
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  couple.i  --  near-neighbor atom connectivity lists   ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     maxn13   maximum number of atoms 1-3 connected to an atom */
/*     maxn14   maximum number of atoms 1-4 connected to an atom */
/*     maxn15   maximum number of atoms 1-5 connected to an atom */

/*     n12      number of atoms directly bonded to each atom */
/*     i12      atom numbers of atoms 1-2 connected to each atom */
/*     n13      number of atoms in a 1-3 relation to each atom */
/*     i13      atom numbers of atoms 1-3 connected to each atom */
/*     n14      number of atoms in a 1-4 relation to each atom */
/*     i14      atom numbers of atoms 1-4 connected to each atom */
/*     n15      number of atoms in a 1-5 relation to each atom */
/*     i15      atom numbers of atoms 1-5 connected to each atom */




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  group.i  --  partitioning of system into atom groups  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     grpmass     total mass of all the atoms in each group */
/*     wgrp        weight for each set of group-group interactions */
/*     ngrp        total number of atom groups in the system */
/*     kgrp        contiguous list of the atoms in each group */
/*     igrp        first and last atom of each group in the list */
/*     grplist     number of the group to which each atom belongs */
/*     use_group   flag to use partitioning of system into groups */
/*     use_intra   flag to include only intragroup interactions */
/*     use_inter   flag to include only intergroup interactions */




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
/*     ##  kgeoms.i  --  parameters for the geometrical restraints  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     xpfix      x-coordinate target for each restrained position */
/*     ypfix      y-coordinate target for each restrained position */
/*     zpfix      z-coordinate target for each restrained position */
/*     pfix       force constant and flat-well range for each position */
/*     dfix       force constant and target range for each distance */
/*     afix       force constant and target range for each angle */
/*     tfix       force constant and target range for each torsion */
/*     gfix       force constant and target range for each group distance */
/*     chir       force constant and target range for chiral centers */
/*     depth      depth of shallow Gaussian basin restraint */
/*     width      exponential width coefficient of Gaussian basin */
/*     rwall      radius of spherical droplet boundary restraint */
/*     npfix      number of position restraints to be applied */
/*     ipfix      atom number involved in each position restraint */
/*     kpfix      flags to use x-, y-, z-coordinate position restraints */
/*     ndfix      number of distance restraints to be applied */
/*     idfix      atom numbers defining each distance restraint */
/*     nafix      number of angle restraints to be applied */
/*     iafix      atom numbers defining each angle restraint */
/*     ntfix      number of torsional restraints to be applied */
/*     itfix      atom numbers defining each torsional restraint */
/*     ngfix      number of group distance restraints to be applied */
/*     igfix      group numbers defining each group distance restraint */
/*     nchir      number of chirality restraints to be applied */
/*     ichir      atom numbers defining each chirality restraint */
/*     use_basin  logical flag governing use of Gaussian basin */
/*     use_wall   logical flag governing use of droplet boundary */




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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  rigid.i  --  rigid body coordinates for atom groups  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     xrb         rigid body reference x-coordinate for each atom */
/*     yrb         rigid body reference y-coordinate for each atom */
/*     zrb         rigid body reference z-coordinate for each atom */
/*     rbc         current rigid body coordinates for each group */
/*     use_rigid   flag to mark use of rigid body coordinate system */




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




/*     set all atoms to be active during energy evaluations */

    usage_1.nuse = atoms_1.n;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	usage_1.use[i__ - 1] = TRUE_;
    }

/*     only geometric restraints will by used in optimization */

    potoff_();
    potent_1.use_geom__ = TRUE_;

/*     set the default values for the restraint variables */

    kgeoms_1.npfix = 0;
    kgeoms_1.ndfix = 0;
    kgeoms_1.ntfix = 0;
    kgeoms_1.ngfix = 0;
    kgeoms_1.nchir = 0;
    kgeoms_1.use_basin__ = TRUE_;
    kgeoms_1.depth = 3.;
    kgeoms_1.width = 1.5;
    kgeoms_1.use_wall__ = FALSE_;

/*     assign each chain to a separate molecule-based group */

    group_1.use_group__ = TRUE_;
    group_1.ngrp = molcul_1.nmol;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	igrp_ref(1, i__) = imol_ref(1, i__);
	igrp_ref(2, i__) = imol_ref(2, i__);
	i__2 = igrp_ref(2, i__);
	for (j = igrp_ref(1, i__); j <= i__2; ++j) {
	    group_1.kgrp[j - 1] = molcul_1.kmol[j - 1];
	    group_1.grplist[group_1.kgrp[j - 1] - 1] = i__;
	}
    }
    i__1 = group_1.ngrp;
    for (i__ = 0; i__ <= i__1; ++i__) {
	i__2 = group_1.ngrp;
	for (j = 0; j <= i__2; ++j) {
	    wgrp_ref(j, i__) = 1.;
	}
	wgrp_ref(i__, i__) = 1.;
    }

/*     assume unit mass for each atom and set group masses */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atmtyp_1.mass[i__ - 1] = 1.;
    }
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	group_1.grpmass[i__ - 1] = (doublereal) (igrp_ref(2, i__) - igrp_ref(
		1, i__) + 1);
    }

/*     set pairwise restraints between the centers of chains */

    i__1 = group_1.ngrp - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = group_1.ngrp;
	for (j = i__ + 1; j <= i__2; ++j) {
	    ++kgeoms_1.ngfix;
	    igfix_ref(1, kgeoms_1.ngfix) = i__;
	    igfix_ref(2, kgeoms_1.ngfix) = j;
	    gfix_ref(1, kgeoms_1.ngfix) = 1.;
	    gfix_ref(2, kgeoms_1.ngfix) = (doublereal) (j - i__) * 11.;
	    gfix_ref(3, kgeoms_1.ngfix) = (doublereal) (j - i__) * 11.;
	}
    }

/*     set position restraints on alpha carbons of each chain */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (katoms_1.atmnum[atoms_1.type__[i__ - 1] - 1] == 6) {
	    i__2 = couple_1.n12[i__ - 1];
	    for (j = 1; j <= i__2; ++j) {
		if (katoms_1.atmnum[atoms_1.type__[i12_ref(j, i__) - 1] - 1] 
			== 7) {
		    i__3 = couple_1.n13[i__ - 1];
		    for (k = 1; k <= i__3; ++k) {
			if (katoms_1.atmnum[atoms_1.type__[i13_ref(k, i__) - 
				1] - 1] == 8) {
			    ++kgeoms_1.npfix;
			    kgeoms_1.ipfix[kgeoms_1.npfix - 1] = i__;
			    kpfix_ref(1, kgeoms_1.npfix) = 1;
			    kpfix_ref(2, kgeoms_1.npfix) = 1;
			    kpfix_ref(3, kgeoms_1.npfix) = 0;
			    kgeoms_1.xpfix[kgeoms_1.npfix - 1] = (doublereal) 
				    (group_1.grplist[i__ - 1] - 1) * 11.;
			    kgeoms_1.ypfix[kgeoms_1.npfix - 1] = 0.;
			    kgeoms_1.zpfix[kgeoms_1.npfix - 1] = 0.;
			    pfix_ref(1, kgeoms_1.npfix) = 1.;
			    pfix_ref(2, kgeoms_1.npfix) = 0.;
			    goto L10;
			}
		    }
		}
	    }
	}
L10:
	;
    }

/*     get rigid body reference coordinates for each chain */

    orient_();

/*     transfer rigid body coordinates to optimization parameters */

    nvar = 0;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    ++nvar;
	    xx[nvar - 1] = rbc_ref(j, i__);
	}
    }

/*     make the call to the optimization routine */

    inform_1.iprint = 0;
    inform_1.iwrite = 0;
    grdmin = .1;
    s_copy(output_1.coordtype, "NONE", (ftnlen)9, (ftnlen)4);
    ocvm_(&nvar, xx, &minimum, &grdmin, (D_fp)pauling1_, (U_fp)optsave_);

/*     transfer optimization parameters to rigid body coordinates */

    nvar = 0;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    ++nvar;
	    rbc_ref(j, i__) = xx[nvar - 1];
	}
    }

/*     convert from rigid body to Cartesian coordinates */

    rigidxyz_();
    return 0;
} /* pauling_ */

#undef kpfix_ref
#undef igfix_ref
#undef wgrp_ref
#undef pfix_ref
#undef igrp_ref
#undef imol_ref
#undef gfix_ref
#undef rbc_ref
#undef i13_ref
#undef i12_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  function pauling1  --  energy and gradient for pauling  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "pauling1" is a service routine that computes the energy */
/*     and gradient for optimally conditioned variable metric */
/*     optimization of rigid bodies */


doublereal pauling1_(doublereal *xx, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    extern /* Subroutine */ int rigidxyz_(void);
    static doublereal e;
    static integer i__, j, nvar;
    static doublereal derivs[6000]	/* was [6][1000] */;
    extern /* Subroutine */ int gradrgd_(doublereal *, doublereal *);


#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]
#define derivs_ref(a_1,a_2) derivs[(a_2)*6 + a_1 - 7]



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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  group.i  --  partitioning of system into atom groups  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     grpmass     total mass of all the atoms in each group */
/*     wgrp        weight for each set of group-group interactions */
/*     ngrp        total number of atom groups in the system */
/*     kgrp        contiguous list of the atoms in each group */
/*     igrp        first and last atom of each group in the list */
/*     grplist     number of the group to which each atom belongs */
/*     use_group   flag to use partitioning of system into groups */
/*     use_intra   flag to include only intragroup interactions */
/*     use_inter   flag to include only intergroup interactions */




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
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  rigid.i  --  rigid body coordinates for atom groups  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     xrb         rigid body reference x-coordinate for each atom */
/*     yrb         rigid body reference y-coordinate for each atom */
/*     zrb         rigid body reference z-coordinate for each atom */
/*     rbc         current rigid body coordinates for each group */
/*     use_rigid   flag to mark use of rigid body coordinate system */




/*     translate optimization parameters to rigid body coordinates */

    /* Parameter adjustments */
    --g;
    --xx;

    /* Function Body */
    nvar = 0;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    ++nvar;
	    rbc_ref(j, i__) = xx[nvar];
	}
    }

/*     compute and store the energy and gradient */

    rigidxyz_();
    gradrgd_(&e, derivs);
    ret_val = e;

/*     translate rigid body gradient to optimization gradient */

    nvar = 0;
    i__1 = group_1.ngrp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    ++nvar;
	    g[nvar] = derivs_ref(j, i__);
	}
    }
    return ret_val;
} /* pauling1_ */

#undef derivs_ref
#undef rbc_ref


/* Main program alias */ int protein_ () { MAIN__ (); return 0; }
