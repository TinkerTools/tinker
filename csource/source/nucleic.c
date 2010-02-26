/* nucleic.f -- translated by f2c (version 20050501).
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
    integer n12[25000], i12[200000]	/* was [8][25000] */, n13[25000], i13[
	    600000]	/* was [24][25000] */, n14[25000], i14[1800000]	/* 
	    was [72][25000] */, n15[25000], i15[5400000]	/* was [216][
	    25000] */;
} couple_;

#define couple_1 couple_

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
    doublereal bkbone[60000]	/* was [6][10000] */, glyco[10000];
    integer pucker[10000];
    logical dblhlx, deoxy[10000];
    char hlxform[1];
} nucleo_;

#define nucleo_1 nucleo_

struct {
    integer ltitle;
    char title[120];
} titles_;

#define titles_1 titles_

struct {
    char amino[93], nuclz[36], amino1[31], nuclz1[12];
} resdue_;

#define resdue_1 resdue_

struct {
    integer nseq, nchain, ichain[20000]	/* was [2][10000] */, seqtyp[10000];
    char seq[30000], chnnam[10000];
} sequen_;

#define sequen_1 sequen_

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
static doublereal c_b108 = 0.;
static integer c__0 = 0;
static doublereal c_b115 = 1.52;
static doublereal c_b122 = 113.;
static doublereal c_b126 = 30.;
static doublereal c_b127 = 150.;
static doublereal c_b128 = 180.;
static integer c_n2 = -2;
static doublereal c_b148 = 1.63;
static doublereal c_b149 = 106.;
static integer c_n1 = -1;
static doublereal c_b161 = 1.44;
static doublereal c_b168 = 110.1;
static doublereal c_b172 = .96;
static doublereal c_b183 = 119.;
static doublereal c_b190 = 1.46;
static doublereal c_b191 = 108.9;
static doublereal c_b193 = 1.42;
static doublereal c_b194 = 109.8;
static doublereal c_b195 = 145.;
static doublereal c_b199 = 107.;
static doublereal c_b203 = 140.;
static doublereal c_b205 = 1.53;
static doublereal c_b206 = 115.9;
static doublereal c_b209 = 102.4;
static doublereal c_b218 = 112.1;
static doublereal c_b229 = 1.43;
static doublereal c_b230 = 109.5;
static doublereal c_b237 = 1.09;
static doublereal c_b251 = 120.;
static doublereal c_b255 = 115.;
static doublereal c_b259 = 90.;
static doublereal c_b281 = 1.6;
static doublereal c_b284 = 1.48;
static doublereal c_b285 = 109.;
static doublereal c_b291 = 101.8;
static integer c__1251 = 1251;
static integer c__1239 = 1239;
static doublereal c_b390 = 60.;
static doublereal c_b394 = -60.;
static integer c__1017 = 1017;
static doublereal c_b505 = 108.1;
static doublereal c_b506 = 113.7;
static integer c__1021 = 1021;
static doublereal c_b509 = 1.37;
static doublereal c_b510 = 128.4;
static integer c__1020 = 1020;
static doublereal c_b513 = 1.3;
static doublereal c_b514 = 113.8;
static integer c__1019 = 1019;
static doublereal c_b518 = 1.39;
static doublereal c_b519 = 104.;
static integer c__1025 = 1025;
static doublereal c_b523 = 1.4;
static doublereal c_b524 = 132.4;
static integer c__1027 = 1027;
static doublereal c_b528 = 1.34;
static doublereal c_b529 = 123.5;
static integer c__1024 = 1024;
static doublereal c_b533 = 1.35;
static doublereal c_b534 = 117.4;
static integer c__1023 = 1023;
static doublereal c_b538 = 1.33;
static doublereal c_b539 = 118.8;
static integer c__1022 = 1022;
static doublereal c_b543 = 1.32;
static doublereal c_b544 = 129.2;
static integer c__1018 = 1018;
static doublereal c_b549 = 110.9;
static integer c__1030 = 1030;
static doublereal c_b565 = 1.08;
static doublereal c_b566 = 123.1;
static integer c__1028 = 1028;
static doublereal c_b570 = 1.;
static integer c__1029 = 1029;
static integer c__1026 = 1026;
static doublereal c_b581 = 115.4;
static integer c__1047 = 1047;
static integer c__1051 = 1051;
static doublereal c_b591 = 1.38;
static integer c__1050 = 1050;
static doublereal c_b595 = 1.31;
static doublereal c_b596 = 114.;
static integer c__1049 = 1049;
static doublereal c_b601 = 103.8;
static integer c__1055 = 1055;
static doublereal c_b606 = 130.1;
static integer c__1060 = 1060;
static doublereal c_b610 = 1.23;
static doublereal c_b611 = 128.8;
static integer c__1054 = 1054;
static doublereal c_b616 = 111.4;
static integer c__1053 = 1053;
static doublereal c_b621 = 125.2;
static integer c__1057 = 1057;
static doublereal c_b626 = 116.1;
static integer c__1052 = 1052;
static doublereal c_b631 = 123.3;
static integer c__1048 = 1048;
static doublereal c_b635 = 1.36;
static doublereal c_b636 = 112.3;
static integer c__1061 = 1061;
static doublereal c_b653 = 123.;
static integer c__1056 = 1056;
static integer c__1058 = 1058;
static integer c__1059 = 1059;
static integer c__1078 = 1078;
static integer c__1079 = 1079;
static doublereal c_b679 = 117.8;
static integer c__1084 = 1084;
static doublereal c_b682 = 1.24;
static doublereal c_b683 = 118.9;
static integer c__1080 = 1080;
static doublereal c_b688 = 118.7;
static integer c__1081 = 1081;
static doublereal c_b693 = 120.6;
static integer c__1085 = 1085;
static doublereal c_b698 = 118.3;
static integer c__1082 = 1082;
static doublereal c_b703 = 121.6;
static integer c__1083 = 1083;
static doublereal c_b708 = 116.9;
static integer c__1086 = 1086;
static integer c__1087 = 1087;
static integer c__1088 = 1088;
static integer c__1089 = 1089;
static doublereal c_b734 = 119.5;
static integer c__1106 = 1106;
static integer c__1107 = 1107;
static doublereal c_b745 = 117.1;
static integer c__1112 = 1112;
static doublereal c_b748 = 1.22;
static doublereal c_b749 = 123.2;
static integer c__1108 = 1108;
static doublereal c_b754 = 114.8;
static integer c__1109 = 1109;
static doublereal c_b759 = 127.;
static integer c__1114 = 1114;
static doublereal c_b764 = 119.8;
static integer c__1110 = 1110;
static doublereal c_b769 = 114.7;
static integer c__1111 = 1111;
static doublereal c_b774 = 119.2;
static integer c__1113 = 1113;
static doublereal c_b785 = 116.5;
static integer c__1115 = 1115;
static doublereal c_b790 = 120.4;
static integer c__1116 = 1116;
static doublereal c_b795 = 118.6;
static integer c__1132 = 1132;
static integer c__1136 = 1136;
static integer c__1135 = 1135;
static integer c__1134 = 1134;
static integer c__1140 = 1140;
static integer c__1142 = 1142;
static integer c__1139 = 1139;
static integer c__1138 = 1138;
static integer c__1137 = 1137;
static integer c__1133 = 1133;
static integer c__1145 = 1145;
static integer c__1143 = 1143;
static integer c__1144 = 1144;
static integer c__1141 = 1141;
static integer c__1161 = 1161;
static integer c__1165 = 1165;
static integer c__1164 = 1164;
static integer c__1163 = 1163;
static integer c__1169 = 1169;
static integer c__1174 = 1174;
static integer c__1168 = 1168;
static integer c__1167 = 1167;
static integer c__1171 = 1171;
static integer c__1166 = 1166;
static integer c__1162 = 1162;
static integer c__1175 = 1175;
static integer c__1170 = 1170;
static integer c__1172 = 1172;
static integer c__1173 = 1173;
static integer c__1191 = 1191;
static integer c__1192 = 1192;
static integer c__1197 = 1197;
static integer c__1193 = 1193;
static integer c__1194 = 1194;
static integer c__1198 = 1198;
static integer c__1195 = 1195;
static integer c__1196 = 1196;
static integer c__1199 = 1199;
static integer c__1200 = 1200;
static integer c__1201 = 1201;
static integer c__1202 = 1202;
static integer c__1218 = 1218;
static integer c__1219 = 1219;
static integer c__1224 = 1224;
static doublereal c_b1045 = 122.9;
static integer c__1220 = 1220;
static integer c__1221 = 1221;
static doublereal c_b1055 = 126.4;
static integer c__1226 = 1226;
static doublereal c_b1060 = 120.5;
static integer c__1222 = 1222;
static doublereal c_b1065 = 114.1;
static integer c__1227 = 1227;
static doublereal c_b1069 = 1.5;
static doublereal c_b1070 = 117.5;
static integer c__1223 = 1223;
static doublereal c_b1075 = 120.8;
static integer c__1225 = 1225;
static doublereal c_b1086 = 116.8;
static integer c__1228 = 1228;
static integer c__1229 = 1229;
static doublereal c_b1106 = 119.4;



/*     ############################################################# */
/*     ##                  COPYRIGHT (C) 1999 by                  ## */
/*     ##  Marina A. Vorobieva, Nina N. Sokolova & Jay W. Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  program nucleic  --  build a nucleic acid from sequence  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "nucleic" builds the internal and Cartesian coordinates */
/*     of a polynucleotide from nucleic acid sequence and torsional */
/*     angle values for the nucleic acid backbone and side chains */


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
    extern /* Subroutine */ int basefile_(char *, ftnlen), nucchain_(void), 
	    molecule_(void);
    extern integer freeunit_(void), trimtext_(char *, ftnlen);
    static integer i__, mode, iseq, izmt, ixyz;
    extern /* Subroutine */ int field_(void);
    static integer natom;
    static logical exist;
    extern /* Subroutine */ int delete_(integer *), getkey_(void), watson_(
	    void), prtseq_(integer *), prtint_(integer *), prtxyz_(integer *);
    static char seqfile[120];
    extern /* Subroutine */ int initial_(void);
    static char intfile[120];
    extern /* Subroutine */ int connect_(void), inertia_(integer *), makeint_(
	    integer *), getseqn_(void), nextarg_(char *, logical *, ftnlen), 
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




/*     ############################################################# */
/*     ##                  COPYRIGHT (C) 1999 by                  ## */
/*     ##  Marina A. Vorobieva, Nina N. Sokolova & Jay W. Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  nucleo.i  --  parameters for nucleic acid structure  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     bkbone    phosphate backbone angles for each nucleotide */
/*     glyco     glycosidic torsional angle for each nucleotide */
/*     pucker    sugar pucker, either 2=2'-endo or 3=3'-endo */
/*     dblhlx    flag to mark system as nucleic acid double helix */
/*     deoxy     flag to mark deoxyribose or ribose sugar units */
/*     hlxform   helix form (A, B or Z) of polynucleotide strands */




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

    getseqn_();
    nucchain_();
    connect_();
    molecule_();
    makexyz_();

/*     perform the alignment of the strands of a double helix */

    if (nucleo_1.dblhlx) {
	watson_();
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

/*     write out a nucleic acid sequence file */

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
    return 0;
} /* MAIN__ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine getseqn  --  nucleic acid sequence and angles  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "getseqn" asks the user for the nucleotide sequence and */
/*     torsional angle values needed to define a nucleic acid */


/* Subroutine */ int getseqn_(void)
{
    /* Initialized data */

    static char ucase[1*26] = "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" 
	    "L" "M" "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X" "Y" "Z";

    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter A-, B- or Z-Form Helix for the Str"
	    "ucture\002,\002 [B] :  \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_30[] = "(/,\002 Enter One Nucleotide per Line, 5' to 3"
	    "': \002,\002 Give PDB Residue Code,\002,/,\002 followed by Backb"
	    "one Torsions (6F) and\002,\002 Glycosidic Torsion (1F)\002,//"
	    ",\002 Use Residue=MOL to Begin a New Strand,\002,\002 Residue=<C"
	    "R> to End Entry\002)";
    static char fmt_40[] = "(/,\002 Enter Residue\002,i4,\002 :  \002,$)";
    static char fmt_50[] = "(a120)";
    static char fmt_70[] = "(/,\002 GETSEQN  --  Nucleotide Type \002,a3,"
	    "\002 is Not Supported\002)";
    static char fmt_80[] = "(/,\002 Build a Double Helix using Complimentary"
	    " Bases\002,\002 [N] :  \002,$)";
    static char fmt_90[] = "(a120)";
    static char fmt_100[] = "(/,\002 Combine the Two Single Strands into Dou"
	    "ble Helix\002,\002 [Y] :  \002,$)";
    static char fmt_110[] = "(a120)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void), s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer i__, j, k;
    static char name__[3];
    static logical done;
    static integer next, stop, start;
    static char record[120];
    static integer length;
    extern /* Subroutine */ int upcase_(char *, ftnlen);
    static logical purine[10000];
    static char answer[1], string[120], resname[3];
    extern /* Subroutine */ int getword_(char *, char *, integer *, ftnlen, 
	    ftnlen), gettext_(char *, char *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___16 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_50, 0 };
    static icilist io___30 = { 1, string, 1, 0, 120, 1 };
    static cilist io___31 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_110, 0 };



#define seq_ref(a_0,a_1) &sequen_1.seq[(a_1)*3 + a_0 - 3]
#define nuclz_ref(a_0,a_1) &resdue_1.nuclz[(a_1)*3 + a_0 - 3]
#define ichain_ref(a_1,a_2) sequen_1.ichain[(a_2)*2 + a_1 - 3]
#define bkbone_ref(a_1,a_2) nucleo_1.bkbone[(a_2)*6 + a_1 - 7]



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




/*     ############################################################# */
/*     ##                  COPYRIGHT (C) 1999 by                  ## */
/*     ##  Marina A. Vorobieva, Nina N. Sokolova & Jay W. Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  nucleo.i  --  parameters for nucleic acid structure  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     bkbone    phosphate backbone angles for each nucleotide */
/*     glyco     glycosidic torsional angle for each nucleotide */
/*     pucker    sugar pucker, either 2=2'-endo or 3=3'-endo */
/*     dblhlx    flag to mark system as nucleic acid double helix */
/*     deoxy     flag to mark deoxyribose or ribose sugar units */
/*     hlxform   helix form (A, B or Z) of polynucleotide strands */




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




/*     choose to generate either an A-, B- or Z-form helix */

    io___16.ciunit = iounit_1.iout;
    s_wsfe(&io___16);
    e_wsfe();
    io___17.ciunit = iounit_1.input;
    s_rsfe(&io___17);
    do_fio(&c__1, record, (ftnlen)120);
    e_rsfe();
    upcase_(record, (ftnlen)120);
    next = 1;
    getword_(record, answer, &next, (ftnlen)120, (ftnlen)1);
    *(unsigned char *)nucleo_1.hlxform = 'B';
    if (*(unsigned char *)answer == 'A') {
	*(unsigned char *)nucleo_1.hlxform = 'A';
    }
    if (*(unsigned char *)answer == 'Z') {
	*(unsigned char *)nucleo_1.hlxform = 'Z';
    }

/*     provide a header to explain the method of sequence input */

    io___21.ciunit = iounit_1.iout;
    s_wsfe(&io___21);
    e_wsfe();

/*     initially, assume that only a single strand is present */

    sequen_1.nchain = 1;
    ichain_ref(1, 1) = 1;
    *(unsigned char *)&sequen_1.chnnam[0] = ' ';

/*     get the nucleotide sequence data and dihedral angle values */

    i__ = 0;
    done = FALSE_;
    while(! done) {
	++i__;
	for (j = 1; j <= 6; ++j) {
	    bkbone_ref(j, i__) = 0.;
	}
	nucleo_1.glyco[i__ - 1] = 0.;
	nucleo_1.pucker[i__ - 1] = 0;
	io___25.ciunit = iounit_1.iout;
	s_wsfe(&io___25);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
	io___26.ciunit = iounit_1.input;
	s_rsfe(&io___26);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	upcase_(record, (ftnlen)120);
	next = 1;
	getword_(record, name__, &next, (ftnlen)120, (ftnlen)3);
	length = trimtext_(name__, (ftnlen)3);
	s_copy(string, record + (next - 1), (ftnlen)120, 120 - (next - 1));
	i__1 = s_rsli(&io___30);
	if (i__1 != 0) {
	    goto L60;
	}
	for (j = 1; j <= 6; ++j) {
	    i__1 = do_lio(&c__5, &c__1, (char *)&bkbone_ref(j, i__), (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L60;
	    }
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&nucleo_1.glyco[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L60;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L60;
	}
L60:

/*     process and store the current nucleotide type */

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
		s_copy(seq_ref(0, i__), nuclz_ref(0, 12), (ftnlen)3, (ftnlen)
			3);
		sequen_1.seqtyp[i__ - 1] = 0;
		if (length == 1) {
		    for (j = 1; j <= 12; ++j) {
			if (*(unsigned char *)name__ == *(unsigned char *)&
				resdue_1.nuclz1[j - 1]) {
			    s_copy(seq_ref(0, i__), nuclz_ref(0, j), (ftnlen)
				    3, (ftnlen)3);
			    sequen_1.seqtyp[i__ - 1] = j;
			}
		    }
		} else {
		    for (j = 1; j <= 12; ++j) {
			if (s_cmp(name__, nuclz_ref(0, j), (ftnlen)3, (ftnlen)
				3) == 0) {
			    s_copy(seq_ref(0, i__), nuclz_ref(0, j), (ftnlen)
				    3, (ftnlen)3);
			    sequen_1.seqtyp[i__ - 1] = j;
			}
		    }
		}
		if (sequen_1.seqtyp[i__ - 1] == 0) {
		    --i__;
		    io___31.ciunit = iounit_1.iout;
		    s_wsfe(&io___31);
		    do_fio(&c__1, name__, (ftnlen)3);
		    e_wsfe();
		}
	    }
	}
    }

/*     offer the option to construct an idealized double helix */

    nucleo_1.dblhlx = FALSE_;
    if (sequen_1.nchain == 1) {
	io___32.ciunit = iounit_1.iout;
	s_wsfe(&io___32);
	e_wsfe();
	io___33.ciunit = iounit_1.input;
	s_rsfe(&io___33);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
	upcase_(answer, (ftnlen)1);
	if (*(unsigned char *)answer == 'Y') {
	    nucleo_1.dblhlx = TRUE_;
	}
    } else if (sequen_1.nchain == 2) {
	io___34.ciunit = iounit_1.iout;
	s_wsfe(&io___34);
	e_wsfe();
	io___35.ciunit = iounit_1.input;
	s_rsfe(&io___35);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	next = 1;
	gettext_(record, answer, &next, (ftnlen)120, (ftnlen)1);
	upcase_(answer, (ftnlen)1);
	if (*(unsigned char *)answer != 'N') {
	    nucleo_1.dblhlx = TRUE_;
	}
    }

/*     build a second strand as the reverse-compliment sequence */

    if (sequen_1.nchain == 1 && nucleo_1.dblhlx) {
	start = 1;
	stop = sequen_1.nseq;
	s_copy(resname, nuclz_ref(0, sequen_1.seqtyp[0]), (ftnlen)3, (ftnlen)
		3);
	if (s_cmp(resname, "MP ", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(resname,
		 "DP ", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(resname, "TP ", (
		ftnlen)3, (ftnlen)3) == 0) {
	    k = sequen_1.nseq + 1;
	    s_copy(seq_ref(0, k), seq_ref(0, 1), (ftnlen)3, (ftnlen)3);
	    sequen_1.seqtyp[k - 1] = sequen_1.seqtyp[0];
	    start = 2;
	}
	s_copy(resname, nuclz_ref(0, sequen_1.seqtyp[sequen_1.nseq - 1]), (
		ftnlen)3, (ftnlen)3);
	if (s_cmp(resname, "MP ", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(resname,
		 "DP ", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(resname, "TP ", (
		ftnlen)3, (ftnlen)3) == 0) {
	    k = sequen_1.nseq << 1;
	    s_copy(seq_ref(0, k), seq_ref(0, sequen_1.nseq), (ftnlen)3, (
		    ftnlen)3);
	    sequen_1.seqtyp[k - 1] = sequen_1.seqtyp[sequen_1.nseq - 1];
	    stop = sequen_1.nseq - 1;
	}
	i__1 = stop;
	for (i__ = start; i__ <= i__1; ++i__) {
	    s_copy(resname, nuclz_ref(0, sequen_1.seqtyp[i__ - 1]), (ftnlen)3,
		     (ftnlen)3);
	    if (s_cmp(resname, "A  ", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(resname, "U  ", (ftnlen)3, (ftnlen)3);
	    } else if (s_cmp(resname, "G  ", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(resname, "C  ", (ftnlen)3, (ftnlen)3);
	    } else if (s_cmp(resname, "C  ", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(resname, "G  ", (ftnlen)3, (ftnlen)3);
	    } else if (s_cmp(resname, "U  ", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(resname, "A  ", (ftnlen)3, (ftnlen)3);
	    } else if (s_cmp(resname, "DA ", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(resname, "DT ", (ftnlen)3, (ftnlen)3);
	    } else if (s_cmp(resname, "DG ", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(resname, "DC ", (ftnlen)3, (ftnlen)3);
	    } else if (s_cmp(resname, "DC ", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(resname, "DG ", (ftnlen)3, (ftnlen)3);
	    } else if (s_cmp(resname, "DT ", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(resname, "DA ", (ftnlen)3, (ftnlen)3);
	    }
	    k = sequen_1.nseq + stop + start - i__;
	    for (j = 1; j <= 12; ++j) {
		if (s_cmp(resname, nuclz_ref(0, j), (ftnlen)3, (ftnlen)3) == 
			0) {
		    s_copy(seq_ref(0, k), nuclz_ref(0, j), (ftnlen)3, (ftnlen)
			    3);
		    sequen_1.seqtyp[k - 1] = j;
		}
	    }
	}
	i__1 = sequen_1.nseq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = sequen_1.nseq + i__;
	    for (j = 1; j <= 6; ++j) {
		bkbone_ref(j, k) = bkbone_ref(j, i__);
	    }
	    nucleo_1.glyco[k - 1] = nucleo_1.glyco[i__ - 1];
	    nucleo_1.pucker[k - 1] = nucleo_1.pucker[i__ - 1];
	}
	sequen_1.nchain = 2;
	sequen_1.nseq <<= 1;
	ichain_ref(1, sequen_1.nchain) = sequen_1.nseq / 2 + 1;
	ichain_ref(2, sequen_1.nchain) = sequen_1.nseq;
    }

/*     set chain identifiers if multiple chains are present */

    if (sequen_1.nchain > 1) {
	i__1 = sequen_1.nchain;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    *(unsigned char *)&sequen_1.chnnam[i__ - 1] = *(unsigned char *)&
		    ucase[i__ - 1];
	}
    }

/*     set the nucleic acid base and sugar structural type */

    i__1 = sequen_1.nseq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(resname, nuclz_ref(0, sequen_1.seqtyp[i__ - 1]), (ftnlen)3, (
		ftnlen)3);
	purine[i__ - 1] = FALSE_;
	if (s_cmp(resname, "A  ", (ftnlen)3, (ftnlen)3) == 0) {
	    purine[i__ - 1] = TRUE_;
	}
	if (s_cmp(resname, "G  ", (ftnlen)3, (ftnlen)3) == 0) {
	    purine[i__ - 1] = TRUE_;
	}
	if (s_cmp(resname, "DA ", (ftnlen)3, (ftnlen)3) == 0) {
	    purine[i__ - 1] = TRUE_;
	}
	if (s_cmp(resname, "DG ", (ftnlen)3, (ftnlen)3) == 0) {
	    purine[i__ - 1] = TRUE_;
	}
	nucleo_1.deoxy[i__ - 1] = FALSE_;
	if (s_cmp(resname, "DA ", (ftnlen)3, (ftnlen)3) == 0) {
	    nucleo_1.deoxy[i__ - 1] = TRUE_;
	}
	if (s_cmp(resname, "DG ", (ftnlen)3, (ftnlen)3) == 0) {
	    nucleo_1.deoxy[i__ - 1] = TRUE_;
	}
	if (s_cmp(resname, "DC ", (ftnlen)3, (ftnlen)3) == 0) {
	    nucleo_1.deoxy[i__ - 1] = TRUE_;
	}
	if (s_cmp(resname, "DT ", (ftnlen)3, (ftnlen)3) == 0) {
	    nucleo_1.deoxy[i__ - 1] = TRUE_;
	}
    }

/*     set the backbone and glycosidic torsions and sugar pucker */

    i__1 = sequen_1.nseq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	done = FALSE_;
	for (j = 1; j <= 6; ++j) {
	    if (bkbone_ref(j, i__) != 0.) {
		done = TRUE_;
	    }
	}
	if (nucleo_1.glyco[i__ - 1] != 0.) {
	    done = TRUE_;
	}
	if (nucleo_1.pucker[i__ - 1] != 0) {
	    done = TRUE_;
	}
	if (! done) {
	    if (*(unsigned char *)nucleo_1.hlxform == 'A') {
		bkbone_ref(1, i__) = -51.8;
		bkbone_ref(2, i__) = 174.8;
		bkbone_ref(3, i__) = 41.7;
		bkbone_ref(4, i__) = 79.1;
		bkbone_ref(5, i__) = -148.;
		bkbone_ref(6, i__) = -75.;
		nucleo_1.glyco[i__ - 1] = -157.2;
		nucleo_1.pucker[i__ - 1] = 3;
	    } else if (*(unsigned char *)nucleo_1.hlxform == 'B') {
		bkbone_ref(1, i__) = -46.1;
		bkbone_ref(2, i__) = -146.5;
		bkbone_ref(3, i__) = 36.4;
		bkbone_ref(4, i__) = 156.5;
		bkbone_ref(5, i__) = 154.7;
		bkbone_ref(6, i__) = -95.6;
		nucleo_1.glyco[i__ - 1] = -97.8;
		nucleo_1.pucker[i__ - 1] = 2;
	    } else if (*(unsigned char *)nucleo_1.hlxform == 'Z') {
		if (purine[i__ - 1]) {
		    bkbone_ref(1, i__) = 48.;
		    bkbone_ref(2, i__) = 179.;
		    bkbone_ref(3, i__) = -170.;
		    bkbone_ref(4, i__) = 100.;
		    bkbone_ref(5, i__) = -104.;
		    bkbone_ref(6, i__) = -69.;
		    nucleo_1.glyco[i__ - 1] = 67.;
		    nucleo_1.pucker[i__ - 1] = 3;
		} else {
		    bkbone_ref(1, i__) = -137.;
		    bkbone_ref(2, i__) = -139.;
		    bkbone_ref(3, i__) = 55.;
		    bkbone_ref(4, i__) = 138.;
		    bkbone_ref(5, i__) = -94.;
		    bkbone_ref(6, i__) = 80.;
		    nucleo_1.glyco[i__ - 1] = -159.;
		    nucleo_1.pucker[i__ - 1] = 1;
		}
	    }
	}
    }
    return 0;
} /* getseqn_ */

#undef bkbone_ref
#undef ichain_ref
#undef nuclz_ref
#undef seq_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine nucchain  --  build polynucleotide backbone  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "nucchain" builds up the internal coordinates for a nucleic */
/*     acid sequence from the sugar type, backbone and glycosidic */
/*     torsional values */


/* Subroutine */ int nucchain_(void)
{
    /* Initialized data */

    static integer o5typ[12] = { 1001,1031,1062,1090,1117,1146,1176,1203,0,0,
	    0,0 };

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, k, m, c1i, c2i, c3i, c4i, c5i, o2i, o3i, o4i, o5i, 
	    poi;
    static logical cap3, cap5;
    static integer ptyp;
    extern /* Subroutine */ int zatom_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *);
    static integer httyp, optyp, ostyp, ottyp;
    static logical single;
    extern /* Subroutine */ int nucbase_(char *, integer *, integer *, 
	    integer *, integer *, ftnlen);
    static char resname[3];


#define nuclz_ref(a_0,a_1) &resdue_1.nuclz[(a_1)*3 + a_0 - 3]
#define ichain_ref(a_1,a_2) sequen_1.ichain[(a_2)*2 + a_1 - 3]
#define bkbone_ref(a_1,a_2) nucleo_1.bkbone[(a_2)*6 + a_1 - 7]



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




/*     ############################################################# */
/*     ##                  COPYRIGHT (C) 1999 by                  ## */
/*     ##  Marina A. Vorobieva, Nina N. Sokolova & Jay W. Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  nucleo.i  --  parameters for nucleic acid structure  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     bkbone    phosphate backbone angles for each nucleotide */
/*     glyco     glycosidic torsional angle for each nucleotide */
/*     pucker    sugar pucker, either 2=2'-endo or 3=3'-endo */
/*     dblhlx    flag to mark system as nucleic acid double helix */
/*     deoxy     flag to mark deoxyribose or ribose sugar units */
/*     hlxform   helix form (A, B or Z) of polynucleotide strands */




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



/*     biopolymer atom types for nucleotide backbone atoms */



/*     initialize the atom counter to the first atom */

    atoms_1.n = 1;

/*     check for single residue and 3'- or 5'-phosphate caps */

    i__1 = sequen_1.nchain;
    for (m = 1; m <= i__1; ++m) {
	single = FALSE_;
	cap5 = FALSE_;
	cap3 = FALSE_;
	if (ichain_ref(1, m) == ichain_ref(2, m)) {
	    single = TRUE_;
	}
	i__ = ichain_ref(1, m);
	k = sequen_1.seqtyp[i__ - 1];
	s_copy(resname, nuclz_ref(0, k), (ftnlen)3, (ftnlen)3);
	if (s_cmp(resname, "MP ", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(resname,
		 "DP ", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(resname, "TP ", (
		ftnlen)3, (ftnlen)3) == 0) {
	    cap5 = TRUE_;
	}
	i__ = ichain_ref(2, m);
	k = sequen_1.seqtyp[i__ - 1];
	s_copy(resname, nuclz_ref(0, k), (ftnlen)3, (ftnlen)3);
	if (s_cmp(resname, "MP ", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(resname,
		 "DP ", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(resname, "TP ", (
		ftnlen)3, (ftnlen)3) == 0) {
	    cap3 = TRUE_;
	}

/*     build the first residue or a phosphate capping group */

	i__ = ichain_ref(1, m);
	k = sequen_1.seqtyp[i__ - 1];
	j = o5typ[k - 1];
	s_copy(resname, nuclz_ref(0, k), (ftnlen)3, (ftnlen)3);
	if (s_cmp(resname, "MP ", (ftnlen)3, (ftnlen)3) == 0) {
	    if (nucleo_1.deoxy[i__]) {
		ostyp = 1246;
		ptyp = 1247;
		optyp = 1248;
	    } else {
		ostyp = 1234;
		ptyp = 1235;
		optyp = 1236;
	    }
	    if (m == 1) {
		o3i = atoms_1.n;
		zatom_(&optyp, &c_b108, &c_b108, &c_b108, &c__0, &c__0, &c__0,
			 &c__0);
		poi = atoms_1.n;
		zatom_(&ptyp, &c_b115, &c_b108, &c_b108, &o3i, &c__0, &c__0, &
			c__0);
		zatom_(&optyp, &c_b115, &c_b122, &c_b108, &poi, &o3i, &c__0, &
			c__0);
	    } else {
		o3i = atoms_1.n;
		i__2 = atoms_1.n - 1;
		i__3 = atoms_1.n - 2;
		i__4 = atoms_1.n - 3;
		zatom_(&optyp, &c_b126, &c_b127, &c_b128, &i__2, &i__3, &i__4,
			 &c__0);
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 1;
		zatom_(&c_n2, &c_b108, &c_b108, &c_b108, &i__2, &i__3, &c__0, 
			&c__0);
		poi = atoms_1.n;
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 3;
		zatom_(&ptyp, &c_b115, &c_b127, &c_b128, &o3i, &i__2, &i__3, &
			c__0);
		i__2 = atoms_1.n - 3;
		zatom_(&optyp, &c_b115, &c_b122, &c_b128, &poi, &o3i, &i__2, &
			c__0);
	    }
	    i__2 = atoms_1.n - 1;
	    zatom_(&optyp, &c_b115, &c_b122, &c_b122, &poi, &o3i, &i__2, &
		    c__1);
	    o5i = atoms_1.n;
	    i__2 = atoms_1.n - 2;
	    zatom_(&ostyp, &c_b148, &c_b149, &c_b149, &poi, &o3i, &i__2, &
		    c_n1);
	} else if (s_cmp(resname, "DP ", (ftnlen)3, (ftnlen)3) == 0) {
	} else if (s_cmp(resname, "TP ", (ftnlen)3, (ftnlen)3) == 0) {
	} else {
	    if (nucleo_1.deoxy[i__ - 1]) {
		ottyp = 1244;
	    } else {
		ottyp = 1232;
	    }
	    if (m == 1) {
		o5i = atoms_1.n;
		zatom_(&ottyp, &c_b108, &c_b108, &c_b108, &c__0, &c__0, &c__0,
			 &c__0);
		c5i = atoms_1.n;
		i__2 = j + 1;
		zatom_(&i__2, &c_b161, &c_b108, &c_b108, &o5i, &c__0, &c__0, &
			c__0);
		c4i = atoms_1.n;
		i__2 = j + 4;
		zatom_(&i__2, &c_b115, &c_b168, &c_b108, &c5i, &o5i, &c__0, &
			c__0);
	    } else {
		o5i = atoms_1.n;
		i__2 = atoms_1.n - 1;
		i__3 = atoms_1.n - 2;
		i__4 = atoms_1.n - 3;
		zatom_(&ottyp, &c_b172, &c_b127, &c_b128, &i__2, &i__3, &i__4,
			 &c__0);
		i__2 = atoms_1.n - 2;
		i__3 = atoms_1.n - 1;
		zatom_(&c_n2, &c_b108, &c_b108, &c_b108, &i__2, &i__3, &c__0, 
			&c__0);
		c5i = atoms_1.n;
		i__2 = j + 1;
		i__3 = atoms_1.n - 2;
		i__4 = atoms_1.n - 3;
		zatom_(&i__2, &c_b161, &c_b183, &c_b128, &o5i, &i__3, &i__4, &
			c__0);
		c4i = atoms_1.n;
		i__2 = j + 4;
		i__3 = atoms_1.n - 3;
		zatom_(&i__2, &c_b115, &c_b168, &c_b128, &c5i, &o5i, &i__3, &
			c__0);
	    }
	    o4i = atoms_1.n;
	    i__2 = j + 6;
	    d__1 = bkbone_ref(3, i__) - 120.;
	    zatom_(&i__2, &c_b190, &c_b191, &d__1, &c4i, &c5i, &o5i, &c__0);
	    c1i = atoms_1.n;
	    if (nucleo_1.pucker[i__ - 1] == 3) {
		i__2 = j + 7;
		zatom_(&i__2, &c_b193, &c_b194, &c_b195, &o4i, &c4i, &c5i, &
			c__0);
	    } else if (nucleo_1.pucker[i__ - 1] == 2) {
		i__2 = j + 7;
		zatom_(&i__2, &c_b193, &c_b194, &c_b199, &o4i, &c4i, &c5i, &
			c__0);
	    } else if (nucleo_1.pucker[i__ - 1] == 1) {
		i__2 = j + 7;
		zatom_(&i__2, &c_b193, &c_b194, &c_b203, &o4i, &c4i, &c5i, &
			c__0);
	    }
	    c3i = atoms_1.n;
	    i__2 = j + 9;
	    zatom_(&i__2, &c_b205, &c_b206, &bkbone_ref(3, i__), &c4i, &c5i, &
		    o5i, &c__0);
	    c2i = atoms_1.n;
	    i__2 = j + 11;
	    d__1 = bkbone_ref(4, i__) + 120.;
	    zatom_(&i__2, &c_b205, &c_b209, &d__1, &c3i, &c4i, &c5i, &c__0);
	    zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &c1i, &c2i, &c__0, &c__0)
		    ;
	    o3i = atoms_1.n;
	    if (nucleo_1.deoxy[i__ - 1]) {
		if (single) {
		    ottyp = 1249;
		    zatom_(&ottyp, &c_b193, &c_b218, &bkbone_ref(4, i__), &
			    c3i, &c4i, &c5i, &c__0);
		} else {
		    i__2 = j + 14;
		    zatom_(&i__2, &c_b193, &c_b218, &bkbone_ref(4, i__), &c3i,
			     &c4i, &c5i, &c__0);
		}
	    } else {
		if (single) {
		    ottyp = 1237;
		    zatom_(&ottyp, &c_b193, &c_b218, &bkbone_ref(4, i__), &
			    c3i, &c4i, &c5i, &c__0);
		} else {
		    i__2 = j + 15;
		    zatom_(&i__2, &c_b193, &c_b218, &bkbone_ref(4, i__), &c3i,
			     &c4i, &c5i, &c__0);
		}
		o2i = atoms_1.n;
		i__2 = j + 13;
		zatom_(&i__2, &c_b229, &c_b230, &c_b230, &c2i, &c3i, &c1i, &
			c__1);
	    }
	    if (nucleo_1.deoxy[i__ - 1]) {
		httyp = 1245;
	    } else {
		httyp = 1233;
	    }
	    zatom_(&httyp, &c_b172, &c_b199, &c_b128, &o5i, &c5i, &c4i, &c__0)
		    ;
	    i__2 = j + 2;
	    zatom_(&i__2, &c_b237, &c_b230, &c_b230, &c5i, &o5i, &c4i, &c__1);
	    i__2 = j + 3;
	    zatom_(&i__2, &c_b237, &c_b230, &c_b230, &c5i, &o5i, &c4i, &c_n1);
	    i__2 = j + 5;
	    zatom_(&i__2, &c_b237, &c_b230, &c_b230, &c4i, &c5i, &c3i, &c_n1);
	    if (nucleo_1.pucker[i__ - 1] == 3) {
		i__2 = j + 8;
		zatom_(&i__2, &c_b237, &c_b230, &c_b251, &c1i, &o4i, &c2i, &
			c_n1);
	    } else if (nucleo_1.pucker[i__ - 1] == 2) {
		i__2 = j + 8;
		zatom_(&i__2, &c_b237, &c_b230, &c_b255, &c1i, &o4i, &c2i, &
			c_n1);
	    } else if (nucleo_1.pucker[i__ - 1] == 1) {
		i__2 = j + 8;
		zatom_(&i__2, &c_b237, &c_b230, &c_b259, &c1i, &o4i, &c2i, &
			c_n1);
	    }
	    i__2 = j + 10;
	    zatom_(&i__2, &c_b237, &c_b230, &c_b230, &c3i, &c4i, &c2i, &c_n1);
	    i__2 = j + 12;
	    zatom_(&i__2, &c_b237, &c_b230, &c_b230, &c2i, &c3i, &c1i, &c_n1);
	    if (nucleo_1.deoxy[i__ - 1]) {
		i__2 = j + 13;
		zatom_(&i__2, &c_b237, &c_b230, &c_b230, &c2i, &c3i, &c1i, &
			c__1);
	    } else {
		i__2 = j + 14;
		zatom_(&i__2, &c_b172, &c_b199, &c_b128, &o2i, &c2i, &c3i, &
			c__0);
	    }
	    if (single) {
		if (nucleo_1.deoxy[i__ - 1]) {
		    httyp = 1250;
		} else {
		    httyp = 1238;
		}
		zatom_(&httyp, &c_b172, &c_b255, &c_b128, &o3i, &c3i, &c4i, &
			c__0);
	    }
	    nucbase_(resname, &i__, &c1i, &o4i, &c2i, (ftnlen)3);
	}

/*     build atoms for residues in the middle of the chain */

	i__2 = ichain_ref(2, m) - 1;
	for (i__ = ichain_ref(1, m) + 1; i__ <= i__2; ++i__) {
	    k = sequen_1.seqtyp[i__ - 1];
	    j = o5typ[k - 1];
	    s_copy(resname, nuclz_ref(0, k), (ftnlen)3, (ftnlen)3);
	    if (cap5) {
		cap5 = FALSE_;
	    } else {
		if (nucleo_1.deoxy[i__ - 1]) {
		    ptyp = 1242;
		    optyp = 1243;
		} else {
		    ptyp = 1230;
		    optyp = 1231;
		}
		poi = atoms_1.n;
		zatom_(&ptyp, &c_b281, &c_b183, &bkbone_ref(5, i__ - 1), &o3i,
			 &c3i, &c4i, &c__0);
		d__1 = bkbone_ref(6, i__ - 1) + 120.;
		zatom_(&optyp, &c_b284, &c_b285, &d__1, &poi, &o3i, &c3i, &
			c__0);
		d__1 = bkbone_ref(6, i__ - 1) - 120.;
		zatom_(&optyp, &c_b284, &c_b285, &d__1, &poi, &o3i, &c3i, &
			c__0);
		o5i = atoms_1.n;
		zatom_(&j, &c_b281, &c_b291, &bkbone_ref(6, i__ - 1), &poi, &
			o3i, &c3i, &c__0);
	    }
	    c5i = atoms_1.n;
	    i__3 = j + 1;
	    zatom_(&i__3, &c_b161, &c_b183, &bkbone_ref(1, i__), &o5i, &poi, &
		    o3i, &c__0);
	    c4i = atoms_1.n;
	    i__3 = j + 4;
	    zatom_(&i__3, &c_b115, &c_b168, &bkbone_ref(2, i__), &c5i, &o5i, &
		    poi, &c__0);
	    o4i = atoms_1.n;
	    i__3 = j + 6;
	    d__1 = bkbone_ref(3, i__) - 120.;
	    zatom_(&i__3, &c_b190, &c_b191, &d__1, &c4i, &c5i, &o5i, &c__0);
	    c1i = atoms_1.n;
	    if (nucleo_1.pucker[i__ - 1] == 3) {
		i__3 = j + 7;
		zatom_(&i__3, &c_b193, &c_b194, &c_b195, &o4i, &c4i, &c5i, &
			c__0);
	    } else if (nucleo_1.pucker[i__ - 1] == 2) {
		i__3 = j + 7;
		zatom_(&i__3, &c_b193, &c_b194, &c_b199, &o4i, &c4i, &c5i, &
			c__0);
	    } else if (nucleo_1.pucker[i__ - 1] == 1) {
		i__3 = j + 7;
		zatom_(&i__3, &c_b193, &c_b194, &c_b203, &o4i, &c4i, &c5i, &
			c__0);
	    }
	    c3i = atoms_1.n;
	    i__3 = j + 9;
	    zatom_(&i__3, &c_b205, &c_b206, &bkbone_ref(3, i__), &c4i, &c5i, &
		    o5i, &c__0);
	    c2i = atoms_1.n;
	    i__3 = j + 11;
	    d__1 = bkbone_ref(4, i__) + 120.;
	    zatom_(&i__3, &c_b205, &c_b209, &d__1, &c3i, &c4i, &c5i, &c__0);
	    zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &c1i, &c2i, &c__0, &c__0)
		    ;
	    o3i = atoms_1.n;
	    if (nucleo_1.deoxy[i__ - 1]) {
		if (cap3) {
		    zatom_(&c__1251, &c_b193, &c_b218, &bkbone_ref(4, i__), &
			    c3i, &c4i, &c5i, &c__0);
		} else {
		    i__3 = j + 14;
		    zatom_(&i__3, &c_b193, &c_b218, &bkbone_ref(4, i__), &c3i,
			     &c4i, &c5i, &c__0);
		}
	    } else {
		if (cap3) {
		    zatom_(&c__1239, &c_b193, &c_b218, &bkbone_ref(4, i__), &
			    c3i, &c4i, &c5i, &c__0);
		} else {
		    i__3 = j + 15;
		    zatom_(&i__3, &c_b193, &c_b218, &bkbone_ref(4, i__), &c3i,
			     &c4i, &c5i, &c__0);
		}
		o2i = atoms_1.n;
		i__3 = j + 13;
		zatom_(&i__3, &c_b229, &c_b230, &c_b230, &c2i, &c3i, &c1i, &
			c__1);
	    }
	    i__3 = j + 2;
	    zatom_(&i__3, &c_b237, &c_b230, &c_b230, &c5i, &o5i, &c4i, &c__1);
	    i__3 = j + 3;
	    zatom_(&i__3, &c_b237, &c_b230, &c_b230, &c5i, &o5i, &c4i, &c_n1);
	    i__3 = j + 5;
	    zatom_(&i__3, &c_b237, &c_b230, &c_b230, &c4i, &c5i, &c3i, &c_n1);
	    if (nucleo_1.pucker[i__ - 1] == 3) {
		i__3 = j + 8;
		zatom_(&i__3, &c_b237, &c_b230, &c_b251, &c1i, &o4i, &c2i, &
			c_n1);
	    } else if (nucleo_1.pucker[i__ - 1] == 2) {
		i__3 = j + 8;
		zatom_(&i__3, &c_b237, &c_b230, &c_b255, &c1i, &o4i, &c2i, &
			c_n1);
	    } else if (nucleo_1.pucker[i__ - 1] == 1) {
		i__3 = j + 8;
		zatom_(&i__3, &c_b237, &c_b230, &c_b259, &c1i, &o4i, &c2i, &
			c_n1);
	    }
	    i__3 = j + 10;
	    zatom_(&i__3, &c_b237, &c_b230, &c_b230, &c3i, &c4i, &c2i, &c_n1);
	    i__3 = j + 12;
	    zatom_(&i__3, &c_b237, &c_b230, &c_b230, &c2i, &c3i, &c1i, &c_n1);
	    if (nucleo_1.deoxy[i__ - 1]) {
		i__3 = j + 13;
		zatom_(&i__3, &c_b237, &c_b230, &c_b230, &c2i, &c3i, &c1i, &
			c__1);
	    } else {
		i__3 = j + 14;
		zatom_(&i__3, &c_b172, &c_b199, &c_b128, &o2i, &c2i, &c3i, &
			c__0);
	    }
	    nucbase_(resname, &i__, &c1i, &o4i, &c2i, (ftnlen)3);
	}

/*     build the last residue or a phosphate capping group */

	i__ = ichain_ref(2, m);
	k = sequen_1.seqtyp[i__ - 1];
	j = o5typ[k - 1];
	s_copy(resname, nuclz_ref(0, k), (ftnlen)3, (ftnlen)3);
	if (single) {
	} else if (s_cmp(resname, "MP ", (ftnlen)3, (ftnlen)3) == 0) {
	    if (nucleo_1.deoxy[i__ - 2]) {
		ptyp = 1252;
		optyp = 1253;
	    } else {
		ptyp = 1240;
		optyp = 1241;
	    }
	    poi = atoms_1.n;
	    zatom_(&ptyp, &c_b148, &c_b183, &bkbone_ref(5, i__ - 1), &o3i, &
		    c3i, &c4i, &c__0);
	    zatom_(&optyp, &c_b115, &c_b149, &c_b390, &poi, &o3i, &c3i, &c__0)
		    ;
	    zatom_(&optyp, &c_b115, &c_b149, &c_b394, &poi, &o3i, &c3i, &c__0)
		    ;
	    zatom_(&optyp, &c_b115, &c_b149, &c_b128, &poi, &o3i, &c3i, &c__0)
		    ;
	} else if (s_cmp(resname, "DP ", (ftnlen)3, (ftnlen)3) == 0) {
	} else if (s_cmp(resname, "TP ", (ftnlen)3, (ftnlen)3) == 0) {
	} else {
	    if (cap5) {
		cap5 = FALSE_;
	    } else {
		if (nucleo_1.deoxy[i__ - 1]) {
		    ptyp = 1242;
		    optyp = 1243;
		} else {
		    ptyp = 1230;
		    optyp = 1231;
		}
		poi = atoms_1.n;
		zatom_(&ptyp, &c_b281, &c_b183, &bkbone_ref(5, i__ - 1), &o3i,
			 &c3i, &c4i, &c__0);
		d__1 = bkbone_ref(6, i__ - 1) + 120.;
		zatom_(&optyp, &c_b284, &c_b285, &d__1, &poi, &o3i, &c3i, &
			c__0);
		d__1 = bkbone_ref(6, i__ - 1) - 120.;
		zatom_(&optyp, &c_b284, &c_b285, &d__1, &poi, &o3i, &c3i, &
			c__0);
		o5i = atoms_1.n;
		zatom_(&j, &c_b281, &c_b291, &bkbone_ref(6, i__ - 1), &poi, &
			o3i, &c3i, &c__0);
	    }
	    c5i = atoms_1.n;
	    i__2 = j + 1;
	    zatom_(&i__2, &c_b161, &c_b183, &bkbone_ref(1, i__), &o5i, &poi, &
		    o3i, &c__0);
	    c4i = atoms_1.n;
	    i__2 = j + 4;
	    zatom_(&i__2, &c_b115, &c_b168, &bkbone_ref(2, i__), &c5i, &o5i, &
		    poi, &c__0);
	    o4i = atoms_1.n;
	    i__2 = j + 6;
	    d__1 = bkbone_ref(3, i__) - 120.;
	    zatom_(&i__2, &c_b190, &c_b191, &d__1, &c4i, &c5i, &o5i, &c__0);
	    c1i = atoms_1.n;
	    if (nucleo_1.pucker[i__ - 1] == 3) {
		i__2 = j + 7;
		zatom_(&i__2, &c_b193, &c_b194, &c_b195, &o4i, &c4i, &c5i, &
			c__0);
	    } else if (nucleo_1.pucker[i__ - 1] == 2) {
		i__2 = j + 7;
		zatom_(&i__2, &c_b193, &c_b194, &c_b199, &o4i, &c4i, &c5i, &
			c__0);
	    } else if (nucleo_1.pucker[i__ - 1] == 1) {
		i__2 = j + 7;
		zatom_(&i__2, &c_b193, &c_b194, &c_b203, &o4i, &c4i, &c5i, &
			c__0);
	    }
	    c3i = atoms_1.n;
	    i__2 = j + 9;
	    zatom_(&i__2, &c_b205, &c_b206, &bkbone_ref(3, i__), &c4i, &c5i, &
		    o5i, &c__0);
	    c2i = atoms_1.n;
	    i__2 = j + 11;
	    d__1 = bkbone_ref(4, i__) + 120.;
	    zatom_(&i__2, &c_b205, &c_b209, &d__1, &c3i, &c4i, &c5i, &c__0);
	    zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &c1i, &c2i, &c__0, &c__0)
		    ;
	    o3i = atoms_1.n;
	    if (nucleo_1.deoxy[i__ - 1]) {
		ottyp = 1249;
		zatom_(&ottyp, &c_b193, &c_b218, &bkbone_ref(4, i__), &c3i, &
			c4i, &c5i, &c__0);
	    } else {
		ottyp = 1237;
		zatom_(&ottyp, &c_b193, &c_b218, &bkbone_ref(4, i__), &c3i, &
			c4i, &c5i, &c__0);
		o2i = atoms_1.n;
		i__2 = j + 13;
		zatom_(&i__2, &c_b229, &c_b230, &c_b230, &c2i, &c3i, &c1i, &
			c__1);
	    }
	    i__2 = j + 2;
	    zatom_(&i__2, &c_b237, &c_b230, &c_b230, &c5i, &o5i, &c4i, &c__1);
	    i__2 = j + 3;
	    zatom_(&i__2, &c_b237, &c_b230, &c_b230, &c5i, &o5i, &c4i, &c_n1);
	    i__2 = j + 5;
	    zatom_(&i__2, &c_b237, &c_b230, &c_b230, &c4i, &c5i, &c3i, &c_n1);
	    if (nucleo_1.pucker[i__ - 1] == 3) {
		i__2 = j + 8;
		zatom_(&i__2, &c_b237, &c_b230, &c_b251, &c1i, &o4i, &c2i, &
			c_n1);
	    } else if (nucleo_1.pucker[i__ - 1] == 2) {
		i__2 = j + 8;
		zatom_(&i__2, &c_b237, &c_b230, &c_b255, &c1i, &o4i, &c2i, &
			c_n1);
	    } else if (nucleo_1.pucker[i__ - 1] == 1) {
		i__2 = j + 8;
		zatom_(&i__2, &c_b237, &c_b230, &c_b259, &c1i, &o4i, &c2i, &
			c_n1);
	    }
	    i__2 = j + 10;
	    zatom_(&i__2, &c_b237, &c_b230, &c_b230, &c3i, &c4i, &c2i, &c_n1);
	    i__2 = j + 12;
	    zatom_(&i__2, &c_b237, &c_b230, &c_b230, &c2i, &c3i, &c1i, &c_n1);
	    if (nucleo_1.deoxy[i__ - 1]) {
		httyp = 1250;
		i__2 = j + 13;
		zatom_(&i__2, &c_b237, &c_b230, &c_b230, &c2i, &c3i, &c1i, &
			c__1);
	    } else {
		httyp = 1238;
		i__2 = j + 14;
		zatom_(&i__2, &c_b172, &c_b199, &c_b128, &o2i, &c2i, &c3i, &
			c__0);
	    }
	    zatom_(&httyp, &c_b172, &c_b255, &c_b128, &o3i, &c3i, &c4i, &c__0)
		    ;
	    nucbase_(resname, &i__, &c1i, &o4i, &c2i, (ftnlen)3);
	}
    }

/*     finally, set the total number of atoms */

    --atoms_1.n;
    return 0;
} /* nucchain_ */

#undef bkbone_ref
#undef ichain_ref
#undef nuclz_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine nucbase  --  build nucleotide base side chain  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "nucbase" builds the side chain for a single nucleotide base */
/*     in terms of internal coordinates */

/*     resname   3-letter name of current nucleotide residue */
/*     i         number of the current nucleotide residue */
/*     c1i       atom number of carbon C1' in residue i */
/*     o4i       atom number of oxygen O4' in residue i */
/*     c2i       atom number of carbon C2' in residue i */

/*     literature references: */

/*     R. Lavery, K. Zakrzewska, "Base and Base Pair Morphologies, */
/*     Helical Parameters, and Definitions" in "Oxford Handbook of */
/*     Nucleic Acid Structure", S. Neidel, Editor, Oxford University */
/*     Press, 1999, pages 40-42 */

/*     W. Saenger, "Principles of Nucleic Acid Structure", Springer- */
/*     Verlag, 1984, page 52 */


/* Subroutine */ int nucbase_(char *resname, integer *i__, integer *c1i, 
	integer *o4i, integer *c2i, ftnlen resname_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int zatom_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *);



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




/*     ############################################################# */
/*     ##                  COPYRIGHT (C) 1999 by                  ## */
/*     ##  Marina A. Vorobieva, Nina N. Sokolova & Jay W. Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  nucleo.i  --  parameters for nucleic acid structure  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     bkbone    phosphate backbone angles for each nucleotide */
/*     glyco     glycosidic torsional angle for each nucleotide */
/*     pucker    sugar pucker, either 2=2'-endo or 3=3'-endo */
/*     dblhlx    flag to mark system as nucleic acid double helix */
/*     deoxy     flag to mark deoxyribose or ribose sugar units */
/*     hlxform   helix form (A, B or Z) of polynucleotide strands */




/*     adenine in adenosine residue  (A) */

    if (s_cmp(resname, "A  ", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__1017, &c_b284, &c_b505, &c_b506, c1i, o4i, c2i, &c__1);
	d__1 = nucleo_1.glyco[*i__ - 1] + 180.;
	i__1 = atoms_1.n - 1;
	zatom_(&c__1021, &c_b509, &c_b510, &d__1, &i__1, c1i, o4i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__1020, &c_b513, &c_b514, &c_b128, &i__1, &i__2, c1i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1019, &c_b518, &c_b519, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1025, &c_b523, &c_b524, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1027, &c_b528, &c_b529, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1024, &c_b533, &c_b534, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1023, &c_b538, &c_b539, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1022, &c_b543, &c_b544, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1018, &c_b533, &c_b549, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 7;
	zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 10;
	zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 7;
	zatom_(&c__1030, &c_b565, &c_b566, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	zatom_(&c__1028, &c_b570, &c_b251, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	zatom_(&c__1029, &c_b570, &c_b251, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1026, &c_b565, &c_b581, &c_b128, &i__1, &i__2, &i__3, &
		c__0);

/*     guanine in guanosine residue  (G) */

    } else if (s_cmp(resname, "G  ", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__1047, &c_b284, &c_b505, &c_b506, c1i, o4i, c2i, &c__1);
	d__1 = nucleo_1.glyco[*i__ - 1] + 180.;
	i__1 = atoms_1.n - 1;
	zatom_(&c__1051, &c_b591, &c_b510, &d__1, &i__1, c1i, o4i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__1050, &c_b595, &c_b596, &c_b128, &i__1, &i__2, c1i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1049, &c_b518, &c_b601, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1055, &c_b523, &c_b606, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1060, &c_b610, &c_b611, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1054, &c_b523, &c_b616, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1053, &c_b591, &c_b621, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1057, &c_b528, &c_b626, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1052, &c_b538, &c_b631, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1048, &c_b635, &c_b636, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 8;
	zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 11;
	zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 8;
	zatom_(&c__1061, &c_b565, &c_b653, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	zatom_(&c__1056, &c_b570, &c_b534, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	zatom_(&c__1058, &c_b570, &c_b251, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	zatom_(&c__1059, &c_b570, &c_b251, &c_b128, &i__1, &i__2, &i__3, &
		c__0);

/*     cytosine in cytidine residue  (C) */

    } else if (s_cmp(resname, "C ", (ftnlen)3, (ftnlen)2) == 0) {
	zatom_(&c__1078, &c_b284, &c_b505, &c_b506, c1i, o4i, c2i, &c__1);
	i__1 = atoms_1.n - 1;
	zatom_(&c__1079, &c_b509, &c_b679, &nucleo_1.glyco[*i__ - 1], &i__1, 
		c1i, o4i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__1084, &c_b682, &c_b683, &c_b108, &i__1, &i__2, c1i, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	zatom_(&c__1080, &c_b591, &c_b688, &c_b128, &i__1, &i__2, c1i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1081, &c_b528, &c_b693, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1085, &c_b543, &c_b698, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 5;
	zatom_(&c__1082, &c_b229, &c_b703, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1083, &c_b635, &c_b708, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 8;
	zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 5;
	zatom_(&c__1086, &c_b570, &c_b251, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 6;
	zatom_(&c__1087, &c_b570, &c_b251, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	zatom_(&c__1088, &c_b565, &c_b703, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 7;
	zatom_(&c__1089, &c_b565, &c_b734, &c_b128, &i__1, &i__2, &i__3, &
		c__0);

/*     uracil in uridine residue  (U) */

    } else if (s_cmp(resname, "U  ", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__1106, &c_b284, &c_b505, &c_b506, c1i, o4i, c2i, &c__1);
	i__1 = atoms_1.n - 1;
	zatom_(&c__1107, &c_b591, &c_b745, &nucleo_1.glyco[*i__ - 1], &i__1, 
		c1i, o4i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__1112, &c_b748, &c_b749, &c_b108, &i__1, &i__2, c1i, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	zatom_(&c__1108, &c_b509, &c_b754, &c_b128, &i__1, &i__2, c1i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1109, &c_b591, &c_b759, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1114, &c_b610, &c_b764, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 5;
	zatom_(&c__1110, &c_b161, &c_b769, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1111, &c_b528, &c_b774, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 8;
	zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	zatom_(&c__1113, &c_b570, &c_b785, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 6;
	zatom_(&c__1115, &c_b565, &c_b790, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 6;
	zatom_(&c__1116, &c_b565, &c_b795, &c_b128, &i__1, &i__2, &i__3, &
		c__0);

/*     adenine in deoxyadenosine residue  (DA) */

    } else if (s_cmp(resname, "DA ", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__1132, &c_b284, &c_b505, &c_b506, c1i, o4i, c2i, &c__1);
	d__1 = nucleo_1.glyco[*i__ - 1] + 180.;
	i__1 = atoms_1.n - 1;
	zatom_(&c__1136, &c_b509, &c_b510, &d__1, &i__1, c1i, o4i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__1135, &c_b513, &c_b514, &c_b128, &i__1, &i__2, c1i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1134, &c_b518, &c_b519, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1140, &c_b523, &c_b524, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1142, &c_b528, &c_b529, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1139, &c_b533, &c_b534, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1138, &c_b538, &c_b539, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1137, &c_b543, &c_b544, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1133, &c_b533, &c_b549, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 7;
	zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 10;
	zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 7;
	zatom_(&c__1145, &c_b565, &c_b566, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	zatom_(&c__1143, &c_b570, &c_b251, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	zatom_(&c__1144, &c_b570, &c_b251, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1141, &c_b565, &c_b581, &c_b128, &i__1, &i__2, &i__3, &
		c__0);

/*     guanine in deoxyguanosine residue  (DG) */

    } else if (s_cmp(resname, "DG ", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__1161, &c_b284, &c_b505, &c_b506, c1i, o4i, c2i, &c__1);
	d__1 = nucleo_1.glyco[*i__ - 1] + 180.;
	i__1 = atoms_1.n - 1;
	zatom_(&c__1165, &c_b591, &c_b510, &d__1, &i__1, c1i, o4i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__1164, &c_b595, &c_b596, &c_b128, &i__1, &i__2, c1i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1163, &c_b518, &c_b601, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1169, &c_b523, &c_b606, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 3;
	zatom_(&c__1174, &c_b610, &c_b611, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1168, &c_b523, &c_b616, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1167, &c_b591, &c_b621, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1171, &c_b528, &c_b626, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1166, &c_b538, &c_b631, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1162, &c_b635, &c_b636, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 8;
	zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 11;
	zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 8;
	zatom_(&c__1175, &c_b565, &c_b653, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	zatom_(&c__1170, &c_b570, &c_b534, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	zatom_(&c__1172, &c_b570, &c_b251, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	zatom_(&c__1173, &c_b570, &c_b251, &c_b128, &i__1, &i__2, &i__3, &
		c__0);

/*     cytosine in deoxycytidine residue  (DC) */

    } else if (s_cmp(resname, "DC ", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__1191, &c_b284, &c_b505, &c_b506, c1i, o4i, c2i, &c__1);
	i__1 = atoms_1.n - 1;
	zatom_(&c__1192, &c_b509, &c_b679, &nucleo_1.glyco[*i__ - 1], &i__1, 
		c1i, o4i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__1197, &c_b682, &c_b683, &c_b108, &i__1, &i__2, c1i, &c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	zatom_(&c__1193, &c_b591, &c_b688, &c_b128, &i__1, &i__2, c1i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1194, &c_b528, &c_b693, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1198, &c_b543, &c_b698, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 5;
	zatom_(&c__1195, &c_b229, &c_b703, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1196, &c_b635, &c_b708, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 8;
	zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 5;
	zatom_(&c__1199, &c_b570, &c_b251, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 6;
	zatom_(&c__1200, &c_b570, &c_b251, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	zatom_(&c__1201, &c_b565, &c_b703, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 7;
	zatom_(&c__1202, &c_b565, &c_b734, &c_b128, &i__1, &i__2, &i__3, &
		c__0);

/*     thymine in deoxythymidine residue  (DT) */

    } else if (s_cmp(resname, "DT ", (ftnlen)3, (ftnlen)3) == 0) {
	zatom_(&c__1218, &c_b284, &c_b505, &c_b506, c1i, o4i, c2i, &c__1);
	i__1 = atoms_1.n - 1;
	zatom_(&c__1219, &c_b509, &c_b745, &nucleo_1.glyco[*i__ - 1], &i__1, 
		c1i, o4i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	zatom_(&c__1224, &c_b748, &c_b1045, &c_b108, &i__1, &i__2, c1i, &c__0)
		;
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	zatom_(&c__1220, &c_b591, &c_b581, &c_b128, &i__1, &i__2, c1i, &c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1221, &c_b591, &c_b1055, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 2;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1226, &c_b610, &c_b1060, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 5;
	zatom_(&c__1222, &c_b161, &c_b1065, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	zatom_(&c__1227, &c_b1069, &c_b1070, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 5;
	zatom_(&c__1223, &c_b528, &c_b1075, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 9;
	zatom_(&c_n1, &c_b108, &c_b108, &c_b108, &i__1, &i__2, &c__0, &c__0);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	zatom_(&c__1225, &c_b570, &c_b1086, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 6;
	zatom_(&c__1228, &c_b237, &c_b230, &c_b108, &i__1, &i__2, &i__3, &
		c__0);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 1;
	zatom_(&c__1228, &c_b237, &c_b230, &c_b230, &i__1, &i__2, &i__3, &
		c__1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 2;
	zatom_(&c__1228, &c_b237, &c_b230, &c_b230, &i__1, &i__2, &i__3, &
		c_n1);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 9;
	zatom_(&c__1229, &c_b565, &c_b1106, &c_b128, &i__1, &i__2, &i__3, &
		c__0);
    }
    return 0;
} /* nucbase_ */



/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine watson  --  align strands of a double helix  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "watson" uses a rigid body optimization to approximately */
/*     align the paired strands of a nucleic acid double helix */


/* Subroutine */ int watson_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int rigidxyz_(void);
    static integer i__, j, ia, ib;
    static doublereal xx[1000];
    static integer kseq, nvar, list[20000]	/* was [2][10000] */;
    static doublereal dist;
    extern /* Subroutine */ int ocvm_(integer *, doublereal *, doublereal *, 
	    doublereal *, D_fp, U_fp);
    static integer root[10000], stop, nbase, nphos, iphos[10000], start, 
	    offset;
    static doublereal grdmin;
    extern /* Subroutine */ int potoff_(void), orient_(void);
    extern doublereal watson1_();
    static char resname[3];
    static doublereal minimum;
    extern /* Subroutine */ int optsave_();


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define rbc_ref(a_1,a_2) rigid_1.rbc[(a_2)*6 + a_1 - 7]
#define dfix_ref(a_1,a_2) kgeoms_1.dfix[(a_2)*3 + a_1 - 4]
#define imol_ref(a_1,a_2) molcul_1.imol[(a_2)*2 + a_1 - 3]
#define igrp_ref(a_1,a_2) group_1.igrp[(a_2)*2 + a_1 - 1]
#define tfix_ref(a_1,a_2) kgeoms_1.tfix[(a_2)*3 + a_1 - 4]
#define list_ref(a_1,a_2) list[(a_2)*2 + a_1 - 3]
#define wgrp_ref(a_1,a_2) group_1.wgrp[(a_2)*1001 + a_1 - 0]
#define idfix_ref(a_1,a_2) kgeoms_1.idfix[(a_2)*2 + a_1 - 3]
#define itfix_ref(a_1,a_2) kgeoms_1.itfix[(a_2)*4 + a_1 - 5]
#define nuclz_ref(a_0,a_1) &resdue_1.nuclz[(a_1)*3 + a_0 - 3]



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




/*     ############################################################# */
/*     ##                  COPYRIGHT (C) 1999 by                  ## */
/*     ##  Marina A. Vorobieva, Nina N. Sokolova & Jay W. Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  nucleo.i  --  parameters for nucleic acid structure  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     bkbone    phosphate backbone angles for each nucleotide */
/*     glyco     glycosidic torsional angle for each nucleotide */
/*     pucker    sugar pucker, either 2=2'-endo or 3=3'-endo */
/*     dblhlx    flag to mark system as nucleic acid double helix */
/*     deoxy     flag to mark deoxyribose or ribose sugar units */
/*     hlxform   helix form (A, B or Z) of polynucleotide strands */




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
    kgeoms_1.use_basin__ = FALSE_;
    kgeoms_1.use_wall__ = FALSE_;

/*     find root atom and hydrogen bond partners for each base */

    kseq = 0;
    nbase = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (katoms_1.atmnum[atoms_1.type__[i__ - 1] - 1] == 6 && couple_1.n12[
		i__ - 1] == 4) {
	    ia = katoms_1.atmnum[atoms_1.type__[i12_ref(1, i__) - 1] - 1];
	    ib = katoms_1.atmnum[atoms_1.type__[i12_ref(4, i__) - 1] - 1];
	    if (ia == 8 && ib == 7) {
		++nbase;
		j = i12_ref(4, i__);
		root[nbase - 1] = j;
		++kseq;
		s_copy(resname, nuclz_ref(0, sequen_1.seqtyp[kseq - 1]), (
			ftnlen)3, (ftnlen)3);
		while(s_cmp(resname, "MP ", (ftnlen)3, (ftnlen)3) == 0 || 
			s_cmp(resname, "DP ", (ftnlen)3, (ftnlen)3) == 0 || 
			s_cmp(resname, "TP ", (ftnlen)3, (ftnlen)3) == 0) {
		    ++kseq;
		    s_copy(resname, nuclz_ref(0, sequen_1.seqtyp[kseq - 1]), (
			    ftnlen)3, (ftnlen)3);
		}
		if (s_cmp(resname, "A  ", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(
			resname, "DA ", (ftnlen)3, (ftnlen)3) == 0) {
		    list_ref(1, nbase) = j + 6;
		    list_ref(2, nbase) = j + 11;
		} else if (s_cmp(resname, "G  ", (ftnlen)3, (ftnlen)3) == 0 ||
			 s_cmp(resname, "DG ", (ftnlen)3, (ftnlen)3) == 0) {
		    list_ref(1, nbase) = j + 12;
		    list_ref(2, nbase) = j + 5;
		} else if (s_cmp(resname, "C  ", (ftnlen)3, (ftnlen)3) == 0 ||
			 s_cmp(resname, "DC ", (ftnlen)3, (ftnlen)3) == 0) {
		    list_ref(1, nbase) = j + 3;
		    list_ref(2, nbase) = j + 8;
		} else if (s_cmp(resname, "U  ", (ftnlen)3, (ftnlen)3) == 0) {
		    list_ref(1, nbase) = j + 8;
		    list_ref(2, nbase) = j + 5;
		} else if (s_cmp(resname, "DT ", (ftnlen)3, (ftnlen)3) == 0) {
		    list_ref(1, nbase) = j + 9;
		    list_ref(2, nbase) = j + 5;
		}
	    }
	}
    }

/*     distance restraints for the base pair hydrogen bonds */

    i__1 = nbase / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = nbase + 1 - i__;
	++kgeoms_1.ndfix;
	idfix_ref(1, kgeoms_1.ndfix) = list_ref(1, i__);
	idfix_ref(2, kgeoms_1.ndfix) = list_ref(1, j);
	dfix_ref(1, kgeoms_1.ndfix) = 50.;
	dfix_ref(2, kgeoms_1.ndfix) = 1.85;
	dfix_ref(3, kgeoms_1.ndfix) = 1.95;
	++kgeoms_1.ndfix;
	idfix_ref(1, kgeoms_1.ndfix) = list_ref(2, i__);
	idfix_ref(2, kgeoms_1.ndfix) = list_ref(2, j);
	dfix_ref(1, kgeoms_1.ndfix) = 50.;
	dfix_ref(2, kgeoms_1.ndfix) = 1.85;
	dfix_ref(3, kgeoms_1.ndfix) = 1.95;
    }

/*     torsional restraints to enforce base pair planarity */

    i__1 = nbase / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = nbase + 1 - i__;
	++kgeoms_1.ntfix;
	itfix_ref(1, kgeoms_1.ntfix) = root[i__ - 1];
	itfix_ref(2, kgeoms_1.ntfix) = list_ref(1, i__);
	itfix_ref(3, kgeoms_1.ntfix) = list_ref(2, i__);
	itfix_ref(4, kgeoms_1.ntfix) = list_ref(1, j);
	tfix_ref(1, kgeoms_1.ntfix) = 2.5;
	tfix_ref(2, kgeoms_1.ntfix) = 180.;
	tfix_ref(3, kgeoms_1.ntfix) = 180.;
	++kgeoms_1.ntfix;
	itfix_ref(1, kgeoms_1.ntfix) = root[i__ - 1];
	itfix_ref(2, kgeoms_1.ntfix) = list_ref(2, i__);
	itfix_ref(3, kgeoms_1.ntfix) = list_ref(1, i__);
	itfix_ref(4, kgeoms_1.ntfix) = list_ref(2, j);
	tfix_ref(1, kgeoms_1.ntfix) = 2.5;
	tfix_ref(2, kgeoms_1.ntfix) = 180.;
	tfix_ref(3, kgeoms_1.ntfix) = 180.;
	++kgeoms_1.ntfix;
	itfix_ref(1, kgeoms_1.ntfix) = root[j - 1];
	itfix_ref(2, kgeoms_1.ntfix) = list_ref(1, j);
	itfix_ref(3, kgeoms_1.ntfix) = list_ref(2, j);
	itfix_ref(4, kgeoms_1.ntfix) = list_ref(1, i__);
	tfix_ref(1, kgeoms_1.ntfix) = 2.5;
	tfix_ref(2, kgeoms_1.ntfix) = 180.;
	tfix_ref(3, kgeoms_1.ntfix) = 180.;
	++kgeoms_1.ntfix;
	itfix_ref(1, kgeoms_1.ntfix) = root[j - 1];
	itfix_ref(2, kgeoms_1.ntfix) = list_ref(2, j);
	itfix_ref(3, kgeoms_1.ntfix) = list_ref(1, j);
	itfix_ref(4, kgeoms_1.ntfix) = list_ref(2, i__);
	tfix_ref(1, kgeoms_1.ntfix) = 2.5;
	tfix_ref(2, kgeoms_1.ntfix) = 180.;
	tfix_ref(3, kgeoms_1.ntfix) = 180.;
    }

/*     distance restraints between interstrand phosphates */

    nphos = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (katoms_1.atmnum[atoms_1.type__[i__ - 1] - 1] == 15) {
	    ++nphos;
	    iphos[nphos - 1] = i__;
	}
    }
    start = 1;
    stop = nphos / 2;
    s_copy(resname, nuclz_ref(0, sequen_1.seqtyp[0]), (ftnlen)3, (ftnlen)3);
    if (s_cmp(resname, "MP ", (ftnlen)3, (ftnlen)3) == 0) {
	++start;
    }
    if (s_cmp(resname, "DP ", (ftnlen)3, (ftnlen)3) == 0) {
	start += 2;
    }
    if (s_cmp(resname, "TP ", (ftnlen)3, (ftnlen)3) == 0) {
	start += 3;
    }
    s_copy(resname, nuclz_ref(0, sequen_1.seqtyp[sequen_1.nseq - 1]), (ftnlen)
	    3, (ftnlen)3);
    if (s_cmp(resname, "MP ", (ftnlen)3, (ftnlen)3) == 0) {
	--stop;
    }
    if (s_cmp(resname, "DP ", (ftnlen)3, (ftnlen)3) == 0) {
	stop += -2;
    }
    if (s_cmp(resname, "TP ", (ftnlen)3, (ftnlen)3) == 0) {
	stop += -3;
    }
    offset = stop + nphos / 2 + 1;
    if (*(unsigned char *)nucleo_1.hlxform == 'A') {
	dist = 17.78;
    }
    if (*(unsigned char *)nucleo_1.hlxform == 'B') {
	dist = 17.46;
    }
    if (*(unsigned char *)nucleo_1.hlxform == 'Z') {
	dist = 13.2;
    }
    i__1 = stop;
    for (i__ = start; i__ <= i__1; ++i__) {
	++kgeoms_1.ndfix;
	idfix_ref(1, kgeoms_1.ndfix) = iphos[i__ - 1];
	idfix_ref(2, kgeoms_1.ndfix) = iphos[offset - i__ - 1];
	dfix_ref(1, kgeoms_1.ndfix) = 100.;
	dfix_ref(2, kgeoms_1.ndfix) = dist;
	dfix_ref(3, kgeoms_1.ndfix) = dist;
    }

/*     assign each strand to a separate molecule-based group */

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
	wgrp_ref(i__, i__) = 0.;
    }

/*     get rigid body reference coordinates for each strand */

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
    ocvm_(&nvar, xx, &minimum, &grdmin, (D_fp)watson1_, (U_fp)optsave_);

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
} /* watson_ */

#undef nuclz_ref
#undef itfix_ref
#undef idfix_ref
#undef wgrp_ref
#undef list_ref
#undef tfix_ref
#undef igrp_ref
#undef imol_ref
#undef dfix_ref
#undef rbc_ref
#undef i12_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  function watson1  --  energy and gradient for watson  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "watson1" is a service routine that computes the energy */
/*     and gradient for optimally conditioned variable metric */
/*     optimization of rigid bodies */


doublereal watson1_(doublereal *xx, doublereal *g)
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
} /* watson1_ */

#undef derivs_ref
#undef rbc_ref


/* Main program alias */ int nucleic_ () { MAIN__ (); return 0; }
