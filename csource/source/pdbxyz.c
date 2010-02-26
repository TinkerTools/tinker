/* pdbxyz.f -- translated by f2c (version 20050501).
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
    doublereal mass[25000];
    integer tag[25000], class__[25000], atomic[25000], valence[25000];
    char name__[75000], story[600000];
} atmtyp_;

#define atmtyp_1 atmtyp_

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
    doublereal xpdb[25000], ypdb[25000], zpdb[25000];
    integer npdb, resnum[25000], npdb12[25000], ipdb12[200000]	/* was [8][
	    25000] */, pdblist[25000];
    char pdbtyp[150000], atmnam[100000], resnam[75000], chntyp[20], altsym[1],
	     instyp[20];
} pdb_;

#define pdb_1 pdb_

struct {
    char amino[93], nuclz[36], amino1[31], nuclz1[12];
} resdue_;

#define resdue_1 resdue_

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
    integer nseq, nchain, ichain[20000]	/* was [2][10000] */, seqtyp[10000];
    char seq[30000], chnnam[10000];
} sequen_;

#define sequen_1 sequen_

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static integer c__25000 = 25000;
static integer c__0 = 0;
static doublereal c_b63 = 1.01;
static doublereal c_b64 = 119.;
static doublereal c_b75 = 109.5;
static doublereal c_b76 = 0.;
static doublereal c_b81 = -120.;
static doublereal c_b87 = -60.;
static doublereal c_b92 = 180.;
static doublereal c_b97 = 60.;
static doublereal c_b107 = 120.9;
static doublereal c_b112 = 120.3;
static doublereal c_b129 = 1.1;
static integer c_n1 = -1;
static doublereal c_b135 = 1.12;
static doublereal c_b182 = 1.25;
static doublereal c_b183 = 117.;
static doublereal c_b184 = 126.;
static integer c__355 = 355;
static integer c__506 = 506;
static integer c__6 = 6;
static integer c__13 = 13;
static integer c__14 = 14;
static doublereal c_b222 = 110.2;
static integer c__21 = 21;
static integer c__23 = 23;
static integer c__25 = 25;
static integer c__22 = 22;
static doublereal c_b247 = 107.;
static doublereal c_b248 = 108.2;
static integer c__24 = 24;
static doublereal c_b253 = 111.6;
static integer c__26 = 26;
static integer c__33 = 33;
static integer c__35 = 35;
static integer c__37 = 37;
static integer c__39 = 39;
static integer c__34 = 34;
static doublereal c_b298 = 107.9;
static doublereal c_b299 = 110.;
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
static integer c__54 = 54;
static integer c__61 = 61;
static integer c__63 = 63;
static integer c__62 = 62;
static doublereal c_b421 = 109.2;
static integer c__64 = 64;
static doublereal c_b432 = .94;
static doublereal c_b433 = 106.9;
static integer c__71 = 71;
static integer c__73 = 73;
static integer c__75 = 75;
static integer c__72 = 72;
static integer c__74 = 74;
static integer c__76 = 76;
static integer c__83 = 83;
static integer c__85 = 85;
static integer c__84 = 84;
static doublereal c_b482 = 107.5;
static integer c__86 = 86;
static doublereal c_b492 = 1.34;
static doublereal c_b493 = 96.;
static integer c__93 = 93;
static integer c__95 = 95;
static integer c__94 = 94;
static integer c__101 = 101;
static integer c__103 = 103;
static integer c__410 = 410;
static integer c__105 = 105;
static integer c__102 = 102;
static doublereal c_b524 = 111.2;
static integer c__104 = 104;
static integer c__411 = 411;
static integer c__106 = 106;
static integer c__113 = 113;
static integer c__115 = 115;
static integer c__116 = 116;
static integer c__118 = 118;
static integer c__120 = 120;
static integer c__114 = 114;
static integer c__117 = 117;
static doublereal c_b598 = 1.09;
static doublereal c_b599 = 120.;
static integer c__119 = 119;
static integer c__121 = 121;
static integer c__128 = 128;
static integer c__130 = 130;
static integer c__131 = 131;
static integer c__133 = 133;
static integer c__135 = 135;
static integer c__136 = 136;
static integer c__129 = 129;
static integer c__132 = 132;
static integer c__134 = 134;
static integer c__137 = 137;
static doublereal c_b681 = .97;
static doublereal c_b682 = 108.;
static integer c__144 = 144;
static integer c__146 = 146;
static integer c__147 = 147;
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
static doublereal c_b727 = 126.3;
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
static doublereal c_b781 = 1.02;
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
static doublereal c_b1126 = 110.9;
static doublereal c_b1127 = 107.3;
static integer c__286 = 286;
static doublereal c_b1137 = 1.04;
static doublereal c_b1138 = 110.5;
static integer c__293 = 293;
static integer c__295 = 295;
static integer c__297 = 297;
static integer c__299 = 299;
static integer c__301 = 301;
static integer c__302 = 302;
static integer c__294 = 294;
static integer c__296 = 296;
static integer c__298 = 298;
static integer c__300 = 300;
static doublereal c_b1207 = 118.5;
static integer c__303 = 303;
static doublereal c_b1213 = 122.5;
static doublereal c_b1219 = 118.8;
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
static doublereal c_b1425 = 1.;
static integer c__1017 = 1017;
static integer c__1021 = 1021;
static integer c__1020 = 1020;
static integer c__1019 = 1019;
static integer c__1025 = 1025;
static integer c__1027 = 1027;
static integer c__1024 = 1024;
static integer c__1023 = 1023;
static integer c__1022 = 1022;
static integer c__1018 = 1018;
static integer c__1030 = 1030;
static doublereal c_b1503 = 1.08;
static doublereal c_b1504 = 123.1;
static integer c__1028 = 1028;
static integer c__1029 = 1029;
static integer c__1026 = 1026;
static doublereal c_b1522 = 115.4;
static integer c__1047 = 1047;
static integer c__1051 = 1051;
static integer c__1050 = 1050;
static integer c__1049 = 1049;
static integer c__1055 = 1055;
static integer c__1060 = 1060;
static integer c__1054 = 1054;
static integer c__1053 = 1053;
static integer c__1057 = 1057;
static integer c__1052 = 1052;
static integer c__1048 = 1048;
static integer c__1061 = 1061;
static doublereal c_b1551 = 123.;
static integer c__1056 = 1056;
static doublereal c_b1557 = 117.4;
static integer c__1058 = 1058;
static integer c__1059 = 1059;
static integer c__1078 = 1078;
static integer c__1079 = 1079;
static integer c__1084 = 1084;
static integer c__1080 = 1080;
static integer c__1081 = 1081;
static integer c__1085 = 1085;
static integer c__1082 = 1082;
static integer c__1083 = 1083;
static integer c__1086 = 1086;
static integer c__1087 = 1087;
static integer c__1088 = 1088;
static doublereal c_b1604 = 121.6;
static integer c__1089 = 1089;
static doublereal c_b1610 = 119.4;
static integer c__1106 = 1106;
static integer c__1107 = 1107;
static integer c__1112 = 1112;
static integer c__1108 = 1108;
static integer c__1109 = 1109;
static integer c__1114 = 1114;
static integer c__1110 = 1110;
static integer c__1111 = 1111;
static integer c__1113 = 1113;
static doublereal c_b1633 = 116.5;
static integer c__1115 = 1115;
static doublereal c_b1639 = 120.4;
static integer c__1116 = 1116;
static doublereal c_b1645 = 118.6;
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
static integer c__1220 = 1220;
static integer c__1221 = 1221;
static integer c__1226 = 1226;
static integer c__1222 = 1222;
static integer c__1227 = 1227;
static integer c__1223 = 1223;
static integer c__1225 = 1225;
static doublereal c_b1803 = 116.8;
static integer c__1228 = 1228;
static integer c__1229 = 1229;
static integer c__2001 = 2001;
static integer c__2002 = 2002;
static doublereal c_b1845 = .96;
static integer c__2003 = 2003;
static integer c__2004 = 2004;
static integer c__2005 = 2005;
static integer c__2006 = 2006;
static integer c__2007 = 2007;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  program pdbxyz  --  Protein Data Bank to XYZ coordinates  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "pdbxyz" takes as input a Protein Data Bank file and then */
/*     converts to and writes out a Cartesian coordinates file and, */
/*     for biopolymers, a sequence file */


/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2], i__3;
    doublereal d__1, d__2, d__3;
    olist o__1;
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_open(olist *), f_rew(alist *);
    double sqrt(doublereal);
    integer f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int ribosome_(void);
    extern integer freeunit_(void);
    static integer i__, j, it;
    static doublereal xi, yi, zi, rij;
    static integer ipdb, iseq, last;
    static doublereal rmax[10];
    static integer size[25000];
    static doublereal rcut;
    static integer next;
    extern /* Subroutine */ int sort_(integer *, integer *);
    static integer ixyz;
    extern /* Subroutine */ int field_(void), final_(void);
    static logical clash;
    extern /* Subroutine */ int delete_(integer *), ligase_(void), getpdb_(
	    void);
    static char letter[1];
    extern /* Subroutine */ int suffix_(char *, char *, ftnlen, ftnlen), 
	    upcase_(char *, ftnlen), chkxyz_(logical *), prtseq_(integer *), 
	    prtxyz_(integer *), readpdb_(integer *);
    static char pdbfile[120];
    static logical nucacid;
    static char seqfile[120];
    extern /* Subroutine */ int initial_(void);
    static logical peptide;
    static char resname[3];
    extern /* Subroutine */ int hetatom_(void), getnumb_(char *, integer *, 
	    integer *, ftnlen);
    static char reslast[3];
    extern /* Subroutine */ int version_(char *, char *, ftnlen, ftnlen);
    static char xyzfile[120];


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define amino_ref(a_0,a_1) &resdue_1.amino[(a_1)*3 + a_0 - 3]
#define nuclz_ref(a_0,a_1) &resdue_1.nuclz[(a_1)*3 + a_0 - 3]
#define atmnam_ref(a_0,a_1) &pdb_1.atmnam[(a_1)*4 + a_0 - 4]
#define resnam_ref(a_0,a_1) &pdb_1.resnam[(a_1)*3 + a_0 - 3]
#define pdbtyp_ref(a_0,a_1) &pdb_1.pdbtyp[(a_1)*6 + a_0 - 6]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  pdb.i  --  definition of a Protein Data Bank structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xpdb      x-coordinate of each atom stored in PDB format */
/*     ypdb      y-coordinate of each atom stored in PDB format */
/*     zpdb      z-coordinate of each atom stored in PDB format */
/*     npdb      number of atoms stored in Protein Data Bank format */
/*     resnum    number of the residue to which each atom belongs */
/*     npdb12    number of atoms directly bonded to each CONECT atom */
/*     ipdb12    atom numbers of atoms connected to each CONECT atom */
/*     pdblist   list of the Protein Data Bank atom number of each atom */
/*     pdbtyp    Protein Data Bank record type assigned to each atom */
/*     atmnam    Protein Data Bank atom name assigned to each atom */
/*     resnam    Protein Data Bank residue name assigned to each atom */
/*     chntyp    string with PDB chain identifiers to be included */
/*     altsym    string with PDB alternate locations to be included */
/*     instyp    string with PDB insertion records to be included */




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




/*     get the Protein Data Bank file and a parameter set */

    initial_();
    getpdb_();
    field_();

/*     decide whether the system contains only polypeptides */

    peptide = FALSE_;
    s_copy(reslast, "***", (ftnlen)3, (ftnlen)3);
    i__1 = pdb_1.npdb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (s_cmp(pdbtyp_ref(0, i__), "ATOM  ", (ftnlen)6, (ftnlen)6) == 0) {
	    s_copy(resname, resnam_ref(0, i__), (ftnlen)3, (ftnlen)3);
	    if (s_cmp(resname, reslast, (ftnlen)3, (ftnlen)3) != 0) {
		s_copy(reslast, resname, (ftnlen)3, (ftnlen)3);
		for (j = 1; j <= 31; ++j) {
		    if (s_cmp(resname, amino_ref(0, j), (ftnlen)3, (ftnlen)3) 
			    == 0) {
			peptide = TRUE_;
			goto L10;
		    }
		}
		peptide = FALSE_;
		goto L20;
L10:
		;
	    }
	}
    }
L20:

/*     decide whether the system contains only nucleic acids */

    nucacid = FALSE_;
    s_copy(reslast, "***", (ftnlen)3, (ftnlen)3);
    i__1 = pdb_1.npdb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (s_cmp(pdbtyp_ref(0, i__), "ATOM  ", (ftnlen)6, (ftnlen)6) == 0) {
	    s_copy(resname, resnam_ref(0, i__), (ftnlen)3, (ftnlen)3);
	    if (s_cmp(resname, reslast, (ftnlen)3, (ftnlen)3) != 0) {
		s_copy(reslast, resname, (ftnlen)3, (ftnlen)3);
		for (j = 1; j <= 12; ++j) {
		    if (s_cmp(resname, nuclz_ref(0, j), (ftnlen)3, (ftnlen)3) 
			    == 0) {
			nucacid = TRUE_;
			goto L30;
		    }
		}
		nucacid = FALSE_;
		goto L40;
L30:
		;
	    }
	}
    }
L40:

/*     open the TINKER coordinates file to be used for output */

    ixyz = freeunit_();
/* Writing concatenation */
    i__2[0] = files_1.leng, a__1[0] = files_1.filename;
    i__2[1] = 4, a__1[1] = ".xyz";
    s_cat(xyzfile, a__1, i__2, &c__2, (ftnlen)120);
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

/*     reopen the PDB file and read the first coordinate set */

    ipdb = freeunit_();
    s_copy(pdbfile, files_1.filename, (ftnlen)120, (ftnlen)120);
    suffix_(pdbfile, "pdb", (ftnlen)120, (ftnlen)3);
    version_(pdbfile, "old", (ftnlen)120, (ftnlen)3);
    o__1.oerr = 0;
    o__1.ounit = ipdb;
    o__1.ofnmlen = 120;
    o__1.ofnm = pdbfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    al__1.aerr = 0;
    al__1.aunit = ipdb;
    f_rew(&al__1);
    readpdb_(&ipdb);

/*     use special translation mechanisms used for biopolymers */

    while(! inform_1.abort) {
	if (peptide || nucacid) {
	    if (peptide) {
		ribosome_();
	    }
	    if (nucacid) {
		ligase_();
	    }
	    hetatom_();
	    last = atoms_1.n;
	    for (i__ = last; i__ >= 1; --i__) {
		if (atoms_1.type__[i__ - 1] == 0) {
		    delete_(&i__);
		}
	    }

/*     get general atom properties for distance-based connectivity */

	} else {
	    atoms_1.n = pdb_1.npdb;
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		atoms_1.x[i__ - 1] = pdb_1.xpdb[i__ - 1];
		atoms_1.y[i__ - 1] = pdb_1.ypdb[i__ - 1];
		atoms_1.z__[i__ - 1] = pdb_1.zpdb[i__ - 1];
		s_copy(name___ref(0, i__), atmnam_ref(1, i__), (ftnlen)3, (
			ftnlen)3);
		couple_1.n12[i__ - 1] = 0;
		next = 1;
		getnumb_(resnam_ref(0, i__), &atoms_1.type__[i__ - 1], &next, 
			(ftnlen)3);
		it = atoms_1.type__[i__ - 1];
		if (it == 0) {
		    *(unsigned char *)letter = *(unsigned char *)name___ref(0,
			     i__);
		    upcase_(letter, (ftnlen)1);
		    if (*(unsigned char *)letter == 'H') {
			size[i__ - 1] = 1;
		    } else if (*(unsigned char *)letter == 'C') {
			size[i__ - 1] = 2;
		    } else if (*(unsigned char *)letter == 'N') {
			size[i__ - 1] = 2;
		    } else if (*(unsigned char *)letter == 'O') {
			size[i__ - 1] = 2;
		    } else if (*(unsigned char *)letter == 'P') {
			size[i__ - 1] = 3;
		    } else if (*(unsigned char *)letter == 'S') {
			size[i__ - 1] = 3;
		    } else {
			size[i__ - 1] = 0;
		    }
		} else if (katoms_1.ligand[it - 1] == 0) {
		    size[i__ - 1] = 0;
		} else if (katoms_1.atmnum[it - 1] <= 2) {
		    size[i__ - 1] = 1;
		} else if (katoms_1.atmnum[it - 1] <= 10) {
		    size[i__ - 1] = 2;
		} else {
		    size[i__ - 1] = 3;
		}
	    }

/*     set the maximum bonded distance between atom type pairs */

	    rmax[0] = -1.;
	    rmax[1] = -1.;
	    rmax[2] = 1.3;
	    rmax[3] = 1.55;
	    rmax[4] = 1.75;
	    rmax[6] = 2.;
	    rmax[9] = 2.2;

/*     find and connect atom pairs within bonding distance */

	    i__1 = atoms_1.n - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		xi = atoms_1.x[i__ - 1];
		yi = atoms_1.y[i__ - 1];
		zi = atoms_1.z__[i__ - 1];
		i__3 = atoms_1.n;
		for (j = i__ + 1; j <= i__3; ++j) {
		    rcut = rmax[size[i__ - 1] * size[j - 1]];
/* Computing 2nd power */
		    d__1 = xi - atoms_1.x[j - 1];
/* Computing 2nd power */
		    d__2 = yi - atoms_1.y[j - 1];
/* Computing 2nd power */
		    d__3 = zi - atoms_1.z__[j - 1];
		    rij = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
		    if (rij <= rcut) {
			++couple_1.n12[i__ - 1];
			i12_ref(couple_1.n12[i__ - 1], i__) = j;
			++couple_1.n12[j - 1];
			i12_ref(couple_1.n12[j - 1], j) = i__;
		    }
		}
	    }
	}

/*     sort the attached atom lists into ascending order */

	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sort_(&couple_1.n12[i__ - 1], &i12_ref(1, i__));
	}

/*     check for atom pairs with identical coordinates */

	clash = FALSE_;
	chkxyz_(&clash);

/*     write the TINKER coordinates and reset the connectivities */

	prtxyz_(&ixyz);
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    couple_1.n12[i__ - 1] = 0;
	}

/*     read the next coordinate set from Protein Data Bank file */

	readpdb_(&ipdb);
    }

/*     write a sequence file for proteins and nucleic acids */

    if (peptide || nucacid) {
	iseq = freeunit_();
/* Writing concatenation */
	i__2[0] = files_1.leng, a__1[0] = files_1.filename;
	i__2[1] = 4, a__1[1] = ".seq";
	s_cat(seqfile, a__1, i__2, &c__2, (ftnlen)120);
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
    }

/*     perform any final tasks before program exit */

    cl__1.cerr = 0;
    cl__1.cunit = ixyz;
    cl__1.csta = 0;
    f_clos(&cl__1);
    final_();
    return 0;
} /* MAIN__ */

#undef pdbtyp_ref
#undef resnam_ref
#undef atmnam_ref
#undef nuclz_ref
#undef amino_ref
#undef name___ref
#undef i12_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine oldatm  --  transfer coordinates from PDB  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "oldatm" get the Cartesian coordinates for an atom from */
/*     the Protein Data Bank file, then assigns the atom type */
/*     and atomic connectivities */


/* Subroutine */ int oldatm_(integer *i__, integer *bionum, integer *i1, 
	integer *ires)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 OLDATM  --  A PDB Atom of Biotype\002,"
	    "i5,\002 is Missing in Residue\002,i5,\002-\002,a3)";

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int fatal_(void), addbond_(integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___25 = { 0, 0, 0, fmt_10, 0 };



#define seq_ref(a_0,a_1) &sequen_1.seq[(a_1)*3 + a_0 - 3]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  pdb.i  --  definition of a Protein Data Bank structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xpdb      x-coordinate of each atom stored in PDB format */
/*     ypdb      y-coordinate of each atom stored in PDB format */
/*     zpdb      z-coordinate of each atom stored in PDB format */
/*     npdb      number of atoms stored in Protein Data Bank format */
/*     resnum    number of the residue to which each atom belongs */
/*     npdb12    number of atoms directly bonded to each CONECT atom */
/*     ipdb12    atom numbers of atoms connected to each CONECT atom */
/*     pdblist   list of the Protein Data Bank atom number of each atom */
/*     pdbtyp    Protein Data Bank record type assigned to each atom */
/*     atmnam    Protein Data Bank atom name assigned to each atom */
/*     resnam    Protein Data Bank residue name assigned to each atom */
/*     chntyp    string with PDB chain identifiers to be included */
/*     altsym    string with PDB alternate locations to be included */
/*     instyp    string with PDB insertion records to be included */




/*     get coordinates, assign atom type, and update connectivities */

    if (*bionum != 0) {
	if (*i__ != 0) {
	    atoms_1.type__[atoms_1.n - 1] = fields_1.biotyp[*bionum - 1];
	    if (atoms_1.type__[atoms_1.n - 1] != 0) {
		s_copy(name___ref(0, atoms_1.n), symbol_ref(0, atoms_1.type__[
			atoms_1.n - 1]), (ftnlen)3, (ftnlen)3);
	    } else {
		s_copy(name___ref(0, atoms_1.n), "   ", (ftnlen)3, (ftnlen)3);
	    }
	    atoms_1.x[atoms_1.n - 1] = pdb_1.xpdb[*i__ - 1];
	    atoms_1.y[atoms_1.n - 1] = pdb_1.ypdb[*i__ - 1];
	    atoms_1.z__[atoms_1.n - 1] = pdb_1.zpdb[*i__ - 1];
	    if (*i1 != 0) {
		addbond_(&atoms_1.n, i1);
	    }
	    ++atoms_1.n;
	} else {
	    io___25.ciunit = iounit_1.iout;
	    s_wsfe(&io___25);
	    do_fio(&c__1, (char *)&(*bionum), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*ires), (ftnlen)sizeof(integer));
	    do_fio(&c__1, seq_ref(0, *ires), (ftnlen)3);
	    e_wsfe();
	    fatal_();
	}
    }
    return 0;
} /* oldatm_ */

#undef symbol_ref
#undef name___ref
#undef seq_ref




/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  subroutine newatm  --  create and define a new atom  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     "newatm" creates and defines an atom needed for the */
/*     Cartesian coordinates file, but which may not present */
/*     in the original Protein Data Bank file */


/* Subroutine */ int newatm_(integer *i__, integer *bionum, integer *ia, 
	doublereal *bond, integer *ib, doublereal *angle1, integer *ic, 
	doublereal *angle2, integer *chiral)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int xyzatm_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *), 
	    addbond_(integer *, integer *);


#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  pdb.i  --  definition of a Protein Data Bank structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xpdb      x-coordinate of each atom stored in PDB format */
/*     ypdb      y-coordinate of each atom stored in PDB format */
/*     zpdb      z-coordinate of each atom stored in PDB format */
/*     npdb      number of atoms stored in Protein Data Bank format */
/*     resnum    number of the residue to which each atom belongs */
/*     npdb12    number of atoms directly bonded to each CONECT atom */
/*     ipdb12    atom numbers of atoms connected to each CONECT atom */
/*     pdblist   list of the Protein Data Bank atom number of each atom */
/*     pdbtyp    Protein Data Bank record type assigned to each atom */
/*     atmnam    Protein Data Bank atom name assigned to each atom */
/*     resnam    Protein Data Bank residue name assigned to each atom */
/*     chntyp    string with PDB chain identifiers to be included */
/*     altsym    string with PDB alternate locations to be included */
/*     instyp    string with PDB insertion records to be included */




/*     set the atom type, compute coordinates, assign */
/*     connectivities and increment the atom counter */

    if (*bionum != 0) {
	atoms_1.type__[atoms_1.n - 1] = fields_1.biotyp[*bionum - 1];
	if (atoms_1.type__[atoms_1.n - 1] != 0) {
	    s_copy(name___ref(0, atoms_1.n), symbol_ref(0, atoms_1.type__[
		    atoms_1.n - 1]), (ftnlen)3, (ftnlen)3);
	} else {
	    s_copy(name___ref(0, atoms_1.n), "   ", (ftnlen)3, (ftnlen)3);
	}
	if (*i__ == 0) {
	    xyzatm_(&atoms_1.n, ia, bond, ib, angle1, ic, angle2, chiral);
	} else {
	    atoms_1.x[atoms_1.n - 1] = pdb_1.xpdb[*i__ - 1];
	    atoms_1.y[atoms_1.n - 1] = pdb_1.ypdb[*i__ - 1];
	    atoms_1.z__[atoms_1.n - 1] = pdb_1.zpdb[*i__ - 1];
	}
	addbond_(&atoms_1.n, ia);
	++atoms_1.n;
    }
    return 0;
} /* newatm_ */

#undef symbol_ref
#undef name___ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine addbond  --  add a bond between two atoms  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "addbond" adds entries to the attached atoms list in */
/*     order to generate a direct connection between two atoms */


/* Subroutine */ int addbond_(integer *i__, integer *j)
{

#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]



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




/*     add connectivity between the two atoms */

    if (*i__ != 0 && *j != 0) {
	++couple_1.n12[*i__ - 1];
	i12_ref(couple_1.n12[*i__ - 1], *i__) = *j;
	++couple_1.n12[*j - 1];
	i12_ref(couple_1.n12[*j - 1], *j) = *i__;
    }
    return 0;
} /* addbond_ */

#undef i12_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine findatm  --  locate PDB atom in a residue  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "findatm" locates a specific PDB atom name type within a */
/*     range of atoms from the PDB file, returns zero if the name */
/*     type was not found */


/* Subroutine */ int findatm_(char *name__, integer *start, integer *stop, 
	integer *ipdb, ftnlen name_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;


#define atmnam_ref(a_0,a_1) &pdb_1.atmnam[(a_1)*4 + a_0 - 4]



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
/*     ##  pdb.i  --  definition of a Protein Data Bank structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xpdb      x-coordinate of each atom stored in PDB format */
/*     ypdb      y-coordinate of each atom stored in PDB format */
/*     zpdb      z-coordinate of each atom stored in PDB format */
/*     npdb      number of atoms stored in Protein Data Bank format */
/*     resnum    number of the residue to which each atom belongs */
/*     npdb12    number of atoms directly bonded to each CONECT atom */
/*     ipdb12    atom numbers of atoms connected to each CONECT atom */
/*     pdblist   list of the Protein Data Bank atom number of each atom */
/*     pdbtyp    Protein Data Bank record type assigned to each atom */
/*     atmnam    Protein Data Bank atom name assigned to each atom */
/*     resnam    Protein Data Bank residue name assigned to each atom */
/*     chntyp    string with PDB chain identifiers to be included */
/*     altsym    string with PDB alternate locations to be included */
/*     instyp    string with PDB insertion records to be included */




/*     search for the specified atom within the residue */

    *ipdb = 0;
    i__1 = *stop;
    for (i__ = *start; i__ <= i__1; ++i__) {
	if (s_cmp(atmnam_ref(0, i__), name__, (ftnlen)4, (ftnlen)4) == 0) {
	    *ipdb = i__;
	    goto L10;
	}
    }
L10:
    return 0;
} /* findatm_ */

#undef atmnam_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine ribosome  --  coordinates from PDB polypeptide  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "ribosome" translates a polypeptide structure in Protein Data */
/*     Bank format to a Cartesian coordinate file and sequence file */


/* Subroutine */ int ribosome_(void)
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
    static char fmt_20[] = "(/,\002 RIBOSOME  --  The Maximum of\002,i8,\002"
	    " Atoms\002,\002 has been Exceeded\002)";
    static char fmt_30[] = "()";
    static char fmt_40[] = "(\002 Peptide Cyclization between Residues : "
	    " \002,2i5)";
    static char fmt_50[] = "()";
    static char fmt_60[] = "(\002 Disulfide Bond between Residues :  \002,2i"
	    "5)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static logical endchain, midchain, newchain;
    static integer i__, j, k, ci[10000], ni[10000], oi[10000], cai[10000];
    static doublereal rik;
    static integer jres, kres, nres, icys[10000], ncys, ityp, stop;
    static doublereal xcys[10000], ycys[10000], zcys[10000];
    extern /* Subroutine */ int fatal_(void);
    static integer start;
    static logical header, cyclic;
    static integer resatm[20000]	/* was [2][10000] */;
    extern /* Subroutine */ int oldatm_(integer *, integer *, integer *, 
	    integer *), newatm_(integer *, integer *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *);
    static integer cyxtyp;
    extern /* Subroutine */ int addbond_(integer *, integer *), addside_(char 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, ftnlen);
    static integer icyclic[10000];
    static char atmname[4];
    extern /* Subroutine */ int findatm_(char *, integer *, integer *, 
	    integer *, ftnlen);
    static char resname[3];
    static integer idisulf[20000]	/* was [2][10000] */, ndisulf;

    /* Fortran I/O blocks */
    static cilist io___69 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___77 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___78 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___79 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___80 = { 0, 0, 0, fmt_60, 0 };



#define seq_ref(a_0,a_1) &sequen_1.seq[(a_1)*3 + a_0 - 3]
#define amino_ref(a_0,a_1) &resdue_1.amino[(a_1)*3 + a_0 - 3]
#define ichain_ref(a_1,a_2) sequen_1.ichain[(a_2)*2 + a_1 - 3]
#define resnam_ref(a_0,a_1) &pdb_1.resnam[(a_1)*3 + a_0 - 3]
#define pdbtyp_ref(a_0,a_1) &pdb_1.pdbtyp[(a_1)*6 + a_0 - 6]
#define resatm_ref(a_1,a_2) resatm[(a_2)*2 + a_1 - 3]
#define idisulf_ref(a_1,a_2) idisulf[(a_2)*2 + a_1 - 3]



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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  pdb.i  --  definition of a Protein Data Bank structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xpdb      x-coordinate of each atom stored in PDB format */
/*     ypdb      y-coordinate of each atom stored in PDB format */
/*     zpdb      z-coordinate of each atom stored in PDB format */
/*     npdb      number of atoms stored in Protein Data Bank format */
/*     resnum    number of the residue to which each atom belongs */
/*     npdb12    number of atoms directly bonded to each CONECT atom */
/*     ipdb12    atom numbers of atoms connected to each CONECT atom */
/*     pdblist   list of the Protein Data Bank atom number of each atom */
/*     pdbtyp    Protein Data Bank record type assigned to each atom */
/*     atmnam    Protein Data Bank atom name assigned to each atom */
/*     resnam    Protein Data Bank residue name assigned to each atom */
/*     chntyp    string with PDB chain identifiers to be included */
/*     altsym    string with PDB alternate locations to be included */
/*     instyp    string with PDB insertion records to be included */




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



/*     set a pointer to the first and last atom of each residue */

    nres = 0;
    k = 0;
    i__1 = pdb_1.npdb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (s_cmp(pdbtyp_ref(0, i__), "ATOM  ", (ftnlen)6, (ftnlen)6) == 0) {
	    if (pdb_1.resnum[i__ - 1] != k) {
		k = pdb_1.resnum[i__ - 1];
		if (nres != 0) {
		    resatm_ref(2, nres) = i__ - 1;
		}
		++nres;
		resatm_ref(1, nres) = i__;
	    }
	}
    }
    if (nres != 0) {
	resatm_ref(2, nres) = pdb_1.npdb;
    }

/*     get the three-letter sequence and code for each residue */

    sequen_1.nseq = nres;
    i__1 = nres;
    for (i__ = 1; i__ <= i__1; ++i__) {
	start = resatm_ref(1, i__);
	s_copy(resname, resnam_ref(0, start), (ftnlen)3, (ftnlen)3);
	s_copy(seq_ref(0, i__), "UNK", (ftnlen)3, (ftnlen)3);
	sequen_1.seqtyp[i__ - 1] = 31;
	for (k = 1; k <= 31; ++k) {
	    if (s_cmp(resname, amino_ref(0, k), (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(seq_ref(0, i__), amino_ref(0, k), (ftnlen)3, (ftnlen)3)
			;
		sequen_1.seqtyp[i__ - 1] = k;
		goto L10;
	    }
	}
L10:
	;
    }

/*     search for the presence of cyclic polypeptide chains */

    i__1 = nres;
    for (i__ = 1; i__ <= i__1; ++i__) {
	icyclic[i__ - 1] = 0;
    }
    i__1 = sequen_1.nchain;
    for (i__ = 1; i__ <= i__1; ++i__) {
	jres = ichain_ref(1, i__);
	kres = ichain_ref(2, i__);
	start = resatm_ref(1, jres);
	stop = resatm_ref(2, jres);
	findatm_(" N  ", &start, &stop, &j, (ftnlen)4);
	start = resatm_ref(1, kres);
	stop = resatm_ref(2, kres);
	findatm_(" C  ", &start, &stop, &k, (ftnlen)4);
	if (j != 0 && k != 0) {
/* Computing 2nd power */
	    d__1 = pdb_1.xpdb[k - 1] - pdb_1.xpdb[j - 1];
/* Computing 2nd power */
	    d__2 = pdb_1.ypdb[k - 1] - pdb_1.ypdb[j - 1];
/* Computing 2nd power */
	    d__3 = pdb_1.zpdb[k - 1] - pdb_1.zpdb[j - 1];
	    rik = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	    if (rik <= 3.) {
		ni[jres - 1] = j;
		ci[kres - 1] = k;
		icyclic[jres - 1] = kres;
		icyclic[kres - 1] = jres;
	    }
	}
    }

/*     search for any potential cystine disulfide bonds */

    for (i__ = 1; i__ <= 31; ++i__) {
	if (s_cmp(amino_ref(0, i__), "CYX", (ftnlen)3, (ftnlen)3) == 0) {
	    cyxtyp = i__;
	}
    }
    ncys = 0;
    i__1 = nres;
    for (i__ = 1; i__ <= i__1; ++i__) {
	start = resatm_ref(1, i__);
	s_copy(resname, resnam_ref(0, start), (ftnlen)3, (ftnlen)3);
	if (s_cmp(resname, "CYS", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(resname,
		 "CYX", (ftnlen)3, (ftnlen)3) == 0) {
	    stop = resatm_ref(2, i__);
	    findatm_(" SG ", &start, &stop, &k, (ftnlen)4);
	    ++ncys;
	    icys[ncys - 1] = i__;
	    xcys[ncys - 1] = pdb_1.xpdb[k - 1];
	    ycys[ncys - 1] = pdb_1.ypdb[k - 1];
	    zcys[ncys - 1] = pdb_1.zpdb[k - 1];
	}
    }
    ndisulf = 0;
    i__1 = ncys - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ncys;
	for (k = i__ + 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    d__1 = xcys[i__ - 1] - xcys[k - 1];
/* Computing 2nd power */
	    d__2 = ycys[i__ - 1] - ycys[k - 1];
/* Computing 2nd power */
	    d__3 = zcys[i__ - 1] - zcys[k - 1];
	    rik = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	    if (rik <= 3.) {
		++ndisulf;
/* Computing MIN */
		i__3 = icys[i__ - 1], i__4 = icys[k - 1];
		idisulf_ref(1, ndisulf) = min(i__3,i__4);
/* Computing MAX */
		i__3 = icys[i__ - 1], i__4 = icys[k - 1];
		idisulf_ref(2, ndisulf) = max(i__3,i__4);
	    }
	}
    }
    i__1 = nres;
    for (i__ = 1; i__ <= i__1; ++i__) {
	icys[i__ - 1] = 0;
    }
    i__1 = ndisulf;
    for (i__ = 1; i__ <= i__1; ++i__) {
	icys[idisulf_ref(1, i__) - 1] = 1;
	icys[idisulf_ref(2, i__) - 1] = 1;
    }

/*     set the current atom to be the first atom */

    atoms_1.n = 1;
    newchain = TRUE_;

/*     locate and assign the atoms that make up each residue */

    i__1 = nres;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ityp = sequen_1.seqtyp[i__ - 1];
	start = resatm_ref(1, i__);
	stop = resatm_ref(2, i__);
	s_copy(resname, resnam_ref(0, start), (ftnlen)3, (ftnlen)3);

/*     check that the maximum allowed atoms is not exceeded */

	if (atoms_1.n + 25 > 25000) {
	    io___69.ciunit = iounit_1.iout;
	    s_wsfe(&io___69);
	    do_fio(&c__1, (char *)&c__25000, (ftnlen)sizeof(integer));
	    e_wsfe();
	    fatal_();
	}

/*     check for a cysteine that should be changed to cystine */

	if (icys[i__ - 1] == 1) {
	    s_copy(resname, "CYX", (ftnlen)3, (ftnlen)3);
	    sequen_1.seqtyp[i__ - 1] = cyxtyp;
	    ityp = cyxtyp;
	}

/*     test for the final residue of a peptide chain */

	endchain = FALSE_;
	i__2 = sequen_1.nchain;
	for (k = 1; k <= i__2; ++k) {
	    if (i__ == ichain_ref(2, k)) {
		endchain = TRUE_;
	    }
	}

/*     residue not at start or end of chain is in the middle */

	if (newchain || endchain) {
	    midchain = FALSE_;
	} else {
	    midchain = TRUE_;
	}

/*     check to see if the current chain is cyclic */

	cyclic = FALSE_;
	if (newchain || endchain) {
	    if (icyclic[i__ - 1] != 0) {
		cyclic = TRUE_;
	    }
	}

/*     build the amide nitrogen of the current residue */

	findatm_(" N  ", &start, &stop, &k, (ftnlen)4);
	if (k != 0) {
	    ni[i__ - 1] = atoms_1.n;
	}
	if (midchain) {
	    j = ntyp[ityp - 1];
	    oldatm_(&k, &j, &ci[i__ - 2], &i__);
	} else if (newchain) {
	    if (cyclic) {
		j = ntyp[ityp - 1];
	    } else {
		j = nntyp[ityp - 1];
	    }
	    oldatm_(&k, &j, &c__0, &i__);
	} else if (endchain) {
	    if (cyclic) {
		j = ntyp[ityp - 1];
	    } else {
		j = nctyp[ityp - 1];
	    }
	    oldatm_(&k, &j, &ci[i__ - 2], &i__);
	}

/*     build the alpha carbon of the current residue */

	s_copy(atmname, " CA ", (ftnlen)4, (ftnlen)4);
	if (s_cmp(resname, "ACE", (ftnlen)3, (ftnlen)3) == 0) {
	    s_copy(atmname, " CH3", (ftnlen)4, (ftnlen)4);
	}
	if (s_cmp(resname, "NME", (ftnlen)3, (ftnlen)3) == 0) {
	    s_copy(atmname, " CH3", (ftnlen)4, (ftnlen)4);
	}
	findatm_(atmname, &start, &stop, &k, (ftnlen)4);
	if (k != 0) {
	    cai[i__ - 1] = atoms_1.n;
	}
	if (midchain || cyclic || nres == 1) {
	    j = catyp[ityp - 1];
	    oldatm_(&k, &j, &ni[i__ - 1], &i__);
	} else if (newchain) {
	    j = cantyp[ityp - 1];
	    oldatm_(&k, &j, &ni[i__ - 1], &i__);
	} else if (endchain) {
	    j = cactyp[ityp - 1];
	    oldatm_(&k, &j, &ni[i__ - 1], &i__);
	}

/*     build the carbonyl carbon of the current residue */

	findatm_(" C  ", &start, &stop, &k, (ftnlen)4);
	if (k != 0) {
	    ci[i__ - 1] = atoms_1.n;
	}
	if (midchain || cyclic) {
	    j = ctyp[ityp - 1];
	    oldatm_(&k, &j, &cai[i__ - 1], &i__);
	} else if (newchain) {
	    j = cntyp[ityp - 1];
	    oldatm_(&k, &j, &cai[i__ - 1], &i__);
	} else if (endchain) {
	    j = cctyp[ityp - 1];
	    oldatm_(&k, &j, &cai[i__ - 1], &i__);
	}

/*     build the carbonyl oxygen of the current residue */

	findatm_(" O  ", &start, &stop, &k, (ftnlen)4);
	if (k != 0) {
	    oi[i__ - 1] = atoms_1.n;
	}
	if (midchain || cyclic) {
	    j = otyp[ityp - 1];
	    oldatm_(&k, &j, &ci[i__ - 1], &i__);
	} else if (newchain) {
	    j = ontyp[ityp - 1];
	    oldatm_(&k, &j, &ci[i__ - 1], &i__);
	} else if (endchain) {
	    j = octyp[ityp - 1];
	    oldatm_(&k, &j, &ci[i__ - 1], &i__);
	}

/*     build the amide hydrogens of the current residue */

	if (midchain || endchain && cyclic) {
	    j = hntyp[ityp - 1];
	    findatm_(" H  ", &start, &stop, &k, (ftnlen)4);
	    newatm_(&k, &j, &ni[i__ - 1], &c_b63, &ci[i__ - 2], &c_b64, &cai[
		    i__ - 1], &c_b64, &c__1);
	} else if (newchain && cyclic) {
	    j = hntyp[ityp - 1];
	    findatm_(" H  ", &start, &stop, &k, (ftnlen)4);
	    newatm_(&k, &j, &ni[i__ - 1], &c_b63, &ci[icyclic[i__ - 1] - 1], &
		    c_b64, &cai[i__ - 1], &c_b64, &c__1);
	} else if (newchain) {
	    j = hnntyp[ityp - 1];
	    if (s_cmp(resname, "PRO", (ftnlen)3, (ftnlen)3) == 0) {
		findatm_(" H2 ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &ni[i__ - 1], &c_b63, &cai[i__ - 1], &c_b75, &
			ci[i__ - 1], &c_b76, &c__0);
		findatm_(" H3 ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &ni[i__ - 1], &c_b63, &cai[i__ - 1], &c_b75, &
			ci[i__ - 1], &c_b81, &c__0);
	    } else if (s_cmp(resname, "PCA", (ftnlen)3, (ftnlen)3) == 0) {
		findatm_(" H  ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &ni[i__ - 1], &c_b63, &cai[i__ - 1], &c_b75, &
			ci[i__ - 1], &c_b87, &c__0);
	    } else {
		findatm_(" H1 ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &ni[i__ - 1], &c_b63, &cai[i__ - 1], &c_b75, &
			ci[i__ - 1], &c_b92, &c__0);
		findatm_(" H2 ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &ni[i__ - 1], &c_b63, &cai[i__ - 1], &c_b75, &
			ci[i__ - 1], &c_b97, &c__0);
		findatm_(" H3 ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &ni[i__ - 1], &c_b63, &cai[i__ - 1], &c_b75, &
			ci[i__ - 1], &c_b87, &c__0);
	    }
	} else if (endchain) {
	    j = hnctyp[ityp - 1];
	    if (s_cmp(resname, "NH2", (ftnlen)3, (ftnlen)3) == 0) {
		findatm_(" H1 ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &ni[i__ - 1], &c_b63, &ci[i__ - 2], &c_b107, &
			cai[i__ - 2], &c_b76, &c__0);
		findatm_(" H2 ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &ni[i__ - 1], &c_b63, &ci[i__ - 2], &c_b112, &
			cai[i__ - 2], &c_b92, &c__0);
	    } else if (s_cmp(resname, "NME", (ftnlen)3, (ftnlen)3) == 0) {
		findatm_(" H  ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &ni[i__ - 1], &c_b63, &ci[i__ - 2], &c_b64, &
			cai[i__ - 1], &c_b64, &c__1);
	    } else {
		findatm_(" H  ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &ni[i__ - 1], &c_b63, &ci[i__ - 2], &c_b64, &
			cai[i__ - 1], &c_b64, &c__1);
	    }
	}

/*     build the alpha hydrogen of the current residue */

	if (s_cmp(resname, "GLY", (ftnlen)3, (ftnlen)3) == 0) {
	    findatm_(" HA2", &start, &stop, &k, (ftnlen)4);
	} else {
	    findatm_(" HA ", &start, &stop, &k, (ftnlen)4);
	}
	if (midchain || cyclic) {
	    j = hatyp[ityp - 1];
	    newatm_(&k, &j, &cai[i__ - 1], &c_b129, &ni[i__ - 1], &c_b75, &ci[
		    i__ - 1], &c_b75, &c_n1);
	} else if (newchain) {
	    j = hantyp[ityp - 1];
	    if (s_cmp(resname, "FOR", (ftnlen)3, (ftnlen)3) == 0) {
		findatm_(" H  ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &ci[i__ - 1], &c_b135, &oi[i__ - 1], &c_b76, &
			c__0, &c_b76, &c__0);
	    } else if (s_cmp(resname, "ACE", (ftnlen)3, (ftnlen)3) == 0) {
		findatm_(" H1 ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &cai[i__ - 1], &c_b129, &ci[i__ - 1], &c_b75, 
			&oi[i__ - 1], &c_b92, &c__0);
		findatm_(" H2 ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &cai[i__ - 1], &c_b129, &ci[i__ - 1], &c_b75, 
			&oi[i__ - 1], &c_b97, &c__0);
		findatm_(" H3 ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &cai[i__ - 1], &c_b129, &ci[i__ - 1], &c_b75, 
			&oi[i__ - 1], &c_b87, &c__0);
	    } else {
		newatm_(&k, &j, &cai[i__ - 1], &c_b129, &ni[i__ - 1], &c_b75, 
			&ci[i__ - 1], &c_b75, &c_n1);
	    }
	} else if (endchain) {
	    j = hactyp[ityp - 1];
	    if (s_cmp(resname, "NME", (ftnlen)3, (ftnlen)3) == 0) {
		findatm_(" H1 ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &cai[i__ - 1], &c_b129, &ni[i__ - 1], &c_b75, 
			&ci[i__ - 2], &c_b92, &c__0);
		findatm_(" H2 ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &cai[i__ - 1], &c_b129, &ni[i__ - 1], &c_b75, 
			&ci[i__ - 2], &c_b97, &c__0);
		findatm_(" H3 ", &start, &stop, &k, (ftnlen)4);
		newatm_(&k, &j, &cai[i__ - 1], &c_b129, &ni[i__ - 1], &c_b75, 
			&ci[i__ - 2], &c_b87, &c__0);
	    } else {
		newatm_(&k, &j, &cai[i__ - 1], &c_b129, &ni[i__ - 1], &c_b75, 
			&ci[i__ - 1], &c_b75, &c_n1);
	    }
	}

/*     build the side chain atoms of the current residue */

	addside_(resname, &i__, &start, &stop, &cai[i__ - 1], &ni[i__ - 1], &
		ci[i__ - 1], &icys[i__ - 1], (ftnlen)3);

/*     build the terminal oxygen at the end of a peptide chain */

	if (endchain && ! cyclic) {
	    findatm_(" OXT", &start, &stop, &k, (ftnlen)4);
	    if (k == 0) {
		findatm_(" OT2", &start, &stop, &k, (ftnlen)4);
	    }
	    j = octyp[ityp - 1];
	    newatm_(&k, &j, &ci[i__ - 1], &c_b182, &cai[i__ - 1], &c_b183, &
		    oi[i__ - 1], &c_b184, &c__1);
	}

/*     next residue starts a new chain if current one just ended */

	if (endchain) {
	    newchain = TRUE_;
	} else {
	    newchain = FALSE_;
	}
    }

/*     connect the terminal residues of a cyclic polypeptide */

    header = TRUE_;
    i__1 = sequen_1.nchain;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = ichain_ref(1, i__);
	k = icyclic[j - 1];
	if (k != 0) {
	    addbond_(&ni[j - 1], &ci[k - 1]);
	    if (inform_1.verbose) {
		if (header) {
		    header = FALSE_;
		    io___77.ciunit = iounit_1.iout;
		    s_wsfe(&io___77);
		    e_wsfe();
		}
		io___78.ciunit = iounit_1.iout;
		s_wsfe(&io___78);
		do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	}
    }

/*     connect the sulfur atoms involved in disulfide bonds */

    header = TRUE_;
    i__1 = ndisulf;
    for (i__ = 1; i__ <= i__1; ++i__) {
	addbond_(&icys[idisulf_ref(1, i__) - 1], &icys[idisulf_ref(2, i__) - 
		1]);
	if (inform_1.verbose) {
	    if (header) {
		header = FALSE_;
		io___79.ciunit = iounit_1.iout;
		s_wsfe(&io___79);
		e_wsfe();
	    }
	    io___80.ciunit = iounit_1.iout;
	    s_wsfe(&io___80);
	    do_fio(&c__1, (char *)&idisulf_ref(1, i__), (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&idisulf_ref(2, i__), (ftnlen)sizeof(
		    integer));
	    e_wsfe();
	}
    }

/*     total number of atoms is one less than the current atom */

    --atoms_1.n;
    return 0;
} /* ribosome_ */

#undef idisulf_ref
#undef resatm_ref
#undef pdbtyp_ref
#undef resnam_ref
#undef ichain_ref
#undef amino_ref
#undef seq_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine addside  --  build the amino acid side chains  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "addside" builds the Cartesian coordinates for a single amino */
/*     acid side chain; coordinates are read from the Protein Data */
/*     Bank file or found from internal coordinates, then atom types */
/*     are assigned and connectivity data generated */


/* Subroutine */ int addside_(char *resname, integer *ires, integer *start, 
	integer *stop, integer *cai, integer *ni, integer *ci, integer *
	cystype, ftnlen resname_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int oldatm_(integer *, integer *, integer *, 
	    integer *), newatm_(integer *, integer *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *), 
	    addbond_(integer *, integer *), findatm_(char *, integer *, 
	    integer *, integer *, ftnlen);



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
	findatm_(" HA3", start, stop, &i__, (ftnlen)4);
	if (*ires == 1) {
	    newatm_(&i__, &c__355, cai, &c_b129, ni, &c_b75, ci, &c_b75, &
		    c__1);
	} else if (*ires == sequen_1.nseq) {
	    newatm_(&i__, &c__506, cai, &c_b129, ni, &c_b75, ci, &c_b75, &
		    c__1);
	} else {
	    newatm_(&i__, &c__6, cai, &c_b129, ni, &c_b75, ci, &c_b75, &c__1);
	}

/*     alanine residue  (ALA) */

    } else if (s_cmp(resname, "ALA", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__13, cai, ires);
	findatm_(" HB1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	newatm_(&i__, &c__14, &i__1, &c_b129, cai, &c_b222, ni, &c_b92, &c__0)
		;
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	newatm_(&i__, &c__14, &i__1, &c_b129, cai, &c_b222, ni, &c_b97, &c__0)
		;
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	newatm_(&i__, &c__14, &i__1, &c_b129, cai, &c_b222, ni, &c_b87, &c__0)
		;

/*     valine residue  (VAL) */

    } else if (s_cmp(resname, "VAL", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__21, cai, ires);
	findatm_(" CG1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__23, &i__1, ires);
	findatm_(" CG2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__25, &i__1, ires);
	findatm_(" HB ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 2;
	newatm_(&i__, &c__22, &i__1, &c_b129, cai, &c_b247, &i__2, &c_b248, &
		c__1);
	findatm_("HG11", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	newatm_(&i__, &c__24, &i__1, &c_b129, &i__2, &c_b253, cai, &c_b92, &
		c__0);
	findatm_("HG12", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	newatm_(&i__, &c__24, &i__1, &c_b129, &i__2, &c_b253, cai, &c_b97, &
		c__0);
	findatm_("HG13", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	newatm_(&i__, &c__24, &i__1, &c_b129, &i__2, &c_b253, cai, &c_b87, &
		c__0);
	findatm_("HG21", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	newatm_(&i__, &c__26, &i__1, &c_b129, &i__2, &c_b253, cai, &c_b92, &
		c__0);
	findatm_("HG22", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	newatm_(&i__, &c__26, &i__1, &c_b129, &i__2, &c_b253, cai, &c_b97, &
		c__0);
	findatm_("HG23", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	newatm_(&i__, &c__26, &i__1, &c_b129, &i__2, &c_b253, cai, &c_b87, &
		c__0);

/*     leucine residue  (LEU) */

    } else if (s_cmp(resname, "LEU", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__33, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__35, &i__1, ires);
	findatm_(" CD1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__37, &i__1, ires);
	findatm_(" CD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__39, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	newatm_(&i__, &c__34, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	newatm_(&i__, &c__34, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_(" HG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 4;
	newatm_(&i__, &c__36, &i__1, &c_b129, &i__2, &c_b247, &i__3, &c_b248, 
		&c__1);
	findatm_("HD11", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__38, &i__1, &c_b129, &i__2, &c_b253, &i__3, &c_b92, &
		c__0);
	findatm_("HD12", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__38, &i__1, &c_b129, &i__2, &c_b253, &i__3, &c_b97, &
		c__0);
	findatm_("HD13", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__38, &i__1, &c_b129, &i__2, &c_b253, &i__3, &c_b87, &
		c__0);
	findatm_("HD21", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 10;
	newatm_(&i__, &c__40, &i__1, &c_b129, &i__2, &c_b253, &i__3, &c_b92, &
		c__0);
	findatm_("HD22", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 11;
	newatm_(&i__, &c__40, &i__1, &c_b129, &i__2, &c_b253, &i__3, &c_b97, &
		c__0);
	findatm_("HD23", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 12;
	newatm_(&i__, &c__40, &i__1, &c_b129, &i__2, &c_b253, &i__3, &c_b87, &
		c__0);

/*     isoleucine residue  (ILE) */

    } else if (s_cmp(resname, "ILE", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__47, cai, ires);
	findatm_(" CG1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__49, &i__1, ires);
	findatm_(" CG2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__51, &i__1, ires);
	findatm_(" CD1", start, stop, &i__, (ftnlen)4);
	if (i__ == 0) {
	    findatm_(" CD ", start, stop, &i__, (ftnlen)4);
	}
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__53, &i__1, ires);
	findatm_(" HB ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	newatm_(&i__, &c__48, &i__1, &c_b129, cai, &c_b247, &i__2, &c_b248, &
		c_n1);
	findatm_("HG12", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 2;
	newatm_(&i__, &c__50, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c__1);
	findatm_("HG13", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 3;
	newatm_(&i__, &c__50, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c_n1);
	findatm_("HG21", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	newatm_(&i__, &c__52, &i__1, &c_b129, &i__2, &c_b253, cai, &c_b92, &
		c__0);
	findatm_("HG22", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	newatm_(&i__, &c__52, &i__1, &c_b129, &i__2, &c_b253, cai, &c_b97, &
		c__0);
	findatm_("HG23", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	newatm_(&i__, &c__52, &i__1, &c_b129, &i__2, &c_b253, cai, &c_b87, &
		c__0);
	findatm_("HD11", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 10;
	newatm_(&i__, &c__54, &i__1, &c_b129, &i__2, &c_b253, &i__3, &c_b92, &
		c__0);
	findatm_("HD12", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 11;
	newatm_(&i__, &c__54, &i__1, &c_b129, &i__2, &c_b253, &i__3, &c_b97, &
		c__0);
	findatm_("HD13", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 12;
	newatm_(&i__, &c__54, &i__1, &c_b129, &i__2, &c_b253, &i__3, &c_b87, &
		c__0);

/*     serine residue  (SER) */

    } else if (s_cmp(resname, "SER", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__61, cai, ires);
	findatm_(" OG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__63, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	newatm_(&i__, &c__62, &i__1, &c_b129, cai, &c_b421, &i__2, &c_b75, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 2;
	newatm_(&i__, &c__62, &i__1, &c_b129, cai, &c_b421, &i__2, &c_b75, &
		c_n1);
	findatm_(" HG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	newatm_(&i__, &c__64, &i__1, &c_b432, &i__2, &c_b433, cai, &c_b92, &
		c__0);

/*     threonine residue  (THR) */

    } else if (s_cmp(resname, "THR", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__71, cai, ires);
	findatm_(" OG1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__73, &i__1, ires);
	findatm_(" CG2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__75, &i__1, ires);
	findatm_(" HB ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 2;
	newatm_(&i__, &c__72, &i__1, &c_b129, cai, &c_b247, &i__2, &c_b248, &
		c_n1);
	findatm_(" HG1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	newatm_(&i__, &c__74, &i__1, &c_b432, &i__2, &c_b433, cai, &c_b92, &
		c__0);
	findatm_("HG21", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 5;
	newatm_(&i__, &c__76, &i__1, &c_b129, &i__2, &c_b253, cai, &c_b92, &
		c__0);
	findatm_("HG22", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 6;
	newatm_(&i__, &c__76, &i__1, &c_b129, &i__2, &c_b253, cai, &c_b97, &
		c__0);
	findatm_("HG23", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	newatm_(&i__, &c__76, &i__1, &c_b129, &i__2, &c_b253, cai, &c_b87, &
		c__0);

/*     cysteine residue  (CYS) */

    } else if (s_cmp(resname, "CYS", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__83, cai, ires);
	findatm_(" SG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__85, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	newatm_(&i__, &c__84, &i__1, &c_b129, cai, &c_b75, &i__2, &c_b482, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 2;
	newatm_(&i__, &c__84, &i__1, &c_b129, cai, &c_b75, &i__2, &c_b482, &
		c_n1);
	findatm_(" HG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	newatm_(&i__, &c__86, &i__1, &c_b492, &i__2, &c_b493, cai, &c_b92, &
		c__0);

/*     cystine residue  (CYX) */

    } else if (s_cmp(resname, "CYX", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__93, cai, ires);
	findatm_(" SG ", start, stop, &i__, (ftnlen)4);
	*cystype = atoms_1.n;
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__95, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	i__2 = atoms_1.n - 1;
	newatm_(&i__, &c__94, &i__1, &c_b129, cai, &c_b75, &i__2, &c_b482, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 2;
	newatm_(&i__, &c__94, &i__1, &c_b129, cai, &c_b75, &i__2, &c_b482, &
		c_n1);

/*     proline residue  (PRO) */

    } else if (s_cmp(resname, "PRO", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__101, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__103, &i__1, ires);
	findatm_(" CD ", start, stop, &i__, (ftnlen)4);
	if (*ires == 1) {
	    i__1 = atoms_1.n - 1;
	    oldatm_(&i__, &c__410, &i__1, ires);
	} else {
	    i__1 = atoms_1.n - 1;
	    oldatm_(&i__, &c__105, &i__1, ires);
	}
	i__1 = atoms_1.n - 1;
	addbond_(&i__1, ni);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 2;
	newatm_(&i__, &c__102, &i__1, &c_b129, cai, &c_b524, &i__2, &c_b524, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	newatm_(&i__, &c__102, &i__1, &c_b129, cai, &c_b524, &i__2, &c_b524, &
		c_n1);
	findatm_(" HG2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 3;
	newatm_(&i__, &c__104, &i__1, &c_b129, &i__2, &c_b524, &i__3, &c_b524,
		 &c__1);
	findatm_(" HG3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 4;
	newatm_(&i__, &c__104, &i__1, &c_b129, &i__2, &c_b524, &i__3, &c_b524,
		 &c_n1);
	if (*ires == 1) {
	    findatm_(" HD2", start, stop, &i__, (ftnlen)4);
	    i__1 = atoms_1.n - 5;
	    i__2 = atoms_1.n - 6;
	    newatm_(&i__, &c__411, &i__1, &c_b129, &i__2, &c_b524, ni, &
		    c_b524, &c__1);
	    findatm_(" HD3", start, stop, &i__, (ftnlen)4);
	    i__1 = atoms_1.n - 6;
	    i__2 = atoms_1.n - 7;
	    newatm_(&i__, &c__411, &i__1, &c_b129, &i__2, &c_b524, ni, &
		    c_b524, &c_n1);
	} else {
	    findatm_(" HD2", start, stop, &i__, (ftnlen)4);
	    i__1 = atoms_1.n - 5;
	    i__2 = atoms_1.n - 6;
	    newatm_(&i__, &c__106, &i__1, &c_b129, &i__2, &c_b524, ni, &
		    c_b524, &c__1);
	    findatm_(" HD3", start, stop, &i__, (ftnlen)4);
	    i__1 = atoms_1.n - 6;
	    i__2 = atoms_1.n - 7;
	    newatm_(&i__, &c__106, &i__1, &c_b129, &i__2, &c_b524, ni, &
		    c_b524, &c_n1);
	}

/*     phenylalanine residue  (PHE) */

    } else if (s_cmp(resname, "PHE", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__113, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__115, &i__1, ires);
	findatm_(" CD1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__116, &i__1, ires);
	findatm_(" CD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__116, &i__1, ires);
	findatm_(" CE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__118, &i__1, ires);
	findatm_(" CE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__118, &i__1, ires);
	findatm_(" CZ ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__120, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	addbond_(&i__1, &i__2);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 6;
	newatm_(&i__, &c__114, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 7;
	newatm_(&i__, &c__114, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_(" HD1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__117, &i__1, &c_b598, &i__2, &c_b599, &i__3, &c_b76, 
		&c__0);
	findatm_(" HD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 10;
	newatm_(&i__, &c__117, &i__1, &c_b598, &i__2, &c_b599, &i__3, &c_b76, 
		&c__0);
	findatm_(" HE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 10;
	newatm_(&i__, &c__119, &i__1, &c_b598, &i__2, &c_b599, &i__3, &c_b92, 
		&c__0);
	findatm_(" HE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 11;
	newatm_(&i__, &c__119, &i__1, &c_b598, &i__2, &c_b599, &i__3, &c_b92, 
		&c__0);
	findatm_(" HZ ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 10;
	newatm_(&i__, &c__121, &i__1, &c_b598, &i__2, &c_b599, &i__3, &c_b92, 
		&c__0);

/*     tyrosine residue  (TYR) */

    } else if (s_cmp(resname, "TYR", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__128, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__130, &i__1, ires);
	findatm_(" CD1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__131, &i__1, ires);
	findatm_(" CD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__131, &i__1, ires);
	findatm_(" CE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__133, &i__1, ires);
	findatm_(" CE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__133, &i__1, ires);
	findatm_(" CZ ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__135, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	addbond_(&i__1, &i__2);
	findatm_(" OH ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__136, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 7;
	newatm_(&i__, &c__129, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 8;
	newatm_(&i__, &c__129, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_(" HD1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 10;
	newatm_(&i__, &c__132, &i__1, &c_b598, &i__2, &c_b599, &i__3, &c_b76, 
		&c__0);
	findatm_(" HD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 11;
	newatm_(&i__, &c__132, &i__1, &c_b598, &i__2, &c_b599, &i__3, &c_b76, 
		&c__0);
	findatm_(" HE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 11;
	newatm_(&i__, &c__134, &i__1, &c_b598, &i__2, &c_b599, &i__3, &c_b92, 
		&c__0);
	findatm_(" HE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 12;
	newatm_(&i__, &c__134, &i__1, &c_b598, &i__2, &c_b599, &i__3, &c_b92, 
		&c__0);
	findatm_(" HH ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__137, &i__1, &c_b681, &i__2, &c_b682, &i__3, &c_b76, 
		&c__0);

/*     tryptophan residue  (TRP) */

    } else if (s_cmp(resname, "TRP", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__144, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__146, &i__1, ires);
	findatm_(" CD1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__147, &i__1, ires);
	findatm_(" CD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__149, &i__1, ires);
	findatm_(" NE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__150, &i__1, ires);
	findatm_(" CE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__152, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	addbond_(&i__1, &i__2);
	findatm_(" CE3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	oldatm_(&i__, &c__153, &i__1, ires);
	findatm_(" CZ2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__155, &i__1, ires);
	findatm_(" CZ3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__157, &i__1, ires);
	findatm_(" CH2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__159, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	addbond_(&i__1, &i__2);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 9;
	newatm_(&i__, &c__145, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 11;
	i__2 = atoms_1.n - 10;
	newatm_(&i__, &c__145, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_(" HD1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 12;
	newatm_(&i__, &c__148, &i__1, &c_b598, &i__2, &c_b184, &i__3, &c_b76, 
		&c__0);
	findatm_(" HE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 12;
	newatm_(&i__, &c__151, &i__1, &c_b63, &i__2, &c_b727, &i__3, &c_b92, &
		c__0);
	findatm_(" HE3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 5;
	newatm_(&i__, &c__154, &i__1, &c_b598, &i__2, &c_b599, &i__3, &c_b92, 
		&c__0);
	findatm_(" HZ2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__156, &i__1, &c_b598, &i__2, &c_b599, &i__3, &c_b92, 
		&c__0);
	findatm_(" HZ3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__158, &i__1, &c_b598, &i__2, &c_b599, &i__3, &c_b92, 
		&c__0);
	findatm_(" HH2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 11;
	newatm_(&i__, &c__160, &i__1, &c_b598, &i__2, &c_b599, &i__3, &c_b92, 
		&c__0);

/*     histidine (HD and HE) residue  (HIS) */

    } else if (s_cmp(resname, "HIS", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__167, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__169, &i__1, ires);
	findatm_(" ND1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__170, &i__1, ires);
	findatm_(" CD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__172, &i__1, ires);
	findatm_(" CE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__174, &i__1, ires);
	findatm_(" NE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__176, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	addbond_(&i__1, &i__2);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	newatm_(&i__, &c__168, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 6;
	newatm_(&i__, &c__168, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_(" HD1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 3;
	newatm_(&i__, &c__171, &i__1, &c_b781, &i__2, &c_b184, &i__3, &c_b92, 
		&c__0);
	findatm_(" HD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 5;
	newatm_(&i__, &c__173, &i__1, &c_b598, &i__2, &c_b184, &i__3, &c_b92, 
		&c__0);
	findatm_(" HE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__175, &i__1, &c_b598, &i__2, &c_b184, &i__3, &c_b92, 
		&c__0);
	findatm_(" HE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__177, &i__1, &c_b781, &i__2, &c_b184, &i__3, &c_b92, 
		&c__0);

/*     histidine (HD only) residue  (HID) */

    } else if (s_cmp(resname, "HID", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__184, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__186, &i__1, ires);
	findatm_(" ND1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__187, &i__1, ires);
	findatm_(" CD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__189, &i__1, ires);
	findatm_(" CE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__191, &i__1, ires);
	findatm_(" NE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__193, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	addbond_(&i__1, &i__2);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	newatm_(&i__, &c__185, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 6;
	newatm_(&i__, &c__185, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_(" HD1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 3;
	newatm_(&i__, &c__188, &i__1, &c_b781, &i__2, &c_b184, &i__3, &c_b92, 
		&c__0);
	findatm_(" HD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 5;
	newatm_(&i__, &c__190, &i__1, &c_b598, &i__2, &c_b184, &i__3, &c_b92, 
		&c__0);
	findatm_(" HE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__192, &i__1, &c_b598, &i__2, &c_b184, &i__3, &c_b92, 
		&c__0);

/*     histidine (HE only) residue  (HIE) */

    } else if (s_cmp(resname, "HIE", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__200, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__202, &i__1, ires);
	findatm_(" ND1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__203, &i__1, ires);
	findatm_(" CD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__204, &i__1, ires);
	findatm_(" CE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__206, &i__1, ires);
	findatm_(" NE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__208, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 3;
	addbond_(&i__1, &i__2);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	newatm_(&i__, &c__201, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 6;
	newatm_(&i__, &c__201, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_(" HD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 3;
	i__3 = atoms_1.n - 4;
	newatm_(&i__, &c__205, &i__1, &c_b598, &i__2, &c_b184, &i__3, &c_b92, 
		&c__0);
	findatm_(" HE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 6;
	newatm_(&i__, &c__207, &i__1, &c_b598, &i__2, &c_b184, &i__3, &c_b92, 
		&c__0);
	findatm_(" HE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__209, &i__1, &c_b781, &i__2, &c_b184, &i__3, &c_b92, 
		&c__0);

/*     aspartic acid residue  (ASP) */

    } else if (s_cmp(resname, "ASP", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__216, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__218, &i__1, ires);
	findatm_(" OD1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__219, &i__1, ires);
	findatm_(" OD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__219, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	newatm_(&i__, &c__217, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	newatm_(&i__, &c__217, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);

/*     asparagine residue  (ASN) */

    } else if (s_cmp(resname, "ASN", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__226, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__228, &i__1, ires);
	findatm_(" OD1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__229, &i__1, ires);
	findatm_(" ND2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__230, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	newatm_(&i__, &c__227, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	newatm_(&i__, &c__227, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_("HD21", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 6;
	newatm_(&i__, &c__231, &i__1, &c_b63, &i__2, &c_b107, &i__3, &c_b76, &
		c__0);
	findatm_("HD22", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__231, &i__1, &c_b63, &i__2, &c_b112, &i__3, &c_b92, &
		c__0);

/*     glutamic acid residue  (GLU) */

    } else if (s_cmp(resname, "GLU", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__238, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__240, &i__1, ires);
	findatm_(" CD ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__242, &i__1, ires);
	findatm_(" OE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__243, &i__1, ires);
	findatm_(" OE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__243, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	newatm_(&i__, &c__239, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	newatm_(&i__, &c__239, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_(" HG2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 5;
	newatm_(&i__, &c__241, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c__1);
	findatm_(" HG3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 6;
	newatm_(&i__, &c__241, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c_n1);

/*     glutamine residue  (GLN) */

    } else if (s_cmp(resname, "GLN", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__250, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__252, &i__1, ires);
	findatm_(" CD ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__254, &i__1, ires);
	findatm_(" OE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__255, &i__1, ires);
	findatm_(" NE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__256, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	newatm_(&i__, &c__251, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	newatm_(&i__, &c__251, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_(" HG2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 5;
	newatm_(&i__, &c__253, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c__1);
	findatm_(" HG3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 6;
	newatm_(&i__, &c__253, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c_n1);
	findatm_("HE21", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__257, &i__1, &c_b63, &i__2, &c_b107, &i__3, &c_b76, &
		c__0);
	findatm_("HE22", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__257, &i__1, &c_b63, &i__2, &c_b112, &i__3, &c_b92, &
		c__0);

/*     methionine residue  (MET) */

    } else if (s_cmp(resname, "MET", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__264, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__266, &i__1, ires);
	findatm_(" SD ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__268, &i__1, ires);
	findatm_(" CE ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__269, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	newatm_(&i__, &c__265, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	newatm_(&i__, &c__265, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_(" HG2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 4;
	newatm_(&i__, &c__267, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c__1);
	findatm_(" HG3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 5;
	newatm_(&i__, &c__267, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c_n1);
	findatm_(" HE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__270, &i__1, &c_b129, &i__2, &c_b222, &i__3, &c_b92, 
		&c__0);
	findatm_(" HE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__270, &i__1, &c_b129, &i__2, &c_b222, &i__3, &c_b97, 
		&c__0);
	findatm_(" HE3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__270, &i__1, &c_b129, &i__2, &c_b222, &i__3, &c_b87, 
		&c__0);

/*     lysine residue  (LYS) */

    } else if (s_cmp(resname, "LYS", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__277, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__279, &i__1, ires);
	findatm_(" CD ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__281, &i__1, ires);
	findatm_(" CE ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__283, &i__1, ires);
	findatm_(" NZ ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__285, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	newatm_(&i__, &c__278, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	newatm_(&i__, &c__278, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_(" HG2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 5;
	newatm_(&i__, &c__280, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c__1);
	findatm_(" HG3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 6;
	newatm_(&i__, &c__280, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c_n1);
	findatm_(" HD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 6;
	newatm_(&i__, &c__282, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c__1);
	findatm_(" HD3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__282, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c_n1);
	findatm_(" HE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__284, &i__1, &c_b129, &i__2, &c_b1126, &i__3, &
		c_b1127, &c__1);
	findatm_(" HE3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__284, &i__1, &c_b129, &i__2, &c_b1126, &i__3, &
		c_b1127, &c_n1);
	findatm_(" HZ1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 11;
	newatm_(&i__, &c__286, &i__1, &c_b1137, &i__2, &c_b1138, &i__3, &
		c_b92, &c__0);
	findatm_(" HZ2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 12;
	newatm_(&i__, &c__286, &i__1, &c_b1137, &i__2, &c_b1138, &i__3, &
		c_b97, &c__0);
	findatm_(" HZ3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 11;
	i__2 = atoms_1.n - 12;
	i__3 = atoms_1.n - 13;
	newatm_(&i__, &c__286, &i__1, &c_b1137, &i__2, &c_b1138, &i__3, &
		c_b87, &c__0);

/*     arginine residue  (ARG) */

    } else if (s_cmp(resname, "ARG", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__293, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__295, &i__1, ires);
	findatm_(" CD ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__297, &i__1, ires);
	findatm_(" NE ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__299, &i__1, ires);
	findatm_(" CZ ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__301, &i__1, ires);
	findatm_(" NH1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__302, &i__1, ires);
	findatm_(" NH2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__302, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 6;
	newatm_(&i__, &c__294, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 7;
	newatm_(&i__, &c__294, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_(" HG2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__296, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c__1);
	findatm_(" HG3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__296, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c_n1);
	findatm_(" HD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__298, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c__1);
	findatm_(" HD3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__298, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c_n1);
	findatm_(" HE ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__300, &i__1, &c_b63, &i__2, &c_b1207, &i__3, &c_b599,
		 &c__1);
	findatm_("HH11", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 11;
	newatm_(&i__, &c__303, &i__1, &c_b63, &i__2, &c_b1213, &i__3, &c_b76, 
		&c__0);
	findatm_("HH12", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 11;
	i__3 = atoms_1.n - 12;
	newatm_(&i__, &c__303, &i__1, &c_b63, &i__2, &c_b1219, &i__3, &c_b92, 
		&c__0);
	findatm_("HH21", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 12;
	i__3 = atoms_1.n - 13;
	newatm_(&i__, &c__303, &i__1, &c_b63, &i__2, &c_b1213, &i__3, &c_b76, 
		&c__0);
	findatm_("HH22", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 11;
	i__2 = atoms_1.n - 13;
	i__3 = atoms_1.n - 14;
	newatm_(&i__, &c__303, &i__1, &c_b63, &i__2, &c_b1219, &i__3, &c_b92, 
		&c__0);

/*     ornithine residue  (ORN) */

    } else if (s_cmp(resname, "ORN", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__310, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__312, &i__1, ires);
	findatm_(" CD ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__314, &i__1, ires);
	findatm_(" NE ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__316, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	newatm_(&i__, &c__311, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	newatm_(&i__, &c__311, &i__1, &c_b129, cai, &c_b298, &i__2, &c_b299, &
		c_n1);
	findatm_(" HG2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 4;
	newatm_(&i__, &c__313, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c__1);
	findatm_(" HG3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 5;
	newatm_(&i__, &c__313, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c_n1);
	findatm_(" HD2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 5;
	newatm_(&i__, &c__315, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c__1);
	findatm_(" HD3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 6;
	newatm_(&i__, &c__315, &i__1, &c_b129, &i__2, &c_b75, &i__3, &c_b75, &
		c_n1);
	findatm_(" HE1", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__317, &i__1, &c_b1137, &i__2, &c_b1138, &i__3, &
		c_b92, &c__0);
	findatm_(" HE2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 8;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 10;
	newatm_(&i__, &c__317, &i__1, &c_b1137, &i__2, &c_b1138, &i__3, &
		c_b97, &c__0);
	findatm_(" HE3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 10;
	i__3 = atoms_1.n - 11;
	newatm_(&i__, &c__317, &i__1, &c_b1137, &i__2, &c_b1138, &i__3, &
		c_b87, &c__0);

/*     methylalanine residue  (AIB) */

    } else if (s_cmp(resname, "AIB", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB1", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__323, cai, ires);
	findatm_(" CB2", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__323, cai, ires);
	findatm_("HB11", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	newatm_(&i__, &c__324, &i__1, &c_b129, cai, &c_b222, ni, &c_b92, &
		c__0);
	findatm_("HB12", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	newatm_(&i__, &c__324, &i__1, &c_b129, cai, &c_b222, ni, &c_b97, &
		c__0);
	findatm_("HB13", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	newatm_(&i__, &c__324, &i__1, &c_b129, cai, &c_b222, ni, &c_b87, &
		c__0);
	findatm_("HB21", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	newatm_(&i__, &c__324, &i__1, &c_b129, cai, &c_b222, ni, &c_b92, &
		c__0);
	findatm_("HB22", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	newatm_(&i__, &c__324, &i__1, &c_b129, cai, &c_b222, ni, &c_b97, &
		c__0);
	findatm_("HB23", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	newatm_(&i__, &c__324, &i__1, &c_b129, cai, &c_b222, ni, &c_b87, &
		c__0);

/*     pyroglutamic acid residue  (PCA) */

    } else if (s_cmp(resname, "PCA", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" CB ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__331, cai, ires);
	findatm_(" CG ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__333, &i__1, ires);
	findatm_(" CD ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__335, &i__1, ires);
	i__1 = atoms_1.n - 1;
	addbond_(&i__1, ni);
	findatm_(" OE ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__336, &i__1, ires);
	findatm_(" HB2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 3;
	newatm_(&i__, &c__332, &i__1, &c_b129, cai, &c_b524, &i__2, &c_b524, &
		c__1);
	findatm_(" HB3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 4;
	newatm_(&i__, &c__332, &i__1, &c_b129, cai, &c_b524, &i__2, &c_b524, &
		c_n1);
	findatm_(" HG2", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 4;
	newatm_(&i__, &c__334, &i__1, &c_b129, &i__2, &c_b524, &i__3, &c_b524,
		 &c__1);
	findatm_(" HG3", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 5;
	newatm_(&i__, &c__334, &i__1, &c_b129, &i__2, &c_b524, &i__3, &c_b524,
		 &c_n1);

/*     unknown residue  (UNK) */

    } else if (s_cmp(resname, "UNK", (ftnlen)3, (ftnlen)3) == 0) {
	if (*ires == 1) {
	    newatm_(&c__0, &c__355, cai, &c_b129, ni, &c_b75, ci, &c_b75, &
		    c__1);
	} else if (*ires == sequen_1.nseq) {
	    newatm_(&c__0, &c__506, cai, &c_b129, ni, &c_b75, ci, &c_b75, &
		    c__1);
	} else {
	    newatm_(&c__0, &c__6, cai, &c_b129, ni, &c_b75, ci, &c_b75, &c__1)
		    ;
	}
    }
    return 0;
} /* addside_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine ligase  --  coordinates from PDB nucleic acid  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "ligase" translates a nucleic acid structure in Protein Data */
/*     Bank format to a Cartesian coordinate file and sequence file */


/* Subroutine */ int ligase_(void)
{
    /* Initialized data */

    static integer h5ttyp[12] = { 1233,1233,1233,1233,1245,1245,1245,1245,0,0,
	    0,0 };
    static integer h3ttyp[12] = { 1238,1238,1238,1238,1250,1250,1250,1250,0,0,
	    0,0 };
    static integer o5typ[12] = { 1001,1031,1062,1090,1117,1146,1176,1203,0,0,
	    0,0 };
    static integer c5typ[12] = { 1002,1032,1063,1091,1118,1147,1177,1204,0,0,
	    0,0 };
    static integer h51typ[12] = { 1003,1033,1064,1092,1119,1148,1178,1205,0,0,
	    0,0 };
    static integer h52typ[12] = { 1004,1034,1065,1093,1120,1149,1179,1206,0,0,
	    0,0 };
    static integer c4typ[12] = { 1005,1035,1066,1094,1121,1150,1180,1207,0,0,
	    0,0 };
    static integer h4typ[12] = { 1006,1036,1067,1095,1122,1151,1181,1208,0,0,
	    0,0 };
    static integer o4typ[12] = { 1007,1037,1068,1096,1123,1152,1182,1209,0,0,
	    0,0 };
    static integer c1typ[12] = { 1008,1038,1069,1097,1124,1153,1183,1210,0,0,
	    0,0 };
    static integer h1typ[12] = { 1009,1039,1070,1098,1125,1154,1184,1211,0,0,
	    0,0 };
    static integer c3typ[12] = { 1010,1040,1071,1099,1126,1155,1185,1212,0,0,
	    0,0 };
    static integer h3typ[12] = { 1011,1041,1072,1100,1127,1156,1186,1213,0,0,
	    0,0 };
    static integer c2typ[12] = { 1012,1042,1073,1101,1128,1157,1187,1214,0,0,
	    0,0 };
    static integer h21typ[12] = { 1013,1043,1074,1102,1129,1158,1188,1215,0,0,
	    0,0 };
    static integer h22typ[12] = { 1015,1045,1076,1104,1130,1159,1189,1216,0,0,
	    0,0 };
    static integer o3typ[12] = { 1016,1046,1077,1105,1131,1160,1190,1217,0,0,
	    0,0 };
    static integer o2typ[12] = { 1014,1044,1075,1103,0,0,0,0,0,0,0,0 };
    static integer ptyp[12] = { 1230,1230,1230,1230,1242,1242,1242,1242,0,0,0,
	    0 };
    static integer optyp[12] = { 1231,1231,1231,1231,1243,1243,1243,1243,0,0,
	    0,0 };

    /* Format strings */
    static char fmt_20[] = "(/,\002 LIGASE  --  The Maximum of\002,i8,\002 A"
	    "toms\002,\002 has been Exceeded\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static logical endchain, newchain;
    static integer i__, j, k, c1i, c2i, c3i, c4i, c5i, o2i, o3i, o4i, o5i, 
	    poi, nres, ityp, stop;
    extern /* Subroutine */ int fatal_(void);
    static logical deoxy[10000];
    static integer start;
    extern /* Subroutine */ int oldatm_(integer *, integer *, integer *, 
	    integer *);
    static integer resatm[20000]	/* was [2][10000] */;
    extern /* Subroutine */ int newatm_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), addbase_(char *, integer *, integer *, integer *, 
	    integer *, ftnlen), addbond_(integer *, integer *), findatm_(char 
	    *, integer *, integer *, integer *, ftnlen);
    static char resname[3];

    /* Fortran I/O blocks */
    static cilist io___112 = { 0, 0, 0, fmt_20, 0 };



#define seq_ref(a_0,a_1) &sequen_1.seq[(a_1)*3 + a_0 - 3]
#define nuclz_ref(a_0,a_1) &resdue_1.nuclz[(a_1)*3 + a_0 - 3]
#define ichain_ref(a_1,a_2) sequen_1.ichain[(a_2)*2 + a_1 - 3]
#define resnam_ref(a_0,a_1) &pdb_1.resnam[(a_1)*3 + a_0 - 3]
#define resatm_ref(a_1,a_2) resatm[(a_2)*2 + a_1 - 3]
#define pdbtyp_ref(a_0,a_1) &pdb_1.pdbtyp[(a_1)*6 + a_0 - 6]



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

/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  pdb.i  --  definition of a Protein Data Bank structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xpdb      x-coordinate of each atom stored in PDB format */
/*     ypdb      y-coordinate of each atom stored in PDB format */
/*     zpdb      z-coordinate of each atom stored in PDB format */
/*     npdb      number of atoms stored in Protein Data Bank format */
/*     resnum    number of the residue to which each atom belongs */
/*     npdb12    number of atoms directly bonded to each CONECT atom */
/*     ipdb12    atom numbers of atoms connected to each CONECT atom */
/*     pdblist   list of the Protein Data Bank atom number of each atom */
/*     pdbtyp    Protein Data Bank record type assigned to each atom */
/*     atmnam    Protein Data Bank atom name assigned to each atom */
/*     resnam    Protein Data Bank residue name assigned to each atom */
/*     chntyp    string with PDB chain identifiers to be included */
/*     altsym    string with PDB alternate locations to be included */
/*     instyp    string with PDB insertion records to be included */




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



/*     biopolymer atom types for nucleic acid backbone atoms */



/*     set a pointer to the first and last atom of each residue */

    nres = 0;
    k = 0;
    i__1 = pdb_1.npdb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (s_cmp(pdbtyp_ref(0, i__), "ATOM  ", (ftnlen)6, (ftnlen)6) == 0) {
	    if (pdb_1.resnum[i__ - 1] != k) {
		k = pdb_1.resnum[i__ - 1];
		if (nres != 0) {
		    resatm_ref(2, nres) = i__ - 1;
		}
		++nres;
		resatm_ref(1, nres) = i__;
	    }
	}
    }
    if (nres != 0) {
	resatm_ref(2, nres) = pdb_1.npdb;
    }

/*     check for deoxyribose and change residue name if necessary */

    i__1 = nres;
    for (i__ = 1; i__ <= i__1; ++i__) {
	deoxy[i__ - 1] = FALSE_;
	start = resatm_ref(1, i__);
	stop = resatm_ref(2, i__);
	s_copy(resname, resnam_ref(0, start), (ftnlen)3, (ftnlen)3);
	findatm_(" O2'", &start, &stop, &k, (ftnlen)4);
	if (k == 0) {
	    deoxy[i__ - 1] = TRUE_;
	    i__2 = stop;
	    for (j = start; j <= i__2; ++j) {
		if (s_cmp(resname, "A  ", (ftnlen)3, (ftnlen)3) == 0) {
		    s_copy(resnam_ref(0, j), "DA ", (ftnlen)3, (ftnlen)3);
		}
		if (s_cmp(resname, "G  ", (ftnlen)3, (ftnlen)3) == 0) {
		    s_copy(resnam_ref(0, j), "DG ", (ftnlen)3, (ftnlen)3);
		}
		if (s_cmp(resname, "C  ", (ftnlen)3, (ftnlen)3) == 0) {
		    s_copy(resnam_ref(0, j), "DC ", (ftnlen)3, (ftnlen)3);
		}
		if (s_cmp(resname, "U  ", (ftnlen)3, (ftnlen)3) == 0) {
		    s_copy(resnam_ref(0, j), "DU ", (ftnlen)3, (ftnlen)3);
		}
		if (s_cmp(resname, "T  ", (ftnlen)3, (ftnlen)3) == 0) {
		    s_copy(resnam_ref(0, j), "DT ", (ftnlen)3, (ftnlen)3);
		}
	    }
	}
    }

/*     get the three-letter sequence and code for each residue */

    sequen_1.nseq = nres;
    i__1 = nres;
    for (i__ = 1; i__ <= i__1; ++i__) {
	start = resatm_ref(1, i__);
	s_copy(resname, resnam_ref(0, start), (ftnlen)3, (ftnlen)3);
	s_copy(seq_ref(0, i__), "UNK", (ftnlen)3, (ftnlen)3);
	sequen_1.seqtyp[i__ - 1] = 12;
	for (k = 1; k <= 12; ++k) {
	    if (s_cmp(resname, nuclz_ref(0, k), (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(seq_ref(0, i__), nuclz_ref(0, k), (ftnlen)3, (ftnlen)3)
			;
		sequen_1.seqtyp[i__ - 1] = k;
		goto L10;
	    }
	}
L10:
	;
    }

/*     set the current atom to be the first atom */

    atoms_1.n = 1;

/*     locate and assign the atoms that make up each residue */

    i__1 = nres;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ityp = sequen_1.seqtyp[i__ - 1];
	start = resatm_ref(1, i__);
	stop = resatm_ref(2, i__);
	s_copy(resname, resnam_ref(0, start), (ftnlen)3, (ftnlen)3);

/*     check that the maximum allowed atoms is not exceeded */

	if (atoms_1.n + 25 > 25000) {
	    io___112.ciunit = iounit_1.iout;
	    s_wsfe(&io___112);
	    do_fio(&c__1, (char *)&c__25000, (ftnlen)sizeof(integer));
	    e_wsfe();
	    fatal_();
	}

/*     test for initial or final residue of a nucleotide chain */

	newchain = FALSE_;
	endchain = FALSE_;
	i__2 = sequen_1.nchain;
	for (j = 1; j <= i__2; ++j) {
	    if (i__ == ichain_ref(1, j)) {
		newchain = TRUE_;
		poi = 0;
		o3i = 0;
	    }
	    if (i__ == ichain_ref(2, j)) {
		endchain = TRUE_;
	    }
	}

/*     build the phosphate atoms of the current residue */

	if (s_cmp(resname, "TP ", (ftnlen)3, (ftnlen)3) == 0) {
	} else if (s_cmp(resname, "DP ", (ftnlen)3, (ftnlen)3) == 0) {
	} else if (s_cmp(resname, "MP ", (ftnlen)3, (ftnlen)3) == 0) {
	} else if (! newchain) {
	    findatm_(" P  ", &start, &stop, &k, (ftnlen)4);
	    if (k != 0) {
		poi = atoms_1.n;
	    }
	    j = ptyp[ityp - 1];
	    oldatm_(&k, &j, &o3i, &i__);
	    findatm_(" OP1", &start, &stop, &k, (ftnlen)4);
	    j = optyp[ityp - 1];
	    i__2 = atoms_1.n - 1;
	    oldatm_(&k, &j, &i__2, &i__);
	    findatm_(" OP2", &start, &stop, &k, (ftnlen)4);
	    j = optyp[ityp - 1];
	    i__2 = atoms_1.n - 2;
	    oldatm_(&k, &j, &i__2, &i__);
	}

/*     build the ribose sugar atoms of the current residue */

	findatm_(" O5'", &start, &stop, &k, (ftnlen)4);
	if (k != 0) {
	    o5i = atoms_1.n;
	}
	j = o5typ[ityp - 1];
	if (newchain) {
	    if (deoxy[i__ - 1]) {
		j = 1244;
	    } else {
		j = 1232;
	    }
	}
	oldatm_(&k, &j, &poi, &i__);
	findatm_(" C5'", &start, &stop, &k, (ftnlen)4);
	if (k != 0) {
	    c5i = atoms_1.n;
	}
	j = c5typ[ityp - 1];
	i__2 = atoms_1.n - 1;
	oldatm_(&k, &j, &i__2, &i__);
	findatm_(" C4'", &start, &stop, &k, (ftnlen)4);
	if (k != 0) {
	    c4i = atoms_1.n;
	}
	j = c4typ[ityp - 1];
	i__2 = atoms_1.n - 1;
	oldatm_(&k, &j, &i__2, &i__);
	findatm_(" O4'", &start, &stop, &k, (ftnlen)4);
	if (k != 0) {
	    o4i = atoms_1.n;
	}
	j = o4typ[ityp - 1];
	i__2 = atoms_1.n - 1;
	oldatm_(&k, &j, &i__2, &i__);
	findatm_(" C1'", &start, &stop, &k, (ftnlen)4);
	if (k != 0) {
	    c1i = atoms_1.n;
	}
	j = c1typ[ityp - 1];
	i__2 = atoms_1.n - 1;
	oldatm_(&k, &j, &i__2, &i__);
	findatm_(" C3'", &start, &stop, &k, (ftnlen)4);
	if (k != 0) {
	    c3i = atoms_1.n;
	}
	j = c3typ[ityp - 1];
	i__2 = atoms_1.n - 3;
	oldatm_(&k, &j, &i__2, &i__);
	findatm_(" C2'", &start, &stop, &k, (ftnlen)4);
	if (k != 0) {
	    c2i = atoms_1.n;
	}
	j = c2typ[ityp - 1];
	i__2 = atoms_1.n - 1;
	oldatm_(&k, &j, &i__2, &i__);
	i__2 = atoms_1.n - 1;
	i__3 = atoms_1.n - 3;
	addbond_(&i__2, &i__3);
	findatm_(" O3'", &start, &stop, &k, (ftnlen)4);
	if (k != 0) {
	    o3i = atoms_1.n;
	}
	j = o3typ[ityp - 1];
	if (endchain) {
	    if (deoxy[i__ - 1]) {
		j = 1249;
	    } else {
		j = 1237;
	    }
	}
	i__2 = atoms_1.n - 2;
	oldatm_(&k, &j, &i__2, &i__);
	if (! deoxy[i__ - 1]) {
	    findatm_(" O2'", &start, &stop, &k, (ftnlen)4);
	    if (k != 0) {
		o2i = atoms_1.n;
	    }
	    j = o2typ[ityp - 1];
	    i__2 = atoms_1.n - 2;
	    oldatm_(&k, &j, &i__2, &i__);
	}

/*     build the hydrogen atoms of the current residue */

	if (newchain) {
	    findatm_(" H5T", &start, &stop, &k, (ftnlen)4);
	    j = h5ttyp[ityp - 1];
	    newatm_(&k, &j, &o5i, &c_b1425, &c5i, &c_b75, &c4i, &c_b92, &c__0)
		    ;
	}
	findatm_(" H5'", &start, &stop, &k, (ftnlen)4);
	j = h51typ[ityp - 1];
	newatm_(&k, &j, &c5i, &c_b598, &o5i, &c_b75, &c4i, &c_b75, &c__1);
	findatm_("H5''", &start, &stop, &k, (ftnlen)4);
	j = h52typ[ityp - 1];
	newatm_(&k, &j, &c5i, &c_b598, &o5i, &c_b75, &c4i, &c_b75, &c_n1);
	findatm_(" H4'", &start, &stop, &k, (ftnlen)4);
	j = h4typ[ityp - 1];
	newatm_(&k, &j, &c4i, &c_b598, &c5i, &c_b75, &c3i, &c_b75, &c_n1);
	findatm_(" H3'", &start, &stop, &k, (ftnlen)4);
	j = h3typ[ityp - 1];
	newatm_(&k, &j, &c3i, &c_b598, &c4i, &c_b75, &c2i, &c_b75, &c_n1);
	if (deoxy[i__ - 1]) {
	    findatm_(" H2'", &start, &stop, &k, (ftnlen)4);
	    j = h21typ[ityp - 1];
	    newatm_(&k, &j, &c2i, &c_b598, &c3i, &c_b75, &c1i, &c_b75, &c_n1);
	    findatm_("H2''", &start, &stop, &k, (ftnlen)4);
	    j = h22typ[ityp - 1];
	    newatm_(&k, &j, &c2i, &c_b598, &c3i, &c_b75, &c1i, &c_b75, &c__1);
	} else {
	    findatm_(" H2'", &start, &stop, &k, (ftnlen)4);
	    j = h21typ[ityp - 1];
	    newatm_(&k, &j, &c2i, &c_b598, &c3i, &c_b75, &c1i, &c_b75, &c_n1);
	    findatm_("HO2'", &start, &stop, &k, (ftnlen)4);
	    j = h22typ[ityp - 1];
	    newatm_(&k, &j, &o2i, &c_b1425, &c2i, &c_b75, &c3i, &c_b92, &c__0)
		    ;
	}
	findatm_(" H1'", &start, &stop, &k, (ftnlen)4);
	j = h1typ[ityp - 1];
	newatm_(&k, &j, &c1i, &c_b598, &o4i, &c_b75, &c2i, &c_b75, &c_n1);
	if (endchain) {
	    findatm_(" H3T", &start, &stop, &k, (ftnlen)4);
	    j = h3ttyp[ityp - 1];
	    newatm_(&k, &j, &o3i, &c_b1425, &c3i, &c_b75, &c4i, &c_b92, &c__0)
		    ;
	}

/*     build the standard base atoms of the current residue */

	addbase_(resname, &i__, &start, &stop, &c1i, (ftnlen)3);
    }

/*     total number of atoms is one less than the current atom */

    --atoms_1.n;
    return 0;
} /* ligase_ */

#undef pdbtyp_ref
#undef resatm_ref
#undef resnam_ref
#undef ichain_ref
#undef nuclz_ref
#undef seq_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine addbase  --  build a single nucleic acid base  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "addbase" builds the Cartesian coordinates for a single nucleic */
/*     acid base; coordinates are read from the Protein Data Bank file */
/*     or found from internal coordinates, then atom types are assigned */
/*     and connectivity data generated */


/* Subroutine */ int addbase_(char *resname, integer *ires, integer *start, 
	integer *stop, integer *c1i, ftnlen resname_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int oldatm_(integer *, integer *, integer *, 
	    integer *), newatm_(integer *, integer *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *), 
	    addbond_(integer *, integer *), findatm_(char *, integer *, 
	    integer *, integer *, ftnlen);



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




/*     adenine in adenosine residue  (A) */

    if (s_cmp(resname, "A ", (ftnlen)3, (ftnlen)2) == 0) {
	findatm_(" N9 ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__1017, c1i, ires);
	findatm_(" C8 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1021, &i__1, ires);
	findatm_(" N7 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1020, &i__1, ires);
	findatm_(" C5 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1019, &i__1, ires);
	findatm_(" C6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1025, &i__1, ires);
	findatm_(" N6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1027, &i__1, ires);
	findatm_(" N1 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1024, &i__1, ires);
	findatm_(" C2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1023, &i__1, ires);
	findatm_(" N3 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1022, &i__1, ires);
	findatm_(" C4 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1018, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 7;
	addbond_(&i__1, &i__2);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 10;
	addbond_(&i__1, &i__2);
	findatm_(" H8 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__1030, &i__1, &c_b1503, &i__2, &c_b1504, &i__3, &
		c_b92, &c__0);
	findatm_(" H61", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__1028, &i__1, &c_b1425, &i__2, &c_b599, &i__3, &
		c_b92, &c__0);
	findatm_(" H62", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__1029, &i__1, &c_b1425, &i__2, &c_b599, &i__3, &
		c_b76, &c__0);
	findatm_(" H2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 4;
	newatm_(&i__, &c__1026, &i__1, &c_b1503, &i__2, &c_b1522, &i__3, &
		c_b92, &c__0);

/*     guanine in guanosine residue  (G) */

    } else if (s_cmp(resname, "G  ", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" N9 ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__1047, c1i, ires);
	findatm_(" C8 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1051, &i__1, ires);
	findatm_(" N7 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1050, &i__1, ires);
	findatm_(" C5 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1049, &i__1, ires);
	findatm_(" C6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1055, &i__1, ires);
	findatm_(" O6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1060, &i__1, ires);
	findatm_(" N1 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1054, &i__1, ires);
	findatm_(" C2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1053, &i__1, ires);
	findatm_(" N2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1057, &i__1, ires);
	findatm_(" N3 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1052, &i__1, ires);
	findatm_(" C4 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1048, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 8;
	addbond_(&i__1, &i__2);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 11;
	addbond_(&i__1, &i__2);
	findatm_(" H8 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__1061, &i__1, &c_b1503, &i__2, &c_b1551, &i__3, &
		c_b92, &c__0);
	findatm_(" H1 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__1056, &i__1, &c_b1425, &i__2, &c_b1557, &i__3, &
		c_b92, &c__0);
	findatm_(" H21", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__1058, &i__1, &c_b1425, &i__2, &c_b599, &i__3, &
		c_b76, &c__0);
	findatm_(" H22", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__1059, &i__1, &c_b1425, &i__2, &c_b599, &i__3, &
		c_b92, &c__0);

/*     cytosine in cytidine residue  (C) */

    } else if (s_cmp(resname, "C  ", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" N1 ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__1078, c1i, ires);
	findatm_(" C2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1079, &i__1, ires);
	findatm_(" O2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1084, &i__1, ires);
	findatm_(" N3 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1080, &i__1, ires);
	findatm_(" C4 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1081, &i__1, ires);
	findatm_(" N4 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1085, &i__1, ires);
	findatm_(" C5 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1082, &i__1, ires);
	findatm_(" C6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1083, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 8;
	addbond_(&i__1, &i__2);
	findatm_(" H41", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 5;
	newatm_(&i__, &c__1086, &i__1, &c_b1425, &i__2, &c_b599, &i__3, &
		c_b76, &c__0);
	findatm_(" H42", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 6;
	newatm_(&i__, &c__1087, &i__1, &c_b1425, &i__2, &c_b599, &i__3, &
		c_b92, &c__0);
	findatm_(" H5 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__1088, &i__1, &c_b1503, &i__2, &c_b1604, &i__3, &
		c_b92, &c__0);
	findatm_(" H6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__1089, &i__1, &c_b1503, &i__2, &c_b1610, &i__3, &
		c_b92, &c__0);

/*     uracil in uridine residue  (U) */

    } else if (s_cmp(resname, "U  ", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" N1 ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__1106, c1i, ires);
	findatm_(" C2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1107, &i__1, ires);
	findatm_(" O2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1112, &i__1, ires);
	findatm_(" N3 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1108, &i__1, ires);
	findatm_(" C4 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1109, &i__1, ires);
	findatm_(" O4 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1114, &i__1, ires);
	findatm_(" C5 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1110, &i__1, ires);
	findatm_(" C6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1111, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 8;
	addbond_(&i__1, &i__2);
	findatm_(" H3 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__1113, &i__1, &c_b1425, &i__2, &c_b1633, &i__3, &
		c_b92, &c__0);
	findatm_(" H5 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 6;
	newatm_(&i__, &c__1115, &i__1, &c_b1503, &i__2, &c_b1639, &i__3, &
		c_b92, &c__0);
	findatm_(" H6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 6;
	newatm_(&i__, &c__1116, &i__1, &c_b1503, &i__2, &c_b1645, &i__3, &
		c_b92, &c__0);

/*     adenine in deoxyadenosine residue  (DA) */

    } else if (s_cmp(resname, "DA ", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" N9 ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__1132, c1i, ires);
	findatm_(" C8 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1136, &i__1, ires);
	findatm_(" N7 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1135, &i__1, ires);
	findatm_(" C5 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1134, &i__1, ires);
	findatm_(" C6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1140, &i__1, ires);
	findatm_(" N6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1142, &i__1, ires);
	findatm_(" N1 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1139, &i__1, ires);
	findatm_(" C2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1138, &i__1, ires);
	findatm_(" N3 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1137, &i__1, ires);
	findatm_(" C4 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1133, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 7;
	addbond_(&i__1, &i__2);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 10;
	addbond_(&i__1, &i__2);
	findatm_(" H8 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 9;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__1145, &i__1, &c_b1503, &i__2, &c_b1504, &i__3, &
		c_b92, &c__0);
	findatm_(" H61", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__1143, &i__1, &c_b1425, &i__2, &c_b599, &i__3, &
		c_b92, &c__0);
	findatm_(" H62", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 7;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__1144, &i__1, &c_b1425, &i__2, &c_b599, &i__3, &
		c_b76, &c__0);
	findatm_(" H2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 4;
	newatm_(&i__, &c__1141, &i__1, &c_b1503, &i__2, &c_b1522, &i__3, &
		c_b92, &c__0);

/*     guanine in deoxyguanosine residue  (DG) */

    } else if (s_cmp(resname, "DG ", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" N9 ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__1161, c1i, ires);
	findatm_(" C8 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1165, &i__1, ires);
	findatm_(" N7 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1164, &i__1, ires);
	findatm_(" C5 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1163, &i__1, ires);
	findatm_(" C6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1169, &i__1, ires);
	findatm_(" O6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1174, &i__1, ires);
	findatm_(" N1 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1168, &i__1, ires);
	findatm_(" C2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1167, &i__1, ires);
	findatm_(" N2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1171, &i__1, ires);
	findatm_(" N3 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1166, &i__1, ires);
	findatm_(" C4 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1162, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 8;
	addbond_(&i__1, &i__2);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 11;
	addbond_(&i__1, &i__2);
	findatm_(" H8 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 10;
	i__2 = atoms_1.n - 9;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__1175, &i__1, &c_b1503, &i__2, &c_b1551, &i__3, &
		c_b92, &c__0);
	findatm_(" H1 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__1170, &i__1, &c_b1425, &i__2, &c_b1557, &i__3, &
		c_b92, &c__0);
	findatm_(" H21", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__1172, &i__1, &c_b1425, &i__2, &c_b599, &i__3, &
		c_b76, &c__0);
	findatm_(" H22", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 8;
	newatm_(&i__, &c__1173, &i__1, &c_b1425, &i__2, &c_b599, &i__3, &
		c_b92, &c__0);

/*     cytosine in deoxycytidine residue  (DC) */

    } else if (s_cmp(resname, "DC ", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" N1 ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__1191, c1i, ires);
	findatm_(" C2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1192, &i__1, ires);
	findatm_(" O2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1197, &i__1, ires);
	findatm_(" N3 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1193, &i__1, ires);
	findatm_(" C4 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1194, &i__1, ires);
	findatm_(" N4 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1198, &i__1, ires);
	findatm_(" C5 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1195, &i__1, ires);
	findatm_(" C6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1196, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 8;
	addbond_(&i__1, &i__2);
	findatm_(" H41", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 5;
	newatm_(&i__, &c__1199, &i__1, &c_b1425, &i__2, &c_b599, &i__3, &
		c_b76, &c__0);
	findatm_(" H42", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 6;
	newatm_(&i__, &c__1200, &i__1, &c_b1425, &i__2, &c_b599, &i__3, &
		c_b92, &c__0);
	findatm_(" H5 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__1201, &i__1, &c_b1503, &i__2, &c_b1604, &i__3, &
		c_b92, &c__0);
	findatm_(" H6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 7;
	newatm_(&i__, &c__1202, &i__1, &c_b1503, &i__2, &c_b1610, &i__3, &
		c_b92, &c__0);

/*     thymine in deoxythymidine residue  (DT) */

    } else if (s_cmp(resname, "DT ", (ftnlen)3, (ftnlen)3) == 0) {
	findatm_(" N1 ", start, stop, &i__, (ftnlen)4);
	oldatm_(&i__, &c__1218, c1i, ires);
	findatm_(" C2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1219, &i__1, ires);
	findatm_(" O2 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1224, &i__1, ires);
	findatm_(" N3 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1220, &i__1, ires);
	findatm_(" C4 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1221, &i__1, ires);
	findatm_(" O4 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1226, &i__1, ires);
	findatm_(" C5 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1222, &i__1, ires);
	findatm_(" C7 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 1;
	oldatm_(&i__, &c__1227, &i__1, ires);
	findatm_(" C6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 2;
	oldatm_(&i__, &c__1223, &i__1, ires);
	i__1 = atoms_1.n - 1;
	i__2 = atoms_1.n - 9;
	addbond_(&i__1, &i__2);
	findatm_(" H3 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 6;
	i__2 = atoms_1.n - 8;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__1225, &i__1, &c_b1425, &i__2, &c_b1803, &i__3, &
		c_b92, &c__0);
	findatm_(" H71", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 3;
	i__2 = atoms_1.n - 4;
	i__3 = atoms_1.n - 6;
	newatm_(&i__, &c__1228, &i__1, &c_b598, &i__2, &c_b75, &i__3, &c_b76, 
		&c__0);
	findatm_(" H72", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 4;
	i__2 = atoms_1.n - 5;
	i__3 = atoms_1.n - 1;
	newatm_(&i__, &c__1228, &i__1, &c_b598, &i__2, &c_b75, &i__3, &c_b75, 
		&c__1);
	findatm_(" H73", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 6;
	i__3 = atoms_1.n - 2;
	newatm_(&i__, &c__1228, &i__1, &c_b598, &i__2, &c_b75, &i__3, &c_b75, 
		&c_n1);
	findatm_(" H6 ", start, stop, &i__, (ftnlen)4);
	i__1 = atoms_1.n - 5;
	i__2 = atoms_1.n - 7;
	i__3 = atoms_1.n - 9;
	newatm_(&i__, &c__1229, &i__1, &c_b1503, &i__2, &c_b1610, &i__3, &
		c_b92, &c__0);
    }
    return 0;
} /* addbase_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine hetatom  --  coordinates of PDB water and ions  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "hetatom" translates water molecules and ions in Protein Data */
/*     Bank format to a Cartesian coordinate file and sequence file */


/* Subroutine */ int hetatom_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int oldatm_(integer *, integer *, integer *, 
	    integer *), newatm_(integer *, integer *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *);


#define atmnam_ref(a_0,a_1) &pdb_1.atmnam[(a_1)*4 + a_0 - 4]
#define resnam_ref(a_0,a_1) &pdb_1.resnam[(a_1)*3 + a_0 - 3]
#define pdbtyp_ref(a_0,a_1) &pdb_1.pdbtyp[(a_1)*6 + a_0 - 6]



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
/*     ##  pdb.i  --  definition of a Protein Data Bank structure  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     xpdb      x-coordinate of each atom stored in PDB format */
/*     ypdb      y-coordinate of each atom stored in PDB format */
/*     zpdb      z-coordinate of each atom stored in PDB format */
/*     npdb      number of atoms stored in Protein Data Bank format */
/*     resnum    number of the residue to which each atom belongs */
/*     npdb12    number of atoms directly bonded to each CONECT atom */
/*     ipdb12    atom numbers of atoms connected to each CONECT atom */
/*     pdblist   list of the Protein Data Bank atom number of each atom */
/*     pdbtyp    Protein Data Bank record type assigned to each atom */
/*     atmnam    Protein Data Bank atom name assigned to each atom */
/*     resnam    Protein Data Bank residue name assigned to each atom */
/*     chntyp    string with PDB chain identifiers to be included */
/*     altsym    string with PDB alternate locations to be included */
/*     instyp    string with PDB insertion records to be included */




/*     find water molecules and ions in PDB HETATM records */

    ++atoms_1.n;
    i__ = 0;
    while(i__ < pdb_1.npdb) {
	++i__;
	if (s_cmp(pdbtyp_ref(0, i__), "HETATM", (ftnlen)6, (ftnlen)6) == 0) {
	    if (s_cmp(resnam_ref(0, i__), "HOH", (ftnlen)3, (ftnlen)3) == 0) {
		if (s_cmp(atmnam_ref(0, i__), " O  ", (ftnlen)4, (ftnlen)4) ==
			 0) {
		    oldatm_(&i__, &c__2001, &c__0, &c__0);
		    if (s_cmp(atmnam_ref(0, i__ + 1), " H  ", (ftnlen)4, (
			    ftnlen)4) == 0 && s_cmp(atmnam_ref(0, i__ + 2), 
			    " H  ", (ftnlen)4, (ftnlen)4) == 0) {
			i__1 = i__ + 1;
			i__2 = atoms_1.n - 1;
			oldatm_(&i__1, &c__2002, &i__2, &c__0);
			i__1 = i__ + 2;
			i__2 = atoms_1.n - 2;
			oldatm_(&i__1, &c__2002, &i__2, &c__0);
			i__ += 2;
		    } else {
			i__1 = atoms_1.n - 1;
			i__2 = atoms_1.n - 2;
			i__3 = atoms_1.n - 3;
			newatm_(&c__0, &c__2002, &i__1, &c_b1845, &i__2, &
				c_b75, &i__3, &c_b599, &c__0);
			i__1 = atoms_1.n - 2;
			i__2 = atoms_1.n - 1;
			i__3 = atoms_1.n - 3;
			newatm_(&c__0, &c__2002, &i__1, &c_b1845, &i__2, &
				c_b75, &i__3, &c_b599, &c__0);
		    }
		}
	    } else if (s_cmp(resnam_ref(0, i__), "NA ", (ftnlen)3, (ftnlen)3) 
		    == 0) {
		oldatm_(&i__, &c__2003, &c__0, &c__0);
	    } else if (s_cmp(resnam_ref(0, i__), "K  ", (ftnlen)3, (ftnlen)3) 
		    == 0) {
		oldatm_(&i__, &c__2004, &c__0, &c__0);
	    } else if (s_cmp(resnam_ref(0, i__), "MG ", (ftnlen)3, (ftnlen)3) 
		    == 0) {
		oldatm_(&i__, &c__2005, &c__0, &c__0);
	    } else if (s_cmp(resnam_ref(0, i__), "CA ", (ftnlen)3, (ftnlen)3) 
		    == 0) {
		oldatm_(&i__, &c__2006, &c__0, &c__0);
	    } else if (s_cmp(resnam_ref(0, i__), "CL ", (ftnlen)3, (ftnlen)3) 
		    == 0) {
		oldatm_(&i__, &c__2007, &c__0, &c__0);
	    }
	}
    }
    --atoms_1.n;
    return 0;
} /* hetatom_ */

#undef pdbtyp_ref
#undef resnam_ref
#undef atmnam_ref


/* Main program alias */ int pdbxyz_ () { MAIN__ (); return 0; }
