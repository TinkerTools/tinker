/* makepdb.f -- translated by f2c (version 20050501).
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
    doublereal molmass[25000], totmass;
    integer nmol, kmol[25000], imol[50000]	/* was [2][25000] */, molcule[
	    25000];
} molcul_;

#define molcul_1 molcul_

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
    integer nseq, nchain, ichain[20000]	/* was [2][10000] */, seqtyp[10000];
    char seq[30000], chnnam[10000];
} sequen_;

#define sequen_1 sequen_

/* Table of constant values */

static integer c__2 = 2;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine makepdb  --  convert Cartesian to PDB format  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "makexyz" converts a set of Cartesian coordinates to Protein */
/*     Data Bank format with special handling for systems consisting */
/*     of polypeptide chains, ligands and water molecules */


/* Subroutine */ int makepdb_(void)
{
    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2], i__3, i__4;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer f_inqu(inlist *), f_open(olist *), f_rew(alist *), f_clos(cllist *
	    );
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern integer freeunit_(void);
    static integer i__, j, k, m, ka, ci[10000], ni[10000], oi[10000], kn, kp, 
	    c1i[10000], c2i[10000], c3i[10000], c4i[10000], c5i[10000], o2i[
	    10000], o3i[10000], o4i[10000], o5i[10000], op1[10000], op2[10000]
	    , op3[10000], cai[10000], cbi, poi[10000], iseq, stop, noxy;
    static logical cbone, nbone, obone, water[25000], exist;
    static integer start;
    extern /* Subroutine */ int attach_(void), pdbatm_(char *, char *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer pdbnum;
    static logical hetmol[25000];
    static integer atmnum, nhydro;
    static char moltyp[7];
    extern /* Subroutine */ int getbase_(char *, integer *, integer *, ftnlen)
	    ;
    static char atmname[4];
    extern /* Subroutine */ int readseq_(integer *), getside_(char *, integer 
	    *, integer *, integer *, integer *, ftnlen);
    static char seqfile[120], resname[3];
    extern /* Subroutine */ int getnuch_(char *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, ftnlen), numeral_(integer *, char *, 
	    integer *, ftnlen), getproh_(char *, integer *, integer *, 
	    integer *, integer *, integer *, ftnlen), version_(char *, char *,
	     ftnlen, ftnlen);
    static integer justify;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
#define i13_ref(a_1,a_2) couple_1.i13[(a_2)*24 + a_1 - 25]
#define i14_ref(a_1,a_2) couple_1.i14[(a_2)*72 + a_1 - 73]
#define seq_ref(a_0,a_1) &sequen_1.seq[(a_1)*3 + a_0 - 3]
#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define imol_ref(a_1,a_2) molcul_1.imol[(a_2)*2 + a_1 - 3]
#define ipdb12_ref(a_1,a_2) pdb_1.ipdb12[(a_2)*8 + a_1 - 9]
#define amino_ref(a_0,a_1) &resdue_1.amino[(a_1)*3 + a_0 - 3]
#define nuclz_ref(a_0,a_1) &resdue_1.nuclz[(a_1)*3 + a_0 - 3]
#define ichain_ref(a_1,a_2) sequen_1.ichain[(a_2)*2 + a_1 - 3]
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




/*     initialize the mapping between TINKER and PDB atoms */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pdb_1.pdblist[i__ - 1] = 0;
    }

/*     read the biopolymer sequence file if one exists */

    iseq = freeunit_();
/* Writing concatenation */
    i__2[0] = files_1.leng, a__1[0] = files_1.filename;
    i__2[1] = 4, a__1[1] = ".seq";
    s_cat(seqfile, a__1, i__2, &c__2, (ftnlen)120);
    version_(seqfile, "old", (ftnlen)120, (ftnlen)3);
    ioin__1.inerr = 0;
    ioin__1.infilen = 120;
    ioin__1.infile = seqfile;
    ioin__1.inex = &exist;
    ioin__1.inopen = 0;
    ioin__1.innum = 0;
    ioin__1.innamed = 0;
    ioin__1.inname = 0;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    if (exist) {
	o__1.oerr = 0;
	o__1.ounit = iseq;
	o__1.ofnmlen = 120;
	o__1.ofnm = seqfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	al__1.aerr = 0;
	al__1.aunit = iseq;
	f_rew(&al__1);
	readseq_(&iseq);
	cl__1.cerr = 0;
	cl__1.cunit = iseq;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/*     assign the molecule type based on sequence information */

    s_copy(moltyp, "GENERIC", (ftnlen)7, (ftnlen)7);
    i__1 = sequen_1.nseq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(resname, seq_ref(0, i__), (ftnlen)3, (ftnlen)3);
	for (j = 1; j <= 31; ++j) {
	    if (s_cmp(resname, amino_ref(0, j), (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(moltyp, "PEPTIDE", (ftnlen)7, (ftnlen)7);
		goto L10;
	    }
	}
	s_copy(moltyp, "GENERIC", (ftnlen)7, (ftnlen)7);
L10:
	;
    }
    if (s_cmp(moltyp, "GENERIC", (ftnlen)7, (ftnlen)7) == 0) {
	i__1 = sequen_1.nseq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s_copy(resname, seq_ref(0, i__), (ftnlen)3, (ftnlen)3);
	    for (j = 1; j <= 12; ++j) {
		if (s_cmp(resname, nuclz_ref(0, j), (ftnlen)3, (ftnlen)3) == 
			0) {
		    s_copy(moltyp, "NUCACID", (ftnlen)7, (ftnlen)7);
		    goto L20;
		}
	    }
	    s_copy(moltyp, "GENERIC", (ftnlen)7, (ftnlen)7);
L20:
	    ;
	}
    }

/*     zero out the backbone atoms for biopolymer sequences */

    if (s_cmp(moltyp, "PEPTIDE", (ftnlen)7, (ftnlen)7) == 0) {
	i__1 = sequen_1.nseq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ni[i__ - 1] = 0;
	    cai[i__ - 1] = 0;
	    ci[i__ - 1] = 0;
	    oi[i__ - 1] = 0;
	}
    }
    if (s_cmp(moltyp, "NUCACID", (ftnlen)7, (ftnlen)7) == 0) {
	i__1 = sequen_1.nseq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    poi[i__ - 1] = 0;
	    op1[i__ - 1] = 0;
	    op2[i__ - 1] = 0;
	    op3[i__ - 1] = 0;
	    c5i[i__ - 1] = 0;
	    o5i[i__ - 1] = 0;
	    c4i[i__ - 1] = 0;
	    o4i[i__ - 1] = 0;
	    c3i[i__ - 1] = 0;
	    o3i[i__ - 1] = 0;
	    c2i[i__ - 1] = 0;
	    o2i[i__ - 1] = 0;
	    c1i[i__ - 1] = 0;
	}
    }

/*     check each molecule to see if it is a water molecule */

    i__1 = molcul_1.nmol;
    for (i__ = 1; i__ <= i__1; ++i__) {
	water[i__ - 1] = FALSE_;
	if (imol_ref(2, i__) - imol_ref(1, i__) == 2) {
	    noxy = 0;
	    nhydro = 0;
	    i__3 = imol_ref(2, i__);
	    for (j = imol_ref(1, i__); j <= i__3; ++j) {
		k = molcul_1.kmol[j - 1];
		if (atmtyp_1.atomic[k - 1] == 8) {
		    ++noxy;
		}
		if (atmtyp_1.atomic[k - 1] == 1) {
		    ++nhydro;
		}
	    }
	    if (noxy == 1 && nhydro == 2) {
		water[i__ - 1] = TRUE_;
	    }
	}
    }

/*     for general structures, transfer each atom to PDB format */

    if (s_cmp(moltyp, "GENERIC", (ftnlen)7, (ftnlen)7) == 0) {
	i__1 = molcul_1.nmol;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = imol_ref(2, i__);
	    for (j = imol_ref(1, i__); j <= i__3; ++j) {
		k = molcul_1.kmol[j - 1];
/* Writing concatenation */
		i__2[0] = 1, a__1[0] = " ";
		i__2[1] = 3, a__1[1] = name___ref(0, k);
		s_cat(atmname, a__1, i__2, &c__2, (ftnlen)4);
		if (water[i__ - 1]) {
		    s_copy(resname, "HOH", (ftnlen)3, (ftnlen)3);
		} else {
		    justify = 0;
		    numeral_(&atoms_1.type__[k - 1], resname, &justify, (
			    ftnlen)3);
		}
		pdbnum = i__;
		pdbatm_(atmname, resname, &pdbnum, &k, (ftnlen)4, (ftnlen)3);
		s_copy(pdbtyp_ref(0, pdb_1.npdb), "HETATM", (ftnlen)6, (
			ftnlen)6);
	    }
	}
	i__1 = molcul_1.nmol;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = imol_ref(2, i__);
	    for (j = imol_ref(1, i__); j <= i__3; ++j) {
		k = molcul_1.kmol[j - 1];
		kp = pdb_1.pdblist[k - 1];
		pdb_1.npdb12[kp - 1] = couple_1.n12[k - 1];
		i__4 = couple_1.n12[k - 1];
		for (m = 1; m <= i__4; ++m) {
		    ipdb12_ref(m, kp) = pdb_1.pdblist[i12_ref(m, k) - 1];
		}
	    }
	}
    }

/*     find the amide nitrogens and other peptide backbone atoms */

    if (s_cmp(moltyp, "PEPTIDE", (ftnlen)7, (ftnlen)7) == 0) {
	attach_();
	m = 1;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s_copy(resname, amino_ref(0, sequen_1.seqtyp[m - 1]), (ftnlen)3, (
		    ftnlen)3);
	    if (s_cmp(resname, "FOR", (ftnlen)3, (ftnlen)3) == 0) {
		if (atmtyp_1.atomic[i__ - 1] == 6) {
		    nbone = FALSE_;
		    obone = FALSE_;
		    i__3 = couple_1.n12[i__ - 1];
		    for (j = 1; j <= i__3; ++j) {
			k = i12_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 7) {
			    nbone = TRUE_;
			} else if (atmtyp_1.atomic[k - 1] == 8) {
			    obone = TRUE_;
			}
		    }
		    if (nbone && obone) {
			cai[m - 1] = i__;
			ci[m - 1] = i__;
			oi[m - 1] = i__ + 1;
			++m;
		    }
		}
	    } else if (s_cmp(resname, "ACE", (ftnlen)3, (ftnlen)3) == 0) {
		if (atmtyp_1.atomic[i__ - 1] == 6) {
		    nbone = FALSE_;
		    obone = FALSE_;
		    i__3 = couple_1.n13[i__ - 1];
		    for (j = 1; j <= i__3; ++j) {
			k = i13_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 7) {
			    nbone = TRUE_;
			} else if (atmtyp_1.atomic[k - 1] == 8) {
			    obone = TRUE_;
			}
		    }
		    if (nbone && obone) {
			cai[m - 1] = i__;
			ci[m - 1] = i__ + 1;
			oi[m - 1] = i__ + 2;
			++m;
		    }
		}
	    } else if (s_cmp(resname, "NH2", (ftnlen)3, (ftnlen)3) == 0) {
		if (atmtyp_1.atomic[i__ - 1] == 7) {
		    nbone = FALSE_;
		    obone = FALSE_;
		    i__3 = couple_1.n13[i__ - 1];
		    for (j = 1; j <= i__3; ++j) {
			k = i13_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 8) {
			    obone = TRUE_;
			}
		    }
		    i__3 = couple_1.n14[i__ - 1];
		    for (j = 1; j <= i__3; ++j) {
			k = i14_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 7) {
			    nbone = TRUE_;
			}
		    }
		    if (nbone && obone) {
			ni[m - 1] = i__;
			++m;
		    }
		}
	    } else if (s_cmp(resname, "NME", (ftnlen)3, (ftnlen)3) == 0) {
		if (atmtyp_1.atomic[i__ - 1] == 7) {
		    nbone = FALSE_;
		    obone = FALSE_;
		    i__3 = couple_1.n13[i__ - 1];
		    for (j = 1; j <= i__3; ++j) {
			k = i13_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 8) {
			    obone = TRUE_;
			}
		    }
		    i__3 = couple_1.n14[i__ - 1];
		    for (j = 1; j <= i__3; ++j) {
			k = i14_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 7) {
			    nbone = TRUE_;
			}
		    }
		    if (nbone && obone) {
			ni[m - 1] = i__;
			cai[m - 1] = i__ + 1;
			++m;
		    }
		}
	    } else {
		if (atmtyp_1.atomic[i__ - 1] == 7) {
		    obone = FALSE_;
		    i__3 = couple_1.n14[i__ - 1];
		    for (j = 1; j <= i__3; ++j) {
			k = i14_ref(j, i__);
			if (atmtyp_1.atomic[k - 1] == 8) {
			    obone = TRUE_;
			}
		    }
		    if (obone) {
			ni[m - 1] = i__;
			cai[m - 1] = i__ + 1;
			ci[m - 1] = i__ + 2;
			oi[m - 1] = i__ + 3;
			++m;
		    }
		}
	    }
	    if (m > sequen_1.nseq) {
		goto L30;
	    }
	}
L30:
	;
    }

/*     get all the atoms for each peptide residue in order */

    if (s_cmp(moltyp, "PEPTIDE", (ftnlen)7, (ftnlen)7) == 0) {
	pdb_1.npdb = 0;
	i__1 = sequen_1.nchain;
	for (m = 1; m <= i__1; ++m) {
	    start = ichain_ref(1, m);
	    stop = ichain_ref(2, m);
	    i__3 = stop;
	    for (i__ = start; i__ <= i__3; ++i__) {
		s_copy(resname, amino_ref(0, sequen_1.seqtyp[i__ - 1]), (
			ftnlen)3, (ftnlen)3);
		if (s_cmp(resname, "FOR", (ftnlen)3, (ftnlen)3) == 0) {
		    pdbatm_(" C  ", resname, &i__, &ci[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" O  ", resname, &i__, &oi[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		} else if (s_cmp(resname, "ACE", (ftnlen)3, (ftnlen)3) == 0) {
		    pdbatm_(" CH3", resname, &i__, &cai[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" C  ", resname, &i__, &ci[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" O  ", resname, &i__, &oi[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		} else if (s_cmp(resname, "NH2", (ftnlen)3, (ftnlen)3) == 0) {
		    pdbatm_(" N  ", resname, &i__, &ni[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		} else if (s_cmp(resname, "NME", (ftnlen)3, (ftnlen)3) == 0) {
		    pdbatm_(" N  ", resname, &i__, &ni[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" CH3", resname, &i__, &cai[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		} else {
		    pdbatm_(" N  ", resname, &i__, &ni[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" CA ", resname, &i__, &cai[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" C  ", resname, &i__, &ci[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" O  ", resname, &i__, &oi[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		}
		getside_(resname, &i__, &ci[i__ - 1], &cai[i__ - 1], &cbi, (
			ftnlen)3);
		if (s_cmp(resname, "CYS", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(
			resname, "CYX", (ftnlen)3, (ftnlen)3) == 0) {
		    s_copy(resname, "CYS", (ftnlen)3, (ftnlen)3);
		    i__4 = couple_1.n13[cbi - 1];
		    for (j = 1; j <= i__4; ++j) {
			if (atmtyp_1.atomic[i13_ref(j, cbi) - 1] == 16) {
			    s_copy(resname, "CYX", (ftnlen)3, (ftnlen)3);
			}
		    }
		}
		if (i__ == stop && ci[i__ - 1] != 0) {
		    i__4 = couple_1.n12[ci[i__ - 1] - 1];
		    for (j = 1; j <= i__4; ++j) {
			k = i12_ref(j, ci[i__ - 1]);
			if (atmtyp_1.atomic[k - 1] == 8 && k != oi[i__ - 1]) {
			    pdbatm_(" OXT", resname, &i__, &k, (ftnlen)4, (
				    ftnlen)3);
			    goto L40;
			}
		    }
L40:
		    ;
		}
		getproh_(resname, &i__, &m, &ni[i__ - 1], &cai[i__ - 1], &cbi,
			 (ftnlen)3);
	    }
	}
    }

/*     find the phosphates and other nucleotide backbone atoms */

    if (s_cmp(moltyp, "NUCACID", (ftnlen)7, (ftnlen)7) == 0) {
	attach_();
	m = 1;
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s_copy(resname, nuclz_ref(0, sequen_1.seqtyp[m - 1]), (ftnlen)3, (
		    ftnlen)3);
	    if (s_cmp(resname, "MP ", (ftnlen)3, (ftnlen)3) == 0) {
		if (atmtyp_1.atomic[i__ - 1] == 15) {
		    poi[m - 1] = i__;
		    ++m;
		}
	    }
	    if (atmtyp_1.atomic[i__ - 1] == 6 && couple_1.n12[i__ - 1] == 4) {
		cbone = FALSE_;
		nbone = FALSE_;
		obone = FALSE_;
		i__3 = couple_1.n12[i__ - 1];
		for (j = 1; j <= i__3; ++j) {
		    k = i12_ref(j, i__);
		    ka = atmtyp_1.atomic[k - 1];
		    kn = couple_1.n12[k - 1];
		    if (ka == 6) {
			cbone = TRUE_;
		    }
		    if (ka == 7 && kn == 3) {
			nbone = TRUE_;
		    }
		    if (ka == 8 && kn == 2) {
			obone = TRUE_;
		    }
		}
		if (cbone && nbone && obone) {
		    c1i[m - 1] = i__;
		    ++m;
		}
	    }
	}
	i__1 = sequen_1.nseq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    m = c1i[i__ - 1];
	    if (m != 0) {
		i__3 = couple_1.n12[m - 1];
		for (j = 1; j <= i__3; ++j) {
		    k = i12_ref(j, m);
		    ka = atmtyp_1.atomic[k - 1];
		    if (ka == 6) {
			c2i[i__ - 1] = k;
		    }
		    if (ka == 7) {
			ni[i__ - 1] = k;
		    }
		    if (ka == 8) {
			o4i[i__ - 1] = k;
		    }
		}
		m = o4i[i__ - 1];
		i__3 = couple_1.n12[m - 1];
		for (j = 1; j <= i__3; ++j) {
		    k = i12_ref(j, m);
		    ka = atmtyp_1.atomic[k - 1];
		    if (ka == 6 && k != c1i[i__ - 1]) {
			c4i[i__ - 1] = k;
		    }
		}
		m = c2i[i__ - 1];
		i__3 = couple_1.n12[m - 1];
		for (j = 1; j <= i__3; ++j) {
		    k = i12_ref(j, m);
		    ka = atmtyp_1.atomic[k - 1];
		    if (ka == 8) {
			o2i[i__ - 1] = k;
		    }
		    if (ka == 6 && k != c1i[i__ - 1]) {
			c3i[i__ - 1] = k;
		    }
		}
		m = c3i[i__ - 1];
		i__3 = couple_1.n12[m - 1];
		for (j = 1; j <= i__3; ++j) {
		    k = i12_ref(j, m);
		    ka = atmtyp_1.atomic[k - 1];
		    if (ka == 8) {
			o3i[i__ - 1] = k;
		    }
		}
		m = c4i[i__ - 1];
		i__3 = couple_1.n12[m - 1];
		for (j = 1; j <= i__3; ++j) {
		    k = i12_ref(j, m);
		    ka = atmtyp_1.atomic[k - 1];
		    if (ka == 6 && k != c3i[i__ - 1]) {
			c5i[i__ - 1] = k;
		    }
		}
		m = c5i[i__ - 1];
		i__3 = couple_1.n12[m - 1];
		for (j = 1; j <= i__3; ++j) {
		    k = i12_ref(j, m);
		    ka = atmtyp_1.atomic[k - 1];
		    if (ka == 8) {
			o5i[i__ - 1] = k;
		    }
		}
		m = o5i[i__ - 1];
		i__3 = couple_1.n12[m - 1];
		for (j = 1; j <= i__3; ++j) {
		    k = i12_ref(j, m);
		    ka = atmtyp_1.atomic[k - 1];
		    if (ka == 15) {
			poi[i__ - 1] = k;
		    }
		}
	    }
	    if (i__ > 1) {
		s_copy(resname, nuclz_ref(0, sequen_1.seqtyp[i__ - 2]), (
			ftnlen)3, (ftnlen)3);
		if (s_cmp(resname, "MP ", (ftnlen)3, (ftnlen)3) == 0) {
		    poi[i__ - 1] = 0;
		}
		if (s_cmp(resname, "DP ", (ftnlen)3, (ftnlen)3) == 0) {
		    poi[i__ - 1] = 0;
		}
		if (s_cmp(resname, "TP ", (ftnlen)3, (ftnlen)3) == 0) {
		    poi[i__ - 1] = 0;
		}
	    }
	    m = poi[i__ - 1];
	    if (m != 0) {
		i__3 = couple_1.n12[m - 1];
		for (j = 1; j <= i__3; ++j) {
		    k = i12_ref(j, m);
		    ka = atmtyp_1.atomic[k - 1];
		    if (ka == 8 && couple_1.n12[k - 1] == 1) {
			if (op1[i__ - 1] == 0) {
			    op1[i__ - 1] = k;
			} else if (op2[i__ - 1] == 0) {
			    op2[i__ - 1] = k;
			} else {
			    op3[i__ - 1] = k;
			}
		    }
		}
	    }
	}
    }

/*     get all the atoms for each nucleotide residue in order */

    if (s_cmp(moltyp, "NUCACID", (ftnlen)7, (ftnlen)7) == 0) {
	pdb_1.npdb = 0;
	i__1 = sequen_1.nchain;
	for (m = 1; m <= i__1; ++m) {
	    start = ichain_ref(1, m);
	    stop = ichain_ref(2, m);
	    i__3 = stop;
	    for (i__ = start; i__ <= i__3; ++i__) {
		s_copy(resname, nuclz_ref(0, sequen_1.seqtyp[i__ - 1]), (
			ftnlen)3, (ftnlen)3);
		if (s_cmp(resname, "MP ", (ftnlen)3, (ftnlen)3) == 0) {
		    pdbatm_(" P  ", resname, &i__, &poi[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" OP1", resname, &i__, &op1[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" OP2", resname, &i__, &op2[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" OP3", resname, &i__, &op3[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		} else if (s_cmp(resname, "DP ", (ftnlen)3, (ftnlen)3) == 0) {
		} else if (s_cmp(resname, "TP ", (ftnlen)3, (ftnlen)3) == 0) {
		} else {
		    pdbatm_(" P  ", resname, &i__, &poi[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" OP1", resname, &i__, &op1[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" OP2", resname, &i__, &op2[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" O5'", resname, &i__, &o5i[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" C5'", resname, &i__, &c5i[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" C4'", resname, &i__, &c4i[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" O4'", resname, &i__, &o4i[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" C3'", resname, &i__, &c3i[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" O3'", resname, &i__, &o3i[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" C2'", resname, &i__, &c2i[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" O2'", resname, &i__, &o2i[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    pdbatm_(" C1'", resname, &i__, &c1i[i__ - 1], (ftnlen)4, (
			    ftnlen)3);
		    getbase_(resname, &i__, &ni[i__ - 1], (ftnlen)3);
		    getnuch_(resname, &i__, &ni[i__ - 1], &c1i[i__ - 1], &c2i[
			    i__ - 1], &o2i[i__ - 1], &c3i[i__ - 1], &o3i[i__ 
			    - 1], &c4i[i__ - 1], &c5i[i__ - 1], &o5i[i__ - 1],
			     (ftnlen)3);
		}
	    }
	}
    }

/*     get any water, ions and ligands following biopolymer chains */

    if (s_cmp(moltyp, "GENERIC", (ftnlen)7, (ftnlen)7) != 0) {
	i__1 = molcul_1.nmol;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    hetmol[i__ - 1] = TRUE_;
	}
	i__1 = atoms_1.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (pdb_1.pdblist[i__ - 1] != 0) {
		hetmol[molcul_1.molcule[i__ - 1] - 1] = FALSE_;
	    }
	}
	i__1 = molcul_1.nmol;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (hetmol[i__ - 1]) {
		i__3 = imol_ref(2, i__);
		for (j = imol_ref(1, i__); j <= i__3; ++j) {
		    k = molcul_1.kmol[j - 1];
		    atmnum = atmtyp_1.atomic[k - 1];
/* Writing concatenation */
		    i__2[0] = 1, a__1[0] = " ";
		    i__2[1] = 3, a__1[1] = name___ref(0, k);
		    s_cat(atmname, a__1, i__2, &c__2, (ftnlen)4);
		    justify = 0;
		    numeral_(&atoms_1.type__[k - 1], resname, &justify, (
			    ftnlen)3);
		    if (water[i__ - 1]) {
			if (atmnum == 1) {
			    s_copy(atmname, " H  ", (ftnlen)4, (ftnlen)4);
			}
			if (atmnum == 6) {
			    s_copy(atmname, " O  ", (ftnlen)4, (ftnlen)4);
			}
			s_copy(resname, "HOH", (ftnlen)3, (ftnlen)3);
		    } else if (atmnum == 11) {
			s_copy(atmname, "NA  ", (ftnlen)4, (ftnlen)4);
			s_copy(resname, " NA", (ftnlen)3, (ftnlen)3);
		    } else if (atmnum == 12) {
			s_copy(atmname, "MG  ", (ftnlen)4, (ftnlen)4);
			s_copy(resname, " MG", (ftnlen)3, (ftnlen)3);
		    } else if (atmnum == 17) {
			s_copy(atmname, "CL  ", (ftnlen)4, (ftnlen)4);
			s_copy(resname, " CL", (ftnlen)3, (ftnlen)3);
		    } else if (atmnum == 19) {
			s_copy(atmname, " K  ", (ftnlen)4, (ftnlen)4);
			s_copy(resname, "  K", (ftnlen)3, (ftnlen)3);
		    } else if (atmnum == 20) {
			s_copy(atmname, "CA  ", (ftnlen)4, (ftnlen)4);
			s_copy(resname, " CA", (ftnlen)3, (ftnlen)3);
		    }
		    pdbnum = sequen_1.nseq + i__ - 1;
		    pdbatm_(atmname, resname, &pdbnum, &k, (ftnlen)4, (ftnlen)
			    3);
		    s_copy(pdbtyp_ref(0, pdb_1.npdb), "HETATM", (ftnlen)6, (
			    ftnlen)6);
		}
	    }
	}
	i__1 = molcul_1.nmol;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (hetmol[i__ - 1]) {
		i__3 = imol_ref(2, i__);
		for (j = imol_ref(1, i__); j <= i__3; ++j) {
		    k = molcul_1.kmol[j - 1];
		    kp = pdb_1.pdblist[k - 1];
		    pdb_1.npdb12[kp - 1] = couple_1.n12[k - 1];
		    i__4 = couple_1.n12[k - 1];
		    for (m = 1; m <= i__4; ++m) {
			ipdb12_ref(m, kp) = pdb_1.pdblist[i12_ref(m, k) - 1];
		    }
		}
	    }
	}
    }
    return 0;
} /* makepdb_ */

#undef pdbtyp_ref
#undef ichain_ref
#undef nuclz_ref
#undef amino_ref
#undef ipdb12_ref
#undef imol_ref
#undef name___ref
#undef seq_ref
#undef i14_ref
#undef i13_ref
#undef i12_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine pdbatm  --  add a single atom to PDB file  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "pdbatm" adds an atom to the Protein Data Bank file */


/* Subroutine */ int pdbatm_(char *atmname, char *resname, integer *ires, 
	integer *icoord, ftnlen atmname_len, ftnlen resname_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);


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




/*     for each atom set the sequential number, record type, atom */
/*     name, residue name, residue number and atomic coordinates */

    if (*icoord != 0) {
	++pdb_1.npdb;
	s_copy(pdbtyp_ref(0, pdb_1.npdb), "ATOM  ", (ftnlen)6, (ftnlen)6);
	s_copy(atmnam_ref(0, pdb_1.npdb), atmname, (ftnlen)4, (ftnlen)4);
	s_copy(resnam_ref(0, pdb_1.npdb), resname, (ftnlen)3, (ftnlen)3);
	pdb_1.resnum[pdb_1.npdb - 1] = *ires;
	pdb_1.xpdb[pdb_1.npdb - 1] = atoms_1.x[*icoord - 1];
	pdb_1.ypdb[pdb_1.npdb - 1] = atoms_1.y[*icoord - 1];
	pdb_1.zpdb[pdb_1.npdb - 1] = atoms_1.z__[*icoord - 1];
	pdb_1.npdb12[pdb_1.npdb - 1] = 0;
	pdb_1.pdblist[*icoord - 1] = pdb_1.npdb;
    }
    return 0;
} /* pdbatm_ */

#undef pdbtyp_ref
#undef resnam_ref
#undef atmnam_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine getside  --  extract the amino acid side chains  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "getside" finds the side chain heavy atoms for a single amino */
/*     acid residue and copies the names and coordinates to the Protein */
/*     Data Bank file */


/* Subroutine */ int getside_(char *resname, integer *ires, integer *ci, 
	integer *cai, integer *cbi, ftnlen resname_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int pdbatm_(char *, char *, integer *, integer *, 
	    ftnlen, ftnlen);


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




/*     if residue is glycine or a cap, there is no side chain */

    *cbi = 0;
    if (s_cmp(resname, "GLY", (ftnlen)3, (ftnlen)3) == 0) {
	return 0;
    }
    if (s_cmp(resname, "UNK", (ftnlen)3, (ftnlen)3) == 0) {
	return 0;
    }
    if (s_cmp(resname, "FOR", (ftnlen)3, (ftnlen)3) == 0) {
	return 0;
    }
    if (s_cmp(resname, "ACE", (ftnlen)3, (ftnlen)3) == 0) {
	return 0;
    }
    if (s_cmp(resname, "NH2", (ftnlen)3, (ftnlen)3) == 0) {
	return 0;
    }
    if (s_cmp(resname, "NME", (ftnlen)3, (ftnlen)3) == 0) {
	return 0;
    }

/*     find the beta carbon atom for the current residue */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ != *ci && atmtyp_1.atomic[i__ - 1] == 6) {
	    for (j = 1; j <= 4; ++j) {
		if (i12_ref(j, i__) == *cai) {
		    *cbi = i__;
		    if (s_cmp(resname, "AIB", (ftnlen)3, (ftnlen)3) != 0) {
			pdbatm_(" CB ", resname, ires, cbi, (ftnlen)4, (
				ftnlen)3);
		    } else {
			pdbatm_(" CB1", resname, ires, cbi, (ftnlen)4, (
				ftnlen)3);
		    }
		    goto L10;
		}
	    }
	}
    }
L10:

/*     glycine residue  (GLY) */

    if (s_cmp(resname, "GLY", (ftnlen)3, (ftnlen)3) == 0) {

/*     alanine residue  (ALA) */

    } else if (s_cmp(resname, "ALA", (ftnlen)3, (ftnlen)3) == 0) {

/*     valine residue  (VAL) */

    } else if (s_cmp(resname, "VAL", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CG2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     leucine residue  (LEU) */

    } else if (s_cmp(resname, "LEU", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" CD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     isoleucine residue  (ILE) */

    } else if (s_cmp(resname, "ILE", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CG2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" CD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     serine residue  (SER) */

    } else if (s_cmp(resname, "SER", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" OG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     threonine residue  (THR) */

    } else if (s_cmp(resname, "THR", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" OG1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CG2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     cysteine residue  (CYS) */

    } else if (s_cmp(resname, "CYS", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" SG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     cysteine residue  (CYX) */

    } else if (s_cmp(resname, "CYX", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" SG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     proline residue  (PRO) */

    } else if (s_cmp(resname, "PRO", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CD ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     phenylalanine residue  (PHE) */

    } else if (s_cmp(resname, "PHE", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" CD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 4;
	pdbatm_(" CE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 5;
	pdbatm_(" CE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 6;
	pdbatm_(" CZ ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     tyrosine residue  (TYR) */

    } else if (s_cmp(resname, "TYR", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" CD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 4;
	pdbatm_(" CE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 5;
	pdbatm_(" CE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 6;
	pdbatm_(" CZ ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 7;
	pdbatm_(" OH ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     tryptophan residue  (TRP) */

    } else if (s_cmp(resname, "TRP", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" CD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 4;
	pdbatm_(" NE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 5;
	pdbatm_(" CE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 6;
	pdbatm_(" CE3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 7;
	pdbatm_(" CZ2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 8;
	pdbatm_(" CZ3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 9;
	pdbatm_(" CH2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     histidine (HD and HE) residue  (HIS) */

    } else if (s_cmp(resname, "HIS", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" ND1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" CD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 4;
	pdbatm_(" CE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 5;
	pdbatm_(" NE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     histidine (HD only) residue  (HID) */

    } else if (s_cmp(resname, "HID", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" ND1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" CD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 4;
	pdbatm_(" CE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 5;
	pdbatm_(" NE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     histidine (HE only) residue  (HIE) */

    } else if (s_cmp(resname, "HIE", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" ND1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" CD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 4;
	pdbatm_(" CE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 5;
	pdbatm_(" NE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     aspartic acid residue  (ASP) */

    } else if (s_cmp(resname, "ASP", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" OD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" OD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     asparagine residue  (ASN) */

    } else if (s_cmp(resname, "ASN", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" OD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" ND2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     glutamic acid residue  (GLU) */

    } else if (s_cmp(resname, "GLU", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CD ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" OE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 4;
	pdbatm_(" OE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     glutamine residue  (GLN) */

    } else if (s_cmp(resname, "GLN", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CD ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" OE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 4;
	pdbatm_(" NE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     methionine residue  (MET) */

    } else if (s_cmp(resname, "MET", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" SD ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" CE ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     lysine residue  (LYS) */

    } else if (s_cmp(resname, "LYS", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CD ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" CE ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 4;
	pdbatm_(" NZ ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     arginine residue  (ARG) */

    } else if (s_cmp(resname, "ARG", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CD ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" NE ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 4;
	pdbatm_(" CZ ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 5;
	pdbatm_(" NH1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 6;
	pdbatm_(" NH2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     ornithine residue  (ORN) */

    } else if (s_cmp(resname, "ORN", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CD ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" NE ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     methylalanine residue  (AIB) */

    } else if (s_cmp(resname, "AIB", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     pyroglutamic acid residue  (PCA) */

    } else if (s_cmp(resname, "PCA", (ftnlen)3, (ftnlen)3) == 0) {
	i__1 = *cbi + 1;
	pdbatm_(" CG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 2;
	pdbatm_(" CD ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *cbi + 3;
	pdbatm_(" OE ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     unknown residue  (UNK) */

    } else if (s_cmp(resname, "UNK", (ftnlen)3, (ftnlen)3) == 0) {
    }
    return 0;
} /* getside_ */

#undef i12_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine getproh  --  extract the amino acid hydrogens  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "getproh" finds the hydrogen atoms for a single amino acid */
/*     residue and copies the names and coordinates to the Protein */
/*     Data Bank file */


/* Subroutine */ int getproh_(char *resname, integer *ires, integer *jchain, 
	integer *ni, integer *cai, integer *cbi, ftnlen resname_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, nh, hca;
    extern /* Subroutine */ int pdbatm_(char *, char *, integer *, integer *, 
	    ftnlen, ftnlen);
    static char atmname[4];
    static logical allatom;


#define i12_ref(a_1,a_2) couple_1.i12[(a_2)*8 + a_1 - 9]
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




/*     get any amide hydrogen atoms for non-N-terminal residues */

    if (*ires != ichain_ref(1, *jchain) || couple_1.n12[*ni - 1] != 4) {
	if (s_cmp(resname, "PRO", (ftnlen)3, (ftnlen)3) != 0) {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (atmtyp_1.atomic[i__ - 1] == 1 && i12_ref(1, i__) == *ni) {
		    if (s_cmp(resname, "NH2", (ftnlen)3, (ftnlen)3) == 0) {
			pdbatm_(" H1 ", resname, ires, &i__, (ftnlen)4, (
				ftnlen)3);
			i__2 = i__ + 1;
			pdbatm_(" H2 ", resname, ires, &i__2, (ftnlen)4, (
				ftnlen)3);
		    } else {
			pdbatm_(" H  ", resname, ires, &i__, (ftnlen)4, (
				ftnlen)3);
		    }
		    goto L10;
		}
	    }
	}

/*     get any amide hydrogen atoms for N-terminal residues */

    } else {
	if (s_cmp(resname, "PRO", (ftnlen)3, (ftnlen)3) == 0) {
	    nh = 0;
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (atmtyp_1.atomic[i__ - 1] == 1 && i12_ref(1, i__) == *ni) {
		    ++nh;
		    if (nh == 1) {
			s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
		    } else if (nh == 2) {
			s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
		    }
		    pdbatm_(atmname, resname, ires, &i__, (ftnlen)4, (ftnlen)
			    3);
		    if (nh == 2) {
			goto L10;
		    }
		}
	    }
	} else if (s_cmp(resname, "PCA", (ftnlen)3, (ftnlen)3) == 0) {
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (atmtyp_1.atomic[i__ - 1] == 1 && i12_ref(1, i__) == *ni) {
		    s_copy(atmname, " H  ", (ftnlen)4, (ftnlen)4);
		    pdbatm_(atmname, resname, ires, &i__, (ftnlen)4, (ftnlen)
			    3);
		    goto L10;
		}
	    }
	} else {
	    nh = 0;
	    i__1 = atoms_1.n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (atmtyp_1.atomic[i__ - 1] == 1 && i12_ref(1, i__) == *ni) {
		    ++nh;
		    if (nh == 1) {
			s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
		    } else if (nh == 2) {
			s_copy(atmname, " H2 ", (ftnlen)4, (ftnlen)4);
		    } else if (nh == 3) {
			s_copy(atmname, " H3 ", (ftnlen)4, (ftnlen)4);
		    }
		    pdbatm_(atmname, resname, ires, &i__, (ftnlen)4, (ftnlen)
			    3);
		    if (nh == 3) {
			goto L10;
		    }
		}
	    }
	}
    }
L10:

/*     get the alpha hydrogen atom for the current residue */

    hca = 0;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (atmtyp_1.atomic[i__ - 1] == 1 && i12_ref(1, i__) == *cai) {
	    hca = i__;
	    if (s_cmp(resname, "GLY", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(atmname, " HA2", (ftnlen)4, (ftnlen)4);
	    } else if (s_cmp(resname, "FOR", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(atmname, " H  ", (ftnlen)4, (ftnlen)4);
	    } else if (s_cmp(resname, "ACE", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	    } else if (s_cmp(resname, "NME", (ftnlen)3, (ftnlen)3) == 0) {
		s_copy(atmname, " H1 ", (ftnlen)4, (ftnlen)4);
	    } else {
		s_copy(atmname, " HA ", (ftnlen)4, (ftnlen)4);
	    }
	    pdbatm_(atmname, resname, ires, &i__, (ftnlen)4, (ftnlen)3);
	    goto L20;
	}
    }
L20:

/*     if no alpha hydrogen, then united atom force field */

    if (hca != 0) {
	allatom = TRUE_;
    } else if (s_cmp(resname, "AIB", (ftnlen)3, (ftnlen)3) == 0) {
	if (couple_1.n12[*cbi - 1] == 1) {
	    allatom = FALSE_;
	} else {
	    allatom = TRUE_;
	}
    } else {
	allatom = FALSE_;
    }

/*     glycine residue  (GLY) */

    if (s_cmp(resname, "GLY", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 1;
	    pdbatm_(" HA3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     alanine residue  (ALA) */

    } else if (s_cmp(resname, "ALA", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 2;
	    pdbatm_(" HB1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 3;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 4;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     valine residue  (VAL) */

    } else if (s_cmp(resname, "VAL", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 4;
	    pdbatm_(" HB ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 5;
	    pdbatm_("HG11", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 6;
	    pdbatm_("HG12", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 7;
	    pdbatm_("HG13", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_("HG21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_("HG22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 10;
	    pdbatm_("HG23", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     leucine residue  (LEU) */

    } else if (s_cmp(resname, "LEU", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 5;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 6;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 7;
	    pdbatm_(" HG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_("HD11", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_("HD12", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 10;
	    pdbatm_("HD13", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 11;
	    pdbatm_("HD21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 12;
	    pdbatm_("HD22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 13;
	    pdbatm_("HD23", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     isoleucine residue  (ILE) */

    } else if (s_cmp(resname, "ILE", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 5;
	    pdbatm_(" HB ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 6;
	    pdbatm_("HG12", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 7;
	    pdbatm_("HG13", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_("HG21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_("HG22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 10;
	    pdbatm_("HG23", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 11;
	    pdbatm_("HD11", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 12;
	    pdbatm_("HD12", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 13;
	    pdbatm_("HD13", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     serine residue  (SER) */

    } else if (s_cmp(resname, "SER", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 3;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 4;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 5;
	    pdbatm_(" HG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *cbi + 2;
	    pdbatm_(" HG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     threonine residue  (THR) */

    } else if (s_cmp(resname, "THR", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 4;
	    pdbatm_(" HB ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 5;
	    pdbatm_(" HG1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 6;
	    pdbatm_("HG21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 7;
	    pdbatm_("HG22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_("HG23", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *cbi + 3;
	    pdbatm_(" HG1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     cysteine residue  (CYS) */

    } else if (s_cmp(resname, "CYS", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 3;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 4;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 5;
	    pdbatm_(" HG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *cbi + 2;
	    pdbatm_(" HG ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     cystine residue  (CYX) */

    } else if (s_cmp(resname, "CYX", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 3;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 4;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     proline residue  (PRO) */

    } else if (s_cmp(resname, "PRO", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 4;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 5;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 6;
	    pdbatm_(" HG2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 7;
	    pdbatm_(" HG3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_(" HD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_(" HD3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     phenylalanine residue  (PHE) */

    } else if (s_cmp(resname, "PHE", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 8;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 10;
	    pdbatm_(" HD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 11;
	    pdbatm_(" HD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 12;
	    pdbatm_(" HE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 13;
	    pdbatm_(" HE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 14;
	    pdbatm_(" HZ ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     tyrosine residue  (TYR) */

    } else if (s_cmp(resname, "TYR", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 9;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 10;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 11;
	    pdbatm_(" HD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 12;
	    pdbatm_(" HD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 13;
	    pdbatm_(" HE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 14;
	    pdbatm_(" HE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 15;
	    pdbatm_(" HH ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *cbi + 12;
	    pdbatm_(" HH ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     tryptophan residue  (TRP) */

    } else if (s_cmp(resname, "TRP", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 11;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 12;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 13;
	    pdbatm_(" HD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 14;
	    pdbatm_(" HE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 15;
	    pdbatm_(" HE3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 16;
	    pdbatm_(" HZ2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 17;
	    pdbatm_(" HZ3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 18;
	    pdbatm_(" HH2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *cbi + 11;
	    pdbatm_(" HE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     histidine (HD and HE) residue  (HIS) */

    } else if (s_cmp(resname, "HIS", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 7;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_(" HD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 10;
	    pdbatm_(" HD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 11;
	    pdbatm_(" HE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 12;
	    pdbatm_(" HE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *cbi + 6;
	    pdbatm_(" HD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 9;
	    pdbatm_(" HE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     histidine (HD only) residue  (HID) */

    } else if (s_cmp(resname, "HID", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 7;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_(" HD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 10;
	    pdbatm_(" HD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 11;
	    pdbatm_(" HE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *cbi + 6;
	    pdbatm_(" HD1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     histidine (HE only) residue  (HIE) */

    } else if (s_cmp(resname, "HIE", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 7;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_(" HD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 10;
	    pdbatm_(" HE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 11;
	    pdbatm_(" HE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *cbi + 8;
	    pdbatm_(" HE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     aspartic acid residue  (ASP) */

    } else if (s_cmp(resname, "ASP", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 5;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 6;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     asparagine residue  (ASN) */

    } else if (s_cmp(resname, "ASN", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 5;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 6;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 7;
	    pdbatm_("HD21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_("HD22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *cbi + 4;
	    pdbatm_("HD21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 5;
	    pdbatm_("HD22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     glutamic acid residue  (GLU) */

    } else if (s_cmp(resname, "GLU", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 6;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 7;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_(" HG2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_(" HG3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     glutamine residue  (GLN) */

    } else if (s_cmp(resname, "GLN", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 6;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 7;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_(" HG2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_(" HG3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 10;
	    pdbatm_("HE21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 11;
	    pdbatm_("HE22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *cbi + 5;
	    pdbatm_("HE21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 6;
	    pdbatm_("HE22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     methionine residue  (MET) */

    } else if (s_cmp(resname, "MET", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 5;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 6;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 7;
	    pdbatm_(" HG2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_(" HG3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_(" HE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 10;
	    pdbatm_(" HE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 11;
	    pdbatm_(" HE3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     lysine residue  (LYS) */

    } else if (s_cmp(resname, "LYS", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 6;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 7;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_(" HG2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_(" HG3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 10;
	    pdbatm_(" HD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 11;
	    pdbatm_(" HD3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 12;
	    pdbatm_(" HE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 13;
	    pdbatm_(" HE3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 14;
	    pdbatm_(" HZ1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 15;
	    pdbatm_(" HZ2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 16;
	    pdbatm_(" HZ3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *cbi + 5;
	    pdbatm_(" HZ1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 6;
	    pdbatm_(" HZ2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 7;
	    pdbatm_(" HZ3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     arginine residue  (ARG) */

    } else if (s_cmp(resname, "ARG", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 8;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 10;
	    pdbatm_(" HG2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 11;
	    pdbatm_(" HG3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 12;
	    pdbatm_(" HD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 13;
	    pdbatm_(" HD3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 14;
	    pdbatm_(" HE ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 15;
	    pdbatm_("HH11", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 16;
	    pdbatm_("HH12", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 17;
	    pdbatm_("HH21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 18;
	    pdbatm_("HH22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *cbi + 7;
	    pdbatm_(" HE ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 8;
	    pdbatm_("HH11", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 9;
	    pdbatm_("HH12", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 10;
	    pdbatm_("HH21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 11;
	    pdbatm_("HH22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     ornithine residue  (ORN) */

    } else if (s_cmp(resname, "ORN", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 5;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 6;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 7;
	    pdbatm_(" HG2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_(" HG3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 9;
	    pdbatm_(" HD2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 10;
	    pdbatm_(" HD3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 11;
	    pdbatm_(" HE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 12;
	    pdbatm_(" HE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 13;
	    pdbatm_(" HE3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *cbi + 4;
	    pdbatm_(" HE1", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 5;
	    pdbatm_(" HE2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 6;
	    pdbatm_(" HE3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     methylalanine residue  (AIB) */

    } else if (s_cmp(resname, "AIB", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = *cbi + 2;
	    pdbatm_("HB11", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 3;
	    pdbatm_("HB12", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 4;
	    pdbatm_("HB13", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 5;
	    pdbatm_("HB21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 6;
	    pdbatm_("HB22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *cbi + 7;
	    pdbatm_("HB23", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     pyroglutamic acid residue  (PCA) */

    } else if (s_cmp(resname, "PCA", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 5;
	    pdbatm_(" HB2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 6;
	    pdbatm_(" HB3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 7;
	    pdbatm_(" HG2", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 8;
	    pdbatm_(" HG3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     unknown residue  (UNK) */

    } else if (s_cmp(resname, "UNK", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 1;
	    pdbatm_(" HA3", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     N-terminal acetyl residue  (ACE) */

    } else if (s_cmp(resname, "ACE", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 1;
	    pdbatm_(" H2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 2;
	    pdbatm_(" H3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     N-terminal formyl residue  (FOR) */

    } else if (s_cmp(resname, "FOR", (ftnlen)3, (ftnlen)3) == 0) {

/*     C-terminal N-methylamide residue  (NME) */

    } else if (s_cmp(resname, "NME", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = hca + 1;
	    pdbatm_(" H2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = hca + 2;
	    pdbatm_(" H3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     C-terminal amide residue  (NH2) */

    } else if (s_cmp(resname, "NH2", (ftnlen)3, (ftnlen)3) == 0) {
    }
    return 0;
} /* getproh_ */

#undef ichain_ref
#undef i12_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine getbase  --  extract the nucleotide side chains  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "getbase" finds the base heavy atoms for a single nucleotide */
/*     residue and copies the names and coordinates to the Protein */
/*     Data Bank file */


/* Subroutine */ int getbase_(char *resname, integer *ires, integer *ni, 
	ftnlen resname_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int pdbatm_(char *, char *, integer *, integer *, 
	    ftnlen, ftnlen);



/*     adenine in adenosine residue  (A) */

    if (s_cmp(resname, "A  ", (ftnlen)3, (ftnlen)3) == 0) {
	pdbatm_(" N9 ", resname, ires, ni, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 1;
	pdbatm_(" C8 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 2;
	pdbatm_(" N7 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 3;
	pdbatm_(" C5 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 4;
	pdbatm_(" C6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 5;
	pdbatm_(" N6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 6;
	pdbatm_(" N1 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 7;
	pdbatm_(" C2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 8;
	pdbatm_(" N3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 9;
	pdbatm_(" C4 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     guanine in guanosine residue  (G) */

    } else if (s_cmp(resname, "G  ", (ftnlen)3, (ftnlen)3) == 0) {
	pdbatm_(" N9 ", resname, ires, ni, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 1;
	pdbatm_(" C8 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 2;
	pdbatm_(" N7 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 3;
	pdbatm_(" C5 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 4;
	pdbatm_(" C6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 5;
	pdbatm_(" O6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 6;
	pdbatm_(" N1 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 7;
	pdbatm_(" C2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 8;
	pdbatm_(" N2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 9;
	pdbatm_(" N3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 10;
	pdbatm_(" C4 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     cytosine in cytidine residue  (C) */

    } else if (s_cmp(resname, "C  ", (ftnlen)3, (ftnlen)3) == 0) {
	pdbatm_(" N1 ", resname, ires, ni, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 1;
	pdbatm_(" C2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 2;
	pdbatm_(" O2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 3;
	pdbatm_(" N3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 4;
	pdbatm_(" C4 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 5;
	pdbatm_(" N4 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 6;
	pdbatm_(" C5 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 7;
	pdbatm_(" C6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     uracil in uridine residue  (U) */

    } else if (s_cmp(resname, "U  ", (ftnlen)3, (ftnlen)3) == 0) {
	pdbatm_(" N1 ", resname, ires, ni, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 1;
	pdbatm_(" C2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 2;
	pdbatm_(" O2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 3;
	pdbatm_(" N3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 4;
	pdbatm_(" C4 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 5;
	pdbatm_(" O4 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 6;
	pdbatm_(" C5 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 7;
	pdbatm_(" C6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     adenine in deoxyadenosine residue  (DA) */

    } else if (s_cmp(resname, "DA ", (ftnlen)3, (ftnlen)3) == 0) {
	pdbatm_(" N9 ", resname, ires, ni, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 1;
	pdbatm_(" C8 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 2;
	pdbatm_(" N7 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 3;
	pdbatm_(" C5 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 4;
	pdbatm_(" C6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 5;
	pdbatm_(" N6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 6;
	pdbatm_(" N1 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 7;
	pdbatm_(" C2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 8;
	pdbatm_(" N3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 9;
	pdbatm_(" C4 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     guanine in deoxyguanosine residue  (DG) */

    } else if (s_cmp(resname, "DG ", (ftnlen)3, (ftnlen)3) == 0) {
	pdbatm_(" N9 ", resname, ires, ni, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 1;
	pdbatm_(" C8 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 2;
	pdbatm_(" N7 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 3;
	pdbatm_(" C5 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 4;
	pdbatm_(" C6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 5;
	pdbatm_(" O6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 6;
	pdbatm_(" N1 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 7;
	pdbatm_(" C2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 8;
	pdbatm_(" N2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 9;
	pdbatm_(" N3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 10;
	pdbatm_(" C4 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     cytosine in deoxycytidine residue  (DC) */

    } else if (s_cmp(resname, "DC ", (ftnlen)3, (ftnlen)3) == 0) {
	pdbatm_(" N1 ", resname, ires, ni, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 1;
	pdbatm_(" C2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 2;
	pdbatm_(" O2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 3;
	pdbatm_(" N3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 4;
	pdbatm_(" C4 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 5;
	pdbatm_(" N4 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 6;
	pdbatm_(" C5 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 7;
	pdbatm_(" C6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);

/*     thymine in deoxythymidine residue  (DT) */

    } else if (s_cmp(resname, "DT ", (ftnlen)3, (ftnlen)3) == 0) {
	pdbatm_(" N1 ", resname, ires, ni, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 1;
	pdbatm_(" C2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 2;
	pdbatm_(" O2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 3;
	pdbatm_(" N3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 4;
	pdbatm_(" C4 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 5;
	pdbatm_(" O4 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 6;
	pdbatm_(" C5 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 7;
	pdbatm_(" C7 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	i__1 = *ni + 8;
	pdbatm_(" C6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
    }
    return 0;
} /* getbase_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine getnuch  --  extract the nucleotide hydrogens  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "getnuch" finds the nucleotide hydrogen atoms for a single */
/*     residue and copies the names and coordinates to the Protein */
/*     Data Bank file */


/* Subroutine */ int getnuch_(char *resname, integer *ires, integer *ni, 
	integer *c1i, integer *c2i, integer *o2i, integer *c3i, integer *o3i, 
	integer *c4i, integer *c5i, integer *o5i, ftnlen resname_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, k;
    static logical done;
    extern /* Subroutine */ int pdbatm_(char *, char *, integer *, integer *, 
	    ftnlen, ftnlen);
    static logical allatom;


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




/*     if no ribose C1 hydrogen, then united atom force field */

    allatom = TRUE_;
    if (couple_1.n12[*c1i - 1] != 4) {
	allatom = FALSE_;
    }

/*     get sugar ring hydrogen atoms for the current residue */

    done = FALSE_;
    i__1 = couple_1.n12[*c5i - 1];
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = i12_ref(i__, *c5i);
	if (atmtyp_1.atomic[k - 1] == 1) {
	    if (! done) {
		pdbatm_(" H5'", resname, ires, &k, (ftnlen)4, (ftnlen)3);
		done = TRUE_;
	    } else {
		pdbatm_("H5''", resname, ires, &k, (ftnlen)4, (ftnlen)3);
	    }
	}
    }
    i__1 = couple_1.n12[*c4i - 1];
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = i12_ref(i__, *c4i);
	if (atmtyp_1.atomic[k - 1] == 1) {
	    pdbatm_(" H4'", resname, ires, &k, (ftnlen)4, (ftnlen)3);
	}
    }
    i__1 = couple_1.n12[*c3i - 1];
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = i12_ref(i__, *c3i);
	if (atmtyp_1.atomic[k - 1] == 1) {
	    pdbatm_(" H3'", resname, ires, &k, (ftnlen)4, (ftnlen)3);
	}
    }
    done = FALSE_;
    i__1 = couple_1.n12[*c2i - 1];
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = i12_ref(i__, *c2i);
	if (atmtyp_1.atomic[k - 1] == 1) {
	    if (! done) {
		pdbatm_(" H2'", resname, ires, &k, (ftnlen)4, (ftnlen)3);
		done = TRUE_;
	    } else {
		pdbatm_("H2''", resname, ires, &k, (ftnlen)4, (ftnlen)3);
	    }
	}
    }
    if (*o2i != 0) {
	i__1 = couple_1.n12[*o2i - 1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = i12_ref(i__, *o2i);
	    if (atmtyp_1.atomic[k - 1] == 1) {
		pdbatm_("HO2'", resname, ires, &k, (ftnlen)4, (ftnlen)3);
	    }
	}
    }
    i__1 = couple_1.n12[*c1i - 1];
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = i12_ref(i__, *c1i);
	if (atmtyp_1.atomic[k - 1] == 1) {
	    pdbatm_(" H1'", resname, ires, &k, (ftnlen)4, (ftnlen)3);
	}
    }

/*     adenine in adenosine residue  (A) */

    if (s_cmp(resname, "A  ", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = *ni + 10;
	    pdbatm_(" H8 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 11;
	    pdbatm_(" H61", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 12;
	    pdbatm_(" H62", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 13;
	    pdbatm_(" H2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *ni + 10;
	    pdbatm_(" H61", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 11;
	    pdbatm_(" H62", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     guanine in guanosine residue  (G) */

    } else if (s_cmp(resname, "G  ", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = *ni + 11;
	    pdbatm_(" H8 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 12;
	    pdbatm_(" H1 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 13;
	    pdbatm_(" H21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 14;
	    pdbatm_(" H22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *ni + 11;
	    pdbatm_(" H1 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 12;
	    pdbatm_(" H21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 13;
	    pdbatm_(" H22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     cytosine in cytidine residue  (C) */

    } else if (s_cmp(resname, "C  ", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = *ni + 8;
	    pdbatm_(" H41", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 9;
	    pdbatm_(" H42", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 10;
	    pdbatm_(" H5 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 11;
	    pdbatm_(" H6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *ni + 8;
	    pdbatm_(" H41", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 9;
	    pdbatm_(" H42", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     uracil in uridine residue  (U) */

    } else if (s_cmp(resname, "U  ", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = *ni + 8;
	    pdbatm_(" H3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 9;
	    pdbatm_(" H5 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 10;
	    pdbatm_(" H6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *ni + 8;
	    pdbatm_(" H3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     adenine in deoxyadenosine residue  (DA) */

    } else if (s_cmp(resname, "DA ", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = *ni + 10;
	    pdbatm_(" H8 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 11;
	    pdbatm_(" H61", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 12;
	    pdbatm_(" H62", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 13;
	    pdbatm_(" H2 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *ni + 10;
	    pdbatm_(" H61", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 11;
	    pdbatm_(" H62", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     guanine in deoxyguanosine residue  (DG) */

    } else if (s_cmp(resname, "DG ", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = *ni + 11;
	    pdbatm_(" H8 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 12;
	    pdbatm_(" H1 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 13;
	    pdbatm_(" H21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 14;
	    pdbatm_(" H22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *ni + 11;
	    pdbatm_(" H1 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 12;
	    pdbatm_(" H21", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 13;
	    pdbatm_(" H22", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     cytosine in deoxycytidine residue  (DC) */

    } else if (s_cmp(resname, "DC ", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = *ni + 8;
	    pdbatm_(" H41", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 9;
	    pdbatm_(" H42", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 10;
	    pdbatm_(" H5 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 11;
	    pdbatm_(" H6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *ni + 8;
	    pdbatm_(" H41", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 9;
	    pdbatm_(" H42", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}

/*     thymine in deoxythymidine residue  (DT) */

    } else if (s_cmp(resname, "DT ", (ftnlen)3, (ftnlen)3) == 0) {
	if (allatom) {
	    i__1 = *ni + 9;
	    pdbatm_(" H3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 10;
	    pdbatm_(" H71", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 11;
	    pdbatm_(" H72", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 12;
	    pdbatm_(" H73", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	    i__1 = *ni + 13;
	    pdbatm_(" H6 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	} else {
	    i__1 = *ni + 9;
	    pdbatm_(" H3 ", resname, ires, &i__1, (ftnlen)4, (ftnlen)3);
	}
    }

/*     get any capping hydrogen atoms for the current residue */

    i__1 = couple_1.n12[*o5i - 1];
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = i12_ref(i__, *o5i);
	if (atmtyp_1.atomic[k - 1] == 1) {
	    pdbatm_(" H5T", resname, ires, &k, (ftnlen)4, (ftnlen)3);
	}
    }
    i__1 = couple_1.n12[*o3i - 1];
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = i12_ref(i__, *o3i);
	if (atmtyp_1.atomic[k - 1] == 1) {
	    pdbatm_(" H3T", resname, ires, &k, (ftnlen)4, (ftnlen)3);
	}
    }
    return 0;
} /* getnuch_ */

#undef i12_ref


