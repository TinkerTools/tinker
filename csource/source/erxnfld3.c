/* erxnfld3.f -- translated by f2c (version 20050501).
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
    integer neb, nea, neba, neub, neaa, neopb, neopd, neid, neit, net, nept, 
	    nebt, nett, nev, nec, necd, ned, nem, nep, new__, ner, nes, nelf, 
	    neg, nex;
} action_;

#define action_1 action_

struct {
    doublereal aesum[25000], aeb[25000], aea[25000], aeba[25000], aeub[25000],
	     aeaa[25000], aeopb[25000], aeopd[25000], aeid[25000], aeit[25000]
	    , aet[25000], aept[25000], aebt[25000], aett[25000], aev[25000], 
	    aec[25000], aecd[25000], aed[25000], aem[25000], aep[25000], aer[
	    25000], aes[25000], aelf[25000], aeg[25000], aex[25000];
} analyz_;

#define analyz_1 analyz_

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
    doublereal electric, dielec, ebuffer, c2scale, c3scale, c4scale, c5scale;
    logical neutnbr, neutcut;
} chgpot_;

#define chgpot_1 chgpot_

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
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_

struct {
    doublereal off, off2, cut, cut2, c0, c1, c2, c3, c4, c5, f0, f1, f2, f3, 
	    f4, f5, f6, f7;
} shunt_;

#define shunt_1 shunt_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

/* Table of constant values */

static integer c__1 = 1;



/*     ############################################################ */
/*     ##  COPYRIGHT (C) 1996 by Yong Kong & Jay William Ponder  ## */
/*     ##                  All Rights Reserved                   ## */
/*     ############################################################ */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine erxnfld3  --  reaction field energy & analysis  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "erxnfld3" calculates the macroscopic reaction field energy, */
/*     and also partitions the energy among the atoms */

/*     literature reference: */

/*     Y. Kong and J. W. Ponder, "Reaction Field Methods for Off-Center */
/*     Multipoles", Journal of Chemical Physics, 107, 481-492 (1997) */


/* Subroutine */ int erxnfld3_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Individual Reaction Field\002,\002 Inter"
	    "actions :\002,//,\002 Type\002,13x,\002Atom Names\002,9x,\002Dis"
	    "t from Origin\002,4x,\002R(1-2)\002,6x,\002Energy\002,/)";
    static char fmt_20[] = "(\002 RxnFld\002,5x,i5,\002-\002,a3,1x,i5,\002"
	    "-\002,a3,2x,3f10.4,f12.4)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    double sqrt(doublereal);
    integer do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k;
    static doublereal r__, r2, di, dk;
    static integer ii, kk, ix, iz, kx, kz;
    static doublereal xr, yr, zr, eik, rpi[13], rpk[13];
    static logical huge__, usei, usek;
    extern /* Subroutine */ int erfik_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *);
    static logical header;
    extern /* Subroutine */ int switch_(char *, ftnlen), ijkpts_(void), 
	    chkpole_(void), rotpole_(void);

    /* Fortran I/O blocks */
    static cilist io___21 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_20, 0 };



#define name___ref(a_0,a_1) &atmtyp_1.name__[(a_1)*3 + a_0 - 3]
#define rpole_ref(a_1,a_2) mpole_1.rpole[(a_2)*13 + a_1 - 14]



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

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  shunt.i  --  polynomial switching function coefficients  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     off    distance at which the potential energy goes to zero */
/*     off2   square of distance at which the potential goes to zero */
/*     cut    distance at which switching of the potential begins */
/*     cut2   square of distance at which the switching begins */
/*     c0     zeroth order coefficient of multiplicative switch */
/*     c1     first order coefficient of multiplicative switch */
/*     c2     second order coefficient of multiplicative switch */
/*     c3     third order coefficient of multiplicative switch */
/*     c4     fourth order coefficient of multiplicative switch */
/*     c5     fifth order coefficient of multiplicative switch */
/*     f0     zeroth order coefficient of additive switch function */
/*     f1     first order coefficient of additive switch function */
/*     f2     second order coefficient of additive switch function */
/*     f3     third order coefficient of additive switch function */
/*     f4     fourth order coefficient of additive switch function */
/*     f5     fifth order coefficient of additive switch function */
/*     f6     sixth order coefficient of additive switch function */
/*     f7     seventh order coefficient of additive switch function */




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




/*     zero out the reaction field energy and partitioning */

    action_1.ner = 0;
    energi_1.er = 0.;
    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	analyz_1.aer[i__ - 1] = 0.;
    }
    header = TRUE_;

/*     set the switching function coefficients */

    switch_("MPOLE", (ftnlen)5);

/*     check the sign of multipole components at chiral sites */

    chkpole_();

/*     rotate the multipole components into the global frame */

    rotpole_();

/*     compute the indices used in reaction field calculations */

    ijkpts_();

/*     calculate the reaction field interaction energy term */

    i__1 = mpole_1.npole;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = mpole_1.ipole[ii - 1];
	iz = mpole_1.zaxis[ii - 1];
	ix = mpole_1.xaxis[ii - 1];
	usei = usage_1.use[i__ - 1] || usage_1.use[iz - 1] || usage_1.use[ix 
		- 1];
	i__2 = mpole_1.polsiz[ii - 1];
	for (j = 1; j <= i__2; ++j) {
	    rpi[j - 1] = rpole_ref(j, ii);
	}
	i__2 = mpole_1.npole;
	for (kk = ii; kk <= i__2; ++kk) {
	    k = mpole_1.ipole[kk - 1];
	    kz = mpole_1.zaxis[kk - 1];
	    kx = mpole_1.xaxis[kk - 1];
	    usek = usage_1.use[k - 1] || usage_1.use[kz - 1] || usage_1.use[
		    kx - 1];
	    if (usei || usek) {
		xr = atoms_1.x[k - 1] - atoms_1.x[i__ - 1];
		yr = atoms_1.y[k - 1] - atoms_1.y[i__ - 1];
		zr = atoms_1.z__[k - 1] - atoms_1.z__[i__ - 1];
		r2 = xr * xr + yr * yr + zr * zr;
		if (r2 <= shunt_1.off2) {
		    i__3 = mpole_1.polsiz[kk - 1];
		    for (j = 1; j <= i__3; ++j) {
			rpk[j - 1] = rpole_ref(j, kk);
		    }
		    erfik_(&ii, &kk, &i__, &k, rpi, rpk, &eik);
		    ++action_1.ner;
		    energi_1.er += eik;
		    analyz_1.aer[i__ - 1] += eik * .5;
		    analyz_1.aer[k - 1] += eik * .5;

/*     print a message if the energy of this interaction is large */

		    huge__ = eik > 10.;
		    if (inform_1.debug || inform_1.verbose && huge__) {
			if (header) {
			    header = FALSE_;
			    io___21.ciunit = iounit_1.iout;
			    s_wsfe(&io___21);
			    e_wsfe();
			}
			r__ = sqrt(r2);
/* Computing 2nd power */
			d__1 = atoms_1.x[i__ - 1];
/* Computing 2nd power */
			d__2 = atoms_1.y[i__ - 1];
/* Computing 2nd power */
			d__3 = atoms_1.z__[i__ - 1];
			di = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
			d__1 = atoms_1.x[k - 1];
/* Computing 2nd power */
			d__2 = atoms_1.y[k - 1];
/* Computing 2nd power */
			d__3 = atoms_1.z__[k - 1];
			dk = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
			io___25.ciunit = iounit_1.iout;
			s_wsfe(&io___25);
			do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
			do_fio(&c__1, name___ref(0, i__), (ftnlen)3);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			do_fio(&c__1, name___ref(0, k), (ftnlen)3);
			do_fio(&c__1, (char *)&di, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&dk, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&r__, (ftnlen)sizeof(doublereal)
				);
			do_fio(&c__1, (char *)&eik, (ftnlen)sizeof(doublereal)
				);
			e_wsfe();
		    }
		}
	    }
	}
    }
    return 0;
} /* erxnfld3_ */

#undef rpole_ref
#undef name___ref


