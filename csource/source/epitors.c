/* epitors.f -- translated by f2c (version 20050501).
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
    doublereal polycut, polycut2;
    logical use_bounds__, use_replica__, use_polymer__;
} bound_;

#define bound_1 bound_

struct {
    doublereal esum, eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, 
	    ebt, ett, ev, ec, ecd, ed, em, ep, er, es, elf, eg, ex;
} energi_;

#define energi_1 energi_

struct {
    doublereal grpmass[1000], wgrp[1002001]	/* was [1001][1001] */;
    integer ngrp, kgrp[25000], igrp[2002]	/* was [2][1001] */, grplist[
	    25000];
    logical use_group__, use_intra__, use_inter__;
} group_;

#define group_1 group_

struct {
    doublereal kpit[100000];
    integer npitors, ipit[600000]	/* was [6][100000] */;
} pitors_;

#define pitors_1 pitors_

struct {
    doublereal idihunit, itorunit, torsunit, ptorunit, storunit, ttorunit;
} torpot_;

#define torpot_1 torpot_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine epitors  --  pi-orbit torsion potential energy  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "epitors" calculates the pi-orbital torsion potential energy */


/* Subroutine */ int epitors_(void)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal e;
    static integer i__;
    static doublereal c2, s2, v2;
    static integer ia, ib, ic, id, ie, ig;
    static doublereal xt, yt, zt, xu, yu, zu, rt2, ru2, rdc, xad, yad, zad, 
	    xbd, ybd, xia, yia, zia, xib, yib, zib, xic, yic, zic, xid, yid, 
	    zid, xie, yie, zie, xig, yig, zig, xip, yip, zip, xiq, yiq, ziq, 
	    zbd, xec, yec, xtu, ytu, ztu, zec, xgc, ygc, zgc, xcp, ycp, zcp, 
	    xdc, ydc, zdc, xqd, yqd, zqd, phi2, fgrp, sine, rtru, sine2;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal cosine;
    extern /* Subroutine */ int groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal cosine2;
    static logical proceed;


#define ipit_ref(a_1,a_2) pitors_1.ipit[(a_2)*6 + a_1 - 7]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     zero out the pi-orbital torsion potential energy */

    energi_1.ept = 0.;

/*     calculate the pi-orbital torsion angle energy term */

    i__1 = pitors_1.npitors;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = ipit_ref(1, i__);
	ib = ipit_ref(2, i__);
	ic = ipit_ref(3, i__);
	id = ipit_ref(4, i__);
	ie = ipit_ref(5, i__);
	ig = ipit_ref(6, i__);

/*     decide whether to compute the current interaction */

	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &ie, &ig);
	}
	if (proceed) {
	    proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
		    usage_1.use[ic - 1] || usage_1.use[id - 1] || usage_1.use[
		    ie - 1] || usage_1.use[ig - 1];
	}

/*     compute the value of the pi-orbital torsion angle */

	if (proceed) {
	    xia = atoms_1.x[ia - 1];
	    yia = atoms_1.y[ia - 1];
	    zia = atoms_1.z__[ia - 1];
	    xib = atoms_1.x[ib - 1];
	    yib = atoms_1.y[ib - 1];
	    zib = atoms_1.z__[ib - 1];
	    xic = atoms_1.x[ic - 1];
	    yic = atoms_1.y[ic - 1];
	    zic = atoms_1.z__[ic - 1];
	    xid = atoms_1.x[id - 1];
	    yid = atoms_1.y[id - 1];
	    zid = atoms_1.z__[id - 1];
	    xie = atoms_1.x[ie - 1];
	    yie = atoms_1.y[ie - 1];
	    zie = atoms_1.z__[ie - 1];
	    xig = atoms_1.x[ig - 1];
	    yig = atoms_1.y[ig - 1];
	    zig = atoms_1.z__[ig - 1];
	    xad = xia - xid;
	    yad = yia - yid;
	    zad = zia - zid;
	    xbd = xib - xid;
	    ybd = yib - yid;
	    zbd = zib - zid;
	    xec = xie - xic;
	    yec = yie - yic;
	    zec = zie - zic;
	    xgc = xig - xic;
	    ygc = yig - yic;
	    zgc = zig - zic;
	    if (bound_1.use_polymer__) {
		image_(&xad, &yad, &zad);
		image_(&xbd, &ybd, &zbd);
		image_(&xec, &yec, &zec);
		image_(&xgc, &ygc, &zgc);
	    }
	    xip = yad * zbd - ybd * zad + xic;
	    yip = zad * xbd - zbd * xad + yic;
	    zip = xad * ybd - xbd * yad + zic;
	    xiq = yec * zgc - ygc * zec + xid;
	    yiq = zec * xgc - zgc * xec + yid;
	    ziq = xec * ygc - xgc * yec + zid;
	    xcp = xic - xip;
	    ycp = yic - yip;
	    zcp = zic - zip;
	    xdc = xid - xic;
	    ydc = yid - yic;
	    zdc = zid - zic;
	    xqd = xiq - xid;
	    yqd = yiq - yid;
	    zqd = ziq - zid;
	    if (bound_1.use_polymer__) {
		image_(&xcp, &ycp, &zcp);
		image_(&xdc, &ydc, &zdc);
		image_(&xqd, &yqd, &zqd);
	    }
	    xt = ycp * zdc - ydc * zcp;
	    yt = zcp * xdc - zdc * xcp;
	    zt = xcp * ydc - xdc * ycp;
	    xu = ydc * zqd - yqd * zdc;
	    yu = zdc * xqd - zqd * xdc;
	    zu = xdc * yqd - xqd * ydc;
	    xtu = yt * zu - yu * zt;
	    ytu = zt * xu - zu * xt;
	    ztu = xt * yu - xu * yt;
	    rt2 = xt * xt + yt * yt + zt * zt;
	    ru2 = xu * xu + yu * yu + zu * zu;
	    rtru = sqrt(rt2 * ru2);
	    if (rtru != 0.) {
		rdc = sqrt(xdc * xdc + ydc * ydc + zdc * zdc);
		cosine = (xt * xu + yt * yu + zt * zu) / rtru;
		sine = (xdc * xtu + ydc * ytu + zdc * ztu) / (rdc * rtru);

/*     set the pi-orbital torsion parameters for this angle */

		v2 = pitors_1.kpit[i__ - 1];
		c2 = -1.;
		s2 = 0.;

/*     compute the multiple angle trigonometry and the phase terms */

		cosine2 = cosine * cosine - sine * sine;
		sine2 = cosine * 2. * sine;
		phi2 = cosine2 * c2 + sine2 * s2 + 1.;

/*     calculate the pi-orbital torsion energy for this angle */

		e = torpot_1.ptorunit * v2 * phi2;

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    e *= fgrp;
		}

/*     increment the total pi-orbital torsion angle energy */

		energi_1.ept += e;
	    }
	}
    }
    return 0;
} /* epitors_ */

#undef ipit_ref


