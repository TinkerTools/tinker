/* etortor.f -- translated by f2c (version 20050501).
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
    integer nbitor, ibitor[500000]	/* was [5][100000] */;
} bitor_;

#define bitor_1 bitor_

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
    doublereal ttx[3000]	/* was [30][100] */, tty[3000]	/* was [30][
	    100] */, tbf[90000]	/* was [900][100] */, tbx[90000]	/* 
	    was [900][100] */, tby[90000]	/* was [900][100] */, tbxy[
	    90000]	/* was [900][100] */;
    integer tnx[100], tny[100];
    char ktt[2000];
} ktrtor_;

#define ktrtor_1 ktrtor_

struct {
    doublereal idihunit, itorunit, torsunit, ptorunit, storunit, ttorunit;
} torpot_;

#define torpot_1 torpot_

struct {
    integer ntortor, itt[300000]	/* was [3][100000] */;
} tortor_;

#define tortor_1 tortor_

struct {
    integer nuse, iuse[25000];
    logical use[25000];
} usage_;

#define usage_1 usage_

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

/* Table of constant values */

static integer c__0 = 0;



/*     ############################################################# */
/*     ##  COPYRIGHT (C) 2003 by Pengyu Ren & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine etortor  --  torsion-torsion cross term energy  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "etortor" calculates the torsion-torsion potential energy */


/* Subroutine */ int etortor_(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal);

    /* Local variables */
    static doublereal e;
    static integer i__, k, ia, ib, ic, id, ie;
    static doublereal xh, yh;
    static integer nt;
    static doublereal xt, yt, zt, xu, yu, zu, xv, yv, zv, ft1[4], ft2[4], x1l,
	     y1l, rt2, ru2, rv2, x1u, y1u, xba;
    static integer nhi;
    static doublereal yba, zba, xia, yia, zia, xib, yib, zib;
    static integer nlo;
    static doublereal xic, yic, zic;
    static integer xlo, ylo;
    static doublereal xid, yid, zid, xie, yie, zie, xdc, ydc, xtu, ytu, ztu, 
	    xuv, yuv, zuv, zdc, xcb, ycb, zcb, xed, yed, zed, ftt[4], ft12[4];
    static integer pos1, pos2;
    static doublereal fgrp, sign, rtru, rurv;
    extern /* Subroutine */ int image_(doublereal *, doublereal *, doublereal 
	    *);
    static doublereal angle1, angle2, value1, value2;
    extern /* Subroutine */ int bcuint_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), groups_(logical *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal cosine1, cosine2;
    static logical proceed;
    extern /* Subroutine */ int chkttor_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    static integer itortor;


#define tbf_ref(a_1,a_2) ktrtor_1.tbf[(a_2)*900 + a_1 - 901]
#define tbx_ref(a_1,a_2) ktrtor_1.tbx[(a_2)*900 + a_1 - 901]
#define tby_ref(a_1,a_2) ktrtor_1.tby[(a_2)*900 + a_1 - 901]
#define itt_ref(a_1,a_2) tortor_1.itt[(a_2)*3 + a_1 - 4]
#define ttx_ref(a_1,a_2) ktrtor_1.ttx[(a_2)*30 + a_1 - 31]
#define tty_ref(a_1,a_2) ktrtor_1.tty[(a_2)*30 + a_1 - 31]
#define tbxy_ref(a_1,a_2) ktrtor_1.tbxy[(a_2)*900 + a_1 - 901]
#define ibitor_ref(a_1,a_2) bitor_1.ibitor[(a_2)*5 + a_1 - 6]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  usage.i  --  atoms active during energy computation  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     nuse   total number of active atoms in energy calculation */
/*     iuse   numbers of the atoms active in energy calculation */
/*     use    true if an atom is active, false if inactive */




/*     zero out the torsion-torsion energy */

    energi_1.ett = 0.;

/*     calculate the torsion-torsion interaction energy term */

    i__1 = tortor_1.ntortor;
    for (itortor = 1; itortor <= i__1; ++itortor) {
	i__ = itt_ref(1, itortor);
	k = itt_ref(2, itortor);
	if (itt_ref(3, itortor) == 1) {
	    ia = ibitor_ref(1, i__);
	    ib = ibitor_ref(2, i__);
	    ic = ibitor_ref(3, i__);
	    id = ibitor_ref(4, i__);
	    ie = ibitor_ref(5, i__);
	} else {
	    ia = ibitor_ref(5, i__);
	    ib = ibitor_ref(4, i__);
	    ic = ibitor_ref(3, i__);
	    id = ibitor_ref(2, i__);
	    ie = ibitor_ref(1, i__);
	}

/*     decide whether to compute the current interaction */

	proceed = TRUE_;
	if (group_1.use_group__) {
	    groups_(&proceed, &fgrp, &ia, &ib, &ic, &id, &ie, &c__0);
	}
	if (proceed) {
	    proceed = usage_1.use[ia - 1] || usage_1.use[ib - 1] || 
		    usage_1.use[ic - 1] || usage_1.use[id - 1] || usage_1.use[
		    ie - 1];
	}

/*     compute the values of the torsional angles */

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
	    xba = xib - xia;
	    yba = yib - yia;
	    zba = zib - zia;
	    xcb = xic - xib;
	    ycb = yic - yib;
	    zcb = zic - zib;
	    xdc = xid - xic;
	    ydc = yid - yic;
	    zdc = zid - zic;
	    xed = xie - xid;
	    yed = yie - yid;
	    zed = zie - zid;
	    if (bound_1.use_polymer__) {
		image_(&xba, &yba, &zba);
		image_(&xcb, &ycb, &zcb);
		image_(&xdc, &ydc, &zdc);
		image_(&xed, &yed, &zed);
	    }
	    xt = yba * zcb - ycb * zba;
	    yt = zba * xcb - zcb * xba;
	    zt = xba * ycb - xcb * yba;
	    xu = ycb * zdc - ydc * zcb;
	    yu = zcb * xdc - zdc * xcb;
	    zu = xcb * ydc - xdc * ycb;
	    xtu = yt * zu - yu * zt;
	    ytu = zt * xu - zu * xt;
	    ztu = xt * yu - xu * yt;
	    rt2 = xt * xt + yt * yt + zt * zt;
	    ru2 = xu * xu + yu * yu + zu * zu;
	    rtru = sqrt(rt2 * ru2);
	    xv = ydc * zed - yed * zdc;
	    yv = zdc * xed - zed * xdc;
	    zv = xdc * yed - xed * ydc;
	    xuv = yu * zv - yv * zu;
	    yuv = zu * xv - zv * xu;
	    zuv = xu * yv - xv * yu;
	    rv2 = xv * xv + yv * yv + zv * zv;
	    rurv = sqrt(ru2 * rv2);
	    if (rtru != 0. && rurv != 0.) {
		cosine1 = (xt * xu + yt * yu + zt * zu) / rtru;
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine1);
		cosine1 = min(d__1,d__2);
		angle1 = acos(cosine1) * 57.29577951308232088;
		sign = xba * xu + yba * yu + zba * zu;
		if (sign < 0.) {
		    angle1 = -angle1;
		}
		value1 = angle1;
		cosine2 = (xu * xv + yu * yv + zu * zv) / rurv;
/* Computing MIN */
		d__1 = 1., d__2 = max(-1.,cosine2);
		cosine2 = min(d__1,d__2);
		angle2 = acos(cosine2) * 57.29577951308232088;
		sign = xcb * xv + ycb * yv + zcb * zv;
		if (sign < 0.) {
		    angle2 = -angle2;
		}
		value2 = angle2;

/*     check for inverted chirality at the central atom */

		chkttor_(&ib, &ic, &id, &sign, &value1, &value2);

/*     use bicubic interpolation to compute spline values */

		nlo = 1;
		nhi = ktrtor_1.tnx[k - 1];
		while(nhi - nlo > 1) {
		    nt = (nhi + nlo) / 2;
		    if (ttx_ref(nt, k) > value1) {
			nhi = nt;
		    } else {
			nlo = nt;
		    }
		}
		xlo = nlo;
		nlo = 1;
		nhi = ktrtor_1.tny[k - 1];
		while(nhi - nlo > 1) {
		    nt = (nhi + nlo) / 2;
		    if (tty_ref(nt, k) > value2) {
			nhi = nt;
		    } else {
			nlo = nt;
		    }
		}
		ylo = nlo;
		xh = ttx_ref(xlo + 1, k) - ttx_ref(xlo, k);
		yh = tty_ref(ylo + 1, k) - tty_ref(ylo, k);
		x1l = ttx_ref(xlo, k);
		x1u = ttx_ref(xlo + 1, k);
		y1l = tty_ref(ylo, k);
		y1u = tty_ref(ylo + 1, k);
		pos2 = ylo * ktrtor_1.tnx[k - 1] + xlo;
		pos1 = pos2 - ktrtor_1.tnx[k - 1];
		ftt[0] = tbf_ref(pos1, k);
		ftt[1] = tbf_ref(pos1 + 1, k);
		ftt[2] = tbf_ref(pos2 + 1, k);
		ftt[3] = tbf_ref(pos2, k);
		ft1[0] = tbx_ref(pos1, k);
		ft1[1] = tbx_ref(pos1 + 1, k);
		ft1[2] = tbx_ref(pos2 + 1, k);
		ft1[3] = tbx_ref(pos2, k);
		ft2[0] = tby_ref(pos1, k);
		ft2[1] = tby_ref(pos1 + 1, k);
		ft2[2] = tby_ref(pos2 + 1, k);
		ft2[3] = tby_ref(pos2, k);
		ft12[0] = tbxy_ref(pos1, k);
		ft12[1] = tbxy_ref(pos1 + 1, k);
		ft12[2] = tbxy_ref(pos2 + 1, k);
		ft12[3] = tbxy_ref(pos2, k);
		bcuint_(ftt, ft1, ft2, ft12, &x1l, &x1u, &y1l, &y1u, &value1, 
			&value2, &e);
		e = torpot_1.ttorunit * e;

/*     scale the interaction based on its group membership */

		if (group_1.use_group__) {
		    e *= fgrp;
		}

/*     increment the total torsion-torsion energy */

		energi_1.ett += e;
	    }
	}
    }
    return 0;
} /* etortor_ */

#undef ibitor_ref
#undef tbxy_ref
#undef tty_ref
#undef ttx_ref
#undef itt_ref
#undef tby_ref
#undef tbx_ref
#undef tbf_ref




/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine chkttor  --  check torsion-torsion chirality  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "chkttor" tests the attached atoms at a torsion-torsion central */
/*     site and inverts the angle values if the site is chiral */

/*     note that the sign convention used in this version is correct */
/*     for phi-psi torsion-torsion interactions as defined in the */
/*     AMOEBA protein force field; the code may need to be altered */
/*     for other uses of the torsion-torsion potential, and will not */
/*     correctly handle enantiomeric sugar rings in nucleic acids */


/* Subroutine */ int chkttor_(integer *ib, integer *ic, integer *id, 
	doublereal *sign, doublereal *value1, doublereal *value2)
{
    static integer i__, j, k, m;
    static doublereal c1, c2, c3;
    static integer ia;
    static doublereal xac, yac, zac, xbc, ybc, zbc, xdc, ydc, zdc, vol;


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




/*     test for chirality at the central torsion-torsion site */

    *sign = 1.;
    if (couple_1.n12[*ic - 1] == 4) {
	j = 0;
	for (i__ = 1; i__ <= 4; ++i__) {
	    m = i12_ref(i__, *ic);
	    if (m != *ib && m != *id) {
		if (j == 0) {
		    j = m;
		} else {
		    k = m;
		}
	    }
	}
	ia = 0;
	if (atoms_1.type__[j - 1] > atoms_1.type__[k - 1]) {
	    ia = j;
	}
	if (atoms_1.type__[k - 1] > atoms_1.type__[j - 1]) {
	    ia = k;
	}
	if (atmtyp_1.atomic[j - 1] > atmtyp_1.atomic[k - 1]) {
	    ia = j;
	}
	if (atmtyp_1.atomic[k - 1] > atmtyp_1.atomic[j - 1]) {
	    ia = k;
	}

/*     compute the signed parallelpiped volume at central site */

	if (ia != 0) {
	    xac = atoms_1.x[ia - 1] - atoms_1.x[*ic - 1];
	    yac = atoms_1.y[ia - 1] - atoms_1.y[*ic - 1];
	    zac = atoms_1.z__[ia - 1] - atoms_1.z__[*ic - 1];
	    xbc = atoms_1.x[*ib - 1] - atoms_1.x[*ic - 1];
	    ybc = atoms_1.y[*ib - 1] - atoms_1.y[*ic - 1];
	    zbc = atoms_1.z__[*ib - 1] - atoms_1.z__[*ic - 1];
	    xdc = atoms_1.x[*id - 1] - atoms_1.x[*ic - 1];
	    ydc = atoms_1.y[*id - 1] - atoms_1.y[*ic - 1];
	    zdc = atoms_1.z__[*id - 1] - atoms_1.z__[*ic - 1];
	    c1 = ybc * zdc - zbc * ydc;
	    c2 = ydc * zac - zdc * yac;
	    c3 = yac * zbc - zac * ybc;
	    vol = xac * c1 + xbc * c2 + xdc * c3;

/*     invert the angle values if chirality has an inverted sign */

	    if (vol < 0.) {
		*sign = -1.;
		*value1 = -(*value1);
		*value2 = -(*value2);
	    }
	}
    }
    return 0;
} /* chkttor_ */

#undef i12_ref


