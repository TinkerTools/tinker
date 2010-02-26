/* torque.f -- translated by f2c (version 20050501).
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
    doublereal desum[75000]	/* was [3][25000] */, deb[75000]	/* 
	    was [3][25000] */, dea[75000]	/* was [3][25000] */, deba[
	    75000]	/* was [3][25000] */, deub[75000]	/* was [3][
	    25000] */, deaa[75000]	/* was [3][25000] */, deopb[75000]	
	    /* was [3][25000] */, deopd[75000]	/* was [3][25000] */, deid[
	    75000]	/* was [3][25000] */, deit[75000]	/* was [3][
	    25000] */, det[75000]	/* was [3][25000] */, dept[75000]	
	    /* was [3][25000] */, debt[75000]	/* was [3][25000] */, dett[
	    75000]	/* was [3][25000] */, dev[75000]	/* was [3][
	    25000] */, dec[75000]	/* was [3][25000] */, decd[75000]	
	    /* was [3][25000] */, ded[75000]	/* was [3][25000] */, dem[
	    75000]	/* was [3][25000] */, dep[75000]	/* was [3][
	    25000] */, der[75000]	/* was [3][25000] */, des[75000]	
	    /* was [3][25000] */, delf[75000]	/* was [3][25000] */, deg[
	    75000]	/* was [3][25000] */, dex[75000]	/* was [3][
	    25000] */;
} deriv_;

#define deriv_1 deriv_

struct {
    doublereal pole[325000]	/* was [13][25000] */, rpole[325000]	/* 
	    was [13][25000] */;
    integer npole, ipole[25000], polsiz[25000], pollist[25000], zaxis[25000], 
	    xaxis[25000], yaxis[25000];
    char polaxe[200000];
} mpole_;

#define mpole_1 mpole_



/*     ############################################################# */
/*     ##  COPYRIGHT (C) 2007 by Pengyu Ren & Jay William Ponder  ## */
/*     ##                   All Rights Reserved                   ## */
/*     ############################################################# */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine torque  --  convert torque to Cartesian force  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "torque" takes the torque values on sites defined by local */
/*     coordinate frames and converts to Cartesian forces on the */
/*     original sites and sites specifying the local frames */

/*     literature reference: */

/*     P. L. Popelier and A. J. Stone, "Formulae for the First and */
/*     Second Derivatives of Aniostropic Potentials with Respect to */
/*     Geometrical Parameters", Molecular Physics, 82, 411-425 (1994) */


/* Subroutine */ int torque_(integer *i__, doublereal *trq1, doublereal *trq2,
	 doublereal *frcx, doublereal *frcy, doublereal *frcz)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer j;
    static doublereal r__[3], s[3], u[3], v[3], w[3], t1[3], t2[3];
    static integer ia, ib, ic, id;
    static doublereal du, dv, dw, ur[3], us[3], vs[3], ws[3], uv[3], uw[3], 
	    vw[3], rsiz, ssiz, usiz, vsiz, wsiz, t1siz, t2siz, urcos, uvcos, 
	    uwcos, vwcos, uscos, vscos, wscos, uvsin, uwsin, vwsin, ursin, 
	    ussin, vssin, wssin, ursiz, ussiz, vssiz, wssiz, uvsiz, uwsiz, 
	    vwsiz, ut1cos, ut2cos, ut1sin, ut2sin, dphidr, dphids, dphidu, 
	    dphidv, dphidw;
    static char axetyp[8];


#define dem_ref(a_1,a_2) deriv_1.dem[(a_2)*3 + a_1 - 4]
#define dep_ref(a_1,a_2) deriv_1.dep[(a_2)*3 + a_1 - 4]
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
/*     ##  deriv.i  --  Cartesian coordinate derivative components  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     desum   total energy Cartesian coordinate derivatives */
/*     deb     bond stretch Cartesian coordinate derivatives */
/*     dea     angle bend Cartesian coordinate derivatives */
/*     deba    stretch-bend Cartesian coordinate derivatives */
/*     deub    Urey-Bradley Cartesian coordinate derivatives */
/*     deaa    angle-angle Cartesian coordinate derivatives */
/*     deopb   out-of-plane bend Cartesian coordinate derivatives */
/*     deopd   out-of-plane distance Cartesian coordinate derivatives */
/*     deid    improper dihedral Cartesian coordinate derivatives */
/*     deit    improper torsion Cartesian coordinate derivatives */
/*     det     torsional Cartesian coordinate derivatives */
/*     dept    pi-orbital torsion Cartesian coordinate derivatives */
/*     debt    stretch-torsion Cartesian coordinate derivatives */
/*     dett    torsion-torsion Cartesian coordinate derivatives */
/*     dev     van der Waals Cartesian coordinate derivatives */
/*     dec     charge-charge Cartesian coordinate derivatives */
/*     decd    charge-dipole Cartesian coordinate derivatives */
/*     ded     dipole-dipole Cartesian coordinate derivatives */
/*     dem     multipole Cartesian coordinate derivatives */
/*     dep     polarization Cartesian coordinate derivatives */
/*     der     reaction field Cartesian coordinate derivatives */
/*     des     solvation Cartesian coordinate derivatives */
/*     delf    metal ligand field Cartesian coordinate derivatives */
/*     deg     geometric restraint Cartesian coordinate derivatives */
/*     dex     extra energy term Cartesian coordinate derivatives */




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




/*     get the local frame type and the frame-defining atoms */

    /* Parameter adjustments */
    --frcz;
    --frcy;
    --frcx;
    --trq2;
    --trq1;

    /* Function Body */
    ia = mpole_1.zaxis[*i__ - 1];
    ib = mpole_1.ipole[*i__ - 1];
    ic = mpole_1.xaxis[*i__ - 1];
    id = mpole_1.yaxis[*i__ - 1];
    s_copy(axetyp, polaxe_ref(0, *i__), (ftnlen)8, (ftnlen)8);

/*     zero out force components on local frame-defining atoms */

    for (j = 1; j <= 3; ++j) {
	frcz[j] = 0.;
	frcx[j] = 0.;
	frcy[j] = 0.;
    }

/*     construct the three rotation axes for the local frame */

    u[0] = atoms_1.x[ia - 1] - atoms_1.x[ib - 1];
    u[1] = atoms_1.y[ia - 1] - atoms_1.y[ib - 1];
    u[2] = atoms_1.z__[ia - 1] - atoms_1.z__[ib - 1];
    v[0] = atoms_1.x[ic - 1] - atoms_1.x[ib - 1];
    v[1] = atoms_1.y[ic - 1] - atoms_1.y[ib - 1];
    v[2] = atoms_1.z__[ic - 1] - atoms_1.z__[ib - 1];
    if (s_cmp(axetyp, "Z-then-X", (ftnlen)8, (ftnlen)8) == 0 || s_cmp(axetyp, 
	    "Bisector", (ftnlen)8, (ftnlen)8) == 0) {
	w[0] = u[1] * v[2] - u[2] * v[1];
	w[1] = u[2] * v[0] - u[0] * v[2];
	w[2] = u[0] * v[1] - u[1] * v[0];
    } else if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0 || s_cmp(
	    axetyp, "3-Fold", (ftnlen)8, (ftnlen)6) == 0) {
	w[0] = atoms_1.x[id - 1] - atoms_1.x[ib - 1];
	w[1] = atoms_1.y[id - 1] - atoms_1.y[ib - 1];
	w[2] = atoms_1.z__[id - 1] - atoms_1.z__[ib - 1];
    }
    usiz = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    vsiz = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    wsiz = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    for (j = 1; j <= 3; ++j) {
	u[j - 1] /= usiz;
	v[j - 1] /= vsiz;
	w[j - 1] /= wsiz;
    }

/*     build some additional axes needed for the Z-bisect method */

    if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	r__[0] = v[0] + w[0];
	r__[1] = v[1] + w[1];
	r__[2] = v[2] + w[2];
	s[0] = u[1] * r__[2] - u[2] * r__[1];
	s[1] = u[2] * r__[0] - u[0] * r__[2];
	s[2] = u[0] * r__[1] - u[1] * r__[0];
	rsiz = sqrt(r__[0] * r__[0] + r__[1] * r__[1] + r__[2] * r__[2]);
	ssiz = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
	for (j = 1; j <= 3; ++j) {
	    r__[j - 1] /= rsiz;
	    s[j - 1] /= ssiz;
	}
    }

/*     find the perpendicular and angle for each pair of axes */

    uv[0] = v[1] * u[2] - v[2] * u[1];
    uv[1] = v[2] * u[0] - v[0] * u[2];
    uv[2] = v[0] * u[1] - v[1] * u[0];
    uw[0] = w[1] * u[2] - w[2] * u[1];
    uw[1] = w[2] * u[0] - w[0] * u[2];
    uw[2] = w[0] * u[1] - w[1] * u[0];
    vw[0] = w[1] * v[2] - w[2] * v[1];
    vw[1] = w[2] * v[0] - w[0] * v[2];
    vw[2] = w[0] * v[1] - w[1] * v[0];
    uvsiz = sqrt(uv[0] * uv[0] + uv[1] * uv[1] + uv[2] * uv[2]);
    uwsiz = sqrt(uw[0] * uw[0] + uw[1] * uw[1] + uw[2] * uw[2]);
    vwsiz = sqrt(vw[0] * vw[0] + vw[1] * vw[1] + vw[2] * vw[2]);
    for (j = 1; j <= 3; ++j) {
	uv[j - 1] /= uvsiz;
	uw[j - 1] /= uwsiz;
	vw[j - 1] /= vwsiz;
    }
    if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	ur[0] = r__[1] * u[2] - r__[2] * u[1];
	ur[1] = r__[2] * u[0] - r__[0] * u[2];
	ur[2] = r__[0] * u[1] - r__[1] * u[0];
	us[0] = s[1] * u[2] - s[2] * u[1];
	us[1] = s[2] * u[0] - s[0] * u[2];
	us[2] = s[0] * u[1] - s[1] * u[0];
	vs[0] = s[1] * v[2] - s[2] * v[1];
	vs[1] = s[2] * v[0] - s[0] * v[2];
	vs[2] = s[0] * v[1] - s[1] * v[0];
	ws[0] = s[1] * w[2] - s[2] * w[1];
	ws[1] = s[2] * w[0] - s[0] * w[2];
	ws[2] = s[0] * w[1] - s[1] * w[0];
	ursiz = sqrt(ur[0] * ur[0] + ur[1] * ur[1] + ur[2] * ur[2]);
	ussiz = sqrt(us[0] * us[0] + us[1] * us[1] + us[2] * us[2]);
	vssiz = sqrt(vs[0] * vs[0] + vs[1] * vs[1] + vs[2] * vs[2]);
	wssiz = sqrt(ws[0] * ws[0] + ws[1] * ws[1] + ws[2] * ws[2]);
	for (j = 1; j <= 3; ++j) {
	    ur[j - 1] /= ursiz;
	    us[j - 1] /= ussiz;
	    vs[j - 1] /= vssiz;
	    ws[j - 1] /= wssiz;
	}
    }

/*     compute the sine of the angle between the rotation axes */

    uvcos = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    uvsin = sqrt(1. - uvcos * uvcos);
    uwcos = u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
    uwsin = sqrt(1. - uwcos * uwcos);
    vwcos = v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
    vwsin = sqrt(1. - vwcos * vwcos);
    if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	urcos = u[0] * r__[0] + u[1] * r__[1] + u[2] * r__[2];
	ursin = sqrt(1. - urcos * urcos);
	uscos = u[0] * s[0] + u[1] * s[1] + u[2] * s[2];
	ussin = sqrt(1. - uscos * uscos);
	vscos = v[0] * s[0] + v[1] * s[1] + v[2] * s[2];
	vssin = sqrt(1. - vscos * vscos);
	wscos = w[0] * s[0] + w[1] * s[1] + w[2] * s[2];
	wssin = sqrt(1. - wscos * wscos);
    }

/*     compute the projection of v and w onto the ru-plane */

    if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	for (j = 1; j <= 3; ++j) {
	    t1[j - 1] = v[j - 1] - s[j - 1] * vscos;
	    t2[j - 1] = w[j - 1] - s[j - 1] * wscos;
	}
	t1siz = sqrt(t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]);
	t2siz = sqrt(t2[0] * t2[0] + t2[1] * t2[1] + t2[2] * t2[2]);
	for (j = 1; j <= 3; ++j) {
	    t1[j - 1] /= t1siz;
	    t2[j - 1] /= t2siz;
	}
	ut1cos = u[0] * t1[0] + u[1] * t1[1] + u[2] * t1[2];
	ut1sin = sqrt(1. - ut1cos * ut1cos);
	ut2cos = u[0] * t2[0] + u[1] * t2[1] + u[2] * t2[2];
	ut2sin = sqrt(1. - ut2cos * ut2cos);
    }

/*     negative of dot product of torque with unit vectors gives */
/*     result of infinitesimal rotation along these vectors */

    dphidu = -trq1[1] * u[0] - trq1[2] * u[1] - trq1[3] * u[2];
    dphidv = -trq1[1] * v[0] - trq1[2] * v[1] - trq1[3] * v[2];
    dphidw = -trq1[1] * w[0] - trq1[2] * w[1] - trq1[3] * w[2];
    if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	dphidr = -trq1[1] * r__[0] - trq1[2] * r__[1] - trq1[3] * r__[2];
	dphids = -trq1[1] * s[0] - trq1[2] * s[1] - trq1[3] * s[2];
    }

/*     force distribution for the Z-then-X local coordinate method */

    if (s_cmp(axetyp, "Z-then-X", (ftnlen)8, (ftnlen)8) == 0) {
	for (j = 1; j <= 3; ++j) {
	    du = uv[j - 1] * dphidv / (usiz * uvsin) + uw[j - 1] * dphidw / 
		    usiz;
	    dv = -uv[j - 1] * dphidu / (vsiz * uvsin);
	    dem_ref(j, ia) = dem_ref(j, ia) + du;
	    dem_ref(j, ic) = dem_ref(j, ic) + dv;
	    dem_ref(j, ib) = dem_ref(j, ib) - du - dv;
	    frcz[j] += du;
	    frcx[j] += dv;
	}

/*     force distribution for the bisector local coordinate method */

    } else if (s_cmp(axetyp, "Bisector", (ftnlen)8, (ftnlen)8) == 0) {
	for (j = 1; j <= 3; ++j) {
	    du = uv[j - 1] * dphidv / (usiz * uvsin) + uw[j - 1] * .5 * 
		    dphidw / usiz;
	    dv = -uv[j - 1] * dphidu / (vsiz * uvsin) + vw[j - 1] * .5 * 
		    dphidw / vsiz;
	    dem_ref(j, ia) = dem_ref(j, ia) + du;
	    dem_ref(j, ic) = dem_ref(j, ic) + dv;
	    dem_ref(j, ib) = dem_ref(j, ib) - du - dv;
	    frcz[j] += du;
	    frcx[j] += dv;
	}

/*     force distribution for the Z-bisect local coordinate method */

    } else if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	for (j = 1; j <= 3; ++j) {
	    du = ur[j - 1] * dphidr / (usiz * ursin) + us[j - 1] * dphids / 
		    usiz;
	    dv = (vssin * s[j - 1] - vscos * t1[j - 1]) * dphidu / (vsiz * (
		    ut1sin + ut2sin));
	    dw = (wssin * s[j - 1] - wscos * t2[j - 1]) * dphidu / (wsiz * (
		    ut1sin + ut2sin));
	    dem_ref(j, ia) = dem_ref(j, ia) + du;
	    dem_ref(j, ic) = dem_ref(j, ic) + dv;
	    dem_ref(j, id) = dem_ref(j, id) + dw;
	    dem_ref(j, ib) = dem_ref(j, ib) - du - dv - dw;
	    frcz[j] += du;
	    frcx[j] += dv;
	    frcy[j] += dw;
	}

/*     force distribution for the 3-fold local coordinate method */
/*        (correct for uv, uw and vw angles all equal to 90) */

    } else if (s_cmp(axetyp, "3-Fold", (ftnlen)8, (ftnlen)6) == 0) {
	for (j = 1; j <= 3; ++j) {
	    du = uw[j - 1] * dphidw / (usiz * uwsin) + uv[j - 1] * dphidv / (
		    usiz * uvsin) - uw[j - 1] * dphidu / (usiz * uwsin) - uv[
		    j - 1] * dphidu / (usiz * uvsin);
	    dv = vw[j - 1] * dphidw / (vsiz * vwsin) - uv[j - 1] * dphidu / (
		    vsiz * uvsin) - vw[j - 1] * dphidv / (vsiz * vwsin) + uv[
		    j - 1] * dphidv / (vsiz * uvsin);
	    dw = -uw[j - 1] * dphidu / (wsiz * uwsin) - vw[j - 1] * dphidv / (
		    wsiz * vwsin) + uw[j - 1] * dphidw / (wsiz * uwsin) + vw[
		    j - 1] * dphidw / (wsiz * vwsin);
	    du /= 3.;
	    dv /= 3.;
	    dw /= 3.;
	    dem_ref(j, ia) = dem_ref(j, ia) + du;
	    dem_ref(j, ic) = dem_ref(j, ic) + dv;
	    dem_ref(j, id) = dem_ref(j, id) + dw;
	    dem_ref(j, ib) = dem_ref(j, ib) - du - dv - dw;
	    frcz[j] += du;
	    frcx[j] += dv;
	    frcy[j] += dw;
	}
    }

/*     negative of dot product of torque with unit vectors gives */
/*     result of infinitesimal rotation along these vectors */

    dphidu = -trq2[1] * u[0] - trq2[2] * u[1] - trq2[3] * u[2];
    dphidv = -trq2[1] * v[0] - trq2[2] * v[1] - trq2[3] * v[2];
    dphidw = -trq2[1] * w[0] - trq2[2] * w[1] - trq2[3] * w[2];
    if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	dphidr = -trq2[1] * r__[0] - trq2[2] * r__[1] - trq2[3] * r__[2];
	dphids = -trq2[1] * s[0] - trq2[2] * s[1] - trq2[3] * s[2];
    }

/*     force distribution for the Z-then-X local coordinate method */

    if (s_cmp(axetyp, "Z-then-X", (ftnlen)8, (ftnlen)8) == 0) {
	for (j = 1; j <= 3; ++j) {
	    du = uv[j - 1] * dphidv / (usiz * uvsin) + uw[j - 1] * dphidw / 
		    usiz;
	    dv = -uv[j - 1] * dphidu / (vsiz * uvsin);
	    dep_ref(j, ia) = dep_ref(j, ia) + du;
	    dep_ref(j, ic) = dep_ref(j, ic) + dv;
	    dep_ref(j, ib) = dep_ref(j, ib) - du - dv;
	    frcz[j] += du;
	    frcx[j] += dv;
	}

/*     force distribution for the bisector local coordinate method */

    } else if (s_cmp(axetyp, "Bisector", (ftnlen)8, (ftnlen)8) == 0) {
	for (j = 1; j <= 3; ++j) {
	    du = uv[j - 1] * dphidv / (usiz * uvsin) + uw[j - 1] * .5 * 
		    dphidw / usiz;
	    dv = -uv[j - 1] * dphidu / (vsiz * uvsin) + vw[j - 1] * .5 * 
		    dphidw / vsiz;
	    dep_ref(j, ia) = dep_ref(j, ia) + du;
	    dep_ref(j, ic) = dep_ref(j, ic) + dv;
	    dep_ref(j, ib) = dep_ref(j, ib) - du - dv;
	    frcz[j] += du;
	    frcx[j] += dv;
	}

/*     force distribution for the Z-bisect local coordinate method */

    } else if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	for (j = 1; j <= 3; ++j) {
	    du = ur[j - 1] * dphidr / (usiz * ursin) + us[j - 1] * dphids / 
		    usiz;
	    dv = (vssin * s[j - 1] - vscos * t1[j - 1]) * dphidu / (vsiz * (
		    ut1sin + ut2sin));
	    dw = (wssin * s[j - 1] - wscos * t2[j - 1]) * dphidu / (wsiz * (
		    ut1sin + ut2sin));
	    dep_ref(j, ia) = dep_ref(j, ia) + du;
	    dep_ref(j, ic) = dep_ref(j, ic) + dv;
	    dep_ref(j, id) = dep_ref(j, id) + dw;
	    dep_ref(j, ib) = dep_ref(j, ib) - du - dv - dw;
	    frcz[j] += du;
	    frcx[j] += dv;
	    frcy[j] += dw;
	}

/*     force distribution for the 3-fold local coordinate method */
/*        (correct for uv, uw and vw angles all equal to 90) */

    } else if (s_cmp(axetyp, "3-Fold", (ftnlen)8, (ftnlen)6) == 0) {
	for (j = 1; j <= 3; ++j) {
	    du = uw[j - 1] * dphidw / (usiz * uwsin) + uv[j - 1] * dphidv / (
		    usiz * uvsin) - uw[j - 1] * dphidu / (usiz * uwsin) - uv[
		    j - 1] * dphidu / (usiz * uvsin);
	    dv = vw[j - 1] * dphidw / (vsiz * vwsin) - uv[j - 1] * dphidu / (
		    vsiz * uvsin) - vw[j - 1] * dphidv / (vsiz * vwsin) + uv[
		    j - 1] * dphidv / (vsiz * uvsin);
	    dw = -uw[j - 1] * dphidu / (wsiz * uwsin) - vw[j - 1] * dphidv / (
		    wsiz * vwsin) + uw[j - 1] * dphidw / (wsiz * uwsin) + vw[
		    j - 1] * dphidw / (wsiz * vwsin);
	    du /= 3.;
	    dv /= 3.;
	    dw /= 3.;
	    dep_ref(j, ia) = dep_ref(j, ia) + du;
	    dep_ref(j, ic) = dep_ref(j, ic) + dv;
	    dep_ref(j, id) = dep_ref(j, id) + dw;
	    dep_ref(j, ib) = dep_ref(j, ib) - du - dv - dw;
	    frcz[j] += du;
	    frcx[j] += dv;
	    frcy[j] += dw;
	}
    }
    return 0;
} /* torque_ */

#undef polaxe_ref
#undef dep_ref
#undef dem_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine torque2  --  convert torque to Cartesian force  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "torque2" takes the torque values on sites defined by local */
/*     coordinate frames and converts to Cartesian forces on the */
/*     original sites and sites specifying the local frames */

/*     literature reference: */

/*     P. L. Popelier and A. J. Stone, "Formulae for the First and */
/*     Second Derivatives of Aniostropic Potentials with Respect to */
/*     Geometrical Parameters", Molecular Physics, 82, 411-425 (1994) */


/* Subroutine */ int torque2_(doublereal *trq, doublereal *derivs)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal r__[3], s[3], u[3], v[3], w[3], t1[3], t2[3];
    static integer ia, ib, ic, id;
    static doublereal du, dv, dw, ur[3], us[3], vs[3], ws[3], uv[3], uw[3], 
	    vw[3], rsiz, ssiz, usiz, vsiz, wsiz, t1siz, t2siz, urcos, uvcos, 
	    uwcos, vwcos, uscos, vscos, wscos, uvsin, uwsin, vwsin, ursin, 
	    ussin, vssin, wssin, ursiz, ussiz, vssiz, wssiz, uvsiz, uwsiz, 
	    vwsiz, ut1cos, ut2cos, ut1sin, ut2sin, dphidr, dphids, dphidu, 
	    dphidv, dphidw;
    static char axetyp[8];


#define trq_ref(a_1,a_2) trq[(a_2)*3 + a_1]
#define polaxe_ref(a_0,a_1) &mpole_1.polaxe[(a_1)*8 + a_0 - 8]
#define derivs_ref(a_1,a_2) derivs[(a_2)*3 + a_1]



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




/*     get the local frame type and the frame-defining atoms */

    /* Parameter adjustments */
    derivs -= 4;
    trq -= 4;

    /* Function Body */
    i__1 = mpole_1.npole;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia = mpole_1.zaxis[i__ - 1];
	ib = mpole_1.ipole[i__ - 1];
	ic = mpole_1.xaxis[i__ - 1];
	id = mpole_1.yaxis[i__ - 1];
	s_copy(axetyp, polaxe_ref(0, i__), (ftnlen)8, (ftnlen)8);

/*     construct the three rotation axes for the local frame */

	u[0] = atoms_1.x[ia - 1] - atoms_1.x[ib - 1];
	u[1] = atoms_1.y[ia - 1] - atoms_1.y[ib - 1];
	u[2] = atoms_1.z__[ia - 1] - atoms_1.z__[ib - 1];
	v[0] = atoms_1.x[ic - 1] - atoms_1.x[ib - 1];
	v[1] = atoms_1.y[ic - 1] - atoms_1.y[ib - 1];
	v[2] = atoms_1.z__[ic - 1] - atoms_1.z__[ib - 1];
	if (s_cmp(axetyp, "Z-then-X", (ftnlen)8, (ftnlen)8) == 0 || s_cmp(
		axetyp, "Bisector", (ftnlen)8, (ftnlen)8) == 0) {
	    w[0] = u[1] * v[2] - u[2] * v[1];
	    w[1] = u[2] * v[0] - u[0] * v[2];
	    w[2] = u[0] * v[1] - u[1] * v[0];
	} else if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0 || 
		s_cmp(axetyp, "3-Fold", (ftnlen)8, (ftnlen)6) == 0) {
	    w[0] = atoms_1.x[id - 1] - atoms_1.x[ib - 1];
	    w[1] = atoms_1.y[id - 1] - atoms_1.y[ib - 1];
	    w[2] = atoms_1.z__[id - 1] - atoms_1.z__[ib - 1];
	}
	usiz = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
	vsiz = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	wsiz = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
	for (j = 1; j <= 3; ++j) {
	    u[j - 1] /= usiz;
	    v[j - 1] /= vsiz;
	    w[j - 1] /= wsiz;
	}

/*     build some additional axes needed for the Z-bisect method */

	if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	    r__[0] = v[0] + w[0];
	    r__[1] = v[1] + w[1];
	    r__[2] = v[2] + w[2];
	    s[0] = u[1] * r__[2] - u[2] * r__[1];
	    s[1] = u[2] * r__[0] - u[0] * r__[2];
	    s[2] = u[0] * r__[1] - u[1] * r__[0];
	    rsiz = sqrt(r__[0] * r__[0] + r__[1] * r__[1] + r__[2] * r__[2]);
	    ssiz = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
	    for (j = 1; j <= 3; ++j) {
		r__[j - 1] /= rsiz;
		s[j - 1] /= ssiz;
	    }
	}

/*     find the perpendicular and angle for each pair of axes */

	uv[0] = v[1] * u[2] - v[2] * u[1];
	uv[1] = v[2] * u[0] - v[0] * u[2];
	uv[2] = v[0] * u[1] - v[1] * u[0];
	uw[0] = w[1] * u[2] - w[2] * u[1];
	uw[1] = w[2] * u[0] - w[0] * u[2];
	uw[2] = w[0] * u[1] - w[1] * u[0];
	vw[0] = w[1] * v[2] - w[2] * v[1];
	vw[1] = w[2] * v[0] - w[0] * v[2];
	vw[2] = w[0] * v[1] - w[1] * v[0];
	uvsiz = sqrt(uv[0] * uv[0] + uv[1] * uv[1] + uv[2] * uv[2]);
	uwsiz = sqrt(uw[0] * uw[0] + uw[1] * uw[1] + uw[2] * uw[2]);
	vwsiz = sqrt(vw[0] * vw[0] + vw[1] * vw[1] + vw[2] * vw[2]);
	for (j = 1; j <= 3; ++j) {
	    uv[j - 1] /= uvsiz;
	    uw[j - 1] /= uwsiz;
	    vw[j - 1] /= vwsiz;
	}
	if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	    ur[0] = r__[1] * u[2] - r__[2] * u[1];
	    ur[1] = r__[2] * u[0] - r__[0] * u[2];
	    ur[2] = r__[0] * u[1] - r__[1] * u[0];
	    us[0] = s[1] * u[2] - s[2] * u[1];
	    us[1] = s[2] * u[0] - s[0] * u[2];
	    us[2] = s[0] * u[1] - s[1] * u[0];
	    vs[0] = s[1] * v[2] - s[2] * v[1];
	    vs[1] = s[2] * v[0] - s[0] * v[2];
	    vs[2] = s[0] * v[1] - s[1] * v[0];
	    ws[0] = s[1] * w[2] - s[2] * w[1];
	    ws[1] = s[2] * w[0] - s[0] * w[2];
	    ws[2] = s[0] * w[1] - s[1] * w[0];
	    ursiz = sqrt(ur[0] * ur[0] + ur[1] * ur[1] + ur[2] * ur[2]);
	    ussiz = sqrt(us[0] * us[0] + us[1] * us[1] + us[2] * us[2]);
	    vssiz = sqrt(vs[0] * vs[0] + vs[1] * vs[1] + vs[2] * vs[2]);
	    wssiz = sqrt(ws[0] * ws[0] + ws[1] * ws[1] + ws[2] * ws[2]);
	    for (j = 1; j <= 3; ++j) {
		ur[j - 1] /= ursiz;
		us[j - 1] /= ussiz;
		vs[j - 1] /= vssiz;
		ws[j - 1] /= wssiz;
	    }
	}

/*     compute the sine of the angle between the rotation axes */

	uvcos = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
	uvsin = sqrt(1. - uvcos * uvcos);
	uwcos = u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
	uwsin = sqrt(1. - uwcos * uwcos);
	vwcos = v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
	vwsin = sqrt(1. - vwcos * vwcos);
	if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	    urcos = u[0] * r__[0] + u[1] * r__[1] + u[2] * r__[2];
	    ursin = sqrt(1. - urcos * urcos);
	    uscos = u[0] * s[0] + u[1] * s[1] + u[2] * s[2];
	    ussin = sqrt(1. - uscos * uscos);
	    vscos = v[0] * s[0] + v[1] * s[1] + v[2] * s[2];
	    vssin = sqrt(1. - vscos * vscos);
	    wscos = w[0] * s[0] + w[1] * s[1] + w[2] * s[2];
	    wssin = sqrt(1. - wscos * wscos);
	}

/*     compute the projection of v and w onto the ru-plane */

	if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	    for (j = 1; j <= 3; ++j) {
		t1[j - 1] = v[j - 1] - s[j - 1] * vscos;
		t2[j - 1] = w[j - 1] - s[j - 1] * wscos;
	    }
	    t1siz = sqrt(t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]);
	    t2siz = sqrt(t2[0] * t2[0] + t2[1] * t2[1] + t2[2] * t2[2]);
	    for (j = 1; j <= 3; ++j) {
		t1[j - 1] /= t1siz;
		t2[j - 1] /= t2siz;
	    }
	    ut1cos = u[0] * t1[0] + u[1] * t1[1] + u[2] * t1[2];
	    ut1sin = sqrt(1. - ut1cos * ut1cos);
	    ut2cos = u[0] * t2[0] + u[1] * t2[1] + u[2] * t2[2];
	    ut2sin = sqrt(1. - ut2cos * ut2cos);
	}

/*     negative of dot product of torque with unit vectors gives */
/*     result of infinitesimal rotation along these vectors */

	dphidu = -trq_ref(1, i__) * u[0] - trq_ref(2, i__) * u[1] - trq_ref(3,
		 i__) * u[2];
	dphidv = -trq_ref(1, i__) * v[0] - trq_ref(2, i__) * v[1] - trq_ref(3,
		 i__) * v[2];
	dphidw = -trq_ref(1, i__) * w[0] - trq_ref(2, i__) * w[1] - trq_ref(3,
		 i__) * w[2];
	if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	    dphidr = -trq_ref(1, i__) * r__[0] - trq_ref(2, i__) * r__[1] - 
		    trq_ref(3, i__) * r__[2];
	    dphids = -trq_ref(1, i__) * s[0] - trq_ref(2, i__) * s[1] - 
		    trq_ref(3, i__) * s[2];
	}

/*     force distribution for the Z-then-X local coordinate method */

	if (s_cmp(axetyp, "Z-then-X", (ftnlen)8, (ftnlen)8) == 0) {
	    for (j = 1; j <= 3; ++j) {
		du = uv[j - 1] * dphidv / (usiz * uvsin) + uw[j - 1] * dphidw 
			/ usiz;
		dv = -uv[j - 1] * dphidu / (vsiz * uvsin);
		derivs_ref(j, ia) = derivs_ref(j, ia) + du;
		derivs_ref(j, ic) = derivs_ref(j, ic) + dv;
		derivs_ref(j, ib) = derivs_ref(j, ib) - du - dv;
	    }

/*     force distribution for the bisector local coordinate method */

	} else if (s_cmp(axetyp, "Bisector", (ftnlen)8, (ftnlen)8) == 0) {
	    for (j = 1; j <= 3; ++j) {
		du = uv[j - 1] * dphidv / (usiz * uvsin) + uw[j - 1] * .5 * 
			dphidw / usiz;
		dv = -uv[j - 1] * dphidu / (vsiz * uvsin) + vw[j - 1] * .5 * 
			dphidw / vsiz;
		derivs_ref(j, ia) = derivs_ref(j, ia) + du;
		derivs_ref(j, ic) = derivs_ref(j, ic) + dv;
		derivs_ref(j, ib) = derivs_ref(j, ib) - du - dv;
	    }

/*     force distribution for the Z-bisect local coordinate method */

	} else if (s_cmp(axetyp, "Z-Bisect", (ftnlen)8, (ftnlen)8) == 0) {
	    for (j = 1; j <= 3; ++j) {
		du = ur[j - 1] * dphidr / (usiz * ursin) + us[j - 1] * dphids 
			/ usiz;
		dv = (vssin * s[j - 1] - vscos * t1[j - 1]) * dphidu / (vsiz *
			 (ut1sin + ut2sin));
		dw = (wssin * s[j - 1] - wscos * t2[j - 1]) * dphidu / (wsiz *
			 (ut1sin + ut2sin));
		derivs_ref(j, ia) = derivs_ref(j, ia) + du;
		derivs_ref(j, ic) = derivs_ref(j, ic) + dv;
		derivs_ref(j, id) = derivs_ref(j, id) + dw;
		derivs_ref(j, ib) = derivs_ref(j, ib) - du - dv - dw;
	    }

/*     force distribution for the 3-fold local coordinate method */
/*        (correct for uv, uw and vw angles all equal to 90) */

	} else if (s_cmp(axetyp, "3-Fold", (ftnlen)8, (ftnlen)6) == 0) {
	    for (j = 1; j <= 3; ++j) {
		du = uw[j - 1] * dphidw / (usiz * uwsin) + uv[j - 1] * dphidv 
			/ (usiz * uvsin) - uw[j - 1] * dphidu / (usiz * uwsin)
			 - uv[j - 1] * dphidu / (usiz * uvsin);
		dv = vw[j - 1] * dphidw / (vsiz * vwsin) - uv[j - 1] * dphidu 
			/ (vsiz * uvsin) - vw[j - 1] * dphidv / (vsiz * vwsin)
			 + uv[j - 1] * dphidv / (vsiz * uvsin);
		dw = -uw[j - 1] * dphidu / (wsiz * uwsin) - vw[j - 1] * 
			dphidv / (wsiz * vwsin) + uw[j - 1] * dphidw / (wsiz *
			 uwsin) + vw[j - 1] * dphidw / (wsiz * vwsin);
		du /= 3.;
		dv /= 3.;
		dw /= 3.;
		derivs_ref(j, ia) = derivs_ref(j, ia) + du;
		derivs_ref(j, ic) = derivs_ref(j, ic) + dv;
		derivs_ref(j, id) = derivs_ref(j, id) + dw;
		derivs_ref(j, ib) = derivs_ref(j, ib) - du - dv - dw;
	    }
	}
    }
    return 0;
} /* torque2_ */

#undef derivs_ref
#undef polaxe_ref
#undef trq_ref


