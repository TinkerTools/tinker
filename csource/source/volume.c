/* volume.f -- translated by f2c (version 20050501).
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
    integer iout, input;
} iounit_;

#define iounit_1 iounit_



/*     ################################################################ */
/*     ##  COPYRIGHT (C) 1990 by Craig Kundrot & Jay William Ponder  ## */
/*     ##                    All Rights Reserved                     ## */
/*     ################################################################ */

/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine volume  --  excluded volume term via Connolly  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "volume" calculates the excluded volume via the Connolly */
/*     analytical volume and surface area algorithm */


/* Subroutine */ int volume_(doublereal *volume_tot__, doublereal *radius, 
	doublereal *exclude)
{
    static doublereal area_tot__;
    extern /* Subroutine */ int connolly_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal probe;



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




/*     make call to the volume and surface area routine */

    /* Parameter adjustments */
    --radius;

    /* Function Body */
    probe = 0.;
    connolly_(volume_tot__, &area_tot__, &radius[1], &probe, exclude);
    return 0;
} /* volume_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine volume1  --  Cartesian excluded volume derivs  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "volume1" calculates first derivatives of the total excluded */
/*     volume with respect to the Cartesian coordinates of each atom */

/*     literature reference: */

/*     C. E. Kundrot, J. W. Ponder and F. M. Richards, "Algorithms for */
/*     Calculating Excluded Volume and Its Derivatives as a Function */
/*     of Molecular Conformation and Their Use in Energy Minimization", */
/*     Journal of Computational Chemistry, 12, 402-409 (1991) */


/* Subroutine */ int volume1_(doublereal *radius, doublereal *probe, 
	doublereal *dex)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 VOLUME1  --  Increase the Value of MAXCU"
	    "BE\002)";
    static char fmt_20[] = "(/,\002 VOLUME1  --  Increase \002,\002 the Valu"
	    "e of MAXARC\002)";
    static char fmt_30[] = "(/,\002 VOLUME1  --  Increase\002,\002 the Value"
	    " of MAXARC\002)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    double sqrt(doublereal), acos(doublereal), atan2(doublereal, doublereal), 
	    sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal phi_term__, d__[1000];
    static integer i__, j, k, m;
    static doublereal aa, ztopshave, bb;
    static integer in, io;
    static doublereal tf;
    static integer ir;
    static doublereal dx[1000], ti, dy[1000], rr;
    static integer nx, ny, nz;
    static doublereal xr, yr, zr, dsq[1000], phi1, phi2, pix2, rrx2, edge, 
	    beta, arcf[1000];
    static integer cube[6750]	/* was [2][15][15][15] */, itab[25000];
    static doublereal arci[1000];
    static integer narc;
    static doublereal temp;
    static logical skip[25000];
    static doublereal rmax;
    static integer inov[1000];
    static doublereal xmin;
    static integer isum;
    static doublereal ymin, zmin, xmax, ymax, zmax, rrsq, ztop, dist2, alpha;
    static integer icube;
    extern /* Subroutine */ int fatal_(void);
    static doublereal rdiff, rsecn;
    static integer itemp;
    static doublereal zgrid, rsecr, rinsq;
    static integer istop, jstop, kstop, mstop;
    static doublereal zstep, theta1, theta2, rsec2n, rsec2r, dtheta, seg_dx__,
	     seg_dy__, seg_dz__, pre_dx__, pre_dy__, pre_dz__, cosine;
    static integer istart, jstart, kstart, mstart;
    static doublereal vdwrad[25000], vdwsum, zstart, cos_phi1__, cos_phi2__;

    /* Fortran I/O blocks */
    static cilist io___19 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___72 = { 0, 0, 0, fmt_30, 0 };



#define dex_ref(a_1,a_2) dex[(a_2)*3 + a_1]
#define cube_ref(a_1,a_2,a_3,a_4) cube[(((a_4)*15 + (a_3))*15 + (a_2))*2 \
+ a_1 - 483]



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




/*     fix the stepsize in the z-direction; this value sets */
/*     the accuracy of the numerical derivatives; zstep=0.06 */
/*     is a good balance between compute time and accuracy */

    /* Parameter adjustments */
    dex -= 4;
    --radius;

    /* Function Body */
    zstep = .0601;

/*     initialize minimum and maximum ranges of atoms */

    pix2 = 6.2831853071795862;
    rmax = 0.;
    xmin = atoms_1.x[0];
    xmax = atoms_1.x[0];
    ymin = atoms_1.y[0];
    ymax = atoms_1.y[0];
    zmin = atoms_1.z__[0];
    zmax = atoms_1.z__[0];

/*     assign van der Waals radii to the atoms; note that */
/*     the radii are incremented by the size of the probe; */
/*     then get the maximum and minimum ranges of atoms */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vdwrad[i__ - 1] = radius[i__];
	if (vdwrad[i__ - 1] == 0.) {
	    skip[i__ - 1] = TRUE_;
	} else {
	    skip[i__ - 1] = FALSE_;
	    vdwrad[i__ - 1] += *probe;
	    if (vdwrad[i__ - 1] > rmax) {
		rmax = vdwrad[i__ - 1];
	    }
	    if (atoms_1.x[i__ - 1] < xmin) {
		xmin = atoms_1.x[i__ - 1];
	    }
	    if (atoms_1.x[i__ - 1] > xmax) {
		xmax = atoms_1.x[i__ - 1];
	    }
	    if (atoms_1.y[i__ - 1] < ymin) {
		ymin = atoms_1.y[i__ - 1];
	    }
	    if (atoms_1.y[i__ - 1] > ymax) {
		ymax = atoms_1.y[i__ - 1];
	    }
	    if (atoms_1.z__[i__ - 1] < zmin) {
		zmin = atoms_1.z__[i__ - 1];
	    }
	    if (atoms_1.z__[i__ - 1] > zmax) {
		zmax = atoms_1.z__[i__ - 1];
	    }
	}
    }

/*     load the cubes based on coarse lattice; first of all */
/*     set edge length to the maximum diameter of any atom */

    edge = rmax * 2.;
    nx = (integer) ((xmax - xmin) / edge) + 1;
    ny = (integer) ((ymax - ymin) / edge) + 1;
    nz = (integer) ((zmax - zmin) / edge) + 1;
/* Computing MAX */
    i__1 = max(nx,ny);
    if (max(i__1,nz) > 15) {
	io___19.ciunit = iounit_1.iout;
	s_wsfe(&io___19);
	e_wsfe();
	fatal_();
    }

/*     initialize the coarse lattice of cubes */

    i__1 = nx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ny;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = nz;
	    for (k = 1; k <= i__3; ++k) {
		cube_ref(1, i__, j, k) = 0;
		cube_ref(2, i__, j, k) = 0;
	    }
	}
    }

/*     find the number of atoms in each cube */

    i__1 = atoms_1.n;
    for (m = 1; m <= i__1; ++m) {
	if (! skip[m - 1]) {
	    i__ = (integer) ((atoms_1.x[m - 1] - xmin) / edge) + 1;
	    j = (integer) ((atoms_1.y[m - 1] - ymin) / edge) + 1;
	    k = (integer) ((atoms_1.z__[m - 1] - zmin) / edge) + 1;
	    cube_ref(1, i__, j, k) = cube_ref(1, i__, j, k) + 1;
	}
    }

/*     determine the highest index in the array "itab" for the */
/*     atoms that fall into each cube; the first cube that has */
/*     atoms defines the first index for "itab"; the final index */
/*     for the atoms in the present cube is the final index of */
/*     the last cube plus the number of atoms in the present cube */

    isum = 0;
    i__1 = nx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ny;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = nz;
	    for (k = 1; k <= i__3; ++k) {
		icube = cube_ref(1, i__, j, k);
		if (icube != 0) {
		    isum += icube;
		    cube_ref(2, i__, j, k) = isum;
		}
	    }
	}
    }

/*     "cube(2,,,)" now contains a pointer to the array "itab" */
/*     giving the position of the last entry for the list of */
/*     atoms in that cube of total number equal to "cube(1,,,)" */

    i__1 = atoms_1.n;
    for (m = 1; m <= i__1; ++m) {
	if (! skip[m - 1]) {
	    i__ = (integer) ((atoms_1.x[m - 1] - xmin) / edge) + 1;
	    j = (integer) ((atoms_1.y[m - 1] - ymin) / edge) + 1;
	    k = (integer) ((atoms_1.z__[m - 1] - zmin) / edge) + 1;
	    icube = cube_ref(2, i__, j, k);
	    itab[icube - 1] = m;
	    cube_ref(2, i__, j, k) = icube - 1;
	}
    }

/*     set "cube(2,,,)" to be the starting index in "itab" */
/*     for atom list of that cube; and "cube(1,,,)" to be */
/*     the stop index */

    isum = 0;
    i__1 = nx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ny;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = nz;
	    for (k = 1; k <= i__3; ++k) {
		icube = cube_ref(1, i__, j, k);
		if (icube != 0) {
		    isum += icube;
		    cube_ref(1, i__, j, k) = isum;
		    cube_ref(2, i__, j, k) = cube_ref(2, i__, j, k) + 1;
		}
	    }
	}
    }

/*     process in turn each atom from the coordinate list; */
/*     first select the potential intersecting atoms */

    i__1 = atoms_1.n;
    for (ir = 1; ir <= i__1; ++ir) {
	pre_dx__ = 0.;
	pre_dy__ = 0.;
	pre_dz__ = 0.;
	if (skip[ir - 1]) {
	    goto L50;
	}
	rr = vdwrad[ir - 1];
	rrx2 = rr * 2.;
	rrsq = rr * rr;
	xr = atoms_1.x[ir - 1];
	yr = atoms_1.y[ir - 1];
	zr = atoms_1.z__[ir - 1];

/*     find cubes to search for overlaps of current atom */

	istart = (integer) ((xr - xmin) / edge);
/* Computing MIN */
	i__2 = istart + 2;
	istop = min(i__2,nx);
	istart = max(istart,1);
	jstart = (integer) ((yr - ymin) / edge);
/* Computing MIN */
	i__2 = jstart + 2;
	jstop = min(i__2,ny);
	jstart = max(jstart,1);
	kstart = (integer) ((zr - zmin) / edge);
/* Computing MIN */
	i__2 = kstart + 2;
	kstop = min(i__2,nz);
	kstart = max(kstart,1);

/*     load all overlapping atoms into "inov" */

	io = 0;
	i__2 = istop;
	for (i__ = istart; i__ <= i__2; ++i__) {
	    i__3 = jstop;
	    for (j = jstart; j <= i__3; ++j) {
		i__4 = kstop;
		for (k = kstart; k <= i__4; ++k) {
		    mstart = cube_ref(2, i__, j, k);
		    if (mstart != 0) {
			mstop = cube_ref(1, i__, j, k);
			i__5 = mstop;
			for (m = mstart; m <= i__5; ++m) {
			    in = itab[m - 1];
			    if (in != ir) {
				++io;
				if (io > 1000) {
				    io___47.ciunit = iounit_1.iout;
				    s_wsfe(&io___47);
				    e_wsfe();
				    fatal_();
				}
				dx[io - 1] = atoms_1.x[in - 1] - xr;
				dy[io - 1] = atoms_1.y[in - 1] - yr;
/* Computing 2nd power */
				d__1 = dx[io - 1];
/* Computing 2nd power */
				d__2 = dy[io - 1];
				dsq[io - 1] = d__1 * d__1 + d__2 * d__2;
/* Computing 2nd power */
				d__1 = atoms_1.z__[in - 1] - zr;
				dist2 = dsq[io - 1] + d__1 * d__1;
/* Computing 2nd power */
				d__1 = rr + vdwrad[in - 1];
				vdwsum = d__1 * d__1;
				if (dist2 > vdwsum || dist2 == 0.) {
				    --io;
				} else {
				    d__[io - 1] = sqrt(dsq[io - 1]);
				    inov[io - 1] = in;
				}
			    }
			}
		    }
		}
	    }
	}

/*     determine resolution along the z-axis */

	if (io != 0) {
	    ztop = zr + rr;
	    ztopshave = ztop - zstep;
	    zgrid = zr - rr;

/*     half of the part not covered by the planes */

	    zgrid += (rrx2 - (integer) (rrx2 / zstep) * zstep) * .5;
	    zstart = zgrid;

/*     section atom spheres perpendicular to the z axis */

	    while(zgrid <= ztop) {

/*     "rsecr" is radius of circle of intersection */
/*     of "ir" sphere on the current sphere */

/* Computing 2nd power */
		d__1 = zgrid - zr;
		rsec2r = rrsq - d__1 * d__1;
		if (rsec2r < 0.) {
		    rsec2r = 1e-6;
		}
		rsecr = sqrt(rsec2r);
		if (zgrid >= ztopshave) {
		    cos_phi1__ = 1.;
		    phi1 = 0.;
		} else {
		    cos_phi1__ = (zgrid + zstep * .5 - zr) / rr;
		    phi1 = acos(cos_phi1__);
		}
		if (zgrid == zstart) {
		    cos_phi2__ = -1.;
		    phi2 = 3.141592653589793238;
		} else {
		    cos_phi2__ = (zgrid - zstep * .5 - zr) / rr;
		    phi2 = acos(cos_phi2__);
		}

/*     check intersections of neighbor circles */

		narc = 0;
		i__2 = io;
		for (k = 1; k <= i__2; ++k) {
		    in = inov[k - 1];
/* Computing 2nd power */
		    d__1 = vdwrad[in - 1];
		    rinsq = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = zgrid - atoms_1.z__[in - 1];
		    rsec2n = rinsq - d__1 * d__1;
		    if (rsec2n > 0.) {
			rsecn = sqrt(rsec2n);
			if (d__[k - 1] < rsecr + rsecn) {
			    rdiff = rsecr - rsecn;
			    if (d__[k - 1] <= abs(rdiff)) {
				if (rdiff < 0.) {
				    narc = 1;
				    arci[narc - 1] = 0.;
				    arcf[narc - 1] = pix2;
				}
				goto L40;
			    }
			    ++narc;
			    if (narc > 1000) {
				io___72.ciunit = iounit_1.iout;
				s_wsfe(&io___72);
				e_wsfe();
				fatal_();
			    }

/*     initial and final arc endpoints are found for intersection */
/*     of "ir" circle with another circle contained in same plane; */
/*     the initial endpoint of the enclosed arc is stored in "arci", */
/*     the final endpoint in "arcf"; get "cosine" via law of cosines */

			    cosine = (dsq[k - 1] + rsec2r - rsec2n) / (d__[k 
				    - 1] * 2. * rsecr);
/* Computing MIN */
			    d__1 = 1., d__2 = max(-1.,cosine);
			    cosine = min(d__1,d__2);

/*     "alpha" is the angle between a line containing either point */
/*     of intersection and the reference circle center and the */
/*     line containing both circle centers; "beta" is the angle */
/*     between the line containing both circle centers and x-axis */

			    alpha = acos(cosine);
			    beta = atan2(dy[k - 1], dx[k - 1]);
			    if (dy[k - 1] < 0.) {
				beta += pix2;
			    }
			    ti = beta - alpha;
			    tf = beta + alpha;
			    if (ti < 0.) {
				ti += pix2;
			    }
			    if (tf > pix2) {
				tf -= pix2;
			    }
			    arci[narc - 1] = ti;

/*     if the arc crosses zero, then it is broken into two segments; */
/*     the first ends at two pi and the second begins at zero */

			    if (tf < ti) {
				arcf[narc - 1] = pix2;
				++narc;
				arci[narc - 1] = 0.;
			    }
			    arcf[narc - 1] = tf;
L40:
			    ;
			}
		    }
		}

/*     find the pre-area and pre-forces on this section (band), */
/*     "pre-" means a multiplicative factor is yet to be applied */

		if (narc == 0) {
/* Computing 2nd power */
		    d__1 = cos_phi1__;
/* Computing 2nd power */
		    d__2 = cos_phi2__;
		    seg_dz__ = pix2 * (d__1 * d__1 - d__2 * d__2);
		    pre_dz__ += seg_dz__;
		} else {

/*     sort the arc endpoint arrays, each with "narc" entries, */
/*     in order of increasing values of the arguments in "arci" */

		    k = 1;
		    while(k < narc) {
			aa = arci[k - 1];
			bb = arcf[k - 1];
			temp = 1e6;
			i__2 = narc;
			for (i__ = k; i__ <= i__2; ++i__) {
			    if (arci[i__ - 1] <= temp) {
				temp = arci[i__ - 1];
				itemp = i__;
			    }
			}
			arci[k - 1] = arci[itemp - 1];
			arcf[k - 1] = arcf[itemp - 1];
			arci[itemp - 1] = aa;
			arcf[itemp - 1] = bb;
			++k;
		    }

/*     consolidate arcs by removing overlapping arc endpoints */

		    temp = arcf[0];
		    j = 1;
		    i__2 = narc;
		    for (k = 2; k <= i__2; ++k) {
			if (temp < arci[k - 1]) {
			    arcf[j - 1] = temp;
			    ++j;
			    arci[j - 1] = arci[k - 1];
			    temp = arcf[k - 1];
			} else if (temp < arcf[k - 1]) {
			    temp = arcf[k - 1];
			}
		    }
		    arcf[j - 1] = temp;
		    narc = j;
		    if (narc == 1) {
			narc = 2;
			arcf[1] = pix2;
			arci[1] = arcf[0];
			arcf[0] = arci[0];
			arci[0] = 0.;
		    } else {
			temp = arci[0];
			i__2 = narc - 1;
			for (k = 1; k <= i__2; ++k) {
			    arci[k - 1] = arcf[k - 1];
			    arcf[k - 1] = arci[k];
			}
			if (temp == 0. && arcf[narc - 1] == pix2) {
			    --narc;
			} else {
			    arci[narc - 1] = arcf[narc - 1];
			    arcf[narc - 1] = temp;
			}
		    }

/*     compute the numerical pre-derivative values */

		    i__2 = narc;
		    for (k = 1; k <= i__2; ++k) {
			theta1 = arci[k - 1];
			theta2 = arcf[k - 1];
			if (theta2 >= theta1) {
			    dtheta = theta2 - theta1;
			} else {
			    dtheta = theta2 + pix2 - theta1;
			}
			phi_term__ = phi2 - phi1 - (sin(phi2 * 2.) - sin(phi1 
				* 2.)) * .5;
			seg_dx__ = (sin(theta2) - sin(theta1)) * phi_term__;
			seg_dy__ = (cos(theta1) - cos(theta2)) * phi_term__;
/* Computing 2nd power */
			d__1 = cos_phi1__;
/* Computing 2nd power */
			d__2 = cos_phi2__;
			seg_dz__ = dtheta * (d__1 * d__1 - d__2 * d__2);
			pre_dx__ += seg_dx__;
			pre_dy__ += seg_dy__;
			pre_dz__ += seg_dz__;
		    }
		}
		zgrid += zstep;
	    }
	}
L50:
	dex_ref(1, ir) = rrsq * .5 * pre_dx__;
	dex_ref(2, ir) = rrsq * .5 * pre_dy__;
	dex_ref(3, ir) = rrsq * .5 * pre_dz__;
    }
    return 0;
} /* volume1_ */

#undef cube_ref
#undef dex_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine volume2  --  Cartesian excluded volume Hessian  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "volume2" calculates second derivatives of the total excluded */
/*     volume with respect to the Cartesian coordinates of the atoms */

/*     literature reference: */

/*     C. E. Kundrot, J. W. Ponder and F. M. Richards, "Algorithms for */
/*     Calculating Excluded Volume and Its Derivatives as a Function */
/*     of Molecular Conformation and Their Use in Energy Minimization", */
/*     Journal of Computational Chemistry, 12, 402-409 (1991) */


/* Subroutine */ int volume2_(integer *iatom, doublereal *radius, doublereal *
	probe, doublereal *xhess, doublereal *yhess, doublereal *zhess)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 VOLUME2  --  Increase\002,\002 the Value"
	    " of MAXARC\002)";
    static char fmt_20[] = "(/,\002 VOLUME2  -- Increase\002,\002 the Value "
	    "of MAXARC\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void);
    double acos(doublereal), atan(doublereal), cos(doublereal), sin(
	    doublereal);

    /* Local variables */
    static doublereal dfdtheta[6]	/* was [3][2] */;
    static integer arcfatom[1000], arciatom[1000];
    static doublereal dthetadx[18]	/* was [2][3][3] */, b, d__[1000];
    static integer i__, j, k, m;
    static doublereal r__[3], u[2], s2, ztopshave, aa, bb;
    static integer id[3], in;
    static doublereal tf, ri, dx[1000], ti, dy[1000], rr, xr, yr, zr;
    static integer iaa, ibb;
    static doublereal r_r__[3], r_s__[2], dsq[1000], phi1, phi2, cos1, cos2, 
	    r_s2__[2], sin1, sin2, pix2, rrx2, beta, arcf[1000], arci[1000];
    static integer narc;
    static doublereal delx[2], dely[2], delz[2], duds[2], dudr[2], dsdx[12]	
	    /* was [2][2][3] */, drdz[6]	/* was [2][3] */, dudx[18]	
	    /* was [2][3][3] */, temp, rrsq, ztop, dist2, rcut2, gamma, alpha;
    extern /* Subroutine */ int fatal_(void);
    static integer inear[1000], nnear;
    static doublereal rsecn, tempf, phi_z__;
    static integer itemp;
    static doublereal zgrid, rsecr, rinsq, zstep, theta1, theta2, rsec2n, 
	    rsec2r;
    static integer iblock, idtemp;
    static doublereal phiold, cosine, firsti, phi_xy__, u_term__[2], vdwrad[
	    25000], zstart, dbetadx[12]	/* was [2][2][3] */, dalphdx[18]	
	    /* was [2][3][3] */;
    static logical covered;
    static integer idfirst;

    /* Fortran I/O blocks */
    static cilist io___108 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___130 = { 0, 0, 0, fmt_20, 0 };



#define dfdtheta_ref(a_1,a_2) dfdtheta[(a_2)*3 + a_1 - 4]
#define dthetadx_ref(a_1,a_2,a_3) dthetadx[((a_3)*3 + (a_2))*2 + a_1 - 3]
#define dsdx_ref(a_1,a_2,a_3) dsdx[((a_3)*2 + (a_2))*2 + a_1 - 3]
#define drdz_ref(a_1,a_2) drdz[(a_2)*2 + a_1 - 1]
#define dudx_ref(a_1,a_2,a_3) dudx[((a_3)*3 + (a_2))*2 + a_1 - 3]
#define xhess_ref(a_1,a_2) xhess[(a_2)*3 + a_1]
#define yhess_ref(a_1,a_2) yhess[(a_2)*3 + a_1]
#define zhess_ref(a_1,a_2) zhess[(a_2)*3 + a_1]
#define dbetadx_ref(a_1,a_2,a_3) dbetadx[((a_3)*2 + (a_2))*2 + a_1 - 3]
#define dalphdx_ref(a_1,a_2,a_3) dalphdx[((a_3)*3 + (a_2))*2 + a_1 - 3]



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




/*     fix the stepsize in the z-direction; this value sets */
/*     the accuracy of the numerical derivatives; zstep=0.06 */
/*     is a good balance between compute time and accuracy */

    /* Parameter adjustments */
    zhess -= 4;
    yhess -= 4;
    xhess -= 4;
    --radius;

    /* Function Body */
    zstep = .0601;

/*     zero out the Hessian elements for current atom */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    xhess_ref(j, i__) = 0.;
	    yhess_ref(j, i__) = 0.;
	    zhess_ref(j, i__) = 0.;
	}
    }
    if (radius[*iatom] == 0.) {
	return 0;
    }
    pix2 = 6.2831853071795862;

/*     assign van der Waals radii to the atoms; note that */
/*     the radii are incremented by the size of the probe */

    i__1 = atoms_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vdwrad[i__ - 1] = radius[i__];
	if (vdwrad[i__ - 1] != 0.) {
	    vdwrad[i__ - 1] += *probe;
	}
    }

/*     set the radius and coordinates for current atom */

    rr = vdwrad[*iatom - 1];
    rrx2 = rr * 2.;
/* Computing 2nd power */
    d__1 = rr;
    rrsq = d__1 * d__1;
    xr = atoms_1.x[*iatom - 1];
    yr = atoms_1.y[*iatom - 1];
    zr = atoms_1.z__[*iatom - 1];

/*     select potential intersecting atoms */

    nnear = 1;
    i__1 = atoms_1.n;
    for (j = 1; j <= i__1; ++j) {
	if (j != *iatom && vdwrad[j - 1] != 0.) {
	    dx[nnear - 1] = atoms_1.x[j - 1] - xr;
	    dy[nnear - 1] = atoms_1.y[j - 1] - yr;
/* Computing 2nd power */
	    d__1 = dx[nnear - 1];
/* Computing 2nd power */
	    d__2 = dy[nnear - 1];
	    dsq[nnear - 1] = d__1 * d__1 + d__2 * d__2;
/* Computing 2nd power */
	    d__1 = atoms_1.z__[j - 1] - zr;
	    dist2 = dsq[nnear - 1] + d__1 * d__1;
/* Computing 2nd power */
	    d__1 = vdwrad[j - 1] + rr;
	    rcut2 = d__1 * d__1;
	    if (dist2 < rcut2) {
		d__[nnear - 1] = sqrt(dsq[nnear - 1]);
		inear[nnear - 1] = j;
		++nnear;
		if (nnear > 1000) {
		    io___108.ciunit = iounit_1.iout;
		    s_wsfe(&io___108);
		    e_wsfe();
		    fatal_();
		}
	    }
	}
    }
    --nnear;

/*     determine the z resolution */

    if (nnear != 0) {
	ztop = zr + rr;
	ztopshave = ztop - zstep;
	zgrid = zr - rr;

/*     half of the part not covered by the planes */

	zgrid += (rrx2 - (integer) (rrx2 / zstep) * zstep) * .5;
	zstart = zgrid;

/*     section atom spheres perpendicular to the z axis */

	while(zgrid <= ztop) {

/*     "rsecr" is radius of current atom sphere on the z-plane */

/* Computing 2nd power */
	    d__1 = zgrid - zr;
	    rsec2r = rrsq - d__1 * d__1;
	    if (rsec2r < 0.) {
		rsec2r = 1e-6;
	    }
	    rsecr = sqrt(rsec2r);
	    if (zgrid >= ztopshave) {
		phi1 = 0.;
	    } else {
		phi1 = acos((zgrid + zstep * .5 - zr) / rr);
	    }
	    if (zgrid == zstart) {
		phi2 = 3.141592653589793238;
	    } else {
		phi2 = phiold;
	    }

/*     check intersections of neighbor circles */

	    k = 0;
	    narc = 0;
	    covered = FALSE_;
	    while(! covered && k < nnear && narc < 1000) {
		++k;
		in = inear[k - 1];
/* Computing 2nd power */
		d__1 = vdwrad[in - 1];
		rinsq = d__1 * d__1;
/* Computing 2nd power */
		d__1 = zgrid - atoms_1.z__[in - 1];
		rsec2n = rinsq - d__1 * d__1;
		if (rsec2n > 0.) {
		    rsecn = sqrt(rsec2n);
		    if (d__[k - 1] < rsecr + rsecn) {
			b = rsecr - rsecn;
			if (d__[k - 1] <= abs(b)) {
			    if (b < 0.) {
				narc = 1;
				arci[narc - 1] = 0.;
				arcf[narc - 1] = pix2;
				arciatom[narc - 1] = in;
				arcfatom[narc - 1] = in;
				covered = TRUE_;
			    }
			} else {
			    ++narc;
			    if (narc > 1000) {
				io___130.ciunit = iounit_1.iout;
				s_wsfe(&io___130);
				e_wsfe();
				fatal_();
			    } else {

/*     initial and final arc endpoints are found for intersection */
/*     of "ir" circle with another circle contained in same plane; */
/*     the initial endpoint of the enclosed arc is stored in "arci", */
/*     the final endpoint in "arcf"; get "cosine" via law of cosines */

				cosine = (dsq[k - 1] + rsec2r - rsec2n) / (
					d__[k - 1] * 2. * rsecr);
/* Computing MIN */
				d__1 = 1., d__2 = max(-1.,cosine);
				cosine = min(d__1,d__2);

/*     "alpha" is the angle between a line containing either point */
/*     of intersection and the reference circle center and the */
/*     line containing both circle centers; "beta" is the angle */
/*     between the line containing both circle centers and x-axis */

				alpha = acos(cosine);
				if (dx[k - 1] == 0.) {
				    gamma = 1.5707963267948966;
				} else {
				    gamma = atan((d__1 = dy[k - 1] / dx[k - 1]
					    , abs(d__1)));
				}
				if (dy[k - 1] > 0.) {
				    if (dx[k - 1] > 0.) {
					beta = gamma;
				    } else {
					beta = 3.141592653589793238 - gamma;
				    }
				} else {
				    if (dx[k - 1] > 0.) {
					beta = pix2 - gamma;
				    } else {
					beta = gamma + 3.141592653589793238;
				    }
				}

/*     finally, the arc endpoints */

				ti = beta - alpha;
				tf = beta + alpha;
				if (ti < 0.) {
				    ti += pix2;
				}
				if (tf > pix2) {
				    tf -= pix2;
				}
				arci[narc - 1] = ti;
				arciatom[narc - 1] = in;
				arcfatom[narc - 1] = in;
				if (tf < ti) {
				    arcf[narc - 1] = pix2;
				    ++narc;
				    arci[narc - 1] = 0.;
				    arciatom[narc - 1] = in;
				    arcfatom[narc - 1] = in;
				}
				arcf[narc - 1] = tf;
			    }
			}
		    }
		}
	    }

/*     find the pre-area and pre-forces on this section (band) */
/*     through sphere "ir"; the "pre-" means a multiplicative */
/*     factor is yet to be applied */

	    if (narc != 0) {

/*     general case; sort arc endpoints */

		k = 1;
		while(k < narc) {
		    aa = arci[k - 1];
		    bb = arcf[k - 1];
		    iaa = arciatom[k - 1];
		    ibb = arcfatom[k - 1];
		    temp = 1e7;
		    i__1 = narc;
		    for (i__ = k; i__ <= i__1; ++i__) {
			if (arci[i__ - 1] <= temp) {
			    temp = arci[i__ - 1];
			    itemp = i__;
			}
		    }
		    arci[k - 1] = arci[itemp - 1];
		    arcf[k - 1] = arcf[itemp - 1];
		    arciatom[k - 1] = arciatom[itemp - 1];
		    arcfatom[k - 1] = arcfatom[itemp - 1];
		    arci[itemp - 1] = aa;
		    arcf[itemp - 1] = bb;
		    arciatom[itemp - 1] = iaa;
		    arcfatom[itemp - 1] = ibb;
		    ++k;
		}

/*     eliminate overlapping arc endpoints; */
/*     first, consolidate the occluded arcs */

		m = 1;
		tempf = arcf[0];
		idtemp = arcfatom[0];
		i__1 = narc;
		for (k = 2; k <= i__1; ++k) {
		    if (tempf < arci[k - 1]) {
			arcf[m - 1] = tempf;
			arcfatom[m - 1] = idtemp;
			++m;
			arci[m - 1] = arci[k - 1];
			arciatom[m - 1] = arciatom[k - 1];
			tempf = arcf[k - 1];
			idtemp = arcfatom[k - 1];
		    } else if (tempf < arcf[k - 1]) {
			tempf = arcf[k - 1];
			idtemp = arcfatom[k - 1];
		    }
		}
		arcf[m - 1] = tempf;
		arcfatom[m - 1] = idtemp;
		narc = m;

/*     change occluded arcs to accessible arcs */

		if (narc == 1) {
		    if (arci[0] == 0. && arcf[0] == pix2) {
			narc = 0;
		    } else {
			firsti = arci[0];
			idfirst = arciatom[0];
			arci[0] = arcf[0];
			arciatom[0] = arcfatom[0];
			arcf[0] = firsti + pix2;
			arcfatom[0] = idfirst;
		    }
		} else {
		    firsti = arci[0];
		    idfirst = arciatom[0];
		    i__1 = narc - 1;
		    for (k = 1; k <= i__1; ++k) {
			arci[k - 1] = arcf[k - 1];
			arciatom[k - 1] = arcfatom[k - 1];
			arcf[k - 1] = arci[k];
			arcfatom[k - 1] = arciatom[k];
		    }

/*     check gap between first and last arcs; if the */
/*     occluded arc crossed zero, then no accessible arc */

		    if (firsti == 0. && arcf[narc - 1] == pix2) {
			--narc;
		    } else {
			arci[narc - 1] = arcf[narc - 1];
			arciatom[narc - 1] = arcfatom[narc - 1];
			arcf[narc - 1] = firsti;
			arcfatom[narc - 1] = idfirst;
		    }
		}

/*     setup prior to application of chain rule */

		i__1 = narc;
		for (k = 1; k <= i__1; ++k) {
/* Computing 2nd power */
		    d__1 = zgrid - zr;
		    ri = sqrt(rrsq - d__1 * d__1);
		    for (i__ = 1; i__ <= 2; ++i__) {
			if (i__ == 1) {
			    id[1] = arciatom[k - 1];
			} else {
			    id[2] = arcfatom[k - 1];
			}
			delx[i__ - 1] = atoms_1.x[id[i__] - 1] - xr;
			dely[i__ - 1] = atoms_1.y[id[i__] - 1] - yr;
			delz[i__ - 1] = zgrid - atoms_1.z__[id[i__] - 1];
/* Computing 2nd power */
			d__1 = delx[i__ - 1];
/* Computing 2nd power */
			d__2 = dely[i__ - 1];
			s2 = d__1 * d__1 + d__2 * d__2;
			r_s__[i__ - 1] = 1. / sqrt(s2);
/* Computing 2nd power */
			d__1 = r_s__[i__ - 1];
			r_s2__[i__ - 1] = d__1 * d__1;
/* Computing 2nd power */
			d__1 = vdwrad[id[i__] - 1];
/* Computing 2nd power */
			d__2 = delz[i__ - 1];
			r__[i__] = sqrt(d__1 * d__1 - d__2 * d__2);
			r_r__[i__] = 1. / r__[i__];
/* Computing 2nd power */
			d__1 = ri;
/* Computing 2nd power */
			d__2 = r__[i__];
			u[i__ - 1] = (d__1 * d__1 + s2 - d__2 * d__2) * (
				r_s__[i__ - 1] * .5 / ri);
		    }

/*     apply the chain rule repeatedly */

		    theta1 = arci[k - 1];
		    theta2 = arcf[k - 1];
		    cos1 = cos(theta1);
		    cos2 = cos(theta2);
		    sin1 = sin(theta1);
		    sin2 = sin(theta2);
		    phi_xy__ = phi2 - phi1 - (sin(phi2 * 2.) - sin(phi1 * 2.))
			     * .5;
/* Computing 2nd power */
		    d__1 = sin(phi2);
/* Computing 2nd power */
		    d__2 = sin(phi1);
		    phi_z__ = d__1 * d__1 - d__2 * d__2;
		    phi_xy__ = rrsq * .5 * phi_xy__;
		    phi_z__ = rrsq * .5 * phi_z__;
		    dfdtheta_ref(1, 1) = -cos1 * phi_xy__;
		    dfdtheta_ref(2, 1) = -sin1 * phi_xy__;
		    dfdtheta_ref(3, 1) = -phi_z__;
		    dfdtheta_ref(1, 2) = cos2 * phi_xy__;
		    dfdtheta_ref(2, 2) = sin2 * phi_xy__;
		    dfdtheta_ref(3, 2) = phi_z__;
		    for (i__ = 1; i__ <= 2; ++i__) {
			dbetadx_ref(i__, 1, 0) = dely[i__ - 1] * r_s2__[i__ - 
				1];
			dbetadx_ref(i__, 2, 0) = -delx[i__ - 1] * r_s2__[i__ 
				- 1];
			dbetadx_ref(i__, 1, i__) = -dbetadx_ref(i__, 1, 0);
			dbetadx_ref(i__, 2, i__) = -dbetadx_ref(i__, 2, 0);
		    }
		    for (i__ = 1; i__ <= 2; ++i__) {
			duds[i__ - 1] = 1. / ri - u[i__ - 1] * r_s__[i__ - 1];
			dsdx_ref(i__, 1, i__) = delx[i__ - 1] * r_s__[i__ - 1]
				;
			dsdx_ref(i__, 2, i__) = dely[i__ - 1] * r_s__[i__ - 1]
				;
			dsdx_ref(i__, 1, 0) = -dsdx_ref(i__, 1, i__);
			dsdx_ref(i__, 2, 0) = -dsdx_ref(i__, 2, i__);
			dudr[i__ - 1] = -r__[i__] * r_s__[i__ - 1] / ri;
			drdz_ref(i__, i__) = delz[i__ - 1] * r_r__[i__];
			drdz_ref(i__, 0) = -drdz_ref(i__, i__);
		    }
		    for (m = 0; m <= 2; ++m) {
			for (i__ = 1; i__ <= 2; ++i__) {
			    dudx_ref(i__, 1, m) = duds[i__ - 1] * dsdx_ref(
				    i__, 1, m);
			    dudx_ref(i__, 2, m) = duds[i__ - 1] * dsdx_ref(
				    i__, 2, m);
			    dudx_ref(i__, 3, m) = dudr[i__ - 1] * drdz_ref(
				    i__, m);
			}
		    }
		    for (i__ = 1; i__ <= 2; ++i__) {
/* Computing 2nd power */
			d__1 = u[i__ - 1];
			u_term__[i__ - 1] = -1. / sqrt(1. - d__1 * d__1);
		    }
		    for (j = 1; j <= 3; ++j) {
			for (m = 0; m <= 2; ++m) {
			    for (i__ = 1; i__ <= 2; ++i__) {
				dalphdx_ref(i__, j, m) = u_term__[i__ - 1] * 
					dudx_ref(i__, j, m);
			    }
			}
		    }
		    for (j = 1; j <= 2; ++j) {
			for (m = 0; m <= 2; ++m) {
			    dthetadx_ref(1, j, m) = dbetadx_ref(1, j, m) + 
				    dalphdx_ref(1, j, m);
			    dthetadx_ref(2, j, m) = dbetadx_ref(2, j, m) - 
				    dalphdx_ref(2, j, m);
			}
		    }
		    for (m = 0; m <= 2; ++m) {
			dthetadx_ref(1, 3, m) = dalphdx_ref(1, 3, m);
			dthetadx_ref(2, 3, m) = -dalphdx_ref(2, 3, m);
		    }

/*     partials with respect to coordinates of serial atom id(m) */

		    id[0] = *iatom;
		    for (m = 0; m <= 2; ++m) {
			iblock = id[m];
			for (j = 1; j <= 3; ++j) {
			    xhess_ref(j, iblock) = xhess_ref(j, iblock) + 
				    dfdtheta_ref(1, 1) * dthetadx_ref(1, j, m)
				     + dfdtheta_ref(1, 2) * dthetadx_ref(2, j,
				     m);
			    yhess_ref(j, iblock) = yhess_ref(j, iblock) + 
				    dfdtheta_ref(2, 1) * dthetadx_ref(1, j, m)
				     + dfdtheta_ref(2, 2) * dthetadx_ref(2, j,
				     m);
			    zhess_ref(j, iblock) = zhess_ref(j, iblock) + 
				    dfdtheta_ref(3, 1) * dthetadx_ref(1, j, m)
				     + dfdtheta_ref(3, 2) * dthetadx_ref(2, j,
				     m);
			}
		    }
		}
	    }
	    zgrid += zstep;
	    phiold = phi1;
	}
    }
    return 0;
} /* volume2_ */

#undef dalphdx_ref
#undef dbetadx_ref
#undef zhess_ref
#undef yhess_ref
#undef xhess_ref
#undef dudx_ref
#undef drdz_ref
#undef dsdx_ref
#undef dthetadx_ref
#undef dfdtheta_ref


