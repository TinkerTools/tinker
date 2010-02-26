/* connolly.f -- translated by f2c (version 20050501).
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
    doublereal a[75000]	/* was [3][25000] */, ar[25000], pr;
    integer na;
} face01_;

#define face01_1 face01_

struct {
    logical skip[25000], nosurf[25000], afree[25000], abur[25000];
} face02_;

#define face02_1 face02_

struct {
    integer acls[50000]	/* was [2][25000] */, cls[1250000], clst[1250000];
} face03_;

#define face03_1 face03_

struct {
    integer ntt, tta[1250000]	/* was [2][625000] */, ttfe[625000], ttle[
	    625000], enext[125000];
    logical ttbur[625000], ttfree[625000];
} face04_;

#define face04_1 face04_

struct {
    doublereal t[225000]	/* was [3][75000] */, tr[75000], tax[225000]	
	    /* was [3][75000] */;
    integer nt, ta[150000]	/* was [2][75000] */, tfe[75000];
    logical tfree[75000];
} face05_;

#define face05_1 face05_

struct {
    doublereal p[150000]	/* was [3][50000] */;
    integer np, pa[150000]	/* was [3][50000] */;
} face06_;

#define face06_1 face06_

struct {
    doublereal v[375000]	/* was [3][125000] */;
    integer nv, va[125000], vp[125000];
} face07_;

#define face07_1 face07_

struct {
    integer nen, env[250000]	/* was [2][125000] */, nfn, fnen[150000]	
	    /* was [3][50000] */;
} face08_;

#define face08_1 face08_

struct {
    doublereal c__[375000]	/* was [3][125000] */, cr[125000];
    integer nc, ca[125000], ct[125000];
} face09_;

#define face09_1 face09_

struct {
    integer nep, epc[125000], epv[250000]	/* was [2][125000] */, afe[
	    25000], ale[25000], epnext[125000];
} face10_;

#define face10_1 face10_

struct {
    integer nfs, fsen[150000]	/* was [2][75000] */, fsep[150000]	/* 
	    was [2][75000] */;
} face11_;

#define face11_1 face11_

struct {
    integer ncy, cynep[25000], cyep[750000]	/* was [30][25000] */;
} face12_;

#define face12_1 face12_

struct {
    integer nfp, fpa[25000], fpcy[250000]	/* was [10][25000] */, fpncy[
	    25000];
} face13_;

#define face13_1 face13_

struct {
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

/* Table of constant values */

static doublereal c_b95 = 0.;
static integer c__1 = 1;
static doublereal c_b154 = -1.;
static doublereal c_b155 = 1.;



/*     ################################################################## */
/*     ##  COPYRIGHT (C) 1990 by Jay William Ponder                    ## */
/*     ##  COPYRIGHT (C) 1985 by Scripps Clinic & Research Foundation  ## */
/*     ##                     All Rights Reserved                      ## */
/*     ################################################################## */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine connolly  --  analytical surface area & volume  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "connolly" uses the algorithms from the AMS/VAM programs of */
/*     Michael Connolly to compute the analytical molecular surface */
/*     area and volume of a collection of spherical atoms; thus */
/*     it implements Fred Richards' molecular surface definition as */
/*     a set of analytically defined spherical and toroidal polygons */

/*     literature references: */

/*     M. L. Connolly, "Analytical Molecular Surface Calculation", */
/*     Journal of Applied Crystallography, 16, 548-558 (1983) */

/*     M. L. Connolly, "Computation of Molecular Volume", Journal */
/*     of the American Chemical Society, 107, 1118-1124 (1985) */

/*     variables only in the Connolly routines: */

/*     na      number of atoms */
/*     ntt     number of temporary tori */
/*     nt      number of tori */
/*     np      number of probe positions */
/*     nv      number of vertices */
/*     nen     number of concave edges */
/*     nfn     number of concave faces */
/*     nc      number of circles */
/*     nep     number of convex edges */
/*     nfs     number of saddle faces */
/*     ncy     number of cycles */
/*     fpncy   number of cycles bounding convex face */
/*     nfp     number of convex faces */
/*     cynep   number of convex edges in cycle */

/*     a       atomic coordinates */
/*     ar      atomic radii */
/*     pr      probe radius */
/*     skip    if true, atom is not used */
/*     nosurf  if true, atom has no free surface */
/*     afree   atom free of neighbors */
/*     abur    atom buried */

/*     acls    begin and end pointers for atoms neighbors */
/*     cls     atom numbers of neighbors */
/*     clst    pointer from neighbor to torus */

/*     tta     torus atom numbers */
/*     ttfe    first edge of each temporary torus */
/*     ttle    last edge of each temporary torus */
/*     enext   pointer to next edge of torus */
/*     ttbur   torus buried */
/*     ttfree  torus free */

/*     t       torus center */
/*     tr      torus radius */
/*     tax     torus axis */
/*     ta      torus atom numbers */
/*     tfe     torus first edge */
/*     tfree   torus free of neighbors */

/*     p       probe coordinates */
/*     pa      probe atom numbers */
/*     v       vertex coordinates */
/*     va      vertex atom number */
/*     vp      vertex probe number */
/*     c       circle center */
/*     cr      circle radius */
/*     ca      circle atom number */
/*     ct      circle torus number */

/*     env     concave edge vertex numbers */
/*     fnen    concave face concave edge numbers */
/*     epc     convex edge circle number */
/*     epv     convex edge vertex numbers */
/*     afe     first convex edge of each atom */
/*     ale     last convex edge of each atom */
/*     epnext  pointer to next convex edge of atom */
/*     fsen    saddle face concave edge numbers */
/*     fsep    saddle face convex edge numbers */
/*     cyep    cycle convex edge numbers */
/*     fpa     atom number of convex face */
/*     fpcy    convex face cycle numbers */


/* Subroutine */ int connolly_(doublereal *volume, doublereal *area, 
	doublereal *radius, doublereal *probe, doublereal *exclude)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int compress_(void);
    static integer i__;
    extern /* Subroutine */ int vam_(doublereal *, doublereal *), place_(void)
	    , torus_(void), nearby_(void), saddles_(void), contact_(void);


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]



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
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




/*     set the probe radius and the number of atoms */

    /* Parameter adjustments */
    --radius;

    /* Function Body */
    face01_1.pr = *probe;
    face01_1.na = atoms_1.n;

/*     set atom coordinates and radii, the excluded buffer */
/*     radius ("exclude") is added to atomic radii */

    i__1 = face01_1.na;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a_ref(1, i__) = atoms_1.x[i__ - 1];
	a_ref(2, i__) = atoms_1.y[i__ - 1];
	a_ref(3, i__) = atoms_1.z__[i__ - 1];
	face01_1.ar[i__ - 1] = radius[i__];
	if (face01_1.ar[i__ - 1] == 0.) {
	    face02_1.skip[i__ - 1] = TRUE_;
	} else {
	    face01_1.ar[i__ - 1] += *exclude;
	    face02_1.skip[i__ - 1] = FALSE_;
	}
    }

/*     find the analytical volume and surface area */

/*     call wiggle */
    nearby_();
    torus_();
    place_();
    compress_();
    saddles_();
    contact_();
    vam_(volume, area);
    return 0;
} /* connolly_ */

#undef a_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine wiggle  --  random perturbation of coordinates  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "wiggle" applies a random perturbation to the atomic coordinates */
/*     to avoid numerical instabilities for symmetrical structures */


/* Subroutine */ int wiggle_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal size;
    extern /* Subroutine */ int ranvec_(doublereal *);
    static doublereal vector[3];


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




/*     apply a small perturbation of fixed magnitude to each atom */

    size = 1e-6;
    i__1 = face01_1.na;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ranvec_(vector);
	a_ref(1, i__) = a_ref(1, i__) + size * vector[0];
	a_ref(2, i__) = a_ref(2, i__) + size * vector[1];
	a_ref(3, i__) = a_ref(3, i__) + size * vector[2];
    }
    return 0;
} /* wiggle_ */

#undef a_ref




/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  subroutine nearby  --  list of neighboring atom pairs  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "nearby" finds all of the through-space neighbors of each */
/*     atom for use in surface area and volume calculations */

/*     local variables : */

/*     ico      integer cube coordinates */
/*     icuptr   pointer to next atom in cube */
/*     comin    minimum atomic coordinates (cube corner) */
/*     icube    pointer to first atom in list for cube */
/*     scube    true if cube contains active atoms */
/*     sscube   true if cube or adjacent cubes have active atoms */
/*     itnl     temporary neighbor list, before sorting */


/* Subroutine */ int nearby_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal d2;
    static integer i1, j1, k1;
    static doublereal r2;
    static integer ici, icj, ick, jci, jcj, jck, ico[75000]	/* was [3][
	    25000] */;
    static doublereal sum;
    static integer clsa[1000], jcls, jmin, ncls, juse, itnl[1000];
    static doublereal sumi;
    static integer iptr;
    static doublereal vect1, vect2, vect3;
    extern doublereal dist2_(doublereal *, doublereal *);
    static integer icube[64000]	/* was [40][40][40] */, nclsa;
    static logical scube[64000]	/* was [40][40][40] */;
    static integer jmold;
    static doublereal comin[3];
    static integer iatom, jatom;
    static doublereal width, radmax;
    static logical sscube[64000]	/* was [40][40][40] */;
    extern /* Subroutine */ int cerror_(char *, ftnlen);
    static integer icuptr[25000], jmincls;


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]
#define ico_ref(a_1,a_2) ico[(a_2)*3 + a_1 - 4]
#define acls_ref(a_1,a_2) face03_1.acls[(a_2)*2 + a_1 - 3]
#define icube_ref(a_1,a_2,a_3) icube[((a_3)*40 + (a_2))*40 + a_1 - 1641]
#define scube_ref(a_1,a_2,a_3) scube[((a_3)*40 + (a_2))*40 + a_1 - 1641]
#define sscube_ref(a_1,a_2,a_3) sscube[((a_3)*40 + (a_2))*40 + a_1 - 1641]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




/*     ignore all atoms that are completely inside another atom; */
/*     may give nonsense results if this step is not taken */

    i__1 = face01_1.na - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (! face02_1.skip[i__ - 1]) {
	    i__2 = face01_1.na;
	    for (j = i__ + 1; j <= i__2; ++j) {
		d2 = dist2_(&a_ref(1, i__), &a_ref(1, j));
/* Computing 2nd power */
		d__1 = face01_1.ar[i__ - 1] - face01_1.ar[j - 1];
		r2 = d__1 * d__1;
		if (! face02_1.skip[j - 1] && d2 < r2) {
		    if (face01_1.ar[i__ - 1] < face01_1.ar[j - 1]) {
			face02_1.skip[i__ - 1] = TRUE_;
		    } else {
			face02_1.skip[j - 1] = TRUE_;
		    }
		}
	    }
	}
    }

/*     check for new coordinate minima and radii maxima */

    radmax = 0.;
    for (k = 1; k <= 3; ++k) {
	comin[k - 1] = a_ref(k, 1);
    }
    i__1 = face01_1.na;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (k = 1; k <= 3; ++k) {
	    if (a_ref(k, i__) < comin[k - 1]) {
		comin[k - 1] = a_ref(k, i__);
	    }
	}
	if (face01_1.ar[i__ - 1] > radmax) {
	    radmax = face01_1.ar[i__ - 1];
	}
    }

/*     calculate width of cube from maximum */
/*     atom radius and probe radius */

    width = (radmax + face01_1.pr) * 2.;

/*     set up cube arrays; first the integer coordinate arrays */

    i__1 = face01_1.na;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (k = 1; k <= 3; ++k) {
	    ico_ref(k, i__) = (integer) ((a_ref(k, i__) - comin[k - 1]) / 
		    width + 1);
	    if (ico_ref(k, i__) < 1) {
		cerror_("Cube Coordinate Too Small", (ftnlen)25);
	    } else if (ico_ref(k, i__) > 40) {
		cerror_("Cube Coordinate Too Large", (ftnlen)25);
	    }
	}
    }

/*     initialize head pointer and srn=2 arrays */

    for (i__ = 1; i__ <= 40; ++i__) {
	for (j = 1; j <= 40; ++j) {
	    for (k = 1; k <= 40; ++k) {
		icube_ref(i__, j, k) = 0;
		scube_ref(i__, j, k) = FALSE_;
		sscube_ref(i__, j, k) = FALSE_;
	    }
	}
    }

/*     initialize linked list pointers */

    i__1 = face01_1.na;
    for (i__ = 1; i__ <= i__1; ++i__) {
	icuptr[i__ - 1] = 0;
    }

/*     set up head and later pointers for each atom */

    i__1 = face01_1.na;
    for (iatom = 1; iatom <= i__1; ++iatom) {

/*     skip atoms with surface request numbers of zero */

	if (face02_1.skip[iatom - 1]) {
	    goto L30;
	}
	i__ = ico_ref(1, iatom);
	j = ico_ref(2, iatom);
	k = ico_ref(3, iatom);
	if (icube_ref(i__, j, k) <= 0) {

/*     first atom in this cube */

	    icube_ref(i__, j, k) = iatom;
	} else {

/*     add to end of linked list */

	    iptr = icube_ref(i__, j, k);
L10:

/*     check for duplicate atoms, turn off one of them */

	    if (dist2_(&a_ref(1, iatom), &a_ref(1, iptr)) <= 0.) {
		face02_1.skip[iatom - 1] = TRUE_;
		goto L30;
	    }

/*     move on down the list */

	    if (icuptr[iptr - 1] <= 0) {
		goto L20;
	    }
	    iptr = icuptr[iptr - 1];
	    goto L10;
L20:

/*     store atom number */

	    icuptr[iptr - 1] = iatom;
	}

/*     check for surfaced atom */

	if (! face02_1.skip[iatom - 1]) {
	    scube_ref(i__, j, k) = TRUE_;
	}
L30:
	;
    }

/*     check if this cube or any adjacent cube has active atoms */

    for (k = 1; k <= 40; ++k) {
	for (j = 1; j <= 40; ++j) {
	    for (i__ = 1; i__ <= 40; ++i__) {
		if (icube_ref(i__, j, k) != 0) {
/* Computing MAX */
		    i__1 = k - 1;
/* Computing MIN */
		    i__3 = k + 1;
		    i__2 = min(i__3,40);
		    for (k1 = max(i__1,1); k1 <= i__2; ++k1) {
/* Computing MAX */
			i__1 = j - 1;
/* Computing MIN */
			i__4 = j + 1;
			i__3 = min(i__4,40);
			for (j1 = max(i__1,1); j1 <= i__3; ++j1) {
/* Computing MAX */
			    i__1 = i__ - 1;
/* Computing MIN */
			    i__5 = i__ + 1;
			    i__4 = min(i__5,40);
			    for (i1 = max(i__1,1); i1 <= i__4; ++i1) {
				if (scube_ref(i1, j1, k1)) {
				    sscube_ref(i__, j, k) = TRUE_;
				}
			    }
			}
		    }
		}
	    }
	}
    }
    ncls = 0;

/*     zero pointers for atom and find its cube */

    i__2 = face01_1.na;
    for (i__ = 1; i__ <= i__2; ++i__) {
	nclsa = 0;
	face02_1.nosurf[i__ - 1] = face02_1.skip[i__ - 1];
	acls_ref(1, i__) = 0;
	acls_ref(2, i__) = 0;
	if (face02_1.skip[i__ - 1]) {
	    goto L70;
	}
	ici = ico_ref(1, i__);
	icj = ico_ref(2, i__);
	ick = ico_ref(3, i__);

/*     skip iatom if its cube and adjoining */
/*     cubes contain only blockers */

	if (! sscube_ref(ici, icj, ick)) {
	    goto L70;
	}
	sumi = face01_1.pr * 2. + face01_1.ar[i__ - 1];

/*     check iatom cube and adjacent cubes for neighboring atoms */

/* Computing MAX */
	i__3 = ick - 1;
/* Computing MIN */
	i__1 = ick + 1;
	i__4 = min(i__1,40);
	for (jck = max(i__3,1); jck <= i__4; ++jck) {
/* Computing MAX */
	    i__3 = icj - 1;
/* Computing MIN */
	    i__5 = icj + 1;
	    i__1 = min(i__5,40);
	    for (jcj = max(i__3,1); jcj <= i__1; ++jcj) {
/* Computing MAX */
		i__3 = ici - 1;
/* Computing MIN */
		i__6 = ici + 1;
		i__5 = min(i__6,40);
		for (jci = max(i__3,1); jci <= i__5; ++jci) {
		    j = icube_ref(jci, jcj, jck);
L40:

/*     check for end of linked list for this cube */

		    if (j <= 0) {
			goto L60;
		    }
		    if (i__ == j) {
			goto L50;
		    }
		    if (face02_1.skip[j - 1]) {
			goto L50;
		    }

/*     distance check */

		    sum = sumi + face01_1.ar[j - 1];
		    vect1 = (d__1 = a_ref(1, j) - a_ref(1, i__), abs(d__1));
		    if (vect1 >= sum) {
			goto L50;
		    }
		    vect2 = (d__1 = a_ref(2, j) - a_ref(2, i__), abs(d__1));
		    if (vect2 >= sum) {
			goto L50;
		    }
		    vect3 = (d__1 = a_ref(3, j) - a_ref(3, i__), abs(d__1));
		    if (vect3 >= sum) {
			goto L50;
		    }
/* Computing 2nd power */
		    d__1 = vect1;
/* Computing 2nd power */
		    d__2 = vect2;
/* Computing 2nd power */
		    d__3 = vect3;
		    d2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
/* Computing 2nd power */
		    d__1 = sum;
		    if (d2 >= d__1 * d__1) {
			goto L50;
		    }

/*     atoms are neighbors, save atom number in temporary array */

		    if (! face02_1.skip[j - 1]) {
			face02_1.nosurf[i__ - 1] = FALSE_;
		    }
		    ++nclsa;
		    if (nclsa > 1000) {
			cerror_("Too many Neighbors for Atom", (ftnlen)27);
		    }
		    itnl[nclsa - 1] = j;
L50:

/*     get number of next atom in cube */

		    j = icuptr[j - 1];
		    goto L40;
L60:
		    ;
		}
	    }
	}
	if (face02_1.nosurf[i__ - 1]) {
	    goto L70;
	}

/*     set up neighbors arrays with jatom in increasing order */

	jmold = 0;
	i__4 = nclsa;
	for (juse = 1; juse <= i__4; ++juse) {
	    jmin = face01_1.na + 1;
	    i__1 = nclsa;
	    for (jcls = 1; jcls <= i__1; ++jcls) {

/*     don't use ones already sorted */

		if (itnl[jcls - 1] > jmold) {
		    if (itnl[jcls - 1] < jmin) {
			jmin = itnl[jcls - 1];
			jmincls = jcls;
		    }
		}
	    }
	    jmold = jmin;
	    jcls = jmincls;
	    jatom = itnl[jcls - 1];
	    clsa[juse - 1] = jatom;
	}

/*     set up pointers to first and last neighbors of atom */

	if (nclsa > 0) {
	    acls_ref(1, i__) = ncls + 1;
	    i__4 = nclsa;
	    for (m = 1; m <= i__4; ++m) {
		++ncls;
		if (ncls > 1250000) {
		    cerror_("Too many Neighboring Atom Pairs", (ftnlen)31);
		}
		face03_1.cls[ncls - 1] = clsa[m - 1];
	    }
	    acls_ref(2, i__) = ncls;
	}
L70:
	;
    }
    return 0;
} /* nearby_ */

#undef sscube_ref
#undef scube_ref
#undef icube_ref
#undef acls_ref
#undef ico_ref
#undef a_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine torus  --  position of each temporary torus  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "torus" sets a list of all of the temporary torus positions */
/*     by testing for a torus between each atom and its neighbors */


/* Subroutine */ int torus_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer ia, ja, jn;
    static doublereal tt[3], ttr;
    static integer ibeg, iend;
    static doublereal ttax[3];
    static logical ttok;
    extern /* Subroutine */ int cerror_(char *, ftnlen), gettor_(integer *, 
	    integer *, logical *, doublereal *, doublereal *, doublereal *);


#define tta_ref(a_1,a_2) face04_1.tta[(a_2)*2 + a_1 - 3]
#define acls_ref(a_1,a_2) face03_1.acls[(a_2)*2 + a_1 - 3]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




/*     no torus is possible if there is only one atom */

    face04_1.ntt = 0;
    i__1 = face01_1.na;
    for (ia = 1; ia <= i__1; ++ia) {
	face02_1.afree[ia - 1] = TRUE_;
    }
    if (face01_1.na <= 1) {
	return 0;
    }

/*     get begin and end pointers to neighbors of this atom */

    i__1 = face01_1.na;
    for (ia = 1; ia <= i__1; ++ia) {
	if (! face02_1.nosurf[ia - 1]) {
	    ibeg = acls_ref(1, ia);
	    iend = acls_ref(2, ia);

/*     check for no neighbors */

	    if (ibeg > 0) {
		i__2 = iend;
		for (jn = ibeg; jn <= i__2; ++jn) {

/*     clear pointer from neighbor to torus */

		    face03_1.clst[jn - 1] = 0;

/*     get atom number of neighbor */

		    ja = face03_1.cls[jn - 1];

/*     don't create torus twice */

		    if (ja >= ia) {

/*     do some solid geometry */

			gettor_(&ia, &ja, &ttok, tt, &ttr, ttax);
			if (ttok) {

/*     we have a temporary torus, set up variables */

			    ++face04_1.ntt;
			    if (face04_1.ntt > 625000) {
				cerror_("Too many Temporary Tori", (ftnlen)23)
					;
			    }

/*     mark both atoms not free */

			    face02_1.afree[ia - 1] = FALSE_;
			    face02_1.afree[ja - 1] = FALSE_;
			    tta_ref(1, face04_1.ntt) = ia;
			    tta_ref(2, face04_1.ntt) = ja;

/*     pointer from neighbor to torus */

			    face03_1.clst[jn - 1] = face04_1.ntt;

/*     initialize torus as both free and buried */

			    face04_1.ttfree[face04_1.ntt - 1] = TRUE_;
			    face04_1.ttbur[face04_1.ntt - 1] = TRUE_;

/*     clear pointers from torus to first and last concave edges */

			    face04_1.ttfe[face04_1.ntt - 1] = 0;
			    face04_1.ttle[face04_1.ntt - 1] = 0;
			}
		    }
		}
	    }
	}
    }
    return 0;
} /* torus_ */

#undef acls_ref
#undef tta_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine place  --  locate positions of the probe sites  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "place" finds the probe sites by putting the probe sphere */
/*     tangent to each triple of neighboring atoms */


/* Subroutine */ int place_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer k, l;
    static doublereal d2;
    static integer l1, l2, ia, ja, ka, la, ke, ik, ip, jk, km, lm;
    static logical tb;
    static integer kv;
    static doublereal bij[3], hij;
    static integer lkf, mnb[500];
    static doublereal det;
    static integer ikt[500], jkt[500];
    static doublereal uij[3];
    static integer itt;
    static doublereal aijk[3];
    static integer iend, jend;
    static doublereal bijk[3], hijk;
    static integer nmnb;
    static doublereal pijk[3], uijk[3];
    static integer iptr, jptr;
    static logical ttok;
    extern doublereal dist2_(doublereal *, doublereal *);
    static integer lkcls[1250000];
    static logical prbok;
    static doublereal tempv[3];
    extern /* Subroutine */ int inedge_(integer *, integer *);
    static doublereal discls[500];
    extern /* Subroutine */ int cerror_(char *, ftnlen);
    static doublereal sumcls[500];
    extern /* Subroutine */ int gettor_(integer *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *), getprb_(integer *, 
	    integer *, integer *, logical *, logical *, doublereal *, 
	    doublereal *, doublereal *);


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]
#define p_ref(a_1,a_2) face06_1.p[(a_2)*3 + a_1 - 4]
#define v_ref(a_1,a_2) face07_1.v[(a_2)*3 + a_1 - 4]
#define pa_ref(a_1,a_2) face06_1.pa[(a_2)*3 + a_1 - 4]
#define tta_ref(a_1,a_2) face04_1.tta[(a_2)*2 + a_1 - 3]
#define env_ref(a_1,a_2) face08_1.env[(a_2)*2 + a_1 - 3]
#define acls_ref(a_1,a_2) face03_1.acls[(a_2)*2 + a_1 - 3]
#define fnen_ref(a_1,a_2) face08_1.fnen[(a_2)*3 + a_1 - 4]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




/*     no possible placement if there are no temporary tori */

    face06_1.np = 0;
    face08_1.nfn = 0;
    face08_1.nen = 0;
    face07_1.nv = 0;
    if (face04_1.ntt <= 0) {
	return 0;
    }

/*     consider each torus in turn */

    i__1 = face04_1.ntt;
    for (itt = 1; itt <= i__1; ++itt) {

/*     get atom numbers */

	ia = tta_ref(1, itt);
	ja = tta_ref(2, itt);

/*     form mutual neighbor list; clear number */
/*     of mutual neighbors of atoms ia and ja */

	nmnb = 0;

/*     get begin and end pointers for each atom's neighbor list */

	iptr = acls_ref(1, ia);
	jptr = acls_ref(1, ja);
	if (iptr <= 0 || jptr <= 0) {
	    goto L130;
	}
	iend = acls_ref(2, ia);
	jend = acls_ref(2, ja);

/*     collect mutual neighbors */

L10:

/*     check for end of loop */

	if (iptr > iend) {
	    goto L40;
	}
	if (jptr > jend) {
	    goto L40;
	}

/*     go move the lagging pointer */

	if (face03_1.cls[iptr - 1] < face03_1.cls[jptr - 1]) {
	    goto L20;
	}
	if (face03_1.cls[jptr - 1] < face03_1.cls[iptr - 1]) {
	    goto L30;
	}

/*     both point at same neighbor; one more mutual neighbor */
/*     save atom number of mutual neighbor */

	++nmnb;
	if (nmnb > 500) {
	    cerror_("Too many Mutual Neighbors", (ftnlen)25);
	}
	mnb[nmnb - 1] = face03_1.cls[iptr - 1];

/*     save pointers to second and third tori */

	ikt[nmnb - 1] = face03_1.clst[iptr - 1];
	jkt[nmnb - 1] = face03_1.clst[jptr - 1];
L20:

/*     increment pointer to ia atom neighbors */

	++iptr;
	goto L10;
L30:

/*     increment pointer to ja atom neighbors */

	++jptr;
	goto L10;
L40:

/*     we have all the mutual neighbors of ia and ja */
/*     if no mutual neighbors, skip to end of loop */

	if (nmnb <= 0) {
	    face04_1.ttbur[itt - 1] = FALSE_;
	    goto L130;
	}
	gettor_(&ia, &ja, &ttok, bij, &hij, uij);
	i__2 = nmnb;
	for (km = 1; km <= i__2; ++km) {
	    ka = mnb[km - 1];
	    discls[km - 1] = dist2_(bij, &a_ref(1, ka));
/* Computing 2nd power */
	    d__1 = face01_1.pr + face01_1.ar[ka - 1];
	    sumcls[km - 1] = d__1 * d__1;

/*     initialize link to next farthest out neighbor */

	    lkcls[km - 1] = 0;
	}

/*     set up a linked list of neighbors in order of */
/*     increasing distance from ia-ja torus center */

	lkf = 1;
	if (nmnb <= 1) {
	    goto L70;
	}

/*     put remaining neighbors in linked list at proper position */

	i__2 = nmnb;
	for (l = 2; l <= i__2; ++l) {
	    l1 = 0;
	    l2 = lkf;
L50:
	    if (discls[l - 1] < discls[l2 - 1]) {
		goto L60;
	    }
	    l1 = l2;
	    l2 = lkcls[l2 - 1];
	    if (l2 != 0) {
		goto L50;
	    }
L60:

/*     add to list */

	    if (l1 == 0) {
		lkf = l;
		lkcls[l - 1] = l2;
	    } else {
		lkcls[l1 - 1] = l;
		lkcls[l - 1] = l2;
	    }
	}
L70:

/*     loop thru mutual neighbors */

	i__2 = nmnb;
	for (km = 1; km <= i__2; ++km) {

/*     get atom number of neighbors */

	    ka = mnb[km - 1];
	    if (face02_1.skip[ia - 1] && face02_1.skip[ja - 1] && 
		    face02_1.skip[ka - 1]) {
		goto L120;
	    }

/*     get tori numbers for neighbor */

	    ik = ikt[km - 1];
	    jk = jkt[km - 1];

/*     possible new triple, do some geometry to */
/*     retrieve saddle center, axis and radius */

	    getprb_(&ia, &ja, &ka, &prbok, &tb, bijk, &hijk, uijk);
	    if (tb) {
		face04_1.ttbur[itt - 1] = TRUE_;
		face04_1.ttfree[itt - 1] = FALSE_;
		goto L120;
	    }

/*     no duplicate triples */

	    if (ka < ja) {
		goto L120;
	    }

/*     check whether any possible probe positions */

	    if (! prbok) {
		goto L120;
	    }

/*     altitude vector */

	    for (k = 1; k <= 3; ++k) {
		aijk[k - 1] = hijk * uijk[k - 1];
	    }

/*     we try two probe placements */

	    for (ip = 1; ip <= 2; ++ip) {
		for (k = 1; k <= 3; ++k) {
		    if (ip == 1) {
			pijk[k - 1] = bijk[k - 1] + aijk[k - 1];
		    } else {
			pijk[k - 1] = bijk[k - 1] - aijk[k - 1];
		    }
		}

/*     mark three tori not free */

		face04_1.ttfree[itt - 1] = FALSE_;
		face04_1.ttfree[ik - 1] = FALSE_;
		face04_1.ttfree[jk - 1] = FALSE_;

/*     check for collisions */

		lm = lkf;
L80:
		if (lm <= 0) {
		    goto L100;
		}

/*     get atom number of mutual neighbor */

		la = mnb[lm - 1];

/*     must not equal third atom */

		if (la == ka) {
		    goto L90;
		}

/*     compare distance to sum of radii */

		d2 = dist2_(pijk, &a_ref(1, la));
		if (d2 <= sumcls[lm - 1]) {
		    goto L110;
		}
L90:
		lm = lkcls[lm - 1];
		goto L80;
L100:

/*     we have a new probe position */

		++face06_1.np;
		if (face06_1.np > 50000) {
		    cerror_("Too many Probe Positions", (ftnlen)24);
		}

/*     mark three tori not buried */

		face04_1.ttbur[itt - 1] = FALSE_;
		face04_1.ttbur[ik - 1] = FALSE_;
		face04_1.ttbur[jk - 1] = FALSE_;

/*     store probe center */

		for (k = 1; k <= 3; ++k) {
		    p_ref(k, face06_1.np) = pijk[k - 1];
		}

/*     calculate vectors from probe to atom centers */

		if (face07_1.nv + 3 > 125000) {
		    cerror_("Too many Vertices", (ftnlen)17);
		}
		for (k = 1; k <= 3; ++k) {
		    v_ref(k, face07_1.nv + 1) = a_ref(k, ia) - p_ref(k, 
			    face06_1.np);
		    v_ref(k, face07_1.nv + 2) = a_ref(k, ja) - p_ref(k, 
			    face06_1.np);
		    v_ref(k, face07_1.nv + 3) = a_ref(k, ka) - p_ref(k, 
			    face06_1.np);
		}

/*     calculate determinant of vectors defining triangle */

		det = v_ref(1, face07_1.nv + 1) * v_ref(2, face07_1.nv + 2) * 
			v_ref(3, face07_1.nv + 3) + v_ref(1, face07_1.nv + 2) 
			* v_ref(2, face07_1.nv + 3) * v_ref(3, face07_1.nv + 
			1) + v_ref(1, face07_1.nv + 3) * v_ref(2, face07_1.nv 
			+ 1) * v_ref(3, face07_1.nv + 2) - v_ref(1, 
			face07_1.nv + 3) * v_ref(2, face07_1.nv + 2) * v_ref(
			3, face07_1.nv + 1) - v_ref(1, face07_1.nv + 2) * 
			v_ref(2, face07_1.nv + 1) * v_ref(3, face07_1.nv + 3) 
			- v_ref(1, face07_1.nv + 1) * v_ref(2, face07_1.nv + 
			3) * v_ref(3, face07_1.nv + 2);

/*     now add probe coordinates to vertices */

		for (k = 1; k <= 3; ++k) {
		    v_ref(k, face07_1.nv + 1) = p_ref(k, face06_1.np) + v_ref(
			    k, face07_1.nv + 1) * face01_1.pr / (face01_1.ar[
			    ia - 1] + face01_1.pr);
		    v_ref(k, face07_1.nv + 2) = p_ref(k, face06_1.np) + v_ref(
			    k, face07_1.nv + 2) * face01_1.pr / (face01_1.ar[
			    ja - 1] + face01_1.pr);
		    v_ref(k, face07_1.nv + 3) = p_ref(k, face06_1.np) + v_ref(
			    k, face07_1.nv + 3) * face01_1.pr / (face01_1.ar[
			    ka - 1] + face01_1.pr);
		}

/*     want the concave face to have counter-clockwise orientation */

		if (det > 0.) {

/*     swap second and third vertices */

		    for (k = 1; k <= 3; ++k) {
			tempv[k - 1] = v_ref(k, face07_1.nv + 2);
			v_ref(k, face07_1.nv + 2) = v_ref(k, face07_1.nv + 3);
			v_ref(k, face07_1.nv + 3) = tempv[k - 1];
		    }

/*     set up pointers from probe to atoms */

		    pa_ref(1, face06_1.np) = ia;
		    pa_ref(2, face06_1.np) = ka;
		    pa_ref(3, face06_1.np) = ja;

/*     set up pointers from vertices to atoms */

		    face07_1.va[face07_1.nv] = ia;
		    face07_1.va[face07_1.nv + 1] = ka;
		    face07_1.va[face07_1.nv + 2] = ja;

/*     insert concave edges into linked lists for appropriate tori */

		    i__3 = face08_1.nen + 1;
		    inedge_(&i__3, &ik);
		    i__3 = face08_1.nen + 2;
		    inedge_(&i__3, &jk);
		    i__3 = face08_1.nen + 3;
		    inedge_(&i__3, &itt);
		} else {

/*     similarly, if face already counter clockwise */

		    pa_ref(1, face06_1.np) = ia;
		    pa_ref(2, face06_1.np) = ja;
		    pa_ref(3, face06_1.np) = ka;
		    face07_1.va[face07_1.nv] = ia;
		    face07_1.va[face07_1.nv + 1] = ja;
		    face07_1.va[face07_1.nv + 2] = ka;
		    i__3 = face08_1.nen + 1;
		    inedge_(&i__3, &itt);
		    i__3 = face08_1.nen + 2;
		    inedge_(&i__3, &jk);
		    i__3 = face08_1.nen + 3;
		    inedge_(&i__3, &ik);
		}

/*     set up pointers from vertices to probe */

		for (kv = 1; kv <= 3; ++kv) {
		    face07_1.vp[face07_1.nv + kv - 1] = face06_1.np;
		}

/*     set up concave edges and concave face */

		if (face08_1.nen + 3 > 125000) {
		    cerror_("Too many Concave Edges", (ftnlen)22);
		}

/*     edges point to vertices */

		env_ref(1, face08_1.nen + 1) = face07_1.nv + 1;
		env_ref(2, face08_1.nen + 1) = face07_1.nv + 2;
		env_ref(1, face08_1.nen + 2) = face07_1.nv + 2;
		env_ref(2, face08_1.nen + 2) = face07_1.nv + 3;
		env_ref(1, face08_1.nen + 3) = face07_1.nv + 3;
		env_ref(2, face08_1.nen + 3) = face07_1.nv + 1;
		if (face08_1.nfn + 1 > 50000) {
		    cerror_("Too many Concave Faces", (ftnlen)22);
		}

/*     face points to edges */

		for (ke = 1; ke <= 3; ++ke) {
		    fnen_ref(ke, face08_1.nfn + 1) = face08_1.nen + ke;
		}

/*     increment counters for number of faces, edges and vertices */

		++face08_1.nfn;
		face08_1.nen += 3;
		face07_1.nv += 3;
L110:
		;
	    }
L120:
	    ;
	}
L130:
	;
    }
    return 0;
} /* place_ */

#undef fnen_ref
#undef acls_ref
#undef env_ref
#undef tta_ref
#undef pa_ref
#undef v_ref
#undef p_ref
#undef a_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine inedge  --  manage linked list of torus edges  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "inedge" inserts a concave edge into the */
/*     linked list for its temporary torus */


/* Subroutine */ int inedge_(integer *ien, integer *itt)
{
    static integer iepen;
    extern /* Subroutine */ int cerror_(char *, ftnlen);



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




/*     check for a serious error in the calling arguments */

    if (*ien <= 0) {
	cerror_("Bad Edge Number in INEDGE", (ftnlen)25);
    }
    if (*itt <= 0) {
	cerror_("Bad Torus Number in INEDGE", (ftnlen)26);
    }

/*     set beginning of list or add to end */

    if (face04_1.ttfe[*itt - 1] == 0) {
	face04_1.ttfe[*itt - 1] = *ien;
	face04_1.enext[*ien - 1] = 0;
	face04_1.ttle[*itt - 1] = *ien;
    } else {
	iepen = face04_1.ttle[*itt - 1];
	face04_1.enext[iepen - 1] = *ien;
	face04_1.enext[*ien - 1] = 0;
	face04_1.ttle[*itt - 1] = *ien;
    }
    return 0;
} /* inedge_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine compress  --  condense temporary to final tori  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "compress" transfers only the non-buried tori from */
/*     the temporary tori arrays to the final tori arrays */


/* Subroutine */ int compress_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer ia, ja, ip1, ip2, iv1, iv2, ned, itt, iptr;
    static logical ttok;
    extern /* Subroutine */ int cerror_(char *, ftnlen), gettor_(integer *, 
	    integer *, logical *, doublereal *, doublereal *, doublereal *);


#define t_ref(a_1,a_2) face05_1.t[(a_2)*3 + a_1 - 4]
#define ta_ref(a_1,a_2) face05_1.ta[(a_2)*2 + a_1 - 3]
#define tta_ref(a_1,a_2) face04_1.tta[(a_2)*2 + a_1 - 3]
#define env_ref(a_1,a_2) face08_1.env[(a_2)*2 + a_1 - 3]
#define tax_ref(a_1,a_2) face05_1.tax[(a_2)*3 + a_1 - 4]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




/*     initialize the number of nonburied tori */

    face05_1.nt = 0;
    if (face04_1.ntt <= 0) {
	return 0;
    }

/*     if torus is free, then it is not buried; */
/*     skip to end of loop if buried torus */

    i__1 = face04_1.ntt;
    for (itt = 1; itt <= i__1; ++itt) {
	if (face04_1.ttfree[itt - 1]) {
	    face04_1.ttbur[itt - 1] = FALSE_;
	}
	if (! face04_1.ttbur[itt - 1]) {

/*     first, transfer information */

	    ++face05_1.nt;
	    if (face05_1.nt > 75000) {
		cerror_("Too many NonBuried Tori", (ftnlen)23);
	    }
	    ia = tta_ref(1, itt);
	    ja = tta_ref(2, itt);
	    gettor_(&ia, &ja, &ttok, &t_ref(1, face05_1.nt), &face05_1.tr[
		    face05_1.nt - 1], &tax_ref(1, face05_1.nt));
	    ta_ref(1, face05_1.nt) = ia;
	    ta_ref(2, face05_1.nt) = ja;
	    face05_1.tfree[face05_1.nt - 1] = face04_1.ttfree[itt - 1];
	    face05_1.tfe[face05_1.nt - 1] = face04_1.ttfe[itt - 1];

/*     special check for inconsistent probes */

	    iptr = face05_1.tfe[face05_1.nt - 1];
	    ned = 0;
	    while(iptr != 0) {
		++ned;
		iptr = face04_1.enext[iptr - 1];
	    }
	    if (ned % 2 != 0) {
		iptr = face05_1.tfe[face05_1.nt - 1];
		while(iptr != 0) {
		    iv1 = env_ref(1, iptr);
		    iv2 = env_ref(2, iptr);
		    ip1 = face07_1.vp[iv1 - 1];
		    ip2 = face07_1.vp[iv2 - 1];
		    cerror_("Odd Torus for Probes IP1 and IP2", (ftnlen)32);
		    iptr = face04_1.enext[iptr - 1];
		}
	    }
	}
    }
    return 0;
} /* compress_ */

#undef tax_ref
#undef env_ref
#undef tta_ref
#undef ta_ref
#undef t_ref




/*     ############################################################## */
/*     ##                                                          ## */
/*     ##  subroutine saddles  --  builds saddle pieces from tori  ## */
/*     ##                                                          ## */
/*     ############################################################## */


/*     "saddles" constructs circles, convex edges and saddle faces */


/* Subroutine */ int saddles_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal);

    /* Local variables */
    static integer k, l1, l2, m1, n1, ia, in, ip;
    static doublereal dt;
    static integer it, iv, ien, ten[500];
    static doublereal tev[1500]	/* was [3][500] */;
    static integer ient;
    static doublereal dtev;
    static integer nent, itwo;
    static doublereal teang[500];
    extern /* Subroutine */ int ipedge_(integer *, integer *);
    static doublereal factor;
    static integer nxtang[500];
    extern doublereal triple_(doublereal *, doublereal *, doublereal *);
    static doublereal atvect[3];
    extern /* Subroutine */ int cerror_(char *, ftnlen);
    static logical sdstrt[500];


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]
#define c___ref(a_1,a_2) face09_1.c__[(a_2)*3 + a_1 - 4]
#define p_ref(a_1,a_2) face06_1.p[(a_2)*3 + a_1 - 4]
#define t_ref(a_1,a_2) face05_1.t[(a_2)*3 + a_1 - 4]
#define ta_ref(a_1,a_2) face05_1.ta[(a_2)*2 + a_1 - 3]
#define env_ref(a_1,a_2) face08_1.env[(a_2)*2 + a_1 - 3]
#define epv_ref(a_1,a_2) face10_1.epv[(a_2)*2 + a_1 - 3]
#define tax_ref(a_1,a_2) face05_1.tax[(a_2)*3 + a_1 - 4]
#define tev_ref(a_1,a_2) tev[(a_2)*3 + a_1 - 4]
#define fsen_ref(a_1,a_2) face11_1.fsen[(a_2)*2 + a_1 - 3]
#define fsep_ref(a_1,a_2) face11_1.fsep[(a_2)*2 + a_1 - 3]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




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




/*     zero the number of circles, convex edges and saddle faces */

    face09_1.nc = 0;
    face10_1.nep = 0;
    face11_1.nfs = 0;
    i__1 = face01_1.na;
    for (ia = 1; ia <= i__1; ++ia) {
	face10_1.afe[ia - 1] = 0;
	face10_1.ale[ia - 1] = 0;
	face02_1.abur[ia - 1] = TRUE_;
    }

/*     no saddle faces if no tori */

    if (face05_1.nt < 1) {
	return 0;
    }

/*     cycle through tori */

    i__1 = face05_1.nt;
    for (it = 1; it <= i__1; ++it) {
	if (face02_1.skip[ta_ref(1, it) - 1] && face02_1.skip[ta_ref(2, it) - 
		1]) {
	    goto L80;
	}

/*     set up two circles */

	for (in = 1; in <= 2; ++in) {
	    ia = ta_ref(in, it);

/*     mark atom not buried */

	    face02_1.abur[ia - 1] = FALSE_;

/*     vector from atom to torus center */

	    for (k = 1; k <= 3; ++k) {
		atvect[k - 1] = t_ref(k, it) - a_ref(k, ia);
	    }
	    factor = face01_1.ar[ia - 1] / (face01_1.ar[ia - 1] + face01_1.pr)
		    ;

/*     one more circle */

	    ++face09_1.nc;
	    if (face09_1.nc > 125000) {
		cerror_("Too many Circles", (ftnlen)16);
	    }

/*     circle center */

	    for (k = 1; k <= 3; ++k) {
		c___ref(k, face09_1.nc) = a_ref(k, ia) + factor * atvect[k - 
			1];
	    }

/*     pointer from circle to atom */

	    face09_1.ca[face09_1.nc - 1] = ia;

/*     pointer from circle to torus */

	    face09_1.ct[face09_1.nc - 1] = it;

/*     circle radius */

	    face09_1.cr[face09_1.nc - 1] = factor * face05_1.tr[it - 1];
	}

/*     skip to special code if free torus */

	if (face05_1.tfree[it - 1]) {
	    goto L70;
	}

/*     now we collect all the concave edges for this torus; */
/*     for each concave edge, calculate vector from torus center */
/*     thru probe center and the angle relative to first such vector */

/*     clear the number of concave edges for torus */

	nent = 0;

/*     pointer to start of linked list */

	ien = face05_1.tfe[it - 1];
L10:

/*     finished if concave edge pointer is zero */

	if (ien <= 0) {
	    goto L20;
	}

/*     one more concave edge */

	++nent;
	if (nent > 500) {
	    cerror_("Too many Edges for Torus", (ftnlen)24);
	}

/*     first vertex of edge */

	iv = env_ref(1, ien);

/*     probe number of vertex */

	ip = face07_1.vp[iv - 1];
	for (k = 1; k <= 3; ++k) {
	    tev_ref(k, nent) = p_ref(k, ip) - t_ref(k, it);
	}
	dtev = 0.;
	for (k = 1; k <= 3; ++k) {
/* Computing 2nd power */
	    d__1 = tev_ref(k, nent);
	    dtev += d__1 * d__1;
	}
	if (dtev <= 0.) {
	    cerror_("Probe on Torus Axis", (ftnlen)19);
	}
	dtev = sqrt(dtev);
	for (k = 1; k <= 3; ++k) {
	    tev_ref(k, nent) = tev_ref(k, nent) / dtev;
	}

/*     store concave edge number */

	ten[nent - 1] = ien;
	if (nent > 1) {

/*     calculate angle between this vector and first vector */

	    dt = 0.;
	    for (k = 1; k <= 3; ++k) {
		dt += tev_ref(k, 1) * tev_ref(k, nent);
	    }

/*     be careful */

	    if (dt > 1.) {
		dt = 1.;
	    }
	    if (dt < -1.) {
		dt = -1.;
	    }

/*     store angle */

	    teang[nent - 1] = acos(dt);

/*     get the sign right */

	    if (triple_(&tev_ref(1, 1), &tev_ref(1, nent), &tax_ref(1, it)) < 
		    0.) {
		teang[nent - 1] = 6.2831853071795862 - teang[nent - 1];
	    }
	} else {
	    teang[0] = 0.;
	}

/*     saddle face starts with this edge if it points parallel */
/*     to torus axis vector (which goes from first to second atom) */

	sdstrt[nent - 1] = face07_1.va[iv - 1] == ta_ref(1, it);

/*     next edge in list */

	ien = face04_1.enext[ien - 1];
	goto L10;
L20:
	if (nent <= 0) {
	    cerror_("No Edges for Non-free Torus", (ftnlen)27);
	}
	itwo = 2;
	if (nent % itwo != 0) {
	    cerror_("Odd Number of Edges for Torus", (ftnlen)29);
	}

/*     set up linked list of concave edges in order */
/*     of increasing angle around the torus axis; */
/*     clear second linked (angle-ordered) list pointers */

	i__2 = nent;
	for (ient = 1; ient <= i__2; ++ient) {
	    nxtang[ient - 1] = 0;
	}
	i__2 = nent;
	for (ient = 2; ient <= i__2; ++ient) {

/*     we have an entry to put into linked list */
/*     search for place to put it */

	    l1 = 0;
	    l2 = 1;
L30:
	    if (teang[ient - 1] < teang[l2 - 1]) {
		goto L40;
	    }

/*     not yet, move along */

	    l1 = l2;
	    l2 = nxtang[l2 - 1];
	    if (l2 != 0) {
		goto L30;
	    }
L40:

/*     we are at end of linked list or between l1 and l2; */
/*     insert edge */

	    if (l1 <= 0) {
		cerror_("Logic Error in SADDLES", (ftnlen)22);
	    }
	    nxtang[l1 - 1] = ient;
	    nxtang[ient - 1] = l2;
	}

/*     collect pairs of concave edges into saddles */
/*     create convex edges while you're at it */

	l1 = 1;
L50:
	if (l1 <= 0) {
	    goto L60;
	}

/*     check for start of saddle */

	if (sdstrt[l1 - 1]) {

/*     one more saddle face */

	    ++face11_1.nfs;
	    if (face11_1.nfs > 75000) {
		cerror_("Too many Saddle Faces", (ftnlen)21);
	    }

/*     get edge number */

	    ien = ten[l1 - 1];

/*     first concave edge of saddle */

	    fsen_ref(1, face11_1.nfs) = ien;

/*     one more convex edge */

	    ++face10_1.nep;
	    if (face10_1.nep > 125000) {
		cerror_("Too many Convex Edges", (ftnlen)21);
	    }

/*     first convex edge points to second circle */

	    face10_1.epc[face10_1.nep - 1] = face09_1.nc;

/*     atom circle lies on */

	    ia = face09_1.ca[face09_1.nc - 1];

/*     insert convex edge into linked list for atom */

	    ipedge_(&face10_1.nep, &ia);

/*     first vertex of convex edge is second vertex of concave edge */

	    epv_ref(1, face10_1.nep) = env_ref(2, ien);

/*     first convex edge of saddle */

	    fsep_ref(1, face11_1.nfs) = face10_1.nep;

/*     one more convex edge */

	    ++face10_1.nep;
	    if (face10_1.nep > 125000) {
		cerror_("Too many Convex Edges", (ftnlen)21);
	    }

/*     second convex edge points to first circle */

	    face10_1.epc[face10_1.nep - 1] = face09_1.nc - 1;
	    ia = face09_1.ca[face09_1.nc - 2];

/*     insert convex edge into linked list for atom */

	    ipedge_(&face10_1.nep, &ia);

/*     second vertex of second convex edge */
/*     is first vertex of first concave edge */

	    epv_ref(2, face10_1.nep) = env_ref(1, ien);
	    l1 = nxtang[l1 - 1];

/*     wrap around */

	    if (l1 <= 0) {
		l1 = 1;
	    }
	    if (sdstrt[l1 - 1]) {
		m1 = nxtang[l1 - 1];
		if (m1 <= 0) {
		    m1 = 1;
		}
		if (sdstrt[m1 - 1]) {
		    cerror_("Three Starts in a Row", (ftnlen)21);
		}
		n1 = nxtang[m1 - 1];

/*     the old switcheroo */

		nxtang[l1 - 1] = n1;
		nxtang[m1 - 1] = l1;
		l1 = m1;
	    }
	    ien = ten[l1 - 1];

/*     second concave edge for saddle face */

	    fsen_ref(2, face11_1.nfs) = ien;

/*     second vertex of first convex edge is */
/*     first vertex of second concave edge */

	    epv_ref(2, face10_1.nep - 1) = env_ref(1, ien);

/*     first vertex of second convex edge is */
/*     second vertex of second concave edge */

	    epv_ref(1, face10_1.nep) = env_ref(2, ien);
	    fsep_ref(2, face11_1.nfs) = face10_1.nep;

/*     quit if we have wrapped around to first edge */

	    if (l1 == 1) {
		goto L60;
	    }
	}

/*     next concave edge */

	l1 = nxtang[l1 - 1];
	goto L50;
L60:
	goto L80;

/*     free torus */

L70:

/*     set up entire circles as convex edges for new saddle surface; */
/*     one more saddle face */

	++face11_1.nfs;
	if (face11_1.nfs > 75000) {
	    cerror_("Too many Saddle Faces", (ftnlen)21);
	}

/*     no concave edges for saddle */

	fsen_ref(1, face11_1.nfs) = 0;
	fsen_ref(2, face11_1.nfs) = 0;

/*     one more convex edge */

	++face10_1.nep;
	ia = face09_1.ca[face09_1.nc - 1];

/*     insert convex edge into linked list for atom */

	ipedge_(&face10_1.nep, &ia);

/*     no vertices for convex edge */

	epv_ref(1, face10_1.nep) = 0;
	epv_ref(2, face10_1.nep) = 0;

/*     pointer from convex edge to second circle */

	face10_1.epc[face10_1.nep - 1] = face09_1.nc;

/*     first convex edge for saddle face */

	fsep_ref(1, face11_1.nfs) = face10_1.nep;

/*     one more convex edge */

	++face10_1.nep;
	ia = face09_1.ca[face09_1.nc - 2];

/*     insert second convex edge into linked list */

	ipedge_(&face10_1.nep, &ia);

/*     no vertices for convex edge */

	epv_ref(1, face10_1.nep) = 0;
	epv_ref(2, face10_1.nep) = 0;

/*     convex edge points to first circle */

	face10_1.epc[face10_1.nep - 1] = face09_1.nc - 1;

/*     second convex edge for saddle face */

	fsep_ref(2, face11_1.nfs) = face10_1.nep;

/*     buried torus; do nothing with it */

L80:
	;
    }
    return 0;
} /* saddles_ */

#undef fsep_ref
#undef fsen_ref
#undef tev_ref
#undef tax_ref
#undef epv_ref
#undef env_ref
#undef ta_ref
#undef t_ref
#undef p_ref
#undef c___ref
#undef a_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine gettor  --  test torus site between two atoms  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "gettor" tests for a possible torus position at the interface */
/*     between two atoms, and finds the torus radius, center and axis */


/* Subroutine */ int gettor_(integer *ia, integer *ja, logical *ttok, 
	doublereal *torcen, doublereal *torad, doublereal *torax)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k;
    static doublereal bij[3], dij, uij[3], vij[3], temp;
    extern doublereal dist2_(doublereal *, doublereal *);
    static doublereal temp1, temp2;


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




/*     get the distance between the two atoms */

    /* Parameter adjustments */
    --torax;
    --torcen;

    /* Function Body */
    *ttok = FALSE_;
    dij = sqrt(dist2_(&a_ref(1, *ia), &a_ref(1, *ja)));

/*     find a unit vector along interatomic (torus) axis */

    for (k = 1; k <= 3; ++k) {
	vij[k - 1] = a_ref(k, *ja) - a_ref(k, *ia);
	uij[k - 1] = vij[k - 1] / dij;
    }

/*     find coordinates of the center of the torus */

/* Computing 2nd power */
    d__1 = face01_1.ar[*ia - 1] + face01_1.pr;
/* Computing 2nd power */
    d__2 = face01_1.ar[*ja - 1] + face01_1.pr;
/* Computing 2nd power */
    d__3 = dij;
    temp = (d__1 * d__1 - d__2 * d__2) / (d__3 * d__3) + 1.;
    for (k = 1; k <= 3; ++k) {
	bij[k - 1] = a_ref(k, *ia) + vij[k - 1] * .5 * temp;
    }

/*     skip if atoms too far apart (should not happen) */

/* Computing 2nd power */
    d__1 = face01_1.ar[*ia - 1] + face01_1.ar[*ja - 1] + face01_1.pr * 2.;
/* Computing 2nd power */
    d__2 = dij;
    temp1 = d__1 * d__1 - d__2 * d__2;
    if (temp1 >= 0.) {

/*     skip if one atom is inside the other */

/* Computing 2nd power */
	d__1 = dij;
/* Computing 2nd power */
	d__2 = face01_1.ar[*ia - 1] - face01_1.ar[*ja - 1];
	temp2 = d__1 * d__1 - d__2 * d__2;
	if (temp2 >= 0.) {

/*     store the torus radius, center and axis */

	    *ttok = TRUE_;
	    *torad = sqrt(temp1 * temp2) / (dij * 2.);
	    for (k = 1; k <= 3; ++k) {
		torcen[k] = bij[k - 1];
		torax[k] = uij[k - 1];
	    }
	}
    }
    return 0;
} /* gettor_ */

#undef a_ref




/*     ################################################################## */
/*     ##                                                              ## */
/*     ##  subroutine getprb  --  test probe site between three atoms  ## */
/*     ##                                                              ## */
/*     ################################################################## */


/*     "getprb" tests for a possible probe position at the interface */
/*     between three neighboring atoms */


/* Subroutine */ int getprb_(integer *ia, integer *ja, integer *ka, logical *
	prbok, logical *tb, doublereal *bijk, doublereal *hijk, doublereal *
	uijk)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double acos(doublereal), sin(doublereal), sqrt(doublereal);

    /* Local variables */
    static integer k;
    static doublereal dba, rad, rij;
    extern doublereal dot_(doublereal *, doublereal *);
    static doublereal rik, uij[3], uik[3], utb[3], tij[3], tik[3];
    static logical tok;
    static doublereal rad2, dat2, rip2, fact, wijk;
    extern doublereal dist2_(doublereal *, doublereal *);
    static doublereal tijik[3], swijk, dotut, dotijk;
    extern /* Subroutine */ int gettor_(integer *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *), vcross_(doublereal *, 
	    doublereal *, doublereal *);


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




/*     initialize, then check torus over atoms "ia" and "ja" */

    /* Parameter adjustments */
    --uijk;
    --bijk;

    /* Function Body */
    *prbok = FALSE_;
    *tb = FALSE_;
    gettor_(ia, ja, &tok, tij, &rij, uij);
    if (! tok) {
	return 0;
    }
    dat2 = dist2_(&a_ref(1, *ka), tij);
/* Computing 2nd power */
    d__1 = face01_1.ar[*ka - 1] + face01_1.pr;
/* Computing 2nd power */
    d__2 = rij;
    rad2 = d__1 * d__1 - d__2 * d__2;

/*     if "ka" less than "ja", then all we care about */
/*     is whether the torus is buried */

    if (*ka < *ja) {
	if (rad2 <= 0.) {
	    return 0;
	}
	if (dat2 > rad2) {
	    return 0;
	}
    }
    gettor_(ia, ka, &tok, tik, &rik, uik);
    if (! tok) {
	return 0;
    }
    dotijk = dot_(uij, uik);
    if (dotijk > 1.) {
	dotijk = 1.;
    }
    if (dotijk < -1.) {
	dotijk = -1.;
    }
    wijk = acos(dotijk);
    swijk = sin(wijk);

/*     if the three atoms are colinear, then there is no */
/*     probe placement; but we still care whether the torus */
/*     is buried by atom "k" */

    if (swijk == 0.) {
	*tb = rad2 > 0. && dat2 <= rad2;
	return 0;
    }
    vcross_(uij, uik, &uijk[1]);
    for (k = 1; k <= 3; ++k) {
	uijk[k] /= swijk;
    }
    vcross_(&uijk[1], uij, utb);
    for (k = 1; k <= 3; ++k) {
	tijik[k - 1] = tik[k - 1] - tij[k - 1];
    }
    dotut = dot_(uik, tijik);
    fact = dotut / swijk;
    for (k = 1; k <= 3; ++k) {
	bijk[k] = tij[k - 1] + utb[k - 1] * fact;
    }
    dba = dist2_(&a_ref(1, *ia), &bijk[1]);
/* Computing 2nd power */
    d__1 = face01_1.ar[*ia - 1] + face01_1.pr;
    rip2 = d__1 * d__1;
    rad = rip2 - dba;
    if (rad < 0.) {
	*tb = rad2 > 0. && dat2 <= rad2;
    } else {
	*prbok = TRUE_;
	*hijk = sqrt(rad);
    }
    return 0;
} /* getprb_ */

#undef a_ref




/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine ipedge  --  manage linked list of convex edges  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "ipedge" inserts convex edge into linked list for atom */


/* Subroutine */ int ipedge_(integer *iep, integer *ia)
{
    static integer iepen;
    extern /* Subroutine */ int cerror_(char *, ftnlen);



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




/*     first, check for an error condition */

    if (*iep <= 0) {
	cerror_("Bad Edge Number in IPEDGE", (ftnlen)25);
    }
    if (*ia <= 0) {
	cerror_("Bad Atom Number in IPEDGE", (ftnlen)25);
    }

/*     set beginning of list or add to end */

    if (face10_1.afe[*ia - 1] == 0) {
	face10_1.afe[*ia - 1] = *iep;
	face10_1.epnext[*iep - 1] = 0;
	face10_1.ale[*ia - 1] = *iep;
    } else {
	iepen = face10_1.ale[*ia - 1];
	face10_1.epnext[iepen - 1] = *iep;
	face10_1.epnext[*iep - 1] = 0;
	face10_1.ale[*ia - 1] = *iep;
    }
    return 0;
} /* ipedge_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine contact  --  builds exposed contact surfaces  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "contact" constructs the contact surface, cycles and convex faces */


/* Subroutine */ int contact_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, k, ia, ic, jc, av[600]	/* was [2][300] */, it, ia2, 
	    aia[300], aic[300], aep[300];
    static doublereal acr[300];
    static integer iep, jcy;
    static doublereal anaa;
    static integer iepa, jepa, nepa, icya, jcya, kcya;
    static doublereal pole[3];
    static logical cycy[10000]	/* was [100][100] */, samef[10000]	/* 
	    was [100][100] */;
    static integer cyepa[3000]	/* was [30][100] */, ncypa, icyep, ncyep, 
	    jcyep, nused;
    extern doublereal anorm_(doublereal *);
    static integer lookv;
    static doublereal aavect[900]	/* was [3][300] */, acvect[900]	/* 
	    was [3][300] */, factor;
    static integer ncyepa[100], ncyold;
    static doublereal unvect[3];
    extern logical ptincy_(doublereal *, doublereal *, integer *);
    static logical epused[300], cyused[100];
    extern /* Subroutine */ int cerror_(char *, ftnlen);


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]
#define c___ref(a_1,a_2) face09_1.c__[(a_2)*3 + a_1 - 4]
#define ta_ref(a_1,a_2) face05_1.ta[(a_2)*2 + a_1 - 3]
#define av_ref(a_1,a_2) av[(a_2)*2 + a_1 - 3]
#define epv_ref(a_1,a_2) face10_1.epv[(a_2)*2 + a_1 - 3]
#define cyep_ref(a_1,a_2) face12_1.cyep[(a_2)*30 + a_1 - 31]
#define fpcy_ref(a_1,a_2) face13_1.fpcy[(a_2)*10 + a_1 - 11]
#define cycy_ref(a_1,a_2) cycy[(a_2)*100 + a_1 - 101]
#define samef_ref(a_1,a_2) samef[(a_2)*100 + a_1 - 101]
#define cyepa_ref(a_1,a_2) cyepa[(a_2)*30 + a_1 - 31]
#define aavect_ref(a_1,a_2) aavect[(a_2)*3 + a_1 - 4]
#define acvect_ref(a_1,a_2) acvect[(a_2)*3 + a_1 - 4]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




/*     zero out the number of cycles and convex faces */

    face12_1.ncy = 0;
    face13_1.nfp = 0;

/*     mark all free atoms not buried */

    i__1 = face01_1.na;
    for (ia = 1; ia <= i__1; ++ia) {
	if (face02_1.afree[ia - 1]) {
	    face02_1.abur[ia - 1] = FALSE_;
	}
    }

/*     go through all atoms */

    i__1 = face01_1.na;
    for (ia = 1; ia <= i__1; ++ia) {
	if (face02_1.skip[ia - 1]) {
	    goto L130;
	}

/*     skip to end of loop if buried atom */

	if (face02_1.abur[ia - 1]) {
	    goto L130;
	}

/*     special code for completely solvent-accessible atom */

	if (face02_1.afree[ia - 1]) {
	    goto L120;
	}

/*     gather convex edges for atom */
/*     clear number of convex edges for atom */

	nepa = 0;

/*     pointer to first edge */

	iep = face10_1.afe[ia - 1];
L10:

/*     check whether finished gathering */

	if (iep <= 0) {
	    goto L20;
	}

/*     one more edge */

	++nepa;
	if (nepa > 300) {
	    cerror_("Too many Convex Edges for Atom", (ftnlen)30);
	}

/*      store vertices of edge */

	av_ref(1, nepa) = epv_ref(1, iep);
	av_ref(2, nepa) = epv_ref(2, iep);

/*     store convex edge number */

	aep[nepa - 1] = iep;
	ic = face10_1.epc[iep - 1];

/*     store circle number */

	aic[nepa - 1] = ic;

/*     get neighboring atom */

	it = face09_1.ct[ic - 1];
	if (ta_ref(1, it) == ia) {
	    ia2 = ta_ref(2, it);
	} else {
	    ia2 = ta_ref(1, it);
	}

/*     store other atom number, we might need it sometime */

	aia[nepa - 1] = ia2;

/*     vector from atom to circle center; also */
/*     vector from atom to center of neighboring atom */
/*     sometimes we use one vector, sometimes the other */

	for (k = 1; k <= 3; ++k) {
	    acvect_ref(k, nepa) = c___ref(k, ic) - a_ref(k, ia);
	    aavect_ref(k, nepa) = a_ref(k, ia2) - a_ref(k, ia);
	}

/*     circle radius */

	acr[nepa - 1] = face09_1.cr[ic - 1];

/*     pointer to next edge */

	iep = face10_1.epnext[iep - 1];
	goto L10;
L20:
	if (nepa <= 0) {
	    cerror_("No Edges for Non-buried, Non-free Atom", (ftnlen)38);
	}


/*     form cycles; initialize all the */
/*     convex edges as not used in cycle */

	i__2 = nepa;
	for (iepa = 1; iepa <= i__2; ++iepa) {
	    epused[iepa - 1] = FALSE_;
	}

/*     save old number of cycles */

	ncyold = face12_1.ncy;
	nused = 0;
	ncypa = 0;
L30:

/*     look for starting edge */

	i__2 = nepa;
	for (iepa = 1; iepa <= i__2; ++iepa) {
	    if (! epused[iepa - 1]) {
		goto L40;
	    }
	}

/*     cannot find starting edge, finished */

	goto L80;
L40:

/*     pointer to edge */

	iep = aep[iepa - 1];

/*     one edge so far for this cycle */

	ncyep = 1;

/*     one more cycle for atom */

	++ncypa;
	if (ncypa > 100) {
	    cerror_("Too many Cycles per Atom", (ftnlen)24);
	}

/*     mark edge used in cycle */

	epused[iepa - 1] = TRUE_;
	++nused;

/*     one more cycle for molecule */

	++face12_1.ncy;
	if (face12_1.ncy > 25000) {
	    cerror_("Too many Cycles", (ftnlen)15);
	}

/*     index of edge in atom cycle array */

	cyepa_ref(ncyep, ncypa) = iepa;

/*     store in molecule cycle array a pointer to edge */

	cyep_ref(ncyep, face12_1.ncy) = iep;

/*     second vertex of this edge is the vertex to look */
/*     for next as the first vertex of another edge */

	lookv = av_ref(2, iepa);

/*     if no vertex, this cycle is finished */

	if (lookv <= 0) {
	    goto L70;
	}
L50:

/*     look for next connected edge */

	i__2 = nepa;
	for (jepa = 1; jepa <= i__2; ++jepa) {
	    if (epused[jepa - 1]) {
		goto L60;
	    }

/*     check second vertex of iepa versus first vertex of jepa */

	    if (av_ref(1, jepa) != lookv) {
		goto L60;
	    }

/*     edges are connected */
/*     pointer to edge */

	    iep = aep[jepa - 1];

/*     one more edge for this cycle */

	    ++ncyep;
	    if (ncyep > 30) {
		cerror_("Too many Edges per Cycle", (ftnlen)24);
	    }
	    epused[jepa - 1] = TRUE_;
	    ++nused;

/*     store index in local edge array */

	    cyepa_ref(ncyep, ncypa) = jepa;

/*     store pointer to edge */

	    cyep_ref(ncyep, face12_1.ncy) = iep;

/*     new vertex to look for */

	    lookv = av_ref(2, jepa);

/*     if no vertex, this cycle is in trouble */

	    if (lookv <= 0) {
		cerror_("Pointer Error in Cycle", (ftnlen)22);
	    }
	    goto L50;
L60:
	    ;
	}

/*     it better connect to first edge of cycle */

	if (lookv != av_ref(1, iepa)) {
	    cerror_("Cycle does not Close", (ftnlen)20);
	}
L70:

/*     this cycle is finished */
/*     store number of edges in cycle */

	ncyepa[ncypa - 1] = ncyep;
	face12_1.cynep[face12_1.ncy - 1] = ncyep;
	if (nused >= nepa) {
	    goto L80;
	}

/*     look for more cycles */

	goto L30;
L80:

/*     compare cycles for inside/outside relation; */
/*     check to see if cycle i is inside cycle j */

	i__2 = ncypa;
	for (icya = 1; icya <= i__2; ++icya) {
	    i__3 = ncypa;
	    for (jcya = 1; jcya <= i__3; ++jcya) {
		jcy = ncyold + jcya;

/*     initialize */

		cycy_ref(icya, jcya) = TRUE_;

/*     check for same cycle */

		if (icya == jcya) {
		    goto L90;
		}

/*     if cycle j has two or fewer edges, nothing can */
/*     lie in its exterior; i is therefore inside j */

		if (ncyepa[jcya - 1] <= 2) {
		    goto L90;
		}

/*     if cycles i and j have a pair of edges belonging */
/*     to the same circle, then they are outside each other */

		i__4 = ncyepa[icya - 1];
		for (icyep = 1; icyep <= i__4; ++icyep) {
		    iepa = cyepa_ref(icyep, icya);
		    ic = aic[iepa - 1];
		    i__5 = ncyepa[jcya - 1];
		    for (jcyep = 1; jcyep <= i__5; ++jcyep) {
			jepa = cyepa_ref(jcyep, jcya);
			jc = aic[jepa - 1];
			if (ic == jc) {
			    cycy_ref(icya, jcya) = FALSE_;
			    goto L90;
			}
		    }
		}
		iepa = cyepa_ref(1, icya);
		anaa = anorm_(&aavect_ref(1, iepa));
		factor = face01_1.ar[ia - 1] / anaa;

/*     north pole and unit vector pointing south */

		for (k = 1; k <= 3; ++k) {
		    pole[k - 1] = factor * aavect_ref(k, iepa) + a_ref(k, ia);
		    unvect[k - 1] = -aavect_ref(k, iepa) / anaa;
		}
		cycy_ref(icya, jcya) = ptincy_(pole, unvect, &jcy);
L90:
		;
	    }
	}

/*     group cycles into faces; direct comparison for i and j */

	i__2 = ncypa;
	for (icya = 1; icya <= i__2; ++icya) {
	    i__3 = ncypa;
	    for (jcya = 1; jcya <= i__3; ++jcya) {

/*     tentatively say that cycles i and j bound */
/*     the same face if they are inside each other */

		samef_ref(icya, jcya) = cycy_ref(icya, jcya) && cycy_ref(jcya,
			 icya);
	    }
	}

/*     if i is in exterior of k, and k is in interior of */
/*     i and j, then i and j do not bound the same face */

	i__2 = ncypa;
	for (icya = 1; icya <= i__2; ++icya) {
	    i__3 = ncypa;
	    for (jcya = 1; jcya <= i__3; ++jcya) {
		if (icya != jcya) {
		    i__4 = ncypa;
		    for (kcya = 1; kcya <= i__4; ++kcya) {
			if (kcya != icya && kcya != jcya) {
			    if (cycy_ref(kcya, icya) && cycy_ref(kcya, jcya) 
				    && ! cycy_ref(icya, kcya)) {
				samef_ref(icya, jcya) = FALSE_;
				samef_ref(jcya, icya) = FALSE_;
			    }
			}
		    }
		}
	    }
	}

/*     fill gaps so that "samef" falls into complete blocks */

	i__2 = ncypa - 2;
	for (icya = 1; icya <= i__2; ++icya) {
	    i__3 = ncypa - 1;
	    for (jcya = icya + 1; jcya <= i__3; ++jcya) {
		if (samef_ref(icya, jcya)) {
		    i__4 = ncypa;
		    for (kcya = jcya + 1; kcya <= i__4; ++kcya) {
			if (samef_ref(jcya, kcya)) {
			    samef_ref(icya, kcya) = TRUE_;
			    samef_ref(kcya, icya) = TRUE_;
			}
		    }
		}
	    }
	}

/*     group cycles belonging to the same face */

	i__2 = ncypa;
	for (icya = 1; icya <= i__2; ++icya) {
	    cyused[icya - 1] = FALSE_;
	}

/*     clear number of cycles used in bounding faces */

	nused = 0;
	i__2 = ncypa;
	for (icya = 1; icya <= i__2; ++icya) {

/*     check for already used */

	    if (cyused[icya - 1]) {
		goto L110;
	    }

/*     one more convex face */

	    ++face13_1.nfp;
	    if (face13_1.nfp > 25000) {
		cerror_("Too many Convex Faces", (ftnlen)21);
	    }

/*     clear number of cycles for face */

	    face13_1.fpncy[face13_1.nfp - 1] = 0;

/*     pointer from face to atom */

	    face13_1.fpa[face13_1.nfp - 1] = ia;

/*     look for all other cycles belonging to same face */

	    i__3 = ncypa;
	    for (jcya = 1; jcya <= i__3; ++jcya) {

/*     check for cycle already used in another face */

		if (cyused[jcya - 1]) {
		    goto L100;
		}

/*     cycles i and j belonging to same face */

		if (! samef_ref(icya, jcya)) {
		    goto L100;
		}

/*     mark cycle used */

		cyused[jcya - 1] = TRUE_;
		++nused;

/*     one more cycle for face */

		++face13_1.fpncy[face13_1.nfp - 1];
		if (face13_1.fpncy[face13_1.nfp - 1] > 10) {
		    cerror_("Too many Cycles bounding Convex Face", (ftnlen)
			    36);
		}
		i__ = face13_1.fpncy[face13_1.nfp - 1];

/*     store cycle number */

		fpcy_ref(i__, face13_1.nfp) = ncyold + jcya;

/*     check for finished */

		if (nused >= ncypa) {
		    goto L130;
		}
L100:
		;
	    }
L110:
	    ;
	}

/*     should not fall through end of do loops */

	cerror_("Not all Cycles grouped into Convex Faces", (ftnlen)40);
L120:

/*     one face for free atom; no cycles */

	++face13_1.nfp;
	if (face13_1.nfp > 25000) {
	    cerror_("Too many Convex Faces", (ftnlen)21);
	}
	face13_1.fpa[face13_1.nfp - 1] = ia;
	face13_1.fpncy[face13_1.nfp - 1] = 0;
L130:
	;
    }
    return 0;
} /* contact_ */

#undef acvect_ref
#undef aavect_ref
#undef cyepa_ref
#undef samef_ref
#undef cycy_ref
#undef fpcy_ref
#undef cyep_ref
#undef epv_ref
#undef av_ref
#undef ta_ref
#undef c___ref
#undef a_ref




/*     ########################################################## */
/*     ##                                                      ## */
/*     ##  subroutine vam  --  volumes and areas of molecules  ## */
/*     ##                                                      ## */
/*     ########################################################## */


/*     "vam" takes the analytical molecular surface defined */
/*     as a collection of spherical and toroidal polygons */
/*     and uses it to compute the volume and surface area */


/* Subroutine */ int vam_(doublereal *volume, doublereal *area)
{
    /* Format strings */
    static char fmt_150[] = "(/,\002 Convex Surface Area for Individual Atom"
	    "s :\002,/)";
    static char fmt_160[] = "(1x,5(i7,f8.3))";
    static char fmt_170[] = "(/,\002 Surface Area and Volume by Geometry Typ"
	    "e :\002)";
    static char fmt_180[] = "(/,\002 Convex Faces :\002,i12,5x,\002Area :"
	    "\002,f13.3,4x,\002Volume :\002,f13.3)";
    static char fmt_190[] = "(\002 Saddle Faces :\002,i12,5x,\002Area :\002,"
	    "f13.3,4x,\002Volume :\002,f13.3)";
    static char fmt_200[] = "(\002 Concave Faces :\002,i11,5x,\002Area :\002"
	    ",f13.3,4x,\002Volume :\002,f13.3)";
    static char fmt_210[] = "(\002 Buried Polyhedra :\002,36x,\002Volume "
	    ":\002,f13.3)";
    static char fmt_220[] = "(/,\002 Spindle Correction :\002,11x,\002Area "
	    ":\002,f13.3,4x,\002Volume :\002,f13.3)";
    static char fmt_230[] = "(\002 Lens Analytical Correction :\002,3x,\002A"
	    "rea :\002,f13.3,4x,\002Volume :\002,f13.3)";
    static char fmt_240[] = "(\002 Lens Numerical Correction :\002,4x,\002Ar"
	    "ea :\002,f13.3,4x,\002Volume :\002,f13.3)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), asin(doublereal), acos(doublereal), sin(
	    doublereal), cos(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer k, ia, ic, ke, ip;
    static doublereal dt, uc[3];
    static integer it, iv;
    static doublereal rm;
    static integer kv;
    static doublereal uq[3];
    static integer ke2;
    static doublereal ds2;
    static integer iv1, iv2;
    static logical ate[100];
    static integer ien, ifn, iep, ifp, ifs, isc, jfn, nop, iop, ivs[3], fnt[
	    150000]	/* was [3][50000] */;
    extern doublereal dot_(doublereal *, doublereal *);
    static doublereal dpp, rat, rsc, rho, stq, tau[3], ppm[3], qij[3], qji[3],
	     upp[3], umq[3], upq[3], uij[3];
    static logical vip[3];
    static doublereal dij2;
    static logical badt[50000], alli, allj;
    static doublereal cora[50000];
    static integer nate, neat, nlap[50000], enfs[125000], idot;
    static doublereal dota, alts[150000]	/* was [3][50000] */, corv[
	    50000], sdot[3], dotv[20], voln, dots[3000]	/* was [3][1000] */, 
	    vint, volp, vols;
    static integer nspt[150000]	/* was [3][50000] */;
    static logical anyi, anyj, case1, case2;
    static doublereal vpyr, vect1[3], vect2[3], vect3[3];
    extern doublereal dist2_(doublereal *, doublereal *);
    static doublereal vect4[3], vect5[3], vect6[3], vect7[3], vect8[3], xpnt1[
	    3], xpnt2[3];
    static logical badav[50000];
    static doublereal arean, areap, fncen[150000]	/* was [3][50000] */, 
	    areas, scinc, alens;
    extern doublereal depth_(integer *, doublereal *);
    static doublereal coran;
    static integer ifnop[100];
    static doublereal vcone, totan, voldo;
    static integer ndots;
    static doublereal vlens, totap, totas, prism, sumsc, corvn, cenop[300]	
	    /* was [3][100] */, vects[9]	/* was [3][3] */, volsp, 
	    tdots[3000]	/* was [3][1000] */;
    static logical fcins[150000]	/* was [3][50000] */, fcint[150000]	
	    /* was [3][50000] */, cinsp, cintp;
    extern /* Subroutine */ int vnorm_(doublereal *, doublereal *);
    static doublereal totvn, totvp, totvs;
    static integer ispnd2[3];
    static doublereal areado;
    extern /* Subroutine */ int measfn_(integer *, doublereal *, doublereal *)
	    ;
    static doublereal areasp;
    extern /* Subroutine */ int measfp_(integer *, doublereal *, doublereal *)
	    , measfs_(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal hedron;
    static integer neatmx, ispind[3];
    static doublereal totasp, alenst, alensn, vlenst, vlensn;
    extern doublereal triple_(doublereal *, doublereal *, doublereal *);
    static doublereal sumsig, sumlam, depths[50000], fnvect[450000]	/* 
	    was [3][3][50000] */, thetaq[3], totvsp, sigmaq[3];
    static logical usenum, spindl, fntrev[150000]	/* was [3][50000] */;
    extern /* Subroutine */ int measpm_(integer *, doublereal *), cerror_(
	    char *, ftnlen), gendot_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), vcross_(doublereal *, 
	    doublereal *, doublereal *), cirpln_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, logical *, logical *, 
	    doublereal *, doublereal *);
    static doublereal atmarea[25000];

    /* Fortran I/O blocks */
    static cilist io___328 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___329 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___330 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___331 = { 0, 0, 0, fmt_180, 0 };
    static cilist io___332 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___333 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___334 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___335 = { 0, 0, 0, fmt_220, 0 };
    static cilist io___336 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___337 = { 0, 0, 0, fmt_240, 0 };



#define p_ref(a_1,a_2) face06_1.p[(a_2)*3 + a_1 - 4]
#define t_ref(a_1,a_2) face05_1.t[(a_2)*3 + a_1 - 4]
#define v_ref(a_1,a_2) face07_1.v[(a_2)*3 + a_1 - 4]
#define ta_ref(a_1,a_2) face05_1.ta[(a_2)*2 + a_1 - 3]
#define env_ref(a_1,a_2) face08_1.env[(a_2)*2 + a_1 - 3]
#define tax_ref(a_1,a_2) face05_1.tax[(a_2)*3 + a_1 - 4]
#define fnt_ref(a_1,a_2) fnt[(a_2)*3 + a_1 - 4]
#define fnen_ref(a_1,a_2) face08_1.fnen[(a_2)*3 + a_1 - 4]
#define fsen_ref(a_1,a_2) face11_1.fsen[(a_2)*2 + a_1 - 3]
#define fsep_ref(a_1,a_2) face11_1.fsep[(a_2)*2 + a_1 - 3]
#define alts_ref(a_1,a_2) alts[(a_2)*3 + a_1 - 4]
#define dots_ref(a_1,a_2) dots[(a_2)*3 + a_1 - 4]
#define nspt_ref(a_1,a_2) nspt[(a_2)*3 + a_1 - 4]
#define fncen_ref(a_1,a_2) fncen[(a_2)*3 + a_1 - 4]
#define cenop_ref(a_1,a_2) cenop[(a_2)*3 + a_1 - 4]
#define vects_ref(a_1,a_2) vects[(a_2)*3 + a_1 - 4]
#define tdots_ref(a_1,a_2) tdots[(a_2)*3 + a_1 - 4]
#define fcins_ref(a_1,a_2) fcins[(a_2)*3 + a_1 - 4]
#define fcint_ref(a_1,a_2) fcint[(a_2)*3 + a_1 - 4]
#define fnvect_ref(a_1,a_2,a_3) fnvect[((a_3)*3 + (a_2))*3 + a_1 - 13]
#define fntrev_ref(a_1,a_2) fntrev[(a_2)*3 + a_1 - 4]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




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




/*     compute the volume of the interior polyhedron */

    hedron = 0.;
    i__1 = face08_1.nfn;
    for (ifn = 1; ifn <= i__1; ++ifn) {
	measpm_(&ifn, &prism);
	hedron += prism;
    }

/*     compute the area and volume due to convex faces */
/*     as well as the area partitioned among the atoms */

    totap = 0.;
    totvp = 0.;
    i__1 = face01_1.na;
    for (ia = 1; ia <= i__1; ++ia) {
	atmarea[ia - 1] = 0.;
    }
    i__1 = face13_1.nfp;
    for (ifp = 1; ifp <= i__1; ++ifp) {
	measfp_(&ifp, &areap, &volp);
	ia = face13_1.fpa[ifp - 1];
	atmarea[ia - 1] += areap;
	totap += areap;
	totvp += volp;
    }

/*     compute the area and volume due to saddle faces */
/*     as well as the spindle correction value */

    totas = 0.;
    totvs = 0.;
    totasp = 0.;
    totvsp = 0.;
    i__1 = face11_1.nfs;
    for (ifs = 1; ifs <= i__1; ++ifs) {
	for (k = 1; k <= 2; ++k) {
	    ien = fsen_ref(k, ifs);
	    if (ien > 0) {
		enfs[ien - 1] = ifs;
	    }
	}
	measfs_(&ifs, &areas, &vols, &areasp, &volsp);
	totas += areas;
	totvs += vols;
	totasp += areasp;
	totvsp += volsp;
	if (areas - areasp < 0.) {
	    cerror_("Negative Area for Saddle Face", (ftnlen)29);
	}
    }

/*     compute the area and volume due to concave faces */

    totan = 0.;
    totvn = 0.;
    i__1 = face08_1.nfn;
    for (ifn = 1; ifn <= i__1; ++ifn) {
	measfn_(&ifn, &arean, &voln);
	totan += arean;
	totvn += voln;
    }

/*     compute the area and volume lens correction values */

    alenst = 0.;
    alensn = 0.;
    vlenst = 0.;
    vlensn = 0.;
    if (face01_1.pr <= 0.) {
	goto L140;
    }
    ndots = 1000;
    gendot_(&ndots, dots, &face01_1.pr, &c_b95, &c_b95, &c_b95);
/* Computing 2nd power */
    d__1 = face01_1.pr;
    dota = d__1 * d__1 * 12.566370614359172 / ndots;
    i__1 = face08_1.nfn;
    for (ifn = 1; ifn <= i__1; ++ifn) {
	nlap[ifn - 1] = 0;
	cora[ifn - 1] = 0.;
	corv[ifn - 1] = 0.;
	badav[ifn - 1] = FALSE_;
	badt[ifn - 1] = FALSE_;
	for (k = 1; k <= 3; ++k) {
	    nspt_ref(k, ifn) = 0;
	}
	ien = fnen_ref(1, ifn);
	iv = env_ref(1, ien);
	ip = face07_1.vp[iv - 1];
	depths[ifn - 1] = depth_(&ip, &alts_ref(1, ifn));
	for (k = 1; k <= 3; ++k) {
	    fncen_ref(k, ifn) = p_ref(k, ip);
	}
	ia = face07_1.va[iv - 1];

/*     get vertices and vectors */

	for (ke = 1; ke <= 3; ++ke) {
	    ien = fnen_ref(ke, ifn);
	    ivs[ke - 1] = env_ref(1, ien);
	    ia = face07_1.va[ivs[ke - 1] - 1];
	    ifs = enfs[ien - 1];
	    iep = fsep_ref(1, ifs);
	    ic = face10_1.epc[iep - 1];
	    it = face09_1.ct[ic - 1];
	    fnt_ref(ke, ifn) = it;
	    fntrev_ref(ke, ifn) = ta_ref(1, it) != ia;
	}
	for (ke = 1; ke <= 3; ++ke) {
	    for (k = 1; k <= 3; ++k) {
		vects_ref(k, ke) = v_ref(k, ivs[ke - 1]) - p_ref(k, ip);
	    }
	}

/*     calculate normal vectors for the three planes */
/*     that cut out the geodesic triangle */

	vcross_(&vects_ref(1, 1), &vects_ref(1, 2), &fnvect_ref(1, 1, ifn));
	vnorm_(&fnvect_ref(1, 1, ifn), &fnvect_ref(1, 1, ifn));
	vcross_(&vects_ref(1, 2), &vects_ref(1, 3), &fnvect_ref(1, 2, ifn));
	vnorm_(&fnvect_ref(1, 2, ifn), &fnvect_ref(1, 2, ifn));
	vcross_(&vects_ref(1, 3), &vects_ref(1, 1), &fnvect_ref(1, 3, ifn));
	vnorm_(&fnvect_ref(1, 3, ifn), &fnvect_ref(1, 3, ifn));
    }
    i__1 = face08_1.nfn - 1;
    for (ifn = 1; ifn <= i__1; ++ifn) {
	i__2 = face08_1.nfn;
	for (jfn = ifn + 1; jfn <= i__2; ++jfn) {
	    dij2 = dist2_(&fncen_ref(1, ifn), &fncen_ref(1, jfn));
/* Computing 2nd power */
	    d__1 = face01_1.pr;
	    if (dij2 > d__1 * d__1 * 4.) {
		goto L90;
	    }
	    if (depths[ifn - 1] > face01_1.pr && depths[jfn - 1] > 
		    face01_1.pr) {
		goto L90;
	    }

/*     these two probes may have intersecting surfaces */

	    dpp = sqrt(dist2_(&fncen_ref(1, ifn), &fncen_ref(1, jfn)));

/*     compute the midpoint */

	    for (k = 1; k <= 3; ++k) {
		ppm[k - 1] = (fncen_ref(k, ifn) + fncen_ref(k, jfn)) / 2.;
		upp[k - 1] = (fncen_ref(k, jfn) - fncen_ref(k, ifn)) / dpp;
	    }
/* Computing 2nd power */
	    d__1 = face01_1.pr;
/* Computing 2nd power */
	    d__2 = dpp / 2.;
	    rm = d__1 * d__1 - d__2 * d__2;
	    if (rm < 0.) {
		rm = 0.;
	    }
	    rm = sqrt(rm);
	    rat = dpp / (face01_1.pr * 2.);
	    if (rat > 1.) {
		rat = 1.;
	    }
	    if (rat < -1.) {
		rat = -1.;
	    }
	    rho = asin(rat);

/*     use circle-plane intersection routine */

	    alli = TRUE_;
	    anyi = FALSE_;
	    spindl = FALSE_;
	    for (k = 1; k <= 3; ++k) {
		ispind[k - 1] = 0;
		ispnd2[k - 1] = 0;
	    }
	    for (ke = 1; ke <= 3; ++ke) {
		thetaq[ke - 1] = 0.;
		sigmaq[ke - 1] = 0.;
		tau[ke - 1] = 0.;
		cirpln_(ppm, &rm, upp, &fncen_ref(1, ifn), &fnvect_ref(1, ke, 
			ifn), &cinsp, &cintp, xpnt1, xpnt2);
		fcins_ref(ke, ifn) = cinsp;
		fcint_ref(ke, ifn) = cintp;
		if (! cinsp) {
		    alli = FALSE_;
		}
		if (cintp) {
		    anyi = TRUE_;
		}
		if (! cintp) {
		    goto L10;
		}
		it = fnt_ref(ke, ifn);
		if (face05_1.tr[it - 1] > face01_1.pr) {
		    goto L10;
		}
		for (ke2 = 1; ke2 <= 3; ++ke2) {
		    if (it == fnt_ref(ke2, jfn)) {
			ispind[ke - 1] = it;
			nspt_ref(ke, ifn) = nspt_ref(ke, ifn) + 1;
			ispnd2[ke2 - 1] = it;
			nspt_ref(ke2, jfn) = nspt_ref(ke2, jfn) + 1;
			spindl = TRUE_;
		    }
		}
		if (ispind[ke - 1] == 0) {
		    goto L10;
		}

/*     check that the two ways of calculating */
/*     intersection points match */

		rat = face05_1.tr[it - 1] / face01_1.pr;
		if (rat > 1.) {
		    rat = 1.;
		}
		if (rat < -1.) {
		    rat = -1.;
		}
		thetaq[ke - 1] = acos(rat);
		stq = sin(thetaq[ke - 1]);
		if (fntrev_ref(ke, ifn)) {
		    for (k = 1; k <= 3; ++k) {
			uij[k - 1] = -tax_ref(k, it);
		    }
		} else {
		    for (k = 1; k <= 3; ++k) {
			uij[k - 1] = tax_ref(k, it);
		    }
		}
		for (k = 1; k <= 3; ++k) {
		    qij[k - 1] = t_ref(k, it) - stq * face01_1.pr * uij[k - 1]
			    ;
		    qji[k - 1] = t_ref(k, it) + stq * face01_1.pr * uij[k - 1]
			    ;
		}
		for (k = 1; k <= 3; ++k) {
		    umq[k - 1] = (qij[k - 1] - ppm[k - 1]) / rm;
		    upq[k - 1] = (qij[k - 1] - fncen_ref(k, ifn)) / 
			    face01_1.pr;
		}
		vcross_(uij, upp, vect1);
		dt = dot_(umq, vect1);
		if (dt > 1.) {
		    dt = 1.;
		}
		if (dt < -1.) {
		    dt = -1.;
		}
		sigmaq[ke - 1] = acos(dt);
		vcross_(upq, &fnvect_ref(1, ke, ifn), vect1);
		vnorm_(vect1, uc);
		vcross_(upp, upq, vect1);
		vnorm_(vect1, uq);
		dt = dot_(uc, uq);
		if (dt > 1.) {
		    dt = 1.;
		}
		if (dt < -1.) {
		    dt = -1.;
		}
		tau[ke - 1] = 3.141592653589793238 - acos(dt);
L10:
		;
	    }
	    allj = TRUE_;
	    anyj = FALSE_;
	    for (ke = 1; ke <= 3; ++ke) {
		cirpln_(ppm, &rm, upp, &fncen_ref(1, jfn), &fnvect_ref(1, ke, 
			jfn), &cinsp, &cintp, xpnt1, xpnt2);
		fcins_ref(ke, jfn) = cinsp;
		fcint_ref(ke, jfn) = cintp;
		if (! cinsp) {
		    allj = FALSE_;
		}
		if (cintp) {
		    anyj = TRUE_;
		}
	    }
	    case1 = alli && allj && ! anyi && ! anyj;
	    case2 = anyi && anyj && spindl;
	    if (! case1 && ! case2) {
		goto L90;
	    }

/*     this kind of overlap can be handled */

	    ++nlap[ifn - 1];
	    ++nlap[jfn - 1];
	    for (ke = 1; ke <= 3; ++ke) {
		ien = fnen_ref(ke, ifn);
		iv1 = env_ref(1, ien);
		iv2 = env_ref(2, ien);
		for (k = 1; k <= 3; ++k) {
		    vect3[k - 1] = v_ref(k, iv1) - fncen_ref(k, ifn);
		    vect4[k - 1] = v_ref(k, iv2) - fncen_ref(k, ifn);
		}
		for (ke2 = 1; ke2 <= 3; ++ke2) {
		    if (ispind[ke - 1] == ispnd2[ke2 - 1]) {
			goto L40;
		    }
		    if (ispind[ke - 1] == 0) {
			goto L40;
		    }
		    cirpln_(&fncen_ref(1, ifn), &face01_1.pr, &fnvect_ref(1, 
			    ke, ifn), &fncen_ref(1, jfn), &fnvect_ref(1, ke2, 
			    jfn), &cinsp, &cintp, xpnt1, xpnt2);
		    if (! cintp) {
			goto L40;
		    }
		    ien = fnen_ref(ke2, jfn);
		    iv1 = env_ref(1, ien);
		    iv2 = env_ref(2, ien);
		    for (k = 1; k <= 3; ++k) {
			vect7[k - 1] = v_ref(k, iv1) - fncen_ref(k, jfn);
			vect8[k - 1] = v_ref(k, iv2) - fncen_ref(k, jfn);
		    }

/*     check whether point lies on spindle arc */

		    for (k = 1; k <= 3; ++k) {
			vect1[k - 1] = xpnt1[k - 1] - fncen_ref(k, ifn);
			vect2[k - 1] = xpnt2[k - 1] - fncen_ref(k, ifn);
			vect5[k - 1] = xpnt1[k - 1] - fncen_ref(k, jfn);
			vect6[k - 1] = xpnt2[k - 1] - fncen_ref(k, jfn);
		    }
		    if (triple_(vect3, vect1, &fnvect_ref(1, ke, ifn)) < 0.) {
			goto L20;
		    }
		    if (triple_(vect1, vect4, &fnvect_ref(1, ke, ifn)) < 0.) {
			goto L20;
		    }
		    if (triple_(vect7, vect5, &fnvect_ref(1, ke2, jfn)) < 0.) 
			    {
			goto L20;
		    }
		    if (triple_(vect5, vect8, &fnvect_ref(1, ke2, jfn)) < 0.) 
			    {
			goto L20;
		    }
		    goto L30;
L20:
		    if (triple_(vect3, vect2, &fnvect_ref(1, ke, ifn)) < 0.) {
			goto L40;
		    }
		    if (triple_(vect2, vect4, &fnvect_ref(1, ke, ifn)) < 0.) {
			goto L40;
		    }
		    if (triple_(vect7, vect6, &fnvect_ref(1, ke2, jfn)) < 0.) 
			    {
			goto L40;
		    }
		    if (triple_(vect6, vect8, &fnvect_ref(1, ke2, jfn)) < 0.) 
			    {
			goto L40;
		    }
L30:
		    badav[ifn - 1] = TRUE_;
L40:
		    ;
		}
	    }
	    for (ke = 1; ke <= 3; ++ke) {
		ien = fnen_ref(ke, ifn);
		iv1 = env_ref(1, ien);
		iv2 = env_ref(2, ien);
		for (k = 1; k <= 3; ++k) {
		    vect3[k - 1] = v_ref(k, iv1) - fncen_ref(k, ifn);
		    vect4[k - 1] = v_ref(k, iv2) - fncen_ref(k, ifn);
		}
		for (ke2 = 1; ke2 <= 3; ++ke2) {
		    if (ispind[ke - 1] == ispnd2[ke2 - 1]) {
			goto L70;
		    }
		    if (ispnd2[ke2 - 1] == 0) {
			goto L70;
		    }
		    cirpln_(&fncen_ref(1, jfn), &face01_1.pr, &fnvect_ref(1, 
			    ke2, jfn), &fncen_ref(1, ifn), &fnvect_ref(1, ke, 
			    ifn), &cinsp, &cintp, xpnt1, xpnt2);
		    if (! cintp) {
			goto L70;
		    }
		    ien = fnen_ref(ke2, jfn);
		    iv1 = env_ref(1, ien);
		    iv2 = env_ref(2, ien);
		    for (k = 1; k <= 3; ++k) {
			vect7[k - 1] = v_ref(k, iv1) - fncen_ref(k, jfn);
			vect8[k - 1] = v_ref(k, iv2) - fncen_ref(k, jfn);
		    }

/*     check whether point lies on spindle arc */

		    for (k = 1; k <= 3; ++k) {
			vect1[k - 1] = xpnt1[k - 1] - fncen_ref(k, ifn);
			vect2[k - 1] = xpnt2[k - 1] - fncen_ref(k, ifn);
			vect5[k - 1] = xpnt1[k - 1] - fncen_ref(k, jfn);
			vect6[k - 1] = xpnt2[k - 1] - fncen_ref(k, jfn);
		    }
		    if (triple_(vect3, vect1, &fnvect_ref(1, ke, ifn)) < 0.) {
			goto L50;
		    }
		    if (triple_(vect1, vect4, &fnvect_ref(1, ke, ifn)) < 0.) {
			goto L50;
		    }
		    if (triple_(vect7, vect5, &fnvect_ref(1, ke2, jfn)) < 0.) 
			    {
			goto L50;
		    }
		    if (triple_(vect5, vect8, &fnvect_ref(1, ke2, jfn)) < 0.) 
			    {
			goto L50;
		    }
		    goto L60;
L50:
		    if (triple_(vect3, vect2, &fnvect_ref(1, ke, ifn)) < 0.) {
			goto L70;
		    }
		    if (triple_(vect2, vect4, &fnvect_ref(1, ke, ifn)) < 0.) {
			goto L70;
		    }
		    if (triple_(vect7, vect6, &fnvect_ref(1, ke2, jfn)) < 0.) 
			    {
			goto L70;
		    }
		    if (triple_(vect6, vect8, &fnvect_ref(1, ke2, jfn)) < 0.) 
			    {
			goto L70;
		    }
L60:
		    badav[jfn - 1] = TRUE_;
L70:
		    ;
		}
	    }
	    sumlam = 0.;
	    sumsig = 0.;
	    sumsc = 0.;
	    for (ke = 1; ke <= 3; ++ke) {
		if (ispind[ke - 1] != 0) {
		    sumlam = sumlam + 3.141592653589793238 - tau[ke - 1];
		    sumsig = sumsig + sigmaq[ke - 1] - 3.141592653589793238;
		    sumsc += sin(sigmaq[ke - 1]) * cos(sigmaq[ke - 1]);
		}
	    }
/* Computing 2nd power */
	    d__1 = face01_1.pr;
	    alens = d__1 * d__1 * 2. * (3.141592653589793238 - sumlam - sin(
		    rho) * (sumsig + 3.141592653589793238));
	    vint = alens * face01_1.pr / 3.;
/* Computing 2nd power */
	    d__1 = rm;
	    vcone = face01_1.pr * (d__1 * d__1) * sin(rho) * (sumsig + 
		    3.141592653589793238) / 3.;
/* Computing 2nd power */
	    d__1 = rm;
	    vpyr = face01_1.pr * (d__1 * d__1) * sin(rho) * sumsc / 3.;
	    vlens = vint - vcone + vpyr;
	    cora[ifn - 1] += alens;
	    cora[jfn - 1] += alens;
	    corv[ifn - 1] += vlens;
	    corv[jfn - 1] += vlens;

/*     check for vertex on opposing probe in face */

	    for (kv = 1; kv <= 3; ++kv) {
		vip[kv - 1] = FALSE_;
		ien = fnen_ref(kv, jfn);
		iv = env_ref(1, ien);
		for (k = 1; k <= 3; ++k) {
		    vect1[k - 1] = v_ref(k, iv) - fncen_ref(k, ifn);
		}
		vnorm_(vect1, vect1);
		for (ke = 1; ke <= 3; ++ke) {
		    dt = dot_(&fnvect_ref(1, ke, ifn), &v_ref(1, iv));
		    if (dt > 0.) {
			goto L80;
		    }
		}
		vip[kv - 1] = TRUE_;
L80:
		;
	    }
L90:
	    ;
	}
    }
    i__1 = face08_1.nfn;
    for (ifn = 1; ifn <= i__1; ++ifn) {
	for (ke = 1; ke <= 3; ++ke) {
	    if (nspt_ref(ke, ifn) > 1) {
		badt[ifn - 1] = TRUE_;
	    }
	}
    }
    i__1 = face08_1.nfn;
    for (ifn = 1; ifn <= i__1; ++ifn) {
	if (nlap[ifn - 1] <= 0) {
	    goto L130;
	}

/*     gather all overlapping probes */

	nop = 0;
	i__2 = face08_1.nfn;
	for (jfn = 1; jfn <= i__2; ++jfn) {
	    if (ifn != jfn) {
		dij2 = dist2_(&fncen_ref(1, ifn), &fncen_ref(1, jfn));
/* Computing 2nd power */
		d__1 = face01_1.pr;
		if (dij2 <= d__1 * d__1 * 4.) {
		    if (depths[jfn - 1] <= face01_1.pr) {
			++nop;
			if (nop > 100) {
			    cerror_("NOP Overflow in VAM", (ftnlen)19);
			}
			ifnop[nop - 1] = jfn;
			for (k = 1; k <= 3; ++k) {
			    cenop_ref(k, nop) = fncen_ref(k, jfn);
			}
		    }
		}
	    }
	}

/*     numerical calculation of the correction */

	areado = 0.;
	voldo = 0.;
	scinc = .050000000000000003;
	for (isc = 1; isc <= 20; ++isc) {
	    rsc = isc - .5;
/* Computing 2nd power */
	    d__1 = rsc;
/* Computing 3rd power */
	    d__2 = scinc;
	    dotv[isc - 1] = face01_1.pr * dota * (d__1 * d__1) * (d__2 * (
		    d__2 * d__2));
	}
	i__2 = nop;
	for (iop = 1; iop <= i__2; ++iop) {
	    ate[iop - 1] = FALSE_;
	}
	neatmx = 0;
	i__2 = ndots;
	for (idot = 1; idot <= i__2; ++idot) {
	    for (ke = 1; ke <= 3; ++ke) {
		dt = dot_(&fnvect_ref(1, ke, ifn), &dots_ref(1, idot));
		if (dt > 0.) {
		    goto L120;
		}
	    }
	    for (k = 1; k <= 3; ++k) {
		tdots_ref(k, idot) = fncen_ref(k, ifn) + dots_ref(k, idot);
	    }
	    i__3 = nop;
	    for (iop = 1; iop <= i__3; ++iop) {
		jfn = ifnop[iop - 1];
		ds2 = dist2_(&tdots_ref(1, idot), &fncen_ref(1, jfn));
/* Computing 2nd power */
		d__1 = face01_1.pr;
		if (ds2 < d__1 * d__1) {
		    areado += dota;
		    goto L100;
		}
	    }
L100:
	    for (isc = 1; isc <= 20; ++isc) {
		rsc = isc - .5;
		for (k = 1; k <= 3; ++k) {
		    sdot[k - 1] = fncen_ref(k, ifn) + rsc * scinc * dots_ref(
			    k, idot);
		}
		neat = 0;
		i__3 = nop;
		for (iop = 1; iop <= i__3; ++iop) {
		    jfn = ifnop[iop - 1];
		    ds2 = dist2_(sdot, &fncen_ref(1, jfn));
/* Computing 2nd power */
		    d__1 = face01_1.pr;
		    if (ds2 < d__1 * d__1) {
			for (k = 1; k <= 3; ++k) {
			    vect1[k - 1] = sdot[k - 1] - fncen_ref(k, jfn);
			}
			for (ke = 1; ke <= 3; ++ke) {
			    dt = dot_(&fnvect_ref(1, ke, jfn), vect1);
			    if (dt > 0.) {
				goto L110;
			    }
			}
			++neat;
			ate[iop - 1] = TRUE_;
L110:
			;
		    }
		}
		if (neat > neatmx) {
		    neatmx = neat;
		}
		if (neat > 0) {
		    voldo += dotv[isc - 1] * (neat / (neat + 1.));
		}
	    }
L120:
	    ;
	}
	coran = areado;
	corvn = voldo;
	nate = 0;
	i__2 = nop;
	for (iop = 1; iop <= i__2; ++iop) {
	    if (ate[iop - 1]) {
		++nate;
	    }
	}

/*     use either the analytical or numerical correction */

	usenum = nate > nlap[ifn - 1] || neatmx > 1 || badt[ifn - 1];
	if (usenum) {
	    cora[ifn - 1] = coran;
	    corv[ifn - 1] = corvn;
	    alensn += cora[ifn - 1];
	    vlensn += corv[ifn - 1];
	} else if (badav[ifn - 1]) {
	    corv[ifn - 1] = corvn;
	    vlensn += corv[ifn - 1];
	}
	alenst += cora[ifn - 1];
	vlenst += corv[ifn - 1];
L130:
	;
    }
L140:

/*     print out the decomposition of the area and volume */

    if (inform_1.debug) {
	io___328.ciunit = iounit_1.iout;
	s_wsfe(&io___328);
	e_wsfe();
	k = 1;
	while(k <= face01_1.na) {
	    io___329.ciunit = iounit_1.iout;
	    s_wsfe(&io___329);
/* Computing MIN */
	    i__2 = k + 4;
	    i__1 = min(i__2,face01_1.na);
	    for (ia = k; ia <= i__1; ++ia) {
		do_fio(&c__1, (char *)&ia, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&atmarea[ia - 1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	    k += 5;
	}
	io___330.ciunit = iounit_1.iout;
	s_wsfe(&io___330);
	e_wsfe();
	io___331.ciunit = iounit_1.iout;
	s_wsfe(&io___331);
	do_fio(&c__1, (char *)&face13_1.nfp, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&totap, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&totvp, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___332.ciunit = iounit_1.iout;
	s_wsfe(&io___332);
	do_fio(&c__1, (char *)&face11_1.nfs, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&totas, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&totvs, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___333.ciunit = iounit_1.iout;
	s_wsfe(&io___333);
	do_fio(&c__1, (char *)&face08_1.nfn, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&totan, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&totvn, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___334.ciunit = iounit_1.iout;
	s_wsfe(&io___334);
	do_fio(&c__1, (char *)&hedron, (ftnlen)sizeof(doublereal));
	e_wsfe();
	if (totasp != 0. || totvsp != 0. || alenst != 0. || vlenst != 0.) {
	    io___335.ciunit = iounit_1.iout;
	    s_wsfe(&io___335);
	    d__1 = -totasp;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    d__2 = -totvsp;
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    io___336.ciunit = iounit_1.iout;
	    s_wsfe(&io___336);
	    d__1 = -alenst - alensn;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    d__2 = vlenst - vlensn;
	    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	if (alensn != 0. || vlensn != 0.) {
	    io___337.ciunit = iounit_1.iout;
	    s_wsfe(&io___337);
	    do_fio(&c__1, (char *)&alensn, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&vlensn, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*     finally, compute the total area and total volume */

    *area = totap + totas + totan - totasp - alenst;
    *volume = totvp + totvs + totvn + hedron - totvsp + vlenst;
    return 0;
} /* vam_ */

#undef fntrev_ref
#undef fnvect_ref
#undef fcint_ref
#undef fcins_ref
#undef tdots_ref
#undef vects_ref
#undef cenop_ref
#undef fncen_ref
#undef nspt_ref
#undef dots_ref
#undef alts_ref
#undef fsep_ref
#undef fsen_ref
#undef fnen_ref
#undef fnt_ref
#undef tax_ref
#undef env_ref
#undef ta_ref
#undef v_ref
#undef t_ref
#undef p_ref




/*     ###################### */
/*     ##                  ## */
/*     ##  function depth  ## */
/*     ##                  ## */
/*     ###################### */


doublereal depth_(integer *ip, doublereal *alt)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer k, ia1, ia2, ia3;
    extern doublereal dot_(doublereal *, doublereal *);
    static doublereal vect1[3], vect2[3], vect3[3], vect4[3];
    extern /* Subroutine */ int vnorm_(doublereal *, doublereal *), vcross_(
	    doublereal *, doublereal *, doublereal *);


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]
#define p_ref(a_1,a_2) face06_1.p[(a_2)*3 + a_1 - 4]
#define pa_ref(a_1,a_2) face06_1.pa[(a_2)*3 + a_1 - 4]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




    /* Parameter adjustments */
    --alt;

    /* Function Body */
    ia1 = pa_ref(1, *ip);
    ia2 = pa_ref(2, *ip);
    ia3 = pa_ref(3, *ip);
    for (k = 1; k <= 3; ++k) {
	vect1[k - 1] = a_ref(k, ia1) - a_ref(k, ia3);
	vect2[k - 1] = a_ref(k, ia2) - a_ref(k, ia3);
	vect3[k - 1] = p_ref(k, *ip) - a_ref(k, ia3);
    }
    vcross_(vect1, vect2, vect4);
    vnorm_(vect4, vect4);
    ret_val = dot_(vect4, vect3);
    for (k = 1; k <= 3; ++k) {
	alt[k] = vect4[k - 1];
    }
    return ret_val;
} /* depth_ */

#undef pa_ref
#undef p_ref
#undef a_ref




/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  subroutine measpm  --  volume of interior polyhedron  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     "measpm" computes the volume of a single prism section of */
/*     the full interior polyhedron */


/* Subroutine */ int measpm_(integer *ifn, doublereal *prism)
{
    static integer k, ia, ke, ip, iv, ien;
    static doublereal pav[9]	/* was [3][3] */, vect1[3], vect2[3], vect3[3]
	    , height;
    extern /* Subroutine */ int vcross_(doublereal *, doublereal *, 
	    doublereal *);


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]
#define p_ref(a_1,a_2) face06_1.p[(a_2)*3 + a_1 - 4]
#define pav_ref(a_1,a_2) pav[(a_2)*3 + a_1 - 4]
#define env_ref(a_1,a_2) face08_1.env[(a_2)*2 + a_1 - 3]
#define fnen_ref(a_1,a_2) face08_1.fnen[(a_2)*3 + a_1 - 4]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




    height = 0.;
    for (ke = 1; ke <= 3; ++ke) {
	ien = fnen_ref(ke, *ifn);
	iv = env_ref(1, ien);
	ia = face07_1.va[iv - 1];
	height += a_ref(3, ia);
	ip = face07_1.vp[iv - 1];
	for (k = 1; k <= 3; ++k) {
	    pav_ref(k, ke) = a_ref(k, ia) - p_ref(k, ip);
	}
    }
    height /= 3.;
    for (k = 1; k <= 3; ++k) {
	vect1[k - 1] = pav_ref(k, 2) - pav_ref(k, 1);
	vect2[k - 1] = pav_ref(k, 3) - pav_ref(k, 1);
    }
    vcross_(vect1, vect2, vect3);
    *prism = height * vect3[2] / 2.;
    return 0;
} /* measpm_ */

#undef fnen_ref
#undef env_ref
#undef pav_ref
#undef p_ref
#undef a_ref




/*     ######################### */
/*     ##                     ## */
/*     ##  subroutine measfp  ## */
/*     ##                     ## */
/*     ######################### */


/* Subroutine */ int measfp_(integer *ifp, doublereal *areap, doublereal *
	volp)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer k, ia, ic, ke;
    static doublereal dt;
    static integer it, ia2, iv1, iv2;
    static doublereal geo;
    static integer iep, icy;
    extern doublereal dot_(doublereal *, doublereal *);
    static doublereal tanv[180]	/* was [3][2][30] */, vect1[3], vect2[3];
    static integer nedge;
    static doublereal angle, gauss;
    extern /* Subroutine */ int vnorm_(doublereal *, doublereal *);
    static doublereal radial[90]	/* was [3][30] */;
    extern doublereal vecang_(doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal aavect[3], acvect[3];
    static integer ncycle, ieuler;
    static doublereal pcurve;
    static integer icyptr;
    static doublereal gcurve;
    extern /* Subroutine */ int cerror_(char *, ftnlen), vcross_(doublereal *,
	     doublereal *, doublereal *);


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]
#define c___ref(a_1,a_2) face09_1.c__[(a_2)*3 + a_1 - 4]
#define v_ref(a_1,a_2) face07_1.v[(a_2)*3 + a_1 - 4]
#define ta_ref(a_1,a_2) face05_1.ta[(a_2)*2 + a_1 - 3]
#define epv_ref(a_1,a_2) face10_1.epv[(a_2)*2 + a_1 - 3]
#define cyep_ref(a_1,a_2) face12_1.cyep[(a_2)*30 + a_1 - 31]
#define fpcy_ref(a_1,a_2) face13_1.fpcy[(a_2)*10 + a_1 - 11]
#define tanv_ref(a_1,a_2,a_3) tanv[((a_3)*2 + (a_2))*3 + a_1 - 10]
#define radial_ref(a_1,a_2) radial[(a_2)*3 + a_1 - 4]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




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




    ia = face13_1.fpa[*ifp - 1];
    pcurve = 0.;
    gcurve = 0.;
    ncycle = face13_1.fpncy[*ifp - 1];
    if (ncycle > 0) {
	ieuler = 2 - ncycle;
    } else {
	ieuler = 2;
    }
    i__1 = ncycle;
    for (icyptr = 1; icyptr <= i__1; ++icyptr) {
	icy = fpcy_ref(icyptr, *ifp);
	nedge = face12_1.cynep[icy - 1];
	i__2 = nedge;
	for (ke = 1; ke <= i__2; ++ke) {
	    iep = cyep_ref(ke, icy);
	    ic = face10_1.epc[iep - 1];
	    it = face09_1.ct[ic - 1];
	    if (ia == ta_ref(1, it)) {
		ia2 = ta_ref(2, it);
	    } else {
		ia2 = ta_ref(1, it);
	    }
	    for (k = 1; k <= 3; ++k) {
		acvect[k - 1] = c___ref(k, ic) - a_ref(k, ia);
		aavect[k - 1] = a_ref(k, ia2) - a_ref(k, ia);
	    }
	    vnorm_(aavect, aavect);
	    dt = dot_(acvect, aavect);
	    geo = -dt / (face01_1.ar[ia - 1] * face09_1.cr[ic - 1]);
	    iv1 = epv_ref(1, iep);
	    iv2 = epv_ref(2, iep);
	    if (iv1 == 0 || iv2 == 0) {
		angle = 6.2831853071795862;
	    } else {
		for (k = 1; k <= 3; ++k) {
		    vect1[k - 1] = v_ref(k, iv1) - c___ref(k, ic);
		    vect2[k - 1] = v_ref(k, iv2) - c___ref(k, ic);
		    radial_ref(k, ke) = v_ref(k, iv1) - a_ref(k, ia);
		}
		vnorm_(&radial_ref(1, ke), &radial_ref(1, ke));
		vcross_(vect1, aavect, &tanv_ref(1, 1, ke));
		vnorm_(&tanv_ref(1, 1, ke), &tanv_ref(1, 1, ke));
		vcross_(vect2, aavect, &tanv_ref(1, 2, ke));
		vnorm_(&tanv_ref(1, 2, ke), &tanv_ref(1, 2, ke));
		angle = vecang_(vect1, vect2, aavect, &c_b154);
	    }
	    gcurve += face09_1.cr[ic - 1] * angle * geo;
	    if (nedge != 1) {
		if (ke > 1) {
		    angle = vecang_(&tanv_ref(1, 2, ke - 1), &tanv_ref(1, 1, 
			    ke), &radial_ref(1, ke), &c_b155);
		    if (angle < 0.) {
			cerror_("Negative Angle in MEASFP", (ftnlen)24);
		    }
		    pcurve += angle;
		}
	    }
	}
	if (nedge > 1) {
	    angle = vecang_(&tanv_ref(1, 2, nedge), &tanv_ref(1, 1, 1), &
		    radial_ref(1, 1), &c_b155);
	    if (angle < 0.) {
		cerror_("Negative Angle in MEASFP", (ftnlen)24);
	    }
	    pcurve += angle;
	}
    }
    gauss = ieuler * 6.2831853071795862 - pcurve - gcurve;
/* Computing 2nd power */
    d__1 = face01_1.ar[ia - 1];
    *areap = gauss * (d__1 * d__1);
    *volp = *areap * face01_1.ar[ia - 1] / 3.;
    return 0;
} /* measfp_ */

#undef radial_ref
#undef tanv_ref
#undef fpcy_ref
#undef cyep_ref
#undef epv_ref
#undef ta_ref
#undef v_ref
#undef c___ref
#undef a_ref




/*     ######################### */
/*     ##                     ## */
/*     ##  subroutine measfs  ## */
/*     ##                     ## */
/*     ######################### */


/* Subroutine */ int measfs_(integer *ifs, doublereal *areas, doublereal *
	vols, doublereal *areasp, doublereal *volsp)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double atan2(doublereal, doublereal), acos(doublereal), sin(doublereal), 
	    cos(doublereal);

    /* Local variables */
    static integer k;
    static doublereal d1, d2, w1, w2;
    static integer ic, it, ia1, ia2, ic1, ic2, iv1, iv2, iep;
    static doublereal phi;
    extern doublereal dot_(doublereal *, doublereal *);
    static doublereal rat, spin;
    static logical cusp;
    static doublereal volt, cone1, cone2, vect1[3], vect2[3], term1, term2, 
	    term3;
    extern /* Subroutine */ int vnorm_(doublereal *, doublereal *);
    static doublereal theta1, theta2;
    extern doublereal vecang_(doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal aavect[3], thetaq;
    extern /* Subroutine */ int cerror_(char *, ftnlen);


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]
#define c___ref(a_1,a_2) face09_1.c__[(a_2)*3 + a_1 - 4]
#define t_ref(a_1,a_2) face05_1.t[(a_2)*3 + a_1 - 4]
#define v_ref(a_1,a_2) face07_1.v[(a_2)*3 + a_1 - 4]
#define ta_ref(a_1,a_2) face05_1.ta[(a_2)*2 + a_1 - 3]
#define epv_ref(a_1,a_2) face10_1.epv[(a_2)*2 + a_1 - 3]
#define fsep_ref(a_1,a_2) face11_1.fsep[(a_2)*2 + a_1 - 3]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




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




    iep = fsep_ref(1, *ifs);
    ic = face10_1.epc[iep - 1];
    it = face09_1.ct[ic - 1];
    ia1 = ta_ref(1, it);
    ia2 = ta_ref(2, it);
    for (k = 1; k <= 3; ++k) {
	aavect[k - 1] = a_ref(k, ia2) - a_ref(k, ia1);
    }
    vnorm_(aavect, aavect);
    iv1 = epv_ref(1, iep);
    iv2 = epv_ref(2, iep);
    if (iv1 == 0 || iv2 == 0) {
	phi = 6.2831853071795862;
    } else {
	for (k = 1; k <= 3; ++k) {
	    vect1[k - 1] = v_ref(k, iv1) - c___ref(k, ic);
	    vect2[k - 1] = v_ref(k, iv2) - c___ref(k, ic);
	}
	phi = vecang_(vect1, vect2, aavect, &c_b155);
    }
    for (k = 1; k <= 3; ++k) {
	vect1[k - 1] = a_ref(k, ia1) - t_ref(k, it);
	vect2[k - 1] = a_ref(k, ia2) - t_ref(k, it);
    }
    d1 = -dot_(vect1, aavect);
    d2 = dot_(vect2, aavect);
    theta1 = atan2(d1, face05_1.tr[it - 1]);
    theta2 = atan2(d2, face05_1.tr[it - 1]);

/*     check for cusps */

    if (face05_1.tr[it - 1] < face01_1.pr && theta1 > 0. && theta2 > 0.) {
	cusp = TRUE_;
	rat = face05_1.tr[it - 1] / face01_1.pr;
	if (rat > 1.) {
	    rat = 1.;
	}
	if (rat < -1.) {
	    rat = -1.;
	}
	thetaq = acos(rat);
    } else {
	cusp = FALSE_;
	thetaq = 0.;
	*areasp = 0.;
	*volsp = 0.;
    }
    term1 = face05_1.tr[it - 1] * face01_1.pr * (theta1 + theta2);
/* Computing 2nd power */
    d__1 = face01_1.pr;
    term2 = d__1 * d__1 * (sin(theta1) + sin(theta2));
    *areas = phi * (term1 - term2);
    if (cusp) {
/* Computing 2nd power */
	d__1 = face01_1.pr;
	spin = face05_1.tr[it - 1] * face01_1.pr * thetaq - d__1 * d__1 * sin(
		thetaq);
	*areasp = phi * 2. * spin;
    }

    iep = fsep_ref(1, *ifs);
    ic2 = face10_1.epc[iep - 1];
    iep = fsep_ref(2, *ifs);
    ic1 = face10_1.epc[iep - 1];
    if (face09_1.ca[ic1 - 1] != ia1) {
	cerror_("IA1 Inconsistency in MEASFS", (ftnlen)27);
    }
    for (k = 1; k <= 3; ++k) {
	vect1[k - 1] = c___ref(k, ic1) - a_ref(k, ia1);
	vect2[k - 1] = c___ref(k, ic2) - a_ref(k, ia2);
    }
    w1 = dot_(vect1, aavect);
    w2 = -dot_(vect2, aavect);
/* Computing 2nd power */
    d__1 = face09_1.cr[ic1 - 1];
    cone1 = phi * (w1 * (d__1 * d__1)) / 6.;
/* Computing 2nd power */
    d__1 = face09_1.cr[ic2 - 1];
    cone2 = phi * (w2 * (d__1 * d__1)) / 6.;
/* Computing 2nd power */
    d__1 = face05_1.tr[it - 1];
    term1 = d__1 * d__1 * face01_1.pr * (sin(theta1) + sin(theta2));
    term2 = sin(theta1) * cos(theta1) + theta1 + sin(theta2) * cos(theta2) + 
	    theta2;
/* Computing 2nd power */
    d__1 = face01_1.pr;
    term2 = face05_1.tr[it - 1] * (d__1 * d__1) * term2;
/* Computing 2nd power */
    d__1 = cos(theta1);
/* Computing 2nd power */
    d__2 = cos(theta2);
    term3 = sin(theta1) * (d__1 * d__1) + sin(theta1) * 2. + sin(theta2) * (
	    d__2 * d__2) + sin(theta2) * 2.;
/* Computing 3rd power */
    d__1 = face01_1.pr;
    term3 = d__1 * (d__1 * d__1) / 3. * term3;
    volt = phi / 2. * (term1 - term2 + term3);
    *vols = volt + cone1 + cone2;
    if (cusp) {
/* Computing 2nd power */
	d__1 = face05_1.tr[it - 1];
	term1 = d__1 * d__1 * face01_1.pr * sin(thetaq);
	term2 = sin(thetaq) * cos(thetaq) + thetaq;
/* Computing 2nd power */
	d__1 = face01_1.pr;
	term2 = face05_1.tr[it - 1] * (d__1 * d__1) * term2;
/* Computing 2nd power */
	d__1 = cos(thetaq);
	term3 = sin(thetaq) * (d__1 * d__1) + sin(thetaq) * 2.;
/* Computing 3rd power */
	d__1 = face01_1.pr;
	term3 = d__1 * (d__1 * d__1) / 3. * term3;
	*volsp = phi * (term1 - term2 + term3);
    }
    return 0;
} /* measfs_ */

#undef fsep_ref
#undef epv_ref
#undef ta_ref
#undef v_ref
#undef t_ref
#undef c___ref
#undef a_ref




/*     ######################### */
/*     ##                     ## */
/*     ##  subroutine measfn  ## */
/*     ##                     ## */
/*     ######################### */


/* Subroutine */ int measfn_(integer *ifn, doublereal *arean, doublereal *
	voln)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static integer k, ia, ke, je, ip, iv, ien;
    static doublereal pav[9]	/* was [3][3] */, pvv[9]	/* was [3][3] 
	    */, angle[3];
    extern /* Subroutine */ int vnorm_(doublereal *, doublereal *);
    static doublereal defect;
    extern doublereal vecang_(doublereal *, doublereal *, doublereal *, 
	    doublereal *), triple_(doublereal *, doublereal *, doublereal *);
    static doublereal planev[9]	/* was [3][3] */;
    extern /* Subroutine */ int cerror_(char *, ftnlen);
    static doublereal simplx;
    extern /* Subroutine */ int vcross_(doublereal *, doublereal *, 
	    doublereal *);


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]
#define p_ref(a_1,a_2) face06_1.p[(a_2)*3 + a_1 - 4]
#define v_ref(a_1,a_2) face07_1.v[(a_2)*3 + a_1 - 4]
#define pav_ref(a_1,a_2) pav[(a_2)*3 + a_1 - 4]
#define env_ref(a_1,a_2) face08_1.env[(a_2)*2 + a_1 - 3]
#define pvv_ref(a_1,a_2) pvv[(a_2)*3 + a_1 - 4]
#define fnen_ref(a_1,a_2) face08_1.fnen[(a_2)*3 + a_1 - 4]
#define planev_ref(a_1,a_2) planev[(a_2)*3 + a_1 - 4]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




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




    for (ke = 1; ke <= 3; ++ke) {
	ien = fnen_ref(ke, *ifn);
	iv = env_ref(1, ien);
	ia = face07_1.va[iv - 1];
	ip = face07_1.vp[iv - 1];
	for (k = 1; k <= 3; ++k) {
	    pvv_ref(k, ke) = v_ref(k, iv) - p_ref(k, ip);
	    pav_ref(k, ke) = a_ref(k, ia) - p_ref(k, ip);
	}
	if (face01_1.pr > 0.) {
	    vnorm_(&pvv_ref(1, ke), &pvv_ref(1, ke));
	}
    }
    if (face01_1.pr <= 0.) {
	*arean = 0.;
    } else {
	for (ke = 1; ke <= 3; ++ke) {
	    je = ke + 1;
	    if (je > 3) {
		je = 1;
	    }
	    vcross_(&pvv_ref(1, ke), &pvv_ref(1, je), &planev_ref(1, ke));
	    vnorm_(&planev_ref(1, ke), &planev_ref(1, ke));
	}
	for (ke = 1; ke <= 3; ++ke) {
	    je = ke - 1;
	    if (je < 1) {
		je = 3;
	    }
	    angle[ke - 1] = vecang_(&planev_ref(1, je), &planev_ref(1, ke), &
		    pvv_ref(1, ke), &c_b154);
	    if (angle[ke - 1] < 0.) {
		cerror_("Negative Angle in MEASFN", (ftnlen)24);
	    }
	}
	defect = 6.2831853071795862 - (angle[0] + angle[1] + angle[2]);
/* Computing 2nd power */
	d__1 = face01_1.pr;
	*arean = d__1 * d__1 * defect;
    }
    simplx = -triple_(&pav_ref(1, 1), &pav_ref(1, 2), &pav_ref(1, 3)) / 6.;
    *voln = simplx - *arean * face01_1.pr / 3.;
    return 0;
} /* measfn_ */

#undef planev_ref
#undef fnen_ref
#undef pvv_ref
#undef env_ref
#undef pav_ref
#undef v_ref
#undef p_ref
#undef a_ref




/*     ######################### */
/*     ##                     ## */
/*     ##  subroutine projct  ## */
/*     ##                     ## */
/*     ######################### */


/* Subroutine */ int projct_(doublereal *pnt, doublereal *unvect, integer *
	icy, integer *ia, doublereal *spv, integer *nedge, logical *fail)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal f;
    static integer k, ke;
    static doublereal dt;
    static integer iv, iep;
    extern doublereal dot_(doublereal *, doublereal *);
    static doublereal polev[3];


#define v_ref(a_1,a_2) face07_1.v[(a_2)*3 + a_1 - 4]
#define epv_ref(a_1,a_2) face10_1.epv[(a_2)*2 + a_1 - 3]
#define spv_ref(a_1,a_2) spv[(a_2)*3 + a_1]
#define cyep_ref(a_1,a_2) face12_1.cyep[(a_2)*30 + a_1 - 31]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




    /* Parameter adjustments */
    spv -= 4;
    --unvect;
    --pnt;

    /* Function Body */
    *fail = FALSE_;
    *nedge = face12_1.cynep[*icy - 1];
    i__1 = face12_1.cynep[*icy - 1];
    for (ke = 1; ke <= i__1; ++ke) {

/*     vertex number (use first vertex of edge) */

	iep = cyep_ref(ke, *icy);
	iv = epv_ref(1, iep);
	if (iv != 0) {

/*     vector from north pole to vertex */

	    for (k = 1; k <= 3; ++k) {
		polev[k - 1] = v_ref(k, iv) - pnt[k];
	    }

/*     calculate multiplication factor */

	    dt = dot_(polev, &unvect[1]);
	    if (dt == 0.) {
		*fail = TRUE_;
		return 0;
	    }
	    f = face01_1.ar[*ia - 1] * 2 / dt;
	    if (f < 1.) {
		*fail = TRUE_;
		return 0;
	    }

/*     projected vertex for this convex edge */

	    for (k = 1; k <= 3; ++k) {
		spv_ref(k, ke) = pnt[k] + f * polev[k - 1];
	    }
	}
    }
    return 0;
} /* projct_ */

#undef cyep_ref
#undef spv_ref
#undef epv_ref
#undef v_ref




/*     ####################### */
/*     ##                   ## */
/*     ##  function ptincy  ## */
/*     ##                   ## */
/*     ####################### */


logical ptincy_(doublereal *pnt, doublereal *unvect, integer *icy)
{
    /* System generated locals */
    integer i__1;
    logical ret_val;

    /* Local variables */
    static integer k, ic, ke, it, iep;
    extern doublereal dot_(doublereal *, doublereal *);
    static doublereal epu[90]	/* was [3][30] */, spv[90]	/* was [3][30]
	     */;
    static logical fail;
    static integer nedge, iaoth, iatom;
    static doublereal acvect[3];
    extern /* Subroutine */ int epuclc_(doublereal *, integer *, doublereal *)
	    ;
    static doublereal cpvect[3];
    extern doublereal rotang_(doublereal *, integer *, doublereal *);
    static doublereal totang;
    extern /* Subroutine */ int projct_(doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, logical *);


#define a_ref(a_1,a_2) face01_1.a[(a_2)*3 + a_1 - 4]
#define c___ref(a_1,a_2) face09_1.c__[(a_2)*3 + a_1 - 4]
#define ta_ref(a_1,a_2) face05_1.ta[(a_2)*2 + a_1 - 3]
#define cyep_ref(a_1,a_2) face12_1.cyep[(a_2)*30 + a_1 - 31]



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

/*     ########################################################### */
/*     ##                                                       ## */
/*     ##  faces.i  --  variables for Connolly area and volume  ## */
/*     ##                                                       ## */
/*     ########################################################### */


/*     maxcls   maximum number of neighboring atom pairs */
/*     maxtt    maximum number of temporary tori */
/*     maxt     maximum number of total tori */
/*     maxp     maximum number of probe positions */
/*     maxv     maximum number of vertices */
/*     maxen    maximum number of concave edges */
/*     maxfn    maximum number of concave faces */
/*     maxc     maximum number of circles */
/*     maxep    maximum number of convex edges */
/*     maxfs    maximum number of saddle faces */
/*     maxcy    maximum number of cycles */
/*     mxcyep   maximum number of cycle convex edges */
/*     maxfp    maximum number of convex faces */
/*     mxfpcy   maximum number of convex face cycles */




/*     a        atomic coordinates */
/*     ar       atomic radii */
/*     pr       probe radius */
/*     na       number of atoms */




/*     skip     if true, atom is not used */
/*     nosurf   if true, atom has no free surface */
/*     afree    atom free of neighbors */
/*     abur     atom buried */




/*     acls     begin and end pointers for atoms neighbors */
/*     cls      atom numbers of neighbors */
/*     clst     pointer from neighbor to torus */




/*     ntt      number of temporary tori */
/*     tta      temporary torus atom numbers */
/*     ttfe     first edge of each temporary torus */
/*     ttle     last edge of each temporary torus */
/*     enext    pointer to next edge of temporary torus */
/*     ttbur    temporary torus buried */
/*     ttfree   temporary torus free */




/*     t        torus center */
/*     tr       torus radius */
/*     tax      torus axis */
/*     nt       number of tori */
/*     ta       torus atom numbers */
/*     tfe      torus first edge */
/*     tfree    torus free of neighbors */




/*     p        probe position coordinates */
/*     np       number of probe positions */
/*     pa       probe position atom numbers */




/*     v        vertex coordinates */
/*     nv       number of vertices */
/*     va       vertex atom number */
/*     vp       vertex probe number */




/*     nen      number of concave edges */
/*     env      vertex numbers for each concave edge */
/*     nfn      number of concave faces */
/*     fnen     concave face concave edge numbers */




/*     c        circle center */
/*     cr       circle radius */
/*     nc       number of circles */
/*     ca       circle atom number */
/*     ct       circle torus number */




/*     nep      number of convex edges */
/*     epc      convex edge circle number */
/*     epv      convex edge vertex numbers */
/*     afe      first convex edge of each atom */
/*     ale      last convex edge of each atom */
/*     epnext   pointer to next convex edge of atom */




/*     nfs      number of saddle faces */
/*     fsen     saddle face concave edge numbers */
/*     fsep     saddle face convex edge numbers */




/*     ncy      number of cycles */
/*     cynep    number of convex edges in cycle */
/*     cyep     cycle convex edge numbers */




/*     nfp      number of convex faces */
/*     fpa      atom number of convex face */
/*     fpcy     convex face cycle numbers */
/*     fpncy    number of cycles bounding convex face */




/*     check for eaten by neighbor */

    /* Parameter adjustments */
    --unvect;
    --pnt;

    /* Function Body */
    i__1 = face12_1.cynep[*icy - 1];
    for (ke = 1; ke <= i__1; ++ke) {
	iep = cyep_ref(ke, *icy);
	ic = face10_1.epc[iep - 1];
	it = face09_1.ct[ic - 1];
	iatom = face09_1.ca[ic - 1];
	if (ta_ref(1, it) == iatom) {
	    iaoth = ta_ref(2, it);
	} else {
	    iaoth = ta_ref(1, it);
	}
	for (k = 1; k <= 3; ++k) {
	    acvect[k - 1] = a_ref(k, iaoth) - a_ref(k, iatom);
	    cpvect[k - 1] = pnt[k] - c___ref(k, ic);
	}
	if (dot_(acvect, cpvect) >= 0.) {
	    ret_val = FALSE_;
	    return ret_val;
	}
    }
    if (face12_1.cynep[*icy - 1] <= 2) {
	ret_val = TRUE_;
	return ret_val;
    }
    projct_(&pnt[1], &unvect[1], icy, &iatom, spv, &nedge, &fail);
    if (fail) {
	ret_val = TRUE_;
	return ret_val;
    }
    epuclc_(spv, &nedge, epu);
    totang = rotang_(epu, &nedge, &unvect[1]);
    ret_val = totang > 0.;
    return ret_val;
} /* ptincy_ */

#undef cyep_ref
#undef ta_ref
#undef c___ref
#undef a_ref




/*     ######################### */
/*     ##                     ## */
/*     ##  subroutine epuclc  ## */
/*     ##                     ## */
/*     ######################### */


/* Subroutine */ int epuclc_(doublereal *spv, integer *nedge, doublereal *epu)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, ke, le, ke2;
    static doublereal epun;
    extern doublereal anorm_(doublereal *);


#define epu_ref(a_1,a_2) epu[(a_2)*3 + a_1]
#define spv_ref(a_1,a_2) spv[(a_2)*3 + a_1]



/*     calculate unit vectors along edges */

    /* Parameter adjustments */
    epu -= 4;
    spv -= 4;

    /* Function Body */
    i__1 = *nedge;
    for (ke = 1; ke <= i__1; ++ke) {

/*     get index of second edge of corner */

	if (ke < *nedge) {
	    ke2 = ke + 1;
	} else {
	    ke2 = 1;
	}

/*     unit vector along edge of cycle */

	for (k = 1; k <= 3; ++k) {
	    epu_ref(k, ke) = spv_ref(k, ke2) - spv_ref(k, ke);
	}
	epun = anorm_(&epu_ref(1, ke));
/*        if (epun .le. 0.0d0)  call cerror ('Null Edge in Cycle') */

/*     normalize */

	if (epun > 0.) {
	    for (k = 1; k <= 3; ++k) {
		epu_ref(k, ke) = epu_ref(k, ke) / epun;
	    }
	} else {
	    for (k = 1; k <= 3; ++k) {
		epu_ref(k, ke) = 0.;
	    }
	}
    }

/*     vectors for null edges come from following or preceding edges */

    i__1 = *nedge;
    for (ke = 1; ke <= i__1; ++ke) {
	if (anorm_(&epu_ref(1, ke)) <= 0.) {
	    le = ke - 1;
	    if (le <= 0) {
		le = *nedge;
	    }
	    for (k = 1; k <= 3; ++k) {
		epu_ref(k, ke) = epu_ref(k, le);
	    }
	}
    }
    return 0;
} /* epuclc_ */

#undef spv_ref
#undef epu_ref




/*     ####################### */
/*     ##                   ## */
/*     ##  function rotang  ## */
/*     ##                   ## */
/*     ####################### */


doublereal rotang_(doublereal *epu, integer *nedge, doublereal *unvect)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double acos(doublereal);

    /* Local variables */
    static integer ke;
    static doublereal dt, ang;
    extern doublereal dot_(doublereal *, doublereal *);
    static doublereal crs[3], totang;
    extern /* Subroutine */ int vcross_(doublereal *, doublereal *, 
	    doublereal *);


#define epu_ref(a_1,a_2) epu[(a_2)*3 + a_1]



    /* Parameter adjustments */
    epu -= 4;
    --unvect;

    /* Function Body */
    totang = 0.;

/*     sum angles at vertices of cycle */

    i__1 = *nedge;
    for (ke = 1; ke <= i__1; ++ke) {
	if (ke < *nedge) {
	    dt = dot_(&epu_ref(1, ke), &epu_ref(1, ke + 1));
	    vcross_(&epu_ref(1, ke), &epu_ref(1, ke + 1), crs);
	} else {

/*     closing edge of cycle */

	    dt = dot_(&epu_ref(1, ke), &epu_ref(1, 1));
	    vcross_(&epu_ref(1, ke), &epu_ref(1, 1), crs);
	}
	if (dt < -1.) {
	    dt = -1.;
	}
	if (dt > 1.) {
	    dt = 1.;
	}
	ang = acos(dt);
	if (dot_(crs, &unvect[1]) > 0.) {
	    ang = -ang;
	}

/*     add to total for cycle */

	totang += ang;
    }
    ret_val = totang;
    return ret_val;
} /* rotang_ */

#undef epu_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine vcross  --  find cross product of two vectors  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "vcross" finds the cross product of two vectors */


/* Subroutine */ int vcross_(doublereal *x, doublereal *y, doublereal *z__)
{


    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    z__[1] = x[2] * y[3] - x[3] * y[2];
    z__[2] = x[3] * y[1] - x[1] * y[3];
    z__[3] = x[1] * y[2] - x[2] * y[1];
    return 0;
} /* vcross_ */



/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  function dot  --  find the dot product of two vectors  ## */
/*     ##                                                         ## */
/*     ############################################################# */


/*     "dot" finds the dot product of two vectors */


doublereal dot_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    doublereal ret_val;



    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    ret_val = x[1] * y[1] + x[2] * y[2] + x[3] * y[3];
    return ret_val;
} /* dot_ */



/*     ####################################################### */
/*     ##                                                   ## */
/*     ##  function anorm  --  find the length of a vector  ## */
/*     ##                                                   ## */
/*     ####################################################### */


/*     "anorm" finds the norm (length) of a vector; used as a */
/*     service routine by the Connolly surface area and volume */
/*     computation */


doublereal anorm_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);



    /* Parameter adjustments */
    --x;

    /* Function Body */
/* Computing 2nd power */
    d__1 = x[1];
/* Computing 2nd power */
    d__2 = x[2];
/* Computing 2nd power */
    d__3 = x[3];
    ret_val = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    if (ret_val < 0.) {
	ret_val = 0.;
    }
    ret_val = sqrt(ret_val);
    return ret_val;
} /* anorm_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine vnorm  --  normalize a vector to unit length  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "vnorm" normalizes a vector to unit length; used as a */
/*     service routine by the Connolly surface area and volume */
/*     computation */


/* Subroutine */ int vnorm_(doublereal *x, doublereal *xn)
{
    static integer k;
    static doublereal ax;
    extern doublereal anorm_(doublereal *);



    /* Parameter adjustments */
    --xn;
    --x;

    /* Function Body */
    ax = anorm_(&x[1]);
    for (k = 1; k <= 3; ++k) {
	xn[k] = x[k] / ax;
    }
    return 0;
} /* vnorm_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  function dist2  --  distance squared between two points  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "dist2" finds the distance squared between two points; used */
/*     as a service routine by the Connolly surface area and volume */
/*     computation */


doublereal dist2_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;



    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
/* Computing 2nd power */
    d__1 = x[1] - y[1];
/* Computing 2nd power */
    d__2 = x[2] - y[2];
/* Computing 2nd power */
    d__3 = x[3] - y[3];
    ret_val = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return ret_val;
} /* dist2_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  function triple  --  form triple product of three vectors  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "triple" finds the triple product of three vectors; used as */
/*     a service routine by the Connolly surface area and volume */
/*     computation */


doublereal triple_(doublereal *x, doublereal *y, doublereal *z__)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal xy[3];
    extern doublereal dot_(doublereal *, doublereal *);
    extern /* Subroutine */ int vcross_(doublereal *, doublereal *, 
	    doublereal *);



    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    vcross_(&x[1], &y[1], xy);
    ret_val = dot_(xy, &z__[1]);
    return ret_val;
} /* triple_ */



/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  function vecang  --  finds the angle between two vectors  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "vecang" finds the angle between two vectors handed with respect */
/*     to a coordinate axis; returns an angle in the range [0,2*pi] */


doublereal vecang_(doublereal *v1, doublereal *v2, doublereal *axis, 
	doublereal *hand)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double acos(doublereal);

    /* Local variables */
    static doublereal a1, a2, a12, dt;
    extern doublereal dot_(doublereal *, doublereal *);
    static doublereal angle;
    extern doublereal anorm_(doublereal *), triple_(doublereal *, doublereal *
	    , doublereal *);



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




    /* Parameter adjustments */
    --axis;
    --v2;
    --v1;

    /* Function Body */
    a1 = anorm_(&v1[1]);
    a2 = anorm_(&v2[1]);
    dt = dot_(&v1[1], &v2[1]);
    a12 = a1 * a2;
    if (abs(a12) != 0.) {
	dt /= a12;
    }
    if (dt < -1.) {
	dt = -1.;
    }
    if (dt > 1.) {
	dt = 1.;
    }
    angle = acos(dt);
    if (*hand * triple_(&v1[1], &v2[1], &axis[1]) < 0.) {
	ret_val = 6.2831853071795862 - angle;
    } else {
	ret_val = angle;
    }
    return ret_val;
} /* vecang_ */



/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine cirpln  --  locate circle-plane intersection  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "cirpln" determines the points of intersection between a */
/*     specified circle and plane */


/* Subroutine */ int cirpln_(doublereal *circen, doublereal *cirrad, 
	doublereal *cirvec, doublereal *plncen, doublereal *plnvec, logical *
	cinsp, logical *cintp, doublereal *xpnt1, doublereal *xpnt2)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k;
    static doublereal dcp, dir;
    extern doublereal dot_(doublereal *, doublereal *);
    static doublereal pnt1[3], rlen, vect1[3], vect2[3];
    extern doublereal anorm_(doublereal *);
    static doublereal ratio;
    extern /* Subroutine */ int vnorm_(doublereal *, doublereal *);
    static doublereal uvect1[3], uvect2[3], cpvect[3];
    extern /* Subroutine */ int vcross_(doublereal *, doublereal *, 
	    doublereal *);



    /* Parameter adjustments */
    --xpnt2;
    --xpnt1;
    --plnvec;
    --plncen;
    --cirvec;
    --circen;

    /* Function Body */
    for (k = 1; k <= 3; ++k) {
	cpvect[k - 1] = plncen[k] - circen[k];
    }
    dcp = dot_(cpvect, &plnvec[1]);
    *cinsp = dcp > 0.;
    vcross_(&plnvec[1], &cirvec[1], vect1);
    if (anorm_(vect1) > 0.) {
	vnorm_(vect1, uvect1);
	vcross_(&cirvec[1], uvect1, vect2);
	if (anorm_(vect2) > 0.) {
	    vnorm_(vect2, uvect2);
	    dir = dot_(uvect2, &plnvec[1]);
	    if (dir != 0.) {
		ratio = dcp / dir;
		if (abs(ratio) <= *cirrad) {
		    for (k = 1; k <= 3; ++k) {
			pnt1[k - 1] = circen[k] + ratio * uvect2[k - 1];
		    }
/* Computing 2nd power */
		    d__1 = *cirrad;
/* Computing 2nd power */
		    d__2 = ratio;
		    rlen = d__1 * d__1 - d__2 * d__2;
		    if (rlen < 0.) {
			rlen = 0.;
		    }
		    rlen = sqrt(rlen);
		    for (k = 1; k <= 3; ++k) {
			xpnt1[k] = pnt1[k - 1] - rlen * uvect1[k - 1];
			xpnt2[k] = pnt1[k - 1] + rlen * uvect1[k - 1];
		    }
		    *cintp = TRUE_;
		    return 0;
		}
	    }
	}
    }
    *cintp = FALSE_;
    return 0;
} /* cirpln_ */



/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine gendot  --  find surface points on unit sphere  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "gendot" finds the coordinates of a specified number of surface */
/*     points for a sphere with the input radius and coordinate center */


/* Subroutine */ int gendot_(integer *ndots, doublereal *dots, doublereal *
	radius, doublereal *xcenter, doublereal *ycenter, doublereal *zcenter)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal x, y, z__, fi, fj, xy;
    static integer nvert, nequat, nhoriz;


#define dots_ref(a_1,a_2) dots[(a_2)*3 + a_1]



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




    /* Parameter adjustments */
    dots -= 4;

    /* Function Body */
    nequat = (integer) sqrt((doublereal) (*ndots) * 3.141592653589793238);
    nvert = (integer) (nequat * .5);
    if (nvert < 1) {
	nvert = 1;
    }
    k = 0;
    i__1 = nvert;
    for (i__ = 0; i__ <= i__1; ++i__) {
	fi = (doublereal) i__ * 3.141592653589793238 / (doublereal) nvert;
	z__ = cos(fi);
	xy = sin(fi);
	nhoriz = (integer) (nequat * xy);
	if (nhoriz < 1) {
	    nhoriz = 1;
	}
	i__2 = nhoriz - 1;
	for (j = 0; j <= i__2; ++j) {
	    fj = (doublereal) j * 6.2831853071795862 / (doublereal) nhoriz;
	    x = cos(fj) * xy;
	    y = sin(fj) * xy;
	    ++k;
	    dots_ref(1, k) = x * *radius + *xcenter;
	    dots_ref(2, k) = y * *radius + *ycenter;
	    dots_ref(3, k) = z__ * *radius + *zcenter;
	    if (k >= *ndots) {
		goto L10;
	    }
	}
    }
L10:
    *ndots = k;
    return 0;
} /* gendot_ */

#undef dots_ref




/*     ################################################################ */
/*     ##                                                            ## */
/*     ##  subroutine cerror  --  surface area-volume error message  ## */
/*     ##                                                            ## */
/*     ################################################################ */


/*     "cerror" is the error handling routine for the Connolly */
/*     surface area and volume computation */


/* Subroutine */ int cerror_(char *string, ftnlen string_len)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 CONNOLLY  --  \002,a)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern integer trimtext_(char *, ftnlen);
    static integer leng;
    extern /* Subroutine */ int fatal_(void);

    /* Fortran I/O blocks */
    static cilist io___489 = { 0, 0, 0, fmt_10, 0 };




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




/*     write out the error message and quit */

    leng = trimtext_(string, string_len);
    io___489.ciunit = iounit_1.iout;
    s_wsfe(&io___489);
    do_fio(&c__1, string, leng);
    e_wsfe();
    fatal_();
    return 0;
} /* cerror_ */

