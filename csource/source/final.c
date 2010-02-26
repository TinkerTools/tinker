/* final.f -- translated by f2c (version 20050501).
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
    integer digits, iprint, iwrite, isend;
    logical verbose, debug, holdup, abort;
} inform_;

#define inform_1 inform_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

struct {
    integer runtyp, cstep;
    doublereal cdt, cenergy, cdx[25000], cdy[25000], cdz[25000];
    logical use_socket__, skt_init__, skt_close__;
} socket_;

#define socket_1 socket_

struct {
    doublereal rsolv[25000], asolv[25000], rborn[25000], drb[25000], drbp[
	    25000], drobc[25000], doffset, p1, p2, p3, p4, p5, gpol[25000], 
	    shct[25000], aobc[25000], bobc[25000], gobc[25000], vsolv[25000], 
	    wace[1000000]	/* was [1000][1000] */, s2ace[1000000]	/* 
	    was [1000][1000] */, uace[1000000]	/* was [1000][1000] */;
    char solvtyp[8], borntyp[8];
} solute_;

#define solute_1 solute_



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  subroutine final  --  final actions before program exit  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     "final" performs any final program actions, prints a status */
/*     message, and then pauses if necessary to avoid closing the */
/*     execution window */


/* Subroutine */ int final_(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 TINKER is Exiting following Normal Termi"
	    "nation\002,\002 of the Program\002,/)";
    static char fmt_20[] = "()";

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_wsfe(cilist *), e_wsfe(
	    void), s_rsfe(cilist *), e_rsfe(void);

    /* Local variables */
    extern /* Subroutine */ int apbsfinal_(void), sktkill_(void);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___2 = { 0, 0, 0, fmt_20, 0 };




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
/*     ##  COPYRIGHT (C)  2002  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  socket.i  --  control parameters for socket communication  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     runtyp      calculation type for passing socket information */
/*     cstep       current optimization or dynamics step number */
/*     cdt         current dynamics cumulative simulation time */
/*     cenergy     current potential energy from simulation */
/*     cdx         current gradient components along the x-axis */
/*     cdy         current gradient components along the y-axis */
/*     cdz         current gradient components along the z-axis */
/*     use_socket  logical flag governing use of external sockets */
/*     skt_init    logical flag to indicate socket initialization */
/*     skt_close   logical flag to indicate socket shutdown */




/*     close any open socket used for external communication */



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################### */
/*     ##                                                           ## */
/*     ##  solute.i  --  parameters for continuum solvation models  ## */
/*     ##                                                           ## */
/*     ############################################################### */


/*     rsolv     atomic radius of each atom for continuum solvation */
/*     asolv     atomic surface area solvation parameters */
/*     rborn     Born radius of each atom for GB/SA solvation */
/*     drb       solvation derivatives with respect to Born radii */
/*     drbp      GK polarization derivatives with respect to Born radii */
/*     drobc     chain rule term for Onufriev-Bashford-Case radii */
/*     doffset   dielectric offset to continuum solvation atomic radii */
/*     p1        single-atom scale factor for analytical Still radii */
/*     p2        1-2 interaction scale factor for analytical Still radii */
/*     p3        1-3 interaction scale factor for analytical Still radii */
/*     p4        nonbonded scale factor for analytical Still radii */
/*     p5        soft cutoff parameter for analytical Still radii */
/*     gpol      polarization self-energy values for each atom */
/*     shct      overlap scale factors for Hawkins-Cramer-Truhlar radii */
/*     aobc      alpha values for Onufriev-Bashford-Case radii */
/*     bobc      beta values for Onufriev-Bashford-Case radii */
/*     gobc      gamma values for Onufriev-Bashford-Case radii */
/*     vsolv     atomic volume of each atom for use with ACE */
/*     wace      "omega" values for atom class pairs for use with ACE */
/*     s2ace     "sigma^2" values for atom class pairs for use with ACE */
/*     uace      "mu" values for atom class pairs for use with ACE */
/*     solvtyp   type of continuum solvation energy model in use */
/*     borntyp   method to be used for the Born radius computation */


    if (socket_1.use_socket__) {
	sktkill_();
    }

/*     free memory used by the APBS Poisson-Boltzmann solver */

    if (s_cmp(solute_1.solvtyp, "PB", (ftnlen)8, (ftnlen)2) == 0) {
	apbsfinal_();
    }

/*     print a final status message before exiting TINKER */

    if (inform_1.debug) {
	io___1.ciunit = iounit_1.iout;
	s_wsfe(&io___1);
	e_wsfe();
    }

/*     may need a pause to avoid closing the execution window */

    if (inform_1.holdup) {
	io___2.ciunit = iounit_1.input;
	s_rsfe(&io___2);
	e_rsfe();
    }
    return 0;
} /* final_ */

