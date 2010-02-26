/* pmpb.f -- translated by f2c (version 20050501).
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



/*     ################################################################ */
/*     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  routines below implement dummy versions of the APBS   ## */
/*     ##  calls required for TINKER to interface with the APBS  ## */
/*     ##  Poisson-Boltzmann solver package from Nathan Baker    ## */
/*     ##                                                        ## */
/*     ############################################################ */

/*     ############################## */
/*     ##                          ## */
/*     ##  subroutine apbsinitial  ## */
/*     ##                          ## */
/*     ############################## */


/* Subroutine */ int apbsinitial_(integer *dime, doublereal *grid, doublereal 
	*gcent, doublereal *cgrid, doublereal *cgcent, doublereal *fgrid, 
	doublereal *fgcent, doublereal *pdie, doublereal *sdie, doublereal *
	srad, doublereal *swin, doublereal *sdens, doublereal *kelvin, 
	integer *ionn, doublereal *ionc, integer *ionq, doublereal *ionr, 
	char *pbtyp, integer *pbtyplen, char *pbsoln, integer *pbsolnlen, 
	char *bcfl, integer *bcfllen, char *chgm, integer *chgmlen, char *
	srfm, integer *srfmlen, ftnlen pbtyp_len, ftnlen pbsoln_len, ftnlen 
	bcfl_len, ftnlen chgm_len, ftnlen srfm_len)
{
    /* Parameter adjustments */
    --ionr;
    --ionq;
    --ionc;
    --fgcent;
    --fgrid;
    --cgcent;
    --cgrid;
    --gcent;
    --grid;
    --dime;

    /* Function Body */
    return 0;
} /* apbsinitial_ */



/*     ############################# */
/*     ##                         ## */
/*     ##  subroutine apbsempole  ## */
/*     ##                         ## */
/*     ############################# */


/* Subroutine */ int apbsempole_(integer *n, doublereal *pos, doublereal *
	rsolv, doublereal *pbpole, doublereal *pbe, doublereal *apbe, 
	doublereal *pbep, doublereal *pbfp, doublereal *pbtp)
{
    /* Parameter adjustments */
    --pbtp;
    --pbfp;
    --pbep;
    --apbe;
    --pbpole;
    --rsolv;
    --pos;

    /* Function Body */
    return 0;
} /* apbsempole_ */



/*     ############################# */
/*     ##                         ## */
/*     ##  subroutine apbsinduce  ## */
/*     ##                         ## */
/*     ############################# */


/* Subroutine */ int apbsinduce_(doublereal *indpole, doublereal *pbeuind)
{
    /* Parameter adjustments */
    --pbeuind;
    --indpole;

    /* Function Body */
    return 0;
} /* apbsinduce_ */



/*     ############################### */
/*     ##                           ## */
/*     ##  subroutine apbsnlinduce  ## */
/*     ##                           ## */
/*     ############################### */


/* Subroutine */ int apbsnlinduce_(doublereal *inppole, doublereal *pbeuinp)
{
    /* Parameter adjustments */
    --pbeuinp;
    --inppole;

    /* Function Body */
    return 0;
} /* apbsnlinduce_ */



/*     ################################### */
/*     ##                               ## */
/*     ##  subroutine pbdirectpolforce  ## */
/*     ##                               ## */
/*     ################################### */


/* Subroutine */ int pbdirectpolforce_(doublereal *indpole, doublereal *
	inppole, doublereal *directf, doublereal *directt)
{
    /* Parameter adjustments */
    --directt;
    --directf;
    --inppole;
    --indpole;

    /* Function Body */
    return 0;
} /* pbdirectpolforce_ */



/*     ################################### */
/*     ##                               ## */
/*     ##  subroutine pbmutualpolforce  ## */
/*     ##                               ## */
/*     ################################### */


/* Subroutine */ int pbmutualpolforce_(doublereal *indpole, doublereal *
	inppole, doublereal *mutualf)
{
    /* Parameter adjustments */
    --mutualf;
    --inppole;
    --indpole;

    /* Function Body */
    return 0;
} /* pbmutualpolforce_ */



/*     ############################ */
/*     ##                        ## */
/*     ##  subroutine apbsfinal  ## */
/*     ##                        ## */
/*     ############################ */


/* Subroutine */ int apbsfinal_(void)
{
    return 0;
} /* apbsfinal_ */

