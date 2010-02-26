/* server.f -- translated by f2c (version 20050501).
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
/*     ##  COPYRIGHT (C) 2003 by Michael Schnieders & Jay W. Ponder  ## */
/*     ##                     All Rights Reserved                    ## */
/*     ################################################################ */

/*     ############################################################# */
/*     ##                                                         ## */
/*     ##  routines below implement dummy versions of the socket  ## */
/*     ##  communication calls required for the transmission of   ## */
/*     ##  information between TINKER and Force Field Explorer;   ## */
/*     ##  functional C code is in "server.c", while the dummy    ## */
/*     ##  calls in this file are written in standard Fortran     ## */
/*     ##                                                         ## */
/*     ############################################################# */

/*     ############################ */
/*     ##                        ## */
/*     ##  subroutine chksocket  ## */
/*     ##                        ## */
/*     ############################ */


/* Subroutine */ int chksocket_(integer *flag__)
{


/*     set flag that will disable socket communications */

    *flag__ = 0;
    return 0;
} /* chksocket_ */



/*     ############################ */
/*     ##                        ## */
/*     ##  subroutine createjvm  ## */
/*     ##                        ## */
/*     ############################ */


/* Subroutine */ int createjvm_(integer *flag__)
{
    return 0;
} /* createjvm_ */



/*     ############################# */
/*     ##                         ## */
/*     ##  subroutine destroyjvm  ## */
/*     ##                         ## */
/*     ############################# */


/* Subroutine */ int destroyjvm_(void)
{
    return 0;
} /* destroyjvm_ */



/*     ############################### */
/*     ##                           ## */
/*     ##  subroutine createserver  ## */
/*     ##                           ## */
/*     ############################### */


/* Subroutine */ int createserver_(integer *flag__)
{
    return 0;
} /* createserver_ */



/*     ################################ */
/*     ##                            ## */
/*     ##  subroutine destroyserver  ## */
/*     ##                            ## */
/*     ################################ */


/* Subroutine */ int destroyserver_(void)
{
    return 0;
} /* destroyserver_ */



/*     ############################### */
/*     ##                           ## */
/*     ##  subroutine createsystem  ## */
/*     ##                           ## */
/*     ############################### */


/* Subroutine */ int createsystem_(integer *n, integer *nkey, integer *flag__)
{
    return 0;
} /* createsystem_ */



/*     ############################# */
/*     ##                         ## */
/*     ##  subroutine getmonitor  ## */
/*     ##                         ## */
/*     ############################# */


/* Subroutine */ int getmonitor_(void)
{
    return 0;
} /* getmonitor_ */



/*     ################################# */
/*     ##                             ## */
/*     ##  subroutine releasemonitor  ## */
/*     ##                             ## */
/*     ################################# */


/* Subroutine */ int releasemonitor_(void)
{
    return 0;
} /* releasemonitor_ */



/*     ############################### */
/*     ##                           ## */
/*     ##  subroutine createupdate  ## */
/*     ##                           ## */
/*     ############################### */


/* Subroutine */ int createupdate_(integer *n, integer *mode, integer *amoeba,
	 integer *flag__)
{
    return 0;
} /* createupdate_ */



/*     ############################# */
/*     ##                         ## */
/*     ##  subroutine needupdate  ## */
/*     ##                         ## */
/*     ############################# */


/* Subroutine */ int needupdate_(integer *flag__)
{
    return 0;
} /* needupdate_ */



/*     ############################# */
/*     ##                         ## */
/*     ##  subroutine setupdated  ## */
/*     ##                         ## */
/*     ############################# */


/* Subroutine */ int setupdated_(void)
{
    return 0;
} /* setupdated_ */



/*     ########################## */
/*     ##                      ## */
/*     ##  subroutine setfile  ## */
/*     ##                      ## */
/*     ########################## */


/* Subroutine */ int setfile_(char *filename, ftnlen filename_len)
{
    return 0;
} /* setfile_ */



/*     ################################ */
/*     ##                            ## */
/*     ##  subroutine setforcefield  ## */
/*     ##                            ## */
/*     ################################ */


/* Subroutine */ int setforcefield_(char *forcefield, ftnlen forcefield_len)
{
    return 0;
} /* setforcefield_ */



/*     ############################# */
/*     ##                         ## */
/*     ##  subroutine setkeyword  ## */
/*     ##                         ## */
/*     ############################# */


/* Subroutine */ int setkeyword_(integer *i__, char *keyline, ftnlen 
	keyline_len)
{
    return 0;
} /* setkeyword_ */



/*     ############################### */
/*     ##                           ## */
/*     ##  subroutine setatomtypes  ## */
/*     ##                           ## */
/*     ############################### */


/* Subroutine */ int setatomtypes_(integer *n, integer *type__)
{
    /* Parameter adjustments */
    --type__;

    /* Function Body */
    return 0;
} /* setatomtypes_ */



/*     ############################ */
/*     ##                        ## */
/*     ##  subroutine setatomic  ## */
/*     ##                        ## */
/*     ############################ */


/* Subroutine */ int setatomic_(integer *n, integer *atomic)
{
    /* Parameter adjustments */
    --atomic;

    /* Function Body */
    return 0;
} /* setatomic_ */



/*     ########################## */
/*     ##                      ## */
/*     ##  subroutine setmass  ## */
/*     ##                      ## */
/*     ########################## */


/* Subroutine */ int setmass_(integer *n, doublereal *mass)
{
    /* Parameter adjustments */
    --mass;

    /* Function Body */
    return 0;
} /* setmass_ */



/*     ############################ */
/*     ##                        ## */
/*     ##  subroutine setcharge  ## */
/*     ##                        ## */
/*     ############################ */


/* Subroutine */ int setcharge_(integer *n, doublereal *charge)
{
    /* Parameter adjustments */
    --charge;

    /* Function Body */
    return 0;
} /* setcharge_ */



/*     ################################## */
/*     ##                              ## */
/*     ##  subroutine setconnectivity  ## */
/*     ##                              ## */
/*     ################################## */


/* Subroutine */ int setconnectivity_(integer *n, integer *b1, integer *b2, 
	integer *b3, integer *b4)
{
    /* Parameter adjustments */
    --b4;
    --b3;
    --b2;
    --b1;

    /* Function Body */
    return 0;
} /* setconnectivity_ */



/*     ########################## */
/*     ##                      ## */
/*     ##  subroutine setname  ## */
/*     ##                      ## */
/*     ########################## */


/* Subroutine */ int setname_(integer *i__, char *name__, ftnlen name_len)
{
    return 0;
} /* setname_ */



/*     ########################### */
/*     ##                       ## */
/*     ##  subroutine setstory  ## */
/*     ##                       ## */
/*     ########################### */


/* Subroutine */ int setstory_(integer *i__, char *story, ftnlen story_len)
{
    return 0;
} /* setstory_ */



/*     ################################# */
/*     ##                             ## */
/*     ##  subroutine setcoordinates  ## */
/*     ##                             ## */
/*     ################################# */


/* Subroutine */ int setcoordinates_(integer *n, doublereal *x, doublereal *y,
	 doublereal *z__)
{
    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    return 0;
} /* setcoordinates_ */



/*     ########################## */
/*     ##                      ## */
/*     ##  subroutine setstep  ## */
/*     ##                      ## */
/*     ########################## */


/* Subroutine */ int setstep_(integer *ncycle)
{
    return 0;
} /* setstep_ */



/*     ########################## */
/*     ##                      ## */
/*     ##  subroutine settime  ## */
/*     ##                      ## */
/*     ########################## */


/* Subroutine */ int settime_(doublereal *time)
{
    return 0;
} /* settime_ */



/*     ############################ */
/*     ##                        ## */
/*     ##  subroutine setenergy  ## */
/*     ##                        ## */
/*     ############################ */


/* Subroutine */ int setenergy_(doublereal *energy)
{
    return 0;
} /* setenergy_ */



/*     ############################### */
/*     ##                           ## */
/*     ##  subroutine setgradients  ## */
/*     ##                           ## */
/*     ############################### */


/* Subroutine */ int setgradients_(integer *n, doublereal *dx, doublereal *dy,
	 doublereal *dz)
{
    /* Parameter adjustments */
    --dz;
    --dy;
    --dx;

    /* Function Body */
    return 0;
} /* setgradients_ */



/*     ############################## */
/*     ##                          ## */
/*     ##  subroutine setvelocity  ## */
/*     ##                          ## */
/*     ############################## */


/* Subroutine */ int setvelocity_(integer *n, doublereal *vx, doublereal *vy, 
	doublereal *vz)
{
    /* Parameter adjustments */
    --vz;
    --vy;
    --vx;

    /* Function Body */
    return 0;
} /* setvelocity_ */



/*     ################################## */
/*     ##                              ## */
/*     ##  subroutine setacceleration  ## */
/*     ##                              ## */
/*     ################################## */


/* Subroutine */ int setacceleration_(integer *n, doublereal *ax, doublereal *
	ay, doublereal *az)
{
    /* Parameter adjustments */
    --az;
    --ay;
    --ax;

    /* Function Body */
    return 0;
} /* setacceleration_ */



/*     ############################# */
/*     ##                         ## */
/*     ##  subroutine setinduced  ## */
/*     ##                         ## */
/*     ############################# */


/* Subroutine */ int setinduced_(integer *n, doublereal *ux, doublereal *uy, 
	doublereal *uz)
{
    /* Parameter adjustments */
    --uz;
    --uy;
    --ux;

    /* Function Body */
    return 0;
} /* setinduced_ */

