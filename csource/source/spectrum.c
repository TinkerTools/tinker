/* spectrum.f -- translated by f2c (version 20050501).
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
    integer nprior, ldir, leng;
    char filename[120], outfile[120];
} files_;

#define files_1 files_

struct {
    integer iout, input;
} iounit_;

#define iounit_1 iounit_

/* Table of constant values */

static integer c__1 = 1;
static integer c__5 = 5;
static integer c__3 = 3;



/*     ################################################### */
/*     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  program spectrum  --  power spectrum from autocorrelation  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "spectrum" computes a power spectrum over a wavelength range */
/*     from the velocity autocorrelation as a function of time */


/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_10[] = "(/,\002 Enter Name of Velocity Autocorrelatio"
	    "n\002,\002 File :  \002,$)";
    static char fmt_20[] = "(a120)";
    static char fmt_40[] = "(/,\002 Enter Time Step for Autocorrelation Dat"
	    "a\002,\002 [1.0 fs] :  \002,$)";
    static char fmt_50[] = "(f20.0)";
    static char fmt_60[] = "(a120)";
    static char fmt_70[] = "()";
    static char fmt_80[] = "(a120)";
    static char fmt_100[] = "(/,\002 Power Spectrum from Velocity Autocorrel"
	    "ation :\002,//,4x,\002Frequency (cm-1)\002,10x,\002Intensity\002"
	    ",/)";
    static char fmt_110[] = "(3x,f12.2,8x,f16.6)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;
    olist o__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *), s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *)
	    , do_fio(integer *, char *, ftnlen), e_rsfe(void), s_rsli(icilist 
	    *), do_lio(integer *, integer *, char *, ftnlen), e_rsli(void), 
	    f_open(olist *), f_rew(alist *), s_cmp(char *, char *, ftnlen, 
	    ftnlen);
    double cos(doublereal);

    /* Local variables */
    extern /* Subroutine */ int basefile_(char *, ftnlen);
    extern integer freeunit_(void);
    static integer i__, k;
    static doublereal vel[100001];
    static logical done;
    static doublereal aver, freq;
    static integer ivel;
    static doublereal time;
    static integer nvel;
    static doublereal norm, step;
    static integer nsamp;
    static logical exist;
    static doublereal factor;
    static char record[120];
    extern /* Subroutine */ int suffix_(char *, char *, ftnlen, ftnlen);
    static char string[120], velfile[120];
    extern /* Subroutine */ int initial_(void);
    static doublereal intense[5000];
    extern /* Subroutine */ int nextarg_(char *, logical *, ftnlen), version_(
	    char *, char *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_20, 0 };
    static icilist io___7 = { 1, string, 1, 0, 120, 1 };
    static cilist io___8 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___16 = { 1, 0, 1, fmt_80, 0 };
    static icilist io___17 = { 0, record, 0, 0, 120, 1 };
    static cilist io___28 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_110, 0 };




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




/*     ################################################### */
/*     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ############################################################ */
/*     ##                                                        ## */
/*     ##  units.i  --  physical constants and unit conversions  ## */
/*     ##                                                        ## */
/*     ############################################################ */


/*     literature reference: */

/*     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended */
/*     Values of the Fundamental Physical Constants: 2006", Reviews of */
/*     Modern Physics, 80, 633-730 (2008) */

/*     The "2006 CODATA Recommended Values" are also available from */
/*     the NIST Reference on Constants, Units, and Uncertainty site */
/*     at http://physics.nist.gov/cuu/index.html */

/*     Most values below are derived from the 2006 CODATA reference */

/*     The conversion from calorie to Joule is the definition of the */
/*     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992) */

/*     The "coulomb" energy conversion factor is found by dimensional */
/*     analysis of Coulomb's Law, ie, by dividing the square of the */
/*     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is */
/*     the permittivity of vacuum (the "electric constant"); note that */
/*     eps0 is typically given in F/m, equivalent to C**2/(J-m) */

/*     The approximate value used for the Debye, 3.33564 x 10-30 C-m, */
/*     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997) */

/*     The value of "prescon" is based on definition of 1 atmosphere */
/*     as 101325 Pa set by the 10th Conference Generale des Poids et */
/*     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3 */

/*     avogadro    Avogadro's number (N) in particles/mole */
/*     lightspd    speed of light in vacuum (c) in cm/ps */
/*     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K */
/*     gasconst    ideal gas constant (R) in kcal/mole/K */
/*     emass       mass of an electron in atomic mass units */
/*     joule       conversion from calories to joules */
/*     convert     conversion from kcal to g*Ang**2/ps**2 */
/*     bohr        conversion from Bohrs to Angstroms */
/*     hartree     conversion from Hartree to kcal/mole */
/*     evolt       conversion from Hartree to electron-volts */
/*     efreq       conversion from Hartree to cm-1 */
/*     coulomb     conversion from electron**2/Ang to kcal/mole */
/*     debye       conversion from electron-Ang to Debyes */
/*     prescon     conversion from kcal/mole/Ang**3 to Atm */




/*     perform the standard initialization functions */

    initial_();

/*     try to get a filename from the command line arguments */

    nextarg_(velfile, &exist, (ftnlen)120);
    if (exist) {
	basefile_(velfile, (ftnlen)120);
	suffix_(velfile, "vel", (ftnlen)120, (ftnlen)3);
	version_(velfile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = velfile;
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
    }

/*     ask for the velocity autocorrelation data filename */

    while(! exist) {
	io___3.ciunit = iounit_1.iout;
	s_wsfe(&io___3);
	e_wsfe();
	io___4.ciunit = iounit_1.input;
	s_rsfe(&io___4);
	do_fio(&c__1, velfile, (ftnlen)120);
	e_rsfe();
	basefile_(velfile, (ftnlen)120);
	suffix_(velfile, "vel", (ftnlen)120, (ftnlen)3);
	version_(velfile, "old", (ftnlen)120, (ftnlen)3);
	ioin__1.inerr = 0;
	ioin__1.infilen = 120;
	ioin__1.infile = velfile;
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
    }

/*     get the time step between autocorrelation data points */

    step = -1.;
    nextarg_(string, &exist, (ftnlen)120);
    if (exist) {
	i__1 = s_rsli(&io___7);
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = do_lio(&c__5, &c__1, (char *)&step, (ftnlen)sizeof(doublereal))
		;
	if (i__1 != 0) {
	    goto L30;
	}
	i__1 = e_rsli();
	if (i__1 != 0) {
	    goto L30;
	}
    }
L30:
    if (step <= 0.) {
	io___8.ciunit = iounit_1.iout;
	s_wsfe(&io___8);
	e_wsfe();
	io___9.ciunit = iounit_1.input;
	s_rsfe(&io___9);
	do_fio(&c__1, (char *)&step, (ftnlen)sizeof(doublereal));
	e_rsfe();
    }
    if (step <= 0.) {
	step = 1.;
    }
    step *= .001;

/*     open the velocity autocorrelation data file */

    ivel = freeunit_();
    o__1.oerr = 0;
    o__1.ounit = ivel;
    o__1.ofnmlen = 120;
    o__1.ofnm = velfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    al__1.aerr = 0;
    al__1.aunit = ivel;
    f_rew(&al__1);

/*     read through file headers to the start of the data */

    done = FALSE_;
    while(! done) {
	io___12.ciunit = ivel;
	s_rsfe(&io___12);
	do_fio(&c__1, record, (ftnlen)120);
	e_rsfe();
	if (s_cmp(record + 3, "Separation", (ftnlen)10, (ftnlen)10) == 0) {
	    done = TRUE_;
	    io___14.ciunit = ivel;
	    s_rsfe(&io___14);
	    e_rsfe();
	}
    }

/*     read the velocity autocorrelation as a function of time */

    for (i__ = 1; i__ <= 100000; ++i__) {
	io___16.ciunit = ivel;
	i__1 = s_rsfe(&io___16);
	if (i__1 != 0) {
	    goto L90;
	}
	i__1 = do_fio(&c__1, record, (ftnlen)120);
	if (i__1 != 0) {
	    goto L90;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L90;
	}
	s_rsli(&io___17);
	do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&nsamp, (ftnlen)sizeof(integer));
	do_lio(&c__5, &c__1, (char *)&aver, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&norm, (ftnlen)sizeof(doublereal));
	e_rsli();
	nvel = k;
	vel[k] = norm;
    }
L90:

/*     compute the power spectrum via discrete Fourier transform */

    factor = .1883651567308853;
    for (i__ = 1; i__ <= 5000; ++i__) {
	freq = factor * (doublereal) i__;
	intense[i__ - 1] = 0.;
	i__1 = nvel;
	for (k = 0; k <= i__1; ++k) {
	    time = step * (doublereal) k;
	    intense[i__ - 1] += vel[k] * cos(freq * time);
	}
	intense[i__ - 1] = step * 1e3 * intense[i__ - 1];
    }

/*     print the power spectrum intensity at each wavelength */

    io___28.ciunit = iounit_1.iout;
    s_wsfe(&io___28);
    e_wsfe();
    for (i__ = 1; i__ <= 5000; ++i__) {
	io___29.ciunit = iounit_1.iout;
	s_wsfe(&io___29);
	d__1 = (doublereal) i__;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&intense[i__ - 1], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* MAIN__ */

/* Main program alias */ int spectrum_ () { MAIN__ (); return 0; }
