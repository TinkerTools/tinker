#
#
#  #############################################################
#  ##                                                         ##
#  ##  linkapbs.make  --  link the TINKER programs with APBS  ##
#  ##            (Intel Fortran for MacOS Version)            ##
#  ##                                                         ##
#  #############################################################
#
#
#  The libraries in APBSLIBS are generated during the standard
#  build procedure for the APBS software package, and should be
#  moved to a directory compatible with the APBSLIB definition
#
#
setenv APBSLIBS "../apbs/macosx/lib/libapbsmainroutines.a ../apbs/macosx/lib/libapbs.a ../apbs/macosx/lib/libmaloc.a ../apbs/macosx/lib/libapbsblas.a"
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o alchemy.x alchemy.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip alchemy.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o analyze.x analyze.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip analyze.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o anneal.x anneal.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip anneal.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o archive.x archive.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip archive.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o bar.x bar.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip bar.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o correlate.x correlate.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip correlate.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o crystal.x crystal.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip crystal.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o diffuse.x diffuse.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip diffuse.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o distgeom.x distgeom.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip distgeom.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o document.x document.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip document.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o dynamic.x dynamic.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip dynamic.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o gda.x gda.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip gda.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o intedit.x intedit.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip intedit.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o intxyz.x intxyz.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip intxyz.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o minimize.x minimize.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip minimize.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o minirot.x minirot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip minirot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o minrigid.x minrigid.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip minrigid.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o molxyz.x molxyz.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip molxyz.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o monte.x monte.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip monte.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o newton.x newton.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip newton.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o newtrot.x newtrot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip newtrot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o nucleic.x nucleic.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip nucleic.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o optimize.x optimize.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip optimize.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o optirot.x optirot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip optirot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o optrigid.x optrigid.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip optrigid.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o path.x path.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip path.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o pdbxyz.x pdbxyz.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip pdbxyz.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o polarize.x polarize.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip polarize.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o poledit.x poledit.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip poledit.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o potential.x potential.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip potential.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o prmedit.x prmedit.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip prmedit.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o protein.x protein.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip protein.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o pss.x pss.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip pss.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o pssrigid.x pssrigid.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip pssrigid.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o pssrot.x pssrot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip pssrot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o radial.x radial.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip radial.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o saddle.x saddle.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip saddle.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o scan.x scan.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip scan.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o sniffer.x sniffer.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip sniffer.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o spacefill.x spacefill.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip spacefill.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o spectrum.x spectrum.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip spectrum.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o superpose.x superpose.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip superpose.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o sybylxyz.x sybylxyz.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip sybylxyz.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o testgrad.x testgrad.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip testgrad.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o testhess.x testhess.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip testhess.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o testpair.x testpair.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip testpair.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o testpol.x testpol.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip testpol.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o testrot.x testrot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip testrot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o timer.x timer.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip timer.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o timerot.x timerot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip timerot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o torsfit.x torsfit.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip torsfit.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o valence.x valence.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip valence.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o vibbig.x vibbig.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip vibbig.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o vibrate.x vibrate.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip vibrate.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o vibrot.x vibrot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip vibrot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o xtalfit.x xtalfit.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip xtalfit.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o xtalmin.x xtalmin.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip xtalmin.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o xyzedit.x xyzedit.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip xyzedit.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o xyzint.x xyzint.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip xyzint.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o xyzpdb.x xyzpdb.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip xyzpdb.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o xyzsybyl.x xyzsybyl.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip xyzsybyl.x
