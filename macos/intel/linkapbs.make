#
#
#  #############################################################
#  ##                                                         ##
#  ##  linkapbs.make  --  link the Tinker programs with APBS  ##
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
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o alchemy.x alchemy.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip alchemy.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o analyze.x analyze.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip analyze.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o anneal.x anneal.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip anneal.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o archive.x archive.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip archive.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o bar.x bar.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip bar.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o correlate.x correlate.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip correlate.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o critical.x critical.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip critical.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o crystal.x crystal.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip crystal.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o diffuse.x diffuse.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip diffuse.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o distgeom.x distgeom.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip distgeom.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o document.x document.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip document.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o dynamic.x dynamic.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip dynamic.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o freefix.x freefix.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip freefix.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o gda.x gda.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip gda.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o intedit.x intedit.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip intedit.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o intxyz.x intxyz.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip intxyz.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o minimize.x minimize.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip minimize.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o minirot.x minirot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip minirot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o minrigid.x minrigid.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip minrigid.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o mol2xyz.x mol2xyz.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip mol2xyz.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o molxyz.x molxyz.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip molxyz.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o monte.x monte.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip monte.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o newton.x newton.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip newton.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o newtrot.x newtrot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip newtrot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o nucleic.x nucleic.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip nucleic.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o optimize.x optimize.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip optimize.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o optirot.x optirot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip optirot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o optrigid.x optrigid.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip optrigid.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o path.x path.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip path.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o pdbxyz.x pdbxyz.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip pdbxyz.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o polarize.x polarize.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip polarize.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o poledit.x poledit.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip poledit.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o potential.x potential.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip potential.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o prmedit.x prmedit.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip prmedit.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o protein.x protein.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip protein.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o pss.x pss.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip pss.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o pssrigid.x pssrigid.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip pssrigid.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o pssrot.x pssrot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip pssrot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o radial.x radial.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip radial.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o saddle.x saddle.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip saddle.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o scan.x scan.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip scan.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o sniffer.x sniffer.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip sniffer.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o spacefill.x spacefill.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip spacefill.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o spectrum.x spectrum.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip spectrum.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o superpose.x superpose.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip superpose.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o testgrad.x testgrad.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip testgrad.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o testhess.x testhess.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip testhess.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o testpair.x testpair.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip testpair.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o testpol.x testpol.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip testpol.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o testrot.x testrot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip testrot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o testvir.x testvir.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip testvir.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o timer.x timer.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip timer.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o timerot.x timerot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip timerot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o torsfit.x torsfit.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip torsfit.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o valence.x valence.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip valence.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o vibbig.x vibbig.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip vibbig.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o vibrate.x vibrate.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip vibrate.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o vibrot.x vibrot.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip vibrot.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o xtalfit.x xtalfit.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip xtalfit.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o xtalmin.x xtalmin.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip xtalmin.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o xyzedit.x xyzedit.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip xyzedit.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o xyzint.x xyzint.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip xyzint.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o xyzmol2.x xyzmol2.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip xyzmol2.x
ifort -O3 -no-ipo -no-prec-div -qopenmp -static-intel -o xyzpdb.x xyzpdb.o libtinker.a libfftw3_threads.a libfftw3.a $APBSLIBS ; strip xyzpdb.x
