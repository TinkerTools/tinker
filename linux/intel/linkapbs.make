#
#
#  #############################################################
#  ##                                                         ##
#  ##  linkapbs.make  --  link the Tinker programs with APBS  ##
#  ##            (Intel Fortran for Linux Version)            ##
#  ##                                                         ##
#  #############################################################
#
#
#  The libraries in APBSLIBS are generated during the standard
#  build procedure for the APBS software package, and should be
#  moved to a directory compatible with the APBSLIB definition
#
#
setenv APBSLIBS "-L../apbs/lib/linux -lapbsmainroutines -lapbs -lmaloc -lapbsblas"
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o alchemy.x alchemy.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip alchemy.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o analyze.x analyze.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip analyze.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o anneal.x anneal.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip anneal.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o arcedit.x arcedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip arcedit.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o bar.x bar.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip bar.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o correlate.x correlate.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip correlate.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o critical.x critical.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip critical.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o crystal.x crystal.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip crystal.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o diffuse.x diffuse.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip diffuse.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o distgeom.x distgeom.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip distgeom.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o document.x document.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip document.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o dynamic.x dynamic.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip dynamic.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o freefix.x freefix.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip freefix.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o gda.x gda.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip gda.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o intedit.x intedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip intedit.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o intxyz.x intxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip intxyz.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o minimize.x minimize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip minimize.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o minirot.x minirot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip minirot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o minrigid.x minrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip minrigid.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o mol2xyz.x mol2xyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip mol2xyz.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o molxyz.x molxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip molxyz.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o monte.x monte.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip monte.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o newton.x newton.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip newton.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o newtrot.x newtrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip newtrot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o nucleic.x nucleic.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip nucleic.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o optimize.x optimize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip optimize.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o optirot.x optirot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip optirot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o optrigid.x optrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip optrigid.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o path.x path.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip path.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o pdbxyz.x pdbxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip pdbxyz.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o polarize.x polarize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip polarize.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o poledit.x poledit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip poledit.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o potential.x potential.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip potential.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o prmedit.x prmedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip prmedit.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o protein.x protein.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip protein.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o pss.x pss.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip pss.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o pssrigid.x pssrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip pssrigid.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o pssrot.x pssrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip pssrot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o radial.x radial.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip radial.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o saddle.x saddle.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip saddle.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o scan.x scan.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip scan.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o sniffer.x sniffer.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip sniffer.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o spacefill.x spacefill.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip spacefill.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o spectrum.x spectrum.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip spectrum.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o superpose.x superpose.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip superpose.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o testgrad.x testgrad.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip testgrad.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o testhess.x testhess.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip testhess.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o testpair.x testpair.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip testpair.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o testpol.x testpol.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip testpol.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o testrot.x testrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip testrot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o testvir.x testvir.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip testvir.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o timer.x timer.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip timer.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o timerot.x timerot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip timerot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o torsfit.x torsfit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip torsfit.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o valence.x valence.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip valence.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o vibbig.x vibbig.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip vibbig.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o vibrate.x vibrate.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip vibrate.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o vibrot.x vibrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip vibrot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o xtalfit.x xtalfit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip xtalfit.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o xtalmin.x xtalmin.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip xtalmin.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o xyzedit.x xyzedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip xyzedit.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o xyzint.x xyzint.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip xyzint.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o xyzmol2.x xyzmol2.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip xyzmol2.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o xyzpdb.x xyzpdb.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 $APBSLIBS ; strip xyzpdb.x
