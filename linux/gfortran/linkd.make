#
#
#  ##############################################################
#  ##                                                          ##
#  ##  linkd.make  --  link the Tinker programs for debugging  ##
#  ##             (GNU gfortran for Linux Version)             ##
#  ##                                                          ##
#  ##############################################################
#
#
#  if the -static-libgcc flag is not enforced, then copy static
#  versions of the main gcc libraries (libgcc.a, libgfortran.a,
#  libgomp.a and libquadmath.a) into the present build directory;
#  with Linux, all four libraries must be present locally, while
#  for MacOS, libquadmath.a is the only one needed
#
#
gfortran -Og -g -fbacktrace -static-libgcc -o alchemy.x alchemy.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip alchemy.x
gfortran -Og -g -fbacktrace -static-libgcc -o analyze.x analyze.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip analyze.x
gfortran -Og -g -fbacktrace -static-libgcc -o anneal.x anneal.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip anneal.x
gfortran -Og -g -fbacktrace -static-libgcc -o arcedit.x arcedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip arcedit.x
gfortran -Og -g -fbacktrace -static-libgcc -o bar.x bar.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip bar.x
gfortran -Og -g -fbacktrace -static-libgcc -o correlate.x correlate.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip correlate.x
gfortran -Og -g -fbacktrace -static-libgcc -o critical.x critical.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip critical.x
gfortran -Og -g -fbacktrace -static-libgcc -o crystal.x crystal.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip crystal.x
gfortran -Og -g -fbacktrace -static-libgcc -o diffuse.x diffuse.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip diffuse.x
gfortran -Og -g -fbacktrace -static-libgcc -o distgeom.x distgeom.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip distgeom.x
gfortran -Og -g -fbacktrace -static-libgcc -o document.x document.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip document.x
gfortran -Og -g -fbacktrace -static-libgcc -o dynamic.x dynamic.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip dynamic.x
gfortran -Og -g -fbacktrace -static-libgcc -o freefix.x freefix.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip freefix.x
gfortran -Og -g -fbacktrace -static-libgcc -o gda.x gda.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip gda.x
gfortran -Og -g -fbacktrace -static-libgcc -o intedit.x intedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip intedit.x
gfortran -Og -g -fbacktrace -static-libgcc -o intxyz.x intxyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip intxyz.x
gfortran -Og -g -fbacktrace -static-libgcc -o minimize.x minimize.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minimize.x
gfortran -Og -g -fbacktrace -static-libgcc -o minirot.x minirot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minirot.x
gfortran -Og -g -fbacktrace -static-libgcc -o minrigid.x minrigid.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minrigid.x
gfortran -Og -g -fbacktrace -static-libgcc -o mol2xyz.x mol2xyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip mol2xyz.x
gfortran -Og -g -fbacktrace -static-libgcc -o molxyz.x molxyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip molxyz.x
gfortran -Og -g -fbacktrace -static-libgcc -o monte.x monte.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip monte.x
gfortran -Og -g -fbacktrace -static-libgcc -o newton.x newton.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip newton.x
gfortran -Og -g -fbacktrace -static-libgcc -o newtrot.x newtrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip newtrot.x
gfortran -Og -g -fbacktrace -static-libgcc -o nucleic.x nucleic.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip nucleic.x
gfortran -Og -g -fbacktrace -static-libgcc -o optimize.x optimize.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optimize.x
gfortran -Og -g -fbacktrace -static-libgcc -o optirot.x optirot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optirot.x
gfortran -Og -g -fbacktrace -static-libgcc -o optrigid.x optrigid.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optrigid.x
gfortran -Og -g -fbacktrace -static-libgcc -o path.x path.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip path.x
gfortran -Og -g -fbacktrace -static-libgcc -o pdbxyz.x pdbxyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pdbxyz.x
gfortran -Og -g -fbacktrace -static-libgcc -o polarize.x polarize.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip polarize.x
gfortran -Og -g -fbacktrace -static-libgcc -o poledit.x poledit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip poledit.x
gfortran -Og -g -fbacktrace -static-libgcc -o potential.x potential.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip potential.x
gfortran -Og -g -fbacktrace -static-libgcc -o prmedit.x prmedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip prmedit.x
gfortran -Og -g -fbacktrace -static-libgcc -o protein.x protein.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip protein.x
gfortran -Og -g -fbacktrace -static-libgcc -o pss.x pss.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pss.x
gfortran -Og -g -fbacktrace -static-libgcc -o pssrigid.x pssrigid.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pssrigid.x
gfortran -Og -g -fbacktrace -static-libgcc -o pssrot.x pssrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pssrot.x
gfortran -Og -g -fbacktrace -static-libgcc -o radial.x radial.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip radial.x
gfortran -Og -g -fbacktrace -static-libgcc -o saddle.x saddle.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip saddle.x
gfortran -Og -g -fbacktrace -static-libgcc -o scan.x scan.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip scan.x
gfortran -Og -g -fbacktrace -static-libgcc -o sniffer.x sniffer.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip sniffer.x
gfortran -Og -g -fbacktrace -static-libgcc -o spacefill.x spacefill.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip spacefill.x
gfortran -Og -g -fbacktrace -static-libgcc -o spectrum.x spectrum.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip spectrum.x
gfortran -Og -g -fbacktrace -static-libgcc -o superpose.x superpose.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip superpose.x
gfortran -Og -g -fbacktrace -static-libgcc -o testgrad.x testgrad.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testgrad.x
gfortran -Og -g -fbacktrace -static-libgcc -o testhess.x testhess.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testhess.x
gfortran -Og -g -fbacktrace -static-libgcc -o testpair.x testpair.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testpair.x
gfortran -Og -g -fbacktrace -static-libgcc -o testpol.x testpol.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testpol.x
gfortran -Og -g -fbacktrace -static-libgcc -o testrot.x testrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testrot.x
gfortran -Og -g -fbacktrace -static-libgcc -o testsurf.x testsurf.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testsurf.x
gfortran -Og -g -fbacktrace -static-libgcc -o testvir.x testvir.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testvir.x
gfortran -Og -g -fbacktrace -static-libgcc -o timer.x timer.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip timer.x
gfortran -Og -g -fbacktrace -static-libgcc -o timerot.x timerot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip timerot.x
gfortran -Og -g -fbacktrace -static-libgcc -o torsfit.x torsfit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip torsfit.x
gfortran -Og -g -fbacktrace -static-libgcc -o valence.x valence.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip valence.x
gfortran -Og -g -fbacktrace -static-libgcc -o vibbig.x vibbig.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibbig.x
gfortran -Og -g -fbacktrace -static-libgcc -o vibrate.x vibrate.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibrate.x
gfortran -Og -g -fbacktrace -static-libgcc -o vibrot.x vibrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibrot.x
gfortran -Og -g -fbacktrace -static-libgcc -o xtalfit.x xtalfit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xtalfit.x
gfortran -Og -g -fbacktrace -static-libgcc -o xtalmin.x xtalmin.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xtalmin.x
gfortran -Og -g -fbacktrace -static-libgcc -o xyzedit.x xyzedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzedit.x
gfortran -Og -g -fbacktrace -static-libgcc -o xyzint.x xyzint.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzint.x
gfortran -Og -g -fbacktrace -static-libgcc -o xyzmol2.x xyzmol2.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzmol2.x
gfortran -Og -g -fbacktrace -static-libgcc -o xyzpdb.x xyzpdb.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzpdb.x
