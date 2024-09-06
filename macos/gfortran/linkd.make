#
#
#  ##############################################################
#  ##                                                          ##
#  ##  linkd.make  --  link the Tinker programs for debugging  ##
#  ##             (GNU gfortran for macOS Version)             ##
#  ##                                                          ##
#  ##############################################################
#
#
gfortran -Og -g -fbacktrace -static-libgcc -o alchemy.x alchemy.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o analyze.x analyze.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o anneal.x anneal.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o arcedit.x arcedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o bar.x bar.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o correlate.x correlate.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o critical.x critical.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o crystal.x crystal.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o diffuse.x diffuse.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o distgeom.x distgeom.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o document.x document.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o dynamic.x dynamic.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o freefix.x freefix.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o gda.x gda.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o intedit.x intedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o intxyz.x intxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o minimize.x minimize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o minirot.x minirot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o minrigid.x minrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o mol2xyz.x mol2xyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o molxyz.x molxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o monte.x monte.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o newton.x newton.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o newtrot.x newtrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o nucleic.x nucleic.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o optimize.x optimize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o optirot.x optirot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o optrigid.x optrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o path.x path.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o pdbxyz.x pdbxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o polarize.x polarize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o poledit.x poledit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o potential.x potential.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o prmedit.x prmedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o protein.x protein.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o pss.x pss.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o pssrigid.x pssrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o pssrot.x pssrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o radial.x radial.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o saddle.x saddle.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o scan.x scan.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o sniffer.x sniffer.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o spacefill.x spacefill.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o spectrum.x spectrum.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o superpose.x superpose.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o testgrad.x testgrad.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o testhess.x testhess.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o testpair.x testpair.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o testpol.x testpol.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o testrot.x testrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o testsurf.x testsurf.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o testvir.x testvir.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o timer.x timer.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o timerot.x timerot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o torsfit.x torsfit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o valence.x valence.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o vibbig.x vibbig.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o vibrate.x vibrate.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o vibrot.x vibrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o xtalfit.x xtalfit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o xtalmin.x xtalmin.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o xyzedit.x xyzedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o xyzint.x xyzint.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o xyzmol2.x xyzmol2.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
gfortran -Og -g -fbacktrace -static-libgcc -o xyzpdb.x xyzpdb.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
