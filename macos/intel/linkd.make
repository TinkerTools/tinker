#
#
#  ##############################################################
#  ##                                                          ##
#  ##  linkd.make  --  link the Tinker programs for debugging  ##
#  ##             (Intel Fortran for macOS Version)            ##
#  ##                                                          ##
#  ##############################################################
#
#
ifort -g -Wl,-no_pie -o alchemy.x alchemy.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o analyze.x analyze.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o anneal.x anneal.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o arcedit.x arcedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o bar.x bar.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o correlate.x correlate.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o critical.x critical.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o crystal.x crystal.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o diffuse.x diffuse.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o distgeom.x distgeom.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o document.x document.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o dynamic.x dynamic.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o freefix.x freefix.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o gda.x gda.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o intedit.x intedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o intxyz.x intxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o minimize.x minimize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o minirot.x minirot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o minrigid.x minrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o mol2xyz.x mol2xyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o molxyz.x molxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o monte.x monte.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o newton.x newton.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o newtrot.x newtrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o nucleic.x nucleic.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o optimize.x optimize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o optirot.x optirot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o optrigid.x optrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o path.x path.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o pdbxyz.x pdbxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o polarize.x polarize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o poledit.x poledit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o potential.x potential.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o prmedit.x prmedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o protein.x protein.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o pss.x pss.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o pssrigid.x pssrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o pssrot.x pssrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o radial.x radial.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o saddle.x saddle.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o scan.x scan.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o sniffer.x sniffer.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o spacefill.x spacefill.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o spectrum.x spectrum.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o superpose.x superpose.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o testgrad.x testgrad.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o testhess.x testhess.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o testpair.x testpair.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o testpol.x testpol.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o testrot.x testrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o testsurf.x testsurf.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o testvir.x testvir.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o timer.x timer.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o timerot.x timerot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o torsfit.x torsfit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o valence.x valence.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o vibbig.x vibbig.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o vibrate.x vibrate.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o vibrot.x vibrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o xtalfit.x xtalfit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o xtalmin.x xtalmin.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o xyzedit.x xyzedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o xyzint.x xyzint.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o xyzmol2.x xyzmol2.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -Wl,-no_pie -o xyzpdb.x xyzpdb.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
