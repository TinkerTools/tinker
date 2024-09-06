#
#
#  ################################################################
#  ##                                                            ##
#  ##  link.make  --  link debug version of the Tinker programs  ##
#  ##              (Intel Fortran for Linux Version)             ##
#  ##                                                            ##
#  ################################################################
#
#
ifort -g -o alchemy.x alchemy.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o analyze.x analyze.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o anneal.x anneal.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o arcedit.x arcedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o bar.x bar.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o correlate.x correlate.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o critical.x critical.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o crystal.x crystal.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o diffuse.x diffuse.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o distgeom.x distgeom.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o document.x document.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o dynamic.x dynamic.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o freefix.x freefix.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o gda.x gda.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o intedit.x intedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o intxyz.x intxyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o minimize.x minimize.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o minirot.x minirot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o minrigid.x minrigid.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o mol2xyz.x mol2xyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o molxyz.x molxyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o monte.x monte.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o newton.x newton.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o newtrot.x newtrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o nucleic.x nucleic.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o optimize.x optimize.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o optirot.x optirot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o optrigid.x optrigid.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o path.x path.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o pdbxyz.x pdbxyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o polarize.x polarize.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o poledit.x poledit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o potential.x potential.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o prmedit.x prmedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o protein.x protein.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o pss.x pss.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o pssrigid.x pssrigid.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o pssrot.x pssrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o radial.x radial.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o saddle.x saddle.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o scan.x scan.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o sniffer.x sniffer.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o spacefill.x spacefill.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o spectrum.x spectrum.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o superpose.x superpose.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o testgrad.x testgrad.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o testhess.x testhess.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o testpair.x testpair.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o testpol.x testpol.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o testrot.x testrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o testsurf.x testsurf.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o testvir.x testvir.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o timer.x timer.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o timerot.x timerot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o torsfit.x torsfit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o valence.x valence.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o vibbig.x vibbig.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o vibrate.x vibrate.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o vibrot.x vibrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o xtalfit.x xtalfit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o xtalmin.x xtalmin.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o xyzedit.x xyzedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o xyzint.x xyzint.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o xyzmol2.x xyzmol2.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
ifort -g -o xyzpdb.x xyzpdb.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3
