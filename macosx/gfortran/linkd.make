#
#
#  ##############################################################
#  ##                                                          ##
#  ##  linkd.make  --  link the TINKER programs for debugging  ##
#  ##            (GNU gfortran for Mac OSX Version)            ##
#  ##                                                          ##
#  ##############################################################
#
#
gfortran -g -o alchemy.x alchemy.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o analyze.x analyze.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o anneal.x anneal.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o archive.x archive.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o bar.x bar.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o correlate.x correlate.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o crystal.x crystal.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o diffuse.x diffuse.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o distgeom.x distgeom.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o document.x document.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o dynamic.x dynamic.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o gda.x gda.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o intedit.x intedit.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o intxyz.x intxyz.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o minimize.x minimize.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o minirot.x minirot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o minrigid.x minrigid.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o molxyz.x molxyz.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o monte.x monte.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o newton.x newton.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o newtrot.x newtrot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o nucleic.x nucleic.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o optimize.x optimize.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o optirot.x optirot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o optrigid.x optrigid.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o path.x path.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o pdbxyz.x pdbxyz.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o polarize.x polarize.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o poledit.x poledit.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o potential.x potential.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o prmedit.x prmedit.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o protein.x protein.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o pss.x pss.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o pssrigid.x pssrigid.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o pssrot.x pssrot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o radial.x radial.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o saddle.x saddle.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o scan.x scan.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o sniffer.x sniffer.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o spacefill.x spacefill.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o spectrum.x spectrum.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o superpose.x superpose.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o sybylxyz.x sybylxyz.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o testgrad.x testgrad.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o testhess.x testhess.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o testpair.x testpair.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o testpol.x testpol.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o testrot.x testrot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o timer.x timer.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o timerot.x timerot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o torsfit.x torsfit.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o valence.x valence.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o vibbig.x vibbig.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o vibrate.x vibrate.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o vibrot.x vibrot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o xtalfit.x xtalfit.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o xtalmin.x xtalmin.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o xyzedit.x xyzedit.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o xyzint.x xyzint.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o xyzpdb.x xyzpdb.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -o xyzsybyl.x xyzsybyl.o libtinker.a libfftw3_threads.a libfftw3.a
