#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##           (Windows/Cygwin/GNU gfortran Version)           ##
#  ##                                                           ##
#  ###############################################################
#
#
gfortran -o alchemy.x alchemy.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o analyze.x analyze.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o anneal.x anneal.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o archive.x archive.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o bar.x bar.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o correlate.x correlate.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o crystal.x crystal.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o diffuse.x diffuse.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o distgeom.x distgeom.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o document.x document.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o dynamic.x dynamic.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o gda.x gda.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o intedit.x intedit.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o intxyz.x intxyz.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o minimize.x minimize.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o minirot.x minirot.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o minrigid.x minrigid.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o molxyz.x molxyz.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o monte.x monte.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o newton.x newton.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o newtrot.x newtrot.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o nucleic.x nucleic.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o optimize.x optimize.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o optirot.x optirot.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o optrigid.x optrigid.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o path.x path.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o pdbxyz.x pdbxyz.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o polarize.x polarize.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o poledit.x poledit.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o potential.x potential.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o prmedit.x prmedit.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o protein.x protein.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o pss.x pss.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o pssrigid.x pssrigid.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o pssrot.x pssrot.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o radial.x radial.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o saddle.x saddle.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o scan.x scan.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o sniffer.x sniffer.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o spacefill.x spacefill.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o spectrum.x spectrum.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o superpose.x superpose.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o sybylxyz.x sybylxyz.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o testgrad.x testgrad.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o testhess.x testhess.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o testpair.x testpair.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o testpol.x testpol.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o testrot.x testrot.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o timer.x timer.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o timerot.x timerot.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o torsfit.x torsfit.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o valence.x valence.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o vibbig.x vibbig.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o vibrate.x vibrate.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o vibrot.x vibrot.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o xtalfit.x xtalfit.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o xtalmin.x xtalmin.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o xyzedit.x xyzedit.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o xyzint.x xyzint.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o xyzpdb.x xyzpdb.o libtinker.a libfftw3.a libfftw3_threads.a
gfortran -o xyzsybyl.x xyzsybyl.o libtinker.a libfftw3.a libfftw3_threads.a
