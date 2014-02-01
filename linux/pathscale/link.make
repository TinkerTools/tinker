#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##           (PathScale pathf90 for Linux Version)           ##
#  ##                                                           ##
#  ###############################################################
#
#
pathf90 -O -msse3 -openmp -o alchemy.x alchemy.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o analyze.x analyze.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o anneal.x anneal.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o archive.x archive.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o bar.x bar.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o correlate.x correlate.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o crystal.x crystal.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o diffuse.x diffuse.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o distgeom.x distgeom.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o document.x document.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o dynamic.x dynamic.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o gda.x gda.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o intedit.x intedit.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o intxyz.x intxyz.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o minimize.x minimize.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o minirot.x minirot.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o minrigid.x minrigid.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o molxyz.x molxyz.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o monte.x monte.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o newton.x newton.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o newtrot.x newtrot.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o nucleic.x nucleic.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o optimize.x optimize.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o optirot.x optirot.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o optrigid.x optrigid.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o path.x path.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o pdbxyz.x pdbxyz.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o polarize.x polarize.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o poledit.x poledit.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o potential.x potential.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o prmedit.x prmedit.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o protein.x protein.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o pss.x pss.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o pssrigid.x pssrigid.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o pssrot.x pssrot.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o radial.x radial.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o saddle.x saddle.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o scan.x scan.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o sniffer.x sniffer.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o spacefill.x spacefill.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o spectrum.x spectrum.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o superpose.x superpose.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o sybylxyz.x sybylxyz.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o testgrad.x testgrad.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o testhess.x testhess.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o testpair.x testpair.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o testpol.x testpol.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o testrot.x testrot.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o timer.x timer.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o timerot.x timerot.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o torsfit.x torsfit.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o valence.x valence.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o vibbig.x vibbig.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o vibrate.x vibrate.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o vibrot.x vibrot.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o xtalfit.x xtalfit.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o xtalmin.x xtalmin.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o xyzedit.x xyzedit.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o xyzint.x xyzint.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o xyzpdb.x xyzpdb.o libtinker.a libfftw3_threads.a libfftw3.a
pathf90 -O -msse3 -openmp -o xyzsybyl.x xyzsybyl.o libtinker.a libfftw3_threads.a libfftw3.a
