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
gfortran -g -static-libgcc -o alchemy.x alchemy.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o analyze.x analyze.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o anneal.x anneal.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o archive.x archive.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o bar.x bar.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o correlate.x correlate.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o crystal.x crystal.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o diffuse.x diffuse.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o distgeom.x distgeom.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o document.x document.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o dynamic.x dynamic.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o gda.x gda.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o intedit.x intedit.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o intxyz.x intxyz.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o minimize.x minimize.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o minirot.x minirot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o minrigid.x minrigid.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o molxyz.x molxyz.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o monte.x monte.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o newton.x newton.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o newtrot.x newtrot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o nucleic.x nucleic.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o optimize.x optimize.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o optirot.x optirot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o optrigid.x optrigid.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o path.x path.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o pdbxyz.x pdbxyz.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o polarize.x polarize.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o poledit.x poledit.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o potential.x potential.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o prmedit.x prmedit.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o protein.x protein.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o pss.x pss.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o pssrigid.x pssrigid.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o pssrot.x pssrot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o radial.x radial.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o saddle.x saddle.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o scan.x scan.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o sniffer.x sniffer.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o spacefill.x spacefill.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o spectrum.x spectrum.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o superpose.x superpose.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o sybylxyz.x sybylxyz.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o testgrad.x testgrad.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o testhess.x testhess.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o testpair.x testpair.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o testpol.x testpol.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o testrot.x testrot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o timer.x timer.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o timerot.x timerot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o torsfit.x torsfit.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o valence.x valence.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o vibbig.x vibbig.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o vibrate.x vibrate.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o vibrot.x vibrot.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o xtalfit.x xtalfit.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o xtalmin.x xtalmin.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o xyzedit.x xyzedit.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o xyzint.x xyzint.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o xyzpdb.x xyzpdb.o libtinker.a libfftw3_threads.a libfftw3.a
gfortran -g -static-libgcc -o xyzsybyl.x xyzsybyl.o libtinker.a libfftw3_threads.a libfftw3.a
