#
#
#  #################################################################
#  ##                                                             ##
#  ##  linkprof.make  --  link the Tinker programs for profiling  ##
#  ##              (Intel Fortran for MacOS Version)              ##
#  ##                                                             ##
#  #################################################################
#
#
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o alchemy.x alchemy.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o analyze.x analyze.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o anneal.x anneal.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o archive.x archive.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o bar.x bar.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o correlate.x correlate.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o crystal.x crystal.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o diffuse.x diffuse.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o distgeom.x distgeom.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o document.x document.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o dynamic.x dynamic.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o gda.x gda.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o intedit.x intedit.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o intxyz.x intxyz.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o minimize.x minimize.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o minirot.x minirot.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o minrigid.x minrigid.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o molxyz.x molxyz.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o monte.x monte.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o newton.x newton.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o newtrot.x newtrot.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o nucleic.x nucleic.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o optimize.x optimize.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o optirot.x optirot.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o optrigid.x optrigid.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o path.x path.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o pdbxyz.x pdbxyz.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o polarize.x polarize.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o poledit.x poledit.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o potential.x potential.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o prmedit.x prmedit.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o protein.x protein.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o pss.x pss.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o pssrigid.x pssrigid.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o pssrot.x pssrot.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o radial.x radial.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o saddle.x saddle.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o scan.x scan.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o sniffer.x sniffer.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o spacefill.x spacefill.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o spectrum.x spectrum.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o superpose.x superpose.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o sybylxyz.x sybylxyz.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o testgrad.x testgrad.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o testhess.x testhess.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o testpair.x testpair.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o testpol.x testpol.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o testrot.x testrot.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o timer.x timer.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o timerot.x timerot.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o torsfit.x torsfit.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o valence.x valence.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o vibbig.x vibbig.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o vibrate.x vibrate.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o vibrot.x vibrot.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o xtalfit.x xtalfit.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o xtalmin.x xtalmin.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o xyzedit.x xyzedit.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o xyzint.x xyzint.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o xyzpdb.x xyzpdb.o libtinker.a libfftw3_threads.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -qopenmp -static-intel -mmacosx-version-min=10.6 -o xyzsybyl.x xyzsybyl.o libtinker.a libfftw3_threads.a libfftw3.a
