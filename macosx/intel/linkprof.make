#
#
#  #################################################################
#  ##                                                             ##
#  ##  linkprof.make  --  link the TINKER programs for profiling  ##
#  ##        (Intel Fortran Compiler for Mac OSX Version)         ##
#  ##                                                             ##
#  #################################################################
#
#
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o alchemy.x alchemy.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o analyze.x analyze.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o anneal.x anneal.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o archive.x archive.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o bar.x bar.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o correlate.x correlate.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o crystal.x crystal.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o diffuse.x diffuse.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o distgeom.x distgeom.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o document.x document.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o dynamic.x dynamic.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o gda.x gda.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o intedit.x intedit.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o intxyz.x intxyz.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o minimize.x minimize.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o minirot.x minirot.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o minrigid.x minrigid.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o molxyz.x molxyz.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o monte.x monte.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o newton.x newton.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o newtrot.x newtrot.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o nucleic.x nucleic.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o optimize.x optimize.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o optirot.x optirot.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o optrigid.x optrigid.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o path.x path.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o pdbxyz.x pdbxyz.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o polarize.x polarize.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o poledit.x poledit.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o potential.x potential.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o prmedit.x prmedit.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o protein.x protein.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o pss.x pss.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o pssrigid.x pssrigid.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o pssrot.x pssrot.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o radial.x radial.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o saddle.x saddle.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o scan.x scan.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o sniffer.x sniffer.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o spacefill.x spacefill.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o spectrum.x spectrum.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o superpose.x superpose.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o sybylxyz.x sybylxyz.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o testgrad.x testgrad.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o testhess.x testhess.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o testpair.x testpair.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o testpol.x testpol.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o testrot.x testrot.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o timer.x timer.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o timerot.x timerot.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o torsfit.x torsfit.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o valence.x valence.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o vibbig.x vibbig.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o vibrate.x vibrate.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o vibrot.x vibrot.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o xtalfit.x xtalfit.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o xtalmin.x xtalmin.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o xyzedit.x xyzedit.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o xyzint.x xyzint.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o xyzpdb.x xyzpdb.o libtinker.a libfftw3_omp.a libfftw3.a
ifort -O3 -g -no-ipo -no-prec-div -openmp -static-intel -mmacosx-version-min=10.5 -o xyzsybyl.x xyzsybyl.o libtinker.a libfftw3_omp.a libfftw3.a
