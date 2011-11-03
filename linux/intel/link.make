#
#
#  #################################################################
#  ##                                                             ##
#  ##  link.make  --  link all of the TINKER programs for OpenMP  ##
#  ##         (Intel Fortran Compiler for Linux Version)          ##
#  ##                                                             ##
#  #################################################################
#
#
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o alchemy.x alchemy.o libtinker.a libfftw3_omp.a libfftw3.a ; strip alchemy.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o analyze.x analyze.o libtinker.a libfftw3_omp.a libfftw3.a ; strip analyze.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o anneal.x anneal.o libtinker.a libfftw3_omp.a libfftw3.a ; strip anneal.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o archive.x archive.o libtinker.a libfftw3_omp.a libfftw3.a ; strip archive.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o correlate.x correlate.o libtinker.a libfftw3_omp.a libfftw3.a ; strip correlate.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o crystal.x crystal.o libtinker.a libfftw3_omp.a libfftw3.a ; strip crystal.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o diffuse.x diffuse.o libtinker.a libfftw3_omp.a libfftw3.a ; strip diffuse.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o distgeom.x distgeom.o libtinker.a libfftw3_omp.a libfftw3.a ; strip distgeom.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o document.x document.o libtinker.a libfftw3_omp.a libfftw3.a ; strip document.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o dynamic.x dynamic.o libtinker.a libfftw3_omp.a libfftw3.a ; strip dynamic.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o gda.x gda.o libtinker.a libfftw3_omp.a libfftw3.a ; strip gda.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o intedit.x intedit.o libtinker.a libfftw3_omp.a libfftw3.a ; strip intedit.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o intxyz.x intxyz.o libtinker.a libfftw3_omp.a libfftw3.a ; strip intxyz.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o minimize.x minimize.o libtinker.a libfftw3_omp.a libfftw3.a ; strip minimize.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o minirot.x minirot.o libtinker.a libfftw3_omp.a libfftw3.a ; strip minirot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o minrigid.x minrigid.o libtinker.a libfftw3_omp.a libfftw3.a ; strip minrigid.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o monte.x monte.o libtinker.a libfftw3_omp.a libfftw3.a ; strip monte.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o newton.x newton.o libtinker.a libfftw3_omp.a libfftw3.a ; strip newton.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o newtrot.x newtrot.o libtinker.a libfftw3_omp.a libfftw3.a ; strip newtrot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o nucleic.x nucleic.o libtinker.a libfftw3_omp.a libfftw3.a ; strip nucleic.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o optimize.x optimize.o libtinker.a libfftw3_omp.a libfftw3.a ; strip optimize.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o optirot.x optirot.o libtinker.a libfftw3_omp.a libfftw3.a ; strip optirot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o optrigid.x optrigid.o libtinker.a libfftw3_omp.a libfftw3.a ; strip optrigid.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o path.x path.o libtinker.a libfftw3_omp.a libfftw3.a ; strip path.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o pdbxyz.x pdbxyz.o libtinker.a libfftw3_omp.a libfftw3.a ; strip pdbxyz.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o polarize.x polarize.o libtinker.a libfftw3_omp.a libfftw3.a ; strip polarize.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o poledit.x poledit.o libtinker.a libfftw3_omp.a libfftw3.a ; strip poledit.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o potential.x potential.o libtinker.a libfftw3_omp.a libfftw3.a ; strip potential.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o prmedit.x prmedit.o libtinker.a libfftw3_omp.a libfftw3.a ; strip prmedit.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o protein.x protein.o libtinker.a libfftw3_omp.a libfftw3.a ; strip protein.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o pss.x pss.o libtinker.a libfftw3_omp.a libfftw3.a ; strip pss.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o pssrigid.x pssrigid.o libtinker.a libfftw3_omp.a libfftw3.a ; strip pssrigid.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o pssrot.x pssrot.o libtinker.a libfftw3_omp.a libfftw3.a ; strip pssrot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o radial.x radial.o libtinker.a libfftw3_omp.a libfftw3.a ; strip radial.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o saddle.x saddle.o libtinker.a libfftw3_omp.a libfftw3.a ; strip saddle.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o scan.x scan.o libtinker.a libfftw3_omp.a libfftw3.a ; strip scan.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o sniffer.x sniffer.o libtinker.a libfftw3_omp.a libfftw3.a ; strip sniffer.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o spacefill.x spacefill.o libtinker.a libfftw3_omp.a libfftw3.a ; strip spacefill.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o spectrum.x spectrum.o libtinker.a libfftw3_omp.a libfftw3.a ; strip spectrum.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o superpose.x superpose.o libtinker.a libfftw3_omp.a libfftw3.a ; strip superpose.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o sybylxyz.x sybylxyz.o libtinker.a libfftw3_omp.a libfftw3.a ; strip sybylxyz.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o testgrad.x testgrad.o libtinker.a libfftw3_omp.a libfftw3.a ; strip testgrad.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o testhess.x testhess.o libtinker.a libfftw3_omp.a libfftw3.a ; strip testhess.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o testpair.x testpair.o libtinker.a libfftw3_omp.a libfftw3.a ; strip testpair.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o testrot.x testrot.o libtinker.a libfftw3_omp.a libfftw3.a ; strip testrot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o timer.x timer.o libtinker.a libfftw3_omp.a libfftw3.a ; strip timer.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o timerot.x timerot.o libtinker.a libfftw3_omp.a libfftw3.a ; strip timerot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o torsfit.x torsfit.o libtinker.a libfftw3_omp.a libfftw3.a ; strip torsfit.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o valence.x valence.o libtinker.a libfftw3_omp.a libfftw3.a ; strip valence.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o vibbig.x vibbig.o libtinker.a libfftw3_omp.a libfftw3.a ; strip vibbig.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o vibrate.x vibrate.o libtinker.a libfftw3_omp.a libfftw3.a ; strip vibrate.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o vibrot.x vibrot.o libtinker.a libfftw3_omp.a libfftw3.a ; strip vibrot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o xtalfit.x xtalfit.o libtinker.a libfftw3_omp.a libfftw3.a ; strip xtalfit.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o xtalmin.x xtalmin.o libtinker.a libfftw3_omp.a libfftw3.a ; strip xtalmin.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o xyzedit.x xyzedit.o libtinker.a libfftw3_omp.a libfftw3.a ; strip xyzedit.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o xyzint.x xyzint.o libtinker.a libfftw3_omp.a libfftw3.a ; strip xyzint.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o xyzpdb.x xyzpdb.o libtinker.a libfftw3_omp.a libfftw3.a ; strip xyzpdb.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto-scalar -static-intel -o xyzsybyl.x xyzsybyl.o libtinker.a libfftw3_omp.a libfftw3.a ; strip xyzsybyl.x
