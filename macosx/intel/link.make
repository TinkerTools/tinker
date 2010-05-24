#
#
#  #################################################################
#  ##                                                             ##
#  ##  link.make  --  link all of the TINKER programs for OpenMP  ##
#  ##        (Intel Fortran Compiler for Mac OSX Version)         ##
#  ##                                                             ##
#  #################################################################
#
#
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o alchemy.x alchemy.o libtinker.a libfftw3_threads.a libfftw3.a ; strip alchemy.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o analyze.x analyze.o libtinker.a libfftw3_threads.a libfftw3.a ; strip analyze.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o anneal.x anneal.o libtinker.a libfftw3_threads.a libfftw3.a ; strip anneal.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o archive.x archive.o libtinker.a libfftw3_threads.a libfftw3.a ; strip archive.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o correlate.x correlate.o libtinker.a libfftw3_threads.a libfftw3.a ; strip correlate.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o crystal.x crystal.o libtinker.a libfftw3_threads.a libfftw3.a ; strip crystal.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o diffuse.x diffuse.o libtinker.a libfftw3_threads.a libfftw3.a ; strip diffuse.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o distgeom.x distgeom.o libtinker.a libfftw3_threads.a libfftw3.a ; strip distgeom.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o document.x document.o libtinker.a libfftw3_threads.a libfftw3.a ; strip document.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o dynamic.x dynamic.o libtinker.a libfftw3_threads.a libfftw3.a ; strip dynamic.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o gda.x gda.o libtinker.a libfftw3_threads.a libfftw3.a ; strip gda.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o intedit.x intedit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip intedit.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o intxyz.x intxyz.o libtinker.a libfftw3_threads.a libfftw3.a ; strip intxyz.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o minimize.x minimize.o libtinker.a libfftw3_threads.a libfftw3.a ; strip minimize.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o minirot.x minirot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip minirot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o minrigid.x minrigid.o libtinker.a libfftw3_threads.a libfftw3.a ; strip minrigid.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o monte.x monte.o libtinker.a libfftw3_threads.a libfftw3.a ; strip monte.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o newton.x newton.o libtinker.a libfftw3_threads.a libfftw3.a ; strip newton.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o newtrot.x newtrot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip newtrot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o nucleic.x nucleic.o libtinker.a libfftw3_threads.a libfftw3.a ; strip nucleic.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o optimize.x optimize.o libtinker.a libfftw3_threads.a libfftw3.a ; strip optimize.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o optirot.x optirot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip optirot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o optrigid.x optrigid.o libtinker.a libfftw3_threads.a libfftw3.a ; strip optrigid.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o path.x path.o libtinker.a libfftw3_threads.a libfftw3.a ; strip path.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o pdbxyz.x pdbxyz.o libtinker.a libfftw3_threads.a libfftw3.a ; strip pdbxyz.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o polarize.x polarize.o libtinker.a libfftw3_threads.a libfftw3.a ; strip polarize.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o poledit.x poledit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip poledit.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o potential.x potential.o libtinker.a libfftw3_threads.a libfftw3.a ; strip potential.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o prmedit.x prmedit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip prmedit.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o protein.x protein.o libtinker.a libfftw3_threads.a libfftw3.a ; strip protein.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o pss.x pss.o libtinker.a libfftw3_threads.a libfftw3.a ; strip pss.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o pssrigid.x pssrigid.o libtinker.a libfftw3_threads.a libfftw3.a ; strip pssrigid.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o pssrot.x pssrot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip pssrot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o radial.x radial.o libtinker.a libfftw3_threads.a libfftw3.a ; strip radial.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o saddle.x saddle.o libtinker.a libfftw3_threads.a libfftw3.a ; strip saddle.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o scan.x scan.o libtinker.a libfftw3_threads.a libfftw3.a ; strip scan.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o sniffer.x sniffer.o libtinker.a libfftw3_threads.a libfftw3.a ; strip sniffer.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o spacefill.x spacefill.o libtinker.a libfftw3_threads.a libfftw3.a ; strip spacefill.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o spectrum.x spectrum.o libtinker.a libfftw3_threads.a libfftw3.a ; strip spectrum.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o superpose.x superpose.o libtinker.a libfftw3_threads.a libfftw3.a ; strip superpose.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o sybylxyz.x sybylxyz.o libtinker.a libfftw3_threads.a libfftw3.a ; strip sybylxyz.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o testgrad.x testgrad.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testgrad.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o testhess.x testhess.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testhess.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o testpair.x testpair.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testpair.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o testrot.x testrot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testrot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o timer.x timer.o libtinker.a libfftw3_threads.a libfftw3.a ; strip timer.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o timerot.x timerot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip timerot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o valence.x valence.o libtinker.a libfftw3_threads.a libfftw3.a ; strip valence.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o vibbig.x vibbig.o libtinker.a libfftw3_threads.a libfftw3.a ; strip vibbig.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o vibrate.x vibrate.o libtinker.a libfftw3_threads.a libfftw3.a ; strip vibrate.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o vibrot.x vibrot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip vibrot.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o xtalfit.x xtalfit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xtalfit.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o xtalmin.x xtalmin.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xtalmin.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o xyzedit.x xyzedit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzedit.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o xyzint.x xyzint.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzint.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o xyzpdb.x xyzpdb.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzpdb.x
ifort -O3 -no-prec-div -fno-omit-frame-pointer -openmp -recursive -auto_scalar -mmacosx-version-min=10.4 -o xyzsybyl.x xyzsybyl.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzsybyl.x
