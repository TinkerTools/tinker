#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the Tinker package programs  ##
#  ##             (Intel Fortran for Linux Version)             ##
#  ##                                                           ##
#  ###############################################################
#
#
ifort -O3 -no-prec-div -recursive -openmp -static -o alchemy.x alchemy.o libtinker.a libfftw3_threads.a libfftw3.a ; strip alchemy.x
ifort -O3 -no-prec-div -recursive -openmp -static -o analyze.x analyze.o libtinker.a libfftw3_threads.a libfftw3.a ; strip analyze.x
ifort -O3 -no-prec-div -recursive -openmp -static -o anneal.x anneal.o libtinker.a libfftw3_threads.a libfftw3.a ; strip anneal.x
ifort -O3 -no-prec-div -recursive -openmp -static -o archive.x archive.o libtinker.a libfftw3_threads.a libfftw3.a ; strip archive.x
ifort -O3 -no-prec-div -recursive -openmp -static -o bar.x bar.o libtinker.a libfftw3_threads.a libfftw3.a ; strip bar.x
ifort -O3 -no-prec-div -recursive -openmp -static -o correlate.x correlate.o libtinker.a libfftw3_threads.a libfftw3.a ; strip correlate.x
ifort -O3 -no-prec-div -recursive -openmp -static -o crystal.x crystal.o libtinker.a libfftw3_threads.a libfftw3.a ; strip crystal.x
ifort -O3 -no-prec-div -recursive -openmp -static -o diffuse.x diffuse.o libtinker.a libfftw3_threads.a libfftw3.a ; strip diffuse.x
ifort -O3 -no-prec-div -recursive -openmp -static -o distgeom.x distgeom.o libtinker.a libfftw3_threads.a libfftw3.a ; strip distgeom.x
ifort -O3 -no-prec-div -recursive -openmp -static -o document.x document.o libtinker.a libfftw3_threads.a libfftw3.a ; strip document.x
ifort -O3 -no-prec-div -recursive -openmp -static -o dynamic.x dynamic.o libtinker.a libfftw3_threads.a libfftw3.a ; strip dynamic.x
ifort -O3 -no-prec-div -recursive -openmp -static -o gda.x gda.o libtinker.a libfftw3_threads.a libfftw3.a ; strip gda.x
ifort -O3 -no-prec-div -recursive -openmp -static -o intedit.x intedit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip intedit.x
ifort -O3 -no-prec-div -recursive -openmp -static -o intxyz.x intxyz.o libtinker.a libfftw3_threads.a libfftw3.a ; strip intxyz.x
ifort -O3 -no-prec-div -recursive -openmp -static -o minimize.x minimize.o libtinker.a libfftw3_threads.a libfftw3.a ; strip minimize.x
ifort -O3 -no-prec-div -recursive -openmp -static -o minirot.x minirot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip minirot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o minrigid.x minrigid.o libtinker.a libfftw3_threads.a libfftw3.a ; strip minrigid.x
ifort -O3 -no-prec-div -recursive -openmp -static -o mol2xyz.x mol2xyz.o libtinker.a libfftw3_threads.a libfftw3.a ; strip mol2xyz.x
ifort -O3 -no-prec-div -recursive -openmp -static -o molxyz.x molxyz.o libtinker.a libfftw3_threads.a libfftw3.a ; strip molxyz.x
ifort -O3 -no-prec-div -recursive -openmp -static -o monte.x monte.o libtinker.a libfftw3_threads.a libfftw3.a ; strip monte.x
ifort -O3 -no-prec-div -recursive -openmp -static -o newton.x newton.o libtinker.a libfftw3_threads.a libfftw3.a ; strip newton.x
ifort -O3 -no-prec-div -recursive -openmp -static -o newtrot.x newtrot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip newtrot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o nucleic.x nucleic.o libtinker.a libfftw3_threads.a libfftw3.a ; strip nucleic.x
ifort -O3 -no-prec-div -recursive -openmp -static -o optimize.x optimize.o libtinker.a libfftw3_threads.a libfftw3.a ; strip optimize.x
ifort -O3 -no-prec-div -recursive -openmp -static -o optirot.x optirot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip optirot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o optrigid.x optrigid.o libtinker.a libfftw3_threads.a libfftw3.a ; strip optrigid.x
ifort -O3 -no-prec-div -recursive -openmp -static -o path.x path.o libtinker.a libfftw3_threads.a libfftw3.a ; strip path.x
ifort -O3 -no-prec-div -recursive -openmp -static -o pdbxyz.x pdbxyz.o libtinker.a libfftw3_threads.a libfftw3.a ; strip pdbxyz.x
ifort -O3 -no-prec-div -recursive -openmp -static -o polarize.x polarize.o libtinker.a libfftw3_threads.a libfftw3.a ; strip polarize.x
ifort -O3 -no-prec-div -recursive -openmp -static -o poledit.x poledit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip poledit.x
ifort -O3 -no-prec-div -recursive -openmp -static -o potential.x potential.o libtinker.a libfftw3_threads.a libfftw3.a ; strip potential.x
ifort -O3 -no-prec-div -recursive -openmp -static -o prmedit.x prmedit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip prmedit.x
ifort -O3 -no-prec-div -recursive -openmp -static -o protein.x protein.o libtinker.a libfftw3_threads.a libfftw3.a ; strip protein.x
ifort -O3 -no-prec-div -recursive -openmp -static -o pss.x pss.o libtinker.a libfftw3_threads.a libfftw3.a ; strip pss.x
ifort -O3 -no-prec-div -recursive -openmp -static -o pssrigid.x pssrigid.o libtinker.a libfftw3_threads.a libfftw3.a ; strip pssrigid.x
ifort -O3 -no-prec-div -recursive -openmp -static -o pssrot.x pssrot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip pssrot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o radial.x radial.o libtinker.a libfftw3_threads.a libfftw3.a ; strip radial.x
ifort -O3 -no-prec-div -recursive -openmp -static -o saddle.x saddle.o libtinker.a libfftw3_threads.a libfftw3.a ; strip saddle.x
ifort -O3 -no-prec-div -recursive -openmp -static -o scan.x scan.o libtinker.a libfftw3_threads.a libfftw3.a ; strip scan.x
ifort -O3 -no-prec-div -recursive -openmp -static -o sniffer.x sniffer.o libtinker.a libfftw3_threads.a libfftw3.a ; strip sniffer.x
ifort -O3 -no-prec-div -recursive -openmp -static -o spacefill.x spacefill.o libtinker.a libfftw3_threads.a libfftw3.a ; strip spacefill.x
ifort -O3 -no-prec-div -recursive -openmp -static -o spectrum.x spectrum.o libtinker.a libfftw3_threads.a libfftw3.a ; strip spectrum.x
ifort -O3 -no-prec-div -recursive -openmp -static -o superpose.x superpose.o libtinker.a libfftw3_threads.a libfftw3.a ; strip superpose.x
ifort -O3 -no-prec-div -recursive -openmp -static -o testgrad.x testgrad.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testgrad.x
ifort -O3 -no-prec-div -recursive -openmp -static -o testhess.x testhess.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testhess.x
ifort -O3 -no-prec-div -recursive -openmp -static -o testpair.x testpair.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testpair.x
ifort -O3 -no-prec-div -recursive -openmp -static -o testpol.x testpol.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testpol.x
ifort -O3 -no-prec-div -recursive -openmp -static -o testrot.x testrot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testrot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o timer.x timer.o libtinker.a libfftw3_threads.a libfftw3.a ; strip timer.x
ifort -O3 -no-prec-div -recursive -openmp -static -o timerot.x timerot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip timerot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o torsfit.x torsfit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip torsfit.x
ifort -O3 -no-prec-div -recursive -openmp -static -o valence.x valence.o libtinker.a libfftw3_threads.a libfftw3.a ; strip valence.x
ifort -O3 -no-prec-div -recursive -openmp -static -o vibbig.x vibbig.o libtinker.a libfftw3_threads.a libfftw3.a ; strip vibbig.x
ifort -O3 -no-prec-div -recursive -openmp -static -o vibrate.x vibrate.o libtinker.a libfftw3_threads.a libfftw3.a ; strip vibrate.x
ifort -O3 -no-prec-div -recursive -openmp -static -o vibrot.x vibrot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip vibrot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o xtalfit.x xtalfit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xtalfit.x
ifort -O3 -no-prec-div -recursive -openmp -static -o xtalmin.x xtalmin.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xtalmin.x
ifort -O3 -no-prec-div -recursive -openmp -static -o xyzedit.x xyzedit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzedit.x
ifort -O3 -no-prec-div -recursive -openmp -static -o xyzint.x xyzint.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzint.x
ifort -O3 -no-prec-div -recursive -openmp -static -o xyzmol2.x xyzmol2.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzmol2.x
ifort -O3 -no-prec-div -recursive -openmp -static -o xyzpdb.x xyzpdb.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzpdb.x
