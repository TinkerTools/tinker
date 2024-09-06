#
#
#  ###################################################################
#  ##                                                               ##
#  ##  linkstatic.make  --  link static version of Tinker programs  ##
#  ##               (Intel Fortran for Linux Version)               ##
#  ##                                                               ##
#  ###################################################################
#
#
ifort -O3 -no-prec-div -recursive -openmp -static -o alchemy.x alchemy.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip alchemy.x
ifort -O3 -no-prec-div -recursive -openmp -static -o analyze.x analyze.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip analyze.x
ifort -O3 -no-prec-div -recursive -openmp -static -o anneal.x anneal.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip anneal.x
ifort -O3 -no-prec-div -recursive -openmp -static -o arcedit.x arcedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip arcedit.x
ifort -O3 -no-prec-div -recursive -openmp -static -o bar.x bar.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip bar.x
ifort -O3 -no-prec-div -recursive -openmp -static -o correlate.x correlate.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip correlate.x
ifort -O3 -no-prec-div -recursive -openmp -static -o critical.x critical.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip critical.x
ifort -O3 -no-prec-div -recursive -openmp -static -o crystal.x crystal.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip crystal.x
ifort -O3 -no-prec-div -recursive -openmp -static -o diffuse.x diffuse.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip diffuse.x
ifort -O3 -no-prec-div -recursive -openmp -static -o distgeom.x distgeom.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip distgeom.x
ifort -O3 -no-prec-div -recursive -openmp -static -o document.x document.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip document.x
ifort -O3 -no-prec-div -recursive -openmp -static -o dynamic.x dynamic.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip dynamic.x
ifort -O3 -no-prec-div -recursive -openmp -static -o freefix.x freefix.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip freefix.x
ifort -O3 -no-prec-div -recursive -openmp -static -o gda.x gda.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip gda.x
ifort -O3 -no-prec-div -recursive -openmp -static -o intedit.x intedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip intedit.x
ifort -O3 -no-prec-div -recursive -openmp -static -o intxyz.x intxyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip intxyz.x
ifort -O3 -no-prec-div -recursive -openmp -static -o minimize.x minimize.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minimize.x
ifort -O3 -no-prec-div -recursive -openmp -static -o minirot.x minirot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minirot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o minrigid.x minrigid.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minrigid.x
ifort -O3 -no-prec-div -recursive -openmp -static -o mol2xyz.x mol2xyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip mol2xyz.x
ifort -O3 -no-prec-div -recursive -openmp -static -o molxyz.x molxyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip molxyz.x
ifort -O3 -no-prec-div -recursive -openmp -static -o monte.x monte.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip monte.x
ifort -O3 -no-prec-div -recursive -openmp -static -o newton.x newton.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip newton.x
ifort -O3 -no-prec-div -recursive -openmp -static -o newtrot.x newtrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip newtrot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o nucleic.x nucleic.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip nucleic.x
ifort -O3 -no-prec-div -recursive -openmp -static -o optimize.x optimize.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optimize.x
ifort -O3 -no-prec-div -recursive -openmp -static -o optirot.x optirot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optirot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o optrigid.x optrigid.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optrigid.x
ifort -O3 -no-prec-div -recursive -openmp -static -o path.x path.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip path.x
ifort -O3 -no-prec-div -recursive -openmp -static -o pdbxyz.x pdbxyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pdbxyz.x
ifort -O3 -no-prec-div -recursive -openmp -static -o polarize.x polarize.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip polarize.x
ifort -O3 -no-prec-div -recursive -openmp -static -o poledit.x poledit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip poledit.x
ifort -O3 -no-prec-div -recursive -openmp -static -o potential.x potential.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip potential.x
ifort -O3 -no-prec-div -recursive -openmp -static -o prmedit.x prmedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip prmedit.x
ifort -O3 -no-prec-div -recursive -openmp -static -o protein.x protein.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip protein.x
ifort -O3 -no-prec-div -recursive -openmp -static -o pss.x pss.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pss.x
ifort -O3 -no-prec-div -recursive -openmp -static -o pssrigid.x pssrigid.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pssrigid.x
ifort -O3 -no-prec-div -recursive -openmp -static -o pssrot.x pssrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pssrot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o radial.x radial.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip radial.x
ifort -O3 -no-prec-div -recursive -openmp -static -o saddle.x saddle.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip saddle.x
ifort -O3 -no-prec-div -recursive -openmp -static -o scan.x scan.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip scan.x
ifort -O3 -no-prec-div -recursive -openmp -static -o sniffer.x sniffer.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip sniffer.x
ifort -O3 -no-prec-div -recursive -openmp -static -o spacefill.x spacefill.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip spacefill.x
ifort -O3 -no-prec-div -recursive -openmp -static -o spectrum.x spectrum.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip spectrum.x
ifort -O3 -no-prec-div -recursive -openmp -static -o superpose.x superpose.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip superpose.x
ifort -O3 -no-prec-div -recursive -openmp -static -o testgrad.x testgrad.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testgrad.x
ifort -O3 -no-prec-div -recursive -openmp -static -o testhess.x testhess.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testhess.x
ifort -O3 -no-prec-div -recursive -openmp -static -o testpair.x testpair.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testpair.x
ifort -O3 -no-prec-div -recursive -openmp -static -o testpol.x testpol.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testpol.x
ifort -O3 -no-prec-div -recursive -openmp -static -o testrot.x testrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testrot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o testvir.x testvir.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testvir.x
ifort -O3 -no-prec-div -recursive -openmp -static -o timer.x timer.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip timer.x
ifort -O3 -no-prec-div -recursive -openmp -static -o timerot.x timerot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip timerot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o torsfit.x torsfit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip torsfit.x
ifort -O3 -no-prec-div -recursive -openmp -static -o valence.x valence.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip valence.x
ifort -O3 -no-prec-div -recursive -openmp -static -o vibbig.x vibbig.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibbig.x
ifort -O3 -no-prec-div -recursive -openmp -static -o vibrate.x vibrate.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibrate.x
ifort -O3 -no-prec-div -recursive -openmp -static -o vibrot.x vibrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibrot.x
ifort -O3 -no-prec-div -recursive -openmp -static -o xtalfit.x xtalfit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xtalfit.x
ifort -O3 -no-prec-div -recursive -openmp -static -o xtalmin.x xtalmin.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xtalmin.x
ifort -O3 -no-prec-div -recursive -openmp -static -o xyzedit.x xyzedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzedit.x
ifort -O3 -no-prec-div -recursive -openmp -static -o xyzint.x xyzint.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzint.x
ifort -O3 -no-prec-div -recursive -openmp -static -o xyzmol2.x xyzmol2.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzmol2.x
ifort -O3 -no-prec-div -recursive -openmp -static -o xyzpdb.x xyzpdb.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzpdb.x
