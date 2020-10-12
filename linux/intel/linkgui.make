#
#
#  ##############################################################
#  ##                                                          ##
#  ##  linkgui.make  --  link the Tinker programs for FFE use  ##
#  ##             (Intel Fortran for Linux Version)            ##
#  ##                                                          ##
#  ##############################################################
#
#
#  Copy the Java libjvm.so library from ../jre/lib/i386/client
#  in the standard Java distribution to use the script below,
#  or add it to the library search path
#
#
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o alchemy.x alchemy.o libtinker.a libjvm.so ; strip alchemy.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o analyze.x analyze.o libtinker.a libjvm.so ; strip analyze.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o anneal.x anneal.o libtinker.a libjvm.so ; strip anneal.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o archive.x archive.o libtinker.a libjvm.so ; strip archive.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o bar.x bar.o libtinker.a libjvm.so ; strip bar.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o correlate.x correlate.o libtinker.a libjvm.so ; strip correlate.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o critical.x critical.o libtinker.a libjvm.so ; strip critical.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o crystal.x crystal.o libtinker.a libjvm.so ; strip crystal.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o diffuse.x diffuse.o libtinker.a libjvm.so ; strip diffuse.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o distgeom.x distgeom.o libtinker.a libjvm.so ; strip distgeom.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o document.x document.o libtinker.a libjvm.so ; strip document.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o dynamic.x dynamic.o libtinker.a libjvm.so ; strip dynamic.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o gda.x gda.o libtinker.a libjvm.so ; strip gda.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o intedit.x intedit.o libtinker.a libjvm.so ; strip intedit.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o intxyz.x intxyz.o libtinker.a libjvm.so ; strip intxyz.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o minimize.x minimize.o libtinker.a libjvm.so ; strip minimize.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o minirot.x minirot.o libtinker.a libjvm.so ; strip minirot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o minrigid.x minrigid.o libtinker.a libjvm.so ; strip minrigid.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o mol2xyz.x mol2xyz.o libtinker.a libjvm.so ; strip mol2xyz.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o molxyz.x molxyz.o libtinker.a libjvm.so ; strip molxyz.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o monte.x monte.o libtinker.a libjvm.so ; strip monte.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o newton.x newton.o libtinker.a libjvm.so ; strip newton.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o newtrot.x newtrot.o libtinker.a libjvm.so ; strip newtrot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o nucleic.x nucleic.o libtinker.a libjvm.so ; strip nucleic.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o optimize.x optimize.o libtinker.a libjvm.so ; strip optimize.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o optirot.x optirot.o libtinker.a libjvm.so ; strip optirot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o optrigid.x optrigid.o libtinker.a libjvm.so ; strip optrigid.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o path.x path.o libtinker.a libjvm.so ; strip path.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o pdbxyz.x pdbxyz.o libtinker.a libjvm.so ; strip pdbxyz.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o polarize.x polarize.o libtinker.a libjvm.so ; strip polarize.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o poledit.x poledit.o libtinker.a libjvm.so ; strip poledit.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o potential.x potential.o libtinker.a libjvm.so ; strip potential.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o prmedit.x prmedit.o libtinker.a libjvm.so ; strip prmedit.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o protein.x protein.o libtinker.a libjvm.so ; strip protein.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o pss.x pss.o libtinker.a libjvm.so ; strip pss.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o pssrigid.x pssrigid.o libtinker.a libjvm.so ; strip pssrigid.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o pssrot.x pssrot.o libtinker.a libjvm.so ; strip pssrot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o radial.x radial.o libtinker.a libjvm.so ; strip radial.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o saddle.x saddle.o libtinker.a libjvm.so ; strip saddle.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o scan.x scan.o libtinker.a libjvm.so ; strip scan.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o sniffer.x sniffer.o libtinker.a libjvm.so ; strip sniffer.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o spacefill.x spacefill.o libtinker.a libjvm.so ; strip spacefill.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o spectrum.x spectrum.o libtinker.a libjvm.so ; strip spectrum.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o superpose.x superpose.o libtinker.a libjvm.so ; strip superpose.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o testgrad.x testgrad.o libtinker.a libjvm.so ; strip testgrad.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o testhess.x testhess.o libtinker.a libjvm.so ; strip testhess.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o testpair.x testpair.o libtinker.a libjvm.so ; strip testpair.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o testpol.x testpol.o libtinker.a libjvm.so ; strip testpol.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o testrot.x testrot.o libtinker.a libjvm.so ; strip testrot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o testvir.x testvir.o libtinker.a libjvm.so ; strip testvir.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o timer.x timer.o libtinker.a libjvm.so ; strip timer.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o timerot.x timerot.o libtinker.a libjvm.so ; strip timerot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o torsfit.x torsfit.o libtinker.a libjvm.so ; strip torsfit.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o valence.x valence.o libtinker.a libjvm.so ; strip valence.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o vibbig.x vibbig.o libtinker.a libjvm.so ; strip vibbig.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o vibrate.x vibrate.o libtinker.a libjvm.so ; strip vibrate.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o vibrot.x vibrot.o libtinker.a libjvm.so ; strip vibrot.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o xtalfit.x xtalfit.o libtinker.a libjvm.so ; strip xtalfit.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o xtalmin.x xtalmin.o libtinker.a libjvm.so ; strip xtalmin.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o xyzedit.x xyzedit.o libtinker.a libjvm.so ; strip xyzedit.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o xyzint.x xyzint.o libtinker.a libjvm.so ; strip xyzint.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o xyzmol2.x xyzmol2.o libtinker.a libjvm.so ; strip xyzmol2.x
ifort -O3 -no-prec-div -recursive -openmp -static-libgcc -static-intel -o xyzpdb.x xyzpdb.o libtinker.a libjvm.so ; strip xyzpdb.x
