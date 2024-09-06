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
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o alchemy.x alchemy.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip alchemy.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o analyze.x analyze.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip analyze.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o anneal.x anneal.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip anneal.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o arcedit.x arcedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip arcedit.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o bar.x bar.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip bar.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o correlate.x correlate.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip correlate.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o critical.x critical.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip critical.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o crystal.x crystal.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip crystal.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o diffuse.x diffuse.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip diffuse.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o distgeom.x distgeom.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip distgeom.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o document.x document.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip document.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o dynamic.x dynamic.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip dynamic.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o freefix.x freefix.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip freefix.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o gda.x gda.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip gda.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o intedit.x intedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip intedit.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o intxyz.x intxyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip intxyz.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o minimize.x minimize.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip minimize.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o minirot.x minirot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip minirot.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o minrigid.x minrigid.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip minrigid.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o mol2xyz.x mol2xyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip mol2xyz.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o molxyz.x molxyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip molxyz.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o monte.x monte.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip monte.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o newton.x newton.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip newton.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o newtrot.x newtrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip newtrot.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o nucleic.x nucleic.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip nucleic.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o optimize.x optimize.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip optimize.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o optirot.x optirot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip optirot.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o optrigid.x optrigid.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip optrigid.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o path.x path.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip path.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o pdbxyz.x pdbxyz.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip pdbxyz.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o polarize.x polarize.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip polarize.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o poledit.x poledit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip poledit.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o potential.x potential.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip potential.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o prmedit.x prmedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip prmedit.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o protein.x protein.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip protein.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o pss.x pss.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip pss.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o pssrigid.x pssrigid.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip pssrigid.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o pssrot.x pssrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip pssrot.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o radial.x radial.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip radial.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o saddle.x saddle.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip saddle.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o scan.x scan.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip scan.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o sniffer.x sniffer.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip sniffer.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o spacefill.x spacefill.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip spacefill.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o spectrum.x spectrum.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip spectrum.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o superpose.x superpose.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip superpose.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o testgrad.x testgrad.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip testgrad.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o testhess.x testhess.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip testhess.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o testpair.x testpair.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip testpair.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o testpol.x testpol.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip testpol.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o testrot.x testrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip testrot.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o testsurf.x testsurf.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip testsurf.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o testvir.x testvir.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip testvir.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o timer.x timer.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip timer.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o timerot.x timerot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip timerot.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o torsfit.x torsfit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip torsfit.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o valence.x valence.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip valence.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o vibbig.x vibbig.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip vibbig.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o vibrate.x vibrate.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip vibrate.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o vibrot.x vibrot.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip vibrot.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o xtalfit.x xtalfit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip xtalfit.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o xtalmin.x xtalmin.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip xtalmin.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o xyzedit.x xyzedit.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip xyzedit.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o xyzint.x xyzint.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip xyzint.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o xyzmol2.x xyzmol2.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip xyzmol2.x
ifort -O3 -no-prec-div -recursive -qopenmp -static-libgcc -static-intel -o xyzpdb.x xyzpdb.o -L. -L../lib/linux -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 libjvm.so ; strip xyzpdb.x
